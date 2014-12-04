#include "../copyright.h"

/*=============================================================================
 * FILE: update_starparticles.c
 *
 * PURPOSE: 1. update particles mass and lifetime
 *          2. update the gas where the particles reside
 *
 * CONTAINS PUBLIC FUNCTIONS:
 *   update_starparticles();
 *   age_starparticles();
 *
 * CONTAINS PRIVATE FUNCTIONS:
 *   mass_inflow();
 *
 * IMPORTANT NOTE:
 *   As star particles are located at [0,1] zone away from the Grid edge,
 *   the star particles' mass and momentum are updated with the 1st-order fluxes
 *   in the ghost zones (VL integrator doesn't calculate the 2nd-order flux
 *   in the ghost zones).  This will lead to a tiny difference between parallel 
 *   run and single processor run in velocity fields.
 *   
 * HISTORY:
 *   Written by Hao Gong, Aug. 2011
 *   Modified by Aaron Skinner, Apr. 2013
 *   Modified by Chang-Goo Kim, Jan. 2014 (for non-isothermal EOS)
 *============================================================================*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../athena.h"
#include "../defs.h"
#include "../prototypes.h"
#include "../globals.h"


#ifdef STAR_PARTICLE
//#define CURRENT

/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 *============================================================================*/
void mass_inflow(GridS *pG, StarParS *pStar, Cons1DS ***x1Flux,
                 Cons1DS ***x2Flux, Cons1DS ***x3Flux, Real *dm, Real *dM1,
                 Real *dM2, Real *dM3);
void outflow_flux(GridS *pG, ConsS ***U, Cons1DS ***x1Flux, Cons1DS ***x2Flux,
                  Cons1DS ***x3Flux);

/*=========================== PUBLIC FUNCTIONS ===============================*/
/*----------------------------------------------------------------------------*/

void update_starparticles(GridS *pG, Cons1DS ***x1Flux, Cons1DS ***x2Flux,
                          Cons1DS ***x3Flux)
{
  StarParListS *pList=NULL;
  StarParS *pStar=NULL;
  int ip,jp,kp;
  Real minv,dm,dM1,dM2,dM3;
  
  pList = pG->Lstars;
  while (pList) {
    pStar = &(pList->starpar);
    cc_ijk(pG,pStar->x1,pStar->x2,pStar->x3,&ip,&jp,&kp);
    
    mass_inflow(pG,pStar,x1Flux,x2Flux,x3Flux,&dm,&dM1,&dM2,&dM3);
    if (dm > 0.0) {
      minv = 1.0/(pStar->m+dm);
      pStar->v1 = (pStar->m*pStar->v1+dM1)*minv;
      pStar->v2 = (pStar->m*pStar->v2+dM2)*minv;
      pStar->v3 = (pStar->m*pStar->v3+dM3)*minv;
      pStar->m += dm;
      pStar->mdot = dm/pG->dt;
    }
    if(pStar->m < 0) ath_error("[update_SP] mass cannot be negative!\n");
    
    pList = pList->next;
  }
  return;
}

void age_starparticles(GridS *pG)
{
  StarParListS *pList=NULL;
  
  pList = pG->Lstars;
  while (pList) {
//    pList->starpar.age += pG->dt;
    pList->starpar.age += pG->dt_old; // it is moved to the integrator loop; after new_dt
    pList = pList->next;
  }
  
  return;
}

/*=========================== PRIVATE FUNCTIONS ==============================*/
/*----------------------------------------------------------------------------*/

void mass_inflow(GridS *pG, StarParS *pStar,
                 Cons1DS ***x1Flux, Cons1DS ***x2Flux, Cons1DS ***x3Flux,
                 Real *dm, Real *dM1, Real *dM2, Real *dM3)
{
  int i,j,k;
  int ip,jp,kp;
  Real x1=pStar->x1,x2=pStar->x2,x3=pStar->x3;
  
  /* update particle mass & momentum using term \delta rho * v */
  int ip_old,jp_old,kp_old;
  Real m=0.0, m_old=0.0;
  Real M1=0.0, M2=0.0, M3=0.0;
  Real M1_old=0.0, M2_old=0.0, M3_old=0.0;
  Real dV, tmp;
  
  dV = pG->dx1*pG->dx2*pG->dx3;

  cc_ijk(pG,pStar->x1,pStar->x2,pStar->x3,&ip,&jp,&kp);
  cc_ijk(pG,pStar->x1_old,pStar->x2_old,pStar->x3_old,&ip_old,&jp_old,&kp_old);
  
  if ((ip != ip_old) || (jp != jp_old) || (kp != kp_old)) {
/* previous implementation. 
    for (k=kp-NSINK_STARP; k<=kp+NSINK_STARP; k++) {
      for (j=jp-NSINK_STARP; j<=jp+NSINK_STARP; j++) {
        for (i=ip-NSINK_STARP; i<=ip+NSINK_STARP; i++) {
          m  += pG->U[k][j][i].d;
          M1 += pG->U[k][j][i].M1;
          M2 += pG->U[k][j][i].M2;
          M3 += pG->U[k][j][i].M3;
        }
      }
    }
    
    for (k=kp_old-NSINK_STARP; k<=kp_old+NSINK_STARP; k++) {
      for (j=jp_old-NSINK_STARP; j<=jp_old+NSINK_STARP; j++) {
        for (i=ip_old-NSINK_STARP; i<=ip_old+NSINK_STARP; i++) {
          m_old  += pG->U[k][j][i].d;
          M1_old += pG->U[k][j][i].M1;
          M2_old += pG->U[k][j][i].M2;
          M3_old += pG->U[k][j][i].M3;
        }
      }
    }
*/
    dV   = pG->dx1*pG->dx2*pG->dx3;
    *dm  = -pStar->mghost*dV;
    *dM1 = -pStar->M1ghost*dV;
    *dM2 = -pStar->M2ghost*dV;
    *dM3 = -pStar->M3ghost*dV;
  } else {
    *dm = 0.0;
    *dM1 = 0.0;
    *dM2 = 0.0;
    *dM3 = 0.0;
  }
 

  /* update the mass of star particle using the flux of the outer face of the
   * ghost zones, rho*\delta \dot v */  
  tmp = pG->dt*pG->dx2*pG->dx3;
  for (k=kp-NSINK_STARP; k<=kp+NSINK_STARP; k++) {
    for (j=jp-NSINK_STARP; j<=jp+NSINK_STARP; j++) {
      *dm  -= tmp*(x1Flux[k][j][ip+NSINK_STARP+1].d  -
                   x1Flux[k][j][ip-NSINK_STARP  ].d );
      *dM1 -= tmp*(x1Flux[k][j][ip+NSINK_STARP+1].Mx -
                   x1Flux[k][j][ip-NSINK_STARP  ].Mx);
      *dM2 -= tmp*(x1Flux[k][j][ip+NSINK_STARP+1].My -
                   x1Flux[k][j][ip-NSINK_STARP  ].My);
      *dM3 -= tmp*(x1Flux[k][j][ip+NSINK_STARP+1].Mz -
                   x1Flux[k][j][ip-NSINK_STARP  ].Mz);
    }
  }
  
  tmp = pG->dt*pG->dx1*pG->dx3;
  for (k=kp-NSINK_STARP; k<=kp+NSINK_STARP; k++) {
    for (i=ip-NSINK_STARP; i<=ip+NSINK_STARP; i++) {
      *dm  -= tmp*(x2Flux[k][jp+NSINK_STARP+1][i].d  -
                   x2Flux[k][jp-NSINK_STARP  ][i].d );
      *dM1 -= tmp*(x2Flux[k][jp+NSINK_STARP+1][i].Mz -
                   x2Flux[k][jp-NSINK_STARP  ][i].Mz);
      *dM2 -= tmp*(x2Flux[k][jp+NSINK_STARP+1][i].Mx -
                   x2Flux[k][jp-NSINK_STARP  ][i].Mx);
      *dM3 -= tmp*(x2Flux[k][jp+NSINK_STARP+1][i].My -
                   x2Flux[k][jp-NSINK_STARP  ][i].My);
    }
  }
  
  tmp = pG->dt*pG->dx1*pG->dx2;
  for (j=jp-NSINK_STARP; j<=jp+NSINK_STARP; j++) {
    for (i=ip-NSINK_STARP; i<=ip+NSINK_STARP; i++) {
      *dm  -= tmp*(x3Flux[kp+NSINK_STARP+1][j][i].d  -
                   x3Flux[kp-NSINK_STARP  ][j][i].d );
      *dM1 -= tmp*(x3Flux[kp+NSINK_STARP+1][j][i].My -
                   x3Flux[kp-NSINK_STARP  ][j][i].My);
      *dM2 -= tmp*(x3Flux[kp+NSINK_STARP+1][j][i].Mz -
                   x3Flux[kp-NSINK_STARP  ][j][i].Mz);
      *dM3 -= tmp*(x3Flux[kp+NSINK_STARP+1][j][i].Mx -
                   x3Flux[kp-NSINK_STARP  ][j][i].Mx);
    }
  }
    
  return;
}

#endif /* STAR_PARTICLE */
