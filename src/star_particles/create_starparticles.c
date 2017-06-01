#include "../copyright.h"

/*=============================================================================
 * FILE: create_starparticles.c
 *
 * PURPOSE: 1. Check if any cell satisfies criteria for star particle creation
 *             a. Non-strict star particle creation criteria:
 *                Exceeding density threshold + local grav minimum
 *             b. Strict star particle creation criteria:
 *                density threshold+ grav minimum + converging flow + grav bound
 *          2. Create star particle if criteria are satisfied
 *
 * CONTAINS PUBLIC FUNCTIONS:
 *   create_starparticles();
 *
 * CONTAINS PRIVATE FUNCTIONS:
 *   converging_flow();
 *   local_phi_min();
 *   grav_bound();
 *   push_to_local_list();
 *
 * HISTORY:
 *   Written by Hao Gong, July 2011
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


#define OVERLAP
#define CONVERGING_ALL_DIRECTION

#ifdef STAR_PARTICLE

/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 *   converging_flow()    - check if flow around cell is converging
 *   local_phi_min()      - check if cell is local potential minimum
 *   grav_bound()         - check if gas in control vol is gravitionally bound
 *   push_to_local_list() - push the particle just created to the local list
 *============================================================================*/
int  density_crit(GridS *pG, int ic, int jc, int kc);
int  converging_flow(GridS *pG, int ic, int jc, int kc);
#ifdef SELF_GRAVITY
int  local_phi_min(GridS *pG, int ic, int jc, int kc);
int  grav_bound();
#endif
void push_to_local_list(GridS *pG, int ic, int jc, int kc, int starpar_id);
Real Stromgren_radius(GridS *pG, const int ip, const int jp, const int kp,
                      StarParS *pStar);


/*=========================== PUBLIC FUNCTIONS ===============================*/
/*----------------------------------------------------------------------------*/

void create_starparticles(DomainS *pD)
{
  GridS *pG=pD->Grid;
  StarParListS *pLstars=NULL;
  StarParListS *pGstars=NULL;
  StarParS *pStar=NULL;
  int i,j,k,is,ie,js,je,ks,ke;
  int new_star_particle_flag;
  int starpar_id;
  Real rdp,v1,v2,v3,dm,drho,vdiff;
  Real dx,dr;
  Real xp,yp,zp;

  is = pG->is; ie = pG->ie;
  js = pG->js; je = pG->je;
  ks = pG->ks; ke = pG->ke;

  dx = MAX(MAX(pG->dx1,pG->dx2),pG->dx3);
  dr = (double)((NGHOST_STARP+0.5)*dx);
  
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        /* a. Check Truelove density criteria */
        if (density_crit(pG,i,j,k)) {
          /* b. Check if local potential minimum */
#ifdef SELF_GRAVITY
          if (local_phi_min(pG,i,j,k)) {
#endif
#ifdef STRICT_STARPARTICLE_CREATION
            /* c. Check if converging flow */
            if (converging_flow(pG,i,j,k)) {
              /*  d. Check if accretion zone is gravitationally bound */
#ifdef SELF_GRAVITY
              if (grav_bound(pG,i,j,k)) {
#endif
#endif
                pLstars = pG->Lstars;
                pGstars = pG->Gstars;
                /* Generate unique identifier (global Domain position) */
                starpar_id = (i-is+pG->Disp[0]) +
                  pD->Nx[0]*((j-js+pG->Disp[1]) + pD->Nx[1]*(k-ks+pG->Disp[2]));
                if (pGstars == NULL) {
                  push_to_local_list(pG,i,j,k,starpar_id);
                }
                else {
                  new_star_particle_flag = 1;
                  /* e. Check if within existing star particle's control volume */
                  while (pGstars && new_star_particle_flag) {
                    pStar = &(pGstars->starpar);
                    cc_pos(pG,i,j,k,&xp,&yp,&zp);
                    rdp = sqrt((double)(SQR(xp-pStar->x1)+SQR(yp-pStar->x2)+SQR(zp-pStar->x3)));
                    
                    if (rdp <= 2*dr) new_star_particle_flag = 0;
                    pGstars = pGstars->next;
                  }
                  if (new_star_particle_flag)
                    push_to_local_list(pG,i,j,k,starpar_id);
                }
#ifdef STRICT_STARPARTICLE_CREATION
#ifdef SELF_GRAVITY
              }
#endif
            }
#endif
#ifdef SELF_GRAVITY
          }
#endif
        }
/* add this line if rarefaction happens in the gas */
/*       if (pG->U[k][j][i].d <= 1e-5) pG->U[k][j][i].d=1e-3; */
      }
    }
  }

  return;
}


/*=========================== PRIVATE FUNCTIONS ==============================*/
/*----------------------------------------------------------------------------*/
//#define TRUELOVE_DENSITY_CRITERION
#define LARSON_PENSTON_DENSITY_CRITERION
int density_crit(GridS *pG, int ic, int jc, int kc)
{
  Real tmp,grid_rho_crit = rho_crit;
  Real dx = MAX(MAX(pG->dx1,pG->dx2),pG->dx3);
  Real csound2,Press;
  ConsS U;
#ifndef SELF_GRAVITY
  Real four_pi_G=4.0*PI*pG->units.G;
#endif

  U = pG->U[kc][jc][ic];
#ifndef BAROTROPIC
  Press = U.E - 0.5*(SQR(U.M1)+SQR(U.M2)+SQR(U.M3))/U.d;
#ifdef MHD
  Press -= (0.5)*(SQR(U.B1c)+SQR(U.B2c)+SQR(U.B3c));
#endif
  Press *= Gamma_1;
  csound2 = Press/U.d;
#else /* BAROTROPIC */
  csound2 = Iso_csound2;
#endif /* BAROTROPIC */
  
#ifdef TRUELOVE_DENSITY_CRITERION
  tmp = SQR(PI)*csound2/(4.0*four_pi_G*SQR(dx));
  grid_rho_crit = MAX(grid_rho_crit,tmp);
#endif
#ifdef LARSON_PENSTON_DENSITY_CRITERION
  tmp = 8.86*csound2*4.0/(four_pi_G*SQR(dx));
  grid_rho_crit = MAX(grid_rho_crit,tmp);
#endif
  
  return (pG->U[kc][jc][ic].d > grid_rho_crit);
}

/*----------------------------------------------------------------------------*/

int converging_flow(GridS *pG, int ic, int jc, int kc)
{
  int i,j,k,is,ie,js,je,ks,ke;
  Real divV=0.0;
  Real divVx=0.0;
  Real divVy=0.0;
  Real divVz=0.0;
  
  /* A.S. NOTE:  CHANGE THE +/- 1 TO SOME GLOBAL CONSTANT?  
   * (AND DEFINE THIS IN GLOBALS.H/MAIN.C)?
   */
  
  is = ic-1; ie = ic+1;
  js = jc-1; je = jc+1;
  ks = kc-1; ke = kc+1;
  
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      divVx += pG->U[k][j][ie].M1/pG->U[k][j][ie].d -
        pG->U[k][j][is].M1/pG->U[k][j][is].d;
    }
  }
  
  for (k=ks; k<=ke; k++) {
    for (i=is; i<=ie; i++) {
      divVy += pG->U[k][je][i].M2/pG->U[k][je][i].d -
        pG->U[k][js][i].M2/pG->U[k][js][i].d;
    }
  }
  
  for (j=js; j<=je; j++) {
    for (i=is; i<=ie; i++) {
      divVz += pG->U[ke][j][i].M3/pG->U[ke][j][i].d -
        pG->U[ks][j][i].M3/pG->U[ks][j][i].d;
    }
  }
  
  divV = divVx + divVy +divVz;
#ifdef CONVERGING_ALL_DIRECTION
  return ((divVx < 0.0) && (divVy < 0.0) && (divVz < 0.0));
#else
  return (divV < 0.0);
#endif
}

/*----------------------------------------------------------------------------*/

#ifdef SELF_GRAVITY
int local_phi_min(GridS *pG, int ic, int jc, int kc)
{
  int i,j,k,is,ie,js,je,ks,ke;
  int iflag=0;
  
  is = ic-1; ie = ic+1;
  js = jc-1; je = jc+1;
  ks = kc-1; ke = kc+1;
    
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        if (pG->Phi[kc][jc][ic] < pG->Phi[k][j][i]) iflag++;
      }
    }
  }
  
  return (iflag >= (ie-is+1)*(je-js+1)*(ke-ks+1) - 1);
}

/*----------------------------------------------------------------------------*/

int grav_bound(GridS *pG, int ic, int jc, int kc)
{
  int i,j,k,is,ie,js,je,ks,ke;
  Real v01,v02,v03,v1,v2,v3;
  Real phi_edge=0.0,phi_total=0.0,ke_total=0.0,te_total=0.0,eng_total=0.0;
  Real dinv;
#ifndef BAROTROPIC
  Real eint;
  ConsS U;
#endif
#ifdef RADIATION
  Real re_total=0.0,cocrad=c/crad;
#endif
  
  is = ic-1; ie = ic+1;
  js = jc-1; je = jc+1;
  ks = kc-1; ke = kc+1;
  
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      phi_edge += pG->Phi[k][j][is-1] + pG->Phi[k][j][ie+1];
    }
  }
  for (k=ks; k<=ke; k++) {
    for (i=is-1; i<=ie+1; i++) {
      phi_edge += pG->Phi[k][js-1][i] + pG->Phi[k][je+1][i];
    }
  }
  for (j=js-1; j<=je+1; j++) {
    for (i=is-1; i<=ie+1; i++) {
      phi_edge += pG->Phi[ks-1][j][i] + pG->Phi[ke+1][j][i];
    }
  }
  phi_edge /= 98.0;
  
  dinv = 1.0/pG->U[kc][jc][ic].d;
  v01 = pG->U[kc][jc][ic].M1*dinv;
  v02 = pG->U[kc][jc][ic].M2*dinv;
  v03 = pG->U[kc][jc][ic].M3*dinv;
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<ie; i++) {
        phi_total += 0.5*pG->U[k][j][i].d*(pG->Phi[k][j][i] - phi_edge);
        dinv = 1.0/pG->U[k][j][i].d;
        v1 = pG->U[k][j][i].M1*dinv;
        v2 = pG->U[k][j][i].M2*dinv;
        v3 = pG->U[k][j][i].M3*dinv;
        ke_total += pG->U[k][j][i].d*(SQR(v1-v01) + SQR(v2-v02) + SQR(v3-v03));
#ifndef BAROTROPIC
        U = pG->U[k][j][i];
        eint = U.E - 0.5*(SQR(U.M1)+SQR(U.M2)+SQR(U.M3))/U.d;
#ifdef MHD
        eint -= (0.5)*(SQR(U.B1c)+SQR(U.B2c)+SQR(U.B3c));
#endif
        te_total += eint;
#else /* BAROTROPIC */
        te_total += 1.5*pG->U[k][j][i].d*Iso_csound2;
#endif /* BAROTROPIC */
#ifdef RADIATION
        re_total += cocrad*pG->Urad[k][j][i].Er;
#endif
      }
    }
  }

  eng_total = phi_total + ke_total + te_total;
#ifdef RADIATION
  eng_total += re_total;
#endif
  
  return (eng_total < 0.0);
}
#endif

/*----------------------------------------------------------------------------*/

void push_to_local_list(GridS *pG, int ip, int jp, int kp, int starpar_id)
{  
  StarParListS *pGstars=NULL;
  StarParListS *pList=NULL;
  StarParS *pStar=NULL;
  StarParS *pStar_tmp=NULL;
  int i,j,k;
  int ip_tmp,jp_tmp,kp_tmp,Noverlap=0;
  Real xp,yp,zp;
  Real rhot,M1t,M2t,M3t,v1,v2,v3,dV;
  
  if ((pList = (StarParListS *) calloc(1,sizeof(StarParListS))) == NULL)
    ath_error("[push_to_local_list]: Error callocing memory for starparlist\n");
  
  pStar = &(pList->starpar);
  pStar->id = starpar_id;
  pStar->merge_history = 0;
  pStar->isnew = 1;

  cc_pos(pG,ip,jp,kp,&xp,&yp,&zp);
  
  rhot = 0.0;
  M1t = 0.0;
  M2t = 0.0;
  M3t = 0.0;
  dV = pG->dx1*pG->dx2*pG->dx3;

  for (k=kp-NSINK_STARP; k<=kp+NSINK_STARP; k++) {
    for (j=jp-NSINK_STARP; j<=jp+NSINK_STARP; j++) {
      for (i=ip-NSINK_STARP; i<=ip+NSINK_STARP; i++) {
        rhot += pG->U[k][j][i].d;
        M1t += pG->U[k][j][i].M1;
        M2t += pG->U[k][j][i].M2;
        M3t += pG->U[k][j][i].M3;
#ifdef OVERLAP
/* exclude overlapped control volume with existing star particles */
        pGstars=pG->Gstars;
        while(pGstars){
          pStar_tmp = &(pGstars->starpar);
          cc_ijk(pG,pStar_tmp->x1,pStar_tmp->x2,pStar_tmp->x3,&ip_tmp,&jp_tmp,&kp_tmp);
          if(((i>=ip_tmp-1) && (i<=ip_tmp+1)) &&
             ((j>=jp_tmp-1) && (j<=jp_tmp+1)) &&
             ((k>=kp_tmp-1) && (k<=kp_tmp+1))){
            rhot -= pG->U[k][j][i].d;
            M1t -= pG->U[k][j][i].M1;
            M2t -= pG->U[k][j][i].M2;
            M3t -= pG->U[k][j][i].M3;
            Noverlap++;
          }
          pGstars = pGstars->next;
        }
#endif
      }
    }
  }

  v1 = M1t/rhot;
  v2 = M2t/rhot;
  v3 = M3t/rhot;
 
  pStar->m = rhot*dV;
  pStar->v1 = v1;
  pStar->v2 = v2;
  pStar->v3 = v3;
  pStar->x1 = xp;
  pStar->x2 = yp;
  pStar->x3 = zp;
  pStar->x1_old = xp;
  pStar->x2_old = yp;
  pStar->x3_old = zp;
  pStar->age = 0.0;
  pStar->mdot = 0.0;
  pStar->mghost  = 0.0;
  pStar->M1ghost = 0.0;
  pStar->M2ghost = 0.0;
  pStar->M3ghost = 0.0;

  set_ghost_region(pG,pStar);  
  rhot = 0.0; M1t = 0.0; M2t = 0.0; M3t = 0.0;
  Noverlap = 0;
  for (k=kp-NSINK_STARP; k<=kp+NSINK_STARP; k++) {
    for (j=jp-NSINK_STARP; j<=jp+NSINK_STARP; j++) {
      for (i=ip-NSINK_STARP; i<=ip+NSINK_STARP; i++) {
        rhot += pG->U[k][j][i].d;
        M1t += pG->U[k][j][i].M1;
        M2t += pG->U[k][j][i].M2;
        M3t += pG->U[k][j][i].M3;

#if defined(MCTRACERS) || defined(VFTRACERS)
       flag_tracer_star(pG, ip, jp, kp, starpar_id);
#endif /* TRACERS */

#ifdef OVERLAP
/* exclude overlapped control volume with existing star particles */
        pGstars=pG->Gstars;
        while(pGstars){
          pStar_tmp = &(pGstars->starpar);
          cc_ijk(pG,pStar_tmp->x1,pStar_tmp->x2,pStar_tmp->x3,&ip_tmp,&jp_tmp,&kp_tmp);
          if(((i>=ip_tmp-1) && (i<=ip_tmp+1)) &&
             ((j>=jp_tmp-1) && (j<=jp_tmp+1)) &&
             ((k>=kp_tmp-1) && (k<=kp_tmp+1))){
            rhot -= pG->U[k][j][i].d;
            M1t -= pG->U[k][j][i].M1;
            M2t -= pG->U[k][j][i].M2;
            M3t -= pG->U[k][j][i].M3;
            Noverlap++;
          }
          pGstars = pGstars->next;
        }
#endif
      }
    }
  }

  pStar->mghost  = rhot;
  pStar->M1ghost = M1t;
  pStar->M2ghost = M2t;
  pStar->M3ghost = M3t;


#ifdef RADIATION
  Stromgren_radius(pG,ip,jp,kp,pStar);
#endif
  
  ath_perr(-1,"\tProc %d is creating a star of mass %f at t=%g i=(%d,%d,%d) with id=%d, nGstars=%d\n",
         myID_Comm_world,pStar->m,pG->time,ip,jp,kp,pStar->id,pG->nGstars);
  printf("\tProc %d is creating a star of mass %f at t=%g i=(%d,%d,%d) with id=%d, nGstars=%d, Noverlap=%d\n",
         myID_Comm_world,pStar->m,pG->time,ip,jp,kp,pStar->id,pG->nGstars,Noverlap);
  /* push it to the local list */
  starpar_push_local(pG,pList);
  
  starpar_printlist(-1, pG);
  
  return;
}


Real Stromgren_radius(GridS *pG, const int ip, const int jp, const int kp,
                      StarParS *pStar)
{
  Real RS, R0 = 3.31;
  Real rho_amb = 0.0;
  int i,j,k,N=0;
  
  for (k=kp-1; k<=kp+1; k++) {
    for (j=jp-1; j<=jp+1; j++) {
      rho_amb += pG->U[k][j][ip-2].d + pG->U[k][j][ip+2].d;
      N++;
    }
  }
  
  for (i=ip-2; i<=ip+2; i++) {
    for (k=kp-1; k<=kp+1; k++) {
      rho_amb += pG->U[k][j-2][i].d + pG->U[k][j+2][ip].d;
      N++;
    }
  }
  
  for (j=jp-2; j<=jp+2; j++) {
    for (i=ip-2; i<=ip+2; i++) {
      rho_amb += pG->U[k-2][j][i].d + pG->U[k+2][j][ip].d;
      N++;
    }
  }
  
  rho_amb /= (Real)N;
  
  RS = R0*pow(pStar->m,ONE_3RD)*pow(rho_amb,-TWO_3RDS);
  
  printf("First guess for star %d is R_S = %e, with rho_amb = %e and m = %e\n",
         pStar->id,RS,rho_amb,pStar->m);

  return RS;
}

#endif /* STAR_PARTICLE && SELF_GRAVITY */
