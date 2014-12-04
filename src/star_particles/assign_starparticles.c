#include "../copyright.h"
/*=============================================================================
 * FILE: assign_starparticles.c
 *
 * PURPOSE: Smooth each star particle's mass to surrounding cells so that the
 *   potential from gas+particle can be calculated by FFT.  
 *   Triangular-shaped-cloud (TSC) (see Hockney & Eastwood, 1981) scheme is 
 *   implemented.
 *
 * CONTAINS PUBLIC FUNCTIONS:
 *   assign_starparticles_3d();
 *
 * HISTORY:
 *   Written by Hao Gong, July 2011
 *   Modified by Aaron Skinner, Apr. 2013
 *============================================================================*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "../defs.h"
#include "../athena.h"
#include "../globals.h"
#include "../prototypes.h"

#if defined(STAR_PARTICLE) && defined(SELF_GRAVITY) && defined(FFT_ENABLED)

#ifdef STATIC_MESH_REFINEMENT
#error self gravity with FFT not yet implemented to work with SMR
#endif

/*=========================== PUBLIC FUNCTIONS ===============================*/
/*----------------------------------------------------------------------------*/

void assign_starparticles_3d(DomainS *pD, ath_fft_data *work)
{
  GridS *pG = pD->Grid;
  StarParListS *pList=NULL;
  StarParS *pStar=NULL;
  int i, is = pG->is, ie = pG->ie;
  int j, js = pG->js, je = pG->je;
  int k, ks = pG->ks, ke = pG->ke;
  int ip,jp,kp;
  Real x1,x2,x3;
  Real tmp,W1[9],W2[9],W3[9],d;

#if defined(SHEARING_BOX) && !defined(SELF_GRAVITY_USING_FFT_OBC)
  ath_error("[assign_SP] This shouldn't be called with Shearing Box\n");
#endif
  /* Interpolate star particle density onto the mesh using CIC or TSC shape
   * functions as described in Hockney & Eastwood. */
  pList = pG->Gstars_fda;
  while (pList) {
    pStar = &(pList->starpar);

    /* If the star particle is within NSINK_STARP+1 zones of an outflow 
     * (bcflag==2) boundary, then do not assign its density to the grid.
     * There is no problem if near a periodic (bcflag==4) boundary.  */
    if (((pStar->x1 < pD->MinX[0] + (NSINK_STARP+1)*pG->dx1) && (pD->bcflag_ix1 == 2)) ||
        ((pStar->x1 > pD->MaxX[0] - (NSINK_STARP+1)*pG->dx1) && (pD->bcflag_ox1 == 2)) ||
        ((pStar->x2 < pD->MinX[1] + (NSINK_STARP+1)*pG->dx2) && (pD->bcflag_ix2 == 2)) ||
        ((pStar->x2 > pD->MaxX[1] - (NSINK_STARP+1)*pG->dx2) && (pD->bcflag_ox2 == 2)) ||
        ((pStar->x3 < pD->MinX[2] + (NSINK_STARP+1)*pG->dx3) && (pD->bcflag_ix3 == 2)) ||
        ((pStar->x3 > pD->MaxX[2] - (NSINK_STARP+1)*pG->dx3) && (pD->bcflag_ox3 == 2))) {
      /* Do nothing! */
    }
    else {
      cc_ijk(pG,pStar->x1,pStar->x2,pStar->x3,&ip,&jp,&kp);
      cc_pos(pG,ip,jp,kp,&x1,&x2,&x3);
      
      /* Compute TSC weights.   */
      tmp = 0.5 + (pStar->x1 - x1)/pG->dx1;
      W1[2] = 0.5*SQR(tmp);  W1[0] = W1[2]-tmp+0.5;  W1[1] = 1.0-W1[0]-W1[2];
      tmp = 0.5 + (pStar->x2 - x2)/pG->dx2;
      W2[2] = 0.5*SQR(tmp);  W2[0] = W2[2]-tmp+0.5;  W2[1] = 1.0-W2[0]-W2[2];
      tmp = 0.5 + (pStar->x3 - x3)/pG->dx3;
      W3[2] = 0.5*SQR(tmp);  W3[0] = W3[2]-tmp+0.5;  W3[1] = 1.0-W3[0]-W3[2];
      //    printf("proc %d has W1[0]=%e, W1[1]=%e, W1[2]=%e\n",myID_Comm_world,W1[0],W1[1],W1[2]);
      //    printf("proc %d has W2[0]=%e, W2[1]=%e, W2[2]=%e\n",myID_Comm_world,W2[0],W2[1],W2[2]);
      //    printf("proc %d has W3[0]=%e, W3[1]=%e, W3[2]=%e\n",myID_Comm_world,W3[0],W3[1],W3[2]);
      
      /* Add weighted particle density to work array 0. */
      d = pStar->m/(pG->dx1*pG->dx2*pG->dx3);
      for (k=0; k<=2; k++) {
        for (j=0; j<=2; j++) {
          for (i=0; i<=2; i++) {
            if ((ip-is+i-1 >= 0) && (ip-is+i-1 <= pG->Nx[0] - 1) &&
                (jp-js+j-1 >= 0) && (jp-js+j-1 <= pG->Nx[1] - 1) &&
                (kp-ks+k-1 >= 0) && (kp-ks+k-1 <= pG->Nx[2] - 1)) {
              if ((abs(i-1) <= NSINK_STARP) ||
                  (abs(j-1) <= NSINK_STARP) ||
                  (abs(k-1) <= NSINK_STARP)) {
                work[F3DI(ip-is+i-1,jp-js+j-1,kp-ks+k-1,pG->Nx[0],pG->Nx[1],pG->Nx[2])][0] = d*W1[i]*W2[j]*W3[k];
              }
              else {
		ath_error("[assign_starparticles]: This should not be the case!\n");
                work[F3DI(ip-is+i-1,jp-js+j-1,kp-ks+k-1,pG->Nx[0],pG->Nx[1],pG->Nx[2])][0] += d*W1[i]*W2[j]*W3[k];
              }
            }
          }
        }
      }
    }
    
    pList = pList->next;
  }

  return;
}

#ifdef SHEARING_BOX
void assign_starparticles_3d_shear(DomainS *pD, Real ***RollDen)
{
  GridS *pG = pD->Grid;
  StarParListS *pList=NULL;
  StarParS *pStar=NULL;
  int i, is = pG->is, ie = pG->ie;
  int j, js = pG->js, je = pG->je;
  int k, ks = pG->ks, ke = pG->ke;
  int ip,jp,kp;
  Real x1,x2,x3;
  Real tmp,W1[9],W2[9],W3[9],d;

  /* Interpolate star particle density onto the mesh using CIC or TSC shape
   * functions as described in Hockney & Eastwood. */
  pList = pG->Gstars_fda;
  while (pList) {
    pStar = &(pList->starpar);

    /* If the star particle is within NSINK_STARP+1 zones of an outflow 
     * (bcflag==2) boundary, then do not assign its density to the grid.
     * There is no problem if near a periodic (bcflag==4) boundary.  */
    if (((pStar->x1 < pD->MinX[0] + (NSINK_STARP+1)*pG->dx1) && (pD->bcflag_ix1 == 2)) ||
        ((pStar->x1 > pD->MaxX[0] - (NSINK_STARP+1)*pG->dx1) && (pD->bcflag_ox1 == 2)) ||
        ((pStar->x2 < pD->MinX[1] + (NSINK_STARP+1)*pG->dx2) && (pD->bcflag_ix2 == 2)) ||
        ((pStar->x2 > pD->MaxX[1] - (NSINK_STARP+1)*pG->dx2) && (pD->bcflag_ox2 == 2)) ||
        ((pStar->x3 < pD->MinX[2] + (NSINK_STARP+1)*pG->dx3) && (pD->bcflag_ix3 == 2)) ||
        ((pStar->x3 > pD->MaxX[2] - (NSINK_STARP+1)*pG->dx3) && (pD->bcflag_ox3 == 2))) {
      /* Do nothing! */
    }
    else {
      cc_ijk(pG,pStar->x1,pStar->x2,pStar->x3,&ip,&jp,&kp);
      cc_pos(pG,ip,jp,kp,&x1,&x2,&x3);
      
      /* Compute TSC weights.   */
      tmp = 0.5 + (pStar->x1 - x1)/pG->dx1;
      W1[2] = 0.5*SQR(tmp);  W1[0] = W1[2]-tmp+0.5;  W1[1] = 1.0-W1[0]-W1[2];
      tmp = 0.5 + (pStar->x2 - x2)/pG->dx2;
      W2[2] = 0.5*SQR(tmp);  W2[0] = W2[2]-tmp+0.5;  W2[1] = 1.0-W2[0]-W2[2];
      tmp = 0.5 + (pStar->x3 - x3)/pG->dx3;
      W3[2] = 0.5*SQR(tmp);  W3[0] = W3[2]-tmp+0.5;  W3[1] = 1.0-W3[0]-W3[2];
      
      /* Add weighted particle density to work array 0. */
      d = pStar->m/(pG->dx1*pG->dx2*pG->dx3);
      for (k=0; k<=2; k++) {
        for (i=0; i<=2; i++) {
          for (j=0; j<=2; j++) {
            if ((ip+i-1 >= is) && (ip+i-1 <= ie) &&
                (jp+j-1 >= is) && (jp+j-1 <= je) &&
                (kp+k-1 >= is) && (kp+k-1 <= ke)) {
                RollDen[kp+k-1][ip+i-1][jp+j-1] = d*W1[i]*W2[j]*W3[k];
            }
          }
        }
      }
    }
    
    pList = pList->next;
  }

  return;
}
#endif /* SHEARING_BOX */
#endif /* STAR_PARTICLE && SELF_GRAVITY && FFT_ENABLED */
