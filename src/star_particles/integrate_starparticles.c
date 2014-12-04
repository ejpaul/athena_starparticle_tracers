#include "../copyright.h"

/*=============================================================================
 * FILE: integrate_starparticles.c
 *
 * PURPOSE: 2nd-order explicit particle integrator
 *
 * CONTAINS PUBLIC FUNCTIONS:
 *   integrate_starparticles();
 *   kick_correct_starparticles();
 *   kick_starparticles();
 *   drift_starparticles();
 *
 * CONTAINS PRIVATE FUNCTIONS:
 *   Get_Force();
 *
 * HISTORY:
 *   Written by Aaron Skinner, Nov. 2010, following Xuening Bai, Mar. 2009
 *   Modified by Hao Gong, Apr. 2011
 *   Modified by Aaron Skinner, Apr. 2013
 *   Modified by Aaron Skinner, Jul. 2013
 *============================================================================*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../defs.h"
#include "../athena.h"
#include "../prototypes.h"
#include "../globals.h"

#ifdef STAR_PARTICLE

/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 *   Get_Force()      - calculate forces other than the drag
 *============================================================================*/
void Get_Force(const GridS *pG, const StarParS *pStar, Real *f1, Real *f2,
               Real *f3);

/*=========================== PUBLIC FUNCTIONS ===============================*/

void integrate_starparticles(DomainS *pD)
{
  GridS *pG=pD->Grid;
  StarParListS *pList=NULL;  /* pointer to the current working particle */
  StarParListS *pListprev=NULL;  /* pointer to the previous particle */
  StarParS *pStar=NULL;
  Real f1,f2,f3;  /* force components */
  Real hdt_old=0.5*pG->dt_old, hdt=0.5*pG->dt;
  Real m,M1,M2,M3,dVol=pG->dx1*pG->dx2*pG->dx3;
  Real mold,M1old,M2old,M3old;
  Real dm,dM1,dM2,dM3;
  int ip,jp,kp,i,j,k;
  int remove_flag;
  Real Lx, Ly, Lz;
#ifdef SHEARING_BOX
  Real Omdt = Omega_0*pG->dt, Py;
  Real qomL, yshear, deltay;
#endif

  Lx = pD->RootMaxX[0] - pD->RootMinX[0];
  Ly = pD->RootMaxX[1] - pD->RootMinX[1];
  Lz = pD->RootMaxX[2] - pD->RootMinX[2];

#ifdef SHEARING_BOX
  qomL = qshear*Omega_0*Lx;
  yshear = qomL*pG->time;
  deltay = fmod(yshear, Ly);
#endif

  pList = pG->Lstars;

  /* loop over all radiation particles */
  while (pList) {
    pStar = &(pList->starpar);
    
    /* Step 1:  Calculate the (specific) force at current position */
    /* If the star particle is near an outflow boundary, reset the force to 0
     * so that it flies off the Grid.  Otherwise, weird things can happen with
     * the extrapolation of the self-grativation potential using OBCs. */
    if (((pStar->x1 < pD->MinX[0] + (NSINK_STARP+1)*pG->dx1) && (pD->bcflag_ix1 == 2)) ||
        ((pStar->x1 > pD->MaxX[0] - (NSINK_STARP+1)*pG->dx1) && (pD->bcflag_ox1 == 2)) ||
        ((pStar->x2 < pD->MinX[1] + (NSINK_STARP+1)*pG->dx2) && (pD->bcflag_ix2 == 2)) ||
        ((pStar->x2 > pD->MaxX[1] - (NSINK_STARP+1)*pG->dx2) && (pD->bcflag_ox2 == 2)) ||
        ((pStar->x3 < pD->MinX[2] + (NSINK_STARP+1)*pG->dx3) && (pD->bcflag_ix3 == 2)) ||
        ((pStar->x3 > pD->MaxX[2] - (NSINK_STARP+1)*pG->dx3) && (pD->bcflag_ox3 == 2))) {
      f1 = 0.0;
      f2 = 0.0;
      f3 = 0.0;
    }
    else {
      Get_Force(pG, pStar, &f1, &f2, &f3);
    }

    cc_ijk(pG,pStar->x1,pStar->x2,pStar->x3,&ip,&jp,&kp);
    mold=0.0; M1old=0.0; M2old=0.0; M3old=0.0;
    for (k=kp-NSINK_STARP; k<=kp+NSINK_STARP; k++) {
      for (j=jp-NSINK_STARP; j<=jp+NSINK_STARP; j++) {
        for (i=ip-NSINK_STARP; i<=ip+NSINK_STARP; i++) {
          mold  += pG->U[k][j][i].d;
          M1old += pG->U[k][j][i].M1;
          M2old += pG->U[k][j][i].M2;
          M3old += pG->U[k][j][i].M3;
        }
      }
    }
    pStar->mghost = mold;
    pStar->M1ghost = M1old;
    pStar->M2ghost = M2old;
    pStar->M3ghost = M3old;
    pStar->x1_old = pStar->x1;
    pStar->x2_old = pStar->x2;
    pStar->x3_old = pStar->x3;
 
    /* Step 2:  (KICK) Update velocity for a half time step */
    /* v^{n-1/2} -> v^n */
    if (pG->time > 0.0) {
      pStar->v1 += hdt_old*f1;
      pStar->v2 += hdt_old*f2;
      pStar->v3 += hdt_old*f3;
    }
   
    /* Step 3:  (KICK) Update velocity for a half time step */
    /* v^n -> v^{n+1/2} */
    pStar->v1 += hdt*f1;
    pStar->v2 += hdt*f2;
    pStar->v3 += hdt*f3;

#ifdef SHEARING_BOX
/* Shearing box source temrs in Quinn et al. (2010) */
/* The order does matter */
    Py = pStar->v2 + 2.0*Omega_0*pStar->x1;
    pStar->v1 += (qshear-2.0)*Omdt*Omega_0*pStar->x1 + Omdt*Py;
    pStar->v2 -= Omdt*pStar->v1;
#endif
   
    /* Step 4:  (DRIFT) Update position for a full time step */
    /* x^n -> x^{n+1/2} */
    if (pG->Nx[0] > 1) pStar->x1 += pG->dt*pStar->v1;
    if (pG->Nx[1] > 1) pStar->x2 += pG->dt*pStar->v2;
    if (pG->Nx[2] > 1) pStar->x3 += pG->dt*pStar->v3;

#ifdef SHEARING_BOX
    /* Step 5:  (KICK) Update velocity for a half time step due to Shear */
    /* v^{n+1/2} -> v^{n+1} */
    pStar->v1 += (qshear-2.0)*Omdt*Omega_0*pStar->x1 + Omdt*Py;
    pStar->v2 = Py - 2.0*Omega_0*pStar->x1;
#endif

/*
    cc_ijk(pG,pStar->x1,pStar->x2,pStar->x3,&ip,&jp,&kp);
    mold=0.0; M1old=0.0; M2old=0.0; M3old=0.0;
    for (k=kp-NSINK_STARP; k<=kp+NSINK_STARP; k++) {
      for (j=jp-NSINK_STARP; j<=jp+NSINK_STARP; j++) {
        for (i=ip-NSINK_STARP; i<=ip+NSINK_STARP; i++) {
          mold  += pG->U[k][j][i].d;
          M1old += pG->U[k][j][i].M1;
          M2old += pG->U[k][j][i].M2;
          M3old += pG->U[k][j][i].M3;
        }
      }
    }
    pStar->mghost -= mold;
    pStar->M1ghost -= M1old;
    pStar->M2ghost -= M2old;
    pStar->M3ghost -= M3old;
*/
 
    remove_flag = 0;

    if (pStar->x1 <  pD->RootMinX[0]) {
      if ((pD->bcflag_ix1 == 4) && (pD->bcflag_ox1 == 4)){
        pStar->x1 += Lx;
#ifdef SHEARING_BOX
        pStar->x2 -= deltay;
#ifndef FARGO
        pStar->v2 -= qomL;
#endif
#endif
      } else remove_flag++;
    }
    if (pStar->x1 > pD->RootMaxX[0]) {
      if ((pD->bcflag_ix1 == 4) && (pD->bcflag_ox1 == 4)){
        pStar->x1 -= Lx;
#ifdef SHEARING_BOX
        pStar->x2 += deltay;
#ifndef FARGO
        pStar->v2 += qomL;
#endif
#endif
      } else remove_flag++;
    }
    if (pStar->x2 <  pD->RootMinX[1]) {
      if ((pD->bcflag_ix2 == 4) && (pD->bcflag_ox2 == 4))
        pStar->x2 += Ly;
      else remove_flag++;
    }
    if (pStar->x2 > pD->RootMaxX[1]) {
      if ((pD->bcflag_ix2 == 4) && (pD->bcflag_ox2 == 4))
        pStar->x2 -= Ly;
      else remove_flag++;
    }
    if (pStar->x3 <  pD->RootMinX[2]) {
      if ((pD->bcflag_ix3 == 4) && (pD->bcflag_ox3 == 4))
        pStar->x3 += Lz;
      else remove_flag++;
    }
    if (pStar->x3 > pD->RootMaxX[2]) {
      if ((pD->bcflag_ix3 == 4) && (pD->bcflag_ox3 == 4))
        pStar->x3 -= Lz;
      else remove_flag++;
    }

    if (remove_flag) {
      printf("Deleting star %d\n", pStar->id);
      pStar->id=-1;
      if (pListprev == NULL) {
        /* the particle is the head of the linklist to be removed */
        pG->Lstars = pList->next;
        free_1d_array(pList);
        pList = pG->Lstars;
        pListprev = NULL;
      }
      else {
        /* the particle is in the middle of the linklist to be removed */
        pListprev->next = pList->next;
        free_1d_array(pList);
        pList = pListprev->next;
      }
      /* modify the corresponding number of star particles */
/*
      pG->nLstars--;
      pG->nGstars--;
*/
    }
    else {
      /* the particle is not to be removed */

      pListprev = pList;
      pList = pList->next;
    }
    
  } /* End of the loop */

  return;
}


/*=========================== PRIVATE FUNCTIONS ==============================*/
/*--------------------------------------------------------------------------- */
/* Interpolate particle forces from the mesh using TSC shape functions as
 *   described in Hockney & Eastwood.  NOTE:  There are more efficient ways
 *   to implement this, but since the radiation particle lists will typically
 *   be small, there is no need.
 * Input:
 *   pG: grid;
 *   x1,x2,x3: particle position;
 * Output:
 *   f1,f2,f3: force components;  */
void Get_Force(const GridS *pG, const StarParS *pStar, 
               Real *f1, Real *f2, Real *f3)
{
  int ip,jp,kp;
  Real x1,x2,x3,tmp,W[9],phi[4];
  Real dx1i=1.0/pG->dx1, dx2i=1.0/pG->dx2, dx3i=1.0/pG->dx3;

  /* initialize forces to zero */
  *f1 = 0.0;  *f2 = 0.0;  *f3 = 0.0;

  /* calculate indices and spatial coordinates of the lower nearest neighbor */
  cc_ijk(pG,pStar->x1,pStar->x2,pStar->x3,&ip,&jp,&kp);
  cc_pos(pG,ip,jp,kp,&x1,&x2,&x3);

  if (pG->Nx[0] > 1) {
    /* Compute TSC weights.   */
    tmp = 0.5 + (pStar->x1 - x1)*dx1i;
    W[2] = 0.5*SQR(tmp);  W[0] = W[2] - tmp + 0.5;  W[1] = 1.0 - W[0] - W[2];
    if (StaticGravPot != NULL){
      phi[0] = (*StaticGravPot)((x1-1.5*pG->dx1),x2,x3);
      phi[1] = (*StaticGravPot)((x1-0.5*pG->dx1),x2,x3);
      phi[2] = (*StaticGravPot)((x1+0.5*pG->dx1),x2,x3);
      phi[3] = (*StaticGravPot)((x1+1.5*pG->dx1),x2,x3);
      *f1 -= ((phi[1]-phi[0])*W[0] +
              (phi[2]-phi[1])*W[1] +
              (phi[3]-phi[2])*W[2])*dx1i;
    }
#ifdef SELF_GRAVITY
    /* compute (specific) force components from self gravity */
    *f1 -= 0.5*((pG->Phi[kp][jp][ip  ] - pG->Phi[kp][jp][ip-2])*W[0]
              + (pG->Phi[kp][jp][ip+1] - pG->Phi[kp][jp][ip-1])*W[1]
              + (pG->Phi[kp][jp][ip+2] - pG->Phi[kp][jp][ip  ])*W[2])*dx1i;
#endif /* SELF_GRAVITY */
  }

  if (pG->Nx[1] > 1) {
    /* Compute TSC weights.   */
    tmp = 0.5 + (pStar->x2 - x2)*dx2i;
    W[2] = 0.5*SQR(tmp);  W[0] = W[2] - tmp + 0.5;  W[1] = 1.0 - W[0] - W[2];
    if (StaticGravPot != NULL){
      phi[0] = (*StaticGravPot)(x1,(x2-1.5*pG->dx2),x3);
      phi[1] = (*StaticGravPot)(x1,(x2-0.5*pG->dx2),x3);
      phi[2] = (*StaticGravPot)(x1,(x2+0.5*pG->dx2),x3);
      phi[3] = (*StaticGravPot)(x1,(x2+1.5*pG->dx2),x3);
      *f2 -= ((phi[1]-phi[0])*W[0] +
              (phi[2]-phi[1])*W[1] +
              (phi[3]-phi[2])*W[2])*dx2i;
    }
#ifdef SELF_GRAVITY
    /* compute (specific) force components from self gravity */
    *f2 -= 0.5*((pG->Phi[kp][jp  ][ip] - pG->Phi[kp][jp-2][ip])*W[0]
              + (pG->Phi[kp][jp+1][ip] - pG->Phi[kp][jp-1][ip])*W[1]
              + (pG->Phi[kp][jp+2][ip] - pG->Phi[kp][jp  ][ip])*W[2])*dx2i;
#endif /* SELF_GRAVITY */
  }

  if (pG->Nx[2] > 1) {
    /* Compute TSC weights.   */
    tmp = 0.5 + (pStar->x3 - x3)*dx3i;
    W[2] = 0.5*SQR(tmp);  W[0] = W[2] - tmp + 0.5;  W[1] = 1.0 - W[0] - W[2];
    if (StaticGravPot != NULL){
      phi[0] = (*StaticGravPot)(x1,x2,(x3-1.5*pG->dx3));
      phi[1] = (*StaticGravPot)(x1,x2,(x3-0.5*pG->dx3));
      phi[2] = (*StaticGravPot)(x1,x2,(x3+0.5*pG->dx3));
      phi[3] = (*StaticGravPot)(x1,x2,(x3+1.5*pG->dx3));
      *f3 -= ((phi[1]-phi[0])*W[0] +
              (phi[2]-phi[1])*W[1] +
              (phi[3]-phi[2])*W[2])*dx3i;
    }
#ifdef SELF_GRAVITY
    /* compute (specific) force components from self gravity */
    *f3 -= 0.5*((pG->Phi[kp  ][jp][ip] - pG->Phi[kp-2][jp][ip])*W[0]
              + (pG->Phi[kp+1][jp][ip] - pG->Phi[kp-1][jp][ip])*W[1]
              + (pG->Phi[kp+2][jp][ip] - pG->Phi[kp  ][jp][ip])*W[2])*dx3i;
#endif /* SELF_GRAVITY */
  }

  return;
}

#endif /* STAR_PARTICLE */
