#include "../copyright.h"

/*=============================================================================
 * FILE: synchro_starparticles.c
 *
 * PURPOSE: keep a copy of all star particle information on each processor
 *
 * CONTAINS PUBLIC FUNCTIONS:
 *   synchro_starparticles();
 *
 * CONTAINS PRIVATE FUNCTIONS:
 *   merge_particles();
 *   BDghost_particle_global_fda();
 *
 * HISTORY:
 *   Written by Hao Gong, May 2011
 *   Modified by Aaron Skinner, Mar. 2013
 *============================================================================*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../defs.h"
#include "../athena.h"
#include "../prototypes.h"
#include "../globals.h"

#define OVERLAP
#ifdef STAR_PARTICLE

/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 *============================================================================*/
void merge_particles(GridS *pG, int nGStars, StarParS *pStarbuff);
void BDghost_particle_global_fda(DomainS *pD, StarParS *pStar);

/*=========================== PUBLIC FUNCTIONS ===============================*/
/*----------------------------------------------------------------------------*/

void synchro_starparticles(DomainS *pD)
{
  GridS *pG=pD->Grid;
  StarParListS *pLstars=NULL;  /* pointer to local star particle list */
  StarParListS *pGstars=NULL;  /* pointer to global star particle list */
  StarParListS *pGstars_fda=NULL;  /* pointer to list for density assignment */
  StarParS *pStar=NULL;  /* pointer to current working star particle */
  StarParS *sendbuf,*recvbuf;
  int i,nStars,starpar_id;
  
#ifdef MPI_PARALLEL
  int *recvbuf_iproc,*recvbuf_displs_iproc;
  int max_nStars;
  int m,gm;
  int nProc;
  int mpierr;
  
  MPI_Comm_size(MPI_COMM_WORLD, &nProc);
  /* find the maximum particle number on each processor */
  m = pG->nLstars;
  mpierr = MPI_Allreduce(&m, &gm, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  if (mpierr != MPI_SUCCESS)
    ath_error("[synchro_starparticles]: Error calling MPI_Allreduce!\n");
  max_nStars = gm;
  nStars = nProc*max_nStars;

  if (max_nStars > 0) {
    /* Create local buffer to store particle information for passing */
    sendbuf = (StarParS *)calloc_1d_array(max_nStars,sizeof(StarParS));
    if (sendbuf == NULL)
      ath_error("[synchro_starparticles]: Error calling calloc_1d_array!\n");
    recvbuf = (StarParS *)calloc_1d_array(nStars,sizeof(StarParS));
    if (recvbuf == NULL)
      ath_error("[synchro_starparticles]: Error calling calloc_1d_array!\n");

    /* Fill send buffer with local star particles. */
    /* NB:  Empty slots in send and receive buffers will have id=-1. */
    pLstars = pG->Lstars;
    for (i=0; i<max_nStars; i++) {
      if (pLstars) {
        sendbuf[i] = pLstars->starpar;
        pLstars = pLstars->next;
      }
      else sendbuf[i].id = -1;
    }

    /* Set up displs & recv_counts for passing array */
    for (i=0; i<nStars; i++) recvbuf[i].id = -1;
    recvbuf_iproc = (int *)calloc(nProc,sizeof(int));
    recvbuf_displs_iproc = (int *)calloc(nProc,sizeof(int));
    for (i=0; i<nProc; i++) {
      recvbuf_iproc[i] = max_nStars;
      recvbuf_displs_iproc[i] = i*max_nStars;
    }
        
    /* Communicate local star particle lists all-to-all */
    mpierr = MPI_Allgatherv(sendbuf, max_nStars, pD->StarParType, recvbuf,
                            recvbuf_iproc, recvbuf_displs_iproc,
                            pD->StarParType, MPI_COMM_WORLD);
    if (mpierr != MPI_SUCCESS)
      ath_error("[synchro_starparticles]: Error calling MPI_Allgatherv!\n");
  }
  
#else /* not MPI_PARALLEL */
  nStars = pG->nLstars;

  if (nStars > 0) {
    recvbuf = (StarParS *)calloc_1d_array(nStars,sizeof(StarParS));
    if (recvbuf == NULL)
      ath_error("[synchro_starparticles]: Error calling calloc_1d_array!\n");

    /* Fill send buffer with local star particles.  Empty slots have id = -1. */
    pLstars = pG->Lstars;
    for (i=0; i<nStars; i++) {
      if (pLstars) {
        recvbuf[i] = pLstars->starpar;
        pLstars = pLstars->next;
      }
      else recvbuf[i].id = -1;
    }
  }
  
#endif  /* MPI_PARALLEL */
  
  if (nStars > 0) {
    
    if (nStars > 1)
      merge_particles(pG,nStars,recvbuf);
    
    /* Destroy all current star particle lists */
    starpar_destruct_local(pG);
    starpar_destruct_global(pG);
    starpar_destruct_global_fda(pG);
    
    for (i=0; i<nStars; i++) {
      pStar = &(recvbuf[i]);
      if (pStar->id != -1) {
        /* Give newly created star particles a unique identifier */
        if (pStar->isnew==1) {
          pStar->id = pG->next_starpar_id++;
          pStar->isnew = 0;
        }
        
        /* Create global star particle list gathering from all local lists */
        pGstars = (StarParListS *)calloc_1d_array(1,sizeof(StarParListS));
        if (pGstars == NULL)
          ath_error("[synchro_starparticles]: Error calling calloc_1d_array!\n");
        pGstars->starpar = *pStar;
        starpar_push_global(pG,pGstars);
        
        /* Create star particle list for density assignment starting from 
         * global list */
        pGstars_fda = (StarParListS *)calloc_1d_array(1,sizeof(StarParListS));
        if (pGstars_fda == NULL)
          ath_error("[synchro_starparticles]: Error calling calloc_1d_array!\n");
        pGstars_fda->starpar = *pStar;
        starpar_push_global_fda(pG,pGstars_fda);
        
        /* Push ghost particles which are close to the border to the list for
         * density assignment if direction is periodic */
        BDghost_particle_global_fda(pD,pStar);
        
        /* Create local star particle list */
        if (starpar_ingrid(pD,pStar)) {
          pLstars = (StarParListS *)calloc_1d_array(1,sizeof(StarParListS));
          if (pLstars == NULL)
            ath_error("[synchro_starparticles]: Error calling calloc_1d_array!\n");
          pLstars->starpar = *pStar;
          starpar_push_local(pG,pLstars);
        }
      }
    }
    
#ifdef MPI_PARALLEL
    free_1d_array(sendbuf);
#endif /* MPI_PARALLEL */
    free_1d_array(recvbuf);
  }
  
  return;
}


/*=========================== PRIVATE FUNCTIONS ==============================*/
/*----------------------------------------------------------------------------*/

void merge_particles(GridS *pG, int nGStars, StarParS *pStarbuff) {
  int i,j;
  int ii,jj,kk,i1i,i2i,i3i,i1j,i2j,i3j;
  int starpar_id;
  Real distance;
  Real x1i,x2i,x3i,x1j,x2j,x3j;
  Real v1i,v2i,v3i,v1j,v2j,v3j;
  Real M1i,M2i,M3i,M1j,M2j,M3j;
  Real mi,mj;
  Real new_m,new_m_inv,new_x1,new_x2,new_x3,new_v1,new_v2,new_v3;
  Real r_control = ((double)NSINK_STARP+0.5)*MAX(pG->dx1,MAX(pG->dx2,pG->dx3));
  
  for (i=0; i<=nGStars-2; i++) {
    for (j=i+1; j<=nGStars-1; j++) {
      if ((pStarbuff[i].id != -1) && (pStarbuff[j].id != -1)) {
        x1i = pStarbuff[i].x1; x2i = pStarbuff[i].x2; x3i = pStarbuff[i].x3;
        x1j = pStarbuff[j].x1; x2j = pStarbuff[j].x2; x3j = pStarbuff[j].x3;
        
        distance = sqrt(SQR(x1i-x1j) + SQR(x2i-x2j) + SQR(x3i-x3j));
        if (distance <= 2.0*r_control) {
          mi  = pStarbuff[i].m;  mj  = pStarbuff[j].m;
          M1i = pStarbuff[i].v1*mi; M2i = pStarbuff[i].v2*mi; M3i = pStarbuff[i].v3*mi;
          M1j = pStarbuff[j].v1*mj; M2j = pStarbuff[j].v2*mj; M3j = pStarbuff[j].v3*mj;
/*          
          if (mi >= mj) {
            mj -= pStarbuff[j].mghost;
            M1j -= pStarbuff[j].M1ghost;
            M2j -= pStarbuff[j].M2ghost;
            M3j -= pStarbuff[j].M3ghost;
            if (mj < 0) ath_error("[syncrho_SP] mass cannot be negative! %g %g %g\n",mi,mj,pStarbuff[j].mghost);
          } else {
            mi -= pStarbuff[i].mghost;
            M1i -= pStarbuff[i].M1ghost;
            M2i -= pStarbuff[i].M2ghost;
            M3i -= pStarbuff[i].M3ghost;
            if (mi < 0) ath_error("[syncrho_SP] mass cannot be negative! %g %g %g\n",mi,mj,pStarbuff[j].mghost);
          }
*/

          new_m = mi + mj;
          new_m_inv = 1.0/new_m;
          new_v1 = (M1i + M1j)*new_m_inv;
          new_v2 = (M2i + M2j)*new_m_inv;
          new_v3 = (M3i + M3j)*new_m_inv;
          new_x1 = (x1i*mi + x1j*mj)*new_m_inv;
          new_x2 = (x2i*mi + x2j*mj)*new_m_inv;
          new_x3 = (x3i*mi + x3j*mj)*new_m_inv;
          
          if (mi >= mj) {
            if (myID_Comm_world==0)
              printf("Merging star %d into star %d\n",
                     pStarbuff[j].id,pStarbuff[i].id);
            pStarbuff[j].id = -1;
            
            pStarbuff[i].m = new_m;
            pStarbuff[i].x1 = new_x1;
            pStarbuff[i].x2 = new_x2;
            pStarbuff[i].x3 = new_x3;
            pStarbuff[i].v1 = new_v1;
            pStarbuff[i].v2 = new_v2;
            pStarbuff[i].v3 = new_v3;
            pStarbuff[i].merge_history++;
          } else {
            if (myID_Comm_world==0)
              printf("Merging star %d into star %d\n",
                     pStarbuff[i].id,pStarbuff[j].id);
            pStarbuff[i].id = -1;
            
            pStarbuff[j].m = new_m;
            pStarbuff[j].x1 = new_x1;
            pStarbuff[j].x2 = new_x2;
            pStarbuff[j].x3 = new_x3;
            pStarbuff[j].v1 = new_v1;
            pStarbuff[j].v2 = new_v2;
            pStarbuff[j].v3 = new_v3;
            pStarbuff[j].merge_history++;
          }
        }
      }
    }
  }
}

/*----------------------------------------------------------------------------*/

/* If the star particle is within NSINK_STARP zone of a periodic boundary, 
 * create a new star particle at the opposite boundary */
/* NB:  There will then be duplicate star particles on the list for density
 * assignment, but only the part of each star particle's shape function that
 * lies within the active Domain will be used.  */
void BDghost_particle_global_fda(DomainS *pD, StarParS *pStar)
{
  GridS *pG = pD->Grid;
  StarParListS *pListNew=NULL;
  Real Lx, Ly, Lz;

#ifdef SHEARING_BOX
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


  
  if ((pD->bcflag_ix1 == 4) && (pD->bcflag_ox1 == 4)) {
//    if ((pStar->x1 >= pD->RootMinX[0]) && (pStar->x1 < pD->RootMinX[0] + NSINK_STARP*pG->dx1)) {
//    if ((pStar->x1 >= pD->RootMinX[0]) && (pStar->x1 < pD->RootMinX[0] + NFEEDBACK*pG->dx1)) {
    if ((pStar->x1 >= pD->RootMinX[0]) && (pStar->x1 < pD->RootMaxX[0])) {
//      printf("Copying star %d periodically in -x direction...\n",pStar->id);
      pListNew = (StarParListS *)calloc_1d_array(1,sizeof(StarParListS));
      if (pListNew == NULL)
        ath_error("[BDghost_particle_global_fda]: Error calling calloc_1d_array!\n");
      pListNew->starpar = *pStar;
      pListNew->starpar.x1 += Lx;
#ifdef SHEARING_BOX
      pListNew->starpar.x2 -= deltay;
#ifndef FARGO
      pListNew->starpar.v2 -= qomL;
#endif
#endif
      starpar_push_global_fda(pG,pListNew);

/* additional copy due to shearing box */
#ifdef SHEARING_BOX
      pListNew = (StarParListS *)calloc_1d_array(1,sizeof(StarParListS));
      if (pListNew == NULL)
        ath_error("[BDghost_particle_global_fda]: Error calling calloc_1d_array!\n");
      pListNew->starpar = *pStar;
      pListNew->starpar.x1 += Lx;
      pListNew->starpar.x2 -= deltay + Ly;
#ifndef FARGO
      pListNew->starpar.v2 -= qomL;
#endif
      starpar_push_global_fda(pG,pListNew);
#endif
    }
    
//    if ((pStar->x1 >= pD->RootMaxX[0] - NSINK_STARP*pG->dx1) && (pStar->x1 < pD->RootMaxX[0])) {
//    if ((pStar->x1 >= pD->RootMaxX[0] - NFEEDBACK*pG->dx1) && (pStar->x1 < pD->RootMaxX[0])) {
    if ((pStar->x1 >= pD->RootMinX[0]) && (pStar->x1 < pD->RootMaxX[0])) {
//      printf("Copying star %d periodically in +x direction...\n",pStar->id);
      pListNew = (StarParListS *)calloc_1d_array(1,sizeof(StarParListS));
      if (pListNew == NULL)
        ath_error("[BDghost_particle_global_fda]: Error calling calloc_1d_array!\n");
      pListNew->starpar = *pStar;
      pListNew->starpar.x1 -= Lx;
#ifdef SHEARING_BOX
      pListNew->starpar.x2 += deltay;
#ifndef FARGO
      pListNew->starpar.v2 += qomL;
#endif
#endif
      starpar_push_global_fda(pG,pListNew);

/* additional copy due to shearing box */
#ifdef SHEARING_BOX
      pListNew = (StarParListS *)calloc_1d_array(1,sizeof(StarParListS));
      if (pListNew == NULL)
        ath_error("[BDghost_particle_global_fda]: Error calling calloc_1d_array!\n");
      pListNew->starpar = *pStar;
      pListNew->starpar.x1 -= Lx;
      pListNew->starpar.x2 += deltay - Ly;
#ifndef FARGO
      pListNew->starpar.v2 += qomL;
#endif
      starpar_push_global_fda(pG,pListNew);
#endif
    }
  }
  
  if ((pD->bcflag_ix2 == 4) && (pD->bcflag_ox2 == 4)) {
//    if ((pStar->x2 >= pD->RootMinX[1]) && (pStar->x2 < pD->RootMinX[1] + NSINK_STARP*pG->dx2)) {
//    if ((pStar->x2 >= pD->RootMinX[1]) && (pStar->x2 < pD->RootMinX[1] + NFEEDBACK*pG->dx2)) {
    if ((pStar->x2 >= pD->RootMinX[1]) && (pStar->x2 < pD->RootMaxX[1])) {
      pListNew = (StarParListS *)calloc_1d_array(1,sizeof(StarParListS));
      if (pListNew == NULL)
        ath_error("[BDghost_particle_global_fda]: Error calling calloc_1d_array!\n");
      pListNew->starpar = *pStar;
      pListNew->starpar.x2 += Ly;
      starpar_push_global_fda(pG,pListNew);
    }
    
//    if ((pStar->x2 >= pD->RootMaxX[1] - NSINK_STARP*pG->dx2) && (pStar->x2 < pD->RootMaxX[1])) {
//    if ((pStar->x2 >= pD->RootMaxX[1] - NFEEDBACK*pG->dx2) && (pStar->x2 < pD->RootMaxX[1])) {
    if ((pStar->x2 >= pD->RootMinX[1]) && (pStar->x2 < pD->RootMaxX[1])) {
      pListNew = (StarParListS *)calloc_1d_array(1,sizeof(StarParListS));
      if (pListNew == NULL)
        ath_error("[BDghost_particle_global_fda]: Error calling calloc_1d_array!\n");
      pListNew->starpar = *pStar;
      pListNew->starpar.x2 -= Ly;
      starpar_push_global_fda(pG,pListNew);
    }
  }
  
#ifndef SELF_GRAVITY_USING_FFT_DISK
  if ((pD->bcflag_ix3 == 4) && (pD->bcflag_ox3 == 4)) {
    if ((pStar->x3 >= pD->RootMinX[2]) && (pStar->x3 < pD->RootMinX[2] + NSINK_STARP*pG->dx3)) {
      pListNew = (StarParListS *)calloc_1d_array(1,sizeof(StarParListS));
      if (pListNew == NULL)
        ath_error("[BDghost_particle_global_fda]: Error calling calloc_1d_array!\n");
      pListNew->starpar = *pStar;
      pListNew->starpar.x3 += Lz;
      starpar_push_global_fda(pG,pListNew);
    }
    
    if ((pStar->x3 >= pD->RootMaxX[2] - NSINK_STARP*pG->dx3) && (pStar->x3 < pD->RootMaxX[2])) {
      pListNew = (StarParListS *)calloc_1d_array(1,sizeof(StarParListS));
      if (pListNew == NULL)
        ath_error("[BDghost_particle_global_fda]: Error calling calloc_1d_array!\n");
      pListNew->starpar = *pStar;
      pListNew->starpar.x3 -= Lz;
      starpar_push_global_fda(pG,pListNew);
    }
  }
#endif

  return;
}

#endif /* STAR_PARTICLE */
