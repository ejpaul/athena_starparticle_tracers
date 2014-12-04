#include "../copyright.h"

/*=============================================================================
 * FILE: modify_ghost_region_starparticles.c
 *
 * PURPOSE: Modify gas density and momentum (+ internal energy) inside the star particle's control 
 *          volume.   Density and momentum (+ internal energy) are reset using the average values 
 *          from the surrounding active zones.
 *
 * CONTAINS PUBLIC FUNCTIONS:
 *   modify_ghost_region_starparticles();
 *
 * HISTORY:
 *   Written by Hao Gong, Aug. 2011
 *   Modified by Aaron Skinner, Apr. 2013
 *   Modified by Chang-Goo Kim, Jan. 2014 for non-isothermal EOS
 *============================================================================*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../athena.h"
#include "../defs.h"
#include "../prototypes.h"
#include "../globals.h"


#ifdef STAR_PARTICLE

#ifdef MHD
#error This is not working with MHD
#endif

#define PRIMITIVE_AVG
#define FACE_ONLY

/*=========================== PUBLIC FUNCTIONS ===============================*/
/*----------------------------------------------------------------------------*/

void modify_ghost_region_starparticles(DomainS *pD) {
  GridS *pG = pD->Grid;
  StarParListS *pList=NULL;
  StarParS *pStar=NULL;
  int is=pG->is, ie=pG->ie, js=pG->js, je=pG->je, ks=pG->ks, ke=pG->ke;
  int i,j,k;
  int ip,jp,kp;
  int ip_old,jp_old,kp_old;
  Real m,M1,M2,M3;


  pList = pG->Lstars;
  while (pList) {
    pStar = &(pList->starpar);
    cc_ijk(pG,pStar->x1_old,pStar->x2_old,pStar->x3_old,&ip_old,&jp_old,&kp_old);
    cc_ijk(pG,pStar->x1,pStar->x2,pStar->x3,&ip,&jp,&kp);
  
    if ((ip != ip_old) || (jp != jp_old) || (kp != kp_old)) {
      m=0.0; M1=0.0; M2=0.0; M3=0.0;
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
      pStar->mghost -= m;
      pStar->M1ghost -= M1;
      pStar->M2ghost -= M2;
      pStar->M3ghost -= M3;   
    }
    pList = pList->next;
  }

  pList = pG->Gstars_fda;
  while (pList) {
    pStar = &(pList->starpar);
    cc_ijk(pG,pStar->x1,pStar->x2,pStar->x3,&ip,&jp,&kp);

    if ((ip >= is-1) && (ip <= ie+1) &&
        (jp >= js-1) && (jp <= je+1) &&
        (kp >= ks-1) && (kp <= ke+1)) {
      set_ghost_region(pG,pStar);  
    }
    pList = pList->next;
  }
}

#ifdef FACE_ONLY
void set_ghost_region(GridS *pG,StarParS *pStar){
  int i, j, k;
  int ip, jp, kp;
  Real fact,w1,w2,w3,davg,M1avg,M2avg,M3avg,Pavg,Eavg;
  ConsS U[2][2][2];
#ifdef PRIMITIVE_AVG
  PrimS W[2][2][2];
#endif

  cc_ijk(pG,pStar->x1,pStar->x2,pStar->x3,&ip,&jp,&kp);

  /* 6 front faces */
  pG->U[kp][jp][ip+1] = pG->U[kp][jp][ip+2];
  pG->U[kp][jp][ip-1] = pG->U[kp][jp][ip-2];
      
  pG->U[kp][jp+1][ip] = pG->U[kp][jp+2][ip];
  pG->U[kp][jp-1][ip] = pG->U[kp][jp-2][ip];
     
  pG->U[kp+1][jp][ip] = pG->U[kp+2][jp][ip];
  pG->U[kp-1][jp][ip] = pG->U[kp-2][jp][ip];
      
  fact=ONE_3RD;
  for(k=-1;k<=1;k+=2){
  for(j=-1;j<=1;j+=2){
  for(i=-1;i<=1;i+=2){
    U[0][0][1] = pG->U[kp+k  ][jp+j  ][ip+i*2];
    U[0][1][0] = pG->U[kp+k  ][jp+j*2][ip+i  ];
    U[1][0][0] = pG->U[kp+k*2][jp+j  ][ip+i  ];

    davg  = fact*(U[0][0][1].d +U[0][1][0].d +U[1][0][0].d);
    M1avg = fact*(U[0][0][1].M1+U[0][1][0].M1+U[1][0][0].M1);
    M2avg = fact*(U[0][0][1].M2+U[0][1][0].M2+U[1][0][0].M2);
    M3avg = fact*(U[0][0][1].M3+U[0][1][0].M3+U[1][0][0].M3);

    pG->U[kp+k][jp+j][ip+i].d = davg;
    pG->U[kp+k][jp+j][ip+i].M1=M1avg;
    pG->U[kp+k][jp+j][ip+i].M2=M2avg;
    pG->U[kp+k][jp+j][ip+i].M3=M3avg;
#ifndef BAROTROPIC
#ifdef PRIMITIVE_AVG
    W[0][0][1] = Cons_to_Prim(&(U[0][0][1]));
    W[0][1][0] = Cons_to_Prim(&(U[0][1][0]));
    W[1][0][0] = Cons_to_Prim(&(U[1][0][0]));

    Pavg = fact*(W[0][0][1].P+W[0][1][0].P+W[1][0][0].P);

    Eavg = Pavg/Gamma_1 + 0.5*(SQR(M1avg)+SQR(M2avg)+SQR(M3avg))/davg;
#else /* PRIMITIVE_AVG */
    Eavg  = fact*(U[0][0][1].E+U[0][1][0].E+U[1][0][0].E);
#endif /* PRIMITIVE_AVG */
    pG->U[kp+k][jp+j][ip+i].E = Eavg;
#endif /* BAROTROPIC */
  }}}
      
/* 4 sides for 3 middle-slices */

  fact=0.5;

/* ip slice */
  for(k=-1;k<=1;k+=2){
  for(j=-1;j<=1;j+=2){
    U[0][1][0] = pG->U[kp+k  ][jp+j*2][ip];
    U[1][0][0] = pG->U[kp+k*2][jp+j  ][ip];

    davg  = fact*(U[0][1][0].d +U[1][0][0].d );
    M1avg = fact*(U[0][1][0].M1+U[1][0][0].M1);
    M2avg = fact*(U[0][1][0].M2+U[1][0][0].M2);
    M3avg = fact*(U[0][1][0].M3+U[1][0][0].M3);

    pG->U[kp+k][jp+j][ip].d = davg;
    pG->U[kp+k][jp+j][ip].M1=M1avg;
    pG->U[kp+k][jp+j][ip].M2=M2avg;
    pG->U[kp+k][jp+j][ip].M3=M3avg;
#ifndef BAROTROPIC
#ifdef PRIMITIVE_AVG
    W[0][1][0] = Cons_to_Prim(&(U[0][1][0]));
    W[1][0][0] = Cons_to_Prim(&(U[1][0][0]));

    Pavg = fact*(W[0][1][0].P+W[1][0][0].P);

    Eavg = Pavg/Gamma_1 + 0.5*(SQR(M1avg)+SQR(M2avg)+SQR(M3avg))/davg;
#else /* PRIMITIVE_AVG */
    Eavg = fact*(U[0][1][0].E+U[1][0][0].E);
#endif /* PRIMITIVE_AVG */
    pG->U[kp+k][jp+j][ip].E = Eavg;
#endif /* BAROTROPIC */
  }}
 
/* jp slice */
  for(k=-1;k<=1;k+=2){
  for(i=-1;i<=1;i+=2){
    U[0][1][0] = pG->U[kp+k  ][jp][ip+i*2];
    U[1][0][0] = pG->U[kp+k*2][jp][ip+i  ];

    davg  = fact*(U[0][1][0].d +U[1][0][0].d );
    M1avg = fact*(U[0][1][0].M1+U[1][0][0].M1);
    M2avg = fact*(U[0][1][0].M2+U[1][0][0].M2);
    M3avg = fact*(U[0][1][0].M3+U[1][0][0].M3);

    pG->U[kp+k][jp][ip+i].d = davg;
    pG->U[kp+k][jp][ip+i].M1=M1avg;
    pG->U[kp+k][jp][ip+i].M2=M2avg;
    pG->U[kp+k][jp][ip+i].M3=M3avg;
#ifndef BAROTROPIC
#ifdef PRIMITIVE_AVG
    W[0][1][0] = Cons_to_Prim(&(U[0][1][0]));
    W[1][0][0] = Cons_to_Prim(&(U[1][0][0]));

    Pavg = fact*(W[0][1][0].P+W[1][0][0].P);

    Eavg = Pavg/Gamma_1 + 0.5*(SQR(M1avg)+SQR(M2avg)+SQR(M3avg))/davg;
#else /* PRIMITIVE_AVG */
    Eavg  = fact*(U[0][1][0].E+U[1][0][0].E);
#endif /* PRIMITIVE_AVG */
    pG->U[kp+k][jp][ip+i].E = Eavg;
#endif /* BAROTROPIC */
  }}

/* kp slice */
  for(j=-1;j<=1;j+=2){
  for(i=-1;i<=1;i+=2){
    U[0][1][0] = pG->U[kp][jp+j  ][ip+i*2];
    U[1][0][0] = pG->U[kp][jp+j*2][ip+i  ];

    davg  = fact*(U[0][1][0].d +U[1][0][0].d );
    M1avg = fact*(U[0][1][0].M1+U[1][0][0].M1);
    M2avg = fact*(U[0][1][0].M2+U[1][0][0].M2);
    M3avg = fact*(U[0][1][0].M3+U[1][0][0].M3);

    pG->U[kp][jp+j][ip+i].d = davg;
    pG->U[kp][jp+j][ip+i].M1=M1avg;
    pG->U[kp][jp+j][ip+i].M2=M2avg;
    pG->U[kp][jp+j][ip+i].M3=M3avg;
#ifndef BAROTROPIC
#ifdef PRIMITIVE_AVG
    W[0][1][0] = Cons_to_Prim(&(U[0][1][0]));
    W[1][0][0] = Cons_to_Prim(&(U[1][0][0]));

    Pavg = fact*(W[0][1][0].P+W[1][0][0].P);

    Eavg = Pavg/Gamma_1 + 0.5*(SQR(M1avg)+SQR(M2avg)+SQR(M3avg))/davg;
#else /* PRIMITIVE_AVG */
    Eavg  = fact*(U[0][1][0].E+U[1][0][0].E);
#endif /* PRIMITIVE_AVG */
    pG->U[kp][jp+j][ip+i].E = Eavg;
#endif /* BAROTROPIC */
  }}
      
  /* Reset density and momentum of the zone where the sink particle resides */
  fact =1/6.0;

  U[0][0][1] = pG->U[kp  ][jp  ][ip+1];
  U[0][1][0] = pG->U[kp  ][jp+1][ip  ];
  U[1][0][0] = pG->U[kp+1][jp  ][ip  ];
  U[0][1][1] = pG->U[kp-1][jp  ][ip  ];
  U[1][0][1] = pG->U[kp  ][jp-1][ip  ];
  U[1][1][0] = pG->U[kp  ][jp  ][ip-1];


  davg  = fact*((U[0][0][1].d+U[0][1][0].d+U[1][0][0].d) +
                (U[0][1][1].d+U[1][1][0].d+U[1][0][1].d));
/*
  if(pStar->isnew){
    M1avg = fact*((U[0][0][1].M1+U[0][1][0].M1+U[1][0][0].M1) +
                  (U[0][1][1].M1+U[1][1][0].M1+U[1][0][1].M1));
    M2avg = fact*((U[0][0][1].M2+U[0][1][0].M2+U[1][0][0].M2) +
                  (U[0][1][1].M2+U[1][1][0].M2+U[1][0][1].M2));
    M3avg = fact*((U[0][0][1].M3+U[0][1][0].M3+U[1][0][0].M3) +
                  (U[0][1][1].M3+U[1][1][0].M3+U[1][0][1].M3));
  } else {
*/
    M1avg = davg*pStar->v1;
    M2avg = davg*pStar->v2;
    M3avg = davg*pStar->v3;
//  }

#ifndef BAROTROPIC
#ifdef PRIMITIVE_AVG
  W[0][0][1] = Cons_to_Prim(&(U[0][0][1]));
  W[0][1][0] = Cons_to_Prim(&(U[0][1][0]));
  W[1][0][0] = Cons_to_Prim(&(U[1][0][0]));
  W[0][1][1] = Cons_to_Prim(&(U[0][1][1]));
  W[1][0][1] = Cons_to_Prim(&(U[1][0][1]));
  W[1][1][0] = Cons_to_Prim(&(U[1][1][0]));

  Pavg = fact*((W[0][0][1].P+W[0][1][0].P+W[1][0][0].P) +
               (W[0][1][1].P+W[1][1][0].P+W[1][0][1].P));

  Eavg = Pavg/Gamma_1 + 0.5*(SQR(M1avg)+SQR(M2avg)+SQR(M3avg))/davg;
#else /* PRIMITIVE_AVG */
  Eavg  = fact*((U[0][0][1].E+U[0][1][0].E+U[1][0][0].E) +
                (U[0][1][1].E+U[1][1][0].E+U[1][0][1].E));
#endif /* PRIMITIVE_AVG */
#endif /* BAROTROPIC */

  pG->U[kp][jp][ip].d = davg;
  pG->U[kp][jp][ip].M1=M1avg;
  pG->U[kp][jp][ip].M2=M2avg;
  pG->U[kp][jp][ip].M3=M3avg;
#ifndef BAROTROPIC
  pG->U[kp][jp][ip].E = Eavg;
#endif
 
}

#else /* FACE_ONLY */

void set_ghost_region(GridS *pG, StarParS *pStar){

  int i, j, k;
  int ip, jp, kp;
  Real fact,w1,w2,w3,davg,M1avg,M2avg,M3avg,Pavg,Eavg;
  ConsS U[2][2][2];
#ifdef PRIMITIVE_AVG
  PrimS W[2][2][2];
#endif

  cc_ijk(pG,pStar->x1,pStar->x2,pStar->x3,&ip,&jp,&kp);

  /* 6 front faces */
  pG->U[kp][jp][ip+1] = pG->U[kp][jp][ip+2];
  pG->U[kp][jp][ip-1] = pG->U[kp][jp][ip-2];
      
  pG->U[kp][jp+1][ip] = pG->U[kp][jp+2][ip];
  pG->U[kp][jp-1][ip] = pG->U[kp][jp-2][ip];
     
  pG->U[kp+1][jp][ip] = pG->U[kp+2][jp][ip];
  pG->U[kp-1][jp][ip] = pG->U[kp-2][jp][ip];
      
  w1=6.0; w2=3.0; w3=2.0; // distance squared weight
  fact=1.0/(w1*3+w2*3+w3);

  for(k=-1;k<=1;k+=2){
  for(j=-1;j<=1;j+=2){
  for(i=-1;i<=1;i+=2){
// weigted average for all nearby 7 zones
    U[0][0][1] = pG->U[kp+k  ][jp+j  ][ip+i*2];
    U[0][1][0] = pG->U[kp+k  ][jp+j*2][ip+i  ];
    U[1][0][0] = pG->U[kp+k*2][jp+j  ][ip+i  ];
    U[0][1][1] = pG->U[kp+k  ][jp+j*2][ip+i*2];
    U[1][0][1] = pG->U[kp+k*2][jp+j  ][ip+i*2];
    U[1][1][0] = pG->U[kp+k*2][jp+j*2][ip+i  ];
    U[1][1][1] = pG->U[kp+k*2][jp+j*2][ip+i*2];

    davg  = fact*(w1*(U[0][0][1].d+U[0][1][0].d+U[1][0][0].d) +
                  w2*(U[0][1][1].d+U[1][1][0].d+U[1][0][1].d) +
                  w3*U[1][1][1].d);
    M1avg = fact*(w1*(U[0][0][1].M1+U[0][1][0].M1+U[1][0][0].M1) +
                  w2*(U[0][1][1].M1+U[1][1][0].M1+U[1][0][1].M1) +
                  w3*U[1][1][1].M1);
    M2avg = fact*(w1*(U[0][0][1].M2+U[0][1][0].M2+U[1][0][0].M2) +
                  w2*(U[0][1][1].M2+U[1][1][0].M2+U[1][0][1].M2) +
                  w3*U[1][1][1].M2);
    M3avg = fact*(w1*(U[0][0][1].M3+U[0][1][0].M3+U[1][0][0].M3) +
                  w2*(U[0][1][1].M3+U[1][1][0].M3+U[1][0][1].M3) +
                  w3*U[1][1][1].M3);

    pG->U[kp+k][jp+j][ip+i].d = davg;
    pG->U[kp+k][jp+j][ip+i].M1=M1avg;
    pG->U[kp+k][jp+j][ip+i].M2=M2avg;
    pG->U[kp+k][jp+j][ip+i].M3=M3avg;
#ifndef BAROTROPIC
#ifdef PRIMITIVE_AVG
    W[0][0][0] = Cons_to_Prim(&(U[0][0][0]));
    W[0][0][1] = Cons_to_Prim(&(U[0][0][1]));
    W[0][1][0] = Cons_to_Prim(&(U[0][1][0]));
    W[1][0][0] = Cons_to_Prim(&(U[1][0][0]));
    W[0][1][1] = Cons_to_Prim(&(U[0][1][1]));
    W[1][0][1] = Cons_to_Prim(&(U[1][0][1]));
    W[1][1][0] = Cons_to_Prim(&(U[1][1][0]));
    W[1][1][1] = Cons_to_Prim(&(U[1][1][1]));

    Pavg = fact*(w1*(W[0][0][1].P+W[0][1][0].P+W[1][0][0].P) +
                 w2*(W[0][1][1].P+W[1][1][0].P+W[1][0][1].P) +
                 w3*W[1][1][1].P);

    Eavg = Pavg/Gamma_1 + 0.5*(SQR(M1avg)+SQR(M2avg)+SQR(M3avg))/davg;
#else /* PRIMITIVE_AVG */
    Eavg  = fact*(w1*(U[0][0][1].E+U[0][1][0].E+U[1][0][0].E) +
                  w2*(U[0][1][1].E+U[1][1][0].E+U[1][0][1].E) +
                  w3*U[1][1][1].E);
#endif /* PRIMITIVE_AVG */
    pG->U[kp+k][jp+j][ip+i].E = Eavg;
#endif /* BAROTROPIC */
  }}}
      
/* 4 sides for 3 middle-slices */
  w1=2.0; w2=1.0; w3=0.0; // distance squared weight
  fact=1.0/(w1*2+w2);

/* ip slice */
  for(k=-1;k<=1;k+=2){
  for(j=-1;j<=1;j+=2){
    U[0][1][0] = pG->U[kp+k  ][jp+j*2][ip];
    U[1][0][0] = pG->U[kp+k*2][jp+j  ][ip];
    U[1][1][0] = pG->U[kp+k*2][jp+j*2][ip];

    davg  = fact*(w1*(U[0][1][0].d+U[1][0][0].d) + w2*U[1][1][0].d);
    M1avg = fact*(w1*(U[0][1][0].M1+U[1][0][0].M1) + w2*U[1][1][0].M1);
    M2avg = fact*(w1*(U[0][1][0].M2+U[1][0][0].M2) + w2*U[1][1][0].M2);
    M3avg = fact*(w1*(U[0][1][0].M3+U[1][0][0].M3) + w2*U[1][1][0].M3);

    pG->U[kp+k][jp+j][ip].d = davg;
    pG->U[kp+k][jp+j][ip].M1=M1avg;
    pG->U[kp+k][jp+j][ip].M2=M2avg;
    pG->U[kp+k][jp+j][ip].M3=M3avg;
#ifndef BAROTROPIC
#ifdef PRIMITIVE_AVG
    W[0][1][0] = Cons_to_Prim(&(U[0][1][0]));
    W[1][0][0] = Cons_to_Prim(&(U[1][0][0]));
    W[1][1][0] = Cons_to_Prim(&(U[1][1][0]));

    Pavg = fact*(w1*(W[0][1][0].P+W[1][0][0].P) +
                 w2*W[1][1][0].P);

    Eavg = Pavg/Gamma_1 + 0.5*(SQR(M1avg)+SQR(M2avg)+SQR(M3avg))/davg;
#else /* PRIMITIVE_AVG */
    Eavg = fact*(w1*(U[0][1][0].E+U[1][0][0].E) +
                 w2*U[1][1][0].E);
#endif /* PRIMITIVE_AVG */
    pG->U[kp+k][jp+j][ip].E = Eavg;
#endif /* BAROTROPIC */
  }}
 
/* jp slice */
  for(k=-1;k<=1;k+=2){
  for(i=-1;i<=1;i+=2){
    U[0][1][0] = pG->U[kp+k  ][jp][ip+i*2];
    U[1][0][0] = pG->U[kp+k*2][jp][ip+i  ];
    U[1][1][0] = pG->U[kp+k*2][jp][ip+i*2];

    davg  = fact*(w1*(U[0][1][0].d+U[1][0][0].d) + w2*U[1][1][0].d);
    M1avg = fact*(w1*(U[0][1][0].M1+U[1][0][0].M1) + w2*U[1][1][0].M1);
    M2avg = fact*(w1*(U[0][1][0].M2+U[1][0][0].M2) + w2*U[1][1][0].M2);
    M3avg = fact*(w1*(U[0][1][0].M3+U[1][0][0].M3) + w2*U[1][1][0].M3);

    pG->U[kp+k][jp][ip+i].d = davg;
    pG->U[kp+k][jp][ip+i].M1=M1avg;
    pG->U[kp+k][jp][ip+i].M2=M2avg;
    pG->U[kp+k][jp][ip+i].M3=M3avg;
#ifndef BAROTROPIC
#ifdef PRIMITIVE_AVG
    W[0][1][0] = Cons_to_Prim(&(U[0][1][0]));
    W[1][0][0] = Cons_to_Prim(&(U[1][0][0]));
    W[1][1][0] = Cons_to_Prim(&(U[1][1][0]));

    Pavg = fact*(w1*(W[0][1][0].P+W[1][0][0].P) +
                 w2*W[1][1][0].P);

    Eavg = Pavg/Gamma_1 + 0.5*(SQR(M1avg)+SQR(M2avg)+SQR(M3avg))/davg;
#else /* PRIMITIVE_AVG */
    Eavg  = fact*(w1*(U[0][1][0].E+U[1][0][0].E) +
                  w2*U[1][1][0].E);
#endif /* PRIMITIVE_AVG */
    pG->U[kp+k][jp][ip+i].E = Eavg;
#endif /* BAROTROPIC */
  }}

/* kp slice */
  for(j=-1;j<=1;j+=2){
  for(i=-1;i<=1;i+=2){
    U[0][1][0] = pG->U[kp][jp+j  ][ip+i*2];
    U[1][0][0] = pG->U[kp][jp+j*2][ip+i  ];
    U[1][1][0] = pG->U[kp][jp+j*2][ip+i*2];

    davg  = fact*(w1*(U[0][1][0].d+U[1][0][0].d) + w2*U[1][1][0].d);
    M1avg = fact*(w1*(U[0][1][0].M1+U[1][0][0].M1) + w2*U[1][1][0].M1);
    M2avg = fact*(w1*(U[0][1][0].M2+U[1][0][0].M2) + w2*U[1][1][0].M2);
    M3avg = fact*(w1*(U[0][1][0].M3+U[1][0][0].M3) + w2*U[1][1][0].M3);

    pG->U[kp][jp+j][ip+i].d = davg;
    pG->U[kp][jp+j][ip+i].M1=M1avg;
    pG->U[kp][jp+j][ip+i].M2=M2avg;
    pG->U[kp][jp+j][ip+i].M3=M3avg;
#ifndef BAROTROPIC
#ifdef PRIMITIVE_AVG
    W[0][1][0] = Cons_to_Prim(&(U[0][1][0]));
    W[1][0][0] = Cons_to_Prim(&(U[1][0][0]));
    W[1][1][0] = Cons_to_Prim(&(U[1][1][0]));

    Pavg = fact*(w1*(W[0][1][0].P+W[1][0][0].P) +
                 w2*W[1][1][0].P);

    Eavg = Pavg/Gamma_1 + 0.5*(SQR(M1avg)+SQR(M2avg)+SQR(M3avg))/davg;
#else /* PRIMITIVE_AVG */
    Eavg  = fact*(w1*(U[0][1][0].E+U[1][0][0].E) +
                  w2*U[1][1][0].E);
#endif /* PRIMITIVE_AVG */
    pG->U[kp][jp+j][ip+i].E = Eavg;
#endif /* BAROTROPIC */
  }}
      
  davg =0.;
  M1avg =0.;
  M2avg =0.;
  M3avg =0.;
  fact = 0.;
#ifndef BAROTROPIC
  Pavg =0.;
  Eavg =0.;
#endif
  for(k=-1; k<=1; k++){
  for(j=-1; j<=1; j++){
  for(i=-1; i<=1; i++){
    w1 =(Real)(SQR(i)+SQR(j)+SQR(k));
    if(w1 != 0) {
      w1 = 6.0/w1;
      fact += w1;
      davg += w1*pG->U[kp+k][jp+j][ip+i].d;
      M1avg += w1*pG->U[kp+k][jp+j][ip+i].M1;
      M2avg += w1*pG->U[kp+k][jp+j][ip+i].M2;
      M3avg += w1*pG->U[kp+k][jp+j][ip+i].M3;
#ifndef BAROTROPIC
#ifdef PRIMITIVE_AVG
      W[0][0][0] = Cons_to_Prim(&(pG->U[kp+k][jp+j][ip+i]));
      Pavg += w1*W[0][0][0].P;
#else
      Eavg += w1*pG->U[kp+k][jp+j][ip+i].E;
#endif
#endif /* BAROTROPIC */
    }
  }}}

  fact = 1./fact;
  davg *= fact;
/*
  M1avg *= fact;
  M2avg *= fact;
  M3avg *= fact;
*/
  M1avg = davg*pStar->v1;
  M2avg = davg*pStar->v2;
  M3avg = davg*pStar->v3;
#ifndef BAROTROPIC
#ifdef PRIMITIVE_AVG
  Pavg *= fact;
  Eavg = Pavg/Gamma_1 + 0.5*(SQR(M1avg)+SQR(M2avg)+SQR(M3avg))/davg;
  if (Pavg < 0.0) ath_error("[modify_ghost_region] Pressure becomes negative at SP ghost region\n");
#else
  Eavg *= fact;
#endif
#endif /* BAROTROPIC */

  pG->U[kp][jp][ip].d = davg;
  pG->U[kp][jp][ip].M1=M1avg;
  pG->U[kp][jp][ip].M2=M2avg;
  pG->U[kp][jp][ip].M3=M3avg;
#ifndef BAROTROPIC
  pG->U[kp][jp][ip].E = Eavg;
#endif
 
}

#endif /* FACE_ONLY */

#if NGHOST_STARP == 2
void cal_avg_nghost2(GridS *pG,StarParS *pStar, const int ip, const int jp, const int kp){

  int i, j, k;
/**************************************************************************/
    /* Outer cells */
    /* ip+-2 face */
    
      for (k=kp-1; k<=kp+1; k++) {
        for (j=jp-1; j<=jp+1; j++) {
          pG->U[k][j][ip-2] = pG->U[k][j][ip-3];
          pG->U[k][j][ip+2] = pG->U[k][j][ip+3];
        }
      }
      
      /* jp+-2 face */
      for (k=kp-1; k<=kp+1; k++) {
        for (i=ip-1; i<=ip+1; i++) {
          pG->U[k][jp-2][i] = pG->U[k][jp-3][i];
          pG->U[k][jp+2][i] = pG->U[k][jp+3][i];
        }
      }
      
      /* kp+-2 face */
      for (j=jp-1; j<=jp+1; j++) {
        for (i=ip-1; i<=ip+1; i++) {
          pG->U[kp-2][j][i] = pG->U[kp-3][j][i];
          pG->U[kp+2][j][i] = pG->U[kp+3][j][i];
        }
      }
      
      /* kp-2, kp+2: jp-2 & jp+2 */
      for (i=ip-1; i<=ip+1; i++) {
        pG->U[kp-2][jp-2][i].d  = ONE_3RD*(pG->U[kp-2][jp-3][i].d  +
                                           pG->U[kp-3][jp-2][i].d  +
                                           pG->U[kp-3][jp-3][i].d );
        pG->U[kp-2][jp-2][i].M1 = ONE_3RD*(pG->U[kp-2][jp-3][i].M1 +
                                           pG->U[kp-3][jp-2][i].M1 +
                                           pG->U[kp-3][jp-3][i].M1);
        pG->U[kp-2][jp-2][i].M2 = ONE_3RD*(pG->U[kp-2][jp-3][i].M2 +
                                           pG->U[kp-3][jp-2][i].M2 +
                                           pG->U[kp-3][jp-3][i].M2);
        pG->U[kp-2][jp-2][i].M3 = ONE_3RD*(pG->U[kp-2][jp-3][i].M3 +
                                           pG->U[kp-3][jp-2][i].M3 +
                                           pG->U[kp-3][jp-3][i].M3);
        
        pG->U[kp-2][jp+2][i].d  = ONE_3RD*(pG->U[kp-2][jp+3][i].d  +
                                           pG->U[kp-3][jp+2][i].d  +
                                           pG->U[kp-3][jp+3][i].d );
        pG->U[kp-2][jp+2][i].M1 = ONE_3RD*(pG->U[kp-2][jp+3][i].M1 +
                                           pG->U[kp-3][jp+2][i].M1 +
                                           pG->U[kp-3][jp+3][i].M1);
        pG->U[kp-2][jp+2][i].M2 = ONE_3RD*(pG->U[kp-2][jp+3][i].M2 +
                                           pG->U[kp-3][jp+2][i].M2 +
                                           pG->U[kp-3][jp+3][i].M2);
        pG->U[kp-2][jp+2][i].M3 = ONE_3RD*(pG->U[kp-2][jp+3][i].M3 +
                                           pG->U[kp-3][jp+2][i].M3 +
                                           pG->U[kp-3][jp+3][i].M3);
        
        pG->U[kp+2][jp-2][i].d  = ONE_3RD*(pG->U[kp+2][jp-3][i].d  +
                                           pG->U[kp+3][jp-2][i].d  +
                                           pG->U[kp+3][jp-3][i].d );
        pG->U[kp+2][jp-2][i].M1 = ONE_3RD*(pG->U[kp+2][jp-3][i].M1 +
                                           pG->U[kp+3][jp-2][i].M1 +
                                           pG->U[kp+3][jp-3][i].M1);
        pG->U[kp+2][jp-2][i].M2 = ONE_3RD*(pG->U[kp+2][jp-3][i].M2 +
                                           pG->U[kp+3][jp-2][i].M2 +
                                           pG->U[kp+3][jp-3][i].M2);
        pG->U[kp+2][jp-2][i].M3 = ONE_3RD*(pG->U[kp+2][jp-3][i].M3 +
                                           pG->U[kp+3][jp-2][i].M3 +
                                           pG->U[kp+3][jp-3][i].M3);
        
        pG->U[kp+2][jp+2][i].d  = ONE_3RD*(pG->U[kp+2][jp+3][i].d  +
                                           pG->U[kp+3][jp+2][i].d  +
                                           pG->U[kp+3][jp+3][i].d );
        pG->U[kp+2][jp+2][i].M1 = ONE_3RD*(pG->U[kp+2][jp+3][i].M1 +
                                           pG->U[kp+3][jp+2][i].M1 +
                                           pG->U[kp+3][jp+3][i].M1);
        pG->U[kp+2][jp+2][i].M2 = ONE_3RD*(pG->U[kp+2][jp+3][i].M2 +
                                           pG->U[kp+3][jp+2][i].M2 +
                                           pG->U[kp+3][jp+3][i].M2);
        pG->U[kp+2][jp+2][i].M3 = ONE_3RD*(pG->U[kp+2][jp+3][i].M3 +
                                           pG->U[kp+3][jp+2][i].M3 +
                                           pG->U[kp+3][jp+3][i].M3);
      }
      
      
      /* kp-2, kp+2: ip-2 & ip+2 */
      for (j=jp-1; j<=jp+1; j++) {
        pG->U[kp-2][j][ip-2].d  = ONE_3RD*(pG->U[kp-2][j][ip-3].d  +
                                           pG->U[kp-3][j][ip-2].d  +
                                           pG->U[kp-3][j][ip-3].d );
        pG->U[kp-2][j][ip-2].M1 = ONE_3RD*(pG->U[kp-2][j][ip-3].M1 +
                                           pG->U[kp-3][j][ip-2].M1 +
                                           pG->U[kp-3][j][ip-3].M1);
        pG->U[kp-2][j][ip-2].M2 = ONE_3RD*(pG->U[kp-2][j][ip-3].M2 +
                                           pG->U[kp-3][j][ip-2].M2 +
                                           pG->U[kp-3][j][ip-3].M2);
        pG->U[kp-2][j][ip-2].M3 = ONE_3RD*(pG->U[kp-2][j][ip-3].M3 +
                                           pG->U[kp-3][j][ip-2].M3 +
                                           pG->U[kp-3][j][ip-3].M3);
        
        pG->U[kp-2][j][ip+2].d  = ONE_3RD*(pG->U[kp-2][j][ip+3].d  +
                                           pG->U[kp-3][j][ip-2].d  +
                                           pG->U[kp-3][j][ip-3].d );
        pG->U[kp-2][j][ip+2].M1 = ONE_3RD*(pG->U[kp-2][j][ip+3].M1 +
                                           pG->U[kp-3][j][ip-2].M1 +
                                           pG->U[kp-3][j][ip-3].M1);
        pG->U[kp-2][j][ip+2].M2 = ONE_3RD*(pG->U[kp-2][j][ip+3].M2 +
                                           pG->U[kp-3][j][ip-2].M2 +
                                           pG->U[kp-3][j][ip-3].M2);
        pG->U[kp-2][j][ip+2].M3 = ONE_3RD*(pG->U[kp-2][j][ip+3].M3 +
                                           pG->U[kp-3][j][ip-2].M3 +
                                           pG->U[kp-3][j][ip-3].M3);
        
        pG->U[kp+2][j][ip-2].d  = ONE_3RD*(pG->U[kp+2][j][ip-3].d  +
                                           pG->U[kp+3][j][ip-2].d  +
                                           pG->U[kp+3][j][ip-3].d );
        pG->U[kp+2][j][ip-2].M1 = ONE_3RD*(pG->U[kp+2][j][ip-3].M1 +
                                           pG->U[kp+3][j][ip-2].M1 +
                                           pG->U[kp+3][j][ip-3].M1);
        pG->U[kp+2][j][ip-2].M2 = ONE_3RD*(pG->U[kp+2][j][ip-3].M2 +
                                           pG->U[kp+3][j][ip-2].M2 +
                                           pG->U[kp+3][j][ip-3].M2);
        pG->U[kp+2][j][ip-2].M3 = ONE_3RD*(pG->U[kp+2][j][ip-3].M3 +
                                           pG->U[kp+3][j][ip-2].M3 +
                                           pG->U[kp+3][j][ip-3].M3);
        
        pG->U[kp+2][j][ip+2].d  = ONE_3RD*(pG->U[kp+2][j][ip+3].d  +
                                           pG->U[kp+3][j][ip-2].d  +
                                           pG->U[kp+3][j][ip-3].d );
        pG->U[kp+2][j][ip+2].M1 = ONE_3RD*(pG->U[kp+2][j][ip+3].M1 +
                                           pG->U[kp+3][j][ip-2].M1 +
                                           pG->U[kp+3][j][ip-3].M1);
        pG->U[kp+2][j][ip+2].M2 = ONE_3RD*(pG->U[kp+2][j][ip+3].M2 +
                                           pG->U[kp+3][j][ip-2].M2 +
                                           pG->U[kp+3][j][ip-3].M2);
        pG->U[kp+2][j][ip+2].M3 = ONE_3RD*(pG->U[kp+2][j][ip+3].M3 +
                                           pG->U[kp+3][j][ip-2].M3 +
                                           pG->U[kp+3][j][ip-3].M3);
      }
      
      /* jp-2,jp+2: ip-2 & ip+2 */
      for (k=kp-1; k<=kp+1; k++) {
        pG->U[k][jp-2][ip-2].d  = ONE_3RD*(pG->U[k][jp-3][ip-2].d  +
                                           pG->U[k][jp-2][ip-3].d  +
                                           pG->U[k][jp-3][ip-3].d );
        pG->U[k][jp-2][ip-2].M1 = ONE_3RD*(pG->U[k][jp-3][ip-2].M1 +
                                           pG->U[k][jp-2][ip-3].M1 +
                                           pG->U[k][jp-3][ip-3].M1);
        pG->U[k][jp-2][ip-2].M2 = ONE_3RD*(pG->U[k][jp-3][ip-2].M2 +
                                           pG->U[k][jp-2][ip-3].M2 +
                                           pG->U[k][jp-3][ip-3].M2);
        pG->U[k][jp-2][ip-2].M3 = ONE_3RD*(pG->U[k][jp-3][ip-2].M3 +
                                           pG->U[k][jp-2][ip-3].M3 +
                                           pG->U[k][jp-3][ip-3].M3);
        
        pG->U[k][jp-2][ip+2].d  = ONE_3RD*(pG->U[k][jp-3][ip+2].d  +
                                           pG->U[k][jp-2][ip+3].d  +
                                           pG->U[k][jp-3][ip+3].d );
        pG->U[k][jp-2][ip+2].M1 = ONE_3RD*(pG->U[k][jp-3][ip+2].M1 +
                                           pG->U[k][jp-2][ip+3].M1 +
                                           pG->U[k][jp-3][ip+3].M1);
        pG->U[k][jp-2][ip+2].M2 = ONE_3RD*(pG->U[k][jp-3][ip+2].M2 +
                                           pG->U[k][jp-2][ip+3].M2 +
                                           pG->U[k][jp-3][ip+3].M2);
        pG->U[k][jp-2][ip+2].M3 = ONE_3RD*(pG->U[k][jp-3][ip+2].M3 +
                                           pG->U[k][jp-2][ip+3].M3 +
                                           pG->U[k][jp-3][ip+3].M3);
        
        pG->U[k][jp+2][ip-2].d  = ONE_3RD*(pG->U[k][jp+3][ip-2].d  +
                                           pG->U[k][jp+2][ip-3].d  +
                                           pG->U[k][jp+3][ip-3].d );
        pG->U[k][jp+2][ip-2].M1 = ONE_3RD*(pG->U[k][jp+3][ip-2].M1 +
                                           pG->U[k][jp+2][ip-3].M1 +
                                           pG->U[k][jp+3][ip-3].M1);
        pG->U[k][jp+2][ip-2].M2 = ONE_3RD*(pG->U[k][jp+3][ip-2].M2 +
                                           pG->U[k][jp+2][ip-3].M2 +
                                           pG->U[k][jp+3][ip-3].M2);
        pG->U[k][jp+2][ip-2].M3 = ONE_3RD*(pG->U[k][jp+3][ip-2].M3 +
                                           pG->U[k][jp+2][ip-3].M3 +
                                           pG->U[k][jp+3][ip-3].M3);
        
        pG->U[k][jp+2][ip+2].d  = ONE_3RD*(pG->U[k][jp+3][ip+2].d  +
                                           pG->U[k][jp+2][ip+3].d  +
                                           pG->U[k][jp+3][ip+3].d );
        pG->U[k][jp+2][ip+2].M1 = ONE_3RD*(pG->U[k][jp+3][ip+2].M1 +
                                           pG->U[k][jp+2][ip+3].M1 +
                                           pG->U[k][jp+3][ip+3].M1);
        pG->U[k][jp+2][ip+2].M2 = ONE_3RD*(pG->U[k][jp+3][ip+2].M2 +
                                           pG->U[k][jp+2][ip+3].M2 +
                                           pG->U[k][jp+3][ip+3].M2);
        pG->U[k][jp+2][ip+2].M3 = ONE_3RD*(pG->U[k][jp+3][ip+2].M3 +
                                           pG->U[k][jp+2][ip+3].M3 +
                                           pG->U[k][jp+3][ip+3].M3);
      }
      
      /* eight corners outflow */
      pG->U[kp-2][jp-2][ip-2].d  = 0.25*(pG->U[kp-3][jp-2][ip-2].d  +
                                         pG->U[kp-2][jp-3][ip-2].d  +
                                         pG->U[kp-2][jp-2][ip-3].d  +
                                         pG->U[kp-3][jp-3][ip-3].d );
      pG->U[kp-2][jp-2][ip-2].M1 = 0.25*(pG->U[kp-3][jp-2][ip-2].M1 +
                                         pG->U[kp-2][jp-3][ip-2].M1 +
                                         pG->U[kp-2][jp-2][ip-3].M1 +
                                         pG->U[kp-3][jp-3][ip-3].M1);
      pG->U[kp-2][jp-2][ip-2].M2 = 0.25*(pG->U[kp-3][jp-2][ip-2].M2 +
                                         pG->U[kp-2][jp-3][ip-2].M2 +
                                         pG->U[kp-2][jp-2][ip-3].M2 +
                                         pG->U[kp-3][jp-3][ip-3].M2);
      pG->U[kp-2][jp-2][ip-2].M3 = 0.25*(pG->U[kp-3][jp-2][ip-2].M3 +
                                         pG->U[kp-2][jp-3][ip-2].M3 +
                                         pG->U[kp-2][jp-2][ip-3].M3 +
                                         pG->U[kp-3][jp-3][ip-3].M3);
      
      pG->U[kp-2][jp-2][ip+2].d  = 0.25*(pG->U[kp-3][jp-2][ip+2].d  +
                                         pG->U[kp-2][jp-3][ip+2].d  +
                                         pG->U[kp-2][jp-2][ip+3].d  +
                                         pG->U[kp-3][jp-3][ip+3].d );
      pG->U[kp-2][jp-2][ip+2].M1 = 0.25*(pG->U[kp-3][jp-2][ip+2].M1 +
                                         pG->U[kp-2][jp-3][ip+2].M1 +
                                         pG->U[kp-2][jp-2][ip+3].M1 +
                                         pG->U[kp-3][jp-3][ip+3].M1);
      pG->U[kp-2][jp-2][ip+2].M2 = 0.25*(pG->U[kp-3][jp-2][ip+2].M2 +
                                         pG->U[kp-2][jp-3][ip+2].M2 +
                                         pG->U[kp-2][jp-2][ip+3].M2 +
                                         pG->U[kp-3][jp-3][ip+3].M2);
      pG->U[kp-2][jp-2][ip+2].M3 = 0.25*(pG->U[kp-3][jp-2][ip+2].M3 +
                                         pG->U[kp-2][jp-3][ip+2].M3 +
                                         pG->U[kp-2][jp-2][ip+3].M3 +
                                         pG->U[kp-3][jp-3][ip+3].M3);
      
      pG->U[kp-2][jp+2][ip+2].d  = 0.25*(pG->U[kp-3][jp+2][ip+2].d  +
                                         pG->U[kp-2][jp+3][ip+2].d  +
                                         pG->U[kp-2][jp+2][ip+3].d  +
                                         pG->U[kp-3][jp+3][ip+3].d );
      pG->U[kp-2][jp+2][ip+2].M1 = 0.25*(pG->U[kp-3][jp+2][ip+2].M1 +
                                         pG->U[kp-2][jp+3][ip+2].M1 +
                                         pG->U[kp-2][jp+2][ip+3].M1 +
                                         pG->U[kp-3][jp+3][ip+3].M1);
      pG->U[kp-2][jp+2][ip+2].M2 = 0.25*(pG->U[kp-3][jp+2][ip+2].M2 +
                                         pG->U[kp-2][jp+3][ip+2].M2 +
                                         pG->U[kp-2][jp+2][ip+3].M2 +
                                         pG->U[kp-3][jp+3][ip+3].M2);
      pG->U[kp-2][jp+2][ip+2].M3 = 0.25*(pG->U[kp-3][jp+2][ip+2].M3 +
                                         pG->U[kp-2][jp+3][ip+2].M3 +
                                         pG->U[kp-2][jp+2][ip+3].M3 +
                                         pG->U[kp-3][jp+3][ip+3].M3);
      
      pG->U[kp-2][jp+2][ip-2].d  = 0.25*(pG->U[kp-3][jp+2][ip-2].d  +
                                         pG->U[kp-2][jp+3][ip-2].d  +
                                         pG->U[kp-2][jp+2][ip-3].d  +
                                         pG->U[kp-3][jp+3][ip-3].d );
      pG->U[kp-2][jp+2][ip-2].M1 = 0.25*(pG->U[kp-3][jp+2][ip-2].M1 +
                                         pG->U[kp-2][jp+3][ip-2].M1 +
                                         pG->U[kp-2][jp+2][ip-3].M1 +
                                         pG->U[kp-3][jp+3][ip-3].M1);
      pG->U[kp-2][jp+2][ip-2].M2 = 0.25*(pG->U[kp-3][jp+2][ip-2].M2 +
                                         pG->U[kp-2][jp+3][ip-2].M2 +
                                         pG->U[kp-2][jp+2][ip-3].M2 +
                                         pG->U[kp-3][jp+3][ip-3].M2);
      pG->U[kp-2][jp+2][ip-2].M3 = 0.25*(pG->U[kp-3][jp+2][ip-2].M3 +
                                         pG->U[kp-2][jp+3][ip-2].M3 +
                                         pG->U[kp-2][jp+2][ip-3].M3 +
                                         pG->U[kp-3][jp+3][ip-3].M3);
      
      pG->U[kp+2][jp-2][ip-2].d  = 0.25*(pG->U[kp+3][jp-2][ip-2].d  +
                                         pG->U[kp+2][jp-3][ip-2].d  +
                                         pG->U[kp+2][jp-2][ip-3].d  +
                                         pG->U[kp+3][jp-3][ip-3].d );
      pG->U[kp+2][jp-2][ip-2].M1 = 0.25*(pG->U[kp+3][jp-2][ip-2].M1 +
                                         pG->U[kp+2][jp-3][ip-2].M1 +
                                         pG->U[kp+2][jp-2][ip-3].M1 +
                                         pG->U[kp+3][jp-3][ip-3].M1);
      pG->U[kp+2][jp-2][ip-2].M2 = 0.25*(pG->U[kp+3][jp-2][ip-2].M2 +
                                         pG->U[kp+2][jp-3][ip-2].M2 +
                                         pG->U[kp+2][jp-2][ip-3].M2 +
                                         pG->U[kp+3][jp-3][ip-3].M2);
      pG->U[kp+2][jp-2][ip-2].M3 = 0.25*(pG->U[kp+3][jp-2][ip-2].M3 +
                                         pG->U[kp+2][jp-3][ip-2].M3 +
                                         pG->U[kp+2][jp-2][ip-3].M3 +
                                         pG->U[kp+3][jp-3][ip-3].M3);
      
      pG->U[kp+2][jp-2][ip+2].d  = 0.25*(pG->U[kp+3][jp-2][ip+2].d  +
                                         pG->U[kp+2][jp-3][ip+2].d  +
                                         pG->U[kp+2][jp-2][ip+3].d  +
                                         pG->U[kp+3][jp-3][ip+3].d );
      pG->U[kp+2][jp-2][ip+2].M1 = 0.25*(pG->U[kp+3][jp-2][ip+2].M1 +
                                         pG->U[kp+2][jp-3][ip+2].M1 +
                                         pG->U[kp+2][jp-2][ip+3].M1 +
                                         pG->U[kp+3][jp-3][ip+3].M1);
      pG->U[kp+2][jp-2][ip+2].M2 = 0.25*(pG->U[kp+3][jp-2][ip+2].M2 +
                                         pG->U[kp+2][jp-3][ip+2].M2 +
                                         pG->U[kp+2][jp-2][ip+3].M2 +
                                         pG->U[kp+3][jp-3][ip+3].M2);
      pG->U[kp+2][jp-2][ip+2].M3 = 0.25*(pG->U[kp+3][jp-2][ip+2].M3 +
                                         pG->U[kp+2][jp-3][ip+2].M3 +
                                         pG->U[kp+2][jp-2][ip+3].M3 +
                                         pG->U[kp+3][jp-3][ip+3].M3);
      
      pG->U[kp+2][jp+2][ip+2].d  = 0.25*(pG->U[kp+3][jp+2][ip+2].d  +
                                         pG->U[kp+2][jp+3][ip+2].d  +
                                         pG->U[kp+2][jp+2][ip+3].d  +
                                         pG->U[kp+3][jp+3][ip+3].d );
      pG->U[kp+2][jp+2][ip+2].M1 = 0.25*(pG->U[kp+3][jp+2][ip+2].M1 +
                                         pG->U[kp+2][jp+3][ip+2].M1 +
                                         pG->U[kp+2][jp+2][ip+3].M1 +
                                         pG->U[kp+3][jp+3][ip+3].M1);
      pG->U[kp+2][jp+2][ip+2].M2 = 0.25*(pG->U[kp+3][jp+2][ip+2].M2 +
                                         pG->U[kp+2][jp+3][ip+2].M2 +
                                         pG->U[kp+2][jp+2][ip+3].M2 +
                                         pG->U[kp+3][jp+3][ip+3].M2);
      pG->U[kp+2][jp+2][ip+2].M3 = 0.25*(pG->U[kp+3][jp+2][ip+2].M3 +
                                         pG->U[kp+2][jp+3][ip+2].M3 +
                                         pG->U[kp+2][jp+2][ip+3].M3 +
                                         pG->U[kp+3][jp+3][ip+3].M3);
      
      pG->U[kp+2][jp+2][ip-2].d  = 0.25*(pG->U[kp+3][jp+2][ip-2].d  +
                                         pG->U[kp+2][jp+3][ip-2].d  +
                                         pG->U[kp+2][jp+2][ip-3].d  +
                                         pG->U[kp+3][jp+3][ip-3].d );
      pG->U[kp+2][jp+2][ip-2].M1 = 0.25*(pG->U[kp+3][jp+2][ip-2].M1 +
                                         pG->U[kp+2][jp+3][ip-2].M1 +
                                         pG->U[kp+2][jp+2][ip-3].M1 +
                                         pG->U[kp+3][jp+3][ip-3].M1);
      pG->U[kp+2][jp+2][ip-2].M2 = 0.25*(pG->U[kp+3][jp+2][ip-2].M2 +
                                         pG->U[kp+2][jp+3][ip-2].M2 +
                                         pG->U[kp+2][jp+2][ip-3].M2 +
                                         pG->U[kp+3][jp+3][ip-3].M2);
      pG->U[kp+2][jp+2][ip-2].M3 = 0.25*(pG->U[kp+3][jp+2][ip-2].M3 +
                                         pG->U[kp+2][jp+3][ip-2].M3 +
                                         pG->U[kp+2][jp+2][ip-3].M3 +
                                         pG->U[kp+3][jp+3][ip-3].M3);
      
      /***********************************************************/
      /*inner cells */
      
      /* x direction outflow */
      pG->U[kp][jp][ip-1] = pG->U[kp][jp][ip-2];
      pG->U[kp][jp][ip+1] = pG->U[kp][jp][ip+2];
      
      /* y direction outflow */
      pG->U[kp][jp-1][ip] = pG->U[kp][jp-2][ip];
      pG->U[kp][jp+1][ip] = pG->U[kp][jp+2][ip];
      
      /* z direction outflow */
      pG->U[kp-1][jp][ip] = pG->U[kp-2][jp][ip];
      pG->U[kp+1][jp][ip] = pG->U[kp+2][jp][ip];
      
      /* eight corners */
      pG->U[kp-1][jp-1][ip-1].d  = 0.25*(pG->U[kp-2][jp-1][ip-1].d  +
                                         pG->U[kp-1][jp-2][ip-1].d  +
                                         pG->U[kp-1][jp-1][ip-2].d  +
                                         pG->U[kp-2][jp-2][ip-2].d );
      pG->U[kp-1][jp-1][ip-1].M1 = 0.25*(pG->U[kp-2][jp-1][ip-1].M1 +
                                         pG->U[kp-1][jp-2][ip-1].M1 +
                                         pG->U[kp-1][jp-1][ip-2].M1 +
                                         pG->U[kp-2][jp-2][ip-2].M1);
      pG->U[kp-1][jp-1][ip-1].M2 = 0.25*(pG->U[kp-2][jp-1][ip-1].M2 +
                                         pG->U[kp-1][jp-2][ip-1].M2 +
                                         pG->U[kp-1][jp-1][ip-2].M2 +
                                         pG->U[kp-2][jp-2][ip-2].M2);
      pG->U[kp-1][jp-1][ip-1].M3 = 0.25*(pG->U[kp-2][jp-1][ip-1].M3 +
                                         pG->U[kp-1][jp-2][ip-1].M3 +
                                         pG->U[kp-1][jp-1][ip-2].M3 +
                                         pG->U[kp-2][jp-2][ip-2].M3);
      
      pG->U[kp-1][jp-1][ip+1].d  = 0.25*(pG->U[kp-2][jp-1][ip+1].d  +
                                         pG->U[kp-1][jp-2][ip+1].d  +
                                         pG->U[kp-1][jp-1][ip+2].d  +
                                         pG->U[kp-2][jp-2][ip+2].d );
      pG->U[kp-1][jp-1][ip+1].M1 = 0.25*(pG->U[kp-2][jp-1][ip+1].M1 +
                                         pG->U[kp-1][jp-2][ip+1].M1 +
                                         pG->U[kp-1][jp-1][ip+2].M1 +
                                         pG->U[kp-2][jp-2][ip+2].M1);
      pG->U[kp-1][jp-1][ip+1].M2 = 0.25*(pG->U[kp-2][jp-1][ip+1].M2 +
                                         pG->U[kp-1][jp-2][ip+1].M2 +
                                         pG->U[kp-1][jp-1][ip+2].M2 +
                                         pG->U[kp-2][jp-2][ip+2].M2);
      pG->U[kp-1][jp-1][ip+1].M3 = 0.25*(pG->U[kp-2][jp-1][ip+1].M3 +
                                         pG->U[kp-1][jp-2][ip+1].M3 +
                                         pG->U[kp-1][jp-1][ip+2].M3 +
                                         pG->U[kp-2][jp-2][ip+2].M3);
      
      pG->U[kp-1][jp+1][ip+1].d  = 0.25*(pG->U[kp-2][jp+1][ip+1].d  +
                                         pG->U[kp-1][jp+2][ip+1].d  +
                                         pG->U[kp-1][jp+1][ip+2].d  +
                                         pG->U[kp-2][jp+2][ip+2].d );
      pG->U[kp-1][jp+1][ip+1].M1 = 0.25*(pG->U[kp-2][jp+1][ip+1].M1 +
                                         pG->U[kp-1][jp+2][ip+1].M1 +
                                         pG->U[kp-1][jp+1][ip+2].M1 +
                                         pG->U[kp-2][jp+2][ip+2].M1);
      pG->U[kp-1][jp+1][ip+1].M2 = 0.25*(pG->U[kp-2][jp+1][ip+1].M2 +
                                         pG->U[kp-1][jp+2][ip+1].M2 +
                                         pG->U[kp-1][jp+1][ip+2].M2 +
                                         pG->U[kp-2][jp+2][ip+2].M2);
      pG->U[kp-1][jp+1][ip+1].M3 = 0.25*(pG->U[kp-2][jp+1][ip+1].M3 +
                                         pG->U[kp-1][jp+2][ip+1].M3 +
                                         pG->U[kp-1][jp+1][ip+2].M3 +
                                         pG->U[kp-2][jp+2][ip+2].M3);
      
      pG->U[kp-1][jp+1][ip-1].d  = 0.25*(pG->U[kp-2][jp+1][ip-1].d  +
                                         pG->U[kp-1][jp+2][ip-1].d  +
                                         pG->U[kp-1][jp+1][ip-2].d  +
                                         pG->U[kp-2][jp+2][ip-2].d );
      pG->U[kp-1][jp+1][ip-1].M1 = 0.25*(pG->U[kp-2][jp+1][ip-1].M1 +
                                         pG->U[kp-1][jp+2][ip-1].M1 +
                                         pG->U[kp-1][jp+1][ip-2].M1 +
                                         pG->U[kp-2][jp+2][ip-2].M1);
      pG->U[kp-1][jp+1][ip-1].M2 = 0.25*(pG->U[kp-2][jp+1][ip-1].M2 +
                                         pG->U[kp-1][jp+2][ip-1].M2 +
                                         pG->U[kp-1][jp+1][ip-2].M2 +
                                         pG->U[kp-2][jp+2][ip-2].M2);
      pG->U[kp-1][jp+1][ip-1].M3 = 0.25*(pG->U[kp-2][jp+1][ip-1].M3 +
                                         pG->U[kp-1][jp+2][ip-1].M3 +
                                         pG->U[kp-1][jp+1][ip-2].M3 +
                                         pG->U[kp-2][jp+2][ip-2].M3);
      
      pG->U[kp+1][jp-1][ip-1].d  = 0.25*(pG->U[kp+2][jp-1][ip-1].d  +
                                         pG->U[kp+1][jp-2][ip-1].d  +
                                         pG->U[kp+1][jp-1][ip-2].d  +
                                         pG->U[kp+2][jp-2][ip-2].d );
      pG->U[kp+1][jp-1][ip-1].M1 = 0.25*(pG->U[kp+2][jp-1][ip-1].M1 +
                                         pG->U[kp+1][jp-2][ip-1].M1 +
                                         pG->U[kp+1][jp-1][ip-2].M1 +
                                         pG->U[kp+2][jp-2][ip-2].M1);
      pG->U[kp+1][jp-1][ip-1].M2 = 0.25*(pG->U[kp+2][jp-1][ip-1].M2 +
                                         pG->U[kp+1][jp-2][ip-1].M2 +
                                         pG->U[kp+1][jp-1][ip-2].M2 +
                                         pG->U[kp+2][jp-2][ip-2].M2);
      pG->U[kp+1][jp-1][ip-1].M3 = 0.25*(pG->U[kp+2][jp-1][ip-1].M3 +
                                         pG->U[kp+1][jp-2][ip-1].M3 +
                                         pG->U[kp+1][jp-1][ip-2].M3 +
                                         pG->U[kp+2][jp-2][ip-2].M3);
      
      pG->U[kp+1][jp-1][ip+1].d  = 0.25*(pG->U[kp+2][jp-1][ip+1].d  +
                                         pG->U[kp+1][jp-2][ip+1].d  +
                                         pG->U[kp+1][jp-1][ip+2].d  +
                                         pG->U[kp+2][jp-2][ip+2].d );
      pG->U[kp+1][jp-1][ip+1].M1 = 0.25*(pG->U[kp+2][jp-1][ip+1].M1 +
                                         pG->U[kp+1][jp-2][ip+1].M1 +
                                         pG->U[kp+1][jp-1][ip+2].M1 +
                                         pG->U[kp+2][jp-2][ip+2].M1);
      pG->U[kp+1][jp-1][ip+1].M2 = 0.25*(pG->U[kp+2][jp-1][ip+1].M2 +
                                         pG->U[kp+1][jp-2][ip+1].M2 +
                                         pG->U[kp+1][jp-1][ip+2].M2 +
                                         pG->U[kp+2][jp-2][ip+2].M2);
      pG->U[kp+1][jp-1][ip+1].M3 = 0.25*(pG->U[kp+2][jp-1][ip+1].M3 +
                                         pG->U[kp+1][jp-2][ip+1].M3 +
                                         pG->U[kp+1][jp-1][ip+2].M3 +
                                         pG->U[kp+2][jp-2][ip+2].M3);
      
      pG->U[kp+1][jp+1][ip+1].d  = 0.25*(pG->U[kp+2][jp+1][ip+1].d  +
                                         pG->U[kp+1][jp+2][ip+1].d  +
                                         pG->U[kp+1][jp+1][ip+2].d  +
                                         pG->U[kp+2][jp+2][ip+2].d );
      pG->U[kp+1][jp+1][ip+1].M1 = 0.25*(pG->U[kp+2][jp+1][ip+1].M1 +
                                         pG->U[kp+1][jp+2][ip+1].M1 +
                                         pG->U[kp+1][jp+1][ip+2].M1 +
                                         pG->U[kp+2][jp+2][ip+2].M1);
      pG->U[kp+1][jp+1][ip+1].M2 = 0.25*(pG->U[kp+2][jp+1][ip+1].M2 +
                                         pG->U[kp+1][jp+2][ip+1].M2 +
                                         pG->U[kp+1][jp+1][ip+2].M2 +
                                         pG->U[kp+2][jp+2][ip+2].M2);
      pG->U[kp+1][jp+1][ip+1].M3 = 0.25*(pG->U[kp+2][jp+1][ip+1].M3 +
                                         pG->U[kp+1][jp+2][ip+1].M3 +
                                         pG->U[kp+1][jp+1][ip+2].M3 +
                                         pG->U[kp+2][jp+2][ip+2].M3);
      
      pG->U[kp+1][jp+1][ip-1].d  = 0.25*(pG->U[kp+2][jp+1][ip-1].d  +
                                         pG->U[kp+1][jp+2][ip-1].d  +
                                         pG->U[kp+1][jp+1][ip-2].d  +
                                         pG->U[kp+2][jp+2][ip-2].d );
      pG->U[kp+1][jp+1][ip-1].M1 = 0.25*(pG->U[kp+2][jp+1][ip-1].M1 +
                                         pG->U[kp+1][jp+2][ip-1].M1 +
                                         pG->U[kp+1][jp+1][ip-2].M1 +
                                         pG->U[kp+2][jp+2][ip-2].M1);
      pG->U[kp+1][jp+1][ip-1].M2 = 0.25*(pG->U[kp+2][jp+1][ip-1].M2 +
                                         pG->U[kp+1][jp+2][ip-1].M2 +
                                         pG->U[kp+1][jp+1][ip-2].M2 +
                                         pG->U[kp+2][jp+2][ip-2].M2);
      pG->U[kp+1][jp+1][ip-1].M3 = 0.25*(pG->U[kp+2][jp+1][ip-1].M3 +
                                         pG->U[kp+1][jp+2][ip-1].M3 +
                                         pG->U[kp+1][jp+1][ip-2].M3 +
                                         pG->U[kp+2][jp+2][ip-2].M3);
      
      /* bottom plane outflow */
      pG->U[kp-1][jp-1][ip].d  = ONE_3RD*(pG->U[kp-2][jp-1][ip].d  +
                                          pG->U[kp-1][jp-2][ip].d  +
                                          pG->U[kp-2][jp-2][ip].d );
      pG->U[kp-1][jp-1][ip].M1 = ONE_3RD*(pG->U[kp-2][jp-1][ip].M1 +
                                          pG->U[kp-1][jp-2][ip].M1 +
                                          pG->U[kp-2][jp-2][ip].M1);
      pG->U[kp-1][jp-1][ip].M2 = ONE_3RD*(pG->U[kp-2][jp-1][ip].M2 +
                                          pG->U[kp-1][jp-2][ip].M2 +
                                          pG->U[kp-2][jp-2][ip].M2);
      pG->U[kp-1][jp-1][ip].M3 = ONE_3RD*(pG->U[kp-2][jp-1][ip].M3 +
                                          pG->U[kp-1][jp-2][ip].M3 +
                                          pG->U[kp-2][jp-2][ip].M3);
      
      pG->U[kp-1][jp][ip+1].d  = ONE_3RD*(pG->U[kp-2][jp][ip+1].d  +
                                          pG->U[kp-1][jp][ip+2].d  +
                                          pG->U[kp-2][jp][ip+2].d );
      pG->U[kp-1][jp][ip+1].M1 = ONE_3RD*(pG->U[kp-2][jp][ip+1].M1 +
                                          pG->U[kp-1][jp][ip+2].M1 +
                                          pG->U[kp-2][jp][ip+2].M1);
      pG->U[kp-1][jp][ip+1].M2 = ONE_3RD*(pG->U[kp-2][jp][ip+1].M2 +
                                          pG->U[kp-1][jp][ip+2].M2 +
                                          pG->U[kp-2][jp][ip+2].M2);
      pG->U[kp-1][jp][ip+1].M3 = ONE_3RD*(pG->U[kp-2][jp][ip+1].M3 +
                                          pG->U[kp-1][jp][ip+2].M3 +
                                          pG->U[kp-2][jp][ip+2].M3);
      
      pG->U[kp-1][jp+1][ip].d  = ONE_3RD*(pG->U[kp-2][jp+1][ip].d  +
                                          pG->U[kp-1][jp+2][ip].d  +
                                          pG->U[kp-2][jp+2][ip].d );
      pG->U[kp-1][jp+1][ip].M1 = ONE_3RD*(pG->U[kp-2][jp+1][ip].M1 +
                                          pG->U[kp-1][jp+2][ip].M1 +
                                          pG->U[kp-2][jp+2][ip].M1);
      pG->U[kp-1][jp+1][ip].M2 = ONE_3RD*(pG->U[kp-2][jp+1][ip].M2 +
                                          pG->U[kp-1][jp+2][ip].M2 +
                                          pG->U[kp-2][jp+2][ip].M2);
      pG->U[kp-1][jp+1][ip].M3 = ONE_3RD*(pG->U[kp-2][jp+1][ip].M3 +
                                          pG->U[kp-1][jp+2][ip].M3 +
                                          pG->U[kp-2][jp+2][ip].M3);
      
      pG->U[kp-1][jp][ip-1].d  = ONE_3RD*(pG->U[kp-2][jp][ip-1].d  +
                                          pG->U[kp-1][jp][ip-2].d  +
                                          pG->U[kp-2][jp][ip-2].d );
      pG->U[kp-1][jp][ip-1].M1 = ONE_3RD*(pG->U[kp-2][jp][ip-1].M1 +
                                          pG->U[kp-1][jp][ip-2].M1 +
                                          pG->U[kp-2][jp][ip-2].M1);
      pG->U[kp-1][jp][ip-1].M2 = ONE_3RD*(pG->U[kp-2][jp][ip-1].M2 +
                                          pG->U[kp-1][jp][ip-2].M2 +
                                          pG->U[kp-2][jp][ip-2].M2);
      pG->U[kp-1][jp][ip-1].M3 = ONE_3RD*(pG->U[kp-2][jp][ip-1].M3 +
                                          pG->U[kp-1][jp][ip-2].M3 +
                                          pG->U[kp-2][jp][ip-2].M3);
      
      /* middle plane outflow */
      pG->U[kp][jp-1][ip-1].d  = ONE_3RD*(pG->U[kp][jp-2][ip-1].d  +
                                          pG->U[kp][jp-1][ip-2].d  +
                                          pG->U[kp][jp-2][ip-2].d );
      pG->U[kp][jp-1][ip-1].M1 = ONE_3RD*(pG->U[kp][jp-2][ip-1].M1 +
                                          pG->U[kp][jp-1][ip-2].M1 +
                                          pG->U[kp][jp-2][ip-2].M1);
      pG->U[kp][jp-1][ip-1].M2 = ONE_3RD*(pG->U[kp][jp-2][ip-1].M2 +
                                          pG->U[kp][jp-1][ip-2].M2 +
                                          pG->U[kp][jp-2][ip-2].M2);
      pG->U[kp][jp-1][ip-1].M3 = ONE_3RD*(pG->U[kp][jp-2][ip-1].M3 +
                                          pG->U[kp][jp-1][ip-2].M3 +
                                          pG->U[kp][jp-2][ip-2].M3);
      
      pG->U[kp][jp-1][ip+1].d  = ONE_3RD*(pG->U[kp][jp-2][ip+1].d  +
                                          pG->U[kp][jp-1][ip+2].d  +
                                          pG->U[kp][jp-2][ip+2].d );
      pG->U[kp][jp-1][ip+1].M1 = ONE_3RD*(pG->U[kp][jp-2][ip+1].M1 +
                                          pG->U[kp][jp-1][ip+2].M1 +
                                          pG->U[kp][jp-2][ip+2].M1);
      pG->U[kp][jp-1][ip+1].M2 = ONE_3RD*(pG->U[kp][jp-2][ip+1].M2 +
                                          pG->U[kp][jp-1][ip+2].M2 +
                                          pG->U[kp][jp-2][ip+2].M2);
      pG->U[kp][jp-1][ip+1].M3 = ONE_3RD*(pG->U[kp][jp-2][ip+1].M3 +
                                          pG->U[kp][jp-1][ip+2].M3 +
                                          pG->U[kp][jp-2][ip+2].M3);
      
      pG->U[kp][jp+1][ip+1].d  = ONE_3RD*(pG->U[kp][jp+2][ip+1].d  +
                                          pG->U[kp][jp+1][ip+2].d  +
                                          pG->U[kp][jp+2][ip+2].d );
      pG->U[kp][jp+1][ip+1].M1 = ONE_3RD*(pG->U[kp][jp+2][ip+1].M1 +
                                          pG->U[kp][jp+1][ip+2].M1 +
                                          pG->U[kp][jp+2][ip+2].M1);
      pG->U[kp][jp+1][ip+1].M2 = ONE_3RD*(pG->U[kp][jp+2][ip+1].M2 +
                                          pG->U[kp][jp+1][ip+2].M2 +
                                          pG->U[kp][jp+2][ip+2].M2);
      pG->U[kp][jp+1][ip+1].M3 = ONE_3RD*(pG->U[kp][jp+2][ip+1].M3 +
                                          pG->U[kp][jp+1][ip+2].M3 +
                                          pG->U[kp][jp+2][ip+2].M3);
      
      pG->U[kp][jp+1][ip-1].d  = ONE_3RD*(pG->U[kp][jp+2][ip-1].d  +
                                          pG->U[kp][jp+1][ip-2].d  +
                                          pG->U[kp][jp+2][ip-2].d );
      pG->U[kp][jp+1][ip-1].M1 = ONE_3RD*(pG->U[kp][jp+2][ip-1].M1 +
                                          pG->U[kp][jp+1][ip-2].M1 +
                                          pG->U[kp][jp+2][ip-2].M1);
      pG->U[kp][jp+1][ip-1].M2 = ONE_3RD*(pG->U[kp][jp+2][ip-1].M2 +
                                          pG->U[kp][jp+1][ip-2].M2 +
                                          pG->U[kp][jp+2][ip-2].M2);
      pG->U[kp][jp+1][ip-1].M3 = ONE_3RD*(pG->U[kp][jp+2][ip-1].M3 +
                                          pG->U[kp][jp+1][ip-2].M3 +
                                          pG->U[kp][jp+2][ip-2].M3);
      
      /* top plane outflow */
      pG->U[kp+1][jp-1][ip].d  = ONE_3RD*(pG->U[kp+2][jp-1][ip].d  +
                                          pG->U[kp+1][jp-2][ip].d  +
                                          pG->U[kp+2][jp-2][ip].d );
      pG->U[kp+1][jp-1][ip].M1 = ONE_3RD*(pG->U[kp+2][jp-1][ip].M1 +
                                          pG->U[kp+1][jp-2][ip].M1 +
                                          pG->U[kp+2][jp-2][ip].M1);
      pG->U[kp+1][jp-1][ip].M2 = ONE_3RD*(pG->U[kp+2][jp-1][ip].M2 +
                                          pG->U[kp+1][jp-2][ip].M2 +
                                          pG->U[kp+2][jp-2][ip].M2);
      pG->U[kp+1][jp-1][ip].M3 = ONE_3RD*(pG->U[kp+2][jp-1][ip].M3 +
                                          pG->U[kp+1][jp-2][ip].M3 +
                                          pG->U[kp+2][jp-2][ip].M3);
      
      pG->U[kp+1][jp][ip+1].d  = ONE_3RD*(pG->U[kp+2][jp][ip+1].d  +
                                          pG->U[kp+1][jp][ip+2].d  +
                                          pG->U[kp+2][jp][ip+2].d );
      pG->U[kp+1][jp][ip+1].M1 = ONE_3RD*(pG->U[kp+2][jp][ip+1].M1 +
                                          pG->U[kp+1][jp][ip+2].M1 +
                                          pG->U[kp+2][jp][ip+2].M1);
      pG->U[kp+1][jp][ip+1].M2 = ONE_3RD*(pG->U[kp+2][jp][ip+1].M2 +
                                          pG->U[kp+1][jp][ip+2].M2 +
                                          pG->U[kp+2][jp][ip+2].M2);
      pG->U[kp+1][jp][ip+1].M3 = ONE_3RD*(pG->U[kp+2][jp][ip+1].M3 +
                                          pG->U[kp+1][jp][ip+2].M3 +
                                          pG->U[kp+2][jp][ip+2].M3);
      
      pG->U[kp+1][jp+1][ip].d  = ONE_3RD*(pG->U[kp+2][jp+1][ip].d  +
                                          pG->U[kp+1][jp+2][ip].d  +
                                          pG->U[kp+2][jp+2][ip].d );
      pG->U[kp+1][jp+1][ip].M1 = ONE_3RD*(pG->U[kp+2][jp+1][ip].M1 +
                                          pG->U[kp+1][jp+2][ip].M1 +
                                          pG->U[kp+2][jp+2][ip].M1);
      pG->U[kp+1][jp+1][ip].M2 = ONE_3RD*(pG->U[kp+2][jp+1][ip].M2 +
                                          pG->U[kp+1][jp+2][ip].M2 +
                                          pG->U[kp+2][jp+2][ip].M2);
      pG->U[kp+1][jp+1][ip].M3 = ONE_3RD*(pG->U[kp+2][jp+1][ip].M3 +
                                          pG->U[kp+1][jp+2][ip].M3 +
                                          pG->U[kp+2][jp+2][ip].M3);
      
      pG->U[kp+1][jp][ip-1].d  = ONE_3RD*(pG->U[kp+2][jp][ip-1].d  +
                                          pG->U[kp+1][jp][ip-2].d  +
                                          pG->U[kp+2][jp][ip-2].d );
      pG->U[kp+1][jp][ip-1].M1 = ONE_3RD*(pG->U[kp+2][jp][ip-1].M1 +
                                          pG->U[kp+1][jp][ip-2].M1 +
                                          pG->U[kp+2][jp][ip-2].M1);
      pG->U[kp+1][jp][ip-1].M2 = ONE_3RD*(pG->U[kp+2][jp][ip-1].M2 +
                                          pG->U[kp+1][jp][ip-2].M2 +
                                          pG->U[kp+2][jp][ip-2].M2);
      pG->U[kp+1][jp][ip-1].M3 = ONE_3RD*(pG->U[kp+2][jp][ip-1].M3 +
                                          pG->U[kp+1][jp][ip-2].M3 +
                                          pG->U[kp+2][jp][ip-2].M3);
      
// need to set pG at ip, jp, kp
}
#endif
#endif /* STAR_PARTICLE */
