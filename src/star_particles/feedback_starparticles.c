#include "../copyright.h"

/*=============================================================================
 * FILE: feedback_starparticles.c
 *
 * PURPOSE: find and reset feedback regions
 *
 * CONTAINS PUBLIC FUNCTIONS:
 *   modify_feedback_region_starparticles();
 *
 * HISTORY:
 *   Written by Chang-Goo Kim, Mar. 2014
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

//#define SEDOV
#define NFEEDBACK_VAR 7

typedef struct FeedbackRegion_s{
  Real Mtot;
  Real M2tot;
  Real Vol;
  Real mom1;
  Real mom2;
  Real mom3;
  Real eint;
}FeedbackRegionS;

#ifdef MPI_PARALLEL
static MPI_Datatype FeedbackRegionType;
#endif

static int iHII, iSN, ifeedback;
static Real rf_min, rf_max=50.; // min and max of feedback region
static Real Rs0=3.13, TII = 64.856; // T=10^4K in the code unit
static Real Q490=3.981, eta=4./7.,c2=8.05,alphaB=7.52e-7; //
static Real Rrad0=23.7, Ush0=188;
static Real hLmax;

void compute_mean_feedback_region(GridS *pG, StarParS *pStar, Real rfeedback);
void assign_mean_feedback_region(GridS *pG, StarParS *pStar, Real rfeedback);
void assign_HII_feedback(GridS *pG, int starpar_id, Real Ri);
void assign_SN_feedback(GridS *pG, StarParS *pGStar);
void get_Sedov_Taylor(Real xi, Real *alpha, Real *v, Real *p);
void remove_starpar(GridS *pG, int starpar_id);
void push_feedback_to_local_list(GridS *pG, StarParS *pStar);
Real compute_HII_radius(GridS *pG, StarParS *pStar);
void compute_SN_radius(GridS *pG, StarParS *pStar);

/*=========================== PUBLIC FUNCTIONS ===============================*/
/*----------------------------------------------------------------------------*/
void feedback_starparticles_init(MeshS *pM){
  GridS *pG;
  int nl,nd,size1=1,size2=1,size3=1,Nx1,Nx2,Nx3;
  int i,j,k;
#ifdef MPI_PARALLEL
  int ierr;
  ierr = MPI_Type_contiguous(NFEEDBACK_VAR, MPI_DOUBLE, &FeedbackRegionType);
  if(ierr != MPI_SUCCESS) ath_error("[feedback_starparticles]: Error calling MPI_Type_contiguous\n");
  ierr = MPI_Type_commit(&FeedbackRegionType);
  if(ierr != MPI_SUCCESS) ath_error("[feedback_starparticles]: Error calling MPI_Type_commit\n");
#endif

  hLmax=MIN(MIN(pM->RootMaxX[0],pM->RootMaxX[1]),pM->RootMaxX[2]);
  rf_min=pM->dx[0]*NHII;
  rf_max=hLmax;

/* Cycle over all Grids on this processor to find maximum Nx1, Nx2, Nx3 */
  for (nl=0; nl<(pM->NLevels); nl++){
    for (nd=0; nd<(pM->DomainsPerLevel[nl]); nd++){
      if (pM->Domain[nl][nd].Grid != NULL) {
        pG=pM->Domain[nl][nd].Grid;
        size1 = pG->Nx[0]+2*nghost;
        size2 = pG->Nx[1]+2*nghost;
        size3 = pG->Nx[2]+2*nghost;
        if ((pG->sfr_ratio = (Real***)calloc_3d_array(size3,size2,size1, sizeof(Real))) == NULL)
          goto on_error;

        for (k=pG->ks;k<=pG->ke;k++){
          for (j=pG->js;j<=pG->je;j++){
            for (i=pG->is;i<=pG->ie;i++){
              pG->sfr_ratio[k][j][i]=1.0;
            }
          }
        }

/*
        synchro_starparticles(&(pM->Domain[nl][nd]));
        feedback_starparticles_assign(&(pM->Domain[nl][nd]));
*/
      }
    }
  }

  if(tHII < 0.0) {
    if(myID_Comm_world==0) ath_pout(0,"NO HII REGION FEEDBACK\n"); 
    iHII = 0;
  } else {
    if(myID_Comm_world==0) ath_pout(0,"HII REGION FEEDBACK for t<%g Myr\n",tHII); 
    iHII = 1;
  }

  if(tSN < 0.0) {
    if(myID_Comm_world==0) ath_pout(0,"NO SN FEEDBACK\n");
    iSN = 0;
  } else {
    if(myID_Comm_world==0) ath_pout(0,"SN FEEDBACK after t>%g Myr\n",tSN); 
    iSN = 1;
  }

  ifeedback = iHII + 2*iSN;
/*
  for (nl=0; nl<(pM->NLevels); nl++){
    for (nd=0; nd<(pM->DomainsPerLevel[nl]); nd++){
      if (pM->Domain[nl][nd].Grid != NULL) {
        pG=pM->Domain[nl][nd].Grid;

        synchro_starparticles(&(pM->Domain[nl][nd]));
        feedback_starparticles(&(pM->Domain[nl][nd]));
      }
    }
  }
*/
  return;

  on_error:
    feedback_starparticles_destruct();
    ath_error("[feedback_starparticle_init]: malloc returned a NULL pointer\n");


  return;
}

void feedback_starparticles_destruct(void)
{
  return;
}

void feedback_starparticles_assign(DomainS *pD){
  Real Vsh;
  GridS *pG = pD->Grid;
  StarParListS *pGstars=NULL;
  StarParS *pGStar=NULL;

  if(ifeedback != 0){
    pGstars = pG->Gstars;
    while (pGstars) {
      pGStar = &(pGstars->starpar);
/* CHECK: are there any new star particles? */
/* for newly formed star particle, calculate initial Stromgren radius, iteratively */
      if(iHII && pGStar->age == 0){
        pGStar->radius = 0.0;
        Vsh=compute_HII_radius(pG,pGStar);
        if(myID_Comm_world==0)
          printf("[HII feedback init] age=%g, navg=%g, radius=%g, Vol=%g\n",
                pGStar->age,pGStar->navg,pGStar->radius,pGStar->Vol);
        push_feedback_to_local_list(pG,pGStar);
        assign_mean_feedback_region(pG,pGStar,pGStar->radius);
      }
      pGstars = pGstars->next;
    }
  } 

  return;
}

void feedback_starparticles(DomainS *pD){
  GridS *pG = pD->Grid;
  StarParListS *pGstars=NULL;
  StarParS *pGStar=NULL;
  Real age,Vsh;
  int is=pG->is, ie=pG->ie, js=pG->js, je=pG->je, ks=pG->ks, ke=pG->ke;
  int i,j,k;

  if(ifeedback != 0){
    for (k=pG->ks;k<=pG->ke;k++){
      for (j=pG->js;j<=pG->je;j++){
        for (i=pG->is;i<=pG->ie;i++){
          pG->sfr_ratio[k][j][i]=1.0;
        }
      }
    }

    pGstars = pG->Gstars;
    while (pGstars) {
      pGStar = &(pGstars->starpar);
      age =(pGStar->age*pG->units.Myr);
      if((ifeedback == 3 && age < tHII) || (ifeedback == 1 && age < tHII)){
        Vsh = compute_HII_radius(pG,pGStar);
/*
        if(Vsh > c2){
          tHII = 0.0;
          break;
        }
*/
        if(pGStar->radius > rf_min) assign_HII_feedback(pG,pGStar->id,pGStar->radius);
      } // HII feedback

      if((ifeedback == 3 && age >= tSN) || (ifeedback == 2 && age >= tSN)){
        compute_SN_radius(pG,pGStar);
        assign_SN_feedback(pG,pGStar);
        tSN += 1.;
        if(tSN > 40.) 
        remove_starpar(pG,pGStar->id);
      } // SN feedback
      push_feedback_to_local_list(pG,pGStar);
      pGstars = pGstars->next;
    } // pGstar
  } // ifeedback

  return;
}

Real compute_HII_radius(GridS *pG, StarParS *pStar){
  Real age;
  Real msp3,Q49,namb2,Q1,Q2,Q3,R1,R2,R3,V0,Vsh;
  Real radius;
  int count=0;

  msp3=pStar->m*pG->units.Msun/1.e3;
  age=pStar->age*pG->units.Myr;
  Q49=age < 3.0 ? Q490*msp3: Q490*msp3*pow(age/3.0,-4.7);
  if(pStar->radius == 0.0){
    compute_mean_feedback_region(pG,pStar,rf_min);
    pStar->radius = rf_min;
  }

  Q1=alphaB*pStar->Vol*pStar->n2avg;
  Q2=Q1;
  R1=pStar->radius;
  R2=pStar->radius;
  V0=pStar->Vol;
  while(Q1 > Q49 && R1 > rf_min){
    R2 = R1;
    Q2 = Q1;
    R1 = R1-0.3*pG->dx1;
    compute_mean_feedback_region(pG,pStar,R1);
    Q1=alphaB*pStar->Vol*pStar->n2avg;
    if(myID_Comm_world == 0) printf("[HII feedback] Q49=%g, Q1=%g, R1=%g, Vol=%g, n2avg=%g\n",Q49,Q1,R1,pStar->Vol,pStar->n2avg);
  }
  while(Q2 <= Q49 && R2 <= rf_max){
    R1 = R2;
    Q1 = Q2;
    R2 = R2+0.3*pG->dx1;
    compute_mean_feedback_region(pG,pStar,R2);
    Q2=alphaB*pStar->Vol*pStar->n2avg;
    if(myID_Comm_world == 0) printf("[HII feedback] Q49=%g, Q2=%g, R2=%g, Vol=%g, n2avg=%g\n",Q49,Q2,R2,pStar->Vol,pStar->n2avg);
  }

  R3 = R2-((Q2-Q49)*(R2-R1)/(Q2-Q1));
  compute_mean_feedback_region(pG,pStar,R3);
  Q3 = alphaB*pStar->Vol*pStar->n2avg;
  if(myID_Comm_world == 0) printf("[HII feedback] final: Q49=%g, (R1,Q1)=(%g,%g), (R3,Q3)=(%g,%g), (R2,Q2)=(%g,%g)\n",Q49,R1,Q1,R3,Q3,R2,Q2);
  if(Q3 < Q49)
    R3 = R3 + (R2-R3)*(Q49-Q3)/(Q2-Q3);
  else
    R3 = R1 + (R3-R1)*(Q49-Q1)/(Q3-Q1);
  Vsh = (R3-pStar->radius)/pG->dt;
  pStar->radius = MIN(R3,rf_max);

  if(myID_Comm_world == 0) printf("[HII feedback] Q49=%g, Q3=%g, Vsh=%g\n",Q49,Q3,Vsh);
  if(myID_Comm_world == 0) printf("[HII feedback] age=%g, navg=%g, radius=%g, Vol=%g, count=%d\n",
                pStar->age,pStar->navg,pStar->radius,pStar->Vol,count);
  return Vsh;
}

void compute_SN_radius(GridS *pG, StarParS *pStar){
  Real Rrad,T9,R1,R2,dR=pG->dx1;
  int iresolved=0,count=0;

  rf_min = pG->dx1*(NSN+1);
  pStar->radius = pG->dx1*(NSN+1);
  while(pStar->radius >= rf_min && count < 5){
    compute_mean_feedback_region(pG,pStar,pStar->radius);
    Rrad = 0.3*Rrad0*pow(pStar->navg,-0.42);
    T9 = 1.0/(pStar->navg*CUBE(pStar->radius/3.3));
#ifdef SEDOV
    if(pStar->radius > Rrad){
      pStar->radius -= pG->dx1;
    } else {
      if((Rrad-pStar->radius) < pG->dx1) {
        iresolved=1;
        break;
      }
      pStar->radius += pG->dx1;
    }
    if(myID_Comm_world == 0) printf("[SN feedback] age=%g, navg=%g, radius=%g, Rrad=%g, ires=%d\n",
          pStar->age,pStar->navg,pStar->radius,Rrad,iresolved);
    count++;
#else
    if(fabs(pStar->radius-Rrad) < pG->dx1 || Rrad < rf_min) break;
    pStar->radius = Rrad;
    count++;
#endif
  }

  if(myID_Comm_world == 0) printf("[SN feedback] age=%g, navg=%g, radius=%g, Rrad=%g, ires=%d\n",
          pStar->age,pStar->navg,pStar->radius,Rrad,iresolved);

  return;
}

void compute_mean_feedback_region(GridS *pG, StarParS *pStar, Real rfeedback) {
  StarParListS *pGList=NULL;
  StarParS *pGStar=NULL;
  int is=pG->is, ie=pG->ie, js=pG->js, je=pG->je, ks=pG->ks, ke=pG->ke;
  int i,j,k;
  int ip,jp,kp;
  int nfeedback=(int)(rfeedback/pG->dx1)+1;
  Real x1,x2,x3,r,dVol=pG->dx1*pG->dx2*pG->dx3;
  
  FeedbackRegionS sendbuf,*recvbuf;
/* Notet that nGstars are set for total number of stars have ever been created */
  int nProc = 1;
#ifdef MPI_PARALLEL
  int mpierr;

  MPI_Comm_size(MPI_COMM_WORLD, &nProc);
#endif

    /* Create local buffer to store feedback region information for passing */
/*
    sendbuf = (FeedbackRegionS *)calloc_1d_array(1,sizeof(FeedbackRegionS));
    if (sendbuf == NULL)
      ath_error("[synchro_starparticles]: Error calling calloc_1d_array!\n");
*/
    recvbuf = (FeedbackRegionS *)calloc_1d_array(nProc,sizeof(FeedbackRegionS));
    if (recvbuf == NULL)
      ath_error("[synchro_starparticles]: Error calling calloc_1d_array!\n");

    /* Initialize sendbuf */
    sendbuf.Mtot = 0.0;
    sendbuf.M2tot = 0.0;
    sendbuf.mom1 = 0.0;
    sendbuf.mom2 = 0.0;
    sendbuf.mom3 = 0.0;
    sendbuf.Vol = 0.0;
    sendbuf.eint = 0.0;

    pGList = pG->Gstars_fda;
    while (pGList) {
      pGStar = &(pGList->starpar);
      if(pGStar->id == pStar->id){
      cc_ijk(pG,pGStar->x1,pGStar->x2,pGStar->x3,&ip,&jp,&kp);

/* Loop over all possible feedback regions. 
 * Further restriction will be checked with spherical radius */
      for (k=kp-nfeedback; k<=kp+nfeedback; k++) {
        for (j=jp-nfeedback; j<=jp+nfeedback; j++) {
          for (i=ip-nfeedback; i<=ip+nfeedback; i++) {
            cc_pos(pG,i,j,k,&x1,&x2,&x3);
/* Calculate distance from the star particle. */
            r = sqrt(SQR(x1-pGStar->x1)+SQR(x2-pGStar->x2)+SQR(x3-pGStar->x3));
//            r = sqrt(SQR(i-ip)+SQR(j-jp)+SQR(k-kp))*pG->dx1;
/* Check whether these zones are in grid and inside feedback radius. */
            if ((i >= is) && (i <= ie) &&
                (j >= js) && (j <= je) &&
                (k >= ks) && (k <= ke) && 
                (r <= rfeedback)) {
              sendbuf.Mtot += pG->U[k][j][i].d*dVol;
              sendbuf.M2tot += SQR(pG->U[k][j][i].d)*dVol;
              sendbuf.Vol  += dVol;
              sendbuf.mom1 += pG->U[k][j][i].M1*dVol;
              sendbuf.mom2 += pG->U[k][j][i].M2*dVol;
              sendbuf.mom3 += pG->U[k][j][i].M3*dVol;
#if defined (SHEARING_BOX) && !defined (FARGO)
              sendbuf.mom2 += qshear*Omega_0*x1*pG->U[k][j][i].d*dVol;
#endif
#ifndef BAROTROPIC
              sendbuf.eint += pG->U[k][j][i].E -0.5*(SQR(pG->U[k][j][i].M1)
                                                   + SQR(pG->U[k][j][i].M2)
                                                   + SQR(pG->U[k][j][i].M3))
                                                   /pG->U[k][j][i].d;
#endif
            }
          }
        }
      }
      } // check starpar_id
      pGList = pGList->next;
    } // while pList

#ifdef MPI_PARALLEL
    /* Communicate local star particle lists all-to-all */
    mpierr = MPI_Allgather(&(sendbuf), 1, FeedbackRegionType, 
                           recvbuf, 1, FeedbackRegionType, MPI_COMM_WORLD);

    if (mpierr != MPI_SUCCESS)
      ath_error("[feedback_starpraticles]: Error calling MPI_Allgather!\n");

    sendbuf.Mtot = 0.0;
    sendbuf.M2tot = 0.0;
    sendbuf.mom1 = 0.0;
    sendbuf.mom2 = 0.0;
    sendbuf.mom3 = 0.0;
    sendbuf.Vol = 0.0;
    sendbuf.eint = 0.0;


    for(i=0;i<nProc; i++){
      sendbuf.Mtot += recvbuf[i].Mtot;
      sendbuf.M2tot += recvbuf[i].M2tot;
      sendbuf.Vol += recvbuf[i].Vol;
      sendbuf.mom1 += recvbuf[i].mom1;
      sendbuf.mom2 += recvbuf[i].mom2;
      sendbuf.mom3 += recvbuf[i].mom3;
      sendbuf.eint += recvbuf[i].eint;
//      printf("[send/recv] time=%g, proc id=%d, %d, star id=%d, Mtot=%g, Mtot=%g\n",
//          pG->time,myID_Comm_world,i,pStar->id,recvbuf[i].Mtot,sendbuf.Mtot);
    }
#endif

    pStar->navg = sendbuf.Mtot/sendbuf.Vol;
    pStar->n2avg = sendbuf.M2tot/sendbuf.Vol;
    pStar->v1avg = sendbuf.mom1/sendbuf.Mtot;
    pStar->v2avg = sendbuf.mom2/sendbuf.Mtot;
    pStar->v3avg = sendbuf.mom3/sendbuf.Mtot;
    pStar->eavg = sendbuf.eint/sendbuf.Vol;
    pStar->Vol = sendbuf.Vol;

  return;
}

void assign_mean_feedback_region(GridS *pG, StarParS *pStar, Real rfeedback) {
  StarParListS *pGList=NULL;
  StarParS *pGStar=NULL;
  int is=pG->is, ie=pG->ie, js=pG->js, je=pG->je, ks=pG->ks, ke=pG->ke;
  int i,j,k;
  int ip,jp,kp;
  int nfeedback=(int)(rfeedback/pG->dx1)+1;
  Real x1,x2,x3,r,dVol=pG->dx1*pG->dx2*pG->dx3;
  
  pGList = pG->Gstars_fda;
  while (pGList) {
    pGStar = &(pGList->starpar);
    if(pGStar->id == pStar->id){
      cc_ijk(pG,pGStar->x1,pGStar->x2,pGStar->x3,&ip,&jp,&kp);

      for (k=kp-nfeedback; k<=kp+nfeedback; k++) {
        for (j=jp-nfeedback; j<=jp+nfeedback; j++) {
          for (i=ip-nfeedback; i<=ip+nfeedback; i++) {
            cc_pos(pG,i,j,k,&x1,&x2,&x3);
/* Calculate distance from the star particle. */
            r = sqrt(SQR(x1-pStar->x1)+SQR(x2-pStar->x2)+SQR(x3-pStar->x3));
//            r = sqrt(SQR(i-ip)+SQR(j-jp)+SQR(k-kp))*pG->dx1;
/* Check whether these zones are in grid and inside feedback radius. */
            if ((i >= is) && (i <= ie) &&
                (j >= js) && (j <= je) &&
                (k >= ks) && (k <= ke) && 
                (r <= rfeedback)) {
              pG->U[k][j][i].d = pStar->navg;

              pG->U[k][j][i].M1 = pStar->v1avg*pStar->navg;
              pG->U[k][j][i].M2 = pStar->v2avg*pStar->navg;
#if defined (SHEARING_BOX) && !defined (FARGO)
              pG->U[k][j][i].M2 -= qshear*Omega_0*x1*pG->U[k][j][i].d;
#endif
              pG->U[k][j][i].M3 = pStar->v3avg*pStar->navg;
#ifndef BAROTROPIC
              pG->U[k][j][i].E = TII*pG->U[k][j][i].d/Gamma_1+0.5/pG->U[k][j][i].d*
                                 (SQR(pG->U[k][j][i].M1)+SQR(pG->U[k][j][i].M2)
                                 +SQR(pG->U[k][j][i].M3));
/*
              pG->U[k][j][i].E = pStar->eavg+0.5/pG->U[k][j][i].d*
                                 (SQR(pG->U[k][j][i].M1)+SQR(pG->U[k][j][i].M2)
                                 +SQR(pG->U[k][j][i].M3));
*/
#endif
              pG->sfr_ratio[k][j][i]=100.0;
            }
          }
        }
      }
    }
    pGList = pGList->next;
  }

  return;
}

void assign_HII_feedback(GridS *pG, int starpar_id, Real Ri){
  StarParListS *pList=NULL;
  StarParS *pStar=NULL;
  int is=pG->is, ie=pG->ie, js=pG->js, je=pG->je, ks=pG->ks, ke=pG->ke;
  int i,j,k;
  int ip,jp,kp;
  int nf,ns;
  Real x1,x2,x3,r;
  Real dv1,dv2,dv3;
  Real Vi,psh,rf;
  
  pList = pG->Gstars_fda;
  while(pList){
    pStar=&(pList->starpar);
    if(pStar->id == starpar_id){
      cc_ijk(pG,pStar->x1,pStar->x2,pStar->x3,&ip,&jp,&kp);

// for log Q=49.6+log(M/10^3), T=10^4. See Draine (2011) Eqs. (15.3) (37.26)
// Rs in units of pc, c2 is in unit of km/s
//  Vi = 4*PI/3.*CUBE(Ri);
//  psh = (15./4.)*Ri/SQR(ts*(1+tau))*pG->dt*Vi/pStar->Vol; // for p(r)\propto psh*r^2
//  psh = (15./4.)*pStar->navg*Ri/SQR(ts*(1+tau))*pG->dt; // for p(r)\propto psh*r^2
//  psh = (25./4.)*pStar->navg*pG->dx1/SQR(ts*(1+tau))*pG->dt;
//    psh = (15./4.)*pStar->navg*Ri/SQR(ts*(1+tau))*tSN*pG->units.Myr; // put total momentum at once: too strong

      rf = Ri;
      nf = (int)(rf/pG->dx1)+1;

//      printf("Rs=%g, Ri=%g, rf=%g, psh=%g\n",Rs,Ri,rf,psh);
/* Loop over all possible feedback regions. 
 * Further restriction will be checked with spherical radius */
      for (k=kp-nf; k<=kp+nf; k++) {
      for (j=jp-nf; j<=jp+nf; j++) {
      for (i=ip-nf; i<=ip+nf; i++) {
        cc_pos(pG,i,j,k,&x1,&x2,&x3);
/* Calculate distance from the star particle. */
        r = sqrt(SQR(x1-pStar->x1)+SQR(x2-pStar->x2)+SQR(x3-pStar->x3));
/* Check whether these zones are in grid and inside feedback radius. */
        if ((i >= is) && (i <= ie) &&
            (j >= js) && (j <= je) &&
            (k >= ks) && (k <= ke) && 
            (r <= rf)) {
/* calculate feedback momentum */
/*
          dv1 = psh*SQR(r/rf)*(x1-pStar->x1)/r/pG->U[k][j][i].d;
          dv2 = psh*SQR(r/rf)*(x2-pStar->x2)/r/pG->U[k][j][i].d;
          dv3 = psh*SQR(r/rf)*(x3-pStar->x3)/r/pG->U[k][j][i].d;
*/

#ifndef BAROTROPIC
              pG->U[k][j][i].E = TII*pG->U[k][j][i].d/Gamma_1+0.5/pG->U[k][j][i].d*
                                   (SQR(pG->U[k][j][i].M1)+SQR(pG->U[k][j][i].M2)
                                   +SQR(pG->U[k][j][i].M3));
/*
          pG->U[k][j][i].E += 0.5*pG->U[k][j][i].d*(SQR(dv1) + SQR(dv2) +SQR(dv3))
                            + pG->U[k][j][i].M1*dv1 + pG->U[k][j][i].M2*dv2
                            + pG->U[k][j][i].M3*dv3;
*/
#endif
/*
          pG->U[k][j][i].M1 += pG->U[k][j][i].d*dv1;
          pG->U[k][j][i].M2 += pG->U[k][j][i].d*dv2;
          pG->U[k][j][i].M3 += pG->U[k][j][i].d*dv3;
*/
          pG->sfr_ratio[k][j][i]=100.0; // assign non-unity value to turn off cooling
        }
      }}}
    }
    pList=pList->next;
  }
  return;
}

void assign_SN_feedback(GridS *pG, StarParS *pGStar){
  StarParListS *pList=NULL;
  StarParS *pStar=NULL;
  int is=pG->is, ie=pG->ie, js=pG->js, je=pG->je, ks=pG->ks, ke=pG->ke;
  int i,j,k;
  int ip,jp,kp;
  int nSN,rind;
  int starpar_id=pGStar->id;
  Real navg=pGStar->navg, v1avg=pGStar->v1avg, v2avg=pGStar->v2avg, v3avg=pGStar->v3avg;
  Real eavg=pGStar->eavg, Vol=pGStar->Vol;
  Real Rrad=pGStar->radius;
  Real x1,x2,x3,r,dr;
  Real dv1,dv2,dv3,vr,Pr,rhor;
  Real rsh,Ush,trad,rho1,xi,rhomin=1.e-3;
  Real tanh,fr=1.0;

  
  pList = pG->Gstars_fda;
  while(pList) {
    pStar = &(pList->starpar);
    if(pStar->id == starpar_id){
      cc_ijk(pG,pStar->x1,pStar->x2,pStar->x3,&ip,&jp,&kp);

      rho1 = navg;
      rsh = Rrad;
      trad = 0.0057/pG->units.Myr*pow(rsh/10.,2.5)*sqrt(navg);
      Ush = 0.4*rsh/trad;
      nSN = (int)(rsh/pG->dx1)+1;

/* Loop over all possible feedback regions. 
 * Further restriction will be checked with spherical radius */
      for (k=kp-nSN; k<=kp+nSN; k++) {
      for (j=jp-nSN; j<=jp+nSN; j++) {
      for (i=ip-nSN; i<=ip+nSN; i++) {
        cc_pos(pG,i,j,k,&x1,&x2,&x3);
/* Calculate distance from the star particle. */
        r = sqrt(SQR(x1-pStar->x1)+SQR(x2-pStar->x2)+SQR(x3-pStar->x3));
/* Check whether these zones are in grid and inside feedback radius. */
        if ((i >= is) && (i <= ie) &&
            (j >= js) && (j <= je) &&
            (k >= ks) && (k <= ke) && 
            (r <= rsh)) {
/* smoothing as in Wise & Abel 2008 */
          xi=r/rsh;
//          tanh = exp(10.*(xi-1))-exp(-10.*(xi-1));
//          tanh = tanh/(exp(10.*(xi-1))+exp(-10.*(xi-1)));
//          fr= 1.28*(0.5-0.5*tanh);

#ifdef SEDOV
/* Using Sedov-Taylor Solution:
*/
          get_Sedov_Taylor(xi,&rhor,&vr,&Pr);

          rhor = MAX(rhor*rho1*(Gamma+1)/Gamma_1,rhomin)*fr;
          vr = vr*2.0/(Gamma+1)*Ush*xi*fr;
          Pr = Pr*2.0/(Gamma+1)*rho1*SQR(Ush)*SQR(xi)*fr;
          if(r > pG->dx1){
            dv1 = vr*(x1-pStar->x1)/r;
            dv2 = vr*(x2-pStar->x2)/r;
            dv3 = vr*(x3-pStar->x3)/r;
          } else {
            dv1 = 0.0;
            dv2 = 0.0;
            dv3 = 0.0;
          }
          pG->U[k][j][i].d = rhor;
          pG->U[k][j][i].M1 = pG->U[k][j][i].d*dv1;
          pG->U[k][j][i].M2 = pG->U[k][j][i].d*dv2;
          pG->U[k][j][i].M3 = pG->U[k][j][i].d*dv3;
#ifndef BAROTROPIC
          pG->U[k][j][i].E = (Pr/Gamma_1+eavg)+0.5*pG->U[k][j][i].d*(SQR(dv1) + SQR(dv2) +SQR(dv3));
#endif
#else
/* Using thermal energy alone:
*/
          pG->U[k][j][i].d = navg;
          pG->U[k][j][i].M1 = navg*v1avg;
          pG->U[k][j][i].M2 = navg*v2avg;
          pG->U[k][j][i].M3 = navg*v3avg;

#ifndef BAROTROPIC
          Pr = Gamma_1*3.4e-5/Vol*fr;
          Pr = Pr/(pG->units.Dcode*SQR(pG->units.Vcode));
          pG->sfr_ratio[k][j][i]=1.0; // assign non-unity value to turn off cooling
          pG->U[k][j][i].E = Pr/Gamma_1+0.5/pG->U[k][j][i].d*(SQR(pG->U[k][j][i].M1)
                           + SQR(pG->U[k][j][i].M2)+SQR(pG->U[k][j][i].M3));
#endif
#endif
        }
      }}}
    }
    pList=pList->next;
  }
  return;
}

void remove_starpar(GridS *pG, int starpar_id){
  StarParListS *pList=NULL; 
  StarParS *pStar= NULL;

  pList = pG->Lstars; 
  while(pList){
    pStar = &(pList->starpar);
    if(pStar->id == starpar_id) pStar->id = -1;
    pList = pList->next;
  }

  pList = pG->Gstars; 
  while(pList){
    pStar = &(pList->starpar);
    if(pStar->id == starpar_id) pStar->id = -1;
    pList = pList->next;
  }

  pList = pG->Gstars_fda;
  while(pList){
    pStar = &(pList->starpar);
    if(pStar->id == starpar_id) pStar->id = -1;
    pList = pList->next;
  }

  return;
}
void get_Sedov_Taylor(Real xi, Real *alpha, Real *v, Real *p){
  int ixi,ixip1;
  Real dxi;
  Real STsol[99*4] = {
   0.02000,   0.00001,   0.01698,   0.30620,
   0.03000,   0.00001,   0.02444,   0.30620,
   0.04000,   0.00001,   0.03225,   0.30620,
   0.05000,   0.00001,   0.04016,   0.30620,
   0.06000,   0.00001,   0.04811,   0.30620,
   0.07000,   0.00001,   0.05608,   0.30620,
   0.08000,   0.00001,   0.06406,   0.30620,
   0.09000,   0.00001,   0.07205,   0.30620,
   0.10000,   0.00001,   0.08004,   0.30620,
   0.11000,   0.00002,   0.08803,   0.30620,
   0.12000,   0.00002,   0.09603,   0.30620,
   0.13000,   0.00004,   0.10402,   0.30620,
   0.14000,   0.00005,   0.11202,   0.30620,
   0.15000,   0.00007,   0.12002,   0.30620,
   0.16000,   0.00009,   0.12802,   0.30620,
   0.17000,   0.00012,   0.13601,   0.30620,
   0.18000,   0.00015,   0.14401,   0.30620,
   0.19000,   0.00019,   0.15201,   0.30620,
   0.20000,   0.00025,   0.16001,   0.30621,
   0.21000,   0.00031,   0.16801,   0.30621,
   0.22000,   0.00038,   0.17601,   0.30621,
   0.23000,   0.00046,   0.18401,   0.30622,
   0.24000,   0.00056,   0.19201,   0.30623,
   0.25000,   0.00067,   0.20001,   0.30624,
   0.26000,   0.00080,   0.20802,   0.30625,
   0.27000,   0.00095,   0.21602,   0.30626,
   0.28000,   0.00112,   0.22402,   0.30628,
   0.29000,   0.00131,   0.23203,   0.30630,
   0.30000,   0.00152,   0.24003,   0.30632,
   0.31000,   0.00176,   0.24804,   0.30635,
   0.32000,   0.00204,   0.25605,   0.30639,
   0.33000,   0.00234,   0.26406,   0.30643,
   0.34000,   0.00267,   0.27207,   0.30648,
   0.35000,   0.00305,   0.28009,   0.30654,
   0.36000,   0.00346,   0.28811,   0.30662,
   0.37000,   0.00391,   0.29614,   0.30670,
   0.38000,   0.00442,   0.30416,   0.30679,
   0.39000,   0.00496,   0.31220,   0.30690,
   0.40000,   0.00557,   0.32024,   0.30703,
   0.41000,   0.00622,   0.32829,   0.30717,
   0.42000,   0.00694,   0.33634,   0.30734,
   0.43000,   0.00772,   0.34441,   0.30753,
   0.44000,   0.00856,   0.35249,   0.30775,
   0.45000,   0.00948,   0.36058,   0.30799,
   0.46000,   0.01048,   0.36868,   0.30827,
   0.47000,   0.01155,   0.37680,   0.30858,
   0.48000,   0.01271,   0.38493,   0.30893,
   0.49000,   0.01396,   0.39309,   0.30933,
   0.50000,   0.01531,   0.40127,   0.30977,
   0.51000,   0.01676,   0.40947,   0.31027,
   0.52000,   0.01832,   0.41770,   0.31082,
   0.53000,   0.02000,   0.42596,   0.31144,
   0.54000,   0.02179,   0.43426,   0.31212,
   0.55000,   0.02372,   0.44259,   0.31288,
   0.56000,   0.02578,   0.45097,   0.31373,
   0.57000,   0.02800,   0.45939,   0.31467,
   0.58000,   0.03036,   0.46786,   0.31570,
   0.59000,   0.03290,   0.47639,   0.31685,
   0.60000,   0.03561,   0.48498,   0.31811,
   0.61000,   0.03850,   0.49363,   0.31950,
   0.62000,   0.04160,   0.50236,   0.32104,
   0.63000,   0.04492,   0.51118,   0.32272,
   0.64000,   0.04846,   0.52008,   0.32458,
   0.65000,   0.05225,   0.52907,   0.32662,
   0.66000,   0.05631,   0.53817,   0.32886,
   0.67000,   0.06065,   0.54738,   0.33132,
   0.68000,   0.06529,   0.55672,   0.33402,
   0.69000,   0.07027,   0.56619,   0.33698,
   0.70000,   0.07561,   0.57581,   0.34023,
   0.71000,   0.08134,   0.58558,   0.34379,
   0.72000,   0.08749,   0.59552,   0.34769,
   0.73000,   0.09411,   0.60564,   0.35197,
   0.74000,   0.10123,   0.61595,   0.35666,
   0.75000,   0.10891,   0.62648,   0.36181,
   0.76000,   0.11720,   0.63723,   0.36746,
   0.77000,   0.12615,   0.64821,   0.37367,
   0.78000,   0.13585,   0.65945,   0.38048,
   0.79000,   0.14638,   0.67097,   0.38797,
   0.80000,   0.15782,   0.68277,   0.39620,
   0.81000,   0.17028,   0.69488,   0.40527,
   0.82000,   0.18388,   0.70731,   0.41526,
   0.83000,   0.19876,   0.72008,   0.42628,
   0.84000,   0.21509,   0.73320,   0.43845,
   0.85000,   0.23305,   0.74670,   0.45190,
   0.86000,   0.25285,   0.76059,   0.46679,
   0.87000,   0.27476,   0.77488,   0.48329,
   0.88000,   0.29906,   0.78958,   0.50161,
   0.89000,   0.32610,   0.80471,   0.52198,
   0.90000,   0.35629,   0.82028,   0.54466,
   0.91000,   0.39010,   0.83629,   0.56997,
   0.92000,   0.42811,   0.85274,   0.59826,
   0.93000,   0.47099,   0.86965,   0.62995,
   0.94000,   0.51954,   0.88700,   0.66552,
   0.95000,   0.57472,   0.90479,   0.70552,
   0.96000,   0.63768,   0.92302,   0.75063,
   0.97000,   0.70981,   0.94167,   0.80162,
   0.98000,   0.79280,   0.96073,   0.85941,
   0.99000,   0.88869,   0.98018,   0.92510,
   1.00000,   1.00000,   1.00000,   1.00000,
  };

  dxi=0.01;
  ixi = MAX((int)(xi/dxi)-2,0);
  ixip1 = MIN(ixi+1,98);

  *alpha = STsol[ixi*4+1]+(STsol[ixip1*4+1]-STsol[ixi*4+1])*(xi-STsol[ixi*4])/dxi;
  *v = STsol[ixi*4+2]+(STsol[ixip1*4+2]-STsol[ixi*4+2])*(xi-STsol[ixi*4])/dxi;
  *p = STsol[ixi*4+3]+(STsol[ixip1*4+3]-STsol[ixi*4+3])*(xi-STsol[ixi*4])/dxi;

  return;
}

void push_feedback_to_local_list(GridS *pG, StarParS *pStar){
  StarParListS *pList=NULL;
  StarParS *pLStar=NULL;

  pList = pG->Lstars;
  while(pList){
    pLStar = &(pList->starpar);
    if(pLStar->id == pStar->id){
      pLStar->navg = pStar->navg;
      pLStar->n2avg = pStar->n2avg;
      pLStar->v1avg = pStar->v1avg;
      pLStar->v2avg = pStar->v2avg;
      pLStar->v3avg = pStar->v3avg;
      pLStar->eavg = pStar->eavg;
      pLStar->Vol = pStar->Vol;
      pLStar->radius = pStar->radius;
    }
    pList = pList->next;
  }
  return;
}
#endif /* STAR_PARTICLE */
