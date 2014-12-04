#include "copyright.h"
/*==============================================================================
 * FILE: msa.c
 *
 * PURPOSE:  Problem generator for MSA test
 *   3D shearing box code. 
 *
 * Code must be configured using --enable-shearing-box
 *
 * REFERENCE: Kim, W.-T and Ostriker, E. C. (2001)
 * Wriiten by Chang-Goo Kim
 *============================================================================*/

#include <float.h>
#include <math.h>

#include <stdlib.h>
#include <string.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"

/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 *============================================================================*/
static Real ShearPot(const Real x1, const Real x2, const Real x3);
static Real VerticalGrav(const Real x1, const Real x2, const Real x3);

static Real hst_sigma(const GridS *pG, const int i, const int j, const int k);
static Real hst_dmax(const GridS *pG, const int i, const int j, const int k);

static void slab_mhd_bc_lower(GridS *pG);
static void slab_mhd_bc_upper(GridS *pG);
static void slab_grav_bc_lower(GridS *pG);
static void slab_grav_bc_upper(GridS *pG);

static void initialize(DomainS *pD);
static int nwx,nwy;
static Real beta, amp, Q, nJ, s0, cs, cs2, Gcons, scaleH, fz;
static Real kx, ky, Lx, Ly, Lz;
static Real dmax;

/*=========================== PUBLIC FUNCTIONS =================================
 * Contains the usual, plus:
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/* problem:  */

void problem(DomainS *pDomain)
{
  GridS *pG = pDomain->Grid;
  int is = pG->is, ie = pG->ie;
  int js = pG->js, je = pG->je;
  int ks = pG->ks, ke = pG->ke;
  int ixs,jxs,kxs,i,j,k;
  int iprob;
  Real x1,x2,x3;
  Real rd, rp, rvx, rvy, rvz, rbx, rby, rbz;
  Real B0,P0;
  Real time0,kxt;
#ifdef SELF_GRAVITY
  Real Gcons;
#endif
  Real ***drho=NULL, std_dev;
  Real zfact = 1.0;
  int nx1,nx2,nx3;

  double rval;

#ifdef MPI_PARALLEL
  Real my_max;
  int ierr;
#endif


  iprob = par_geti_def("problem","iprob",1);

  printf("initializing...\n");
  initialize(pDomain);
  printf("initializing PS...\n");
  if(iprob == 2){
    /* Get local grid size */
    nx1 = pG->Nx[0];
    nx2 = pG->Nx[1];
    nx3 = pG->Nx[2];
 
 
    if ((drho=(Real***)calloc_3d_array(nx3+2*nghost,nx2+2*nghost,nx1+2*nghost,sizeof(Real)))==NULL)
      ath_error("[problem]: Error allocating memory for vel pert\n");
 
    initializePS(pDomain,drho,&std_dev);
    printf("drho[ks][js][is]=%g,std_dev=%g\n",drho[ks][js][is],std_dev);

  }

  time0=par_getd_def("problem","time0",0.0);

  B0 = cs/sqrt(beta);
#ifndef BAROTROPIC
  P0 = cs2/Gamma;
#endif

  kxt = kx+qshear*Omega_0*ky*time0;

  pG->time=time0;

  dmax=1.e-30;
  printf("assigning initial condition...\n");
  for (k=ks-nghost; k<=ke+nghost; k++) {
  for (j=js; j<=je; j++) {
    for (i=is; i<=ie; i++) {
      cc_pos(pG,i,j,k,&x1,&x2,&x3);
      if(iprob == 1){
      /* eigen mode perturbation */
        rd  = 1.0+amp*cos(kxt*x1+ky*x2);
        rvx = amp*kx/ky*sin(kxt*x1+ky*x2);
        rvy = amp*sin(kxt*x1+ky*x2);
        rvz = 0.0;
        rp  = cs2*(rd-1.0);
 
        rbx = amp*nwy*cos(kxt*(x1-0.5*pG->dx1)+ky*x2);
        rby = -amp*nwx*cos(kxt*x1+ky*(x2-0.5*pG->dx2));
        rbz = 0.0;
      }

      if(iprob == 2){
        rd = 1.0+amp*(drho[k][j][i]/std_dev);
        rvx = 0.0;
        rvy = 0.0;
        rvz = 0.0;
        rp  = cs2*(rd-1.0);

        zfact = exp(-0.5*SQR(x3/scaleH*fz));
        rd = rd*zfact;
      }
      /* gaussian random field */
      pG->U[k][j][i].d  = rd;
      pG->U[k][j][i].M1 = rd*rvx;
      pG->U[k][j][i].M2 = rd*rvy;
#ifndef FARGO
      pG->U[k][j][i].M2 -= rd*(qshear*Omega_0*x1);
#endif
      pG->U[k][j][i].M3 = rd*rvz;
#ifdef ADIABATIC
      pG->U[k][j][i].E = (P0+rp)/Gamma_1*zfact
         + 0.5*(SQR(pG->U[k][j][i].M1) + SQR(pG->U[k][j][i].M2) 
         + SQR(pG->U[k][j][i].M3))/rd;
#endif

#ifdef MHD
        pG->B1i[k][j][i] = rbx;
        pG->B2i[k][j][i] = B0+rby;
        pG->B3i[k][j][i] = 0.0;

        if (i==ie) cc_pos(pG,ie+1,j,k,&x1,&x2,&x3);
        rbx = amp*nwy*cos(kx*(x1-0.5*pG->dx1)+ky*x2);
        if (j==je) cc_pos(pG,i,je+1,k,&x1,&x2,&x3);
        rby = -amp*nwx*cos(kx*x1+ky*(x2-0.5*pG->dx2));
        if (i==ie) pG->B1i[k][j][ie+1] = rbx;
        if (j==je) pG->B2i[k][je+1][i] = B0+rby;
        if (pG->Nx[2] > 1 && k==ke) pG->B3i[ke+1][j][i] = 0.0;
#endif /* MHD */
      dmax = MAX(dmax,pG->U[k][j][i].d);
    }
  }}
#ifdef MHD
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        pG->U[k][j][i].B1c = 0.5*(pG->B1i[k][j][i]+pG->B1i[k][j][i+1]);
        pG->U[k][j][i].B2c = 0.5*(pG->B2i[k][j][i]+pG->B2i[k][j+1][i]);
        if (pG->Nx[2] >1) pG->U[k][j][i].B3c = 0.5*(pG->B3i[k][j][i]+pG->B3i[k+1][j][i]); else pG->U[k][j][i].B3c =pG->B3i[k][j][i];
#ifdef ADIABATIC
        pG->U[k][j][i].E += 0.5*(SQR(pG->U[k][j][i].B1c)
         + SQR(pG->U[k][j][i].B2c) + SQR(pG->U[k][j][i].B3c));
#endif
      }
    }
  }
#endif /* MHD */

#ifdef MPI_PARALLEL
  my_max = dmax;
  ierr = MPI_Allreduce(&my_max, &dmax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#endif


  if(iprob == 2) free_3d_array(drho);
  printf("=== end of problem setting ===\n");
  return;
}

/*==============================================================================
 * PUBLIC PROBLEM USER FUNCTIONS:
 * problem_write_restart() - writes problem-specific user data to restart files
 * problem_read_restart()  - reads problem-specific user data from restart files
 * get_usr_expr()          - sets pointer to expression for special output data
 * get_usr_out_fun()       - returns a user defined output function pointer
 * Userwork_in_loop        - problem specific work IN     main loop
 * Userwork_after_loop     - problem specific work AFTER  main loop
 *----------------------------------------------------------------------------*/

void problem_write_restart(MeshS *pM, FILE *fp)
{
  return;
}

/*
 * 'problem_read_restart' must enroll gravity on restarts
 */

void problem_read_restart(MeshS *pM, FILE *fp)
{
  DomainS *pD=NULL;
  GridS *pG=NULL;
  int nl,nd;

  for (nl=0; nl<pM->NLevels; nl++){
    for (nd=0; nd<pM->DomainsPerLevel[nl]; nd++){
      if (pM->Domain[nl][nd].Grid != NULL) {
        pD = &(pM->Domain[nl][nd]);
        pG = pM->Domain[nl][nd].Grid;

        initialize(pD);
      }
    }
  }

  return;
}

ConsFun_t get_usr_expr(const char *expr)
{
  return NULL;
}

VOutFun_t get_usr_out_fun(const char *name){
  return NULL;
}

void Userwork_in_loop(MeshS *pM)
{
  DomainS *pD=NULL;
  GridS *pG=NULL;
  Real dVol;
  int i,j,k;
  int nl,nd;

#ifdef MPI_PARALLEL
  Real my_max;
  int ierr;
#endif

  for (nl=0; nl<pM->NLevels; nl++){
    for (nd=0; nd<pM->DomainsPerLevel[nl]; nd++){
      if (pM->Domain[nl][nd].Grid != NULL) {
        pD = &(pM->Domain[nl][nd]);
        pG = pM->Domain[nl][nd].Grid;

        dmax = 1.e-30;
	dVol = pG->dx1*pG->dx2*pG->dx3;
        for (k=pG->ks; k<=pG->ke; k++) {
          for (j=pG->js; j<=pG->je; j++) {
            for (i=pG->is; i<=pG->ie; i++) {
              dmax = MAX(dmax,pG->U[k][j][i].d);
            }
          }
        }
      }
    }
  }

#ifdef MPI_PARALLEL
  my_max = dmax;
  ierr = MPI_Allreduce(&my_max, &dmax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#endif

}

void Userwork_after_loop(MeshS *pM)
{
}

/*==============================================================================
 * PRIVATE FUNCTION:
 *============================================================================*/
void initialize(DomainS* pD){
  GridS *pG = pD->Grid;

  int bc_ix3,bc_ox3;
  if(pG->Nx[2] == 1) ShBoxCoord = xy; /* 2D xy-plane */

  Lx = pD->MaxX[0]-pD->MinX[0];
  Ly = pD->MaxX[1]-pD->MinX[1];
  Lz = pD->MaxX[2]-pD->MinX[2];

/* Read problem parameters. */
  Omega_0 = par_getd("problem","omega");
  qshear = par_getd("problem","qshear");

  amp = par_getd("problem","amp");

/* Read parameters for magnetic field */
  beta = par_getd("problem","beta"); 

/* Read parameters for self gravity */
  Q=par_getd("problem","Q");
  nJ= par_getd("problem","nJ");
  s0 = par_getd("problem","s0");

#ifdef SELF_GRAVITY
  fz = sqrt(PI*s0)*(1.0+sqrt(1.0+16.0/PI/s0))/4.0;
#else
  fz = 1.0;
#endif

  scaleH = sqrt(s0/2.0)/PI/nJ;
  cs=sqrt(4.0-2.0*qshear)/PI/nJ/Q;
  cs2=SQR(cs);
#ifdef ISOTHERMAL
  Iso_csound=cs;
  Iso_csound2=cs2;
#endif

#ifdef SELF_GRAVITY
  Gcons = nJ*cs2/sqrt(2*PI)/scaleH;
  grav_mean_rho = 0.0;
#ifndef SELF_GRAVITY_USING_FFT_DISK
  if(pD->Nx[2] >1) grav_mean_rho = 1.0;
#endif

/* Set gravity constant*/
  four_pi_G = 4.0*PI*Gcons;
#endif /* SELF_GRAVITY */

/* initialize wavenumbers, given input number of waves per L */
  nwx = par_geti_def("problem","nwx",-6);
  nwy = par_geti_def("problem","nwy",1);

  ky = nwy*2.0*PI;
  kx = nwx*2.0*PI;

/* enroll gravitational potential function */

  ShearingBoxPot = ShearPot;
  StaticGravPot = VerticalGrav;

/* history dump for linear perturbation amplitude. See Kim & Ostriker 2001 */
  dump_history_enroll(hst_sigma, "<sigma>");
  dump_history_enroll(hst_dmax, "dmax");

#ifdef SELF_GRAVITY_USING_FFT_DISK
  bc_ix3=par_geti("domain1","bc_ix3");
  bc_ox3=par_geti("domain1","bc_ox3");
  bvals_grav_fun(pD, left_x3, slab_grav_bc_lower);
  bvals_grav_fun(pD, right_x3, slab_grav_bc_upper);
  if(bc_ix3 == 2) bvals_mhd_fun(pD, left_x3, slab_mhd_bc_lower);
  if(bc_ox3 == 2) bvals_mhd_fun(pD, right_x3, slab_mhd_bc_upper); 
#endif

}


static Real ShearPot(const Real x1, const Real x2, const Real x3)
{
  Real phi=0.0;
#ifndef FARGO
  phi -= qshear*Omega_0*Omega_0*x1*x1;
#endif
  return phi;
}

static Real VerticalGrav(const Real x1, const Real x2, const Real x3)
{
  Real phi=0.0;
  phi += 2.0*(2.0-qshear)/SQR(Q)/s0*SQR(x3);

  return phi;
}

/* Volume average of SQUARED perturbed quantity. 
 * int_0^L q^2 cos^2(kx*x+ky*y)dxdydz = 0.5 LxLyLz q^2
 * To get q, one should take square root of history dumps.
 */
static Real hst_sigma(const GridS *pG, const int i, const int j, const int k)
{
  return 2*SQR(pG->U[k][j][i].d-1);
}

static Real hst_dmax(const GridS *pG, const int i, const int j, const int k)
{
  return dmax;
}

/*--------------------------------------------------------------------------*/
/*  Two functions to set BCs of U   begin here                              */
/*--------------------------------------------------------------------------*/
/* slab_mhd_bc_upper:
 *   Applies open bc in vertical direction to obtain values of MHD variables 
 *   in ghost zones at: 
 *
 *      j>je (2D)
 *      k>ke (3D)
 * 
 *   Switches for 1D, 2D, 3D
 *   logarithmic extrapolation for density and pressure
 *   vertical velocity and other variables are extrapolated with zero slope
*/
static void slab_mhd_bc_upper(GridS *pG)
{
  int i, j, k;

  int is = pG->is, ie = pG->ie;
  int js = pG->js, je = pG->je;
  int ks = pG->ks, ke = pG->ke;

  int iu, il; /* i-upper, i-lower  needed for 2D and 3D */
  int ju, jl; /* j-upper, j-lower  needed for 3D */

  int dim;

  PrimS W_1, W_2;
  Real T, Press;

  iu = ie + nghost;
  il = is - nghost;

  dim = 0;
  if(pG->Nx[0] > 1) dim++;
  if(pG->Nx[1] > 1) dim++;
  if(pG->Nx[2] > 1) dim++;


  switch(dim){
    case 1:
      ath_error("Slab BC only works with 2D or 3D\n");
      break;
    case 2:

/* 2D case -- apply BC in x2 direction */
#ifdef MHD
      ath_error("Slab BC only works with HD for 2D\n");
#endif /* MHD */

    for (j=1; j<=nghost; j++) {
      for (i=il; i<=iu; i++) {
/* get primitive variables at last two active zones */
        W_1 = Cons_to_Prim(&(pG->U[ks][je][i]));
        W_2 = Cons_to_Prim(&(pG->U[ks][je-1][i]));
/* constant logarithmic gradient for density if decreasing outward; constant density otherwise*/
        pG->U[ks][je+j][i].d = W_1.d;
        if ((W_1.d-W_2.d) < 0) pG->U[ks][je+j][i].d *= pow(W_1.d/W_2.d,j);
/* ZERO GRADIENT FOR VELOCITY (zero force) and outflow for V2*/
        pG->U[ks][je+j][i].M1 = W_1.V1*pG->U[ks][je+j][i].d;
        pG->U[ks][je+j][i].M2 = MAX(W_1.V2*pG->U[ks][je+j][i].d, 0.0);
        pG->U[ks][je+j][i].M3 = W_1.V3*pG->U[ks][je+j][i].d;
#ifndef BAROTROPIC
        T = W_1.P/W_1.d;
        Press = W_1.P;
        Press *= pow(W_1.P/W_2.P,j);
        pG->U[ks][je+j][i].E = T*pG->U[ks][je+j][i].d/Gamma_1 + 0.5/pG->U[ks][je+j][i].d
//        pG->U[ks][je+j][i].E = Press/Gamma_1 + 0.5/pG->U[ks][je+j][i].d
          * (SQR(pG->U[ks][je+j][i].M1) + SQR(pG->U[ks][je+j][i].M2) +SQR(pG->U[ks][je+j][i].M3));
#endif /* BAROTROPIC */
      }
    }
#ifdef MHD
    for (j=1; j<=nghost; j++) {
      for (i=il; i<=iu; i++) {
        pG->B1i[ks][je+j][i] = pG->B1i[ks][je][i];
        pG->B3i[ks][je+j][i] = pG->B3i[ks][je][i];
        pG->U[ks][je+j][i].B1c = pG->U[ks][je][i].B1c;
        pG->U[ks][je+j][i].B2c = pG->U[ks][je][i].B2c;
        pG->U[ks][je+j][i].B3c = pG->U[ks][je][i].B3c;
#ifndef BAROTROPIC
        pG->U[ks][je+j][i].E += 0.5*(SQR(pG->U[ks][je+j][i].B1c)
                             +  SQR(pG->U[ks][je+j][i].B2c)+SQR(pG->U[ks][je+j][i].B3c));
#endif
      }
    }
/* Note that k=ke+1 is not a boundary condition for the interface field B3i */
    for (j=2; j<=nghost; j++) {
      for (i=il; i<=iu; i++) {
        pG->B2i[ks][je+j][i] = pG->B2i[ks][je+1][i];
      }
    }
#endif /* MHD */

    break;
    case 3:
/* 3D case -- apply BC in x3 direction */
    ju = je + nghost;
    jl = js - nghost;
    for (k=1; k<=nghost; k++) {
      for (j=jl; j<=ju; j++) {
        for (i=il; i<=iu; i++) {
/* get primitive variables at last two active zones */
          W_1 = Cons_to_Prim(&(pG->U[ke  ][j][i]));
          W_2 = Cons_to_Prim(&(pG->U[ke-1][j][i]));
/* constant logarithmic gradient for density if decreasing outward; constant density otherwise*/
          pG->U[ke+k][j][i].d = W_1.d;
          if ((W_1.d-W_2.d) < 0) pG->U[ke+k][j][i].d *= pow(W_1.d/W_2.d,k);
/* ZERO GRADIENT FOR VELOCITY (zero force) and outflow for V2*/
          pG->U[ke+k][j][i].M1 = W_1.V1*pG->U[ke+k][j][i].d;
          pG->U[ke+k][j][i].M2 = W_1.V2*pG->U[ke+k][j][i].d;
          pG->U[ke+k][j][i].M3 = MAX(W_1.V3*pG->U[ke+k][j][i].d, 0.0);
#ifndef BAROTROPIC
          T = W_1.P/W_1.d;
          Press = W_1.P;
          Press *= pow(W_1.P/W_2.P,k);
          pG->U[ke+k][j][i].E = T*pG->U[ke+k][j][i].d/Gamma_1 + 0.5/pG->U[ke+k][j][i].d
//        pG->U[ks][je+j][i].E = Press/Gamma_1 + 0.5/pG->U[ks][je+j][i].d
          * (SQR(pG->U[ke+k][j][i].M1) + SQR(pG->U[ke+k][j][i].M2) +SQR(pG->U[ke+k][j][i].M3));
#endif /* BAROTROPIC */
        }
      }
    }
#ifdef MHD
    for (k=1; k<=nghost; k++) {
      for (j=jl; j<=ju; j++) {
        for (i=il; i<=iu; i++) {
          pG->B1i[ke+k][j][i] = pG->B1i[ke][j][i];
          pG->B2i[ke+k][j][i] = pG->B2i[ke][j][i];
          pG->U[ke+k][j][i].B1c = pG->U[ke][j][i].B1c;
          pG->U[ke+k][j][i].B2c = pG->U[ke][j][i].B2c;
          pG->U[ke+k][j][i].B3c = pG->U[ke][j][i].B3c;
#ifndef BAROTROPIC
          pG->U[ke+k][j][i].E += 0.5*(SQR(pG->U[ke+k][j][i].B1c)
                              +  SQR(pG->U[ke+k][j][i].B2c)+SQR(pG->U[ke+k][j][i].B3c));
#endif
        }
      }
    }
/* Note that k=ke+1 is not a boundary condition for the interface field B3i */
    for (k=2; k<=nghost; k++) {
      for (j=jl; j<=ju; j++) {
        for (i=il; i<=iu; i++) {
          pG->B3i[ke+k][j][i] = pG->B3i[ke+1][j][i];
        }
      }
    }
#endif /* MHD */
  }
  return;
}
/* slab_mhd_bc_lower:
 *   Applies open bc in vertical direction to obtain values of MHD variables 
 *   in ghost zones at: 
 *
 *      j<js (2D)
 *      k<ks (3D)
 * 
 *   Switches for 2D, 3D
 *   logarithmic extrapolation for density and pressure
 *   vertical velocity and other variables are extrapolated with zero slope
*/
static void slab_mhd_bc_lower(GridS *pG)
{
  int i, j, k;

  int is = pG->is, ie = pG->ie;
  int js = pG->js, je = pG->je;
  int ks = pG->ks, ke = pG->ke;

  int dim;

  int iu, il; /* i-upper, i-lower  needed for 2D and 3D */
  int ju, jl; /* j-upper, j-lower  needed for 3D */
  PrimS W_1, W_2;
  Real T,Press;

  iu = ie + nghost;
  il = is - nghost;

  dim = 0;
  if(pG->Nx[0] > 1) dim++;
  if(pG->Nx[1] > 1) dim++;
  if(pG->Nx[2] > 1) dim++;

  switch(dim){
    case 1:
      ath_error("Slab BC only works with 2D or 3D\n");
      break;
    case 2:

/* 2D case -- apply BC in x2 direction */
#ifdef MHD
      ath_error("Slab BC only works with HD for 2D\n");
#endif /* MHD */

    for (j=1; j<=nghost; j++) {
      for (i=il; i<=iu; i++) {
/* get primitive variables at last two active zones */
        W_1 = Cons_to_Prim(&(pG->U[ks][js][i]));
        W_2 = Cons_to_Prim(&(pG->U[ks][js+1][i]));
/* constant logarithmic gradient for density if decreasing outward; constant density otherwise*/
        pG->U[ks][js-j][i].d = W_1.d;
        if ((W_1.d-W_2.d) < 0) pG->U[ks][js-j][i].d *= pow(W_1.d/W_2.d,j);
/* ZERO GRADIENT FOR VELOCITY (zero force) and outflow for V2*/
        pG->U[ks][js-j][i].M1 = W_1.V1*pG->U[ks][js-j][i].d;
        pG->U[ks][js-j][i].M2 = MIN(W_1.V2*pG->U[ks][js-j][i].d, 0.0);
        pG->U[ks][js-j][i].M3 = W_1.V3*pG->U[ks][js-j][i].d;
#ifndef BAROTROPIC
        T = W_1.P/W_1.d;
        Press = W_1.P;
        Press *= pow(W_1.P/W_2.P,j);
        pG->U[ks][js-j][i].E = T*pG->U[ks][js-j][i].d/Gamma_1 + 0.5/pG->U[ks][js-j][i].d
//        pG->U[ks][js-j][i].E = Press/Gamma_1 + 0.5/pG->U[ks][js-j][i].d
          * (SQR(pG->U[ks][js-j][i].M1) + SQR(pG->U[ks][js-j][i].M2) +SQR(pG->U[ks][js-j][i].M3));
#endif /* BAROTROPIC */
      }
    }
#ifdef MHD
    for (j=1; j<=nghost; j++) {
      for (i=il; i<=iu; i++) {
        pG->B1i[ks][js-j][i]   = pG->B1i[ks][js][i];
        pG->B2i[ks][js-j][i]   = pG->B2i[ks][js][i];
        pG->B3i[ks][js-j][i]   = pG->B3i[ks][js][i];
        pG->U[ks][js-j][i].B1c = pG->U[ks][js][i].B1c;
        pG->U[ks][js-j][i].B2c = pG->U[ks][js][i].B2c;
        pG->U[ks][js-j][i].B3c = pG->U[ks][js][i].B3c;
#ifndef BAROTROPIC
        pG->U[ks][js-j][i].E += 0.5*(SQR(pG->U[ks][js-j][i].B1c)
                             +  SQR(pG->U[ks][js-j][i].B2c)+SQR(pG->U[ks][js-j][i].B3c));
#endif
      }
    }
#endif /* MHD */
    break;
    case 3:
/* 3D case -- apply BC in x3 direction */
    ju = je + nghost;
    jl = js - nghost;
    for (k=1; k<=nghost; k++) {
      for (j=jl; j<=ju; j++) {
        for (i=il; i<=iu; i++) {
/* get primitive variables at last two active zones */
          W_1 = Cons_to_Prim(&(pG->U[ks  ][j][i]));
          W_2 = Cons_to_Prim(&(pG->U[ks+1][j][i]));
/* constant logarithmic gradient for density if decreasing outward; constant density otherwise*/
          pG->U[ks-k][j][i].d = W_1.d;
          if ((W_1.d-W_2.d) < 0) pG->U[ks-k][j][i].d *= pow(W_1.d/W_2.d,k);
/* ZERO GRADIENT FOR VELOCITY (zero force) and outflow for V2*/
          pG->U[ks-k][j][i].M1 = W_1.V1*pG->U[ks-k][j][i].d;
          pG->U[ks-k][j][i].M2 = W_1.V2*pG->U[ks-k][j][i].d;
          pG->U[ks-k][j][i].M3 = MAX(W_1.V3*pG->U[ks-k][j][i].d, 0.0);
#ifndef BAROTROPIC
          T = W_1.P/W_1.d;
          Press = W_1.P;
          Press *= pow(W_1.P/W_2.P,k);
          pG->U[ks-k][j][i].E = T*pG->U[ks-k][j][i].d/Gamma_1 + 0.5/pG->U[ks-k][j][i].d
//        pG->U[ks][je+j][i].E = Press/Gamma_1 + 0.5/pG->U[ks][je+j][i].d
          * (SQR(pG->U[ks-k][j][i].M1) + SQR(pG->U[ks-k][j][i].M2) +SQR(pG->U[ks-k][j][i].M3));
#endif /* BAROTROPIC */
        }
      }
    }
#ifdef MHD
    for (k=1; k<=nghost; k++) {
      for (j=jl; j<=ju; j++) {
        for (i=il; i<=iu; i++) {
          pG->B1i[ks-k][j][i]   = pG->B1i[ks][j][i];
          pG->B2i[ks-k][j][i]   = pG->B2i[ks][j][i];
          pG->B3i[ks-k][j][i]   = pG->B3i[ks][j][i];
          pG->U[ks-k][j][i].B1c = pG->U[ks][j][i].B1c;
          pG->U[ks-k][j][i].B2c = pG->U[ks][j][i].B2c;
          pG->U[ks-k][j][i].B3c = pG->U[ks][j][i].B3c;
#ifndef BAROTROPIC
          pG->U[ks-k][j][i].E += 0.5*(SQR(pG->U[ks-k][j][i].B1c)
                              +  SQR(pG->U[ks-k][j][i].B2c)+SQR(pG->U[ks-k][j][i].B3c));
#endif
        }
      }
    }
#endif /* MHD */
  }
  return;
}


/*--------------------------------------------------------------------------*/
/*  Two functions to set BCs of Phi begin here                              */
/*--------------------------------------------------------------------------*/
/* slab_grav_bc_upper:
 *   Applies open bc in vertical direction to obtain values of potential in 
 *   ghost zones at: 
 *
 *      j>je (2D)
 *      k>ke (3D)
 *
 *   Assumes BCs have already been set in the in-plane direction (periodic)
 * 
 *   Switches for 1D, 2D, 3D
*/
static void slab_grav_bc_upper(GridS *pG)
{


#ifdef SELF_GRAVITY_USING_FFT_DISK

  int is = pG->is, ie = pG->ie;
  int js = pG->js, je = pG->je;
  int ks = pG->ks, ke = pG->ke;

  int i, j, k;
  int il,iu; /* i-lower/upper */
  int jl,ju; /* j-lower/upper */
  Real dx1sq=(pG->dx1*pG->dx1),dx2sq=(pG->dx2*pG->dx2);
  Real dx3sq=(pG->dx3*pG->dx3); 

  int dim;

  dim = 0;
  if(pG->Nx[0] > 1) dim++;
  if(pG->Nx[1] > 1) dim++;
  if(pG->Nx[2] > 1) dim++;


  iu = pG->ie + nghost;
  il = pG->is - nghost;

  switch(dim){
    case 1:
      /* ghost zones already calculated in selfg_by_fft_slab_1d ! */
      return;
      break;
    case 2:
/* 2D case -- apply open bc in x2 direction */
/* assumes periodic BC have already been applied in x1 direction to set 
   all ghost zones from is-nghost to ie+nghost
*/
   for (j=1; j<=nghost; j++){
     for (i=il+j; i<=iu-j; i++){
       pG->Phi[ks][je+j][i] = 2.0*pG->Phi[ks][je+(j-1)][i] - pG->Phi[ks][je+(j-2)][i] 
	 -(dx2sq/dx1sq)*(pG->Phi[ks][je+(j-1)][i+1] - 2.0*pG->Phi[ks][je+(j-1)][i]
                        +pG->Phi[ks][je+(j-1)][i-1])
         +dx2sq*four_pi_G*pG->U[ks][je+(j-1)][i].d;

      }
   }
/*
       for (j=je+nghost;j>=je;j--){
       for (i=1;i<=nghost;i++){
         pG->Phi[ks][j][is-i]=pG->Phi[ks][j][ie-(i-1)];
         pG->Phi[ks][j][ie+i]=pG->Phi[ks][j][is+(i-1)];
       }}
*/
    break;
    case 3:
/* 3D case -- apply open bc in x3 direction */
/* assumes periodic BC have already been applied in x1 and x2 directions 
   to set all ghost zones between is-nghost and ie+nghost 
  and js-nghost to je+nghost
*/
  ju = pG->je + nghost;
  jl = pG->js - nghost;
  for (k=1; k<=nghost; k++) {
    for (j=jl+k; j<=ju-k; j++){
      for (i=il+k; i<=iu-k; i++){
        pG->Phi[ke+k][j][i] =  2.0 * pG->Phi[ke+(k-1)][j][i] - pG->Phi[ke+(k-2)][j][i] 
          - (dx3sq/dx1sq)*(pG->Phi[ke+(k-1)][j][i+1]-2.0*pG->Phi[ke+(k-1)][j][i]+pG->Phi[ke+(k-1)][j][i-1])
          - (dx3sq/dx2sq)*(pG->Phi[ke+(k-1)][j+1][i]-2.0*pG->Phi[ke+(k-1)][j][i]+pG->Phi[ke+(k-1)][j-1][i])
          + dx3sq*four_pi_G*pG->U[ke+(k-1)][j][i].d;    
      }
    }
  }
   /* end of 3D case */
  }
#endif /* SELF_GRAVITY_USING_FFT_DISK */
  return;
}
/*--------------------------------------------------------------------------*/
/* slab_grav_bc_lower:
 *   Applies open bc in vertical direction to obtain values of potential in 
 *   ghost zones at: 
 *
 *      j<js (2D)
 *      k<ks (3D)
 *
 *   Assumes BCs have already been set in the in-plane direction (periodic)
 * 
 *   Switches for 1D, 2D, 3D
*/
static void slab_grav_bc_lower(GridS *pG)
{
#ifdef SELF_GRAVITY_USING_FFT_DISK

  int is = pG->is, ie = pG->ie;
  int js = pG->js, je = pG->je;
  int ks = pG->ks, ke = pG->ke;

  int i, j, k;
  int il,iu; /* i-lower/upper */
  int jl,ju; /* j-lower/upper */
  Real dx1sq=(pG->dx1*pG->dx1),dx2sq=(pG->dx2*pG->dx2);
  Real dx3sq=(pG->dx3*pG->dx3); 

  int dim;

  dim = 0;
  if(pG->Nx[0] > 1) dim++;
  if(pG->Nx[1] > 1) dim++;
  if(pG->Nx[2] > 1) dim++;


  iu = pG->ie + nghost;
  il = pG->is - nghost;

  switch(dim){
    case 1:
      /* ghost zones already calculated in selfg_by_fft_slab_1d ! */
      return;
      break;
    case 2:
/* 2D case -- apply open bc in x2 direction */
/* assumes periodic BC have already been applied in x1 direction to set 
   all ghost zones from is-nghost to ie+nghost
*/
    for (j=1; j<=nghost;j++){
      for (i=il+j; i<=iu-j; i++){  
        pG->Phi[ks][js-j][i] = 2.0*pG->Phi[ks][js-(j-1)][i] - pG->Phi[ks][js-(j-2)][i] 
          -(dx2sq/dx1sq)*(pG->Phi[ks][js-(j-1)][i+1] - 2.0*pG->Phi[ks][js-(j-1)][i]
                         +pG->Phi[ks][js-(j-1)][i-1])
          +dx2sq*four_pi_G*pG->U[ks][js-(j-1)][i].d;
      }
    }
/*
       for (j=js-nghost;j<=js;j++){
       for (i=1;i<=nghost;i++){
         pG->Phi[ks][j][is-i]=pG->Phi[ks][j][ie-(i-1)];
         pG->Phi[ks][j][ie+i]=pG->Phi[ks][j][is+(i-1)];
       }}
*/
    break;
    case 3:
/* 3D case -- apply open bc in x3 direction */
/* assumes periodic BC have already been applied in x1 and x2 directions 
   to set all ghost zones between is-nghost and ie+nghost and 
   js-nghost to je+nghost
*/
  ju = pG->je + nghost;
  jl = pG->js - nghost;
  for (k=1; k<=nghost; k++) {
    for (j=jl+k; j<=ju-k; j++){
      for (i=il+k; i<=iu-k; i++){
      pG->Phi[ks-k][j][i] =  2.0*pG->Phi[ks-(k-1)][j][i] - pG->Phi[ks-(k-2)][j][i] 
        -(dx3sq/dx1sq)*(pG->Phi[ks-(k-1)][j][i+1]-2.0*pG->Phi[ks-(k-1)][j][i]+pG->Phi[ks-(k-1)][j][i-1])
        -(dx3sq/dx2sq)*(pG->Phi[ks-(k-1)][j+1][i]-2.0*pG->Phi[ks-(k-1)][j][i]+pG->Phi[ks-(k-1)][j-1][i])
        +dx3sq*four_pi_G*pG->U[ks-(k-1)][j][i].d;
      }
    }
  }
   /* end of 3D case */
  }
#endif /* SELF_GRAVITY_USING_FFT_DISK */
  return;
}

/* Two functions to set BC of Phi end here */
/*--------------------------------------------------------------------------*/


