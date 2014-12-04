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
static void slab_grav_bc_lower(GridS *pG);
static void slab_grav_bc_upper(GridS *pG);
static void slab_mhd_bc_lower(GridS *pG);
static void slab_mhd_bc_upper(GridS *pG);
static Real UnstratifiedDisk(const Real x1, const Real x2, const Real x3);
static Real hst_sigma(const GridS *pG, const int i, const int j, const int k);
static Real hst_ux(const GridS *pG, const int i, const int j, const int k);
static Real hst_uy(const GridS *pG, const int i, const int j, const int k);
#ifdef MHD
static Real hst_m1(const GridS *pG, const int i, const int j, const int k);
static Real hst_m2(const GridS *pG, const int i, const int j, const int k);
#endif
#ifdef SELF_GRAVITY
static Real hst_dPhi(const GridS *pG, const int i, const int j, const int k);
#endif
#ifndef BAROTROPIC
static Real hst_dP(const GridS *pG, const int i, const int j, const int k);
#endif

static void initialize(DomainS *pD);
static int nwx,nwy;
static Real beta, amp, Q, nJ, cs, cs2, Gcons;
static Real kx, ky, Lx, Ly, Lz, Phi0;

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
  long int iseed;
  Real x1,x2,x3;
  Real rd, rp, rvx, rvy, rvz, rbx, rby, rbz;
  Real B0,P0;
  Real time0,kxt;
#ifdef SELF_GRAVITY
  Real Gcons;
#endif

  double rval;

  initialize(pDomain);

  time0=par_getd_def("problem","time0",0.0);

  B0 = cs/sqrt(beta);
#ifndef BAROTROPIC
  P0 = cs2/Gamma;
#endif

/* Ensure a different initial random seed for each process in an MPI calc. */
  ixs = pG->Disp[0];
  jxs = pG->Disp[1];
  kxs = pG->Disp[2];
  iseed = -1 - (ixs + pDomain->Nx[0]*(jxs + pDomain->Nx[1]*kxs));

  kxt = kx+qshear*Omega_0*ky*time0;

  pG->time=time0;

  for (k=ks; k<=ke; k++) {
  for (j=js; j<=je; j++) {
    for (i=is; i<=ie; i++) {
      cc_pos(pG,i,j,k,&x1,&x2,&x3);

      rd  = 1.0+amp*cos(kxt*x1+ky*x2);
      rvx = amp*kx/ky*sin(kxt*x1+ky*x2);
      rvy = amp*sin(kxt*x1+ky*x2);
      rvz = 0.0;
      rp  = cs2*(rd-1.0);

      rbx = amp*nwy*cos(kxt*(x1-0.5*pG->dx1)+ky*x2);
      rby = -amp*nwx*cos(kxt*x1+ky*(x2-0.5*pG->dx2));
      rbz = 0.0;

      pG->U[k][j][i].d  = rd;
      pG->U[k][j][i].M1 = rd*rvx;
      pG->U[k][j][i].M2 = rd*rvy;
#ifndef FARGO
      pG->U[k][j][i].M2 -= rd*(qshear*Omega_0*x1);
#endif
      pG->U[k][j][i].M3 = rd*rvz;
#ifdef ADIABATIC
      pG->U[k][j][i].E = (P0+rp)/Gamma_1
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
  Real my_Phi;
  int ierr;
#endif

  if(pM->time == 0.0){
  for (nl=0; nl<pM->NLevels; nl++){
    for (nd=0; nd<pM->DomainsPerLevel[nl]; nd++){
      if (pM->Domain[nl][nd].Grid != NULL) {
        pD = &(pM->Domain[nl][nd]);
        pG = pM->Domain[nl][nd].Grid;

        Phi0 = 0.0;
	dVol = pG->dx1*pG->dx2*pG->dx3;
        for (k=pG->ks; k<=pG->ke; k++) {
          for (j=pG->js; j<=pG->je; j++) {
            for (i=pG->is; i<=pG->ie; i++) {
              Phi0 += pG->Phi[k][j][i]*dVol;
            }
          }
        }

      }
    }
  }

#ifdef MPI_PARALLEL
  my_Phi = Phi0;
  ierr = MPI_Allreduce(&my_Phi, &Phi0, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
  Phi0 /= Lx*Ly*Lz;
  }



}

void Userwork_after_loop(MeshS *pM)
{
}

/*==============================================================================
 * PRIVATE FUNCTION:
 *============================================================================*/
void initialize(DomainS* pD){
  GridS *pG = pD->Grid;

  if(pG->Nx[2] == 1) ShBoxCoord = xy; /* 2D xy-plane */

/* Read problem parameters. */
  Omega_0 = par_getd("problem","omega");
  qshear = par_getd("problem","qshear");

  amp = par_getd("problem","amp");

/* Read parameters for magnetic field */
  beta = par_getd("problem","beta"); 

/* Read parameters for self gravity */
  Q=par_getd("problem","Q");
  nJ= par_getd("problem","nJ");

  cs=sqrt(4.0-2.0*qshear)/PI/nJ/Q;
  cs2=SQR(cs);

#ifdef SELF_GRAVITY
  Gcons = nJ*cs2;
  grav_mean_rho = 1.0;
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

  Lx = pD->MaxX[0]-pD->MinX[0];
  Ly = pD->MaxX[1]-pD->MinX[1];
  Lz = pD->MaxX[2]-pD->MinX[2];

/* enroll gravitational potential function */

  ShearingBoxPot = UnstratifiedDisk;

/* history dump for linear perturbation amplitude. See Kim & Ostriker 2001 */
  dump_history_enroll(hst_sigma, "<sigma>");
  dump_history_enroll(hst_ux, "<ux>");
  dump_history_enroll(hst_uy, "<uy>");
#ifdef MHD
  dump_history_enroll(hst_m1, "<m1>");
  dump_history_enroll(hst_m2, "<m2>");
#endif

#ifndef BAROTORPIC
  dump_history_enroll(hst_dP, "<dP>");
#endif

#ifdef SELF_GRAVITY
  dump_history_enroll(hst_dPhi, "<dPhi>");
#endif

}


static Real UnstratifiedDisk(const Real x1, const Real x2, const Real x3)
{
  Real phi=0.0;
#ifndef FARGO
  phi -= qshear*Omega_0*Omega_0*x1*x1;
#endif
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

static Real hst_ux(const GridS *pG, const int i, const int j, const int k)
{
  return 2*SQR(pG->U[k][j][i].M1/pG->U[k][j][i].d);
}

static Real hst_uy(const GridS *pG, const int i, const int j, const int k)
{
  Real vy0=0;
  Real x1,x2,x3;
  cc_pos(pG,i,j,k,&x1,&x2,&x3);
#ifndef FARGO
  vy0 -= qshear*Omega_0*x1;
#endif
  return 2*SQR(pG->U[k][j][i].M2/pG->U[k][j][i].d-vy0);
}

#ifdef MHD
static Real hst_m1(const GridS *pG, const int i, const int j, const int k)
{
  Real B0 = cs/sqrt(beta);
  return 2*SQR(pG->U[k][j][i].B1c/B0);
}

static Real hst_m2(const GridS *pG, const int i, const int j, const int k)
{
  Real B0 = cs/sqrt(beta);

  return 2*SQR(1-pG->U[k][j][i].B2c/B0);
}
#endif

#ifndef BAROTROPIC
static Real hst_dP(const GridS *pG, const int i, const int j, const int k)
{
  Real Eint,Ekin,Emag=0.0;

  Ekin = 0.5*(SQR(pG->U[k][j][i].M1)+SQR(pG->U[k][j][i].M2)
             +SQR(pG->U[k][j][i].M3))/pG->U[k][j][i].d;
#ifdef MHD
  Emag = 0.5*(SQR(pG->U[k][j][i].B1c)+SQR(pG->U[k][j][i].B2c)+SQR(pG->U[k][j][i].B3c));
#endif
  Eint = pG->U[k][j][i].E - Ekin - Emag;

  return 2*SQR(Eint/Gamma_1-cs2/Gamma);
}
#endif

#ifdef SELF_GRAVITY
static Real hst_dPhi(const GridS *pG, const int i, const int j, const int k)
{
  Real dPhi;
  Real kxt,k2,k0;
  Real x1,x2,x3;
  Real gz=1.0,gz2=1.0,hLz,kLz;
  int ks = pG->ks, ke = pG->ke;

  cc_pos(pG,i,j,k,&x1,&x2,&x3);

  kxt = kx+qshear*ky*pG->time;
  k2 =kxt*kxt+ky*ky;
  k0 =sqrt(kx*kx+ky*ky);

#ifdef SELF_GRAVITY_USING_FFT_DISK
  hLz=0.5*Lz;
  kLz=sqrt(k2)*Lz;
  gz =1-0.5*exp(-k0*hLz)*(exp(k0*fabs(x3))+exp(-k0*fabs(x3)));
  gz2 = 1+0.5*exp(-kLz)-2.0*(1-exp(-kLz))/kLz+0.25*(1-exp(-2*kLz))/kLz;
#endif


  dPhi = pG->Phi[k][j][i]-Phi0;//*k2/(4.0*PI*nJ*cs2);
/*
#ifdef SELF_GRAVITY_USING_FFT_DISK
  dPhi /=  1-0.5*(exp(-sqrt(k2)*(zmax-fabs(x3)))+exp(-sqrt(k2)*(zmax+fabs(x3))));
#endif
*/

  return 2*SQR(dPhi);
}
#endif
