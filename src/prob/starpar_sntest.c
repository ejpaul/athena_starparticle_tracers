#include "copyright.h"
/*==============================================================================
 * FILE: starpar_test.c
 *
 * PURPOSE: Problem generator for thermal instability with operator split cooling.
 *          iprob=1 for eigenmode test
 *          iprob=2 for white noise
 *          iprob=3 for gaussian random field
 *============================================================================*/

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"

#ifdef STAR_PARTICLE
#define TEST_STAR
#else
//#define EJECTA
//#define SEDOV
#endif

//#define THERMAL_EQUILIBRIUM

#ifdef MHD
#error This problem does not work with MHD
#endif

/* FFT indexing Nfast=k, Nmid=j, Nslow=i (opposite to Athena)
 *  * For OFST, i,j,k,nx2,nx3 reference the local grid */
#define OFST(i, j, k) ((k) + nx3*((j) + nx2*(i)))
/* KWVM: magnitude of wavenumber k in units of dkx */
#define KWVM(i, j, k) (sqrt(SQR(KCOMP(i,gis,gnx1))+ \
                            SQR(KCOMP(j,gjs,gnx2))+SQR(KCOMP(k,gks,gnx3))))


/* variables for generate gaussan random field */
/* Number of cells in local grid, number of cells in global grid */
static int nx1,nx2,nx3,gnx1,gnx2,gnx3;
/* Starting indices for global grid */
static int gis,gjs,gks;
/* Cutoff wavenumbers, G&O spect peak, power law spect exponent, 2 pi/L */
static Real klow,khigh,kpeak,expo,dkx,Lx;
/* Driving properties */
static int ispect;
/* Seed for random number generator */
static long int rseed;

/* FFTW - Variables, Plan, etc. */
/* These are made static global variables so that they need not be
 *    allocated AND destroyed with each call to pspect! */
static struct ath_3d_fft_plan *plan;
/* Between calls to generate(), these have unshifted, unnormalized
 *  * velocity perturbations. */
static ath_fft_data *frho=NULL;
static Real ***drho=NULL, std_dev;
static Real ***dv1=NULL, ***dv2=NULL, ***dv3=NULL;

static ConstS consts;
static Real Punit, dmax, rin, ton, toff;
#ifndef STAR_PARTICLE
static Real tSN,fm,ESN;
#endif

/*==============================================================================
 * History function:
 *============================================================================*/
static Real hst_dmax(const GridS *pG, const int i, const int j, const int k);
static Real hst_radial_mom(const GridS *pG, const int i, const int j, const int k);
static Real hst_Mcold(const GridS *pG, const int i, const int j, const int k);
static Real hst_Minter(const GridS *pG, const int i, const int j, const int k);
static Real hst_Mwarm(const GridS *pG, const int i, const int j, const int k);
static Real hst_Mhot(const GridS *pG, const int i, const int j, const int k);
static Real hst_Rshell(const GridS *pG, const int i, const int j, const int k);
static Real hst_Mshell(const GridS *pG, const int i, const int j, const int k);
static Real hst_Mshell2(const GridS *pG, const int i, const int j, const int k);
static Real hst_Vcenter(const GridS *pG, const int i, const int j, const int k);
static Real hst_Tcenter(const GridS *pG, const int i, const int j, const int k);
static Real hst_ncenter(const GridS *pG, const int i, const int j, const int k);
static Real hst_Pcenter(const GridS *pG, const int i, const int j, const int k);
static Real hst_Vhot(const GridS *pG, const int i, const int j, const int k);
static Real hst_Thot(const GridS *pG, const int i, const int j, const int k);
static Real hst_nhot(const GridS *pG, const int i, const int j, const int k);
static Real hst_Phot(const GridS *pG, const int i, const int j, const int k);
static Real hst_Pshell(const GridS *pG, const int i, const int j, const int k);
static Real hst_rmom_hot(const GridS *pG, const int i, const int j, const int k);
static Real hst_rmom_shell(const GridS *pG, const int i, const int j, const int k);

#ifdef STAR_PARTICLE
static Real hst_age_sp(const GridS *pG, const int i, const int j, const int k);
static Real hst_mass_in_sp(const GridS *pG, const int i, const int j, const int k);
static Real hst_mass_in_sp_ghost(const GridS *pG, const int i, const int j, const int k);
static Real hst_mass_in_sp_ghost2(const GridS *pG, const int i, const int j, const int k);
#endif
/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 * ran2()
 *============================================================================*/

Real Teq(const Real Pok, const Real temp);
Real neq(const Real temp);
#ifdef SHEARING_BOX
static Real ShearPot(const Real x1, const Real x2, const Real x3);
#endif
static void initialize_PS(DomainS *pDomain);
static void pspect(ath_fft_data *ampl);

static void initialize(DomainS *pD);
static double ran2(long int *idum);
static Real logd(const GridS *pG, const int i, const int j, const int k);

void correct_bad_zones(DomainS* pD);

#ifdef TEST_STAR
static void create_teststar(GridS *pG, Real m0, Real x0, Real y0, Real z0, Real v0, int starpar_id);
static Real GravPot(const Real x1, const Real x2, const Real x3);
static Real GM;
static int tag=0;
#endif
#ifndef STAR_PARTICLE
void get_Sedov_Taylor(Real xi, Real *alpha, Real *v, Real *p);
void assign_SN(GridS *pGrid);
#endif
/*=========================== PUBLIC FUNCTIONS ===============================*/
/*----------------------------------------------------------------------------*/
/* problem:  */

void problem(DomainS *pDomain)
{
  GridS *pGrid=(pDomain->Grid);
  int i=0,j=0,k=0;
  int is,ie,js,je,ks,ke;
  int iprob,converged;
  Real n0,P0,v0,r,rout,n1,Tc,Tw;
  long int iseed = -1-myID_Comm_world;

#ifdef MPI_PARALLEL
  Real my_max;
  int ierr;
#endif

  Real x1,x2,x3;
  Real x1min,x1max,kwx,amp;

  is = pGrid->is; ie = pGrid->ie;
  js = pGrid->js; je = pGrid->je;
  ks = pGrid->ks; ke = pGrid->ke;

/* Get box size and set wavenumber */
  x1min = pDomain->MinX[0];
  x1max = pDomain->MaxX[0];
  Lx = x1max - x1min;
  kwx = 2.0*PI/Lx;
  amp   = par_getd("problem","amp");
  rout = par_getd("problem","rout");

  initialize(pDomain);
/* Read problem parameters */

  P0    = par_getd("problem","P0")*pGrid->heat_ratio;
#ifdef THERMAL_EQUILIBRIUM
  converged = bisection(Teq,1.,184.,P0,&Tc);
  converged = bisection(Teq,5000.,10000.,P0,&Tw);
  n0 = P0/Tc/1.1;
  n1 = P0/Tw/1.1;
#else
  n0    = par_getd("problem","n0");
  Tc	= P0/1.1/n0;
  n1    = par_getd("problem","n1");
#endif
  P0    *= Punit;

  iprob = par_geti_def("problem","iprob",2);
  if(myID_Comm_world == 0) printf("P0/k0 = %g, P0 in code unit = %g\n",1.1*n0*Tc,P0);
  if(iprob == 3 || iprob == 4){
    rseed = -(unsigned)time(NULL);
    if(myID_Comm_world==0) printf("seed using time = %ld\n",rseed);
    rseed -= (gis + pDomain->Nx[0]*(gjs + pDomain->Nx[1]*gks));
#ifndef FFT_ENABLED
    ath_error("[problem]: --enable-fft is needed to generate PS\n");
#endif
  }

  if(iprob == 3){
    if ((drho=(Real***)calloc_3d_array(nx3+2*nghost,nx2+2*nghost,nx1+2*nghost,sizeof(Real)))==NULL)      ath_error("[problem]: Error allocating memory for vel pert\n");
    initialize_drho(pDomain,drho,&std_dev,rseed);
  }

  if(iprob == 4){
    if ((dv1=(Real***)calloc_3d_array(nx3+2*nghost,nx2+2*nghost,nx1+2*nghost,sizeof(Real)))==NULL)      ath_error("[problem]: Error allocating memory for vel pert\n");
    if ((dv2=(Real***)calloc_3d_array(nx3+2*nghost,nx2+2*nghost,nx1+2*nghost,sizeof(Real)))==NULL)      ath_error("[problem]: Error allocating memory for vel pert\n");
    if ((dv3=(Real***)calloc_3d_array(nx3+2*nghost,nx2+2*nghost,nx1+2*nghost,sizeof(Real)))==NULL)      ath_error("[problem]: Error allocating memory for vel pert\n");
    initialize_dv(pDomain,dv1,dv2,dv3,rseed,1); // 1 for projection to Div V =0
  }

/* Constant density and temperature initially */
  dmax=1.e-30;

  for (k=ks; k<=ke; k++) {
  for (j=js; j<=je; j++) {
    for (i=is; i<=ie; i++) {
      cc_pos(pGrid,i,j,k,&x1,&x2,&x3);
      r=sqrt(SQR(x1)+SQR(x2)+SQR(x3));

      if(r<=rin) pGrid->U[k][j][i].d = n0;
      else pGrid->U[k][j][i].d = n0*SQR(rin/r);

      if(r>rout) pGrid->U[k][j][i].d = n1;

      if(iprob == 1) pGrid->U[k][j][i].d *= (1+amp*sin(kwx*x1));
      if(iprob == 2) pGrid->U[k][j][i].d *= (1+amp*(ran2(&iseed)-0.5));
      if(iprob == 3) pGrid->U[k][j][i].d *= (1+amp*(drho[k][j][i]/std_dev));
      pGrid->U[k][j][i].M1 = 0.0;
      pGrid->U[k][j][i].M2 = 0.0;
      pGrid->U[k][j][i].M3 = 0.0;
#ifdef SHEARING_BOX
#ifndef FARGO
      pGrid->U[k][j][i].M2 -= qshear*Omega_0*x1*pGrid->U[k][j][i].d;
#endif
#endif
#ifndef BAROTROPIC
      pGrid->U[k][j][i].E = P0/Gamma_1
             + 0.5*(SQR(pGrid->U[k][j][i].M1) + SQR(pGrid->U[k][j][i].M2)
             + SQR(pGrid->U[k][j][i].M3))/pGrid->U[k][j][i].d;

#endif
      dmax = MAX(dmax,pGrid->U[k][j][i].d);
    }
  }}

#ifdef MPI_PARALLEL
  my_max = dmax;
  ierr = MPI_Allreduce(&my_max, &dmax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#endif

#ifndef STAR_PARTICLE
  if(ton == 0.0){
    assign_SN(pGrid);
    ton += tSN;
  }
#endif

#ifdef TEST_STAR
  if(myID_Comm_world==0 && ton==0) {
    create_teststar(pGrid,1.e3/pGrid->units.Msun,0.5*pGrid->dx1,0.5*pGrid->dx2,0.5*pGrid->dx3,0.0,0);
    tag=1;
    printf("creating star particles at t=%g\n",pGrid->time);
  }
#endif

  if(iprob == 3) free_3d_array(drho);
  if(iprob == 4){
    shift_normalize(pGrid,dv1,dv2,dv3,amp);
    free_3d_array(dv1);
    free_3d_array(dv2);
    free_3d_array(dv3);
  }

}

/*==============================================================================
 * PROBLEM USER FUNCTIONS:
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

void problem_read_restart(MeshS *pM, FILE *fp)
{

  DomainS *pD=NULL;
  GridS *pG=NULL;
  int i,j,k;
  int nl,nd;
  Real eint;

#ifdef MPI_PARALLEL
  Real my_max;
  int ierr;
#endif

  pM->time = par_getd("time","time");
  pM->nstep = par_getd("time","nstep");
  for (nl=0; nl<pM->NLevels; nl++){
    for (nd=0; nd<pM->DomainsPerLevel[nl]; nd++){
      if (pM->Domain[nl][nd].Grid != NULL) {
        pD = &(pM->Domain[nl][nd]);
        pG = pM->Domain[nl][nd].Grid;

        initialize(pD);

        dmax = 1.e-30;
        for (k=pG->ks; k<=pG->ke; k++) {
          for (j=pG->js; j<=pG->je; j++) {
            for (i=pG->is; i<=pG->ie; i++) {
              dmax = MAX(dmax,pG->U[k][j][i].d);
#ifndef BAROTROPIC
              eint = pG->U[k][j][i].E - 0.5/pG->U[k][j][i].d*(SQR(pG->U[k][j][i].M1)
                    + SQR(pG->U[k][j][i].M2)+ SQR(pG->U[k][j][i].M3));
#endif
              pG->U[k][j][i].M1=0.0;
              pG->U[k][j][i].M2=0.0;
              pG->U[k][j][i].M3=0.0;
#ifndef BAROTROPIC
              pG->U[k][j][i].E = eint;
#endif
            }
          }
        }

#ifndef STAR_PARTICLE
        pG->time = pM->time;
        if(pG->time >= ton && ton >= 0.0){
          assign_SN(pG);
          ton += tSN;
        }
#endif
      }
    }
  }

  new_dt(pM);
//  data_output(pM,1);

#ifdef MPI_PARALLEL
  my_max = dmax;
  ierr = MPI_Allreduce(&my_max, &dmax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#endif



  return;
}

static Real logd(const GridS *pG, const int i, const int j, const int k)
{
  return log10(pG->U[k][j][i].d);
}


ConsFun_t get_usr_expr(const char *expr)
{ 
  if(strcmp(expr,"logd")==0) return logd;
  return NULL;
}

VOutFun_t get_usr_out_fun(const char *name){
  return NULL;
}

void Userwork_in_loop(MeshS *pM)
{

  DomainS *pD=NULL;
  GridS *pG=NULL;
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

#ifdef TEST_STAR
        if(myID_Comm_world==0 && pG->time > ton && tag == 0) {
          create_teststar(pG,1.e3/pG->units.Msun,0.5*pG->dx1,0.5*pG->dx2,0.5*pG->dx3,0.0,0);
          printf("creating star particles at t=%g\n",pM->time);
          tag=1;
        }
#endif
#ifndef STAR_PARTICLE
        if(pG->time >= toff) ton = -1;

        if(pG->time >= ton && ton >= 0.0){
          assign_SN(pG);
          //data_output(pM, 1);
          ton += tSN; 
        }
#endif
	correct_bad_zones(pD);
      }
    }
  }



  return;
}

void Userwork_after_loop(MeshS *pM)
{
  return;
}

/*=========================== PRIVATE FUNCTIONS ==============================*/

Real Teq(const Real Pok, const Real temp){
  Real nden;
 
  nden=Pok/1.1/temp;

  return (nden/neq(temp)-1.0);
}

Real neq(const Real temp){
  Real c1=1.e7, c2=1.4e-2, t1=1.184e5, t2=92.;
 
  return 1.0/(c2*sqrt(temp)*exp(-t2/temp)+c1*exp(-t1/(temp+1.e3)));
}


void correct_bad_zones(DomainS* pD){
  GridS *pG = (pD->Grid);
  int is = pG->is, ie = pG->ie;
  int js = pG->js, je = pG->je;
  int ks = pG->ks, ke = pG->ke;
  int i,j,k;
  int ii,jj,kk;
  int ibad,Ncells;

  Real nden,cs2,P,norignal;
  Real eint,v1,v2,v3;

  ConsS U;

  for(k=ks;k<=ke;k++){ 
    for(j=js;j<=je;j++){ 
      for(i=is;i<=ie;i++){ 
        ibad = 0;
        U=pG->U[k][j][i];

        eint = (U.E - (0.5/U.d)*(SQR(U.M1)+SQR(U.M2)+SQR(U.M3)));
#ifdef MHD
        eint -= (0.5)*(SQR(U.B1c)+SQR(U.B2c)+SQR(U.B3c));
#endif
        nden = U.d;
        cs2 = Gamma_1*eint/U.d;

        if(nden < 0 || nden != nden) ibad=1;
        if(cs2 < 0.06) ibad=2;

	if(ibad != 0){
	  if(ibad == 1){
            ath_perr(-1,"[Userwork in Loop] Bad cell is detected at id=%d [%d][%d][%d]\n",
                    k,j,i,myID_Comm_world);
            ath_perr(-1,"[Userwork in Loop] Original: n = %g, eint= %g\n",nden,eint);
          }

          eint = 0.0;
          v1 = 0.0;
          v2 = 0.0;
          v3 = 0.0;
          nden = 0.0;
          Ncells=0;

/* Using six zones
// for [k][j][i-1] and [k][j][i+1] 
          jj = j;
          kk = k;
          for(ii=i-1;ii<=i+1;ii+=2){ 
            U=pG->U[kk][jj][ii];
            P = (U.E - (0.5/U.d)*(SQR(U.M1)+SQR(U.M2)+SQR(U.M3)));
#ifdef MHD
            P -= (0.5)*(SQR(U.B1c)+SQR(U.B2c)+SQR(U.B3c));
#endif
            P *= Gamma_1;

            if((U.d > 0) && (P/U.d > 0.06)){
              eint += P/Gamma_1;
              v1 += U.M1;
              v2 += U.M2;
              v3 += U.M3;
              nden += U.d; 
              Ncells++;
            }
          }

// for [k][j-1][i] and [k][j+1][i] 
          ii = i;
          kk = k;
          for(jj=j-1;jj<=j+1;jj+=2){ 
            U=pG->U[kk][jj][ii];
            P = (U.E - (0.5/U.d)*(SQR(U.M1)+SQR(U.M2)+SQR(U.M3)));
#ifdef MHD
            P -= (0.5)*(SQR(U.B1c)+SQR(U.B2c)+SQR(U.B3c));
#endif
            P *= Gamma_1;

            if((U.d > 0) && (P/U.d > 0.06)){
              eint += P/Gamma_1;
              v1 += U.M1;
              v2 += U.M2;
              v3 += U.M3;
              nden += U.d; 
              Ncells++;
            }
          }

// for [k-1][j][i] and [k+1][j][i] 
          ii = i;
          jj = j;
          for(kk=k-1;kk<=k+1;kk+=2){ 
            U=pG->U[kk][jj][ii];
            P = (U.E - (0.5/U.d)*(SQR(U.M1)+SQR(U.M2)+SQR(U.M3)));
#ifdef MHD
            P -= (0.5)*(SQR(U.B1c)+SQR(U.B2c)+SQR(U.B3c));
#endif
            P *= Gamma_1;

            if((U.d > 0) && (P/U.d > 0.06)){
              eint += P/Gamma_1;
              v1 += U.M1;
              v2 += U.M2;
              v3 += U.M3;
              nden += U.d; 
              Ncells++;
            }
          }
*/
/* Searching for 27 zones */

          for(kk=k-1;kk<=k+1;kk++){ 
          for(jj=j-1;jj<=j+1;jj++){ 
          for(ii=i-1;ii<=i+1;ii++){ 
            U=pG->U[kk][jj][ii];
            P = (U.E - (0.5/U.d)*(SQR(U.M1)+SQR(U.M2)+SQR(U.M3)));
#ifdef MHD
            P -= (0.5)*(SQR(U.B1c)+SQR(U.B2c)+SQR(U.B3c));
#endif
            P *= Gamma_1;

            if((U.d > 0) && (P/U.d > 0.06)){
              eint += P/Gamma_1;
              v1 += U.M1;
              v2 += U.M2;
              v3 += U.M3;
              nden += U.d; 
              Ncells++;
            }
          }}}


          eint = eint/(Real)Ncells;
          v1 = v1/nden;
          v2 = v2/nden;
          v3 = v3/nden;
          nden = nden/(Real)Ncells;
	  if(ibad == 1){
            ath_perr(-1,"[Userwork in Loop] Averaged: n = %g, eint= %g, Ncells=%d\n",
                     nden,eint,Ncells);
            pG->U[k][j][i].d = nden;
            pG->U[k][j][i].M1 = nden*v1;
            pG->U[k][j][i].M2 = nden*v2;
            pG->U[k][j][i].M3 = nden*v3;
          }
          pG->U[k][j][i].E = eint + 0.5/pG->U[k][j][i].d*(SQR(pG->U[k][j][i].M1)+SQR(pG->U[k][j][i].M2)+SQR(pG->U[k][j][i].M3));

        }
      }
    }
  }

  return;

}
static void initialize(DomainS* pD){
  GridS *pG = (pD->Grid);

  Real kappa;

  /* Get local grid size */
  nx1 = pG->Nx[0];
  nx2 = pG->Nx[1];
  nx3 = pG->Nx[2];

  /* Get global grid size */
  gnx1 = pD->Nx[0];
  gnx2 = pD->Nx[1];
  gnx3 = pD->Nx[2];

  /* Get extents of local FFT grid in global coordinates */
  gis=pG->Disp[0];
  gjs=pG->Disp[1];
  gks=pG->Disp[2];

  init_consts(&consts);
  pG->units.Dcode = 1.4*consts.mH;
  pG->units.Vcode = consts.kms;
  pG->units.Lcode = consts.pc;
  pG->units.Mcode = pG->units.Dcode*CUBE(pG->units.Lcode);
  pG->units.Tcode = pG->units.Lcode/pG->units.Vcode; 
  Punit = consts.kB/pG->units.Dcode/SQR(pG->units.Vcode);
  pG->units.Msun = pG->units.Mcode/consts.Msun;
  pG->units.Myr = pG->units.Tcode/consts.Myr;
  if(myID_Comm_world==0) printf("Myr in code = %g, Msun in code= %g\n",pG->units.Myr,pG->units.Msun);

#ifdef SELF_GRAVITY
  pG->units.G = consts.G*pG->units.Dcode*SQR(pG->units.Lcode)/SQR(pG->units.Vcode);
#endif

#ifdef OPERATOR_SPLIT_COOLING
  pG->heat0 = 2.e-26;
  pG->heat_ratio = par_getd_def("problem","heat_ratio",1.0);
#endif

#ifdef THERMAL_CONDUCTION
  kappa = par_getd("problem","kappa");
  kappa_iso = kappa/(consts.kB*pG->units.Lcode*pG->units.Vcode);
  if(myID_Comm_world==0) printf("kappa in c.g.s. = %g, in code unit = %g\n",kappa,kappa_iso);
#endif

#ifdef SELF_GRAVITY
  four_pi_G = 4.0*PI*pG->units.G;
  grav_mean_rho = 0.0;
#endif

#ifdef SHEARING_BOX
  Omega_0 = par_getd_def("problem","Omega",28.e-3);
  qshear  = par_getd_def("problem","qshear",1.0);

  ShearingBoxPot = ShearPot;
#endif

#ifdef STAR_PARTICLE
  tHII = par_getd_def("problem","tHII",3.0);
#endif
  tSN = par_getd_def("problem","tSN",0.0)/pG->units.Myr; // in units of Myr
  ton = par_getd("problem","ton")/pG->units.Myr; // in units of Myr
  toff = par_getd_def("problem","toff",1.0)/pG->units.Myr; // in units of Myr
  rin = par_getd("problem","rin");
  fm = CUBE(par_getd_def("problem","fm",0.3));
  ESN = par_getd_def("problem","ESN",1.0);

  dump_history_enroll(hst_dmax, "dmax");
  dump_history_enroll(hst_radial_mom, "r_mom");
  dump_history_enroll(hst_Mcold, "Mcold");
  dump_history_enroll(hst_Minter, "Minter");
  dump_history_enroll(hst_Mwarm, "Mwarm");
  dump_history_enroll(hst_Mhot, "Mhot");
  dump_history_enroll(hst_Mshell, "Msh");
  dump_history_enroll(hst_Mshell2, "Msh2");
  dump_history_enroll(hst_Rshell, "Rsh");
  dump_history_enroll(hst_ncenter, "nc");
  dump_history_enroll(hst_Pcenter, "Pc");
  dump_history_enroll(hst_Tcenter, "Tc");
  dump_history_enroll(hst_Vcenter, "Vc");
  dump_history_enroll(hst_nhot, "nhot");
  dump_history_enroll(hst_Phot, "Phot");
  dump_history_enroll(hst_Thot, "Thot");
  dump_history_enroll(hst_Vhot, "Vhot");
  dump_history_enroll(hst_Pshell, "Pshell");
  dump_history_enroll(hst_rmom_hot, "rmom_hot");
  dump_history_enroll(hst_rmom_shell, "rmom_shell");
#ifdef STAR_PARTICLE
  dump_history_enroll(hst_age_sp, "age");
  dump_history_enroll(hst_mass_in_sp, "msp");
  dump_history_enroll(hst_mass_in_sp_ghost, "mghost");
  dump_history_enroll(hst_mass_in_sp_ghost, "mghost2");
#endif
}

#ifdef TEST_STAR
static void create_teststar(GridS *pG, Real m0, Real x0, Real y0, Real z0, Real v0, int starpar_id)
{
  StarParListS *pList=NULL;
  StarParS *pStar=NULL;

  Real xp,yp,zp;

  if ((pList = (StarParListS *) calloc(1,sizeof(StarParListS))) == NULL)
    ath_error("[push_to_local_list]: Error callocing memory for starparlist\n");

  pStar = &(pList->starpar);
  pStar->id = starpar_id;
  pStar->merge_history = 0;
  pStar->isnew = 1;


  pStar->m = m0;
  pStar->x1 = x0;
  pStar->v1 = v0;
#ifdef SHEARING_BOX
  pStar->x2 = v0/Omega_0/(2-qshear)+y0;
  pStar->v2 = -qshear*Omega_0*x0;
#else
  pStar->x2 = y0;
  pStar->v2 = 0.0;
#endif
  pStar->x3 = z0;
  pStar->v3 = 0.0;
  pStar->age = 0.0;
  pStar->mdot = 0.0;
  pStar->mghost = 0.0;
  pStar->M1ghost = 0.0;
  pStar->M2ghost = 0.0;
  pStar->M3ghost = 0.0;
  pStar->navg = 0.0;
  pStar->v1avg = 0.0;
  pStar->v2avg = 0.0;
  pStar->v3avg = 0.0;
  pStar->eavg = 0.0;

  starpar_push_local(pG,pList);
  
  starpar_printlist(-1, pG);

  return;
}
#endif

static Real hst_dmax(const GridS *pG, const int i, const int j, const int k)
{
  return dmax;
}

static Real hst_radial_mom(const GridS *pG, const int i, const int j, const int k)
{
  Real x1,x2,x3;
  Real x10,x20,x30;
  Real rmom,r;
#ifdef STAR_PARTICLE
  StarParListS *pList=NULL;
  StarParS *pStar=NULL;

  pList = pG->Gstars;
  if(pList) {
    pStar = &(pList->starpar);
    x10=pStar->x1;
    x20=pStar->x2;
    x30=pStar->x3;
  }
#else
  x10=0.0;
  x20=0.0;
  x30=0.0;
#endif


  cc_pos(pG,i,j,k,&x1,&x2,&x3);
  r=sqrt(SQR(x1)+SQR(x2)+SQR(x3));
  rmom = pG->U[k][j][i].M1*(x1-x10)/r
        +pG->U[k][j][i].M2*(x2-x20)/r
        +pG->U[k][j][i].M3*(x3-x30)/r;
  return rmom;
}

static Real hst_Mcold(const GridS *pG, const int i, const int j, const int k)
{
  Real Press,Temp;

  Press = pG->U[k][j][i].E - 0.5/pG->U[k][j][i].d*(SQR(pG->U[k][j][i].M1)
                           + SQR(pG->U[k][j][i].M2)+ SQR(pG->U[k][j][i].M3));
  Press = Press*Gamma_1;
  Press = Press*pG->units.Dcode*SQR(pG->units.Vcode)/consts.kB;
  Temp = Press/1.1/pG->U[k][j][i].d;
  
  if(Temp <= 184.) return pG->U[k][j][i].d;
  else return 0.;
}


static Real hst_Minter(const GridS *pG, const int i, const int j, const int k)
{
  Real Press,Temp;

  Press = pG->U[k][j][i].E - 0.5/pG->U[k][j][i].d*(SQR(pG->U[k][j][i].M1)
                           + SQR(pG->U[k][j][i].M2)+ SQR(pG->U[k][j][i].M3));
  Press = Press*Gamma_1;
  Press = Press*pG->units.Dcode*SQR(pG->units.Vcode)/consts.kB;
  Temp = Press/1.1/pG->U[k][j][i].d;
  
  if(Temp > 184. && Temp < 5050) return pG->U[k][j][i].d;
  else return 0.;
}


static Real hst_Mwarm(const GridS *pG, const int i, const int j, const int k)
{
  Real Press,Temp;

  Press = pG->U[k][j][i].E - 0.5/pG->U[k][j][i].d*(SQR(pG->U[k][j][i].M1)
                           + SQR(pG->U[k][j][i].M2)+ SQR(pG->U[k][j][i].M3));
  Press = Press*Gamma_1;
  Press = Press*pG->units.Dcode*SQR(pG->units.Vcode)/consts.kB;
  Temp = Press/1.1/pG->U[k][j][i].d;
  
  if(Temp >= 5050. && Temp <2.e4 ) return pG->U[k][j][i].d;
  else return 0.;
}


static Real hst_Mhot(const GridS *pG, const int i, const int j, const int k)
{
  Real Press,Temp;

  Press = pG->U[k][j][i].E - 0.5/pG->U[k][j][i].d*(SQR(pG->U[k][j][i].M1)
                           + SQR(pG->U[k][j][i].M2)+ SQR(pG->U[k][j][i].M3));
  Press = Press*Gamma_1;
  Press = Press*pG->units.Dcode*SQR(pG->units.Vcode)/consts.kB;
  Temp = Press/1.1/pG->U[k][j][i].d;
  
  if(Temp >= 2.e4) return pG->U[k][j][i].d;
  else return 0.;
}

static Real hst_Rshell(const GridS *pG, const int i, const int j, const int k)
{
  Real r,x1,x2,x3;
  Real Press,Temp;
  Real vr=0.0;

  cc_pos(pG,i,j,k,&x1,&x2,&x3);
  r=sqrt(SQR(x1)+SQR(x2)+SQR(x3));
  if(r != 0.0) vr = (pG->U[k][j][i].M1*x1+pG->U[k][j][i].M2*x2+pG->U[k][j][i].M3*x3)/(pG->U[k][j][i].d*r);

  Press = pG->U[k][j][i].E - 0.5/pG->U[k][j][i].d*(SQR(pG->U[k][j][i].M1)
                           + SQR(pG->U[k][j][i].M2)+ SQR(pG->U[k][j][i].M3));
  Press = Press*Gamma_1;
  Press = Press*pG->units.Dcode*SQR(pG->units.Vcode)/consts.kB;
  Temp = Press/1.1/pG->U[k][j][i].d;
  
  if(Temp < 2.e4 && vr > 1.0) return r*pG->U[k][j][i].d;
  else return 0.;
}


static Real hst_Mshell(const GridS *pG, const int i, const int j, const int k)
{
  Real r,x1,x2,x3;
  Real Press,Temp;
  Real vr=0.0;

  cc_pos(pG,i,j,k,&x1,&x2,&x3);
  r=sqrt(SQR(x1)+SQR(x2)+SQR(x3));
  if(r != 0.0) vr = (pG->U[k][j][i].M1*x1+pG->U[k][j][i].M2*x2+pG->U[k][j][i].M3*x3)/(pG->U[k][j][i].d*r);

  Press = pG->U[k][j][i].E - 0.5/pG->U[k][j][i].d*(SQR(pG->U[k][j][i].M1)
                           + SQR(pG->U[k][j][i].M2)+ SQR(pG->U[k][j][i].M3));
  Press = Press*Gamma_1;
  Press = Press*pG->units.Dcode*SQR(pG->units.Vcode)/consts.kB;
  Temp = Press/1.1/pG->U[k][j][i].d;
  
  if(Temp < 2.e4 && vr > 1.0) return pG->U[k][j][i].d;
  else return 0.;
}

static Real hst_Mshell2(const GridS *pG, const int i, const int j, const int k)
{
  Real r,x1,x2,x3;
  Real Press,Temp;

  cc_pos(pG,i,j,k,&x1,&x2,&x3);

  Press = pG->U[k][j][i].E - 0.5/pG->U[k][j][i].d*(SQR(pG->U[k][j][i].M1)
                           + SQR(pG->U[k][j][i].M2)+ SQR(pG->U[k][j][i].M3));
  Press = Press*Gamma_1;
  Press = Press*pG->units.Dcode*SQR(pG->units.Vcode)/consts.kB;
  Temp = Press/1.1/pG->U[k][j][i].d;
  
  if((Temp < 2.e4) && (Press > 5000.*pG->heat_ratio)) return pG->U[k][j][i].d;
  else return 0.;
}

static Real hst_Tcenter(const GridS *pG, const int i, const int j, const int k)
{
  Real r,x1,x2,x3;
  Real Press,Temp;

  cc_pos(pG,i,j,k,&x1,&x2,&x3);
  r=sqrt(SQR(x1)+SQR(x2)+SQR(x3));
  
  if(r <= rin){
    Press = pG->U[k][j][i].E - 0.5/pG->U[k][j][i].d*(SQR(pG->U[k][j][i].M1)
                             + SQR(pG->U[k][j][i].M2)+ SQR(pG->U[k][j][i].M3));
    Press = Press*Gamma_1;
    Press = Press*pG->units.Dcode*SQR(pG->units.Vcode)/consts.kB;
    Temp = Press/1.1/pG->U[k][j][i].d;

    return Temp;
  } else return 0.;
}

static Real hst_Pcenter(const GridS *pG, const int i, const int j, const int k)
{
  Real r,x1,x2,x3;
  Real Press,Temp;
  cc_pos(pG,i,j,k,&x1,&x2,&x3);

  r=sqrt(SQR(x1)+SQR(x2)+SQR(x3));
  
  if(r <= rin){
//  if((fabs(x1)<pG->dx1) && (fabs(x2)<pG->dx2) &&(fabs(x3)<pG->dx3)){
    Press = pG->U[k][j][i].E - 0.5/pG->U[k][j][i].d*(SQR(pG->U[k][j][i].M1)
                             + SQR(pG->U[k][j][i].M2)+ SQR(pG->U[k][j][i].M3));
    Press = Press*Gamma_1;
    Press = Press*pG->units.Dcode*SQR(pG->units.Vcode)/consts.kB;

    return Press;
  } else return 0.;
}

static Real hst_ncenter(const GridS *pG, const int i, const int j, const int k)
{
  Real r,x1,x2,x3;
  cc_pos(pG,i,j,k,&x1,&x2,&x3);

  r=sqrt(SQR(x1)+SQR(x2)+SQR(x3));
  
  if(r <= rin)
//  if((fabs(x1)<pG->dx1) && (fabs(x2)<pG->dx2) &&(fabs(x3)<pG->dx3))
    return pG->U[k][j][i].d;
  else return 0.;
}

static Real hst_Vcenter(const GridS *pG, const int i, const int j, const int k)
{
  Real r,x1,x2,x3;
  cc_pos(pG,i,j,k,&x1,&x2,&x3);

  r=sqrt(SQR(x1)+SQR(x2)+SQR(x3));
  
  if(r <= rin)
//  if((fabs(x1)<pG->dx1) && (fabs(x2)<pG->dx2) &&(fabs(x3)<pG->dx3))
    return 1.;
  else return 0.;
}

static Real hst_Thot(const GridS *pG, const int i, const int j, const int k)
{
  Real r,x1,x2,x3;
  Real Press,Temp;

  cc_pos(pG,i,j,k,&x1,&x2,&x3);
  r=sqrt(SQR(x1)+SQR(x2)+SQR(x3));
  
  Press = pG->U[k][j][i].E - 0.5/pG->U[k][j][i].d*(SQR(pG->U[k][j][i].M1)
                           + SQR(pG->U[k][j][i].M2)+ SQR(pG->U[k][j][i].M3));
  Press = Press*Gamma_1;
  Press = Press*pG->units.Dcode*SQR(pG->units.Vcode)/consts.kB;
  Temp = Press/1.1/pG->U[k][j][i].d;

  if(Temp > 2.e4) return Temp;
  else return 0.;
}

static Real hst_Phot(const GridS *pG, const int i, const int j, const int k)
{
  Real r,x1,x2,x3;
  Real Press,Temp;
  cc_pos(pG,i,j,k,&x1,&x2,&x3);

  r=sqrt(SQR(x1)+SQR(x2)+SQR(x3));
  
  Press = pG->U[k][j][i].E - 0.5/pG->U[k][j][i].d*(SQR(pG->U[k][j][i].M1)
                           + SQR(pG->U[k][j][i].M2)+ SQR(pG->U[k][j][i].M3));
  Press = Press*Gamma_1;
  Press = Press*pG->units.Dcode*SQR(pG->units.Vcode)/consts.kB;
  Temp = Press/1.1/pG->U[k][j][i].d;

  if(Temp >= 2.e4) return Press;
  else return 0.;
}

static Real hst_Pshell(const GridS *pG, const int i, const int j, const int k)
{
  Real r,x1,x2,x3;
  Real Press,Temp;
  Real vr = 0.0;

  cc_pos(pG,i,j,k,&x1,&x2,&x3);
  r=sqrt(SQR(x1)+SQR(x2)+SQR(x3));
  if(r != 0.0) vr = (pG->U[k][j][i].M1*x1+pG->U[k][j][i].M2*x2+pG->U[k][j][i].M3*x3)/(pG->U[k][j][i].d*r);

  
  Press = pG->U[k][j][i].E - 0.5/pG->U[k][j][i].d*(SQR(pG->U[k][j][i].M1)
                           + SQR(pG->U[k][j][i].M2)+ SQR(pG->U[k][j][i].M3));
  Press = Press*Gamma_1;
  Press = Press*pG->units.Dcode*SQR(pG->units.Vcode)/consts.kB;
  Temp = Press/1.1/pG->U[k][j][i].d;

  if(Temp < 2.e4 && vr > 1.0) return Press;
  else return 0.;
}

static Real hst_rmom_hot(const GridS *pG, const int i, const int j, const int k)
{
  Real x1,x2,x3;
  Real x10,x20,x30;
  Real rmom,r;

  Real Press,Temp;
  Real vr = 0.0;

 #ifdef STAR_PARTICLE
  StarParListS *pList=NULL;
  StarParS *pStar=NULL;

  pList = pG->Gstars;
  if(pList) {
    pStar = &(pList->starpar);
    x10=pStar->x1;
    x20=pStar->x2;
    x30=pStar->x3;
  }
#else
  x10=0.0;
  x20=0.0;
  x30=0.0;
#endif


  cc_pos(pG,i,j,k,&x1,&x2,&x3);
  r=sqrt(SQR(x1)+SQR(x2)+SQR(x3));
  rmom = pG->U[k][j][i].M1*(x1-x10)/r
        +pG->U[k][j][i].M2*(x2-x20)/r
        +pG->U[k][j][i].M3*(x3-x30)/r;
  if(r != 0.0) vr = (pG->U[k][j][i].M1*x1+pG->U[k][j][i].M2*x2+pG->U[k][j][i].M3*x3)/(pG->U[k][j][i].d*r);

  Press = pG->U[k][j][i].E - 0.5/pG->U[k][j][i].d*(SQR(pG->U[k][j][i].M1)
                           + SQR(pG->U[k][j][i].M2)+ SQR(pG->U[k][j][i].M3));
  Press = Press*Gamma_1;
  Press = Press*pG->units.Dcode*SQR(pG->units.Vcode)/consts.kB;
  Temp = Press/1.1/pG->U[k][j][i].d;

  if(Temp >= 2.e4) return rmom; 
  else return 0.;
}

static Real hst_rmom_shell(const GridS *pG, const int i, const int j, const int k)
{
  Real x1,x2,x3;
  Real x10,x20,x30;
  Real rmom,r;

  Real Press,Temp;
  Real vr = 0.0;

 #ifdef STAR_PARTICLE
  StarParListS *pList=NULL;
  StarParS *pStar=NULL;

  pList = pG->Gstars;
  if(pList) {
    pStar = &(pList->starpar);
    x10=pStar->x1;
    x20=pStar->x2;
    x30=pStar->x3;
  }
#else
  x10=0.0;
  x20=0.0;
  x30=0.0;
#endif


  cc_pos(pG,i,j,k,&x1,&x2,&x3);
  r=sqrt(SQR(x1)+SQR(x2)+SQR(x3));
  rmom = pG->U[k][j][i].M1*(x1-x10)/r
        +pG->U[k][j][i].M2*(x2-x20)/r
        +pG->U[k][j][i].M3*(x3-x30)/r;
  if(r != 0.0) vr = (pG->U[k][j][i].M1*x1+pG->U[k][j][i].M2*x2+pG->U[k][j][i].M3*x3)/(pG->U[k][j][i].d*r);

  Press = pG->U[k][j][i].E - 0.5/pG->U[k][j][i].d*(SQR(pG->U[k][j][i].M1)
                           + SQR(pG->U[k][j][i].M2)+ SQR(pG->U[k][j][i].M3));
  Press = Press*Gamma_1;
  Press = Press*pG->units.Dcode*SQR(pG->units.Vcode)/consts.kB;
  Temp = Press/1.1/pG->U[k][j][i].d;

  if(Temp < 2.e4 && vr > 1.0) return rmom; 
  else return 0.;
}



static Real hst_nhot(const GridS *pG, const int i, const int j, const int k)
{
  Real r,x1,x2,x3;
  Real Press,Temp;
  cc_pos(pG,i,j,k,&x1,&x2,&x3);

  r=sqrt(SQR(x1)+SQR(x2)+SQR(x3));
  
  Press = pG->U[k][j][i].E - 0.5/pG->U[k][j][i].d*(SQR(pG->U[k][j][i].M1)
                           + SQR(pG->U[k][j][i].M2)+ SQR(pG->U[k][j][i].M3));
  Press = Press*Gamma_1;
  Press = Press*pG->units.Dcode*SQR(pG->units.Vcode)/consts.kB;
  Temp = Press/1.1/pG->U[k][j][i].d;

  if(Temp > 2.e4) return pG->U[k][j][i].d;
  else return 0.;
}

static Real hst_Vhot(const GridS *pG, const int i, const int j, const int k)
{
  Real r,x1,x2,x3;
  Real Press,Temp;
  cc_pos(pG,i,j,k,&x1,&x2,&x3);

  r=sqrt(SQR(x1)+SQR(x2)+SQR(x3));
  
  Press = pG->U[k][j][i].E - 0.5/pG->U[k][j][i].d*(SQR(pG->U[k][j][i].M1)
                           + SQR(pG->U[k][j][i].M2)+ SQR(pG->U[k][j][i].M3));
  Press = Press*Gamma_1;
  Press = Press*pG->units.Dcode*SQR(pG->units.Vcode)/consts.kB;
  Temp = Press/1.1/pG->U[k][j][i].d;

  if(Temp > 2.e4) return 1.;
  else return 0.;
}


#ifdef STAR_PARTICLE
static Real hst_age_sp(const GridS *pG, const int i, const int j, const int k)
{
  StarParListS *pList=NULL;
  StarParS *pStar=NULL;

  pList = pG->Gstars;
  while (pList) {
    pStar = &(pList->starpar);
    return pStar->age;
    pList = pList->next;
  }
  return 0.;

}


static Real hst_mass_in_sp(const GridS *pG, const int i, const int j, const int k)
{
  StarParListS *pList=NULL;
  StarParS *pStar=NULL;
  int ip,jp,kp;
  Real dVol = pG->dx1*pG->dx2*pG->dx3;
  Real msp=0.;

  pList = pG->Gstars;
  while (pList) {
    pStar = &(pList->starpar);
    
    cc_ijk(pG,pStar->x1,pStar->x2,pStar->x3,&ip,&jp,&kp);

    if((i==ip) && (j==jp) && (k==kp)) msp = pStar->m/dVol;

    pList = pList->next;
  }

  return msp;
}

static Real hst_mass_in_sp_ghost(const GridS *pG, const int i, const int j, const int k)
{
  StarParListS *pList=NULL;
  StarParS *pStar=NULL;
  int ip,jp,kp;
  Real dVol = pG->dx1*pG->dx2*pG->dx3;
  Real msp=0.;

  pList = pG->Gstars_fda;
  while (pList) {
    pStar = &(pList->starpar);
    
    cc_ijk(pG,pStar->x1,pStar->x2,pStar->x3,&ip,&jp,&kp);
    if ((i >= ip-1) && (i <= ip+1) &&
        (j >= jp-1) && (j <= jp+1) &&
        (k >= kp-1) && (k <= kp+1)) {

      msp = pG->U[k][j][i].d;
    }
    pList = pList->next;
  }

  return msp;
}

static Real hst_mass_in_sp_ghost2(const GridS *pG, const int i, const int j, const int k)
{
  StarParListS *pList=NULL;
  StarParS *pStar=NULL;
  int ip,jp,kp;
  Real dVol = pG->dx1*pG->dx2*pG->dx3;
  Real msp=0.;

  pList = pG->Gstars;
  while (pList) {
    pStar = &(pList->starpar);
    
    cc_ijk(pG,pStar->x1,pStar->x2,pStar->x3,&ip,&jp,&kp);

    if((i==ip) && (j==jp) && (k==kp)) msp = pStar->mghost;

    pList = pList->next;
  }

  return msp;
}
#endif

/*----------------------------------------------------------------------------*/
/*! \fn static Real ShearPot(const Real x1, const Real x2,const Real x3)
 *  \brief tidal potential in 3D shearing box
 */
#ifdef SHEARING_BOX
static Real ShearPot(const Real x1, const Real x2, const Real x3)
{
  Real phi=0.0;
#ifndef FARGO
  phi -= qshear*Omega_0*Omega_0*x1*x1;
#endif
  return phi;
}
#endif
#ifdef TEST_STAR
static Real GravPot(const Real x1, const Real x2, const Real x3)
{
  Real phi=0.0,r;
  r=sqrt(SQR(x1)+SQR(x2)+SQR(x3));
  phi -= GM/r;
  return phi;
}
#endif
/*------------------------------------------------------------------------------
 * initialize perturbation.
 * Gaussian random perturbation in Fourier space with power law PS
 */

static void initialize_PS(DomainS *pDomain)
{
  GridS *pGrid=pDomain->Grid;
  int i,j,k;
  int is=pGrid->is, ie=pGrid->ie;
  int js=pGrid->js, je=pGrid->je;
  int ks=pGrid->ks, ke=pGrid->ke;
  int ind;

  Real gdvol, sum1, sum2;
#ifdef MPI_PARALLEL
  int err;
  Real my_sum1,my_sum2;
#endif

  rseed = -1 - (gis + pDomain->Nx[0]*(gjs + pDomain->Nx[1]*gks));

  if ((drho=(Real***)calloc_3d_array(nx3+2*nghost,nx2+2*nghost,nx1+2*nghost,sizeof(Real)))==NULL) {
    ath_error("[problem]: Error allocating memory for vel pert\n");
  }

  /* parameters for spectrum */
  ispect = par_geti("problem","ispect");
  if (ispect == 1) {
    expo = par_getd("problem","expo");
  } else if (ispect == 2) {
    kpeak = par_getd("problem","kpeak")*2.0*PI;
  } else if (ispect > 3) {
    ath_error("Invalid value for ispect\n");
  }
  /* Cutoff wavenumbers of spectrum */
  klow = par_getd("problem","klow"); /* in integer units */
  khigh = par_getd("problem","khigh"); /* in integer units */
  dkx = 2.0*PI/Lx; /* convert k from integer */

  /* Initialize the FFT plan */
  plan = ath_3d_fft_quick_plan(pDomain, NULL, ATH_FFT_BACKWARD);

  /* Allocate memory for FFTs */
  frho = ath_3d_fft_malloc(plan);

  /* Generate new perturbations following appropriate power spectrum */
  pspect(frho);

  /* Transform velocities from k space to physical space */
  ath_3d_fft(plan, frho);

  /* Set the density in real space */
  sum1 = 0.;
  sum2 = 0.;
  gdvol = 1.0/((Real)(gnx1*gnx2*gnx3));
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        ind = OFST(i-is,j-js,k-ks);
        drho[k][j][i] = frho[ind][0]*gdvol;
        sum1 += drho[k][j][i];
        sum2 += SQR(drho[k][j][i]);
      }
    }
  }

#ifdef MPI_PARALLEL
  my_sum1=sum1;
  my_sum2=sum2;
  err = MPI_Allreduce(&my_sum1, &sum1, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  if(err) ath_error("[pspect]: MPI_Allreduce returned error code %d\n",err);
  err = MPI_Allreduce(&my_sum2, &sum2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  if(err) ath_error("[pspect]: MPI_Allreduce returned error code %d\n",err);
#endif
  sum1 = sum1*gdvol;
  sum2 = sum2*gdvol;

  std_dev = sqrt(sum2-sum1*sum1);
}


/*------------------------------------------------------------------------------
 *  pspect -- computes component of velocity with specific power
 *  spectrum in Fourier space determined by ispect
 *
 *  Velocity power spectrum returned in ampl
 *    klow   = multiple of 2 pi/L for cut-off at low  wavenumbers
 *    khigh  = multiple of 2 pi/L for cut-off at high wavenumbers
 *    expo   = exponent of power law
 *    ispect = integer flag which specifies spectrum
 *  Note that the fourier amplitudes are stored in an array with no
 *  ghost zones
 */


static void pspect(ath_fft_data *ampl)
{
  int i,j,k;
  double q1,q2,q3;

  /* set random amplitudes with gaussian deviation */
  for (k=0; k<nx3; k++) {
    for (j=0; j<nx2; j++) {
      for (i=0; i<nx1; i++) {
        q1 = ran2(&rseed);
        q2 = ran2(&rseed);
        q3 = sqrt(-2.0*log(q1+1.0e-20))*cos(2.0*PI*q2);
        q1 = ran2(&rseed);
        ampl[OFST(i,j,k)][0] = q3*cos(2.0*PI*q1);
        ampl[OFST(i,j,k)][1] = q3*sin(2.0*PI*q1);
      }
    }
  }

  /* set power spectrum
   *   ispect=1: power law - original form
   *   ispect=2: form from Gammie&Ostriker
   */
  for (k=0; k<nx3; k++) {
    for (j=0; j<nx2; j++) {
      for (i=0; i<nx1; i++) {
        /* compute k/dkx */
        q3 = KWVM(i,j,k);
        if ((q3 > klow) && (q3 < khigh)) {
          q3 *= dkx; /* multiply by 2 pi/L */
          if (ispect == 1) {
            /* decreasing power law */
            ampl[OFST(i,j,k)][0] /= pow(q3,(expo+2.0)/2.0);
            ampl[OFST(i,j,k)][1] /= pow(q3,(expo+2.0)/2.0);
          } else if (ispect == 2) {
            /* G&O form */
            ampl[OFST(i,j,k)][0] *= pow(q3,3.0)*exp(-4.0*q3/kpeak);
            ampl[OFST(i,j,k)][1] *= pow(q3,3.0)*exp(-4.0*q3/kpeak);
          }
        } else {
          /* introduce cut-offs at klow and khigh */
          ampl[OFST(i,j,k)][0] = 0.0;
          ampl[OFST(i,j,k)][1] = 0.0;
        }
      }
    }
  }
  ampl[0][0] = 0.0;
  ampl[0][1] = 0.0;

  return;
}
/*------------------------------------------------------------------------------
 * ran2: extracted from the Numerical Recipes in C (version 2) code.  Modified
 *   to use doubles instead of floats. -- T. A. Gardiner -- Aug. 12, 2003
 */

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define RNMX (1.0-DBL_EPSILON)

/* Long period (> 2 x 10^{18}) random number generator of L'Ecuyer
 * with Bays-Durham shuffle and added safeguards.  Returns a uniform
 * random deviate between 0.0 and 1.0 (exclusive of the endpoint
 * values).  Call with idum = a negative integer to initialize;
 * thereafter, do not alter idum between successive deviates in a
 * sequence.  RNMX should appriximate the largest floating point value
 * that is less than 1. 

 */
double ran2(long int *idum){
  int j;
  long int k;
  static long int idum2=123456789;
  static long int iy=0;
  static long int iv[NTAB];
  double temp;

  if (*idum <= 0) { /* Initialize */
    if (-(*idum) < 1) *idum=1; /* Be sure to prevent idum = 0 */
    else *idum = -(*idum);
    idum2=(*idum);
    for (j=NTAB+7;j>=0;j--) { /* Load the shuffle table (after 8 warm-ups) */
      k=(*idum)/IQ1;
      *idum=IA1*(*idum-k*IQ1)-k*IR1;
      if (*idum < 0) *idum += IM1;
      if (j < NTAB) iv[j] = *idum;
    }
    iy=iv[0];
  }
  k=(*idum)/IQ1;                 /* Start here when not initializing */
  *idum=IA1*(*idum-k*IQ1)-k*IR1; /* Compute idum=(IA1*idum) % IM1 without */
  if (*idum < 0) *idum += IM1;   /* overflows by Schrage's method */
  k=idum2/IQ2;
  idum2=IA2*(idum2-k*IQ2)-k*IR2; /* Compute idum2=(IA2*idum) % IM2 likewise */
  if (idum2 < 0) idum2 += IM2;
  j=(int)(iy/NDIV);              /* Will be in the range 0...NTAB-1 */
  iy=iv[j]-idum2;                /* Here idum is shuffled, idum and idum2 */
  iv[j] = *idum;                 /* are combined to generate output */
  if (iy < 1) iy += IMM1;
  if ((temp=AM*iy) > RNMX) return RNMX; /* No endpoint values */
  else return temp;
}

#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NTAB
#undef NDIV
#undef RNMX

#ifndef STAR_PARTICLE
void assign_SN(GridS *pGrid){
  Real r,x1,x2,x3,rin0=rin;
  int i,j,k;
  int is=pGrid->is, ie=pGrid->ie;
  int js=pGrid->js, je=pGrid->je;
  int ks=pGrid->ks, ke=pGrid->ke;
  Real Psn,dVol=pGrid->dx1*pGrid->dx2*pGrid->dx3;
  Real SNvol=0.0,Mtot=0.0,M1tot=0.0,M2tot=0.0,M3tot=0.0,einttot=0.0,eint;
  Real Mrad = 1679.;

#ifdef SEDOV
  Real n0,v1r,v2r,v3r,vr,Pr,rhor,xi0=1.15167;
  Real Ush,trad,xi,rhomin=1.e-3;
  Real tanh,fr=1.0,fth=1.0,fv=1.0;
  Real Etot=0., Ektot=0., Ethtot=0.;
#endif
#ifdef EJECTA
  Real Mej=10.,vej,rej,vr,v1r,v2r,v3r;
  vej = 3.17e3*sqrt(10.*ESN/Mej);
#endif


#ifdef MPI_PARALLEL
  Real my_val;
  int ierr;
#endif
  Real sendbuf[6],recvbuf[6];

  while(1){
    sendbuf[0]=0.0;
    sendbuf[1]=0.0;
    sendbuf[2]=0.0;
    sendbuf[3]=0.0;
    sendbuf[4]=0.0;
    sendbuf[5]=0.0;
    for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
    for (i=is; i<=ie; i++) {
      cc_pos(pGrid,i,j,k,&x1,&x2,&x3);
      r=sqrt(SQR(x1)+SQR(x2)+SQR(x3));

      if(r<=rin){
        sendbuf[0] += dVol;
        sendbuf[1] += pGrid->U[k][j][i].d*dVol;
        sendbuf[2] += pGrid->U[k][j][i].M1*dVol;
        sendbuf[3] += pGrid->U[k][j][i].M2*dVol;
        sendbuf[4] += pGrid->U[k][j][i].M3*dVol;
#ifndef BAROTROPIC
        eint = pGrid->U[k][j][i].E-0.5/pGrid->U[k][j][i].d*(SQR(pGrid->U[k][j][i].M1)
                                  +SQR(pGrid->U[k][j][i].M2)+SQR(pGrid->U[k][j][i].M3));
        sendbuf[5] += eint*dVol;
#endif
#ifdef SEDOV
        xi=r/rin;
        get_Sedov_Taylor(xi,&rhor,&vr,&Pr);
        if(r<pGrid->dx1) vr = 0.0; 
        Ektot += 8.0/25.0/(SQR(Gamma)-1)*SQR(xi0)*CUBE(xi0/rin)*rhor*SQR(vr)*dVol;
        Ethtot += 8.0/25.0/(SQR(Gamma)-1)*SQR(xi0)*CUBE(xi0/rin)*Pr*dVol; 
#endif
        }
    }}}

#ifdef MPI_PARALLEL
    ierr = MPI_Allreduce(&(sendbuf[0]), &(recvbuf[0]), 6, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
    &(recvbuf[0]) = &(sendbuf[0]);
#endif
    if(rin == rin0){
      SNvol=recvbuf[0];
      Mtot =recvbuf[1];
      M1tot=recvbuf[2];
      M2tot=recvbuf[3];
      M3tot=recvbuf[4];
      einttot=recvbuf[5];
    }

    Mrad = 1679./pGrid->units.Msun*pow(recvbuf[1]/recvbuf[0],-0.26)*pow(ESN,0.87);
    if (recvbuf[1]/Mrad > fm || fm == 0.0){
      if(rin != rin0) rin -= 0.5*pGrid->dx1;
      break;
    }

    if(myID_Comm_world==0) printf("[problem] Enlarge SN size to rin=%g, Mrad=%g, Mtot=%g\n",rin,fm*Mrad,recvbuf[1]);
 
    SNvol=recvbuf[0];
    Mtot =recvbuf[1];
    M1tot=recvbuf[2];
    M2tot=recvbuf[3];
    M3tot=recvbuf[4];
    einttot=recvbuf[5];
    rin += 0.5*pGrid->dx1;
  }
  Mrad = 1679./pGrid->units.Msun*pow(Mtot/SNvol,-0.26)*pow(ESN,0.87);
  if(myID_Comm_world==0) printf("[problem] SNR is set to rin=%g, Mrad=%g, Mtot=%g, navg=%g\n",rin,fm*Mrad,Mtot,Mtot/SNvol);

#ifdef SEDOV
#ifdef MPI_PARALLEL
  my_val = Ektot;
  ierr = MPI_Allreduce(&my_val, &Ektot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  my_val = Ethtot;
  ierr = MPI_Allreduce(&my_val, &Ethtot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif

  n0 = Mtot/SNvol;
  Etot = (Ektot+Ethtot)/ESN;
  fth = 0.717/Ethtot*ESN;
  fv = sqrt(0.283/Ektot*ESN);
  printf("Etot=%g Ektot=%g Ethtot=%g fth=%g fv=%g\n",Etot,Ektot,Ethtot,fth,fv);
  trad = 0.0057/pGrid->units.Myr*pow(rin/10.,2.5)*sqrt(n0);
  Ush = 0.4*rin/trad;
#else
  Psn = Gamma_1*ESN*1.e51/SNvol/CUBE(pGrid->units.Lcode);
  Psn = Psn/(pGrid->units.Dcode*SQR(pGrid->units.Vcode));
#endif

  for (k=ks; k<=ke; k++) {
  for (j=js; j<=je; j++) {
    for (i=is; i<=ie; i++) {
      cc_pos(pGrid,i,j,k,&x1,&x2,&x3);
      r=sqrt(SQR(x1)+SQR(x2)+SQR(x3));
      if(r<=rin) {
#ifdef SEDOV
        xi=r/rin;
        get_Sedov_Taylor(xi,&rhor,&vr,&Pr);
        rhor = MAX(rhor*n0*(Gamma+1)/Gamma_1,rhomin);
        vr = fv*vr*2.0/(Gamma+1)*Ush;
        Pr = fth*Pr*2.0/(Gamma+1)*n0*SQR(Ush);
        if(r<pGrid->dx1) vr = 0.0; 
        v1r = vr*(x1)/r;
        v2r = vr*(x2)/r;
        v3r = vr*(x3)/r;
        pGrid->U[k][j][i].d = rhor;
        pGrid->U[k][j][i].M1 = pGrid->U[k][j][i].d*v1r;
        pGrid->U[k][j][i].M2 = pGrid->U[k][j][i].d*v2r;
        pGrid->U[k][j][i].M3 = pGrid->U[k][j][i].d*v3r;
#else
        pGrid->U[k][j][i].d = Mtot/SNvol;
        pGrid->U[k][j][i].M1 = pGrid->U[k][j][i].d*M1tot/Mtot;
        pGrid->U[k][j][i].M2 = pGrid->U[k][j][i].d*M2tot/Mtot;
        pGrid->U[k][j][i].M3 = pGrid->U[k][j][i].d*M3tot/Mtot;
#ifdef EJECTA
        pGrid->U[k][j][i].d = Mej/pGrid->units.Msun/SNvol;
        vr = sqrt(5./3.)*vej*(r/rin);
	v1r = vr*x1/rin;
	v2r = vr*x2/rin;
	v3r = vr*x3/rin;
        pGrid->U[k][j][i].M1 += pGrid->U[k][j][i].d*v1r;
        pGrid->U[k][j][i].M2 += pGrid->U[k][j][i].d*v2r;
        pGrid->U[k][j][i].M3 += pGrid->U[k][j][i].d*v3r;
        pGrid->U[k][j][i].d += Mtot/SNvol;
#endif /* EJECTA */
#endif

#ifndef BAROTROPIC
	pGrid->U[k][k][i].E = einttot/SNvol;
#ifdef SEDOV
        pGrid->U[k][j][i].E += Pr/Gamma_1;
#else
#ifndef EJECTA
        pGrid->U[k][j][i].E += Psn/Gamma_1;
#endif
#endif
        pGrid->U[k][j][i].E += 0.5/pGrid->U[k][j][i].d
                             *(SQR(pGrid->U[k][j][i].M1) 
                             + SQR(pGrid->U[k][j][i].M2)
                             + SQR(pGrid->U[k][j][i].M3));

#endif
      }
    }
  }}
  return;
}

void get_Sedov_Taylor(Real xi, Real *alpha, Real *v, Real *p){
  int ixi,ixip1;
  Real dxi;
/* ST solution from Shu's book II, Eqs. (17.10)-(17.12) */
/*    xi       alpha       xi*v      xi^2*p             */
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


#endif
