#include "copyright.h"
/*==============================================================================
 * FILE: starpar_ti.c
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
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"

#ifndef OPERATOR_SPLIT_COOLING
#error This problem only works with --enable-cooling option
#endif

#ifdef MHD
#error This problem does not work with MHD
#endif

/* FFT indexing Nfast=k, Nmid=j, Nslow=i (opposite to Athena)
 *  * For OFST, i,j,k,nx2,nx3 reference the local grid */
#define OFST(i, j, k) ((k) + nx3*((j) + nx2*(i)))
/* KWVM: magnitude of wavenumber k in units of dkx */
#define KWVM(i, j, k) (sqrt(SQR(KCOMP(i,gis,gnx1))+ \
                            SQR(KCOMP(j,gjs,gnx2))+SQR(KCOMP(k,gks,gnx3))))

#define VERTICAL_GRAVITY

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
static Real Punit, dmax;
static Real n0,T0,P0;
static Real ngbc=100.;
static Real ton;

/*==============================================================================
 * History function:
 *============================================================================*/
static Real hst_dmax(const GridS *pG, const int i, const int j, const int k);
#ifdef STAR_PARTICLE
static Real hst_mass_in_sp(const GridS *pG, const int i, const int j, const int k);
static Real hst_mass_in_sp_ghost(const GridS *pG, const int i, const int j, const int k);
static Real hst_mass_diff_in_sp_ghost(const GridS *pG, const int i, const int j, const int k);
#endif
static Real hst_x2_dke(const GridS *pG, const int i, const int j, const int k);
static Real hst_scaleH2(const GridS *pG, const int i, const int j, const int k);
static Real hst_Mw(const GridS *pG, const int i, const int j, const int k);
static Real hst_Mu(const GridS *pG, const int i, const int j, const int k);
static Real hst_Mc(const GridS *pG, const int i, const int j, const int k);
static Real hst_Mgbc(const GridS *pG, const int i, const int j, const int k);

/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 * ran2()
 *============================================================================*/
#ifdef SHEARING_BOX
static Real ShearPot(const Real x1, const Real x2, const Real x3);
#endif
static void initialize_PS(DomainS *pDomain);
static void pspect(ath_fft_data *ampl);

static void initialize(DomainS *pD);
static double ran2(long int *idum);
static Real logd(const GridS *pG, const int i, const int j, const int k);

#ifdef VERTICAL_GRAVITY
static Real VertGravPot(const Real x1, const Real x2, const Real x3);
static Real gvert, scaleH;
#endif
/*=========================== PUBLIC FUNCTIONS ===============================*/
/*----------------------------------------------------------------------------*/
/* problem:  */

void problem(DomainS *pDomain)
{
  GridS *pGrid=(pDomain->Grid);
  int i=0,j=0,k=0;
  int is,ie,js,je,ks,ke;
  int iprob;
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

  initialize(pDomain);
/* Read problem parameters */


  iprob = par_geti_def("problem","iprob",2);
  if(myID_Comm_world == 0) printf("P0/k0 = %g, P0 in code unit = %g\n",1.1*n0*T0,P0);
  if(iprob == 3){
#ifdef FFT_ENABLED
    if ((drho=(Real***)calloc_3d_array(nx3+2*nghost,nx2+2*nghost,nx1+2*nghost,sizeof(Real)))==NULL)      ath_error("[problem]: Error allocating memory for vel pert\n");
    rseed = -1 - (gis + pDomain->Nx[0]*(gjs + pDomain->Nx[1]*gks));
    initialize_drho(pDomain,drho,&std_dev,rseed);
#else
    ath_error("[problem]: --enable-fft is needed to generate PS\n");
#endif
  }

  if(iprob == 4){
#ifdef FFT_ENABLED
    if ((dv1=(Real***)calloc_3d_array(nx3+2*nghost,nx2+2*nghost,nx1+2*nghost,sizeof(Real)))==NULL)      ath_error("[problem]: Error allocating memory for vel pert\n");
    if ((dv2=(Real***)calloc_3d_array(nx3+2*nghost,nx2+2*nghost,nx1+2*nghost,sizeof(Real)))==NULL)      ath_error("[problem]: Error allocating memory for vel pert\n");
    if ((dv3=(Real***)calloc_3d_array(nx3+2*nghost,nx2+2*nghost,nx1+2*nghost,sizeof(Real)))==NULL)      ath_error("[problem]: Error allocating memory for vel pert\n");
    rseed = -1 - (gis + pDomain->Nx[0]*(gjs + pDomain->Nx[1]*gks));
    initialize_dv(pDomain,dv1,dv2,dv3,rseed,1); // 1 for projection to Div V =0
#else
    ath_error("[problem]: --enable-fft is needed to generate PS\n");
#endif
  }

/* Constant density and temperature initially */
  dmax=1.e-30;

  for (k=ks; k<=ke; k++) {
  for (j=js; j<=je; j++) {
    for (i=is; i<=ie; i++) {
      cc_pos(pGrid,i,j,k,&x1,&x2,&x3);

      if(iprob == 1) pGrid->U[k][j][i].d = n0*(1+amp*sin(kwx*x1));
      if(iprob == 2) pGrid->U[k][j][i].d = n0*(1+amp*(ran2(&iseed)-0.5));
      if(iprob == 3) pGrid->U[k][j][i].d = n0*(1+amp*(drho[k][j][i]/std_dev));
      if(iprob == 4) pGrid->U[k][j][i].d = n0;
#ifdef VERTICAL_GRAVITY
//      pGrid->U[k][j][i].d *= exp(-0.5*SQR(x3/scaleH));
#endif
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

#ifdef MPI_PARALLEL
  Real my_max;
  int ierr;
#endif

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

  Real gfact;

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
        for (k=pG->ks; k<=pG->ke; k++) {
          for (j=pG->js; j<=pG->je; j++) {
            for (i=pG->is; i<=pG->ie; i++) {
              dmax = MAX(dmax,pG->U[k][j][i].d);
            }
          }
        }

#ifdef SELF_GRAVITY
        gfact = 1.0/(1.0+exp(-4.0*(pM->time*Omega_0-(PI+ton*2*PI))));
        if((pM->time*Omega_0 > (ton+1)*2*PI) || (ton < 0)) 
          four_pi_G = 4.0*PI*pG->units.G; 
        else 
          four_pi_G = 4.0*PI*pG->units.G*gfact;
#endif

      }
    }
  }

#ifdef MPI_PARALLEL
  my_max = dmax;
  ierr = MPI_Allreduce(&my_max, &dmax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#endif

  return;
}

void Userwork_after_loop(MeshS *pM)
{
  return;
}

/*=========================== PRIVATE FUNCTIONS ==============================*/

static void initialize(DomainS* pD){
  GridS *pG = (pD->Grid);

  Real kappa,rhosd,gfact;

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

  T0    = par_getd("problem","T0");
  n0    = par_getd("problem","n0");
  P0    = 1.1*n0*T0*Punit;

  pG->units.G = consts.G*pG->units.Dcode*SQR(pG->units.Lcode)/SQR(pG->units.Vcode);

#ifdef OPERATOR_SPLIT_COOLING
  pG->heat0 = 2.e-26;
  pG->heat_ratio = 1.0;
#endif

#ifdef THERMAL_CONDUCTION
  kappa = par_getd("problem","kappa");
  kappa_iso = kappa/(consts.kB*pG->units.Lcode*pG->units.Vcode);
  if(myID_Comm_world==0) printf("kappa in c.g.s. = %g, in code unit = %g\n",kappa,kappa_iso);
#endif

#ifdef SHEARING_BOX
  Omega_0 = par_getd_def("problem","Omega",28.e-3);
  qshear  = par_getd_def("problem","qshear",1.0);

  ShearingBoxPot = ShearPot;
#endif

#ifdef SELF_GRAVITY
  ton = par_getd("problem","ton"); // in units of orbit time
  gfact = 1.0/(1.0+exp(-4.0*((pG->time-pG->dt)*Omega_0-(PI+ton*2*PI))));
  if((pG->time*Omega_0 > (ton+1)*2*PI) || (ton < 0)) 
    four_pi_G = 4.0*PI*pG->units.G; 
  else 
   four_pi_G = 4.0*PI*pG->units.G*gfact;
  grav_mean_rho = 0.0;
  if(myID_Comm_world==0) printf("time=%g, Gcons,gfact = %g %g\n",pG->time,pG->units.G,gfact);
#endif

#ifdef VERTICAL_GRAVITY
  rhosd = par_getd("problem","rhosd");
  rhosd *= consts.Msun/pG->units.Mcode;
  gvert = 4.0*PI*pG->units.G*rhosd;
  scaleH = sqrt(P0/n0/gvert);
  if(myID_Comm_world==0) printf("rhosd, gvert in code unit = %g %g, scaleH = %g\n",rhosd,gvert,scaleH);
  StaticGravPot = VertGravPot; 
#else
  par_setd("problem","rhosd","%.15e",0.,"No vertical gravity"); 
#endif

#ifdef STAR_PARTICLE
  tHII = par_getd_def("problem","tHII",3.0);
  tSN = par_getd_def("problem","tSN",4.0);
#endif

/*
  dump_history_enroll(hst_dmax, "dmax");
  dump_history_enroll(hst_x2_dke, "x2dke");
  dump_history_enroll(hst_scaleH2, "H2");
  dump_history_enroll(hst_Mw, "Mw");
  dump_history_enroll(hst_Mu, "Mu");
  dump_history_enroll(hst_Mc, "Mc");
  dump_history_enroll(hst_Mgbc, "Mgbc");

#ifdef STAR_PARTICLE
  dump_history_enroll(hst_mass_in_sp, "msp");
  dump_history_enroll(hst_mass_in_sp_ghost, "mghost");
  dump_history_enroll(hst_mass_diff_in_sp_ghost, "dmghost");
#endif

*/
}

static Real hst_dmax(const GridS *pG, const int i, const int j, const int k)
{
  return dmax;
}

static Real hst_x2_dke(const GridS *pG, const int i, const int j, const int k)
{
  Real x1,x2,x3,dVy;
  cc_pos(pG,i,j,k,&x1,&x2,&x3);
  dVy=pG->U[k][j][i].M2/pG->U[k][j][i].d;
#if defined (SHEARING_BOX) && !defined(FARGO)
  dVy += qshear*Omega_0*x1;
#endif
  return 0.5*SQR(dVy)*pG->U[k][j][i].d;
}

static Real hst_scaleH2(const GridS *pG, const int i, const int j, const int k)
{ 
  Real x1, x2, x3;
  cc_pos(pG,i,j,k,&x1,&x2,&x3);

  return pG->U[k][j][i].d*x3*x3;
}

static Real hst_Mw(const GridS *pG, const int i, const int j, const int k)
{ /* The density fraction if 1> d: warm gas */

  if ( pG->U[k][j][i].d < (pG->heat_ratio)) return pG->U[k][j][i].d;
  else return 0.;
}

static Real hst_Mu(const GridS *pG, const int i, const int j, const int k)
{ /* The density fraction if 1< d < 8.67: unstable gas */

  if ( pG->U[k][j][i].d >= (pG->heat_ratio) && pG->U[k][j][i].d < (8.67*pG->heat_ratio) ) return pG->U[k][j][i].d;
  else return 0.;
}

static Real hst_Mc(const GridS *pG, const int i, const int j, const int k)
{ /* The density fraction if 8.67 <= d < 130: cold gas */

  if ( pG->U[k][j][i].d >= (8.67*pG->heat_ratio) ) return pG->U[k][j][i].d;
  else return 0.;
}

static Real hst_Mgbc(const GridS *pG, const int i, const int j, const int k)
{ /* The density fraction if 130 <= d: self-gravitating gas  */

  if ( pG->U[k][j][i].d >= ngbc ) return pG->U[k][j][i].d;
  else return 0.;
}

#ifdef STAR_PARTICLE
static Real hst_mass_in_sp(const GridS *pG, const int i, const int j, const int k)
{
  StarParListS *pList=NULL;
  StarParS *pStar=NULL;
  int ip,jp,kp;
  Real dVol = pG->dx1*pG->dx2*pG->dx3;
  Real msp = 0;

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
  Real mghost=0.;
  int ip,jp,kp;

  pList = pG->Gstars_fda;
  while (pList) {
    pStar = &(pList->starpar);
    
    cc_ijk(pG,pStar->x1,pStar->x2,pStar->x3,&ip,&jp,&kp);
    if ((i >= ip-NSINK_STARP) && (i <= ip+NSINK_STARP) &&
        (j >= jp-NSINK_STARP) && (j <= jp+NSINK_STARP) &&
        (k >= kp-NSINK_STARP) && (k <= kp+NSINK_STARP)) {
      mghost = pG->U[k][j][i].d;
    }
    pList = pList->next;
  }

  return mghost;
}

static Real hst_mass_diff_in_sp_ghost(const GridS *pG, const int i, const int j, const int k)
{
  StarParListS *pList=NULL;
  StarParS *pStar=NULL;
  int ip,jp,kp;
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
#ifdef VERTICAL_GRAVITY
static Real VertGravPot(const Real x1, const Real x2, const Real x3)
{
  Real phi=0.0;
  phi += 0.5*gvert*x3*x3;
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
