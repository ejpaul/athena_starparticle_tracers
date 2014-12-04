#include "copyright.h"
/*==============================================================================
 * FILE: ti_test.c
 *
 * PURPOSE: Problem generator for thermal instability with operator split cooling.
 *          iprob=1 for eigenmode test
 *          iprob=2 for random perturbation
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

static ConstS consts;
static Real Punit, dmax;

/*==============================================================================
 * History function:
 *============================================================================*/
static Real hst_dmax(const GridS *pG, const int i, const int j, const int k);
#ifdef STAR_PARTICLE
static Real hst_mass_in_sp(const GridS *pG, const int i, const int j, const int k);
static Real hst_mass_in_sp_ghost(const GridS *pG, const int i, const int j, const int k);
#endif
/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 * ran2()
 *============================================================================*/

void initialize(DomainS *pD, GridS *pG);
static double ran2(long int *idum);
static Real logd(const GridS *pG, const int i, const int j, const int k);

/*=========================== PUBLIC FUNCTIONS ===============================*/
/*----------------------------------------------------------------------------*/
/* problem:  */

void problem(DomainS *pDomain)
{
  GridS *pGrid=(pDomain->Grid);
  int i=0,j=0,k=0;
  int is,ie,js,je,ks,ke;
  int iprob;
  Real n0,T0,P0;
  long int iseed = -1-myID_Comm_world;

#ifdef MPI_PARALLEL
  Real my_max;
  int ierr;
#endif



  Real x1,x2,x3;
  Real x1min,x1max,Lx,kwx,amp;

  is = pGrid->is; ie = pGrid->ie;
  js = pGrid->js; je = pGrid->je;
  ks = pGrid->ks; ke = pGrid->ke;

/* Get box size and set wavenumber */
  x1min = pDomain->MinX[0];
  x1max = pDomain->MaxX[0];
  Lx = x1max - x1min;
  kwx = 2.0*PI/Lx;
  amp   = par_getd("problem","amp");

  initialize(pDomain,pGrid);
/* Read problem parameters */

  T0    = par_getd("problem","T0");
  n0    = par_getd("problem","n0");
  P0    = 1.1*n0*T0*Punit;

  iprob = par_geti_def("problem","iprob",2);
  if(myID_Comm_world == 0) printf("P0/k0 = %g, P0 in code unit = %g\n",1.1*n0*T0,P0);

/* Constant density and temperature initially */
  dmax=1.e-30;

  for (k=ks; k<=ke; k++) {
  for (j=js; j<=je; j++) {
    for (i=is; i<=ie; i++) {
      cc_pos(pGrid,i,j,k,&x1,&x2,&x3);

      if(iprob == 1) pGrid->U[k][j][i].d = n0*(1+amp*sin(kwx*x1));
      if(iprob == 2) pGrid->U[k][j][i].d = n0*(1+amp*(ran2(&iseed)-0.5));
      pGrid->U[k][j][i].M1 = 0.0;
      pGrid->U[k][j][i].M2 = 0.0;
      pGrid->U[k][j][i].M3 = 0.0;
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

        initialize(pD,pG);

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
  return;
}

void Userwork_after_loop(MeshS *pM)
{
  return;
}

/*=========================== PRIVATE FUNCTIONS ==============================*/

void initialize(DomainS* pD, GridS *pG){

  Real kappa;

  init_consts(&consts);
  pG->units.Dcode = 1.4*consts.mH;
  pG->units.Vcode = consts.kms;
  pG->units.Lcode = consts.pc;
  pG->units.Mcode = pG->units.Dcode*CUBE(pG->units.Lcode);
  pG->units.Tcode = pG->units.Lcode/pG->units.Vcode; 
  Punit = consts.kB/pG->units.Dcode/SQR(pG->units.Vcode);

#ifdef SELF_GRAVITY
  pG->units.G = consts.G*pG->units.Dcode*SQR(pG->units.Lcode)/SQR(pG->units.Vcode);
#endif

#ifdef OPERATOR_SPLIT_COOLING
  pG->heat0 = 2.e-26;
  pG->heat_ratio = 1.0;
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

  dump_history_enroll(hst_dmax, "dmax");
  dump_history_enroll(hst_mass_in_sp, "msp");
  dump_history_enroll(hst_mass_in_sp_ghost, "mghost");
}

static Real hst_dmax(const GridS *pG, const int i, const int j, const int k)
{
  return dmax;
}

#ifdef STAR_PARTICLE
static Real hst_mass_in_sp(const GridS *pG, const int i, const int j, const int k)
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
#endif

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
