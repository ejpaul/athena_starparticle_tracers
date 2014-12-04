#include "copyright.h"
/*==============================================================================
 * FILE: ti_collapse.c
 *
 * PURPOSE: Problem generator for spherical collapse with cooling.
 *============================================================================*/

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"

#ifdef MHD
#error This problem does not work with MHD
#endif

#define FIXED_BC

static Real nnow,dmax,Pmax,nLP,Punit;
static ConstS consts;
static Real rho_small;

/*==============================================================================
 * History function:
 *============================================================================*/
static Real hst_dmax(const GridS *pG, const int i, const int j, const int k);
static Real hst_Pmax(const GridS *pG, const int i, const int j, const int k);
static Real hst_Mcore(const GridS *pG, const int i, const int j, const int k);

/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 * ran2()
 *============================================================================*/

static double ran2(long int *idum);
static Real logd(const GridS *pG, const int i, const int j, const int k);
Real Teq(const Real Pok, const Real temp);
Real neq(const Real temp);
void BEprof(const Real Tcenter,const Real Pedge, Real *dx, Real *nprof, Real *Pprof);
void derivs(Real x, Real y[], Real dy[]);
Real TLP(const Real Pok, const Real temp);
 

/*=========================== PUBLIC FUNCTIONS ===============================*/
/*----------------------------------------------------------------------------*/
/* problem:  */

void problem(DomainS *pDomain)
{
  GridS *pGrid=(pDomain->Grid);
  int i=0,j=0,k=0;
  int is,ie,js,je,ks,ke;
  int iprob;
  Real Tc,Tw,kappa,pok0,Press;
  Real *nprof=NULL, *Pprof=NULL;
  int ii,NSOLN=100000,converged;
  Real Pedge,Tcenter,dr,f;
  Real TLPeq;
#ifdef MPI_PARALLEL
  Real my_Mtot,my_max;
  int ierr;
#endif


  long int iseed= -1 - myID_Comm_world;


  Real x1,x2,x3;
  Real x1min,x1max,Lx,kwx,amp;
  Real r,rin,Mtot,dvol,M0;
  Real v0,v1=0.,v2=0.,v3=0.,tanh;
  Real rd,rP;

  is = pGrid->is; ie = pGrid->ie;
  js = pGrid->js; je = pGrid->je;
  ks = pGrid->ks; ke = pGrid->ke;

  rho_small=1.e-4;
/* Get box size and set wavenumber */
  x1min = pDomain->MinX[0];
  x1max = pDomain->MaxX[0];
  Lx = x1max - x1min;
  kwx = 2.0*PI/Lx;
  amp   = par_getd("problem","amp");

/* Set Units */
  init_consts(&consts);
  pGrid->units.Dcode = 1.4*consts.mH;
  pGrid->units.Vcode = consts.kms;
  pGrid->units.Lcode = consts.pc;
  pGrid->units.Mcode = pGrid->units.Dcode*CUBE(pGrid->units.Lcode);
// this is time unit not temperature unit
  pGrid->units.Tcode = pGrid->units.Lcode/pGrid->units.Vcode; 
//  pGrid->units.Pcode = pGrid->units.Dcode*SQR(pGrid->units.Vcode);
#ifdef OPERATOR_SPLIT_COOLING
  pGrid->heat0 = 2.e-26;
  pGrid->heat_ratio = 1.0;
#endif

/* Read problem parameters */
  iprob = par_geti("problem","iprob");

  pok0 = par_getd("problem","pok");
  Punit = consts.kB/pGrid->units.Dcode/SQR(pGrid->units.Vcode);
  Press = pok0*Punit;

#ifdef SELF_GRAVITY
  pGrid->units.G = consts.G*pGrid->units.Dcode*SQR(pGrid->units.Lcode)/SQR(pGrid->units.Vcode);
  four_pi_G = 4.0*PI*pGrid->units.G;
  grav_mean_rho = 0.0;
#endif

  nLP = 8.86*1.1*consts.kB/SQR(1.4*consts.mH)/PI/consts.G/SQR(pGrid->dx1*consts.pc);
  converged = bisection(TLP,1.,184.,pok0,&TLPeq);
  nLP = neq(TLPeq);
  if(myID_Comm_world == 0) printf("nLP = %g for dx = %g\n",nLP,pGrid->dx1);

  converged = bisection(Teq,1.,184.,pok0,&Tc);
  converged = bisection(Teq,5000.,10000.,pok0,&Tw);
  if(myID_Comm_world == 0) printf("T_c, T_w = %g %g, M_0 = %g at P/k = %g\n",Tc,Tw,M0,pok0);

  M0 = SQR(1.1*consts.kB*Tc/(pGrid->units.Dcode))/sqrt(4.0*PI*CUBE(consts.G)*pok0*consts.kB)/consts.Msun;
  if(myID_Comm_world == 0) printf("Mcode = %g Msun\n",pGrid->units.Mcode/consts.Msun);

#ifdef ISOTHERMAL
    Tcenter = par_getd("problem","Tcenter");
    Iso_csound2 = (1.1*consts.kB/(1.4*consts.mH)*Tcenter)/SQR(pGrid->units.Vcode);
    Iso_csound = sqrt(Iso_csound2);
    if(myID_Comm_world == 0) printf("Tcenter = %g, c_s = %g\n",Tcenter, Iso_csound);
#endif

  if(iprob == 2){
    Pedge = pok0;
    Tcenter = par_getd("problem","Tcenter");

    if ((nprof = (Real*)calloc_1d_array(NSOLN,sizeof(Real))) == NULL)
      ath_error("[ti_collapse]: Error allocating memory for 1D analytic solution\n");
    if ((Pprof = (Real*)calloc_1d_array(NSOLN,sizeof(Real))) == NULL)
      ath_error("[ti_collapse]: Error allocating memory for 1D analytic solution\n");

    dr = 1.e-3;
    BEprof(Tcenter,Pedge,&dr,nprof,Pprof);
  }

  if(iprob == 3){
    v0 = par_getd("problem","v0");
  }

  rin = par_getd_def("problem","rin",0.0);
  Mtot = 0.0;
  dvol = pGrid->dx1*pGrid->dx2*pGrid->dx3;

#ifdef THERMAL_CONDUCTION
  kappa = par_getd("problem","kappa");
  kappa_iso = kappa/(consts.kB*pGrid->units.Lcode*pGrid->units.Vcode);
  if(myID_Comm_world==0) printf("kappa in c.g.s. = %g, in code unit = %g\n",kappa,kappa_iso);
#endif

/* Constant density and temperature initially */

  dmax=1.e-30;
  Pmax=1.e-30;
  for (k=ks-nghost; k<=ke+nghost; k++) {
  for (j=js-nghost; j<=je+nghost; j++) {
    for (i=is-nghost; i<=ie+nghost; i++) {
      cc_pos(pGrid,i,j,k,&x1,&x2,&x3);

      r=sqrt(SQR(x1)+SQR(x2)+SQR(x3));

      if(iprob == 1){ //uniform cold cloud in TE

        rP=Press;
        if(r<rin){
          rd = pok0/1.1/Tc*amp;//(1+amp*(ran2(&iseed)-0.5));
          if((i>=is) && (i<=ie) && 
             (j>=js) && (j<=je) && 
             (k>=ks) && (k<=ke)) Mtot = Mtot+rd*dvol;
        } else {
          rd = pok0/1.1/Tw;
        }
      }
      if(iprob == 2){ //non-isothermal BE sphere
        ii = MIN((int)(r/dr),NSOLN-1);
        f = r/dr-ii;
        rd = (1-f)*nprof[ii]+f*nprof[ii+1];
        rP = ((1-f)*Pprof[ii]+f*Pprof[ii+1])*Punit;

        ii = MIN((int)rint(r/dr),NSOLN-1);
        if(Pprof[ii] > Pedge) rd = nprof[ii]*amp; else rd = nprof[ii];
        rP = Pprof[ii]*Punit;
#ifdef ISOTHERMAL
        if(r>rin) rd = rho_small;
#endif

        if((Pprof[ii] > Pedge) && 
           (i>=is) && (i<=ie) && 
           (j>=js) && (j<=je) && 
           (k>=ks) && (k<=ke)) Mtot = Mtot+rd*dvol;
      }

      if(iprob == 3){ //uniform warm medium with conversing velocity
        rP=Press;
        rd = pok0/1.1/Tw;
	//tanh = 0.5*((exp(0.5*(r-x1max))-exp(-0.5*(r-x1max)))/(exp(0.5*(r-x1max))+exp(-0.5*(r-x1max)))+1);
	tanh = (exp(0.5*(r))-exp(-0.5*(r)))/(exp(0.5*(r))+exp(-0.5*(r)));
	tanh = 1.0;
	if(r!=0){
	  v1 = -v0*tanh*x1/r;
	  v2 = -v0*tanh*x2/r;
	  v3 = -v0*tanh*x3/r;
        } else{
	  v1 = 0.;
	  v2 = 0.;
	  v3 = 0.;
        }
      }

      pGrid->U[k][j][i].d = rd;
      pGrid->U[k][j][i].M1 = rd*v1;
      pGrid->U[k][j][i].M2 = rd*v2;
      pGrid->U[k][j][i].M3 = rd*v3;
#ifndef BAROTROPIC
      pGrid->U[k][j][i].E = rP/Gamma_1
             + 0.5*(SQR(pGrid->U[k][j][i].M1) + SQR(pGrid->U[k][j][i].M2)
             + SQR(pGrid->U[k][j][i].M3))/pGrid->U[k][j][i].d;
#endif

      dmax = MAX(dmax,rd);
      Pmax = MAX(Pmax,rP/Punit);
    }
  }}

#ifdef MPI_PARALLEL
  my_Mtot = Mtot;
  ierr = MPI_Allreduce(&my_Mtot, &Mtot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  my_max = dmax;
  ierr = MPI_Allreduce(&my_max, &dmax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  my_max = Pmax;
  ierr = MPI_Allreduce(&my_max, &Pmax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#endif

  Mtot=Mtot*pGrid->units.Mcode/consts.Msun;

  if(myID_Comm_world==0) printf("Mtot = %g Msun %g M0\n",Mtot,Mtot/M0);
  dump_history_enroll(hst_dmax, "dmax");
  dump_history_enroll(hst_Pmax, "Pmax");
  dump_history_enroll(hst_Mcore, "Mcore");

#ifdef FIXED_BC
  bvals_mhd_fun(pDomain,left_x1,do_nothing_bc);
  bvals_mhd_fun(pDomain,right_x1,do_nothing_bc);
  bvals_mhd_fun(pDomain,left_x2,do_nothing_bc);
  bvals_mhd_fun(pDomain,right_x2,do_nothing_bc);
  bvals_mhd_fun(pDomain,left_x3,do_nothing_bc);
  bvals_mhd_fun(pDomain,right_x3,do_nothing_bc);
#endif


  free_1d_array(nprof);
  free_1d_array(Pprof);

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
  return;
}

static Real logd(const GridS *pG, const int i, const int j, const int k)
{
  return log10(pG->U[k][j][i].d);
}

static Real hst_Pmax(const GridS *pG, const int i, const int j, const int k)
{
  return Pmax;
}

static Real hst_dmax(const GridS *pG, const int i, const int j, const int k)
{
  return dmax;
}

static Real hst_Mcore(const GridS *pG, const int i, const int j, const int k)
{
   
  if (pG->U[k][j][i].d > 5.0) return pG->U[k][j][i].d*pG->units.Mcode/consts.Msun;
  return 0.0;
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
  PrimS W;
  Real x1,x2,x3,r;
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
        Pmax = 1.e-30;
        for (k=pG->ks; k<=pG->ke; k++) {
          for (j=pG->js; j<=pG->je; j++) {
            for (i=pG->is; i<=pG->ie; i++) {
              cc_pos(pG,i,j,k,&x1,&x2,&x3);
              r=sqrt(SQR(x1)+SQR(x2)+SQR(x3));
              dmax = MAX(dmax,pG->U[k][j][i].d);
#ifndef BAROTROPIC
              W = Cons_to_Prim(&pG->U[k][j][i]);
              Pmax = MAX(Pmax,W.P/Punit);
#else
              Pmax = dmax*Iso_csound2/Punit;
#endif

            }
          }
        }

#ifdef ISOTHERMAL
        for (k=pG->ks-nghost; k<=pG->ke+nghost; k++) {
          for (j=pG->js-nghost; j<=pG->je+nghost; j++) {
            for (i=pG->is-nghost; i<=pG->ie+nghost; i++) {
              if (pG->U[k][j][i].d < rho_small) {
                pG->U[k][j][i].d = rho_small;
                pG->U[k][j][i].M1 = 0.0;
                pG->U[k][j][i].M2 = 0.0;
                pG->U[k][j][i].M3 = 0.0;
              }
            }
          }
        }
#endif
      }
    }
  }
  

#ifdef MPI_PARALLEL
  my_max = dmax;
  ierr = MPI_Allreduce(&my_max, &dmax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  my_max = Pmax;
  ierr = MPI_Allreduce(&my_max, &Pmax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#endif

  if(dmax > 2.0*nLP){
#ifndef STAR_PARTICLE
    ath_error("[problem] Singularity has been reached!\n");
#else
    data_output(pM, 1); 
    ath_pout(0,"[problem] additional data dump has been made at t=%g\n",pM->time);
    if(myID_Comm_world == 0) printf("[problem] additional data dump has been made at t=%g\n",pM->time);
#endif
  }

  return;
}

void Userwork_after_loop(MeshS *pM)
{
  return;
}

/*=========================== PRIVATE FUNCTIONS ==============================*/
/*------------------------------------------------------------------------------
 * Teq: return thermal equilibrium density and temperature for a given pressure.
 * cooling: Koyama & Inutsuka cooling function for low temperature.
 * 12/30/2013 by Chang-Goo Kim
 */
Real Teq(const Real Pok, const Real temp){
  Real nden;
 
  nden=Pok/1.1/temp;

  return (nden/neq(temp)-1.0);
}

Real TLP(const Real Pok, const Real temp){
  return neq(temp) - nLP*temp;
}

Real neq(const Real temp){
  Real c1=1.e7, c2=1.4e-2, t1=1.184e5, t2=92.;
 
  return 1.0/(c2*sqrt(temp)*exp(-t2/temp)+c1*exp(-t1/(temp+1.e3)));
}

void BEprof(const Real Tcenter,const Real Pedge, Real *dx, Real *nprof, Real *Pprof){
  Real ncenter,pokcenter,rcenter;
  Real y[2],eps,r0,r1,dr;
  int converged,i,NSOLN=100000;
  Real Tnow,pok;

  ncenter = neq(Tcenter);
  pokcenter = 1.1*ncenter*Tcenter;
  rcenter = sqrt(1.1*consts.kB/(4.0*PI*consts.G))/(1.4*consts.mH)*sqrt(Tcenter/ncenter)/consts.pc;
  if(myID_Comm_world == 0) printf("T, n, P/k, r/pc = %g %g %g %g\n",Tcenter,ncenter,pokcenter,rcenter);

  dr = (*dx);
  eps = 1.e-8;
  r0 = eps;
  r1 = r0 + dr;

  y[0] = 1.0;
  y[1] = 0.0;
  for (i=0; i<NSOLN; i++) {
    pok = y[0]*pokcenter;
#ifndef BAROTROPIC
    if(pok < Pedge){
      converged = bisection(Teq,5000.,10000.,pok,&Tnow);
    }else{
      converged = bisection(Teq,1.,184.,pok,&Tnow);
    }
#else
    Tnow = Tcenter;
#endif
    nnow= pok/Tnow/1.1/ncenter;
    Pprof[i] = pok;
    nprof[i] = pok/Tnow/1.1;

    odeint_lite(y, 2, r0, r1, eps, dr, 0.0, derivs);
    r0 = r1;
    r1 = r0 + dr;

//    if(pok < Pedge) break;
//    else{
//      odeint_lite(y, 2, r0, r1, eps, dr, 0.0, derivs);
//      r0 = r1;
//      r1 = r0 + dr;
//    }
  }
  *dx = dr*rcenter;
}

void derivs(Real x, Real y[], Real dy[]){
  dy[0] = nnow*y[1]/SQR(x);
  dy[1] = -nnow*SQR(x);
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
