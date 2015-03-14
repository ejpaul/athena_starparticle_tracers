#include "copyright.h"
/*==============================================================================
 * FILE: starpar_besphere.c
 * Collapse of a Bonnor-Ebert sphere (Bonnor 1956)
 * By Hao Gong
 * Adapted for testing of tracer following of sink particle active 
 * zone.
*============================================================================*/

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "defs.h"
#include "athena.h"
#include "prototypes.h"
#include "globals.h"

#ifndef ISOTHERMAL
#error Problem generator only works for isothermal EOS
#endif

#ifndef SELF_GRAVITY
#error Problem generator only works for self-gravity
#endif

#ifdef MPI_PARALLEL
#include "mpi.h"
#endif

/* Initial density (will be average density throughout simulation) */
static Real rho_small;
void derivs(Real x, Real y[], Real dy[]);

/* ========================================================================== */
/*
 *  Function problem
 *
 *  Set up initial conditions, allocate memory, and initialize FFT plans
 */

void problem(DomainS *pD)
{
  GridS *pG = pD->Grid;
  int i, is=pG->is, ie = pG->ie;
  int j, js=pG->js, je = pG->je;
  int k, ks=pG->ks, ke = pG->ke;
  int ixs,jxs,kxs;
  Real Gcons;
  Real xx,yy,zz;
  Real xc,yc,zc,r_out,r,dr;
  int ii,NSOLN=10000;
  Real y[2],eps,x0,x1,dx;
  Real dratio,xiunit,rho_crt;
  Real *r_soln=NULL,*rho_soln=NULL;

  /* Parse input file */
  Gcons=par_getd("problem","Gcons");
  xc = par_getd("problem","xc");
  yc = par_getd("problem","yc");
  zc = par_getd("problem","zc");
  r_out = par_getd("problem","r_out"); // edge radius (!= r_crit)
  rho_crt = par_getd("problem","rho_crt"); // density enhancement factor
  dratio=par_getd_def("problem","dratio",17.75); // ratio of central to edge density
  
  /* Set gravity constant*/
  four_pi_G = 4.0*PI*Gcons;
  grav_mean_rho = 0.0;
  rho_small = 1.0e-4;
  
  /* Allocate memory for semi-analytic solution */
  if ((r_soln = (Real*)calloc_1d_array(NSOLN,sizeof(Real))) == NULL)
    ath_error("[starpar_besphere]: Error allocating memory for 1D semi-analytic solution\n");
  if ((rho_soln = (Real*)calloc_1d_array(NSOLN,sizeof(Real))) == NULL)
    ath_error("[starpar_besphere]: Error allocating memory for 1D semi-analytic solution\n");
  
  /* Solve for 1D radial solution by integrating ODE */
  /* x is r, y[0] is exp(-rho/rho_c), y[1] is dy[0]/dx */

/* CGK added and modfied this part to clarify the original code by Hao.
 * BE sphere equation typically adopts a unit system
 * 	[rho] = rho_c
 * 	[r] = c_s/(4*PI*G*rho_c)^(1/2)
 * 	[v] = c_s
 * 	[4*PI*G] = 1
 * Gong & Ostriker (2013) adopts 
 * 	[rho] = rho_e
 * 	[r] = c_s(PI/G/rho_e)^(1/2)
 * 	[v] = c_s
 * 	[4*PI*G] = 4*PI*PI
 * BE dimensionless variables can be obtained by following factors to GO code unit
 * 	dratio = (rho_c/rho_e)
 * 	xiunit = 2*PI*sqrt(rho_c/rho_e)
 * NOTE: The results presented in GO paper can be reproduced when commented out 
 *       #define TRUE_GREENS_FUNCTION_KERNEL in selfg_fft_obc.c, which resluts in
 *       expanding corners via weaker force at large radii. 
 *       Otherwise, corners are collapsing.
 */

  eps = 1.0e-8;
  xiunit=2.0*PI*sqrt(dratio); 
  x0 = eps;  x1 = r_out*xiunit;
  dx = (x1-x0)/(Real)(NSOLN-1);
  y[0] = SQR(x0)/6.0;
  y[1] = x0/3.0;
  r_soln[0] = x0/xiunit;
  rho_soln[0] = rho_crt*exp(-y[0])*dratio;
  for (i=1; i<NSOLN; i++) {
    x1 = x0 + dx;
    odeint_lite(y, 2, x0, x1, eps, dx, 0.0, derivs);
    r_soln[i] = x1/xiunit;
    rho_soln[i] = rho_crt*exp(-y[0])*dratio;
    x0 = x1;
  }

  dr = dx/xiunit;
  /* Initialize density and momenta */
  for (k=ks-nghost; k<=ke+nghost; k++) {
    for (j=js-nghost; j<=je+nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        memset(&(pG->U[k][j][i]),0.0,sizeof(ConsS));
        cc_pos(pG,i,j,k,&xx,&yy,&zz);
        r = sqrt(SQR(xx-xc) + SQR(yy-yc) + SQR(zz-zc));

        if (r <= r_out) {
          ii = MIN((int)rint(r/dr),NSOLN-1);
          pG->U[k][j][i].d = rho_soln[ii];
          pG->U[k][j][i].M1 = rho_soln[ii]*2;
        }
        else {
          pG->U[k][j][i].d = rho_small;
          pG->U[k][j][i].M1 = rho_soln[ii]*2;
        }
      }
    }
  }
  
#ifdef MHD
  /* Initialize uniform magnetic field */
  for (k=ks-nghost; k<=ke+nghost; k++) {
    for (j=js-nghost; j<=je+nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pG->U[k][j][i].B1c  = B0;
        pG->U[k][j][i].B2c  = 0.0;
        pG->U[k][j][i].B3c  = 0.0;
        pG->B1i[k][j][i] = B0;
        pG->B2i[k][j][i] = 0.0;
        pG->B3i[k][j][i] = 0.0;
      }
    }
  }
#endif /* MHD */

  free_1d_array(r_soln);
  free_1d_array(rho_soln);
  
#if defined(MCTRACERS) || defined(VFTRACERS)
  if (strcmp(par_gets("problem","distribution"), "uniform") == 0)
    tracer_init_unif(pGrid);
  else if (strcmp(par_gets("problem","distribution"), "prop") == 0)
    tracer_init_proportional(pGrid);
#endif // TRACERS //

  return;
}


/* ========================================================================== */

/*
 *  Function Userwork_in_loop
 *
 *  Drive velocity field for turbulence in GMC problems
 */

void Userwork_in_loop(MeshS *pM)
{
  DomainS *pD=NULL;
  GridS *pG=NULL;
  int i,j,k;
  int nl,nd;
  
#ifdef STAR_PARTICLE
  for (nl=0; nl<pM->NLevels; nl++){
    for (nd=0; nd<pM->DomainsPerLevel[nl]; nd++){
      if (pM->Domain[nl][nd].Grid != NULL) {
        pD = &(pM->Domain[nl][nd]);
        pG = pM->Domain[nl][nd].Grid;

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
      }
    }
  }
    
#endif /* STAR_PARTICLE */
  
  return;
}

/* ========================================================================== */

void Userwork_after_loop(MeshS *pM)
{
  /* Don't free memory here if doing any analysis because final
   * output hasn't been written yet!! */
  return;
}

void problem_write_restart(MeshS *pM, FILE *fp)
{  return;  }

void problem_read_restart(MeshS *pM, FILE *fp)
{
  return;
}

ConsFun_t get_usr_expr(const char *expr)
{
  return NULL;
}

VOutFun_t get_usr_out_fun(const char *name){
  return NULL;
}

/* ========================================================================== */
/*
 *  Function hst_*
 *
 *  Dumps to history file
 */

/* Dump kinetic energy in perturbations */
static Real hst_dEk(const GridS *pG, const int i, const int j, const int k)
{ /* The kinetic energy in perturbations is 0.5*d*V^2 */
  return 0.5*(pG->U[k][j][i].M1*pG->U[k][j][i].M1 +
              pG->U[k][j][i].M2*pG->U[k][j][i].M2 +
              pG->U[k][j][i].M3*pG->U[k][j][i].M3)/pG->U[k][j][i].d;
}

/* Dump magnetic energy in perturbations */
static Real hst_dEb(const GridS *pG, const int i, const int j, const int k)
{ /* The magnetic energy in perturbations is 0.5*B^2 - 0.5*B0^2 */
#ifdef MHD
  return 0.5*((pG->U[k][j][i].B1c*pG->U[k][j][i].B1c +
               pG->U[k][j][i].B2c*pG->U[k][j][i].B2c +
               pG->U[k][j][i].B3c*pG->U[k][j][i].B3c)-B0*B0);
#else /* MHD */
  return 0.0;
#endif /* MHD */
}

/* -------------------------------------------------------------------------- */

void derivs(Real x, Real y[], Real dy[]) {
  dy[0] = y[1];
  dy[1] = -2.0*y[1]/x + exp(-y[0]);
}


