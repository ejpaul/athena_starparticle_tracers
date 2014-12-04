#include "copyright.h"
/*==============================================================================
 * FILE: peturb.c
 *
 * PURPOSE: generator for gaussian random field
 *============================================================================*/

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"

#ifdef FFT_ENABLED /* ENDIF AT END OF FILE */

/* FFT indexing Nfast=k, Nmid=j, Nslow=i (opposite to Athena)
 *  * For OFST, i,j,k,nx2,nx3 reference the local grid */
#define OFST(i, j, k) ((k) + nx3*((j) + nx2*(i)))
/* KWVM: magnitude of wavenumber k in units of dkx */
#define KWVM(i, j, k) (sqrt(SQR(KCOMP(i,gis,gnx1))+ \
                            SQR(KCOMP(j,gjs,gnx2))+SQR(KCOMP(k,gks,gnx3))))

void pspect(ath_fft_data *ampl,long int *rseed);
void project(ath_fft_data *fv1,ath_fft_data *fv2,ath_fft_data *fv3);
void initialize(DomainS *pD);
double ran2(long int *idum);

/* variables for generate gaussan random field */
/* Number of cells in local grid, number of cells in global grid */
static int nx1,nx2,nx3,gnx1,gnx2,gnx3;
/* Starting indices for global grid */
static int gis,gjs,gks;
/* Seed for random number generator */
/* Driving properties */
static int ispect;

/* Cutoff wavenumbers, G&O spect peak, power law spect exponent, 2 pi/L */
static Real klow,khigh,kpeak,expo,dkx;

/* FFTW - Variables, Plan, etc. */
/* These are made static global variables so that they need not be
 *    allocated AND destroyed with each call to pspect! */
static struct ath_3d_fft_plan *plan;
/* Between calls to generate(), these have unshifted, unnormalized
 *  * velocity perturbations. */

/*------------------------------------------------------------------------------
 * initialize perturbation.
 * Gaussian random perturbation in Fourier space with power law PS
 */
void initialize(DomainS *pD)
{
  GridS *pG=pD->Grid;

  int i,j,k;
  int ind;

  Real x1min, x1max, Lx;

#ifndef FFT_ENABLED
  ath_error("[problem]: --enable-fft is needed to generate PS\n");
#endif


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

  x1min = pD->MinX[0];
  x1max = pD->MaxX[0];
  Lx = x1max - x1min;

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

  return;
}


void initialize_drho(DomainS *pD, Real ***drho, Real *std_dev, long int rseed)
{
  GridS *pG=pD->Grid;

  int i,j,k;
  int is=pG->is, ie=pG->ie;
  int js=pG->js, je=pG->je;
  int ks=pG->ks, ke=pG->ke;
  int ind;

  ath_fft_data *frho=NULL;

  Real gdvol, sum1, sum2;
#ifdef MPI_PARALLEL
  int err;
  Real my_sum1,my_sum2;
#endif

#ifndef FFT_ENABLED
  ath_error("[problem]: --enable-fft is needed to generate PS\n");
#endif

  initialize(pD);

  /* Initialize the FFT plan */
  plan = ath_3d_fft_quick_plan(pD, NULL, ATH_FFT_BACKWARD);

  /* Allocate memory for FFTs */
  frho = ath_3d_fft_malloc(plan);

  /* Generate new perturbations following appropriate power spectrum */
  pspect(frho,&rseed);

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

  *std_dev = sqrt(sum2-sum1*sum1);

  free_1d_array(frho);
}

/*------------------------------------------------------------------------------
 * initialize perturbation.
 * Gaussian random perturbation in Fourier space with power law PS
 */

void initialize_dv(DomainS *pD, Real ***dv1, Real ***dv2, Real ***dv3,
                   long int rseed, int iproj)
{
  GridS *pG=pD->Grid;

  int i,j,k;
  int is=pG->is, ie=pG->ie;
  int js=pG->js, je=pG->je;
  int ks=pG->ks, ke=pG->ke;

  Real dvol;
  int ind;

  ath_fft_data *fv1=NULL, *fv2=NULL, *fv3=NULL;

#ifndef FFT_ENABLED
  ath_error("[problem]: --enable-fft is needed to generate PS\n");
#endif

  initialize(pD);

  /* Initialize the FFT plan */
  plan = ath_3d_fft_quick_plan(pD, NULL, ATH_FFT_BACKWARD);

  /* Allocate memory for FFTs */
  fv1 = ath_3d_fft_malloc(plan);
  fv2 = ath_3d_fft_malloc(plan);
  fv3 = ath_3d_fft_malloc(plan);

  /* Generate new perturbations following appropriate power spectrum */
  pspect(fv1,&rseed);
  pspect(fv2,&rseed);
  pspect(fv3,&rseed);

  /* Require div V = 0 */
  if(iproj) project(fv1,fv2,fv3); 

  /* Transform velocities from k space to physical space */
  ath_3d_fft(plan, fv1);
  ath_3d_fft(plan, fv2);
  ath_3d_fft(plan, fv3);

  /* Set the velocities in real space */
  dvol = 1.0/((Real)(gnx1*gnx2*gnx3));
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        ind = OFST(i-is,j-js,k-ks);
        dv1[k][j][i] = fv1[ind][0]*dvol;
        dv2[k][j][i] = fv2[ind][0]*dvol;
        dv3[k][j][i] = fv3[ind][0]*dvol;
      }
    }
  }

  free_1d_array(fv1);
  free_1d_array(fv2);
  free_1d_array(fv3);

  return;
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


void pspect(ath_fft_data *ampl,long int *rseed)
{
  int i,j,k;
  double q1,q2,q3;

  /* set random amplitudes with gaussian deviation */
  for (k=0; k<nx3; k++) {
    for (j=0; j<nx2; j++) {
      for (i=0; i<nx1; i++) {
        q1 = ran2(rseed);
        q2 = ran2(rseed);
        q3 = sqrt(-2.0*log(q1+1.0e-20))*cos(2.0*PI*q2);
        q1 = ran2(rseed);
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

/* ========================================================================== */
/*! \fn static void project()
 *  \brief Makes velocity perturbations divergence free
 */
void project(ath_fft_data *fv1, ath_fft_data *fv2, ath_fft_data *fv3)
{
  int i,j,k,m,ind;
  double kap[3], kapn[3], mag;
  ath_fft_data dot;
  
  /* Project off non-solenoidal component of velocity */
  for (k=0; k<nx3; k++) {
    kap[2] = sin(2.0*PI*(gks+k)/gnx3);
    for (j=0; j<nx2; j++) {
      kap[1] = sin(2.0*PI*(gjs+j)/gnx2);
      for (i=0; i<nx1; i++) {
        if (((gis+i)+(gjs+j)+(gks+k)) != 0) {
          kap[0] = sin(2.0*PI*(gis+i)/gnx1);
          ind = OFST(i,j,k);

          /* make kapn a unit vector */
          mag = sqrt(SQR(kap[0]) + SQR(kap[1]) + SQR(kap[2]));
          for (m=0; m<3; m++) kapn[m] = kap[m] / mag;

          /* find fv_0 dot kapn */
          dot[0] = fv1[ind][0]*kapn[0]+fv2[ind][0]*kapn[1]+fv3[ind][0]*kapn[2];
          dot[1] = fv1[ind][1]*kapn[0]+fv2[ind][1]*kapn[1]+fv3[ind][1]*kapn[2];

          /* fv = fv_0 - (fv_0 dot kapn) * kapn */
          fv1[ind][0] -= dot[0]*kapn[0];
          fv2[ind][0] -= dot[0]*kapn[1];
          fv3[ind][0] -= dot[0]*kapn[2];

          fv1[ind][1] -= dot[1]*kapn[0];
          fv2[ind][1] -= dot[1]*kapn[1];
          fv3[ind][1] -= dot[1]*kapn[2];
        }
      }
    }
  }

  return;
}


void shift_normalize(GridS *pG, Real ***dv1, Real ***dv2, Real ***dv3, Real amp){
  int i,j,k;

  int is=pG->is, ie=pG->ie;
  int js=pG->js, je=pG->je;
  int ks=pG->ks, ke=pG->ke;

  Real t0, t0ij, t0i, t1, t1ij, t1i;
  Real t2, t2ij, t2i, t3, t3ij, t3i;
#ifdef MPI_PARALLEL
  int mpierr;
  Real m[4], gm[4];
#endif

  /* Calculate net momentum pertubation components t1, t2, t3 */
  t0 = 0.0;  t1 = 0.0;  t2 = 0.0;  t3 = 0.0;
  for (k=ks; k<=ke; k++) {
    t0ij = 0.0;  t1ij = 0.0;  t2ij = 0.0;  t3ij = 0.0;
    for (j=js; j<=je; j++) {
      t0i = 0.0;  t1i = 0.0;  t2i = 0.0;  t3i = 0.0;
      for (i=is; i<=ie; i++) {
        t0i += pG->U[k][j][i].d;

	/* The net momentum perturbation */
        t1i += pG->U[k][j][i].d * dv1[k][j][i];
        t2i += pG->U[k][j][i].d * dv2[k][j][i];
        t3i += pG->U[k][j][i].d * dv3[k][j][i];
      }
      t0ij += t0i;  t1ij += t1i;  t2ij += t2i;  t3ij += t3i;
    }
    t0 += t0ij;  t1 += t1ij;  t2 += t2ij;  t3 += t3ij;
  }

#ifdef MPI_PARALLEL
  /* Sum the perturbations over all processors */
  m[0] = t0;  m[1] = t1;  m[2] = t2;  m[3] = t3;
  mpierr = MPI_Allreduce(m, gm, 4, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  if (mpierr) ath_error("[normalize]: MPI_Allreduce error = %d\n", mpierr);
  t0 = gm[0];  t1 = gm[1];  t2 = gm[2];  t3 = gm[3];
#endif /* MPI_PARALLEL */

  /* Subtract the mean velocity perturbation so that the net momentum
   * perturbation is zero. */
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        dv1[k][j][i] -= t1/t0;
        dv2[k][j][i] -= t2/t0;
        dv3[k][j][i] -= t3/t0;
      }
    }
  }

  /* Calculate unscaled velocity dispersions */
  t0 = 0.0;  t1 = 0.0;  t2 = 0.0;  t3 = 0.0;
  for (k=ks; k<=ke; k++) {
    t0ij = 0.0;  t1ij = 0.0;  t2ij = 0.0;  t3ij = 0.0;
    for (j=js; j<=je; j++) {
      t0i = 0.0;  t1i = 0.0;  t2i = 0.0;  t3i = 0.0;
      for (i=is; i<=ie; i++) {
        t0i += pG->U[k][j][i].d;

	/* The net momentum perturbation */
        t1i += pG->U[k][j][i].d * SQR(dv1[k][j][i]);
        t2i += pG->U[k][j][i].d * SQR(dv2[k][j][i]);
        t3i += pG->U[k][j][i].d * SQR(dv3[k][j][i]);
      }
      t0ij += t0i;  t1ij += t1i;  t2ij += t2i;  t3ij += t3i;
    }
    t0 += t0ij;  t1 += t1ij;  t2 += t2ij;  t3 += t3ij;
  }

#ifdef MPI_PARALLEL
  /* Sum the perturbations over all processors */
  m[0] = t0;  m[1] = t1;  m[2] = t2;  m[3] = t3;
  mpierr = MPI_Allreduce(m, gm, 4, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  if (mpierr) ath_error("[normalize]: MPI_Allreduce error = %d\n", mpierr);
  t0 = gm[0];  t1 = gm[1];  t2 = gm[2];  t3 = gm[3];
#endif /* MPI_PARALLEL */

  /* Normalize perturbation amplitude to zero */
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        dv1[k][j][i] *= sqrt(t0/t1);
        dv2[k][j][i] *= sqrt(t0/t2);
        dv3[k][j][i] *= sqrt(t0/t3);
      }
    }
  }

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        pG->U[k][j][i].M1 += amp*pG->U[k][j][i].d*dv1[k][j][i];
        pG->U[k][j][i].M2 += amp*pG->U[k][j][i].d*dv2[k][j][i];
        pG->U[k][j][i].M3 += amp*pG->U[k][j][i].d*dv3[k][j][i];
#ifndef BAROTROPIC
        pG->U[k][j][i].E += 0.5*pG->U[k][j][i].d*SQR(amp*dv1[k][j][i]);
        pG->U[k][j][i].E += 0.5*pG->U[k][j][i].d*SQR(amp*dv2[k][j][i]);
        pG->U[k][j][i].E += 0.5*pG->U[k][j][i].d*SQR(amp*dv3[k][j][i]);
#endif
      }
    }
  }


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

#endif /* FFT_ENABLED */
