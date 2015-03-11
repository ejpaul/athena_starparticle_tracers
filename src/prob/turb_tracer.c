#include "copyright.h"
/*==============================================================================
 * FILE: turb_tracer.c
 *
 * Momentum injection to the ISM from radiation-driven shells around star
 *   clusters.  This is the shell formation problem.  This version reads in a
 *   user-defined value of the efficiency, epsilon.  This version begins with
 *   force balance and adds in perturbations.
 *
 * REFERENCES: E. Ostriker, & R. Shetty, "Maximally star-forming galactic
 *  disks I. Starbust regulation via feedback-driven turbulence," accepted
 *  by ApJ. (2011)
 *============================================================================*/

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"
#include "units.h"

#ifndef ISOTHERMAL
#error This problem generator requires --with-eos=isothermal
#endif /* ISOTHERMAL */

static Real rho_small,rho,M_GMC;
static Real M0,tau0,tau_shell,H,Lstar,rstar,rcloud,rcrit,fcrit,tau_cloud,v_turb;
static Real eps_min,eps_max,rho_cloud,L_Jeans,Egrav;
static Real Lx,Ly,Lz;

/* Function prototypes for analysis and outputs */
static Real hst_dEk(const GridS *pG,const int i,const int j,const int k);
static Real hst_dEb(const GridS *pG,const int i,const int j,const int k);
static Real hst_mass_stars(const GridS *pG,const int i,const int j,const int k);
static Real hst_Egrav_tot(const GridS *pG,const int i,const int j,const int k);
static Real hst_Egrav_gas(const GridS *pG,const int i,const int j,const int k);
static Real hst_Egrav_stars(const GridS *pG,const int i,const int j,const int k);
static Real hst_Wgrav_gas(const GridS *pG,const int i,const int j,const int k);
static Real hst_Wgrav_stars(const GridS *pG,const int i,const int j,const int k);
static Real hst_Mr_out(const GridS *pG,const int i,const int j,const int k);
static Real hst_Efree_out(const GridS *pG,const int i,const int j,const int k);
static Real log10d(const GridS *pG,const int i,const int j,const int k);
static Real usr_Sigma1(const GridS *pG,const int i,const int j,const int k);
static Real usr_Sigma2(const GridS *pG,const int i,const int j,const int k);
static Real usr_Sigma3(const GridS *pG,const int i,const int j,const int k);
static Real dEk(const GridS *pG, const int i, const int j, const int k);

/* Uncomment the following define to drive the flow in an impulsive manner
 as was done originally.  Restarts for this mode not yet implemented! */
/* #define IMPULSIVE_DRIVING */

/* KEEP SEMI-COLONS OUT OF THESE PRE-PROCESSOR DIRECTIVES! */
/* FFT indexing Nfast=k, Nmid=j, Nslow=i (opposite to Athena)
 * For OFST, i,j,k,nx2,nx3 reference the local grid */
#define OFST(i, j, k) ((k) + nx3*((j) + nx2*(i)))
/* KWVM: magnitude of wavenumber k in units of dkx */
#define KWVM(i, j, k) (sqrt(SQR(KCOMP(i,gis,gnx1)) + \
SQR(KCOMP(j,gjs,gnx2)) + \
SQR(KCOMP(k,gks,gnx3))))

/* FFTW - Variables, Plan, etc. */
/* These are made static global variables so that they need not be
 allocated AND destroyed with each call to pspect! */
static struct ath_3d_fft_plan *plan;
/* Between calls to generate(), these have unshifted, unnormalized
 * velocity perturbations. */
static ath_fft_data *fv1=NULL, *fv2=NULL, *fv3=NULL;

/* Normalized, shifted velocity perturbations */
static Real ***dv1=NULL, ***dv2=NULL, ***dv3=NULL;
/* Cutoff wavenumbers, G&O spect peak, power law spect exponent, 2 pi/L */
static Real klow,khigh,kpeak,expo,dkx;
/* Energy injection rate, last planned driving time, driving interval.
 * If not using impulsive driving, then the time quantities above are for
 * computing a new spectrum, not driving */
static Real dedt;
/* Number of cells in local grid (with and without ghost zones),
 * number of cells in global grid */
static int nx1,nx2,nx3,nx1gh,nx2gh,nx3gh,gnx1,gnx2,gnx3;
/* Starting and ending indices for local grid */
static int is,ie,il,iu,js,je,jl,ju,ks,ke,kl,ku;
/* Starting and ending indices for global grid */
static int gis,gie,gjs,gje,gks,gke;
/* Spatial extents of local grid */
static Real x1min,x1max,x2min,x2max,x3min,x3max;
/* Length and volume elements */
static Real dx,dV;

#ifdef MCTRACERS
/* Boolean- have tracers been deposited ? */
static int MC_deposit = 0;
#endif /* MCTRACERS */

#ifdef VFTRACERS
/* Boolean- have tracers been deposited ? */
static int VF_deposit = 0;
#endif /* VFTRACERS */

#if defined(MCTRACERS) || defined(VFTRACERS)
/* Time to deposit tracers */
static Real t_dep;
/* Tracer density per cell */
static Real d_prop;
#endif

/* Seed for random number generator */
long int rseed;
#ifdef MHD
/* beta = isothermal pressure / magnetic pressure
 * B0 = sqrt(2.0*Iso_csound2*rhobar/beta) is init magnetic field strength */
static Real beta,B0;
#endif /* MHD */
/* Initial density (will be average density throughout simulation) */
static const Real rhobar = 1.0;

/* Functions appear in this file in the same order that they appear in the
 * prototypes below */

/* Function prototypes for generating velocity perturbations */
static void pspect(ath_fft_data *ampl);
static void project();
static inline void transform();
static inline void generate();
static void perturb(GridS *pG);

/* Function prototypes for initializing and interfacing with Athena */
static void initialize(DomainS *pD);

/* Function prototypes for Numerical Recipes functions */
static Real ran2(long int *idum);

/*=========================== PUBLIC FUNCTIONS ===============================*/
/*----------------------------------------------------------------------------*/
/* problem:  */
void problem(DomainS *pD)
{
    GridS *pG = pD->Grid;
    int i,j,k,ii;
    int in,jn,kn;
    Real x1,x2,x3;
    Real tmp;
    Real rmin,rmax,rstart,eps,F_soln,rijk,dr,sigma,sum;
#ifdef STARPARTICLE
    StarParListS *pList=NULL;
    StarParS *pStar=NULL;
#endif
    Real dfdrcrit,A,B,C;
    Real *r_soln=NULL,*f_soln=NULL;
    Real phi,theta,Ylm_max=1.0,M_sum,M_sum_tot;
    
    rseed = par_geti("problem","rseed");
    if (rseed > 0.0) ath_error("[radpargrav]:  rseed must be <= 0\n");
    initialize(pD);
    
    /* Initialize gas density */
    M_sum = 0.0;
    for (k=kl; k<=ku; k++) {
        for (j=jl; j<=ju; j++) {
            for (i=il; i<=iu; i++) {
                memset(&(pG->U[k][j][i]),0.0,sizeof(ConsS));
                
                pG->U[k][j][i].d = rho;            }
        }
    }
    
    /* Set the initial perturbations.  Note that we're putting in too much
     * energy this time.  This is okay since we're only interested in the
     * saturated state. */
    generate();
    perturb(pG);
    
    /* If decaying turbulence, no longer need the driving memory */
    ath_pout(0,"De-allocating driving memory.\n");
        
    /* Free Athena-style arrays */
    free_3d_array(dv1);
    free_3d_array(dv2);
    free_3d_array(dv3);
        
    /* Free FFTW-style arrays */
    ath_3d_fft_free(fv1);
    ath_3d_fft_free(fv2);
    ath_3d_fft_free(fv3);
    
    return;
}


/*==============================================================================
 * PROBLEM USER FUNCTIONS:
 * problem_write_restart() - writes problem-specific user data to restart files
 * problem_read_restart()  - reads problem-specific user data from restart files
 * get_usr_expr()          - sets pointer to expression for special output data
 * get_usr_out_fun()       - returns a user defined output function pointer
 * Userwork_in_loop        - problem specific work IN     main loop
 * Userwork_after_loop     - problem specific work AFTER  main loop
 *============================================================================*/

void problem_write_restart(MeshS *pM, FILE *fp)
{
    return;
}

void problem_read_restart(MeshS *pM, FILE *fp)
{
    DomainS *pD=NULL;
    GridS *pG=NULL;
    int nl,nd;
    
    for (nl=0; nl<pM->NLevels; nl++) {
        for (nd=0; nd<pM->DomainsPerLevel[nl]; nd++) {
            if (pM->Domain[nl][nd].Grid != NULL) {
                pD = &(pM->Domain[nl][nd]);
                initialize(pD);
            }
        }
    }
    
#ifdef STAR_PARTICLE
    for (nl=0; nl<pM->NLevels; nl++) {
        for (nd=0; nd<pM->DomainsPerLevel[nl]; nd++) {
            if (pM->Domain[nl][nd].Grid != NULL) {
                pD = &(pM->Domain[nl][nd]);
                pG = pD->Grid;
                
                if (myID_Comm_world == 0) {
                    starpar_printlist(0, pG);
                }
//                if (pG->Gstars)
//                    source_exists = 1;  // PUT THIS LINE IN RESTART.C
            }
        }
    }
#endif /* STAR_PARTICLE */
    
    return;
}

ConsFun_t get_usr_expr(const char *expr)
{
    if(strcmp(expr,"Sigma1")==0) return usr_Sigma1;
    if(strcmp(expr,"Sigma2")==0) return usr_Sigma2;
    if(strcmp(expr,"Sigma3")==0) return usr_Sigma3;
    if(strcmp(expr, "dEk")==0) return dEk;
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
    Real newtime;
    const Real rho_floor = 0.01*rho_small;
    
    /* If rho_floor is non-zero, enforce density/momentum floor */
    if (rho_floor) {
        for (nl=0; nl<pM->NLevels; nl++) {
            for (nd=0; nd<pM->DomainsPerLevel[nl]; nd++) {
                if (pM->Domain[nl][nd].Grid != NULL) {
                    pD = &(pM->Domain[nl][nd]);
                    pG = pD->Grid;
                    
                    for (k=kl; k<=ku; k++) {
                        for (j=kl; j<=ju; j++) {
                            for (i=il; i<=iu; i++) {
                                if (pG->U[k][j][i].d < rho_floor) {
                                    pG->U[k][j][i].d = rho_floor;
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
    }
    
#ifdef STAR_PARTICLE
    StarParListS *pGstars = NULL;
    StarParS *pStar = NULL;
    char filename[16];
    char *mode;
    FILE *hp;
    int ip,jp,kp;
    int tmp,mpierr;
    
    for (nl=0; nl<pM->NLevels; nl++) {
        for (nd=0; nd<pM->DomainsPerLevel[nl]; nd++) {
            if (pM->Domain[nl][nd].Grid != NULL) {
                pD = &(pM->Domain[nl][nd]);
                pG = pD->Grid;
                
                /* Write star particle status to output file */
                if (myID_Comm_world == 0) {
                    pGstars = pG->Gstars;
                    
                    while (pGstars) {
                        pStar = &(pGstars->starpar);
                        mode = (pStar->age == 0.0) ? "w" : "a";
#ifdef MPI_PARALLEL
                        sprintf(filename,"../star%4.4d.dat",pStar->id);
#else
                        sprintf(filename,"star%4.4d.dat",pStar->id);
#endif
                        if((hp = fopen(filename,mode)) == NULL) {
                            ath_error("[radpargrav]: Unable to open starpar dump file\n");
                            return;
                        }
                        fprintf(hp,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d\n",
                                pG->time,pStar->m,pStar->x1,pStar->x2,pStar->x3,
                                pStar->v1,pStar->v2,pStar->v3,pStar->age,pStar->mdot,
                                pStar->merge_history);
                        fclose(hp);
                        pGstars = pGstars->next;
                    }
                }
            }
        }
    }
    
#endif /* STAR_PARTICLE */
    
#if defined(MCTRACERS) || defined(VFTRACERS)
    /* if t >= T_f */
    newtime = pG->time + pG->dt;
    if (newtime >= t_dep) {
#ifdef MCTRACERS
        if (MC_deposit == 0) {
            mc_init_threshold(pG, rho*d_prop);
            MC_deposit = 1;
        }
#endif /* MCTRACERS */
#ifdef VFTRACERS
        if (VF_deposit == 0) {
            vf_init_threshold(pG, rho*d_prop);
            VF_deposit = 1;
        }
    }
#endif
    
    if (isnan(pG->dt)) ath_error("[radpargrav]:  Time step is NaN!");
    
    return;
}

void Userwork_in_rad_loop(MeshS *pM)
{
    return;
}

void Userwork_after_loop(MeshS *pM)
{
    Userwork_in_loop(pM);
    
    return;
}

/*==============================================================================
 * PRIVATE FUNCTIONS
 *============================================================================*/

/*------------------------------------------------------------------------------
 *  Function hst_*
 *
 *  Dumps to history file
 *  "History" variables are usually volume integrals: e.g.,
 *    S = \Sum_{ijk} q[k][j][i]*dx1[i]*dx2[j]*dx3[k].  Note that history
 *    variables are TOTAL, VOLUME INTEGRATED quantities, NOT volume averages
 *    (just divide by the Domain volume to compute the latter).  With MPI, the
 *    sum is performed over all Grids in the Domain.
 */

/* Dump kinetic energy in perturbations */
static Real hst_dEk(const GridS *pG,const int i,const int j,const int k)
{
    /* The kinetic energy in perturbations is 0.5*d*V^2 */
    return 0.5*(pG->U[k][j][i].M1*pG->U[k][j][i].M1 +
                pG->U[k][j][i].M2*pG->U[k][j][i].M2 +
                pG->U[k][j][i].M3*pG->U[k][j][i].M3)/pG->U[k][j][i].d;
}

/* Dump magnetic energy in perturbations */
static Real hst_dEb(const GridS *pG,const int i,const int j,const int k)
{
    /* The magnetic energy in perturbations is 0.5*B^2 - 0.5*B0^2 */
#ifdef MHD
    return 0.5*((pG->U[k][j][i].B1c*pG->U[k][j][i].B1c +
                 pG->U[k][j][i].B2c*pG->U[k][j][i].B2c +
                 pG->U[k][j][i].B3c*pG->U[k][j][i].B3c)-B0*B0);
#else /* MHD */
    return 0.0;
#endif /* MHD */
}

/* Dump mass in star particles */
static Real hst_mass_stars(const GridS *pG,const int i,const int j,const int k)
{
#ifdef STARPARTICLE
      StarParListS *pList=NULL;
      StarParS *pStar=NULL;
    
      int ip,jp,kp;
      Real x1,x2,x3;
    
      pList = pG->Lstars;
      while (pList) {
        pStar = &(pList->starpar);
        cc_ijk(pG,pStar->x1,pStar->x2,pStar->x3,&ip,&jp,&kp);
        /* If particle in this cell */
        if (i==ip && j==jp && k==kp) {
          return pStar->m/(pG->dx1*pG->dx2*pG->dx3);
        }
        pList = pList->next;
      }
#endif
    return 0.0;
}

/* Dump gravitational potential energy in gas plus stars */
static Real hst_Egrav_tot(const GridS *pG,const int i,const int j,const int k)
{
#ifdef STARPARTICLE
    StarParListS *pList=pG->Gstars;
    StarParS *pStar=NULL;
#endif
    
    int ip,jp,kp;
    Real tmp,x1,x2,x3,W1[3],W2[3],W3[3];
    
#ifdef STARPARTICLE
    /* The gravitational energy is d*Phi, except where particles are, so check to
     * see if there is a particle. */
    pList = pG->Gstars;
    while (pList) {
        pStar = &(pList->starpar);
        cc_ijk(pG,pStar->x1,pStar->x2,pStar->x3,&ip,&jp,&kp);
        /* If within 1 zone of a star particle */
        if (abs(i-ip)<=1 && abs(j-jp)<=1 && abs(k-kp)<=1) {
            cc_pos(pG,i,j,k,&x1,&x2,&x3);
            /* Compute TSC weights   */
            tmp = 0.5 + (pStar->x1 - x1)/pG->dx1;
            W1[2] = 0.5*SQR(tmp);  W1[0] = W1[2]-tmp+0.5;  W1[1] = 1.0-W1[0]-W1[2];
            tmp = 0.5 + (pStar->x2 - x2)/pG->dx2;
            W2[2] = 0.5*SQR(tmp);  W2[0] = W2[2]-tmp+0.5;  W2[1] = 1.0-W2[0]-W2[2];
            tmp = 0.5 + (pStar->x3 - x3)/pG->dx3;
            W3[2] = 0.5*SQR(tmp);  W3[0] = W3[2]-tmp+0.5;  W3[1] = 1.0-W3[0]-W3[2];
            
            /* Return Phi times weighted particle density */
            tmp = pStar->m/(pG->dx1*pG->dx2*pG->dx3);
            return tmp*pG->Phi[k][j][i]*W1[i-(ip-1)]*W2[j-(jp-1)]*W3[k-(kp-1)];
        }
        pList = pList->next;
    }
#endif
    
    return pG->U[k][j][i].d*pG->Phi[k][j][i];
}

/* Dump gravitational potential energy in gas */
static Real hst_Egrav_gas(const GridS *pG,const int i,const int j,const int k)
{
#ifdef STARPARTICLE
    StarParListS *pList=pG->Gstars;
    StarParS *pStar=NULL;
    int ip,jp,kp;
    
    /* The gravitational energy is d*Phi, except where particles are, so check to
     * see if there is a particle. */
    pList = pG->Gstars;
    while (pList) {
        pStar = &(pList->starpar);
        cc_ijk(pG,pStar->x1,pStar->x2,pStar->x3,&ip,&jp,&kp);
        if (abs(i-ip)<=1 && abs(j-jp)<=1 && abs(k-kp)<=1)
            return 0.0;
        pList = pList->next;
    }
#endif
    return pG->U[k][j][i].d*pG->Phi[k][j][i];
}

/* Dump gravitational potential energy in stars */
static Real hst_Egrav_stars(const GridS *pG,const int i,const int j,const int k)
{
#ifdef STARPARTICLES
    StarParListS *pList=pG->Gstars;
    StarParS *pStar=NULL;
    int ip,jp,kp;
    Real tmp,x1,x2,x3,W1[3],W2[3],W3[3];
    
    /* The gravitational energy is d*Phi, except where particles are, so check to
     * see if there is a particle. */
    pList = pG->Gstars;
    while (pList) {
        pStar = &(pList->starpar);
        cc_ijk(pG,pStar->x1,pStar->x2,pStar->x3,&ip,&jp,&kp);
        /* If within 1 zone of a star particle */
        if (abs(i-ip)<=1 && abs(j-jp)<=1 && abs(k-kp)<=1) {
            cc_pos(pG,i,j,k,&x1,&x2,&x3);
            /* Compute TSC weights   */
            tmp = 0.5 + (pStar->x1 - x1)/pG->dx1;
            W1[2] = 0.5*SQR(tmp);  W1[0] = W1[2]-tmp+0.5;  W1[1] = 1.0-W1[0]-W1[2];
            tmp = 0.5 + (pStar->x2 - x2)/pG->dx2;
            W2[2] = 0.5*SQR(tmp);  W2[0] = W2[2]-tmp+0.5;  W2[1] = 1.0-W2[0]-W2[2];
            tmp = 0.5 + (pStar->x3 - x3)/pG->dx3;
            W3[2] = 0.5*SQR(tmp);  W3[0] = W3[2]-tmp+0.5;  W3[1] = 1.0-W3[0]-W3[2];
            
            /* Return Phi times weighted particle density */
            tmp = pStar->m/(pG->dx1*pG->dx2*pG->dx3);
            return tmp*pG->Phi[k][j][i]*W1[i-(ip-1)]*W2[j-(jp-1)]*W3[k-(kp-1)];
        }
        pList = pList->next;
    }
#endif
    return 0.0;
}

/* Dump rate of work done by gravitational field on gas */
static Real hst_Wgrav_gas(const GridS *pG,const int i,const int j,const int k)
{
    Real x1,x2,x3,f1,f2,f3;
#ifdef STARPARTICLES
    StarParListS *pList=pG->Gstars;
    StarParS *pStar=NULL;
    int ip,jp,kp;
    
    /* The gravitational work is -d*(v dot Phi), except where particles are, so
     * first check to see if there is a particle. */
    pList = pG->Gstars;
    while (pList) {
        pStar = &(pList->starpar);
        cc_ijk(pG,pStar->x1,pStar->x2,pStar->x3,&ip,&jp,&kp);
        if (abs(i-ip)<=1 && abs(j-jp)<=1 && abs(k-kp)<=1)
            return 0.0;
        pList = pList->next;
    }
#endif
    
    /* Otherwise, calculate force using centered-difference approximation */
    f1 = -0.5*(pG->Phi[k][j][i+1]-pG->Phi[k][j][i-1])/pG->dx1;
    f2 = -0.5*(pG->Phi[k][j+1][i]-pG->Phi[k][j-1][i])/pG->dx2;
    f3 = -0.5*(pG->Phi[k+1][j][i]-pG->Phi[k-1][j][i])/pG->dx3;
    return (pG->U[k][j][i].M1*f1 +
            pG->U[k][j][i].M2*f2 +
            pG->U[k][j][i].M3*f3)/pG->U[k][j][i].d;
}

/* Dump rate of work done by gravitational field on stars */
static Real hst_Wgrav_stars(const GridS *pG,const int i,const int j,const int k)
{
#ifdef STARPARTICLES
    StarParListS *pList=pG->Gstars;
    StarParS *pStar=NULL;
    int ip,jp,kp;
    Real tmp,x1,x2,x3,f1,f2,f3,W1[3],W2[3],W3[3];
    
    pList = pG->Gstars;
    while (pList) {
        pStar = &(pList->starpar);
        cc_ijk(pG,pStar->x1,pStar->x2,pStar->x3,&ip,&jp,&kp);
        /* If particle in this cell */
        if (i==ip && j==jp && k==kp) {
            cc_pos(pG,i,j,k,&x1,&x2,&x3);
            /* Compute TSC weights   */
            tmp = 0.5 + (pStar->x1 - x1)/pG->dx1;
            W1[2] = 0.5*SQR(tmp);  W1[0] = W1[2]-tmp+0.5;  W1[1] = 1.0-W1[0]-W1[2];
            tmp = 0.5 + (pStar->x2 - x2)/pG->dx2;
            W2[2] = 0.5*SQR(tmp);  W2[0] = W2[2]-tmp+0.5;  W2[1] = 1.0-W2[0]-W2[2];
            tmp = 0.5 + (pStar->x3 - x3)/pG->dx3;
            W3[2] = 0.5*SQR(tmp);  W3[0] = W3[2]-tmp+0.5;  W3[1] = 1.0-W3[0]-W3[2];
            
            /* Return v dot grad(Phi) times weighted particle density */
            tmp = pStar->m/(pG->dx1*pG->dx2*pG->dx3);
            f1 = -0.5*((pG->Phi[kp][jp][ip  ] - pG->Phi[kp][jp][ip-2])*W1[0] +
                       (pG->Phi[kp][jp][ip+1] - pG->Phi[kp][jp][ip-1])*W1[1] +
                       (pG->Phi[kp][jp][ip+2] - pG->Phi[kp][jp][ip  ])*W1[2])/pG->dx1;
            f2 = -0.5*((pG->Phi[kp][jp  ][ip] - pG->Phi[kp][jp-2][ip])*W2[0] +
                       (pG->Phi[kp][jp+1][ip] - pG->Phi[kp][jp-1][ip])*W2[1] +
                       (pG->Phi[kp][jp+2][ip] - pG->Phi[kp][jp  ][ip])*W2[2])/pG->dx2;
            f3 = -0.5*((pG->Phi[kp  ][jp][ip] - pG->Phi[kp-2][jp][ip])*W3[0] +
                       (pG->Phi[kp+1][jp][ip] - pG->Phi[kp-1][jp][ip])*W3[1] +
                       (pG->Phi[kp+2][jp][ip] - pG->Phi[kp  ][jp][ip])*W3[2])/pG->dx3;
            
            return tmp*(f1*pStar->v1 + f2*pStar->v2 + f3*pStar->v3);
        }
        pList = pList->next;
    }
#endif
    return 0.0;
}

/* Dump radial component of the outward momentum flux */
static Real hst_Mr_out(const GridS *pG,const int i,const int j,const int k)
{
    Real x1,x2,x3,tmp,rmf=0.0;
    
    cc_pos(pG,i,j,k,&x1,&x2,&x3);
    tmp = (x1*pG->U[k][j][i].M1 +
           x2*pG->U[k][j][i].M2 +
           x3*pG->U[k][j][i].M3)/sqrt(SQR(x1) + SQR(x2) + SQR(x3));
    
    /* Inner x1-boundary */
    if (i==pG->is && pG->lx1_id==-1)
        rmf -= tmp*(pG->U[k][j][i].M1/pG->U[k][j][i].d)/pG->dx1;
    
    /* Outer x1-boundary */
    if (i==pG->ie && pG->rx1_id==-1)
        rmf += tmp*(pG->U[k][j][i].M1/pG->U[k][j][i].d)/pG->dx1;
    
    /* Inner x2-boundary */
    if (j==pG->js && pG->lx2_id==-1)
        rmf -= tmp*(pG->U[k][j][i].M2/pG->U[k][j][i].d)/pG->dx2;
    
    /* Outer x2-boundary */
    if (j==pG->je && pG->rx2_id==-1)
        rmf += tmp*(pG->U[k][j][i].M2/pG->U[k][j][i].d)/pG->dx2;
    
    /* Inner x3-boundary */
    if (k==pG->ks && pG->lx3_id==-1)
        rmf -= tmp*(pG->U[k][j][i].M3/pG->U[k][j][i].d)/pG->dx3;
    
    /* Outer x3-boundary */
    if (k==pG->ke && pG->rx3_id==-1)
        rmf += tmp*(pG->U[k][j][i].M3/pG->U[k][j][i].d)/pG->dx3;
    
    return rmf;
}

/* Dump outward free energy flux */
/* NOTE:  This is actually an energy density flux, since history dumps are
 *   volume integrals */
static Real hst_Efree_out(const GridS *pG,const int i,const int j,const int k)
{
    Real x1,x2,x3,tmp,ef=0.0;
    
    /* Free Energy = Kinetic + Thermal + Gravitational */
    cc_pos(pG,i,j,k,&x1,&x2,&x3);
    tmp = 0.5*(SQR(pG->U[k][j][i].M1) +
               SQR(pG->U[k][j][i].M2) +
               SQR(pG->U[k][j][i].M3))/SQR(pG->U[k][j][i].d) +
    1.5*Iso_csound2 + pG->Phi[k][j][i];
    
    /* Inner x1-boundary */
    if (i==pG->is && pG->lx1_id==-1)
        ef -= tmp*pG->U[k][j][i].M1/pG->dx1;
    
    /* Outer x1-boundary */
    if (i==pG->ie && pG->rx1_id==-1)
        ef += tmp*pG->U[k][j][i].M1/pG->dx1;
    
    /* Inner x2-boundary */
    if (j==pG->js && pG->lx2_id==-1)
        ef -= tmp*pG->U[k][j][i].M2/pG->dx2;
    
    /* Outer x2-boundary */
    if (j==pG->je && pG->rx2_id==-1)
        ef += tmp*pG->U[k][j][i].M2/pG->dx2;
    
    /* Inner x3-boundary */
    if (k==pG->ks && pG->lx3_id==-1)
        ef -= tmp*pG->U[k][j][i].M3/pG->dx3;
    
    /* Outer x3-boundary */
    if (k==pG->ke && pG->rx3_id==-1)
        ef += tmp*pG->U[k][j][i].M3/pG->dx3;
    
    return ef;
}

/*------------------------------------------------------------------------------
 *  Function usr_*
 *
 *  User-defined output expressions.  Surface densities are produced by using
 *  the ':' reduction operator in the input file, which produces a global
 *  average along a given dimension.  For example, to obtain Sigma1, multiplying
 *  by Lx and averaging is the same as multiplying the sum over all i of
 *  rho(i,j,k) by dx=Lx/Nx.
 */

/* Dump kinetic energy in perturbations */
static Real dEk(const GridS *pG, const int i, const int j, const int k)
{ /* The kinetic energy in perturbations is 0.5*d*V^2 */
    return 0.5*(pG->U[k][j][i].M1*pG->U[k][j][i].M1 +
                pG->U[k][j][i].M2*pG->U[k][j][i].M2 +
                pG->U[k][j][i].M3*pG->U[k][j][i].M3)/pG->U[k][j][i].d;
}

static Real log10d(const GridS *pG,const int i,const int j,const int k)
{
    return log10(pG->U[k][j][i].d);
}

static Real usr_Sigma1(const GridS *pG,const int i,const int j,const int k)
{
    return pG->U[k][j][i].d*Lx;
}

static Real usr_Sigma2(const GridS *pG,const int i,const int j,const int k)
{
    return pG->U[k][j][i].d*Ly;
}

static Real usr_Sigma3(const GridS *pG,const int i,const int j,const int k)
{
    return pG->U[k][j][i].d*Lz;
}

/*------------------------------------------------------------------------------
 *  TURBULENCE FUNCTIONS
 *----------------------------------------------------------------------------*/

/*------------------------------------------------------------------------------
 *  Function pspect -- computes component of velocity with specific
 *  power spectrum in Fourier space
 *
 *  Velocity power spectrum returned in ampl
 *    klow   = multiple of 2 pi/L for cut-off at low  wavenumbers
 *    khigh  = multiple of 2 pi/L for cut-off at high wavenumbers
 *    expo   = exponent of power law
 *  Note that the fourier amplitudes are stored in an array with no
 *  ghost zones
 */
static void pspect(ath_fft_data *ampl)
{
    int i,j,k;
    Real q1,q2,q3,q4,q5;
    
    for (k=0; k<gnx3; k++) {
        for (j=0; j<gnx2; j++) {
            for (i=0; i<gnx1; i++) {
                /* Calculate the wave vector magnitude, k, in global coordinates,
                 * even if (i-gis,j-gjs,k-gks) is not on the local Grid. */
                q4 = KWVM(i-gis,j-gjs,k-gks);
                /* If k is within the cutoff range, generate 3 random numbers.
                 * This is done so that different processor numbers/topologies will
                 * produce the exact same random number sequences. */
                if ((q4 > klow) && (q4 < khigh)) {
                    q1 = ran2(&rseed);
                    q2 = ran2(&rseed);
                    q3 = ran2(&rseed);
                }
                
                /* If (i-gis,j-gjs,k-gks) is on the local Grid (in global coordinates),
                 * assign a wave amplitude (possibly 0). */
                if ((i>=gis) && (i<=gie) &&
                    (j>=gjs) && (j<=gje) &&
                    (k>=gks) && (k<=gke)) {
                    /* If k is within the cutoff range, compute an amplitude. */
                    if ((q4 > klow) && (q4 < khigh)) {
                        q5 = sqrt(-2.0*log(q1+1.0e-20))*cos(2.0*PI*q2);
                        ampl[OFST(i-gis,j-gjs,k-gks)][0] = q5*cos(2.0*PI*q3);
                        ampl[OFST(i-gis,j-gjs,k-gks)][1] = q5*sin(2.0*PI*q3);
                        
                        // Set power spectrum
                        q4 *= dkx; /* Multiply by 2 pi/L */

                        /* Decreasing power law */
                        ampl[OFST(i-gis,j-gjs,k-gks)][0] /= pow(q4,(expo+2.0)/2.0);
                        ampl[OFST(i-gis,j-gjs,k-gks)][1] /= pow(q4,(expo+2.0)/2.0);
                    } else {
                        /* Otherwise, introduce cut-offs at klow and khigh */
                        ampl[OFST(i-gis,j-gjs,k-gks)][0] = 0.0;
                        ampl[OFST(i-gis,j-gjs,k-gks)][1] = 0.0;
                    }
                }
            }
        }
    }
    if (gis==0 && gie==0 && gjs==0 && gje==0 && gks==0 && gke==0) {
        ampl[0][0] = 0.0;
        ampl[0][1] = 0.0;
    }
    
    return;
}

/*------------------------------------------------------------------------------
 *  Function project
 *
 *  Makes velocity perturbations divergence free
 */
static void project()
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

/*------------------------------------------------------------------------------
 *  Function transform
 *
 *  Generate velocities from fourier transform
 */
static inline void transform()
{
    /* Transform velocities from k space to physical space */
    ath_3d_fft(plan, fv1);
    ath_3d_fft(plan, fv2);
    ath_3d_fft(plan, fv3);
    
    /* Should technically renormalize (divide by gnx1*gnx2*gnx3) here, but
     * since we're going to renormalize to get the desired energy injection
     * rate anyway, there's no point */
    
    return;
}

/*------------------------------------------------------------------------------
 *  Function generate
 *
 *  Generate the velocity perturbations
 */
static inline void generate()
{
    /* Generate new perturbations following appropriate power spectrum */
    pspect(fv1);
    pspect(fv2);
    pspect(fv3);
    
    /* Transform perturbations to real space, but don't normalize until
     * just before we apply them in perturb() */
    transform();
    
    return;
}

/*------------------------------------------------------------------------------
 *  Function perturb
 *
 *  Shifts velocities so no net momentum change, normalizes to keep
 *  dedt fixed, and then sets velocities (momenta)
 */
static void perturb(GridS *pG)
{
    int i,j,k;
    int ind, mpierr;
    Real dvol, s, de, qa, v1, v2, v3;
    Real t0,t1,t2,t3,aa,bb,cc;
    Real m[4], gm[4];
    
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
    
    /* Calculate net momentum pertubation components t1, t2, t3 */
    t0 = 0.0;  t1 = 0.0;  t2 = 0.0;  t3 = 0.0;
    for (k=ks; k<=ke; k++) {
        for (j=js; j<=je; j++) {
            for (i=is; i<=ie; i++) {
                t0 += pG->U[k][j][i].d;
                
                /* The net momentum perturbation */
                t1 += pG->U[k][j][i].d * dv1[k][j][i];
                t2 += pG->U[k][j][i].d * dv2[k][j][i];
                t3 += pG->U[k][j][i].d * dv3[k][j][i];
            }
        }
    }
    
#ifdef MPI_PARALLEL
    /* Sum the perturbations over all processors */
    m[0] = t0;  m[1] = t1;  m[2] = t2;  m[3] = t3;
    mpierr = MPI_Allreduce(m, gm, 4, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    if (mpierr) ath_error("[radpargrav]: MPI_Allreduce error = %d\n", mpierr);
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
    
    /* Calculate unscaled energy of perturbations */
    aa = 0.0;  bb = 0.0;  cc = 0.0;
    for (k=ks; k<=ke; k++) {
        for (j=js; j<=je; j++) {
            for (i=is; i<=ie; i++) {
                /* Calculate velocity pertubation at cell center from
                 * perturbations at cell faces */
                v1 = dv1[k][j][i];
                v2 = dv2[k][j][i];
                v3 = dv3[k][j][i];
                
                aa += 0.5*(pG->U[k][j][i].d)*(SQR(v1) + SQR(v2) + SQR(v3));
                bb += (pG->U[k][j][i].M1)*v1 + (pG->U[k][j][i].M2)*v2 +
                (pG->U[k][j][i].M3)*v3;
                cc += 0.5*(SQR(pG->U[k][j][i].M1) + SQR(pG->U[k][j][i].M2) +
                           SQR(pG->U[k][j][i].M3))/pG->U[k][j][i].d;
            }
        }
    }
    
#ifdef MPI_PARALLEL
    /* Sum the perturbations over all processors */
    m[0] = aa;  m[1] = bb;  m[2] = cc;
    mpierr = MPI_Allreduce(m, gm, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    if (mpierr) ath_error("[radpargrav]: MPI_Allreduce error = %d\n", mpierr);
    aa = gm[0];  bb = gm[1];  cc = gm[2];
#endif /* MPI_PARALLEL */
    
    /* Rescale to give the correct energy injection rate */
    dvol = pG->dx1*pG->dx2*pG->dx3;
    aa = MAX(aa,1.0e-20);
    /* decaying turbulence (all in one shot) */
    /* NOTE:  In this case, dedt is really the desired total kinetic energy */
    
    cc -= dedt/dvol;
    s = (-bb + sqrt(SQR(bb) - 4.0*aa*cc))/(2.0*aa);
    if (isnan(s)) ath_error("[radpargrav]: s is NaN!\n");
    
    /* Apply momentum pertubations */
    t1 = 0.0;
    for (k=ks; k<=ke; k++) {
        for (j=js; j<=je; j++) {
            for (i=is; i<=ie; i++) {
                qa = s*pG->U[k][j][i].d;
                pG->U[k][j][i].M1 += qa*dv1[k][j][i];
                pG->U[k][j][i].M2 += qa*dv2[k][j][i];
                pG->U[k][j][i].M3 += qa*dv3[k][j][i];
            }
        }
    }
    
    return;
}

/*  */
void tracer_stars(const GridS *pG)
{
#ifdef STARPARTICLE
    StarParListS *pList=NULL;
    StarParS *pStar=NULL;
    
    int ip,jp,kp;
    Real x1,x2,x3;
    
    pList = pG->Lstars;
    while (pList) {
        pStar = &(pList->starpar);
        cc_ijk(pG,pStar->x1,pStar->x2,pStar->x3,&ip,&jp,&kp);
        /* If particle in this cell */
        if (i==ip && j==jp && k==kp) {
            return pStar->m/(pG->dx1*pG->dx2*pG->dx3);
        }
        pList = pList->next;
    }
#endif
    return 0.0;
}

/*------------------------------------------------------------------------------
 *  Function initialize
 *
 *  Allocate memory and initialize FFT plans
 *----------------------------------------------------------------------------*/
static void initialize(DomainS *pD)
{
    GridS *pG = pD->Grid;
    int i,j,k;
    int nbuf, mpierr;
    Real kwv, kpara, kperp;
    char donedrive = 0;
    Real pc_cgs=3.085678e18,kms_cgs=1.0e5,mH_cgs=1.6733e-24;
    
    /*----------------------------------------------------------------------------
     * Variables within this block are stored globally, and used
     * within preprocessor macros.  Don't create variables with
     * these names within your function if you are going to use
     * OFST(), KCOMP(), or KWVM() within the function! */
    
    /* Get local grid size */
    nx1 = pG->Nx[0];
    nx2 = pG->Nx[1];
    nx3 = pG->Nx[2];
    
    /* Get extents of local grid */
    is = pG->is;  ie = pG->ie;
    js = pG->js;  je = pG->je;
    ks = pG->ks;  ke = pG->ke;
    
    /* Get global grid size */
    gnx1 = pD->Nx[0];
    gnx2 = pD->Nx[1];
    gnx3 = pD->Nx[2];
    
    /* Get extents of local FFT grid in global coordinates */
    gis = pG->Disp[0];  gie = pG->Disp[0]+pG->Nx[0]-1;
    gjs = pG->Disp[1];  gje = pG->Disp[1]+pG->Nx[1]-1;
    gks = pG->Disp[2];  gke = pG->Disp[2]+pG->Nx[2]-1;
    
    /* Get size of arrays with ghost cells */
    il = is-nghost*(pG->Nx[0]>1);  iu = ie+nghost*(pG->Nx[0]>1);  nx1gh = iu-il+1;
    jl = js-nghost*(pG->Nx[1]>1);  ju = je+nghost*(pG->Nx[1]>1);  nx2gh = ju-jl+1;
    kl = ks-nghost*(pG->Nx[2]>1);  ku = ke+nghost*(pG->Nx[2]>1);  nx3gh = ku-kl+1;
    
    /* Get spatial extents of local grid */
    cc_pos(pG,il,jl,kl,&x1min,&x2min,&x3min);
    cc_pos(pG,iu,ju,ku,&x1max,&x2max,&x3max);
    
    /* Get local grid length, volume */
    dx = MAX(MAX(pG->dx1,pG->dx2),pG->dx3);
    dV = pG->dx1*pG->dx2*pG->dx3;
    
    //  printf("proc %d:  nx1=%d, nx2=%d, nx3=%d\n", myID_Comm_world,nx1,nx2,nx3);
    //  printf("proc %d:  gis=%d, gjs=%d, gks=%d\n", myID_Comm_world,gis,gjs,gks);
    //  ath_pout(0,"gnx1=%d, gnx2=%d, gnx3=%d\n", gnx1,gnx2,gnx3);
    //  ath_pout(0,"nx1gh=%d, nx2gh=%d, nx3gh=%d\n", nx1gh,nx2gh,nx3gh);
    
    /*--------------------------------------------------------------------------*/
    /* Get input parameters */
    
    /* Parse input file */
    rho_small   = par_getd("problem", "rho_small");
    rho         = par_getd("problem", "rho");
    rcloud      = par_getd("problem", "rcloud");
#ifdef MCTRACERS
    t_dep       = par_getd("problem", "t_dep");
    d_prop      = par_getd("problem", "d_prop");
#endif
    
    /* Initialize units
     *
     * This assumes inputs are converted such that code units are:
     * length unit = pc,
     * velocity unit = km/s
     * density unit = 1.4 m_H * 1 /cm^3  [ i.e. # density of H is in cm^-3]
     *
     * Thus:
     *
     * code length unit = 3.0856 e18 cm = 1 pc
     * code time unit   = pc/(km s^-1)=3.0856 e13 s = 0.978 Myr
     * code mass unit   = 1.4 mH * (pc/cm)^3 = 6.88 e31 g = 0.035 Msun
     */
    UnitS units;
    
    /* Code length unit */
    units.Lcode = pc_cgs;
    ath_pout(0,"length unit [cm]:  %e\n",units.Lcode);
    
    /* Code mass unit */
    units.Mcode = 1.4*mH_cgs*CUBE(pc_cgs);
    ath_pout(0,"mass unit [g]:  %e\n",units.Mcode);
    
    /* Code velocity unit */
    units.Vcode = kms_cgs;
    ath_pout(0,"velocity unit [cm s^-1]:  %e\n",units.Vcode);
    
    /* Code time unit */
    units.Tcode = units.Lcode/units.Vcode;
    ath_pout(0,"time unit [s]:  %e\n",units.Tcode);
    
    init_units(&units);
    ath_pout(1,"units.cm = %e\n",units.cm);
    ath_pout(1,"units.g  = %e\n",units.g);
    ath_pout(1,"units.s  = %e\n",units.s);
    
    Lx = (pD->MaxX[0] - pD->MinX[0])*pc_cgs*units.cm;
    Ly = (pD->MaxX[1] - pD->MinX[1])*pc_cgs*units.cm;
    Lz = (pD->MaxX[2] - pD->MinX[2])*pc_cgs*units.cm;
    
    /* Code density unit */
    ath_pout(0,"density unit [g cm^-3]:  %e\n",units.Dcode);
    
    /* Code energy unit */
    ath_pout(0,"1 erg in code units = %e\n",units.erg);
    
    /* Code temperature unit */
    /* Chosen so that aR=1 in code units */
    ath_pout(0,"1 K in code units = %e\n",units.K);
    
    /* Scale rho to code units */
    rho *= 1.4*mH_cgs; // convert to g/cm^3
    rho *= units.g/CUBE(units.cm);
    ath_pout(0,"rho in code units = %e\n",rho);
    
    /* Scale rcloud to code units */
    rcloud *= units.pc;
    ath_pout(0,"rcloud in code units = %e\n",rcloud);
    
    /* scale rho_small by density maximum */
    rho_small *= rho;
    ath_pout(0,"rho_small in code units = %e\n",rho_small);
    
    /* MGMC */
    M_GMC = rho*CUBE(rcloud);
    ath_pout(0,"M_GMC in code units = %e\n", M_GMC);
    
#ifdef SELF_GRAVITY
    /* Set gravity constant */
    /* G = G_cgs [cm^3 g^-1 s^-2] */
    four_pi_G = 4.0*PI*units.G;
    grav_mean_rho = rho;
    ath_pout(0,"4*pi*G in code units = %e\n",four_pi_G);
    
    ath_pout(0,"rho_Truelove in code units = %e\n",
             SQR(PI)*Iso_csound2/(4.0*four_pi_G*SQR(dx)));
    ath_pout(0,"rho_LP in code units = %e\n",
             8.86*Iso_csound2*4.0/(four_pi_G*SQR(dx)));
    
    /* Ensure the Jeans length is resolved by at least 4 zones */
    L_Jeans = Iso_csound*2.0*PI/sqrt(four_pi_G*rho);
    ath_pout(0,"L_Jeans in code units ~ %e, resolved over %1.1f zones\n",
             L_Jeans,L_Jeans/dx);
    if (L_Jeans/dx < 4.0)
        ath_error("[radpargrav]:  The Jeans length is not resolved!\n");
    ath_pout(0,"t_ff in code units = %e\n",sqrt(3.0*PI*PI/(8.0*four_pi_G*rho)));
#endif
    
    /*--------------------------------------------------------------------------*/
    /* Set up perturbation parameters */
    
    /* interval for generating new driving spectrum; also interval for
     * driving when IMPULSIVE_DRIVING is used */
#ifdef MHD
    /* magnetic field strength */
    beta = par_getd("problem","beta");
    /* beta = isothermal pressure/magnetic pressure */
    B0 = sqrt(2.0*Iso_csound2*rhobar/beta);
#endif /* MHD */
    /* energy injection rate */
    //  dedt = par_getd("problem","dedt");
    //  dedt *= units.erg;
    /* Set dedt(=Ekin,tot) such that Ekin = 0.5*alpha_vir*Egrav, where
     * Egrav = (3/5)*G*M^2/R  */
    Egrav = 0.6*(four_pi_G/(4.0*PI))*SQR(M_GMC)/rcloud;
    dedt = 0.5*par_getd("problem","alpha_vir")*Egrav;
    ath_pout(0,"Ekin in code units = %e\n",dedt);
    ath_pout(0,"Egrav in code units = %e\n",Egrav);
    v_turb = sqrt(2.0*dedt/M_GMC);
    ath_pout(0,"v_turb in code units = %e\n",v_turb);
    
    /* parameters for spectrum */
    expo = par_getd("problem","expo");
        
    /* Cutoff wavenumbers of spectrum */
    klow = par_getd("problem","klow"); /* in integer units */
    khigh = par_getd("problem","khigh"); /* in integer units */
    dkx = 2.0*PI/MAX(MAX(pG->dx1*gnx1,pG->dx2*gnx2),pG->dx3*gnx3); /* convert k from integer */
    
    /* If restarting with decaying turbulence, no driving necessary. */
    if (pG->time > 0) {
        donedrive = 1;
    }
    
    if (donedrive == 0) {
        /* Allocate memory for components of velocity perturbation */
        if ((dv1=(Real***)calloc_3d_array(nx3gh,nx2gh,nx1gh,sizeof(Real)))==NULL) {
            ath_error("[radpargrav]: Error allocating memory for vel pert\n");
        }
        if ((dv2=(Real***)calloc_3d_array(nx3gh,nx2gh,nx1gh,sizeof(Real)))==NULL) {
            ath_error("[radpargrav]: Error allocating memory for vel pert\n");
        }
        if ((dv3=(Real***)calloc_3d_array(nx3gh,nx2gh,nx1gh,sizeof(Real)))==NULL) {
            ath_error("[radpargrav]: Error allocating memory for vel pert\n");
        }
    }
    
    /* Initialize the FFT plan */
    plan = ath_3d_fft_quick_plan(pD, NULL, ATH_FFT_BACKWARD);
    
    /* Allocate memory for FFTs */
    if (donedrive == 0) {
        fv1 = ath_3d_fft_malloc(plan);
        fv2 = ath_3d_fft_malloc(plan);
        fv3 = ath_3d_fft_malloc(plan);
    }
    
    /* Enroll outputs */
    dump_history_enroll(hst_dEk,"<hst_dEk>");
    //dump_history_enroll(hst_mass_stars,"<mass_stars>");
    dump_history_enroll(hst_Egrav_tot,"<dEgrav_tot>");
    dump_history_enroll(hst_Egrav_gas,"<dEgrav_gas>");
    //dump_history_enroll(hst_Egrav_stars,"<dEgrav_stars>");
    dump_history_enroll(hst_Wgrav_gas,"<dWgrav_gas>");
    //dump_history_enroll(hst_Wgrav_stars,"<dWgrav_stars>");
    dump_history_enroll(hst_Mr_out,"<Mr_out>");
    dump_history_enroll(hst_Efree_out,"<Efree_out>");
#ifdef STARPARTICLES
    dump_history_enroll(hst_Er,"\\int Er dV");
#endif
    
    return;
}


/*------------------------------------------------------------------------------
 *  Function ran2
 *
 *  The routine ran2() is extracted from the Numerical Recipes in C
 *  (version 2) code.  I've modified it to use doubles instead of
 *  floats. -- T. A. Gardiner -- Aug. 12, 2003
 *
 *  Long period (> 2 x 10^{18}) random number generator of L'Ecuyer
 *  with Bays-Durham shuffle and added safeguards.  Returns a uniform
 *  random deviate between 0.0 and 1.0 (exclusive of the endpoint
 *  values).  Call with idum = a negative integer to initialize;
 *  thereafter, do not alter idum between successive deviates in a
 *  sequence.  RNMX should appriximate the largest floating point 
 *  value that is less than 1.
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
#define NDIV (1+IMM1/NTAB)
#define RNMX (1.0-DBL_EPSILON)
#define NTAB 32

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

#undef OFST
#undef KCOMP
#undef KWVM
