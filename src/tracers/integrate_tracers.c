#include "../copyright.h"
/*===========================================================================*/
/* \file integrate_mctracers.c                                               */
/* \brief Initialize Monte Carlo tracer related structures and functions.    */
/*									     */
/* PURPOSE: 								     */
/*									     */
/* CONTAINS PUBLIC FUNCTIONS:						     */
/*============================================================================*/
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "../defs.h"
#include "../athena.h"
#include "../prototypes.h"
#include "prototypes.h"
#include "tracers.h"
#include "../globals.h"
#include <float.h>
#include <gsl/gsl_rng.h>
#include <time.h>

#if defined(MCTRACERS) || defined(VFTRACERS)  /* endif at the end of the file */

/*----------------------------------------------------------------------------*/
/* Calculates pseudo random number based on seed */
/*----------------------------------------------------------------------------*/
static long int seed;
/* Current random array index */
static gsl_rng *rng;

/* ========================================================================== */

/*----------------------------------------------------------------------------*/
/* Fill pG->random with ran2() generated doubles                              */
/*----------------------------------------------------------------------------*/
void ran_gen(GridS *pG) {

    long s = time(NULL);
    int pid = myID_Comm_world;
    
    seed = abs(((s*181)*((pid-83)*359))%104729);
    rng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rng, seed);
    
//    n1z = pG->Nx[0] + 2*nghost;
//    n2z = pG->Nx[1] + 2*nghost;
//    ran2(&seed);
//    pG->rand_list = (double **)calloc_2d_array(n1z, n2z, sizeof(double));
    
    return;
}

/*----------------------------------------------------------------------------*/
/* Fill pG->random with ran2() generated doubles                              */
/*----------------------------------------------------------------------------*/
void ran_gen_list(GridS *pG) {
    
//    seed = 0;
//    rng = gsl_rng_alloc(gsl_rng_mt19937);
//    gsl_rng_set(rng, seed);
//    
//    for (j = pG->js-nghost; j<=pG->je+nghost; j++) {
//        for (i = pG->is-nghost; i<=pG->ie+nghost; i++) {
//            pG->rand_list[i][j] = gsl_rng_uniform(rng);
//        }
//    }
    
    return;
}

/*----------------------------------------------------------------------------*/
/* Sweeps through list and moves nodes flagged for removal                    */
/*----------------------------------------------------------------------------*/
void Tracerlist_sweep(TracerListS *list, GridS *pG)
{
	TracerS *tracer = list->Head;
    TracerS *next;
    
	while (tracer) {
        /* If have reached nodes added during current cycle, return */
        //        if (list->currTail) {
        //            if (pnode == list->currTail) return;
        //        }
        /* If pnode flagged, move pointers around in current list */
        next = tracer->Next;
        if (tracer->newList) {
            //          pnode->newList->currTail = pnode;
            if (tracer->newList != list) {
                tracer_list_remove(list, tracer);
                Tracerlist_add(tracer->newList, tracer);
            }
            else {
                tracer->newList = NULL;
            }
        }
        tracer = next;
	}
    return;
}

/*----------------------------------------------------------------------------*/
/* Sweeps through list and moves nodes flagged for removal at boundary        */
/*----------------------------------------------------------------------------*/
void Tracerlist_sweep_bc(TracerListS *list)
{
	TracerS *pnode = list->Head;
    TracerS *next;
    
	while (pnode) {
        next = pnode->Next;
#ifdef DEBUG
        assert(pnode->newList != list);
#endif /* DEBUG */
        /* Remove from old list */
        tracer_list_remove(list, pnode);
        /* Move to end of new list */
        Tracerlist_add(pnode->newList, pnode);
        pnode->newList = NULL;
        /* Iterate through list */
		pnode = next;
	}
    
    return;
}

/*----------------------------------------------------------------------------*/
/* Iterate through list and determine tracers to be moved */
/*----------------------------------------------------------------------------*/
void prob_iterate_x1(TracerListS *list, double pflux, GridS *pG)
{
    TracerS *pnode;
    pnode = list->Head;
    double rand;
    int k = list->k;
    int j = list->j;
    int i = list->i;
    
    while(pnode) {
        /* If pnode has not been flagged earlier in current loop */
        if (pnode->newList == NULL) {
            rand = gsl_rng_uniform(rng);
            
            if (pflux > 0) {
                if (rand <= pflux) {
                    pnode->newList = &((pG->GridLists)[k][j][i+1]);
                }
            }
            if (pflux < 0) {
                if (rand <= -pflux) {
                    pnode->newList = &((pG->GridLists)[k][j][i-1]);
                }
            }
        }
        pnode = pnode->Next;
    }
}

/*----------------------------------------------------------------------------*/
/* Iterate through list and determine tracers to be moved */
/*----------------------------------------------------------------------------*/
void prob_iterate_x2(TracerListS *list, double pflux, GridS *pG)
{
    TracerS *pnode;
    pnode = list->Head;
    double rand;
    int k = list->k;
    int j = list->j;
    int i = list->i;
    
    while(pnode) {
        /* If pnode has not been flagged earlier in current loop */
        if (pnode->newList == NULL) {
            rand = gsl_rng_uniform(rng);
//            rand = rand_list[curr_index];
//            curr_index++;
//            if (curr_index == rand_len) ran_gen_list(pG);
            if (pflux > 0) {
                if (rand <= pflux) {
                    pnode->newList = &((pG->GridLists)[k][j+1][i]);
                }
            }
            if (pflux < 0) {
                if (rand <= -pflux) {
                    pnode->newList = &((pG->GridLists)[k][j-1][i]);
                }
            }
        }
        pnode = pnode->Next;
    }
}

/*----------------------------------------------------------------------------*/
/* Iterate through list and determine tracers to be moved */
/*----------------------------------------------------------------------------*/

void prob_iterate_x3(TracerListS *list, double pflux, GridS *pG)
{
    TracerS *pnode;
    pnode = list->Head;
    double rand;
    int i = list->i;
    int j = list->j;
    int k = list->k;
    
    while(pnode) {
        /* If pnode has not been flagged earlier in current loop */
        if (pnode->newList == NULL) {
            rand = gsl_rng_uniform(rng);
            if (pflux > 0) {
                if (rand <= pflux) {
                    pnode->newList = &((pG->GridLists)[k+1][j][i]);
                }
            }
            if (pflux < 0) {
                if (rand <= -pflux) {
                    pnode->newList = &((pG->GridLists)[k-1][j][i]);
                }
            }
        }
        pnode = pnode->Next;
    }
}

/*----------------------------------------------------------------------------*/
/* Mark tracers within control volume when star particle created   */
/* ic, jc, kc are the indices of the new star particle */
/*----------------------------------------------------------------------------*/
#ifdef STAR_PARTICLE
void flag_tracer_grid(GridS *pG)
{
    Real3Vect cell1;
    Real a, b, c;
    int ic, jc, kc;
    StarParListS *pLstars = NULL;
    StarParS *pStar = NULL;
    
    /* cell1 is a shortcut expressions as well as dimension indicator */
    if (pG->Nx[0] > 1)  cell1.x1 = 1.0/pG->dx1;  else cell1.x1 = 0.0;
    if (pG->Nx[1] > 1)  cell1.x2 = 1.0/pG->dx2;  else cell1.x2 = 0.0;
    if (pG->Nx[2] > 1)  cell1.x3 = 1.0/pG->dx3;  else cell1.x3 = 0.0;
    
    pLstars = pG->Lstars;
    while (pLstars) {
        pStar = &(pLstars->starpar);
        celli(pG, pStar->x1, cell1.x1, &ic, &a);
        cellj(pG, pStar->x2, cell1.x2, &jc, &b);
        cellk(pG, pStar->x3, cell1.x3, &kc, &c);
        flag_tracer_star(pG, ic, jc, kc, pStar->id);
        pLstars = pLstars->next;
    }
}


/*----------------------------------------------------------------------------*/
/* Mark tracers within control volume when star particle created   */
/* ic, jc, kc are the indices of the new star particle */
/*----------------------------------------------------------------------------*/
void flag_tracer_star(GridS *pG, int ic, int jc, int kc, int star_id)
{
    int i, j, k;
    TracerListS *list;
    TracerS *tracer;
    
    for (k = kc-NSINK_STARP; k <= kc+NSINK_STARP; k++) {
        for (j = jc-NSINK_STARP; j <= jc+NSINK_STARP; j++) {
            for (i = ic-NSINK_STARP; i <= ic+NSINK_STARP; i++) {
                list = &((pG->GridLists)[k][j][i]);
                tracer = list->Head;
                while(tracer) {
                    tracer->prop->star_id = star_id;
                    tracer = tracer->Next;
                }
            }
        }
    }
}

/*----------------------------------------------------------------------------*/
/* Output number of tracers in each threshold category within starpar  */
/*----------------------------------------------------------------------------*/
void output_tracer_star(GridS *pG, double thresh)
{
    TracerListS *list;
    TracerS *tracer;
    double thresh1 = thresh;
    double thresh2 = thresh*10;
    double thresh3 = thresh*10;
    double thresh4 = thresh*100;
    double thresh5 = thresh*100;
    int thresh1_cnt = 0;
    int thresh2_cnt = 0;
    int thresh3_cnt = 0;
    int thresh4_cnt = 0;
    int thresh5_cnt = 0;
    int thresh1_star = 0;
    int thresh2_star = 0;
    int thresh3_star = 0;
    int thresh4_star = 0;
    int thresh5_star = 0;
    int is = pG->is;
    int ie = pG->ie;
    int js = pG->js;
    int je = pG->je;
    int ks = pG->ks;
    int ke = pG->ke;
    int i, j, k;
    FILE *hp;
    
    if((hp = fopen("tracer_star.txt","a")) == NULL) {
        ath_error("[output_tracer_star]: Unable to open dump file\n");
        return;
    }
    for (k=ks; k<=ke; k++) {
        for (j=js; j<=je; j++) {
            for (i=is; i<=ie; i++) {
                list = &((pG->GridLists)[k][j][i]);
                tracer = list->Head;
                while (tracer) {
                    if (tracer->prop->d_init > thresh1) {
                        thresh1_cnt++;
                        if (tracer->prop->star_id != -1) {
                            thresh1_star++;
                        }
                    }
                    if (tracer->prop->d_init > thresh2) {
                        thresh2_cnt++;
                        if (tracer->prop->star_id != -1) {
                            thresh2_star++;
                        }
                    }
                    if (tracer->prop->d_init > thresh3) {
                        thresh3_cnt++;
                        if (tracer->prop->star_id != -1) {
                            thresh3_star++;
                        }
                    }
                    if (tracer->prop->d_init > thresh4) {
                        thresh4_cnt++;
                        if (tracer->prop->star_id != -1) {
                            thresh4_star++;
                        }
                    }
                    if (tracer->prop->d_init > thresh5) {
                        thresh5_cnt++;
                        if (tracer->prop->star_id != -1) {
                            thresh5_star++;
                        }
                    }
                    tracer = tracer->Next;
                }
            }
        }
    }
    
    fprintf(hp,"%lf %d %d %d %d %d %d %d %d %d %d\n",
            pG->time, thresh1_cnt, thresh1_star, thresh2_cnt, thresh2_star,
            thresh3_cnt, thresh3_star, thresh4_cnt, thresh4_star, thresh5_cnt, thresh5_star);
    fclose(hp);
}

#endif /* STAR_PARTICLE */

/*----------------------------------------------------------------------------*/
/* Implement top hat algorithm */
/*----------------------------------------------------------------------------*/
#ifdef TOPHAT
void mc_tophat(GridS *pG)
{
    int i, ii, j, jj, k, kk;
    int il, ir, jl, jr, kl, kr;
    int is, ie, js, je, ks, ke;
    int tot, N, nx, ny, nz;
    int rfilter = 2;
    MCtophatS *list;
    
    is = pG->is;
    ie = pG->ie;
    js = pG->js;
    je = pG->je;
    ks = pG->ks;
    ke = pG->ke;
    
    nx = pG->Nx[0];
    ny = pG->Nx[1];
    nz = pG->Nx[2];
        
    for (k=ks; k<=ke; k++) {
    	for (j=js; j<=je; j++) {
            for (i=is; i<=ie; i++) {
                list = &((pG->TopHatGrid)[k][j][i]);
                il = i-rfilter;
                ir = i+rfilter+1;
                jl = j-rfilter;
                jr = j+rfilter+1;
                kl = k-rfilter;
                kr = k+rfilter+1;
                N = 0;
                tot = 0;
                for (kk=kl; kk<=kr; kk++) {
                    for (jj=jl; jj<=jr; jj++) {
                        for (ii=il; ii<=ir; ii++) {
                            if ((ii-i)*(ii-i) + (jj-j)*(jj-j) + (kk-k)*(kk-k) <= rfilter*rfilter) {
                                tot += pG->U[kk%nz][jj%ny][ii%nx].d;
                                N++;
                            }
                        }
                    }
                }
                list->count = tot/N;
            }
        }
    }
}
#endif /* TOPHAT */
#endif /* TRACERS */
