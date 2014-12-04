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

    int i, j, n1z, n2z;
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
    
    int i, j, n1z, n2z;
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
            //          assert(tracer->newList != list);
            //          pnode->newList->currTail = pnode;
            tracer_list_remove(list, tracer);
            Tracerlist_add(tracer->newList, tracer);
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
        assert(pnode->newList != list);
        /* Remove from old list */
        tracer_list_remove(list, pnode);
        /* Move to end of new list */
        Tracerlist_add(pnode->newList, pnode);
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
#endif