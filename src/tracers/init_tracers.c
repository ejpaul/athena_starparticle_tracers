#include "../copyright.h"
/*===========================================================================*/
/* \file init_tracers.c                                                    */
/* \brief Initialize Monte Carlo tracer related structures and functions.    */
/*									     */
/* PURPOSE: 								     */
/*									     */
/* CONTAINS PUBLIC FUNCTIONS:						     */
/*============================================================================*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "../defs.h"
#include "../athena.h"
#include "../prototypes.h"
#include "prototypes.h"
#include "tracers.h"
#include "../globals.h"
#include <float.h>
#include <gsl/gsl_rng.h>
#include <time.h>

#if defined(MCTRACERS) || defined(VFTRACERS)/* endif at the end of the file */

/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 *============================================================================*/

/* Current id */
static double current_tracer;
/* Current random array index */

/*=========================== PUBLIC FUNCTIONS ===============================*/
/*----------------------------------------------------------------------------*/
/* Initializes tracer grid */
/*----------------------------------------------------------------------------*/
void init_tracer_grid(MeshS *pM)
{
    current_tracer = 0;
    int n1z, n2z, n3z, i, j, k, il, ir, jl, jr, kl, kr;
    int ks, ke, js, je, is, ie;
    TracerListS *list;
    
    DomainS *pD;
    GridS   *pG;
    pD = (DomainS*)&(pM->Domain[0][0]);  /* set ptr to Domain */
    pG = pD->Grid;          /* set ptr to Grid */
    
    is = pG->is;
    ie = pG->ie;
    js = pG->js;
    je = pG->je;
    ks = pG->ks;
    ke = pG->ke;
    
    /* ---------  Allocate 3D arrays for MClistS based on grid size ------ */
    /* pG->GridLists same size as pG->U */
    
    if (pG->Nx[0] > 1) {
        n1z = pG->Nx[0] + 2*nghost;
        il = is - nghost;
        ir = ie + nghost;
    }
    else {
        n1z = 1;
        il = is;
        ir = ie;
    }
    
    if (pG->Nx[1] > 1) {
        n2z = pG->Nx[1] + 2*nghost;
        jl = js - nghost;
        jr = je + nghost;
    }
    else {
        n2z = 1;
        jl = js;
        jr = je;
    }
    
    if (pG->Nx[2] > 1) {
        n3z = pG->Nx[2] + 2*nghost;
        kl = ks - nghost;
        kr = ke + nghost;
    }
    else {
        n3z = 1;
        kl = ks;
        kr = ke;
    }
	
    /* Build a 3D array of type TracerListS */
    pG->GridLists = (TracerListS***)calloc_3d_array(n3z, n2z, n1z, sizeof(TracerListS));
    if (pG->GridLists == NULL) goto on_error;
    
#ifdef TOPHAT
    /* Build a 3D array of type MCtophatS */
    pG->TopHatGrid = (TracertophatS***)calloc_3d_array(n3z, n2z, n1z, sizeof(TracertophatS));
    if (pG->TopHatGrid == NULL) goto on_error;
#endif /* TOPHAT */
    
    /* Initialize PRNG */
    ran_gen(pG);
    
    /* Loop over grid, initialize TracerlistS */
    for (k=kl; k<=kr; k++) {
    	for (j=jl; j<=jr; j++) {
            for (i=il; i<=ir; i++) {
                list = &((pG->GridLists)[k][j][i]);
                /* Initialize count to 0, increment when adding */
                list->count = 0;
                /* Initialize reduced mass */
#ifdef MCTRACERS
                list->Rmass = pG->U[k][j][i].d;
                list->currTail = NULL;
#endif
                /* Initialize position of grid cell */
                list->i = i;
                list->j = j;
                list->k = k;
                list->Head = NULL;
                list->Tail = NULL;
            }
        }
    }
    return;
on_error:
    ath_error("[init_mctracers]: Error allocating memory.\n");
}

/*----------------------------------------------------------------------------*/
/* Initialize uniform density mc tracers */
/*----------------------------------------------------------------------------*/
void tracer_init_unif(GridS *pG) {
    int ks, ke, js, je, is, ie;
    int i, j, k;
    TracerListS *list;
    int n = par_geti("problem","N_tracers");
    is = pG->is;
    ie = pG->ie;
    js = pG->js;
    je = pG->je;
    ks = pG->ks;
    ke = pG->ke;
    
    /* Only initialize tracers in active cells */
    for (k=ks; k<=ke; k++) {
        for (j=js; j<=je; j++) {
            for (i=is; i<=ie; i++) {
                list = &((pG->GridLists)[k][j][i]);
                list->count = 0;
                /* Initialize position, pointers of list */
                list->i = i;
                list->j = j;
                list->k = k;
                list->Head = NULL;
                list->Tail = NULL;
#ifdef MCTRACERS
                list->currTail = NULL;
#endif /* MCTRACERS */
                init_tracer_list(pG, list, n);
            }
        }
    }
    return;
}

void tracer_init_threshold(GridS *pG, Real rho) {
    int ks, ke, js, je, is, ie;
    int i, j, k;
    TracerListS *list;
    Real d;
    int n = par_geti("problem","N_proportional");
    is = pG->is;
    ie = pG->ie;
    js = pG->js;
    je = pG->je;
    ks = pG->ks;
    ke = pG->ke;
    
    /* Only initialize tracers in active cells */
    for (k=ks; k<=ke; k++) {
        for (j=js; j<=je; j++) {
            for (i=is; i<=ie; i++) {
                list = &((pG->GridLists)[k][j][i]);
                list->count = 0;
                /* Initialize position, pointers of list */
                list->i = i;
                list->j = j;
                list->k = k;
                list->Head = NULL;
                list->Tail = NULL;
#ifdef MCTRACERS
                list->currTail = NULL;
#endif /* MCTRACERS */
                d = pG->U[k][j][i].d;
                if (d > rho) {
                    init_tracer_list(pG, list, d*n);
                }
            }
        }
    }
    return;
}

/*----------------------------------------------------------------------------*/
/* Initialize density mc tracers proportional to fluid density */
/*----------------------------------------------------------------------------*/
void tracer_init_proportional(GridS *pG) {
    int ks, ke, js, je, is, ie;
    int i, j, k;
    TracerListS *list;
    Real d;
    int m;
    int n = par_geti("problem","N_proportional");
    is = pG->is;
    ie = pG->ie;
    js = pG->js;
    je = pG->je;
    ks = pG->ks;
    ke = pG->ke;
    
    /* Only initialize tracers in active cells */
    for (k=ks; k<=ke; k++) {
        for (j=js; j<=je; j++) {
            for (i=is; i<=ie; i++) {
                list = &((pG->GridLists)[k][j][i]);
                list->count = 0;
                /* Initialize position, pointers of list */
                list->i = i;
                list->j = j;
                list->k = k;
                list->Head = NULL;
                list->Tail = NULL;
                d = pG->U[k][j][i].d;
                m = d*n;
#ifdef MCTRACERS
                list->currTail = NULL;
#endif /* MCTRACERS */
                init_tracer_list(pG, list, m);
            }
        }
    }
#ifdef DEBUG
    tracer_debug(pG);
#endif /* DEBUG */
    return;
}

/*----------------------------------------------------------------------------*/
/* Initialize mc tracers in x_l ghost zone (1 layer) for outflow problem      */
/* e.g. shk_cloud                                                             */
/*----------------------------------------------------------------------------*/
void tracer_init_xlinflow(GridS *pG) {
    int ks, ke, js, je, is, ie;
    int i, j, k;
    int n = par_geti("problem","N_proportional");
    TracerListS *list;
    Real d;
    is = pG->is;
    ie = pG->ie;
    js = pG->js;
    je = pG->je;
    ks = pG->ks;
    ke = pG->ke;
    
    /* Initilize tracers in l ghotst zone */
    for (k=ks; k<=ke; k++) {
        for (j=js; j<=je; j++) {
            i = is-1;
                list = &((pG->GridLists)[k][j][i]);
                list->count = 0;
                /* Initialize position, pointers of list */
                list->i = i;
                list->j = j;
                list->k = k;
                list->Head = NULL;
                list->Tail = NULL;
                d = pG->U[k][j][i].d;
#ifdef MCTRACERS
            list->currTail = NULL;
#endif /* MCTRACERS */
            init_tracer_list(pG, list, n*d);
        }
    }
    return;
}

/*----------------------------------------------------------------------------*/
/* Debugs MCLIST through assert statements */
/*----------------------------------------------------------------------------*/
void tracer_debug(GridS *pG) {
    int i, j, k, count, tot;
    int is = pG->is;
    int ie = pG->ie;
    int js = pG->js;
    int je = pG->je;
    int ks = pG->ks;
    int ke = pG->ke;
    TracerListS *list;
    TracerS *pnode;
    
#ifdef VFTRACERS 
    int i1, i2, i3;
    Real a, b, c;
    Real3Vect cell1;              /* one over dx1, dx2, dx3 */
    /* cell1 is a shortcut expressions as well as dimension indicator */
    if (pG->Nx[0] > 1)  cell1.x1 = 1.0/pG->dx1;  else cell1.x1 = 0.0;
    if (pG->Nx[1] > 1)  cell1.x2 = 1.0/pG->dx2;  else cell1.x2 = 0.0;
    if (pG->Nx[2] > 1)  cell1.x3 = 1.0/pG->dx3;  else cell1.x3 = 0.0;
#endif /* VFTRACERS */
    
    tot = 0;
    for (k=ks; k<=ke; k++) {
    	for (j=js; j<=je; j++) {
            for (i=is; i<=ie; i++) {
                list = &((pG->GridLists)[k][j][i]);
                count = 0;
                pnode = list->Head;
                if (pnode) {
                    assert(pnode->Prev == NULL);
                }
                if (list->count == 0) {
                    assert(list->Head == NULL);
                    assert(list->Tail == NULL);
                }
#ifdef MCTRACERS
                assert(list->Rmass >= 0);
#endif /* MCTRACERS */
                while(pnode) {
                    assert(pnode->newList==NULL);
                    if (pnode->Prev == NULL) {
                        assert(pnode == list->Head);
                    }
                    if (pnode->Next == NULL) {
                        assert(pnode == list->Tail);
                    }
                    if (pnode->Next) {
                        assert(pnode->Next->Prev == pnode);
                    }
                    if (pnode->Prev) {
                        assert(pnode->Prev->Next == pnode);
                    }
#ifdef VFTRACERS
// Check that tracer is in the right cell
                    celli(pG, pnode->x1, cell1.x1, &i1, &a);
                    assert(i == i1);
                    cellj(pG, pnode->x2, cell1.x2, &i2, &b);
                    assert(j == i2);
                    cellk(pG, pnode->x3, cell1.x3, &i3, &c);
                    assert(k == i3);
#endif /* VFTRACERS */
                    count++;
                    pnode = pnode->Next;
                }
                assert(count == list->count);
                tot += count;
            }
        }
    }
    printf("Total # tracers: %d\n", tot);
    fflush(0);
}

/*----------------------------------------------------------------------------*/
/* Calculates unique id number for mc tracer */
/* In the form of (subsequent integers).(myID_Comm_world) */
/*----------------------------------------------------------------------------*/

double get_tracer_id() {
    current_tracer++;
    int len;
    double ext;
    if (myID_Comm_world <= 0) {
        ext = 0;
    }
    else {
        len = floor(log10(abs(myID_Comm_world))) + 1;
        ext = (double)myID_Comm_world/((double)len*10.0);
    }
    return current_tracer + ext;
}

/*----------------------------------------------------------------------------*/
/* Initializes nodes of given list */
/*----------------------------------------------------------------------------*/
void init_tracer_list(GridS *pG, TracerListS *list, int n)
{
    Real x1, x2, x3;
    int i = list->i;
    int j = list->j;
    int k = list->k;
    TracerS *tracer = NULL;
    double rand;
    
    long s = time(NULL);
    int pid = myID_Comm_world;
    long seed = abs(((s*181)*((pid-83)*359))%104729);
    static gsl_rng *rng;
    rng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rng, seed);
    
#ifdef MCTRACERS
    list->Rmass = pG->U[k][j][i].d;
#endif /* MCTRACERS */
    list->Tail = NULL;
    list->Head = NULL;
    
    // Calculates cell-centered x,y,z given i,j,k
    cc_pos(pG, i, j, k, &x1, &x2, &x3);
	for (int i = 1; i <= n; i++) {
        /* Initialize tracer, add to tail of list */
		tracer = init_tracer();
        tracer->prop = (TracerPropS *)malloc(sizeof(TracerPropS));
        tracer->prop->id = get_tracer_id();
        tracer->prop->d_init = pG->U[k][j][i].d;
        tracer->prop->i_init = (int)list->i;
        tracer->prop->j_init = (int)list->j;
        tracer->prop->k_init = (int)list->k;
#ifdef STAR_PARTICLE
        tracer->prop->star_id = -1;
#endif /* STAR_PARTICLE */
#ifdef VFTRACERS 
        // Randomly distribute ntracers in coordinate space
        rand = gsl_rng_uniform(rng);
        tracer->x1 = x1 - 0.5*pG->dx1 + pG->dx1*rand;
        rand = gsl_rng_uniform(rng);
        tracer->x2 = x2 - 0.5*pG->dx2 + pG->dx2*rand;
        rand = gsl_rng_uniform(rng);
        tracer->x3 = x3 - 0.5*pG->dx3 + pG->dx3*rand;
#endif /* VFTRACERS */
		Tracerlist_add(list, tracer);
	}
	return;
}

/*----------------------------------------------------------------------------*/
/* Add tracer to end of list */
/*----------------------------------------------------------------------------*/
void Tracerlist_add(TracerListS *list, TracerS *tracer)
{
    /* If nonempty list add to end */
	if (list->Tail != NULL) {
        assert(list->Head);
        assert(list->count > 0);
        assert(list->Tail->Next == NULL);
        
		list->Tail->Next = tracer;
        tracer->Prev = list->Tail;
	}
    /* If empty list */
    else {
        assert(list->count == 0);
        assert(list->Tail == NULL);
        assert(list->Head == NULL);
        
        list->Head = tracer;
        tracer->Prev = NULL;
    }
    
    list->Tail = tracer;
    tracer->Next = NULL;
    tracer->newList = NULL;
	list->count++;
}

/*----------------------------------------------------------------------------*/
/* Returns initialized tracer for given timestep and id */
/*----------------------------------------------------------------------------*/
TracerS *init_tracer()
{
	TracerS *tracer = NULL;
	tracer = (TracerS*) malloc(sizeof(TracerS));
	if (tracer == NULL) goto on_error;
	tracer->Prev = NULL;
    tracer->Next = NULL;
    tracer->prop = (TracerPropS*) malloc(sizeof(TracerPropS));
	tracer->newList = NULL;
	return tracer;
    
on_error:
    ath_error("[init_mctracers]: Error allocating memory.\n");
    return NULL;
}

/*----------------------------------------------------------------------------*/
/* Remove tracer from list  */
/*----------------------------------------------------------------------------*/
void tracer_list_remove(TracerListS *list, TracerS *pnode)
{
    /* if pnode is only in list */
    if (!pnode->Prev && !pnode->Next) {
        assert(list->count == 1);
        assert(list->Head == pnode);
        assert(list->Tail == pnode);
        
        list->Head = NULL;
        list->Tail = NULL;
    }
    /* if pnode is last in list (but not first) */
    else if (!pnode->Next) {
        assert(list->Tail == pnode);
        
        list->Tail = pnode->Prev;
        list->Tail->Next = NULL;
    }
    /* if pnode is first in list (but not last) */
    else if (!pnode->Prev) {
        assert(list->Head == pnode);
        
        list->Head = pnode->Next;
        list->Head->Prev = NULL;
    }
    /* If in middle of list */
    else {
        assert(pnode->Prev->Next == pnode);
        assert(pnode->Next->Prev == pnode);
        
        pnode->Next->Prev = pnode->Prev;
        pnode->Prev->Next = pnode->Next;
    }
    list->count--;
    pnode->Next = NULL;
    pnode->Prev = NULL;
}

/*----------------------------------------------------------------------------*/
/* Iterate through all grids, free all tracers and tracer lists */
/*----------------------------------------------------------------------------*/
void tracer_destruct(MeshS *Mesh)
{
    GridS *pG;
    int nl, nd;
    int n1z, n2z, n3z, i, j, k, il, ir, jl, jr, kl, kr;
    int ks, ke, js, je, is, ie;
    TracerS *tracer;
    TracerS *Next;
    TracerListS *list;
    
    for (nl=0; nl<(Mesh->NLevels); nl++){
        for (nd=0; nd<(Mesh->DomainsPerLevel[nl]); nd++){
            if (Mesh->Domain[nl][nd].Grid != NULL){
                pG=(Mesh->Domain[nl][nd].Grid);
                is = pG->is;
                ie = pG->ie;
                js = pG->js;
                je = pG->je;
                ks = pG->ks;
                ke = pG->ke;
                if (pG->Nx[0] > 1) {
                    n1z = pG->Nx[0] + 2*nghost;
                    il = is - nghost;
                    ir = ie + nghost;
                }
                else {
                    n1z = 1;
                    il = is;
                    ir = ie;
                }
                if (pG->Nx[1] > 1) {
                    n2z = pG->Nx[1] + 2*nghost;
                    jl = js - nghost;
                    jr = je + nghost;
                }
                else {
                    n2z = 1;
                    jl = js;
                    jr = je;
                }
                if (pG->Nx[2] > 1) {
                    n3z = pG->Nx[2] + 2*nghost;
                    kl = ks - nghost;
                    kr = ke + nghost;
                }
                else {
                    n3z = 1;
                    kl = ks;
                    kr = ke;
                }
                for (k=kl; k<=kr; k++) {
                    for (j=jl; j<=jr; j++) {
                        for (i=il; i<=ir; i++) {
                            list = &((pG->GridLists)[k][j][i]);
                            tracer = list->Head;
                            while(tracer) {
                                Next = tracer->Next;
                                tracer_list_remove(list, tracer);
                                free(tracer);
                                tracer = Next;
                            }
                        }
                    }
                }
                free_3d_array(pG->GridLists);
            }
        }
    }
    
}

#endif /* TRACERS */
