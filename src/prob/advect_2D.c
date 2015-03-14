#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"

void problem(DomainS *pDomain)
{
    GridS *pGrid=(pDomain->Grid);
    int i, is = pGrid->is, ie = pGrid->ie;
    int j, js = pGrid->js, je = pGrid->je;
    int k, ks = pGrid->ks, ke = pGrid->ke;
    
    for (k=ks; k<=ke; k++) {
        for (j=js; j<=je; j++) {
            for (i=is; i<=ie; i++) {
                if (((i-50)*(i-50) + (j-50)*(j-50)) <= 100) {
                    pGrid->U[k][j][i].d  = 2;
                    pGrid->U[k][j][i].M1 = 0;
                    pGrid->U[k][j][i].M2 = 2;
                    pGrid->U[k][j][i].M3 = 0;
                    pGrid->U[k][j][i].E = 1.0;
                    pGrid->U[k][j][i].B1c = 0;
                    pGrid->U[k][j][i].B2c = 0;
                    pGrid->U[k][j][i].B3c = 0;
                }
                else {
                    pGrid->U[k][j][i].d  = 1;
                    pGrid->U[k][j][i].M1 = 0;
                    pGrid->U[k][j][i].M2 = 1;
                    pGrid->U[k][j][i].M3 = 0;
                    pGrid->U[k][j][i].E = 0.5;
                    pGrid->U[k][j][i].B1c = 0;
                    pGrid->U[k][j][i].B2c = 0;
                    pGrid->U[k][j][i].B3c = 0;
                }
            }
        }
    }
    
#if defined(MCTRACERS) || defined(VFTRACERS)
    if (strcmp(par_gets("problem","distribution"), "uniform") == 0)
        tracer_init_unif(pGrid);
    else if (strcmp(par_gets("problem","distribution"), "prop") == 0)
        tracer_init_proportional(pGrid);
#endif // TRACERS //
    
}

/*==============================================================================
 * PROBLEM USER FUNCTIONS:
 * problem_write_restart() - writes problem-specific user data to restart files
 * problem_read_restart()  - reads problem-specific user data from restart files
 * get_usr_expr()          - sets pointer to expression for special output data
 * get_usr_out_fun()       - returns a user defined output function pointer
 * get_usr_par_prop()      - returns a user defined particle selection function
 * Userwork_in_loop        - problem specific work IN     main loop
 * Userwork_after_loop     - problem specific work AFTER  main loop
 *----------------------------------------------------------------------------*/

#if defined(MCTRACERS) || defined(VFTRACERS)
static Real num_density(const GridS *pG, const int i, const int j, const int k)
{
    Real count = (Real) (pG->GridLists)[k][j][i].count;
    return count;
}
#endif /* TRACERS */


#ifdef MCTRACERS
#ifdef TOPHAT
static Real top_hat(const GridS *pG, const int i, const int j, const int k)
{
    Real count = (Real) (pG->TopHatGrid)[k][j][i].count;
    return count;
}
#endif /* TOPHAT */
#endif // MCTRACERS //

#if defined(VFTRACERS) || defined(MCTRACERS)
static Real ratio_map(const GridS *pG, const int i, const int j, const int k)
{
    Real count = (Real) (pG->GridLists)[k][j][i].count;
    Real d = (Real) pG->U[k][j][i].d;
    Real ratio = (Real) count/d;
    return ratio;
}
#endif // TRACERS //

#if defined(MCTRACERS) || defined(VFTRACERS)
static Real d_init(const GridS *pG, const int i, const int j, const int k)
{
    Real d = 0.0;
    Real n = 0.0;
    TracerS *tracer = pG->GridLists[k][j][i].Head;
    while (tracer) {
        n++;
        d += tracer->prop->d_init;
        tracer = tracer->Next;
    }
    if (n != 0){
        return d/n;
    }
    else return 0;
}
#endif // MCTRACERS //

#if defined(MCTRACERS) || defined(VFTRACERS)
static Real i_init(const GridS *pG, const int i, const int j, const int k)
{
    Real d = 0.0;
    Real n = 0.0;
    TracerS *tracer = pG->GridLists[k][j][i].Head;
    while (tracer) {
        n++;
        d += tracer->prop->d_init;
        tracer = tracer->Next;
    }
    if (n != 0){
        return d/n;
    }
    else return 0;
}
#endif // TRACERS //

void problem_write_restart(MeshS *pM, FILE *fp)
{
    return;
}

void problem_read_restart(MeshS *pM, FILE *fp)
{
    return;
}

ConsFun_t get_usr_expr(const char *expr)
{
#if defined(MCTRACERS) || defined(VFTRACERS)
    if(strcmp(expr, "num_density")==0) return num_density;
    if(strcmp(expr, "d_init")==0) return d_init;
    if(strcmp(expr, "ratio_map")==0) return d_init;
#endif /* TRACERS */
    return NULL;
}

VOutFun_t get_usr_out_fun(const char *name){
    return NULL;
}

void Userwork_in_loop(MeshS *pM)
{
}

void Userwork_after_loop(MeshS *pM)
{
}


