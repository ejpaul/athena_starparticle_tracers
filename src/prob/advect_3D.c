#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"

/*----------------------------------------------------------------------------*/
/* problem:  */

void problem(DomainS *pDomain)
{
    GridS *pGrid=(pDomain->Grid);
#ifdef MCTRACERS
    MClistS *list;
#endif
    int i, is = pGrid->is, ie = pGrid->ie;
    int j, js = pGrid->js, je = pGrid->je;
    int k, ks = pGrid->ks, ke = pGrid->ke;
    
    for (k=ks; k<=ke; k++) {
        for (j=js; j<=je; j++) {
            for (i=is; i<=ie; i++) {
                if ((i-45)*(i-45) + (j-45)*(j-45) + (k-45)*(k-45) <= 1600) {
                    pGrid->U[k][j][i].d  = 2;
                    pGrid->U[k][j][i].M1 = 2;
                    pGrid->U[k][j][i].M2 = 0;
                    pGrid->U[k][j][i].M3 = 0;
                    pGrid->U[k][j][i].E = 1.0;
                    pGrid->U[k][j][i].B1c = 0;
                    pGrid->U[k][j][i].B2c = 0;
                    pGrid->U[k][j][i].B3c = 0;
#ifdef MCTRACERS
                    list = &((pGrid->GridLists)[k][j][i]);
                    init_mctracer_list(list, 40, 2);
#endif
                }
                else {
                    pGrid->U[k][j][i].d  = 1;
                    pGrid->U[k][j][i].M1 = 1;
                    pGrid->U[k][j][i].M2 = 0;
                    pGrid->U[k][j][i].M3 = 0;
                    pGrid->U[k][j][i].E = 0.5;
                    pGrid->U[k][j][i].B1c = 0;
                    pGrid->U[k][j][i].B2c = 0;
                    pGrid->U[k][j][i].B3c = 0;
                }
            }
        }
    }
    
#ifdef MCTRACERS
//    mc_init_proportional(pGrid);
#endif
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

#ifdef MCTRACERS
static Real num_density(const GridS *pG, const int i, const int j, const int k)
{
    Real count = (Real) (pG->GridLists)[k][j][i].count;
    return count;
}

static Real top_hat(const GridS *pG, const int i, const int j, const int k)
{
    Real count = (Real) (pG->TopHatGrid)[k][j][i].count;
    return count;
}

static Real ratio_map(const GridS *pG, const int i, const int j, const int k)
{
    Real count = (Real) (pG->GridLists)[k][j][i].count;
    Real d = (Real) pG->U[k][j][i].d;
    Real ratio = (Real) count/d;
    return ratio;
}

static Real d_init(const GridS *pG, const int i, const int j, const int k)
{
    Real d = 0;
    int n = 0;
    MCtracerS *tracer = pG->GridLists[k][j][i].Head;
    while (tracer) {
        n++;
        d += tracer->prop->d_init;
        tracer = tracer->Next;
    }
    return d/n;
}
#endif

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
#ifdef MCTRACERS
    if(strcmp(expr, "MC_num")==0) return num_density;
    if(strcmp(expr, "MC_tophat")==0) return top_hat;
    if(strcmp(expr, "ratio_map")==0) return ratio_map;
    if(strcmp(expr, "d_init")==0) return d_init;
#endif
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


