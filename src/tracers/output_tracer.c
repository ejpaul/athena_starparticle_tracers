#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../defs.h"
#include "../athena.h"
#include "../prototypes.h"
#include "prototypes.h"
#include "../globals.h"
#include <assert.h>
#ifdef VFTRACERS
/*----------------------------------------------------------------------------*/
/*! \fn void particle_to_grid(Grid *pG, PropFun_t par_prop)
 *  \brief Bin the particles to grid cells
 */
void tracer_to_grid(DomainS *pD)
{
    GridS *pG = pD->Grid;
    int i,j,k, is,js,ks, i0,j0,k0, i1,j1,k1, i2,j2,k2;
    long p;
    Real drho;
    Real weight[3][3][3];
    Real3Vect cell1;
    TracerS *tracer;
    int n0 = 2;
    /* Get grid limit related quantities */
    if (pG->Nx[0] > 1)  cell1.x1 = 1.0/pG->dx1;
    else                cell1.x1 = 0.0;
    
    if (pG->Nx[1] > 1)  cell1.x2 = 1.0/pG->dx2;
    else                cell1.x2 = 0.0;
    
    if (pG->Nx[2] > 1)  cell1.x3 = 1.0/pG->dx3;
    else                cell1.x3 = 0.0;

    /* initialization */
    for (k=klp; k<=kup; k++) {
        for (j=jlp; j<=jup; j++) {
            for (i=ilp; i<=iup; i++) {
                pG->TracerGrid[k][j][i].out_d = 0.0;
            }
        }
    }
    
    /* bin the particles */
    for (p=0; p<pG->ntracers; p++) {
        tracer = &(pG->tracers[p]);
        
        getwei_TSC(pG, tracer->x1, tracer->x2, tracer->x3, cell1, weight, &is, &js, &ks);
        
        /* distribute particles */
        k1 = MAX(ks, klp);    k2 = MIN(ks+n0, kup);
        j1 = MAX(js, jlp);    j2 = MIN(js+n0, jup);
        i1 = MAX(is, ilp);    i2 = MIN(is+n0, iup);

        for (k=k1; k<=k2; k++) {
            k0 = k-k1;
            for (j=j1; j<=j2; j++) {
                j0 = j-j1;
                for (i=i1; i<=i2; i++) {
                    i0 = i-i1;
                    pG->TracerGrid[k][j][i].out_d  += weight[k0][j0][i0];
                    pG->TracerGrid[k][j][i].d_init += weight[k0][j0][i0]*(tracer->d_init);
                }
            }
        }
    }
    return;
}
#endif