#include "../copyright.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../defs.h"
#include "../athena.h"
#include "../prototypes.h"
#include "prototypes.h"
#include "tracers.h"
#include "../globals.h"
#include <assert.h>

#ifdef VFTRACERS

/*=========================== PUBLIC FUNCTIONS ===============================*/
void Integrate_vftracers(DomainS *pD)
{
    GridS *pG = pD->Grid;         /* set ptr to Grid */

    TracerS *curr;                /* pointer of the current working position */
    TracerListS *list;
    Real3Vect cell1;              /* one over dx1, dx2, dx3 */
    int inew, jnew, knew;
    Real weight[3][3][3];		  /* weight function */
    /* interpolated current velocity */
    Real v1 = 0;
    Real v2 = 0;
    Real v3 = 0;
    /* interpolated predicted velocity */
    Real u1 = 0;
    Real u2 = 0;
    Real u3 = 0;
    /* predicted position */
    Real xp1 = 0;
    Real xp2 = 0;
    Real xp3 = 0;

    const Real dx1_1 = 1./pG->dx1;
    const Real dx2_1 = 1./pG->dx2;
    const Real dx3_1 = 1./pG->dx3;
    
    int is = pG->is;
    int ie = pG->ie;
    int js = pG->js;
    int je = pG->je;
    int ks = pG->ks;
    int ke = pG->ke;
    int i_s, j_s, k_s;
    int i1, i2, i3;
    Real a, b, c;
    
    /* cell1 is a shortcut expressions as well as dimension indicator */
    if (pG->Nx[0] > 1)  cell1.x1 = 1.0/pG->dx1;  else cell1.x1 = 0.0;
    if (pG->Nx[1] > 1)  cell1.x2 = 1.0/pG->dx2;  else cell1.x2 = 0.0;
    if (pG->Nx[2] > 1)  cell1.x3 = 1.0/pG->dx3;  else cell1.x3 = 0.0;
    
    for (int k=ks; k<=ke; k++) {
        for (int j=js; j<=je; j++) {
            for (int i=is; i<=ie; i++) {
                list = &((pG->GridLists)[k][j][i]);
                curr = list->Head;
                while(curr) {
                    /* Interpolate velocity at current position */
                    getwei_TSC(pG, curr->x1, curr->x2, curr->x3, cell1, weight, &i_s, &j_s, &k_s);
                    interp(pG, weight, i_s, j_s, k_s, &v1, &v2, &v3);
                    
                    /* Predicted positions */
                    xp1 = curr->x1 + pG->dt*v1;
                    xp2 = curr->x2 + pG->dt*v2;
                    xp3 = curr->x3 + pG->dt*v3;
                                        
                    /* Interpolate velocity at predicted position */
                    getwei_TSC(pG, xp1, xp2, xp3, cell1, weight, &i_s, &j_s, &k_s);
                    interp(pG, weight, i_s, j_s, k_s, &u1, &u2, &u3);
                    
                    /* Improved Euler method */
                    if (pG->Nx[0] > 1) {
                        curr->x1 += 0.5*pG->dt*v1 + 0.5*pG->dt*u1;
                    }
                    if (pG->Nx[1] > 1) {
                        curr->x2 += 0.5*pG->dt*v2 + 0.5*pG->dt*u2;
                    }
                    if (pG->Nx[2] > 1) {
                        curr->x3 += 0.5*pG->dt*v3 + 0.5*pG->dt*u3;
                    }
                    curr = curr->Next;
                }
            }
        }
    }
}

void Integrate_vf_2nd(DomainS *pD)
{
    GridS *pG = pD->Grid;         /* set ptr to Grid */
    
    TracerS *curr;                /* pointer of the current working position */
    TracerListS *list;
    Real3Vect cell1;              /* one over dx1, dx2, dx3 */
    int inew, jnew, knew;
    Real weight[3][3][3];		  /* weight function */
    /* interpolated current velocity */
    Real v1 = 0;
    Real v2 = 0;
    Real v3 = 0;
    /* interpolated previous velocity */
    Real u1 = 0;
    Real u2 = 0;
    Real u3 = 0;
    /* predicted position */
    Real xp1 = 0;
    Real xp2 = 0;
    Real xp3 = 0;
    
    const Real dx1_1 = 1./pG->dx1;
    const Real dx2_1 = 1./pG->dx2;
    const Real dx3_1 = 1./pG->dx3;
    
    int is = pG->is;
    int ie = pG->ie;
    int js = pG->js;
    int je = pG->je;
    int ks = pG->ks;
    int ke = pG->ke;
    int i_s, j_s, k_s;
    int i1, i2, i3;
    Real a, b, c;
    
    /* cell1 is a shortcut expressions as well as dimension indicator */
    if (pG->Nx[0] > 1)  cell1.x1 = 1.0/pG->dx1;  else cell1.x1 = 0.0;
    if (pG->Nx[1] > 1)  cell1.x2 = 1.0/pG->dx2;  else cell1.x2 = 0.0;
    if (pG->Nx[2] > 1)  cell1.x3 = 1.0/pG->dx3;  else cell1.x3 = 0.0;
    
    for (int k=ks; k<=ke; k++) {
        for (int j=js; j<=je; j++) {
            for (int i=is; i<=ie; i++) {
                list = &((pG->GridLists)[k][j][i]);
                curr = list->Head;
                while(curr) {
                    /* Predicted positions */
                    xp1 = curr->x1 + pG->dt*v1;
                    xp2 = curr->x2 + pG->dt*v2;
                    xp3 = curr->x3 + pG->dt*v3;

                    /* Interpolate Uprev to predicted position */
                    getwei_TSC(pG, xp1, xp2, xp3, cell1, weight, &i_s, &j_s, &k_s);
                    interp(pG, weight, i_s, j_s, k_s, &v1, &v2, &v3);
                    
                    /* Interpolate U to predicted position */
                    interp_prev(pG, weight, i_s, j_s, k_s, &u1, &u2, &u3);
                    
                    /* 2nd order method */
                    if (pG->Nx[0] > 1) {
                        curr->x1 += 0.5*pG->dt*v1 + 0.5*pG->dt*u1;
                    }
                    if (pG->Nx[1] > 1) {
                        curr->x2 += 0.5*pG->dt*v2 + 0.5*pG->dt*u2;
                    }
                    if (pG->Nx[2] > 1) {
                        curr->x3 += 0.5*pG->dt*v3 + 0.5*pG->dt*u3;
                    }
                    curr = curr->Next;
                }
            }
        }
    }
}

void Integrate_vf_reconstruction(DomainS *pD)
{
    GridS *pG = pD->Grid;         /* set ptr to Grid */
    TracerS *curr;                /* pointer of the current working position */
    TracerListS *list;
    Real3Vect cell1;              /* one over dx1, dx2, dx3 */
    int inew, jnew, knew;
    Real weight[3][3][3];		  /* weight function */
    /* interpolated current velocity */
    Real v1 = 0;
    Real v2 = 0;
    Real v3 = 0;
    /* interpolated predicted velocity */
    Real u1 = 0;
    Real u2 = 0;
    Real u3 = 0;
    /* predicted position */
    Real xp1 = 0;
    Real xp2 = 0;
    Real xp3 = 0;
    
    const Real dx1_1 = 1./pG->dx1;
    const Real dx2_1 = 1./pG->dx2;
    const Real dx3_1 = 1./pG->dx3;
    
    int is = pG->is;
    int ie = pG->ie;
    int js = pG->js;
    int je = pG->je;
    int ks = pG->ks;
    int ke = pG->ke;
    int i_s, j_s, k_s;
    int i1, i2, i3;
    Real a, b, c;
    
    /* cell1 is a shortcut expressions as well as dimension indicator */
    if (pG->Nx[0] > 1)  cell1.x1 = 1.0/pG->dx1;  else cell1.x1 = 0.0;
    if (pG->Nx[1] > 1)  cell1.x2 = 1.0/pG->dx2;  else cell1.x2 = 0.0;
    if (pG->Nx[2] > 1)  cell1.x3 = 1.0/pG->dx3;  else cell1.x3 = 0.0;
    
    for (int k=ks; k<=ke; k++) {
        for (int j=js; j<=je; j++) {
            for (int i=is; i<=ie; i++) {
                list = &((pG->GridLists)[k][j][i]);
                curr = list->Head;
                while(curr) {
                    /* Predicted positions */
                    xp1 = curr->x1 + pG->dt*v1;
                    xp2 = curr->x2 + pG->dt*v2;
                    xp3 = curr->x3 + pG->dt*v3;
                    
                    /* Interpolate Uprev to predicted position */
                    getwei_TSC(pG, xp1, xp2, xp3, cell1, weight, &i_s, &j_s, &k_s);
                    interp(pG, weight, i_s, j_s, k_s, &v1, &v2, &v3);
                    
                    /* Interpolate U to predicted position */
                    interp_prev(pG, weight, i_s, j_s, k_s, &u1, &u2, &u3);
                    
                    /* 2nd order method */
                    if (pG->Nx[0] > 1) {
                        curr->x1 += 0.5*pG->dt*v1 + 0.5*pG->dt*u1;
                    }
                    if (pG->Nx[1] > 1) {
                        curr->x2 += 0.5*pG->dt*v2 + 0.5*pG->dt*u2;
                    }
                    if (pG->Nx[2] > 1) {
                        curr->x3 += 0.5*pG->dt*v3 + 0.5*pG->dt*u3;
                    }
                    curr = curr->Next;
                }
            }
        }
    }
}


void vf_newijk(DomainS* pD) {
    GridS *pG = (pD->Grid);
    int is = pG->is;
    int ie = pG->ie;
    int js = pG->js;
    int je = pG->je;
    int ks = pG->ks;
    int ke = pG->ke;
    TracerS *tracer;
    TracerListS *list;
    Real3Vect cell1;              /* one over dx1, dx2, dx3 */
    int inew, jnew, knew;
    Real a, b, c;
    /* cell1 is a shortcut expressions as well as dimension indicator */
    if (pG->Nx[0] > 1)  cell1.x1 = 1.0/pG->dx1;  else cell1.x1 = 0.0;
    if (pG->Nx[1] > 1)  cell1.x2 = 1.0/pG->dx2;  else cell1.x2 = 0.0;
    if (pG->Nx[2] > 1)  cell1.x3 = 1.0/pG->dx3;  else cell1.x3 = 0.0;
    
    for (int k=ks; k<=ke; k++) {
        for (int j=js; j<=je; j++) {
            for (int i=is; i<=ie; i++) {
                list = &((pG->GridLists)[k][j][i]);
                tracer = list->Head;
                while(tracer) {
                    celli(pG, tracer->x1, cell1.x1, &inew, &a);
                    cellj(pG, tracer->x2, cell1.x2, &jnew, &b);
                    cellk(pG, tracer->x3, cell1.x3, &knew, &c);
                    tracer->newList = &((pG->GridLists)[knew][jnew][inew]);
                    tracer = tracer->Next;
                }
                Tracerlist_sweep(list, pG);
            }
        }
    }
}

void vf_newpos(GridS* pG, TracerListS *list, TracerListS *newList) {
    TracerS *tracer = list->Head;
    double x1, x2, x3, dx1, dx2, dx3, x1new, x2new, x3new;
    int inew = newList->i;
    int i = list->i;
    int jnew = newList->j;
    int j = list->j;
    int knew = newList->k;
    int k = list->k;
    int i1, i2, i3;
    
    Real a, b, c;
    Real3Vect cell1;              /* one over dx1, dx2, dx3 */
    /* cell1 is a shortcut expressions as well as dimension indicator */
    if (pG->Nx[0] > 1)  cell1.x1 = 1.0/pG->dx1;  else cell1.x1 = 0.0;
    if (pG->Nx[1] > 1)  cell1.x2 = 1.0/pG->dx2;  else cell1.x2 = 0.0;
    if (pG->Nx[2] > 1)  cell1.x3 = 1.0/pG->dx3;  else cell1.x3 = 0.0;
    
    cc_pos(pG, inew, jnew, knew, &x1new, &x2new, &x3new);
    cc_pos(pG, i, j, k, &x1, &x2, &x3);
    
    while(tracer) {
        dx1 = tracer->x1 - x1;
        dx2 = tracer->x2 - x2;
        dx3 = tracer->x3 - x3;
        tracer->x1 = x1new + dx1;
        tracer->x2 = x2new + dx2;
        tracer->x3 = x3new + dx3;
        tracer = tracer->Next;
    }
}

int interp(GridS *pG, Real weight[3][3][3], int is, int js, int ks, Real *u1, Real *u2, Real *u3)
{
    int n0,i,j,k,i0,j0,k0,i1,j1,k1,i2,j2,k2;
    Real v1, v2, v3;
    Real totwei, totwei1;		/* total weight (in case of edge cells) */
    
    /* linear interpolation */
    v1 = 0.0; v2 = 0.0; v3 = 0.0;
    totwei = 0.0; totwei1 = 1.0;
    n0 = 2;
    
    int ie = pG->ie;
    int je = pG->je;
    int ke = pG->ke;
    /* Interpolate density, velocity and sound speed */
    /* Note: in lower dimensions only wei[0] is non-zero */
    k1 = ks;	k2 = MIN(ks+n0, ke);
    j1 = js;	j2 = MIN(js+n0, je);
    i1 = je;	i2 = MIN(is+n0, ie);
    
    for (k=k1; k<=k2; k++) {
        k0=k-k1;
        for (j=j1; j<=j2; j++) {
            j0=j-j1;
            for (i=i1; i<=i2; i++) {
                i0=i-i1;
                v1 += weight[k0][j0][i0] * (pG->U[k][j][i].M1)/(pG->U[k][j][i].d);
                v2 += weight[k0][j0][i0] * (pG->U[k][j][i].M2)/(pG->U[k][j][i].d);
                v3 += weight[k0][j0][i0] * (pG->U[k][j][i].M3)/(pG->U[k][j][i].d);
                totwei += weight[k0][j0][i0];
            }
        }
    }
    if (totwei < TINY_NUMBER) /* particle lies out of the grid, warning! */
        return -1;
    
    totwei1 = 1.0/totwei;
    *u1 = v1*totwei1;	*u2 = v2*totwei1;	*u3 = v3*totwei1;

    return 0;
}

int interp_prev(GridS *pG, Real weight[3][3][3], int is, int js, int ks, Real *u1, Real *u2, Real *u3)
{
    int n0,i,j,k,i0,j0,k0,i1,j1,k1,i2,j2,k2;
    Real v1, v2, v3;
    Real totwei, totwei1;		/* total weight (in case of edge cells) */
    
    /* linear interpolation */
    v1 = 0.0; v2 = 0.0; v3 = 0.0;
    totwei = 0.0; totwei1 = 1.0;
    n0 = 2;
    
    int ie = pG->ie;
    int je = pG->je;
    int ke = pG->ke;
    /* Interpolate density, velocity and sound speed */
    /* Note: in lower dimensions only wei[0] is non-zero */
    k1 = ks;	k2 = MIN(ks+n0, ke);
    j1 = js;	j2 = MIN(js+n0, je);
    i1 = je;	i2 = MIN(is+n0, ie);
    
    for (k=k1; k<=k2; k++) {
        k0=k-k1;
        for (j=j1; j<=j2; j++) {
            j0=j-j1;
            for (i=i1; i<=i2; i++) {
                i0=i-i1;
                v1 += weight[k0][j0][i0] * (pG->U_prev[k][j][i].M1)/(pG->U[k][j][i].d);
                v2 += weight[k0][j0][i0] * (pG->U_prev[k][j][i].M2)/(pG->U[k][j][i].d);
                v3 += weight[k0][j0][i0] * (pG->U_prev[k][j][i].M3)/(pG->U[k][j][i].d);
                totwei += weight[k0][j0][i0];
            }
        }
    }
    if (totwei < TINY_NUMBER) /* particle lies out of the grid, warning! */
        return -1;
    
    totwei1 = 1.0/totwei;
    *u1 = v1*totwei1;	*u2 = v2*totwei1;	*u3 = v3*totwei1;
    
    return 0;
}

#endif /* VFTRACERS */