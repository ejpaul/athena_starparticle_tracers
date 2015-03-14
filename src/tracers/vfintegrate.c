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
void Integrate_vf_2nd(DomainS *pD)
{
    GridS *pG = pD->Grid;         /* set ptr to Grid */
    
    TracerS *curr;                /* pointer of the current working position */
    TracerListS *list;
    Real3Vect cell1;              /* one over dx1, dx2, dx3 */
    
    Real weight[3][3][3];		  /* weight function */
    /* position x^n */
    Real x1, x2, x3;
    /* U^(n+1)[x^n] */
    Real v1 = 0;
    Real v2 = 0;
    Real v3 = 0;
    /* U^n[x^n] */
    Real v1_prev = 0;
    Real v2_prev = 0;
    Real v3_prev = 0;
    /* U^n[x^(n+1/2)] */
    Real vp1_prev = 0;
    Real vp2_prev = 0;
    Real vp3_prev = 0;
    /* U^(n+1)[x^(n+1/2)] */
    Real vp1 = 0;
    Real vp2 = 0;
    Real vp3 = 0;
    /* predicted position */
    Real xp1 = 0;
    Real xp2 = 0;
    Real xp3 = 0;
    
    int is = pG->is;
    int ie = pG->ie;
    int js = pG->js;
    int je = pG->je;
    int ks = pG->ks;
    int ke = pG->ke;
    int i_s, j_s, k_s;
    
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
                    /* Interpolate U^n to x^n */
                    x1 = curr->x1;
                    x2 = curr->x2;
                    x3 = curr->x3;
                    getwei_TSC(pG, x1, x2, x3, cell1, weight, &i_s, &j_s, &k_s);
                    interp_prev(pG, weight, i_s, j_s, k_s, &v1_prev, &v2_prev, &v3_prev);
                    
                    /* Interpolate U^(n+1) to x^n */
                    interp(pG, weight, i_s, j_s, k_s, &v1, &v2, &v3);
                    
                    /* Predicted positions */
                    /* v^(n+1/4) = 3/4*v^n + 1/4*v^(n+1) */
                    xp1 = x1 + 0.5*pG->dt*(0.75*v1_prev + 0.25*v1);
                    xp2 = x2 + 0.5*pG->dt*(0.75*v2_prev + 0.25*v2);
                    xp3 = x3 + 0.5*pG->dt*(0.75*v3_prev + 0.25*v3);

                    /* Interpolate U^n to x^(n+1/2) */
                    getwei_TSC(pG, xp1, xp2, xp3, cell1, weight, &i_s, &j_s, &k_s);
                    interp_prev(pG, weight, i_s, j_s, k_s, &vp1_prev, &vp2_prev, &vp3_prev);
                    
                    /* Interpolate U^(n+1) to x^(n+1/2) */
                    interp(pG, weight, i_s, j_s, k_s, &vp1, &vp2, &vp3);
                    
                    /* 2nd order method */
                    if (pG->Nx[0] > 1) {
                        curr->x1 += 0.5*pG->dt*(0.25*vp1_prev + 0.75*vp1);
                    }
                    if (pG->Nx[1] > 1) {
                        curr->x2 += 0.5*pG->dt*(0.25*vp2_prev + 0.75*vp2);
                    }
                    if (pG->Nx[2] > 1) {
                        curr->x3 += 0.5*pG->dt*(0.25*vp3_prev + 0.75*vp3);
                    }
                    curr = curr->Next;
                }
            }
        }
    }
}

void Integrate_vf_2nd_lower(DomainS *pD) {
    GridS *pG = pD->Grid;         /* set ptr to Grid */
    
    TracerS *curr;                /* pointer of the current working position */
    TracerListS *list;
    Real3Vect cell1;              /* one over dx1, dx2, dx3 */
    
    Real weight[3][3][3];		  /* weight function */
    /* position x^n */
    Real x1, x2, x3;
    /* U^(n)[x^n] */
    Real v1_prev = 0;
    Real v2_prev = 0;
    Real v3_prev = 0;
    /* U^(n)[x^(n+1/2)] */
    Real vp1_prev = 0;
    Real vp2_prev = 0;
    Real vp3_prev = 0;
    /* U^(n+1)[x^(n+1/2)] */
    Real vp1 = 0;
    Real vp2 = 0;
    Real vp3 = 0;
    /* predicted position */
    Real xp1 = 0;
    Real xp2 = 0;
    Real xp3 = 0;
    
    int is = pG->is;
    int ie = pG->ie;
    int js = pG->js;
    int je = pG->je;
    int ks = pG->ks;
    int ke = pG->ke;
    int i_s, j_s, k_s;
    
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
                    /* Interpolate U^n to x^n */
                    x1 = curr->x1;
                    x2 = curr->x2;
                    x3 = curr->x3;
                    getwei_TSC(pG, x1, x2, x3, cell1, weight, &i_s, &j_s, &k_s);
                    interp_prev(pG, weight, i_s, j_s, k_s, &v1_prev, &v2_prev, &v3_prev);
                    
                    /* Predicted positions */ 
                    xp1 = x1 + 0.5*pG->dt*(v1_prev);
                    xp2 = x2 + 0.5*pG->dt*(v2_prev);
                    xp3 = x3 + 0.5*pG->dt*(v3_prev);
                    
                    /* Interpolate U^n to x^(n+1/2) */
                    getwei_TSC(pG, xp1, xp2, xp3, cell1, weight, &i_s, &j_s, &k_s);
                    interp_prev(pG, weight, i_s, j_s, k_s, &vp1_prev, &vp2_prev, &vp3_prev);
                    
                    /* Interpolate U^(n+1) to x^(n+1/2) */
                    interp(pG, weight, i_s, j_s, k_s, &vp1, &vp2, &vp3);
                    
                    /* 2nd order method */
                    if (pG->Nx[0] > 1) {
                        curr->x1 += pG->dt*(0.5*vp1_prev + 0.5*vp1);
                    }
                    if (pG->Nx[1] > 1) {
                        curr->x2 += pG->dt*(0.5*vp2_prev + 0.5*vp2);
                    }
                    if (pG->Nx[2] > 1) {
                        curr->x3 += pG->dt*(0.5*vp3_prev + 0.5*vp3);
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

/* Function to move vftracer from ghost zone position in list to 
   active zone position in newList */

void vf_newpos(GridS* pG, TracerListS *list, TracerListS *newList) {

    TracerS *tracer = list->Head;
    double x1, x2, x3, dx1, dx2, dx3, x1new, x2new, x3new;
    int inew = newList->i;
    int i = list->i;
    int jnew = newList->j;
    int j = list->j;
    int knew = newList->k;
    int k = list->k;
    int ie = pG->ie;
    int je = pG->je;
    int ke = pG->ke;
    int is = pG->is;
    int js = pG->js;
    int ks = pG->ks;
    
    Real3Vect cell1;              /* one over dx1, dx2, dx3 */
    /* cell1 is a shortcut expressions as well as dimension indicator */
    if (pG->Nx[0] > 1)  cell1.x1 = 1.0/pG->dx1;  else cell1.x1 = 0.0;
    if (pG->Nx[1] > 1)  cell1.x2 = 1.0/pG->dx2;  else cell1.x2 = 0.0;
    if (pG->Nx[2] > 1)  cell1.x3 = 1.0/pG->dx3;  else cell1.x3 = 0.0;
    
    /* Given i, j, k returns cell-centered x1, x2, x3 position */
    cc_pos(pG, inew, jnew, knew, &x1new, &x2new, &x3new);
    cc_pos(pG, i, j, k, &x1, &x2, &x3);

    assert(i < is || i > ie || j < js || j > je || k < ks || k > ke);
    assert(inew >= is && jnew >= js || knew >= ks);
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
    
    /* TSC interpolation */
    v1 = 0.0; v2 = 0.0; v3 = 0.0;
    totwei = 0.0; totwei1 = 1.0;
    n0 = 3;
    
//    int is = pG->is;
//    int js = pG->js;
//    int ks = pG->ks;
    int ie = pG->ie;
    int je = pG->je;
    int ke = pG->ke;
    
    /* Note: in lower dimensions only wei[0] is non-zero */
    k1 = ks;	k2 = MIN(ks+n0, ke);
    j1 = js;	j2 = MIN(js+n0, je);
    i1 = is;	i2 = MIN(is+n0, ie);
    
//    if (k1 == (pG->ks-1)) k1 = pG->ke;
//    else k1 = ks;
//    k2 = ks+1;
//    if (ks+2 == (pG->ke+1)) k3 = pG->ks;
//    else k3 = k2+2;
    
    for (k=k1; k<=k2; k++) {
        k0=k-k1;
        for (j=j1; j<=j2; j++) {
            j0=j-j1;
            for (i=i1; i<=i2; i++) {
                i0 = i-i1;

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
    n0 = 3;
    
//    int is = pG->is;
//    int js = pG->js;
//    int ks = pG->ks;
    int ie = pG->ie;
    int je = pG->je;
    int ke = pG->ke;
    /* Interpolate density, velocity and sound speed */
    /* Note: in lower dimensions only wei[0] is non-zero */
    k1 = ks;	k2 = MIN(ks+n0, ke);
    j1 = js;	j2 = MIN(js+n0, je);
    i1 = is;	i2 = MIN(is+n0, ie);
    
    for (k=k1; k<=k2; k++) {
        k0=k-k1;
        for (j=j1; j<=j2; j++) {
            j0=j-j1;
            for (i=i1; i<=i2; i++) {
                i0=i-i1;
//                if (i == (ie+1)) i = is;
//                if (i == (is-1)) i = ie;
//                if (j == (je+1)) j = is;
//                if (j == (js-1)) j = je;
//                if (k == (ke+1)) k = ks;
//                if (k == (ks-1)) k = ke;
                v1 += weight[k0][j0][i0] * (pG->U_prev[k][j][i].M1)/(pG->U_prev[k][j][i].d);
                v2 += weight[k0][j0][i0] * (pG->U_prev[k][j][i].M2)/(pG->U_prev[k][j][i].d);
                v3 += weight[k0][j0][i0] * (pG->U_prev[k][j][i].M3)/(pG->U_prev[k][j][i].d);
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
