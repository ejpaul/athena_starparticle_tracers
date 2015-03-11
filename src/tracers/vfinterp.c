#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../defs.h"
#include "../athena.h"
#include "../prototypes.h"
#include "prototypes.h"
#include "tracers.h"
#include "../globals.h"

#ifdef VFTRACERS
/*----------------------------------------------------------------------------*/
/*! \fn void getwei_linear(GridS *pG, Real x1, Real x2, Real x3,
 *           Real3Vect cell1, Real weight[3][3][3], int *is, int *js, int *ks)
 *  \brief Get weight using linear interpolation
 *
 * Input: pG: grid; x1,x2,x3: global coordinate; cell1: 1 over dx1,dx2,dx3
 * Output: weight: weight function; is,js,ks: starting cell indices in the grid.
 * Note: this interpolation works in any 1-3 dimensions.
 */

void getwei_linear(GridS *pG, Real x1, Real x2, Real x3, Real3Vect cell1,
                   Real weight[3][3][3], int *is, int *js, int *ks)
{
    int i, j, k, i1, j1, k1;
    Real a, b, c;		/* grid coordinate for the position (x1,x2,x3) */
    Real wei1[2], wei2[2], wei3[2];/* weight function in x1,x2,x3 directions */
    
    /* find cell locations and calculate 1D weight */
    /* x1 direction */
    if (cell1.x1 > 0.0) {
        /* Given x, returns containing cell i coordinate */
        /* x: global x coordinate;
         * i1: i-index containing x; a: grid index coordinate of x;
         * Return: 0: x is on the left of the ith cell;
         *         1: x is on the right of the ith cell;
         */
        i = celli(pG, x1, cell1.x1, &i1, &a);	/* x1 index */
        i1 = i+i1-1;	*is = i1;		/* starting x1 index */
        wei1[1] = a - i1 - 0.5;			/* one direction weight */
        wei1[0] = 1.0 - wei1[1];			/* 0: left; 1: right */
    }
    else { /* x1 dimension collapses */
        *is = pG->is;
        wei1[1] = 0.0;
        wei1[0] = 1.0;
    }
    
    /* x2 direction */
    if (cell1.x2 > 0.0) {
        j = cellj(pG, x2, cell1.x2, &j1, &b);	/* x2 index */
        j1 = j+j1-1;	*js = j1;		/* starting x2 index */
        wei2[1] = b - j1 - 0.5;			/* one direction weight */
        wei2[0] = 1.0 - wei2[1];			/* 0: left; 1: right */
    }
    else { /* x2 dimension collapses */
        *js = pG->js;
        wei2[1] = 0.0;
        wei2[0] = 1.0;
    }
    
    /* x3 direction */
    if (cell1.x3 > 0.0) {
        k = cellk(pG, x3, cell1.x3, &k1, &c);	/* x3 index */
        k1 = k+k1-1;	*ks = k1;		/* starting x3 index */
        wei3[1] = c - k1 - 0.5;			/* one direction weight */
        wei3[0] = 1.0 - wei3[1];			/* 0: left; 1: right */
    }
    else { /* x3 dimension collapses */
        *ks = pG->ks;
        wei3[1] = 0.0;
        wei3[0] = 1.0;
    }
    
    /* calculate 3D weight */
    for (k=0; k<2; k++)
        for (j=0; j<2; j++)
            for (i=0; i<2; i++)
                weight[k][j][i] = wei1[i] * wei2[j] * wei3[k];
    
    return;
}

/*----------------------------------------------------------------------------*/
/*! \fn void getwei_TSC(GridS *pG, Real x1, Real x2, Real x3, Real3Vect cell1,
 *                      Real weight[3][3][3], int *is, int *js, int *ks)
 *
 *  \brief Get weight using Triangular Shaped Cloud (TSC) interpolation
 * Input: pG: grid; x1,x2,x3: global coordinate; cell1: 1 over dx1,dx2,dx3
 * Output: weight: weight function; is,js,ks: starting cell indices in the grid.
 * Note: this interpolation works in any 1-3 dimensions.
 */
void getwei_TSC(GridS *pG, Real x1, Real x2, Real x3, Real3Vect cell1,
                Real weight[3][3][3], int *is, int *js, int *ks)
{
    int i, j, k;
    Real a, b, c, d;	/* grid coordinate for the position (x1,x2,x3) */
    Real wei1[3], wei2[3], wei3[3];/* weight function in x1,x2,x3 directions */
    
    /* find cell locations and calculate 1D weight */
    /* x1 direction */
    if (cell1.x1 > 0.0) {
        celli(pG, x1, cell1.x1, &i, &a);		/* x1 index */
        *is = i - 1;				/* starting x1 index, wei[0] */
        d = a - i;
        wei1[0] = 0.5*SQR(1.0-d);			/* 0: left; 2: right */
        wei1[1] = 0.75-SQR(d-0.5);			/* one direction weight */
        wei1[2] = 0.5*SQR(d);
    }
    else { /* x1 dimension collapses */
        *is = pG->is;
        wei1[1] = 0.0;	wei1[2] = 0.0;
        wei1[0] = 1.0;
    }
    
    /* x2 direction */
    if (cell1.x2 > 0.0) {
        cellj(pG, x2, cell1.x2, &j, &b);		/* x2 index */
        *js = j - 1;				/* starting x2 index */
        d = b - j;
        wei2[0] = 0.5*SQR(1.0-d);			/* 0: left; 2: right */
        wei2[1] = 0.75-SQR(d-0.5);			/* one direction weight */
        wei2[2] = 0.5*SQR(d);
    }
    else { /* x2 dimension collapses */
        *js = pG->js;
        wei2[1] = 0.0;	wei2[2] = 0.0;
        wei2[0] = 1.0;
    }
    
    /* x3 direction */
    if (cell1.x3 > 0.0) {
        cellk(pG, x3, cell1.x3, &k, &c);		/* x3 index */
        *ks = k - 1;				/* starting x3 index */
        d = c - k;
        wei3[0] = 0.5*SQR(1.0-d);			/* 0: left; 2: right */
        wei3[1] = 0.75-SQR(d-0.5);			/* one direction weight */
        wei3[2] = 0.5*SQR(d);
    }
    else { /* x3 dimension collapses */
        *ks = pG->ks;
        wei3[1] = 0.0;	wei3[2] = 0.0;
        wei3[0] = 1.0;
    }
    
    /* calculate 3D weight */
    for (k=0; k<3; k++) {
        for (j=0; j<3; j++) {
            for (i=0; i<3; i++) {
                weight[k][j][i] = wei1[i] * wei2[j] * wei3[k];
//                printf("weight[%d][%d][%d] = %f\n", k,j,i,weight[k][j][i]);
            }
        }
    }
    return;
}

/*----------------------------------------------------------------------------*/
/*! \fn void getwei_QP(GridS *pG, Real x1, Real x2, Real x3, Real3Vect cell1,
 *                       Real weight[3][3][3], int *is, int *js, int *ks)
 *  \brief Get weight using quadratic polynomial interpolation
 *
 * Input: pG: grid; x1,x2,x3: global coordinate; cell1: 1 over dx1,dx2,dx3
 * Output: weight: weight function; is,js,ks: starting cell indices in the grid.
 * Note: this interpolation works in any 1-3 dimensions.
 */
void getwei_QP(GridS *pG, Real x1, Real x2, Real x3, Real3Vect cell1,
               Real weight[3][3][3], int *is, int *js, int *ks)
{
    int i, j, k;
    Real a, b, c, d;	/* grid coordinate for the position (x1,x2,x3) */
    Real wei1[3], wei2[3], wei3[3];/* weight function in x1,x2,x3 directions */
    
    /* find cell locations and calculate 1D weight */
    /* x1 direction */
    if (cell1.x1 > 0.0) {
        celli(pG, x1, cell1.x1, &i, &a);		/* x1 index */
        *is = i - 1;				/* starting x1 index, wei[0] */
        d = a - i;
        wei1[0] = 0.5*(0.5-d)*(1.5-d);		/* 0: left; 2: right */
        wei1[1] = 1.0-SQR(d-0.5);			/* one direction weight */
        wei1[2] = 0.5*(d-0.5)*(d+0.5);
    }
    else { /* x1 dimension collapses */
        *is = pG->is;
        wei1[1] = 0.0;	wei1[2] = 0.0;
        wei1[0] = 1.0;
    }
    
    /* x2 direction */
    if (cell1.x2 > 0.0) {
        cellj(pG, x2, cell1.x2, &j, &b);		/* x2 index */
        *js = j - 1;				/* starting x2 index */
        d = b - j;
        wei2[0] = 0.5*(0.5-d)*(1.5-d);		/* 0: left; 2: right */
        wei2[1] = 1.0-SQR(d-0.5);			/* one direction weight */
        wei2[2] = 0.5*(d-0.5)*(d+0.5);
    }
    else { /* x2 dimension collapses */
        *js = pG->js;
        wei2[1] = 0.0;	wei2[2] = 0.0;
        wei2[0] = 1.0;
    }
    
    /* x3 direction */
    if (cell1.x3 > 0.0) {
        cellk(pG, x3, cell1.x3, &k, &c);		/* x3 index */
        *ks = k - 1;				/* starting x3 index */
        d = c - k;
        wei3[0] = 0.5*(0.5-d)*(1.5-d);		/* 0: left; 2: right */
        wei3[1] = 1.0-SQR(d-0.5);			/* one direction weight */
        wei3[2] = 0.5*(d-0.5)*(d+0.5);
    }
    else { /* x3 dimension collapses */
        *ks = pG->ks;
        wei3[1] = 0.0;	wei3[2] = 0.0;
        wei3[0] = 1.0;
    }
    
    /* calculate 3D weight */
    for (k=0; k<3; k++)
        for (j=0; j<3; j++)
            for (i=0; i<3; i++)
                weight[k][j][i] = wei1[i] * wei2[j] * wei3[k];
    
    return;
}
#endif /* VFTRACERS */
