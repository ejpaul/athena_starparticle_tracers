#ifndef TRACERS_PROTOTYPES_H
#define TRACERS_PROTOTYPES_H
#include "../copyright.h"
/*============================================================================*/
/*! \file prototypes.h
 *  \brief Prototypes for all public functions in the /src/mctracers dir */
/*============================================================================*/

#include <stdio.h>
#include <stdarg.h>
#include "../athena.h"
#include "../defs.h"
#include "../config.h"
#include "tracers.h"

#if defined(VFTRACERS) || defined(MCTRACERS)

/* init_tracers.c */
void init_tracer_grid(MeshS *pM);
void tracer_init_unif(GridS *pG);
void tracer_init_proportional(GridS *pG);
void tracer_init_xlinflow(GridS *pG);
void tracer_debug(GridS *pG);
double get_tracer_id();
void init_tracer_list(TracerListS *list, int n, Real d);
void Tracerlist_add(TracerListS *list, TracerS *tracer);
TracerS *init_tracer();
void tracer_list_remove(TracerListS *list, TracerS *pnode);
void tracer_destruct(MeshS *Mesh);
void init_vftracer_list(GridS *pG, TracerListS *list, int N);

/* integrate_tracers.c */
void Tracerlist_sweep(TracerListS *list, GridS *pG);
void Tracerlist_sweep_bc(TracerListS *list);
void prob_iterate_x1(TracerListS *list, double pflux, GridS *pG);
void prob_iterate_x2(TracerListS *list, double pflux, GridS *pG);
void prob_iterate_x3(TracerListS *list, double pflux, GridS *pG);
void ran_gen_list(GridS *pG);
void ran_gen(GridS *pG);

/* vfintegrate.c */
void vf_newijk(DomainS* pD);
void Integrate_vftracers(DomainS *pD);
void Integrate_vf_2nd(DomainS *pD);
void vf_newpos(GridS *pG, TracerListS *list, TracerListS *newList);
int  interp(GridS *pG, Real weight[3][3][3], int is, int js, int ks, Real *u1, Real *u2, Real *u3);

/* bvals_tracer.c */
void bvals_tracer(DomainS *pD);
void bvals_tracer_init(MeshS *pM);

/* vfinterp.c */
void getwei_TSC(GridS *pG, Real x1, Real x2, Real x3, Real3Vect cell1,
                Real weight[3][3][3], int *is, int *js, int *ks);

#endif /* TRACERS */

#endif /* TRACERS_PROTYPES_H */


