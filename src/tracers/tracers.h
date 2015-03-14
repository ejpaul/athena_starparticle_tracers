#ifndef TRACERS_H
#define TRACERS_H
#include "../copyright.h"
#include <stdio.h>
#include <stdarg.h>
#include "../athena.h"
#include "../defs.h"
#include "../config.h"

#if defined(MCTRACERS) || defined(VFTRACERS)

typedef struct TracerListS TracerListS;
typedef struct TracerPropS TracerPropS;
typedef struct TracerS TracerS;
typedef struct TracertophatS TracertophatS;
typedef struct TracerList_MPI TracerList_MPI;
typedef struct Tracer_MPI Tracer_MPI;

/* Monte Carlo tracer particle structure */
struct TracerS{
    TracerListS* newList; /*!< pointer to new list, if !=NULL flags for movement */
    TracerS* Next;  /*!< pointer to next tracer in list */
    TracerS* Prev;  /*!< pointer to previous tracer in list */
    TracerPropS* prop; /*!< data structure storing properties of tracer */
#ifdef VFTRACERS
    double x1,x2,x3;	/*!< coordinate in X,Y,Z */
#endif /* VFTRACERS */
};

/* Doubly linked list for MCtracer Grid */
struct TracerListS{
    int count;  /*!< number of tracers in list */
    int i,j,k;  /*!< current location in grid */
    TracerS *Head;	/*!< pointer to first tracer in grid's list */
    TracerS *Tail;  /*!< pointer to last tracer in grid's list */
#ifdef MCTRACERS
    double Rmass;     /*!< reduced mass of grid cell */
    TracerS *currTail; /*!< pointer to last tracer in grid's list (not added in current loop */
#endif /* MCTRACERS */
    };

/* Data structure for storing tracer properties properties */
struct TracerPropS{
    double id;
    double d_init; /*!< initial density of underlying fluid */
    int i_init,j_init,k_init; /*!< initial coordinates of tracer */
#ifdef STAR_PARTICLE
    int star_id; /*!< id of associated star particle */
#endif /* STAR_PARTICLE */
};

/* Top hat algorithm structure */
struct TracertophatS{
    int count;  /*!< number of tracers in list */
};

#ifdef MPI_PARALLEL
/* Structure for passing list information via MPI */
struct TracerList_MPI{
    int count;  /*!< number of tracers in list */
};

/* Structure for passing tracer information via MPI */
struct Tracer_MPI{
    double id;         /*!< unique tracer identifier */
    double d_init; /*!< initial density of underlying fluid */
    int i_init;
    int j_init;
    int k_init;
    int star_id; /*!< id of associated star particle */
    double x1, x2, x3;
};
#endif /* MPI */

#endif /* TRACERS */
#endif /* TRACERS_H */
