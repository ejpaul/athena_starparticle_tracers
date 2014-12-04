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
typedef struct TracerHistS TracerHistS;
typedef struct TracerS TracerS;
typedef struct TracertophatS TracertophatS;
typedef struct TracerList_MPI TracerList_MPI;
typedef struct Tracer_MPI Tracer_MPI;

/* Monte Carlo tracer particle structure */
struct TracerS{
    TracerListS* newList; /*!< pointer to new list, if !=NULL flags for movement */
    TracerS* Next;  /*!< pointer to next tracer in list */
    TracerS* Prev;  /*!< pointer to previous tracer in list */
    TracerHistS* hist; /*!< data structure including history of underlying fluid */
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

/* Data structure for tracing history underlying fluid properties */
struct TracerHistS{
    double id;
    double d_init; /*!< initial density of underlying fluid */
    int i_init,j_init,k_init; /*!< initial coordinates of tracer */
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
    double i_init;
    double j_init;
    double k_init;
//    int i_init,j_init,k_init; /*!< initial coordinates of tracer */
};
#endif /* MPI */

#endif /* TRACERS */
#endif /* TRACERS_H */
