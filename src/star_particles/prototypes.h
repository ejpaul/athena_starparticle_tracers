#ifndef STAR_PARTICLE_PROTOTYPES_H
#define STAR_PARTICLE_PROTOTYPES_H
#ifdef STAR_PARTICLE
#include "../copyright.h"
/*==============================================================================
 * FILE: prototypes.h
 *
 * PURPOSE: Prototypes for all public functions in the /src/star_particle directory.
 *============================================================================*/
#include <stdio.h>
#include <stdarg.h>
#include "../athena.h"
#include "../defs.h"

#include "../config.h"

#ifdef SELF_GRAVITY
/* assign_starparticles.c */
void assign_starparticles_3d(DomainS *pD, ath_fft_data *work);
#ifdef SHEARING_BOX
void assign_starparticles_3d_shear(DomainS *pD, Real ***RoolDen);
#endif
#endif

/* create_starparticles.c */
void create_starparticles(DomainS *pD);

/* integrate_starparticles.c */
void integrate_starparticles(DomainS *pD);
void kick_correct_starparticles(DomainS *pD);

/* modify_ghost_region_starparticles.c */
void modify_ghost_region_starparticles(DomainS *pD);
void set_ghost_region(GridS *pG, StarParS *pStar);

/* synchro_starparticles.c */
void synchro_starparticles(DomainS *pD);

/* update_starparticles.c */
void update_starparticles(GridS *pG, Cons1DS ***x1Flux, Cons1DS ***x2Flux, Cons1DS ***x3Flux);
void age_starparticles(GridS *pG);

/* utils_starparticles.c */
void starpar_push_local(GridS *pG, StarParListS *pList);
void starpar_push_global(GridS *pG, StarParListS *pList);
void starpar_push_global_fda(GridS *pG, StarParListS *pList);
void starpar_destruct_local(GridS *pG);
void starpar_destruct_global(GridS *pG);
void starpar_destruct_global_fda(GridS *pG);
int  starpar_ingrid(DomainS *pD, StarParS *pStar);
void starpar_printlist(const int level, const GridS *pG);

/* dump_starpar_history.c */
/*! \fn Real (*StarParFun_t)(const GridS *pG, const StarParS *pStar)
 *  \brief user define Star Particle history function */
typedef Real (*StarParFun_t)(const GridS *pG, const StarParS *pStar);

void dump_starpar_history(MeshS *pM, OutputS *pOut);
void dump_starpar_history_enroll(const StarParFun_t pfun, const char *label);

/* output_starpar_vtk.c */
void output_starpar_vtk(MeshS *pM, OutputS *pOut);

/* feedback_starparticles.c */
void feedback_starparticles_init(MeshS *pM);
void feedback_starparticles_destruct(void);
void feedback_starparticles_assign(DomainS *pD);
void feedback_starparticles(DomainS *pD);

#endif /* STAR_PARTICLE */
#endif /* STAR_PARTICLE_PROTOTYPES_H */
