#include "../copyright.h"
/*=============================================================================
 * FILE: utils_starpar.c
 *
 * PURPOSE: Contains most of the utilities for star particles:
 *
 * CONTAINS PUBLIC FUNCTIONS:
 *   starpar_push_local()
 *   starpar_push_global()
 *   starpar_destruct_local()
 *   starpar_destruct_global()
 *
 * History:
 *   Written by Aaron Skinner, May 2010
 *   Modified by Hao Gong, Apr. 2011
 *   Modified by Aaron Skinner, Apr. 2013
 *============================================================================*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../defs.h"
#include "../athena.h"
#include "../prototypes.h"
#include "../globals.h"

#ifdef STAR_PARTICLE

/*-----------------------------PUSH/POP---------------------------------------*/
/* Push a new star particle pointed to by pList onto the head of the source list
 * of the Grid pointed to by pG.  Note that this function does not create the
 * new source struct. */
void starpar_push_local(GridS *pG, StarParListS *pList)
{
  pList->next = pG->Lstars;
  pG->Lstars = pList;
  pG->nLstars++;
  return;
}

void starpar_push_global(GridS *pG, StarParListS *pList)
{
  pList->next = pG->Gstars;
  pG->Gstars = pList;
  pG->nGstars++;
  return;
}

void starpar_push_global_fda(GridS *pG, StarParListS *pList)
{
  pList->next = pG->Gstars_fda;
  pG->Gstars_fda = pList;
  return;
}

/*-----------------------------DESTRUCTOR-------------------------------------*/

void starpar_destruct_local(GridS *pG)
{
  StarParListS *pList = pG->Lstars;

  while (pList) {
    pList = pList->next;
    free_1d_array(pG->Lstars);
    pG->Lstars = pList;
  }

  pG->nLstars = 0;
  pG->Lstars = NULL;
  return;
}

void starpar_destruct_global(GridS *pG)
{
  StarParListS *pList = pG->Gstars;

  while (pList) {
    pList = pList->next;
    free_1d_array(pG->Gstars);
    pG->Gstars = pList;
  }
  
  pG->nGstars = 0;
  pG->Gstars = NULL;
  return;
}

void starpar_destruct_global_fda(GridS *pG)
{
  StarParListS *pList = pG->Gstars_fda;

  while (pList) {
    pList = pList->next;
    free_1d_array(pG->Gstars_fda);
    pG->Gstars_fda = pList;
  }
  
  pG->Gstars_fda = NULL;
  return;
}

int starpar_ingrid(DomainS *pD, StarParS *pStar){
  GridS *pG = pD->Grid;
  
  if (pStar->x1 == pD->MaxX[0])
    if ((pD->bcflag_ix1 == 4) && (pD->bcflag_ox1 == 4)) pStar->x1 = pD->MinX[0];
    else return 0;
  if (pStar->x2 == pD->MaxX[1])
    if ((pD->bcflag_ix2 == 4) && (pD->bcflag_ox2 == 4)) pStar->x2 = pD->MinX[1];
    else return 0;
  if (pStar->x3 == pD->MaxX[2])
    if ((pD->bcflag_ix3 == 4) && (pD->bcflag_ox3 == 4)) pStar->x3 = pD->MinX[2];
    else return 0;
  
  return ((pStar->x1 >= pG->MinX[0]) && (pStar->x1 < pG->MaxX[0])
       && (pStar->x2 >= pG->MinX[1]) && (pStar->x2 < pG->MaxX[1])
       && (pStar->x3 >= pG->MinX[2]) && (pStar->x3 < pG->MaxX[2]));
}

/*-----------------------------OUTPUT/DEBUG-----------------------------------*/

void starpar_printlist(const int level, const GridS *pG)
{
  StarParListS *pLstars = pG->Lstars;
  StarParListS *pGstars = pG->Gstars_fda;
  StarParS *pStar = NULL;
  
  if (pLstars)
    ath_perr(level, "\nPrinting list of local star particles...\n");
  else
    ath_perr(level, "\nNo local star particles!\n");
  
  while (pLstars) {
    pStar = &(pLstars->starpar);
    ath_perr(level, "\nLOCAL STAR PARTICLE %d:\n",pStar->id);
    ath_perr(level, "\tMass:  m = %f\n", pStar->m);
    ath_perr(level, "\tAccretion Rate:  mdot = %f\n", pStar->mdot);
    ath_perr(level, "\tLocation (%f,%f,%f):\n",pStar->x1,pStar->x2,pStar->x3);
    ath_perr(level, "\tVelocity (%f,%f,%f):\n",pStar->v1,pStar->v2,pStar->v3);
    ath_perr(level, "\tAge:  %f\n", pStar->age);
    pLstars = pLstars->next;
  }
  
  if (pGstars)
    ath_perr(level, "\nPrinting list of global star particles...\n");
  else
    ath_perr(level, "\nNo global star particles!\n");
  
  while (pGstars) {
    pStar = &(pGstars->starpar);
    ath_perr(level, "\nGLOBAL STAR PARTICLE %d:\n",pStar->id);
    ath_perr(level, "\tMass:  m = %f\n", pStar->m);
    ath_perr(level, "\tAccretion Rate:  mdot = %f\n", pStar->mdot);
    ath_perr(level, "\tLocation (%f,%f,%f):\n",pStar->x1,pStar->x2,pStar->x3);
    ath_perr(level, "\tVelocity (%f,%f,%f):\n",pStar->v1,pStar->v2,pStar->v3);
    ath_perr(level, "\tAge:  %f\n", pStar->age);
    pGstars = pGstars->next;
  }
  
  return;
}

#endif /*STAR_PARTICLE */
