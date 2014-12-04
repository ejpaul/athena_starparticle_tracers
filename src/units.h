#ifndef RADIATION_PROTOTYPES_H
#define RADIATION_PROTOTYPES_H 
#ifdef RADIATION
#include "../copyright.h"
/*==============================================================================
 * FILE: units.h
 *
 * PURPOSE: Definitions of various physical units, code units, etc.
 *============================================================================*/
#include "athena.h"
#include "defs.h"

typedef struct Unit_S{
  Real cm;
  Real g;
  Real s;
  Real Lcode;
  Real Mcode;
  Real Tcode;
  Real Dcode;
  Real dyne;
  Real erg;
  Real G;
  Real Msun;
  Real Lsun;
  Real Myr;
  Real pc;
  Real kpc;
  Real kms;
  Real mH;
  Real aR;
  Real kB;
  Real c;
}UnitS;

#endif /* RADIATION */
#endif /* RADIATION_PROTOTYPES_H */
