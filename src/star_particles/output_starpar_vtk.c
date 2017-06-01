#include "../copyright.h"
/*==============================================================================
 * FILE: output_starpar_vtk.c
 *
 * PURPOSE: Function to write star particle data in VTK "legacy" format.
 *
 * CONTAINS PUBLIC FUNCTIONS: 
 *   output_starpar_vtk() - writes VTK file containing star particle data.
 *============================================================================*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../defs.h"
#include "../athena.h"
#include "../prototypes.h"
#include "../globals.h"

#ifdef STAR_PARTICLE

/*=========================== PUBLIC FUNCTIONS ===============================*/
/*----------------------------------------------------------------------------*/
/* output_starpar_vtk:   */

void output_starpar_vtk(MeshS *pM, OutputS *pOut)
{
  GridS *pGrid = pM->Domain[0][0].Grid;
  FILE *pfile;
  char *fname;
  StarParListS *pGstars = NULL;
  StarParS *pStar = NULL;
  int i, ierr, nstars = 0;
  int *idata = NULL;
  float *data = NULL;
  float time = (float)pM->time;
  int step = pM->nstep;
  int big_end = ath_big_endian();

  if (myID_Comm_world == 0) {
    /* construct output filename.  pOut->id will either be name of variable,
     * if 'id=...' was included in <ouput> block, or 'outN' where N is number of
     * <output> block.  */
    sprintf(pOut->id,"starpar");
    fname = ath_fname(NULL, pM->outfilename, NULL, NULL, num_digit,
                      pOut->num, pOut->id, "vtk");
    if (fname == NULL)
      ath_error("[output_starpar_vtk]: Error constructing filename\n");
    
    /* open output file */
    pfile = fopen(fname,"w");
    if (pfile == NULL)
      ath_error("[output_starpar_vtk]: Unable to open vtk file %s\n",fname);
    
    /* Count the number of star particles, compute global MIN/MAX mass */
    pGstars = pGrid->Gstars;
    pOut->gmin = (pGstars) ? HUGE_NUMBER : 0.0;
    pOut->gmax = 0.0;
    nstars = 0;
    while (pGstars) {
      pStar = &(pGstars->starpar);
      pOut->gmin = MIN(pOut->gmin,pStar->m);
      pOut->gmax = MAX(pOut->gmax,pStar->m);
      nstars++;
      pGstars = pGstars->next;
    }
    
    /* Allocate memory for temporary array of floats */
    data = (float *)calloc_1d_array(3*nstars,sizeof(float));
    if (data == NULL) {
      ath_error("[output_starpar_vtk]: Error calling calloc_1d_array!\n");
      return;
    }
    
    /* Allocate memory for temporary array of ints */
    idata = (int *)calloc_1d_array(2*nstars,sizeof(int));
    if (idata == NULL) {
      ath_error("[output_starpar_vtk]: Error calling calloc_1d_array!\n");
      return;
    }
    
    /* There are five basic parts to the VTK "legacy" file format.  */
    /*  1. Write file version and identifier */
    fprintf(pfile,"# vtk DataFile Version 2.0\n");
    
    /*  2. Header */
    fprintf(pfile,"STAR PARTICLES at time= %e\n", pGrid->time);
    
    /*  3. File format */
    fprintf(pfile,"BINARY\n");
    
    /*  4. Dataset structure (star particle positions) */
    fprintf(pfile,"DATASET UNSTRUCTURED_GRID\n");

    fprintf(pfile,"POINTS %d float \n", nstars);
      pGstars = pGrid->Gstars;
      i = 0;
      while (pGstars) {
          pStar = &(pGstars->starpar);
          data[3*i  ] = (float)pStar->x1;
          data[3*i+1] = (float)pStar->x2;
          data[3*i+2] = (float)pStar->x3;
          i++;
          pGstars = pGstars->next;
      }
      if (!big_end) ath_bswap(data,sizeof(float),3*nstars);
      if (nstars) fwrite(data,sizeof(float),3*nstars,pfile);
    
    /* Treat all star particles as VTK_VERTEX (type=1) data structures.
     * There is 1 "cell" of type 1 containing each vertex, and it
     * takes 2*nstars numbers to describe all the cells. */
    fprintf(pfile,"CELLS %d %d \n", nstars,2*nstars);
    for (i=0; i<nstars; i++) {
      idata[2*i  ] = 1;
      idata[2*i+1] = i;
    }
    if (!big_end) ath_bswap(idata,sizeof(int),2*nstars);
    if (nstars) fwrite(idata,sizeof(int),2*nstars,pfile);
    fprintf(pfile,"\n");
    
    fprintf(pfile,"CELL_TYPES %d \n", nstars);
    for (i=0; i<nstars; i++) {
      idata[i] = 1;
    }
    if (!big_end) ath_bswap(idata,sizeof(int),nstars);
    if (nstars) fwrite(idata,sizeof(int),nstars,pfile);
    fprintf(pfile,"\n");

    /* Write time and step */
    fprintf(pfile,"FIELD FieldData 2\n");
    fprintf(pfile,"TIME 1 1 float \n");
    if (!big_end) ath_bswap(&time,sizeof(float),1);
    fwrite(&time,sizeof(float),1,pfile);
    fprintf(pfile,"\n");
    fprintf(pfile,"CYCLE 1 1 int \n");
    if (!big_end) ath_bswap(&step,sizeof(int),1);
    fwrite(&step,sizeof(int),1,pfile);
    fprintf(pfile,"\n");

    /* Write data  */
    fprintf(pfile,"CELL_DATA %d\n",nstars);
    fprintf(pfile,"POINT_DATA %d\n",nstars);

      
    /* Write star particle IDs */
    fprintf(pfile,"SCALARS star_particle_id int\n");
    fprintf(pfile,"LOOKUP_TABLE default\n");
    pGstars = pGrid->Gstars;
    i = 0;
    while (pGstars) {
      pStar = &(pGstars->starpar);
      idata[i] = pStar->id;
      i++;
      pGstars = pGstars->next;
    }
    if (!big_end) ath_bswap(idata,sizeof(int),nstars);
    if (nstars) fwrite(idata,sizeof(int),nstars,pfile);
    fprintf(pfile,"\n");
    
    /* Write star particle masses */
    fprintf(pfile,"SCALARS star_particle_mass float\n");
    fprintf(pfile,"LOOKUP_TABLE default\n");
    pGstars = pGrid->Gstars;
    i = 0;
    while (pGstars) {
      pStar = &(pGstars->starpar);
      data[i] = (float)pStar->m;
      i++;
      pGstars = pGstars->next;
    }
    if (!big_end) ath_bswap(data,sizeof(float),nstars);
    if (nstars) fwrite(data,sizeof(float),nstars,pfile);
    fprintf(pfile,"\n");
    
    /* Write star particle velocities */
    fprintf(pfile,"VECTORS star_particle_velocity float\n");
    pGstars = pGrid->Gstars;
    i = 0;
    while (pGstars) {
      pStar = &(pGstars->starpar);
      data[3*i  ] = (float)pStar->v1;
      data[3*i+1] = (float)pStar->v2;
      data[3*i+2] = (float)pStar->v3;
      i++;
      pGstars = pGstars->next;
    }
    if (!big_end) ath_bswap(data,sizeof(float),3*nstars);
    if (nstars) fwrite(data,sizeof(float),3*nstars,pfile);
    fprintf(pfile,"\n");
    
    fclose(pfile);
    free_1d_array(data);
    free_1d_array(idata);
  }

  return;
}
#endif /* STAR_PARTICLE */
