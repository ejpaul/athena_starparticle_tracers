#include "../copyright.h"
/*==============================================================================
 * FILE: output_tracer_vtk.c
 *
 * PURPOSE: Function to write tracer data in VTK "legacy" format.
 *
 * CONTAINS PUBLIC FUNCTIONS: 
 *   output_tracer_vtk() - writes VTK file containing tracer data.
 *============================================================================*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../defs.h"
#include "../athena.h"
#include "../prototypes.h"
#include "../globals.h"
#include "prototypes.h"

#if defined(MCTRACERS) || defined(VFTRACERS)
#ifdef STAR_PARTICLE

/*=========================== PUBLIC FUNCTIONS ===============================*/
/*----------------------------------------------------------------------------*/
/* output_tracer_vtk:   */

void output_tracer_vtk(MeshS *pM, OutputS *pOut)
{
  GridS *pGrid = pM->Domain[0][0].Grid;
  FILE *pfile;
  char *fname;
  int ierr;
  int *idata = NULL;
  float *data = NULL;
  float time = (float)pM->time;
  int step = pM->nstep;
  int big_end = ath_big_endian();
    int ntracers;
    int ks = pGrid->ks;
    int js = pGrid->js;
    int is = pGrid->is;
    int ke = pGrid->ke;
    int je = pGrid->je;
    int ie = pGrid->ie;
    int i, j, k;
    int count, m;
    Real x1, x2, x3;
    TracerListS *list;
    TracerS *tracer;

  if (myID_Comm_world == 0) {
    /* construct output filename.  pOut->id will either be name of variable,
     * if 'id=...' was included in <ouput> block, or 'outN' where N is number of
     * <output> block.  */
    sprintf(pOut->id,"tracer");
    fname = ath_fname(NULL, pM->outfilename, NULL, NULL, num_digit,
                      pOut->num, pOut->id, "vtk");
    if (fname == NULL)
      ath_error("[output_tracer_vtk]: Error constructing filename\n");
    
    /* open output file */
    pfile = fopen(fname,"w");
    if (pfile == NULL)
      ath_error("[output_tracer_vtk]: Unable to open vtk file %s\n",fname);
    
    /* Count the number of tracers within sink particles */
    ntracers = 0;
      for (k=ks; k<=ke; k++) {
          for (j=js; j<=je; j++) {
              for (i=is; i<=ie; i++) {
                  count = 0;
                  list = &((pGrid->GridLists)[k][j][i]);
                  tracer = list->Head;
                  while(tracer) {
                      if (tracer->prop->star_id != -1) {
                          count++;
                      }
                      tracer = tracer->Next;
                  }
                  ntracers += count;
              }
          }
      }
    
    /* Allocate memory for temporary array of floats */
    data = (float *)calloc_1d_array(3*ntracers,sizeof(float));
    if (data == NULL) {
      ath_error("[output_tracer_vtk]: Error calling calloc_1d_array!\n");
      return;
    }
    
    /* Allocate memory for temporary array of ints */
    idata = (int *)calloc_1d_array(2*ntracers,sizeof(int));
    if (idata == NULL) {
      ath_error("[output_tracer_vtk]: Error calling calloc_1d_array!\n");
      return;
    }
    
    /* There are five basic parts to the VTK "legacy" file format.  */
    /*  1. Write file version and identifier */
    fprintf(pfile,"# vtk DataFile Version 2.0\n");
    
    /*  2. Header */
    fprintf(pfile,"TRACERS at time= %e\n", pGrid->time);
    
    /*  3. File format */
    fprintf(pfile,"BINARY\n");
    
    /*  4. Dataset structure (star particle positions) */
    fprintf(pfile,"DATASET UNSTRUCTURED_GRID\n");

    fprintf(pfile,"NTRACERS %d \n", ntracers);
    
    /* Treat all tracers as VTK_VERTEX (type=1) data structures.
     * There is 1 "cell" of type 1 containing each vertex, and it
     * takes 2*ntracers numbers to describe all the cells. */
    fprintf(pfile,"CELLS %d %d \n", ntracers,2*ntracers);
    for (i=0; i<ntracers; i++) {
      idata[2*i  ] = 1;
      idata[2*i+1] = i;
    }
    if (!big_end) ath_bswap(idata,sizeof(int),2*ntracers);
    if (ntracers) fwrite(idata,sizeof(int),2*ntracers,pfile);
    fprintf(pfile,"\n");
    
    fprintf(pfile,"CELL_TYPES %d \n", ntracers);
    for (i=0; i<ntracers; i++) {
      idata[i] = 1;
    }
    if (!big_end) ath_bswap(idata,sizeof(int),ntracers);
    if (ntracers) fwrite(idata,sizeof(int),ntracers,pfile);
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
    fprintf(pfile,"CELL_DATA %d\n",ntracers);
    fprintf(pfile,"POINT_DATA %d\n",ntracers);

    /* Write tracer IDs */
    fprintf(pfile,"SCALARS tracer_id float\n");
    fprintf(pfile,"LOOKUP_TABLE default\n");
      m = 0;
      for (k=ks; k<=ke; k++) {
          for (j=js; j<=je; j++) {
              for (i=is; i<=ie; i++) {
                  list = &((pGrid->GridLists)[k][j][i]);
                  tracer = list->Head;
                  while(tracer) {
                      if (tracer->prop->star_id != -1) {
                          data[m] = (float)tracer->prop->id;
                          m++;
                      }
                      tracer = tracer->Next;
                  }
              }
          }
      }
    if (!big_end) ath_bswap(data,sizeof(float),ntracers);
    if (ntracers) fwrite(data,sizeof(float),ntracers,pfile);
    fprintf(pfile,"\n");
      
      /* Write associated star particle IDs */
      fprintf(pfile,"SCALARS star_id int\n");
      fprintf(pfile,"LOOKUP_TABLE default\n");
      m = 0;
      for (k=ks; k<=ke; k++) {
          for (j=js; j<=je; j++) {
              for (i=is; i<=ie; i++) {
                  list = &((pGrid->GridLists)[k][j][i]);
                  tracer = list->Head;
                  while(tracer) {
                      if (tracer->prop->star_id != -1) {
                          idata[m] = (int)tracer->prop->star_id;
                          m++;
                      }
                      tracer = tracer->Next;
                  }
              }
          }
      }
      if (!big_end) ath_bswap(idata,sizeof(int),ntracers);
      if (ntracers) fwrite(idata,sizeof(int),ntracers,pfile);
      fprintf(pfile,"\n");
    
    /* Write tracer initial density */
    fprintf(pfile,"SCALARS d_init float\n");
    fprintf(pfile,"LOOKUP_TABLE default\n");
      m = 0;
      for (k=ks; k<=ke; k++) {
          for (j=js; j<=je; j++) {
              for (i=is; i<=ie; i++) {
                  list = &((pGrid->GridLists)[k][j][i]);
                  tracer = list->Head;
                  while(tracer) {
                      if (tracer->prop->star_id != -1) {
                          data[m] = (float)tracer->prop->d_init;
                          m++;
                      }
                      tracer = tracer->Next;
                  }
              }
          }
      }

    if (!big_end) ath_bswap(data,sizeof(float),ntracers);
    if (ntracers) fwrite(data,sizeof(float),ntracers, pfile);
    fprintf(pfile,"\n");
      
      /* Write tracer initial i */
      fprintf(pfile,"SCALARS i_init float\n");
      fprintf(pfile,"LOOKUP_TABLE default\n");
      m = 0;
      for (k=ks; k<=ke; k++) {
          for (j=js; j<=je; j++) {
              for (i=is; i<=ie; i++) {
                  list = &((pGrid->GridLists)[k][j][i]);
                  tracer = list->Head;
                  while(tracer) {
                      if (tracer->prop->star_id != -1) {
                          data[m] = (float)tracer->prop->i_init;
                          m++;
                      }
                      tracer = tracer->Next;
                  }
              }
          }
      }
      
      if (!big_end) ath_bswap(data,sizeof(float),ntracers);
      if (ntracers) fwrite(data,sizeof(float),ntracers, pfile);
      fprintf(pfile,"\n");
      
      /* Write tracer initial j */
      fprintf(pfile,"SCALARS j_init float\n");
      fprintf(pfile,"LOOKUP_TABLE default\n");
      m = 0;
      for (k=ks; k<=ke; k++) {
          for (j=js; j<=je; j++) {
              for (i=is; i<=ie; i++) {
                  list = &((pGrid->GridLists)[k][j][i]);
                  tracer = list->Head;
                  while(tracer) {
                      if (tracer->prop->star_id != -1) {
                          data[m] = (float)tracer->prop->j_init;
                          m++;
                      }
                      tracer = tracer->Next;
                  }
              }
          }
      }
      
      if (!big_end) ath_bswap(data,sizeof(float),ntracers);
      if (ntracers) fwrite(data,sizeof(float),ntracers, pfile);
      fprintf(pfile,"\n");
      
      /* Write tracer initial k */
      fprintf(pfile,"SCALARS k_init float\n");
      fprintf(pfile,"LOOKUP_TABLE default\n");
      m = 0;
      for (k=ks; k<=ke; k++) {
          for (j=js; j<=je; j++) {
              for (i=is; i<=ie; i++) {
                  list = &((pGrid->GridLists)[k][j][i]);
                  tracer = list->Head;
                  while(tracer) {
                      if (tracer->prop->star_id != -1) {
                          data[m] = (float)tracer->prop->j_init;
                          m++;
                      }
                      tracer = tracer->Next;
                  }
              }
          }
      }
      
      if (!big_end) ath_bswap(data,sizeof(float),ntracers);
      if (ntracers) fwrite(data,sizeof(float),ntracers, pfile);
      fprintf(pfile,"\n");

     /* Write tracer positions */
    fprintf(pfile,"VECTORS tracer_position float\n");
    m = 0;
      for (k=ks; k<=ke; k++) {
          for (j=js; j<=je; j++) {
              for (i=is; i<=ie; i++) {
                  list = &((pGrid->GridLists)[k][j][i]);
                  tracer = list->Head;
                  while(tracer) {
                      if (tracer->prop->star_id != -1) {
#ifdef VFTRACERS
                          data[3*m  ] = (float)tracer->x1;
                          data[3*m+1] = (float)tracer->x2;
                          data[3*m+2] = (float)tracer->x3;
#endif /* VFTRACERS */
    /* In case of MCTRACERS use cell-centered position */
#ifdef MCTRACERS
                          cc_pos(pGrid, i, j, k, &x1, &x2, &x3);
                          data[3*m  ] = (float)x1;
                          data[3*m+1] = (float)x2;
                          data[3*m+2] = (float)x3;
#endif /* MCTRACERS */
                          m++;
                      }
                      tracer = tracer->Next;
                  }
              }
          }
      }
      
    if (!big_end) ath_bswap(data,sizeof(float),3*ntracers);
    if (ntracers) fwrite(data,sizeof(float),3*ntracers,pfile);
    fprintf(pfile,"\n");
    
    fclose(pfile);
    free_1d_array(data);
    free_1d_array(idata);
  }

  return;
}
#endif /* STAR_PARTICLE */
#endif /* TRACERS */
