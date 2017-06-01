#include "../copyright.h"
/*===========================================================================
 * \file bvals_mc.c
 * \brief
 *
 * PURPOSE:
 *============================================================================*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "../defs.h"
#include "../athena.h"
#include "../prototypes.h"
#include "prototypes.h"
#include "tracers.h"
#include "../globals.h"
#include <float.h>

#if defined(MCTRACERS) || defined(VFTRACERS)  /* endif at the end of the file */

/* boundary condition function pointers. local to this function  */
static VGFun_t apply_mc_ix1 = NULL, apply_mc_ox1 = NULL;
static VGFun_t apply_mc_ix2 = NULL, apply_mc_ox2 = NULL;
static VGFun_t apply_mc_ix3 = NULL, apply_mc_ox3 = NULL;

#ifdef MPI_PARALLEL
/* MPI send and receive buffers */
static int *send_buf = NULL, *recv_buf = NULL;
static Tracer_MPI *send_buf_mc0x1 = NULL, *send_buf_mc1x1 = NULL;
static Tracer_MPI *recv_buf_mc0x1 = NULL, *recv_buf_mc1x1 = NULL;
static Tracer_MPI *send_buf_mc0x2 = NULL, *send_buf_mc1x2 = NULL;
static Tracer_MPI *recv_buf_mc0x2 = NULL, *recv_buf_mc1x2 = NULL;
static Tracer_MPI *send_buf_mc0x3 = NULL, *send_buf_mc1x3 = NULL;
static Tracer_MPI *recv_buf_mc0x3 = NULL, *recv_buf_mc1x3 = NULL;

static TracerList_MPI *send_buf_list0x1 = NULL, *send_buf_list1x1 = NULL;
static TracerList_MPI *recv_buf_list0x1 = NULL, *recv_buf_list1x1 = NULL;
static TracerList_MPI *send_buf_list0x2 = NULL, *send_buf_list1x2 = NULL;
static TracerList_MPI *recv_buf_list0x2 = NULL, *recv_buf_list1x2 = NULL;
static TracerList_MPI *send_buf_list0x3 = NULL, *send_buf_list1x3 = NULL;
static TracerList_MPI *recv_buf_list0x3 = NULL, *recv_buf_list1x3 = NULL;

static MPI_Request *recv_rq = NULL, *send_rq = NULL;
#endif /* MPI_PARALLEL */

/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 *   reflect_???()  - reflecting BCs at boundary ???
 *   outflow_???()  - outflow BCs at boundary ???
 *   periodic_???() - periodic BCs at boundary ???
 *   conduct_???()  - conducting BCs at boundary ???
 *   pack_???()     - pack data for MPI non-blocking send at ??? boundary
 *   unpack_???()   - unpack data for MPI non-blocking receive at ??? boundary
 *============================================================================*/

static void reflect_ix1_mc(GridS *pG);
static void reflect_ox1_mc(GridS *pG);
static void reflect_ix2_mc(GridS *pG);
static void reflect_ox2_mc(GridS *pG);
static void reflect_ix3_mc(GridS *pG);
static void reflect_ox3_mc(GridS *pG);

static void outflow_ix1_mc(GridS *pG);
static void outflow_ox1_mc(GridS *pG);
static void outflow_ix2_mc(GridS *pG);
static void outflow_ox2_mc(GridS *pG);
static void outflow_ix3_mc(GridS *pG);
static void outflow_ox3_mc(GridS *pG);

static void periodic_ix1_mc(GridS *pG);
static void periodic_ox1_mc(GridS *pG);
static void periodic_ix2_mc(GridS *pG);
static void periodic_ox2_mc(GridS *pG);
static void periodic_ix3_mc(GridS *pG);
static void periodic_ox3_mc(GridS *pG);

#ifdef MPI_PARALLEL
static int pack_ix1_mc(GridS *pG);
static int pack_ox1_mc(GridS *pG);
static int pack_ix2_mc(GridS *pG);
static int pack_ox2_mc(GridS *pG);
static int pack_ix3_mc(GridS *pG);
static int pack_ox3_mc(GridS *pG);

static void pack_ix1_tracers(GridS *pG);
static void pack_ox1_tracers(GridS *pG);
static void pack_ix2_tracers(GridS *pG);
static void pack_ox2_tracers(GridS *pG);
static void pack_ix3_tracers(GridS *pG);
static void pack_ox3_tracers(GridS *pG);

static void unpack_ix1_mc(GridS *pG);
static void unpack_ox1_mc(GridS *pG);
static void unpack_ix2_mc(GridS *pG);
static void unpack_ox2_mc(GridS *pG);
static void unpack_ix3_mc(GridS *pG);
static void unpack_ox3_mc(GridS *pG);
#endif /* MPI_PARALLEL */

/*=========================== PUBLIC FUNCTIONS ===============================*/

void bvals_tracer(DomainS *pD) {
    GridS *pGrid = (pD->Grid);
#ifdef MPI_PARALLEL
    int err, mIndex;
    int cnt;
    
    if((recv_buf = (int*)calloc_1d_array(2,sizeof(int))) == NULL)
        ath_error("[bvals_init]: Failed to allocate recv buffer\n");
    if((send_buf = (int*)calloc_1d_array(2,sizeof(int))) == NULL)
        ath_error("[bvals_init]: Failed to allocate send buffer\n");
    if((recv_rq = (MPI_Request*) calloc_1d_array(2,sizeof(MPI_Request))) == NULL)
        ath_error("[bvals_init]: Failed to allocate recv MPI_Request array\n");
    if((send_rq = (MPI_Request*) calloc_1d_array(2,sizeof(MPI_Request))) == NULL)
        ath_error("[bvals_init]: Failed to allocate send MPI_Request array\n");
    
#endif /* MPI_PARALLEL */
    
    /*--- Step 1. ------------------------------------------------------------------
     * Boundary Conditions in x1-direction */
    if (pGrid->Nx[0] > 1){
                
#ifdef MPI_PARALLEL
        cnt = nghost*(pGrid->Nx[1])*(pGrid->Nx[2]); // Number of cells
        
        if((send_buf_list0x1 = (TracerList_MPI*)calloc_1d_array(cnt,sizeof(TracerList_MPI))) == NULL)
            ath_error("[bvals_init]: Failed to allocate send buffer\n");
        if((send_buf_list1x1 = (TracerList_MPI*)calloc_1d_array(cnt,sizeof(TracerList_MPI))) == NULL)
            ath_error("[bvals_init]: Failed to allocate send buffer\n");
        if((recv_buf_list0x1 = (TracerList_MPI*)calloc_1d_array(cnt,sizeof(TracerList_MPI))) == NULL)
            ath_error("[bvals_init]: Failed to allocate send buffer\n");
        if((recv_buf_list1x1 = (TracerList_MPI*)calloc_1d_array(cnt,sizeof(TracerList_MPI
                                                                       ))) == NULL)
            ath_error("[bvals_init]: Failed to allocate send buffer\n");
        
        /* MPI blocks to both left and right */
        if (pGrid->rx1_id >= 0 && pGrid->lx1_id >= 0) {
            
            /* Post non-blocking receives for data size from L and R Grids */
            err = MPI_Irecv(&(recv_buf[0]), 1, MPI_INTEGER, pGrid->lx1_id, LtoR_tag, pD->Comm_Domain, &(recv_rq[0]));
            err = MPI_Irecv(&(recv_buf[1]), 1, MPI_INTEGER, pGrid->rx1_id, RtoL_tag, pD->Comm_Domain, &(recv_rq[1]));
            
            /* pack and send data size to L and R */
            send_buf[0] = pack_ix1_mc(pGrid);

            /* calloc tracer buffer for send to l */
            if((send_buf_mc0x1 = (Tracer_MPI *)calloc_1d_array(send_buf[0], sizeof(Tracer_MPI))) == NULL)
                ath_error("[bvals_init]: Failed to allocate recv buffer\n");
            
            err = MPI_Isend(&(send_buf[0]), 1, MPI_INTEGER, pGrid->lx1_id, RtoL_tag, pD->Comm_Domain, &(send_rq[0]));
            
            send_buf[1] = pack_ox1_mc(pGrid);
            
             // calloc tracer buffer for send to r
            if((send_buf_mc1x1 = (Tracer_MPI *)calloc_1d_array(send_buf[1], sizeof(Tracer_MPI))) == NULL)
                ath_error("[bvals_init]: Failed to allocate recv buffer\n");
            
            err = MPI_Isend(&(send_buf[1]), 1, MPI_INTEGER, pGrid->rx1_id, LtoR_tag, pD->Comm_Domain, &(send_rq[1]));

            /* check non-blocking sends have completed. */
            err = MPI_Waitall(2, send_rq, MPI_STATUS_IGNORE);
            
            /* check non-blocking data size receives */
            err = MPI_Waitany(2, recv_rq, &mIndex, MPI_STATUS_IGNORE);
            
            if (mIndex == 0) {
                /* calloc buffer for tracer receives from l */
                if((recv_buf_mc0x1 = (Tracer_MPI *)calloc_1d_array(recv_buf[0], sizeof(Tracer_MPI))) == NULL)
                    ath_error("[bvals_init]: Failed to allocate recv buffer\n");
            }
            if (mIndex == 1) {
                /* calloc buffer for tracer recevies from r */
                if((recv_buf_mc1x1 = (Tracer_MPI *)calloc_1d_array(recv_buf[1], sizeof(Tracer_MPI))) == NULL)
                    ath_error("[bvals_init]: Failed to allocate recv buffer\n");
            }
            err = MPI_Waitany(2, recv_rq, &mIndex,MPI_STATUS_IGNORE);
            if (mIndex == 0) {
                /* calloc buffer for tracer receives from l */
                if((recv_buf_mc0x1 = (Tracer_MPI *)calloc_1d_array(recv_buf[0], sizeof(Tracer_MPI))) == NULL)
                    ath_error("[bvals_init]: Failed to allocate recv buffer\n");
            }
            if (mIndex == 1) {
                /* calloc buffer for tracer recevies from r */
                if((recv_buf_mc1x1 = (Tracer_MPI *)calloc_1d_array(recv_buf[1], sizeof(Tracer_MPI))) == NULL)
                    ath_error("[bvals_init]: Failed to allocate recv buffer\n");
            }

//             Post non-blocking receives for tracers from L and R Grids 
            err = MPI_Irecv(&(recv_buf_mc0x1[0]),recv_buf[0],pD->TRACERTYPE,pGrid->lx1_id,LtoR_tag,
                             pD->Comm_Domain, &(recv_rq[0]));
            err = MPI_Irecv(&(recv_buf_mc1x1[0]),recv_buf[1],pD->TRACERTYPE,pGrid->rx1_id,RtoL_tag,
                             pD->Comm_Domain, &(recv_rq[1]));
            
            /* pack and send data L and R */
            pack_ix1_tracers(pGrid);
            
            err = MPI_Isend(&(send_buf_mc0x1[0]),send_buf[0],pD->TRACERTYPE,pGrid->lx1_id,RtoL_tag, pD->Comm_Domain, &(send_rq[0]));

            pack_ox1_tracers(pGrid);

            err = MPI_Isend(&(send_buf_mc1x1[0]),send_buf[1],pD->TRACERTYPE,pGrid->rx1_id,LtoR_tag, pD->Comm_Domain, &(send_rq[1]));

            /* check non-blocking sends have completed. */
            err = MPI_Waitall(2,send_rq, MPI_STATUS_IGNORE);

            /* check non-blocking receives and unpack data in any order. */
            err = MPI_Waitany(2,recv_rq,&mIndex,MPI_STATUS_IGNORE);
            err = MPI_Waitany(2,recv_rq,&mIndex,MPI_STATUS_IGNORE);
            
            /* Post non-blocking receives for list structures from L and R Grids */
            err = MPI_Irecv(&(recv_buf_list0x1[0]), cnt, pD->LISTTYPE, pGrid->lx1_id,LtoR_tag, pD->Comm_Domain, &(recv_rq[0]));
            err = MPI_Irecv(&(recv_buf_list1x1[0]), cnt, pD->LISTTYPE, pGrid->rx1_id,RtoL_tag, pD->Comm_Domain, &(recv_rq[1]));
            
            /* Send list structures L and R */
            err = MPI_Isend(&(send_buf_list0x1[0]), cnt, pD->LISTTYPE,pGrid->lx1_id,RtoL_tag, pD->Comm_Domain, &(send_rq[0]));
            err = MPI_Isend(&(send_buf_list1x1[0]), cnt, pD->LISTTYPE,pGrid->rx1_id,LtoR_tag, pD->Comm_Domain, &(send_rq[1]));
            
            /* check non-blocking sends have completed. */
            err = MPI_Waitall(2, send_rq, MPI_STATUS_IGNORE);
            
            /* check non-blocking receives and unpack data in any order. */
            err = MPI_Waitany(2, recv_rq,&mIndex,MPI_STATUS_IGNORE);
            if (mIndex == 0) unpack_ix1_mc(pGrid);
            if (mIndex == 1) unpack_ox1_mc(pGrid);
            err = MPI_Waitany(2, recv_rq,&mIndex,MPI_STATUS_IGNORE);
            if (mIndex == 0) unpack_ix1_mc(pGrid);
            if (mIndex == 1) unpack_ox1_mc(pGrid);
            
            if (send_buf_list0x1) free_1d_array(send_buf_list0x1);
            if (send_buf_list1x1) free_1d_array(send_buf_list1x1);
            if (recv_buf_list0x1) free_1d_array(recv_buf_list0x1);
            if (recv_buf_list1x1) free_1d_array(recv_buf_list1x1);
            if (recv_buf_mc0x1) free_1d_array(recv_buf_mc0x1);
            if (recv_buf_mc1x1) free_1d_array(recv_buf_mc1x1);
            if (send_buf_mc0x1) free_1d_array(send_buf_mc0x1);
            if (send_buf_mc1x1) free_1d_array(send_buf_mc1x1);
        
        }
    
        /* Physical boundary on left, MPI block on right */
        if (pGrid->rx1_id >= 0 && pGrid->lx1_id < 0) {
            
            /* Post non-blocking receives for data size from R Grids */
            err = MPI_Irecv(&(recv_buf[1]), 1, MPI_INTEGER, pGrid->rx1_id, RtoL_tag, pD->Comm_Domain, &(recv_rq[1]));
            
            /* pack and send data size to R */
            send_buf[1] = pack_ox1_mc(pGrid);
            
            /* calloc tracer buffer for send to r */
            if((send_buf_mc1x1 = (Tracer_MPI *)calloc_1d_array(send_buf[1], sizeof(Tracer_MPI))) == NULL)
                ath_error("[bvals_init]: Failed to allocate recv buffer\n");

            err = MPI_Isend(&(send_buf[1]), 1, MPI_INTEGER, pGrid->rx1_id, LtoR_tag, pD->Comm_Domain, &(send_rq[1]));
            
            /* set physical boundary */
            (*(apply_mc_ix1))(pGrid);
            
            /* check non-blocking sends have completed. */
            err = MPI_Wait(&(send_rq[1]), MPI_STATUS_IGNORE);
            
            /* check non-blocking data size received from R */
            err = MPI_Wait(&(recv_rq[1]), MPI_STATUS_IGNORE);

            /* calloc buffer for tracer recevies from r */
            if((recv_buf_mc1x1 = (Tracer_MPI *)calloc_1d_array(recv_buf[1], sizeof(Tracer_MPI))) == NULL)
                ath_error("[bvals_init]: Failed to allocate recv buffer\n");

            /* Post non-blocking receives for tracers from R Grids */
            err = MPI_Irecv(&(recv_buf_mc1x1[0]),recv_buf[1],pD->TRACERTYPE,pGrid->rx1_id,RtoL_tag, pD->Comm_Domain, &(recv_rq[1]));

            /* pack and send tracers R */
            pack_ox1_tracers(pGrid);

            err = MPI_Isend(&(send_buf_mc1x1[0]),send_buf[1],pD->TRACERTYPE,pGrid->rx1_id,LtoR_tag, pD->Comm_Domain, &(send_rq[1]));

            /* check non-blocking sends have completed. */
            err = MPI_Wait(&(send_rq[1]), MPI_STATUS_IGNORE);
            
            /* check non-blocking data size received from R */
            err = MPI_Wait(&(recv_rq[1]), MPI_STATUS_IGNORE);

            /* Post non-blocking receives for list structures from R Grid */
            err = MPI_Irecv(&(recv_buf_list1x1[0]),cnt, pD->LISTTYPE, pGrid->rx1_id,RtoL_tag, pD->Comm_Domain, &(recv_rq[1]));
            
            /* Send list structures R */
            err = MPI_Isend(&(send_buf_list1x1[0]),cnt, pD->LISTTYPE, pGrid->rx1_id,LtoR_tag, pD->Comm_Domain, &(send_rq[1]));
            
            /* check non-blocking sends have completed. */
            err = MPI_Wait(&(send_rq[1]), MPI_STATUS_IGNORE);
            
            /* check non-blocking data size received from R */
            err = MPI_Wait(&(recv_rq[1]), MPI_STATUS_IGNORE);
            
            /* Unpack */
            unpack_ox1_mc(pGrid);
            
            if (send_buf_list1x1) free_1d_array(send_buf_list1x1);
            if (recv_buf_list1x1) free_1d_array(recv_buf_list1x1);
            if (recv_buf_mc1x1) free_1d_array(recv_buf_mc1x1);
            if (send_buf_mc1x1) free_1d_array(send_buf_mc1x1);
            
        }
        
        /* MPI block on left, Physical boundary on right */
        if (pGrid->rx1_id < 0 && pGrid->lx1_id >= 0) {
            
            /* Post non-blocking receives for data size from L Grid */
            err = MPI_Irecv(&(recv_buf[0]), 1, MPI_INTEGER, pGrid->lx1_id, LtoR_tag, pD->Comm_Domain, &(recv_rq[0]));
            
            /* pack and send data size to L */
            send_buf[0] = pack_ix1_mc(pGrid);
            
            /* calloc tracer buffer for send to L */
            if((send_buf_mc0x1 = (Tracer_MPI *)calloc_1d_array(send_buf[0], sizeof(Tracer_MPI))) == NULL)
                ath_error("[bvals_init]: Failed to allocate recv buffer\n");
            
            err = MPI_Isend(&(send_buf[0]), 1, MPI_INTEGER, pGrid->lx1_id, RtoL_tag, pD->Comm_Domain, &(send_rq[0]));
            
            /* set physical boundary */
            (*(apply_mc_ox1))(pGrid);
            
            /* check non-blocking sends have completed. */
            err = MPI_Wait(&(send_rq[0]), MPI_STATUS_IGNORE);
            
            /* check non-blocking data size received from L */
            err = MPI_Wait(&(recv_rq[0]), MPI_STATUS_IGNORE);
            
            /* calloc buffer for tracer recevies from L */
            if((recv_buf_mc0x1 = (Tracer_MPI *)calloc_1d_array(recv_buf[0], sizeof(Tracer_MPI))) == NULL)
                ath_error("[bvals_init]: Failed to allocate recv buffer\n");
            
            /* Post non-blocking receives for tracers from L Grids */
            err = MPI_Irecv(&(recv_buf_mc0x1[0]),recv_buf[0],pD->TRACERTYPE,pGrid->lx1_id,LtoR_tag, pD->Comm_Domain, &(recv_rq[0]));

            /* pack and send tracers L */
            pack_ix1_tracers(pGrid);
            
            err = MPI_Isend(&(send_buf_mc0x1[0]),send_buf[0],pD->TRACERTYPE,pGrid->lx1_id,RtoL_tag, pD->Comm_Domain, &(send_rq[0]));

            /* check non-blocking sends have completed. */
            err = MPI_Wait(&(send_rq[0]), MPI_STATUS_IGNORE);
            
            /* check non-blocking data size received from L */
            err = MPI_Wait(&(recv_rq[0]), MPI_STATUS_IGNORE);

            /* Post non-blocking receives for list structures from L Grid */
            err = MPI_Irecv(&(recv_buf_list0x1[0]),cnt, pD->LISTTYPE, pGrid->lx1_id, LtoR_tag, pD->Comm_Domain, &(recv_rq[0]));
            
            /* Send list structures L */
            err = MPI_Isend(&(send_buf_list0x1[0]),cnt, pD->LISTTYPE, pGrid->lx1_id, RtoL_tag, pD->Comm_Domain, &(send_rq[0]));
            
            /* check non-blocking sends have completed. */
            err = MPI_Wait(&(send_rq[0]), MPI_STATUS_IGNORE);
            
            /* check non-blocking data size received from R */
            err = MPI_Wait(&(recv_rq[0]), MPI_STATUS_IGNORE);
            unpack_ix1_mc(pGrid);
            
            if (send_buf_list0x1) free_1d_array(send_buf_list0x1);
            if (recv_buf_list0x1) free_1d_array(recv_buf_list0x1);
            if (recv_buf_mc0x1) free_1d_array(recv_buf_mc0x1);
            if (send_buf_mc0x1) free_1d_array(send_buf_mc0x1);
        }
        
#endif /* MPI_PARALLEL */
        
        /* Physical boundaries on both left and right */
        if (pGrid->rx1_id < 0 && pGrid->lx1_id < 0) {
            (*(apply_mc_ix1))(pGrid);
            (*(apply_mc_ox1))(pGrid);
        }
    }
    /*--- Step 2. ------------------------------------------------------------------
     * Boundary Conditions in x2-direction */
    
    if (pGrid->Nx[1] > 1){
        
#ifdef MPI_PARALLEL
        cnt = (pGrid->Nx[0] + 2*nghost)*nghost*(pGrid->Nx[2]);

        if((send_buf_list0x2 = (TracerList_MPI *)calloc_1d_array(cnt,sizeof(TracerList_MPI))) == NULL)
            ath_error("[bvals_init]: Failed to allocate send buffer\n");
        if((send_buf_list1x2 = (TracerList_MPI*)calloc_1d_array(cnt,sizeof(TracerList_MPI))) == NULL)
            ath_error("[bvals_init]: Failed to allocate send buffer\n");
        if((recv_buf_list0x2 = (TracerList_MPI*)calloc_1d_array(cnt,sizeof(TracerList_MPI))) == NULL)
            ath_error("[bvals_init]: Failed to allocate send buffer\n");
        if((recv_buf_list1x2 = (TracerList_MPI*)calloc_1d_array(cnt,sizeof(TracerList_MPI))) == NULL)
            ath_error("[bvals_init]: Failed to allocate send buffer\n");

        if (pGrid->rx2_id >= 0 && pGrid->lx2_id >= 0) {
            
            /* Post non-blocking receives for data size from L and R Grids */
            err = MPI_Irecv(&(recv_buf[0]), 1, MPI_INTEGER, pGrid->lx2_id, LtoR_tag, pD->Comm_Domain, &(recv_rq[0]));
            err = MPI_Irecv(&(recv_buf[1]), 1, MPI_INTEGER, pGrid->rx2_id, RtoL_tag, pD->Comm_Domain, &(recv_rq[1]));
            
            /* pack and send data size to L and R */
            send_buf[0] = pack_ix2_mc(pGrid);
            
            /* calloc tracer buffer for send to l */
            if((send_buf_mc0x2 = (Tracer_MPI *)calloc_1d_array(send_buf[0], sizeof(Tracer_MPI))) == NULL)
                ath_error("[bvals_init]: Failed to allocate recv buffer\n");
            err = MPI_Isend(&(send_buf[0]), 1, MPI_INTEGER, pGrid->lx2_id, RtoL_tag, pD->Comm_Domain, &(send_rq[0]));

            /* pack data size */
            send_buf[1] = pack_ox2_mc(pGrid);
            
            /* calloc tracer buffer for send to r */
            if((send_buf_mc1x2 = (Tracer_MPI *)calloc_1d_array(send_buf[1], sizeof(Tracer_MPI))) == NULL)
                ath_error("[bvals_init]: Failed to allocate recv buffer\n");
            err = MPI_Isend(&(send_buf[1]), 1, MPI_INTEGER, pGrid->rx2_id, LtoR_tag, pD->Comm_Domain, &(send_rq[1]));
            
            /* check non-blocking sends have completed. */
            err = MPI_Waitall(2, send_rq, MPI_STATUS_IGNORE);
            
            /* check non-blocking data size receives */
            err = MPI_Waitany(2, recv_rq, &mIndex, MPI_STATUS_IGNORE);
            if (mIndex == 0) {
                /* calloc buffer for tracer receives from l */
                if((recv_buf_mc0x2 = (Tracer_MPI *)calloc_1d_array(recv_buf[0], sizeof(Tracer_MPI))) == NULL)
                    ath_error("[bvals_init]: Failed to allocate recv buffer\n");
            }
            if (mIndex == 1) {
                /* calloc buffer for tracer recevies from r */
                if((recv_buf_mc1x2 = (Tracer_MPI *)calloc_1d_array(recv_buf[1], sizeof(Tracer_MPI))) == NULL)
                    ath_error("[bvals_init]: Failed to allocate recv buffer\n");
            }
            err = MPI_Waitany(2, recv_rq, &mIndex,MPI_STATUS_IGNORE);
            if (mIndex == 0) {
                /* calloc buffer for tracer receives from l */
                if((recv_buf_mc0x2 = (Tracer_MPI *)calloc_1d_array(recv_buf[0], sizeof(Tracer_MPI))) == NULL)
                    ath_error("[bvals_init]: Failed to allocate recv buffer\n");
            }
            if (mIndex == 1) {
                /* calloc buffer for tracer recevies from r */
                if((recv_buf_mc1x2 = (Tracer_MPI *)calloc_1d_array(recv_buf[1], sizeof(Tracer_MPI))) == NULL)
                    ath_error("[bvals_init]: Failed to allocate recv buffer\n");
            }
            
            /* Post non-blocking receives for tracers from L and R Grids */
            err = MPI_Irecv(&(recv_buf_mc0x2[0]),recv_buf[0],pD->TRACERTYPE,pGrid->lx2_id,LtoR_tag,
                            pD->Comm_Domain, &(recv_rq[0]));
            err = MPI_Irecv(&(recv_buf_mc1x2[0]),recv_buf[1],pD->TRACERTYPE,pGrid->rx2_id,RtoL_tag,
                            pD->Comm_Domain, &(recv_rq[1]));
            
            /* pack and send data L and R */
            pack_ix2_tracers(pGrid);
            err = MPI_Isend(&(send_buf_mc0x2[0]),send_buf[0],pD->TRACERTYPE,pGrid->lx2_id,RtoL_tag, pD->Comm_Domain, &(send_rq[0]));
            
            pack_ox2_tracers(pGrid);
            err = MPI_Isend(&(send_buf_mc1x2[0]),send_buf[1],pD->TRACERTYPE,pGrid->rx2_id,LtoR_tag, pD->Comm_Domain, &(send_rq[1]));
            
            /* check non-blocking sends have completed. */
            err = MPI_Waitall(2,send_rq, MPI_STATUS_IGNORE);
            
            /* check non-blocking receives and unpack data in any order. */
            err = MPI_Waitany(2,recv_rq,&mIndex,MPI_STATUS_IGNORE);
            err = MPI_Waitany(2,recv_rq,&mIndex,MPI_STATUS_IGNORE);
        
            /* Post non-blocking receives for list structures from L and R Grids */
            
            err = MPI_Irecv(&(recv_buf_list0x2[0]),cnt, pD->LISTTYPE, pGrid->lx2_id,LtoR_tag, pD->Comm_Domain, &(recv_rq[0]));
            err = MPI_Irecv(&(recv_buf_list1x2[0]),cnt, pD->LISTTYPE, pGrid->rx2_id,RtoL_tag, pD->Comm_Domain, &(recv_rq[1]));
            
            /* Send list structures L and R */
            err = MPI_Isend(&(send_buf_list0x2[0]),cnt, pD->LISTTYPE,pGrid->lx2_id,RtoL_tag, pD->Comm_Domain, &(send_rq[0]));
            err = MPI_Isend(&(send_buf_list1x2[0]),cnt, pD->LISTTYPE,pGrid->rx2_id,LtoR_tag, pD->Comm_Domain, &(send_rq[1]));
            
            /* check non-blocking sends have completed. */
            err = MPI_Waitall(2, send_rq, MPI_STATUS_IGNORE);
            
            /* check non-blocking receives and unpack data in any order. */
            err = MPI_Waitany(2,recv_rq,&mIndex,MPI_STATUS_IGNORE);
            if (mIndex == 0) unpack_ix2_mc(pGrid);
            if (mIndex == 1) unpack_ox2_mc(pGrid);
            err = MPI_Waitany(2,recv_rq,&mIndex,MPI_STATUS_IGNORE);
            if (mIndex == 0) unpack_ix2_mc(pGrid);
            if (mIndex == 1) unpack_ox2_mc(pGrid);
            
            if (send_buf_list0x2) free_1d_array(send_buf_list0x2);
            if (send_buf_list1x2) free_1d_array(send_buf_list1x2);
            if (recv_buf_list0x2) free_1d_array(recv_buf_list0x2);
            if (recv_buf_list1x2) free_1d_array(recv_buf_list1x2);
            if (recv_buf_mc0x2) free_1d_array(recv_buf_mc0x2);
            if (recv_buf_mc1x2) free_1d_array(recv_buf_mc1x2);
            if (send_buf_mc0x2) free_1d_array(send_buf_mc0x2);
            if (send_buf_mc1x2) free_1d_array(send_buf_mc1x2);

        }

        /* Physical boundary on left, MPI block on right */
        if (pGrid->rx2_id >= 0 && pGrid->lx2_id < 0) {
            
            /* Post non-blocking receives for data size from R Grids */
            err = MPI_Irecv(&(recv_buf[1]), 1, MPI_INTEGER, pGrid->rx2_id, RtoL_tag, pD->Comm_Domain, &(recv_rq[1]));
            
            /* pack and send data size to R */
            send_buf[1] = pack_ox2_mc(pGrid);
            /* calloc tracer buffer for send to r */
            if((send_buf_mc1x2 = (Tracer_MPI *)calloc_1d_array(send_buf[1], sizeof(Tracer_MPI))) == NULL)
                ath_error("[bvals_init]: Failed to allocate recv buffer\n");
            err = MPI_Isend(&(send_buf[1]), 1, MPI_INTEGER, pGrid->rx2_id, LtoR_tag, pD->Comm_Domain, &(send_rq[1]));
            
            /* set physical boundary */
            (*(apply_mc_ix2))(pGrid);
            
            /* check non-blocking sends have completed. */
            err = MPI_Wait(&(send_rq[1]), MPI_STATUS_IGNORE);
            
            /* check non-blocking data size received from R */
            err = MPI_Wait(&(recv_rq[1]), MPI_STATUS_IGNORE);
            
            /* calloc buffer for tracer recevies from r */
            if((recv_buf_mc1x2 = (Tracer_MPI *)calloc_1d_array(recv_buf[1], sizeof(Tracer_MPI))) == NULL)
                ath_error("[bvals_init]: Failed to allocate recv buffer\n");
            
            /* Post non-blocking receives for tracers from R Grids */
            err = MPI_Irecv(&(recv_buf_mc1x2[0]),recv_buf[1],pD->TRACERTYPE,pGrid->rx2_id,RtoL_tag,
                            pD->Comm_Domain, &(recv_rq[1]));
            
            /* pack and send tracers R */
            pack_ox2_tracers(pGrid);
            err = MPI_Isend(&(send_buf_mc1x2[0]),send_buf[1],pD->TRACERTYPE,pGrid->rx2_id,LtoR_tag, pD->Comm_Domain, &(send_rq[1]));
            
            /* check non-blocking sends have completed. */
            err = MPI_Wait(&(send_rq[1]), MPI_STATUS_IGNORE);
            
            /* check non-blocking data size received from R */
            err = MPI_Wait(&(recv_rq[1]), MPI_STATUS_IGNORE);
            
            /* Post non-blocking receives for list structures from R Grid */
            err = MPI_Irecv(&(recv_buf_list1x2[0]),cnt, pD->LISTTYPE, pGrid->rx2_id,RtoL_tag, pD->Comm_Domain, &(recv_rq[1]));
            
            /* Send list structures R */
            err = MPI_Isend(&(send_buf_list1x2[0]),cnt, pD->LISTTYPE,pGrid->rx2_id,LtoR_tag, pD->Comm_Domain, &(send_rq[1]));
            
            /* check non-blocking sends have completed. */
            err = MPI_Wait(&(send_rq[1]), MPI_STATUS_IGNORE);
            
            /* check non-blocking data size received from R */
            err = MPI_Wait(&(recv_rq[1]), MPI_STATUS_IGNORE);
            unpack_ox2_mc(pGrid);
            
            if (send_buf_list1x2) free_1d_array(send_buf_list1x2);
            if (recv_buf_list1x2) free_1d_array(recv_buf_list1x2);
            if (recv_buf_mc1x2) free_1d_array(recv_buf_mc1x2);
            if (send_buf_mc1x2) free_1d_array(send_buf_mc1x2);

        }
        
        /* MPI block on left, Physical boundary on right */
        if (pGrid->rx2_id < 0 && pGrid->lx2_id >= 0) {
            
            /* Post non-blocking receives for data size from L Grid */
            err = MPI_Irecv(&(recv_buf[0]), 1, MPI_INTEGER, pGrid->lx2_id, LtoR_tag, pD->Comm_Domain, &(recv_rq[0]));
            
            /* pack and send data size to L */
            send_buf[0] = pack_ix2_mc(pGrid);
            /* calloc tracer buffer for send to L */
            if((send_buf_mc0x2 = (Tracer_MPI *)calloc_1d_array(send_buf[0], sizeof(Tracer_MPI))) == NULL)
                ath_error("[bvals_init]: Failed to allocate recv buffer\n");
            err = MPI_Isend(&(send_buf[0]), 1, MPI_INTEGER, pGrid->lx2_id, RtoL_tag, pD->Comm_Domain, &(send_rq[0]));
            
            /* set physical boundary */
            (*(apply_mc_ox2))(pGrid);
            
            /* check non-blocking sends have completed. */
            err = MPI_Wait(&(send_rq[0]), MPI_STATUS_IGNORE);
            
            /* check non-blocking data size received from L */
            err = MPI_Wait(&(recv_rq[0]), MPI_STATUS_IGNORE);
            
            /* calloc buffer for tracer recevies from L */
            if((recv_buf_mc0x2 = (Tracer_MPI *)calloc_1d_array(recv_buf[0], sizeof(Tracer_MPI))) == NULL)
                ath_error("[bvals_init]: Failed to allocate recv buffer\n");
            
            /* Post non-blocking receives for tracers from L Grids */
            err = MPI_Irecv(&(recv_buf_mc0x2[0]),recv_buf[0],pD->TRACERTYPE,pGrid->lx2_id,LtoR_tag,pD->Comm_Domain, &(recv_rq[0]));
            
            /* pack and send tracers L */
            pack_ix2_tracers(pGrid);
            err = MPI_Isend(&(send_buf_mc0x2[0]),send_buf[0],pD->TRACERTYPE,pGrid->lx2_id,RtoL_tag, pD->Comm_Domain, &(send_rq[0]));
            
            /* check non-blocking sends have completed. */
            err = MPI_Wait(&(send_rq[0]), MPI_STATUS_IGNORE);
            
            /* check non-blocking data size received from L */
            err = MPI_Wait(&(recv_rq[0]), MPI_STATUS_IGNORE);
            
            /* Post non-blocking receives for list structures from L Grid */
            err = MPI_Irecv(&(recv_buf_list0x2[0]),cnt, pD->LISTTYPE, pGrid->lx2_id,LtoR_tag, pD->Comm_Domain, &(recv_rq[0]));
            
            /* Send list structures L */
            err = MPI_Isend(&(send_buf_list0x2[0]),cnt, pD->LISTTYPE, pGrid->lx2_id,RtoL_tag, pD->Comm_Domain, &(send_rq[0]));
            
            /* check non-blocking sends have completed. */
            err = MPI_Wait(&(send_rq[0]), MPI_STATUS_IGNORE);
            
            /* check non-blocking data size received from R */
            err = MPI_Wait(&(recv_rq[0]), MPI_STATUS_IGNORE);
            unpack_ix2_mc(pGrid);
            
            if (send_buf_list0x2) free_1d_array(send_buf_list0x2);
            if (recv_buf_list0x2) free_1d_array(recv_buf_list0x2);
            if (recv_buf_mc0x2) free_1d_array(recv_buf_mc0x2);
            if (send_buf_mc0x2) free_1d_array(send_buf_mc0x2);
            
        }
        
#endif /* MPI_PARALLEL */
        
        /* Physical boundaries on both left and right */
        if (pGrid->rx2_id < 0 && pGrid->lx2_id < 0) {
            (*(apply_mc_ix2))(pGrid);
            (*(apply_mc_ox2))(pGrid);
        }
        
    }
    /*--- Step 3. ------------------------------------------------------------------
     * Boundary Conditions in x3-direction */
    
    if (pGrid->Nx[2] > 1){
 #ifdef MPI_PARALLEL
        
        cnt = (pGrid->Nx[0] + 2*nghost)*(pGrid->Nx[1] + 2*nghost)*nghost;
        
        if((send_buf_list0x3 = (TracerList_MPI*)calloc_1d_array(cnt,sizeof(TracerList_MPI))) == NULL)
            ath_error("[bvals_init]: Failed to allocate send buffer\n");
        if((send_buf_list1x3 = (TracerList_MPI*)calloc_1d_array(cnt,sizeof(TracerList_MPI))) == NULL)
            ath_error("[bvals_init]: Failed to allocate send buffer\n");
        if((recv_buf_list0x3 = (TracerList_MPI*)calloc_1d_array(cnt,sizeof(TracerList_MPI))) == NULL)
            ath_error("[bvals_init]: Failed to allocate send buffer\n");
        if((recv_buf_list1x3 = (TracerList_MPI*)calloc_1d_array(cnt,sizeof(TracerList_MPI))) == NULL)
            ath_error("[bvals_init]: Failed to allocate send buffer\n");
        
        /* MPI blocks to both left and right */
        if (pGrid->rx3_id >= 0 && pGrid->lx3_id >= 0) {
            
            /* Post non-blocking receives for data size from L and R Grids */
            err = MPI_Irecv(&(recv_buf[0]), 1, MPI_INTEGER, pGrid->lx3_id, LtoR_tag, pD->Comm_Domain, &(recv_rq[0]));
            err = MPI_Irecv(&(recv_buf[1]), 1, MPI_INTEGER, pGrid->rx3_id, RtoL_tag, pD->Comm_Domain, &(recv_rq[1]));
            
            /* pack and send data size to L and R */
            
            send_buf[0] = pack_ix3_mc(pGrid);
            
            /* calloc tracer buffer for send to l */
            if((send_buf_mc0x3 = (Tracer_MPI *)calloc_1d_array(send_buf[0], sizeof(Tracer_MPI))) == NULL)
                ath_error("[bvals_init]: Failed to allocate recv buffer\n");
            err = MPI_Isend(&(send_buf[0]), 1, MPI_INTEGER, pGrid->lx3_id, RtoL_tag, pD->Comm_Domain, &(send_rq[0]));
            
            send_buf[1] = pack_ox3_mc(pGrid);
            
            /* calloc tracer buffer for send to r */
            if((send_buf_mc1x3 = (Tracer_MPI *)calloc_1d_array(send_buf[1], sizeof(Tracer_MPI))) == NULL)
                ath_error("[bvals_init]: Failed to allocate recv buffer\n");
            
            err = MPI_Isend(&(send_buf[1]), 1, MPI_INTEGER, pGrid->rx3_id, LtoR_tag, pD->Comm_Domain, &(send_rq[1]));
            
            /* check non-blocking sends have completed. */
            err = MPI_Waitall(2, send_rq, MPI_STATUS_IGNORE);
            
            /* check non-blocking data size receives */
            err = MPI_Waitany(2, recv_rq, &mIndex,MPI_STATUS_IGNORE);
            
            if (mIndex == 0) {
                /* calloc buffer for tracer receives from l */
                if((recv_buf_mc0x3 = (Tracer_MPI *)calloc_1d_array(recv_buf[0], sizeof(Tracer_MPI))) == NULL)
                    ath_error("[bvals_init]: Failed to allocate recv buffer\n");
            }
            if (mIndex == 1) {
                /* calloc buffer for tracer recevies from r */
                if((recv_buf_mc1x3 = (Tracer_MPI *)calloc_1d_array(recv_buf[1], sizeof(Tracer_MPI))) == NULL)
                    ath_error("[bvals_init]: Failed to allocate recv buffer\n");
            }
            err = MPI_Waitany(2, recv_rq, &mIndex,MPI_STATUS_IGNORE);
            if (mIndex == 0) {
                /* calloc buffer for tracer receives from l */
                if((recv_buf_mc0x3 = (Tracer_MPI *)calloc_1d_array(recv_buf[0], sizeof(Tracer_MPI))) == NULL)
                    ath_error("[bvals_init]: Failed to allocate recv buffer\n");
            }
            if (mIndex == 1) {
                /* calloc buffer for tracer recevies from r */
                if((recv_buf_mc1x3 = (Tracer_MPI *)calloc_1d_array(recv_buf[1], sizeof(Tracer_MPI))) == NULL)
                    ath_error("[bvals_init]: Failed to allocate recv buffer\n");
            }
            
            /* Post non-blocking receives for tracers from L and R Grids */
            err = MPI_Irecv(&(recv_buf_mc0x3[0]),recv_buf[0],pD->TRACERTYPE,pGrid->lx3_id,LtoR_tag,
                            pD->Comm_Domain, &(recv_rq[0]));
            err = MPI_Irecv(&(recv_buf_mc1x3[0]),recv_buf[1],pD->TRACERTYPE,pGrid->rx3_id,RtoL_tag,
                            pD->Comm_Domain, &(recv_rq[1]));
            
            /* pack and send data L and R */
            pack_ix3_tracers(pGrid);
            err = MPI_Isend(&(send_buf_mc0x3[0]),send_buf[0],pD->TRACERTYPE,pGrid->lx3_id,RtoL_tag, pD->Comm_Domain, &(send_rq[0]));
            
            pack_ox3_tracers(pGrid);
            err = MPI_Isend(&(send_buf_mc1x3[0]),send_buf[1],pD->TRACERTYPE,pGrid->rx3_id,LtoR_tag, pD->Comm_Domain, &(send_rq[1]));
            
            /* check non-blocking sends have completed. */
            err = MPI_Waitall(2,send_rq, MPI_STATUS_IGNORE);
            
            /* check non-blocking receives and unpack data in any order. */
            err = MPI_Waitany(2,recv_rq,&mIndex,MPI_STATUS_IGNORE);
            err = MPI_Waitany(2,recv_rq,&mIndex,MPI_STATUS_IGNORE);
            
            /* Post non-blocking receives for list structures from L and R Grids */
            
            err = MPI_Irecv(&(recv_buf_list0x3[0]),cnt, pD->LISTTYPE, pGrid->lx3_id,LtoR_tag, pD->Comm_Domain, &(recv_rq[0]));
            err = MPI_Irecv(&(recv_buf_list1x3[0]),cnt, pD->LISTTYPE, pGrid->rx3_id,RtoL_tag, pD->Comm_Domain, &(recv_rq[1]));
            
            /* Send list structures L and R */
            err = MPI_Isend(&(send_buf_list0x3[0]),cnt, pD->LISTTYPE,pGrid->lx3_id,RtoL_tag, pD->Comm_Domain, &(send_rq[0]));
            err = MPI_Isend(&(send_buf_list1x3[0]),cnt, pD->LISTTYPE,pGrid->rx3_id,LtoR_tag, pD->Comm_Domain, &(send_rq[1]));
            
            /* check non-blocking sends have completed. */
            err = MPI_Waitall(2, send_rq, MPI_STATUS_IGNORE);
            
            /* check non-blocking receives and unpack data in any order. */
            err = MPI_Waitany(2,recv_rq,&mIndex,MPI_STATUS_IGNORE);
            if (mIndex == 0) unpack_ix3_mc(pGrid);
            if (mIndex == 1) unpack_ox3_mc(pGrid);
            err = MPI_Waitany(2,recv_rq,&mIndex,MPI_STATUS_IGNORE);
            if (mIndex == 0) unpack_ix3_mc(pGrid);
            if (mIndex == 1) unpack_ox3_mc(pGrid);
            
            if (send_buf_list0x3) free_1d_array(send_buf_list0x3);
            if (send_buf_list1x3) free_1d_array(send_buf_list1x3);
            if (recv_buf_list0x3) free_1d_array(recv_buf_list0x3);
            if (recv_buf_list1x3) free_1d_array(recv_buf_list1x3);
            if (recv_buf_mc0x3) free_1d_array(recv_buf_mc0x3);
            if (recv_buf_mc1x3) free_1d_array(recv_buf_mc1x3);
            if (send_buf_mc0x3) free_1d_array(send_buf_mc0x3);
            if (send_buf_mc1x3) free_1d_array(send_buf_mc1x3);
        }
        
        /* Physical boundary on left, MPI block on right */
        if (pGrid->rx3_id >= 0 && pGrid->lx3_id < 0) {

            /* Post non-blocking receives for data size from R Grids */
            err = MPI_Irecv(&(recv_buf[1]), 1, MPI_INTEGER, pGrid->rx3_id, RtoL_tag, pD->Comm_Domain, &(recv_rq[1]));
            
            /* pack and send data size to R */
            send_buf[1] = pack_ox3_mc(pGrid);
            /* calloc tracer buffer for send to r */
            if((send_buf_mc1x3 = (Tracer_MPI *)calloc_1d_array(send_buf[1], sizeof(Tracer_MPI))) == NULL)
                ath_error("[bvals_init]: Failed to allocate recv buffer\n");

            err = MPI_Isend(&(send_buf[1]), 1, MPI_INTEGER, pGrid->rx3_id, LtoR_tag, pD->Comm_Domain, &(send_rq[1]));
            
            /* set physical boundary */
            (*(apply_mc_ix3))(pGrid);
            
            /* check non-blocking sends have completed. */
            err = MPI_Wait(&(send_rq[1]), MPI_STATUS_IGNORE);
            
            /* check non-blocking data size received from R */
            err = MPI_Wait(&(recv_rq[1]), MPI_STATUS_IGNORE);
            
            /* calloc buffer for tracer recevies from r */
            if((recv_buf_mc1x3 = (Tracer_MPI *)calloc_1d_array(recv_buf[1], sizeof(Tracer_MPI))) == NULL)
                ath_error("[bvals_init]: Failed to allocate recv buffer\n");

            /* Post non-blocking receives for tracers from R Grids */
            err = MPI_Irecv(&(recv_buf_mc1x3[0]),recv_buf[1],pD->TRACERTYPE,pGrid->rx3_id,RtoL_tag, pD->Comm_Domain, &(recv_rq[1]));
            
            /* pack and send tracers R */
            pack_ox3_tracers(pGrid);
            
            err = MPI_Isend(&(send_buf_mc1x3[0]),send_buf[1],pD->TRACERTYPE,pGrid->rx3_id,LtoR_tag, pD->Comm_Domain, &(send_rq[1]));

            /* check non-blocking sends have completed. */
            err = MPI_Wait(&(send_rq[1]), MPI_STATUS_IGNORE);

            /* check non-blocking data size received from R */
            err = MPI_Wait(&(recv_rq[1]), MPI_STATUS_IGNORE);
            
            /* Post non-blocking receives for list structures from R Grid */
            err = MPI_Irecv(&(recv_buf_list1x3[0]),cnt, pD->LISTTYPE, pGrid->rx3_id,RtoL_tag, pD->Comm_Domain, &(recv_rq[1]));

            /* Send list structures R */
            err = MPI_Isend(&(send_buf_list1x3[0]),cnt, pD->LISTTYPE,pGrid->rx3_id,LtoR_tag, pD->Comm_Domain, &(send_rq[1]));
            
            /* check non-blocking sends have completed. */
            err = MPI_Wait(&(send_rq[1]), MPI_STATUS_IGNORE);
            
            /* check non-blocking data size received from R */
            err = MPI_Wait(&(recv_rq[1]), MPI_STATUS_IGNORE);
            unpack_ox3_mc(pGrid);
            
            if (send_buf_list1x3) free_1d_array(send_buf_list1x3);
            if (recv_buf_list1x3) free_1d_array(recv_buf_list1x3);
            if (recv_buf_mc1x3) free_1d_array(recv_buf_mc1x3);
            if (send_buf_mc1x3) free_1d_array(send_buf_mc1x3);
        }

        /* MPI block on left, Physical boundary on right */
        if (pGrid->rx3_id < 0 && pGrid->lx3_id >= 0) {

            /* Post non-blocking receives for data size from L Grid */
            err = MPI_Irecv(&(recv_buf[0]), 1, MPI_INTEGER, pGrid->lx3_id, LtoR_tag, pD->Comm_Domain, &(recv_rq[0]));
            
            /* pack and send data size to L */
            send_buf[0] = pack_ix3_mc(pGrid);
            /* calloc tracer buffer for send to L */
            if((send_buf_mc0x3 = (Tracer_MPI *)calloc_1d_array(send_buf[0], sizeof(Tracer_MPI))) == NULL)
                ath_error("[bvals_init]: Failed to allocate recv buffer\n");

            err = MPI_Isend(&(send_buf[0]), 1, MPI_INTEGER, pGrid->lx3_id, RtoL_tag, pD->Comm_Domain, &(send_rq[0]));
            
            /* set physical boundary */
            (*(apply_mc_ox3))(pGrid);

            /* check non-blocking sends have completed. */
            err = MPI_Wait(&(send_rq[0]), MPI_STATUS_IGNORE);
            
            /* check non-blocking data size received from L */
            err = MPI_Wait(&(recv_rq[0]), MPI_STATUS_IGNORE);

            /* calloc buffer for tracer recevies from L */
            if((recv_buf_mc0x3 = (Tracer_MPI *)calloc_1d_array(recv_buf[0], sizeof(Tracer_MPI))) == NULL)
                ath_error("[bvals_init]: Failed to allocate recv buffer\n");
            
            /* Post non-blocking receives for tracers from L Grids */
            err = MPI_Irecv(&(recv_buf_mc0x3[0]),recv_buf[0],pD->TRACERTYPE,pGrid->lx3_id,LtoR_tag, pD->Comm_Domain, &(recv_rq[0]));

            /* pack and send tracers L */
            pack_ix3_tracers(pGrid);

            err = MPI_Isend(&(send_buf_mc0x3[0]),send_buf[0],pD->TRACERTYPE,pGrid->lx3_id,RtoL_tag, pD->Comm_Domain, &(send_rq[0]));

            /* check non-blocking sends have completed. */
            err = MPI_Wait(&(send_rq[0]), MPI_STATUS_IGNORE);
            
            /* check non-blocking data size received from L */
            err = MPI_Wait(&(recv_rq[0]), MPI_STATUS_IGNORE);
            
            /* Post non-blocking receives for list structures from L Grid */
            err = MPI_Irecv(&(recv_buf_list0x3[0]),cnt, pD->LISTTYPE, pGrid->lx3_id, LtoR_tag, pD->Comm_Domain, &(recv_rq[0]));

            /* Send list structures L */
            err = MPI_Isend(&(send_buf_list0x3[0]),cnt, pD->LISTTYPE, pGrid->lx3_id, RtoL_tag, pD->Comm_Domain, &(send_rq[0]));
            
            /* check non-blocking sends have completed. */
            err = MPI_Wait(&(send_rq[0]), MPI_STATUS_IGNORE);
            
            /* check non-blocking data size received from R */
            err = MPI_Wait(&(recv_rq[0]), MPI_STATUS_IGNORE);
            unpack_ix3_mc(pGrid);
            
            if (send_buf_list0x3) free_1d_array(send_buf_list0x3);
            if (recv_buf_list0x3) free_1d_array(recv_buf_list0x3);
            if (recv_buf_mc0x3) free_1d_array(recv_buf_mc0x3);
            if (send_buf_mc0x3) free_1d_array(send_buf_mc0x3);
            
        }
      
#endif /* MPI_PARALLEL */
        
        /* Physical boundaries on both left and right */
        if (pGrid->rx3_id < 0 && pGrid->lx3_id < 0) {
            
            (*(apply_mc_ix3))(pGrid);
            (*(apply_mc_ox3))(pGrid);
        }
    }
    
#ifdef MPI_PARALLEL
    
    if (recv_buf) free_1d_array(recv_buf);
    if (send_buf) free_1d_array(send_buf);
    if (recv_rq) free_1d_array(recv_rq);
    if (send_rq) free_1d_array(send_rq);
    
#endif /* MPI_PARALLEL */

    return;
}

/*----------------------------------------------------------------------------*/
/*! \fn void bvals_mc_init(MeshS *pM)
 *  \brief Sets function pointers for physical boundaries during
 *   initialization, allocates memory for send/receive buffers with MPI
 */

void bvals_tracer_init(MeshS *pM)
{
    GridS *pG;
    DomainS *pD;
    
    pD = &(pM->Domain[0][0]);
    pG = pD->Grid;
    
    /* Set function pointers for physical boundaries in x1-direction */
    if(pG->Nx[0] > 1) {
        if(apply_mc_ix1 == NULL){
            switch(pM->BCFlag_ix1){
                    
                case 1: /* Reflecting */
                    apply_mc_ix1 = reflect_ix1_mc;
                    break;
                case 5: /* Reflecting */
                    apply_mc_ix1 = reflect_ix1_mc;
                    break;
                case 2: /* Outflow */
                    apply_mc_ix1 = outflow_ix1_mc;
                    break;
                case 4: /* Periodic */
                    apply_mc_ix1 = periodic_ix1_mc;
                    break;
                default:
                    ath_perr(-1,"[set_bvals_mc_init]: bc_ix1 = %d unknown\n",
                             pM->BCFlag_ix1);
                    exit(EXIT_FAILURE);
            }
        }
        
        if(apply_mc_ox1 == NULL){
            switch(pM->BCFlag_ox1){
                case 1: /* Reflecting */
                    apply_mc_ox1 = reflect_ox1_mc;
                    break;
                case 5: /* Reflecting */
                    apply_mc_ox1 = reflect_ox1_mc;
                    break;
                case 2: /* Outflow */
                    apply_mc_ox1 = outflow_ox1_mc;
                    break;
                case 4: /* Periodic */
                    apply_mc_ox1 = periodic_ox1_mc;
                    break;
                    
                default:
                    ath_perr(-1,"[set_bvals_mc_init]: bc_ox1 = %d unknown\n",
                             pM->BCFlag_ox1);
                    exit(EXIT_FAILURE);
            }
            
        }
    }
        /* Set function pointers for physical boundaries in x2-direction */
        
        if(pG->Nx[1] > 1) {
            if(apply_mc_ix2 == NULL){
                switch(pM->BCFlag_ix2){
                    case 1: /* Reflecting */
                        apply_mc_ix2 = reflect_ix2_mc;
                        break;
                    case 5: /* Reflecting */
                        apply_mc_ix2 = reflect_ix2_mc;
                        break;
                    case 2: /* Outflow */
                        apply_mc_ix2 = outflow_ix2_mc;
                        break;
                    case 4: /* Periodic */
                        apply_mc_ix2 = periodic_ix2_mc;
                        break;
                    default:
                        ath_perr(-1,"[set_bvals_mc_init]: bc_ix2 = %d unknown\n",
                                 pM->BCFlag_ix2);
                        exit(EXIT_FAILURE);
                }
            }
            
            if(apply_mc_ox2 == NULL){
                switch(pM->BCFlag_ox2){
                    case 1: /* Reflecting */
                        apply_mc_ox2 = reflect_ox2_mc;
                        break;
                    case 5: /* Reflecting */
                        apply_mc_ox2 = reflect_ox2_mc;
                        break;
                    case 2: /* Outflow */
                        apply_mc_ox2 = outflow_ox2_mc;
                        break;
                    case 4: /* Periodic */
                        apply_mc_ox2 = periodic_ox2_mc;
                        break;
                    default:
                        ath_perr(-1,"[set_bvals_mc_init]: bc_ox2 = %d unknown\n",
                                 pM->BCFlag_ox2);
                        exit(EXIT_FAILURE);
                }
                
            }
        }
    /* Set function pointers for physical boundaries in x3-direction */
    if(pG->Nx[2] > 1) {
        if(apply_mc_ix3 == NULL){
            switch(pM->BCFlag_ix3){
                case 1: /* Reflecting */
                    apply_mc_ix3 = reflect_ix3_mc;
                    break;
                case 5: /* Reflecting */
                    apply_mc_ix3 = reflect_ix3_mc;
                    break;
                case 2: /* Outflow */
                    apply_mc_ix3 = outflow_ix3_mc;
                    break;
                case 4: /* Periodic */
                    apply_mc_ix3 = periodic_ix3_mc;
                    break;
                default:
                    ath_perr(-1,"[set_bvals_mc_init]: bc_ix3 = %d unknown\n",
                             pM->BCFlag_ix3);
                    exit(EXIT_FAILURE);
            }
            
        }
        if(apply_mc_ox3 == NULL){
            switch(pM->BCFlag_ox3){
                case 1: /* Reflecting */
                    apply_mc_ox3 = reflect_ox3_mc;
                    break;
                case 5: /* Reflecting */
                    apply_mc_ox3 = reflect_ox3_mc;
                    break;
                case 2: /* Outflow */
                    apply_mc_ox3 = outflow_ox3_mc;
                    break;
                case 4: /* Periodic */
                    apply_mc_ox3 = periodic_ox3_mc;
                    break;
                default:
                    ath_perr(-1,"[set_bvals_mc_init]: bc_ox3 = %d unknown\n",
                             pM->BCFlag_ox3);
                    exit(EXIT_FAILURE);
            }
            
        }
    }
    return;
}

static void reflect_ix1_mc(GridS *pG)
{
//    int is = pG->is;
//    int js = pG->js, je = pG->je;
//    int ks = pG->ks, ke = pG->ke;
//    int i,j,k;
//    TracerS *pnode;
//    TracerListS *list, *newList;
//    
//    for (k=ks; k<=ke; k++) {
//        for (j=js; j<=je; j++) {
//            for (i=1; i<=nghost; i++) {
//                
//                /* Tracer list to move (in ghost zone) */
//                list = &((pG->GridLists)[k][j][is-i]);
//                /* New list in active zone */
//                newList = &((pG->GridLists)[k][j][is+(i-1)]);
//                pnode = list->Head;
//                
//                /* flag tracers for movement */
//                while(pnode) {
//                    pnode->newList = newList;
//                    pnode = pnode->Next;
//                }
//                /* Move tracers */
//                Tracerlist_sweep_bc(list);
//            }
//        }
//    }
    return;
}

static void reflect_ox1_mc(GridS *pG)
{
//    int  ie = pG->ie;
//    int js = pG->js, je = pG->je;
//    int ks = pG->ks, ke = pG->ke;
//    int i,j,k;
//    TracerS *pnode;
//    TracerListS *list, *newList;
//
//    for (k=ks; k<=ke; k++) {
//        for (j=js; j<=je; j++) {
//            for (i=1; i<=nghost; i++) {
//                /* Tracer list to move (in ghost zone) */
//                list = &((pG->GridLists)[k][j][ie+i]);
//                /* New list in active zone */
//                newList = &((pG->GridLists)[k][j][ie-(i-1)]);
//                pnode = list->Head;
//                /* flag tracers for movement */
//                while(pnode) {
//                    pnode->newList = newList;
//                    pnode = pnode->Next;
//                }
//                /* Move tracers */
//                Tracerlist_sweep_bc(list);
//            }
//        }
//    }
    return;
}

static void reflect_ix2_mc(GridS *pG)
{
    int is = pG->is, ie = pG->ie;
    int js = pG->js;
    int ks = pG->ks, ke = pG->ke;
    int i,j,k;
    TracerS *pnode;
    TracerListS *list, *newList;
    
    for (k=ks; k<=ke; k++) {
        for (j=1; j<=nghost; j++) {
            for (i=is-nghost; i<=ie+nghost; i++) {
                /* Tracer list to move (in ghost zone) */
                list = &((pG->GridLists)[k][js-j][i]);
                /* New list in active zone */
                newList = &((pG->GridLists)[k][js+(j-1)][i]);
                pnode = list->Head;
                /* flag tracers for movement */
                while(pnode) {
                    pnode->newList = newList;
                    pnode = pnode->Next;
                }
                /* Move tracers */
                Tracerlist_sweep_bc(list);
            }
        }
    }
    return;
}

static void reflect_ox2_mc(GridS *pG)
{
    int is = pG->is, ie = pG->ie;
    int je = pG->je;
    int ks = pG->ks, ke = pG->ke;
    int i,j,k;
    TracerS *pnode;
    TracerListS *list, *newList;

    for (k=ks; k<=ke; k++) {
        for (j=1; j<=nghost; j++) {
            for (i=is-nghost; i<=ie+nghost; i++) {
                /* Tracer list to move (in ghost zone) */
                list = &((pG->GridLists)[k][je+j][i]);
                /* New list in active zone */
                newList = &((pG->GridLists)[k][je-(j-1)][i]);
                pnode = list->Head;
                /* flag tracers for movement */
                while(pnode) {
                    pnode->newList = newList;
                    pnode = pnode->Next;
                }
                /* Move tracers */
                Tracerlist_sweep_bc(list);
            }
        }
    }
    return;
}

static void reflect_ix3_mc(GridS *pG)
{
    int is = pG->is, ie = pG->ie;
    int js = pG->js, je = pG->je;
    int ks = pG->ks;
    int i,j,k;
    TracerS *pnode;
    TracerListS *list, *newList;

    for (k=1; k<=nghost; k++) {
        for (j=js-nghost; j<=je+nghost; j++) {
            for (i=is-nghost; i<=ie+nghost; i++) {
                /* Tracer list to move (in ghost zone) */
                list = &((pG->GridLists)[ks-k][j][i]);
                /* New list in active zone */
                newList = &((pG->GridLists)[ks+(k-1)][j][i]);
                pnode = list->Head;
                /* flag tracers for movement */
                while(pnode) {
                    pnode->newList = newList;
                    pnode = pnode->Next;
                }
                /* Move tracers */
                Tracerlist_sweep_bc(list);
            }
        }
    }
    return;
}

static void reflect_ox3_mc(GridS *pG)
{
    int is = pG->is, ie = pG->ie;
    int js = pG->js, je = pG->je;
    int ke = pG->ke;
    int i,j,k;
    TracerS *pnode;
    TracerListS *list, *newList;
    
    for (k=1; k<=nghost; k++) {
        for (j=js-nghost; j<=je+nghost; j++) {
            for (i=is-nghost; i<=ie+nghost; i++) {
                /* Tracer list to move (in ghost zone) */
                list = &((pG->GridLists)[ke+k][j][i]);
                /* New list in active zone */
                newList = &((pG->GridLists)[ke-(k-1)][j][i]);
                pnode = list->Head;
                /* flag tracers for movement */
                while(pnode) {
                    pnode->newList = newList;
                    pnode = pnode->Next;
                }
                /* Move tracers */
                Tracerlist_sweep_bc(list);
            }
        }
    }
}

static void outflow_ix1_mc(GridS *pG)
{
    return;
}
static void outflow_ox1_mc(GridS *pG)
{
    return;
}
static void outflow_ix2_mc(GridS *pG)
{
    return;
}
static void outflow_ox2_mc(GridS *pG)
{
    return;
}
static void outflow_ix3_mc(GridS *pG)
{
    return;
}
static void outflow_ox3_mc(GridS *pG)
{
    return;
}

#ifdef MPI_PARALLEL
static int pack_ix1_mc(GridS *pG)
{
    int is = pG->is;
    int js = pG->js, je = pG->je;
    int ks = pG->ks, ke = pG->ke;
    int i,j,k,n;
    TracerListS *list;
    
    n = 0;
    for (k=ks; k<=ke; k++){
        for (j=js; j<=je; j++){
            i = is-1;
            list = &((pG->GridLists)[k][j][i]);
            n += list->count;
        }
    }
    return n;
}

static void pack_ix1_tracers(GridS *pG)
{
    int is = pG->is;
    int js = pG->js, je = pG->je;
    int ks = pG->ks, ke = pG->ke;
    int i,j,k,l,m;
    TracerListS *list;
    TracerS *tracer;
    TracerS *Next;
    
    l = 0;
    m = 0;
    for (k=ks; k<=ke; k++){
        for (j=js; j<=je; j++){
            i = is-1;
            list = &((pG->GridLists)[k][j][i]);
            send_buf_list0x1[m].count = list->count;
            m++;
            tracer = list->Head;
            while(tracer) {
                send_buf_mc0x1[l].id = (double)tracer->prop->id;
                send_buf_mc0x1[l].d_init = (double)tracer->prop->d_init;
                send_buf_mc0x1[l].i_init = (int)tracer->prop->i_init;
                send_buf_mc0x1[l].j_init = (int)tracer->prop->j_init;
                send_buf_mc0x1[l].k_init = (int)tracer->prop->k_init;
#ifdef STAR_PARTICLE
                send_buf_mc0x1[l].star_id = (int)tracer->prop->star_id;
#endif /* STAR_PARTICLE */
#ifdef VFTRACERS
                send_buf_mc0x1[l].x1 = (double)tracer->x1;
                send_buf_mc0x1[l].x2 = (double)tracer->x2;
                send_buf_mc0x1[l].x3 = (double)tracer->x3;
#endif /* VFTRACERS */
                l++;
                Next = tracer->Next;
                tracer = Next;
            }
            tracer = list->Head;
            while(tracer) {
                Next = tracer->Next;
                tracer_list_remove(list, tracer);
                free(tracer);
                tracer = Next;
            }
#ifdef DEBUG
            assert(list->count == 0);
#endif /* DEBUG */
    }
}
#ifdef DEBUG
    assert(l == send_buf[0]);
#endif /* DEBUG */
}


static int pack_ox1_mc(GridS *pG)
{
    int is = pG->is, ie = pG->ie;
    int js = pG->js, je = pG->je;
    int ks = pG->ks, ke = pG->ke;
    int i,j,k,n;
    TracerListS *list;
    n = 0;
    for (k=ks; k<=ke; k++) {
        for (j=js; j<=je; j++) {
            i = ie+1;
            list = &((pG->GridLists)[k][j][i]);
            n += list->count;
        }
    }
    return n;
}

static void pack_ox1_tracers(GridS *pG)
{
    int is = pG->is, ie = pG->ie;
    int js = pG->js, je = pG->je;
    int ks = pG->ks, ke = pG->ke;
    int i,j,k,l,m;
    TracerListS *list;
    TracerS *tracer;
    TracerS *Next;

    l=0;
    m=0;
    for (k=ks; k<=ke; k++) {
        for (j=js; j<=je; j++) {
            i = ie+1;
            list = &((pG->GridLists)[k][j][i]);
            send_buf_list1x1[m].count = list->count;
            m++;
            tracer = list->Head;
            while(tracer) {
                send_buf_mc1x1[l].id = (double)tracer->prop->id;
                send_buf_mc1x1[l].d_init = (double)tracer->prop->d_init;
                send_buf_mc1x1[l].i_init = (int)tracer->prop->i_init;
                send_buf_mc1x1[l].j_init = (int)tracer->prop->j_init;
                send_buf_mc1x1[l].k_init = (int)tracer->prop->k_init;
#ifdef STAR_PARTICLE
                send_buf_mc1x1[l].star_id = (int)tracer->prop->star_id;
#endif /* STAR_PARTICLE */
#ifdef VFTRACERS
                send_buf_mc1x1[l].x1 = (double)tracer->x1;
                send_buf_mc1x1[l].x2 = (double)tracer->x2;
                send_buf_mc1x1[l].x3 = (double)tracer->x3;
#endif /* VFTRACERS */
                l++;
                Next = tracer->Next;
                tracer = Next;
            }
            tracer = list->Head;
            while(tracer) {
                Next = tracer->Next;
                tracer_list_remove(list, tracer);
                free(tracer);
                tracer = Next;
            }
#ifdef DEBUG
            assert(list->count == 0);
#endif /* DEBUG */
        }
    }
#ifdef DEBUG
    assert(l == send_buf[1]);
#endif /* DEBUG */
}


static int pack_ix2_mc(GridS *pG)
{
    int is = pG->is, ie = pG->ie;
    int js = pG->js, je = pG->je;
    int ks = pG->ks, ke = pG->ke;
    int i,j,k,n;
    TracerListS *list;

    n = 0;
    for (k=ks; k<=ke; k++) {
        j = js-1;
        for (i=is-nghost; i<=ie+nghost; i++) {
            list = &((pG->GridLists)[k][j][i]);
            n += list->count;
        }
    }
    return n;
}

static void pack_ix2_tracers(GridS *pG)
{
    int is = pG->is, ie = pG->ie;
    int js = pG->js, je = pG->je;
    int ks = pG->ks, ke = pG->ke;
    int i,j,k,l,m;
    TracerListS *list;
    TracerS *tracer;
    TracerS *Next;
    
    l = 0;
    m = 0;
    for (k=ks; k<=ke; k++) {
        j = js-1;
        for (i=is-nghost; i<=ie+nghost; i++) {
            list = &((pG->GridLists)[k][j][i]);
            send_buf_list0x2[m].count = list->count;
            m++;
            tracer = list->Head;
            while(tracer) {
                send_buf_mc0x2[l].id = (double)tracer->prop->id;
                send_buf_mc0x2[l].d_init = (double)tracer->prop->d_init;
                send_buf_mc0x2[l].i_init = (int)tracer->prop->i_init;
                send_buf_mc0x2[l].j_init = (int)tracer->prop->j_init;
                send_buf_mc0x2[l].k_init = (int)tracer->prop->k_init;
#ifdef STAR_PARTICLE
                send_buf_mc0x2[l].star_id = (int)tracer->prop->star_id;
#endif /* STAR_PARTICLE */
#ifdef VFTRACERS
                send_buf_mc0x2[l].x1 = (double)tracer->x1;
                send_buf_mc0x2[l].x2 = (double)tracer->x2;
                send_buf_mc0x2[l].x3 = (double)tracer->x3;
#endif /* VFTRACERS */
                l++;
                Next = tracer->Next;
                tracer = Next;
            }
            tracer = list->Head;
            while(tracer) {
                Next = tracer->Next;
                tracer_list_remove(list, tracer);
                free(tracer);
                tracer = Next;
            }
        }
    }
}

static int pack_ox2_mc(GridS *pG)
{
    int is = pG->is, ie = pG->ie;
    int js = pG->js, je = pG->je;
    int ks = pG->ks, ke = pG->ke;
    int i,j,k,n;
    TracerListS *list;
    
    n = 0;
    for (k=ks; k<=ke; k++) {
        j = je+1;
        for (i=is-nghost; i<=ie+nghost; i++) {
            list = &((pG->GridLists)[k][j][i]);
            n += list->count;
        }
    }
    return n;
}

static void pack_ox2_tracers(GridS *pG)
{
    int is = pG->is, ie = pG->ie;
    int js = pG->js, je = pG->je;
    int ks = pG->ks, ke = pG->ke;
    int i,j,k,l,m;
    TracerListS *list;
    TracerS *tracer;
    TracerS *Next;
    
    l = 0;
    m = 0;
    for (k=ks; k<=ke; k++) {
        j = je+1;
        for (i=is-nghost; i<=ie+nghost; i++) {
            list = &((pG->GridLists)[k][j][i]);
            send_buf_list1x2[m].count = list->count;
            m++;
            tracer = list->Head;
            while(tracer) {
                send_buf_mc1x2[l].id = (double)tracer->prop->id;
                send_buf_mc1x2[l].d_init = (double)tracer->prop->d_init;
                send_buf_mc1x2[l].i_init = (int)tracer->prop->i_init;
                send_buf_mc1x2[l].j_init = (int)tracer->prop->j_init;
                send_buf_mc1x2[l].k_init = (int)tracer->prop->k_init;
#ifdef STAR_PARTICLE
                send_buf_mc1x2[l].star_id = (int)tracer->prop->star_id;
#endif /* STAR_PARTICLE */
#ifdef VFTRACERS
                send_buf_mc1x2[l].x1 = (double)tracer->x1;
                send_buf_mc1x2[l].x2 = (double)tracer->x2;
                send_buf_mc1x2[l].x3 = (double)tracer->x3;
#endif /* VFTRACERS */
                l++;
                Next = tracer->Next;
                tracer = Next;
            }
            tracer = list->Head;
            while(tracer) {
                Next = tracer->Next;
                tracer_list_remove(list, tracer);
                free(tracer);
                tracer = Next;
            }
        }
    }
}

static int pack_ix3_mc(GridS *pG)
{
    int is = pG->is, ie = pG->ie;
    int js = pG->js, je = pG->je;
    int ks = pG->ks, ke = pG->ke;
    int i,j,k,n;
    TracerListS *list;
    
    n = 0;
    k = ks-1;
    for (j=js-nghost; j<=je+nghost; j++) {
        for (i=is-nghost; i<=ie+nghost; i++) {
            list = &((pG->GridLists)[k][j][i]);
            n += list->count;
        }
    }
    return n;
}

static void pack_ix3_tracers(GridS *pG)
{
    int is = pG->is, ie = pG->ie;
    int js = pG->js, je = pG->je;
    int ks = pG->ks, ke = pG->ke;
    int i,j,k,l,m;
    TracerListS *list;
    TracerS *tracer;
    TracerS *Next;
    
    l = 0;
    m = 0;
    k = ks-1;
    for (j=js-nghost; j<=je+nghost; j++) {
        for (i=is-nghost; i<=ie+nghost; i++) {
            list = &((pG->GridLists)[k][j][i]);
            send_buf_list0x3[m].count = list->count;
            m++;
            tracer = list->Head;
            while(tracer) {
                send_buf_mc0x3[l].x1 = (double)tracer->prop->id;
                send_buf_mc0x3[l].d_init = (double)tracer->prop->d_init;
                send_buf_mc0x3[l].i_init = (int)tracer->prop->i_init;
                send_buf_mc0x3[l].j_init = (int)tracer->prop->j_init;
                send_buf_mc0x3[l].k_init = (int)tracer->prop->k_init;
#ifdef STAR_PARTICLE
                send_buf_mc0x3[l].star_id = (int)tracer->prop->star_id;
#endif /* STAR_PARTICLE */
#ifdef VFTRACERS
                send_buf_mc0x3[l].x1 = (double)tracer->x1;
                send_buf_mc0x3[l].x2 = (double)tracer->x2;
                send_buf_mc0x3[l].x3 = (double)tracer->x3;
#endif /* VFTRACERS */
                l++;
                Next = tracer->Next;
                tracer = Next;
            }
            tracer = list->Head;
            while(tracer) {
                Next = tracer->Next;
                tracer_list_remove(list, tracer);
                if(tracer) free(tracer);
                tracer = Next;
            }
#ifdef DEBUG
            assert(list->count == 0);
#endif /* DEBUG */
        }
    }
#ifdef DEBUG
    assert(l == send_buf[0]);
#endif /* DEBUG */
}

static void pack_ox3_tracers(GridS *pG)
{
    int is = pG->is, ie = pG->ie;
    int js = pG->js, je = pG->je;
    int ks = pG->ks, ke = pG->ke;
    int i,j,k,l,m;
    TracerListS *list;
    TracerS *tracer;
    TracerS *Next;
    
    l = 0;
    m = 0;
    k = ke+1;
    for (j=js-nghost; j<=je+nghost; j++) {
        for (i=is-nghost; i<=ie+nghost; i++) {
            list = &((pG->GridLists)[k][j][i]);
            send_buf_list1x3[m].count = list->count;
            m++;
            tracer = list->Head;
            while(tracer) {
                send_buf_mc1x3[l].id = (double)tracer->prop->id;
                send_buf_mc1x3[l].d_init = (double)tracer->prop->d_init;
                send_buf_mc1x3[l].i_init = (int)tracer->prop->i_init;
                send_buf_mc1x3[l].j_init = (int)tracer->prop->j_init;
                send_buf_mc1x3[l].k_init = (int)tracer->prop->k_init;
#ifdef STAR_PARTICLE
                send_buf_mc1x3[l].star_id = (int)tracer->prop->star_id;
#endif /* STAR_PARTICLE */
#ifdef VFTRACERS
                send_buf_mc1x3[l].x1 = (double)tracer->x1;
                send_buf_mc1x3[l].x2 = (double)tracer->x2;
                send_buf_mc1x3[l].x3 = (double)tracer->x3;
#endif /* VFTRACERS */
                l++;
                Next = tracer->Next;
                tracer = Next;
            }
            tracer = list->Head;
            while(tracer) {
                Next = tracer->Next;
                tracer_list_remove(list, tracer);
                if(tracer) free(tracer);
                tracer = Next;
            }
#ifdef DEBUG
            assert(list->count == 0);
#endif /* DEBUG */
        }
    }
#ifdef DEBUG
    assert(l == send_buf[1]);
#endif /* DEBUG */
}

static int pack_ox3_mc(GridS *pG)
{
    int is = pG->is, ie = pG->ie;
    int js = pG->js, je = pG->je;
    int ks = pG->ks, ke = pG->ke;
    int i,j,k,n;
    TracerListS *list;
    
    n = 0;
    k = ke+1;
    for (j=js-nghost; j<=je+nghost; j++) {
        for (i=is-nghost; i<=ie+nghost; i++) {
            list = &((pG->GridLists)[k][j][i]);
            n += list->count;
        }
    }
    return n;
}

#endif /* MPI_PARALLEL */

static void periodic_ix1_mc(GridS *pG)
{
    int is = pG->is, ie = pG->ie;
    int js = pG->js, je = pG->je;
    int ks = pG->ks, ke = pG->ke;
    int i,j,k;
    TracerS *pnode;
    TracerListS *list, *newList;

    /* Move tracers from ghost zone into active zone */
    for (k=ks; k<=ke; k++) {
        for (j=js; j<=je; j++) {
            for (i=1; i<=nghost; i++) {
                /* Tracer list to move */
                list = &((pG->GridLists)[k][j][is-1]);
                /* List in active zone */
                newList = &((pG->GridLists)[k][j][ie-(i-1)]);
                pnode = list->Head;
                /* flag tracers for movement */
                while(pnode) {
                    pnode->newList = newList;
                    pnode = pnode->Next;
                }
                /* Changes x position of tracer corresponding to new list */
#ifdef VFTRACERS
                vf_newpos(pG, list, newList);
#endif /* VFTRACERS */
                Tracerlist_sweep_bc(list);
            }
        }
    }
    return;
}

static void periodic_ox1_mc(GridS *pG)
{
    int is = pG->is, ie = pG->ie;
    int js = pG->js, je = pG->je;
    int ks = pG->ks, ke = pG->ke;
    int i,j,k;
    TracerS *pnode;
    TracerListS *list, *newList;
    
    for (k=ks; k<=ke; k++) {
        for (j=js; j<=je; j++) {
            for (i=1; i<=nghost; i++) {
                /* Tracer list to move */
                list = &((pG->GridLists)[k][j][ie+i]);
                /* List in active zone */
                newList = &((pG->GridLists)[k][j][is+(i-1)]);
                pnode = list->Head;
                
                /* flag tracers for movement */
                while(pnode) {
                    pnode->newList = newList;
                    pnode = pnode->Next;
                }
#ifdef VFTRACERS
                vf_newpos(pG, list, newList);
#endif /* VFTRACERS */
                /* Move tracers */
                Tracerlist_sweep_bc(list);
            }
        }
    }
    return;
}

static void periodic_ix2_mc(GridS *pG)
{
    int is = pG->is, ie = pG->ie;
    int js = pG->js, je = pG->je;
    int ks = pG->ks, ke = pG->ke;
    int i,j,k;
    TracerS *pnode;
    TracerListS *list, *newList;
    
    for (k=ks; k<=ke; k++) {
        for (j=1; j<=nghost; j++) {
            for (i=is-nghost; i<=ie+nghost; i++) {
                /* Tracer list to move (from ghost zone) */
                list = &((pG->GridLists)[k][js-j][i]);
                /* List in active zone */
                newList = &((pG->GridLists)[k][je-(j-1)][i]);
                pnode = list->Head;
                /* flag tracers for movement */
                while(pnode) {
                    pnode->newList = newList;
                    pnode = pnode->Next;
                }
#ifdef VFTRACERS
                vf_newpos(pG, list, newList);
#endif /* VFTRACERS */
                /* Move tracers */
                Tracerlist_sweep_bc(list);
            }
        }
    }
    return;
}

static void periodic_ox2_mc(GridS *pG)
{
    int is = pG->is, ie = pG->ie;
    int js = pG->js, je = pG->je;
    int ks = pG->ks, ke = pG->ke;
    int i,j,k;
    TracerS *pnode;
    TracerListS *list, *newList;
    
    for (k=ks; k<=ke; k++) {
        for (j=1; j<=nghost; j++) {
            for (i=is-nghost; i<=ie+nghost; i++) {
                
                /* Tracer list to move (from ghost zone) */
                list = &((pG->GridLists)[k][je+j][i]);
                /* List in active zone */
                newList = &((pG->GridLists)[k][js+(j-1)][i]);
                pnode = list->Head;
                
                /* flag tracers for movement */
                while(pnode) {
                    pnode->newList = newList;
                    pnode = pnode->Next;
                }
#ifdef VFTRACERS
                vf_newpos(pG, list, newList);
#endif /* VFTRACERS */
                /* Move tracers */
                Tracerlist_sweep_bc(list);
            }
        }
    }
    return;
}

static void periodic_ix3_mc(GridS *pG)
{
    int is = pG->is, ie = pG->ie;
    int js = pG->js, je = pG->je;
    int ks = pG->ks, ke = pG->ke;
    int i,j,k;
    TracerS *pnode;
    TracerListS *list, *newList;
    
    for (k=1; k<=nghost; k++) {
        for (j=js-nghost; j<=je+nghost; j++) {
            for (i=is-nghost; i<=ie+nghost; i++) {
                /* Tracer list to move (from ghost zone) */
                list = &((pG->GridLists)[ks-k][j][i]);
                /* List in active zone */
                newList = &((pG->GridLists)[ke-(k-1)][j][i]);
                pnode = list->Head;
                
                /* flag tracers for movement */
                while(pnode) {
                    pnode->newList = newList;
                    pnode = pnode->Next;
                }
#ifdef VFTRACERS
                vf_newpos(pG, list, newList);
#endif /* VFTRACERS */
                /* Move tracers */
                Tracerlist_sweep_bc(list);
            }
        }
    }
    return;
}

static void periodic_ox3_mc(GridS *pG)
{
    int is = pG->is, ie = pG->ie;
    int js = pG->js, je = pG->je;
    int ks = pG->ks, ke = pG->ke;
    int i,j,k;
    TracerS *pnode;
    TracerListS *list, *newList;
    
    for (k=1; k<=nghost; k++) {
        for (j=js-nghost; j<=je+nghost; j++) {
            for (i=is-nghost; i<=ie+nghost; i++) {
                /* Tracer list to move (from ghost zone) */
                list = &((pG->GridLists)[ke+k][j][i]);
                /* List in active zone */
                newList = &((pG->GridLists)[ks+(k-1)][j][i]);
                pnode = list->Head;
                /* flag tracers for movement */
                while(pnode) {
                    pnode->newList = newList;
                    pnode = pnode->Next;
                }
#ifdef VFTRACERS
                vf_newpos(pG, list, newList);
#endif /* VFTRACERS */
                /* Move tracers */
                Tracerlist_sweep_bc(list);
            }
        }
    }
    return;
}

#ifdef MPI_PARALLEL

static void unpack_ix1_mc(GridS *pG)
{
    int is = pG->is, ie = pG->ie;
    int js = pG->js, je = pG->je;
    int ks = pG->ks, ke = pG->ke;
    int i,j,k,l,m,n, count;
    TracerListS *list;
    TracerS *tracer;
    Tracer_MPI *mpi;
    TracerPropS *prop;
    
    l = 0;
    m = 0;
    for (k=ks; k<=ke; k++){
        for (j=js; j<=je; j++){
            i = is;
            list = &((pG->GridLists)[k][j][i]);
            count = recv_buf_list0x1[l].count;
            for (n=1; n<=count; n++) {
                tracer = init_tracer();
                prop = (TracerPropS *)malloc(sizeof(TracerPropS));
                prop->id = (double)recv_buf_mc0x1[m].id;
                prop->d_init = (double)recv_buf_mc0x1[m].d_init;
                prop->i_init = (int)recv_buf_mc0x1[m].i_init;
                prop->j_init = (int)recv_buf_mc0x1[m].j_init;
                prop->k_init = (int)recv_buf_mc0x1[m].k_init;
#ifdef STAR_PARTICLE
                prop->star_id = (int)recv_buf_mc0x1[m].star_id;
#endif /* STAR_PARTICLE */
#ifdef VFTRACERS
                tracer->x1 = (double)recv_buf_mc0x1[m].x1;
                tracer->x2 = (double)recv_buf_mc0x1[m].x2;
                tracer->x3 = (double)recv_buf_mc0x1[m].x3;
#endif /* VFTRACERS */
                tracer->prop = prop;
                Tracerlist_add(list,tracer);
                m++;
            }
#ifdef MCTRACERS
            list->Rmass = pG->U[k][j][i].d;
#endif /* MCTRACERS */
            l++;
        }
    }
}

static void unpack_ox1_mc(GridS *pG)
{
    int is = pG->is, ie = pG->ie;
    int js = pG->js, je = pG->je;
    int ks = pG->ks, ke = pG->ke;
    int i,j,k,l,m,n,count;
    TracerListS *list;
    TracerS *tracer;
    Tracer_MPI *mpi;
    TracerPropS *prop;
    
    l = 0;
    m = 0;
    for (k=ks; k<=ke; k++){
        for (j=js; j<=je; j++){
            i = ie;
            list = &((pG->GridLists)[k][j][i]);
            count = recv_buf_list1x1[l].count;
            for (n=1; n<=count; n++) {
                tracer = init_tracer();
                prop = (TracerPropS *)malloc(sizeof(TracerPropS));
#ifdef DEBUG
                assert(recv_buf_mc1x1 != NULL);
#endif /* DEBUG */
                prop->id = (double)recv_buf_mc1x1[m].id;
                prop->d_init = (double)recv_buf_mc1x1[m].d_init;
                prop->i_init = (int)recv_buf_mc1x1[m].i_init;
                prop->j_init = (int)recv_buf_mc1x1[m].j_init;
                prop->k_init = (int)recv_buf_mc1x1[m].k_init;
#ifdef STAR_PARTICLE
                prop->star_id = (int)recv_buf_mc1x1[m].star_id;
#endif /* STAR_PARTICLE */
#ifdef VFTRACERS
                tracer->x1 = (double)recv_buf_mc1x1[m].x1;
                tracer->x2 = (double)recv_buf_mc1x1[m].x2;
                tracer->x3 = (double)recv_buf_mc1x1[m].x3;
#endif /* VFTRACERS */
                tracer->prop = prop;
                Tracerlist_add(list,tracer);
                m++;
            }
#ifdef MCTRACERS
            list->Rmass = pG->U[k][j][i].d;
#endif /* MCTRACERS */
            l++;
        }
    }
}

static void unpack_ix2_mc(GridS *pG)
{
    int is = pG->is, ie = pG->ie;
    int js = pG->js, je = pG->je;
    int ks = pG->ks, ke = pG->ke;
    int i,j,k,l,m,n,count;
    TracerListS *list;
    TracerS *tracer;
    Tracer_MPI *mpi;
    TracerPropS *prop;
    l = 0;
    m = 0;
    for (k=ks; k<=ke; k++) {
        j = js;
        for (i=is-nghost; i<=ie+nghost; i++) {
            list = &((pG->GridLists)[k][j][i]);
            count = recv_buf_list0x2[l].count;
            for (n=1; n<=count; n++) {
                tracer = init_tracer();
                prop = (TracerPropS *)malloc(sizeof(TracerPropS));
                prop->id = (double)recv_buf_mc0x2[m].id;
                prop->d_init = (double)recv_buf_mc0x2[m].d_init;
                prop->i_init = (int)recv_buf_mc0x2[m].i_init;
                prop->j_init = (int)recv_buf_mc0x2[m].j_init;
                prop->k_init = (int)recv_buf_mc0x2[m].k_init;
#ifdef STAR_PARTICLE
                prop->star_id = (int)recv_buf_mc0x2[m].star_id;
#endif /* STAR_PARTICLE */
#ifdef VFTRACERS
                tracer->x1 = (double)recv_buf_mc0x2[m].x1;
                tracer->x2 = (double)recv_buf_mc0x2[m].x2;
                tracer->x3 = (double)recv_buf_mc0x2[m].x3;
#endif /* VFTRACERS */
                tracer->prop = prop;
                Tracerlist_add(list,tracer);
                m++;
            }
#ifdef MCTRACERS
            list->Rmass = pG->U[k][j][i].d;
#endif /* MCTRACERS */
            l++;
        }
    }
}

static void unpack_ox2_mc(GridS *pG)
{
    int is = pG->is, ie = pG->ie;
    int je = pG->je, js = pG->js;
    int ks = pG->ks, ke = pG->ke;
    int i,j,k,l,m,n,count;
    TracerListS *list;
    TracerS *tracer;
    Tracer_MPI *mpi;
    TracerPropS *prop;
    l = 0;
    m = 0;
    for (k=ks; k<=ke; k++){
        j = je;
        for (i=is-nghost; i<=ie+nghost; i++){
            list = &((pG->GridLists)[k][j][i]);
            count = recv_buf_list1x2[l].count;
            for (n=1; n<=count; n++) {
                tracer = init_tracer();
                prop = (TracerPropS *)malloc(sizeof(TracerPropS));
                prop->id = (double)recv_buf_mc1x2[m].id;
                prop->d_init = (double)recv_buf_mc1x2[m].d_init;
                prop->i_init = (int)recv_buf_mc1x2[m].i_init;
                prop->j_init = (int)recv_buf_mc1x2[m].j_init;
                prop->k_init = (int)recv_buf_mc1x2[m].k_init;
#ifdef STAR_PARTICLE
                prop->star_id = (int)recv_buf_mc1x2[m].star_id;
#endif /* STAR_PARTICLE */
#ifdef VFTRACERS
                tracer->x1 = (double)recv_buf_mc1x2[m].x1;
                tracer->x2 = (double)recv_buf_mc1x2[m].x2;
                tracer->x3 = (double)recv_buf_mc1x2[m].x3;
#endif /* VFTRACERS */
                tracer->prop = prop;
                Tracerlist_add(list,tracer);
                m++;
            }
#ifdef MCTRACERS
            list->Rmass = pG->U[k][j][i].d;
#endif /* MCTRACERS */
            l++;
        }
    }
}

static void unpack_ix3_mc(GridS *pG)
{
    int is = pG->is, ie = pG->ie;
    int js = pG->js, je = pG->je;
    int ks = pG->ks, ke = pG->ke;
    int i,j,k,l,m,n,count;
    TracerListS *list;
    TracerS *tracer;
    Tracer_MPI *mpi;
    TracerPropS *prop;
    l = 0;
    m = 0;
    k = ks;
    for (j=js-nghost; j<=je+nghost; j++) {
        for (i=is-nghost; i<=ie+nghost; i++) {
            list = &((pG->GridLists)[k][j][i]);
            count = recv_buf_list0x3[l].count;
            for (n=1; n<=count; n++) {
                tracer = init_tracer();
                prop = (TracerPropS *)malloc(sizeof(TracerPropS));
                prop->id = (double)recv_buf_mc0x3[m].id;
                prop->d_init = (double)recv_buf_mc0x3[m].d_init;
                prop->i_init = (int)recv_buf_mc0x3[m].i_init;
                prop->j_init = (int)recv_buf_mc0x3[m].j_init;
                prop->k_init = (int)recv_buf_mc0x3[m].k_init;
#ifdef STAR_PARTICLE
                prop->star_id = (int)recv_buf_mc0x3[m].star_id;
#endif /* STAR_PARTICLE */
#ifdef VFTRACERS
                tracer->x1 = (double)recv_buf_mc0x3[m].x1;
                tracer->x2 = (double)recv_buf_mc0x3[m].x2;
                tracer->x3 = (double)recv_buf_mc0x3[m].x3;
#endif /* VFtRACERS */
                tracer->prop = prop;
                Tracerlist_add(list,tracer);
                m++;
            }
#ifdef MCTRACERS
            list->Rmass = pG->U[k][j][i].d;
#endif /* MCTRACERS */
            l++;
        }
    }
}

static void unpack_ox3_mc(GridS *pG)
{
    int is = pG->is, ie = pG->ie;
    int js = pG->js, je = pG->je;
    int ke = pG->ke, ks = pG->ks;
    int i,j,k,l,m,n,count;
    TracerListS *list;
    TracerS *tracer;
    Tracer_MPI *mpi;
    TracerPropS *prop;
    l = 0;
    m = 0;
    k = ke;
    for (j=js-nghost; j<=je+nghost; j++) {
        for (i=is-nghost; i<=ie+nghost; i++) {
            list = &((pG->GridLists)[k][j][i]);
            count = recv_buf_list1x3[l].count;
            for (n=1; n<=count; n++) {
                tracer = init_tracer();
                prop = (TracerPropS *)malloc(sizeof(TracerPropS));
                prop->id = (double)recv_buf_mc1x3[m].id;
                prop->d_init = (double)recv_buf_mc1x3[m].d_init;
                prop->i_init = (int)recv_buf_mc1x3[m].i_init;
                prop->j_init = (int)recv_buf_mc1x3[m].j_init;
                prop->k_init = (int)recv_buf_mc1x3[m].k_init;
#ifdef STAR_PARTICLE
                prop->star_id = (int)recv_buf_mc1x3[m].star_id;
#endif /* STAR_PARTICLE */
#ifdef VFTRACERS
                tracer->x1 = (double)recv_buf_mc1x3[m].x1;
                tracer->x2 = (double)recv_buf_mc1x3[m].x2;
                tracer->x3 = (double)recv_buf_mc1x3[m].x3;
#endif /* VFTRACERS */
                tracer->prop = prop;
                Tracerlist_add(list, tracer);
                m++;
            }
#ifdef MCTRACERS
            list->Rmass = pG->U[k][j][i].d;
#endif /* MCTRACERS */
            l++;
        }
    }
}

#endif /* MPI_PARALLEL */

#endif /* VFTRACERS */
