/*----------------------------------------------------------------------------*/
/*  \brief Function to write dumps of mc tracer "history" in a
/*  formatted table. */
/*  This functio is currently not being called and probably has not been tested */ 
void mctracer_out(MeshS *pM, OutputS *pOut)
{
    GridS *pG;
    DomainS *pD;
    
#ifdef STAR_PARTICLE
    StarParListS *pGstars = NULL;
    StarParS *pStar = NULL;
#endif
    
    int nl,nd;
    double scal[NSCAL + MAX_USR_H_COUNT];
    FILE *pfile;
    char *fname,*plev=NULL,*pdom=NULL,*pdir=NULL,fmt[80];
    char levstr[8],domstr[8],dirstr[20];
    int n, total_hst_cnt, mhst, myID_Comm_Domain=1;
#ifdef MPI_PARALLEL
    int ierr;
#endif
    
    total_hst_cnt = NSCAL + usr_hst_cnt;
    
    /* Add a white space to the format */
    if(pOut->dat_fmt == NULL){
        sprintf(fmt," %%14.6e"); /* Use a default format */
    }
    else{
        sprintf(fmt," %s",pOut->dat_fmt);
    }
    
    /* Loop over all Domains in Mesh, and output Grid data */
    
    for (nl=0; nl<(pM->NLevels); nl++){
        for (nd=0; nd<(pM->DomainsPerLevel[nl]); nd++){
            if (pM->Domain[nl][nd].Grid != NULL){
                pG = pM->Domain[nl][nd].Grid;
                pD = (DomainS*)&(pM->Domain[nl][nd]);
                
                /* Only the parent (rank=0) process computes the average and writes output.
                 * For single-processor jobs, myID_Comm_world is always zero. */
                
#ifdef MPI_PARALLEL
                ierr = MPI_Comm_rank(pD->Comm_Domain, &myID_Comm_Domain);
#endif
                if((myID_Comm_Domain==0) || (myID_Comm_world==0)){  /* I'm the parent */
                    
                    pGtracers = pG->GridLists;
                    for (k=kl; k<=kr; k++) {
                        for (j=jl; j<=jr; j++) {
                            for (i=il; i<=ir; i++) {
                                list = &((pG->GridLists)[k][j][i]);
                                tracer = list->Head;
                                while(tracer) {
                                    mhst = 0;
                                    scal[mhst++] = pG->time;
                                    scal[mhst++] = pG->U[j][k][i].d;
                                    scal[mhst++] = tracer->d_init;
                                    scal[mhst++] = tracer->i_init;
                                    scal[mhst++] = tracer->j_init;
                                    scal[mhst++] = tracer->k_init;

                                }
                    while(pGstars) {
                        pStar = &(pGstars->starpar);
                        
                        mhst = 0;
                        scal[mhst++] = pG->time;
                        scal[mhst++] = pStar->m;
                        scal[mhst++] = pStar->x1;
                        scal[mhst++] = pStar->x2;
                        scal[mhst++] = pStar->x3;
                        scal[mhst++] = pStar->v1;
                        scal[mhst++] = pStar->v2;
                        scal[mhst++] = pStar->v3;
                        scal[mhst++] = pStar->age;
                        scal[mhst++] = pStar->mdot;
                        scal[mhst++] = pStar->merge_history;
                        scal[mhst++] = pStar->id;
                        /* Calculate the user defined history variables */
                        for(n=0; n<usr_hst_cnt; n++){
                            scal[mhst++] = (*phst_fun[n])(pG, pStar);
                        }
                        
                        /* Create filename and open file.  History files are always written in lev#
                         * directories of root process (rank=0 in MPI_COMM_WORLD) */
#ifdef MPI_PARALLEL
                        if (nl>0) {
                            plev = &levstr[0];
                            sprintf(plev,"lev%d",nl);
                            pdir = &dirstr[0];
                            sprintf(pdir,"../id0/lev%d",nl);
                        }
#else
                        if (nl>0) {
                            plev = &levstr[0];
                            sprintf(plev,"lev%d",nl);
                            pdir = &dirstr[0];
                            sprintf(pdir,"lev%d",nl);
                        }
#endif
                        if (nd>0) {
                            pdom = &domstr[0];
                            sprintf(pdom,"dom%d",nd);
                        }
                        
                        fname = ath_fname(pdir,pM->outfilename,plev,pdom,num_digit,pStar->id,NULL,"star");
                        if(fname == NULL){
                            ath_perr(-1,"[dump_starpar_history]: Unable to create history filename\n");
                        }
                        if(pStar->age == 0.0) pfile = fopen(fname,"w");
                        else pfile = fopen(fname,"a");
                        if(pfile == NULL){
                            ath_perr(-1,"[dump_starparhistory]: Unable to open the history file\n");
                        }
                        free(fname);
                        
                        /* Write out column headers, but only for first dump */
                        mhst = 1;
                        if(pStar->age == 0){
                            fprintf(pfile,"# Athena star particle history dump for level=%i domain=%i id=%i\n",nl,nd,pStar->id);
                            fprintf(pfile,"#   [%i]=time   ",mhst++);
                            fprintf(pfile,"   [%i]=mass    ",mhst++);
                            fprintf(pfile,"   [%i]=x1      ",mhst++);
                            fprintf(pfile,"   [%i]=x2      ",mhst++);
                            fprintf(pfile,"   [%i]=x3      ",mhst++);
                            fprintf(pfile,"   [%i]=v1      ",mhst++);
                            fprintf(pfile,"   [%i]=v2      ",mhst++);
                            fprintf(pfile,"   [%i]=v3      ",mhst++);
                            fprintf(pfile,"   [%i]=age     ",mhst++);
                            fprintf(pfile,"   [%i]=mdot    ",mhst++);
                            fprintf(pfile,"   [%i]=mhist   ",mhst++);
                            fprintf(pfile,"   [%i]=id   ",mhst++);
                            for(n=0; n<usr_hst_cnt; n++){
                                fprintf(pfile,"  [%i]=%s",mhst++,usr_label[n]);
                            }
                            fprintf(pfile,"\n#\n");
                        }
                        
                        /* Write out data, and close file */
                        for (n=0; n<total_hst_cnt; n++) {
                            fprintf(pfile,fmt,scal[n]);
                        }
                        fprintf(pfile,"\n");
                        
                        fclose(pfile);
                        /* Move to next star particle */
                        pGstars = pGstars->next;
                    }
                }
            }
        }
    }
    
    return;
}
