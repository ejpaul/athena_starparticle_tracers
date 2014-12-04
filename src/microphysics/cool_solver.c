#include "../copyright.h"
/*==============================================================================
 * FILE: cool_solver.c
 *
 * PURPOSE: Add cooling to energy.
 *   Koyama & Inutsuka (2002)'s cooling function is adopted.  Assuming local
 *   density is a constant within cooling time scale.  Temperature is updated
 *   by using either integration using Simpson's rule or fully implicit method.
 *   Timestep is limited by this routine to ensure stability.
 *
 *   originally written by R. Piontek  
 *   modified in C by Chang-Goo Kim at 2006-11-10 
 *
 * CONTAINS PUBLIC FUNCTIONS:
 *   cooling_solver() 
 *   cooling_solver_init() - allocates memory needed
 *   cooling_solver_destruct() - frees memory used
 *
 * CONTAINS PRIVATE FUNCTIONS:
 * Real Lam(Real temp);
 * Real dLam(Real temp);
 * Real temp_next(Real t_old, Real nden, Real *dt, int update);
 *============================================================================*/

#include <math.h>
#include "../defs.h"
#include "../athena.h"
#include "../globals.h"
#include "../prototypes.h"
#include "prototypes.h"

#ifdef OPERATOR_SPLIT_COOLING
#ifdef BAROTROPIC
#error : cooling solver does not work with Isothermal EOS
#endif

/* uncomment to use sub_cycle */
//#define SUB_CYCLE
//#define SIMPSON
#ifdef STAR_PARTICLE
#define STAR_PARTICLE_MASK
#endif

#define SD93
//#define GF12
Real Lam(Real temp);
Real dLam(Real temp);
Real temp_next(Real t_old, Real nden, Real heat, Real *dt, int update);

#ifdef SUBCYCLE
static Real ***dt_sub=NULL;
#endif

#ifdef STAR_PARTICLE_MASK
static int ***mask=NULL;
void set_starpar_mask(GridS *pG);
#endif

// constraint for maximum temperature change at each time step and tolerance for NR convergence
static Real dtempmax=1.00, toler=0.01, dfloor=1.e-5; 
static Real unitT, unitC; // units for temperature and cooling rate

//=======================================================================
void cooling_solver(GridS *pG)
{

  int i, is = pG->is, ie = pG->ie;
  int j, js = pG->js, je = pG->je;
  int k, ks = pG->ks, ke = pG->ke;

#ifdef MPI_PARALLEL
  int err;
#endif

  int cycle, ncycle=1, ncycle_max=0, cycle_max=1;
  int reset=0;

  Real temp, t_old, nden;
  Real my_dt,dt;
  Real my_dt_min=HUGE_NUMBER;
  Real t_bt, dt_off;
  Real x1,x2,x3;

  Real Tmin = 10.;
  Real heat = pG->heat0*pG->heat_ratio;
  Real tcool;

  Real eint,v1,v2,v3;
  int Ncells,kk,jj,ii;

  ConsS U;


#ifdef STAR_PARTICLE_MASK
  set_starpar_mask(pG);
#endif
 
  for(k=ks;k<=ke;k++){ 
    for(j=js;j<=je;j++){ 
      for(i=is;i<=ie;i++){ 
#ifdef STAR_PARTICLE_MASK
      if(mask[k][j][i]){
#endif
        my_dt=pG->dt;

        U=pG->U[k][j][i];

        t_old = U.E - (0.5/U.d)*(SQR(U.M1)+SQR(U.M2)+SQR(U.M3));
#ifdef MHD
        t_old -= (0.5)*(SQR(U.B1c)+SQR(U.B2c)+SQR(U.B3c));
#endif
        t_old *= (Gamma_1/U.d);
        t_old = MAX(unitT*t_old,Tmin); // in unit of Kelvin
        //nden  = MAX(U.d,dfloor);  // in unit of mbar
        nden  = U.d;  // in unit of mbar

//	tcool = t_old/nden/(Lam(t_old)*unitC);
//	my_dt = MIN(2.0*tcool,pG->dt);
	
#ifdef STAR_PARTICLE
        heat = pG->heat0*pG->sfr_ratio[k][j][i];
#endif
        temp=temp_next(t_old,nden,heat,&my_dt,0);
        if( temp != temp){
#ifdef STAR_PARTICLE_MASK
          cc_pos(pG,i,j,k,&x1,&x2,&x3);
          starpar_printlist(-1, pG);
          printf("[cool_solver.c] Temperature is NaN: id=%d; x,y,z=%g %g %g\n",myID_Comm_world,x1,x2,x3);
#endif
          printf("[cool_solver.c] Temperature is NaN: id=%d; i,j,k=%d %d %d; t_old=%g temp=%g nden=%g; dt=%g %g tcool=%g\n",myID_Comm_world,i,j,k,t_old,temp,nden,pG->dt,my_dt,tcool);
          ath_error("[cool_solver.c] Temperature is NaN: id=%d; i,j,k=%d %d %d; t_old=%g temp=%g nden=%g; dt=%g %g tcool=%g\n",myID_Comm_world,i,j,k,t_old,temp,nden,pG->dt,my_dt,tcool);
        }
//        if(myID_Comm_world == 0) printf("[cool_solver.c] Temperature is NaN: id=%d; i,j,k=%d %d %d; t_old=%g temp=%g nden=%g; dt=%g %g tcool=%g\n",myID_Comm_world,i,j,k,t_old,temp,nden,pG->dt,my_dt,tcool);
#ifdef SUB_CYCLE
        dt_sub[k-ks][j-js][i-is]=my_dt;
#endif
        my_dt_min=MIN(my_dt_min,my_dt);
#ifdef STAR_PARTICLE_MASK
      } /* if mask */
#endif
      }
    }
  }

/* Set a time step which guarantees NR convergency and 
 * small relative change of temperature for whole grids
 * and recalculate temperature change with the time step */

#ifdef MPI_PARALLEL
  my_dt=my_dt_min;
  err = MPI_Allreduce(&my_dt, &dt, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  if(err) ath_error("[sync_dt]: MPI_Allreduce returned error code %d\n",err);
#else
  dt = my_dt_min;
#endif // MPI_PARALLEL
//  printf("[cool_solver.c] dt is changed from %g to %g at t=%g, id=%d\n",pG->dt,dt,pG->time,myID_Comm_world);


#ifndef SUB_CYCLE
  pG->dt=dt;
#endif
/* If NR converged and relative change of temperature is not too large for entire grid,
 * update energy with temperature floor Tmin=10K.
 */
  for(k=ks;k<=ke;k++){ 
    for(j=js;j<=je;j++){ 
      for(i=is;i<=ie;i++){ 
#ifdef STAR_PARTICLE_MASK
      if(mask[k][j][i]){
#endif
//        reset_dt:

        U=pG->U[k][j][i];

        t_old = U.E - (0.5/U.d)*(SQR(U.M1)+SQR(U.M2)+SQR(U.M3));
#ifdef MHD
        t_old -= (0.5)*(SQR(U.B1c)+SQR(U.B2c)+SQR(U.B3c));
#endif
        t_old *= (Gamma_1/U.d);
        t_old = MAX(unitT*t_old,Tmin); // in unit of Kelvin
        //nden  = MAX(U.d,dfloor);  // in unit of mbar
        nden  = U.d;  // in unit of mbar

#ifdef SUB_CYCLE
          reset=0;
          dt = dt_sub[k-ks][j-js][i-is];
          ncycle = (int)(pG->dt/dt);
          ncycle_max = MAX(ncycle_max,ncycle);
          cycle_max=ncycle;
#endif
        temp = t_old;

        for(cycle = 0; cycle < cycle_max; cycle++){ 
          my_dt=dt;
#ifdef STAR_PARTICLE
          heat = pG->heat0*pG->sfr_ratio[k][j][i];
#endif
          temp=temp_next(t_old,nden,heat,&dt,1);
          if( my_dt != dt ) {
/*
#ifdef SUB_CYCLE
            dt_sub[k-ks][j-js][i-is] = dt; 
#endif
            reset++;
            if(reset > 3) ath_error("[cool_solver.c] time step becomes too small reset=%d at t=%g\nfor i,j,k=%d %d %d, t_old=%g nden=%g",reset,pG->time,i,j,k,t_old,U.d);
            cycle_max = (int)(my_dt/dt);
            ath_perr(-1,"[cool_solver.c] dt is changed from %g to %g ncycle=%d reset=%d\n",my_dt,dt,cycle_max,reset);
            goto reset_dt;
*/
            ath_perr(-1,"[cool_solver.c] Solution didn't converge. Simply take average for [k][j][i]=[%d][%d][%d] of id=%d\n",k,j,i,myID_Comm_world);
          }
/* If NR converged, check relative change of temperature. */
          if( temp != temp) ath_error("[cool_solver.c] Temperature is NaN cycle=%d ncycle=%d\n",cycle,ncycle);
          t_old=temp;
        }

        if(dt == pG->dt) {
          temp=MAX(temp,Tmin)/unitT;
//        pG->U[k][j][i].d = MAX(U.d,dfloor);  // in unit of mbar
          pG->U[k][j][i].E = temp*pG->U[k][j][i].d/Gamma_1;
          pG->U[k][j][i].E += (0.5/pG->U[k][j][i].d)*
             (SQR(pG->U[k][j][i].M1)+SQR(pG->U[k][j][i].M2)+SQR(pG->U[k][j][i].M3));
#ifdef MHD
          pG->U[k][j][i].E += (0.5)*
             (SQR(pG->U[k][j][i].B1c)+SQR(pG->U[k][j][i].B2c)+SQR(pG->U[k][j][i].B3c));
#endif
        } else {
#ifdef STAR_PARTICLE_MASK
// mark bad cell to mask array
          mask[k][j][i] = -1;      
#endif
        }
        dt=pG->dt;
        reset=0;
#ifdef STAR_PARTICLE_MASK
      } /* if mask */
#endif
      }
    }
  }

  if (ncycle_max > 1 ) ath_perr(-1,"subcycling %d times at %g\n",ncycle_max,pG->time);


#ifdef STAR_PARTICLE_MASK
/* correct bad cells by taking averages */
  for(k=ks;k<=ke;k++){ 
    for(j=js;j<=je;j++){ 
      for(i=is;i<=ie;i++){ 
        if(mask[k][j][i] == -1){
          eint +=(U.E - (0.5/U.d)*(SQR(U.M1)+SQR(U.M2)+SQR(U.M3)));
#ifdef MHD
          eint -= (0.5)*(SQR(U.B1c)+SQR(U.B2c)+SQR(U.B3c));
#endif
          ath_perr(-1,"[cool_solver.c] Bad cell is detected: [k][j][i]=[%d][%d][%d], id=%d\n",
                  k,j,i,myID_Comm_world);
          ath_perr(-1,"[cool_solver.c] Original: n = %g, eint= %g\n",pG->U[k][j][i].d,eint);

          eint = 0.0;
          v1 = 0.0;
          v2 = 0.0;
          v3 = 0.0;
          nden = 0.0;
          Ncells=0;
          for(kk=k-1;kk<=k+1;kk++){ 
          for(jj=j-1;jj<=j+1;jj++){ 
          for(ii=i-1;ii<=i+1;ii++){ 
            if(mask[kk][jj][ii] == 1){
              U=pG->U[kk][jj][ii];
              eint +=(U.E - (0.5/U.d)*(SQR(U.M1)+SQR(U.M2)+SQR(U.M3)));
#ifdef MHD
              eint -= (0.5)*(SQR(U.B1c)+SQR(U.B2c)+SQR(U.B3c));
#endif
              v1 += U.M1;
              v2 += U.M2;
              v3 += U.M3;
              nden += U.d; 
              Ncells++;
            }
          }}}
          eint = eint/(Real)Ncells;
          v1 = v1/nden;
          v2 = v2/nden;
          v3 = v3/nden;
          nden = nden/(Real)Ncells;
          ath_perr(-1,"[cool_solver.c] Average: n = %g, eint= %g, Ncells=%d\n",nden,eint,Ncells);
/*
          pG->U[k][j][i].d = nden;
          pG->U[k][j][i].M1 = nden*v1;
          pG->U[k][j][i].M2 = nden*v2;
          pG->U[k][j][i].M3 = nden*v3;
*/
          pG->U[k][j][i].E = eint + 0.5/pG->U[k][j][i].d*(SQR(pG->U[k][j][i].M1)+SQR(pG->U[k][j][i].M2)+SQR(pG->U[k][j][i].M3));
        }
      }
    }
  }
#endif

  return;
}

/*----------------------------------------------------------------------------*/
/* cooling_solver_init: 
*/

void cooling_solver_init(MeshS *pM){
  int nl,nd,size1=1,size2=1,size3=1,Nx1,Nx2,Nx3;
  UnitS units;
  ConstS consts;

/* Cycle over all Grids on this processor to find maximum Nx1, Nx2, Nx3 */
  for (nl=0; nl<(pM->NLevels); nl++){
    for (nd=0; nd<(pM->DomainsPerLevel[nl]); nd++){
      if (pM->Domain[nl][nd].Grid != NULL) {
          size1 = MAX(size1,pM->Domain[nl][nd].Grid->Nx[0]);
          size2 = MAX(size2,pM->Domain[nl][nd].Grid->Nx[1]);
          size3 = MAX(size3,pM->Domain[nl][nd].Grid->Nx[2]);
          units = pM->Domain[nl][nd].Grid->units;
      }
    }
  }

  init_consts(&consts);
  unitT = SQR(units.Vcode)*units.Dcode/1.1/consts.kB;
  unitC = units.Lcode/units.Vcode*Gamma_1/consts.kB/1.1; // unit for cooling rate

#ifdef SUB_CYCLE
  if ((dt_sub = (Real***)calloc_3d_array(size3,size2,size1, sizeof(Real))) == NULL)
    goto on_error;
#endif

  size1 = size1 + 2*nghost;
  size2 = size2 + 2*nghost;
  size3 = size3 + 2*nghost;

#ifdef STAR_PARTICLE_MASK
  if ((mask = (int***)calloc_3d_array(size3,size2,size1, sizeof(int))) == NULL)
    goto on_error;
#endif

  return;

  on_error:
  cooling_solver_destruct();
  ath_error("[cooling_solver_init]: malloc returned a NULL pointer\n");

}

void cooling_solver_destruct(void)
{
#ifdef STAR_PARTICLE_MASK
  if (mask != NULL) free_3d_array(mask);
#endif

#ifdef SUB_CYCLE
  if (dt_sub != NULL) free_3d_array(dt_sub);
#endif
  return;
}

/*----------------------------------------------------------------------------*/
/* PRIVATE FUCNTION                                                           */
/*----------------------------------------------------------------------------*/

#ifdef STAR_PARTICLE_MASK
void set_starpar_mask(GridS *pG){
  StarParListS *pList=NULL;
  StarParS *pStar=NULL;
  int i, j, k;
  int ip, jp, kp;
  int is=pG->is, ie=pG->ie, js=pG->js, je=pG->je, ks=pG->ks, ke=pG->ke;

  for(k=ks-nghost;k<=ke+nghost;k++){
  for(j=js-nghost;j<=je+nghost;j++){
  for(i=is-nghost;i<=ie+nghost;i++){
    mask[k][j][i]=1;
    if(pG->sfr_ratio[k][j][i] != 1) mask[k][j][i]=0;
  }}}

/*
  pList = pG->Gstars_fda;
  while (pList) {
    pStar = &(pList->starpar);
    
    cc_ijk(pG,pStar->x1,pStar->x2,pStar->x3,&ip,&jp,&kp);

    if ((ip >= is-1) && (ip <= ie+1) &&
        (jp >= js-1) && (jp <= je+1) &&
        (kp >= ks-1) && (kp <= ke+1)) {
      for(k=kp-1;k<=kp+1;k++){
      for(j=jp-1;j<=jp+1;j++){
      for(i=ip-1;i<=ip+1;i++){
        mask[k][j][i]=0;
      }}}
    }
    pList = pList->next;
  }
*/
}
#endif
/*----------------------------------------------------------------------------*/
/* temp_next:
*/


Real temp_next(Real t_old, Real nden, Real heat, Real *dt, int update)
{
  int itmax=20;

#ifdef SIMPSON
  Real t1,t2,L1,L2,L3,L4,dL2,dL3,dL4;
  Real h2,h3,h4;
  Real Lsimp,dLsimp;
#endif

  int count=0, cmax=16;
  int iter;
  int err=0;

  Real temp,dtemp,h1;
  Real epse,epst;

  Real my_dt;

  my_dt=(*dt);

  if( nden != nden) ath_error("[temp_next] Density is NaN\n");
  if( t_old != t_old) ath_error("[temp_next] Temperature is NaN\n");

/* Newton-Rhapson iteration to add cooling terms 
 * using fully implicit discretization. */

  start_iter:

  temp  = t_old;
  for(iter=1;iter<=itmax;iter++) {
    h1=t_old < 1.e4 ? heat: 0.;
#ifdef SIMPSON
/* using Simpson's rule */
    t1=(2.*t_old+temp)/3.;
    t2=(t_old+2.*temp)/3.;
    h2=t1 < 1.e4 ? heat: 0.;
    h3=t2 < 1.e4 ? heat: 0.;
    h4=temp < 1.e4 ? heat: 0.;
    L1=nden*Lam(t_old)-h1;
    L2=nden*Lam(t1)-h2;
    L3=nden*Lam(t2)-h3;
    L4=nden*Lam(temp)-h4;

    dL2=nden*dLam(t1);
    dL3=nden*dLam(t2);
    dL4=nden*dLam(temp);

//    Lsimp=8./(1./L1+3./L2+3./L3+1./L4);
//    dLsimp=0.125*SQR(Lsimp)*(dL2/SQR(L2)+2*dL3/SQR(L3)+dL4/SQR(L4));
    Lsimp=0.125*(L1+3.*L2+3.*L3+L4);
    dLsimp=0.125*(dL2+2*dL3+dL4);

    dtemp=(temp+unitC*Lsimp*my_dt-t_old)/
          (1.0+unitC*dLsimp*my_dt);
#else

/* using fully implicit method */
    dtemp=(temp+unitC*(nden*Lam(temp)-h1)*my_dt-t_old)/
          (1.0+unitC*dLam(temp)*nden*my_dt);
#endif

    temp =(temp - dtemp);
    epse=fabs(dtemp/temp);
    if(epse < toler) goto check_dT; // check convergency of NR iteration
  }

/* If NR didn't converge, recalculate with halved time step */
  if(count <= cmax) {
    if(update)
    ath_perr(-1,"Newton didn't converge %d %d: epse=%g temp=%g dtemp=%g rho=%g t_old=%g dt=%g my_dt=%g\n",
                      update,myID_Comm_world,epse,temp,dtemp,nden,t_old,*dt,my_dt);
    my_dt *= 0.5;
    count++;
    goto start_iter;
  }
      
  check_dT:
/* If NR converged, check relative change of temperature. */
//  epst=fabs(t_old-temp)/t_old;
  epst=fabs(log10(t_old)-log10(temp));
  if (epst > dtempmax && count <= cmax ) { 
    if(update)
    ath_perr(-1,"Temperature change is too large %d: epst=%g temp=%g t_old=%g dt=%g my_dt=%g\n",
                      update,epst,temp,t_old,*dt,my_dt);
    my_dt *= 0.5;
    count++;
    goto start_iter;
  }
/* If NR didn't converge and relative change of temperature is too large
 * after halved dt count times, stop calculation with error message
 */
  if (epst > dtempmax && count > cmax ) 
    ath_error("[cool_solver.c] Not gonna do it %d: told=%e tnew=%e epst=%g count=%d\n",
              update,t_old,temp,epst,count);

  *dt = my_dt;

  return temp;
}

//=======================================================================
// cooling function in unit of erg cm^3 s^{-1}
//=======================================================================
Real Lam(Real temp)
{
  Real c1=1.e7, c2=1.4e-2, t1=1.184e5, t2=92.;
  Real lamtmp;
  Real logt;
  Real beta,C;
 
  logt=log10(temp);

  if(logt < 4.2) {
  lamtmp = c1*exp(-t1/(temp+1.e3))+
           c2*sqrt(temp)*exp(-t2/temp);
  lamtmp *= 2.e-26;
  return (lamtmp);
  } 
#ifdef SD93
  else if(logt < 4.35){
  beta=-1.0;
  C=-17.55;
  }
  else if(logt < 4.90){
  beta=1.63636;
  C=-29.0182;
  }
  else if(logt < 5.40){
  beta=0.0;
  C=-21.0;
  }
  else if(logt < 5.90){
  beta=-1.92;
  C=-10.632;
  }
  else if(logt < 6.25){
  beta=0.0;
  C=-21.96;
  }
  else if(logt < 6.50){
  beta=-1.92;
  C=-9.96;
  }
  else if(logt < 7.50){
  beta=-0.34;
  C=-20.23;
  }
  else {
  beta=0.40;
  C=-25.78;
  }
#endif

#ifdef GF12
  else if(logt < 4.5){
  beta=-1.0;
  C=-17.55;
  }
  else if(logt < 4.90){
  beta=1.6;
  C=-29.25;
  }
  else if(logt < 5.40){
  beta=0.0;
  C=-21.41;
  }
  else if(logt < 6.25){
  beta=-0.647;
  C=-17.9162;
  }
  else if(logt < 6.5){
  beta=-1.92;
  C=-9.96;
  }
  else if(logt < 7.50){
  beta=-0.34;
  C=-20.23;
  }
  else {
  beta=0.40;
  C=-25.78;
  }
#endif

  lamtmp = C+beta*logt;
  lamtmp = pow(10.,lamtmp);
  return (lamtmp);
  
}
//=======================================================================
// derivative of cooling function
//=======================================================================
Real dLam(Real temp)
{
  Real c1=1.e7, c2=1.4e-2, t1=1.184e5, t2=92.;
  Real dlamtmp,beta,C;
  Real logt,lamtmp;

  logt=log10(temp);

  if(logt < 4.2) {
  dlamtmp = c1*t1*exp(-t1/(temp+1.e3))/pow(temp+1.e3,2.)+
            c2*(0.5/sqrt(temp)+t2/pow(temp,1.5))*exp(-t2/temp);
  dlamtmp *= 2.e-26;
  return (dlamtmp);
  } 
#ifdef SD93
  else if(logt < 4.35){
  beta=-1.0;
  C=-17.55;
  }
  else if(logt < 4.90){
  beta=1.63636;
  C=-29.0182;
  }
  else if(logt < 5.40){
  beta=0.0;
  C=-21.0;
  }
  else if(logt < 5.90){
  beta=-1.92;
  C=-10.632;
  }
  else if(logt < 6.25){
  beta=0.0;
  C=-21.96;
  }
  else if(logt < 6.50){
  beta=-1.92;
  C=-9.96;
  }
  else if(logt < 7.50){
  beta=-0.34;
  C=-20.23;
  }
  else {
  beta=0.40;
  C=-25.78;
  }
#endif
#ifdef GF12
  else if(logt < 4.5){
  beta=-1.0;
  C=-17.55;
  }
  else if(logt < 4.90){
  beta=1.6;
  C=-29.25;
  }
  else if(logt < 5.40){
  beta=0.0;
  C=-21.41;
  }
  else if(logt < 6.25){
  beta=-0.647;
  C=-17.9162;
  }
  else if(logt < 6.5){
  beta=-1.92;
  C=-9.96;
  }
  else if(logt < 7.50){
  beta=-0.34;
  C=-20.23;
  }
  else {
  beta=0.40;
  C=-25.78;
  }
#endif 

  lamtmp = C+beta*logt;
  dlamtmp=beta*pow(10.,lamtmp)/temp;
  return (dlamtmp);
}

#endif
