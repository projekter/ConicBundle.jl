/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/QPDirectKKTSolver.cxx
    This file is part of ConciBundle, a C/C++ library for convex optimization.

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.

***************************************************************************** */



#include <iomanip>
#include <sstream>
#include <fstream>
#include "QPDirectKKTSolver.hxx"

using namespace CH_Matrix_Classes;

namespace ConicBundle {

  QPDirectKKTSolver::~QPDirectKKTSolver()
  {}
 
// *************************************************************************
//                             clear
// *************************************************************************

void QPDirectKKTSolver::clear()
{  
  dim=0;
  Anr=0;
  bsz=0;
  csz=0;

  Vp=0;

  Hfactor=1.;

  Diag_inv.init(0,0,0.);
  Qchol.init(0,0.);
  Schur_complement.init(0,0.);
  AQiAt_inv.init(0,0.);
  CABinvCt_inv.init(0,0.);
  LinvC.init(0,0,0.);
  piv.init(0,0,Integer(0));
  
  QPKKTSolverObject::clear();
}


// *************************************************************************
//                            QPinit_KKTdata
// *************************************************************************

int QPDirectKKTSolver::QPinit_KKTdata(QPSolverProxObject* in_Hp,
				     QPModelBlockObject* in_model, 
				     const Sparsemat* in_A, 
				     const Indexmatrix* in_eq_indices)
{
  assert(in_Hp);
  clear();
  Hp=in_Hp;
  model=in_model;
  A=in_A;
  eq_indices=in_eq_indices;

  if (A)
    blockA_norm=max(1.,norm2(*A));
  
  dim=-1; //not yet available
  Anr=(A==0)?0:A->rowdim();
  bsz=0;   //bundle size
  csz=0;   //constraint size
  if (model){
    bsz=model->dim_model();
    csz=model->dim_constraints();
  }

  return 0;
}

// *************************************************************************
//                         compute_DLR_Schur_complement
// *************************************************************************

int QPDirectKKTSolver::compute_DLR_Schur_complement(Symmatrix& sc)
{
  sc.init(Anr+bsz+csz,0.);
  
  //for low rank compute a lowrank representation of the inverse
  Matrix VtDi;
  if (Vp){    
    VtDi.init(*Vp,std::sqrt(Hfactor),1); //=V^t
    VtDi.scale_cols(Diag_inv);
  }

  //now inverse(Q+KKTdiagx) = Diag_inv - transpose(VtDi)* Qcholsolve * Qcholsolve * VtDi
      
  //std::cout<<"A="<<A<<std::endl;    
  Symmatrix tmpsym;
  Matrix LiVtDiAt;
  
  if (Anr>0){
    scaledrankadd(*A,Diag_inv,tmpsym);
    //lowrank part
    if (Vp){
      genmult(VtDi,*A,LiVtDiAt,1.,0.,0,1);
      Qchol.Chol_Lsolve(LiVtDiAt);
      rankadd(LiVtDiAt,tmpsym,-1.,1.,1);
    }
    //copy to position
    for(Integer i=0;i<Anr;i++){
      sc(i,i)-=tmpsym(i,i);
      for (Integer j=i+1;j<Anr;j++)
	sc(i,j)-=tmpsym(i,j);
    }
    /// sc now holds   - [ A inverse(Q+diagx) A']  in the A block
  } 
    
	
  if (bsz>0){
    //compute A*Diag_inv*Bundle  and LiVtDiBt (if Vp)
    Matrix tmpB(Anr,bsz); chk_set_init(tmpB,1);
    
    if (Anr>0){
      Matrix tmpA(*A);
      tmpA.scale_cols(Diag_inv);
      model->times_B(tmpA,tmpB,1.,0.,0,1);
      /// tmpB =  A* Diag(Q+diagx)^-1 * B'
    }
    
    if (Vp){
      Matrix LiVtDiBt(Vp->coldim(),bsz); chk_set_init(LiVtDiBt,1);
      model->times_B(VtDi,LiVtDiBt,1.,0.,0,1);
      Qchol.Chol_Lsolve(LiVtDiBt);
      /// LiVtDiBt  now holds  Li * Vt * Diag(Q+diagx)^{-1} * Bt (Inverse grampartB)
      
      //compute transpose(LiVtDiBt)*LiVtDiBt into tmpsym
      rankadd(LiVtDiBt,tmpsym,-1.,0.,1);
      // tmpsym = - grampartB' * grampartB
      
      //copy tmpsym into position
      for (Integer i=0;i<bsz;i++){
	for(Integer j=i;j<bsz;j++){
	  sc(i+Anr,j+Anr)-=tmpsym(i,j);
	}
      }
      /// the B block now holds the lowrank, the diagonal part comes later
      
      //add transpose(LiVtDiAt)*LiVtDiBt
      if (Anr>0)
	genmult(LiVtDiAt,LiVtDiBt,tmpB,-1.,1.,1,0);
      
      /// tmpB =  A * Diag(Q+diagx)^-1 * B' - grampartA' * grampartB
      
    }
    
    //copy tmpB into its position
    for (Integer j=0;j<Anr;j++){
      for (Integer i=0;i<bsz;i++){
	sc(j,Anr+i)-=tmpB(j,i);
	}
    }
    
    //compute and subtract BDiBt, the diagonal part for the B block
    /*
      for (Integer i=0;i<bsz;i++){
      const MinorantPointer& p=model->get_bundle()[unsigned(i)];
      for(Integer j=i;j<bsz;j++){
      AQiAt_inv(i+Anr,j+Anr)-=p.ip(model->get_bundle()[unsigned(j)],0,&Diag_inv);
      }
      }
    */
      model->add_BDBt(Diag_inv,sc,true,Anr);
  } //endif bsz>0
  
  return 0;
}

// *************************************************************************
//                         compute_dense_Schur_complement
// *************************************************************************

int QPDirectKKTSolver::compute_dense_Schur_complement(Symmatrix& sc)
{
  sc.init(Anr+bsz+csz,0.);

  Matrix tmpmat(dim,Anr+bsz+csz);
  tmpmat.init(dim,0,0.);
  if (Anr)
    tmpmat.concat_right(*A,1);
  if (bsz)
    model->get_Bt(tmpmat,Anr);
  Qchol.Chol_Lsolve(tmpmat);
  if (csz)
    tmpmat.enlarge_right(csz,0.);
  rankadd(tmpmat,sc,-1.,0.,1);

  return 0;
}


// *************************************************************************
//                             QPinit_KKTsystem
// *************************************************************************

// set up the KKT System

// KKTdiagx and KTTdiagy are both the positive compelementary systems,
// but the one of KKTdiagy is on the dual side to KKTdiagx and taken with minus 
  
  
int QPDirectKKTSolver::QPinit_KKTsystem(const Matrix& KKTdiagx,
					const Matrix& KKTdiagy,
					Real Hfac,
					Real /* prec */,
					QPSolverParameters* /* params */)
{
  dim=KKTdiagx.rowdim();
  assert((A==0)||(A->rowdim()==KKTdiagy.dim()));
  assert((A==0)||(A->coldim()==dim));
	 
  int status=0;
  
  //------  prepare inverse of KKT Q if necessary
  if ((norm2(KKTdiagx)!=0)
      ||((Hp->is_DLR())&&(Diag_inv.rowdim()==0))
      ||((!Hp->is_DLR())&&(Qchol.rowdim()==0))
      ||(Hfac!=1.)
      ||(Hfactor!=1.)
      ){
    Schur_complement.init(0,0.);
    if (Hp->is_DLR()) {
      //-- low rank case
      Hp->get_precond(Diag_inv,Vp);
      Diag_inv*=Hfac;
      Diag_inv+=KKTdiagx;
      Diag_inv.inv();
      if (Vp){
	scaledrankadd(*Vp,Diag_inv,Qchol,Hfac,0.,1);
	for (Integer i=0;i<Qchol.rowdim();i++){
	  Qchol(i,i)+=1.;
	}
	status=Qchol.Chol_factor(1e-20);
	if (status){
	  if (cb_out()){
	    get_out()<<"**** ERROR: QPDirectKKTSolver::QPinit_KKTsystem(): Chol_factor() failed for low rank inversion and returned "<<status<<std::endl;
	  }
	}
      }
    }
    else {
      //-- dense case
      Qchol.init(dim,0.);
      Hp->add_H(Qchol);
      Qchol*=Hfac;
      for(Integer i=0;i<dim;i++)
	Qchol(i,i)=KKTdiagx(i);
      status=Qchol.Chol_factor(1e-20);
      if (status){
	if (cb_out()){
	  get_out()<<"**** ERROR: QPDirectKKTSolver::QPinit_KKTsystem(): Chol_factor() failed and returned "<<status<<std::endl;
	}
      }
    }
  }

  Hfactor=Hfac;
  
  //--------  form Schur complement 
  if (Anr+bsz+csz>0){
    if ((norm2(KKTdiagx)!=0.)||(Schur_complement.rowdim()!=Anr+bsz+csz)||(Hfactor!=1.)){
      if (Hp->is_DLR()) {
	status=compute_DLR_Schur_complement(Schur_complement);
	if (status){
	  if (cb_out()){
	    get_out()<<"*** ERROR: QPDirectKKTSolver::QPinit_KKTsystem(): compute_DLR_Schur_complement(...): Chol_factor() failed for quadratic term and returned "<<status<<std::endl;
	  }
	}
      }
      else {
	status=compute_dense_Schur_complement(Schur_complement);
	if (status){
	  if (cb_out()){
	    get_out()<<"*** ERROR: QPDirectKKTSolver::QPinit_KKTsystem(): compute_dense_Schur_complement(...): Chol_factor() failed for quadratic term and returned "<<status<<std::endl;
	  }
	}
      }
    }
  }


  //--------- build and factorize the AB or the ABC system matrix
  if ((!status)&&(Schur_complement.rowdim()>0)){

    /// if not explicitly requested otherwise, first try to build the system for AB
    if (!factorize_ABC){
      
      AQiAt_inv=Schur_complement;
	
      // add the diagonal blocks of the complementarity parts
      if (Anr>0){
	for(Integer i=0;i<KKTdiagy.dim();i++){
	  AQiAt_inv(i,i)-=KKTdiagy(i);
	}
      }
      if (model){
	status=model->add_localsys(AQiAt_inv,Anr,Anr+bsz);
	if (status){
	  if (cb_out()){
	    get_out()<<"**** ERROR in QPDirectKKTSolver::QPinit_KKTsystem(...): model->add_localsys(...) failed and returned "<<status<<std::endl;
	  }
	}
      }

      LinvC.newsize(Anr+bsz,csz); chk_set_init(LinvC,1);      
      CABinvCt_inv.newsize(csz); chk_set_init(CABinvCt_inv,1);
      if (csz>0){
	Symmatrix tmpsym;
	swap(AQiAt_inv,tmpsym);
	for(Integer i=0;i<LinvC.rowdim();i++){
	  for(Integer j=0;j<LinvC.coldim();j++){
	    LinvC(i,j)=tmpsym(i,Anr+bsz+j);
	  }
	}
	AQiAt_inv.newsize(Anr+bsz); chk_set_init(AQiAt_inv,1);
	for(Integer i=0;i<AQiAt_inv.rowdim();i++){
	  for(Integer j=i;j<AQiAt_inv.rowdim();j++){
	    AQiAt_inv(i,j)=-tmpsym(i,j);
	  }
	}
	// std::cout<<"AQiAt="<<AQiAt_inv;  //TEST output
	status=AQiAt_inv.Chol_factor();
	if (status){
	  if (cb_out())
	    get_out()<<"**** ERROR QPDirectKKTSolver::QPinit_KKTsystem(...): AQiAt_inv.Chol_factor(.) failed and returned "<<status<<std::endl;
	}
	else {
	  status=AQiAt_inv.Chol_Lsolve(LinvC);
	  if (status){
	    if ( cb_out())
	      get_out()<<"**** ERROR QPDirectKKTSolver::QPinit_KKTsystem(...): AQiAt_inv.Chol_Lsolve(.) failed and returned "<<status<<std::endl;
	    return status;
	  }
	  else {
	    for(Integer i=0;i<CABinvCt_inv.rowdim();i++){
	      for(Integer j=i;j<CABinvCt_inv.rowdim();j++){
		CABinvCt_inv(i,j)=tmpsym(Anr+bsz+i,Anr+bsz+j);
	      }
	    }
	    rankadd(LinvC,CABinvCt_inv,1.,1.,1);
	    // std::cout<<"CABinvCt="<<CABinvCt_inv;   //TEST output
	    status=CABinvCt_inv.Chol_factor();
	    if (status){
	      if (cb_out())
		get_out()<<"**** ERROR QPDirectKKTSolver::QPinit_KKTsystem(...): CABinvCt_inv.Chol_factor(.) failed and returned "<<status<<std::endl;
	    }
	  }
	}
      }
      else {
	AQiAt_inv*=-1.;
	// std::cout<<"AQiAt="<<AQiAt_inv;  //TEST output
	status=AQiAt_inv.Chol_factor();
	if ((status)&& cb_out()){
	  get_out()<<"**** ERROR QPDirectKKTSolver::QPinit_KKTsystem(...): AQiAt_inv.Chol_factor(.) failed and returned "<<status<<std::endl;
	}	    
      }
      //if this failed try ABC
      if (status){
	if (cb_out()){
	  get_out()<<"**** WARNING in QPDirectKKTSolver::QPinit_KKTsystem(...): solving AB with Cholesky failed, switching to factorize_ABC with Aasen "<<std::endl;
	}
	factorize_ABC=true;
	status=0;
      }
    }
    
    
    if ((!status)&&(factorize_ABC)) {
	
      AQiAt_inv=Schur_complement;
      
      // add the diagonal blocks of the complementarity parts
      if (Anr>0){
	for(Integer i=0;i<KKTdiagy.dim();i++){
	  AQiAt_inv(i,i)-=KKTdiagy(i);
	}
      }
      if (model){
	status=model->add_localsys(AQiAt_inv,Anr,Anr+bsz);
	if (status){
	  if (cb_out()){
	    get_out()<<"**** ERROR in LowRank_QPDirectKKTSolver::QPinit_KKTsystem(...): model->add_localsys(...) failed and returned "<<status<<std::endl;
	  }
	}
      }
      	
      if (!status){
	status=AQiAt_inv.Aasen_factor(piv);
	if (status){
	  if (cb_out())
	    get_out()<<"**** ERROR QPDirectKKTSolver::QPinit_KKTsystem(...): AQiAt_inv.Aasen_factor(.) failed and returned "<<status<<std::endl;
	  return status;
	}
      }
    }
  }

  return status;
}

// *************************************************************************
//                             QPsolve_KKTsystem
// *************************************************************************

// solve the KKT System

/// on input
/// dualrhs = -(Qx+A'y+G'modelx+c) +mu (...)
/// primalrhs = -(Ax+s) + mu ()...
int QPDirectKKTSolver::QPsolve_KKTsystem(Matrix& solx,Matrix& soly,
					 const Matrix& primalrhs,
					 const Matrix& dualrhs,
					 Real rhsmu,
					 Real rhscorr,
					 Real /* prec */,
					 QPSolverParameters* /* params */)
{
  assert(dualrhs.dim()==dim);
  assert(primalrhs.dim()==Anr);
  int status=0;

  solx.init(dualrhs);
  soly.init(primalrhs);
  
  ///apply inverse of (Q+KKTdiagx) to solx
  if (Hp->is_DLR()){
    solx%=Diag_inv;
    Matrix tmpmat;
    if (Vp){
      genmult(*Vp,solx,tmpmat,std::sqrt(Hfactor),0.,1);
      status=Qchol.Chol_solve(tmpmat);
      solx.init(dualrhs);
      genmult(*Vp,tmpmat,solx,-std::sqrt(Hfactor),1.);
      solx%=Diag_inv;
    }
  }
  else {
    if (Qchol.rowdim()){
      status=Qchol.Chol_solve(solx);
    }
    else {
      Hp->apply_Hinv(solx);
    }
  }
  //solx is the solution unless Anr+bsz+csz>0

  if (Anr+bsz+csz>0){ 
    //----------- form the right hand side for Schur complement with the Q block
    /// subtract A * (inv solx) from soly  
    if (Anr>0){
      genmult(*A,solx,soly,-1.,1.);
    }
    /// subtract B * (inv solx) from the model rhs side appended to soly (add it here)
    if (model){
      /*
	Integer bsz=model->dim_model();
	for (Integer i=0;i<bsz;i++){
	soly(A.rowdim()+i)-=model->get_bundle()[unsigned(i)].ip(solx);
	}
      */
      Matrix tmpvec(bsz,1); chk_set_init(tmpvec,1);
      model->B_times(solx,tmpvec,-1.,0.);
      soly.concat_below(tmpvec);
      soly.enlarge_below(csz,0.);
      model->add_localrhs(soly,rhsmu,rhscorr,Anr,Anr+bsz,true);
    }
    
    // std::cout<<" rhs="<<soly<<std::endl;   /// TEST
    // Matrix tmprhs=soly;                         ///TEST
    
    
    ///-------------------       solve for the rhs 
    if (factorize_ABC) {
      status=AQiAt_inv.Aasen_solve(soly,piv);
    }
    else {
      if (csz>0){
	Matrix soltr(csz,1,soly.get_store()+Anr+bsz);
	soly.reduce_length(Anr+bsz);
	status=AQiAt_inv.Chol_Lsolve(soly);
	
	genmult(LinvC,soly,soltr,1.,1.,1,0);
	status=CABinvCt_inv.Chol_solve(soltr);
	genmult(LinvC,soltr,soly,-1.,1.);
	status=AQiAt_inv.Chol_Ltsolve(soly);
	
	soly*=-1;
	soly.concat_below(soltr);
      }
      else {
	status=AQiAt_inv.Chol_solve(soly);
	
      soly*=-1;
      }
    }
  
    ///-------------------     extract the solution  
    
    // inform the model block about the solution
    
    solx.init(dualrhs);
    if (model){
      Matrix solmodelx(bsz,1,soly.get_store()+Anr);
      Matrix solmodelconstr(csz,1,soly.get_store()+Anr+bsz);
      soly.reduce_length(Anr);
      model->B_times(solmodelx,solx,-1.,1.,1,0);
      model->computed_step(solmodelx,solmodelconstr);
    }
    
    // recover the Q block solution
    if (Anr>0)
      genmult(*A,soly,solx,-1.,1.,1,0);
    
    if (Hp->is_DLR()){
      solx%=Diag_inv;
      if (Vp){
	Matrix tmpmat;
	genmult(*Vp,solx,tmpmat,std::sqrt(Hfactor),0.,1);
        Qchol.Chol_solve(tmpmat);
	Matrix tmp2;
	genmult(*Vp,tmpmat,tmp2,std::sqrt(Hfactor));
	tmp2%=Diag_inv;
	solx-=tmp2;
      }
    }
    else {
      if (Qchol.rowdim()==dim){
	status=Qchol.Chol_solve(solx);
      }
      else {
	Hp->apply_Hinv(solx);
      }
    }
  }
   
  return status;
}



}

