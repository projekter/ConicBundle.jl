/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/QPModelBlock.cxx
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



#include "QPModelBlock.hxx"
#include "QPSumModelBlock.hxx"
#include "QPConeModelBlock.hxx"

using namespace CH_Matrix_Classes;

namespace ConicBundle {

  QPModelBlock::~QPModelBlock()
  {}

// *****************************************************************************
//                                push_aft
// *****************************************************************************

  int QPModelBlock::push_aft(const AffineFunctionTransformation* aft,
			      const Indexmatrix* global_indices,
			      const Indexmatrix* local_indices,
			      std::map<MinorantPointer,MinorantPointer>* precomputed)
  {
    assert(constant_minorant.size()>0);
    assert(bundle.size()>0);
    constant_minorant.push_back(constant_minorant.back());
    bundle.push_back(bundle.back());
    if (aft==0){
      return 0;
    }
    modelx_aggregate.clear();
    Bt.init(0,0,0.);
    
    MinorantPointer& cm=get_constant_minorant();
    MinorantBundle& bun=get_bundle();

    MinorantPointer tmpm;
    int err=0;
    if (precomputed==0){
      if (!cm.empty()){
	if (aft->transform_minorant(tmpm,cm,1.,true,local_indices,global_indices)){
	  err++;
	  if (cb_out(0)){
	    get_out()<<"\n**** ERROR: QPModelBlock::apply_aft(..): transform_minorant failed for constant_minorant"<<std::endl;
	  }
	}
	cm=tmpm;
      }
      for(unsigned int i=0;i<bun.size();i++){
	tmpm.clear();
	if (aft->transform_minorant(tmpm,bun[i],1.,false,local_indices,global_indices)){
	  err++;
	  if (cb_out(0)){
	    get_out()<<"\n**** ERROR: QPModelBlock::apply_aft(..): transform_minorant failed for minorant "<<i<<std::endl;
	  }
	}
	bun[i]=tmpm;
      }
      return err;
    }

    // precomputed is available, deal with it
    std::map<MinorantPointer,MinorantPointer>::iterator mapit;
    if (!cm.empty()){
      mapit=precomputed->find(cm);
      if ((mapit==precomputed->end())||(! mapit->second.valid())){
	if (aft->transform_minorant(tmpm,cm,1.,false,local_indices,global_indices)){
	  if (cb_out(0)){
	    get_out()<<"\n**** ERROR: QPModelBlock::apply_aft(..): transform_minorant failed for constant_minorant"<<std::endl;
	  }
	  err++;
	}
	(*precomputed)[cm]=tmpm;
	cm=tmpm;
      }
      else {
	cm=mapit->second;
      }
    }
    for(unsigned int i=0;i<bun.size();i++){
      mapit=precomputed->find(bun[i]);
      if ((mapit==precomputed->end())||(!mapit->second.valid())){
	tmpm.clear();
	if (aft->transform_minorant(tmpm,bun[i],1.,false,local_indices,global_indices)){
	  if (cb_out(0)){
	    get_out()<<"\n**** ERROR: QPModelBlock::apply_aft(..): transform_minorant failed for minorant "<<i<<std::endl;
	  }
	  err++;
	}
	(*precomputed)[bun[i]]=tmpm;
	bun[i]=tmpm;
      }
      else {
	bun[i]=mapit->second;
      }
    }

    return err;
    
  }

// *****************************************************************************
//                                pop_aft
// *****************************************************************************

  int QPModelBlock::pop_aft()
  {
    assert(constant_minorant.size()==bundle.size());
    if (bundle.size()<=1)
      return 1;
    bundle.pop_back();
    constant_minorant.pop_back();
    modelx_aggregate.clear();
    Bt.init(0,0,0.);
    return 0;
  }
    

// *****************************************************************************
//                                globalx_cost
// *****************************************************************************
  Real  QPModelBlock::globalx_cost(const Matrix& globalx)
  {
    if (! get_constant_minorant().empty())
      return get_constant_minorant().evaluate(-1,globalx);
    return 0.;
  }

  /// computes and returns C=alpha*B*A+beta*C where B and A may be transposed
  Matrix&  QPModelBlock::B_times(const Matrix& A,
				   Matrix& C,
				   Real alpha,
				   Real beta,
				   int Btrans,
				   int Atrans)
  {
    if (dim_model()==Integer(get_bundle().size()))
      return genmult(get_bundle(),A,C,alpha,beta,(Btrans==0),Atrans);

    return B_times(A,C,alpha,beta,Btrans,Atrans,0,get_bundle(),0);
  }

   /// computes and returns C=alpha*A*(B)+beta*C where A and B may be transposed
   Matrix&  QPModelBlock::times_B(const Matrix& A,
					   Matrix& C,
					   Real alpha,
					   Real beta,
					   int Atrans,
					   int Btrans)
  {
    if (dim_model()==Integer(get_bundle().size()))
      return genmult(A,get_bundle(),C,alpha,beta,Atrans,(Btrans==0));

    return times_B(A,C,alpha,beta,Atrans,Btrans,0,get_bundle(),0);
  }
   

  /// add B*Diag(diagvec)*Bt to S  in the principal block starting at startindex
  Symmatrix& QPModelBlock::add_BDBt(const Matrix& diagvec,
				      Symmatrix& S,
				      bool minus,
				      Integer startindex)
  {
    if (dim_model()==Integer(get_bundle().size())){
      const Integer bsz=Integer(get_bundle().size());
      if (minus) {
	for (Integer i=0;i<bsz;i++){
	  const MinorantPointer& p=get_bundle()[unsigned(i)];
	  for(Integer j=i;j<bsz;j++){
	    S(i+startindex,j+startindex)-=p.ip(get_bundle()[unsigned(j)],0,&diagvec);
	  }
	}
      }
      else {
	for (Integer i=0;i<bsz;i++){
	  const MinorantPointer& p=get_bundle()[unsigned(i)];
	  for(Integer j=i;j<bsz;j++){
	    S(i+startindex,j+startindex)+=p.ip(get_bundle()[unsigned(j)],0,&diagvec);
	  }
	}
      }
      return S;
    }

    if (Bt.coldim()!=dim_model()){
      Bt.newsize(diagvec.rowdim(),dim_model());chk_set_init(Bt,1);
      get_Bt(Bt,0,get_bundle(),0);
    }
    return add_BDBt(diagvec,S,minus,startindex,Bt,0,get_bundle(),0);
  }


  /// get the current matrix for the coupling matrix Bt in the first row
  Matrix&  QPModelBlock::get_Bt(Matrix& inBt,Integer start_col)
  {
    if (Bt.coldim()!=dim_model()){
      Bt.newsize(inBt.rowdim(),dim_model());chk_set_init(Bt,1);
      get_Bt(Bt,0,get_bundle(),0);
    }
    assert(Bt.rowdim()==inBt.rowdim());
    assert(start_col>=0);
    if (inBt.coldim()-start_col<Bt.coldim()){
      inBt.enlarge_right(Bt.coldim()+start_col-inBt.coldim(),0.);
    }
    if (inBt.coldim()==Bt.coldim())
      inBt=Bt;
    else {
      mat_xey(Bt.dim(),inBt.get_store()+start_col*inBt.rowdim(),Bt.get_store());
    }
    return inBt;
  }
  
  /// get the vector formed by all model x variables
  Matrix&  QPModelBlock::get_x()
  {
    if (modelx.dim()!=dim_model()){
      modelx.newsize(dim_model(),1); chk_set_init(modelx,1);
      get_modelx(modelx,0);
    }
    return modelx;
  }
 
  /// get the vector formed by all delta model x variables
  Matrix&  QPModelBlock::get_dx()
  {
    if (modeldx.dim()!=dim_model()){
      modeldx.newsize(dim_model(),1); chk_set_init(modeldx,1);
      if (get_modeldx(modeldx,0)){
	modeldx.init(0,0,0.);
      }
    }
    return modeldx;   
  }

  /// get the vector formed by all delta model x variables
  Matrix&  QPModelBlock::get_dcstr()
  {
    if (modeldcstr.dim()!=dim_constraints()){
      modeldcstr.newsize(dim_constraints(),1); chk_set_init(modeldcstr,1);
      if (get_modeldcstr(modeldcstr,0)){
	modeldcstr.init(0,0,0.);
      }
    }
    return modeldcstr;   
  }

  /// adds opB transposed times modelx and constant affine term to the arguments
  int  QPModelBlock::add_modelx_aggregate(Real& val,
				   Matrix& vec)
  {
    if (! get_constant_minorant().empty())
      get_constant_minorant().get_minorant(val,vec,0,1.,true);
    add_modelx_aggregate(val,vec,get_bundle(),0);
    return 0;
  }

    
    
  /// get the model violation for the current system solution
  Matrix&  QPModelBlock::get_sysviol_model(const Matrix& dy) 
  {
    sysviol_model.newsize(dim_model(),1); chk_set_init(sysviol_model,1);
    get_sysviol_model(sysviol_model,0,dy,get_bundle(),0);
    return sysviol_model;
  }

  
  /// get the constraint violation for the current system solution
  Matrix&  QPModelBlock::get_sysviol_constraints()
  {
    sysviol_constraints.newsize(dim_constraints(),1); chk_set_init(sysviol_constraints,1);
    if (get_sysviol_constraints(sysviol_constraints,0))
      sysviol_constraints.init(0,0,0.);
    return sysviol_constraints;    
  }

  void QPModelBlock::display_model_values(const Matrix& y,std::ostream& out)
  {
    display_model_values(y,get_bundle(),0,out);
  }

  /// initialize the model variables to a strictly feasible "central" starting point; this is the first call when the next QP problem is solved, so other initialization steps may be appropriate as well here. 
  int  QPModelBlock::reset_starting_point(const Matrix& y,Real mu)
  {
    return reset_starting_point(y,mu,get_bundle(),0);
  }
   

  /// compute the step in the model space given the step in the design space
  int  QPModelBlock::compute_step(const Matrix& ystep)
  {
    QPModelBlock::modelx_changed();
    return compute_step(ystep,get_bundle(),0);
  }

  /// store the computed step and compute the missing dual step information
  int  QPModelBlock::computed_step(const Matrix& modelxstep,
				const Matrix& modelconstrstep)
  {return computed_step(modelxstep,0,modelconstrstep,0);}

  /// move in the last computed step direction by a step of length alpha
  int  QPModelBlock::do_step(Real alpha,const Matrix& nexty)
  {
    return do_step(alpha,nexty,get_bundle(),0);
  }

  /// if mu is not zero, always add the centering term for this mu as well, if append is false, add the Schur complement rhs for add_BtinvsysB, if append is true, append the rhs of the local system
  int  QPModelBlock::add_localrhs(Matrix& globalrhs,
		   Real rhsmu,
		   Real rhscorr,
		   Integer startindex_model,
		   Integer startindex_constraints,
		   bool append) 
  {
    return add_localrhs(globalrhs,rhsmu,rhscorr,startindex_model,startindex_constraints,append,get_bundle(),0);
  }

  int QPModelBlock::solve_constrsys(const Symmatrix& ABchol,
			      Matrix& ABCrhs_and_sol,
			      Integer startindex_model,
			      Integer startindex_constraints)
  {
    assert(startindex_constraints=ABCrhs_and_sol.dim()-dim_constraints());
    Matrix Crhs_and_sol(dim_constraints(),1,ABCrhs_and_sol.get_store()+startindex_constraints);
    ABCrhs_and_sol.reduce_length(startindex_constraints);
    int status=ABchol.Chol_Lsolve(ABCrhs_and_sol);
    if (status){
      if (cb_out())
	get_out()<<"**** WARNING QPModelBlock::solve_constrsys(....): ABchol.Chol_Lsolve failed and returned "<<status<<std::endl;
      return status;
    }
    Matrix ABrhs(ABCrhs_and_sol);
    status=solve_constrsys(ABchol,ABrhs,ABCrhs_and_sol,startindex_model,Crhs_and_sol,0);
    if (status){
      if (cb_out())
	get_out()<<"**** WARNING QPModelBlock::solve_constrsys(....): call to inernal solve_constrsys(......) failed and returned "<<status<<std::endl;
      return status;
    }
    status=ABchol.Chol_Ltsolve(ABCrhs_and_sol);
    if (status){
      if (cb_out())
	get_out()<<"**** WARNING QPModelBlock::solve_constrsys(....): ABchol.Chol_Ltsolve failed and returned "<<status<<std::endl;
    }
    ABCrhs_and_sol*=-1;
    ABCrhs_and_sol.concat_below(Crhs_and_sol);
    return status;
  }

 
  /** @brief add the diagonal of the Schur complemented blocks belonging to bundle and local constraints (used for diagonal preconditioning)
  */
  
  int QPModelBlock::add_BCSchur_diagonal(Matrix& diagonal)
  {
    return add_BCSchur_diagonal(diagonal,get_bundle(),0);
  }


  int QPModelBlock::propose_BCSchur_pcsubspace(Matrix& lowrank,
					 Matrix& sigma_guess,
					 const Matrix& Diag_inv,
					 Real minval,
					 Real diaginvval)
  {
    return propose_BCSchur_pcsubspace(lowrank,sigma_guess,Diag_inv,minval,diaginvval,get_bundle(),0);
  }


  /** @brief compute the preconditioning low-rank representation of the Schur complementd blocks belonging to bundle and local constraints by adding a Johnson-Lindenstrauss projection onto the given subspace to glob_lowrank   
  
  */ 
  
  int QPModelBlock::prepare_BCSchur_JLprecond(Matrix& glob_lowrank,
					       Matrix& subspace,
					       bool append_subspace)

  {
    return prepare_BCSchur_JLprecond(glob_lowrank,subspace,append_subspace,get_bundle(),0);
  }
    
  /// add the contributions to glob_diagonal and glob_rhs of the Schur complemented parts, and return local_rhs, local_globblock, local_diagblock of the non complemented parts 
  int QPModelBlock::add_Schur_rhs(Matrix& glob_rhs,
				   Matrix* local_rhs,
				   Real rhsmu,
				   Real rhs_corr)
  {
    if (local_rhs)
      local_rhs->init(dim_constraints(),1,0.);
    return add_Schur_rhs(glob_rhs,local_rhs,rhsmu,rhs_corr,
			 0,get_bundle(),0);
  }

  /// multiply in_vec with the local contribution to the global main block and add it to out_vec; the other local multiplications are carried out externally with the information provide in prepare_Schur_precond and are not done here.
  int QPModelBlock::add_Schur_mult(const Matrix& in_vec,
				    Matrix& out_vec,
				    const Matrix* in_cvec,
				    Matrix* out_cvec)
  {
    return add_Schur_mult(in_vec,out_vec,in_cvec,out_cvec,0,get_bundle(),0);
  }

  /// use the computed step information to also compute the steps of the complemented parts
  int QPModelBlock::computed_Schur_step(const Matrix& xstep,
				  const Matrix& local_step)
  {
    return computed_Schur_step(xstep,local_step,0,get_bundle(),0);
  }
   
  


  int  QPModelBlock::add_BtinvsysB(Symmatrix& globalsys) 
  {
    return add_BtinvsysB(globalsys,get_bundle(),0);
  }

  QPModelPointer::~QPModelPointer()
  {}

  QPSumModelDataObject* QPModelPointer::generate_summodel_data(BundleModel* )
  { return new QPSumModelBlock(this); }
  
  QPConeModelDataObject*  QPModelPointer::generate_conemodel_data(BundleModel*)
  { return new QPConeModelBlock(this); }

 
}
