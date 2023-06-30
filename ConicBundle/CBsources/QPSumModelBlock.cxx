/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/QPSumModelBlock.cxx
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



#include "QPSumModelBlock.hxx"

using namespace CH_Matrix_Classes;

namespace ConicBundle {


  void QPSumModelBlock::clear()
  {
    blocks.clear();
    modeldim=-1;
    constrdim=-1;
    QPModelBlock::clear();
  }
  
  QPSumModelBlock::QPSumModelBlock(CBout* cb,int cbinc):CBout(cb,cbinc),QPModelBlock(cb,cbinc),QPSumModelDataObject(cb,cbinc){clear();}
  
  QPSumModelBlock::~QPSumModelBlock()
  {}

  void QPSumModelBlock::recursive_delete_and_clear()
  {
    for (unsigned i=0;i<blocks.size();i++){
      blocks[i]->recursive_delete_and_clear();
      delete blocks[i];
    }
    clear();
  }

  int QPSumModelBlock::recursive_copy_data_of(QPModelBlockObject* p)
  {
    QPSumModelBlock* sp=dynamic_cast<QPSumModelBlock*>(p);
    if ((sp==0)||(sp->blocks.size()!=blocks.size()))
      return 1;
    
    constant_minorant=sp->constant_minorant;
    bundle=sp->bundle;
    modelx=sp->modelx;
    Bt=sp->Bt;
    modeldx=sp->modeldx;
    modeldcstr=sp->modeldcstr;
    sysviol_model=sp->sysviol_model;
    sysviol_constraints=sp->sysviol_constraints;
    modelx_aggregate=sp->modelx_aggregate;

    modeldim=sp->modeldim;
    constrdim=sp->constrdim;

    int err=0;
    for (unsigned i=0;i<blocks.size();i++){
      err+=blocks[i]->recursive_copy_data_of(sp->blocks[i]);
    }
    return err;
  }

  /// return a cloned object on the heap
  QPModelBlockObject* QPSumModelBlock::clone(){
    QPSumModelBlock* p=new QPSumModelBlock(this,0);
    
    p->constant_minorant=constant_minorant;
    p->bundle=bundle;
    p->modelx=modelx;
    p->Bt=Bt;
    p->modeldx=modeldx;
    p->modeldcstr=modeldcstr;
    p->sysviol_model=sysviol_model;
    p->sysviol_constraints=sysviol_constraints;
    p->modelx_aggregate=modelx_aggregate;

    p->modeldim=modeldim;
    p->constrdim=constrdim;

    p->blocks.resize(blocks.size());
    for (unsigned i=0;i<blocks.size();i++){
      p->blocks[i]=dynamic_cast<QPModelBlock*>(blocks[i]->clone());
      assert(p->blocks[i]);
    }

    return p;
  }


  
  int QPSumModelBlock::append(QPModelDataObject* inblock)
  {
    
    if (inblock==0)
      return 0;
    QPModelBlock* qpbp=dynamic_cast<QPModelBlock*>(inblock);
    if (qpbp==0){
      if (cb_out())
	get_out()<<"**** ERROR in QPSumModelBlock::append(): not a derived QPModelBlock*"<<std::endl;
      return 1;
    }
    blocks.push_back(qpbp);
    modeldim=-1;
    constrdim=-1;
    modelx_changed();

    if (bundle.size()==0){
      //initialize
      bundle.push_back(inblock->get_bundle());
      constant_minorant.push_back(inblock->get_constant_minorant());
    }
    else {
      if (!inblock->get_constant_minorant().empty())
      inblock->get_constant_minorant().get_minorant(get_constant_minorant());
      get_bundle().insert(get_bundle().end(),inblock->get_bundle().begin(),inblock->get_bundle().end());
    }

    return 0;
  }   
  
  //-----------  QPBlock routines
  
   /// returns the dimension of the model set (here the same as the bundle size)
  Integer QPSumModelBlock::dim_model() 
  {
    if (modeldim<0){
      modeldim=0;
      for(unsigned i=0;i<blocks.size();i++)
	modeldim+=blocks[i]->dim_model();
    }
    assert(modeldim<=Integer(get_bundle().size()));
    return modeldim;
  }
  
  /// returns the dimension of the system describing the model set (may contain further constraints)
  Integer QPSumModelBlock::dim_constraints()
  {
    if (constrdim<0){
      constrdim=0;
      for(unsigned i=0;i<blocks.size();i++)
	constrdim+=blocks[i]->dim_constraints();
    }
    return constrdim;
  }
  
  ///returns the dual upper bound to the model value (the trace weighted sum of the dual trace variables); it returns 0. if no model is contained
  Real QPSumModelBlock::constraints_cost(void)
  {
    Real val=0.;
    for (unsigned i=0;i<blocks.size();i++)
      val+=blocks[i]->constraints_cost();
    return val;
  }
  
  
  Matrix& QPSumModelBlock::B_times(const Matrix& A,
			      Matrix& C,
			      Real alpha,
			      Real beta,
			      int Btrans,
			      int Atrans,
			      Integer startindex_model,
			      MinorantBundle& globalbundle,
			      Integer startindex_bundle)
  {
    for(unsigned i=0;i<blocks.size();i++){
      blocks[i]->B_times(A,C,alpha,beta,Btrans,Atrans,startindex_model,
			 globalbundle,startindex_bundle);
      startindex_model+=blocks[i]->dim_model();
      startindex_bundle+=blocks[i]->dim_bundle();
    }
    return C;
  }
    
  Matrix& QPSumModelBlock::times_B(const Matrix& A,
			      Matrix& C,
			      Real alpha,
			      Real beta,
			      int Atrans,
			      int Btrans,
			      Integer startindex_model,
			      MinorantBundle& globalbundle,
			      Integer startindex_bundle)
  {
    for(unsigned i=0;i<blocks.size();i++){
      blocks[i]->times_B(A,C,alpha,beta,Atrans,Btrans,startindex_model,
			 globalbundle,startindex_bundle);
      startindex_model+=blocks[i]->dim_model();
      startindex_bundle+=blocks[i]->dim_bundle();
    }
    return C;
  }

  ///add the main diagonal block tranpose(projection)*diagvec*projection to bigS starting at startindex
  Symmatrix& QPSumModelBlock::add_BDBt(const Matrix& diagvec,
				  Symmatrix& bigS,
				  bool minus,
				  Integer startindex,
				  Matrix& Bt,
				  Integer startindex_model,
				  MinorantBundle& globalbundle,
				  Integer startindex_bundle)
  {
    for(unsigned i=0;i<blocks.size();i++){
      blocks[i]->add_BDBt(diagvec,bigS,minus,startindex,
			  Bt,startindex_model,globalbundle,startindex_bundle);
      startindex_model+=blocks[i]->dim_model();
      startindex_bundle+=blocks[i]->dim_bundle();
    }
    return bigS;
  }

  
  /// get the current matrix for the coupling matrix Bt in the first row of blocks
  Matrix& QPSumModelBlock::get_Bt(Matrix& Bt,
			     Integer startindex_model,
			     MinorantBundle& globalbundle,
			     Integer startindex_bundle)
  {
    for(unsigned i=0;i<blocks.size();i++){
      blocks[i]->get_Bt(Bt,startindex_model,
			 globalbundle,startindex_bundle);
      startindex_model+=blocks[i]->dim_model();
      startindex_bundle+=blocks[i]->dim_bundle();
    }
    return Bt;
  }


  /// set the local modelx value in modelx beginning with startindex (initialize it, do not add)
  int QPSumModelBlock::get_modelx(Matrix& inmodelx, Integer startindex)
  {
    for(unsigned i=0;i<blocks.size();i++){
      blocks[i]->get_modelx(inmodelx,startindex);
      startindex+=blocks[i]->dim_model();
    }
    return 0;
  }
  
  int QPSumModelBlock::get_modeldx(Matrix& inmodeldx, Integer startindex)
  {
    for(unsigned i=0;i<blocks.size();i++){
      blocks[i]->get_modeldx(inmodeldx,startindex);
      startindex+=blocks[i]->dim_model();
    }
    return 0;
  }
  
  int QPSumModelBlock::get_modeldcstr(Matrix& inmodeldcstr, Integer startindex)
  {
    for(unsigned i=0;i<blocks.size();i++){
      blocks[i]->get_modeldcstr(inmodeldcstr,startindex);
      startindex+=blocks[i]->dim_constraints();
    }
    return 0;
  }
  
  /// set the local modelx value in modelx beginning with startindex (initialize it, do not add)
  int QPSumModelBlock::add_modelx_aggregate(Real& val,
				       Matrix& vec,
				       MinorantBundle& global_bundle,
				       Integer startindex_bundle)
  {
    for(unsigned i=0;i<blocks.size();i++){
      blocks[i]->add_modelx_aggregate(val,vec,global_bundle,startindex_bundle);
      startindex_bundle+=blocks[i]->dim_bundle();
    }
    return 0;
  }
  
  int QPSumModelBlock::get_sysviol_model(Matrix& insysviol_model,
				    Integer startindex_model,
				    const Matrix& dy,
				    MinorantBundle& global_bundle,
				    Integer startindex_bundle)
  {
    for(unsigned i=0;i<blocks.size();i++){
      blocks[i]->get_sysviol_model(insysviol_model,startindex_model,dy,global_bundle,startindex_bundle);
      startindex_model+=blocks[i]->dim_model();
      startindex_bundle+=blocks[i]->dim_bundle();
    }
    return 0;
  }
  
  int QPSumModelBlock::get_sysviol_constraints(Matrix& insysviol_constraints, Integer startindex)
  {
    for(unsigned i=0;i<blocks.size();i++){
      blocks[i]->get_sysviol_constraints(insysviol_constraints,startindex);
      startindex+=blocks[i]->dim_constraints();
    }
    return 0;
  }

  void QPSumModelBlock::display_model_values(const Matrix& y,
					MinorantBundle& globalbundle,
					Integer startindex_bundle,
					std::ostream& out)
  {
    for(unsigned i=0;i<blocks.size();i++){
      out<<" block["<<i<<"]:";
      blocks[i]->display_model_values(y,globalbundle,startindex_bundle,out);
      startindex_bundle+=blocks[i]->dim_bundle();
    }
  }
  
  
  // bundlevalues holds the evaluation of the bundle for the current y 
  int QPSumModelBlock::reset_starting_point(const Matrix& y,
					    Real mu,
				       MinorantBundle& globalbundle,
				       Integer startindex_bundle)
  {
    modelx_changed();
    for(unsigned i=0;i<blocks.size();i++){
      blocks[i]->reset_starting_point(y,mu,globalbundle,startindex_bundle);
      startindex_bundle+=blocks[i]->dim_bundle();
    }
    modeldim=-1;
    constrdim=-1;
    return 0;
  }
  

  ///Euclidean norm of constraint violation of modelx
  Real QPSumModelBlock::primalviol_2normsqr()
  {
    Real val=0;
    for(unsigned i=0;i<blocks.size();i++){
      val+=blocks[i]->primalviol_2normsqr();
    }
    return val;
  }
  
  
  ///Euclidean norm of constraint violation on the dual model side
  Real QPSumModelBlock::dualviol_2normsqr()
  {
    Real val=0;
    for(unsigned i=0;i<blocks.size();i++){
      val+=blocks[i]->dualviol_2normsqr();
    }
    return val;
  }
  
  
  int QPSumModelBlock::get_mu_info(Integer& mudim,
				   Real& tr_xz,
				   Real& tr_xdzpdxz,
				   Real& tr_dxdz,
				   Real& min_xz,
				   Real& max_xz) const
  {
    for(unsigned i=0;i<blocks.size();i++){
      blocks[i]->get_mu_info(mudim,tr_xz,tr_xdzpdxz,tr_dxdz,min_xz,max_xz);
    }
    return 0;
  }

  int QPSumModelBlock::get_nbh_info(Integer mudim,
				    Real tr_xz,
				    Real tr_xdzpdxz,
				    Real tr_dxdz,
				    Real nbh_ubnd,
				    Real& alpha,
				    Real& max_nbh,
				    Real& nrmsqr_xz,
				    Real& nrmsqr_xdzpdxz,
				    Real& nrmsqr_dxdz,
				    Real& ip_xz_xdzpdxz,
				    Real& ip_xz_dxdz,
				    Real& ip_dxdz_xdzpdxz) const
  {
    for(unsigned i=0;i<blocks.size();i++){
      blocks[i]->get_nbh_info(mudim,tr_xz,tr_xdzpdxz,tr_dxdz,nbh_ubnd,alpha,max_nbh,
			      nrmsqr_xz,nrmsqr_xdzpdxz,nrmsqr_dxdz,
			      ip_xz_xdzpdxz,ip_xz_dxdz,ip_dxdz_xdzpdxz);
    }
    return 0;
  }

  
  
  /// if necessary, reduce alpha to the biggest value so that feasibility is maintained with this step size
  int QPSumModelBlock::linesearch(Real& alpha) const
  {
    for(unsigned i=0;i<blocks.size();i++){
      blocks[i]->linesearch(alpha);
    }
    return 0;
  }
  
  /// compute the step in the model space given the step in the design space
  int QPSumModelBlock::compute_step(const Matrix& ystep,
			       MinorantBundle& globalbundle,
			       Integer startindex_bundle)
  {
    for(unsigned i=0;i<blocks.size();i++){
      blocks[i]->compute_step(ystep,globalbundle,startindex_bundle);
      startindex_bundle+=blocks[i]->dim_bundle();
    }
    modeldx.init(0,0,0.);
    return 0;
  }
   
   
  /// store this computed step locally and compute the missing local dual step information
  int QPSumModelBlock::computed_step(const Matrix& modelxstep, Integer startindex_model,
				const Matrix& modelconstrstep, Integer startindex_constr)
  {
    for(unsigned i=0;i<blocks.size();i++){
      blocks[i]->computed_step(modelxstep,startindex_model,modelconstrstep,startindex_constr);
      startindex_model+=blocks[i]->dim_model();
      startindex_constr+=blocks[i]->dim_constraints();
    }
    modeldx.init(0,0,0.);
    return 0;
  }
  
  
  /// move in the last computed step direction by a step of length alpha
  int QPSumModelBlock::do_step(Real alpha,
			  const Matrix& y,
			  MinorantBundle& globalbundle,
			  Integer startindex_bundle)
  {
    modelx_changed();
    for(unsigned i=0;i<blocks.size();i++){
      blocks[i]->do_step(alpha,y,globalbundle,startindex_bundle);
      startindex_bundle+=blocks[i]->dim_bundle();
    }
    modeldim=-1;
    constrdim=-1;
    return 0;
  }



  int QPSumModelBlock::add_localrhs(Matrix& globalrhs, 
			       Real rhsmu,
			       Real rhscorr,
			       Integer startindex_model,
			       Integer startindex_constraints,
			       bool append,
			       MinorantBundle& globalbundle,
			       Integer startindex_bundle)
  {
    for(unsigned i=0;i<blocks.size();i++){
      blocks[i]->add_localrhs(globalrhs,rhsmu,rhscorr,startindex_model,startindex_constraints,append,globalbundle,startindex_bundle);
      startindex_model+=blocks[i]->dim_model();
      startindex_bundle+=blocks[i]->dim_bundle();
      startindex_constraints+=blocks[i]->dim_constraints();
    }
    return 0;
  }
  
  
  
  ///add the "scaled" minorant outer products to globalsys, where the correct minorants start at the given index
  int QPSumModelBlock::add_BtinvsysB(Symmatrix& globalsys,
				MinorantBundle& bundle,
				Integer startindex_bundle) 
  {
    for(unsigned i=0;i<blocks.size();i++){
      blocks[i]->add_BtinvsysB(globalsys,bundle,startindex_bundle);
      startindex_bundle+=blocks[i]->dim_bundle();
    }
    return 0;
  }

    


   int QPSumModelBlock::add_localsys(Symmatrix& globalsys,
			   Integer startindex_model,
			   Integer startindex_constraints)
  {
    for(unsigned i=0;i<blocks.size();i++){
      blocks[i]->add_localsys(globalsys,startindex_model,startindex_constraints);
      if (startindex_model>=0)
	startindex_model+=blocks[i]->dim_model();
      if (startindex_constraints>=0)
	startindex_constraints+=blocks[i]->dim_constraints();
    }
    return 0;
  }


  int QPSumModelBlock::localsys_mult(const Matrix& in_vec,
			    Matrix& out_vec,
			    Integer startindex_model,
			    Integer startindex_constraints)
  {
    for(unsigned i=0;i<blocks.size();i++){
      blocks[i]->localsys_mult(in_vec,out_vec,startindex_model,startindex_constraints);
      startindex_model+=blocks[i]->dim_model();
      startindex_constraints+=blocks[i]->dim_constraints();
    }
    return 0;
  }
  

  /** @brief  given the Cholesky factorization LL' of *minus* the blocks A and B (contraints on design variables and Bundle-modelx) and LinvABrhs, solve for the local constraints C and add the new contribution of tracedual*LinvTrace to LinvABsol; store the tracedual in Crhs_and_sol but not yet locally (this will be done by computed_step() ). 
   */
  int QPSumModelBlock::solve_constrsys(const Symmatrix& ABchol,
				  const Matrix& LinvABrhs,
				  Matrix& LinvABsol,
				  Integer startindex_model,
				  Matrix& Crhs_and_sol,
				  Integer startindex_constraints)
  {
    int status=0;
    for(unsigned i=0;i<blocks.size();i++){
      status=blocks[i]->solve_constrsys(ABchol,LinvABrhs,LinvABsol,startindex_model,
				 Crhs_and_sol,startindex_constraints);
      if (status){
	if (cb_out())
	  get_out()<<"***** WARNING QPSumModelBlock::solve_constrsys(......): recursive call failed for block "<<i<<" and returned "<<status<<std::endl;
	break;
      }
      startindex_model+=blocks[i]->dim_model();
      startindex_constraints+=blocks[i]->dim_constraints();
    }
    return status;
  }

  
  /** @brief add the diagonal of the Schur complemented blocks belonging to bundle and local constraints (used for diagonal preconditioning)
  */ 
  int QPSumModelBlock::add_BCSchur_diagonal(Matrix& diagonal,
				   MinorantBundle& globalbundle,
				   Integer startindex_bundle)
  {
    int status=0;
    for(unsigned i=0;i<blocks.size();i++){
      status|=blocks[i]->add_BCSchur_diagonal(diagonal,globalbundle,startindex_bundle);
      startindex_bundle+=blocks[i]->dim_bundle();
    }
    return status;
  }
  

  
  /** @brief append to lowrank "large" columns that should serve well for generating a low rank projection of the Schur complemented model part. For each column i the coordinate sigma_guess(i) gives the Diag_inv-norm for this column. The parameter minval asks to ignore columns whose norms are smaller than minval. If diaginvval is positive, the vector Diag_inv is this value times the all ones vector. 

    On input lowrank must have the correct number of rows already but may
    have 0 columns.  
  */ 
  int QPSumModelBlock::propose_BCSchur_pcsubspace(Matrix& lowrank,
					 Matrix& sigma_guess,
					 const Matrix& Diag_inv,
					 Real minval,
					 Real diaginvval,
					 MinorantBundle& globalbundle,
					 Integer startindex_bundle)
  {
    int status=0;
    for(unsigned i=0;i<blocks.size();i++){
      status|=blocks[i]->propose_BCSchur_pcsubspace(lowrank,sigma_guess,
						    Diag_inv,minval,diaginvval,
						    globalbundle,startindex_bundle);
      startindex_bundle+=blocks[i]->dim_bundle();
    }
    return status;
  }


  /** @brief compute the preconditioning low-rank representation of the Schur complementd blocks belonging to bundle and local constraints by adding a Johnson-Lindenstrauss projection onto the given subspace to glob_lowrank   
  
  */ 
  int QPSumModelBlock::prepare_BCSchur_JLprecond(Matrix& glob_lowrank,
					    Matrix& subspace,
					    bool append_subspace,
					    MinorantBundle& globalbundle,
					    Integer startindex_bundle)
  {
    int status=0;
    for(unsigned i=0;i<blocks.size();i++){
      status|=blocks[i]->prepare_BCSchur_JLprecond(glob_lowrank,subspace,append_subspace,
						   globalbundle,startindex_bundle);
      startindex_bundle+=blocks[i]->dim_bundle();
    }
    return status;
  }
  

  /// add the contributions to glob_rhs of the Schur complemented model block, and return local_rhs of the non complemented constraint block in the rows/columns/diagonal block starting at startindex_constraints
   int QPSumModelBlock::add_Schur_rhs(Matrix& glob_rhs,
				 Matrix* local_rhs,
				 Real rhsmu,
				 Real rhscorr,
				 Integer startindex_constraints,
				 MinorantBundle& globalbundle,
				 Integer startindex_bundle)
  {
    for(unsigned i=0;i<blocks.size();i++){
      blocks[i]->add_Schur_rhs(glob_rhs,local_rhs,rhsmu,rhscorr,
				       startindex_constraints,
				       globalbundle,startindex_bundle);
      startindex_constraints+=blocks[i]->dim_constraints();
      startindex_bundle+=blocks[i]->dim_bundle();
    }
    return 0;
  }

  /// multiply in_vec with the local contribution to the global main block and add it to out_vec; the other local multiplications are carried out externally with the information provide in prepare_Schur_precond and are not done here.
  int QPSumModelBlock::add_Schur_mult(const Matrix& in_vec,
				 Matrix& out_vec,
				 const Matrix* in_cvec,
				 Matrix* out_cvec,
				 Integer startindex_constraints,
				 MinorantBundle& globalbundle,
				 Integer startindex_bundle)
  {
    for(unsigned i=0;i<blocks.size();i++){
      blocks[i]->add_Schur_mult(in_vec,out_vec,in_cvec,out_cvec,startindex_constraints,
				globalbundle,startindex_bundle);
      startindex_constraints+=blocks[i]->dim_constraints();
      startindex_bundle+=blocks[i]->dim_bundle();
    }
    return 0;
  }

  int QPSumModelBlock::computed_Schur_step(const Matrix& xstep,
			  const Matrix& local_step,
			  Integer startindex_constraints,
			  MinorantBundle& globalbundle,
			  Integer startindex_bundle)
  {
    for(unsigned i=0;i<blocks.size();i++){
      blocks[i]->computed_Schur_step(xstep,local_step,
				     startindex_constraints,
				     globalbundle,startindex_bundle);
      startindex_constraints+=blocks[i]->dim_constraints();
      startindex_bundle+=blocks[i]->dim_bundle();
    }
    return 0;
  }

    

}
