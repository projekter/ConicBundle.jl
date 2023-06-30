/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/PSCAffineFunction.cxx
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



#include "PSCAffineFunction.hxx"
#include "CMsymsparse.hxx"
#include "lanczpol.hxx"
#include "LanczMaxEig.hxx"
#include <queue>
#include <utility>
#include <string.h>
 

 
using namespace CH_Matrix_Classes;


namespace ConicBundle{


//*****************************************************************************
//                            PSCAffineMinorantExtender
//*****************************************************************************

  PSCAffineMinorantExtender::PSCAffineMinorantExtender(PSCAffineFunction* in_amf):
    CBout()
  {
    assert(in_amf);
    amf=in_amf;
    set_cbout(amf,0);
  }

  PSCAffineMinorantExtender::~PSCAffineMinorantExtender()
  {}

  int PSCAffineMinorantExtender::extend(Minorant& mnrt,int nc,const int* indices) 
  {
    int err=0;
    if (nc>0) {
      const PSCPrimal* psd=dynamic_cast<PSCPrimal*>(mnrt.get_primal());
      if (psd==0){
	err++;
	if (cb_out()){
	  get_out()<<"**** WARNING: PSCAffineMinorantExtender::extend(...): no primal or not a PSCPrimal"<<std::endl; 
	}
      }
      else {
	for(int i=0;i<nc;i++){
	  Integer ind=Integer(indices[i]);
	  Real value;
	  if (psd->primal_ip(value,amf->get_opAt(),ind)){
	    err++;
	    if (cb_out()){
	      get_out()<<"**** WARNING: PSCAffineMinorantExtender::extend(...): psd->primal_ip(...) failed for index="<<ind<<std::endl; 
	    }
	  }
	  else {
	    if (mnrt.add_coeff(ind,value)){
	      err++;
	      if (cb_out()){
		get_out()<<"**** WARNING: PSCAffineMinorantExtender::extend(...): adding a new coefficint to the minorant failed for index="<<ind<<std::endl; 
	      }
	    }
	  }
	}
      }
    }
    return err;
  }

//*****************************************************************************
//                               AMFMaxEigSolver
//*****************************************************************************

  class AMFMaxEigSolver:public CBout
  {
  private:
    
    Integer dim;
    Bigmatrix bigmat;
    bool bigmat_init;
    
    Lanczos* lanczos;
    
    /** if zero, then imprecise eigenvalue computation is allowed;
        for nonegative values the eigenvalue solver has to deliver
        at least this number of eigenvalues and eigenvectors */
    Integer exacteigs;
    
    Integer dense_limit;
        
    
  public:
    void clear()
    { 
      dim=-1;
      bigmat.clear();
      bigmat_init=false;
    }
    
    AMFMaxEigSolver(const CBout* cb=0,int incr=-1):
      CBout(cb,incr)
    {
      dense_limit=50;
      //dense_limit=500;
      lanczos=new Lanczpol;
      //lanczos=new LanczMaxEig;
      assert(lanczos);
      exacteigs=0;
      clear();
    }
    
    ~AMFMaxEigSolver()
    { 
      delete lanczos;
    }
    
    bool is_init(){return bigmat_init;}
    
    void set_exacteigs(Integer i){ exacteigs=max(i,Integer(0));}
    
    void set_dense_limit(Integer lim){dense_limit=max(lim,Integer(0));}

    const Bigmatrix& get_bigmat() const {return bigmat;}
    
    int init(const Matrix& y,const Integer indim, const CoeffmatPointer C,
	     const SparseCoeffmatVector* opAt, bool dense=false)
    {
      clear();
      dim=indim;
      int status=bigmat.init(y,indim,C,opAt,dense);
      if (status==0) {
	bigmat_init=true;
      }
      return status;
    }
    
    
    int evaluate_projection (const Matrix& P, 
			     const double /* relprec */,
			     Matrix& projected_Ritz_vectors,
			     Matrix& projected_Ritz_values)
    {
      assert(bigmat_init);
      //----  compute tmpsym=P'*bigm*P   for P=bundlevecs (same span as bundle)
      bigmat.lanczosmult(P,projected_Ritz_vectors);
      Integer k=P.coldim();
      Symmatrix tmpsym; tmpsym.newsize(k);
      chk_set_init(tmpsym,1);
      tmpsym.xetriu_yza(P,projected_Ritz_vectors);
      
      //---- determine maximum eigenvalue of tmpsym
      int ret_val=tmpsym.eig(projected_Ritz_vectors,projected_Ritz_values,false);
      if (ret_val){
	if (cb_out()){
	  get_out()<<"**** WARNING: AMFMaxEigSolver::evaluate_projection(.....): eig returned "<<ret_val<<std::endl;
	}
	return ret_val;
      }
      if (cb_out(10)){
	get_out().precision(8);
	get_out()<<"  PSCAF eigmod="<<max(projected_Ritz_values)<<std::endl;
      }
      
      return 0;
    }
    
    int compute_projection (const Matrix& P, Symmatrix& S)
    {
      assert(bigmat_init);
      //----  compute tmpsym=P'*bigm*P   for P=bundlevecs (same span as bundle)
      Matrix tmpmat;
      bigmat.lanczosmult(P,tmpmat);
      Integer k=P.coldim();
      S.newsize(k);
      chk_set_init(S,1);
      S.xetriu_yza(P,tmpmat);
      
      return 0;
    }
    
    int evaluate(const Matrix& bundlevecs, 
		 const double relprec, const double Ritz_bound,
		 Matrix& Ritz_vectors, Matrix&  Ritz_values)
    {
      assert(bigmat_init);
      /* for testing:
      (bigmat.get_symrep()).eig(Ritz_vectors,Ritz_values,false);
      std::cout<<"eval: lmax="<<Ritz_values(0)<<std::endl;
      */
      int status=0;
      Integer nreig=0;
      if (bigmat.lanczosdim()<dense_limit){
	(bigmat.get_symrep()).eig(Ritz_vectors,Ritz_values,false);
	nreig=dim;
	if (cb_out(0)){
	  get_out()<<" PSCAF dense: ";
	  get_out()<<" eigmax=";get_out().precision(8);get_out()<<Ritz_values(0)<<std::endl;
	}
      }
      //-------------------- sparse routine for lmax computation
      else {
	lanczos->set_relprec(relprec);
	Integer blocksz=1;
	if (exacteigs==0) {
          /*
	  if (relprec<1e-7){
	    lanczos->enable_stop_above(CB_plus_infinity);
	    nreig=5;
	    blocksz=5;
	  }
	  else{
	  */
	  lanczos->enable_stop_above(Ritz_bound);
	  nreig=1;
	  /*
	  if (bundlevecs.coldim()>5){
	    blocksz=nreig=bundlevecs.coldim()+1;
	    //Integer ncols=min(bigmat.lanczosdim()/3,max(300,3*blocksz));
	    Integer ncols=min(max(bigmat.lanczosdim()/3,2*blocksz),4*blocksz);
	    Integer nblockmult=max(Integer(ncols/blocksz),Integer(2));
	    Integer nchebit=7;
	    Integer multflops=bigmat.lanczosflops();
	    //if (multflops<nblockmult*blocksz*bigmat.lanczosdim()){
	    //  nchebit=10;
	    //}
	    lanczos->set_nblockmult(nblockmult);
	    lanczos->set_nchebit(nchebit);
	  }
	  else {//default
	    blocksz=nreig=1;
	    lanczos->set_nblockmult(-1);
	    lanczos->set_nchebit(-1);
	  }
	  */
    	  /* } */
	}
	else {
	  lanczos->enable_stop_above(CB_plus_infinity);
	  nreig=exacteigs;
	}
	
	bigmat.reset_nmult();
	
	
	Ritz_values.init(0,0,0.);
        if (bundlevecs.dim()!=0) {
	  Ritz_vectors=bundlevecs;
	}
	else {
	  Ritz_vectors.init(0,0,0.);
	}
	status=lanczos->compute(&bigmat,Ritz_values,Ritz_vectors,nreig,blocksz);
       

	/*
	Ritz_values.init(0,0,0.);
	if ((Ritz_vectors.rowdim()!=bigmat.lanczosdim())||(Ritz_vectors.coldim()==0))
	  Ritz_vectors.init(0,0,0.);
	Integer blocksz=1;
	status=lanczos->compute(&bigmat,Ritz_values,Ritz_vectors,nreig,blocksz);
	*/
	  

	nreig=Ritz_values.dim();
	if (cb_out(0)){
	  get_out()<<"  PSCAF Lanczos {"<<lanczos->get_iter()<<",";
	  get_out()<<bigmat.get_nmult()<<","<<nreig<<"} ";
	}
	
	lanczos->get_lanczosvecs(Ritz_values,Ritz_vectors);
	
	if (cb_out(0)){
	  if (nreig>0) {get_out()<<" eigmax=";get_out().precision(8);get_out()<<Ritz_values(0);}
	  else {get_out()<<" Ritz_val=";get_out().precision(8);get_out()<<Ritz_values(0);}
	  get_out()<<" rp=";get_out().precision(2);get_out()<<relprec;
	  get_out()<<" vecs="<<Ritz_values.dim();
	  get_out().precision(8);get_out()<<" bd="<<Ritz_bound<<std::endl;
	}
      }
      return status;
    }
    
    
    void  set_out(std::ostream* o=0,int pril=1)
    {CBout::set_out(o,pril);if (lanczos) lanczos->set_out(o,pril-1);}
    
  };  //end of class AMFMaxEigSolver
  
  
  
//*****************************************************************************
//                            PSCAffineFunction
//*****************************************************************************

  

  int PSCAffineFunction::form_bigmatrix(const Matrix& current_point)
  {
    //initialize the maxeigsolver if this is necessary 
    if (opAt.coldim()!=current_point.dim()){
      if (cb_out())
	get_out()<<"**** ERROR: PSCAffineFunction::form_bigmatrix(.): input dimension="<<current_point.dim()<<" does not match argument dimension="<<opAt.coldim()<<std::endl;
      return 1;
    }

    bool same_y=(last_bigmat_y.coldim()==1);
    if ((same_y)&&(current_point.dim()==last_bigmat_y.dim())){
      for(Integer i=0;i<current_point.dim();++i){      
	if (fabs(current_point(i)-last_bigmat_y(i))>eps_Real*max(fabs(last_bigmat_y(i)),1.)) {
	  same_y=false;
	  break;
	}
      }
    }
    else {
      same_y=false;
    }
    bool must_init=!same_y;
    if (!must_init){
      for (Integer i=0;i<opAt.blockdim().dim();++i){
	if (!maxeigsolver[(unsigned long)(i)]->is_init()){
	  must_init=true;
	  break;
	}
      }
    }
    
    if (must_init){
      for (Integer i=0;i<opAt.blockdim().dim();++i){
	if ((same_y)&&(maxeigsolver[(unsigned long)(i)]->is_init())) continue;
	CoeffmatPointer cp=C(i,0);
	const SparseCoeffmatVector* opAp=opAt.block(i);
	if (maxeigsolver[(unsigned long)(i)]->init(current_point,opAt.blockdim(i),cp,opAp,(C.get_dense_cnt(i)+opAt.get_dense_cnt(i)!=0))){
	  if (cb_out()) {
	    get_out()<<"**** ERROR: PSCAffineFunction::evaluate(...): initialization of Eigenvaluesolver failed for matrix "<<i<<std::endl;
	  }
	  return 1;
	}
      }
      last_bigmat_y=current_point;
    } 

    return 0;
  }

    
  PSCAffineFunction::PSCAffineFunction(const CBout* cb,int incr):
    CBout(cb,incr)
  {
    generating_primal=0;
    check_correctness_flag=true;
    clear();
  }  
  
  PSCAffineFunction::PSCAffineFunction(const SparseCoeffmatMatrix& in_C,
					     const SparseCoeffmatMatrix& in_opAt,
					     PSCPrimal* gen_prim,
					     const CBout* cb,int incr):
    CBout(cb,incr)
  {
    generating_primal=0;
    check_correctness_flag=true;
    clear();
    generating_primal=gen_prim;
    PSCAffineModification amfmod(opAt.coldim(),opAt.blockdim(),this);
    if (amfmod.add_append_blocks(in_C.blockdim(),&in_C,&in_opAt)){
      if (cb_out())
	get_out()<<"**** ERROR: PSCAffineFunction::PSCAffineFunction(.....): modfication with add_append_blocks(..) failed"<<std::endl;
    }
    if (apply_modification(amfmod)){
      if (cb_out())
	get_out()<<"**** ERROR: PSCAffineFunction::PSCAffineFunction(.....): executing the modfication failed"<<std::endl;	
    }    
  }  
  
  PSCAffineFunction::~PSCAffineFunction()
  {
    clear();
  }

  void PSCAffineFunction::clear()
  {
    C.init(Indexmatrix(0,1,Integer(0)),1);
    opAt.clear();
    C.set_cbout(this,0);
    opAt.set_cbout(this,0);
     
    for(unsigned int i=0;i<maxeigsolver.size();++i){
      delete maxeigsolver[i];
    }
    maxeigsolver.clear();

    delete generating_primal;
    generating_primal=0;

    last_bigmat_y.init(0,0,0.); //0 columns means not initialized
    maxvecs=5;
  }
  

  Minorant* PSCAffineFunction::generate_minorant(const Matrix& P)
  {
    if (P.coldim()==0)
      return new Minorant;
    assert(P.rowdim()==sum(opAt.blockdim()));
    Real cval;
    C.Gram_ip(cval,P,0);
    Matrix tmpvec;
    opAt.Gram_ip(tmpvec,P);
    
    PrimalData* p=0;
    if (generating_primal){
      PSCPrimal* sdpp=dynamic_cast<PSCPrimal*>(generating_primal->clone_primal_data());
      if (sdpp!=0) {
	sdpp->assign_Gram_matrix(P);
	p=sdpp;
      }
    }
    return new Minorant(true,cval,tmpvec.rowdim(),tmpvec.get_store(),0,1.,p);
  }

  int  PSCAffineFunction::svec_projection(Matrix& svec_offset,
					     Matrix& svec_coeffs,
					     const Matrix& P,
					     const Indexmatrix* index_subset)
  {
    int err=0;

    assert(P.rowdim()==sum(opAt.blockdim()));
    Symmatrix S;
    Integer s2dim=(P.coldim()*(P.coldim()+1))/2;
    if (index_subset==0){
      svec_coeffs.newsize(opAt.coldim(),s2dim); chk_set_init(svec_coeffs,1);
      for(Integer i=0;i<opAt.coldim();i++){
	err+=opAt.project(S,P,i);
	svec(S,svec_offset);
	mat_xey(s2dim,svec_coeffs.get_store()+i,opAt.coldim(),svec_offset.get_store(),1);
      }
    }
    else {
      Integer dim=index_subset->dim();
      svec_coeffs.newsize(dim,s2dim); chk_set_init(svec_coeffs,1);
      for(Integer i=0;i<dim;i++){
	err+=opAt.project(S,P,(*index_subset)(i));
	svec(S,svec_offset);
	mat_xey(s2dim,svec_coeffs.get_store()+i,dim,svec_offset.get_store(),1);
      }
    }
    err+=C.project(S,P,0);
    svec(S,svec_offset);
    
    return err;
  }
  
  int PSCAffineFunction::evaluate(const Matrix& current_point,
				     const Matrix& bundlevecs, 
				     const double relprec,
				     const double Ritz_bound,
				     Matrix& Ritz_vectors,
				     Matrix&  Ritz_values,
				     PSCPrimalExtender*& primal_extender)
  {
    int err=0;
    primal_extender=0;

    if (form_bigmatrix(current_point)){
      if (cb_out()) {
	get_out()<<"**** ERROR: PSCAffineFunction::evaluate(.......): form_bigmatrix() failed"<<std::endl;
      }
      err++;
      return err;
    }
    
    //evaluate
    int retval=0;
    if (opAt.blockdim().dim()==1){
      retval=maxeigsolver[0]->evaluate(bundlevecs,relprec,Ritz_bound,
				       Ritz_vectors,Ritz_values);
      if (retval){
	if (cb_out()) {
	  get_out()<<"**** ERROR: PSCAffineFunction::evaluate(...): Eigenvaluesolver failed with code "<<retval<<std::endl;
	}
      }
    }
    else {
      typedef std::pair<Real,Integer> Sol;
      std::priority_queue<Sol> priq;
      Integer start_row=0;
      Integer dim=sum(opAt.blockdim());
      Matrix tmp_sol_vecs(dim,maxvecs); chk_set_init(tmp_sol_vecs,1);
      Matrix tmp_sol_vals(maxvecs,1); chk_set_init(tmp_sol_vals,1);
      for (Integer i=0;i<opAt.blockdim().dim();++i){
	Indexmatrix tmp_ind(Range(start_row,start_row+opAt.blockdim(i)-1));
	Matrix tmp_bundle;
        if (bundlevecs.rowdim()>=dim){
	  tmp_bundle=bundlevecs.rows(tmp_ind);
	}
	Matrix tmp_Ritz_vec;
        if (Ritz_vectors.rowdim()>=dim){
	  tmp_Ritz_vec=Ritz_vectors.rows(tmp_ind);
	}
	Matrix tmp_Ritz_val;
	int lretval=maxeigsolver[(unsigned long)(i)]->evaluate(tmp_bundle,relprec,Ritz_bound,
					      tmp_Ritz_vec,tmp_Ritz_val);
	if (lretval){
	  retval|=lretval;
	  if (cb_out()) {
	    get_out()<<"**** ERROR: PSCAffineFunction::evaluate(...): Eigenvaluesolver failed for matrix "<<i<<" with code"<<lretval<<std::endl;
	  }
	}
	else {
	  Indexmatrix ind;
	  sortindex(-tmp_Ritz_val,ind);
	  Integer j=0;
	  while(j<ind.dim()){
	    Integer jj=ind(j);
	    j++;
	    Real val=tmp_Ritz_val(jj);
	    Integer ii;
	    if (Integer(priq.size())<maxvecs){
	      ii=Integer(priq.size());
	    }
	    else if (-priq.top().first<val){
	      ii=priq.top().second;
	      priq.pop();
	    }
	    else {
	      break;
	    }
	    mat_xea(start_row,tmp_sol_vecs.get_store()+ii*dim,0.);
	    mat_xey(opAt.blockdim(i),tmp_sol_vecs.get_store()+ii*dim+start_row,tmp_Ritz_vec.get_store()+jj*opAt.blockdim(i));
	    mat_xea(dim-opAt.blockdim(i)-start_row,tmp_sol_vecs.get_store()+ii*dim+start_row+opAt.blockdim(i),0.);
	    tmp_sol_vals(ii)=val;
	    priq.push(Sol(-val,ii));
	  }  
	}
	start_row+=opAt.blockdim(i);
      } //end for
      tmp_sol_vecs.delete_cols(Range(Integer(priq.size()),maxvecs-1));
      tmp_sol_vals.delete_rows(Range(Integer(priq.size()),maxvecs-1));      
      Indexmatrix sind;
      sortindex(-tmp_sol_vals,sind);
      Ritz_vectors=tmp_sol_vecs.cols(sind);
      Ritz_values=tmp_sol_vals.rows(sind);  
      if (cb_out(1)){
	get_out()<<" PSCAF lmax="; get_out().precision(6);
	for(Integer i=0;i<min(Integer(5),Ritz_values.dim());i++){
	  get_out()<<" "<<Ritz_values(i);
	}
	get_out()<<std::endl;
      }
    }
    
    return err+retval;
  }

  
  int PSCAffineFunction::evaluate_projection(const Matrix& current_point,
					     const Matrix& P, 
					     const double relprec,
					     Matrix& projected_Ritz_vectors,
					     Matrix& projected_Ritz_values)
  {
    int err=0;

    if (form_bigmatrix(current_point)){
      if (cb_out())
	get_out()<<"**** ERROR: PSCAffineFunction::evaluate_projection(.....): form_bigmatrix() failed"<<std::endl;
      err++;
      return err;
    }
    
    //evaluate
    int retval=0;
    if (opAt.blockdim().dim()==1){
      retval=maxeigsolver[0]->evaluate_projection(P,relprec,
				    projected_Ritz_vectors,projected_Ritz_values);
      if (retval){
	if (cb_out()) {
	  get_out()<<"**** ERROR: PSCAffineFunction::evaluate_projection(.....): Eigenvaluesolver failed with code "<<retval<<std::endl;
	}
      }
    }
    else {
      Integer start_row=0;
      Real maxval=CB_minus_infinity;
      Matrix tmp_sol_vecs;
      Matrix tmp_sol_vals;
      for (Integer i=0;i<opAt.blockdim().dim();++i){
	Indexmatrix tmp_ind(Range(start_row,start_row+opAt.blockdim(i)-1));
	Matrix tmp_P=P.rows(tmp_ind);
	Matrix tmp_Ritz_vec=projected_Ritz_vectors;
	Matrix tmp_Ritz_val=projected_Ritz_values;
	int lretval=maxeigsolver[(unsigned long)(i)]->evaluate_projection(tmp_P,relprec,
					   tmp_Ritz_vec,tmp_Ritz_val);
	if (lretval){
	  retval|=lretval;
	  if (cb_out()) {
	    get_out()<<"**** ERROR: PSCAffineFunction::evaluate_projection(.....): Eigenvaluesolver failed for matrix "<<i<<" with code"<<lretval<<std::endl;
	  }
	  continue;
	}
	if (maxval<max(tmp_Ritz_val)){
	  maxval=max(tmp_Ritz_val);
	  tmp_sol_vecs=tmp_Ritz_vec;
	  tmp_sol_vals=tmp_Ritz_val;
	}
      } //end for
      
      projected_Ritz_vectors=tmp_sol_vecs;
      projected_Ritz_values=tmp_sol_vals;      
    }
    
    return err+retval;
  }


  int PSCAffineFunction::left_right_product(int i, const Matrix& E,
			 const Matrix& F, Matrix& G)
{
  int err=0;

  if ((sum(opAt.blockdim())!=F.rowdim())||(sum(opAt.blockdim())!=E.rowdim())){
    if (cb_out())
	get_out()<<"**** ERROR: PSCAffineFunction::left_right_prod(...): dimension do not match"<<std::endl;
    err++;
    return err;
  }
  if (i<0){
    if (cb_out())
      get_out()<<"**** ERROR: PSCAffineFunction::left_right_prod(...): row index i="<<i<<" negative"<<std::endl;
    err++;
    return err;
  }
    
  if (opAt.blockdim().dim()==1){
    CoeffmatPointer cp=opAt(0,i);
    if (cp==0) {
      G.init(E.coldim(),F.coldim(),0.);
      return err;
    }
    cp->left_right_prod(E,F,G);
    return err;
  }
  Integer maxsz=max(opAt.blockdim());
  Matrix tmpE(maxsz,E.coldim());
  Matrix tmpF(maxsz,F.coldim());
  Matrix tmpmat;
  G.init(E.coldim(),F.coldim(),0.);
  Integer firstrow=0;
  for(Integer j=0;j<opAt.blockdim().dim();j++){
    Integer dim=opAt.blockdim(j);
    CoeffmatPointer cp=opAt(j,i);
    if (cp==0) {
      firstrow+=dim;
      continue;
    }
    Indexmatrix rows(Range(firstrow,firstrow+dim));
    tmpE=E.rows(rows);
    tmpF=F.rows(rows);
    cp->left_right_prod(tmpE,tmpF,tmpmat);
    G+=tmpmat;
    firstrow+=dim;
  }
  return err;
} 

 
  int PSCAffineFunction::apply_modification(const PSCAffineModification& amfmod)
  {
    int err=0;

    if (amfmod.no_modification()){
      return err;
    }
    
    if (amfmod.apply_to_PSCAffine(&C,&opAt)){
      err++;
      if (cb_out()){
	get_out()<<"**** ERROR in PSCAffineFunction::apply_modification(.): modification failed"<<std::endl;
      }
    }
  
    if (amfmod.get_reset_primal()){
      delete generating_primal;
      if (amfmod.get_generating_primal()){
	generating_primal=dynamic_cast<PSCPrimal*>(amfmod.get_generating_primal()->clone_primal_data());
	assert(generating_primal);
      }
      else {
	generating_primal=0;
      }
    }

    if (amfmod.variable_modifications()||amfmod.block_modifications()){
      last_bigmat_y.init(0,0,0.);
    }

    if (amfmod.appended_blockdim().dim()>0){
      maxeigsolver.resize((unsigned long)(amfmod.old_blockdim().dim()+amfmod.appended_blockdim().dim()),0);
      for(unsigned int i=(unsigned int)(amfmod.old_blockdim().dim());i<maxeigsolver.size();i++){
	maxeigsolver[i]=new AMFMaxEigSolver(this);
	assert(maxeigsolver[i]);
      }
    }

    if (amfmod.map_to_old_blocks()){
      std::vector<AMFMaxEigSolver*> tmp(unsigned(amfmod.map_to_old_blocks()->dim()));
      for (unsigned int i=0;i<tmp.size();i++){
	unsigned int ind=unsigned((*(amfmod.map_to_old_blocks()))(Integer(i)));
	tmp[i]=maxeigsolver[ind];
	maxeigsolver[ind]=0;
      }
      for (unsigned int i=0;i<maxeigsolver.size();i++){
	delete maxeigsolver[i];
      }
      swap(maxeigsolver,tmp);
    }
    
    return err;
  }
 
  int PSCAffineFunction::apply_modification(const OracleModification& omod,
					       const Matrix* new_center,
					       const Matrix* old_center,
					       bool& discard_obj_in_center,
					       bool& discard_model,
					       bool& discard_aggregates,
					       MinorantExtender*& extender)
  {
    int err=0;

    PSCAffineModification amfmod(0,Indexmatrix(0,0,Integer(0)),this,0);
    const PSCAffineModification* amfmodp=dynamic_cast<const PSCAffineModification*>(&omod);
    if (amfmodp==0){
      amfmod.clear(opAt.coldim(),opAt.blockdim());
      amfmodp=&amfmod;
      if (amfmod.incorporate(omod)){
	err++;
	if (cb_out()){
	  get_out()<<"**** ERROR in PSCAffineFunction::apply_modification(...): inocorporating a general oraclemodifcation failed"<<std::endl;
	}
      }
    }
    
    discard_obj_in_center=false;
    discard_model=false;
    discard_aggregates=false;
    extender=0;

    if (amfmodp->no_modification()){
      return err;
    }

    //--- check whether deletions concern only zero variables/zero matrices
    if (amfmodp->block_modifications()){
      discard_obj_in_center=true;
      discard_model=true;
      discard_aggregates=true;
    }

    if (amfmodp->get_reset_primal()){
      discard_aggregates=true;
    }

    if ((old_center==0)||(new_center==0)) {
      discard_obj_in_center=true;
    }
    else {
      if((!amfmodp->deleted_variables_are_zero(*old_center,opAt))||
	 (!amfmodp->mapped_variables_are_equal(*new_center,*old_center))){
	discard_obj_in_center=true;
      }
    }
    
    //--- apply the modification
    if (apply_modification(*amfmodp)){
      err++;
      if (cb_out()){
	get_out()<<"**** ERROR in PSCAffineFunction::apply_modification(...): modifcation failed for the affine matrix function"<<std::endl;
      }
    }
      
    //--- check whether additions have zero variables and require/allow extensions
    if (!discard_obj_in_center){
      if(!amfmodp->new_variables_are_zero(*new_center,opAt)){
	discard_obj_in_center=true;
      }
    }

    //provide an extender if needed
    if ((!discard_model)&&(!discard_aggregates)&&(!amfmodp->get_skip_extension())){
      extender=new PSCAffineMinorantExtender(this);
    }

    return err;
  }
  
 
  void  PSCAffineFunction::set_out(std::ostream* o,int pril)
  {
    CBout::set_out(o,pril);
    for(unsigned int i=0;i<maxeigsolver.size();i++){
      maxeigsolver[i]->set_cbout(this);
    }
  }

  void  PSCAffineFunction::set_cbout(const CBout* cb,int incr)
  {
    CBout::set_cbout(cb,incr);
    C.set_cbout(cb,incr);
    opAt.set_cbout(cb,incr);
    for(unsigned int i=0;i<maxeigsolver.size();i++){
      maxeigsolver[i]->set_cbout(this);
    }
  }

  std::ostream& PSCAffineFunction::print_problem_data(std::ostream& o) const
  {

    o<<"\nBEGIN_AFFINEMATRIXFUNCTION\n";
    o<<"\nBLOCKS\n";
    for(Integer i=0;i<opAt.blockdim().dim();i++){
      o<<" "<<opAt.blockdim(i);
    }
    o<<"\n";
    o<<"\nVARIABLES\n";
    o<<" "<<opAt.coldim();
    o<<"\n";
    o<<"\nCOST_MATRICES\n";
    for(Integer i=0;i<opAt.blockdim().dim();i++){
      const CoeffmatPointer cm=C(i,0);
      if (cm!=0) {
	o<<"\n"<<i<<"\n";
	cm->out(o);
      }
    }
    o<<"\nCONSTRAINT_MATRICES\n"; 
    for(Integer i=0;i<opAt.blockdim().dim();i++){
      const SparseCoeffmatVector* row=opAt.block(i);
      if (row){
	for(SparseCoeffmatVector::const_iterator it=row->begin();it!=row->end();it++){
	  o<<"\n"<<i<<" "<<it->first<<"\n";
	  it->second->out(o);
	}
      }
    }
    o<<"\nEND_AFFINEMATRIXFUNCTION"<<std::endl;
    return o;
  }

  std::istream& PSCAffineFunction::read_problem_data(std::istream& in)
  {
    clear();
    char next_word[80];
    char next_char;
    if (! in.good()){
      if (cb_out()){
	get_out()<<"*** ERROR in PSCAffineFunction::read_problem_data(): ";
        get_out()<<" instream is not good";
      }
      return in;
    }
    in>>next_word;
    if (strcmp(next_word,"BEGIN_AFFINEMATRIXFUNCTION")!=0){
      if (cb_out()){
	get_out()<<"*** ERROR in PSCAffineFunction::read_problem_data(): ";
        get_out()<<"expected BEGIN_AFFINEMATRIXFUNCTION but got "<<next_word<<std::endl;
      }
      in.clear(in.rdstate()|std::ios::failbit);
      return in;
    }
    in>>next_word;
    if (strcmp(next_word,"BLOCKS")!=0){
      if (cb_out()){
	get_out()<<"*** ERROR in PSCAffineFunction::read_problem_data(): ";
        get_out()<<"expected BLOCKS but got "<<next_word<<std::endl;
      }
      in.clear(in.rdstate()|std::ios::failbit);
      return in;
    }
    Indexmatrix blockdim(0,1,Integer(0));
    do {
      in>>std::ws;
      next_char=char(in.peek());
      if ((next_char<'0')||(next_char>'9'))
	break;
      Integer next_dim;
      in>>next_dim;
      blockdim.concat_below(next_dim);
    } while (in.good());
    in>>next_word;
    if (strcmp(next_word,"VARIABLES")!=0){
      if (cb_out()){
	get_out()<<"*** ERROR in PSCAffineFunction::read_problem_data(): ";
        get_out()<<"expected VARIABLES but got "<<next_word<<std::endl;
      }
      in.clear(in.rdstate()|std::ios::failbit);
      return in;
    }
    in>>std::ws;
    next_char=char(in.peek());
    if ((next_char<'0')||(next_char>'9')){
      if (cb_out()){
	get_out()<<"*** ERROR in PSCAffineFunction::read_problem_data(): ";
        get_out()<<"expected nonnegative integer variable but got next character "<<next_char<<std::endl;
      }
      in.clear(in.rdstate()|std::ios::failbit);
      return in;
    }
    Integer ydim;
    in>>ydim;
    in>>next_word;
    if (strcmp(next_word,"COST_MATRICES")!=0){
      if (cb_out()){
	get_out()<<"*** ERROR in PSCAffineFunction::read_problem_data(): ";
        get_out()<<"expected COST_MATRICES but got "<<next_word<<std::endl;
      }
      in.clear(in.rdstate()|std::ios::failbit);
      return in;
    }
    Indexmatrix indi(0,1,Integer(0));
    Indexmatrix indj(0,1,Integer(0));
    CoeffmatVector cmvec;
    do {
      in>>std::ws;
      next_char=char(in.peek());
      if ((next_char<'0')||(next_char>'9'))
	break;
      Integer next_i;
      in>>next_i;
      CoeffmatPointer cmp=coeffmat_read(in);
      if (cmp==0){
	if (cb_out()){
	  get_out()<<"*** ERROR in PSCAffineFunction::read_problem_data(): ";
	  get_out()<<"coeffmat_read failed in reading the cost coefficient matrix for "<<next_i<<std::endl;
	}
	in.clear(in.rdstate()|std::ios::failbit);
	return in;
      }
      indi.concat_below(next_i);
      indj.concat_below(0);
      cmvec.push_back(cmp);
    } while (in.good());
    if (C.init(blockdim,1,&indi,&indj,&cmvec)){
      if (cb_out()){
	get_out()<<"*** ERROR in PSCAffineFunction::read_problem_data(): ";
	get_out()<<"C.init() failed"<<std::endl;
      }
      in.clear(in.rdstate()|std::ios::failbit);
      return in;
    }
    in>>next_word;
    if (strcmp(next_word,"CONSTRAINT_MATRICES")!=0){
      if (cb_out()){
	get_out()<<"*** ERROR in PSCAffineFunction::read_problem_data(): ";
        get_out()<<"expected CONSTRAINT_MATRICES but got "<<next_word<<std::endl;
      }
      in.clear(in.rdstate()|std::ios::failbit);
      return in;
    }
    indi.init(0,1,Integer(0));
    indj.init(0,1,Integer(0));
    cmvec.clear();
    do {
      in>>std::ws;
      next_char=char(in.peek());
      if ((next_char<'0')||(next_char>'9'))
	break;
      Integer next_i;
      in>>next_i;
      Integer next_j;
      in>>next_j;
      if (next_j<0) {
	if (cb_out()){
	  get_out()<<"*** ERROR in PSCAffineFunction::read_problem_data(): ";
	  get_out()<<"coeffmat_read failed in reading the column index to row "<<next_i<<": "<<next_j<<std::endl;
	}
	in.clear(in.rdstate()|std::ios::failbit);
	return in;
      }
      CoeffmatPointer cmp=coeffmat_read(in);
      if (cmp==0){
	if (cb_out()){
	  get_out()<<"*** ERROR in PSCAffineFunction::read_problem_data(): ";
	  get_out()<<"coeffmat_read failed in reading the constraint coefficient matrix for ("<<next_i<<","<<next_j<<")"<<std::endl;
	}
	in.clear(in.rdstate()|std::ios::failbit);
	return in;
      }
      indi.concat_below(next_i);
      indj.concat_below(next_j);
      cmvec.push_back(cmp);
    } while (in.good());
    if (opAt.init(blockdim,1,&indi,&indj,&cmvec)){
      if (cb_out()){
	get_out()<<"*** ERROR in PSCAffineFunction::read_problem_data(): ";
	get_out()<<"opAt.init failed"<<std::endl;
      }
      in.clear(in.rdstate()|std::ios::failbit);
      return in;
    }
    in>>next_word;
    if (strcmp(next_word,"END_AFFINEMATRIXFUNCTION")!=0){
      if (cb_out()){
	get_out()<<"*** ERROR in PSCAffineFunction::read_problem_data(): ";
        get_out()<<"expected END_AFFINEMATRIXFUNCTION but got "<<next_word<<std::endl;
      }
      in.clear(in.rdstate()|std::ios::failbit);
      return in;
    }
    return in;
  }


  std::ostream& PSCAffineFunction::print_problem_data_to_mfile(std::ostream& o,Integer block_nr) const
  {

    o<<"\n% BEGIN_AFFINEMATRIXFUNCTION "<<block_nr<<"\n";
    o<<"\n% BLOCKS within this block\n";
    o<<"blk{"<<block_nr<<",1} = 's';\n";
    o<<"blk{"<<block_nr<<",2} = [";
    for(Integer i=0;i<opAt.blockdim().dim();i++){
      o<<" "<<opAt.blockdim(i);
    }
    o<<"];\n";
    Integer dimsum=sum(opAt.blockdim());
    Integer dim2sum=(dimsum*(dimsum+1))/2;
    o<<"\n% VARIABLES, number of\n";
    o<<" ydim="<<opAt.coldim();
    o<<"\n% COST_MATRIX\n";
    Indexmatrix indi;
    Indexmatrix indj;
    Matrix val;
    Symmatrix tmpsym;
    Sparsesym tmpssym;
    Sparsemat tmpsmat;
    for(Integer i=0;i<C.blockdim().dim();i++){
      const CoeffmatPointer cm=C(i,0);
      if (cm==0)
	continue;
      Integer joffset=1;
      for(Integer j=0;j<i;j++){
	joffset+=opAt.blockdim(i);
      }
      if (opAt.blockdim(i)>5000){
	o<<"error('matrix dimension exceeds 5000, aborting to avoid allocation problems');"<<std::endl;
	o.fail();
	return o;
      } 
      Indexmatrix ii;
      Indexmatrix ij;
      Matrix v;
      cm->make_symmatrix(tmpsym);
      tmpsmat.init(tmpsym);
      tmpsmat.get_edge_rep(ii,ij,v);
      ii+=joffset;
      ij+=joffset;
      indi.concat_below(ii);
      indj.concat_below(ij);
      val.concat_below(v);
    }
    o<<"indi=[";
    for(Integer j=0;j<indi.dim();j++){
      o<<indi(j);
      if (j<indi.dim()-1)
	o<<"\n";
    }
    o<<"];\n";
    o<<"indj=[";
    for(Integer j=0;j<indj.dim();j++){
      o<<indj(j);
      if (j<indj.dim()-1)
	o<<"\n";
    }
    o<<"];\n";
    o<<"val=[";
    for(Integer j=0;j<val.dim();j++){
      o<<std::setprecision(16)<<val(j);
      if (j<val.dim()-1)
	o<<"\n";
    }
    o<<"];\n";
    o<<"C{"<<block_nr<<"}=sparse(indi,indj,val,"<<dimsum<<","<<dimsum<<");\n";

    o<<"\n% CONSTRAINT_MATRICES\n"; 
    indi.init(0,0,0.);
    indj.init(0,0,0.);
    val.init(0,0,0.);
    for(Integer colj=0;colj<opAt.coldim();colj++){
      const SparseCoeffmatVector* col=opAt.column(colj);
      if (col==0)
	continue;
      for(SparseCoeffmatVector::const_iterator rowit=col->begin();rowit!=col->end();rowit++){
	Integer j=rowit->first;
	Integer joffset=1;
	for(Integer i=0;i<j;i++){
	  joffset+=opAt.blockdim(i);
	}
	Indexmatrix ii;
	Indexmatrix ij;
	Matrix v;
	rowit->second->make_symmatrix(tmpsym);
	tmpssym.init(tmpsym);
	tmpssym.get_edge_rep(ii,ij,v);
	ii+=joffset;
	ij+=joffset;
 	for(Integer k=0;k<ii.dim();k++){
	  indi.concat_below((ij(k)*(ij(k)-1))/2+ii(k));
	  indj.concat_below(colj+1);
	  if (ij(k)==ii(k))
	    val.concat_below(v(k));
	  else
	    val.concat_below(sqrt(2)*v(k));
	}
      }
    }
    o<<"indi=[";
    for(Integer i=0;i<indi.dim();i++){
      o<<indi(i);
      if (i<indi.dim()-1)
	o<<"\n";
    }
    o<<"];\n";
    o<<"indj=[";
    for(Integer i=0;i<indj.dim();i++){
      o<<indj(i);
      if (i<indj.dim()-1)
	o<<"\n";
    }
    o<<"];\n";
    o<<"val=[";
    for(Integer i=0;i<val.dim();i++){
      o<<std::setprecision(16)<<val(i);
      if (i<val.dim()-1)
	o<<"\n";
    }
    o<<"];\n";
    o<<"At{"<<block_nr<<"}=sparse(indi,indj,val,"<<dim2sum<<","<<opAt.coldim()<<");\n";
    
    o<<"clear indi indj val;\n";
    o<<"\n% END_AFFINEMATRIXFUNCTION"<<std::endl;
    return o;
  }




}

