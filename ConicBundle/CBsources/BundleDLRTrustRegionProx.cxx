/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/BundleDLRTrustRegionProx.cxx
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



#include <math.h>
#include <stdlib.h>
#include "BundleDLRTrustRegionProx.hxx"

 
using namespace CH_Matrix_Classes;

namespace ConicBundle {

// *****************************************************************************
//                       BundleDLRTrustRegionProx::provide_inverse_data
// *****************************************************************************

void BundleDLRTrustRegionProx::provide_inverse_data(bool skip_old_fixed_ind) const
{
  if ((!skip_old_fixed_ind)||(old_fixed_ind.dim()==0)){
    if (LinvHt.coldim()!=vecH.rowdim()){
  
      Dinv=D;
      Dinv.inv();
      
      LinvHt.init(vecH,1.,1);
      
      Symmatrix S;
      scaledrankadd(LinvHt,Dinv,S);
      if ((S.rowdim()!=0)&&(norm2(S-LinvHt*Diag(Dinv)*vecH)>=1e-10*(1.+max(diag(S))))){
	if (cb_out()){
	  get_out()<<" precision-problems in BundleDLRTrustRegionProx::provide_inverse_data"<<std::endl;
	  get_out()<<"S=";S.mfile_output(get_out());
	  get_out()<<"LinvHt=";LinvHt.mfile_output(get_out());
	  get_out()<<"D=";D.mfile_output(get_out());
	  get_out()<<"vecH=";vecH.mfile_output(get_out());
	  get_out()<<"LinvHt*Diag(Dinv)*vecH=";(LinvHt*Diag(Dinv)*vecH).mfile_output(get_out());
	  get_out()<<"S-LinvHt*Diag(Dinv)*vecH=";(S-LinvHt*Diag(Dinv)*vecH).mfile_output(get_out());
	  get_out()<<std::endl;
	}
      }
      assert((S.rowdim()==0)||(norm2(S-LinvHt*Diag(Dinv)*vecH)<1e-10*(1.+max(diag(S)))));
      for(Integer i=0;i<S.rowdim();i++)
	S(i,i)+=1.;
      
      if (S.Chol_factor()){
	if (cb_out())
	  get_out()<<"**** WARNING: provide_inverse_data(.): S.Chol_factor() failed";
      }
      LinvHt.scale_cols(Dinv);
      S.Chol_Lsolve(LinvHt);
    
    }
    if (skip_old_fixed_ind){
      indDinv=Dinv;
      indLinvHt=LinvHt;
    }
  }
  else {
    indDinv=D;
    indDinv.delete_rows(old_fixed_ind,true);
    indDinv.inv();
    
    indLinvHt.init(vecH,1.,1);
    indLinvHt.delete_cols(old_fixed_ind,true);
    
    Symmatrix S;
    scaledrankadd(indLinvHt,indDinv,S);
    assert((S.rowdim()==0)||(norm2(S-indLinvHt*Diag(indDinv)*transpose(indLinvHt))<1e-10*max(diag(S))));
    for(Integer i=0;i<S.rowdim();i++)
      S(i,i)+=1.;
      
    if (S.Chol_factor()){
      if (cb_out())
	get_out()<<"**** WARNING: provide_inverse_data(.): S.Chol_factor() failed when skipping fixed indices";
    }
    
    indLinvHt.scale_cols(indDinv);
    S.Chol_Lsolve(indLinvHt);
    
  }
}


// *****************************************************************************
//                       BundleDLRTrustRegionProx::provide_inverse_data
// *****************************************************************************

void BundleDLRTrustRegionProx::clear_inverse_data()
{
  Dinv.init(0,0,0.);
  LinvHt.init(0,0,0.);

  old_fixed_ind.init(0,0,Integer(0));
  indDinv.init(0,0,0.);
  indLinvHt.init(0,0,0.);
}


// *****************************************************************************
//                       BundleDLRTrustRegionProx::clean
// *****************************************************************************

void BundleDLRTrustRegionProx::clean()
{
  if (vecH.coldim()==0)
    return;

  //compute the SVD
  Symmatrix S;
  rankadd(vecH,S,1.,0.,1);
  int retval=S.eig(indLinvHt,Dinv,false);
  if (retval){
    if (cb_out())
      get_out()<<"**** WARNING: BundleDLRTrustRegionProx::clean(...): S.eig failed and returned "<<retval<<" (order="<<S.rowdim()<<")"<<std::endl;
  }
  Real lmax=Dinv(0);
  Integer cnt=1;
  Integer ubcnt=min(max_columns,Dinv.dim());
  while ((cnt<ubcnt)&&(Dinv(cnt)>1e-4*lmax))
    cnt++;

  if (cb_out(2)){
    get_out()<<" BDLRclean("<<cnt;
    if (Dinv.dim()>0)
      get_out()<<","<<lmax<<","<<Dinv(cnt-1);
    get_out()<<")";
  }

  genmult(vecH,indLinvHt,LinvHt);
  swap(vecH,LinvHt);
  
  for(Integer j=cnt;j<Dinv.dim();j++){
    if (Dinv(j)<1e-10*lmax)
      break;
    for(Integer i=0;i<D.dim();i++){
      D(i)+=sqr(vecH(i,j));
    }
  }

  vecH.delete_cols(Range(cnt,vecH.coldim()-1));

  needs_cleaning=false;
}

// *****************************************************************************
//                       BundleDLRTrustRegionProx::init
// *****************************************************************************

void BundleDLRTrustRegionProx::init(const CH_Matrix_Classes::Matrix& in_D, 
	    const CH_Matrix_Classes::Matrix& in_vecH)
{
  //assert column vector
  assert(in_D.coldim()==1);
  //assert matching dimsions
  assert(vecH.rowdim()==in_D.rowdim());
  
  D=in_D;
  D+=weightu;
  vecH=in_vecH; 

  clear_inverse_data();

}


// *****************************************************************************
//                       BundleDLRTrustRegionProx::set_weightu
// *****************************************************************************

void BundleDLRTrustRegionProx::set_weightu(CH_Matrix_Classes::Real in_weightu)
{
  assert(in_weightu>0.);
  if (fabs(weightu-in_weightu)>1e-10*fabs(weightu)){
    D-=weightu;
    weightu=in_weightu;
    D+=weightu;

    clear_inverse_data();
    compute_corr();
  }
}

// *****************************************************************************
//                       BundleDLRTrustRegionProx::norm_sqr
// *****************************************************************************

Real BundleDLRTrustRegionProx::norm_sqr(const Matrix& B) const
{
  Real d=normDsquared(B,D);
  if (vecH.coldim()>0){
    Matrix tmpmat;
    genmult(vecH,B,tmpmat,1.,0.,1);
    d+=ip(tmpmat,tmpmat);
  }
  return d;
}


// *****************************************************************************
//                       BundleDLRTrustRegionProx::dnorm_sqr
// *****************************************************************************

Real BundleDLRTrustRegionProx::dnorm_sqr(const MinorantPointer& B) const
{
  Real d=B.dual_norm_squared(&D);
  if (vecH.coldim()>0){
    provide_inverse_data();
    Real dummy;
    Matrix tmpvec(vecH.rowdim(),1); chk_set_init(tmpvec,1);
    B.get_minorant(dummy,tmpvec,0);
    Matrix tmpmat;
    genmult(LinvHt,tmpvec,tmpmat);
    d-=ip(tmpmat,tmpmat);
  }
  return d;
}

// *****************************************************************************
//                                add_H
// *****************************************************************************
  
int BundleDLRTrustRegionProx::add_H(Symmatrix& big_sym,
					Integer start_index) const
{
  assert((start_index>=0)&&(big_sym.rowdim()-start_index>=D.rowdim()));
  for(Integer i=0;i<D.rowdim();i++)
    big_sym(i+start_index,i+start_index)+=D(i);
  if (vecH.coldim()>0){
    if (big_sym.rowdim()==vecH.rowdim()){
      rankadd(vecH,big_sym,1.,1.);
    }
    else {
      Symmatrix tmpsym;
      rankadd(vecH,tmpsym);
      for(Integer i=0;i<tmpsym.rowdim();i++){
	for(Integer j=i;j<tmpsym.rowdim();j++)
	  big_sym(i+start_index,j+start_index)+=tmpsym(i,j);
      }
    }
  }
  return 0;
}

// *****************************************************************************
//                       BundleDLRTrustRegionProx::add_Hx
// *****************************************************************************

Matrix& BundleDLRTrustRegionProx::add_Hx(const Matrix& x, Matrix& outplusHx,Real alpha) const
{ 
  outplusHx.xpeya(D%x,alpha);
  if (vecH.coldim()>0){
    Matrix tmpmat;
    genmult(vecH,x,tmpmat,alpha,0.,1);
    genmult(vecH,tmpmat,outplusHx,1.,1.);
  }
  return outplusHx;
}

// *****************************************************************************
//                       BundleDLRTrustRegionProx::get_precond
// *****************************************************************************

void BundleDLRTrustRegionProx::get_precond(Matrix& inD,
					   const Matrix*& Vp) const
{ 
  inD.init(D);
  if (vecH.rowdim()>0)
    Vp=&vecH;
  else 
    Vp=0;
}

// *****************************************************************************
//                       BundleDLRTrustRegionProx::compute_QP_costs
// *****************************************************************************


int BundleDLRTrustRegionProx::compute_QP_costs(Symmatrix& Q,
				      Matrix& d,
				      Real& offset,
				      const MinorantPointer& constant_minorant,
				      const MinorantBundle& bundle, 
				      const Matrix& y,
				      const MinorantPointer& groundset_minorant,
				      Indexmatrix* yfixed)
{
  //--- determine the fixed indices and values
  Indexmatrix ind(y.dim(),Integer(1));
  ind.init(0,0,0.);
  Matrix val(y.dim(),Integer(1));
  val.init(0,0,0.);
  _y.newsize(y.dim(),Integer(1)); chk_set_init(_y,1);
  Integer ydim=0;		  
  if (yfixed){
    for(Integer i=0;i<yfixed->dim();i++){
      if ((*yfixed)(i)>0){
	(*yfixed)(i)=1;
	ind.concat_below(i);
	val.concat_below(y(i));
      }
      else {
	_y(ydim++)=y(i);
      }
    }
    _y.reduce_length(ydim);
  }
  else {
    _y=y;
    ydim=y.dim();
  }
  assert(ydim==y.dim()-ind.dim());
  bool fixed_values=(ind.dim()>0);
  bool fixed_changed=!equal(old_fixed_ind,ind);
  if (fixed_changed){
    old_fixed_ind=ind;
  }

  //--- compute the factorization on the current subset of indices
  if (!fixed_values){
    if ((fixed_changed)||(indLinvHt.coldim()!=vecH.rowdim()))
      provide_inverse_data(true);
  }
  else {
    if ((fixed_changed)||(indLinvHt.coldim()+old_fixed_ind.dim()!=vecH.rowdim())){
      provide_inverse_data(true);
    }
  }
    
  //--- get the bundle data into matrix form
  Integer xdim=Integer(bundle.size());
  _A.newsize(ydim,xdim); chk_set_init(_A,1);
  _b.newsize(ydim,Integer(1)); chk_set_init(_b,1); 
  _c.newsize(xdim,Integer(1)); chk_set_init(_c,1);
  _delta=0.;

  if (groundset_minorant.get_minorant(_delta,_b,0,1.,false,fixed_values?&ind:0,fixed_values?&val:0)){
    if (cb_out())
      get_out()<<"*** ERROR in BundleDiagonalProx::compute_QP_costs(...): groundset_minorant.get_minorant failed"<<std::endl;
    return 1;    
  }
  if ((!constant_minorant.empty())&&(constant_minorant.get_minorant(_delta,_b,0,1.,true,fixed_values?&ind:0,fixed_values?&val:0))){
    if (cb_out())
      get_out()<<"*** ERROR in BundleDiagonalProx::compute_QP_costs(...): groundset_minorant.get_minorant failed"<<std::endl;
    return 1;    
  }
      
  for (Integer i=0;i<xdim;i++){
    if (bundle[unsigned(i)].get_minorant(_c(i),_A,i,1,false,fixed_values?&ind:0,fixed_values?&val:0)){
      if (cb_out())
	get_out()<<"*** ERROR in BundleDiagonalProx::compute_QP_costs(...): groundset_minorant.get_minorant failed"<<std::endl;
      return 1;    
    } 
  }

  //compute the linear and part of the constant coefficient
  scaledrankadd(_A,indDinv,Q,1.,0.,1);
  assert(norm2(Q-transpose(_A)*Diag(indDinv)*_A)<1e-10*max(diag(Q)));
  Matrix tmpmat;
  genmult(indLinvHt,_A,tmpmat);
  rankadd(tmpmat,Q,-1.,1.,1);
  old_LinvQ=Q;

  genmult(indLinvHt,_b,oldd);
  offset=_delta+ip(_b,_y)-(normDsquared(_b,indDinv)-ip(oldd,oldd))/2.;

  tmpmat.init(_b,-1.);
  tmpmat%=indDinv;
  tmpmat+=_y;
  genmult(indLinvHt,oldd,tmpmat,1.,1.,1);

  d=_c;
  genmult(_A,tmpmat,d,1.,1.,1);

  oldd=d;
  oldoffset=offset;


  return 0;
}


// *****************************************************************************
//                       BundleDLRTrustRegionProx::update_QP_costs
// *****************************************************************************

int BundleDLRTrustRegionProx::update_QP_costs(Symmatrix& delta_Q, 
							Matrix& delta_d,  
							Real& delta_offset,
							const MinorantPointer& /* constant_minorant */,
							const MinorantBundle& /* bundle */,
							const Matrix& /* y */,
							const MinorantPointer& /* subg */,
							const MinorantPointer& delta_subg,
							const Indexmatrix& delta_index,
							Indexmatrix* yfixed)
{

  //--- determine the changes in fixed indices
  Integer xdim=_A.coldim();
  Integer change_dim=delta_index.dim();
  Indexmatrix new_fixed_ind(change_dim,Integer(1));
  new_fixed_ind.init(0,0,0.);
  Indexmatrix new_fixed_newind(change_dim,Integer(1));
  new_fixed_newind.init(0,0,0.);
  
  _delta+=delta_subg.offset();  
    
  Integer corr_cnt=0;
  for(Integer j=0;j<delta_index.dim();j++){
    Integer ind=delta_index[j];
    if (yfixed){
      switch((*yfixed)(ind)){
      case 0: {
	//was and stayed not fixed, but the subgradient changed
	//treat this later, if no case 2 occured
	Real sgval=delta_subg.coeff(ind);
	while ((corr_cnt<old_fixed_ind.dim())&&(old_fixed_ind(corr_cnt)<ind))
	  corr_cnt++;
	ind-=corr_cnt;
	_b[ind]+=sgval;
	break;
      }
      case 2: {
	new_fixed_ind.concat_below(ind);
	while ((corr_cnt<old_fixed_ind.dim())&&(old_fixed_ind(corr_cnt)<ind))
	  corr_cnt++;
	ind-=corr_cnt;
	new_fixed_newind.concat_below(ind);
	(*yfixed)(ind)=1;
	break;
      }
      default:
	if (cb_out())
	  get_out()<<"*** ERROR in BundleDLRTrustRegionProx::update_QP_costs(...):  internal error, yfixed("<<ind<<")="<<(*yfixed)(ind)<<" should not occur here"<<std::endl;
	std::abort();
      }
    }
    else {
	Real sgval=delta_subg.coeff(ind);
	while ((corr_cnt<old_fixed_ind.dim())&&(old_fixed_ind(corr_cnt)<ind))
	  corr_cnt++;
	ind-=corr_cnt;
	_b[ind]+=sgval;
    }
  }
  
  //--- check whether Hind_chol changed and recompute all if so
  Matrix tmpmat;

  if (new_fixed_ind.dim()>0){
    //prepare indLinvHt
    old_fixed_ind.concat_below(new_fixed_ind);
    Indexmatrix sind;
    sortindex(old_fixed_ind,sind);
    old_fixed_ind=old_fixed_ind(sind);
    provide_inverse_data(true);

    //update the values of _delta and _c and partly Q by the newly fixed values
    _delta+=ip(_b(new_fixed_ind),_y(new_fixed_ind));
    tmpmat=_A.rows(new_fixed_newind);
    genmult(tmpmat,_y(new_fixed_newind),_c,1.,1.,1); 
    rankadd(tmpmat,delta_Q,-1/weightu,0.,1);

    //eliminate the newly fixed rows
    _A.delete_rows(new_fixed_ind);
    _b.delete_rows(new_fixed_ind);
    _y.delete_rows(new_fixed_ind);

    //recompute Q
    delta_Q-=old_LinvQ;
    scaledrankadd(_A,indDinv,old_LinvQ,1.,0.,1);
    assert(norm2(old_LinvQ-transpose(_A)*Diag(indDinv)*_A)<1e-10*max(diag(old_LinvQ)));
    genmult(indLinvHt,_A,tmpmat);
    rankadd(tmpmat,old_LinvQ,-1.,1.,1);
    delta_Q+=old_LinvQ;

    
  }//endif (new_fixed_ind)
  else {
    delta_Q.init(xdim,0.);
  }
  
  delta_offset=-oldoffset;
  delta_d=-oldd;

  genmult(indLinvHt,_b,oldd);
  oldoffset=_delta+ip(_b,_y)-(normDsquared(_b,indDinv)-ip(oldd,oldd))/2.;
  delta_offset+=oldoffset;

  tmpmat.init(_b,-1.);
  tmpmat%=indDinv;
  tmpmat+=_y;
  genmult(indLinvHt,oldd,tmpmat,1.,1.,1);

  oldd=_c;
  genmult(_A,tmpmat,oldd,1.,1.,1);
  delta_d+=oldd;
  
  return 0;
}


// *****************************************************************************
//                       BundleDLRTrustRegionProx::apply_modification
// *****************************************************************************

int BundleDLRTrustRegionProx::apply_modification(const GroundsetModification& gsmdf)
{
  Integer olddim=vecH.rowdim();
  if (gsmdf.old_vardim()!=olddim) {
    if (cb_out())
      get_out()<<"**** ERROR BundleDLRTrustRegionProx::apply_modification: dim="<<olddim<<" but modification assumes "<<gsmdf.old_vardim()<<std::endl;
    return 1;
  }
  
  D.enlarge_below(gsmdf.appended_vardim(),weightu);
  vecH.enlarge_below(gsmdf.appended_vardim(),0.);

  if (gsmdf.map_to_old_variables()){
    vecH=vecH.rows(*(gsmdf.map_to_old_variables()));
  }
  compute_corr();

  clear_inverse_data();

  return 0;
}


// *****************************************************************************
//                       BundleLowRankTrustRegionProx::projected_clone
// *****************************************************************************

      /** @brief in order to allow for fixed variables, this generates a clone restricted to the given indices
   */
BundleProxObject* BundleDLRTrustRegionProx::projected_clone(const CH_Matrix_Classes::Indexmatrix& indices ) 
  {
    Matrix tmpmat=vecH.rows(indices);
    Indexmatrix piv;
    Integer r=tmpmat.QR_factor(piv);
    Matrix tmpvecH=tmpmat.rows(Range(0,r));
    Symmatrix S;
    rankadd(tmpvecH,S,1.,0.,1);
    Real maxval=1.;
    for(Integer i=0;i<S.rowdim();i++)
      maxval=max(maxval,S(i,i));
    S/=maxval;
    Matrix tmplamH;
    if (S.eig(tmpvecH,tmplamH)){
      if (cb_out())
	get_out()<<"**** WARNING BundleLowRankTrustRegionProx::projected_clone(): eig failed"<<std::endl;
    }
    tmpvecH.enlarge_below(tmpmat.rowdim(),0.);
    tmpmat.Q_times(tmpvecH,r);
    tmplamH*=maxval;
    tmplamH.sqrt();
    tmpvecH.scale_cols(tmplamH);
    tmpmat=D(indices);
    tmpmat-=weightu;
 
    BundleDLRTrustRegionProx* pp=new BundleDLRTrustRegionProx(tmpmat,tmpvecH,0,get_use_local_metric(),this,0);
    pp->set_weightu(weightu);
    pp->apply_factor(factor);
    if (get_variable_metric_selection())
      pp->set_variable_metric_selection(get_variable_metric_selection()->clone_VariableMetricSelection());
    return pp;
  }

// *****************************************************************************
//         BundleDLRTrustRegionProx::apply_variable_metric
// *****************************************************************************

int BundleDLRTrustRegionProx::apply_variable_metric(VariableMetricModel* groundset,
						    VariableMetricModel* model,
						    const Matrix& /* aggr */,
						    Integer y_id,
						    const Matrix& y,
						    bool descent_step,
						    Real& current_weight,
						    Real model_maxviol,
						    const Indexmatrix* in_new_indices)
{
  assert(aft_stack.size()==0);

  if(in_new_indices){
    if (cb_out()) 
      get_out()<<"**** WARNING BundleDLRTrustRegionProx::apply_variable_metric(): this implementation does not support updates on subsets of indices -> updating everything"<<std::endl;
  }
    

  new_indices=in_new_indices;

  if (y.dim()==0)
    return 0;
  aft_stack.clear();


  if ((weightu<=0.)||(current_weight<=0.)){
    //not initialized, initialize
    if (current_weight>0.){
      weightu=current_weight;
    }
    else {
      if (weightu<=0.)
	weightu=1.;
      current_weight=weightu;
    }
    vecH.init(0,1,0.);
    D.init(0,1,0.);
  }

  bool add_weightu=false;
  if ((descent_step)||(vecH.rowdim()!=y.dim())){
    add_weightu=true;
    weightu=max(current_weight,1e-10);
    vecH.init(y.dim(),0,0.);
    D.init(y.dim(),1,0.);
    new_indices=0;
  }


  max_columns=min(max(y.dim()/5,30),y.dim());
  //max_columns=min(y.dim(),50);
  needs_cleaning=false;

  int err=0;
  if (groundset->variable_metric_transform()->add_variable_metric(*this,y_id,y,
								  descent_step,
								  weightu,
								  model_maxviol,
								  in_new_indices)){
    if (cb_out()) 
      get_out()<<"**** WARNING BundleDLRTrustRegionProx::apply_variable_metric(): groundset->add_variable_metric(...) failed "<<std::endl;
    err++;
  }
  if (model->variable_metric_transform()->add_variable_metric(*this,y_id,y,
							      descent_step,
							      weightu,
							      model_maxviol,
							      in_new_indices)){
    if (cb_out()) 
      get_out()<<"**** WARNING BundleDLRTrustRegionProx::apply_variable_metric(): model->transform()->add_variable_metric(...) failed "<<std::endl;
    err++;
  }

  if ((needs_cleaning)||(vecH.coldim()>max_columns))
    clean();

  if (cb_out(2)){
    get_out()<<" BDLRTRS(";
    get_out()<<max(D)<<","<<min(D);
    get_out()<<","<<vecH.coldim()<<")";
  }

  current_weight=weightu;
  if (add_weightu) 
    D+=weightu;

  assert(aft_stack.size()==0);
  new_indices=0;
  aft_stack.clear();

  compute_corr();

  clear_inverse_data();
  
  return err;
}


// *****************************************************************************
//         BundleDLRTrustRegionProx::add_variable_metric
// *****************************************************************************

int BundleDLRTrustRegionProx::add_variable_metric(Matrix& in_diagH,
						     Matrix& in_vecH)
{
  if ((in_vecH.coldim()==0)&&(in_diagH.dim()==0))
    return 0;


  //apply the transposed trafos down the stack
  Real fun_factor=1.;
  bool no_trafo=true;
  Matrix tmpmat;
  for(Integer i=Integer(aft_stack.size());--i>=0;){
    fun_factor*=aft_stack[unsigned(i)]->get_fun_coeff();
    if (aft_stack[unsigned(i)]->get_arg_trafo()){
      if (in_vecH.coldim()>0){
	genmult(*aft_stack[unsigned(i)]->get_arg_trafo(),in_vecH,tmpmat,1.,0.,1);
	swap(tmpmat,in_vecH);
	needs_cleaning=true;
      }
      if (in_diagH.dim()>0){
	if (no_trafo){
	  no_trafo=false;
	  tmpmat.init(*aft_stack[unsigned(i)]->get_arg_trafo(),1.,1);
	  in_diagH.sqrt();
	  tmpmat.scale_cols(in_diagH);
	}
	else {
	  genmult(*aft_stack[unsigned(i)]->get_arg_trafo(),in_diagH,tmpmat,1.,0.,1);
	} 
	swap(tmpmat,in_diagH);
	needs_cleaning=true;
      }
    }
  }
  if (fun_factor!=1.){
    in_vecH*=fun_factor;
    in_diagH*=fun_factor;
  }
  
  //add vecH*transpose(vecH) to H
  if ((vecH.coldim()>0)&&(in_vecH.coldim()>0))
    needs_cleaning=true;
  if (in_vecH.coldim()>0)
    vecH.concat_right(in_vecH);
  if (in_diagH.rowdim()>0){
    if (no_trafo)
      D+=in_diagH;
    else 
      vecH.concat_right(in_diagH);
  }
      
  if ((needs_cleaning)&&(vecH.coldim()>2*max_columns))
    clean();
    
  return 0;
}


// *****************************************************************************
//         BundleDLRTrustRegionProx::apply_Hinv
// *****************************************************************************

Matrix& BundleDLRTrustRegionProx::apply_Hinv(Matrix& x) const
{
  //first transform x if necessary
  Real fun_coeff=1.;
  Matrix tmpmat;
  for(Integer i=Integer(aft_stack.size());--i>=0;){
    fun_coeff*=aft_stack[unsigned(i)]->get_fun_coeff();
    if (aft_stack[unsigned(i)]->get_arg_trafo()){
      genmult(*aft_stack[unsigned(i)]->get_arg_trafo(),x,tmpmat,1.,0.,1);
      swap(tmpmat,x);
    }
  }

  //apply the inverse
  if (vecH.coldim()==0)
    x/=D;
  else {
    provide_inverse_data();
    Matrix tmpmat;
    genmult(LinvHt,x,tmpmat);
    x/=D;
    genmult(LinvHt,tmpmat,x,-1.,1.,1);
  }
    
  //apply the trafos up the stack  
  for(unsigned i=0;i<aft_stack.size();i++){
    fun_coeff*=aft_stack[unsigned(i)]->get_fun_coeff();
    if (aft_stack[unsigned(i)]->get_arg_trafo()){
      genmult(*aft_stack[unsigned(i)]->get_arg_trafo(),x,tmpmat);
      swap(tmpmat,x);
    }
  }
 
  if (fun_coeff!=1.)
    x*=fun_coeff;
 
  return x;
}

// *****************************************************************************
//                       BundleDLRTrustRegionProx::mfile_data
// *****************************************************************************


int BundleDLRTrustRegionProx::mfile_data(std::ostream& out) const
{
  out<<"clear weightu qp_diagHscale qp_vecHscale;\n";
  out<<"weightu=";
  out.precision(16);
  out<<weightu<<"\n";
  out<<"qp_diagHscale="; D.mfile_output(out);
  out<<"qp_vecHscale="; vecH.mfile_output(out);
  return 0;
}





}  //namespace ConicBundle
