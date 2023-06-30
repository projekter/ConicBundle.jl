/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/AffineFunctionTransformation.cxx
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


#include <set>
#include "AffineFunctionTransformation.hxx"

 
using namespace CH_Matrix_Classes;

namespace ConicBundle {

// *****************************************************************************
//                           ~AffineFunctionTransformation()
// *****************************************************************************


  AffineFunctionTransformation::~AffineFunctionTransformation()
  {
    delete linear_cost;
    delete arg_offset;
    delete arg_trafo;
  }

// *****************************************************************************
//                        init
// *****************************************************************************


int AffineFunctionTransformation::init(
				       Real in_fun_coeff,
				       Real in_fun_offset,
				       Matrix* in_linear_cost,
				       Matrix* in_arg_offset,
				       Sparsemat* in_arg_trafo,
				       bool in_model_calls_delete)
{
  int err=0;
  if ((in_linear_cost)&&(in_linear_cost->coldim()!=1)){
    err=1;
    if (cb_out()){
      get_out()<<"**** ERROR AffineFunctionTransformation::init: linear_cost must be a columnt vector but instead of 1 it has column dimension ="<<in_linear_cost->coldim()<<std::endl;
    }
  }
  if ((in_arg_offset)&&(in_arg_offset->coldim()!=1)){
    err=1;
    if (cb_out()){
      get_out()<<"**** ERROR AffineFunctionTransformation::init: arg_offset must be a columnt vector but instead of 1 it has column dimension ="<<in_arg_offset->coldim()<<std::endl;
    }
  }
  if ((in_arg_trafo==0)&&(in_linear_cost)&&(in_arg_offset)&&
      (in_linear_cost->rowdim()!=in_arg_offset->rowdim())){
    err=1;
    if (cb_out()){
      get_out()<<"**** ERROR AffineFunctionTransformation::init: identity transformation but rowdim(argument_offset)="<<in_arg_offset->rowdim()<<" differs from rowdim(linear_cost)="<<in_linear_cost->rowdim()<<std::endl;
    }
  }
  if ((in_arg_trafo)&&(in_linear_cost)&&
      (in_linear_cost->rowdim()!=in_arg_trafo->coldim())){
    err=1;
    if (cb_out()){
      get_out()<<"**** ERROR AffineFunctionTransformation::init: coldim(arg_trafo)="<<in_arg_trafo->coldim()<<" differs from rowdim(linear_cost)"<<in_linear_cost->rowdim()<<std::endl;
    }
  }
  if ((in_arg_trafo)&&(in_arg_offset)&&
      (in_arg_trafo->rowdim()!=in_arg_offset->rowdim())){
    err=1;
    if (cb_out()){
      get_out()<<"**** ERROR AffineFunctionTransformation::init: rowdim(arg_trafo)="<<in_arg_trafo->rowdim()<<" differs from rowdim(arg_offset)"<<arg_offset->rowdim()<<std::endl;
    }
  }
  if (err)
    return err;
 
  delete linear_cost; 
  delete arg_offset; 
  delete arg_trafo; 
  
  fun_coeff=in_fun_coeff;
  fun_offset=in_fun_offset;
  linear_cost=in_linear_cost;
  arg_offset=in_arg_offset;
  arg_trafo=in_arg_trafo;
  model_calls_delete=in_model_calls_delete;

  Minorant* mnrt=new Minorant(true,fun_offset);
  if (linear_cost){
    mnrt->add_coeffs(linear_cost->dim(),linear_cost->get_store());
    //mnrt->sparsify(1e-12*norm2(*linear_cost));
  }
  constant_minorant.init(mnrt,0);
  
  return 0;
}
  
// *****************************************************************************
//                              scaled_index_subset
// *****************************************************************************

  bool AffineFunctionTransformation::scaled_index_subset(const Indexmatrix* col_ind,
							 const Indexmatrix* row_ind) const

    {
      if (arg_offset!=0) return false;
      if (arg_trafo==0) return true;
      if (col_ind==0){
	if (max((arg_trafo->get_colinfo()).col(1))>1)
	  return false;
      }
      else {
	Integer nind=col_ind->dim();
	Integer i=0;
	while ((i<nind)&&(arg_trafo->col_nonzeros((*col_ind)(i))<=1))
	  i++;
	if (i<nind)
	  return false;
      }
      if (row_ind==0){
	if (max((arg_trafo->get_rowinfo()).col(1))>1)
	  return false;
      }
      else {
	Integer nind=row_ind->dim();
	Integer i=0;
	while ((i<nind)&&(arg_trafo->row_nonzeros((*row_ind)(i))<=1))
	  i++;
	if (i<nind)
	  return false;
      }
      return true;
    }

// *****************************************************************************
//                              scaled_index
// *****************************************************************************

  bool AffineFunctionTransformation::scaled_index(CH_Matrix_Classes::Integer& mapped_index,
						  CH_Matrix_Classes::Real& coeff,
						  CH_Matrix_Classes::Integer index) const
  {
    if (arg_offset!=0) return false;
    if (arg_trafo==0) {mapped_index=index; coeff=1.; return true;}
    Integer i;
    Integer nz=arg_trafo->col_nonzeros(index,&i);
    if (nz==0) {mapped_index=-1; return true;}
    if (nz==1) {
      mapped_index=arg_trafo->get_colindex()(i);
      coeff=arg_trafo->get_colval()(i);
      return true;
    }
    return false;
  }



// *****************************************************************************
//                           copy_traforows
// *****************************************************************************

int 
AffineFunctionTransformation::copy_traforows(CH_Matrix_Classes::Matrix& copy_to,
					     const CH_Matrix_Classes::Matrix& copy_from) const
{
  if (arg_trafo==0){
    copy_to=copy_from;
    return 0;
  }
  assert((copy_to.dim()==arg_trafo->rowdim())&&(copy_from.dim()==arg_trafo->rowdim()));
  const Indexmatrix& rowinfo=arg_trafo->get_rowinfo();
  for (Integer i=0;i<rowinfo.rowdim();i++){
    Integer ind=rowinfo(i,0);
    copy_to(ind)=copy_from(ind);
  }
  return 0;
}



// *****************************************************************************
//                           transform_argument
// *****************************************************************************

const Matrix&
AffineFunctionTransformation::transform_argument(Matrix& out_y, 
						 Real& offset,
						 const Matrix& in_y) const
{
  offset=fun_offset;
  if (linear_cost)
    offset+=ip(*linear_cost,in_y);

  if (arg_trafo==0){ //identity matrix
    if (arg_offset==0) {
      //no transformation at all
      out_y.init(0,0,0.);
      return in_y;
    }
    
    out_y.init(*arg_offset);
    out_y+=in_y;
    return out_y;
  }

  if ((out_y.rowdim()==arg_trafo->rowdim()) &&
      (out_y.coldim()==1)&&
      (sparse_argument_changes())){
    Integer nrows=arg_trafo->get_rowinfo().rowdim();
    const Integer *rowip=arg_trafo->get_rowinfo().get_store();
    Real* out_yp=out_y.get_store();
    if (arg_offset) {
      const Real* argop=arg_offset->get_store();
      for (Integer i=nrows;--i>=0;rowip++){
	*(out_yp+*rowip)=*(argop+*rowip);
      }
    }
    else {
      for (Integer i=nrows;--i>=0;){
        *(out_yp+*rowip++)=0.;
      }
    }
    genmult(*arg_trafo,in_y,out_y,1.,1.);
  }
  else {
    if (arg_offset){
      out_y.init(*arg_offset);
      genmult(*arg_trafo,in_y,out_y,1.,1.);
    }
    else
      genmult(*arg_trafo,in_y,out_y,1.,0.);
  }
  
  return out_y;
}



// *****************************************************************************
//                           modified_transform_argument
// *****************************************************************************

const Matrix&
AffineFunctionTransformation::modified_transform_argument(Matrix& out_y,
							  const Matrix& in_y,
							  const AFTModification* aftmdf,
							  const GroundsetModification& gsmdf) const
{
  /*
  AffineFunctionTransformation* test_aft=0;
  assert((test_aft=new AffineFunctionTransformation(1.,0.,0,arg_offset?new Matrix(*arg_offset):0,arg_trafo?new Sparsemat(* arg_trafo):0))!=0);
  assert(test_aft->apply_modification(aftmdf,gsmdf)==0);
  Matrix testoy,testdummy;
  Real testoffset=0.;
  assert(&(testoy=test_aft->transform_argument(testdummy,testoffset,in_y)));
  delete test_aft;
  */

  if ((aftmdf==0)&&(gsmdf.no_modification())){
    Real dummy=0.;
    const Matrix& mat=transform_argument(out_y,dummy,in_y);
    //assert(norm2(out_y-testoy)<1e-10*(1.+norm2(testoy)));
    return mat;
  }
  if (aftmdf){
    const Matrix& mat=aftmdf->apply_modified_transform(out_y,in_y,arg_trafo,arg_offset);
    //assert(norm2(out_y-testoy)<1e-10*(1.+norm2(testoy)));
    return mat;
  }

  //initialize out_y to the offset
  if (arg_offset)
    out_y=*arg_offset;
  else
    out_y.init(arg_trafo?arg_trafo->rowdim():gsmdf.old_vardim(),1,0.);
     
  //spread out the input vector on the old part to tmpvec
  if (gsmdf.map_to_old_variables()){
    if (arg_trafo==0) {
      for (Integer i=0;i<in_y.dim();i++){
	Integer ind=(*(gsmdf.map_to_old_variables()))(i);
	if (ind<gsmdf.old_vardim())
	  out_y(ind)+=in_y(i);
      }
    }
    else {
      Matrix tmpvec(gsmdf.old_vardim(),1,0.);
      for (Integer i=0;i<in_y.dim();i++){
	Integer ind=(*(gsmdf.map_to_old_variables()))(i);
	if (ind<gsmdf.old_vardim())
	  tmpvec(ind)+=in_y(i);
      }
      genmult(*arg_trafo,tmpvec,out_y,1.,1.);
    }
  }
  else {
    if (arg_trafo==0)
      out_y+=in_y;
    else 
      genmult(*arg_trafo,in_y,out_y,1.,1.);
  }
  //assert(norm2(out_y-testoy)<1e-10*(1.+norm2(testoy)));
  return out_y;
}


// *****************************************************************************
//                             transform_minorant
// *****************************************************************************

  int AffineFunctionTransformation::transform_minorant(MinorantPointer& out_minorant,
						       const MinorantPointer& in_minorant,
						       Real alpha,
						       bool add_trafo_minorant,
						       const Indexmatrix* provided_row_indices,
						       const Indexmatrix* needed_col_indices) const
  {
    if (add_trafo_minorant) {
      if (constant_minorant.get_minorant(out_minorant,alpha,0,needed_col_indices,needed_col_indices)){
	if (cb_out()){
	  get_out()<<"**** ERROR AffineFunctionTransformation::transform_minorant(...): constant_minorant.get_minorant(..) failed"<<std::endl;
	}
	return 1;
      }
    }

    alpha*=fun_coeff;

    if (in_minorant.get_minorant(out_minorant,alpha,arg_trafo,provided_row_indices,needed_col_indices,true)){
      if (cb_out()){
	get_out()<<"**** ERROR AffineFunctionTransformation::transform_minorant(...): in_minorant.get_minorant(..) failed"<<std::endl;
      }
      return 1;
    }
      
    if (arg_offset!=0){
      out_minorant.add_offset(alpha*in_minorant.ip(*arg_offset));
    }

    return 0;
  }
	

// *****************************************************************************
//                             transform_minorants
// *****************************************************************************

  int AffineFunctionTransformation::transform_minorants(MinorantBundle& out_minorants, 
						  const MinorantBundle& in_minorants, 
						  Real alpha) const
  {
    assert(out_minorants.size()==in_minorants.size());
    int err=0;
    for (unsigned int i=0;i<in_minorants.size();i++){
      int lerr=transform_minorant(out_minorants[i],in_minorants[i],alpha);
      if (lerr) {
	if (cb_out()){
	  get_out()<<"**** ERROR AffineFunctionTransformation::get_minorants(...): transform_minorant(..) failed for "<<i<<" and returned "<<lerr<<std::endl;
	}
	err+=lerr;
      }
    }
    return err;
  }
	

  
// *****************************************************************************
//                           scaling_indices
// *****************************************************************************

int
AffineFunctionTransformation::scaling_indices(Indexmatrix& row_indices,const Indexmatrix& col_indices) const
{
  if ((col_indices.dim()==0)||(arg_trafo==0)){
    row_indices=col_indices;
    return 0;
  }
  
  std::set<Integer> visited;
  const Integer n=col_indices.dim();
  for(Integer j=0;j<n;j++){
    Integer k=col_indices(j);
    Integer i;
    Integer nz=arg_trafo->col_nonzeros(k,&i);
    Integer endi=i+nz;
    for (;i<endi;i++){
      visited.insert(arg_trafo->get_colindex()(i));
    }
  }
  row_indices.newsize(Integer(visited.size()),1); chk_set_init(row_indices,1);
  Integer i=0;
  for (std::set<Integer>::const_iterator it=visited.begin();it!=visited.end();it++,i++){
    row_indices(i)=*it;
  }
    
  return 0;
}

// *****************************************************************************
//                           qp_cost_indices
// *****************************************************************************

int
AffineFunctionTransformation::qp_cost_indices(Indexmatrix& row_indices,const Indexmatrix* col_indices) const
{
  if ((col_indices!=0)&&((arg_offset==0)&&(arg_trafo==0))){
    row_indices=*col_indices;
    return 0;
  }
  if (col_indices==0){
    if ((arg_offset==0)&&(arg_trafo!=0)&&(sparse_argument_changes())){
      row_indices.init(arg_trafo->get_rowinfo().col(0));
    }
    else {
      row_indices.init(0,0,Integer(0));
    }
    return 0;
  } 
 
  std::set<Integer> visited;
  if (arg_offset){
    for(Integer i=0;i<arg_offset->dim();i++){
      if ((*arg_offset)(i)!=0.)
	visited.insert(i);
    }
  }
  if (Integer(visited.size())<to_dim()) {
    if (arg_trafo){
      const Integer n=col_indices->dim();
      for(Integer j=0;j<n;j++){
	Integer k=(*col_indices)(j);
	Integer i;
	Integer nz=arg_trafo->col_nonzeros(k,&i);
	Integer endi=i+nz;
	for (;i<endi;i++){
	  visited.insert(arg_trafo->get_colindex()(i));
	}
	if (Integer(visited.size())==to_dim()) break;
      }
    }
    else {
      const Integer n=col_indices->dim();
      for(Integer j=0;j<n;j++){
	Integer k=(*col_indices)(j);
	visited.insert(k);
      }
    }
  }

  if (Integer(visited.size())==to_dim()){
    row_indices.init(Range(0,to_dim()-1));
    return 0;
  }
  row_indices.newsize(Integer(visited.size()),1); chk_set_init(row_indices,1);
  Integer i=0;
  for (std::set<Integer>::const_iterator it=visited.begin();it!=visited.end();it++,i++){
    row_indices(i)=*it;
  }
    
  return 0;
}

// *****************************************************************************
//                          add_diagonal_scaling
// *****************************************************************************

int 
AffineFunctionTransformation::add_diagonal_scaling(
						   Matrix& diagscale,
						   const Indexmatrix* indices,
						   Real alpha,
						   const Matrix& in_diagscale
						   ) const 
{
  Real gamma=alpha*fun_coeff;
  if (gamma==0.)
    return 0;

  if (indices==0){//do it for all indices
    if (arg_trafo==0){
      diagscale.xpeya(in_diagscale,gamma);
      return 0;
    }
   
    Integer nz_cols=arg_trafo->get_colinfo().rowdim();
    for (Integer j=0;j<nz_cols;j++){
      Integer i=arg_trafo->get_colinfo()(j,2);
      Integer endi=i+arg_trafo->get_colinfo()(j,1);
      Real d=0.;
      for (;i<endi;i++){
	d+=sqr(arg_trafo->get_colval()(i))*in_diagscale(arg_trafo->get_colindex()(i));
      }
      diagscale(arg_trafo->get_colinfo()(j,0))+=d*gamma;
    }
    return 0;
  }

  if (arg_trafo==0){
    Integer k=indices->dim();
    for (Integer i=0;i<k;i++){
      Integer indi=(*indices)(i);
      diagscale(indi)+=gamma*in_diagscale(indi);
    }
    return 0;
  }

  Integer k=indices->dim();
  for(Integer j=0;j<k;j++){
    Integer indj=(*indices)(j);
    Integer i;
    Integer nz=arg_trafo->col_nonzeros(indj,&i);
    Integer endi=i+nz;
    for (;i<endi;i++){
      Real d=0.;
      for (;i<endi;i++){
	d+=sqr(arg_trafo->get_colval()(i))*in_diagscale(arg_trafo->get_colindex()(i));
      }
      diagscale(indj)+=d*gamma;
    }
  }

  return 0;
}

// *****************************************************************************
//                            apply_modification
// *****************************************************************************

int AffineFunctionTransformation::apply_modification(
						     const AFTModification* aftmdf,
						     const GroundsetModification& gsmdf
						     )
{
  if (aftmdf){
    aftmdf->apply_to_factor(fun_coeff);
    aftmdf->apply_to_offset(fun_offset);
    constant_minorant.add_offset(aftmdf->get_additional_offset());
  }

  if ((aftmdf==0)||((aftmdf->only_scalars_change())&&(!(aftmdf->ignore_groundset_modification())))){

    //-----------  default actions based on information in gsmdf 
    if (gsmdf.appended_vardim()>0){
      if (linear_cost){
	linear_cost->concat_below(Matrix(gsmdf.appended_vardim(),1,0.));
      }
      if (arg_trafo==0){  //treat as identity
	if (arg_offset)
	  arg_offset->concat_below(Matrix(gsmdf.appended_vardim(),1,0.));
      }
      else {
	arg_trafo->concat_right(Sparsemat(arg_trafo->rowdim(),gsmdf.appended_vardim()));
      }
    }
    if (gsmdf.map_to_old_variables()){
      if (linear_cost){
	*linear_cost=linear_cost->rows(*(gsmdf.map_to_old_variables()));
	Minorant* mnrt=new Minorant(true,fun_offset);
	mnrt->add_coeffs(linear_cost->dim(),linear_cost->get_store());
	//mnrt->sparsify(1e-12*norm2(*linear_cost));
	constant_minorant.init(mnrt,0);
      }
      if (arg_trafo==0){ //treat as identity
	if (arg_offset)
	  *arg_offset=arg_offset->rows(*(gsmdf.map_to_old_variables()));
      }   
      else {
	*arg_trafo=arg_trafo->cols(*(gsmdf.map_to_old_variables()));
      }
    }
    
    return 0;
  }

  int err=0;
  if ((from_dim()>=0)&&(from_dim()!=aftmdf->old_vardim())){
    if (cb_out())
      get_out()<<"**** ERROR AffineFunctionTransformation::apply_modification: this has "<<from_dim()<<" columns and modification assumes "<<aftmdf->old_vardim()<<std::endl;
    err++;
  }
  if ((to_dim()>=0)&&(to_dim()!=aftmdf->old_rowdim())){
    if (cb_out())
      get_out()<<"**** ERROR AffineFunctionTransformation::apply_modification: this has "<<to_dim()<<" rows and modification assumes "<<aftmdf->old_rowdim()<<std::endl;
    err++;
  }
  if (gsmdf.new_vardim()!=aftmdf->new_vardim()){
    if (cb_out())
      get_out()<<"**** ERROR AffineFunctionTransformation::apply_modification: local groundset will have "<<from_dim()<<" variables while modification leads to "<<aftmdf->new_vardim()<<" variables "<<std::endl;
    err++;
  }
  if (err)
    return err;

  aftmdf->apply_to_rows(arg_trafo,arg_offset);
  aftmdf->apply_to_costs(linear_cost);
  Minorant* mnrt=new Minorant(true,fun_offset);
  if (linear_cost){
    mnrt->add_coeffs(linear_cost->dim(),linear_cost->get_store());
    //mnrt->sparsify(1e-12*norm2(*linear_cost));
  }
  constant_minorant.init(mnrt,0);

  return 0;
}




// *****************************************************************************
//                            analyze_modification
// *****************************************************************************

const GroundsetModification* 
AffineFunctionTransformation::analyze_modification(
						   bool& minorant_trafo_differs,
						   const AFTModification* aftmdf,
						   const GroundsetModification& gsmdf
						   ) const
{
  minorant_trafo_differs=false;
  if (aftmdf){
    if (aftmdf->get_additional_factor()!=1.)
      minorant_trafo_differs=true;
    if (aftmdf->get_additional_offset()!=0.)
      minorant_trafo_differs=true;
  }

  if ((aftmdf==0)||((aftmdf->only_scalars_change())&&(!(aftmdf->ignore_groundset_modification())))){

    //-----------  default actions based on information in gsmdf 
    bool identity_lost=false;
    
    if (gsmdf.appended_vardim()>0){
      if (arg_trafo==0){  //treat as identity
	identity_lost=true;
      }
      //else zero columns will be appended, minorant image is preserved
    }
    if (gsmdf.map_to_old_variables()){
      if (arg_trafo==0){ //treat as identity
	minorant_trafo_differs=true;
	identity_lost=true;
      }   
      else {
	//it differs if a nonzero column is mapped somewhere else or deleted
	if (!minorant_trafo_differs) {
	  //find largest nonzero column
	  Integer maxind=-1;
	  if (arg_trafo->get_colinfo().rowdim()>0)
	    maxind=arg_trafo->get_colinfo()(arg_trafo->get_colinfo().rowdim()-1,0);
	  Integer newdim=(gsmdf.map_to_old_variables())->dim();
	  for (Integer i=0;i<newdim;i++){
	    Integer ind=(*gsmdf.map_to_old_variables())(i);
	    if ((ind!=i)&&
		(
		 ((linear_cost)&&(ind<linear_cost->dim())&&((*linear_cost)(ind)!=0.))
		 ||
		 ((linear_cost)&&(i<linear_cost->dim())&&((*linear_cost)(i)!=0.))
		 ||
		 ((ind<=maxind)&&(arg_trafo->col_nonzeros(ind)>0))
                 ||
		 ((i<=maxind)&&(arg_trafo->col_nonzeros(i)>0))
		 )){
	      minorant_trafo_differs=true;
	      break;
	    }
	  }
	}	
      }
    }
    
    if (arg_trafo==0){
      if (!identity_lost)
	return &gsmdf;
      else {
	local_gsmdf.clear(gsmdf.old_vardim());
	return &local_gsmdf;
      }
    }
  
    local_gsmdf.clear(arg_trafo->rowdim());
    return &local_gsmdf;
  }

  int err=0;
  if ((from_dim()>=0)&&(from_dim()!=aftmdf->old_vardim())){
    if (cb_out())
      get_out()<<"**** ERROR AffineFunctionTransformation::apply_modification: this has "<<from_dim()<<" columns and modification assumes "<<aftmdf->old_vardim()<<std::endl;
    err++;
  }
  if ((to_dim()>=0)&&(to_dim()!=aftmdf->old_rowdim())){
    if (cb_out())
      get_out()<<"**** ERROR AffineFunctionTransformation::apply_modification: this has "<<to_dim()<<" rows and modification assumes "<<aftmdf->old_rowdim()<<std::endl;
    err++;
  }
  if (gsmdf.new_vardim()!=aftmdf->new_vardim()){
    if (cb_out())
      get_out()<<"**** ERROR AffineFunctionTransformation::apply_modification: local groundset will have "<<from_dim()<<" variables while modification leads to "<<aftmdf->new_vardim()<<" variables "<<std::endl;
    err++;
  }
  if (err)
    return 0;

  bool gsmdf_determines_groundset=false;
  if ((arg_trafo==0)&&(aftmdf->preserves_identity())){   
    gsmdf_determines_groundset=true;
    if (aftmdf->appended_vardim()!=gsmdf.appended_vardim()){
      gsmdf_determines_groundset=false;
      if (cb_out(1))
	get_out()<<"**** WARNING AffineFunctionTransformation::apply_modification: arg_trafo stays NULL but modification does not behave like the identity: modification appends "<<aftmdf->appended_vardim()<<" variables to AFT and "<<aftmdf->appended_rowdim()<<" rows to AFT and local groundset appended "<<gsmdf.appended_vardim()<<std::endl;
    }
    if ((aftmdf->map_to_old_variables()==0)!=(gsmdf.map_to_old_variables()==0)){
      gsmdf_determines_groundset=false;
      if (cb_out(1))
	get_out()<<"**** WARNING AffineFunctionTransformation::apply_modification: arg_trafo stays NULL but modification does not behave like the identity: modification reorders variables differently than in the local groundset"<<std::endl;
    }
    else if ((gsmdf_determines_groundset)&&(aftmdf->map_to_old_variables())){
      Integer n=aftmdf->map_to_old_variables()->dim();
      Integer i=0;
      while ((i<n)&&((*(gsmdf.map_to_old_variables()))(i)==(*(aftmdf->map_to_old_variables()))(i)))
	  i++;
      if (i<n){
	gsmdf_determines_groundset=false;
	if (cb_out(1))
	  get_out()<<"**** WARNING AffineFunctionTransformation::apply_modification: arg_trafo stays NULL but modification does not behave like the identity: modification reorders variables differently than in the local groundset"<<std::endl;
      }
    }
  }	

    
  if (gsmdf_determines_groundset){
    if(!(gsmdf.no_additions_or_deletions_in_vars()))
      minorant_trafo_differs=true;
  }

  
  //transformation changes when appending nonzero columns
  if ((!minorant_trafo_differs)&&((aftmdf->get_append_cols())||(aftmdf->get_append_costs()))){
    minorant_trafo_differs=true;
  }

  //check wether reordering of the variables may change the nonzeros
  if ((!minorant_trafo_differs)&&(aftmdf->map_to_old_variables())){
    //nothing is appended, because this was checked before already
    if (arg_trafo==0){ //treat as identity 
      minorant_trafo_differs=true;
    }   
    else {
      //it differs if a nonzero column is mapped somewhere else or deleted
      //find largest nonzero column
      Integer maxind=-1;
      if (arg_trafo->get_colinfo().rowdim()>0)
	maxind=arg_trafo->get_colinfo()(arg_trafo->get_colinfo().rowdim()-1,0);
      Integer newdim=(gsmdf.map_to_old_variables())->dim();
      for (Integer i=0;i<newdim;i++){
	Integer ind=(*(aftmdf->map_to_old_variables()))(i);
	if ((ind!=i)&&
	    (
	     ((linear_cost)&&(ind<linear_cost->dim())&&((*linear_cost)(ind)!=0.))
	     ||
	     ((linear_cost)&&(i<linear_cost->dim())&&((*linear_cost)(i)!=0.))
	     ||
	     ((ind<=maxind)&&(arg_trafo->col_nonzeros(ind)>0))
	     ||
	     ((i<=maxind)&&(arg_trafo->col_nonzeros(i)>0))
	     )){
	  minorant_trafo_differs=true;
	  break;
	}
      }
    }	
  }

  //transformation changes when appending nonzero rows
  if ((!minorant_trafo_differs)&&((aftmdf->get_append_rows())||(aftmdf->get_append_rhs()))){
    minorant_trafo_differs=true;
  }

  //check wether reordering of the rows may change the nonzeros
  if ((!minorant_trafo_differs)&&(aftmdf->map_to_old_rows())){
    //nothing is appended, because this was checked before already
    if (arg_trafo==0){ //treat as identity 
      minorant_trafo_differs=true;
    }   
    else {
      //it differs if a nonzero column is mapped somewhere else or deleted
      //find largest nonzero column
      Integer maxind=-1;
      if (arg_trafo->get_rowinfo().rowdim()>0)
	maxind=arg_trafo->get_rowinfo()(arg_trafo->get_rowinfo().rowdim()-1,0);
      Integer newdim=aftmdf->map_to_old_rows()->dim();
      for (Integer i=0;i<newdim;i++){
	Integer ind=(*(aftmdf->map_to_old_rows()))(i);
	if ((ind!=i)&&
	    (
	     ((arg_offset)&&(ind<arg_offset->dim())&&((*arg_offset)(ind)!=0.))
	     ||
	     ((arg_offset)&&(i<arg_offset->dim())&&((*arg_offset)(i)!=0.))
	     ||
	     ((ind<=maxind)&&(arg_trafo->row_nonzeros(ind)>0))
	     ||
	     ((i<=maxind)&&(arg_trafo->row_nonzeros(i)>0))
	     )){
	  minorant_trafo_differs=true;
	  break;
	}
      }
    }	
  }


    
  if (gsmdf_determines_groundset){
    return &gsmdf;
  }
  
  local_gsmdf.clear(aftmdf->old_rowdim());
  if (aftmdf->appended_rowdim()>=0){
    local_gsmdf.add_append_vars(aftmdf->appended_rowdim());
  }
  if (aftmdf->map_to_old_rows()){
    local_gsmdf.add_reassign_vars(*(aftmdf->map_to_old_rows()));
  }
  return &local_gsmdf;
}

// *****************************************************************************
//                               output_aft_data
// *****************************************************************************

std::ostream& AffineFunctionTransformation::output_aft_data(std::ostream& out) const
{ 
  out<<"\n%(begin aft)\n";
  out.precision(16);
  out.width(18);
  out<<" fun_coeff="<<fun_coeff<<";";
  out.width(18);
  out<<"\n fun_offset="<<fun_offset<<";";
  out<<"\n linear_cost=[";
    if (linear_cost){
      for (Integer i=0;i<linear_cost->dim();i++){
	out.width(18);
	out<<" "<<(*linear_cost)(i);
      }
    }
  out<<"]';";
  out<<"\n arg_offset=[";
    if (linear_cost){
      for (Integer i=0;i<arg_offset->dim();i++){
	out.width(18);
	out<<" "<<(*arg_offset)(i);
      }
    }
  out<<"]';";
  out<<"\n% arg_trafo\n";
  Indexmatrix indi,indj;
  Matrix val;
  if (arg_trafo){
    arg_trafo->get_edge_rep(indi,indj,val);
    out<<"\n indi=[";
    for (Integer i=0;i<indi.dim();i++){
      out.width(18);
      out<<" "<<indi(i);
    }
    out<<"\n indj=[";
    for (Integer i=0;i<indj.dim();i++){
      out.width(18);
      out<<" "<<indj(i);
    }
    out<<"]';";
    out<<"\n val=[";
    for (Integer i=0;i<val.dim();i++){
      out.width(18);
      out<<" "<<val(i);
    }
    out<<"]';";
    out<<"\n arg_trafo=sparse(indi,indj,val,";
    out<<arg_trafo->rowdim()<<","<<arg_trafo->coldim()<<");";
  }
  else {
    out<<"\n arg_trafo=[];";
  }
 
  out<<"\n%(end aft)\n";
  return out;
}


}
