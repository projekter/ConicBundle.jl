/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/AFTModification.cxx
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



#include "CBSolver.hxx"
#include "AFTModification.hxx"

using namespace CH_Matrix_Classes;

namespace ConicBundle {


  AFTModification::~AFTModification()
  {}

  AFTModification::AFTModification(Integer var_olddim,Integer row_olddim,bool ignore_groundset_mdf):mdf(var_olddim,row_olddim),ignore_gs_mdf(ignore_groundset_mdf)
  {factor=1;offset=0.;preserves_id_flag=-1;}

  void AFTModification::clear(Integer var_olddim,Integer row_olddim)
  {mdf.clear(var_olddim,row_olddim);factor=1.;offset=0.;preserves_id_flag=-1;ignore_gs_mdf=false;}

  int AFTModification::incorporate(const AFTModification& m)
  {
    factor*=m.factor;
    offset+=m.offset;
    ignore_gs_mdf = ignore_gs_mdf || m.ignore_gs_mdf;
    preserves_id_flag=-1;
    return mdf.incorporate(m.mdf);
  }
  
  int AFTModification::apply_to_costs(CH_Matrix_Classes::Matrix*& linear_cost) const
  {
    if ((linear_cost==0)&&(mdf.get_var_append_costs()==0))
      return 0;

    if (linear_cost==0)
      linear_cost=new Matrix(mdf.old_vardim(),1,0.);

    return mdf.apply_to_vars(0,0,0,linear_cost);
  }

  int AFTModification::apply_to_rows(CH_Matrix_Classes::Sparsemat*& rows,
		    CH_Matrix_Classes::Matrix*& rhs) const
  {
    if ((rhs==0)&&(mdf.get_row_append_rhslb()))
      rhs=new Matrix(mdf.old_rowdim(),1,0.);

    if ((rows==0)&&(!preserves_identity())){
      Integer dim=mdf.old_rowdim();
      Indexmatrix indi(Range(0,dim-1));
      rows=new Sparsemat(dim,dim,dim,indi,indi,Matrix(dim,1,1.));
    }
  
    int retval= mdf.apply_to_rows(rows,rhs,0);

    if ((retval==0)&&(rhs)&&(mdf.new_row_indices())){
      for (Integer i=0;i<mdf.new_row_indices()->rowdim();i++){
	Integer ind=(*(mdf.new_row_indices()))(i);
	if ((*rhs)(ind)<=CB_minus_infinity)
	  (*rhs)(ind)=0.;
      }
    }

    return retval;
  }

  const Matrix& AFTModification::apply_modified_transform(Matrix& out_y,
						const Matrix& in_y,
						const Sparsemat* arg_trafo,
						const Matrix* arg_offset) const
  {
    //--------- the case of the identity transformation: find the right offset
    if ((arg_trafo==0)&&(preserves_identity())){
      assert(in_y.dim()==new_vardim());
      if ((arg_offset==0)&&(get_append_rhs()==0)){
	return in_y;
      }
      if (arg_offset){
	//initialize by arg_offset and append further rows (or zero)
	assert(arg_offset->dim()==old_rowdim());
	out_y.init(*arg_offset);
	if (get_append_rhs()){
	  out_y.concat_below(*(get_append_rhs()));
	}
	else {
	  out_y.concat_below(Matrix(appended_rowdim(),1,0.));
	}
      }
      else {
	if (get_append_rhs()){
	  out_y.init(old_rowdim(),1,0.);
	  out_y.concat_below(*(get_append_rhs()));
	}
	else 
	  out_y.init(old_rowdim()+appended_rowdim(),1,0.);
      }
      if (map_to_old_rows())
	out_y=out_y(*(map_to_old_rows()));
      out_y+=in_y;
      return out_y;
    }
    
    //---------- the non identity case
    //spread out the input vector to invec1 (old columns) and invec2 (appended ones)
    Matrix invec1,invec2;
    if (map_to_old_variables()){
      invec1.init(old_vardim(),1,0.);
      invec2.init(appended_vardim(),1,0.);
      for (Integer i=0;i<in_y.dim();i++){
	Integer ind=(*(map_to_old_variables()))(i);
	if (ind>=old_vardim())
	  invec2(ind-old_vardim())=in_y(i);
	else
	  invec1(ind)=in_y(i);
      }
    }
    else {
      invec1=in_y(Range(0,old_vardim()-1));
      invec2=in_y(Range(old_vardim(),old_vardim()+appended_vardim()-1));
    }
    //tansform the part of the old columns in out_y
    if (arg_trafo==0) 
      out_y=invec1;
    else 
      genmult(*arg_trafo,invec1,out_y);
    //add the part of the new columns to out_y
    if (get_append_cols())
      genmult(*get_append_cols(),invec2,out_y,1.,1.);
    //add the old offset to out_y
    if (arg_offset)
      out_y+=*arg_offset;
    //form the rows to be appended to out_y in invec1
    if (get_append_rhs())
      invec1=*(get_append_rhs());
    else 
      invec1.init(appended_rowdim(),1,0.);
    //add the part of the new rows (we may use in_y now)
    if (get_append_rows()){
      genmult(*get_append_rows(),in_y,invec1,1.,1.);
    }
    out_y.concat_below(invec1);
    //reorder out_y
    if (map_to_old_rows())
      out_y=out_y(*(map_to_old_rows()));
    return out_y;
  }
  
  bool AFTModification::groundset_changes_suffice_for_identity()
  {
    if ((factor!=1.)||(offset!=0.))
      return false;
    if (get_append_costs()!=0)
      return false;
    if (get_append_rhs()!=0)
      return false;
    return preserves_identity();
  }
    

  bool AFTModification::preserves_identity() const
  {
    if (preserves_id_flag>=0){
      return (preserves_id_flag==1);
    }
    if ((mdf.old_vardim()!=mdf.old_rowdim())||(mdf.new_vardim()!=mdf.new_rowdim())){
      preserves_id_flag=0;
      return false;
    }

    if (mdf.appended_vardim()!=mdf.appended_rowdim()){
      preserves_id_flag=0;
      return false;
    }

    if ((mdf.map_to_old_variables()==0)!=(mdf.map_to_old_rows()==0)){
      preserves_id_flag=0;
      return false;
    }

    if (mdf.map_to_old_variables()){
      Integer n=mdf.map_to_old_variables()->dim();
      Integer i=0;
      while ((i<n)&&((*(mdf.map_to_old_variables()))(i)==(*(mdf.map_to_old_rows()))(i)))
	  i++;
      if (i<n){
	preserves_id_flag=0;
	return false;
      }
    }	

    if ((mdf.get_var_append_cols()!=0)||(mdf.get_row_append_mat()!=0)){
      if (mdf.get_var_append_cols()->nonzeros()!=0){
	preserves_id_flag=0;
	return false;
      }
      if (mdf.get_row_append_mat()->nonzeros()!=mdf.appended_rowdim()){
	preserves_id_flag=0;
	return false;
      }
      for (Integer i=0;i<mdf.appended_rowdim();i++){
	Integer startind;
	if (mdf.get_row_append_mat()->row_nonzeros(i,&startind)!=1){
	  preserves_id_flag=0;
	  return false;
	}
	if (mdf.get_row_append_mat()->get_rowval()(startind)!=1.){
	  preserves_id_flag=0;
	  return false;
	}
	Integer colind=mdf.get_row_append_mat()->get_rowindex()(startind);
	if (mdf.map_to_old_rows()==0){
	  if (colind!=mdf.old_rowdim()+i){
	    preserves_id_flag=0;
	    return false;
	  }
	}
	else if ((*(mdf.map_to_old_rows()))(colind)!=mdf.old_rowdim()+1)
	{
	  preserves_id_flag=0;
	  return false;
	}
      }
    }

    preserves_id_flag=1;
    return true;
  }

  bool AFTModification::no_additions_or_deletions_in_vars() const
  {
    return ((mdf.appended_vardim()==0)&&((mdf.map_to_old_variables()==0)||(mdf.map_to_old_variables()->dim()==mdf.old_vardim())));
  }

  bool AFTModification::no_additions_or_deletions_in_rows() const
  {
    return ((mdf.appended_rowdim()==0)&&((mdf.map_to_old_rows()==0)||(mdf.map_to_old_rows()->dim()==mdf.old_rowdim())));
  }

}

