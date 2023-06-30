/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/SparseCoeffmatMatrix.cxx
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



#include <string.h>
#include <stdlib.h>
#include <cctype>

#include "SparseCoeffmatMatrix.hxx"
#include "PSCPrimal.hxx"

using namespace CH_Matrix_Classes;

namespace ConicBundle {


//******************************************************************************
//                  SparseCoeffmatMatrix::SparseCoeffmatMatrix
//******************************************************************************

void SparseCoeffmatMatrix::form_colrep() const
{
  if (colrep!=0)
    return;
  
  colrep=new SCMcolrep;
  for (unsigned i=0;i<blockrep.size();i++){
    const SparseCoeffmatVector& sv=blockrep[i];
    for (SparseCoeffmatVector::const_iterator it=sv.begin();it!=sv.end();it++){
      assert(it->second!=0);
      (*colrep)[it->first][Integer(i)]=it->second;
    }
  }
  
  return;
}

//******************************************************************************
//                  SparseCoeffmatMatrix::~SparseCoeffmatMatrix
//******************************************************************************
    
  SparseCoeffmatMatrix::~SparseCoeffmatMatrix()
  {
    clear();
  }

//******************************************************************************
//                  SparseCoeffmatMatrix::clear
//******************************************************************************
  
void SparseCoeffmatMatrix::clear()
{
  block_dim.init(0,1,Integer(0));
  dense_cnt.init(0,1,Integer(0));
  col_dim=0;

  delete colrep;
  colrep=0;

  blockrep.clear();
}

//******************************************************************************
//                  SparseCoeffmatMatrix::operator=
//******************************************************************************
    
  SparseCoeffmatMatrix& SparseCoeffmatMatrix::operator=(const SparseCoeffmatMatrix& S)
  {
    clear();
    block_dim=S.block_dim;
    dense_cnt=S.dense_cnt;
    col_dim=S.col_dim;
    blockrep=S.blockrep;
    return *this;
  }

//******************************************************************************
//                  SparseCoeffmatMatrix::init
//******************************************************************************
  
int SparseCoeffmatMatrix::init(const Indexmatrix& in_block_dim,
			       Integer in_col_dim,
			       const Indexmatrix* block_ind,
			       const Indexmatrix* col_ind,
			       const CoeffmatVector* coeff_vec)
{
  clear();
  int err=0;
  if (in_block_dim.coldim()!=1){
    if (cb_out())
      get_out()<<"**** ERROR: SparseCoeffmatMatrix::init(.....): in_block_dim.coldim()="<<block_dim.coldim()<<" but should equal to 1 (even with 0 rows), restructering as column vector"<<std::endl;
    err++;
  }
  block_dim.init(in_block_dim.dim(),1,in_block_dim.get_store());
  if (in_col_dim<0){
    if (cb_out())
      get_out()<<"**** ERROR: SparseCoeffmatMatrix::init(.....): in_col_dim="<<in_col_dim<<" but should be nonnegative, setting to zero"<<std::endl;
    err++;
  }
  col_dim=max(0,in_col_dim);
  for (Integer i=0;i<block_dim.dim();i++){
    if (block_dim(i)<0){
      if (cb_out())
	get_out()<<"**** ERROR: SparseCoeffmatMatrix::init(.....): block_dim(i)="<<block_dim(i)<<" but should be nonnegative, setting to 0"<<std::endl;
      err++;
      block_dim(i)=0;
    }
  }

  dense_cnt.init(in_block_dim.dim(),1,Integer(0));
  blockrep.resize(unsigned(block_dim.dim()));
  if (((block_ind==0)!=(col_ind==0))||
      ((block_ind==0)!=(coeff_vec==0))){
      if (cb_out())
	get_out()<<"**** ERROR: SparseCoeffmatMatrix::init(.....): in the sparse representation of block_ind, col_ind, coeff_vec some are zero, some not; either all or none must be given, skipping all"<<std::endl;
      err++;
      return err;
  }

  if (block_ind!=0){
    if ((block_ind->dim()!=col_ind->dim())||(block_ind->dim()!=Integer(coeff_vec->size()))){
      if (cb_out())
	get_out()<<"**** ERROR: SparseCoeffmatMatrix::init(.....): in the sparse representation the size of block_ind, col_ind, coeff_vec is not identical but should be, skipping all"<<std::endl;
      err++;
      return err;
    }
     
    for (Integer i=0;i<block_ind->dim();i++){
      Integer indi=(*block_ind)(i);
      Integer indj=(*col_ind)(i);
      CoeffmatPointer cm=(*coeff_vec)[unsigned(i)];
      if (set(indi,indj,cm)){
	if (cb_out())
	  get_out()<<"**** ERROR: SparseCoeffmatMatrix::init(.....): inserting item "<<i<<" failed"<<std::endl;
	err++;
      }
    }
  }
  return err;
}
  
//******************************************************************************
//                  SparseCoeffmatMatrix::operator()
//******************************************************************************
  
Integer SparseCoeffmatMatrix::blockdim(Integer i) const
{
  if ((i<0)||(i>=block_dim.dim())){
    if (cb_out())
      get_out()<<"**** ERROR: SparseCoeffmatMatrix::blockdim(.): block index i="<<i<<" exceeds the range [0,"<<block_dim.dim()-1<<"]"<<std::endl;
    return -1;
  }
  return block_dim(i);
}
  
//******************************************************************************
//                  SparseCoeffmatMatrix::set
//******************************************************************************
  
  int SparseCoeffmatMatrix::set(Integer i,Integer j,const CoeffmatPointer& cm) 
{
  int err=0;
  if ((i<0)||(i>=block_dim.dim())){
      if (cb_out())
	get_out()<<"**** ERROR: SparseCoeffmatMatrix::set(i,j,cm): block index i="<<i<<" exceeds the range [0,"<<block_dim.dim()-1<<"]"<<std::endl;
      err++;
  }
  if ((j<0)||(j>=col_dim)){
    if (cb_out())
      get_out()<<"**** ERROR: SparseCoeffmatMatrix::set(i,j,cm): block index j="<<j<<" exceeds the range [0,"<<col_dim-1<<"]"<<std::endl;
    err++;
  }
  if (err)
    return err;
  
  delete colrep;
  colrep=0;

  if (cm==0){
    SparseCoeffmatVector::iterator it=blockrep[unsigned(i)].find(j);
    if (it==blockrep[unsigned(i)].end())
      return err;
    if (it->second->dense()){
      assert(dense_cnt(i)>0);
      dense_cnt(i)--;
    }
    blockrep[unsigned(i)].erase(it);
    return err;
  }

  if (cm->dim()!=block_dim(i)){
    if (cb_out())
      get_out()<<"**** ERROR: SparseCoeffmatMatrix::set(.....): the order of block index i="<<i<<" is "<<block_dim(i)<<", but the coefficient matrix has order "<<cm->dim()<<";  replacing anyway"<<std::endl;
    err++;
  } 
 
  if (cm->dense()){
    dense_cnt(i)++;
  }

  SparseCoeffmatVector::iterator it=blockrep[unsigned(i)].find(j);
  if (it!=blockrep[unsigned(i)].end()){
    if (it->second->dense()){
      assert(dense_cnt(i)>0);
      dense_cnt(i)--;
    }
    it->second=cm;
  }
  else {
    (blockrep[unsigned(i)])[j]=cm;
  }
  
  return err;
}

//******************************************************************************
//                  SparseCoeffmatMatrix::operator()
//******************************************************************************
  
CoeffmatPointer SparseCoeffmatMatrix::operator()(Integer i,Integer j) const
{
  int err=0;
  if ((i<0)||(i>=block_dim.dim())){
    if (cb_out())
      get_out()<<"**** ERROR: SparseCoeffmatMatrix::operator(i,j): block index i="<<i<<" exceeds the range [0,"<<block_dim.dim()-1<<"]"<<std::endl;
    err++;
  }
  if ((j<0)||(j>=col_dim)){
    if (cb_out())
      get_out()<<"**** ERROR: SparseCoeffmatMatrix::operator(i,j): block index j="<<i<<" exceeds the range [0,"<<col_dim-1<<"]"<<std::endl;
    err++;
  }
  if ((err)||(blockrep[unsigned(i)].size()==0))
    return 0;
  SparseCoeffmatVector::const_iterator it=blockrep[unsigned(i)].find(j);
  if (it==blockrep[unsigned(i)].end())
    return 0;
  return it->second;
}

//******************************************************************************
//                  SparseCoeffmatMatrix::block
//******************************************************************************
 
const SparseCoeffmatVector* SparseCoeffmatMatrix::block(Integer i) const
{
  if ((i<0)||(i>=block_dim.dim())){
    if (cb_out())
      get_out()<<"**** ERROR: SparseCoeffmatMatrix::block(.): block index i="<<i<<" exceeds the range [0,"<<block_dim.dim()-1<<"]"<<std::endl;
    return 0;
  }
  return (blockrep[unsigned(i)].size()==0)?0:&(blockrep[unsigned(i)]);
}

//******************************************************************************
//                  SparseCoeffmatMatrix::column
//******************************************************************************
 
const SparseCoeffmatVector* SparseCoeffmatMatrix::column(Integer i) const
{
  if ((i<0)||(i>=col_dim)){
    if (cb_out())
      get_out()<<"**** ERROR: SparseCoeffmatMatrix::block(.): block index i="<<i<<" exceeds the range [0,"<<col_dim-1<<"]"<<std::endl;
    return 0;
  }
  form_colrep();
  SCMcolrep::const_iterator it=colrep->find(i);
  if (it!=colrep->end())
    return &(it->second);
  return 0;
}

//******************************************************************************
//                  SparseCoeffmatMatrix::append_blocks
//******************************************************************************
 
int SparseCoeffmatMatrix::append_blocks(const SparseCoeffmatMatrix& append_mat,
					const Indexmatrix* blocks,
					const Indexmatrix* cols)
{
  int err=0;
  if (blocks){
    if (max(*blocks)>append_mat.blockdim().dim()){
      if (cb_out())
	get_out()<<"**** ERROR: SparseCoeffmatMatrix::append_blocks(...): indices within *blocks exceed range of append_mat"<<std::endl;
      err++;
    }
  }
  if (cols){
    if (max(*cols)>append_mat.coldim()){
      if (cb_out())
	get_out()<<"**** ERROR: SparseCoeffmatMatrix::append_blocks(...): indices within *cols exceed range of append_mat"<<std::endl;
      err++;
    }
    if (((block_dim.dim()>0)||(col_dim>0))&&(cols->dim()!=col_dim)){
      if (cb_out())
	get_out()<<"**** ERROR: SparseCoeffmatMatrix::append_blocks(...): dimensoin of *cols ="<<cols->dim()<<" != "<<col_dim<<"= number of columns of this"<<std::endl;
      err++;
    }
  }
  else {
    if (((block_dim.dim()>0)||(col_dim>0))&&(col_dim!=append_mat.col_dim)){
      if (cb_out())
	get_out()<<"**** ERROR: SparseCoeffmatMatrix::append_blocks(...): col_dim="<<col_dim<<" does not match appended col_dim="<<append_mat.col_dim<<", skipping append"<<std::endl;
      err++;
    }
  }

  if (err){
    if (cb_out())
      get_out()<<"**** ERROR: SparseCoeffmatMatrix::append_blocks(...): fatal errors on input paramters, skipping append; return value="<<err<<std::endl;
    return err;
  }

  delete colrep;
  colrep=0;

  if (col_dim==0){
    if (cols)
      col_dim=cols->dim();
    else
      col_dim=append_mat.col_dim;
  }

  Integer olddim=block_dim.dim();
  if (blocks) {
    if (block_dim.coldim()==0){
      block_dim.init(blocks->dim(),1,Integer(0));
      dense_cnt.init(blocks->dim(),1,Integer(0));
    }
    else {
      block_dim.enlarge_below(blocks->dim(),Integer(0));
      dense_cnt.enlarge_below(blocks->dim(),Integer(0));
    }
    unsigned s=unsigned(blockrep.size());
    blockrep.resize(unsigned(block_dim.dim()));
    for(Integer i=0;i<blocks->dim();i++){
      Integer iind=(*blocks)(i);
      if (iind<0){
	block_dim(olddim+i)=-iind;
      }
      else {
	block_dim(olddim+i)=append_mat.block_dim(iind);
	if (append_mat.blockrep[unsigned(iind)].size()>0){
	  if (cols){
	    for(Integer j=0;j<col_dim;j++){
	      Integer jind=(*cols)(j);
	      if (jind>=0){
		SparseCoeffmatVector::const_iterator it=append_mat.blockrep[unsigned(iind)].find(jind);
		if (it!=append_mat.blockrep[unsigned(iind)].end()){
		  (blockrep[s+unsigned(i)])[jind]=it->second;
		  if (it->second->dense()){
		    dense_cnt(Integer(s)+i)++;
		  }
		}
	      }
	    }
	  }
	  else {
	    blockrep[s+unsigned(i)] = append_mat.blockrep[unsigned(iind)];
	    dense_cnt(Integer(s)+i)+= append_mat.get_dense_cnt(iind);
	  }
	}
      }
    }
  }
  else { //blocks==0
    block_dim.concat_below(append_mat.block_dim);
    dense_cnt.enlarge_below(append_mat.rowdim(),Integer(0));
    unsigned s=unsigned(blockrep.size());
    blockrep.resize(s+append_mat.blockrep.size());
    for(unsigned i=0;i<append_mat.blockrep.size();i++){
      if (append_mat.blockrep[i].size()>0){
	if (cols){
	  for(Integer j=0;j<col_dim;j++){
	    Integer jind=(*cols)(j);
	    if (jind>=0){
	      SparseCoeffmatVector::const_iterator it=append_mat.blockrep[i].find(jind);
	      if (it!=append_mat.blockrep[i].end()){
		(blockrep[s+i])[jind]=it->second;
		if (it->second->dense())
		  dense_cnt(Integer(s+i))++;
	      }
	    }
	  }
	}
	else {
	  blockrep[s+i] = append_mat.blockrep[i];
	  dense_cnt(Integer(s+i))+=append_mat.get_dense_cnt(Integer(i));
	}
      }
    }
  }
 
  return 0;
}

//******************************************************************************
//                  SparseCoeffmatMatrix::append_columns
//******************************************************************************
 
int SparseCoeffmatMatrix::append_columns(const SparseCoeffmatMatrix& append_mat,
					 const Indexmatrix* blocks,
					 const Indexmatrix* cols)
{
  int err=0;
  if (cols){
    if (max(*cols)>append_mat.coldim()){
      if (cb_out())
	get_out()<<"**** ERROR: SparseCoeffmatMatrix::append_columns(...): indices within *cols exceed range of append_mat"<<std::endl;
      err++;
    }
  }
  if (blocks){
    if (max(*blocks)>append_mat.blockdim().dim()){
      if (cb_out())
	get_out()<<"**** ERROR: SparseCoeffmatMatrix::append_columns(...): indices within *blocks exceed range of append_mat"<<std::endl;
      err++;
    }
    if ((block_dim.dim()>0)||(col_dim>0)){
      if (blocks->dim()!=block_dim.dim()){
	if (cb_out())
	  get_out()<<"**** ERROR: SparseCoeffmatMatrix::append_columns(...): blocks->dim()="<<blocks->dim()<<" != "<<block_dim.dim()<<"= number of blocks in this"<<std::endl;
	err++;
      }
      else {
	for(Integer i=0;i<block_dim.dim();i++){
	  Integer d=(*blocks)(i);
	  if (d<0){
	    if (block_dim(i)!=-d){
	      if (cb_out())
		get_out()<<"**** ERROR: SparseCoeffmatMatrix::append_columns(...): -(*blocks)("<<i<<")="<<-d<<" != "<<block_dim(i)<<"= order of block "<<i<<" in this"<<std::endl;
	      err++;
	    }
	  }
	  else {
	    if (block_dim(i)!=append_mat.block_dim(d)){
	      if (cb_out())
		get_out()<<"**** ERROR: SparseCoeffmatMatrix::append_columns(...): append_mat.block_dim((*blocks)("<<i<<"))="<<append_mat.block_dim(d)<<" != "<<block_dim(i)<<"= order of block "<<i<<" in this"<<std::endl;
	      err++;
	    }
	  }
	}
      }
    }
  }
  else {
    if ((block_dim.dim()>0)||(col_dim>0)){
      if (append_mat.block_dim.dim()!=block_dim.dim()){
	if (cb_out())
	  get_out()<<"**** ERROR: SparseCoeffmatMatrix::append_columns(...): append_mat.block_dim.dim()="<<append_mat.block_dim.dim()<<" != "<<block_dim.dim()<<"= number of blocks in this"<<std::endl;
	err++;
      }
      else {
	for(Integer i=0;i<block_dim.dim();i++){
	  if (block_dim(i)!=append_mat.block_dim(i)){
	    if (cb_out())
	      get_out()<<"**** ERROR: SparseCoeffmatMatrix::append_columns(...): append_mat.block_dim("<<i<<")="<<append_mat.block_dim(i)<<" != "<<block_dim(i)<<"= order of block "<<i<<" in this"<<std::endl;
	    err++;
	  }
	}
      }
    }
  }

  if (err){
    if (cb_out())
      get_out()<<"**** ERROR: SparseCoeffmatMatrix::append_columns(...): fatal errors on input paramters, skipping append; return value="<<err<<std::endl;
    return err;
  }

  delete colrep;
  colrep=0;

  if (block_dim.dim()==0){
    if (blocks){
      block_dim.init(blocks->dim(),1,Integer(0));
      dense_cnt.init(blocks->dim(),1,Integer(0));
    }
    else {
      block_dim=append_mat.block_dim;
      dense_cnt.init(block_dim.dim(),1,Integer(0));
    }
    blockrep.resize(unsigned(block_dim.dim()));
  }

  Integer olddim=col_dim;
  col_dim+=cols? cols->dim():append_mat.col_dim;
  if (blocks){
    for (Integer i=0;i<block_dim.dim();i++){
      Integer d=(*blocks)(i);
      if (d<0){
	block_dim(i)=-d;
      }
      else {
	block_dim(i)=append_mat.block_dim(d);
	const SparseCoeffmatVector& sv=append_mat.blockrep[unsigned(d)];
	if (sv.size()>0){
	  if (cols){
	    for(Integer j=0;j<cols->dim();j++){
	      Integer jind=(*cols)(j);
	      if (jind>=0){
		SparseCoeffmatVector::const_iterator it=sv.find(jind);
		if (it!=sv.end())
		  err+=set(i,olddim+j,it->second);
	      }
	    }
	  }
	  else {
	    for(SparseCoeffmatVector::const_iterator it=sv.begin();it!=sv.end();it++){
	      err+=set(i,olddim+it->first,it->second);
	    }
	  }
	}
      }
    }
  }
  else {
    for (unsigned i=0;i<append_mat.blockrep.size();i++){
      const SparseCoeffmatVector& sv=append_mat.blockrep[i];
      if (sv.size()>0){
	if (cols){
	  for(Integer j=0;j<cols->dim();j++){
	    Integer jind=(*cols)(j);
	    if (jind>=0){
	      SparseCoeffmatVector::const_iterator it=sv.find((*cols)(j));
	      if (it!=sv.end())
		err+=set(Integer(i),olddim+j,it->second);
	    }
	  }
	}
	else {
	  for(SparseCoeffmatVector::const_iterator it=sv.begin();it!=sv.end();it++){
	    err+=set(Integer(i),olddim+it->first,it->second);
	  }
	}
      }
    }
  }
  if (err)
    if (cb_out())
      get_out()<<"**** ERROR: SparseCoeffmatMatrix::append_columns(.): errors occured in appending the coefficient matrices by set()"<<std::endl;
 
  return err;
}

//******************************************************************************
//                  SparseCoeffmatMatrix::reassign_blocks
//******************************************************************************
 
int SparseCoeffmatMatrix::reassign_blocks(const Indexmatrix& map_to_old)
{
  if ((map_to_old.dim()>0)&&((min(map_to_old)<0)||(max(map_to_old)>=block_dim.dim()))){
    if (cb_out())
      get_out()<<"**** ERROR: SparseCoeffmatMatrix::reassign_blocks(.): some map entries exceed the range [0,"<<block_dim.dim()-1<<"]"<<std::endl;
    return 1;
  }

  delete colrep;
  colrep=0;

  Indexmatrix new_blocks(block_dim(map_to_old));
  swap(block_dim,new_blocks);
  new_blocks.init(new_blocks.dim(),1,-1);
  SCMblockrep new_rep(unsigned(map_to_old.dim()));
  
  int err=0;
  for (Integer i=0;i<map_to_old.dim();i++){
    Integer ind=map_to_old(i);
    if (new_blocks(ind)>=0){
      if (cb_out())
	get_out()<<"**** ERROR: SparseCoeffmatMatrix::reassign_blocks(.): map entry "<<i<<" refers again to "<<ind<<" (first occurence was map entry "<<new_blocks(ind)<<"), assigning zero instead"<<std::endl;
      err++;
      continue;
    }
    new_blocks(ind)=i;
    swap(new_rep[unsigned(i)],blockrep[unsigned(ind)]);
  }

  swap(new_rep,blockrep);
  dense_cnt=dense_cnt(map_to_old);
  return err;
}


//******************************************************************************
//                  SparseCoeffmatMatrix::reassign_columns
//******************************************************************************
 
int SparseCoeffmatMatrix::reassign_columns(const Indexmatrix& map_to_old)
{
  if ((map_to_old.dim()>0)&&((min(map_to_old)<0)||(max(map_to_old)>=col_dim))){
    if (cb_out())
      get_out()<<"**** ERROR: SparseCoeffmatMatrix::reassign_columns(.): some map entries exceed the range [0,"<<col_dim-1<<"]"<<std::endl;
    return 1;
  }

  delete colrep;
  colrep=0;

  int err=0;
  Indexmatrix map_to_new(col_dim,1,-1);
  for(Integer i=0;i<map_to_old.dim();i++){
    Integer ind=map_to_old(i);
    if (map_to_new(ind)<0)
      map_to_new(ind)=i;
    else {
      if (cb_out())
	get_out()<<"**** ERROR: SparseCoeffmatMatrix::reassign_columns(.): map entry "<<i<<" refers again to "<<ind<<" (first occurence was map entry "<<map_to_new(ind)<<"), assigning zero instead"<<std::endl;
      err++;
    }
  }
  col_dim=map_to_old.dim();
  
  dense_cnt.init(dense_cnt.dim(),1,Integer(0));
  for(unsigned i=0;i<blockrep.size();i++){
    if (blockrep[i].size()>0){
      SparseCoeffmatVector tmp;
      for (SparseCoeffmatVector::iterator it=blockrep[i].begin();
	   it!=blockrep[i].end();
	   it++){
	Integer ind=map_to_new(it->first);
	if (ind>=0) {
	  tmp[ind]=it->second;
	  if (it->second->dense())
	    dense_cnt(Integer(i))++;
	}
      }
      swap(tmp,blockrep[i]);
    }
  }
  
  return err;
}

//******************************************************************************
//                  SparseCoeffmatMatrix::delete_blocks
//******************************************************************************
 
  int SparseCoeffmatMatrix::delete_blocks(const Indexmatrix& delete_indices, Indexmatrix* map_to_old)
{ 
  if (delete_indices.dim()==0){
    if (map_to_old)
      map_to_old->init(Range(0,block_dim.dim()-1));
    return 0;
  }

  if ((min(delete_indices)<0)||(max(delete_indices)>=block_dim.dim())){
    if (cb_out())
      get_out()<<"**** ERROR: SparseCoeffmatMatrix::delete_blocks(..): some index entries exceed the range [0,"<<block_dim.dim()-1<<"]"<<std::endl;
    return 1;
  }

  Indexmatrix mto;
  if (map_to_old==0)
    map_to_old=&mto;

  map_to_old->init(Range(0,block_dim.dim()-1));
  int err=0;
  for(Integer i=0;i<delete_indices.dim();i++){
    Integer ind=delete_indices(i);
    if ((*map_to_old)(ind)>=0)
      (*map_to_old)(ind)=-(i+1);
    else {
      err++;
      if (cb_out())
	get_out()<<"**** ERROR: SparseCoeffmatMatrix::delete_blocks(..): deletion index "<<i<<" asks for deleting block "<<ind<<" previously deleted due to index "<<-(*map_to_old)(ind)-1<<", ignoring this"<<std::endl;
    }
  }
  Integer cnt=0;
  for(Integer i=0;i<map_to_old->dim();i++){
    if ((*map_to_old)(i)>=0){
      (*map_to_old)(cnt++)=(*map_to_old)(i);
    }
  }
  map_to_old->reduce_length(cnt);
  err+=reassign_blocks(*map_to_old);

  return err;
}

//******************************************************************************
//                  SparseCoeffmatMatrix::delete_columns
//******************************************************************************
 
  int SparseCoeffmatMatrix::delete_columns(const Indexmatrix& delete_indices, Indexmatrix* map_to_old)
{ 
  if (delete_indices.dim()==0){
    if (map_to_old)
      map_to_old->init(Range(0,col_dim-1));
    return 0;
  }

  if ((min(delete_indices)<0)||(max(delete_indices)>=col_dim)){
    if (cb_out())
      get_out()<<"**** ERROR: SparseCoeffmatMatrix::delete_columns(..): some index entries exceed the range [0,"<<col_dim-1<<"]"<<std::endl;
    return 1;
  }

  Indexmatrix mto;
  if (map_to_old==0)
    map_to_old=&mto;

  map_to_old->init(Range(0,col_dim-1));
  int err=0;
  for(Integer i=0;i<delete_indices.dim();i++){
    Integer ind=delete_indices(i);
    if ((*map_to_old)(ind)>=0)
      (*map_to_old)(ind)=-(i+1);
    else {
      err++;
      if (cb_out())
	get_out()<<"**** ERROR: SparseCoeffmatMatrix::delete_columns(..): deletion index "<<i<<" asks for deleting block "<<ind<<" previously deleted due to index "<<-(*map_to_old)(ind)-1<<", ignoring this"<<std::endl;
    }
  }
  Integer cnt=0;
  for(Integer i=0;i<map_to_old->dim();i++){
    if ((*map_to_old)(i)>=0){
      (*map_to_old)(cnt++)=(*map_to_old)(i);
    }
  }
  map_to_old->reduce_length(cnt);
  err+=reassign_columns(*map_to_old);

  return err;
}

//******************************************************************************
//                  SparseCoeffmatMatrix::delete_columns
//******************************************************************************
 
bool SparseCoeffmatMatrix::operator==(const SparseCoeffmatMatrix& mat) const
{
  if ((block_dim.dim()!=mat.block_dim.dim())||(col_dim!=mat.col_dim))
    return false;
  for(Integer i=0;i<block_dim.dim();i++){
    if (block_dim(i)!=mat.block_dim(i))
      return false;
    if (dense_cnt(i)!=mat.dense_cnt(i))
      return false;
    if (blockrep[unsigned(i)].size()!=mat.blockrep[unsigned(i)].size())
      return false;
    for(Integer j=0;j<col_dim;j++)
      if ((*this)(i,j)!=mat(i,j))
	return false;
  }
	  
  return true;
}

//******************************************************************************
//                  SparseCoeffmatMatrix::gramip
//******************************************************************************

int SparseCoeffmatMatrix::Gram_ip(Matrix& ipvec,
				  const Matrix& P,
				  const Matrix* Lam,
				  const Indexmatrix* ind) const
{
  int err=0;
  if (sum(block_dim)!=P.rowdim()){
    err++;
    if (cb_out()){
      get_out()<<"**** ERROR: SparseCoeffmatMatrix::Gram_ip(...): order of the block diagonal matrix ="<<sum(block_dim)<<" is not the same as the number of rows of the Gram matrix ="<<P.rowdim()<<std::endl;
    }
  }
  if ((ind)&&(ind->dim()>0)&&((min(*ind)<0)||(max(*ind)>=col_dim))){
    err++;
    if (cb_out()){      
      get_out()<<"**** ERROR: SparseCoeffmatMatrix::Gram_ip(...): some indices exceed the column range [0,"<<col_dim-1<<"]"<<std::endl;
    }
  }
  if (err)
    return err;
  
  if (ind)
    ipvec.init(ind->dim(),1,0.); 
  else 
    ipvec.init(col_dim,1,0.);
  
  if (block_dim.rowdim()==1){
    const SparseCoeffmatVector& row=blockrep[0];
    if (row.size()>0){
      if (ind){
	for(Integer i=0;i<ipvec.rowdim();i++){
	  SparseCoeffmatVector::const_iterator it=row.find((*ind)(i));
	  if (it!=row.end()) {
	    ipvec(i)=it->second->gramip(P,0,Lam);
	  }
	}
      }
      else {
	for(SparseCoeffmatVector::const_iterator it=row.begin();it!=row.end();it++){
	  ipvec(it->first)=it->second->gramip(P,0,Lam);
	}
      }
    }
  }
  else { 
    Integer start=0;
    for(Integer j=0;j<block_dim.rowdim();j++){
      const SparseCoeffmatVector& row=blockrep[unsigned(j)];
      if (row.size()>0){
	if (ind){
	  for(Integer i=0;i<ipvec.rowdim();i++){
	    SparseCoeffmatVector::const_iterator it=row.find((*ind)(i));
	    if (it!=row.end()) {
	      ipvec(i)+=it->second->gramip(P,start,Lam);
	    }
	  }
	}
	else {
	  for(SparseCoeffmatVector::const_iterator it=row.begin();it!=row.end();it++){
	    ipvec(it->first)=it->second->gramip(P,start,Lam);
	  }
	}
      }
      start+=block_dim(j);
    }
  }
  
  return err;
}
  
//******************************************************************************
//                  SparseCoeffmatMatrix::Gram_ip
//******************************************************************************

int SparseCoeffmatMatrix::Gram_ip(Real& ipval,const Matrix& P,Integer j) const
{
  int err=0;
  if (sum(block_dim)!=P.rowdim()){
    err++;
    if (cb_out()){
      get_out()<<"**** ERROR: SparseCoeffmatMatrix::Gram_ip(...): order of the block diagonal matrix ="<<sum(block_dim)<<" is not the same as the number of rows of the Gram matrix ="<<P.rowdim()<<std::endl;
    }
  }
  if ((j<0)||(j>=col_dim)){
    err++;
    if (cb_out()){      
      get_out()<<"**** ERROR: SparseCoeffmatMatrix::Gram_ip(...): index ="<<j<<"exceeds the column range [0,"<<col_dim-1<<"]"<<std::endl;
    }
  }
  if (err)
    return err;

  ipval=0.;

  if (block_dim.rowdim()==1){
    SparseCoeffmatVector::const_iterator colj=blockrep[0].find(j);
    if (colj!=blockrep[0].end())
      ipval=colj->second->gramip(P);
  }
  else {
    if (colrep==0)
      form_colrep();
    SCMcolrep::const_iterator colit=colrep->find(j);
    if (colit!=colrep->end()){
      Integer start=0;
      Integer cnt=0;
      for(SparseCoeffmatVector::const_iterator colj=colit->second.begin();colj!=colit->second.end(); colj++){
	assert(cnt<=colj->first);
	while (cnt<colj->first){
	  start+=block_dim(cnt);
	  cnt++;
	}
	ipval+=colj->second->gramip(P,start);
      }
    }
  }
  
  return err;
}

//******************************************************************************
//                  SparseCoeffmatMatrix::primal_ip
//******************************************************************************

int SparseCoeffmatMatrix::primal_ip(Matrix& ipvec,const PSCPrimal* primal,const Indexmatrix* ind) const
{
  int err=0;
  if ((ind)&&(ind->dim()>0)&&((min(*ind)<0)||(max(*ind)>=col_dim))){
    err++;
    if (cb_out()){      
      get_out()<<"**** ERROR: SparseCoeffmatMatrix::primal_ip(...): some indices exceed the column range [0,"<<col_dim-1<<"]"<<std::endl;
    }
  }
  if (err)
    return err;

  if (ind)
    ipvec.init(ind->dim(),1,0.); 
  else 
    ipvec.init(col_dim,1,0.);

  if (primal==0){ //is ok if this has only zero entries corresponding to ind
    if (ind) {
      for(Integer i=0;i<ipvec.rowdim();i++){
	if (column((*ind)(i))!=0){ 
	  err++;
	  break;
	}
      }
    }
    else {
      for (unsigned int i=0;i<blockrep.size();i++){
	if (blockrep[i].size()>0){
	  err++;
	  break;
	}
      }
    }
  }
  else {  //primal!=0
    for(Integer i=0;i<ipvec.rowdim();i++){
      Integer indi=i;
      if (ind)
	indi=(*ind)(i);
      if (primal->primal_ip(ipvec(i),*this,indi)){
	err++;
	if (cb_out()){
	  get_out()<<"**** ERROR in SparseCoeffmatMatrix::primal_ip(...): for multiple indices, primal->primal_ip(..) failed for index="<<indi<<std::endl;
	}
      }
    }
  }
  return err;
}

  
//******************************************************************************
//                  SparseCoeffmatMatrix::primal_ip
//******************************************************************************

int SparseCoeffmatMatrix::primal_ip(Real& ipval,const PSCPrimal* primal,Integer j) const
{
  int err=0;
  if ((j<0)||(j>=col_dim)){
    err++;
    if (cb_out()){      
      get_out()<<"**** ERROR: SparseCoeffmatMatrix::primal_ip(...): index ="<<j<<"exceeds the column range [0,"<<col_dim-1<<"]"<<std::endl;
    }
  }
  if (err)
    return err;

  if (primal==0){ //is ok if this has only zero entries corresponding to ind
    ipval=0.;
    if (column(j)!=0){ 
      err++;
    }
  }
  else {  //primal!=0
    if (primal->primal_ip(ipval,*this,j)){
      err++;
      if (cb_out()){
	get_out()<<"**** ERROR in SparseCoeffmatMatrix::primal_ip(...): failed for index="<<j<<std::endl;
      }
    }
  }

  return err;
}


//******************************************************************************
//                  SparseCoeffmatMatrix::project
//******************************************************************************

  
    
int SparseCoeffmatMatrix::project(Symmatrix& S,const Matrix& P,Integer j) const
{
  int err=0;
  if (sum(block_dim)!=P.rowdim()){
    err++;
    if (cb_out()){
      get_out()<<"**** ERROR: SparseCoeffmatMatrix::Gram_ip(...): order of the block diagonal matrix ="<<sum(block_dim)<<" is not the same as the number of rows of the Gram matrix ="<<P.rowdim()<<std::endl;
    }
  }
  if ((j<0)||(j>=col_dim)){
    err++;
    if (cb_out()){      
      get_out()<<"**** ERROR: SparseCoeffmatMatrix::Gram_ip(...): index ="<<j<<"exceeds the column range [0,"<<col_dim-1<<"]"<<std::endl;
    }
  }
  if (err)
    return err;

  if (block_dim.rowdim()==1){
    const SparseCoeffmatVector::const_iterator it=blockrep[0].find(j);
    if (it!=blockrep[0].end()){
      assert(it->second!=0);
      it->second->project(S,P);
    }
    else {
      S.init(P.coldim(),0.);
    }
  }
  else {
    const SparseCoeffmatVector* colj=column(j);
    if (colj==0){ //no nonzero matrix
      S.init(P.coldim(),0.);
    }
    else {
      SparseCoeffmatVector::const_iterator it=colj->begin();
      S.init(P.coldim(),0.);
      Integer start=0;
      Integer cnt=0;
      for (;it!=colj->end();it++){
	assert(cnt<=it->first);
	while(cnt<it->first){
	  start+=block_dim(cnt);
	  ++cnt;
	}
	assert(it->second!=0);
	it->second->add_projection(S,P,1.,start);
      }
    }
  }

  return err;
}



}

