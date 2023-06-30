/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/ModificationBase.cxx
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



#include <typeinfo>
#include "ModificationBase.hxx"

using namespace CH_Matrix_Classes;

namespace ConicBundle {


ModificationBase::~ModificationBase()
{}



// *****************************************************************************
//                           ModificationBase::adapt_map_to_old
// *****************************************************************************

int ModificationBase::adapt_map_to_old(Indexmatrix*& map_to_old,
				   Indexmatrix*& del_ind,
				   Indexmatrix*& new_ind,
				   Indexmatrix& append_del_ind,
				   const Indexmatrix& input_map,
				   Integer& append_dim,
				   Integer olddim,
				   Integer newdim) const
{
  Indexmatrix inv_map(olddim+append_dim,1,-1);
  int err=0;
  Integer cnt_append=0;
  for (Integer i=0;i<input_map.dim();i++){
    Integer ind=input_map(i);
    if ((ind<0)||(ind>=newdim)){
      if (cb_out())
	get_out()<<"**** ERROR: ModificationBase::adapt_map_to_old(...): input_map("<<i<<")="<<ind<<" exceeds range ["<<0<<", "<<newdim-1<<"]"<<std::endl;
      err++;
      continue;
    }
    if(map_to_old)
      ind=(*map_to_old)(ind);
    if (inv_map(ind)>=0){
      if (cb_out())
	get_out()<<"**** ERROR: ModificationBase::adapt_map_to_old(...): input_map("<<inv_map(ind)<<")="<<input_map(i)<<"=input_map("<<i<<") assigned twice"<<std::endl;
      err++;
      continue;
    }
    inv_map(ind)=i;
    if (ind>=olddim){
      cnt_append++;
    }
  }
  if (err){
    return err;
  }
  
  if (map_to_old==0)
    map_to_old=new Indexmatrix(0,1,Integer(0));
  map_to_old->init(input_map.dim(),1,-1);
  if (cnt_append>0){
    if (new_ind==0)
      new_ind=new Indexmatrix;
    new_ind->init(cnt_append,1,Integer(0));
  }
  else { 
    delete new_ind;
    new_ind=0;
  } 
  append_del_ind.init(append_dim-cnt_append,1,0);
  if (del_ind==0)
    del_ind=new Indexmatrix;
  del_ind->init(olddim+cnt_append-input_map.dim(),1,0);


  Integer vdii=0;
  for(Integer i=0;i<olddim;i++){
    Integer ind=inv_map(i);
    if (ind<0){
      (*del_ind)(vdii)=i;
      vdii++;
    }
    else {
      assert((*map_to_old)(ind)==-1);
      (*map_to_old)(ind)=i;
    }
  }
  assert(vdii==del_ind->dim());
  Integer adii=0;
  for (Integer i=0;i<append_dim;i++){
    Integer ind=inv_map(i+olddim);
    if (ind<0){
      append_del_ind(adii)=i;
      adii++;
    }
    else {
      assert((*map_to_old)(ind)==-1);
      (*map_to_old)(ind)=i-adii+olddim;
      (*new_ind)(i-adii)=ind;
    }
  }
  assert(adii==append_del_ind.dim());
  assert(min(*map_to_old)>=0);
  assert(max(*map_to_old)<olddim+cnt_append);
  assert(map_to_old->coldim()==1);

  append_dim=cnt_append;

  return 0;
}

// *****************************************************************************
//                           ModificationBase::form_map_to_old
// *****************************************************************************

int ModificationBase::form_map_to_old(Indexmatrix& map_to_old,
		       const Indexmatrix& del_ind,
		       Integer dim) const
{
  map_to_old.init(Range(0,dim-1));
  if (del_ind.dim()==0){
    return 0;
  }
  int err=0;
  for(Integer i=0;i<del_ind.dim();i++){
    Integer ind=del_ind(i);
    if ((ind<0)||(ind>=dim)){
      if (cb_out())
	get_out()<<"**** ERROR: ModificationBase::form_map_to_old(...): delete_index("<<i<<")="<<ind<<" exceeds range [0,"<<dim-1<<"]"<<std::endl;
      err++;
      continue;
    }
    if (map_to_old(ind)<0){
      if (cb_out())
	get_out()<<"**** ERROR: ModificationBase::form_map_to_okd(...): delete_index("<<i<<")="<<ind<<"=delete_index("<<-map_to_old(ind)-1<<"), deleting the same variable more than once"<<std::endl;
      err++;
      continue;
    }
    map_to_old(ind)=-(i+1);
  }

  if (err)
    return err;

  Integer j=0;
  for(Integer i=0;i<dim;i++){
    Integer ind=map_to_old(i);
    if (ind>=0){
      map_to_old[j++]=ind;
    }
  }

  map_to_old.reduce_length(j);

  assert(map_to_old.coldim()==1);

  assert(j==dim-del_ind.dim());
  return 0;
}


}

