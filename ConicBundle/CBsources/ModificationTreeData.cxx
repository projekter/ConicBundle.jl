/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/ModificationTreeData.cxx
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




#include "ModificationTreeData.hxx"
#include "SumBlockModel.hxx"
#include "AFTModel.hxx"

using namespace CH_Matrix_Classes;

namespace ConicBundle {
  
ModificationTreeData::ModificationTreeData(const FunctionObject* fo,
					   FunctionObject* wr,
					   SumBlockModel* sbm,
					   Integer gs_dim,
					   Integer fixed_dim,
					   AffineFunctionTransformation* aft,
					   const CBout* cb):
  CBout(cb),oraclemod(0)
{  
  assert(fo);
  assert(sbm);
  
  if (aft!=0){
    if ((aft->from_dim()>=0)&&(aft->from_dim()!=gs_dim)){
      if (cb_out())
	get_out()<<"**** ERROR ModificationTreeData::ModificationTreeData(......): column dimension of affine function transformation ="<<aft->from_dim()<<" != "<<gs_dim<<" = input groundset dimension, dimensions do not match"<<std::endl;
    }
    if ((aft->to_dim()>=0)&&(fixed_dim>=0)&&(aft->to_dim()!=fixed_dim)){
      if (cb_out())
	get_out()<<"**** ERROR ModificationTreeData::ModificationTreeData(......): row dimension of affine function transformation to_dim="<<aft->to_dim()<<" != "<<fixed_dim<<" = fixed_dim of function, dimensions do not match"<<std::endl;
    }	
  }
  if (((aft==0)||(aft->get_arg_trafo()==0))&&(fixed_dim>=0)&&(gs_dim!=fixed_dim)){
    if (cb_out())
      get_out()<<"**** ERROR ModificationTreeData::ModificationTreeData(......): no or identiy transformation, yet input groundset dimension ="<<gs_dim<<" != "<<fixed_dim<<" = fixed_dim of function, dimensions do not match"<<std::endl;
  }	
  
  funobject=fo;
  wrapper=wr;
  model=sbm;
  parent=0;
  
  fixed_dimension=fixed_dim;
  
  aftmod.clear(gs_dim,gs_dim);
  if ((aft)&&(aft->to_dim()>=0))
    aftmod.clear(gs_dim,aft->to_dim());
}
  
  
ModificationTreeData::~ModificationTreeData()
{
  assert(modification_pending()==false);
  assert(children.size()==0);
  unlink_subtree();
  delete model;
  delete wrapper;
}


bool ModificationTreeData::modification_pending() const
{
  return (oraclemod!=0) || (!aftmod.no_modification());
}
  
bool ModificationTreeData::pending_oracle_modification(Integer& old_dim,
						       Integer& new_dim,
						       Integer& append_dim,
						       const Indexmatrix*& map_to_old,
						       const Indexmatrix*& deleted_indices,
						       const Indexmatrix*& new_indices,
						       const OracleModification*& oracle_modification) const
{
  bool modify= (((oraclemod)&&(!oraclemod->no_modification()))||
		(aftmod.appended_rowdim()>0)||
		(aftmod.map_to_old_rows())
		);
  old_dim=aftmod.old_rowdim();
  new_dim=aftmod.new_rowdim();
  append_dim=aftmod.appended_rowdim();
  map_to_old=aftmod.map_to_old_rows();
  deleted_indices=aftmod.deleted_row_indices();
  new_indices=aftmod.new_row_indices();
  oracle_modification=oraclemod;
  return modify;
}
    
int ModificationTreeData::add_child(ModificationTreeData* fmd)
{
  assert(fmd!=this);
  assert(modification_pending()==false);
  int err =0;
  if (fmd->aftmod.old_vardim()!=aftmod.old_rowdim()){
    if (cb_out())
      get_out()<<"**** ERROR ModificationTreeData::add_child(......): child from_dim="<<fmd->aftmod.old_vardim()<<" != "<<aftmod.old_rowdim()<<" = to_dim of this function, dimensions do not match"<<std::endl;
    err++;
  }		
  if ((err==0)&& model->add_model(fmd->model)){
    if (cb_out()){
      get_out()<<"**** ERROR ModificationTreeData::add_child(......): model->add_model() failed for child model"<<std::endl;
    }		
    err++;
  }
  
  if (err==0){
    children[fmd->funobject]=fmd;
    fmd->parent=this;
  }

  return err;
}
    

int ModificationTreeData::unlink_subtree()
{
  assert(modification_pending()==false);
  if (parent==0)
    return 0;
  
  int err=0;
  
  SumBlockModel* sbm=parent->model->remove_model(model);
  if (sbm==0){
    if (cb_out())
      get_out()<<"**** ERROR ModificationTreeData::unlink_me(): parent->model->remove_model(model) failed"<<std::endl;
    err++;
  }		
  
  FunctionMap::iterator it=parent->children.find(funobject);
  if (it==parent->children.end()){
    if (cb_out())
      get_out()<<"**** ERROR ModificationTreeData::unlink_me(): parent->model->remove_model(model) failed"<<std::endl;
    err++;
  }
  else
    parent->children.erase(it);
  
  parent=0;
  
  return err;
}
  
    
int ModificationTreeData::delete_descendants(FunctionMap& funmap)
{
  int err=0;
  for(FunctionMap::iterator it=children.begin();it!=children.end();it++){
    if (it->second->delete_descendants(funmap)){
      if (cb_out())
	get_out()<<"**** ERROR ModificationTreeData::delete_subtree(.): delete_subtree failed for a child"<<std::endl;
      err++;
    }
    it->second->parent=0;
    SumBlockModel* sbm=model->remove_model(it->second->model);
    if (sbm==0){
      if (cb_out())
	get_out()<<"**** ERROR ModificationTreeData::delete_descendatns(): model->remove_model(model) failed"<<std::endl;
      err++;
    }		
    FunctionMap::iterator delit=funmap.find(it->first);
    if (delit==funmap.end()){
      if (cb_out())
	get_out()<<"**** ERROR ModificationTreeData::delete_subtree(.): child not in input function map"<<std::endl;
      err++;
    }
    else {
      funmap.erase(delit);
    }
    delete it->second;
  }
  children.clear();
  return err;
}
    
  
int ModificationTreeData::add_parents_to_map(FunObjModMap& funmap)
{
  if (parent==0)
    return 0;
  FunObjModMap::const_iterator it=funmap.find(parent->funobject);
  if (it!=funmap.end())
    return 0;
  funmap[parent->funobject]=FunctionObjectModification();
  return parent->add_parents_to_map(funmap);
}
  
  
int ModificationTreeData::update_subtree_modification(Integer add_dim,const Indexmatrix* map_to_old,const FunObjModMap* modmap)
{
  int err=0;

  Integer child_add_dim=0;   ///number of new parameters added to this oracle and its children (on top of old ones) in consequence of the current modifications (may e.g. be fewer than in add_dim due to resassignments)
  const Indexmatrix* child_map_to_old=0;  ///the reassignment of the variables after the addition (on top of old ones)
  Indexmatrix tmp_map_to_old;

  const OracleModification* in_om=0;
  const AFTModification* in_aftm=0;
  bool in_map=false;
  if (modmap!=0){
    FunObjModMap::const_iterator it=modmap->find(funobject);
    if (it!=modmap->end()) {
      in_map=true;
      in_om = it->second.get_oracle_modification();
      in_aftm= it->second.get_aft_modification();
    }
  }
  if (in_om!=0) {
    if (oraclemod==0){
      oraclemod=in_om->new_initial_oraclemodification(aftmod.old_rowdim());
      if (aftmod.appended_rowdim()>0)
	oraclemod->add_append_variables(aftmod.appended_rowdim());
      if (aftmod.map_to_old_variables())
	oraclemod->add_reassign_variables(int(aftmod.map_to_old_variables()->dim()),(const int *)(aftmod.map_to_old_variables()->get_store()));
    }
    if (oraclemod->incorporate(*in_om)){
      if (cb_out())
	get_out()<<"**** ERROR ModificationTreeData::update_subtree_modification(...): incorporating oracle modification failed"<<std::endl;
      err++;      
    }
  }
  if (in_aftm==0){
    //default action
    if ((fixed_dimension>=0)||((modmap!=0)&&(!in_map))){
      //map new variables to zero
      if (add_dim>0){
	Sparsemat s(aftmod.new_rowdim(),add_dim);
	if (aftmod.add_append_vars(add_dim,&s,0)){
	  if (cb_out())
	    get_out()<<"**** ERROR ModificationTreeData::update_subtree_modification(...): default add_append_vars for avoiding changes failed for the AFT modification"<<std::endl;
	  err++;      
	}
	//no action needed for oraclemod
      }
      if (map_to_old){
	if (aftmod.add_reassign_vars(*map_to_old)){
	  if (cb_out())
	    get_out()<<"**** ERROR ModificationTreeData::update_subtree_modification(...): default add_reassign_vars for avoiding changes failed for the AFT modification"<<std::endl;
	  err++;      
	}
      }
      //no action needed for oraclemod
      if ((oraclemod)&&(oraclemod->get_new_vardim()!=aftmod.new_rowdim())){
	if (cb_out())
	  get_out()<<"**** ERROR ModificationTreeData::update_subtree_modification(...): after default actions on the AFT for avoiding argument changes the number of arguments of the oracle="<<oraclemod->get_new_vardim()<<" does not match that resulting from the AFT="<<aftmod.new_rowdim()<<std::endl;
	err++;      
      }
    }
    else {// in_aftm==0 with fixed_dimension<0 && ((modmap==0)||in_map)
      //append new variables
      if (add_dim>0){
	child_add_dim=add_dim;
	Sparsemat s(aftmod.new_rowdim(),add_dim);
	if (aftmod.add_append_vars(add_dim,&s,0)){
	  if (cb_out())
	    get_out()<<"**** ERROR ModificationTreeData::update_subtree_modification(...): add_append_vars failed for passing new variables on via the AFT modification"<<std::endl;
	  err++;      
	}
	Indexmatrix indi(Range(0,add_dim-1));
	Indexmatrix indj(Range(aftmod.new_vardim()-add_dim,aftmod.new_vardim()-1));
	Matrix val(add_dim,1,1.);
	s.init(add_dim,aftmod.new_vardim(),add_dim,indi,indj,val);
	if (aftmod.add_append_rows(add_dim,&s,0)){
	  if (cb_out())
	    get_out()<<"**** ERROR ModificationTreeData::update_subtree_modification(...): add_append_rows failed for passing new variables on via the AFT modification"<<std::endl;
	  err++;      
	}
	if ((oraclemod)&&(in_om==0)){
	  if (oraclemod->add_append_variables(add_dim)){
	    if (cb_out())
	      get_out()<<"**** ERROR ModificationTreeData::update_subtree_modification(...): add_append_variables failed for appending new variables via the oracle modification"<<std::endl;
	    err++;      
	  }
	}
      }
      if (map_to_old){
	const AffineFunctionTransformation* aft=model->get_aftmodel()?model->get_aftmodel()->get_aft():0;
	bool id_case=((aft?(aft->get_arg_trafo()==0):true)&&(aftmod.preserves_identity()));

	if (aftmod.add_reassign_vars(*map_to_old)){
	  if (cb_out())
	    get_out()<<"**** ERROR ModificationTreeData::update_subtree_modification(...): add_reassign_vars failed for reassigning variables via the AFT modification"<<std::endl;
	  err++;      
	}

	if (id_case){
	  //perform the same changes on the function
	  child_map_to_old=map_to_old;
	  if (aftmod.add_reassign_rows(*map_to_old)){
	    if (cb_out())
	      get_out()<<"**** ERROR ModificationTreeData::update_subtree_modification(...): add_reassign_rows failed for passing on reassignments via the AFT modification"<<std::endl;
	    err++;      
	  }
	  if ((oraclemod)&&(in_om==0)){
	    if (oraclemod->add_reassign_variables(map_to_old->dim(),map_to_old->get_store())){
	      if (cb_out())
		get_out()<<"**** ERROR ModificationTreeData::update_subtree_modification(...): add_reassign_variables failed for reassigning variables via the oracle modification"<<std::endl;
	      err++;      
	    }
	  }
	}
	else { //non id case, check for zero rows after deletion
	  //for this identify all deleted columns by setting those entries to 0 
	  //that are mapped by map_to_old_variables
	  Indexmatrix delvar;
	  const Indexmatrix* varmap=aftmod.map_to_old_variables();
	  if (varmap){
	    delvar.init(aftmod.old_vardim()+aftmod.appended_vardim(),1,1);
	    for(Integer i=0;i<varmap->dim();i++){
	      delvar((*varmap)(i))=0;
	    }
	  }
	  else {
	    delvar.init(aftmod.old_vardim()+aftmod.appended_vardim(),1,0);
	  }
          //test the rows that have to be tested
	  Indexmatrix delrow(0,0,Integer(0));
	  const Indexmatrix* rowmap=aftmod.map_to_old_rows();
	  const Matrix* aftoffset=aft?aft->get_arg_offset():0;
	  const Sparsemat* afttrafo=aft?aft->get_arg_trafo():0;
	  const Matrix* modoffset=aftmod.get_append_rhs();
	  const Sparsemat* modcoltrafo=aftmod.get_append_cols();
	  const Sparsemat* modrowtrafo=aftmod.get_append_rows();
	  Integer dim=rowmap?rowmap->dim():aftmod.old_rowdim()+aftmod.appended_rowdim();
	  for(Integer i=0;i<dim;i++){
	    Integer ind=rowmap?(*rowmap)(i):i;
	    if (ind<aftmod.old_rowdim()){
	      //check for offset ==0
	      if ((aftoffset)&&((*aftoffset)(ind)!=0.))
		continue;
	      //check for trafo part ==0
	      if(afttrafo){
		Integer startind;
		Integer nz=afttrafo->row_nonzeros(ind,&startind);
		if (nz>0)
		  nz+=startind;
		else
		  startind=0;
		while((startind<nz)&&
		      (delvar(afttrafo->get_rowindex()(startind))==1))
		  startind++;
		if (startind<nz)
		  continue;
	      }
	      //check for modcoltrafo part ==0
	      if ((modcoltrafo)&&(modcoltrafo->row_nonzeros(ind)>0))
		continue;
	      //zero map, delete this row
	    }
	    else {  //ind >= aftmod.old_rowdim(), newly append rows
	      Integer appendind=ind-aftmod.old_rowdim();
	      //check for offset ==0
	      if ((modoffset)&&((*modoffset)(appendind)!=0.))
		continue;
	      //check for modrowtrafo 
	      if ((modrowtrafo)&&(modrowtrafo->row_nonzeros(appendind)>0))
		continue;
	      //zero map, delete this row
	    }
	    delrow.concat_below(i);
	  }  //end for
	  if (delrow.dim()>0){
	    tmp_map_to_old.init(Range(0,aftmod.new_rowdim()-1));
	    tmp_map_to_old.delete_rows(delrow);
	    child_map_to_old=&tmp_map_to_old;
	    //perform the same changes on the function
	    if (aftmod.add_reassign_rows(*child_map_to_old)){
	      if (cb_out())
		get_out()<<"**** ERROR ModificationTreeData::update_subtree_modification(...): add_reassign_rows failed for passing on deletion reassignments via the AFT modification"<<std::endl;
	      err++;      
	    }
	    if ((oraclemod)&&(in_om==0)){
	      if (oraclemod->add_reassign_variables(child_map_to_old->dim(),child_map_to_old->get_store())){
		if (cb_out())
		  get_out()<<"**** ERROR ModificationTreeData::update_subtree_modification(...): add_reassign_variables failed for reassigning variables via the oracle modification after deletions"<<std::endl;
		err++;      
	      }
	    }
	  } //endif delrow.dim()>0
	} 
      }
      if ((oraclemod)&&(oraclemod->get_new_vardim()!=aftmod.new_rowdim())){
	if (cb_out())
	  get_out()<<"**** ERROR ModificationTreeData::update_subtree_modification(...): after default actions on the AFT with argument changes the number of arguments of the oracle="<<oraclemod->get_new_vardim()<<" does not match that resulting from the AFT="<<aftmod.new_rowdim()<<std::endl;
	err++;      
      }
 
    }
  }
  else { //in_aftm!=0
    if (aftmod.incorporate(*in_aftm)){
      if (cb_out()){
	get_out()<<"**** ERROR ModificationTreeData::append_variables(...): incorporating the given AFT modification failed"<<std::endl;
      }
      err++;
    }
    if (parent&&parent->aftmod.new_rowdim()!=aftmod.new_vardim()){
      if (cb_out()){
	get_out()<<"**** ERROR ModificationTreeData::append_variables(...): incorporating the given AFT modification resulted in a new incoming dimension = "<<aftmod.new_vardim()<<" != "<<parent->aftmod.new_rowdim()<<" the dimension of the input to this function"<<std::endl;
      }
      err++;
    }
    if ((oraclemod)&&(in_om==0)){
      if (in_aftm->appended_rowdim()>0){
	if (oraclemod->add_append_variables(in_aftm->appended_rowdim())){
	  if (cb_out())
	    get_out()<<"**** ERROR ModificationTreeData::update_subtree_modification(...): add_append_variables failed for incorporated appending of variables via the oracle modification"<<std::endl;
	  err++;      
	}
      }
      if (in_aftm->map_to_old_rows()){
	if (oraclemod->add_reassign_variables(in_aftm->map_to_old_rows()->dim(),in_aftm->map_to_old_rows()->get_store())){
	  if (cb_out())
	    get_out()<<"**** ERROR ModificationTreeData::update_subtree_modification(...): add_reassign_variables failed for incorporated reassigning of variables via the oracle modification"<<std::endl;
	  err++;      
	}
      }
    }
    child_add_dim=in_aftm->appended_rowdim();
    child_map_to_old=in_aftm->map_to_old_rows();

    if ((oraclemod)&&(oraclemod->get_new_vardim()!=aftmod.new_rowdim())){
      if (cb_out())
	get_out()<<"**** ERROR ModificationTreeData::update_subtree_modification(...): after user changes on the AFT the number of arguments of the oracle="<<oraclemod->get_new_vardim()<<" does not match that resulting from the AFT="<<aftmod.new_rowdim()<<std::endl;
      err++;      
    }
  }
  
  for(FunctionMap::iterator it=children.begin();it!=children.end();it++){
    if (it->second->update_subtree_modification(child_add_dim,child_map_to_old,modmap)){
      if (cb_out())
	get_out()<<"**** ERROR ModificationTreeData::update_subtree_modification(...): update_subtree_modification failed for a child"<<std::endl;
      err++;      
    }
  }

  return err;
}


int ModificationTreeData::collect_subtree_modification(FunObjModMap& modmap)
{
  int err=0;
  AFTModification* aftm=0;
  if (!aftmod.no_modification()){
    bool aft_stays_id=((model->get_aftmodel()==0)||(model->get_aftmodel()->get_aft()->get_arg_trafo()==0));
    aft_stays_id =  (aft_stays_id && aftmod.groundset_changes_suffice_for_identity());
    if (! aft_stays_id)
      aftm=&aftmod;
  }
  if (aftm||oraclemod){
    if (wrapper)
      modmap[wrapper]=FunctionObjectModification(oraclemod,aftm);
    else
      modmap[funobject]=FunctionObjectModification(oraclemod,aftm);
  }
  for(FunctionMap::iterator it=children.begin();it!=children.end();it++){
    if (it->second->collect_subtree_modification(modmap)){
      if (cb_out())
	get_out()<<"**** ERROR ModificationTreeData::collect_subtree_modification(.): collect_subtree_modification failed for a child"<<std::endl;
      err++;      
    }
  }
  return err;
}

int ModificationTreeData::subtree_modification_performed()
{
  aftmod.clear(aftmod.new_vardim(),aftmod.new_rowdim());
  delete oraclemod;
  oraclemod=0;
  
  int err=0;
  for(FunctionMap::iterator it=children.begin();it!=children.end();it++){
    if (it->second->subtree_modification_performed()){
      if (cb_out())
	get_out()<<"**** ERROR ModificationTreeData::subtree_modification_performed(): collect_subtree_modification failed for a child"<<std::endl;
      err++;      
    }
  }
  return err;
}

int ModificationTreeData::clear_subtree_modification()
{
  aftmod.clear(aftmod.old_vardim(),aftmod.old_rowdim());
  delete oraclemod;
  oraclemod=0;
  
  int err=0;
  for(FunctionMap::iterator it=children.begin();it!=children.end();it++){
    if (it->second->clear_subtree_modification()){
      if (cb_out())
	get_out()<<"**** ERROR ModificationTreeData::clear_subtree_modification: collect_subtree_modification failed for a child"<<std::endl;
      err++;      
    }
  }
  return err;
}

  void ModificationTreeData::set_cbout(const CBout* cb,int /* incr */)
{
  CBout::set_cbout(cb,0);
  if (model)
    model->set_cbout(this,0);
  aftmod.set_cbout(this,0);
}


} //end namespace ConicBundle
