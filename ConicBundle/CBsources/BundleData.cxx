/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/BundleData.cxx
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



#include "mymath.hxx"
#include "BundleData.hxx"

 
using namespace CH_Matrix_Classes;

namespace ConicBundle {

// *****************************************************************************
//                              BundleData
// *****************************************************************************

BundleData::BundleData()
{
  function_factor=1.;
  function_task=ObjectiveFunction;
  max_old_minorants=0;
  prex_id=0;
  clear();
}

// *****************************************************************************
//                              ~BundleData
// *****************************************************************************

BundleData::~BundleData()
{}

// *****************************************************************************
//                                clear
// *****************************************************************************

void BundleData::clear(CH_Matrix_Classes::Integer start_modification_id)
{
  dim=-1;

  modification_id=start_modification_id;


  center_ub_mid=-1;
  center_ub=CB_plus_infinity;
  center_id=-1;
  center_relprec=1e-3;

  cand_ub_mid=-1;
  cand_ub=CB_plus_infinity;
  cand_id=-1;
  cand_relprec=1e-3;

  local_aggregate.clear();
  aggregate.clear();
  aggregate_id=0;

  old_minorants.clear();
  minorant_nexti=0;
  set_lowrankH().init(0,0,0.);
  set_diagH().init(0,0,0.);
  set_denseH().init(0,0.);
  model_curvature=0.;

}

// *****************************************************************************
//                                init
// *****************************************************************************

int BundleData::init(const BundleData* bd)
{
  dim=bd->dim;

  function_factor=bd->function_factor;
  function_task=bd->function_task;

  modification_id=bd->modification_id;

  center_ub_mid=bd->center_ub_mid;
  center_id=bd->center_id;
  center_ub=bd->center_ub;
  center_relprec=bd->center_relprec;

  cand_ub_mid=bd->cand_ub_mid;
  cand_id=bd->cand_id;
  cand_ub=bd->cand_ub;
  cand_relprec=bd->cand_relprec;

  local_aggregate=bd->local_aggregate;
  aggregate=bd->aggregate;
  aggregate_id=bd->aggregate_id;

  max_old_minorants=bd->max_old_minorants;
  old_minorants=bd->old_minorants;
  set_lowrankH()=bd->get_lowrankH();
  set_diagH()=bd->get_diagH();
  set_denseH()=bd->get_denseH();
  minorant_nexti=bd->minorant_nexti;
  model_curvature=bd->model_curvature;

  sumbundle.init(bd->sumbundle);

  prex_id=bd->prex_id;


  return 0;
}

// *****************************************************************************
//                                clone
// *****************************************************************************

BundleData* BundleData::clone() const
{
  BundleData* bd=new BundleData;
  bd->init(this);
  return bd;
}

// *****************************************************************************
//                              do_step
// *****************************************************************************


int BundleData::do_step(Integer point_id)
{
  if (cb_out(10)){
    get_out()<<"\n  entering  BundleData::do_step"<<std::endl;
  }
  assert(modification_id>=0);

  if ((cand_ub_mid!=modification_id)||(point_id!=cand_id)) {
    if (cb_out())
      get_out()<<"\n**** ERROR BundleData::do_step(): unavailable or inconsistent candidate information (cand_ub_mid="<<cand_ub_mid<<", modificatin_id="<<modification_id<<", point_id="<<point_id<<", cand_id="<<cand_id<<"), return 1"<<std::endl;
    return 1;
  }

  center_ub_mid=cand_ub_mid;
  center_ub=cand_ub;
  center_relprec=cand_relprec;

  center_id=cand_id;

  if (cb_out(10)){
    get_out()<<"\n  leaving  BundleData::do_step"<<std::endl;
  }

  return 0;
}

// *****************************************************************************
//                          store_old_minorant
// *****************************************************************************

int BundleData::store_old_minorant(MinorantPointer& minorant)
{
  if (max_old_minorants<=0)
    return 0;
  
  //--- make sure sizes are ok with max_old_minorants requirements
  if (max_old_minorants<Integer(old_minorants.size())){
    //MinorantBundle tmp; tmp.resize(unsigned(max_old_minorants));
    MinorantBundle tmp((unsigned(max_old_minorants)));
    for(Integer i=1;i<=max_old_minorants;i++){
      if (--minorant_nexti<0)
	minorant_nexti=Integer(old_minorants.size())-1;
      tmp[unsigned(max_old_minorants-i)]=old_minorants[unsigned(minorant_nexti)];
    }
    swap(tmp,old_minorants);
    minorant_nexti=0;
  }
  else if ((max_old_minorants>Integer(old_minorants.size()))&&(minorant_nexti<Integer(old_minorants.size()))){
    MinorantBundle tmp(old_minorants.size());
    for(Integer i=0;i<=Integer(tmp.size());i++,minorant_nexti++){
      if (minorant_nexti==Integer(old_minorants.size()))
	minorant_nexti=0;
      tmp[unsigned(i)]=old_minorants[unsigned(minorant_nexti)];
    }
    swap(tmp,old_minorants);
    minorant_nexti=Integer(old_minorants.size());
  }

  //--- if the MinorantPointer is not empty, store it
  if (!minorant.empty()){
    minorant_nexti%=max_old_minorants;
    if (minorant_nexti==Integer(old_minorants.size())){
      old_minorants.push_back(minorant);
    }
    else {
      old_minorants[unsigned(minorant_nexti)]=minorant;
    }
    minorant_nexti++;
  }
  return 0;
}

// *****************************************************************************
//                          get_latest_minorants
// *****************************************************************************

int BundleData::get_latest_minorants(MinorantBundle& latest_minorants,
				     Integer max_number)
{
  max_old_minorants=max(max_old_minorants,max_number);
  if (Integer(old_minorants.size())<max_number)
    latest_minorants=old_minorants;
  else {
    latest_minorants.resize(unsigned(max_number));
    unsigned ind=unsigned(minorant_nexti);
    for(unsigned i=0;i<unsigned(max_number);i++){
      if (ind==0){
	ind=unsigned(old_minorants.size());
      }
      --ind;
      latest_minorants[i]=old_minorants[ind];
    }
  }
  return 0;
}

// *****************************************************************************
//                          get_model_data
// *****************************************************************************

int BundleData::get_model_data(MinorantBundle& model_minorants,
			       Matrix& modelcoeff) const
{
  model_minorants.clear();
  modelcoeff.init(0,0,0.);
  return 0;
}
  
// *****************************************************************************
//                            synchronize_ids
// *****************************************************************************

int BundleData::synchronize_ids(Integer& new_center_ub_fid, 
				Integer new_center_id,
				Integer old_center_id,
				Integer& new_cand_ub_fid, 
				Integer new_cand_id,
				Integer old_cand_id,
				Integer& new_aggregate_id,
			        Integer new_prex_id)
{
  assert(modification_id>=0);
  assert((new_center_ub_fid<0)||(new_center_id>=0));
  assert((new_center_id<0)||(old_center_id>=0));
  assert((new_cand_ub_fid<0)||(new_cand_id>=0));
  assert((new_cand_id<0)||(old_cand_id>=0));
  assert((new_center_id!=new_cand_id)||(new_center_id<0)||(new_cand_id<0)||(old_center_id==old_cand_id));
  prex_id=new_prex_id;

  if ((new_center_ub_fid>=0)&&
      ((old_center_id!=center_id)||
       (center_ub_mid!=modification_id)||
       (new_center_ub_fid!=modification_id)
       )){
    if (cb_out(1))
      get_out()<<"\n**** WARNING BundleData::synchronize_ids(): new_center_ub_fid="<<new_center_ub_fid<<" with old_center_id="<<old_center_id<<", but center_id="<<center_id<<", center_ub_mid="<<center_ub_mid<<", and modification_id="<<modification_id<<",  so center_ub is not up to date and return value for center_ub_fid will be -1"<<std::endl;
    new_center_ub_fid=-1;
  }
  if ((new_center_id>=0)&&(old_center_id!=center_id)){
    if (cb_out(1))
      get_out()<<"\n**** WARNING BundleData::synchronize_ids(): new_center_id="<<new_center_id<<" with old_center_id="<<old_center_id<<", but center_id="<<center_id<<",  so center cannot by synchronized, discarding all center data"<<std::endl;
    new_center_id=-1;
  }
  if ((new_cand_ub_fid>=0)&&
      ((old_cand_id!=cand_id)||
       (cand_ub_mid!=modification_id)||
       (new_cand_ub_fid!=modification_id)
       )){
    if (cb_out(1))
      get_out()<<"\n**** WARNING BundleData::synchronize_ids(): new_cand_ub_fid="<<new_cand_ub_fid<<" with old_cand_id="<<old_cand_id<<", but cand_id="<<cand_id<<", cand_ub_mid="<<cand_ub_mid<<", and modification_id="<<modification_id<<",  so cand_ub is not up to date and return value for cand_ub_fid will be -1"<<std::endl;
    new_cand_ub_fid=-1;
  }
  if ((new_cand_id>=0)&&(old_cand_id!=cand_id)){
    if (cb_out())
      get_out()<<"\n**** WARNING BundleData::synchronize_ids(): new_cand_id="<<new_cand_id<<" with old_cand_id="<<old_cand_id<<", but cand_id="<<cand_id<<",  so cand cannot by synchronized, discarding all candidate data"<<std::endl;
    new_cand_id=-1;
  }
  if ((new_aggregate_id>=0)&&(new_aggregate_id!=aggregate_id)){
    if (cb_out())
      get_out()<<"\n**** WARNING BundleData::synchronize_ids(): new_aggregate_id="<<new_aggregate_id<<", but aggreagte_id="<<aggregate_id<<",  so aggregate is not up to date and return value for aggregate_id will be -1"<<std::endl;
    new_aggregate_id=-1;
  }
    
    
  if (new_center_id>=0){
    center_id=new_center_id;
    if (center_ub_mid==modification_id){
      center_ub_mid=0;
    }
    else {
      center_ub_mid=-1;
    }
    if (new_center_ub_fid>=0){
      new_center_ub_fid=center_ub_mid;
    }
    else {
      new_center_ub_fid=-1;
    }
  }
  else {
    center_id=center_ub_mid=new_center_ub_fid=-1;
  }
  
  if (new_cand_id>=0){
    cand_id=new_cand_id;
    if (cand_ub_mid==modification_id){
      cand_ub_mid=0;
    }
    else {
      cand_ub_mid=-1;
    }
    if (new_cand_ub_fid>=0){
      new_cand_ub_fid=cand_ub_mid;
    }
    else {
      new_cand_ub_fid=-1;
    }
  }
  else {
    cand_id=cand_ub_mid=new_cand_ub_fid=-1;
  }
  
  modification_id=0;  
  aggregate_id=0;

 
  if ((aggregate.valid())&&(new_aggregate_id>=0)){
    if (!local_aggregate.empty())
      local_aggregate.synchronize_ids(modification_id,
				      new_center_id,old_center_id,
				      new_cand_id,old_cand_id,prex_id);
    aggregate.synchronize_ids(modification_id,
			      new_center_id,old_center_id,
			      new_cand_id,old_cand_id,prex_id);
    new_aggregate_id=aggregate_id;
  }
  else {
    new_aggregate_id=-1;
  }
  sumbundle.synchronize_ids(modification_id,
			new_center_id,old_center_id,
			new_cand_id,old_cand_id,
			prex_id);   

  for(unsigned i=0;i<old_minorants.size();i++){
    old_minorants[i].synchronize_ids(modification_id,
			new_center_id,old_center_id,
			new_cand_id,old_cand_id,
			prex_id);
  }

  return 0;
}
  
  
// *****************************************************************************
//                             center_modified
// *****************************************************************************


  bool BundleData::center_modified(Integer& check_center_ub_fid, Integer check_center_id)
 {
  if (cb_out(10)){
    get_out()<<"\n  entering  BundleData::center_modified"<<std::endl;
  }

  if ((center_ub_mid!=modification_id)||
      (check_center_ub_fid!=center_ub_mid)||
      (check_center_id!=center_id)) {
    if (cb_out(10)){
      get_out()<<"\n  leaving  BundleData::center_modified with return value true"<<std::endl;
    }
    check_center_ub_fid=modification_id;
    return true;
  }
  
  if (cb_out(10)){
    get_out()<<"\n  leaving  BundleData::center_modified with return value false"<<std::endl;
  }
  return false;
 }




// *****************************************************************************
//                           model_aggregate_modified
// *****************************************************************************


  bool BundleData::model_aggregate_modified(Integer last_aggr_id)
{
  if (cb_out(10)){
    get_out()<<"\n  entering  BundleData::model_aggregate_modified"<<std::endl;
  }
  
  if ((!aggregate.valid())||(last_aggr_id!=aggregate_id)) {
    if (cb_out(10)){
      get_out()<<"\n  leaving  BundleData::model_aggregate_modified with return value true"<<std::endl;
    }
    return true;
  }
  
  if (cb_out(10)){
    get_out()<<"\n  leaving  BundleData::model_aggregate_modified with return value false"<<std::endl;
  }
  return false;
}

// *****************************************************************************
//                               apply_modification
// *****************************************************************************

int BundleData::apply_modification(const GroundsetModification& gsmdf,MinorantExtender* mex)
{
  if (gsmdf.no_modification())
    return 0;

  int err=0;

  if ((gsmdf.appended_vardim()>0)&&(mex==0)) {

    aggregate.clear();
    local_aggregate.clear();
    old_minorants.clear();
    minorant_nexti=0;
    //sumbundle manages itself except for the lowest level, which is not dealt with here 

  }
  else {

    if ((aggregate.valid())&&(aggregate.apply_modification(gsmdf,modification_id,mex))){
      if (cb_out()) 
	get_out()<<"**** WARNING: BundleData::apply_modification(..): apply_modification failed for aggregate"<<std::endl;
      err++;
    }
    else {
      if ((aggregate.valid())&&(!gsmdf.no_modification()))
	aggregate_id++;
    }

    if ((local_aggregate.valid())&&(local_aggregate.apply_modification(gsmdf,modification_id,mex))){
      if (cb_out()) 
	get_out()<<"**** WARNING: BundleData::apply_modification(..): apply_modification failed for local_aggregate"<<std::endl;
      err++;
    }

    for(unsigned int i=0;i<old_minorants.size();i++){
      if (old_minorants[i].apply_modification(gsmdf,modification_id,mex)){
	if (cb_out())
	  get_out()<<"**** WARNING: BundleData::apply_modification(..): apply_modification failed for old_minorants number "<<i<<std::endl;
	err++;
      }
    }

    //sumbundle manages itself except for the lowest level, which is not dealt with here 
  }
  
  if (set_lowrankH().coldim()>0){
    if(set_lowrankH().rowdim()!=gsmdf.old_vardim()){
      if (cb_out()) 
	get_out()<<"**** WARNING: BundleData::apply_modification(..): set_lowrankH().rowdim()="<<set_lowrankH().rowdim()<<" while GroundsetModification assumes dim="<<gsmdf.old_vardim()<<std::endl;
      err++;
      set_lowrankH().init(0,0,0.);
    }
    else {
      set_lowrankH().enlarge_below(gsmdf.appended_vardim(),0.);
      if (gsmdf.map_to_old_variables())
	set_lowrankH()=set_lowrankH().rows(*gsmdf.map_to_old_variables());
    }
  }

  if (set_diagH().coldim()>0){
    if(set_diagH().rowdim()!=gsmdf.old_vardim()){
      if (cb_out()) 
	get_out()<<"**** WARNING: BundleData::apply_modification(..): set_diagH().rowdim()="<<set_diagH().rowdim()<<" while GroundsetModification assumes dim="<<gsmdf.old_vardim()<<std::endl;
      err++;
      set_diagH().init(0,0,0.);
    }
    else {
      set_diagH().enlarge_below(gsmdf.appended_vardim(),0.);
      if (gsmdf.map_to_old_variables())
	set_diagH()=set_diagH().rows(*gsmdf.map_to_old_variables());
    }
  }

  if (set_denseH().rowdim()>0){
    if(set_denseH().rowdim()!=gsmdf.old_vardim()){
      if (cb_out()) 
	get_out()<<"**** WARNING: BundleData::apply_modification(..): set_denseH().rowdim()="<<set_denseH().rowdim()<<" while GroundsetModification assumes dim="<<gsmdf.old_vardim()<<std::endl;
      err++;
      set_denseH().init(0,0.);
    }
    else {
      set_denseH().enlarge_below(gsmdf.appended_vardim(),0.);
      if (gsmdf.map_to_old_variables())
	set_denseH()=set_denseH().principal_submatrix(*gsmdf.map_to_old_variables());
    }
  }

  return err; 

}

// *****************************************************************************
//                                   clear_model
// *****************************************************************************

 void BundleData::clear_model(bool discard_minorants_only)
 {
   
   old_minorants.clear();
   minorant_nexti=0;
   local_aggregate.clear();
   aggregate.clear();

   set_lowrankH().init(0,0,0.);
   set_diagH().init(0,0,0.);
   set_denseH().init(0,0.);

   if (!discard_minorants_only)
     center_ub_mid=cand_ub_mid=-1;

 }

// *****************************************************************************
//                                   clear_aggregates
// *****************************************************************************

  void BundleData::clear_aggregates()
 {
   local_aggregate.clear();
   aggregate.clear();
   
   unsigned cnt=0;
   for(unsigned i=0;i<old_minorants.size();i++){
     if (! old_minorants[i].aggregate()){
       if (cnt<i) {
	 old_minorants[cnt]=old_minorants[i];
       }
       cnt++;
     }
   }
   if (minorant_nexti>Integer(cnt)){
     minorant_nexti=Integer(cnt);
   }
   old_minorants.resize(cnt);
   
 }

// *****************************************************************************
//                               call_primal_extender
// *****************************************************************************

  int BundleData::call_primal_extender(PrimalExtender& prex,bool /*include_candidates*/)
{
  prex_id++;
  int err=0;
  if (aggregate.call_primal_extender(prex,prex_id)){
    if (cb_out()) 
      get_out()<<"**** WARNING: NNCData::call_primal_extender(..): PrimalExtender::extend failed for aggregate"<<std::endl;
    err++;
  }
  else {
    if (aggregate.valid())
      aggregate_id++;
  }
 
  for(unsigned int i=0;i<old_minorants.size();i++){
    if (old_minorants[i].call_primal_extender(prex,prex_id)){
      if (cb_out())
	get_out()<<"**** WARNING: NNCData::call_primal_extender(..): PrimalExtender::extend failed for old_minorants number "<<i<<std::endl;
      err++;
    }
  }
 
  if (sumbundle.has_bundle_for(function_task)){
    if (sumbundle.call_primal_extender(prex,prex_id,function_task)){
      if (cb_out()) 
	get_out()<<"**** WARNING: NNCData::call_primal_extender(..):  sumbundle.call_primal_extender(..) failed"<<std::endl;
      err++;
    }
  }

  return err; 
}

  

}

