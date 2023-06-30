/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/AFTData.cxx
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



#include <sstream>
#include <fstream>
#include "mymath.hxx"
#include "AFTData.hxx"

 
using namespace CH_Matrix_Classes;

namespace ConicBundle {

// *****************************************************************************
//                              ~AFTData
// *****************************************************************************

AFTData::~AFTData()
{}

// *****************************************************************************
//                                clear
// *****************************************************************************

void AFTData::clear(Integer start_modification_id)
{
  BundleData::clear(start_modification_id);

  aftpoint_id=-1;
  function_id=-1;

  center_ub_fid=-1;
  center_mid=-1;
  center_pid=-1;
  aft_center.init(0,0,0.);
  center_offset=0.;

  cand_ub_fid=-1;
  cand_mid=-1;
  cand_pid=-1;
  aft_cand.init(0,0,0.);
  cand_offset=0.;

  aggregate_fid=-1;
  aggregate_available=false;

}


// *****************************************************************************
//                                init
// *****************************************************************************

int AFTData::init(const BundleData* bd)
{
  const AFTData* aftd=dynamic_cast<const AFTData*>(bd);
  if (aftd==0){
    if (cb_out())
      get_out()<<"**** ERROR AFTData::init(): dynamic cast failed, argument is not of type const AFTData*"<<std::endl; 
    return 1;
  }
  clear();
  BundleData::init(bd);

  function_id=aftd->function_id;
  aftpoint_id=aftd->aftpoint_id;

  center_ub_fid=aftd->center_ub_fid;
  center_mid=aftd->center_mid;
  center_pid=aftd->center_pid;
  aft_center=aftd->aft_center;
  center_offset=aftd->center_offset;

  cand_ub_fid=aftd->cand_ub_fid;
  cand_mid=aftd->cand_mid;
  cand_pid=aftd->cand_pid;
  aft_cand=aftd->aft_cand;
  cand_offset=aftd->cand_offset;

  aggregate_fid=aftd->aggregate_fid;
  aggregate_available=aftd->aggregate_available;

  return 0;
}

// *****************************************************************************
//                                clone
// *****************************************************************************

BundleData* AFTData::clone() const
{
  AFTData* bd=new AFTData;
  bd->init(this);
  return bd;
}

// *****************************************************************************
//                              do_step
// *****************************************************************************


int AFTData::do_step(Integer point_id)
{
  if (cb_out(10)){
    get_out()<<"\n  entering  AFTData::do_step"<<std::endl;
  }

  int retval=BundleData::do_step(point_id);
  if (retval){
    if (cb_out())
      get_out()<<"\n**** ERROR AFTData::do_step(.): BundleData::do_step(.) returned "<<std::endl;
    if (cb_out(10)){
      get_out()<<"\n  leaving  AFTData::do_step with return value "<<retval<<std::endl;
    }
    return retval;
  }
  
  assert(cand_mid==modification_id);
  assert(cand_pid>=0);
  
  center_ub_fid=cand_ub_fid;  
  center_mid=cand_mid;  
  center_pid=cand_pid;
  aft_center=aft_cand;    // it would be better to use aft->copy_traforows here
  center_offset=cand_offset;

  if (cb_out(10)){
    get_out()<<"\n  leaving  AFTData::do_step with return value "<<retval<<std::endl;
  }

  return retval;
}

// *****************************************************************************
//                            synchronize_ids
// *****************************************************************************

  int AFTData::synchronize_ids(Integer& new_center_ub_fid, 
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

  function_id=max(function_id,center_ub_fid);
  function_id=max(function_id,cand_ub_fid);
  prex_id=new_prex_id;

  if ((new_center_ub_fid>=0)&&
      ((old_center_id!=center_id)||
       (center_ub_fid!=function_id)||
       (center_ub_mid!=modification_id)||
       (new_center_ub_fid!=modification_id)
       )){
    if (cb_out())
      get_out()<<"\n**** ERROR AFTData::synchronize_ids(): new_center_ub_fid="<<new_center_ub_fid<<" with old_center_id="<<old_center_id<<", but center_id="<<center_id<<", center_ub_mid="<<center_ub_mid<<", and modification_id="<<modification_id<<",  so center_ub is not up to date and return value for center_ub_fid will be -1"<<std::endl;
    new_center_ub_fid=-1;
  }
  if ((new_center_id>=0)&&
      ((old_center_id!=center_id)||
       (center_mid!=modification_id)
       )){
    if (cb_out())
      get_out()<<"\n**** ERROR AFTData::synchronize_ids(): new_center_id="<<new_center_id<<" with old_center_id="<<old_center_id<<", but center_id="<<center_id<<", center_ub_fid="<<center_ub_fid<<", function_id="<<function_id<<",  so center cannot by synchronized, discarding all center data"<<std::endl;
    new_center_id=-1;
  }
  if ((new_cand_ub_fid>=0)&&
      ((old_cand_id!=cand_id)||
       (cand_ub_fid!=function_id)||
       (cand_ub_mid!=modification_id)||
       (new_cand_ub_fid!=modification_id)
       )){
    if (cb_out())
      get_out()<<"\n**** ERROR AFTData::synchronize_ids(): new_cand_ub_fid="<<new_cand_ub_fid<<" with old_cand_id="<<old_cand_id<<", but cand_id="<<cand_id<<", cand_ub_fid="<<cand_ub_fid<<", function_id="<<function_id<<", cand_ub_mid="<<cand_ub_mid<<", and modification_id="<<modification_id<<",  so cand_ub is not up to date nad return value for cand_ub_fid will be -1"<<std::endl;
    new_cand_ub_fid=-1;
  }
  if ((new_cand_id>=0)&&
      ((old_cand_id!=cand_id)||
       (cand_mid!=modification_id)
       )){
    if (cb_out())
      get_out()<<"\n**** ERROR AFTData::synchronize_ids(): new_cand_id="<<new_cand_id<<" with old_cand_id="<<old_cand_id<<", but cand_id="<<cand_id<<",  so cand cannot by synchronized, discarding all candidate data"<<std::endl;
    new_cand_id=-1;
  }
  if ((new_aggregate_id>=0)&&(new_aggregate_id!=aggregate_id)){
    if (cb_out())
      get_out()<<"\n**** ERROR AFTData::synchronize_ids(): new_aggregate_id="<<new_aggregate_id<<", but aggreagte_id="<<aggregate_id<<",  so aggregate is not up to date and return value for aggregate_id will be -1"<<std::endl;
    new_aggregate_id=-1;
  }
  

  aftpoint_id=-1;
  //modificaton_id will be set to 0 afterwards
  //function_id will be synchronized aftwerwards but outside this routine 

  if (new_center_id>=0){
    center_id=new_center_id;
    center_pid=++aftpoint_id;
    center_mid=0;
    if ((center_ub_mid==modification_id)&&(center_ub_fid==function_id)){
      center_ub_mid=0;
    }
    else {
      center_ub_mid=center_ub_fid=-1;
    }
    if (new_center_ub_fid>=0){
      new_center_ub_fid=center_ub_mid;
    }
    else {
      new_center_ub_fid=-1;
    }
  }
  else {
    center_id=center_mid=center_pid=-1;
    new_center_ub_fid=center_ub_mid=center_ub_fid=-1;
  }
  
  if (new_cand_id>=0){
    cand_id=new_cand_id;
    cand_pid=++aftpoint_id;
    cand_mid=0;
    if ((cand_ub_mid==modification_id)&&(cand_ub_fid==function_id)){
      cand_ub_mid=0;
    }
    else {
      cand_ub_mid=cand_ub_fid=-1;
    }
    if (new_cand_ub_fid>=0){
      new_cand_ub_fid=cand_ub_mid;
    }
    else {
      new_cand_ub_fid=-1;
    }
  }
  else {
    cand_id=cand_mid=cand_pid=-1;
    new_cand_ub_fid=cand_ub_mid=cand_ub_fid=-1;
  }
  
  modification_id=0;
  
  aggregate_id=0;
  //aggregate_fid will be set in AFTModel::synchronize_ids
  if ((aggregate_available)&&(new_aggregate_id>=0)){
    if (aggregate.valid())
      aggregate.synchronize_ids(modification_id,
				new_center_id,old_center_id,
				new_cand_id,old_cand_id,new_prex_id);
    new_aggregate_id=aggregate_id;
  }
  else {
    new_aggregate_id=-1;
  }

  sumbundle.synchronize_ids(modification_id,
			new_center_id,old_center_id,
			new_cand_id,old_cand_id,
			new_prex_id);   
 
  //aggregate_fid will be set in AFTModel::synchronize_ids

  assert(old_minorants.size()==0);
    
  return 0;
}


// *****************************************************************************
//                             center_modified
// *****************************************************************************


  bool AFTData::center_modified(Integer& in_center_ub_mid,Integer check_center_id)
 {
  if (cb_out(10)){
    get_out()<<"\n  entering  AFTData::center_modified"<<std::endl;
  }
  assert(modification_id>=0);

  if ((check_center_id!=center_id)||
      (center_mid!=modification_id)||
      (center_ub_mid!=modification_id)||
      (center_ub_mid!=in_center_ub_mid)||
      (center_ub_fid!=function_id)) {
    if (cb_out(10)){
      get_out()<<"\n  leaving  AFTData::center_modified with return value true"<<std::endl;
    }
    in_center_ub_mid=modification_id;
    return true;
  }
  
  if (cb_out(10)){
    get_out()<<"\n  leaving  AFTData::center_modified with return value false"<<std::endl;
  }
  return false;
 }

// *****************************************************************************
//                           model_aggregate_modified
// *****************************************************************************


bool AFTData::model_aggregate_modified(Integer last_aggr_id)
{
  if (cb_out(10)){
    get_out()<<"\n  entering  AFTData::model_aggregate_modified"<<std::endl;
  }
  
  if ((!aggregate_available)||(last_aggr_id!=aggregate_id)) {
    if (cb_out(10)){
      get_out()<<"\n  leaving  AFTData::model_aggregate_modified with return value true"<<std::endl;
    }
    return true;
  }
  
  if (cb_out(10)){
    get_out()<<"\n  leaving  AFTData::model_aggregate_modified with return value false"<<std::endl;
  }
  return false;
}

// *****************************************************************************
//                                   clear_model
// *****************************************************************************

 void AFTData::clear_model(bool discard_minorants_only)
 {
   BundleData::clear_model(discard_minorants_only);

   //here all can be regenerated at negligible cost if possible at all, so rather do that
   center_ub_mid=-1;
   cand_ub_mid=-1;
   clear_aggregates();
 }

// *****************************************************************************
//                                   clear_aggregates
// *****************************************************************************

  void AFTData::clear_aggregates()
 {
   aggregate_fid=-1;
   aggregate_available=false;
   aggregate.clear();
 }


  

}

