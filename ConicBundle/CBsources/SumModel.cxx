/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/SumModel.cxx
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
#include "SumModelParameters.hxx"
#include "BundleIdProx.hxx"

using namespace CH_Matrix_Classes;

namespace ConicBundle {

// *****************************************************************************
//                              SumModel::ModelData()
// *****************************************************************************

  SumModel::ModelData::ModelData(SumBlockModel* om):
    oracle_model(om),function_id(-1),center_ub_fid(-1),cand_ub_fid(-1),aggregate_id(-1),model_deviation(0.),deviation_contribution(0.),model_curvature(0.),weighted_curvature(0.),suggested_mode(SumBundle::inactive)
{
  assert(om);
}

// *****************************************************************************
//                              SumModel::~ModelData()
// *****************************************************************************

  SumModel::ModelData::~ModelData()
{
}

// *****************************************************************************
//                              SumModel()
// *****************************************************************************

  SumModel::SumModel(CBout* cb):SumBlockModel(cb),ncalls(0),block(0),model_selection(0)
{
  clear();
}

// *****************************************************************************
//                              ~SumModel()
// *****************************************************************************

SumModel::~SumModel()
{
  clear();
  delete model_selection;
}

// *****************************************************************************
//                              clear()
// *****************************************************************************

void SumModel::clear()
{
  if (bundlehandler){
    if (bundlehandler->remove_contributions()){
      if (cb_out()) 
	get_out()<<"**** ERROR: SumModel::clear(...): bundlehandler->remove_contribution failed"<<std::endl;
    }   
    for (ModelMap::iterator it=modelmap.begin();it!=modelmap.end();it++){
      SumBundle::Mode unavail=SumBundle::unavailable;
      it->second->model()->sumbundle_mode(unavail);
    }
  }
  data.clear();
  
  SumBlockModel::clear();
  
  for (ModelMap::iterator it=modelmap.begin();it!=modelmap.end();it++){
    delete it->second;
    it->second=0;
  }
  
  modelmap.clear();  
  delete block;
  block=0;
  ncalls=0;

  delete model_selection;
  model_selection=new  SumModelParameters(this);
  
  cand_minorant.clear();
  center_minorant.clear();
}
 
// *****************************************************************************
//                            model
// *****************************************************************************
 
const SumBlockModel* SumModel::model(const FunctionObject* fo) const
{
  if (fo==0) 
    return 0;
  ModelMap::const_iterator it=modelmap.find(fo);
  if (it==modelmap.end())
    return 0;
  return it->second->oracle_model;
}
 
// *****************************************************************************
//                            add_model
// *****************************************************************************

int SumModel::add_model(SumBlockModel* model)
{
  if (cb_out(10)){
    get_out()<<"\n  entering  SumModel::add_model"<<std::endl;
  }    
  if (model==0){
    if (cb_out()) get_out()<<"**** ERROR: SumModel::add_model(...): model pointer is NULL"<<std::endl;
    if (cb_out(10)){
      get_out()<<"\n  leaving  SumModel::add_model"<<std::endl;
    }    
    return 1;
  }
    
  const FunctionObject* fo=model->get_oracle_object();
  const ModelMap::iterator it=modelmap.find(fo);
  if (it!=modelmap.end()) {
    if (cb_out()) get_out()<<"**** ERROR: SumModel::add_model(...): model pointer holds a function object that is added already; function object must be unique"<<std::endl;
    if (cb_out(10)){
      get_out()<<"\n  leaving  SumModel::add_model"<<std::endl;
    }    
    return 1;
  }
  modelmap[fo]=new ModelData(model);
  model->set_cbout(this,0);   //output at the same level for member problems

  clear_model();
  data.modification_id++;

  if (cb_out(10)){
    get_out()<<"\n  leaving  SumModel::add_model"<<std::endl;
  }    
  return 0;
}

// *****************************************************************************
//                            remove_model
// *****************************************************************************

SumBlockModel* SumModel::remove_model(const FunctionObject* fo)
{
  if (cb_out(10)){
    get_out()<<"\n  entering  SumModel::remove_model"<<std::endl;
  }    
  if (fo==0){
    if (cb_out()) get_out()<<"**** WARNING: SumModel::remove_model(...): argument pointer is NULL"<<std::endl;
    if (cb_out(10)){
      get_out()<<"\n  leaving  SumModel::remove_model"<<std::endl;
    }    
   return 0;
  }
    
  const ModelMap::iterator it=modelmap.find(fo);
  if (it==modelmap.end()) {
    if (cb_out()) get_out()<<"**** WARNING: SumModel::remove_model(...): no model for the input function object"<<std::endl;
    if (cb_out(10)){
      get_out()<<"\n  leaving  SumModel::remove_model"<<std::endl;
    }    
    return 0;
  }

  SumBundle::Mode unavail=SumBundle::unavailable;
  if (it->second->model()->sumbundle_mode(unavail)){
      if (cb_out()) 
	get_out()<<"**** ERROR: SumModel::remove_model(...): sumbundle_mode() failed"<<std::endl;
      if (cb_out(10)){
	get_out()<<"\n  leaving  SumModel::remove_model"<<std::endl;
      }    
  }

  SumBlockModel* mp = it->second->oracle_model;

  delete it->second;
  modelmap.erase(it);

  clear_model();
  data.modification_id++;
  
  if (cb_out(10)){
    get_out()<<"\n  leaving  SumModel::remove_model"<<std::endl;
  }    
  return mp;
}
  
// *****************************************************************************
//                            eval_function
// *****************************************************************************

int SumModel::eval_function(Integer& ub_fid,
			    Real& ub,
			    Integer y_id,
			    const Matrix& y,
			    Real nullstep_bound,
			    Real relprec)
{
  if (cb_out(10)){
    get_out()<<"\n  entering  SumModel::eval_function"<<std::endl;
  }    
  CH_Tools::Microseconds start_eval=clock.time();

  if (data.dim<0)
    data.dim=y.dim();
  assert(data.dim==y.dim());
  assert(relprec>0.);
  assert(data.function_factor==1.);

  ncalls++;

  data.cand_ub_mid=-1;
  data.cand_id=y_id;
  data.cand_relprec=relprec;

  cand_minorant.clear();

  if (bundlehandler)
    bundlehandler->clear_cand_minorants();
  
  //compute a quick lower bound on each subfunction
  Matrix fun_lb(Integer(modelmap.size()),1,CB_minus_infinity);
  Integer i=0;
  for (ModelMap::iterator it=modelmap.begin();it!=modelmap.end();it++,i++){
    fun_lb(i)=it->second->model()->lb_function(y_id,y);
  }
  Real sum_lb=sum(fun_lb);
  
  //compute values for submodels
  int sumretval=0;
  data.cand_ub=0.;
  bool fid_increased=false;
  i=0;
  preeval_time+=clock.time()-start_eval;
  for (ModelMap::iterator it=modelmap.begin();it!=modelmap.end();it++,i++){
    sum_lb-=fun_lb(i);
    Real local_nullstep_bound=nullstep_bound-data.cand_ub-sum_lb;

    Real ub_val;
    start_eval=clock.time();
    int retval=it->second->model()->eval_function(it->second->cand_ub_fid,ub_val,y_id,y,
						 local_nullstep_bound,
						 relprec);
    
    eval_time+=clock.time()-start_eval;
    if (retval<0){
      if (cb_out()){
	get_out()<<"**** WARNING SumModel::eval_function: eval_function returned"<<retval<<" for function "<<"i"<<std::endl;
      }
      sumretval+=retval;
    }    
    if (retval>0){
      if (cb_out()){
	get_out()<<"**** ERROR SumModel::eval_function: eval_function failed for functionc "<<i<<" and returned"<<retval<<std::endl;
      }
      eval_time+=clock.time()-start_eval;
      if (cb_out(10)){
	get_out()<<"\n  leaving  SumModel::eval_function with return value "<<retval<<std::endl;
      }    
      return retval;
    }    
    if (it->second->cand_ub_fid>it->second->function_id){
      fid_increased=true;
      it->second->function_id=it->second->cand_ub_fid;
    }
    data.cand_ub+=ub_val;
    
  }
  start_eval=clock.time();
  
  //all information collected
  if (fid_increased){
    data.modification_id++;
    center_minorant.clear();
    cand_minorant.clear();
  }
  ub_fid=data.cand_ub_mid=data.modification_id;
  ub=data.cand_ub;

  //update contribution to the parent sumbundle if necessary
  if (data.sumbundle.has_contributions()){
    assert(bundlehandler);
    if (bundlehandler->contribute_new_minorants()){
	if (cb_out())
	  get_out()<<"**** ERROR SumModel::eval_function: bundlehandler->contribute_new_minorants failed"<<std::endl;
    }
  }

  //if needed collect information for dynamic scaling
   if (data.max_old_minorants>0){
    MinorantPointer tmp;
    if (get_function_minorant(tmp)){
      if (cb_out())
	get_out()<<"**** WARNING SumModel::eval_function: get_function_minorant(.) failed for adding olld minorant"<<std::endl;
    }
    else if (data.store_old_minorant(tmp)){
      if (cb_out())
	get_out()<<"**** WARNING SumModel::eval_function(): data.store_old_minorant() failed"<<std::endl;
    }
  }

  eval_time+=clock.time()-start_eval;
 
  if (cb_out(10)){
    get_out()<<"\n  leaving  SumModel::eval_function"<<std::endl;
  }    

  posteval_time+=clock.time()-start_eval;
  return sumretval;
}

// *****************************************************************************
//                            eval_model
// *****************************************************************************

int SumModel::eval_model(Real& lb,
			 Integer y_id,
			 const Matrix& y,
			 Real relprec)
{
  if (cb_out(10)){
    get_out()<<"\n  entering  SumModel::eval_model"<<std::endl;
  }    

  assert(relprec>0.);

  CH_Tools::Microseconds start_eval=clock.time();

  lb=0.;
 
  if ((bundlehandler)&&(data.sumbundle.active())){
    if (data.sumbundle.has_working_roots()){
      if (bundlehandler->eval_model(lb,y_id,y)){
	if (cb_out()){
	  get_out()<<"**** ERROR SumModel::eval_model(.....): bundlehandler->eval_model() failed"<<std::endl;
	}
      }
    }

    int sumretval=0;
    Integer i=0;
    for(std::list<const ModelData*>::iterator it=local_augmodel.begin();
	it!=local_augmodel.end();it++,i++){
      
      Real lb_val;
      int retval=(*it)->model()->eval_model(lb_val,y_id,y,
						 relprec);
      if (retval<0){
	if (cb_out()){
	  get_out()<<"**** WARNING SumModel::eval_model: eval_model returned"<<retval<<" for local_augmodel "<<"i"<<std::endl;
	}
	sumretval+=retval;
      }    
      if (retval>0){
	if (cb_out()){
	  get_out()<<"**** ERROR SumModel::eval_model: eval_model failed for model "<<i<<" and returned"<<retval<<std::endl;
	}
	eval_time+=clock.time()-start_eval;
	if (cb_out(10)){
	  get_out()<<"\n  leaving  SumModel::eval_model with return value "<<retval<<std::endl;
	}    
	return retval;
      }    
      lb+=lb_val;
    }
    
    if (cb_out(10)){
      get_out()<<"\n  leaving  SumModel::eval_model with retunr value "<<sumretval<<std::endl;
    }    
    
    return sumretval;
    
  }

  
  //compute values for submodels
  int sumretval=0;
  Integer i=0;
  for (ModelMap::iterator it=modelmap.begin();it!=modelmap.end();it++,i++){

    Real lb_val;
    int retval=it->second->model()->eval_model(lb_val,y_id,y,
					       relprec);
    if (retval<0){
      if (cb_out()){
	get_out()<<"**** WARNING SumModel::eval_model: eval_model returned"<<retval<<" for model "<<"i"<<std::endl;
      }
      sumretval+=retval;
    }    
    if (retval>0){
      if (cb_out()){
	get_out()<<"**** ERROR SumModel::eval_model: eval_model failed for model "<<i<<" and returned"<<retval<<std::endl;
      }
      eval_time+=clock.time()-start_eval;
      if (cb_out(10)){
	get_out()<<"\n  leaving  SumModel::eval_model with return value "<<retval<<std::endl;
      }    
      return retval;
    }    
    lb+=lb_val;
  }

  eval_time+=clock.time()-start_eval;
 
  if (cb_out(10)){
    get_out()<<"\n  leaving  SumModel::eval_model with return value "<<sumretval<<std::endl;
  }    

  return sumretval;
}


// *****************************************************************************
//                                update_model
// *****************************************************************************

  int SumModel::update_model(ModelUpdate model_update,
			     Integer center_id,
			     const Matrix& center_y,
			     Integer y_id,
			     const Matrix& y,
			     Real model_maxviol,
			     BundleProxObject& H)
{
  if (bundlehandler!=0){
    Real model_deviation=0;
    Real model_curvature=0;
    return update_model(model_update,center_id,center_y,y_id,y,model_maxviol,H,model_deviation,model_curvature);
  }

  if (cb_out(10)){
    get_out()<<"\n  entering  SumModel::update_model"<<std::endl;
  }


  assert((y_id==data.center_id)||(y_id==data.cand_id));

  if (((data.cand_id!=y_id)||(data.cand_ub_mid!=data.modification_id))&&
      ((data.center_id!=y_id)||(data.center_ub_mid!=data.modification_id))
      ) {
    if (cb_out(0)){
      get_out()<<"\n**** WARNING: SumModel::update_model(...) called without any available subgradient information"<<std::endl;
    }
    if (cb_out(10)){
      get_out()<<"\n  leaving  SumModel::update_model with return value 1"<<std::endl;
    }
    return 1;
  }


  CH_Tools::Microseconds update_start=clock.time();

  Real factor=Real(modelmap.size())*data.function_factor;
  Integer i=0;
  H.apply_factor(1./data.function_factor);
  for (ModelMap::iterator it=modelmap.begin();it!=modelmap.end();it++,i++){
    int retval=it->second->model()->update_model(model_update,center_id,center_y,y_id,y,model_maxviol/factor,H);
    if (retval){
      if (cb_out()) 
	get_out()<<"*** ERROR: SumModel::update_model: calling update_model for submodel "<<i<<" returned "<<retval<<std::endl;
      if (cb_out(10)){
	get_out()<<"\n leaving SumModel::update_model: return value 1"<<std::endl;
      }    
      H.apply_factor(data.function_factor);
      return 1;
    }  
    it->second->model_deviation=0.;
    it->second->model_curvature=0.;
    if (model_update==descent_step){
      it->second->center_ub_fid=it->second->cand_ub_fid;
    }
  }
  H.apply_factor(data.function_factor);

  if (model_update==descent_step){
    center_minorant.clear();
    data.do_step(y_id);
  }

  updatemodel_time+=clock.time()-update_start;

  if (cb_out(10)){
    get_out()<<"\n  leaving  SumModel::update_model with return value 0"<<std::endl;
  }
  return 0;
}
  
// *****************************************************************************
//                                update_model
// *****************************************************************************

int SumModel::update_model(ModelUpdate model_update,
			   Integer center_id,
			   const Matrix& center_y,
			   Integer y_id,
			   const Matrix& y,
			   Real model_maxviol,
			   BundleProxObject& H,
			   Real& model_deviation,
			   Real& model_curvature)
{
  if (cb_out(10)){
    get_out()<<"\n  entering  SumModel::update_model"<<std::endl;
  }
  assert((!data.sumbundle.has_contributions())|| bundlehandler);
  assert(!(data.sumbundle.has_roots()&&(bundlehandler==0)));
  assert((y_id==data.center_id)||(y_id==data.cand_id));

  if (((data.cand_id!=y_id)||(data.cand_ub_mid!=data.modification_id))&&
      ((data.center_id!=y_id)||(data.center_ub_mid!=data.modification_id))
      ) {
    if (cb_out(0)){
      get_out()<<"\n**** WARNING: SumModel::update_model(...) called without any available subgradient information"<<std::endl;
    }
    if (cb_out(10)){
      get_out()<<"\n  leaving  SumModel::update_model with return value 1"<<std::endl;
    }
    return 1;
  }

  CH_Tools::Microseconds update_start=clock.time();

  if ((!data.aggregate.valid())||((y_id!=data.center_id)&&(y_id!=data.cand_id))){
    model_deviation=0.;
  }
  else {
    if (data.cand_id==y_id){
      model_deviation=data.cand_ub;
    }
    else {
      model_deviation=data.center_ub;
    }
    model_deviation-=data.aggregate.evaluate(y_id,y);
  }


  //depending on the mode of the sumbundle and on global bh update the sumbundle

  Real factor=Real(modelmap.size())*data.function_factor;
  if (bundlehandler){
    if (bundlehandler->update_model(model_update,center_id,center_y,y_id,y,(Real(modelmap.size())-Real(local_augmodel.size()))*model_maxviol/factor,H)){
      if (cb_out(0)){
	get_out()<<"\n**** WARNING: SumModel::update_model(...): bundlehandler->update_model(.) failed"<<std::endl;
      }
      if (cb_out(10)){
	get_out()<<"\n  leaving  SumModel::update_model with return value 1"<<std::endl;
      }
      return 1;
    }
  }

  model_curvature=0;
  Integer i=0;
  H.apply_factor(1./data.function_factor);
  for (ModelMap::iterator it=modelmap.begin();it!=modelmap.end();it++,i++){
    Real md=0.;
    Real cv=0.;
    int retval=0;
    retval=it->second->model()->update_model(model_update,center_id,center_y,y_id,y,model_maxviol/factor,H,md,cv);
    assert(cv>=0.);
    if (retval){
      if (cb_out()) 
	get_out()<<"*** ERROR: SumModel::update_model: calling update_model for submodel "<<i<<" returned "<<retval<<std::endl;
      if (cb_out(10)){
	get_out()<<"\n leaving SumModel::update_model: return value 1"<<std::endl;
      }    
      H.apply_factor(data.function_factor);
      return 1;
    }  
    it->second->model_deviation=md;
    it->second->model_curvature=cv;
    if (cb_out(2)){
      get_out()<<"["<<i<<","<<md<<","<<cv<<"]"<<std::endl;
    }
    model_curvature+=cv;
    if (model_update==descent_step){
      it->second->center_ub_fid=it->second->cand_ub_fid;
    }
  }
  H.apply_factor(data.function_factor);

  if (model_update==descent_step){
    center_minorant.clear();
    data.do_step(y_id);
  }

  //--- if use_local_model and the new subgradient is already part of the model,
  //    it should be deleted again,  this is not yet implemented but should be done here

  updatemodel_time+=clock.time()-update_start;

  if (cb_out(10)){
    get_out()<<"\n  leaving  SumModel::update_model with return value 0"<<std::endl;
  }
  return 0;
}
  
// *****************************************************************************
//                             synchronize_ids
// *****************************************************************************


int SumModel::synchronize_ids(Integer& new_center_ub_fid, 
			      Integer new_center_id,
			      Integer old_center_id,
			      Integer& new_cand_ub_fid, 
			      Integer new_cand_id,
			      Integer old_cand_id,
			      Integer& new_aggregate_id)
{
  if (cb_out(10)){
    get_out()<<"\n  entering  SumModel::synchronize_ids"<<std::endl;
  }

  int retval=data.synchronize_ids(new_center_ub_fid,
				  new_center_id,
				  old_center_id,
				  new_cand_ub_fid,
				  new_cand_id,
				  old_cand_id,
				  new_aggregate_id,0);
  if (retval){
    if (cb_out()) 
      get_out()<<"*** ERROR: SumModel::synchronize_ids: calling data.synchronize_ids returned "<<retval<<std::endl;
    if (cb_out(10)){
      get_out()<<"\n  leaving  SumModel::synchronize_ids with return value 1"<<std::endl;
    }
    return 1;
  }

  bool center_changed=false;
  bool cand_changed=false;
  bool aggr_changed=false;

  Integer i=0;
  for (ModelMap::iterator it=modelmap.begin();it!=modelmap.end();it++,i++){
    ModelData& mdata=*(it->second);
    if ((data.center_id<0)||(mdata.center_ub_fid!=mdata.function_id))
      mdata.center_ub_fid=-1;
    if ((data.cand_id<0)||(mdata.cand_ub_fid!=mdata.function_id))
      mdata.cand_ub_fid=-1;
    int retval=mdata.model()->synchronize_ids(mdata.center_ub_fid,
					      data.center_id,
					      old_center_id,
					      mdata.cand_ub_fid,
					      data.cand_id,
					      old_cand_id,
					      mdata.aggregate_id);
    if (retval){
      if (cb_out()) 
	get_out()<<"*** ERROR: SumModel::synchronize_ids: calling synchronize_ids for submodel "<<i<<" returned "<<retval<<std::endl;
      if (cb_out(10)){
	get_out()<<"\n leaving SumModel::synchronize_ids: return value 1"<<std::endl;
      }    
      return 1;
    }
    mdata.function_id=max(-1,max(mdata.center_ub_fid,mdata.cand_ub_fid));
    if ((mdata.function_id<0)||(mdata.center_ub_fid!=mdata.function_id)){
      center_changed=true;
    }
    if ((mdata.function_id<0)||(mdata.cand_ub_fid!=mdata.function_id)){
      cand_changed=true;
    }
    if (mdata.aggregate_id<0)
      aggr_changed=true;
  }

  if (center_changed){
    data.center_ub_mid=-1;
    new_center_ub_fid=-1;
    center_minorant.clear();
  }
  if (cand_changed){
    data.cand_ub_mid=-1;
    new_cand_ub_fid=-1;
    cand_minorant.clear();
  }
  if (aggr_changed){
    data.local_aggregate.clear();
    data.aggregate.clear();
    new_aggregate_id=-1;
  }

  if (cb_out(10)){
    get_out()<<"\n  leaving  SumModel::synchronize_ids with return value 0"<<std::endl;
  }
  return 0;
}


// *****************************************************************************
//                             center_modified
// *****************************************************************************


  bool SumModel::center_modified(Integer& function_id,Integer center_id)
 {
  if (cb_out(10)){
    get_out()<<"\n  entering  SumModel::center_modified"<<std::endl;
  }

  bool modified=false;
  bool function_changed=false;
  Integer i=0;
  for (ModelMap::iterator it=modelmap.begin();it!=modelmap.end();it++,i++){
    //all have to be called recursively so that each can announce the need to recompute to its parent
    ModelData& mdata=*(it->second);
    Integer fid=mdata.center_ub_fid;
    if (mdata.model()->center_modified(fid,center_id)){
      modified=true;
      mdata.center_ub_fid=-1;
    }
    else {
      mdata.center_ub_fid=fid;
    }
    if (fid!=mdata.function_id){
      assert(fid>mdata.function_id);
      mdata.function_id=fid;
      function_changed=true;
    }
  }
  if (function_changed){
    data.modification_id++;
  }

  function_id=data.modification_id;

  if ((modified)||
      ((data.center_ub_mid!=data.modification_id)||(center_id!=data.center_id))
      ){
    data.center_ub_mid=-1;
    center_minorant.clear();
    cand_minorant.clear();
    if (cb_out(10)){
      get_out()<<"\n  leaving  SumModel::center_modified with return value true"<<std::endl;
    }
    return true;
  }
  
  if (cb_out(10)){
    get_out()<<"\n  leaving  SumModel::center_modified with return value false"<<std::endl;
  }
  return false;
 }


// *****************************************************************************
//                             recompute_center
// *****************************************************************************


  int SumModel::recompute_center(Integer& center_ub_fid,
				 Real& new_ub,
			       Integer center_id,
			       const Matrix& center_y,
			       bool accept_only_higher,
			       Real relprec)
 {
  if (cb_out(10)){
    get_out()<<"\n  entering  SumModel::recompute_center"<<std::endl;
  }
  if (data.dim<0)
    data.dim=center_y.dim();
  assert(data.dim==center_y.dim());

  if ((relprec<=0)&&((center_id!=data.center_id)||(data.center_relprec<=0))){
    Real rp=1e-3;
    if ((data.center_relprec>0)&&(rp>data.center_relprec))
      rp=data.center_relprec;
    if (cb_out()){
      get_out()<<"\n**** WARNING SumModel::recompute_center: called with relprec="<<relprec<<"<=0 but no previous precision available for center_id="<<center_id<<", using default choice "<<rp<<std::endl;
    }
    relprec=rp;
  }
  if (relprec<=0)
    relprec=data.center_relprec;

  if ((accept_only_higher)&&
      ((data.center_id!=center_id)||(data.center_ub_mid!=data.modification_id))
      ){
    if (cb_out()){
      get_out()<<"\n**** WARNING SumModel::recompute_center: called with accept_only_higher=true but there is no valid previous value or data.center_id="<<data.center_id<<" differs from center_id="<<center_id<<"; this should never happen, setting accept_only_higher to false"<<std::endl;
    }
    accept_only_higher=false;
  }
  
  if ((data.center_id==center_id)&&
      (data.center_ub_mid==data.modification_id)&&
      (center_ub_fid==data.modification_id)&&
      (relprec>=data.center_relprec)
      ){
    new_ub=data.center_ub;
    center_ub_fid=data.modification_id;
    return 0;
  }

  CH_Tools::Microseconds start_eval=clock.time();
  ncalls++;
  
  new_ub=0.;
  int sumretval=0;
  bool fid_increased=false;
  Integer i=0;
  for (ModelMap::iterator it=modelmap.begin();it!=modelmap.end();it++,i++){
    ModelData& mdata=*(it->second);
    Real ub;
    int retval=mdata.model()->recompute_center(mdata.center_ub_fid,ub,center_id,center_y,accept_only_higher,relprec);
    if (retval<0){
      if (cb_out(0)){
	get_out()<<"\n**** ERROR: SumModel::recompute_center(...): recompute_center failed to obtain the desired precision for submodel i="<<i<<" and returned "<<retval<<std::endl;
      }
      sumretval+=retval;
    }
    if (retval>0){
      if (cb_out(0)){
	get_out()<<"\n**** ERROR: SumModel::recompute_center(...): recompute_center failed for submodel i="<<i<<" and returned "<<retval<<std::endl;
      }
      if (cb_out(10)){
	get_out()<<"\n  leaving  SumModel::recompute_center with return value "<<retval<<std::endl;
      }
      eval_time+=clock.time()-start_eval;
      return retval;
    }
    if (mdata.center_ub_fid>mdata.function_id){
      fid_increased=true;
      mdata.function_id=mdata.center_ub_fid;
    }
    
    new_ub+=ub;
  }
  if (fid_increased){
    data.modification_id++;
    accept_only_higher=false;
    center_minorant.clear();
    cand_minorant.clear();
  }
  if ((!accept_only_higher)||(new_ub>data.center_ub)){
    data.center_ub=new_ub;
    data.center_relprec=relprec;
    center_minorant.clear();
  }
  data.center_id=center_id;
  data.center_ub_mid=data.modification_id;
  center_ub_fid=data.modification_id;

  eval_time+=clock.time()-start_eval;

  new_ub=data.center_ub;

  if (cb_out(10)){
    get_out()<<"\n  leaving  SumModel::recompute_center with return value "<<sumretval<<std::endl;
  }

  return sumretval;
 }



// *****************************************************************************
//                             model_aggregate_modified
// *****************************************************************************


  bool SumModel::model_aggregate_modified(Integer last_aggr_id)
 {
  if (cb_out(10)){
    get_out()<<"\n  entering  SumModel::model_aggregate_modified"<<std::endl;
  }

  for (ModelMap::iterator it=modelmap.begin();it!=modelmap.end();it++){
    if (it->second->model()->model_aggregate_modified(it->second->aggregate_id)){
      data.local_aggregate.clear();
      data.aggregate.clear();
    }    
  }
  
  if ((!data.aggregate.valid())||(last_aggr_id!=data.aggregate_id)) {
    if (cb_out(10)){
      get_out()<<"\n  leaving  SumModel::model_aggreagte_modified with return value true"<<std::endl;
    }
    return true;
  }
  
  if (cb_out(10)){
    get_out()<<"\n  leaving  SumModel::model_aggreagte_modified with return value false"<<std::endl;
  }
  return false;
 }


// *****************************************************************************
//                             provide_model_aggregate
// *****************************************************************************


  int SumModel::provide_model_aggregate(Integer y_id,const Matrix& y)
 {
  if (cb_out(10)){
    get_out()<<"\n  entering  SumModel::provide_model_aggregate"<<std::endl;
  }

  if (data.aggregate.valid()) {
    if (cb_out(10)){
      get_out()<<"\n  leaving  SumModel::provide_model_aggreagte with return value 0"<<std::endl;
    }
    return 0;
  }

  //--- first ensure that a bundle exists, initializing it if possible
  if ((bundlehandler)&&(bundlehandler->initialization_needed())){
    BundleIdProx H;
    int retval=update_model(new_subgradient,y_id,y,y_id,y,max(1e-6,1e-3*std::fabs(data.cand_ub)),H);
    if (retval){
      if (cb_out())
	get_out()<<"**** ERROR SumModel::provide_model_aggregate: update_model failed and returned "<<retval<<std::endl;
      return 1;
    }
  }

  data.local_aggregate.clear();
  data.aggregate.clear();

  //collect the local part
  Integer i=0;
  for (ModelMap::iterator it=modelmap.begin();it!=modelmap.end();it++,i++){
    int retval=it->second->model()->provide_model_aggregate(y_id,y);
    if (retval){
      if (cb_out(0)){
	get_out()<<"\n**** ERROR:  SumModel::provide_model_aggregate(...): provide_model_aggregate failed for submodel i="<<i<<" and returned "<<retval<<std::endl;
      }
      if (cb_out(10)){
	get_out()<<"\n  leaving  SumModel::provide_model_aggregate with return value 1"<<std::endl;
      }
      return 1;
    }
    retval=it->second->model()->get_model_aggregate(it->second->aggregate_id,data.local_aggregate,false);
    if (retval){
      if (cb_out(0)){
	get_out()<<"\n**** ERROR:  SumModel::provide_model_aggregate(...): get_model_aggregate failed for submodel i="<<i<<" and returned "<<retval<<std::endl;
      }
      if (cb_out(10)){
	get_out()<<"\n  leaving  SumModel::provide_model_aggregate with return value 1"<<std::endl;
      }
      
      return 1;
    }
  }
  if ((bundlehandler)&&(data.sumbundle.active())){
    if (bundlehandler->provide_model_aggregate()){
      if (cb_out())
	get_out()<<"**** ERROR SumModel::provide_model_aggregate: bundlehandler->ürovide_model_aggregate() failed"<<std::endl;
      return 1;
    }
  }
  if (data.sumbundle.get_local_model_aggregate(data.local_aggregate)){
    if (cb_out(0)){
      get_out()<<"\n**** ERROR:  SumModel::provide_model_aggregate(...): sumbundle.get_local_model_aggregate failed"<<std::endl;
    }
    if (cb_out(10)){
      get_out()<<"\n  leaving  SumModel::provide_model_aggregate with return value 1"<<std::endl;
    }  
    return 1;
  }

  //add the contributed part
  data.aggregate=data.local_aggregate;
  if (data.sumbundle.get_contributed_model_aggregate(data.aggregate)){
    if (cb_out(0)){
      get_out()<<"\n**** ERROR:  SumModel::provide_model_aggregate(...): sumbundle.get_local_model_aggregate failed"<<std::endl;
    }
    if (cb_out(10)){
      get_out()<<"\n  leaving  SumModel::provide_model_aggregate with return value 1"<<std::endl;
    }  
    return 1;
  }

  data.aggregate_id++;
   
  if (cb_out(10)){
    get_out()<<"\n  leaving  SumModel::provide_model_aggregate with return value 0"<<std::endl;
  }

  return 0;
 }


// *****************************************************************************
//                         apply_modification
// *****************************************************************************

int SumModel::apply_modification(bool& no_changes,
				 const GroundsetModification& gsmdf,
				 const FunObjModMap& funmdfmap,
				 Integer new_center_id,
				 const Matrix& new_center,
				 Integer old_center_id,
				 const Matrix& old_center)
{
  if (call_aftmodel_first(funmdfmap)){
    return sbm_transform()->apply_modification(no_changes,
					       gsmdf, funmdfmap,
					       new_center_id,new_center,
					       old_center_id,old_center);
  }

  assert ((data.dim<0)||(old_center_id<0)||(data.dim==old_center.dim()));
  assert((data.dim<0)||(data.dim==gsmdf.old_vardim()));

  // ------ if the contribution changes, remove it first
  int err=0;
  if ((!no_changes)||(!gsmdf.no_modification())){
    if (bundlehandler){
      if (bundlehandler->remove_contributions()){
	if (cb_out()){
	  get_out()<<"**** ERROR FunctionModel::apply_modification(.......): sumbundle->remove_contributions failed"<<std::endl;
	}
	err++;
      }
    }
  }

  no_changes=true;    // will be set to false if there are changes

  Integer i=0;
  bool and_no_changes=true;
  for (ModelMap::iterator it=modelmap.begin();it!=modelmap.end();it++,i++){
    bool local_no_changes=no_changes;
    int retval=it->second->model()->apply_modification(local_no_changes,
						       gsmdf,funmdfmap,
						       new_center_id,new_center,
						       old_center_id,old_center);
    if (retval){
      err++;
      if (cb_out(0)){
	get_out()<<"\n**** WARNING:  SumModel::apply_modification(..): apply_modification(..) failed for submodel i="<<i<<" and returned "<<retval<<std::endl;
      }
    }
    and_no_changes = and_no_changes && local_no_changes;
  }

  if ((!gsmdf.no_modification())&& and_no_changes){
    if (data.apply_modification(gsmdf,0)){
      if (cb_out(0))
	get_out()<<"\n**** WARNING:  SumModel::apply_modification(..): apply_modification(..) failed for data"<<std::endl;
      err++;
      no_changes=false;
    }
    else {
      local_augmodel.clear();
      cand_minorant.clear();
      center_minorant.clear();
      data.modification_id++;
    }
  }

  if ((!no_changes) || (!and_no_changes)){
    clear_model();
    data.modification_id++;
    no_changes=false;
  }
  
  data.dim=gsmdf.new_vardim();

  if (cb_out(10)){
    get_out()<<"\n   leaving SumModel::apply_modification with return value "<<err<<std::endl;
  }
  return err;  
}

// *****************************************************************************
//                              add_variable_metric
// *****************************************************************************


int SumModel::add_variable_metric(VariableMetric& H,
				  Integer y_id,
				  const Matrix& y,
				  bool descent_step,
				  Real weightu,
				  Real model_maxviol,
				  const Indexmatrix* indices)
{
  if ((!H.employ_variable_metric())||(!H.get_use_local_metric())){
    if (SumBlockModel::add_variable_metric(H,y_id,y,descent_step,weightu,
					   model_maxviol,indices)){
      if (cb_out())
	get_out()<<"**** WARNING: FunctionModel::add_variable_metric(): SumBlockModel::add_variable_metric() failed"<<std::endl;
      return 1;
    }
    return 0;
  }

  if (cb_out(10)){
    get_out()<<"\n  entering  SumModel::add_variable_metric"<<std::endl;
  }

  int err=0;
  CH_Tools::Microseconds start_metric=clock.time();


  Real factor=Real(modelmap.size())*data.function_factor;
  if ((bundlehandler)&&(data.sumbundle.has_working_roots())){
    if (bundlehandler->add_variable_metric(H,y_id,y,descent_step,
					   weightu,(Real(modelmap.size())-Real(local_augmodel.size()))*model_maxviol/factor,indices)){
      err++;
      if (cb_out())
	get_out()<<"**** WARNING: SumModel::add_variable_metric(...) bundlehandler->add_variable_metric(...) failed"<<std::endl;
    }
  }
  
  if (err==0) {
    Integer i=0; 
    for (ModelMap::iterator it=modelmap.begin();it!=modelmap.end();it++,i++){
      int retval=it->second->model()->add_variable_metric(H,y_id,y,
							  descent_step,
							  weightu/data.function_factor,
							  model_maxviol/factor,
							  indices);
      if (retval){
	if (cb_out(0)){
	  get_out()<<"\n**** WARNING:  SumModel::add_variable_metric(...): add_variable_metric failed for submodel i="<<i<<" and returned "<<retval<<std::endl;
	}
	err++;
	break;
      }
    }
  }

  metric_time+=clock.time()-start_metric;


  if (cb_out(10)){
    get_out()<<"\n leaving  SumModel::add_variable_metric with return value "<<err<<std::endl;
  }
  return err;
}

// *****************************************************************************
//                         check_center_validity_by_candidate
// *****************************************************************************


  int SumModel::check_center_validity_by_candidate(bool& cand_subg_is_below,
						   Integer center_id,
						   const Matrix& center_y)
{
  if (cb_out(10)){
    get_out()<<"\n  entering  SumModel::check_center_validity_by_candidate"<<std::endl;
  }    
  
  if ((data.center_id!=center_id)||
      (data.center_ub_mid!=data.modification_id)||
      (data.cand_ub_mid!=data.modification_id)
      ){
    if (cb_out()) 
      get_out()<<"**** WARNING: FunctionModel::check_center_validity_by_candidate: no center or no candidate minorant available for checking validity"<<std::endl;
    return 1;
  }  

  cand_subg_is_below=true;
  Integer i=0; 
  for (ModelMap::iterator it=modelmap.begin();it!=modelmap.end();it++,i++){
    bool is_below=true;
    int status=it->second->model()->check_center_validity_by_candidate(is_below,center_id,center_y);
    if (status){
      if (cb_out()) 
	get_out()<<"*** ERROR: SumModel::check_center_validity_by_candidate: calling check_center_validity_by_candidate for submodel "<<i<<" returned "<<status<<std::endl;
      if (cb_out(10)){
	get_out()<<"\n leaving SumModel::check_center_validity_by_candidate: return value 1"<<std::endl;
      }    
      return 1;
    }  
    if ((is_below==false)&&(cb_out(0))){
      get_out()<<"\n SumModel::check_center_validity_by_candidate: center_ub lies below cand_minorant for submodel "<<i<<std::endl;
    }
    cand_subg_is_below &= is_below;
  }

  if (cb_out(10)){
    get_out()<<"\n leaving SumModel::check_center_validity_by_candidate: cand_subg_is_below="<<cand_subg_is_below<<", return value 0"<<std::endl;
  } 
   
  return 0;
}
  
// *****************************************************************************
//                        sumbundle_mode
// *****************************************************************************

int SumModel::sumbundle_mode(SumBundle::Mode& mode,
			      SumBundleHandler* bh,
			      AffineFunctionTransformation* aft)
{
  assert((mode!=SumBundle::child)||(bh!=0));

  SumBundle::Mode next_mode=sumbundle_parameters->get_acceptable_mode();
  if ((next_mode==SumBundle::root)||
      ((next_mode==SumBundle::inactive)&&(mode==SumBundle::root))
      ){
    next_mode=SumBundle::root;
  }
  else if((next_mode==SumBundle::unavailable)||
	  ((bh==0)&&((mode!=SumBundle::root)||(next_mode==SumBundle::child)))
	  ){
    next_mode=SumBundle::unavailable;
  }
  else {   //note: bh!=0 (all cases==0 are covered)
    assert(bh);
    //next_mode is either inactive or child
    if ((next_mode==SumBundle::child)&&(mode==SumBundle::root)) {
      next_mode=SumBundle::inactive;
    }
    else if ((mode!=SumBundle::unavailable)&&(!bh->handles(get_data()->function_task))){
      next_mode=SumBundle::unavailable;
    }
    else { //mode = inactive or child or unavailable
      next_mode=mode;
    }
  }
  
  mode=next_mode;

  SumBundle::Mode current_mode=SumBundle::unavailable;
  if (get_data()->sumbundle.has_bundle_for(get_data()->function_task)){
    current_mode=get_data()->sumbundle.get_mode(get_data()->function_task);
  }
  assert(!(((current_mode==SumBundle::root)||(current_mode==SumBundle::child))&&(get_data()->get_local_model_active()))); //note: initially both might be inactive


  int err=0;
  // normalize sumbundle before adding or removing anything here or later
  if ((bundlehandler)&&(bundlehandler->normalize_sumbundle())){
    if (cb_out()){
      get_out()<<"*** ERROR: SumModel::sumbundle_mode(): bundlehandler->normalize_sumbundle() failed"<<std::endl;
    }
    err++;
  }

  //--- if no longer a child, remove contributions
  if ((current_mode==SumBundle::child)&&(mode!=SumBundle::child)){
    assert(bundlehandler);
    assert(bundlehandler->handles(get_data()->function_task));
    bundlehandler->remove_contributions();
  }
  
  //carry out the decision:
  //--- if mode==SumBundle::unavailable or no parent handler, remove an existing bundlehandler;
  //    if mode==SumBundle::inactive with parent handler, set it to inactive but keep following 
  if ((mode==SumBundle::unavailable)||(mode==SumBundle::inactive)){
    get_data()->local_model_active=true;
    if (mode==SumBundle::unavailable){
      bh=0;
    }
    if (bundlehandler!=0){
      if (bundlehandler->set_parent_information(bh,aft,SumBundle::inactive)){
	if (cb_out()){
	  get_out()<<"*** ERROR: ConeModel::sumbundle_mode(): SumBundleHandler::set_parent_information() failed"<<std::endl;
	}
	err++;
      }
      if (bh==0){
	bundlehandler->clear_model();
	delete bundlehandler;
	bundlehandler=0;
	mode=SumBundle::unavailable;
      }
    }
  }
  else {
    //initialize the local handler if necessary
    if (bundlehandler==0){
      std::map<FunctionTask,Real> factor_map;
      if (bh==0) {
	factor_map[ObjectiveFunction]=1.;
	factor_map[ConstantPenaltyFunction]=1.;
	factor_map[AdaptivePenaltyFunction]=1.;
      }
      else {
	if (bh->handles(ObjectiveFunction))
	  factor_map[ObjectiveFunction]=1.;
	if (bh->handles(ConstantPenaltyFunction))
	  factor_map[ConstantPenaltyFunction]=1.;
	if (bh->handles(AdaptivePenaltyFunction))
	  factor_map[AdaptivePenaltyFunction]=1.;      
      }
      bundlehandler=new SumBundleHandler(data.sumbundle,factor_map,sumbundle_parameters);
      bundlehandler->set_cbout(this,0);
      
      if (bh==0){
	//no bundle sizes are prespecified by a parent,
	//do a pseudo-update so that there is at least one virtual entry in each bundle that the others can contribute their aggregate to to start with
	Matrix dummy(0,0,0.);
	BundleIdProx H;
	H.set_weightu(1.);
	if (bundlehandler->update_model(BundleModel::null_step,0,dummy,0,dummy,1.,H)){
	  if (cb_out()){
	    get_out()<<"*** ERROR: SumModel::sumbundle_mode(): the pseudo call to bundlehandler->update_model() failed"<<std::endl;
	  }
	  err++;
	}
      }
    }
    
    if (bundlehandler->set_parent_information(bh,aft,mode)){
      if (cb_out()){
	get_out()<<"*** ERROR: SumModel::sumbundle_mode(): SumBundleHandler::set_parent_information() failed"<<std::endl;
      }
      err++;
    }
    
    //--- if it is a new root, set the local sumbundle_parameters
    if (mode==SumBundle::root){
      if (bundlehandler->set_bundle_parameters(*sumbundle_parameters)){
	if (cb_out()){
	  get_out()<<"*** ERROR: ConeModel::sumbundle_mode(): set_bundle_parameters() failed"<<std::endl;
	}
	err++;
      }
    }
  }
  
  return err;
}
  

// *****************************************************************************
//                                start_sumaugmodel
// *****************************************************************************

int SumModel::start_sumaugmodel(QPModelDataPointer& qpblockp,
				Integer cand_id,
				const Matrix& cand_y,
				const Indexmatrix* indices,
				SumBundleHandler* bh,
				SumBundle::Mode mode,
				AffineFunctionTransformation* aft)
{
  if (cb_out(10)){
    get_out()<<"\n  entering  SumModel::start_sumaugmodel"<<std::endl;
  }

  data.local_aggregate.clear();
  data.aggregate.clear();

  if (qpblockp.set_model_data(block)){
    delete block;
    block=qpblockp.generate_summodel_data();
    block->set_cbout(this);
  }
  assert(block);
  block->clear();
  
  int err=0;
  bool bundlehandler_existed=(bundlehandler!=0);
  if (sumbundle_mode(mode,bh,aft)){
    if (cb_out()){
      get_out()<<"**** ERROR ConeModel::start_sumaugmodel(): sumbundle_mode() failed"<<std::endl;
    }
    err++;
  }
  else if (((bundlehandler_existed)||(mode!=SumBundle::unavailable))
	   &&(modelmap.size()>0)){
    if (mode==SumBundle::root){
      // set which submodels should be asked to contribute how
      if (model_selection->select_models(modelmap)){
	if (cb_out()){
	  get_out()<<"*** WARNING: SumModel::start_sumaugmodel(): model_selection->select_models() failed"<<std::endl;
	}
	err++;
      }
    }
    else {
      // ask all to follow the own mode
      for (SumModel::ModelMap::iterator it=modelmap.begin();it!=modelmap.end();it++){
	it->second->set_suggested_mode()=mode;
      }
    }
  }

  // call all submodels to start their augmented models or to contribute
  local_augmodel.clear();
  Integer i=0; 
  for (ModelMap::iterator it=modelmap.begin();it!=modelmap.end();it++,i++){
    int retval=it->second->model()->start_sumaugmodel(qpblockp,cand_id,cand_y,indices,bundlehandler,it->second->set_suggested_mode());
    if (retval){
      if (cb_out(0)){
	get_out()<<"\n**** WARNING: SumModel::start_sumaugmodel(...): start_sumaugmodel failed for submodel i="<<i<<" and returned "<<retval<<std::endl;
      }
      err++;
    }
    else if (qpblockp.get_model_data_ptr()){
      local_augmodel.push_back(it->second);
      if (block->append(qpblockp.get_model_data_ptr())){
	if (cb_out()){
	  get_out()<<"*** ERROR: SumModel::start_sumaugmodel(): block->append() failed for submodel i="<<i<<std::endl;
	}
	err++;
      }
    }
  }

  //--- check if a decision between local and global models is required, and if so decide
  if ((bundlehandler)&&(data.sumbundle.has_working_roots())){
    if (bundlehandler->start_augmodel(qpblockp,*block,cand_id,cand_y,indices)){
      if (cb_out()){
	get_out()<<"*** ERROR: SumModel::start_sumaugmodel(): bundlehandler->start_augmodel() failed"<<std::endl;
      }
      err++;
    }
  }
  else if ((mode==SumBundle::child)&&(bundlehandler->add_contributions())){
    if (cb_out()){
      get_out()<<"*** ERROR: SumModel::start_sumaugmodel(): add_contributions() failed"<<std::endl;
    }
    err++;
  }  

  //--- check if there was some some model at all
  if ((block->get_bundle().size()==0)&&(block->get_constant_minorant().empty())) {
    qpblockp.clear_model_data_ptr();
  }
  else {
    if (cb_out(2)){
      get_out()<<" SMxdim="<<block->get_bundle().size();
    }
    qpblockp.set_model_data(block);
    assert(qpblockp.get_model_data_ptr()==block);
  }
    
  if (cb_out(10)){
    get_out()<<"\n  leaving  SumModel::start_sumaugmodel with return value "<<err<<std::endl;
  }

  return err;
}
 
// *****************************************************************************
//                                make_model_aggregate
// *****************************************************************************

// SumBundleHandler is NULL if the SumBundle is not used by the parent
// otherwise it holds information about the bundle_coefficients or mutliplier adjustments

int SumModel::make_model_aggregate(bool& increased, bool fixed)  
{
  if (cb_out(10)){
    get_out()<<"\n  entering  SumModel::make_model_aggregate"<<std::endl;
  }
  
  data.local_aggregate.clear();
  data.aggregate.clear();

  if (bundlehandler) {
    bool local_increased=increased;
    if (bundlehandler->make_model_aggregate(local_increased,fixed)){
      if (cb_out(0)){
	get_out()<<"\n**** ERROR:  SumModel::make_model_aggregate(...): make_model_aggregate failed for bundlehandler"<<std::endl;
      }
      if (cb_out(10)){
	get_out()<<"\n  leaving  SumModel::make_model_aggregate with return value 1"<<std::endl;
      }
      
      return 1;
    }
    increased |=local_increased;
  }


  Integer i=0; 
  for (ModelMap::iterator it=modelmap.begin();it!=modelmap.end();it++,i++){
    bool local_increased=increased;
    int retval=it->second->model()->make_model_aggregate(local_increased,fixed);
    if (retval){
      if (cb_out(0)){
	get_out()<<"\n**** ERROR:  SumModel::make_model_aggregate(..): make_model_aggregate failed for submodel i="<<i<<" and returned "<<retval<<std::endl;
	}
      if (cb_out(10)){
	get_out()<<"\n  leaving  SumModel::make_model_aggregate with return value 1"<<std::endl;
      }
      
      return 1;
    }
    increased |=local_increased;
  }


  if (increased) {
    data.center_ub_mid=-1;
  }
  else { //collect the aggregate

    //first the local part
    if (data.sumbundle.get_local_model_aggregate(data.local_aggregate)){
      if (cb_out(0)){
	get_out()<<"\n**** ERROR:  SumModel::make_model_aggregate(...): get_model_aggregate failed for bundlehandler"<<std::endl;
      }
      if (cb_out(10)){
	get_out()<<"\n  leaving  SumModel::make_model_aggregate with return value 1"<<std::endl;
      }
      
      return 1;
      
    }
    
    Integer i=0; 
    for (ModelMap::iterator it=modelmap.begin();it!=modelmap.end();it++,i++){
      int retval=it->second->model()->get_model_aggregate(it->second->aggregate_id,data.local_aggregate,false);
      if (retval){
	if (cb_out(0)){
	  get_out()<<"\n**** ERROR:  SumModel::make_model_aggregate(...): get_model_aggregate failed for submodel i="<<i<<" and returned "<<retval<<std::endl;
	}
	if (cb_out(10)){
	  get_out()<<"\n  leaving  SumModel::make_model_aggregate with return value 1"<<std::endl;
	}
	
	return 1;
      }
    }
    
    //then add the contributed part
    data.aggregate=data.local_aggregate;
    if (data.sumbundle.get_contributed_model_aggregate(data.aggregate)){
      if (cb_out(0)){
	get_out()<<"\n**** ERROR:  SumModel::make_model_aggregate(...): get_model_aggregate failed for bundlehandler"<<std::endl;
	}
      if (cb_out(10)){
	get_out()<<"\n  leaving  SumModel::make_model_aggregate with return value 1"<<std::endl;
      }
      
      return 1;
      
    }

    data.aggregate_id++;

  }
      
  if (cb_out(10)){
    get_out()<<"\n  leaving  SumModel::make_aug_linmodel with return value 0"<<std::endl;
  }
  return 0;
} 


// *****************************************************************************
//                                 lb_function
// *****************************************************************************

//returns a *quick* lower bound for the function value at y
//(eg by a previous subgradient)

Real SumModel::lb_function(Integer y_id,const Matrix& y)
{
  if (cb_out(10)){
    get_out()<<"\n  entering  SumModel::lb_function"<<std::endl;
  }
  if (data.dim<0)
    data.dim=y.dim();
  assert(data.dim==y.dim());

  if (data.aggregate.valid()){
    if (cb_out(10)){
      get_out()<<"\n  leaving  SumModel::lb_function, using aggregate"<<std::endl;
    }
    return data.aggregate.evaluate(y_id,y);
  }

  Real lb=0.;
  for (ModelMap::iterator it=modelmap.begin();it!=modelmap.end();it++){
    lb+=it->second->model()->lb_function(y_id,y);
  }

  if (cb_out(10)){
    get_out()<<"\n  leaving  SumModel::lb_function"<<std::endl;
  }
  return lb;
}


// *****************************************************************************
//                            get_function_minorant
// *****************************************************************************


int SumModel::get_function_minorant(MinorantPointer& minorant,
				    const AffineFunctionTransformation* aft)
{
  if (cb_out(10)){
    get_out()<<"\n  entering  SumModel::get_function_minorant"<<std::endl;
  }    

  if (!cand_minorant.valid()){
    cand_minorant.clear();
    Integer i=0; 
    for (ModelMap::iterator it=modelmap.begin();it!=modelmap.end();it++,i++){
      int retval=it->second->model()->get_function_minorant(cand_minorant);
      if (retval){
	if (cb_out()){
	  get_out()<<"**** ERROR SumModel::get_function_minorant: get_function_minorant failed for function "<<i<<" and returned "<<retval<<std::endl;
	}
	if (cb_out(10)){
	  get_out()<<"\n  leaving  SumModel::get_function_minorant with return value "<<retval<<std::endl;
	}    
	return retval;
      }
    }
  }
    
  if (aft) {
    int retval=aft->transform_minorant(minorant,cand_minorant,1.);
    if (retval){
      if (cb_out()) 
	get_out()<<"**** ERROR: SumModel::get_function_minorant(...): transforming minorant information failed and returned "<<retval<<std::endl;
      if (cb_out(10)){
	get_out()<<"\n  leaving  SumModel::get_function_minorant with return value "<<1<<std::endl;
      }    
      return 1;
    }    
  }
  else {
    assert(cand_minorant.valid());
    cand_minorant.get_minorant(minorant);
  }
   
  if (cb_out(10)){
    get_out()<<"\n  leaving  SumModel::get_function_minorant with return value "<<0<<std::endl;
  }
 
  return 0;
}


// *****************************************************************************
//                            get_center_minorant
// *****************************************************************************


  int SumModel::get_center_minorant(MinorantPointer& minorant,
				    const AffineFunctionTransformation* aft)
{
  if (cb_out(10)){
    get_out()<<"\n  entering  SumModel::get_center_minorant"<<std::endl;
  }    

  if (!center_minorant.valid()){
    center_minorant.clear();
    Integer i=0; 
    for (ModelMap::iterator it=modelmap.begin();it!=modelmap.end();it++,i++){
      int retval=it->second->model()->get_center_minorant(center_minorant);
      if (retval){
	if (cb_out()){
	  get_out()<<"**** ERROR SumModel::get_center_minorant: get_center_minorant failed for function "<<i<<" and returned"<<retval<<std::endl;
	}
	if (cb_out(10)){
	  get_out()<<"\n  leaving  SumModel::get_center_minorant with return value "<<retval<<std::endl;
	}    
	return retval;
      }
    }
  }

  if (aft){
    int retval=aft->transform_minorant(minorant,center_minorant,1.);
    if (retval){
      if (cb_out()) 
	get_out()<<"**** ERROR: SumModel::get_center_minorant(...): transforming minorant information failed and returned "<<retval<<std::endl;
      if (cb_out(10)){
	get_out()<<"\n  leaving  SumModel::get_center_minorant with return value "<<1<<std::endl;
      }    
      return 1;
    }
  }
  else {
    assert(center_minorant.valid());
    center_minorant.get_minorant(minorant);
  }
   
  if (cb_out(10)){
    get_out()<<"\n  leaving  SumModel::get_center_minorant with return value "<<0<<std::endl;
  }

  return 0;
}


// *****************************************************************************
//                                adjust_multiplier
// *****************************************************************************

  int SumModel::adjust_multiplier(bool& values_may_have_changed)  
{
  if (cb_out(10)){
    get_out()<<"\n  entering  SumModel::adjust_multplier"<<std::endl;
  }

  if (bundlehandler){
    bool maybe_changes=values_may_have_changed;
    if (bundlehandler->adjust_multiplier(maybe_changes)){
      if (cb_out()){
	get_out()<<"**** ERROR SumModel::adjust_multiplier(.): bundlehandler->adjust_multilier(.) failed"<<std::endl;
      }
      if (cb_out(10)){
	get_out()<<"\n  leaving  SumModel::adjust_multiplier(.) with return value 1"<<std::endl;
      }
      return 1;
    }
    values_may_have_changed |= maybe_changes;    
  }

  Integer i=0; 
  for (ModelMap::iterator it=modelmap.begin();it!=modelmap.end();it++,i++){
    bool maybe_changes=values_may_have_changed;
    int retval=it->second->model()->adjust_multiplier(maybe_changes);
    if (retval){
      if (cb_out(0)){
	get_out()<<"\n**** WARNING: SumModel::adjust_multplier(.): adjust_multiplier(.) failed for submodel i="<<i<<" and returned "<<retval<<std::endl;
      }
      clear_model();
      values_may_have_changed=true;
      
      if (cb_out(10)){
	get_out()<<"\n  entering  SumModel::addjust_multplier with return value 1"<<std::endl;
      }
      return 1;
    }
    values_may_have_changed |= maybe_changes;
  }

  if (values_may_have_changed){
    clear_model();
  }

  if (cb_out(10)){
    get_out()<<"\n  leaving  SumModel::adjust_multplier with return value 0"<<std::endl;
  }
  return 0;
}
 


// *****************************************************************************
//                                clear_model
// *****************************************************************************

void SumModel::clear_model(bool discard_minorants_only)  
{
  if (cb_out(10)){
    get_out()<<"\n  entering  SumModel::clear_model"<<std::endl;
  }

  Integer i=0;
  for (ModelMap::iterator it=modelmap.begin();it!=modelmap.end();it++,i++){
    SumBundle::Mode inact=SumBundle::inactive;
    if (it->second->model()->sumbundle_mode(inact,bundlehandler)){
      if (cb_out()) 
	get_out()<<"**** ERROR: SumModel::clear_model(): sumbundle_mode failed for model "<<i<<std::endl;
    }      
    
    it->second->aggregate_id=-1;
  }
  if (bundlehandler){
    bundlehandler->clear_model();
  }

  local_augmodel.clear();
  data.clear_model(discard_minorants_only);
  center_minorant.clear();
  cand_minorant.clear();

  if (cb_out(10)){
    get_out()<<"\n  leaving  SumModel::clear_model"<<std::endl;
  }

}

// *****************************************************************************
//                                clear_aggregates
// *****************************************************************************

void SumModel::clear_aggregates()  
{
  if (cb_out(10)){
    get_out()<<"\n  entering  SumModel::clear_aggregates"<<std::endl;
  }

  Integer i=0;
  for (ModelMap::iterator it=modelmap.begin();it!=modelmap.end();it++,i++){
    SumBundle::Mode inact=SumBundle::inactive;
    if (it->second->model()->sumbundle_mode(inact,bundlehandler)){
      if (cb_out()) 
	get_out()<<"**** ERROR: SumModel::clear_aggregates(): sumbundle_mode failed for model "<<i<<std::endl;
    }      
    
    it->second->aggregate_id=-1;
  }

  if (bundlehandler){
    bundlehandler->clear_aggregates();
  }
  data.clear_aggregates();
  local_augmodel.clear();

  if (cb_out(10)){
    get_out()<<"\n  leaving  SumModel::clear_aggregates"<<std::endl;
  }

}

// *****************************************************************************
//                         set_bundle_parameters      
// *****************************************************************************

int SumModel::set_bundle_parameters(const BundleParameters& bp)
{
  const SumModelParametersObject* smp=dynamic_cast<const SumModelParametersObject*>(&bp);
  if (smp){
    delete model_selection;
    model_selection=dynamic_cast<SumModelParametersObject*>(smp->clone_BundleParameters());
    assert(model_selection);
  }
  else {
    assert(model_selection);
    return model_selection->init(bp);
  }
  return 0;
}

// *****************************************************************************
//                         get_bundle_parameters      
// *****************************************************************************

///  see SumBlockModel::get_bundle_parameters()  
const BundleParameters* SumModel::get_bundle_parameters() const 
{return model_selection;}

// *****************************************************************************
//                         get_bundle_parameters      
// *****************************************************************************

/// set output and outputlevel of warnings and errors recursively, see CBout
void SumModel::set_out(std::ostream* o,int pril)
  {
    SumBlockModel::set_out(o,pril);
    for (ModelMap::iterator it=modelmap.begin();it!=modelmap.end();it++){
      it->second->model()->set_out(o,pril);
    }
    if (block)
      block->set_cbout(this,-1);
    if (model_selection)
      model_selection->set_cbout(this,0);
  }


}

