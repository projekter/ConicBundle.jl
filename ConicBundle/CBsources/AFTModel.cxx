/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/AFTModel.cxx
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
#include "AFTModel.hxx"

using namespace CH_Matrix_Classes;

namespace ConicBundle {

// *****************************************************************************
//                             generate_aft_cand()
// *****************************************************************************

// this should be the only place where a new point number can be assigned
// and the transformed coordinates are actually computed

//given on input: center_mid, center_id, center_y
//returned on output: center_mid, center_pid, 
//       pointer to center_y (arg changes not) or to center_aft (arg changes)

const CH_Matrix_Classes::Matrix& 
AFTModel::generate_aft_cand(
			    CH_Matrix_Classes::Integer cand_id,
			    const CH_Matrix_Classes::Matrix& cand_y
			    )
{
  //check if the point exists already as candidate for the current AFT
  if ((cand_id==data.cand_id)&&(data.cand_mid==data.modification_id)){
    //then all is already set correctly
    if (aft->argument_changes()){
      return data.aft_cand;
    }
    else
      return cand_y;
  }
  //check if the point exists already as center for the current AFT
  if ((cand_id==data.center_id)&&(data.center_mid==data.modification_id)){
    data.cand_id=cand_id;
    data.cand_pid=data.center_pid;
    data.cand_offset=data.center_offset;
    data.cand_mid=data.modification_id;

    data.cand_ub_mid=data.center_ub_mid;
    data.cand_ub_fid=data.center_ub_fid;
    data.cand_ub=data.center_ub;
    data.cand_relprec=data.center_relprec;

    cand_minorant.clear();

    if (aft->argument_changes()){
      aft->copy_traforows(data.aft_cand,data.aft_center);
      return data.aft_cand;
    }		
    data.aft_cand.init(0,0,0.);
    return cand_y;
  }

  //we assume it gets new coordinates and compute these
  cand_minorant.clear();
  data.cand_id=cand_id;
  data.cand_pid=++data.aftpoint_id;
  data.cand_ub_mid=-1;   //no evaluation information yet
  data.cand_ub_fid=-1;   //no evaluation information yet

  if (aft->get_fun_coeff()==0.){
    data.cand_ub=data.cand_offset;
    data.cand_ub_mid=data.modification_id;
    data.cand_relprec=1e-12;
  }

  data.cand_mid=data.modification_id;      
  return  aft->transform_argument(data.aft_cand,data.cand_offset,cand_y); 
}

// *****************************************************************************
//                       check_function_id_change()
// *****************************************************************************

int AFTModel::check_function_id_change(Integer function_id)
{
  if (function_id<data.function_id){
    if (cb_out())
      get_out()<<"**** ERROR AFTModel::check_function_id_change: input function_id="<<function_id<<" is smaller than stored data.function_id="<<data.function_id<<std::endl;
    return 1;
  }
  if (data.function_id<function_id){
    data.function_id=function_id;
    if (data.center_mid==data.modification_id)
      data.center_mid++;
    if (data.center_ub_mid==data.modification_id)
      data.center_ub_mid++;
    if (data.cand_mid==data.modification_id)
      data.cand_mid++;
    if (data.cand_ub_mid==data.modification_id)
      data.cand_ub_mid++;
    data.modification_id++;
    center_minorant.clear();
    cand_minorant.clear();
  }
  return 0;
}
  
  // *****************************************************************************
//                              AFTModel()
// *****************************************************************************

  AFTModel::AFTModel(SumBlockModel* in_model, 
		     AffineFunctionTransformation* inaft,
		     Integer start_modification_id,
		     bool in_model_is_owner,
		     CBout* cb):
    CBout(cb),SumBlockModel(cb),aft(0),block(0)
{
  model=in_model;
  model_is_owner=in_model? in_model_is_owner:false;
  clear(inaft,start_modification_id);

  CBout::set_cbout(in_model);
  aft->set_cbout(this);
  block=0;
}
  
// *****************************************************************************
//                              ~AFTModel()
// *****************************************************************************

AFTModel::~AFTModel()
{
  if (bundlehandler){
    SumBundle::Mode inact=SumBundle::inactive;
    if (cond_model_trafo()->sumbundle_mode(inact)){
      if (cb_out()) 
	get_out()<<"**** ERROR: AFTModel::~AFTModel(...): sumbundle_mode failed"<<std::endl;
    }      
    assert(!data.sumbundle.active());
    delete bundlehandler;
    bundlehandler=0;
  }
 
  if ((aft)&&(aft->get_model_calls_delete()))
    delete aft;
  aft=0;
  delete block;
  block=0;

  SumBlockModel::clear();
}

// *****************************************************************************
//                              clear()
// *****************************************************************************

  void AFTModel::clear(AffineFunctionTransformation* inaft,Integer stmodid)
{
  
  if (bundlehandler){
    SumBundle::Mode inact=SumBundle::inactive;
    if (cond_model_trafo()->sumbundle_mode(inact)){
      if (cb_out()) 
	get_out()<<"**** ERROR: AFTModel::clear(...): sumbundle_mode failed"<<std::endl;
    }      
    assert(!data.sumbundle.active());
    delete bundlehandler;
    bundlehandler=0;
  }

  SumBlockModel::clear();
 
  if (aft==0){
    if (inaft){
      aft=inaft;
    }
    else {
      aft=new AffineFunctionTransformation;
    }
  }
  else if (inaft) {
    if (aft->get_model_calls_delete())
      delete aft;
    aft=inaft;
  }
  data.clear(stmodid);

  ncalls=0;

  cand_minorant.clear();
  center_minorant.clear();

  delete block;
  block=0;

  precomputed.clear();
  qp_indices.init(0,0,Integer(0));
}

// *****************************************************************************
//                            eval_function
// *****************************************************************************

//evaluates the objective function in $y$
//if evaluated by an iterative method that provides upper and lower bounds,
//it may stop when the lower bound (lb) is above the nullstep_bound
//evalutation may also stop, if (ub-lb)<relprec*(|ub|+1) is satisfied
//returns:  0 ... if all is ok, use get_function_sol or do_step
//          1 ... if solution could not be computed to desired precision

  int AFTModel::eval_function(Integer& ub_fid,
			      Real& ub,
			      Integer y_id,
			      const Matrix& iny,
			      Real nullstep_bound,
			      Real relprec)
{
  if (cb_out(10)){
    get_out()<<"\n  entering  AFTModel::eval_function"<<std::endl;
  }    
  assert(relprec>0.);

  CH_Tools::Microseconds start_eval=clock.time();

  if (data.dim<0)
    data.dim=iny.dim();
  assert(data.dim==iny.dim());

  const Matrix& cand_y=generate_aft_cand(y_id,iny);

  if ((aft->get_fun_coeff()==0.)||(model==0)){
    //trivial solution, don't do the recursive call
    data.cand_ub=0.;
    data.cand_relprec=1e-10;
    data.cand_ub_fid=-1;
    data.cand_ub_mid=data.modification_id;
    ub_fid=data.modification_id;
    ub=data.cand_offset;
    //no need to consider the sumbundle if this just mimicks a linear function,
    //because a linear function is never included in a sumbundle
    if (cb_out(10)){
      get_out()<<"\n  leaving  AFTModel::eval_function (old computation results suffice)"<<std::endl;
    }    
    preeval_time+=clock.time()-start_eval;
    return 0;
  }

  if (bundlehandler)
    bundlehandler->clear_cand_minorants();
    
  nullstep_bound=(nullstep_bound-data.cand_offset)/aft->get_fun_coeff();
  
  if ((data.cand_ub_fid==data.function_id)&&
      (data.cand_ub_mid==data.modification_id)&&
      (data.cand_relprec<=relprec)){
    //data was already computed once for this point with same or higher precision
    ub= aft->objective_value(data.cand_offset,data.cand_ub);
    if (cb_out(10)){
      get_out()<<"\n  leaving  AFTModel::eval_function (old computation results suffice)"<<std::endl;
    }    
    preeval_time+=clock.time()-start_eval;
    return 0;
  }
    
  data.cand_relprec=relprec;
 
  ncalls++;
  preeval_time+=clock.time()-start_eval;
  start_eval=clock.time();
  int retval=cond_model_trafo()->eval_function(data.cand_ub_fid,data.cand_ub,data.cand_pid,cand_y,nullstep_bound,relprec);
  eval_time+=clock.time()-start_eval;
  start_eval=clock.time();

  if (retval<0){
    if (cb_out()){
      get_out()<<"**** WARNING AFTModel::eval_function: eval_function returned"<<retval<<std::endl;
    }
  }    
  if (data.cand_ub_fid<0){
    if (cb_out()){
      get_out()<<"**** WARNING AFTModel::eval_function: eval_function returned negative function_id="<<data.cand_ub_fid<<std::endl;
    }
    data.cand_ub_fid=-1;
  }    
  else {
    check_function_id_change(data.cand_ub_fid);
  }
  if (retval>0){
    if (cb_out()){
      get_out()<<"**** ERROR AFTModel::eval_function: eval_function failed and returned"<<retval<<std::endl;
    }
    if (cb_out(10)){
      get_out()<<"\n  leaving  AFTModel::eval_function with return value "<<retval<<std::endl;
    }    
    posteval_time+=clock.time()-start_eval;
    return retval;
  }    
  
  //all information collected
  ub=aft->objective_value(data.cand_offset,data.cand_ub);
  ub_fid=data.cand_ub_mid=data.modification_id;

  //update old minorants if necessary
  if (data.max_old_minorants>0){
    MinorantPointer tmp;
    if (get_function_minorant(tmp)){
      if (cb_out())
	get_out()<<"**** WARNING AFTModel::eval_function: get_function_minorant(.) failed for adding olld minorant"<<std::endl;
    }
    else if (data.store_old_minorant(tmp)){
      if (cb_out())
	get_out()<<"**** WARNING FunctionModel::eval_function(): data.store_old_minorant() failed"<<std::endl;
    }
  }

  //update contribution to the parent sumbundle if necessary
  if ((bundlehandler)&&(data.sumbundle.has_contributions())){
    if (bundlehandler->contribute_new_minorants()){
	if (cb_out())
	  get_out()<<"**** ERROR AFTModel::eval_function: bundlehandler->contribute_new_minorants failed"<<std::endl;
    }
  }

  posteval_time+=clock.time()-start_eval;
  if (cb_out(10)){
    get_out()<<"\n  leaving  AFTModel::eval_function with return value "<<retval<<std::endl;
  }    

  return retval;
}

// *****************************************************************************
//                            eval_model
// *****************************************************************************

//evaluates the objective model in $y$
//if evaluated by an iterative method that provides upper and lower bounds,
//it may stop when the lower bound (lb) is above the nullstep_bound
//evalutation may also stop, if (ub-lb)<relprec*(|ub|+1) is satisfied
//returns:  0 ... if all is ok, use get_model_sol or do_step
//          1 ... if solution could not be computed to desired precision

  int AFTModel::eval_model(Real& lb,
			   Integer y_id,
			   const Matrix& iny,
			   Real relprec)
{
  if (cb_out(10)){
    get_out()<<"\n  entering  AFTModel::eval_model"<<std::endl;
  }    
  assert(relprec>0.);

  assert(data.dim==iny.dim());
  const Matrix& cand_y=generate_aft_cand(y_id,iny);

  if ((aft->get_fun_coeff()==0.)||(model==0)){
    //trivial solution, don't do the recursive call
    lb=data.cand_offset;
    if (cb_out(10)){
      get_out()<<"\n  leaving  AFTModel::eval_model (function coefficient 0.)"<<std::endl;
    }    
    return 0;
  }

  assert((bundlehandler==0)||(!data.sumbundle.has_roots()));
    
 
  Real cand_lb;
  CH_Tools::Microseconds start_eval=clock.time();
  int retval=cond_model_trafo()->eval_model(cand_lb,data.cand_pid,cand_y,relprec);
  evalmodel_time+=clock.time()-start_eval;

  if (retval<0){
    if (cb_out()){
      get_out()<<"**** WARNING AFTModel::eval_model: eval_model returned"<<retval<<std::endl;
    }
  }    
  if (retval>0){
    if (cb_out()){
      get_out()<<"**** ERROR AFTModel::eval_model: eval_model failed and returned"<<retval<<std::endl;
    }
    if (cb_out(10)){
      get_out()<<"\n  leaving  AFTModel::eval_model with return value "<<retval<<std::endl;
    }    
    return retval;
  }    
  
  //all information collected
  lb=aft->objective_value(data.cand_offset,cand_lb);
  
  if (cb_out(3)){
    get_out().precision(12);
    get_out()<<" AFTModval="<<lb;
  }

  if (cb_out(10)){
    get_out()<<"\n  leaving  AFTModel::eval_model with return value "<<retval<<std::endl;
  }    

  return retval;
}


// *****************************************************************************
//                                update_model
// *****************************************************************************

int AFTModel::update_model(ModelUpdate model_update,
			   Integer center_id,
			   const Matrix& center_y,
			   Integer y_id,
			   const Matrix& iny,
			   Real model_maxviol,
			   BundleProxObject& H)
{
  if (cb_out(10)){
    get_out()<<"\n  entering  AFTModel::update_model"<<std::endl;
  }
  assert((y_id==data.center_id)||(y_id==data.cand_id));
  assert(bundlehandler==0);

  CH_Tools::Microseconds update_start=clock.time();

  if ((aft->get_fun_coeff()==0.)||(model==0)){
    //function is linear in y
    int retval=0;
     if (model_update==descent_step){
      retval=data.do_step(y_id);
      center_minorant.clear();
      if ((retval)&&(cb_out())) 
	get_out()<<"*** ERROR: AFTModel::update_model: function coefficient 0. and calling data.do_step() returned "<<retval<<std::endl;
    }
    updatemodel_time+=clock.time()-update_start;
    if (cb_out(10)){
      get_out()<<"\n  leaving  AFTModel::update_model with return value "<<retval<<std::endl;
    }
    return retval;
  }


  if ((data.cand_ub_mid!=data.modification_id)&&(data.center_ub_mid!=data.modification_id)){
    if (cb_out(0)){
      get_out()<<"\n**** WARNING: AFTModel::update_model(...) called without nay available subgradient information"<<std::endl;
    }
    if (cb_out(10)){
      get_out()<<"\n  leaving  AFTModel::update_model with return value 1"<<std::endl;
    }
    return 1;
  }
  Integer cenid;
  const Matrix* my_ceny;
  if (center_id==data.center_id){
    cenid=data.center_pid;
    if (aft->argument_changes())
      my_ceny=&(data.aft_center);
    else
      my_ceny=&center_y;
  }
  else {
    cenid=data.cand_pid;
    if (aft->argument_changes())
      my_ceny=&(data.aft_cand);
    else
      my_ceny=&iny;
  }
  Integer canid;
  const Matrix* my_cany;
  if (y_id==data.center_id){
    canid=data.center_pid;
    if (aft->argument_changes())
      my_cany=&(data.aft_center);
    else
      my_cany=&center_y;
  }
  else {
    canid=data.cand_pid;
    if (aft->argument_changes())
      my_cany=&(data.aft_cand);
    else
      my_cany=&iny;
  }

  Real factor=data.function_factor*aft->get_fun_coeff();
  H.apply_factor(1./factor);
  int retval1=cond_model_trafo()->update_model(model_update,cenid,*my_ceny,canid,*my_cany,model_maxviol/factor,H);
  if (retval1){
    if (cb_out()) 
      get_out()<<"*** ERROR: AFTModel::update_model: calling update_model returned "<<retval1<<std::endl;
  }  
  H.apply_factor(factor);

  int retval2=0;
  if (model_update==descent_step){
    retval2=data.do_step(y_id);
    center_minorant.clear();
    if ((retval2)&&(cb_out())) 
      get_out()<<"*** ERROR: AFTModel::update_model: calling data.do_step() returned "<<retval2<<std::endl;
  }

  updatemodel_time+=clock.time()-update_start;

  if (cb_out(10)){
    get_out()<<"\n  leaving  AFTModel::update_model with return value "<<retval1+retval2<<std::endl;
  }
  return retval1+retval2;
}

// *****************************************************************************
//                                update_model
// *****************************************************************************

int AFTModel::update_model(ModelUpdate model_update,
			   Integer center_id,
			   const Matrix& center_y,
			   Integer y_id,
			   const Matrix& iny,
			   Real model_maxviol,
			   BundleProxObject& H,
			   Real & model_deviation,
			   Real & model_curvature)
{
  if (cb_out(10)){
    get_out()<<"\n  entering  AFTModel::update_model"<<std::endl;
  }
  assert((y_id==data.center_id)||(y_id==data.cand_id));

  CH_Tools::Microseconds update_start=clock.time();

  if ((aft->get_fun_coeff()==0.)||(model==0)){
    //function is linear in y
    model_deviation=0;
    model_curvature=0;
    int retval=0;
    if (model_update==descent_step){
      retval=data.do_step(y_id);
      center_minorant.clear();
      if ((retval)&&(cb_out())) 
	get_out()<<"*** ERROR: AFTModel::update_model: function coefficient 0. and calling data.do_step() returned "<<retval<<std::endl;
    }
    updatemodel_time+=clock.time()-update_start;
    if (cb_out(10)){
      get_out()<<"\n  leaving  AFTModel::update_model with return value "<<retval<<std::endl;
    }
    return retval;
  }

  if (bundlehandler){
    if (bundlehandler->update_model(model_update,center_id,center_y,y_id,iny,model_maxviol,H)){
      if (cb_out(0)){
	get_out()<<"\n**** WARNING: AFTModel::update_model(...): bundlehandler->update_model(.) failed"<<std::endl;
      }
      updatemodel_time+=clock.time()-update_start;
      if (cb_out(10)){
	get_out()<<"\n  leaving  AFTModel::update_model with return value 1"<<std::endl;
      }
      return 1;
    }
  }

  if ((data.cand_ub_mid!=data.modification_id)&&(data.center_ub_mid!=data.modification_id)){
    if (cb_out(0)){
      get_out()<<"\n**** WARNING: AFTModel::update_model(...) called without any available subgradient information"<<std::endl;
    }
    updatemodel_time+=clock.time()-update_start;
    if (cb_out(10)){
      get_out()<<"\n  leaving  AFTModel::update_model with return value 1"<<std::endl;
    }
    return 1;
  }
  Integer cenid;
  const Matrix* my_ceny;
  if (center_id==data.center_id){
    cenid=data.center_pid;
    if (aft->argument_changes())
      my_ceny=&(data.aft_center);
    else
      my_ceny=&center_y;
  }
  else {
    cenid=data.cand_pid;
    if (aft->argument_changes())
      my_ceny=&(data.aft_cand);
    else
      my_ceny=&iny;
  }
  Integer canid;
  const Matrix* my_cany;
  if (y_id==data.center_id){
    canid=data.center_pid;
    if (aft->argument_changes())
      my_cany=&(data.aft_center);
    else
      my_cany=&center_y;
  }
  else {
    canid=data.cand_pid;
    if (aft->argument_changes())
      my_cany=&(data.aft_cand);
    else
      my_cany=&iny;
  }

  Real factor=data.function_factor*aft->get_fun_coeff();
  H.apply_factor(1./factor);
  int retval1=cond_model_trafo()->update_model(model_update,cenid,*my_ceny,canid,*my_cany,model_maxviol/factor,H,model_deviation,model_curvature);
  model_deviation*=factor;
  model_curvature*=factor;
  if (retval1){
    if (cb_out()) 
      get_out()<<"*** ERROR: AFTModel::update_model: calling update_model returned "<<retval1<<std::endl;
  }  
  H.apply_factor(factor);

  int retval2=0;
  if (model_update==descent_step){
    retval2=data.do_step(y_id);
    center_minorant.clear();
    if ((retval2)&&(cb_out())) 
      get_out()<<"*** ERROR: AFTModel::update_model: calling data.do_step() returned "<<retval2<<std::endl;
  }
  
  updatemodel_time+=clock.time()-update_start;

  if (cb_out(10)){
    get_out()<<"\n  leaving  AFTModel::update_model with return value "<<retval1+retval2<<std::endl;
  }
  return retval1+retval2;
}
  
// *****************************************************************************
//                             synchronize_ids
// *****************************************************************************


int AFTModel::synchronize_ids(Integer& new_center_ub_fid, 
			      Integer new_center_id,
			      Integer old_center_id,
			      Integer& new_cand_ub_fid, 
			      Integer new_cand_id,
			      Integer old_cand_id,
                              Integer& new_aggregate_id)

 {
  if (cb_out(10)){
    get_out()<<"\n  entering  AFTModel::synchronize_ids"<<std::endl;
  }

  Integer old_center_pid=data.center_pid;
  Integer old_cand_pid=data.cand_pid;

  int retval1=data.synchronize_ids(new_center_ub_fid,
				   new_center_id,
				   old_center_id,
				   new_cand_ub_fid,
				   new_cand_id,
				   old_cand_id,
                                   new_aggregate_id,0);

  if (retval1){
    if (cb_out()) 
      get_out()<<"*** ERROR: AFTModel::synchronize_ids: calling data.synchronize_ids returned "<<retval1<<std::endl;
  }

  int retval2=0;
  if (model!=0){ 
    if (aft->get_fun_coeff()==0.){
      Integer dummy_center_fid=-1;
      Integer dummy_cand_fid=-1;
      Integer dummy_aggr_id=-1;

      retval2=cond_model_trafo()->synchronize_ids(dummy_center_fid,Integer(-1),Integer(-1),
						  dummy_cand_fid,Integer(-1),Integer(-1),dummy_aggr_id);
      if (retval2){
	if (cb_out()) 
	  get_out()<<"*** ERROR: AFTModel::synchronize_ids: calling synchronize_ids returned "<<retval2<<std::endl;
      }
      data.function_id=-1;
      new_center_ub_fid=data.center_ub_fid=-1;
      new_cand_ub_fid=data.cand_ub_fid=-1;
    }
    else {
      retval2=cond_model_trafo()->synchronize_ids(data.center_ub_fid,
						  data.center_pid,
						  old_center_pid,
						  data.cand_ub_fid,
						  data.cand_pid,
						  old_cand_pid,
						  data.aggregate_fid);
    }
    if (retval2){
      if (cb_out()) 
	get_out()<<"*** ERROR: AFTModel::synchronize_ids: calling synchronize_ids returned "<<retval2<<std::endl;
    }

    if (data.center_ub_fid<0){
      new_center_ub_fid=data.center_ub_mid=data.center_ub_fid=-1;
    }
    if (data.cand_ub_fid<0){
      new_cand_ub_fid=data.cand_ub_mid=data.cand_ub_fid=-1;
    }
    data.function_id=max(-1,max(data.center_ub_fid,data.cand_ub_fid));
    if (data.aggregate_fid<0){
      data.aggregate_available=false;
      new_aggregate_id=data.aggregate_fid=-1;
    }
  }
  else {
    data.function_id=-1;
    data.center_ub_fid=-1;
    data.cand_ub_fid=-1;
  }

  if ((data.center_ub_fid!=data.function_id)||
      (data.center_ub_mid!=data.modification_id))
    center_minorant.clear();

  if ((data.cand_ub_fid!=data.function_id)||
      (data.cand_ub_mid!=data.modification_id))
    cand_minorant.clear();

  if (cb_out(10)){
    get_out()<<"\n  leaving  AFTModel::synchronize_ids with return value "<<retval1+retval2<<std::endl;
  }
  return retval1+retval2;
 }


// *****************************************************************************
//                             center_modified
// *****************************************************************************


  bool AFTModel::center_modified(Integer& function_modification_id,
				 Integer last_center_id)
 {
  if (cb_out(10)){
    get_out()<<"\n  entering  AFTModel::center_modified"<<std::endl;
  }
  bool modified=false;
  
  if (data.center_modified(function_modification_id,last_center_id)){
    modified=true;
  }  

  if ((aft->get_fun_coeff()!=0.)&&(model!=0)){
    Integer fid=data.center_ub_fid;
    bool model_modified = cond_model_trafo()->center_modified(fid,data.center_pid) ;
    modified |= model_modified;
    if (check_function_id_change(fid)){
      if (cb_out()){
	get_out()<<"**** ERROR in AFTModel::center_modified(): check_function_id_change failed"<<std::endl;
      }
    }
    if (model_modified){
      data.center_ub_mid=data.center_ub_fid=-1;
    }
    else{
      if (data.center_ub_fid>=0){
	data.center_ub_mid=data.modification_id;
	data.center_ub_fid=data.function_id;
      }
    }
  }
  else {
    assert(data.center_ub_fid==-1);
    assert(data.cand_ub_fid==-1);
    assert(data.aggregate_fid==-1);
    data.function_id=-1;
  }

  if (modified){
    center_minorant.clear();
    cand_minorant.clear();
  }
  function_modification_id=data.modification_id;

  if (cb_out(10)){
    get_out()<<"\n  leaving  AFTModel::center_modified with return value "<<modified;
  }
  return modified;
 }


// *****************************************************************************
//                             recompute_center
// *****************************************************************************


  int AFTModel::recompute_center(Integer& new_ub_fid,
				 Real& new_ub,
				 Integer center_id,
				 const Matrix& center_y,
				 bool accept_only_higher,
				 Real relprec)
 {
  if (cb_out(10)){
    get_out()<<"\n  entering  AFTModel::recompute_center"<<std::endl;
  }
  if (data.dim<0)
    data.dim=center_y.dim();
  assert(data.dim==center_y.dim());

  if ((aft->get_fun_coeff()==0.)||(model==0)){
    if((data.center_id!=center_id)||
       (data.center_mid!=data.modification_id)){
      generate_aft_cand(center_id,center_y);
      //here this also fills up the candidate with the function evaluation data
      data.do_step(center_id);
      center_minorant.clear();
    }
    new_ub=aft->objective_value(data.center_offset,0.);
    new_ub_fid=data.modification_id;
    if (cb_out(10)){
      get_out()<<"\n  leaving  AFTModel::recompute_center with available center and return value 0"<<std::endl;
    }
    return 0.;
  }

  if ((relprec<=0)&&((center_id!=data.center_id)||(data.center_relprec<=0))){
    Real rp=1e-3;
    if ((data.center_relprec>0)&&(rp>data.center_relprec))
      rp=data.center_relprec;
    if (cb_out()){
      get_out()<<"\n**** WARNING AFTModel::recompute_center: called with relprec="<<relprec<<"<=0 but no previous precision available for center_id="<<center_id<<", using default choice "<<rp<<std::endl;
    }
    relprec=rp;
  }
  if (relprec<=0)
    relprec=data.center_relprec;

  if ((accept_only_higher)&&
      ((data.center_id!=center_id)||(data.center_ub_mid!=data.modification_id))){
    if (cb_out()){
      get_out()<<"\n**** WARNING AFTModel::recompute_center: called with accept_only_higher=true but there is no valid previous value or data.center_id="<<data.center_id<<" differs from center_id="<<center_id<<"; this should never happen, setting accept_only_higher to false"<<std::endl;
    }
    accept_only_higher=false;
  }
  
  if ((data.center_id==center_id)&&
      (data.center_ub_fid==data.function_id)&&
      (data.center_ub_mid==data.modification_id)&&
      (data.center_mid==data.modification_id)&&
      (relprec<=data.center_relprec)
      ){
    //new computation should not change anything to the better
    new_ub=aft->objective_value(data.center_offset,data.center_ub);
    new_ub_fid=data.modification_id;
    if (cb_out(10)){
      get_out()<<"\n  leaving  AFTModel::recompute_center with available center and return value 0"<<std::endl;
    }
    return 0;
  }


  const Matrix& my_y=generate_aft_cand(center_id,center_y);
  
  Real ub=CB_minus_infinity;
  Integer ub_fid=-1;
  ncalls++;
  CH_Tools::Microseconds start_eval=clock.time();
  int retval=cond_model_trafo()->recompute_center(ub_fid,ub,data.cand_pid,my_y,accept_only_higher,relprec);
  eval_time+=clock.time()-start_eval;
  if (retval>0){
    if (cb_out(0)){
      get_out()<<"\n**** ERROR:  AFTModel::recompute_center(...): recompute_center failed and returned"<<retval<<std::endl;
    }
    if (cb_out(10)){
      get_out()<<"\n  leaving  AFTModel::recompute_center with return value 1"<<std::endl;
    }
    return 1;
  }
  if (check_function_id_change(ub_fid)){
    if (cb_out(0)){
      get_out()<<"\n**** ERROR:  AFTModel::recompute_center(...): check_function_id_change failed"<<std::endl;
    }
  }
  
  if ((!accept_only_higher)||(ub>data.center_ub)){
    data.cand_ub=ub;
    data.cand_relprec=relprec;
    data.cand_ub_mid=data.modification_id;
    data.cand_ub_fid=ub_fid;
    data.do_step(center_id);
    center_minorant.clear();
  }
  new_ub=aft->objective_value(data.center_offset,data.center_ub);
  new_ub_fid=data.modification_id;

  if (cb_out(10)){
    get_out()<<"\n  leaving  AFTModel::recompute_center with return value "<<retval<<std::endl;
  }

  return retval;
 }



// *****************************************************************************
//                             model_aggregate_modified
// *****************************************************************************


bool AFTModel::model_aggregate_modified(Integer last_aggr_id)
 {
  if (cb_out(10)){
    get_out()<<"\n  entering  AFTModel::model_aggregate_modified"<<std::endl;
  }

  if ((aft->get_fun_coeff()==0.)||(model==0)){
    if ((!data.aggregate_available)||(last_aggr_id<data.aggregate_id)) {
      if (cb_out(10)){
	get_out()<<"\n  leaving  AFTModel::model_aggreagte_modified with return value true"<<std::endl;
      }
      return true;
    }
    else {
      if (cb_out(10)){
	get_out()<<"\n  leaving  AFTModel::model_aggreagte_modified with return value false"<<std::endl;
      }
      return false;
    }
  }
    
  if (cond_model_trafo()->model_aggregate_modified(data.aggregate_fid)){
    data.aggregate_available=false;
  }

  if ((!data.aggregate_available)||(last_aggr_id<data.aggregate_id)) {
    if (cb_out(10)){
      get_out()<<"\n  leaving  AFTModel::model_aggreagte_modified with return value true"<<std::endl;
    }
    return true;
  }
  
  if (cb_out(10)){
    get_out()<<"\n  leaving  AFTModel::model_aggreagte_modified with return value false"<<std::endl;
  }
  return false;
 }


// *****************************************************************************
//                             provide_model_aggregate
// *****************************************************************************


  int AFTModel::provide_model_aggregate(Integer y_id,const Matrix& y)
 {
  if (cb_out(10)){
    get_out()<<"\n  entering  AFTModel::provide_model_aggregate"<<std::endl;
  }

  assert(data.dim==y.dim());

  if (data.aggregate_available) {
    if (cb_out(10)){
      get_out()<<"\n  leaving  AFTModel::provide_model_aggreagte with return value 0"<<std::endl;
    }
    return 0;
  }

  if ((aft->get_fun_coeff()==0.)||(model==0)){
    //it is stored directly in aft, nothing to do
    if (!data.aggregate_available)
      data.aggregate_id++;
    data.aggregate_available=true;
    data.local_aggregate.clear();
    data.aggregate.clear();
    if (cb_out(10)){
      get_out()<<"\n  leaving  AFTModel::provide_model_aggreagte with return value 0"<<std::endl;
    }
    return 0;
  }

  const Matrix& cand_y=generate_aft_cand(y_id,y);

  int retval=cond_model_trafo()->provide_model_aggregate(data.cand_pid,cand_y);
  if (retval){
    if (cb_out(0)){
      get_out()<<"\n**** ERROR:  AFTModel::provide_model_aggregate(...): provide_model_aggregate failed and returned "<<retval<<std::endl;
    }
    if (cb_out(10)){
      get_out()<<"\n  leaving  AFTModel::provide_model_aggregate with return value 1"<<std::endl;
    }
    return 1;
  }
  data.aggregate_id++;
  data.aggregate_available=true;
  data.local_aggregate.clear();
  data.aggregate.clear();
   
  if (cb_out(10)){
    get_out()<<"\n  leaving  AFTModel::provide_model_aggreagte with return value 0"<<std::endl;
  }
  return 0;
 }


// *****************************************************************************
//                         apply_modification
// *****************************************************************************

int AFTModel::apply_modification(bool& no_changes,
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

  if (cb_out(10)){
    get_out()<<"\n   entering AFTModel::apply_modification(.......) "<<std::endl;
  }

  assert ((data.dim<0)||(old_center_id<0)||(data.dim==old_center.dim()));
  assert((data.dim<0)||(data.dim==gsmdf.old_vardim()));
  assert((new_center_id<0)==(old_center_id<0));

  // ----- retrieve explicit requests -----
  int err=0;
  const AFTModification* aftmdf=0;

  FunObjModMap::const_iterator mapit;
  if (model_is_owner){
    mapit=funmdfmap.find(model->get_oracle_object());
    if (mapit!=funmdfmap.end()){
      aftmdf=mapit->second.get_aft_modification();
    }
  }
  else {
    mapit=funmdfmap.find(aft);
    if (mapit!=funmdfmap.end()){
      const OracleModification* omdf=mapit->second.get_oracle_modification();
      if (omdf){
	aftmdf=dynamic_cast<const AFTModification*>(omdf);
	if (aftmdf==0){
	  if (cb_out())
	    get_out()<<"**** ERROR AFTModel::apply_modification(.......): casting OracleModification* to AFTModification* failed"<<std::endl;
	  err++;
	}
      }
    }
  }

  // ----- check whether the transposed aft-map changes -----

  bool  minorant_trafo_differs;
  const GroundsetModification* my_gsmdf=aft->analyze_modification(minorant_trafo_differs,
								 aftmdf,
								 gsmdf);
  if (my_gsmdf==0){
    if (cb_out()){
      get_out()<<"**** ERROR FunctionModel::apply_modification(.......): aft->analyze_modification(...) failed"<<std::endl;
    }
    return 1;
  }

  // ----- determine the transformed center points ------
  const Matrix* my_old_center=&old_center;
  Integer my_old_center_pid=-1;
  const Matrix* my_new_center=&new_center;
  Integer my_new_center_pid=-1;
  Matrix new_aft_center;
  bool aft_changes=((aftmdf)&&(!aftmdf->no_modification()))||
    ((aftmdf==0)&&(!gsmdf.no_modification()));
  bool center_exists=!((data.center_ub_mid<0)||(old_center_id!=data.center_id)||(new_center_id<0));
  if (center_exists){
    my_old_center_pid=data.center_pid;
    if (aft->argument_changes()){
      my_old_center=&data.aft_center;
    }
    if (new_center_id!=old_center_id)
      my_new_center_pid=data.aftpoint_id+1;
    else 
      my_new_center_pid=my_old_center_pid;
    if ((aft->argument_changes())||aft_changes){
      my_new_center = & aft->modified_transform_argument(new_aft_center,new_center,aftmdf,gsmdf);
      if (my_new_center_pid==my_old_center_pid)
	my_new_center_pid=data.aftpoint_id+1;
    }
  }
    
  // ------ if the contribution changes, remove it first
  if ((!no_changes)|| aft_changes){
    if (bundlehandler){
      if (bundlehandler->remove_contributions()){
	if (cb_out()){
	  get_out()<<"**** ERROR AFTModel::apply_modification(.......): sumbundle->remove_contributions failed"<<std::endl;
	}
	err++;
      }
    }
  }


  // ----- call the transformed function
 
  int retval=0;
  if (model!=0){
    retval=cond_model_trafo()->apply_modification(
						  no_changes,
						  *my_gsmdf,funmdfmap,
						  my_new_center_pid,*my_new_center,
						  my_old_center_pid,*my_old_center
						  );
    if (retval){
      if (cb_out(0)){
	get_out()<<"\n**** WARNING:  AFTModel::apply_modification(..): apply_modification(..) failed and returned "<<retval<<std::endl;
      }
    }
  }

  if ((!gsmdf.no_modification())&& no_changes){
    if (data.apply_modification(gsmdf,0)){
      if (cb_out(0))
	get_out()<<"\n**** WARNING:  AFTModel::apply_modification(..): apply_modification(..) failed for data"<<std::endl;
      err++;
      no_changes=false;
    }
    data.aggregate_available=false;
    data.modification_id++;
    cand_minorant.clear();
    center_minorant.clear();
    precomputed.clear();
  }

  bool cl_model=false;
  if (!no_changes){
    cl_model=true;
    clear_model();
  }


  // ----- execute the changes in the aft and install the new center  
  if ((aft_changes)||(!no_changes)){
    data.aggregate_available=false;
    data.modification_id++;
  }

  if (aft_changes) {
    if (aft->apply_modification(aftmdf,gsmdf)){
      if (cb_out()){
	get_out()<<"**** ERROR AFTModel::apply_modification(.......):  aft->apply_modification(..) failed"<<std::endl;
      }
      err++;
    }
    cand_minorant.clear();
    center_minorant.clear();
    precomputed.clear();
  }

  if ((center_exists)&&(!cl_model)){
    if (my_old_center_pid!=my_new_center_pid){
      data.center_ub_fid=-1;
      data.center_id=-1;
      data.cand_id=-1;
      generate_aft_cand(new_center_id,new_center);
      assert(data.cand_pid==my_new_center_pid);
      data.center_id=data.cand_id;  // indeed center_id and not data.center_id!
      data.center_mid=data.modification_id;
      data.center_pid=data.cand_pid;
      data.aft_center=data.aft_cand;
      data.center_offset=data.cand_offset;
    }
    else {
      //nothing changed except possibly the modification_id 
      data.center_ub_fid=data.modification_id;
      data.center_mid=data.modification_id;
    }
  }
  else {
    data.center_ub_fid=-1;
    data.center_mid=-1;
    data.cand_mid=-1;
    //make sure a later new center is not mixed up with the virtual one of just now
    if (my_new_center_pid==data.aftpoint_id) {
      data.aftpoint_id++;
    }
  }

  data.dim=gsmdf.new_vardim();
  
  if (cb_out(10)){
    get_out()<<"\n   leaving AFTModel::apply_modification with return value "<<((retval)?retval:err)<<std::endl;
  }
  return (retval)?retval:err;  
}



// *****************************************************************************
//                              add_variable_metric
// *****************************************************************************


int AFTModel::add_variable_metric(VariableMetric& H,
				  Integer center_id,
				  const Matrix& center_y,
				  bool descent_step,
				  Real weightu,
				  Real model_maxviol,
				  const Indexmatrix* indices)
{
  if ((!H.employ_variable_metric())||(!H.get_use_local_metric())){
    if (SumBlockModel::add_variable_metric(H,center_id,center_y,descent_step,weightu,
					   model_maxviol,indices)){
      if (cb_out())
	get_out()<<"**** WARNING: AFTModel::add_variable_metric(): SumBlockModel::add_variable_metric() failed"<<std::endl;
      return 1;
    }
    return 0;
  }
  
  if (cb_out(10)){
    get_out()<<"\n  entering  AFTModel::add_variable_metric"<<std::endl;
  }
  assert(data.dim==center_y.dim());


  if ((aft->get_fun_coeff()==0.)||(model==0)){
    if (cb_out(10)){
      get_out()<<"\n  leaving  AFTModel::add_variable_metric with return value 0"<<std::endl;
    }
    return 0;
  }

  CH_Tools::Microseconds start_metric=clock.time();

  assert((center_id==data.center_id)&&(data.modification_id==data.center_ub_mid));

  if ((bundlehandler)&&(data.sumbundle.has_working_roots())){
    if (bundlehandler->add_variable_metric(H,center_id,center_y,
					   descent_step,weightu,model_maxviol,indices)){
      if (cb_out())
	get_out()<<"**** WARNING: AFTModel::add_variable_metric(...) bundlehandler->add_variable_metric(...) failed"<<std::endl;

      metric_time+=clock.time()-start_metric;

      if (cb_out(10)){
	get_out()<<"\n leaving  AFTModel::add_variable_metric with return value 1"<<std::endl;
      }
      return 1;
    }
  }

  H.push_aft(aft);

  Integer id=data.center_pid;

  const Matrix* my_y;
  if (aft->argument_changes())
    my_y=&(data.aft_center);
  else
    my_y=&center_y;

  int retval=cond_model_trafo()->add_variable_metric(H,id,*my_y,
						     descent_step,
						     weightu/aft->get_fun_coeff(),
						     model_maxviol/aft->get_fun_coeff(),
						     indices);

  H.pop_aft();

  metric_time+=clock.time()-start_metric;

  if (retval){
    if (cb_out(0)){
      get_out()<<"\n**** WARNING:  AFTModel::add_variable_metric(...): add_variable_metric failed and returned "<<retval<<std::endl;
    }

    if (cb_out(10)){
      get_out()<<"\n leaving  AFTModel::add_variable_metric with return value 1"<<std::endl;
    }
    return 1;
  }

  if (cb_out(10)){
    get_out()<<"\n leaving  AFTModel::add_variable_metric with return value 0"<<std::endl;
  }
  return 0;
}

// *****************************************************************************
//                         check_center_validity_by_candidate
// *****************************************************************************


  int AFTModel::check_center_validity_by_candidate(bool& cand_minorant_is_below,
                                                   Integer center_id,
						     const Matrix& center_y)
{
  if (cb_out(10)){
    get_out()<<"\n  entering  AFTModel::check_center_validity_by_candidate"<<std::endl;
  }    
  assert(data.dim==center_y.dim());
  
  if ((center_id!=data.center_id)||
      (data.cand_ub_mid!=data.modification_id)||
      (data.center_ub_mid!=data.modification_id)){
    if (cb_out()) 
      get_out()<<"*** ERROR: AFTModel::check_center_validity_by_candidate: no such center or candidate available for checking validity"<<std::endl;
    if (cb_out(10)){
      get_out()<<"\n leaving AFTModel::check_center_validity_by_candidate: return value 1"<<std::endl;
    }    
    return 1;
  }  

  if ((aft->get_fun_coeff()==0.)||(model==0)){
    //function is linear in y, don't even check
    cand_minorant_is_below=true;
    if (cb_out(10)){
      get_out()<<"\n leaving AFTModel::check_center_validity_by_candidate: return value 0"<<std::endl;
    }    
    return 0;
  }

  int status; 
  cand_minorant_is_below=true;
  if (aft->argument_changes())
    status=cond_model_trafo()->check_center_validity_by_candidate(cand_minorant_is_below,data.center_pid,data.aft_center);
  else 
    status=cond_model_trafo()->check_center_validity_by_candidate(cand_minorant_is_below,data.center_pid,center_y);
  if (status){
    if (cb_out()) 
      get_out()<<"*** ERROR: AFTModel::check_center_validity_by_candidate: calling check_center_validity_by_candidate returned "<<status<<std::endl;
    if (cb_out(10)){
      get_out()<<"\n leaving AFTModel::check_center_validity_by_candidate: return value 1"<<std::endl;
    }    
    return 1;
  }  
  if ((cand_minorant_is_below==false)&&(cb_out(0))){
    get_out()<<"\n AFTModel::check_center_validity_by_candidate: center_ub lies below cand_minorant"<<std::endl;
  }

  if (cb_out(10)){
    get_out()<<"\n leaving AFTModel::check_center_validity_by_candidate: cand_minorant_is_below="<<cand_minorant_is_below<<", return value 0"<<std::endl;
  } 
   
  return 0;
}


// *****************************************************************************
//                                sumbundle_mode
// *****************************************************************************
 

int AFTModel::sumbundle_mode(SumBundle::Mode& mode,
			      SumBundleHandler* bh,
			      AffineFunctionTransformation* pre_aft)
{
  if (cb_out(10)){
    get_out()<<"\n  entering  AFTModel::sumbundle_mode"<<std::endl;
  }

  int err=0;
  if ((mode==SumBundle::unavailable)||
      (bh==0)||(pre_aft==0)||
      (aft->get_fun_coeff()==0)||(model==0)) {
    if (bundlehandler){
      // there is a local bundlehandler, but we do not need it here
      if (bundlehandler->remove_contributions()){
	if (cb_out()){
	  get_out()<<"*** ERROR: AFTModel::sumbundle_mode(): bundlehandler->normalize_sumbundle() failed"<<std::endl;
	}
	err++;
      }
      int retval=bundlehandler->set_parent_information(0,0,SumBundle::inactive);
      if (retval){
	if (cb_out()){
	  get_out()<<"*** ERROR: AFTModel::sumbundle_mode(): bundlehandler->set_parent_information() failed and returned "<<retval<<std::endl;
	}
	err++;
      }
    }
    mode=SumBundle::unavailable;
  }

  else {
    //parent handler bh is not NULL, there is a model and one aft before, so copy the next sumbundle here if needed

    SumBundle::Mode local_mode=mode;
    if (bundlehandler==0){
      std::map<FunctionTask,Real> factor_map;
      if (bh->handles(ObjectiveFunction))
	factor_map[ObjectiveFunction]=1.;
      if (bh->handles(ConstantPenaltyFunction))
	factor_map[ConstantPenaltyFunction]=1.;
      if (bh->handles(AdaptivePenaltyFunction))
	factor_map[AdaptivePenaltyFunction]=1.;      
      bundlehandler=new SumBundleHandler(data.sumbundle,factor_map);
      bundlehandler->set_cbout(this,0);
      if (mode==SumBundle::child){
	//for add_contribution to work the mode has to be root
	local_mode=SumBundle::root;  
      }
    }

    if (bundlehandler->normalize_sumbundle()){
      if (cb_out()){
	get_out()<<"*** ERROR: AFTModel::sumbundle_mode(): bundlehandler->normalize_sumbundle() failed"<<std::endl;
      }
      err++;
    }

    if (bundlehandler->set_parent_information(bh,pre_aft,local_mode)){
      if (cb_out()){
	get_out()<<"*** ERROR: AFTModel::sumbundle_mode(): SumBundleHandler::initialize_subhandler() failed"<<std::endl;
      }
      err++;
    }
  }

  if (cb_out(10)){
    get_out()<<"\n  leaving  AFTModel::sumbundle_mode with return value "<<err<<std::endl;
  }
  return err;
}


// *****************************************************************************
//                                start_sumaugmodel
// *****************************************************************************

int AFTModel::start_sumaugmodel(QPModelDataPointer& qpblockp,
			     Integer in_id,
			     const Matrix& in_y,
			     const Indexmatrix* indices,
			     SumBundleHandler* bh,
			     SumBundle::Mode mode,
			     AffineFunctionTransformation* pre_aft)
{
  if (cb_out(10)){
    get_out()<<"\n  entering  AFTModel::start_sumaugmodel"<<std::endl;
  }
  
  int err=0;
  SumBundle::Mode this_mode=mode;
  if (sumbundle_mode(this_mode,bh,pre_aft)){
    if (cb_out()){
      get_out()<<"**** ERROR AFTModel::start_sumaugmodel(): sumbundle_mode() failed"<<std::endl;
    }
    err++;
  }
  else {
    
    //==== call start_sumaugmodel of the model

    if ((aft->get_fun_coeff()==0.)||(model==0)){
      //--- there is no model and no contribution to sumbundle
      assert(this_mode==SumBundle::unavailable);
      if (qpblockp.set_model_data(block)){
	delete block;
	block=qpblockp.generate_conemodel_data();
	block->set_cbout(this);	  
      }
      assert(block);
      MinorantBundle dummy_bundle;
      Matrix dummy_mat(0,0,0.);
      block->init(aft->get_constant_minorant(),dummy_bundle,
		  0,Indexmatrix(0,0,Integer(0)),Indexmatrix(0,0,Integer(0)),
		  dummy_mat,dummy_mat,1.,ObjectiveFunction);
      assert(dummy_bundle.size()==0);
      qpblockp.set_model_data(block);
      assert(qpblockp.get_model_data_ptr()==block);
    }

    else {
      
      //--- model needs to be called

      //determine sumbundle arguments
      SumBundleHandler* call_bh=bh;
      SumBundle::Mode call_mode=mode;
      AffineFunctionTransformation* call_aft=aft;
      if (bh!=0){
	if (mode==SumBundle::unavailable){
	  assert(this_mode==SumBundle::unavailable);
	  call_bh=0;
	}
	else if((pre_aft!=0)&&(bundlehandler)){
	  call_bh=bundlehandler;
	}
      }

      // compute the indices
      aft->qp_cost_indices(qp_indices,indices);
      const Indexmatrix* l_ind=&qp_indices;
      if ((qp_indices.coldim()==0)||
	((aft->to_dim()>=0)&&(aft->to_dim()==qp_indices.rowdim()))||
	  ((aft->to_dim()<0)&&(qp_indices.rowdim()==data.dim))
	  ){  
	l_ind=0;
	qp_indices.init(0,0,Integer(0));
      }

      // clean up the old precomputed map
      std::map<MinorantPointer,MinorantPointer>::iterator mapit;
      if (((l_ind==0)&&(old_indices.dim()>0))||((l_ind)&&(!equal(*l_ind,old_indices)))){
	precomputed.clear();
      }
      else {
	mapit=precomputed.begin();
	while(mapit!=precomputed.end()) {
	  if (mapit->first.one_user()){
	    //this is the only copy left
	    std::map<MinorantPointer,MinorantPointer>::iterator delit=mapit;
	    mapit++;
	    precomputed.erase(delit);
	  }
	  else {
	    mapit++;
	  }
	}
      }
      if (l_ind==0)
	old_indices.init(0,0,Integer(0));
      else 
	old_indices=*l_ind;

      // get the transformed point
      const Matrix& local_cand_y=generate_aft_cand(in_id,in_y);

      // call the model
      int retval=cond_model_trafo()->start_sumaugmodel(qpblockp,data.cand_pid,local_cand_y,l_ind,call_bh,call_mode,call_aft);
      if (retval){
	if (cb_out(0)){
	  get_out()<<"\n**** ERROR: AFTModel::start_sumaugmodel(...): start_sumaugmodel failed and returned "<<retval<<std::endl;
	}
	err++;
      }
      else if (qpblockp.get_model_data_ptr()){
	//transform the model
	if (qpblockp.get_model_data_ptr()->push_aft(aft,indices,l_ind,&precomputed)){
	  if (cb_out(0)){
	    get_out()<<"\n**** ERROR: AFTModel::start_sumaugmodel(...): ->push_aft failed for the QP model"<<std::endl;
	  }
	  err++;
	}
      }
      
      data.aggregate.clear();
      data.aggregate_fid=-1;
      data.aggregate_available=false;

    }//end else   call the mode
  } //end else sumbundle_model did not fail
  
  //--- bundlehandler actions
  if (bundlehandler){
    if (this_mode==SumBundle::unavailable){
      //everything should be removed already by now
      assert(data.sumbundle.has_contributions()==false);
      assert(data.sumbundle.has_working_roots()==false);
      bundlehandler->clear_model();
      delete bundlehandler;
      bundlehandler=0;
    }
    else {
      //add all new roots, if there are any; does not hurt to call repeatedly
      if (bundlehandler->add_contributions()){
	if (cb_out(0)){
	  get_out()<<"\n**** ERROR: AFTModel::start_sumaugmodel(...): ->apply_aft failed for the QP model"<<std::endl;
	}
	err++;
      }
    }
  }
      
  assert((bundlehandler==0)||(!data.sumbundle.has_roots()));

  if (cb_out(10)){
    get_out()<<"\n  leaving  AFTModel::start_sumaugmodel with return value "<<err<<std::endl;
  }

  return err;
}



// *****************************************************************************
//                                make_model_aggregate
// *****************************************************************************

int AFTModel::make_model_aggregate(bool& penalty_increased,bool penalty_fixed)  
{
  if (cb_out(10)){
    get_out()<<"\n  entering  AFTModel::make_model_aggregate"<<std::endl;
  }

  penalty_increased=false;
  data.aggregate_available=false;
  data.local_aggregate.clear();
  data.aggregate.clear();

  if ((aft->get_fun_coeff()==0.)||(model==0)){
    //it is stored directly in aft, nothing to do
    if (!data.aggregate_available)
      data.aggregate_id++;
    data.aggregate_available=true;
    if (cb_out(10)){
      get_out()<<"\n  leaving  AFTModel::make_model_aggreagte with return value 0"<<std::endl;
    }
    return 0;
  }

  if (bundlehandler){
    bool local_increased=penalty_increased;
    if (bundlehandler->make_model_aggregate(local_increased,penalty_fixed)){
      if (cb_out(0)){
	get_out()<<"\n**** ERROR:  AFTModel::make_model_aggregate(...): init_make_model_aggregate failed for bundlehandler"<<std::endl;
      }
      if (cb_out(10)){
	get_out()<<"\n  leaving  AFTModel::make_model_aggregate with return value 1"<<std::endl;
      }
      return 1;
    }
    penalty_increased |=local_increased;
  }


  data.aggregate_fid=-1;
    
  int retval=cond_model_trafo()->make_model_aggregate(penalty_increased,penalty_fixed);
  if (retval){
    if (cb_out(0)){
      get_out()<<"\n**** ERROR:  AFTModel::make_model_aggregate(...): make_model_aggregate failed  and returned "<<retval<<std::endl;
    }
    if (cb_out(10)){
      get_out()<<"\n  leaving  AFTModel::make_model_aggregate with return value 1"<<std::endl;
    }
    
    return 1;
  }

  if (penalty_increased) {
    data.center_ub_mid=-1;
  }
  else {
    data.aggregate_id++;
    data.aggregate_available=true;
  }

  if (cb_out(10)){
    get_out()<<"\n  leaving  AFTModel::make_aug_linmodel with return value 0"<<std::endl;
  }
  return 0;
} 
 

// *****************************************************************************
//                              get_model_aggregate
// *****************************************************************************

  int AFTModel::get_model_aggregate(Integer& o_aggr_id,MinorantPointer& aggr, bool all_parts,const AffineFunctionTransformation* in_aft)
{
  if (cb_out(10)){
    get_out()<<"\n  entering  AFTModel::get_model_aggregate"<<std::endl;
  }

  if (!data.aggregate_available) {
    if (cb_out(0)){
      get_out()<<"\n**** WARNING: AFTModel::get_model_aggregate(...) called without an available aggregate"<<std::endl;
    }
    if (cb_out(10)){
      get_out()<<"\n  leaving  AFTModel::get_model_aggregate with return value 1"<<std::endl;
    }
    return 1;
  }

  if (in_aft!=0){
    // treat the case of *this just representing a linear function that is added
    if ((aft->get_fun_coeff()==0.)||(model==0)){
      int retval=in_aft->transform_minorant(aggr,aft->get_constant_minorant(),1,true);
      if (retval){
	if (cb_out()){
	  get_out()<<"\n**** WARNING: AFTModel::get_model_aggregate(...): for in_aft!=0 and aft->get_fun_coeff==0. routine aft->transform_minorant returned "<<retval<<std::endl;
	}
	if (cb_out(10)){
	  get_out()<<"\n  leaving  AFTModel::get_model_aggregate with return value "<<retval<<std::endl;
	}
	return retval;
      }
      o_aggr_id=data.aggregate_id;
      if (cb_out(10)){
	get_out()<<"\n  leaving  AFTModel::get_model_aggregate with return value "<<0<<std::endl;
      }    
      return 0;
    }

    // treat the case that only the local part of the transformed model without the sumbundle part is to be added
    if (!all_parts) {
      if (!data.local_aggregate.valid()){
	data.local_aggregate.clear();
	int retval=cond_model_trafo()->get_model_aggregate(data.aggregate_fid,data.local_aggregate,all_parts,aft);
	if (retval){
	  if (cb_out()){
	    get_out()<<"\n**** WARNING: AFTModel::get_model_aggregate(...): cond_model_trafo()->get_model_aggregate returned "<<retval<<std::endl;
	  }
	  if (cb_out(10)){
	    get_out()<<"\n  leaving  AFTModel::get_model_aggregate with return value "<<retval<<std::endl;
	  }    
	  return retval;
	}
      }
      if (!data.local_aggregate.empty()){
	int retval=in_aft->transform_minorant(aggr,data.local_aggregate,1.,true);
	if (retval){
	  if (cb_out()){
	    get_out()<<"\n**** WARNING: AFTModel::get_model_aggregate(...): for in_aft!=0 routine aft->transform_minorant returned "<<retval<<std::endl;
	  }
	  if (cb_out(10)){
	    get_out()<<"\n  leaving  AFTModel::get_model_aggregate with return value "<<retval<<std::endl;
	  }    
	  return retval;
	}
      }
      o_aggr_id=data.aggregate_id;
      if (cb_out(10)){
	get_out()<<"\n  leaving  AFTModel::get_model_aggregate with return value "<<0<<std::endl;
      }    
      return 0;
    }

    //treat the usual aggregate 
    assert(all_parts==true);
    if (!data.aggregate.valid()){
      int retval=cond_model_trafo()->get_model_aggregate(data.aggregate_fid,data.aggregate,all_parts,aft);
      if (retval){
        if (cb_out()){
          get_out()<<"\n**** WARNING: AFTModel::get_model_aggregate(...): cond_model_trafo()->get_model_aggregate returned "<<retval<<std::endl;
        }
        if (cb_out(10)){
          get_out()<<"\n  leaving  AFTModel::get_model_aggregate with return value "<<retval<<std::endl;
        }    
        return retval;
      }
    } 
    int retval=in_aft->transform_minorant(aggr,data.aggregate,1.,true);
    if (retval){
      if (cb_out()){
	get_out()<<"\n**** WARNING: AFTModel::get_model_aggregate(...): for in_aft!=0 routine aft->transfrom_minorant returned "<<retval<<std::endl;
      }
      if (cb_out(10)){
	get_out()<<"\n  leaving  AFTModel::get_model_aggregate with return value "<<retval<<std::endl;
      }    
      return retval;
    }
    o_aggr_id=data.aggregate_id;
    if (cb_out(10)){
      get_out()<<"\n  leaving  AFTModel::get_model_aggregate with return value "<<0<<std::endl;
    }    
    return 0;
  }

  assert(in_aft==0);
  
  if ((aft->get_fun_coeff()==0.)||(model==0)){
    //it is stored directly in aft, nothing to do
    if (aft->get_constant_minorant().get_minorant(aggr)){
      if (cb_out()){
	get_out()<<"\n**** WARNING: AFTModel::get_model_aggregate(...): for in_aft==0 routine aft->get_constant_minorant().get_minorant() failed "<<std::endl;
      }
      if (cb_out(10)){
	get_out()<<"\n  leaving  AFTModel::get_model_aggregate with return value "<<1<<std::endl;
      }    
      return 1;
    }
    o_aggr_id=data.aggregate_id;

    if (cb_out(10)){
      get_out()<<"\n  leaving  AFTModel::get_model_aggregate with return value "<<0<<std::endl;
    }    
    return 0;
  }

  if ((all_parts)&&(data.aggregate.valid())){
    if(data.aggregate.get_minorant(aggr)){
      if (cb_out()){
	get_out()<<"\n**** WARNING: AFTModel::get_model_aggregate(...): for in_aft==0 routine data.aggregate.get_minorant() failed "<<std::endl;
      }
      if (cb_out(10)){
	get_out()<<"\n  leaving  AFTModel::get_model_aggregate with return value "<<1<<std::endl;
      }    
      return 1;
    }
  }
  else if ((!all_parts)&&(data.local_aggregate.valid())){
    if(data.local_aggregate.get_minorant(aggr)){
      if (cb_out()){
	get_out()<<"\n**** WARNING: AFTModel::get_model_aggregate(...): for in_aft==0 routine data.local_aggregate.get_minorant() failed "<<std::endl;
      }
      if (cb_out(10)){
	get_out()<<"\n  leaving  AFTModel::get_model_aggregate with return value "<<1<<std::endl;
      }    
      return 1;
    }
  }
  else {
    int retval=cond_model_trafo()->get_model_aggregate(data.aggregate_fid,aggr,all_parts,aft);
    if (retval){
      if (cb_out()){
	get_out()<<"\n**** WARNING: AFTModel::get_model_aggregate(...): cond_model_trafo()->get_model_aggregate returned "<<retval<<std::endl;
      }
      if (cb_out(10)){
	get_out()<<"\n  leaving  AFTModel::get_model_aggregate with return value "<<retval<<std::endl;
      }    
      return retval;
    }
  }
  o_aggr_id=data.aggregate_id;

  if (cb_out(10)){
    get_out()<<"\n  leaving  AFTModel::get_model_aggregate with return value "<<0<<std::endl;
  }    
  return 0;
}


// *****************************************************************************
//                                 lb_function
// *****************************************************************************

//returns a *quick* lower bound for the function value at y
//(eg by a previous subgradient)

Real AFTModel::lb_function(Integer y_id,const Matrix& y)
{
  if (cb_out(10)){
    get_out()<<"\n  entering  AFTModel::lb_function"<<std::endl;
  }
  if (data.dim<0)
    data.dim=y.dim();
  assert(data.dim==y.dim());

  if ((aft->get_fun_coeff()==0.)||(model==0)){
    if (cb_out(10)){
      get_out()<<"\n  leaving  AFTModel::lb_function, fun_coeff==0., using linear part"<<std::endl;
    }
    return aft->get_constant_minorant().evaluate(y_id,y);
  }

  if ((data.aggregate_available)&&(data.aggregate.valid())){
    if (cb_out(10)){
      get_out()<<"\n  leaving  AFTModel::lb_function, using aggregate"<<std::endl;
    }
    return data.aggregate.evaluate(y_id,y);
  }

  const Matrix& my_y=generate_aft_cand(y_id,y);

  Real lb=cond_model_trafo()->lb_function(data.cand_pid,my_y);

  if (cb_out(10)){
    get_out()<<"\n  leaving  AFTModel::lb_function"<<std::endl;
  }
  return aft->objective_value(data.cand_offset,lb);
}


// *****************************************************************************
//                            get_function_minorant
// *****************************************************************************


int AFTModel::get_function_minorant(MinorantPointer& minorant,const AffineFunctionTransformation* in_aft)
{
  if (cb_out(10)){
    get_out()<<"\n  entering  AFTModel::get_function_minorant"<<std::endl;
  }  

  if (in_aft!=0){
    if (!cand_minorant.valid()){
      if ((aft->get_fun_coeff()==0.)||(model==0)){
	int retval=in_aft->transform_minorant(minorant,aft->get_constant_minorant(),1.,true);
	if (retval){
	  if (cb_out()){
	    get_out()<<"\n**** WARNING: AFTModel::get_function_minorant(....): for in_aft!=0 and aft->get_fun_coeff==0. routine in_aft->transform_minorant returned "<<retval<<std::endl;
	  }
	  if (cb_out(10)){
	    get_out()<<"\n  leaving  AFTModel::get_function_minorant with return value "<<retval<<std::endl;
	  }
	  return retval;
	}
	if (cb_out(10)){
	  get_out()<<"\n  leaving  AFTModel::get_function_minorant with return value "<<0<<std::endl;
	}    
	return 0;
      }

      assert(aft->get_fun_coeff()!=0);

      cand_minorant.clear();
      int retval=cond_model_trafo()->get_function_minorant(cand_minorant,aft);
      if (retval){
	if (cb_out()){
	  get_out()<<"\n**** WARNING: AFTModel::get_function_minorant(....): cond_model_trafo()->get_function_minorant returned "<<retval<<std::endl;
	}
	if (cb_out(10)){
	  get_out()<<"\n  leaving  AFTModel::get_function_minorant with return value "<<retval<<std::endl;
	}    
	return retval;
      }
    } 
    int retval=in_aft->transform_minorant(minorant,cand_minorant,1.,true);
    if (retval){
      if (cb_out()){
	get_out()<<"\n**** WARNING: AFTModel::get_function_minorant(....): for in_aft!=0 routine aft->transform_minorant returned "<<retval<<std::endl;
      }
      if (cb_out(10)){
	get_out()<<"\n  leaving  AFTModel::get_function_minorant with return value "<<retval<<std::endl;
      }    
      return retval;
    }
    if (cb_out(10)){
      get_out()<<"\n  leaving  AFTModel::get_function_minorant with return value "<<0<<std::endl;
    }    
    return 0;
  }

  assert(in_aft==0);
  
  if ((aft->get_fun_coeff()==0.)||(model==0)){
    //it is stored directly in aft, nothing to do
    if (aft->get_constant_minorant().get_minorant(minorant)){
      if (cb_out()){
	get_out()<<"\n**** WARNING: AFTModel::get_function_minorant(....): for in_aft==0 routine aft->get_constant_minorant().get_minorant() failed "<<std::endl;
      }
      if (cb_out(10)){
	get_out()<<"\n  leaving  AFTModel::get_function_minorant with return value "<<1<<std::endl;
      }    
      return 1;
    }
    if (cb_out(10)){
      get_out()<<"\n  leaving  AFTModel::get_function_minorant with return value "<<0<<std::endl;
    }    
    return 0;
  }

  if (cand_minorant.valid()){
    if (cand_minorant.get_minorant(minorant)){
      if (cb_out()){
	get_out()<<"\n**** WARNING: AFTModel::get_function_minorant(....): for in_aft==0 routine cand_minorant.get_minorant() failed "<<std::endl;
      }
      if (cb_out(10)){
	get_out()<<"\n  leaving  AFTModel::get_function_minorant with return value "<<1<<std::endl;
      }    
      return 1;
    }
  }
  else {
    int retval=cond_model_trafo()->get_function_minorant(minorant,aft);
    if (retval){
      if (cb_out()){
	get_out()<<"\n**** WARNING: AFTModel::get_function_minorant(....): cond_model_trafo()->get_function_minorant returned "<<retval<<std::endl;
      }
      if (cb_out(10)){
	get_out()<<"\n  leaving  AFTModel::get_function_minorant with return value "<<retval<<std::endl;
      }    
      return retval;
    }
  }

  if (cb_out(10)){
    get_out()<<"\n  leaving  AFTModel::get_function_minorant with return value "<<0<<std::endl;
  }

  return 0;
}


// *****************************************************************************
//                            get_center_minorant
// *****************************************************************************


int AFTModel::get_center_minorant(MinorantPointer& minorant,const AffineFunctionTransformation* in_aft)
{
  if (cb_out(10)){
    get_out()<<"\n  entering  AFTModel::get_center_minorant"<<std::endl;
  }  

  if (in_aft!=0){
    if (!center_minorant.valid()){
      if ((aft->get_fun_coeff()==0.)||(model==0)){
	int retval=in_aft->transform_minorant(minorant,aft->get_constant_minorant(),1.,true);
	if (retval){
	  if (cb_out()){
	    get_out()<<"\n**** WARNING: AFTModel::get_center_minorant(....): for in_aft!=0 and aft->get_fun_coeff==0. routine in_aft->transform_minorant returned "<<retval<<std::endl;
	  }
	  if (cb_out(10)){
	    get_out()<<"\n  leaving  AFTModel::get_center_minorant with return value "<<retval<<std::endl;
	  }
	  return retval;
	}
	if (cb_out(10)){
	  get_out()<<"\n  leaving  AFTModel::get_center_minorant with return value "<<0<<std::endl;
	}    
	return 0;
      }

      assert(aft->get_fun_coeff()!=0);

      center_minorant.clear();
      int retval=cond_model_trafo()->get_center_minorant(center_minorant,aft);
      if (retval){
	if (cb_out()){
	  get_out()<<"\n**** WARNING: AFTModel::get_center_minorant(....): cond_model_trafo()->get_center_minorant returned "<<retval<<std::endl;
	}
	if (cb_out(10)){
	  get_out()<<"\n  leaving  AFTModel::get_center_minorant with return value "<<retval<<std::endl;
	}    
	return retval;
      }
    } 
    int retval=in_aft->transform_minorant(minorant,center_minorant,1.,true);
    if (retval){
      if (cb_out()){
	get_out()<<"\n**** WARNING: AFTModel::get_center_minorant(....): for in_aft!=0 routine aft->transform_minorant returned "<<retval<<std::endl;
      }
      if (cb_out(10)){
	get_out()<<"\n  leaving  AFTModel::get_center_minorant with return value "<<retval<<std::endl;
      }    
      return retval;
    }
    if (cb_out(10)){
      get_out()<<"\n  leaving  AFTModel::get_center_minorant with return value "<<0<<std::endl;
    }    
    return 0;
  }

  assert(in_aft==0);
  
  if ((aft->get_fun_coeff()==0.)||(model==0)){
    //it is stored directly in aft, nothing to do
    if (aft->get_constant_minorant().get_minorant(minorant)){
      if (cb_out()){
	get_out()<<"\n**** WARNING: AFTModel::get_center_minorant(....): for in_aft==0 routine aft->get_constant_minorant().get_minorant() failed "<<std::endl;
      }
      if (cb_out(10)){
	get_out()<<"\n  leaving  AFTModel::get_center_minorant with return value "<<1<<std::endl;
      }    
      return 1;
    }
    if (cb_out(10)){
      get_out()<<"\n  leaving  AFTModel::get_center_minorant with return value "<<0<<std::endl;
    }    
    return 0;
  }

  if (center_minorant.valid()){
    if (center_minorant.get_minorant(minorant)){
      if (cb_out()){
	get_out()<<"\n**** WARNING: AFTModel::get_center_minorant(....): for in_aft==0 routine center_minorant.get_minorant() failed "<<std::endl;
      }
      if (cb_out(10)){
	get_out()<<"\n  leaving  AFTModel::get_center_minorant with return value "<<1<<std::endl;
      }    
      return 1;
    }
  }
  else {
    int retval=cond_model_trafo()->get_center_minorant(minorant,aft);
    if (retval){
      if (cb_out()){
	get_out()<<"\n**** WARNING: AFTModel::get_center_minorant(....): cond_model_trafo()->get_center_minorant returned "<<retval<<std::endl;
      }
      if (cb_out(10)){
	get_out()<<"\n  leaving  AFTModel::get_center_minorant with return value "<<retval<<std::endl;
      }    
      return retval;
    }
  }

  if (cb_out(10)){
    get_out()<<"\n  leaving  AFTModel::get_center_minorant with return value "<<0<<std::endl;
  }

  return 0;
}



// *****************************************************************************
//                                adjust_multiplier
// *****************************************************************************

int AFTModel::adjust_multiplier(bool& values_may_have_changed)  
{
  if (cb_out(10)){
    get_out()<<"\n  entering  AFTModel::adjust_multplier"<<std::endl;
  }

  if ((aft->get_fun_coeff()==0)||(model==0)){
    if (cb_out(10)){
      get_out()<<"\n  leaving  AFTModel::adjust_multplier with return value 0"<<std::endl;
    }
    return 0;
  }

  if (bundlehandler){
    bool maybe_changes=values_may_have_changed;
    if (bundlehandler->adjust_multiplier(maybe_changes)){
      if (cb_out()){
	get_out()<<"**** ERROR AFTModel::adjust_multiplier(.): bundlehandler->adjust_multilier(.) failed"<<std::endl;
      }
      if (cb_out(10)){
	get_out()<<"\n  leaving  AFTModel::adjust_multiplier(.) with return value 1"<<std::endl;
      }
      return 1;
    }
    values_may_have_changed |= maybe_changes;    
  }


  bool maybe_changes=false;
  int retval=cond_model_trafo()->adjust_multiplier(maybe_changes);
  if (retval){
    if (cb_out(0)){
      get_out()<<"\n**** WARNING: AFTModel::adjust_multplier(): adjust_multiplier failed and returned "<<retval<<std::endl;
    }
    maybe_changes=true;
  }

  if (maybe_changes){
    values_may_have_changed=true;
    clear_model();
  }

  if (cb_out(10)){
    get_out()<<"\n  leaving  AFTModel::adjust_multplier with return value "<<retval<<std::endl;
  }
  return 0;
}

// *****************************************************************************
//                                clear_model
// *****************************************************************************

void AFTModel::clear_model(bool discard_minorants_only)  
{
  if (cb_out(10)){
    get_out()<<"\n  entering  AFTModel::clear_model"<<std::endl;
  }

  if (bundlehandler){
    SumBundle::Mode inact=SumBundle::inactive;
    if (cond_model_trafo()->sumbundle_mode(inact)){
      if (cb_out()) 
	get_out()<<"**** ERROR: AFTModel::clear_model(): sumbundle_mode failed"<<std::endl;
    } 
    bundlehandler->clear_model();
    assert(!data.sumbundle.active());
    delete bundlehandler;
    bundlehandler=0;
  }

  data.clear_model(discard_minorants_only);
  center_minorant.clear();
  cand_minorant.clear();
  precomputed.clear();

  if (cb_out(10)){
    get_out()<<"\n  leaving  AFTModel::clear_model"<<std::endl;
  }

}

// *****************************************************************************
//                                clear_aggregates
// *****************************************************************************

void AFTModel::clear_aggregates()  
{
  if (cb_out(10)){
    get_out()<<"\n  entering  AFTModel::clear_aggregates"<<std::endl;
  }

  if (data.sumbundle.has_bundle_data()){
    if (bundlehandler){
      bundlehandler->clear_aggregates();
    }
    else {
      std::map<FunctionTask,Real> fmap;
      bundlehandler=new SumBundleHandler(data.sumbundle,fmap);
      bundlehandler->set_cbout(this,0);
      bundlehandler->clear_aggregates();
      delete bundlehandler;
      bundlehandler=0;
    }
  }    

  data.clear_aggregates();

  if (cb_out(10)){
    get_out()<<"\n  leaving  AFTModel::clear_aggregates"<<std::endl;
  }
}

// *****************************************************************************
//                               output_bundle_data
// *****************************************************************************

/*
std::ostream& AFTModel::output_bundle_data(std::ostream& out) const
{ 
  out<<"\n%(begin AFTModel)\n";
  assert(aft);
    
  aft->output_aft_data(out);
  if (model){
      out<<"\n%(begin AFTSubModel)\n";
      model->output_bundle_data(out);
      out<<"\n%(end AFTSubModel)\n";
  }
  
  out<<"\n%(end AFTModel)\n";
  return out;
}
*/

}

