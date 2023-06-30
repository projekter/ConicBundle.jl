/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/SumBlockModel.cxx
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
#include "SumBlockModel.hxx"
#include "AFTModel.hxx"
#include "SumBundleParameters.hxx"

#include "BundleDiagonalTrustRegionProx.hxx"
 
using namespace CH_Matrix_Classes;

namespace ConicBundle {

// *****************************************************************************
//                              ~SumBlockModel
// *****************************************************************************

SumBlockModel::~SumBlockModel()
{
  SumBlockModel::clear();
  delete sumbundle_parameters;
}

// *****************************************************************************
//                             initialize_aft
// *****************************************************************************

int SumBlockModel::initialize_aft(AffineFunctionTransformation* aft)
{
  if ((aftmodel==0)&&(aft==0)){
    return 0;
  }

  Integer old_modification_id=get_data()->modification_id;
  Integer aft_mod_id=-1;
  bool argument_changes=false;

  if (bundlehandler)
    bundlehandler->remove_contributions();

  if (aftmodel!=0){
    aft_mod_id=aftmodel->get_data()->modification_id;
    argument_changes=aftmodel->get_aft()->argument_changes();
    if (aft) { 
      aftmodel->clear(aft,aft_mod_id+1);
      argument_changes = argument_changes || aft->argument_changes();
    }
    else {
      delete aftmodel;
      aftmodel=0;
    }
  }
  else {
    aftmodel=new AFTModel(this,aft,aft_mod_id+1,true);
    argument_changes = aft->argument_changes();
  }
   
  if (aftmodel==0){
    //to the outside world it must become visible that the function changed
    Integer next_modification_id=max(old_modification_id,aft_mod_id)+1;
    if (!argument_changes){
      //the points, where the function was evaluated are not affected
      if (get_data()->center_ub_mid==get_data()->modification_id)
	get_data()->center_ub_mid=next_modification_id;
      if (get_data()->cand_ub_mid==get_data()->modification_id)
	get_data()->cand_ub_mid=next_modification_id;
    }
    get_data()->modification_id=next_modification_id;
  }
  
  if (argument_changes){
    //make sure that the evaluations are marked as not intialized
    get_data()->center_ub_mid=-1;
    get_data()->center_id=-1;
    get_data()->cand_ub_mid=-1;
    get_data()->cand_id=-1;
  }

  //the aggregate and the model stay valid
  return 0;
}

// *****************************************************************************
//                             call_aftmodel_first()
// *****************************************************************************

bool SumBlockModel::call_aftmodel_first(const FunObjModMap& funmdfmap)
{
  if (aftmodel) 
    return false;
  FunObjModMap::const_iterator mapit=funmdfmap.find(get_oracle_object());
  if (mapit==funmdfmap.end())
    return false;
  if (mapit->second.get_aft_modification()==0)
    return false;
  if (mapit->second.get_aft_modification()->no_modification())
    return false;
  AffineFunctionTransformation* aftp=new AffineFunctionTransformation;
  initialize_aft(aftp);

  return true;
}
    

// *****************************************************************************
//                              SumBlockModel
// *****************************************************************************

  SumBlockModel::SumBlockModel(CBout* cb,int cbinc):
    BundleModel(cb,cbinc),aftmodel(0),vm_selection(0),bundlehandler(0),sumbundle_parameters(0)
{  
  SumBlockModel::clear();
}

// *****************************************************************************
//                                  clear()
// *****************************************************************************

void SumBlockModel::clear()
{
  delete aftmodel; 
  aftmodel=0;
  
  clock.start();
  evalmodel_time=0;
  updatemodel_time=0;
  eval_time=0;
  preeval_time=0;
  posteval_time=0;
  metric_time=0;

  delete vm_selection;
  vm_selection=0;
  
  delete bundlehandler;
  bundlehandler= 0; //new SumBundleHandler(data.sumbundle);

  delete sumbundle_parameters;
  sumbundle_parameters=new SumBundleParameters;
}



// *****************************************************************************
//                             synchronize_ids
// *****************************************************************************


  int SumBlockModel::synchronize_ids(Integer& new_center_ub_fid, 
			      Integer new_center_id,
			      Integer old_center_id,
			      Integer& new_cand_ub_fid, 
			      Integer new_cand_id,
			      Integer old_cand_id,
			      Integer& new_aggregate_id)
 {
  if (cb_out(10)){
    get_out()<<"\n  entering  SumBlockModel::synchronize_ids"<<std::endl;
  }

  if (get_data()->synchronize_ids(new_center_ub_fid,new_center_id,old_center_id,new_cand_ub_fid,new_cand_id,old_cand_id,new_aggregate_id,0)){
    if (cb_out(10)){
      get_out()<<"\n  leaving  SumBlockModel::synchronize_ids with return value 1"<<std::endl;
    }
    return 1;
  }
  
  if (cb_out(10)){
    get_out()<<"\n  leaving  SumBlockModel::synchronize_ids with return value 0"<<std::endl;
  }
  return 0;
 }


// *****************************************************************************
//                             center_modified
// *****************************************************************************


  bool SumBlockModel::center_modified(Integer& center_ub_fid,Integer check_center_id)
 {
  if (cb_out(10)){
    get_out()<<"\n  entering  SumBlockModel::center_modified"<<std::endl;
  }

  if (get_data()->center_modified(center_ub_fid,check_center_id)){
    if (cb_out(10)){
      get_out()<<"\n  leaving  SumBlockModel::center_modified with return value true"<<std::endl;
    }
    return true;
  }
  
  if (cb_out(10)){
    get_out()<<"\n  leaving  SumBlockModel::center_modified with return value false"<<std::endl;
  }
  return false;
 }



// *****************************************************************************
//                              start_augmodel
// *****************************************************************************

int SumBlockModel::start_augmodel(QPModelDataPointer& blockp,
				  Integer cand_id,
				  const Matrix& cand_y,
				  const Indexmatrix* indices)
{
  return start_sumaugmodel(blockp,cand_id,cand_y,indices);
}

// *****************************************************************************
//                              get_model_aggregate
// *****************************************************************************

int SumBlockModel::get_model_aggregate(Integer& o_aggr_id,MinorantPointer& o_aggr)
{
  return get_model_aggregate(o_aggr_id,o_aggr,true);
}

// *****************************************************************************
//                              get_model_aggregate
// *****************************************************************************

  int SumBlockModel::get_model_aggregate(Integer& o_aggr_id,MinorantPointer& o_aggr, bool all_parts, const AffineFunctionTransformation* aft)
{
  if (cb_out(10)){
    get_out()<<"\n  entering  SumBlockModel::get_model_aggregate"<<std::endl;
  }
  
  if (!get_data()->aggregate.valid()) {
    if (cb_out(0)){
      get_out()<<"\n**** WARNING: SumBlockModel::get_model_aggregate(...) called without an available aggregate"<<std::endl;
    }
    if (cb_out(10)){
      get_out()<<"\n  leaving  SumBlockModel::get_model_aggregate with return value 1"<<std::endl;
    }
    return 1;
  }
  
  if (all_parts){
    if (aft) {
      if (aft->transform_minorant(o_aggr,get_data()->aggregate,1.,true)){
	if (cb_out()) 
	  get_out()<<"**** WARNING: SumBlockModel::get_model_aggregate(...): collecting minorant information failed"<<std::endl;
	if (cb_out(10)){
	  get_out()<<"\n  leaving  SumBlockModel::get_model_aggregate with return value "<<1<<std::endl;
	}    
	return 1;
      }
    }
    else {
      if (get_data()->aggregate.get_minorant(o_aggr)){
	if (cb_out()) 
	  get_out()<<"**** WARNING: SumBlockModel::get_model_aggregate(...): aggregate.get_minorant failed"<<std::endl;
	if (cb_out(10)){
	  get_out()<<"\n  leaving  SumBlockModel::get_model_aggregate with return value "<<1<<std::endl;
	}    
	return 1;
      }
    }
  }
  else { //get only the local part
    if (!get_data()->local_aggregate.empty()){
      if (aft) {
	if (aft->transform_minorant(o_aggr,get_data()->local_aggregate,1.,true)){
	  if (cb_out()) 
	    get_out()<<"**** WARNING: SumBlockModel::get_model_aggregate(...): collecting minorant information failed"<<std::endl;
	  if (cb_out(10)){
	    get_out()<<"\n  leaving  SumBlockModel::get_model_aggregate with return value "<<1<<std::endl;
	  }    
	  return 1;
	}
      }
      else {
	if (get_data()->local_aggregate.get_minorant(o_aggr)){
	  if (cb_out()) 
	    get_out()<<"**** WARNING: SumBlockModel::get_model_aggregate(...): aggregate.get_minorant failed"<<std::endl;
	  if (cb_out(10)){
	    get_out()<<"\n  leaving  SumBlockModel::get_model_aggregate with return value "<<1<<std::endl;
	  }    
	  return 1;
	}
      }
    }
  }

  o_aggr_id=get_data()->aggregate_id;

  if (cb_out(10)){
    get_out()<<"\n  leaving  SumBlockModel::get_model_aggregate with return value "<<0<<std::endl;
  }    
  return 0;
}

// *****************************************************************************
//                         apply_modification
// *****************************************************************************

int SumBlockModel::apply_modification(bool& no_changes,
				    const GroundsetModification& gsmdf,
				    const FunObjModMap& funmdfmap,
				    Integer new_center_id,
				    const Matrix& new_center,
				    Integer old_center_id,
				    const Matrix& old_center)
{
  if (cb_out(10)){
    get_out()<<"\n  entering  SumBlockModel::apply_modifiction"<<std::endl;
  }

  if (call_aftmodel_first(funmdfmap)){
    int retval=sbm_transform()->apply_modification(no_changes,
					       gsmdf, funmdfmap,
					       new_center_id,new_center,
					       old_center_id,old_center);
    if (cb_out(10)){
      get_out()<<"\n  leaving  SumBlockModel::apply_modification with return value "<<retval<<std::endl;
    }
    return retval;
  }

  assert ((get_data()->dim<0)||(old_center_id<0)||(get_data()->dim==old_center.dim()));
  assert((get_data()->dim<0)||(get_data()->dim==gsmdf.old_vardim()));
  assert((new_center_id<0)==(old_center_id<0));

  // ---- check whether the center is even considered preservable and update its id if so
  if ((old_center_id<0)||(old_center_id!=get_data()->center_id)){
    if ((old_center_id>=0)&&(get_data()->center_ub_mid>=0)&&(cb_out()))
      get_out()<<"**** WARNING SumBlockModel::apply_modification(.......): old_center_id="<<old_center_id<<" does not match stored center_id="<<get_data()->center_id<<", discarding center"<<std::endl;
    get_data()->center_ub_mid=-1;
  }  
  
  // ---- ask the PSCOracle for consequences of its changes (maybe by changing itself at the same time)

  int err=0;
  bool discard_objective_in_center=true;
  bool discard_model=true;
  bool discard_aggregate=true;
  MinorantExtender* mex=0;
  const OracleModification* omdf=0;
  const FunObjModMap::const_iterator mapit=funmdfmap.find(get_oracle_object());
  if (mapit!=funmdfmap.end())
    omdf=mapit->second.get_oracle_modification();
  
  if (omdf==0){
    if (gsmdf.no_modification()){
      discard_objective_in_center=false;
      discard_model=false;
      discard_aggregate=false;
    }
    else {
      omdf=&gsmdf;
    }
  }
  else {
    //check consistency of modification information
    if (gsmdf.get_old_vardim()!=omdf->get_old_vardim()){
      if (cb_out())
	get_out()<<"**** ERROR SumBlockModel::apply_modification(.......): old_vardim="<<gsmdf.get_old_vardim()<<" of the ground set does not match old_vardim="<<omdf->get_old_vardim()<<" of the OracleModification"<<std::endl;
      err++;
    }
    if (gsmdf.get_new_vardim()!=omdf->get_new_vardim()){
      if (cb_out())
	get_out()<<"**** ERROR SumBlockModel::apply_modification(.......): new_vardim="<<gsmdf.get_new_vardim()<<" of the ground set does not match new_vardim="<<omdf->get_new_vardim()<<" of the OracleModification"<<std::endl;
      err++;
    }
    if (gsmdf.get_appended_vardim()!=omdf->get_appended_vardim()){
      if (cb_out())
	get_out()<<"**** ERROR SumBlockModel::apply_modification(.......): append_vardim="<<gsmdf.get_appended_vardim()<<" of the ground set does not match append_vardim="<<omdf->get_appended_vardim()<<" of the OracleModification"<<std::endl;
      err++;
    }
    if (gsmdf.get_map_to_old_variables()||omdf->get_map_to_old_variables()){
      if ((gsmdf.get_map_to_old_variables()==0)
	  ||(omdf->get_map_to_old_variables()==0)
	  ||(gsmdf.get_new_vardim()!=omdf->get_new_vardim())
	  ||(sum(*(gsmdf.map_to_old_variables())-Indexmatrix(omdf->get_new_vardim(),1,gsmdf.get_map_to_old_variables()))!=0)){
	if (cb_out())
	  get_out()<<"**** ERROR SumBlockModel::apply_modification(.......): map_to_old_variables of the ground set does not match map_to_old_variables of the OracleModification"<<std::endl;
	err++;
      }
    }
  }

  if (omdf){
    get_data()->modification_id++;
    int retval=0;
    if (get_data()->center_ub_mid>=0){
      retval=get_oracle_object()->apply_modification(*omdf,&new_center,&old_center,
						     discard_objective_in_center,
						     discard_model,discard_aggregate,
						     mex);
    }
    else {
      retval=get_oracle_object()->apply_modification(*omdf,0,0,
						     discard_objective_in_center,
						     discard_model,discard_aggregate,
						     mex);
    }
    if (retval) {
      if (cb_out()){
	get_out()<<"**** ERROR SumBlockModel::apply_modification(.......): oracle->apply_modfication(.......) failed for MatrixPSCOracle and returned "<<retval<<std::endl;
      }
      err++;
    }
  }
  
  // ------ if the contribution changes, remove it first
  if ((!no_changes)||(discard_model)||(discard_aggregate)||(!gsmdf.no_modification())){
    if (bundlehandler){
      if (bundlehandler->remove_contributions()){
	if (cb_out()){
	  get_out()<<"**** ERROR SumBlockModel::apply_modification(.......): sumbundle->remove_contributions failed"<<std::endl;
	}
	err++;
      }
    }
  }
  
  
  // ----- execute the changes -----
  no_changes=true;   //will be set to false if externally some data needs recomputation
  bool discard_minorants_only=false;

  if (discard_objective_in_center){
      no_changes=false;
    get_data()->center_ub_mid=-1;
  }
  
  if (!discard_model){
    if (discard_aggregate){
      no_changes=false;
      clear_aggregates();  
    }

    if (!gsmdf.no_modification()){
      if ((gsmdf.appended_vardim()>0)&&(mex==0)){
	discard_model=true;
	discard_minorants_only=true;
      }
      else if (get_data()->apply_modification(gsmdf,mex)){
	if (cb_out()){
	  get_out()<<"**** ERROR SumBlockModel::apply_modification(.......): modifying/extending the minorants failed"<<std::endl;
	}
	err++;
	discard_model=true;
      }
    }
  }
  delete mex; mex=0;

  if (discard_model){
    no_changes=false;
    clear_model(discard_minorants_only);
  }
      
  get_data()->dim=gsmdf.new_vardim();
  if (get_data()->center_ub_mid>=0) {
    //value in the center is preserved
    get_data()->center_ub_mid=get_data()->modification_id;
    get_data()->center_id=new_center_id;
  }

  if (cb_out(10)){
    get_out()<<"\n  leaving  SumBlockModel::apply_modification with return value "<<err<<std::endl;
  }
  
  return err;
}


// *****************************************************************************
//                                add_variable_metric
// *****************************************************************************

int SumBlockModel::add_variable_metric(VariableMetric& H,
				       Integer y_id,
				       const Matrix& y,
				       bool descent_step,
				       Real weightu,
				       Real model_maxviol,
				       const Indexmatrix* indices)
{
  if (cb_out(10)){
    get_out()<<"\n  entering  SumBlockModel::add_variable_metric"<<std::endl;
  }

  CH_Tools::Microseconds start_metric=clock.time();

  
  int err=0;
  if (!H.employ_variable_metric()){
    get_data()->set_lowrankH().init(0,0,0.);
    get_data()->set_diagH().init(0,0,0.);
    get_data()->set_denseH().init(0,0.);
  }
  else {
    assert(get_data());
    bool sumbundle_active=((bundlehandler)&&(get_data()->sumbundle.has_working_roots()));
    if ((!H.get_use_local_metric())||(!sumbundle_active)){
      VariableMetricSelection* vms=vm_selection;
      if ((vm_selection==0)
	  ||((!H.get_use_local_metric())&&(H.get_variable_metric_selection()))
	  ){
	vms=H.get_variable_metric_selection();
      }
      if (vms){
	if (vms->add_variable_metric(H,
				     y_id,y,descent_step,
				     weightu,model_maxviol,indices,
				     *get_data())){
	  err++;
	  if (cb_out())
	    get_out()<<"**** WARNING: SumBlockModel::add_variable_metric(): vms->add_variable_metric(.........) failed"<<std::endl;
	}
      }
    }
    else {
      // use local metric and sumbundle is active
      assert(get_data()->sumbundle.active());
      assert(bundlehandler);
      //never use the history here
      get_data()->set_lowrankH().init(0,0,0.);
      get_data()->set_diagH().init(0,0,0.);
      get_data()->set_denseH().init(0,0.);
      if (bundlehandler->add_variable_metric(H,y_id,y,descent_step,weightu,model_maxviol,indices)){
	err++;
	if (cb_out())
	  get_out()<<"**** WARNING: SumBlockModel::add_variable_metric(): bundlehandler->add_variable_metric() failed"<<std::endl;	
      }
    }    
  }

  metric_time+=clock.time()-start_metric;

  if (cb_out(10)){
    get_out()<<"\n  leaving SumBlockModel::add_variable_metric with return value "<<err<<std::endl;
  }
  return err;
}

// *****************************************************************************
//                               get_function_minorant()
// *****************************************************************************

  int SumBlockModel::get_function_minorant(Integer& fmid,MinorantPointer& mnrt)
{
  if (cb_out(10)){
    get_out()<<"\n  entering  SumBlockModel::get_function_minorant()"<<std::endl;
  }
  int err=0;
  fmid=get_data()->get_cand_modification_id();
  if ((fmid<0)||(get_data()->get_modification_id()!=fmid)){
    if (cb_out())
      get_out()<<"\n **** WARNING: SumBlockModel::get_function_minorant(): no minorant is availalbe that is up to date"<<std::endl;
    mnrt.clear();
    fmid=-1;
    err++;
  }
  else {
    err=get_function_minorant(mnrt,0);
    if (err) {
      if (cb_out())
	get_out()<<"\n **** WARNING: SumBlockModel::get_function_minorant(): call to get_function_minorant(mnrt,0) returned "<<err<<std::endl;
      mnrt.clear();
    }
  }
  if (cb_out(10)){
    get_out()<<"\n leaving  SumBlockModel::get_function_minorant() with return value "<<err<<std::endl;
  }
  return err;
}

// *****************************************************************************
//                               get_center_minorant()
// *****************************************************************************

  int SumBlockModel::get_center_minorant(Integer& fmid,MinorantPointer& mnrt)
{
  if (cb_out(10)){
    get_out()<<"\n  entering  SumBlockModel::get_center_minorant()"<<std::endl;
  }
  int err=0;
  fmid=get_data()->get_center_modification_id();
  if ((fmid<0)||(get_data()->get_modification_id()!=fmid)){
    if (cb_out())
      get_out()<<"\n **** WARNING: SumBlockModel::get_center_minorant(): no minorant is availalbe that is up to date"<<std::endl;
    mnrt.clear();
    fmid=-1;
    err++;
  }
  else {
    err=get_center_minorant(mnrt,0);
    if (err) {
      if (cb_out())
	get_out()<<"\n **** WARNING: SumBlockModel::get_center_minorant(): call to get_center_minorant(mnrt,0) returned "<<err<<std::endl;
      mnrt.clear();
    }
  }
  if (cb_out(10)){
    get_out()<<"\n leaving  SumBlockModel::get_center_minorant() with return value "<<err<<std::endl;
  }
  return err;
}




// *****************************************************************************
//                           model_aggregate_modified
// *****************************************************************************


bool SumBlockModel::model_aggregate_modified(Integer last_aggr_id)
 {
  if (cb_out(10)){
    get_out()<<"\n  entering  SumBlockModel::model_aggregate_modified"<<std::endl;
  }

  if (get_data()->model_aggregate_modified(last_aggr_id)) {
    if (cb_out(10)){
      get_out()<<"\n  leaving  SumBlockModel::model_aggregate_modified with return value true"<<std::endl;
    }
    return true;
  }
  
  if (cb_out(10)){
    get_out()<<"\n  leaving  SumBlockModel::model_aggregate_modified with return value false"<<std::endl;
  }
  return false;
 }



// *****************************************************************************
//                                 lb_function
// *****************************************************************************

//returns a *quick* lower bound for the function value at y
//(eg by a previous subgradient)

  Real SumBlockModel::lb_function(Integer y_id,const Matrix& y)
{
  if (cb_out(10)){
    get_out()<<"\n  entering  SumBlockModel::lb_function"<<std::endl;
  }

  if (get_data()->aggregate.valid()){
    if (cb_out(10)){
      get_out()<<"\n  leaving  SumBlockModel::lb_function, using aggregate"<<std::endl;
    }
    return get_data()->aggregate.evaluate(y_id,y);
  }

  if (cb_out(10)){
    get_out()<<"\n  leaving  SumBlockModel::lb_function"<<std::endl;
  }
  return CB_minus_infinity;
}

// *****************************************************************************
//                              transform
// *****************************************************************************

BundleModel* SumBlockModel::transform()
{
  return (aftmodel)? aftmodel : this;
}

// *****************************************************************************
//                              sbm_transform
// *****************************************************************************

SumBlockModel* SumBlockModel::sbm_transform()
{
  return (aftmodel)? aftmodel : this;
}

// *****************************************************************************
//                         set_bundle_parameters      
// *****************************************************************************

int SumBlockModel::set_sumbundle_parameters(const BundleParameters& bp)
{
  const SumBundleParametersObject* spo=dynamic_cast<const SumBundleParametersObject*>(&bp);
  int retval=0;
  if (spo){
    delete sumbundle_parameters;
    sumbundle_parameters=dynamic_cast<SumBundleParametersObject*>(spo->clone_BundleParameters());
    assert(sumbundle_parameters);
  }
  else {
    assert(sumbundle_parameters);
    retval=sumbundle_parameters->init(bp);
  }
  return retval;
}

// *****************************************************************************
//                              print_statistics
// *****************************************************************************

std::ostream& SumBlockModel::print_statistics(std::ostream& out) const
{
  out<<" preeval "<<get_preeval_time();
  out<<" eval "<<get_eval_time();
  out<<" posteval "<<get_posteval_time();
  out<<" model "<<get_evalmodel_time();
  out<<" update "<<get_updatemodel_time();
  out<<" scale "<<get_metric_time()<<"\n";

  return out;
}
  

}

