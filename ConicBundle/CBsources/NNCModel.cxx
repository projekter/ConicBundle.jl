/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/NNCModel.cxx
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


#include <fstream>
#include "mymath.hxx"
#include "NNCModel.hxx"
#include "NNCModelParameters.hxx"
#include "BundleIdProx.hxx"

using namespace CH_Tools;
using namespace CH_Matrix_Classes;

namespace ConicBundle {

// *****************************************************************************
//                            NNCModel
// *****************************************************************************

void NNCModel::clear()
{
  ConeModel::clear();
  data.clear();
  delete model_selection;
  model_selection=new NNCModelParameters(this);
  
  ret_code=0;
  nr_eval=0;
  delete block;
  block=0;
}

// *****************************************************************************
//                            NNCModel
// *****************************************************************************

  NNCModel::NNCModel(MatrixFunctionOracle* fo,
		     Real fun_factor,
		     FunctionTask fun_task,
		     CBout* cb,int cbinc):
    CBout(cb,cbinc),ConeModel(cb,cbinc),oracle(fo),data(fun_factor,fun_task),model_selection(0),block(0)
{
  assert(fo);
  clear();
}

// *****************************************************************************
//                            NNCModel
// *****************************************************************************

NNCModel::~NNCModel()
{
  delete model_selection;
  delete block;
}

// *****************************************************************************
//                            NNCModel
// *****************************************************************************

int NNCModel::adapt_bundle_to_parameters(Integer cand_id,const Matrix& cand_y)
{
  data.outside_bundle=false;
  assert(data.bundlecoeff.rowdim()==Integer(data.bundle.size()));
  Integer kmax=(model_selection->get_max_model_size()<2)?Integer(data.bundle.size()):Integer(model_selection->get_max_model_size());
  assert(kmax>=2);
  if (Integer(data.bundle.size())<=kmax) //nothing to do
    return 0;
  
  MinorantBundle bundle=data.bundle;
  Matrix coeff=data.bundlecoeff;
  data.bundlecoeff.init(kmax,1,0.);
  data.bundle.resize(unsigned(kmax));
  
  // first find the best subgradient and ensure it is contained
  Matrix objval(Integer(bundle.size()),1,0.);
  Integer bestind=-1;
  Real bestval=min_Real;
  for(Integer i=0;i<objval.rowdim();i++){
    Real d=bundle[unsigned(i)].evaluate(cand_id,cand_y);
    objval(i)=d;
    if (d>bestval){
      bestval=d;
      bestind=i;
    }
  }
  Integer k=0;
  if (bestind>=0){
    data.bundle[unsigned(k)]=bundle[unsigned(bestind)];
    data.bundlecoeff(k)=coeff(bestind);
    k++;
    objval(bestind)=max_Real; //mark as treated
  }
            
  // next collect all with possitive coefficients or aggregate them if too many
  Indexmatrix sind;
  sortindex(coeff,sind,false);
  MinorantBundle aggrbundle;
  Matrix aggrcoeff;
  Integer i=0;
  while((i<coeff.rowdim())&&(coeff(i)>eps_Real*data.get_function_factor())){
    Integer bi=sind(i++);
    if (bi==bestind)
      continue;
    objval(bi)=max_Real; //mark as treated
    const MinorantPointer& mp=bundle[unsigned(bi)];
    Integer idind=0;
    while ((idind<k)&&(! mp.equals(data.bundle[unsigned(idind)],eps_Real)))
      idind++;
    if (idind<k) {
      data.bundlecoeff(idind)+=coeff(bi);
    }
    else {
      if (k<kmax){
	data.bundle[unsigned(k)]=mp;
	data.bundlecoeff(k)=coeff(bi);
	k++;
      }
      else {
	Real d=coeff(bi);
	if (d>1e-10*data.get_function_factor()){
	  aggrcoeff.concat_below(d);
	  aggrbundle.push_back(mp);
	}
      }
    }
  }
  if (aggrcoeff.rowdim()>0){
    assert(Integer(k)==kmax);
    MinorantPointer& aggr=data.bundle[unsigned(kmax-1)];
    Real d=data.bundlecoeff(kmax-1)+sum(aggrcoeff);
    aggr.scale(1./d);
    aggr.aggregate(aggrbundle,aggrcoeff,1./d);
    data.bundlecoeff(kmax-1)=d;
    aggrbundle.clear();
  }
   
  //if not all have been treated and there is still room, sort by objective value
  if ((i<coeff.rowdim())&&(k<kmax)){
    sortindex(objval,sind,false);
    while ((i<coeff.rowdim())&&(k<kmax)){
      Integer bi=sind(i++);
      if (bi==bestind)
	continue;
      const MinorantPointer& mp=bundle[unsigned(bi)];
      Integer idind=0;
      while ((idind<k)&&(! mp.equals(data.bundle[unsigned(idind)],eps_Real)))
	idind++;
      if (idind<k) //identical minorant 
	continue;
      data.bundle[unsigned(k)]=mp;
      data.bundlecoeff(k)=coeff(bi);
      k++;
    }
  }

  if (k<kmax){
    data.bundlecoeff.reduce_length(k);
    data.bundle.resize(unsigned(k));
  }

  data.primal_tapia=(data.bundlecoeff>1e-8*data.get_function_factor());

  return 0;
}


// *****************************************************************************
//                            evaluate_oracle
// *****************************************************************************


  int NNCModel::evaluate_oracle(Integer y_id,
				const Matrix& y,
				Real nullstep_bound,
				Real relprec)
{
  assert(relprec>0.);

  CH_Tools::Microseconds start_eval=clock.time();

  data.cand_minorants.clear();
  data.cand_ub=nullstep_bound;
  data.cand_id=y_id;
  data.cand_ub_mid=-1;  //signals no valid computation yet

  //---- evaluate  	
  int err=0;
  std::vector<Minorant*> minorants;
  PrimalExtender* pep=0;			
  nr_eval++;
  preeval_time+=clock.time()-start_eval;  
  start_eval=clock.time();
  ret_code = oracle->evaluate(y,relprec,data.cand_ub,
			      minorants,pep);
  eval_time+=clock.time()-start_eval;
  start_eval=clock.time();

  if (ret_code){
    if (cb_out()) get_out()<<"**** WARNING: NNCModel::eval_function(): function returned code = "<<ret_code<<std::endl;
  }

  if (minorants.size()==0) {
    if (cb_out()) 
      get_out()<<"**** ERROR: NNCModel::eval_function(): function returned no subgradients"<<std::endl;
    posteval_time+=clock.time()-start_eval;
    return 1;
  }

  data.cand_ub_mid=data.get_modification_id();
 
  //--- first apply the modifications to the old data
  if (pep){
    if (call_primal_extender(*pep))
      if (cb_out()) 
	get_out()<<"**** WARNING: NNCModel::eval_function(): data.call_primal_extender failed"<<std::endl;
    delete pep;
  }

  //---- restructure (and check) the minorants
  bool check_data=oracle->check_correctness();
  Real cand_lb=CB_minus_infinity;
  Real max_lb=cand_lb;
  
  data.cand_minorants.resize(minorants.size());
  for(unsigned int i=0;i<minorants.size();i++){
 
    Real check_value=0.;
    bool check_value_init=false;

    if ((check_data)&&(! minorants[i]->offset_gives_value_at_origin())){
      check_value=minorants[i]->offset();
      check_value_init=true;
    } 

    data.cand_minorants[i].init(minorants[i],data.get_modification_id());
    if (!minorants[i]->offset_gives_value_at_origin()){
      minorants[i]->add_offset(-data.cand_minorants[i].ip(y));
      minorants[i]->offset_gives_value_at_origin()=true;
    }

    if (check_data){
      if (!check_value_init) {
	check_value=data.cand_minorants[i].evaluate(y_id,y);
      }
      if (check_value>data.cand_ub+0.1*relprec*(fabs(data.cand_ub)+1.)){
	if (cb_out()){ 
	  get_out().precision(12);
	  get_out()<<"**** ERROR: NNCModel::eval_function(): the oracle returned an upper bound="<<data.cand_ub<<" on the function value that is smaller than the lower bound="<<check_value<<" resulting from minorant "<<i<<std::endl;
	}
	err++;
      }
      max_lb=max(max_lb,check_value);
      if (cand_lb<check_value){
	if (i==0)
	  cand_lb=check_value;
	else {
	  if (cb_out()){ 
	    get_out().precision(12);
	    get_out()<<"**** ERROR: NNCModel::eval_function(): the value="<<check_value<<" of minorant "<<i<<" exceeds the value="<<cand_lb<<" of minorant 0 which should be largest"<<std::endl;
	  }
	  err++;
	}
      }	  
    } 
  }
  data.cand_minorant=data.cand_minorants[0];

  if (check_data){
    if ((max_lb<nullstep_bound)&&
	(data.cand_ub-max_lb>relprec*(fabs(data.cand_ub)+1.))){
      if (cb_out()) 
	get_out()<<"**** WARNING: NNCModel::eval_function(): insufficient precision in evaluation routine"<<std::endl;
      if (err==0) 
	err=-1;
    }
  }

  if (data.cand_ub>=CB_plus_infinity){
    if (cb_out()) 
      get_out()<<"**** WARNING: NNCModel::eval_function(): function value cand_ub>=CB_plus_infinity"<<std::endl;
  }    
  
  posteval_time+=clock.time()-start_eval;
  return (err>0)?err:ret_code;
}

// *****************************************************************************
//                              evalaluate_local_model
// *****************************************************************************

//evaluate the current cutting model in $y$ 
//if evaluated by an iterative method that provides upper and lower bounds
//it may stop when the lower bound (lb) is above the nullstep_bound
//evalutation may also stop, if (ub-lb)<relprec*(|ub|+1) is satisfied

int NNCModel::evaluate_local_model(Real& lb,
				   Integer y_id,
				   const Matrix& y,
				   Real /* relprec */)
{
  assert(data.get_local_model_active());
  assert(data.bundle.size()>0);
  assert(data.bundle[0].valid());

  lb=data.bundle[0].evaluate(y_id,y);
  for (unsigned i=1;i<data.bundle.size();i++){
    lb=max(lb,data.bundle[i].evaluate(y_id,y));
  }

  return 0;
}


// *****************************************************************************
//                                update_local_model
// *****************************************************************************


int NNCModel::update_local_model(ModelUpdate model_update,
				 Integer center_id,
				 const Matrix& center_y,
				 Integer y_id,
				 const Matrix& y,
				 Real model_maxviol,
				 BundleProxObject& H)
{
  assert(unsigned(data.bundlecoeff.dim())==data.bundle.size());

  if ((model_update!=new_subgradient)&&(data.cand_ub_mid!=data.get_modification_id())){   
    if (cb_out())
      get_out()<<"**** ERROR NNCModel::update_model: candidate required but its values are not up to date"<<std::endl;
    return 1;
  }

  if (model_update==descent_step){
    data.do_step(y_id);
  }

  //initialize the bundle if necessary
  if (data.bundle.size()==0){
    data.outside_bundle=false;
    if (data.cand_ub_mid!=data.get_modification_id()){
      //input_id must be center_id
      if ((y_id!=data.get_center_id())||(data.get_center_modification_id()!=data.get_modification_id())){
	if (cb_out())
	  get_out()<<"**** ERROR NNCModel::update_model: point_id should match recomputed center but it does not or center is not up to date"<<std::endl;
	return 1;	
      }
      data.bundle.push_back(data.center_minorant);
      data.bundlecoeff.init(1,1,data.get_function_factor());
      data.primal_tapia.init(1,1,1.);
    }
    else if ((model_selection->get_max_model_size()<2)||(unsigned(model_selection->get_max_model_size())>=data.cand_minorants.size())){
      data.bundle=data.cand_minorants;
      data.bundlecoeff.init(Integer(data.bundle.size()),1,0.);
      data.bundlecoeff(0)=data.get_function_factor();
      data.primal_tapia=(data.bundlecoeff>1e-8*data.get_function_factor());
    }
    else {
      Matrix tmpvec(Integer(data.cand_minorants.size()),1); chk_set_init(tmpvec,1);
      for (unsigned i=0;i<data.cand_minorants.size();i++){
	tmpvec(Integer(i))=data.cand_minorants[i].evaluate(y_id,y);
      }
      Indexmatrix tmpind;
      sortindex(tmpvec,tmpind,false);
      tmpind.reduce_length(model_selection->get_max_model_size());
      data.bundle.resize(unsigned(model_selection->get_max_model_size()));
      for (Integer i=0;i<model_selection->get_max_model_size();i++)
	data.bundle[unsigned(i)]=data.cand_minorants[unsigned(tmpind(i))];
      data.bundlecoeff.init(model_selection->get_max_model_size(),1,0.);
      data.bundlecoeff(0)=data.get_function_factor();
      data.primal_tapia=(data.bundlecoeff>1e-8*data.get_function_factor());
    }
    
    return 0;
  }


  if ((!data.get_local_model_active())||(model_update==new_subgradient)) {
    //it is a root or a child in sumbundle, all is taken care of by sumbundle
    //but if desired here the local bundle could be updated along, maybe with more
    
    return 0;
  }

  //-- bring the bundle into shape if it comes from somwhere else 
  if (data.outside_bundle)
    adapt_bundle_to_parameters(y_id,y);

  //-- prepare for calling the model selection routine
  MinorantBundle latest_minorants;
  const MinorantBundle* omp;
  if ((model_selection->get_max_bundle_size()<2)||
      (model_selection->get_max_bundle_size()>=int(data.get_old_minorants().size()))
      ){
    omp=&(data.get_old_minorants());
  }
  else {
    omp=&latest_minorants;
    data.get_latest_minorants(latest_minorants,model_selection->get_max_bundle_size());
  }

  //-- call the model selection routine
  int retval=model_selection->select_model(data.bundle,
					   data.bundlecoeff,
					   data.primal_tapia,
					   data.local_aggregate,
					   data.center_minorant,
					   data.cand_minorants,
					   *omp,
					   oracle,
					   data.get_function_task(),
					   data.get_function_factor(),
					   model_update,
					   center_id,center_y,
					   y_id,y,
					   model_maxviol,
					   H);
  assert(data.bundle.size()>0);
  assert(Integer(data.bundle.size())==data.bundlecoeff.rowdim());
  assert(sum(data.bundlecoeff)<=data.get_function_factor()*(1+1e-8));

  if ((retval)&&(cb_out())){
    get_out()<<"**** WARNING NNCModel::update_model: model_selection->select_model() returned "<<retval<<std::endl;
  }
    
  if (cb_out(2)){
    get_out()<<" NNCModel::update_model: bundle_size="<<data.bundle.size()<<std::endl;
  }
  
  //Test on the quality of the aggregate
  // if (cb_out(2)){
  //   Matrix tmpvec(y.rowdim(),1,0.);
  //   Real val=0.;
  //   data.aggregate.get_minorant(val,tmpvec,1,1.);
  //   for(Integer i=0;i<data.bundlecoeff.rowdim();i++)
  //     data.bundle[unsigned(i)].get_minorant(val,tmpvec,1,-data.bundlecoeff(i),true);
  //   get_out()<<" Fun update aggrdevnorm="<<norm2(tmpvec)<<" valdev="<<std::fabs(val)<<std::endl;
  // }
  
  
  return retval;

}
  

 
// *****************************************************************************
//                        provide_local_aggregate
// *****************************************************************************


int NNCModel::provide_local_aggregate(Integer /* y_id */ ,const Matrix& /* iny */)
{
  if (data.local_aggregate.valid()) {
    return 0;
  }
  
  //--- form an aggregate
  data.local_aggregate.clear();

  assert(data.get_local_model_active());
  assert(data.bundle.size()>0);
  assert((data.get_function_task()!=ObjectiveFunction)||(std::fabs(sum(data.bundlecoeff)-data.get_function_factor())<1e-6*data.get_function_factor()));

  if (data.local_aggregate.aggregate(data.bundle,data.bundlecoeff)){
    if (cb_out())
      get_out()<<"**** ERROR NNCModel::provide_model_aggregate: forming the aggregate failed"<<std::endl;
    return 1;
  }

  return 0;
}


// *****************************************************************************
//                         switch_to_local
// *****************************************************************************

int NNCModel::switch_to_local( const MinorantBundle& current_model,
			     const Matrix& current_coeff)
 {
   assert(Integer(current_model.size())==current_coeff.rowdim());
   data.bundle=current_model;
   data.bundlecoeff=current_coeff;
   data.outside_bundle=true;
   return 0;
 }

// *****************************************************************************
//                        get_aggregate_and_coefficient
// *****************************************************************************

int NNCModel::get_aggregate_and_coefficient(MinorantPointer& sumbundle_aggr,
					    Real& aggr_coeff)
 {
   //set new model if required
   assert(data.bundle.size()>0);
   aggr_coeff=sum(data.bundlecoeff);
   if (aggr_coeff>eps_Real){
     if (data.local_aggregate.valid()){
       sumbundle_aggr.init(data.local_aggregate,1./aggr_coeff);
     }
     else {
       sumbundle_aggr.clear();
       sumbundle_aggr.aggregate(data.bundle,data.bundlecoeff,1./aggr_coeff);
     }
   }
   else {
     sumbundle_aggr=data.bundle[0];
   }
   return 0;
 }
 
  

// *****************************************************************************
//                                start_local_augmodel
// *****************************************************************************

//return a pointer to the variables/constraints generating the cutting model
//returns 0 on success, 1 on failure

  int NNCModel::start_local_augmodel(QPModelDataPointer& qpblockp,
				     Integer cand_id,
				     const Matrix& cand_y,
				     const Indexmatrix* /* needed_indexsubset */)
{
  if (qpblockp.set_model_data(block)){
    delete block;
    block=qpblockp.generate_conemodel_data();
    block->set_cbout(this);
  }
  assert(block);
  
  if (data.outside_bundle){
    adapt_bundle_to_parameters(cand_id,cand_y);
  }

  Integer xdim=Integer(data.bundle.size());
  Indexmatrix tmpind(0,1,Integer(0));
  Matrix dummymat(0,1,0.);
  int err=0;
  if (block->init(MinorantPointer(),data.bundle,xdim,tmpind,tmpind,dummymat,dummymat,data.get_function_factor(),data.get_function_task())){
    if (cb_out()){
      get_out()<<"*** ERROR: NNCModel::start_local_augmodel(): initializing the QP data with block->init() failed"<<std::endl;
    }
    err++;
  }
    
  qpblockp.set_model_data(block);
  assert(qpblockp.get_model_data_ptr()==block);

  if (cb_out(2)){
    get_out()<<" NNCxdim="<<xdim;
  }

  return 0;
}
 
// *****************************************************************************
//                        make_local_aggregate
// *****************************************************************************

  int NNCModel::make_local_aggregate(Real& new_function_factor,
				     bool keep_penalty_fixed)  
{
  assert(new_function_factor==data.get_function_factor());
  
  bool penalty_changed=false;
  data.local_aggregate.clear();
  assert(block);
  
  if ((!keep_penalty_fixed) && (data.get_function_task()==AdaptivePenaltyFunction)){
    // increasing the penalty only makes sense if the constraint is active
    // and the model does not promise a strong reduction anyways.
    Real prec=0;
    Real ubmodelval=block->tracedual(&prec);
    Real traceval=block->evaluate_trace();
    if ((ubmodelval>100.*prec)&&
	(data.get_center_modification_id()==data.get_modification_id())&&
	(ubmodelval>.9*data.get_center_ub())
	){
      new_function_factor=max(2*data.get_function_factor(),data.get_function_factor()+ubmodelval);
      penalty_changed=true;
    }
    else if ((ubmodelval<0.01*prec)&&
	     (data.get_center_modification_id()==data.get_modification_id())&&
	     (data.get_function_factor()>1.)&&
	     (10.*traceval<data.get_function_factor())
	     ){
      new_function_factor=max(1.,max(0.01*data.get_function_factor(),5.*traceval));
      penalty_changed=true;
    }

    if (penalty_changed) {
      block->adjust_trace(new_function_factor);
    }
  }

  
  if (!penalty_changed){

    //get coefficients and activity
    if (block->get_nncx(data.bundlecoeff,&data.primal_tapia)) {
      if (cb_out()){
	get_out()<<"**** ERROR NNCModel::make_model_aggregate(): retrieving the QP solution failed"<<std::endl;
    }
      return 1;
    }
    
    if (data.bundlecoeff.dim()!=Integer(data.bundle.size())){
      if (cb_out()){
	get_out()<<"**** ERROR NNCModel::make_model_aggregate(): the retrieved QP solution does not match the bundle size"<<std::endl;
      }
      return 1;
    }
    
    //form aggregate
    if (data.local_aggregate.aggregate(data.bundle,data.bundlecoeff)){
      if (cb_out()){
	get_out()<<"**** ERROR NNCModel::make_model_aggregate(): forming the aggregate failed"<<std::endl;
      }
      return 1;
    }
  }
  
  return 0;
} 

// *****************************************************************************
//                                 lb_function
// *****************************************************************************

//returns a *quick* lower bound for the function value at y
//(eg by a previous subgradient)

  Real NNCModel::lb_function(Integer y_id,const Matrix& y)
{

  Real lb=CB_minus_infinity;

  if (data.get_aggregate().valid()){
    lb=data.get_aggregate().evaluate(y_id,y);
  }
  else if (data.bundle.size()>0){
    lb=data.bundle[0].evaluate(y_id,y);
    for (unsigned i=1;i<data.bundle.size();i++){
      lb=max(lb,data.bundle[i].evaluate(y_id,y));
    }
    lb=data.get_function_factor()*lb;
  }
  else if ((data.cand_ub_mid==data.get_modification_id())&&
	   (data.cand_minorants.size()>0)){
    lb=data.cand_minorants[0].evaluate(y_id,y);
    for (unsigned i=1;i<data.cand_minorants.size();i++){
      lb=max(lb,data.cand_minorants[i].evaluate(y_id,y));
    }
    lb=data.get_function_factor()*lb;
  }

  if (data.get_function_task()!=ObjectiveFunction)
    lb=max(0.,lb);

  return lb;
}

// *****************************************************************************
//                           provide_cand__minorant
// *****************************************************************************


MinorantPointer& NNCModel::provide_cand_minorant()
{
  if (data.cand_ub_mid!=data.get_modification_id())
    data.cand_minorant.clear();
  return data.cand_minorant;
}

// *****************************************************************************
//                            provide_center_minorant
// *****************************************************************************


MinorantPointer& NNCModel::provide_center_minorant()
{
  if (data.get_center_modification_id()!=data.get_modification_id())
    data.center_minorant.clear();

  return data.center_minorant;
}

// *****************************************************************************
//                           local_adjust_multiplier
// *****************************************************************************


int NNCModel::local_adjust_multiplier(Real& new_function_factor)
{
  assert(new_function_factor==data.get_function_factor());

  if (data.get_function_task()!=AdaptivePenaltyFunction){
    return 0;
  }

  assert(block);
  Real oldval=data.get_function_factor();
  Real traceval=block->evaluate_trace();
  new_function_factor=min(data.get_function_factor(),max(1.,1.5*traceval));
  if (new_function_factor<oldval){
    block->adjust_trace(new_function_factor);
  }

  return 0;
}


// *****************************************************************************
//                         set_bundle_parameters      
// *****************************************************************************

int NNCModel::set_bundle_parameters(const BundleParameters& bp)
{
  const NNCModelParametersObject* mpo=dynamic_cast<const NNCModelParametersObject*>(&bp);
  int retval=0;
  if (mpo){
    delete model_selection;
    model_selection=dynamic_cast<NNCModelParametersObject*>(mpo->clone_BundleParameters());
    assert(model_selection);
  }
  else {
    assert(model_selection);
    retval=model_selection->init(bp);
  }
  return retval;
}


} 
