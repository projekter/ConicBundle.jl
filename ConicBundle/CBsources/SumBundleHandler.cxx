/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/SumBundleHandler.cxx
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



#include "SumBundleHandler.hxx"
#include "SumBundleParameters.hxx"
#include "BundleIdProx.hxx"

using namespace CH_Matrix_Classes;

namespace ConicBundle {

    ///
    SumBundleHandler::BundleInformation::BundleInformation():block(0)
    {sbp=new SumBundleParameters(2,2,0);bundle_size=aggr_index=new_index=0;
      increase_factor=1.;}
  
    ///
    SumBundleHandler::BundleInformation::~BundleInformation()
    {delete sbp;delete block;}

    ///only copies the information regarding the bundle size and initializes some members
    SumBundleHandler::BundleInformation&
    SumBundleHandler::BundleInformation::operator=(const SumBundleHandler::BundleInformation& bi)
    {
      bundle_size=bi.bundle_size;
      delete sbp;
      sbp=dynamic_cast<SumBundleParametersObject*>(bi.sbp->clone_BundleParameters());
      delete block;
      block=0; //simply rebuild it next time
      assert(sbp);
      aggr_index=bi.aggr_index;
      new_index=bi.new_index;
      map_to_old=bi.map_to_old;
      increase_factor=1.;
      return *this;
    }


// *****************************************************************************
//                              
// *****************************************************************************

  int SumBundleHandler::init(FunctionTask ft,Real funfactor,const BundleParameters* bp)
{
  if (!sumbundle->has_bundle_for(ft)){
    sumbundle->init(ft,funfactor);
    if (sumbundle->get_mode(ft)==SumBundle::unavailable)
      sumbundle->set_mode(ft)=SumBundle::inactive;
  }
  BundleInformation* bi;
  if (bundleinfo[ft]==0){
    bi=new BundleInformation;
    bundleinfo[ft]=bi;
    bi->bundle_size=sumbundle->bundle_size(ft);
  }
  else {
    bi=bundleinfo[ft];
    assert((sumbundle->get_n_contributors(ft)==0)||(bi->bundle_size==sumbundle->bundle_size(ft)));
  }

  sumbundle->set_function_factor(ft)=funfactor;
  store_aggregate(ft);

  Integer max_model_size=bp? bp->get_max_model_size():0;
  Integer max_bundle_size=bp? bp->get_max_bundle_size():0;

  if (max_model_size<2){
    max_model_size=max(bi->sbp->get_max_model_size(),2);
  }

  if (max_bundle_size<max_model_size){
    max_bundle_size=max(bi->bundle_size,max(50,2*max_model_size));
  }

  const SumBundleParametersObject* insbp=dynamic_cast<const SumBundleParametersObject*>(bp);
  if (insbp){
    delete bi->sbp;
    bi->sbp=dynamic_cast<SumBundleParametersObject*>(insbp->clone_BundleParameters());
  }
    
  bi->sbp->set_max_model_size(max_model_size);  
  bi->sbp->set_max_bundle_size(max_bundle_size);
  bi->sbp->set_cbout(this,0);

  if (bi->bundle_size>max_bundle_size){
    assert((bi->aggr_index>=0)&&(bi->aggr_index<bi->bundle_size));
    assert((bi->new_index>=0)&&(bi->new_index<bi->bundle_size));
    bi->map_to_old.init(max_bundle_size,1,Integer(0));
    //usually bi->aggr_index==0; if not, set this initially
    Integer aggr_index=0;
    Integer new_index=min(bi->new_index,max_bundle_size-1);
    if (bi->aggr_index!=aggr_index){
      bi->map_to_old(aggr_index)=bi->aggr_index;
    }
    //assuming cyclic increase of new_index keep the most recent vectors
    Integer cnt=max_bundle_size-1;
    Integer pos=new_index;
    Integer from=bi->new_index;
    while (cnt>0){
      if (from==bi->aggr_index){
	--from;
      }
      if (from<0){
	from=bi->bundle_size-1;
	continue;
      }
      if (pos==0){
	pos=bi->sbp->get_max_bundle_size()-1;
      }
      bi->map_to_old(pos)=from;
      --cnt;
      --from;
      --pos;
    }
    reassign_minorants(bi->map_to_old,ft);
    bi->aggr_index=aggr_index;
    bi->new_index=new_index;
    assert(bi->bundle_size==max_bundle_size);
  } 
  else {
    bi->map_to_old.init(0,0,Integer(0));
  }
  return 0;
}

// *****************************************************************************
//                              
// *****************************************************************************

  int SumBundleHandler::set_bundle_parameters(FunctionTask ft,const BundleParameters& bp)
{
  if(!sumbundle->has_bundle_for(ft))
    return 0;

  Real f=sumbundle->get_function_factor(ft);

  BundleParameters *p=bp.clone_BundleParameters();
  
  if (p->get_max_model_size()<2){
    p->set_max_model_size(2);
  }

  if (p->get_max_bundle_size()<p->get_max_model_size()){
    p->set_max_bundle_size(max(max(50,2*p->get_max_model_size()),sumbundle->bundle_size(ft)));
  }
  if ((parent_handler!=0)&&(parent_handler->bundleinfo[ft]!=0)){
    BundleInformation* bhi=parent_handler->bundleinfo[ft];
    p->set_max_model_size(max(p->get_max_model_size(),bhi->sbp->get_max_model_size()));
    p->set_max_bundle_size(max(p->get_max_bundle_size(),bhi->sbp->get_max_bundle_size()));
  }
  int retval=init(ft,f,p);
  delete p;
  return retval;
}

// *****************************************************************************
//                              
// *****************************************************************************

int SumBundleHandler::set_bundle_parameters(const BundleParameters& bp)
{
  int err=0;
  err += set_bundle_parameters(ObjectiveFunction,bp);
  err += set_bundle_parameters(ConstantPenaltyFunction,bp);
  err += set_bundle_parameters(AdaptivePenaltyFunction,bp);
  return err;
}

// *****************************************************************************
//                              
// *****************************************************************************

  SumBundleHandler::SumBundleHandler(SumBundle& sb,const std::map<FunctionTask,Real>& factor_map,const BundleParameters* bp)
{
  sumbundle=&sb;
  parent_handler=0;
  aft=0;
  bundleinfo.resize(3,0);
  
  if(factor_map.find(ObjectiveFunction)==factor_map.end()){
    sumbundle->clear(ObjectiveFunction);
    sumbundle->set_mode(ObjectiveFunction)=SumBundle::unavailable;
  }
  if(factor_map.find(ConstantPenaltyFunction)==factor_map.end()){
    sumbundle->clear(ConstantPenaltyFunction);
    sumbundle->set_mode(ConstantPenaltyFunction)=SumBundle::unavailable;
  }
  if(factor_map.find(AdaptivePenaltyFunction)==factor_map.end()){
    sumbundle->clear(AdaptivePenaltyFunction);
    sumbundle->set_mode(AdaptivePenaltyFunction)=SumBundle::unavailable;
  }

  for (std::map<FunctionTask,Real>::const_iterator it=factor_map.begin();
       it!=factor_map.end();it++){
    assert(it->second>0.);
    init(it->first,it->second,bp);
  }

  
}

// *****************************************************************************
//                              
// *****************************************************************************

SumBundleHandler::~SumBundleHandler()
{
  for (unsigned long i=0;i<bundleinfo.size();i++)
    delete bundleinfo[i];
  bundleinfo.clear();
}

// *****************************************************************************
//                              
// *****************************************************************************

int SumBundleHandler::store_aggregate(FunctionTask ft)
{
  if ((!sumbundle->has_bundle_for(ft))||(sumbundle->get_n_contributors(ft)==0)||(sumbundle->bundle_size(ft)==0))
    return 0;

  BundleInformation* bi=bundleinfo[ft];

  assert(bi->bundle_size==sumbundle->bundle_size(ft));
  assert(bi->bundle_size==sumbundle->get_coeff(ft).dim());
  assert(bi->sbp->get_max_bundle_size()>=2);
  assert((bi->aggr_index>=0)&&(bi->aggr_index<bi->bundle_size));

  //if the aggregate is already in its place and at its bound, do nothing
  if (CH_Matrix_Classes::abs(sumbundle->get_coeff(ft)(bi->aggr_index)-sumbundle->get_function_factor(ft))<1e-10*sumbundle->get_function_factor(ft)){
    return 0;
  }

  //if the aggregate is not yet available, form it
  if (!sumbundle->get_aggregate(ft).valid()){
    sumbundle->set_aggregate(ft).aggregate(sumbundle->get_bundle(ft),sumbundle->get_coeff(ft));
  }

  

  //compute the coefficients for "normalizing" the aggregate
  Real factor=sum(sumbundle->get_coeff(ft))/sumbundle->get_function_factor(ft);
  sumbundle->set_coeff(ft)*=factor;

  //copy the aggregate into its column
  sumbundle->set_bundle(ft)[unsigned(bi->aggr_index)].init(sumbundle->get_aggregate(ft),1./sumbundle->get_function_factor(ft));

  //set the new coefficient vector
  sumbundle->set_coeff(ft).init(bi->bundle_size,1,0.);
  sumbundle->set_coeff(ft)(bi->aggr_index)=sumbundle->get_function_factor(ft);

  return 0;
}


// *****************************************************************************
//                              
// *****************************************************************************

int SumBundleHandler::reassign_minorants(const Indexmatrix& map_to_old,FunctionTask ft)
{
  if (!sumbundle->has_bundle_for(ft))
    return 0;

  if ((sumbundle->get_n_contributors(ft)==0)||(sumbundle->bundle_size(ft)==0)){
    //no contributors, just organizational for potential children
    BundleInformation* bi=bundleinfo[ft];
    bi->bundle_size=map_to_old.dim();
    return 0;
  }

  BundleInformation* bi=bundleinfo[ft];

  assert(bi->bundle_size==sumbundle->bundle_size(ft));
  assert(bi->bundle_size==sumbundle->get_coeff(ft).dim());
  assert((0<=min(map_to_old))&&(max(map_to_old)<bi->bundle_size));
  assert(bi->sbp->get_max_bundle_size()>=2);
  assert((bi->aggr_index>=0)&&(bi->aggr_index<bi->bundle_size));
  assert(! sumbundle->get_bundle(ft)[unsigned(bi->new_index)].empty());

  //check correctness and whether a copy is needed to do the job
  Indexmatrix map_to(bi->bundle_size,1,Integer(-1));
  bool need_copy=false;
  for (Integer i=0;i<map_to_old.dim();i++){
    Integer ind=map_to_old(i);
    if (map_to(ind)<0){
      //this is the first occurence of ind and i will get the value of coeff(ind)
      map_to(ind)=i;
    }
    if ((ind<i)&&(map_to(ind)!=ind)){
      need_copy=true;
    }
  }


  int err=0;
  for (Integer i=0;i<bi->bundle_size;i++){
    if (map_to(i)<0){
      if ((map_to_old.dim()>0)&&(sumbundle->get_coeff(ft)(i)!=0.)){
	if (cb_out()){
	  get_out()<<"**** ERROR SumBundleHandler::reassign_minorants(..): minorant "<<i<<" is deleted but has coeff("<<i<<")="<<sumbundle->get_coeff(ft)(i)<<std::endl;
	}
	err++;
      }
    }
  }

  if (!need_copy){

    if (sumbundle->bundle_size(ft)<map_to_old.dim()){
      Integer nnew=map_to_old.dim()-sumbundle->bundle_size(ft);
      sumbundle->set_bundle(ft).resize((unsigned long)(map_to_old.dim()));
      sumbundle->set_coeff(ft).concat_right(Matrix(nnew,1,0.));
    }

    for(Integer i=0;i<map_to_old.dim();i++){
      Integer ind=map_to_old(i);
      if (i==ind)
	continue;
      sumbundle->set_bundle(ft)[unsigned(i)]=sumbundle->get_bundle(ft)[unsigned(ind)];
      if(i==map_to(ind)){
	sumbundle->set_coeff(ft)(i)=sumbundle->get_coeff(ft)(ind);
      }
      else {
	sumbundle->set_coeff(ft)(i)=0.;
      }
    }

    if (bi->bundle_size>map_to_old.dim()){
      sumbundle->set_bundle(ft).resize((unsigned long)(map_to_old.dim()),0);
      sumbundle->set_coeff(ft).reduce_length(map_to_old.dim());
    }
  }
  else {
    MinorantBundle tmpbundle((unsigned long)(map_to_old.dim()));
    Matrix coeff(sumbundle->get_coeff(ft)(map_to_old));
    for (Integer i=0;i<map_to_old.dim();i++){
      tmpbundle[unsigned(i)]=sumbundle->get_bundle(ft)[unsigned(i)];
      if (i!=map_to(map_to_old(i))){
	coeff(i)=0.;
      }
    }
    assert(sum(coeff)-sum(sumbundle->get_coeff(ft))<1e-6*sum(coeff));
    swap(sumbundle->set_coeff(ft),coeff);
    swap(sumbundle->set_bundle(ft),tmpbundle);
  }
  
  assert(sumbundle->get_coeff(ft).dim()==map_to_old.dim());
  assert(sumbundle->get_bundle(ft).size()==unsigned(map_to_old.dim()));
  bi->bundle_size=map_to_old.dim();

  return err;
}

// *****************************************************************************
//                              
// *****************************************************************************

bool SumBundleHandler::initialization_needed(FunctionTask ft) const
{
  if ((!sumbundle->has_bundle_for(ft))||
      (sumbundle->get_mode(ft)!=SumBundle::root)||
      (sumbundle->get_n_contributors(ft)==0)||
      (sumbundle->bundle_size(ft)>0))
    return false;
  return true;
}

// *****************************************************************************
//                              
// *****************************************************************************

bool SumBundleHandler::initialization_needed() const
{
  bool needed=false;
  needed|=initialization_needed(ObjectiveFunction);
  needed|=initialization_needed(ConstantPenaltyFunction);
  needed|=initialization_needed(AdaptivePenaltyFunction);
  return needed;
}

// *****************************************************************************
//                              
// *****************************************************************************

  int SumBundleHandler::set_parent_information(SumBundleHandler *parent_sbh,const AffineFunctionTransformation* in_aft,SumBundle::Mode mode)
{
  bool use_parent_update=false;
  if ((parent_handler)&&(parent_sbh)){
    //in this case the next parent is assumed to be compatible with the old one
    use_parent_update=true;
  }
    
  parent_handler=parent_sbh;
  aft=in_aft;

  if (parent_handler==0){
    sumbundle->take_control(mode);
  }
  
  if (align_bundle(use_parent_update)){
    if (cb_out()){
      get_out()<<"*** ERROR: SumBundleHandler::set_parent_information(): align_bundle() failed"<<std::endl;
    }
    return 1;
  }

  if (mode!=SumBundle::child){
    if (remove_contributions()){
      if (cb_out()){
	get_out()<<"*** ERROR: SumBundleHandler::set_parent_information(): remove_contributions() failed"<<std::endl;
      }
      return 1;
    }
    sumbundle->take_control(mode);
  }
  else if (!sumbundle->active()) {
    sumbundle->take_control(mode);
  }    

  return 0;
}

// *****************************************************************************
//                              
// *****************************************************************************
  
  int SumBundleHandler::reset_function_factor(FunctionTask ft,Real factor)
{
  if (! sumbundle->has_bundle_for(ft)){
    return 0;
  }

  if (remove_contributions(ft)){
    if (cb_out()){
      get_out()<<"**** ERROR SumBundleHandler::reset_function_factor(..): remove_contribution failed for FunctionTask"<<ft<<std::endl;
    }
    return 1;
  }

  assert(factor>0.);
  sumbundle->set_function_factor(ft)=factor;
  sumbundle->set_aggregate(ft).clear();

  return 0;
}

// *****************************************************************************
//                              
// *****************************************************************************

Real SumBundleHandler::get_increase_factor() const
{ 
  if (sumbundle->has_bundle_for(AdaptivePenaltyFunction)){
    assert(bundleinfo[AdaptivePenaltyFunction]);
    return bundleinfo[AdaptivePenaltyFunction]->increase_factor; 
  }
  return 1.;
}


// *****************************************************************************
//                              
// *****************************************************************************

int SumBundleHandler::update_model(BundleModel::ModelUpdate model_update,
				   Integer /* center_id */,
				   const Matrix& /* center_y */,
				   Integer y_id,
				   const Matrix& y,
				   Real model_maxviol,
				   BundleProxObject& H,
				   FunctionTask ft)
{
  if (!sumbundle->has_bundle_for(ft))
    return 0;

  assert(bundleinfo[ft]);
  BundleInformation* bi=bundleinfo[ft];

  bi->model_indices.init(0,0,Integer(0));

  //Integer n=y.dim();

  if (sumbundle->get_n_contributors(ft)==0){
    sumbundle->init(ft,sumbundle->get_function_factor(ft));
    //here is no bundle to update but the rules have to be communicated
    //if there is a parent copy this only concerns the size informations 
    if ((parent_handler!=0)&&(parent_handler->bundleinfo[ft])){
      BundleInformation* bhi=parent_handler->bundleinfo[ft];
      assert(bhi->bundle_size<=bhi->sbp->get_max_bundle_size());
      assert(bhi->new_index<bhi->bundle_size);
      assert(bhi->aggr_index<bhi->bundle_size);
      bi->bundle_size=bhi->bundle_size;
      bi->aggr_index=bhi->aggr_index;
      bi->new_index=bhi->new_index;
      bi->map_to_old=bhi->map_to_old;
      return 0;
    }
    //without parent, use default values
    bi->aggr_index=0;
    bi->map_to_old.init(0,0,Integer(0));
    if (bi->bundle_size<bi->sbp->get_max_bundle_size()){
      bi->new_index=bi->bundle_size;
      bi->bundle_size++;
    }
    else {
      bi->new_index++;
      if (bi->new_index>=bi->sbp->get_max_bundle_size())
	bi->new_index=1;
    }
    return 0;
  }

  assert(bi->bundle_size==sumbundle->bundle_size(ft));
  assert(bi->bundle_size==sumbundle->get_coeff(ft).dim());
  assert(bi->sbp->get_max_bundle_size()>=2);

  //if not yet the case, store the aggregate in bi->aggr_index
  if (store_aggregate(ft)){
    if (cb_out()){
      get_out()<<"**** ERROR SumBundleHandler::update_model(..):  store_aggregate("<<ft<<") failed"<<std::endl;
    }
    return 1;
  }
 
  //the parent handler's choice has priority; if it exists, apply it 
  if ((parent_handler!=0)&&(parent_handler->bundleinfo[ft])){
    BundleInformation* bhi=parent_handler->bundleinfo[ft];
    assert(bhi->bundle_size<=bhi->sbp->get_max_bundle_size());
    assert(bhi->new_index<bhi->bundle_size);
    assert(bhi->aggr_index<bhi->bundle_size);

    bi->sbp->set_max_bundle_size(max(bi->sbp->get_max_bundle_size(),bhi->sbp->get_max_bundle_size()));
    bi->aggr_index=bhi->aggr_index;
    bi->new_index=bhi->new_index;

    //no excess bundle currently, so simply take over the rules of the parent
    bi->map_to_old.init(bhi->map_to_old);
    if (bi->map_to_old.dim()>0){
      if (reassign_minorants(bi->map_to_old,ft)){
	if (cb_out()){
	  get_out()<<"**** ERROR SumBundleHandler::update_model(..):  reassign_minorants failed for parents order"<<std::endl;
	}
	return 1;
      }
    }
  }
  else { //not parents_choice, use cyclic order 
    Integer aggr_index=0;
    //reorder and shrink if required assuming cyclic order for age of minorants
    if ((aggr_index!=bi->aggr_index)||(Integer(sumbundle->get_bundle(ft).size())>bi->sbp->get_max_bundle_size())){
      Integer ind=bi->new_index; //newest entry in the bundle
      Integer mapind=min(Integer(sumbundle->get_bundle(ft).size()),bi->sbp->get_max_bundle_size());
      bi->map_to_old.init(mapind,1,Integer(0));
      bi->new_index=--mapind;
      for(Integer i=mapind;i>=0;i--,ind--){
	if (ind<0)
	  ind=Integer(sumbundle->get_bundle(ft).size()-1);
	if (ind==bi->aggr_index){
	  bi->map_to_old(aggr_index)=ind;
	  continue;
	}
	bi->map_to_old(mapind)=ind;
	mapind--;
      }
      if (reassign_minorants(bi->map_to_old,ft)){
	if (cb_out()){
	  get_out()<<"**** ERROR SumBundleHandler::update_model(..):  reassign_minorants failed"<<std::endl;
	}
	return 1;
      }
      bi->aggr_index=aggr_index;
    }
    bi->new_index=min(bi->new_index+1,Integer(sumbundle->get_bundle(ft).size()));
    if (bi->new_index==bi->sbp->get_max_bundle_size())
      bi->new_index=1;
  }

  //place the new candidate at bi->new_index
  if (bi->new_index==bi->bundle_size){
    assert(bi->bundle_size<bi->sbp->get_max_bundle_size());
    bi->bundle_size++;
    sumbundle->set_bundle(ft).push_back(sumbundle->get_cand_minorant(ft));
    if ((bi->bundle_size==1)&&(ft==ObjectiveFunction))
      sumbundle->set_coeff(ft).concat_below(sumbundle->get_function_factor(ft));
    else
      sumbundle->set_coeff(ft).concat_below(0.);
  }
  else {
    assert(bi->new_index<bi->bundle_size);
    sumbundle->set_bundle(ft)[unsigned(bi->new_index)]=sumbundle->get_cand_minorant(ft);
    sumbundle->set_coeff(ft)(bi->new_index)=0.;
  }


  //if *this is a root, decide on which bundle minorants should be part of the model
  if (sumbundle->get_mode(ft)==SumBundle::root){
    bi->tmodel_maxviol=model_maxviol;
    bi->tweightu=H.get_weightu()*H.get_factor();
    bi->tmodel_update=model_update;
    bi->sbp->set_cbout(this,0);
    H.apply_factor(1./sumbundle->get_function_factor(ft));
    bi->sbp->select_model(bi->model_indices,y_id,y,
			  sumbundle->get_bundle(ft),0,
			  model_maxviol/sumbundle->get_function_factor(ft),
			  H,
			  model_update);
    H.apply_factor(sumbundle->get_function_factor(ft));
  }

  return 0;
}

// *****************************************************************************
//                              
// *****************************************************************************

int SumBundleHandler::update_model(BundleModel::ModelUpdate model_update,
				   Integer center_id,
				   const Matrix& center_y,
				   Integer y_id,
				   const Matrix& y,
				   Real model_maxviol,
				   BundleProxObject& H)
{
  int err=0;
  err+=update_model(model_update,center_id,center_y,y_id,y,model_maxviol,H,ObjectiveFunction);
  err+=update_model(model_update,center_id,center_y,y_id,y,model_maxviol,H,ConstantPenaltyFunction);
  err+=update_model(model_update,center_id,center_y,y_id,y,model_maxviol,H,AdaptivePenaltyFunction);
  return err;  
}

// *****************************************************************************
//                              
// *****************************************************************************

  int SumBundleHandler::eval_model(Real& lb,Integer y_id,const Matrix& y,FunctionTask ft) const
{
  if ((!sumbundle->has_bundle_for(ft))||
      (sumbundle->get_mode(ft)!=SumBundle::root)||
      (sumbundle->get_n_contributors(ft)==0)||
      (sumbundle->bundle_size(ft)==0))
    return 1;
  
  if (bundleinfo[ft]->model_indices.dim()>0){
    const Indexmatrix& model_indices=bundleinfo[ft]->model_indices;
    const MinorantBundle& bundle=sumbundle->get_bundle(ft);
    lb=bundle[unsigned(model_indices(0))].evaluate(y_id,y);
    for (Integer i=1;i<model_indices.dim();i++){
      Real v=bundle[unsigned(model_indices(i))].evaluate(y_id,y);
      if (v>lb)
	lb=v;
    }
  }
  else {
    const MinorantBundle& bundle=sumbundle->get_bundle(ft);
    lb=bundle[0].evaluate(y_id,y);
    for (unsigned i=1;i<bundle.size();i++){
      Real v=bundle[i].evaluate(y_id,y);
      if (v>lb)
	lb=v;
    }
  }
    
  if ((ft!=ObjectiveFunction)&&(lb<0.))
    lb=0.;
  else 
    lb*=sumbundle->get_function_factor(ft);

  return 0;
}

// *****************************************************************************
//                              
// *****************************************************************************

int SumBundleHandler::eval_model(Real& lb,Integer y_id,const Matrix& y) const
{
  bool initialized=false;
  int err=0;
  lb=0.;

  FunctionTask ft=ObjectiveFunction;
  if ((sumbundle->get_mode(ft)==SumBundle::root)&&
      (sumbundle->get_n_contributors(ft)>0)){
    err+=eval_model(lb,y_id,y,ft);
    initialized=true;
  }

  ft=ConstantPenaltyFunction;
  if ((sumbundle->get_mode(ft)==SumBundle::root)&&
      (sumbundle->get_n_contributors(ft)>0)){
    Real plb;
    err+=eval_model(plb,y_id,y,ft);
    lb+=plb;
    initialized=true;
  }

  ft=AdaptivePenaltyFunction;
  if ((sumbundle->get_mode(ft)==SumBundle::root)&&
      (sumbundle->get_n_contributors(ft)>0)){
    Real plb;
    err+=eval_model(plb,y_id,y,ft);
    lb+=plb;
    initialized=true;
  }

  if (!initialized)
    err++;

  return err;
}



// *****************************************************************************
//                              
// *****************************************************************************

Real SumBundleHandler::lb_model(Integer y_id,const Matrix& y,FunctionTask ft) const
{
  assert(sumbundle->has_bundle_for(ft));
  assert(sumbundle->get_mode(ft)==SumBundle::root);
  assert(sumbundle->get_n_contributors(ft)>0);

  Real lb;
  if (sumbundle->get_aggregate(ft).valid())
    lb=sumbundle->get_aggregate(ft).evaluate(y_id,y);
  else {
    if (bundleinfo[ft]->model_indices.dim()>0){
      const Indexmatrix& model_indices=bundleinfo[ft]->model_indices;
      const MinorantBundle& bundle=sumbundle->get_bundle(ft);
      lb=bundle[unsigned(model_indices(0))].evaluate(y_id,y);
      for (Integer i=1;i<model_indices.dim();i++){
	Real v=bundle[unsigned(model_indices(i))].evaluate(y_id,y);
	if (v>lb)
	  lb=v;
      }
    }
    else {
      const MinorantBundle& bundle=sumbundle->get_bundle(ft);
      lb=bundle[0].evaluate(y_id,y);
      for (unsigned i=1;i<bundle.size();i++){
	Real v=bundle[i].evaluate(y_id,y);
	if (v>lb)
	  lb=v;
      }
    }
    lb*=sumbundle->get_function_factor(ft);
  }

  if ((ft!=ObjectiveFunction)&&(lb<0.))
    lb=0.;

  return lb;
}

// *****************************************************************************
//                              
// *****************************************************************************

Real SumBundleHandler::lb_model(Integer y_id,const Matrix& y) const
{
  assert(sumbundle->has_working_roots());

  Real lb=0.;

  FunctionTask ft=ObjectiveFunction;
  if ((sumbundle->get_mode(ft)==SumBundle::root)&&
      (sumbundle->get_n_contributors(ft)>0)){
    lb=lb_model(y_id,y,ft);
  }

  ft=ConstantPenaltyFunction;
  if ((sumbundle->get_mode(ft)==SumBundle::root)&&
      (sumbundle->get_n_contributors(ft)>0)){
    lb+=lb_model(y_id,y,ft);
  }

  ft=AdaptivePenaltyFunction;
  if ((sumbundle->get_mode(ft)==SumBundle::root)&&
      (sumbundle->get_n_contributors(ft)>0)){
    lb+=lb_model(y_id,y,ft);
  }

  return lb;
}





// *****************************************************************************
//                              
// *****************************************************************************

int SumBundleHandler::contribute_new_minorants(FunctionTask ft)
{
  if ((!sumbundle->has_bundle_for(ft))
      ||(sumbundle->get_mode(ft)!=SumBundle::child)
      ||(sumbundle->get_n_contributors(ft)==0))
    return 0;

  assert(parent_handler);

  Real factor=sumbundle->get_function_factor(ft)/parent_handler->sumbundle->get_function_factor(ft);

  if (aft) {
    int retval=aft->transform_minorant(parent_handler->sumbundle->set_cand_minorant(ft),sumbundle->get_cand_minorant(ft),factor);
    if (retval){
      if (cb_out()){
	get_out()<<"**** ERROR SumBundleHandler::contribute_new_minorants(..):  aft->transform_minorant failed with return value "<<retval<<std::endl;
      }
      return 1;
    }
  }
  else {
    int retval=sumbundle->get_cand_minorant(ft).get_minorant(parent_handler->sumbundle->set_cand_minorant(ft),factor);
    if (retval){
      if (cb_out()){
	get_out()<<"**** ERROR SumBundleHandler::contribute_new_minorants(..):  get_minorant failed with return value "<<retval<<std::endl;
      }
      return 1;
    }
  }

  return 0;
}


// *****************************************************************************
//                              
// *****************************************************************************

int SumBundleHandler::contribute_new_minorants()
{
  int err=0;
  err+=contribute_new_minorants(ObjectiveFunction);
  err+=contribute_new_minorants(ConstantPenaltyFunction);
  err+=contribute_new_minorants(AdaptivePenaltyFunction);
  return err;  
}


// *****************************************************************************
//                              
// *****************************************************************************

  ///remove valid contributions and set the states correspondingly
int SumBundleHandler::normalize_sumbundle()
{
  int err=0;
  err+=store_aggregate(ObjectiveFunction);
  err+=store_aggregate(ConstantPenaltyFunction);
  err+=store_aggregate(AdaptivePenaltyFunction);
  return err;    
}


// *****************************************************************************
//                              
// *****************************************************************************

  ///remove valid contributions and set the states correspondingly
int SumBundleHandler::remove_contributions(FunctionTask ft)
{
  if ((!sumbundle->has_bundle_for(ft))
      ||(sumbundle->get_mode(ft)!=SumBundle::child)){
    return 0;
  }

  if (sumbundle->get_n_contributors(ft)==0){
    sumbundle->set_mode(ft)=SumBundle::inactive;
    return 0;
  }
  
  // let the parent handler make sure that the parent bundle is the root or inactive by maybe removing it first
  assert(parent_handler);
  if (parent_handler->remove_contributions(ft)){
    if (cb_out()){
      get_out()<<"**** ERROR SumBundleHandler::remove_contributions(..):  parent_handler->remove_contributions failed"<<std::endl;
    }
    return 1;
  } 

  assert(bundleinfo[ft]);
  BundleInformation* bi=bundleinfo[ft];
  assert(bi->bundle_size==sumbundle->bundle_size(ft));
  
  assert(parent_handler->bundleinfo[ft]);
  assert(parent_handler->bundleinfo[ft]->bundle_size==parent_handler->sumbundle->bundle_size(ft));
  assert(bi->bundle_size==parent_handler->bundleinfo[ft]->bundle_size);

  Real factor=-sumbundle->get_function_factor(ft)/parent_handler->sumbundle->get_function_factor(ft);

  int err=0;
  if (aft) {
    err=aft->transform_minorants(parent_handler->sumbundle->set_bundle(ft),
			   sumbundle->get_bundle(ft),
			   factor);
    if (err) {
      if (cb_out()){
	get_out()<<"*** ERROR: SumBundleHandler::remove_contributions("<<ft<<"): aft->transform_minorants failed and returned "<<err<<std::endl;
      }
    }
  }
  else {
    for (unsigned int i=0;i<unsigned(bi->bundle_size);i++){
      int lerr=sumbundle->get_bundle(ft)[i].get_minorant(parent_handler->sumbundle->set_bundle(ft)[i],factor);
      if (lerr) {
	if (cb_out()){
	  get_out()<<"*** ERROR: SumBundleHandler::remove_contributions("<<ft<<"): get_minorant failed for minorant "<<i<<" and returned "<<lerr<<std::endl;
	}
	err+=lerr;
      }
    }
  }

  assert(parent_handler->sumbundle->get_n_contributors(ft)>=1);
  parent_handler->sumbundle->set_n_contributors(ft)--;
  parent_handler->sumbundle->set_aggregate(ft).clear();
  sumbundle->set_mode(ft)=SumBundle::root;

  return err;
}

// *****************************************************************************
//                              
// *****************************************************************************

  ///remove valid contributions and set the states correspondingly
int SumBundleHandler::remove_contributions()
{
  int err=0;
  err+=remove_contributions(ObjectiveFunction);
  err+=remove_contributions(ConstantPenaltyFunction);
  err+=remove_contributions(AdaptivePenaltyFunction);
  return err;  
}

// *****************************************************************************
//                              
// *****************************************************************************

/// read the update rules of the parent handler bh so that contributing to it is possible
  int SumBundleHandler::align_bundle(bool use_parent_update,FunctionTask ft)
{
  if (!sumbundle->has_bundle_for(ft))
    return 0;

  BundleInformation* bi=bundleinfo[ft];
  assert(bi);
  bi->map_to_old.init(0,0,Integer(0));
  bi->model_indices.init(0,0,Integer(0));

  if ((parent_handler==0)||(parent_handler->bundleinfo[ft]==0)){
    //nothing to read, nothing to do
    if (sumbundle->get_mode(ft)==SumBundle::child)
      sumbundle->set_mode(ft)=SumBundle::root;
    return 0;
  }

  assert((bi->new_index>=bi->bundle_size)||(! sumbundle->get_bundle(ft)[unsigned(bi->new_index)].empty()));

  const BundleInformation* bhi=parent_handler->bundleinfo[ft];
  assert(bhi);

  if (sumbundle->get_n_contributors(ft)>0){
    assert(bi->bundle_size==sumbundle->bundle_size(ft));
    assert(bi->bundle_size==sumbundle->get_coeff(ft).dim());
    
    if (store_aggregate(ft)){
      if (cb_out()){
	get_out()<<"**** ERROR SumBundleHandler::align_bundle(..):  store_aggregate("<<ft<<") failed"<<std::endl;
      }
      return 1;
    }
    
    if((bi->bundle_size!=bhi->bundle_size)||(use_parent_update && (bhi->map_to_old.dim()!=0))){
      //resize the bundle to match that of bh
      if ((bhi->bundle_size>0)&&(use_parent_update)){
	assert(max(bhi->map_to_old)<bi->bundle_size);
	assert(bhi->map_to_old(bhi->aggr_index)==bi->aggr_index);
	bi->map_to_old=bhi->map_to_old;
      }
      else {
	bi->map_to_old.init(Range(0,bhi->bundle_size-1));
	if (bhi->bundle_size<bi->bundle_size){
	  if (bhi->bundle_size>0){
	    //move aggregate to right position and delete columns exceeding bhi->bundle_size
	    bi->map_to_old(bhi->aggr_index)=bi->aggr_index;
	    if (bi->aggr_index<bhi->bundle_size){
	      bi->map_to_old(bi->aggr_index)=bhi->aggr_index;
	    }
	  }
	}
	else {
	  //bhi->bundelesize>bi->bundle_size, enlarge by copying the aggregate
	  for (Integer i=bi->bundle_size;i<bhi->bundle_size;i++){
	    bi->map_to_old(i)=bi->aggr_index;
	  }
	  bi->map_to_old(bhi->aggr_index)=bi->aggr_index;
	  if (bhi->aggr_index<bi->bundle_size)
	    bi->map_to_old(bi->aggr_index)=bhi->aggr_index;
	}
      }
      if (reassign_minorants(bi->map_to_old,ft)){
	if (cb_out()){
	  get_out()<<"**** ERROR SumBundleHandler::align_bundle(..):  reassign_minorants failed"<<std::endl;
	}
	return 1;
      }
    }
  }

  bi->bundle_size=bhi->bundle_size;
  bi->aggr_index=bhi->aggr_index;
  bi->new_index=bhi->new_index;
  bi->sbp->set_max_bundle_size(max(bi->sbp->get_max_bundle_size(),bhi->sbp->get_max_bundle_size()));
  assert(bi->sbp->get_max_model_size()<=bi->sbp->get_max_bundle_size());
  bi->map_to_old=bhi->map_to_old;

  return 0;
}
      

// *****************************************************************************
//                              
// *****************************************************************************

int SumBundleHandler::align_bundle(bool use_parent_update)
{
  int err=0;
  err+=align_bundle(use_parent_update,ObjectiveFunction);
  err+=align_bundle(use_parent_update,ConstantPenaltyFunction);
  err+=align_bundle(use_parent_update,AdaptivePenaltyFunction);
  return err;  
}

// *****************************************************************************
//                              
// *****************************************************************************

  ///add root contributions to the parent and set the states correspondingly
  int SumBundleHandler::add_contributions(FunctionTask ft)
{
  if ((!sumbundle->has_bundle_for(ft))||
      (sumbundle->get_mode(ft)!=SumBundle::root)){
    return 0;
  }

  if (sumbundle->get_n_contributors(ft)==0){
    sumbundle->set_mode(ft)=SumBundle::inactive;
    return 0;
  }

  if ((parent_handler==0)||(parent_handler->bundleinfo[ft]==0)){
    return 0;
  }
  
  BundleInformation* bi=bundleinfo[ft];
  assert(bi);
  assert(bi->bundle_size==sumbundle->bundle_size(ft));
  assert((bi->new_index>=bi->bundle_size)||(! sumbundle->get_bundle(ft)[unsigned(bi->new_index)].empty()));

  bi->old_diagonal.init(0,0,0.);
  bi->old_lowrank.init(0,0,0.);
  bi->old_sym.init(0,0.);
  
  BundleInformation* bhi=parent_handler->bundleinfo[ft];
  assert(bi->bundle_size==bhi->bundle_size);
  assert(bi->aggr_index==bhi->aggr_index);

  if (parent_handler->sumbundle->get_mode(ft)==SumBundle::child){
    //unlink the parent's sumbundle from its parent's sumbundle
    if (parent_handler->remove_contributions(ft)){
      if (cb_out()){
	get_out()<<"*** ERROR: SumBundleHandler::add_contributions("<<ft<<"): parent_handler->remove_contributions("<<ft<<") failed"<<std::endl;
      }
      return 1;
    }
  }

  //the aggregate of this bundle should have the only nonzero coefficient
  //and this coefficient should be at its upper bound
  assert((bi->bundle_size==0)||(std::fabs(sumbundle->get_coeff(ft)(bi->aggr_index)-sumbundle->get_function_factor(ft))<1e-8*sumbundle->get_function_factor(ft)));

  if (parent_handler->sumbundle->get_n_contributors(ft)==0){
    //initialize this sumbundle
    parent_handler->sumbundle->init(ft,parent_handler->sumbundle->get_function_factor(ft));
    parent_handler->sumbundle->set_coeff(ft).init(bhi->bundle_size,1,0.);
    if (bhi->bundle_size>0){
      parent_handler->sumbundle->set_coeff(ft)(bhi->aggr_index)=parent_handler->sumbundle->get_function_factor(ft);
    }
    parent_handler->sumbundle->set_bundle(ft).clear();
    parent_handler->sumbundle->set_bundle(ft).resize(unsigned(bhi->bundle_size));
  }

  assert(bhi->bundle_size==parent_handler->sumbundle->bundle_size(ft));

  Real factor=sumbundle->get_function_factor(ft)/parent_handler->sumbundle->get_function_factor(ft);

  int err=0;
  if (aft) {
    err=aft->transform_minorants(parent_handler->sumbundle->set_bundle(ft),
			   sumbundle->get_bundle(ft),
			   factor);
    if (err){
      if (cb_out()){
	get_out()<<"*** ERROR: SumBundleHandler::add_contributions("<<ft<<"): aft->transform_minorants failed"<<std::endl;
      }
    }
  }
  else {
    for (unsigned int i=0;i<unsigned(bi->bundle_size);i++){
      int lerr=sumbundle->get_bundle(ft)[i].get_minorant(parent_handler->sumbundle->set_bundle(ft)[i],factor);
      if (lerr){
	if (cb_out()){
	  get_out()<<"*** ERROR: SumBundleHandler::add_contributions("<<ft<<"): get_minorant failed for minorant "<<i<<" and returned "<<lerr<<std::endl;
	}
	err+=lerr;
      }
    }
  }

  parent_handler->sumbundle->set_n_contributors(ft)++;
  parent_handler->sumbundle->set_aggregate(ft).clear();
  parent_handler->sumbundle->set_mode(ft)=SumBundle::root;
  sumbundle->set_mode(ft)=SumBundle::child;
 
  return err;
}

// *****************************************************************************
//                              
// *****************************************************************************

  ///add valid contributions and set the states correspondingly
int SumBundleHandler::add_contributions()
{
  int err=0;
  err+=add_contributions(ObjectiveFunction);
  err+=add_contributions(ConstantPenaltyFunction);
  err+=add_contributions(AdaptivePenaltyFunction);
  return err;  
}

// *****************************************************************************
//                              
// *****************************************************************************

int SumBundleHandler::start_augmodel(QPModelDataPointer& qpblockp,
				     Integer cand_id,
				     const Matrix& cand_y,
				     const Indexmatrix* /* indices */,
				     FunctionTask ft)
{
  if ((! sumbundle->has_bundle_for(ft))
      ||(sumbundle->get_mode(ft)!=SumBundle::root)
      ||(sumbundle->get_n_contributors(ft)==0)
      ||(sumbundle->bundle_size(ft)==0))
    return 0;
  BundleInformation* bi=bundleinfo[ft];
  assert(bi);
  assert((bi->new_index>=bi->bundle_size)||(! sumbundle->get_bundle(ft)[unsigned(bi->new_index)].empty()));

  //because the contributors might have changed it is worth to redetermine the model
  bi->sbp->set_cbout(this,0);
  BundleIdProx H;
  H.set_weightu(bi->tweightu/sumbundle->get_function_factor(ft));
  bi->sbp->select_model(bi->model_indices,cand_id,cand_y,
			sumbundle->get_bundle(ft),0,
			bi->tmodel_maxviol/sumbundle->get_function_factor(ft),
			H,
			bi->tmodel_update);

  //set the block size and copy the selected bundle to the output 
  Integer lxdim=bi->model_indices.dim();
  Indexmatrix tmpind(0,0,Integer(0));

  if (qpblockp.set_model_data(bi->block)){
    delete bi->block;
    bi->block=qpblockp.generate_conemodel_data();
  }
  assert(bi->block);
  MinorantBundle out_bundle;
  out_bundle.resize(unsigned(lxdim));
  for(Integer i=0;i<lxdim;i++){
    out_bundle[unsigned(i)]=(sumbundle->get_bundle(ft))[unsigned(bi->model_indices(i))];
  }
  Matrix dummymat(0,0,0.);
  bi->block->init(MinorantPointer(),out_bundle,lxdim,tmpind,tmpind,dummymat,dummymat,sumbundle->get_function_factor(ft),ft);
  qpblockp.set_model_data(bi->block);
  assert(qpblockp.get_model_data_ptr()==bi->block);
  bi->increase_factor=1.;
  if (cb_out(2)){
    get_out()<<" SBHxdim="<<lxdim;
  }
  return 0;
}

// *****************************************************************************
//                              
// *****************************************************************************

int SumBundleHandler::start_augmodel(QPModelDataPointer& bp,
				     QPSumModelDataObject& sumblock,
				     Integer cand_id,
				     const Matrix& cand_y,
				     const Indexmatrix* indices)
{
  int err=0;
  bp.clear_model_data_ptr();
  err+=start_augmodel(bp,cand_id,cand_y,indices,ObjectiveFunction);
  sumblock.append(bp.get_model_data_ptr());
  bp.clear_model_data_ptr();
  err+=start_augmodel(bp,cand_id,cand_y,indices,ConstantPenaltyFunction);
  sumblock.append(bp.get_model_data_ptr());
  bp.clear_model_data_ptr();
  err+=start_augmodel(bp,cand_id,cand_y,indices,AdaptivePenaltyFunction);
  sumblock.append(bp.get_model_data_ptr());
  bp.clear_model_data_ptr();
  return err;
}


// *****************************************************************************
//                              
// *****************************************************************************

int SumBundleHandler::make_model_aggregate(FunctionTask ft)
{
  if (! sumbundle->has_bundle_for(ft))
    return 0;

  BundleInformation* bi=bundleinfo[ft];
  assert(bi);

  if (sumbundle->get_mode(ft)==SumBundle::inactive){
    if ((parent_handler)&&(parent_handler->sumbundle->has_bundle_for(ft))){
      //distribute the coefficients
      Real factor=sumbundle->get_function_factor(ft)/parent_handler->sumbundle->get_function_factor(ft);
      sumbundle->set_coeff(ft).init(parent_handler->sumbundle->get_coeff(ft),factor);
    }
    else {
      sumbundle->set_coeff(ft).init(bi->bundle_size,1,0.);
      if (bi->bundle_size>0){
	sumbundle->set_coeff(ft)(bi->aggr_index)=sumbundle->get_function_factor(ft);
      }
    }
    assert((ft!=ObjectiveFunction)||(std::fabs(sum(sumbundle->get_coeff(ft))-sumbundle->get_function_factor(ft))<1e-8*sumbundle->get_function_factor(ft)));
    assert((ft==ObjectiveFunction)||(sum(sumbundle->get_coeff(ft))<(1+1e-8)*sumbundle->get_function_factor(ft)));
    sumbundle->set_aggregate(ft).clear();
    return 0;
  }
  
  //store the coefficients in sumbundle
  if (sumbundle->get_mode(ft)==SumBundle::root) {
    if (sumbundle->get_n_contributors(ft)>0){
      Matrix coeff;
      if (bi->block->get_nncx(coeff)){
	if (cb_out()){
	  get_out()<<"**** ERROR SumBundleHandler::make_model_aggregate(....): retrieving the QP solution failed for block "<<ft<<std::endl;
	}
	return 1;
      }
      sumbundle->set_coeff(ft).init(sumbundle->bundle_size(ft),1,0.);
      sumbundle->set_coeff(ft).subassign(bi->model_indices,coeff);
    }
    else {
      sumbundle->set_coeff(ft).init(bi->bundle_size,1,0.);
      if (bi->bundle_size>0){
	sumbundle->set_coeff(ft)(bi->aggr_index)=sumbundle->get_function_factor(ft);
      }
      sumbundle->set_aggregate(ft).clear();
      return 0;
    }
  }
  else {
    assert((parent_handler!=0)&&(sumbundle->get_mode(ft)==SumBundle::child));
    Real factor=sumbundle->get_function_factor(ft)/parent_handler->sumbundle->get_function_factor(ft);
    sumbundle->set_coeff(ft).init(parent_handler->sumbundle->get_coeff(ft),factor);
  }
  assert((ft!=ObjectiveFunction)||(std::fabs(sum(sumbundle->get_coeff(ft))-sumbundle->get_function_factor(ft))<1e-8*sumbundle->get_function_factor(ft)));
  assert((ft==ObjectiveFunction)||(sum(sumbundle->get_coeff(ft))<(1+1e-8)*sumbundle->get_function_factor(ft)));

  sumbundle->set_aggregate(ft).clear();
  sumbundle->set_aggregate(ft).aggregate(sumbundle->get_bundle(ft),sumbundle->get_coeff(ft));
  
  return 0;
}

// *****************************************************************************
//                              
// *****************************************************************************

int SumBundleHandler::make_model_aggregate(bool& increased, bool factor_fixed)
{
  bool incr=false;
  if (sumbundle->has_bundle_for(AdaptivePenaltyFunction)){
    BundleInformation* bi=bundleinfo[AdaptivePenaltyFunction];
    assert(bi);
    bi->increase_factor=1.;
    if(!factor_fixed){
      if (sumbundle->get_mode(AdaptivePenaltyFunction)==SumBundle::root) {
	if (sumbundle->get_n_contributors(AdaptivePenaltyFunction)>0){
	  Real funfactor=sumbundle->get_function_factor(AdaptivePenaltyFunction);
	  Real prec=0;
	  Real ubmodelval=bi->block->tracedual(&prec);
	  Real traceval=bi->block->evaluate_trace();
	  if (ubmodelval>100.*prec){
	    bi->increase_factor=max(2*funfactor,funfactor+ubmodelval)/funfactor;
	    incr=true;
	  }
	  else if ((ubmodelval<0.01*prec)&&
		   (funfactor>1.)&&
		   (10.*traceval<funfactor)
		   ){
	    bi->increase_factor=max(1.,max(0.01*funfactor,5.*traceval))/funfactor;
	    incr=true;
	  }

	  if (incr) {
	    bi->block->adjust_trace(bi->increase_factor*funfactor);
	    if (cb_out(1)){
	      get_out()<<" increase_factor="<<bi->increase_factor<<" ";
	    }
	  }
	    
	}	  
      }
      else if ((sumbundle->get_mode(AdaptivePenaltyFunction)==SumBundle::child)&&(parent_handler->get_increase_factor()!=1.)){
	bi->increase_factor=parent_handler->get_increase_factor();
	incr=true;
      }
    }
    if (incr) {
      sumbundle->set_function_factor(AdaptivePenaltyFunction)*=bi->increase_factor;
    }
  }
  
  increased |= incr;
 

  int err=0;

  if (increased){
    //no need to compute the aggregate
    return err;
  }

  err+=make_model_aggregate(ObjectiveFunction);
  err+=make_model_aggregate(ConstantPenaltyFunction);
  err+=make_model_aggregate(AdaptivePenaltyFunction);

  return err;
}

// *****************************************************************************
//                              
// *****************************************************************************

int SumBundleHandler::provide_model_aggregate(FunctionTask ft)
{
  if (! sumbundle->has_bundle_for(ft))
    return 0;

  if ((sumbundle->get_mode(ft)==SumBundle::inactive)||
      (sumbundle->get_n_contributors(ft)==0) ||
      (sumbundle->get_aggregate(ft).valid()))
    return 0;

  assert(sumbundle->bundle_size(ft)>0);
  assert(min(sumbundle->get_coeff(ft))>=0.);
  assert(sum(sumbundle->get_coeff(ft))<(1.+1e-8)*sumbundle->get_function_factor(ft));
  assert((ft!=ObjectiveFunction)||(sum(sumbundle->get_coeff(ft))>(1.-1e-8)*sumbundle->get_function_factor(ft)));

  
  sumbundle->set_aggregate(ft).aggregate(sumbundle->get_bundle(ft),sumbundle->get_coeff(ft));

  return 0;
}

// *****************************************************************************
//                              
// *****************************************************************************

int SumBundleHandler::provide_model_aggregate()
{
  int err=0;
  err+=provide_model_aggregate(ObjectiveFunction);
  err+=provide_model_aggregate(ConstantPenaltyFunction);
  err+=provide_model_aggregate(AdaptivePenaltyFunction);
  return err;
}

// *****************************************************************************
//                              
// *****************************************************************************

/// see SumBlockModel::adjust_multiplier
int SumBundleHandler::adjust_multiplier(bool& values_may_have_changed)
{
  bool changed=false;
  if (sumbundle->has_bundle_for(AdaptivePenaltyFunction)){
    BundleInformation* bi=bundleinfo[AdaptivePenaltyFunction];
    assert(bi);
    bi->increase_factor=1.;
    if (sumbundle->get_mode(AdaptivePenaltyFunction)==SumBundle::root) {        
      if (sumbundle->get_n_contributors(AdaptivePenaltyFunction)>0){
	Matrix tmpvec;
	if (bi->block->get_nncx(tmpvec)){
	  if (cb_out()){
	    get_out()<<"**** ERROR SumBundleHandler::adjust_multplier(....): retrieving the QP solution failed for AdaptivePenaltyFunction block"<<std::endl;
	  }
	  return 1;
	}
	Real oldval=sumbundle->get_function_factor(AdaptivePenaltyFunction);
	sumbundle->set_function_factor(AdaptivePenaltyFunction)=min(oldval,max(1.,1.5*sum(tmpvec)));
	if (sumbundle->get_function_factor(AdaptivePenaltyFunction)<oldval){
	  changed=true;
	  bi->increase_factor=sumbundle->get_function_factor(AdaptivePenaltyFunction)/oldval;
	  bi->block->adjust_trace(sumbundle->get_function_factor(AdaptivePenaltyFunction));
	}
      }
    }
    else if (parent_handler->get_increase_factor()!=1.){
      bi->increase_factor=parent_handler->get_increase_factor();
      changed=true;
    }
    if (changed) {
      sumbundle->set_function_factor(AdaptivePenaltyFunction)*=bi->increase_factor;
      sumbundle->set_aggregate(AdaptivePenaltyFunction).clear();
    }
  }
  values_may_have_changed |= changed;

  return 0;
}

// *****************************************************************************
//                              
// *****************************************************************************


int SumBundleHandler::contribute_initial_bundle(FunctionTask ft,
						const MinorantBundle& bundle,
						const Matrix& coeff)
{
  if (!sumbundle->has_bundle_for(ft)){
    if (cb_out()){
      get_out()<<"**** ERROR SumBundleHandler::contribute_initial_bundle(.....): no bundle for FunctionTask "<<ft<<std::endl;
    }
    return 1;
  }
  if (sumbundle->get_mode(ft)==SumBundle::child){
    if (remove_contributions(ft)){
      if (cb_out()){
	get_out()<<"**** ERROR SumBundleHandler::contribute_initial_bundle(.....): remove_contribution failed for FunctionTask "<<ft<<std::endl;
      }
      return 1;
    }
  }

  assert(Integer(bundle.size())==coeff.dim());
  assert(min(coeff)>=0.);
  assert(sum(coeff)<=(1.+1e-8)*sumbundle->get_function_factor(ft));

  sumbundle->init(ft,sumbundle->get_function_factor(ft));
  BundleInformation* bi=bundleinfo[ft];
  assert(bi);
  
  sumbundle->set_bundle(ft)=bundle;
  sumbundle->set_coeff(ft).init(coeff);
  sumbundle->set_n_contributors(ft)=1;
  bi->aggr_index=0;
  bi->bundle_size=sumbundle->bundle_size(ft);
  store_aggregate(ft);
  if (bi->bundle_size>bi->sbp->get_max_bundle_size()){
    Indexmatrix map_to_old(Range(0,bi->sbp->get_max_bundle_size()-1));
    reassign_minorants(map_to_old,ft);
  } 

  return 0;
}

// *****************************************************************************
//                              
// *****************************************************************************


int SumBundleHandler::install_external_aggregate(FunctionTask ft,
						 const MinorantPointer& mp,
						 Real aggr_coeff)
{
  if (!sumbundle->has_bundle_for(ft)){
    if (cb_out()){
      get_out()<<"**** ERROR SumBundleHandler::install_external_aggregate(.....): no bundle for FunctionTask "<<ft<<std::endl;
    }
    return 1;
  }
  if (sumbundle->get_mode(ft)!=SumBundle::inactive){
    if (cb_out()){
      get_out()<<"**** ERROR SumBundleHandler::install_external_aggregate(.....): may only be called for FunctionTask "<<ft<<" if inactive, but it is active" <<std::endl;
    }
    return 1;
  }
  if (sumbundle->get_n_contributors(ft)==0){
    if (cb_out()){
      get_out()<<"**** ERROR SumBundleHandler::install_external_aggregate(.....): may not be called for FunctionTask "<<ft<<" whithout contributors but has none" <<std::endl;
    }
    return 1;
  }
  assert(aggr_coeff>=0.);
  assert((ft!=ObjectiveFunction)||(std::fabs(aggr_coeff-sumbundle->get_function_factor(ft))<1e-8*sumbundle->get_function_factor(ft)));
  assert((ft==ObjectiveFunction)||(aggr_coeff<=sumbundle->get_function_factor(ft)));

  BundleInformation* bi=bundleinfo[ft];
  assert(bi);

  if (bi->bundle_size==0){
    sumbundle->set_bundle(ft).resize(1,mp);
    sumbundle->set_coeff(ft).init(1,1,aggr_coeff);
    bi->aggr_index=0;
    bi->bundle_size=sumbundle->bundle_size(ft);
  }
  else {
    assert((bi->aggr_index>=0)&&(bi->aggr_index<bi->bundle_size));
    sumbundle->set_bundle(ft)[unsigned(bi->aggr_index)]=mp;
    sumbundle->set_coeff(ft).init(bi->bundle_size,1,0.);
    sumbundle->set_coeff(ft)(bi->aggr_index)=aggr_coeff;
  }
  sumbundle->set_aggregate(ft).clear();

  store_aggregate(ft);

  return 0;
}

// *****************************************************************************
//                              
// *****************************************************************************

  int SumBundleHandler::set_cand_minorant(FunctionTask ft,
					  const MinorantPointer& mp)
{
  if (!sumbundle->has_bundle_for(ft)){
    if (cb_out()){
      get_out()<<"**** ERROR SumBundleHandler::add_new_minorant(....): no bundle for FunctionTask "<<ft<<std::endl;
    }
    return 1;
  }
  if (sumbundle->get_n_contributors(ft)==0){
    if (cb_out()){
      get_out()<<"**** ERROR SumBundleHandler::add_new_minorant(....): no contributions announced for FunctionTask "<<ft<<std::endl;
    }
    return 1;
  }

  sumbundle->set_cand_minorant(ft)=mp;

  return 0;
}
  
 
// *****************************************************************************
//                              
// *****************************************************************************

void SumBundleHandler::clear_model(FunctionTask ft)
{
  if ((! sumbundle->has_bundle_for(ft))||(sumbundle->bundle_size(ft)==0))
    return;

  if (remove_contributions(ft)){
    if (cb_out()){
      get_out()<<"**** ERROR SumBundleHandler::clear_model(.): remove_contribution failed for FunctionTask"<<ft<<std::endl;
    }
    return;
  }

  sumbundle->init(ft,sumbundle->get_function_factor(ft));
  BundleInformation* bi=bundleinfo[ft];
  assert(bi);
  bi->bundle_size=0;
  bi->new_index=0;
  bi->aggr_index=0;
  bi->increase_factor=1.;
  if (bi->block)
    bi->block->clear();
  bi->model_indices.init(0,0,0.);
  bi->old_diagonal.init(0,0,0.);
  bi->old_lowrank.init(0,0,0.);
  bi->old_sym.init(0,0.);
}

// *****************************************************************************
//                              
// *****************************************************************************

void SumBundleHandler::clear_model()
{
  clear_model(ObjectiveFunction);
  clear_model(ConstantPenaltyFunction);
  clear_model(AdaptivePenaltyFunction);  
}

// *****************************************************************************
//                              
// *****************************************************************************

void SumBundleHandler::clear_aggregates(FunctionTask ft)
{
  if (! sumbundle->has_bundle_for(ft))
    return;

  sumbundle->set_aggregate(ft).clear();

  
  BundleInformation* bi=bundleinfo[ft];
  assert(bi);
  bi->model_indices.init(0,0,0.);
  bool no_aggrs=true;
  for (Integer i=0;(no_aggrs)&&(i<bi->bundle_size);i++){
    if (sumbundle->get_bundle(ft)[unsigned(i)].aggregate()){
      no_aggrs=false;
    }
  }
 
  if (remove_contributions(ft)){
    if (cb_out()){
      get_out()<<"**** ERROR SumBundleHandler::clear_aggregates(.): remove_contribution failed for FunctionTask"<<ft<<std::endl;
    }
    return;
  }

  if (sumbundle->get_n_contributors(ft)==0){
    clear_model(ft);
    return;
  }
    
  Indexmatrix map_to_old(bi->bundle_size,1,Integer(0));
  Integer cnt=0;
  Real sumcoeff=0.; //collects the values of the deleted nonzero coefficients
  Integer maxi=-1; //for searching the largest remaining coefficient
  Real maxv=-1.;
  for (Integer i=0;i<bi->bundle_size;i++){
    if (sumbundle->get_bundle(ft)[unsigned(i)].aggregate()){
      sumcoeff+=sumbundle->get_coeff(ft)(i);
      sumbundle->set_coeff(ft)(i)=0.;
      continue;
    }
    map_to_old(cnt)=i;
    if (sumbundle->get_coeff(ft)(i)>maxv){
      maxv=sumbundle->get_coeff(ft)(i);
      maxi=i;
    }
    cnt++;
  }
  map_to_old.reduce_length(cnt);
  if (maxi>=0) {
    //add the deleted coeff values to a remainig sugradient having maximum coefficient
    sumbundle->set_coeff(ft)(maxi)+=sumcoeff;
  }

  if (reassign_minorants(map_to_old,ft)){
    if (cb_out()){
      get_out()<<"**** ERROR SumBundleHandler::clear_aggregates(.): remove_contribution failed for FunctionTask"<<ft<<std::endl;
    }
  }

  if (
      (map_to_old(bi->aggr_index)!=bi->aggr_index) ||
      (CH_Matrix_Classes::abs(sum(sumbundle->get_coeff(ft))-sumbundle->get_function_factor(ft))>1e-8*sumbundle->get_function_factor(ft))
      ){
    bi->aggr_index=0;
    if (bi->bundle_size>0){
      sumbundle->set_coeff(ft).init(bi->bundle_size,1,0.);
      sumbundle->set_coeff(ft)(bi->aggr_index)=sumbundle->get_function_factor(ft);
    }
  }
  bi->new_index=0;
}

// *****************************************************************************
//                              
// *****************************************************************************

void SumBundleHandler::clear_aggregates()
{
  clear_aggregates(ObjectiveFunction);
  clear_aggregates(ConstantPenaltyFunction);
  clear_aggregates(AdaptivePenaltyFunction);  
}

// *****************************************************************************
//                              
// *****************************************************************************

void SumBundleHandler::clear_cand_minorants()
{
  if (sumbundle->has_bundle_for(ObjectiveFunction))
      sumbundle->set_cand_minorant(ObjectiveFunction).clear();
  if (sumbundle->has_bundle_for(ConstantPenaltyFunction))
      sumbundle->set_cand_minorant(ConstantPenaltyFunction).clear();
  if (sumbundle->has_bundle_for(AdaptivePenaltyFunction))
      sumbundle->set_cand_minorant(AdaptivePenaltyFunction).clear();
}

// *****************************************************************************
//                        add_variable_metric      
// *****************************************************************************

int SumBundleHandler::add_variable_metric(FunctionTask ft,
					  VariableMetric& H,
					  Integer y_id ,
					  const Matrix& y,
					  bool descent_step,
					  Real weightu,
					  Real model_maxviol,
					  const Indexmatrix* indices) 
{
  if (!H.employ_variable_metric())
    return 0;

  if ((! sumbundle->has_bundle_for(ft))||
      (!(sumbundle->get_mode(ft)==SumBundle::root))||
      (sumbundle->get_n_contributors(ft)==0)||
      (sumbundle->get_aggregate(ft).empty())||
      (sumbundle->bundle_size(ft)<=bundleinfo[ft]->model_indices.dim())){
    return 0;
  }

  class VariableMetricSumBundleData: public VariableMetricBundleData
  {
  public:
    BundleInformation* bi;
    SumBundle* sumbundle;
    FunctionTask ft;

    Real get_function_factor() const
    { return sumbundle->get_function_factor(ft);}

    int get_latest_minorants(MinorantBundle& latest_minorants,
			     Integer max_number)
    {
      assert(bi->sbp);
      bi->sbp->set_max_bundle_size(max(bi->sbp->get_max_bundle_size(),max_number));
      if (bi->bundle_size<=max_number)
	latest_minorants=sumbundle->get_bundle(ft);
      else if (bi->bundle_size>0){
	latest_minorants.resize(unsigned(max_number));
	long unsigned int ind=unsigned(bi->new_index);
	for(unsigned int i=0;i<latest_minorants.size();i++){
	  if (ind==0)
	    ind=sumbundle->get_bundle(ft).size();
	  --ind;
	  latest_minorants[i]=(sumbundle->get_bundle(ft))[ind];
	}
      }
      else {
	latest_minorants.clear();
      }
      return 0;
    }

    int get_model_data(MinorantBundle& model_minorants,Matrix& model_coeff) const
    {
      Integer msz=bi->model_indices.dim();
      model_minorants.resize(unsigned(msz));
      for(Integer i=0;i<msz;i++){
	model_minorants[unsigned(i)]=(sumbundle->get_bundle(ft))[unsigned(bi->model_indices(i))];
      }
      model_coeff=sumbundle->get_coeff(ft);
      assert(model_coeff.dim()==msz);
      return 0;
    }

    const MinorantPointer& get_aggregate() const
    { return sumbundle->get_aggregate(ft);}

    Symmatrix& set_denseH()
    { return bi->old_sym; }

    Matrix& set_lowrankH()
    { return bi->old_lowrank; }

    Matrix& set_diagH()
    { return bi->old_diagonal; }

    const Symmatrix& get_denseH() const
    { return bi->old_sym; }

    const Matrix& set_lowrankH() const
    { return bi->old_lowrank; }

    const Matrix& set_diagH() const
    { return bi->old_diagonal; }
  };

  BundleInformation* bi=bundleinfo[ft];
  assert(bi);


  VariableMetricSelection* vm_selection=bi->sbp->get_variable_metric_selection();
  if (vm_selection==0)
    vm_selection=H.get_variable_metric_selection();
  if (vm_selection!=0){
    VariableMetricSumBundleData vmsbd;
    vmsbd.bi=bi;
    vmsbd.sumbundle=sumbundle;
    vmsbd.ft=ft;

    if (vm_selection->add_variable_metric(H,y_id,y,descent_step,
					  weightu,model_maxviol,
					  indices,vmsbd)){
      if (cb_out())
	get_out()<<"**** WARNING: SumBundleHandler::add_variable_metric(ft,.........): H.add_variable_metric(.........) failed"<<std::endl;
      return 1;
    }
  }
  
  return 0;
}

// *****************************************************************************
//                        add_variable_metric      
// *****************************************************************************

int SumBundleHandler::add_variable_metric(VariableMetric& H,
					  Integer yid,
					  const Matrix& y,
					  bool descent_step,
					  Real weightu,
					  Real model_maxviol,
					  const Indexmatrix* indices) 
{
  //currently this is not a good strategy in combination with frequent switching
  //so it is turned off and no old_minorants are collected
  if (!H.employ_variable_metric()) 
    return 0;

  int err=0;
  err+=add_variable_metric(ObjectiveFunction,H,yid,y,descent_step,weightu,model_maxviol,indices);
  err+=add_variable_metric(ConstantPenaltyFunction,H,yid,y,descent_step,weightu,model_maxviol,indices);
  err+=add_variable_metric(AdaptivePenaltyFunction,H,yid,y,descent_step,weightu,model_maxviol,indices);
  return err;  
}

// *****************************************************************************
//                        guess_curvature      
// *****************************************************************************

Real SumBundleHandler::guess_curvature(const MinorantBundle& mnrts,
				       const Indexmatrix& selected_indices,
				       Integer cand_id,
				       const Matrix& cand_y,
				       Real model_maxviol) const
{
  assert(model_maxviol>0.);
  
  Integer n=selected_indices.dim();
  if (n<2)
    return 0.;

  //-- compute function values 
  Matrix val(n,1); chk_set_init(val,1);
  for (Integer i=0;i<selected_indices.dim();i++){
    assert(selected_indices(i)<=Integer(mnrts.size()));
    const MinorantPointer& mnrti(mnrts[unsigned(selected_indices(i))]);
    val(i)=mnrti.evaluate(cand_id,cand_y);
  }
  
  //-- sort them nonincreasingly and determine cutt off bounds
  Indexmatrix sind;
  sortindex(val,sind,false);  //largest index first
  Real minval=val(sind(0))-5.*model_maxviol; //value considered too small for curvature

  //-- determine smallest pairwise curvature; skip reappearing minorants
  Matrix normsqr(n,1); chk_set_init(normsqr,1); 
  Matrix mincurve(n,1,0.);
  Indexmatrix skip(n,1,Integer(0));
  Integer cnt=0;
  Integer stopi=n;
  normsqr(0)=mnrts[unsigned(selected_indices(sind(0)))].dual_norm_squared();
  for (Integer i=1;i<n;i++){
    Integer indi=sind(i);
    Real vali=val(indi);
    if (vali<minval) {//value considered too small for curvature
      break;
      stopi=i;
    }
    
    const MinorantPointer& mnrti(mnrts[unsigned(selected_indices(indi))]);
    Real normsqri=normsqr(i)=mnrti.dual_norm_squared();
    Real mincurvei=max_Real;
    for(Integer j=i;--j>=0;){
      if (skip(j))
	continue;
      Integer indj=sind(j);
      Real ipval=mnrts[unsigned(selected_indices(indj))].ip(mnrti);   
      Real diffnsqr=normsqri-2.*ipval+normsqr(j);
      assert(diffnsqr>=-1e-8*max(1.,normsqri));
      diffnsqr=max(0.,diffnsqr);
      Real valj=val(indj);
      //check whether the minorants are identical and skip i if so
      if ((valj-vali<1e-8*(1.+std::fabs(valj)))
	  &&(diffnsqr<1e-8*(1+normsqri))
	  &&(mnrts[unsigned(selected_indices(j))].equals(mnrti))){
	skip(i)=1;
	break;
      }
      //determine curvature and check for improvement on previous best
      Real curvei=.5*diffnsqr/(valj+model_maxviol-vali);
      assert(curvei>=0.);
      if (curvei<mincurvei)
	mincurvei=curvei;
    }
    if (skip(i)==0){
      cnt++;
      mincurve(i)=mincurvei;
    }
  }

  Real avgcurve=0.;
  if (cnt>0){
    avgcurve=sum(mincurve)/Real(cnt);
  }
  
  if (cb_out(2))
    get_out()<<" {"<<n<<","<<stopi<<","<<cnt<<","<<avgcurve<<","<<model_maxviol<<"}";
  
  return avgcurve;
}

}  //end namespace ConicBundle;
