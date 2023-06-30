/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/NNCModelParameters.cxx
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



#include "NNCModelParameters.hxx"

using namespace CH_Matrix_Classes;

namespace ConicBundle {


NNCModelParametersObject::~NNCModelParametersObject()
{}
 
NNCModelParameters::~NNCModelParameters()
{}
  
//******************************************************************************
//                              select_model
//******************************************************************************

int NNCModelParameters::select_model(MinorantBundle& model,
				     Matrix& coefficients,
				     Matrix& indicators,
				     const MinorantPointer& aggregate,
				     const MinorantPointer& /* center_minorant */,
				     const MinorantBundle& cand_minorants,
				     const MinorantBundle& /* old_minorants */,
				     MatrixFunctionOracle* /* oracle */,
				     FunctionTask function_task,
				     Real function_factor,
				     BundleModel::ModelUpdate model_update,
				     Integer /* center_id */,
				     const Matrix& /* center_y */,
				     Integer cand_id,
				     const Matrix& cand_y,
				     Real /* model_maxviol */,
				     BundleProxObject& /* H */)
{
  assert(Integer(model.size())==coefficients.dim());
  assert(indicators.dim()==coefficients.dim());

  
  
  switch(update_rule){
  case 0: default:
    return select_active_model(model,coefficients,indicators,cand_id,cand_y,cand_minorants,false);
  case 1:
    return select_active_model(model,coefficients,indicators,cand_id,cand_y,cand_minorants,(model_update==BundleModel::null_step));
  case 2:
    return select_model_heuristic(model,coefficients,indicators,cand_id,cand_y,cand_minorants,aggregate,function_task,function_factor);
  }
  return 0;
}
  
//******************************************************************************
//                              select_active_model
//******************************************************************************

int NNCModelParameters::select_active_model(MinorantBundle& model,
					    Matrix& modelcoeff,
					    Matrix& indicators,
					    Integer y_id,
					    const Matrix& y,
					    const MinorantBundle& minorants,
					    bool grow)
{
  //assert((0<=aggr_index)&&(aggr_index<Integer(minorants.size())));
  assert(Integer(model.size())==modelcoeff.rowdim());
  assert(modelcoeff.rowdim()==indicators.dim());
  assert(min(modelcoeff)>=0.);

  if ((!grow)||((max_model_size>=2)&&(modelcoeff.rowdim()>=max_model_size))){
    
    //determine which coefficient/indicator values should be aggregated
    bool use_indicators=(modelcoeff.dim()==indicators.dim());
    Real sumcoeffs=sum(modelcoeff);
    Real coeffzero=max(eps_Real,sumcoeffs*1e-10);
    Real sumindicators=sum(indicators);
    Real indicatorzero=max(eps_Real,1e-3*sumindicators/indicators.rowdim());
    
    //collect the new bundle and coefficients here and swap afterwards
    MinorantBundle newmodel;
    Matrix newcoeff(modelcoeff.dim()+1,1,0.);
    newcoeff.init(0,1,0.);
    Matrix newindicator(modelcoeff.dim()+1,1,0.);
    newindicator.init(0,1,0.);

    //collect those that should be aggregated here and aggregate afterwards
    MinorantBundle aggrmodel;
    Matrix aggrcoeff(modelcoeff.dim()+1,1,0.);
    aggrcoeff.init(0,1,0.);
    Matrix aggrindicator(modelcoeff.dim()+1,1,0.);
    aggrindicator.init(0,1,0.);

    //go through the old model, add the unaggregated ones with active indicator, aggregate the others if necessary
    for (unsigned ui=0;ui<model.size();ui++){
      if ((model[ui].aggregate())&&
	  ((use_indicators&&(indicators(Integer(ui))>indicatorzero))||(modelcoeff(Integer(ui))>coeffzero))
	  ){
	aggrmodel.push_back(model[ui]);
	aggrcoeff.concat_below(modelcoeff(Integer(ui)));
	aggrindicator.concat_below(indicators(Integer(ui)));
      }
      else if (use_indicators && (indicators(Integer(ui))>indicatorzero)){
	newmodel.push_back(model[ui]);
	newcoeff.concat_below(modelcoeff(Integer(ui)));
	newindicator.concat_below(indicators(Integer(ui)));
      }
      else if (modelcoeff(Integer(ui))>coeffzero){
	aggrmodel.push_back(model[ui]);
	aggrcoeff.concat_below(modelcoeff(Integer(ui)));
	aggrindicator.concat_below(0.);
      }
    }

    //if the size of the newmodel is still too big, sort by coefficient values
    if ((max_model_size>=2)&&(newcoeff.rowdim()>=max_model_size-1-(aggrcoeff.rowdim()>0))){
      Indexmatrix sind;
      sortindex(newcoeff,sind,false);
      model.clear();
      modelcoeff.init(0,1,0.);
      indicators.init(0,1,0.);
      for(Integer i=0;i<max_model_size-2;i++){
	Integer ind=sind(i);
	model.push_back(newmodel[unsigned(ind)]);
	modelcoeff.concat_below(newcoeff(ind));
	indicators.concat_below(newindicator(ind));
      }
      for(Integer i=max_model_size-2;i<newcoeff.rowdim();i++){
	Integer ind=sind(i);
	aggrmodel.push_back(newmodel[unsigned(ind)]);
	aggrcoeff.concat_below(newcoeff(ind));
	aggrindicator.concat_below(0.);
      }
    }
    else {
      swap(model,newmodel);
      swap(modelcoeff,newcoeff);
      swap(indicators,newindicator);
    }
    
    //form a new aggregate (if there are candidates) and append it
    if (aggrmodel.size()>0){
      Real sumaggrcoeff=sum(aggrcoeff);
      MinorantPointer aggr;
      aggr.aggregate(aggrmodel,aggrcoeff,1./sumaggrcoeff);
      model.push_back(aggr);
      modelcoeff.concat_below(sumaggrcoeff);
      sumaggrcoeff=sum(aggrindicator);
      indicators.concat_below(sumaggrcoeff);
    }
  }

  //compute value of the new cutting planes at y and find the maximizer
  Real maxval=min_Real;
  Integer maxind=-1;
  for (unsigned i=0;i<minorants.size();i++){
    Real val=minorants[i].evaluate(y_id,y);;
    if (val>maxval){
      maxval=val;
      maxind=Integer(i);
    }		     
  }
  assert(maxind>=0);

  //if the new maximizing minorant is not yet in the new bundle, add it
  if (maxind>=0){
    MinorantPointer newmnrt=minorants[unsigned(maxind)];
    bool append=true;
    for (unsigned ui=0;ui<model.size();ui++){
      if (newmnrt.equals(model[ui])){
	append=false;
	break;
      }
    }
    if (append){
      model.push_back(newmnrt);
      modelcoeff.concat_below(0.);
      indicators.concat_below(0.);
    }
  }

  //set the new model and the coefficients by swapping

  if (cb_out(2)){
    get_out()<<" select_active "<<model.size()<<" sum(coeff)="<<sum(modelcoeff)<<" sum(indic)="<<sum(indicators)<<std::endl;
    get_out()<<transpose(modelcoeff);
    get_out()<<transpose(indicators);
  }
  
  assert(modelcoeff.dim()==Integer(model.size()));
  
  return 0;
}


//******************************************************************************
//                              select_model_heuristic
//******************************************************************************

int NNCModelParameters::select_model_heuristic(MinorantBundle& model,
					       Matrix& modelcoeff,
					       Matrix& indicators,
					       Integer y_id,
					       const Matrix& y,
					       const MinorantBundle& minorants,
					       const MinorantPointer& aggregate,
					       FunctionTask function_task,
					       Real function_factor)
{
  assert(modelcoeff.dim()==Integer(model.size()));
  assert(modelcoeff.dim()==indicators.dim());
  assert(min(modelcoeff)>=0.);

  //compute value of the new cutting planes at y
  Matrix tmpvec2(Integer(minorants.size()),1); chk_set_init(tmpvec2,1);
  for (unsigned i=0;i<minorants.size();i++){
    tmpvec2(Integer(i))=minorants[i].evaluate(y_id,y);
  }
  
  Real aggrval;
  if (aggregate.valid())
    aggrval=aggregate.evaluate(y_id,y)/function_factor;
  else
    aggrval=sum(tmpvec2)/tmpvec2.rowdim();
  
  //if the function is of penalty type AND the value in the candidate is 
  //negative for all new subgradients, do not change the model at all
  if ((function_task!=ObjectiveFunction)&&(max(tmpvec2)<=0.)){
    if (cb_out(1)){
      get_out()<<" penaltyskip model_update "<<std::endl;
    }
    return 0;
  }

  // determine a maximal size of the model
  Integer max_heu_model_size=max_model_size;
  if (max_heu_model_size<2){
    Real cutoff=max(eps_Real,1e-6*sum(modelcoeff)/modelcoeff.dim());
    Integer cnt=0;
    for(Integer i=0;i<modelcoeff.dim();i++){
      if(modelcoeff(i)>cutoff)
	cnt++;
    }
    max_heu_model_size=min(cnt+5,y.dim());
  }

  //compute value of model cutting planes at y
  Matrix tmpvec(Integer(model.size()),1); chk_set_init(tmpvec,1);
  for (unsigned i=0;i<model.size();i++){
    tmpvec(Integer(i))=model[i].evaluate(y_id,y);
  }
  //append new values
  //(after a null step one of these should have higher value than model)
  tmpvec.concat_below(tmpvec2);
  //sort by this value and use this sorting for generating the next model
  //(we must keep at least one from the model and 
  // try to eliminate identical cuts)
  Indexmatrix tmpind;
  sortindex(tmpvec,tmpind);
  Integer sind=tmpind.dim();
  Indexmatrix modeldel;
  Indexmatrix modelkeep;
  Indexmatrix subgkeep;
  Real maxval=tmpvec(tmpind(sind-1));
  
  while (--sind>=0){
    Integer ind=tmpind(sind);
    bool bunvec=(ind<Integer(model.size()));
    if ((modelkeep.dim()+subgkeep.dim()>=max_heu_model_size)||
        ((!bunvec)&&(subgkeep.dim()>=max_heu_model_size-1))
	){//there are already enough new vectors
      if (bunvec) { 
	modeldel.concat_below(ind);
      }
      continue;
    }
    
    //skip new ones that are worse than aggrval minus a bit
    if ((!bunvec)&&(tmpvec(ind)<aggrval-max(1e-3*(maxval-aggrval),1e-8*(1.+std::fabs(maxval)))))
      continue;
    
    //--- is there a cluster of cutting planes with identical value in cand_y?
    Real lbval=tmpvec(ind);
    Real tol=(fabs(lbval)+1.)*(100.*eps_Real);
    lbval-=tol;
    
    if ((sind==0)||(tmpvec(tmpind(sind-1))<=lbval)){ 
      //no cluster, only one element
      if (bunvec){
	modelkeep.concat_below(ind);
      }
      else {
	subgkeep.concat_below(ind-Integer(model.size()));
      }
      continue;
    }
	
    //--- find the cluster and split it into model vectors and new subgradients
    Integer si=sind;
    Indexmatrix bunind;
    Indexmatrix subgind;
    if (bunvec) bunind.concat_below(ind);
    else subgind.concat_below(ind-Integer(model.size()));
    while ((--si>=0)&&(tmpvec(tmpind(si))>lbval)){
      Integer i=tmpind(si);
      if (i<Integer(model.size())){
	bunind.concat_below(i);
      }
      else {
	subgind.concat_below(i-Integer(model.size()));
      }	   
    }
    sind=si+1;

    //--- first treat cluster of model vectors
    Indexmatrix scoeffind;
    sortindex(modelcoeff(bunind),scoeffind);
    si=scoeffind.dim();
    Integer newbunstart=modelkeep.dim();
    while(--si>=0){
      Integer ind=bunind(scoeffind(si));
      bool deleteind=false;
      if (modelkeep.dim()+subgkeep.dim()>=max_heu_model_size) {
	deleteind=1;
      }
      else {
	for(Integer i=newbunstart;i<modelkeep.dim();i++){
	  unsigned j=unsigned(modelkeep(i)); 
	  if (model[unsigned(ind)].equals(model[j],tol)){
            deleteind=true;
	    break;
	  }
	} //end for
      } //end else
      if (deleteind){
	modeldel.concat_below(ind);
      }
      else {
        modelkeep.concat_below(ind);
      }
    }  
    
    //--- now treat cluster of subgradient vectors
    Integer newsubgstart=subgkeep.dim();
    for(si=0;si<subgind.dim();si++){
      Integer ind=subgind(si);
      bool deleteind=false;
      if ((modelkeep.dim()+subgkeep.dim()>=max_heu_model_size)||
          ((modelkeep.dim()==0)&&(subgkeep.dim()==max_heu_model_size-1))){
	deleteind=true;
      }
      else {
	for(Integer i=newbunstart;i<modelkeep.dim();i++){
	  unsigned j=unsigned(modelkeep(i)); 
	  if (minorants[unsigned(ind)].equals(model[j],tol)){
            deleteind=true;
	    break;
	  }
	} //end for
	{for(Integer i=newsubgstart;(i<subgkeep.dim())&&(!deleteind);i++){
	   unsigned j=unsigned(subgkeep(i)); 
	  if (minorants[unsigned(ind)].equals(minorants[j],tol)){
            deleteind=true;
	    break;
	  }
	}} //end for
      } //end else
      if (!deleteind){
	subgkeep.concat_below(ind);
      }
    }//end for cluster of new subgradients

  }//end while(cutting planes)

  assert(modelkeep.dim()+modeldel.dim()==Integer(model.size()));

  if (cb_out(1)){
    get_out().precision(8);
  get_out()<<" select3: keep="<<modelkeep.dim()<<"("<<max(tmpvec(modelkeep))<<","<<min(tmpvec(modelkeep))<<")";
  get_out()<<" del="<<modeldel.dim()<<"("<<max(tmpvec(modeldel))<<","<<min(tmpvec(modeldel))<<")";
  get_out()<<" new="<<subgkeep.dim()<<"("<<max(tmpvec2(subgkeep))<<","<<min(tmpvec2(subgkeep))<<")";
  get_out()<<" ignore="<<Integer(minorants.size())-subgkeep.dim()<<"("<<min(tmpvec2)<<")"<<std::endl;
  }
    
  //--- aggregate deleted model columns into smallest model column kept
  if (modeldel.dim()>0){
    Integer aggregind;
    Real maxcoeff=max(modelcoeff);
    min(modelcoeff(modelkeep),&aggregind);
    aggregind=modelkeep(aggregind);
    Real s=sum(modelcoeff(modeldel))+modelcoeff(aggregind);
    if (s>max(1.,maxcoeff)*eps_Real) {
      model[unsigned(aggregind)].scale(modelcoeff(aggregind)/s);
      for(Integer i=0;i<modeldel.dim();i++){
	Integer delind=modeldel(i);
	Real a=modelcoeff(delind)/s;
	if (a>max(1.,maxcoeff)*eps_Real) {
	  model[unsigned(aggregind)].aggregate(model[unsigned(delind)],a);
	}
      }
    }
    modelcoeff(aggregind)=s;
    indicators(aggregind)+=sum(indicators(modeldel));
    //modelvecs.delete_cols(modeldel);
    sortindex(modeldel,tmpind);
    Integer cnt=modeldel(tmpind(0));
    Integer delcnt=0;
    Integer nextdel;
    if (++delcnt<modeldel.dim()){
      nextdel=modeldel(tmpind(delcnt));
    }
    else {
      nextdel=Integer(model.size());
    }
    for(unsigned long i=(unsigned long)(cnt)+1;i<model.size();i++){
      if (Integer(i)==nextdel){
	if (++delcnt<modeldel.dim()){
	  nextdel=modeldel(tmpind(delcnt));
	}
	else {
	  nextdel=Integer(model.size());
	}
	continue;
      } 
      model[(unsigned long)(cnt)]=model[i];
      cnt++;
    }
    model.resize((unsigned long)(cnt));
    modelcoeff.delete_rows(modeldel);
    indicators.delete_rows(modeldel);
  }
  
  //--- add new subgradient columns to model
  if (subgkeep.dim()>0){
    if (subgkeep.dim()==Integer(minorants.size())){
      unsigned long olddim=model.size();
      model.resize(olddim+minorants.size());
      for(unsigned long i=0;i<(unsigned long)(subgkeep.dim());i++){
	model[olddim+i]=minorants[i];
      }
      modelcoeff.enlarge_below(Integer(minorants.size()),0.);
      indicators.enlarge_below(Integer(minorants.size()),0.);
    }
    else {
      unsigned long olddim=model.size();
      model.resize(olddim+(unsigned long)(subgkeep.dim()));
      for(unsigned long i=0;i<(unsigned long)(subgkeep.dim());i++){
	model[olddim+i]=minorants[(unsigned long)(subgkeep(Integer(i)))];
      }
      modelcoeff.enlarge_below(subgkeep.dim(),0.);
      indicators.enlarge_below(subgkeep.dim(),0.);
    }    
  }
  
  assert(modelcoeff.dim()==Integer(model.size()));
  assert(modelcoeff.dim()==indicators.dim());
  
  return 0;

}




} // end namespace ConicBundle
