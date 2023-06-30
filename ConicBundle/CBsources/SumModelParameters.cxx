/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/SumModelParameters.cxx
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



#include "SumModelParameters.hxx"

using namespace CH_Matrix_Classes;

namespace ConicBundle {

  SumModelParametersObject::~SumModelParametersObject()
  {}

  SumModelParameters::~SumModelParameters()
  {}

// *****************************************************************************
//                             select_models
// *****************************************************************************

int SumModelParameters::select_models(SumModel::ModelMap& modelmap)
{  
  Integer n_models=Integer(modelmap.size());
  if (n_models==0) //nothing to do
    return 0;

  //---  perpare the relevant decision data 
  Integer at_most_k_local_models=max_local_models;
  Real maxdev=0.;
  Real sumdev=0.;
  for (SumModel::ModelMap::iterator it=modelmap.begin();it!=modelmap.end();it++){
    Real devval=it->second->get_model_deviation();
    sumdev+=devval;
    maxdev=max(maxdev,devval);
  }
  if (maxdev<eps_Real){
    maxdev=1.;
    sumdev=1.;
  }
  Matrix model_dev(n_models,1,0.);
  Real maxdc=0.;
  Matrix dev_contrib(n_models,1,0.);
  Real maxwc=0.;
  Matrix weighted_curvature(n_models,1,0.);
  std::vector<SumBundle::Mode*> model_mode((unsigned)n_models,0);
  Integer i=0;
  for (SumModel::ModelMap::iterator it=modelmap.begin();it!=modelmap.end();it++,i++){
    Real devval=it->second->get_model_deviation();
    model_dev(i)=devval;
    Real dc=it->second->set_deviation_contribution();
    dc=0.9*dc+0.1*model_dev(i)/sumdev;
    it->second->set_deviation_contribution()=dc;
    dev_contrib(i)=dc;
    maxdc=max(maxdc,dc);
    Real wc=.9*it->second->set_weighted_curvature();
    wc+=.1*it->second->get_model_curvature();
    it->second->set_weighted_curvature()=wc;
    weighted_curvature(i)=wc;
    maxwc=max(maxwc,wc);
    model_mode[unsigned(i)]=&(it->second->set_suggested_mode());
  }
  Real avgval=min(maxdev,sumdev/Real(n_models));
  if (cb_out(2)){
    get_out()<<" modeldev="<<sumdev<<" avgval="<<avgval;
  }

  //--- decide on which submodels should have their own model
  
  SumBundle::Mode submode=SumBundle::inactive;
  Integer strategy=update_rule;
  if (strategy>=10){
    submode=SumBundle::root;
    strategy-=10;
  }
  
  switch(strategy){
    
  case 0: default: {
    // select at_most_k_local_models sorted by current absolute deviation
    // but skip all those whose deviation is almost zero
    Indexmatrix sind;
    if (at_most_k_local_models>0){
      sortindex(model_dev,sind,false);
    }
    else {
      sind.init(Range(0,model_dev.rowdim()-1));
    }
    Real threshold=1e-3*maxdev;
    if (at_most_k_local_models>sind.dim())
      threshold=max(threshold,(1+1e-3)*avgval);
    if (cb_out(2)){
      get_out()<<" devthrsh="<<threshold;
    }
    Integer i=0;
    Integer ind;
    while ((i<at_most_k_local_models)&&
	   (i<sind.dim())&&
	   (model_dev((ind=sind(i)))>threshold)){
      if (cb_out(2)){
	get_out()<<" "<<ind<<"("<<model_dev(ind)<<","<<dev_contrib(ind)<<","<<weighted_curvature(ind)<<")";
      }
      *model_mode[unsigned(ind)]=submode;
      i++;
    }
    while(i<sind.dim()){
      *model_mode[unsigned(sind(i))]=SumBundle::child;
      i++;      
    }
    break;
  }
    
  case 1: {
    // select at_most_k_local_models sorted by average deviation contribution
    // but skip all those whose deviation contribution is almost zero
    Indexmatrix sind;
    if (at_most_k_local_models>0){
      sortindex(dev_contrib,sind,false);
    }
    else {
      sind.init(Range(0,model_dev.rowdim()-1));
    }
    Real threshold=max(1e-3*maxdc,0.1/Real(n_models));
    if (at_most_k_local_models>sind.dim())
      threshold=max(threshold,1./Real(n_models));
    if (cb_out(2)){
      get_out()<<" dcthrsh="<<threshold;
    }
    Integer i=0;
    Integer ind;
    while ((i<at_most_k_local_models)&&
	   (i<sind.dim())&&
	   (dev_contrib((ind=sind(i)))>threshold)){
      if (cb_out(2)){
	get_out()<<" "<<ind<<"("<<model_dev(ind)<<","<<dev_contrib(ind)<<","<<weighted_curvature(ind)<<")";
      }
      *model_mode[unsigned(ind)]=submode;
      i++;
    }
    while(i<sind.dim()){
      *model_mode[unsigned(sind(i))]=SumBundle::child;
      i++;      
    }
    break;
  }
    
  case 2: {
    // select at_most_k_local_models sorted by average deviation contribution
    // but skip all those whose deviation contribution is almost zero
    Real avgcurv=max(eps_Real,sum(weighted_curvature)/weighted_curvature.dim());
    if (cb_out(2)){
      get_out()<<" avgcurv="<<avgcurv;
    }
    Indexmatrix sind;
    if (at_most_k_local_models>0){
      sortindex(weighted_curvature,sind,false);
    }
    else {
      sind.init(Range(0,model_dev.rowdim()-1));
    }
    Real threshold=max(1e-6*maxwc,1e-3*avgcurv);
    if (at_most_k_local_models>sind.dim())
      threshold=max(threshold,0.1*avgcurv);
    if (cb_out(2)){
      get_out()<<" wcthrsh="<<threshold;
    }
    Integer i=0;
    Integer ind;
    while ((i<at_most_k_local_models)&&
	   (i<sind.dim())&&
	   (weighted_curvature((ind=sind(i)))>threshold)){
      if (cb_out(2)){
	get_out()<<" "<<ind<<"("<<model_dev(ind)<<","<<dev_contrib(ind)<<","<<weighted_curvature(ind)<<")";
      }
      *model_mode[unsigned(ind)]=submode;
      i++;
    }
    while(i<sind.dim()){
      *model_mode[unsigned(sind(i))]=SumBundle::child;
      i++;      
    }
    break;
  }  

  } //end switch(strategy)

  return 0;
}



}

