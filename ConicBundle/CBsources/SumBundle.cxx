/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/SumBundle.cxx
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
#include "SumBundle.hxx"

 
using namespace CH_Matrix_Classes;

namespace ConicBundle {

// *****************************************************************************
//                              BData
// *****************************************************************************

SumBundle::BData::BData()
{
  mode=unavailable;
  clear();
}

// *****************************************************************************
//                              BData(const BData&)  
// *****************************************************************************

  SumBundle::BData::BData(const SumBundle::BData& bd):CBout()
{
  init(bd);
}

// *****************************************************************************
//                              ~BData
// *****************************************************************************

SumBundle::BData::~BData()
{
  clear();
}
  
// *****************************************************************************
//                              BData::clear
// *****************************************************************************

void SumBundle::BData::clear(Real fun_factor)
{
  // mode is not influenced by this
  n_contributors=0;
  function_factor=fun_factor;
  
  bundle.clear();
  cand_minorant.clear();

  aggregate.clear();
  coeff.init(0,1,0.);			     

}
  
// *****************************************************************************
//                              BData::init
// *****************************************************************************

void SumBundle::BData::init(const SumBundle::BData& bd)
{
  clear();
  mode=bd.mode;
  n_contributors=bd.n_contributors;
  function_factor=bd.function_factor;

  bundle=bd.bundle;
  aggregate=bd.aggregate;
  coeff.init(bd.coeff);			     

}

// *****************************************************************************
//                              BData::synchronize_ids
// *****************************************************************************

void SumBundle::BData::synchronize_ids(Integer new_modification_id,
				       Integer new_center_id,
				       Integer old_center_id,
				       Integer new_cand_id,
				       Integer old_cand_id,
				       Integer new_prex_id)
{
  aggregate.synchronize_ids(new_modification_id,
			new_center_id,old_center_id,
			new_cand_id,old_cand_id,
			new_prex_id);
  for(unsigned i=0;i<bundle.size();i++){
    bundle[i].synchronize_ids(new_modification_id,
			new_center_id,old_center_id,
			new_cand_id,old_cand_id,
			new_prex_id);
  }

}
  
// *****************************************************************************
//                              BData::operator=
// *****************************************************************************

SumBundle::BData&  SumBundle::BData::operator=(const SumBundle::BData& bd)
{
  this->init(bd);
  return *this;
}


// *****************************************************************************
//                           BData::get_model_aggregate
// *****************************************************************************

int SumBundle::BData::get_model_aggregate(MinorantPointer& aggr,
					  Real factor,
					  const AffineFunctionTransformation* aft) const
{
  assert(aggregate.valid());
  if (aft){
    int retval=aft->transform_minorant(aggr,aggregate,factor,false);
    return retval;
  }

  int retval=aggregate.get_minorant(aggr,factor);
  return retval;
}

// *****************************************************************************
//                             call_primal_extender
// *****************************************************************************

int SumBundle::BData::call_primal_extender(PrimalExtender& prex,Integer prex_id)
{
  int err=0;
  if (aggregate.call_primal_extender(prex,prex_id)){
    if (cb_out()) 
      get_out()<<"**** WARNING: SumBundle::Bdata::call_primal_extender(..): PrimalExtender::extend failed for aggregate"<<std::endl;
    err++;
  }
   
  for(unsigned int i=0;i<bundle.size();i++){
    if (bundle[i].call_primal_extender(prex,prex_id)){
      if (cb_out())
	get_out()<<"**** WARNING: SumBundle::Bdata::call_primal_extender(..): PrimalExtender::extend failed for bundle minorant number "<<i<<std::endl;
      err++;
    }
  }
 

  return err;
}

// *****************************************************************************
//                             BData::apply_modification
// *****************************************************************************

  int SumBundle::BData::apply_modification(const GroundsetModification& gsmdf,
			 CH_Matrix_Classes::Integer mod_id,
			 MinorantExtender* mex)
{
  int err=0;
  if ((aggregate.valid())&&(aggregate.apply_modification(gsmdf,mod_id,mex))){
    if (cb_out()) 
      get_out()<<"**** WARNING: SumBundle::BData::call_primal_extender(..): PrimalExtender::extend failed for aggregate"<<std::endl;
    err++;
  }
   
  for(unsigned int i=0;i<bundle.size();i++){
    if (bundle[i].apply_modification(gsmdf,mod_id,mex)){
      if (cb_out())
	get_out()<<"**** WARNING: SumBundle::BData::call_primal_extender(..): PrimalExtender::extend failed for bundle minorant number "<<i<<std::endl;
      err++;
    }
  }
 

  return err;
}

// *****************************************************************************
//                              SumBundle
// *****************************************************************************

  SumBundle::SumBundle():CBout(),bdata(3)
{
  clear();
}

// *****************************************************************************
//                              SumBundle(const SumBundle&)
// *****************************************************************************

  SumBundle::SumBundle(const SumBundle& sb):CBout(),bdata(3)
{
  init(sb);
}

// *****************************************************************************
//                              ~SumBundle
// *****************************************************************************

SumBundle::~SumBundle()
{
  clear();
  bdata.clear();
}

// *****************************************************************************
//                              clear
// *****************************************************************************

void SumBundle::clear()
{
  clear(ObjectiveFunction);
  clear(ConstantPenaltyFunction);
  clear(AdaptivePenaltyFunction);
}

// *****************************************************************************
//                              clear(FunctionTask)
// *****************************************************************************

void SumBundle::clear(FunctionTask ft)
{
  bdata[ft].clear();
}

// *****************************************************************************
//                                operator=
// *****************************************************************************

SumBundle& SumBundle::operator=(const SumBundle& sb)
{
  init(sb);
  return *this;
}

// *****************************************************************************
//                                init
// *****************************************************************************

void SumBundle::init(const SumBundle& sb)
{
  clear();
  bdata=sb.bdata;
}

// *****************************************************************************
//                                init
// *****************************************************************************

  void SumBundle::init(FunctionTask ft,Real fun_factor)
{
  assert(fun_factor>0.);
  bdata[ft].clear(fun_factor);
}

// *****************************************************************************
//                              synchronize_modification_id
// *****************************************************************************

void SumBundle::synchronize_ids(Integer new_modification_id,
				Integer new_center_id,
				Integer old_center_id,
				Integer new_cand_id,
				Integer old_cand_id,
				Integer new_prex_id)
{
  bdata[ObjectiveFunction].synchronize_ids(new_modification_id,
			new_center_id,old_center_id,
			new_cand_id,old_cand_id,
			new_prex_id);
  bdata[ConstantPenaltyFunction].synchronize_ids(new_modification_id,
			new_center_id,old_center_id,
			new_cand_id,old_cand_id,
			new_prex_id);
  bdata[AdaptivePenaltyFunction].synchronize_ids(new_modification_id,
			new_center_id,old_center_id,
			new_cand_id,old_cand_id,
			new_prex_id);
}

// *****************************************************************************
//                                has_bundle_for
// *****************************************************************************

bool SumBundle::has_bundle_for(FunctionTask ft) const
{
  return bdata[ft].mode!=unavailable;
}
  
// *****************************************************************************
//                                has_bundle_data
// *****************************************************************************

bool SumBundle::has_bundle_data() const
{
  return has_bundle_for(ObjectiveFunction)||
    has_bundle_for(ConstantPenaltyFunction)||
    has_bundle_for(AdaptivePenaltyFunction);
}

// *****************************************************************************
//                                bundle_size
// *****************************************************************************

Integer SumBundle::bundle_size(FunctionTask ft) const
{
  return Integer(bdata[ft].bundle.size());
}

  
// *****************************************************************************
//                                has_roots
// *****************************************************************************

bool SumBundle::has_roots() const
{
  return (bdata[ObjectiveFunction].mode==root)||
    (bdata[ConstantPenaltyFunction].mode==root)||
    (bdata[AdaptivePenaltyFunction].mode==root);
}
  
// *****************************************************************************
//                                has_working_roots
// *****************************************************************************

bool SumBundle::has_working_roots() const
{
  return ((bdata[ObjectiveFunction].mode==root)&&(bdata[ObjectiveFunction].n_contributors>0))||
    ((bdata[ConstantPenaltyFunction].mode==root)&&(bdata[ConstantPenaltyFunction].n_contributors>0))||
    ((bdata[AdaptivePenaltyFunction].mode==root)&&(bdata[AdaptivePenaltyFunction].n_contributors>0));
}
  
// *****************************************************************************
//                                active
// *****************************************************************************

bool SumBundle::active() const
{
  return (bdata[ObjectiveFunction].mode==root)||
    (bdata[ObjectiveFunction].mode==child)||
    (bdata[ConstantPenaltyFunction].mode==root)||
    (bdata[ConstantPenaltyFunction].mode==child)||
    (bdata[AdaptivePenaltyFunction].mode==root)||
    (bdata[AdaptivePenaltyFunction].mode==child);
}
  
// *****************************************************************************
//                              has_contributions
// *****************************************************************************

bool SumBundle::has_contributions() const
{
  return (bdata[ObjectiveFunction].mode==child)||
    (bdata[ConstantPenaltyFunction].mode==child)||
    (bdata[AdaptivePenaltyFunction].mode==child);
}
  
// *****************************************************************************
//                                take_control
// *****************************************************************************

  void SumBundle::take_control(SumBundle::Mode mode) 
{
  for(BDataVector::iterator it=bdata.begin();it!=bdata.end();it++){
    if (it->mode==unavailable)
      continue;
    if (mode!=inactive)
      it->mode=root;
    else
      it->mode=inactive;
  }
}
  

// *****************************************************************************
//                                get_valid
// *****************************************************************************

//bool SumBundle::get_valid(FunctionTask ft) const
//{
//  assert(bdata[ft].mode!=unavailable);
//  return bdata[ft].aggregate.valid();
//}
  
// *****************************************************************************
//                                get_mode
// *****************************************************************************

  SumBundle::Mode SumBundle::get_mode(FunctionTask ft) const
{
  return bdata[ft].mode;
}
  
// *****************************************************************************
//                                set_mode
// *****************************************************************************

  SumBundle::Mode& SumBundle::set_mode(FunctionTask ft) 
{
  return bdata[ft].mode;
}
  
// *****************************************************************************
//                                get_function_factor
// *****************************************************************************

Real SumBundle::get_function_factor(FunctionTask ft) const
{
  assert(bdata[ft].mode!=unavailable);
  return bdata[ft].function_factor;
}
  
// *****************************************************************************
//                                set_function_factor
// *****************************************************************************

Real& SumBundle::set_function_factor(FunctionTask ft) 
{
  assert(bdata[ft].mode!=unavailable);
  return bdata[ft].function_factor;
}
  
// *****************************************************************************
//                                get_n_contributors
// *****************************************************************************

Integer SumBundle::get_n_contributors(FunctionTask ft) const
{
  assert(bdata[ft].mode!=unavailable);
  return bdata[ft].n_contributors;
}
  
// *****************************************************************************
//                                set_n_contributors
// *****************************************************************************

Integer& SumBundle::set_n_contributors(FunctionTask ft) 
{
  assert(bdata[ft].mode!=unavailable);
  return bdata[ft].n_contributors;
}
  
// *****************************************************************************
//                                get_bundle
// *****************************************************************************

const MinorantBundle& SumBundle::get_bundle(FunctionTask ft) const
{
  assert(bdata[ft].mode!=unavailable);
  return bdata[ft].bundle;
}
  
// *****************************************************************************
//                                set_bundle
// *****************************************************************************

MinorantBundle& SumBundle::set_bundle(FunctionTask ft) 
{
  assert(bdata[ft].mode!=unavailable);
  return bdata[ft].bundle;
}
  
  
// *****************************************************************************
//                                get_coeff
// *****************************************************************************

const Matrix& SumBundle::get_coeff(FunctionTask ft) const
{
  assert(bdata[ft].mode!=unavailable);
  return bdata[ft].coeff;
}
  
// *****************************************************************************
//                                set_coeff
// *****************************************************************************

Matrix& SumBundle::set_coeff(FunctionTask ft) 
{
  assert(bdata[ft].mode!=unavailable);
  return bdata[ft].coeff;
}
  
// *****************************************************************************
//                                get_aggregate
// *****************************************************************************

const MinorantPointer& SumBundle::get_aggregate(FunctionTask ft) const
{
  assert(bdata[ft].mode!=unavailable);
  return bdata[ft].aggregate;
}
  
// *****************************************************************************
//                                set_aggregate
// *****************************************************************************

MinorantPointer& SumBundle::set_aggregate(FunctionTask ft) 
{
  assert(bdata[ft].mode!=unavailable);
  return bdata[ft].aggregate;
}
  
// *****************************************************************************
//                                get_cand_minorant
// *****************************************************************************

const MinorantPointer& SumBundle::get_cand_minorant(FunctionTask ft) const
{
  assert(bdata[ft].mode!=unavailable);
  return bdata[ft].cand_minorant;
}
  
// *****************************************************************************
//                                set_cand_minorant
// *****************************************************************************

MinorantPointer& SumBundle::set_cand_minorant(FunctionTask ft) 
{
  assert(bdata[ft].mode!=unavailable);
  return bdata[ft].cand_minorant;
}
  

// *****************************************************************************
//                             get_local_model_aggregate
// *****************************************************************************

int SumBundle::get_local_model_aggregate(MinorantPointer& aggregate,
					 Real factor,
					 const AffineFunctionTransformation* aft) const
{
  int err=0;
  FunctionTask ft=ObjectiveFunction;
  if ((bdata[ft].mode==root)&&
      (bdata[ft].n_contributors>0)){
    if (!bdata[ft].aggregate.valid()){
      err++;
      if (cb_out()) 
	get_out()<<"**** WARNING: SumBundle::get_local_model_aggregate(...): aggregate for ObjectiveFunction part is not valid"<<std::endl;
    }
    else if (bdata[ft].get_model_aggregate(aggregate,factor,aft)){
      err++;
      if (cb_out()) 
	get_out()<<"**** WARNING: SumBundle::get_local_model_aggregate(...): collecting aggregate for ObjectiveFunction part failed"<<std::endl;
    }
  }

  ft=ConstantPenaltyFunction;
  if ((bdata[ft].mode==root)&&
      (bdata[ft].n_contributors>0)){
    if (!bdata[ft].aggregate.valid()){
      err++;
      if (cb_out()) 
	get_out()<<"**** WARNING: SumBundle::get_local_model_aggregate(...): aggregate for ConstantPenaltyFunction part is not valid"<<std::endl;
    }
    else if (bdata[ft].get_model_aggregate(aggregate,factor,aft)){
      err++;
      if (cb_out()) 
	get_out()<<"**** WARNING: SumBundle::get_local_model_aggregate(...): collecting aggregate for ConstantPenaltyFunction part failed"<<std::endl;
    }
  }

  ft=AdaptivePenaltyFunction;
  if ((bdata[ft].mode==root)&&
      (bdata[ft].n_contributors>0)){
    if (!bdata[ft].aggregate.valid()){
      err++;
      if (cb_out()) 
	get_out()<<"**** WARNING: SumBundle::get_local_model_aggregate(...): aggregate for AdaptivePenaltyFunction part is not valid"<<std::endl;
    }
    else if (bdata[ft].get_model_aggregate(aggregate,factor,aft)){
      err++;
      if (cb_out()) 
	get_out()<<"**** WARNING: SumBundle::get_local_model_aggregate(...): collecting aggregate for AdaptivePenaltyFunction part failed"<<std::endl;
    }
  }

  return err;
}

// *****************************************************************************
//                             get_contributed_model_aggregate
// *****************************************************************************

int SumBundle::get_contributed_model_aggregate(MinorantPointer& aggregate,
					 Real factor,
					 const AffineFunctionTransformation* aft) const
{
  int err=0;
  FunctionTask ft=ObjectiveFunction;
  if ((bdata[ft].mode==child)&&
      (bdata[ft].n_contributors>0)){
    if (!bdata[ft].aggregate.valid()){
      err++;
      if (cb_out()) 
	get_out()<<"**** WARNING: SumBundle::get_contributed_model_aggregate(...): aggregate for ObjectiveFunction part is not valid"<<std::endl;
    }
    else if (bdata[ft].get_model_aggregate(aggregate,factor,aft)){
      err++;
      if (cb_out()) 
	get_out()<<"**** WARNING: SumBundle::get_contributed_model_aggregate(...): collecting aggregate for ObjectiveFunction part failed"<<std::endl;
    }
  }

  ft=ConstantPenaltyFunction;
  if ((bdata[ft].mode==child)&&
      (bdata[ft].n_contributors>0)){
    if (!bdata[ft].aggregate.valid()){
      err++;
      if (cb_out()) 
	get_out()<<"**** WARNING: SumBundle::get_contributed_model_aggregate(...): aggregate for ConstantPenaltyFunction part is not valid"<<std::endl;
    }
    else if (bdata[ft].get_model_aggregate(aggregate,factor,aft)){
      err++;
      if (cb_out()) 
	get_out()<<"**** WARNING: SumBundle::get_contributed_model_aggregate(...): collecting aggregate for ConstantPenaltyFunction part failed"<<std::endl;
    }
  }

  ft=AdaptivePenaltyFunction;
  if ((bdata[ft].mode==child)&&
      (bdata[ft].n_contributors>0)){
    if (!bdata[ft].aggregate.valid()){
      err++;
      if (cb_out()) 
	get_out()<<"**** WARNING: SumBundle::get_contributed_model_aggregate(...): aggregate for AdaptivePenaltyFunction part is not valid"<<std::endl;
    }
    else if (bdata[ft].get_model_aggregate(aggregate,factor,aft)){
      err++;
      if (cb_out()) 
	get_out()<<"**** WARNING: SumBundle::get_contributed_model_aggregate(...): collecting aggregate for AdaptivePenaltyFunction part failed"<<std::endl;
    }
  }

  return err;
}

// *****************************************************************************
//                             get_approximate_primal
// *****************************************************************************

const PrimalData* SumBundle::get_approximate_primal(FunctionTask ft) const
{
  if ((bdata[ft].mode==unavailable)||(bdata[ft].mode==inactive))
    return 0;

  if (!(bdata[ft].aggregate.valid())){
    if (cb_out()) 
      get_out()<<"**** WARNING: SumBundle::get_approximate_primal(...): corresponding aggregate is not valid for FunctionTask "<<ft<<std::endl;
    return 0;
  }
 
  return const_cast<MinorantPointer*>(&(bdata[ft].aggregate))->get_primal();
}

// *****************************************************************************
//                             call_primal_extender
// *****************************************************************************

int SumBundle::call_primal_extender(PrimalExtender& prex,Integer prex_id,FunctionTask ft)
{
  if (bdata[ft].mode==unavailable)
    return 0;
  
  return bdata[ft].call_primal_extender(prex,prex_id);
}
      
// *****************************************************************************
//                             apply_modification
// *****************************************************************************

int SumBundle::apply_modification(const GroundsetModification& gsmdf,Integer mod_id,MinorantExtender* mex,FunctionTask ft)
{
  if (bdata[ft].mode==unavailable)
    return 0;
  
  return bdata[ft].apply_modification(gsmdf,mod_id,mex);
}
      
    
  
 
}

