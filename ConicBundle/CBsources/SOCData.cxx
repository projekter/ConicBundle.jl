/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/SOCData.cxx
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
#include "SOCData.hxx"
#include "SOCOracle.hxx"

using namespace CH_Matrix_Classes;

namespace ConicBundle {

// *****************************************************************************
//                              SOCData
// *****************************************************************************

  SOCData::SOCData(Real fun_factor,FunctionTask fun_task)
{
  assert(fun_factor>0.);
  set_function_factor()=fun_factor;
  set_function_task()=fun_task;
  SOCvecs_max=10;
  clear();
}

// *****************************************************************************
//                              ~SOCData
// *****************************************************************************

SOCData::~SOCData()
{  
  clear();
}

// *****************************************************************************
//                                clear
// *****************************************************************************

void SOCData::clear(Integer start_modification_id)
{
  BundleData::clear(start_modification_id);

  center_SOCval=0.;
  center_SOCvec.init(0,0,0.);
  center_minorant.clear();

  cand_SOCval=0.;
  cand_SOCvec.init(0,0,0.);
  cand_minorant.clear();

  bundlevecs.init(0,0,0.);
 
  SOCvecs.init(0,0,0.);
  SOCvecs_nexti=0;

}


// *****************************************************************************
//                                init
// *****************************************************************************

int SOCData::init(const BundleData* bd)
{
  const SOCData* mbd=dynamic_cast<const SOCData*>(bd);
  if (mbd==0){
    if (cb_out())
      get_out()<<"**** ERROR SOCData::init(): dynamic cast failed, argument is not of type const SOCData*"<<std::endl; 
    return 1;
  }
  clear();
  BundleData::init(bd);

  center_SOCval=mbd->center_SOCval;
  center_SOCvec=mbd->center_SOCvec;
  center_minorant=mbd->center_minorant;
  
  cand_SOCval=mbd->cand_SOCval;
  cand_SOCvec=mbd->cand_SOCvec;
  cand_minorant=mbd->cand_minorant;

  primalvec=mbd->primalvec;
  bundlevecs=mbd->bundlevecs;

  SOCvecs=mbd->SOCvecs;
  SOCvecs_nexti=mbd->SOCvecs_nexti;
  SOCvecs_max=mbd->SOCvecs_max;

  return 0;
}

// *****************************************************************************
//                                clone
// *****************************************************************************

BundleData* SOCData::clone() const
{
  SOCData* mbd=new SOCData;
  mbd->init(this);
  return mbd;
}

// *****************************************************************************
//                          store_SOCvec
// *****************************************************************************

int SOCData::store_SOCvec(const Matrix& socvec)
{
  if (SOCvecs_max<=0)
    return 0;
  Integer n=socvec.rowdim()-1;
  if (SOCvecs_max<SOCvecs.coldim()){
    Matrix tmpvecs(n,SOCvecs_max); chk_set_init(tmpvecs,1);
    SOCvecs_nexti-=(SOCvecs_max-1);
    if (SOCvecs_nexti<0)
      SOCvecs_nexti+=SOCvecs.coldim();
    for(Integer i=0;i<SOCvecs_max-1;i++){
      mat_xey(n,tmpvecs.get_store()+i*n,SOCvecs.get_store()+SOCvecs_nexti*n);
      SOCvecs_nexti++;
      if (SOCvecs_nexti>=SOCvecs.coldim())
	SOCvecs_nexti=0;
    }
    SOCvecs_nexti=SOCvecs_max-1;
    swap(tmpvecs,SOCvecs);
  }

  SOCvecs_nexti%=SOCvecs_max;

  Real d=std::sqrt(mat_ip(n,socvec.get_store()+1));
  if (d<1e-10){
    //store e_1
    if (SOCvecs_nexti==SOCvecs.coldim()){
      if (SOCvecs.coldim()==0)
	SOCvecs.init(n,1,0.);
      else 
	SOCvecs.enlarge_right(1,0.);
    }
    else {
      mat_xea(n,SOCvecs.get_store()+n*SOCvecs_nexti,0.);
    }
    SOCvecs(0,SOCvecs_nexti)=1.;
  }
  else {
    //store normalized part
    if (SOCvecs_nexti==SOCvecs.coldim()){
      if (SOCvecs.coldim()==0)
	SOCvecs.init(n,1,socvec.get_store()+1,1,1./d);
      else 
	SOCvecs.enlarge_right(1,socvec.get_store()+1,1./d);
    }
    else {
      assert((n==SOCvecs.rowdim())&&(0<=SOCvecs_nexti)&&(SOCvecs_nexti<SOCvecs.coldim()));
      mat_xeya(n,SOCvecs.get_store()+n*SOCvecs_nexti,socvec.get_store()+1,1./d);
    }
  }
  SOCvecs_nexti++;
  return 0;
}

// *****************************************************************************
//                             form_bundlevecs
// *****************************************************************************

int SOCData::form_bundlevecs(Integer max_columns)
{
  assert((primalvec.coldim()==1)||(bundlevecs.coldim()==0));
  assert(cand_SOCvec.rowdim()==SOCvecs.rowdim()+1);

  Matrix tmpvec;
  bundlevecs.init(0,0,0.);
  Real tol=1e-10;

  //-- insert the aggregate
  if (primalvec.coldim()==1) {
    tmpvec.init(primalvec.rowdim()-1,1,primalvec.get_store()+1);
    Real n2=norm2(tmpvec);
    if (n2<tol){
      tmpvec.init(primalvec.rowdim()-1,1,0.);
      tmpvec(0)=1.;
    }
    else
      tmpvec/=n2;
    bundlevecs.concat_right(tmpvec);
  }
  
  //-- insert the candidate
  if (cand_ub_mid==get_modification_id()){
    tmpvec.init(cand_SOCvec.rowdim()-1,1,cand_SOCvec.get_store()+1);
    if (bundlevecs.coldim()==1){
      //orthogonalize to the aggregate
      tmpvec.xpeya(bundlevecs,-ip(tmpvec,bundlevecs));
    } 
    Real n2=norm2(tmpvec);
    if (n2<tol){
      tmpvec.init(cand_SOCvec.rowdim()-1,1,0.);
      tmpvec(0)=1.;
    }
    else
      tmpvec/=n2;
    bundlevecs.concat_right(tmpvec);

    if (bundlevecs.coldim()==1){
      //generate a new primalvec
      if ((get_function_task()!=ObjectiveFunction)&&(cand_SOCval<0)){
	primalvec.init(cand_SOCvec.rowdim(),1,0.);
      }
      else {
	primalvec.init(cand_SOCvec,get_function_factor());
      }	
    }
  }

  assert(bundlevecs.coldim()>0);

  //add further columns if bundle_parameters suggest this
  if ((max_columns>bundlevecs.coldim())&&(SOCvecs.coldim()>1)){
    // number of vectors to add
    Integer k=max(max_columns-bundlevecs.coldim(),SOCvecs.coldim()-1);
    //orthogonalize the vectors to the existing ones
    Matrix tmpmat=SOCvecs;
    tmpmat-=bundlevecs*genmult(bundlevecs,SOCvecs,tmpvec,1.,0.,1);
    //compute basis to k largest singular values
    Symmatrix S;
    if (tmpmat.rowdim()<=tmpmat.coldim()){
      rankadd(tmpmat,S);
      S.eig(tmpmat,tmpvec,false);
      Integer i=0;
      while ((i<k+1)&&(tmpvec(i)>tol))
	i++;
      tmpmat.delete_cols(Range(i,tmpmat.coldim()-1),true);
      bundlevecs.concat_right(tmpmat);
    }
    else {
      rankadd(tmpmat,S,1.,0.,1);
      Matrix P;
      S.eig(P,tmpvec,false);
      Integer i=0;
      while ((i<k+1)&&(tmpvec(i)>tol))
	i++;
      P.delete_cols(Range(i,P.coldim()-1));
      genmult(tmpmat,P,tmpvec);
      bundlevecs.concat_right(tmpvec);
    }
    //orthogonalize once more
    bundlevecs.QR_factor();
    tmpmat.init(bundlevecs.rowdim(),min(bundlevecs.coldim(),max(2,max_columns)),0.);
    for (Integer i=0;i<tmpmat.coldim();i++)
      tmpmat(i,i)=1.;
    bundlevecs.Q_times(tmpmat,bundlevecs.coldim());
    swap(bundlevecs,tmpmat);
  }

  return 0;
} 

// *****************************************************************************
//                            synchronize_ids
// *****************************************************************************

int SOCData::synchronize_ids(Integer& new_center_ub_fid, 
				  Integer new_center_id,
				  Integer old_center_id,
				  Integer& new_cand_ub_fid, 
				  Integer new_cand_id,
				  Integer old_cand_id,
				  Integer& new_aggregate_id,
				  Integer new_prex_id)
{
  int err=BundleData::synchronize_ids(new_center_ub_fid,new_center_id,old_center_id,
				      new_cand_ub_fid,new_cand_id,old_cand_id,new_aggregate_id,new_prex_id);
  
  center_minorant.synchronize_ids(get_modification_id(),
			new_center_id,old_center_id,
			new_cand_id,old_cand_id,
			get_prex_id());
  cand_minorant.synchronize_ids(get_modification_id(),
			new_center_id,old_center_id,
			new_cand_id,old_cand_id,
			get_prex_id());
  
  return err;  
}

// *****************************************************************************
//                              do_step
// *****************************************************************************


int SOCData::do_step(Integer point_id)
{

  int retval=BundleData::do_step(point_id);
  if (retval){
    if (cb_out())
      get_out()<<"\n**** ERROR SOCData::do_step(.): BundleData::do_step(.) returned "<<std::endl;
    return retval;
  }

  center_SOCval=cand_SOCval;
  center_SOCvec=cand_SOCvec;
  center_minorant=cand_minorant;
  
  return retval;
}

// *****************************************************************************
//                          get_latest_minorants
// *****************************************************************************

int SOCData::get_latest_minorants(MinorantBundle& latest_minorants,
				  Integer max_number)
{
  SOCvecs_max=max(SOCvecs_max,max_number);
  return BundleData::get_latest_minorants(latest_minorants,max_number);
}

// *****************************************************************************
//                                   clear_model
// *****************************************************************************

void SOCData::clear_model(bool discard_minorants_only)
{

  BundleData::clear_model(discard_minorants_only);
  cand_minorant.clear();
  center_minorant.clear();

  if (!discard_minorants_only){
    clear_aggregates();
    
    center_SOCval=0.;
    center_SOCvec.init(0,0,0.);
  
    cand_SOCval=0.;
    cand_SOCvec.init(0,0,0.);

    primalvec.init(0,0,0.);

    bundlevecs.init(0,0,0.);
    SOCvecs.init(0,0,0.);
    SOCvecs_nexti=0;
  }
  
}

// *****************************************************************************
//                                   clear_aggregates
// *****************************************************************************

void SOCData::clear_aggregates()
{
  BundleData::clear_aggregates();
}

// *****************************************************************************
//                               call_primal_extender
// *****************************************************************************

int SOCData::call_primal_extender(PrimalExtender& prex,bool include_candidates)
{
  int err=BundleData::call_primal_extender(prex,include_candidates);

  if (local_aggregate.call_primal_extender(prex,get_prex_id())){
    if (cb_out()) 
      get_out()<<"**** WARNING: FunctionData::call_primal_extender(..): PrimalExtender::extend failed for local_aggregate"<<std::endl;
    err++;
  }

  center_minorant.clear(); //can be reconstructed any time 
  
  cand_minorant.clear(); //can be reconstructed any time

  SOCPrimalExtender* meprex=dynamic_cast<SOCPrimalExtender*>(&prex);
  if (meprex==0){
    if (cb_out()) 
      get_out()<<"**** WARNING: SOCData::call_primal_extender(..): cast to SOCPrimalExtender failed"<<std::endl;
    err++;
  }
  else {
    if ((get_center_modification_id()==get_modification_id())&&(meprex->extend_SOC(center_SOCvec))){
      if (cb_out()) 
	get_out()<<"**** WARNING: SOCData::call_primal_extender(..):  extending center_SOCvec failed"<<std::endl;
      err++;
      set_center_ub_mid()=-1;
    }

    if ((cand_ub_mid==get_modification_id())&&(meprex->extend_SOC(cand_SOCvec))){
      if (cb_out()) 
	get_out()<<"**** WARNING: SOCData::call_primal_extender(..):  extending cand_SOCvec failed"<<std::endl;
      err++;
      cand_SOCvec.init(0,0,0.);
      cand_ub_mid=-1;
    }

    if ((primalvec.coldim()==1)&&(meprex->extend_SOC(primalvec))){
      if (cb_out()) 
	get_out()<<"**** WARNING: SOCData::call_primal_extender(..):  extending primalvec failed"<<std::endl;
      err++;
      primalvec.init(0,0,0.);
    }

    if (SOCvecs.coldim()>0){
      Matrix tmpvec;
      Matrix tmpmat(0,0,0.);
      Integer n=SOCvecs.rowdim();
      for(Integer i=0;i<SOCvecs.coldim();i++){
	tmpvec.newsize(n+1,1); chk_set_init(tmpvec,1);
	tmpvec(0)=1.;
	mat_xey(n,tmpvec.get_store()+1,SOCvecs.get_store()+i*n);
	if (meprex->extend_SOC(tmpvec)){
	  if (i<SOCvecs_nexti)
	    SOCvecs_nexti--;
	}
	else {
	  assert((tmpvec.coldim()==1)&&(tmpvec.rowdim()>=1));
	  if (tmpvec.dim()>1){
	    Real d=std::sqrt(mat_ip(tmpvec.dim()-1,tmpvec.get_store()-1));
	    if (d<1e-10){
	      tmpvec.init(tmpvec.dim(),1,0.);
	      tmpvec(1)=1.;
	    }
	    else 
	      tmpvec/=d;
	  }
	  if (tmpmat.coldim()==0)
	    tmpmat.init(tmpvec.dim()-1,1,tmpvec.get_store()+1);
	  else {
	    assert(tmpmat.rowdim()==tmpvec.dim()-1);
	    tmpmat.enlarge_right(1,tmpvec.get_store()+1);
	  }
	}
      }
      swap(SOCvecs,tmpmat);
    }

    if ((bundlevecs.coldim()>0)&&((primalvec.coldim()==1)||(cand_SOCvec.coldim()==1))){
      form_bundlevecs(bundlevecs.coldim());
    }
    else {
      bundlevecs.init(0,0,0.);
    }
  }

  return err; 

}

// *****************************************************************************
//                               apply_modification
// *****************************************************************************

int SOCData::apply_modification(const GroundsetModification& gsmdf,MinorantExtender* mex)
{
  ///modification_id has already been increased
  int err=BundleData::apply_modification(gsmdf,mex);

  if ((center_minorant.valid())&&(center_minorant.apply_modification(gsmdf,get_modification_id(),mex))){
    if (cb_out()) 
      get_out()<<"**** WARNING: SOCData::apply_modification(..): apply_modification failed for center_minorant"<<std::endl;
    err++;
  }
  
  if ((cand_minorant.valid())&&(cand_minorant.apply_modification(gsmdf,get_modification_id(),mex))){
    if (cb_out()) 
      get_out()<<"**** WARNING: SOCData::apply_modification(..): apply_modification failed for center_minorant"<<std::endl;
    err++;
  }
  
  if (get_sumbundle().has_bundle_for(get_function_task())){
    if (set_sumbundle().apply_modification(gsmdf,get_modification_id(),mex,get_function_task())){
      if (cb_out()) 
	get_out()<<"**** WARNING: SOCData::apply_modification(..):  sumbundle.apply_modificaiton(..) failed"<<std::endl;
      err++;
    }
  }

  return err;
}


// *****************************************************************************
//                                get_approximate_primal
// *****************************************************************************

const PrimalData* SOCData::get_approximate_primal() const 
{
  return const_cast<MinorantPointer*>(&get_aggregate())->get_primal();
}

// *****************************************************************************
//                                get_center_primal
// *****************************************************************************

const PrimalData* SOCData::get_center_primal() const 
{
  return const_cast<MinorantPointer*>(&center_minorant)->get_primal();
}

// *****************************************************************************
//                                get_candidate_primal
// *****************************************************************************

const PrimalData* SOCData::get_candidate_primal() const 
{
  return const_cast<MinorantPointer*>(&cand_minorant)->get_primal();
}



}

