/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/SumBundleParameters.cxx
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



#include "SumBundleParameters.hxx"

using namespace CH_Matrix_Classes;

namespace ConicBundle {

//******************************************************************************
//                              select_model
//******************************************************************************

int SumBundleParameters::select_model(Indexmatrix& model_indices,
				   Integer cand_id,
				   const Matrix& cand_y,
				   const MinorantBundle& minorants,
				   Integer aggr_index,
				   Real model_maxviol,
			           BundleProxObject& H,
				   BundleModel::ModelUpdate model_update)
{
  switch(update_rule){
  case 1: return select_model1(model_indices,cand_id,cand_y,minorants,aggr_index,model_maxviol,H,model_update);
  default: return select_model0(model_indices,cand_id,cand_y,minorants,aggr_index,model_maxviol,H,model_update);
  }
}
  
//******************************************************************************
//                              select_model0
//******************************************************************************

int SumBundleParameters::select_model0(Indexmatrix& model_indices,
					    Integer cand_id,
					    const Matrix& cand_y,
					    const MinorantBundle& minorants,
					    Integer aggr_index,
					    Real model_maxviol,
					    BundleProxObject& H,
					    BundleModel::ModelUpdate /* model_update */)
{
  assert((0<=aggr_index)&&(aggr_index<Integer(minorants.size())));

  if (max_model_size<2){
    max_model_size=2;
  }
  if (max_bundle_size<max_model_size){
    max_bundle_size=max(Integer(minorants.size()),max(max_model_size,50));
  }

  if (cb_out(1))
    get_out()<<" select0: mnrtsz="<<minorants.size()<<" msz="<<max_model_size<<" ["<<aggr_index<<"]";


  //first the aggregate is included if it exists
  model_indices.newsize(max_model_size,1);
  model_indices.init(1,1,aggr_index);

  //--- include the minorant with the largest value at y and initialize mnrtval
  Integer bsz=Integer(minorants.size());
  Matrix mnrtval(bsz,1); chk_set_init(mnrtval,1);
  Real maxval=min_Real;
  Integer besti=-1;
  for (Integer i=0;i<bsz;i++){
    Real v=minorants[unsigned(i)].evaluate(cand_id,cand_y);
    mnrtval(i)=v;
    if (maxval<v){
      maxval=v;
      besti=i;
    }
  }
  Real aggrval=mnrtval(aggr_index);
  assert(besti>=0);
  if ((besti!=aggr_index)
      &&(
	 (std::fabs(maxval-aggrval)>1e-10*(std::fabs(maxval)+1.))
	 ||(! minorants[unsigned(aggr_index)].equals(minorants[unsigned(besti)]))
	 )
      ){
    model_indices.concat_below(besti);
    if (cb_out(1))
      get_out()<<" ["<<besti<<"]";
  }

  if (model_indices.rowdim()>=max_model_size){
    if (cb_out(1)){
      get_out()<<" bsz="<<model_indices.rowdim();
      get_out().precision(12);
      get_out()<<" mviol="<<model_maxviol<<" weightu="<<H.get_weightu();
      get_out()<<" av="<<aggrval;
      get_out()<<std::endl;
    }
    return 0;
  }

  Indexmatrix skip(bsz,1,Integer(0));
  Matrix mnrtnormsqr(bsz,1,0.);

  skip(aggr_index)=1;
  Real aggrnorm=std::sqrt(mnrtnormsqr(aggr_index)=minorants[unsigned(aggr_index)].dual_norm_squared());
  

  if (model_indices.rowdim()>=2){
    skip(besti)=1;
    mnrtnormsqr(besti)=minorants[unsigned(besti)].dual_norm_squared();
  }

  //------------ determine the parameters for selecting more elements
  Real weightu=H.get_weightu();
  assert(model_maxviol>=0.);  
  assert(maxval>=aggrval-(1.+std::fabs(aggrval))*eps_Real);
  assert(weightu>=0.);
  Real violation_eps=max(.9*model_maxviol,(1+std::fabs(aggrval))*eps_Real);
  //violation_eps=max(violation_eps,max(maxval-aggrval,1e-10*(1.+std::fabs(maxval))));
  Real maxmult=1000*max(weightu,eps_Real);
  //Real normbnd=max(1e-6,0.01*mnrtnormsqr(aggr_index));
  Real normbnd=1e-4*max(1.,mnrtnormsqr(aggr_index));
  Real mindiff=.01*violation_eps;

  //------------ initialize bestnorm and find the largest with respect to the current bundle 
  Matrix violoffset(bsz,1,0.);
  Matrix minnormsqu(bsz,1,max_Real);
  Real maxnormsqu=min_Real;
  Real sumnorm=0;
  int too_bad=0;
  int too_small=0;
  besti=-1;
  for (Integer i=0;i<bsz;i++){
    if (skip(i))
      continue;

    Real mnsqri=mnrtnormsqr(i)=minorants[unsigned(i)].dual_norm_squared();    
    if (mnrtval(i)+.1*aggrnorm*std::sqrt(mnsqri)/weightu<aggrval-violation_eps){
      too_bad++;
      skip(i)=1;
      continue;
    }

    Real voff=violoffset(i)=violation_eps-mnrtval(i);

    for (Integer j=0;j<model_indices.rowdim();j++){
      Integer indj=model_indices(j);
      Real dv=voff+mnrtval(indj);
      if (dv<mindiff)
	//continue;
	dv=mindiff;
      Real ipval=minorants[unsigned(i)].ip(minorants[unsigned(indj)]);   
      Real diffnormsqu=mnrtnormsqr(indj)-2.*ipval+mnsqri;
      Real normsqu=.5*diffnormsqu/dv;
      if ((diffnormsqu<normbnd)||(normsqu<maxmult)){
	too_small++;
	minnormsqu(i)=normsqu;
	skip(i)=1;
	break;
      }
      if (normsqu<minnormsqu(i)){
	minnormsqu(i)=normsqu;
      }
    }
    assert(minnormsqu(i)<.1*max_Real);
    if ((minnormsqu(i)>maxnormsqu)&&(skip(i)==0)){
      maxnormsqu=minnormsqu(i);
      besti=i;
    }
  }

  if (maxnormsqu>maxmult){
    sumnorm+=std::sqrt(maxnormsqu);
    model_indices.concat_below(besti);
    if (cb_out(1))
      get_out()<<" "<<maxnormsqu<<"["<<besti<<"]";
  }
  else {
    besti=-1;
  }

  //------ check whether more bundle elements are allowed and worth it

  while((besti>=0)&&(model_indices.rowdim()<max_model_size)){

    skip(besti)=1;
    const MinorantPointer& bun(minorants[unsigned(besti)]);
    Real bunnormsqu=mnrtnormsqr(besti);
    Real bunval=mnrtval(besti);
    besti=-1;
    Real maxnormsqu=min_Real;
     
    for (Integer i=0;i<bsz;i++){
      if (skip(i))
	continue;      
      Real dv=violoffset(i)+bunval;
      if (dv>mindiff){
	Real ipval=minorants[unsigned(i)].ip(bun);   
	Real diffnormsqu=bunnormsqu-2.*ipval+mnrtnormsqr(i);
	Real normsqu=.5*diffnormsqu/dv;
	if ((diffnormsqu<normbnd)||(normsqu<maxmult)){
	  too_small++;
	  minnormsqu(i)=normsqu;
	  skip(i)=1;
	  continue;
	}
	if (normsqu<minnormsqu(i)){
	  minnormsqu(i)=normsqu;
	}
      }
      if (minnormsqu(i)>maxnormsqu){
	maxnormsqu=minnormsqu(i);
	besti=i;
      }
    }

    if (maxnormsqu>maxmult){
      sumnorm+=std::sqrt(maxnormsqu);
      model_indices.concat_below(besti);
      if (cb_out(1))
	get_out()<<" "<<maxnormsqu<<"["<<besti<<"]";
    }
    else {
      besti=-1;
    }
  } // end while((besti>=0)&& still room)

  if (cb_out(1)){
    get_out()<<" bsz="<<model_indices.rowdim();
    get_out().precision(12);
    get_out()<<" viol="<<violation_eps<<" mviol="<<model_maxviol<<" weightu="<<weightu;
    get_out()<<" av="<<aggrval<<" snorm="<<sumnorm;
    get_out()<<" bad="<<too_bad<<" small="<<too_small;
    get_out()<<std::endl;
  }

  return 0;        	       
}



//******************************************************************************
//                              select_model
//******************************************************************************

int SumBundleParameters::select_model1(Indexmatrix& model_indices,
				       Integer cand_id,
				       const Matrix& cand_y,
				       const MinorantBundle& minorants,
				       Integer aggr_index,
				       Real /* model_maxviol */,
				       BundleProxObject& /* H */,
				       BundleModel::ModelUpdate /* model_update */)
{
  assert((0<=aggr_index)&&(aggr_index<Integer(minorants.size())));
  
  if (max_model_size<2){
    max_model_size=2;
  }
  if (max_bundle_size<max_model_size){
    max_bundle_size=max(Integer(minorants.size()),max(max_model_size,50));
  }

  //collect the indices of the model minorants, the aggregate is always included
  model_indices.newsize(max_model_size,1);
  model_indices.init(1,1,aggr_index);
  
  //determine the minorants' values at y and sort according to this, remove duplicates
  Integer bsz=Integer(minorants.size());
  Matrix val(bsz,1); chk_set_init(val,1);
  Real maxval=min_Real;
  for (Integer i=0;i<bsz;i++){
    Real v=minorants[unsigned(i)].evaluate(cand_id,cand_y);
    val(i)=v;
    if (maxval<v)
      maxval=v;
  }
  Real aggrval=val(aggr_index);    
  Indexmatrix sind;
  sortindex(val,sind,false);
  Real maxdiff=max(1e-6,10.*(maxval-aggrval));
  
  for(Integer i=0;(i<bsz)&&(model_indices.dim()<max_model_size);i++){
    Integer ind=sind(i);
    if (ind==aggr_index)
      continue;
    //check whether its value is still good enough
    Real v=val(ind);
    if (maxval-v>=maxdiff)
      break;
    //check whether it is identical to a previous one
    bool included=false;
    for (Integer j=model_indices.dim();--j>=0;){
      if (CH_Matrix_Classes::abs(val(model_indices(j))-v)>1e-10*(1.+CH_Matrix_Classes::abs(v)))
	break;
      if (minorants[unsigned(ind)].equals(minorants[unsigned(model_indices(j))])){
	  included=true;
	  break;
      }
    }
    if (!included)
      model_indices.concat_below(ind);
  }

  if (cb_out(1)){
    get_out()<<" select1 bsz="<<model_indices.rowdim();
    get_out().precision(12);
    get_out()<<" av="<<aggrval;
    if (model_indices.dim()>0)
      get_out()<<"["<<val(model_indices(0))<<","<<val(model_indices(model_indices.dim()-1))<<"]";
    get_out()<<std::endl;
  }


  return 0;

}

} // end namespace ConicBundle
