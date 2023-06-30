/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/BoxModelParameters.cxx
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



#include "BoxModelParameters.hxx"

using namespace CH_Matrix_Classes;

namespace ConicBundle {


BoxModelParametersObject::~BoxModelParametersObject()
{}
 
BoxModelParameters::~BoxModelParameters()
{}
  
//******************************************************************************
//                              select_model
//******************************************************************************

int BoxModelParameters::select_model(MinorantBundle& box_model, 
				     Matrix& box_coeff,  
				     Matrix& /* box_indicators */, 
				     Indexmatrix& box_coords, 
				     Matrix& box_complvalues, 
				     MinorantBundle& nnc_model, 
				     Matrix& nnc_coeff, 
				     Matrix& /* nnc_indicators */,
				     const Matrix& coord_switching,
				     const MinorantBundle& /* minorants */,
				     const MinorantPointer cand_minorant,
				     const PrimalMatrix cand_boxvec,
				     const PrimalMatrix aggr_boxvec,
				     Real aggr_scaleval,
				     BoxOracle* oracle,
				     Integer modification_id,
				     FunctionTask function_task,
				     Real function_factor,
				     BundleModel::ModelUpdate /* model_update */,
				     Integer /* center_id */,
				     const Matrix& /* center_y */,
				     Integer /* y_id */,
				     const Matrix& y,
				     Real /* model_maxviol */,
				     BundleProxObject& /* H */)
{
  assert(coord_switching.rowdim()==y.rowdim());
  Integer dim=y.rowdim();
  
  //-----  initialize to no coordinates selected, all in box_complvalues 
  if (aggr_boxvec.dim()==0){
    box_complvalues=cand_boxvec;
    aggr_scaleval=function_factor;
  }
  else {
    box_complvalues=aggr_boxvec;
  }
  box_coords.init(0,1,0.);
  box_coeff.init(0,1,0.);
  box_model.clear();
  nnc_model.clear();
  nnc_coeff.init(0,1,0.);

  //----- select the coordinates according to their switching rate
  Indexmatrix sind;
  sortindex(coord_switching,sind,false);
  Real val=1.;
  Integer i=0;
  Integer ind;
  Integer maxk=min(coord_switching.dim(),max_model_size>1?max_model_size-2:coord_switching.dim());
  //while((i<maxk)&&((coord_switching(ind=sind(i)))>100)){
  //while((i<maxk)&&((coord_switching(ind=sind(i)))>3.)){
  while((i<maxk)&&((coord_switching(ind=sind(i)))>1e-3)){
    if (oracle->get_upper_bounds()(ind)-oracle->get_lower_bounds()(ind)<1e-10*function_factor){
      maxk=min(maxk+1,coord_switching.dim());
      i++;
      continue;
    }
    box_coords.concat_below(ind);
    box_model.push_back(MinorantPointer(new Minorant(true,0.,1,&val,&ind,1.,0),modification_id));
    box_coeff.concat_below(box_complvalues(ind)*aggr_scaleval);
    box_complvalues(ind)=0.;
    i++;
  }

  //--- install the part not covered directly by coordinates
  if ((box_coords.dim()<dim)||(function_task!=ObjectiveFunction)||(i==0)){
    if (box_coords.dim()>0){
      // some coordinates selected
      box_model.push_back(MinorantPointer(new MatrixMinorant(0,box_complvalues,0,true),modification_id));
      box_coeff.concat_below(aggr_scaleval);
    }
    else {
      //no coordinates selected (box_coords.dim() == 0)
      if (aggr_boxvec.dim()!=0){
	nnc_model.push_back(MinorantPointer(new MatrixMinorant(0,aggr_boxvec,0,true),modification_id));
	nnc_coeff.concat_below(aggr_scaleval);
      }	
    }
    if ((box_coords.dim()<dim)||(box_coords.dim()==0)){
      nnc_model.push_back(cand_minorant);
      nnc_coeff.concat_below(nnc_model.size()==1?function_factor:0.);
    }
  }

  return 0;
}




} // end namespace ConicBundle
