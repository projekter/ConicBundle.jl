/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/BoxModelParameters.hxx
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



#ifndef CONICBUNDLE_BOXMODELPARAMETERS_HXX
#define CONICBUNDLE_BOXMODELPARAMETERS_HXX

/**  @file BoxModelParameters.hxx
    @brief Header declaring the class ConicBundle::BoxModelParameters
    @version 1.0
    @date 2020-03-28
    @author Christoph Helmberg
*/

#include "BoxModelParametersObject.hxx"

namespace ConicBundle {
/** @ingroup InternalBundleModel

*/
//@{

/** @brief abstract interface for BoxModel for the model selection routine


*/

  class BoxModelParameters: public BoxModelParametersObject
{
protected:
  //int max_model_size;       // maximum number of minorants to be selected for the cutting model (numbers<=1 for no limit, numbers >=2 impose a strict limit)
  //int max_bundle_size;      // suggested maximum number of latest minorants stored for use in a model, for constructing variable metric information etc. (negative numbers give no preference; the size may be increased internally in case of confliciting requirements, eg. in n_model_size or by variable metric routines) 
  //int update_rule;      // in case several update rules are available

  //bool enforce_separate_model; // when forming combined models for a sum of functions, enforce a separate model for this function
  
public:

  /// default constructor
  BoxModelParameters(const CBout* cb=0,int cbinc=-1):
    BoxModelParametersObject(cb,cbinc)
  {}

  /// copy constructor for BundleParameters
  BoxModelParameters(const BundleParameters& bp,const CBout* cb=0,int cbinc=-1):
    BoxModelParametersObject(bp,cb,cbinc)
  {}

  ///copy constructor
  BoxModelParameters(const BoxModelParameters& sms):
    CBout(sms),BoxModelParametersObject(sms)
  {}

  ///virtual destructor, implemented in BoxModelParameters.cxx
  virtual ~BoxModelParameters();
  
  /// clone
  BundleParameters* clone_BundleParameters() const 
  { return new BoxModelParameters(*this); }

  /** @brief BoxModel calls this for selecting the next coordinates 
       for a specialized polyhedral model with a box part and an
       nnc part for aggregates, see the general explanation of the
       class.

       There is little experience on how to do this. 
       The current routine is minimalistic and simply uses
       the weighted number of switches in each coordinate
       between lower and upper bounds to select coordinates.

       To allow maybe better choices in future implementations
       the arguments try to pass all potentially relevant
       information items. 

       On putput the coefficient values of the new model must be
       feasible and have to generate the same aggregate (the aggregate
       is maintained), and the new candidate minorant must be in the
       feasible set.

    @param[in,out] box_model
       the boxmodel holds the minorants describing the BoxBlock part of
       the model for selected coordinates and, unless exact or empty,
       in the last position the complement coordinates of a feasible
       point (e.g. aggr_boxvec) 

    @param[in,out] box_coeff       
       - on input coefficients of BoxBlock determined in last 
         BundleMethod::eval_augmodel giving rise (together with nnc_coeff)
         to the aggrgate,
       - on output they match the new model and still give rise to
         the same aggregate
       
    @param[in,out] box_indicators 
       indicators for activity of box minorants, indicators may but need
       not be maintained

    @param[in,out] box_coords
       the coordinates selected to have their respective interval range 
       in the model 
			   
    @param[in,out] box_complvalues
       a point with 0 in the box_coords and feasible
       coordinate values in the complement (if not empty)

    @param[in,out] nnc_model 
       if box_model is empty, this spans at least the aggregate (if available) 
       and the candidate (always); if box_model is not empty but not 
       the entire box, nnc_model typically holds one of the candidate or 
       the aggregate; if box_model is the entire box, nnc_model is empty
			   
    @param[in,out] nnc_coeff  
       - on input coefficients of BoxBlock determined in last 
         BundleMethod::eval_augmodel giving rise (together with box_coeff)
         to the aggrgate,
       - on output they match the new model and still give rise to
         the same aggregate
			   
    @param[in,out] nnc_indicators 
       indicators for activity of minorants, indicators may but need
       not be maintained
			   
    @param[in] coord_switching
       keeps track of which coordinates where changing the most in the past
       by forming a weighted average in BoxModel::eval_function()
			   
    @param[in] minorants
        the vector of MinorantPointer gives additional minorants 
	collected over time (some may be duplicates also of those in
        nnc_model)

   @param[in] cand_minorant
       the (eps)sugradient linear minorant returned  
       by BoxModel::eval_function for the candidate (without function factor)
      
   @param[in] cand_boxvec
       the maximizer over the box for the current candidate

   @param[in] aggr_boxvec
       the primal aggregate vector in the box (without function_factor) 
       giving rise to the aggregate; not initialized if zerodimensional
			  
   @param[in] aggr_scaleval
       0<= aggr_scalevale <= function_factor, ==function_factor if 
       function_task==Objective_Function; the aggregate with 
       function_factor is aggr_boxvec*aggr_scaleval;

    @param[in] oracle
        gives access to lower and upper bounds of the box

    @param[in] modification_id
        the identifier of the current function version to be used in generating
        specialized minorants corresponding to the coordinate vectors

    @param[in] function_task
        see FunctionTask

    @param[in] function_factor
         interpreted according to function_task and the coefficients sum up to at most this value

    @param[in] model_update
        informs about whether cand_y is the result of a null_step or descent_step or a new set up.

    @param[in] center_id
        the identifier of the center point  

    @param[in] center_y
        the center point

    @param[in] cand_id
        the identifier of the candidate point 

    @param[in] cand_y
        the candidate (mostly differnt from the center), close to it the model should be good

    @param[in] model_maxviol
        a minorant violated by this would have caused a null step

    @param[in] H
        the variable metric used in the proximal term (function_factor is already removed in this) 

    @return 
     - 0 on success
     - 1 on failure    
  */
  virtual int select_model(MinorantBundle& box_model, 
			   CH_Matrix_Classes::Matrix& box_coeff,  
			   CH_Matrix_Classes::Matrix& box_indicators, 
			   CH_Matrix_Classes::Indexmatrix& box_coords, 
			   CH_Matrix_Classes::Matrix& box_complvalues, 
			   MinorantBundle& nnc_model, 
			   CH_Matrix_Classes::Matrix& nnc_coeff, 
			   CH_Matrix_Classes::Matrix& nnc_indicators,
			   const CH_Matrix_Classes::Matrix& coord_switching,
			   const MinorantBundle& minorants,
			   const MinorantPointer cand_minorant,
			   const PrimalMatrix cand_boxvec,
			   const PrimalMatrix aggr_boxvec,
			   CH_Matrix_Classes::Real aggr_scaleval,
			   BoxOracle* oracle,
			   CH_Matrix_Classes::Integer modification_id,
			   FunctionTask function_task,
			   CH_Matrix_Classes::Real function_factor,
			   BundleModel::ModelUpdate model_update,
			   CH_Matrix_Classes::Integer center_id,
			   const CH_Matrix_Classes::Matrix& center_y,
			   CH_Matrix_Classes::Integer cand_id,
			   const CH_Matrix_Classes::Matrix& cand_y,
			   CH_Matrix_Classes::Real model_maxviol,
			   BundleProxObject& H);
  
  

};


  //@}

}

#endif

