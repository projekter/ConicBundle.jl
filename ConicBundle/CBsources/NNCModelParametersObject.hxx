/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/NNCModelParametersObject.hxx
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



#ifndef CONICBUNDLE_NNCMODELPARAMETERSOBJECT_HXX
#define CONICBUNDLE_NNCMODELPARAMETERSOBJECT_HXX

/**  @file NNCModelParametersObject.hxx
    @brief Header declaring the class ConicBundle::NNCModelParametersObject
    @version 1.0
    @date 2020-03-28
    @author Christoph Helmberg
*/

#include "MatrixCBSolver.hxx"
#include "SumBlockModel.hxx"

namespace ConicBundle {
/** @ingroup InternalBundleModel

*/
//@{

/** @brief abstract interface for NNCModel for the model selection routine select_model()


*/

  class NNCModelParametersObject: public virtual CBout, public BundleParameters
{
protected:
  //int max_model_size;       // maximum number of minorants to be selected for the cutting model (numbers<=1 for no limit, numbers >=2 impose a strict limit)
  //int max_bundle_size;      // suggested maximum number of latest minorants stored for use in a model, for constructing variable metric information etc. (negative numbers give no preference; the size may be increased internally in case of confliciting requirements, eg. in n_model_size or by variable metric routines) 
  //int update_rule;      // in case several update rules are available

  //bool enforce_separate_model; // when forming combined models for a sum of functions, enforce a separate model for this function
  
public:
  /// initialize BundleParameters to given values, if bp is a SumModleParametersObject alos set n_local_models to its values, otherwise leave it unchanged
  int init(const BundleParameters& bp){
    BundleParameters::init(bp);
    const NNCModelParametersObject* mpo=dynamic_cast<const NNCModelParametersObject*>(&bp);
    if (mpo){
      set_cbout(mpo,0);
    }
    return 0;
  }

  /// default constructor
  NNCModelParametersObject(const CBout* cb=0,int cbinc=-1):
    CBout(cb,cbinc),BundleParameters()
  {}

  /// copy constructor for BundleParameters
  NNCModelParametersObject(const BundleParameters& bp,const CBout* cb=0,int cbinc=-1):
    CBout(cb,cbinc),BundleParameters(bp)
  {}

  ///copy constructor
  NNCModelParametersObject(const NNCModelParametersObject& sms):
    CBout(sms),BundleParameters(sms)
  {}

  ///virtual destructor, implemented in NNCModelParameters.cxx
  virtual ~NNCModelParametersObject();
  
  /** @brief NNCModel calls this for selecting the next minorants for a polyhedral model

    @param[in,out] model
        contains the model of the previous subproblem or a minimal starting choice

    @param[in,out]  coefficients
        contains the coefficients resulting from the last bundel subproblem (on input and on output they should generate the aggregate, if the aggregate is valid) 

    @param[in,out]  activity_indicators
        contains activity_indicators resulting from the last bundel subproblem (currently 1 if considered active and 0 otherwise). Those indicated active in the bundle on input will also be indicated active 

    @param[in] aggregate
        if valid, it holds the current aggregate. It arises as nonnegative combinations of the model by the coefficients (so it inculde the current function_factor). This property of model and coefficients should be maintained on output

    @param[in] center_minorant
        if valid, it holds the current aggregate. It arises as nonnegative combinations of the model by the coefficients, and that should be maintained on output

    @param[in] cand_minorants
        holds the vector of minorants returned by evaluation in the candidate point

    @param[in] old_minorants
        the vector of MinorantPointer gives additional minorants 
	collected over time (some may be duplicates, some are
	most certainly already contained in the bundle on input)

    @param[in] oracle
        gives access to the evaluation oracle

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
  virtual int select_model(MinorantBundle& model,
			   CH_Matrix_Classes::Matrix& coefficients,
			   CH_Matrix_Classes::Matrix& activity_indicators,
			   const MinorantPointer& aggregate,
			   const MinorantPointer& center_minorant,
			   const MinorantBundle& cand_minorants,
			   const MinorantBundle& old_minorants,
			   MatrixFunctionOracle* oracle,
			   FunctionTask function_task,
			   CH_Matrix_Classes::Real function_factor,
			   BundleModel::ModelUpdate model_update,
			   CH_Matrix_Classes::Integer center_id,
			   const CH_Matrix_Classes::Matrix& center_y,
			   CH_Matrix_Classes::Integer cand_id,
			   const CH_Matrix_Classes::Matrix& cand_y,
			   CH_Matrix_Classes::Real model_maxviol,
			   BundleProxObject& H)=0;
  
  

};


  //@}

}

#endif

