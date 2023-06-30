/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/SOCModelParameters.hxx
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



#ifndef CONICBUNDLE_SOCMODELPARAMETERS_HXX
#define CONICBUNDLE_SOCMODELPARAMETERS_HXX


/**  @file SOCModelParameters.hxx
    @brief Header declaring the class ConicBundle::SOCModelParameters
    @version 1.0
    @date 2020-03-28
    @author Christoph Helmberg
*/


#include "SOCModelParametersObject.hxx"

namespace ConicBundle {

/** @ingroup InternalBundleModel

*/
//@{


/** @brief default model selection routine for SOCModel

  There is not yet much experience with this and it
  is certainly worth to try to improve this.

*/
  class SOCModelParameters : public SOCModelParametersObject
{
private:

public:
  /// default constructor with the possibility to set the output
  SOCModelParameters(const CBout* cb=0,int incr=-1):
    SOCModelParametersObject(cb,incr)
  {}

  /// constructor for size parameters
  SOCModelParameters(int modelsize,int bundlesize=10,int updaterule=0,const CBout* cb=0,int incr=-1):
    SOCModelParametersObject(BundleParameters(modelsize,bundlesize,updaterule),cb,incr)
  {}

  /// copy constructor for BundleParameters
  SOCModelParameters(const BundleParameters& bp,const CBout* cb=0,int incr=-1):
    SOCModelParametersObject(bp,cb,incr)
  {}

  ///copy constructor
  SOCModelParameters(const SOCModelParameters& sms):
    CBout(sms),SOCModelParametersObject(sms)
  {}

  /// destructor
  virtual ~SOCModelParameters();

  /// clone
  BundleParameters* clone_BundleParameters() const 
  { return new SOCModelParameters(*this); }
  
  
    /** @brief SOCModel calls this for selecting the next minorants for a polyhedral model

    @param[in,out] modelvecs
        the columns span the barx subspace of the SOC face (the model is initialized if modelvecs.coldim()>0). On output it has to span at least the subspace spanned by aggrvec and cand_SOCvec

    @param[in,out] aggrvec  
        current aggregate soc vector (includes all coordinates and @a function_factor according to @a function_task); if not valid it has coldim==0.

    @param[in] cand_SOCval
        lower bound on the candidate value

    @param[in] cand_SOCvec
        the SOCvector generating cand_SOCval (includes all coordinates but no function_factor)

    @param[in] center_SOCval
        lower bound on the center value

    @param[in] center_SOCvec
        the SOCvector generating center_SOCval (includes all coordinates but no function_factor)

    @param[in] SOCvecs
        collects the barx parts of the old SOCvecs  

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
  virtual int select_model(CH_Matrix_Classes::Matrix& modelvecs,
			   const CH_Matrix_Classes::Matrix& aggrvec,
			   CH_Matrix_Classes::Real cand_SOCval,
			   const CH_Matrix_Classes::Matrix& cand_SOCvec,
			   CH_Matrix_Classes::Real center_SOCval,
			   const CH_Matrix_Classes::Matrix& center_SOCvec,
			   const CH_Matrix_Classes::Matrix& SOCvecs,
			   SOCOracle* oracle,
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

