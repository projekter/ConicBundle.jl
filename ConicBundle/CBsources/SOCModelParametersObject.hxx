/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/SOCModelParametersObject.hxx
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



#ifndef CONICBUNDLE_SOCMODELPARAMETERSOBJECT_HXX
#define CONICBUNDLE_SOCMODELPARAMETERSOBJECT_HXX

/**  @file SOCModelParametersObject.hxx
    @brief Header declaring the class ConicBundle::SOCModelParametersObject
    @version 1.0
    @date 2020-03-28
    @author Christoph Helmberg
*/

#include "SumBlockModel.hxx"
#include "SOCOracle.hxx"
//#include "SOCData.hxx"

namespace ConicBundle {
/** @ingroup InternalBundleModel

*/
//@{

/** @brief abstract interface for SOCModel for the model selection routine select_model()
    
  Coordinate zero of a vector x in a second order cone of dimension dim
  gives the maxium Euclidean norm of the vector of the coordinates 1 to dim-1
  (the barx part).

  Because SOCModel implements a model for an AffineMatrixTransformation
  given support function over the Second order cone, all subgradients 
  are described by points inside this cone. For generating the model
  it suffices to know the extreme rays in the barx space or some 
  appropriate basis of the barx-space containing them. There
  is no need to consider the generated minorants at all as long as
  the generating primal points are covered by the subspace.

  The model is always a scond order cone of size @a
  modelvecs.coldim()+1 and the columns of @a modelvecs form an
  orthonormal basis of some subspace of the barx space (coordinates 1
  to dim-1). This subspace must always span the aggregate stored in
  @a aggrvec and the candidate vector @a cand_SOCvec, so typically 
  it is a least two dimensional but may be one-dimensional in special 
  cases (e.g. if there is only one minorant to start with). Older
  candidate vectors may be available in @a SOCvecs in order to
  select new @a modelvecs.
    
*/

  class SOCModelParametersObject: public virtual CBout, public BundleParameters
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
    const SOCModelParametersObject* mpo=dynamic_cast<const SOCModelParametersObject*>(&bp);
    if (mpo){
      set_cbout(mpo,0);
    }
    return 0;
  }

  /// default constructor
  SOCModelParametersObject(const CBout* cb=0,int cbinc=-1):
    CBout(cb,cbinc),BundleParameters()
  {}

  /// copy constructor for BundleParameters
  SOCModelParametersObject(const BundleParameters& bp,const CBout* cb=0,int cbinc=-1):
    CBout(cb,cbinc),BundleParameters(bp)
  {}

  ///copy constructor
  SOCModelParametersObject(const SOCModelParametersObject& sms):
    CBout(sms),BundleParameters(sms)
  {}

  ///virtual destructor, implemented in SOCModelParameters.cxx
  virtual ~SOCModelParametersObject();
  
    /** @brief SOCModel calls this for selecting the next minorants for a polyhedral model

    @param[in,out] modelvecs
        the columns span the barx subspace of the SOC face (the model is initialized if modelvecs.coldim()>0). On output it has to span at least the subspace spanned by aggrvec and cand_SOCvec

    @param[in] aggrvec  
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
			   BundleProxObject& H)=0;

};


  //@}

}

#endif

