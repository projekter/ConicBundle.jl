/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/SumBundleParameters.hxx
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



#ifndef CONICBUNDLE_SUMBUNDLEPARAMETERS_HXX
#define CONICBUNDLE_SUMBUNDLEPARAMETERS_HXX


/**  @file SumBundleParameters.hxx
    @brief Header declaring the class ConicBundle::SumBundleParameters
    @version 1.0
    @date 2020-03-12
    @author Christoph Helmberg
*/


#include "SumBundleParametersObject.hxx"

namespace ConicBundle {

/** @ingroup InternalBundleModel

*/
//@{


/** @brief implements two versions of select_model() (update_rule in {0,1})

   See the general description of SumBundleParametersObject
   regarding design issues for select_model() routines.

*/
  class SumBundleParameters : public SumBundleParametersObject 
{
private:

  /** @brief default heuristic by curvatuer guess; selected by update_rule 0

      first it includes the aggregate and the minorant with highest value in
      the candidate (if there is one that beats the aggregate).

      After that it includes those minorants in sequence 
      that incur the largest curvature (a rough estiamte of it) 
      that the quadratic model would need in order to avoid
      a violation of more thant model_maxviol of this minorant.
      
  */
  int select_model0(CH_Matrix_Classes::Indexmatrix& model_indices,
		    CH_Matrix_Classes::Integer cand_id,
		    const CH_Matrix_Classes::Matrix& cand_y,
		    const MinorantBundle& minorants,
		    CH_Matrix_Classes::Integer aggr_index,
		    CH_Matrix_Classes::Real model_maxviol,
		    BundleProxObject& H,
		    BundleModel::ModelUpdate model_update);

  /** @brief simple heuristic by function value; selected by update_rule 1

      Include the aggragte.

      Sort the others by their value in candidate, and include
      them in this order (if not a duplicate) until the model_size 
      is reached.
  */
  int select_model1(CH_Matrix_Classes::Indexmatrix& model_indices,
		    CH_Matrix_Classes::Integer cand_id,
		    const CH_Matrix_Classes::Matrix& cand_y,
		    const MinorantBundle& minorants,
		    CH_Matrix_Classes::Integer aggr_index,
		    CH_Matrix_Classes::Real model_maxviol,
		    BundleProxObject& H,
		    BundleModel::ModelUpdate model_update);


public:
  /// default constructor with the possibility to set the output
  SumBundleParameters(const CBout* cb=0,int incr=-1):
    SumBundleParametersObject(cb,incr)
  {}

  /// constructor for customized parameters
  SumBundleParameters(int modelsize,int bundlesize=10,int updaterule=0,const CBout* cb=0,int incr=-1):
    SumBundleParametersObject(BundleParameters(modelsize,bundlesize,updaterule),cb,incr)
  {}

  /// copy constructor
  SumBundleParameters(const SumBundleParameters& sbp):
    CBout(sbp),SumBundleParametersObject(sbp)
  {}

    ///
  virtual ~SumBundleParameters()
  {}

  /// clone
  virtual BundleParameters* clone_BundleParameters() const
  { return new SumBundleParameters(*this); }

  
  /** @brief FunctionModel and SumBundleHandler call this for selecting the next minorants for a polyhedral model

    @param[out] model_indices 
        the indices of minorants selected for the model; index 0 is always the aggregate indicated by the input index aggr_index

    @param[in] cand_id
        the identifier of the candidate point supplied next 

    @param[in] cand_y
        the candidate (differnte from the center), close to it the model should be good

    @param[in] minorants
        the vector of MinorantPointer gives the minorants 
	out of which the model should be selected.

    @param[in] aggr_index
        the index of the aggregate within the minorants 

    @param[in] model_maxviol
        a minorant violated by this would have caused a null step

    @param[in] H
        the proximal term used for determining the given cand_y

    @param[in] model_update
        informs about whether cand_y is the result of a null_step or descent_step or aomw other a new set up.

    @return 
     - 0 on success
     - 1 on failure    
  */
  virtual int select_model(CH_Matrix_Classes::Indexmatrix& model_indices,
			   CH_Matrix_Classes::Integer cand_id,
			   const CH_Matrix_Classes::Matrix& cand_y,
			   const MinorantBundle& minorants,
			   CH_Matrix_Classes::Integer aggr_index,
			   CH_Matrix_Classes::Real model_maxviol,
			   BundleProxObject& H,
			   BundleModel::ModelUpdate model_update);
  
 
};

  //@}

}

#endif

