/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/PSCModelParameters.hxx
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



#ifndef CONICBUNDLE_PSCMODELPARAMETERS_HXX
#define CONICBUNDLE_PSCMODELPARAMETERS_HXX


/**  @file PSCModelParameters.hxx
    @brief Header declaring the class ConicBundle::PSCModelParameters
    @version 1.0
    @date 2020-03-28
    @author Christoph Helmberg
*/


#include "PSCModelParametersObject.hxx"

namespace ConicBundle {

/** @ingroup InternalBundleModel

*/
//@{


/** @brief default model selection routine for PSCModel

  The current routine is quite advanced and stable for a reasonalby
  wide range of applications. It clearly favors stability and
  precision over fast initial convergence. It ignores all parameter
  settings and picks its own sizes.

*/
  class PSCModelParameters : public PSCModelParametersObject
{
private:
  CH_Matrix_Classes::Real cutoffval; ///< an averaged bound on Ritzvalues for being included in the model
  CH_Matrix_Classes::Real gapsz; ///< a weighted average of the gaps between maximum Ritz_value and the Ritz_value of the largest active index
  CH_Matrix_Classes::Real Ritz_skipval; ///< a weighted average for Ritz_values to be discarded
  CH_Matrix_Classes::Integer n_nullsteps; ///< counts the number of null steps since last descent step 

public:
  /// default constructor with the possibility to set the output
  PSCModelParameters(const CBout* cb=0,int incr=-1):
    PSCModelParametersObject(cb,incr)
  {cutoffval=CH_Matrix_Classes::max_Real;gapsz=0.; Ritz_skipval=CH_Matrix_Classes::max_Real;n_nullsteps=0;}

  /// constructor for size parameters
  PSCModelParameters(int modelsize,int bundlesize=10,int updaterule=0,const CBout* cb=0,int incr=-1):
    PSCModelParametersObject(BundleParameters(modelsize,bundlesize,updaterule),cb,incr)
  {cutoffval=CH_Matrix_Classes::min_Real;gapsz=0.; Ritz_skipval=CH_Matrix_Classes::max_Real;n_nullsteps=0;}

  /// copy constructor for BundleParameters
  PSCModelParameters(const BundleParameters& bp,const CBout* cb=0,int incr=-1):
    PSCModelParametersObject(bp,cb,incr)
  {cutoffval=CH_Matrix_Classes::max_Real;gapsz=0.; Ritz_skipval=CH_Matrix_Classes::max_Real;n_nullsteps=0;}

  ///copy constructor
  PSCModelParameters(const PSCModelParameters& sms):
    CBout(sms),PSCModelParametersObject(sms)
  {cutoffval=sms.cutoffval;gapsz=sms.gapsz; Ritz_skipval=sms.Ritz_skipval;n_nullsteps=sms.n_nullsteps;}

  /// destructor
  virtual ~PSCModelParameters();

  /// clone
  BundleParameters* clone_BundleParameters() const 
  { return new PSCModelParameters(*this); }
  
  /** @brief PSCModel calls this for selecting the next positive semidefinite model

    @param[in,out] modelvecs
        * on input: orthonormal basis of the subspace used in the last 
          semidefinite model
        * on output: orthonormal subspace basis to be used in the next 
          semidefinite model, it includes the subspace of the primalvecs of 
          the keepsize largest primaleigs. In particular for null steps
	  the remaining columns of primalvecs with nonzero primaleigs
          have to be included with the current aggregate in a new aggregate
          of the model; this will be taken care of outside this routine 
          aftwards
    
    @param[in,out] model_aggregate (MinorantPointer)
        aggregate in use in the last and then the next model. 

    @param[in,out] topvecs
        * on input: orthonormal basis of the collected subspace that is
          supposed to approximate the eigenspace to the largest eigenvalues,
	  see Ritz_values
        * on output: same thing but maybe reduced in size
    
    @param[in,out] Ritz_values
        * on input: Ritz_values in cand_y for the vectors in topvecs
        * on output: same thing but maybe reduced in size as in topvecs
    
    @param[in,out] activedim
        * on input: the dimension of the subspace (first columns in topvecs) 
           regarded as active in the last iterations of the bundle subproblem solution 
	* on output: the dimension of the subspace (first columns of topvecs) 
	  regarded as (maybe weakly) active now
           (typically it will increase during null steps and may be tightened
           at descent steps) 

    @param[out] keepsize
        columns 0..keepsize-1 of primalvecs (corresponding to the keepsize 
	largest primaleigs) are included in modelvecs. The remaining
        columns need to be aggregated afterwards into the aggregate of the model

    @param[out] skippedsize
        columns activedim..activedim+skippedsize-1 of topvecs should be used 
        for setting up the scaling matrix after descent steps.

    @param[in] primal_Ritzval
        the (common) Ritz value of the active subspace of the model 
        (if not available use some guess like Ritz_values(0))

    @param[in] primaleigs
        eigenvalues of the last primal semdifinite model matrix
        (sorted nonincreasingly)

    @param[in] primalvecs
        corresponding (orthonormal) eigenvectors to primaleigs

    @param[in] primal_aggregate (MinorantPointer)
        aggregate in use in the last model. 

    @param[in] primal_aggregate_coeff
        coefficient on how strongly the aggregated was used in the last 
        primal solution to the model

    @param[in] growthrate (Real)
        factor <X,Z>/<X^-,Z^->, where X^- and Z⁻- are the last but one 
        iterates of the interior point method
 
    @param[in] primalgrowth (Matrix)
        factor by which primaleigs changed in the last interior point iteration 

    @param[in] dualgrowth (Matrix)
        factor by which the dual Ritz values to primalvecs changed 
	during the last interior point iteration 

    @param[in] cand_Ritzvec (const Matrix&)
        the (orthonormal) vectors returned by the evaluation call to the oracle

    @param[in] cand_Ritzval (const Matrix&)
	the Ritz values of cand_Ritzvec returned by the evaluation call to the oracle

    @param[in] oracle
        gives access to the evaluation oracle

    @param[in] modification_id
       the identifier for the current version of the function accounting for dynamic modifications

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

    @param[in] diffval_center_aggregate
        difference of center value to aggregate value (nonnegative, without function_factor) 

    @param[in] H
        the variable metric used in the proximal term (function_factor is already removed in this) 

    @return 
     - 0 on success
     - 1 on failure    
  */
  virtual int select_model(CH_Matrix_Classes::Matrix& modelvecs,
			   MinorantPointer& model_aggregate,
			   CH_Matrix_Classes::Matrix& topvecs,
			   CH_Matrix_Classes::Matrix& Ritz_values,
			   CH_Matrix_Classes::Integer& activedim,
			   CH_Matrix_Classes::Integer& keepsize,
			   CH_Matrix_Classes::Integer& skippedsize,
			   const CH_Matrix_Classes::Real primal_Ritzval,
			   const CH_Matrix_Classes::Matrix& primaleigs,
			   const CH_Matrix_Classes::Matrix& primalvecs,
			   const MinorantPointer& primal_aggregate,
			   CH_Matrix_Classes::Real primal_aggregate_coeff,
			   CH_Matrix_Classes::Real growthrate,
			   const CH_Matrix_Classes::Matrix& primalgrowth,
			   const CH_Matrix_Classes::Matrix& dualgrowth,
			   const CH_Matrix_Classes::Matrix& cand_Ritzvec,
			   const CH_Matrix_Classes::Matrix& cand_Ritzval,
			   PSCOracle* oracle,
			   CH_Matrix_Classes::Integer modification_id,
			   FunctionTask function_task,
			   CH_Matrix_Classes::Real function_factor,
			   BundleModel::ModelUpdate model_update,
			   CH_Matrix_Classes::Integer center_id,
			   const CH_Matrix_Classes::Matrix& center_y,
			   CH_Matrix_Classes::Integer cand_id,
			   const CH_Matrix_Classes::Matrix& cand_y,
			   CH_Matrix_Classes::Real model_maxviol,
			   CH_Matrix_Classes::Real diffval_center_aggregate,
			   BundleProxObject& H);
  

 
};

  //@}

}

#endif

