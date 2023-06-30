/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/PSCModelParametersObject.hxx
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



#ifndef CONICBUNDLE_PSCMODELPARAMETERSOBJECT_HXX
#define CONICBUNDLE_PSCMODELPARAMETERSOBJECT_HXX

/**  @file PSCModelParametersObject.hxx
    @brief Header declaring the class ConicBundle::PSCModelParametersObject
    @version 1.0
    @date 2020-03-28
    @author Christoph Helmberg
*/

#include "SumBlockModel.hxx"
#include "PSCOracle.hxx"
//#include "PSCData.hxx"

namespace ConicBundle {
/** @ingroup InternalBundleModel

*/
//@{

/** @brief abstract interface for PSCModel for the model selection routine select_model()
    
  Because PSCModel implements a model for an AffineMatrixFunction
  given support function over the positive semidefinite cone, all
  subgradients are described by points inside this cone. For
  generating the model it suffices to know the extreme rays or
  faces. These may described by an orthonormal basis of the subspace
  of the eigenspaces of nonzero eigenvales. There is no need to
  consider the generated minorants at all as long as the generating
  primal matrices are covered by the subspace.

  The model is a face of the semidefinite cone spanned by the
  orthonormal columns of @a modelvecs -- this face includes the
  candidate and maybe part of the aggregate -- possibly in convex
  combination with a particular positive semidefinite matrix that is
  needed to span the @a aggregate; in slight abuse of notation this
  additional matrix is called the @a model_aggregate and its
  conribution to the full aggregate is the @a
  model_aggregate_coeff. The eigenvector-eigenvalue pair @a primalvecs
  and @a primaleigs (sorted nonincreasingly) describe the part of the
  aggregate in the face spanned by @a modelvecs on input; the columns
  of @a primalvecs are also orthonormal and span the same space as @a
  modelvecs on input. @a growthrate, @a primalgrowth and @a dualgrowth
  help to estimate how many of the eigenvalues could actually be
  active in the final optimal solution (the interior point code never
  gets to zero eigenvalues). For good convergence it seems important
  to keep at least the part that is reduced slower than the dual side
  and maybe more in addition to the new subspace due to the new
  eigenvalue computation.

  The results of the oracle call in the candidate point @a cand_y are
  stored in @a cand_Ritzvecs and @a cand_Ritzvals (sorted
  nonincreasingly).  Due to numerical imprecision or early termination
  they need not deliver the actual maximal eigenvectors and
  eigenvalues. In order to obtain a more stable guess of those, they
  are used together with the previously important subspace to first
  form the common subspace of both and to determine the spectral
  decomposition of the projection of the AffineMatrixFunction at
  cand_y. The results for the projected matrix are stored in @a
  topvecs (orthonormal columns) and @a Ritz_values (sorted
  nondecreasingly). Quite often the maximum Ritz_value is a bit larger
  than the value returned by the eigenvalue computation and the first
  column of @a topvecs is typically a better estimate tan that of 
  @a cand_Ritzvecs, but the latter brings in more new direction 
  infomation.
  
  The task of this routine is to 

  - form the new @a modelvecs basis from this information
  
  - decide how many of the topmost @a Ritz_values and @a topvecs 
    should be regarded as active (returned in @a activedim),

  - decide how many of the topmost @a primaleigs and @a primalvecs are
    actually included and return this number in @a modelvecs (all
    smaller ones will be aggregated in @a model_aggregate afterwards
    outside this routine),

  - decide how many additional elements of the topmost @a Ritz_values 
    and @a topvecs should be kept on top of @a activedim for the
    purpose of using this in variable metric heuristics. The 
    dimension of this additional subspace is returned in 
    @a skippedsize. The subspace of @a topvecs (and the
    corresponding elements of @a Ritz_values) beyond the
    first activedim+skippedsize elemnts should be deleted here,
    so that the matrix does not grow arbitrarily (usually it gets
    much bigger than hoped for anyways)
    
*/

  class PSCModelParametersObject: public virtual CBout, public BundleParameters
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
    const PSCModelParametersObject* mpo=dynamic_cast<const PSCModelParametersObject*>(&bp);
    if (mpo){
      set_cbout(mpo,0);
    }
    return 0;
  }

  /// default constructor
  PSCModelParametersObject(const CBout* cb=0,int cbinc=-1):
    CBout(cb,cbinc),BundleParameters()
  {}

  /// copy constructor for BundleParameters
  PSCModelParametersObject(const BundleParameters& bp,const CBout* cb=0,int cbinc=-1):
    CBout(cb,cbinc),BundleParameters(bp)
  {}

  ///copy constructor
  PSCModelParametersObject(const PSCModelParametersObject& sms):
    CBout(sms),BundleParameters(sms)
  {}

  ///virtual destructor, implemented in PSCModelParameters.cxx
  virtual ~PSCModelParametersObject();
  
  /// Compute the aggregate minornat corresponding to the matrix P*Diag(d)*P'/sum(d) for the given oracle and aggregate it to (or, if empty or mp_coeff=0., store it in) mp; P has orthonormal columns and d is a nonengative vector of appropriate size.
  int get_minorant(MinorantPointer& mp,
		   CH_Matrix_Classes::Real& mp_coeff,
		   const CH_Matrix_Classes::Matrix& P,
		   const CH_Matrix_Classes::Matrix& d,
		   PSCOracle* oracle,
		   CH_Matrix_Classes::Integer modification_id);
  
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
			   CH_Matrix_Classes::Real primal_Ritzval,
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
			   BundleProxObject& H)=0;
  

};


  //@}

}

#endif

