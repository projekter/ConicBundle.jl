/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/PSCData.hxx
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



#ifndef CONICBUNDLE_PSCDATA_HXX
#define CONICBUNDLE_PSCDATA_HXX


/**  @file PSCData.hxx
    @brief Header declaring the class ConicBundle::PSCData (needed for ConicBundle::FunctionModel)
    @version 1.0
    @date 2015-09-28
    @author Christoph Helmberg
*/


#include "BundleData.hxx"
#include "PSCOracle.hxx"
//#include "PSCAffineFunction.hxx"

namespace ConicBundle {

/** @ingroup InternalBundleModel 

*/
//@{

/** @brief extends BundleData for use with PSCModel (see there for the description of the function) by adding information on Ritz_vectors generating the subgradients and the cutting model

    In difference to the original meaning of the variables in BundleData, 
    the function value and the subgradient information is not stored 
    explicitly, but instead the linear part and the Ritzvalue/Ritzvector 
    information are stored so that function value and subgradient can be 
    (re)computed upon need. In particular
    - center_ub, cand_ub   give the (hopefully) upper bound on the 
      maximum eigenvalue (Ritzvalue) in center_y (with center_id), 
      in cand_y (with cand_id) 
    - subg_offset and subgrad are recomputed for the respective Ritzvector
      v whenever subgrad.dim()==0 and the corresponding center or cand 
      information is up to date. In particular also for cand only the
      subgradient to the respective Ritzvector will be computed, and for 
      center as well as for cand they will be computed to the values 
      subg_offset=function_factor* v^TCv and 
      subgrad_i=rhs_i-function_factor*(v^TA_iv) (i=1,...,ydim)
 */

class PSCData:public BundleData
{
private:
  friend class PSCModel;
  friend class QPPSCOracleData;
  
  //===== variables =====
  //CH_Matrix_Classes::Integer dim; ///< the incoming dimension, once it is detected, -1 if not yet clear (some need it, but it is mainly for consistency checks and should be initialized in eval_function, lb_function, recompute_center, and updated in apply_modification)
  //FunctionTask function_task;   
  //CH_Matrix_Classes::Real function_factor;  ///< >0, interpreted according to function_task

  //CH_Matrix_Classes::Integer modification_id; ///< increased on any modification applied to the function

  //CH_Matrix_Classes::Integer center_ub_mid; ///< modification id for which center_ub was computed, <0 if not initialized
  //CH_Matrix_Classes::Real center_ub;     ///< upper bound on function value for modification id center_mid and center with point id center_id
  //CH_Matrix_Classes::Real center_relprec; ///< relative precision used in computing center_ub
  //CH_Matrix_Classes::Integer center_id;  ///< point id of the stored center
  CH_Matrix_Classes::Real center_Ritzval; ///< lower bound on the center value
  CH_Matrix_Classes::Matrix center_Ritzvec; ///< a Ritzvector generating center_Ritzval
  MinorantPointer center_minorant; ///< if not empty it matches center_Ritzvec


  //CH_Matrix_Classes::Integer cand_ub_mid; ///< modification id for which cand_ub was computed, <0 if not initialized
  //CH_Matrix_Classes::Real cand_ub;     ///< upper bound on function value for modification id cand_mid and candidate with point id cand_id
  //CH_Matrix_Classes::Real cand_relprec; ///< relative precision used in computing cand_ub
  //CH_Matrix_Classes::Integer cand_id;  ///< point id of the stored candidate
  CH_Matrix_Classes::Matrix cand_Ritzval; ///< Ritz values returned by the oracle
  CH_Matrix_Classes::Matrix cand_Ritzvec; ///< the Ritzvectors generating cand_Ritzval
  MinorantPointer cand_minorant; ///< if not empty, it matches topvecs.col(0) (!)

  //CH_Matrix_Classes::Integer aggregate_id;
  //MinorantPointer local_aggregate; ///< the local (possibly empty) aggregate linear minorant of the model handled by this and its descendants (including @a function_factor and @a function_task)
  //MinorantPointer aggregate; ///< the aggregate linear minorant of the model (including @a function_factor and @a function_task)
  CH_Matrix_Classes::Matrix primaleigs;  ///< eigs of eigval-decomposition of last eval_augmodel (if not empty, sum of this and aggrcoeffs is feasible for function_task) 
  CH_Matrix_Classes::Matrix primalvecs;  ///< vecs of eigval-decomposition of last eval_augmodel (always the same span as bundlevecs)
  MinorantPointer primal_aggregate; ///< for getting inside the semidifinite cone
  CH_Matrix_Classes::Real primal_aggregate_coeff; ///< coefficient determined in last eval_augmodel, maybe with some aggregated parts of primalvecs if not all were kept
  CH_Matrix_Classes::Real primal_Ritzval; ///< the (common) Ritz_value of the active subspace in primaleigs (= cutting model value)

  CH_Matrix_Classes::Real growthrate; ///< factor <X,Z>/<X^-,Z^->, where X^- and Z⁻- are the last but one iterates of the interior point method   
  CH_Matrix_Classes::Matrix primalgrowth; ///< factor by which primaleigs changed in the last interior point iteration
  CH_Matrix_Classes::Matrix dualgrowth; ///< factor by which the dual Ritz values to primalvecs changed during the last interior point iteration 

  //the bundle is initialized if bundlevecs.coldim()>0 (then it must have rowdim>0!) 
  CH_Matrix_Classes::Matrix bundlevecs;   ///< the columns span the subspace of the semidefinte face (the bundle is initialized if bundlevecs.coldim()>0, then it must have rowdim>0!) 
  MinorantPointer bundle_aggregate; ///< for getting inside the semidifinite cone
  

  CH_Matrix_Classes::Matrix topvecs; ///< orthonormal columns giving the largest Ritz_values; the range of this includes the range of bundlevecs 
  CH_Matrix_Classes::Matrix Ritz_values; ///< of the vectors in topvecs and skippedvecs
  CH_Matrix_Classes::Integer activedim; ///< identified size of the active subspace
  CH_Matrix_Classes::Integer keepsize; ///< extended space dimension kept in bundlevecs
  CH_Matrix_Classes::Integer skippedsize; ///< number of columns in topvecs exceeding activedim 

  //SumBundle sumbundle;

  //CH_Matrix_Classes::Integer prex_id; // counter for keeping track of PriamlExtender calls

  //--- data needed for generic scaling
  //CH_Matrix_Classes::Integer max_old_minorants; 
  //MinorantBundle old_minorants; 
  //CH_Matrix_Classes::Matrix old_lowrank;
  //CH_Matrix_Classes::Matrix old_diagonal;
  //CH_Matrix_Classes::Symmatrix old_sym; 
  //CH_Matrix_Classes::Integer minorant_nexti;

  //--- data needed for collecting specific scaling relevant information
  CH_Matrix_Classes::Real scaling_gap;  ///< Ritz threshold for consideration in scaling
  CH_Matrix_Classes::Real center_gap;  ///< the gap observed for the current center

public:
  
  //===== methods =====
  /// calls clear()
  ~PSCData();
 
  /// reset to initial state (also used by the default constructor) 
  void clear(CH_Matrix_Classes::Integer start_modification_id=0);  

  /// calls clear()
  PSCData(CH_Matrix_Classes::Real fun_factor=1.,FunctionTask fun_task=ObjectiveFunction);

  /// returns primalvecs
  const CH_Matrix_Classes::Matrix& get_primalvecs() const
  { return primalvecs; }
    
  /// returns primaleigs
  const CH_Matrix_Classes::Matrix& get_primaleigs() const
  { return primaleigs; }

  /// returns topvecs
  const CH_Matrix_Classes::Matrix& get_topvecs() const
  { return topvecs; }
    
  /// returns Ritz_values
  const CH_Matrix_Classes::Matrix& get_Ritz_values() const
  { return Ritz_values; }
    


  /// returns keepsize
  CH_Matrix_Classes::Integer get_keepsize() const
  { return keepsize; }
    
  /// returns acitvedim
  CH_Matrix_Classes::Integer get_activedim() const
  { return activedim; }
    
  /// returns skippedsize
  CH_Matrix_Classes::Integer get_skippedsize() const
  { return skippedsize; }
    
  
  /// if @a bd is of type PSCData, initialize to this data  
  int init(const BundleData* bd);

  ///return a pointer to a clone of this 
  BundleData* clone() const;
  
  ///if the candidate information is available and consitent for point_id, copy it from cand to center and return 0, otherwise return 1
  int do_step(CH_Matrix_Classes::Integer point_id);

  ///synchronize ids in any case by discarding inconsistent parts but return number of errors
  int synchronize_ids(CH_Matrix_Classes::Integer& new_center_ub_fid, 
			      CH_Matrix_Classes::Integer new_center_id,
			      CH_Matrix_Classes::Integer old_center_id,
			      CH_Matrix_Classes::Integer& new_cand_ub_fid, 
			      CH_Matrix_Classes::Integer new_cand_id,
			      CH_Matrix_Classes::Integer old_cand_id,
			      CH_Matrix_Classes::Integer& new_aggregate_id,
			      CH_Matrix_Classes::Integer new_prex_id=0);

  //int center_modified()            // as in BundleData

  //int model_aggregate_modified()   // as in BundleData

  ///clear all function evaluations and aggregates, but not the collected subspace information
  void clear_model_except_bundlevecs();   

  ///clear the cutting model and all function evaluations
  void clear_model(bool discard_minorants_only=false);   

  ///delete all kinds of aggregates but keep explicit parts of the cutting model
  void clear_aggregates();

  ///see the last argument of FunctionOracle::evaluate() 
  int call_primal_extender(PrimalExtender&,bool include_candidates=true);

  ///rearrange/extend the minorants according to the given groundset modifications 
  int apply_modification(const GroundsetModification&,MinorantExtender* mex);

  ///return the PrimalData corresponding to the aggregate
  const PrimalData* get_approximate_primal() const;

  ///return the PrimalData corresponding to the aggregate
  const PrimalData* get_center_primal() const;

  ///return the PrimalData delivered by the last call of FunctionOracle::evaluate()
  const PrimalData* get_candidate_primal() const;

};


  //@}

}

#endif

