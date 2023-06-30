/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/SOCData.hxx
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



#ifndef CONICBUNDLE_SOCDATA_HXX
#define CONICBUNDLE_SOCDATA_HXX


/**  @file SOCData.hxx
    @brief Header declaring the class ConicBundle::SOCData (needed for ConicBundle::FunctionModel)
    @version 1.0
    @date 2017-09-28
    @author Christoph Helmberg
*/


#include "BundleData.hxx"

namespace ConicBundle {

/** @ingroup InternalBundleModel 

*/
//@{

/** @brief extends BundleData for use with SOCModel (see there for the description of the function) by adding information on SOC_vectors generating the subgradients and the cutting model

    In difference to the original meaning of the variables in BundleData, 
    the function value and the subgradient information is not stored 
    explicitly, but instead the linear part and the SOCvalue/SOCvector 
    information are stored so that function value and subgradient can be 
    (re)computed upon need. In particular
    - center_ub, cand_ub   give the (hopefully) upper bound on the 
      maximum eigenvalue (SOCvalue) in center_y (with center_id), 
      in cand_y (with cand_id) 
    - subg_offset and subgrad are recomputed for the respective SOCvector
      v whenever subgrad.dim()==0 and the corresponding center or cand 
      information is up to date. In particular also for cand only the
      subgradient to the respective SOCvector will be computed, and for 
      center as well as for cand they will be computed to the values 
      subg_offset=function_factor* v^TCv and 
      subgrad_i=rhs_i-function_factor*(v^TA_iv) (i=1,...,ydim)
 */

class SOCData:public BundleData
{
private:
  friend class SOCModel;

  //===== variables =====
  //CH_Matrix_Classes::Integer dim; ///< the incoming dimension, once it is detected, -1 if not yet clear (some need it, but it is mainly for consistency checks and should be initialized in eval_function, lb_function, recompute_center, and updated in apply_modification)
  //FunctionTask function_task;   
  //CH_Matrix_Classes::Real function_factor;  ///< >0, interpreted according to function_task

  //CH_Matrix_Classes::Integer modification_id; ///< increased on any modification applied to the function

  //CH_Matrix_Classes::Integer center_ub_mid; ///< modification id for which center_ub was computed, <0 if not initialized
  //CH_Matrix_Classes::Real center_ub;     ///< upper bound on function value for modification id center_mid and center with point id center_id
  //CH_Matrix_Classes::Real center_relprec; ///< relative precision used in computing center_ub
  //CH_Matrix_Classes::Integer center_id;  ///< point id of the stored center
  CH_Matrix_Classes::Real center_SOCval; ///< lower bound on the center value
  CH_Matrix_Classes::Matrix center_SOCvec; ///< the SOCvector generating center_SOCval
  MinorantPointer center_minorant; ///< if not empty it matches center_SOCvec


  //CH_Matrix_Classes::Integer cand_ub_mid; ///< modification id for which cand_ub was computed, <0 if not initialized
  //CH_Matrix_Classes::Real cand_ub;     ///< upper bound on function value for modification id cand_mid and candidate with point id cand_id
  //CH_Matrix_Classes::Real cand_relprec; ///< relative precision used in computing cand_ub
  //CH_Matrix_Classes::Integer cand_id;  ///< point id of the stored candidate
  CH_Matrix_Classes::Real cand_SOCval; ///< lower bound on the candidate value
  CH_Matrix_Classes::Matrix cand_SOCvec; ///< the SOCvector generating cand_SOCval
  MinorantPointer cand_minorant; ///< if not empty, it matches cand_SOCvec

  //CH_Matrix_Classes::Integer aggregate_id;
  //MinorantPointer local_aggregate; ///< the local (possibly empty) aggregate linear minorant of the model handled by this and its descendants (including @a function_factor and @a function_task)
  //MinorantPointer aggregate; ///< the aggregate linear minorant of the model (including @a function_factor and @a function_task)
  CH_Matrix_Classes::Matrix primalvec;  ///< current aggregate soc vector (includes @a function_factor and @a function_task); if not valid it has coldim==0

  //the bundle is initialized if bundlevecs.coldim()>0 (then it must have rowdim>0!) 
  CH_Matrix_Classes::Matrix bundlevecs;   ///< the columns span the subspace of the SOC face (the bundle is initialized if bundlevecs.coldim()>0) 
  

  CH_Matrix_Classes::Matrix SOCvecs; ///< collects the barx parts of the old SOCvecs  
  CH_Matrix_Classes::Integer SOCvecs_max; ///< upper bound on stored SOCvecs
  CH_Matrix_Classes::Integer SOCvecs_nexti; ///< points at next "free" column

  //SumBundle sumbundle;


  //--- data needed for generic scaling
  //CH_Matrix_Classes::Integer max_old_minorants; 
  //MinorantBundle old_minorants; 
  //CH_Matrix_Classes::Matrix old_lowrank;
  //CH_Matrix_Classes::Matrix old_diagonal;
  //CH_Matrix_Classes::Symmatrix old_sym; 
  //CH_Matrix_Classes::Integer minorant_nexti; 

public:
  
  //===== methods =====
  /// calls clear()
  ~SOCData();
 
  /// reset to initial state (also used by the default constructor) 
  void clear(CH_Matrix_Classes::Integer start_modification_id=0);  

  /// calls clear()
  SOCData(CH_Matrix_Classes::Real fun_factor=1.,FunctionTask fun_task=ObjectiveFunction);

  /// if @a bd is of type SOCData, initialize to this data  
  int init(const BundleData* bd);

  ///return a pointer to a clone of this 
  BundleData* clone() const;
  
  ///if the candidate information is available and consitent for point_id, copy it from cand to center and return 0, otherwise return 1
  int do_step(CH_Matrix_Classes::Integer point_id);

  //int store_old_minorant(MinorantPointer& minorant);  //as in Bundledata

  ///if max_old_minorants > 0, it adds the SOCvec cyclically to SOCvecs keeping max_old_minorants of them
  virtual int store_SOCvec(const CH_Matrix_Classes::Matrix& SOCvec);

  /// starting with aggregate and cand_SOCvec add further ones as needed
  int form_bundlevecs(CH_Matrix_Classes::Integer max_columns);

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

  /// return the max_number latest minorants if available; 
  virtual int get_latest_minorants(MinorantBundle& latest_minorants,
				   CH_Matrix_Classes::Integer max_number);

};


  //@}

}

#endif

