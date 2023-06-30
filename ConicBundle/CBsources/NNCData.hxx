/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/NNCData.hxx
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



#ifndef CONICBUNDLE_NNCDATA_HXX
#define CONICBUNDLE_NNCDATA_HXX


/**  @file NNCData.hxx
    @brief Header declaring the class ConicBundle::NNCData (needed for ConicBundle::NNCModel)
    @version 1.0
    @date 2014-07-24
    @author Christoph Helmberg
*/


#include "BundleData.hxx"

namespace ConicBundle {

/** @ingroup InternalBundleModel 

*/
//@{

/** @brief extends BundleData for use with FunctionModel by adding information on subgradients and the cutting plane model
 */

class NNCData:public BundleData
{
private:
  friend class NNCModel;

  //===== variables =====
  //CH_Matrix_Classes::Integer dim; ///< the incoming dimension, once it is detected, -1 if not yet clear (some need it, but it is mainly for consistency checks and should be initialized in eval_function, lb_function, recompute_center, and updated in apply_modification)
  //FunctionTask function_task;   
  //CH_Matrix_Classes::Real function_factor;  ///< >0, interpreted according to function_task

  //CH_Matrix_Classes::Integer modification_id; ///< increased on any modification applied to the function or @a function_factor

  //CH_Matrix_Classes::Integer center_ub_mid; ///< modification id for which center_ub was computed, <0 if not initialized
  //CH_Matrix_Classes::Real center_ub;     ///< upper bound on function value for modification id center_mid and center with point id center_id
  //CH_Matrix_Classes::Real center_relprec; ///< relative precision used in computing center_ub
  MinorantPointer center_minorant;   ///< the (eps-)subgradient linear minorant in the center (available if center_ub_mid is the modification id)


  //CH_Matrix_Classes::Integer cand_ub_mid; ///< modification id for which cand_ub was computed, <0 if not initialized
  //CH_Matrix_Classes::Real cand_ub;     ///< upper bound on function value for modification id cand_mid and candidate with point id cand_id
  //CH_Matrix_Classes::Real cand_relprec; ///< relative precision used in computing cand_ub
  //CH_Matrix_Classes::Integer cand_id;  ///< point id of the stored candidate
  MinorantPointer cand_minorant;    ///< the "best" (eps)sugradient linear minorant returned  by eval_function (available if valid and cand_ub_mid is the modification id)
  MinorantBundle cand_minorants;    ///< the (eps)sugradient linear minorants returned  by eval_function for the candidate (available if cand_ub_mid is the modification id)

  //CH_Matrix_Classes::Integer aggregate_id;  ///< if available and ids match, the aggregate need not be recomputed
  //MinorantPointer aggregate; ///< the aggregate linear minorant of the model (including @a function_factor and @a function_task)

  //SumBundle sumbundle;

  //CH_Matrix_Classes::Integer max_old_minorants; ///< ==0 by default, otherwise it gives the maximum number of the most recent candidate minorants to be stored 
  //MinorantBundle old_minorants;  ///< minorants collected for scaling and no other purposes, holds the at most max_old_minorants most recent candidate minorants

  //CH_Matrix_Classes::Matrix old_lowrank; ///< used by the dynamic scaling low rank heuristic (without function_factor)
  //CH_Matrix_Classes::Matrix old_diagonal; ///< used by the dynamic scaling low rank heuristic (without function_factor)
  //CH_Matrix_Classes::Symmatrix old_sym; ///< used by the dynamic scaling heuristic 
  //CH_Matrix_Classes::Integer minorant_nexti; ///< used by the dynamic scaling low rank heuristic 
 

  //the bundle is initialized if bundle.size()>0 (each subg might be of length 0!) 

  CH_Matrix_Classes::Matrix bundlecoeff;  ///< coefficients determined in last BundleMethod::eval_augmodel
  CH_Matrix_Classes::Matrix primal_tapia; ///< indicators for activity of minorants


  MinorantBundle bundle;  ///< the bundle holds the minorants describing the model

  bool outside_bundle; ///< true if bundle and bundlecoeff are initialized from an external source, e.g by a request to switch from sumbundle to local bundle; these may not yet satisfy the usual requirements specified by the bundle parameters


public:
  
  //===== methods =====
  /// calls clear()
  ~NNCData();
 
  /// reset to initial state (also used by the default constructor) 
  void clear(CH_Matrix_Classes::Integer start_modification_id=0);  

  /// initializes BundleData, sets center_primal to NULL and calls clear()
  NNCData(CH_Matrix_Classes::Real fun_factor=1.,FunctionTask fun_task=ObjectiveFunction);

  /// if @a bd is of type NNCData, initialize to this data  
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

  ///clear the cutting model and all function evaluations
  void clear_model(bool discard_minorants_only=false);   

  ///remove all aggregate minorants
  void clear_aggregates();          

  ///see the last argument of FunctionOracle::evaluate() 
  int call_primal_extender(PrimalExtender& prex,bool include_candidates=true);

  ///rearrange/extend the minorants according to the given groundset modifications 
  int apply_modification(const GroundsetModification&,MinorantExtender* mex);

  ///return the PrimalData corresponding to the aggregate
  const PrimalData* get_approximate_primal() const;

  ///return the PrimalData corresponding to the aggregate
  const PrimalData* get_center_primal() const;

  ///return the PrimalData delivered by the last call of FunctionOracle::evaluate()
  const PrimalData* get_candidate_primal() const;

  /// the minorants currently used in the model; the list may be empty or max contain other minorants than returned in get_latest_minorants(); the minorants still need to be mutliplied by function_factor
  virtual int get_model_data(MinorantBundle& model_minorants,
			     CH_Matrix_Classes::Matrix& model_coeff) const
  {
    model_minorants=bundle;
    model_coeff=bundlecoeff;
    return 0;
  }


};


  //@}

}

#endif

