/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/BundleData.hxx
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



#ifndef CONICBUNDLE_BUNDLEDATA_HXX
#define CONICBUNDLE_BUNDLEDATA_HXX


/**  @file BundleData.hxx
    @brief Header declaring the class ConicBundle::BundleData (see ConicBundle::SumBlockModel)
    @version 1.0
    @date 2014-07-22
    @author Christoph Helmberg
*/


#include "SumBundle.hxx"

namespace ConicBundle {

  /** @ingroup InternalBundleModel

  */
  //@{


  /** @brief base class for use with SumBlockModel for storing and managing
      essential data of evaluations and the bundle describing the model

      This and all derived classes serve to store the essential part of the
      model, so that saving and restoring this should allow to continue from
      exactly the same model again. This is for supporting branch
      and bound approaches and the like.

      It also provides some essential common data and routines that are needed
      for consistent handling of problem changes or contributions to a common
      bundle model for sums of convex functions as required by the abstract
      SumBlockModel extending the BundleModel.

      In order to support changes on the fly while keeping the number of
      recomputations at a minimum, all modifications of the function or the
      model lead to appropriate changes in certain id counters.

      In particular, a nonnegative @a modification_id will be increased whenever the
      function changes. In order to check validity of the function values of the
      candidate and the center, appropriate ids of the arguments (a point id -
      this id being assigned by BundleSolver and not here - and the modification id)
      are stored at evaluation together with the relative precision, so that
      validity of the stored results can be checked by comparing the ids to the
      current values.  A negative id will always signal an uninitiliazed or
      outdated number. This is used, e.g., if information is discarded by
      clear_model().

      The handling of the model aggregate is a bit different, because its
      validity and availablility depends on the changes of the cutting
      model, which may have many different causes. Furthermore the aggregate
      may consist of different parts due to contriubtions in the SumBundle.
      The @a aggregate_id is increased whenever a new @a aggregate is
      computed (all parts will be computed at the same time).
      The overall aggregate can be marked invalid by setting @a aggregate to
      empty (no need to increase the aggregate id at such an occasion).

      The part stored in @a local_aggregate represents the model parts
      handled by the model itself or by its descendants, i.e., it includes the
      local part, the root parts of its SumBundleHandler and the local
      parts of its descendants. It does not contain the contributions of
      itself and its descendants to the SumBundle of a parent. In particular,
      if the entire model is contributed to the SumBundle of a parent,
      @a local_aggregate will be empty.

      The full @a aggregate contains all parts, i.e., on top of the
      local_aggregate it includes the aggregate of the contributions
      of itself and its descendants. In essence, the @a aggregate
      gives the full aggregate of the entire model representing this
      funciton with all its descendants.

  */

  class BundleData : public VariableMetricBundleData, public virtual CBout {
  private:

    friend class SumBlockModel;
    friend class ConeModel;
    friend class SumModel;
    friend class AFTModel;
    friend class AFTData;

    //===== variables =====
    CH_Matrix_Classes::Integer dim; ///< the incoming dimension, once it is detected, -1 if not yet clear (some need it, but it is mainly for consistency checks and should be initialized in eval_function, lb_function, recompute_center, and updated in apply_modification)
    FunctionTask function_task;   ///< see FunctionTask
    CH_Matrix_Classes::Real function_factor;  ///< >0, interpreted according to function_task
    CH_Matrix_Classes::Integer modification_id; ///< increased on any modification applied to the function  or @a function_factor

    CH_Matrix_Classes::Integer center_ub_mid; ///< modification id for which center_ub was computed, <0 if not initialized
    CH_Matrix_Classes::Real center_ub;     ///< upper bound on function value for modification id center_mid and center with point id center_id (does not include the effect of @a function_factor and @a function_task)
    CH_Matrix_Classes::Real center_relprec; ///< relative precision used in computing center_ub
    CH_Matrix_Classes::Integer center_id;  ///< point id of the stored center

    CH_Matrix_Classes::Integer aggregate_id;  ///< if available and ids match, the aggregate need not be recomputed
    MinorantPointer aggregate; ///< the aggregate linear minorant of the entire model (including @a function_factor and @a function_task)

    /// true if the local model is initialized and contributes to the model (maybe in addition to the sumbundle) 
    bool local_model_active;

    /// for switching/splitting between a local model and a common SumModel this stores the (potential) contribution to a common model
    SumBundle sumbundle;

    CH_Matrix_Classes::Integer max_old_minorants; ///< ==0 by default, otherwise it gives the maximum number of the most recent candidate minorants to be stored 
    MinorantBundle old_minorants;  ///< minorants collected e.g. for variable metric, holds the at most max_old_minorants most recent candidate minorants
    CH_Matrix_Classes::Integer minorant_nexti; ///< points to next position for storing in old_minorants (in a cyclic fashion)  

    CH_Matrix_Classes::Real model_curvature;  ///< this value >= 0 gives an estimate of the "model curvature" in cand_y (may be zero)

    ///if max_old_minorants > 0, it adds the minorant cyclically to old_minorants keeping max_old_minorants of them
    virtual int store_old_minorant(MinorantPointer& minorant);

    CH_Matrix_Classes::Integer prex_id; ///< counter for keeping track of PriamlExtender calls

  protected:

    CH_Matrix_Classes::Integer cand_ub_mid; ///< modification id for which cand_ub was computed, <0 if not initialized
    CH_Matrix_Classes::Real cand_ub;     ///< upper bound on function value for modification id cand_mid and candidate with point id cand_id (does not include the effect of @a function_factor and @a function_task)
    CH_Matrix_Classes::Real cand_relprec; ///< relative precision used in computing cand_ub
    CH_Matrix_Classes::Integer cand_id;  ///< point id of the stored candidate


    MinorantPointer local_aggregate; ///< the local (possibly empty) aggregate linear minorant of the model handled by this and its descendants. It includes @a function_factor and @a function_task, and comprises the local_model parts and the active sumbundle parts


    /// the purpose or mode the function is used in, see FunctionTask
    virtual FunctionTask& set_function_task() {
      return function_task;
    }

    /// the factor by which the minorants (except for the aggregate) need to be multiplied in order to match the current function scaling
    virtual CH_Matrix_Classes::Real& set_function_factor() {
      return function_factor;
    }

    /// the modification_id when center was computed; -1 indicates not initialized
    virtual CH_Matrix_Classes::Integer& set_center_ub_mid() {
      return center_ub_mid;
    }

    /// the aggregate minorant as currently in use in the bundle method; it already includes the function_factor
    virtual MinorantPointer& set_aggregate() {
      return aggregate;
    }

    /// the identification number of the aggregate minorant
    virtual CH_Matrix_Classes::Integer& set_aggregate_id() {
      return aggregate_id;
    }

    ///return the sumbundle
    virtual SumBundle& set_sumbundle() {
      return sumbundle;
    }

    ///return old_minorants
    virtual MinorantBundle& set_old_minorants() {
      return old_minorants;
    }

    ///for keeping track of changes by primal extenders
    virtual CH_Matrix_Classes::Integer get_prex_id() const {
      return prex_id;
    }



  public:

    //===== methods =====
    /// calls clear()
    BundleData();

    /// 
    virtual ~BundleData();

    /// reset to initial state (also used by the default constructor) 
    virtual void clear(CH_Matrix_Classes::Integer start_modification_id = 0);

    ///initialize from other BundleData
    virtual int init(const BundleData* bd);

    ///return a pointer to a clone of this
    BundleData* clone() const;

    ///returns the observed dimension; intially this is negative until the first point is worked on 
    virtual CH_Matrix_Classes::Integer get_dim() const {
      return dim;
    }

    /// the purpose or mode the function is used in, see FunctionTask
    virtual FunctionTask get_function_task() const {
      return function_task;
    }

    /// the factor by which the minorants (except for the aggregate) need to be multiplied in order to match the current function scaling
    virtual CH_Matrix_Classes::Real get_function_factor() const {
      return function_factor;
    }

    ///the value is negative if not initialized and otherwise gives the number of modifications applied to the function
    virtual CH_Matrix_Classes::Integer get_modification_id() const {
      return modification_id;
    }

    ///the value is negative if not initialized and otherwise gives the number of modifications applied to the function
    virtual CH_Matrix_Classes::Integer get_center_modification_id() const {
      return center_ub_mid;
    }

    ///the value of the upper bound stored for center (without function_factor) if its modification id is up to date
    virtual CH_Matrix_Classes::Real get_center_ub() const {
      return center_ub;
    }

    ///the point id for which center was computed if its modification id is up to date
    virtual CH_Matrix_Classes::Integer get_center_id() const {
      return center_id;
    }

    ///the value is negative if not initialized and otherwise gives the number of modifications applied to the function
    virtual CH_Matrix_Classes::Integer get_cand_modification_id() const {
      return cand_ub_mid;
    }

    ///the value of the upper bound stored for the candidate (without function_factor) if its modification id is up to date
    virtual CH_Matrix_Classes::Real get_cand_ub() const {
      return cand_ub;
    }

    ///the point id for which cand was computed if its modification id is up to date
    virtual CH_Matrix_Classes::Integer get_cand_id() const {
      return cand_id;
    }

    ///returns the bundle collected in store_old_minorants()
    const MinorantBundle& get_old_minorants() const {
      return old_minorants;
    }

    ///if the candidate information is available and consitent for point_id, copy it from cand to center and return 0, otherwise return 1
    virtual int do_step(CH_Matrix_Classes::Integer point_id);

    ///synchronize ids in any case by discarding inconsistent parts but return number of errors
    virtual int synchronize_ids(CH_Matrix_Classes::Integer& new_center_ub_fid,
      CH_Matrix_Classes::Integer new_center_id,
      CH_Matrix_Classes::Integer old_center_id,
      CH_Matrix_Classes::Integer& new_cand_ub_fid,
      CH_Matrix_Classes::Integer new_cand_id,
      CH_Matrix_Classes::Integer old_cand_id,
      CH_Matrix_Classes::Integer& new_aggregate_id,
      CH_Matrix_Classes::Integer new_prex_id = 0);

    ///check whether center computation is still valid for this point and modification id
    virtual bool center_modified(CH_Matrix_Classes::Integer& center_ub_fid,
      CH_Matrix_Classes::Integer check_center_id);

    ///check whether aggregate is available and has the same id
    virtual bool model_aggregate_modified(CH_Matrix_Classes::Integer model_aggregate_id);

    ///clear the cutting model and all function evaluations; for special cutting models of support functions it may sometimes be possible and sufficient to keep the model and its values but to regenerate the minorants when needed; if in this case only the current minorants should be deleted, set discard_minorants_only to true.  
    virtual void clear_model(bool discard_minorants_only = false);

    ///delete all kinds of aggregates but keep explicit parts of the cutting model
    virtual void clear_aggregates();

    ///see the last argument of FunctionOracle::evaluate() 
    virtual int call_primal_extender(PrimalExtender&, bool include_candidates = true);

    ///rearrange/extend the minorants according to the given groundset modifications 
    virtual int apply_modification(const GroundsetModification&, MinorantExtender* mex);

    ///return the PrimalData corresponding to the aggregate
    virtual const PrimalData* get_approximate_primal() const {
      return 0;
    }

    ///return the PrimalData corresponding to the aggregate
    virtual const PrimalData* get_center_primal() const {
      return 0;
    }

    ///return the PrimalData delivered by the last call of FunctionOracle::evaluate()
    virtual const PrimalData* get_candidate_primal() const {
      return 0;
    }

    ///return the sumbundle
    virtual const SumBundle& get_sumbundle() const {
      return sumbundle;
    }

    ///return the sumbundle
    virtual bool get_local_model_active() const {
      return local_model_active;
    }

    ///return a rough curvature estimate to judge the relevance of providing a good model for this oracle 
    virtual CH_Matrix_Classes::Real get_model_curvature() const {
      return model_curvature;
    }

    /// the latest minorants available; the number may fall below or exceed the one requested in request_n_latest_minorants(); the minorants still need to be mutliplied by function_factor
    virtual int get_latest_minorants(MinorantBundle& latest_minorants,
      CH_Matrix_Classes::Integer max_number);

    /// the minorants currently used in the model; the list may be empty or max contain other minorants than returned in get_latest_minorants(); the minorants still need to be mutliplied by function_factor
    virtual int get_model_data(MinorantBundle& model_minorants,
      CH_Matrix_Classes::Matrix& model_coeff) const;

    /// the aggregate minorant as currently in use in the bundle method; it already includes the function_factor
    virtual const MinorantPointer& get_aggregate() const {
      return aggregate;
    }


  };


  //@}

}

#endif

