/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/AFTData.hxx
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



#ifndef CONICBUNDLE_AFTDATA_HXX
#define CONICBUNDLE_AFTDATA_HXX


/**  @file AFTData.hxx
    @brief Header declaring the class ConicBundle::AFTData (needed for ConicBundle::ConvexProblem)
    @version 1.0
    @date 2014-07-23
    @author Christoph Helmberg
*/


#include "BundleData.hxx"

namespace ConicBundle {

  /** @ingroup InternalBundleModel

  */
  //@{

  /** @brief specialized BundleData class for use with AFTModel  for storing and managing essential data of evaluations and the bundle describing the model

      The general philosophy regarding ids is the same as described in
      BundleData, but this time the @a modification_id is increased
      whenever the AFT is modified or it detects that the function it
      transforms was modified.  Furthermore, the AFTModel also
      generates new points, namly the transformed arguments arising from
      the input points.  In fact, whenever the AFT is modified,
      typically a new argument needs to be generated for the same old
      point, so the AFTModle has to assign its own @a aftpoint_id
      whenever it generates a new point by way of transforming an
      argument. Therefore the data has to store for each input point_id
      the transformed aft_point together with a point_mid (storing what
      AFT it was generated for) and a point_pid (storing the new point
      identifier of the AFTModel). Whenever a point is transformed, the
      affine part dependending on this point only is stored in
      point_offset. If the function is then evaluated, the returned
      value is called point_ub with point_ub_mid giving the
      modification_id when it was determined.

      Likewise the situation with the aggregate is a bit more involved.  The
      validity of the aggregate now depends on changes in the function as well
      as that of the AFT. So whenever the AFT changes but the function does not,
      the aggregate is in principle still available but has changed, so an extra
      @a aggregate_id is needed here again and the id assigned by the function
      model is stored in @a aggregate_mid.  The role of @a aggregate_available
      is as before.  Typically there is no need to have an explicit copy of the
      aggregate available here, it suffices to transform the one available in
      the function model. In some cases it is hard to avoid storing it
      explicitly and if it is done, this is marked by setting @a
      aggregate_stored to true.

  */

  class AFTData :public BundleData {
  public:
    //===== variables =====

    //FunctionTask function_task;   ///< see FunctionData::FunctionTask
    //CH_Matrix_Classes::Real function_factor;  ///< >0, interpreted according to function_t //CH_Matrix_Classes::Integer modification_id; ///< increased on any modification applied to the AffineFunctionTransformation  
    CH_Matrix_Classes::Integer aftpoint_id; ///< increased whenever "a new point" (neither cand nor center) enters from outside or is created by changes in aft
    CH_Matrix_Classes::Integer function_id; ///< the most recent (largest) function id known (-1 on initializaton or if model==0)

    //CH_Matrix_Classes::Integer center_ub_mid; ///< modification_id for which center_ub was computed, <0 if not initialized
    //CH_Matrix_Classes::Real center_ub;     ///< (untransformed) upper bound on function value returned for aft_center with aftpoint_id aft_center_id
    //CH_Matrix_Classes::Real center_relprec;     ///< relative precision by which it was computed
    CH_Matrix_Classes::Integer center_ub_fid; ///< the function id returned on its computation  
    CH_Matrix_Classes::Integer center_mid; ///< modification_id for which center_offset and aft_center were computed, <0 if not initialized
    CH_Matrix_Classes::Integer center_pid; ///< aftpoint_id of the AFT-transformed center
    //CH_Matrix_Classes::Integer center_id;  ///< point id of the generating center
    CH_Matrix_Classes::Matrix aft_center; ///< AFT-transformed center
    CH_Matrix_Classes::Real center_offset;     ///< constant plus linear y term for modification id center_mid and center with point id center_id

    //CH_Matrix_Classes::Integer cand_ub_mid; ///< modification_id for which cand_ub was computed, <0 if not initialized
    //CH_Matrix_Classes::Real cand_ub;     ///< (untransformed) upper bound on function value returned for aft_cand with aftpoint_id aft_cand_id
    //CH_Matrix_Classes::Real cand_relprec;     ///< relative precision by which it was computed
    CH_Matrix_Classes::Integer cand_ub_fid; ///< the function id returned on its computation  
    CH_Matrix_Classes::Integer cand_mid; ///< modification_id for which cand_offset and aft_cand were computed, <0 if not initialized
    CH_Matrix_Classes::Integer cand_pid; ///< aftpoint_id of the AFT-transformed cand
    //CH_Matrix_Classes::Integer cand_id;  ///< point id of the generating cand
    CH_Matrix_Classes::Matrix aft_cand; ///< AFT-transformed cand
    CH_Matrix_Classes::Real cand_offset;     ///< constant plus linear y term for modification id center_mid and cand with point id cand_id


    CH_Matrix_Classes::Integer aggregate_fid; ///< the function's aggregate id for which this aggregate was computed (this was returned by model for get_aggregate_id)
    //MinorantPointer aggregate; ///< the aggregate linear minorant of the model (including @a function_factor and @a function_task)
    //CH_Matrix_Classes::Integer aggregate_id;  ///< the aggregate id assigned when computed (but not nec. stored) here, at the same time as aggregate_available is set to true
    bool aggregate_available; ///< true if still valid, false if outdated or not initialized


    //===== methods =====
    ///
    virtual ~AFTData();

    /// reset to initial state (also used by the default constructor) 
    virtual void clear(CH_Matrix_Classes::Integer start_modification_id = 0);

    /// calls clear()
    AFTData(CH_Matrix_Classes::Integer start_modification_id = 0) {
      clear(start_modification_id);
    }

    ///initialize from other BundleData
    virtual int init(const BundleData* bd);

    ///return a pointer to a clone of this
    BundleData* clone() const;

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
    virtual bool center_modified(CH_Matrix_Classes::Integer& function_id,
      CH_Matrix_Classes::Integer center_id);

    ///check whether aggregate is available and has the same id
    virtual bool model_aggregate_modified(CH_Matrix_Classes::Integer last_aggr_id);
    ///clear the cutting modle and all function evaluatoins
    virtual void clear_model(bool discard_minorants_only = false);

    ///delete all kinds of aggregates but keep explicit parts of the cutting model
    virtual void clear_aggregates();

    // as in BundleData 
    //virtual int call_primal_extender(PrimalExtender&,bool /*include_candidates*/=true){return 1;}

    // as in BundleData
    //virtual int apply_modification(const GroundsetModification&,MinorantExtender* mex){return 1;}

    // as in BundleData
    //virtual int get_approximate_primal(PrimalData& ) const {return 1;}

    // as in BundleData
    //virtual int get_center_primal(PrimalData& ) const {return 1;}

    //as in BundleData
    //virtual int get_candidate_primal(PrimalData& ) const {return 1;}

    //as in BundleData
    //virtual const SumBundle& get_sumbundle() const {return sumbundle;}

  };


  //@}

}

#endif

