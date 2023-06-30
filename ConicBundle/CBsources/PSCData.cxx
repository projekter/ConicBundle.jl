/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/PSCData.cxx
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



#include "mymath.hxx"
#include "PSCData.hxx"

using namespace CH_Matrix_Classes;

namespace ConicBundle {

  // *****************************************************************************
  //                              PSCData
  // *****************************************************************************

  PSCData::PSCData(Real fun_factor, FunctionTask fun_task) {
    assert(fun_factor > 0.);
    set_function_factor() = fun_factor;
    set_function_task() = fun_task;
    clear();
  }

  // *****************************************************************************
  //                              ~PSCData
  // *****************************************************************************

  PSCData::~PSCData() {
    clear();
  }

  // *****************************************************************************
  //                                clear
  // *****************************************************************************

  void PSCData::clear(Integer start_modification_id) {
    BundleData::clear(start_modification_id);

    center_Ritzval = 0.;
    center_Ritzvec.init(0, 0, 0.);
    center_minorant.clear();

    cand_Ritzval.init(0, 0, 0.);
    cand_Ritzvec.init(0, 0, 0.);
    cand_minorant.clear();

    primaleigs.init(0, 0, 0.);
    primalvecs.init(0, 0, 0.);
    primal_Ritzval = 0;
    primal_aggregate.clear();
    primal_aggregate_coeff = 0.;
    bundle_aggregate.clear();
    growthrate = 0;
    primalgrowth.init(0, 0, 0.);
    dualgrowth.init(0, 0, 0.);

    bundlevecs.init(0, 0, 0.);

    topvecs.init(0, 0, 0.);
    Ritz_values.init(0, 0, 0.);
    activedim = 0;
    keepsize = 0;
    skippedsize = 0;

    scaling_gap = 0.;
    center_gap = 0.;
  }


  // *****************************************************************************
  //                                init
  // *****************************************************************************

  int PSCData::init(const BundleData* bd) {
    const PSCData* mbd = dynamic_cast<const PSCData*>(bd);
    if (mbd == 0) {
      if (cb_out())
        get_out() << "**** ERROR PSCData::init(): dynamic cast failed, argument is not of type const PSCData*" << std::endl;
      return 1;
    }
    clear();
    BundleData::init(bd);

    center_Ritzval = mbd->center_Ritzval;
    center_Ritzvec = mbd->center_Ritzvec;
    center_minorant = mbd->center_minorant;

    cand_Ritzval = mbd->cand_Ritzval;
    cand_Ritzvec = mbd->cand_Ritzvec;
    cand_minorant = mbd->cand_minorant;


    primaleigs = mbd->primaleigs;
    primalvecs = mbd->primalvecs;
    primal_Ritzval = mbd->primal_Ritzval;
    primal_aggregate = mbd->primal_aggregate;
    primal_aggregate_coeff = mbd->primal_aggregate_coeff;
    bundle_aggregate = mbd->bundle_aggregate;
    growthrate = mbd->growthrate;
    primalgrowth = mbd->primalgrowth;
    dualgrowth = mbd->dualgrowth;

    bundlevecs = mbd->bundlevecs;

    topvecs = mbd->topvecs;
    Ritz_values = mbd->Ritz_values;
    activedim = mbd->activedim;
    keepsize = mbd->keepsize;
    skippedsize = mbd->skippedsize;

    scaling_gap = mbd->scaling_gap;
    center_gap = mbd->center_gap;

    return 0;
  }

  // *****************************************************************************
  //                                init
  // *****************************************************************************

  BundleData* PSCData::clone() const {
    PSCData* mbd = new PSCData;
    mbd->init(this);
    return mbd;
  }

  // *****************************************************************************
  //                            synchronize_ids
  // *****************************************************************************

  int PSCData::synchronize_ids(Integer& new_center_ub_fid,
    Integer new_center_id,
    Integer old_center_id,
    Integer& new_cand_ub_fid,
    Integer new_cand_id,
    Integer old_cand_id,
    Integer& new_aggregate_id,
    Integer new_prex_id) {
    int err = BundleData::synchronize_ids(new_center_ub_fid, new_center_id, old_center_id,
      new_cand_ub_fid, new_cand_id, old_cand_id, new_aggregate_id, new_prex_id);

    center_minorant.synchronize_ids(get_modification_id(),
      new_center_id, old_center_id,
      new_cand_id, old_cand_id,
      get_prex_id());
    cand_minorant.synchronize_ids(get_modification_id(),
      new_center_id, old_center_id,
      new_cand_id, old_cand_id,
      get_prex_id());
    primal_aggregate.synchronize_ids(get_modification_id(),
      new_center_id, old_center_id,
      new_cand_id, old_cand_id,
      get_prex_id());
    bundle_aggregate.synchronize_ids(get_modification_id(),
      new_center_id, old_center_id,
      new_cand_id, old_cand_id,
      get_prex_id());

    return err;
  }

  // *****************************************************************************
  //                              do_step
  // *****************************************************************************


  int PSCData::do_step(Integer point_id) {

    int retval = BundleData::do_step(point_id);
    if (retval) {
      if (cb_out())
        get_out() << "\n**** ERROR PSCData::do_step(.): BundleData::do_step(.) returned " << std::endl;
      return retval;
    }

    assert(topvecs.coldim() > 0);
    assert(cand_Ritzvec.coldim() > 0);
    center_Ritzval = Ritz_values(0);
    center_Ritzvec = topvecs.col(0);
    center_minorant = cand_minorant;

    center_gap = 0.;

    return retval;
  }

  // *****************************************************************************
  //                                   clear_model
  // *****************************************************************************

  void PSCData::clear_model(bool discard_minorants_only) {
    center_minorant.clear();
    cand_minorant.clear();
    set_old_minorants().clear();
    clear_aggregates();

    BundleData::clear_model(discard_minorants_only);

    if (!discard_minorants_only) {

      center_Ritzval = 0.;
      center_Ritzvec.init(0, 0, 0.);

      cand_Ritzval.init(0, 0, 0.);
      cand_Ritzvec.init(0, 0, 0.);

      primaleigs.init(0, 0, 0.);
      primalvecs.init(0, 0, 0.);
      primal_Ritzval = 0;
      primal_aggregate_coeff = 0.;
      growthrate = 0.;
      primalgrowth.init(0, 0, 0.);
      dualgrowth.init(0, 0, 0.);

      bundlevecs.init(0, 0, 0.);
      topvecs.init(0, 0, 0.);
      Ritz_values.init(0, 0, 0.);
      activedim = 0;
      keepsize = 0;
      skippedsize = 0;

      scaling_gap = 0.;
    }
  }

  // *****************************************************************************
  //                                   clear_aggregates
  // *****************************************************************************

  void PSCData::clear_aggregates() {
    BundleData::clear_aggregates();

    //before deleting the aggregates, keep the trace at the same value
    if ((primal_aggregate_coeff > 1e-10 * get_function_factor()) && (primaleigs.dim() > 0)) {
      primaleigs(0) += primal_aggregate_coeff; //this is usually the most important vector
    }

    primal_aggregate_coeff = 0.;
    primal_aggregate.clear();
    bundle_aggregate.clear();

  }

  // *****************************************************************************
  //                               call_primal_extender
  // *****************************************************************************

  int PSCData::call_primal_extender(PrimalExtender& prex, bool include_candidates) {
    int err = BundleData::call_primal_extender(prex, include_candidates);

    if (local_aggregate.call_primal_extender(prex, get_prex_id())) {
      if (cb_out())
        get_out() << "**** WARNING: FunctionData::call_primal_extender(..): PrimalExtender::extend failed for local_aggregate" << std::endl;
      err++;
    }

    center_minorant.clear(); //can be reconstructed any time 

    cand_minorant.clear(); //can be reconstructed any time

    if (primal_aggregate.call_primal_extender(prex, get_prex_id())) {
      if (cb_out())
        get_out() << "**** WARNING: PSCData::call_primal_extender(..): PrimalExtender::extend failed for primal_aggregate" << std::endl;
      err++;
    }

    if (bundle_aggregate.call_primal_extender(prex, get_prex_id())) {
      if (cb_out())
        get_out() << "**** WARNING: PSCData::call_primal_extender(..): PrimalExtender::extend failed for bundle_aggregate" << std::endl;
      err++;
    }

    PSCPrimalExtender* meprex = dynamic_cast<PSCPrimalExtender*>(&prex);
    if (meprex == 0) {
      if (cb_out())
        get_out() << "**** WARNING: PSCData::call_primal_extender(..): cast to PSCPrimalExtender failed" << std::endl;
      err++;
    } else {
      if ((get_center_modification_id() == get_modification_id()) && (meprex->extend_Ritz(center_Ritzvec))) {
        if (cb_out())
          get_out() << "**** WARNING: PSCData::call_primal_extender(..):  extending center_Ritzvec failed" << std::endl;
        err++;
      }

      if ((cand_ub_mid == get_modification_id()) && (meprex->extend_Ritz(cand_Ritzvec))) {
        if (cb_out())
          get_out() << "**** WARNING: PSCData::call_primal_extender(..):  extending cand_Ritzvec failed" << std::endl;
        err++;
      }

      if ((bundlevecs.coldim() > 0) && (meprex->extend_Ritz(bundlevecs))) {
        if (cb_out())
          get_out() << "**** WARNING: PSCData::call_primal_extender(..):  extending bundlevecs failed" << std::endl;
        err++;
      }

      if ((topvecs.coldim() > 0) && (meprex->extend_Ritz(topvecs))) {
        if (cb_out())
          get_out() << "**** WARNING: PSCData::call_primal_extender(..):  extending topvecs failed" << std::endl;
        err++;
      }
    }

    return err;

  }

  // *****************************************************************************
  //                               apply_modification
  // *****************************************************************************

  int PSCData::apply_modification(const GroundsetModification& gsmdf, MinorantExtender* mex) {
    ///modification_id has already been increased
    int err = BundleData::apply_modification(gsmdf, mex);

    if ((center_minorant.valid()) && (center_minorant.apply_modification(gsmdf, get_modification_id(), mex))) {
      if (cb_out())
        get_out() << "**** WARNING: PSCData::apply_modification(..): apply_modification failed for center_minorant" << std::endl;
      err++;
    }

    if ((cand_minorant.valid()) && (cand_minorant.apply_modification(gsmdf, get_modification_id(), mex))) {
      if (cb_out())
        get_out() << "**** WARNING: PSCData::apply_modification(..): apply_modification failed for center_minorant" << std::endl;
      err++;
    }

    if ((primal_aggregate.valid()) && (primal_aggregate.apply_modification(gsmdf, get_modification_id(), mex))) {
      if (cb_out())
        get_out() << "**** WARNING: PSCData::apply_modification(..): apply_modification failed for primal_aggregate" << std::endl;
      err++;
    }

    if ((bundle_aggregate.valid()) && (bundle_aggregate.apply_modification(gsmdf, get_modification_id(), mex))) {
      if (cb_out())
        get_out() << "**** WARNING: PSCData::apply_modification(..): apply_modification failed for bundle_aggregate" << std::endl;
      err++;
    }

    if (get_sumbundle().has_bundle_for(get_function_task())) {
      if (set_sumbundle().apply_modification(gsmdf, get_modification_id(), mex, get_function_task())) {
        if (cb_out())
          get_out() << "**** WARNING: PSCData::apply_modification(..):  sumbundle.apply_modificaiton(..) failed" << std::endl;
        err++;
      }
    }

    return err;
  }


  // *****************************************************************************
  //                                get_approximate_primal
  // *****************************************************************************

  const PrimalData* PSCData::get_approximate_primal() const {
    return const_cast<MinorantPointer*>(&get_aggregate())->get_primal();
  }

  // *****************************************************************************
  //                                get_center_primal
  // *****************************************************************************

  const PrimalData* PSCData::get_center_primal() const {
    return const_cast<MinorantPointer*>(&center_minorant)->get_primal();
  }

  // *****************************************************************************
  //                                get_candidate_primal
  // *****************************************************************************

  const PrimalData* PSCData::get_candidate_primal() const {
    return const_cast<MinorantPointer*>(&cand_minorant)->get_primal();
  }



}

