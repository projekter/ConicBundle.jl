/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/NNCData.cxx
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
#include "NNCData.hxx"
//for testing:
//#include "MatrixCBSolver.hxx"

using namespace CH_Matrix_Classes;

namespace ConicBundle {

  // *****************************************************************************
  //                              NNCData
  // *****************************************************************************

  NNCData::NNCData(Real fun_factor, FunctionTask fun_task) {
    assert(fun_factor > 0.);
    set_function_factor() = fun_factor;
    set_function_task() = fun_task;
    clear();
  }

  // *****************************************************************************
  //                              ~NNCData
  // *****************************************************************************

  NNCData::~NNCData() {
    clear();
  }

  // *****************************************************************************
  //                                clear
  // *****************************************************************************

  void NNCData::clear(Integer start_modification_id) {
    BundleData::clear(start_modification_id);

    center_minorant.clear();
    cand_minorant.clear();
    cand_minorants.clear();

    bundlecoeff.init(0, 0, 0.);
    primal_tapia.init(0, 0, 0.);

    bundle.clear();
    outside_bundle = false;
  }


  // *****************************************************************************
  //                                init
  // *****************************************************************************

  int NNCData::init(const BundleData* bd) {
    const NNCData* fbd = dynamic_cast<const NNCData*>(bd);
    if (fbd == 0) {
      if (cb_out())
        get_out() << "**** ERROR NNCData::init(): dynamic cast failed, argument is not of type const NNCData*" << std::endl;
      return 1;
    }
    clear();
    BundleData::init(bd);

    center_minorant = fbd->center_minorant;
    cand_minorant = fbd->cand_minorant;
    cand_minorants = fbd->cand_minorants;
    bundlecoeff = fbd->bundlecoeff;
    primal_tapia = fbd->primal_tapia;
    bundle = fbd->bundle;
    outside_bundle = fbd->outside_bundle;

    return 0;
  }

  // *****************************************************************************
  //                                init
  // *****************************************************************************

  BundleData* NNCData::clone() const {
    NNCData* fbd = new NNCData;
    fbd->init(this);
    return fbd;
  }

  // *****************************************************************************
  //                            synchronize_ids
  // *****************************************************************************

  int NNCData::synchronize_ids(Integer& new_center_ub_fid,
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
    for (unsigned i = 0; i < cand_minorants.size(); i++) {
      cand_minorants[i].synchronize_ids(get_modification_id(),
        new_center_id, old_center_id,
        new_cand_id, old_cand_id,
        get_prex_id());
    }
    for (unsigned i = 0; i < bundle.size(); i++) {
      bundle[i].synchronize_ids(get_modification_id(),
        new_center_id, old_center_id,
        new_cand_id, old_cand_id,
        get_prex_id());
    }
    return err;
  }

  // *****************************************************************************
  //                              do_step
  // *****************************************************************************


  int NNCData::do_step(Integer point_id) {

    int retval = BundleData::do_step(point_id);
    if (retval) {
      if (cb_out())
        get_out() << "\n**** ERROR NNCData::do_step(.): BundleData::do_step(.) returned " << std::endl;
      return retval;
    }
    assert(cand_minorant.valid());

    center_minorant = cand_minorant;

    return retval;
  }

  // *****************************************************************************
  //                                   clear_model
  // *****************************************************************************

  void NNCData::clear_model(bool discard_minorants_only) {
    BundleData::clear_model(discard_minorants_only);

    center_minorant.clear();
    cand_minorant.clear();
    cand_minorants.clear();

    bundlecoeff.init(0, 0, 0.);
    primal_tapia.init(0, 0, 0.);

    bundle.clear();
    outside_bundle = false;
  }


  // *****************************************************************************
  //                                   clear_aggregates
  // *****************************************************************************

  void NNCData::clear_aggregates() {
    BundleData::clear_aggregates();

    unsigned cnt = 0;
    Real coeffval = 0.;
    Indexmatrix delind;
    for (unsigned i = 0; i < bundle.size(); i++) {
      if (!bundle[i].aggregate()) {
        if (cnt < i) {
          bundle[cnt] = bundle[i];
        }
        cnt++;
      } else {
        delind.concat_below(Integer(i));
        coeffval += bundlecoeff(Integer(i));
      }
    }
    bundle.resize(cnt);
    bundlecoeff.delete_rows(delind, true);
    if (cnt > 0) {
      bundlecoeff(0) += coeffval;
    }

    cnt = 0;

  }


  // *****************************************************************************
  //                               call_primal_extender
  // *****************************************************************************

  int NNCData::call_primal_extender(PrimalExtender& prex, bool include_candidates) {
    int err = 0;
    err = BundleData::call_primal_extender(prex, include_candidates);

    if (local_aggregate.call_primal_extender(prex, get_prex_id())) {
      if (cb_out())
        get_out() << "**** WARNING: NNCData::call_primal_extender(..): PrimalExtender::extend failed for local_aggregate" << std::endl;
      err++;
    }

    if (center_minorant.call_primal_extender(prex, get_prex_id())) {
      if (cb_out())
        get_out() << "**** WARNING: NNCData::call_primal_extender(..): PrimalExtender::extend failed for center_minorant" << std::endl;
      err++;
    }

    if (cand_minorant.call_primal_extender(prex, get_prex_id())) {
      if (cb_out())
        get_out() << "**** WARNING: NNCData::call_primal_extender(..): PrimalExtender::extend failed for cand_minorant" << std::endl;
      err++;
    }

    for (unsigned int i = 0; i < cand_minorants.size(); i++) {
      if (cand_minorants[i].call_primal_extender(prex, get_prex_id())) {
        if (cb_out())
          get_out() << "**** WARNING: NNCData::call_primal_extender(..): PrimalExtender::extend failed for cand_minorants number " << i << std::endl;
        err++;
      }
    }

    for (unsigned int i = 0; i < bundle.size(); i++) {
      if (bundle[i].call_primal_extender(prex, get_prex_id())) {
        if (cb_out())
          get_out() << "**** WARNING: NNCData::call_primal_extender(..): PrimalExtender::extend failed for bundle minorant number " << i << std::endl;
        err++;
      }
    }

    return err;
  }

  // *****************************************************************************
  //                               apply_modification
  // *****************************************************************************

  int NNCData::apply_modification(const GroundsetModification& gsmdf, MinorantExtender* mex) {
    ///modification_id has already been increased
    int err = BundleData::apply_modification(gsmdf, mex);

    if ((center_minorant.valid()) && (center_minorant.apply_modification(gsmdf, get_modification_id(), mex))) {
      if (cb_out())
        get_out() << "**** WARNING: NNCData::apply_modification(..): apply_modification failed for center_minorant" << std::endl;
      err++;
    }

    if ((cand_minorant.valid()) && (cand_minorant.apply_modification(gsmdf, get_modification_id(), mex))) {
      if (cb_out())
        get_out() << "**** WARNING: SOCData::apply_modification(..): apply_modification failed for cand_minorant" << std::endl;
      err++;
    }

    for (unsigned int i = 0; i < cand_minorants.size(); i++) {
      if (cand_minorants[i].apply_modification(gsmdf, get_modification_id(), mex)) {
        if (cb_out())
          get_out() << "**** WARNING: NNCData::apply_modification(..): apply_modification failed for cand_minorants number " << i << std::endl;
        err++;
      }
    }

    for (unsigned int i = 0; i < bundle.size(); i++) {
      if (bundle[i].apply_modification(gsmdf, get_modification_id(), mex)) {
        if (cb_out())
          get_out() << "**** WARNING: NNCData::apply_modification(..): apply_modification failed for bundle minorant number " << i << std::endl;
        err++;
      }
    }

    if (get_sumbundle().has_bundle_for(get_function_task())) {
      if (set_sumbundle().apply_modification(gsmdf, get_modification_id(), mex, get_function_task())) {
        if (cb_out())
          get_out() << "**** WARNING: NNCData::apply_modification(..):  sumbundle.apply_modificaiton(..) failed" << std::endl;
        err++;
      }
    }

    //----------- testing
    /*
    if ((aggregate.valid())&&(bundle.size()>0)){

      Matrix aggrv(gsmdf.new_vardim(),1,0.);
      Real offset;
      aggregate.get_minorant(offset,aggrv,0);
      std::cout<<"offset="<<offset<<" aggr="<<aggrv;
      PrimalMatrix prim;
      if (aggregate.get_primal(prim)==0)
        std::cout<<"prim="<<prim;

      MinorantPointer mnrtp;
      if (mnrtp.aggregate(bundle,bundlecoeff)){
        if (cb_out())
    get_out()<<"**** WARNING: NNCData::apply_modification(..):  test aggregation failed"<<std::endl;
      }
      if (!mnrtp.equals(aggregate)){
        if (cb_out())
    get_out()<<"**** WARNING: NNCData::apply_modification(..):  test aggregation delivered different aggregate"<<std::endl;

        Matrix taggrv(gsmdf.new_vardim(),1,0.);
        Real toffset;
        mnrtp.get_minorant(toffset,taggrv,0);
        std::cout<<"toffset="<<toffset<<" taggr="<<taggrv;
        PrimalMatrix tprim;
        if (mnrtp.get_primal(tprim)==0)
    std::cout<<"tprim="<<tprim;

        std::cout<<"coeff="<<bundlecoeff;

        for (unsigned i=0;i<bundle.size();i++){
    Matrix bunv(gsmdf.new_vardim(),1,0.);
    Real bunoff;
    bundle[i].get_minorant(bunoff,bunv,0);
    std::cout<<"bunoff["<<i<<"]="<<bunoff<<" bunv["<<i<<"]="<<bunv;
    PrimalMatrix bunprim;
    if (bundle[i].get_primal(bunprim)==0)
      std::cout<<"bunprim["<<i<<"]="<<bunprim;

        }
      }
    }
    */

    return err;

  }

  // *****************************************************************************
  //                                get_approximate_primal
  // *****************************************************************************

  const PrimalData* NNCData::get_approximate_primal() const {
    return const_cast<MinorantPointer*>(&get_aggregate())->get_primal();
  }

  // *****************************************************************************
  //                                get_center_primal
  // *****************************************************************************

  const PrimalData* NNCData::get_center_primal() const {
    return const_cast<MinorantPointer*>(&center_minorant)->get_primal();
  }

  // *****************************************************************************
  //                                get_candidate_primal
  // *****************************************************************************

  const PrimalData* NNCData::get_candidate_primal() const {
    assert(cand_minorant.valid());
    return const_cast<MinorantPointer*>(&(cand_minorant))->get_primal();
  }


}

