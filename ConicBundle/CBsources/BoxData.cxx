/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/BoxData.cxx
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
#include "BoxData.hxx"
//for testing:
#include "MatrixCBSolver.hxx"

using namespace CH_Matrix_Classes;

namespace ConicBundle {

  // *****************************************************************************
  //                              BoxData
  // *****************************************************************************

  BoxData::BoxData(Real fun_factor, FunctionTask fun_task) {
    assert(fun_factor > 0.);
    set_function_factor() = fun_factor;
    set_function_task() = fun_task;
    prex_id = 0;
    clear();
  }

  // *****************************************************************************
  //                              ~BoxData
  // *****************************************************************************

  BoxData::~BoxData() {
    clear();
  }

  // *****************************************************************************
  //                                clear
  // *****************************************************************************

  void BoxData::clear(Integer start_modification_id) {
    BundleData::clear(start_modification_id);

    center_minorant.clear();
    center_boxvec.init(0, 1, 0.);

    cand_minorant.clear();
    cand_boxvec.init(0, 1, 0.);

    aggr_boxvec.init(0, 1, 0.);
    aggr_scaleval = 0.;

    boxbundle.clear();
    boxbundle_coeff.init(0, 1, 0.);
    boxbundle_coords.init(0, 1, Integer(0));
    boxbundle_complvalues.init(0, 1, 0.);
    primal_boxtapia.init(0, 1, 0.);

    nncbundle.clear();
    nncbundle_coeff.init(0, 1, 0.);
    primal_nnctapia.init(0, 1, 0.);
    outside_bundle = false;

    coord_switching.init(0, 1, 0.);
  }


  // *****************************************************************************
  //                                init
  // *****************************************************************************

  int BoxData::init(const BundleData* bd) {
    const BoxData* fbd = dynamic_cast<const BoxData*>(bd);
    if (fbd == 0) {
      if (cb_out())
        get_out() << "**** ERROR BoxData::init(): dynamic cast failed, argument is not of type const BoxData*" << std::endl;
      return 1;
    }
    clear();
    BundleData::init(bd);

    center_minorant = fbd->center_minorant;
    center_boxvec.init(fbd->center_boxvec);

    cand_minorant = fbd->cand_minorant;
    cand_boxvec.init(fbd->cand_boxvec);

    aggr_boxvec.init(fbd->aggr_boxvec);
    aggr_scaleval = fbd->aggr_scaleval;

    boxbundle = fbd->boxbundle;
    boxbundle_coeff = fbd->boxbundle_coeff;
    boxbundle_coords = fbd->boxbundle_coords;
    boxbundle_complvalues = fbd->boxbundle_complvalues;
    primal_boxtapia = fbd->primal_boxtapia;

    nncbundle = fbd->nncbundle;
    nncbundle_coeff = fbd->nncbundle_coeff;
    primal_nnctapia = fbd->primal_nnctapia;
    outside_bundle = fbd->outside_bundle;

    coord_switching = fbd->coord_switching;


    return 0;
  }

  // *****************************************************************************
  //                                init
  // *****************************************************************************

  BundleData* BoxData::clone() const {
    BoxData* fbd = new BoxData;
    fbd->init(this);
    return fbd;
  }

  // *****************************************************************************
  //                            synchronize_ids
  // *****************************************************************************

  int BoxData::synchronize_ids(Integer& new_center_ub_fid,
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
      prex_id);

    cand_minorant.synchronize_ids(get_modification_id(),
      new_center_id, old_center_id,
      new_cand_id, old_cand_id,
      prex_id);

    for (unsigned i = 0; i < boxbundle.size(); i++) {
      boxbundle[i].synchronize_ids(get_modification_id(),
        new_center_id, old_center_id,
        new_cand_id, old_cand_id,
        prex_id);
    }
    for (unsigned i = 0; i < nncbundle.size(); i++) {
      nncbundle[i].synchronize_ids(get_modification_id(),
        new_center_id, old_center_id,
        new_cand_id, old_cand_id,
        prex_id);
    }
    return err;
  }

  // *****************************************************************************
  //                              do_step
  // *****************************************************************************


  int BoxData::do_step(Integer point_id) {

    int retval = BundleData::do_step(point_id);
    if (retval) {
      if (cb_out())
        get_out() << "\n**** ERROR BoxData::do_step(.): BundleData::do_step(.) returned " << std::endl;
      return retval;
    }
    assert(cand_minorant.valid());

    center_minorant = cand_minorant;
    center_boxvec.init(cand_boxvec);

    return retval;
  }

  // *****************************************************************************
  //                                   clear_model
  // *****************************************************************************

  void BoxData::clear_model(bool discard_minorants_only) {
    /// here all data is unaffected as all data relates to the box
    if (!discard_minorants_only) {
      BundleData::clear_model();

      center_boxvec.init(0, 1, 0.);
      center_minorant.clear();

      cand_boxvec.init(0, 1, 0.);
      cand_minorant.clear();

      aggr_boxvec.init(0, 1, 0.);
      aggr_scaleval = 0.;

      boxbundle.clear();
      boxbundle_coeff.init(0, 1, 0.);
      boxbundle_coords.init(0, 1, Integer(0));
      boxbundle_complvalues.init(0, 1, 0.);
      primal_boxtapia.init(0, 1, 0.);

      nncbundle.clear();
      nncbundle_coeff.init(0, 1, 0.);
      primal_nnctapia.init(0, 1, 0.);
      outside_bundle = false;

      coord_switching.init(0, 1, 0.);
    }
  }


  // *****************************************************************************
  //                                   clear_aggregates
  // *****************************************************************************

  void BoxData::clear_aggregates() {
    BundleData::clear_aggregates();

    unsigned cnt = 0;
    Real coeffval = 0.;
    Indexmatrix delind;
    for (unsigned i = 0; i < nncbundle.size(); i++) {
      if (!nncbundle[i].aggregate()) {
        if (cnt < i) {
          nncbundle[cnt] = nncbundle[i];
        }
        cnt++;
      } else {
        delind.concat_below(Integer(i));
        coeffval += nncbundle_coeff(Integer(i));
      }
    }
    nncbundle.resize(cnt);
    nncbundle_coeff.delete_rows(delind, true);
    if (cnt > 0) {
      nncbundle_coeff(0) += coeffval;
    }

    aggr_boxvec.init(0, 1, 0.);
    aggr_scaleval = 0.;

    cnt = 0;

  }


  // *****************************************************************************
  //                               call_primal_extender
  // *****************************************************************************

  int BoxData::call_primal_extender(PrimalExtender& prex, bool include_candidates) {
    int err = 0;
    err = BundleData::call_primal_extender(prex, include_candidates);

    if (local_aggregate.call_primal_extender(prex, prex_id)) {
      if (cb_out())
        get_out() << "**** WARNING: BoxData::call_primal_extender(..): PrimalExtender::extend failed for local_aggregate" << std::endl;
      err++;
    }

    if (center_minorant.call_primal_extender(prex, prex_id)) {
      if (cb_out())
        get_out() << "**** WARNING: BoxData::call_primal_extender(..): PrimalExtender::extend failed for center_minorant" << std::endl;
      err++;
    }

    if (cand_minorant.call_primal_extender(prex, prex_id)) {
      if (cb_out())
        get_out() << "**** WARNING: BoxData::call_primal_extender(..): PrimalExtender::extend failed for cand_minorant" << std::endl;
      err++;
    }

    for (unsigned int i = 0; i < boxbundle.size(); i++) {
      if (boxbundle[i].call_primal_extender(prex, prex_id)) {
        if (cb_out())
          get_out() << "**** WARNING: BoxData::call_primal_extender(..): PrimalExtender::extend failed for bundle minorant number " << i << std::endl;
        err++;
      }
    }

    for (unsigned int i = 0; i < nncbundle.size(); i++) {
      if (nncbundle[i].call_primal_extender(prex, prex_id)) {
        if (cb_out())
          get_out() << "**** WARNING: BoxData::call_primal_extender(..): PrimalExtender::extend failed for bundle minorant number " << i << std::endl;
        err++;
      }
    }

    BoxPrimalExtender* meprex = dynamic_cast<BoxPrimalExtender*>(&prex);
    if (meprex == 0) {
      if (cb_out())
        get_out() << "**** WARNING: BoxData::call_primal_extender(..): cast to BoxPrimalExtender failed" << std::endl;
      err++;
    } else {
      if ((get_center_modification_id() == get_modification_id()) && (meprex->extend_Box(center_boxvec))) {
        if (cb_out())
          get_out() << "**** WARNING: BoxData::call_primal_extender(..):  extending center_boxvec failed" << std::endl;
        err++;
      }

      if ((get_cand_modification_id() == get_modification_id()) && (meprex->extend_Box(cand_boxvec))) {
        if (cb_out())
          get_out() << "**** WARNING: BoxData::call_primal_extender(..):  extending cand_boxvec failed" << std::endl;
        err++;
      }

      if ((aggr_boxvec.dim() > 0) && (meprex->extend_Box(aggr_boxvec))) {
        if (cb_out())
          get_out() << "**** WARNING: BoxData::call_primal_extender(..):  extending aggr_boxvec failed" << std::endl;
        err++;
      }

    }

    return err;
  }

  // *****************************************************************************
  //                               apply_modification
  // *****************************************************************************

  int BoxData::apply_modification(const GroundsetModification& gsmdf, MinorantExtender* mex) {
    ///modification_id has already been increased
    int err = BundleData::apply_modification(gsmdf, mex);
    Integer modification_id = get_modification_id();

    if ((center_minorant.valid()) && (center_minorant.apply_modification(gsmdf, modification_id, mex))) {
      if (cb_out())
        get_out() << "**** WARNING: BoxData::apply_modification(..): apply_modification failed for center_minorant" << std::endl;
      err++;
    }

    if (cand_minorant.apply_modification(gsmdf, modification_id, mex)) {
      if (cb_out())
        get_out() << "**** WARNING: BoxData::apply_modification(..): apply_modification failed for cand_minorant" << std::endl;
      err++;
    }

    for (unsigned int i = 0; i < nncbundle.size(); i++) {
      if (nncbundle[i].apply_modification(gsmdf, modification_id, mex)) {
        if (cb_out())
          get_out() << "**** WARNING: BoxData::apply_modification(..): apply_modification failed for bundle minorant number " << i << std::endl;
        err++;
      }
    }

    if (get_sumbundle().has_bundle_for(get_function_task())) {
      if (set_sumbundle().apply_modification(gsmdf, modification_id, mex, get_function_task())) {
        if (cb_out())
          get_out() << "**** WARNING: BoxData::apply_modification(..):  sumbundle.apply_modificaiton(..) failed" << std::endl;
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
    get_out()<<"**** WARNING: BoxData::apply_modification(..):  test aggregation failed"<<std::endl;
      }
      if (!mnrtp.equals(aggregate)){
        if (cb_out())
    get_out()<<"**** WARNING: BoxData::apply_modification(..):  test aggregation delivered different aggregate"<<std::endl;

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

  const PrimalData* BoxData::get_approximate_primal() const {
    return &aggr_boxvec;
  }

  // *****************************************************************************
  //                                get_center_primal
  // *****************************************************************************

  const PrimalData* BoxData::get_center_primal() const {
    return &center_boxvec;
  }

  // *****************************************************************************
  //                                get_candidate_primal
  // *****************************************************************************

  const PrimalData* BoxData::get_candidate_primal() const {
    assert(cand_minorant.valid());
    return &cand_boxvec;
  }


}

