/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/SOCModel.cxx
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
#include "SOCModel.hxx"
#include "BundleIdProx.hxx"
#include "SOCModelParameters.hxx"

using namespace CH_Tools;
using namespace CH_Matrix_Classes;

namespace ConicBundle {

  // *****************************************************************************
  //                                get_minorant
  // *****************************************************************************

  //Compute the subgradient corresponding to the matrix P*Diag(d)*P'/sum(d)
  //and store it in a new column. P has orthonormal columns and d is 
  //a nonengative vector of appropriate size.

  int SOCModel::get_minorant(MinorantPointer& mp,
    const Matrix& socvec) {
    assert((socvec.coldim() == 1) && (socvec.rowdim() >= 1));

    int err = 0;
    Minorant* mnrt = oracle->generate_minorant(socvec);
    if ((mnrt == 0) || (!mnrt->offset_gives_value_at_origin())) {
      if (cb_out())
        get_out() << "**** ERROR SOCModel::get_minorant(....): oracle->generate_minorant failed to generate a minorant or returned one with offset_gives_value_at_origin()==false" << std::endl;
      return 1;
    }
    mp.init(mnrt, data.get_modification_id());

    return err;
  }


  // *****************************************************************************
  //                            clear()
  // *****************************************************************************

  void SOCModel::clear() {
    ConeModel::clear();
    data.clear();
    delete model_selection;
    model_selection = new SOCModelParameters;

    ret_code = 0;
    nr_eval = 0;
    delete block;
    block = 0;
  }

  // *****************************************************************************
  //                            SOCModel
  // *****************************************************************************

  SOCModel::SOCModel(SOCOracle* fo,
    Real fun_factor,
    FunctionTask fun_task,
    CBout* cb, int cbinc) :
    CBout(cb, cbinc), ConeModel(cb, cbinc), oracle(fo), data(fun_factor, fun_task), model_selection(0), block(0) {
    assert(fo);
    clear();
  }

  // *****************************************************************************
  //                            SOCModel
  // *****************************************************************************

  SOCModel::~SOCModel() {
    delete model_selection;
    delete block;
  }

  // *****************************************************************************
  //                            evaluate_oracle
  // *****************************************************************************


  int SOCModel::evaluate_oracle(Integer y_id,
    const Matrix& y,
    Real nullstep_bound,
    Real relprec) {
    assert(relprec > 0.);

    CH_Tools::Microseconds start_eval = clock.time();

    //delete old info
    Real cmpval = ((data.cand_ub_mid >= 0) ? std::fabs(data.cand_ub) : 1.);
    data.cand_minorant.clear();
    data.cand_relprec = relprec;
    data.cand_id = y_id;
    data.cand_ub_mid = -1;  //signals no valid computation yet


    //---- evaluate  	
    int err = 0;
    Real SOCval_bound = nullstep_bound / data.get_function_factor();
    Real SOCval_relprec = 0.1 * relprec * (cmpval + 1. / data.get_function_factor()) / (cmpval + 1.);
    data.cand_SOCval = SOCval_bound;

    Matrix tmpmat;
    SOCPrimalExtender* pep = 0;
    nr_eval++;

    preeval_time += clock.time() - start_eval;
    start_eval = clock.time();

    ret_code = oracle->evaluate(y, SOCval_relprec, data.cand_SOCval,
      data.cand_SOCvec, pep);
    eval_time += clock.time() - start_eval;
    start_eval = clock.time();


    if (data.cand_SOCvec.coldim() != 1) {
      if (cb_out())
        get_out() << "**** ERROR: SOCModel::eval_function(): function has to return one SOC_vector but returned " << data.cand_SOCvec.coldim() << std::endl;
      err++;
    } else if (data.cand_SOCvec.rowdim() == 0) {
      if (cb_out())
        get_out() << "**** ERROR: SOCModel::eval_function(): the returned SOC_vector has 0 rows but it must have at least 1" << std::endl;
      err++;
    } else if (data.cand_SOCvec(0) != 1.) {
      if (cb_out())
        get_out() << "**** ERROR: SOCModel::eval_function(): SOC_vector(0)=" << data.cand_SOCvec(0) << "!= 1." << std::endl;
      err++;
    }

    //--- if needed, apply the modifications to the old data
    if (pep) {
      if (call_primal_extender(*pep)) {
        if (cb_out())
          get_out() << "**** WARNING: SOCModel::eval_function(): call_primal_extender failed" << std::endl;
        err++;
      }
      delete pep;
    }

    if ((data.get_center_modification_id() == data.get_modification_id()) && (data.center_SOCvec.coldim() == 1) && (data.center_SOCvec.rowdim() != data.cand_SOCvec.rowdim())) {
      if (cb_out())
        get_out() << "**** ERROR: SOCModel::eval_function(): order of the matrix changed without having called SOCModel::clear_model first" << std::endl;
      err++;
    }

    if (err) {
      posteval_time += clock.time() - start_eval;
      return err;
    }

    //--- update the most important subspace (thereby maybe correcting the function value)
    data.cand_ub = data.cand_SOCval;
    data.cand_ub_mid = data.get_modification_id();

    //update old SOCvecs 
    data.SOCvecs_max = max(data.SOCvecs_max, model_selection->get_max_bundle_size());
    if (data.store_SOCvec(data.cand_SOCvec)) {
      if (cb_out())
        get_out() << "**** WARNING SOCModel::eval_function(): data.store_SOCvec() failed" << std::endl;
    }

    if (ret_code) {
      if (cb_out()) get_out() << "**** WARNING: SOCModel::eval_function(): function returned code = " << ret_code << std::endl;
      posteval_time += clock.time() - start_eval;
      return ret_code;
    }

    posteval_time += clock.time() - start_eval;
    return err;
  }

  // *****************************************************************************
  //                              evaluate_local_model
  // *****************************************************************************

  //evaluate the current cutting model in $y$ 
  //if evaluated by an iterative method that provides upper and lower bounds
  //it may stop when the lower bound (lb) is above the nullstep_bound
  //evalutation may also stop, if (ub-lb)<relprec*(|ub|+1) is satisfied

  int SOCModel::evaluate_local_model(Real& lb,
    Integer /* y_id */,
    const Matrix& y,
    Real relprec) {

    int model_ret_code = oracle->evaluate_projection(y, data.bundlevecs, relprec, lb);
    if (model_ret_code) {
      if (cb_out())
        get_out() << "**** WARNING SOCModel::eval_model: oracle->evaluate_projection failed and returned " << model_ret_code << std::endl;
    }

    return model_ret_code;
  }


  // *****************************************************************************
  //                                update_local_model
  // *****************************************************************************

  int SOCModel::update_local_model(ModelUpdate model_update,
    Integer center_id,
    const Matrix& center_y,
    Integer y_id,
    const Matrix& y,
    Real model_maxviol,
    BundleProxObject& H) {
    assert((data.primalvec.coldim() == 1) || (data.bundlevecs.coldim() == 0));
    assert((data.SOCvecs.coldim() == 0) || (data.cand_SOCvec.rowdim() == data.SOCvecs.rowdim() + 1));

    if ((model_update != new_subgradient) && (data.cand_ub_mid != data.get_modification_id())) {
      if (cb_out())
        get_out() << "**** ERROR SOCModel::update_local_model: candidate required but its values are not up to date" << std::endl;
      return 1;
    }

    if (model_update == descent_step) {
      data.do_step(y_id);
    }


    //initialize the bundle if necessary
    if (data.bundlevecs.coldim() == 0) {
      data.form_bundlevecs(max(2, model_selection->get_max_model_size() - 1));
      return 0;
    }

    if ((!data.get_local_model_active()) || (model_update == new_subgradient)) {
      //it is a root or a child, all is taken care of by sumbundle
      //but if desired here the local bundle could be updated along, maybe with more 
      return 0;
    }

    //--- if the trace is not fixed and there is no active vector, do nothing
    if ((data.get_function_task() != ObjectiveFunction) && (data.cand_SOCval < 0.)) {
      if (cb_out(1)) {
        get_out() << " SOC: penaltyskip bundle_update " << std::endl;
      }
      return 0;
    }

    int retval = model_selection->select_model(data.bundlevecs,
      data.primalvec,
      data.cand_SOCval,
      data.cand_SOCvec,
      data.center_SOCval,
      data.center_SOCvec,
      data.SOCvecs,
      oracle,
      data.get_function_task(),
      data.get_function_factor(),
      model_update,
      center_id,
      center_y,
      y_id,
      y,
      model_maxviol,
      H);

    if (data.primalvec.coldim() == 0) {
      //generate a new one
      if ((data.get_function_task() != ObjectiveFunction) && (data.cand_SOCval < 0)) {
        data.primalvec.init(data.cand_SOCvec.rowdim(), 1, 0.);
      } else {
        data.primalvec.init(data.cand_SOCvec, data.get_function_factor());
      }
    }

    if ((retval) && (cb_out())) {
      get_out() << "**** WARNING SOCModel::update_local_model: model_selection->select_model() returned " << retval << std::endl;
    }

    if (cb_out(1))
      get_out() << "  SOCModel::update_local_model(...): bundlesize=" << data.bundlevecs.coldim() << std::endl;

    return retval;
  }


  // *****************************************************************************
  //                        provide_local_aggregate
  // *****************************************************************************


  int SOCModel::provide_local_aggregate(Integer /* y_id */, const Matrix& /* iny */) {
    if (data.local_aggregate.valid()) {
      return 0;
    }

    //--- form an aggregate 
    data.local_aggregate.clear();

    assert(data.get_local_model_active());
    assert(data.primalvec.coldim() == 1);
    assert((data.get_function_task() != ObjectiveFunction) || (std::fabs(data.primalvec(0) - data.get_function_factor()) < 1e-6 * data.get_function_factor()));

    if (get_minorant(data.local_aggregate, data.primalvec)) {
      if (cb_out())
        get_out() << "**** ERROR SOCModel::provide_local_aggregate: get_minorant failed for the aggregate" << std::endl;
      return 1;
    }

    return 0;
  }


  // *****************************************************************************
  //                         switch_to_local
  // *****************************************************************************

  int SOCModel::switch_to_local(const MinorantBundle& current_model,
    const Matrix& current_coeff) {
    assert(Integer(current_model.size()) == current_coeff.rowdim());
    data.local_aggregate.clear();
    if (data.local_aggregate.aggregate(current_model, current_coeff)) {
      if (cb_out()) {
        get_out() << "*** WARNING: SOCModel::sumbundle_mode(): getting the bundle_aggregate failed" << std::endl;
      }
    }
    if (oracle->extract_SOCvector(data.primalvec, data.local_aggregate.get_minorant())) {
      if (cb_out()) {
        get_out() << "*** WARNING: SOCModel::sumbundle_mode(): getting the bundle_aggregate failed" << std::endl;
      }
    }
    assert((data.get_function_task() != ObjectiveFunction) || (std::fabs(data.primalvec(0) - data.get_function_factor()) < 1e-6 * data.get_function_factor()));
    data.form_bundlevecs(max(2, model_selection->get_max_model_size() - 1));
    return 0;
  }

  // *****************************************************************************
  //                        get_aggregate_and_coefficient
  // *****************************************************************************

  int SOCModel::get_aggregate_and_coefficient(MinorantPointer& sumbundle_aggr,
    Real& aggr_coeff) {
    //set new model if required
    sumbundle_aggr.clear();
    if (data.primalvec.coldim() == 1) {
      aggr_coeff = data.primalvec(0);
      if (aggr_coeff > 1e-10) {
        get_minorant(sumbundle_aggr, data.primalvec);
        sumbundle_aggr.scale(1. / aggr_coeff);
      } else {
        assert(data.bundlevecs.coldim() > 1);
        Integer k = data.bundlevecs.rowdim();
        Matrix tmpvec(k + 1, 1); chk_set_init(tmpvec, 1);
        tmpvec(0) = 1.;
        mat_xey(k, tmpvec.get_store() + 1, data.bundlevecs.get_store());
        get_minorant(sumbundle_aggr, tmpvec);
      }
    } else {
      assert(data.cand_SOCvec.coldim() == 1);
      aggr_coeff = data.get_function_factor();
      get_minorant(sumbundle_aggr, data.cand_SOCvec);
    }
    return 0;
  }


  // *****************************************************************************
  //                                start_local_augmodel
  // *****************************************************************************


  int SOCModel::start_local_augmodel(QPModelDataPointer& qpblockp,
    Integer /* cand_id */,
    const Matrix& /* cand_y */,
    const Indexmatrix* needed_indexsubset) {
    if (qpblockp.set_model_data(block)) {
      delete block;
      block = qpblockp.generate_conemodel_data();
      block->set_cbout(this);
    }
    assert(block);

    //use the local model
    Integer xdim = data.bundlevecs.coldim() + 1;

    MinorantBundle bundle;
    bundle.resize(unsigned(xdim));
    unsigned cnt = 0;
    Matrix tmpvec;
    Matrix tmpmat;
    if (oracle->projection(tmpvec, tmpmat, data.bundlevecs, needed_indexsubset)) {
      if (cb_out())
        get_out() << "**** ERROR: SOCModel::start_augmodel(): oracle->projection() failed" << std::endl;
      return 1;
    }
    if (needed_indexsubset == 0) {
      //assert(tmpmat.rowdim()==cand_y.rowdim());
      for (Integer i = 0; i < xdim; i++) {
        bundle[cnt++].init(new Minorant(true, tmpvec(i), tmpmat.rowdim(), tmpmat.get_store() + i * tmpmat.rowdim()), data.get_modification_id());
      }
    } else {
      Integer dim = needed_indexsubset->dim();
      assert(tmpmat.rowdim() == dim);
      for (Integer i = 0; i < xdim; i++) {
        bundle[cnt++].init(new Minorant(true, tmpvec(i), dim, tmpmat.get_store() + i * dim, needed_indexsubset->get_store()), data.get_modification_id());
      }
    }
    assert(Integer(cnt) == xdim);

    Indexmatrix tmpind(0, 1, Integer(0));
    Matrix dummymat(0, 1, 0.);
    block->init(MinorantPointer(), bundle, 0, Indexmatrix(1, 1, xdim), tmpind, dummymat, dummymat, data.get_function_factor(), data.get_function_task());

    qpblockp.set_model_data(block);
    assert(qpblockp.get_model_data_ptr() == block);

    if (cb_out(2)) {
      get_out() << " SOCxdim=" << xdim;
    }

    return 0;
  }


  // *****************************************************************************
  //                                make_local_aggregate
  // *****************************************************************************

  int SOCModel::make_local_aggregate(Real& new_function_factor,
    bool keep_penalty_fixed) {
    assert(new_function_factor == data.get_function_factor());

    bool penalty_changed = false;
    data.local_aggregate.clear();
    assert(block);

    if ((!keep_penalty_fixed) && (data.get_function_task() == AdaptivePenaltyFunction)) {
      // increasing the penalty only makes sense if the constraint is active
      // and the model does not promise a strong reduction anyways.
      Real prec = 0;
      Real ubmodelval = block->tracedual(&prec);
      Real traceval = block->evaluate_trace();
      if ((ubmodelval > 100. * prec) &&
        (data.get_center_modification_id() == data.get_modification_id()) &&
        (ubmodelval > .9 * data.get_center_ub())
        ) {
        new_function_factor = max(2 * data.get_function_factor(), data.get_function_factor() + ubmodelval);
        penalty_changed = true;
      } else if ((ubmodelval < 0.01 * prec) &&
        (data.get_center_modification_id() == data.get_modification_id()) &&
        (data.get_function_factor() > 1.) &&
        (10. * traceval < data.get_function_factor())
        ) {
        new_function_factor = max(1., max(0.01 * data.get_function_factor(), 5. * traceval));
        penalty_changed = true;
      }

      if (penalty_changed) {
        block->adjust_trace(new_function_factor);
      }
    }

    if (!penalty_changed) {
      Matrix tmpvec;
      Real soc_activity;
      if (block->get_socx(0, tmpvec, &soc_activity)) {
        if (cb_out()) {
          get_out() << "**** ERROR SOCModel::make_local_aggregate(): retrieving the SOC part of the QP solution failed" << std::endl;
        }
        return 1;
      }

      // generate data.primalvec
      data.primalvec.newsize(data.bundlevecs.rowdim() + 1, 1); chk_set_init(data.primalvec, 1);
      data.primalvec(0) = tmpvec(0);
      tmpvec.delete_rows(Indexmatrix(1, 1, 0), true);
      Matrix tmpv;
      genmult(data.bundlevecs, tmpvec, tmpv);
      mat_xey(tmpv.rowdim(), data.primalvec.get_store() + 1, tmpv.get_store());
      // generate the corresponding minorant
      if (get_minorant(data.local_aggregate, data.primalvec)) {
        if (cb_out()) {
          get_out() << "*** WARNING: SOCModel::make_local_aggregate(): get_minorant failed when forming the aggregate" << std::endl;
        }
      }

      // if (cb_out(2)){
      //   get_out()<<" x0="<<data.primalvec(0)<<std::endl;
      // }

    }
    return 0;
  }

  // *****************************************************************************
  //                                 lb_function
  // *****************************************************************************

  //returns a *quick* lower bound for the function value at y
  //(eg by a previous subgradient)

  Real SOCModel::lb_function(Integer y_id, const Matrix& y) {
    Real lb = CB_minus_infinity;

    if (data.get_aggregate().valid()) {
      lb = data.get_aggregate().evaluate(y_id, y);
    } else if ((data.get_center_modification_id() == data.get_modification_id()) && (provide_center_minorant().valid())) {
      lb = data.get_function_factor() * data.center_minorant.evaluate(y_id, y);
    } else if ((data.cand_ub_mid == data.get_modification_id()) && (provide_cand_minorant().valid())) {
      lb = data.get_function_factor() * data.cand_minorant.evaluate(y_id, y);
    }

    if (data.get_function_task() != ObjectiveFunction)
      lb = max(0., lb);

    return lb;
  }

  // *****************************************************************************
  //                             provide_center_minorant
  // *****************************************************************************

  MinorantPointer& SOCModel::provide_center_minorant() {
    if (data.get_center_modification_id() != data.get_modification_id()) {
      if (cb_out())
        get_out() << "**** WARNING: SOCModel::provide_center_minorant(...): no center available" << std::endl;
      data.center_minorant.clear();
    } else if (!data.center_minorant.valid()) {
      data.center_minorant.clear();
      if (get_minorant(data.center_minorant, data.center_SOCvec)) {
        if (cb_out())
          get_out() << "**** WARNING: SOCModel::provide_center_minorant(...): get_minorant failed" << std::endl;
        data.center_minorant.clear();
      }
    }
    return data.center_minorant;
  }

  // *****************************************************************************
  //                             provide_cand_minorant
  // *****************************************************************************

  MinorantPointer& SOCModel::provide_cand_minorant() {
    if (data.cand_ub_mid != data.get_modification_id()) {
      if (cb_out())
        get_out() << "**** WARNING: SOCModel::provide_cand_minorant(...): no cand available" << std::endl;
      data.cand_minorant.clear();
    } else if (!data.cand_minorant.valid()) {
      data.cand_minorant.clear();
      if (get_minorant(data.cand_minorant, data.cand_SOCvec)) {
        if (cb_out())
          get_out() << "**** WARNING: SOCModel::provide_cand_minorant(...): get_minorant failed" << std::endl;
        data.cand_minorant.clear();
      }
    }
    return data.cand_minorant;
  }


  // *****************************************************************************
  //                            local_adjust_multiplier
  // *****************************************************************************


  int SOCModel::local_adjust_multiplier(Real& new_function_factor) {
    assert(new_function_factor == data.get_function_factor());

    //no multipliers of relevance in this class
    if (data.get_function_task() != AdaptivePenaltyFunction) {
      return 0;
    }

    assert(block);
    Real traceval = block->evaluate_trace();
    Real oldval = data.get_function_factor();
    new_function_factor = min(data.get_function_factor(), max(1., 1.5 * traceval));
    if (new_function_factor < oldval) {
      block->adjust_trace(new_function_factor);
    }

    return 0;
  }

  // *****************************************************************************
  //                               output_bundle_data
  // *****************************************************************************

  std::ostream& SOCModel::output_bundle_data(std::ostream& out) const {
    out << "\n%(begin SOCModel)\n";

    out << "\n n_soc=" << data.bundlevecs.rowdim() << ";";
    Integer xdim = data.bundlevecs.coldim() + 1;
    out << "\n xdim=" << xdim << ";";
    Integer ydim = data.get_dim();
    out << "\n ydim=" << ydim << ";";

    Matrix bigvec(xdim, 1, 0.);
    Matrix bigmat(ydim, xdim, 0.);
    oracle->projection(bigvec, bigmat, data.bundlevecs);

    out << "\n costs=[";
    for (Integer i = 0; i < xdim; i++) {
      out.width(18);
      out << " " << (bigvec)(i);
    }
    out << "]';";

    out << "\n subgmat=[";
    for (Integer i = 0; i < ydim; i++) {
      for (Integer j = 0; j < xdim; j++) {
        out.width(18);
        out << " " << (bigmat)(i, j);
      }
      if (i < ydim - 1)
        out << "\n";
    }
    out << "];";

    out << "\n%(end SOCModel)\n";
    return out;
  }

  // *****************************************************************************
  //                         set_bundle_parameters      
  // *****************************************************************************

  int SOCModel::set_bundle_parameters(const BundleParameters& bp) {
    const SOCModelParametersObject* mpo = dynamic_cast<const SOCModelParametersObject*>(&bp);
    int retval = 0;
    if (mpo) {
      delete model_selection;
      model_selection = dynamic_cast<SOCModelParametersObject*>(mpo->clone_BundleParameters());
      assert(model_selection);
    } else {
      assert(model_selection);
      retval = model_selection->init(bp);
    }
    if (model_selection->get_max_bundle_size() >= 2) {
      data.SOCvecs_max = max(data.SOCvecs_max, model_selection->get_max_bundle_size());
    }
    return retval;
  }


}
