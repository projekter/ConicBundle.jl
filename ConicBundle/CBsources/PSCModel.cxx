/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/PSCModel.cxx
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
#include "PSCModel.hxx"
#include "PSCModelParameters.hxx"
#include "PSCVariableMetricSelection.hxx"
#include "BundleIdProx.hxx"

using namespace CH_Tools;
using namespace CH_Matrix_Classes;

namespace ConicBundle {

  // *****************************************************************************
  //                            clear()
  // *****************************************************************************

  void PSCModel::clear() {
    ConeModel::clear();
    data.clear();
    delete model_selection;
    model_selection = new PSCModelParameters;

    delete vm_selection;
    PSCVariableMetricSelection* pscvms = new PSCVariableMetricSelection(this, 0);
    pscvms->set_oracle(oracle);
    vm_selection = pscvms;

    ret_code = 0;
    nr_eval = 0;
    delete block;
    block = 0;

  }

  // *****************************************************************************
  //                            PSCModel
  // *****************************************************************************


  int PSCModel::set_variable_metric_selection(VariableMetricSelection* vms) {
    delete vm_selection;
    vm_selection = vms;
    PSCVariableMetricSelection* pscvms = dynamic_cast<PSCVariableMetricSelection*>(vms);
    if (pscvms)
      pscvms->set_oracle(oracle);
    return 0;
  }


  // *****************************************************************************
  //                            PSCModel
  // *****************************************************************************

  PSCModel::PSCModel(PSCOracle* fo,
    Real fun_factor,
    FunctionTask fun_task,
    CBout* cb, int cbinc) :
    CBout(cb, cbinc), ConeModel(cb, cbinc), oracle(fo), data(fun_factor, fun_task), model_selection(0), block(0) {
    assert(fo);
    clear();

  }

  // *****************************************************************************
  //                            PSCModel
  // *****************************************************************************

  PSCModel::~PSCModel() {
    delete model_selection;
    delete block;
  }

  // *****************************************************************************
  //                            evaluate_oracle
  // *****************************************************************************


  int PSCModel::evaluate_oracle(Integer y_id,
    const Matrix& y,
    Real nullstep_bound,
    Real relprec) {
    assert(relprec > 0.);

    CH_Tools::Microseconds start_eval = clock.time();

    //delete old info

    data.cand_minorant.clear();
    data.cand_id = y_id;
    data.cand_relprec = relprec;
    data.cand_ub_mid = -1;


    //---- evaluate  	
    int err = 0;
    Real eigval_bound = nullstep_bound / data.get_function_factor();
    Real eigval_relprec = 0.1 * relprec;
    //Real cmpval=((data.cand_ub_mid>=0)?0.1*std::fabs(data.cand_ub):1.);
    //Real eigval_relprec=max(10*eps_Real,0.5*relprec/min(1.,(1.-relprec))*(cmpval+1./data.get_function_factor())/(cmpval+1.));
    nr_eval++;

    Matrix tmpmat;

    if (data.topvecs.coldim() > data.bundlevecs.coldim()) {
      if (data.activedim <= 0)
        tmpmat = data.topvecs;
      else
        tmpmat.init(data.topvecs.rowdim(), min(data.topvecs.coldim(), data.activedim), data.topvecs.get_store());
    } else {
      if (data.activedim <= 0)
        tmpmat = data.bundlevecs;
      else
        tmpmat.init(data.bundlevecs.rowdim(), min(data.bundlevecs.coldim(), data.activedim), data.bundlevecs.get_store());
    }

    PSCPrimalExtender* pep = 0;

    preeval_time += clock.time() - start_eval;
    start_eval = clock.time();

    ret_code = oracle->evaluate(y, tmpmat, eigval_relprec, eigval_bound,
      data.cand_Ritzvec, data.cand_Ritzval, pep);
    eval_time += clock.time() - start_eval;
    start_eval = clock.time();

    if (data.cand_Ritzval.dim() == 0) {
      if (cb_out())
        get_out() << "**** ERROR: PSCModel::eval_function(): function returned no Ritz_values" << std::endl;
      err++;
    }

    if (data.cand_Ritzval.dim() != data.cand_Ritzvec.coldim()) {
      if (cb_out())
        get_out() << "**** ERROR: PSCModel::eval_function(): function returned different number of Ritz_values and Ritz_vectors" << std::endl;
      err++;
    }

    //--- if needed, apply the modifications to the old data
    if (pep) {
      if (call_primal_extender(*pep)) {
        if (cb_out())
          get_out() << "**** WARNING: PSCModel::eval_function(): call_primal_extender failed" << std::endl;
        err++;
      }
      delete pep;
    }

    if (err) {
      posteval_time += clock.time() - start_eval;
      return err;
    }


    //--- update the most important subspace (thereby maybe correcting the function value)
    Real lmax = max(data.cand_Ritzval);

    //compute a basis of the subspace spanned by topvecs/bundlevecs and cand_Ritzvec into cand_Ritzvec
    if (data.topvecs.coldim() > 0)
      tmpmat = data.topvecs;
    else
      tmpmat = data.bundlevecs;
    tmpmat.concat_right(data.cand_Ritzvec);
    Indexmatrix piv;
    Integer r1 = tmpmat.QR_factor(piv);
    tmpmat.delete_cols(Range(r1, tmpmat.coldim() - 1));
    piv.init(Range(0, r1 - 1));
    data.topvecs.init(tmpmat.rowdim(), tmpmat.coldim(), 0.);
    for (Integer i = 0; i < data.topvecs.coldim(); i++)
      data.topvecs(i, i) = 1.;
    tmpmat.Q_times(data.topvecs, r1);

    //compute the eigenvalue decomposition of the projected matrix
    Matrix Ritz_vecs;
    int ret_val = oracle->evaluate_projection(y, data.topvecs, 1e-10, Ritz_vecs, data.Ritz_values);
    if (ret_val) {
      if (cb_out()) {
        get_out() << "**** WARNING: PSCModel::eval_function(): evaluate_projection returned code = " << ret_val << std::endl;
      }
    }
    genmult(data.topvecs, Ritz_vecs, tmpmat);
    swap(data.topvecs, tmpmat);

    if (cb_out(2))
      get_out() << " lmax=" << lmax << " Ritz_values=" << transpose(data.Ritz_values);

    lmax = data.Ritz_values(0);

    data.cand_ub = lmax + eigval_relprec * (fabs(lmax) + Real(Ritz_vecs.rowdim()) * eps_Real);
    data.cand_ub_mid = data.get_modification_id();

    if (ret_code) {
      if (cb_out()) get_out() << "**** WARNING: PSCModel::eval_function(): function returned code = " << ret_code << std::endl;
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

  int PSCModel::evaluate_local_model(Real& lb,
    Integer y_id,
    const Matrix& y,
    Real relprec) {
    assert(data.bundlevecs.coldim() > 0);

    //--- determine model value for bundle vectors

    Matrix tmpvec, tmpmat;
    int model_ret_code = oracle->evaluate_projection(y, data.bundlevecs, relprec, tmpmat, tmpvec);
    if (model_ret_code) {
      if (cb_out())
        get_out() << "**** WARNING PSCModel::eval_model: oracle->evaluate_projection failed and returned " << model_ret_code << std::endl;
    }
    lb = max(tmpvec);

    //--- determine model values for aggregates
    if (data.bundle_aggregate.valid()) {
      lb = max(lb, data.bundle_aggregate.evaluate(y_id, y));
    }

    if (cb_out(3)) {
      get_out().precision(12);
      get_out() << " PSCModval=" << lb;
    }

    return model_ret_code;
  }


  // *****************************************************************************
  //                                update_local_model
  // *****************************************************************************

  int PSCModel::update_local_model(ModelUpdate model_update,
    Integer center_id,
    const Matrix& center_y,
    Integer y_id,
    const Matrix& y,
    Real model_maxviol,
    BundleProxObject& H) {
    assert((data.primalvecs.coldim() > 0) || (data.bundlevecs.coldim() == 0));
    assert(data.primaleigs.dim() == data.primalvecs.coldim());
    assert(data.primalvecs.rowdim() == data.bundlevecs.rowdim());
    assert(data.cand_Ritzvec.rowdim() == data.topvecs.rowdim());
    assert(data.primalgrowth.dim() == data.primaleigs.dim());

    if ((model_update != new_subgradient) && (data.cand_ub_mid != data.get_modification_id())) {
      if (cb_out())
        get_out() << "**** ERROR PSCModel::update_local_model: candidate required but its values are not up to date" << std::endl;
      return 1;
    }

    if (!data.cand_minorant.empty()) {
      Real val = data.cand_minorant.evaluate(center_id, center_y);
      //if (data.get_center_ub()-val>data.center_gap)
      data.center_gap = data.get_center_ub() - val;
      if (data.scaling_gap > data.center_gap)
        data.scaling_gap = .9 * data.scaling_gap + .1 * data.center_gap;
      else
        data.scaling_gap = .5 * data.scaling_gap + .5 * data.center_gap;
      if (cb_out(3)) {
        get_out() << " PSC center_gap=" << data.center_gap << " scaling_gap=" << data.scaling_gap << " model_maxviol=" << model_maxviol / data.get_function_factor() << std::endl;
      }
    }

    if (model_update == descent_step) {
      data.center_gap = 0;
      data.do_step(y_id);
    }

    //topvecs holds the Ritz-vectors sorted nonincreasingly by Ritz-values

    //initialize the bundle if necessary
    if (data.bundlevecs.coldim() == 0) {
      assert(data.bundle_aggregate.empty());
      Integer nc = min(data.topvecs.coldim() - 1, max(3, model_selection->get_max_model_size()));
      Indexmatrix ind;
      data.bundlevecs.init(data.topvecs.cols(Range(0, nc - 1)));
      data.primal_Ritzval = data.Ritz_values(0);
      data.primalvecs = data.bundlevecs;
      data.primaleigs.init(nc, 1, 0.);
      data.growthrate = .1;
      data.primalgrowth.init(nc, 1, .1);
      data.dualgrowth.init(nc, 1, 1.);
      if ((data.get_function_task() == ObjectiveFunction) || (data.Ritz_values(0) >= 0.)) {
        data.primaleigs(0) = data.get_function_factor();
        data.primalgrowth(0) = 1.;
        data.dualgrowth(0) = .1;
      }
      data.activedim = data.keepsize = 1;

      return 0;
    }

    if ((!data.get_local_model_active()) || (model_update == new_subgradient)) {
      //it is a root or a child, all is taken care of by sumbundle
      //but if desired here the local bundle could be updated along, maybe with more 
      return 0;
    }

    //--- if the trace is not fixed and there is no active vector, do nothing
    if ((data.get_function_task() != ObjectiveFunction) && (data.Ritz_values(0) < 0.)) {
      if (cb_out(1)) {
        get_out() << " PSC penaltyskip bundle_update " << std::endl;
      }
      return 0;
    }

    //---- call the update routine
    if (data.primaleigs.coldim() == 0) {
      assert(data.Ritz_values.rowdim() > 0);
      data.primal_Ritzval = data.Ritz_values(0);
    }
    Real gap = max(model_maxviol, max(data.scaling_gap, data.center_gap) * data.get_function_factor());

    Real diffval;
    if (model_update == descent_step) {
      diffval = data.get_center_ub() - data.get_aggregate().evaluate(y_id, y) / data.get_function_factor();
    } else {
      diffval = data.get_center_ub() - data.get_aggregate().evaluate(center_id, center_y) / data.get_function_factor();
    }
    //assert(diffval>-1e-10*data.get_function_factor());
    diffval = max(diffval, 1e-10 * data.get_function_factor());

    int retval = model_selection->select_model(data.bundlevecs,
      data.bundle_aggregate,
      data.topvecs, data.Ritz_values,
      data.activedim, data.keepsize,
      data.skippedsize,
      data.primal_Ritzval,
      data.primaleigs,
      data.primalvecs,
      data.primal_aggregate,
      data.primal_aggregate_coeff,
      data.growthrate,
      data.primalgrowth,
      data.dualgrowth,
      data.cand_Ritzvec,
      data.cand_Ritzval,
      oracle,
      data.get_modification_id(),
      data.get_function_task(),
      data.get_function_factor(),
      model_update,
      center_id, center_y,
      y_id, y,
      //model_maxviol,
      gap,
      diffval,
      H);
    if ((retval) && (cb_out())) {
      get_out() << "**** ERROR PSCModel::update_local_model: model_selection->select_modle failed and returned " << retval << std::endl;
    }

    //---- aggregate the remaining primal vectors
    //check whether some of the primal vectors are not well reperesented in the new bundle
    if (data.keepsize < data.primaleigs.dim()) {
      Indexmatrix ind(Range(data.keepsize, data.primaleigs.dim() - 1));
      if (cb_out(2)) {
        get_out() << " PSC: aggrdim=" << ind.dim();
        get_out().precision(4);
      }

      //aggregate 
      if (ind.dim() < data.primalvecs.coldim()) {  //only a subset of the columns
        if (model_selection->get_minorant(data.primal_aggregate, data.primal_aggregate_coeff, data.primalvecs.cols(ind), data.primaleigs(ind), oracle, data.get_modification_id())) {
          if (cb_out())
            get_out() << "**** ERROR PSCModel::update_local_model: get_minorant failed for primal_aggregate" << std::endl;
        }
      } else { //aggregate all columns
        if (model_selection->get_minorant(data.primal_aggregate, data.primal_aggregate_coeff, data.primalvecs, data.primaleigs, oracle, data.get_modification_id())) {
          if (cb_out())
            get_out() << "**** ERROR PSCModel::update_model: get_minorant failed for primal_aggregate" << std::endl;
        }
      }
      data.primalvecs.delete_cols(ind);
      data.primaleigs.delete_rows(ind);
      data.primalgrowth.delete_rows(ind);
      data.dualgrowth.delete_rows(ind);
    }

    if (cb_out(1))
      get_out() << "  PSCModel::update_local_model(...): bundlesize=" << data.bundlevecs.coldim() << std::endl;


    //Test on the quality of the aggregate
    if (cb_out(2)) {
      Matrix tmpmat, tmp2;
      genmult(data.bundlevecs, data.primalvecs, tmp2, 1., 0., 1);
      genmult(data.bundlevecs, tmp2, tmpmat);
      MinorantPointer testmp;
      Real coeff = 0.;
      if (model_selection->get_minorant(testmp, coeff, tmpmat, data.primaleigs, oracle, data.get_modification_id())) {
        if (cb_out()) {
          get_out() << "*** WARNING: PSCModel::update_model(): get_minorant failed when forming the test aggregate" << std::endl;
        }
      }
      get_out() << " projcoeff=" << coeff;
      get_out() << " primalaggrcoeff=" << data.primal_aggregate_coeff;
      testmp.scale(coeff);
      if (!data.primal_aggregate.empty()) {
        if (testmp.aggregate(data.primal_aggregate, data.primal_aggregate_coeff))
          if (cb_out()) {
            get_out() << "*** WARNING: PSCModel::make_model_aggregate(): aggregation failed when forming the test aggregate" << std::endl;
          }
      }
      Matrix tmpvec(data.get_dim(), 1, 0.);
      Real val = 0.;
      data.get_aggregate().get_minorant(val, tmpvec, 1, 1.);
      testmp.get_minorant(val, tmpvec, 1, -1., true);
      get_out() << " updatePSCM aggrdevnorm=" << norm2(tmpvec) << " valdev=" << std::fabs(val) << std::endl;
    }


    assert((data.primalvecs.coldim() > 0) || (data.bundlevecs.coldim() == 0));
    assert(data.primaleigs.dim() == data.primalvecs.coldim());
    assert(data.primalvecs.rowdim() == data.bundlevecs.rowdim());
    assert(data.primalvecs.rowdim() == data.topvecs.rowdim());
    assert(data.primalgrowth.dim() == data.primaleigs.dim());

    return retval;
  }


  // *****************************************************************************
  //                        provide_local_aggregate
  // *****************************************************************************


  int PSCModel::provide_local_aggregate(Integer /* y_id */, const Matrix& /* iny */) {
    if (data.local_aggregate.valid()) {
      return 0;
    }

    //--- form an aggregate 
    data.local_aggregate.clear();


    assert((data.get_function_task() != ObjectiveFunction) || (std::fabs(sum(data.primaleigs) + data.primal_aggregate_coeff - data.get_function_factor()) < 1e-6 * data.get_function_factor()));

    Real coeff = 0.;
    if (model_selection->get_minorant(data.local_aggregate, coeff, data.primalvecs, data.primaleigs, oracle, data.get_modification_id())) {
      if (cb_out())
        get_out() << "**** ERROR PSCModel::provide_local_aggregate: get_minorant failed for the aggregate" << std::endl;
      return 1;
    }
    data.local_aggregate.scale(coeff);
    if ((!data.primal_aggregate.empty()) && (data.primal_aggregate_coeff > eps_Real * data.get_function_factor())) {
      if (data.local_aggregate.aggregate(data.primal_aggregate, data.primal_aggregate_coeff))
        return 1;
    }

    return 0;
  }

  // *****************************************************************************
  //                         switch_to_local
  // *****************************************************************************

  int PSCModel::switch_to_local(const MinorantBundle& current_model,
    const Matrix& current_coeff) {
    assert(Integer(current_model.size()) == current_coeff.rowdim());
    data.primal_aggregate.clear();
    data.primal_aggregate_coeff = sum(current_coeff);
    if ((data.primal_aggregate_coeff < eps_Real * data.get_function_factor()) &&
      (current_model.size() > 0)) {
      data.primal_aggregate = current_model[0];
    } else {
      if (data.primal_aggregate.aggregate(current_model, current_coeff, 1. / data.get_function_factor())) {
        if (cb_out()) {
          get_out() << "*** WARNING: PSCModel::sumbundle_mode(): getting the bundle_aggregate failed" << std::endl;
        }
      }
    }
    data.bundle_aggregate = data.primal_aggregate;
    data.primaleigs.init(0, 1, 0.);   //signals that no solution information is available
    data.growthrate = 0.;
    data.primalgrowth.init(0, 1, 0.);
    data.dualgrowth.init(0, 1, 0.);
    return 0;
  }

  // *****************************************************************************
  //                        get_aggregate_and_coefficient
  // *****************************************************************************

  int PSCModel::get_aggregate_and_coefficient(MinorantPointer& sumbundle_aggr,
    Real& aggr_coeff) {
    //set new model if required
    if (!data.primal_aggregate.empty()) {
      sumbundle_aggr = data.primal_aggregate;
      aggr_coeff = data.primal_aggregate_coeff;
      model_selection->get_minorant(sumbundle_aggr, aggr_coeff, data.primalvecs, data.primaleigs, oracle, data.get_modification_id());
    } else if (provide_cand_minorant().valid()) {
      sumbundle_aggr = data.cand_minorant;
      aggr_coeff = data.get_function_factor();
    } else if (provide_center_minorant().valid()) {
      sumbundle_aggr = data.center_minorant;
      aggr_coeff = data.get_function_factor();
    } else
      return 1;

    return 0;
  }


  // *****************************************************************************
  //                                start_local_augmodel
  // *****************************************************************************


  int PSCModel::start_local_augmodel(QPModelDataPointer& qpblockp,
    Integer /* cand_id */,
    const Matrix& /* cand_y */,
    const Indexmatrix* needed_indexsubset) {
    if (qpblockp.set_model_data(block)) {
      delete block;
      block = qpblockp.generate_conemodel_data();
      block->set_cbout(this);
    }
    assert(block);
    QPPSCOracleData* oracle_info = &qp_oracle_data;
    if (qp_oracle_data.init(oracle, &data))
      oracle_info = 0;

    //use the local model
    Integer s2dim = (data.bundlevecs.coldim() * (data.bundlevecs.coldim() + 1)) / 2;
    Integer xdim = s2dim + (data.bundle_aggregate.empty() ? 0 : 1);
    MinorantBundle bundle;
    bundle.resize(unsigned(xdim));
    unsigned cnt = 0;
    if (!data.bundle_aggregate.empty()) {
      bundle[cnt++] = data.bundle_aggregate;
      // if (cb_out(2)){
      // 	get_out()<<" bundle["<<cnt-1<<"]=";
      // 	bundle[cnt-1].display(get_out());
      // }
    }

    Matrix tmpvec;
    Matrix tmpmat;
    if (oracle->svec_projection(tmpvec, tmpmat, data.bundlevecs, needed_indexsubset)) {
      if (cb_out())
        get_out() << "**** ERROR: PSCModel::start_augmodel(): oracle->svec_projection() failed" << std::endl;
      return 1;
    }
    if (needed_indexsubset == 0) {
      //assert(tmpmat.rowdim()==cand_y.rowdim());
      assert(tmpmat.coldim() == s2dim);
      for (Integer i = 0; i < s2dim; i++) {
        bundle[cnt++].init(new Minorant(true, tmpvec(i), tmpmat.rowdim(), tmpmat.get_store() + i * tmpmat.rowdim()), data.get_modification_id());
        // if (cb_out(2)){
        //   get_out()<<" bundle["<<cnt-1<<"]=";
        //   bundle[cnt-1].display(get_out());
        // }
      }
    } else {
      Integer dim = needed_indexsubset->dim();
      assert((tmpmat.rowdim() == dim) && (tmpmat.coldim() == s2dim));
      for (Integer i = 0; i < s2dim; i++) {
        bundle[cnt++].init(new Minorant(true, tmpvec(i), dim, tmpmat.get_store() + i * dim, needed_indexsubset->get_store()), data.get_modification_id());
        // if (cb_out(2)){
        //   get_out()<<" bundle["<<cnt<<"]=";
        //   bundle[cnt-1].display(get_out());
        // }
      }
    }
    assert(Integer(cnt) == xdim);

    Indexmatrix tmpind(0, 1, Integer(0));
    Matrix dummymat(0, 1, 0.);
    if (!data.bundle_aggregate.empty()) {
      block->init(MinorantPointer(), bundle, 1, tmpind, Indexmatrix(1, 1, data.bundlevecs.coldim()), dummymat, dummymat, data.get_function_factor(), data.get_function_task(), oracle_info);
    } else {
      block->init(MinorantPointer(), bundle, 0, tmpind, Indexmatrix(1, 1, data.bundlevecs.coldim()), dummymat, dummymat, data.get_function_factor(), data.get_function_task(), oracle_info);
    }

    qpblockp.set_model_data(block);
    assert(qpblockp.get_model_data_ptr() == block);

    if (cb_out(2)) {
      get_out() << " PSCxdim=" << xdim;
    }

    return 0;
  }


  // *****************************************************************************
  //                                make_local_aggregate
  // *****************************************************************************

  int PSCModel::make_local_aggregate(Real& new_function_factor,
    bool keep_penalty_fixed) {
    assert(new_function_factor == data.get_function_factor());

    bool penalty_changed = false;
    data.local_aggregate.clear();

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
      Matrix linact;
      if (block->get_nncx(tmpvec, &linact)) {
        if (cb_out()) {
          get_out() << "**** ERROR PSCModel::make_local_aggregate(): retrieving the linear part of the QP solution failed" << std::endl;
        }
        return 1;
      }
      if ((tmpvec.dim() == 0) != data.bundle_aggregate.empty()) {
        if (cb_out()) {
          get_out() << "**** ERROR PSCModel::make_local_aggregate(): the retrieved linear QP solution does not match the number of aggregates" << std::endl;
        }
        return 1;
      }
      if (!data.bundle_aggregate.empty()) {
        data.primal_aggregate = data.bundle_aggregate;
        data.primal_aggregate_coeff = tmpvec(0);
      } else {
        data.primal_aggregate.clear();
        data.primal_aggregate_coeff = 0.;
      }

      Real cutvalub = block->tracedual();
      data.primal_Ritzval = cutvalub;
      bool cautious = false;
      //bool cautious=((data.get_function_task()!=ObjectiveFunction)&&(cutvalub<1e-6*data.get_function_factor()));
      if (cb_out(2)) {
        get_out().precision(8);
        get_out() << " PSC cutvalub=" << cutvalub;
        if (cautious)
          get_out() << " cautious" << std::endl;
      }

      Matrix tmpmat;
      if (block->get_pscx(0, data.primaleigs, tmpmat, data.growthrate, data.primalgrowth, data.dualgrowth)) {
        if (cb_out()) {
          get_out() << "**** ERROR PSCModel::make_local_aggregate(): retrieving the semidefinite part of the QP solution failed" << std::endl;
        }
        return 1;
      }
      if (data.primaleigs.rowdim() != data.bundlevecs.coldim()) {
        if (cb_out()) {
          get_out() << "**** ERROR PSCModel::make_local_aggregate(): the retrieved semidefinite QP solution does not match the number of bundle vectors" << std::endl;
        }
        return 1;
      }
      genmult(data.bundlevecs, tmpmat, data.primalvecs);

      Real coeff = 0.;
      data.local_aggregate.clear();
      if (model_selection->get_minorant(data.local_aggregate, coeff, data.primalvecs, data.primaleigs, oracle, data.get_modification_id())) {
        if (cb_out()) {
          get_out() << "*** WARNING: PSCModel::make_local_aggregate(): get_minorant failed when forming the aggregate" << std::endl;
        }
      }
      data.local_aggregate.scale(coeff);
      if ((!data.primal_aggregate.empty()) && (data.primal_aggregate_coeff > eps_Real * data.get_function_factor())) {
        if (data.local_aggregate.aggregate(data.primal_aggregate, data.primal_aggregate_coeff))
          return 1;
      }

      // //check consistency of this construction
      // if (cb_out(2)){
      //   Real val1=0.;
      //   Matrix grad1(data.get_dim(),1,0.);
      //   data.local_aggregate.get_minorant(val1,grad1,0,1.);
      //   Real val2=0.;
      //   Matrix grad2(data.get_dim(),1,0.);
      //   block->add_modelx_aggregate(val2,grad2);
      //   get_out().precision(6);
      //   get_out()<<" aggrconsisteny: valdev="<<val1-val2<<" normdev="<<norm2(grad1-grad2)<<std::endl;
      //   if ((std::fabs(val1-val2)>1e-10*(std::fabs(val1)+1.))||(norm2(grad1-grad2)>1e-10*(norm2(grad1)+1.))) {
      // 	data.local_aggregate.display(get_out());
      // 	if (!data.primal_aggregate.empty())
      // 	  get_out()<<" aggr_coeff="<<data.primal_aggregate_coeff;
      // 	get_out()<<" X="<<tmpmat*Diag(data.primaleigs)*transpose(tmpmat);
      // 	get_out()<<"svec="<<transpose(svec(Symmatrix(tmpmat*Diag(data.primaleigs)*transpose(tmpmat))));
      // 	get_out()<<" locval="<<val1<<" locgrad="<<transpose(grad1);
      // 	get_out()<<" blockval="<<val2<<" blockgrad="<<transpose(grad2);
      // 	for (unsigned int i=0;i<block->get_bundle().size();i++){
      // 	  get_out()<<" blockbundle["<<i<<"]=";
      // 	  block->get_bundle()[i].display(get_out());
      // 	}
      //   }
      // }

    }
    return 0;
  }

  // *****************************************************************************
  //                                 lb_function
  // *****************************************************************************

  //returns a *quick* lower bound for the function value at y
  //(eg by a previous subgradient)

  Real PSCModel::lb_function(Integer y_id, const Matrix& y) {
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

  MinorantPointer& PSCModel::provide_center_minorant() {
    if (data.get_center_modification_id() != data.get_modification_id()) {
      if (cb_out())
        get_out() << "**** WARNING: PSCModel::provide_center_minorant(...): no center available" << std::endl;
      data.center_minorant.clear();
    } else if (!data.center_minorant.valid()) {
      Real coeff = 0.;
      data.center_minorant.clear();
      if (model_selection->get_minorant(data.center_minorant, coeff, data.center_Ritzvec, Matrix(1, 1, 1.), oracle, data.get_modification_id())) {
        if (cb_out())
          get_out() << "**** WARNING: PSCModel::provide_center_minorant(...): get_minorant failed" << std::endl;
        data.center_minorant.clear();
      }
    }

    return data.center_minorant;
  }

  // *****************************************************************************
  //                             provide_cand_minorant
  // *****************************************************************************

  MinorantPointer& PSCModel::provide_cand_minorant() {
    if (data.cand_ub_mid != data.get_modification_id()) {
      if (cb_out())
        get_out() << "**** WARNING: PSCModel::provide_cand_minorant(...): no cand available" << std::endl;
      data.cand_minorant.clear();
    } else if (!data.cand_minorant.valid()) {
      Real coeff = 0.;
      if (model_selection->get_minorant(data.cand_minorant, coeff, data.topvecs.col(0), Matrix(1, 1, 1.), oracle, data.get_modification_id())) {
        if (cb_out())
          get_out() << "**** WARNING: PSCModel::provide_cand_minorant(...): get_minorant failed" << std::endl;
        data.cand_minorant.clear();
      }
    }

    return data.cand_minorant;
  }

  // *****************************************************************************
  //                            local_adjust_multiplier
  // *****************************************************************************


  int PSCModel::local_adjust_multiplier(Real& new_function_factor) {
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

  std::ostream& PSCModel::output_bundle_data(std::ostream& out) const {
    out << "\n%(begin PSCModel)\n";

    out << " n_aggr=" << data.bundle_aggregate.valid() << ";";
    out << "\n n_sym=" << data.bundlevecs.rowdim() << ";";
    Integer svecdim = (data.bundlevecs.coldim() * (data.bundlevecs.coldim() + 1)) / 2;
    Integer xdim = svecdim + (data.bundle_aggregate.empty() ? 0 : 1);
    out << "\n xdim=" << xdim << ";";
    Integer ydim = data.get_dim();
    out << "\n ydim=" << ydim << ";";


    Matrix bigvec(svecdim, 1, 0.);
    Matrix bigmat(ydim, svecdim, 0.);
    oracle->svec_projection(bigvec, bigmat, data.bundlevecs);

    out << "\n costs=[";
    if (data.bundle_aggregate.valid()) {
      out.width(18);
      out << " " << data.bundle_aggregate.offset();
    }
    for (Integer i = 0; i < svecdim; i++) {
      out.width(18);
      out << " " << (bigvec)(i);
    }
    out << "]';";

    out << "\n subgmat=[";
    for (Integer i = 0; i < ydim; i++) {
      if (data.bundle_aggregate.valid()) {
        out.width(18);
        out << " " << data.bundle_aggregate.coeff(i);
      }
      for (Integer j = 0; j < svecdim; j++) {
        out.width(18);
        out << " " << (bigmat)(i, j);
      }
      if (i < ydim - 1)
        out << "\n";
    }
    out << "];";

    out << "\n%(end PSCModel)\n";
    return out;
  }

  // *****************************************************************************
  //                         set_bundle_parameters      
  // *****************************************************************************

  int PSCModel::set_bundle_parameters(const BundleParameters& bp) {
    const PSCModelParametersObject* mpo = dynamic_cast<const PSCModelParametersObject*>(&bp);
    int retval = 0;
    if (mpo) {
      delete model_selection;
      model_selection = dynamic_cast<PSCModelParametersObject*>(mpo->clone_BundleParameters());
      assert(model_selection);
    } else {
      assert(model_selection);
      retval = model_selection->init(bp);
    }
    return retval;
  }



}
