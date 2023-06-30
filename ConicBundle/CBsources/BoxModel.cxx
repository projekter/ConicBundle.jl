/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/BoxModel.cxx
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


#include <fstream>
#include "mymath.hxx"
#include "BoxModel.hxx"
#include "BoxModelParameters.hxx"
#include "BundleIdProx.hxx"

using namespace CH_Tools;
using namespace CH_Matrix_Classes;

namespace ConicBundle {

  // *****************************************************************************
  //                            BoxModel
  // *****************************************************************************

  void BoxModel::clear() {
    ConeModel::clear();
    data.clear();
    delete model_selection;
    model_selection = new BoxModelParameters;

    ret_code = 0;
    nr_eval = 0;
    delete block;
    block = 0;
  }

  // *****************************************************************************
  //                            BoxModel
  // *****************************************************************************

  BoxModel::BoxModel(BoxOracle* bo,
    Real fun_factor, FunctionTask fun_task, CBout* cb, int cbinc) :
    CBout(cb, cbinc), ConeModel(cb), oracle(bo), data(fun_factor, fun_task), model_selection(0), block(0) {
    assert(bo);
    if (min(bo->get_upper_bounds() - bo->get_lower_bounds()) < -eps_Real * (1. + std::fabs(max(bo->get_upper_bounds())))) {
      if (cb_out()) {
        get_out() << "**** WARNING BoxModel::BoxModel(): box oracle bounds do not satisfy lower bounds <= upper bounds" << std::endl;
      }
    }
    clear();
  }

  // *****************************************************************************
  //                            BoxModel
  // *****************************************************************************

  BoxModel::~BoxModel() {
    delete model_selection;
    delete block;
  }

  // *****************************************************************************
  //                            BoxModel
  // *****************************************************************************

  int BoxModel::adapt_bundle_to_parameters(Integer cand_id, const Matrix& cand_y) {
    data.outside_bundle = false;
    assert(data.nncbundle_coeff.rowdim() == Integer(data.nncbundle.size()));
    Integer kmax = 2; //do not use model_selection->n_model_size;
    assert(kmax >= 2);
    if (Integer(data.nncbundle.size()) <= kmax) //nothing to do
      return 0;

    MinorantBundle bundle = data.nncbundle;
    Matrix coeff = data.nncbundle_coeff;
    data.nncbundle_coeff.init(kmax, 1, 0.);
    data.nncbundle.resize(unsigned(kmax));

    // first find the best subgradient and ensure it is contained
    Matrix objval(Integer(bundle.size()), 1, 0.);
    Integer bestind = -1;
    Real bestval = min_Real;
    for (Integer i = 0; i < objval.rowdim(); i++) {
      Real d = bundle[unsigned(i)].evaluate(cand_id, cand_y);
      objval(i) = d;
      if (d > bestval) {
        bestval = d;
        bestind = i;
      }
    }
    Integer k = 0;
    if (bestind >= 0) {
      data.nncbundle[unsigned(k)] = bundle[unsigned(bestind)];
      data.nncbundle_coeff(k) = coeff(bestind);
      k++;
      objval(bestind) = max_Real; //mark as treated
    }

    // next aggregate all with possitive coefficients
    Indexmatrix sind;
    sortindex(coeff, sind, false);
    MinorantBundle aggrbundle;
    Matrix aggrcoeff;
    Integer i = 0;
    while ((i < coeff.rowdim()) && (coeff(i) > eps_Real * data.get_function_factor())) {
      Integer bi = sind(i++);
      if (bi == bestind)
        continue;
      objval(bi) = max_Real; //mark as treated
      const MinorantPointer& mp = bundle[unsigned(bi)];
      Integer idind = 0;
      while ((idind < k) && (!mp.equals(data.nncbundle[unsigned(idind)], eps_Real)))
        idind++;
      if (idind < k) {
        data.nncbundle_coeff(idind) += coeff(bi);
      } else {
        if (k < kmax) {
          data.nncbundle[unsigned(k)] = mp;
          data.nncbundle_coeff(k) = coeff(bi);
          k++;
        } else {
          Real d = coeff(bi);
          if (d > 1e-10 * data.get_function_factor()) {
            aggrcoeff.concat_below(d);
            aggrbundle.push_back(mp);
          }
        }
      }
    }
    if (aggrcoeff.rowdim() > 0) {
      assert(Integer(k) == kmax);
      MinorantPointer& aggr = data.nncbundle[unsigned(kmax - 1)];
      Real d = data.nncbundle_coeff(kmax - 1) + sum(aggrcoeff);
      aggr.scale(1. / d);
      aggr.aggregate(aggrbundle, aggrcoeff, 1. / d);
      data.nncbundle_coeff(kmax - 1) = d;
      aggrbundle.clear();
    }

    //if not all have been treated and there is still room, sort by objective value
    if ((i < coeff.rowdim()) && (k < kmax)) {
      sortindex(objval, sind, false);
      while ((i < coeff.rowdim()) && (k < kmax)) {
        Integer bi = sind(i++);
        if (bi == bestind)
          continue;
        const MinorantPointer& mp = bundle[unsigned(bi)];
        Integer idind = 0;
        while ((idind < k) && (!mp.equals(data.nncbundle[unsigned(idind)], eps_Real)))
          idind++;
        if (idind < k) //identical minorant 
          continue;
        data.nncbundle[unsigned(k)] = mp;
        data.nncbundle_coeff(k) = coeff(bi);
        k++;
      }
    }

    if (k < kmax) {
      data.nncbundle_coeff.reduce_length(k);
      data.nncbundle.resize(unsigned(k));
    }

    data.primal_nnctapia = (data.nncbundle_coeff > 1e-8 * data.get_function_factor());

    return 0;
  }


  // *****************************************************************************
  //                            evaluate_oracle
  // *****************************************************************************


  int BoxModel::evaluate_oracle(Integer y_id,
    const Matrix& y,
    Real nullstep_bound,
    Real relprec) {
    assert(relprec > 0.);

    CH_Tools::Microseconds start_eval = clock.time();


    int err = 0;
    Integer dim = y.rowdim();
    if (dim != oracle->get_lower_bounds().rowdim()) {
      if (cb_out())
        get_out() << "**** ERROR: BoxModel::eval_function(): dimension of lower bound vector lb = " << oracle->get_lower_bounds().rowdim() << " does not match to function space dimension = " << y.rowdim() << std::endl;
      err++;
    }
    if (dim != oracle->get_upper_bounds().rowdim()) {
      if (cb_out())
        get_out() << "**** ERROR: BoxModel::eval_function(): dimension of upper bound vector ub = " << oracle->get_upper_bounds().rowdim() << " does not match to function space dimension = " << y.rowdim() << std::endl;
      err++;
    }
    if (err)
      return err;

    //delete old info

    data.cand_minorant.clear();
    data.cand_ub = nullstep_bound;
    data.cand_id = y_id;
    data.cand_ub_mid = -1;//signals no valid computation yet

    //---- evaluate  	
    std::vector<Minorant*> minorants;
    //PrimalExtender* pep=0;			

    nr_eval++;
    preeval_time += clock.time() - start_eval;
    start_eval = clock.time();

    //ret_code = oracle->evaluate(y,relprec,data.cand_ub, minorants,pep);

    const Matrix& olb = oracle->get_lower_bounds();
    const Matrix& oub = oracle->get_upper_bounds();
    Matrix& cand_boxvec = data.cand_boxvec;
    Matrix& coord_switching = data.coord_switching;

    bool trace_switches = true;
    if (coord_switching.rowdim() != dim) {
      coord_switching.init(dim, 1, 0.);
      trace_switches = false;
    } else {
      coord_switching *= .9;
    }

    if (cand_boxvec.dim() != dim) {
      cand_boxvec.init(olb);
      trace_switches = false;
    }

    Real val = 0;
    if (trace_switches) {
      for (Integer i = 0; i < dim; i++) {
        if (y(i) < 0) {
          Real d = olb(i);
          val += d * y(i);
          if (cand_boxvec(i) != d) {
            cand_boxvec(i) = d;
            coord_switching(i) += 1.;
          }
        } else {
          Real d = oub(i);
          val += d * y(i);
          if (cand_boxvec(i) != d) {
            cand_boxvec(i) = d;
            coord_switching(i) += 1.;
          }
        }
      }
    } else {
      for (Integer i = 0; i < dim; i++) {
        if (y(i) < 0) {
          Real d = olb(i);
          val += d * y(i);
          cand_boxvec(i) = d;
        } else {
          Real d = oub(i);
          val += d * y(i);
          cand_boxvec(i) = d;
        }
      }
    }

    data.cand_ub_mid = data.get_modification_id();
    data.cand_ub = val;
    data.cand_relprec = relprec;
    data.cand_minorant = MinorantPointer(new MatrixMinorant(0, cand_boxvec, 0, true), data.get_modification_id());

    eval_time += clock.time() - start_eval;
    start_eval = clock.time();


    //--- maybe required later, as a reminder: first apply here the modifications to the old data
    // if (pep){
    //   if (call_primal_extender(*pep))
    //     if (cb_out()) 
    // 	get_out()<<"**** WARNING: BoxModel::eval_function(): data.call_primal_extender failed"<<std::endl;
    //   delete pep;
    // }

    if (data.cand_ub >= CB_plus_infinity) {
      if (cb_out())
        get_out() << "**** WARNING: BoxModel::eval_function(): function value cand_ub>=CB_plus_infinity" << std::endl;
    }

    posteval_time += clock.time() - start_eval;
    return (err > 0) ? err : ret_code;
  }

  // *****************************************************************************
  //                              evaluate_local_model
  // *****************************************************************************

  //evaluate the current cutting model in $y$ 
  //if evaluated by an iterative method that provides upper and lower bounds
  //it may stop when the lower bound (lb) is above the nullstep_bound
  //evalutation may also stop, if (ub-lb)<relprec*(|ub|+1) is satisfied

  int BoxModel::evaluate_local_model(Real& lb,
    Integer y_id,
    const Matrix& y,
    Real /* relprec */) {
    assert(data.nncbundle.size() + data.boxbundle.size() > 0);

    if (data.boxbundle.size() > 0) {
      const Matrix& olb = oracle->get_lower_bounds();
      const Matrix& oub = oracle->get_upper_bounds();
      lb = ip(y, data.boxbundle_complvalues);
      for (Integer i = 1; i < data.boxbundle_coords.dim(); i++) {
        Integer ind = data.boxbundle_coords(i);
        if (y(ind) < 0.)
          lb += y(ind) * olb(ind);
        else
          lb += y(ind) * oub(ind);
      }
    } else {
      lb = min_Real;
    }
    for (unsigned i = 0; i < data.nncbundle.size(); i++)
      lb = max(lb, data.nncbundle[i].evaluate(y_id, y));

    return 0;
  }


  // *****************************************************************************
  //                                update_model
  // *****************************************************************************


  int BoxModel::update_local_model(ModelUpdate model_update,
    Integer center_id,
    const Matrix& center_y,
    Integer y_id,
    const Matrix& y,
    Real model_maxviol,
    BundleProxObject& H) {
    assert(unsigned(data.nncbundle_coeff.dim()) == data.nncbundle.size());

    if ((model_update != new_subgradient) && (data.cand_ub_mid != data.get_modification_id())) {
      if (cb_out())
        get_out() << "**** ERROR BoxModel::update_model: candidate required but its values are not up to date" << std::endl;
      return 1;
    }

    if (model_update == descent_step) {
      data.do_step(y_id);
    }

    //initialize the bundle if necessary
    if (data.boxbundle.size() + data.nncbundle.size() == 0) {
      data.outside_bundle = false;
      data.boxbundle_coords.init(0, 1, Integer(0));
      data.boxbundle_complvalues.init(0, 1, 0.);
      if (data.cand_ub_mid != data.get_modification_id()) {
        //input_id must be center_id
        if ((y_id != data.get_center_id()) || (data.get_center_modification_id() != data.get_modification_id())) {
          if (cb_out())
            get_out() << "**** ERROR BoxModel::update_model: point_id should match recomputed center but it does not or center is not up to date" << std::endl;
          return 1;
        }
        data.nncbundle.push_back(data.center_minorant);
        data.nncbundle_coeff.init(1, 1, data.get_function_factor());
        data.primal_nnctapia.init(1, 1, 1.);
      } else {
        data.nncbundle.push_back(data.cand_minorant);
        data.nncbundle_coeff.init(1, 1, data.get_function_factor());
        data.primal_nnctapia.init(1, 1, 1.);
      }


      return 0;
    }

    if ((!data.get_local_model_active()) || (model_update == new_subgradient)) {
      //it is a root or a child in sumbundle, all is taken care of by sumbundle
      //but if desired here the local bundle could be updated along, maybe with more

      return 0;
    }

    //=========== the model selection routine
    data.outside_bundle = false;
    assert(data.coord_switching.dim() == y.dim());

    assert(model_selection);
    int retval = model_selection->select_model(data.boxbundle,
      data.boxbundle_coeff,
      data.primal_boxtapia,
      data.boxbundle_coords,
      data.boxbundle_complvalues,
      data.nncbundle,
      data.nncbundle_coeff,
      data.primal_nnctapia,
      data.coord_switching,
      data.get_old_minorants(),
      data.cand_minorant,
      data.cand_boxvec,
      data.aggr_boxvec,
      data.aggr_scaleval,
      oracle,
      data.get_modification_id(),
      data.get_function_task(),
      data.get_function_factor(),
      model_update,
      center_id,
      center_y,
      y_id,
      y,
      model_maxviol,
      H);
    if ((retval) && (cb_out())) {
      get_out() << "**** WARNING BoxModel::update_model: model_selection->select_model() returned " << retval << std::endl;
    }

    if (cb_out(2)) {
      get_out() << " BoxModel::update_model: boxbundle_size=" << data.boxbundle.size() << " linbundle_size=" << data.nncbundle.size() << std::endl;
    }

    return retval;

  }


  // *****************************************************************************
  //                        provide_local_aggregate
  // *****************************************************************************


  int BoxModel::provide_local_aggregate(Integer /* y_id */, const Matrix& iny) {
    if (data.local_aggregate.valid()) {
      return 0;
    }

    //--- form an aggregate
    data.local_aggregate.clear();

    if (data.aggr_boxvec.rowdim() == iny.rowdim()) {
      data.local_aggregate = MinorantPointer(new MatrixMinorant(0, data.aggr_boxvec, 0, true), data.get_modification_id(), data.aggr_scaleval);
    } else if (data.cand_boxvec.rowdim() == iny.rowdim()) {
      data.local_aggregate = MinorantPointer(new MatrixMinorant(0, data.cand_boxvec, 0, true), data.get_modification_id(), data.get_function_factor());
    } else {
      assert(data.center_boxvec.rowdim() == iny.rowdim());
      data.local_aggregate = MinorantPointer(new MatrixMinorant(0, data.center_boxvec, 0, true), data.get_modification_id(), data.get_function_factor());
    }

    return 0;
  }



  // *****************************************************************************
  //                         switch_to_local
  // *****************************************************************************

  int BoxModel::switch_to_local(const MinorantBundle& current_model,
    const Matrix& current_coeff) {
    assert(Integer(current_model.size()) == current_coeff.rowdim());
    data.nncbundle = current_model;
    data.nncbundle_coeff = current_coeff;
    if (data.aggr_boxvec.rowdim() == 0) {
      data.aggr_scaleval = sum(data.nncbundle_coeff);
      data.aggr_boxvec.init(oracle->get_lower_bounds().rowdim(), 1, 0.);
      if (data.aggr_scaleval > eps_Real * data.get_function_factor())
        genmult(data.nncbundle, data.nncbundle_coeff, data.aggr_boxvec, 1. / data.aggr_scaleval);
    }
    data.outside_bundle = true;
    return 0;
  }

  // *****************************************************************************
  //                        get_aggregate_and_coefficient
  // *****************************************************************************

  int BoxModel::get_aggregate_and_coefficient(MinorantPointer& sumbundle_aggr,
    Real& aggr_coeff) {
    //set new model if required
    if (data.aggr_boxvec.rowdim() > 0) {
      sumbundle_aggr.init(new MatrixMinorant(0, data.aggr_boxvec, 0, true), data.get_modification_id());
      aggr_coeff = data.aggr_scaleval;
    } else if (provide_cand_minorant().valid()) {
      sumbundle_aggr = data.cand_minorant;
      aggr_coeff = data.get_function_factor();
    } else
      return 1;
    return 0;
  }


  // *****************************************************************************
  //                                start_local_augmodel
  // *****************************************************************************

  //return a pointer to the variables/constraints generating the cutting model
  //returns 0 on success, 1 on failure

  int BoxModel::start_local_augmodel(QPModelDataPointer& qpblockp,
    Integer cand_id,
    const Matrix& cand_y,
    const Indexmatrix* /* needed_indexsubset */) {

    if (qpblockp.set_model_data(block)) {
      delete block;
      block = qpblockp.generate_conemodel_data();
      block->set_cbout(this);
    }
    assert(block);

    if (data.outside_bundle) {
      adapt_bundle_to_parameters(cand_id, cand_y);
    }
    MinorantBundle bundle;
    bundle.reserve(data.nncbundle.size() + data.boxbundle.size());
    for (unsigned i = 0; i < data.nncbundle.size(); i++) {
      bundle.push_back(data.nncbundle[i]);
    }
    for (unsigned i = 0; i < data.boxbundle.size(); i++) {
      bundle.push_back(data.boxbundle[i]);
    }
    Integer xdim = Integer(data.nncbundle.size());
    bool scale_box = (xdim > 0) || (Integer(data.boxbundle.size()) == data.boxbundle_coords.dim() + 1);
    assert((scale_box) || ((data.get_function_task() == ObjectiveFunction) && (Integer(data.boxbundle.size()) == oracle->get_lower_bounds().dim())));
    Indexmatrix tmpind(0, 0, Integer(0));
    int status = block->init(MinorantPointer(), bundle, xdim, tmpind, tmpind, oracle->get_lower_bounds()(data.boxbundle_coords), oracle->get_upper_bounds()(data.boxbundle_coords), data.get_function_factor(), data.get_function_task(), 0, scale_box);
    if (status) {
      if (cb_out()) {
        get_out() << "*** ERROR: BoxModel::start_augmodel(): initializing the QP block failed and returned" << status << std::endl;
      }
      return status;
    }

    qpblockp.set_model_data(block);
    assert(qpblockp.get_model_data_ptr() == block);

    if (cb_out(2)) {
      get_out() << " BMsz=" << bundle.size();
    }

    return status;
  }

  // *****************************************************************************
  //                        make_local_aggregate
  // *****************************************************************************

  int BoxModel::make_local_aggregate(Real& new_function_factor,
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

      Real scaleval = 0.;
      data.aggr_boxvec.init(oracle->get_lower_bounds().dim(), 1, 0.);

      //get coefficients and activity for linbundle
      if (data.nncbundle.size() > 0) {
        if (block->get_nncx(data.nncbundle_coeff, &data.primal_nnctapia)) {
          if (cb_out()) {
            get_out() << "**** ERROR BoxModel::make_model_aggregate(): retrieving the QP solution failed for the linear part" << std::endl;
          }
          return 1;
        }

        if (data.nncbundle_coeff.dim() != Integer(data.nncbundle.size())) {
          if (cb_out()) {
            get_out() << "**** ERROR BoxModel::make_model_aggregate(): the retrieved QP solution does not match the bundle size" << std::endl;
          }
          return 1;
        }

        //form aggregate contribution of linbundle into data.aggr_boxvec;
        scaleval = sum(data.nncbundle_coeff);
        if (scaleval < 1e-10 * data.get_function_factor()) {
          //use the old point scaled to 0
          Real dummy;
          data.nncbundle[0].get_minorant(dummy, data.aggr_boxvec, 0, 1.);
        } else {
          Real dummy = 0.;
          for (Integer i = 0; i < data.nncbundle_coeff.dim(); i++) {
            if (data.nncbundle[unsigned(i)].get_minorant(dummy, data.aggr_boxvec, 0, data.nncbundle_coeff(i) / scaleval, true)) {
              if (cb_out()) {
                get_out() << "**** ERROR BoxModel::make_model_aggregate(): forming the aggregate failed" << std::endl;
              }
              return 1;
            }
            assert(dummy == 0.);
          }
        }

      }

      //get coefficients and activity for boxbundle
      if (data.boxbundle.size() > 0) {
        if (block->get_boxx(data.boxbundle_coeff, &data.primal_boxtapia)) {
          if (cb_out()) {
            get_out() << "**** ERROR BoxModel::make_model_aggregate(): retrieving the QP solution failed for the box part" << std::endl;
          }
          return 1;
        }

        if (data.boxbundle_coeff.dim() != Integer(data.boxbundle.size())) {
          if (cb_out()) {
            get_out() << "**** ERROR BoxModel::make_model_aggregate(): the retrieved QP solution does not match the bundle size" << std::endl;
          }
          return 1;
        }

        Matrix boxaggr;
        Real boxscal;
        if (data.boxbundle_coeff.dim() == data.boxbundle_coords.dim()) {
          assert(data.boxbundle_coeff.dim() == oracle->get_upper_bounds().dim());
          assert(data.get_function_task() == ObjectiveFunction);
          assert(data.nncbundle.size() == 0);
          boxscal = data.get_function_factor();
          boxaggr.newsize(data.boxbundle_coeff.dim(), 1); chk_set_init(boxaggr, 1);
        } else {
          assert(data.boxbundle_coeff.dim() == data.boxbundle_coords.dim() + 1);
          boxscal = data.boxbundle_coeff(data.boxbundle_coeff.dim() - 1);
          if (boxscal < 1e-10 * data.get_function_factor()) {
            //scaled point is zero, avoid numerical difficulties in feasibility
            if (data.cand_boxvec.dim() == oracle->get_lower_bounds().dim())
              boxaggr = data.cand_boxvec;
            else {
              assert(data.center_boxvec.dim() == oracle->get_lower_bounds().dim());
              boxaggr = data.center_boxvec;
            }
          } else {
            boxaggr.init(data.boxbundle_complvalues);
            for (Integer i = 0; i < data.boxbundle_coords.dim(); i++) {
              boxaggr(data.boxbundle_coords(i)) = data.boxbundle_coeff(i) / boxscal;
            }
          }
        }

        scaleval += boxscal;
        Real alpha = boxscal / scaleval;
        xbpeya(data.aggr_boxvec, boxaggr, alpha, 1. - alpha);
      }

      //check if the constructed aggregate is box feasible
      for (Integer i = 0; i < data.aggr_boxvec.dim(); i++) {
        if (data.aggr_boxvec(i) < oracle->get_lower_bounds()(i)) {
          if ((data.aggr_boxvec(i) - oracle->get_lower_bounds()(i)) < -eps_Real * (1. + std::fabs(oracle->get_lower_bounds()(i))))
            if (cb_out())
              get_out() << "**** WARNING BoxModel::make_model_aggregate(): aggregate(" << i << ")=" << data.aggr_boxvec(i) << " < lower bound = " << oracle->get_lower_bounds()(i) << ", setting to lower bound" << std::endl;
          data.aggr_boxvec(i) = oracle->get_lower_bounds()(i);
        }
        if (data.aggr_boxvec(i) > oracle->get_upper_bounds()(i)) {
          if ((oracle->get_upper_bounds()(i) - data.aggr_boxvec(i)) < -eps_Real * (1. + std::fabs(oracle->get_upper_bounds()(i))))
            if (cb_out())
              get_out() << "**** WARNING BoxModel::make_model_aggregate(): aggregate(" << i << ")=" << data.aggr_boxvec(i) << " > upper bound = " << oracle->get_upper_bounds()(i) << ", setting to upper bound" << std::endl;
          data.aggr_boxvec(i) = oracle->get_upper_bounds()(i);
        }
      }

      scaleval = min(max(0., scaleval), data.get_function_factor());

      data.aggr_scaleval = scaleval;
      data.local_aggregate = MinorantPointer(new MatrixMinorant(0., data.aggr_boxvec, 0, true), data.get_modification_id());
      data.local_aggregate.scale(scaleval);
    }

    return 0;
  }

  // *****************************************************************************
  //                                 lb_function
  // *****************************************************************************

  //returns a *quick* lower bound for the function value at y
  //(eg by a previous subgradient)

  Real BoxModel::lb_function(Integer y_id, const Matrix& y) {

    Real lb = CB_minus_infinity;

    if (data.get_aggregate().valid()) {
      lb = data.get_aggregate().evaluate(y_id, y);
    } else if (data.nncbundle.size() > 0) {
      lb = data.nncbundle[0].evaluate(y_id, y);
      for (unsigned i = 1; i < data.nncbundle.size(); i++) {
        lb = max(lb, data.nncbundle[i].evaluate(y_id, y));
      }
      lb = data.get_function_factor() * lb;
    } else if ((data.cand_ub_mid == data.get_modification_id()) &&
      (data.cand_minorant.valid())) {
      lb = data.cand_minorant.evaluate(y_id, y);
      lb = data.get_function_factor() * lb;
    }

    if (data.get_function_task() != ObjectiveFunction)
      lb = max(0., lb);

    return lb;
  }

  // *****************************************************************************
  //                         provide_cand_minorant
  // *****************************************************************************


  MinorantPointer& BoxModel::provide_cand_minorant() {
    if (data.cand_ub_mid != data.get_modification_id())
      data.cand_minorant.clear();
    return data.cand_minorant;
  }


  // *****************************************************************************
  //                            provide_center_minorant
  // *****************************************************************************


  MinorantPointer& BoxModel::provide_center_minorant() {
    if (data.get_center_modification_id() != data.get_modification_id())
      data.center_minorant.clear();
    return data.center_minorant;
  }

  // *****************************************************************************
  //                            local_adjust_multiplier
  // *****************************************************************************


  int BoxModel::local_adjust_multiplier(Real& new_function_factor) {
    assert(new_function_factor == data.get_function_factor());

    if (data.get_function_task() != AdaptivePenaltyFunction) {
      return 0;
    }

    assert(block);
    Real oldval = data.get_function_factor();
    Real traceval = block->evaluate_trace();
    new_function_factor = min(data.get_function_factor(), max(1., 1.5 * traceval));
    if (new_function_factor < oldval) {
      block->adjust_trace(new_function_factor);
      data.aggr_boxvec.init(0, 0, 0.);
    }

    return 0;
  }


  // *****************************************************************************
  //                         set_bundle_parameters      
  // *****************************************************************************

  int BoxModel::set_bundle_parameters(const BundleParameters& bp) {
    const BoxModelParametersObject* mpo = dynamic_cast<const BoxModelParametersObject*>(&bp);
    int retval = 0;
    if (mpo) {
      delete model_selection;
      model_selection = dynamic_cast<BoxModelParametersObject*>(mpo->clone_BundleParameters());
      assert(model_selection);
    } else {
      assert(model_selection);
      retval = model_selection->init(bp);
    }
    return retval;
  }



}
