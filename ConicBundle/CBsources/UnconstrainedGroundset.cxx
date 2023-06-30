/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/UnconstrainedGroundset.cxx
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



#include <typeinfo>
#include "mymath.hxx"
#include "UnconstrainedGroundset.hxx"
#include "BundleIdProx.hxx"
#include "UQPSolver.hxx"

using namespace CH_Matrix_Classes;

namespace ConicBundle {

  // *****************************************************************************
  //                              clear()
  // *****************************************************************************

  void UnconstrainedGroundset::clear(Integer indim, Integer gs_id) {
    indim = (indim > 0) ? indim : 0;
    assert(gs_id >= 0);
    groundset_id = gs_id;
    dim = 0;
    starting_point.init(0, 1, 0.);
    gs_aggregate.init(new Minorant, groundset_id);
    use_yfixing = false;
    yfixed.init(0, 1, Integer(0));
    Hp = 0;
    c.init(0, 1, 0.);
    gamma = 0;
    delete qp_solver;
    UQPSolver* qps = new UQPSolver;
    qp_solver = qps;

    GroundsetModification mdf(0);
    mdf.add_append_vars(indim);
    groundset_id--;  //next operation increments the counter
    apply_modification(mdf);
  }

  // *****************************************************************************
  //                              UnconstrainedGroundset()
  // *****************************************************************************

  UnconstrainedGroundset::UnconstrainedGroundset(CH_Matrix_Classes::Integer indim,
    const CH_Matrix_Classes::Matrix* start_val,
    const CH_Matrix_Classes::Matrix* costs,
    const CH_Matrix_Classes::Real offset,
    CH_Matrix_Classes::Integer in_groundset_id) {
    qp_solver = 0;
    indim = (indim > 0) ? indim : 0;
    clear(0, in_groundset_id);
    GroundsetModification mdf(0);
    mdf.add_append_vars(indim, start_val, costs);
    mdf.add_offset(offset);
    groundset_id--; //next operation increments the counter
    apply_modification(mdf);
  }


  // *****************************************************************************
  //                              ~UnconstrainedGroundset()
  // *****************************************************************************

  UnconstrainedGroundset::~UnconstrainedGroundset(void) {
  }

  // *****************************************************************************
  //                           UnconstrainedGroundset::is_feasible()
  // *****************************************************************************

  bool UnconstrainedGroundset::is_feasible(Integer& gs_id, const Matrix& y, Real) {
    if (y.dim() != dim) {
      return false;
    }
    gs_id = groundset_id;
    return true;
  }

  // *****************************************************************************
  //                           UnconstrainedGroundset::ensure_feasibility()
  // *****************************************************************************

  int UnconstrainedGroundset::ensure_feasibility(Integer& gs_id,
    Matrix& y,
    bool& ychanged,
    BundleProxObject*,
    Real) {
    gs_id = groundset_id;
    if (y.dim() != dim) {
      y.init(dim, 1, 0.);
      ychanged = true;
    }
    return 0;
  }

  // *****************************************************************************
  //                           UnconstrainedGroundset::get_qp_solver()
  // *****************************************************************************

  QPSolverObject* UnconstrainedGroundset::get_qp_solver(bool& solves_model_without_gs,
    BundleProxObject* /* Hp */) {
    assert(qp_solver);
    solves_model_without_gs = false;
    return qp_solver;
  }

  // *****************************************************************************
  //                              UnconstrainedGroundset::candidate
  // *****************************************************************************

  int UnconstrainedGroundset::candidate(Integer& gs_id,
    Matrix& newy,
    Real& cand_gs_val,
    Real& linval,
    Real& augval_lb,
    Real& augval_ub,
    Real& normsubg2,
    const Matrix& center_y,
    Real /* center_value */,
    const MinorantPointer& model_subg,
    BundleProxObject* Hp,
    MinorantPointer* delta_groundset_aggregate,
    Indexmatrix* delta_index,
    Real /* relprec */) {
    assert(Hp);
    assert((delta_groundset_aggregate == 0) == (delta_index == 0));

    gs_id = groundset_id;

    Real dummy;
    newy.newsize(dim, 1); chk_set_init(newy, 1);
    model_subg.get_minorant(dummy, newy, 0, -1, false);
    gs_aggregate.get_minorant(dummy, newy, 0, -1., true);
    Hp->apply_Hinv(newy);
    normsubg2 = -model_subg.ip(newy) - gs_aggregate.ip(newy);
    newy += center_y;
    cand_gs_val = gs_aggregate.evaluate(-1, newy);
    linval = model_subg.evaluate(-1, newy) + cand_gs_val;
    //augval should not change
    assert(augval_lb < linval + normsubg2 / 2. + 1e-8 * (1. + std::fabs(linval + normsubg2 / 2.)));
    //for consistency with the subgradients we adapt it anyway
    augval_lb = augval_ub = linval + normsubg2 / 2.;

    if (delta_groundset_aggregate) {
      delta_groundset_aggregate->init(new Minorant, 0);
      delta_index->init(0, 0, Integer(0));
    }

    return 0;
  }


  // *****************************************************************************
  //                              UnconstrainedGroundset::apply_modification
  // *****************************************************************************

  int UnconstrainedGroundset::apply_modification(const GroundsetModification& mdf) {
    if (dim != mdf.old_vardim()) {
      if (cb_out())
        get_out() << "**** ERROR: UnconstrainedGroundset::apply_modification(.): there are " << dim << " variables, but modification assumes " << mdf.old_vardim() << " variables" << std::endl;
      return 1;
    }
    if (mdf.no_modification() && (mdf.get_add_offset() == 0.))
      return 0;
    int err = 0;
    dim = mdf.new_vardim();
    yfixed.init(dim, 1, 0);
    groundset_id++;
    if (gs_aggregate.apply_modification(mdf, groundset_id, 0, true)) {
      if (cb_out())
        get_out() << "**** ERROR: UnconstrainedGroundset::apply_modification(.): modification of the groundset aggregate failed" << std::endl;
      err++;
    }
    if (mdf.apply_to_vars(starting_point)) {
      if (cb_out())
        get_out() << "**** ERROR: UnconstrainedGroundset::apply_modification(.): modification of the starting point failed" << std::endl;
      err++;
    }

    // if (qp_solver->apply_modification(mdf)){
    //   if (cb_out())
    //     get_out()<<"**** ERROR: UnconstrainedGroundset::apply_modification(.): modification of qp_solver failed"<<std::endl;
    //   err++;
    // }

    int dummy = -1;
    if (!is_feasible(dummy, starting_point)) {
      if (cb_out())
        get_out() << "**** WARNING: UnconstrainedGroundset::apply_modification(.): starting point is not feasible (but this is allowed)" << std::endl;
      err++;
    }

    return err;
  }


  // *************************************************************************
  //                       UnconstrainedGroundset::mfile_data
  // *************************************************************************

  int UnconstrainedGroundset::mfile_data(std::ostream& out) const {
    out << "clear gs_subg gs_sugb_offset yfixed G rhs lby uby\n";
    out << "gs_subg=[";
    for (Integer i = 0; i < dim; i++) {
      out.precision(16);
      out.width(18);
      out << gs_aggregate.coeff(i);
      if (i < dim - 1)
        out << "\n";
    }
    out << "];\n";
    out << "gs_subg_offset=" << gs_aggregate.offset() << ";\n";
    out << "yfixed=[";
    for (Integer i = 0; i < yfixed.dim(); i++) {
      out.precision(16);
      out.width(18);
      out << yfixed(i);
      if (i < yfixed.dim() - 1)
        out << "\n";
    }
    out << "];\n";
    return 0;
  }






}

