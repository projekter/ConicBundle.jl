/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/CBSolver.cxx
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



#include "MatrixCBSolver.hxx"
#include <algorithm>
#include <map>

//------------------------------------------------------------


using namespace CH_Matrix_Classes;

namespace ConicBundle {


  PrimalData::~PrimalData() {
  }

  PrimalExtender::~PrimalExtender() {
  }

  MinorantExtender::~MinorantExtender() {
  }

  OracleModification::~OracleModification() {
  }

  FunctionObject::~FunctionObject() {
  }

  int FunctionOracle::apply_modification(
    const OracleModification& /* oracle_modification */,
    const double* /* new_center */,
    const double* /* old_center */,
    bool& /* discard_objective_in_center */,
    bool& /* discard_model */,
    bool& /* discard_aggregates */,
    MinorantExtender*& /* minorant_extender */
  ) {
    return 1;
  }

  BundleParameters::~BundleParameters() {
  }

  //------------------------------------------------------------
  // CBmethod implementation - mostly passing data to MatrixCBSolver
  //------------------------------------------------------------

  //--------------------
  CBSolver::CBSolver(std::ostream* out, int print_level) {
    solver = new MatrixCBSolver(out, print_level);
    assert(solver);
  }

  //--------------------
  CBSolver::~CBSolver() {
    assert(solver);
    delete solver;
  }

  //--------------------
  void CBSolver::clear() {
    assert(solver);
    solver->clear();
  }

  //--------------------
  void CBSolver::set_defaults() {
    assert(solver);
    solver->set_defaults();
  }

  //--------------------
  int CBSolver::init_problem(int dim,
    const DVector* lbounds,
    const DVector* ubounds) {
    assert(solver);
    Matrix lb;
    Matrix* lbp;
    if (lbounds == 0) lbp = 0;
    else {
      lb.init(*lbounds);
      lbp = &lb;
    }
    Matrix ub;
    Matrix* ubp;
    if (ubounds == 0) ubp = 0;
    else {
      ub.init(*ubounds);
      ubp = &ub;
    }
    return solver->init_problem(dim, lbp, ubp);
  }

  //--------------------
  int CBSolver::add_function(FunctionObject& function) {
    assert(solver);
    return solver->add_function(function);
  }


  //----------------------------------------
  // append new variables (always in last postions in this order)
  int CBSolver::append_variables(int n_append,
    const DVector* lbounds,
    const DVector* ubounds) {
    assert(solver);
    Matrix lb;
    Matrix* lbp;
    if (lbounds == 0) lbp = 0;
    else {
      lb.init(*lbounds);
      lbp = &lb;
    }
    Matrix ub;
    Matrix* ubp;
    if (ubounds == 0) ubp = 0;
    else {
      ub.init(*ubounds);
      ubp = &ub;
    }
    return solver->append_variables(n_append, lbp, ubp);
  }

  //----------------------------------------
  // delete variables
  int CBSolver::delete_variables(const IVector& del_indices,
    IVector& map_to_old) {
    assert(solver);
    Indexmatrix mapold;
    int ret_code = solver->delete_variables(Indexmatrix(del_indices), mapold);
    assign(map_to_old, mapold);
    return ret_code;
  }


  //----------------------------------------
  // reassign variables
  int CBSolver::reassign_variables(const IVector& aind) {
    assert(solver);
    return solver->reassign_variables(Indexmatrix(aind));
  }

  //----------------------------------------
  int CBSolver::set_lower_bound(int i, double lb) {
    assert(solver);
    return solver->set_lower_bound(i, lb);
  }

  //----------------------------------------
  int CBSolver::set_upper_bound(int i, double ub) {
    assert(solver);
    return solver->set_upper_bound(i, ub);
  }


  //--------------------
  int CBSolver::solve(int maxsteps, bool stop_at_descent) {
    assert(solver);
    return solver->solve(maxsteps, stop_at_descent);
  }

  //--------------------
  int CBSolver::termination_code() const {
    assert(solver);
    return solver->termination_code();
  }

  //--------------------
  std::ostream& CBSolver::print_termination_code(std::ostream& out) const {
    assert(solver);
    return solver->print_termination_code(out);
  }

  //--------------------
  double CBSolver::get_objval() const {
    assert(solver);
    return solver->get_objval();
  }

  //--------------------
  int CBSolver::get_center(DVector& center) const {
    assert(solver);
    Matrix cen;
    int ret_code = solver->get_center(cen);
    assign(center, cen);
    return ret_code;
  }

  //--------------------
  double CBSolver::get_candidate_value() const {
    assert(solver);
    return solver->get_candidate_value();
  }

  //--------------------
  int CBSolver::get_candidate(DVector& center) const {
    assert(solver);
    Matrix cen;
    int ret_code = solver->get_candidate(cen);
    assign(center, cen);
    return ret_code;
  }

  //--------------------
  int CBSolver::get_approximate_slacks(DVector& slacks)const {
    assert(solver);
    Matrix sl;
    int ret_code = solver->get_approximate_slacks(sl);
    assign(slacks, sl);
    return ret_code;
  }

  //--------------------
  const PrimalData* CBSolver::get_approximate_primal(const FunctionObject& function) const {
    assert(solver);
    return solver->get_approximate_primal(function);
  }


  //--------------------
  const PrimalData* CBSolver::get_center_primal(const FunctionObject& function) const {
    assert(solver);
    return solver->get_center_primal(function);
  }

  //--------------------
  const PrimalData* CBSolver::get_candidate_primal(const FunctionObject& function) const {
    assert(solver);
    return solver->get_candidate_primal(function);
  }

  //--------------------
  double CBSolver::get_sgnorm() const {
    assert(solver);
    return solver->get_sgnorm();
  }

  //--------------------
  int CBSolver::get_subgradient(DVector& subgradient) const {
    assert(solver);
    Matrix subg;
    int ret_code = solver->get_subgradient(subg);
    assign(subgradient, subg);
    return ret_code;
  }

  //--------------------
  int CBSolver::get_function_status(const FunctionObject& function) const {
    assert(solver);
    return solver->get_function_status(function);
  }

  //--------------------
  int CBSolver::set_max_modelsize(const FunctionObject& function,
    int ms) {
    assert(solver);
    return solver->set_max_modelsize(ms, &function);
  }

  //--------------------
  int CBSolver::set_max_bundlesize(const FunctionObject& function,
    int mb) {
    assert(solver);
    return solver->set_max_bundlesize(mb, &function);
  }

  //--------------------
  int CBSolver::set_bundle_parameters(const FunctionObject& function,
    const BundleParameters& bp) {
    assert(solver);
    return solver->set_bundle_parameters(bp, &function);
  }

  //--------------------
  const BundleParameters* CBSolver::get_bundle_parameters(const FunctionObject& function) const {
    assert(solver);
    return solver->get_bundle_parameters(&function);
  }

  //--------------------
  int CBSolver::reinit_function_model(const FunctionObject& function) {
    assert(solver);
    return solver->reinit_function_model(&function);
  }

  //--------------------
  int CBSolver::call_primal_extender(const FunctionObject& function, PrimalExtender& primal_extender) {
    assert(solver);
    return solver->call_primal_extender(function, primal_extender);
  }

  //--------------------
  int CBSolver::set_term_relprec(const double term_relprec) {
    assert(solver);
    solver->set_term_relprec(term_relprec);
    return 0;
  }

  //--------------------
  double CBSolver::get_last_weight() const {
    assert(solver);
    return solver->get_last_weight();
  }

  //--------------------
  int CBSolver::set_next_weight(const double weight) {
    assert(solver);
    return solver->set_next_weight(weight);
  }

  //--------------------
  int CBSolver::set_min_weight(const double weight) {
    assert(solver);
    return solver->set_min_weight(weight);
  }

  //--------------------
  int CBSolver::set_max_weight(const double weight) {
    assert(solver);
    return solver->set_max_weight(weight);
  }

  //--------------------
  int CBSolver::set_new_center_point(const DVector& center_point) {
    assert(solver);
    return solver->set_new_center_point(Matrix(center_point));
  }

  //--------------------
  int CBSolver::set_variable_metric(int do_scaling) {
    assert(solver);
    return solver->set_variable_metric(do_scaling);
  }

  //--------------------
  void CBSolver::set_active_bounds_fixing(bool allow_fixing) {
    assert(solver);
    solver->set_active_bounds_fixing(allow_fixing);
  }


  void CBSolver::clear_fail_counts(void) {
    assert(solver);
    return solver->clear_fail_counts();
  }

  void CBSolver::set_eval_limit(int eval_limit) {
    assert(solver);
    return solver->set_eval_limit(eval_limit);
  }

  void CBSolver::set_inner_update_limit(int update_limit) {
    assert(solver);
    return solver->set_inner_update_limit(update_limit);
  }

  int CBSolver::get_dim() {
    assert(solver);
    return solver->get_dim();
  }

  int CBSolver::get_n_functions() {
    assert(solver);
    return solver->get_n_functions();
  }

  int CBSolver::get_fixed_active_bounds(IVector& vind) const {
    assert(solver);
    const Indexmatrix* ind = solver->get_fixed_active_bounds();
    if (ind) {
      vind.resize(unsigned(ind->rowdim()));
      for (Integer i = 0; i < ind->rowdim(); i++) {
        vind[unsigned(i)] = (*ind)[i];
      }
    } else {
      vind.resize(unsigned(solver->get_dim()));
      for (unsigned i = 0; i < vind.size(); i++) {
        vind[i] = 0;
      }
    }
    return 0;
  }


  void CBSolver::set_out(std::ostream* o, int pril) {
    assert(solver);
    solver->set_out(o, pril);
  }

}

