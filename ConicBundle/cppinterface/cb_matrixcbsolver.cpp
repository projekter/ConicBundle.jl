dll void cb_matrixcbsolver_destroy(MatrixCBSolver* self) {
  delete self;
}

dll MatrixCBSolver* cb_matrixcbsolver_new(int print_level = 0) {
  return new MatrixCBSolver(&std::cout, print_level);
}

dll void cb_matrixcbsolver_clear(MatrixCBSolver* self) {
  self->clear();
}

dll void cb_matrixcbsolver_set_defaults(MatrixCBSolver* self) {
  self->set_defaults();
}

dll int cb_matrixcbsolver_init_problem(MatrixCBSolver* self, int dim, const Matrix* lbounds = 0, const Matrix* ubounds = 0, const Matrix* startval = 0, const Matrix* costs = 0, Real offset = 0.) {
  return self->init_problem(dim, lbounds, ubounds, startval, costs, offset);
}

dll int cb_matrixcbsolver_add_function(MatrixCBSolver* self, FunctionObject* function_, Real fun_factor = 1., int fun_task = ObjectiveFunction, AffineFunctionTransformation* aft = 0, int argument_list_may_change_dynamically = 0) {
  return self->add_function(*function_, fun_factor, (FunctionTask)fun_task, aft, (bool)argument_list_may_change_dynamically);
}

dll int cb_matrixcbsolver_set_lower_bound(MatrixCBSolver* self, int i, double lb) {
  return self->set_lower_bound(i, lb);
}

dll int cb_matrixcbsolver_set_upper_bound(MatrixCBSolver* self, int i, double ub) {
  return self->set_upper_bound(i, ub);
}

dll int cb_matrixcbsolver_append_variables(MatrixCBSolver* self, int n_append, const Matrix* lbounds = 0, const Matrix* ubounds = 0, const Sparsemat* constraint_columns = 0, const Matrix* startval = 0, const Matrix* costs = 0) {
  return self->append_variables(n_append, lbounds, ubounds, constraint_columns, startval, costs, (FunObjModMap*)0);
}

dll int cb_matrixcbsolver_delete_variables(MatrixCBSolver* self, const Indexmatrix* delete_indices, Indexmatrix* map_to_old) {
  return self->delete_variables(*delete_indices, *map_to_old, (FunObjModMap*)0);
}

dll int cb_matrixcbsolver_reassign_variables(MatrixCBSolver* self, const Indexmatrix* assign_new_from_old) {
  return self->reassign_variables(*assign_new_from_old, (FunObjModMap*)0);
}

dll int cb_matrixcbsolver_append_constraints(MatrixCBSolver* self, Integer append_n_rows, const Sparsemat* append_rows = 0, const Matrix* append_rhslb = 0, const Matrix* append_rhsub = 0) {
  return self->append_constraints(append_n_rows, append_rows, append_rhslb, append_rhsub);
}

dll int cb_matrixcbsolver_solve(MatrixCBSolver* self, int maxsteps = 0, int stop_at_descent_steps = 0) {
  return self->solve(maxsteps, (bool)stop_at_descent_steps);
}

dll int cb_matrixcbsolver_termination_code(const MatrixCBSolver* self) {
  return self->termination_code();
}

dll void cb_matrixcbsolver_print_termination_code(const MatrixCBSolver* self) {
  self->print_termination_code(std::cout);
}

dll double cb_matrixcbsolver_get_objval(const MatrixCBSolver* self) {
  return self->get_objval();
}

dll int cb_matrixcbsolver_get_center(const MatrixCBSolver* self, Matrix* center) {
  return self->get_center(*center);
}

dll double cb_matrixcbsolver_get_sgnorm(const MatrixCBSolver* self) {
  return self->get_sgnorm();
}

dll int cb_matrixcbsolver_get_subgradient(const MatrixCBSolver* self, Matrix* subgradient) {
  return self->get_subgradient(*subgradient);
}

dll double cb_matrixcbsolver_get_cutval(const MatrixCBSolver* self) {
  return self->get_cutval();
}

dll double cb_matrixcbsolver_get_candidate_value(const MatrixCBSolver* self) {
  return self->get_candidate_value();
}

dll int cb_matrixcbsolver_get_candidate(const MatrixCBSolver* self, Matrix* center) {
  return self->get_candidate(*center);
}

dll int cb_matrixcbsolver_set_term_relprec(MatrixCBSolver* self, const double term_relprec) {
  return self->set_term_relprec(term_relprec);
}

dll int cb_matrixcbsolver_set_new_center_point(MatrixCBSolver* self, const Matrix* center_point) {
  return self->set_new_center_point(*center_point);
}

dll int cb_matrixcbsolver_get_function_status(const MatrixCBSolver* self, const FunctionObject* function_) {
  return self->get_function_status(*function_);
}

dll int cb_matrixcbsolver_get_approximate_slacks(const MatrixCBSolver* self, Matrix* param0) {
  return self->get_approximate_slacks(*param0);
}

dll const PrimalData* cb_matrixcbsolver_get_approximate_primal(const MatrixCBSolver* self, const FunctionObject* function_) {
  return self->get_approximate_primal(*function_);
}

dll const PrimalData* cb_matrixcbsolver_get_center_primal(const MatrixCBSolver* self, const FunctionObject* function_) {
  return self->get_center_primal(*function_);
}

dll const PrimalData* cb_matrixcbsolver_get_candidate_primal(const MatrixCBSolver* self, const FunctionObject* function_) {
  return self->get_candidate_primal(*function_);
}

dll int cb_matrixcbsolver_set_sumbundle(MatrixCBSolver* self, int use_sumbundle, int n_local_models = -1, const BundleParameters* bundle_parameters = 0, int strategy = 1) {
  return self->set_sumbundle((bool)use_sumbundle, n_local_models, bundle_parameters, strategy);
}

dll int cb_matrixcbsolver_set_max_modelsize(MatrixCBSolver* self, int max_modelsize, const FunctionObject* function_ = 0) {
  return self->set_max_modelsize(max_modelsize, function_);
}

dll int cb_matrixcbsolver_set_max_bundlesize(MatrixCBSolver* self, int max_bundlesize, const FunctionObject* function_ = 0) {
  return self->set_max_bundlesize(max_bundlesize, function_);
}

dll int cb_matrixcbsolver_set_bundle_parameters(MatrixCBSolver* self, const BundleParameters* params, const FunctionObject* function_ = 0) {
  return self->set_bundle_parameters(*params, function_);
}

dll const BundleParameters* cb_matrixcbsolver_get_bundle_parameters(const MatrixCBSolver* self, const FunctionObject* function_ = 0) {
  return self->get_bundle_parameters(function_);
}

dll int cb_matrixcbsolver_set_sumbundle_parameters(MatrixCBSolver* self, const SumBundleParametersObject* params, const FunctionObject* function_ = 0) {
  return self->set_sumbundle_parameters(*params, function_);
}

dll const BundleData* cb_matrixcbsolver_get_bundle_data(const MatrixCBSolver* self, const FunctionObject* function_ = 0) {
  return self->get_bundle_data(function_);
}

dll int cb_matrixcbsolver_reinit_function_model(MatrixCBSolver* self, const FunctionObject* function_ = 0) {
  return self->reinit_function_model(function_);
}

dll int cb_matrixcbsolver_clear_aggregates(MatrixCBSolver* self, const FunctionObject* function_ = 0) {
  return self->clear_aggregates(function_);
}

dll int cb_matrixcbsolver_call_primal_extender(MatrixCBSolver* self, const FunctionObject* function_, PrimalExtender* primal_extender) {
  return self->call_primal_extender(*function_, *primal_extender);
}

dll double cb_matrixcbsolver_get_last_weight(const MatrixCBSolver* self) {
  return self->get_last_weight();
}

dll double cb_matrixcbsolver_get_next_weight(const MatrixCBSolver* self) {
  return self->get_next_weight();
}

dll int cb_matrixcbsolver_set_next_weight(MatrixCBSolver* self, const double weight) {
  return self->set_next_weight(weight);
}

dll int cb_matrixcbsolver_set_min_weight(MatrixCBSolver* self, const double min_weight) {
  return self->set_min_weight(min_weight);
}

dll int cb_matrixcbsolver_set_max_weight(MatrixCBSolver* self, const double max_weight) {
  return self->set_max_weight(max_weight);
}

dll int cb_matrixcbsolver_set_weight_update(MatrixCBSolver* self, BundleWeight* bw) {
  return self->set_weight_update(bw);
}

dll int cb_matrixcbsolver_adjust_multiplier(MatrixCBSolver* self) {
  return self->adjust_multiplier();
}

dll int cb_matrixcbsolver_set_variable_metric(MatrixCBSolver* self, int do_variable_metric) {
  return self->set_variable_metric(do_variable_metric);
}

dll int cb_matrixcbsolver_set_prox(MatrixCBSolver* self, BundleProxObject* proxp) {
  return self->set_prox(proxp);
}

dll void cb_matrixcbsolver_set_active_bounds_fixing(MatrixCBSolver* self, int allow_fixing) {
  self->set_active_bounds_fixing((bool)allow_fixing);
}

dll void cb_matrixcbsolver_clear_fail_counts(MatrixCBSolver* self) {
  self->clear_fail_counts();
}

dll void cb_matrixcbsolver_set_eval_limit(MatrixCBSolver* self, Integer eval_limit) {
  self->set_eval_limit(eval_limit);
}

dll void cb_matrixcbsolver_set_inner_update_limit(MatrixCBSolver* self, Integer update_limit) {
  self->set_inner_update_limit(update_limit);
}

dll void cb_matrixcbsolver_set_time_limit(MatrixCBSolver* self, Integer time_limit) {
  self->set_time_limit(time_limit);
}

dll int cb_matrixcbsolver_set_qp_solver(MatrixCBSolver* self, QPSolverParametersObject* qpparams, QPSolverObject* newqpsolver = 0) {
  return self->set_qp_solver(qpparams, newqpsolver);
}

dll int cb_matrixcbsolver_get_dim(const MatrixCBSolver* self) {
  return self->get_dim();
}

dll int cb_matrixcbsolver_get_n_functions(const MatrixCBSolver* self) {
  return self->get_n_functions();
}

dll int cb_matrixcbsolver_get_n_oracle_calls(const MatrixCBSolver* self) {
  return self->get_n_oracle_calls();
}

dll int cb_matrixcbsolver_get_n_descent_steps(const MatrixCBSolver* self) {
  return self->get_n_descent_steps();
}

dll int cb_matrixcbsolver_get_n_inner_iterations(const MatrixCBSolver* self) {
  return self->get_n_inner_iterations();
}

dll int cb_matrixcbsolver_get_n_inner_updates(const MatrixCBSolver* self) {
  return self->get_n_inner_updates();
}

dll int cb_matrixcbsolver_get_descent_step(const MatrixCBSolver* self) {
  return self->get_descent_step();
}

dll int cb_matrixcbsolver_get_null_step(const MatrixCBSolver* self) {
  return self->get_null_step();
}

dll int cb_matrixcbsolver_get_costs(const MatrixCBSolver* self, Matrix* costs) {
  return self->get_costs(*costs);
}

dll const Matrix* cb_matrixcbsolver_get_lbounds(const MatrixCBSolver* self) {
  return self->get_lbounds();
}

dll const Matrix* cb_matrixcbsolver_get_ubounds(const MatrixCBSolver* self) {
  return self->get_ubounds();
}

dll const Indexmatrix* cb_matrixcbsolver_get_fixed_active_bounds(const MatrixCBSolver* self) {
  return self->get_fixed_active_bounds();
}

dll BundleProxObject* cb_matrixcbsolver_get_prox(const MatrixCBSolver* self) {
  return self->get_prox();
}

dll void cb_matrixcbsolver_set_out(MatrixCBSolver* self, int print_level = 1) {
  self->set_out(&std::cout, print_level);
}

dll void cb_matrixcbsolver_print_line_summary(const MatrixCBSolver* self) {
  self->print_line_summary(std::cout);
}

dll void cb_matrixcbsolver_print_statistics(const MatrixCBSolver* self) {
  self->print_statistics(std::cout);
}

dll const BundleSolver* cb_matrixcbsolver_get_solver(const MatrixCBSolver* self) {
  return self->get_solver();
}

