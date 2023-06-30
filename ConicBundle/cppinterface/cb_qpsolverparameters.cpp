dll void cb_qpsolverparameters_destroy(QPSolverParameters* self) {
  delete self;
}

dll QPSolverParameters* cb_qpsolverparameters_new(int incr = -1) {
  return new QPSolverParameters(0, incr);
}

dll Real cb_qpsolverparameters_qpget_min_objective_relprec(const QPSolverParameters* self) {
  return self->QPget_min_objective_relprec();
}

dll Real cb_qpsolverparameters_qpget_objective_gap_eps(const QPSolverParameters* self) {
  return self->QPget_objective_gap_eps();
}

dll Real cb_qpsolverparameters_qpget_primal_infeasibility_eps(const QPSolverParameters* self) {
  return self->QPget_primal_infeasibility_eps();
}

dll Real cb_qpsolverparameters_qpget_dual_infeasibility_eps(const QPSolverParameters* self) {
  return self->QPget_dual_infeasibility_eps();
}

dll Real cb_qpsolverparameters_qpget_lower_bound(const QPSolverParameters* self) {
  return self->QPget_lower_bound();
}

dll Real cb_qpsolverparameters_qpget_upper_bound(const QPSolverParameters* self) {
  return self->QPget_upper_bound();
}

dll Real cb_qpsolverparameters_qpget_lower_bound_gap_eps(const QPSolverParameters* self) {
  return self->QPget_lower_bound_gap_eps();
}

dll Real cb_qpsolverparameters_qpget_upper_bound_gap_eps(const QPSolverParameters* self) {
  return self->QPget_upper_bound_gap_eps();
}

dll Integer cb_qpsolverparameters_qpget_maxiter(const QPSolverParameters* self) {
  return self->QPget_maxiter();
}

dll QPKKTSolverObject* cb_qpsolverparameters_qpget_kktsolver(QPSolverParameters* self) {
  return self->QPget_KKTsolver();
}

dll int cb_qpsolverparameters_qpget_use_predictor_corrector(const QPSolverParameters* self) {
  return self->QPget_use_predictor_corrector();
}

dll int cb_qpsolverparameters_qpget_use_neighborhood(const QPSolverParameters* self) {
  return self->QPget_use_neighborhood();
}

dll Real cb_qpsolverparameters_qpget_nbh_ub(const QPSolverParameters* self) {
  return self->QPget_nbh_ub();
}

dll Real cb_qpsolverparameters_qpget_nbh_lb(const QPSolverParameters* self) {
  return self->QPget_nbh_lb();
}

dll int cb_qpsolverparameters_qpget_use_socqp(const QPSolverParameters* self) {
  return self->QPget_use_socqp();
}

dll int cb_qpsolverparameters_qpset_min_objective_relprec(QPSolverParameters* self, Real eps) {
  return self->QPset_min_objective_relprec(eps);
}

dll int cb_qpsolverparameters_qpset_objective_gap_eps(QPSolverParameters* self, Real eps) {
  return self->QPset_objective_gap_eps(eps);
}

dll int cb_qpsolverparameters_qpset_primal_infeasibility_eps(QPSolverParameters* self, Real eps) {
  return self->QPset_primal_infeasibility_eps(eps);
}

dll int cb_qpsolverparameters_qpset_dual_infeasibility_eps(QPSolverParameters* self, Real eps) {
  return self->QPset_dual_infeasibility_eps(eps);
}

dll int cb_qpsolverparameters_qpset_lower_and_upper_bounds(QPSolverParameters* self, Real lb, Real ub) {
  return self->QPset_lower_and_upper_bounds(lb, ub);
}

dll int cb_qpsolverparameters_qpset_lower_bound_gap_eps(QPSolverParameters* self, Real eps) {
  return self->QPset_lower_bound_gap_eps(eps);
}

dll int cb_qpsolverparameters_qpset_upper_bound_gap_eps(QPSolverParameters* self, Real eps) {
  return self->QPset_upper_bound_gap_eps(eps);
}

dll int cb_qpsolverparameters_qpset_maxiter(QPSolverParameters* self, Integer mi) {
  return self->QPset_maxiter(mi);
}

dll int cb_qpsolverparameters_qpset_kktsolver(QPSolverParameters* self, QPKKTSolverObject* in_KKTsolver) {
  return self->QPset_KKTsolver(in_KKTsolver);
}

dll int cb_qpsolverparameters_qpset_use_predictor_corrector(QPSolverParameters* self, int upc) {
  return self->QPset_use_predictor_corrector((bool)upc);
}

dll int cb_qpsolverparameters_qpset_use_neighborhood(QPSolverParameters* self, int nbh) {
  return self->QPset_use_neighborhood((bool)nbh);
}

dll int cb_qpsolverparameters_qpset_nbh_bounds(QPSolverParameters* self, Real nbhlb, Real nbhub) {
  return self->QPset_nbh_bounds(nbhlb, nbhub);
}

dll int cb_qpsolverparameters_qpset_use_socqp(QPSolverParameters* self, int s) {
  return self->QPset_use_socqp((bool)s);
}

dll int cb_qpsolverparameters_qpset_allow_uqpsolver(QPSolverParameters* self, int allow) {
  return self->QPset_allow_UQPSolver((bool)allow);
}

dll int cb_qpsolverparameters_qpallow_uqpsolver(QPSolverParameters* self) {
  return self->QPallow_UQPSolver();
}

