dll void cb_uqpsolver_destroy(UQPSolver* self) {
  delete self;
}

dll void cb_uqpsolver_clear(UQPSolver* self) {
  self->clear();
}

dll void cb_uqpsolver_set_defaults(UQPSolver* self) {
  self->set_defaults();
}

dll UQPSolver* cb_uqpsolver_new(int cbinc = -1) {
  return new UQPSolver(0, cbinc);
}

dll void cb_uqpsolver_init_size(UQPSolver* self, Integer maxdim) {
  self->init_size(maxdim);
}

dll const Symmatrix* cb_uqpsolver_get_q(const UQPSolver* self) {
  return &self->get_Q();
}

dll const Matrix* cb_uqpsolver_get_c(const UQPSolver* self) {
  return &self->get_c();
}

dll Real cb_uqpsolver_get_offset(const UQPSolver* self) {
  return self->get_offset();
}

dll void cb_uqpsolver_set_termbounds(UQPSolver* self, Real lb, Real ub) {
  self->set_termbounds(lb, ub);
}

dll void cb_uqpsolver_set_termeps(UQPSolver* self, Real te) {
  self->set_termeps(te);
}

dll void cb_uqpsolver_set_maxiter(UQPSolver* self, Integer mi) {
  self->set_maxiter(mi);
}

dll Integer cb_uqpsolver_get_iter(const UQPSolver* self) {
  return self->get_iter();
}

dll Integer cb_uqpsolver_get_status(const UQPSolver* self) {
  return self->get_status();
}

dll Real cb_uqpsolver_get_termeps(const UQPSolver* self) {
  return self->get_termeps();
}

dll Integer cb_uqpsolver_get_maxiter(const UQPSolver* self) {
  return self->get_maxiter();
}

dll Real cb_uqpsolver_get_primalval(const UQPSolver* self) {
  return self->get_primalval();
}

dll Real cb_uqpsolver_get_dualval(const UQPSolver* self) {
  return self->get_dualval();
}

dll const Matrix* cb_uqpsolver_get_x(UQPSolver* self) {
  return &self->get_x();
}

dll const Matrix* cb_uqpsolver_get_y(UQPSolver* self) {
  return &self->get_y();
}

dll int cb_uqpsolver_solve(UQPSolver* self, const Symmatrix* Q, const Matrix* c, Real offset) {
  return self->solve(*Q, *c, offset);
}

dll int cb_uqpsolver_resolve(UQPSolver* self) {
  return self->resolve();
}

dll int cb_uqpsolver_update(UQPSolver* self, const Symmatrix* dQ, const Matrix* dc, Real doffset) {
  return self->update(*dQ, *dc, doffset);
}

dll void cb_uqpsolver_print_statistics(const UQPSolver* self) {
  self->print_statistics(std::cout);
}

dll void cb_uqpsolver_save(const UQPSolver* self) {
  self->save(std::cout);
}

dll void cb_uqpsolver_qpclear(UQPSolver* self) {
  self->QPclear();
}

dll int cb_uqpsolver_qpset_parameters(UQPSolver* self, QPSolverParametersObject* param0) {
  return self->QPset_parameters(param0);
}

dll int cb_uqpsolver_apply_modification(UQPSolver* self, const GroundsetModification* param0) {
  return self->apply_modification(*param0);
}

dll int cb_uqpsolver_qpsupports_yfixing(UQPSolver* self) {
  return self->QPsupports_yfixing();
}

dll int cb_uqpsolver_qpsupports_updates(UQPSolver* self) {
  return self->QPsupports_updates();
}

dll Real cb_uqpsolver_qpget_lower_bound(UQPSolver* self) {
  return self->QPget_lower_bound();
}

dll int cb_uqpsolver_qpsolve(UQPSolver* self, const Matrix* center_y, Real lower_bound, Real upper_bound, Real relprec, QPSolverProxObject* Hp, const MinorantPointer* gs_aggr, Indexmatrix* yfixed) {
  return self->QPsolve(*center_y, lower_bound, upper_bound, relprec, Hp, *gs_aggr, yfixed);
}

dll int cb_uqpsolver_qpupdate(UQPSolver* self, const Matrix* center_y, Real lower_bound, Real upper_bound, Real relprec, QPSolverProxObject* Hp, const MinorantPointer* gs_aggr, Indexmatrix* yfixed, const MinorantPointer* delta_gs_subg, const Indexmatrix* delta_index) {
  return self->QPupdate(*center_y, lower_bound, upper_bound, relprec, Hp, *gs_aggr, yfixed, *delta_gs_subg, *delta_index);
}

dll int cb_uqpsolver_qpresolve(UQPSolver* self, Real lower_bound, Real upper_bound, Real relprec) {
  return self->QPresolve(lower_bound, upper_bound, relprec);
}

dll int cb_uqpsolver_qpget_solution(UQPSolver* self, Real* augval_lb, Real* augval_ub, Matrix* new_point, Real* gsaggr_offset, Matrix* gsaggr_gradient) {
  return self->QPget_solution(*augval_lb, *augval_ub, *new_point, *gsaggr_offset, *gsaggr_gradient);
}

dll int cb_uqpsolver_qpprefer_uqpsolver(const UQPSolver* self, QPSolverProxObject* param0) {
  return self->QPprefer_UQPSolver(param0);
}

dll int cb_uqpsolver_qpconstrained(const UQPSolver* self) {
  return self->QPconstrained();
}

dll GroundsetModification* cb_uqpsolver_qpstart_modification(UQPSolver* self) {
  return self->QPstart_modification();
}

dll int cb_uqpsolver_qpapply_modification(UQPSolver* self, const GroundsetModification* param0) {
  return self->QPapply_modification(*param0);
}

dll int cb_uqpsolver_qpis_feasible(UQPSolver* self, const Matrix* param0, Real param1 = 1e-10) {
  return self->QPis_feasible(*param0, param1);
}

dll int cb_uqpsolver_qpensure_feasibility(UQPSolver* self, Matrix* param0, int* ychanged, QPSolverProxObject* param2, Real param3 = 1e-10) {
  ychanged = 0;
  return self->QPensure_feasibility(*param0, *(bool*)ychanged, param2, param3);
}

dll void cb_uqpsolver_qpprint_statistics(UQPSolver* self, int param1 = 0) {
  self->QPprint_statistics(std::cout, param1);
}

