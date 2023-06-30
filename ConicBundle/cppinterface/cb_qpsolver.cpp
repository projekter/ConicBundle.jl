dll void cb_qpsolver_destroy(QPSolver* self) {
  delete self;
}

dll void cb_qpsolver_qpclear(QPSolver* self) {
  self->QPclear();
}

dll QPSolver* cb_qpsolver_new(int cbinc = -1) {
  return new QPSolver(0, cbinc);
}

dll int cb_qpsolver_qpset_parameters(QPSolver* self, QPSolverParametersObject* params) {
  return self->QPset_parameters(params);
}

dll int cb_qpsolver_qpsupports_yfixing(QPSolver* self) {
  return self->QPsupports_yfixing();
}

dll int cb_qpsolver_qpsupports_updates(QPSolver* self) {
  return self->QPsupports_updates();
}

dll Real cb_qpsolver_qpget_lower_bound(QPSolver* self) {
  return self->QPget_lower_bound();
}

dll int cb_qpsolver_qpsolve(QPSolver* self, const Matrix* center_y, Real lower_bound, Real upper_bound, Real relprec, QPSolverProxObject* Hp, const MinorantPointer* gs_aggr, Indexmatrix* yfixed) {
  return self->QPsolve(*center_y, lower_bound, upper_bound, relprec, Hp, *gs_aggr, yfixed);
}

dll int cb_qpsolver_qpupdate(QPSolver* self, const Matrix* center_y, Real lower_bound, Real upper_bound, Real relprec, QPSolverProxObject* Hp, const MinorantPointer* gs_aggr, Indexmatrix* yfixed, const MinorantPointer* delta_gs_aggr, const Indexmatrix* delta_index) {
  return self->QPupdate(*center_y, lower_bound, upper_bound, relprec, Hp, *gs_aggr, yfixed, *delta_gs_aggr, *delta_index);
}

dll int cb_qpsolver_qpresolve(QPSolver* self, Real lower_bound, Real upper_bound, Real relprec) {
  return self->QPresolve(lower_bound, upper_bound, relprec);
}

dll int cb_qpsolver_qpget_solution(QPSolver* self, Real* augval_lb, Real* augval_ub, Matrix* new_point, Real* gsaggr_offset, Matrix* gsaggr_gradient) {
  return self->QPget_solution(*augval_lb, *augval_ub, *new_point, *gsaggr_offset, *gsaggr_gradient);
}

dll void cb_qpsolver_qpprint_statistics(QPSolver* self, int param1 = 0) {
  self->QPprint_statistics(std::cout, param1);
}

dll GroundsetModification* cb_qpsolver_qpstart_modification(QPSolver* self) {
  return self->QPstart_modification();
}

dll int cb_qpsolver_qpapply_modification(QPSolver* self, const GroundsetModification* mdf) {
  return self->QPapply_modification(*mdf);
}

dll int cb_qpsolver_qpis_feasible(QPSolver* self, const Matrix* y, Real relprec = 1e-10) {
  return self->QPis_feasible(*y, relprec);
}

dll int cb_qpsolver_qpensure_feasibility(QPSolver* self, Matrix* y, int* ychanged, QPSolverProxObject* inHp, Real relprec = 1e-10) {
  ychanged = 0;
  return self->QPensure_feasibility(*y, *(bool*)ychanged, inHp, relprec);
}

dll int cb_qpsolver_solve(QPSolver* self, BundleProxObject* Hp, const Matrix* c, Real gamma, Real lowerbound, Real upperbound, Real relprec, Real skip_factor) {
  return self->solve(Hp, *c, gamma, lowerbound, upperbound, relprec, skip_factor);
}

dll int cb_qpsolver_qpprefer_uqpsolver(const QPSolver* self, QPSolverProxObject* param0) {
  return self->QPprefer_UQPSolver(param0);
}

dll int cb_qpsolver_qpconstrained(const QPSolver* self) {
  return self->QPconstrained();
}

dll Integer cb_qpsolver_rowdim(const QPSolver* self) {
  return self->rowdim();
}

dll const Matrix* cb_qpsolver_get_lby(const QPSolver* self) {
  return &self->get_lby();
}

dll const Matrix* cb_qpsolver_get_uby(const QPSolver* self) {
  return &self->get_uby();
}

dll const Indexmatrix* cb_qpsolver_get_lbindex(const QPSolver* self) {
  return &self->get_lbindex();
}

dll const Indexmatrix* cb_qpsolver_get_ubindex(const QPSolver* self) {
  return &self->get_ubindex();
}

dll const Sparsemat* cb_qpsolver_get_a(const QPSolver* self) {
  return &self->get_A();
}

dll const Matrix* cb_qpsolver_get_rhslb(const QPSolver* self) {
  return &self->get_rhslb();
}

dll const Matrix* cb_qpsolver_get_rhsub(const QPSolver* self) {
  return &self->get_rhsub();
}

dll const Indexmatrix* cb_qpsolver_get_rhslbind(const QPSolver* self) {
  return &self->get_rhslbind();
}

dll const Indexmatrix* cb_qpsolver_get_rhsubind(const QPSolver* self) {
  return &self->get_rhsubind();
}

dll int cb_qpsolver_mfile_data(const QPSolver* self) {
  return self->mfile_data(std::cout);
}

dll void cb_qpsolver_set_cbout(QPSolver* self, int incr = -1) {
  self->set_cbout(0, incr);
}

