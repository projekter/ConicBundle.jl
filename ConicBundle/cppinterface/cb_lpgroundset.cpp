dll void cb_lpgroundset_destroy(LPGroundset* self) {
  delete self;
}

dll void cb_lpgroundset_clear(LPGroundset* self, Integer indim = 0, Integer in_groundset_id = 0) {
  self->clear(indim, in_groundset_id);
}

dll LPGroundset* cb_lpgroundset_new() {
  return new LPGroundset(0);
}

dll LPGroundset* cb_lpgroundset_new2(Integer dim, const Matrix* lbyp = 0, const Matrix* ubyp = 0, const Sparsemat* Gp = 0, const Matrix* rhslbp = 0, const Matrix* rhsubp = 0, const Matrix* start_val = 0, const Matrix* costs = 0, const Real offset = 0., Integer in_groundset_id = 0) {
  return new LPGroundset(dim, lbyp, ubyp, Gp, rhslbp, rhsubp, start_val, costs, offset, in_groundset_id, 0);
}

dll Integer cb_lpgroundset_get_groundset_id(const LPGroundset* self) {
  return self->get_groundset_id();
}

dll void cb_lpgroundset_set_groundset_id(LPGroundset* self, Integer gsid) {
  self->set_groundset_id(gsid);
}

dll Integer cb_lpgroundset_get_dim(const LPGroundset* self) {
  return self->get_dim();
}

dll int cb_lpgroundset_set_qpsolver(LPGroundset* self, QPSolverParametersObject* qpparams, QPSolverObject* qpsolver = 0) {
  return self->set_qpsolver(qpparams, qpsolver);
}

dll int cb_lpgroundset_constrained(const LPGroundset* self) {
  return self->constrained();
}

dll const Matrix* cb_lpgroundset_get_lby(const LPGroundset* self) {
  return self->get_lby();
}

dll const Matrix* cb_lpgroundset_get_uby(const LPGroundset* self) {
  return self->get_uby();
}

dll int cb_lpgroundset_is_feasible(LPGroundset* self, Integer* in_groundset_id, const Matrix* y, Real relprec = 1e-10) {
  return self->is_feasible(*in_groundset_id, *y, relprec);
}

dll int cb_lpgroundset_ensure_feasibility(LPGroundset* self, Integer* in_groundset_id, Matrix* y, int* ychanged, BundleProxObject* Hp, Real relprec = 1e-10) {
  ychanged = 0;
  return self->ensure_feasibility(*in_groundset_id, *y, *(bool*)ychanged, Hp, relprec);
}

dll QPSolverObject* cb_lpgroundset_get_qp_solver(LPGroundset* self, int* solves_model_without_gs, BundleProxObject* Hp) {
  solves_model_without_gs = 0;
  return self->get_qp_solver(*(bool*)solves_model_without_gs, Hp);
}

dll int cb_lpgroundset_set_qp_solver_parameters(LPGroundset* self, QPSolverParametersObject* param0) {
  return self->set_qp_solver_parameters(param0);
}

dll const Matrix* cb_lpgroundset_get_starting_point(const LPGroundset* self) {
  return &self->get_starting_point();
}

dll int cb_lpgroundset_set_starting_point(LPGroundset* self, const Matrix* vec) {
  return self->set_starting_point(*vec);
}

dll int cb_lpgroundset_candidate(LPGroundset* self, Integer* gs_id, Matrix* newy, Real* cand_gs_val, Real* linval, Real* augval_lb, Real* augval_ub, Real* subgnorm2, const Matrix* center_y, Real center_value, const MinorantPointer* model_minorant, BundleProxObject* Hp, MinorantPointer* delta_groundset_minorant = 0, Indexmatrix* delta_index = 0, Real relprec = 1e-2) {
  return self->candidate(*gs_id, *newy, *cand_gs_val, *linval, *augval_lb, *augval_ub, *subgnorm2, *center_y, center_value, *model_minorant, Hp, delta_groundset_minorant, delta_index, relprec);
}

dll const MinorantPointer* cb_lpgroundset_get_gs_aggregate(const LPGroundset* self) {
  return &self->get_gs_aggregate();
}

dll const MinorantPointer* cb_lpgroundset_get_gs_minorant(const LPGroundset* self) {
  return &self->get_gs_minorant();
}

dll const Indexmatrix* cb_lpgroundset_get_yfixed(const LPGroundset* self) {
  return self->get_yfixed();
}

dll Indexmatrix* cb_lpgroundset_set_yfixed(LPGroundset* self) {
  return self->set_yfixed();
}

dll int cb_lpgroundset_get_use_yfixing(const LPGroundset* self) {
  return self->get_use_yfixing();
}

dll void cb_lpgroundset_set_use_yfixing(LPGroundset* self, int uyf) {
  self->set_use_yfixing((bool)uyf);
}

dll int cb_lpgroundset_set_variable_metric_selection(LPGroundset* self, VariableMetricSelection* vms = 0) {
  return self->set_variable_metric_selection(vms);
}

dll VariableMetricSelection* cb_lpgroundset_get_variable_metric_selection(const LPGroundset* self) {
  return self->get_variable_metric_selection();
}

dll int cb_lpgroundset_add_variable_metric(LPGroundset* self, VariableMetric* H, Integer y_id, const Matrix* y, int descent_step, Real weightu, Real model_maxviol, const Indexmatrix* indices = 0) {
  return self->add_variable_metric(*H, y_id, *y, (bool)descent_step, weightu, model_maxviol, indices);
}

dll GroundsetModification* cb_lpgroundset_start_modification(LPGroundset* self) {
  return self->start_modification();
}

dll int cb_lpgroundset_apply_modification(LPGroundset* self, const GroundsetModification* mdf) {
  return self->apply_modification(*mdf);
}

dll int cb_lpgroundset_mfile_data(const LPGroundset* self) {
  return self->mfile_data(std::cout);
}

dll void cb_lpgroundset_set_cbout(LPGroundset* self, int incr = -1) {
  self->set_cbout(0, incr);
}

