dll void cb_unconstrainedgroundset_destroy(UnconstrainedGroundset* self) {
  delete self;
}

dll void cb_unconstrainedgroundset_clear(UnconstrainedGroundset* self, Integer indim = 0, Integer in_groundset_id = 0) {
  self->clear(indim, in_groundset_id);
}

dll UnconstrainedGroundset* cb_unconstrainedgroundset_new(Integer indim = 0, const Matrix* start_val = 0, const Matrix* costs = 0, const Real offset = 0., Integer in_groundset_id = 0) {
  return new UnconstrainedGroundset(indim, start_val, costs, offset, in_groundset_id);
}

dll Integer cb_unconstrainedgroundset_get_groundset_id(const UnconstrainedGroundset* self) {
  return self->get_groundset_id();
}

dll void cb_unconstrainedgroundset_set_groundset_id(UnconstrainedGroundset* self, Integer gsid) {
  self->set_groundset_id(gsid);
}

dll Integer cb_unconstrainedgroundset_get_dim(const UnconstrainedGroundset* self) {
  return self->get_dim();
}

dll int cb_unconstrainedgroundset_constrained(const UnconstrainedGroundset* self) {
  return self->constrained();
}

dll int cb_unconstrainedgroundset_is_feasible(UnconstrainedGroundset* self, Integer* in_groundset_id, const Matrix* y, Real relprec = 1e-10) {
  return self->is_feasible(*in_groundset_id, *y, relprec);
}

dll int cb_unconstrainedgroundset_ensure_feasibility(UnconstrainedGroundset* self, Integer* in_groundset_id, Matrix* y, int* ychanged, BundleProxObject* Hp = 0, Real relprec = 1e-10) {
  ychanged = 0;
  return self->ensure_feasibility(*in_groundset_id, *y, *(bool*)ychanged, Hp, relprec);
}

dll QPSolverObject* cb_unconstrainedgroundset_get_qp_solver(UnconstrainedGroundset* self, int* solves_model_without_gs, BundleProxObject* Hp) {
  solves_model_without_gs = 0;
  return self->get_qp_solver(*(bool*)solves_model_without_gs, Hp);
}

dll int cb_unconstrainedgroundset_set_qp_solver_parameters(UnconstrainedGroundset* self, QPSolverParametersObject* param0) {
  return self->set_qp_solver_parameters(param0);
}

dll const Matrix* cb_unconstrainedgroundset_get_starting_point(const UnconstrainedGroundset* self) {
  return &self->get_starting_point();
}

dll int cb_unconstrainedgroundset_set_starting_point(UnconstrainedGroundset* self, const Matrix* vec) {
  return self->set_starting_point(*vec);
}

dll int cb_unconstrainedgroundset_candidate(UnconstrainedGroundset* self, Integer* gs_id, Matrix* newy, Real* cand_gs_val, Real* linval, Real* augval_lb, Real* augval_ub, Real* subgnorm2, const Matrix* center_y, Real center_value, const MinorantPointer* model_minorant, BundleProxObject* Hp, MinorantPointer* delta_groundset_minorant = 0, Indexmatrix* delta_index = 0, Real relprec = 1e-2) {
  return self->candidate(*gs_id, *newy, *cand_gs_val, *linval, *augval_lb, *augval_ub, *subgnorm2, *center_y, center_value, *model_minorant, Hp, delta_groundset_minorant, delta_index, relprec);
}

dll const MinorantPointer* cb_unconstrainedgroundset_get_gs_aggregate(const UnconstrainedGroundset* self) {
  return &self->get_gs_aggregate();
}

dll const MinorantPointer* cb_unconstrainedgroundset_get_gs_minorant(const UnconstrainedGroundset* self) {
  return &self->get_gs_minorant();
}

dll const Indexmatrix* cb_unconstrainedgroundset_get_yfixed(const UnconstrainedGroundset* self) {
  return self->get_yfixed();
}

dll Indexmatrix* cb_unconstrainedgroundset_set_yfixed(UnconstrainedGroundset* self) {
  return self->set_yfixed();
}

dll int cb_unconstrainedgroundset_get_use_yfixing(const UnconstrainedGroundset* self) {
  return self->get_use_yfixing();
}

dll void cb_unconstrainedgroundset_set_use_yfixing(UnconstrainedGroundset* self, int uyf) {
  self->set_use_yfixing((bool)uyf);
}

dll GroundsetModification* cb_unconstrainedgroundset_start_modification(UnconstrainedGroundset* self) {
  return self->start_modification();
}

dll int cb_unconstrainedgroundset_apply_modification(UnconstrainedGroundset* self, const GroundsetModification* mdf) {
  return self->apply_modification(*mdf);
}

dll int cb_unconstrainedgroundset_mfile_data(const UnconstrainedGroundset* self) {
  return self->mfile_data(std::cout);
}

dll void cb_unconstrainedgroundset_set_cbout(UnconstrainedGroundset* self, int incr = -1) {
  self->set_cbout(0, incr);
}

