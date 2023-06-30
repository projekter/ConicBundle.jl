dll void cb_boxdata_destroy(BoxData* self) {
  delete self;
}

dll void cb_boxdata_clear(BoxData* self, Integer start_modification_id = 0) {
  self->clear(start_modification_id);
}

dll BoxData* cb_boxdata_new(Real fun_factor = 1., int fun_task = ObjectiveFunction) {
  return new BoxData(fun_factor, (FunctionTask)fun_task);
}

dll int cb_boxdata_init(BoxData* self, const BundleData* bd) {
  return self->init(bd);
}

dll BundleData* cb_boxdata_clone(const BoxData* self) {
  return self->clone();
}

dll int cb_boxdata_do_step(BoxData* self, Integer point_id) {
  return self->do_step(point_id);
}

dll int cb_boxdata_synchronize_ids(BoxData* self, Integer* new_center_ub_fid, Integer new_center_id, Integer old_center_id, Integer* new_cand_ub_fid, Integer new_cand_id, Integer old_cand_id, Integer* new_aggregate_id, Integer new_prex_id = 0) {
  return self->synchronize_ids(*new_center_ub_fid, new_center_id, old_center_id, *new_cand_ub_fid, new_cand_id, old_cand_id, *new_aggregate_id, new_prex_id);
}

dll void cb_boxdata_clear_model(BoxData* self, int discard_minorants_only = 0) {
  self->clear_model((bool)discard_minorants_only);
}

dll void cb_boxdata_clear_aggregates(BoxData* self) {
  self->clear_aggregates();
}

dll int cb_boxdata_call_primal_extender(BoxData* self, PrimalExtender* prex, int include_candidates = 1) {
  return self->call_primal_extender(*prex, (bool)include_candidates);
}

dll int cb_boxdata_apply_modification(BoxData* self, const GroundsetModification* param0, MinorantExtender* mex) {
  return self->apply_modification(*param0, mex);
}

dll const PrimalData* cb_boxdata_get_approximate_primal(const BoxData* self) {
  return self->get_approximate_primal();
}

dll const PrimalData* cb_boxdata_get_center_primal(const BoxData* self) {
  return self->get_center_primal();
}

dll const PrimalData* cb_boxdata_get_candidate_primal(const BoxData* self) {
  return self->get_candidate_primal();
}

