dll void cb_socdata_destroy(SOCData* self) {
  delete self;
}

dll void cb_socdata_clear(SOCData* self, Integer start_modification_id = 0) {
  self->clear(start_modification_id);
}

dll SOCData* cb_socdata_new(Real fun_factor = 1., int fun_task = ObjectiveFunction) {
  return new SOCData(fun_factor, (FunctionTask)fun_task);
}

dll int cb_socdata_init(SOCData* self, const BundleData* bd) {
  return self->init(bd);
}

dll BundleData* cb_socdata_clone(const SOCData* self) {
  return self->clone();
}

dll int cb_socdata_do_step(SOCData* self, Integer point_id) {
  return self->do_step(point_id);
}

dll int cb_socdata_store_socvec(SOCData* self, const Matrix* SOCvec) {
  return self->store_SOCvec(*SOCvec);
}

dll int cb_socdata_form_bundlevecs(SOCData* self, Integer max_columns) {
  return self->form_bundlevecs(max_columns);
}

dll int cb_socdata_synchronize_ids(SOCData* self, Integer* new_center_ub_fid, Integer new_center_id, Integer old_center_id, Integer* new_cand_ub_fid, Integer new_cand_id, Integer old_cand_id, Integer* new_aggregate_id, Integer new_prex_id = 0) {
  return self->synchronize_ids(*new_center_ub_fid, new_center_id, old_center_id, *new_cand_ub_fid, new_cand_id, old_cand_id, *new_aggregate_id, new_prex_id);
}

dll void cb_socdata_clear_model(SOCData* self, int discard_minorants_only = 0) {
  self->clear_model((bool)discard_minorants_only);
}

dll void cb_socdata_clear_aggregates(SOCData* self) {
  self->clear_aggregates();
}

dll int cb_socdata_call_primal_extender(SOCData* self, PrimalExtender* param0, int include_candidates = 1) {
  return self->call_primal_extender(*param0, (bool)include_candidates);
}

dll int cb_socdata_apply_modification(SOCData* self, const GroundsetModification* param0, MinorantExtender* mex) {
  return self->apply_modification(*param0, mex);
}

dll const PrimalData* cb_socdata_get_approximate_primal(const SOCData* self) {
  return self->get_approximate_primal();
}

dll const PrimalData* cb_socdata_get_center_primal(const SOCData* self) {
  return self->get_center_primal();
}

dll const PrimalData* cb_socdata_get_candidate_primal(const SOCData* self) {
  return self->get_candidate_primal();
}

dll int cb_socdata_get_latest_minorants(SOCData* self, MinorantBundle* latest_minorants, Integer max_number) {
  return self->get_latest_minorants(*latest_minorants, max_number);
}

