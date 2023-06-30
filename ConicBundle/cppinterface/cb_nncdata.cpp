dll void cb_nncdata_destroy(NNCData* self) {
  delete self;
}

dll void cb_nncdata_clear(NNCData* self, Integer start_modification_id = 0) {
  self->clear(start_modification_id);
}

dll NNCData* cb_nncdata_new(Real fun_factor = 1., int fun_task = ObjectiveFunction) {
  return new NNCData(fun_factor, (FunctionTask)fun_task);
}

dll int cb_nncdata_init(NNCData* self, const BundleData* bd) {
  return self->init(bd);
}

dll BundleData* cb_nncdata_clone(const NNCData* self) {
  return self->clone();
}

dll int cb_nncdata_do_step(NNCData* self, Integer point_id) {
  return self->do_step(point_id);
}

dll int cb_nncdata_synchronize_ids(NNCData* self, Integer* new_center_ub_fid, Integer new_center_id, Integer old_center_id, Integer* new_cand_ub_fid, Integer new_cand_id, Integer old_cand_id, Integer* new_aggregate_id, Integer new_prex_id = 0) {
  return self->synchronize_ids(*new_center_ub_fid, new_center_id, old_center_id, *new_cand_ub_fid, new_cand_id, old_cand_id, *new_aggregate_id, new_prex_id);
}

dll void cb_nncdata_clear_model(NNCData* self, int discard_minorants_only = 0) {
  self->clear_model((bool)discard_minorants_only);
}

dll void cb_nncdata_clear_aggregates(NNCData* self) {
  self->clear_aggregates();
}

dll int cb_nncdata_call_primal_extender(NNCData* self, PrimalExtender* prex, int include_candidates = 1) {
  return self->call_primal_extender(*prex, (bool)include_candidates);
}

dll int cb_nncdata_apply_modification(NNCData* self, const GroundsetModification* param0, MinorantExtender* mex) {
  return self->apply_modification(*param0, mex);
}

dll const PrimalData* cb_nncdata_get_approximate_primal(const NNCData* self) {
  return self->get_approximate_primal();
}

dll const PrimalData* cb_nncdata_get_center_primal(const NNCData* self) {
  return self->get_center_primal();
}

dll const PrimalData* cb_nncdata_get_candidate_primal(const NNCData* self) {
  return self->get_candidate_primal();
}

dll int cb_nncdata_get_model_data(const NNCData* self, MinorantBundle* model_minorants, Matrix* model_coeff) {
  return self->get_model_data(*model_minorants, *model_coeff);
}

