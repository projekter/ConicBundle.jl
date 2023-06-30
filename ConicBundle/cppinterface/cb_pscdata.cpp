dll void cb_pscdata_destroy(PSCData* self) {
  delete self;
}

dll void cb_pscdata_clear(PSCData* self, Integer start_modification_id = 0) {
  self->clear(start_modification_id);
}

dll PSCData* cb_pscdata_new(Real fun_factor = 1., int fun_task = ObjectiveFunction) {
  return new PSCData(fun_factor, (FunctionTask)fun_task);
}

dll const Matrix* cb_pscdata_get_primalvecs(const PSCData* self) {
  return &self->get_primalvecs();
}

dll const Matrix* cb_pscdata_get_primaleigs(const PSCData* self) {
  return &self->get_primaleigs();
}

dll const Matrix* cb_pscdata_get_topvecs(const PSCData* self) {
  return &self->get_topvecs();
}

dll const Matrix* cb_pscdata_get_ritz_values(const PSCData* self) {
  return &self->get_Ritz_values();
}

dll Integer cb_pscdata_get_keepsize(const PSCData* self) {
  return self->get_keepsize();
}

dll Integer cb_pscdata_get_activedim(const PSCData* self) {
  return self->get_activedim();
}

dll Integer cb_pscdata_get_skippedsize(const PSCData* self) {
  return self->get_skippedsize();
}

dll int cb_pscdata_init(PSCData* self, const BundleData* bd) {
  return self->init(bd);
}

dll BundleData* cb_pscdata_clone(const PSCData* self) {
  return self->clone();
}

dll int cb_pscdata_do_step(PSCData* self, Integer point_id) {
  return self->do_step(point_id);
}

dll int cb_pscdata_synchronize_ids(PSCData* self, Integer* new_center_ub_fid, Integer new_center_id, Integer old_center_id, Integer* new_cand_ub_fid, Integer new_cand_id, Integer old_cand_id, Integer* new_aggregate_id, Integer new_prex_id = 0) {
  return self->synchronize_ids(*new_center_ub_fid, new_center_id, old_center_id, *new_cand_ub_fid, new_cand_id, old_cand_id, *new_aggregate_id, new_prex_id);
}

dll void cb_pscdata_clear_model(PSCData* self, int discard_minorants_only = 0) {
  self->clear_model((bool)discard_minorants_only);
}

dll void cb_pscdata_clear_aggregates(PSCData* self) {
  self->clear_aggregates();
}

dll int cb_pscdata_call_primal_extender(PSCData* self, PrimalExtender* param0, int include_candidates = 1) {
  return self->call_primal_extender(*param0, (bool)include_candidates);
}

dll int cb_pscdata_apply_modification(PSCData* self, const GroundsetModification* param0, MinorantExtender* mex) {
  return self->apply_modification(*param0, mex);
}

dll const PrimalData* cb_pscdata_get_approximate_primal(const PSCData* self) {
  return self->get_approximate_primal();
}

dll const PrimalData* cb_pscdata_get_center_primal(const PSCData* self) {
  return self->get_center_primal();
}

dll const PrimalData* cb_pscdata_get_candidate_primal(const PSCData* self) {
  return self->get_candidate_primal();
}

