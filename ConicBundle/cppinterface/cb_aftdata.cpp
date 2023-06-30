dll void cb_aftdata_destroy(AFTData* self) {
  delete self;
}

dll void cb_aftdata_clear(AFTData* self, Integer start_modification_id = 0) {
  self->clear(start_modification_id);
}

dll AFTData* cb_aftdata_new(Integer start_modification_id = 0) {
  return new AFTData(start_modification_id);
}

dll int cb_aftdata_init(AFTData* self, const BundleData* bd) {
  return self->init(bd);
}

dll BundleData* cb_aftdata_clone(const AFTData* self) {
  return self->clone();
}

dll int cb_aftdata_do_step(AFTData* self, Integer point_id) {
  return self->do_step(point_id);
}

dll int cb_aftdata_synchronize_ids(AFTData* self, Integer* new_center_ub_fid, Integer new_center_id, Integer old_center_id, Integer* new_cand_ub_fid, Integer new_cand_id, Integer old_cand_id, Integer* new_aggregate_id, Integer new_prex_id = 0) {
  return self->synchronize_ids(*new_center_ub_fid, new_center_id, old_center_id, *new_cand_ub_fid, new_cand_id, old_cand_id, *new_aggregate_id, new_prex_id);
}

dll int cb_aftdata_center_modified(AFTData* self, Integer* function_id, Integer center_id) {
  return self->center_modified(*function_id, center_id);
}

dll int cb_aftdata_model_aggregate_modified(AFTData* self, Integer last_aggr_id) {
  return self->model_aggregate_modified(last_aggr_id);
}

dll void cb_aftdata_clear_model(AFTData* self, int discard_minorants_only = 0) {
  self->clear_model((bool)discard_minorants_only);
}

dll void cb_aftdata_clear_aggregates(AFTData* self) {
  self->clear_aggregates();
}

