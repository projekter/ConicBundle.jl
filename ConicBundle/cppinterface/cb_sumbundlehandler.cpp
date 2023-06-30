dll void cb_sumbundlehandler_destroy(SumBundleHandler* self) {
  delete self;
}

dll const SumBundle* cb_sumbundlehandler_get_sumbundle(const SumBundleHandler* self) {
  return self->get_sumbundle();
}

dll int cb_sumbundlehandler_handles(SumBundleHandler* self, int ft) {
  return self->handles((FunctionTask)ft);
}

dll Integer cb_sumbundlehandler_get_new_index(const SumBundleHandler* self, int ft) {
  return self->get_new_index((FunctionTask)ft);
}

dll int cb_sumbundlehandler_initialization_needed(const SumBundleHandler* self, int ft) {
  return self->initialization_needed((FunctionTask)ft);
}

dll int cb_sumbundlehandler_initialization_needed2(const SumBundleHandler* self) {
  return self->initialization_needed();
}

dll int cb_sumbundlehandler_set_parent_information(SumBundleHandler* self, SumBundleHandler* parent_sbh, const AffineFunctionTransformation* aft, int in_mode) {
  return self->set_parent_information(parent_sbh, aft, (SumBundle::Mode)in_mode);
}

dll int cb_sumbundlehandler_reset_function_factor(SumBundleHandler* self, int ft, Real factor) {
  return self->reset_function_factor((FunctionTask)ft, factor);
}

dll int cb_sumbundlehandler_set_bundle_parameters(SumBundleHandler* self, const BundleParameters* bp) {
  return self->set_bundle_parameters(*bp);
}

dll Real cb_sumbundlehandler_get_increase_factor(const SumBundleHandler* self) {
  return self->get_increase_factor();
}

dll int cb_sumbundlehandler_update_model(SumBundleHandler* self, BundleModel::ModelUpdate* model_update, Integer center_id, const Matrix* center_y, Integer cand_id, const Matrix* cand_y, Real model_maxviol, BundleProxObject* H) {
  return self->update_model(*model_update, center_id, *center_y, cand_id, *cand_y, model_maxviol, *H);
}

dll int cb_sumbundlehandler_eval_model(const SumBundleHandler* self, Real* lb, Integer yid, const Matrix* y) {
  return self->eval_model(*lb, yid, *y);
}

dll Real cb_sumbundlehandler_lb_model(const SumBundleHandler* self, Integer yid, const Matrix* y) {
  return self->lb_model(yid, *y);
}

dll int cb_sumbundlehandler_normalize_sumbundle(SumBundleHandler* self) {
  return self->normalize_sumbundle();
}

dll int cb_sumbundlehandler_contribute_new_minorants(SumBundleHandler* self) {
  return self->contribute_new_minorants();
}

dll int cb_sumbundlehandler_remove_contributions(SumBundleHandler* self) {
  return self->remove_contributions();
}

dll int cb_sumbundlehandler_add_contributions(SumBundleHandler* self) {
  return self->add_contributions();
}

dll int cb_sumbundlehandler_start_augmodel(SumBundleHandler* self, QPModelDataPointer* bp, Integer cand_id, const Matrix* cand_y, const Indexmatrix* indices, int ft) {
  return self->start_augmodel(*bp, cand_id, *cand_y, indices, (FunctionTask)ft);
}

dll int cb_sumbundlehandler_start_augmodel2(SumBundleHandler* self, QPModelDataPointer* bp, QPSumModelDataObject* sumblock, Integer cand_id, const Matrix* cand_y, const Indexmatrix* indices = 0) {
  return self->start_augmodel(*bp, *sumblock, cand_id, *cand_y, indices);
}

dll int cb_sumbundlehandler_make_model_aggregate(SumBundleHandler* self, int* increased, int fixed) {
  increased = 0;
  return self->make_model_aggregate(*(bool*)increased, (bool)fixed);
}

dll int cb_sumbundlehandler_provide_model_aggregate(SumBundleHandler* self) {
  return self->provide_model_aggregate();
}

dll int cb_sumbundlehandler_adjust_multiplier(SumBundleHandler* self, int* values_may_have_changed) {
  values_may_have_changed = 0;
  return self->adjust_multiplier(*(bool*)values_may_have_changed);
}

dll int cb_sumbundlehandler_contribute_initial_bundle(SumBundleHandler* self, int ft, const MinorantBundle* bundle_minorants, const Matrix* coeff) {
  return self->contribute_initial_bundle((FunctionTask)ft, *bundle_minorants, *coeff);
}

dll int cb_sumbundlehandler_install_external_aggregate(SumBundleHandler* self, int ft, const MinorantPointer* aggr, Real aggr_coeff) {
  return self->install_external_aggregate((FunctionTask)ft, *aggr, aggr_coeff);
}

dll int cb_sumbundlehandler_set_cand_minorant(SumBundleHandler* self, int ft, const MinorantPointer* minorant) {
  return self->set_cand_minorant((FunctionTask)ft, *minorant);
}

dll void cb_sumbundlehandler_clear_model(SumBundleHandler* self) {
  self->clear_model();
}

dll void cb_sumbundlehandler_clear_aggregates(SumBundleHandler* self) {
  self->clear_aggregates();
}

dll void cb_sumbundlehandler_clear_cand_minorants(SumBundleHandler* self) {
  self->clear_cand_minorants();
}

dll int cb_sumbundlehandler_add_variable_metric(SumBundleHandler* self, int ft, VariableMetric* H, Integer yid, const Matrix* y, int descent_step, Real weightu, Real model_maxviol, const Indexmatrix* indices = 0) {
  return self->add_variable_metric((FunctionTask)ft, *H, yid, *y, (bool)descent_step, weightu, model_maxviol, indices);
}

dll int cb_sumbundlehandler_add_variable_metric2(SumBundleHandler* self, VariableMetric* H, Integer yid, const Matrix* center_y, int descent_step, Real weightu, Real model_maxviol, const Indexmatrix* indices = 0) {
  return self->add_variable_metric(*H, yid, *center_y, (bool)descent_step, weightu, model_maxviol, indices);
}

dll Real cb_sumbundlehandler_guess_curvature(const SumBundleHandler* self, const MinorantBundle* mnrts, const Indexmatrix* selected_indices, Integer cand_id, const Matrix* cand_y, Real model_maxviol) {
  return self->guess_curvature(*mnrts, *selected_indices, cand_id, *cand_y, model_maxviol);
}

