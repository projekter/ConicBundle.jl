dll void cb_summodel_destroy(SumModel* self) {
  delete self;
}

dll void cb_summodel_clear(SumModel* self) {
  self->clear();
}

dll SumModel* cb_summodel_new() {
  return new SumModel(0);
}

dll const SumBlockModel* cb_summodel_model(const SumModel* self, const FunctionObject* fo) {
  return self->model(fo);
}

dll Integer cb_summodel_nsubmodels(const SumModel* self) {
  return self->nsubmodels();
}

dll int cb_summodel_add_model(SumModel* self, SumBlockModel* model) {
  return self->add_model(model);
}

dll SumBlockModel* cb_summodel_remove_model(SumModel* self, const FunctionObject* fo) {
  return self->remove_model(fo);
}

dll int cb_summodel_eval_function(SumModel* self, Integer* ub_fid, Real* ub, Integer y_id, const Matrix* y, Real nullstep_bound, Real relprec) {
  return self->eval_function(*ub_fid, *ub, y_id, *y, nullstep_bound, relprec);
}

dll int cb_summodel_eval_model(SumModel* self, Real* lb, Integer y_id, const Matrix* y, Real relprec) {
  return self->eval_model(*lb, y_id, *y, relprec);
}

dll int cb_summodel_update_model(SumModel* self, int model_update, Integer center_id, const Matrix* center_y, Integer cand_id, const Matrix* cand_y, Real model_maxviol, BundleProxObject* H) {
  return self->update_model((ModelUpdate)model_update, center_id, *center_y, cand_id, *cand_y, model_maxviol, *H);
}

dll int cb_summodel_synchronize_ids(SumModel* self, Integer* new_center_ub_fid, Integer new_center_id, Integer old_center_id, Integer* new_cand_ub_fid, Integer new_cand_id, Integer old_cand_id, Integer* new_aggregate_id) {
  return self->synchronize_ids(*new_center_ub_fid, new_center_id, old_center_id, *new_cand_ub_fid, new_cand_id, old_cand_id, *new_aggregate_id);
}

dll int cb_summodel_center_modified(SumModel* self, Integer* function_modification_id, Integer center_id) {
  return self->center_modified(*function_modification_id, center_id);
}

dll int cb_summodel_recompute_center(SumModel* self, Integer* new_center_ub_fid, Real* new_center_ub, Integer center_id, const Matrix* y, int accept_only_higher_values = 0, Real relprec = -1.) {
  return self->recompute_center(*new_center_ub_fid, *new_center_ub, center_id, *y, (bool)accept_only_higher_values, relprec);
}

dll int cb_summodel_model_aggregate_modified(SumModel* self, Integer old_model_aggregate_id) {
  return self->model_aggregate_modified(old_model_aggregate_id);
}

dll int cb_summodel_provide_model_aggregate(SumModel* self, Integer y_id, const Matrix* y) {
  return self->provide_model_aggregate(y_id, *y);
}

dll int cb_summodel_add_variable_metric(SumModel* self, VariableMetric* H, Integer y_id, const Matrix* y, int descent_step, Real weightu, Real model_maxviol, const Indexmatrix* indices = 0) {
  return self->add_variable_metric(*H, y_id, *y, (bool)descent_step, weightu, model_maxviol, indices);
}

dll int cb_summodel_check_center_validity_by_candidate(SumModel* self, int* cand_minorant_is_below, Integer center_id, const Matrix* center_y) {
  cand_minorant_is_below = 0;
  return self->check_center_validity_by_candidate(*(bool*)cand_minorant_is_below, center_id, *center_y);
}

dll ModifiableOracleObject* cb_summodel_get_oracle_object(SumModel* self) {
  return self->get_oracle_object();
}

dll int cb_summodel_make_model_aggregate(SumModel* self, int* penalty_parameter_increased, int keep_penalty_fixed) {
  penalty_parameter_increased = 0;
  return self->make_model_aggregate(*(bool*)penalty_parameter_increased, (bool)keep_penalty_fixed);
}

dll Real cb_summodel_lb_function(SumModel* self, Integer y_id, const Matrix* y) {
  return self->lb_function(y_id, *y);
}

dll int cb_summodel_get_function_minorant(SumModel* self, MinorantPointer* minorant, const AffineFunctionTransformation* aft = 0) {
  return self->get_function_minorant(*minorant, aft);
}

dll int cb_summodel_get_center_minorant(SumModel* self, MinorantPointer* minorant, const AffineFunctionTransformation* aft = 0) {
  return self->get_center_minorant(*minorant, aft);
}

dll int cb_summodel_adjust_multiplier(SumModel* self, int* values_may_have_changed) {
  values_may_have_changed = 0;
  return self->adjust_multiplier(*(bool*)values_may_have_changed);
}

dll int cb_summodel_sumbundle_mode(SumModel* self, int* mode, SumBundleHandler* bh = 0, AffineFunctionTransformation* aft = 0) {
  return self->sumbundle_mode(*(SumBundle::Mode*)mode, bh, aft);
}

dll int cb_summodel_start_sumaugmodel(SumModel* self, QPModelDataPointer* blockp, Integer cand_id, const Matrix* cand_y, const Indexmatrix* indices = 0, SumBundleHandler* bh = 0, int mode = SumBundle::inactive, AffineFunctionTransformation* aft = 0) {
  return self->start_sumaugmodel(*blockp, cand_id, *cand_y, indices, bh, (SumBundle::Mode)mode, aft);
}

dll int cb_summodel_update_model2(SumModel* self, int model_update, Integer center_id, const Matrix* center_y, Integer cand_id, const Matrix* cand_y, Real null_step_factor, BundleProxObject* H, Real* model_deviation, Real* model_curvature) {
  return self->update_model((ModelUpdate)model_update, center_id, *center_y, cand_id, *cand_y, null_step_factor, *H, *model_deviation, *model_curvature);
}

dll BundleData* cb_summodel_get_data(SumModel* self) {
  return self->get_data();
}

dll const BundleData* cb_summodel_get_data2(const SumModel* self) {
  return self->get_data();
}

dll int cb_summodel_set_data(SumModel* self, BundleData* bd) {
  return self->set_data(bd);
}

dll const PrimalData* cb_summodel_get_approximate_primal(const SumModel* self) {
  return self->get_approximate_primal();
}

dll const PrimalData* cb_summodel_get_center_primal(const SumModel* self) {
  return self->get_center_primal();
}

dll const PrimalData* cb_summodel_get_candidate_primal(const SumModel* self) {
  return self->get_candidate_primal();
}

dll int cb_summodel_call_primal_extender(SumModel* self, PrimalExtender* param0) {
  return self->call_primal_extender(*param0);
}

dll int cb_summodel_set_bundle_parameters(SumModel* self, const BundleParameters* bp) {
  return self->set_bundle_parameters(*bp);
}

dll const BundleParameters* cb_summodel_get_bundle_parameters(const SumModel* self) {
  return self->get_bundle_parameters();
}

dll void cb_summodel_clear_model(SumModel* self, int discard_minorants_only = 0) {
  self->clear_model((bool)discard_minorants_only);
}

dll void cb_summodel_clear_aggregates(SumModel* self) {
  self->clear_aggregates();
}

dll void cb_summodel_set_out(SumModel* self, int pril = 1) {
  self->set_out(&std::cout, pril);
}

dll CH_Tools::Microseconds* cb_summodel_new_get_preeval_time(const SumModel* self) {
  return new CH_Tools::Microseconds(self->get_preeval_time());
}

dll CH_Tools::Microseconds* cb_summodel_new_get_eval_time(const SumModel* self) {
  return new CH_Tools::Microseconds(self->get_eval_time());
}

dll CH_Tools::Microseconds* cb_summodel_new_get_posteval_time(const SumModel* self) {
  return new CH_Tools::Microseconds(self->get_posteval_time());
}

