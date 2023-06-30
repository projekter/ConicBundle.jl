dll void cb_aftmodel_destroy(AFTModel* self) {
  delete self;
}

dll void cb_aftmodel_clear(AFTModel* self, AffineFunctionTransformation* inaft, Integer start_modification_id = 0) {
  self->clear(inaft, start_modification_id);
}

dll void cb_aftmodel_clear2(AFTModel* self) {
  self->clear();
}

dll AFTModel* cb_aftmodel_new(SumBlockModel* in_model, AffineFunctionTransformation* inaft = 0, Integer start_modification_id = 0, int in_model_is_owner = 0) {
  return new AFTModel(in_model, inaft, start_modification_id, (bool)in_model_is_owner, 0);
}

dll const AffineFunctionTransformation* cb_aftmodel_get_aft(const AFTModel* self) {
  return self->get_aft();
}

dll int cb_aftmodel_eval_function(AFTModel* self, Integer* ub_fid, Real* ub, Integer y_id, const Matrix* y, Real nullstep_bound, Real relprec) {
  return self->eval_function(*ub_fid, *ub, y_id, *y, nullstep_bound, relprec);
}

dll int cb_aftmodel_eval_model(AFTModel* self, Real* lb, Integer y_id, const Matrix* y, Real relprec) {
  return self->eval_model(*lb, y_id, *y, relprec);
}

dll int cb_aftmodel_update_model(AFTModel* self, int model_update, Integer center_id, const Matrix* center_y, Integer cand_id, const Matrix* cand_y, Real model_maxviol, BundleProxObject* H) {
  return self->update_model((ModelUpdate)model_update, center_id, *center_y, cand_id, *cand_y, model_maxviol, *H);
}

dll int cb_aftmodel_synchronize_ids(AFTModel* self, Integer* new_center_ub_fid, Integer new_center_id, Integer old_center_id, Integer* new_cand_ub_fid, Integer new_cand_id, Integer old_cand_id, Integer* new_aggregate_id) {
  return self->synchronize_ids(*new_center_ub_fid, new_center_id, old_center_id, *new_cand_ub_fid, new_cand_id, old_cand_id, *new_aggregate_id);
}

dll int cb_aftmodel_center_modified(AFTModel* self, Integer* center_ub_fid, Integer center_id) {
  return self->center_modified(*center_ub_fid, center_id);
}

dll int cb_aftmodel_recompute_center(AFTModel* self, Integer* new_center_ub_fid, Real* new_center_ub, Integer center_id, const Matrix* y, int accept_only_higher_values = 0, Real relprec = -1.) {
  return self->recompute_center(*new_center_ub_fid, *new_center_ub, center_id, *y, (bool)accept_only_higher_values, relprec);
}

dll int cb_aftmodel_model_aggregate_modified(AFTModel* self, Integer old_model_aggregate_id) {
  return self->model_aggregate_modified(old_model_aggregate_id);
}

dll int cb_aftmodel_provide_model_aggregate(AFTModel* self, Integer y_id, const Matrix* y) {
  return self->provide_model_aggregate(y_id, *y);
}

dll int cb_aftmodel_add_variable_metric(AFTModel* self, VariableMetric* H, Integer center_id, const Matrix* y, int descent_step, Real weightu, Real model_maxviol, const Indexmatrix* indices = 0) {
  return self->add_variable_metric(*H, center_id, *y, (bool)descent_step, weightu, model_maxviol, indices);
}

dll int cb_aftmodel_check_center_validity_by_candidate(AFTModel* self, int* cand_minorant_is_below, Integer center_id, const Matrix* center_y) {
  cand_minorant_is_below = 0;
  return self->check_center_validity_by_candidate(*(bool*)cand_minorant_is_below, center_id, *center_y);
}

dll ModifiableOracleObject* cb_aftmodel_get_oracle_object(AFTModel* self) {
  return self->get_oracle_object();
}

dll int cb_aftmodel_make_model_aggregate(AFTModel* self, int* penalty_parameter_increased, int keep_penalty_fixed) {
  penalty_parameter_increased = 0;
  return self->make_model_aggregate(*(bool*)penalty_parameter_increased, (bool)keep_penalty_fixed);
}

dll int cb_aftmodel_get_model_aggregate(AFTModel* self, Integer* model_aggregate_id, MinorantPointer* model_aggregate, int all_parts = 1, const AffineFunctionTransformation* aft = 0) {
  return self->get_model_aggregate(*model_aggregate_id, *model_aggregate, (bool)all_parts, aft);
}

dll Real cb_aftmodel_lb_function(AFTModel* self, Integer y_id, const Matrix* y) {
  return self->lb_function(y_id, *y);
}

dll int cb_aftmodel_get_function_minorant(AFTModel* self, MinorantPointer* minorant, const AffineFunctionTransformation* aft = 0) {
  return self->get_function_minorant(*minorant, aft);
}

dll int cb_aftmodel_get_center_minorant(AFTModel* self, MinorantPointer* minorant, const AffineFunctionTransformation* aft = 0) {
  return self->get_center_minorant(*minorant, aft);
}

dll int cb_aftmodel_adjust_multiplier(AFTModel* self, int* values_may_have_changed) {
  values_may_have_changed = 0;
  return self->adjust_multiplier(*(bool*)values_may_have_changed);
}

dll int cb_aftmodel_sumbundle_mode(AFTModel* self, int* mode, SumBundleHandler* bh = 0, AffineFunctionTransformation* aft = 0) {
  return self->sumbundle_mode(*(SumBundle::Mode*)mode, bh, aft);
}

dll int cb_aftmodel_start_sumaugmodel(AFTModel* self, QPModelDataPointer* blockp, Integer cand_id, const Matrix* cand_y, const Indexmatrix* indices = 0, SumBundleHandler* bh = 0, int mode = SumBundle::inactive, AffineFunctionTransformation* aft = 0) {
  return self->start_sumaugmodel(*blockp, cand_id, *cand_y, indices, bh, (SumBundle::Mode)mode, aft);
}

dll int cb_aftmodel_update_model2(AFTModel* self, int model_update, Integer center_id, const Matrix* center_y, Integer cand_id, const Matrix* cand_y, Real model_maxviol, BundleProxObject* H, Real* model_deviation, Real* model_curvature) {
  return self->update_model((ModelUpdate)model_update, center_id, *center_y, cand_id, *cand_y, model_maxviol, *H, *model_deviation, *model_curvature);
}

dll BundleData* cb_aftmodel_get_data(AFTModel* self) {
  return self->get_data();
}

dll const BundleData* cb_aftmodel_get_data2(const AFTModel* self) {
  return self->get_data();
}

dll int cb_aftmodel_set_data(AFTModel* self, BundleData* bd) {
  return self->set_data(bd);
}

dll const PrimalData* cb_aftmodel_get_approximate_primal(const AFTModel* self) {
  return self->get_approximate_primal();
}

dll const PrimalData* cb_aftmodel_get_center_primal(const AFTModel* self) {
  return self->get_center_primal();
}

dll const PrimalData* cb_aftmodel_get_candidate_primal(const AFTModel* self) {
  return self->get_candidate_primal();
}

dll int cb_aftmodel_call_primal_extender(AFTModel* self, PrimalExtender* param0) {
  return self->call_primal_extender(*param0);
}

dll int cb_aftmodel_set_bundle_parameters(AFTModel* self, const BundleParameters* param0) {
  return self->set_bundle_parameters(*param0);
}

dll BundleParameters* cb_aftmodel_get_bundle_parameters(const AFTModel* self) {
  return self->get_bundle_parameters();
}

dll void cb_aftmodel_clear_model(AFTModel* self, int discard_minorants_only = 0) {
  self->clear_model((bool)discard_minorants_only);
}

dll void cb_aftmodel_clear_aggregates(AFTModel* self) {
  self->clear_aggregates();
}

dll CH_Tools::Microseconds* cb_aftmodel_new_get_preeval_time(const AFTModel* self) {
  return new CH_Tools::Microseconds(self->get_preeval_time());
}

dll CH_Tools::Microseconds* cb_aftmodel_new_get_eval_time(const AFTModel* self) {
  return new CH_Tools::Microseconds(self->get_eval_time());
}

dll CH_Tools::Microseconds* cb_aftmodel_new_get_posteval_time(const AFTModel* self) {
  return new CH_Tools::Microseconds(self->get_posteval_time());
}

dll void cb_aftmodel_set_out(AFTModel* self, int pril = 1) {
  self->set_out(&std::cout, pril);
}

