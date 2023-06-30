dll void cb_sumbundle_destroy(SumBundle* self) {
  delete self;
}

dll SumBundle* cb_sumbundle_new() {
  return new SumBundle();
}

dll SumBundle* cb_sumbundle_new2(const SumBundle* sb) {
  return new SumBundle(*sb);
}

dll void cb_sumbundle_init(SumBundle* self, const SumBundle* sb) {
  self->init(*sb);
}

dll void cb_sumbundle_synchronize_ids(SumBundle* self, Integer new_modification_id, Integer new_center_id, Integer old_center_id, Integer new_cand_id, Integer old_cand_id, Integer new_prex_id = 0) {
  self->synchronize_ids(new_modification_id, new_center_id, old_center_id, new_cand_id, old_cand_id, new_prex_id);
}

dll int cb_sumbundle_has_bundle_for(const SumBundle* self, int ft) {
  return self->has_bundle_for((FunctionTask)ft);
}

dll int cb_sumbundle_has_bundle_data(const SumBundle* self) {
  return self->has_bundle_data();
}

dll int cb_sumbundle_bundle_size(const SumBundle* self, int ft) {
  return self->bundle_size((FunctionTask)ft);
}

dll int cb_sumbundle_has_roots(const SumBundle* self) {
  return self->has_roots();
}

dll int cb_sumbundle_has_working_roots(const SumBundle* self) {
  return self->has_working_roots();
}

dll int cb_sumbundle_active(const SumBundle* self) {
  return self->active();
}

dll int cb_sumbundle_has_contributions(const SumBundle* self) {
  return self->has_contributions();
}

dll int cb_sumbundle_get_mode(const SumBundle* self, int ft) {
  return (int)self->get_mode((FunctionTask)ft);
}

dll Real cb_sumbundle_get_function_factor(const SumBundle* self, int ft) {
  return self->get_function_factor((FunctionTask)ft);
}

dll Integer cb_sumbundle_get_n_contributors(const SumBundle* self, int ft) {
  return self->get_n_contributors((FunctionTask)ft);
}

dll const MinorantBundle* cb_sumbundle_get_bundle(const SumBundle* self, int ft) {
  return &self->get_bundle((FunctionTask)ft);
}

dll const Matrix* cb_sumbundle_get_coeff(const SumBundle* self, int ft) {
  return &self->get_coeff((FunctionTask)ft);
}

dll const MinorantPointer* cb_sumbundle_get_aggregate(const SumBundle* self, int ft) {
  return &self->get_aggregate((FunctionTask)ft);
}

dll const MinorantPointer* cb_sumbundle_get_cand_minorant(const SumBundle* self, int ft) {
  return &self->get_cand_minorant((FunctionTask)ft);
}

dll int cb_sumbundle_get_local_model_aggregate(const SumBundle* self, MinorantPointer* aggregate, Real factor = 1., const AffineFunctionTransformation* aft = 0) {
  return self->get_local_model_aggregate(*aggregate, factor, aft);
}

dll int cb_sumbundle_get_contributed_model_aggregate(const SumBundle* self, MinorantPointer* aggregate, Real factor = 1., const AffineFunctionTransformation* aft = 0) {
  return self->get_contributed_model_aggregate(*aggregate, factor, aft);
}

dll int cb_sumbundle_call_primal_extender(SumBundle* self, PrimalExtender* prex, Integer prex_id, int ft) {
  return self->call_primal_extender(*prex, prex_id, (FunctionTask)ft);
}

dll int cb_sumbundle_apply_modification(SumBundle* self, const GroundsetModification* gsmdf, Integer mod_id, MinorantExtender* mex, int ft) {
  return self->apply_modification(*gsmdf, mod_id, mex, (FunctionTask)ft);
}

