dll void cb_bundlelowranktrustregionprox_destroy(BundleLowRankTrustRegionProx* self) {
  delete self;
}

dll BundleLowRankTrustRegionProx* cb_bundlelowranktrustregionprox_new(Integer dim = 0, VariableMetricSelection* vp = 0, int local_metric = 0, int inc = -1) {
  return new BundleLowRankTrustRegionProx(dim, vp, (bool)local_metric, 0, inc);
}

dll BundleLowRankTrustRegionProx* cb_bundlelowranktrustregionprox_new2(const Matrix* in_vecH, const Matrix* in_lamH, VariableMetricSelection* vp = 0, int local_metric = 0, int inc = -1) {
  return new BundleLowRankTrustRegionProx(*in_vecH, *in_lamH, vp, (bool)local_metric, 0, inc);
}

dll void cb_bundlelowranktrustregionprox_set_weightu(BundleLowRankTrustRegionProx* self, Real in_weightu) {
  self->set_weightu(in_weightu);
}

dll Real cb_bundlelowranktrustregionprox_get_weightu(const BundleLowRankTrustRegionProx* self) {
  return self->get_weightu();
}

dll Real cb_bundlelowranktrustregionprox_get_term_corr(const BundleLowRankTrustRegionProx* self) {
  return self->get_term_corr();
}

dll void cb_bundlelowranktrustregionprox_init(BundleLowRankTrustRegionProx* self, const Matrix* in_vecH, const Matrix* in_lamH) {
  self->init(*in_vecH, *in_lamH);
}

dll Real cb_bundlelowranktrustregionprox_norm_sqr(const BundleLowRankTrustRegionProx* self, const Matrix* B) {
  return self->norm_sqr(*B);
}

dll Real cb_bundlelowranktrustregionprox_dnorm_sqr(const BundleLowRankTrustRegionProx* self, const MinorantPointer* B) {
  return self->dnorm_sqr(*B);
}

dll int cb_bundlelowranktrustregionprox_is_dlr(const BundleLowRankTrustRegionProx* self) {
  return self->is_DLR();
}

dll int cb_bundlelowranktrustregionprox_add_h(const BundleLowRankTrustRegionProx* self, Symmatrix* big_sym, Integer start_index = 0) {
  return self->add_H(*big_sym, start_index);
}

dll Matrix* cb_bundlelowranktrustregionprox_add_hx(const BundleLowRankTrustRegionProx* self, const Matrix* x, Matrix* outplusHx, Real alpha = 1.) {
  return &self->add_Hx(*x, *outplusHx, alpha);
}

dll Matrix* cb_bundlelowranktrustregionprox_apply_hinv(const BundleLowRankTrustRegionProx* self, Matrix* x) {
  return &self->apply_Hinv(*x);
}

dll int cb_bundlelowranktrustregionprox_compute_qp_costs(BundleLowRankTrustRegionProx* self, Symmatrix* Q, Matrix* d, Real* offset, const MinorantPointer* constant_minorant, const MinorantBundle* bundle, const Matrix* y, const MinorantPointer* groundset_minorant, Indexmatrix* yfixed) {
  return self->compute_QP_costs(*Q, *d, *offset, *constant_minorant, *bundle, *y, *groundset_minorant, yfixed);
}

dll int cb_bundlelowranktrustregionprox_update_qp_costs(BundleLowRankTrustRegionProx* self, Symmatrix* delta_Q, Matrix* delta_d, Real* delta_offset, const MinorantPointer* constant_minorant, const MinorantBundle* bundle, const Matrix* center_y, const MinorantPointer* groundset_minorant, const MinorantPointer* delta_groundset_minorant, const Indexmatrix* delta_index, Indexmatrix* yfixed) {
  return self->update_QP_costs(*delta_Q, *delta_d, *delta_offset, *constant_minorant, *bundle, *center_y, *groundset_minorant, *delta_groundset_minorant, *delta_index, yfixed);
}

dll int cb_bundlelowranktrustregionprox_apply_modification(BundleLowRankTrustRegionProx* self, const GroundsetModification* gsmdf) {
  return self->apply_modification(*gsmdf);
}

dll BundleProxObject* cb_bundlelowranktrustregionprox_projected_clone(BundleLowRankTrustRegionProx* self, const Indexmatrix* indices) {
  return self->projected_clone(*indices);
}

dll int cb_bundlelowranktrustregionprox_supports_diagonal_bounds_scaling(const BundleLowRankTrustRegionProx* self) {
  return self->supports_diagonal_bounds_scaling();
}

dll int cb_bundlelowranktrustregionprox_diagonal_scaling_heuristic_update(BundleLowRankTrustRegionProx* self, const Matrix* param0) {
  return self->diagonal_scaling_heuristic_update(*param0);
}

dll int cb_bundlelowranktrustregionprox_supports_dense_variable_metric(const BundleLowRankTrustRegionProx* self) {
  return self->supports_dense_variable_metric();
}

dll int cb_bundlelowranktrustregionprox_supports_lowrank_variable_metric(const BundleLowRankTrustRegionProx* self) {
  return self->supports_lowrank_variable_metric();
}

dll int cb_bundlelowranktrustregionprox_supports_diagonal_variable_metric(const BundleLowRankTrustRegionProx* self) {
  return self->supports_diagonal_variable_metric();
}

dll int cb_bundlelowranktrustregionprox_apply_variable_metric(BundleLowRankTrustRegionProx* self, VariableMetricModel* groundset, VariableMetricModel* model, const Matrix* aggr, Integer y_id, const Matrix* y, int descent_step, Real* current_weight, Real model_maxviol, const Indexmatrix* new_indices = 0) {
  return self->apply_variable_metric(groundset, model, *aggr, y_id, *y, (bool)descent_step, *current_weight, model_maxviol, new_indices);
}

dll int cb_bundlelowranktrustregionprox_add_variable_metric(BundleLowRankTrustRegionProx* self, Matrix* diagH, Matrix* vecH) {
  return self->add_variable_metric(*diagH, *vecH);
}

dll int cb_bundlelowranktrustregionprox_push_aft(BundleLowRankTrustRegionProx* self, const AffineFunctionTransformation* aft) {
  return self->push_aft(aft);
}

dll int cb_bundlelowranktrustregionprox_pop_aft(BundleLowRankTrustRegionProx* self) {
  return self->pop_aft();
}

dll int cb_bundlelowranktrustregionprox_mfile_data(const BundleLowRankTrustRegionProx* self) {
  return self->mfile_data(std::cout);
}

