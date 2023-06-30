dll void cb_bundledlrtrustregionprox_destroy(BundleDLRTrustRegionProx* self) {
  delete self;
}

dll BundleDLRTrustRegionProx* cb_bundledlrtrustregionprox_new(Integer dim = 0, VariableMetricSelection* vp = 0, int local_metric = 0, int inc = -1) {
  return new BundleDLRTrustRegionProx(dim, vp, (bool)local_metric, 0, inc);
}

dll BundleDLRTrustRegionProx* cb_bundledlrtrustregionprox_new2(const Matrix* in_D, const Matrix* in_vecH, VariableMetricSelection* vp = 0, int local_metric = 0, int inc = -1) {
  return new BundleDLRTrustRegionProx(*in_D, *in_vecH, vp, (bool)local_metric, 0, inc);
}

dll void cb_bundledlrtrustregionprox_set_weightu(BundleDLRTrustRegionProx* self, Real in_weightu) {
  self->set_weightu(in_weightu);
}

dll Real cb_bundledlrtrustregionprox_get_weightu(const BundleDLRTrustRegionProx* self) {
  return self->get_weightu();
}

dll Real cb_bundledlrtrustregionprox_get_term_corr(const BundleDLRTrustRegionProx* self) {
  return self->get_term_corr();
}

dll void cb_bundledlrtrustregionprox_init(BundleDLRTrustRegionProx* self, const Matrix* in_D, const Matrix* in_vecH) {
  self->init(*in_D, *in_vecH);
}

dll Real cb_bundledlrtrustregionprox_norm_sqr(const BundleDLRTrustRegionProx* self, const Matrix* B) {
  return self->norm_sqr(*B);
}

dll Real cb_bundledlrtrustregionprox_dnorm_sqr(const BundleDLRTrustRegionProx* self, const MinorantPointer* B) {
  return self->dnorm_sqr(*B);
}

dll int cb_bundledlrtrustregionprox_is_dlr(const BundleDLRTrustRegionProx* self) {
  return self->is_DLR();
}

dll int cb_bundledlrtrustregionprox_add_h(const BundleDLRTrustRegionProx* self, Symmatrix* big_sym, Integer start_index = 0) {
  return self->add_H(*big_sym, start_index);
}

dll Matrix* cb_bundledlrtrustregionprox_add_hx(const BundleDLRTrustRegionProx* self, const Matrix* x, Matrix* outplusHx, Real alpha = 1.) {
  return &self->add_Hx(*x, *outplusHx, alpha);
}

dll Matrix* cb_bundledlrtrustregionprox_apply_hinv(const BundleDLRTrustRegionProx* self, Matrix* x) {
  return &self->apply_Hinv(*x);
}

dll int cb_bundledlrtrustregionprox_compute_qp_costs(BundleDLRTrustRegionProx* self, Symmatrix* Q, Matrix* d, Real* offset, const MinorantPointer* constant_minorant, const MinorantBundle* bundle, const Matrix* y, const MinorantPointer* groundset_minorant, Indexmatrix* yfixed) {
  return self->compute_QP_costs(*Q, *d, *offset, *constant_minorant, *bundle, *y, *groundset_minorant, yfixed);
}

dll int cb_bundledlrtrustregionprox_update_qp_costs(BundleDLRTrustRegionProx* self, Symmatrix* delta_Q, Matrix* delta_d, Real* delta_offset, const MinorantPointer* constant_minorant, const MinorantBundle* bundle, const Matrix* center_y, const MinorantPointer* groundset_minorant, const MinorantPointer* delta_groundset_minorant, const Indexmatrix* delta_index, Indexmatrix* yfixed) {
  return self->update_QP_costs(*delta_Q, *delta_d, *delta_offset, *constant_minorant, *bundle, *center_y, *groundset_minorant, *delta_groundset_minorant, *delta_index, yfixed);
}

dll int cb_bundledlrtrustregionprox_apply_modification(BundleDLRTrustRegionProx* self, const GroundsetModification* gsmdf) {
  return self->apply_modification(*gsmdf);
}

dll BundleProxObject* cb_bundledlrtrustregionprox_projected_clone(BundleDLRTrustRegionProx* self, const Indexmatrix* indices) {
  return self->projected_clone(*indices);
}

dll int cb_bundledlrtrustregionprox_supports_diagonal_bounds_scaling(const BundleDLRTrustRegionProx* self) {
  return self->supports_diagonal_bounds_scaling();
}

dll int cb_bundledlrtrustregionprox_diagonal_bounds_scaling_update(BundleDLRTrustRegionProx* self, const Matrix* param0) {
  return self->diagonal_bounds_scaling_update(*param0);
}

dll int cb_bundledlrtrustregionprox_supports_dense_variable_metric(const BundleDLRTrustRegionProx* self) {
  return self->supports_dense_variable_metric();
}

dll int cb_bundledlrtrustregionprox_supports_lowrank_variable_metric(const BundleDLRTrustRegionProx* self) {
  return self->supports_lowrank_variable_metric();
}

dll int cb_bundledlrtrustregionprox_supports_diagonal_variable_metric(const BundleDLRTrustRegionProx* self) {
  return self->supports_diagonal_variable_metric();
}

dll int cb_bundledlrtrustregionprox_apply_variable_metric(BundleDLRTrustRegionProx* self, VariableMetricModel* groundset, VariableMetricModel* model, const Matrix* aggr, Integer y_id, const Matrix* y, int descent_step, Real* current_weight, Real model_maxviol, const Indexmatrix* new_indices = 0) {
  return self->apply_variable_metric(groundset, model, *aggr, y_id, *y, (bool)descent_step, *current_weight, model_maxviol, new_indices);
}

dll int cb_bundledlrtrustregionprox_add_variable_metric(BundleDLRTrustRegionProx* self, Matrix* diagH, Matrix* vecH) {
  return self->add_variable_metric(*diagH, *vecH);
}

dll int cb_bundledlrtrustregionprox_push_aft(BundleDLRTrustRegionProx* self, const AffineFunctionTransformation* aft) {
  return self->push_aft(aft);
}

dll int cb_bundledlrtrustregionprox_pop_aft(BundleDLRTrustRegionProx* self) {
  return self->pop_aft();
}

dll int cb_bundledlrtrustregionprox_mfile_data(const BundleDLRTrustRegionProx* self) {
  return self->mfile_data(std::cout);
}

