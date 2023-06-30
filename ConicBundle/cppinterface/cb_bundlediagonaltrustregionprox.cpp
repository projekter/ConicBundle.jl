dll void cb_bundlediagonaltrustregionprox_destroy(BundleDiagonalTrustRegionProx* self) {
  delete self;
}

dll BundleDiagonalTrustRegionProx* cb_bundlediagonaltrustregionprox_new(const Matrix* Din, VariableMetricSelection* vp = 0, int local_scaling = 0, int bounds_scaling = 0, int cbinc = -1) {
  return new BundleDiagonalTrustRegionProx(*Din, vp, (bool)local_scaling, (bool)bounds_scaling, 0, cbinc);
}

dll BundleDiagonalTrustRegionProx* cb_bundlediagonaltrustregionprox_new2(Integer dim, Real d, VariableMetricSelection* vp = 0, int local_scaling = 0, int bounds_scaling = 0, int cbinc = -1) {
  return new BundleDiagonalTrustRegionProx(dim, d, vp, (bool)local_scaling, (bool)bounds_scaling, 0, cbinc);
}

dll BundleDiagonalTrustRegionProx* cb_bundlediagonaltrustregionprox_new3(Integer dim = 0, VariableMetricSelection* vp = 0, int local_scaling = 0, int bounds_scaling = 0, int cbinc = -1) {
  return new BundleDiagonalTrustRegionProx(dim, vp, (bool)local_scaling, (bool)bounds_scaling, 0, cbinc);
}

dll void cb_bundlediagonaltrustregionprox_set_weightu(BundleDiagonalTrustRegionProx* self, Real in_weightu) {
  self->set_weightu(in_weightu);
}

dll Real cb_bundlediagonaltrustregionprox_get_weightu(const BundleDiagonalTrustRegionProx* self) {
  return self->get_weightu();
}

dll Real cb_bundlediagonaltrustregionprox_get_term_corr(const BundleDiagonalTrustRegionProx* self) {
  return self->get_term_corr();
}

dll const Matrix* cb_bundlediagonaltrustregionprox_get_d(const BundleDiagonalTrustRegionProx* self) {
  return &self->get_D();
}

dll void cb_bundlediagonaltrustregionprox_set_d(BundleDiagonalTrustRegionProx* self, Matrix* in_D) {
  self->set_D(*in_D);
}

dll Integer cb_bundlediagonaltrustregionprox_dim(const BundleDiagonalTrustRegionProx* self) {
  return self->dim();
}

dll Real cb_bundlediagonaltrustregionprox_get(const BundleDiagonalTrustRegionProx* self, Integer i) {
  return (*self)(i);
}

dll Real cb_bundlediagonaltrustregionprox_norm_sqr(const BundleDiagonalTrustRegionProx* self, const Matrix* B) {
  return self->norm_sqr(*B);
}

dll Real cb_bundlediagonaltrustregionprox_dnorm_sqr(const BundleDiagonalTrustRegionProx* self, const MinorantPointer* B) {
  return self->dnorm_sqr(*B);
}

dll int cb_bundlediagonaltrustregionprox_is_dlr(const BundleDiagonalTrustRegionProx* self) {
  return self->is_DLR();
}

dll int cb_bundlediagonaltrustregionprox_add_h(const BundleDiagonalTrustRegionProx* self, Symmatrix* big_sym, Integer start_index = 0) {
  return self->add_H(*big_sym, start_index);
}

dll Matrix* cb_bundlediagonaltrustregionprox_add_hx(const BundleDiagonalTrustRegionProx* self, const Matrix* x, Matrix* outplusHx, Real alpha = 1.) {
  return &self->add_Hx(*x, *outplusHx, alpha);
}

dll Matrix* cb_bundlediagonaltrustregionprox_apply_hinv(const BundleDiagonalTrustRegionProx* self, Matrix* x) {
  return &self->apply_Hinv(*x);
}

dll int cb_bundlediagonaltrustregionprox_compute_qp_costs(BundleDiagonalTrustRegionProx* self, Symmatrix* Q, Matrix* d, Real* offset, const MinorantPointer* constant_minorant, const MinorantBundle* bundle, const Matrix* y, const MinorantPointer* groundset_minorant, Indexmatrix* yfixed) {
  return self->compute_QP_costs(*Q, *d, *offset, *constant_minorant, *bundle, *y, *groundset_minorant, yfixed);
}

dll int cb_bundlediagonaltrustregionprox_update_qp_costs(BundleDiagonalTrustRegionProx* self, Symmatrix* delta_Q, Matrix* delta_d, Real* delta_offset, const MinorantPointer* constant_minorant, const MinorantBundle* bundle, const Matrix* center_y, const MinorantPointer* groundset_minorant, const MinorantPointer* delta_groundset_minorant, const Indexmatrix* delta_index, Indexmatrix* yfixed) {
  return self->update_QP_costs(*delta_Q, *delta_d, *delta_offset, *constant_minorant, *bundle, *center_y, *groundset_minorant, *delta_groundset_minorant, *delta_index, yfixed);
}

dll int cb_bundlediagonaltrustregionprox_apply_modification(BundleDiagonalTrustRegionProx* self, const GroundsetModification* gsmdf) {
  return self->apply_modification(*gsmdf);
}

dll BundleProxObject* cb_bundlediagonaltrustregionprox_projected_clone(BundleDiagonalTrustRegionProx* self, const Indexmatrix* indices) {
  return self->projected_clone(*indices);
}

dll int cb_bundlediagonaltrustregionprox_supports_diagonal_bounds_scaling(const BundleDiagonalTrustRegionProx* self) {
  return self->supports_diagonal_bounds_scaling();
}

dll int cb_bundlediagonaltrustregionprox_diagonal_bounds_scaling_update(BundleDiagonalTrustRegionProx* self, const Matrix* D_update) {
  return self->diagonal_bounds_scaling_update(*D_update);
}

dll int cb_bundlediagonaltrustregionprox_supports_lowrank_variable_metric(const BundleDiagonalTrustRegionProx* self) {
  return self->supports_lowrank_variable_metric();
}

dll int cb_bundlediagonaltrustregionprox_supports_diagonal_variable_metric(const BundleDiagonalTrustRegionProx* self) {
  return self->supports_diagonal_variable_metric();
}

dll int cb_bundlediagonaltrustregionprox_apply_variable_metric(BundleDiagonalTrustRegionProx* self, VariableMetricModel* groundset, VariableMetricModel* model, const Matrix* aggr, Integer y_id, const Matrix* y, int descent_step, Real* current_weight, Real model_maxviol, const Indexmatrix* new_indices = 0) {
  return self->apply_variable_metric(groundset, model, *aggr, y_id, *y, (bool)descent_step, *current_weight, model_maxviol, new_indices);
}

dll int cb_bundlediagonaltrustregionprox_add_variable_metric(BundleDiagonalTrustRegionProx* self, Matrix* diagH, Matrix* vecH) {
  return self->add_variable_metric(*diagH, *vecH);
}

dll int cb_bundlediagonaltrustregionprox_push_aft(BundleDiagonalTrustRegionProx* self, const AffineFunctionTransformation* aft) {
  return self->push_aft(aft);
}

dll int cb_bundlediagonaltrustregionprox_pop_aft(BundleDiagonalTrustRegionProx* self) {
  return self->pop_aft();
}

dll int cb_bundlediagonaltrustregionprox_mfile_data(const BundleDiagonalTrustRegionProx* self) {
  return self->mfile_data(std::cout);
}

