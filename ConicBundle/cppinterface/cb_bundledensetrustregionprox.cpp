dll void cb_bundledensetrustregionprox_destroy(BundleDenseTrustRegionProx* self) {
  delete self;
}

dll BundleDenseTrustRegionProx* cb_bundledensetrustregionprox_new(const Symmatrix* Hin, VariableMetricSelection* vp = 0, int local_metric = 0, int bounds_scaling = 0, int cbinc = -1) {
  return new BundleDenseTrustRegionProx(*Hin, vp, (bool)local_metric, (bool)bounds_scaling, 0, cbinc);
}

dll BundleDenseTrustRegionProx* cb_bundledensetrustregionprox_new2(Integer dim = 0, VariableMetricSelection* vp = 0, int local_metric = 0, int bounds_scaling = 0, int cbinc = -1) {
  return new BundleDenseTrustRegionProx(dim, vp, (bool)local_metric, (bool)bounds_scaling, 0, cbinc);
}

dll void cb_bundledensetrustregionprox_set_weightu(BundleDenseTrustRegionProx* self, Real in_weightu) {
  self->set_weightu(in_weightu);
}

dll Real cb_bundledensetrustregionprox_get_weightu(const BundleDenseTrustRegionProx* self) {
  return self->get_weightu();
}

dll Real cb_bundledensetrustregionprox_get_term_corr(const BundleDenseTrustRegionProx* self) {
  return self->get_term_corr();
}

dll const Symmatrix* cb_bundledensetrustregionprox_init(BundleDenseTrustRegionProx* self, Symmatrix* in_H) {
  return &self->init(*in_H);
}

dll int cb_bundledensetrustregionprox_get_factored(const BundleDenseTrustRegionProx* self) {
  return self->get_factored();
}

dll const Symmatrix* cb_bundledensetrustregionprox_get_h(const BundleDenseTrustRegionProx* self) {
  return &self->get_H();
}

dll const Symmatrix* cb_bundledensetrustregionprox_get_hchol(const BundleDenseTrustRegionProx* self) {
  return &self->get_Hchol();
}

dll Integer cb_bundledensetrustregionprox_dim(const BundleDenseTrustRegionProx* self) {
  return self->dim();
}

dll const Real cb_bundledensetrustregionprox_get(const BundleDenseTrustRegionProx* self, Integer i, Integer j) {
  return (*self)(i, j);
}

dll Real cb_bundledensetrustregionprox_norm_sqr(const BundleDenseTrustRegionProx* self, const Matrix* B) {
  return self->norm_sqr(*B);
}

dll Real cb_bundledensetrustregionprox_dnorm_sqr(const BundleDenseTrustRegionProx* self, const MinorantPointer* B) {
  return self->dnorm_sqr(*B);
}

dll int cb_bundledensetrustregionprox_is_dlr(const BundleDenseTrustRegionProx* self) {
  return self->is_DLR();
}

dll int cb_bundledensetrustregionprox_add_h(const BundleDenseTrustRegionProx* self, Symmatrix* big_sym, Integer start_index = 0) {
  return self->add_H(*big_sym, start_index);
}

dll Matrix* cb_bundledensetrustregionprox_add_hx(const BundleDenseTrustRegionProx* self, const Matrix* x, Matrix* outplusHx, Real alpha = 1.) {
  return &self->add_Hx(*x, *outplusHx, alpha);
}

dll Matrix* cb_bundledensetrustregionprox_apply_hinv(const BundleDenseTrustRegionProx* self, Matrix* x) {
  return &self->apply_Hinv(*x);
}

dll int cb_bundledensetrustregionprox_compute_qp_costs(BundleDenseTrustRegionProx* self, Symmatrix* Q, Matrix* d, Real* offset, const MinorantPointer* constant_minorant, const MinorantBundle* bundle, const Matrix* y, const MinorantPointer* groundset_minorant, Indexmatrix* yfixed) {
  return self->compute_QP_costs(*Q, *d, *offset, *constant_minorant, *bundle, *y, *groundset_minorant, yfixed);
}

dll int cb_bundledensetrustregionprox_update_qp_costs(BundleDenseTrustRegionProx* self, Symmatrix* delta_Q, Matrix* delta_d, Real* delta_offset, const MinorantPointer* constant_minorant, const MinorantBundle* bundle, const Matrix* center_y, const MinorantPointer* groundset_minorant, const MinorantPointer* delta_groundset_minorant, const Indexmatrix* delta_index, Indexmatrix* yfixed) {
  return self->update_QP_costs(*delta_Q, *delta_d, *delta_offset, *constant_minorant, *bundle, *center_y, *groundset_minorant, *delta_groundset_minorant, *delta_index, yfixed);
}

dll int cb_bundledensetrustregionprox_apply_modification(BundleDenseTrustRegionProx* self, const GroundsetModification* gsmdf) {
  return self->apply_modification(*gsmdf);
}

dll BundleProxObject* cb_bundledensetrustregionprox_projected_clone(BundleDenseTrustRegionProx* self, const Indexmatrix* indices) {
  return self->projected_clone(*indices);
}

dll int cb_bundledensetrustregionprox_supports_diagonal_bounds_scaling(const BundleDenseTrustRegionProx* self) {
  return self->supports_diagonal_bounds_scaling();
}

dll int cb_bundledensetrustregionprox_diagonal_bounds_scaling_update(BundleDenseTrustRegionProx* self, const Matrix* param0) {
  return self->diagonal_bounds_scaling_update(*param0);
}

dll int cb_bundledensetrustregionprox_supports_dense_variable_metric(const BundleDenseTrustRegionProx* self) {
  return self->supports_dense_variable_metric();
}

dll int cb_bundledensetrustregionprox_supports_lowrank_variable_metric(const BundleDenseTrustRegionProx* self) {
  return self->supports_lowrank_variable_metric();
}

dll int cb_bundledensetrustregionprox_supports_diagonal_variable_metric(const BundleDenseTrustRegionProx* self) {
  return self->supports_diagonal_variable_metric();
}

dll int cb_bundledensetrustregionprox_apply_variable_metric(BundleDenseTrustRegionProx* self, VariableMetricModel* groundset, VariableMetricModel* model, const Matrix* aggr, Integer y_id, const Matrix* y, int descent_step, Real* current_weight, Real model_maxviol, const Indexmatrix* new_indices = 0) {
  return self->apply_variable_metric(groundset, model, *aggr, y_id, *y, (bool)descent_step, *current_weight, model_maxviol, new_indices);
}

dll int cb_bundledensetrustregionprox_add_variable_metric(BundleDenseTrustRegionProx* self, Symmatrix* addH) {
  return self->add_variable_metric(*addH);
}

dll int cb_bundledensetrustregionprox_add_variable_metric2(BundleDenseTrustRegionProx* self, Matrix* diagH, Matrix* vecH) {
  return self->add_variable_metric(*diagH, *vecH);
}

dll int cb_bundledensetrustregionprox_push_aft(BundleDenseTrustRegionProx* self, const AffineFunctionTransformation* aft) {
  return self->push_aft(aft);
}

dll int cb_bundledensetrustregionprox_pop_aft(BundleDenseTrustRegionProx* self) {
  return self->pop_aft();
}

dll int cb_bundledensetrustregionprox_mfile_data(const BundleDenseTrustRegionProx* self) {
  return self->mfile_data(std::cout);
}

