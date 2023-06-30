dll void cb_bundleidprox_destroy(BundleIdProx* self) {
  delete self;
}

dll BundleIdProx* cb_bundleidprox_new(Integer d = 0, Real w = 1., int cbinc = -1) {
  return new BundleIdProx(d, w, 0, cbinc);
}

dll void cb_bundleidprox_set_weightu(BundleIdProx* self, Real in_weightu) {
  self->set_weightu(in_weightu);
}

dll Real cb_bundleidprox_get_weightu(const BundleIdProx* self) {
  return self->get_weightu();
}

dll Real cb_bundleidprox_get_term_corr(const BundleIdProx* self) {
  return self->get_term_corr();
}

dll Real cb_bundleidprox_norm_sqr(const BundleIdProx* self, const Matrix* B) {
  return self->norm_sqr(*B);
}

dll Real cb_bundleidprox_dnorm_sqr(const BundleIdProx* self, const MinorantPointer* B) {
  return self->dnorm_sqr(*B);
}

dll int cb_bundleidprox_is_dlr(const BundleIdProx* self) {
  return self->is_DLR();
}

dll int cb_bundleidprox_add_h(const BundleIdProx* self, Symmatrix* big_sym, Integer start_index = 0) {
  return self->add_H(*big_sym, start_index);
}

dll Matrix* cb_bundleidprox_add_hx(const BundleIdProx* self, const Matrix* x, Matrix* outplusHx, Real alpha = 1.) {
  return &self->add_Hx(*x, *outplusHx, alpha);
}

dll Matrix* cb_bundleidprox_apply_hinv(const BundleIdProx* self, Matrix* x) {
  return &self->apply_Hinv(*x);
}

dll int cb_bundleidprox_compute_qp_costs(BundleIdProx* self, Symmatrix* Q, Matrix* d, Real* offset, const MinorantPointer* constant_minorant, const MinorantBundle* bundle, const Matrix* y, const MinorantPointer* groundset_minorant, Indexmatrix* yfixed) {
  return self->compute_QP_costs(*Q, *d, *offset, *constant_minorant, *bundle, *y, *groundset_minorant, yfixed);
}

dll int cb_bundleidprox_update_qp_costs(BundleIdProx* self, Symmatrix* delta_Q, Matrix* delta_d, Real* delta_offset, const MinorantPointer* constant_minorant, const MinorantBundle* bundle, const Matrix* center_y, const MinorantPointer* groundset_minorant, const MinorantPointer* delta_groundset_minorant, const Indexmatrix* delta_index, Indexmatrix* yfixed) {
  return self->update_QP_costs(*delta_Q, *delta_d, *delta_offset, *constant_minorant, *bundle, *center_y, *groundset_minorant, *delta_groundset_minorant, *delta_index, yfixed);
}

dll int cb_bundleidprox_apply_modification(BundleIdProx* self, const GroundsetModification* gsmdf) {
  return self->apply_modification(*gsmdf);
}

dll BundleProxObject* cb_bundleidprox_projected_clone(BundleIdProx* self, const Indexmatrix* indices) {
  return self->projected_clone(*indices);
}

dll int cb_bundleidprox_mfile_data(const BundleIdProx* self) {
  return self->mfile_data(std::cout);
}

