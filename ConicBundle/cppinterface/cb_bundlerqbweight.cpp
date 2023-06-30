dll void cb_bundlerqbweight_destroy(BundleRQBWeight* self) {
  delete self;
}

dll BundleRQBWeight* cb_bundlerqbweight_new(BundleWeight* bwp = 0, int incr = -1) {
  return new BundleRQBWeight(bwp, 0, incr);
}

dll BundleRQBWeight* cb_bundlerqbweight_new2(Real m1, Real m2 = .2, Real m3 = 1., Real eta = 1e-6, BundleWeight* bwp = 0, int incr = -1) {
  return new BundleRQBWeight(m1, m2, m3, eta, bwp, 0, incr);
}

dll void cb_bundlerqbweight_set_defaults(BundleRQBWeight* self) {
  self->set_defaults();
}

dll void cb_bundlerqbweight_clear(BundleRQBWeight* self) {
  self->clear();
}

dll int cb_bundlerqbweight_init(BundleRQBWeight* self, Real aggr_dnmormsqr, Groundset* groundset, BundleModel* model) {
  return self->init(aggr_dnmormsqr, groundset, model);
}

dll void cb_bundlerqbweight_set_next_weight(BundleRQBWeight* self, Real u) {
  self->set_next_weight(u);
}

dll void cb_bundlerqbweight_set_minweight(BundleRQBWeight* self, Real mw) {
  self->set_minweight(mw);
}

dll int cb_bundlerqbweight_get_next_weight_set(const BundleRQBWeight* self) {
  return self->get_next_weight_set();
}

dll Real cb_bundlerqbweight_get_minweight(const BundleRQBWeight* self) {
  return self->get_minweight();
}

dll void cb_bundlerqbweight_set_maxweight(BundleRQBWeight* self, Real mw) {
  self->set_maxweight(mw);
}

dll Real cb_bundlerqbweight_get_maxweight(const BundleRQBWeight* self) {
  return self->get_maxweight();
}

dll Real cb_bundlerqbweight_get_weight(const BundleRQBWeight* self) {
  return self->get_weight();
}

dll int cb_bundlerqbweight_weight_changed(const BundleRQBWeight* self) {
  return self->weight_changed();
}

dll int cb_bundlerqbweight_descent_update(BundleRQBWeight* self, Real newval, Real oldval, Real modelval, const Matrix* y, const Matrix* newy, Real normsubg2, BundleProxObject* Hp) {
  return self->descent_update(newval, oldval, modelval, *y, *newy, normsubg2, Hp);
}

dll int cb_bundlerqbweight_nullstep_update(BundleRQBWeight* self, Real newval, Real oldval, Real modelval, const Matrix* y, const Matrix* newy, MinorantPointer* new_minorant, MinorantPointer* aggregate, Real nullstep_bound, Real normsubg2, BundleProxObject* Hp) {
  return self->nullstep_update(newval, oldval, modelval, *y, *newy, *new_minorant, *aggregate, nullstep_bound, normsubg2, Hp);
}

dll int cb_bundlerqbweight_apply_modification(BundleRQBWeight* self, const GroundsetModification* gsmdf) {
  return self->apply_modification(*gsmdf);
}

