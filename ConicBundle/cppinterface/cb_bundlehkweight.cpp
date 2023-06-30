dll void cb_bundlehkweight_destroy(BundleHKWeight* self) {
  delete self;
}

dll BundleHKWeight* cb_bundlehkweight_new(Real mRin = .5, BundleWeight* bwp = 0, int incr = -1) {
  return new BundleHKWeight(mRin, bwp, 0, incr);
}

dll void cb_bundlehkweight_set_defaults(BundleHKWeight* self) {
  self->set_defaults();
}

dll void cb_bundlehkweight_set_nullstep_updates(BundleHKWeight* self, int nu = 0) {
  self->set_nullstep_updates(nu);
}

dll void cb_bundlehkweight_clear(BundleHKWeight* self) {
  self->clear();
}

dll int cb_bundlehkweight_init(BundleHKWeight* self, Real aggr_dnmormsqr, Groundset* groundset, BundleModel* model) {
  return self->init(aggr_dnmormsqr, groundset, model);
}

dll void cb_bundlehkweight_set_next_weight(BundleHKWeight* self, Real u) {
  self->set_next_weight(u);
}

dll void cb_bundlehkweight_set_minweight(BundleHKWeight* self, Real mw) {
  self->set_minweight(mw);
}

dll int cb_bundlehkweight_get_next_weight_set(const BundleHKWeight* self) {
  return self->get_next_weight_set();
}

dll Real cb_bundlehkweight_get_minweight(const BundleHKWeight* self) {
  return self->get_minweight();
}

dll void cb_bundlehkweight_set_maxweight(BundleHKWeight* self, Real mw) {
  self->set_maxweight(mw);
}

dll Real cb_bundlehkweight_get_maxweight(const BundleHKWeight* self) {
  return self->get_maxweight();
}

dll Real cb_bundlehkweight_get_weight(const BundleHKWeight* self) {
  return self->get_weight();
}

dll int cb_bundlehkweight_weight_changed(const BundleHKWeight* self) {
  return self->weight_changed();
}

dll int cb_bundlehkweight_descent_update(BundleHKWeight* self, Real newval, Real oldval, Real modelval, const Matrix* y, const Matrix* newy, Real normsubg2, BundleProxObject* Hp) {
  return self->descent_update(newval, oldval, modelval, *y, *newy, normsubg2, Hp);
}

dll int cb_bundlehkweight_nullstep_update(BundleHKWeight* self, Real newval, Real oldval, Real modelval, const Matrix* y, const Matrix* newy, MinorantPointer* new_minorant, MinorantPointer* aggregate, Real nullstep_bound, Real normsubg2, BundleProxObject* Hp) {
  return self->nullstep_update(newval, oldval, modelval, *y, *newy, *new_minorant, *aggregate, nullstep_bound, normsubg2, Hp);
}

dll int cb_bundlehkweight_apply_modification(BundleHKWeight* self, const GroundsetModification* gsmdf) {
  return self->apply_modification(*gsmdf);
}

