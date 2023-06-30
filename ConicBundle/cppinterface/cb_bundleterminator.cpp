dll void cb_bundleterminator_destroy(BundleTerminator* self) {
  delete self;
}

dll void cb_bundleterminator_set_defaults(BundleTerminator* self) {
  self->set_defaults();
}

dll void cb_bundleterminator_clear(BundleTerminator* self) {
  self->clear();
}

dll BundleTerminator* cb_bundleterminator_new(int incr = -1) {
  return new BundleTerminator(0, incr);
}

dll void cb_bundleterminator_set_termeps(BundleTerminator* self, Real teps) {
  self->set_termeps(teps);
}

dll Real cb_bundleterminator_get_termeps(const BundleTerminator* self) {
  return self->get_termeps();
}

dll void cb_bundleterminator_set_aggr_dnormsqr(BundleTerminator* self, Real sg) {
  self->set_aggr_dnormsqr(sg);
}

dll Real cb_bundleterminator_get_aggr_dnormsqr(const BundleTerminator* self) {
  return self->get_aggr_dnormsqr();
}

dll void cb_bundleterminator_set_timelimit(BundleTerminator* self, const CH_Tools::Clock* cp, CH_Tools::Microseconds* tl) {
  self->set_timelimit(cp, *tl);
}

dll CH_Tools::Microseconds* cb_bundleterminator_new_get_timelimit(const BundleTerminator* self) {
  return new CH_Tools::Microseconds(self->get_timelimit());
}

dll void cb_bundleterminator_set_recomplimit(BundleTerminator* self, Integer rl) {
  self->set_recomplimit(rl);
}

dll Integer cb_bundleterminator_get_recomplimit(const BundleTerminator* self) {
  return self->get_recomplimit();
}

dll void cb_bundleterminator_set_qpfailslimit(BundleTerminator* self, Integer ql) {
  self->set_qpfailslimit(ql);
}

dll Integer cb_bundleterminator_get_qpfailslimit(const BundleTerminator* self) {
  return self->get_qpfailslimit();
}

dll void cb_bundleterminator_set_modelfailslimit(BundleTerminator* self, Integer ml) {
  self->set_modelfailslimit(ml);
}

dll Integer cb_bundleterminator_get_modelfailslimit(const BundleTerminator* self) {
  return self->get_modelfailslimit();
}

dll void cb_bundleterminator_set_augvalfailslimit(BundleTerminator* self, Integer al) {
  self->set_augvalfailslimit(al);
}

dll Integer cb_bundleterminator_get_augvalfailslimit(const BundleTerminator* self) {
  return self->get_augvalfailslimit();
}

dll void cb_bundleterminator_set_objevallimit(BundleTerminator* self, Integer ol) {
  self->set_objevallimit(ol);
}

dll Integer cb_bundleterminator_get_objevallimit(const BundleTerminator* self) {
  return self->get_objevallimit();
}

dll void cb_bundleterminator_set_oraclefailslimit(BundleTerminator* self, Integer ol) {
  self->set_oraclefailslimit(ol);
}

dll Integer cb_bundleterminator_get_oraclefailslimit(const BundleTerminator* self) {
  return self->get_oraclefailslimit();
}

dll int cb_bundleterminator_get_terminated(const BundleTerminator* self) {
  return self->get_terminated();
}

dll void cb_bundleterminator_clear_terminated(BundleTerminator* self) {
  self->clear_terminated();
}

dll void cb_bundleterminator_print_status(const BundleTerminator* self) {
  self->print_status(std::cout);
}

dll void cb_bundleterminator_save(const BundleTerminator* self) {
  self->save(std::cout);
}

