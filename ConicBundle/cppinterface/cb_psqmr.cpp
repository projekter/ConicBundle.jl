dll void cb_psqmr_destroy(Psqmr* self) {
  delete self;
}

dll Psqmr* cb_psqmr_new(int pril = -1) {
  return new Psqmr(&std::cout, pril);
}

dll void cb_psqmr_set_maxit(Psqmr* self, Integer in_maxit) {
  self->set_maxit(in_maxit);
}

dll Integer cb_psqmr_get_maxit(const Psqmr* self) {
  return self->get_maxit();
}

dll int cb_psqmr_get_err(const Psqmr* self) {
  return self->get_err();
}

dll Integer cb_psqmr_get_nmult(const Psqmr* self) {
  return self->get_nmult();
}

dll Real cb_psqmr_get_residual_norm(const Psqmr* self) {
  return self->get_residual_norm();
}

dll Real cb_psqmr_get_avg_reduction(const Psqmr* self) {
  return self->get_avg_reduction();
}

dll Real cb_psqmr_get_termprec(const Psqmr* self) {
  return self->get_termprec();
}

dll int cb_psqmr_compute(Psqmr* self, IterativeSystemObject* system, Matrix* x, Real termprec, Matrix* storex = 0, Integer storestep = 0) {
  return self->compute(*system, *x, termprec, storex, storestep);
}

dll void cb_psqmr_set_out(Psqmr* self, int in_print_level = 1) {
  self->set_out(&std::cout, in_print_level);
}

