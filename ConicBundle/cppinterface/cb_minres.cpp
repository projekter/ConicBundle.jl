dll void cb_minres_destroy(MinRes* self) {
  delete self;
}

dll MinRes* cb_minres_new(int pril = -1) {
  return new MinRes(&std::cout, pril);
}

dll void cb_minres_set_maxit(MinRes* self, Integer in_maxit) {
  self->set_maxit(in_maxit);
}

dll Integer cb_minres_get_maxit(const MinRes* self) {
  return self->get_maxit();
}

dll int cb_minres_get_err(const MinRes* self) {
  return self->get_err();
}

dll Integer cb_minres_get_nmult(const MinRes* self) {
  return self->get_nmult();
}

dll Real cb_minres_get_residual_norm(const MinRes* self) {
  return self->get_residual_norm();
}

dll Real cb_minres_get_avg_reduction(const MinRes* self) {
  return self->get_avg_reduction();
}

dll Real cb_minres_get_termprec(const MinRes* self) {
  return self->get_termprec();
}

dll int cb_minres_compute(MinRes* self, IterativeSystemObject* system, Matrix* x, Real termprec, Matrix* storex = 0, Integer storestep = 0) {
  return self->compute(*system, *x, termprec, storex, storestep);
}

dll void cb_minres_set_out(MinRes* self, int in_print_level = 1) {
  self->set_out(&std::cout, in_print_level);
}

