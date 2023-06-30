dll void cb_pcg_destroy(PCG* self) {
  delete self;
}

dll PCG* cb_pcg_new(int pril = -1) {
  return new PCG(&std::cout, pril);
}

dll void cb_pcg_set_maxit(PCG* self, Integer in_maxit) {
  self->set_maxit(in_maxit);
}

dll Integer cb_pcg_get_maxit(const PCG* self) {
  return self->get_maxit();
}

dll int cb_pcg_get_err(const PCG* self) {
  return self->get_err();
}

dll Integer cb_pcg_get_nmult(const PCG* self) {
  return self->get_nmult();
}

dll Real cb_pcg_get_residual_norm(const PCG* self) {
  return self->get_residual_norm();
}

dll Real cb_pcg_get_avg_reduction(const PCG* self) {
  return self->get_avg_reduction();
}

dll Real cb_pcg_get_termprec(const PCG* self) {
  return self->get_termprec();
}

dll int cb_pcg_compute(PCG* self, IterativeSystemObject* system, Matrix* x, Real termprec, Matrix* storex = 0, Integer storestep = 0) {
  return self->compute(*system, *x, termprec, storex, storestep);
}

dll void cb_pcg_set_out(PCG* self, int in_print_level = 1) {
  self->set_out(&std::cout, in_print_level);
}

