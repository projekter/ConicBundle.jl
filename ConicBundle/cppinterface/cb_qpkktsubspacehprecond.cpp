dll void cb_qpkktsubspacehprecond_destroy(QPKKTSubspaceHPrecond* self) {
  delete self;
}

dll QPKKTSubspaceHPrecond* cb_qpkktsubspacehprecond_new(Integer inmethod = 0, int cbinc = -1) {
  return new QPKKTSubspaceHPrecond(inmethod, 0, cbinc);
}

dll int cb_qpkktsubspacehprecond_init_data(QPKKTSubspaceHPrecond* self, QPSolverProxObject* Hp, QPModelBlockObject* model, const Sparsemat* A, const Indexmatrix* eq_indices, int SchurComplAineq) {
  return self->init_data(Hp, model, A, eq_indices, (bool)SchurComplAineq);
}

dll int cb_qpkktsubspacehprecond_init_system(QPKKTSubspaceHPrecond* self, const Matrix* KKTdiagx, const Matrix* KKTdiagy, Real Hfactor, Real prec, QPSolverParameters* params) {
  return self->init_system(*KKTdiagx, *KKTdiagy, Hfactor, prec, params);
}

dll int cb_qpkktsubspacehprecond_precondm1(QPKKTSubspaceHPrecond* self, Matrix* vec) {
  return self->precondM1(*vec);
}

dll int cb_qpkktsubspacehprecond_set_subspace(QPKKTSubspaceHPrecond* self, const Matrix* insubspace) {
  return self->set_subspace(*insubspace);
}

dll Real cb_qpkktsubspacehprecond_get_lmin_invm1(QPKKTSubspaceHPrecond* self) {
  return self->get_lmin_invM1();
}

dll int cb_qpkktsubspacehprecond_precond_invg1(QPKKTSubspaceHPrecond* self, Matrix* vec) {
  return self->precond_invG1(*vec);
}

dll int cb_qpkktsubspacehprecond_precond_invg1tran(QPKKTSubspaceHPrecond* self, Matrix* vec) {
  return self->precond_invG1tran(*vec);
}

dll Integer cb_qpkktsubspacehprecond_precond_size(QPKKTSubspaceHPrecond* self) {
  return self->precond_size();
}

dll int cb_qpkktsubspacehprecond_cond_number_mult(QPKKTSubspaceHPrecond* self, Matrix* vec, const Matrix* KKTdiagx, const Matrix* KKTdiagy) {
  return self->cond_number_mult(*vec, *KKTdiagx, *KKTdiagy);
}

dll Integer cb_qpkktsubspacehprecond_get_precond_rank(QPKKTSubspaceHPrecond* self) {
  return self->get_precond_rank();
}

dll CH_Tools::Microseconds* cb_qpkktsubspacehprecond_new_get_t_precond_mult(QPKKTSubspaceHPrecond* self) {
  return new CH_Tools::Microseconds(self->get_t_precond_mult());
}

dll void cb_qpkktsubspacehprecond_reset_t_precond_mult(QPKKTSubspaceHPrecond* self) {
  self->reset_t_precond_mult();
}

