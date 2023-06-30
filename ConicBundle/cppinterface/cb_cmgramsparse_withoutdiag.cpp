dll void cb_cmgramsparse_withoutdiag_destroy(CMgramsparse_withoutdiag* self) {
  delete self;
}

dll CMgramsparse_withoutdiag* cb_cmgramsparse_withoutdiag_new(const Sparsemat* Ain, int pos = 1, CoeffmatInfo* cip = 0) {
  return new CMgramsparse_withoutdiag(*Ain, (bool)pos, cip);
}

dll Coeffmat* cb_cmgramsparse_withoutdiag_clone(const CMgramsparse_withoutdiag* self) {
  return self->clone();
}

dll Integer cb_cmgramsparse_withoutdiag_dim(const CMgramsparse_withoutdiag* self) {
  return self->dim();
}

dll Real cb_cmgramsparse_withoutdiag_get(const CMgramsparse_withoutdiag* self, Integer i, Integer j) {
  return (*self)(i, j);
}

dll void cb_cmgramsparse_withoutdiag_make_symmatrix(const CMgramsparse_withoutdiag* self, Symmatrix* S) {
  self->make_symmatrix(*S);
}

dll Real cb_cmgramsparse_withoutdiag_norm(const CMgramsparse_withoutdiag* self) {
  return self->norm();
}

dll Coeffmat* cb_cmgramsparse_withoutdiag_subspace(const CMgramsparse_withoutdiag* self, const Matrix* P) {
  return self->subspace(*P);
}

dll void cb_cmgramsparse_withoutdiag_multiply(CMgramsparse_withoutdiag* self, Real d) {
  self->multiply(d);
}

dll Real cb_cmgramsparse_withoutdiag_ip(const CMgramsparse_withoutdiag* self, const Symmatrix* S) {
  return self->ip(*S);
}

dll Real cb_cmgramsparse_withoutdiag_gramip(const CMgramsparse_withoutdiag* self, const Matrix* P) {
  return self->gramip(*P);
}

dll Real cb_cmgramsparse_withoutdiag_gramip2(const CMgramsparse_withoutdiag* self, const Matrix* P, Integer start_row, const Matrix* Lam = 0) {
  return self->gramip(*P, start_row, Lam);
}

dll void cb_cmgramsparse_withoutdiag_addmeto(const CMgramsparse_withoutdiag* self, Symmatrix* S, Real d = 1.) {
  self->addmeto(*S, d);
}

dll void cb_cmgramsparse_withoutdiag_addprodto(const CMgramsparse_withoutdiag* self, Matrix* B, const Matrix* C, Real d = 1.) {
  self->addprodto(*B, *C, d);
}

dll void cb_cmgramsparse_withoutdiag_addprodto2(const CMgramsparse_withoutdiag* self, Matrix* B, const Sparsemat* C, Real d = 1.) {
  self->addprodto(*B, *C, d);
}

dll void cb_cmgramsparse_withoutdiag_left_right_prod(const CMgramsparse_withoutdiag* self, const Matrix* P, const Matrix* Q, Matrix* R) {
  self->left_right_prod(*P, *Q, *R);
}

dll Integer cb_cmgramsparse_withoutdiag_prodvec_flops(const CMgramsparse_withoutdiag* self) {
  return self->prodvec_flops();
}

dll int cb_cmgramsparse_withoutdiag_dense(const CMgramsparse_withoutdiag* self) {
  return self->dense();
}

dll int cb_cmgramsparse_withoutdiag_sparse(const CMgramsparse_withoutdiag* self) {
  return self->sparse();
}

dll int cb_cmgramsparse_withoutdiag_sparse2(const CMgramsparse_withoutdiag* self, Indexmatrix* param0, Indexmatrix* param1, Matrix* param2, Real param3) {
  return self->sparse(*param0, *param1, *param2, param3);
}

dll int cb_cmgramsparse_withoutdiag_support_in(const CMgramsparse_withoutdiag* self, const Sparsesym* param0) {
  return self->support_in(*param0);
}

dll Real cb_cmgramsparse_withoutdiag_ip2(const CMgramsparse_withoutdiag* self, const Sparsesym* S) {
  return self->ip(*S);
}

dll void cb_cmgramsparse_withoutdiag_project(const CMgramsparse_withoutdiag* self, Symmatrix* S, const Matrix* P) {
  self->project(*S, *P);
}

dll void cb_cmgramsparse_withoutdiag_add_projection(const CMgramsparse_withoutdiag* self, Symmatrix* S, const Matrix* P, Real alpha = 1., Integer start_row = 0) {
  self->add_projection(*S, *P, alpha, start_row);
}

dll const Matrix* cb_cmgramsparse_withoutdiag_postgenmult(const CMgramsparse_withoutdiag* self, const Matrix* B, Matrix* C, Real alpha = 1., Real beta = 0., int btrans = 0) {
  return &self->postgenmult(*B, *C, alpha, beta, btrans);
}

dll const Matrix* cb_cmgramsparse_withoutdiag_pregenmult(const CMgramsparse_withoutdiag* self, const Matrix* B, Matrix* C, Real alpha = 1., Real beta = 0., int btrans = 0) {
  return &self->pregenmult(*B, *C, alpha, beta, btrans);
}

dll int cb_cmgramsparse_withoutdiag_equal(const CMgramsparse_withoutdiag* self, const Coeffmat* p, double tol = 1e-6) {
  return self->equal(p, tol);
}

dll void cb_cmgramsparse_withoutdiag_display(const CMgramsparse_withoutdiag* self) {
  self->display(std::cout);
}

dll void cb_cmgramsparse_withoutdiag_out(const CMgramsparse_withoutdiag* self) {
  self->out(std::cout);
}

dll const Sparsemat* cb_cmgramsparse_withoutdiag_get_a(const CMgramsparse_withoutdiag* self) {
  return &self->get_A();
}

dll int cb_cmgramsparse_withoutdiag_get_positive(const CMgramsparse_withoutdiag* self) {
  return self->get_positive();
}

