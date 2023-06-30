dll void cb_cmgramsparse_destroy(CMgramsparse* self) {
  delete self;
}

dll CMgramsparse* cb_cmgramsparse_new(const Sparsemat* Ain, int pos = 1, CoeffmatInfo* cip = 0) {
  return new CMgramsparse(*Ain, (bool)pos, cip);
}

dll Coeffmat* cb_cmgramsparse_clone(const CMgramsparse* self) {
  return self->clone();
}

dll Integer cb_cmgramsparse_dim(const CMgramsparse* self) {
  return self->dim();
}

dll Real cb_cmgramsparse_get(const CMgramsparse* self, Integer i, Integer j) {
  return (*self)(i, j);
}

dll void cb_cmgramsparse_make_symmatrix(const CMgramsparse* self, Symmatrix* S) {
  self->make_symmatrix(*S);
}

dll Real cb_cmgramsparse_norm(const CMgramsparse* self) {
  return self->norm();
}

dll Coeffmat* cb_cmgramsparse_subspace(const CMgramsparse* self, const Matrix* P) {
  return self->subspace(*P);
}

dll void cb_cmgramsparse_multiply(CMgramsparse* self, Real d) {
  self->multiply(d);
}

dll Real cb_cmgramsparse_ip(const CMgramsparse* self, const Symmatrix* S) {
  return self->ip(*S);
}

dll Real cb_cmgramsparse_gramip(const CMgramsparse* self, const Matrix* P) {
  return self->gramip(*P);
}

dll Real cb_cmgramsparse_gramip2(const CMgramsparse* self, const Matrix* P, Integer start_row, const Matrix* Lam = 0) {
  return self->gramip(*P, start_row, Lam);
}

dll void cb_cmgramsparse_addmeto(const CMgramsparse* self, Symmatrix* S, Real d = 1.) {
  self->addmeto(*S, d);
}

dll void cb_cmgramsparse_addprodto(const CMgramsparse* self, Matrix* B, const Matrix* C, Real d = 1.) {
  self->addprodto(*B, *C, d);
}

dll void cb_cmgramsparse_addprodto2(const CMgramsparse* self, Matrix* B, const Sparsemat* C, Real d = 1.) {
  self->addprodto(*B, *C, d);
}

dll void cb_cmgramsparse_left_right_prod(const CMgramsparse* self, const Matrix* P, const Matrix* Q, Matrix* R) {
  self->left_right_prod(*P, *Q, *R);
}

dll Integer cb_cmgramsparse_prodvec_flops(const CMgramsparse* self) {
  return self->prodvec_flops();
}

dll int cb_cmgramsparse_dense(const CMgramsparse* self) {
  return self->dense();
}

dll int cb_cmgramsparse_sparse(const CMgramsparse* self) {
  return self->sparse();
}

dll int cb_cmgramsparse_sparse2(const CMgramsparse* self, Indexmatrix* param0, Indexmatrix* param1, Matrix* param2, Real param3) {
  return self->sparse(*param0, *param1, *param2, param3);
}

dll int cb_cmgramsparse_support_in(const CMgramsparse* self, const Sparsesym* param0) {
  return self->support_in(*param0);
}

dll Real cb_cmgramsparse_ip2(const CMgramsparse* self, const Sparsesym* S) {
  return self->ip(*S);
}

dll void cb_cmgramsparse_project(const CMgramsparse* self, Symmatrix* S, const Matrix* P) {
  self->project(*S, *P);
}

dll void cb_cmgramsparse_add_projection(const CMgramsparse* self, Symmatrix* S, const Matrix* P, Real alpha = 1., Integer start_row = 0) {
  self->add_projection(*S, *P, alpha, start_row);
}

dll const Matrix* cb_cmgramsparse_postgenmult(const CMgramsparse* self, const Matrix* B, Matrix* C, Real alpha = 1., Real beta = 0., int btrans = 0) {
  return &self->postgenmult(*B, *C, alpha, beta, btrans);
}

dll const Matrix* cb_cmgramsparse_pregenmult(const CMgramsparse* self, const Matrix* B, Matrix* C, Real alpha = 1., Real beta = 0., int btrans = 0) {
  return &self->pregenmult(*B, *C, alpha, beta, btrans);
}

dll int cb_cmgramsparse_equal(const CMgramsparse* self, const Coeffmat* p, double tol = 1e-6) {
  return self->equal(p, tol);
}

dll void cb_cmgramsparse_display(const CMgramsparse* self) {
  self->display(std::cout);
}

dll void cb_cmgramsparse_out(const CMgramsparse* self) {
  self->out(std::cout);
}

dll const Sparsemat* cb_cmgramsparse_get_a(const CMgramsparse* self) {
  return &self->get_A();
}

dll int cb_cmgramsparse_get_positive(const CMgramsparse* self) {
  return self->get_positive();
}

