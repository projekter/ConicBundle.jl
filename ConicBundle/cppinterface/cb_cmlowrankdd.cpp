dll void cb_cmlowrankdd_destroy(CMlowrankdd* self) {
  delete self;
}

dll CMlowrankdd* cb_cmlowrankdd_new(const Matrix* Ain, const Matrix* Bin, CoeffmatInfo* cip = 0) {
  return new CMlowrankdd(*Ain, *Bin, cip);
}

dll Coeffmat* cb_cmlowrankdd_clone(const CMlowrankdd* self) {
  return self->clone();
}

dll Integer cb_cmlowrankdd_dim(const CMlowrankdd* self) {
  return self->dim();
}

dll Real cb_cmlowrankdd_get(const CMlowrankdd* self, Integer i, Integer j) {
  return (*self)(i, j);
}

dll void cb_cmlowrankdd_make_symmatrix(const CMlowrankdd* self, Symmatrix* S) {
  self->make_symmatrix(*S);
}

dll Real cb_cmlowrankdd_norm(const CMlowrankdd* self) {
  return self->norm();
}

dll Coeffmat* cb_cmlowrankdd_subspace(const CMlowrankdd* self, const Matrix* P) {
  return self->subspace(*P);
}

dll void cb_cmlowrankdd_multiply(CMlowrankdd* self, Real d) {
  self->multiply(d);
}

dll Real cb_cmlowrankdd_ip(const CMlowrankdd* self, const Symmatrix* S) {
  return self->ip(*S);
}

dll Real cb_cmlowrankdd_gramip(const CMlowrankdd* self, const Matrix* P) {
  return self->gramip(*P);
}

dll Real cb_cmlowrankdd_gramip2(const CMlowrankdd* self, const Matrix* P, Integer start_row, const Matrix* Lam = 0) {
  return self->gramip(*P, start_row, Lam);
}

dll void cb_cmlowrankdd_addmeto(const CMlowrankdd* self, Symmatrix* S, Real d = 1.) {
  self->addmeto(*S, d);
}

dll void cb_cmlowrankdd_addprodto(const CMlowrankdd* self, Matrix* D, const Matrix* C, Real d = 1.) {
  self->addprodto(*D, *C, d);
}

dll void cb_cmlowrankdd_addprodto2(const CMlowrankdd* self, Matrix* D, const Sparsemat* C, Real d = 1.) {
  self->addprodto(*D, *C, d);
}

dll void cb_cmlowrankdd_left_right_prod(const CMlowrankdd* self, const Matrix* P, const Matrix* Q, Matrix* R) {
  self->left_right_prod(*P, *Q, *R);
}

dll Integer cb_cmlowrankdd_prodvec_flops(const CMlowrankdd* self) {
  return self->prodvec_flops();
}

dll int cb_cmlowrankdd_dense(const CMlowrankdd* self) {
  return self->dense();
}

dll int cb_cmlowrankdd_sparse(const CMlowrankdd* self) {
  return self->sparse();
}

dll int cb_cmlowrankdd_sparse2(const CMlowrankdd* self, Indexmatrix* param0, Indexmatrix* param1, Matrix* param2, Real param3) {
  return self->sparse(*param0, *param1, *param2, param3);
}

dll int cb_cmlowrankdd_support_in(const CMlowrankdd* self, const Sparsesym* param0) {
  return self->support_in(*param0);
}

dll Real cb_cmlowrankdd_ip2(const CMlowrankdd* self, const Sparsesym* S) {
  return self->ip(*S);
}

dll void cb_cmlowrankdd_project(const CMlowrankdd* self, Symmatrix* S, const Matrix* P) {
  self->project(*S, *P);
}

dll void cb_cmlowrankdd_add_projection(const CMlowrankdd* self, Symmatrix* S, const Matrix* P, Real alpha = 1., Integer start_row = 0) {
  self->add_projection(*S, *P, alpha, start_row);
}

dll const Matrix* cb_cmlowrankdd_postgenmult(const CMlowrankdd* self, const Matrix* D, Matrix* C, Real alpha = 1., Real beta = 0., int dtrans = 0) {
  return &self->postgenmult(*D, *C, alpha, beta, dtrans);
}

dll const Matrix* cb_cmlowrankdd_pregenmult(const CMlowrankdd* self, const Matrix* D, Matrix* C, Real alpha = 1., Real beta = 0., int dtrans = 0) {
  return &self->pregenmult(*D, *C, alpha, beta, dtrans);
}

dll int cb_cmlowrankdd_equal(const CMlowrankdd* self, const Coeffmat* p, double tol = 1e-6) {
  return self->equal(p, tol);
}

dll void cb_cmlowrankdd_display(const CMlowrankdd* self) {
  self->display(std::cout);
}

dll void cb_cmlowrankdd_out(const CMlowrankdd* self) {
  self->out(std::cout);
}

