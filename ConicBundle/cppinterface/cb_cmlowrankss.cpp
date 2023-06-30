dll void cb_cmlowrankss_destroy(CMlowrankss* self) {
  delete self;
}

dll CMlowrankss* cb_cmlowrankss_new(const Sparsemat* Ain, const Sparsemat* Bin, CoeffmatInfo* cip = 0) {
  return new CMlowrankss(*Ain, *Bin, cip);
}

dll Coeffmat* cb_cmlowrankss_clone(const CMlowrankss* self) {
  return self->clone();
}

dll Integer cb_cmlowrankss_dim(const CMlowrankss* self) {
  return self->dim();
}

dll Real cb_cmlowrankss_get(const CMlowrankss* self, Integer i, Integer j) {
  return (*self)(i, j);
}

dll void cb_cmlowrankss_make_symmatrix(const CMlowrankss* self, Symmatrix* S) {
  self->make_symmatrix(*S);
}

dll Real cb_cmlowrankss_norm(const CMlowrankss* self) {
  return self->norm();
}

dll Coeffmat* cb_cmlowrankss_subspace(const CMlowrankss* self, const Matrix* P) {
  return self->subspace(*P);
}

dll void cb_cmlowrankss_multiply(CMlowrankss* self, Real d) {
  self->multiply(d);
}

dll Real cb_cmlowrankss_ip(const CMlowrankss* self, const Symmatrix* S) {
  return self->ip(*S);
}

dll Real cb_cmlowrankss_gramip(const CMlowrankss* self, const Matrix* P) {
  return self->gramip(*P);
}

dll Real cb_cmlowrankss_gramip2(const CMlowrankss* self, const Matrix* P, Integer start_row, const Matrix* Lam = 0) {
  return self->gramip(*P, start_row, Lam);
}

dll void cb_cmlowrankss_addmeto(const CMlowrankss* self, Symmatrix* S, Real d = 1.) {
  self->addmeto(*S, d);
}

dll void cb_cmlowrankss_addprodto(const CMlowrankss* self, Matrix* D, const Matrix* C, Real d = 1.) {
  self->addprodto(*D, *C, d);
}

dll void cb_cmlowrankss_addprodto2(const CMlowrankss* self, Matrix* D, const Sparsemat* C, Real d = 1.) {
  self->addprodto(*D, *C, d);
}

dll void cb_cmlowrankss_left_right_prod(const CMlowrankss* self, const Matrix* P, const Matrix* Q, Matrix* R) {
  self->left_right_prod(*P, *Q, *R);
}

dll Integer cb_cmlowrankss_prodvec_flops(const CMlowrankss* self) {
  return self->prodvec_flops();
}

dll int cb_cmlowrankss_dense(const CMlowrankss* self) {
  return self->dense();
}

dll int cb_cmlowrankss_sparse(const CMlowrankss* self) {
  return self->sparse();
}

dll int cb_cmlowrankss_sparse2(const CMlowrankss* self, Indexmatrix* param0, Indexmatrix* param1, Matrix* param2, Real param3) {
  return self->sparse(*param0, *param1, *param2, param3);
}

dll int cb_cmlowrankss_support_in(const CMlowrankss* self, const Sparsesym* param0) {
  return self->support_in(*param0);
}

dll Real cb_cmlowrankss_ip2(const CMlowrankss* self, const Sparsesym* S) {
  return self->ip(*S);
}

dll void cb_cmlowrankss_project(const CMlowrankss* self, Symmatrix* S, const Matrix* P) {
  self->project(*S, *P);
}

dll void cb_cmlowrankss_add_projection(const CMlowrankss* self, Symmatrix* S, const Matrix* P, Real alpha = 1., Integer start_row = 0) {
  self->add_projection(*S, *P, alpha, start_row);
}

dll const Matrix* cb_cmlowrankss_postgenmult(const CMlowrankss* self, const Matrix* D, Matrix* C, Real alpha = 1., Real beta = 0., int dtrans = 0) {
  return &self->postgenmult(*D, *C, alpha, beta, dtrans);
}

dll const Matrix* cb_cmlowrankss_pregenmult(const CMlowrankss* self, const Matrix* D, Matrix* C, Real alpha = 1., Real beta = 0., int dtrans = 0) {
  return &self->pregenmult(*D, *C, alpha, beta, dtrans);
}

dll int cb_cmlowrankss_equal(const CMlowrankss* self, const Coeffmat* p, double tol = 1e-6) {
  return self->equal(p, tol);
}

dll void cb_cmlowrankss_display(const CMlowrankss* self) {
  self->display(std::cout);
}

dll void cb_cmlowrankss_out(const CMlowrankss* self) {
  self->out(std::cout);
}

