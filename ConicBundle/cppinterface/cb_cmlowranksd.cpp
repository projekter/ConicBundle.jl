dll void cb_cmlowranksd_destroy(CMlowranksd* self) {
  delete self;
}

dll CMlowranksd* cb_cmlowranksd_new(const Sparsemat* Ain, const Matrix* Bin, CoeffmatInfo* cip = 0) {
  return new CMlowranksd(*Ain, *Bin, cip);
}

dll Coeffmat* cb_cmlowranksd_clone(const CMlowranksd* self) {
  return self->clone();
}

dll Integer cb_cmlowranksd_dim(const CMlowranksd* self) {
  return self->dim();
}

dll Real cb_cmlowranksd_get(const CMlowranksd* self, Integer i, Integer j) {
  return (*self)(i, j);
}

dll void cb_cmlowranksd_make_symmatrix(const CMlowranksd* self, Symmatrix* S) {
  self->make_symmatrix(*S);
}

dll Real cb_cmlowranksd_norm(const CMlowranksd* self) {
  return self->norm();
}

dll Coeffmat* cb_cmlowranksd_subspace(const CMlowranksd* self, const Matrix* P) {
  return self->subspace(*P);
}

dll void cb_cmlowranksd_multiply(CMlowranksd* self, Real d) {
  self->multiply(d);
}

dll Real cb_cmlowranksd_ip(const CMlowranksd* self, const Symmatrix* S) {
  return self->ip(*S);
}

dll Real cb_cmlowranksd_gramip(const CMlowranksd* self, const Matrix* P) {
  return self->gramip(*P);
}

dll Real cb_cmlowranksd_gramip2(const CMlowranksd* self, const Matrix* P, Integer start_row, const Matrix* Lam = 0) {
  return self->gramip(*P, start_row, Lam);
}

dll void cb_cmlowranksd_addmeto(const CMlowranksd* self, Symmatrix* S, Real d = 1.) {
  self->addmeto(*S, d);
}

dll void cb_cmlowranksd_addprodto(const CMlowranksd* self, Matrix* D, const Matrix* C, Real d = 1.) {
  self->addprodto(*D, *C, d);
}

dll void cb_cmlowranksd_addprodto2(const CMlowranksd* self, Matrix* D, const Sparsemat* C, Real d = 1.) {
  self->addprodto(*D, *C, d);
}

dll void cb_cmlowranksd_left_right_prod(const CMlowranksd* self, const Matrix* P, const Matrix* Q, Matrix* R) {
  self->left_right_prod(*P, *Q, *R);
}

dll Integer cb_cmlowranksd_prodvec_flops(const CMlowranksd* self) {
  return self->prodvec_flops();
}

dll int cb_cmlowranksd_dense(const CMlowranksd* self) {
  return self->dense();
}

dll int cb_cmlowranksd_sparse(const CMlowranksd* self) {
  return self->sparse();
}

dll int cb_cmlowranksd_sparse2(const CMlowranksd* self, Indexmatrix* param0, Indexmatrix* param1, Matrix* param2, Real param3) {
  return self->sparse(*param0, *param1, *param2, param3);
}

dll int cb_cmlowranksd_support_in(const CMlowranksd* self, const Sparsesym* param0) {
  return self->support_in(*param0);
}

dll Real cb_cmlowranksd_ip2(const CMlowranksd* self, const Sparsesym* S) {
  return self->ip(*S);
}

dll void cb_cmlowranksd_project(const CMlowranksd* self, Symmatrix* S, const Matrix* P) {
  self->project(*S, *P);
}

dll void cb_cmlowranksd_add_projection(const CMlowranksd* self, Symmatrix* S, const Matrix* P, Real alpha = 1., Integer start_row = 0) {
  self->add_projection(*S, *P, alpha, start_row);
}

dll const Matrix* cb_cmlowranksd_postgenmult(const CMlowranksd* self, const Matrix* D, Matrix* C, Real alpha = 1., Real beta = 0., int dtrans = 0) {
  return &self->postgenmult(*D, *C, alpha, beta, dtrans);
}

dll const Matrix* cb_cmlowranksd_pregenmult(const CMlowranksd* self, const Matrix* D, Matrix* C, Real alpha = 1., Real beta = 0., int dtrans = 0) {
  return &self->pregenmult(*D, *C, alpha, beta, dtrans);
}

dll int cb_cmlowranksd_equal(const CMlowranksd* self, const Coeffmat* p, double tol = 1e-6) {
  return self->equal(p, tol);
}

dll void cb_cmlowranksd_display(const CMlowranksd* self) {
  self->display(std::cout);
}

dll void cb_cmlowranksd_out(const CMlowranksd* self) {
  self->out(std::cout);
}

