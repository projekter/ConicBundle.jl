dll void cb_cmsymsparse_destroy(CMsymsparse* self) {
  delete self;
}

dll CMsymsparse* cb_cmsymsparse_new(const Sparsesym* Ain, CoeffmatInfo* cip = 0) {
  return new CMsymsparse(*Ain, cip);
}

dll Coeffmat* cb_cmsymsparse_clone(const CMsymsparse* self) {
  return self->clone();
}

dll Integer cb_cmsymsparse_dim(const CMsymsparse* self) {
  return self->dim();
}

dll Real cb_cmsymsparse_get(const CMsymsparse* self, Integer i, Integer j) {
  return (*self)(i, j);
}

dll void cb_cmsymsparse_make_symmatrix(const CMsymsparse* self, Symmatrix* S) {
  self->make_symmatrix(*S);
}

dll Real cb_cmsymsparse_norm(const CMsymsparse* self) {
  return self->norm();
}

dll Coeffmat* cb_cmsymsparse_subspace(const CMsymsparse* self, const Matrix* P) {
  return self->subspace(*P);
}

dll void cb_cmsymsparse_multiply(CMsymsparse* self, Real d) {
  self->multiply(d);
}

dll Real cb_cmsymsparse_ip(const CMsymsparse* self, const Symmatrix* S) {
  return self->ip(*S);
}

dll Real cb_cmsymsparse_gramip(const CMsymsparse* self, const Matrix* P) {
  return self->gramip(*P);
}

dll Real cb_cmsymsparse_gramip2(const CMsymsparse* self, const Matrix* P, Integer start_row, const Matrix* Lam = 0) {
  return self->gramip(*P, start_row, Lam);
}

dll void cb_cmsymsparse_addmeto(const CMsymsparse* self, Symmatrix* S, Real d = 1.) {
  self->addmeto(*S, d);
}

dll void cb_cmsymsparse_addprodto(const CMsymsparse* self, Matrix* B, const Matrix* C, Real d = 1.) {
  self->addprodto(*B, *C, d);
}

dll void cb_cmsymsparse_addprodto2(const CMsymsparse* self, Matrix* B, const Sparsemat* C, Real d = 1.) {
  self->addprodto(*B, *C, d);
}

dll void cb_cmsymsparse_left_right_prod(const CMsymsparse* self, const Matrix* P, const Matrix* Q, Matrix* R) {
  self->left_right_prod(*P, *Q, *R);
}

dll Integer cb_cmsymsparse_prodvec_flops(const CMsymsparse* self) {
  return self->prodvec_flops();
}

dll int cb_cmsymsparse_dense(const CMsymsparse* self) {
  return self->dense();
}

dll int cb_cmsymsparse_sparse(const CMsymsparse* self) {
  return self->sparse();
}

dll int cb_cmsymsparse_sparse2(const CMsymsparse* self, Indexmatrix* I, Indexmatrix* J, Matrix* val, Real d = 1.) {
  return self->sparse(*I, *J, *val, d);
}

dll int cb_cmsymsparse_support_in(const CMsymsparse* self, const Sparsesym* S) {
  return self->support_in(*S);
}

dll Real cb_cmsymsparse_ip2(const CMsymsparse* self, const Sparsesym* S) {
  return self->ip(*S);
}

dll void cb_cmsymsparse_project(const CMsymsparse* self, Symmatrix* S, const Matrix* P) {
  self->project(*S, *P);
}

dll void cb_cmsymsparse_add_projection(const CMsymsparse* self, Symmatrix* S, const Matrix* P, Real alpha = 1., Integer start_row = 0) {
  self->add_projection(*S, *P, alpha, start_row);
}

dll const Matrix* cb_cmsymsparse_postgenmult(const CMsymsparse* self, const Matrix* B, Matrix* C, Real alpha = 1., Real beta = 0., int btrans = 0) {
  return &self->postgenmult(*B, *C, alpha, beta, btrans);
}

dll const Matrix* cb_cmsymsparse_pregenmult(const CMsymsparse* self, const Matrix* B, Matrix* C, Real alpha = 1., Real beta = 0., int btrans = 0) {
  return &self->pregenmult(*B, *C, alpha, beta, btrans);
}

dll int cb_cmsymsparse_equal(const CMsymsparse* self, const Coeffmat* p, double tol = 1e-6) {
  return self->equal(p, tol);
}

dll void cb_cmsymsparse_display(const CMsymsparse* self) {
  self->display(std::cout);
}

dll void cb_cmsymsparse_out(const CMsymsparse* self) {
  self->out(std::cout);
}

dll const Sparsesym* cb_cmsymsparse_get_a(const CMsymsparse* self) {
  return &self->get_A();
}

