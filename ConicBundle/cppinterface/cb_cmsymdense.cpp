dll void cb_cmsymdense_destroy(CMsymdense* self) {
  delete self;
}

dll CMsymdense* cb_cmsymdense_new(const Symmatrix* Ain, CoeffmatInfo* cip = 0) {
  return new CMsymdense(*Ain, cip);
}

dll Coeffmat* cb_cmsymdense_clone(const CMsymdense* self) {
  return self->clone();
}

dll Integer cb_cmsymdense_dim(const CMsymdense* self) {
  return self->dim();
}

dll Real cb_cmsymdense_get(const CMsymdense* self, Integer i, Integer j) {
  return (*self)(i, j);
}

dll void cb_cmsymdense_make_symmatrix(const CMsymdense* self, Symmatrix* S) {
  self->make_symmatrix(*S);
}

dll Real cb_cmsymdense_norm(const CMsymdense* self) {
  return self->norm();
}

dll Coeffmat* cb_cmsymdense_subspace(const CMsymdense* self, const Matrix* P) {
  return self->subspace(*P);
}

dll void cb_cmsymdense_multiply(CMsymdense* self, Real d) {
  self->multiply(d);
}

dll Real cb_cmsymdense_ip(const CMsymdense* self, const Symmatrix* S) {
  return self->ip(*S);
}

dll Real cb_cmsymdense_gramip(const CMsymdense* self, const Matrix* P) {
  return self->gramip(*P);
}

dll Real cb_cmsymdense_gramip2(const CMsymdense* self, const Matrix* P, Integer start_row, const Matrix* Lam = 0) {
  return self->gramip(*P, start_row, Lam);
}

dll void cb_cmsymdense_addmeto(const CMsymdense* self, Symmatrix* S, Real d = 1.) {
  self->addmeto(*S, d);
}

dll void cb_cmsymdense_addprodto(const CMsymdense* self, Matrix* B, const Matrix* C, Real d = 1.) {
  self->addprodto(*B, *C, d);
}

dll void cb_cmsymdense_addprodto2(const CMsymdense* self, Matrix* B, const Sparsemat* C, Real d = 1.) {
  self->addprodto(*B, *C, d);
}

dll void cb_cmsymdense_left_right_prod(const CMsymdense* self, const Matrix* P, const Matrix* Q, Matrix* R) {
  self->left_right_prod(*P, *Q, *R);
}

dll Integer cb_cmsymdense_prodvec_flops(const CMsymdense* self) {
  return self->prodvec_flops();
}

dll int cb_cmsymdense_dense(const CMsymdense* self) {
  return self->dense();
}

dll int cb_cmsymdense_sparse(const CMsymdense* self) {
  return self->sparse();
}

dll int cb_cmsymdense_sparse2(const CMsymdense* self, Indexmatrix* param0, Indexmatrix* param1, Matrix* param2, Real param3) {
  return self->sparse(*param0, *param1, *param2, param3);
}

dll int cb_cmsymdense_support_in(const CMsymdense* self, const Sparsesym* param0) {
  return self->support_in(*param0);
}

dll Real cb_cmsymdense_ip2(const CMsymdense* self, const Sparsesym* S) {
  return self->ip(*S);
}

dll void cb_cmsymdense_project(const CMsymdense* self, Symmatrix* S, const Matrix* P) {
  self->project(*S, *P);
}

dll void cb_cmsymdense_add_projection(const CMsymdense* self, Symmatrix* S, const Matrix* P, Real alpha = 1., Integer start_row = 0) {
  self->add_projection(*S, *P, alpha, start_row);
}

dll const Matrix* cb_cmsymdense_postgenmult(const CMsymdense* self, const Matrix* B, Matrix* C, Real alpha = 1., Real beta = 0., int btrans = 0) {
  return &self->postgenmult(*B, *C, alpha, beta, btrans);
}

dll const Matrix* cb_cmsymdense_pregenmult(const CMsymdense* self, const Matrix* B, Matrix* C, Real alpha = 1., Real beta = 0., int btrans = 0) {
  return &self->pregenmult(*B, *C, alpha, beta, btrans);
}

dll int cb_cmsymdense_equal(const CMsymdense* self, const Coeffmat* p, double tol = 1e-6) {
  return self->equal(p, tol);
}

dll void cb_cmsymdense_display(const CMsymdense* self) {
  self->display(std::cout);
}

dll void cb_cmsymdense_out(const CMsymdense* self) {
  self->out(std::cout);
}

dll const Symmatrix* cb_cmsymdense_get_a(const CMsymdense* self) {
  return &self->get_A();
}

