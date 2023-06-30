dll void cb_cmsingleton_destroy(CMsingleton* self) {
  delete self;
}

dll CMsingleton* cb_cmsingleton_new(Integer innr, Integer ini, Integer inj, Real inval, CoeffmatInfo* cip = 0) {
  return new CMsingleton(innr, ini, inj, inval, cip);
}

dll Coeffmat* cb_cmsingleton_clone(const CMsingleton* self) {
  return self->clone();
}

dll Integer cb_cmsingleton_dim(const CMsingleton* self) {
  return self->dim();
}

dll Real cb_cmsingleton_get(const CMsingleton* self, Integer i, Integer j) {
  return (*self)(i, j);
}

dll void cb_cmsingleton_make_symmatrix(const CMsingleton* self, Symmatrix* S) {
  self->make_symmatrix(*S);
}

dll Real cb_cmsingleton_norm(const CMsingleton* self) {
  return self->norm();
}

dll Coeffmat* cb_cmsingleton_subspace(const CMsingleton* self, const Matrix* P) {
  return self->subspace(*P);
}

dll void cb_cmsingleton_multiply(CMsingleton* self, Real d) {
  self->multiply(d);
}

dll Real cb_cmsingleton_ip(const CMsingleton* self, const Symmatrix* S) {
  return self->ip(*S);
}

dll Real cb_cmsingleton_gramip(const CMsingleton* self, const Matrix* P) {
  return self->gramip(*P);
}

dll Real cb_cmsingleton_gramip2(const CMsingleton* self, const Matrix* P, Integer start_row, const Matrix* Lam = 0) {
  return self->gramip(*P, start_row, Lam);
}

dll void cb_cmsingleton_addmeto(const CMsingleton* self, Symmatrix* S, Real d = 1.) {
  self->addmeto(*S, d);
}

dll void cb_cmsingleton_addprodto(const CMsingleton* self, Matrix* B, const Matrix* C, Real d = 1.) {
  self->addprodto(*B, *C, d);
}

dll void cb_cmsingleton_addprodto2(const CMsingleton* self, Matrix* B, const Sparsemat* C, Real d = 1.) {
  self->addprodto(*B, *C, d);
}

dll void cb_cmsingleton_left_right_prod(const CMsingleton* self, const Matrix* P, const Matrix* Q, Matrix* R) {
  self->left_right_prod(*P, *Q, *R);
}

dll Integer cb_cmsingleton_prodvec_flops(const CMsingleton* self) {
  return self->prodvec_flops();
}

dll int cb_cmsingleton_dense(const CMsingleton* self) {
  return self->dense();
}

dll int cb_cmsingleton_sparse(const CMsingleton* self) {
  return self->sparse();
}

dll int cb_cmsingleton_sparse2(const CMsingleton* self, Indexmatrix* I, Indexmatrix* J, Matrix* v, Real d = 1.) {
  return self->sparse(*I, *J, *v, d);
}

dll int cb_cmsingleton_support_in(const CMsingleton* self, const Sparsesym* S) {
  return self->support_in(*S);
}

dll Real cb_cmsingleton_ip2(const CMsingleton* self, const Sparsesym* S) {
  return self->ip(*S);
}

dll void cb_cmsingleton_project(const CMsingleton* self, Symmatrix* S, const Matrix* P) {
  self->project(*S, *P);
}

dll void cb_cmsingleton_add_projection(const CMsingleton* self, Symmatrix* S, const Matrix* P, Real alpha = 1., Integer start_row = 0) {
  self->add_projection(*S, *P, alpha, start_row);
}

dll const Matrix* cb_cmsingleton_postgenmult(const CMsingleton* self, const Matrix* B, Matrix* C, Real alpha = 1., Real beta = 0., int btrans = 0) {
  return &self->postgenmult(*B, *C, alpha, beta, btrans);
}

dll const Matrix* cb_cmsingleton_pregenmult(const CMsingleton* self, const Matrix* B, Matrix* C, Real alpha = 1., Real beta = 0., int btrans = 0) {
  return &self->pregenmult(*B, *C, alpha, beta, btrans);
}

dll int cb_cmsingleton_equal(const CMsingleton* self, const Coeffmat* p, double tol = 1e-6) {
  return self->equal(p, tol);
}

dll void cb_cmsingleton_display(const CMsingleton* self) {
  self->display(std::cout);
}

dll void cb_cmsingleton_out(const CMsingleton* self) {
  self->out(std::cout);
}

dll int cb_cmsingleton_get_ijval(const CMsingleton* self, Integer* i, Integer* j, Real* v) {
  return self->get_ijval(*i, *j, *v);
}

