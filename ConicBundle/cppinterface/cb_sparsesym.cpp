dll void cb_sparsesym_destroy(Sparsesym* self) {
  delete self;
}

dll Sparsesym* cb_sparsesym_new() {
  return new Sparsesym();
}

dll Sparsesym* cb_sparsesym_new2(const Sparsesym* A, Real d = 1.) {
  return new Sparsesym(*A, d);
}

dll Sparsesym* cb_sparsesym_new3(Integer nr) {
  return new Sparsesym(nr);
}

dll Sparsesym* cb_sparsesym_new4(Integer nr, Integer nz, const Integer* ini, const Integer* inj, const Real* va) {
  return new Sparsesym(nr, nz, ini, inj, va);
}

dll Sparsesym* cb_sparsesym_new5(Integer nr, Integer nz, const Indexmatrix* ini, const Indexmatrix* inj, const Matrix* va) {
  return new Sparsesym(nr, nz, *ini, *inj, *va);
}

dll Sparsesym* cb_sparsesym_init(Sparsesym* self, const Sparsesym* param0, Real d = 1.) {
  return &self->init(*param0, d);
}

dll Sparsesym* cb_sparsesym_init2(Sparsesym* self, const Matrix* param0, Real d = 1.) {
  return &self->init(*param0, d);
}

dll Sparsesym* cb_sparsesym_init3(Sparsesym* self, const Indexmatrix* param0, Real d = 1.) {
  return &self->init(*param0, d);
}

dll Sparsesym* cb_sparsesym_init4(Sparsesym* self, const Symmatrix* param0, Real d = 1.) {
  return &self->init(*param0, d);
}

dll Sparsesym* cb_sparsesym_init5(Sparsesym* self, const Sparsemat* param0, Real d = 1.) {
  return &self->init(*param0, d);
}

dll Sparsesym* cb_sparsesym_init6(Sparsesym* self, Integer nr) {
  return &self->init(nr);
}

dll Sparsesym* cb_sparsesym_init7(Sparsesym* self, Integer nr, Integer nz, const Integer* ini, const Integer* inj, const Real* va) {
  return &self->init(nr, nz, ini, inj, va);
}

dll Sparsesym* cb_sparsesym_init8(Sparsesym* self, Integer nr, Integer nz, const Indexmatrix* ini, const Indexmatrix* inj, const Matrix* va) {
  return &self->init(nr, nz, *ini, *inj, *va);
}

dll Sparsesym* cb_sparsesym_init_support(Sparsesym* self, const Sparsesym* A, Real d = 0.) {
  return &self->init_support(*A, d);
}

dll void cb_sparsesym_set_tol(Sparsesym* self, Real t) {
  self->set_tol(t);
}

dll Sparsesym* cb_sparsesym_new6(const Matrix* param0, Real d = 1.) {
  return new Sparsesym(*param0, d);
}

dll Sparsesym* cb_sparsesym_new7(const Indexmatrix* param0, Real d = 1.) {
  return new Sparsesym(*param0, d);
}

dll Sparsesym* cb_sparsesym_new8(const Symmatrix* param0, Real d = 1.) {
  return new Sparsesym(*param0, d);
}

dll Sparsesym* cb_sparsesym_new9(const Sparsemat* param0, Real d = 1.) {
  return new Sparsesym(*param0, d);
}

dll void cb_sparsesym_dim(const Sparsesym* self, Integer* r, Integer* c) {
  self->dim(*r, *c);
}

dll Integer cb_sparsesym_dim2(const Sparsesym* self) {
  return self->dim();
}

dll Integer cb_sparsesym_rowdim(const Sparsesym* self) {
  return self->rowdim();
}

dll Integer cb_sparsesym_coldim(const Sparsesym* self) {
  return self->coldim();
}

dll Integer cb_sparsesym_nonzeros(const Sparsesym* self) {
  return self->nonzeros();
}

dll Real cb_sparsesym_get(const Sparsesym* self, Integer i, Integer j) {
  return (*self)(i, j);
}

dll Real cb_sparsesym_get2(const Sparsesym* self, Integer i) {
  return (*self)(i);
}

dll const Indexmatrix* cb_sparsesym_get_colinfo(const Sparsesym* self) {
  return &self->get_colinfo();
}

dll const Indexmatrix* cb_sparsesym_get_colindex(const Sparsesym* self) {
  return &self->get_colindex();
}

dll const Matrix* cb_sparsesym_get_colval(const Sparsesym* self) {
  return &self->get_colval();
}

dll const Indexmatrix* cb_sparsesym_get_suppind(const Sparsesym* self) {
  return &self->get_suppind();
}

dll const Indexmatrix* cb_sparsesym_get_suppcol(const Sparsesym* self) {
  return &self->get_suppcol();
}

dll void cb_sparsesym_get_edge_rep(const Sparsesym* self, Indexmatrix* I, Indexmatrix* J, Matrix* val) {
  self->get_edge_rep(*I, *J, *val);
}

dll int cb_sparsesym_contains_support(const Sparsesym* self, const Sparsesym* A) {
  return self->contains_support(*A);
}

dll int cb_sparsesym_check_support(const Sparsesym* self, Integer i, Integer j) {
  return self->check_support(i, j);
}

dll Matrix* cb_sparsesym_new_diag(const Sparsesym* A) {
  return new Matrix(diag(*A));
}

dll Sparsesym* cb_sparsesym_new_sparsediag(const Matrix* A, Real tol) {
  return new Sparsesym(sparseDiag(*A, tol));
}

dll void cb_sparsesym_swap(Sparsesym* A, Sparsesym* B) {
  swap(*A, *B);
}

dll Sparsesym* cb_sparsesym_xeya(Sparsesym* self, const Sparsesym* A, Real d = 1.) {
  return &self->xeya(*A, d);
}

dll Sparsesym* cb_sparsesym_xeya2(Sparsesym* self, const Matrix* A, Real d = 1.) {
  return &self->xeya(*A, d);
}

dll Sparsesym* cb_sparsesym_xeya3(Sparsesym* self, const Indexmatrix* A, Real d = 1.) {
  return &self->xeya(*A, d);
}

dll Sparsesym* cb_sparsesym_support_xbpeya(Sparsesym* self, const Sparsesym* y, Real alpha = 1., Real beta = 0.) {
  return &self->support_xbpeya(*y, alpha, beta);
}

dll Sparsesym* cb_sparsesym_xbpeya(Sparsesym* x, const Sparsesym* y, Real alpha, Real beta) {
  return &xbpeya(*x, *y, alpha, beta);
}

dll Sparsesym* cb_sparsesym_xeyapzb(Sparsesym* x, const Sparsesym* y, const Sparsesym* z, Real alpha, Real beta) {
  return &xeyapzb(*x, *y, *z, alpha, beta);
}

dll Sparsesym* cb_sparsesym_support_rankadd(const Matrix* A, Sparsesym* C, Real alpha, Real beta, int trans) {
  return &support_rankadd(*A, *C, alpha, beta, trans);
}

dll Matrix* cb_sparsesym_genmult(const Sparsesym* A, const Matrix* B, Matrix* C, Real alpha, Real beta, int btrans) {
  return &genmult(*A, *B, *C, alpha, beta, btrans);
}

dll Matrix* cb_sparsesym_genmult2(const Matrix* A, const Sparsesym* B, Matrix* C, Real alpha, Real beta, int atrans) {
  return &genmult(*A, *B, *C, alpha, beta, atrans);
}

dll Sparsesym* cb_sparsesym_assign(Sparsesym* self, const Sparsesym* A) {
  return &(*self = *A);
}

dll Sparsesym* cb_sparsesym_plus(Sparsesym* self, const Sparsesym* v) {
  return &(*self += *v);
}

dll Sparsesym* cb_sparsesym_minus(Sparsesym* self, const Sparsesym* v) {
  return &(*self -= *v);
}

dll Sparsesym* cb_sparsesym_new_minus(const Sparsesym* self) {
  return new Sparsesym(-(*self));
}

dll Sparsesym* cb_sparsesym_times(Sparsesym* self, Real d) {
  return &(*self *= d);
}

dll Sparsesym* cb_sparsesym_divide(Sparsesym* self, Real d) {
  return &(*self /= d);
}

dll Sparsesym* cb_sparsesym_assign2(Sparsesym* self, const Matrix* A) {
  return &(*self = *A);
}

dll Sparsesym* cb_sparsesym_assign3(Sparsesym* self, const Indexmatrix* A) {
  return &(*self = *A);
}

dll Sparsesym* cb_sparsesym_transpose(Sparsesym* self) {
  return &self->transpose();
}

dll Sparsesym* cb_sparsesym_new_plus(const Sparsesym* A, const Sparsesym* B) {
  return new Sparsesym(*A + *B);
}

dll Sparsesym* cb_sparsesym_new_minus2(const Sparsesym* A, const Sparsesym* B) {
  return new Sparsesym(*A - *B);
}

dll Sparsesym* cb_sparsesym_new_times(const Sparsesym* A, Real d) {
  return new Sparsesym(*A * d);
}

dll Sparsesym* cb_sparsesym_new_times2(Real d, const Sparsesym* A) {
  return new Sparsesym(d * *A);
}

dll Sparsesym* cb_sparsesym_new_divide(const Sparsesym* A, Real d) {
  return new Sparsesym(*A / d);
}

dll Matrix* cb_sparsesym_new_times3(const Matrix* A, const Sparsesym* B) {
  return new Matrix(*A * *B);
}

dll Matrix* cb_sparsesym_new_times4(const Sparsesym* A, const Matrix* B) {
  return new Matrix(*A * *B);
}

dll Matrix* cb_sparsesym_new_plus2(const Matrix* A, const Sparsesym* B) {
  return new Matrix(*A + *B);
}

dll Matrix* cb_sparsesym_new_plus3(const Sparsesym* A, const Matrix* B) {
  return new Matrix(*A + *B);
}

dll Matrix* cb_sparsesym_new_minus3(const Matrix* A, const Sparsesym* B) {
  return new Matrix(*A - *B);
}

dll Matrix* cb_sparsesym_new_minus4(const Sparsesym* A, const Matrix* B) {
  return new Matrix(*A - *B);
}

dll Sparsesym* cb_sparsesym_new_transpose(const Sparsesym* A) {
  return new Sparsesym(transpose(*A));
}

dll Sparsesym* cb_sparsesym_xeya4(Sparsesym* self, const Symmatrix* A, Real d = 1.) {
  return &self->xeya(*A, d);
}

dll Sparsesym* cb_sparsesym_xeya5(Sparsesym* self, const Sparsemat* A, Real d = 1.) {
  return &self->xeya(*A, d);
}

dll Sparsesym* cb_sparsesym_assign4(Sparsesym* self, const Symmatrix* A) {
  return &(*self = *A);
}

dll Sparsesym* cb_sparsesym_assign5(Sparsesym* self, const Sparsemat* A) {
  return &(*self = *A);
}

dll Sparsemat* cb_sparsesym_new_sparsemult(const Sparsesym* self, const Matrix* A) {
  return new Sparsemat(self->sparsemult(*A));
}

dll Symmatrix* cb_sparsesym_new_plus4(const Sparsesym* A, const Symmatrix* B) {
  return new Symmatrix(*A + *B);
}

dll Symmatrix* cb_sparsesym_new_plus5(const Symmatrix* A, const Sparsesym* B) {
  return new Symmatrix(*A + *B);
}

dll Symmatrix* cb_sparsesym_new_minus5(const Sparsesym* A, const Symmatrix* B) {
  return new Symmatrix(*A - *B);
}

dll Symmatrix* cb_sparsesym_new_minus6(const Symmatrix* A, const Sparsesym* B) {
  return new Symmatrix(*A - *B);
}

dll Real cb_sparsesym_ip(const Symmatrix* A, const Sparsesym* B) {
  return ip(*A, *B);
}

dll Real cb_sparsesym_ip2(const Sparsesym* A, const Symmatrix* B) {
  return ip(*A, *B);
}

dll Matrix* cb_sparsesym_genmult3(const Sparsesym* A, const Sparsemat* B, Matrix* C, Real alpha, Real beta, int btrans) {
  return &genmult(*A, *B, *C, alpha, beta, btrans);
}

dll Matrix* cb_sparsesym_genmult4(const Sparsemat* A, const Sparsesym* B, Matrix* C, Real alpha, Real beta, int atrans) {
  return &genmult(*A, *B, *C, alpha, beta, atrans);
}

dll Matrix* cb_sparsesym_new_times5(const Sparsesym* A, const Sparsemat* B) {
  return new Matrix(*A * *B);
}

dll Matrix* cb_sparsesym_new_times6(const Sparsemat* A, const Sparsesym* B) {
  return new Matrix(*A * *B);
}

dll Sparsesym* cb_sparsesym_new_abs(const Sparsesym* A) {
  return new Sparsesym(abs(*A));
}

dll Real cb_sparsesym_trace(const Sparsesym* A) {
  return trace(*A);
}

dll Real cb_sparsesym_ip3(const Sparsesym* A, const Sparsesym* B) {
  return ip(*A, *B);
}

dll Real cb_sparsesym_ip4(const Matrix* A, const Sparsesym* B) {
  return ip(*A, *B);
}

dll Real cb_sparsesym_ip5(const Sparsesym* A, const Matrix* B) {
  return ip(*A, *B);
}

dll Real cb_sparsesym_norm2(const Sparsesym* A) {
  return norm2(*A);
}

dll Matrix* cb_sparsesym_new_sumrows(const Sparsesym* A) {
  return new Matrix(sumrows(*A));
}

dll Matrix* cb_sparsesym_new_sumcols(const Sparsesym* A) {
  return new Matrix(sumcols(*A));
}

dll Real cb_sparsesym_sum(const Sparsesym* A) {
  return sum(*A);
}

dll int cb_sparsesym_equal(const Sparsesym* A, const Sparsesym* B, Real eqtol) {
  return equal(*A, *B, eqtol);
}

dll void cb_sparsesym_display(const Sparsesym* self, int precision = 0, int width = 0, int screenwidth = 0) {
  self->display(std::cout, precision, width, screenwidth);
}

