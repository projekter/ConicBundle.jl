dll void cb_sparsemat_destroy(Sparsemat* self) {
  delete self;
}

dll Sparsemat* cb_sparsemat_new() {
  return new Sparsemat();
}

dll Sparsemat* cb_sparsemat_new2(const Sparsemat* A, Real d = 1.) {
  return new Sparsemat(*A, d);
}

dll Sparsemat* cb_sparsemat_new3(Integer nr, Integer nc) {
  return new Sparsemat(nr, nc);
}

dll Sparsemat* cb_sparsemat_new4(Integer nr, Integer nc, Integer nz, const Integer* ini, const Integer* inj, const Real* va) {
  return new Sparsemat(nr, nc, nz, ini, inj, va);
}

dll Sparsemat* cb_sparsemat_new5(Integer nr, Integer nc, Integer nz, const Indexmatrix* ini, const Indexmatrix* inj, const Matrix* va) {
  return new Sparsemat(nr, nc, nz, *ini, *inj, *va);
}

dll Sparsemat* cb_sparsemat_init(Sparsemat* self, const Sparsemat* A, Real d = 1.) {
  return &self->init(*A, d);
}

dll Sparsemat* cb_sparsemat_init2(Sparsemat* self, const Matrix* A, Real d = 1.) {
  return &self->init(*A, d);
}

dll Sparsemat* cb_sparsemat_init3(Sparsemat* self, const Indexmatrix* A, Real d = 1.) {
  return &self->init(*A, d);
}

dll Sparsemat* cb_sparsemat_init4(Sparsemat* self, const Symmatrix* param0, Real d = 1.) {
  return &self->init(*param0, d);
}

dll Sparsemat* cb_sparsemat_init5(Sparsemat* self, const Sparsesym* param0, Real d = 1.) {
  return &self->init(*param0, d);
}

dll Sparsemat* cb_sparsemat_init6(Sparsemat* self, Integer nr, Integer nc) {
  return &self->init(nr, nc);
}

dll Sparsemat* cb_sparsemat_init7(Sparsemat* self, Integer nr, Integer nc, Integer nz, const Integer* ini, const Integer* inj, const Real* va) {
  return &self->init(nr, nc, nz, ini, inj, va);
}

dll Sparsemat* cb_sparsemat_init8(Sparsemat* self, Integer nr, Integer nc, Integer nz, const Indexmatrix* ini, const Indexmatrix* inj, const Matrix* va) {
  return &self->init(nr, nc, nz, *ini, *inj, *va);
}

dll void cb_sparsemat_set_tol(Sparsemat* self, Real t) {
  self->set_tol(t);
}

dll Sparsemat* cb_sparsemat_new6(const Matrix* A, Real d = 1.) {
  return new Sparsemat(*A, d);
}

dll Sparsemat* cb_sparsemat_new7(const Indexmatrix* A, Real d = 1.) {
  return new Sparsemat(*A, d);
}

dll Sparsemat* cb_sparsemat_new8(const Symmatrix* A, Real d = 1.) {
  return new Sparsemat(*A, d);
}

dll Sparsemat* cb_sparsemat_new9(const Sparsesym* A, Real d = 1.) {
  return new Sparsemat(*A, d);
}

dll void cb_sparsemat_dim(const Sparsemat* self, Integer* in_nr, Integer* in_nc) {
  self->dim(*in_nr, *in_nc);
}

dll Integer cb_sparsemat_dim2(const Sparsemat* self) {
  return self->dim();
}

dll Integer cb_sparsemat_rowdim(const Sparsemat* self) {
  return self->rowdim();
}

dll Integer cb_sparsemat_coldim(const Sparsemat* self) {
  return self->coldim();
}

dll Integer cb_sparsemat_nonzeros(const Sparsemat* self) {
  return self->nonzeros();
}

dll Integer cb_sparsemat_col_nonzeros(const Sparsemat* self, Integer i, Integer* startind = 0) {
  return self->col_nonzeros(i, startind);
}

dll Integer cb_sparsemat_row_nonzeros(const Sparsemat* self, Integer i, Integer* startind = 0) {
  return self->row_nonzeros(i, startind);
}

dll Real cb_sparsemat_get(const Sparsemat* self, Integer i, Integer j) {
  return (*self)(i, j);
}

dll Real cb_sparsemat_get2(const Sparsemat* self, Integer i) {
  return (*self)(i);
}

dll Sparsemat* cb_sparsemat_new_col(const Sparsemat* self, Integer i) {
  return new Sparsemat(self->col(i));
}

dll Sparsemat* cb_sparsemat_new_row(const Sparsemat* self, Integer i) {
  return new Sparsemat(self->row(i));
}

dll Sparsemat* cb_sparsemat_new_cols(const Sparsemat* self, const Indexmatrix* ind) {
  return new Sparsemat(self->cols(*ind));
}

dll Sparsemat* cb_sparsemat_new_rows(const Sparsemat* self, const Indexmatrix* ind) {
  return new Sparsemat(self->rows(*ind));
}

dll Sparsemat* cb_sparsemat_delete_rows(Sparsemat* self, const Indexmatrix* ind) {
  return &self->delete_rows(*ind);
}

dll Sparsemat* cb_sparsemat_delete_cols(Sparsemat* self, const Indexmatrix* ind) {
  return &self->delete_cols(*ind);
}

dll Sparsemat* cb_sparsemat_insert_row(Sparsemat* self, Integer i, const Sparsemat* v) {
  return &self->insert_row(i, *v);
}

dll Sparsemat* cb_sparsemat_insert_col(Sparsemat* self, Integer i, const Sparsemat* v) {
  return &self->insert_col(i, *v);
}

dll Sparsemat* cb_sparsemat_concat_right(Sparsemat* self, const Sparsemat* A) {
  return &self->concat_right(*A);
}

dll Sparsemat* cb_sparsemat_concat_below(Sparsemat* self, const Sparsemat* A) {
  return &self->concat_below(*A);
}

dll const Indexmatrix* cb_sparsemat_get_colinfo(const Sparsemat* self) {
  return &self->get_colinfo();
}

dll const Indexmatrix* cb_sparsemat_get_colindex(const Sparsemat* self) {
  return &self->get_colindex();
}

dll const Matrix* cb_sparsemat_get_colval(const Sparsemat* self) {
  return &self->get_colval();
}

dll const Indexmatrix* cb_sparsemat_get_rowinfo(const Sparsemat* self) {
  return &self->get_rowinfo();
}

dll const Indexmatrix* cb_sparsemat_get_rowindex(const Sparsemat* self) {
  return &self->get_rowindex();
}

dll const Matrix* cb_sparsemat_get_rowval(const Sparsemat* self) {
  return &self->get_rowval();
}

dll void cb_sparsemat_get_edge_rep(const Sparsemat* self, Indexmatrix* I, Indexmatrix* J, Matrix* val) {
  self->get_edge_rep(*I, *J, *val);
}

dll int cb_sparsemat_get_edge(const Sparsemat* self, Integer i, Integer* indi, Integer* indj, Real* val) {
  return self->get_edge(i, *indi, *indj, *val);
}

dll int cb_sparsemat_contains_support(const Sparsemat* self, const Sparsemat* A) {
  return self->contains_support(*A);
}

dll Sparsemat* cb_sparsemat_new_concat_right(const Sparsemat* A, const Sparsemat* B) {
  return new Sparsemat(concat_right(*A, *B));
}

dll Sparsemat* cb_sparsemat_new_concat_below(const Sparsemat* A, const Sparsemat* B) {
  return new Sparsemat(concat_below(*A, *B));
}

dll void cb_sparsemat_swap(Sparsemat* A, Sparsemat* B) {
  swap(*A, *B);
}

dll Sparsemat* cb_sparsemat_xeya(Sparsemat* self, const Sparsemat* A, Real d = 1.) {
  return &self->xeya(*A, d);
}

dll Sparsemat* cb_sparsemat_xeya2(Sparsemat* self, const Matrix* A, Real d = 1.) {
  return &self->xeya(*A, d);
}

dll Sparsemat* cb_sparsemat_xeya3(Sparsemat* self, const Indexmatrix* A, Real d = 1.) {
  return &self->xeya(*A, d);
}

dll Sparsemat* cb_sparsemat_xbpeya(Sparsemat* x, const Sparsemat* y, Real alpha, Real beta, int ytrans) {
  return &xbpeya(*x, *y, alpha, beta, ytrans);
}

dll Matrix* cb_sparsemat_genmult(const Sparsemat* A, const Matrix* B, Matrix* C, Real alpha, Real beta, int atrans, int btrans) {
  return &genmult(*A, *B, *C, alpha, beta, atrans, btrans);
}

dll Matrix* cb_sparsemat_genmult2(const Sparsemat* A, const Matrix* B, Integer colB, Matrix* C, Integer colC, Real alpha, Real beta, int atrans, int btrans) {
  return &genmult(*A, *B, colB, *C, colC, alpha, beta, atrans, btrans);
}

dll Matrix* cb_sparsemat_genmult3(const Matrix* A, const Sparsemat* B, Matrix* C, Real alpha, Real beta, int atrans, int btrans) {
  return &genmult(*A, *B, *C, alpha, beta, atrans, btrans);
}

dll Matrix* cb_sparsemat_genmult4(const Sparsemat* A, const Sparsemat* B, Matrix* C, Real alpha, Real beta, int atrans, int btrans) {
  return &genmult(*A, *B, *C, alpha, beta, atrans, btrans);
}

dll Sparsemat* cb_sparsemat_assign(Sparsemat* self, const Sparsemat* A) {
  return &(*self = *A);
}

dll Sparsemat* cb_sparsemat_plus(Sparsemat* self, const Sparsemat* A) {
  return &(*self += *A);
}

dll Sparsemat* cb_sparsemat_minus(Sparsemat* self, const Sparsemat* A) {
  return &(*self -= *A);
}

dll Sparsemat* cb_sparsemat_new_minus(const Sparsemat* self) {
  return new Sparsemat(-(*self));
}

dll Sparsemat* cb_sparsemat_times(Sparsemat* self, Real d) {
  return &(*self *= d);
}

dll Sparsemat* cb_sparsemat_divide(Sparsemat* self, Real d) {
  return &(*self /= d);
}

dll Sparsemat* cb_sparsemat_assign2(Sparsemat* self, const Matrix* A) {
  return &(*self = *A);
}

dll Sparsemat* cb_sparsemat_assign3(Sparsemat* self, const Indexmatrix* A) {
  return &(*self = *A);
}

dll Sparsemat* cb_sparsemat_transpose(Sparsemat* self) {
  return &self->transpose();
}

dll Sparsemat* cb_sparsemat_new_times(const Sparsemat* A, const Sparsemat* B) {
  return new Sparsemat(*A * *B);
}

dll Sparsemat* cb_sparsemat_new_plus(const Sparsemat* A, const Sparsemat* B) {
  return new Sparsemat(*A + *B);
}

dll Sparsemat* cb_sparsemat_new_minus2(const Sparsemat* A, const Sparsemat* B) {
  return new Sparsemat(*A - *B);
}

dll Sparsemat* cb_sparsemat_new_times2(const Sparsemat* A, Real d) {
  return new Sparsemat(*A * d);
}

dll Sparsemat* cb_sparsemat_new_times3(Real d, const Sparsemat* A) {
  return new Sparsemat(d * *A);
}

dll Sparsemat* cb_sparsemat_new_divide(const Sparsemat* A, Real d) {
  return new Sparsemat(*A / d);
}

dll Matrix* cb_sparsemat_new_times4(const Sparsemat* A, const Matrix* B) {
  return new Matrix(*A * *B);
}

dll Matrix* cb_sparsemat_new_times5(const Matrix* A, const Sparsemat* B) {
  return new Matrix(*A * *B);
}

dll Matrix* cb_sparsemat_new_plus2(const Sparsemat* A, const Matrix* B) {
  return new Matrix(*A + *B);
}

dll Matrix* cb_sparsemat_new_plus3(const Matrix* A, const Sparsemat* B) {
  return new Matrix(*A + *B);
}

dll Matrix* cb_sparsemat_new_minus3(const Sparsemat* A, const Matrix* B) {
  return new Matrix(*A - *B);
}

dll Matrix* cb_sparsemat_new_minus4(const Matrix* A, const Sparsemat* B) {
  return new Matrix(*A - *B);
}

dll Sparsemat* cb_sparsemat_new_transpose(const Sparsemat* A) {
  return new Sparsemat(transpose(*A));
}

dll Sparsemat* cb_sparsemat_xeya4(Sparsemat* self, const Sparsesym* A, Real d = 1.) {
  return &self->xeya(*A, d);
}

dll Sparsemat* cb_sparsemat_assign4(Sparsemat* self, const Symmatrix* A) {
  return &(*self = *A);
}

dll Sparsemat* cb_sparsemat_assign5(Sparsemat* self, const Sparsesym* A) {
  return &(*self = *A);
}

dll Matrix* cb_sparsemat_genmult5(const Symmatrix* A, const Sparsemat* B, Matrix* C, Real alpha, Real beta, int btrans) {
  return &genmult(*A, *B, *C, alpha, beta, btrans);
}

dll Matrix* cb_sparsemat_genmult6(const Sparsemat* A, const Symmatrix* B, Matrix* C, Real alpha, Real beta, int atrans) {
  return &genmult(*A, *B, *C, alpha, beta, atrans);
}

dll Matrix* cb_sparsemat_genmult7(const Sparsesym* A, const Sparsemat* B, Matrix* C, Real alpha, Real beta, int btrans) {
  return &genmult(*A, *B, *C, alpha, beta, btrans);
}

dll Matrix* cb_sparsemat_genmult8(const Sparsemat* A, const Sparsesym* B, Matrix* C, Real alpha, Real beta, int atrans) {
  return &genmult(*A, *B, *C, alpha, beta, atrans);
}

dll Symmatrix* cb_sparsemat_rankadd(const Sparsemat* A, Symmatrix* C, Real alpha, Real beta, int trans) {
  return &rankadd(*A, *C, alpha, beta, trans);
}

dll Symmatrix* cb_sparsemat_scaledrankadd(const Sparsemat* A, const Matrix* D, Symmatrix* C, Real alpha, Real beta, int trans) {
  return &scaledrankadd(*A, *D, *C, alpha, beta, trans);
}

dll Symmatrix* cb_sparsemat_rank2add(const Sparsemat* A, const Matrix* B, Symmatrix* C, Real alpha, Real beta, int trans) {
  return &rank2add(*A, *B, *C, alpha, beta, trans);
}

dll Matrix* cb_sparsemat_new_times6(const Symmatrix* A, const Sparsemat* B) {
  return new Matrix(*A * *B);
}

dll Matrix* cb_sparsemat_new_times7(const Sparsemat* A, const Symmatrix* B) {
  return new Matrix(*A * *B);
}

dll Sparsemat* cb_sparsemat_new_abs(const Sparsemat* A) {
  return new Sparsemat(abs(*A));
}

dll Sparsemat* cb_sparsemat_scale_rows(Sparsemat* self, const Matrix* vec) {
  return &self->scale_rows(*vec);
}

dll Sparsemat* cb_sparsemat_scale_cols(Sparsemat* self, const Matrix* vec) {
  return &self->scale_cols(*vec);
}

dll Real cb_sparsemat_trace(const Sparsemat* A) {
  return trace(*A);
}

dll Real cb_sparsemat_ip(const Sparsemat* A, const Sparsemat* B) {
  return ip(*A, *B);
}

dll Real cb_sparsemat_ip2(const Sparsemat* A, const Matrix* B) {
  return ip(*A, *B);
}

dll Real cb_sparsemat_ip3(const Matrix* A, const Sparsemat* B) {
  return ip(*A, *B);
}

dll Real cb_sparsemat_colip(const Sparsemat* A, Integer j, const Matrix* scaling) {
  return colip(*A, j, scaling);
}

dll Real cb_sparsemat_rowip(const Sparsemat* A, Integer i, const Matrix* scaling) {
  return rowip(*A, i, scaling);
}

dll Sparsemat* cb_sparsemat_new_colsip(const Sparsemat* A, const Matrix* scaling) {
  return new Sparsemat(colsip(*A, scaling));
}

dll Sparsemat* cb_sparsemat_new_rowsip(const Sparsemat* A, const Matrix* scaling) {
  return new Sparsemat(rowsip(*A, scaling));
}

dll Real cb_sparsemat_norm2(const Sparsemat* A) {
  return norm2(*A);
}

dll Matrix* cb_sparsemat_new_sumrows(const Sparsemat* A) {
  return new Matrix(sumrows(*A));
}

dll Matrix* cb_sparsemat_new_sumcols(const Sparsemat* A) {
  return new Matrix(sumcols(*A));
}

dll Real cb_sparsemat_sum(const Sparsemat* A) {
  return sum(*A);
}

dll int cb_sparsemat_equal(const Sparsemat* A, const Sparsemat* B, Real eqtol) {
  return equal(*A, *B, eqtol);
}

dll void cb_sparsemat_display(const Sparsemat* self, int precision = 0, int width = 0, int screenwidth = 0) {
  self->display(std::cout, precision, width, screenwidth);
}

