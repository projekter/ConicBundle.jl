dll void cb_matrix_destroy(Matrix* self) {
  delete self;
}

dll Matrix* cb_matrix_new() {
  return new Matrix();
}

dll Matrix* cb_matrix_new2(const Matrix* param0, Real d = 1., int atrans = 0) {
  return new Matrix(*param0, d, atrans);
}

dll Matrix* cb_matrix_new3(Real param0_from, Real param0_to, Real param0_step = 1, Real param0_tol = 1e-8) {
  return new Matrix(Realrange(param0_from, param0_to, param0_step, param0_tol));
}

dll Matrix* cb_matrix_new4(Integer nr, Integer nc) {
  return new Matrix(nr, nc);
}

dll Matrix* cb_matrix_new5(Integer nr, Integer nc, Real d) {
  return new Matrix(nr, nc, d);
}

dll Matrix* cb_matrix_new6(Integer nr, Integer nc, const Real* dp, Integer incr = 1, Real d = 1.) {
  return new Matrix(nr, nc, dp, incr, d);
}

dll Matrix* cb_matrix_init(Matrix* self, const Matrix* A, Real d = 1., int atrans = 0) {
  return &self->init(*A, d, atrans);
}

dll Matrix* cb_matrix_init2(Matrix* self, const Indexmatrix* A, Real d = 1.) {
  return &self->init(*A, d);
}

dll Matrix* cb_matrix_init3(Matrix* self, const Sparsemat* A, Real d = 1.) {
  return &self->init(*A, d);
}

dll Matrix* cb_matrix_init4(Matrix* self, const Symmatrix* S, Real d = 1.) {
  return &self->init(*S, d);
}

dll Matrix* cb_matrix_init5(Matrix* self, const Sparsesym* param0, Real d = 1.) {
  return &self->init(*param0, d);
}

dll Matrix* cb_matrix_init6(Matrix* self, Real param0_from, Real param0_to, Real param0_step = 1, Real param0_tol = 1e-8) {
  return &self->init(Realrange(param0_from, param0_to, param0_step, param0_tol));
}

dll Matrix* cb_matrix_init7(Matrix* self, Integer nr, Integer nc, Real d) {
  return &self->init(nr, nc, d);
}

dll Matrix* cb_matrix_init8(Matrix* self, Integer nr, Integer nc, const Real* dp, Integer incr = 1, Real d = 1.) {
  return &self->init(nr, nc, dp, incr, d);
}

dll Matrix* cb_matrix_init_diag(Matrix* self, int nr, Real d = 1.) {
  return &self->init_diag(nr, d);
}

dll Matrix* cb_matrix_init_diag2(Matrix* self, const Matrix* vec, Real d = 1.) {
  return &self->init_diag(*vec, d);
}

dll Matrix* cb_matrix_init_diag3(Matrix* self, const Indexmatrix* vec, Real d = 1.) {
  return &self->init_diag(*vec, d);
}

dll void cb_matrix_newsize(Matrix* self, Integer nr, Integer nc) {
  self->newsize(nr, nc);
}

dll Matrix* cb_matrix_new7(const Indexmatrix* A, Real d = 1.) {
  return new Matrix(*A, d);
}

dll Matrix* cb_matrix_new8(const Sparsemat* A, Real d = 1.) {
  return new Matrix(*A, d);
}

dll Matrix* cb_matrix_new9(const Symmatrix* S, Real d = 1.) {
  return new Matrix(*S, d);
}

dll Matrix* cb_matrix_new10(const Sparsesym* param0, Real d = 1.) {
  return new Matrix(*param0, d);
}

dll void cb_matrix_dim(const Matrix* self, Integer* _nr, Integer* _nc) {
  self->dim(*_nr, *_nc);
}

dll Integer cb_matrix_dim2(const Matrix* self) {
  return self->dim();
}

dll Integer cb_matrix_rowdim(const Matrix* self) {
  return self->rowdim();
}

dll Integer cb_matrix_coldim(const Matrix* self) {
  return self->coldim();
}

dll void cb_matrix_set(Matrix* self, Integer i, Integer j, Real value) {
  (*self)(i, j) = value;
}

dll void cb_matrix_set2(Matrix* self, Integer i, Real value) {
  (*self)(i) = value;
}

dll Real cb_matrix_get(const Matrix* self, Integer i, Integer j) {
  return (*self)(i, j);
}

dll Real cb_matrix_get2(const Matrix* self, Integer i) {
  return (*self)(i);
}

dll Matrix* cb_matrix_new_get(const Matrix* self, const Indexmatrix* vecrow, const Indexmatrix* veccol) {
  return new Matrix((*self)(*vecrow, *veccol));
}

dll Matrix* cb_matrix_new_get2(const Matrix* self, const Indexmatrix* A) {
  return new Matrix((*self)(*A));
}

dll Matrix* cb_matrix_new_col(const Matrix* self, Integer i) {
  return new Matrix(self->col(i));
}

dll Matrix* cb_matrix_new_row(const Matrix* self, Integer i) {
  return new Matrix(self->row(i));
}

dll Matrix* cb_matrix_new_cols(const Matrix* self, const Indexmatrix* vec) {
  return new Matrix(self->cols(*vec));
}

dll Matrix* cb_matrix_new_rows(const Matrix* self, const Indexmatrix* vec) {
  return new Matrix(self->rows(*vec));
}

dll Matrix* cb_matrix_swap_rowsij(Matrix* self, Integer i, Integer j) {
  return &self->swap_rowsij(i, j);
}

dll Matrix* cb_matrix_swap_colsij(Matrix* self, Integer i, Integer j) {
  return &self->swap_colsij(i, j);
}

dll Matrix* cb_matrix_pivot_permute_rows(Matrix* self, const Indexmatrix* piv, int inverse = 0) {
  return &self->pivot_permute_rows(*piv, (bool)inverse);
}

dll Matrix* cb_matrix_pivot_permute_cols(Matrix* self, const Indexmatrix* piv, int inverse = 0) {
  return &self->pivot_permute_cols(*piv, (bool)inverse);
}

dll Matrix* cb_matrix_triu(Matrix* self, Integer d = 0) {
  return &self->triu(d);
}

dll Matrix* cb_matrix_tril(Matrix* self, Integer d = 0) {
  return &self->tril(d);
}

dll Matrix* cb_matrix_subassign(Matrix* self, const Indexmatrix* vecrow, const Indexmatrix* veccol, const Matrix* A) {
  return &self->subassign(*vecrow, *veccol, *A);
}

dll Matrix* cb_matrix_subassign2(Matrix* self, const Indexmatrix* vec, const Matrix* A) {
  return &self->subassign(*vec, *A);
}

dll Matrix* cb_matrix_delete_rows(Matrix* self, const Indexmatrix* ind, int sorted_increasingly = 0) {
  return &self->delete_rows(*ind, (bool)sorted_increasingly);
}

dll Matrix* cb_matrix_delete_cols(Matrix* self, const Indexmatrix* ind, int sorted_increasingly = 0) {
  return &self->delete_cols(*ind, (bool)sorted_increasingly);
}

dll Matrix* cb_matrix_insert_row(Matrix* self, Integer i, const Matrix* v) {
  return &self->insert_row(i, *v);
}

dll Matrix* cb_matrix_insert_col(Matrix* self, Integer i, const Matrix* v) {
  return &self->insert_col(i, *v);
}

dll Matrix* cb_matrix_reduce_length(Matrix* self, Integer n) {
  return &self->reduce_length(n);
}

dll Matrix* cb_matrix_concat_right(Matrix* self, const Matrix* A, int Atrans = 0) {
  return &self->concat_right(*A, Atrans);
}

dll Matrix* cb_matrix_concat_below(Matrix* self, const Matrix* A) {
  return &self->concat_below(*A);
}

dll Matrix* cb_matrix_concat_right2(Matrix* self, Real d) {
  return &self->concat_right(d);
}

dll Matrix* cb_matrix_concat_below2(Matrix* self, Real d) {
  return &self->concat_below(d);
}

dll Matrix* cb_matrix_enlarge_right(Matrix* self, Integer addnc) {
  return &self->enlarge_right(addnc);
}

dll Matrix* cb_matrix_enlarge_below(Matrix* self, Integer addnr) {
  return &self->enlarge_below(addnr);
}

dll Matrix* cb_matrix_enlarge_right2(Matrix* self, Integer addnc, Real d) {
  return &self->enlarge_right(addnc, d);
}

dll Matrix* cb_matrix_enlarge_below2(Matrix* self, Integer addnr, Real d) {
  return &self->enlarge_below(addnr, d);
}

dll Matrix* cb_matrix_enlarge_right3(Matrix* self, Integer addnc, const Real* dp, Real d = 1.) {
  return &self->enlarge_right(addnc, dp, d);
}

dll Matrix* cb_matrix_enlarge_below3(Matrix* self, Integer addnr, const Real* dp, Real d = 1.) {
  return &self->enlarge_below(addnr, dp, d);
}

dll const Real* cb_matrix_get_store2(const Matrix* self) {
  return self->get_store();
}

dll Matrix* cb_matrix_new_diag(const Matrix* A) {
  return new Matrix(diag(*A));
}

dll Matrix* cb_matrix_new_triu(const Matrix* A, Integer i) {
  return new Matrix(triu(*A, i));
}

dll Matrix* cb_matrix_new_tril(const Matrix* A, Integer i) {
  return new Matrix(tril(*A, i));
}

dll Matrix* cb_matrix_new_concat_right(const Matrix* A, const Matrix* B) {
  return new Matrix(concat_right(*A, *B));
}

dll Matrix* cb_matrix_new_concat_below(const Matrix* A, const Matrix* B) {
  return new Matrix(concat_below(*A, *B));
}

dll void cb_matrix_swap(Matrix* A, Matrix* B) {
  swap(*A, *B);
}

dll Matrix* cb_matrix_xeya(Matrix* self, const Matrix* A, Real d = 1., int atrans = 0) {
  return &self->xeya(*A, d, atrans);
}

dll Matrix* cb_matrix_xpeya(Matrix* self, const Matrix* A, Real d = 1.) {
  return &self->xpeya(*A, d);
}

dll Matrix* cb_matrix_xeya2(Matrix* self, const Indexmatrix* A, Real d = 1.) {
  return &self->xeya(*A, d);
}

dll Matrix* cb_matrix_xpeya2(Matrix* self, const Indexmatrix* A, Real d = 1.) {
  return &self->xpeya(*A, d);
}

dll Matrix* cb_matrix_xbpeya(Matrix* x, const Matrix* y, Real alpha, Real beta, int ytrans) {
  return &xbpeya(*x, *y, alpha, beta, ytrans);
}

dll Matrix* cb_matrix_xeyapzb(Matrix* x, const Matrix* y, const Matrix* z, Real alpha, Real beta) {
  return &xeyapzb(*x, *y, *z, alpha, beta);
}

dll Matrix* cb_matrix_genmult(const Matrix* A, const Matrix* B, Matrix* C, Real alpha, Real beta, int atrans, int btrans) {
  return &genmult(*A, *B, *C, alpha, beta, atrans, btrans);
}

dll Matrix* cb_matrix_assign(Matrix* self, const Matrix* A) {
  return &(*self = *A);
}

dll Matrix* cb_matrix_times(Matrix* self, const Matrix* s) {
  return &(*self *= *s);
}

dll Matrix* cb_matrix_plus(Matrix* self, const Matrix* v) {
  return &(*self += *v);
}

dll Matrix* cb_matrix_minus(Matrix* self, const Matrix* v) {
  return &(*self -= *v);
}

dll Matrix* cb_matrix_rem(Matrix* self, const Matrix* A) {
  return &(*self %= *A);
}

dll Matrix* cb_matrix_divide(Matrix* self, const Matrix* A) {
  return &(*self /= *A);
}

dll Matrix* cb_matrix_new_minus(const Matrix* self) {
  return new Matrix(-(*self));
}

dll Matrix* cb_matrix_times2(Matrix* self, Real d) {
  return &(*self *= d);
}

dll Matrix* cb_matrix_divide2(Matrix* self, Real d) {
  return &(*self /= d);
}

dll Matrix* cb_matrix_plus2(Matrix* self, Real d) {
  return &(*self += d);
}

dll Matrix* cb_matrix_minus2(Matrix* self, Real d) {
  return &(*self -= d);
}

dll Matrix* cb_matrix_transpose(Matrix* self) {
  return &self->transpose();
}

dll Matrix* cb_matrix_new_times(const Matrix* A, const Matrix* B) {
  return new Matrix(*A * *B);
}

dll Matrix* cb_matrix_new_plus(const Matrix* A, const Matrix* B) {
  return new Matrix(*A + *B);
}

dll Matrix* cb_matrix_new_minus2(const Matrix* A, const Matrix* B) {
  return new Matrix(*A - *B);
}

dll Matrix* cb_matrix_new_rem(const Matrix* A, const Matrix* B) {
  return new Matrix(*A % *B);
}

dll Matrix* cb_matrix_new_divide(const Matrix* A, const Matrix* B) {
  return new Matrix(*A / *B);
}

dll Matrix* cb_matrix_new_times2(const Matrix* A, Real d) {
  return new Matrix(*A * d);
}

dll Matrix* cb_matrix_new_times3(Real d, const Matrix* A) {
  return new Matrix(d * *A);
}

dll Matrix* cb_matrix_new_divide2(const Matrix* A, Real d) {
  return new Matrix(*A / d);
}

dll Matrix* cb_matrix_new_plus2(const Matrix* A, Real d) {
  return new Matrix(*A + d);
}

dll Matrix* cb_matrix_new_plus3(Real d, const Matrix* A) {
  return new Matrix(d + *A);
}

dll Matrix* cb_matrix_new_minus3(const Matrix* A, Real d) {
  return new Matrix(*A - d);
}

dll Matrix* cb_matrix_new_minus4(Real d, const Matrix* A) {
  return new Matrix(d - *A);
}

dll Matrix* cb_matrix_new_transpose(const Matrix* A) {
  return new Matrix(transpose(*A));
}

dll Matrix* cb_matrix_xeya3(Matrix* self, const Symmatrix* A, Real d = 1.) {
  return &self->xeya(*A, d);
}

dll Matrix* cb_matrix_xpeya3(Matrix* self, const Symmatrix* A, Real d = 1.) {
  return &self->xpeya(*A, d);
}

dll Matrix* cb_matrix_assign2(Matrix* self, const Symmatrix* S) {
  return &(*self = *S);
}

dll Matrix* cb_matrix_times3(Matrix* self, const Symmatrix* S) {
  return &(*self *= *S);
}

dll Matrix* cb_matrix_plus3(Matrix* self, const Symmatrix* S) {
  return &(*self += *S);
}

dll Matrix* cb_matrix_minus3(Matrix* self, const Symmatrix* S) {
  return &(*self -= *S);
}

dll Matrix* cb_matrix_xeya4(Matrix* self, const Sparsesym* A, Real d = 1.) {
  return &self->xeya(*A, d);
}

dll Matrix* cb_matrix_xpeya4(Matrix* self, const Sparsesym* A, Real d = 1.) {
  return &self->xpeya(*A, d);
}

dll Matrix* cb_matrix_assign3(Matrix* self, const Sparsesym* param0) {
  return &(*self = *param0);
}

dll Matrix* cb_matrix_times4(Matrix* self, const Sparsesym* S) {
  return &(*self *= *S);
}

dll Matrix* cb_matrix_plus4(Matrix* self, const Sparsesym* S) {
  return &(*self += *S);
}

dll Matrix* cb_matrix_minus4(Matrix* self, const Sparsesym* S) {
  return &(*self -= *S);
}

dll Matrix* cb_matrix_xeya5(Matrix* self, const Sparsemat* A, Real d = 1.) {
  return &self->xeya(*A, d);
}

dll Matrix* cb_matrix_xpeya5(Matrix* self, const Sparsemat* A, Real d = 1.) {
  return &self->xpeya(*A, d);
}

dll Matrix* cb_matrix_assign4(Matrix* self, const Sparsemat* A) {
  return &(*self = *A);
}

dll Matrix* cb_matrix_times5(Matrix* self, const Sparsemat* A) {
  return &(*self *= *A);
}

dll Matrix* cb_matrix_plus5(Matrix* self, const Sparsemat* A) {
  return &(*self += *A);
}

dll Matrix* cb_matrix_minus5(Matrix* self, const Sparsemat* A) {
  return &(*self -= *A);
}

dll Matrix* cb_matrix_genmult2(const Symmatrix* A, const Matrix* B, Matrix* C, Real alpha, Real beta, int btrans) {
  return &genmult(*A, *B, *C, alpha, beta, btrans);
}

dll Matrix* cb_matrix_genmult3(const Matrix* A, const Symmatrix* B, Matrix* C, Real alpha, Real beta, int atrans) {
  return &genmult(*A, *B, *C, alpha, beta, atrans);
}

dll Matrix* cb_matrix_genmult4(const Sparsesym* A, const Matrix* B, Matrix* C, Real alpha, Real beta, int btrans) {
  return &genmult(*A, *B, *C, alpha, beta, btrans);
}

dll Matrix* cb_matrix_genmult5(const Matrix* A, const Sparsesym* B, Matrix* C, Real alpha, Real beta, int atrans) {
  return &genmult(*A, *B, *C, alpha, beta, atrans);
}

dll Matrix* cb_matrix_genmult6(const Sparsemat* A, const Matrix* B, Matrix* C, Real alpha, Real beta, int atrans, int btrans) {
  return &genmult(*A, *B, *C, alpha, beta, atrans, btrans);
}

dll Matrix* cb_matrix_genmult7(const Sparsemat* A, const Matrix* B, Integer colB, Matrix* C, Integer colC, Real alpha, Real beta, int atrans, int btrans) {
  return &genmult(*A, *B, colB, *C, colC, alpha, beta, atrans, btrans);
}

dll Matrix* cb_matrix_genmult8(const Matrix* A, const Sparsemat* B, Matrix* C, Real alpha, Real beta, int atrans, int btrans) {
  return &genmult(*A, *B, *C, alpha, beta, atrans, btrans);
}

dll Matrix* cb_matrix_rand(Matrix* self, Integer nr, Integer nc, CH_Tools::GB_rand* random_generator = 0) {
  return &self->rand(nr, nc, random_generator);
}

dll Matrix* cb_matrix_rand_normal(Matrix* self, Integer nr, Integer nc, Real mean = 0., Real variance = 1., int generator_type = 0) {
  return &self->rand_normal(nr, nc, mean, variance, generator_type);
}

dll Matrix* cb_matrix_shuffle(Matrix* self, CH_Tools::GB_rand* random_generator = 0) {
  return &self->shuffle(random_generator);
}

dll Matrix* cb_matrix_inv(Matrix* self) {
  return &self->inv();
}

dll Matrix* cb_matrix_sqrt(Matrix* self) {
  return &self->sqrt();
}

dll Matrix* cb_matrix_sqr(Matrix* self) {
  return &self->sqr();
}

dll Matrix* cb_matrix_sign(Matrix* self, Real tol = 1e-12) {
  return &self->sign(tol);
}

dll Matrix* cb_matrix_floor(Matrix* self) {
  return &self->floor();
}

dll Matrix* cb_matrix_ceil(Matrix* self) {
  return &self->ceil();
}

dll Matrix* cb_matrix_rint(Matrix* self) {
  return &self->rint();
}

dll Matrix* cb_matrix_round(Matrix* self) {
  return &self->round();
}

dll Matrix* cb_matrix_abs(Matrix* self) {
  return &self->abs();
}

dll Matrix* cb_matrix_new_rand(Integer nr, Integer nc, CH_Tools::GB_rand* random_generator) {
  return new Matrix(rand(nr, nc, random_generator));
}

dll Matrix* cb_matrix_new_inv(const Matrix* A) {
  return new Matrix(inv(*A));
}

dll Matrix* cb_matrix_new_sqrt(const Matrix* A) {
  return new Matrix(sqrt(*A));
}

dll Matrix* cb_matrix_new_sqr(const Matrix* A) {
  return new Matrix(sqr(*A));
}

dll Matrix* cb_matrix_new_sign(const Matrix* A, Real tol) {
  return new Matrix(sign(*A, tol));
}

dll Matrix* cb_matrix_new_floor(const Matrix* A) {
  return new Matrix(floor(*A));
}

dll Matrix* cb_matrix_new_ceil(const Matrix* A) {
  return new Matrix(ceil(*A));
}

dll Matrix* cb_matrix_new_rint(const Matrix* A) {
  return new Matrix(rint(*A));
}

dll Matrix* cb_matrix_new_round(const Matrix* A) {
  return new Matrix(round(*A));
}

dll Matrix* cb_matrix_new_abs(const Matrix* A) {
  return new Matrix(abs(*A));
}

dll int cb_matrix_contains_nan(Matrix* self) {
  return self->contains_nan();
}

dll Matrix* cb_matrix_scale_rows(Matrix* self, const Matrix* vec) {
  return &self->scale_rows(*vec);
}

dll Matrix* cb_matrix_scale_cols(Matrix* self, const Matrix* vec) {
  return &self->scale_cols(*vec);
}

dll int cb_matrix_triu_solve(Matrix* self, Matrix* rhs, Real tol = 1e-10) {
  return self->triu_solve(*rhs, tol);
}

dll int cb_matrix_tril_solve(Matrix* self, Matrix* rhs, Real tol = 1e-10) {
  return self->tril_solve(*rhs, tol);
}

dll int cb_matrix_qr_factor(Matrix* self, Real tol = 1e-10) {
  return self->QR_factor(tol);
}

dll int cb_matrix_qr_factor2(Matrix* self, Matrix* Q, Real tol = 1e-10) {
  return self->QR_factor(*Q, tol);
}

dll int cb_matrix_qr_factor3(const Matrix* self, Matrix* Q, Matrix* R, Real tol) {
  return self->QR_factor(*Q, *R, tol);
}

dll int cb_matrix_qr_factor4(Matrix* self, Indexmatrix* piv, Real tol = 1e-10) {
  return self->QR_factor(*piv, tol);
}

dll int cb_matrix_qr_factor_relpiv(Matrix* self, Indexmatrix* piv, Real tol = 1e-10) {
  return self->QR_factor_relpiv(*piv, tol);
}

dll int cb_matrix_qr_factor5(Matrix* self, Matrix* Q, Indexmatrix* piv, Real tol = 1e-10) {
  return self->QR_factor(*Q, *piv, tol);
}

dll int cb_matrix_qr_factor6(const Matrix* self, Matrix* Q, Matrix* R, Indexmatrix* piv, Real tol) {
  return self->QR_factor(*Q, *R, *piv, tol);
}

dll int cb_matrix_qt_times(const Matrix* self, Matrix* A, Integer r) {
  return self->Qt_times(*A, r);
}

dll int cb_matrix_q_times(const Matrix* self, Matrix* A, Integer r) {
  return self->Q_times(*A, r);
}

dll int cb_matrix_times_q(const Matrix* self, Matrix* A, Integer r) {
  return self->times_Q(*A, r);
}

dll int cb_matrix_qr_solve(Matrix* self, Matrix* rhs, Real tol = 1e-10) {
  return self->QR_solve(*rhs, tol);
}

dll int cb_matrix_qr_concat_right(Matrix* self, const Matrix* A, Indexmatrix* piv, Integer r, Real tol = 1e-10) {
  return self->QR_concat_right(*A, *piv, r, tol);
}

dll int cb_matrix_ls(Matrix* self, Matrix* rhs, Real tol) {
  return self->ls(*rhs, tol);
}

dll int cb_matrix_nnls(const Matrix* self, Matrix* rhs, Matrix* dual = 0, Real tol = 1e-10) {
  return self->nnls(*rhs, dual, tol);
}

dll Real cb_matrix_trace(const Matrix* A) {
  return trace(*A);
}

dll Real cb_matrix_ip(const Matrix* A, const Matrix* B) {
  return ip(*A, *B);
}

dll Real cb_matrix_ip_min_max(const Matrix* A, const Matrix* B, Real* minval, Real* maxval) {
  return ip_min_max(*A, *B, *minval, *maxval);
}

dll Real cb_matrix_colip(const Matrix* A, Integer j, const Matrix* scaling) {
  return colip(*A, j, scaling);
}

dll Real cb_matrix_rowip(const Matrix* A, Integer i, const Matrix* scaling) {
  return rowip(*A, i, scaling);
}

dll Matrix* cb_matrix_new_colsip(const Matrix* A) {
  return new Matrix(colsip(*A));
}

dll Matrix* cb_matrix_new_rowsip(const Matrix* A) {
  return new Matrix(rowsip(*A));
}

dll Real cb_matrix_norm2(const Matrix* A) {
  return norm2(*A);
}

dll Real cb_matrix_normdsquared(const Matrix* A, const Matrix* d, int atrans, int dinv) {
  return normDsquared(*A, *d, atrans, dinv);
}

dll Matrix* cb_matrix_new_sumrows(const Matrix* A) {
  return new Matrix(sumrows(*A));
}

dll Matrix* cb_matrix_new_sumcols(const Matrix* A) {
  return new Matrix(sumcols(*A));
}

dll Real cb_matrix_sum(const Matrix* A) {
  return sum(*A);
}

dll Matrix* cb_matrix_new_house(const Matrix* A, Integer i, Integer j, Real tol) {
  return new Matrix(house(*A, i, j, tol));
}

dll int cb_matrix_rowhouse(Matrix* A, const Matrix* v, Integer i, Integer j) {
  return rowhouse(*A, *v, i, j);
}

dll int cb_matrix_colhouse(Matrix* A, const Matrix* v, Integer i, Integer j) {
  return colhouse(*A, *v, i, j);
}

dll int cb_matrix_qr_factor7(const Matrix* A, Matrix* Q, Matrix* R, Real tol) {
  return QR_factor(*A, *Q, *R, tol);
}

dll int cb_matrix_qr_factor8(const Matrix* A, Matrix* Q, Matrix* R, Indexmatrix* piv, Real tol) {
  return QR_factor(*A, *Q, *R, *piv, tol);
}

dll Indexmatrix* cb_matrix_new_find(const Matrix* self, Real tol = 1e-10) {
  return new Indexmatrix(self->find(tol));
}

dll Indexmatrix* cb_matrix_new_find_number(const Matrix* self, Real num = 0., Real tol = 1e-10) {
  return new Indexmatrix(self->find_number(num, tol));
}

dll Matrix* cb_matrix_new_less(const Matrix* A, const Matrix* B) {
  return new Matrix(*A < *B);
}

dll Matrix* cb_matrix_new_greater(const Matrix* A, const Matrix* B) {
  return new Matrix(*A > *B);
}

dll Matrix* cb_matrix_new_lessequal(const Matrix* A, const Matrix* B) {
  return new Matrix(*A <= *B);
}

dll Matrix* cb_matrix_new_greaterequal(const Matrix* A, const Matrix* B) {
  return new Matrix(*A >= *B);
}

dll Matrix* cb_matrix_new_equal(const Matrix* A, const Matrix* B) {
  return new Matrix(*A == *B);
}

dll Matrix* cb_matrix_new_inequal(const Matrix* A, const Matrix* B) {
  return new Matrix(*A != *B);
}

dll Matrix* cb_matrix_new_less2(const Matrix* A, Real d) {
  return new Matrix(*A < d);
}

dll Matrix* cb_matrix_new_greater2(const Matrix* A, Real d) {
  return new Matrix(*A > d);
}

dll Matrix* cb_matrix_new_lessequal2(const Matrix* A, Real d) {
  return new Matrix(*A <= d);
}

dll Matrix* cb_matrix_new_greaterequal2(const Matrix* A, Real d) {
  return new Matrix(*A >= d);
}

dll Matrix* cb_matrix_new_equal2(const Matrix* A, Real d) {
  return new Matrix(*A == d);
}

dll Matrix* cb_matrix_new_inequal2(const Matrix* A, Real d) {
  return new Matrix(*A != d);
}

dll Matrix* cb_matrix_new_less3(Real d, const Matrix* A) {
  return new Matrix(d < *A);
}

dll Matrix* cb_matrix_new_greater3(Real d, const Matrix* A) {
  return new Matrix(d > *A);
}

dll Matrix* cb_matrix_new_lessequal3(Real d, const Matrix* A) {
  return new Matrix(d <= *A);
}

dll Matrix* cb_matrix_new_greaterequal3(Real d, const Matrix* A) {
  return new Matrix(d >= *A);
}

dll Matrix* cb_matrix_new_equal3(Real d, const Matrix* A) {
  return new Matrix(d == *A);
}

dll Matrix* cb_matrix_new_inequal3(Real d, const Matrix* A) {
  return new Matrix(d != *A);
}

dll int cb_matrix_equal(const Matrix* A, const Matrix* B) {
  return equal(*A, *B);
}

dll Matrix* cb_matrix_new_minrows(const Matrix* A) {
  return new Matrix(minrows(*A));
}

dll Matrix* cb_matrix_new_mincols(const Matrix* A) {
  return new Matrix(mincols(*A));
}

dll Real cb_matrix_min(const Matrix* A, Integer* iindex, Integer* jindex) {
  return min(*A, iindex, jindex);
}

dll Matrix* cb_matrix_new_maxrows(const Matrix* A) {
  return new Matrix(maxrows(*A));
}

dll Matrix* cb_matrix_new_maxcols(const Matrix* A) {
  return new Matrix(maxcols(*A));
}

dll Real cb_matrix_max(const Matrix* A, Integer* iindex, Integer* jindex) {
  return max(*A, iindex, jindex);
}

dll Indexmatrix* cb_matrix_new_sortindex(const Matrix* vec, int nondecreasing) {
  return new Indexmatrix(sortindex(*vec, (bool)nondecreasing));
}

dll void cb_matrix_sortindex(const Matrix* vec, Indexmatrix* ind, int nondecreasing) {
  sortindex(*vec, *ind, (bool)nondecreasing);
}

dll void cb_matrix_display(const Matrix* self, int precision = 0, int width = 0, int screenwidth = 0) {
  self->display(std::cout, precision, width, screenwidth);
}

dll void cb_matrix_mfile_output(const Matrix* self, int precision = 16, int width = 0) {
  self->mfile_output(std::cout, precision, width);
}

