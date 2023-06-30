dll void cb_symmatrix_destroy(Symmatrix* self) {
  delete self;
}

dll Symmatrix* cb_symmatrix_new() {
  return new Symmatrix();
}

dll Symmatrix* cb_symmatrix_new2(const Symmatrix* A, double d = 1.) {
  return new Symmatrix(*A, d);
}

dll Symmatrix* cb_symmatrix_new3(Integer nr) {
  return new Symmatrix(nr);
}

dll Symmatrix* cb_symmatrix_new4(Integer nr, Real d) {
  return new Symmatrix(nr, d);
}

dll Symmatrix* cb_symmatrix_new5(Integer nr, const Real* dp) {
  return new Symmatrix(nr, dp);
}

dll Symmatrix* cb_symmatrix_init(Symmatrix* self, const Symmatrix* A, double d = 1.) {
  return &self->init(*A, d);
}

dll Symmatrix* cb_symmatrix_init2(Symmatrix* self, const Matrix* A, double d = 1.) {
  return &self->init(*A, d);
}

dll Symmatrix* cb_symmatrix_init3(Symmatrix* self, const Indexmatrix* A, double d = 1.) {
  return &self->init(*A, d);
}

dll Symmatrix* cb_symmatrix_init4(Symmatrix* self, const Sparsesym* A, Real d = 1.) {
  return &self->init(*A, d);
}

dll Symmatrix* cb_symmatrix_init5(Symmatrix* self, Integer nr, Real d) {
  return &self->init(nr, d);
}

dll Symmatrix* cb_symmatrix_init6(Symmatrix* self, Integer nr, const Real* dp) {
  return &self->init(nr, dp);
}

dll void cb_symmatrix_newsize(Symmatrix* self, Integer n) {
  self->newsize(n);
}

dll Symmatrix* cb_symmatrix_new6(const Matrix* param0, double d = 1.) {
  return new Symmatrix(*param0, d);
}

dll Symmatrix* cb_symmatrix_new7(const Indexmatrix* param0, double d = 1.) {
  return new Symmatrix(*param0, d);
}

dll Symmatrix* cb_symmatrix_new8(const Sparsesym* A, Real d = 1.) {
  return new Symmatrix(*A, d);
}

dll void cb_symmatrix_dim(const Symmatrix* self, Integer* _nr, Integer* _nc) {
  self->dim(*_nr, *_nc);
}

dll Integer cb_symmatrix_dim2(const Symmatrix* self) {
  return self->dim();
}

dll Integer cb_symmatrix_rowdim(const Symmatrix* self) {
  return self->rowdim();
}

dll Integer cb_symmatrix_coldim(const Symmatrix* self) {
  return self->coldim();
}

dll void cb_symmatrix_set(Symmatrix* self, Integer i, Integer j, Real value) {
  (*self)(i, j) = value;
}

dll void cb_symmatrix_set2(Symmatrix* self, Integer i, Real value) {
  (*self)(i) = value;
}

dll Real cb_symmatrix_get(const Symmatrix* self, Integer i, Integer j) {
  return (*self)(i, j);
}

dll Real cb_symmatrix_get2(const Symmatrix* self, Integer i) {
  return (*self)(i);
}

dll Matrix* cb_symmatrix_new_col(const Symmatrix* self, Integer i) {
  return new Matrix(self->col(i));
}

dll Matrix* cb_symmatrix_new_row(const Symmatrix* self, Integer i) {
  return new Matrix(self->row(i));
}

dll Matrix* cb_symmatrix_new_cols(const Symmatrix* self, const Indexmatrix* vec) {
  return new Matrix(self->cols(*vec));
}

dll Matrix* cb_symmatrix_new_rows(const Symmatrix* self, const Indexmatrix* vec) {
  return new Matrix(self->rows(*vec));
}

dll Symmatrix* cb_symmatrix_swapij(Symmatrix* self, Integer i, Integer j) {
  return &self->swapij(i, j);
}

dll Symmatrix* cb_symmatrix_pivot_permute(Symmatrix* self, const Indexmatrix* piv, int inverse = 0) {
  return &self->pivot_permute(*piv, (bool)inverse);
}

dll Symmatrix* cb_symmatrix_principal_submatrix(const Symmatrix* self, const Indexmatrix* ind, Symmatrix* S) {
  return &self->principal_submatrix(*ind, *S);
}

dll Symmatrix* cb_symmatrix_new_principal_submatrix(const Symmatrix* self, const Indexmatrix* ind) {
  return new Symmatrix(self->principal_submatrix(*ind));
}

dll Symmatrix* cb_symmatrix_delete_principal_submatrix(Symmatrix* self, const Indexmatrix* ind, int sorted_increasingly = 0) {
  return &self->delete_principal_submatrix(*ind, (bool)sorted_increasingly);
}

dll Symmatrix* cb_symmatrix_enlarge_below(Symmatrix* self, Integer addn) {
  return &self->enlarge_below(addn);
}

dll Symmatrix* cb_symmatrix_enlarge_below2(Symmatrix* self, Integer addn, Real d) {
  return &self->enlarge_below(addn, d);
}

dll const Real* cb_symmatrix_get_store2(const Symmatrix* self) {
  return self->get_store();
}

dll Matrix* cb_symmatrix_new_diag(const Symmatrix* A) {
  return new Matrix(diag(*A));
}

dll Symmatrix* cb_symmatrix_new_diag2(const Matrix* A) {
  return new Symmatrix(Diag(*A));
}

dll void cb_symmatrix_swap(Symmatrix* A, Symmatrix* B) {
  swap(*A, *B);
}

dll Symmatrix* cb_symmatrix_xeya(Symmatrix* self, const Symmatrix* A, Real d = 1.) {
  return &self->xeya(*A, d);
}

dll Symmatrix* cb_symmatrix_xpeya(Symmatrix* self, const Symmatrix* A, Real d = 1.) {
  return &self->xpeya(*A, d);
}

dll Symmatrix* cb_symmatrix_rankadd(const Matrix* A, Symmatrix* C, Real alpha, Real beta, int trans) {
  return &rankadd(*A, *C, alpha, beta, trans);
}

dll Symmatrix* cb_symmatrix_scaledrankadd(const Matrix* A, const Matrix* D, Symmatrix* C, Real alpha, Real beta, int trans) {
  return &scaledrankadd(*A, *D, *C, alpha, beta, trans);
}

dll Symmatrix* cb_symmatrix_rank2add(const Matrix* A, const Matrix* B, Symmatrix* C, Real alpha, Real beta, int trans) {
  return &rank2add(*A, *B, *C, alpha, beta, trans);
}

dll Symmatrix* cb_symmatrix_xbpeya(Symmatrix* x, const Symmatrix* y, Real alpha, Real beta) {
  return &xbpeya(*x, *y, alpha, beta);
}

dll Symmatrix* cb_symmatrix_xeyapzb(Symmatrix* x, const Symmatrix* y, const Symmatrix* z, Real alpha, Real beta) {
  return &xeyapzb(*x, *y, *z, alpha, beta);
}

dll Matrix* cb_symmatrix_genmult(const Symmatrix* A, const Matrix* B, Matrix* C, Real alpha, Real beta, int btrans) {
  return &genmult(*A, *B, *C, alpha, beta, btrans);
}

dll Matrix* cb_symmatrix_genmult2(const Matrix* A, const Symmatrix* B, Matrix* C, Real alpha, Real beta, int atrans) {
  return &genmult(*A, *B, *C, alpha, beta, atrans);
}

dll Symmatrix* cb_symmatrix_assign(Symmatrix* self, const Symmatrix* A) {
  return &(*self = *A);
}

dll Symmatrix* cb_symmatrix_plus(Symmatrix* self, const Symmatrix* A) {
  return &(*self += *A);
}

dll Symmatrix* cb_symmatrix_minus(Symmatrix* self, const Symmatrix* A) {
  return &(*self -= *A);
}

dll Symmatrix* cb_symmatrix_rem(Symmatrix* self, const Symmatrix* A) {
  return &(*self %= *A);
}

dll Symmatrix* cb_symmatrix_new_minus(const Symmatrix* self) {
  return new Symmatrix(-(*self));
}

dll Symmatrix* cb_symmatrix_times(Symmatrix* self, Real d) {
  return &(*self *= d);
}

dll Symmatrix* cb_symmatrix_divide(Symmatrix* self, Real d) {
  return &(*self /= d);
}

dll Symmatrix* cb_symmatrix_plus2(Symmatrix* self, Real d) {
  return &(*self += d);
}

dll Symmatrix* cb_symmatrix_minus2(Symmatrix* self, Real d) {
  return &(*self -= d);
}

dll Symmatrix* cb_symmatrix_transpose(Symmatrix* self) {
  return &self->transpose();
}

dll Matrix* cb_symmatrix_new_times(const Symmatrix* A, const Symmatrix* B) {
  return new Matrix(*A * *B);
}

dll Symmatrix* cb_symmatrix_new_rem(const Symmatrix* A, const Symmatrix* B) {
  return new Symmatrix(*A % *B);
}

dll Symmatrix* cb_symmatrix_new_plus(const Symmatrix* A, const Symmatrix* B) {
  return new Symmatrix(*A + *B);
}

dll Symmatrix* cb_symmatrix_new_minus2(const Symmatrix* A, const Symmatrix* B) {
  return new Symmatrix(*A - *B);
}

dll Matrix* cb_symmatrix_new_times2(const Symmatrix* A, const Matrix* B) {
  return new Matrix(*A * *B);
}

dll Matrix* cb_symmatrix_new_times3(const Matrix* A, const Symmatrix* B) {
  return new Matrix(*A * *B);
}

dll Matrix* cb_symmatrix_new_plus2(const Symmatrix* A, const Matrix* B) {
  return new Matrix(*A + *B);
}

dll Matrix* cb_symmatrix_new_plus3(const Matrix* A, const Symmatrix* B) {
  return new Matrix(*A + *B);
}

dll Matrix* cb_symmatrix_new_minus3(const Symmatrix* A, const Matrix* B) {
  return new Matrix(*A - *B);
}

dll Matrix* cb_symmatrix_new_minus4(const Matrix* A, const Symmatrix* B) {
  return new Matrix(*A - *B);
}

dll Symmatrix* cb_symmatrix_new_times4(const Symmatrix* A, Real d) {
  return new Symmatrix(*A * d);
}

dll Symmatrix* cb_symmatrix_new_times5(Real d, const Symmatrix* A) {
  return new Symmatrix(d * *A);
}

dll Symmatrix* cb_symmatrix_new_divide(const Symmatrix* A, Real d) {
  return new Symmatrix(*A / d);
}

dll Symmatrix* cb_symmatrix_new_plus4(const Symmatrix* A, Real d) {
  return new Symmatrix(*A + d);
}

dll Symmatrix* cb_symmatrix_new_plus5(Real d, const Symmatrix* A) {
  return new Symmatrix(d + *A);
}

dll Symmatrix* cb_symmatrix_new_minus5(const Symmatrix* A, Real d) {
  return new Symmatrix(*A - d);
}

dll Symmatrix* cb_symmatrix_new_minus6(Real d, const Symmatrix* A) {
  return new Symmatrix(d - *A);
}

dll Symmatrix* cb_symmatrix_new_transpose(const Symmatrix* A) {
  return new Symmatrix(transpose(*A));
}

dll Symmatrix* cb_symmatrix_xeya2(Symmatrix* self, const Matrix* A, Real d = 1.) {
  return &self->xeya(*A, d);
}

dll Symmatrix* cb_symmatrix_xpeya2(Symmatrix* self, const Matrix* A, Real d = 1.) {
  return &self->xpeya(*A, d);
}

dll Symmatrix* cb_symmatrix_xeya3(Symmatrix* self, const Indexmatrix* A, Real d = 1.) {
  return &self->xeya(*A, d);
}

dll Symmatrix* cb_symmatrix_xpeya3(Symmatrix* self, const Indexmatrix* A, Real d = 1.) {
  return &self->xpeya(*A, d);
}

dll Symmatrix* cb_symmatrix_xeya4(Symmatrix* self, const Sparsesym* A, Real d = 1.) {
  return &self->xeya(*A, d);
}

dll Symmatrix* cb_symmatrix_xpeya4(Symmatrix* self, const Sparsesym* A, Real d = 1.) {
  return &self->xpeya(*A, d);
}

dll Symmatrix* cb_symmatrix_xetriu_yza(Symmatrix* self, const Matrix* A, const Matrix* B, Real d = 1.) {
  return &self->xetriu_yza(*A, *B, d);
}

dll Symmatrix* cb_symmatrix_xpetriu_yza(Symmatrix* self, const Matrix* A, const Matrix* B, Real d = 1.) {
  return &self->xpetriu_yza(*A, *B, d);
}

dll Symmatrix* cb_symmatrix_xetriu_yza2(Symmatrix* self, const Sparsemat* A, const Matrix* B, Real d = 1.) {
  return &self->xetriu_yza(*A, *B, d);
}

dll Symmatrix* cb_symmatrix_xpetriu_yza2(Symmatrix* self, const Sparsemat* A, const Matrix* B, Real d = 1.) {
  return &self->xpetriu_yza(*A, *B, d);
}

dll Symmatrix* cb_symmatrix_assign2(Symmatrix* self, const Sparsesym* A) {
  return &(*self = *A);
}

dll Symmatrix* cb_symmatrix_plus3(Symmatrix* self, const Sparsesym* A) {
  return &(*self += *A);
}

dll Symmatrix* cb_symmatrix_minus3(Symmatrix* self, const Sparsesym* A) {
  return &(*self -= *A);
}

dll Matrix* cb_symmatrix_genmult3(const Symmatrix* A, const Sparsemat* B, Matrix* C, Real alpha, Real beta, int btrans) {
  return &genmult(*A, *B, *C, alpha, beta, btrans);
}

dll Matrix* cb_symmatrix_genmult4(const Sparsemat* A, const Symmatrix* B, Matrix* C, Real alpha, Real beta, int atrans) {
  return &genmult(*A, *B, *C, alpha, beta, atrans);
}

dll Symmatrix* cb_symmatrix_rankadd2(const Sparsemat* A, Symmatrix* C, Real alpha, Real beta, int trans) {
  return &rankadd(*A, *C, alpha, beta, trans);
}

dll Symmatrix* cb_symmatrix_scaledrankadd2(const Sparsemat* A, const Matrix* D, Symmatrix* C, Real alpha, Real beta, int trans) {
  return &scaledrankadd(*A, *D, *C, alpha, beta, trans);
}

dll Symmatrix* cb_symmatrix_rank2add2(const Sparsemat* A, const Matrix* B, Symmatrix* C, Real alpha, Real beta, int trans) {
  return &rank2add(*A, *B, *C, alpha, beta, trans);
}

dll Symmatrix* cb_symmatrix_new_abs(const Symmatrix* A) {
  return new Symmatrix(abs(*A));
}

dll Symmatrix* cb_symmatrix_shift_diag(Symmatrix* self, Real s) {
  return &self->shift_diag(s);
}

dll int cb_symmatrix_ldlfactor(Symmatrix* self, Real tol = 1e-10) {
  return self->LDLfactor(tol);
}

dll int cb_symmatrix_ldlsolve(const Symmatrix* self, Matrix* x) {
  return self->LDLsolve(*x);
}

dll int cb_symmatrix_ldlinverse(const Symmatrix* self, Symmatrix* S) {
  return self->LDLinverse(*S);
}

dll int cb_symmatrix_chol_factor(Symmatrix* self, Real tol = 1e-10) {
  return self->Chol_factor(tol);
}

dll int cb_symmatrix_chol_solve(const Symmatrix* self, Matrix* x) {
  return self->Chol_solve(*x);
}

dll int cb_symmatrix_chol_inverse(const Symmatrix* self, Symmatrix* S) {
  return self->Chol_inverse(*S);
}

dll int cb_symmatrix_chol_lsolve(const Symmatrix* self, Matrix* rhs) {
  return self->Chol_Lsolve(*rhs);
}

dll int cb_symmatrix_chol_ltsolve(const Symmatrix* self, Matrix* rhs) {
  return self->Chol_Ltsolve(*rhs);
}

dll int cb_symmatrix_chol_scaleli(const Symmatrix* self, Symmatrix* S) {
  return self->Chol_scaleLi(*S);
}

dll int cb_symmatrix_chol_scalelt(const Symmatrix* self, Symmatrix* S) {
  return self->Chol_scaleLt(*S);
}

dll int cb_symmatrix_chol_lmult(const Symmatrix* self, Matrix* rhs) {
  return self->Chol_Lmult(*rhs);
}

dll int cb_symmatrix_chol_ltmult(const Symmatrix* self, Matrix* rhs) {
  return self->Chol_Ltmult(*rhs);
}

dll int cb_symmatrix_chol_factor2(Symmatrix* self, Indexmatrix* piv, Real tol = 1e-10) {
  return self->Chol_factor(*piv, tol);
}

dll int cb_symmatrix_chol_solve2(const Symmatrix* self, Matrix* x, const Indexmatrix* piv) {
  return self->Chol_solve(*x, *piv);
}

dll int cb_symmatrix_chol_inverse2(const Symmatrix* self, Symmatrix* S, const Indexmatrix* piv) {
  return self->Chol_inverse(*S, *piv);
}

dll int cb_symmatrix_aasen_factor(Symmatrix* self, Indexmatrix* piv) {
  return self->Aasen_factor(*piv);
}

dll int cb_symmatrix_aasen_lsolve(const Symmatrix* self, Matrix* x) {
  return self->Aasen_Lsolve(*x);
}

dll int cb_symmatrix_aasen_ltsolve(const Symmatrix* self, Matrix* x) {
  return self->Aasen_Ltsolve(*x);
}

dll int cb_symmatrix_aasen_tridiagsolve(const Symmatrix* self, Matrix* x) {
  return self->Aasen_tridiagsolve(*x);
}

dll int cb_symmatrix_aasen_solve(const Symmatrix* self, Matrix* x, const Indexmatrix* piv) {
  return self->Aasen_solve(*x, *piv);
}

dll Integer cb_symmatrix_eig(const Symmatrix* self, Matrix* P, Matrix* d, int sort_non_decreasingly = 1) {
  return self->eig(*P, *d, (bool)sort_non_decreasingly);
}

dll Real cb_symmatrix_trace(const Symmatrix* A) {
  return trace(*A);
}

dll Real cb_symmatrix_ip(const Symmatrix* A, const Symmatrix* B) {
  return ip(*A, *B);
}

dll Real cb_symmatrix_ip2(const Matrix* A, const Symmatrix* B) {
  return ip(*A, *B);
}

dll Real cb_symmatrix_ip3(const Symmatrix* A, const Matrix* B) {
  return ip(*A, *B);
}

dll Real cb_symmatrix_norm2(const Symmatrix* A) {
  return norm2(*A);
}

dll Matrix* cb_symmatrix_new_sumrows(const Symmatrix* A) {
  return new Matrix(sumrows(*A));
}

dll Matrix* cb_symmatrix_new_sumcols(const Symmatrix* A) {
  return new Matrix(sumcols(*A));
}

dll Real cb_symmatrix_sum(const Symmatrix* A) {
  return sum(*A);
}

dll void cb_symmatrix_svec(const Symmatrix* A, Matrix* v, Real a, int add, Integer startindex_vec, Integer startindex_A, Integer blockdim) {
  svec(*A, *v, a, (bool)add, startindex_vec, startindex_A, blockdim);
}

dll Matrix* cb_symmatrix_new_svec(const Symmatrix* A) {
  return new Matrix(svec(*A));
}

dll void cb_symmatrix_sveci(const Matrix* v, Symmatrix* A, Real a, int add, Integer startindex_vec, Integer startindex_A, Integer blockdim) {
  sveci(*v, *A, a, (bool)add, startindex_vec, startindex_A, blockdim);
}

dll void cb_symmatrix_init_svec(Symmatrix* self, Integer nr, const Real* dp, Integer incr = 1, Real d = 1.) {
  self->init_svec(nr, dp, incr, d);
}

dll void cb_symmatrix_store_svec(const Symmatrix* self, Real* dp, Integer incr = 1, Real d = 1.) {
  self->store_svec(dp, incr, d);
}

dll Symmatrix* cb_symmatrix_new_skron(const Symmatrix* A, const Symmatrix* B, Real alpha, int add, Integer startindex_S) {
  return new Symmatrix(skron(*A, *B, alpha, (bool)add, startindex_S));
}

dll void cb_symmatrix_skron(const Symmatrix* A, const Symmatrix* B, Symmatrix* S, Real a, int add, Integer startindex_S) {
  skron(*A, *B, *S, a, (bool)add, startindex_S);
}

dll void cb_symmatrix_symscale(const Symmatrix* A, const Matrix* B, Symmatrix* S, Real a, Real b, int btrans) {
  symscale(*A, *B, *S, a, b, btrans);
}

dll Matrix* cb_symmatrix_new_minrows(const Symmatrix* A) {
  return new Matrix(minrows(*A));
}

dll Matrix* cb_symmatrix_new_mincols(const Symmatrix* A) {
  return new Matrix(mincols(*A));
}

dll Real cb_symmatrix_min(const Symmatrix* A) {
  return min(*A);
}

dll Matrix* cb_symmatrix_new_maxrows(const Symmatrix* A) {
  return new Matrix(maxrows(*A));
}

dll Matrix* cb_symmatrix_new_maxcols(const Symmatrix* A) {
  return new Matrix(maxcols(*A));
}

dll Real cb_symmatrix_max(const Symmatrix* A) {
  return max(*A);
}

dll void cb_symmatrix_display(const Symmatrix* self, int precision = 0, int width = 0, int screenwidth = 0) {
  self->display(std::cout, precision, width, screenwidth);
}

dll void cb_symmatrix_mfile_output(const Symmatrix* self, int precision = 16, int width = 0) {
  self->mfile_output(std::cout, precision, width);
}

