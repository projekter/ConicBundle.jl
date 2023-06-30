dll void cb_indexmatrix_destroy(Indexmatrix* self) {
  delete self;
}

dll Indexmatrix* cb_indexmatrix_new() {
  return new Indexmatrix();
}

dll Indexmatrix* cb_indexmatrix_new2(const Indexmatrix* A, Integer d = 1) {
  return new Indexmatrix(*A, d);
}

dll Indexmatrix* cb_indexmatrix_new3(Integer param0_from, Integer param0_to, Integer param0_step = 1) {
  return new Indexmatrix(Range(param0_from, param0_to, param0_step));
}

dll Indexmatrix* cb_indexmatrix_new4(Integer nr, Integer nc) {
  return new Indexmatrix(nr, nc);
}

dll Indexmatrix* cb_indexmatrix_new5(Integer nr, Integer nc, Integer d) {
  return new Indexmatrix(nr, nc, d);
}

dll Indexmatrix* cb_indexmatrix_new6(Integer nr, Integer nc, const Integer* dp, Integer incr = 1) {
  return new Indexmatrix(nr, nc, dp, incr);
}

dll Indexmatrix* cb_indexmatrix_init(Indexmatrix* self, const Indexmatrix* A, Integer d = 1) {
  return &self->init(*A, d);
}

dll Indexmatrix* cb_indexmatrix_init2(Indexmatrix* self, Integer param0_from, Integer param0_to, Integer param0_step = 1) {
  return &self->init(Range(param0_from, param0_to, param0_step));
}

dll Indexmatrix* cb_indexmatrix_init3(Indexmatrix* self, Integer nr, Integer nc, Integer d) {
  return &self->init(nr, nc, d);
}

dll Indexmatrix* cb_indexmatrix_init4(Indexmatrix* self, Integer nr, Integer nc, const Integer* dp, Integer incr = 1) {
  return &self->init(nr, nc, dp, incr);
}

dll void cb_indexmatrix_newsize(Indexmatrix* self, Integer nr, Integer nc) {
  self->newsize(nr, nc);
}

dll Indexmatrix* cb_indexmatrix_new7(const Matrix* param0) {
  return new Indexmatrix(*param0);
}

dll Indexmatrix* cb_indexmatrix_new9(const Sparsemat* param0) {
  return new Indexmatrix(*param0);
}

dll void cb_indexmatrix_dim(const Indexmatrix* self, Integer* _nr, Integer* _nc) {
  self->dim(*_nr, *_nc);
}

dll Integer cb_indexmatrix_dim2(const Indexmatrix* self) {
  return self->dim();
}

dll Integer cb_indexmatrix_rowdim(const Indexmatrix* self) {
  return self->rowdim();
}

dll Integer cb_indexmatrix_coldim(const Indexmatrix* self) {
  return self->coldim();
}

dll void cb_indexmatrix_set(Indexmatrix* self, Integer i, Integer j, Integer value) {
  (*self)(i, j) = value;
}

dll void cb_indexmatrix_set2(Indexmatrix* self, Integer i, Integer value) {
  (*self)(i) = value;
}

dll Integer cb_indexmatrix_get(const Indexmatrix* self, Integer i, Integer j) {
  return (*self)(i, j);
}

dll Integer cb_indexmatrix_get2(const Indexmatrix* self, Integer i) {
  return (*self)(i);
}

dll Indexmatrix* cb_indexmatrix_new_get(const Indexmatrix* self, const Indexmatrix* vecrow, const Indexmatrix* veccol) {
  return new Indexmatrix((*self)(*vecrow, *veccol));
}

dll Indexmatrix* cb_indexmatrix_new_get2(const Indexmatrix* self, const Indexmatrix* A) {
  return new Indexmatrix((*self)(*A));
}

dll Indexmatrix* cb_indexmatrix_new_col(const Indexmatrix* self, Integer i) {
  return new Indexmatrix(self->col(i));
}

dll Indexmatrix* cb_indexmatrix_new_row(const Indexmatrix* self, Integer i) {
  return new Indexmatrix(self->row(i));
}

dll Indexmatrix* cb_indexmatrix_new_cols(const Indexmatrix* self, const Indexmatrix* vec) {
  return new Indexmatrix(self->cols(*vec));
}

dll Indexmatrix* cb_indexmatrix_new_rows(const Indexmatrix* self, const Indexmatrix* vec) {
  return new Indexmatrix(self->rows(*vec));
}

dll Indexmatrix* cb_indexmatrix_triu(Indexmatrix* self, Integer d = 0) {
  return &self->triu(d);
}

dll Indexmatrix* cb_indexmatrix_tril(Indexmatrix* self, Integer d = 0) {
  return &self->tril(d);
}

dll Indexmatrix* cb_indexmatrix_subassign(Indexmatrix* self, const Indexmatrix* vecrow, const Indexmatrix* veccol, const Indexmatrix* A) {
  return &self->subassign(*vecrow, *veccol, *A);
}

dll Indexmatrix* cb_indexmatrix_subassign2(Indexmatrix* self, const Indexmatrix* vec, const Indexmatrix* A) {
  return &self->subassign(*vec, *A);
}

dll Indexmatrix* cb_indexmatrix_delete_rows(Indexmatrix* self, const Indexmatrix* ind, int sorted_increasingly = 0) {
  return &self->delete_rows(*ind, (bool)sorted_increasingly);
}

dll Indexmatrix* cb_indexmatrix_delete_cols(Indexmatrix* self, const Indexmatrix* ind, int sorted_increasingly = 0) {
  return &self->delete_cols(*ind, (bool)sorted_increasingly);
}

dll Indexmatrix* cb_indexmatrix_insert_row(Indexmatrix* self, Integer i, const Indexmatrix* v) {
  return &self->insert_row(i, *v);
}

dll Indexmatrix* cb_indexmatrix_insert_col(Indexmatrix* self, Integer i, const Indexmatrix* v) {
  return &self->insert_col(i, *v);
}

dll Indexmatrix* cb_indexmatrix_reduce_length(Indexmatrix* self, Integer n) {
  return &self->reduce_length(n);
}

dll Indexmatrix* cb_indexmatrix_concat_right(Indexmatrix* self, const Indexmatrix* A) {
  return &self->concat_right(*A);
}

dll Indexmatrix* cb_indexmatrix_concat_below(Indexmatrix* self, const Indexmatrix* A) {
  return &self->concat_below(*A);
}

dll Indexmatrix* cb_indexmatrix_concat_below2(Indexmatrix* self, Integer d) {
  return &self->concat_below(d);
}

dll Indexmatrix* cb_indexmatrix_concat_right2(Indexmatrix* self, Integer d) {
  return &self->concat_right(d);
}

dll Indexmatrix* cb_indexmatrix_enlarge_right(Indexmatrix* self, Integer addnc) {
  return &self->enlarge_right(addnc);
}

dll Indexmatrix* cb_indexmatrix_enlarge_below(Indexmatrix* self, Integer addnr) {
  return &self->enlarge_below(addnr);
}

dll Indexmatrix* cb_indexmatrix_enlarge_right2(Indexmatrix* self, Integer addnc, Integer d) {
  return &self->enlarge_right(addnc, d);
}

dll Indexmatrix* cb_indexmatrix_enlarge_below2(Indexmatrix* self, Integer addnr, Integer d) {
  return &self->enlarge_below(addnr, d);
}

dll Indexmatrix* cb_indexmatrix_enlarge_right3(Indexmatrix* self, Integer addnc, const Integer* dp, Integer d = 1) {
  return &self->enlarge_right(addnc, dp, d);
}

dll Indexmatrix* cb_indexmatrix_enlarge_below3(Indexmatrix* self, Integer addnr, const Integer* dp, Integer d = 1) {
  return &self->enlarge_below(addnr, dp, d);
}

dll const Integer* cb_indexmatrix_get_store2(const Indexmatrix* self) {
  return self->get_store();
}

dll Indexmatrix* cb_indexmatrix_new_diag(const Indexmatrix* A) {
  return new Indexmatrix(diag(*A));
}

dll Indexmatrix* cb_indexmatrix_new_triu(const Indexmatrix* A, Integer d) {
  return new Indexmatrix(triu(*A, d));
}

dll Indexmatrix* cb_indexmatrix_new_tril(const Indexmatrix* A, Integer d) {
  return new Indexmatrix(tril(*A, d));
}

dll Indexmatrix* cb_indexmatrix_new_concat_right(const Indexmatrix* A, const Indexmatrix* B) {
  return new Indexmatrix(concat_right(*A, *B));
}

dll Indexmatrix* cb_indexmatrix_new_concat_below(const Indexmatrix* A, const Indexmatrix* B) {
  return new Indexmatrix(concat_below(*A, *B));
}

dll void cb_indexmatrix_swap(Indexmatrix* A, Indexmatrix* B) {
  swap(*A, *B);
}

dll Indexmatrix* cb_indexmatrix_xeya(Indexmatrix* self, const Indexmatrix* A, Integer d = 1) {
  return &self->xeya(*A, d);
}

dll Indexmatrix* cb_indexmatrix_xpeya(Indexmatrix* self, const Indexmatrix* A, Integer d = 1) {
  return &self->xpeya(*A, d);
}

dll Indexmatrix* cb_indexmatrix_xbpeya(Indexmatrix* x, const Indexmatrix* y, Integer alpha, Integer beta, int ytrans) {
  return &xbpeya(*x, *y, alpha, beta, ytrans);
}

dll Indexmatrix* cb_indexmatrix_xeyapzb(Indexmatrix* x, const Indexmatrix* y, const Indexmatrix* z, Integer alpha, Integer beta) {
  return &xeyapzb(*x, *y, *z, alpha, beta);
}

dll Indexmatrix* cb_indexmatrix_genmult(const Indexmatrix* A, const Indexmatrix* B, Indexmatrix* C, Integer alpha, Integer beta, int atrans, int btrans) {
  return &genmult(*A, *B, *C, alpha, beta, atrans, btrans);
}

dll Indexmatrix* cb_indexmatrix_assign(Indexmatrix* self, const Indexmatrix* A) {
  return &(*self = *A);
}

dll Indexmatrix* cb_indexmatrix_times(Indexmatrix* self, const Indexmatrix* s) {
  return &(*self *= *s);
}

dll Indexmatrix* cb_indexmatrix_plus(Indexmatrix* self, const Indexmatrix* v) {
  return &(*self += *v);
}

dll Indexmatrix* cb_indexmatrix_minus(Indexmatrix* self, const Indexmatrix* v) {
  return &(*self -= *v);
}

dll Indexmatrix* cb_indexmatrix_rem(Indexmatrix* self, const Indexmatrix* A) {
  return &(*self %= *A);
}

dll Indexmatrix* cb_indexmatrix_new_minus(const Indexmatrix* self) {
  return new Indexmatrix(-(*self));
}

dll Indexmatrix* cb_indexmatrix_times2(Indexmatrix* self, Integer d) {
  return &(*self *= d);
}

dll Indexmatrix* cb_indexmatrix_divide(Indexmatrix* self, Integer d) {
  return &(*self /= d);
}

dll Indexmatrix* cb_indexmatrix_rem2(Indexmatrix* self, Integer d) {
  return &(*self %= d);
}

dll Indexmatrix* cb_indexmatrix_plus2(Indexmatrix* self, Integer d) {
  return &(*self += d);
}

dll Indexmatrix* cb_indexmatrix_minus2(Indexmatrix* self, Integer d) {
  return &(*self -= d);
}

dll Indexmatrix* cb_indexmatrix_transpose(Indexmatrix* self) {
  return &self->transpose();
}

dll Indexmatrix* cb_indexmatrix_new_times(const Indexmatrix* A, const Indexmatrix* B) {
  return new Indexmatrix(*A * *B);
}

dll Indexmatrix* cb_indexmatrix_new_plus(const Indexmatrix* A, const Indexmatrix* B) {
  return new Indexmatrix(*A + *B);
}

dll Indexmatrix* cb_indexmatrix_new_minus2(const Indexmatrix* A, const Indexmatrix* B) {
  return new Indexmatrix(*A - *B);
}

dll Indexmatrix* cb_indexmatrix_new_rem(const Indexmatrix* A, const Indexmatrix* B) {
  return new Indexmatrix(*A % *B);
}

dll Indexmatrix* cb_indexmatrix_new_times2(const Indexmatrix* A, Integer d) {
  return new Indexmatrix(*A * d);
}

dll Indexmatrix* cb_indexmatrix_new_times3(Integer d, const Indexmatrix* A) {
  return new Indexmatrix(d * *A);
}

dll Indexmatrix* cb_indexmatrix_new_divide(const Indexmatrix* A, Integer d) {
  return new Indexmatrix(*A / d);
}

dll Indexmatrix* cb_indexmatrix_new_rem2(const Indexmatrix* A, Integer d) {
  return new Indexmatrix(*A % d);
}

dll Indexmatrix* cb_indexmatrix_new_plus2(const Indexmatrix* A, Integer d) {
  return new Indexmatrix(*A + d);
}

dll Indexmatrix* cb_indexmatrix_new_plus3(Integer d, const Indexmatrix* A) {
  return new Indexmatrix(d + *A);
}

dll Indexmatrix* cb_indexmatrix_new_minus3(const Indexmatrix* A, Integer d) {
  return new Indexmatrix(*A - d);
}

dll Indexmatrix* cb_indexmatrix_new_minus4(Integer d, const Indexmatrix* A) {
  return new Indexmatrix(d - *A);
}

dll Indexmatrix* cb_indexmatrix_new_transpose(const Indexmatrix* A) {
  return new Indexmatrix(transpose(*A));
}

dll Indexmatrix* cb_indexmatrix_rand(Indexmatrix* self, Integer nr, Integer nc, Integer lowerb, Integer upperb, CH_Tools::GB_rand* random_generator = 0) {
  return &self->rand(nr, nc, lowerb, upperb, random_generator);
}

dll Indexmatrix* cb_indexmatrix_shuffle(Indexmatrix* self, CH_Tools::GB_rand* random_generator = 0) {
  return &self->shuffle(random_generator);
}

dll Indexmatrix* cb_indexmatrix_sign(Indexmatrix* self) {
  return &self->sign();
}

dll Indexmatrix* cb_indexmatrix_abs(Indexmatrix* self) {
  return &self->abs();
}

dll Indexmatrix* cb_indexmatrix_new_rand(Integer nr, Integer nc, Integer lb, Integer ub, CH_Tools::GB_rand* random_generator) {
  return new Indexmatrix(rand(nr, nc, lb, ub, random_generator));
}

dll Indexmatrix* cb_indexmatrix_new_sign(const Indexmatrix* A) {
  return new Indexmatrix(sign(*A));
}

dll Indexmatrix* cb_indexmatrix_new_abs(const Indexmatrix* A) {
  return new Indexmatrix(abs(*A));
}

dll Integer cb_indexmatrix_trace(const Indexmatrix* A) {
  return trace(*A);
}

dll Integer cb_indexmatrix_ip(const Indexmatrix* A, const Indexmatrix* B) {
  return ip(*A, *B);
}

dll Real cb_indexmatrix_norm2(const Indexmatrix* A) {
  return norm2(*A);
}

dll Indexmatrix* cb_indexmatrix_new_sumrows(const Indexmatrix* A) {
  return new Indexmatrix(sumrows(*A));
}

dll Indexmatrix* cb_indexmatrix_new_sumcols(const Indexmatrix* A) {
  return new Indexmatrix(sumcols(*A));
}

dll Integer cb_indexmatrix_sum(const Indexmatrix* A) {
  return sum(*A);
}

dll Indexmatrix* cb_indexmatrix_new_find(const Indexmatrix* self) {
  return new Indexmatrix(self->find());
}

dll Indexmatrix* cb_indexmatrix_new_find_number(const Indexmatrix* self, Integer num = 0) {
  return new Indexmatrix(self->find_number(num));
}

dll Indexmatrix* cb_indexmatrix_new_less(const Indexmatrix* A, const Indexmatrix* B) {
  return new Indexmatrix(*A < *B);
}

dll Indexmatrix* cb_indexmatrix_new_greater(const Indexmatrix* A, const Indexmatrix* B) {
  return new Indexmatrix(*A > *B);
}

dll Indexmatrix* cb_indexmatrix_new_lessequal(const Indexmatrix* A, const Indexmatrix* B) {
  return new Indexmatrix(*A <= *B);
}

dll Indexmatrix* cb_indexmatrix_new_greaterequal(const Indexmatrix* A, const Indexmatrix* B) {
  return new Indexmatrix(*A >= *B);
}

dll Indexmatrix* cb_indexmatrix_new_equal(const Indexmatrix* A, const Indexmatrix* B) {
  return new Indexmatrix(*A == *B);
}

dll Indexmatrix* cb_indexmatrix_new_inequal(const Indexmatrix* A, const Indexmatrix* B) {
  return new Indexmatrix(*A != *B);
}

dll Indexmatrix* cb_indexmatrix_new_less2(const Indexmatrix* A, Integer d) {
  return new Indexmatrix(*A < d);
}

dll Indexmatrix* cb_indexmatrix_new_greater2(const Indexmatrix* A, Integer d) {
  return new Indexmatrix(*A > d);
}

dll Indexmatrix* cb_indexmatrix_new_lessequal2(const Indexmatrix* A, Integer d) {
  return new Indexmatrix(*A <= d);
}

dll Indexmatrix* cb_indexmatrix_new_greaterequal2(const Indexmatrix* A, Integer d) {
  return new Indexmatrix(*A >= d);
}

dll Indexmatrix* cb_indexmatrix_new_equal2(const Indexmatrix* A, Integer d) {
  return new Indexmatrix(*A == d);
}

dll Indexmatrix* cb_indexmatrix_new_inequal2(const Indexmatrix* A, Integer d) {
  return new Indexmatrix(*A != d);
}

dll Indexmatrix* cb_indexmatrix_new_less3(Integer d, const Indexmatrix* A) {
  return new Indexmatrix(d < *A);
}

dll Indexmatrix* cb_indexmatrix_new_greater3(Integer d, const Indexmatrix* A) {
  return new Indexmatrix(d > *A);
}

dll Indexmatrix* cb_indexmatrix_new_lessequal3(Integer d, const Indexmatrix* A) {
  return new Indexmatrix(d <= *A);
}

dll Indexmatrix* cb_indexmatrix_new_greaterequal3(Integer d, const Indexmatrix* A) {
  return new Indexmatrix(d >= *A);
}

dll Indexmatrix* cb_indexmatrix_new_equal3(Integer d, const Indexmatrix* A) {
  return new Indexmatrix(d == *A);
}

dll Indexmatrix* cb_indexmatrix_new_inequal3(Integer d, const Indexmatrix* A) {
  return new Indexmatrix(d != *A);
}

dll int cb_indexmatrix_equal(const Indexmatrix* A, const Indexmatrix* b) {
  return equal(*A, *b);
}

dll Indexmatrix* cb_indexmatrix_new_minrows(const Indexmatrix* A) {
  return new Indexmatrix(minrows(*A));
}

dll Indexmatrix* cb_indexmatrix_new_mincols(const Indexmatrix* A) {
  return new Indexmatrix(mincols(*A));
}

dll Integer cb_indexmatrix_min(const Indexmatrix* A, Integer* iindex, Integer* jindex) {
  return min(*A, iindex, jindex);
}

dll Indexmatrix* cb_indexmatrix_new_maxrows(const Indexmatrix* A) {
  return new Indexmatrix(maxrows(*A));
}

dll Indexmatrix* cb_indexmatrix_new_maxcols(const Indexmatrix* A) {
  return new Indexmatrix(maxcols(*A));
}

dll Integer cb_indexmatrix_max(const Indexmatrix* A, Integer* iindex, Integer* jindex) {
  return max(*A, iindex, jindex);
}

dll Indexmatrix* cb_indexmatrix_new_sortindex(const Indexmatrix* vec, int nondecreasing) {
  return new Indexmatrix(sortindex(*vec, (bool)nondecreasing));
}

dll void cb_indexmatrix_sortindex(const Indexmatrix* vec, Indexmatrix* ind, int nondecreasing) {
  sortindex(*vec, *ind, (bool)nondecreasing);
}

dll void cb_indexmatrix_display(const Indexmatrix* self, int precision = 0, int width = 0, int screenwidth = 0) {
  self->display(std::cout, precision, width, screenwidth);
}

dll void cb_indexmatrix_mfile_output(const Indexmatrix* self, int precision = 16, int width = 0) {
  self->mfile_output(std::cout, precision, width);
}

