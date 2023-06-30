dll void cb_sparsecoeffmatmatrix_destroy(SparseCoeffmatMatrix* self) {
  delete self;
}

dll SparseCoeffmatMatrix* cb_sparsecoeffmatmatrix_assign(SparseCoeffmatMatrix* self, const SparseCoeffmatMatrix* param0) {
  return &(*self = *param0);
}

dll void cb_sparsecoeffmatmatrix_clear(SparseCoeffmatMatrix* self) {
  self->clear();
}

dll int cb_sparsecoeffmatmatrix_init(SparseCoeffmatMatrix* self, const Indexmatrix* block_dim, Integer col_dim, const Indexmatrix* block_ind = 0, const Indexmatrix* col_ind = 0, const CoeffmatVector* coeff_vec = 0) {
  return self->init(*block_dim, col_dim, block_ind, col_ind, coeff_vec);
}

dll SparseCoeffmatMatrix* cb_sparsecoeffmatmatrix_new(int incr = -1) {
  return new SparseCoeffmatMatrix(0, incr);
}

dll SparseCoeffmatMatrix* cb_sparsecoeffmatmatrix_new2(const SparseCoeffmatMatrix* S, int incr = -1) {
  return new SparseCoeffmatMatrix(*S, 0, incr);
}

dll SparseCoeffmatMatrix* cb_sparsecoeffmatmatrix_new3(const Indexmatrix* in_block_dim, Integer in_col_dim, const Indexmatrix* block_ind = 0, const Indexmatrix* col_ind = 0, const CoeffmatVector* coeff_vec = 0, int incr = -1) {
  return new SparseCoeffmatMatrix(*in_block_dim, in_col_dim, block_ind, col_ind, coeff_vec, 0, incr);
}

dll const Indexmatrix* cb_sparsecoeffmatmatrix_blockdim(const SparseCoeffmatMatrix* self) {
  return &self->blockdim();
}

dll Integer cb_sparsecoeffmatmatrix_blockdim2(const SparseCoeffmatMatrix* self, Integer i) {
  return self->blockdim(i);
}

dll Integer cb_sparsecoeffmatmatrix_coldim(const SparseCoeffmatMatrix* self) {
  return self->coldim();
}

dll Integer cb_sparsecoeffmatmatrix_rowdim(const SparseCoeffmatMatrix* self) {
  return self->rowdim();
}

dll Integer cb_sparsecoeffmatmatrix_nzcoldim(const SparseCoeffmatMatrix* self) {
  return self->nzcoldim();
}

dll int cb_sparsecoeffmatmatrix_set(SparseCoeffmatMatrix* self, Integer i, Integer j, const CoeffmatPointer* cm) {
  return self->set(i, j, *cm);
}

dll int cb_sparsecoeffmatmatrix_set2(SparseCoeffmatMatrix* self, Integer i, Integer j, Coeffmat* cm) {
  return self->set(i, j, cm);
}

dll CoeffmatPointer* cb_sparsecoeffmatmatrix_new_get(const SparseCoeffmatMatrix* self, Integer i, Integer j) {
  return new CoeffmatPointer((*self)(i, j));
}

dll int cb_sparsecoeffmatmatrix_append_blocks(SparseCoeffmatMatrix* self, const SparseCoeffmatMatrix* append_mat, const Indexmatrix* blocks = 0, const Indexmatrix* cols = 0) {
  return self->append_blocks(*append_mat, blocks, cols);
}

dll int cb_sparsecoeffmatmatrix_append_columns(SparseCoeffmatMatrix* self, const SparseCoeffmatMatrix* append_mat, const Indexmatrix* blocks = 0, const Indexmatrix* cols = 0) {
  return self->append_columns(*append_mat, blocks, cols);
}

dll int cb_sparsecoeffmatmatrix_reassign_blocks(SparseCoeffmatMatrix* self, const Indexmatrix* map_to_old) {
  return self->reassign_blocks(*map_to_old);
}

dll int cb_sparsecoeffmatmatrix_reassign_columns(SparseCoeffmatMatrix* self, const Indexmatrix* map_to_old) {
  return self->reassign_columns(*map_to_old);
}

dll int cb_sparsecoeffmatmatrix_delete_blocks(SparseCoeffmatMatrix* self, const Indexmatrix* delete_indices, Indexmatrix* map_to_old = 0) {
  return self->delete_blocks(*delete_indices, map_to_old);
}

dll int cb_sparsecoeffmatmatrix_delete_columns(SparseCoeffmatMatrix* self, const Indexmatrix* delete_indices, Indexmatrix* map_to_old = 0) {
  return self->delete_columns(*delete_indices, map_to_old);
}

dll Integer cb_sparsecoeffmatmatrix_get_dense_cnt(const SparseCoeffmatMatrix* self, Integer i) {
  return self->get_dense_cnt(i);
}

dll bool cb_sparsecoeffmatmatrix_new_equal(const SparseCoeffmatMatrix* self, const SparseCoeffmatMatrix* mat) {
  return (*self == *mat);
}

dll bool cb_sparsecoeffmatmatrix_new_inequal(const SparseCoeffmatMatrix* self, const SparseCoeffmatMatrix* mat) {
  return (*self != *mat);
}

dll int cb_sparsecoeffmatmatrix_gram_ip(const SparseCoeffmatMatrix* self, Matrix* ipvec, const Matrix* P, const Matrix* Lam = 0, const Indexmatrix* ind = 0) {
  return self->Gram_ip(*ipvec, *P, Lam, ind);
}

dll int cb_sparsecoeffmatmatrix_gram_ip2(const SparseCoeffmatMatrix* self, Real* ipval, const Matrix* P, Integer j) {
  return self->Gram_ip(*ipval, *P, j);
}

dll int cb_sparsecoeffmatmatrix_primal_ip(const SparseCoeffmatMatrix* self, Matrix* ipvec, const PSCPrimal* primal, const Indexmatrix* ind = 0) {
  return self->primal_ip(*ipvec, primal, ind);
}

dll int cb_sparsecoeffmatmatrix_primal_ip2(const SparseCoeffmatMatrix* self, Real* value, const PSCPrimal* primal, Integer j) {
  return self->primal_ip(*value, primal, j);
}

dll int cb_sparsecoeffmatmatrix_project(const SparseCoeffmatMatrix* self, Symmatrix* S, const Matrix* P, const Integer j) {
  return self->project(*S, *P, j);
}

