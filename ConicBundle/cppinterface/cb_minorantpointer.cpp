dll void cb_minorantpointer_destroy(MinorantPointer* self) {
  delete self;
}

dll MinorantPointer* cb_minorantpointer_assign(MinorantPointer* self, const MinorantPointer* mp) {
  return &(*self = *mp);
}

dll void cb_minorantpointer_init(MinorantPointer* self, const MinorantPointer* mp, Real factor = 1., int enforce_copy = 0) {
  self->init(*mp, factor, (bool)enforce_copy);
}

dll void cb_minorantpointer_init2(MinorantPointer* self, Minorant* mp, Integer modification_id = -1, Real factor = 1.) {
  self->init(mp, modification_id, factor);
}

dll MinorantPointer* cb_minorantpointer_new() {
  return new MinorantPointer();
}

dll MinorantPointer* cb_minorantpointer_new2(MinorantUseData* in_md) {
  return new MinorantPointer(in_md);
}

dll MinorantPointer* cb_minorantpointer_new3(const MinorantPointer* mp) {
  return new MinorantPointer(*mp);
}

dll MinorantPointer* cb_minorantpointer_new4(const MinorantPointer* mp, Real factor) {
  return new MinorantPointer(*mp, factor);
}

dll MinorantPointer* cb_minorantpointer_new5(Minorant* mnrt, Integer modification_id, Real factor = 1.) {
  return new MinorantPointer(mnrt, modification_id, factor);
}

dll void cb_minorantpointer_clear(MinorantPointer* self) {
  self->clear();
}

dll int cb_minorantpointer_synchronize_ids(MinorantPointer* self, Integer new_modification_id, Integer new_center_id, Integer old_center_id, Integer new_cand_id, Integer old_cand_id, Integer new_prex_id = 0) {
  return self->synchronize_ids(new_modification_id, new_center_id, old_center_id, new_cand_id, old_cand_id, new_prex_id);
}

dll int cb_minorantpointer_empty(const MinorantPointer* self) {
  return self->empty();
}

dll int cb_minorantpointer_valid(const MinorantPointer* self) {
  return self->valid();
}

dll int cb_minorantpointer_zero(const MinorantPointer* self) {
  return self->zero();
}

dll int cb_minorantpointer_aggregate(const MinorantPointer* self) {
  return self->aggregate();
}

dll int cb_minorantpointer_one_user(const MinorantPointer* self) {
  return self->one_user();
}

dll const Minorant* cb_minorantpointer_get_minorant(const MinorantPointer* self) {
  return self->get_minorant();
}

dll const PrimalData* cb_minorantpointer_get_primal(MinorantPointer* self) {
  return self->get_primal();
}

dll Real cb_minorantpointer_offset(const MinorantPointer* self) {
  return self->offset();
}

dll Real cb_minorantpointer_coeff(const MinorantPointer* self, Integer i) {
  return self->coeff(i);
}

dll Integer cb_minorantpointer_nonzeros(const MinorantPointer* self) {
  return self->nonzeros();
}

dll int cb_minorantpointer_scale(MinorantPointer* self, Real val) {
  return self->scale(val);
}

dll int cb_minorantpointer_call_primal_extender(MinorantPointer* self, PrimalExtender* prex, Integer prex_id) {
  return self->call_primal_extender(*prex, prex_id);
}

dll int cb_minorantpointer_apply_modification(MinorantPointer* self, const GroundsetModification* gsmdf, Integer mod_id, MinorantExtender* mex, int apply_gsmdf_costs = 0) {
  return self->apply_modification(*gsmdf, mod_id, mex, (bool)apply_gsmdf_costs);
}

dll Real cb_minorantpointer_evaluate(const MinorantPointer* self, Integer yid, const Matrix* y, int with_constant = 1) {
  return self->evaluate(yid, *y, (bool)with_constant);
}

dll int cb_minorantpointer_add_offset(MinorantPointer* self, Real offset) {
  return self->add_offset(offset);
}

dll int cb_minorantpointer_get_minorant2(const MinorantPointer* self, Real* offset, Matrix* mat, Integer column, Real alpha = 1., int add = 0, const Indexmatrix* skip_fixed = 0, const Matrix* fixed_vals = 0) {
  return self->get_minorant(*offset, *mat, column, alpha, (bool)add, skip_fixed, fixed_vals);
}

dll int cb_minorantpointer_get_minorant3(const MinorantPointer* self, MinorantPointer* mp, Real alpha = 1.) {
  return self->get_minorant(*mp, alpha);
}

dll int cb_minorantpointer_get_minorant4(const MinorantPointer* self, MinorantPointer* mp, Real alpha, const Sparsemat* sp, const Indexmatrix* provided_row_indices = 0, const Indexmatrix* needed_col_indices = 0, int enforce_copy = 0) {
  return self->get_minorant(*mp, alpha, sp, provided_row_indices, needed_col_indices, (bool)enforce_copy);
}

dll int cb_minorantpointer_aggregate2(MinorantPointer* self, const MinorantBundle* minorants, const Matrix* coeff, Real factor = 1.) {
  return self->aggregate(*minorants, *coeff, factor);
}

dll int cb_minorantpointer_aggregate3(MinorantPointer* self, const MinorantPointer* minorant, double itsfactor = 1.) {
  return self->aggregate(*minorant, itsfactor);
}

dll Real cb_minorantpointer_ip(const MinorantPointer* self, const MinorantPointer* mp, const Indexmatrix* skip_fixed = 0, const Matrix* ipdiag = 0) {
  return self->ip(*mp, skip_fixed, ipdiag);
}

dll Real cb_minorantpointer_ip2(const MinorantPointer* self, const Matrix* m, const Matrix* ipdiag = 0, Integer startindex_m = 0) {
  return self->ip(*m, ipdiag, startindex_m);
}

dll bool cb_minorantpointer_new_less(const MinorantPointer* self, const MinorantPointer* mp) {
  return self->operator<(*mp);
}

dll bool cb_minorantpointer_new_equal(const MinorantPointer* self, const MinorantPointer* mp) {
  return (*self == *mp);
}

dll bool cb_minorantpointer_new_greater(const MinorantPointer* self, const MinorantPointer* mp) {
  return self->operator>(*mp);
}

dll int cb_minorantpointer_equals(const MinorantPointer* self, const MinorantPointer* mp, Real tol = 1e-10) {
  return self->equals(*mp, tol);
}

dll Real cb_minorantpointer_norm_squared(const MinorantPointer* self, const Matrix* D = 0) {
  return self->norm_squared(D);
}

dll Real cb_minorantpointer_dual_norm_squared(const MinorantPointer* self, const Matrix* D = 0) {
  return self->dual_norm_squared(D);
}

dll Matrix* cb_minorantpointer_left_genmult(const MinorantPointer* self, const Matrix* B, Matrix* C, Real alpha = 1., Real beta = 0., int thistrans = 0, int btrans = 0, Integer thisindex = 0) {
  return &self->left_genmult(*B, *C, alpha, beta, thistrans, btrans, thisindex);
}

dll Matrix* cb_minorantpointer_right_genmult(const MinorantPointer* self, const Matrix* A, Matrix* C, Real alpha = 1., Real beta = 0., int atrans = 0, int thistrans = 0, Integer thisindex = 0) {
  return &self->right_genmult(*A, *C, alpha, beta, atrans, thistrans, thisindex);
}

dll void cb_minorantpointer_display(const MinorantPointer* self, int precision = 8) {
  self->display(std::cout, precision);
}

