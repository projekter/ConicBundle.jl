dll void cb_minorant_destroy(Minorant* self) {
  delete self;
}

dll Minorant* cb_minorant_new(double offset, const DVector* subg, PrimalData* primal = 0, int offset_at_origin = 0) {
  return new Minorant(offset, *subg, primal, (bool)offset_at_origin);
}

dll Minorant* cb_minorant_new2(double offset, const DVector* subg_val, const IVector* subg_ind, PrimalData* primal = 0, int offset_at_origin = 0) {
  return new Minorant(offset, *subg_val, *subg_ind, primal, (bool)offset_at_origin);
}

dll Minorant* cb_minorant_new3(int offset_at_origin = 1, double offset = 0., int n_elementes = 0, const double* coeffs = 0, const int* indices = 0, double scale_val = 1., PrimalData* primal = 0) {
  return new Minorant((bool)offset_at_origin, offset, n_elementes, coeffs, indices, scale_val, primal);
}

dll Minorant* cb_minorant_new4(const Minorant* mnrt, double factor = 1., int with_primal = 0) {
  return new Minorant(mnrt, factor, (bool)with_primal);
}

dll double cb_minorant_offset(const Minorant* self) {
  return self->offset();
}

dll int cb_minorant_add_offset(Minorant* self, double value) {
  return self->add_offset(value);
}

dll bool* cb_minorant_offset_gives_value_at_origin(Minorant* self) {
  return &self->offset_gives_value_at_origin();
}

dll int cb_minorant_offset_gives_value_at_origin2(const Minorant* self) {
  return self->offset_gives_value_at_origin();
}

dll double cb_minorant_coeff(Minorant* self, int i) {
  return self->coeff(i);
}

dll int cb_minorant_add_coeff(Minorant* self, int i, double value) {
  return self->add_coeff(i, value);
}

dll int cb_minorant_add_coeffs(Minorant* self, int n_elements, const double* values, double factor = 1., int start_pos = 0) {
  return self->add_coeffs(n_elements, values, factor, start_pos);
}

dll int cb_minorant_add_coeffs2(Minorant* self, int n_elements, const double* values, const int* indices, double factor = 1.) {
  return self->add_coeffs(n_elements, values, indices, factor);
}

dll int cb_minorant_sparsify(Minorant* self, double tol = CB_minorant_zero_tolerance, double sparsity_ratio = CB_minorant_sparsity_ratio) {
  return self->sparsify(tol, sparsity_ratio);
}

dll int cb_minorant_nonzeros(Minorant* self) {
  return self->nonzeros();
}

dll double* cb_minorant_get_dense_coeff_store(Minorant* self, int n_elements) {
  return self->get_dense_coeff_store(n_elements);
}

dll int cb_minorant_reassign_coeffs(Minorant* self, int n_elements, const int* map_to_old_coeffs) {
  return self->reassign_coeffs(n_elements, map_to_old_coeffs);
}

dll void cb_minorant_set_primal(Minorant* self, PrimalData* param0) {
  self->set_primal(param0);
}

dll PrimalData* cb_minorant_get_primal(Minorant* self) {
  return self->get_primal();
}

dll const PrimalData* cb_minorant_get_primal2(const Minorant* self) {
  return self->get_primal();
}

dll Minorant* cb_minorant_clone_minorant(const Minorant* self, double factor = 1., int with_primal = 1) {
  return self->clone_minorant(factor, (bool)with_primal);
}

dll int cb_minorant_aggregate(Minorant* self, const Minorant* minorant, double factor = 1.) {
  return self->aggregate(*minorant, factor);
}

dll int cb_minorant_number_aggregated(const Minorant* self) {
  return self->number_aggregated();
}

dll int cb_minorant_scale_minorant(Minorant* self, double scale_val) {
  return self->scale_minorant(scale_val);
}

dll double cb_minorant_norm_squared(const Minorant* self) {
  return self->norm_squared();
}

