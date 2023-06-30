dll void cb_pscaffinemodification_destroy(PSCAffineModification* self) {
  delete self;
}

dll PSCAffineModification* cb_pscaffinemodification_new(Integer var_olddim, const Indexmatrix* block_olddim, int incr = 0) {
  return new PSCAffineModification(var_olddim, *block_olddim, 0, incr);
}

dll int cb_pscaffinemodification_clear(PSCAffineModification* self, Integer var_olddim, const Indexmatrix* block_olddim) {
  return self->clear(var_olddim, *block_olddim);
}

dll int cb_pscaffinemodification_add_append_vars(PSCAffineModification* self, Integer append_dim, const SparseCoeffmatMatrix* append_cols) {
  return self->add_append_vars(append_dim, append_cols);
}

dll int cb_pscaffinemodification_add_reassign_vars(PSCAffineModification* self, const Indexmatrix* map_to_old) {
  return self->add_reassign_vars(*map_to_old);
}

dll int cb_pscaffinemodification_add_delete_vars(PSCAffineModification* self, const Indexmatrix* del_ind, Indexmatrix* map_to_old) {
  return self->add_delete_vars(*del_ind, *map_to_old);
}

dll int cb_pscaffinemodification_add_append_blocks(PSCAffineModification* self, const Indexmatrix* append_dim, const SparseCoeffmatMatrix* append_offsets, const SparseCoeffmatMatrix* append_blocks) {
  return self->add_append_blocks(*append_dim, append_offsets, append_blocks);
}

dll int cb_pscaffinemodification_add_reassign_blocks(PSCAffineModification* self, const Indexmatrix* map_to_old) {
  return self->add_reassign_blocks(*map_to_old);
}

dll int cb_pscaffinemodification_add_delete_blocks(PSCAffineModification* self, const Indexmatrix* del_ind, Indexmatrix* map_to_old) {
  return self->add_delete_blocks(*del_ind, *map_to_old);
}

dll int cb_pscaffinemodification_add_reset_generating_primal(PSCAffineModification* self, PSCPrimal* new_generating_primal) {
  return self->add_reset_generating_primal(new_generating_primal);
}

dll int cb_pscaffinemodification_set_skip_extension(PSCAffineModification* self, int skip) {
  return self->set_skip_extension((bool)skip);
}

dll int cb_pscaffinemodification_apply_to_pscaffine(const PSCAffineModification* self, SparseCoeffmatMatrix* offset, SparseCoeffmatMatrix* matrix) {
  return self->apply_to_PSCAffine(offset, matrix);
}

dll int cb_pscaffinemodification_no_modification(const PSCAffineModification* self) {
  return self->no_modification();
}

dll int cb_pscaffinemodification_set_append_to_old(PSCAffineModification* self, int append_only) {
  return self->set_append_to_old((bool)append_only);
}

dll int cb_pscaffinemodification_append_to_old(const PSCAffineModification* self) {
  return self->append_to_old();
}

dll int cb_pscaffinemodification_deleted_variables_are_zero(const PSCAffineModification* self, const Matrix* oldpoint, const SparseCoeffmatMatrix* oldmat) {
  return self->deleted_variables_are_zero(*oldpoint, *oldmat);
}

dll int cb_pscaffinemodification_new_variables_are_zero(const PSCAffineModification* self, const Matrix* newpoint, const SparseCoeffmatMatrix* newmat) {
  return self->new_variables_are_zero(*newpoint, *newmat);
}

dll int cb_pscaffinemodification_mapped_variables_are_equal(const PSCAffineModification* self, const Matrix* newpoint, const Matrix* oldpoint) {
  return self->mapped_variables_are_equal(*newpoint, *oldpoint);
}

dll int cb_pscaffinemodification_variable_modifications(const PSCAffineModification* self) {
  return self->variable_modifications();
}

dll int cb_pscaffinemodification_block_modifications(const PSCAffineModification* self) {
  return self->block_modifications();
}

dll Integer cb_pscaffinemodification_old_vardim(const PSCAffineModification* self) {
  return self->old_vardim();
}

dll Integer cb_pscaffinemodification_new_vardim(const PSCAffineModification* self) {
  return self->new_vardim();
}

dll Integer cb_pscaffinemodification_appended_vardim(const PSCAffineModification* self) {
  return self->appended_vardim();
}

dll const Indexmatrix* cb_pscaffinemodification_old_blockdim(const PSCAffineModification* self) {
  return &self->old_blockdim();
}

dll const Indexmatrix* cb_pscaffinemodification_new_blockdim(const PSCAffineModification* self) {
  return &self->new_blockdim();
}

dll const Indexmatrix* cb_pscaffinemodification_appended_blockdim(const PSCAffineModification* self) {
  return &self->appended_blockdim();
}

dll const Indexmatrix* cb_pscaffinemodification_map_to_old_variables(const PSCAffineModification* self) {
  return self->map_to_old_variables();
}

dll const Indexmatrix* cb_pscaffinemodification_deleted_var_indices(const PSCAffineModification* self) {
  return self->deleted_var_indices();
}

dll const Indexmatrix* cb_pscaffinemodification_new_var_indices(const PSCAffineModification* self) {
  return self->new_var_indices();
}

dll const Indexmatrix* cb_pscaffinemodification_map_to_old_blocks(const PSCAffineModification* self) {
  return self->map_to_old_blocks();
}

dll const Indexmatrix* cb_pscaffinemodification_deleted_block_indices(const PSCAffineModification* self) {
  return self->deleted_block_indices();
}

dll const Indexmatrix* cb_pscaffinemodification_new_block_indices(const PSCAffineModification* self) {
  return self->new_block_indices();
}

dll const SparseCoeffmatMatrix* cb_pscaffinemodification_get_var_append(const PSCAffineModification* self) {
  return &self->get_var_append();
}

dll const SparseCoeffmatMatrix* cb_pscaffinemodification_get_block_append(const PSCAffineModification* self) {
  return &self->get_block_append();
}

dll const SparseCoeffmatMatrix* cb_pscaffinemodification_get_offset_append(const PSCAffineModification* self) {
  return &self->get_offset_append();
}

dll int cb_pscaffinemodification_get_reset_primal(const PSCAffineModification* self) {
  return self->get_reset_primal();
}

dll const PSCPrimal* cb_pscaffinemodification_get_generating_primal(const PSCAffineModification* self) {
  return self->get_generating_primal();
}

dll int cb_pscaffinemodification_get_skip_extension(const PSCAffineModification* self) {
  return self->get_skip_extension();
}

dll int cb_pscaffinemodification_get_old_vardim(const PSCAffineModification* self) {
  return self->get_old_vardim();
}

dll int cb_pscaffinemodification_get_new_vardim(const PSCAffineModification* self) {
  return self->get_new_vardim();
}

dll int cb_pscaffinemodification_get_appended_vardim(const PSCAffineModification* self) {
  return self->get_appended_vardim();
}

dll const int* cb_pscaffinemodification_get_map_to_old_variables(const PSCAffineModification* self) {
  return self->get_map_to_old_variables();
}

dll int cb_pscaffinemodification_incorporate(PSCAffineModification* self, const OracleModification* m) {
  return self->incorporate(*m);
}

dll OracleModification* cb_pscaffinemodification_new_initial_oraclemodification(const PSCAffineModification* self, int old_var_dim) {
  return self->new_initial_oraclemodification(old_var_dim);
}

dll int cb_pscaffinemodification_add_append_variables(PSCAffineModification* self, int append_dim) {
  return self->add_append_variables(append_dim);
}

dll int cb_pscaffinemodification_add_reassign_variables(PSCAffineModification* self, int new_dim, const int* map_to_old_indices) {
  return self->add_reassign_variables(new_dim, map_to_old_indices);
}

dll void cb_pscaffinemodification_set_cbout(PSCAffineModification* self, int incr = -1) {
  self->set_cbout(0, incr);
}

dll void cb_pscaffinemodification_set_out(PSCAffineModification* self, int print_level = 1) {
  self->set_out(&std::cout, print_level);
}

