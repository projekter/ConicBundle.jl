dll void cb_aftmodification_destroy(AFTModification* self) {
  delete self;
}

dll AFTModification* cb_aftmodification_new(Integer var_olddim = 0, Integer row_olddim = 0, int ignore_groundset_modification = 0) {
  return new AFTModification(var_olddim, row_olddim, (bool)ignore_groundset_modification);
}

dll void cb_aftmodification_clear(AFTModification* self, Integer var_olddim, Integer row_olddim) {
  self->clear(var_olddim, row_olddim);
}

dll int cb_aftmodification_add_append_vars(AFTModification* self, Integer append_dim, const Sparsemat* append_cols, const Matrix* linear_costs) {
  return self->add_append_vars(append_dim, append_cols, linear_costs);
}

dll int cb_aftmodification_add_reassign_vars(AFTModification* self, const Indexmatrix* map_to_old) {
  return self->add_reassign_vars(*map_to_old);
}

dll int cb_aftmodification_add_delete_vars(AFTModification* self, const Indexmatrix* del_ind, Indexmatrix* map_to_old) {
  return self->add_delete_vars(*del_ind, *map_to_old);
}

dll int cb_aftmodification_add_append_rows(AFTModification* self, Integer append_dim, const Sparsemat* append_rows, const Matrix* append_rhs) {
  return self->add_append_rows(append_dim, append_rows, append_rhs);
}

dll int cb_aftmodification_add_reassign_rows(AFTModification* self, const Indexmatrix* map_to_old) {
  return self->add_reassign_rows(*map_to_old);
}

dll int cb_aftmodification_add_delete_rows(AFTModification* self, const Indexmatrix* rows_del_ind, Indexmatrix* rows_map_to_old) {
  return self->add_delete_rows(*rows_del_ind, *rows_map_to_old);
}

dll int cb_aftmodification_add_apply_factor(AFTModification* self, Real times_factor) {
  return self->add_apply_factor(times_factor);
}

dll int cb_aftmodification_add_offset(AFTModification* self, Real delta) {
  return self->add_offset(delta);
}

dll int cb_aftmodification_incorporate(AFTModification* self, const AFTModification* m) {
  return self->incorporate(*m);
}

dll int cb_aftmodification_apply_to_factor(const AFTModification* self, Real* f) {
  return self->apply_to_factor(*f);
}

dll int cb_aftmodification_apply_to_offset(const AFTModification* self, Real* o) {
  return self->apply_to_offset(*o);
}

dll const Matrix* cb_aftmodification_apply_modified_transform(const AFTModification* self, Matrix* out_y, const Matrix* in_y, const Sparsemat* arg_trafo, const Matrix* arg_offset) {
  return &self->apply_modified_transform(*out_y, *in_y, arg_trafo, arg_offset);
}

dll int cb_aftmodification_no_modification(const AFTModification* self) {
  return self->no_modification();
}

dll int cb_aftmodification_set_append_to_old(AFTModification* self, int append_only) {
  return self->set_append_to_old((bool)append_only);
}

dll int cb_aftmodification_append_to_old(const AFTModification* self) {
  return self->append_to_old();
}

dll int cb_aftmodification_only_scalars_change(const AFTModification* self) {
  return self->only_scalars_change();
}

dll int cb_aftmodification_ignore_groundset_modification(const AFTModification* self) {
  return self->ignore_groundset_modification();
}

dll int cb_aftmodification_groundset_changes_suffice_for_identity(AFTModification* self) {
  return self->groundset_changes_suffice_for_identity();
}

dll int cb_aftmodification_preserves_identity(const AFTModification* self) {
  return self->preserves_identity();
}

dll int cb_aftmodification_no_additions_or_deletions_in_vars(const AFTModification* self) {
  return self->no_additions_or_deletions_in_vars();
}

dll int cb_aftmodification_no_additions_or_deletions_in_rows(const AFTModification* self) {
  return self->no_additions_or_deletions_in_rows();
}

dll int cb_aftmodification_deleted_variables_are_zero(const AFTModification* self, const Matrix* oldpoint) {
  return self->deleted_variables_are_zero(*oldpoint);
}

dll int cb_aftmodification_new_variables_are_zero(const AFTModification* self, const Matrix* newpoint) {
  return self->new_variables_are_zero(*newpoint);
}

dll int cb_aftmodification_mapped_variables_are_equal(const AFTModification* self, const Matrix* newpoint, const Matrix* oldpoint) {
  return self->mapped_variables_are_equal(*newpoint, *oldpoint);
}

dll int cb_aftmodification_get_old_vardim(const AFTModification* self) {
  return self->get_old_vardim();
}

dll int cb_aftmodification_get_new_vardim(const AFTModification* self) {
  return self->get_new_vardim();
}

dll int cb_aftmodification_get_appended_vardim(const AFTModification* self) {
  return self->get_appended_vardim();
}

dll const int* cb_aftmodification_get_map_to_old_variables(const AFTModification* self) {
  return self->get_map_to_old_variables();
}

dll int cb_aftmodification_incorporate2(AFTModification* self, const OracleModification* m) {
  return self->incorporate(*m);
}

dll OracleModification* cb_aftmodification_new_initial_oraclemodification(const AFTModification* self, int old_var_dim) {
  return self->new_initial_oraclemodification(old_var_dim);
}

dll int cb_aftmodification_add_append_variables(AFTModification* self, int append_dim) {
  return self->add_append_variables(append_dim);
}

dll int cb_aftmodification_add_reassign_variables(AFTModification* self, int new_dim, const int* map_to_old_indices) {
  return self->add_reassign_variables(new_dim, map_to_old_indices);
}

dll Integer cb_aftmodification_old_vardim(const AFTModification* self) {
  return self->old_vardim();
}

dll Integer cb_aftmodification_new_vardim(const AFTModification* self) {
  return self->new_vardim();
}

dll Integer cb_aftmodification_appended_vardim(const AFTModification* self) {
  return self->appended_vardim();
}

dll Integer cb_aftmodification_old_rowdim(const AFTModification* self) {
  return self->old_rowdim();
}

dll Integer cb_aftmodification_new_rowdim(const AFTModification* self) {
  return self->new_rowdim();
}

dll Integer cb_aftmodification_appended_rowdim(const AFTModification* self) {
  return self->appended_rowdim();
}

dll const Indexmatrix* cb_aftmodification_map_to_old_variables(const AFTModification* self) {
  return self->map_to_old_variables();
}

dll const Indexmatrix* cb_aftmodification_deleted_var_indices(const AFTModification* self) {
  return self->deleted_var_indices();
}

dll const Indexmatrix* cb_aftmodification_new_var_indices(const AFTModification* self) {
  return self->new_var_indices();
}

dll const Indexmatrix* cb_aftmodification_map_to_old_rows(const AFTModification* self) {
  return self->map_to_old_rows();
}

dll const Indexmatrix* cb_aftmodification_deleted_row_indices(const AFTModification* self) {
  return self->deleted_row_indices();
}

dll const Indexmatrix* cb_aftmodification_new_row_indices(const AFTModification* self) {
  return self->new_row_indices();
}

dll Real cb_aftmodification_get_additional_offset(const AFTModification* self) {
  return self->get_additional_offset();
}

dll Real cb_aftmodification_get_additional_factor(const AFTModification* self) {
  return self->get_additional_factor();
}

dll const Sparsemat* cb_aftmodification_get_append_cols(const AFTModification* self) {
  return self->get_append_cols();
}

dll const Matrix* cb_aftmodification_get_append_costs(const AFTModification* self) {
  return self->get_append_costs();
}

dll const Sparsemat* cb_aftmodification_get_append_rows(const AFTModification* self) {
  return self->get_append_rows();
}

dll const Matrix* cb_aftmodification_get_append_rhs(const AFTModification* self) {
  return self->get_append_rhs();
}

dll void cb_aftmodification_set_out(AFTModification* self, int print_level = 1) {
  self->set_out(&std::cout, print_level);
}

