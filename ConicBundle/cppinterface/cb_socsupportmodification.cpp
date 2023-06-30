dll void cb_socsupportmodification_destroy(SOCSupportModification* self) {
  delete self;
}

dll SOCSupportModification* cb_socsupportmodification_new(Integer var_olddim = 0, int incr = -1) {
  return new SOCSupportModification(var_olddim, 0, incr);
}

dll void cb_socsupportmodification_clear(SOCSupportModification* self, Integer var_olddim) {
  self->clear(var_olddim);
}

dll int cb_socsupportmodification_add_append_vars(SOCSupportModification* self, Integer append_dim) {
  return self->add_append_vars(append_dim);
}

dll int cb_socsupportmodification_add_reassign_vars(SOCSupportModification* self, const Indexmatrix* map_to_old) {
  return self->add_reassign_vars(*map_to_old);
}

dll int cb_socsupportmodification_add_delete_vars(SOCSupportModification* self, const Indexmatrix* del_ind, Indexmatrix* map_to_old) {
  return self->add_delete_vars(*del_ind, *map_to_old);
}

dll int cb_socsupportmodification_no_modification(const SOCSupportModification* self) {
  return self->no_modification();
}

dll int cb_socsupportmodification_set_append_to_old(SOCSupportModification* self, int append_only) {
  return self->set_append_to_old((bool)append_only);
}

dll int cb_socsupportmodification_append_to_old(const SOCSupportModification* self) {
  return self->append_to_old();
}

dll int cb_socsupportmodification_no_additions_or_deletions_in_vars(const SOCSupportModification* self) {
  return self->no_additions_or_deletions_in_vars();
}

dll int cb_socsupportmodification_deleted_variables_are_zero(const SOCSupportModification* self, const Matrix* oldpoint) {
  return self->deleted_variables_are_zero(*oldpoint);
}

dll int cb_socsupportmodification_new_variables_are_zero(const SOCSupportModification* self, const Matrix* newpoint) {
  return self->new_variables_are_zero(*newpoint);
}

dll int cb_socsupportmodification_mapped_variables_are_equal(const SOCSupportModification* self, const Matrix* newpoint, const Matrix* oldpoint) {
  return self->mapped_variables_are_equal(*newpoint, *oldpoint);
}

dll int cb_socsupportmodification_get_old_vardim(const SOCSupportModification* self) {
  return self->get_old_vardim();
}

dll int cb_socsupportmodification_get_new_vardim(const SOCSupportModification* self) {
  return self->get_new_vardim();
}

dll int cb_socsupportmodification_get_appended_vardim(const SOCSupportModification* self) {
  return self->get_appended_vardim();
}

dll const int* cb_socsupportmodification_get_map_to_old_variables(const SOCSupportModification* self) {
  return self->get_map_to_old_variables();
}

dll int cb_socsupportmodification_incorporate(SOCSupportModification* self, const OracleModification* m) {
  return self->incorporate(*m);
}

dll OracleModification* cb_socsupportmodification_new_initial_oraclemodification(const SOCSupportModification* self, int old_var_dim) {
  return self->new_initial_oraclemodification(old_var_dim);
}

dll int cb_socsupportmodification_add_append_variables(SOCSupportModification* self, int append_dim) {
  return self->add_append_variables(append_dim);
}

dll int cb_socsupportmodification_add_reassign_variables(SOCSupportModification* self, int new_dim, const int* map_to_old_indices) {
  return self->add_reassign_variables(new_dim, map_to_old_indices);
}

dll Integer cb_socsupportmodification_old_vardim(const SOCSupportModification* self) {
  return self->old_vardim();
}

dll Integer cb_socsupportmodification_new_vardim(const SOCSupportModification* self) {
  return self->new_vardim();
}

dll Integer cb_socsupportmodification_appended_vardim(const SOCSupportModification* self) {
  return self->appended_vardim();
}

dll const Indexmatrix* cb_socsupportmodification_map_to_old_variables(const SOCSupportModification* self) {
  return self->map_to_old_variables();
}

dll const Indexmatrix* cb_socsupportmodification_deleted_var_indices(const SOCSupportModification* self) {
  return self->deleted_var_indices();
}

dll const Indexmatrix* cb_socsupportmodification_new_var_indices(const SOCSupportModification* self) {
  return self->new_var_indices();
}

dll void cb_socsupportmodification_set_out(SOCSupportModification* self, int print_level = 1) {
  self->set_out(&std::cout, print_level);
}

dll void cb_socsupportmodification_set_cbout(SOCSupportModification* self, int incr = -1) {
  self->set_cbout(0, incr);
}

