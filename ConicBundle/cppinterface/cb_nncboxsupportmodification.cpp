dll void cb_nncboxsupportmodification_destroy(NNCBoxSupportModification* self) {
  delete self;
}

dll NNCBoxSupportModification* cb_nncboxsupportmodification_new(Integer var_olddim = 0, int incr = -1) {
  return new NNCBoxSupportModification(var_olddim, 0, incr);
}

dll void cb_nncboxsupportmodification_clear(NNCBoxSupportModification* self, Integer var_olddim) {
  self->clear(var_olddim);
}

dll int cb_nncboxsupportmodification_add_append_vars(NNCBoxSupportModification* self, Integer append_dim, const Matrix* append_lb, const Matrix* append_ub) {
  return self->add_append_vars(append_dim, append_lb, append_ub);
}

dll int cb_nncboxsupportmodification_add_reassign_vars(NNCBoxSupportModification* self, const Indexmatrix* map_to_old) {
  return self->add_reassign_vars(*map_to_old);
}

dll int cb_nncboxsupportmodification_add_delete_vars(NNCBoxSupportModification* self, const Indexmatrix* del_ind, Indexmatrix* map_to_old) {
  return self->add_delete_vars(*del_ind, *map_to_old);
}

dll int cb_nncboxsupportmodification_apply_to_bounds(const NNCBoxSupportModification* self, Matrix* lb, Matrix* ub) {
  return self->apply_to_bounds(*lb, *ub);
}

dll int cb_nncboxsupportmodification_no_modification(const NNCBoxSupportModification* self) {
  return self->no_modification();
}

dll int cb_nncboxsupportmodification_set_append_to_old(NNCBoxSupportModification* self, int append_only) {
  return self->set_append_to_old((bool)append_only);
}

dll int cb_nncboxsupportmodification_append_to_old(const NNCBoxSupportModification* self) {
  return self->append_to_old();
}

dll int cb_nncboxsupportmodification_no_additions_or_deletions_in_vars(const NNCBoxSupportModification* self) {
  return self->no_additions_or_deletions_in_vars();
}

dll int cb_nncboxsupportmodification_deleted_variables_are_zero(const NNCBoxSupportModification* self, const Matrix* oldpoint) {
  return self->deleted_variables_are_zero(*oldpoint);
}

dll int cb_nncboxsupportmodification_new_variables_are_zero(const NNCBoxSupportModification* self, const Matrix* newpoint) {
  return self->new_variables_are_zero(*newpoint);
}

dll int cb_nncboxsupportmodification_mapped_variables_are_equal(const NNCBoxSupportModification* self, const Matrix* newpoint, const Matrix* oldpoint) {
  return self->mapped_variables_are_equal(*newpoint, *oldpoint);
}

dll int cb_nncboxsupportmodification_get_old_vardim(const NNCBoxSupportModification* self) {
  return self->get_old_vardim();
}

dll int cb_nncboxsupportmodification_get_new_vardim(const NNCBoxSupportModification* self) {
  return self->get_new_vardim();
}

dll int cb_nncboxsupportmodification_get_appended_vardim(const NNCBoxSupportModification* self) {
  return self->get_appended_vardim();
}

dll const int* cb_nncboxsupportmodification_get_map_to_old_variables(const NNCBoxSupportModification* self) {
  return self->get_map_to_old_variables();
}

dll int cb_nncboxsupportmodification_incorporate(NNCBoxSupportModification* self, const OracleModification* m) {
  return self->incorporate(*m);
}

dll OracleModification* cb_nncboxsupportmodification_new_initial_oraclemodification(const NNCBoxSupportModification* self, int old_var_dim) {
  return self->new_initial_oraclemodification(old_var_dim);
}

dll int cb_nncboxsupportmodification_add_append_variables(NNCBoxSupportModification* self, int append_dim) {
  return self->add_append_variables(append_dim);
}

dll int cb_nncboxsupportmodification_add_reassign_variables(NNCBoxSupportModification* self, int new_dim, const int* map_to_old_indices) {
  return self->add_reassign_variables(new_dim, map_to_old_indices);
}

dll Integer cb_nncboxsupportmodification_old_vardim(const NNCBoxSupportModification* self) {
  return self->old_vardim();
}

dll Integer cb_nncboxsupportmodification_new_vardim(const NNCBoxSupportModification* self) {
  return self->new_vardim();
}

dll Integer cb_nncboxsupportmodification_appended_vardim(const NNCBoxSupportModification* self) {
  return self->appended_vardim();
}

dll const Indexmatrix* cb_nncboxsupportmodification_map_to_old_variables(const NNCBoxSupportModification* self) {
  return self->map_to_old_variables();
}

dll const Indexmatrix* cb_nncboxsupportmodification_deleted_var_indices(const NNCBoxSupportModification* self) {
  return self->deleted_var_indices();
}

dll const Indexmatrix* cb_nncboxsupportmodification_new_var_indices(const NNCBoxSupportModification* self) {
  return self->new_var_indices();
}

dll void cb_nncboxsupportmodification_set_out(NNCBoxSupportModification* self, int print_level = 1) {
  self->set_out(&std::cout, print_level);
}

dll void cb_nncboxsupportmodification_set_cbout(NNCBoxSupportModification* self, int incr = -1) {
  self->set_cbout(0, incr);
}

