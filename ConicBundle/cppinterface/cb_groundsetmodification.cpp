dll void cb_groundsetmodification_destroy(GroundsetModification* self) {
  delete self;
}

dll GroundsetModification* cb_groundsetmodification_new(Integer var_olddim = 0, int incr = -1) {
  return new GroundsetModification(var_olddim, 0, incr);
}

dll void cb_groundsetmodification_clear(GroundsetModification* self, Integer var_olddim) {
  self->clear(var_olddim);
}

dll int cb_groundsetmodification_add_append_vars(GroundsetModification* self, Integer append_dim, const Matrix* start_val = 0, const Matrix* costs = 0) {
  return self->add_append_vars(append_dim, start_val, costs);
}

dll int cb_groundsetmodification_add_reassign_vars(GroundsetModification* self, const Indexmatrix* map_to_old) {
  return self->add_reassign_vars(*map_to_old);
}

dll int cb_groundsetmodification_add_delete_vars(GroundsetModification* self, const Indexmatrix* del_ind, Indexmatrix* map_to_old) {
  return self->add_delete_vars(*del_ind, *map_to_old);
}

dll int cb_groundsetmodification_add_offset(GroundsetModification* self, Real delta) {
  return self->add_offset(delta);
}

dll int cb_groundsetmodification_apply_to_vars(const GroundsetModification* self, Matrix* vars) {
  return self->apply_to_vars(*vars);
}

dll int cb_groundsetmodification_apply_to_costs(const GroundsetModification* self, Matrix* costs, Real* in_offset) {
  return self->apply_to_costs(*costs, *in_offset);
}

dll int cb_groundsetmodification_no_modification(const GroundsetModification* self) {
  return self->no_modification();
}

dll int cb_groundsetmodification_set_append_to_old(GroundsetModification* self, int append_only) {
  return self->set_append_to_old((bool)append_only);
}

dll int cb_groundsetmodification_append_to_old(const GroundsetModification* self) {
  return self->append_to_old();
}

dll int cb_groundsetmodification_no_additions_or_deletions_in_vars(const GroundsetModification* self) {
  return self->no_additions_or_deletions_in_vars();
}

dll int cb_groundsetmodification_deleted_variables_are_zero(const GroundsetModification* self, const Matrix* oldpoint) {
  return self->deleted_variables_are_zero(*oldpoint);
}

dll int cb_groundsetmodification_new_variables_are_zero(const GroundsetModification* self, const Matrix* newpoint) {
  return self->new_variables_are_zero(*newpoint);
}

dll int cb_groundsetmodification_mapped_variables_are_equal(const GroundsetModification* self, const Matrix* newpoint, const Matrix* oldpoint) {
  return self->mapped_variables_are_equal(*newpoint, *oldpoint);
}

dll int cb_groundsetmodification_get_old_vardim(const GroundsetModification* self) {
  return self->get_old_vardim();
}

dll int cb_groundsetmodification_get_new_vardim(const GroundsetModification* self) {
  return self->get_new_vardim();
}

dll int cb_groundsetmodification_get_appended_vardim(const GroundsetModification* self) {
  return self->get_appended_vardim();
}

dll const int* cb_groundsetmodification_get_map_to_old_variables(const GroundsetModification* self) {
  return self->get_map_to_old_variables();
}

dll Real cb_groundsetmodification_get_add_offset(const GroundsetModification* self) {
  return self->get_add_offset();
}

dll const Matrix* cb_groundsetmodification_get_append_costs(const GroundsetModification* self) {
  return self->get_append_costs();
}

dll int cb_groundsetmodification_incorporate(GroundsetModification* self, const OracleModification* m) {
  return self->incorporate(*m);
}

dll OracleModification* cb_groundsetmodification_new_initial_oraclemodification(const GroundsetModification* self, int old_var_dim) {
  return self->new_initial_oraclemodification(old_var_dim);
}

dll int cb_groundsetmodification_add_append_variables(GroundsetModification* self, int append_dim) {
  return self->add_append_variables(append_dim);
}

dll int cb_groundsetmodification_add_reassign_variables(GroundsetModification* self, int new_dim, const int* map_to_old_indices) {
  return self->add_reassign_variables(new_dim, map_to_old_indices);
}

dll Integer cb_groundsetmodification_old_vardim(const GroundsetModification* self) {
  return self->old_vardim();
}

dll Integer cb_groundsetmodification_new_vardim(const GroundsetModification* self) {
  return self->new_vardim();
}

dll Integer cb_groundsetmodification_appended_vardim(const GroundsetModification* self) {
  return self->appended_vardim();
}

dll const Indexmatrix* cb_groundsetmodification_map_to_old_variables(const GroundsetModification* self) {
  return self->map_to_old_variables();
}

dll const Indexmatrix* cb_groundsetmodification_deleted_var_indices(const GroundsetModification* self) {
  return self->deleted_var_indices();
}

dll const Indexmatrix* cb_groundsetmodification_new_var_indices(const GroundsetModification* self) {
  return self->new_var_indices();
}

dll void cb_groundsetmodification_set_out(GroundsetModification* self, int print_level = 1) {
  self->set_out(&std::cout, print_level);
}

dll void cb_groundsetmodification_set_cbout(GroundsetModification* self, int incr = -1) {
  self->set_cbout(0, incr);
}

