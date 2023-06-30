dll void cb_nncmodel_destroy(NNCModel* self) {
  delete self;
}

dll void cb_nncmodel_clear(NNCModel* self) {
  self->clear();
}

dll NNCModel* cb_nncmodel_new(MatrixFunctionOracle* fo, Real fun_factor = 1., int fun_task = ObjectiveFunction, int cbinc = -1) {
  return new NNCModel(fo, fun_factor, (FunctionTask)fun_task, 0, cbinc);
}

dll ModifiableOracleObject* cb_nncmodel_get_oracle_object(NNCModel* self) {
  return self->get_oracle_object();
}

dll Real cb_nncmodel_lb_function(NNCModel* self, Integer y_id, const Matrix* y) {
  return self->lb_function(y_id, *y);
}

dll BundleData* cb_nncmodel_get_data(NNCModel* self) {
  return self->get_data();
}

dll const BundleData* cb_nncmodel_get_data2(const NNCModel* self) {
  return self->get_data();
}

dll int cb_nncmodel_set_data(NNCModel* self, BundleData* bd) {
  return self->set_data(bd);
}

dll const PrimalData* cb_nncmodel_get_approximate_primal(const NNCModel* self) {
  return self->get_approximate_primal();
}

dll const PrimalData* cb_nncmodel_get_center_primal(const NNCModel* self) {
  return self->get_center_primal();
}

dll const PrimalData* cb_nncmodel_get_candidate_primal(const NNCModel* self) {
  return self->get_candidate_primal();
}

dll int cb_nncmodel_call_primal_extender(NNCModel* self, PrimalExtender* prex) {
  return self->call_primal_extender(*prex);
}

dll int cb_nncmodel_set_bundle_parameters(NNCModel* self, const BundleParameters* bp) {
  return self->set_bundle_parameters(*bp);
}

dll const BundleParameters* cb_nncmodel_get_bundle_parameters(const NNCModel* self) {
  return self->get_bundle_parameters();
}

dll int cb_nncmodel_get_ret_code(const NNCModel* self) {
  return self->get_ret_code();
}

dll void cb_nncmodel_set_out(NNCModel* self, int pril = 1) {
  self->set_out(&std::cout, pril);
}

