dll void cb_socmodel_destroy(SOCModel* self) {
  delete self;
}

dll void cb_socmodel_clear(SOCModel* self) {
  self->clear();
}

dll SOCModel* cb_socmodel_new(SOCOracle* fo, Real fun_factor = 1., int fun_task = ObjectiveFunction, int cbinc = -1) {
  return new SOCModel(fo, fun_factor, (FunctionTask)fun_task, 0, cbinc);
}

dll ModifiableOracleObject* cb_socmodel_get_oracle_object(SOCModel* self) {
  return self->get_oracle_object();
}

dll Real cb_socmodel_lb_function(SOCModel* self, Integer y_id, const Matrix* y) {
  return self->lb_function(y_id, *y);
}

dll BundleData* cb_socmodel_get_data(SOCModel* self) {
  return self->get_data();
}

dll const BundleData* cb_socmodel_get_data2(const SOCModel* self) {
  return self->get_data();
}

dll int cb_socmodel_set_data(SOCModel* self, BundleData* bd) {
  return self->set_data(bd);
}

dll const PrimalData* cb_socmodel_get_approximate_primal(const SOCModel* self) {
  return self->get_approximate_primal();
}

dll const PrimalData* cb_socmodel_get_center_primal(const SOCModel* self) {
  return self->get_center_primal();
}

dll const PrimalData* cb_socmodel_get_candidate_primal(const SOCModel* self) {
  return self->get_candidate_primal();
}

dll int cb_socmodel_call_primal_extender(SOCModel* self, PrimalExtender* prex) {
  return self->call_primal_extender(*prex);
}

dll int cb_socmodel_set_bundle_parameters(SOCModel* self, const BundleParameters* bp) {
  return self->set_bundle_parameters(*bp);
}

dll const BundleParameters* cb_socmodel_get_bundle_parameters(const SOCModel* self) {
  return self->get_bundle_parameters();
}

dll int cb_socmodel_get_ret_code(const SOCModel* self) {
  return self->get_ret_code();
}

dll void cb_socmodel_output_bundle_data(const SOCModel* self) {
  self->output_bundle_data(std::cout);
}

dll void cb_socmodel_set_out(SOCModel* self, int pril = 1) {
  self->set_out(&std::cout, pril);
}

