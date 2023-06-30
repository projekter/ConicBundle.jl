dll void cb_pscmodel_destroy(PSCModel* self) {
  delete self;
}

dll void cb_pscmodel_clear(PSCModel* self) {
  self->clear();
}

dll int cb_pscmodel_set_variable_metric_selection(PSCModel* self, VariableMetricSelection* vms = 0) {
  return self->set_variable_metric_selection(vms);
}

dll PSCModel* cb_pscmodel_new(PSCOracle* fo, Real fun_factor = 1., int fun_task = ObjectiveFunction, int cbinc = -1) {
  return new PSCModel(fo, fun_factor, (FunctionTask)fun_task, 0, cbinc);
}

dll ModifiableOracleObject* cb_pscmodel_get_oracle_object(PSCModel* self) {
  return self->get_oracle_object();
}

dll Real cb_pscmodel_lb_function(PSCModel* self, Integer y_id, const Matrix* y) {
  return self->lb_function(y_id, *y);
}

dll BundleData* cb_pscmodel_get_data(PSCModel* self) {
  return self->get_data();
}

dll const BundleData* cb_pscmodel_get_data2(const PSCModel* self) {
  return self->get_data();
}

dll int cb_pscmodel_set_data(PSCModel* self, BundleData* bd) {
  return self->set_data(bd);
}

dll const PrimalData* cb_pscmodel_get_approximate_primal(const PSCModel* self) {
  return self->get_approximate_primal();
}

dll const PrimalData* cb_pscmodel_get_center_primal(const PSCModel* self) {
  return self->get_center_primal();
}

dll const PrimalData* cb_pscmodel_get_candidate_primal(const PSCModel* self) {
  return self->get_candidate_primal();
}

dll int cb_pscmodel_call_primal_extender(PSCModel* self, PrimalExtender* prex) {
  return self->call_primal_extender(*prex);
}

dll int cb_pscmodel_set_bundle_parameters(PSCModel* self, const BundleParameters* bp) {
  return self->set_bundle_parameters(*bp);
}

dll const BundleParameters* cb_pscmodel_get_bundle_parameters(const PSCModel* self) {
  return self->get_bundle_parameters();
}

dll int cb_pscmodel_get_ret_code(const PSCModel* self) {
  return self->get_ret_code();
}

dll void cb_pscmodel_output_bundle_data(const PSCModel* self) {
  self->output_bundle_data(std::cout);
}

dll void cb_pscmodel_set_out(PSCModel* self, int pril = 1) {
  self->set_out(&std::cout, pril);
}

