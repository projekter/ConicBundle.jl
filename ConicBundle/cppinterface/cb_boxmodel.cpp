dll void cb_boxmodel_destroy(BoxModel* self) {
  delete self;
}

dll void cb_boxmodel_clear(BoxModel* self) {
  self->clear();
}

dll BoxModel* cb_boxmodel_new(BoxOracle* fo, Real fun_factor = 1., int fun_task = ObjectiveFunction, int cbinc = -1) {
  return new BoxModel(fo, fun_factor, (FunctionTask)fun_task, 0, cbinc);
}

dll ModifiableOracleObject* cb_boxmodel_get_oracle_object(BoxModel* self) {
  return self->get_oracle_object();
}

dll Real cb_boxmodel_lb_function(BoxModel* self, Integer y_id, const Matrix* y) {
  return self->lb_function(y_id, *y);
}

dll BundleData* cb_boxmodel_get_data(BoxModel* self) {
  return self->get_data();
}

dll const BundleData* cb_boxmodel_get_data2(const BoxModel* self) {
  return self->get_data();
}

dll int cb_boxmodel_set_data(BoxModel* self, BundleData* bd) {
  return self->set_data(bd);
}

dll const PrimalData* cb_boxmodel_get_approximate_primal(const BoxModel* self) {
  return self->get_approximate_primal();
}

dll const PrimalData* cb_boxmodel_get_center_primal(const BoxModel* self) {
  return self->get_center_primal();
}

dll const PrimalData* cb_boxmodel_get_candidate_primal(const BoxModel* self) {
  return self->get_candidate_primal();
}

dll int cb_boxmodel_call_primal_extender(BoxModel* self, PrimalExtender* prex) {
  return self->call_primal_extender(*prex);
}

dll int cb_boxmodel_set_bundle_parameters(BoxModel* self, const BundleParameters* bp) {
  return self->set_bundle_parameters(*bp);
}

dll BundleParameters* cb_boxmodel_get_bundle_parameters(const BoxModel* self) {
  return self->get_bundle_parameters();
}

dll int cb_boxmodel_get_ret_code(const BoxModel* self) {
  return self->get_ret_code();
}

dll void cb_boxmodel_set_out(BoxModel* self, int pril = 1) {
  self->set_out(&std::cout, pril);
}

