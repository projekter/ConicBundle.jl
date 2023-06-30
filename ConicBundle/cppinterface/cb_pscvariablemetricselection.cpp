dll void cb_pscvariablemetricselection_destroy(PSCVariableMetricSelection* self) {
  delete self;
}

dll PSCVariableMetricSelection* cb_pscvariablemetricselection_new(int cbincr = -1) {
  return new PSCVariableMetricSelection(0, cbincr);
}

dll Integer cb_pscvariablemetricselection_get_selection_method(const PSCVariableMetricSelection* self) {
  return self->get_selection_method();
}

dll void cb_pscvariablemetricselection_set_selection_method(PSCVariableMetricSelection* self, Integer sm) {
  self->set_selection_method(sm);
}

dll Real cb_pscvariablemetricselection_get_oldfactor(const PSCVariableMetricSelection* self) {
  return self->get_oldfactor();
}

dll void cb_pscvariablemetricselection_set_oldfactor(PSCVariableMetricSelection* self, Real of) {
  self->set_oldfactor(of);
}

dll Real cb_pscvariablemetricselection_get_maxeigval_factor(const PSCVariableMetricSelection* self) {
  return self->get_maxeigval_factor();
}

dll void cb_pscvariablemetricselection_set_maxeigval_factor(PSCVariableMetricSelection* self, Real ef) {
  self->set_maxeigval_factor(ef);
}

dll Real cb_pscvariablemetricselection_get_mineigval_factor(const PSCVariableMetricSelection* self) {
  return self->get_mineigval_factor();
}

dll void cb_pscvariablemetricselection_set_mineigval_factor(PSCVariableMetricSelection* self, Real ef) {
  self->set_mineigval_factor(ef);
}

dll void cb_pscvariablemetricselection_set_oracle(PSCVariableMetricSelection* self, PSCOracle* psco) {
  self->set_oracle(psco);
}

dll int cb_pscvariablemetricselection_add_variable_metric(PSCVariableMetricSelection* self, VariableMetric* H, Integer y_id, const Matrix* y, int descent_step, Real weightu, Real model_maxviol, const Indexmatrix* indices, VariableMetricBundleData* bundle_data) {
  return self->add_variable_metric(*H, y_id, *y, (bool)descent_step, weightu, model_maxviol, indices, *bundle_data);
}

dll VariableMetricSelection* cb_pscvariablemetricselection_clone_variablemetricselection(PSCVariableMetricSelection* self) {
  return self->clone_VariableMetricSelection();
}

