dll void cb_variablemetricsvdselection_destroy(VariableMetricSVDSelection* self) {
  delete self;
}

dll VariableMetricSVDSelection* cb_variablemetricsvdselection_new(int cbincr = -1) {
  return new VariableMetricSVDSelection(0, cbincr);
}

dll VariableMetricSVDSelection* cb_variablemetricsvdselection_new2(Integer in_n_latest_minorants, Integer in_selection_method, Real in_oldfactor = 0., int cbincr = -1) {
  return new VariableMetricSVDSelection(in_n_latest_minorants, in_selection_method, in_oldfactor, 0, cbincr);
}

dll void cb_variablemetricsvdselection_get_uvlambda(VariableMetricSVDSelection* self, Matrix* U, Matrix* V, Matrix* lam, Matrix* cand) {
  self->get_UVlambda(*U, *V, *lam, *cand);
}

dll Integer cb_variablemetricsvdselection_get_n_latest_minorants(const VariableMetricSVDSelection* self) {
  return self->get_n_latest_minorants();
}

dll void cb_variablemetricsvdselection_set_n_latest_minorants(VariableMetricSVDSelection* self, Integer nlm) {
  self->set_n_latest_minorants(nlm);
}

dll Integer cb_variablemetricsvdselection_get_selection_method(const VariableMetricSVDSelection* self) {
  return self->get_selection_method();
}

dll void cb_variablemetricsvdselection_set_selection_method(VariableMetricSVDSelection* self, Integer sm) {
  self->set_selection_method(sm);
}

dll int cb_variablemetricsvdselection_add_variable_metric(VariableMetricSVDSelection* self, VariableMetric* H, Integer y_id, const Matrix* y, int descent_step, Real weightu, Real model_maxviol, const Indexmatrix* indices, VariableMetricBundleData* bundle_data) {
  return self->add_variable_metric(*H, y_id, *y, (bool)descent_step, weightu, model_maxviol, indices, *bundle_data);
}

dll VariableMetricSelection* cb_variablemetricsvdselection_clone_variablemetricselection(VariableMetricSVDSelection* self) {
  return self->clone_VariableMetricSelection();
}

