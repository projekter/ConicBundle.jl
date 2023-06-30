dll void cb_affinefunctiontransformation_destroy(AffineFunctionTransformation* self) {
  delete self;
}

dll int cb_affinefunctiontransformation_init(AffineFunctionTransformation* self, Real fun_coeff = 1., Real fun_offset = 0., Matrix* linear_cost = 0, Matrix* arg_offset = 0, Sparsemat* arg_trafo = 0, int model_calls_delete = 1) {
  return self->init(fun_coeff, fun_offset, linear_cost, arg_offset, arg_trafo, (bool)model_calls_delete);
}

dll AffineFunctionTransformation* cb_affinefunctiontransformation_new(Real in_fun_coeff = 1., Real in_fun_offset = 0., Matrix* in_linear_cost = 0, Matrix* in_arg_offset = 0, Sparsemat* in_arg_trafo = 0, int in_model_calls_delete = 1, int incr = -1) {
  return new AffineFunctionTransformation(in_fun_coeff, in_fun_offset, in_linear_cost, in_arg_offset, in_arg_trafo, (bool)in_model_calls_delete, 0, incr);
}

dll int cb_affinefunctiontransformation_get_model_calls_delete(AffineFunctionTransformation* self) {
  return self->get_model_calls_delete();
}

dll void cb_affinefunctiontransformation_set_model_calls_delete(AffineFunctionTransformation* self, int mcd) {
  self->set_model_calls_delete((bool)mcd);
}

dll int cb_affinefunctiontransformation_argument_changes(const AffineFunctionTransformation* self) {
  return self->argument_changes();
}

dll int cb_affinefunctiontransformation_sparse_argument_changes(const AffineFunctionTransformation* self) {
  return self->sparse_argument_changes();
}

dll int cb_affinefunctiontransformation_scaled_index_subset(const AffineFunctionTransformation* self, const Indexmatrix* col_ind, const Indexmatrix* row_ind) {
  return self->scaled_index_subset(col_ind, row_ind);
}

dll int cb_affinefunctiontransformation_scaled_index(const AffineFunctionTransformation* self, Integer* mapped_index, Real* coeff, Integer index) {
  return self->scaled_index(*mapped_index, *coeff, index);
}

dll Real cb_affinefunctiontransformation_get_fun_coeff(const AffineFunctionTransformation* self) {
  return self->get_fun_coeff();
}

dll Real* cb_affinefunctiontransformation_set_fun_coeff(AffineFunctionTransformation* self) {
  return &self->set_fun_coeff();
}

dll Real cb_affinefunctiontransformation_get_fun_offset(const AffineFunctionTransformation* self) {
  return self->get_fun_offset();
}

dll Real* cb_affinefunctiontransformation_set_fun_offset(AffineFunctionTransformation* self) {
  return &self->set_fun_offset();
}

dll const Matrix* cb_affinefunctiontransformation_get_linear_cost(const AffineFunctionTransformation* self) {
  return self->get_linear_cost();
}

dll const Matrix* cb_affinefunctiontransformation_get_arg_offset(const AffineFunctionTransformation* self) {
  return self->get_arg_offset();
}

dll const Sparsemat* cb_affinefunctiontransformation_get_arg_trafo(const AffineFunctionTransformation* self) {
  return self->get_arg_trafo();
}

dll Real cb_affinefunctiontransformation_get_linear_cost2(const AffineFunctionTransformation* self, Integer i) {
  return self->get_linear_cost(i);
}

dll const MinorantPointer* cb_affinefunctiontransformation_get_constant_minorant(const AffineFunctionTransformation* self) {
  return &self->get_constant_minorant();
}

dll Integer cb_affinefunctiontransformation_from_dim(const AffineFunctionTransformation* self) {
  return self->from_dim();
}

dll Integer cb_affinefunctiontransformation_to_dim(const AffineFunctionTransformation* self) {
  return self->to_dim();
}

dll int cb_affinefunctiontransformation_copy_traforows(const AffineFunctionTransformation* self, Matrix* copy_to, const Matrix* copy_from) {
  return self->copy_traforows(*copy_to, *copy_from);
}

dll const Matrix* cb_affinefunctiontransformation_transform_argument(const AffineFunctionTransformation* self, Matrix* transformed_y, Real* transformed_offset, const Matrix* input_y) {
  return &self->transform_argument(*transformed_y, *transformed_offset, *input_y);
}

dll const Matrix* cb_affinefunctiontransformation_modified_transform_argument(const AffineFunctionTransformation* self, Matrix* transformed_y, const Matrix* input_y, const AFTModification* aftmdf, const GroundsetModification* gsmdf) {
  return &self->modified_transform_argument(*transformed_y, *input_y, aftmdf, *gsmdf);
}

dll Real cb_affinefunctiontransformation_objective_value(const AffineFunctionTransformation* self, Real offset, Real function_value) {
  return self->objective_value(offset, function_value);
}

dll int cb_affinefunctiontransformation_transform_minorant(const AffineFunctionTransformation* self, MinorantPointer* out_minorant, const MinorantPointer* in_minorant, Real alpha = 1., int add_trafo_minorant = 0, const Indexmatrix* provided_row_indices = 0, const Indexmatrix* needed_column_indices = 0) {
  return self->transform_minorant(*out_minorant, *in_minorant, alpha, (bool)add_trafo_minorant, provided_row_indices, needed_column_indices);
}

dll int cb_affinefunctiontransformation_transform_minorants(const AffineFunctionTransformation* self, MinorantBundle* out_minorants, const MinorantBundle* in_minorants, Real alpha = 1.) {
  return self->transform_minorants(*out_minorants, *in_minorants, alpha);
}

dll int cb_affinefunctiontransformation_qp_cost_indices(const AffineFunctionTransformation* self, Indexmatrix* provide_row_indices, const Indexmatrix* needed_col_indices) {
  return self->qp_cost_indices(*provide_row_indices, needed_col_indices);
}

dll int cb_affinefunctiontransformation_scaling_indices(const AffineFunctionTransformation* self, Indexmatrix* row_indices, const Indexmatrix* col_indices) {
  return self->scaling_indices(*row_indices, *col_indices);
}

dll int cb_affinefunctiontransformation_add_diagonal_scaling(const AffineFunctionTransformation* self, Matrix* diagscale, const Indexmatrix* indices, Real alpha, const Matrix* in_diagscale) {
  return self->add_diagonal_scaling(*diagscale, indices, alpha, *in_diagscale);
}

dll const GroundsetModification* cb_affinefunctiontransformation_analyze_modification(const AffineFunctionTransformation* self, int* minorant_trafo_differs, const AFTModification* aftmdf, const GroundsetModification* gsmdf) {
  minorant_trafo_differs = 0;
  return self->analyze_modification(*(bool*)minorant_trafo_differs, aftmdf, *gsmdf);
}

dll int cb_affinefunctiontransformation_apply_modification(AffineFunctionTransformation* self, const AFTModification* aftmdf, const GroundsetModification* gsmdf) {
  return self->apply_modification(aftmdf, *gsmdf);
}

dll void cb_affinefunctiontransformation_output_aft_data(const AffineFunctionTransformation* self) {
  self->output_aft_data(std::cout);
}

