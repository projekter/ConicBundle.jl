dll void cb_boxmodelparameters_destroy(BoxModelParameters* self) {
  delete self;
}

dll BoxModelParameters* cb_boxmodelparameters_new(int cbinc = -1) {
  return new BoxModelParameters(0, cbinc);
}

dll BoxModelParameters* cb_boxmodelparameters_new2(const BundleParameters* bp, int cbinc = -1) {
  return new BoxModelParameters(*bp, 0, cbinc);
}

dll BoxModelParameters* cb_boxmodelparameters_new3(const BoxModelParameters* sms) {
  return new BoxModelParameters(*sms);
}

dll BundleParameters* cb_boxmodelparameters_clone_bundleparameters(const BoxModelParameters* self) {
  return self->clone_BundleParameters();
}

dll int cb_boxmodelparameters_select_model(BoxModelParameters* self, MinorantBundle* box_model, Matrix* box_coeff, Matrix* box_indicators, Indexmatrix* box_coords, Matrix* box_complvalues, MinorantBundle* nnc_model, Matrix* nnc_coeff, Matrix* nnc_indicators, const Matrix* coord_switching, const MinorantBundle* minorants, const MinorantPointer* cand_minorant, const PrimalMatrix* cand_boxvec, const PrimalMatrix* aggr_boxvec, Real aggr_scaleval, BoxOracle* oracle, Integer modification_id, int function_task, Real function_factor, BundleModel::ModelUpdate* model_update, Integer center_id, const Matrix* center_y, Integer cand_id, const Matrix* cand_y, Real model_maxviol, BundleProxObject* H) {
  return self->select_model(*box_model, *box_coeff, *box_indicators, *box_coords, *box_complvalues, *nnc_model, *nnc_coeff, *nnc_indicators, *coord_switching, *minorants, *cand_minorant, *cand_boxvec, *aggr_boxvec, aggr_scaleval, oracle, modification_id, (FunctionTask)function_task, function_factor, *model_update, center_id, *center_y, cand_id, *cand_y, model_maxviol, *H);
}

