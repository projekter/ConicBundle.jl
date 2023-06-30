dll void cb_nncmodelparameters_destroy(NNCModelParameters* self) {
  delete self;
}

dll NNCModelParameters* cb_nncmodelparameters_new2(int modelsize, int bundlesize = 10, int updaterule = 0, int incr = -1) {
  return new NNCModelParameters(modelsize, bundlesize, updaterule, 0, incr);
}

dll NNCModelParameters* cb_nncmodelparameters_new3(const BundleParameters* bp, int incr = -1) {
  return new NNCModelParameters(*bp, 0, incr);
}

dll NNCModelParameters* cb_nncmodelparameters_new4(const NNCModelParameters* sms) {
  return new NNCModelParameters(*sms);
}

dll BundleParameters* cb_nncmodelparameters_clone_bundleparameters(const NNCModelParameters* self) {
  return self->clone_BundleParameters();
}

dll int cb_nncmodelparameters_select_model(NNCModelParameters* self, MinorantBundle* model, Matrix* coefficients, Matrix* activity_indicators, const MinorantPointer* aggregate, const MinorantPointer* center_minorant, const MinorantBundle* cand_minorants, const MinorantBundle* old_minorants, MatrixFunctionOracle* oracle, int function_task, Real function_factor, BundleModel::ModelUpdate* model_update, Integer center_id, const Matrix* center_y, Integer cand_id, const Matrix* cand_y, Real model_maxviol, BundleProxObject* H) {
  return self->select_model(*model, *coefficients, *activity_indicators, *aggregate, *center_minorant, *cand_minorants, *old_minorants, oracle, (FunctionTask)function_task, function_factor, *model_update, center_id, *center_y, cand_id, *cand_y, model_maxviol, *H);
}

