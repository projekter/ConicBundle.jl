dll void cb_socmodelparameters_destroy(SOCModelParameters* self) {
  delete self;
}

dll SOCModelParameters* cb_socmodelparameters_new2(int modelsize, int bundlesize = 10, int updaterule = 0, int incr = -1) {
  return new SOCModelParameters(modelsize, bundlesize, updaterule, 0, incr);
}

dll SOCModelParameters* cb_socmodelparameters_new3(const BundleParameters* bp, int incr = -1) {
  return new SOCModelParameters(*bp, 0, incr);
}

dll SOCModelParameters* cb_socmodelparameters_new4(const SOCModelParameters* sms) {
  return new SOCModelParameters(*sms);
}

dll BundleParameters* cb_socmodelparameters_clone_bundleparameters(const SOCModelParameters* self) {
  return self->clone_BundleParameters();
}

dll int cb_socmodelparameters_select_model(SOCModelParameters* self, Matrix* modelvecs, const Matrix* aggrvec, Real cand_SOCval, const Matrix* cand_SOCvec, Real center_SOCval, const Matrix* center_SOCvec, const Matrix* SOCvecs, SOCOracle* oracle, int function_task, Real function_factor, BundleModel::ModelUpdate* model_update, Integer center_id, const Matrix* center_y, Integer cand_id, const Matrix* cand_y, Real model_maxviol, BundleProxObject* H) {
  return self->select_model(*modelvecs, *aggrvec, cand_SOCval, *cand_SOCvec, center_SOCval, *center_SOCvec, *SOCvecs, oracle, (FunctionTask)function_task, function_factor, *model_update, center_id, *center_y, cand_id, *cand_y, model_maxviol, *H);
}

