dll void cb_pscmodelparameters_destroy(PSCModelParameters* self) {
  delete self;
}

dll PSCModelParameters* cb_pscmodelparameters_new2(int modelsize, int bundlesize = 10, int updaterule = 0, int incr = -1) {
  return new PSCModelParameters(modelsize, bundlesize, updaterule, 0, incr);
}

dll PSCModelParameters* cb_pscmodelparameters_new3(const BundleParameters* bp, int incr = -1) {
  return new PSCModelParameters(*bp, 0, incr);
}

dll PSCModelParameters* cb_pscmodelparameters_new4(const PSCModelParameters* sms) {
  return new PSCModelParameters(*sms);
}

dll BundleParameters* cb_pscmodelparameters_clone_bundleparameters(const PSCModelParameters* self) {
  return self->clone_BundleParameters();
}

dll int cb_pscmodelparameters_select_model(PSCModelParameters* self, Matrix* modelvecs, MinorantPointer* model_aggregate, Matrix* topvecs, Matrix* Ritz_values, Integer* activedim, Integer* keepsize, Integer* skippedsize, const Real primal_Ritzval, const Matrix* primaleigs, const Matrix* primalvecs, const MinorantPointer* primal_aggregate, Real primal_aggregate_coeff, Real growthrate, const Matrix* primalgrowth, const Matrix* dualgrowth, const Matrix* cand_Ritzvec, const Matrix* cand_Ritzval, PSCOracle* oracle, Integer modification_id, int function_task, Real function_factor, BundleModel::ModelUpdate* model_update, Integer center_id, const Matrix* center_y, Integer cand_id, const Matrix* cand_y, Real model_maxviol, Real diffval_center_aggregate, BundleProxObject* H) {
  return self->select_model(*modelvecs, *model_aggregate, *topvecs, *Ritz_values, *activedim, *keepsize, *skippedsize, primal_Ritzval, *primaleigs, *primalvecs, *primal_aggregate, primal_aggregate_coeff, growthrate, *primalgrowth, *dualgrowth, *cand_Ritzvec, *cand_Ritzval, oracle, modification_id, (FunctionTask)function_task, function_factor, *model_update, center_id, *center_y, cand_id, *cand_y, model_maxviol, diffval_center_aggregate, *H);
}

