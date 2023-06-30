dll void cb_sumbundleparameters_destroy(SumBundleParameters* self) {
  delete self;
}

dll SumBundleParameters* cb_sumbundleparameters_new2(int modelsize, int bundlesize = 10, int updaterule = 0, int incr = -1) {
  return new SumBundleParameters(modelsize, bundlesize, updaterule, 0, incr);
}

dll SumBundleParameters* cb_sumbundleparameters_new3(const SumBundleParameters* sbp) {
  return new SumBundleParameters(*sbp);
}

dll BundleParameters* cb_sumbundleparameters_clone_bundleparameters(const SumBundleParameters* self) {
  return self->clone_BundleParameters();
}

dll int cb_sumbundleparameters_select_model(SumBundleParameters* self, Indexmatrix* model_indices, Integer cand_id, const Matrix* cand_y, const MinorantBundle* minorants, Integer aggr_index, Real model_maxviol, BundleProxObject* H, BundleModel::ModelUpdate* model_update) {
  return self->select_model(*model_indices, cand_id, *cand_y, *minorants, aggr_index, model_maxviol, *H, *model_update);
}

