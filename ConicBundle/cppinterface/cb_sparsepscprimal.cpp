dll void cb_sparsepscprimal_destroy(SparsePSCPrimal* self) {
  delete self;
}

dll SparsePSCPrimal* cb_sparsepscprimal_new(const Sparsesym* sps, double factor = 1.) {
  return new SparsePSCPrimal(*sps, factor);
}

dll SparsePSCPrimal* cb_sparsepscprimal_new2(const SparsePSCPrimal* pr, double factor = 1.) {
  return new SparsePSCPrimal(*pr, factor);
}

dll SparsePSCPrimal* cb_sparsepscprimal_assign(SparsePSCPrimal* self, const Sparsesym* sdp) {
  return &(*self = *sdp);
}

dll int cb_sparsepscprimal_assign_gram_matrix(SparsePSCPrimal* self, const Matrix* P) {
  return self->assign_Gram_matrix(*P);
}

dll int cb_sparsepscprimal_aggregate_primal_data(SparsePSCPrimal* self, const PrimalData* it, double factor = 1.) {
  return self->aggregate_primal_data(*it, factor);
}

dll int cb_sparsepscprimal_aggregate_gram_matrix(SparsePSCPrimal* self, const Matrix* P, double factor = 1.) {
  return self->aggregate_Gram_matrix(*P, factor);
}

dll int cb_sparsepscprimal_scale_primal_data(SparsePSCPrimal* self, double factor) {
  return self->scale_primal_data(factor);
}

dll int cb_sparsepscprimal_primal_ip(const SparsePSCPrimal* self, Real* value, const SparseCoeffmatMatrix* A, Integer column) {
  return self->primal_ip(*value, *A, column);
}

