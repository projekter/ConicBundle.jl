dll void cb_gramsparsepscprimal_destroy(GramSparsePSCPrimal* self) {
  delete self;
}

dll GramSparsePSCPrimal* cb_gramsparsepscprimal_new(const Sparsesym* sps, double factor = 1.) {
  return new GramSparsePSCPrimal(*sps, factor);
}

dll GramSparsePSCPrimal* cb_gramsparsepscprimal_new2(const GramSparsePSCPrimal* pr, double factor = 1.) {
  return new GramSparsePSCPrimal(*pr, factor);
}

dll GramSparsePSCPrimal* cb_gramsparsepscprimal_assign(GramSparsePSCPrimal* self, const Sparsesym* sdp) {
  return &(*self = *sdp);
}

dll GramSparsePSCPrimal* cb_gramsparsepscprimal_assign2(GramSparsePSCPrimal* self, const GramSparsePSCPrimal* sdp) {
  return &(*self = *sdp);
}

dll const Matrix* cb_gramsparsepscprimal_get_grammatrix(const GramSparsePSCPrimal* self) {
  return &self->get_grammatrix();
}

dll int cb_gramsparsepscprimal_assign_gram_matrix(GramSparsePSCPrimal* self, const Matrix* P) {
  return self->assign_Gram_matrix(*P);
}

dll int cb_gramsparsepscprimal_aggregate_primal_data(GramSparsePSCPrimal* self, const PrimalData* it, double factor = 1.) {
  return self->aggregate_primal_data(*it, factor);
}

dll int cb_gramsparsepscprimal_aggregate_gram_matrix(GramSparsePSCPrimal* self, const Matrix* P, double factor = 1.) {
  return self->aggregate_Gram_matrix(*P, factor);
}

dll int cb_gramsparsepscprimal_scale_primal_data(GramSparsePSCPrimal* self, double factor) {
  return self->scale_primal_data(factor);
}

dll int cb_gramsparsepscprimal_primal_ip(const GramSparsePSCPrimal* self, Real* value, const SparseCoeffmatMatrix* A, Integer column) {
  return self->primal_ip(*value, *A, column);
}

