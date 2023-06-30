dll void cb_densepscprimal_destroy(DensePSCPrimal* self) {
  delete self;
}

dll DensePSCPrimal* cb_densepscprimal_new() {
  return new DensePSCPrimal();
}

dll DensePSCPrimal* cb_densepscprimal_new2(const DensePSCPrimal* symmat, double factor = 1.) {
  return new DensePSCPrimal(*symmat, factor);
}

dll DensePSCPrimal* cb_densepscprimal_new3(const Symmatrix* symmat, double factor = 1.) {
  return new DensePSCPrimal(*symmat, factor);
}

dll DensePSCPrimal* cb_densepscprimal_new4(Integer n) {
  return new DensePSCPrimal(n);
}

dll const DensePSCPrimal* cb_densepscprimal_assign(DensePSCPrimal* self, const Symmatrix* symmat) {
  return &(*self = *symmat);
}

dll int cb_densepscprimal_assign_gram_matrix(DensePSCPrimal* self, const Matrix* P) {
  return self->assign_Gram_matrix(*P);
}

dll int cb_densepscprimal_aggregate_primal_data(DensePSCPrimal* self, const PrimalData* it, double factor = 1.) {
  return self->aggregate_primal_data(*it, factor);
}

dll int cb_densepscprimal_aggregate_gram_matrix(DensePSCPrimal* self, const Matrix* P, double factor = 1.) {
  return self->aggregate_Gram_matrix(*P, factor);
}

dll int cb_densepscprimal_scale_primal_data(DensePSCPrimal* self, double factor) {
  return self->scale_primal_data(factor);
}

dll int cb_densepscprimal_primal_ip(const DensePSCPrimal* self, Real* value, const SparseCoeffmatMatrix* A, Integer column) {
  return self->primal_ip(*value, *A, column);
}

