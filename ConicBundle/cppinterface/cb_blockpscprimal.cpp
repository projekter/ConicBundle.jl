dll void cb_blockpscprimal_destroy(BlockPSCPrimal* self) {
  delete self;
}

dll BlockPSCPrimal* cb_blockpscprimal_new(const BlockPSCPrimal* pr, double factor = 1.) {
  return new BlockPSCPrimal(*pr, factor);
}

dll int cb_blockpscprimal_assign_gram_matrix(BlockPSCPrimal* self, const Matrix* P) {
  return self->assign_Gram_matrix(*P);
}

dll int cb_blockpscprimal_aggregate_primal_data(BlockPSCPrimal* self, const PrimalData* it, double factor = 1.) {
  return self->aggregate_primal_data(*it, factor);
}

dll int cb_blockpscprimal_aggregate_gram_matrix(BlockPSCPrimal* self, const Matrix* P, double factor = 1.) {
  return self->aggregate_Gram_matrix(*P, factor);
}

dll Integer cb_blockpscprimal_get_nblocks(const BlockPSCPrimal* self) {
  return self->get_nblocks();
}

dll Integer cb_blockpscprimal_blockdim(const BlockPSCPrimal* self, Integer i) {
  return self->blockdim(i);
}

dll PSCPrimal* cb_blockpscprimal_block(const BlockPSCPrimal* self, Integer i) {
  return self->block(i);
}

dll int cb_blockpscprimal_scale_primal_data(BlockPSCPrimal* self, double factor) {
  return self->scale_primal_data(factor);
}

dll int cb_blockpscprimal_primal_ip(const BlockPSCPrimal* self, Real* value, const SparseCoeffmatMatrix* A, Integer column) {
  return self->primal_ip(*value, *A, column);
}

