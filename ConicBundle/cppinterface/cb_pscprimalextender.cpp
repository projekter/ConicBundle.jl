dll void cb_pscprimalextender_destroy(PSCPrimalExtender* self) {
  delete self;
}

dll int cb_pscprimalextender_extend(PSCPrimalExtender* self, PrimalData* param0) {
  return self->extend(*param0);
}

dll int cb_pscprimalextender_extend_ritz(PSCPrimalExtender* self, Matrix* param0) {
  return self->extend_Ritz(*param0);
}

