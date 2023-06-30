dll void cb_socprimalextender_destroy(SOCPrimalExtender* self) {
  delete self;
}

dll int cb_socprimalextender_extend(SOCPrimalExtender* self, PrimalData* param0) {
  return self->extend(*param0);
}

dll int cb_socprimalextender_extend_soc(SOCPrimalExtender* self, Matrix* param0) {
  return self->extend_SOC(*param0);
}

