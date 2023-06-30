dll void cb_boxprimalextender_destroy(BoxPrimalExtender* self) {
  delete self;
}

dll int cb_boxprimalextender_extend(BoxPrimalExtender* self, PrimalData* param0) {
  return self->extend(*param0);
}

dll int cb_boxprimalextender_extend_box(BoxPrimalExtender* self, Matrix* param0) {
  return self->extend_Box(*param0);
}

