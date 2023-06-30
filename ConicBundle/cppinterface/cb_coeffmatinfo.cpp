dll void cb_coeffmatinfo_destroy(CoeffmatInfo* self) {
  delete self;
}

dll CoeffmatInfo* cb_coeffmatinfo_new(Real sf = 1.) {
  return new CoeffmatInfo(sf);
}

dll CoeffmatInfo* cb_coeffmatinfo_clone(const CoeffmatInfo* self) {
  return self->clone();
}

dll Real cb_coeffmatinfo_get_scalefactor(const CoeffmatInfo* self) {
  return self->get_scalefactor();
}

dll void cb_coeffmatinfo_set_scalefactor(CoeffmatInfo* self, Real sf) {
  self->set_scalefactor(sf);
}

dll void cb_coeffmatinfo_multiply(CoeffmatInfo* self, Real sf) {
  self->multiply(sf);
}

dll void cb_coeffmatinfo_print_id(const CoeffmatInfo* self) {
  self->print_id(std::cout);
}

