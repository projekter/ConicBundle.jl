dll void cb_pscaffineminorantextender_destroy(PSCAffineMinorantExtender* self) {
  delete self;
}

dll PSCAffineMinorantExtender* cb_pscaffineminorantextender_new(PSCAffineFunction* amf) {
  return new PSCAffineMinorantExtender(amf);
}

dll int cb_pscaffineminorantextender_extend(PSCAffineMinorantExtender* self, Minorant* minorant, int n_coords, const int* indices) {
  return self->extend(*minorant, n_coords, indices);
}

