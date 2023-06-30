dll void cb_socsupportminorantextender_destroy(SOCSupportMinorantExtender* self) {
  delete self;
}

dll SOCSupportMinorantExtender* cb_socsupportminorantextender_new(SOCSupportFunction* fun) {
  return new SOCSupportMinorantExtender(fun);
}

dll int cb_socsupportminorantextender_extend(SOCSupportMinorantExtender* self, Minorant* minorant, int n_coords, const int* indices) {
  return self->extend(*minorant, n_coords, indices);
}

