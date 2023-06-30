dll void cb_nncboxsupportminorantextender_destroy(NNCBoxSupportMinorantExtender* self) {
  delete self;
}

dll NNCBoxSupportMinorantExtender* cb_nncboxsupportminorantextender_new(NNCBoxSupportFunction* fun) {
  return new NNCBoxSupportMinorantExtender(fun);
}

dll int cb_nncboxsupportminorantextender_extend(NNCBoxSupportMinorantExtender* self, Minorant* minorant, int n_coords, const int* indices) {
  return self->extend(*minorant, n_coords, indices);
}

