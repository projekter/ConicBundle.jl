dll void cb_cfunctionminorantextender_destroy(CFunctionMinorantExtender* self) {
  delete self;
}

dll CFunctionMinorantExtender* cb_cfunctionminorantextender_new(void* fk, cb_subgextp se) {
  return new CFunctionMinorantExtender(fk, se);
}

dll int cb_cfunctionminorantextender_extend(CFunctionMinorantExtender* self, Minorant* minorant, int n_coords, const int* indices) {
  return self->extend(*minorant, n_coords, indices);
}

