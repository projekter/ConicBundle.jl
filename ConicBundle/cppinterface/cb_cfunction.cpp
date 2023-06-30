dll void cb_cfunction_destroy(CFunction* self) {
  delete self;
}

dll CFunction* cb_cfunction_new(void* fk, cb_functionp fp, cb_subgextp se = 0, int prdim = 0) {
  return new CFunction(fk, fp, se, prdim);
}

dll void cb_cfunction_set_max_new(CFunction* self, Integer mn) {
  self->set_max_new(mn);
}

