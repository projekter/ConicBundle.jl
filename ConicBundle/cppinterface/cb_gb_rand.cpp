dll void cb_gb_rand_destroy(GB_rand* self) {
  delete self;
}

dll void cb_gb_rand_init(GB_rand* self, long seed = 1) {
  self->init(seed);
}

dll GB_rand* cb_gb_rand_new(long seed = 1) {
  return new GB_rand(seed);
}

dll long cb_gb_rand_unif_long(GB_rand* self, long m) {
  return self->unif_long(m);
}

dll double cb_gb_rand_next(GB_rand* self) {
  return self->next();
}

dll void cb_gb_rand_save(const GB_rand* self) {
  self->save(std::cout);
}

