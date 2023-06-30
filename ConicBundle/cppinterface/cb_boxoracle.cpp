dll void cb_boxoracle_destroy(BoxOracle* self) {
  delete self;
}

dll BoxOracle* cb_boxoracle_new(const Matrix* in_lb, const Matrix* in_ub) {
  return new BoxOracle(*in_lb, *in_ub);
}

dll const Matrix* cb_boxoracle_get_lower_bounds(BoxOracle* self) {
  return &self->get_lower_bounds();
}

dll const Matrix* cb_boxoracle_get_upper_bounds(BoxOracle* self) {
  return &self->get_upper_bounds();
}

dll int cb_boxoracle_check_correctness(const BoxOracle* self) {
  return self->check_correctness();
}

