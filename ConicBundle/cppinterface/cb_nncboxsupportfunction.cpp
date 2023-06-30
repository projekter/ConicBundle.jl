dll void cb_nncboxsupportfunction_destroy(NNCBoxSupportFunction* self) {
  delete self;
}

dll NNCBoxSupportFunction* cb_nncboxsupportfunction_new(const Matrix* lb, const Matrix* ub, int incr = -1) {
  return new NNCBoxSupportFunction(*lb, *ub, 0, incr);
}

dll int cb_nncboxsupportfunction_check_correctness(const NNCBoxSupportFunction* self) {
  return self->check_correctness();
}

dll const Matrix* cb_nncboxsupportfunction_get_lower_bounds(NNCBoxSupportFunction* self) {
  return &self->get_lower_bounds();
}

dll const Matrix* cb_nncboxsupportfunction_get_upper_bounds(NNCBoxSupportFunction* self) {
  return &self->get_upper_bounds();
}

dll void cb_nncboxsupportfunction_set_out(NNCBoxSupportFunction* self, int pril = 1) {
  self->set_out(&std::cout, pril);
}

dll void cb_nncboxsupportfunction_set_cbout(NNCBoxSupportFunction* self, int incr = -1) {
  self->set_cbout(0, incr);
}

dll void cb_nncboxsupportfunction_print_problem_data(const NNCBoxSupportFunction* self) {
  self->print_problem_data(std::cout);
}

dll void cb_nncboxsupportfunction_print_problem_data_to_mfile(const NNCBoxSupportFunction* self, Integer blocknr) {
  self->print_problem_data_to_mfile(std::cout, blocknr);
}

