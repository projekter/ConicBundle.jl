dll void cb_socsupportfunction_destroy(SOCSupportFunction* self) {
  delete self;
}

dll SOCSupportFunction* cb_socsupportfunction_new(Integer socdim, int incr = -1) {
  return new SOCSupportFunction(socdim, 0, incr);
}

dll Minorant* cb_socsupportfunction_generate_minorant(SOCSupportFunction* self, const Matrix* SOCvec) {
  return self->generate_minorant(*SOCvec);
}

dll int cb_socsupportfunction_extract_socvector(SOCSupportFunction* self, Matrix* SOCvec, const Minorant* SOCminorant) {
  return self->extract_SOCvector(*SOCvec, SOCminorant);
}

dll int cb_socsupportfunction_projection(SOCSupportFunction* self, Matrix* offset, Matrix* coeffs, const Matrix* bar_P, const Indexmatrix* index_subset = 0) {
  return self->projection(*offset, *coeffs, *bar_P, index_subset);
}

dll int cb_socsupportfunction_evaluate_projection(SOCSupportFunction* self, const Matrix* current_point, const Matrix* P, const Real relprec, Real* projected_SOC_value) {
  return self->evaluate_projection(*current_point, *P, relprec, *projected_SOC_value);
}

dll int cb_socsupportfunction_check_correctness(const SOCSupportFunction* self) {
  return self->check_correctness();
}

dll Integer cb_socsupportfunction_get_socdim(SOCSupportFunction* self) {
  return self->get_socdim();
}

dll void cb_socsupportfunction_set_out(SOCSupportFunction* self, int pril = 1) {
  self->set_out(&std::cout, pril);
}

dll void cb_socsupportfunction_set_cbout(SOCSupportFunction* self, int incr = -1) {
  self->set_cbout(0, incr);
}

dll void cb_socsupportfunction_print_problem_data(const SOCSupportFunction* self) {
  self->print_problem_data(std::cout);
}

dll void cb_socsupportfunction_print_problem_data_to_mfile(const SOCSupportFunction* self, Integer blocknr) {
  self->print_problem_data_to_mfile(std::cout, blocknr);
}

