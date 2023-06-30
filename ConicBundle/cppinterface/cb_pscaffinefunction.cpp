dll void cb_pscaffinefunction_destroy(PSCAffineFunction* self) {
  delete self;
}

dll PSCAffineFunction* cb_pscaffinefunction_new(int incr = -1) {
  return new PSCAffineFunction(0, incr);
}

dll PSCAffineFunction* cb_pscaffinefunction_new2(const SparseCoeffmatMatrix* C, const SparseCoeffmatMatrix* opAt, PSCPrimal* generating_primal = 0, int incr = -1) {
  return new PSCAffineFunction(*C, *opAt, generating_primal, 0, incr);
}

dll void cb_pscaffinefunction_set_check_correctness(PSCAffineFunction* self, int chk) {
  self->set_check_correctness((bool)chk);
}

dll void cb_pscaffinefunction_set_max_ritzvecs(PSCAffineFunction* self, Integer maxv) {
  self->set_max_Ritzvecs(maxv);
}

dll Minorant* cb_pscaffinefunction_generate_minorant(PSCAffineFunction* self, const Matrix* P) {
  return self->generate_minorant(*P);
}

dll int cb_pscaffinefunction_svec_projection(PSCAffineFunction* self, Matrix* svec_offset, Matrix* svec_coeffs, const Matrix* P, const Indexmatrix* index_subset = 0) {
  return self->svec_projection(*svec_offset, *svec_coeffs, *P, index_subset);
}

dll int cb_pscaffinefunction_evaluate_projection(PSCAffineFunction* self, const Matrix* current_point, const Matrix* P, const double relprec, Matrix* projected_Ritz_vectors, Matrix* projected_Ritz_values) {
  return self->evaluate_projection(*current_point, *P, relprec, *projected_Ritz_vectors, *projected_Ritz_values);
}

dll int cb_pscaffinefunction_left_right_product(PSCAffineFunction* self, int i, const Matrix* E, const Matrix* F, Matrix* G) {
  return self->left_right_product(i, *E, *F, *G);
}

dll int cb_pscaffinefunction_check_correctness(const PSCAffineFunction* self) {
  return self->check_correctness();
}

dll const SparseCoeffmatMatrix* cb_pscaffinefunction_get_opat(PSCAffineFunction* self) {
  return &self->get_opAt();
}

dll const SparseCoeffmatMatrix* cb_pscaffinefunction_get_c(PSCAffineFunction* self) {
  return &self->get_C();
}

dll const PSCPrimal* cb_pscaffinefunction_get_generating_primal(PSCAffineFunction* self) {
  return self->get_generating_primal();
}

dll void cb_pscaffinefunction_set_out(PSCAffineFunction* self, int pril = 1) {
  self->set_out(&std::cout, pril);
}

dll void cb_pscaffinefunction_set_cbout(PSCAffineFunction* self, int incr = -1) {
  self->set_cbout(0, incr);
}

dll void cb_pscaffinefunction_print_problem_data(const PSCAffineFunction* self) {
  self->print_problem_data(std::cout);
}

dll void cb_pscaffinefunction_print_problem_data_to_mfile(const PSCAffineFunction* self, Integer blocknr) {
  self->print_problem_data_to_mfile(std::cout, blocknr);
}

