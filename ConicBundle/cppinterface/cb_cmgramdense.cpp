dll void cb_cmgramdense_destroy(CMgramdense* self) {
  delete self;
}

dll CMgramdense* cb_cmgramdense_new(const Matrix* Ain, int pos = 1, CoeffmatInfo* cip = 0) {
  return new CMgramdense(*Ain, (bool)pos, cip);
}

dll Coeffmat* cb_cmgramdense_clone(const CMgramdense* self) {
  return self->clone();
}

dll Integer cb_cmgramdense_dim(const CMgramdense* self) {
  return self->dim();
}

dll Real cb_cmgramdense_get(const CMgramdense* self, Integer i, Integer j) {
  return (*self)(i, j);
}

dll void cb_cmgramdense_make_symmatrix(const CMgramdense* self, Symmatrix* S) {
  self->make_symmatrix(*S);
}

dll Real cb_cmgramdense_norm(const CMgramdense* self) {
  return self->norm();
}

dll Coeffmat* cb_cmgramdense_subspace(const CMgramdense* self, const Matrix* P) {
  return self->subspace(*P);
}

dll void cb_cmgramdense_multiply(CMgramdense* self, Real d) {
  self->multiply(d);
}

dll Real cb_cmgramdense_ip(const CMgramdense* self, const Symmatrix* S) {
  return self->ip(*S);
}

dll Real cb_cmgramdense_gramip(const CMgramdense* self, const Matrix* P) {
  return self->gramip(*P);
}

dll Real cb_cmgramdense_gramip2(const CMgramdense* self, const Matrix* P, Integer start_row, const Matrix* Lam = 0) {
  return self->gramip(*P, start_row, Lam);
}

dll void cb_cmgramdense_addmeto(const CMgramdense* self, Symmatrix* S, Real d = 1.) {
  self->addmeto(*S, d);
}

dll void cb_cmgramdense_addprodto(const CMgramdense* self, Matrix* B, const Matrix* C, Real d = 1.) {
  self->addprodto(*B, *C, d);
}

dll void cb_cmgramdense_addprodto2(const CMgramdense* self, Matrix* B, const Sparsemat* C, Real d = 1.) {
  self->addprodto(*B, *C, d);
}

dll void cb_cmgramdense_left_right_prod(const CMgramdense* self, const Matrix* P, const Matrix* Q, Matrix* R) {
  self->left_right_prod(*P, *Q, *R);
}

dll Integer cb_cmgramdense_prodvec_flops(const CMgramdense* self) {
  return self->prodvec_flops();
}

dll int cb_cmgramdense_dense(const CMgramdense* self) {
  return self->dense();
}

dll int cb_cmgramdense_sparse(const CMgramdense* self) {
  return self->sparse();
}

dll int cb_cmgramdense_sparse2(const CMgramdense* self, Indexmatrix* param0, Indexmatrix* param1, Matrix* param2, Real param3) {
  return self->sparse(*param0, *param1, *param2, param3);
}

dll int cb_cmgramdense_support_in(const CMgramdense* self, const Sparsesym* param0) {
  return self->support_in(*param0);
}

dll Real cb_cmgramdense_ip2(const CMgramdense* self, const Sparsesym* S) {
  return self->ip(*S);
}

dll void cb_cmgramdense_project(const CMgramdense* self, Symmatrix* S, const Matrix* P) {
  self->project(*S, *P);
}

dll void cb_cmgramdense_add_projection(const CMgramdense* self, Symmatrix* S, const Matrix* P, Real alpha = 1., Integer start_row = 0) {
  self->add_projection(*S, *P, alpha, start_row);
}

dll const Matrix* cb_cmgramdense_postgenmult(const CMgramdense* self, const Matrix* B, Matrix* C, Real alpha = 1., Real beta = 0., int btrans = 0) {
  return &self->postgenmult(*B, *C, alpha, beta, btrans);
}

dll const Matrix* cb_cmgramdense_pregenmult(const CMgramdense* self, const Matrix* B, Matrix* C, Real alpha = 1., Real beta = 0., int btrans = 0) {
  return &self->pregenmult(*B, *C, alpha, beta, btrans);
}

dll int cb_cmgramdense_equal(const CMgramdense* self, const Coeffmat* p, double tol = 1e-6) {
  return self->equal(p, tol);
}

dll void cb_cmgramdense_display(const CMgramdense* self) {
  self->display(std::cout);
}

dll void cb_cmgramdense_out(const CMgramdense* self) {
  self->out(std::cout);
}

dll const Matrix* cb_cmgramdense_get_a(const CMgramdense* self) {
  return &self->get_A();
}

dll int cb_cmgramdense_get_positive(const CMgramdense* self) {
  return self->get_positive();
}

