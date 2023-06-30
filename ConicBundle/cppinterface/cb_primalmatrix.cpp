dll void cb_primalmatrix_destroy(PrimalMatrix* self) {
  delete self;
}

dll PrimalMatrix* cb_primalmatrix_new() {
  return new PrimalMatrix();
}

dll PrimalMatrix* cb_primalmatrix_new2(Integer nr, Integer nc) {
  return new PrimalMatrix(nr, nc);
}

dll PrimalMatrix* cb_primalmatrix_new3(Integer r, Integer c, Real d) {
  return new PrimalMatrix(r, c, d);
}

dll PrimalMatrix* cb_primalmatrix_new4(const PrimalMatrix* pm) {
  return new PrimalMatrix(*pm);
}

dll PrimalMatrix* cb_primalmatrix_new5(const Matrix* pm) {
  return new PrimalMatrix(*pm);
}

dll PrimalMatrix* cb_primalmatrix_assign(PrimalMatrix* self, const Matrix* pd) {
  return &(*self = *pd);
}

dll int cb_primalmatrix_aggregate_primal_data(PrimalMatrix* self, const PrimalData* it, double itsfactor) {
  return self->aggregate_primal_data(*it, itsfactor);
}

dll int cb_primalmatrix_scale_primal_data(PrimalMatrix* self, double myfactor) {
  return self->scale_primal_data(myfactor);
}

