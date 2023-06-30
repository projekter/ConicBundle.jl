dll void cb_qpiterativekkthasolver_destroy(QPIterativeKKTHASolver* self) {
  delete self;
}

dll QPIterativeKKTHASolver* cb_qpiterativekkthasolver_new(IterativeSolverObject* insolver, QPKKTPrecondObject* inprecond = 0, int cbinc = -1) {
  return new QPIterativeKKTHASolver(insolver, inprecond, 0, cbinc);
}

dll int cb_qpiterativekkthasolver_qpsolve_kktsystem(QPIterativeKKTHASolver* self, Matrix* solx, Matrix* soly, const Matrix* primalrhs, const Matrix* dualrhs, Real rhsmu, Real rhscorr, Real prec, QPSolverParameters* params) {
  return self->QPsolve_KKTsystem(*solx, *soly, *primalrhs, *dualrhs, rhsmu, rhscorr, prec, params);
}

dll int cb_qpiterativekkthasolver_itsys_mult(QPIterativeKKTHASolver* self, const Matrix* in_vec, Matrix* out_vec) {
  return self->ItSys_mult(*in_vec, *out_vec);
}

dll Integer cb_qpiterativekkthasolver_qpget_system_size(QPIterativeKKTHASolver* self) {
  return self->QPget_system_size();
}

