dll void cb_qpiterativekkthaeqsolver_destroy(QPIterativeKKTHAeqSolver* self) {
  delete self;
}

dll QPIterativeKKTHAeqSolver* cb_qpiterativekkthaeqsolver_new(IterativeSolverObject* insolver, QPKKTPrecondObject* inprecond = 0, int cbinc = -1) {
  return new QPIterativeKKTHAeqSolver(insolver, inprecond, 0, cbinc);
}

dll int cb_qpiterativekkthaeqsolver_qpinit_kktdata(QPIterativeKKTHAeqSolver* self, QPSolverProxObject* Hp, QPModelBlockObject* model, const Sparsemat* A, const Indexmatrix* eq_indices) {
  return self->QPinit_KKTdata(Hp, model, A, eq_indices);
}

dll int cb_qpiterativekkthaeqsolver_qpsolve_kktsystem(QPIterativeKKTHAeqSolver* self, Matrix* solx, Matrix* soly, const Matrix* primalrhs, const Matrix* dualrhs, Real rhsmu, Real rhscorr, Real prec, QPSolverParameters* params) {
  return self->QPsolve_KKTsystem(*solx, *soly, *primalrhs, *dualrhs, rhsmu, rhscorr, prec, params);
}

dll int cb_qpiterativekkthaeqsolver_itsys_mult(QPIterativeKKTHAeqSolver* self, const Matrix* in_vec, Matrix* out_vec) {
  return self->ItSys_mult(*in_vec, *out_vec);
}

dll Integer cb_qpiterativekkthaeqsolver_qpget_system_size(QPIterativeKKTHAeqSolver* self) {
  return self->QPget_system_size();
}

