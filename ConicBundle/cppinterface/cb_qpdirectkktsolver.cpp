dll void cb_qpdirectkktsolver_destroy(QPDirectKKTSolver* self) {
  delete self;
}

dll void cb_qpdirectkktsolver_clear(QPDirectKKTSolver* self) {
  self->clear();
}

dll QPDirectKKTSolver* cb_qpdirectkktsolver_new(int in_factorize_ABC = 0, int cbinc = -1) {
  return new QPDirectKKTSolver((bool)in_factorize_ABC, 0, cbinc);
}

dll int cb_qpdirectkktsolver_qpinit_kktdata(QPDirectKKTSolver* self, QPSolverProxObject* Hp, QPModelBlockObject* model, const Sparsemat* A, const Indexmatrix* eq_indices) {
  return self->QPinit_KKTdata(Hp, model, A, eq_indices);
}

dll int cb_qpdirectkktsolver_qpinit_kktsystem(QPDirectKKTSolver* self, const Matrix* KKTdiagx, const Matrix* KKTdiagy, Real Hfactor, Real prec, QPSolverParameters* params) {
  return self->QPinit_KKTsystem(*KKTdiagx, *KKTdiagy, Hfactor, prec, params);
}

dll int cb_qpdirectkktsolver_qpsolve_kktsystem(QPDirectKKTSolver* self, Matrix* solx, Matrix* soly, const Matrix* primalrhs, const Matrix* dualrhs, Real rhsmu, Real rhscorr, Real prec, QPSolverParameters* params) {
  return self->QPsolve_KKTsystem(*solx, *soly, *primalrhs, *dualrhs, rhsmu, rhscorr, prec, params);
}

