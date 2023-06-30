dll void cb_qpkktsolvercomparison_destroy(QPKKTSolverComparison* self) {
  delete self;
}

dll void cb_qpkktsolvercomparison_clear(QPKKTSolverComparison* self) {
  self->clear();
}

dll QPKKTSolverComparison* cb_qpkktsolvercomparison_new(int cbinc = -1) {
  return new QPKKTSolverComparison(0, cbinc);
}

dll int cb_qpkktsolvercomparison_add_solver(QPKKTSolverComparison* self, QPKKTSolverObject* solver, const char* name) {
  return self->add_solver(solver, name);
}

dll int cb_qpkktsolvercomparison_qpinit_kktdata(QPKKTSolverComparison* self, QPSolverProxObject* Hp, QPModelBlockObject* model, const Sparsemat* A, const Indexmatrix* eq_indices) {
  return self->QPinit_KKTdata(Hp, model, A, eq_indices);
}

dll int cb_qpkktsolvercomparison_qpinit_kktsystem(QPKKTSolverComparison* self, const Matrix* KKTdiagx, const Matrix* KKTdiagy, Real Hfactor, Real prec, QPSolverParameters* params) {
  return self->QPinit_KKTsystem(*KKTdiagx, *KKTdiagy, Hfactor, prec, params);
}

dll int cb_qpkktsolvercomparison_qpsolve_kktsystem(QPKKTSolverComparison* self, Matrix* solx, Matrix* soly, const Matrix* primalrhs, const Matrix* dualrhs, Real rhsmu, Real rhscorr, Real prec, QPSolverParameters* params) {
  return self->QPsolve_KKTsystem(*solx, *soly, *primalrhs, *dualrhs, rhsmu, rhscorr, prec, params);
}

dll Real cb_qpkktsolvercomparison_qpget_blockh_norm(QPKKTSolverComparison* self) {
  return self->QPget_blockH_norm();
}

dll Real cb_qpkktsolvercomparison_qpget_blocka_norm(QPKKTSolverComparison* self) {
  return self->QPget_blockA_norm();
}

dll int cb_qpkktsolvercomparison_get_mu_stats(QPKKTSolverComparison* self, Real lbmu, Real ubmu, Indexmatrix* dims, Matrix* mu, Matrix* prepsecs, Matrix* predsecs, Matrix* corrsecs, Indexmatrix* predcalls, Indexmatrix* corrcalls, Matrix* cond, Indexmatrix* pccols, Matrix* sysviol) {
  return self->get_mu_stats(lbmu, ubmu, *dims, *mu, *prepsecs, *predsecs, *corrsecs, *predcalls, *corrcalls, *cond, *pccols, *sysviol);
}

dll int cb_qpkktsolvercomparison_get_prob_stats(QPKKTSolverComparison* self, Indexmatrix* dims, Indexmatrix* iterations, Matrix* lastmu, Matrix* prepsecs, Matrix* predsecs, Matrix* corrsecs, Indexmatrix* predcalls, Indexmatrix* corrcalls) {
  return self->get_prob_stats(*dims, *iterations, *lastmu, *prepsecs, *predsecs, *corrsecs, *predcalls, *corrcalls);
}

