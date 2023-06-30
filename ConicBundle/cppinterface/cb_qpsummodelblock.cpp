dll void cb_qpsummodelblock_destroy(QPSumModelBlock* self) {
  delete self;
}

dll void cb_qpsummodelblock_clear(QPSumModelBlock* self) {
  self->clear();
}

dll QPSumModelBlock* cb_qpsummodelblock_new(int cbincr = -1) {
  return new QPSumModelBlock(0, cbincr);
}

dll QPModelBlockObject* cb_qpsummodelblock_clone(QPSumModelBlock* self) {
  return self->clone();
}

dll void cb_qpsummodelblock_recursive_delete_and_clear(QPSumModelBlock* self) {
  self->recursive_delete_and_clear();
}

dll int cb_qpsummodelblock_recursive_copy_data_of(QPSumModelBlock* self, QPModelBlockObject* param0) {
  return self->recursive_copy_data_of(param0);
}

dll int cb_qpsummodelblock_append(QPSumModelBlock* self, QPModelDataObject* inblock) {
  return self->append(inblock);
}

dll Matrix* cb_qpsummodelblock_b_times(QPSumModelBlock* self, const Matrix* A, Matrix* C, Real alpha, Real beta, int Btrans, int Atrans, Integer startindex_model, MinorantBundle* globalbundle, Integer startindex_bundle) {
  return &self->B_times(*A, *C, alpha, beta, Btrans, Atrans, startindex_model, *globalbundle, startindex_bundle);
}

dll Matrix* cb_qpsummodelblock_times_b(QPSumModelBlock* self, const Matrix* A, Matrix* C, Real alpha, Real beta, int Atrans, int Btrans, Integer startindex_model, MinorantBundle* globalbundle, Integer startindex_bundle) {
  return &self->times_B(*A, *C, alpha, beta, Atrans, Btrans, startindex_model, *globalbundle, startindex_bundle);
}

dll Symmatrix* cb_qpsummodelblock_add_bdbt(QPSumModelBlock* self, const Matrix* diagvec, Symmatrix* bigS, int minus, Integer startindex, Matrix* Bt, Integer startindex_model, MinorantBundle* globalbundle, Integer startindex_bundle) {
  return &self->add_BDBt(*diagvec, *bigS, (bool)minus, startindex, *Bt, startindex_model, *globalbundle, startindex_bundle);
}

dll Matrix* cb_qpsummodelblock_get_bt(QPSumModelBlock* self, Matrix* Bt, Integer startindex_model, MinorantBundle* global_bundle, Integer startindex_bundle) {
  return &self->get_Bt(*Bt, startindex_model, *global_bundle, startindex_bundle);
}

dll int cb_qpsummodelblock_get_modelx(QPSumModelBlock* self, Matrix* modelx, Integer startindex_model) {
  return self->get_modelx(*modelx, startindex_model);
}

dll int cb_qpsummodelblock_get_modeldx(QPSumModelBlock* self, Matrix* modeldx, Integer startindex_model) {
  return self->get_modeldx(*modeldx, startindex_model);
}

dll int cb_qpsummodelblock_get_modeldcstr(QPSumModelBlock* self, Matrix* modeldcstr, Integer startindex_constraints) {
  return self->get_modeldcstr(*modeldcstr, startindex_constraints);
}

dll int cb_qpsummodelblock_add_modelx_aggregate(QPSumModelBlock* self, Real* val, Matrix* vec, MinorantBundle* global_bundle, Integer startindex_bundle) {
  return self->add_modelx_aggregate(*val, *vec, *global_bundle, startindex_bundle);
}

dll int cb_qpsummodelblock_get_sysviol_model(QPSumModelBlock* self, Matrix* modelvec, Integer startindex_model, const Matrix* dy, MinorantBundle* global_bundle, Integer startindex_bundle) {
  return self->get_sysviol_model(*modelvec, startindex_model, *dy, *global_bundle, startindex_bundle);
}

dll int cb_qpsummodelblock_get_sysviol_constraints(QPSumModelBlock* self, Matrix* constrvec, Integer startindex_constr) {
  return self->get_sysviol_constraints(*constrvec, startindex_constr);
}

dll void cb_qpsummodelblock_display_model_values(QPSumModelBlock* self, const Matrix* y, MinorantBundle* global_bundle, Integer startindex_bundle) {
  self->display_model_values(*y, *global_bundle, startindex_bundle, std::cout);
}

dll int cb_qpsummodelblock_reset_starting_point(QPSumModelBlock* self, const Matrix* y, Real mu, MinorantBundle* global_bundle, Integer startindex_bundle) {
  return self->reset_starting_point(*y, mu, *global_bundle, startindex_bundle);
}

dll int cb_qpsummodelblock_compute_step(QPSumModelBlock* self, const Matrix* ystep, MinorantBundle* global_bundle, Integer startindex_bundle) {
  return self->compute_step(*ystep, *global_bundle, startindex_bundle);
}

dll int cb_qpsummodelblock_computed_step(QPSumModelBlock* self, const Matrix* modelxstep, Integer startindex_model, const Matrix* modelconstrstep, Integer startindex_constr) {
  return self->computed_step(*modelxstep, startindex_model, *modelconstrstep, startindex_constr);
}

dll int cb_qpsummodelblock_do_step(QPSumModelBlock* self, Real alpha, const Matrix* y, MinorantBundle* global_bundle, Integer startindex_bundle) {
  return self->do_step(alpha, *y, *global_bundle, startindex_bundle);
}

dll int cb_qpsummodelblock_add_localrhs(QPSumModelBlock* self, Matrix* globalrhs, Real rhsmu, Real rhscorr, Integer startindex_model, Integer startindex_constraints, int append, MinorantBundle* bundle, Integer startindex_bundel) {
  return self->add_localrhs(*globalrhs, rhsmu, rhscorr, startindex_model, startindex_constraints, (bool)append, *bundle, startindex_bundel);
}

dll int cb_qpsummodelblock_add_btinvsysb(QPSumModelBlock* self, Symmatrix* globalsys, MinorantBundle* bundle, Integer startindex_bundle) {
  return self->add_BtinvsysB(*globalsys, *bundle, startindex_bundle);
}

dll int cb_qpsummodelblock_solve_constrsys(QPSumModelBlock* self, const Symmatrix* ABchol, const Matrix* LinvABrhs, Matrix* LinvABsol, Integer startindex_model, Matrix* Crhs_and_sol, Integer startindex_constraints) {
  return self->solve_constrsys(*ABchol, *LinvABrhs, *LinvABsol, startindex_model, *Crhs_and_sol, startindex_constraints);
}

dll int cb_qpsummodelblock_add_bcschur_diagonal(QPSumModelBlock* self, Matrix* diagonal, MinorantBundle* globalbundle, Integer startindex_bundle) {
  return self->add_BCSchur_diagonal(*diagonal, *globalbundle, startindex_bundle);
}

dll int cb_qpsummodelblock_propose_bcschur_pcsubspace(QPSumModelBlock* self, Matrix* lowrank, Matrix* sigma_guess, const Matrix* Diag_inv, Real minval, Real diaginvval, MinorantBundle* globalbundle, Integer startindex_bundle) {
  return self->propose_BCSchur_pcsubspace(*lowrank, *sigma_guess, *Diag_inv, minval, diaginvval, *globalbundle, startindex_bundle);
}

dll int cb_qpsummodelblock_prepare_bcschur_jlprecond(QPSumModelBlock* self, Matrix* glob_lowrank, Matrix* subspace, int append_globtransp_times_mat_to_subspace, MinorantBundle* globalbundle, Integer startindex_bundle) {
  return self->prepare_BCSchur_JLprecond(*glob_lowrank, *subspace, (bool)append_globtransp_times_mat_to_subspace, *globalbundle, startindex_bundle);
}

dll int cb_qpsummodelblock_add_schur_rhs(QPSumModelBlock* self, Matrix* glob_rhs, Matrix* local_rhs, Real rhsmu, Real rhscorr, Integer startindex_constraints, MinorantBundle* globalbundle, Integer startindex_bundle) {
  return self->add_Schur_rhs(*glob_rhs, local_rhs, rhsmu, rhscorr, startindex_constraints, *globalbundle, startindex_bundle);
}

dll int cb_qpsummodelblock_add_schur_mult(QPSumModelBlock* self, const Matrix* in_vec, Matrix* out_vec, const Matrix* in_cvec, Matrix* out_cvec, Integer startindex_constraints, MinorantBundle* globalbundle, Integer startindex_bundle) {
  return self->add_Schur_mult(*in_vec, *out_vec, in_cvec, out_cvec, startindex_constraints, *globalbundle, startindex_bundle);
}

dll int cb_qpsummodelblock_computed_schur_step(QPSumModelBlock* self, const Matrix* xstep, const Matrix* local_step, Integer startindex_model, MinorantBundle* globalbundle, Integer startindex_bundle) {
  return self->computed_Schur_step(*xstep, *local_step, startindex_model, *globalbundle, startindex_bundle);
}

dll Integer cb_qpsummodelblock_dim_model(QPSumModelBlock* self) {
  return self->dim_model();
}

dll Integer cb_qpsummodelblock_dim_constraints(QPSumModelBlock* self) {
  return self->dim_constraints();
}

dll Real cb_qpsummodelblock_constraints_cost(QPSumModelBlock* self) {
  return self->constraints_cost();
}

dll Real cb_qpsummodelblock_primalviol_2normsqr(QPSumModelBlock* self) {
  return self->primalviol_2normsqr();
}

dll Real cb_qpsummodelblock_dualviol_2normsqr(QPSumModelBlock* self) {
  return self->dualviol_2normsqr();
}

dll int cb_qpsummodelblock_get_mu_info(const QPSumModelBlock* self, Integer* mudim, Real* tr_xz, Real* tr_xdzpdxz, Real* tr_dxdz, Real* min_xz, Real* max_xz) {
  return self->get_mu_info(*mudim, *tr_xz, *tr_xdzpdxz, *tr_dxdz, *min_xz, *max_xz);
}

dll int cb_qpsummodelblock_get_nbh_info(const QPSumModelBlock* self, Integer mudim, Real tr_xz, Real tr_xdzpdxz, Real tr_dxdz, Real nbh_ubnd, Real* alpha, Real* max_nbh, Real* nrmsqr_xz, Real* nrmsqr_xdzpdxz, Real* nrmsqr_dxdz, Real* ip_xz_xdzpdxz, Real* ip_xz_dxdz, Real* ip_dxdz_xdzpdxz) {
  return self->get_nbh_info(mudim, tr_xz, tr_xdzpdxz, tr_dxdz, nbh_ubnd, *alpha, *max_nbh, *nrmsqr_xz, *nrmsqr_xdzpdxz, *nrmsqr_dxdz, *ip_xz_xdzpdxz, *ip_xz_dxdz, *ip_dxdz_xdzpdxz);
}

dll int cb_qpsummodelblock_linesearch(const QPSumModelBlock* self, Real* alpha) {
  return self->linesearch(*alpha);
}

dll int cb_qpsummodelblock_add_localsys(QPSumModelBlock* self, Symmatrix* globalsys, Integer startindex_model, Integer startindex_constraints) {
  return self->add_localsys(*globalsys, startindex_model, startindex_constraints);
}

dll int cb_qpsummodelblock_localsys_mult(QPSumModelBlock* self, const Matrix* in_vec, Matrix* out_vec, Integer startindex_model, Integer startindex_constraints) {
  return self->localsys_mult(*in_vec, *out_vec, startindex_model, startindex_constraints);
}

