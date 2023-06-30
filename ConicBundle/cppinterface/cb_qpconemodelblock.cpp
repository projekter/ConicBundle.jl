dll void cb_qpconemodelblock_destroy(QPConeModelBlock* self) {
  delete self;
}

dll Matrix* cb_qpconemodelblock_b_times(QPConeModelBlock* self, const Matrix* A, Matrix* C, Real alpha, Real beta, int Btrans, int Atrans, Integer startindex_model, MinorantBundle* globalbundle, Integer startindex_bundle) {
  return &self->B_times(*A, *C, alpha, beta, Btrans, Atrans, startindex_model, *globalbundle, startindex_bundle);
}

dll Matrix* cb_qpconemodelblock_times_b(QPConeModelBlock* self, const Matrix* A, Matrix* C, Real alpha, Real beta, int Atrans, int Btrans, Integer startindex_model, MinorantBundle* globalbundle, Integer startindex_bundle) {
  return &self->times_B(*A, *C, alpha, beta, Atrans, Btrans, startindex_model, *globalbundle, startindex_bundle);
}

dll Symmatrix* cb_qpconemodelblock_add_bdbt(QPConeModelBlock* self, const Matrix* diagvec, Symmatrix* bigS, int minus, Integer startindex, Matrix* Bt, Integer startindex_model, MinorantBundle* globalbundle, Integer startindex_bundle) {
  return &self->add_BDBt(*diagvec, *bigS, (bool)minus, startindex, *Bt, startindex_model, *globalbundle, startindex_bundle);
}

dll Matrix* cb_qpconemodelblock_get_bt(QPConeModelBlock* self, Matrix* Bt, Integer startindex_model, MinorantBundle* global_bundle, Integer startindex_bundle) {
  return &self->get_Bt(*Bt, startindex_model, *global_bundle, startindex_bundle);
}

dll int cb_qpconemodelblock_get_modelx(QPConeModelBlock* self, Matrix* modelx, Integer startindex_model) {
  return self->get_modelx(*modelx, startindex_model);
}

dll int cb_qpconemodelblock_get_modeldx(QPConeModelBlock* self, Matrix* modeldx, Integer startindex_model) {
  return self->get_modeldx(*modeldx, startindex_model);
}

dll int cb_qpconemodelblock_get_modeldcstr(QPConeModelBlock* self, Matrix* modeldcstr, Integer startindex_constraints) {
  return self->get_modeldcstr(*modeldcstr, startindex_constraints);
}

dll int cb_qpconemodelblock_add_modelx_aggregate(QPConeModelBlock* self, Real* val, Matrix* vec, MinorantBundle* global_bundle, Integer startindex_bundle) {
  return self->add_modelx_aggregate(*val, *vec, *global_bundle, startindex_bundle);
}

dll int cb_qpconemodelblock_get_sysviol_model(QPConeModelBlock* self, Matrix* modelvec, Integer startindex_model, const Matrix* y_plus_dy, MinorantBundle* global_bundle, Integer startindex_bundle) {
  return self->get_sysviol_model(*modelvec, startindex_model, *y_plus_dy, *global_bundle, startindex_bundle);
}

dll int cb_qpconemodelblock_get_sysviol_constraints(QPConeModelBlock* self, Matrix* constrvec, Integer startindex_constr) {
  return self->get_sysviol_constraints(*constrvec, startindex_constr);
}

dll void cb_qpconemodelblock_display_model_values(QPConeModelBlock* self, const Matrix* y, MinorantBundle* global_bundle, Integer startindex_bundle) {
  self->display_model_values(*y, *global_bundle, startindex_bundle, std::cout);
}

dll int cb_qpconemodelblock_reset_starting_point(QPConeModelBlock* self, const Matrix* y, Real mu, MinorantBundle* global_bundle, Integer startindex_bundle) {
  return self->reset_starting_point(*y, mu, *global_bundle, startindex_bundle);
}

dll int cb_qpconemodelblock_compute_step(QPConeModelBlock* self, const Matrix* ystep, MinorantBundle* global_bundle, Integer startindex_bundle) {
  return self->compute_step(*ystep, *global_bundle, startindex_bundle);
}

dll int cb_qpconemodelblock_computed_step(QPConeModelBlock* self, const Matrix* modelxstep, Integer startindex_model, const Matrix* modelconstrstep, Integer startindex_constr) {
  return self->computed_step(*modelxstep, startindex_model, *modelconstrstep, startindex_constr);
}

dll int cb_qpconemodelblock_do_step(QPConeModelBlock* self, Real alpha, const Matrix* y, MinorantBundle* global_bundle, Integer startindex_bundle) {
  return self->do_step(alpha, *y, *global_bundle, startindex_bundle);
}

dll int cb_qpconemodelblock_add_localrhs(QPConeModelBlock* self, Matrix* globalrhs, Real rhsmu, Real rhscorr, Integer startindex_model, Integer startindex_constraints, int append, MinorantBundle* bundle, Integer startindex_bundel) {
  return self->add_localrhs(*globalrhs, rhsmu, rhscorr, startindex_model, startindex_constraints, (bool)append, *bundle, startindex_bundel);
}

dll int cb_qpconemodelblock_add_btinvsysb(QPConeModelBlock* self, Symmatrix* globalsys, MinorantBundle* bundle, Integer startindex_bundle) {
  return self->add_BtinvsysB(*globalsys, *bundle, startindex_bundle);
}

dll int cb_qpconemodelblock_solve_constrsys(QPConeModelBlock* self, const Symmatrix* ABchol, const Matrix* LinvABrhs, Matrix* LinvABsol, Integer startindex_model, Matrix* Crhs_and_sol, Integer startindex_constraints) {
  return self->solve_constrsys(*ABchol, *LinvABrhs, *LinvABsol, startindex_model, *Crhs_and_sol, startindex_constraints);
}

dll void cb_qpconemodelblock_clear(QPConeModelBlock* self) {
  self->clear();
}

dll QPConeModelBlock* cb_qpconemodelblock_new(int cbinc = -1) {
  return new QPConeModelBlock(0, cbinc);
}

dll QPModelBlockObject* cb_qpconemodelblock_clone(QPConeModelBlock* self) {
  return self->clone();
}

dll void cb_qpconemodelblock_recursive_delete_and_clear(QPConeModelBlock* self) {
  self->recursive_delete_and_clear();
}

dll int cb_qpconemodelblock_recursive_copy_data_of(QPConeModelBlock* self, QPModelBlockObject* param0) {
  return self->recursive_copy_data_of(param0);
}

dll int cb_qpconemodelblock_init(QPConeModelBlock* self, const MinorantPointer* constant_minorant, const MinorantBundle* bundle, Integer nnc_dim, const Indexmatrix* soc_dim, const Indexmatrix* psc_dim, const Matrix* box_lb, const Matrix* box_ub, Real b, int ft, QPModelOracleDataObject* oracle_data = 0, int scale_box = 1) {
  return self->init(*constant_minorant, *bundle, nnc_dim, *soc_dim, *psc_dim, *box_lb, *box_ub, b, (FunctionTask)ft, oracle_data, (bool)scale_box);
}

dll int cb_qpconemodelblock_adjust_trace(QPConeModelBlock* self, Real b) {
  return self->adjust_trace(b);
}

dll Real cb_qpconemodelblock_evaluate_trace(const QPConeModelBlock* self) {
  return self->evaluate_trace();
}

dll Real cb_qpconemodelblock_get_trace(QPConeModelBlock* self) {
  return self->get_trace();
}

dll int cb_qpconemodelblock_get_nncx(QPConeModelBlock* self, Matrix* nncx, Matrix* nncx_activity = 0, int cautious = 0) {
  return self->get_nncx(*nncx, nncx_activity, (bool)cautious);
}

dll int cb_qpconemodelblock_get_socx(QPConeModelBlock* self, Integer i, Matrix* socx, Real* socx_activity, int cautious = 0) {
  return self->get_socx(i, *socx, socx_activity, (bool)cautious);
}

dll int cb_qpconemodelblock_get_pscx(QPConeModelBlock* self, Integer i, Matrix* pscx_eigs, Matrix* pscx_vecs, Real* pscx_growthrate, Matrix* pscx_primalgrowth, Matrix* pscx_dualgrowth) {
  return self->get_pscx(i, *pscx_eigs, *pscx_vecs, *pscx_growthrate, *pscx_primalgrowth, *pscx_dualgrowth);
}

dll int cb_qpconemodelblock_get_boxx(QPConeModelBlock* self, Matrix* boxx, Matrix* linx_activity = 0, int cautious = 0) {
  return self->get_boxx(*boxx, linx_activity, (bool)cautious);
}

dll Real cb_qpconemodelblock_tracedual(const QPConeModelBlock* self, Real* prec = 0) {
  return self->tracedual(prec);
}

dll Integer cb_qpconemodelblock_dim_model(QPConeModelBlock* self) {
  return self->dim_model();
}

dll Integer cb_qpconemodelblock_dim_constraints(QPConeModelBlock* self) {
  return self->dim_constraints();
}

dll Real cb_qpconemodelblock_constraints_cost(QPConeModelBlock* self) {
  return self->constraints_cost();
}

dll Real cb_qpconemodelblock_primalviol_2normsqr(QPConeModelBlock* self) {
  return self->primalviol_2normsqr();
}

dll Real cb_qpconemodelblock_dualviol_2normsqr(QPConeModelBlock* self) {
  return self->dualviol_2normsqr();
}

dll int cb_qpconemodelblock_get_mu_info(const QPConeModelBlock* self, Integer* mudim, Real* tr_xz, Real* tr_xdzpdxz, Real* tr_dxdz, Real* min_xz, Real* max_xz) {
  return self->get_mu_info(*mudim, *tr_xz, *tr_xdzpdxz, *tr_dxdz, *min_xz, *max_xz);
}

dll int cb_qpconemodelblock_get_nbh_info(const QPConeModelBlock* self, Integer mudim, Real tr_xz, Real tr_xdzpdxz, Real tr_dxdz, Real nbh_ubnd, Real* alpha, Real* max_nbh, Real* nrmsqr_xz, Real* nrmsqr_xdzpdxz, Real* nrmsqr_dxdz, Real* ip_xz_xdzpdxz, Real* ip_xz_dxdz, Real* ip_dxdz_xdzpdxz) {
  return self->get_nbh_info(mudim, tr_xz, tr_xdzpdxz, tr_dxdz, nbh_ubnd, *alpha, *max_nbh, *nrmsqr_xz, *nrmsqr_xdzpdxz, *nrmsqr_dxdz, *ip_xz_xdzpdxz, *ip_xz_dxdz, *ip_dxdz_xdzpdxz);
}

dll int cb_qpconemodelblock_linesearch(const QPConeModelBlock* self, Real* alpha) {
  return self->linesearch(*alpha);
}

dll int cb_qpconemodelblock_add_localsys(QPConeModelBlock* self, Symmatrix* globalsys, Integer startindex_model, Integer startindex_constraints) {
  return self->add_localsys(*globalsys, startindex_model, startindex_constraints);
}

dll int cb_qpconemodelblock_localsys_mult(QPConeModelBlock* self, const Matrix* in_vec, Matrix* out_vec, Integer startindex_model, Integer startindex_constraints) {
  return self->localsys_mult(*in_vec, *out_vec, startindex_model, startindex_constraints);
}

dll int cb_qpconemodelblock_add_bcschur_diagonal(QPConeModelBlock* self, Matrix* diagonal, MinorantBundle* globalbundle, Integer startindex_bundle) {
  return self->add_BCSchur_diagonal(*diagonal, *globalbundle, startindex_bundle);
}

dll int cb_qpconemodelblock_propose_bcschur_pcsubspace(QPConeModelBlock* self, Matrix* lowrank, Matrix* sigma_guess, const Matrix* Diag_inv, Real minval, Real diaginvval, MinorantBundle* globalbundle, Integer startindex_bundle) {
  return self->propose_BCSchur_pcsubspace(*lowrank, *sigma_guess, *Diag_inv, minval, diaginvval, *globalbundle, startindex_bundle);
}

dll int cb_qpconemodelblock_prepare_bcschur_jlprecond(QPConeModelBlock* self, Matrix* glob_lowrank, Matrix* subspace, int append_globtransp_times_mat_to_subspace, MinorantBundle* globalbundle, Integer startindex_bundle) {
  return self->prepare_BCSchur_JLprecond(*glob_lowrank, *subspace, (bool)append_globtransp_times_mat_to_subspace, *globalbundle, startindex_bundle);
}

dll int cb_qpconemodelblock_add_schur_rhs(QPConeModelBlock* self, Matrix* glob_rhs, Matrix* local_rhs, Real rhsmu, Real rhscorr, Integer startindex_constraints, MinorantBundle* globalbundle, Integer startindex_bundle) {
  return self->add_Schur_rhs(*glob_rhs, local_rhs, rhsmu, rhscorr, startindex_constraints, *globalbundle, startindex_bundle);
}

dll int cb_qpconemodelblock_add_schur_mult(QPConeModelBlock* self, const Matrix* in_vec, Matrix* out_vec, const Matrix* in_cvec, Matrix* out_cvec, Integer startindex_constraints, MinorantBundle* globalbundle, Integer startindex_bundle) {
  return self->add_Schur_mult(*in_vec, *out_vec, in_cvec, out_cvec, startindex_constraints, *globalbundle, startindex_bundle);
}

dll int cb_qpconemodelblock_computed_schur_step(QPConeModelBlock* self, const Matrix* xstep, const Matrix* local_step, Integer startindex_model, MinorantBundle* globalbundle, Integer startindex_bundle) {
  return self->computed_Schur_step(*xstep, *local_step, startindex_model, *globalbundle, startindex_bundle);
}

dll void cb_qpconemodelblock_set_cbout(QPConeModelBlock* self, int incr = -1) {
  self->set_cbout(0, incr);
}

