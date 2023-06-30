dll void cb_bundlesolver_destroy(BundleSolver* self) {
  delete self;
}

dll BundleSolver* cb_bundlesolver_new2(Integer dim, BundleModel* bp = 0, int incr = -1) {
  return new BundleSolver(dim, bp, 0, incr);
}

dll BundleSolver* cb_bundlesolver_new3(Groundset* gs, BundleModel* bp = 0, int incr = -1) {
  return new BundleSolver(gs, bp, 0, incr);
}

dll void cb_bundlesolver_set_defaults(BundleSolver* self) {
  self->set_defaults();
}

dll void cb_bundlesolver_clear(BundleSolver* self) {
  self->clear();
}

dll void cb_bundlesolver_clear_fails(BundleSolver* self) {
  self->clear_fails();
}

dll int cb_bundlesolver_initialize(BundleSolver* self, Integer dim, BundleModel* bp = 0) {
  return self->initialize(dim, bp);
}

dll int cb_bundlesolver_initialize2(BundleSolver* self, Groundset* gs, BundleModel* bp = 0) {
  return self->initialize(gs, bp);
}

dll void cb_bundlesolver_set_out(BundleSolver* self, int pril = 1) {
  self->set_out(&std::cout, pril);
}

dll void cb_bundlesolver_set_cbout(BundleSolver* self, int incr) {
  self->set_cbout(0, incr);
}

dll int cb_bundlesolver_set_model(BundleSolver* self, BundleModel* bp) {
  return self->set_model(bp);
}

dll int cb_bundlesolver_set_new_center(BundleSolver* self, const Matrix* yp = 0) {
  return self->set_new_center(yp);
}

dll int cb_bundlesolver_solve(BundleSolver* self, int maxsteps = 0, int stop_at_descent_steps = 0) {
  return self->solve(maxsteps, (bool)stop_at_descent_steps);
}

dll void cb_bundlesolver_print_line_summary(const BundleSolver* self) {
  self->print_line_summary(std::cout);
}

dll int cb_bundlesolver_apply_modification2(BundleSolver* self, const GroundsetModification* gsmdf) {
  return self->apply_modification(*gsmdf);
}

dll void cb_bundlesolver_set_terminator(BundleSolver* self, BundleTerminator* bt) {
  self->set_terminator(bt);
}

dll void cb_bundlesolver_set_bundleweight(BundleSolver* self, BundleWeight* bw) {
  self->set_bundleweight(bw);
}

dll void cb_bundlesolver_set_modeleps(BundleSolver* self, Real in_eps) {
  self->set_modeleps(in_eps);
}

dll void cb_bundlesolver_set_ml(BundleSolver* self, Real in_mL) {
  self->set_mL(in_mL);
}

dll void cb_bundlesolver_set_mn(BundleSolver* self, Real in_mN) {
  self->set_mN(in_mN);
}

dll void cb_bundlesolver_set_use_linval(BundleSolver* self, int ul) {
  self->set_use_linval((bool)ul);
}

dll void cb_bundlesolver_set_do_yfixing(BundleSolver* self, int dofix) {
  self->set_do_yfixing((bool)dofix);
}

dll int cb_bundlesolver_set_variable_metric(BundleSolver* self, int ds) {
  return self->set_variable_metric(ds);
}

dll int cb_bundlesolver_set_prox_diagonal(BundleSolver* self, const Matrix* insc) {
  return self->set_prox_diagonal(*insc);
}

dll int cb_bundlesolver_set_prox(BundleSolver* self, BundleProxObject* bsp) {
  return self->set_prox(bsp);
}

dll void cb_bundlesolver_set_clock(BundleSolver* self, const CH_Tools::Clock* myclock) {
  self->set_clock(*myclock);
}

dll void cb_bundlesolver_set_max_updates(BundleSolver* self, Integer mu) {
  self->set_max_updates(mu);
}

dll int cb_bundlesolver_get_terminate(const BundleSolver* self) {
  return self->get_terminate();
}

dll Real cb_bundlesolver_get_center_objval(const BundleSolver* self) {
  return self->get_center_objval();
}

dll Real cb_bundlesolver_get_center_ub(const BundleSolver* self) {
  return self->get_center_ub();
}

dll Real cb_bundlesolver_get_center_gs_val(const BundleSolver* self) {
  return self->get_center_gs_val();
}

dll const Matrix* cb_bundlesolver_get_center_y(const BundleSolver* self) {
  return &self->get_center_y();
}

dll Real cb_bundlesolver_get_cand_objval(const BundleSolver* self) {
  return self->get_cand_objval();
}

dll Real cb_bundlesolver_get_cand_ub(const BundleSolver* self) {
  return self->get_cand_ub();
}

dll Real cb_bundlesolver_get_cand_gs_val(const BundleSolver* self) {
  return self->get_cand_gs_val();
}

dll const Matrix* cb_bundlesolver_get_cand_y(const BundleSolver* self) {
  return &self->get_cand_y();
}

dll Real cb_bundlesolver_get_aggr_dnormsqr(const BundleSolver* self) {
  return self->get_aggr_dnormsqr();
}

dll Real cb_bundlesolver_get_aggregate_offset(const BundleSolver* self) {
  return self->get_aggregate_offset();
}

dll void cb_bundlesolver_get_aggregate(const BundleSolver* self, Matrix* aggregate) {
  self->get_aggregate(*aggregate);
}

dll const MinorantPointer* cb_bundlesolver_get_gs_aggregate(const BundleSolver* self) {
  return &self->get_gs_aggregate();
}

dll const MinorantPointer* cb_bundlesolver_get_model_aggregate(const BundleSolver* self) {
  return &self->get_model_aggregate();
}

dll Real cb_bundlesolver_get_modelval(const BundleSolver* self) {
  return self->get_modelval();
}

dll Real cb_bundlesolver_get_weight(const BundleSolver* self) {
  return self->get_weight();
}

dll Real cb_bundlesolver_get_term_corr(const BundleSolver* self) {
  return self->get_term_corr();
}

dll int cb_bundlesolver_get_descent_step(const BundleSolver* self) {
  return self->get_descent_step();
}

dll int cb_bundlesolver_get_null_step(const BundleSolver* self) {
  return self->get_null_step();
}

dll const Indexmatrix* cb_bundlesolver_get_yfixed(const BundleSolver* self) {
  return self->get_yfixed();
}

dll Real cb_bundlesolver_get_modeleps(const BundleSolver* self) {
  return self->get_modeleps();
}

dll int cb_bundlesolver_get_use_linval(const BundleSolver* self) {
  return self->get_use_linval();
}

dll int cb_bundlesolver_get_do_variable_metric(const BundleSolver* self) {
  return self->get_do_variable_metric();
}

dll int cb_bundlesolver_get_use_variable_metric(const BundleSolver* self) {
  return self->get_use_variable_metric();
}

dll Integer cb_bundlesolver_get_cntobjeval(const BundleSolver* self) {
  return self->get_cntobjeval();
}

dll Integer cb_bundlesolver_get_descent_steps(const BundleSolver* self) {
  return self->get_descent_steps();
}

dll Integer cb_bundlesolver_get_innerit(const BundleSolver* self) {
  return self->get_innerit();
}

dll Integer cb_bundlesolver_get_suminnerit(const BundleSolver* self) {
  return self->get_suminnerit();
}

dll Integer cb_bundlesolver_get_sumupdatecnt(const BundleSolver* self) {
  return self->get_sumupdatecnt();
}

dll Integer cb_bundlesolver_get_recomp(const BundleSolver* self) {
  return self->get_recomp();
}

dll Integer cb_bundlesolver_get_sumrecomp(const BundleSolver* self) {
  return self->get_sumrecomp();
}

dll Integer cb_bundlesolver_get_qpfails(const BundleSolver* self) {
  return self->get_qpfails();
}

dll Integer cb_bundlesolver_get_sumqpfails(const BundleSolver* self) {
  return self->get_sumqpfails();
}

dll Integer cb_bundlesolver_get_modelfails(const BundleSolver* self) {
  return self->get_modelfails();
}

dll Integer cb_bundlesolver_get_summodelfails(const BundleSolver* self) {
  return self->get_summodelfails();
}

dll Integer cb_bundlesolver_get_augvalfails(const BundleSolver* self) {
  return self->get_augvalfails();
}

dll Integer cb_bundlesolver_get_sumaugvalfails(const BundleSolver* self) {
  return self->get_sumaugvalfails();
}

dll Integer cb_bundlesolver_get_oraclefails(const BundleSolver* self) {
  return self->get_oraclefails();
}

dll Integer cb_bundlesolver_get_sumoraclefails(const BundleSolver* self) {
  return self->get_sumoraclefails();
}

dll Integer cb_bundlesolver_get_shallowcut(const BundleSolver* self) {
  return self->get_shallowcut();
}

dll const Groundset* cb_bundlesolver_get_groundset(const BundleSolver* self) {
  return self->get_groundset();
}

dll const BundleModel* cb_bundlesolver_get_model(const BundleSolver* self) {
  return self->get_model();
}

dll BundleProxObject* cb_bundlesolver_get_prox(const BundleSolver* self) {
  return self->get_prox();
}

dll BundleTerminator* cb_bundlesolver_get_terminator(const BundleSolver* self) {
  return self->get_terminator();
}

dll BundleWeight* cb_bundlesolver_get_bundleweight(const BundleSolver* self) {
  return self->get_bundleweight();
}

dll CH_Tools::Microseconds* cb_bundlesolver_new_get_qpcoeff_time(const BundleSolver* self) {
  return new CH_Tools::Microseconds(self->get_QPcoeff_time());
}

dll CH_Tools::Microseconds* cb_bundlesolver_new_get_qpsolve_time(const BundleSolver* self) {
  return new CH_Tools::Microseconds(self->get_QPsolve_time());
}

dll CH_Tools::Microseconds* cb_bundlesolver_new_get_make_aggr_time(const BundleSolver* self) {
  return new CH_Tools::Microseconds(self->get_make_aggr_time());
}

dll CH_Tools::Microseconds* cb_bundlesolver_new_get_evalaugmodel_time(const BundleSolver* self) {
  return new CH_Tools::Microseconds(self->get_evalaugmodel_time());
}

dll void cb_bundlesolver_print_statistics(const BundleSolver* self) {
  self->print_statistics(std::cout);
}

dll int cb_bundlesolver_qp_mfile_data(const BundleSolver* self, const Matrix* center_y, const BundleProxObject* Hp, const MinorantPointer* gs_subg, const Symmatrix* Q, const Matrix* c, Real offset, const Indexmatrix* yfixed) {
  return self->qp_mfile_data(*center_y, Hp, *gs_subg, *Q, *c, offset, *yfixed);
}

