dll void cb_uqpsummodelblock_destroy(UQPSumModelBlock* self) {
  delete self;
}

dll void cb_uqpsummodelblock_clear(UQPSumModelBlock* self) {
  self->clear();
}

dll UQPSumModelBlock* cb_uqpsummodelblock_new(int incr = -1) {
  return new UQPSumModelBlock(0, incr);
}

dll int cb_uqpsummodelblock_append(UQPSumModelBlock* self, QPModelDataObject* inblock) {
  return self->append(inblock);
}

dll Integer cb_uqpsummodelblock_xdim(const UQPSumModelBlock* self) {
  return self->xdim();
}

dll Integer cb_uqpsummodelblock_ydim(const UQPSumModelBlock* self) {
  return self->ydim();
}

dll int cb_uqpsummodelblock_set_qp_xstart(UQPSumModelBlock* self, Integer x_start_index) {
  return self->set_qp_xstart(x_start_index);
}

dll int cb_uqpsummodelblock_set_qp_ystart(UQPSumModelBlock* self, Integer y_start_index) {
  return self->set_qp_ystart(y_start_index);
}

dll int cb_uqpsummodelblock_starting_x(UQPSumModelBlock* self, Matrix* qp_x) {
  return self->starting_x(*qp_x);
}

dll int cb_uqpsummodelblock_starting_y(UQPSumModelBlock* self, Matrix* qp_y, const Matrix* qp_Qx, const Matrix* qp_c) {
  return self->starting_y(*qp_y, *qp_Qx, *qp_c);
}

dll Real cb_uqpsummodelblock_get_local_primalcost(const UQPSumModelBlock* self) {
  return self->get_local_primalcost();
}

dll Real cb_uqpsummodelblock_get_local_dualcost(const UQPSumModelBlock* self) {
  return self->get_local_dualcost();
}

dll int cb_uqpsummodelblock_get_ab(const UQPSumModelBlock* self, Matrix* qp_A, Matrix* qp_b) {
  return self->get_Ab(*qp_A, *qp_b);
}

dll int cb_uqpsummodelblock_restart_x(UQPSumModelBlock* self, Matrix* qp_x, const Matrix* qp_c, const Matrix* qp_dc) {
  return self->restart_x(*qp_x, *qp_c, *qp_dc);
}

dll int cb_uqpsummodelblock_restart_y(UQPSumModelBlock* self, Matrix* qp_y, const Matrix* qp_Qx, const Matrix* qp_c, const Matrix* qp_dc) {
  return self->restart_y(*qp_y, *qp_Qx, *qp_c, *qp_dc);
}

dll int cb_uqpsummodelblock_add_xinv_kron_z(UQPSumModelBlock* self, Symmatrix* barQ) {
  return self->add_xinv_kron_z(*barQ);
}

dll int cb_uqpsummodelblock_add_local_sys(UQPSumModelBlock* self, Symmatrix* sysdy, Matrix* rhs) {
  return self->add_local_sys(*sysdy, *rhs);
}

dll int cb_uqpsummodelblock_suggest_mu(UQPSumModelBlock* self, Real* ip_xz, Integer* mu_dim, Real* sigma, const Matrix* qp_dx, const Matrix* qp_dy, const Matrix* rhs_residual) {
  return self->suggest_mu(*ip_xz, *mu_dim, *sigma, *qp_dx, *qp_dy, *rhs_residual);
}

dll int cb_uqpsummodelblock_get_corr(UQPSumModelBlock* self, Matrix* xcorr, Matrix* rhs, Real mu) {
  return self->get_corr(*xcorr, *rhs, mu);
}

dll int cb_uqpsummodelblock_line_search(UQPSumModelBlock* self, Real* alpha, const Matrix* qp_dx, const Matrix* qp_dy, const Matrix* rhs_residual) {
  return self->line_search(*alpha, *qp_dx, *qp_dy, *rhs_residual);
}

dll int cb_uqpsummodelblock_set_point(UQPSumModelBlock* self, const Matrix* qp_x, const Matrix* qp_y, Real alpha) {
  return self->set_point(*qp_x, *qp_y, alpha);
}

dll int cb_uqpsummodelblock_add_modelx_aggregate(UQPSumModelBlock* self, Real* offset, Matrix* gradient) {
  return self->add_modelx_aggregate(*offset, *gradient);
}

dll void cb_uqpsummodelblock_set_out(UQPSumModelBlock* self, int pril = 1) {
  self->set_out(&std::cout, pril);
}

dll void cb_uqpsummodelblock_set_cbout(UQPSumModelBlock* self, int incr = -1) {
  self->set_cbout(0, incr);
}

dll Matrix* cb_uqpsummodelblock_add_bs(const UQPSumModelBlock* self, Matrix* qp_vec) {
  return &self->add_Bs(*qp_vec);
}

dll Matrix* cb_uqpsummodelblock_subtract_z(const UQPSumModelBlock* self, Matrix* dual_residual, int with_step = 0) {
  return &self->subtract_z(*dual_residual, (bool)with_step);
}

