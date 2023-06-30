dll void cb_uqpconemodelblock_destroy(UQPConeModelBlock* self) {
  delete self;
}

dll UQPConeModelBlock* cb_uqpconemodelblock_new(int cbinc = -1) {
  return new UQPConeModelBlock(0, cbinc);
}

dll void cb_uqpconemodelblock_clear(UQPConeModelBlock* self) {
  self->clear();
}

dll int cb_uqpconemodelblock_init(UQPConeModelBlock* self, const MinorantPointer* constant_minorant, const MinorantBundle* bundle, Integer nnc_dim, const Indexmatrix* soc_dim, const Indexmatrix* sdp_dim, const Matrix* box_lb, const Matrix* box_ub, Real b, int ft, QPModelOracleDataObject* oracle_data = 0, int scale_box = 1) {
  return self->init(*constant_minorant, *bundle, nnc_dim, *soc_dim, *sdp_dim, *box_lb, *box_ub, b, (FunctionTask)ft, oracle_data, (bool)scale_box);
}

dll int cb_uqpconemodelblock_adjust_trace(UQPConeModelBlock* self, Real b) {
  return self->adjust_trace(b);
}

dll Real cb_uqpconemodelblock_evaluate_trace(const UQPConeModelBlock* self) {
  return self->evaluate_trace();
}

dll Real cb_uqpconemodelblock_get_trace(UQPConeModelBlock* self) {
  return self->get_trace();
}

dll int cb_uqpconemodelblock_get_nncx(UQPConeModelBlock* self, Matrix* nncx, Matrix* nncx_activity = 0, int cautious = 0) {
  return self->get_nncx(*nncx, nncx_activity, (bool)cautious);
}

dll int cb_uqpconemodelblock_get_socx(UQPConeModelBlock* self, Integer i, Matrix* socx, Real* socx_activity, int cautious = 0) {
  return self->get_socx(i, *socx, socx_activity, (bool)cautious);
}

dll int cb_uqpconemodelblock_get_pscx(UQPConeModelBlock* self, Integer i, Matrix* pscx_eigs, Matrix* pscx_vecs, Real* pscx_growthrate, Matrix* pscx_primalgrowth, Matrix* pscx_dualgrowth) {
  return self->get_pscx(i, *pscx_eigs, *pscx_vecs, *pscx_growthrate, *pscx_primalgrowth, *pscx_dualgrowth);
}

dll int cb_uqpconemodelblock_get_boxx(UQPConeModelBlock* self, Matrix* param0, Matrix* param1 = 0, int cautious = 0) {
  return self->get_boxx(*param0, param1, (bool)cautious);
}

dll int cb_uqpconemodelblock_add_modelx_aggregate(UQPConeModelBlock* self, Real* offset, Matrix* gradient) {
  return self->add_modelx_aggregate(*offset, *gradient);
}

dll Real cb_uqpconemodelblock_tracedual(const UQPConeModelBlock* self, Real* prec = 0) {
  return self->tracedual(prec);
}

dll int cb_uqpconemodelblock_get_nncz(UQPConeModelBlock* self, Matrix* vecz) {
  return self->get_nncz(*vecz);
}

dll int cb_uqpconemodelblock_get_socz(UQPConeModelBlock* self, Integer i, Matrix* vecz) {
  return self->get_socz(i, *vecz);
}

dll int cb_uqpconemodelblock_get_x(UQPConeModelBlock* self, Integer i, Symmatrix* X) {
  return self->get_X(i, *X);
}

dll int cb_uqpconemodelblock_get_z(UQPConeModelBlock* self, Integer i, Symmatrix* Z) {
  return self->get_Z(i, *Z);
}

dll Real cb_uqpconemodelblock_get_y(UQPConeModelBlock* self) {
  return self->get_y();
}

dll Real cb_uqpconemodelblock_get_s(UQPConeModelBlock* self) {
  return self->get_s();
}

dll Real cb_uqpconemodelblock_get_mu(UQPConeModelBlock* self) {
  return self->get_mu();
}

dll int cb_uqpconemodelblock_get_old_nncx(UQPConeModelBlock* self, Matrix* vecx) {
  return self->get_old_nncx(*vecx);
}

dll int cb_uqpconemodelblock_get_old_nncz(UQPConeModelBlock* self, Matrix* vecz) {
  return self->get_old_nncz(*vecz);
}

dll int cb_uqpconemodelblock_get_old_socx(UQPConeModelBlock* self, Integer i, Matrix* vecx) {
  return self->get_old_socx(i, *vecx);
}

dll int cb_uqpconemodelblock_get_old_socz(UQPConeModelBlock* self, Integer i, Matrix* vecz) {
  return self->get_old_socz(i, *vecz);
}

dll int cb_uqpconemodelblock_get_old_x(UQPConeModelBlock* self, Integer i, Symmatrix* X) {
  return self->get_old_X(i, *X);
}

dll int cb_uqpconemodelblock_get_old_z(UQPConeModelBlock* self, Integer i, Symmatrix* Z) {
  return self->get_old_Z(i, *Z);
}

dll Real cb_uqpconemodelblock_get_old_y(UQPConeModelBlock* self) {
  return self->get_old_y();
}

dll Real cb_uqpconemodelblock_get_old_s(UQPConeModelBlock* self) {
  return self->get_old_s();
}

dll Real cb_uqpconemodelblock_get_old_mu(UQPConeModelBlock* self) {
  return self->get_old_mu();
}

dll Real cb_uqpconemodelblock_get_last_alpha(UQPConeModelBlock* self) {
  return self->get_last_alpha();
}

dll int cb_uqpconemodelblock_compute_local_directions(UQPConeModelBlock* self, const Matrix* qp_dx, const Matrix* qp_dy, const Matrix* rhs_resid, Matrix* dz, Matrix* duz, Real* box_ds, Real* ds) {
  return self->compute_local_directions(*qp_dx, *qp_dy, *rhs_resid, *dz, *duz, *box_ds, *ds);
}

dll int cb_uqpconemodelblock_inner_line_search(UQPConeModelBlock* self, Real* alpha, const Matrix* qp_dx, const Matrix* qp_dy, const Matrix* dz, const Matrix* duz, Real box_ds, Real ds) {
  return self->inner_line_search(*alpha, *qp_dx, *qp_dy, *dz, *duz, box_ds, ds);
}

dll Integer cb_uqpconemodelblock_xdim(const UQPConeModelBlock* self) {
  return self->xdim();
}

dll Integer cb_uqpconemodelblock_ydim(const UQPConeModelBlock* self) {
  return self->ydim();
}

dll int cb_uqpconemodelblock_set_qp_xstart(UQPConeModelBlock* self, Integer x_start_index) {
  return self->set_qp_xstart(x_start_index);
}

dll int cb_uqpconemodelblock_set_qp_ystart(UQPConeModelBlock* self, Integer y_start_index) {
  return self->set_qp_ystart(y_start_index);
}

dll int cb_uqpconemodelblock_starting_x(UQPConeModelBlock* self, Matrix* qp_x) {
  return self->starting_x(*qp_x);
}

dll int cb_uqpconemodelblock_starting_y(UQPConeModelBlock* self, Matrix* qp_y, const Matrix* qp_Qx, const Matrix* qp_c) {
  return self->starting_y(*qp_y, *qp_Qx, *qp_c);
}

dll int cb_uqpconemodelblock_get_ab(const UQPConeModelBlock* self, Matrix* qp_A, Matrix* qp_b) {
  return self->get_Ab(*qp_A, *qp_b);
}

dll int cb_uqpconemodelblock_restart_x(UQPConeModelBlock* self, Matrix* qp_x, const Matrix* qp_c, const Matrix* qp_dc) {
  return self->restart_x(*qp_x, *qp_c, *qp_dc);
}

dll int cb_uqpconemodelblock_restart_y(UQPConeModelBlock* self, Matrix* qp_y, const Matrix* qp_Qx, const Matrix* qp_c, const Matrix* qp_dc) {
  return self->restart_y(*qp_y, *qp_Qx, *qp_c, *qp_dc);
}

dll int cb_uqpconemodelblock_add_xinv_kron_z(UQPConeModelBlock* self, Symmatrix* barQ) {
  return self->add_xinv_kron_z(*barQ);
}

dll int cb_uqpconemodelblock_add_local_sys(UQPConeModelBlock* self, Symmatrix* sysdy, Matrix* rhs) {
  return self->add_local_sys(*sysdy, *rhs);
}

dll Real cb_uqpconemodelblock_get_local_primalcost(const UQPConeModelBlock* self) {
  return self->get_local_primalcost();
}

dll Real cb_uqpconemodelblock_get_local_dualcost(const UQPConeModelBlock* self) {
  return self->get_local_dualcost();
}

dll int cb_uqpconemodelblock_suggest_mu(UQPConeModelBlock* self, Real* ip_xz, Integer* mu_dim, Real* sigma, const Matrix* qp_dx, const Matrix* qp_dy, const Matrix* rhs_residual) {
  return self->suggest_mu(*ip_xz, *mu_dim, *sigma, *qp_dx, *qp_dy, *rhs_residual);
}

dll int cb_uqpconemodelblock_get_corr(UQPConeModelBlock* self, Matrix* xcorr, Matrix* rhs, Real mu) {
  return self->get_corr(*xcorr, *rhs, mu);
}

dll int cb_uqpconemodelblock_line_search(UQPConeModelBlock* self, Real* alpha, const Matrix* qp_dx, const Matrix* qp_dy, const Matrix* rhs_residual) {
  return self->line_search(*alpha, *qp_dx, *qp_dy, *rhs_residual);
}

dll int cb_uqpconemodelblock_set_point(UQPConeModelBlock* self, const Matrix* qp_x, const Matrix* qp_y, Real alpha) {
  return self->set_point(*qp_x, *qp_y, alpha);
}

dll Matrix* cb_uqpconemodelblock_add_bs(const UQPConeModelBlock* self, Matrix* qp_vec) {
  return &self->add_Bs(*qp_vec);
}

dll Matrix* cb_uqpconemodelblock_subtract_z(const UQPConeModelBlock* self, Matrix* dual_residual, int with_step = 0) {
  return &self->subtract_z(*dual_residual, (bool)with_step);
}

