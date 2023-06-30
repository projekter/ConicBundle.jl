dll void cb_minorantusedata_destroy(MinorantUseData* self) {
  delete self;
}

dll MinorantUseData* cb_minorantusedata_new(Minorant* mp, Real sval, Integer modif_id) {
  return new MinorantUseData(mp, sval, modif_id);
}

dll MinorantUseData* cb_minorantusedata_new2(MinorantUseData* mdp, Real sval) {
  return new MinorantUseData(mdp, sval);
}

dll int cb_minorantusedata_valid(const MinorantUseData* self) {
  return self->valid();
}

dll Minorant* cb_minorantusedata_get_minorant(const MinorantUseData* self) {
  return self->get_minorant();
}

dll Integer* cb_minorantusedata_set_modification_id(MinorantUseData* self) {
  return &self->set_modification_id();
}

dll Integer cb_minorantusedata_get_modification_id(const MinorantUseData* self) {
  return self->get_modification_id();
}

dll int cb_minorantusedata_synchronize_ids(MinorantUseData* self, Integer new_modification_id, Integer new_center_id, Integer old_center_id, Integer new_cand_id, Integer old_cand_id, Integer new_prex_id) {
  return self->synchronize_ids(new_modification_id, new_center_id, old_center_id, new_cand_id, old_cand_id, new_prex_id);
}

dll int cb_minorantusedata_aggregate(const MinorantUseData* self) {
  return self->aggregate();
}

dll int cb_minorantusedata_aggregated(MinorantUseData* self, Integer n) {
  return self->aggregated(n);
}

dll Real cb_minorantusedata_offset(const MinorantUseData* self) {
  return self->offset();
}

dll Real cb_minorantusedata_coeff(const MinorantUseData* self, Integer i) {
  return self->coeff(i);
}

dll int cb_minorantusedata_scale(MinorantUseData* self, Real factor) {
  return self->scale(factor);
}

dll int cb_minorantusedata_call_primal_extender(MinorantUseData* self, PrimalExtender* prex, Integer in_prex_id) {
  return self->call_primal_extender(*prex, in_prex_id);
}

dll Real cb_minorantusedata_evaluate(const MinorantUseData* self, Integer yid, const Matrix* y, int with_constant = 1) {
  return self->evaluate(yid, *y, (bool)with_constant);
}

dll int cb_minorantusedata_one_user(const MinorantUseData* self) {
  return self->one_user();
}

