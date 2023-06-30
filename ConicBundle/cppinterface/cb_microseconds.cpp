dll void cb_microseconds_destroy(Microseconds* self) {
  delete self;
}

dll Microseconds* cb_microseconds_new() {
  return new Microseconds();
}

dll Microseconds* cb_microseconds_new2(int infty) {
  return new Microseconds((bool)infty);
}

dll Microseconds* cb_microseconds_new3(const Microseconds* m) {
  return new Microseconds(*m);
}

dll Microseconds* cb_microseconds_new4(long secs, long msecs = 0) {
  return new Microseconds(secs, msecs);
}

dll Microseconds* cb_microseconds_new5(int secs, int msecs = 0) {
  return new Microseconds(secs, msecs);
}

dll Microseconds* cb_microseconds_new6(long hours, long minutes, long secs, long micros) {
  return new Microseconds(hours, minutes, secs, micros);
}

dll Microseconds* cb_microseconds_new7(int hours, int minutes, int secs, int micros) {
  return new Microseconds(hours, minutes, secs, micros);
}

dll Microseconds* cb_microseconds_assign(Microseconds* self, const Microseconds* m) {
  return &(*self = *m);
}

dll Microseconds* cb_microseconds_plus(Microseconds* self, const Microseconds* m) {
  return &(*self += *m);
}

dll Microseconds* cb_microseconds_minus(Microseconds* self, const Microseconds* m) {
  return &(*self -= *m);
}

dll bool cb_microseconds_new_less(const Microseconds* self, const Microseconds* m) {
  return self->operator<(*m);
}

dll bool cb_microseconds_new_greater(const Microseconds* self, const Microseconds* m) {
  return self->operator>(*m);
}

dll bool cb_microseconds_new_lessequal(const Microseconds* self, const Microseconds* m) {
  return (*self <= *m);
}

dll bool cb_microseconds_new_greaterequal(const Microseconds* self, const Microseconds* m) {
  return (*self >= *m);
}

dll bool cb_microseconds_new_equal(const Microseconds* self, const Microseconds* m) {
  return (*self == *m);
}

dll void cb_microseconds_set_infinity(Microseconds* self, int infty) {
  self->set_infinity((bool)infty);
}

dll int cb_microseconds_get_infinity(const Microseconds* self) {
  return self->get_infinity();
}

dll void cb_microseconds_hhmmss(const Microseconds* self, long* hours, long* minutes, long* secs) {
  self->hhmmss(*hours, *minutes, *secs);
}

dll void cb_microseconds_hhmmssdd(const Microseconds* self, long* hours, long* minutes, long* secs, long* hund) {
  self->hhmmssdd(*hours, *minutes, *secs, *hund);
}

dll long cb_microseconds_roundsecs(const Microseconds* self) {
  return self->roundsecs();
}

dll long cb_microseconds_roundhundredths(const Microseconds* self) {
  return self->roundhundredths();
}

dll void cb_microseconds_print_time(const Microseconds* m, int secondsonly = 0) {
  print_time(std::cout, *m, secondsonly);
}

