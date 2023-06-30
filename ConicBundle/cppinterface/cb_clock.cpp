dll void cb_clock_destroy(Clock* self) {
  delete self;
}

dll void cb_clock_start(Clock* self) {
  self->start();
}

dll Clock* cb_clock_new() {
  return new Clock();
}

dll void cb_clock_set_offset(Clock* self, Microseconds* offs) {
  self->set_offset(*offs);
}

dll Microseconds* cb_clock_new_time(const Clock* self) {
  return new Microseconds(self->time());
}

dll Microseconds* cb_clock_new_wall_time(const Clock* self) {
  return new Microseconds(self->wall_time());
}

dll void cb_clock_elapsed_time(const Clock* self) {
  self->elapsed_time(std::cout);
}

