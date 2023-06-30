@doc raw"""
    cb_start!(self::CBClock)

read current time, all further time measurements will be in relation to this time
"""
cb_start!(self::CBClock) = @ccall libcb.cb_clock_start(self.data::Ptr{Cvoid})::Cvoid

@doc raw"""
    CBClock()

calls start()
"""
CBClock() = CBClock(@ccall libcb.cb_clock_new()::Ptr{Cvoid})

@doc raw"""
    cb_set_offset!(self::CBClock, offs::CBMicroseconds)

allows to specify an offset, that will furtheron be added to all time measurements
"""
cb_set_offset!(self::CBClock, offs::CBMicroseconds) = @ccall libcb.cb_clock_set_offset(self.data::Ptr{Cvoid}, offs.data::Ptr{Cvoid})::Cvoid

@doc raw"""
    cb_time(self::CBClock)

return time elapsed since last call to start() in Microseconds (possibly adding an optional offset)
"""
cb_time(self::CBClock) = CBMicroseconds(@ccall libcb.cb_clock_new_time(self.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_wall_time(self::CBClock)

returns the Microseconds passed on the wall clock sind initialization
"""
cb_wall_time(self::CBClock) = CBMicroseconds(@ccall libcb.cb_clock_new_wall_time(self.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_elapsed_time(self::CBClock)

call time() and print the result in format "hh:mm:ss", togehter with current date and time, to out
"""
cb_elapsed_time(self::CBClock) = @ccall libcb.cb_clock_elapsed_time(self.data::Ptr{Cvoid})::Cvoid

