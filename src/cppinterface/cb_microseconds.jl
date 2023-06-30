@doc raw"""
    CBMicroseconds()

default constructor, value 0
"""
CBMicroseconds() = CBMicroseconds(@ccall libcb.cb_microseconds_new()::Ptr{Cvoid})

@doc raw"""
    CBMicroseconds(infty::Bool)

constructor for setting value to 0 or infinity
"""
CBMicroseconds(infty::Bool) = CBMicroseconds(@ccall libcb.cb_microseconds_new2(infty::Cint)::Ptr{Cvoid})

@doc raw"""
    CBMicroseconds(m::CBMicroseconds)

copy constructor
"""
CBMicroseconds(m::CBMicroseconds) = CBMicroseconds(@ccall libcb.cb_microseconds_new3(m.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    CBMicroseconds(secs::Integer, msecs::Integer = 0)

specify directly seconds and microseconds (in [0,10^6], no range check!)
"""
CBMicroseconds(secs::Integer, msecs::Integer = 0) = CBMicroseconds(@ccall libcb.cb_microseconds_new4(secs::Clong, msecs::Clong)::Ptr{Cvoid})

@doc raw"""
    CBMicroseconds(hours::Integer, minutes::Integer, secs::Integer, micros::Integer)

convert hours, minutes, secs, micros to Microseconds (no range check!)
"""
CBMicroseconds(hours::Integer, minutes::Integer, secs::Integer, micros::Integer) = CBMicroseconds(@ccall libcb.cb_microseconds_new6(hours::Clong, minutes::Clong, secs::Clong, micros::Clong)::Ptr{Cvoid})

Base.copy!(self::CBMicroseconds, m::CBMicroseconds) = (@ccall libcb.cb_microseconds_assign(self.data::Ptr{Cvoid}, m.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

MA.operate!(::typeof(+), self::CBMicroseconds, m::CBMicroseconds) = (@ccall libcb.cb_microseconds_plus(self.data::Ptr{Cvoid}, m.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

MA.promote_operation(::typeof(+), ::Type{<:CBMicroseconds}, ::Type{<:CBMicroseconds}) = CBMicroseconds

MA.operate!(::typeof(-), self::CBMicroseconds, m::CBMicroseconds) = (@ccall libcb.cb_microseconds_minus(self.data::Ptr{Cvoid}, m.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

MA.promote_operation(::typeof(-), ::Type{<:CBMicroseconds}, ::Type{<:CBMicroseconds}) = CBMicroseconds

Base.:<(self::CBMicroseconds, m::CBMicroseconds) = Bool(@ccall libcb.cb_microseconds_new_less(self.data::Ptr{Cvoid}, m.data::Ptr{Cvoid})::Bool)

Base.:>(self::CBMicroseconds, m::CBMicroseconds) = Bool(@ccall libcb.cb_microseconds_new_greater(self.data::Ptr{Cvoid}, m.data::Ptr{Cvoid})::Bool)

Base.:<=(self::CBMicroseconds, m::CBMicroseconds) = Bool(@ccall libcb.cb_microseconds_new_lessequal(self.data::Ptr{Cvoid}, m.data::Ptr{Cvoid})::Bool)

Base.:>=(self::CBMicroseconds, m::CBMicroseconds) = Bool(@ccall libcb.cb_microseconds_new_greaterequal(self.data::Ptr{Cvoid}, m.data::Ptr{Cvoid})::Bool)

Base.:(==)(self::CBMicroseconds, m::CBMicroseconds) = Bool(@ccall libcb.cb_microseconds_new_equal(self.data::Ptr{Cvoid}, m.data::Ptr{Cvoid})::Bool)

@doc raw"""
    cb_set_infinity!(self::CBMicroseconds, infty::Bool)

use true to regard value as infinity
"""
cb_set_infinity!(self::CBMicroseconds, infty::Bool) = @ccall libcb.cb_microseconds_set_infinity(self.data::Ptr{Cvoid}, infty::Cint)::Cvoid

@doc raw"""
    cb_get_infinity(self::CBMicroseconds)

if true, value should be regarded as infinity
"""
cb_get_infinity(self::CBMicroseconds) = Bool(@ccall libcb.cb_microseconds_get_infinity(self.data::Ptr{Cvoid})::Cint)

@doc raw"""
    cb_hhmmss(self::CBMicroseconds)

convert and store the value of (*this) to hours, minutes, seconds
"""
function cb_hhmmss(self::CBMicroseconds)
    secs = Ref{Integer}()
    minutes = Ref{Integer}()
    hours = Ref{Integer}()
    @ccall libcb.cb_microseconds_hhmmss(self.data::Ptr{Cvoid}, hours::Ref{Integer}, minutes::Ref{Integer}, secs::Ref{Integer})::Cvoid
    return hours[], minutes[], secs[]
end

@doc raw"""
    cb_hhmmssdd(self::CBMicroseconds)

convert and store the value of (*this) to hours, minutes, seconds, hundredths
"""
function cb_hhmmssdd(self::CBMicroseconds)
    hund = Ref{Integer}()
    secs = Ref{Integer}()
    minutes = Ref{Integer}()
    hours = Ref{Integer}()
    @ccall libcb.cb_microseconds_hhmmssdd(self.data::Ptr{Cvoid}, hours::Ref{Integer}, minutes::Ref{Integer}, secs::Ref{Integer}, hund::Ref{Integer})::Cvoid
    return hours[], minutes[], secs[], hund[]
end

@doc raw"""
    cb_roundsecs(self::CBMicroseconds)

round the value to seconds
"""
cb_roundsecs(self::CBMicroseconds) = @ccall libcb.cb_microseconds_roundsecs(self.data::Ptr{Cvoid})::Clong

@doc raw"""
    cb_roundhundredths(self::CBMicroseconds)

round the value to hundredths
"""
cb_roundhundredths(self::CBMicroseconds) = @ccall libcb.cb_microseconds_roundhundredths(self.data::Ptr{Cvoid})::Clong

@doc raw"""
    cb_print_time(m::CBMicroseconds, secondsonly::Integer = 0)

print Microseconds in the format "hh:mm:ss.dd" or "hh:mm:ss"
"""
cb_print_time(m::CBMicroseconds, secondsonly::Integer = 0) = @ccall libcb.cb_microseconds_print_time(m.data::Ptr{Cvoid}, secondsonly::Cint)::Cvoid

