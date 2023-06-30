@doc raw"""
    CBCoeffmatInfo(sf::Real = 1.)

default value is 1 for no scaling
"""
CBCoeffmatInfo(sf::Real = 1.) = CBCoeffmatInfo(@ccall libcb.cb_coeffmatinfo_new(sf::Cdouble)::Ptr{Cvoid})

@doc raw"""
    cb_clone(self::CBCoeffmatInfo)

generates a new copy of itself on the heap
"""
cb_clone(self::CBCoeffmatInfo) = CBCoeffmatInfo(@ccall libcb.cb_coeffmatinfo_clone(self.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_get_scalefactor(self::CBCoeffmatInfo)

returns the scale factor
"""
cb_get_scalefactor(self::CBCoeffmatInfo) = @ccall libcb.cb_coeffmatinfo_get_scalefactor(self.data::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_set_scalefactor!(self::CBCoeffmatInfo, sf::Real)

sets the scale factor
"""
cb_set_scalefactor!(self::CBCoeffmatInfo, sf::Real) = @ccall libcb.cb_coeffmatinfo_set_scalefactor(self.data::Ptr{Cvoid}, sf::Cdouble)::Cvoid

@doc raw"""
    cb_multiply!(self::CBCoeffmatInfo, sf::Real)

scales the scale factor
"""
cb_multiply!(self::CBCoeffmatInfo, sf::Real) = @ccall libcb.cb_coeffmatinfo_multiply(self.data::Ptr{Cvoid}, sf::Cdouble)::Cvoid

@doc raw"""
    cb_print_id(self::CBCoeffmatInfo)

output a name if there is one for recognizing the type
"""
cb_print_id(self::CBCoeffmatInfo) = @ccall libcb.cb_coeffmatinfo_print_id(self.data::Ptr{Cvoid})::Cvoid

