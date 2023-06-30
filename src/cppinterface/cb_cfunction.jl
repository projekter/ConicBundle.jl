@doc raw"""
    CBCFunction(fk::Union{<:AbstractVector{Nothing},Nothing}, fp::Ptr{Cvoid}, se::Ptr{Cvoid} = 0, prdim::Integer = 0)

constructor
"""
CBCFunction(fk::Union{<:AbstractVector{Nothing},Nothing}, fp::Ptr{Cvoid}, se::Ptr{Cvoid} = 0, prdim::Integer = 0) = GC.@preserve fk begin
    (LinearAlgebra.chkstride1(fk); CBCFunction(@ccall libcb.cb_cfunction_new(fk::Ptr{Cvoid}, fp::Ptr{Cvoid}, se::Ptr{Cvoid}, prdim::Cint)::Ptr{Cvoid}))
end

@doc raw"""
    cb_set_max_new!(self::CBCFunction, mn::Integer)

set the maximum number of new subgardients per evaluations
"""
cb_set_max_new!(self::CBCFunction, mn::Integer) = @ccall libcb.cb_cfunction_set_max_new(self.data::Ptr{Cvoid}, mn::Cint)::Cvoid

