@doc raw"""
    CBCFunctionMinorantExtender(fk::Union{<:AbstractVector{Nothing},Nothing}, se::Ptr{Cvoid})

constructor
"""
CBCFunctionMinorantExtender(fk::Union{<:AbstractVector{Nothing},Nothing}, se::Ptr{Cvoid}) = GC.@preserve fk begin
    (LinearAlgebra.chkstride1(fk); CBCFunctionMinorantExtender(@ccall libcb.cb_cfunctionminorantextender_new(fk::Ptr{Cvoid}, se::Ptr{Cvoid})::Ptr{Cvoid}))
end

@doc raw"""
    cb_extend!(self::CBCFunctionMinorantExtender, minorant::CBMinorant, n_coords::Integer, indices::Union{<:AbstractVector{Integer},Nothing})

see MinorantExtender::extend() for explanations
"""
cb_extend!(self::CBCFunctionMinorantExtender, minorant::CBMinorant, n_coords::Integer, indices::Union{<:AbstractVector{Integer},Nothing}) = GC.@preserve indices begin
    (LinearAlgebra.chkstride1(indices); @ccall libcb.cb_cfunctionminorantextender_extend(self.data::Ptr{Cvoid}, minorant.data::Ptr{Cvoid}, n_coords::Cint, indices::Ptr{Cint})::Cint)
end

