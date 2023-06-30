@doc raw"""
    CBDensePSCPrimal()

initialize to a symmetric matrix of size 0
"""
CBDensePSCPrimal() = CBDensePSCPrimal(@ccall libcb.cb_densepscprimal_new()::Ptr{Cvoid})

@doc raw"""
    CBDensePSCPrimal(symmat::CBDensePSCPrimal, factor::Real = 1.)

copy constructor
"""
CBDensePSCPrimal(symmat::CBDensePSCPrimal, factor::Real = 1.) = CBDensePSCPrimal(@ccall libcb.cb_densepscprimal_new2(symmat.data::Ptr{Cvoid}, factor::Cdouble)::Ptr{Cvoid})

@doc raw"""
    CBDensePSCPrimal(symmat::CBSymmatrix, factor::Real = 1.)

copy constructor from a CH_Matrix_Classes::Symmatrix
"""
CBDensePSCPrimal(symmat::CBSymmatrix, factor::Real = 1.) = CBDensePSCPrimal(@ccall libcb.cb_densepscprimal_new3(symmat.data::Ptr{Cvoid}, factor::Cdouble)::Ptr{Cvoid})

@doc raw"""
    CBDensePSCPrimal(n::Integer)

direct size initialization to a zero matrix
"""
CBDensePSCPrimal(n::Integer) = CBDensePSCPrimal(@ccall libcb.cb_densepscprimal_new4(n::Cint)::Ptr{Cvoid})

@doc raw"""
    Base.copy!(self::CBDensePSCPrimal, symmat::CBSymmatrix)

assigns a symmetric matrix
"""
Base.copy!(self::CBDensePSCPrimal, symmat::CBSymmatrix) = (@ccall libcb.cb_densepscprimal_assign(self.data::Ptr{Cvoid}, symmat.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

@doc raw"""
    cb_assign_Gram_matrix!(self::CBDensePSCPrimal, P::CBMatrix)

assign P*P^T to this
"""
cb_assign_Gram_matrix!(self::CBDensePSCPrimal, P::CBMatrix) = @ccall libcb.cb_densepscprimal_assign_gram_matrix(self.data::Ptr{Cvoid}, P.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_aggregate_primal_data!(self::CBDensePSCPrimal, it::CBPrimalData, factor::Real = 1.)

add factor*it to this (it must also be a DensePSCPrimal)
"""
cb_aggregate_primal_data!(self::CBDensePSCPrimal, it::CBPrimalData, factor::Real = 1.) = @ccall libcb.cb_densepscprimal_aggregate_primal_data(self.data::Ptr{Cvoid}, it.data::Ptr{Cvoid}, factor::Cdouble)::Cint

@doc raw"""
    cb_aggregate_Gram_matrix!(self::CBDensePSCPrimal, P::CBMatrix, factor::Real = 1.)

add factor*P*P^T to this
"""
cb_aggregate_Gram_matrix!(self::CBDensePSCPrimal, P::CBMatrix, factor::Real = 1.) = @ccall libcb.cb_densepscprimal_aggregate_gram_matrix(self.data::Ptr{Cvoid}, P.data::Ptr{Cvoid}, factor::Cdouble)::Cint

@doc raw"""
    cb_scale_primal_data!(self::CBDensePSCPrimal, factor::Real)

multiply/scale *this with a nonnegative factor
"""
cb_scale_primal_data!(self::CBDensePSCPrimal, factor::Real) = @ccall libcb.cb_densepscprimal_scale_primal_data(self.data::Ptr{Cvoid}, factor::Cdouble)::Cint

@doc raw"""
    cb_primal_ip(self::CBDensePSCPrimal, A::CBSparseCoeffmatMatrix, column::Integer)

if compatible evaluate value=ip(*this,A.column[i])
"""
function cb_primal_ip(self::CBDensePSCPrimal, A::CBSparseCoeffmatMatrix, column::Integer)
    value = Ref{Float64}()
    @ccall libcb.cb_densepscprimal_primal_ip(self.data::Ptr{Cvoid}, value::Ref{Float64}, A.data::Ptr{Cvoid}, column::Cint)::Cint
    return value[]
end

