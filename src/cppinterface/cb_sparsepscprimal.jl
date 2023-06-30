@doc raw"""
    CBSparsePSCPrimal(sps::CBSparsesym, factor::Real = 1.)

copy constructor from a CH_Matrix_Classes::Sparssym, only the support of this matrix  will be used in all Gram operations
"""
CBSparsePSCPrimal(sps::CBSparsesym, factor::Real = 1.) = CBSparsePSCPrimal(@ccall libcb.cb_sparsepscprimal_new(sps.data::Ptr{Cvoid}, factor::Cdouble)::Ptr{Cvoid})

@doc raw"""
    CBSparsePSCPrimal(pr::CBSparsePSCPrimal, factor::Real = 1.)

copy constructor, only the same support will be used in all Gram operations
"""
CBSparsePSCPrimal(pr::CBSparsePSCPrimal, factor::Real = 1.) = CBSparsePSCPrimal(@ccall libcb.cb_sparsepscprimal_new2(pr.data::Ptr{Cvoid}, factor::Cdouble)::Ptr{Cvoid})

@doc raw"""
    Base.copy!(self::CBSparsePSCPrimal, sdp::CBSparsesym)

assigns this Sparsesym, only the same support will be used in all Gram operations
"""
Base.copy!(self::CBSparsePSCPrimal, sdp::CBSparsesym) = (@ccall libcb.cb_sparsepscprimal_assign(self.data::Ptr{Cvoid}, sdp.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

@doc raw"""
    cb_assign_Gram_matrix!(self::CBSparsePSCPrimal, P::CBMatrix)

for each element aij in the support set aij=<P.row(i),P.row(j)>
"""
cb_assign_Gram_matrix!(self::CBSparsePSCPrimal, P::CBMatrix) = @ccall libcb.cb_sparsepscprimal_assign_gram_matrix(self.data::Ptr{Cvoid}, P.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_aggregate_primal_data!(self::CBSparsePSCPrimal, it::CBPrimalData, factor::Real = 1.)

if it is a SparseSDPRimal, add factor*it to this on the support of this
"""
cb_aggregate_primal_data!(self::CBSparsePSCPrimal, it::CBPrimalData, factor::Real = 1.) = @ccall libcb.cb_sparsepscprimal_aggregate_primal_data(self.data::Ptr{Cvoid}, it.data::Ptr{Cvoid}, factor::Cdouble)::Cint

@doc raw"""
    cb_aggregate_Gram_matrix!(self::CBSparsePSCPrimal, P::CBMatrix, factor::Real = 1.)

add factor*P*P^T on the support to this
"""
cb_aggregate_Gram_matrix!(self::CBSparsePSCPrimal, P::CBMatrix, factor::Real = 1.) = @ccall libcb.cb_sparsepscprimal_aggregate_gram_matrix(self.data::Ptr{Cvoid}, P.data::Ptr{Cvoid}, factor::Cdouble)::Cint

@doc raw"""
    cb_scale_primal_data!(self::CBSparsePSCPrimal, factor::Real)

multiply/scale *this with a nonnegative factor
"""
cb_scale_primal_data!(self::CBSparsePSCPrimal, factor::Real) = @ccall libcb.cb_sparsepscprimal_scale_primal_data(self.data::Ptr{Cvoid}, factor::Cdouble)::Cint

@doc raw"""
    cb_primal_ip(self::CBSparsePSCPrimal, A::CBSparseCoeffmatMatrix, column::Integer)

if compatible evaluate value=ip(*this,A.column[i])
"""
function cb_primal_ip(self::CBSparsePSCPrimal, A::CBSparseCoeffmatMatrix, column::Integer)
    value = Ref{Float64}()
    @ccall libcb.cb_sparsepscprimal_primal_ip(self.data::Ptr{Cvoid}, value::Ref{Float64}, A.data::Ptr{Cvoid}, column::Cint)::Cint
    return value[]
end

