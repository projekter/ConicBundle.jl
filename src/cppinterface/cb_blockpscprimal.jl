@doc raw"""
    CBBlockPSCPrimal(pr::CBBlockPSCPrimal, factor::Real = 1.)

copy constructor
"""
CBBlockPSCPrimal(pr::CBBlockPSCPrimal, factor::Real = 1.) = CBBlockPSCPrimal(@ccall libcb.cb_blockpscprimal_new(pr.data::Ptr{Cvoid}, factor::Cdouble)::Ptr{Cvoid})

@doc raw"""
    cb_assign_Gram_matrix!(self::CBBlockPSCPrimal, P::CBMatrix)

*@brief for each element aij in the support set aij=<P.row(i),P.row(j)>

       The rows of P corresponding to each block are passed on the
       with a corresponding call to the PSCPrimal of this block
    
"""
cb_assign_Gram_matrix!(self::CBBlockPSCPrimal, P::CBMatrix) = @ccall libcb.cb_blockpscprimal_assign_gram_matrix(self.data::Ptr{Cvoid}, P.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_aggregate_primal_data!(self::CBBlockPSCPrimal, it::CBPrimalData, factor::Real = 1.)

* add factor*it to this
        (it must also be a BlockPSCPrimal with the same block partition) 
"""
cb_aggregate_primal_data!(self::CBBlockPSCPrimal, it::CBPrimalData, factor::Real = 1.) = @ccall libcb.cb_blockpscprimal_aggregate_primal_data(self.data::Ptr{Cvoid}, it.data::Ptr{Cvoid}, factor::Cdouble)::Cint

@doc raw"""
    cb_aggregate_Gram_matrix!(self::CBBlockPSCPrimal, P::CBMatrix, factor::Real = 1.)

*@brief add factor*P*P^T to this

       The rows of P corresponding to each block are passed on the
       with a corresponding call to the PSCPrimal of this block
    
"""
cb_aggregate_Gram_matrix!(self::CBBlockPSCPrimal, P::CBMatrix, factor::Real = 1.) = @ccall libcb.cb_blockpscprimal_aggregate_gram_matrix(self.data::Ptr{Cvoid}, P.data::Ptr{Cvoid}, factor::Cdouble)::Cint

@doc raw"""
    cb_get_nblocks(self::CBBlockPSCPrimal)

returns the number of blocks
"""
cb_get_nblocks(self::CBBlockPSCPrimal) = @ccall libcb.cb_blockpscprimal_get_nblocks(self.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_blockdim(self::CBBlockPSCPrimal, i::Integer)

returns the size of block i
"""
cb_blockdim(self::CBBlockPSCPrimal, i::Integer) = @ccall libcb.cb_blockpscprimal_blockdim(self.data::Ptr{Cvoid}, i::Cint)::Cint

@doc raw"""
    cb_block(self::CBBlockPSCPrimal, i::Integer)

returns the PSCPrimal of block i if there is one, 0 otherwise
"""
cb_block(self::CBBlockPSCPrimal, i::Integer) = CBPSCPrimal(@ccall libcb.cb_blockpscprimal_block(self.data::Ptr{Cvoid}, i::Cint)::Ptr{Cvoid})

@doc raw"""
    cb_scale_primal_data!(self::CBBlockPSCPrimal, factor::Real)

multiply/scale *this with a nonnegative factor
"""
cb_scale_primal_data!(self::CBBlockPSCPrimal, factor::Real) = @ccall libcb.cb_blockpscprimal_scale_primal_data(self.data::Ptr{Cvoid}, factor::Cdouble)::Cint

@doc raw"""
    cb_primal_ip(self::CBBlockPSCPrimal, A::CBSparseCoeffmatMatrix, column::Integer)

if compatible evaluate value=ip(*this,A.column[i])
"""
function cb_primal_ip(self::CBBlockPSCPrimal, A::CBSparseCoeffmatMatrix, column::Integer)
    value = Ref{Float64}()
    @ccall libcb.cb_blockpscprimal_primal_ip(self.data::Ptr{Cvoid}, value::Ref{Float64}, A.data::Ptr{Cvoid}, column::Cint)::Cint
    return value[]
end

