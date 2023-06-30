@doc raw"""
    CBGramSparsePSCPrimal(sps::CBSparsesym, factor::Real = 1.)

initialize to the given sparse symmetric matrix, the gram part is zero
"""
CBGramSparsePSCPrimal(sps::CBSparsesym, factor::Real = 1.) = CBGramSparsePSCPrimal(@ccall libcb.cb_gramsparsepscprimal_new(sps.data::Ptr{Cvoid}, factor::Cdouble)::Ptr{Cvoid})

@doc raw"""
    CBGramSparsePSCPrimal(pr::CBGramSparsePSCPrimal, factor::Real = 1.)

copy constructor
"""
CBGramSparsePSCPrimal(pr::CBGramSparsePSCPrimal, factor::Real = 1.) = CBGramSparsePSCPrimal(@ccall libcb.cb_gramsparsepscprimal_new2(pr.data::Ptr{Cvoid}, factor::Cdouble)::Ptr{Cvoid})

@doc raw"""
    Base.copy!(self::CBGramSparsePSCPrimal, sdp::CBSparsesym)

assign the sparse symmetric matrix to this and set the gram part to zero
"""
Base.copy!(self::CBGramSparsePSCPrimal, sdp::CBSparsesym) = (@ccall libcb.cb_gramsparsepscprimal_assign(self.data::Ptr{Cvoid}, sdp.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

@doc raw"""
    Base.copy!(self::CBGramSparsePSCPrimal, sdp::CBGramSparsePSCPrimal)

copy the information
"""
Base.copy!(self::CBGramSparsePSCPrimal, sdp::CBGramSparsePSCPrimal) = (@ccall libcb.cb_gramsparsepscprimal_assign2(self.data::Ptr{Cvoid}, sdp.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

@doc raw"""
    cb_get_grammatrix(self::CBGramSparsePSCPrimal)

returns the matrix \f$P\f$ giving rise to the Gram matrix \f$PP^T\f$
"""
cb_get_grammatrix(self::CBGramSparsePSCPrimal) = (@ccall libcb.cb_gramsparsepscprimal_get_grammatrix(self.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

@doc raw"""
    cb_assign_Gram_matrix!(self::CBGramSparsePSCPrimal, P::CBMatrix)

Set the grammatrix part to \f$P\f$ and set all values on the sparse support to zero (but keep this support even if it is zero now!)
"""
cb_assign_Gram_matrix!(self::CBGramSparsePSCPrimal, P::CBMatrix) = @ccall libcb.cb_gramsparsepscprimal_assign_gram_matrix(self.data::Ptr{Cvoid}, P.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_aggregate_primal_data!(self::CBGramSparsePSCPrimal, it::CBPrimalData, factor::Real = 1.)

*@brief if it is a GramSparsePSCPrimal or SparseSDPRimal, add factor*it to this on only the support of this sparsematrix

       Even if it is a GramSparsePSCPrimal and it has a nontirival Gram matrix part,
       this part is only added to the sparse part of this on the support of the
       sparse part of this. No attempt is made to enlarge the Gram part.
     
"""
cb_aggregate_primal_data!(self::CBGramSparsePSCPrimal, it::CBPrimalData, factor::Real = 1.) = @ccall libcb.cb_gramsparsepscprimal_aggregate_primal_data(self.data::Ptr{Cvoid}, it.data::Ptr{Cvoid}, factor::Cdouble)::Cint

@doc raw"""
    cb_aggregate_Gram_matrix!(self::CBGramSparsePSCPrimal, P::CBMatrix, factor::Real = 1.)

*@brief add factor*P*P^T to this, collecting
       all available information only in the sparse part, even the own Gram part.

       This operation more or less converts this to a SparsePSCPrimal. The point
       is that this routine is typically only called by PSCModel when aggregating
       the information in a single aggregate matrix over several steps and it is
       pointless to try to keep any gram information in this case.
     
"""
cb_aggregate_Gram_matrix!(self::CBGramSparsePSCPrimal, P::CBMatrix, factor::Real = 1.) = @ccall libcb.cb_gramsparsepscprimal_aggregate_gram_matrix(self.data::Ptr{Cvoid}, P.data::Ptr{Cvoid}, factor::Cdouble)::Cint

@doc raw"""
    cb_scale_primal_data!(self::CBGramSparsePSCPrimal, factor::Real)

multiply/scale *this with a nonnegative factor
"""
cb_scale_primal_data!(self::CBGramSparsePSCPrimal, factor::Real) = @ccall libcb.cb_gramsparsepscprimal_scale_primal_data(self.data::Ptr{Cvoid}, factor::Cdouble)::Cint

@doc raw"""
    cb_primal_ip(self::CBGramSparsePSCPrimal, A::CBSparseCoeffmatMatrix, column::Integer)

if compatible evaluate value=ip(*this,A.column[i])
"""
function cb_primal_ip(self::CBGramSparsePSCPrimal, A::CBSparseCoeffmatMatrix, column::Integer)
    value = Ref{Float64}()
    @ccall libcb.cb_gramsparsepscprimal_primal_ip(self.data::Ptr{Cvoid}, value::Ref{Float64}, A.data::Ptr{Cvoid}, column::Cint)::Cint
    return value[]
end

