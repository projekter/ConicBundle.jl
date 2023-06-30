@doc raw"""
    CBPsqmr(pril::Integer = -1)

default constructor
"""
CBPsqmr(pril::Integer = -1) = CBPsqmr(@ccall libcb.cb_psqmr_new(pril::Cint)::Ptr{Cvoid})

@doc raw"""
    cb_set_maxit!(self::CBPsqmr, in_maxit::Integer)

set maximum number of iterations
"""
cb_set_maxit!(self::CBPsqmr, in_maxit::Integer) = @ccall libcb.cb_psqmr_set_maxit(self.data::Ptr{Cvoid}, in_maxit::Cint)::Cvoid

@doc raw"""
    cb_get_maxit(self::CBPsqmr)

get maximum number of iterations
"""
cb_get_maxit(self::CBPsqmr) = @ccall libcb.cb_psqmr_get_maxit(self.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_get_err(self::CBPsqmr)

returns the error code of the last call
"""
cb_get_err(self::CBPsqmr) = @ccall libcb.cb_psqmr_get_err(self.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_get_nmult(self::CBPsqmr)

returns the number of matrix-vector multiplications of the last call
"""
cb_get_nmult(self::CBPsqmr) = @ccall libcb.cb_psqmr_get_nmult(self.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_get_residual_norm(self::CBPsqmr)

returns the residual norm of last call
"""
cb_get_residual_norm(self::CBPsqmr) = @ccall libcb.cb_psqmr_get_residual_norm(self.data::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_get_avg_reduction(self::CBPsqmr)

returns the average of the achieved reduction factor per iteration
"""
cb_get_avg_reduction(self::CBPsqmr) = @ccall libcb.cb_psqmr_get_avg_reduction(self.data::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_get_termprec(self::CBPsqmr)

return the (absolute) precision requirement for termination used in the last call
"""
cb_get_termprec(self::CBPsqmr) = @ccall libcb.cb_psqmr_get_termprec(self.data::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_compute!(self::CBPsqmr, system::CBIterativeSystemObject, x::CBMatrix, termprec::Real, storex::Union{<:CBMatrix,Nothing} = nothing, storestep::Integer = 0)

compute the solution for system into x with (absolute) residual precision termprec
"""
cb_compute!(self::CBPsqmr, system::CBIterativeSystemObject, x::CBMatrix, termprec::Real, storex::Union{<:CBMatrix,Nothing} = nothing, storestep::Integer = 0) = @ccall libcb.cb_psqmr_compute(self.data::Ptr{Cvoid}, system.data::Ptr{Cvoid}, x.data::Ptr{Cvoid}, termprec::Cdouble, (isnothing(storex) ? C_NULL : storex.data)::Ptr{Cvoid}, storestep::Cint)::Cint

@doc raw"""
    cb_set_out!(self::CBPsqmr, in_print_level::Integer = 1)

set output stream and level of detail of log output (for debugging)
"""
cb_set_out!(self::CBPsqmr, in_print_level::Integer = 1) = @ccall libcb.cb_psqmr_set_out(self.data::Ptr{Cvoid}, in_print_level::Cint)::Cvoid

