@doc raw"""
    CBPCG(pril::Integer = -1)

default constructor
"""
CBPCG(pril::Integer = -1) = CBPCG(@ccall libcb.cb_pcg_new(pril::Cint)::Ptr{Cvoid})

@doc raw"""
    cb_set_maxit!(self::CBPCG, in_maxit::Integer)

set maximum number of iterations
"""
cb_set_maxit!(self::CBPCG, in_maxit::Integer) = @ccall libcb.cb_pcg_set_maxit(self.data::Ptr{Cvoid}, in_maxit::Cint)::Cvoid

@doc raw"""
    cb_get_maxit(self::CBPCG)

get maximum number of iterations
"""
cb_get_maxit(self::CBPCG) = @ccall libcb.cb_pcg_get_maxit(self.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_get_err(self::CBPCG)

returns the error code of the last call
"""
cb_get_err(self::CBPCG) = @ccall libcb.cb_pcg_get_err(self.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_get_nmult(self::CBPCG)

returns the number of matrix-vector multiplications of the last call
"""
cb_get_nmult(self::CBPCG) = @ccall libcb.cb_pcg_get_nmult(self.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_get_residual_norm(self::CBPCG)

returns the residual norm of last call
"""
cb_get_residual_norm(self::CBPCG) = @ccall libcb.cb_pcg_get_residual_norm(self.data::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_get_avg_reduction(self::CBPCG)

returns the average of the achieved reduction factor per iteration
"""
cb_get_avg_reduction(self::CBPCG) = @ccall libcb.cb_pcg_get_avg_reduction(self.data::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_get_termprec(self::CBPCG)

return the (absolute) precision requirement for termination used in the last call
"""
cb_get_termprec(self::CBPCG) = @ccall libcb.cb_pcg_get_termprec(self.data::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_compute!(self::CBPCG, system::CBIterativeSystemObject, x::CBMatrix, termprec::Real, storex::Union{<:CBMatrix,Nothing} = nothing, storestep::Integer = 0)

compute the solution for system into x with (absolute) residual precision termprec
"""
cb_compute!(self::CBPCG, system::CBIterativeSystemObject, x::CBMatrix, termprec::Real, storex::Union{<:CBMatrix,Nothing} = nothing, storestep::Integer = 0) = @ccall libcb.cb_pcg_compute(self.data::Ptr{Cvoid}, system.data::Ptr{Cvoid}, x.data::Ptr{Cvoid}, termprec::Cdouble, (isnothing(storex) ? C_NULL : storex.data)::Ptr{Cvoid}, storestep::Cint)::Cint

@doc raw"""
    cb_set_out!(self::CBPCG, in_print_level::Integer = 1)

set output stream and level of detail of log output (for debugging)
"""
cb_set_out!(self::CBPCG, in_print_level::Integer = 1) = @ccall libcb.cb_pcg_set_out(self.data::Ptr{Cvoid}, in_print_level::Cint)::Cvoid

