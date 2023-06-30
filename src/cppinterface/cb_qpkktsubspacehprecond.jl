@doc raw"""
    CBQPKKTSubspaceHPrecond(inmethod::Integer = 0, cbinc::Integer = -1)

default constructor
"""
CBQPKKTSubspaceHPrecond(inmethod::Integer = 0, cbinc::Integer = -1) = CBQPKKTSubspaceHPrecond(@ccall libcb.cb_qpkktsubspacehprecond_new(inmethod::Cint, cbinc::Cint)::Ptr{Cvoid})

@doc raw"""
    cb_init_data!(self::CBQPKKTSubspaceHPrecond, Hp::Union{<:CBQPSolverProxObject,Nothing}, model::Union{<:CBQPModelBlockObject,Nothing}, A::Union{<:CBSparsemat,Nothing}, eq_indices::Union{<:CBIndexmatrix,Nothing}, SchurComplAineq::Bool)

returns 1 if this class is not applicable in the current data situation, otherwise it stores the data pointers and these need to stay valid throught the use of the other routines but are not deleted here
"""
cb_init_data!(self::CBQPKKTSubspaceHPrecond, Hp::Union{<:CBQPSolverProxObject,Nothing}, model::Union{<:CBQPModelBlockObject,Nothing}, A::Union{<:CBSparsemat,Nothing}, eq_indices::Union{<:CBIndexmatrix,Nothing}, SchurComplAineq::Bool) = @ccall libcb.cb_qpkktsubspacehprecond_init_data(self.data::Ptr{Cvoid}, (isnothing(Hp) ? C_NULL : Hp.data)::Ptr{Cvoid}, (isnothing(model) ? C_NULL : model.data)::Ptr{Cvoid}, (isnothing(A) ? C_NULL : A.data)::Ptr{Cvoid}, (isnothing(eq_indices) ? C_NULL : eq_indices.data)::Ptr{Cvoid}, SchurComplAineq::Cint)::Cint

@doc raw"""
    cb_init_system!(self::CBQPKKTSubspaceHPrecond, KKTdiagx::CBMatrix, KKTdiagy::CBMatrix, Hfactor::Real, prec::Real, params::Union{<:CBQPSolverParameters,Nothing})

set up the primal dual KKT system for being solved for predictor and corrector rhs; the input objects KKTdiagx and KKTdiagy will not change during use of the preconditioner, so it suffices to store the address if they are need during application of the preconditioner
"""
cb_init_system!(self::CBQPKKTSubspaceHPrecond, KKTdiagx::CBMatrix, KKTdiagy::CBMatrix, Hfactor::Real, prec::Real, params::Union{<:CBQPSolverParameters,Nothing}) = @ccall libcb.cb_qpkktsubspacehprecond_init_system(self.data::Ptr{Cvoid}, KKTdiagx.data::Ptr{Cvoid}, KKTdiagy.data::Ptr{Cvoid}, Hfactor::Cdouble, prec::Cdouble, (isnothing(params) ? C_NULL : params.data)::Ptr{Cvoid})::Cint

@doc raw"""
    cb_precondM1!(self::CBQPKKTSubspaceHPrecond, vec::CBMatrix)

returns M1^{-1}*vec; default: M1=I
"""
cb_precondM1!(self::CBQPKKTSubspaceHPrecond, vec::CBMatrix) = @ccall libcb.cb_qpkktsubspacehprecond_precondm1(self.data::Ptr{Cvoid}, vec.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_set_subspace!(self::CBQPKKTSubspaceHPrecond, insubspace::CBMatrix)

if the method admits this, let the subspace be chosen externally
"""
cb_set_subspace!(self::CBQPKKTSubspaceHPrecond, insubspace::CBMatrix) = @ccall libcb.cb_qpkktsubspacehprecond_set_subspace(self.data::Ptr{Cvoid}, insubspace.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_get_lmin_invM1!(self::CBQPKKTSubspaceHPrecond)

return (an estimate of) the minimum eigenvalue of the preconditioner M1^{-1}; this is used, e.g., to correct the precission in MINRES
"""
cb_get_lmin_invM1!(self::CBQPKKTSubspaceHPrecond) = @ccall libcb.cb_qpkktsubspacehprecond_get_lmin_invm1(self.data::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_precond_invG1!(self::CBQPKKTSubspaceHPrecond, vec::CBMatrix)

for estimating the condition number with M1=G*G^T this returns G^{-1}*vec; default: G=I
"""
cb_precond_invG1!(self::CBQPKKTSubspaceHPrecond, vec::CBMatrix) = @ccall libcb.cb_qpkktsubspacehprecond_precond_invg1(self.data::Ptr{Cvoid}, vec.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_precond_invG1tran!(self::CBQPKKTSubspaceHPrecond, vec::CBMatrix)

for estimating the condition number with M1=G*G^T this returns G^{-T}*vec; default: G=I
"""
cb_precond_invG1tran!(self::CBQPKKTSubspaceHPrecond, vec::CBMatrix) = @ccall libcb.cb_qpkktsubspacehprecond_precond_invg1tran(self.data::Ptr{Cvoid}, vec.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_precond_size!(self::CBQPKKTSubspaceHPrecond)

for estimating the condition number directly for the preconditioned part only; negative numbers indicate that the routine is not implemented
"""
cb_precond_size!(self::CBQPKKTSubspaceHPrecond) = @ccall libcb.cb_qpkktsubspacehprecond_precond_size(self.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_cond_number_mult!(self::CBQPKKTSubspaceHPrecond, vec::CBMatrix, KKTdiagx::CBMatrix, KKTdiagy::CBMatrix)

for estimating the condition number directly for the preconditioned part only
"""
cb_cond_number_mult!(self::CBQPKKTSubspaceHPrecond, vec::CBMatrix, KKTdiagx::CBMatrix, KKTdiagy::CBMatrix) = @ccall libcb.cb_qpkktsubspacehprecond_cond_number_mult(self.data::Ptr{Cvoid}, vec.data::Ptr{Cvoid}, KKTdiagx.data::Ptr{Cvoid}, KKTdiagy.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_get_precond_rank!(self::CBQPKKTSubspaceHPrecond)

for evaluation purposes with iterative solvers, return the rank of the precondiontioner used (or the number of n-vector multiplications per call)
"""
cb_get_precond_rank!(self::CBQPKKTSubspaceHPrecond) = @ccall libcb.cb_qpkktsubspacehprecond_get_precond_rank(self.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_get_t_precond_mult!(self::CBQPKKTSubspaceHPrecond)

for evaluation purposes with iterative solvers, return the time spent in the multiplication with
"""
cb_get_t_precond_mult!(self::CBQPKKTSubspaceHPrecond) = CBCH_Tools::Microseconds(@ccall libcb.cb_qpkktsubspacehprecond_new_get_t_precond_mult(self.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_reset_t_precond_mult!(self::CBQPKKTSubspaceHPrecond)

for evaluation purposes with iterative solvers, reset the time spent in the multiplication with the preconditioner to zero
"""
cb_reset_t_precond_mult!(self::CBQPKKTSubspaceHPrecond) = @ccall libcb.cb_qpkktsubspacehprecond_reset_t_precond_mult(self.data::Ptr{Cvoid})::Cvoid

