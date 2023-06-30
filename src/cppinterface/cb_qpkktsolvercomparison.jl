@doc raw"""
    cb_clear!(self::CBQPKKTSolverComparison)

reset data to empty
"""
cb_clear!(self::CBQPKKTSolverComparison) = @ccall libcb.cb_qpkktsolvercomparison_clear(self.data::Ptr{Cvoid})::Cvoid

@doc raw"""
    CBQPKKTSolverComparison(cbinc::Integer = -1)

default constructor
"""
CBQPKKTSolverComparison(cbinc::Integer = -1) = CBQPKKTSolverComparison(@ccall libcb.cb_qpkktsolvercomparison_new(cbinc::Cint)::Ptr{Cvoid})

@doc raw"""
    cb_add_solver!(self::CBQPKKTSolverComparison, solver::Union{<:CBQPKKTSolverObject,Nothing}, name::Union{<:AbstractVector{UInt8},Nothing})

the first solver added is the reference solver
"""
cb_add_solver!(self::CBQPKKTSolverComparison, solver::Union{<:CBQPKKTSolverObject,Nothing}, name::Union{<:AbstractVector{UInt8},Nothing}) = GC.@preserve name begin
    (LinearAlgebra.chkstride1(name); @ccall libcb.cb_qpkktsolvercomparison_add_solver(self.data::Ptr{Cvoid}, (isnothing(solver) ? C_NULL : solver.data)::Ptr{Cvoid}, name::Ptr{Cchar})::Cint)
end

@doc raw"""
    cb_QPinit_KKTdata!(self::CBQPKKTSolverComparison, Hp::Union{<:CBQPSolverProxObject,Nothing}, model::Union{<:CBQPModelBlockObject,Nothing}, A::Union{<:CBSparsemat,Nothing}, eq_indices::Union{<:CBIndexmatrix,Nothing})

returns 1 if this class is not applicable in the current data situation, otherwise it stores the data pointers and these need to stay valid throught the use of the other routines but are not deleted here
"""
cb_QPinit_KKTdata!(self::CBQPKKTSolverComparison, Hp::Union{<:CBQPSolverProxObject,Nothing}, model::Union{<:CBQPModelBlockObject,Nothing}, A::Union{<:CBSparsemat,Nothing}, eq_indices::Union{<:CBIndexmatrix,Nothing}) = @ccall libcb.cb_qpkktsolvercomparison_qpinit_kktdata(self.data::Ptr{Cvoid}, (isnothing(Hp) ? C_NULL : Hp.data)::Ptr{Cvoid}, (isnothing(model) ? C_NULL : model.data)::Ptr{Cvoid}, (isnothing(A) ? C_NULL : A.data)::Ptr{Cvoid}, (isnothing(eq_indices) ? C_NULL : eq_indices.data)::Ptr{Cvoid})::Cint

@doc raw"""
    cb_QPinit_KKTsystem!(self::CBQPKKTSolverComparison, KKTdiagx::CBMatrix, KKTdiagy::CBMatrix, Hfactor::Real, prec::Real, params::Union{<:CBQPSolverParameters,Nothing})

set up the primal dual KKT system for being solved for predictor and corrector rhs in QPsolve_KKTsystem
"""
cb_QPinit_KKTsystem!(self::CBQPKKTSolverComparison, KKTdiagx::CBMatrix, KKTdiagy::CBMatrix, Hfactor::Real, prec::Real, params::Union{<:CBQPSolverParameters,Nothing}) = @ccall libcb.cb_qpkktsolvercomparison_qpinit_kktsystem(self.data::Ptr{Cvoid}, KKTdiagx.data::Ptr{Cvoid}, KKTdiagy.data::Ptr{Cvoid}, Hfactor::Cdouble, prec::Cdouble, (isnothing(params) ? C_NULL : params.data)::Ptr{Cvoid})::Cint

@doc raw"""
    cb_QPsolve_KKTsystem!(self::CBQPKKTSolverComparison, solx::CBMatrix, soly::CBMatrix, primalrhs::CBMatrix, dualrhs::CBMatrix, rhsmu::Real, rhscorr::Real, prec::Real, params::Union{<:CBQPSolverParameters,Nothing})

solve the KKTsystem to precision prec for the given right hand sides that have been computed for the value rhsmu of the barrier parameter and in which a rhscorr fraction (out of [0,1] of the corrector term have been included; in iterative solvers solx and soly may be used as starting points
"""
cb_QPsolve_KKTsystem!(self::CBQPKKTSolverComparison, solx::CBMatrix, soly::CBMatrix, primalrhs::CBMatrix, dualrhs::CBMatrix, rhsmu::Real, rhscorr::Real, prec::Real, params::Union{<:CBQPSolverParameters,Nothing}) = @ccall libcb.cb_qpkktsolvercomparison_qpsolve_kktsystem(self.data::Ptr{Cvoid}, solx.data::Ptr{Cvoid}, soly.data::Ptr{Cvoid}, primalrhs.data::Ptr{Cvoid}, dualrhs.data::Ptr{Cvoid}, rhsmu::Cdouble, rhscorr::Cdouble, prec::Cdouble, (isnothing(params) ? C_NULL : params.data)::Ptr{Cvoid})::Cint

@doc raw"""
    cb_QPget_blockH_norm!(self::CBQPKKTSolverComparison)

for judging violation this returns (an estimate of) the norm of the H-row in the latest system
"""
cb_QPget_blockH_norm!(self::CBQPKKTSolverComparison) = @ccall libcb.cb_qpkktsolvercomparison_qpget_blockh_norm(self.data::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_QPget_blockA_norm!(self::CBQPKKTSolverComparison)

for judging violation this returns (an estimate of) the norm of the A-row in the latest system
"""
cb_QPget_blockA_norm!(self::CBQPKKTSolverComparison) = @ccall libcb.cb_qpkktsolvercomparison_qpget_blocka_norm(self.data::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_get_mu_stats!(self::CBQPKKTSolverComparison, lbmu::Real, ubmu::Real, dims::CBIndexmatrix, mu::CBMatrix, prepsecs::CBMatrix, predsecs::CBMatrix, corrsecs::CBMatrix, predcalls::CBIndexmatrix, corrcalls::CBIndexmatrix, cond::CBMatrix, pccols::CBIndexmatrix, sysviol::CBMatrix)

return those data columns (each a KKT system; columns are more efficient to append than lines) that fall into the given lower and upper bounds on mu
"""
cb_get_mu_stats!(self::CBQPKKTSolverComparison, lbmu::Real, ubmu::Real, dims::CBIndexmatrix, mu::CBMatrix, prepsecs::CBMatrix, predsecs::CBMatrix, corrsecs::CBMatrix, predcalls::CBIndexmatrix, corrcalls::CBIndexmatrix, cond::CBMatrix, pccols::CBIndexmatrix, sysviol::CBMatrix) = @ccall libcb.cb_qpkktsolvercomparison_get_mu_stats(self.data::Ptr{Cvoid}, lbmu::Cdouble, ubmu::Cdouble, dims.data::Ptr{Cvoid}, mu.data::Ptr{Cvoid}, prepsecs.data::Ptr{Cvoid}, predsecs.data::Ptr{Cvoid}, corrsecs.data::Ptr{Cvoid}, predcalls.data::Ptr{Cvoid}, corrcalls.data::Ptr{Cvoid}, cond.data::Ptr{Cvoid}, pccols.data::Ptr{Cvoid}, sysviol.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_get_prob_stats!(self::CBQPKKTSolverComparison, dims::CBIndexmatrix, iterations::CBIndexmatrix, lastmu::CBMatrix, prepsecs::CBMatrix, predsecs::CBMatrix, corrsecs::CBMatrix, predcalls::CBIndexmatrix, corrcalls::CBIndexmatrix)

return one data column per subproblem (more efficient to append than lines) with the sum of the time/calls/etc.
"""
cb_get_prob_stats!(self::CBQPKKTSolverComparison, dims::CBIndexmatrix, iterations::CBIndexmatrix, lastmu::CBMatrix, prepsecs::CBMatrix, predsecs::CBMatrix, corrsecs::CBMatrix, predcalls::CBIndexmatrix, corrcalls::CBIndexmatrix) = @ccall libcb.cb_qpkktsolvercomparison_get_prob_stats(self.data::Ptr{Cvoid}, dims.data::Ptr{Cvoid}, iterations.data::Ptr{Cvoid}, lastmu.data::Ptr{Cvoid}, prepsecs.data::Ptr{Cvoid}, predsecs.data::Ptr{Cvoid}, corrsecs.data::Ptr{Cvoid}, predcalls.data::Ptr{Cvoid}, corrcalls.data::Ptr{Cvoid})::Cint

