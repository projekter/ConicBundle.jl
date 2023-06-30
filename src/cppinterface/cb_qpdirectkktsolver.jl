@doc raw"""
    cb_clear!(self::CBQPDirectKKTSolver)

reset data to empty
"""
cb_clear!(self::CBQPDirectKKTSolver) = @ccall libcb.cb_qpdirectkktsolver_clear(self.data::Ptr{Cvoid})::Cvoid

@doc raw"""
    CBQPDirectKKTSolver(in_factorize_ABC::Bool = false, cbinc::Integer = -1)

default constructor
"""
CBQPDirectKKTSolver(in_factorize_ABC::Bool = false, cbinc::Integer = -1) = CBQPDirectKKTSolver(@ccall libcb.cb_qpdirectkktsolver_new(in_factorize_ABC::Cint, cbinc::Cint)::Ptr{Cvoid})

@doc raw"""
    cb_QPinit_KKTdata!(self::CBQPDirectKKTSolver, Hp::Union{<:CBQPSolverProxObject,Nothing}, model::Union{<:CBQPModelBlockObject,Nothing}, A::Union{<:CBSparsemat,Nothing}, eq_indices::Union{<:CBIndexmatrix,Nothing})

returns 1 if this class is not applicable in the current data situation, otherwise it stores the data pointers and these need to stay valid throught the use of the other routines but are not deleted here
"""
cb_QPinit_KKTdata!(self::CBQPDirectKKTSolver, Hp::Union{<:CBQPSolverProxObject,Nothing}, model::Union{<:CBQPModelBlockObject,Nothing}, A::Union{<:CBSparsemat,Nothing}, eq_indices::Union{<:CBIndexmatrix,Nothing}) = @ccall libcb.cb_qpdirectkktsolver_qpinit_kktdata(self.data::Ptr{Cvoid}, (isnothing(Hp) ? C_NULL : Hp.data)::Ptr{Cvoid}, (isnothing(model) ? C_NULL : model.data)::Ptr{Cvoid}, (isnothing(A) ? C_NULL : A.data)::Ptr{Cvoid}, (isnothing(eq_indices) ? C_NULL : eq_indices.data)::Ptr{Cvoid})::Cint

@doc raw"""
    cb_QPinit_KKTsystem!(self::CBQPDirectKKTSolver, KKTdiagx::CBMatrix, KKTdiagy::CBMatrix, Hfactor::Real, prec::Real, params::Union{<:CBQPSolverParameters,Nothing})

set up the primal dual KKT system for being solved for predictor and corrector rhs in QPsolve_KKTsystem
"""
cb_QPinit_KKTsystem!(self::CBQPDirectKKTSolver, KKTdiagx::CBMatrix, KKTdiagy::CBMatrix, Hfactor::Real, prec::Real, params::Union{<:CBQPSolverParameters,Nothing}) = @ccall libcb.cb_qpdirectkktsolver_qpinit_kktsystem(self.data::Ptr{Cvoid}, KKTdiagx.data::Ptr{Cvoid}, KKTdiagy.data::Ptr{Cvoid}, Hfactor::Cdouble, prec::Cdouble, (isnothing(params) ? C_NULL : params.data)::Ptr{Cvoid})::Cint

@doc raw"""
    cb_QPsolve_KKTsystem!(self::CBQPDirectKKTSolver, solx::CBMatrix, soly::CBMatrix, primalrhs::CBMatrix, dualrhs::CBMatrix, rhsmu::Real, rhscorr::Real, prec::Real, params::Union{<:CBQPSolverParameters,Nothing})

solve the KKTsystem to precision prec for the given right hand sides that have been computed for the value rhsmu of the barrier parameter and in which a rhscorr fraction (out of [0,1] of the corrector term have been included; in iterative solvers solx and soly may be used as starting points
"""
cb_QPsolve_KKTsystem!(self::CBQPDirectKKTSolver, solx::CBMatrix, soly::CBMatrix, primalrhs::CBMatrix, dualrhs::CBMatrix, rhsmu::Real, rhscorr::Real, prec::Real, params::Union{<:CBQPSolverParameters,Nothing}) = @ccall libcb.cb_qpdirectkktsolver_qpsolve_kktsystem(self.data::Ptr{Cvoid}, solx.data::Ptr{Cvoid}, soly.data::Ptr{Cvoid}, primalrhs.data::Ptr{Cvoid}, dualrhs.data::Ptr{Cvoid}, rhsmu::Cdouble, rhscorr::Cdouble, prec::Cdouble, (isnothing(params) ? C_NULL : params.data)::Ptr{Cvoid})::Cint

