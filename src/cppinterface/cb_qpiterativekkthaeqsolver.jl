@doc raw"""
    CBQPIterativeKKTHAeqSolver(insolver::Union{<:CBIterativeSolverObject,Nothing}, inprecond::Union{<:CBQPKKTPrecondObject,Nothing} = nothing, cbinc::Integer = -1)

default constructor
"""
CBQPIterativeKKTHAeqSolver(insolver::Union{<:CBIterativeSolverObject,Nothing}, inprecond::Union{<:CBQPKKTPrecondObject,Nothing} = nothing, cbinc::Integer = -1) = CBQPIterativeKKTHAeqSolver(@ccall libcb.cb_qpiterativekkthaeqsolver_new((isnothing(insolver) ? C_NULL : insolver.data)::Ptr{Cvoid}, (isnothing(inprecond) ? C_NULL : inprecond.data)::Ptr{Cvoid}, cbinc::Cint)::Ptr{Cvoid})

@doc raw"""
    cb_QPinit_KKTdata!(self::CBQPIterativeKKTHAeqSolver, Hp::Union{<:CBQPSolverProxObject,Nothing}, model::Union{<:CBQPModelBlockObject,Nothing}, A::Union{<:CBSparsemat,Nothing}, eq_indices::Union{<:CBIndexmatrix,Nothing})

returns 1 if this class is not applicable in the current data situation, otherwise it stores the data pointers and these need to stay valid throught the use of the other routines but are not deleted here
"""
cb_QPinit_KKTdata!(self::CBQPIterativeKKTHAeqSolver, Hp::Union{<:CBQPSolverProxObject,Nothing}, model::Union{<:CBQPModelBlockObject,Nothing}, A::Union{<:CBSparsemat,Nothing}, eq_indices::Union{<:CBIndexmatrix,Nothing}) = @ccall libcb.cb_qpiterativekkthaeqsolver_qpinit_kktdata(self.data::Ptr{Cvoid}, (isnothing(Hp) ? C_NULL : Hp.data)::Ptr{Cvoid}, (isnothing(model) ? C_NULL : model.data)::Ptr{Cvoid}, (isnothing(A) ? C_NULL : A.data)::Ptr{Cvoid}, (isnothing(eq_indices) ? C_NULL : eq_indices.data)::Ptr{Cvoid})::Cint

@doc raw"""
    cb_QPsolve_KKTsystem!(self::CBQPIterativeKKTHAeqSolver, solx::CBMatrix, soly::CBMatrix, primalrhs::CBMatrix, dualrhs::CBMatrix, rhsmu::Real, rhscorr::Real, prec::Real, params::Union{<:CBQPSolverParameters,Nothing})

solve the KKTsystem to precision prec for the given right hand sides that have been computed for the value rhsmu of the barrier parameter and in which a rhscorr fraction (out of [0,1] of the corrector term have been included; in iterative solvers solx and soly may be used as starting points
"""
cb_QPsolve_KKTsystem!(self::CBQPIterativeKKTHAeqSolver, solx::CBMatrix, soly::CBMatrix, primalrhs::CBMatrix, dualrhs::CBMatrix, rhsmu::Real, rhscorr::Real, prec::Real, params::Union{<:CBQPSolverParameters,Nothing}) = @ccall libcb.cb_qpiterativekkthaeqsolver_qpsolve_kktsystem(self.data::Ptr{Cvoid}, solx.data::Ptr{Cvoid}, soly.data::Ptr{Cvoid}, primalrhs.data::Ptr{Cvoid}, dualrhs.data::Ptr{Cvoid}, rhsmu::Cdouble, rhscorr::Cdouble, prec::Cdouble, (isnothing(params) ? C_NULL : params.data)::Ptr{Cvoid})::Cint

@doc raw"""
    cb_ItSys_mult!(self::CBQPIterativeKKTHAeqSolver, in_vec::CBMatrix, out_vec::CBMatrix)

returns out_vec=(system matrix)*in_vec
"""
cb_ItSys_mult!(self::CBQPIterativeKKTHAeqSolver, in_vec::CBMatrix, out_vec::CBMatrix) = @ccall libcb.cb_qpiterativekkthaeqsolver_itsys_mult(self.data::Ptr{Cvoid}, in_vec.data::Ptr{Cvoid}, out_vec.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_QPget_system_size!(self::CBQPIterativeKKTHAeqSolver)

for evaluation purposes with iterative solvers, return the size of the system matrix
"""
cb_QPget_system_size!(self::CBQPIterativeKKTHAeqSolver) = @ccall libcb.cb_qpiterativekkthaeqsolver_qpget_system_size(self.data::Ptr{Cvoid})::Cint

