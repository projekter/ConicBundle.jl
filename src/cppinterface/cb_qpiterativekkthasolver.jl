@doc raw"""
    CBQPIterativeKKTHASolver(insolver::Union{<:CBIterativeSolverObject,Nothing}, inprecond::Union{<:CBQPKKTPrecondObject,Nothing} = nothing, cbinc::Integer = -1)

default constructor
"""
CBQPIterativeKKTHASolver(insolver::Union{<:CBIterativeSolverObject,Nothing}, inprecond::Union{<:CBQPKKTPrecondObject,Nothing} = nothing, cbinc::Integer = -1) = CBQPIterativeKKTHASolver(@ccall libcb.cb_qpiterativekkthasolver_new((isnothing(insolver) ? C_NULL : insolver.data)::Ptr{Cvoid}, (isnothing(inprecond) ? C_NULL : inprecond.data)::Ptr{Cvoid}, cbinc::Cint)::Ptr{Cvoid})

@doc raw"""
    cb_QPsolve_KKTsystem!(self::CBQPIterativeKKTHASolver, solx::CBMatrix, soly::CBMatrix, primalrhs::CBMatrix, dualrhs::CBMatrix, rhsmu::Real, rhscorr::Real, prec::Real, params::Union{<:CBQPSolverParameters,Nothing})

solve the KKTsystem to precision prec for the given right hand sides that have been computed for the value rhsmu of the barrier parameter and in which a rhscorr fraction (out of [0,1] of the corrector term have been included; in iterative solvers solx and soly may be used as starting points
"""
cb_QPsolve_KKTsystem!(self::CBQPIterativeKKTHASolver, solx::CBMatrix, soly::CBMatrix, primalrhs::CBMatrix, dualrhs::CBMatrix, rhsmu::Real, rhscorr::Real, prec::Real, params::Union{<:CBQPSolverParameters,Nothing}) = @ccall libcb.cb_qpiterativekkthasolver_qpsolve_kktsystem(self.data::Ptr{Cvoid}, solx.data::Ptr{Cvoid}, soly.data::Ptr{Cvoid}, primalrhs.data::Ptr{Cvoid}, dualrhs.data::Ptr{Cvoid}, rhsmu::Cdouble, rhscorr::Cdouble, prec::Cdouble, (isnothing(params) ? C_NULL : params.data)::Ptr{Cvoid})::Cint

@doc raw"""
    cb_ItSys_mult!(self::CBQPIterativeKKTHASolver, in_vec::CBMatrix, out_vec::CBMatrix)

returns out_vec=(system matrix)*in_vec
"""
cb_ItSys_mult!(self::CBQPIterativeKKTHASolver, in_vec::CBMatrix, out_vec::CBMatrix) = @ccall libcb.cb_qpiterativekkthasolver_itsys_mult(self.data::Ptr{Cvoid}, in_vec.data::Ptr{Cvoid}, out_vec.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_QPget_system_size!(self::CBQPIterativeKKTHASolver)

for evaluation purposes with iterative solvers, return the size of the system matrix
"""
cb_QPget_system_size!(self::CBQPIterativeKKTHASolver) = @ccall libcb.cb_qpiterativekkthasolver_qpget_system_size(self.data::Ptr{Cvoid})::Cint

