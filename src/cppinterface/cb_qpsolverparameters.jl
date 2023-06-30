@doc raw"""
    CBQPSolverParameters(incr::Integer = -1)

default constructor
"""
CBQPSolverParameters(incr::Integer = -1) = CBQPSolverParameters(@ccall libcb.cb_qpsolverparameters_new(incr::Cint)::Ptr{Cvoid})

@doc raw"""
    cb_QPget_min_objective_relprec(self::CBQPSolverParameters)

get this variable value
"""
cb_QPget_min_objective_relprec(self::CBQPSolverParameters) = @ccall libcb.cb_qpsolverparameters_qpget_min_objective_relprec(self.data::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_QPget_objective_gap_eps(self::CBQPSolverParameters)

get this variable value
"""
cb_QPget_objective_gap_eps(self::CBQPSolverParameters) = @ccall libcb.cb_qpsolverparameters_qpget_objective_gap_eps(self.data::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_QPget_primal_infeasibility_eps(self::CBQPSolverParameters)

get this variable value
"""
cb_QPget_primal_infeasibility_eps(self::CBQPSolverParameters) = @ccall libcb.cb_qpsolverparameters_qpget_primal_infeasibility_eps(self.data::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_QPget_dual_infeasibility_eps(self::CBQPSolverParameters)

get this variable value
"""
cb_QPget_dual_infeasibility_eps(self::CBQPSolverParameters) = @ccall libcb.cb_qpsolverparameters_qpget_dual_infeasibility_eps(self.data::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_QPget_lower_bound(self::CBQPSolverParameters)

get this variable value
"""
cb_QPget_lower_bound(self::CBQPSolverParameters) = @ccall libcb.cb_qpsolverparameters_qpget_lower_bound(self.data::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_QPget_upper_bound(self::CBQPSolverParameters)

get this variable value
"""
cb_QPget_upper_bound(self::CBQPSolverParameters) = @ccall libcb.cb_qpsolverparameters_qpget_upper_bound(self.data::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_QPget_lower_bound_gap_eps(self::CBQPSolverParameters)

get this variable value
"""
cb_QPget_lower_bound_gap_eps(self::CBQPSolverParameters) = @ccall libcb.cb_qpsolverparameters_qpget_lower_bound_gap_eps(self.data::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_QPget_upper_bound_gap_eps(self::CBQPSolverParameters)

get this variable value
"""
cb_QPget_upper_bound_gap_eps(self::CBQPSolverParameters) = @ccall libcb.cb_qpsolverparameters_qpget_upper_bound_gap_eps(self.data::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_QPget_maxiter(self::CBQPSolverParameters)

get this variable value
"""
cb_QPget_maxiter(self::CBQPSolverParameters) = @ccall libcb.cb_qpsolverparameters_qpget_maxiter(self.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_QPget_KKTsolver!(self::CBQPSolverParameters)

get this variable value
"""
cb_QPget_KKTsolver!(self::CBQPSolverParameters) = CBQPKKTSolverObject(@ccall libcb.cb_qpsolverparameters_qpget_kktsolver(self.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_QPget_use_predictor_corrector(self::CBQPSolverParameters)

get this variable value
"""
cb_QPget_use_predictor_corrector(self::CBQPSolverParameters) = Bool(@ccall libcb.cb_qpsolverparameters_qpget_use_predictor_corrector(self.data::Ptr{Cvoid})::Cint)

@doc raw"""
    cb_QPget_use_neighborhood(self::CBQPSolverParameters)

get this variable value
"""
cb_QPget_use_neighborhood(self::CBQPSolverParameters) = Bool(@ccall libcb.cb_qpsolverparameters_qpget_use_neighborhood(self.data::Ptr{Cvoid})::Cint)

@doc raw"""
    cb_QPget_nbh_ub(self::CBQPSolverParameters)

get this variable value
"""
cb_QPget_nbh_ub(self::CBQPSolverParameters) = @ccall libcb.cb_qpsolverparameters_qpget_nbh_ub(self.data::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_QPget_nbh_lb(self::CBQPSolverParameters)

get this variable value
"""
cb_QPget_nbh_lb(self::CBQPSolverParameters) = @ccall libcb.cb_qpsolverparameters_qpget_nbh_lb(self.data::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_QPget_use_socqp(self::CBQPSolverParameters)

get this variable value
"""
cb_QPget_use_socqp(self::CBQPSolverParameters) = Bool(@ccall libcb.cb_qpsolverparameters_qpget_use_socqp(self.data::Ptr{Cvoid})::Cint)

@doc raw"""
    cb_QPset_min_objective_relprec!(self::CBQPSolverParameters, eps::Real)

set this variable value
"""
cb_QPset_min_objective_relprec!(self::CBQPSolverParameters, eps::Real) = @ccall libcb.cb_qpsolverparameters_qpset_min_objective_relprec(self.data::Ptr{Cvoid}, eps::Cdouble)::Cint

@doc raw"""
    cb_QPset_objective_gap_eps!(self::CBQPSolverParameters, eps::Real)

set this variable value
"""
cb_QPset_objective_gap_eps!(self::CBQPSolverParameters, eps::Real) = @ccall libcb.cb_qpsolverparameters_qpset_objective_gap_eps(self.data::Ptr{Cvoid}, eps::Cdouble)::Cint

@doc raw"""
    cb_QPset_primal_infeasibility_eps!(self::CBQPSolverParameters, eps::Real)

set this variable value
"""
cb_QPset_primal_infeasibility_eps!(self::CBQPSolverParameters, eps::Real) = @ccall libcb.cb_qpsolverparameters_qpset_primal_infeasibility_eps(self.data::Ptr{Cvoid}, eps::Cdouble)::Cint

@doc raw"""
    cb_QPset_dual_infeasibility_eps!(self::CBQPSolverParameters, eps::Real)

set this variable value
"""
cb_QPset_dual_infeasibility_eps!(self::CBQPSolverParameters, eps::Real) = @ccall libcb.cb_qpsolverparameters_qpset_dual_infeasibility_eps(self.data::Ptr{Cvoid}, eps::Cdouble)::Cint

@doc raw"""
    cb_QPset_lower_and_upper_bounds!(self::CBQPSolverParameters, lb::Real, ub::Real)

set this variable value
"""
cb_QPset_lower_and_upper_bounds!(self::CBQPSolverParameters, lb::Real, ub::Real) = @ccall libcb.cb_qpsolverparameters_qpset_lower_and_upper_bounds(self.data::Ptr{Cvoid}, lb::Cdouble, ub::Cdouble)::Cint

@doc raw"""
    cb_QPset_lower_bound_gap_eps!(self::CBQPSolverParameters, eps::Real)

set this variable value
"""
cb_QPset_lower_bound_gap_eps!(self::CBQPSolverParameters, eps::Real) = @ccall libcb.cb_qpsolverparameters_qpset_lower_bound_gap_eps(self.data::Ptr{Cvoid}, eps::Cdouble)::Cint

@doc raw"""
    cb_QPset_upper_bound_gap_eps!(self::CBQPSolverParameters, eps::Real)

set this variable value
"""
cb_QPset_upper_bound_gap_eps!(self::CBQPSolverParameters, eps::Real) = @ccall libcb.cb_qpsolverparameters_qpset_upper_bound_gap_eps(self.data::Ptr{Cvoid}, eps::Cdouble)::Cint

@doc raw"""
    cb_QPset_maxiter!(self::CBQPSolverParameters, mi::Integer)

set this variable value
"""
cb_QPset_maxiter!(self::CBQPSolverParameters, mi::Integer) = @ccall libcb.cb_qpsolverparameters_qpset_maxiter(self.data::Ptr{Cvoid}, mi::Cint)::Cint

@doc raw"""
    cb_QPset_KKTsolver!(self::CBQPSolverParameters, in_KKTsolver::Union{<:CBQPKKTSolverObject,Nothing})

delete previous solver and replace by the new one (should not be zero when calling the solver)
"""
cb_QPset_KKTsolver!(self::CBQPSolverParameters, in_KKTsolver::Union{<:CBQPKKTSolverObject,Nothing}) = @ccall libcb.cb_qpsolverparameters_qpset_kktsolver(self.data::Ptr{Cvoid}, (isnothing(in_KKTsolver) ? C_NULL : in_KKTsolver.data)::Ptr{Cvoid})::Cint

@doc raw"""
    cb_QPset_use_predictor_corrector!(self::CBQPSolverParameters, upc::Bool)

if set to true (=default), a predictor corrector approach is used for solving the KKT system (solve for barrier parameter mu=0, guess mu, solve again for this mu including some bilinear perturbation) otherwise the barrier parameter is set apriori and the step computed in one solve
"""
cb_QPset_use_predictor_corrector!(self::CBQPSolverParameters, upc::Bool) = @ccall libcb.cb_qpsolverparameters_qpset_use_predictor_corrector(self.data::Ptr{Cvoid}, upc::Cint)::Cint

@doc raw"""
    cb_QPset_use_neighborhood!(self::CBQPSolverParameters, nbh::Bool)

if set to true (default: false), the lines search is carried out with respect to the neighborhood polynomial ensuring the afte this step the point is again inside the neighborhood of the central path
"""
cb_QPset_use_neighborhood!(self::CBQPSolverParameters, nbh::Bool) = @ccall libcb.cb_qpsolverparameters_qpset_use_neighborhood(self.data::Ptr{Cvoid}, nbh::Cint)::Cint

@doc raw"""
    cb_QPset_nbh_bounds!(self::CBQPSolverParameters, nbhlb::Real, nbhub::Real)

set the upper bound on the neighborhood that should be ensured in curve searches; ensures eps_Real<=nbhlb<=nbhub (nbhub should be < 1. and <=.35 is safe)
"""
cb_QPset_nbh_bounds!(self::CBQPSolverParameters, nbhlb::Real, nbhub::Real) = @ccall libcb.cb_qpsolverparameters_qpset_nbh_bounds(self.data::Ptr{Cvoid}, nbhlb::Cdouble, nbhub::Cdouble)::Cint

@doc raw"""
    cb_QPset_use_socqp!(self::CBQPSolverParameters, s::Bool)

if set to true (default: false), the quadratic term is modelled via a second order cone approach
"""
cb_QPset_use_socqp!(self::CBQPSolverParameters, s::Bool) = @ccall libcb.cb_qpsolverparameters_qpset_use_socqp(self.data::Ptr{Cvoid}, s::Cint)::Cint

@doc raw"""
    cb_QPset_allow_UQPSolver!(self::CBQPSolverParameters, allow::Bool)

set to true/false if switching to the unconstrained solver is allowed or not
"""
cb_QPset_allow_UQPSolver!(self::CBQPSolverParameters, allow::Bool) = @ccall libcb.cb_qpsolverparameters_qpset_allow_uqpsolver(self.data::Ptr{Cvoid}, allow::Cint)::Cint

@doc raw"""
    cb_QPallow_UQPSolver!(self::CBQPSolverParameters)

set to true/false if switching to the unconstrained solver is allowed or not
"""
cb_QPallow_UQPSolver!(self::CBQPSolverParameters) = Bool(@ccall libcb.cb_qpsolverparameters_qpallow_uqpsolver(self.data::Ptr{Cvoid})::Cint)

