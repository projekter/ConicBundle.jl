@doc raw"""
    cb_QPclear!(self::CBQPSolver)

(re)initialize to empty
"""
cb_QPclear!(self::CBQPSolver) = @ccall libcb.cb_qpsolver_qpclear(self.data::Ptr{Cvoid})::Cvoid

@doc raw"""
    CBQPSolver(cbinc::Integer = -1)

default constructor
"""
CBQPSolver(cbinc::Integer = -1) = CBQPSolver(@ccall libcb.cb_qpsolver_new(cbinc::Cint)::Ptr{Cvoid})

@doc raw"""
    cb_QPset_parameters!(self::CBQPSolver, params::Union{<:CBQPSolverParametersObject,Nothing})

check whether the parameters are QPSolverParameters and set them if so
"""
cb_QPset_parameters!(self::CBQPSolver, params::Union{<:CBQPSolverParametersObject,Nothing}) = @ccall libcb.cb_qpsolver_qpset_parameters(self.data::Ptr{Cvoid}, (isnothing(params) ? C_NULL : params.data)::Ptr{Cvoid})::Cint

@doc raw"""
    cb_QPsupports_yfixing!(self::CBQPSolver)

yfixing is currently not supported, this returns false.
"""
cb_QPsupports_yfixing!(self::CBQPSolver) = Bool(@ccall libcb.cb_qpsolver_qpsupports_yfixing(self.data::Ptr{Cvoid})::Cint)

@doc raw"""
    cb_QPsupports_updates!(self::CBQPSolver)

* @brief return true iff the code supports QPupdate(), i.e., it supports external updates of the groundset aggregate in order to model constraints not included explicitly in the QP's model
    
"""
cb_QPsupports_updates!(self::CBQPSolver) = Bool(@ccall libcb.cb_qpsolver_qpsupports_updates(self.data::Ptr{Cvoid})::Cint)

@doc raw"""
    cb_QPget_lower_bound!(self::CBQPSolver)

returns the current lower bound on the optimal value (if feasibility is good enough)
"""
cb_QPget_lower_bound!(self::CBQPSolver) = @ccall libcb.cb_qpsolver_qpget_lower_bound(self.data::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_QPsolve!(self::CBQPSolver, center_y::CBMatrix, lower_bound::Real, upper_bound::Real, relprec::Real, Hp::Union{<:CBQPSolverProxObject,Nothing}, gs_aggr::CBMinorantPointer, yfixed::Union{<:CBIndexmatrix,Nothing})

solve the current bundle subproblem so that precision requirements are met (see \ref InternalQPSolverInterface)
"""
cb_QPsolve!(self::CBQPSolver, center_y::CBMatrix, lower_bound::Real, upper_bound::Real, relprec::Real, Hp::Union{<:CBQPSolverProxObject,Nothing}, gs_aggr::CBMinorantPointer, yfixed::Union{<:CBIndexmatrix,Nothing}) = @ccall libcb.cb_qpsolver_qpsolve(self.data::Ptr{Cvoid}, center_y.data::Ptr{Cvoid}, lower_bound::Cdouble, upper_bound::Cdouble, relprec::Cdouble, (isnothing(Hp) ? C_NULL : Hp.data)::Ptr{Cvoid}, gs_aggr.data::Ptr{Cvoid}, (isnothing(yfixed) ? C_NULL : yfixed.data)::Ptr{Cvoid})::Cint

@doc raw"""
    cb_QPupdate!(self::CBQPSolver, center_y::CBMatrix, lower_bound::Real, upper_bound::Real, relprec::Real, Hp::Union{<:CBQPSolverProxObject,Nothing}, gs_aggr::CBMinorantPointer, yfixed::Union{<:CBIndexmatrix,Nothing}, delta_gs_aggr::CBMinorantPointer, delta_index::CBIndexmatrix)

solve the bundle subproblem for updated box multipliers so that precision requirements are met (see \ref InternalQPSolverInterface). This routine is typically not called for this solver, because box constraints are included explicitly.
"""
cb_QPupdate!(self::CBQPSolver, center_y::CBMatrix, lower_bound::Real, upper_bound::Real, relprec::Real, Hp::Union{<:CBQPSolverProxObject,Nothing}, gs_aggr::CBMinorantPointer, yfixed::Union{<:CBIndexmatrix,Nothing}, delta_gs_aggr::CBMinorantPointer, delta_index::CBIndexmatrix) = @ccall libcb.cb_qpsolver_qpupdate(self.data::Ptr{Cvoid}, center_y.data::Ptr{Cvoid}, lower_bound::Cdouble, upper_bound::Cdouble, relprec::Cdouble, (isnothing(Hp) ? C_NULL : Hp.data)::Ptr{Cvoid}, gs_aggr.data::Ptr{Cvoid}, (isnothing(yfixed) ? C_NULL : yfixed.data)::Ptr{Cvoid}, delta_gs_aggr.data::Ptr{Cvoid}, delta_index.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_QPresolve!(self::CBQPSolver, lower_bound::Real, upper_bound::Real, relprec::Real)

resolve the bundle subproblem (usually because of modified penalty parameters) so that precision requirements are met (see \ref InternalQPSolverInterface)
"""
cb_QPresolve!(self::CBQPSolver, lower_bound::Real, upper_bound::Real, relprec::Real) = @ccall libcb.cb_qpsolver_qpresolve(self.data::Ptr{Cvoid}, lower_bound::Cdouble, upper_bound::Cdouble, relprec::Cdouble)::Cint

@doc raw"""
    cb_QPget_solution!(self::CBQPSolver, new_point::CBMatrix, gsaggr_gradient::CBMatrix)

retrieve the solution produced (see \ref InternalQPSolverInterface)
"""
function cb_QPget_solution!(self::CBQPSolver, new_point::CBMatrix, gsaggr_gradient::CBMatrix)
    gsaggr_offset = Ref{Float64}()
    augval_ub = Ref{Float64}()
    augval_lb = Ref{Float64}()
    @ccall libcb.cb_qpsolver_qpget_solution(self.data::Ptr{Cvoid}, augval_lb::Ref{Float64}, augval_ub::Ref{Float64}, new_point.data::Ptr{Cvoid}, gsaggr_offset::Ref{Float64}, gsaggr_gradient.data::Ptr{Cvoid})::Cint
    return augval_lb[], augval_ub[], gsaggr_offset[]
end

@doc raw"""
    cb_QPprint_statistics!(self::CBQPSolver, param1::Integer = 0)

currently it does nothing
"""
cb_QPprint_statistics!(self::CBQPSolver, param1::Integer = 0) = @ccall libcb.cb_qpsolver_qpprint_statistics(self.data::Ptr{Cvoid}, param1::Cint)::Cvoid

@doc raw"""
    cb_QPstart_modification!(self::CBQPSolver)

return a new modification object on the heap that is initialized for modification of *this
"""
cb_QPstart_modification!(self::CBQPSolver) = CBGroundsetModification(@ccall libcb.cb_qpsolver_qpstart_modification(self.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_QPapply_modification!(self::CBQPSolver, mdf::CBGroundsetModification)

groundset changes are communicated to the solver here
"""
cb_QPapply_modification!(self::CBQPSolver, mdf::CBGroundsetModification) = @ccall libcb.cb_qpsolver_qpapply_modification(self.data::Ptr{Cvoid}, mdf.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_QPis_feasible!(self::CBQPSolver, y::CBMatrix, relprec::Real = 1e-10)

check feasiblity of y for the current groundset constraints
"""
cb_QPis_feasible!(self::CBQPSolver, y::CBMatrix, relprec::Real = 1e-10) = Bool(@ccall libcb.cb_qpsolver_qpis_feasible(self.data::Ptr{Cvoid}, y.data::Ptr{Cvoid}, relprec::Cdouble)::Cint)

@doc raw"""
    cb_QPensure_feasibility!(self::CBQPSolver, y::CBMatrix, ychanged::Bool, inHp::Union{<:CBQPSolverProxObject,Nothing}, relprec::Real = 1e-10)

makes y feasible if not so, see Groundset::ensure_feasibility()
"""
cb_QPensure_feasibility!(self::CBQPSolver, y::CBMatrix, ychanged::Bool, inHp::Union{<:CBQPSolverProxObject,Nothing}, relprec::Real = 1e-10) = @ccall libcb.cb_qpsolver_qpensure_feasibility(self.data::Ptr{Cvoid}, y.data::Ptr{Cvoid}, ychanged::Ref{Cint}, (isnothing(inHp) ? C_NULL : inHp.data)::Ptr{Cvoid}, relprec::Cdouble)::Cint

@doc raw"""
    cb_solve!(self::CBQPSolver, Hp::Union{<:CBBundleProxObject,Nothing}, c::CBMatrix, gamma::Real, lowerbound::Real, upperbound::Real, relprec::Real, skip_factor::Real)

solve the quadratic problem for the given cost function and precision (without cutting model, usually for finding feasible starting points)
"""
cb_solve!(self::CBQPSolver, Hp::Union{<:CBBundleProxObject,Nothing}, c::CBMatrix, gamma::Real, lowerbound::Real, upperbound::Real, relprec::Real, skip_factor::Real) = @ccall libcb.cb_qpsolver_solve(self.data::Ptr{Cvoid}, (isnothing(Hp) ? C_NULL : Hp.data)::Ptr{Cvoid}, c.data::Ptr{Cvoid}, gamma::Cdouble, lowerbound::Cdouble, upperbound::Cdouble, relprec::Cdouble, skip_factor::Cdouble)::Cint

@doc raw"""
    cb_QPprefer_UQPSolver(self::CBQPSolver, param0::Union{<:CBQPSolverProxObject,Nothing})

returns true if, for the current constraints and the requested ProxObject, it might be better to use the internal unconstrained QP solver (which can deal with box constraints by a work-around)
"""
cb_QPprefer_UQPSolver(self::CBQPSolver, param0::Union{<:CBQPSolverProxObject,Nothing}) = Bool(@ccall libcb.cb_qpsolver_qpprefer_uqpsolver(self.data::Ptr{Cvoid}, (isnothing(param0) ? C_NULL : param0.data)::Ptr{Cvoid})::Cint)

@doc raw"""
    cb_QPconstrained(self::CBQPSolver)

returns false if the feasible set is the entire space (unconstrained optimization), true otherwise.
"""
cb_QPconstrained(self::CBQPSolver) = Bool(@ccall libcb.cb_qpsolver_qpconstrained(self.data::Ptr{Cvoid})::Cint)

@doc raw"""
    cb_rowdim(self::CBQPSolver)

number of linear constraints
"""
cb_rowdim(self::CBQPSolver) = @ccall libcb.cb_qpsolver_rowdim(self.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_get_lby(self::CBQPSolver)

returns the lower bounds on y
"""
cb_get_lby(self::CBQPSolver) = (@ccall libcb.cb_qpsolver_get_lby(self.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

@doc raw"""
    cb_get_uby(self::CBQPSolver)

returns the upper bounds on y
"""
cb_get_uby(self::CBQPSolver) = (@ccall libcb.cb_qpsolver_get_uby(self.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

@doc raw"""
    cb_get_lbindex(self::CBQPSolver)

returns the indices of variable lower bounds > ConicBundle::CB_minus_infinity
"""
cb_get_lbindex(self::CBQPSolver) = (@ccall libcb.cb_qpsolver_get_lbindex(self.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

@doc raw"""
    cb_get_ubindex(self::CBQPSolver)

returns the indices of variable lower bounds < ConicBundle::CB_plus_infinity
"""
cb_get_ubindex(self::CBQPSolver) = (@ccall libcb.cb_qpsolver_get_ubindex(self.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

@doc raw"""
    cb_get_A(self::CBQPSolver)

returns the constraint matrix of the feasible set
"""
cb_get_A(self::CBQPSolver) = (@ccall libcb.cb_qpsolver_get_a(self.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

@doc raw"""
    cb_get_rhslb(self::CBQPSolver)

returns the constraint lower bounds
"""
cb_get_rhslb(self::CBQPSolver) = (@ccall libcb.cb_qpsolver_get_rhslb(self.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

@doc raw"""
    cb_get_rhsub(self::CBQPSolver)

returns the constraint upper bounds
"""
cb_get_rhsub(self::CBQPSolver) = (@ccall libcb.cb_qpsolver_get_rhsub(self.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

@doc raw"""
    cb_get_rhslbind(self::CBQPSolver)

returns the indices with constraint lower bound slacks
"""
cb_get_rhslbind(self::CBQPSolver) = (@ccall libcb.cb_qpsolver_get_rhslbind(self.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

@doc raw"""
    cb_get_rhsubind(self::CBQPSolver)

returns the indices with constraint upper bound slacks
"""
cb_get_rhsubind(self::CBQPSolver) = (@ccall libcb.cb_qpsolver_get_rhsubind(self.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

@doc raw"""
    cb_mfile_data(self::CBQPSolver)

output the data describing the QP in m-file style
"""
cb_mfile_data(self::CBQPSolver) = @ccall libcb.cb_qpsolver_mfile_data(self.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_set_cbout!(self::CBQPSolver, incr::Integer = -1)

set output settings
"""
cb_set_cbout!(self::CBQPSolver, incr::Integer = -1) = @ccall libcb.cb_qpsolver_set_cbout(self.data::Ptr{Cvoid}, incr::Cint)::Cvoid

