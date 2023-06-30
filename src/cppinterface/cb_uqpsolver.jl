@doc raw"""
    cb_clear!(self::CBUQPSolver)

reset to "empty/no" model
"""
cb_clear!(self::CBUQPSolver) = @ccall libcb.cb_uqpsolver_clear(self.data::Ptr{Cvoid})::Cvoid

@doc raw"""
    cb_set_defaults!(self::CBUQPSolver)

reset parameters to default values
"""
cb_set_defaults!(self::CBUQPSolver) = @ccall libcb.cb_uqpsolver_set_defaults(self.data::Ptr{Cvoid})::Cvoid

@doc raw"""
    CBUQPSolver(cbinc::Integer = -1)

default constructor
"""
CBUQPSolver(cbinc::Integer = -1) = CBUQPSolver(@ccall libcb.cb_uqpsolver_new(cbinc::Cint)::Ptr{Cvoid})

@doc raw"""
    cb_init_size!(self::CBUQPSolver, maxdim::Integer)

reserve memory for this size
"""
cb_init_size!(self::CBUQPSolver, maxdim::Integer) = @ccall libcb.cb_uqpsolver_init_size(self.data::Ptr{Cvoid}, maxdim::Cint)::Cvoid

@doc raw"""
    cb_get_Q(self::CBUQPSolver)

returns the quadratic cost matrix
"""
cb_get_Q(self::CBUQPSolver) = (@ccall libcb.cb_uqpsolver_get_q(self.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

@doc raw"""
    cb_get_c(self::CBUQPSolver)

returns the linear cost vector
"""
cb_get_c(self::CBUQPSolver) = (@ccall libcb.cb_uqpsolver_get_c(self.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

@doc raw"""
    cb_get_offset(self::CBUQPSolver)

returns the constant cost offset value
"""
cb_get_offset(self::CBUQPSolver) = @ccall libcb.cb_uqpsolver_get_offset(self.data::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_set_termbounds!(self::CBUQPSolver, lb::Real, ub::Real)

sets the termination lower and upper bounds
"""
cb_set_termbounds!(self::CBUQPSolver, lb::Real, ub::Real) = @ccall libcb.cb_uqpsolver_set_termbounds(self.data::Ptr{Cvoid}, lb::Cdouble, ub::Cdouble)::Cvoid

@doc raw"""
    cb_set_termeps!(self::CBUQPSolver, te::Real)

sets the termination precision
"""
cb_set_termeps!(self::CBUQPSolver, te::Real) = @ccall libcb.cb_uqpsolver_set_termeps(self.data::Ptr{Cvoid}, te::Cdouble)::Cvoid

@doc raw"""
    cb_set_maxiter!(self::CBUQPSolver, mi::Integer)

sets the upper bound on the number of interior point iterations (<0 means no bound)
"""
cb_set_maxiter!(self::CBUQPSolver, mi::Integer) = @ccall libcb.cb_uqpsolver_set_maxiter(self.data::Ptr{Cvoid}, mi::Cint)::Cvoid

@doc raw"""
    cb_get_iter(self::CBUQPSolver)

return the number of iterations of the last solve
"""
cb_get_iter(self::CBUQPSolver) = @ccall libcb.cb_uqpsolver_get_iter(self.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_get_status(self::CBUQPSolver)

return the status of the last solve
"""
cb_get_status(self::CBUQPSolver) = @ccall libcb.cb_uqpsolver_get_status(self.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_get_termeps(self::CBUQPSolver)

return the termination precision
"""
cb_get_termeps(self::CBUQPSolver) = @ccall libcb.cb_uqpsolver_get_termeps(self.data::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_get_maxiter(self::CBUQPSolver)

return the upper bound on interior point interations
"""
cb_get_maxiter(self::CBUQPSolver) = @ccall libcb.cb_uqpsolver_get_maxiter(self.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_get_primalval(self::CBUQPSolver)

return the primal objective value (lower bound) of the last solve
"""
cb_get_primalval(self::CBUQPSolver) = @ccall libcb.cb_uqpsolver_get_primalval(self.data::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_get_dualval(self::CBUQPSolver)

return the dual objective value (upper bound) of the last solve
"""
cb_get_dualval(self::CBUQPSolver) = @ccall libcb.cb_uqpsolver_get_dualval(self.data::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_get_x!(self::CBUQPSolver)

return the joint model vector (primal solution) produced by the last solve
"""
cb_get_x!(self::CBUQPSolver) = (@ccall libcb.cb_uqpsolver_get_x(self.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

@doc raw"""
    cb_get_y!(self::CBUQPSolver)

return the joint dual vector (dual solution) produced by the last solve
"""
cb_get_y!(self::CBUQPSolver) = (@ccall libcb.cb_uqpsolver_get_y(self.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

@doc raw"""
    cb_solve!(self::CBUQPSolver, Q::CBSymmatrix, c::CBMatrix, offset::Real)

solve the QP for this cost function from scratch
"""
cb_solve!(self::CBUQPSolver, Q::CBSymmatrix, c::CBMatrix, offset::Real) = @ccall libcb.cb_uqpsolver_solve(self.data::Ptr{Cvoid}, Q.data::Ptr{Cvoid}, c.data::Ptr{Cvoid}, offset::Cdouble)::Cint

@doc raw"""
    cb_resolve!(self::CBUQPSolver)

resolve the QP for the same cost function as last time with slightly modified feasible set
"""
cb_resolve!(self::CBUQPSolver) = @ccall libcb.cb_uqpsolver_resolve(self.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_update!(self::CBUQPSolver, dQ::CBSymmatrix, dc::CBMatrix, doffset::Real)

resolve the QP for the same feasible set but update the cost terms by the given argumetns first
"""
cb_update!(self::CBUQPSolver, dQ::CBSymmatrix, dc::CBMatrix, doffset::Real) = @ccall libcb.cb_uqpsolver_update(self.data::Ptr{Cvoid}, dQ.data::Ptr{Cvoid}, dc.data::Ptr{Cvoid}, doffset::Cdouble)::Cint

@doc raw"""
    cb_print_statistics(self::CBUQPSolver)

output some statistical information on performance
"""
cb_print_statistics(self::CBUQPSolver) = @ccall libcb.cb_uqpsolver_print_statistics(self.data::Ptr{Cvoid})::Cvoid

@doc raw"""
    cb_save(self::CBUQPSolver)

save the current settings and values
"""
cb_save(self::CBUQPSolver) = @ccall libcb.cb_uqpsolver_save(self.data::Ptr{Cvoid})::Cvoid

@doc raw"""
    cb_QPclear!(self::CBUQPSolver)

calls clear(), i.e. reinitialize completely
"""
cb_QPclear!(self::CBUQPSolver) = @ccall libcb.cb_uqpsolver_qpclear(self.data::Ptr{Cvoid})::Cvoid

@doc raw"""
    cb_QPset_parameters!(self::CBUQPSolver, param0::Union{<:CBQPSolverParametersObject,Nothing})

does nothing here
"""
cb_QPset_parameters!(self::CBUQPSolver, param0::Union{<:CBQPSolverParametersObject,Nothing}) = @ccall libcb.cb_uqpsolver_qpset_parameters(self.data::Ptr{Cvoid}, (isnothing(param0) ? C_NULL : param0.data)::Ptr{Cvoid})::Cint

@doc raw"""
    cb_apply_modification!(self::CBUQPSolver, param0::CBGroundsetModification)

does nothing here (unconstrained case)
"""
cb_apply_modification!(self::CBUQPSolver, param0::CBGroundsetModification) = @ccall libcb.cb_uqpsolver_apply_modification(self.data::Ptr{Cvoid}, param0.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_QPsupports_yfixing!(self::CBUQPSolver)

no difficulty if the BundleProxObject does it
"""
cb_QPsupports_yfixing!(self::CBUQPSolver) = Bool(@ccall libcb.cb_uqpsolver_qpsupports_yfixing(self.data::Ptr{Cvoid})::Cint)

@doc raw"""
    cb_QPsupports_updates!(self::CBUQPSolver)

* @brief return true iff the code supports QPupdate(), i.e., it supports external updates of the groundset aggregate in order to model constraints not included explicitly in the QP's model
    
"""
cb_QPsupports_updates!(self::CBUQPSolver) = Bool(@ccall libcb.cb_uqpsolver_qpsupports_updates(self.data::Ptr{Cvoid})::Cint)

@doc raw"""
    cb_QPget_lower_bound!(self::CBUQPSolver)

return the lower bound on the objective value of the bundle subproblem
"""
cb_QPget_lower_bound!(self::CBUQPSolver) = @ccall libcb.cb_uqpsolver_qpget_lower_bound(self.data::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_QPsolve!(self::CBUQPSolver, center_y::CBMatrix, lower_bound::Real, upper_bound::Real, relprec::Real, Hp::Union{<:CBQPSolverProxObject,Nothing}, gs_aggr::CBMinorantPointer, yfixed::Union{<:CBIndexmatrix,Nothing})

see QPSolverObject::QPsolve() and  \ref UnconstrainedQPSolver
"""
cb_QPsolve!(self::CBUQPSolver, center_y::CBMatrix, lower_bound::Real, upper_bound::Real, relprec::Real, Hp::Union{<:CBQPSolverProxObject,Nothing}, gs_aggr::CBMinorantPointer, yfixed::Union{<:CBIndexmatrix,Nothing}) = @ccall libcb.cb_uqpsolver_qpsolve(self.data::Ptr{Cvoid}, center_y.data::Ptr{Cvoid}, lower_bound::Cdouble, upper_bound::Cdouble, relprec::Cdouble, (isnothing(Hp) ? C_NULL : Hp.data)::Ptr{Cvoid}, gs_aggr.data::Ptr{Cvoid}, (isnothing(yfixed) ? C_NULL : yfixed.data)::Ptr{Cvoid})::Cint

@doc raw"""
    cb_QPupdate!(self::CBUQPSolver, center_y::CBMatrix, lower_bound::Real, upper_bound::Real, relprec::Real, Hp::Union{<:CBQPSolverProxObject,Nothing}, gs_aggr::CBMinorantPointer, yfixed::Union{<:CBIndexmatrix,Nothing}, delta_gs_subg::CBMinorantPointer, delta_index::CBIndexmatrix)

see QPSolverObject::QPupdate() and  \ref UnconstrainedQPSolver
"""
cb_QPupdate!(self::CBUQPSolver, center_y::CBMatrix, lower_bound::Real, upper_bound::Real, relprec::Real, Hp::Union{<:CBQPSolverProxObject,Nothing}, gs_aggr::CBMinorantPointer, yfixed::Union{<:CBIndexmatrix,Nothing}, delta_gs_subg::CBMinorantPointer, delta_index::CBIndexmatrix) = @ccall libcb.cb_uqpsolver_qpupdate(self.data::Ptr{Cvoid}, center_y.data::Ptr{Cvoid}, lower_bound::Cdouble, upper_bound::Cdouble, relprec::Cdouble, (isnothing(Hp) ? C_NULL : Hp.data)::Ptr{Cvoid}, gs_aggr.data::Ptr{Cvoid}, (isnothing(yfixed) ? C_NULL : yfixed.data)::Ptr{Cvoid}, delta_gs_subg.data::Ptr{Cvoid}, delta_index.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_QPresolve!(self::CBUQPSolver, lower_bound::Real, upper_bound::Real, relprec::Real)

see QPSolverObject::QPresolve() and  \ref UnconstrainedQPSolver
"""
cb_QPresolve!(self::CBUQPSolver, lower_bound::Real, upper_bound::Real, relprec::Real) = @ccall libcb.cb_uqpsolver_qpresolve(self.data::Ptr{Cvoid}, lower_bound::Cdouble, upper_bound::Cdouble, relprec::Cdouble)::Cint

@doc raw"""
    cb_QPget_solution!(self::CBUQPSolver, new_point::CBMatrix, gsaggr_gradient::CBMatrix)

the unconstrained solver can only provide this information if the references to the input data of QPsolve() and QPudate() are still available and unchanged, otherwise the behavior is undefined and will hopefully return 1 if not valid
"""
function cb_QPget_solution!(self::CBUQPSolver, new_point::CBMatrix, gsaggr_gradient::CBMatrix)
    gsaggr_offset = Ref{Float64}()
    augval_ub = Ref{Float64}()
    augval_lb = Ref{Float64}()
    @ccall libcb.cb_uqpsolver_qpget_solution(self.data::Ptr{Cvoid}, augval_lb::Ref{Float64}, augval_ub::Ref{Float64}, new_point.data::Ptr{Cvoid}, gsaggr_offset::Ref{Float64}, gsaggr_gradient.data::Ptr{Cvoid})::Cint
    return augval_lb[], augval_ub[], gsaggr_offset[]
end

@doc raw"""
    cb_QPprefer_UQPSolver(self::CBUQPSolver, param0::Union{<:CBQPSolverProxObject,Nothing})

as this IS the internal UQPSolver, the answer doesn't matter, but the solver would certainly say yes
"""
cb_QPprefer_UQPSolver(self::CBUQPSolver, param0::Union{<:CBQPSolverProxObject,Nothing}) = Bool(@ccall libcb.cb_uqpsolver_qpprefer_uqpsolver(self.data::Ptr{Cvoid}, (isnothing(param0) ? C_NULL : param0.data)::Ptr{Cvoid})::Cint)

@doc raw"""
    cb_QPconstrained(self::CBUQPSolver)

it is always unconstrained
"""
cb_QPconstrained(self::CBUQPSolver) = Bool(@ccall libcb.cb_uqpsolver_qpconstrained(self.data::Ptr{Cvoid})::Cint)

@doc raw"""
    cb_QPstart_modification!(self::CBUQPSolver)

returns 0 because no modifications are applicable
"""
cb_QPstart_modification!(self::CBUQPSolver) = CBGroundsetModification(@ccall libcb.cb_uqpsolver_qpstart_modification(self.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_QPapply_modification!(self::CBUQPSolver, param0::CBGroundsetModification)

no modifications need to be carried out here as there is no data to be modified, so any modification succeeds
"""
cb_QPapply_modification!(self::CBUQPSolver, param0::CBGroundsetModification) = @ccall libcb.cb_uqpsolver_qpapply_modification(self.data::Ptr{Cvoid}, param0.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_QPis_feasible!(self::CBUQPSolver, param0::CBMatrix, param1::Real = 1e-10)

for the unconstrained solver any point is feasible, because it cannot even check the dimension of the design space
"""
cb_QPis_feasible!(self::CBUQPSolver, param0::CBMatrix, param1::Real = 1e-10) = Bool(@ccall libcb.cb_uqpsolver_qpis_feasible(self.data::Ptr{Cvoid}, param0.data::Ptr{Cvoid}, param1::Cdouble)::Cint)

@doc raw"""
    cb_QPensure_feasibility!(self::CBUQPSolver, param0::CBMatrix, ychanged::Bool, param2::Union{<:CBQPSolverProxObject,Nothing}, param3::Real = 1e-10)

for the unconstrained solver any point is feasible, because it cannot even check the dimension of the design space, so any y is left unchanged
"""
cb_QPensure_feasibility!(self::CBUQPSolver, param0::CBMatrix, ychanged::Bool, param2::Union{<:CBQPSolverProxObject,Nothing}, param3::Real = 1e-10) = @ccall libcb.cb_uqpsolver_qpensure_feasibility(self.data::Ptr{Cvoid}, param0.data::Ptr{Cvoid}, ychanged::Ref{Cint}, (isnothing(param2) ? C_NULL : param2.data)::Ptr{Cvoid}, param3::Cdouble)::Cint

@doc raw"""
    cb_QPprint_statistics!(self::CBUQPSolver, param1::Integer = 0)

outputs some statistical data about solver performance
"""
cb_QPprint_statistics!(self::CBUQPSolver, param1::Integer = 0) = @ccall libcb.cb_uqpsolver_qpprint_statistics(self.data::Ptr{Cvoid}, param1::Cint)::Cvoid

