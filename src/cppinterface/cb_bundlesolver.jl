@doc raw"""
    CBBundleSolver(dim::Integer, bp::Union{<:CBBundleModel,Nothing} = nothing, incr::Integer = -1)

calls initialize(CH_Matrix_Classes::Integer, BundleModel*)
"""
CBBundleSolver(dim::Integer, bp::Union{<:CBBundleModel,Nothing} = nothing, incr::Integer = -1) = CBBundleSolver(@ccall libcb.cb_bundlesolver_new2(dim::Cint, (isnothing(bp) ? C_NULL : bp.data)::Ptr{Cvoid}, incr::Cint)::Ptr{Cvoid})

@doc raw"""
    CBBundleSolver(gs::Union{<:CBGroundset,Nothing}, bp::Union{<:CBBundleModel,Nothing} = nothing, incr::Integer = -1)

calls initialize(Groundset*,BundleModel*)
"""
CBBundleSolver(gs::Union{<:CBGroundset,Nothing}, bp::Union{<:CBBundleModel,Nothing} = nothing, incr::Integer = -1) = CBBundleSolver(@ccall libcb.cb_bundlesolver_new3((isnothing(gs) ? C_NULL : gs.data)::Ptr{Cvoid}, (isnothing(bp) ? C_NULL : bp.data)::Ptr{Cvoid}, incr::Cint)::Ptr{Cvoid})

@doc raw"""
    cb_set_defaults!(self::CBBundleSolver)

resets all parameters to default values and calls BundleTerminator::set_defaults() for *terminator and BundleWeight::set_defaults() for *bundleweight
"""
cb_set_defaults!(self::CBBundleSolver) = @ccall libcb.cb_bundlesolver_set_defaults(self.data::Ptr{Cvoid})::Cvoid

@doc raw"""
    cb_clear!(self::CBBundleSolver)

resets all variables and pointers to classes to initial state and calls set_defaults()
"""
cb_clear!(self::CBBundleSolver) = @ccall libcb.cb_bundlesolver_clear(self.data::Ptr{Cvoid})::Cvoid

@doc raw"""
    cb_clear_fails!(self::CBBundleSolver)

resets all fail counts to zero (call this to resume computations
"""
cb_clear_fails!(self::CBBundleSolver) = @ccall libcb.cb_bundlesolver_clear_fails(self.data::Ptr{Cvoid})::Cvoid

@doc raw"""
    cb_initialize!(self::CBBundleSolver, dim::Integer, bp::Union{<:CBBundleModel,Nothing} = nothing)

calls clear(), initializes an unconstrained groundset to this dimension and sets the bundle model to bp
"""
cb_initialize!(self::CBBundleSolver, dim::Integer, bp::Union{<:CBBundleModel,Nothing} = nothing) = @ccall libcb.cb_bundlesolver_initialize(self.data::Ptr{Cvoid}, dim::Cint, (isnothing(bp) ? C_NULL : bp.data)::Ptr{Cvoid})::Cint

@doc raw"""
    cb_initialize!(self::CBBundleSolver, gs::Union{<:CBGroundset,Nothing}, bp::Union{<:CBBundleModel,Nothing} = nothing)

calls clear(), initializes the groundset to gs and the bundle model to bp
"""
cb_initialize!(self::CBBundleSolver, gs::Union{<:CBGroundset,Nothing}, bp::Union{<:CBBundleModel,Nothing} = nothing) = @ccall libcb.cb_bundlesolver_initialize2(self.data::Ptr{Cvoid}, (isnothing(gs) ? C_NULL : gs.data)::Ptr{Cvoid}, (isnothing(bp) ? C_NULL : bp.data)::Ptr{Cvoid})::Cint

@doc raw"""
    cb_set_out!(self::CBBundleSolver, pril::Integer = 1)

set output and outputlevel of warnings and errors recursively, see CBout
"""
cb_set_out!(self::CBBundleSolver, pril::Integer = 1) = @ccall libcb.cb_bundlesolver_set_out(self.data::Ptr{Cvoid}, pril::Cint)::Cvoid

@doc raw"""
    cb_set_cbout!(self::CBBundleSolver, incr::Integer)

set output and outputlevel of warnings and errors recursively with CBout
"""
cb_set_cbout!(self::CBBundleSolver, incr::Integer) = @ccall libcb.cb_bundlesolver_set_cbout(self.data::Ptr{Cvoid}, incr::Cint)::Cvoid

@doc raw"""
    cb_set_model!(self::CBBundleSolver, bp::Union{<:CBBundleModel,Nothing})

set/change the model that should be optimized over (for the existing groundset and starting point)
"""
cb_set_model!(self::CBBundleSolver, bp::Union{<:CBBundleModel,Nothing}) = @ccall libcb.cb_bundlesolver_set_model(self.data::Ptr{Cvoid}, (isnothing(bp) ? C_NULL : bp.data)::Ptr{Cvoid})::Cint

@doc raw"""
    cb_set_new_center!(self::CBBundleSolver, yp::Union{<:CBMatrix,Nothing} = nothing)

replace the current center by *yp or, if yp==0, by the default starting point
"""
cb_set_new_center!(self::CBBundleSolver, yp::Union{<:CBMatrix,Nothing} = nothing) = @ccall libcb.cb_bundlesolver_set_new_center(self.data::Ptr{Cvoid}, (isnothing(yp) ? C_NULL : yp.data)::Ptr{Cvoid})::Cint

@doc raw"""
    cb_solve!(self::CBBundleSolver, maxsteps::Integer = 0, stop_at_descent_steps::Bool = false)

execute at most maxsteps iterations of the bundle method stopping before if termination occurs or stop_at_descent_steps==true and a descent step occurs; maxsteps<=0 indicates no bound on the steps
"""
cb_solve!(self::CBBundleSolver, maxsteps::Integer = 0, stop_at_descent_steps::Bool = false) = @ccall libcb.cb_bundlesolver_solve(self.data::Ptr{Cvoid}, maxsteps::Cint, stop_at_descent_steps::Cint)::Cint

@doc raw"""
    cb_print_line_summary(self::CBBundleSolver)

* @brief print a one line summary about the current state of progress of the algorithm to @a out

        Assuming that an external clock has been set, the line looks as follows

        hh:mm:ss.hh endit dd ii uu weightu aggr_dnorm modelval center_ub (N)

        These are
        - hh:mm:ss.hh hours, mintues, seconds, hundreth of secondes
        - "endit" this string is short for "end of iteration" and is
          convenient for the unix command grep; in the case the code
          has terminated it shows "_endit" instead
        - dd gives the number of descent steps as in @a descent_steps
        - ii gives the total number of iterations of the bundle method,
          counting descent steps and null steps
        - uu gives the total number of quadratic bundle subproblem
          evaluations as counted by @a sumupdatecnt
        - weight gives the value of the weight factor for the proximal
          term used in the last quadratic bundle subproblem
        - aggr_dnorm gives the dual norm of the aggregate (dual with
          respect to the quadratic proximal term)
        - modelval gives the model value used for deciding on null or
          descent step; this is @a linval if @a use_linval ==true and
    @a cutval otherwise
        - center_ub gives the upper bound computed by the oracle on the
           function value in the center
        - "N" is only shown if the code has not terminated and a null
           step just occured before returning (null_step is true)
     
"""
cb_print_line_summary(self::CBBundleSolver) = @ccall libcb.cb_bundlesolver_print_line_summary(self.data::Ptr{Cvoid})::Cvoid

@doc raw"""
    cb_apply_modification!(self::CBBundleSolver, gsmdf::CBGroundsetModification)

modify the groundset as described by GroundsetModification and inform the oracle function(s) about this change (calls the other apply_modification for an empty FunObjModMap)
"""
cb_apply_modification!(self::CBBundleSolver, gsmdf::CBGroundsetModification) = @ccall libcb.cb_bundlesolver_apply_modification2(self.data::Ptr{Cvoid}, gsmdf.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_set_terminator!(self::CBBundleSolver, bt::Union{<:CBBundleTerminator,Nothing})

replace the previous BundleTerminator by bt; bt will be deleted when replaced or on destruction of this
"""
cb_set_terminator!(self::CBBundleSolver, bt::Union{<:CBBundleTerminator,Nothing}) = @ccall libcb.cb_bundlesolver_set_terminator(self.data::Ptr{Cvoid}, (isnothing(bt) ? C_NULL : bt.data)::Ptr{Cvoid})::Cvoid

@doc raw"""
    cb_set_bundleweight!(self::CBBundleSolver, bw::Union{<:CBBundleWeight,Nothing})

replace the previous BundleWeight by bw; bw will be deleted when replaced or on destruction of this
"""
cb_set_bundleweight!(self::CBBundleSolver, bw::Union{<:CBBundleWeight,Nothing}) = @ccall libcb.cb_bundlesolver_set_bundleweight(self.data::Ptr{Cvoid}, (isnothing(bw) ? C_NULL : bw.data)::Ptr{Cvoid})::Cvoid

@doc raw"""
    cb_set_modeleps!(self::CBBundleSolver, in_eps::Real)

set the required model precision by @a in_eps (if it is positive)
"""
cb_set_modeleps!(self::CBBundleSolver, in_eps::Real) = @ccall libcb.cb_bundlesolver_set_modeleps(self.data::Ptr{Cvoid}, in_eps::Cdouble)::Cvoid

@doc raw"""
    cb_set_mL!(self::CBBundleSolver, in_mL::Real)

set the acceptance level for descent steps (rather don't change this!)
"""
cb_set_mL!(self::CBBundleSolver, in_mL::Real) = @ccall libcb.cb_bundlesolver_set_ml(self.data::Ptr{Cvoid}, in_mL::Cdouble)::Cvoid

@doc raw"""
    cb_set_mN!(self::CBBundleSolver, in_mN::Real)

set the acceptance level for null steps (mL<=in_mN<1., rather don't change this!)
"""
cb_set_mN!(self::CBBundleSolver, in_mN::Real) = @ccall libcb.cb_bundlesolver_set_mn(self.data::Ptr{Cvoid}, in_mN::Cdouble)::Cvoid

@doc raw"""
    cb_set_use_linval!(self::CBBundleSolver, ul::Bool)

if set to true, the value of the aggregate in the candidate is used for deciding on null or descent step, otherwise the model value
"""
cb_set_use_linval!(self::CBBundleSolver, ul::Bool) = @ccall libcb.cb_bundlesolver_set_use_linval(self.data::Ptr{Cvoid}, ul::Cint)::Cvoid

@doc raw"""
    cb_set_do_yfixing!(self::CBBundleSolver, dofix::Bool)

if set to true, the groundset may use a heuristic to decide whether a variable is fixed to one of its bounds (often helps to reduce inner update iterations)
"""
cb_set_do_yfixing!(self::CBBundleSolver, dofix::Bool) = @ccall libcb.cb_bundlesolver_set_do_yfixing(self.data::Ptr{Cvoid}, dofix::Cint)::Cvoid

@doc raw"""
    cb_set_variable_metric!(self::CBBundleSolver, ds::Integer)

0 ... use no scaling, 1 ... use a scaling heuristic, 2 ... also allow groundset to influence the scaling so as to favor feasibility,
"""
cb_set_variable_metric!(self::CBBundleSolver, ds::Integer) = @ccall libcb.cb_bundlesolver_set_variable_metric(self.data::Ptr{Cvoid}, ds::Cint)::Cint

@doc raw"""
    cb_set_prox_diagonal!(self::CBBundleSolver, insc::CBMatrix)

set the prox term to the given diagonal matrix
"""
cb_set_prox_diagonal!(self::CBBundleSolver, insc::CBMatrix) = @ccall libcb.cb_bundlesolver_set_prox_diagonal(self.data::Ptr{Cvoid}, insc.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_set_prox!(self::CBBundleSolver, bsp::Union{<:CBBundleProxObject,Nothing})

set the proximal term by replacing BundleProxObject with bsp, the latter is deleted when replaced or on destruction of this
"""
cb_set_prox!(self::CBBundleSolver, bsp::Union{<:CBBundleProxObject,Nothing}) = @ccall libcb.cb_bundlesolver_set_prox(self.data::Ptr{Cvoid}, (isnothing(bsp) ? C_NULL : bsp.data)::Ptr{Cvoid})::Cint

@doc raw"""
    cb_set_clock!(self::CBBundleSolver, myclock::CBClock)

set the external clock to be used for output
"""
cb_set_clock!(self::CBBundleSolver, myclock::CBClock) = @ccall libcb.cb_bundlesolver_set_clock(self.data::Ptr{Cvoid}, myclock.data::Ptr{Cvoid})::Cvoid

@doc raw"""
    cb_set_max_updates!(self::CBBundleSolver, mu::Integer)

set the maximum number of Gauss-Seidel iterations until the next evaluations for descent/null step, use negative numbers for infinite, 0 or 1 for at most 1
"""
cb_set_max_updates!(self::CBBundleSolver, mu::Integer) = @ccall libcb.cb_bundlesolver_set_max_updates(self.data::Ptr{Cvoid}, mu::Cint)::Cvoid

@doc raw"""
    cb_get_terminate(self::CBBundleSolver)

returns the value of the last call to BundleTerminator::check_termination()
"""
cb_get_terminate(self::CBBundleSolver) = @ccall libcb.cb_bundlesolver_get_terminate(self.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_get_center_objval(self::CBBundleSolver)

returns the upper bound on the objective value in @a center_y
"""
cb_get_center_objval(self::CBBundleSolver) = @ccall libcb.cb_bundlesolver_get_center_objval(self.data::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_get_center_ub(self::CBBundleSolver)

returns the upper bound on the objective in @a center_y returned by the oracle
"""
cb_get_center_ub(self::CBBundleSolver) = @ccall libcb.cb_bundlesolver_get_center_ub(self.data::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_get_center_gs_val(self::CBBundleSolver)

returns the groundset objective in @a center_y r
"""
cb_get_center_gs_val(self::CBBundleSolver) = @ccall libcb.cb_bundlesolver_get_center_gs_val(self.data::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_get_center_y(self::CBBundleSolver)

returns the current center of stability @a center_y (after a descent step this is the same as the candidate)
"""
cb_get_center_y(self::CBBundleSolver) = (@ccall libcb.cb_bundlesolver_get_center_y(self.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

@doc raw"""
    cb_get_cand_objval(self::CBBundleSolver)

returns the upper bound on the objective in @a cand_y
"""
cb_get_cand_objval(self::CBBundleSolver) = @ccall libcb.cb_bundlesolver_get_cand_objval(self.data::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_get_cand_ub(self::CBBundleSolver)

returns the upper bound on the objective in @a cand_y returned by the oracle
"""
cb_get_cand_ub(self::CBBundleSolver) = @ccall libcb.cb_bundlesolver_get_cand_ub(self.data::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_get_cand_gs_val(self::CBBundleSolver)

returns the groundset objective value in @a cand_y
"""
cb_get_cand_gs_val(self::CBBundleSolver) = @ccall libcb.cb_bundlesolver_get_cand_gs_val(self.data::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_get_cand_y(self::CBBundleSolver)

returns the most recent candidate @a cand_y (after a descent step this is the same as the candidate)
"""
cb_get_cand_y(self::CBBundleSolver) = (@ccall libcb.cb_bundlesolver_get_cand_y(self.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

@doc raw"""
    cb_get_aggr_dnormsqr(self::CBBundleSolver)

returns the dual norm squared of the current aggregate (dual w.r.t. the quadratic proximal term)
"""
cb_get_aggr_dnormsqr(self::CBBundleSolver) = @ccall libcb.cb_bundlesolver_get_aggr_dnormsqr(self.data::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_get_aggregate_offset(self::CBBundleSolver)

returns the offset of the current aggregate linear minorant (should be called before any modifications, otherwise this may no longer be correct)
"""
cb_get_aggregate_offset(self::CBBundleSolver) = @ccall libcb.cb_bundlesolver_get_aggregate_offset(self.data::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_get_aggregate(self::CBBundleSolver, aggregate::CBMatrix)

returns the linear term of the current aggregate linear minorant (this should be called before any modifications, otherwise this may no longer be correct or may even cause an error)
"""
cb_get_aggregate(self::CBBundleSolver, aggregate::CBMatrix) = @ccall libcb.cb_bundlesolver_get_aggregate(self.data::Ptr{Cvoid}, aggregate.data::Ptr{Cvoid})::Cvoid

@doc raw"""
    cb_get_gs_aggregate(self::CBBundleSolver)

returns the linear term of the current groundset aggregate linear minorant
"""
cb_get_gs_aggregate(self::CBBundleSolver) = (@ccall libcb.cb_bundlesolver_get_gs_aggregate(self.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

@doc raw"""
    cb_get_model_aggregate(self::CBBundleSolver)

returns the linear term of the latest model aggregate linear minorant
"""
cb_get_model_aggregate(self::CBBundleSolver) = (@ccall libcb.cb_bundlesolver_get_model_aggregate(self.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

@doc raw"""
    cb_get_modelval(self::CBBundleSolver)

returns the model value in the candidate that was used for deciding on null/descent step
"""
cb_get_modelval(self::CBBundleSolver) = @ccall libcb.cb_bundlesolver_get_modelval(self.data::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_get_weight(self::CBBundleSolver)

returns the weight for the proximal term used in the last quadratic subproblem
"""
cb_get_weight(self::CBBundleSolver) = @ccall libcb.cb_bundlesolver_get_weight(self.data::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_get_term_corr(self::CBBundleSolver)

returns the correction factor used in the termination criterion to compensate the strength of the proximal term
"""
cb_get_term_corr(self::CBBundleSolver) = @ccall libcb.cb_bundlesolver_get_term_corr(self.data::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_get_descent_step(self::CBBundleSolver)

returns true if the latest iteration resulted in a descent step (note, get_descent_step() and get_null_step() may both return false e.g. if termination occurs)
"""
cb_get_descent_step(self::CBBundleSolver) = Bool(@ccall libcb.cb_bundlesolver_get_descent_step(self.data::Ptr{Cvoid})::Cint)

@doc raw"""
    cb_get_null_step(self::CBBundleSolver)

returns true if the latest iteration resulted in a null step (note, get_descent_step() and get_null_step() may both return false e.g. if termination occurs)
"""
cb_get_null_step(self::CBBundleSolver) = Bool(@ccall libcb.cb_bundlesolver_get_null_step(self.data::Ptr{Cvoid})::Cint)

@doc raw"""
    cb_get_yfixed(self::CBBundleSolver)

if the groundset has constraints and set_do_yfixing was set with true enty i of the returned matrix is !=0 if the coordinate was fixed at one of its bounds and 0 if the coordinate is still free
"""
cb_get_yfixed(self::CBBundleSolver) = CBIndexmatrix(@ccall libcb.cb_bundlesolver_get_yfixed(self.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_get_modeleps(self::CBBundleSolver)

returns the model precision @a modeleps
"""
cb_get_modeleps(self::CBBundleSolver) = @ccall libcb.cb_bundlesolver_get_modeleps(self.data::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_get_use_linval(self::CBBundleSolver)

returns true if the aggregate linear minorant is used for the model value
"""
cb_get_use_linval(self::CBBundleSolver) = @ccall libcb.cb_bundlesolver_get_use_linval(self.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_get_do_variable_metric(self::CBBundleSolver)

returns true if the proximal term is not of the type BundleIdProx, i.e. if it is not simply the squared Euclidean norm
"""
cb_get_do_variable_metric(self::CBBundleSolver) = @ccall libcb.cb_bundlesolver_get_do_variable_metric(self.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_get_use_variable_metric(self::CBBundleSolver)

returns the value of the variable @a do_dynamic_scaling, see there
"""
cb_get_use_variable_metric(self::CBBundleSolver) = Bool(@ccall libcb.cb_bundlesolver_get_use_variable_metric(self.data::Ptr{Cvoid})::Cint)

@doc raw"""
    cb_get_cntobjeval(self::CBBundleSolver)

returns the number of calls to the oracle since the last clear()
"""
cb_get_cntobjeval(self::CBBundleSolver) = @ccall libcb.cb_bundlesolver_get_cntobjeval(self.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_get_descent_steps(self::CBBundleSolver)

returns the number of descent steps since the last clear()
"""
cb_get_descent_steps(self::CBBundleSolver) = @ccall libcb.cb_bundlesolver_get_descent_steps(self.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_get_innerit(self::CBBundleSolver)

returns the number of bundle method iterations since the last descent step
"""
cb_get_innerit(self::CBBundleSolver) = @ccall libcb.cb_bundlesolver_get_innerit(self.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_get_suminnerit(self::CBBundleSolver)

returns the number of bundle method itrations since the last clear()
"""
cb_get_suminnerit(self::CBBundleSolver) = @ccall libcb.cb_bundlesolver_get_suminnerit(self.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_get_sumupdatecnt(self::CBBundleSolver)

returns the number of model qp subproblems since the last clear()
"""
cb_get_sumupdatecnt(self::CBBundleSolver) = @ccall libcb.cb_bundlesolver_get_sumupdatecnt(self.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_get_recomp(self::CBBundleSolver)

returns the number of oracle reevaluations for the center due to numerical problems since the last descent step
"""
cb_get_recomp(self::CBBundleSolver) = @ccall libcb.cb_bundlesolver_get_recomp(self.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_get_sumrecomp(self::CBBundleSolver)

returns the number of oracle reevaluations for the center due to numerical problems since the last clear() or clear_fails()
"""
cb_get_sumrecomp(self::CBBundleSolver) = @ccall libcb.cb_bundlesolver_get_sumrecomp(self.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_get_qpfails(self::CBBundleSolver)

returns the number of fails in qp subproblems since the last null/descent step
"""
cb_get_qpfails(self::CBBundleSolver) = @ccall libcb.cb_bundlesolver_get_qpfails(self.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_get_sumqpfails(self::CBBundleSolver)

returns the number of fails in qp subproblems since the last clear() or clear_fails()
"""
cb_get_sumqpfails(self::CBBundleSolver) = @ccall libcb.cb_bundlesolver_get_sumqpfails(self.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_get_modelfails(self::CBBundleSolver)

returns the number of fails in model evaluatoins since the last null/descent step
"""
cb_get_modelfails(self::CBBundleSolver) = @ccall libcb.cb_bundlesolver_get_modelfails(self.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_get_summodelfails(self::CBBundleSolver)

returns the number of fails in model evaluations since the last clear() or clear_fails()
"""
cb_get_summodelfails(self::CBBundleSolver) = @ccall libcb.cb_bundlesolver_get_summodelfails(self.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_get_augvalfails(self::CBBundleSolver)

returns the number of failures to increase the augmented model value since the last null/descent step
"""
cb_get_augvalfails(self::CBBundleSolver) = @ccall libcb.cb_bundlesolver_get_augvalfails(self.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_get_sumaugvalfails(self::CBBundleSolver)

returns the number of failures to increase the augmented model value since the last clear() or clear_fails()
"""
cb_get_sumaugvalfails(self::CBBundleSolver) = @ccall libcb.cb_bundlesolver_get_sumaugvalfails(self.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_get_oraclefails(self::CBBundleSolver)

returns the number of fails in oracle evaluations since the last descent step
"""
cb_get_oraclefails(self::CBBundleSolver) = @ccall libcb.cb_bundlesolver_get_oraclefails(self.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_get_sumoraclefails(self::CBBundleSolver)

returns the number of fails in oracle evaluations since the last clear() or clear_fails()
"""
cb_get_sumoraclefails(self::CBBundleSolver) = @ccall libcb.cb_bundlesolver_get_sumoraclefails(self.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_get_shallowcut(self::CBBundleSolver)

returns the number of oracle evaluations that returned an epsilon subgradient that improved the model by a dangerously small amount (mostly this is due to solving the QP subproblems only approximately)
"""
cb_get_shallowcut(self::CBBundleSolver) = @ccall libcb.cb_bundlesolver_get_shallowcut(self.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_get_groundset(self::CBBundleSolver)

return a pointer to the groundset
"""
cb_get_groundset(self::CBBundleSolver) = CBGroundset(@ccall libcb.cb_bundlesolver_get_groundset(self.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_get_model(self::CBBundleSolver)

return a pointer to the cutting model
"""
cb_get_model(self::CBBundleSolver) = CBBundleModel(@ccall libcb.cb_bundlesolver_get_model(self.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_get_prox(self::CBBundleSolver)

return a pointer to the quadratic term of the proximal term
"""
cb_get_prox(self::CBBundleSolver) = CBBundleProxObject(@ccall libcb.cb_bundlesolver_get_prox(self.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_get_terminator(self::CBBundleSolver)

return a pointer to the termination criterion
"""
cb_get_terminator(self::CBBundleSolver) = CBBundleTerminator(@ccall libcb.cb_bundlesolver_get_terminator(self.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_get_bundleweight(self::CBBundleSolver)

return a pointer to the class for updating the weightu of the proximal term
"""
cb_get_bundleweight(self::CBBundleSolver) = CBBundleWeight(@ccall libcb.cb_bundlesolver_get_bundleweight(self.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_get_QPcoeff_time(self::CBBundleSolver)

return time spent in computing the cost coefficients of the quadratic bundle subproblem
"""
cb_get_QPcoeff_time(self::CBBundleSolver) = CBCH_Tools::Microseconds(@ccall libcb.cb_bundlesolver_new_get_qpcoeff_time(self.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_get_QPsolve_time(self::CBBundleSolver)

return time spent in solving the quadratic bundle subproblem
"""
cb_get_QPsolve_time(self::CBBundleSolver) = CBCH_Tools::Microseconds(@ccall libcb.cb_bundlesolver_new_get_qpsolve_time(self.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_get_make_aggr_time(self::CBBundleSolver)

return time spent in forming the model aggregate
"""
cb_get_make_aggr_time(self::CBBundleSolver) = CBCH_Tools::Microseconds(@ccall libcb.cb_bundlesolver_new_get_make_aggr_time(self.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_get_evalaugmodel_time(self::CBBundleSolver)

return time spent in total for the quadratic bundle subproblem
"""
cb_get_evalaugmodel_time(self::CBBundleSolver) = CBCH_Tools::Microseconds(@ccall libcb.cb_bundlesolver_new_get_evalaugmodel_time(self.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_print_statistics(self::CBBundleSolver)

output some time statistic paramters
"""
cb_print_statistics(self::CBBundleSolver) = @ccall libcb.cb_bundlesolver_print_statistics(self.data::Ptr{Cvoid})::Cvoid

@doc raw"""
    cb_qp_mfile_data(self::CBBundleSolver, center_y::CBMatrix, Hp::Union{<:CBBundleProxObject,Nothing}, gs_subg::CBMinorantPointer, Q::CBSymmatrix, c::CBMatrix, offset::Real, yfixed::CBIndexmatrix)

output the data of the Gauss-Seidel qp to in an m file format
"""
cb_qp_mfile_data(self::CBBundleSolver, center_y::CBMatrix, Hp::Union{<:CBBundleProxObject,Nothing}, gs_subg::CBMinorantPointer, Q::CBSymmatrix, c::CBMatrix, offset::Real, yfixed::CBIndexmatrix) = @ccall libcb.cb_bundlesolver_qp_mfile_data(self.data::Ptr{Cvoid}, center_y.data::Ptr{Cvoid}, (isnothing(Hp) ? C_NULL : Hp.data)::Ptr{Cvoid}, gs_subg.data::Ptr{Cvoid}, Q.data::Ptr{Cvoid}, c.data::Ptr{Cvoid}, offset::Cdouble, yfixed.data::Ptr{Cvoid})::Cint

