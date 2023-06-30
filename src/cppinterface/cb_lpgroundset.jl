@doc raw"""
    cb_clear!(self::CBLPGroundset, indim::Integer = 0, in_groundset_id::Integer = 0)

resets all values as described in Groundset::clear()
"""
cb_clear!(self::CBLPGroundset, indim::Integer = 0, in_groundset_id::Integer = 0) = @ccall libcb.cb_lpgroundset_clear(self.data::Ptr{Cvoid}, indim::Cint, in_groundset_id::Cint)::Cvoid

@doc raw"""
    CBLPGroundset()

calls clear() with the same parameters
"""
CBLPGroundset() = CBLPGroundset(@ccall libcb.cb_lpgroundset_new()::Ptr{Cvoid})

@doc raw"""
    CBLPGroundset(dim::Integer, lbyp::Union{<:CBMatrix,Nothing} = nothing, ubyp::Union{<:CBMatrix,Nothing} = nothing, Gp::Union{<:CBSparsemat,Nothing} = nothing, rhslbp::Union{<:CBMatrix,Nothing} = nothing, rhsubp::Union{<:CBMatrix,Nothing} = nothing, start_val::Union{<:CBMatrix,Nothing} = nothing, costs::Union{<:CBMatrix,Nothing} = nothing, offset::Real = 0., in_groundset_id::Integer = 0)

allows to specify the groundset in the constructor, zero is allowed everywhere
"""
CBLPGroundset(dim::Integer, lbyp::Union{<:CBMatrix,Nothing} = nothing, ubyp::Union{<:CBMatrix,Nothing} = nothing, Gp::Union{<:CBSparsemat,Nothing} = nothing, rhslbp::Union{<:CBMatrix,Nothing} = nothing, rhsubp::Union{<:CBMatrix,Nothing} = nothing, start_val::Union{<:CBMatrix,Nothing} = nothing, costs::Union{<:CBMatrix,Nothing} = nothing, offset::Real = 0., in_groundset_id::Integer = 0) = CBLPGroundset(@ccall libcb.cb_lpgroundset_new2(dim::Cint, (isnothing(lbyp) ? C_NULL : lbyp.data)::Ptr{Cvoid}, (isnothing(ubyp) ? C_NULL : ubyp.data)::Ptr{Cvoid}, (isnothing(Gp) ? C_NULL : Gp.data)::Ptr{Cvoid}, (isnothing(rhslbp) ? C_NULL : rhslbp.data)::Ptr{Cvoid}, (isnothing(rhsubp) ? C_NULL : rhsubp.data)::Ptr{Cvoid}, (isnothing(start_val) ? C_NULL : start_val.data)::Ptr{Cvoid}, (isnothing(costs) ? C_NULL : costs.data)::Ptr{Cvoid}, offset::Cdouble, in_groundset_id::Cint)::Ptr{Cvoid})

@doc raw"""
    cb_get_groundset_id(self::CBLPGroundset)

returns the current groundset_id, increased values indicate changes in the ground set
"""
cb_get_groundset_id(self::CBLPGroundset) = @ccall libcb.cb_lpgroundset_get_groundset_id(self.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_set_groundset_id!(self::CBLPGroundset, gsid::Integer)

sets the groundset_id to the desired value, increasing it is safer here because this is used to indicate changes
"""
cb_set_groundset_id!(self::CBLPGroundset, gsid::Integer) = @ccall libcb.cb_lpgroundset_set_groundset_id(self.data::Ptr{Cvoid}, gsid::Cint)::Cvoid

@doc raw"""
    cb_get_dim(self::CBLPGroundset)

returns the dimension of the ground set, i.e., the length of the variables vector y
"""
cb_get_dim(self::CBLPGroundset) = @ccall libcb.cb_lpgroundset_get_dim(self.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_set_qpsolver!(self::CBLPGroundset, qpparams::Union{<:CBQPSolverParametersObject,Nothing}, qpsolver::Union{<:CBQPSolverObject,Nothing} = nothing)

Set the qp solver's parameters to qpparams (if not null); if the second argument qpsolver is also given, the old solver is first discarded and replaced by this new solver and then the parameters are set (if given).  Any object passed here will be owned and deleted by *this. For correct continuaton a new qpsolver needs to have the same feasible set as the current solver but this must be ensured by the caller.
"""
cb_set_qpsolver!(self::CBLPGroundset, qpparams::Union{<:CBQPSolverParametersObject,Nothing}, qpsolver::Union{<:CBQPSolverObject,Nothing} = nothing) = @ccall libcb.cb_lpgroundset_set_qpsolver(self.data::Ptr{Cvoid}, (isnothing(qpparams) ? C_NULL : qpparams.data)::Ptr{Cvoid}, (isnothing(qpsolver) ? C_NULL : qpsolver.data)::Ptr{Cvoid})::Cint

@doc raw"""
    cb_constrained(self::CBLPGroundset)

returns false if the feasible set is the entire space (unconstrained optimization), true otherwise.
"""
cb_constrained(self::CBLPGroundset) = Bool(@ccall libcb.cb_lpgroundset_constrained(self.data::Ptr{Cvoid})::Cint)

@doc raw"""
    cb_get_lby(self::CBLPGroundset)

returns the lower bounds vector on y if it exists
"""
cb_get_lby(self::CBLPGroundset) = CBMatrix(@ccall libcb.cb_lpgroundset_get_lby(self.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_get_uby(self::CBLPGroundset)

returns the upper bounds vector on y if it exists
"""
cb_get_uby(self::CBLPGroundset) = CBMatrix(@ccall libcb.cb_lpgroundset_get_uby(self.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_is_feasible!(self::CBLPGroundset, y::CBMatrix, relprec::Real = 1e-10)

returns true if still feasible, see Groundset::is_feasible()
"""
function cb_is_feasible!(self::CBLPGroundset, y::CBMatrix, relprec::Real = 1e-10)
    in_groundset_id = Ref{Int}()
    Bool(@ccall libcb.cb_lpgroundset_is_feasible(self.data::Ptr{Cvoid}, in_groundset_id::Ref{Int}, y.data::Ptr{Cvoid}, relprec::Cdouble)::Cint)
    return in_groundset_id[]
end

@doc raw"""
    cb_ensure_feasibility!(self::CBLPGroundset, y::CBMatrix, ychanged::Bool, Hp::Union{<:CBBundleProxObject,Nothing}, relprec::Real = 1e-10)

makes y feasible if not so, see Groundset::ensure_feasibility()
"""
function cb_ensure_feasibility!(self::CBLPGroundset, y::CBMatrix, ychanged::Bool, Hp::Union{<:CBBundleProxObject,Nothing}, relprec::Real = 1e-10)
    in_groundset_id = Ref{Int}()
    @ccall libcb.cb_lpgroundset_ensure_feasibility(self.data::Ptr{Cvoid}, in_groundset_id::Ref{Int}, y.data::Ptr{Cvoid}, ychanged::Ref{Cint}, (isnothing(Hp) ? C_NULL : Hp.data)::Ptr{Cvoid}, relprec::Cdouble)::Cint
    return in_groundset_id[]
end

@doc raw"""
    cb_get_qp_solver!(self::CBLPGroundset, solves_model_without_gs::Bool, Hp::Union{<:CBBundleProxObject,Nothing})

returns a pointer to an internal QPSolverObject that is able to solve bundle suproblems efficiently for this kind of groundset and scaling; if solves_model_without_gs == true the qp solver does not include the groundset and the groundset has to be dealt with by the Gauss Seidel approach
"""
cb_get_qp_solver!(self::CBLPGroundset, solves_model_without_gs::Bool, Hp::Union{<:CBBundleProxObject,Nothing}) = CBQPSolverObject(@ccall libcb.cb_lpgroundset_get_qp_solver(self.data::Ptr{Cvoid}, solves_model_without_gs::Ref{Cint}, (isnothing(Hp) ? C_NULL : Hp.data)::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_set_qp_solver_parameters!(self::CBLPGroundset, param0::Union{<:CBQPSolverParametersObject,Nothing})

set parameters for the QP_Solver
"""
cb_set_qp_solver_parameters!(self::CBLPGroundset, param0::Union{<:CBQPSolverParametersObject,Nothing}) = @ccall libcb.cb_lpgroundset_set_qp_solver_parameters(self.data::Ptr{Cvoid}, (isnothing(param0) ? C_NULL : param0.data)::Ptr{Cvoid})::Cint

@doc raw"""
    cb_get_starting_point(self::CBLPGroundset)

returns a stored starting point, note: this need not be feasible; if generated automatically, its dimension is correct.
"""
cb_get_starting_point(self::CBLPGroundset) = (@ccall libcb.cb_lpgroundset_get_starting_point(self.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

@doc raw"""
    cb_set_starting_point!(self::CBLPGroundset, vec::CBMatrix)

stores the a new starting point irrespective of whether it is feasible or not and returns 0 if it feasible, 1 if it is infeasible
"""
cb_set_starting_point!(self::CBLPGroundset, vec::CBMatrix) = @ccall libcb.cb_lpgroundset_set_starting_point(self.data::Ptr{Cvoid}, vec.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_candidate!(self::CBLPGroundset, newy::CBMatrix, center_y::CBMatrix, center_value::Real, model_minorant::CBMinorantPointer, Hp::Union{<:CBBundleProxObject,Nothing}, delta_groundset_minorant::Union{<:CBMinorantPointer,Nothing} = nothing, delta_index::Union{<:CBIndexmatrix,Nothing} = nothing, relprec::Real = 1e-2)

computes the next ground set minorant and candidate, see Groundset::candidate()
"""
function cb_candidate!(self::CBLPGroundset, newy::CBMatrix, center_y::CBMatrix, center_value::Real, model_minorant::CBMinorantPointer, Hp::Union{<:CBBundleProxObject,Nothing}, delta_groundset_minorant::Union{<:CBMinorantPointer,Nothing} = nothing, delta_index::Union{<:CBIndexmatrix,Nothing} = nothing, relprec::Real = 1e-2)
    subgnorm2 = Ref{Float64}()
    augval_ub = Ref{Float64}()
    augval_lb = Ref{Float64}()
    linval = Ref{Float64}()
    cand_gs_val = Ref{Float64}()
    gs_id = Ref{Int}()
    @ccall libcb.cb_lpgroundset_candidate(self.data::Ptr{Cvoid}, gs_id::Ref{Int}, newy.data::Ptr{Cvoid}, cand_gs_val::Ref{Float64}, linval::Ref{Float64}, augval_lb::Ref{Float64}, augval_ub::Ref{Float64}, subgnorm2::Ref{Float64}, center_y.data::Ptr{Cvoid}, center_value::Cdouble, model_minorant.data::Ptr{Cvoid}, (isnothing(Hp) ? C_NULL : Hp.data)::Ptr{Cvoid}, (isnothing(delta_groundset_minorant) ? C_NULL : delta_groundset_minorant.data)::Ptr{Cvoid}, (isnothing(delta_index) ? C_NULL : delta_index.data)::Ptr{Cvoid}, relprec::Cdouble)::Cint
    return gs_id[], cand_gs_val[], linval[], augval_lb[], augval_ub[], subgnorm2[]
end

@doc raw"""
    cb_get_gs_aggregate(self::CBLPGroundset)

returns the groundset aggregate computed in candidate()
"""
cb_get_gs_aggregate(self::CBLPGroundset) = (@ccall libcb.cb_lpgroundset_get_gs_aggregate(self.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

@doc raw"""
    cb_get_gs_minorant(self::CBLPGroundset)

returns the linear minorant valid on the entire ground set (e.g. a linear cost funciton)
"""
cb_get_gs_minorant(self::CBLPGroundset) = (@ccall libcb.cb_lpgroundset_get_gs_minorant(self.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

@doc raw"""
    cb_get_yfixed(self::CBLPGroundset)

if not NULL (iff get_use_yfixing()==false) it returns the vector yfixed with yfixed(i)=0 if not fixed, =1 is fixed already, =2 if newly fixed
"""
cb_get_yfixed(self::CBLPGroundset) = CBIndexmatrix(@ccall libcb.cb_lpgroundset_get_yfixed(self.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_set_yfixed!(self::CBLPGroundset)

if not NULL (iff get_use_yfixing()==false) returns the vector yfixed with yfixed(i)=0 if not fixed, =1 is fixed already, =2 if newly fixed
"""
cb_set_yfixed!(self::CBLPGroundset) = CBIndexmatrix(@ccall libcb.cb_lpgroundset_set_yfixed(self.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_get_use_yfixing(self::CBLPGroundset)

true if the cooridinate fixing heuristic is switched on (only constrained cases)
"""
cb_get_use_yfixing(self::CBLPGroundset) = Bool(@ccall libcb.cb_lpgroundset_get_use_yfixing(self.data::Ptr{Cvoid})::Cint)

@doc raw"""
    cb_set_use_yfixing!(self::CBLPGroundset, uyf::Bool)

set to true to switch on the cooridinate fixing heuristic (only constrained cases)
"""
cb_set_use_yfixing!(self::CBLPGroundset, uyf::Bool) = @ccall libcb.cb_lpgroundset_set_use_yfixing(self.data::Ptr{Cvoid}, uyf::Cint)::Cvoid

@doc raw"""
    cb_set_variable_metric_selection!(self::CBLPGroundset, vms::Union{<:CBVariableMetricSelection,Nothing} = nothing)

delete old selector and set a new one (0 is allowed resulting in no local selector)
"""
cb_set_variable_metric_selection!(self::CBLPGroundset, vms::Union{<:CBVariableMetricSelection,Nothing} = nothing) = @ccall libcb.cb_lpgroundset_set_variable_metric_selection(self.data::Ptr{Cvoid}, (isnothing(vms) ? C_NULL : vms.data)::Ptr{Cvoid})::Cint

@doc raw"""
    cb_get_variable_metric_selection(self::CBLPGroundset)

delete old selector and set a new one (0 is allowed resulting in no local selector)
"""
cb_get_variable_metric_selection(self::CBLPGroundset) = CBVariableMetricSelection(@ccall libcb.cb_lpgroundset_get_variable_metric_selection(self.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_add_variable_metric!(self::CBLPGroundset, H::CBVariableMetric, y_id::Integer, y::CBMatrix, descent_step::Bool, weightu::Real, model_maxviol::Real, indices::Union{<:CBIndexmatrix,Nothing} = nothing)

see VariableMetric
"""
cb_add_variable_metric!(self::CBLPGroundset, H::CBVariableMetric, y_id::Integer, y::CBMatrix, descent_step::Bool, weightu::Real, model_maxviol::Real, indices::Union{<:CBIndexmatrix,Nothing} = nothing) = @ccall libcb.cb_lpgroundset_add_variable_metric(self.data::Ptr{Cvoid}, H.data::Ptr{Cvoid}, y_id::Cint, y.data::Ptr{Cvoid}, descent_step::Cint, weightu::Cdouble, model_maxviol::Cdouble, (isnothing(indices) ? C_NULL : indices.data)::Ptr{Cvoid})::Cint

@doc raw"""
    cb_start_modification!(self::CBLPGroundset)

propagates the call to QPSolverObject::QPstart_modification() of the current qpsolver
"""
cb_start_modification!(self::CBLPGroundset) = CBGroundsetModification(@ccall libcb.cb_lpgroundset_start_modification(self.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_apply_modification!(self::CBLPGroundset, mdf::CBGroundsetModification)

change the groundset description as specified by the argument
"""
cb_apply_modification!(self::CBLPGroundset, mdf::CBGroundsetModification) = @ccall libcb.cb_lpgroundset_apply_modification(self.data::Ptr{Cvoid}, mdf.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_mfile_data(self::CBLPGroundset)

m-file output routine for debugging or testing in Matlab (not yet working)
"""
cb_mfile_data(self::CBLPGroundset) = @ccall libcb.cb_lpgroundset_mfile_data(self.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_set_cbout!(self::CBLPGroundset, incr::Integer = -1)

output settings
"""
cb_set_cbout!(self::CBLPGroundset, incr::Integer = -1) = @ccall libcb.cb_lpgroundset_set_cbout(self.data::Ptr{Cvoid}, incr::Cint)::Cvoid

