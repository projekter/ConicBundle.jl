@doc raw"""
    module ConicBundle

Solve ``\min_{y\in\mathbf{R}^m}  f_0(y) + f_1(y) + ... + f_k(y)`` for convex functions ``f_i``, the y-variables may be
bounded or box constrained. The most important steps are explained here.
Internal details are sketched in internal_cinterface.

# Setting up the Problem, the Functions, and the Main Loop
First open a new problem by constructing a [`CBProblem`](@ref) object [use `CBProblem(1)` in order to employ a minimal bundle
solver with just one aggregate and one new subgradient in each iteration; this is an attractive choice, if fast iterations
and/or little memory consumption are of special importance]. The object will be needed for every manipulation of this problem.
Cleanup is performed automatically.

Next, set the dimension of the design variables/argument as well as possible box constraints on these by the function
[`cb_init_problem!`](@ref cb_init_problem!(::CBProblem, ::Integer)).

Now set up your functions ``f_i`` as functions that respect the call signature detailed in [`CBFunction`](@ref). Via these
functions you will supply, for a given argument, the function value and a subgradient (=the gradient if the function is
differentiable) to the solver.

The callbacks have to be added to the solver using the routine
[`cb_add_function!`](@ref cb_add_function!(::CBProblem{T}, ::T, ::Function, ::Nothing) where {T}).

Once all functions are added, the optimization process can be started. If you know a good starting point then set it with
[`cb_set_new_center_point!`](@ref cb_set_new_center_point!(::CBProblem, ::Vector{Float64})) now, otherwise the method will pick
the zero vector or, in the case of box constraints, the point closest to zero as starting point.

Finally, set up a loop that calls [`cb_solve!`](@ref cb_solve!(::CBProblem, ::Integer, ::Bool)) until
[`cb_termination_code`](@ref cb_termination_code(::CBProblem)) is nonzero.

After the first call to [`cb_solve!`](@ref cb_solve!(::CBProblem, ::Integer, ::Bool)) you can retrieve, at any time, the
current objective value by [`cb_get_objval`](@ref cb_get_objval(::CBProblem)) and the argument leading to this value by
[`cb_get_center`](@ref cb_get_center(::CBProblem)). For some screen output, use
[`cb_set_print_level!`](@ref cb_set_print_level!(::CBProblem, ::Integer)).

# Lagrangean Relaxation, Primal Approximations, and Cutting Planes
If you are optimizing a Lagrangean relaxation, you might be interested in getting an approximation to your primal optimal
solution. This can be done by specifying in each function for each (epsilon) subgradient the corresponding primal vectors that
generate it, see [`CBFunction`](@ref) and
[`cb_add_function!`](@ref cb_add_function!(::CBProblem{T}, ::T, ::Function, ::Nothing) where {T}) as a start. Then for each of
your functions, you can retrieve the current primal approximation using
[`cb_get_approximate_primal!`](@ref cb_get_approximate_primal!(::Vector{Float64}, ::CBProblem{T}, ::T) where {T}).

If, in addition, you plan to improve your primal relaxation via cutting planes, that are strongly violated by the current
primal approximation, you should have a look at [`cb_append_variables!`](@ref cb_append_variables!(::CBProblem, ::Integer)),
[`CBSubgExt`](@ref), [`cb_reinit_function_model!`](@ref cb_reinit_function_model!(::CBProblem{T}, ::T) where {T}),
[`cb_get_approximate_slacks`](@ref cb_get_approximate_slacks(::CBProblem)), and
[`cb_delete_variables!`](@ref cb_delete_variables!(::CBProblem, ::Vector{Cint})).
"""
module ConicBundle

import LinearAlgebra
import MutableArithmetics as MA
using Pkg.Artifacts

export CBProblem, cb_clear!, cb_set_default!, cb_init_problem!, cb_add_function!, cb_set_lower_bound!, cb_set_upper_bound!,
    cb_append_variables!, cb_delete_variables!, cb_reassign_variables!, cb_solve!, cb_termination_code,
    cb_print_termination_code, cb_get_objval, cb_get_center, cb_get_center!, cb_get_sgnorm, cb_get_subgradient,
    cb_get_subgradient!, cb_get_candidate_value, cb_get_candidate, cb_get_candidate!, cb_set_term_relprec!,
    cb_set_new_center_point!, cb_get_function_status, cb_get_approximate_slacks, cb_get_approximate_slacks!,
    cb_get_approximate_primal!, cb_get_center_primal!, cb_get_candidate_primal!, cb_set_max_modelsize!, cb_set_max_bundlesize!,
    cb_set_max_new_subgradients!, cb_get_bundle_parameters, cb_reinit_function_model!, cb_get_last_weight, cb_set_next_weight!,
    cb_set_min_weight!, cb_set_max_weight!, cb_set_variable_metric!, cb_get_dim, cb_get_n_functions, cb_get_minus_infinity,
    cb_get_plus_infinity, cb_clear_fail_counts!, cb_set_eval_limit!, cb_set_inner_update_limit!, cb_set_active_bounds_fixing!,
    cb_get_fixed_active_bounds, cb_get_fixed_active_bounds!, cb_set_print_level!, cb_minus_infinity, cb_plus_infinity,
    cbft_objective_function, cbft_constant_penalty_function, cbft_adaptive_penalty_function,
    cbt_not_terminated, cbt_relprec_satisfied, cbt_time_limit, cbt_feval_limit, cbt_subfailure_limit, cbt_mfailure_limit,
    cbt_incfailure_limit, cbt_ocall_limit, cbt_ofailure_limit,
    cbvm_no_scaling, cbvm_diagonal_scaling, cbvm_diagonal_scaling_with_bounds

const libcb = joinpath(artifact"ConicBundle", Sys.iswindows() ? "ConicBundle.dll" : "ConicBundle.so")

mutable struct CBCSolver end;

const CBHandle = Ptr{CBCSolver};

"""
    CBProblem{T}(no_bundle::Bool=false)
    CBProblem(no_bundle::Bool=false)

Creates a a new problem object and returns a pointer to it.

Arguments:
- `T`

  type of the function keys (`Any` if omitted)
- `no_bundle::Bool`

  if `true`, then the minimal bundle consisting of just one new and one aggregate gradient is used so that there is no real
  bundle available and bundle size options are then meaningless.
"""
struct CBProblem{T}
    handle::CBHandle
    functions::Dict{T,Base.CFunction}

    function CBProblem{T}(no_bundle::Bool=false) where {T}
        handle = @ccall libcb.cb_construct_problem(no_bundle::Cint)::CBHandle
        handle === C_NULL && error("Construction of problem failed")
        functions = Dict{T,Base.CFunction}()
        finalizer(functions) do _
            @ccall libcb.cb_destruct_problem(Ref(handle)::Ref{CBHandle})::Cint
        end
        new{T}(handle, functions)
    end

    CBProblem(no_bundle::Bool=false) = CBProblem{Any}(no_bundle)
end

include("cppinterface/cb_classes.jl")

Base.convert(::Type{P}, ::Nothing) where {P<:Ptr} = P(C_NULL)
Base.unsafe_convert(::Type{P}, ::Nothing) where {P<:Ptr} = P(C_NULL)

"""
    cb_clear!(p::CBProblem)

Clears all data structures and problem information but keeps ouptut settings and algorithmic parameter settings.
"""
cb_clear!(p::CBProblem) = @ccall libcb.cb_clear(p.handle::CBHandle)::Cvoid

"""
    cb_set_default!(p::CBProblem)

Sets default values for algorithmic parameters that are not function specific (e.g., relative precision, weight and weight
bounds for the augmented problem, etc.)
"""
cb_set_default!(p::CBProblem) = @ccall libcb.cb_set_defaults(p.handle::CBHandle)::Cvoid

@doc raw"""
    cb_init_problem!(p::CBProblem, m::Integer; lowerb=nothing, upperb=nothing)

Initializes the problem by setting the design space (the dimension and possible box constraints of the variables)

Clears all data structures and sets the dimension ``m`` for a new problem.
for solving ``\min_{y in \mathbb R^m}  f_0(y) + f_1(y) + ...``
Box constraints may be specified for ``y``, the functions ``f_i`` must be added by
  [`cb_add_function!`](@ref cb_add_function!(::CBProblem{T}, ::T, ::Function, ::Nothing) where {T}).

Lower and/or upper bounds must be specified for all variables or for none of them. To specify no bounds at all, give `nothing`.
Otherwise use [`cb_get_minus_infinity`](@ref) for unbounded below and [`cb_get_plus_infinity`](@ref) for unbounded above.
For `nothing`, unbounded will be used as default for all variables. Specifying bounds selectively is also possible
by [`cb_set_lower_bound!`](@ref cb_set_lower_bound!(::CBProblem, ::Integer, ::Float64)) or
[`cb_set_upper_bound!`](@ref cb_set_upper_bound!(::CBProblem, ::Integer, ::Float64)).

Arguments:
- `p::CBProblem`

  the current problem
- `m::Integer`

  the dimension of the argument/design space/the number of Lagrange multipliers
- `lowerb::Union{Vector{Float64},Nothing}` (either double array of length `m` or `nothing`)

  If `nothing`, all variables are considered unbounded below, otherwise `lowerb[i]` gives the minimum feasible value for
  variable `y[i]`, use [`cb_get_minus_infinity`](@ref) for unbounded below.
- `upperb::Union{Vector{Float64},Nothing}` (either double array of length `m` or `nothing`)

  If `nothing`, all variables are considered unbounded above, otherwise `upperb[i]` gives the maximum feasible value for
  variable `y[i]`, use [`cb_get_plus_infinity`](@ref) for unbounded above.
"""
Base.@propagate_inbounds function cb_init_problem!(p::CBProblem, m::Integer; lowerb::Union{Vector{Float64},Nothing}=nothing,
    upperb::Union{Vector{Float64},Nothing}=nothing, offset::Float64=0.)
    @boundscheck isnothing(lowerb) || length(lowerb) == m || error("lowerb must be of the specified length m")
    @boundscheck isnothing(upperb) || length(upperb) == m || error("upperb must be of the specified length m")
    GC.@preserve lowerb upperb begin
        iszero(@ccall libcb.cb_init_problem(
            p.handle::CBHandle,
            m::Cint,
            lowerb::Ptr{Cdouble},
            upperb::Ptr{Cdouble},
            offset::Cdouble
        )::Cint) || error("Initialization of problem failed")
    end
    return
end

@doc raw"""
    CBFunction{T}(callback, m, n)

function oracle; describe your function as a function of this type to pass it to the solver

The oracle interface is used to describe a convex function. The dimension of the argument vector of the function must be set in
[`cb_init_problem!`](@ref cb_init_problem!(::CBProblem, ::Integer)), let it be `m` in the following.

If the sum of several such functions is to be minimized, it is the task of the user to guarantee that all dimensions match.

In many applications, computing the function value is an iterative process that approaches the true function value from below.
The code offers a bound for the function value, above which it is certain that the code will reject the current point. If in
the iterative process a lower bound on the function value exceeds this bound, then it is sufficient to return, instead of the
true function value and a subgradient, the current lower bound and a vector so that together they describe a supporting
hyperplane (an epsilon subgradient, lying completely below the function) to the function at this point.

If the function corresponds to Lagrangean relaxation of a primal maximization problem one may want to generate a primal
approximate solution. In this case, set in
[`cb_add_function!`](@ref cb_add_function!(::CBProblem{T}, ::T, ::Function, ::Nothing) where {T}) the desired primal dimension.
Then the solver will provide memory in primal for returning in the function the generating primal vectors for each subgradient.
If the primal dimension is set to zero, primal will be `nothing` and no aggregation takes place.

If primal aggregation is used then it is possible to implement a primal cutting plane framework. This requires the introduction
of new (dual) variables in the design space of the function. In this case a function of signature [`CBSubgExt`](@ref) should be
provided for filling in the missing coordinates in existing subgradients. This function need not be specified even if
constraints are added, but then the cutting model of the objective is lost at each addition of constraints.

Arguments of the callback:
- `function_key::T`

  supplied by the user on entering the function. May be useful if multiple copies of the same function are used with parameters
- `arg::Vector{Float64}` (double array of length `m`)

  argument of the function (e.g. the Lagrange multipliers)
- `relprec::Float64`

  relative precision requirement for objective values that may lead to descent steps (this precision is not required if it is
  already certain that the function value will be too poor)
- `max_subg::Cint`

  at most `max_subg` epsilon-`subgradients` and `subg_values` may be returned, but at least one must be returned!
- `objective_threshold::Float64`

  value gives the threshold for a null step; you may stop, if a cutting plane yields at least this;
- `subg_values::Vector{Float64}` (caller-allocated, length `max_subg`)

  store for each epsilon subgradient the value at the argument
- `subgradients::Matrix{Float64}` (caller-allocated, size `(m, max_subg)`)

  `subgradients[i, j]` = coefficient of subgradient `j` at y-coordinate `i`
- `primal::Union{Matrix{Float64},Nothing}` (caller-allocated, size `(n, max_subg)` or `nothing`)

  If the function arises from Lagrangean relaxation and a primal approximation is desired then set the primal dimension in
  [`cb_add_function!`](@ref cb_add_function!(::CBProblem{T}, ::T, ::Function, ::Nothing) where {T}) and return the primal
  solutions corresponding to the eps-subgradients in the array pointed to by primal.

Callback must return a tuple:
- upper bound on the true function value within ``\mathrm{relprec}*(\operatorname{abs}(\mathrm{objval})+1.)``, if there is no
  hyperplane cutting above the threshold specified in `objective_value` on input.
  Otherwise the max of `subg_values`.
- the number of epsilon-`subgradients` returned. Termination is forced if no new subgradient is returned.
"""
struct CBFunction{T}
    callback::Function
    key::T
    m::Cint
    n::Cint
end

function (cbf::CBFunction)(::Ptr{Cvoid}, arg::Ptr{Cdouble}, relprec::Cdouble, max_subg::Cint,
    objective_value::Ptr{Cdouble}, n_subgrads::Ptr{Cint}, subg_values::Ptr{Cdouble}, subgradients::Ptr{Cdouble},
    primal::Ptr{Cdouble})::Cint
    out_objective, out_nsubgrads = cbf.callback(
        cbf.key,
        unsafe_wrap(Array, arg, cbf.m, own=false),
        relprec,
        max_subg,
        unsafe_load(objective_value),
        unsafe_wrap(Array, subg_values, max_subg, own=false),
        unsafe_wrap(Array, subgradients, (cbf.m, max_subg), own=false),
        primal === C_NULL ? nothing : unsafe_wrap(Array, primal, (cbf.n, max_subg), own=false)
    )
    unsafe_store!(objective_value, out_objective)
    unsafe_store!(n_subgrads, out_nsubgrads)
    return zero(Cint)
end

"""
    CBSubgExt{T}(callback, n)

This routine is not needed unless variables (constraints in Lagrangean relaxation) are added on the fly.

The solver calls this routine whenever new variables have been added on the fly in order to extend old subgradients to the new
coordinates. If primal data was supplied for the subgradients then `generating_primal` holds a pointer to this (possibly
aggregated) data, otherwise it is `nothing`.

In the presence of primal data, the new coordinates correspond to the violation of the new primal constraints. These have to be
returned in the array `new_subgradient_values`; more precisely, for i=1 to `n_indices` the element `new_subgradient_values[i]`
has to hold the subgradient information of constraint `variable_indices[i]`;

If `generating_primal` is `nothing`, then the routine can only successfully extend the subgradients, if the new coordinates
have no influence on the function; then the new subgradient coordinates are all zero and the components of
`new_subgradient_values` have to be initialized to zero.

If you do indeed need this, you have to provide one such function with each evaluation function.

Arguments of the callback:
- `function_key::T`

  supplied by the user on entering the function. May be useful if multiple copies of the same function are used with parameters
- `generating_primal::Union{Nothing,Vector{Float64}}` (double array of primal length `n` or `nothing`)
  if not `nothing` it holds the (possibly aggregated) primal solution that generated the subgradient that needs to be extendend
- `n_indices::Cint`

  gives the number of indices for which the subgradient value has to be computed
- `variable_indices::Vector{Int}` (pointer to int array of length `n_indices`)

  for the `y` variables with indices `variable_indices[i]`, `i=1,..,n_indices` the subgradient coefficient has to be computed
- `new_subgradient_values::Vector{Float64}` (caller-allocated, pointer to double array of length `n_indices`)

  store the the subgradient coefficient of `y` variable with index `variable_indices[i]` at `new_subgradient_values[i]` for
  `i=1,..,n_indices`

Callback must return 0 on success or 1 if extension is impossible
"""
struct CBSubgExt{T}
    callback::Function
    key::T
    n::Cint
end

function (cbs::CBSubgExt)(::Ptr{Cvoid}, generating_primal::Ptr{Cdouble}, n_indices::Cint,
    variable_indices::Ptr{Cint}, new_subgradient_values::Ptr{Cdouble})::Cint
    return cbs.callback(
        cbs.key,
        generating_primal === C_NULL ? nothing : unsafe_wrap(Array, generating_primal, cbs.n, own=false),
        n_indices,
        unsafe_wrap(Array, variable_indices, n_indices, own=false),
        unsafe_wrap(Array, new_subgradient_values, n_indices, own=false)
    )
end

@doc """
    enum CBFunctionTask

- `cbft_objective_function`
- `cbft_constant_penalty_function`
- `cbft_adaptive_penalty_function`
"""
@enum CBFunctionTask::Cint cbft_objective_function=0 cbft_constant_penalty_function=1 cbft_adaptive_penalty_function=2

"""
    cb_add_function!(p::CBProblem{T}, function_key::T, f::Function, se::Union{Function,Nothing}=nothing; primaldim::Integer=0,
        fun_factor::Float64=1., fun_task::CBFunctionTask=cbft_objective_function,
        aft::Union{CBAffineFunctionTransformation,Nothing}=nothing)

Adds the function, the sum of which should be minimized, to the problem description.

Each function added must be given a unique `function_key`, `f` supplies the evaluation function and must not be zero, `se` can
be used to specify a routine for extending subgradients, but it may be `nothing`.
`primaldim` can be used if an approximate primal solution should  be aggregated (In this case storage will be supplied in the
call to the evaluation function for storing for each subgradient the generating primal vector). It may be zero if this is not
needed.

Arguments:
- `p::CBProblem{T}`

  the problem to which the function should be added
- `function_key::T`

  The value of the funciton_key must UNIQUELY identify the function, (it may be the address of the function [if unique], or
  give the address of a struct holding additional user parameters)
- `f::Function`

  the function (see [`CBFunction`](@ref) for parameters)
- `se::Union{Function,Nothing}`

  This parameter may be `nothing`, otherwise the respective function will be called in order to compute coefficients for new
  subgradient coordinates resulting from added variables. See [`CBSubgExt`](@ref) for parameters.
- `primaldim::Integer`

  May be zero, otherwise in each call to `f` enough storage will be provided to store a primal generating vector for each
  subgradient returned. The primal solutions will be aggregated along with the subgradients. This allows to generate
  approximate primal optimal solutions, e.g., in Lagrangean relaxation.
- `fun_factor::Float64`
  Allows to specify a scaling factor for the function. `fun_factor` must be a strictly positive number.
- `fun_task::CBFunctionTask`
  Specifies whether the function is to be used as a an ObjectiveFunction, a ConstantPenaltyFunction with `fun_factor` as
  maximum penalty factor, or as an `AdaptivePenaltyFunction` with `fun_factor` at initial penalty guess that might be increased
  or decreased over time.
- `aft::Union{CBAffineFunctionTransformation,Nothing}`
  May be used to modify the argument and give an additional affine term (linear term plus offset). For adding an affine term
  there are several other possibilities, e.g. in [`cb_init_problem!`](@ref cb_init_problem!(::CBProblem, ::Integer)), so there
  is no need to do so here. If, however, an existing function implementation requires only some subset of the variables, it is
  more convenient to supply a corresponding `aft` instead of reimplementing the function.

See also [`CBFunction`](@ref), [`CBFunctionTask`](@ref), [`CBAffineFunctionTransformation`](@ref).
"""
function cb_add_function!(p::CBProblem{T}, function_key::T, f::Function, ::Nothing=nothing; primaldim::Integer=0,
    fun_factor::Float64=1., fun_task::CBFunctionTask=cbft_objective_function,
    aft::Union{CBAffineFunctionTransformation,Nothing}=nothing) where {T}
    haskey(p.functions, function_key) && error("function_key must be unique")
    fun_factor ≤ 0. && error("fun_factor must be positive")
    callback = @cfunction(
        $(CBFunction{T}(f, function_key, cb_get_dim(p), primaldim)),
        Cint,
        (Ptr{Cvoid}, Ptr{Cdouble}, Cdouble, Cint, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble})
    )
    p.functions[function_key] = callback # need to preserve this from GC
    iszero(@ccall libcb.cb_add_function(
        p.handle::CBHandle,
        callback::Ptr{Cvoid},
        callback::Ptr{Cvoid},
        C_NULL::Ptr{Cvoid},
        primaldim::Cint,
        fun_factor::Cdouble,
        fun_task::Cint,
        (isnothing(aft) ? C_NULL : aft.data)::Ptr{Cvoid}
    )::Cint) || error("Adding a function to the problem failed")
    return
end

struct SEFunctionKey{T}
    key::T
end

function cb_add_function!(p::CBProblem{T}, function_key::T, f::Function, se::Function; primaldim::Integer=0,
    fun_factor::Float64=1., fun_task::CBFunctionTask=cbft_objective_function,
    aft::Union{CBAffineFunctionTransformation,Nothing}=nothing) where {T}
    haskey(p.functions, function_key) && error("function_key must be unique")
    fun_factor ≤ 0. && error("fun_factor must be positive")
    @assert(!haskey(p.functions, SEFunctionKey{T}(function_key)))
    callback = @cfunction(
        $(CBFunction{T}(f, function_key, cb_get_dim(p), primaldim)),
        Cint,
        (Ptr{Cvoid}, Ptr{Cdouble}, Cdouble, Cint, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble})
    )
    p.functions[function_key] = callback
    se_callback = @cfunction(
        $(CBSubgExt{T}(se, function_key, primaldim)),
        Cint,
        (Ptr{Cvoid}, Ptr{Cdouble}, Cint, Ptr{Cint}, Ptr{Cdouble})
    )
    p.functions[SEFunctionKey{T}(function_key)] = se_callback # need to preserve these from GC
    iszero(@ccall libcb.cb_add_function(
        p.handle::CBHandle,
        callback::Ptr{Cvoid},
        callback::Ptr{Cvoid},
        se_callback::Ptr{Cvoid},
        primaldim::Cint,
        fun_factor::Cdouble,
        fun_task::Cint,
        (isnothing(aft) ? C_NULL : aft.data)::Ptr{Cvoid}
    )::Cint) || error("Adding a function to the problem failed")
    return
end

"""
    cb_set_lower_bound!(p::CBProblem, i::Integer, lower_bound::Float64)

set lower bound for variable `i`, use [`cb_get_minus_infinity`](@ref) for unbounded from below.

The algorithm may have to adapt the center point aftwards. In this case the old function values will be marked as outdated and
will be recomputed at the next call to e.g. [`cb_solve!`](@ref cb_solve!(::CBProblem, ::Integer, ::Bool)).

Arguments:
- `p::CBProblem`

  the problem
- `i::Integer`

  index of the variable
- `lower_bound::Float64`

  value of the lower bound on variable `i`
"""
function cb_set_lower_bound!(p::CBProblem, i::Integer, lower_bound::Float64)
    iszero(@ccall libcb.cb_set_lower_bound(
        p.handle::CBHandle,
        i::Cint,
        lower_bound::Cdouble
    )::Cint) || error("Setting the lower bound failed")
    return
end

"""
    cb_set_upper_bound!(p::CBProblem, i::Integer, lower_bound::Float64)

set upper bound for variable `i`, use [`cb_get_plus_infinity`](@ref) for unbounded from below.

The algorithm may have to adapt the center point aftwards. In this case the old function values will be marked as outdated and
will be recomputed at the next call to e.g. [`cb_solve!`](@ref cb_solve!(::CBProblem, ::Integer, ::Bool)).

Arguments:
- `p::CBProblem`

  the problem
- `i::Integer`

  index of the variable
- `lower_bound::Float64`

  value of the upper bound on variable `i`
"""
function cb_set_upper_bound!(p::CBProblem, i::Integer, upper_bound::Float64)
    iszero(@ccall libcb.cb_set_upper_bound(
        p.handle::CBHandle,
        i::Cint,
        upper_bound::Cdouble
    )::Cint) || error("Setting the upper bound failed")
    return
end

"""
    cb_append_variables!(p::CBProblem, n_append::Integer; lowerb=nothing, upperb=nothing)

Append new variables (always in last postions in this order).

If 0 is feasible for the new coordinates then this is selected as starting value for the new coordinates; otherwise, the number
closest to zero is used. If all new coordinates can be set to zero then it is assumed that for an existing center point the
function values need not be recomputed (this is e.g. the case in Lagrangean relaxation; if this is not correct call
[`cb_reinit_function_model!`](@ref cb_reinit_function_model!(::CBProblem{T}, ::T) where {T}) below). Otherwise the old function
values will be marked as outdated and will be recomputed at the next call to e.g.
[`cb_solve!`](@ref cb_solve!(::CBProblem, ::Integer, ::Bool)).

Be sure to update your objective functions so that they can handle the new variables before you call this and any further
ConicBundle routines that require function evaluations. Also, these operations may lead to inavailability of certain other data
such as subgradients and primal approximations.

Arguments:
- `p::CBProblem`

  the problem
- `n_append::Integer`

  number of variables to append (always in last position in the same order)
- `lower_bound::Union{Vector{Float64},Nothing}` (either double array of size `n_append` or `nothing`)

  If `nothing`, all appended variables are considered unbounded below, otherwise `lowerb[i]` gives the minimum feasible value
  for variable `y[i]`, use [`cb_get_minus_infinity`](@ref) for unbounded below.
- `upper_bound::Union{Vector{Float64},Nothing}` (either double array of size `n_append` or `nothing`)

  If `nothing`, all appended variables are considered unbounded above, otherwise `upperb[i]` gives the maximum feasible value
  for variable `y[i]`, use [`cb_get_plus_infinity`](@ref) for unbounded above.
"""
function cb_append_variables!(p::CBProblem, n_append::Integer; lowerb::Union{Vector{Float64},Nothing}=nothing,
    upperb::Union{Vector{Float64},Nothing}=nothing)
    GC.@preserve lowerb upperb begin
        iszero(@ccall libcb.cb_append_variables(
            p.handle::CBHandle,
            n_append::Cint,
            lowerb::Ptr{Cdouble},
            upperb::Ptr{Cdouble}
        )::Cint) || error("Appending variables to problem failed")
    end
    return
end

"""
    cb_delete_variables!(p::CBProblem, delete_indices::Vector{Cint})
    cb_delete_variables!(map_to_old::Vector{Cint}, p::CBProblem, delete_indices::Vector{Cint})

Deletes variables corresponding to the specified indices.

The indices of the remaining variables are reassigned so that they are consecutive again, the routine returns a vector giving
for each new index of these remaining variables the old coordinate.

If all of the deleted variables are zero, function values are assumed to remain correct (if this is not so, call
[`cb_reinit_function_model!`](@ref cb_reinit_function_model!(::CBProblem{T}, ::T) where {T}) below) Otherwise the old function
values will be marked as outdated and will be recomputed at the next call to e.g.
[`cb_solve!`](@ref cb_solve!(::CBProblem, ::Integer, ::Bool)).

Be sure to update your objective functions so that they can handle the new variables before you call any further ConicBundle
routines that require function evaluations. Also, these operations may lead to inavailability of certain other data such as
subgradients and primal approximations.

Arguments:
- `p::CBProblem`

  the problem
- `delete_indices::Vector{Cint}`

  the entries in delete_indices specify the indices of the variables to be deleted

Returns a `Vector{Cint}` of length `cb_get_dim(p)-length(delete_indices)` whose element `i` contains the old index (before the
call) of the variable that now has index position `i`.
Use the three-argument version to avoid allocation.
"""
cb_delete_variables!(p::CBProblem, delete_indices::Vector{Cint}) =
    @inbounds cb_delete_variables!(Vector{Cint}(undef, cb_get_dim(p) - length(delete_indices)), p, delete_indices)

Base.@propagate_inbounds function cb_delete_variables!(map_to_old::Vector{Cint}, p::CBProblem, delete_indices::Vector{Cint})
    n_del::Cint = length(delete_indices)
    @boundscheck length(map_to_old) == cb_get_dim(p) - n_del || error("Length of delete_indices is not appropriate")
    iszero(@ccall libcb.cb_delete_variables(
        p.handle::CBHandle,
        n_del::Cint,
        delete_indices::Ref{Cint},
        map_to_old::Ref{Cint}
    )::Cint) || error("Deleting variables from problem failed")
    return map_to_old
end

"""
    cb_reassign_variables!(p::CBProblem, assign_new_from_old::Vector{Cint})

Reassigns variables to new index positions by mapping to position `i` the variable that previously had index
`assign_new_from_old[i]`.

Old variables, that are not mapped to any position will be deleted. It is allowed to generate several copies of old variables.

If all of the deleted variables as well as new multiple copies are zero, function values are assumed to remain correct (if this
is not so, call [`cb_reinit_function_model!`](@ref cb_reinit_function_model!(::CBProblem{T}, ::T) where {T}) below). Otherwise
the old function values will be marked as outdated and will be recomputed at the next call to e.g.
[`cb_solve!`](@ref cb_solve!(::CBProblem, ::Integer, ::Bool)).

Be sure to update your objective functions so that they can handle the new variables before you call any further ConicBundle
routines that require function evaluations. Also, these operations may lead to inavailability of certain other data such as
subgradients and primal approximations.

Arguments:
- `p::CBProblem`

  the problem
- `assign_new_from_old::Vector{Cint}`

  entry `assign_new_from_old[i]` specifies the old index of the variable, that has to be copied to index position `i`.
"""
function cb_reassign_variables!(p::CBProblem, assign_new_from_old::Vector{Cint})
    iszero(@ccall libcb.cb_reassign_variables(
        p.handle::CBHandle,
        length(assign_new_from_old)::Cint,
        assign_new_from_old::Ref{Cint}
    )::Cint) || error("Reassigned variables in problem failed")
    return
end

"""
    cb_solve!(p::CBProblem, maxsteps::Integer=0, stop_at_descent_steps::Bool=false)

Bundle methods solve a problem by a sequence of so called descent steps that actually bring progress by moving from the current
"center point" to a new center with better objective. A descent step may consist of several function evaluations (null steps),
that lead to no immediate progress but mainly improve a cutting model of the objective function close to the current center
point. A minimizer to the model is accepted as descent step if the function value at this point satisfies a sufficient decrease
criterion in comparison to the decrease predicted by the model. Having found a descent step, the next center is automatically
shifted to this successful candidate. Termination criteria may stop the process of seeking for a descent step, in which case
the current center is kept and the routine [`cb_termination_code`](@ref cb_termination_code(::CBProblem)) returns the
termination code.

Restarting, after each descent step, the bundle method from scratch with the new center as starting point does not endanger
convergence. Therefore, a descent step is the smallest unit, after which user interaction can take place safely. To allow this
there is a flag `stop_at_descent_steps` that will cause the code to return after the next descent step.

If you know what your are doing, you may also use the input parameter `maxsteps` to force the algorithm to return after at most
`maxsteps` null steps. Calling solve again without any intermediate problem configurations will then simply continue the
process where it stopped and convergence is save. During null steps one may not decrease the weight or delete nonzero variables
of the center or the current candidate!

In a Lagrangean relaxation cutting plane approach one may want to separate and enlarge the dimension after a certain number of
null steps. In this case the code will try to preserve the model, given appropriate subgradient extension routines have been
provided. If the model cannot be extended, it has to be discarded (if subgradient extension is not successful this is done
automatically), and the algorithm will be restarted from the current center point.

Arguments:
- `p::CBProblem`

  the problem
- `maxsteps::Integer`

  use value `=0` as default (anything <= serves as infinite upper bound), if `maxsteps>0` the code returns after at most so
  many null steps
- `stop_at_descent_steps::Integer`

  if `true` the code also returns whenever a descent step occured, if `false` it only stops after maxsteps or when a
  termination criterion is met
"""
function cb_solve!(p::CBProblem, maxsteps::Integer=0, stop_at_descent_steps::Bool=false)
    iszero(@ccall libcb.cb_solve(
        p.handle::CBHandle,
        maxsteps::Cint,
        stop_at_descent_steps::Cint
    )::Cint) || error("Solver step failed")
    return
end

@doc """
    enum CBTerminationCode

- `cbt_not_terminated`: Not terminated.

  (Continue with the next [`cb_solve!`](@ref cb_solve!(::CBProblem, ::Integer, ::Bool)))
- `cbt_relprec_satisfied`: Relative precision criterion satisfied.

  (See [`cb_set_term_relprec!`](@ref cb_set_term_relprec!(::CBProblem, ::Float64)))
- `cbt_time_limit`: Timelimit exceeded.

  (Currently the C interface does not offer a timelimit.)
- `cbt_feval_limit`: Maximum number of function reevaluations exceeded.

  (Indicates that there is a problem with one of the function oracles that seems to deliver no valid upper bounds on the true
  function value for descent steps)
- `cbt_subfailure_limit`: Maximum number of quadratic subproblem failures exceeded.

  (Indicates that the numerical limits of the inner quadratic programming solver are reached, no further progress expected)
- `cbt_mfailure_limit`: maximum number of model evaluation failures exceeded

  (Indicates that the numerical limits of the setup of the subproblem are reached, no further progress expected)
- `cbt_incfailure_limit`: maximum number of failures to increase the augmented model value exceeded

  (Indicates that the numerical limits of the interplay between subproblem and quadratic programming solver are reached, no
  further progress expected)
- `cbt_ocall_limit`: maximum number of oracle calls (function evaluations) exceeded, see
  [`cb_set_eval_limit!`](@ref cb_set_eval_limit!(::CBProblem, ::Integer))
- `cbt_ofailure_limit`: maximum number of oracle failures exceeded.

  This refers to function evaluations that terminate with insufficient precision but still provide a new approximate
  subgradient. A failure typically indicates numerical difficulties with the precision requirements.
  (Currently the C interface does not allow to manipulate the limit, it is set to 10)
"""
@enum(CBTerminationCode::Cint, cbt_not_terminated=0, cbt_relprec_satisfied=1, cbt_time_limit=2, cbt_feval_limit=4,
    cbt_subfailure_limit=8, cbt_mfailure_limit=16, cbt_incfailure_limit=32, cbt_ocall_limit=64, cbt_ofailure_limit=128)

"""
    cb_termination_code(p::CBProblem)

Returns the termination code of the bundle algorithm for the latest descent step

For resetting all counters relevant for termination see [`cb_clear_fail_counts!`](@ref cb_clear_fail_counts!(::CBProblem)).

See also [`CBTerminationCode`](@ref).
"""
cb_termination_code(p::CBProblem) = @ccall libcb.cb_termination_code(p.handle::CBHandle)::CBTerminationCode

"""
    cb_print_termination_code(p::CBProblem)

Outputs a text version of termination code, see cb_termination_code().
"""
function cb_print_termination_code(p::CBProblem)
    iszero(@ccall libcb.cb_print_termination_code(p.handle::CBHandle)::Cint) || error("Printing termination code failed")
    return
end

"""
    cb_get_objval(p::CBProblem)

Returns the objective value resulting from last descent step (initially undefined).

If no problem modification routines were called since then, it is the objective value at the point returned by
[`cb_get_center`](@ref).
"""
cb_get_objval(p::CBProblem) = @ccall libcb.cb_get_objval(p.handle::CBHandle)::Cdouble

"""
    cb_get_center(p::CBProblem)

Returns the next center point that was produced by the last call to
[`cb_solve!`](@ref cb_solve!(::CBProblem, ::Integer, ::Bool)) (in some problem modification routines the center point may be
updated immediately, in others the center point will be corrected automatically directly before starting the next descent step
and its values may be infeasible till then).

Returns a double vector of length [`cb_get_dim`](@ref cb_get_dim(::CBProblem)). Element `i` will be the value of design
variable ``y_i`` in the next center point (mostly the result of the latest descent step).
Use [`cb_get_center!`](@ref cb_get_center!(::Vector{Float64}, ::CBProblem)) to avoid allocation.
"""
cb_get_center(p::CBProblem) = @inbounds cb_get_center!(Vector{Float64}(undef, cb_get_dim(p)), p)

"""
    cb_get_center!(center::Vector{Float64}, p::CBProblem)

Mutating variant of [`cb_get_center`](@ref cb_get_center(::CBProblem))
"""
Base.@propagate_inbounds function cb_get_center!(center::Vector{Float64}, p::CBProblem)
    @boundscheck length(center) == cb_get_dim(p) || error("Length of center is not appropriate")
    iszero(@ccall libcb.cb_get_center(p.handle::CBHandle, center::Ref{Cdouble})::Cint) || error("Getting center point failed")
    return center
end

"""
    cb_get_sgnorm(p::CBProblem)

Returns Euclidean norm of the latest aggregate subgradient.
"""
cb_get_sgnorm(p::CBProblem) = @ccall libcb.cb_get_sgnorm(p.handle::CBHandle)::Cdouble

"""
    cb_get_subgradient(p::CBProblem)

Returns the latest aggregate subgradient.

Returns a double vector of length [`cb_get_dim`](@ref cb_get_dim(::CBProblem)). Element `i` will be filled with the coordinate
value `i`.
Use [`cb_get_subgradient!`](@ref cb_get_subgradient!(::Vector{Float64}, ::CBProblem)) to avoid allocation.
"""
cb_get_subgradient(p::CBProblem) = @inbounds cb_get_subgradient!(Vector{Float64}(undef, cb_get_dim(p)), p)

"""
    cb_get_subgradient!(subgradient::Vector{Float64}, p::CBProblem)

Mutating variant of [`cb_get_subgradient`](@ref cb_get_subgradient(::CBProblem))
"""
Base.@propagate_inbounds function cb_get_subgradient!(subgradient::Vector{Float64}, p::CBProblem)
    @boundscheck length(subgradient) == cb_get_dim(p) || error("Length of subgradient is not appropriate")
    iszero(@ccall libcb.cb_get_subgradient(p.handle::CBHandle, subgradient::Ref{Cdouble})::Cint) ||
        error("Getting subgradient failed")
    return subgradient
end

"""
    cb_get_candidate_value(p::CBProblem)

Returns the objective value computed in the last step of [`cb_solve!`](@ref cb_solve!(::CBProblem, ::Integer, ::Bool)),
independent of whether this was a descent step or a null step (initially undefined).

If no problem modification routines were called since then, it is the objective value at the point returned by
[`cb_get_candidate`](@ref cb_get_candidate(::CBProblem)). If this last evaluation led to a descent step, then it is the same
value as in [`cb_get_objval`](@ref cb_get_objval(::CBProblem)).
"""
cb_get_candidate_value(p::CBProblem) = @ccall libcb.cb_get_candidate_value(p.handle::CBHandle)::Cdouble

"""
    cb_get_candidate(p::CBProblem)

Returns the last point, the "candidate", at which the function was evaluated in
[`cb_solve!`](@ref cb_solve!(::CBProblem, ::Integer, ::Bool)).

If this evaluation lead to a descent step, it is the same point as in [`cb_get_center`](@ref cb_get_center(::CBProblem)).

Returns a double array of length [`cb_get_dim`](@ref cb_get_dim(::CBProblem)). Element `i` will be the value of design variable
``y_i`` of the point.
Use [`cb_get_candidate!`](@ref cb_get_candidate!(::Vector{Float64}, ::CBProblem)) to avoid allocation.
"""
cb_get_candidate(p::CBProblem) = @inbounds cb_get_candidate!(Vector{Float64}(undef, cb_get_dim(p)), p)

"""
    cb_get_candidate!(candidate::Vector{Float64}, p::CBProblem)

Mutating variant of [`cb_get_candidate`](@ref cb_get_candidate(::CBProblem))
"""
Base.@propagate_inbounds function cb_get_candidate!(candidate::Vector{Float64}, p::CBProblem)
    @boundscheck length(candidate) == cb_get_dim(p) || error("Length of candidate is not appropriate")
    iszero(@ccall libcb.cb_get_candidate(p.handle::CBHandle, candidate::Ref{Cdouble})::Cint) ||
        error("Getting candidate failed")
    return candidate
end

"""
    cb_set_term_relprec!(p::CBProblem, term_relprec::Float64)

Sets the relative precision requirements for successful termination (default `1e-5`).
The algorithm stops with termination code `cbt_relprec_satisfied`, if predicted progress for the next step is less than
`term_relprec` times absolute function value plus one.
"""
function cb_set_term_relprec!(p::CBProblem, term_relprec::Float64)
    iszero(@ccall libcb.cb_set_term_relprec(p.handle::CBHandle, term_relprec::Cdouble)::Cint) ||
        error("Setting relative precision requirements failed")
    return
end

"""
    cb_set_new_center_point!(p::CBProblem, center::Vector{Float64})

Set the starting point/center that will be used in the next call to
[`cb_solve!`](@ref cb_solve!(::CBProblem, ::Integer, ::Bool)). Each call to this routine causes an immediate evaluation of all
oracles.

Arguments:
- `p::CBProblem`

  the problem
- `center::Vector{Float64}` (length [`cb_get_dim`](@ref cb_get_dim(::CBProblem)))

  `center[i]` holds the value of design variable ``y_i``
"""
Base.@propagate_inbounds function cb_set_new_center_point!(p::CBProblem, center::Vector{Float64})
    @boundscheck length(center) == cb_get_dim(p) || error("Length of center point is not appropriate")
    iszero(@ccall libcb.cb_set_new_center_point(p.handle::CBHandle, center::Ref{Cdouble})::Cint) ||
        error("Setting center point failed")
    return
end

"""
    cb_get_function_status(p::CBProblem{T}, function_key::T)

Returns the return value of the latest evaluation call to the function with this `function_key`

Remember, a unique `function_key` must be specified in
[`cb_add_function!`](@ref cb_add_function!(::CBProblem{T}, ::T, ::Function, ::Nothing) where {T}).

Arguments:
- `p::CBProblem{T}`

  the problem
- `function_key::T`

  unique identifier as set in [`cb_add_function!`](@ref cb_add_function!(::CBProblem{T}, ::T, ::Function, ::Nothing) where {T}).

Returns value of latest call to the function having this `function_key`
"""
cb_get_function_status(p::CBProblem{T}, function_key::T) where {T} =
    @ccall libcb.cb_get_function_status(p.handle::CBHandle, p.functions[function_key]::Ptr{Cvoid})::Cint

"""
    cb_get_approximate_slacks(p::CBProblem)

Returns the multipliers for the bound constraints on the design variables; in Lagrangean relaxation they may be interpreted as
primal slacks for inequality constraints.

Returns a double array of length [`cb_get_dim`](@ref cb_get_dim(::CBProblem)). Element `i` will be filled with the coordinate
value `i`
Use [`cb_get_approximate_slacks!`](@ref cb_get_approximate_slacks!(::Vector{Float64}, ::CBProblem)) to avoid allocation.
"""
cb_get_approximate_slacks(p::CBProblem) = @inbounds cb_get_approximate_slacks!(Vector{Float64}(undef, cb_get_dim(p)), p)

"""
    cb_get_approximate_slacks!(slacks::Vector{Float64}, p::CBProblem)

Mutating variant of [`cb_get_approximate_slacks`](@ref cb_get_approximate_slacks(::CBProblem))
"""
Base.@propagate_inbounds function cb_get_approximate_slacks!(slacks::Vector{Float64}, p::CBProblem)
    @boundscheck length(slacks) == cb_get_dim(p) || error("Length of slacks is not appropriate")
    iszero(@ccall libcb.cb_get_approximate_slacks(p.handle::CBHandle, slacks::Ref{Cdouble})::Cint) ||
        error("Getting slacks failed")
    return slacks
end

"""
    cb_get_approximate_primal!(primal::Vector{Float64}, p::CBProblem{T}, function_key::T)

Returns the current approximate primal solution for the function having this `function_key`

The `function_key` must match the one specified in
[`cb_add_function!`](@ref cb_add_function!(::CBProblem{T}, ::T, ::Function, ::Nothing) where {T}). Likewise, the routine is
meaningful only if `primaldim` was set in
  [`cb_add_function!`](@ref cb_add_function!(::CBProblem{T}, ::T, ::Function, ::Nothing) where {T}) and primal vectors were
  returned along with the subgradients in all calls to [`CBFunction`](@ref) with this `function_key`.
In this case it returns the current approximate primal solution aggregated alongside with the aggregate subgradient.
A primal solution may not be available after addition of constraints, if extension of the aggregate subgradient to the new
coordinates failed.

If no primal dimension was set for this function, the routine does nothing.

Arguments:
- `primal::Vector{Float64}` (caller-allocated of length `primaldim`)

  `primal[i]` will be filled with the coordinate value `i`
- `p::CBProblem{T}`

  the problem
- `function_key::T`

  unique identifier as set in [`cb_add_function!`](@ref cb_add_function!(::CBProblem{T}, ::T, ::Function, ::Nothing) where {T})
"""
function cb_get_approximate_primal!(primal::Vector{Float64}, p::CBProblem{T}, function_key::T) where {T}
    iszero(@ccall libcb.cb_get_approximate_primal(
        p.handle::CBHandle,
        p.functions[function_key]::Ptr{Cvoid},
        primal::Ref{Cdouble}
    )::Cint) || error("Getting approximate primal solution failed")
    return primal
end

"""
    cb_get_center_primal!(primal::Vector{Float64}, p::CBProblem{T}, function_key::T)

Returns the best primal solution obtained in the current center point in evaluating the function having this `function_key`

The `function_key` must match the one specified in
[`cb_add_function!`](@ref cb_add_function!(::CBProblem{T}, ::T, ::Function, ::Nothing) where {T}). Likewise, the routine is
meaningful only if `primaldim` was set in
[`cb_add_function!`](@ref cb_add_function!(::CBProblem{T}, ::T, ::Function, ::Nothing) where {T}) and primal vectors were
returned along with the subgradients in all calls to [`CBFunction`](@ref) with this `function_key`.
In this case it returns the current approximate primal solution aggregated alongside with the aggregate subgradient.
A primal solution may not be available after addition of constraints, if extension of the aggregate subgradient to the new
coordinates failed.

If no primal dimension was set for this function, the routine does nothing.

Arguments:
- `primal::Vector{Float64}` (caller-allocated of length `primaldim`)

  `primal[i]` will be filled with the coordinate value `i`
- `p::CBProblem{T}`

  the problem
- `function_key::T`

  unique identifier as set in [`cb_add_function!`](@ref cb_add_function!(::CBProblem{T}, ::T, ::Function, ::Nothing) where {T})
"""
function cb_get_center_primal!(primal::Vector{Float64}, p::CBProblem{T}, function_key::T) where {T}
    iszero(@ccall libcb.cb_get_center_primal(
        p.handle::CBHandle,
        p.functions[function_key]::Ptr{Cvoid},
        primal::Ref{Cdouble}
    )::Cint) || error("Getting center primal solution failed")
    return primal
end

"""
    cb_get_candidate_primal!(primal::Vector{Float64}, p::CBProblem{T}, function_key::T)

Returns the best primal solution returned by the last evaluation of the function having this `function_key` in the point
[`cb_get_candidate`](@ref cb_get_candidate(::CBProblem)).

The `function_key` must match the one specified in
[`cb_add_function!`](@ref cb_add_function!(::CBProblem{T}, ::T, ::Function, ::Nothing) where {T}). Likewise, the routine is
meaningful only if `primaldim` was set in
  [`cb_add_function!`](@ref cb_add_function!(::CBProblem{T}, ::T, ::Function, ::Nothing) where {T}) and primal vectors were
  returned along with the subgradients in all calls to [`CBFunction`](@ref) with this `function_key`.
In this case it returns the current approximate primal solution aggregated alongside with the aggregate subgradient.
A primal solution may not be available after addition of constraints, if extension of the aggregate subgradient to the new
coordinates failed.

If no primal dimension was set for this function, the routine does nothing.

Arguments:
- `primal::Vector{Float64}` (caller-allocated of length `primaldim`)

  `primal[i]` will be filled with the coordinate value `i`
- `p::CBProblem{T}`

  the problem
- `function_key::T`

  unique identifier as set in [`cb_add_function!`](@ref cb_add_function!(::CBProblem{T}, ::T, ::Function, ::Nothing) where {T})
"""
function cb_get_candidate_primal!(primal::Vector{Float64}, p::CBProblem{T}, function_key::T) where {T}
    iszero(@ccall libcb.cb_get_candidate_primal(
        p.handle::CBHandle,
        p.functions[function_key]::Ptr{Cvoid},
        primal::Ref{Cdouble}
    )::Cint) || error("Getting primal candidate failed")
    return primal
end

"""
    cb_set_max_modelsize!(p::CBProblem{T}, function_key::T, modelsize::Integer)

Sets the maximum number of subgradients used in forming the cutting model of the function having this `function_key`

The `function_key` must match the one specified in
[`cb_add_function!`](@ref cb_add_function!(::CBProblem{T}, ::T, ::Function, ::Nothing) where {T}).

Quite often a very small model, e.g., 2, yields very fast iterations and good progress in time (sometimes at the cost of more
evaluations). By limited numerical experience, a significant reduction in the number of evaluations can only be expected if the
bundle is large enough to wrap the function rather tightly. Quite frequently, unfortunately, this entails that solving the
quadratic subproblems is more expensive than function evaluation.

Arguments:
- `p::CBProblem{T}`

  the problem
- `function_key::T`

  unique identifier as set in [`cb_add_function!`](@ref cb_add_function!(::CBProblem{T}, ::T, ::Function, ::Nothing) where {T})
- `modelsize::Integer`

  maximum number of subgradients to be used in forming the cutting model
"""
function cb_set_max_modelsize!(p::CBProblem{T}, function_key::T, modelsize::Integer) where {T}
    iszero(@ccall libcb.cb_set_max_modelsize(
        p.handle::CBHandle,
        p.functions[function_key]::Ptr{Cvoid},
        modelsize::Cint
    )::Cint) || error("Setting maximal model size failed")
    return
end

"""
    cb_set_max_bundlesize!(p::CBProblem{T}, function_key::T, bundlesize::Integer)

Sets the maximum number of subgradients stored for use in forming the model or determining metric information of the function
having this `function_key`. it must be as least as large as `max_modelsize` (and is increased to this if not)

The `function_key` must match the one specified in
[`cb_add_function!`](@ref cb_add_function!(::CBProblem{T}, ::T, ::Function, ::Nothing) where {T}).

Arguments:
- `p::CBProblem{T}`

  the problem
- `function_key::T`

  unique identifier as set in [`cb_add_function!`](@ref cb_add_function!(::CBProblem{T}, ::T, ::Function, ::Nothing) where {T})
- `bundlesize::Integer`

  maximum number of subgradients to be used in forming the cutting model
"""
function cb_set_max_bundlesize!(p::CBProblem{T}, function_key::T, bundlesize::Integer) where {T}
    iszero(@ccall libcb.cb_set_max_bundlesize(
        p.handle::CBHandle,
        p.functions[function_key]::Ptr{Cvoid},
        bundlesize::Cint
    )::Cint) || error("Setting maximal bundle size failed")
    return
end

"""
    cb_set_max_new_subgradients!(p::CBProblem{T}, function_key::T, max_new_subg::Integer)

Sets the maximum number of epsilon subgradients that can be returned in one call to the function having this `function_key`.

The `function_key` must match the one specified in
[`cb_add_function!`](@ref cb_add_function!(::CBProblem{T}, ::T, ::Function, ::Nothing) where {T}).

The parameter `max_new_subg` corresponds directly to the parameter `max_subg` in [`CBFunction`](@ref).

Arguments:
- `p::CBProblem{T}`

  the problem
- `function_key::T`

  unique identifier as set in [`cb_add_function!`](@ref cb_add_function!(::CBProblem{T}, ::T, ::Function, ::Nothing) where {T})
- `max_new_subg::Integer`

  maximum number of subgradients to be used in forming the cutting model
"""
function cb_set_max_new_subgradients!(p::CBProblem{T}, function_key::T, max_new_subg::Integer) where {T}
    iszero(@ccall libcb.cb_set_max_new_subgradients(
        p.handle::CBHandle,
        p.functions[function_key]::Ptr{Cvoid},
        max_new_subg::Cint
    )::Cint) || error("Setting maximal number of subgradients failed")
    return
end

"""
     cb_get_bundle_parameters(p::CBProblem{T}, function_key::T)

Retrieves the two bundle parameters specified in the routines
[`cb_set_max_modelsize!`](@ref cb_set_max_modelsize!(::CBProblem{T}, ::T, ::Integer) where {T}) and
[`cb_set_max_bundlesize!`](@ref cb_set_max_bundlesize!(::CBProblem{T}, ::T, ::Integer) where {T}) for the function having this
`function_key`.

The `function_key` must match the one specified in
[`cb_add_function!`](@ref cb_add_function!(::CBProblem{T}, ::T, ::Function, ::Nothing) where {T}).

Arguments:
- `p::CBProblem{T}`

  the problem
- `function_key::T`

  unique identifier as set in [`cb_add_function!`](@ref cb_add_function!(::CBProblem{T}, ::T, ::Function, ::Nothing) where {T})

Returns two integers:
- the maximum number of subgradients to be used in forming the cutting model
- the maximum number of subgradients stored for use in forming the cutting model
"""
function cb_get_bundle_parameters(p::CBProblem{T}, function_key::T) where {T}
    max_modelsize = Ref{Cint}()
    max_bundlesize = Ref{Cint}()
    iszero(@ccall libcb.cb_get_bundle_parameters(
        p.handle::CBHandle,
        p.functions[function_key]::Ptr{Cvoid},
        max_modelsize::Ref{Cint},
        max_bundlesize::Ref{Cint}
    )::Cint) || error("Getting bundle parameters failed")
    return max_modelsize[], max_bundlesize[]
end

"""
    cb_reinit_function_model!(p::CBProblem{T}, function_key::T)

Clears cutting model, subgradients and stored function values for the function with this `function_key`

This has to be called whenever the specified function was modified so that the old subgradients and/or primal generators are no
longer valid.

The `function_key` must match the one specified in
[`cb_add_function!`](@ref cb_add_function!(::CBProblem{T}, ::T, ::Function, ::Nothing) where {T}).
"""
function cb_reinit_function_model!(p::CBProblem{T}, function_key::T) where {T}
    iszero(@ccall libcb.cb_reinit_function_model(p.handle::CBHandle, p.functions[function_key]::Ptr{Cvoid},)::Cint) ||
        error("Reinitializing function model failed")
    return
end

"""
    cb_get_last_weight(p::CBProblem)

Returns the current weight for the quadratic term in the augmented subproblem (may be interpreted as 1./step_size or
1./trustregion-radius).
"""
cb_get_last_weight(p::CBProblem) = @ccall libcb.cb_get_last_weight(p.handle::CBHandle)::Cdouble

"""
    cb_set_next_weight!(p::CBProblem, weight::Float64)

Sets the  weight (>0) to be used in the quadratic term of the next augmented subproblem (may be interpreted as 1./step_size or
1./trustregion-radius).

Independent of whether the weight violates current min- and max-bounds set in
[`cb_set_min_weight!`](@ref cb_set_min_weight!(::CBProblem, ::Float64)) and
[`cb_set_max_weight!`](@ref cb_set_max_weight!(::CBProblem, ::Float64)), the next model will be computed for this value.
Thereafter, however, it will be updated as usual; in particular, it may be truncated by min and max bounds immediately after
the first subproblem.

In order to guarantee a constant weight (e.g. 1 is frequently a reasonable choice if the automatic default heuristic performs
poorly), set the min and max bounds to the same value, too.
"""
function cb_set_next_weight!(p::CBProblem, weight::Float64)
    iszero(@ccall libcb.cb_set_next_weight(p.handle::CBHandle, weight::Cdouble)::Cint) || error("Setting next weight failed")
    return
end

"""
    cb_set_min_weight!(p::CBProblem, min_weight::Float64)

Sets a lower bound on the weight for the quadratic term of the augmented subproblem.

Nonpositive values indicate no bound.
The new value shows its effect only at first dynamic change of the weight.
"""
function cb_set_min_weight!(p::CBProblem, min_weight::Float64)
    iszero(@ccall libcb.cb_set_min_weight(p.handle::CBHandle, min_weight::Cdouble)::Cint) ||
        error("Setting minimum weight failed")
    return
end

"""
    cb_set_max_weight!(p::CBProblem, max_weight::Float64)

Sets a upper bound on the weight for the quadratic term of the augmented subproblem.

Nonpositive values indicate no bound.
The new value shows its effect only at first dynamic change of the weight.
"""
function cb_set_max_weight!(p::CBProblem, max_weight::Float64)
    iszero(@ccall libcb.cb_set_max_weight(p.handle::CBHandle, max_weight::Cdouble)::Cint) ||
        error("Setting maximum weight failed")
    return
end

@doc """
    CBVariableMetric

- `cbvm_no_scaling`
- `cbvm_diagonal_scaling`
- `cbvm_diagonal_scaling_with_bounds`
"""
@enum CBVariableMetric::Cint cbvm_no_scaling=0 cbvm_diagonal_scaling=1 cbvm_diagonal_scaling_with_bounds=2

"""
    cb_set_variable_metric!(p::CBProblem, do_variable_metric::CBVariableMetric)

Sets a upper bound on the weight for the quadratic term of the augmented subproblem.

Nonpositive values indicate no bound.
The new value shows its effect only at first dynamic change of the weight.
"""
function cb_set_variable_metric!(p::CBProblem, do_variable_metric::CBVariableMetric)
    iszero(@ccall libcb.cb_set_variable_metric(p.handle::CBHandle, do_variable_metric::CBVariableMetric)::Cint) ||
        error("Setting variable metric failed")
    return
end

"""
    cb_get_dim(p::CBProblem)

Returns the current dimension of the design space/argument or -1 if no dimension is set.
"""
cb_get_dim(p::CBProblem) = @ccall libcb.cb_get_dim(p.handle::CBHandle)::Cint

"""
    cb_get_n_functions(p::CBProblem)

Returns the current number of functions in the problem.
"""
cb_get_n_functions(p::CBProblem) = @ccall libcb.cb_get_n_functions(p.handle::CBHandle)::Cint

"""
    cb_get_minus_infinity()

Returns the value "minus infinity", i.e., all bounds <= this value are set to this value and are regarded as minus infinity
Use the constant [`cb_minus_infinity`](@ref) instead of calling this function.
"""
cb_get_minus_infinity() = @ccall libcb.cb_get_minus_infinity()::Cdouble

"""
    const cb_minus_infinity::Float64

Contains the value "minus infinity" as returned by [`cb_get_minus_infinity`](@ref)
"""
const cb_minus_infinity = cb_get_minus_infinity()

"""
    cb_get_plus_infinity()

Returns the value "plus infinity", i.e., all bounds >= this value are set to this value and are regarded as plus infinity
Use the constant [`cb_plus_infinity`](@ref) instead of calling this function.
"""
cb_get_plus_infinity() = @ccall libcb.cb_get_plus_infinity()::Cdouble

"""
    const cb_plus_infinity::Float64

Contains the value "plus infinity" as returned by [`cb_get_plus_infinity`](@ref)
"""
const cb_plus_infinity = cb_get_plus_infinity()

"""
    cb_clear_fail_counts!(p::CBProblem)

clears all fail counts on numerical function oder model failures, may be useful if this caused premature termination.
"""
cb_clear_fail_counts!(p::CBProblem) = @ccall libcb.cb_clear_fail_counts(p.handle::CBHandle)::Cvoid

"""
    cb_set_eval_limit!(p::CBProblem, eval_limit::Integer)

Sets an upper bound on the number of calls to the oracle (use negative numbers for no limit).

If this number is reached, the algorithm will terminate independently of whether the last step was a descent or a null step.
A negative number will be interepreted as no limit.
"""
cb_set_eval_limit!(p::CBProblem, eval_limit::Integer) =
    @ccall libcb.cb_set_eval_limit(p.handle::CBHandle, eval_limit::Cint)::Cvoid

"""
    cb_set_inner_update_limit!(p::CBProblem, update_limit::Integer)

Set an upper bound on the number of inner updates for the cutting model with primal slacks within one null step (use negative
numbers for no limit).

A negative number will be interepreted as no limit, i.e., the updates will be done till a certain precision of the cutting
model is achieved.
"""
cb_set_inner_update_limit!(p::CBProblem, update_limit::Integer) =
    @ccall libcb.cb_set_inner_update_limit(p.handle::CBHandle, update_limit::Cint)::Cvoid

"""
    cb_set_active_bounds_fixing!(p::CBProblem, allow_fixing::Bool)

If set to `true` (the default is `false`), some variables will be fixed automatically to the center value if their bounds are
strongly active (i.e., the corresponding multipliers are big).

The coordinates to be fixed are redetermined in each call following a descent step or a change of the function. An indicator
vector of the variables fixed in the last call can be obtained via the routine
[`cb_get_fixed_active_bounds`](@ref cb_get_fixed_active_bounds(::CBProblem)).

Setting this value to `true` might improve the performance of the algorithm in some instances but there is no convergence
theory. It might be particularly helpful within Lagrangian relaxation if a primal cutting plane approach is used and non-tight
inequalities should be eliminated quickly (fixing then indicates large primal slack values).
"""
cb_set_active_bounds_fixing!(p::CBProblem, allow_fixing::Bool) =
    @ccall libcb.cb_set_active_bounds_fixing(p.handle::CBHandle, allow_fixing::Cint)::Cvoid

"""
    cb_get_fixed_active_bounds(p::CBProblem)

Returns the indicator vector of variables temporarily fixed to the center value due to significantly positive multipliers for
the box constraints, see [`cb_set_active_bounds_fixing!`](@ref cb_set_active_bounds_fixing!(::CBProblem, ::Bool)).

Such a fixing indicates that the corresponding variables would like to stay at their bounds. If no variables were fixed, the
dimension of the vector is zero.

Returns a int array of length [`cb_get_dim`](@ref cb_get_dim(::CBProblem)). Element `i` will be `1` if the variable `i`
was fixed to the bound and `0` otherwise
Use [`cb_get_fixed_active_bounds!`](@ref cb_get_fixed_active_bounds!(::CBProblem, ::Vector{Cint})) to avoid allocation.
"""
cb_get_fixed_active_bounds(p::CBProblem) = @inbounds cb_get_fixed_active_bounds!(p, Vector{Cint}(undef, cb_get_dim(p)))

"""
    cb_get_fixed_active_bounds!(p::CBProblem, indicator::Vector{Cint})

Mutating variant of [`cb_get_fixed_active_bounds`](@ref cb_get_fixed_active_bounds!(::CBProblem, ::Vector{Cint}))
"""
Base.@propagate_inbounds function cb_get_fixed_active_bounds!(p::CBProblem, indicator::Vector{Cint})
    @boundscheck length(indicator) == cb_get_dim(p) || error("Length of indicator is not appropriate")
    iszero(@ccall libcb.cb_get_fixed_active_bounds(p.handle::CBHandle, indicator::Ref{Cint})::Cint) ||
        error("Getting fixed active bounds failed")
    return indicator
end

"""
    cb_set_print_level!(p::CBProblem, pril::Integer)

Specifies the output level (<0 no output at all, =0 errors and warnings, >0 increasingly detailed information)

Output levels:
- `<0` ... no output, not even errors or warnings
- `0` ... no output except for errors and warnings
- `1` ... line summary after each descent step
- `>1` ... undocumented and increasingly detailed log information.
          These higher levels should only be used if requested
          for debugging purposes.

# Example for level 1:

```
00:00:00.00 endit  1   1   1   563.    563.  39041.188  39043.162
00:00:00.00 endit  2   2   2   563.    559.  38488.165  38490.200
00:00:00.00 endit  3   3   3   56.3    555.  33014.533  33211.856
00:00:00.00 endit  4   4   4   5.63    517. -14306.459  2738.0343
00:00:00.00 endit  5   5   5   4.04    148. -2692.1131  2.2150883
00:00:00.00 endit  6   6   6   4.01    1.29  1.7908952  2.0000581
00:00:00.00 endit  7   7   7   3.95  0.0213  1.9999387  2.0000000
00:00:00.00 _endit  8   8   8   3.95 2.94e-05  2.0000000  2.0000000

Column 1      2     3   4   5    6       7       8          9
```
- Column 1: computation time in hh:mm:ss.dd,
- Column 2: "endit" is convenient for grep and stands for "end of iteration".
    Iterations with [`cb_termination_code`](@ref cb_termination_code(::CBProblem)) != 0 are marked with "_endit".
- Column 3: number of descent steps (= calls to [`cb_solve!`](@ref cb_solve!(::CBProblem, ::Integer, ::Bool)))
- Column 4: number of descent and null steps. Up to initialization calls and reevaluations, this is the number of
    evaluation calls to the function oracles from within the bundle method. In the example all calls led to descent steps.
- Column 5: number of innermost iterations. It differs from column 5 only in the case of variables with bounds in which
    case it gives the number of updates of the multipliers for the bounds (or primal slacks in Lagrangean relaxation).
    Exceedingly high numbers in this column indicate that some variables are constantly at their bounds and it might be
    possible to improve convergence by deleting them (i.e. set them as constants to this bound and remove the variable).
- Column 6: the weight of the quadratic term in the augmented problem.
- Column 7: the norm of the aggregate subgradient. If it is small, say below 0.1, then mostly this is good indication that
    the objective value is close to optimal.
- Column 8: the value of the cutting model in the last candidate point. It is always a lower bound on the true function
    value in this point
- Column 9: the objective value in the latest point that led to a descent step, i.e., the point returend by
    [`cb_get_center`](@ref cb_get_center(::CBProblem)). Whenever [`cb_termination_code`](@ref cb_termination_code(::CBProblem))
    returns 0 this is also the objective value of the latest evaluation call to the function oracles and the value in the
    center point of the next iteration.
"""
cb_set_print_level!(p::CBProblem, pril::Integer) = @ccall libcb.cb_set_print_level(p.handle::CBHandle, pril::Cint)::Cvoid

Base.show(io::IO, m::MIME"text/plain", p::CBProblem) =
  print(io, "ConicBundle problem (using C interface) of dimension ", cb_get_dim(p))

include("cppinterface/ConicBundle_cpp.jl")

end