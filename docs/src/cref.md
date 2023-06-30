```@meta
CurrentModule = ConicBundle
```
# Reference of C interface
The C interface is manually ported to Julia. In particular, the comments were curated and some logical coding choices were
made. This is a good interface.
Some of the functions will also be available in the C++ interface; but do not mix them with objects created by the C interface.

```@docs
CBProblem
cb_clear!(::CBProblem)
cb_set_default!(::CBProblem)
cb_init_problem!(::CBProblem, ::Integer)
CBFunction
CBSubgExt
CBFunctionTask
cb_add_function!(::CBProblem{T}, ::T, ::Function, ::Nothing) where {T}
cb_set_lower_bound!(::CBProblem, ::Integer, ::Float64)
cb_set_upper_bound!(::CBProblem, ::Integer, ::Float64)
cb_append_variables!(::CBProblem, ::Integer)
cb_delete_variables!(::CBProblem, ::Vector{Cint})
cb_reassign_variables!(::CBProblem, ::Vector{Cint})
cb_solve!(::CBProblem, ::Integer, ::Bool)
CBTerminationCode
cb_termination_code(::CBProblem)
cb_print_termination_code(::CBProblem)
cb_get_objval(::CBProblem)
cb_get_center(::CBProblem)
cb_get_center!(::Vector{Float64}, ::CBProblem)
cb_get_sgnorm(::CBProblem)
cb_get_subgradient(::CBProblem)
cb_get_subgradient!(::Vector{Float64}, ::CBProblem)
cb_get_candidate_value(::CBProblem)
cb_get_candidate(::CBProblem)
cb_get_candidate!(::Vector{Float64}, ::CBProblem)
cb_set_term_relprec!(::CBProblem, ::Float64)
cb_set_new_center_point!(::CBProblem, ::Vector{Float64})
cb_get_function_status(::CBProblem{T}, ::T) where {T}
cb_get_approximate_slacks(::CBProblem)
cb_get_approximate_slacks!(::Vector{Float64}, ::CBProblem)
cb_get_approximate_primal!(::Vector{Float64}, ::CBProblem{T}, ::T) where {T}
cb_get_center_primal!(::Vector{Float64}, ::CBProblem{T}, ::T) where {T}
cb_get_candidate_primal!(::Vector{Float64}, ::CBProblem{T}, ::T) where {T}
cb_set_max_modelsize!(::CBProblem{T}, ::T, ::Integer) where {T}
cb_set_max_bundlesize!(::CBProblem{T}, ::T, ::Integer) where {T}
cb_set_max_new_subgradients!(::CBProblem{T}, ::T, ::Integer) where {T}
cb_get_bundle_parameters(::CBProblem{T}, ::T) where {T}
cb_reinit_function_model!(::CBProblem{T}, ::T) where {T}
cb_get_last_weight(::CBProblem)
cb_set_next_weight!(::CBProblem, ::Float64)
cb_set_min_weight!(::CBProblem, ::Float64)
cb_set_max_weight!(::CBProblem, ::Float64)
CBVariableMetric
cb_set_variable_metric!(::CBProblem, ::CBVariableMetric)
cb_get_dim(::CBProblem)
cb_get_n_functions(::CBProblem)
cb_get_minus_infinity
cb_get_plus_infinity
cb_clear_fail_counts!(::CBProblem)
cb_set_eval_limit!(::CBProblem, ::Integer)
cb_set_inner_update_limit!(::CBProblem, ::Integer)
cb_set_active_bounds_fixing!(::CBProblem, ::Bool)
cb_get_fixed_active_bounds(::CBProblem)
cb_get_fixed_active_bounds!(::CBProblem, ::Vector{Cint})
cb_set_print_level!(::CBProblem, ::Integer)
cb_minus_infinity
cb_plus_infinity
```
