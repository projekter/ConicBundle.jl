@doc raw"""
    cb_get_sumbundle(self::CBSumBundleHandler)

returns the sumbundle
"""
cb_get_sumbundle(self::CBSumBundleHandler) = CBSumBundle(@ccall libcb.cb_sumbundlehandler_get_sumbundle(self.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_handles!(self::CBSumBundleHandler, ft::CBFunctionTask)

returns true if this FunctionTask is handled
"""
cb_handles!(self::CBSumBundleHandler, ft::CBFunctionTask) = Bool(@ccall libcb.cb_sumbundlehandler_handles(self.data::Ptr{Cvoid}, ft::Cint)::Cint)

@doc raw"""
    cb_get_new_index(self::CBSumBundleHandler, ft::CBFunctionTask)

returns the index of the newest subgradient in the bundle
"""
cb_get_new_index(self::CBSumBundleHandler, ft::CBFunctionTask) = @ccall libcb.cb_sumbundlehandler_get_new_index(self.data::Ptr{Cvoid}, ft::Cint)::Cint

@doc raw"""
    cb_initialization_needed(self::CBSumBundleHandler, ft::CBFunctionTask)

returns true if the corresponding part is root with contributions but has bundle_size 0
"""
cb_initialization_needed(self::CBSumBundleHandler, ft::CBFunctionTask) = Bool(@ccall libcb.cb_sumbundlehandler_initialization_needed(self.data::Ptr{Cvoid}, ft::Cint)::Cint)

@doc raw"""
    cb_initialization_needed(self::CBSumBundleHandler)

returns true if one of the parts is root with contributions but has bundle_size 0
"""
cb_initialization_needed(self::CBSumBundleHandler) = Bool(@ccall libcb.cb_sumbundlehandler_initialization_needed2(self.data::Ptr{Cvoid})::Cint)

@doc raw"""
    cb_set_parent_information!(self::CBSumBundleHandler, parent_sbh::Union{<:CBSumBundleHandler,Nothing}, aft::Union{<:CBAffineFunctionTransformation,Nothing}, in_mode::CBMode)

* @brief sets a parent handler, an aft and prepares for the next in_mode. If the respective pointers are null and in_mode is active, a part in child mode is changed to root.

    If the parent pointer is not null, the routine calls align_bundle in
    order to match the existing bundle to the parent's.  Therefore, if
    the parent pointer is not null and the parent has a positive bundle
    size, set_parent_information may only be called if all parts that
    have contributors also have a bundle containing an aggregate,
    i.e. they need a bundlesize of at least one.

    If in_mode==root, the handler has to perpare the sumbundle to serve as root.

    If in_mode==inactive, the handler has to remove any relevant
    contributions and to deactivate the sumbundle

    If in_mode==child, the action depends on the the current mode. If
    the mode is child or root already, no changes are required. In
    particular, if the mode is root and the in_mode child is desired,
    this needs to be achieved later by a call to add_contributions()
    (once all bundle components have been collected). If the current
    mode is inactive, the sumbundle's mode is changed to root and again
    an add_contributions() will be needed to change that to child.

    
"""
cb_set_parent_information!(self::CBSumBundleHandler, parent_sbh::Union{<:CBSumBundleHandler,Nothing}, aft::Union{<:CBAffineFunctionTransformation,Nothing}, in_mode::CBMode) = @ccall libcb.cb_sumbundlehandler_set_parent_information(self.data::Ptr{Cvoid}, (isnothing(parent_sbh) ? C_NULL : parent_sbh.data)::Ptr{Cvoid}, (isnothing(aft) ? C_NULL : aft.data)::Ptr{Cvoid}, in_mode::Cint)::Cint

@doc raw"""
    cb_reset_function_factor!(self::CBSumBundleHandler, ft::CBFunctionTask, factor::Real)

resets the value of the function factor for this part of sumbundle
"""
cb_reset_function_factor!(self::CBSumBundleHandler, ft::CBFunctionTask, factor::Real) = @ccall libcb.cb_sumbundlehandler_reset_function_factor(self.data::Ptr{Cvoid}, ft::Cint, factor::Cdouble)::Cint

@doc raw"""
    cb_set_bundle_parameters!(self::CBSumBundleHandler, bp::CBBundleParameters)

sets max_bundle_size and max_model_size for all parts; this may be increased internally if the mode and/or the parent handler require this
"""
cb_set_bundle_parameters!(self::CBSumBundleHandler, bp::CBBundleParameters) = @ccall libcb.cb_sumbundlehandler_set_bundle_parameters(self.data::Ptr{Cvoid}, bp.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_get_increase_factor(self::CBSumBundleHandler)

returns the increase factor for the unbounded part of sumbundle, otherwise 1.
"""
cb_get_increase_factor(self::CBSumBundleHandler) = @ccall libcb.cb_sumbundlehandler_get_increase_factor(self.data::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_update_model!(self::CBSumBundleHandler, model_update::CBModelUpdate, center_id::Integer, center_y::CBMatrix, cand_id::Integer, cand_y::CBMatrix, model_maxviol::Real, H::CBBundleProxObject)

updates the sumbundle maybe even if not active; in this it tries to do the same as bh, storing the aggregate at the same place and providing new room at the same place; the new subgradient is *not* entered here but in set_cand_minorant()
"""
cb_update_model!(self::CBSumBundleHandler, model_update::CBModelUpdate, center_id::Integer, center_y::CBMatrix, cand_id::Integer, cand_y::CBMatrix, model_maxviol::Real, H::CBBundleProxObject) = @ccall libcb.cb_sumbundlehandler_update_model(self.data::Ptr{Cvoid}, model_update.data::Ptr{Cvoid}, center_id::Cint, center_y.data::Ptr{Cvoid}, cand_id::Cint, cand_y.data::Ptr{Cvoid}, model_maxviol::Cdouble, H.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_eval_model(self::CBSumBundleHandler, yid::Integer, y::CBMatrix)

evaluate the model value
"""
function cb_eval_model(self::CBSumBundleHandler, yid::Integer, y::CBMatrix)
    lb = Ref{Float64}()
    @ccall libcb.cb_sumbundlehandler_eval_model(self.data::Ptr{Cvoid}, lb::Ref{Float64}, yid::Cint, y.data::Ptr{Cvoid})::Cint
    return lb[]
end

@doc raw"""
    cb_lb_model(self::CBSumBundleHandler, yid::Integer, y::CBMatrix)

returns a *quick* lower bound for the model value
"""
cb_lb_model(self::CBSumBundleHandler, yid::Integer, y::CBMatrix) = @ccall libcb.cb_sumbundlehandler_lb_model(self.data::Ptr{Cvoid}, yid::Cint, y.data::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_normalize_sumbundle!(self::CBSumBundleHandler)

* brief bring the bundle into a normal form so that contributions may be added and subtracted consistently

    The normal form has the aggregate w.r.t. to the last bundle cofficients
    stored in column aggr->index in scaled form so that in coeff all weight is
    set to it. No other columns are modified.

    When add_contributions is called, it is assumed that parent bundle and
    contributing bundle are in this form. Some care has to be taken that this
    normalization is happening for parent and children in coordinated form so
    that remove_contribution() does not cause havoc. For this, the sumbundle
    should always be normalized in (or right before) calling
    SumBlockModel::sumbundle_contribution() before starting any interaction with
    the parent or the children.
    
"""
cb_normalize_sumbundle!(self::CBSumBundleHandler) = @ccall libcb.cb_sumbundlehandler_normalize_sumbundle(self.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_contribute_new_minorants!(self::CBSumBundleHandler)

once all new minorants have been collected at this level, they are passed to the next if required
"""
cb_contribute_new_minorants!(self::CBSumBundleHandler) = @ccall libcb.cb_sumbundlehandler_contribute_new_minorants(self.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_remove_contributions!(self::CBSumBundleHandler)

remove own contributions to the parent and set the states correspondingly
"""
cb_remove_contributions!(self::CBSumBundleHandler) = @ccall libcb.cb_sumbundlehandler_remove_contributions(self.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_add_contributions!(self::CBSumBundleHandler)

* @brief add root parts of this sumbundle to the parent's sumbundle and set the states correspondingly

        If the current mode is child, nothing is done, because it is
        assumed that the sumbundle is already part of the parent. Thus,
        only root parts are considerd new and are added. This is only
        done, if the parenhandler accepts these kind of parts.

        When a contribution is added to a parent that is currently
        inactive or child, any respective parent's contributions to the
        parent's parent are first removed and then the mode of the
        parent is set to root. Thus, the parent has to call
        add_contributions afterwards, if it still wants to contribute to
        its own parent.

        add_contribution should only be called from
        SumBlockModel::sumbundle_contribution() with a normalized bundle,
        i.e., normalize_sumbundle() should have been called before.
     
"""
cb_add_contributions!(self::CBSumBundleHandler) = @ccall libcb.cb_sumbundlehandler_add_contributions(self.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_start_augmodel!(self::CBSumBundleHandler, bp::CBQPModelDataPointer, cand_id::Integer, cand_y::CBMatrix, indices::Union{<:CBIndexmatrix,Nothing}, ft::CBFunctionTask)

start the augmented model block for FunctionTask ft if to be handled here and increase xdim correspondingly; if there are several, use start_augmodel(QP_SumBlock&,Integer&) instead
"""
cb_start_augmodel!(self::CBSumBundleHandler, bp::CBQPModelDataPointer, cand_id::Integer, cand_y::CBMatrix, indices::Union{<:CBIndexmatrix,Nothing}, ft::CBFunctionTask) = @ccall libcb.cb_sumbundlehandler_start_augmodel(self.data::Ptr{Cvoid}, bp.data::Ptr{Cvoid}, cand_id::Cint, cand_y.data::Ptr{Cvoid}, (isnothing(indices) ? C_NULL : indices.data)::Ptr{Cvoid}, ft::Cint)::Cint

@doc raw"""
    cb_start_augmodel!(self::CBSumBundleHandler, bp::CBQPModelDataPointer, sumblock::CBQPSumModelDataObject, cand_id::Integer, cand_y::CBMatrix, indices::Union{<:CBIndexmatrix,Nothing} = nothing)

add augmented model blocks to the sumblock for parts to be handled here and increase xdim correspondingly
"""
cb_start_augmodel!(self::CBSumBundleHandler, bp::CBQPModelDataPointer, sumblock::CBQPSumModelDataObject, cand_id::Integer, cand_y::CBMatrix, indices::Union{<:CBIndexmatrix,Nothing} = nothing) = @ccall libcb.cb_sumbundlehandler_start_augmodel2(self.data::Ptr{Cvoid}, bp.data::Ptr{Cvoid}, sumblock.data::Ptr{Cvoid}, cand_id::Cint, cand_y.data::Ptr{Cvoid}, (isnothing(indices) ? C_NULL : indices.data)::Ptr{Cvoid})::Cint

@doc raw"""
    cb_make_model_aggregate!(self::CBSumBundleHandler, increased::Bool, fixed::Bool)

see SumBlockModel::make_model_aggregate
"""
cb_make_model_aggregate!(self::CBSumBundleHandler, increased::Bool, fixed::Bool) = @ccall libcb.cb_sumbundlehandler_make_model_aggregate(self.data::Ptr{Cvoid}, increased::Ref{Cint}, fixed::Cint)::Cint

@doc raw"""
    cb_provide_model_aggregate!(self::CBSumBundleHandler)

see SumBlockModel::provide_model_aggregate
"""
cb_provide_model_aggregate!(self::CBSumBundleHandler) = @ccall libcb.cb_sumbundlehandler_provide_model_aggregate(self.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_adjust_multiplier!(self::CBSumBundleHandler, values_may_have_changed::Bool)

see SumBlockModel::adjust_multiplier
"""
cb_adjust_multiplier!(self::CBSumBundleHandler, values_may_have_changed::Bool) = @ccall libcb.cb_sumbundlehandler_adjust_multiplier(self.data::Ptr{Cvoid}, values_may_have_changed::Ref{Cint})::Cint

@doc raw"""
    cb_contribute_initial_bundle!(self::CBSumBundleHandler, ft::CBFunctionTask, bundle_minorants::CBMinorantBundle, coeff::CBMatrix)

* @brief (re)initialize the bundle of the respective part

       This is only possible, if the handler handles a bundle for this.
       If not this causes an error. The main purpose is really
       to start the bundle on the fly if a parent starts a sumbundle
       in the middle of the computation. The information in coeff
       is assumed to yield the current aggregate.

       If this bundle part has been contributed berfore (its mode is child), the
       contribution is first removed from the parent. Then any existing
       information except for the function factor is discarded and replaced by
       the new bundle information with the number of contributions set to 1.

       If the size of primals is not zero, it must have one nonzero entry per
       column of minorants and all future calls via set_cand_minorant() also have
       to provide exactly one primal for each update.

    
"""
cb_contribute_initial_bundle!(self::CBSumBundleHandler, ft::CBFunctionTask, bundle_minorants::CBMinorantBundle, coeff::CBMatrix) = @ccall libcb.cb_sumbundlehandler_contribute_initial_bundle(self.data::Ptr{Cvoid}, ft::Cint, bundle_minorants.data::Ptr{Cvoid}, coeff.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_install_external_aggregate!(self::CBSumBundleHandler, ft::CBFunctionTask, aggr::CBMinorantPointer, aggr_coeff::Real)

* @brief replace the aggregate by one from outside

       This is only possible, if the handler handles a bundle for this and
       the bundle is not active. The main purpose is to switch from an
       external model to a newly contributing sumbundle that has been
       updated all along but not been in use. Installing the external
       aggregate then ensures convergence.

       If primal is not zero, it must already have one nonzero entry per
       column of existing minorants and all future calls via set_cand_minorant()
       also have to provide exactly one primal for each update.

    
"""
cb_install_external_aggregate!(self::CBSumBundleHandler, ft::CBFunctionTask, aggr::CBMinorantPointer, aggr_coeff::Real) = @ccall libcb.cb_sumbundlehandler_install_external_aggregate(self.data::Ptr{Cvoid}, ft::Cint, aggr.data::Ptr{Cvoid}, aggr_coeff::Cdouble)::Cint

@doc raw"""
    cb_set_cand_minorant!(self::CBSumBundleHandler, ft::CBFunctionTask, minorant::CBMinorantPointer)

* @brief set the new minorant information of the candidate

       Only one minorant can take this position.
       This is mainly due to that only one primal information can be associated
       with each minorant. Thus adding must involve compatible primals and
       this is easiest to guarantee if there is only one inital type at the
       lowest level. Contributions to parents will not care about the primals.

    
"""
cb_set_cand_minorant!(self::CBSumBundleHandler, ft::CBFunctionTask, minorant::CBMinorantPointer) = @ccall libcb.cb_sumbundlehandler_set_cand_minorant(self.data::Ptr{Cvoid}, ft::Cint, minorant.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_clear_model!(self::CBSumBundleHandler)

first calls remove_contributions(), then discards all current models
"""
cb_clear_model!(self::CBSumBundleHandler) = @ccall libcb.cb_sumbundlehandler_clear_model(self.data::Ptr{Cvoid})::Cvoid

@doc raw"""
    cb_clear_aggregates!(self::CBSumBundleHandler)

first calls remove_contributions(), then discards all aggregate minorants in the current model
"""
cb_clear_aggregates!(self::CBSumBundleHandler) = @ccall libcb.cb_sumbundlehandler_clear_aggregates(self.data::Ptr{Cvoid})::Cvoid

@doc raw"""
    cb_clear_cand_minorants!(self::CBSumBundleHandler)

before contributing the new evaluation results, the candidates have to be cleared
"""
cb_clear_cand_minorants!(self::CBSumBundleHandler) = @ccall libcb.cb_sumbundlehandler_clear_cand_minorants(self.data::Ptr{Cvoid})::Cvoid

@doc raw"""
    cb_add_variable_metric!(self::CBSumBundleHandler, ft::CBFunctionTask, H::CBVariableMetric, yid::Integer, y::CBMatrix, descent_step::Bool, weightu::Real, model_maxviol::Real, indices::Union{<:CBIndexmatrix,Nothing} = nothing)

see DynamicScaling
"""
cb_add_variable_metric!(self::CBSumBundleHandler, ft::CBFunctionTask, H::CBVariableMetric, yid::Integer, y::CBMatrix, descent_step::Bool, weightu::Real, model_maxviol::Real, indices::Union{<:CBIndexmatrix,Nothing} = nothing) = @ccall libcb.cb_sumbundlehandler_add_variable_metric(self.data::Ptr{Cvoid}, ft::Cint, H.data::Ptr{Cvoid}, yid::Cint, y.data::Ptr{Cvoid}, descent_step::Cint, weightu::Cdouble, model_maxviol::Cdouble, (isnothing(indices) ? C_NULL : indices.data)::Ptr{Cvoid})::Cint

@doc raw"""
    cb_add_variable_metric!(self::CBSumBundleHandler, H::CBVariableMetric, yid::Integer, center_y::CBMatrix, descent_step::Bool, weightu::Real, model_maxviol::Real, indices::Union{<:CBIndexmatrix,Nothing} = nothing)

see DynamicScaling
"""
cb_add_variable_metric!(self::CBSumBundleHandler, H::CBVariableMetric, yid::Integer, center_y::CBMatrix, descent_step::Bool, weightu::Real, model_maxviol::Real, indices::Union{<:CBIndexmatrix,Nothing} = nothing) = @ccall libcb.cb_sumbundlehandler_add_variable_metric2(self.data::Ptr{Cvoid}, H.data::Ptr{Cvoid}, yid::Cint, center_y.data::Ptr{Cvoid}, descent_step::Cint, weightu::Cdouble, model_maxviol::Cdouble, (isnothing(indices) ? C_NULL : indices.data)::Ptr{Cvoid})::Cint

@doc raw"""
    cb_guess_curvature(self::CBSumBundleHandler, mnrts::CBMinorantBundle, selected_indices::CBIndexmatrix, cand_id::Integer, cand_y::CBMatrix, model_maxviol::Real)

computes an estimate of the current curvature to support SumModel in the selection of submodles
"""
cb_guess_curvature(self::CBSumBundleHandler, mnrts::CBMinorantBundle, selected_indices::CBIndexmatrix, cand_id::Integer, cand_y::CBMatrix, model_maxviol::Real) = @ccall libcb.cb_sumbundlehandler_guess_curvature(self.data::Ptr{Cvoid}, mnrts.data::Ptr{Cvoid}, selected_indices.data::Ptr{Cvoid}, cand_id::Cint, cand_y.data::Ptr{Cvoid}, model_maxviol::Cdouble)::Cdouble

