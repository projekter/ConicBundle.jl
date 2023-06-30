@doc raw"""
    cb_clear!(self::CBSumModel)

resets the SumModel to its initial state, in particular it removes but does not delete any Models added in add_model() (their AFTs are deleted if not explicitly denied)
"""
cb_clear!(self::CBSumModel) = @ccall libcb.cb_summodel_clear(self.data::Ptr{Cvoid})::Cvoid

@doc raw"""
    CBSumModel()

constructor
"""
CBSumModel() = CBSumModel(@ccall libcb.cb_summodel_new()::Ptr{Cvoid})

@doc raw"""
    cb_model(self::CBSumModel, fo::Union{<:CBFunctionObject,Nothing})

returns the submodel for FunctionObject fo if it in this model, otherwise 0
"""
cb_model(self::CBSumModel, fo::Union{<:CBFunctionObject,Nothing}) = CBSumBlockModel(@ccall libcb.cb_summodel_model(self.data::Ptr{Cvoid}, (isnothing(fo) ? C_NULL : fo.data)::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_nsubmodels(self::CBSumModel)

returns the number of submodels in this model
"""
cb_nsubmodels(self::CBSumModel) = @ccall libcb.cb_summodel_nsubmodels(self.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_add_model!(self::CBSumModel, model::Union{<:CBSumBlockModel,Nothing})

adds the @a model as submodel to this model (if this model may have submodels); if any error occurs the model is not added
"""
cb_add_model!(self::CBSumModel, model::Union{<:CBSumBlockModel,Nothing}) = @ccall libcb.cb_summodel_add_model(self.data::Ptr{Cvoid}, (isnothing(model) ? C_NULL : model.data)::Ptr{Cvoid})::Cint

@doc raw"""
    cb_remove_model!(self::CBSumModel, fo::Union{<:CBFunctionObject,Nothing})

remove the submodel identified by @a fo from this model, this does NOT destruct the model. It returns the pointer to the model if there is one, otherwise 0
"""
cb_remove_model!(self::CBSumModel, fo::Union{<:CBFunctionObject,Nothing}) = CBSumBlockModel(@ccall libcb.cb_summodel_remove_model(self.data::Ptr{Cvoid}, (isnothing(fo) ? C_NULL : fo.data)::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_eval_function!(self::CBSumModel, y_id::Integer, y::CBMatrix, nullstep_bound::Real, relprec::Real)

see BundleModel::eval_function
"""
function cb_eval_function!(self::CBSumModel, y_id::Integer, y::CBMatrix, nullstep_bound::Real, relprec::Real)
    ub = Ref{Float64}()
    ub_fid = Ref{Int}()
    @ccall libcb.cb_summodel_eval_function(self.data::Ptr{Cvoid}, ub_fid::Ref{Int}, ub::Ref{Float64}, y_id::Cint, y.data::Ptr{Cvoid}, nullstep_bound::Cdouble, relprec::Cdouble)::Cint
    return ub_fid[], ub[]
end

@doc raw"""
    cb_eval_model!(self::CBSumModel, y_id::Integer, y::CBMatrix, relprec::Real)

see BundleModel::eval_model
"""
function cb_eval_model!(self::CBSumModel, y_id::Integer, y::CBMatrix, relprec::Real)
    lb = Ref{Float64}()
    @ccall libcb.cb_summodel_eval_model(self.data::Ptr{Cvoid}, lb::Ref{Float64}, y_id::Cint, y.data::Ptr{Cvoid}, relprec::Cdouble)::Cint
    return lb[]
end

@doc raw"""
    cb_update_model!(self::CBSumModel, model_update::CBModelUpdate, center_id::Integer, center_y::CBMatrix, cand_id::Integer, cand_y::CBMatrix, model_maxviol::Real, H::CBBundleProxObject)

see BundleModel::update_model
"""
cb_update_model!(self::CBSumModel, model_update::CBModelUpdate, center_id::Integer, center_y::CBMatrix, cand_id::Integer, cand_y::CBMatrix, model_maxviol::Real, H::CBBundleProxObject) = @ccall libcb.cb_summodel_update_model(self.data::Ptr{Cvoid}, model_update::Cint, center_id::Cint, center_y.data::Ptr{Cvoid}, cand_id::Cint, cand_y.data::Ptr{Cvoid}, model_maxviol::Cdouble, H.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_synchronize_ids!(self::CBSumModel, new_center_id::Integer, old_center_id::Integer, new_cand_id::Integer, old_cand_id::Integer)

see BundleModel::synchronize_ids
"""
function cb_synchronize_ids!(self::CBSumModel, new_center_id::Integer, old_center_id::Integer, new_cand_id::Integer, old_cand_id::Integer)
    new_aggregate_id = Ref{Int}()
    new_cand_ub_fid = Ref{Int}()
    new_center_ub_fid = Ref{Int}()
    @ccall libcb.cb_summodel_synchronize_ids(self.data::Ptr{Cvoid}, new_center_ub_fid::Ref{Int}, new_center_id::Cint, old_center_id::Cint, new_cand_ub_fid::Ref{Int}, new_cand_id::Cint, old_cand_id::Cint, new_aggregate_id::Ref{Int})::Cint
    return new_center_ub_fid[], new_cand_ub_fid[], new_aggregate_id[]
end

@doc raw"""
    cb_center_modified!(self::CBSumModel, center_id::Integer)

see BundleModel::center_modified
"""
function cb_center_modified!(self::CBSumModel, center_id::Integer)
    function_modification_id = Ref{Int}()
    Bool(@ccall libcb.cb_summodel_center_modified(self.data::Ptr{Cvoid}, function_modification_id::Ref{Int}, center_id::Cint)::Cint)
    return function_modification_id[]
end

@doc raw"""
    cb_recompute_center!(self::CBSumModel, center_id::Integer, y::CBMatrix, accept_only_higher_values::Bool = false, relprec::Real = -1.)

see BundleModel::recompute_center
"""
function cb_recompute_center!(self::CBSumModel, center_id::Integer, y::CBMatrix, accept_only_higher_values::Bool = false, relprec::Real = -1.)
    new_center_ub = Ref{Float64}()
    new_center_ub_fid = Ref{Int}()
    @ccall libcb.cb_summodel_recompute_center(self.data::Ptr{Cvoid}, new_center_ub_fid::Ref{Int}, new_center_ub::Ref{Float64}, center_id::Cint, y.data::Ptr{Cvoid}, accept_only_higher_values::Cint, relprec::Cdouble)::Cint
    return new_center_ub_fid[], new_center_ub[]
end

@doc raw"""
    cb_model_aggregate_modified!(self::CBSumModel, old_model_aggregate_id::Integer)

see BundleModel::model_aggregate_modified
"""
cb_model_aggregate_modified!(self::CBSumModel, old_model_aggregate_id::Integer) = Bool(@ccall libcb.cb_summodel_model_aggregate_modified(self.data::Ptr{Cvoid}, old_model_aggregate_id::Cint)::Cint)

@doc raw"""
    cb_provide_model_aggregate!(self::CBSumModel, y_id::Integer, y::CBMatrix)

see BundleModel::provide_model_aggregate
"""
cb_provide_model_aggregate!(self::CBSumModel, y_id::Integer, y::CBMatrix) = @ccall libcb.cb_summodel_provide_model_aggregate(self.data::Ptr{Cvoid}, y_id::Cint, y.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_add_variable_metric!(self::CBSumModel, H::CBVariableMetric, y_id::Integer, y::CBMatrix, descent_step::Bool, weightu::Real, model_maxviol::Real, indices::Union{<:CBIndexmatrix,Nothing} = nothing)

see BundleModel::add_dynamic_scaling
"""
cb_add_variable_metric!(self::CBSumModel, H::CBVariableMetric, y_id::Integer, y::CBMatrix, descent_step::Bool, weightu::Real, model_maxviol::Real, indices::Union{<:CBIndexmatrix,Nothing} = nothing) = @ccall libcb.cb_summodel_add_variable_metric(self.data::Ptr{Cvoid}, H.data::Ptr{Cvoid}, y_id::Cint, y.data::Ptr{Cvoid}, descent_step::Cint, weightu::Cdouble, model_maxviol::Cdouble, (isnothing(indices) ? C_NULL : indices.data)::Ptr{Cvoid})::Cint

@doc raw"""
    cb_check_center_validity_by_candidate!(self::CBSumModel, cand_minorant_is_below::Bool, center_id::Integer, center_y::CBMatrix)

see BundleModel::check_center_validity_by_candidate
"""
cb_check_center_validity_by_candidate!(self::CBSumModel, cand_minorant_is_below::Bool, center_id::Integer, center_y::CBMatrix) = @ccall libcb.cb_summodel_check_center_validity_by_candidate(self.data::Ptr{Cvoid}, cand_minorant_is_below::Ref{Cint}, center_id::Cint, center_y.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_get_oracle_object!(self::CBSumModel)

as AFTModel has no oracle of its own, this returns the dummy oracle
"""
cb_get_oracle_object!(self::CBSumModel) = CBModifiableOracleObject(@ccall libcb.cb_summodel_get_oracle_object(self.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_make_model_aggregate!(self::CBSumModel, penalty_parameter_increased::Bool, keep_penalty_fixed::Bool)

see SumBlockModel::make_model_aggregate
"""
cb_make_model_aggregate!(self::CBSumModel, penalty_parameter_increased::Bool, keep_penalty_fixed::Bool) = @ccall libcb.cb_summodel_make_model_aggregate(self.data::Ptr{Cvoid}, penalty_parameter_increased::Ref{Cint}, keep_penalty_fixed::Cint)::Cint

@doc raw"""
    cb_lb_function!(self::CBSumModel, y_id::Integer, y::CBMatrix)

see SumBlockModel::lb_function()
"""
cb_lb_function!(self::CBSumModel, y_id::Integer, y::CBMatrix) = @ccall libcb.cb_summodel_lb_function(self.data::Ptr{Cvoid}, y_id::Cint, y.data::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_get_function_minorant!(self::CBSumModel, minorant::CBMinorantPointer, aft::Union{<:CBAffineFunctionTransformation,Nothing} = nothing)

see SumBlockModel::get_function_minorant()
"""
cb_get_function_minorant!(self::CBSumModel, minorant::CBMinorantPointer, aft::Union{<:CBAffineFunctionTransformation,Nothing} = nothing) = @ccall libcb.cb_summodel_get_function_minorant(self.data::Ptr{Cvoid}, minorant.data::Ptr{Cvoid}, (isnothing(aft) ? C_NULL : aft.data)::Ptr{Cvoid})::Cint

@doc raw"""
    cb_get_center_minorant!(self::CBSumModel, minorant::CBMinorantPointer, aft::Union{<:CBAffineFunctionTransformation,Nothing} = nothing)

see SumBlockModel::get_center_minorant()
"""
cb_get_center_minorant!(self::CBSumModel, minorant::CBMinorantPointer, aft::Union{<:CBAffineFunctionTransformation,Nothing} = nothing) = @ccall libcb.cb_summodel_get_center_minorant(self.data::Ptr{Cvoid}, minorant.data::Ptr{Cvoid}, (isnothing(aft) ? C_NULL : aft.data)::Ptr{Cvoid})::Cint

@doc raw"""
    cb_adjust_multiplier!(self::CBSumModel, values_may_have_changed::Bool)

see SumBlockModel::adjust_multiplier()
"""
cb_adjust_multiplier!(self::CBSumModel, values_may_have_changed::Bool) = @ccall libcb.cb_summodel_adjust_multiplier(self.data::Ptr{Cvoid}, values_may_have_changed::Ref{Cint})::Cint

@doc raw"""
    cb_sumbundle_mode!(self::CBSumModel, bh::Union{<:CBSumBundleHandler,Nothing} = nothing, aft::Union{<:CBAffineFunctionTransformation,Nothing} = nothing)

see SumBlockModel::sumbundle_mode
"""
function cb_sumbundle_mode!(self::CBSumModel, bh::Union{<:CBSumBundleHandler,Nothing} = nothing, aft::Union{<:CBAffineFunctionTransformation,Nothing} = nothing)
    mode = Ref{CBMode}()
    @ccall libcb.cb_summodel_sumbundle_mode(self.data::Ptr{Cvoid}, mode::Ref{CBMode}, (isnothing(bh) ? C_NULL : bh.data)::Ptr{Cvoid}, (isnothing(aft) ? C_NULL : aft.data)::Ptr{Cvoid})::Cint
    return mode[]
end

@doc raw"""
    cb_start_sumaugmodel!(self::CBSumModel, blockp::CBQPModelDataPointer, cand_id::Integer, cand_y::CBMatrix, indices::Union{<:CBIndexmatrix,Nothing} = nothing, bh::Union{<:CBSumBundleHandler,Nothing} = nothing, mode::CBMode = cbm_inactive, aft::Union{<:CBAffineFunctionTransformation,Nothing} = nothing)

see BundleModel::start_augmodel() for the first four parameters; for the others see sumbundle_contribution()
"""
cb_start_sumaugmodel!(self::CBSumModel, blockp::CBQPModelDataPointer, cand_id::Integer, cand_y::CBMatrix, indices::Union{<:CBIndexmatrix,Nothing} = nothing, bh::Union{<:CBSumBundleHandler,Nothing} = nothing, mode::CBMode = cbm_inactive, aft::Union{<:CBAffineFunctionTransformation,Nothing} = nothing) = @ccall libcb.cb_summodel_start_sumaugmodel(self.data::Ptr{Cvoid}, blockp.data::Ptr{Cvoid}, cand_id::Cint, cand_y.data::Ptr{Cvoid}, (isnothing(indices) ? C_NULL : indices.data)::Ptr{Cvoid}, (isnothing(bh) ? C_NULL : bh.data)::Ptr{Cvoid}, mode::Cint, (isnothing(aft) ? C_NULL : aft.data)::Ptr{Cvoid})::Cint

@doc raw"""
    cb_get_data!(self::CBSumModel)

see SumBlockModel::get_data()
"""
cb_get_data!(self::CBSumModel) = CBBundleData(@ccall libcb.cb_summodel_get_data(self.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_get_data(self::CBSumModel)

see SumBlockModel::get_data()
"""
cb_get_data(self::CBSumModel) = CBBundleData(@ccall libcb.cb_summodel_get_data2(self.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_set_data!(self::CBSumModel, bd::Union{<:CBBundleData,Nothing})

see SumBlockModel::set_data()
"""
cb_set_data!(self::CBSumModel, bd::Union{<:CBBundleData,Nothing}) = @ccall libcb.cb_summodel_set_data(self.data::Ptr{Cvoid}, (isnothing(bd) ? C_NULL : bd.data)::Ptr{Cvoid})::Cint

@doc raw"""
    cb_get_approximate_primal(self::CBSumModel)

this has no primals, so it returns 1, see SumBlockModel::get_approximate_primal
"""
cb_get_approximate_primal(self::CBSumModel) = CBPrimalData(@ccall libcb.cb_summodel_get_approximate_primal(self.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_get_center_primal(self::CBSumModel)

this has no primals, so it returns 1, see SumBlockModel::get_center_primal
"""
cb_get_center_primal(self::CBSumModel) = CBPrimalData(@ccall libcb.cb_summodel_get_center_primal(self.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_get_candidate_primal(self::CBSumModel)

this has no primals, so it returns 1, see SumBlockModel::get_candidate_primal
"""
cb_get_candidate_primal(self::CBSumModel) = CBPrimalData(@ccall libcb.cb_summodel_get_candidate_primal(self.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_call_primal_extender!(self::CBSumModel, param0::CBPrimalExtender)

this has no primals, so it returns 1, see SumBlockModel::call_primal_extender
"""
cb_call_primal_extender!(self::CBSumModel, param0::CBPrimalExtender) = @ccall libcb.cb_summodel_call_primal_extender(self.data::Ptr{Cvoid}, param0.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_set_bundle_parameters!(self::CBSumModel, bp::CBBundleParameters)

if bp is a SumBundleParametersObject, it is cloned, otherwise the paraemters are copied
"""
cb_set_bundle_parameters!(self::CBSumModel, bp::CBBundleParameters) = @ccall libcb.cb_summodel_set_bundle_parameters(self.data::Ptr{Cvoid}, bp.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_get_bundle_parameters(self::CBSumModel)

see SumBlockModel::get_bundle_parameters()
"""
cb_get_bundle_parameters(self::CBSumModel) = CBBundleParameters(@ccall libcb.cb_summodel_get_bundle_parameters(self.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_clear_model!(self::CBSumModel, discard_minorants_only::Bool = false)

see SumBlockModel::clear_model
"""
cb_clear_model!(self::CBSumModel, discard_minorants_only::Bool = false) = @ccall libcb.cb_summodel_clear_model(self.data::Ptr{Cvoid}, discard_minorants_only::Cint)::Cvoid

@doc raw"""
    cb_clear_aggregates!(self::CBSumModel)

see SumBlockModel::clear_aggregates
"""
cb_clear_aggregates!(self::CBSumModel) = @ccall libcb.cb_summodel_clear_aggregates(self.data::Ptr{Cvoid})::Cvoid

@doc raw"""
    cb_set_out!(self::CBSumModel, pril::Integer = 1)

set output and outputlevel of warnings and errors recursively, see CBout
"""
cb_set_out!(self::CBSumModel, pril::Integer = 1) = @ccall libcb.cb_summodel_set_out(self.data::Ptr{Cvoid}, pril::Cint)::Cvoid

@doc raw"""
    cb_get_preeval_time(self::CBSumModel)

see SumBlockModel::get_preeval_time()
"""
cb_get_preeval_time(self::CBSumModel) = CBCH_Tools::Microseconds(@ccall libcb.cb_summodel_new_get_preeval_time(self.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_get_eval_time(self::CBSumModel)

see SumBlockModel::get_eval_time()
"""
cb_get_eval_time(self::CBSumModel) = CBCH_Tools::Microseconds(@ccall libcb.cb_summodel_new_get_eval_time(self.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_get_posteval_time(self::CBSumModel)

see SumBlockModel::get_posteval_time()
"""
cb_get_posteval_time(self::CBSumModel) = CBCH_Tools::Microseconds(@ccall libcb.cb_summodel_new_get_posteval_time(self.data::Ptr{Cvoid})::Ptr{Cvoid})

