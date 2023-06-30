@doc raw"""
    cb_clear!(self::CBNNCData, start_modification_id::Integer = 0)

reset to initial state (also used by the default constructor)
"""
cb_clear!(self::CBNNCData, start_modification_id::Integer = 0) = @ccall libcb.cb_nncdata_clear(self.data::Ptr{Cvoid}, start_modification_id::Cint)::Cvoid

@doc raw"""
    CBNNCData(fun_factor::Real = 1., fun_task::CBFunctionTask = cbft_objective_function)

initializes BundleData, sets center_primal to NULL and calls clear()
"""
CBNNCData(fun_factor::Real = 1., fun_task::CBFunctionTask = cbft_objective_function) = CBNNCData(@ccall libcb.cb_nncdata_new(fun_factor::Cdouble, fun_task::Cint)::Ptr{Cvoid})

@doc raw"""
    cb_init!(self::CBNNCData, bd::Union{<:CBBundleData,Nothing})

if @a bd is of type NNCData, initialize to this data
"""
cb_init!(self::CBNNCData, bd::Union{<:CBBundleData,Nothing}) = @ccall libcb.cb_nncdata_init(self.data::Ptr{Cvoid}, (isnothing(bd) ? C_NULL : bd.data)::Ptr{Cvoid})::Cint

@doc raw"""
    cb_clone(self::CBNNCData)

return a pointer to a clone of this
"""
cb_clone(self::CBNNCData) = CBBundleData(@ccall libcb.cb_nncdata_clone(self.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_do_step!(self::CBNNCData, point_id::Integer)

if the candidate information is available and consitent for point_id, copy it from cand to center and return 0, otherwise return 1
"""
cb_do_step!(self::CBNNCData, point_id::Integer) = @ccall libcb.cb_nncdata_do_step(self.data::Ptr{Cvoid}, point_id::Cint)::Cint

@doc raw"""
    cb_synchronize_ids!(self::CBNNCData, new_center_id::Integer, old_center_id::Integer, new_cand_id::Integer, old_cand_id::Integer, new_prex_id::Integer = 0)

synchronize ids in any case by discarding inconsistent parts but return number of errors
"""
function cb_synchronize_ids!(self::CBNNCData, new_center_id::Integer, old_center_id::Integer, new_cand_id::Integer, old_cand_id::Integer, new_prex_id::Integer = 0)
    new_aggregate_id = Ref{Int}()
    new_cand_ub_fid = Ref{Int}()
    new_center_ub_fid = Ref{Int}()
    @ccall libcb.cb_nncdata_synchronize_ids(self.data::Ptr{Cvoid}, new_center_ub_fid::Ref{Int}, new_center_id::Cint, old_center_id::Cint, new_cand_ub_fid::Ref{Int}, new_cand_id::Cint, old_cand_id::Cint, new_aggregate_id::Ref{Int}, new_prex_id::Cint)::Cint
    return new_center_ub_fid[], new_cand_ub_fid[], new_aggregate_id[]
end

@doc raw"""
    cb_clear_model!(self::CBNNCData, discard_minorants_only::Bool = false)

clear the cutting model and all function evaluations
"""
cb_clear_model!(self::CBNNCData, discard_minorants_only::Bool = false) = @ccall libcb.cb_nncdata_clear_model(self.data::Ptr{Cvoid}, discard_minorants_only::Cint)::Cvoid

@doc raw"""
    cb_clear_aggregates!(self::CBNNCData)

remove all aggregate minorants
"""
cb_clear_aggregates!(self::CBNNCData) = @ccall libcb.cb_nncdata_clear_aggregates(self.data::Ptr{Cvoid})::Cvoid

@doc raw"""
    cb_call_primal_extender!(self::CBNNCData, prex::CBPrimalExtender, include_candidates::Bool = true)

see the last argument of FunctionOracle::evaluate()
"""
cb_call_primal_extender!(self::CBNNCData, prex::CBPrimalExtender, include_candidates::Bool = true) = @ccall libcb.cb_nncdata_call_primal_extender(self.data::Ptr{Cvoid}, prex.data::Ptr{Cvoid}, include_candidates::Cint)::Cint

@doc raw"""
    cb_apply_modification!(self::CBNNCData, param0::CBGroundsetModification, mex::Union{<:CBMinorantExtender,Nothing})

rearrange/extend the minorants according to the given groundset modifications
"""
cb_apply_modification!(self::CBNNCData, param0::CBGroundsetModification, mex::Union{<:CBMinorantExtender,Nothing}) = @ccall libcb.cb_nncdata_apply_modification(self.data::Ptr{Cvoid}, param0.data::Ptr{Cvoid}, (isnothing(mex) ? C_NULL : mex.data)::Ptr{Cvoid})::Cint

@doc raw"""
    cb_get_approximate_primal(self::CBNNCData)

return the PrimalData corresponding to the aggregate
"""
cb_get_approximate_primal(self::CBNNCData) = CBPrimalData(@ccall libcb.cb_nncdata_get_approximate_primal(self.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_get_center_primal(self::CBNNCData)

return the PrimalData corresponding to the aggregate
"""
cb_get_center_primal(self::CBNNCData) = CBPrimalData(@ccall libcb.cb_nncdata_get_center_primal(self.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_get_candidate_primal(self::CBNNCData)

return the PrimalData delivered by the last call of FunctionOracle::evaluate()
"""
cb_get_candidate_primal(self::CBNNCData) = CBPrimalData(@ccall libcb.cb_nncdata_get_candidate_primal(self.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_get_model_data(self::CBNNCData, model_minorants::CBMinorantBundle, model_coeff::CBMatrix)

the minorants currently used in the model; the list may be empty or max contain other minorants than returned in get_latest_minorants(); the minorants still need to be mutliplied by function_factor
"""
cb_get_model_data(self::CBNNCData, model_minorants::CBMinorantBundle, model_coeff::CBMatrix) = @ccall libcb.cb_nncdata_get_model_data(self.data::Ptr{Cvoid}, model_minorants.data::Ptr{Cvoid}, model_coeff.data::Ptr{Cvoid})::Cint

