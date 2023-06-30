@doc raw"""
    cb_clear!(self::CBSOCData, start_modification_id::Integer = 0)

reset to initial state (also used by the default constructor)
"""
cb_clear!(self::CBSOCData, start_modification_id::Integer = 0) = @ccall libcb.cb_socdata_clear(self.data::Ptr{Cvoid}, start_modification_id::Cint)::Cvoid

@doc raw"""
    CBSOCData(fun_factor::Real = 1., fun_task::CBFunctionTask = cbft_objective_function)

calls clear()
"""
CBSOCData(fun_factor::Real = 1., fun_task::CBFunctionTask = cbft_objective_function) = CBSOCData(@ccall libcb.cb_socdata_new(fun_factor::Cdouble, fun_task::Cint)::Ptr{Cvoid})

@doc raw"""
    cb_init!(self::CBSOCData, bd::Union{<:CBBundleData,Nothing})

if @a bd is of type SOCData, initialize to this data
"""
cb_init!(self::CBSOCData, bd::Union{<:CBBundleData,Nothing}) = @ccall libcb.cb_socdata_init(self.data::Ptr{Cvoid}, (isnothing(bd) ? C_NULL : bd.data)::Ptr{Cvoid})::Cint

@doc raw"""
    cb_clone(self::CBSOCData)

return a pointer to a clone of this
"""
cb_clone(self::CBSOCData) = CBBundleData(@ccall libcb.cb_socdata_clone(self.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_do_step!(self::CBSOCData, point_id::Integer)

if the candidate information is available and consitent for point_id, copy it from cand to center and return 0, otherwise return 1
"""
cb_do_step!(self::CBSOCData, point_id::Integer) = @ccall libcb.cb_socdata_do_step(self.data::Ptr{Cvoid}, point_id::Cint)::Cint

@doc raw"""
    cb_store_SOCvec!(self::CBSOCData, SOCvec::CBMatrix)

if max_old_minorants > 0, it adds the SOCvec cyclically to SOCvecs keeping max_old_minorants of them
"""
cb_store_SOCvec!(self::CBSOCData, SOCvec::CBMatrix) = @ccall libcb.cb_socdata_store_socvec(self.data::Ptr{Cvoid}, SOCvec.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_form_bundlevecs!(self::CBSOCData, max_columns::Integer)

starting with aggregate and cand_SOCvec add further ones as needed
"""
cb_form_bundlevecs!(self::CBSOCData, max_columns::Integer) = @ccall libcb.cb_socdata_form_bundlevecs(self.data::Ptr{Cvoid}, max_columns::Cint)::Cint

@doc raw"""
    cb_synchronize_ids!(self::CBSOCData, new_center_id::Integer, old_center_id::Integer, new_cand_id::Integer, old_cand_id::Integer, new_prex_id::Integer = 0)

synchronize ids in any case by discarding inconsistent parts but return number of errors
"""
function cb_synchronize_ids!(self::CBSOCData, new_center_id::Integer, old_center_id::Integer, new_cand_id::Integer, old_cand_id::Integer, new_prex_id::Integer = 0)
    new_aggregate_id = Ref{Int}()
    new_cand_ub_fid = Ref{Int}()
    new_center_ub_fid = Ref{Int}()
    @ccall libcb.cb_socdata_synchronize_ids(self.data::Ptr{Cvoid}, new_center_ub_fid::Ref{Int}, new_center_id::Cint, old_center_id::Cint, new_cand_ub_fid::Ref{Int}, new_cand_id::Cint, old_cand_id::Cint, new_aggregate_id::Ref{Int}, new_prex_id::Cint)::Cint
    return new_center_ub_fid[], new_cand_ub_fid[], new_aggregate_id[]
end

@doc raw"""
    cb_clear_model!(self::CBSOCData, discard_minorants_only::Bool = false)

clear the cutting model and all function evaluations
"""
cb_clear_model!(self::CBSOCData, discard_minorants_only::Bool = false) = @ccall libcb.cb_socdata_clear_model(self.data::Ptr{Cvoid}, discard_minorants_only::Cint)::Cvoid

@doc raw"""
    cb_clear_aggregates!(self::CBSOCData)

delete all kinds of aggregates but keep explicit parts of the cutting model
"""
cb_clear_aggregates!(self::CBSOCData) = @ccall libcb.cb_socdata_clear_aggregates(self.data::Ptr{Cvoid})::Cvoid

@doc raw"""
    cb_call_primal_extender!(self::CBSOCData, param0::CBPrimalExtender, include_candidates::Bool = true)

see the last argument of FunctionOracle::evaluate()
"""
cb_call_primal_extender!(self::CBSOCData, param0::CBPrimalExtender, include_candidates::Bool = true) = @ccall libcb.cb_socdata_call_primal_extender(self.data::Ptr{Cvoid}, param0.data::Ptr{Cvoid}, include_candidates::Cint)::Cint

@doc raw"""
    cb_apply_modification!(self::CBSOCData, param0::CBGroundsetModification, mex::Union{<:CBMinorantExtender,Nothing})

rearrange/extend the minorants according to the given groundset modifications
"""
cb_apply_modification!(self::CBSOCData, param0::CBGroundsetModification, mex::Union{<:CBMinorantExtender,Nothing}) = @ccall libcb.cb_socdata_apply_modification(self.data::Ptr{Cvoid}, param0.data::Ptr{Cvoid}, (isnothing(mex) ? C_NULL : mex.data)::Ptr{Cvoid})::Cint

@doc raw"""
    cb_get_approximate_primal(self::CBSOCData)

return the PrimalData corresponding to the aggregate
"""
cb_get_approximate_primal(self::CBSOCData) = CBPrimalData(@ccall libcb.cb_socdata_get_approximate_primal(self.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_get_center_primal(self::CBSOCData)

return the PrimalData corresponding to the aggregate
"""
cb_get_center_primal(self::CBSOCData) = CBPrimalData(@ccall libcb.cb_socdata_get_center_primal(self.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_get_candidate_primal(self::CBSOCData)

return the PrimalData delivered by the last call of FunctionOracle::evaluate()
"""
cb_get_candidate_primal(self::CBSOCData) = CBPrimalData(@ccall libcb.cb_socdata_get_candidate_primal(self.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_get_latest_minorants!(self::CBSOCData, latest_minorants::CBMinorantBundle, max_number::Integer)

return the max_number latest minorants if available;
"""
cb_get_latest_minorants!(self::CBSOCData, latest_minorants::CBMinorantBundle, max_number::Integer) = @ccall libcb.cb_socdata_get_latest_minorants(self.data::Ptr{Cvoid}, latest_minorants.data::Ptr{Cvoid}, max_number::Cint)::Cint

