@doc raw"""
    cb_clear!(self::CBPSCData, start_modification_id::Integer = 0)

reset to initial state (also used by the default constructor)
"""
cb_clear!(self::CBPSCData, start_modification_id::Integer = 0) = @ccall libcb.cb_pscdata_clear(self.data::Ptr{Cvoid}, start_modification_id::Cint)::Cvoid

@doc raw"""
    CBPSCData(fun_factor::Real = 1., fun_task::CBFunctionTask = cbft_objective_function)

calls clear()
"""
CBPSCData(fun_factor::Real = 1., fun_task::CBFunctionTask = cbft_objective_function) = CBPSCData(@ccall libcb.cb_pscdata_new(fun_factor::Cdouble, fun_task::Cint)::Ptr{Cvoid})

@doc raw"""
    cb_get_primalvecs(self::CBPSCData)

returns primalvecs
"""
cb_get_primalvecs(self::CBPSCData) = (@ccall libcb.cb_pscdata_get_primalvecs(self.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

@doc raw"""
    cb_get_primaleigs(self::CBPSCData)

returns primaleigs
"""
cb_get_primaleigs(self::CBPSCData) = (@ccall libcb.cb_pscdata_get_primaleigs(self.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

@doc raw"""
    cb_get_topvecs(self::CBPSCData)

returns topvecs
"""
cb_get_topvecs(self::CBPSCData) = (@ccall libcb.cb_pscdata_get_topvecs(self.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

@doc raw"""
    cb_get_Ritz_values(self::CBPSCData)

returns Ritz_values
"""
cb_get_Ritz_values(self::CBPSCData) = (@ccall libcb.cb_pscdata_get_ritz_values(self.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

@doc raw"""
    cb_get_keepsize(self::CBPSCData)

returns keepsize
"""
cb_get_keepsize(self::CBPSCData) = @ccall libcb.cb_pscdata_get_keepsize(self.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_get_activedim(self::CBPSCData)

returns acitvedim
"""
cb_get_activedim(self::CBPSCData) = @ccall libcb.cb_pscdata_get_activedim(self.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_get_skippedsize(self::CBPSCData)

returns skippedsize
"""
cb_get_skippedsize(self::CBPSCData) = @ccall libcb.cb_pscdata_get_skippedsize(self.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_init!(self::CBPSCData, bd::Union{<:CBBundleData,Nothing})

if @a bd is of type PSCData, initialize to this data
"""
cb_init!(self::CBPSCData, bd::Union{<:CBBundleData,Nothing}) = @ccall libcb.cb_pscdata_init(self.data::Ptr{Cvoid}, (isnothing(bd) ? C_NULL : bd.data)::Ptr{Cvoid})::Cint

@doc raw"""
    cb_clone(self::CBPSCData)

return a pointer to a clone of this
"""
cb_clone(self::CBPSCData) = CBBundleData(@ccall libcb.cb_pscdata_clone(self.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_do_step!(self::CBPSCData, point_id::Integer)

if the candidate information is available and consitent for point_id, copy it from cand to center and return 0, otherwise return 1
"""
cb_do_step!(self::CBPSCData, point_id::Integer) = @ccall libcb.cb_pscdata_do_step(self.data::Ptr{Cvoid}, point_id::Cint)::Cint

@doc raw"""
    cb_synchronize_ids!(self::CBPSCData, new_center_id::Integer, old_center_id::Integer, new_cand_id::Integer, old_cand_id::Integer, new_prex_id::Integer = 0)

synchronize ids in any case by discarding inconsistent parts but return number of errors
"""
function cb_synchronize_ids!(self::CBPSCData, new_center_id::Integer, old_center_id::Integer, new_cand_id::Integer, old_cand_id::Integer, new_prex_id::Integer = 0)
    new_aggregate_id = Ref{Int}()
    new_cand_ub_fid = Ref{Int}()
    new_center_ub_fid = Ref{Int}()
    @ccall libcb.cb_pscdata_synchronize_ids(self.data::Ptr{Cvoid}, new_center_ub_fid::Ref{Int}, new_center_id::Cint, old_center_id::Cint, new_cand_ub_fid::Ref{Int}, new_cand_id::Cint, old_cand_id::Cint, new_aggregate_id::Ref{Int}, new_prex_id::Cint)::Cint
    return new_center_ub_fid[], new_cand_ub_fid[], new_aggregate_id[]
end

@doc raw"""
    cb_clear_model!(self::CBPSCData, discard_minorants_only::Bool = false)

clear the cutting model and all function evaluations
"""
cb_clear_model!(self::CBPSCData, discard_minorants_only::Bool = false) = @ccall libcb.cb_pscdata_clear_model(self.data::Ptr{Cvoid}, discard_minorants_only::Cint)::Cvoid

@doc raw"""
    cb_clear_aggregates!(self::CBPSCData)

delete all kinds of aggregates but keep explicit parts of the cutting model
"""
cb_clear_aggregates!(self::CBPSCData) = @ccall libcb.cb_pscdata_clear_aggregates(self.data::Ptr{Cvoid})::Cvoid

@doc raw"""
    cb_call_primal_extender!(self::CBPSCData, param0::CBPrimalExtender, include_candidates::Bool = true)

see the last argument of FunctionOracle::evaluate()
"""
cb_call_primal_extender!(self::CBPSCData, param0::CBPrimalExtender, include_candidates::Bool = true) = @ccall libcb.cb_pscdata_call_primal_extender(self.data::Ptr{Cvoid}, param0.data::Ptr{Cvoid}, include_candidates::Cint)::Cint

@doc raw"""
    cb_apply_modification!(self::CBPSCData, param0::CBGroundsetModification, mex::Union{<:CBMinorantExtender,Nothing})

rearrange/extend the minorants according to the given groundset modifications
"""
cb_apply_modification!(self::CBPSCData, param0::CBGroundsetModification, mex::Union{<:CBMinorantExtender,Nothing}) = @ccall libcb.cb_pscdata_apply_modification(self.data::Ptr{Cvoid}, param0.data::Ptr{Cvoid}, (isnothing(mex) ? C_NULL : mex.data)::Ptr{Cvoid})::Cint

@doc raw"""
    cb_get_approximate_primal(self::CBPSCData)

return the PrimalData corresponding to the aggregate
"""
cb_get_approximate_primal(self::CBPSCData) = CBPrimalData(@ccall libcb.cb_pscdata_get_approximate_primal(self.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_get_center_primal(self::CBPSCData)

return the PrimalData corresponding to the aggregate
"""
cb_get_center_primal(self::CBPSCData) = CBPrimalData(@ccall libcb.cb_pscdata_get_center_primal(self.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_get_candidate_primal(self::CBPSCData)

return the PrimalData delivered by the last call of FunctionOracle::evaluate()
"""
cb_get_candidate_primal(self::CBPSCData) = CBPrimalData(@ccall libcb.cb_pscdata_get_candidate_primal(self.data::Ptr{Cvoid})::Ptr{Cvoid})

