@doc raw"""
    cb_clear!(self::CBAFTData, start_modification_id::Integer = 0)

reset to initial state (also used by the default constructor)
"""
cb_clear!(self::CBAFTData, start_modification_id::Integer = 0) = @ccall libcb.cb_aftdata_clear(self.data::Ptr{Cvoid}, start_modification_id::Cint)::Cvoid

@doc raw"""
    CBAFTData(start_modification_id::Integer = 0)

calls clear()
"""
CBAFTData(start_modification_id::Integer = 0) = CBAFTData(@ccall libcb.cb_aftdata_new(start_modification_id::Cint)::Ptr{Cvoid})

@doc raw"""
    cb_init!(self::CBAFTData, bd::Union{<:CBBundleData,Nothing})

initialize from other BundleData
"""
cb_init!(self::CBAFTData, bd::Union{<:CBBundleData,Nothing}) = @ccall libcb.cb_aftdata_init(self.data::Ptr{Cvoid}, (isnothing(bd) ? C_NULL : bd.data)::Ptr{Cvoid})::Cint

@doc raw"""
    cb_clone(self::CBAFTData)

return a pointer to a clone of this
"""
cb_clone(self::CBAFTData) = CBBundleData(@ccall libcb.cb_aftdata_clone(self.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_do_step!(self::CBAFTData, point_id::Integer)

if the candidate information is available and consitent for point_id, copy it from cand to center and return 0, otherwise return 1
"""
cb_do_step!(self::CBAFTData, point_id::Integer) = @ccall libcb.cb_aftdata_do_step(self.data::Ptr{Cvoid}, point_id::Cint)::Cint

@doc raw"""
    cb_synchronize_ids!(self::CBAFTData, new_center_id::Integer, old_center_id::Integer, new_cand_id::Integer, old_cand_id::Integer, new_prex_id::Integer = 0)

synchronize ids in any case by discarding inconsistent parts but return number of errors
"""
function cb_synchronize_ids!(self::CBAFTData, new_center_id::Integer, old_center_id::Integer, new_cand_id::Integer, old_cand_id::Integer, new_prex_id::Integer = 0)
    new_aggregate_id = Ref{Int}()
    new_cand_ub_fid = Ref{Int}()
    new_center_ub_fid = Ref{Int}()
    @ccall libcb.cb_aftdata_synchronize_ids(self.data::Ptr{Cvoid}, new_center_ub_fid::Ref{Int}, new_center_id::Cint, old_center_id::Cint, new_cand_ub_fid::Ref{Int}, new_cand_id::Cint, old_cand_id::Cint, new_aggregate_id::Ref{Int}, new_prex_id::Cint)::Cint
    return new_center_ub_fid[], new_cand_ub_fid[], new_aggregate_id[]
end

@doc raw"""
    cb_center_modified!(self::CBAFTData, center_id::Integer)

check whether center computation is still valid for this point and modification id
"""
function cb_center_modified!(self::CBAFTData, center_id::Integer)
    function_id = Ref{Int}()
    Bool(@ccall libcb.cb_aftdata_center_modified(self.data::Ptr{Cvoid}, function_id::Ref{Int}, center_id::Cint)::Cint)
    return function_id[]
end

@doc raw"""
    cb_model_aggregate_modified!(self::CBAFTData, last_aggr_id::Integer)

check whether aggregate is available and has the same id
"""
cb_model_aggregate_modified!(self::CBAFTData, last_aggr_id::Integer) = Bool(@ccall libcb.cb_aftdata_model_aggregate_modified(self.data::Ptr{Cvoid}, last_aggr_id::Cint)::Cint)

@doc raw"""
    cb_clear_model!(self::CBAFTData, discard_minorants_only::Bool = false)

clear the cutting modle and all function evaluatoins
"""
cb_clear_model!(self::CBAFTData, discard_minorants_only::Bool = false) = @ccall libcb.cb_aftdata_clear_model(self.data::Ptr{Cvoid}, discard_minorants_only::Cint)::Cvoid

@doc raw"""
    cb_clear_aggregates!(self::CBAFTData)

delete all kinds of aggregates but keep explicit parts of the cutting model
"""
cb_clear_aggregates!(self::CBAFTData) = @ccall libcb.cb_aftdata_clear_aggregates(self.data::Ptr{Cvoid})::Cvoid

