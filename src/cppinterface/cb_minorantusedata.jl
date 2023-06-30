@doc raw"""
    CBMinorantUseData(mp::Union{<:CBMinorant,Nothing}, sval::Real, modif_id::Integer)

constructor for a new minorant with scaling @a sval and modification id @a modfi_id
"""
CBMinorantUseData(mp::Union{<:CBMinorant,Nothing}, sval::Real, modif_id::Integer) = CBMinorantUseData(@ccall libcb.cb_minorantusedata_new((isnothing(mp) ? C_NULL : mp.data)::Ptr{Cvoid}, sval::Cdouble, modif_id::Cint)::Ptr{Cvoid})

@doc raw"""
    CBMinorantUseData(mdp::Union{<:CBMinorantUseData,Nothing}, sval::Real)

constructor for a recursively containing further MinorantUseData with an additional scaling factor @a sval
"""
CBMinorantUseData(mdp::Union{<:CBMinorantUseData,Nothing}, sval::Real) = CBMinorantUseData(@ccall libcb.cb_minorantusedata_new2((isnothing(mdp) ? C_NULL : mdp.data)::Ptr{Cvoid}, sval::Cdouble)::Ptr{Cvoid})

@doc raw"""
    cb_valid(self::CBMinorantUseData)

check validity recursively
"""
cb_valid(self::CBMinorantUseData) = Bool(@ccall libcb.cb_minorantusedata_valid(self.data::Ptr{Cvoid})::Cint)

@doc raw"""
    cb_get_minorant(self::CBMinorantUseData)

return the final minorant (by a recursive call) or 0 if there is none
"""
cb_get_minorant(self::CBMinorantUseData) = CBMinorant(@ccall libcb.cb_minorantusedata_get_minorant(self.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_set_modification_id!(self::CBMinorantUseData)

returns the modification id also for overwriting
"""
cb_set_modification_id!(self::CBMinorantUseData) = @ccall libcb.cb_minorantusedata_set_modification_id(self.data::Ptr{Cvoid})::Ptr{Cint}

@doc raw"""
    cb_get_modification_id(self::CBMinorantUseData)

returns the modification id
"""
cb_get_modification_id(self::CBMinorantUseData) = @ccall libcb.cb_minorantusedata_get_modification_id(self.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_synchronize_ids!(self::CBMinorantUseData, new_modification_id::Integer, new_center_id::Integer, old_center_id::Integer, new_cand_id::Integer, old_cand_id::Integer, new_prex_id::Integer)

sets the modification_id to its new id and reinitializes the evaluation map
"""
cb_synchronize_ids!(self::CBMinorantUseData, new_modification_id::Integer, new_center_id::Integer, old_center_id::Integer, new_cand_id::Integer, old_cand_id::Integer, new_prex_id::Integer) = @ccall libcb.cb_minorantusedata_synchronize_ids(self.data::Ptr{Cvoid}, new_modification_id::Cint, new_center_id::Cint, old_center_id::Cint, new_cand_id::Cint, old_cand_id::Cint, new_prex_id::Cint)::Cint

@doc raw"""
    cb_aggregate(self::CBMinorantUseData)

returns true if the minorant is an aggregate of several minorants
"""
cb_aggregate(self::CBMinorantUseData) = Bool(@ccall libcb.cb_minorantusedata_aggregate(self.data::Ptr{Cvoid})::Cint)

@doc raw"""
    cb_aggregated!(self::CBMinorantUseData, n::Integer)

add @a n to the number couting the aggregations
"""
cb_aggregated!(self::CBMinorantUseData, n::Integer) = @ccall libcb.cb_minorantusedata_aggregated(self.data::Ptr{Cvoid}, n::Cint)::Cint

@doc raw"""
    cb_offset(self::CBMinorantUseData)

return the offset of the minorant
"""
cb_offset(self::CBMinorantUseData) = @ccall libcb.cb_minorantusedata_offset(self.data::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_coeff(self::CBMinorantUseData, i::Integer)

return coefficient @a i of the minorant
"""
cb_coeff(self::CBMinorantUseData, i::Integer) = @ccall libcb.cb_minorantusedata_coeff(self.data::Ptr{Cvoid}, i::Cint)::Cdouble

@doc raw"""
    cb_scale!(self::CBMinorantUseData, factor::Real)

carries through the scaling for the underlying minorant, afterwards scaleval==1., may only be carried out if use_cnt==1
"""
cb_scale!(self::CBMinorantUseData, factor::Real) = @ccall libcb.cb_minorantusedata_scale(self.data::Ptr{Cvoid}, factor::Cdouble)::Cint

@doc raw"""
    cb_call_primal_extender!(self::CBMinorantUseData, prex::CBPrimalExtender, in_prex_id::Integer)

apply the primal extender to the minorant if @a in_prex_id indicates a new extension
"""
cb_call_primal_extender!(self::CBMinorantUseData, prex::CBPrimalExtender, in_prex_id::Integer) = @ccall libcb.cb_minorantusedata_call_primal_extender(self.data::Ptr{Cvoid}, prex.data::Ptr{Cvoid}, in_prex_id::Cint)::Cint

@doc raw"""
    cb_evaluate(self::CBMinorantUseData, yid::Integer, y::CBMatrix, with_constant::Bool = true)

evaluate the minorant for @a y unluess @a yid allows to retrieve a previous evaluation
"""
cb_evaluate(self::CBMinorantUseData, yid::Integer, y::CBMatrix, with_constant::Bool = true) = @ccall libcb.cb_minorantusedata_evaluate(self.data::Ptr{Cvoid}, yid::Cint, y.data::Ptr{Cvoid}, with_constant::Cint)::Cdouble

@doc raw"""
    cb_one_user(self::CBMinorantUseData)

returns true if not valid or use_cnt==1 recursively
"""
cb_one_user(self::CBMinorantUseData) = Bool(@ccall libcb.cb_minorantusedata_one_user(self.data::Ptr{Cvoid})::Cint)

