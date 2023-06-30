CBSumBundle() = CBSumBundle(@ccall libcb.cb_sumbundle_new()::Ptr{Cvoid})

@doc raw"""
    CBSumBundle(sb::CBSumBundle)

copy constructor
"""
CBSumBundle(sb::CBSumBundle) = CBSumBundle(@ccall libcb.cb_sumbundle_new2(sb.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_init!(self::CBSumBundle, sb::CBSumBundle)

initialize
"""
cb_init!(self::CBSumBundle, sb::CBSumBundle) = @ccall libcb.cb_sumbundle_init(self.data::Ptr{Cvoid}, sb.data::Ptr{Cvoid})::Cvoid

@doc raw"""
    cb_synchronize_ids!(self::CBSumBundle, new_modification_id::Integer, new_center_id::Integer, old_center_id::Integer, new_cand_id::Integer, old_cand_id::Integer, new_prex_id::Integer = 0)

sets the modification_id to id
"""
cb_synchronize_ids!(self::CBSumBundle, new_modification_id::Integer, new_center_id::Integer, old_center_id::Integer, new_cand_id::Integer, old_cand_id::Integer, new_prex_id::Integer = 0) = @ccall libcb.cb_sumbundle_synchronize_ids(self.data::Ptr{Cvoid}, new_modification_id::Cint, new_center_id::Cint, old_center_id::Cint, new_cand_id::Cint, old_cand_id::Cint, new_prex_id::Cint)::Cvoid

@doc raw"""
    cb_has_bundle_for(self::CBSumBundle, ft::CBFunctionTask)

returns true if BData exists for this mode
"""
cb_has_bundle_for(self::CBSumBundle, ft::CBFunctionTask) = Bool(@ccall libcb.cb_sumbundle_has_bundle_for(self.data::Ptr{Cvoid}, ft::Cint)::Cint)

@doc raw"""
    cb_has_bundle_data(self::CBSumBundle)

returns true if BData exists
"""
cb_has_bundle_data(self::CBSumBundle) = Bool(@ccall libcb.cb_sumbundle_has_bundle_data(self.data::Ptr{Cvoid})::Cint)

@doc raw"""
    cb_bundle_size(self::CBSumBundle, ft::CBFunctionTask)

if BData exists for this mode it returns its bundle_size (possibly 0), otherwise 0
"""
cb_bundle_size(self::CBSumBundle, ft::CBFunctionTask) = @ccall libcb.cb_sumbundle_bundle_size(self.data::Ptr{Cvoid}, ft::Cint)::Cint

@doc raw"""
    cb_has_roots(self::CBSumBundle)

returns true if one of its parts is a root
"""
cb_has_roots(self::CBSumBundle) = Bool(@ccall libcb.cb_sumbundle_has_roots(self.data::Ptr{Cvoid})::Cint)

@doc raw"""
    cb_has_working_roots(self::CBSumBundle)

returns true if one of its parts is a root with n_contributors>0
"""
cb_has_working_roots(self::CBSumBundle) = Bool(@ccall libcb.cb_sumbundle_has_working_roots(self.data::Ptr{Cvoid})::Cint)

@doc raw"""
    cb_active(self::CBSumBundle)

returns true if one of its parts is not inactive
"""
cb_active(self::CBSumBundle) = Bool(@ccall libcb.cb_sumbundle_active(self.data::Ptr{Cvoid})::Cint)

@doc raw"""
    cb_has_contributions(self::CBSumBundle)

returns true if one of its parts is a child
"""
cb_has_contributions(self::CBSumBundle) = Bool(@ccall libcb.cb_sumbundle_has_contributions(self.data::Ptr{Cvoid})::Cint)

@doc raw"""
    cb_get_mode(self::CBSumBundle, ft::CBFunctionTask)

gets the corresponding mode (call only if a has_bundle_for(ft)==true)
"""
cb_get_mode(self::CBSumBundle, ft::CBFunctionTask) = @ccall libcb.cb_sumbundle_get_mode(self.data::Ptr{Cvoid}, ft::Cint)::CBMode

@doc raw"""
    cb_get_function_factor(self::CBSumBundle, ft::CBFunctionTask)

gets the corresponding function factor (call only if a has_bundle_for(ft)==true)
"""
cb_get_function_factor(self::CBSumBundle, ft::CBFunctionTask) = @ccall libcb.cb_sumbundle_get_function_factor(self.data::Ptr{Cvoid}, ft::Cint)::Cdouble

@doc raw"""
    cb_get_n_contributors(self::CBSumBundle, ft::CBFunctionTask)

gets the corresponding n_contributors (call only if a has_bundle_for(ft)==true)
"""
cb_get_n_contributors(self::CBSumBundle, ft::CBFunctionTask) = @ccall libcb.cb_sumbundle_get_n_contributors(self.data::Ptr{Cvoid}, ft::Cint)::Cint

@doc raw"""
    cb_get_bundle(self::CBSumBundle, ft::CBFunctionTask)

gets the corresponding minorants (call only if a has_bundle_for(ft)==true)
"""
cb_get_bundle(self::CBSumBundle, ft::CBFunctionTask) = (@ccall libcb.cb_sumbundle_get_bundle(self.data::Ptr{Cvoid}, ft::Cint)::Ptr{Cvoid}; return self)

@doc raw"""
    cb_get_coeff(self::CBSumBundle, ft::CBFunctionTask)

gets the corresponding aggregation coefficients (call only if a has_bundle_for(ft)==true)
"""
cb_get_coeff(self::CBSumBundle, ft::CBFunctionTask) = (@ccall libcb.cb_sumbundle_get_coeff(self.data::Ptr{Cvoid}, ft::Cint)::Ptr{Cvoid}; return self)

@doc raw"""
    cb_get_aggregate(self::CBSumBundle, ft::CBFunctionTask)

gets the corresponding aggregate (call only if a has_bundle_for(ft)==true)
"""
cb_get_aggregate(self::CBSumBundle, ft::CBFunctionTask) = (@ccall libcb.cb_sumbundle_get_aggregate(self.data::Ptr{Cvoid}, ft::Cint)::Ptr{Cvoid}; return self)

@doc raw"""
    cb_get_cand_minorant(self::CBSumBundle, ft::CBFunctionTask)

gets the corresponding candidate minorant (call only if a has_bundle_for(ft)==true)
"""
cb_get_cand_minorant(self::CBSumBundle, ft::CBFunctionTask) = (@ccall libcb.cb_sumbundle_get_cand_minorant(self.data::Ptr{Cvoid}, ft::Cint)::Ptr{Cvoid}; return self)

@doc raw"""
    cb_get_local_model_aggregate(self::CBSumBundle, aggregate::CBMinorantPointer, factor::Real = 1., aft::Union{<:CBAffineFunctionTransformation,Nothing} = nothing)

get the aggregate that is due to root sumbundle parts handled here
"""
cb_get_local_model_aggregate(self::CBSumBundle, aggregate::CBMinorantPointer, factor::Real = 1., aft::Union{<:CBAffineFunctionTransformation,Nothing} = nothing) = @ccall libcb.cb_sumbundle_get_local_model_aggregate(self.data::Ptr{Cvoid}, aggregate.data::Ptr{Cvoid}, factor::Cdouble, (isnothing(aft) ? C_NULL : aft.data)::Ptr{Cvoid})::Cint

@doc raw"""
    cb_get_contributed_model_aggregate(self::CBSumBundle, aggregate::CBMinorantPointer, factor::Real = 1., aft::Union{<:CBAffineFunctionTransformation,Nothing} = nothing)

get the aggregate that is due to child parts of the sumbundle, which are contributed to parents
"""
cb_get_contributed_model_aggregate(self::CBSumBundle, aggregate::CBMinorantPointer, factor::Real = 1., aft::Union{<:CBAffineFunctionTransformation,Nothing} = nothing) = @ccall libcb.cb_sumbundle_get_contributed_model_aggregate(self.data::Ptr{Cvoid}, aggregate.data::Ptr{Cvoid}, factor::Cdouble, (isnothing(aft) ? C_NULL : aft.data)::Ptr{Cvoid})::Cint

@doc raw"""
    cb_call_primal_extender!(self::CBSumBundle, prex::CBPrimalExtender, prex_id::Integer, ft::CBFunctionTask)

call this primal extender for the primals; if there minorants but no primals or if it fails, return 1
"""
cb_call_primal_extender!(self::CBSumBundle, prex::CBPrimalExtender, prex_id::Integer, ft::CBFunctionTask) = @ccall libcb.cb_sumbundle_call_primal_extender(self.data::Ptr{Cvoid}, prex.data::Ptr{Cvoid}, prex_id::Cint, ft::Cint)::Cint

@doc raw"""
    cb_apply_modification!(self::CBSumBundle, gsmdf::CBGroundsetModification, mod_id::Integer, mex::Union{<:CBMinorantExtender,Nothing}, ft::CBFunctionTask)

rearrange/extend the minorants according to the given groundset modifications
"""
cb_apply_modification!(self::CBSumBundle, gsmdf::CBGroundsetModification, mod_id::Integer, mex::Union{<:CBMinorantExtender,Nothing}, ft::CBFunctionTask) = @ccall libcb.cb_sumbundle_apply_modification(self.data::Ptr{Cvoid}, gsmdf.data::Ptr{Cvoid}, mod_id::Cint, (isnothing(mex) ? C_NULL : mex.data)::Ptr{Cvoid}, ft::Cint)::Cint

