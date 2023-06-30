@doc raw"""
    CBBundleHKWeight(mRin::Real = .5, bwp::Union{<:CBBundleWeight,Nothing} = nothing, incr::Integer = -1)

the parameter mRin gets the value for accepting descent steps, bwp may be used to communicate the previous values use by another routine
"""
CBBundleHKWeight(mRin::Real = .5, bwp::Union{<:CBBundleWeight,Nothing} = nothing, incr::Integer = -1) = CBBundleHKWeight(@ccall libcb.cb_bundlehkweight_new(mRin::Cdouble, (isnothing(bwp) ? C_NULL : bwp.data)::Ptr{Cvoid}, incr::Cint)::Ptr{Cvoid})

@doc raw"""
    cb_set_defaults!(self::CBBundleHKWeight)

set default values for 'constant' parameters, e.g. minweight and maxweight
"""
cb_set_defaults!(self::CBBundleHKWeight) = @ccall libcb.cb_bundlehkweight_set_defaults(self.data::Ptr{Cvoid})::Cvoid

@doc raw"""
    cb_set_nullstep_updates!(self::CBBundleHKWeight, nu::Integer = 0)

set nullstep update strategy (0 ... original, 1 ... none, 2 ... enlarge if subsequence of three norm increases is found
"""
cb_set_nullstep_updates!(self::CBBundleHKWeight, nu::Integer = 0) = @ccall libcb.cb_bundlehkweight_set_nullstep_updates(self.data::Ptr{Cvoid}, nu::Cint)::Cvoid

@doc raw"""
    cb_clear!(self::CBBundleHKWeight)

reset all adaptive variables and parameters
"""
cb_clear!(self::CBBundleHKWeight) = @ccall libcb.cb_bundlehkweight_clear(self.data::Ptr{Cvoid})::Cvoid

@doc raw"""
    cb_init!(self::CBBundleHKWeight, aggr_dnmormsqr::Real, groundset::Union{<:CBGroundset,Nothing}, model::Union{<:CBBundleModel,Nothing})

compute first weight and set some parameters
"""
cb_init!(self::CBBundleHKWeight, aggr_dnmormsqr::Real, groundset::Union{<:CBGroundset,Nothing}, model::Union{<:CBBundleModel,Nothing}) = @ccall libcb.cb_bundlehkweight_init(self.data::Ptr{Cvoid}, aggr_dnmormsqr::Cdouble, (isnothing(groundset) ? C_NULL : groundset.data)::Ptr{Cvoid}, (isnothing(model) ? C_NULL : model.data)::Ptr{Cvoid})::Cint

@doc raw"""
    cb_set_next_weight!(self::CBBundleHKWeight, u::Real)

<=0 leaves everything unchanged and does nothing
"""
cb_set_next_weight!(self::CBBundleHKWeight, u::Real) = @ccall libcb.cb_bundlehkweight_set_next_weight(self.data::Ptr{Cvoid}, u::Cdouble)::Cvoid

@doc raw"""
    cb_set_minweight!(self::CBBundleHKWeight, mw::Real)

<=0 means no bound
"""
cb_set_minweight!(self::CBBundleHKWeight, mw::Real) = @ccall libcb.cb_bundlehkweight_set_minweight(self.data::Ptr{Cvoid}, mw::Cdouble)::Cvoid

@doc raw"""
    cb_get_next_weight_set(self::CBBundleHKWeight)

true if the next weight was prespecified externally
"""
cb_get_next_weight_set(self::CBBundleHKWeight) = Bool(@ccall libcb.cb_bundlehkweight_get_next_weight_set(self.data::Ptr{Cvoid})::Cint)

cb_get_minweight(self::CBBundleHKWeight) = @ccall libcb.cb_bundlehkweight_get_minweight(self.data::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_set_maxweight!(self::CBBundleHKWeight, mw::Real)

<=0 means no bound
"""
cb_set_maxweight!(self::CBBundleHKWeight, mw::Real) = @ccall libcb.cb_bundlehkweight_set_maxweight(self.data::Ptr{Cvoid}, mw::Cdouble)::Cvoid

cb_get_maxweight(self::CBBundleHKWeight) = @ccall libcb.cb_bundlehkweight_get_maxweight(self.data::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_get_weight(self::CBBundleHKWeight)

returns current value of the weight
"""
cb_get_weight(self::CBBundleHKWeight) = @ccall libcb.cb_bundlehkweight_get_weight(self.data::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_weight_changed(self::CBBundleHKWeight)

returns true if last call of *_update modified current value of tau, else 0
"""
cb_weight_changed(self::CBBundleHKWeight) = Bool(@ccall libcb.cb_bundlehkweight_weight_changed(self.data::Ptr{Cvoid})::Cint)

@doc raw"""
    cb_descent_update!(self::CBBundleHKWeight, newval::Real, oldval::Real, modelval::Real, y::CBMatrix, newy::CBMatrix, normsubg2::Real, Hp::Union{<:CBBundleProxObject,Nothing})

determine next weight after a descent step
"""
cb_descent_update!(self::CBBundleHKWeight, newval::Real, oldval::Real, modelval::Real, y::CBMatrix, newy::CBMatrix, normsubg2::Real, Hp::Union{<:CBBundleProxObject,Nothing}) = @ccall libcb.cb_bundlehkweight_descent_update(self.data::Ptr{Cvoid}, newval::Cdouble, oldval::Cdouble, modelval::Cdouble, y.data::Ptr{Cvoid}, newy.data::Ptr{Cvoid}, normsubg2::Cdouble, (isnothing(Hp) ? C_NULL : Hp.data)::Ptr{Cvoid})::Cint

@doc raw"""
    cb_nullstep_update!(self::CBBundleHKWeight, newval::Real, oldval::Real, modelval::Real, y::CBMatrix, newy::CBMatrix, new_minorant::CBMinorantPointer, aggregate::CBMinorantPointer, nullstep_bound::Real, normsubg2::Real, Hp::Union{<:CBBundleProxObject,Nothing})

determine next weight after a null step
"""
cb_nullstep_update!(self::CBBundleHKWeight, newval::Real, oldval::Real, modelval::Real, y::CBMatrix, newy::CBMatrix, new_minorant::CBMinorantPointer, aggregate::CBMinorantPointer, nullstep_bound::Real, normsubg2::Real, Hp::Union{<:CBBundleProxObject,Nothing}) = @ccall libcb.cb_bundlehkweight_nullstep_update(self.data::Ptr{Cvoid}, newval::Cdouble, oldval::Cdouble, modelval::Cdouble, y.data::Ptr{Cvoid}, newy.data::Ptr{Cvoid}, new_minorant.data::Ptr{Cvoid}, aggregate.data::Ptr{Cvoid}, nullstep_bound::Cdouble, normsubg2::Cdouble, (isnothing(Hp) ? C_NULL : Hp.data)::Ptr{Cvoid})::Cint

@doc raw"""
    cb_apply_modification!(self::CBBundleHKWeight, gsmdf::CBGroundsetModification)

reinitialize after modifications
"""
cb_apply_modification!(self::CBBundleHKWeight, gsmdf::CBGroundsetModification) = @ccall libcb.cb_bundlehkweight_apply_modification(self.data::Ptr{Cvoid}, gsmdf.data::Ptr{Cvoid})::Cint

