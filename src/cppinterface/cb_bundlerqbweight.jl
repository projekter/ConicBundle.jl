@doc raw"""
    CBBundleRQBWeight(bwp::Union{<:CBBundleWeight,Nothing} = nothing, incr::Integer = -1)

bwp may be used to communicate the previous values used by another routine
"""
CBBundleRQBWeight(bwp::Union{<:CBBundleWeight,Nothing} = nothing, incr::Integer = -1) = CBBundleRQBWeight(@ccall libcb.cb_bundlerqbweight_new((isnothing(bwp) ? C_NULL : bwp.data)::Ptr{Cvoid}, incr::Cint)::Ptr{Cvoid})

@doc raw"""
    CBBundleRQBWeight(m1::Real, m2::Real = .2, m3::Real = 1., eta::Real = 1e-6, bwp::Union{<:CBBundleWeight,Nothing} = nothing, incr::Integer = -1)

bwp may be used to communicate the previous values used by another routine
"""
CBBundleRQBWeight(m1::Real, m2::Real = .2, m3::Real = 1., eta::Real = 1e-6, bwp::Union{<:CBBundleWeight,Nothing} = nothing, incr::Integer = -1) = CBBundleRQBWeight(@ccall libcb.cb_bundlerqbweight_new2(m1::Cdouble, m2::Cdouble, m3::Cdouble, eta::Cdouble, (isnothing(bwp) ? C_NULL : bwp.data)::Ptr{Cvoid}, incr::Cint)::Ptr{Cvoid})

@doc raw"""
    cb_set_defaults!(self::CBBundleRQBWeight)

set default values for 'constant' parameters, e.g. minweight and maxweight
"""
cb_set_defaults!(self::CBBundleRQBWeight) = @ccall libcb.cb_bundlerqbweight_set_defaults(self.data::Ptr{Cvoid})::Cvoid

@doc raw"""
    cb_clear!(self::CBBundleRQBWeight)

reset all adaptive variables and parameters
"""
cb_clear!(self::CBBundleRQBWeight) = @ccall libcb.cb_bundlerqbweight_clear(self.data::Ptr{Cvoid})::Cvoid

@doc raw"""
    cb_init!(self::CBBundleRQBWeight, aggr_dnmormsqr::Real, groundset::Union{<:CBGroundset,Nothing}, model::Union{<:CBBundleModel,Nothing})

compute first weight and set some parameters
"""
cb_init!(self::CBBundleRQBWeight, aggr_dnmormsqr::Real, groundset::Union{<:CBGroundset,Nothing}, model::Union{<:CBBundleModel,Nothing}) = @ccall libcb.cb_bundlerqbweight_init(self.data::Ptr{Cvoid}, aggr_dnmormsqr::Cdouble, (isnothing(groundset) ? C_NULL : groundset.data)::Ptr{Cvoid}, (isnothing(model) ? C_NULL : model.data)::Ptr{Cvoid})::Cint

@doc raw"""
    cb_set_next_weight!(self::CBBundleRQBWeight, u::Real)

<=0 leaves everything unchanged and does nothing
"""
cb_set_next_weight!(self::CBBundleRQBWeight, u::Real) = @ccall libcb.cb_bundlerqbweight_set_next_weight(self.data::Ptr{Cvoid}, u::Cdouble)::Cvoid

@doc raw"""
    cb_set_minweight!(self::CBBundleRQBWeight, mw::Real)

<=0 means no bound
"""
cb_set_minweight!(self::CBBundleRQBWeight, mw::Real) = @ccall libcb.cb_bundlerqbweight_set_minweight(self.data::Ptr{Cvoid}, mw::Cdouble)::Cvoid

@doc raw"""
    cb_get_next_weight_set(self::CBBundleRQBWeight)

true if the next weight was prespecified externally
"""
cb_get_next_weight_set(self::CBBundleRQBWeight) = Bool(@ccall libcb.cb_bundlerqbweight_get_next_weight_set(self.data::Ptr{Cvoid})::Cint)

cb_get_minweight(self::CBBundleRQBWeight) = @ccall libcb.cb_bundlerqbweight_get_minweight(self.data::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_set_maxweight!(self::CBBundleRQBWeight, mw::Real)

<=0 means no bound
"""
cb_set_maxweight!(self::CBBundleRQBWeight, mw::Real) = @ccall libcb.cb_bundlerqbweight_set_maxweight(self.data::Ptr{Cvoid}, mw::Cdouble)::Cvoid

cb_get_maxweight(self::CBBundleRQBWeight) = @ccall libcb.cb_bundlerqbweight_get_maxweight(self.data::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_get_weight(self::CBBundleRQBWeight)

returns current value of the weight
"""
cb_get_weight(self::CBBundleRQBWeight) = @ccall libcb.cb_bundlerqbweight_get_weight(self.data::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_weight_changed(self::CBBundleRQBWeight)

returns true if last call of *_update modified current value of tau, else 0
"""
cb_weight_changed(self::CBBundleRQBWeight) = Bool(@ccall libcb.cb_bundlerqbweight_weight_changed(self.data::Ptr{Cvoid})::Cint)

@doc raw"""
    cb_descent_update!(self::CBBundleRQBWeight, newval::Real, oldval::Real, modelval::Real, y::CBMatrix, newy::CBMatrix, normsubg2::Real, Hp::Union{<:CBBundleProxObject,Nothing})

determine next weight after a descent step
"""
cb_descent_update!(self::CBBundleRQBWeight, newval::Real, oldval::Real, modelval::Real, y::CBMatrix, newy::CBMatrix, normsubg2::Real, Hp::Union{<:CBBundleProxObject,Nothing}) = @ccall libcb.cb_bundlerqbweight_descent_update(self.data::Ptr{Cvoid}, newval::Cdouble, oldval::Cdouble, modelval::Cdouble, y.data::Ptr{Cvoid}, newy.data::Ptr{Cvoid}, normsubg2::Cdouble, (isnothing(Hp) ? C_NULL : Hp.data)::Ptr{Cvoid})::Cint

@doc raw"""
    cb_nullstep_update!(self::CBBundleRQBWeight, newval::Real, oldval::Real, modelval::Real, y::CBMatrix, newy::CBMatrix, new_minorant::CBMinorantPointer, aggregate::CBMinorantPointer, nullstep_bound::Real, normsubg2::Real, Hp::Union{<:CBBundleProxObject,Nothing})

determine next weight after a null step
"""
cb_nullstep_update!(self::CBBundleRQBWeight, newval::Real, oldval::Real, modelval::Real, y::CBMatrix, newy::CBMatrix, new_minorant::CBMinorantPointer, aggregate::CBMinorantPointer, nullstep_bound::Real, normsubg2::Real, Hp::Union{<:CBBundleProxObject,Nothing}) = @ccall libcb.cb_bundlerqbweight_nullstep_update(self.data::Ptr{Cvoid}, newval::Cdouble, oldval::Cdouble, modelval::Cdouble, y.data::Ptr{Cvoid}, newy.data::Ptr{Cvoid}, new_minorant.data::Ptr{Cvoid}, aggregate.data::Ptr{Cvoid}, nullstep_bound::Cdouble, normsubg2::Cdouble, (isnothing(Hp) ? C_NULL : Hp.data)::Ptr{Cvoid})::Cint

@doc raw"""
    cb_apply_modification!(self::CBBundleRQBWeight, gsmdf::CBGroundsetModification)

reinitialize after modifications
"""
cb_apply_modification!(self::CBBundleRQBWeight, gsmdf::CBGroundsetModification) = @ccall libcb.cb_bundlerqbweight_apply_modification(self.data::Ptr{Cvoid}, gsmdf.data::Ptr{Cvoid})::Cint

