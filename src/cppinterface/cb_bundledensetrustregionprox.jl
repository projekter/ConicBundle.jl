@doc raw"""
    CBBundleDenseTrustRegionProx(Hin::CBSymmatrix, vp::Union{<:CBVariableMetricSelection,Nothing} = nothing, local_metric::Bool = false, bounds_scaling::Bool = false, cbinc::Integer = -1)

initialize to this Matrix and set the variable_metric option (false by default)
"""
CBBundleDenseTrustRegionProx(Hin::CBSymmatrix, vp::Union{<:CBVariableMetricSelection,Nothing} = nothing, local_metric::Bool = false, bounds_scaling::Bool = false, cbinc::Integer = -1) = CBBundleDenseTrustRegionProx(@ccall libcb.cb_bundledensetrustregionprox_new(Hin.data::Ptr{Cvoid}, (isnothing(vp) ? C_NULL : vp.data)::Ptr{Cvoid}, local_metric::Cint, bounds_scaling::Cint, cbinc::Cint)::Ptr{Cvoid})

@doc raw"""
    CBBundleDenseTrustRegionProx(dim::Integer = 0, vp::Union{<:CBVariableMetricSelection,Nothing} = nothing, local_metric::Bool = false, bounds_scaling::Bool = false, cbinc::Integer = -1)

initialize H to the zero Matrix of this dimension (on the diagonal the weight will be added)  and set the variable_metric option (false by default)
"""
CBBundleDenseTrustRegionProx(dim::Integer = 0, vp::Union{<:CBVariableMetricSelection,Nothing} = nothing, local_metric::Bool = false, bounds_scaling::Bool = false, cbinc::Integer = -1) = CBBundleDenseTrustRegionProx(@ccall libcb.cb_bundledensetrustregionprox_new2(dim::Cint, (isnothing(vp) ? C_NULL : vp.data)::Ptr{Cvoid}, local_metric::Cint, bounds_scaling::Cint, cbinc::Cint)::Ptr{Cvoid})

@doc raw"""
    cb_set_weightu!(self::CBBundleDenseTrustRegionProx, in_weightu::Real)

set the weight of the proximal term
"""
cb_set_weightu!(self::CBBundleDenseTrustRegionProx, in_weightu::Real) = @ccall libcb.cb_bundledensetrustregionprox_set_weightu(self.data::Ptr{Cvoid}, in_weightu::Cdouble)::Cvoid

@doc raw"""
    cb_get_weightu(self::CBBundleDenseTrustRegionProx)

returns the current weight of the proximal term
"""
cb_get_weightu(self::CBBundleDenseTrustRegionProx) = @ccall libcb.cb_bundledensetrustregionprox_get_weightu(self.data::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_get_term_corr(self::CBBundleDenseTrustRegionProx)

returns a correction factor for termination precision if the quadratic term is strong
"""
cb_get_term_corr(self::CBBundleDenseTrustRegionProx) = @ccall libcb.cb_bundledensetrustregionprox_get_term_corr(self.data::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_init!(self::CBBundleDenseTrustRegionProx, in_H::CBSymmatrix)

set H with the information, whether it is factored
"""
cb_init!(self::CBBundleDenseTrustRegionProx, in_H::CBSymmatrix) = (@ccall libcb.cb_bundledensetrustregionprox_init(self.data::Ptr{Cvoid}, in_H.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

@doc raw"""
    cb_get_factored(self::CBBundleDenseTrustRegionProx)

returns true iff get_Hchol() returns the factord matrix of H with weightu
"""
cb_get_factored(self::CBBundleDenseTrustRegionProx) = Bool(@ccall libcb.cb_bundledensetrustregionprox_get_factored(self.data::Ptr{Cvoid})::Cint)

@doc raw"""
    cb_get_H(self::CBBundleDenseTrustRegionProx)

returns the metric matrix without weightu
"""
cb_get_H(self::CBBundleDenseTrustRegionProx) = (@ccall libcb.cb_bundledensetrustregionprox_get_h(self.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

@doc raw"""
    cb_get_Hchol(self::CBBundleDenseTrustRegionProx)

returns the stored factorization of H with weightu (up to date if get_factored()==true)
"""
cb_get_Hchol(self::CBBundleDenseTrustRegionProx) = (@ccall libcb.cb_bundledensetrustregionprox_get_hchol(self.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

@doc raw"""
    cb_dim(self::CBBundleDenseTrustRegionProx)

returns the order of the matrix
"""
cb_dim(self::CBBundleDenseTrustRegionProx) = @ccall libcb.cb_bundledensetrustregionprox_dim(self.data::Ptr{Cvoid})::Cint

@doc raw"""
    Base.getindex(self::CBBundleDenseTrustRegionProx, i::Integer, j::Integer)

returns H(i,j) (without including weightu)
"""
Base.getindex(self::CBBundleDenseTrustRegionProx, i::Integer, j::Integer) = @ccall libcb.cb_bundledensetrustregionprox_get(self.data::Ptr{Cvoid}, i::Cint, j::Cint)::Cdouble

@doc raw"""
    cb_norm_sqr(self::CBBundleDenseTrustRegionProx, B::CBMatrix)

returns $\|B\|^2_H$ (with weightu included in H)
"""
cb_norm_sqr(self::CBBundleDenseTrustRegionProx, B::CBMatrix) = @ccall libcb.cb_bundledensetrustregionprox_norm_sqr(self.data::Ptr{Cvoid}, B.data::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_dnorm_sqr(self::CBBundleDenseTrustRegionProx, B::CBMinorantPointer)

returns $\|B\|^2_{H^{-1}}$ (with weightu included in H)
"""
cb_dnorm_sqr(self::CBBundleDenseTrustRegionProx, B::CBMinorantPointer) = @ccall libcb.cb_bundledensetrustregionprox_dnorm_sqr(self.data::Ptr{Cvoid}, B.data::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_is_DLR(self::CBBundleDenseTrustRegionProx)

return true if H is of the form diagonal matrix plus Gram matrix of a low rank matrix
"""
cb_is_DLR(self::CBBundleDenseTrustRegionProx) = Bool(@ccall libcb.cb_bundledensetrustregionprox_is_dlr(self.data::Ptr{Cvoid})::Cint)

@doc raw"""
    cb_add_H(self::CBBundleDenseTrustRegionProx, big_sym::CBSymmatrix, start_index::Integer = 0)

add H to the dense symmetric matrix as a principal submatrix starting at position start_index
"""
cb_add_H(self::CBBundleDenseTrustRegionProx, big_sym::CBSymmatrix, start_index::Integer = 0) = @ccall libcb.cb_bundledensetrustregionprox_add_h(self.data::Ptr{Cvoid}, big_sym.data::Ptr{Cvoid}, start_index::Cint)::Cint

@doc raw"""
    cb_add_Hx(self::CBBundleDenseTrustRegionProx, x::CBMatrix, outplusHx::CBMatrix, alpha::Real = 1.)

adds \f$alpha*Hx\f$ to outplusHx and returns this
"""
cb_add_Hx(self::CBBundleDenseTrustRegionProx, x::CBMatrix, outplusHx::CBMatrix, alpha::Real = 1.) = (@ccall libcb.cb_bundledensetrustregionprox_add_hx(self.data::Ptr{Cvoid}, x.data::Ptr{Cvoid}, outplusHx.data::Ptr{Cvoid}, alpha::Cdouble)::Ptr{Cvoid}; return self)

@doc raw"""
    cb_apply_Hinv(self::CBBundleDenseTrustRegionProx, x::CBMatrix)

returns \f$H^{-1}x\f$
"""
cb_apply_Hinv(self::CBBundleDenseTrustRegionProx, x::CBMatrix) = (@ccall libcb.cb_bundledensetrustregionprox_apply_hinv(self.data::Ptr{Cvoid}, x.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

@doc raw"""
    cb_compute_QP_costs!(self::CBBundleDenseTrustRegionProx, Q::CBSymmatrix, d::CBMatrix, constant_minorant::CBMinorantPointer, bundle::CBMinorantBundle, y::CBMatrix, groundset_minorant::CBMinorantPointer, yfixed::Union{<:CBIndexmatrix,Nothing})

computes the dual QP costs Q, d, and the constant offset to the bundle subproblem, see BundleProxObject::compute_QP_costs
"""
function cb_compute_QP_costs!(self::CBBundleDenseTrustRegionProx, Q::CBSymmatrix, d::CBMatrix, constant_minorant::CBMinorantPointer, bundle::CBMinorantBundle, y::CBMatrix, groundset_minorant::CBMinorantPointer, yfixed::Union{<:CBIndexmatrix,Nothing})
    offset = Ref{Float64}()
    @ccall libcb.cb_bundledensetrustregionprox_compute_qp_costs(self.data::Ptr{Cvoid}, Q.data::Ptr{Cvoid}, d.data::Ptr{Cvoid}, offset::Ref{Float64}, constant_minorant.data::Ptr{Cvoid}, bundle.data::Ptr{Cvoid}, y.data::Ptr{Cvoid}, groundset_minorant.data::Ptr{Cvoid}, (isnothing(yfixed) ? C_NULL : yfixed.data)::Ptr{Cvoid})::Cint
    return offset[]
end

@doc raw"""
    cb_update_QP_costs!(self::CBBundleDenseTrustRegionProx, delta_Q::CBSymmatrix, delta_d::CBMatrix, constant_minorant::CBMinorantPointer, bundle::CBMinorantBundle, center_y::CBMatrix, groundset_minorant::CBMinorantPointer, delta_groundset_minorant::CBMinorantPointer, delta_index::CBIndexmatrix, yfixed::Union{<:CBIndexmatrix,Nothing})

updates the dual QP costs Q, d, and the constant offset to the bundle subproblem, see BundleProxObject::update_QP_costs
"""
function cb_update_QP_costs!(self::CBBundleDenseTrustRegionProx, delta_Q::CBSymmatrix, delta_d::CBMatrix, constant_minorant::CBMinorantPointer, bundle::CBMinorantBundle, center_y::CBMatrix, groundset_minorant::CBMinorantPointer, delta_groundset_minorant::CBMinorantPointer, delta_index::CBIndexmatrix, yfixed::Union{<:CBIndexmatrix,Nothing})
    delta_offset = Ref{Float64}()
    @ccall libcb.cb_bundledensetrustregionprox_update_qp_costs(self.data::Ptr{Cvoid}, delta_Q.data::Ptr{Cvoid}, delta_d.data::Ptr{Cvoid}, delta_offset::Ref{Float64}, constant_minorant.data::Ptr{Cvoid}, bundle.data::Ptr{Cvoid}, center_y.data::Ptr{Cvoid}, groundset_minorant.data::Ptr{Cvoid}, delta_groundset_minorant.data::Ptr{Cvoid}, delta_index.data::Ptr{Cvoid}, (isnothing(yfixed) ? C_NULL : yfixed.data)::Ptr{Cvoid})::Cint
    return delta_offset[]
end

@doc raw"""
    cb_apply_modification!(self::CBBundleDenseTrustRegionProx, gsmdf::CBGroundsetModification)

when BundleSolver is called to modify the groundset it also calls this
"""
cb_apply_modification!(self::CBBundleDenseTrustRegionProx, gsmdf::CBGroundsetModification) = @ccall libcb.cb_bundledensetrustregionprox_apply_modification(self.data::Ptr{Cvoid}, gsmdf.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_projected_clone!(self::CBBundleDenseTrustRegionProx, indices::CBIndexmatrix)

* @brief in order to allow for fixed variables, this generates a clone restricted to the given indices
   
"""
cb_projected_clone!(self::CBBundleDenseTrustRegionProx, indices::CBIndexmatrix) = CBBundleProxObject(@ccall libcb.cb_bundledensetrustregionprox_projected_clone(self.data::Ptr{Cvoid}, indices.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_supports_diagonal_bounds_scaling(self::CBBundleDenseTrustRegionProx)

* @brief this implementation does not support a diagonal scaling heuristic,
        therefore the following routine has to return true.
     
"""
cb_supports_diagonal_bounds_scaling(self::CBBundleDenseTrustRegionProx) = Bool(@ccall libcb.cb_bundledensetrustregionprox_supports_diagonal_bounds_scaling(self.data::Ptr{Cvoid})::Cint)

@doc raw"""
    cb_diagonal_bounds_scaling_update!(self::CBBundleDenseTrustRegionProx, param0::CBMatrix)

* @brief if supported, D_update has to contain nonnegative numbers that
        are permanently added to the diagonal here. It is important to keep
        track of this change only if afterwards update_QP_costs is called before
        compute_QP_costs. In this case the only nonzero enries in D_update must
        be those of delta_index
     
"""
cb_diagonal_bounds_scaling_update!(self::CBBundleDenseTrustRegionProx, param0::CBMatrix) = @ccall libcb.cb_bundledensetrustregionprox_diagonal_bounds_scaling_update(self.data::Ptr{Cvoid}, param0.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_supports_dense_variable_metric(self::CBBundleDenseTrustRegionProx)

returns true if dynamic scaling with dense symmetric matrices is supported
"""
cb_supports_dense_variable_metric(self::CBBundleDenseTrustRegionProx) = Bool(@ccall libcb.cb_bundledensetrustregionprox_supports_dense_variable_metric(self.data::Ptr{Cvoid})::Cint)

@doc raw"""
    cb_supports_lowrank_variable_metric(self::CBBundleDenseTrustRegionProx)

returns true if dynamic scaling with low rank structure is supported
"""
cb_supports_lowrank_variable_metric(self::CBBundleDenseTrustRegionProx) = Bool(@ccall libcb.cb_bundledensetrustregionprox_supports_lowrank_variable_metric(self.data::Ptr{Cvoid})::Cint)

@doc raw"""
    cb_supports_diagonal_variable_metric(self::CBBundleDenseTrustRegionProx)

returns true if dynamic scaling with diagonal matrices is supported
"""
cb_supports_diagonal_variable_metric(self::CBBundleDenseTrustRegionProx) = Bool(@ccall libcb.cb_bundledensetrustregionprox_supports_diagonal_variable_metric(self.data::Ptr{Cvoid})::Cint)

@doc raw"""
    cb_apply_variable_metric!(self::CBBundleDenseTrustRegionProx, groundset::Union{<:CBVariableMetricModel,Nothing}, model::Union{<:CBVariableMetricModel,Nothing}, aggr::CBMatrix, y_id::Integer, y::CBMatrix, descent_step::Bool, model_maxviol::Real, new_indices::Union{<:CBIndexmatrix,Nothing} = nothing)

see DynamicScaling
"""
function cb_apply_variable_metric!(self::CBBundleDenseTrustRegionProx, groundset::Union{<:CBVariableMetricModel,Nothing}, model::Union{<:CBVariableMetricModel,Nothing}, aggr::CBMatrix, y_id::Integer, y::CBMatrix, descent_step::Bool, model_maxviol::Real, new_indices::Union{<:CBIndexmatrix,Nothing} = nothing)
    current_weight = Ref{Float64}()
    @ccall libcb.cb_bundledensetrustregionprox_apply_variable_metric(self.data::Ptr{Cvoid}, (isnothing(groundset) ? C_NULL : groundset.data)::Ptr{Cvoid}, (isnothing(model) ? C_NULL : model.data)::Ptr{Cvoid}, aggr.data::Ptr{Cvoid}, y_id::Cint, y.data::Ptr{Cvoid}, descent_step::Cint, current_weight::Ref{Float64}, model_maxviol::Cdouble, (isnothing(new_indices) ? C_NULL : new_indices.data)::Ptr{Cvoid})::Cint
    return current_weight[]
end

@doc raw"""
    cb_add_variable_metric!(self::CBBundleDenseTrustRegionProx, addH::CBSymmatrix)

see BundleProxObject::add_variable_metric()
"""
cb_add_variable_metric!(self::CBBundleDenseTrustRegionProx, addH::CBSymmatrix) = @ccall libcb.cb_bundledensetrustregionprox_add_variable_metric(self.data::Ptr{Cvoid}, addH.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_add_variable_metric!(self::CBBundleDenseTrustRegionProx, diagH::CBMatrix, vecH::CBMatrix)

see BundleProxObject::add_lowrank_variable_metric()
"""
cb_add_variable_metric!(self::CBBundleDenseTrustRegionProx, diagH::CBMatrix, vecH::CBMatrix) = @ccall libcb.cb_bundledensetrustregionprox_add_variable_metric2(self.data::Ptr{Cvoid}, diagH.data::Ptr{Cvoid}, vecH.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_push_aft!(self::CBBundleDenseTrustRegionProx, aft::Union{<:CBAffineFunctionTransformation,Nothing})

see  BundleProxObject::push_aft();
"""
cb_push_aft!(self::CBBundleDenseTrustRegionProx, aft::Union{<:CBAffineFunctionTransformation,Nothing}) = @ccall libcb.cb_bundledensetrustregionprox_push_aft(self.data::Ptr{Cvoid}, (isnothing(aft) ? C_NULL : aft.data)::Ptr{Cvoid})::Cint

@doc raw"""
    cb_pop_aft!(self::CBBundleDenseTrustRegionProx)

see  BundleProxObject::pop_aft();
"""
cb_pop_aft!(self::CBBundleDenseTrustRegionProx) = @ccall libcb.cb_bundledensetrustregionprox_pop_aft(self.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_mfile_data(self::CBBundleDenseTrustRegionProx)

output the description of the scaling in mfile-suitable format
"""
cb_mfile_data(self::CBBundleDenseTrustRegionProx) = @ccall libcb.cb_bundledensetrustregionprox_mfile_data(self.data::Ptr{Cvoid})::Cint

