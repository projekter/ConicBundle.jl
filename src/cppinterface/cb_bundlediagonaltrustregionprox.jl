@doc raw"""
    CBBundleDiagonalTrustRegionProx(Din::CBMatrix, vp::Union{<:CBVariableMetricSelection,Nothing} = nothing, local_scaling::Bool = false, bounds_scaling::Bool = false, cbinc::Integer = -1)

initialize to this diagonal matrix
"""
CBBundleDiagonalTrustRegionProx(Din::CBMatrix, vp::Union{<:CBVariableMetricSelection,Nothing} = nothing, local_scaling::Bool = false, bounds_scaling::Bool = false, cbinc::Integer = -1) = CBBundleDiagonalTrustRegionProx(@ccall libcb.cb_bundlediagonaltrustregionprox_new(Din.data::Ptr{Cvoid}, (isnothing(vp) ? C_NULL : vp.data)::Ptr{Cvoid}, local_scaling::Cint, bounds_scaling::Cint, cbinc::Cint)::Ptr{Cvoid})

@doc raw"""
    CBBundleDiagonalTrustRegionProx(dim::Integer, d::Real, vp::Union{<:CBVariableMetricSelection,Nothing} = nothing, local_scaling::Bool = false, bounds_scaling::Bool = false, cbinc::Integer = -1)

initialize to a diaognal matrix d*identity of dimesion dim
"""
CBBundleDiagonalTrustRegionProx(dim::Integer, d::Real, vp::Union{<:CBVariableMetricSelection,Nothing} = nothing, local_scaling::Bool = false, bounds_scaling::Bool = false, cbinc::Integer = -1) = CBBundleDiagonalTrustRegionProx(@ccall libcb.cb_bundlediagonaltrustregionprox_new2(dim::Cint, d::Cdouble, (isnothing(vp) ? C_NULL : vp.data)::Ptr{Cvoid}, local_scaling::Cint, bounds_scaling::Cint, cbinc::Cint)::Ptr{Cvoid})

@doc raw"""
    CBBundleDiagonalTrustRegionProx(dim::Integer = 0, vp::Union{<:CBVariableMetricSelection,Nothing} = nothing, local_scaling::Bool = false, bounds_scaling::Bool = false, cbinc::Integer = -1)

initialize to a zero diaognal matrix of dimesion dim
"""
CBBundleDiagonalTrustRegionProx(dim::Integer = 0, vp::Union{<:CBVariableMetricSelection,Nothing} = nothing, local_scaling::Bool = false, bounds_scaling::Bool = false, cbinc::Integer = -1) = CBBundleDiagonalTrustRegionProx(@ccall libcb.cb_bundlediagonaltrustregionprox_new3(dim::Cint, (isnothing(vp) ? C_NULL : vp.data)::Ptr{Cvoid}, local_scaling::Cint, bounds_scaling::Cint, cbinc::Cint)::Ptr{Cvoid})

@doc raw"""
    cb_set_weightu!(self::CBBundleDiagonalTrustRegionProx, in_weightu::Real)

set the weight of the proximal term
"""
cb_set_weightu!(self::CBBundleDiagonalTrustRegionProx, in_weightu::Real) = @ccall libcb.cb_bundlediagonaltrustregionprox_set_weightu(self.data::Ptr{Cvoid}, in_weightu::Cdouble)::Cvoid

@doc raw"""
    cb_get_weightu(self::CBBundleDiagonalTrustRegionProx)

returns the current weight of the proximal term
"""
cb_get_weightu(self::CBBundleDiagonalTrustRegionProx) = @ccall libcb.cb_bundlediagonaltrustregionprox_get_weightu(self.data::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_get_term_corr(self::CBBundleDiagonalTrustRegionProx)

returns a correction factor for termination precision if the quadratic term is strong
"""
cb_get_term_corr(self::CBBundleDiagonalTrustRegionProx) = @ccall libcb.cb_bundlediagonaltrustregionprox_get_term_corr(self.data::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_get_D(self::CBBundleDiagonalTrustRegionProx)

returns the diagonal D of the diagonal scaling matrix
"""
cb_get_D(self::CBBundleDiagonalTrustRegionProx) = (@ccall libcb.cb_bundlediagonaltrustregionprox_get_d(self.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

@doc raw"""
    cb_set_D!(self::CBBundleDiagonalTrustRegionProx, in_D::CBMatrix)

set the diagonal (it needs to be >=0 but this is not checked)
"""
cb_set_D!(self::CBBundleDiagonalTrustRegionProx, in_D::CBMatrix) = @ccall libcb.cb_bundlediagonaltrustregionprox_set_d(self.data::Ptr{Cvoid}, in_D.data::Ptr{Cvoid})::Cvoid

@doc raw"""
    cb_dim(self::CBBundleDiagonalTrustRegionProx)

returns the dimension of the diagonal
"""
cb_dim(self::CBBundleDiagonalTrustRegionProx) = @ccall libcb.cb_bundlediagonaltrustregionprox_dim(self.data::Ptr{Cvoid})::Cint

@doc raw"""
    Base.getindex(self::CBBundleDiagonalTrustRegionProx, i::Integer)

return the i-th element of the diagonal matrix D
"""
Base.getindex(self::CBBundleDiagonalTrustRegionProx, i::Integer) = @ccall libcb.cb_bundlediagonaltrustregionprox_get(self.data::Ptr{Cvoid}, i::Cint)::Cdouble

@doc raw"""
    cb_norm_sqr(self::CBBundleDiagonalTrustRegionProx, B::CBMatrix)

returns \f$\|B\|^2_H\f$
"""
cb_norm_sqr(self::CBBundleDiagonalTrustRegionProx, B::CBMatrix) = @ccall libcb.cb_bundlediagonaltrustregionprox_norm_sqr(self.data::Ptr{Cvoid}, B.data::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_dnorm_sqr(self::CBBundleDiagonalTrustRegionProx, B::CBMinorantPointer)

returns \f$\|B\|^2_{H^{-1}}\f$
"""
cb_dnorm_sqr(self::CBBundleDiagonalTrustRegionProx, B::CBMinorantPointer) = @ccall libcb.cb_bundlediagonaltrustregionprox_dnorm_sqr(self.data::Ptr{Cvoid}, B.data::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_is_DLR(self::CBBundleDiagonalTrustRegionProx)

return true if H is of the form diagonal matrix plus Gram matrix of a low rank matrix
"""
cb_is_DLR(self::CBBundleDiagonalTrustRegionProx) = Bool(@ccall libcb.cb_bundlediagonaltrustregionprox_is_dlr(self.data::Ptr{Cvoid})::Cint)

@doc raw"""
    cb_add_H(self::CBBundleDiagonalTrustRegionProx, big_sym::CBSymmatrix, start_index::Integer = 0)

add H to the dense symmetric matrix as a principal submatrix starting at position start_index
"""
cb_add_H(self::CBBundleDiagonalTrustRegionProx, big_sym::CBSymmatrix, start_index::Integer = 0) = @ccall libcb.cb_bundlediagonaltrustregionprox_add_h(self.data::Ptr{Cvoid}, big_sym.data::Ptr{Cvoid}, start_index::Cint)::Cint

@doc raw"""
    cb_add_Hx(self::CBBundleDiagonalTrustRegionProx, x::CBMatrix, outplusHx::CBMatrix, alpha::Real = 1.)

adds \f$alpha*Hx\f$ to outplusHx and returns this
"""
cb_add_Hx(self::CBBundleDiagonalTrustRegionProx, x::CBMatrix, outplusHx::CBMatrix, alpha::Real = 1.) = (@ccall libcb.cb_bundlediagonaltrustregionprox_add_hx(self.data::Ptr{Cvoid}, x.data::Ptr{Cvoid}, outplusHx.data::Ptr{Cvoid}, alpha::Cdouble)::Ptr{Cvoid}; return self)

@doc raw"""
    cb_apply_Hinv(self::CBBundleDiagonalTrustRegionProx, x::CBMatrix)

returns \f$H^{-1}x\f$
"""
cb_apply_Hinv(self::CBBundleDiagonalTrustRegionProx, x::CBMatrix) = (@ccall libcb.cb_bundlediagonaltrustregionprox_apply_hinv(self.data::Ptr{Cvoid}, x.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

@doc raw"""
    cb_compute_QP_costs!(self::CBBundleDiagonalTrustRegionProx, Q::CBSymmatrix, d::CBMatrix, constant_minorant::CBMinorantPointer, bundle::CBMinorantBundle, y::CBMatrix, groundset_minorant::CBMinorantPointer, yfixed::Union{<:CBIndexmatrix,Nothing})

computes the dual QP costs Q, d, and the constant offset to the bundle subproblem, see BundleProxObject::compute_QP_costs
"""
function cb_compute_QP_costs!(self::CBBundleDiagonalTrustRegionProx, Q::CBSymmatrix, d::CBMatrix, constant_minorant::CBMinorantPointer, bundle::CBMinorantBundle, y::CBMatrix, groundset_minorant::CBMinorantPointer, yfixed::Union{<:CBIndexmatrix,Nothing})
    offset = Ref{Float64}()
    @ccall libcb.cb_bundlediagonaltrustregionprox_compute_qp_costs(self.data::Ptr{Cvoid}, Q.data::Ptr{Cvoid}, d.data::Ptr{Cvoid}, offset::Ref{Float64}, constant_minorant.data::Ptr{Cvoid}, bundle.data::Ptr{Cvoid}, y.data::Ptr{Cvoid}, groundset_minorant.data::Ptr{Cvoid}, (isnothing(yfixed) ? C_NULL : yfixed.data)::Ptr{Cvoid})::Cint
    return offset[]
end

@doc raw"""
    cb_update_QP_costs!(self::CBBundleDiagonalTrustRegionProx, delta_Q::CBSymmatrix, delta_d::CBMatrix, constant_minorant::CBMinorantPointer, bundle::CBMinorantBundle, center_y::CBMatrix, groundset_minorant::CBMinorantPointer, delta_groundset_minorant::CBMinorantPointer, delta_index::CBIndexmatrix, yfixed::Union{<:CBIndexmatrix,Nothing})

updates the dual QP costs Q, d, and the constant offset to the bundle subproblem, see BundleProxObject::update_QP_costs
"""
function cb_update_QP_costs!(self::CBBundleDiagonalTrustRegionProx, delta_Q::CBSymmatrix, delta_d::CBMatrix, constant_minorant::CBMinorantPointer, bundle::CBMinorantBundle, center_y::CBMatrix, groundset_minorant::CBMinorantPointer, delta_groundset_minorant::CBMinorantPointer, delta_index::CBIndexmatrix, yfixed::Union{<:CBIndexmatrix,Nothing})
    delta_offset = Ref{Float64}()
    @ccall libcb.cb_bundlediagonaltrustregionprox_update_qp_costs(self.data::Ptr{Cvoid}, delta_Q.data::Ptr{Cvoid}, delta_d.data::Ptr{Cvoid}, delta_offset::Ref{Float64}, constant_minorant.data::Ptr{Cvoid}, bundle.data::Ptr{Cvoid}, center_y.data::Ptr{Cvoid}, groundset_minorant.data::Ptr{Cvoid}, delta_groundset_minorant.data::Ptr{Cvoid}, delta_index.data::Ptr{Cvoid}, (isnothing(yfixed) ? C_NULL : yfixed.data)::Ptr{Cvoid})::Cint
    return delta_offset[]
end

@doc raw"""
    cb_apply_modification!(self::CBBundleDiagonalTrustRegionProx, gsmdf::CBGroundsetModification)

when BundleSolver is called to modify the groundset it also calls this
"""
cb_apply_modification!(self::CBBundleDiagonalTrustRegionProx, gsmdf::CBGroundsetModification) = @ccall libcb.cb_bundlediagonaltrustregionprox_apply_modification(self.data::Ptr{Cvoid}, gsmdf.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_projected_clone!(self::CBBundleDiagonalTrustRegionProx, indices::CBIndexmatrix)

* @brief in order to allow for fixed variables, this generates a clone restricted to the given indices
   
"""
cb_projected_clone!(self::CBBundleDiagonalTrustRegionProx, indices::CBIndexmatrix) = CBBundleProxObject(@ccall libcb.cb_bundlediagonaltrustregionprox_projected_clone(self.data::Ptr{Cvoid}, indices.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_supports_diagonal_bounds_scaling(self::CBBundleDiagonalTrustRegionProx)

* @brief this implementation supports a diagonal scaling heuristic for
        bounds in the groundset, therefore the following routine has to return true.
     
"""
cb_supports_diagonal_bounds_scaling(self::CBBundleDiagonalTrustRegionProx) = Bool(@ccall libcb.cb_bundlediagonaltrustregionprox_supports_diagonal_bounds_scaling(self.data::Ptr{Cvoid})::Cint)

@doc raw"""
    cb_diagonal_bounds_scaling_update!(self::CBBundleDiagonalTrustRegionProx, D_update::CBMatrix)

* @brief if supported, D_update has to contain nonnegative numbers that
        are permanently added to the diagonal here. It is important to keep
        track of this change only if afterwards update_QP_costs is called before
        compute_QP_costs. In this case the only nonzero enries in D_update must
        be those of delta_index
     
"""
cb_diagonal_bounds_scaling_update!(self::CBBundleDiagonalTrustRegionProx, D_update::CBMatrix) = @ccall libcb.cb_bundlediagonaltrustregionprox_diagonal_bounds_scaling_update(self.data::Ptr{Cvoid}, D_update.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_supports_lowrank_variable_metric(self::CBBundleDiagonalTrustRegionProx)

returns true if dynamic scaling with low rank matrices is supported
"""
cb_supports_lowrank_variable_metric(self::CBBundleDiagonalTrustRegionProx) = Bool(@ccall libcb.cb_bundlediagonaltrustregionprox_supports_lowrank_variable_metric(self.data::Ptr{Cvoid})::Cint)

@doc raw"""
    cb_supports_diagonal_variable_metric(self::CBBundleDiagonalTrustRegionProx)

returns true if dynamic scaling with diagonal matrices is supported
"""
cb_supports_diagonal_variable_metric(self::CBBundleDiagonalTrustRegionProx) = Bool(@ccall libcb.cb_bundlediagonaltrustregionprox_supports_diagonal_variable_metric(self.data::Ptr{Cvoid})::Cint)

@doc raw"""
    cb_apply_variable_metric!(self::CBBundleDiagonalTrustRegionProx, groundset::Union{<:CBVariableMetricModel,Nothing}, model::Union{<:CBVariableMetricModel,Nothing}, aggr::CBMatrix, y_id::Integer, y::CBMatrix, descent_step::Bool, model_maxviol::Real, new_indices::Union{<:CBIndexmatrix,Nothing} = nothing)

see VariableMetric
"""
function cb_apply_variable_metric!(self::CBBundleDiagonalTrustRegionProx, groundset::Union{<:CBVariableMetricModel,Nothing}, model::Union{<:CBVariableMetricModel,Nothing}, aggr::CBMatrix, y_id::Integer, y::CBMatrix, descent_step::Bool, model_maxviol::Real, new_indices::Union{<:CBIndexmatrix,Nothing} = nothing)
    current_weight = Ref{Float64}()
    @ccall libcb.cb_bundlediagonaltrustregionprox_apply_variable_metric(self.data::Ptr{Cvoid}, (isnothing(groundset) ? C_NULL : groundset.data)::Ptr{Cvoid}, (isnothing(model) ? C_NULL : model.data)::Ptr{Cvoid}, aggr.data::Ptr{Cvoid}, y_id::Cint, y.data::Ptr{Cvoid}, descent_step::Cint, current_weight::Ref{Float64}, model_maxviol::Cdouble, (isnothing(new_indices) ? C_NULL : new_indices.data)::Ptr{Cvoid})::Cint
    return current_weight[]
end

@doc raw"""
    cb_add_variable_metric!(self::CBBundleDiagonalTrustRegionProx, diagH::CBMatrix, vecH::CBMatrix)

see BundleProxObject::add_dynamic_scaling()
"""
cb_add_variable_metric!(self::CBBundleDiagonalTrustRegionProx, diagH::CBMatrix, vecH::CBMatrix) = @ccall libcb.cb_bundlediagonaltrustregionprox_add_variable_metric(self.data::Ptr{Cvoid}, diagH.data::Ptr{Cvoid}, vecH.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_push_aft!(self::CBBundleDiagonalTrustRegionProx, aft::Union{<:CBAffineFunctionTransformation,Nothing})

see  BundleProxObject::push_aft();
"""
cb_push_aft!(self::CBBundleDiagonalTrustRegionProx, aft::Union{<:CBAffineFunctionTransformation,Nothing}) = @ccall libcb.cb_bundlediagonaltrustregionprox_push_aft(self.data::Ptr{Cvoid}, (isnothing(aft) ? C_NULL : aft.data)::Ptr{Cvoid})::Cint

@doc raw"""
    cb_pop_aft!(self::CBBundleDiagonalTrustRegionProx)

see  BundleProxObject::pop_aft();
"""
cb_pop_aft!(self::CBBundleDiagonalTrustRegionProx) = @ccall libcb.cb_bundlediagonaltrustregionprox_pop_aft(self.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_mfile_data(self::CBBundleDiagonalTrustRegionProx)

output the description of the scaling in mfile-suitable format
"""
cb_mfile_data(self::CBBundleDiagonalTrustRegionProx) = @ccall libcb.cb_bundlediagonaltrustregionprox_mfile_data(self.data::Ptr{Cvoid})::Cint

