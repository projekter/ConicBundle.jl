@doc raw"""
    CBBundleDLRTrustRegionProx(dim::Integer = 0, vp::Union{<:CBVariableMetricSelection,Nothing} = nothing, local_metric::Bool = false, inc::Integer = -1)

default constructor with empty H (equal to zero) and the dimension as argument
"""
CBBundleDLRTrustRegionProx(dim::Integer = 0, vp::Union{<:CBVariableMetricSelection,Nothing} = nothing, local_metric::Bool = false, inc::Integer = -1) = CBBundleDLRTrustRegionProx(@ccall libcb.cb_bundledlrtrustregionprox_new(dim::Cint, (isnothing(vp) ? C_NULL : vp.data)::Ptr{Cvoid}, local_metric::Cint, inc::Cint)::Ptr{Cvoid})

@doc raw"""
    CBBundleDLRTrustRegionProx(in_D::CBMatrix, in_vecH::CBMatrix, vp::Union{<:CBVariableMetricSelection,Nothing} = nothing, local_metric::Bool = false, inc::Integer = -1)

constructs H=(D+weightu)+vecH*transpose(vecH) with weightu=1., so in_D must be a column vector with same dimension as rows in in_vecH
"""
CBBundleDLRTrustRegionProx(in_D::CBMatrix, in_vecH::CBMatrix, vp::Union{<:CBVariableMetricSelection,Nothing} = nothing, local_metric::Bool = false, inc::Integer = -1) = CBBundleDLRTrustRegionProx(@ccall libcb.cb_bundledlrtrustregionprox_new2(in_D.data::Ptr{Cvoid}, in_vecH.data::Ptr{Cvoid}, (isnothing(vp) ? C_NULL : vp.data)::Ptr{Cvoid}, local_metric::Cint, inc::Cint)::Ptr{Cvoid})

@doc raw"""
    cb_set_weightu!(self::CBBundleDLRTrustRegionProx, in_weightu::Real)

sets the next weight
"""
cb_set_weightu!(self::CBBundleDLRTrustRegionProx, in_weightu::Real) = @ccall libcb.cb_bundledlrtrustregionprox_set_weightu(self.data::Ptr{Cvoid}, in_weightu::Cdouble)::Cvoid

@doc raw"""
    cb_get_weightu(self::CBBundleDLRTrustRegionProx)

returns the current weight in use
"""
cb_get_weightu(self::CBBundleDLRTrustRegionProx) = @ccall libcb.cb_bundledlrtrustregionprox_get_weightu(self.data::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_get_term_corr(self::CBBundleDLRTrustRegionProx)

returns a correction factor for termination precision if the quadratic term is strong
"""
cb_get_term_corr(self::CBBundleDLRTrustRegionProx) = @ccall libcb.cb_bundledlrtrustregionprox_get_term_corr(self.data::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_init!(self::CBBundleDLRTrustRegionProx, in_D::CBMatrix, in_vecH::CBMatrix)

reset the prox information; the diagaonal part in_D must be a nonnegative column vector, the low rank part in_vecH must have the same number of rows as in_D
"""
cb_init!(self::CBBundleDLRTrustRegionProx, in_D::CBMatrix, in_vecH::CBMatrix) = @ccall libcb.cb_bundledlrtrustregionprox_init(self.data::Ptr{Cvoid}, in_D.data::Ptr{Cvoid}, in_vecH.data::Ptr{Cvoid})::Cvoid

@doc raw"""
    cb_norm_sqr(self::CBBundleDLRTrustRegionProx, B::CBMatrix)

returns $\|B\|^2_H$ (with weight included)
"""
cb_norm_sqr(self::CBBundleDLRTrustRegionProx, B::CBMatrix) = @ccall libcb.cb_bundledlrtrustregionprox_norm_sqr(self.data::Ptr{Cvoid}, B.data::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_dnorm_sqr(self::CBBundleDLRTrustRegionProx, B::CBMinorantPointer)

returns $\|B\|^2_{H^{-1}}$ (with weight included)
"""
cb_dnorm_sqr(self::CBBundleDLRTrustRegionProx, B::CBMinorantPointer) = @ccall libcb.cb_bundledlrtrustregionprox_dnorm_sqr(self.data::Ptr{Cvoid}, B.data::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_is_DLR(self::CBBundleDLRTrustRegionProx)

return true if H is of the form diagonal matrix plus Gram matrix of a low rank matrix
"""
cb_is_DLR(self::CBBundleDLRTrustRegionProx) = Bool(@ccall libcb.cb_bundledlrtrustregionprox_is_dlr(self.data::Ptr{Cvoid})::Cint)

@doc raw"""
    cb_add_H(self::CBBundleDLRTrustRegionProx, big_sym::CBSymmatrix, start_index::Integer = 0)

add H to the dense symmetric matrix as a principal submatrix starting at position start_index
"""
cb_add_H(self::CBBundleDLRTrustRegionProx, big_sym::CBSymmatrix, start_index::Integer = 0) = @ccall libcb.cb_bundledlrtrustregionprox_add_h(self.data::Ptr{Cvoid}, big_sym.data::Ptr{Cvoid}, start_index::Cint)::Cint

@doc raw"""
    cb_add_Hx(self::CBBundleDLRTrustRegionProx, x::CBMatrix, outplusHx::CBMatrix, alpha::Real = 1.)

adds \f$alpha*Hx\f$ to outplusHx and returns this
"""
cb_add_Hx(self::CBBundleDLRTrustRegionProx, x::CBMatrix, outplusHx::CBMatrix, alpha::Real = 1.) = (@ccall libcb.cb_bundledlrtrustregionprox_add_hx(self.data::Ptr{Cvoid}, x.data::Ptr{Cvoid}, outplusHx.data::Ptr{Cvoid}, alpha::Cdouble)::Ptr{Cvoid}; return self)

@doc raw"""
    cb_apply_Hinv(self::CBBundleDLRTrustRegionProx, x::CBMatrix)

returns \f$H^{-1}x\f$ where \f$H^{-1}=\frac1u(I-V(\Lambda/(\Lambda+u))V^\top)\f$
"""
cb_apply_Hinv(self::CBBundleDLRTrustRegionProx, x::CBMatrix) = (@ccall libcb.cb_bundledlrtrustregionprox_apply_hinv(self.data::Ptr{Cvoid}, x.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

@doc raw"""
    cb_compute_QP_costs!(self::CBBundleDLRTrustRegionProx, Q::CBSymmatrix, d::CBMatrix, constant_minorant::CBMinorantPointer, bundle::CBMinorantBundle, y::CBMatrix, groundset_minorant::CBMinorantPointer, yfixed::Union{<:CBIndexmatrix,Nothing})

computes the dual QP costs Q, d, and the constant offset to the bundle subproblem, see BundleProxObject::compute_QP_costs
"""
function cb_compute_QP_costs!(self::CBBundleDLRTrustRegionProx, Q::CBSymmatrix, d::CBMatrix, constant_minorant::CBMinorantPointer, bundle::CBMinorantBundle, y::CBMatrix, groundset_minorant::CBMinorantPointer, yfixed::Union{<:CBIndexmatrix,Nothing})
    offset = Ref{Float64}()
    @ccall libcb.cb_bundledlrtrustregionprox_compute_qp_costs(self.data::Ptr{Cvoid}, Q.data::Ptr{Cvoid}, d.data::Ptr{Cvoid}, offset::Ref{Float64}, constant_minorant.data::Ptr{Cvoid}, bundle.data::Ptr{Cvoid}, y.data::Ptr{Cvoid}, groundset_minorant.data::Ptr{Cvoid}, (isnothing(yfixed) ? C_NULL : yfixed.data)::Ptr{Cvoid})::Cint
    return offset[]
end

@doc raw"""
    cb_update_QP_costs!(self::CBBundleDLRTrustRegionProx, delta_Q::CBSymmatrix, delta_d::CBMatrix, constant_minorant::CBMinorantPointer, bundle::CBMinorantBundle, center_y::CBMatrix, groundset_minorant::CBMinorantPointer, delta_groundset_minorant::CBMinorantPointer, delta_index::CBIndexmatrix, yfixed::Union{<:CBIndexmatrix,Nothing})

updates the dual QP costs Q, d, and the constant offset to the bundle subproblem, see BundleProxObject::update_QP_costs
"""
function cb_update_QP_costs!(self::CBBundleDLRTrustRegionProx, delta_Q::CBSymmatrix, delta_d::CBMatrix, constant_minorant::CBMinorantPointer, bundle::CBMinorantBundle, center_y::CBMatrix, groundset_minorant::CBMinorantPointer, delta_groundset_minorant::CBMinorantPointer, delta_index::CBIndexmatrix, yfixed::Union{<:CBIndexmatrix,Nothing})
    delta_offset = Ref{Float64}()
    @ccall libcb.cb_bundledlrtrustregionprox_update_qp_costs(self.data::Ptr{Cvoid}, delta_Q.data::Ptr{Cvoid}, delta_d.data::Ptr{Cvoid}, delta_offset::Ref{Float64}, constant_minorant.data::Ptr{Cvoid}, bundle.data::Ptr{Cvoid}, center_y.data::Ptr{Cvoid}, groundset_minorant.data::Ptr{Cvoid}, delta_groundset_minorant.data::Ptr{Cvoid}, delta_index.data::Ptr{Cvoid}, (isnothing(yfixed) ? C_NULL : yfixed.data)::Ptr{Cvoid})::Cint
    return delta_offset[]
end

@doc raw"""
    cb_apply_modification!(self::CBBundleDLRTrustRegionProx, gsmdf::CBGroundsetModification)

when BundleSolver is called to modify the groundset it also calls this
"""
cb_apply_modification!(self::CBBundleDLRTrustRegionProx, gsmdf::CBGroundsetModification) = @ccall libcb.cb_bundledlrtrustregionprox_apply_modification(self.data::Ptr{Cvoid}, gsmdf.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_projected_clone!(self::CBBundleDLRTrustRegionProx, indices::CBIndexmatrix)

* @brief in order to allow for fixed variables, this generates a clone restricted to the given indices
   
"""
cb_projected_clone!(self::CBBundleDLRTrustRegionProx, indices::CBIndexmatrix) = CBBundleProxObject(@ccall libcb.cb_bundledlrtrustregionprox_projected_clone(self.data::Ptr{Cvoid}, indices.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_supports_diagonal_bounds_scaling(self::CBBundleDLRTrustRegionProx)

* @brief this implementation does not support a diagonal scaling heuristic,
        therefore the following routine has to return true.
     
"""
cb_supports_diagonal_bounds_scaling(self::CBBundleDLRTrustRegionProx) = Bool(@ccall libcb.cb_bundledlrtrustregionprox_supports_diagonal_bounds_scaling(self.data::Ptr{Cvoid})::Cint)

@doc raw"""
    cb_diagonal_bounds_scaling_update!(self::CBBundleDLRTrustRegionProx, param0::CBMatrix)

* @brief if supported, D_update has to contain nonnegative numbers that
        are permanently added to the diagonal here. It is important to keep
        track of this change only if afterwards update_QP_costs is called before
        compute_QP_costs. In this case the only nonzero enries in D_update must
        be those of delta_index
     
"""
cb_diagonal_bounds_scaling_update!(self::CBBundleDLRTrustRegionProx, param0::CBMatrix) = @ccall libcb.cb_bundledlrtrustregionprox_diagonal_bounds_scaling_update(self.data::Ptr{Cvoid}, param0.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_supports_dense_variable_metric(self::CBBundleDLRTrustRegionProx)

returns true if variable metric with dense symmetric matrices is supported
"""
cb_supports_dense_variable_metric(self::CBBundleDLRTrustRegionProx) = Bool(@ccall libcb.cb_bundledlrtrustregionprox_supports_dense_variable_metric(self.data::Ptr{Cvoid})::Cint)

@doc raw"""
    cb_supports_lowrank_variable_metric(self::CBBundleDLRTrustRegionProx)

returns true if variable metric with low rank structure is supported
"""
cb_supports_lowrank_variable_metric(self::CBBundleDLRTrustRegionProx) = Bool(@ccall libcb.cb_bundledlrtrustregionprox_supports_lowrank_variable_metric(self.data::Ptr{Cvoid})::Cint)

@doc raw"""
    cb_supports_diagonal_variable_metric(self::CBBundleDLRTrustRegionProx)

returns true if variable metric with diagonal matrices is supported
"""
cb_supports_diagonal_variable_metric(self::CBBundleDLRTrustRegionProx) = Bool(@ccall libcb.cb_bundledlrtrustregionprox_supports_diagonal_variable_metric(self.data::Ptr{Cvoid})::Cint)

@doc raw"""
    cb_apply_variable_metric!(self::CBBundleDLRTrustRegionProx, groundset::Union{<:CBVariableMetricModel,Nothing}, model::Union{<:CBVariableMetricModel,Nothing}, aggr::CBMatrix, y_id::Integer, y::CBMatrix, descent_step::Bool, model_maxviol::Real, new_indices::Union{<:CBIndexmatrix,Nothing} = nothing)

see DynamicScaling
"""
function cb_apply_variable_metric!(self::CBBundleDLRTrustRegionProx, groundset::Union{<:CBVariableMetricModel,Nothing}, model::Union{<:CBVariableMetricModel,Nothing}, aggr::CBMatrix, y_id::Integer, y::CBMatrix, descent_step::Bool, model_maxviol::Real, new_indices::Union{<:CBIndexmatrix,Nothing} = nothing)
    current_weight = Ref{Float64}()
    @ccall libcb.cb_bundledlrtrustregionprox_apply_variable_metric(self.data::Ptr{Cvoid}, (isnothing(groundset) ? C_NULL : groundset.data)::Ptr{Cvoid}, (isnothing(model) ? C_NULL : model.data)::Ptr{Cvoid}, aggr.data::Ptr{Cvoid}, y_id::Cint, y.data::Ptr{Cvoid}, descent_step::Cint, current_weight::Ref{Float64}, model_maxviol::Cdouble, (isnothing(new_indices) ? C_NULL : new_indices.data)::Ptr{Cvoid})::Cint
    return current_weight[]
end

@doc raw"""
    cb_add_variable_metric!(self::CBBundleDLRTrustRegionProx, diagH::CBMatrix, vecH::CBMatrix)

see BundleProxObject::add_variable_metric()
"""
cb_add_variable_metric!(self::CBBundleDLRTrustRegionProx, diagH::CBMatrix, vecH::CBMatrix) = @ccall libcb.cb_bundledlrtrustregionprox_add_variable_metric(self.data::Ptr{Cvoid}, diagH.data::Ptr{Cvoid}, vecH.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_push_aft!(self::CBBundleDLRTrustRegionProx, aft::Union{<:CBAffineFunctionTransformation,Nothing})

see  BundleProxObject::push_aft();
"""
cb_push_aft!(self::CBBundleDLRTrustRegionProx, aft::Union{<:CBAffineFunctionTransformation,Nothing}) = @ccall libcb.cb_bundledlrtrustregionprox_push_aft(self.data::Ptr{Cvoid}, (isnothing(aft) ? C_NULL : aft.data)::Ptr{Cvoid})::Cint

@doc raw"""
    cb_pop_aft!(self::CBBundleDLRTrustRegionProx)

see  BundleProxObject::pop_aft();
"""
cb_pop_aft!(self::CBBundleDLRTrustRegionProx) = @ccall libcb.cb_bundledlrtrustregionprox_pop_aft(self.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_mfile_data(self::CBBundleDLRTrustRegionProx)

output the description of the scaling in mfile-suitable format
"""
cb_mfile_data(self::CBBundleDLRTrustRegionProx) = @ccall libcb.cb_bundledlrtrustregionprox_mfile_data(self.data::Ptr{Cvoid})::Cint

