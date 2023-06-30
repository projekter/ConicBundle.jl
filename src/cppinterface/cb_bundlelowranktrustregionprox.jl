@doc raw"""
    CBBundleLowRankTrustRegionProx(dim::Integer = 0, vp::Union{<:CBVariableMetricSelection,Nothing} = nothing, local_metric::Bool = false, inc::Integer = -1)

default constructor with empty H (equal to zero) and the dimension as argument
"""
CBBundleLowRankTrustRegionProx(dim::Integer = 0, vp::Union{<:CBVariableMetricSelection,Nothing} = nothing, local_metric::Bool = false, inc::Integer = -1) = CBBundleLowRankTrustRegionProx(@ccall libcb.cb_bundlelowranktrustregionprox_new(dim::Cint, (isnothing(vp) ? C_NULL : vp.data)::Ptr{Cvoid}, local_metric::Cint, inc::Cint)::Ptr{Cvoid})

@doc raw"""
    CBBundleLowRankTrustRegionProx(in_vecH::CBMatrix, in_lamH::CBMatrix, vp::Union{<:CBVariableMetricSelection,Nothing} = nothing, local_metric::Bool = false, inc::Integer = -1)

constructs H=vecH*Diag(lamH)*transpose(vecH), so in_lamH must be a column vector with same dimension as columns in in_vecH and the row dimension of in_vecH must match the design space, in_vecH is assumed orthogonal
"""
CBBundleLowRankTrustRegionProx(in_vecH::CBMatrix, in_lamH::CBMatrix, vp::Union{<:CBVariableMetricSelection,Nothing} = nothing, local_metric::Bool = false, inc::Integer = -1) = CBBundleLowRankTrustRegionProx(@ccall libcb.cb_bundlelowranktrustregionprox_new2(in_vecH.data::Ptr{Cvoid}, in_lamH.data::Ptr{Cvoid}, (isnothing(vp) ? C_NULL : vp.data)::Ptr{Cvoid}, local_metric::Cint, inc::Cint)::Ptr{Cvoid})

@doc raw"""
    cb_set_weightu!(self::CBBundleLowRankTrustRegionProx, in_weightu::Real)

sets the next weight
"""
cb_set_weightu!(self::CBBundleLowRankTrustRegionProx, in_weightu::Real) = @ccall libcb.cb_bundlelowranktrustregionprox_set_weightu(self.data::Ptr{Cvoid}, in_weightu::Cdouble)::Cvoid

@doc raw"""
    cb_get_weightu(self::CBBundleLowRankTrustRegionProx)

returns the current weight in use
"""
cb_get_weightu(self::CBBundleLowRankTrustRegionProx) = @ccall libcb.cb_bundlelowranktrustregionprox_get_weightu(self.data::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_get_term_corr(self::CBBundleLowRankTrustRegionProx)

returns a correction factor for termination precision if the quadratic term is strong
"""
cb_get_term_corr(self::CBBundleLowRankTrustRegionProx) = @ccall libcb.cb_bundlelowranktrustregionprox_get_term_corr(self.data::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_init!(self::CBBundleLowRankTrustRegionProx, in_vecH::CBMatrix, in_lamH::CBMatrix)

reset the prox information; in_vecH must be an orthogonal matrix with in_vecH.rowdim() matching the dimension but maybe with zero columns; in_lamH must be a column vector with row dimension matching the column dimension of in_vecH and all entries positive
"""
cb_init!(self::CBBundleLowRankTrustRegionProx, in_vecH::CBMatrix, in_lamH::CBMatrix) = @ccall libcb.cb_bundlelowranktrustregionprox_init(self.data::Ptr{Cvoid}, in_vecH.data::Ptr{Cvoid}, in_lamH.data::Ptr{Cvoid})::Cvoid

@doc raw"""
    cb_norm_sqr(self::CBBundleLowRankTrustRegionProx, B::CBMatrix)

returns $\|B\|^2_H$ (with weight included)
"""
cb_norm_sqr(self::CBBundleLowRankTrustRegionProx, B::CBMatrix) = @ccall libcb.cb_bundlelowranktrustregionprox_norm_sqr(self.data::Ptr{Cvoid}, B.data::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_dnorm_sqr(self::CBBundleLowRankTrustRegionProx, B::CBMinorantPointer)

returns $\|B\|^2_{H^{-1}}$ (with weight included)
"""
cb_dnorm_sqr(self::CBBundleLowRankTrustRegionProx, B::CBMinorantPointer) = @ccall libcb.cb_bundlelowranktrustregionprox_dnorm_sqr(self.data::Ptr{Cvoid}, B.data::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_is_DLR(self::CBBundleLowRankTrustRegionProx)

return true if H is of the form diagonal matrix plus Gram matrix of a low rank matrix
"""
cb_is_DLR(self::CBBundleLowRankTrustRegionProx) = Bool(@ccall libcb.cb_bundlelowranktrustregionprox_is_dlr(self.data::Ptr{Cvoid})::Cint)

@doc raw"""
    cb_add_H(self::CBBundleLowRankTrustRegionProx, big_sym::CBSymmatrix, start_index::Integer = 0)

add H to the dense symmetric matrix as a principal submatrix starting at position start_index
"""
cb_add_H(self::CBBundleLowRankTrustRegionProx, big_sym::CBSymmatrix, start_index::Integer = 0) = @ccall libcb.cb_bundlelowranktrustregionprox_add_h(self.data::Ptr{Cvoid}, big_sym.data::Ptr{Cvoid}, start_index::Cint)::Cint

@doc raw"""
    cb_add_Hx(self::CBBundleLowRankTrustRegionProx, x::CBMatrix, outplusHx::CBMatrix, alpha::Real = 1.)

adds \f$alpha*Hx\f$ to outplusHx and returns this
"""
cb_add_Hx(self::CBBundleLowRankTrustRegionProx, x::CBMatrix, outplusHx::CBMatrix, alpha::Real = 1.) = (@ccall libcb.cb_bundlelowranktrustregionprox_add_hx(self.data::Ptr{Cvoid}, x.data::Ptr{Cvoid}, outplusHx.data::Ptr{Cvoid}, alpha::Cdouble)::Ptr{Cvoid}; return self)

@doc raw"""
    cb_apply_Hinv(self::CBBundleLowRankTrustRegionProx, x::CBMatrix)

returns \f$H^{-1}x\f$ where \f$H^{-1}=\frac1u(I-V(\Lambda/(\Lambda+u))V^\top)\f$
"""
cb_apply_Hinv(self::CBBundleLowRankTrustRegionProx, x::CBMatrix) = (@ccall libcb.cb_bundlelowranktrustregionprox_apply_hinv(self.data::Ptr{Cvoid}, x.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

@doc raw"""
    cb_compute_QP_costs!(self::CBBundleLowRankTrustRegionProx, Q::CBSymmatrix, d::CBMatrix, constant_minorant::CBMinorantPointer, bundle::CBMinorantBundle, y::CBMatrix, groundset_minorant::CBMinorantPointer, yfixed::Union{<:CBIndexmatrix,Nothing})

computes the dual QP costs Q, d, and the constant offset to the bundle subproblem, see BundleProxObject::compute_QP_costs
"""
function cb_compute_QP_costs!(self::CBBundleLowRankTrustRegionProx, Q::CBSymmatrix, d::CBMatrix, constant_minorant::CBMinorantPointer, bundle::CBMinorantBundle, y::CBMatrix, groundset_minorant::CBMinorantPointer, yfixed::Union{<:CBIndexmatrix,Nothing})
    offset = Ref{Float64}()
    @ccall libcb.cb_bundlelowranktrustregionprox_compute_qp_costs(self.data::Ptr{Cvoid}, Q.data::Ptr{Cvoid}, d.data::Ptr{Cvoid}, offset::Ref{Float64}, constant_minorant.data::Ptr{Cvoid}, bundle.data::Ptr{Cvoid}, y.data::Ptr{Cvoid}, groundset_minorant.data::Ptr{Cvoid}, (isnothing(yfixed) ? C_NULL : yfixed.data)::Ptr{Cvoid})::Cint
    return offset[]
end

@doc raw"""
    cb_update_QP_costs!(self::CBBundleLowRankTrustRegionProx, delta_Q::CBSymmatrix, delta_d::CBMatrix, constant_minorant::CBMinorantPointer, bundle::CBMinorantBundle, center_y::CBMatrix, groundset_minorant::CBMinorantPointer, delta_groundset_minorant::CBMinorantPointer, delta_index::CBIndexmatrix, yfixed::Union{<:CBIndexmatrix,Nothing})

updates the dual QP costs Q, d, and the constant offset to the bundle subproblem, see BundleProxObject::update_QP_costs
"""
function cb_update_QP_costs!(self::CBBundleLowRankTrustRegionProx, delta_Q::CBSymmatrix, delta_d::CBMatrix, constant_minorant::CBMinorantPointer, bundle::CBMinorantBundle, center_y::CBMatrix, groundset_minorant::CBMinorantPointer, delta_groundset_minorant::CBMinorantPointer, delta_index::CBIndexmatrix, yfixed::Union{<:CBIndexmatrix,Nothing})
    delta_offset = Ref{Float64}()
    @ccall libcb.cb_bundlelowranktrustregionprox_update_qp_costs(self.data::Ptr{Cvoid}, delta_Q.data::Ptr{Cvoid}, delta_d.data::Ptr{Cvoid}, delta_offset::Ref{Float64}, constant_minorant.data::Ptr{Cvoid}, bundle.data::Ptr{Cvoid}, center_y.data::Ptr{Cvoid}, groundset_minorant.data::Ptr{Cvoid}, delta_groundset_minorant.data::Ptr{Cvoid}, delta_index.data::Ptr{Cvoid}, (isnothing(yfixed) ? C_NULL : yfixed.data)::Ptr{Cvoid})::Cint
    return delta_offset[]
end

@doc raw"""
    cb_apply_modification!(self::CBBundleLowRankTrustRegionProx, gsmdf::CBGroundsetModification)

when BundleSolver is called to modify the groundset it also calls this
"""
cb_apply_modification!(self::CBBundleLowRankTrustRegionProx, gsmdf::CBGroundsetModification) = @ccall libcb.cb_bundlelowranktrustregionprox_apply_modification(self.data::Ptr{Cvoid}, gsmdf.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_projected_clone!(self::CBBundleLowRankTrustRegionProx, indices::CBIndexmatrix)

* @brief in order to allow for fixed variables, this generates a clone restricted to the given indices
   
"""
cb_projected_clone!(self::CBBundleLowRankTrustRegionProx, indices::CBIndexmatrix) = CBBundleProxObject(@ccall libcb.cb_bundlelowranktrustregionprox_projected_clone(self.data::Ptr{Cvoid}, indices.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_supports_diagonal_bounds_scaling(self::CBBundleLowRankTrustRegionProx)

* @brief this implementation does not support a diagonal scaling heuristic,
        therefore the following routine has to return true.
     
"""
cb_supports_diagonal_bounds_scaling(self::CBBundleLowRankTrustRegionProx) = Bool(@ccall libcb.cb_bundlelowranktrustregionprox_supports_diagonal_bounds_scaling(self.data::Ptr{Cvoid})::Cint)

@doc raw"""
    cb_diagonal_scaling_heuristic_update!(self::CBBundleLowRankTrustRegionProx, param0::CBMatrix)

* @brief if supported, D_update has to contain nonnegative numbers that
        are permanently added to the diagonal here. It is important to keep
        track of this change only if afterwards update_QP_costs is called before
        compute_QP_costs. In this case the only nonzero enries in D_update must
        be those of delta_index
     
"""
cb_diagonal_scaling_heuristic_update!(self::CBBundleLowRankTrustRegionProx, param0::CBMatrix) = @ccall libcb.cb_bundlelowranktrustregionprox_diagonal_scaling_heuristic_update(self.data::Ptr{Cvoid}, param0.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_supports_dense_variable_metric(self::CBBundleLowRankTrustRegionProx)

returns true if dynamic scaling with dense symmetric matrices is supported
"""
cb_supports_dense_variable_metric(self::CBBundleLowRankTrustRegionProx) = Bool(@ccall libcb.cb_bundlelowranktrustregionprox_supports_dense_variable_metric(self.data::Ptr{Cvoid})::Cint)

@doc raw"""
    cb_supports_lowrank_variable_metric(self::CBBundleLowRankTrustRegionProx)

returns true if dynamic scaling with low rank structure is supported
"""
cb_supports_lowrank_variable_metric(self::CBBundleLowRankTrustRegionProx) = Bool(@ccall libcb.cb_bundlelowranktrustregionprox_supports_lowrank_variable_metric(self.data::Ptr{Cvoid})::Cint)

@doc raw"""
    cb_supports_diagonal_variable_metric(self::CBBundleLowRankTrustRegionProx)

returns true if dynamic scaling with diagonal matrices is supported
"""
cb_supports_diagonal_variable_metric(self::CBBundleLowRankTrustRegionProx) = Bool(@ccall libcb.cb_bundlelowranktrustregionprox_supports_diagonal_variable_metric(self.data::Ptr{Cvoid})::Cint)

@doc raw"""
    cb_apply_variable_metric!(self::CBBundleLowRankTrustRegionProx, groundset::Union{<:CBVariableMetricModel,Nothing}, model::Union{<:CBVariableMetricModel,Nothing}, aggr::CBMatrix, y_id::Integer, y::CBMatrix, descent_step::Bool, model_maxviol::Real, new_indices::Union{<:CBIndexmatrix,Nothing} = nothing)

see DynamicProx
"""
function cb_apply_variable_metric!(self::CBBundleLowRankTrustRegionProx, groundset::Union{<:CBVariableMetricModel,Nothing}, model::Union{<:CBVariableMetricModel,Nothing}, aggr::CBMatrix, y_id::Integer, y::CBMatrix, descent_step::Bool, model_maxviol::Real, new_indices::Union{<:CBIndexmatrix,Nothing} = nothing)
    current_weight = Ref{Float64}()
    @ccall libcb.cb_bundlelowranktrustregionprox_apply_variable_metric(self.data::Ptr{Cvoid}, (isnothing(groundset) ? C_NULL : groundset.data)::Ptr{Cvoid}, (isnothing(model) ? C_NULL : model.data)::Ptr{Cvoid}, aggr.data::Ptr{Cvoid}, y_id::Cint, y.data::Ptr{Cvoid}, descent_step::Cint, current_weight::Ref{Float64}, model_maxviol::Cdouble, (isnothing(new_indices) ? C_NULL : new_indices.data)::Ptr{Cvoid})::Cint
    return current_weight[]
end

@doc raw"""
    cb_add_variable_metric!(self::CBBundleLowRankTrustRegionProx, diagH::CBMatrix, vecH::CBMatrix)

see BundleProxObject::add_dynamic_scaling()
"""
cb_add_variable_metric!(self::CBBundleLowRankTrustRegionProx, diagH::CBMatrix, vecH::CBMatrix) = @ccall libcb.cb_bundlelowranktrustregionprox_add_variable_metric(self.data::Ptr{Cvoid}, diagH.data::Ptr{Cvoid}, vecH.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_push_aft!(self::CBBundleLowRankTrustRegionProx, aft::Union{<:CBAffineFunctionTransformation,Nothing})

see  BundleProxObject::push_aft();
"""
cb_push_aft!(self::CBBundleLowRankTrustRegionProx, aft::Union{<:CBAffineFunctionTransformation,Nothing}) = @ccall libcb.cb_bundlelowranktrustregionprox_push_aft(self.data::Ptr{Cvoid}, (isnothing(aft) ? C_NULL : aft.data)::Ptr{Cvoid})::Cint

@doc raw"""
    cb_pop_aft!(self::CBBundleLowRankTrustRegionProx)

see  BundleProxObject::pop_aft();
"""
cb_pop_aft!(self::CBBundleLowRankTrustRegionProx) = @ccall libcb.cb_bundlelowranktrustregionprox_pop_aft(self.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_mfile_data(self::CBBundleLowRankTrustRegionProx)

output the description of the prox term in mfile-suitable format
"""
cb_mfile_data(self::CBBundleLowRankTrustRegionProx) = @ccall libcb.cb_bundlelowranktrustregionprox_mfile_data(self.data::Ptr{Cvoid})::Cint

