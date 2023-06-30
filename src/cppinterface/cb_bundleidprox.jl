@doc raw"""
    CBBundleIdProx(d::Integer = 0, w::Real = 1., cbinc::Integer = -1)

initialize with dimension and weight
"""
CBBundleIdProx(d::Integer = 0, w::Real = 1., cbinc::Integer = -1) = CBBundleIdProx(@ccall libcb.cb_bundleidprox_new(d::Cint, w::Cdouble, cbinc::Cint)::Ptr{Cvoid})

@doc raw"""
    cb_set_weightu!(self::CBBundleIdProx, in_weightu::Real)

set the weight of the proximal term
"""
cb_set_weightu!(self::CBBundleIdProx, in_weightu::Real) = @ccall libcb.cb_bundleidprox_set_weightu(self.data::Ptr{Cvoid}, in_weightu::Cdouble)::Cvoid

@doc raw"""
    cb_get_weightu(self::CBBundleIdProx)

returns the current weight of the proximal term
"""
cb_get_weightu(self::CBBundleIdProx) = @ccall libcb.cb_bundleidprox_get_weightu(self.data::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_get_term_corr(self::CBBundleIdProx)

returns the correction factor for the termination criterion, here min(1,1/weight)
"""
cb_get_term_corr(self::CBBundleIdProx) = @ccall libcb.cb_bundleidprox_get_term_corr(self.data::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_norm_sqr(self::CBBundleIdProx, B::CBMatrix)

returns \f$\|B\|^2_H\f$
"""
cb_norm_sqr(self::CBBundleIdProx, B::CBMatrix) = @ccall libcb.cb_bundleidprox_norm_sqr(self.data::Ptr{Cvoid}, B.data::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_dnorm_sqr(self::CBBundleIdProx, B::CBMinorantPointer)

returns \f$\|B\|^2_{H^{-1}}\f$
"""
cb_dnorm_sqr(self::CBBundleIdProx, B::CBMinorantPointer) = @ccall libcb.cb_bundleidprox_dnorm_sqr(self.data::Ptr{Cvoid}, B.data::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_is_DLR(self::CBBundleIdProx)

return true if H is of the form diagonal matrix plus Gram matrix of a low rank matrix
"""
cb_is_DLR(self::CBBundleIdProx) = Bool(@ccall libcb.cb_bundleidprox_is_dlr(self.data::Ptr{Cvoid})::Cint)

@doc raw"""
    cb_add_H(self::CBBundleIdProx, big_sym::CBSymmatrix, start_index::Integer = 0)

add H to the dense symmetric matrix as a principal submatrix starting at position start_index
"""
cb_add_H(self::CBBundleIdProx, big_sym::CBSymmatrix, start_index::Integer = 0) = @ccall libcb.cb_bundleidprox_add_h(self.data::Ptr{Cvoid}, big_sym.data::Ptr{Cvoid}, start_index::Cint)::Cint

@doc raw"""
    cb_add_Hx(self::CBBundleIdProx, x::CBMatrix, outplusHx::CBMatrix, alpha::Real = 1.)

adds \f$alpha*Hx\f$ to outplusHx and returns this
"""
cb_add_Hx(self::CBBundleIdProx, x::CBMatrix, outplusHx::CBMatrix, alpha::Real = 1.) = (@ccall libcb.cb_bundleidprox_add_hx(self.data::Ptr{Cvoid}, x.data::Ptr{Cvoid}, outplusHx.data::Ptr{Cvoid}, alpha::Cdouble)::Ptr{Cvoid}; return self)

@doc raw"""
    cb_apply_Hinv(self::CBBundleIdProx, x::CBMatrix)

returns \f$H^{-1}x\f$
"""
cb_apply_Hinv(self::CBBundleIdProx, x::CBMatrix) = (@ccall libcb.cb_bundleidprox_apply_hinv(self.data::Ptr{Cvoid}, x.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

@doc raw"""
    cb_compute_QP_costs!(self::CBBundleIdProx, Q::CBSymmatrix, d::CBMatrix, constant_minorant::CBMinorantPointer, bundle::CBMinorantBundle, y::CBMatrix, groundset_minorant::CBMinorantPointer, yfixed::Union{<:CBIndexmatrix,Nothing})

computes the dual QP costs Q, d, and the constant offset to the bundle subproblem, see BundleProxObject::compute_QP_costs
"""
function cb_compute_QP_costs!(self::CBBundleIdProx, Q::CBSymmatrix, d::CBMatrix, constant_minorant::CBMinorantPointer, bundle::CBMinorantBundle, y::CBMatrix, groundset_minorant::CBMinorantPointer, yfixed::Union{<:CBIndexmatrix,Nothing})
    offset = Ref{Float64}()
    @ccall libcb.cb_bundleidprox_compute_qp_costs(self.data::Ptr{Cvoid}, Q.data::Ptr{Cvoid}, d.data::Ptr{Cvoid}, offset::Ref{Float64}, constant_minorant.data::Ptr{Cvoid}, bundle.data::Ptr{Cvoid}, y.data::Ptr{Cvoid}, groundset_minorant.data::Ptr{Cvoid}, (isnothing(yfixed) ? C_NULL : yfixed.data)::Ptr{Cvoid})::Cint
    return offset[]
end

@doc raw"""
    cb_update_QP_costs!(self::CBBundleIdProx, delta_Q::CBSymmatrix, delta_d::CBMatrix, constant_minorant::CBMinorantPointer, bundle::CBMinorantBundle, center_y::CBMatrix, groundset_minorant::CBMinorantPointer, delta_groundset_minorant::CBMinorantPointer, delta_index::CBIndexmatrix, yfixed::Union{<:CBIndexmatrix,Nothing})

updates the dual QP costs Q, d, and the constant offset to the bundle subproblem, see BundleProxObject::update_QP_costs
"""
function cb_update_QP_costs!(self::CBBundleIdProx, delta_Q::CBSymmatrix, delta_d::CBMatrix, constant_minorant::CBMinorantPointer, bundle::CBMinorantBundle, center_y::CBMatrix, groundset_minorant::CBMinorantPointer, delta_groundset_minorant::CBMinorantPointer, delta_index::CBIndexmatrix, yfixed::Union{<:CBIndexmatrix,Nothing})
    delta_offset = Ref{Float64}()
    @ccall libcb.cb_bundleidprox_update_qp_costs(self.data::Ptr{Cvoid}, delta_Q.data::Ptr{Cvoid}, delta_d.data::Ptr{Cvoid}, delta_offset::Ref{Float64}, constant_minorant.data::Ptr{Cvoid}, bundle.data::Ptr{Cvoid}, center_y.data::Ptr{Cvoid}, groundset_minorant.data::Ptr{Cvoid}, delta_groundset_minorant.data::Ptr{Cvoid}, delta_index.data::Ptr{Cvoid}, (isnothing(yfixed) ? C_NULL : yfixed.data)::Ptr{Cvoid})::Cint
    return delta_offset[]
end

@doc raw"""
    cb_apply_modification!(self::CBBundleIdProx, gsmdf::CBGroundsetModification)

when BundleSolver is called to modify the groundset it also calls this
"""
cb_apply_modification!(self::CBBundleIdProx, gsmdf::CBGroundsetModification) = @ccall libcb.cb_bundleidprox_apply_modification(self.data::Ptr{Cvoid}, gsmdf.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_projected_clone!(self::CBBundleIdProx, indices::CBIndexmatrix)

* @brief in order to allow for fixed variables, this generates a clone restricted to the given indices
   
"""
cb_projected_clone!(self::CBBundleIdProx, indices::CBIndexmatrix) = CBBundleProxObject(@ccall libcb.cb_bundleidprox_projected_clone(self.data::Ptr{Cvoid}, indices.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_mfile_data(self::CBBundleIdProx)

output the description of the prox term in mfile-suitable format
"""
cb_mfile_data(self::CBBundleIdProx) = @ccall libcb.cb_bundleidprox_mfile_data(self.data::Ptr{Cvoid})::Cint

