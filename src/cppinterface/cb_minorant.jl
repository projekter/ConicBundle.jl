@doc raw"""
    CBMinorant(offset::Real, subg::CBDVector, primal::Union{<:CBPrimalData,Nothing} = nothing, offset_at_origin::Bool = false)

* @brief Minorant constructor for a dense subgradient specified via a DVector

      Suppose in the evaluation of your oracle at the current point \f$y\f$ you determine a subgradient \f$s\f$ for the subgradient inequality

      \f[ f(z)\ge f(y)+\langle s,z-y\rangle\quad\forall z\in\mathbf{R}^m, \f]

      then use \f$f(y)\f$ for @a offset and \f$s\f$ in the form of a DVector in @a subg and keep the default value @a offset_at_origin == false.

      If, on the other hand, your oracle implements a support function for some compact set \f$\mathcal{X}\subset\mathbf{R}^m\f$ like

      \f[ f(y) = \max_{x\in\mathcal{X}} x^\top y\f]

      then it is actually more efficient to return 0 for @a offset,
      a maximizing \f$x\f$ in @ subg and to put @a offset_at_origin = true.

      You may also pass over a PrimalData object in primal that will be then be
      owned and later deleted by the minorant and aggregated along with the minornat
     
"""
CBMinorant(offset::Real, subg::CBDVector, primal::Union{<:CBPrimalData,Nothing} = nothing, offset_at_origin::Bool = false) = CBMinorant(@ccall libcb.cb_minorant_new(offset::Cdouble, subg.data::Ptr{Cvoid}, (isnothing(primal) ? C_NULL : primal.data)::Ptr{Cvoid}, offset_at_origin::Cint)::Ptr{Cvoid})

@doc raw"""
    CBMinorant(offset::Real, subg_val::CBDVector, subg_ind::CBIVector, primal::Union{<:CBPrimalData,Nothing} = nothing, offset_at_origin::Bool = false)

* @brief Minorant constructor for a sparse subgradient where the nonzero values are specified by DVector and the corresponding indices by an IVector

      Suppose in the evaluation of your oracle at the current point \f$y\f$ you determine a subgradient \f$s\f$ with few nonzeros for the subgradient inequality

      \f[ f(z)\ge f(y)+\langle s,z-y\rangle\quad\forall z\in\mathbf{R}^m, \f]

      then use \f$f(y)\f$ for @a offset and pass \f$s\f$ by giving the nonzeros in the DVector @a subg_val, the corresponding indices in a IVector of the same length in @a subg_ind and keep the default value @a offset_at_origin == false.

      If, on the other hand, your oracle implements a support function for some compact set \f$\mathcal{X}\subset\mathbf{R}^m\f$ like

      \f[ f(y) = \max_{x\in\mathcal{X}} x^\top y\f]

      then it is actually more efficient to return 0 for @a offset,
      a sparse maximizing \f$x\f$ in @a subg_val and @a subg_ind  and to put @a offset_at_origin = true.

      You may also pass over a PrimalData object in primal that will be then be
      owned and later deleted by the minorant and aggregated along with the minornat
     
"""
CBMinorant(offset::Real, subg_val::CBDVector, subg_ind::CBIVector, primal::Union{<:CBPrimalData,Nothing} = nothing, offset_at_origin::Bool = false) = CBMinorant(@ccall libcb.cb_minorant_new2(offset::Cdouble, subg_val.data::Ptr{Cvoid}, subg_ind.data::Ptr{Cvoid}, (isnothing(primal) ? C_NULL : primal.data)::Ptr{Cvoid}, offset_at_origin::Cint)::Ptr{Cvoid})

@doc raw"""
    CBMinorant(offset_at_origin::Bool = true, offset::Real = 0., n_elementes::Integer = 0, coeffs::Union{<:AbstractVector{Real},Nothing} = nothing, indices::Union{<:AbstractVector{Integer},Nothing} = nothing, scale_val::Real = 1., primal::Union{<:CBPrimalData,Nothing} = nothing)

* @brief default initializes a zero minorant, see full explanation otherwise; NOTE: if the offset supplied by the minorant refers to the evaluation point (i.e., if it is the function value at the point of evaluation), set @a offset_at_origin = false !

    In many applications, in particular for Lagrangian relaxation,
    the offset is more naturally given for the origin, and giving it this way
    also reduces computational cost a bit, so @a offset_at_origin = true
    is the suggested default. If, however, the minorant arises from a subgradient
    inequality by evaluation in the current point, you may as well give the function
    value as offset directly and set @a offset_at_origin=false;

    The data specifying the minorant may be set here or (part of) it may be
    entered/added later. The meaning of the other parameters is as follows.

    @a offset gives the constant value (if offset_at_origin = false then relative to
    the evaluation point)

    @a n_elements gives the number of elements specified by @a coeffs
    (and possibly indices), but if @a coeffs == NULL it just asks to reserve space
    for that many coefficients

    If @a coeffs is not NULL, it points to an array of doubles of size
    at least @a n_elements. If @a indices == NULL then coeff[i] belongs to position
    i for i=0 to n_elements-1, otherwise coeff[i] belongs to position indices[i].
    All data is copied, the arrays are not modified, not used later and not
    deleted here.

    The entire input data (offset and coefficients) is multplied by @a scale_val
    to give the final minorant (scale_val is not memorized internally but executed
    immediately).

    If the minorant arises from @a PrimalData and the primal data should be aggregated
    along, it may be entered here or in the routine Minorant::set_primal(). The
    object pointed to is then owned by Minorant and will be deleted by Minorant on its
    destruction.

    
"""
CBMinorant(offset_at_origin::Bool = true, offset::Real = 0., n_elementes::Integer = 0, coeffs::Union{<:AbstractVector{Real},Nothing} = nothing, indices::Union{<:AbstractVector{Integer},Nothing} = nothing, scale_val::Real = 1., primal::Union{<:CBPrimalData,Nothing} = nothing) = GC.@preserve coeffs indices begin
    (LinearAlgebra.chkstride1(indices); LinearAlgebra.chkstride1(coeffs); CBMinorant(@ccall libcb.cb_minorant_new3(offset_at_origin::Cint, offset::Cdouble, n_elementes::Cint, coeffs::Ptr{Cdouble}, indices::Ptr{Cint}, scale_val::Cdouble, (isnothing(primal) ? C_NULL : primal.data)::Ptr{Cvoid})::Ptr{Cvoid}))
end

@doc raw"""
    CBMinorant(mnrt::Union{<:CBMinorant,Nothing}, factor::Real = 1., with_primal::Bool = false)

they main purpose of this constructor is to allow easy cloning for derived classes
"""
CBMinorant(mnrt::Union{<:CBMinorant,Nothing}, factor::Real = 1., with_primal::Bool = false) = CBMinorant(@ccall libcb.cb_minorant_new4((isnothing(mnrt) ? C_NULL : mnrt.data)::Ptr{Cvoid}, factor::Cdouble, with_primal::Cint)::Ptr{Cvoid})

@doc raw"""
    cb_offset(self::CBMinorant)

returns the current offset value
"""
cb_offset(self::CBMinorant) = @ccall libcb.cb_minorant_offset(self.data::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_add_offset!(self::CBMinorant, value::Real)

adds this value to the current offset value
"""
cb_add_offset!(self::CBMinorant, value::Real) = @ccall libcb.cb_minorant_add_offset(self.data::Ptr{Cvoid}, value::Cdouble)::Cint

@doc raw"""
    cb_offset_gives_value_at_origin!(self::CBMinorant)

allows to specify/modify whether the offset refers to origin or to the point of evaluation at which this minorant was supplied
"""
cb_offset_gives_value_at_origin!(self::CBMinorant) = @ccall libcb.cb_minorant_offset_gives_value_at_origin(self.data::Ptr{Cvoid})::Ptr{Bool}

@doc raw"""
    cb_offset_gives_value_at_origin(self::CBMinorant)

true if the offset refers to origin and false, if the offset refers to the value of the  minorant in the point of the oracle evaluation at which this minorant was supplied
"""
cb_offset_gives_value_at_origin(self::CBMinorant) = Bool(@ccall libcb.cb_minorant_offset_gives_value_at_origin2(self.data::Ptr{Cvoid})::Cint)

@doc raw"""
    cb_coeff!(self::CBMinorant, i::Integer)

returns the value of the coefficient in coordinate i
"""
cb_coeff!(self::CBMinorant, i::Integer) = @ccall libcb.cb_minorant_coeff(self.data::Ptr{Cvoid}, i::Cint)::Cdouble

@doc raw"""
    cb_add_coeff!(self::CBMinorant, i::Integer, value::Real)

adds the value to the coefficient in coordinate i
"""
cb_add_coeff!(self::CBMinorant, i::Integer, value::Real) = @ccall libcb.cb_minorant_add_coeff(self.data::Ptr{Cvoid}, i::Cint, value::Cdouble)::Cint

@doc raw"""
    cb_add_coeffs!(self::CBMinorant, n_elements::Integer, values::Union{<:AbstractVector{Real},Nothing}, factor::Real = 1., start_pos::Integer = 0)

adds the n values (possibly multiplied by factor) to consecutive coefficients starting at start_pos (by default 0); this always converts the minroant into a dense vector first and adds then
"""
cb_add_coeffs!(self::CBMinorant, n_elements::Integer, values::Union{<:AbstractVector{Real},Nothing}, factor::Real = 1., start_pos::Integer = 0) = GC.@preserve values begin
    (LinearAlgebra.chkstride1(values); @ccall libcb.cb_minorant_add_coeffs(self.data::Ptr{Cvoid}, n_elements::Cint, values::Ptr{Cdouble}, factor::Cdouble, start_pos::Cint)::Cint)
end

@doc raw"""
    cb_add_coeffs!(self::CBMinorant, n_elements::Integer, values::Union{<:AbstractVector{Real},Nothing}, indices::Union{<:AbstractVector{Integer},Nothing}, factor::Real = 1.)

adds the n values (possibly multiplied by factor) to the coefficients indicated by indices (if zero, this calls the other add_coeffs for the consecutive version)
"""
cb_add_coeffs!(self::CBMinorant, n_elements::Integer, values::Union{<:AbstractVector{Real},Nothing}, indices::Union{<:AbstractVector{Integer},Nothing}, factor::Real = 1.) = GC.@preserve values indices begin
    (LinearAlgebra.chkstride1(indices); LinearAlgebra.chkstride1(values); @ccall libcb.cb_minorant_add_coeffs2(self.data::Ptr{Cvoid}, n_elements::Cint, values::Ptr{Cdouble}, indices::Ptr{Cint}, factor::Cdouble)::Cint)
end

@doc raw"""
    cb_sparsify!(self::CBMinorant, tol::Real = CB_minorant_zero_tolerance, sparsity_ratio::Real = CB_minorant_sparsity_ratio)

converts to sparse format with zeros of absolut value at most tol*(fabs(offset)+1) if the given ratio of elements is zero in relation to the maximum nonzero index
"""
cb_sparsify!(self::CBMinorant, tol::Real = CB_minorant_zero_tolerance, sparsity_ratio::Real = CB_minorant_sparsity_ratio) = @ccall libcb.cb_minorant_sparsify(self.data::Ptr{Cvoid}, tol::Cdouble, sparsity_ratio::Cdouble)::Cint

@doc raw"""
    cb_nonzeros!(self::CBMinorant)

returns the number of nonzero coefficients
"""
cb_nonzeros!(self::CBMinorant) = @ccall libcb.cb_minorant_nonzeros(self.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_get_dense_coeff_store!(self::CBMinorant, n_elements::Integer)

*@brief If the returned pointer is not NULL it gives direct access to the array of current coefficient values with indices 0 up to n_elements-1;

      If the return value is NULL, the representation may not be
      available or access to the store may not be granted; in this case
      other routines like get_coeffs and add_coeffs have to be used.

      This routine is mainly intended for increasing efficiency in
      some internal computations; the validity of the pointer returned
      may get lost with any call to any other routine of this Minorant,
      so during manipulations of the stored values no other routines
      should be called. Needless to say, this routine should only be
      used by experts.
    
"""
cb_get_dense_coeff_store!(self::CBMinorant, n_elements::Integer) = @ccall libcb.cb_minorant_get_dense_coeff_store(self.data::Ptr{Cvoid}, n_elements::Cint)::Ptr{Cdouble}

@doc raw"""
    cb_reassign_coeffs!(self::CBMinorant, n_elements::Integer, map_to_old_coeffs::Union{<:AbstractVector{Integer},Nothing})

resorts (and deletes) coefficients so that afterwards it has n_elements and the new coeff(i) has the previous value of coeff(map_to_old_coeff(i))
"""
cb_reassign_coeffs!(self::CBMinorant, n_elements::Integer, map_to_old_coeffs::Union{<:AbstractVector{Integer},Nothing}) = GC.@preserve map_to_old_coeffs begin
    (LinearAlgebra.chkstride1(map_to_old_coeffs); @ccall libcb.cb_minorant_reassign_coeffs(self.data::Ptr{Cvoid}, n_elements::Cint, map_to_old_coeffs::Ptr{Cint})::Cint)
end

@doc raw"""
    cb_set_primal!(self::CBMinorant, param0::Union{<:CBPrimalData,Nothing})

if the minorant is generated by PrimalData and this should be aggregated along, insert a heap object of it here (see PrimalData)
"""
cb_set_primal!(self::CBMinorant, param0::Union{<:CBPrimalData,Nothing}) = @ccall libcb.cb_minorant_set_primal(self.data::Ptr{Cvoid}, (isnothing(param0) ? C_NULL : param0.data)::Ptr{Cvoid})::Cvoid

@doc raw"""
    cb_get_primal!(self::CBMinorant)

returns NULL if there is no primal data and otherwise a pointer to it
"""
cb_get_primal!(self::CBMinorant) = CBPrimalData(@ccall libcb.cb_minorant_get_primal(self.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_get_primal(self::CBMinorant)

returns NULL if there is no primal data and otherwise a pointer to it (const version)
"""
cb_get_primal(self::CBMinorant) = CBPrimalData(@ccall libcb.cb_minorant_get_primal2(self.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_clone_minorant(self::CBMinorant, factor::Real = 1., with_primal::Bool = true)

generates a full copy (multiplied by factor) on the heap and returns a pointer to it (it also includes a clone of the primal data if with_primal==true, otherwise the copy will have no primal data)
"""
cb_clone_minorant(self::CBMinorant, factor::Real = 1., with_primal::Bool = true) = CBMinorant(@ccall libcb.cb_minorant_clone_minorant(self.data::Ptr{Cvoid}, factor::Cdouble, with_primal::Cint)::Ptr{Cvoid})

@doc raw"""
    cb_aggregate!(self::CBMinorant, minorant::CBMinorant, factor::Real = 1.)

adds factor*minorant to this and does this also for the primal if it is availabe
"""
cb_aggregate!(self::CBMinorant, minorant::CBMinorant, factor::Real = 1.) = @ccall libcb.cb_minorant_aggregate(self.data::Ptr{Cvoid}, minorant.data::Ptr{Cvoid}, factor::Cdouble)::Cint

@doc raw"""
    cb_number_aggregated(self::CBMinorant)

returns the number of minorants aggregated in this one, value 1 thus means not aggregated
"""
cb_number_aggregated(self::CBMinorant) = @ccall libcb.cb_minorant_number_aggregated(self.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_scale_minorant!(self::CBMinorant, scale_val::Real)

mutliply offset and coefficients (and PrimalData, if given) by scale_val
"""
cb_scale_minorant!(self::CBMinorant, scale_val::Real) = @ccall libcb.cb_minorant_scale_minorant(self.data::Ptr{Cvoid}, scale_val::Cdouble)::Cint

@doc raw"""
    cb_norm_squared(self::CBMinorant)

return the squared Euclidean norm
"""
cb_norm_squared(self::CBMinorant) = @ccall libcb.cb_minorant_norm_squared(self.data::Ptr{Cvoid})::Cdouble

