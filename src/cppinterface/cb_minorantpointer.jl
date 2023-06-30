@doc raw"""
    Base.copy!(self::CBMinorantPointer, mp::CBMinorantPointer)

calls new_data with the data of mp and copies the modification_id
"""
Base.copy!(self::CBMinorantPointer, mp::CBMinorantPointer) = (@ccall libcb.cb_minorantpointer_assign(self.data::Ptr{Cvoid}, mp.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

@doc raw"""
    cb_init!(self::CBMinorantPointer, mp::CBMinorantPointer, factor::Real = 1., enforce_copy::Bool = false)

if factor!=1 it generates another MinorantUseData referring to the one of mp, otherwise it simply uses the same one (empty stays empty)
"""
cb_init!(self::CBMinorantPointer, mp::CBMinorantPointer, factor::Real = 1., enforce_copy::Bool = false) = @ccall libcb.cb_minorantpointer_init(self.data::Ptr{Cvoid}, mp.data::Ptr{Cvoid}, factor::Cdouble, enforce_copy::Cint)::Cvoid

@doc raw"""
    cb_init!(self::CBMinorantPointer, mp::Union{<:CBMinorant,Nothing}, modification_id::Integer = -1, factor::Real = 1.)

if mp==0 it becomes empty, otherwise it creates and then points to a MinorantUseData for holding mp with this modification_id and factor
"""
cb_init!(self::CBMinorantPointer, mp::Union{<:CBMinorant,Nothing}, modification_id::Integer = -1, factor::Real = 1.) = @ccall libcb.cb_minorantpointer_init2(self.data::Ptr{Cvoid}, (isnothing(mp) ? C_NULL : mp.data)::Ptr{Cvoid}, modification_id::Cint, factor::Cdouble)::Cvoid

@doc raw"""
    CBMinorantPointer()

declares the pointer _empty_
"""
CBMinorantPointer() = CBMinorantPointer(@ccall libcb.cb_minorantpointer_new()::Ptr{Cvoid})

@doc raw"""
    CBMinorantPointer(in_md::Union{<:CBMinorantUseData,Nothing})

initialize this to point to in_md
"""
CBMinorantPointer(in_md::Union{<:CBMinorantUseData,Nothing}) = CBMinorantPointer(@ccall libcb.cb_minorantpointer_new2((isnothing(in_md) ? C_NULL : in_md.data)::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    CBMinorantPointer(mp::CBMinorantPointer)

calls new_data
"""
CBMinorantPointer(mp::CBMinorantPointer) = CBMinorantPointer(@ccall libcb.cb_minorantpointer_new3(mp.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    CBMinorantPointer(mp::CBMinorantPointer, factor::Real)

calls init(const MinorantPointer&,CH_Matrix_Classes::Real)
"""
CBMinorantPointer(mp::CBMinorantPointer, factor::Real) = CBMinorantPointer(@ccall libcb.cb_minorantpointer_new4(mp.data::Ptr{Cvoid}, factor::Cdouble)::Ptr{Cvoid})

@doc raw"""
    CBMinorantPointer(mnrt::Union{<:CBMinorant,Nothing}, modification_id::Integer, factor::Real = 1.)

calls init(Minorant*,CH_Matrix_Classes::Integer,CH_Matrix_Classes::Real)
"""
CBMinorantPointer(mnrt::Union{<:CBMinorant,Nothing}, modification_id::Integer, factor::Real = 1.) = CBMinorantPointer(@ccall libcb.cb_minorantpointer_new5((isnothing(mnrt) ? C_NULL : mnrt.data)::Ptr{Cvoid}, modification_id::Cint, factor::Cdouble)::Ptr{Cvoid})

@doc raw"""
    cb_clear!(self::CBMinorantPointer)

afterwards the pointer is _empty_ (calls delete_data())
"""
cb_clear!(self::CBMinorantPointer) = @ccall libcb.cb_minorantpointer_clear(self.data::Ptr{Cvoid})::Cvoid

@doc raw"""
    cb_synchronize_ids!(self::CBMinorantPointer, new_modification_id::Integer, new_center_id::Integer, old_center_id::Integer, new_cand_id::Integer, old_cand_id::Integer, new_prex_id::Integer = 0)

if not empty it sets the modification_id to its new id and reinitializes the evaluation map
"""
cb_synchronize_ids!(self::CBMinorantPointer, new_modification_id::Integer, new_center_id::Integer, old_center_id::Integer, new_cand_id::Integer, old_cand_id::Integer, new_prex_id::Integer = 0) = @ccall libcb.cb_minorantpointer_synchronize_ids(self.data::Ptr{Cvoid}, new_modification_id::Cint, new_center_id::Cint, old_center_id::Cint, new_cand_id::Cint, old_cand_id::Cint, new_prex_id::Cint)::Cint

@doc raw"""
    cb_empty(self::CBMinorantPointer)

returns true if the pointer is _empty_
"""
cb_empty(self::CBMinorantPointer) = Bool(@ccall libcb.cb_minorantpointer_empty(self.data::Ptr{Cvoid})::Cint)

@doc raw"""
    cb_valid(self::CBMinorantPointer)

returns true if the pointer is not empty and the data is valid
"""
cb_valid(self::CBMinorantPointer) = Bool(@ccall libcb.cb_minorantpointer_valid(self.data::Ptr{Cvoid})::Cint)

@doc raw"""
    cb_zero(self::CBMinorantPointer)

returns true if the pointer is not empty but all entrys (also the offset) are zero
"""
cb_zero(self::CBMinorantPointer) = Bool(@ccall libcb.cb_minorantpointer_zero(self.data::Ptr{Cvoid})::Cint)

@doc raw"""
    cb_aggregate(self::CBMinorantPointer)

*@brief returns ture if it points to a combination of minorants from the same function;

       A minorant is considered an aggregate if it is obtained from an aggregate
       or it was computed by a call to aggregate(MinorantBundle&,Real);
       It is not an aggregate if it is the sum of non aggregate minorants from
       different functions.

       An _empty_ pointer is not an aggregate.
    
"""
cb_aggregate(self::CBMinorantPointer) = Bool(@ccall libcb.cb_minorantpointer_aggregate(self.data::Ptr{Cvoid})::Cint)

@doc raw"""
    cb_one_user(self::CBMinorantPointer)

returns true if not valid or this is the only active pointer to the minorant
"""
cb_one_user(self::CBMinorantPointer) = Bool(@ccall libcb.cb_minorantpointer_one_user(self.data::Ptr{Cvoid})::Cint)

@doc raw"""
    cb_get_minorant(self::CBMinorantPointer)

returns the Minorant *this points to or 0 if _empty_
"""
cb_get_minorant(self::CBMinorantPointer) = CBMinorant(@ccall libcb.cb_minorantpointer_get_minorant(self.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_get_primal!(self::CBMinorantPointer)

returns the primal of the Minorant *this points to or 0 if _empty_ (first carrying out any pending scalings on the minorant)
"""
cb_get_primal!(self::CBMinorantPointer) = CBPrimalData(@ccall libcb.cb_minorantpointer_get_primal(self.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_offset(self::CBMinorantPointer)

returns the offset of the minorant (including the internal scalings) or CB_minus_infinity if _empty_
"""
cb_offset(self::CBMinorantPointer) = @ccall libcb.cb_minorantpointer_offset(self.data::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_coeff(self::CBMinorantPointer, i::Integer)

returns coefficient i of the minorant (including the internal scalings) or 0. if _empty_
"""
cb_coeff(self::CBMinorantPointer, i::Integer) = @ccall libcb.cb_minorantpointer_coeff(self.data::Ptr{Cvoid}, i::Cint)::Cdouble

@doc raw"""
    cb_nonzeros(self::CBMinorantPointer)

returns the number of nonzero coefficients
"""
cb_nonzeros(self::CBMinorantPointer) = @ccall libcb.cb_minorantpointer_nonzeros(self.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_scale!(self::CBMinorantPointer, val::Real)

multiply the MinorantPointer by val (this is an external factor for the minorant and possibly its primal information, but it will be used in aggreagation and when retrieving the approximate primal)
"""
cb_scale!(self::CBMinorantPointer, val::Real) = @ccall libcb.cb_minorantpointer_scale(self.data::Ptr{Cvoid}, val::Cdouble)::Cint

@doc raw"""
    cb_call_primal_extender!(self::CBMinorantPointer, prex::CBPrimalExtender, prex_id::Integer)

executes prex on the primal if its id is smaller then prex_id; returns != if this fails or not primal data is available
"""
cb_call_primal_extender!(self::CBMinorantPointer, prex::CBPrimalExtender, prex_id::Integer) = @ccall libcb.cb_minorantpointer_call_primal_extender(self.data::Ptr{Cvoid}, prex.data::Ptr{Cvoid}, prex_id::Cint)::Cint

@doc raw"""
    cb_apply_modification!(self::CBMinorantPointer, gsmdf::CBGroundsetModification, mod_id::Integer, mex::Union{<:CBMinorantExtender,Nothing}, apply_gsmdf_costs::Bool = false)

if valid and the local modification_id is smaller than mod_id, it modifies the coefficients as described by GroundsetModification; if successful, the mod_id is assigend and the function returns 0; otherwise it is invalidated and returns !=0; the costs of gsmdf are used only if apply_gsmdf_costs=false in whic case mex must be NULL
"""
cb_apply_modification!(self::CBMinorantPointer, gsmdf::CBGroundsetModification, mod_id::Integer, mex::Union{<:CBMinorantExtender,Nothing}, apply_gsmdf_costs::Bool = false) = @ccall libcb.cb_minorantpointer_apply_modification(self.data::Ptr{Cvoid}, gsmdf.data::Ptr{Cvoid}, mod_id::Cint, (isnothing(mex) ? C_NULL : mex.data)::Ptr{Cvoid}, apply_gsmdf_costs::Cint)::Cint

@doc raw"""
    cb_evaluate(self::CBMinorantPointer, yid::Integer, y::CBMatrix, with_constant::Bool = true)

negative ids are allowed and indicate there is no need to memorize this result, returns CB_minus_infinity if empty, otherwise offset+ip(minorant,y)
"""
cb_evaluate(self::CBMinorantPointer, yid::Integer, y::CBMatrix, with_constant::Bool = true) = @ccall libcb.cb_minorantpointer_evaluate(self.data::Ptr{Cvoid}, yid::Cint, y.data::Ptr{Cvoid}, with_constant::Cint)::Cdouble

@doc raw"""
    cb_add_offset!(self::CBMinorantPointer, offset::Real)

add the offset; if empty, initialize to offset with zero linear part, if not the only user (use_cnt>1), clone it first
"""
cb_add_offset!(self::CBMinorantPointer, offset::Real) = @ccall libcb.cb_minorantpointer_add_offset(self.data::Ptr{Cvoid}, offset::Cdouble)::Cint

@doc raw"""
    cb_get_minorant(self::CBMinorantPointer, mat::CBMatrix, column::Integer, alpha::Real = 1., add::Bool = false, skip_fixed::Union{<:CBIndexmatrix,Nothing} = nothing, fixed_vals::Union{<:CBMatrix,Nothing} = nothing)

*@brief store/add the minorant in offset and a matrix column, possibly skipping indices skip_fixed. The coordinats of the latter ones are multiplied by the given values or 0 and added to offset

       The dimensions of the matrix must already fit the requirements on input.
       If add is false (as by default), the full length of the column
       is initialized to the gradient and filled up with zeros where needed.
       skip_fixed and fixed_vals must both be given or both not be given.
       If they are both given, both must have the same length.
       If skip_fixed is given, it must have strictly increasing indices.
       The skip part of the routine is used in implementations of BundleScaling::get_QP_costs.
     
"""
function cb_get_minorant(self::CBMinorantPointer, mat::CBMatrix, column::Integer, alpha::Real = 1., add::Bool = false, skip_fixed::Union{<:CBIndexmatrix,Nothing} = nothing, fixed_vals::Union{<:CBMatrix,Nothing} = nothing)
    offset = Ref{Float64}()
    @ccall libcb.cb_minorantpointer_get_minorant2(self.data::Ptr{Cvoid}, offset::Ref{Float64}, mat.data::Ptr{Cvoid}, column::Cint, alpha::Cdouble, add::Cint, (isnothing(skip_fixed) ? C_NULL : skip_fixed.data)::Ptr{Cvoid}, (isnothing(fixed_vals) ? C_NULL : fixed_vals.data)::Ptr{Cvoid})::Cint
    return offset[]
end

@doc raw"""
    cb_get_minorant(self::CBMinorantPointer, mp::CBMinorantPointer, alpha::Real = 1.)

*@brief store/add the minorant in/to mp, possibly scaled by alpha

       If mp is empty, store it there, if mp is not empty, add it.

       if *this is empty, it causes an error.
     
"""
cb_get_minorant(self::CBMinorantPointer, mp::CBMinorantPointer, alpha::Real = 1.) = @ccall libcb.cb_minorantpointer_get_minorant3(self.data::Ptr{Cvoid}, mp.data::Ptr{Cvoid}, alpha::Cdouble)::Cint

@doc raw"""
    cb_get_minorant(self::CBMinorantPointer, mp::CBMinorantPointer, alpha::Real, sp::Union{<:CBSparsemat,Nothing}, provided_row_indices::Union{<:CBIndexmatrix,Nothing} = nothing, needed_col_indices::Union{<:CBIndexmatrix,Nothing} = nothing, enforce_copy::Bool = false)

*@brief store/add the minorant in/to mp, possibly scaled by alpha and transformed by sp, which possibly requires only the indices of provided_row_indices to compute possibly only the indices in needed_col_indices

       If mp is empty, store it there, if mp is not empty, add it.

       if sp==0, it is treated as the identity

       If provided_row_indices or needed_col_indices is given, its indices must
       be in strictly increasing order.  For sp==0 both must coincide and will
       probably be ignored by just returning a scaled reference to *this in mp.

       if *this is empty, it causes an error.
     
"""
cb_get_minorant(self::CBMinorantPointer, mp::CBMinorantPointer, alpha::Real, sp::Union{<:CBSparsemat,Nothing}, provided_row_indices::Union{<:CBIndexmatrix,Nothing} = nothing, needed_col_indices::Union{<:CBIndexmatrix,Nothing} = nothing, enforce_copy::Bool = false) = @ccall libcb.cb_minorantpointer_get_minorant4(self.data::Ptr{Cvoid}, mp.data::Ptr{Cvoid}, alpha::Cdouble, (isnothing(sp) ? C_NULL : sp.data)::Ptr{Cvoid}, (isnothing(provided_row_indices) ? C_NULL : provided_row_indices.data)::Ptr{Cvoid}, (isnothing(needed_col_indices) ? C_NULL : needed_col_indices.data)::Ptr{Cvoid}, enforce_copy::Cint)::Cint

@doc raw"""
    cb_aggregate!(self::CBMinorantPointer, minorants::CBMinorantBundle, coeff::CBMatrix, factor::Real = 1.)

*@brief collect the aggregate as the nonnegative linear combination of the bundle minorants (here *this may not be part of the bundle!). If *this is not empty, add this aggregate to *this. Only aggregation carries along primal data!

       All minorants in the bundle need to be valid. Aggregating over no minorants
       is treated as a zero minorant without any primal information.
     
"""
cb_aggregate!(self::CBMinorantPointer, minorants::CBMinorantBundle, coeff::CBMatrix, factor::Real = 1.) = @ccall libcb.cb_minorantpointer_aggregate2(self.data::Ptr{Cvoid}, minorants.data::Ptr{Cvoid}, coeff.data::Ptr{Cvoid}, factor::Cdouble)::Cint

@doc raw"""
    cb_aggregate!(self::CBMinorantPointer, minorant::CBMinorantPointer, itsfactor::Real = 1.)

*@brief aggregate the minorant*itsfactor to *this, both must be valid for this. Only aggregation carries along primal data!
    
"""
cb_aggregate!(self::CBMinorantPointer, minorant::CBMinorantPointer, itsfactor::Real = 1.) = @ccall libcb.cb_minorantpointer_aggregate3(self.data::Ptr{Cvoid}, minorant.data::Ptr{Cvoid}, itsfactor::Cdouble)::Cint

@doc raw"""
    cb_ip(self::CBMinorantPointer, mp::CBMinorantPointer, skip_fixed::Union{<:CBIndexmatrix,Nothing} = nothing, ipdiag::Union{<:CBMatrix,Nothing} = nothing)

*@brief computes the inner product of the two minorants; if skip_fixed!=NULL the corrsponding indices are not considered, if ipdiag!=0 the inner product is taken with respect to this diagonal matrix, i.e.  sum_i mp(i)*(*this)(i)*(*ipdiag)(i)
     
"""
cb_ip(self::CBMinorantPointer, mp::CBMinorantPointer, skip_fixed::Union{<:CBIndexmatrix,Nothing} = nothing, ipdiag::Union{<:CBMatrix,Nothing} = nothing) = @ccall libcb.cb_minorantpointer_ip(self.data::Ptr{Cvoid}, mp.data::Ptr{Cvoid}, (isnothing(skip_fixed) ? C_NULL : skip_fixed.data)::Ptr{Cvoid}, (isnothing(ipdiag) ? C_NULL : ipdiag.data)::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_ip(self::CBMinorantPointer, m::CBMatrix, ipdiag::Union{<:CBMatrix,Nothing} = nothing, startindex_m::Integer = 0)

*@brief computes the inner product with m; if ipdiag!=0 the inner product is taken with respect to this diagonal matrix, i.e.  sum_i mp(startindex_m+i)*(*this)(i)*(*ipdiag)(i) 
"""
cb_ip(self::CBMinorantPointer, m::CBMatrix, ipdiag::Union{<:CBMatrix,Nothing} = nothing, startindex_m::Integer = 0) = @ccall libcb.cb_minorantpointer_ip2(self.data::Ptr{Cvoid}, m.data::Ptr{Cvoid}, (isnothing(ipdiag) ? C_NULL : ipdiag.data)::Ptr{Cvoid}, startindex_m::Cint)::Cdouble

@doc raw"""
    Base.:<(self::CBMinorantPointer, mp::CBMinorantPointer)

they are equal if they point to the same object, compares the addresses of this objects
"""
Base.:<(self::CBMinorantPointer, mp::CBMinorantPointer) = Bool(@ccall libcb.cb_minorantpointer_new_less(self.data::Ptr{Cvoid}, mp.data::Ptr{Cvoid})::Bool)

@doc raw"""
    Base.:(==)(self::CBMinorantPointer, mp::CBMinorantPointer)

they are equal if they point to the same object, compares the addresses of this objects
"""
Base.:(==)(self::CBMinorantPointer, mp::CBMinorantPointer) = Bool(@ccall libcb.cb_minorantpointer_new_equal(self.data::Ptr{Cvoid}, mp.data::Ptr{Cvoid})::Bool)

@doc raw"""
    Base.:>(self::CBMinorantPointer, mp::CBMinorantPointer)

they are equal if they point to the same object, compares the addresses of this objects
"""
Base.:>(self::CBMinorantPointer, mp::CBMinorantPointer) = Bool(@ccall libcb.cb_minorantpointer_new_greater(self.data::Ptr{Cvoid}, mp.data::Ptr{Cvoid})::Bool)

@doc raw"""
    cb_equals(self::CBMinorantPointer, mp::CBMinorantPointer, tol::Real = 1e-10)

they are equal if they point to the same object or are both 0. If not, they differ if their matrix representations differ; if not, they differ if the entries differ by at least tol*(1.+fabs(this->offset())
"""
cb_equals(self::CBMinorantPointer, mp::CBMinorantPointer, tol::Real = 1e-10) = Bool(@ccall libcb.cb_minorantpointer_equals(self.data::Ptr{Cvoid}, mp.data::Ptr{Cvoid}, tol::Cdouble)::Cint)

@doc raw"""
    cb_norm_squared(self::CBMinorantPointer, D::Union{<:CBMatrix,Nothing} = nothing)

Compute the norm squared of this for the given diagonal matrix D (identity if not given), i.e. \f$\|(*this)\|^2_{D}\f$
"""
cb_norm_squared(self::CBMinorantPointer, D::Union{<:CBMatrix,Nothing} = nothing) = @ccall libcb.cb_minorantpointer_norm_squared(self.data::Ptr{Cvoid}, (isnothing(D) ? C_NULL : D.data)::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_dual_norm_squared(self::CBMinorantPointer, D::Union{<:CBMatrix,Nothing} = nothing)

Compute the dual norm squared of this for the given diagonal matrix D (identity if not given), i.e. \f$\|(*this)\|^2_{D^{-1}}\f$
"""
cb_dual_norm_squared(self::CBMinorantPointer, D::Union{<:CBMatrix,Nothing} = nothing) = @ccall libcb.cb_minorantpointer_dual_norm_squared(self.data::Ptr{Cvoid}, (isnothing(D) ? C_NULL : D.data)::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_left_genmult(self::CBMinorantPointer, B::CBMatrix, C::CBMatrix, alpha::Real = 1., beta::Real = 0., thistrans::Integer = 0, btrans::Integer = 0, thisindex::Integer = 0)

computes and returns C=alpha*(*this)*B+beta*C where B and *this may be transposed and *this is considered to be a column with thisindex in a bigger matrix
"""
cb_left_genmult(self::CBMinorantPointer, B::CBMatrix, C::CBMatrix, alpha::Real = 1., beta::Real = 0., thistrans::Integer = 0, btrans::Integer = 0, thisindex::Integer = 0) = (@ccall libcb.cb_minorantpointer_left_genmult(self.data::Ptr{Cvoid}, B.data::Ptr{Cvoid}, C.data::Ptr{Cvoid}, alpha::Cdouble, beta::Cdouble, thistrans::Cint, btrans::Cint, thisindex::Cint)::Ptr{Cvoid}; return self)

@doc raw"""
    cb_right_genmult(self::CBMinorantPointer, A::CBMatrix, C::CBMatrix, alpha::Real = 1., beta::Real = 0., atrans::Integer = 0, thistrans::Integer = 0, thisindex::Integer = 0)

computes and returns C=alpha*A*(*this)+beta*C where A and *this may be transposed and *this is considered to be a column with thisindex in a bigger matrix
"""
cb_right_genmult(self::CBMinorantPointer, A::CBMatrix, C::CBMatrix, alpha::Real = 1., beta::Real = 0., atrans::Integer = 0, thistrans::Integer = 0, thisindex::Integer = 0) = (@ccall libcb.cb_minorantpointer_right_genmult(self.data::Ptr{Cvoid}, A.data::Ptr{Cvoid}, C.data::Ptr{Cvoid}, alpha::Cdouble, beta::Cdouble, atrans::Cint, thistrans::Cint, thisindex::Cint)::Ptr{Cvoid}; return self)

@doc raw"""
    cb_display(self::CBMinorantPointer, precision::Integer = 8)

output the Minorant in a nice format
"""
cb_display(self::CBMinorantPointer, precision::Integer = 8) = @ccall libcb.cb_minorantpointer_display(self.data::Ptr{Cvoid}, precision::Cint)::Cvoid

