@doc raw"""
    cb_init!(self::CBAffineFunctionTransformation, fun_coeff::Real = 1., fun_offset::Real = 0., linear_cost::Union{<:CBMatrix,Nothing} = nothing, arg_offset::Union{<:CBMatrix,Nothing} = nothing, arg_trafo::Union{<:CBSparsemat,Nothing} = nothing, model_calls_delete::Bool = true)

sets the parameters of the transformation. The ownership of objects pointed to is passed to *this (they will be deleted here). If *this is entered into an AFTModel, model_calls_delete==true tells the AFTModel to delete this AffineFunctionTransformation at the end.
"""
cb_init!(self::CBAffineFunctionTransformation, fun_coeff::Real = 1., fun_offset::Real = 0., linear_cost::Union{<:CBMatrix,Nothing} = nothing, arg_offset::Union{<:CBMatrix,Nothing} = nothing, arg_trafo::Union{<:CBSparsemat,Nothing} = nothing, model_calls_delete::Bool = true) = @ccall libcb.cb_affinefunctiontransformation_init(self.data::Ptr{Cvoid}, fun_coeff::Cdouble, fun_offset::Cdouble, (isnothing(linear_cost) ? C_NULL : linear_cost.data)::Ptr{Cvoid}, (isnothing(arg_offset) ? C_NULL : arg_offset.data)::Ptr{Cvoid}, (isnothing(arg_trafo) ? C_NULL : arg_trafo.data)::Ptr{Cvoid}, model_calls_delete::Cint)::Cint

@doc raw"""
    CBAffineFunctionTransformation(in_fun_coeff::Real = 1., in_fun_offset::Real = 0., in_linear_cost::Union{<:CBMatrix,Nothing} = nothing, in_arg_offset::Union{<:CBMatrix,Nothing} = nothing, in_arg_trafo::Union{<:CBSparsemat,Nothing} = nothing, in_model_calls_delete::Bool = true, incr::Integer = -1)

calls init()
"""
CBAffineFunctionTransformation(in_fun_coeff::Real = 1., in_fun_offset::Real = 0., in_linear_cost::Union{<:CBMatrix,Nothing} = nothing, in_arg_offset::Union{<:CBMatrix,Nothing} = nothing, in_arg_trafo::Union{<:CBSparsemat,Nothing} = nothing, in_model_calls_delete::Bool = true, incr::Integer = -1) = CBAffineFunctionTransformation(@ccall libcb.cb_affinefunctiontransformation_new(in_fun_coeff::Cdouble, in_fun_offset::Cdouble, (isnothing(in_linear_cost) ? C_NULL : in_linear_cost.data)::Ptr{Cvoid}, (isnothing(in_arg_offset) ? C_NULL : in_arg_offset.data)::Ptr{Cvoid}, (isnothing(in_arg_trafo) ? C_NULL : in_arg_trafo.data)::Ptr{Cvoid}, in_model_calls_delete::Cint, incr::Cint)::Ptr{Cvoid})

@doc raw"""
    cb_get_model_calls_delete!(self::CBAffineFunctionTransformation)

retruns true if the model has to delete this
"""
cb_get_model_calls_delete!(self::CBAffineFunctionTransformation) = Bool(@ccall libcb.cb_affinefunctiontransformation_get_model_calls_delete(self.data::Ptr{Cvoid})::Cint)

@doc raw"""
    cb_set_model_calls_delete!(self::CBAffineFunctionTransformation, mcd::Bool)

set to true if the model has to delete this, to false if it is destructed elsewhere
"""
cb_set_model_calls_delete!(self::CBAffineFunctionTransformation, mcd::Bool) = @ccall libcb.cb_affinefunctiontransformation_set_model_calls_delete(self.data::Ptr{Cvoid}, mcd::Cint)::Cvoid

@doc raw"""
    cb_argument_changes(self::CBAffineFunctionTransformation)

returns true if not the identity
"""
cb_argument_changes(self::CBAffineFunctionTransformation) = Bool(@ccall libcb.cb_affinefunctiontransformation_argument_changes(self.data::Ptr{Cvoid})::Cint)

@doc raw"""
    cb_sparse_argument_changes(self::CBAffineFunctionTransformation)

* @brief returns true if arg_trafo influences at most two thirds of the
  entries of local_argument
     
"""
cb_sparse_argument_changes(self::CBAffineFunctionTransformation) = Bool(@ccall libcb.cb_affinefunctiontransformation_sparse_argument_changes(self.data::Ptr{Cvoid})::Cint)

@doc raw"""
    cb_scaled_index_subset(self::CBAffineFunctionTransformation, col_ind::Union{<:CBIndexmatrix,Nothing}, row_ind::Union{<:CBIndexmatrix,Nothing})

* @brief returns true if the transformation maps each index (in col_ind
  if !=0) onto at most one index and vice versa (out of row_ind if !=0)
     
"""
cb_scaled_index_subset(self::CBAffineFunctionTransformation, col_ind::Union{<:CBIndexmatrix,Nothing}, row_ind::Union{<:CBIndexmatrix,Nothing}) = Bool(@ccall libcb.cb_affinefunctiontransformation_scaled_index_subset(self.data::Ptr{Cvoid}, (isnothing(col_ind) ? C_NULL : col_ind.data)::Ptr{Cvoid}, (isnothing(row_ind) ? C_NULL : row_ind.data)::Ptr{Cvoid})::Cint)

@doc raw"""
    cb_scaled_index(self::CBAffineFunctionTransformation, index::Integer)

* @brief returns false if index is mapped to more than one index,
  otherwise true with mapped_index==-1 if mapped to zero, else >=0 and
  @a coeff gives the coefficient
     
"""
function cb_scaled_index(self::CBAffineFunctionTransformation, index::Integer)
    coeff = Ref{Float64}()
    mapped_index = Ref{Int}()
    Bool(@ccall libcb.cb_affinefunctiontransformation_scaled_index(self.data::Ptr{Cvoid}, mapped_index::Ref{Int}, coeff::Ref{Float64}, index::Cint)::Cint)
    return mapped_index[], coeff[]
end

@doc raw"""
    cb_get_fun_coeff(self::CBAffineFunctionTransformation)

returns the factor for the function
"""
cb_get_fun_coeff(self::CBAffineFunctionTransformation) = @ccall libcb.cb_affinefunctiontransformation_get_fun_coeff(self.data::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_set_fun_coeff!(self::CBAffineFunctionTransformation)

allows to set the factor for the function
"""
cb_set_fun_coeff!(self::CBAffineFunctionTransformation) = @ccall libcb.cb_affinefunctiontransformation_set_fun_coeff(self.data::Ptr{Cvoid})::Ptr{Cdouble}

@doc raw"""
    cb_get_fun_offset(self::CBAffineFunctionTransformation)

returns the constant offset for the funciton
"""
cb_get_fun_offset(self::CBAffineFunctionTransformation) = @ccall libcb.cb_affinefunctiontransformation_get_fun_offset(self.data::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_set_fun_offset!(self::CBAffineFunctionTransformation)

allows to set the constant offset for the funciton
"""
cb_set_fun_offset!(self::CBAffineFunctionTransformation) = @ccall libcb.cb_affinefunctiontransformation_set_fun_offset(self.data::Ptr{Cvoid})::Ptr{Cdouble}

@doc raw"""
    cb_get_linear_cost(self::CBAffineFunctionTransformation)

returns the pointer to the linear term added to the funciton
"""
cb_get_linear_cost(self::CBAffineFunctionTransformation) = CBMatrix(@ccall libcb.cb_affinefunctiontransformation_get_linear_cost(self.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_get_arg_offset(self::CBAffineFunctionTransformation)

returns the pointer to the constant offset added to the argument (not neeeded in the code)
"""
cb_get_arg_offset(self::CBAffineFunctionTransformation) = CBMatrix(@ccall libcb.cb_affinefunctiontransformation_get_arg_offset(self.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_get_arg_trafo(self::CBAffineFunctionTransformation)

returns the pointer to the linear transformation of the argument (not neeeded in the code)
"""
cb_get_arg_trafo(self::CBAffineFunctionTransformation) = CBSparsemat(@ccall libcb.cb_affinefunctiontransformation_get_arg_trafo(self.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_get_linear_cost(self::CBAffineFunctionTransformation, i::Integer)

returns the value of the linear cost coefficient for @a i>=0 and for i==-1 the constant offset
"""
cb_get_linear_cost(self::CBAffineFunctionTransformation, i::Integer) = @ccall libcb.cb_affinefunctiontransformation_get_linear_cost2(self.data::Ptr{Cvoid}, i::Cint)::Cdouble

@doc raw"""
    cb_get_constant_minorant(self::CBAffineFunctionTransformation)

return the constant linear minorant corresponding to fun_offset+<*linear_cost,.>
"""
cb_get_constant_minorant(self::CBAffineFunctionTransformation) = (@ccall libcb.cb_affinefunctiontransformation_get_constant_minorant(self.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

@doc raw"""
    cb_from_dim(self::CBAffineFunctionTransformation)

returns the dimension of the input argument or -1 if it is unknown
"""
cb_from_dim(self::CBAffineFunctionTransformation) = @ccall libcb.cb_affinefunctiontransformation_from_dim(self.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_to_dim(self::CBAffineFunctionTransformation)

returns the dimension of the output argument or -1 if it is unknown
"""
cb_to_dim(self::CBAffineFunctionTransformation) = @ccall libcb.cb_affinefunctiontransformation_to_dim(self.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_copy_traforows(self::CBAffineFunctionTransformation, copy_to::CBMatrix, copy_from::CBMatrix)

only copies the elements that are effected by the image of arg_trafo
"""
cb_copy_traforows(self::CBAffineFunctionTransformation, copy_to::CBMatrix, copy_from::CBMatrix) = @ccall libcb.cb_affinefunctiontransformation_copy_traforows(self.data::Ptr{Cvoid}, copy_to.data::Ptr{Cvoid}, copy_from.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_transform_argument(self::CBAffineFunctionTransformation, transformed_y::CBMatrix, input_y::CBMatrix)

* @brief computes @a transformed_offset and, if @a this is not the identity,@a transformed_y and returns either @a input_y (id) or @a transformed_y

  On input @a transformed_y is supposed to be of dimension zero, unless
  it was already the result of a previous transformation for exactly
  this AFT with this data (no intermediate modifications). This is
  important, because if @a transformed_y has the correct size, only the
  data changed by the transformation is overwritten and the constant
  parts are assumed to be already initialized. On output it is again of
  dimension 0 if the transformation is the identity.

        The @a transformed_offset = fun_offset+ip(linear_cost,@a input_y) is
        the value needed in objective_value() together with the result of the
        evaluation of the function for @a transformed_y in order to compute
        the transformed objective value.
    
"""
function cb_transform_argument(self::CBAffineFunctionTransformation, transformed_y::CBMatrix, input_y::CBMatrix)
    transformed_offset = Ref{Float64}()
    (@ccall libcb.cb_affinefunctiontransformation_transform_argument(self.data::Ptr{Cvoid}, transformed_y.data::Ptr{Cvoid}, transformed_offset::Ref{Float64}, input_y.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)
    return transformed_offset[]
end

@doc raw"""
    cb_modified_transform_argument(self::CBAffineFunctionTransformation, transformed_y::CBMatrix, input_y::CBMatrix, aftmdf::Union{<:CBAFTModification,Nothing}, gsmdf::CBGroundsetModification)

* @brief given the modification aftmdf or if 0, gsmdf, compute the
        transformed argument that would arise after this modification as in
  transform_argument()

  On input @a transformed_y is supposed to be of dimension zero.
  The matrix returned is transformed_y unless the modification
  preserves the identity transformation. In this case, the
        returned matrix is input_y.
    
"""
cb_modified_transform_argument(self::CBAffineFunctionTransformation, transformed_y::CBMatrix, input_y::CBMatrix, aftmdf::Union{<:CBAFTModification,Nothing}, gsmdf::CBGroundsetModification) = (@ccall libcb.cb_affinefunctiontransformation_modified_transform_argument(self.data::Ptr{Cvoid}, transformed_y.data::Ptr{Cvoid}, input_y.data::Ptr{Cvoid}, (isnothing(aftmdf) ? C_NULL : aftmdf.data)::Ptr{Cvoid}, gsmdf.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

@doc raw"""
    cb_objective_value(self::CBAffineFunctionTransformation, offset::Real, function_value::Real)

* @brief  if @a offset is the value computed in transform_argument and
       @a function_value results form an evaluation in the respective point,
       the routine returns the objective value obtained by the affine transformation
     
"""
cb_objective_value(self::CBAffineFunctionTransformation, offset::Real, function_value::Real) = @ccall libcb.cb_affinefunctiontransformation_objective_value(self.data::Ptr{Cvoid}, offset::Cdouble, function_value::Cdouble)::Cdouble

@doc raw"""
    cb_transform_minorant(self::CBAffineFunctionTransformation, out_minorant::CBMinorantPointer, in_minorant::CBMinorantPointer, alpha::Real = 1., add_trafo_minorant::Bool = false, provided_row_indices::Union{<:CBIndexmatrix,Nothing} = nothing, needed_column_indices::Union{<:CBIndexmatrix,Nothing} = nothing)

* @brief transforms the in linear minorant (scaled by alpha) and
        initializes or adds it to the out linear minorant

        if out_minorant.empty()==true, the out_minorant is initialized,
        otherwise the information is added.

        the constant_minorant of the transformation is only added if
        add_trafo_minorant is set to true.

  If given, provided_row_indices and needed_column_indices must be as
        specified in qp_cost_indices, i.e. of in_minorant only the
        provided_row_indices (all if NULL) will be used to compute only the
        needed_column_indices (as rows, all if NULL) of out_minorant.
     
"""
cb_transform_minorant(self::CBAffineFunctionTransformation, out_minorant::CBMinorantPointer, in_minorant::CBMinorantPointer, alpha::Real = 1., add_trafo_minorant::Bool = false, provided_row_indices::Union{<:CBIndexmatrix,Nothing} = nothing, needed_column_indices::Union{<:CBIndexmatrix,Nothing} = nothing) = @ccall libcb.cb_affinefunctiontransformation_transform_minorant(self.data::Ptr{Cvoid}, out_minorant.data::Ptr{Cvoid}, in_minorant.data::Ptr{Cvoid}, alpha::Cdouble, add_trafo_minorant::Cint, (isnothing(provided_row_indices) ? C_NULL : provided_row_indices.data)::Ptr{Cvoid}, (isnothing(needed_column_indices) ? C_NULL : needed_column_indices.data)::Ptr{Cvoid})::Cint

@doc raw"""
    cb_transform_minorants(self::CBAffineFunctionTransformation, out_minorants::CBMinorantBundle, in_minorants::CBMinorantBundle, alpha::Real = 1.)

* @brief transforms several linear minorants (scaled by alpha) and
        initializes or adds them to the out linear minorants

        for out_minorant[i].empty()==true the out_minorant is initialized,
        otherwise the information is added.
     
"""
cb_transform_minorants(self::CBAffineFunctionTransformation, out_minorants::CBMinorantBundle, in_minorants::CBMinorantBundle, alpha::Real = 1.) = @ccall libcb.cb_affinefunctiontransformation_transform_minorants(self.data::Ptr{Cvoid}, out_minorants.data::Ptr{Cvoid}, in_minorants.data::Ptr{Cvoid}, alpha::Cdouble)::Cint

@doc raw"""
    cb_qp_cost_indices(self::CBAffineFunctionTransformation, provide_row_indices::CBIndexmatrix, needed_col_indices::Union{<:CBIndexmatrix,Nothing})

* @brief if the algorithm only requires the indices given in @a
       needed_col_indices (NULL means all indices) then it suffices to supply
       the @a provide_row_indices as @a indices in transform_minorant (a 0x0
       return matrix again means all indices, while a 0x1 matrix means no
       indices). Any input or output indices must be in strictly increasing
       order.
     
"""
cb_qp_cost_indices(self::CBAffineFunctionTransformation, provide_row_indices::CBIndexmatrix, needed_col_indices::Union{<:CBIndexmatrix,Nothing}) = @ccall libcb.cb_affinefunctiontransformation_qp_cost_indices(self.data::Ptr{Cvoid}, provide_row_indices.data::Ptr{Cvoid}, (isnothing(needed_col_indices) ? C_NULL : needed_col_indices.data)::Ptr{Cvoid})::Cint

@doc raw"""
    cb_scaling_indices(self::CBAffineFunctionTransformation, row_indices::CBIndexmatrix, col_indices::CBIndexmatrix)

* @brief if AFTModel::add_diagonal_scaling() is called with @a indices
  specified by @a col_indices, then AFTmodel has to provide a
  diagonal scaling matrix of dimension to_dim() as @a in_diagscale
  in add_diagonal_scaling() with the entries in the output vector
  @a row_indices computed correctly.
     
"""
cb_scaling_indices(self::CBAffineFunctionTransformation, row_indices::CBIndexmatrix, col_indices::CBIndexmatrix) = @ccall libcb.cb_affinefunctiontransformation_scaling_indices(self.data::Ptr{Cvoid}, row_indices.data::Ptr{Cvoid}, col_indices.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_add_diagonal_scaling(self::CBAffineFunctionTransformation, diagscale::CBMatrix, indices::Union{<:CBIndexmatrix,Nothing}, alpha::Real, in_diagscale::CBMatrix)

* @brief transform the scaling matrix in_diagscale of the untransformed
  function model and add it as described in
  BundleMethod::add_diagonal_scaling.

        If indices are given, they must be sorted in striclty increasing order.
    
"""
cb_add_diagonal_scaling(self::CBAffineFunctionTransformation, diagscale::CBMatrix, indices::Union{<:CBIndexmatrix,Nothing}, alpha::Real, in_diagscale::CBMatrix) = @ccall libcb.cb_affinefunctiontransformation_add_diagonal_scaling(self.data::Ptr{Cvoid}, diagscale.data::Ptr{Cvoid}, (isnothing(indices) ? C_NULL : indices.data)::Ptr{Cvoid}, alpha::Cdouble, in_diagscale.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_analyze_modification(self::CBAffineFunctionTransformation, minorant_trafo_differs::Bool, aftmdf::Union{<:CBAFTModification,Nothing}, gsmdf::CBGroundsetModification)

* @brief returns information on the changes in the ground set of the
        transformed arguments and checks whether after the modification the
  nonzero image of the transformed minorants will not differ from before;
        returns NULL on errors.

  If arg_trafo==NULL (act as identity) and if aftmdf!=NULL adds
  any explicit matrix parts or has modifications not consistent with
  maintaining the identity, arg_trafo is virtually first set to an explicit
  identity before executing the modification. If arg_trafo==NULL and
  aftmdf==NULL, the transformations of gsmdf are applied to linear_cost
  and arg_offset.

  If the input @a gsmdf and the AFT modifications in @a aftmdf result in
  a different GroundsetModification of the transformed space, this is
  computed into the member @a local_gsmdf and the returned pointer
  points to this @a local_gsmdf. If all transformations are passed on
  exactly as in @a gsmdf (because the trafo still acts as the identity),
  the returned pointer points to @a gsmdf.

        If the changes affect the transformation at all in its effect
  on the nonzero image of transformed minorants,
        @a minorant_trafo_differs will be set to true otherwise to false.

     
"""
cb_analyze_modification(self::CBAffineFunctionTransformation, minorant_trafo_differs::Bool, aftmdf::Union{<:CBAFTModification,Nothing}, gsmdf::CBGroundsetModification) = CBGroundsetModification(@ccall libcb.cb_affinefunctiontransformation_analyze_modification(self.data::Ptr{Cvoid}, minorant_trafo_differs::Ref{Cint}, (isnothing(aftmdf) ? C_NULL : aftmdf.data)::Ptr{Cvoid}, gsmdf.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_apply_modification!(self::CBAffineFunctionTransformation, aftmdf::Union{<:CBAFTModification,Nothing}, gsmdf::CBGroundsetModification)

* @brief if arg_trafo==NULL (act as identity) and if aftmdf!=NULL adds
  any explicit matrix parts or has modifications not consistent with
  maintaining the identity, arg_trafo is first set to an explicit
  identity before executing the modification. If arg_trafo==NULL and
  aftmdf==NULL, the transformations of gsmdf are applied to linear_cost
  and arg_offset.
     
"""
cb_apply_modification!(self::CBAffineFunctionTransformation, aftmdf::Union{<:CBAFTModification,Nothing}, gsmdf::CBGroundsetModification) = @ccall libcb.cb_affinefunctiontransformation_apply_modification(self.data::Ptr{Cvoid}, (isnothing(aftmdf) ? C_NULL : aftmdf.data)::Ptr{Cvoid}, gsmdf.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_output_aft_data(self::CBAffineFunctionTransformation)

for testing purposes this outputs the data in mfile readable form
"""
cb_output_aft_data(self::CBAffineFunctionTransformation) = @ccall libcb.cb_affinefunctiontransformation_output_aft_data(self.data::Ptr{Cvoid})::Cvoid

