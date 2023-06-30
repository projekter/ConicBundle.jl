@doc raw"""
    CBAFTModification(var_olddim::Integer = 0, row_olddim::Integer = 0, ignore_groundset_modification::Bool = false)

initialize and reset to an unmodified object currently having var_olddim columns and row_olddim rows
"""
CBAFTModification(var_olddim::Integer = 0, row_olddim::Integer = 0, ignore_groundset_modification::Bool = false) = CBAFTModification(@ccall libcb.cb_aftmodification_new(var_olddim::Cint, row_olddim::Cint, ignore_groundset_modification::Cint)::Ptr{Cvoid})

@doc raw"""
    cb_clear!(self::CBAFTModification, var_olddim::Integer, row_olddim::Integer)

reset modifications to an unmodified object currently having var_olddim columns and row_olddim rows,  calls Modification::clear
"""
cb_clear!(self::CBAFTModification, var_olddim::Integer, row_olddim::Integer) = @ccall libcb.cb_aftmodification_clear(self.data::Ptr{Cvoid}, var_olddim::Cint, row_olddim::Cint)::Cvoid

@doc raw"""
    cb_add_append_vars!(self::CBAFTModification, append_dim::Integer, append_cols::Union{<:CBSparsemat,Nothing}, linear_costs::Union{<:CBMatrix,Nothing})

append @a append_dim new variables/columns with values specified by @a append_cols (NULL: default value) and @a linear_costs coefficients (NULL: default value) by calling Modification::add_append_vars
"""
cb_add_append_vars!(self::CBAFTModification, append_dim::Integer, append_cols::Union{<:CBSparsemat,Nothing}, linear_costs::Union{<:CBMatrix,Nothing}) = @ccall libcb.cb_aftmodification_add_append_vars(self.data::Ptr{Cvoid}, append_dim::Cint, (isnothing(append_cols) ? C_NULL : append_cols.data)::Ptr{Cvoid}, (isnothing(linear_costs) ? C_NULL : linear_costs.data)::Ptr{Cvoid})::Cint

@doc raw"""
    cb_add_reassign_vars!(self::CBAFTModification, map_to_old::CBIndexmatrix)

reassign the variables/columns as given in @a map_to_old, calls Modification::add_reassign_vars
"""
cb_add_reassign_vars!(self::CBAFTModification, map_to_old::CBIndexmatrix) = @ccall libcb.cb_aftmodification_add_reassign_vars(self.data::Ptr{Cvoid}, map_to_old.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_add_delete_vars!(self::CBAFTModification, del_ind::CBIndexmatrix, map_to_old::CBIndexmatrix)

delete the variables indexed by @a del_ind; for each new index @a map_to_old returns the old one; calls Modification::add_delete_vars
"""
cb_add_delete_vars!(self::CBAFTModification, del_ind::CBIndexmatrix, map_to_old::CBIndexmatrix) = @ccall libcb.cb_aftmodification_add_delete_vars(self.data::Ptr{Cvoid}, del_ind.data::Ptr{Cvoid}, map_to_old.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_add_append_rows!(self::CBAFTModification, append_dim::Integer, append_rows::Union{<:CBSparsemat,Nothing}, append_rhs::Union{<:CBMatrix,Nothing})

append @a append_dim new rows as in append_rows (if NULL, use default value) with affine offset append_rhs (if NULL use default value 0.), calls Modification::add_append_rows
"""
cb_add_append_rows!(self::CBAFTModification, append_dim::Integer, append_rows::Union{<:CBSparsemat,Nothing}, append_rhs::Union{<:CBMatrix,Nothing}) = @ccall libcb.cb_aftmodification_add_append_rows(self.data::Ptr{Cvoid}, append_dim::Cint, (isnothing(append_rows) ? C_NULL : append_rows.data)::Ptr{Cvoid}, (isnothing(append_rhs) ? C_NULL : append_rhs.data)::Ptr{Cvoid})::Cint

@doc raw"""
    cb_add_reassign_rows!(self::CBAFTModification, map_to_old::CBIndexmatrix)

reassign the rows as given in @a map_to_old, calls Modification::add_reassign_vars
"""
cb_add_reassign_rows!(self::CBAFTModification, map_to_old::CBIndexmatrix) = @ccall libcb.cb_aftmodification_add_reassign_rows(self.data::Ptr{Cvoid}, map_to_old.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_add_delete_rows!(self::CBAFTModification, rows_del_ind::CBIndexmatrix, rows_map_to_old::CBIndexmatrix)

delete the rows indexed by @a rows_del_ind; for each new index @a rows_map_to_old returns the old one; calls Modification::add_delete_rows
"""
cb_add_delete_rows!(self::CBAFTModification, rows_del_ind::CBIndexmatrix, rows_map_to_old::CBIndexmatrix) = @ccall libcb.cb_aftmodification_add_delete_rows(self.data::Ptr{Cvoid}, rows_del_ind.data::Ptr{Cvoid}, rows_map_to_old.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_add_apply_factor!(self::CBAFTModification, times_factor::Real)

multiply the current factor by the value @a times_factor (>=0, returns 1 if <0.)
"""
cb_add_apply_factor!(self::CBAFTModification, times_factor::Real) = @ccall libcb.cb_aftmodification_add_apply_factor(self.data::Ptr{Cvoid}, times_factor::Cdouble)::Cint

@doc raw"""
    cb_add_offset!(self::CBAFTModification, delta::Real)

add to the currecnt offset the value @a delta
"""
cb_add_offset!(self::CBAFTModification, delta::Real) = @ccall libcb.cb_aftmodification_add_offset(self.data::Ptr{Cvoid}, delta::Cdouble)::Cint

@doc raw"""
    cb_incorporate!(self::CBAFTModification, m::CBAFTModification)

incorporate the AFTModification @a m into this one; after factor and offset are dealt with it calls Modification::incorporate
"""
cb_incorporate!(self::CBAFTModification, m::CBAFTModification) = @ccall libcb.cb_aftmodification_incorporate(self.data::Ptr{Cvoid}, m.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_apply_to_factor(self::CBAFTModification)

multiply @a f by the factor
"""
function cb_apply_to_factor(self::CBAFTModification)
    f = Ref{Float64}()
    @ccall libcb.cb_aftmodification_apply_to_factor(self.data::Ptr{Cvoid}, f::Ref{Float64})::Cint
    return f[]
end

@doc raw"""
    cb_apply_to_offset(self::CBAFTModification)

add to @a o the offset
"""
function cb_apply_to_offset(self::CBAFTModification)
    o = Ref{Float64}()
    @ccall libcb.cb_aftmodification_apply_to_offset(self.data::Ptr{Cvoid}, o::Ref{Float64})::Cint
    return o[]
end

@doc raw"""
    cb_apply_modified_transform(self::CBAFTModification, out_y::CBMatrix, in_y::CBMatrix, arg_trafo::Union{<:CBSparsemat,Nothing}, arg_offset::Union{<:CBMatrix,Nothing})

transform the vector as if the Modification had been carried out
"""
cb_apply_modified_transform(self::CBAFTModification, out_y::CBMatrix, in_y::CBMatrix, arg_trafo::Union{<:CBSparsemat,Nothing}, arg_offset::Union{<:CBMatrix,Nothing}) = (@ccall libcb.cb_aftmodification_apply_modified_transform(self.data::Ptr{Cvoid}, out_y.data::Ptr{Cvoid}, in_y.data::Ptr{Cvoid}, (isnothing(arg_trafo) ? C_NULL : arg_trafo.data)::Ptr{Cvoid}, (isnothing(arg_offset) ? C_NULL : arg_offset.data)::Ptr{Cvoid})::Ptr{Cvoid}; return self)

@doc raw"""
    cb_no_modification(self::CBAFTModification)

returns true if no modifications need to be executed
"""
cb_no_modification(self::CBAFTModification) = Bool(@ccall libcb.cb_aftmodification_no_modification(self.data::Ptr{Cvoid})::Cint)

@doc raw"""
    cb_set_append_to_old!(self::CBAFTModification, append_only::Bool)

if set to true, no deletions/reassignments may be present or specified in the future, only appensions are allowed
"""
cb_set_append_to_old!(self::CBAFTModification, append_only::Bool) = @ccall libcb.cb_aftmodification_set_append_to_old(self.data::Ptr{Cvoid}, append_only::Cint)::Cint

@doc raw"""
    cb_append_to_old(self::CBAFTModification)

returns true if this only contains appending operations and incorporating this is done with respect to the old dimension
"""
cb_append_to_old(self::CBAFTModification) = Bool(@ccall libcb.cb_aftmodification_append_to_old(self.data::Ptr{Cvoid})::Cint)

@doc raw"""
    cb_only_scalars_change(self::CBAFTModification)

returns true if linear_cost, matrix and affine rhs offset are not changed
"""
cb_only_scalars_change(self::CBAFTModification) = Bool(@ccall libcb.cb_aftmodification_only_scalars_change(self.data::Ptr{Cvoid})::Cint)

@doc raw"""
    cb_ignore_groundset_modification(self::CBAFTModification)

returns true if groundset_modifications should be ignored
"""
cb_ignore_groundset_modification(self::CBAFTModification) = Bool(@ccall libcb.cb_aftmodification_ignore_groundset_modification(self.data::Ptr{Cvoid})::Cint)

@doc raw"""
    cb_groundset_changes_suffice_for_identity!(self::CBAFTModification)

returns true if for an AFT with argtrafo==0 the changes in the ground set reflect all modifications
"""
cb_groundset_changes_suffice_for_identity!(self::CBAFTModification) = Bool(@ccall libcb.cb_aftmodification_groundset_changes_suffice_for_identity(self.data::Ptr{Cvoid})::Cint)

@doc raw"""
    cb_preserves_identity(self::CBAFTModification)

returns true if the modifications are consistent with the AffineFunctionTransformation matrix staying the identity
"""
cb_preserves_identity(self::CBAFTModification) = Bool(@ccall libcb.cb_aftmodification_preserves_identity(self.data::Ptr{Cvoid})::Cint)

@doc raw"""
    cb_no_additions_or_deletions_in_vars(self::CBAFTModification)

returns true if no columns/variables were added or deleted (allows permutations), false otherwise
"""
cb_no_additions_or_deletions_in_vars(self::CBAFTModification) = Bool(@ccall libcb.cb_aftmodification_no_additions_or_deletions_in_vars(self.data::Ptr{Cvoid})::Cint)

@doc raw"""
    cb_no_additions_or_deletions_in_rows(self::CBAFTModification)

returns true if no rows were added or deleted (allows permutations), false otherwise
"""
cb_no_additions_or_deletions_in_rows(self::CBAFTModification) = Bool(@ccall libcb.cb_aftmodification_no_additions_or_deletions_in_rows(self.data::Ptr{Cvoid})::Cint)

@doc raw"""
    cb_deleted_variables_are_zero(self::CBAFTModification, oldpoint::CBMatrix)

returns true if all entries deleted in @a oldpoint (must be a vector of length old_vardim()) are 0,  false otherwise
"""
cb_deleted_variables_are_zero(self::CBAFTModification, oldpoint::CBMatrix) = Bool(@ccall libcb.cb_aftmodification_deleted_variables_are_zero(self.data::Ptr{Cvoid}, oldpoint.data::Ptr{Cvoid})::Cint)

@doc raw"""
    cb_new_variables_are_zero(self::CBAFTModification, newpoint::CBMatrix)

returns true if all entries in newpoint (must be a vector of length new_vardim()) that correspond to new variables have value 0 and false otherwise
"""
cb_new_variables_are_zero(self::CBAFTModification, newpoint::CBMatrix) = Bool(@ccall libcb.cb_aftmodification_new_variables_are_zero(self.data::Ptr{Cvoid}, newpoint.data::Ptr{Cvoid})::Cint)

@doc raw"""
    cb_mapped_variables_are_equal(self::CBAFTModification, newpoint::CBMatrix, oldpoint::CBMatrix)

returns true if the values in newpoint (must be a vector of length new_vardim()) that correspond to old variables match the old values stored in oldpoint (must be a vector of length old_vardim()) and false otherwise
"""
cb_mapped_variables_are_equal(self::CBAFTModification, newpoint::CBMatrix, oldpoint::CBMatrix) = Bool(@ccall libcb.cb_aftmodification_mapped_variables_are_equal(self.data::Ptr{Cvoid}, newpoint.data::Ptr{Cvoid}, oldpoint.data::Ptr{Cvoid})::Cint)

@doc raw"""
    cb_get_old_vardim(self::CBAFTModification)

returns the number of variables before modification
"""
cb_get_old_vardim(self::CBAFTModification) = @ccall libcb.cb_aftmodification_get_old_vardim(self.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_get_new_vardim(self::CBAFTModification)

returns the number of variables once all stored modifications have been performed
"""
cb_get_new_vardim(self::CBAFTModification) = @ccall libcb.cb_aftmodification_get_new_vardim(self.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_get_appended_vardim(self::CBAFTModification)

returns the number of variables that are appended (due to later reassignmentds they may no longer be located at the end)
"""
cb_get_appended_vardim(self::CBAFTModification) = @ccall libcb.cb_aftmodification_get_appended_vardim(self.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_get_map_to_old_variables(self::CBAFTModification)

returns null if there are no index changes, otherwise the Indexmatrix pointed to is a vector whose i-th entry holds the old index of the new i-th variable (injective!), index values exceeding old_vardim() refer to newly appended variables
"""
cb_get_map_to_old_variables(self::CBAFTModification) = @ccall libcb.cb_aftmodification_get_map_to_old_variables(self.data::Ptr{Cvoid})::Ptr{Cint}

@doc raw"""
    cb_incorporate!(self::CBAFTModification, m::CBOracleModification)

incorporate the OracleModification @a m (it should only contain variable changes, but this is not checked!) into this one; calls Modification::incorporate
"""
cb_incorporate!(self::CBAFTModification, m::CBOracleModification) = @ccall libcb.cb_aftmodification_incorporate2(self.data::Ptr{Cvoid}, m.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_new_initial_oraclemodification(self::CBAFTModification, old_var_dim::Integer)

returns a new object on the heap, that allows to incorporate this but starts off from a function whose input argument dimension is old_var_dim
"""
cb_new_initial_oraclemodification(self::CBAFTModification, old_var_dim::Integer) = CBOracleModification(@ccall libcb.cb_aftmodification_new_initial_oraclemodification(self.data::Ptr{Cvoid}, old_var_dim::Cint)::Ptr{Cvoid})

@doc raw"""
    cb_add_append_variables!(self::CBAFTModification, append_dim::Integer)

append append_dim further variables at the end of the argument vector (without specifying their effect, they may e.g. be ignored by the function)
"""
cb_add_append_variables!(self::CBAFTModification, append_dim::Integer) = @ccall libcb.cb_aftmodification_add_append_variables(self.data::Ptr{Cvoid}, append_dim::Cint)::Cint

@doc raw"""
    cb_add_reassign_variables!(self::CBAFTModification, new_dim::Integer, map_to_old_indices::Union{<:AbstractVector{Integer},Nothing})

reorder and resize the variables as given by the first new_dim entries of map_to_old_indices; each former index may appear at most once in this list, so new_dim < get_new_vardim()
"""
cb_add_reassign_variables!(self::CBAFTModification, new_dim::Integer, map_to_old_indices::Union{<:AbstractVector{Integer},Nothing}) = GC.@preserve map_to_old_indices begin
    (LinearAlgebra.chkstride1(map_to_old_indices); @ccall libcb.cb_aftmodification_add_reassign_variables(self.data::Ptr{Cvoid}, new_dim::Cint, map_to_old_indices::Ptr{Cint})::Cint)
end

@doc raw"""
    cb_old_vardim(self::CBAFTModification)

returns the number of variables before modification (given on initialization)
"""
cb_old_vardim(self::CBAFTModification) = @ccall libcb.cb_aftmodification_old_vardim(self.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_new_vardim(self::CBAFTModification)

returns the number of variables once all stored modifications have been performed
"""
cb_new_vardim(self::CBAFTModification) = @ccall libcb.cb_aftmodification_new_vardim(self.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_appended_vardim(self::CBAFTModification)

returns the number of variables that are appended (due to later reassignmentds they may no longer be located at the end)
"""
cb_appended_vardim(self::CBAFTModification) = @ccall libcb.cb_aftmodification_appended_vardim(self.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_old_rowdim(self::CBAFTModification)

returns the number of rows before modification (given on initialization)
"""
cb_old_rowdim(self::CBAFTModification) = @ccall libcb.cb_aftmodification_old_rowdim(self.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_new_rowdim(self::CBAFTModification)

returns the number of rows once all stored modifications have been performed
"""
cb_new_rowdim(self::CBAFTModification) = @ccall libcb.cb_aftmodification_new_rowdim(self.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_appended_rowdim(self::CBAFTModification)

returns the number of rows that are appended (due to later reassignments they may no longer be located at the end)
"""
cb_appended_rowdim(self::CBAFTModification) = @ccall libcb.cb_aftmodification_appended_rowdim(self.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_map_to_old_variables(self::CBAFTModification)

returns null if there are no index changes, otherwise the Indexmatrix pointed to is a vector whose i-th entry holds the old index of the new i-th variable (injective!), index values exceeding old_vardim() refer to newly appended variables
"""
cb_map_to_old_variables(self::CBAFTModification) = CBIndexmatrix(@ccall libcb.cb_aftmodification_map_to_old_variables(self.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_deleted_var_indices(self::CBAFTModification)

returns null if there were no deletions, otherwise the Indexmatrix pointed to is a vector holding the deleted old variable indices in increasing order
"""
cb_deleted_var_indices(self::CBAFTModification) = CBIndexmatrix(@ccall libcb.cb_aftmodification_deleted_var_indices(self.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_new_var_indices(self::CBAFTModification)

returns null if no variables were added, otherwise the Indexmatrix pointed to is a vector holding the new indices of the new variables in increasing order
"""
cb_new_var_indices(self::CBAFTModification) = CBIndexmatrix(@ccall libcb.cb_aftmodification_new_var_indices(self.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_map_to_old_rows(self::CBAFTModification)

returns null if there are no index changes, otherwise the Indexmatrix pointed to is a vector whose i-th entry holds the old index of the new i-th row (injective!), index values exceeding old_rowdim() refer to newly appended rows
"""
cb_map_to_old_rows(self::CBAFTModification) = CBIndexmatrix(@ccall libcb.cb_aftmodification_map_to_old_rows(self.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_deleted_row_indices(self::CBAFTModification)

returns null if there were no deletions, otherwise the Indexmatrix pointed to is a vector holding the deleted old row indices in increasing order
"""
cb_deleted_row_indices(self::CBAFTModification) = CBIndexmatrix(@ccall libcb.cb_aftmodification_deleted_row_indices(self.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_new_row_indices(self::CBAFTModification)

returns null if no variables were added, otherwise the Indexmatrix pointed to is a vector holding the new indices of the new rows in increasing order
"""
cb_new_row_indices(self::CBAFTModification) = CBIndexmatrix(@ccall libcb.cb_aftmodification_new_row_indices(self.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_get_additional_offset(self::CBAFTModification)

returns the value to be added to the offset
"""
cb_get_additional_offset(self::CBAFTModification) = @ccall libcb.cb_aftmodification_get_additional_offset(self.data::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_get_additional_factor(self::CBAFTModification)

returns null if nothing or default values have to be appended, otherwise it points to a sparse matrix whose columns need to be appended to the matrix
"""
cb_get_additional_factor(self::CBAFTModification) = @ccall libcb.cb_aftmodification_get_additional_factor(self.data::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_get_append_cols(self::CBAFTModification)

returns null if nothing or default values have to be appended, otherwise it points to a sparse matrix whose columns need to be appended to the matrix
"""
cb_get_append_cols(self::CBAFTModification) = CBSparsemat(@ccall libcb.cb_aftmodification_get_append_cols(self.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_get_append_costs(self::CBAFTModification)

returns null if nothing or default values have to be appended, otherwise it points to a matrix whose entries need to be appended to the cost vector
"""
cb_get_append_costs(self::CBAFTModification) = CBMatrix(@ccall libcb.cb_aftmodification_get_append_costs(self.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_get_append_rows(self::CBAFTModification)

returns null if nothing or default values have to be appended, otherwise it points to a sparse matrix whose rows need to be appended to the matrix
"""
cb_get_append_rows(self::CBAFTModification) = CBSparsemat(@ccall libcb.cb_aftmodification_get_append_rows(self.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_get_append_rhs(self::CBAFTModification)

returns null if nothing or default values have to be appended, otherwise it points to a matrix whose rows need to be appended to the argument offset
"""
cb_get_append_rhs(self::CBAFTModification) = CBMatrix(@ccall libcb.cb_aftmodification_get_append_rhs(self.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_set_out!(self::CBAFTModification, print_level::Integer = 1)

see CBout::set_out
"""
cb_set_out!(self::CBAFTModification, print_level::Integer = 1) = @ccall libcb.cb_aftmodification_set_out(self.data::Ptr{Cvoid}, print_level::Cint)::Cvoid

