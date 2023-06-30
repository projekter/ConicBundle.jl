@doc raw"""
    CBPSCAffineModification(var_olddim::Integer, block_olddim::CBIndexmatrix, incr::Integer = 0)

calls clear() with these parameters
"""
CBPSCAffineModification(var_olddim::Integer, block_olddim::CBIndexmatrix, incr::Integer = 0) = CBPSCAffineModification(@ccall libcb.cb_pscaffinemodification_new(var_olddim::Cint, block_olddim.data::Ptr{Cvoid}, incr::Cint)::Ptr{Cvoid})

@doc raw"""
    cb_clear!(self::CBPSCAffineModification, var_olddim::Integer, block_olddim::CBIndexmatrix)

* @brief resets all variables so that the object to be modified has
        starting size var_olddim (number of variables) and block_olddim
        (number of rows) and no modifications

        The actual old data is not needed at this point,
        the changes on it will be collected and excuted in the
        routines apply_to_vars and apply_to_blocks

        Setting the parameter ensure_start_val_box_feasibility to true
        will cause the algorithm to check in add_append_vars() whether
        the input values are within the given bounds and in
        apply_to_vars() it will project all start_values onto the bounds
        for all, old and new, indices (which might have been changed by
        then).  If it is false (default), all values will be accepted as
        given.

        Setting the parameter ensure_bounds_consistency to true
        (default) will raise errors in add_append_vars() and in
        apply_to_vars() whenever there are lower bounds greater than the
        respective upper bounds so as to avoid trivial
        infeasibilities. This check is omitted if set to false.

        The remaining values give the values of plus and minus infinity
        the no bounds should exceed. These are the default values at the
        same time (maybe it might be good to have a separate default value,
        but this is not implemented here).
    
"""
cb_clear!(self::CBPSCAffineModification, var_olddim::Integer, block_olddim::CBIndexmatrix) = @ccall libcb.cb_pscaffinemodification_clear(self.data::Ptr{Cvoid}, var_olddim::Cint, block_olddim.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_add_append_vars!(self::CBPSCAffineModification, append_dim::Integer, append_cols::Union{<:CBSparseCoeffmatMatrix,Nothing})

* @brief append information on new variables at the respective ends

        @param append_dim
               number of variables (or columns ) to
         be appended

        @param append_cols
               if NULL, append zero matrices, otherwise it must point to
         a sparse matrix of size new_rowdim() times @a append_dim
               that is to be appended to the constraint matrix on the right

        @return number of errors; if errors occured, none of the new changes are performed
     
"""
cb_add_append_vars!(self::CBPSCAffineModification, append_dim::Integer, append_cols::Union{<:CBSparseCoeffmatMatrix,Nothing}) = @ccall libcb.cb_pscaffinemodification_add_append_vars(self.data::Ptr{Cvoid}, append_dim::Cint, (isnothing(append_cols) ? C_NULL : append_cols.data)::Ptr{Cvoid})::Cint

@doc raw"""
    cb_add_reassign_vars!(self::CBPSCAffineModification, map_to_old::CBIndexmatrix)

* @brief reassign the current variable indices (with modifications) as specified by @a map_to_old

        @a map_to_old must specify an injective map (no two values match)
        into indices 0 up to new_vardim()-1 (not all need to appear). The variable getting index
        i (for i=0 to map_to_old.dim()-1) is the variable with current
        index map_to_old(i) (current refers to considering all previous
        modifications as having been carried out already). The return
        value is the number of errors in @a map_to_old. If such occured, this
        reassign is not performed.
     
"""
cb_add_reassign_vars!(self::CBPSCAffineModification, map_to_old::CBIndexmatrix) = @ccall libcb.cb_pscaffinemodification_add_reassign_vars(self.data::Ptr{Cvoid}, map_to_old.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_add_delete_vars!(self::CBPSCAffineModification, del_ind::CBIndexmatrix, map_to_old::CBIndexmatrix)

* @brief delete the variables indexed by the vector del_ind and
        return the index changes of the others in a vector map_to_old

        @a del_ind need not be ordered in any way, but any index in 0 to
        new_vardim()-1 may appear at most once. On output the dimension
        of the column vector @a map_to_old gives the number of remaining
        variables and its entry i holds the index the variable with new
        index i had before the deletion. The return value is the number
        of errors in @a del_ind. If such occured, this deletion is
        not performed and @a map_to_old may contain garbage.
    
"""
cb_add_delete_vars!(self::CBPSCAffineModification, del_ind::CBIndexmatrix, map_to_old::CBIndexmatrix) = @ccall libcb.cb_pscaffinemodification_add_delete_vars(self.data::Ptr{Cvoid}, del_ind.data::Ptr{Cvoid}, map_to_old.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_add_append_blocks!(self::CBPSCAffineModification, append_dim::CBIndexmatrix, append_offsets::Union{<:CBSparseCoeffmatMatrix,Nothing}, append_blocks::Union{<:CBSparseCoeffmatMatrix,Nothing})

* @brief append information on new rows at the respective ends

       @param append_dim
              the dimension gives the number of (diagonal) blocks and each
        entry the order of the block to be appended

       @param append_blocks
              if NULL, append zero matrics, otherwise it must point
              to a sparse matrix of size @a append_dim times new_vardim()
              that is to be appended to the constraint matrix below.

       @param append_offsets
              if NULL, append zero matrices, otherwise it must point to
        a column vector of size @a append_dim that is to be appended
              to the vector of offsets


       @return number of errors; if errors occured, none of the new changes are performed
    
"""
cb_add_append_blocks!(self::CBPSCAffineModification, append_dim::CBIndexmatrix, append_offsets::Union{<:CBSparseCoeffmatMatrix,Nothing}, append_blocks::Union{<:CBSparseCoeffmatMatrix,Nothing}) = @ccall libcb.cb_pscaffinemodification_add_append_blocks(self.data::Ptr{Cvoid}, append_dim.data::Ptr{Cvoid}, (isnothing(append_offsets) ? C_NULL : append_offsets.data)::Ptr{Cvoid}, (isnothing(append_blocks) ? C_NULL : append_blocks.data)::Ptr{Cvoid})::Cint

@doc raw"""
    cb_add_reassign_blocks!(self::CBPSCAffineModification, map_to_old::CBIndexmatrix)

* @brief reassign the current row indices (with modifications) as specified by @a map_to_old

        @a map_to_old must specify an injective map (no two values
        match) into indices 0 up to new_rowdim()-1 (not all need to
        appear).  The row getting index i (for i=0 to
        map_to_old.dim()-1) is the row with current index map_to_old(i)
        (current refers to considering all previous modifications as
        having been carried out already). The return value is the number
        of errors in @a map_to_old. If such occured, this reassign is
        not performed.
     
"""
cb_add_reassign_blocks!(self::CBPSCAffineModification, map_to_old::CBIndexmatrix) = @ccall libcb.cb_pscaffinemodification_add_reassign_blocks(self.data::Ptr{Cvoid}, map_to_old.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_add_delete_blocks!(self::CBPSCAffineModification, del_ind::CBIndexmatrix, map_to_old::CBIndexmatrix)

* @brief delete the rows indexed by the vector del_ind and
        return the index changes of the others in a vector map_to_old

        @a del_ind need not be ordered in any way, but any index in 0 to
        new_rowdim()-1 may appear at most once. On output the dimension
        of the column vector @a map_to_old gives the number of remaining
        rows and its entry i holds the index the row with new index i
        had before the deletion. The return value is the number of
        errors in @a del_ind. If such occured, this deletion is not
        performed and @a map_to_old may contain garbage.
    
"""
cb_add_delete_blocks!(self::CBPSCAffineModification, del_ind::CBIndexmatrix, map_to_old::CBIndexmatrix) = @ccall libcb.cb_pscaffinemodification_add_delete_blocks(self.data::Ptr{Cvoid}, del_ind.data::Ptr{Cvoid}, map_to_old.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_add_reset_generating_primal!(self::CBPSCAffineModification, new_generating_primal::Union{<:CBPSCPrimal,Nothing})

* @brief replace the current generating primal by new_generating_primal (on the heap, will be deleted here). It may be NULL to switch of generating primals.

        Any changes here cause the deletion of all current aggregate minorants.
        If not NULL, the object pointed to is on the heap and control over it
        is passed over to *this, so *this will make sure it is deleted at due
        time.

        PSCPrimal gives no information on the order of the matrices involved, so no consistency checks are done here.
    
"""
cb_add_reset_generating_primal!(self::CBPSCAffineModification, new_generating_primal::Union{<:CBPSCPrimal,Nothing}) = @ccall libcb.cb_pscaffinemodification_add_reset_generating_primal(self.data::Ptr{Cvoid}, (isnothing(new_generating_primal) ? C_NULL : new_generating_primal.data)::Ptr{Cvoid})::Cint

@doc raw"""
    cb_set_skip_extension!(self::CBPSCAffineModification, skip::Bool)

* @brief if this time no extension is possible for newly added variables
        with the availabel generating primal, set this to true;
     
"""
cb_set_skip_extension!(self::CBPSCAffineModification, skip::Bool) = @ccall libcb.cb_pscaffinemodification_set_skip_extension(self.data::Ptr{Cvoid}, skip::Cint)::Cint

@doc raw"""
    cb_apply_to_PSCAffine(self::CBPSCAffineModification, offset::Union{<:CBSparseCoeffmatMatrix,Nothing}, matrix::Union{<:CBSparseCoeffmatMatrix,Nothing})

* @brief carry out the collected modifications on the data describing the PSCAffineFunction

         If a specific parameter is NULL, no changes are performed on it,
         if it is not null, it must point to a SparseCoeffmatMatrix
         whose sizes correspond to the "old" data. Then the following
         operations will be performed on it in this sequence:

         1. new variables (columns) are appended to the matrix

         2. if reassignment information is given, the columns are
            mapped/reorderd as given by *map_to_old_variables()

         3. new blocks (rows) are appended to the offset and the matrix

         4. if reassignment information is given, the rows are
            mappen/reorderd as given by *map_to_old_blocks()

         @param[in,out] offset
            if not NULL, this points to the offset vector.

         @param[in,out] matrix
            if not NULL, this points to the old matrix.

         @return the number of dimension errors of non NULL inputs,
            if any, no modifications are made to any inputs.
     
"""
cb_apply_to_PSCAffine(self::CBPSCAffineModification, offset::Union{<:CBSparseCoeffmatMatrix,Nothing}, matrix::Union{<:CBSparseCoeffmatMatrix,Nothing}) = @ccall libcb.cb_pscaffinemodification_apply_to_pscaffine(self.data::Ptr{Cvoid}, (isnothing(offset) ? C_NULL : offset.data)::Ptr{Cvoid}, (isnothing(matrix) ? C_NULL : matrix.data)::Ptr{Cvoid})::Cint

@doc raw"""
    cb_no_modification(self::CBPSCAffineModification)

returns true if no modifications need to be executed
"""
cb_no_modification(self::CBPSCAffineModification) = Bool(@ccall libcb.cb_pscaffinemodification_no_modification(self.data::Ptr{Cvoid})::Cint)

@doc raw"""
    cb_set_append_to_old!(self::CBPSCAffineModification, append_only::Bool)

if set to true, no deletions/reassignments may be present or specified in the future, only appensions are allowed
"""
cb_set_append_to_old!(self::CBPSCAffineModification, append_only::Bool) = @ccall libcb.cb_pscaffinemodification_set_append_to_old(self.data::Ptr{Cvoid}, append_only::Cint)::Cint

@doc raw"""
    cb_append_to_old(self::CBPSCAffineModification)

returns true if this only contains appending operations and incorporating this is done with respect to the old dimension
"""
cb_append_to_old(self::CBPSCAffineModification) = Bool(@ccall libcb.cb_pscaffinemodification_append_to_old(self.data::Ptr{Cvoid})::Cint)

@doc raw"""
    cb_deleted_variables_are_zero(self::CBPSCAffineModification, oldpoint::CBMatrix, oldmat::CBSparseCoeffmatMatrix)

returns true if all entries deleted in @a oldpoint (must be a vector of length old_vardim()) or the corresponding entries in the old matrix are 0 and false otherwise
"""
cb_deleted_variables_are_zero(self::CBPSCAffineModification, oldpoint::CBMatrix, oldmat::CBSparseCoeffmatMatrix) = Bool(@ccall libcb.cb_pscaffinemodification_deleted_variables_are_zero(self.data::Ptr{Cvoid}, oldpoint.data::Ptr{Cvoid}, oldmat.data::Ptr{Cvoid})::Cint)

@doc raw"""
    cb_new_variables_are_zero(self::CBPSCAffineModification, newpoint::CBMatrix, newmat::CBSparseCoeffmatMatrix)

returns true if for all indices of new variables the entries in newpoint (must be a vector of length new_vardim()) or the matrices in newmat are 0 and false otherwise
"""
cb_new_variables_are_zero(self::CBPSCAffineModification, newpoint::CBMatrix, newmat::CBSparseCoeffmatMatrix) = Bool(@ccall libcb.cb_pscaffinemodification_new_variables_are_zero(self.data::Ptr{Cvoid}, newpoint.data::Ptr{Cvoid}, newmat.data::Ptr{Cvoid})::Cint)

@doc raw"""
    cb_mapped_variables_are_equal(self::CBPSCAffineModification, newpoint::CBMatrix, oldpoint::CBMatrix)

returns true if the values in newpoint (must be a vector of length new_vardim()) that correspond to old variables match the old values stored in oldpoint (must be a vector of length old_vardim()) and false otherwise
"""
cb_mapped_variables_are_equal(self::CBPSCAffineModification, newpoint::CBMatrix, oldpoint::CBMatrix) = Bool(@ccall libcb.cb_pscaffinemodification_mapped_variables_are_equal(self.data::Ptr{Cvoid}, newpoint.data::Ptr{Cvoid}, oldpoint.data::Ptr{Cvoid})::Cint)

@doc raw"""
    cb_variable_modifications(self::CBPSCAffineModification)

returns true if some modifications are performed on the block structure
"""
cb_variable_modifications(self::CBPSCAffineModification) = Bool(@ccall libcb.cb_pscaffinemodification_variable_modifications(self.data::Ptr{Cvoid})::Cint)

@doc raw"""
    cb_block_modifications(self::CBPSCAffineModification)

returns true if some modifications are performed on the block structure
"""
cb_block_modifications(self::CBPSCAffineModification) = Bool(@ccall libcb.cb_pscaffinemodification_block_modifications(self.data::Ptr{Cvoid})::Cint)

@doc raw"""
    cb_old_vardim(self::CBPSCAffineModification)

returns the number of variables before modification (given on initialization)
"""
cb_old_vardim(self::CBPSCAffineModification) = @ccall libcb.cb_pscaffinemodification_old_vardim(self.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_new_vardim(self::CBPSCAffineModification)

returns the number of variables once all stored modifications have been performed
"""
cb_new_vardim(self::CBPSCAffineModification) = @ccall libcb.cb_pscaffinemodification_new_vardim(self.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_appended_vardim(self::CBPSCAffineModification)

returns the number of variables that are appended (due to later reassignmentds they may no longer be located at the end)
"""
cb_appended_vardim(self::CBPSCAffineModification) = @ccall libcb.cb_pscaffinemodification_appended_vardim(self.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_old_blockdim(self::CBPSCAffineModification)

returns the number of rows before modification (given on initialization)
"""
cb_old_blockdim(self::CBPSCAffineModification) = (@ccall libcb.cb_pscaffinemodification_old_blockdim(self.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

@doc raw"""
    cb_new_blockdim(self::CBPSCAffineModification)

returns the number of rows once all stored modifications have been performed
"""
cb_new_blockdim(self::CBPSCAffineModification) = (@ccall libcb.cb_pscaffinemodification_new_blockdim(self.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

@doc raw"""
    cb_appended_blockdim(self::CBPSCAffineModification)

returns the number of rows that are appended (due to later reassignments they may no longer be located at the end)
"""
cb_appended_blockdim(self::CBPSCAffineModification) = (@ccall libcb.cb_pscaffinemodification_appended_blockdim(self.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

@doc raw"""
    cb_map_to_old_variables(self::CBPSCAffineModification)

returns null if there are no index changes, otherwise the Indexmatrix pointed to is a vector whose i-th entry holds the old index of the new i-th variable (injective!), index values exceeding old_vardim() refer to newly appended variables
"""
cb_map_to_old_variables(self::CBPSCAffineModification) = CBIndexmatrix(@ccall libcb.cb_pscaffinemodification_map_to_old_variables(self.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_deleted_var_indices(self::CBPSCAffineModification)

returns null if there were no deletions, otherwise the Indexmatrix pointed to is a vector holding the deleted old variable indices in increasing order
"""
cb_deleted_var_indices(self::CBPSCAffineModification) = CBIndexmatrix(@ccall libcb.cb_pscaffinemodification_deleted_var_indices(self.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_new_var_indices(self::CBPSCAffineModification)

returns null if no variables were added, otherwise the Indexmatrix pointed to is a vector holding the new indices of the new variables in increasing order
"""
cb_new_var_indices(self::CBPSCAffineModification) = CBIndexmatrix(@ccall libcb.cb_pscaffinemodification_new_var_indices(self.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_map_to_old_blocks(self::CBPSCAffineModification)

returns null if there are index changes, otherwise the Indexmatrix pointed to is a vector whose i-th entry holds the old index of the new i-th row (injective!), index values exceeding old_rowdim() refer to newly appended rows
"""
cb_map_to_old_blocks(self::CBPSCAffineModification) = CBIndexmatrix(@ccall libcb.cb_pscaffinemodification_map_to_old_blocks(self.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_deleted_block_indices(self::CBPSCAffineModification)

returns null if there were no deletions, otherwise the Indexmatrix pointed to is a vector holding the deleted old row indices in increasing order
"""
cb_deleted_block_indices(self::CBPSCAffineModification) = CBIndexmatrix(@ccall libcb.cb_pscaffinemodification_deleted_block_indices(self.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_new_block_indices(self::CBPSCAffineModification)

returns null if no rows were added, otherwise the Indexmatrix pointed ato is a vector holding the new indices of the new rows in increasing order
"""
cb_new_block_indices(self::CBPSCAffineModification) = CBIndexmatrix(@ccall libcb.cb_pscaffinemodification_new_block_indices(self.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_get_var_append(self::CBPSCAffineModification)

returns null if nothing or default values have to be appended, otherwise it points to a sparse matrix whose columns need to be appended to the matrix
"""
cb_get_var_append(self::CBPSCAffineModification) = (@ccall libcb.cb_pscaffinemodification_get_var_append(self.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

@doc raw"""
    cb_get_block_append(self::CBPSCAffineModification)

returns null if nothing or default values have to be appended, otherwise it points to a sparse matrix whose rows need to be appended to the matrix
"""
cb_get_block_append(self::CBPSCAffineModification) = (@ccall libcb.cb_pscaffinemodification_get_block_append(self.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

@doc raw"""
    cb_get_offset_append(self::CBPSCAffineModification)

returns null if nothing or default values have to be appended, otherwise it points to a matrix whose entries need to be appended to the right hand side lower bounds vector
"""
cb_get_offset_append(self::CBPSCAffineModification) = (@ccall libcb.cb_pscaffinemodification_get_offset_append(self.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

@doc raw"""
    cb_get_reset_primal(self::CBPSCAffineModification)

returns true if the generating primal is to be replaced by the one stored here
"""
cb_get_reset_primal(self::CBPSCAffineModification) = Bool(@ccall libcb.cb_pscaffinemodification_get_reset_primal(self.data::Ptr{Cvoid})::Cint)

@doc raw"""
    cb_get_generating_primal(self::CBPSCAffineModification)

returns the generating primal pointer stored here (may be NULL); if get_reset_primal() is true, the PSCAffineFunction should either delete its generating primal (NULL) or replace its generating primal by a clone of this one
"""
cb_get_generating_primal(self::CBPSCAffineModification) = CBPSCPrimal(@ccall libcb.cb_pscaffinemodification_get_generating_primal(self.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_get_skip_extension(self::CBPSCAffineModification)

returns true if the generating primal is to be replaced by the one stored here
"""
cb_get_skip_extension(self::CBPSCAffineModification) = Bool(@ccall libcb.cb_pscaffinemodification_get_skip_extension(self.data::Ptr{Cvoid})::Cint)

@doc raw"""
    cb_get_old_vardim(self::CBPSCAffineModification)

returns the number of variables before modification
"""
cb_get_old_vardim(self::CBPSCAffineModification) = @ccall libcb.cb_pscaffinemodification_get_old_vardim(self.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_get_new_vardim(self::CBPSCAffineModification)

returns the number of variables once all stored modifications have been performed
"""
cb_get_new_vardim(self::CBPSCAffineModification) = @ccall libcb.cb_pscaffinemodification_get_new_vardim(self.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_get_appended_vardim(self::CBPSCAffineModification)

returns the number of variables that are appended (due to later reassignmentds they may no longer be located at the end)
"""
cb_get_appended_vardim(self::CBPSCAffineModification) = @ccall libcb.cb_pscaffinemodification_get_appended_vardim(self.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_get_map_to_old_variables(self::CBPSCAffineModification)

returns null if there are no index changes, otherwise the Indexmatrix pointed to is a vector whose i-th entry holds the old index of the new i-th variable (injective!), index values exceeding old_vardim() refer to newly appended variables
"""
cb_get_map_to_old_variables(self::CBPSCAffineModification) = @ccall libcb.cb_pscaffinemodification_get_map_to_old_variables(self.data::Ptr{Cvoid})::Ptr{Cint}

@doc raw"""
    cb_incorporate!(self::CBPSCAffineModification, m::CBOracleModification)

* @brief add the modification specified in @a m on top of
        the modifications collected so far

        If m is in fact an PSCAffineModification,
        the old_vardim() of modification @a m must be
        identical to new_vardim() of this and
        old_rowdim() of modification @a m must be identical to
        new_rodim() of this.  The return value is the number of
        errors in this respect. If such occured, this incorporation
        is not performed.

        A general OracleModification @a m should only contain variable
        changes (this is not checked) and appending variables appends
        zero blocks.
    
"""
cb_incorporate!(self::CBPSCAffineModification, m::CBOracleModification) = @ccall libcb.cb_pscaffinemodification_incorporate(self.data::Ptr{Cvoid}, m.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_new_initial_oraclemodification(self::CBPSCAffineModification, old_var_dim::Integer)

returns a new object on the heap, that allows to incorporate this but starts off from a function whose input argument dimension is old_var_dim
"""
cb_new_initial_oraclemodification(self::CBPSCAffineModification, old_var_dim::Integer) = CBOracleModification(@ccall libcb.cb_pscaffinemodification_new_initial_oraclemodification(self.data::Ptr{Cvoid}, old_var_dim::Cint)::Ptr{Cvoid})

@doc raw"""
    cb_add_append_variables!(self::CBPSCAffineModification, append_dim::Integer)

append append_dim further variables at the end of the argument vector (without specifying their effect, they may e.g. be ignored by the function)
"""
cb_add_append_variables!(self::CBPSCAffineModification, append_dim::Integer) = @ccall libcb.cb_pscaffinemodification_add_append_variables(self.data::Ptr{Cvoid}, append_dim::Cint)::Cint

@doc raw"""
    cb_add_reassign_variables!(self::CBPSCAffineModification, new_dim::Integer, map_to_old_indices::Union{<:AbstractVector{Integer},Nothing})

reorder and resize the variables as given by the first new_dim entries of map_to_old_indices; each former index may appear at most once in this list, so new_dim < get_new_vardim()
"""
cb_add_reassign_variables!(self::CBPSCAffineModification, new_dim::Integer, map_to_old_indices::Union{<:AbstractVector{Integer},Nothing}) = GC.@preserve map_to_old_indices begin
    (LinearAlgebra.chkstride1(map_to_old_indices); @ccall libcb.cb_pscaffinemodification_add_reassign_variables(self.data::Ptr{Cvoid}, new_dim::Cint, map_to_old_indices::Ptr{Cint})::Cint)
end

@doc raw"""
    cb_set_cbout!(self::CBPSCAffineModification, incr::Integer = -1)

see CBout::set_cbout
"""
cb_set_cbout!(self::CBPSCAffineModification, incr::Integer = -1) = @ccall libcb.cb_pscaffinemodification_set_cbout(self.data::Ptr{Cvoid}, incr::Cint)::Cvoid

@doc raw"""
    cb_set_out!(self::CBPSCAffineModification, print_level::Integer = 1)

see CBout::set_out
"""
cb_set_out!(self::CBPSCAffineModification, print_level::Integer = 1) = @ccall libcb.cb_pscaffinemodification_set_out(self.data::Ptr{Cvoid}, print_level::Cint)::Cvoid

