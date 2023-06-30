@doc raw"""
    CBSOCSupportModification(var_olddim::Integer = 0, incr::Integer = -1)

constructor, calls modification constructor
"""
CBSOCSupportModification(var_olddim::Integer = 0, incr::Integer = -1) = CBSOCSupportModification(@ccall libcb.cb_socsupportmodification_new(var_olddim::Cint, incr::Cint)::Ptr{Cvoid})

@doc raw"""
    cb_clear!(self::CBSOCSupportModification, var_olddim::Integer)

reset modifications to an unmodified object currently having var_olddim variables, calls Modification::clear
"""
cb_clear!(self::CBSOCSupportModification, var_olddim::Integer) = @ccall libcb.cb_socsupportmodification_clear(self.data::Ptr{Cvoid}, var_olddim::Cint)::Cvoid

@doc raw"""
    cb_add_append_vars!(self::CBSOCSupportModification, append_dim::Integer)

append @a append_dim new variables to the box function
"""
cb_add_append_vars!(self::CBSOCSupportModification, append_dim::Integer) = @ccall libcb.cb_socsupportmodification_add_append_vars(self.data::Ptr{Cvoid}, append_dim::Cint)::Cint

@doc raw"""
    cb_add_reassign_vars!(self::CBSOCSupportModification, map_to_old::CBIndexmatrix)

reassign the variables as given in @a map_to_old, calls Modification::add_reassign_vars
"""
cb_add_reassign_vars!(self::CBSOCSupportModification, map_to_old::CBIndexmatrix) = @ccall libcb.cb_socsupportmodification_add_reassign_vars(self.data::Ptr{Cvoid}, map_to_old.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_add_delete_vars!(self::CBSOCSupportModification, del_ind::CBIndexmatrix, map_to_old::CBIndexmatrix)

delete the variables indexed by @a del_ind; for each new index @a map_to_old returns the old one; calls Modification::add_delete_vars
"""
cb_add_delete_vars!(self::CBSOCSupportModification, del_ind::CBIndexmatrix, map_to_old::CBIndexmatrix) = @ccall libcb.cb_socsupportmodification_add_delete_vars(self.data::Ptr{Cvoid}, del_ind.data::Ptr{Cvoid}, map_to_old.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_no_modification(self::CBSOCSupportModification)

returns true if no modifications need to be executed except possibly an offset change for the ground set minorant
"""
cb_no_modification(self::CBSOCSupportModification) = Bool(@ccall libcb.cb_socsupportmodification_no_modification(self.data::Ptr{Cvoid})::Cint)

@doc raw"""
    cb_set_append_to_old!(self::CBSOCSupportModification, append_only::Bool)

if set to true, no deletions/reassignments may be present or specified in the future, only appensions are allowed
"""
cb_set_append_to_old!(self::CBSOCSupportModification, append_only::Bool) = @ccall libcb.cb_socsupportmodification_set_append_to_old(self.data::Ptr{Cvoid}, append_only::Cint)::Cint

@doc raw"""
    cb_append_to_old(self::CBSOCSupportModification)

returns true if this only contains appending operations and incorporating this is done with respect to the old dimension
"""
cb_append_to_old(self::CBSOCSupportModification) = Bool(@ccall libcb.cb_socsupportmodification_append_to_old(self.data::Ptr{Cvoid})::Cint)

@doc raw"""
    cb_no_additions_or_deletions_in_vars(self::CBSOCSupportModification)

returns true if no variables were added or deleted (allows permutations), false otherwise
"""
cb_no_additions_or_deletions_in_vars(self::CBSOCSupportModification) = Bool(@ccall libcb.cb_socsupportmodification_no_additions_or_deletions_in_vars(self.data::Ptr{Cvoid})::Cint)

@doc raw"""
    cb_deleted_variables_are_zero(self::CBSOCSupportModification, oldpoint::CBMatrix)

returns true if all entries deleted in @a oldpoint (must be a vector of length old_vardim()) are 0,  false otherwise
"""
cb_deleted_variables_are_zero(self::CBSOCSupportModification, oldpoint::CBMatrix) = Bool(@ccall libcb.cb_socsupportmodification_deleted_variables_are_zero(self.data::Ptr{Cvoid}, oldpoint.data::Ptr{Cvoid})::Cint)

@doc raw"""
    cb_new_variables_are_zero(self::CBSOCSupportModification, newpoint::CBMatrix)

returns true if all entries in newpoint (must be a vector of length new_vardim()) that correspond to new variables have value 0 and false otherwise
"""
cb_new_variables_are_zero(self::CBSOCSupportModification, newpoint::CBMatrix) = Bool(@ccall libcb.cb_socsupportmodification_new_variables_are_zero(self.data::Ptr{Cvoid}, newpoint.data::Ptr{Cvoid})::Cint)

@doc raw"""
    cb_mapped_variables_are_equal(self::CBSOCSupportModification, newpoint::CBMatrix, oldpoint::CBMatrix)

returns true if the values in newpoint (must be a vector of length new_vardim()) that correspond to old variables match the old values stored in oldpoint (must be a vector of length old_vardim()) and false otherwise
"""
cb_mapped_variables_are_equal(self::CBSOCSupportModification, newpoint::CBMatrix, oldpoint::CBMatrix) = Bool(@ccall libcb.cb_socsupportmodification_mapped_variables_are_equal(self.data::Ptr{Cvoid}, newpoint.data::Ptr{Cvoid}, oldpoint.data::Ptr{Cvoid})::Cint)

@doc raw"""
    cb_get_old_vardim(self::CBSOCSupportModification)

returns the number of variables before modification
"""
cb_get_old_vardim(self::CBSOCSupportModification) = @ccall libcb.cb_socsupportmodification_get_old_vardim(self.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_get_new_vardim(self::CBSOCSupportModification)

returns the number of variables once all stored modifications have been performed
"""
cb_get_new_vardim(self::CBSOCSupportModification) = @ccall libcb.cb_socsupportmodification_get_new_vardim(self.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_get_appended_vardim(self::CBSOCSupportModification)

returns the number of variables that are appended (due to later reassignmentds they may no longer be located at the end)
"""
cb_get_appended_vardim(self::CBSOCSupportModification) = @ccall libcb.cb_socsupportmodification_get_appended_vardim(self.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_get_map_to_old_variables(self::CBSOCSupportModification)

returns null if there are no index changes, otherwise the Indexmatrix pointed to is a vector whose i-th entry holds the old index of the new i-th variable (injective!), index values exceeding old_vardim() refer to newly appended variables
"""
cb_get_map_to_old_variables(self::CBSOCSupportModification) = @ccall libcb.cb_socsupportmodification_get_map_to_old_variables(self.data::Ptr{Cvoid})::Ptr{Cint}

@doc raw"""
    cb_incorporate!(self::CBSOCSupportModification, m::CBOracleModification)

incorporate the OracleModification @a m (it should only contain variable changes, but this is not checked!) into this one; calls Modification::incorporate
"""
cb_incorporate!(self::CBSOCSupportModification, m::CBOracleModification) = @ccall libcb.cb_socsupportmodification_incorporate(self.data::Ptr{Cvoid}, m.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_new_initial_oraclemodification(self::CBSOCSupportModification, old_var_dim::Integer)

returns a new object on the heap, that allows to incorporate this but starts off from a function whose input argument dimension is old_var_dim
"""
cb_new_initial_oraclemodification(self::CBSOCSupportModification, old_var_dim::Integer) = CBOracleModification(@ccall libcb.cb_socsupportmodification_new_initial_oraclemodification(self.data::Ptr{Cvoid}, old_var_dim::Cint)::Ptr{Cvoid})

@doc raw"""
    cb_add_append_variables!(self::CBSOCSupportModification, append_dim::Integer)

append append_dim further variables at the end of the argument vector (without effect on the box function, i.e., lower and upper bounds are set to 0)
"""
cb_add_append_variables!(self::CBSOCSupportModification, append_dim::Integer) = @ccall libcb.cb_socsupportmodification_add_append_variables(self.data::Ptr{Cvoid}, append_dim::Cint)::Cint

@doc raw"""
    cb_add_reassign_variables!(self::CBSOCSupportModification, new_dim::Integer, map_to_old_indices::Union{<:AbstractVector{Integer},Nothing})

reorder and resize the variables as given by the first new_dim entries of map_to_old_indices; each former index may appear at most once in this list, so new_dim < get_new_vardim()
"""
cb_add_reassign_variables!(self::CBSOCSupportModification, new_dim::Integer, map_to_old_indices::Union{<:AbstractVector{Integer},Nothing}) = GC.@preserve map_to_old_indices begin
    (LinearAlgebra.chkstride1(map_to_old_indices); @ccall libcb.cb_socsupportmodification_add_reassign_variables(self.data::Ptr{Cvoid}, new_dim::Cint, map_to_old_indices::Ptr{Cint})::Cint)
end

@doc raw"""
    cb_old_vardim(self::CBSOCSupportModification)

returns the number of variables before modification (given on initialization)
"""
cb_old_vardim(self::CBSOCSupportModification) = @ccall libcb.cb_socsupportmodification_old_vardim(self.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_new_vardim(self::CBSOCSupportModification)

returns the number of variables once all stored modifications have been performed
"""
cb_new_vardim(self::CBSOCSupportModification) = @ccall libcb.cb_socsupportmodification_new_vardim(self.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_appended_vardim(self::CBSOCSupportModification)

returns the number of variables that are appended (due to later reassignmentds they may no longer be located at the end)
"""
cb_appended_vardim(self::CBSOCSupportModification) = @ccall libcb.cb_socsupportmodification_appended_vardim(self.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_map_to_old_variables(self::CBSOCSupportModification)

returns null if there are no index changes, otherwise the Indexmatrix pointed to is a vector whose i-th entry holds the old index of the new i-th variable (injective!), index values exceeding old_vardim() refer to newly appended variables
"""
cb_map_to_old_variables(self::CBSOCSupportModification) = CBIndexmatrix(@ccall libcb.cb_socsupportmodification_map_to_old_variables(self.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_deleted_var_indices(self::CBSOCSupportModification)

returns null if there were no deletions, otherwise the Indexmatrix pointed to is a vector holding the deleted old variable indices in increasing order
"""
cb_deleted_var_indices(self::CBSOCSupportModification) = CBIndexmatrix(@ccall libcb.cb_socsupportmodification_deleted_var_indices(self.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_new_var_indices(self::CBSOCSupportModification)

returns null if no variables were added, otherwise the Indexmatrix pointed ato is a vector holding the new indices of the new variables in increasing order
"""
cb_new_var_indices(self::CBSOCSupportModification) = CBIndexmatrix(@ccall libcb.cb_socsupportmodification_new_var_indices(self.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_set_out!(self::CBSOCSupportModification, print_level::Integer = 1)

see CBout::set_out
"""
cb_set_out!(self::CBSOCSupportModification, print_level::Integer = 1) = @ccall libcb.cb_socsupportmodification_set_out(self.data::Ptr{Cvoid}, print_level::Cint)::Cvoid

@doc raw"""
    cb_set_cbout!(self::CBSOCSupportModification, incr::Integer = -1)

see CBout::set_out
"""
cb_set_cbout!(self::CBSOCSupportModification, incr::Integer = -1) = @ccall libcb.cb_socsupportmodification_set_cbout(self.data::Ptr{Cvoid}, incr::Cint)::Cvoid

