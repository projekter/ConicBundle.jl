@doc raw"""
    CBGroundsetModification(var_olddim::Integer = 0, incr::Integer = -1)

constructor, calls modification constructor
"""
CBGroundsetModification(var_olddim::Integer = 0, incr::Integer = -1) = CBGroundsetModification(@ccall libcb.cb_groundsetmodification_new(var_olddim::Cint, incr::Cint)::Ptr{Cvoid})

@doc raw"""
    cb_clear!(self::CBGroundsetModification, var_olddim::Integer)

reset modifications to an unmodified object currently having var_olddim variables, calls Modification::clear
"""
cb_clear!(self::CBGroundsetModification, var_olddim::Integer) = @ccall libcb.cb_groundsetmodification_clear(self.data::Ptr{Cvoid}, var_olddim::Cint)::Cvoid

@doc raw"""
    cb_add_append_vars!(self::CBGroundsetModification, append_dim::Integer, start_val::Union{<:CBMatrix,Nothing} = nothing, costs::Union{<:CBMatrix,Nothing} = nothing)

append @a append_dim new variables with start_val as initial values (if NULL, use default value), calls Modification::add_append_vars
"""
cb_add_append_vars!(self::CBGroundsetModification, append_dim::Integer, start_val::Union{<:CBMatrix,Nothing} = nothing, costs::Union{<:CBMatrix,Nothing} = nothing) = @ccall libcb.cb_groundsetmodification_add_append_vars(self.data::Ptr{Cvoid}, append_dim::Cint, (isnothing(start_val) ? C_NULL : start_val.data)::Ptr{Cvoid}, (isnothing(costs) ? C_NULL : costs.data)::Ptr{Cvoid})::Cint

@doc raw"""
    cb_add_reassign_vars!(self::CBGroundsetModification, map_to_old::CBIndexmatrix)

reassign the variables as given in @a map_to_old, calls Modification::add_reassign_vars
"""
cb_add_reassign_vars!(self::CBGroundsetModification, map_to_old::CBIndexmatrix) = @ccall libcb.cb_groundsetmodification_add_reassign_vars(self.data::Ptr{Cvoid}, map_to_old.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_add_delete_vars!(self::CBGroundsetModification, del_ind::CBIndexmatrix, map_to_old::CBIndexmatrix)

delete the variables indexed by @a del_ind; for each new index @a map_to_old returns the old one; calls Modification::add_delete_vars
"""
cb_add_delete_vars!(self::CBGroundsetModification, del_ind::CBIndexmatrix, map_to_old::CBIndexmatrix) = @ccall libcb.cb_groundsetmodification_add_delete_vars(self.data::Ptr{Cvoid}, del_ind.data::Ptr{Cvoid}, map_to_old.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_add_offset!(self::CBGroundsetModification, delta::Real)

add to the current offset the value @a delta
"""
cb_add_offset!(self::CBGroundsetModification, delta::Real) = @ccall libcb.cb_groundsetmodification_add_offset(self.data::Ptr{Cvoid}, delta::Cdouble)::Cint

@doc raw"""
    cb_apply_to_vars(self::CBGroundsetModification, vars::CBMatrix)

carry out the collected modifications on the given vector by calling Modification::apply_to_vars
"""
cb_apply_to_vars(self::CBGroundsetModification, vars::CBMatrix) = @ccall libcb.cb_groundsetmodification_apply_to_vars(self.data::Ptr{Cvoid}, vars.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_apply_to_costs(self::CBGroundsetModification, costs::CBMatrix)

carry out the collected modifications on the given vector by calling Modification::apply_to_vars
"""
function cb_apply_to_costs(self::CBGroundsetModification, costs::CBMatrix)
    in_offset = Ref{Float64}()
    @ccall libcb.cb_groundsetmodification_apply_to_costs(self.data::Ptr{Cvoid}, costs.data::Ptr{Cvoid}, in_offset::Ref{Float64})::Cint
    return in_offset[]
end

@doc raw"""
    cb_no_modification(self::CBGroundsetModification)

returns true if no modifications need to be executed except possibly an offset change for the ground set minorant
"""
cb_no_modification(self::CBGroundsetModification) = Bool(@ccall libcb.cb_groundsetmodification_no_modification(self.data::Ptr{Cvoid})::Cint)

@doc raw"""
    cb_set_append_to_old!(self::CBGroundsetModification, append_only::Bool)

if set to true, no deletions/reassignments may be present or specified in the future, only appensions are allowed
"""
cb_set_append_to_old!(self::CBGroundsetModification, append_only::Bool) = @ccall libcb.cb_groundsetmodification_set_append_to_old(self.data::Ptr{Cvoid}, append_only::Cint)::Cint

@doc raw"""
    cb_append_to_old(self::CBGroundsetModification)

returns true if this only contains appending operations and incorporating this is done with respect to the old dimension
"""
cb_append_to_old(self::CBGroundsetModification) = Bool(@ccall libcb.cb_groundsetmodification_append_to_old(self.data::Ptr{Cvoid})::Cint)

@doc raw"""
    cb_no_additions_or_deletions_in_vars(self::CBGroundsetModification)

returns true if no variables were added or deleted (allows permutations), false otherwise
"""
cb_no_additions_or_deletions_in_vars(self::CBGroundsetModification) = Bool(@ccall libcb.cb_groundsetmodification_no_additions_or_deletions_in_vars(self.data::Ptr{Cvoid})::Cint)

@doc raw"""
    cb_deleted_variables_are_zero(self::CBGroundsetModification, oldpoint::CBMatrix)

returns true if all entries deleted in @a oldpoint (must be a vector of length old_vardim()) are 0,  false otherwise
"""
cb_deleted_variables_are_zero(self::CBGroundsetModification, oldpoint::CBMatrix) = Bool(@ccall libcb.cb_groundsetmodification_deleted_variables_are_zero(self.data::Ptr{Cvoid}, oldpoint.data::Ptr{Cvoid})::Cint)

@doc raw"""
    cb_new_variables_are_zero(self::CBGroundsetModification, newpoint::CBMatrix)

returns true if all entries in newpoint (must be a vector of length new_vardim()) that correspond to new variables have value 0 and false otherwise
"""
cb_new_variables_are_zero(self::CBGroundsetModification, newpoint::CBMatrix) = Bool(@ccall libcb.cb_groundsetmodification_new_variables_are_zero(self.data::Ptr{Cvoid}, newpoint.data::Ptr{Cvoid})::Cint)

@doc raw"""
    cb_mapped_variables_are_equal(self::CBGroundsetModification, newpoint::CBMatrix, oldpoint::CBMatrix)

returns true if the values in newpoint (must be a vector of length new_vardim()) that correspond to old variables match the old values stored in oldpoint (must be a vector of length old_vardim()) and false otherwise
"""
cb_mapped_variables_are_equal(self::CBGroundsetModification, newpoint::CBMatrix, oldpoint::CBMatrix) = Bool(@ccall libcb.cb_groundsetmodification_mapped_variables_are_equal(self.data::Ptr{Cvoid}, newpoint.data::Ptr{Cvoid}, oldpoint.data::Ptr{Cvoid})::Cint)

@doc raw"""
    cb_get_old_vardim(self::CBGroundsetModification)

returns the number of variables before modification
"""
cb_get_old_vardim(self::CBGroundsetModification) = @ccall libcb.cb_groundsetmodification_get_old_vardim(self.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_get_new_vardim(self::CBGroundsetModification)

returns the number of variables once all stored modifications have been performed
"""
cb_get_new_vardim(self::CBGroundsetModification) = @ccall libcb.cb_groundsetmodification_get_new_vardim(self.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_get_appended_vardim(self::CBGroundsetModification)

returns the number of variables that are appended (due to later reassignmentds they may no longer be located at the end)
"""
cb_get_appended_vardim(self::CBGroundsetModification) = @ccall libcb.cb_groundsetmodification_get_appended_vardim(self.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_get_map_to_old_variables(self::CBGroundsetModification)

returns null if there are no index changes, otherwise the Indexmatrix pointed to is a vector whose i-th entry holds the old index of the new i-th variable (injective!), index values exceeding old_vardim() refer to newly appended variables
"""
cb_get_map_to_old_variables(self::CBGroundsetModification) = @ccall libcb.cb_groundsetmodification_get_map_to_old_variables(self.data::Ptr{Cvoid})::Ptr{Cint}

@doc raw"""
    cb_get_add_offset(self::CBGroundsetModification)

returns the change in the offste value of the groundset minorant
"""
cb_get_add_offset(self::CBGroundsetModification) = @ccall libcb.cb_groundsetmodification_get_add_offset(self.data::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_get_append_costs(self::CBGroundsetModification)

returns the change in the offste value of the groundset minorant
"""
cb_get_append_costs(self::CBGroundsetModification) = CBMatrix(@ccall libcb.cb_groundsetmodification_get_append_costs(self.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_incorporate!(self::CBGroundsetModification, m::CBOracleModification)

incorporate the OracleModification @a m (it should only contain variable changes, but this is not checked!) into this one; calls Modification::incorporate
"""
cb_incorporate!(self::CBGroundsetModification, m::CBOracleModification) = @ccall libcb.cb_groundsetmodification_incorporate(self.data::Ptr{Cvoid}, m.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_new_initial_oraclemodification(self::CBGroundsetModification, old_var_dim::Integer)

returns a new object on the heap, that allows to incorporate this but starts off from a function whose input argument dimension is old_var_dim
"""
cb_new_initial_oraclemodification(self::CBGroundsetModification, old_var_dim::Integer) = CBOracleModification(@ccall libcb.cb_groundsetmodification_new_initial_oraclemodification(self.data::Ptr{Cvoid}, old_var_dim::Cint)::Ptr{Cvoid})

@doc raw"""
    cb_add_append_variables!(self::CBGroundsetModification, append_dim::Integer)

append append_dim further variables at the end of the argument vector (without specifying their effect, they may e.g. be ignored by the function)
"""
cb_add_append_variables!(self::CBGroundsetModification, append_dim::Integer) = @ccall libcb.cb_groundsetmodification_add_append_variables(self.data::Ptr{Cvoid}, append_dim::Cint)::Cint

@doc raw"""
    cb_add_reassign_variables!(self::CBGroundsetModification, new_dim::Integer, map_to_old_indices::Union{<:AbstractVector{Integer},Nothing})

reorder and resize the variables as given by the first new_dim entries of map_to_old_indices; each former index may appear at most once in this list, so new_dim < get_new_vardim()
"""
cb_add_reassign_variables!(self::CBGroundsetModification, new_dim::Integer, map_to_old_indices::Union{<:AbstractVector{Integer},Nothing}) = GC.@preserve map_to_old_indices begin
    (LinearAlgebra.chkstride1(map_to_old_indices); @ccall libcb.cb_groundsetmodification_add_reassign_variables(self.data::Ptr{Cvoid}, new_dim::Cint, map_to_old_indices::Ptr{Cint})::Cint)
end

@doc raw"""
    cb_old_vardim(self::CBGroundsetModification)

returns the number of variables before modification (given on initialization)
"""
cb_old_vardim(self::CBGroundsetModification) = @ccall libcb.cb_groundsetmodification_old_vardim(self.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_new_vardim(self::CBGroundsetModification)

returns the number of variables once all stored modifications have been performed
"""
cb_new_vardim(self::CBGroundsetModification) = @ccall libcb.cb_groundsetmodification_new_vardim(self.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_appended_vardim(self::CBGroundsetModification)

returns the number of variables that are appended (due to later reassignmentds they may no longer be located at the end)
"""
cb_appended_vardim(self::CBGroundsetModification) = @ccall libcb.cb_groundsetmodification_appended_vardim(self.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_map_to_old_variables(self::CBGroundsetModification)

returns null if there are no index changes, otherwise the Indexmatrix pointed to is a vector whose i-th entry holds the old index of the new i-th variable (injective!), index values exceeding old_vardim() refer to newly appended variables
"""
cb_map_to_old_variables(self::CBGroundsetModification) = CBIndexmatrix(@ccall libcb.cb_groundsetmodification_map_to_old_variables(self.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_deleted_var_indices(self::CBGroundsetModification)

returns null if there were no deletions, otherwise the Indexmatrix pointed to is a vector holding the deleted old variable indices in increasing order
"""
cb_deleted_var_indices(self::CBGroundsetModification) = CBIndexmatrix(@ccall libcb.cb_groundsetmodification_deleted_var_indices(self.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_new_var_indices(self::CBGroundsetModification)

returns null if no variables were added, otherwise the Indexmatrix pointed ato is a vector holding the new indices of the new variables in increasing order
"""
cb_new_var_indices(self::CBGroundsetModification) = CBIndexmatrix(@ccall libcb.cb_groundsetmodification_new_var_indices(self.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_set_out!(self::CBGroundsetModification, print_level::Integer = 1)

see CBout::set_out
"""
cb_set_out!(self::CBGroundsetModification, print_level::Integer = 1) = @ccall libcb.cb_groundsetmodification_set_out(self.data::Ptr{Cvoid}, print_level::Cint)::Cvoid

@doc raw"""
    cb_set_cbout!(self::CBGroundsetModification, incr::Integer = -1)

see CBout::set_out
"""
cb_set_cbout!(self::CBGroundsetModification, incr::Integer = -1) = @ccall libcb.cb_groundsetmodification_set_cbout(self.data::Ptr{Cvoid}, incr::Cint)::Cvoid

