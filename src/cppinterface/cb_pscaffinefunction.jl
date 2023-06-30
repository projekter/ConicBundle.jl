CBPSCAffineFunction(incr::Integer = -1) = CBPSCAffineFunction(@ccall libcb.cb_pscaffinefunction_new(incr::Cint)::Ptr{Cvoid})

@doc raw"""
    CBPSCAffineFunction(C::CBSparseCoeffmatMatrix, opAt::CBSparseCoeffmatMatrix, generating_primal::Union{<:CBPSCPrimal,Nothing} = nothing, incr::Integer = -1)

* @brief initialize the PSCAffineFunction with its matrices and possible a generating_primal

  C and opAt define the constant (block-)offset and the linear
  (block-)matrix function as described in the general text of
  PSCAffineFunction

  generating_primal defines in what form primal matrices should be
  aggregated. If the argument is NULL then no primal aggregation will
  take place.  The control over the generating primal is passed over to
  this. This will delete an existing generating primal whenever a new
  generating primal is set or upon destruction.

  The final two arguments allow to set the output, see CBout.
    
"""
CBPSCAffineFunction(C::CBSparseCoeffmatMatrix, opAt::CBSparseCoeffmatMatrix, generating_primal::Union{<:CBPSCPrimal,Nothing} = nothing, incr::Integer = -1) = CBPSCAffineFunction(@ccall libcb.cb_pscaffinefunction_new2(C.data::Ptr{Cvoid}, opAt.data::Ptr{Cvoid}, (isnothing(generating_primal) ? C_NULL : generating_primal.data)::Ptr{Cvoid}, incr::Cint)::Ptr{Cvoid})

@doc raw"""
    cb_set_check_correctness!(self::CBPSCAffineFunction, chk::Bool)

if set to true, ConicBundle employs some additional consistency checks
"""
cb_set_check_correctness!(self::CBPSCAffineFunction, chk::Bool) = @ccall libcb.cb_pscaffinefunction_set_check_correctness(self.data::Ptr{Cvoid}, chk::Cint)::Cvoid

@doc raw"""
    cb_set_max_Ritzvecs!(self::CBPSCAffineFunction, maxv::Integer)

set the maximum number of new Ritzvectors returned by evaluate(); values<1 default to 5
"""
cb_set_max_Ritzvecs!(self::CBPSCAffineFunction, maxv::Integer) = @ccall libcb.cb_pscaffinefunction_set_max_ritzvecs(self.data::Ptr{Cvoid}, maxv::Cint)::Cvoid

@doc raw"""
    cb_generate_minorant!(self::CBPSCAffineFunction, P::CBMatrix)

see PSCOracle::generate_minorant()
"""
cb_generate_minorant!(self::CBPSCAffineFunction, P::CBMatrix) = CBMinorant(@ccall libcb.cb_pscaffinefunction_generate_minorant(self.data::Ptr{Cvoid}, P.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_svec_projection!(self::CBPSCAffineFunction, svec_offset::CBMatrix, svec_coeffs::CBMatrix, P::CBMatrix, index_subset::Union{<:CBIndexmatrix,Nothing} = nothing)

see PSCOracle::svec_projection()
"""
cb_svec_projection!(self::CBPSCAffineFunction, svec_offset::CBMatrix, svec_coeffs::CBMatrix, P::CBMatrix, index_subset::Union{<:CBIndexmatrix,Nothing} = nothing) = @ccall libcb.cb_pscaffinefunction_svec_projection(self.data::Ptr{Cvoid}, svec_offset.data::Ptr{Cvoid}, svec_coeffs.data::Ptr{Cvoid}, P.data::Ptr{Cvoid}, (isnothing(index_subset) ? C_NULL : index_subset.data)::Ptr{Cvoid})::Cint

@doc raw"""
    cb_evaluate_projection!(self::CBPSCAffineFunction, current_point::CBMatrix, P::CBMatrix, relprec::Real, projected_Ritz_vectors::CBMatrix, projected_Ritz_values::CBMatrix)

see PSCOracle::evaluate_projection()
"""
cb_evaluate_projection!(self::CBPSCAffineFunction, current_point::CBMatrix, P::CBMatrix, relprec::Real, projected_Ritz_vectors::CBMatrix, projected_Ritz_values::CBMatrix) = @ccall libcb.cb_pscaffinefunction_evaluate_projection(self.data::Ptr{Cvoid}, current_point.data::Ptr{Cvoid}, P.data::Ptr{Cvoid}, relprec::Cdouble, projected_Ritz_vectors.data::Ptr{Cvoid}, projected_Ritz_values.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_left_right_product!(self::CBPSCAffineFunction, i::Integer, E::CBMatrix, F::CBMatrix, G::CBMatrix)

see PSCOracle::left_right_product()
"""
cb_left_right_product!(self::CBPSCAffineFunction, i::Integer, E::CBMatrix, F::CBMatrix, G::CBMatrix) = @ccall libcb.cb_pscaffinefunction_left_right_product(self.data::Ptr{Cvoid}, i::Cint, E.data::Ptr{Cvoid}, F.data::Ptr{Cvoid}, G.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_check_correctness(self::CBPSCAffineFunction)

see PSCOracle::check_correctness()
"""
cb_check_correctness(self::CBPSCAffineFunction) = Bool(@ccall libcb.cb_pscaffinefunction_check_correctness(self.data::Ptr{Cvoid})::Cint)

@doc raw"""
    cb_get_opAt!(self::CBPSCAffineFunction)

* returns the row representation of the coefficient matrices
     (each entry of the map represents a row by a SparseCoeffmatVector).
    
"""
cb_get_opAt!(self::CBPSCAffineFunction) = (@ccall libcb.cb_pscaffinefunction_get_opat(self.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

@doc raw"""
    cb_get_C!(self::CBPSCAffineFunction)

* returns the block representation of the coefficient matrices
     (each entry of the map represents a block by a SparseCoeffmatVector).
    
"""
cb_get_C!(self::CBPSCAffineFunction) = (@ccall libcb.cb_pscaffinefunction_get_c(self.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

@doc raw"""
    cb_get_generating_primal!(self::CBPSCAffineFunction)

returns the current setting concerning the generation of an PSCPrimal (0 for none)
"""
cb_get_generating_primal!(self::CBPSCAffineFunction) = CBPSCPrimal(@ccall libcb.cb_pscaffinefunction_get_generating_primal(self.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_set_out!(self::CBPSCAffineFunction, pril::Integer = 1)

see ConicBundle::CBout
"""
cb_set_out!(self::CBPSCAffineFunction, pril::Integer = 1) = @ccall libcb.cb_pscaffinefunction_set_out(self.data::Ptr{Cvoid}, pril::Cint)::Cvoid

@doc raw"""
    cb_set_cbout!(self::CBPSCAffineFunction, incr::Integer = -1)

see ConicBundle::CBout
"""
cb_set_cbout!(self::CBPSCAffineFunction, incr::Integer = -1) = @ccall libcb.cb_pscaffinefunction_set_cbout(self.data::Ptr{Cvoid}, incr::Cint)::Cvoid

@doc raw"""
    cb_print_problem_data(self::CBPSCAffineFunction)

write the problem description to out so that it can be read again by read_problem_data()
"""
cb_print_problem_data(self::CBPSCAffineFunction) = @ccall libcb.cb_pscaffinefunction_print_problem_data(self.data::Ptr{Cvoid})::Cvoid

@doc raw"""
    cb_print_problem_data_to_mfile(self::CBPSCAffineFunction, blocknr::Integer)

undocumented highly volatile variant for external testing
"""
cb_print_problem_data_to_mfile(self::CBPSCAffineFunction, blocknr::Integer) = @ccall libcb.cb_pscaffinefunction_print_problem_data_to_mfile(self.data::Ptr{Cvoid}, blocknr::Cint)::Cvoid

