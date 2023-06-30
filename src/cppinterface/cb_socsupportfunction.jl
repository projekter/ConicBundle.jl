@doc raw"""
    CBSOCSupportFunction(socdim::Integer, incr::Integer = -1)

initialize with dimension >= 1 (and output options)
"""
CBSOCSupportFunction(socdim::Integer, incr::Integer = -1) = CBSOCSupportFunction(@ccall libcb.cb_socsupportfunction_new(socdim::Cint, incr::Cint)::Ptr{Cvoid})

@doc raw"""
    cb_generate_minorant!(self::CBSOCSupportFunction, SOCvec::CBMatrix)

see SOCOracle::generate_minorant()
"""
cb_generate_minorant!(self::CBSOCSupportFunction, SOCvec::CBMatrix) = CBMinorant(@ccall libcb.cb_socsupportfunction_generate_minorant(self.data::Ptr{Cvoid}, SOCvec.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_extract_SOCvector!(self::CBSOCSupportFunction, SOCvec::CBMatrix, SOCminorant::Union{<:CBMinorant,Nothing})

see SOCOracle::extract_SOCvector()
"""
cb_extract_SOCvector!(self::CBSOCSupportFunction, SOCvec::CBMatrix, SOCminorant::Union{<:CBMinorant,Nothing}) = @ccall libcb.cb_socsupportfunction_extract_socvector(self.data::Ptr{Cvoid}, SOCvec.data::Ptr{Cvoid}, (isnothing(SOCminorant) ? C_NULL : SOCminorant.data)::Ptr{Cvoid})::Cint

@doc raw"""
    cb_projection!(self::CBSOCSupportFunction, offset::CBMatrix, coeffs::CBMatrix, bar_P::CBMatrix, index_subset::Union{<:CBIndexmatrix,Nothing} = nothing)

see SOCOracle::projection()
"""
cb_projection!(self::CBSOCSupportFunction, offset::CBMatrix, coeffs::CBMatrix, bar_P::CBMatrix, index_subset::Union{<:CBIndexmatrix,Nothing} = nothing) = @ccall libcb.cb_socsupportfunction_projection(self.data::Ptr{Cvoid}, offset.data::Ptr{Cvoid}, coeffs.data::Ptr{Cvoid}, bar_P.data::Ptr{Cvoid}, (isnothing(index_subset) ? C_NULL : index_subset.data)::Ptr{Cvoid})::Cint

@doc raw"""
    cb_evaluate_projection!(self::CBSOCSupportFunction, current_point::CBMatrix, P::CBMatrix, relprec::Real)

see SOCOracle::evaluate()
"""
function cb_evaluate_projection!(self::CBSOCSupportFunction, current_point::CBMatrix, P::CBMatrix, relprec::Real)
    projected_SOC_value = Ref{Float64}()
    @ccall libcb.cb_socsupportfunction_evaluate_projection(self.data::Ptr{Cvoid}, current_point.data::Ptr{Cvoid}, P.data::Ptr{Cvoid}, relprec::Cdouble, projected_SOC_value::Ref{Float64})::Cint
    return projected_SOC_value[]
end

@doc raw"""
    cb_check_correctness(self::CBSOCSupportFunction)

see SOCOracle::check_correctness() (true only needed for debugging)
"""
cb_check_correctness(self::CBSOCSupportFunction) = Bool(@ccall libcb.cb_socsupportfunction_check_correctness(self.data::Ptr{Cvoid})::Cint)

@doc raw"""
    cb_get_socdim!(self::CBSOCSupportFunction)

returns the dimension of the second order cone
"""
cb_get_socdim!(self::CBSOCSupportFunction) = @ccall libcb.cb_socsupportfunction_get_socdim(self.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_set_out!(self::CBSOCSupportFunction, pril::Integer = 1)

see ConicBundle::CBout
"""
cb_set_out!(self::CBSOCSupportFunction, pril::Integer = 1) = @ccall libcb.cb_socsupportfunction_set_out(self.data::Ptr{Cvoid}, pril::Cint)::Cvoid

@doc raw"""
    cb_set_cbout!(self::CBSOCSupportFunction, incr::Integer = -1)

see ConicBundle::CBout
"""
cb_set_cbout!(self::CBSOCSupportFunction, incr::Integer = -1) = @ccall libcb.cb_socsupportfunction_set_cbout(self.data::Ptr{Cvoid}, incr::Cint)::Cvoid

@doc raw"""
    cb_print_problem_data(self::CBSOCSupportFunction)

write the problem description to out so that it can be read again by read_problem_data()
"""
cb_print_problem_data(self::CBSOCSupportFunction) = @ccall libcb.cb_socsupportfunction_print_problem_data(self.data::Ptr{Cvoid})::Cvoid

@doc raw"""
    cb_print_problem_data_to_mfile(self::CBSOCSupportFunction, blocknr::Integer)

undocumented highly volatile variant for external testing
"""
cb_print_problem_data_to_mfile(self::CBSOCSupportFunction, blocknr::Integer) = @ccall libcb.cb_socsupportfunction_print_problem_data_to_mfile(self.data::Ptr{Cvoid}, blocknr::Cint)::Cvoid

