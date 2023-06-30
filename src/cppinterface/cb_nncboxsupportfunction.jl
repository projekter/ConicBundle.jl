@doc raw"""
    CBNNCBoxSupportFunction(lb::CBMatrix, ub::CBMatrix, incr::Integer = -1)

intialize with lower bounds vector lb and upper bounds vector ub (and output options), both must be column vectors of the same length and lb<=ub componentwise, length 0 results in objective value 0
"""
CBNNCBoxSupportFunction(lb::CBMatrix, ub::CBMatrix, incr::Integer = -1) = CBNNCBoxSupportFunction(@ccall libcb.cb_nncboxsupportfunction_new(lb.data::Ptr{Cvoid}, ub.data::Ptr{Cvoid}, incr::Cint)::Ptr{Cvoid})

@doc raw"""
    cb_check_correctness(self::CBNNCBoxSupportFunction)

see MatrixFunctionOracle::check_correctness()   (true only needed for debugging)
"""
cb_check_correctness(self::CBNNCBoxSupportFunction) = Bool(@ccall libcb.cb_nncboxsupportfunction_check_correctness(self.data::Ptr{Cvoid})::Cint)

@doc raw"""
    cb_get_lower_bounds!(self::CBNNCBoxSupportFunction)

returns the column vector of lower bounds
"""
cb_get_lower_bounds!(self::CBNNCBoxSupportFunction) = (@ccall libcb.cb_nncboxsupportfunction_get_lower_bounds(self.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

@doc raw"""
    cb_get_upper_bounds!(self::CBNNCBoxSupportFunction)

retunrs the column vector of upper bounds
"""
cb_get_upper_bounds!(self::CBNNCBoxSupportFunction) = (@ccall libcb.cb_nncboxsupportfunction_get_upper_bounds(self.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

@doc raw"""
    cb_set_out!(self::CBNNCBoxSupportFunction, pril::Integer = 1)

see ConicBundle::CBout
"""
cb_set_out!(self::CBNNCBoxSupportFunction, pril::Integer = 1) = @ccall libcb.cb_nncboxsupportfunction_set_out(self.data::Ptr{Cvoid}, pril::Cint)::Cvoid

@doc raw"""
    cb_set_cbout!(self::CBNNCBoxSupportFunction, incr::Integer = -1)

see ConicBundle::CBout
"""
cb_set_cbout!(self::CBNNCBoxSupportFunction, incr::Integer = -1) = @ccall libcb.cb_nncboxsupportfunction_set_cbout(self.data::Ptr{Cvoid}, incr::Cint)::Cvoid

@doc raw"""
    cb_print_problem_data(self::CBNNCBoxSupportFunction)

write the problem description to out so that it can be read again by read_problem_data()
"""
cb_print_problem_data(self::CBNNCBoxSupportFunction) = @ccall libcb.cb_nncboxsupportfunction_print_problem_data(self.data::Ptr{Cvoid})::Cvoid

@doc raw"""
    cb_print_problem_data_to_mfile(self::CBNNCBoxSupportFunction, blocknr::Integer)

undocumented highly volatile variant for external testing
"""
cb_print_problem_data_to_mfile(self::CBNNCBoxSupportFunction, blocknr::Integer) = @ccall libcb.cb_nncboxsupportfunction_print_problem_data_to_mfile(self.data::Ptr{Cvoid}, blocknr::Cint)::Cvoid

