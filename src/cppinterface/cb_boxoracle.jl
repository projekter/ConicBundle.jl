@doc raw"""
    CBBoxOracle(in_lb::CBMatrix, in_ub::CBMatrix)

constructor initializing lower and upper bounds (must have the same dimesnion, not checked)
"""
CBBoxOracle(in_lb::CBMatrix, in_ub::CBMatrix) = CBBoxOracle(@ccall libcb.cb_boxoracle_new(in_lb.data::Ptr{Cvoid}, in_ub.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_get_lower_bounds!(self::CBBoxOracle)

* @brief returns the lower bounds vector of the box
     
"""
cb_get_lower_bounds!(self::CBBoxOracle) = (@ccall libcb.cb_boxoracle_get_lower_bounds(self.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

@doc raw"""
    cb_get_upper_bounds!(self::CBBoxOracle)

* @brief returns the upper bounds vector of the box
     
"""
cb_get_upper_bounds!(self::CBBoxOracle) = (@ccall libcb.cb_boxoracle_get_upper_bounds(self.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

@doc raw"""
    cb_check_correctness(self::CBBoxOracle)

*@brief switch on/off some correctnes checks on the oracle 
"""
cb_check_correctness(self::CBBoxOracle) = Bool(@ccall libcb.cb_boxoracle_check_correctness(self.data::Ptr{Cvoid})::Cint)

