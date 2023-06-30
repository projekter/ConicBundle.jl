@doc raw"""
    CBPrimalMatrix()

empty matrix
"""
CBPrimalMatrix() = CBPrimalMatrix(@ccall libcb.cb_primalmatrix_new()::Ptr{Cvoid})

@doc raw"""
    CBPrimalMatrix(nr::Integer, nc::Integer)

* @brief generate a matrix of size nr x nc but WITHOUT initializing the memory

        If initializing the memory externally and CONICBUNDLE_DEBUG is defined, please use
        set_init() via matrix.set_init(true) in order to avoid warnings concerning improper
        initialization
    
"""
CBPrimalMatrix(nr::Integer, nc::Integer) = CBPrimalMatrix(@ccall libcb.cb_primalmatrix_new2(nr::Cint, nc::Cint)::Ptr{Cvoid})

@doc raw"""
    CBPrimalMatrix(r::Integer, c::Integer, d::Real)

generate a matrix of size nr x nc initializing all elements to the value d
"""
CBPrimalMatrix(r::Integer, c::Integer, d::Real) = CBPrimalMatrix(@ccall libcb.cb_primalmatrix_new3(r::Cint, c::Cint, d::Cdouble)::Ptr{Cvoid})

@doc raw"""
    CBPrimalMatrix(pm::CBPrimalMatrix)

copy constructor, *this=pm
"""
CBPrimalMatrix(pm::CBPrimalMatrix) = CBPrimalMatrix(@ccall libcb.cb_primalmatrix_new4(pm.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    CBPrimalMatrix(pm::CBMatrix)

copy constructor, *this=pm
"""
CBPrimalMatrix(pm::CBMatrix) = CBPrimalMatrix(@ccall libcb.cb_primalmatrix_new5(pm.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    Base.copy!(self::CBPrimalMatrix, pd::CBMatrix)

copy operator, *this=pm
"""
Base.copy!(self::CBPrimalMatrix, pd::CBMatrix) = (@ccall libcb.cb_primalmatrix_assign(self.data::Ptr{Cvoid}, pd.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

@doc raw"""
    cb_aggregate_primal_data!(self::CBPrimalMatrix, it::CBPrimalData, itsfactor::Real)

multiply *this Matrix with myfactor and add itsfactor*it (it must dynamic_cast to a PrimalMatrix)
"""
cb_aggregate_primal_data!(self::CBPrimalMatrix, it::CBPrimalData, itsfactor::Real) = @ccall libcb.cb_primalmatrix_aggregate_primal_data(self.data::Ptr{Cvoid}, it.data::Ptr{Cvoid}, itsfactor::Cdouble)::Cint

@doc raw"""
    cb_scale_primal_data!(self::CBPrimalMatrix, myfactor::Real)

multiply/scale *this with a nonnegative myfactor
"""
cb_scale_primal_data!(self::CBPrimalMatrix, myfactor::Real) = @ccall libcb.cb_primalmatrix_scale_primal_data(self.data::Ptr{Cvoid}, myfactor::Cdouble)::Cint

