@doc raw"""
    cb_init!(self::CBPSCBundleParameters, bp::CBBundleParameters)

*@brief initialize to given values 
"""
cb_init!(self::CBPSCBundleParameters, bp::CBBundleParameters) = @ccall libcb.cb_pscbundleparameters_init(self.data::Ptr{Cvoid}, bp.data::Ptr{Cvoid})::Cint

@doc raw"""
    CBPSCBundleParameters()

default constructor
"""
CBPSCBundleParameters() = CBPSCBundleParameters(@ccall libcb.cb_pscbundleparameters_new()::Ptr{Cvoid})

@doc raw"""
    CBPSCBundleParameters(bp::CBBundleParameters)

"copy" constructor
"""
CBPSCBundleParameters(bp::CBBundleParameters) = CBPSCBundleParameters(@ccall libcb.cb_pscbundleparameters_new2(bp.data::Ptr{Cvoid})::Ptr{Cvoid})

