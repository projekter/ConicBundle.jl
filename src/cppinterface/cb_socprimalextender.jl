@doc raw"""
    cb_extend!(self::CBSOCPrimalExtender, param0::CBPrimalData)

like in PrimalExtender, called by ConicBundle to update internal PrimalData objects, has to return 0 on success
"""
cb_extend!(self::CBSOCPrimalExtender, param0::CBPrimalData) = @ccall libcb.cb_socprimalextender_extend(self.data::Ptr{Cvoid}, param0.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_extend_SOC!(self::CBSOCPrimalExtender, param0::CBMatrix)

called by ConicBundle to update internal SOC vectors, has to return 0 on success
"""
cb_extend_SOC!(self::CBSOCPrimalExtender, param0::CBMatrix) = @ccall libcb.cb_socprimalextender_extend_soc(self.data::Ptr{Cvoid}, param0.data::Ptr{Cvoid})::Cint

