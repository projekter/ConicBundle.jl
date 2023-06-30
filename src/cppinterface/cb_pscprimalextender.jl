@doc raw"""
    cb_extend!(self::CBPSCPrimalExtender, param0::CBPrimalData)

like in PrimalExtender, called by ConicBundle to update internal PrimalData objects, has to return 0 on success
"""
cb_extend!(self::CBPSCPrimalExtender, param0::CBPrimalData) = @ccall libcb.cb_pscprimalextender_extend(self.data::Ptr{Cvoid}, param0.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_extend_Ritz!(self::CBPSCPrimalExtender, param0::CBMatrix)

called by ConicBundle to update internal Ritz_vectors, has to return 0 on success
"""
cb_extend_Ritz!(self::CBPSCPrimalExtender, param0::CBMatrix) = @ccall libcb.cb_pscprimalextender_extend_ritz(self.data::Ptr{Cvoid}, param0.data::Ptr{Cvoid})::Cint

