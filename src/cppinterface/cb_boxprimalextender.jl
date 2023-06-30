@doc raw"""
    cb_extend!(self::CBBoxPrimalExtender, param0::CBPrimalData)

like in PrimalExtender, called by ConicBundle to update internal PrimalData objects, has to return 0 on success
"""
cb_extend!(self::CBBoxPrimalExtender, param0::CBPrimalData) = @ccall libcb.cb_boxprimalextender_extend(self.data::Ptr{Cvoid}, param0.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_extend_Box!(self::CBBoxPrimalExtender, param0::CBMatrix)

called by ConicBundle to update internal Ritz_vectors, has to return 0 on success
"""
cb_extend_Box!(self::CBBoxPrimalExtender, param0::CBMatrix) = @ccall libcb.cb_boxprimalextender_extend_box(self.data::Ptr{Cvoid}, param0.data::Ptr{Cvoid})::Cint

