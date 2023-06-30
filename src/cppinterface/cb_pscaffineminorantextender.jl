@doc raw"""
    CBPSCAffineMinorantExtender(amf::Union{<:CBPSCAffineFunction,Nothing})

the PSCAffineFunction pointed to has to be valid for the livetime of this object
"""
CBPSCAffineMinorantExtender(amf::Union{<:CBPSCAffineFunction,Nothing}) = CBPSCAffineMinorantExtender(@ccall libcb.cb_pscaffineminorantextender_new((isnothing(amf) ? C_NULL : amf.data)::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_extend!(self::CBPSCAffineMinorantExtender, minorant::CBMinorant, n_coords::Integer, indices::Union{<:AbstractVector{Integer},Nothing})

@brief called by ConicBundle to update internal Minorant objects, has to return 0 on success

        @param[in,out] minorant  (Minorant&)
            it holds a (possibly aggregated) minorant that was generated
            from minorants returned by oracle calls, e.g. as in
      FunctionOracle::evaluate() If PrimalData was provided in these
      minorants, this will be aggregated along and will also be
      available in this minorant.

        @param[in] n_coords (int)
            the number of coordinate positions that have to be filled in

        @param[out] new_subgradient_values  (DVector &)
      the indices of these coordinate positions (sorted in
      strictly increasing order)

  @return
           -  0 on success,
           -  1 if extension/update is impossible

     
"""
cb_extend!(self::CBPSCAffineMinorantExtender, minorant::CBMinorant, n_coords::Integer, indices::Union{<:AbstractVector{Integer},Nothing}) = GC.@preserve indices begin
    (LinearAlgebra.chkstride1(indices); @ccall libcb.cb_pscaffineminorantextender_extend(self.data::Ptr{Cvoid}, minorant.data::Ptr{Cvoid}, n_coords::Cint, indices::Ptr{Cint})::Cint)
end

