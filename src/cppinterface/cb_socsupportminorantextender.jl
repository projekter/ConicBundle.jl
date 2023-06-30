@doc raw"""
    CBSOCSupportMinorantExtender(fun::Union{<:CBSOCSupportFunction,Nothing})

the SOCSupportFunction pointed to has to be valid for the livetime of this object
"""
CBSOCSupportMinorantExtender(fun::Union{<:CBSOCSupportFunction,Nothing}) = CBSOCSupportMinorantExtender(@ccall libcb.cb_socsupportminorantextender_new((isnothing(fun) ? C_NULL : fun.data)::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_extend!(self::CBSOCSupportMinorantExtender, minorant::CBMinorant, n_coords::Integer, indices::Union{<:AbstractVector{Integer},Nothing})

@brief called by ConicBundle to update internal Minorant objects, has to return 0 on success

        for each relevant index i the minorant value is set to the projection
  of 0 onto the interval [lower_bound(i),upper_bound(i)]

        @param[in,out] minorant  (Minorant&)
            it holds a (possibly aggregated) minorant that was generated
            from minorants returned by oracle calls, e.g. as in
      FunctionOracle::evaluate()

        @param[in] n_coords (int)
            the number of coordinate positions that have to be filled in

        @param[out] new_subgradient_values  (DVector &)
      the indices of these coordinate positions (sorted in
      strictly increasing order)

  @return
           -  0 on success,
           -  1 if extension/update is impossible

     
"""
cb_extend!(self::CBSOCSupportMinorantExtender, minorant::CBMinorant, n_coords::Integer, indices::Union{<:AbstractVector{Integer},Nothing}) = GC.@preserve indices begin
    (LinearAlgebra.chkstride1(indices); @ccall libcb.cb_socsupportminorantextender_extend(self.data::Ptr{Cvoid}, minorant.data::Ptr{Cvoid}, n_coords::Cint, indices::Ptr{Cint})::Cint)
end

