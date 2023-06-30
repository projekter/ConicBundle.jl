@doc raw"""
    cb_init!(self::CBGB_rand, seed::Integer = 1)

restart generator with seed
"""
cb_init!(self::CBGB_rand, seed::Integer = 1) = @ccall libcb.cb_gb_rand_init(self.data::Ptr{Cvoid}, seed::Clong)::Cvoid

@doc raw"""
    CBGB_rand(seed::Integer = 1)

calls  init(seed)
"""
CBGB_rand(seed::Integer = 1) = CBGB_rand(@ccall libcb.cb_gb_rand_new(seed::Clong)::Ptr{Cvoid})

@doc raw"""
    cb_unif_long!(self::CBGB_rand, m::Integer)

returns a random integer number "uniformly distributed" in {0,..,m-1}
"""
cb_unif_long!(self::CBGB_rand, m::Integer) = @ccall libcb.cb_gb_rand_unif_long(self.data::Ptr{Cvoid}, m::Clong)::Clong

@doc raw"""
    cb_next!(self::CBGB_rand)

returns a random double number "uniformly distributed" in (0,1)
"""
cb_next!(self::CBGB_rand) = @ccall libcb.cb_gb_rand_next(self.data::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_save(self::CBGB_rand)

save current configuration to out so as to continue identically after restore
"""
cb_save(self::CBGB_rand) = @ccall libcb.cb_gb_rand_save(self.data::Ptr{Cvoid})::Cvoid

