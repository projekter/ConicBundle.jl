@doc raw"""
    cb_set_defaults!(self::CBBundleTerminator)

sets the default parameter values
"""
cb_set_defaults!(self::CBBundleTerminator) = @ccall libcb.cb_bundleterminator_set_defaults(self.data::Ptr{Cvoid})::Cvoid

@doc raw"""
    cb_clear!(self::CBBundleTerminator)

resets @a terminated
"""
cb_clear!(self::CBBundleTerminator) = @ccall libcb.cb_bundleterminator_clear(self.data::Ptr{Cvoid})::Cvoid

@doc raw"""
    CBBundleTerminator(incr::Integer = -1)

calls set_defaults() and clear()
"""
CBBundleTerminator(incr::Integer = -1) = CBBundleTerminator(@ccall libcb.cb_bundleterminator_new(incr::Cint)::Ptr{Cvoid})

@doc raw"""
    cb_set_termeps!(self::CBBundleTerminator, teps::Real)

set the termination precision (>0!)
"""
cb_set_termeps!(self::CBBundleTerminator, teps::Real) = @ccall libcb.cb_bundleterminator_set_termeps(self.data::Ptr{Cvoid}, teps::Cdouble)::Cvoid

@doc raw"""
    cb_get_termeps(self::CBBundleTerminator)

returns the current termination precision
"""
cb_get_termeps(self::CBBundleTerminator) = @ccall libcb.cb_bundleterminator_get_termeps(self.data::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_set_aggr_dnormsqr!(self::CBBundleTerminator, sg::Real)

set an upper bound for the dual norm squared of the aggregate at termination (or <=0 if no such bound is desired)
"""
cb_set_aggr_dnormsqr!(self::CBBundleTerminator, sg::Real) = @ccall libcb.cb_bundleterminator_set_aggr_dnormsqr(self.data::Ptr{Cvoid}, sg::Cdouble)::Cvoid

@doc raw"""
    cb_get_aggr_dnormsqr(self::CBBundleTerminator)

returns the current bound for the dual norm squared of the aggregate
"""
cb_get_aggr_dnormsqr(self::CBBundleTerminator) = @ccall libcb.cb_bundleterminator_get_aggr_dnormsqr(self.data::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_set_timelimit!(self::CBBundleTerminator, cp::Union{<:CBClock,Nothing}, tl::CBMicroseconds)

set cp==0 for no timelimit, otherwise specify clock and microseconds
"""
cb_set_timelimit!(self::CBBundleTerminator, cp::Union{<:CBClock,Nothing}, tl::CBMicroseconds) = @ccall libcb.cb_bundleterminator_set_timelimit(self.data::Ptr{Cvoid}, (isnothing(cp) ? C_NULL : cp.data)::Ptr{Cvoid}, tl.data::Ptr{Cvoid})::Cvoid

@doc raw"""
    cb_get_timelimit(self::CBBundleTerminator)

returns the timelimit value
"""
cb_get_timelimit(self::CBBundleTerminator) = CBCH_Tools::Microseconds(@ccall libcb.cb_bundleterminator_new_get_timelimit(self.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_set_recomplimit!(self::CBBundleTerminator, rl::Integer)

set upper bound on the value returned by BundleTerminatorData::get_sumrecomp, <0 if no limit
"""
cb_set_recomplimit!(self::CBBundleTerminator, rl::Integer) = @ccall libcb.cb_bundleterminator_set_recomplimit(self.data::Ptr{Cvoid}, rl::Cint)::Cvoid

@doc raw"""
    cb_get_recomplimit(self::CBBundleTerminator)

returns the current value of this parameter
"""
cb_get_recomplimit(self::CBBundleTerminator) = @ccall libcb.cb_bundleterminator_get_recomplimit(self.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_set_qpfailslimit!(self::CBBundleTerminator, ql::Integer)

set upper bound on the value returned by BundleTerminatorData::get_sumqpfails, <0 if no limit
"""
cb_set_qpfailslimit!(self::CBBundleTerminator, ql::Integer) = @ccall libcb.cb_bundleterminator_set_qpfailslimit(self.data::Ptr{Cvoid}, ql::Cint)::Cvoid

@doc raw"""
    cb_get_qpfailslimit(self::CBBundleTerminator)

returns the current value of this parameter
"""
cb_get_qpfailslimit(self::CBBundleTerminator) = @ccall libcb.cb_bundleterminator_get_qpfailslimit(self.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_set_modelfailslimit!(self::CBBundleTerminator, ml::Integer)

set upper bound on the value returned by BundleTerminatorData::get_summodelfails, <0 if no limit
"""
cb_set_modelfailslimit!(self::CBBundleTerminator, ml::Integer) = @ccall libcb.cb_bundleterminator_set_modelfailslimit(self.data::Ptr{Cvoid}, ml::Cint)::Cvoid

@doc raw"""
    cb_get_modelfailslimit(self::CBBundleTerminator)

returns the current value of this parameter
"""
cb_get_modelfailslimit(self::CBBundleTerminator) = @ccall libcb.cb_bundleterminator_get_modelfailslimit(self.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_set_augvalfailslimit!(self::CBBundleTerminator, al::Integer)

set upper bound on the value returned by BundleTerminatorData::get_sumaugvalfails(), <0 if no limit
"""
cb_set_augvalfailslimit!(self::CBBundleTerminator, al::Integer) = @ccall libcb.cb_bundleterminator_set_augvalfailslimit(self.data::Ptr{Cvoid}, al::Cint)::Cvoid

@doc raw"""
    cb_get_augvalfailslimit(self::CBBundleTerminator)

returns the current value of this parameter
"""
cb_get_augvalfailslimit(self::CBBundleTerminator) = @ccall libcb.cb_bundleterminator_get_augvalfailslimit(self.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_set_objevallimit!(self::CBBundleTerminator, ol::Integer)

set upper bound on the value returned by BundleTerminatorData::get_cntobjeval, <0 if no limit
"""
cb_set_objevallimit!(self::CBBundleTerminator, ol::Integer) = @ccall libcb.cb_bundleterminator_set_objevallimit(self.data::Ptr{Cvoid}, ol::Cint)::Cvoid

@doc raw"""
    cb_get_objevallimit(self::CBBundleTerminator)

returns the current value of this parameter
"""
cb_get_objevallimit(self::CBBundleTerminator) = @ccall libcb.cb_bundleterminator_get_objevallimit(self.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_set_oraclefailslimit!(self::CBBundleTerminator, ol::Integer)

set upper bound on the value returned by BundleTerminatorData::get_sumoraclefails, <0 if no limit
"""
cb_set_oraclefailslimit!(self::CBBundleTerminator, ol::Integer) = @ccall libcb.cb_bundleterminator_set_oraclefailslimit(self.data::Ptr{Cvoid}, ol::Cint)::Cvoid

@doc raw"""
    cb_get_oraclefailslimit(self::CBBundleTerminator)

returns the current value of this parameter
"""
cb_get_oraclefailslimit(self::CBBundleTerminator) = @ccall libcb.cb_bundleterminator_get_oraclefailslimit(self.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_get_terminated(self::CBBundleTerminator)

return the termination code returned in the last call to check_termination()
"""
cb_get_terminated(self::CBBundleTerminator) = @ccall libcb.cb_bundleterminator_get_terminated(self.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_clear_terminated!(self::CBBundleTerminator)

reset the termination code to zero (this does not remove the reason for the termination; for this, set other bounds or clear the numbers supplied by BundleTerminationData)
"""
cb_clear_terminated!(self::CBBundleTerminator) = @ccall libcb.cb_bundleterminator_clear_terminated(self.data::Ptr{Cvoid})::Cvoid

@doc raw"""
    cb_print_status(self::CBBundleTerminator)

output an explanation string for the current termination code
"""
cb_print_status(self::CBBundleTerminator) = @ccall libcb.cb_bundleterminator_print_status(self.data::Ptr{Cvoid})::Cvoid

@doc raw"""
    cb_save(self::CBBundleTerminator)

output current parameter settings
"""
cb_save(self::CBBundleTerminator) = @ccall libcb.cb_bundleterminator_save(self.data::Ptr{Cvoid})::Cvoid

