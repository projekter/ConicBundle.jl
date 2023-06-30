@doc raw"""
    cb_clear!(self::CBUQPSumModelBlock)

reset to "empty/no" model
"""
cb_clear!(self::CBUQPSumModelBlock) = @ccall libcb.cb_uqpsummodelblock_clear(self.data::Ptr{Cvoid})::Cvoid

@doc raw"""
    CBUQPSumModelBlock(incr::Integer = -1)

default constructor
"""
CBUQPSumModelBlock(incr::Integer = -1) = CBUQPSumModelBlock(@ccall libcb.cb_uqpsummodelblock_new(incr::Cint)::Ptr{Cvoid})

@doc raw"""
    cb_append!(self::CBUQPSumModelBlock, inblock::Union{<:CBQPModelDataObject,Nothing})

add another model at the end of the list
"""
cb_append!(self::CBUQPSumModelBlock, inblock::Union{<:CBQPModelDataObject,Nothing}) = @ccall libcb.cb_uqpsummodelblock_append(self.data::Ptr{Cvoid}, (isnothing(inblock) ? C_NULL : inblock.data)::Ptr{Cvoid})::Cint

@doc raw"""
    cb_xdim(self::CBUQPSumModelBlock)

sum of all xdim of the sublocks
"""
cb_xdim(self::CBUQPSumModelBlock) = @ccall libcb.cb_uqpsummodelblock_xdim(self.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_ydim(self::CBUQPSumModelBlock)

sum of all ydim of the sublocks
"""
cb_ydim(self::CBUQPSumModelBlock) = @ccall libcb.cb_uqpsummodelblock_ydim(self.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_set_qp_xstart!(self::CBUQPSumModelBlock, x_start_index::Integer)

set the starting index of x also for all subblocks
"""
cb_set_qp_xstart!(self::CBUQPSumModelBlock, x_start_index::Integer) = @ccall libcb.cb_uqpsummodelblock_set_qp_xstart(self.data::Ptr{Cvoid}, x_start_index::Cint)::Cint

@doc raw"""
    cb_set_qp_ystart!(self::CBUQPSumModelBlock, y_start_index::Integer)

set the starting index of y also for all subblocks
"""
cb_set_qp_ystart!(self::CBUQPSumModelBlock, y_start_index::Integer) = @ccall libcb.cb_uqpsummodelblock_set_qp_ystart(self.data::Ptr{Cvoid}, y_start_index::Cint)::Cint

@doc raw"""
    cb_starting_x!(self::CBUQPSumModelBlock, qp_x::CBMatrix)

get the starting x of all subblocks
"""
cb_starting_x!(self::CBUQPSumModelBlock, qp_x::CBMatrix) = @ccall libcb.cb_uqpsummodelblock_starting_x(self.data::Ptr{Cvoid}, qp_x.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_starting_y!(self::CBUQPSumModelBlock, qp_y::CBMatrix, qp_Qx::CBMatrix, qp_c::CBMatrix)

get the starting y information of all subblocks
"""
cb_starting_y!(self::CBUQPSumModelBlock, qp_y::CBMatrix, qp_Qx::CBMatrix, qp_c::CBMatrix) = @ccall libcb.cb_uqpsummodelblock_starting_y(self.data::Ptr{Cvoid}, qp_y.data::Ptr{Cvoid}, qp_Qx.data::Ptr{Cvoid}, qp_c.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_get_local_primalcost(self::CBUQPSumModelBlock)

get joint primalcost of all subblocks
"""
cb_get_local_primalcost(self::CBUQPSumModelBlock) = @ccall libcb.cb_uqpsummodelblock_get_local_primalcost(self.data::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_get_local_dualcost(self::CBUQPSumModelBlock)

get joint dualcost of all subblocks
"""
cb_get_local_dualcost(self::CBUQPSumModelBlock) = @ccall libcb.cb_uqpsummodelblock_get_local_dualcost(self.data::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_get_Ab(self::CBUQPSumModelBlock, qp_A::CBMatrix, qp_b::CBMatrix)

get the A matrix of all subblocks and store it consistently
"""
cb_get_Ab(self::CBUQPSumModelBlock, qp_A::CBMatrix, qp_b::CBMatrix) = @ccall libcb.cb_uqpsummodelblock_get_ab(self.data::Ptr{Cvoid}, qp_A.data::Ptr{Cvoid}, qp_b.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_restart_x!(self::CBUQPSumModelBlock, qp_x::CBMatrix, qp_c::CBMatrix, qp_dc::CBMatrix)

get a good restarting x of all subblocks for this change in costs
"""
cb_restart_x!(self::CBUQPSumModelBlock, qp_x::CBMatrix, qp_c::CBMatrix, qp_dc::CBMatrix) = @ccall libcb.cb_uqpsummodelblock_restart_x(self.data::Ptr{Cvoid}, qp_x.data::Ptr{Cvoid}, qp_c.data::Ptr{Cvoid}, qp_dc.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_restart_y!(self::CBUQPSumModelBlock, qp_y::CBMatrix, qp_Qx::CBMatrix, qp_c::CBMatrix, qp_dc::CBMatrix)

get a good restarting y  of all subblocks for this change in costs
"""
cb_restart_y!(self::CBUQPSumModelBlock, qp_y::CBMatrix, qp_Qx::CBMatrix, qp_c::CBMatrix, qp_dc::CBMatrix) = @ccall libcb.cb_uqpsummodelblock_restart_y(self.data::Ptr{Cvoid}, qp_y.data::Ptr{Cvoid}, qp_Qx.data::Ptr{Cvoid}, qp_c.data::Ptr{Cvoid}, qp_dc.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_add_xinv_kron_z!(self::CBUQPSumModelBlock, barQ::CBSymmatrix)

add this for all subblocks
"""
cb_add_xinv_kron_z!(self::CBUQPSumModelBlock, barQ::CBSymmatrix) = @ccall libcb.cb_uqpsummodelblock_add_xinv_kron_z(self.data::Ptr{Cvoid}, barQ.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_add_local_sys!(self::CBUQPSumModelBlock, sysdy::CBSymmatrix, rhs::CBMatrix)

add this for all subblocks
"""
cb_add_local_sys!(self::CBUQPSumModelBlock, sysdy::CBSymmatrix, rhs::CBMatrix) = @ccall libcb.cb_uqpsummodelblock_add_local_sys(self.data::Ptr{Cvoid}, sysdy.data::Ptr{Cvoid}, rhs.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_suggest_mu!(self::CBUQPSumModelBlock, qp_dx::CBMatrix, qp_dy::CBMatrix, rhs_residual::CBMatrix)

get this from all subblocks
"""
function cb_suggest_mu!(self::CBUQPSumModelBlock, qp_dx::CBMatrix, qp_dy::CBMatrix, rhs_residual::CBMatrix)
    sigma = Ref{Float64}()
    mu_dim = Ref{Int}()
    ip_xz = Ref{Float64}()
    @ccall libcb.cb_uqpsummodelblock_suggest_mu(self.data::Ptr{Cvoid}, ip_xz::Ref{Float64}, mu_dim::Ref{Int}, sigma::Ref{Float64}, qp_dx.data::Ptr{Cvoid}, qp_dy.data::Ptr{Cvoid}, rhs_residual.data::Ptr{Cvoid})::Cint
    return ip_xz[], mu_dim[], sigma[]
end

@doc raw"""
    cb_get_corr!(self::CBUQPSumModelBlock, xcorr::CBMatrix, rhs::CBMatrix, mu::Real)

get this from all subblocks
"""
cb_get_corr!(self::CBUQPSumModelBlock, xcorr::CBMatrix, rhs::CBMatrix, mu::Real) = @ccall libcb.cb_uqpsummodelblock_get_corr(self.data::Ptr{Cvoid}, xcorr.data::Ptr{Cvoid}, rhs.data::Ptr{Cvoid}, mu::Cdouble)::Cint

@doc raw"""
    cb_line_search!(self::CBUQPSumModelBlock, qp_dx::CBMatrix, qp_dy::CBMatrix, rhs_residual::CBMatrix)

get/do this from/for all subblocks
"""
function cb_line_search!(self::CBUQPSumModelBlock, qp_dx::CBMatrix, qp_dy::CBMatrix, rhs_residual::CBMatrix)
    alpha = Ref{Float64}()
    @ccall libcb.cb_uqpsummodelblock_line_search(self.data::Ptr{Cvoid}, alpha::Ref{Float64}, qp_dx.data::Ptr{Cvoid}, qp_dy.data::Ptr{Cvoid}, rhs_residual.data::Ptr{Cvoid})::Cint
    return alpha[]
end

@doc raw"""
    cb_set_point!(self::CBUQPSumModelBlock, qp_x::CBMatrix, qp_y::CBMatrix, alpha::Real)

do this for all subblocks
"""
cb_set_point!(self::CBUQPSumModelBlock, qp_x::CBMatrix, qp_y::CBMatrix, alpha::Real) = @ccall libcb.cb_uqpsummodelblock_set_point(self.data::Ptr{Cvoid}, qp_x.data::Ptr{Cvoid}, qp_y.data::Ptr{Cvoid}, alpha::Cdouble)::Cint

@doc raw"""
    cb_add_modelx_aggregate!(self::CBUQPSumModelBlock, gradient::CBMatrix)

do this for all subblocks
"""
function cb_add_modelx_aggregate!(self::CBUQPSumModelBlock, gradient::CBMatrix)
    offset = Ref{Float64}()
    @ccall libcb.cb_uqpsummodelblock_add_modelx_aggregate(self.data::Ptr{Cvoid}, offset::Ref{Float64}, gradient.data::Ptr{Cvoid})::Cint
    return offset[]
end

@doc raw"""
    cb_set_out!(self::CBUQPSumModelBlock, pril::Integer = 1)

do this for all subblocks
"""
cb_set_out!(self::CBUQPSumModelBlock, pril::Integer = 1) = @ccall libcb.cb_uqpsummodelblock_set_out(self.data::Ptr{Cvoid}, pril::Cint)::Cvoid

@doc raw"""
    cb_set_cbout!(self::CBUQPSumModelBlock, incr::Integer = -1)

do this for all subblocks
"""
cb_set_cbout!(self::CBUQPSumModelBlock, incr::Integer = -1) = @ccall libcb.cb_uqpsummodelblock_set_cbout(self.data::Ptr{Cvoid}, incr::Cint)::Cvoid

@doc raw"""
    cb_add_Bs(self::CBUQPSumModelBlock, qp_vec::CBMatrix)

do this for all subblocks
"""
cb_add_Bs(self::CBUQPSumModelBlock, qp_vec::CBMatrix) = (@ccall libcb.cb_uqpsummodelblock_add_bs(self.data::Ptr{Cvoid}, qp_vec.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

@doc raw"""
    cb_subtract_z(self::CBUQPSumModelBlock, dual_residual::CBMatrix, with_step::Bool = false)

do this for all subblocks
"""
cb_subtract_z(self::CBUQPSumModelBlock, dual_residual::CBMatrix, with_step::Bool = false) = (@ccall libcb.cb_uqpsummodelblock_subtract_z(self.data::Ptr{Cvoid}, dual_residual.data::Ptr{Cvoid}, with_step::Cint)::Ptr{Cvoid}; return self)

