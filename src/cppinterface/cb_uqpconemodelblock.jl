@doc raw"""
    CBUQPConeModelBlock(cbinc::Integer = -1)

default constructor
"""
CBUQPConeModelBlock(cbinc::Integer = -1) = CBUQPConeModelBlock(@ccall libcb.cb_uqpconemodelblock_new(cbinc::Cint)::Ptr{Cvoid})

@doc raw"""
    cb_clear!(self::CBUQPConeModelBlock)

reinitialize to "empty/no" model
"""
cb_clear!(self::CBUQPConeModelBlock) = @ccall libcb.cb_uqpconemodelblock_clear(self.data::Ptr{Cvoid})::Cvoid

@doc raw"""
    cb_init!(self::CBUQPConeModelBlock, constant_minorant::CBMinorantPointer, bundle::CBMinorantBundle, nnc_dim::Integer, soc_dim::CBIndexmatrix, sdp_dim::CBIndexmatrix, box_lb::CBMatrix, box_ub::CBMatrix, b::Real, ft::CBFunctionTask, oracle_data::Union{<:CBQPModelOracleDataObject,Nothing} = nothing, scale_box::Bool = true)

sets up the model with bundle information and how to combine it, see QPConeModelDataObject::init() for a detailed description
"""
cb_init!(self::CBUQPConeModelBlock, constant_minorant::CBMinorantPointer, bundle::CBMinorantBundle, nnc_dim::Integer, soc_dim::CBIndexmatrix, sdp_dim::CBIndexmatrix, box_lb::CBMatrix, box_ub::CBMatrix, b::Real, ft::CBFunctionTask, oracle_data::Union{<:CBQPModelOracleDataObject,Nothing} = nothing, scale_box::Bool = true) = @ccall libcb.cb_uqpconemodelblock_init(self.data::Ptr{Cvoid}, constant_minorant.data::Ptr{Cvoid}, bundle.data::Ptr{Cvoid}, nnc_dim::Cint, soc_dim.data::Ptr{Cvoid}, sdp_dim.data::Ptr{Cvoid}, box_lb.data::Ptr{Cvoid}, box_ub.data::Ptr{Cvoid}, b::Cdouble, ft::Cint, (isnothing(oracle_data) ? C_NULL : oracle_data.data)::Ptr{Cvoid}, scale_box::Cint)::Cint

@doc raw"""
    cb_adjust_trace!(self::CBUQPConeModelBlock, b::Real)

change the right hand side of the trace constraint to b
"""
cb_adjust_trace!(self::CBUQPConeModelBlock, b::Real) = @ccall libcb.cb_uqpconemodelblock_adjust_trace(self.data::Ptr{Cvoid}, b::Cdouble)::Cint

@doc raw"""
    cb_evaluate_trace(self::CBUQPConeModelBlock)

evaluate the left hand side of the trace constraint for modelx
"""
cb_evaluate_trace(self::CBUQPConeModelBlock) = @ccall libcb.cb_uqpconemodelblock_evaluate_trace(self.data::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_get_trace!(self::CBUQPConeModelBlock)

get the right hand side of the trace constraint
"""
cb_get_trace!(self::CBUQPConeModelBlock) = @ccall libcb.cb_uqpconemodelblock_get_trace(self.data::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_get_nncx!(self::CBUQPConeModelBlock, nncx::CBMatrix, nncx_activity::Union{<:CBMatrix,Nothing} = nothing, cautious::Bool = false)

get the linear part of modelx (and a guess, which of them are active, in {0.,1.})
"""
cb_get_nncx!(self::CBUQPConeModelBlock, nncx::CBMatrix, nncx_activity::Union{<:CBMatrix,Nothing} = nothing, cautious::Bool = false) = @ccall libcb.cb_uqpconemodelblock_get_nncx(self.data::Ptr{Cvoid}, nncx.data::Ptr{Cvoid}, (isnothing(nncx_activity) ? C_NULL : nncx_activity.data)::Ptr{Cvoid}, cautious::Cint)::Cint

@doc raw"""
    cb_get_socx!(self::CBUQPConeModelBlock, i::Integer, socx::CBMatrix, socx_activity::Union{<:AbstractVector{Cdouble},Nothing}, cautious::Bool = false)

get the SOC part of modelx (and a guess whether the entire cone is active
"""
cb_get_socx!(self::CBUQPConeModelBlock, i::Integer, socx::CBMatrix, socx_activity::Union{<:AbstractVector{Cdouble},Nothing}, cautious::Bool = false) = GC.@preserve socx_activity begin
    (LinearAlgebra.chkstride1(socx_activity); @ccall libcb.cb_uqpconemodelblock_get_socx(self.data::Ptr{Cvoid}, i::Cint, socx.data::Ptr{Cvoid}, socx_activity::Ptr{Cdouble}, cautious::Cint)::Cint)
end

@doc raw"""
    cb_get_pscx!(self::CBUQPConeModelBlock, i::Integer, pscx_eigs::CBMatrix, pscx_vecs::CBMatrix, pscx_primalgrowth::CBMatrix, pscx_dualgrowth::CBMatrix)

get the PSC part of modelx (and a guess on the rank of the active part)
"""
function cb_get_pscx!(self::CBUQPConeModelBlock, i::Integer, pscx_eigs::CBMatrix, pscx_vecs::CBMatrix, pscx_primalgrowth::CBMatrix, pscx_dualgrowth::CBMatrix)
    pscx_growthrate = Ref{Float64}()
    @ccall libcb.cb_uqpconemodelblock_get_pscx(self.data::Ptr{Cvoid}, i::Cint, pscx_eigs.data::Ptr{Cvoid}, pscx_vecs.data::Ptr{Cvoid}, pscx_growthrate::Ref{Float64}, pscx_primalgrowth.data::Ptr{Cvoid}, pscx_dualgrowth.data::Ptr{Cvoid})::Cint
    return pscx_growthrate[]
end

@doc raw"""
    cb_get_boxx!(self::CBUQPConeModelBlock, param0::CBMatrix, param1::Union{<:CBMatrix,Nothing} = nothing, cautious::Bool = false)

get the box part of modelx (and a guess, which of the bounds are active, in {0.,1.})
"""
cb_get_boxx!(self::CBUQPConeModelBlock, param0::CBMatrix, param1::Union{<:CBMatrix,Nothing} = nothing, cautious::Bool = false) = @ccall libcb.cb_uqpconemodelblock_get_boxx(self.data::Ptr{Cvoid}, param0.data::Ptr{Cvoid}, (isnothing(param1) ? C_NULL : param1.data)::Ptr{Cvoid}, cautious::Cint)::Cint

@doc raw"""
    cb_add_modelx_aggregate!(self::CBUQPConeModelBlock, gradient::CBMatrix)

adds opB transposed times modelx (with offsets but without constant affine term) to the arguments
"""
function cb_add_modelx_aggregate!(self::CBUQPConeModelBlock, gradient::CBMatrix)
    offset = Ref{Float64}()
    @ccall libcb.cb_uqpconemodelblock_add_modelx_aggregate(self.data::Ptr{Cvoid}, offset::Ref{Float64}, gradient.data::Ptr{Cvoid})::Cint
    return offset[]
end

@doc raw"""
    cb_tracedual(self::CBUQPConeModelBlock, prec::Union{<:AbstractVector{Cdouble},Nothing} = nothing)

return the value of the dual variable to the trace consrat == support function value
"""
cb_tracedual(self::CBUQPConeModelBlock, prec::Union{<:AbstractVector{Cdouble},Nothing} = nothing) = GC.@preserve prec begin
    (LinearAlgebra.chkstride1(prec); @ccall libcb.cb_uqpconemodelblock_tracedual(self.data::Ptr{Cvoid}, prec::Ptr{Cdouble})::Cdouble)
end

@doc raw"""
    cb_get_nncz!(self::CBUQPConeModelBlock, vecz::CBMatrix)

get the current dual non negative cone point (of the solution)
"""
cb_get_nncz!(self::CBUQPConeModelBlock, vecz::CBMatrix) = @ccall libcb.cb_uqpconemodelblock_get_nncz(self.data::Ptr{Cvoid}, vecz.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_get_socz!(self::CBUQPConeModelBlock, i::Integer, vecz::CBMatrix)

get the current dual second order cone point to cone i (of the solution)
"""
cb_get_socz!(self::CBUQPConeModelBlock, i::Integer, vecz::CBMatrix) = @ccall libcb.cb_uqpconemodelblock_get_socz(self.data::Ptr{Cvoid}, i::Cint, vecz.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_get_X!(self::CBUQPConeModelBlock, i::Integer, X::CBSymmatrix)

get the current primal positive semidefinite cone point to cone i (of the solution)
"""
cb_get_X!(self::CBUQPConeModelBlock, i::Integer, X::CBSymmatrix) = @ccall libcb.cb_uqpconemodelblock_get_x(self.data::Ptr{Cvoid}, i::Cint, X.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_get_Z!(self::CBUQPConeModelBlock, i::Integer, Z::CBSymmatrix)

get the current dual positive semidefinite cone point to cone i (of the solution)
"""
cb_get_Z!(self::CBUQPConeModelBlock, i::Integer, Z::CBSymmatrix) = @ccall libcb.cb_uqpconemodelblock_get_z(self.data::Ptr{Cvoid}, i::Cint, Z.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_get_y!(self::CBUQPConeModelBlock)

get the current dual value of the trace constraint
"""
cb_get_y!(self::CBUQPConeModelBlock) = @ccall libcb.cb_uqpconemodelblock_get_y(self.data::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_get_s!(self::CBUQPConeModelBlock)

get the current slack value of the trace constraint
"""
cb_get_s!(self::CBUQPConeModelBlock) = @ccall libcb.cb_uqpconemodelblock_get_s(self.data::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_get_mu!(self::CBUQPConeModelBlock)

get the current barrier parameter
"""
cb_get_mu!(self::CBUQPConeModelBlock) = @ccall libcb.cb_uqpconemodelblock_get_mu(self.data::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_get_old_nncx!(self::CBUQPConeModelBlock, vecx::CBMatrix)

get the previous primal non negative cone point (of the solution)
"""
cb_get_old_nncx!(self::CBUQPConeModelBlock, vecx::CBMatrix) = @ccall libcb.cb_uqpconemodelblock_get_old_nncx(self.data::Ptr{Cvoid}, vecx.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_get_old_nncz!(self::CBUQPConeModelBlock, vecz::CBMatrix)

get the previous dual non negative cone point (of the solution)
"""
cb_get_old_nncz!(self::CBUQPConeModelBlock, vecz::CBMatrix) = @ccall libcb.cb_uqpconemodelblock_get_old_nncz(self.data::Ptr{Cvoid}, vecz.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_get_old_socx!(self::CBUQPConeModelBlock, i::Integer, vecx::CBMatrix)

get the previous primal second order cone point to cone i (of the solution)
"""
cb_get_old_socx!(self::CBUQPConeModelBlock, i::Integer, vecx::CBMatrix) = @ccall libcb.cb_uqpconemodelblock_get_old_socx(self.data::Ptr{Cvoid}, i::Cint, vecx.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_get_old_socz!(self::CBUQPConeModelBlock, i::Integer, vecz::CBMatrix)

get the previous dual second order cone point to cone i (of the solution)
"""
cb_get_old_socz!(self::CBUQPConeModelBlock, i::Integer, vecz::CBMatrix) = @ccall libcb.cb_uqpconemodelblock_get_old_socz(self.data::Ptr{Cvoid}, i::Cint, vecz.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_get_old_X!(self::CBUQPConeModelBlock, i::Integer, X::CBSymmatrix)

get the previous primal positive semidefinite cone point to cone i (of the solution)
"""
cb_get_old_X!(self::CBUQPConeModelBlock, i::Integer, X::CBSymmatrix) = @ccall libcb.cb_uqpconemodelblock_get_old_x(self.data::Ptr{Cvoid}, i::Cint, X.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_get_old_Z!(self::CBUQPConeModelBlock, i::Integer, Z::CBSymmatrix)

get the previous dual positive semidefinite cone point to cone i (of the solution)
"""
cb_get_old_Z!(self::CBUQPConeModelBlock, i::Integer, Z::CBSymmatrix) = @ccall libcb.cb_uqpconemodelblock_get_old_z(self.data::Ptr{Cvoid}, i::Cint, Z.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_get_old_y!(self::CBUQPConeModelBlock)

get the previous dual value of the trace constraint
"""
cb_get_old_y!(self::CBUQPConeModelBlock) = @ccall libcb.cb_uqpconemodelblock_get_old_y(self.data::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_get_old_s!(self::CBUQPConeModelBlock)

get the previous slack value of the trace constraint
"""
cb_get_old_s!(self::CBUQPConeModelBlock) = @ccall libcb.cb_uqpconemodelblock_get_old_s(self.data::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_get_old_mu!(self::CBUQPConeModelBlock)

get the previous barrier parameter
"""
cb_get_old_mu!(self::CBUQPConeModelBlock) = @ccall libcb.cb_uqpconemodelblock_get_old_mu(self.data::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_get_last_alpha!(self::CBUQPConeModelBlock)

get the most recent step size
"""
cb_get_last_alpha!(self::CBUQPConeModelBlock) = @ccall libcb.cb_uqpconemodelblock_get_last_alpha(self.data::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_compute_local_directions!(self::CBUQPConeModelBlock, qp_dx::CBMatrix, qp_dy::CBMatrix, rhs_resid::CBMatrix, dz::CBMatrix, duz::CBMatrix)

given the steps of the global part, compute the step of the local part
"""
function cb_compute_local_directions!(self::CBUQPConeModelBlock, qp_dx::CBMatrix, qp_dy::CBMatrix, rhs_resid::CBMatrix, dz::CBMatrix, duz::CBMatrix)
    ds = Ref{Float64}()
    box_ds = Ref{Float64}()
    @ccall libcb.cb_uqpconemodelblock_compute_local_directions(self.data::Ptr{Cvoid}, qp_dx.data::Ptr{Cvoid}, qp_dy.data::Ptr{Cvoid}, rhs_resid.data::Ptr{Cvoid}, dz.data::Ptr{Cvoid}, duz.data::Ptr{Cvoid}, box_ds::Ref{Float64}, ds::Ref{Float64})::Cint
    return box_ds[], ds[]
end

@doc raw"""
    cb_inner_line_search!(self::CBUQPConeModelBlock, qp_dx::CBMatrix, qp_dy::CBMatrix, dz::CBMatrix, duz::CBMatrix, box_ds::Real, ds::Real)

perform a line search for the given direction and return a feasible step length
"""
function cb_inner_line_search!(self::CBUQPConeModelBlock, qp_dx::CBMatrix, qp_dy::CBMatrix, dz::CBMatrix, duz::CBMatrix, box_ds::Real, ds::Real)
    alpha = Ref{Float64}()
    @ccall libcb.cb_uqpconemodelblock_inner_line_search(self.data::Ptr{Cvoid}, alpha::Ref{Float64}, qp_dx.data::Ptr{Cvoid}, qp_dy.data::Ptr{Cvoid}, dz.data::Ptr{Cvoid}, duz.data::Ptr{Cvoid}, box_ds::Cdouble, ds::Cdouble)::Cint
    return alpha[]
end

@doc raw"""
    cb_xdim(self::CBUQPConeModelBlock)

dimension of externally visible primal variables
"""
cb_xdim(self::CBUQPConeModelBlock) = @ccall libcb.cb_uqpconemodelblock_xdim(self.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_ydim(self::CBUQPConeModelBlock)

dimension of externally visible dual variables
"""
cb_ydim(self::CBUQPConeModelBlock) = @ccall libcb.cb_uqpconemodelblock_ydim(self.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_set_qp_xstart!(self::CBUQPConeModelBlock, x_start_index::Integer)

the indices of the local variables correspond to the indices of the qp variables x and z starting with this index; returns 0 on success, 1 on failure
"""
cb_set_qp_xstart!(self::CBUQPConeModelBlock, x_start_index::Integer) = @ccall libcb.cb_uqpconemodelblock_set_qp_xstart(self.data::Ptr{Cvoid}, x_start_index::Cint)::Cint

@doc raw"""
    cb_set_qp_ystart!(self::CBUQPConeModelBlock, y_start_index::Integer)

the indices of the local variables correspond to the indices of the qp variables y starting with this index; returns 0 on success, 1 on failure
"""
cb_set_qp_ystart!(self::CBUQPConeModelBlock, y_start_index::Integer) = @ccall libcb.cb_uqpconemodelblock_set_qp_ystart(self.data::Ptr{Cvoid}, y_start_index::Cint)::Cint

@doc raw"""
    cb_starting_x!(self::CBUQPConeModelBlock, qp_x::CBMatrix)

generate a strictly feasible primal starting point, store it in the qpx_range of x; returns 0 on success, 1 on failure
"""
cb_starting_x!(self::CBUQPConeModelBlock, qp_x::CBMatrix) = @ccall libcb.cb_uqpconemodelblock_starting_x(self.data::Ptr{Cvoid}, qp_x.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_starting_y!(self::CBUQPConeModelBlock, qp_y::CBMatrix, qp_Qx::CBMatrix, qp_c::CBMatrix)

generate a strictly feasible dual starting point, store it in the qpy_range of y,  x is fixed already by a previous call to starting_x and Qx=Q*x; returns 0 on success, 1 on failure
"""
cb_starting_y!(self::CBUQPConeModelBlock, qp_y::CBMatrix, qp_Qx::CBMatrix, qp_c::CBMatrix) = @ccall libcb.cb_uqpconemodelblock_starting_y(self.data::Ptr{Cvoid}, qp_y.data::Ptr{Cvoid}, qp_Qx.data::Ptr{Cvoid}, qp_c.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_get_Ab(self::CBUQPConeModelBlock, qp_A::CBMatrix, qp_b::CBMatrix)

store the local coefficients of matrices A and b in the positions corresponding to qpy_range (rows) and qpx_range (columns); returns 0 on success, 1 on failure
"""
cb_get_Ab(self::CBUQPConeModelBlock, qp_A::CBMatrix, qp_b::CBMatrix) = @ccall libcb.cb_uqpconemodelblock_get_ab(self.data::Ptr{Cvoid}, qp_A.data::Ptr{Cvoid}, qp_b.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_restart_x!(self::CBUQPConeModelBlock, qp_x::CBMatrix, qp_c::CBMatrix, qp_dc::CBMatrix)

* @brief it is assumed that the problem was solved already once and is now
     resolved for a new linear cost term qp_c that resulted from the old
     one by adding qp_dc.

     on input qp_x holds the old optimal solution and on output
     the coorespoind qpx_range should be replaced by a reasonable
     strictly feasible solution for x suitable for restarting
     (see also restart_yz)

     returns 0 on success, 1 on failure
   
"""
cb_restart_x!(self::CBUQPConeModelBlock, qp_x::CBMatrix, qp_c::CBMatrix, qp_dc::CBMatrix) = @ccall libcb.cb_uqpconemodelblock_restart_x(self.data::Ptr{Cvoid}, qp_x.data::Ptr{Cvoid}, qp_c.data::Ptr{Cvoid}, qp_dc.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_restart_y!(self::CBUQPConeModelBlock, qp_y::CBMatrix, qp_Qx::CBMatrix, qp_c::CBMatrix, qp_dc::CBMatrix)

* @brief this is called after restart_x (see there)

       on input qp_y and qp_z hold the old optimal solution and on output
       the coorespoind qpy/qpx_range should be replaced by a reasonable
       strictly feasible solution for y/z suitable for restarting

       returns 0 on success, 1 on failure
    
"""
cb_restart_y!(self::CBUQPConeModelBlock, qp_y::CBMatrix, qp_Qx::CBMatrix, qp_c::CBMatrix, qp_dc::CBMatrix) = @ccall libcb.cb_uqpconemodelblock_restart_y(self.data::Ptr{Cvoid}, qp_y.data::Ptr{Cvoid}, qp_Qx.data::Ptr{Cvoid}, qp_c.data::Ptr{Cvoid}, qp_dc.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_add_xinv_kron_z!(self::CBUQPConeModelBlock, barQ::CBSymmatrix)

add the system term corresponding to (xinv kron z) (that arises from solving the linearized perturbed complementarity system x*z =0 or =mu*I for dx in the preferred search direction) to the diagonal block corresponding to qpx_range x qpx_range
"""
cb_add_xinv_kron_z!(self::CBUQPConeModelBlock, barQ::CBSymmatrix) = @ccall libcb.cb_uqpconemodelblock_add_xinv_kron_z(self.data::Ptr{Cvoid}, barQ.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_add_local_sys!(self::CBUQPConeModelBlock, sysdy::CBSymmatrix, rhs::CBMatrix)

* @brief add the local system informatoin

      on input:
               sysdy= A*barQ^{-1}*A^T    (barQ as returned in add_xinv_kron_z)
               rhs= A*barQ^{-1}*(c-Q*x-A^T*y)-(b-A*x)

      if the block uses additional internal variables
      (like an additional term + B*s with s>=0 in the primal feasibility constr)
      then the corresponding block terms have now to be added to sysdy and rhs, eg,
         sysdy +=  B*(t^{-1} kron s)*B^T     (if t is the dual variable to s)
         rhs   +=  B*s - B*(t^{-1} kron s)*B^T*y
    
"""
cb_add_local_sys!(self::CBUQPConeModelBlock, sysdy::CBSymmatrix, rhs::CBMatrix) = @ccall libcb.cb_uqpconemodelblock_add_local_sys(self.data::Ptr{Cvoid}, sysdy.data::Ptr{Cvoid}, rhs.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_get_local_primalcost(self::CBUQPConeModelBlock)

returns the current local primal cost contribution <d,s>
"""
cb_get_local_primalcost(self::CBUQPConeModelBlock) = @ccall libcb.cb_uqpconemodelblock_get_local_primalcost(self.data::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_get_local_dualcost(self::CBUQPConeModelBlock)

returns the current local dual cost contribution
"""
cb_get_local_dualcost(self::CBUQPConeModelBlock) = @ccall libcb.cb_uqpconemodelblock_get_local_dualcost(self.data::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_suggest_mu!(self::CBUQPConeModelBlock, qp_dx::CBMatrix, qp_dy::CBMatrix, rhs_residual::CBMatrix)

* @brief supply the information for the choice of the next barrier parameter value

       dx, dy is the predictor direction giving rise to the
       rhs_residual -(c-At(y+dy)-Q(x+dx)). Compute the direction dz and
       local step and based on the predictor (x+dx,y+dy,z+dz) suggest a
       value for mu by specifying the inner product of the dual cone
       variables ip_xz=ip(x,z)+ip(s,t), the dimension of the conic
       variable space mu_dim= cone_x.dim+cone_s.dim a value for the
       factor on mu to obtain the new target
   
"""
function cb_suggest_mu!(self::CBUQPConeModelBlock, qp_dx::CBMatrix, qp_dy::CBMatrix, rhs_residual::CBMatrix)
    sigma = Ref{Float64}()
    mu_dim = Ref{Int}()
    ip_xz = Ref{Float64}()
    @ccall libcb.cb_uqpconemodelblock_suggest_mu(self.data::Ptr{Cvoid}, ip_xz::Ref{Float64}, mu_dim::Ref{Int}, sigma::Ref{Float64}, qp_dx.data::Ptr{Cvoid}, qp_dy.data::Ptr{Cvoid}, rhs_residual.data::Ptr{Cvoid})::Cint
    return ip_xz[], mu_dim[], sigma[]
end

@doc raw"""
    cb_get_corr!(self::CBUQPConeModelBlock, xcorr::CBMatrix, rhs::CBMatrix, mu::Real)

* @brief supply the information for the corrector

     on input (w.r.t. corresponding positions)
          xcorr = 0
          rhs as on output of add_local_sys

     on output the corresponding positions of xcorr should hold the corrector
     term of the search direction, eg,  xcorr = mu*x^{-1} - x^{-1}*dx*dz,
     and if the block holds additional local variables as in add_local_sys then

             rhs += B*(mu * t^{-1}- t^{-1}*dt*ds)

     has to be called after suggest_mu which computes the other directions
   
"""
cb_get_corr!(self::CBUQPConeModelBlock, xcorr::CBMatrix, rhs::CBMatrix, mu::Real) = @ccall libcb.cb_uqpconemodelblock_get_corr(self.data::Ptr{Cvoid}, xcorr.data::Ptr{Cvoid}, rhs.data::Ptr{Cvoid}, mu::Cdouble)::Cint

@doc raw"""
    cb_line_search!(self::CBUQPConeModelBlock, qp_dx::CBMatrix, qp_dy::CBMatrix, rhs_residual::CBMatrix)

* @brief perform a line search for the block variables

      dx,dy give the final step direction, alpha is on input
      an upper bound on the step size.

      On output alpha has to be no larger than on input and
      has to guarantee strict feasibility of the primal/dual step on
      the local variables.

      The block has to compute the step direction dz as well as
      for additional internal variables now and to choose alpha so
      that strict feasibility is guaranteed for the internal
      variables as well
   
"""
function cb_line_search!(self::CBUQPConeModelBlock, qp_dx::CBMatrix, qp_dy::CBMatrix, rhs_residual::CBMatrix)
    alpha = Ref{Float64}()
    @ccall libcb.cb_uqpconemodelblock_line_search(self.data::Ptr{Cvoid}, alpha::Ref{Float64}, qp_dx.data::Ptr{Cvoid}, qp_dy.data::Ptr{Cvoid}, rhs_residual.data::Ptr{Cvoid})::Cint
    return alpha[]
end

@doc raw"""
    cb_set_point!(self::CBUQPConeModelBlock, qp_x::CBMatrix, qp_y::CBMatrix, alpha::Real)

x,y,z is the new point and has to be stored, alpha is the step size used in the step, it is passed so thatthe block can take the same step for internal variables if needed.
"""
cb_set_point!(self::CBUQPConeModelBlock, qp_x::CBMatrix, qp_y::CBMatrix, alpha::Real) = @ccall libcb.cb_uqpconemodelblock_set_point(self.data::Ptr{Cvoid}, qp_x.data::Ptr{Cvoid}, qp_y.data::Ptr{Cvoid}, alpha::Cdouble)::Cint

@doc raw"""
    cb_add_Bs(self::CBUQPConeModelBlock, qp_vec::CBMatrix)

add the local product of matrices B and s in the positions corresponding to qpy_range (rows) and return qp_vec; returns 0 on success, 1 on failure
"""
cb_add_Bs(self::CBUQPConeModelBlock, qp_vec::CBMatrix) = (@ccall libcb.cb_uqpconemodelblock_add_bs(self.data::Ptr{Cvoid}, qp_vec.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

@doc raw"""
    cb_subtract_z(self::CBUQPConeModelBlock, dual_residual::CBMatrix, with_step::Bool = false)

add the contributions of the dual slacks and return dual_residual returns 0 on success, 1 on failure
"""
cb_subtract_z(self::CBUQPConeModelBlock, dual_residual::CBMatrix, with_step::Bool = false) = (@ccall libcb.cb_uqpconemodelblock_subtract_z(self.data::Ptr{Cvoid}, dual_residual.data::Ptr{Cvoid}, with_step::Cint)::Ptr{Cvoid}; return self)

