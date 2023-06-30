@doc raw"""
    cb_clear!(self::CBQPSumModelBlock)

reset to "empty/no" model
"""
cb_clear!(self::CBQPSumModelBlock) = @ccall libcb.cb_qpsummodelblock_clear(self.data::Ptr{Cvoid})::Cvoid

@doc raw"""
    CBQPSumModelBlock(cbincr::Integer = -1)

default constructor
"""
CBQPSumModelBlock(cbincr::Integer = -1) = CBQPSumModelBlock(@ccall libcb.cb_qpsummodelblock_new(cbincr::Cint)::Ptr{Cvoid})

@doc raw"""
    cb_clone!(self::CBQPSumModelBlock)

return a cloned object on the heap
"""
cb_clone!(self::CBQPSumModelBlock) = CBQPModelBlockObject(@ccall libcb.cb_qpsummodelblock_clone(self.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_recursive_delete_and_clear!(self::CBQPSumModelBlock)

usually the objects of the recursive block structure and not deleted in a clear. If needed, this can be invoked explicitly here, e.g., in order to clean up clones
"""
cb_recursive_delete_and_clear!(self::CBQPSumModelBlock) = @ccall libcb.cb_qpsummodelblock_recursive_delete_and_clear(self.data::Ptr{Cvoid})::Cvoid

@doc raw"""
    cb_recursive_copy_data_of!(self::CBQPSumModelBlock, param0::Union{<:CBQPModelBlockObject,Nothing})

sofar this is only needed for some comparative evaluations
"""
cb_recursive_copy_data_of!(self::CBQPSumModelBlock, param0::Union{<:CBQPModelBlockObject,Nothing}) = @ccall libcb.cb_qpsummodelblock_recursive_copy_data_of(self.data::Ptr{Cvoid}, (isnothing(param0) ? C_NULL : param0.data)::Ptr{Cvoid})::Cint

@doc raw"""
    cb_append!(self::CBQPSumModelBlock, inblock::Union{<:CBQPModelDataObject,Nothing})

add another model at the end of the list
"""
cb_append!(self::CBQPSumModelBlock, inblock::Union{<:CBQPModelDataObject,Nothing}) = @ccall libcb.cb_qpsummodelblock_append(self.data::Ptr{Cvoid}, (isnothing(inblock) ? C_NULL : inblock.data)::Ptr{Cvoid})::Cint

@doc raw"""
    cb_B_times!(self::CBQPSumModelBlock, A::CBMatrix, C::CBMatrix, alpha::Real, beta::Real, Btrans::Integer, Atrans::Integer, startindex_model::Integer, globalbundle::CBMinorantBundle, startindex_bundle::Integer)

C=beta*C+alpha*B*A where B and A may be transposed; carry out the model part of this beginning at startindex_model and beta for the part, that is added to (the calling routine has to make sure beta is not executed repeatedly if the same part is affected by other models as well)
"""
cb_B_times!(self::CBQPSumModelBlock, A::CBMatrix, C::CBMatrix, alpha::Real, beta::Real, Btrans::Integer, Atrans::Integer, startindex_model::Integer, globalbundle::CBMinorantBundle, startindex_bundle::Integer) = (@ccall libcb.cb_qpsummodelblock_b_times(self.data::Ptr{Cvoid}, A.data::Ptr{Cvoid}, C.data::Ptr{Cvoid}, alpha::Cdouble, beta::Cdouble, Btrans::Cint, Atrans::Cint, startindex_model::Cint, globalbundle.data::Ptr{Cvoid}, startindex_bundle::Cint)::Ptr{Cvoid}; return self)

@doc raw"""
    cb_times_B!(self::CBQPSumModelBlock, A::CBMatrix, C::CBMatrix, alpha::Real, beta::Real, Atrans::Integer, Btrans::Integer, startindex_model::Integer, globalbundle::CBMinorantBundle, startindex_bundle::Integer)

C=beta*C+alpha*A*B where A and B may be transposed; carry out the model part of this beginning at startindex_model
"""
cb_times_B!(self::CBQPSumModelBlock, A::CBMatrix, C::CBMatrix, alpha::Real, beta::Real, Atrans::Integer, Btrans::Integer, startindex_model::Integer, globalbundle::CBMinorantBundle, startindex_bundle::Integer) = (@ccall libcb.cb_qpsummodelblock_times_b(self.data::Ptr{Cvoid}, A.data::Ptr{Cvoid}, C.data::Ptr{Cvoid}, alpha::Cdouble, beta::Cdouble, Atrans::Cint, Btrans::Cint, startindex_model::Cint, globalbundle.data::Ptr{Cvoid}, startindex_bundle::Cint)::Ptr{Cvoid}; return self)

@doc raw"""
    cb_add_BDBt!(self::CBQPSumModelBlock, diagvec::CBMatrix, bigS::CBSymmatrix, minus::Bool, startindex::Integer, Bt::CBMatrix, startindex_model::Integer, globalbundle::CBMinorantBundle, startindex_bundle::Integer)

add the main diagonal block tranpose(projection)*diagvec*projection to bigS starting at startindex
"""
cb_add_BDBt!(self::CBQPSumModelBlock, diagvec::CBMatrix, bigS::CBSymmatrix, minus::Bool, startindex::Integer, Bt::CBMatrix, startindex_model::Integer, globalbundle::CBMinorantBundle, startindex_bundle::Integer) = (@ccall libcb.cb_qpsummodelblock_add_bdbt(self.data::Ptr{Cvoid}, diagvec.data::Ptr{Cvoid}, bigS.data::Ptr{Cvoid}, minus::Cint, startindex::Cint, Bt.data::Ptr{Cvoid}, startindex_model::Cint, globalbundle.data::Ptr{Cvoid}, startindex_bundle::Cint)::Ptr{Cvoid}; return self)

@doc raw"""
    cb_get_Bt!(self::CBQPSumModelBlock, Bt::CBMatrix, startindex_model::Integer, global_bundle::CBMinorantBundle, startindex_bundle::Integer)

get the current matrix for the coupling matrix Bt in the first row of blocks
"""
cb_get_Bt!(self::CBQPSumModelBlock, Bt::CBMatrix, startindex_model::Integer, global_bundle::CBMinorantBundle, startindex_bundle::Integer) = (@ccall libcb.cb_qpsummodelblock_get_bt(self.data::Ptr{Cvoid}, Bt.data::Ptr{Cvoid}, startindex_model::Cint, global_bundle.data::Ptr{Cvoid}, startindex_bundle::Cint)::Ptr{Cvoid}; return self)

@doc raw"""
    cb_get_modelx!(self::CBQPSumModelBlock, modelx::CBMatrix, startindex_model::Integer)

set the local modelx value in modelx beginning with startindex (initialize it, do not add)
"""
cb_get_modelx!(self::CBQPSumModelBlock, modelx::CBMatrix, startindex_model::Integer) = @ccall libcb.cb_qpsummodelblock_get_modelx(self.data::Ptr{Cvoid}, modelx.data::Ptr{Cvoid}, startindex_model::Cint)::Cint

@doc raw"""
    cb_get_modeldx!(self::CBQPSumModelBlock, modeldx::CBMatrix, startindex_model::Integer)

set the local modeldx value in modeldx beginning with startindex (initialize it, do not add)
"""
cb_get_modeldx!(self::CBQPSumModelBlock, modeldx::CBMatrix, startindex_model::Integer) = @ccall libcb.cb_qpsummodelblock_get_modeldx(self.data::Ptr{Cvoid}, modeldx.data::Ptr{Cvoid}, startindex_model::Cint)::Cint

@doc raw"""
    cb_get_modeldcstr!(self::CBQPSumModelBlock, modeldcstr::CBMatrix, startindex_constraints::Integer)

set the local modeldcstr value in modeldcstr beginning with startindex (initialize it, do not add)
"""
cb_get_modeldcstr!(self::CBQPSumModelBlock, modeldcstr::CBMatrix, startindex_constraints::Integer) = @ccall libcb.cb_qpsummodelblock_get_modeldcstr(self.data::Ptr{Cvoid}, modeldcstr.data::Ptr{Cvoid}, startindex_constraints::Cint)::Cint

@doc raw"""
    cb_add_modelx_aggregate!(self::CBQPSumModelBlock, vec::CBMatrix, global_bundle::CBMinorantBundle, startindex_bundle::Integer)

adds opB transposed times modelx (without constant affine term) to the arguments
"""
function cb_add_modelx_aggregate!(self::CBQPSumModelBlock, vec::CBMatrix, global_bundle::CBMinorantBundle, startindex_bundle::Integer)
    val = Ref{Float64}()
    @ccall libcb.cb_qpsummodelblock_add_modelx_aggregate(self.data::Ptr{Cvoid}, val::Ref{Float64}, vec.data::Ptr{Cvoid}, global_bundle.data::Ptr{Cvoid}, startindex_bundle::Cint)::Cint
    return val[]
end

@doc raw"""
    cb_get_sysviol_model!(self::CBQPSumModelBlock, modelvec::CBMatrix, startindex_model::Integer, dy::CBMatrix, global_bundle::CBMinorantBundle, startindex_bundle::Integer)

set the model violation for the current system solution
"""
cb_get_sysviol_model!(self::CBQPSumModelBlock, modelvec::CBMatrix, startindex_model::Integer, dy::CBMatrix, global_bundle::CBMinorantBundle, startindex_bundle::Integer) = @ccall libcb.cb_qpsummodelblock_get_sysviol_model(self.data::Ptr{Cvoid}, modelvec.data::Ptr{Cvoid}, startindex_model::Cint, dy.data::Ptr{Cvoid}, global_bundle.data::Ptr{Cvoid}, startindex_bundle::Cint)::Cint

@doc raw"""
    cb_get_sysviol_constraints!(self::CBQPSumModelBlock, constrvec::CBMatrix, startindex_constr::Integer)

set the constraint violation for the current system solution starting at this index
"""
cb_get_sysviol_constraints!(self::CBQPSumModelBlock, constrvec::CBMatrix, startindex_constr::Integer) = @ccall libcb.cb_qpsummodelblock_get_sysviol_constraints(self.data::Ptr{Cvoid}, constrvec.data::Ptr{Cvoid}, startindex_constr::Cint)::Cint

@doc raw"""
    cb_display_model_values!(self::CBQPSumModelBlock, y::CBMatrix, global_bundle::CBMinorantBundle, startindex_bundle::Integer)

output for debbuging purposes
"""
cb_display_model_values!(self::CBQPSumModelBlock, y::CBMatrix, global_bundle::CBMinorantBundle, startindex_bundle::Integer) = @ccall libcb.cb_qpsummodelblock_display_model_values(self.data::Ptr{Cvoid}, y.data::Ptr{Cvoid}, global_bundle.data::Ptr{Cvoid}, startindex_bundle::Cint)::Cvoid

@doc raw"""
    cb_reset_starting_point!(self::CBQPSumModelBlock, y::CBMatrix, mu::Real, global_bundle::CBMinorantBundle, startindex_bundle::Integer)

reset the starting point for this value of the design variables y
"""
cb_reset_starting_point!(self::CBQPSumModelBlock, y::CBMatrix, mu::Real, global_bundle::CBMinorantBundle, startindex_bundle::Integer) = @ccall libcb.cb_qpsummodelblock_reset_starting_point(self.data::Ptr{Cvoid}, y.data::Ptr{Cvoid}, mu::Cdouble, global_bundle.data::Ptr{Cvoid}, startindex_bundle::Cint)::Cint

@doc raw"""
    cb_compute_step!(self::CBQPSumModelBlock, ystep::CBMatrix, global_bundle::CBMinorantBundle, startindex_bundle::Integer)

compute the step in the model space given the step in the design space
"""
cb_compute_step!(self::CBQPSumModelBlock, ystep::CBMatrix, global_bundle::CBMinorantBundle, startindex_bundle::Integer) = @ccall libcb.cb_qpsummodelblock_compute_step(self.data::Ptr{Cvoid}, ystep.data::Ptr{Cvoid}, global_bundle.data::Ptr{Cvoid}, startindex_bundle::Cint)::Cint

@doc raw"""
    cb_computed_step!(self::CBQPSumModelBlock, modelxstep::CBMatrix, startindex_model::Integer, modelconstrstep::CBMatrix, startindex_constr::Integer)

store this computed step locally and compute the missing local dual step information
"""
cb_computed_step!(self::CBQPSumModelBlock, modelxstep::CBMatrix, startindex_model::Integer, modelconstrstep::CBMatrix, startindex_constr::Integer) = @ccall libcb.cb_qpsummodelblock_computed_step(self.data::Ptr{Cvoid}, modelxstep.data::Ptr{Cvoid}, startindex_model::Cint, modelconstrstep.data::Ptr{Cvoid}, startindex_constr::Cint)::Cint

@doc raw"""
    cb_do_step!(self::CBQPSumModelBlock, alpha::Real, y::CBMatrix, global_bundle::CBMinorantBundle, startindex_bundle::Integer)

move in the last computed step direction by a step of length alpha and compute and store the violation in this point for later use in
"""
cb_do_step!(self::CBQPSumModelBlock, alpha::Real, y::CBMatrix, global_bundle::CBMinorantBundle, startindex_bundle::Integer) = @ccall libcb.cb_qpsummodelblock_do_step(self.data::Ptr{Cvoid}, alpha::Cdouble, y.data::Ptr{Cvoid}, global_bundle.data::Ptr{Cvoid}, startindex_bundle::Cint)::Cint

@doc raw"""
    cb_add_localrhs!(self::CBQPSumModelBlock, globalrhs::CBMatrix, rhsmu::Real, rhscorr::Real, startindex_model::Integer, startindex_constraints::Integer, append::Bool, bundle::CBMinorantBundle, startindex_bundel::Integer)

If mu is not zero, always add the centering term for this mu as well, if append is false, add the Schur complement rhs for add_BtinvsysB, if append is true, fill in the rhs of the local system starting at startindex for the model and at startindex_constraints for the modelconstraints
"""
cb_add_localrhs!(self::CBQPSumModelBlock, globalrhs::CBMatrix, rhsmu::Real, rhscorr::Real, startindex_model::Integer, startindex_constraints::Integer, append::Bool, bundle::CBMinorantBundle, startindex_bundel::Integer) = @ccall libcb.cb_qpsummodelblock_add_localrhs(self.data::Ptr{Cvoid}, globalrhs.data::Ptr{Cvoid}, rhsmu::Cdouble, rhscorr::Cdouble, startindex_model::Cint, startindex_constraints::Cint, append::Cint, bundle.data::Ptr{Cvoid}, startindex_bundel::Cint)::Cint

@doc raw"""
    cb_add_BtinvsysB!(self::CBQPSumModelBlock, globalsys::CBSymmatrix, bundle::CBMinorantBundle, startindex_bundle::Integer)

add the "scaled" minorant outer products to globalsys, where the correct minroants start at the given index
"""
cb_add_BtinvsysB!(self::CBQPSumModelBlock, globalsys::CBSymmatrix, bundle::CBMinorantBundle, startindex_bundle::Integer) = @ccall libcb.cb_qpsummodelblock_add_btinvsysb(self.data::Ptr{Cvoid}, globalsys.data::Ptr{Cvoid}, bundle.data::Ptr{Cvoid}, startindex_bundle::Cint)::Cint

@doc raw"""
    cb_solve_constrsys!(self::CBQPSumModelBlock, ABchol::CBSymmatrix, LinvABrhs::CBMatrix, LinvABsol::CBMatrix, startindex_model::Integer, Crhs_and_sol::CBMatrix, startindex_constraints::Integer)

* @brief  given the Cholesky factorization LL' of *minus* the blocks A and B (contraints on design variables and Bundle-modelx) and LinvABrhs, solve for the local constraints C and add the new contribution of tracedual*LinvTrace to LinvABsol; store the tracedual in Crhs_and_sol but not yet locally (this will be done by computed_step() ).
     
"""
cb_solve_constrsys!(self::CBQPSumModelBlock, ABchol::CBSymmatrix, LinvABrhs::CBMatrix, LinvABsol::CBMatrix, startindex_model::Integer, Crhs_and_sol::CBMatrix, startindex_constraints::Integer) = @ccall libcb.cb_qpsummodelblock_solve_constrsys(self.data::Ptr{Cvoid}, ABchol.data::Ptr{Cvoid}, LinvABrhs.data::Ptr{Cvoid}, LinvABsol.data::Ptr{Cvoid}, startindex_model::Cint, Crhs_and_sol.data::Ptr{Cvoid}, startindex_constraints::Cint)::Cint

@doc raw"""
    cb_add_BCSchur_diagonal!(self::CBQPSumModelBlock, diagonal::CBMatrix, globalbundle::CBMinorantBundle, startindex_bundle::Integer)

* @brief add the diagonal of the Schur complemented blocks belonging to bundle and local constraints (used for diagonal preconditioning)
    
"""
cb_add_BCSchur_diagonal!(self::CBQPSumModelBlock, diagonal::CBMatrix, globalbundle::CBMinorantBundle, startindex_bundle::Integer) = @ccall libcb.cb_qpsummodelblock_add_bcschur_diagonal(self.data::Ptr{Cvoid}, diagonal.data::Ptr{Cvoid}, globalbundle.data::Ptr{Cvoid}, startindex_bundle::Cint)::Cint

@doc raw"""
    cb_propose_BCSchur_pcsubspace!(self::CBQPSumModelBlock, lowrank::CBMatrix, sigma_guess::CBMatrix, Diag_inv::CBMatrix, minval::Real, diaginvval::Real, globalbundle::CBMinorantBundle, startindex_bundle::Integer)

* @brief append to lowrank "large" columns that should serve well for generating a low rank projection of the Schur complemented model part. For each column i the coordinate sigma_guess(i) gives the Diag_inv-norm for this column. The parameter minval asks to ignore columns whose norms are smaller than minval. If diaginvval is positive, the vector Diag_inv is this value times the all ones vector.

      On input lowrank must have the correct number of rows already but may
      have 0 columns.
    
"""
cb_propose_BCSchur_pcsubspace!(self::CBQPSumModelBlock, lowrank::CBMatrix, sigma_guess::CBMatrix, Diag_inv::CBMatrix, minval::Real, diaginvval::Real, globalbundle::CBMinorantBundle, startindex_bundle::Integer) = @ccall libcb.cb_qpsummodelblock_propose_bcschur_pcsubspace(self.data::Ptr{Cvoid}, lowrank.data::Ptr{Cvoid}, sigma_guess.data::Ptr{Cvoid}, Diag_inv.data::Ptr{Cvoid}, minval::Cdouble, diaginvval::Cdouble, globalbundle.data::Ptr{Cvoid}, startindex_bundle::Cint)::Cint

@doc raw"""
    cb_prepare_BCSchur_JLprecond!(self::CBQPSumModelBlock, glob_lowrank::CBMatrix, subspace::CBMatrix, append_globtransp_times_mat_to_subspace::Bool, globalbundle::CBMinorantBundle, startindex_bundle::Integer)

* @brief compute the preconditioning low-rank representation of the Schur complementd blocks belonging to bundle and local constraints by adding a Johnson-Lindenstrauss projection onto the given subspace to glob_lowrank

    
"""
cb_prepare_BCSchur_JLprecond!(self::CBQPSumModelBlock, glob_lowrank::CBMatrix, subspace::CBMatrix, append_globtransp_times_mat_to_subspace::Bool, globalbundle::CBMinorantBundle, startindex_bundle::Integer) = @ccall libcb.cb_qpsummodelblock_prepare_bcschur_jlprecond(self.data::Ptr{Cvoid}, glob_lowrank.data::Ptr{Cvoid}, subspace.data::Ptr{Cvoid}, append_globtransp_times_mat_to_subspace::Cint, globalbundle.data::Ptr{Cvoid}, startindex_bundle::Cint)::Cint

@doc raw"""
    cb_add_Schur_rhs!(self::CBQPSumModelBlock, glob_rhs::CBMatrix, local_rhs::Union{<:CBMatrix,Nothing}, rhsmu::Real, rhscorr::Real, startindex_constraints::Integer, globalbundle::CBMinorantBundle, startindex_bundle::Integer)

add the contributions to glob_rhs of the Schur complemented model block, and return local_rhs of the non complemented constraint block in the rows/columns/diagonal block starting at startindex_constraints
"""
cb_add_Schur_rhs!(self::CBQPSumModelBlock, glob_rhs::CBMatrix, local_rhs::Union{<:CBMatrix,Nothing}, rhsmu::Real, rhscorr::Real, startindex_constraints::Integer, globalbundle::CBMinorantBundle, startindex_bundle::Integer) = @ccall libcb.cb_qpsummodelblock_add_schur_rhs(self.data::Ptr{Cvoid}, glob_rhs.data::Ptr{Cvoid}, (isnothing(local_rhs) ? C_NULL : local_rhs.data)::Ptr{Cvoid}, rhsmu::Cdouble, rhscorr::Cdouble, startindex_constraints::Cint, globalbundle.data::Ptr{Cvoid}, startindex_bundle::Cint)::Cint

@doc raw"""
    cb_add_Schur_mult!(self::CBQPSumModelBlock, in_vec::CBMatrix, out_vec::CBMatrix, in_cvec::Union{<:CBMatrix,Nothing}, out_cvec::Union{<:CBMatrix,Nothing}, startindex_constraints::Integer, globalbundle::CBMinorantBundle, startindex_bundle::Integer)

multiply in_vec with the local contribution to the global main block and add it to out_vec; the other local multiplications are carried out externally with the information provide in prepare_Schur_precond and are not done here.
"""
cb_add_Schur_mult!(self::CBQPSumModelBlock, in_vec::CBMatrix, out_vec::CBMatrix, in_cvec::Union{<:CBMatrix,Nothing}, out_cvec::Union{<:CBMatrix,Nothing}, startindex_constraints::Integer, globalbundle::CBMinorantBundle, startindex_bundle::Integer) = @ccall libcb.cb_qpsummodelblock_add_schur_mult(self.data::Ptr{Cvoid}, in_vec.data::Ptr{Cvoid}, out_vec.data::Ptr{Cvoid}, (isnothing(in_cvec) ? C_NULL : in_cvec.data)::Ptr{Cvoid}, (isnothing(out_cvec) ? C_NULL : out_cvec.data)::Ptr{Cvoid}, startindex_constraints::Cint, globalbundle.data::Ptr{Cvoid}, startindex_bundle::Cint)::Cint

@doc raw"""
    cb_computed_Schur_step!(self::CBQPSumModelBlock, xstep::CBMatrix, local_step::CBMatrix, startindex_model::Integer, globalbundle::CBMinorantBundle, startindex_bundle::Integer)

for passing on the solution information after solvin the system
"""
cb_computed_Schur_step!(self::CBQPSumModelBlock, xstep::CBMatrix, local_step::CBMatrix, startindex_model::Integer, globalbundle::CBMinorantBundle, startindex_bundle::Integer) = @ccall libcb.cb_qpsummodelblock_computed_schur_step(self.data::Ptr{Cvoid}, xstep.data::Ptr{Cvoid}, local_step.data::Ptr{Cvoid}, startindex_model::Cint, globalbundle.data::Ptr{Cvoid}, startindex_bundle::Cint)::Cint

@doc raw"""
    cb_dim_model!(self::CBQPSumModelBlock)

returns the dimension of the model set (here the same as the bundle size)
"""
cb_dim_model!(self::CBQPSumModelBlock) = @ccall libcb.cb_qpsummodelblock_dim_model(self.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_dim_constraints!(self::CBQPSumModelBlock)

returns the dimension of the system describing the model set (may contain further constraints)
"""
cb_dim_constraints!(self::CBQPSumModelBlock) = @ccall libcb.cb_qpsummodelblock_dim_constraints(self.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_constraints_cost!(self::CBQPSumModelBlock)

returns the dual upper bound to the model value (the trace weighted sum of the dual trace variables); it returns 0. if no model is contained
"""
cb_constraints_cost!(self::CBQPSumModelBlock) = @ccall libcb.cb_qpsummodelblock_constraints_cost(self.data::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_primalviol_2normsqr!(self::CBQPSumModelBlock)

return the squared Euclidean norm of constraint violation of modelx
"""
cb_primalviol_2normsqr!(self::CBQPSumModelBlock) = @ccall libcb.cb_qpsummodelblock_primalviol_2normsqr(self.data::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_dualviol_2normsqr!(self::CBQPSumModelBlock)

return the squared Euclidean norm of the dual model violation
"""
cb_dualviol_2normsqr!(self::CBQPSumModelBlock) = @ccall libcb.cb_qpsummodelblock_dualviol_2normsqr(self.data::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_get_mu_info(self::CBQPSumModelBlock)

add dimensions of the primal-dual pairs to mudim and add the "trace" (the inner product with center) of the respective primal-dual pair products for the current step; update the min and max values of x_i*z_i
"""
function cb_get_mu_info(self::CBQPSumModelBlock)
    max_xz = Ref{Float64}()
    min_xz = Ref{Float64}()
    tr_dxdz = Ref{Float64}()
    tr_xdzpdxz = Ref{Float64}()
    tr_xz = Ref{Float64}()
    mudim = Ref{Int}()
    @ccall libcb.cb_qpsummodelblock_get_mu_info(self.data::Ptr{Cvoid}, mudim::Ref{Int}, tr_xz::Ref{Float64}, tr_xdzpdxz::Ref{Float64}, tr_dxdz::Ref{Float64}, min_xz::Ref{Float64}, max_xz::Ref{Float64})::Cint
    return mudim[], tr_xz[], tr_xdzpdxz[], tr_dxdz[], min_xz[], max_xz[]
end

@doc raw"""
    cb_get_nbh_info(self::CBQPSumModelBlock, mudim::Integer, tr_xz::Real, tr_xdzpdxz::Real, tr_dxdz::Real, nbh_ubnd::Real)

for limiting the stepsize with respect to the neighborhood this information about norms and inner products of x(.)*z-tr_xz-tr_xz/mudim(.*)1, x.()*dz+dx(.)*z-tr_xdzpdxz/mudim(.*)1, and dx(.)*dz-tr_dxdz/mudim(.)*1 is required, each block *adds* its contribution to the numbers
"""
function cb_get_nbh_info(self::CBQPSumModelBlock, mudim::Integer, tr_xz::Real, tr_xdzpdxz::Real, tr_dxdz::Real, nbh_ubnd::Real)
    ip_dxdz_xdzpdxz = Ref{Float64}()
    ip_xz_dxdz = Ref{Float64}()
    ip_xz_xdzpdxz = Ref{Float64}()
    nrmsqr_dxdz = Ref{Float64}()
    nrmsqr_xdzpdxz = Ref{Float64}()
    nrmsqr_xz = Ref{Float64}()
    max_nbh = Ref{Float64}()
    alpha = Ref{Float64}()
    @ccall libcb.cb_qpsummodelblock_get_nbh_info(self.data::Ptr{Cvoid}, mudim::Cint, tr_xz::Cdouble, tr_xdzpdxz::Cdouble, tr_dxdz::Cdouble, nbh_ubnd::Cdouble, alpha::Ref{Float64}, max_nbh::Ref{Float64}, nrmsqr_xz::Ref{Float64}, nrmsqr_xdzpdxz::Ref{Float64}, nrmsqr_dxdz::Ref{Float64}, ip_xz_xdzpdxz::Ref{Float64}, ip_xz_dxdz::Ref{Float64}, ip_dxdz_xdzpdxz::Ref{Float64})::Cint
    return alpha[], max_nbh[], nrmsqr_xz[], nrmsqr_xdzpdxz[], nrmsqr_dxdz[], ip_xz_xdzpdxz[], ip_xz_dxdz[], ip_dxdz_xdzpdxz[]
end

@doc raw"""
    cb_linesearch(self::CBQPSumModelBlock)

if necessary, reduce alpha to the biggest value so that feasibility is maintained with this step size
"""
function cb_linesearch(self::CBQPSumModelBlock)
    alpha = Ref{Float64}()
    @ccall libcb.cb_qpsummodelblock_linesearch(self.data::Ptr{Cvoid}, alpha::Ref{Float64})::Cint
    return alpha[]
end

@doc raw"""
    cb_add_localsys!(self::CBQPSumModelBlock, globalsys::CBSymmatrix, startindex_model::Integer, startindex_constraints::Integer)

virtual int add_BtinvsysB(CH_Matrix_Classes::Symmatrix& globalsys);
"""
cb_add_localsys!(self::CBQPSumModelBlock, globalsys::CBSymmatrix, startindex_model::Integer, startindex_constraints::Integer) = @ccall libcb.cb_qpsummodelblock_add_localsys(self.data::Ptr{Cvoid}, globalsys.data::Ptr{Cvoid}, startindex_model::Cint, startindex_constraints::Cint)::Cint

@doc raw"""
    cb_localsys_mult!(self::CBQPSumModelBlock, in_vec::CBMatrix, out_vec::CBMatrix, startindex_model::Integer, startindex_constraints::Integer)

* @brief multiply the local system diagonal block consisting of the model and local contraints rows and columns by in_vec[startindex_model+0,...,+dim_model(),startindex_constraints+0,...,+dim_constraints]  into the same coordinates of out_vec. 
"""
cb_localsys_mult!(self::CBQPSumModelBlock, in_vec::CBMatrix, out_vec::CBMatrix, startindex_model::Integer, startindex_constraints::Integer) = @ccall libcb.cb_qpsummodelblock_localsys_mult(self.data::Ptr{Cvoid}, in_vec.data::Ptr{Cvoid}, out_vec.data::Ptr{Cvoid}, startindex_model::Cint, startindex_constraints::Cint)::Cint

