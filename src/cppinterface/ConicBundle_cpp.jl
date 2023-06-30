function cb_destroy!(obj, rest...)
    cb_destroy!(obj)
    cb_destroy!(rest...)
end

include("cb_matrix.jl")
include("cb_indexmatrix.jl")
include("cb_sparsemat.jl")
include("cb_symmatrix.jl")
include("cb_sparsesym.jl")
include("cb_cmgramdense.jl")
include("cb_cmgramsparse.jl")
include("cb_cmgramsparse_withoutdiag.jl")
include("cb_cmlowrankdd.jl")
include("cb_cmlowranksd.jl")
include("cb_cmlowrankss.jl")
include("cb_cmsingleton.jl")
include("cb_cmsymdense.jl")
include("cb_cmsymsparse.jl")
include("cb_sparsecoeffmatmatrix.jl")
include("cb_gb_rand.jl")
include("cb_coeffmatinfo.jl")
include("cb_primalmatrix.jl")
include("cb_matrixcbsolver.jl")
include("cb_blockpscprimal.jl")
include("cb_densepscprimal.jl")
include("cb_gramsparsepscprimal.jl")
include("cb_sparsepscprimal.jl")
include("cb_aftmodification.jl")
include("cb_groundsetmodification.jl")
include("cb_nncboxsupportmodification.jl")
include("cb_pscaffinemodification.jl")
include("cb_socsupportmodification.jl")
include("cb_affinefunctiontransformation.jl")
include("cb_minorantpointer.jl")
include("cb_minorantusedata.jl")
include("cb_minorant.jl")
include("cb_boxprimalextender.jl")
include("cb_boxoracle.jl")
include("cb_pscprimalextender.jl")
include("cb_pscbundleparameters.jl")
include("cb_socprimalextender.jl")
include("cb_socbundleparameters.jl")
include("cb_cfunctionminorantextender.jl")
include("cb_cfunction.jl")
include("cb_nncboxsupportminorantextender.jl")
include("cb_nncboxsupportfunction.jl")
include("cb_pscaffineminorantextender.jl")
include("cb_pscaffinefunction.jl")
include("cb_socsupportminorantextender.jl")
include("cb_socsupportfunction.jl")
include("cb_boxmodelparameters.jl")
include("cb_nncmodelparameters.jl")
include("cb_pscmodelparameters.jl")
include("cb_socmodelparameters.jl")
include("cb_sumbundleparameters.jl")
include("cb_aftdata.jl")
include("cb_boxdata.jl")
include("cb_nncdata.jl")
include("cb_pscdata.jl")
include("cb_socdata.jl")
include("cb_bundlehkweight.jl")
include("cb_bundlerqbweight.jl")
include("cb_bundledensetrustregionprox.jl")
include("cb_bundlediagonaltrustregionprox.jl")
include("cb_bundledlrtrustregionprox.jl")
include("cb_bundleidprox.jl")
include("cb_bundlelowranktrustregionprox.jl")
include("cb_qpsolverparameters.jl")
include("cb_qpsolver.jl")
include("cb_uqpsolver.jl")
include("cb_lpgroundset.jl")
include("cb_unconstrainedgroundset.jl")
include("cb_aftmodel.jl")
include("cb_boxmodel.jl")
include("cb_nncmodel.jl")
include("cb_pscmodel.jl")
include("cb_socmodel.jl")
include("cb_summodel.jl")
include("cb_pscvariablemetricselection.jl")
include("cb_variablemetricsvdselection.jl")
include("cb_qpdirectkktsolver.jl")
include("cb_qpiterativekkthaeqsolver.jl")
include("cb_qpiterativekkthasolver.jl")
include("cb_qpkktsolvercomparison.jl")
include("cb_sumbundlehandler.jl")
include("cb_qpconemodelblock.jl")
include("cb_qpsummodelblock.jl")
include("cb_uqpconemodelblock.jl")
include("cb_uqpsummodelblock.jl")
include("cb_minres.jl")
include("cb_pcg.jl")
include("cb_psqmr.jl")
include("cb_qpkktsubspacehprecond.jl")
include("cb_sumbundle.jl")
include("cb_bundlesolver.jl")
include("cb_bundleterminator.jl")
include("cb_clock.jl")
include("cb_microseconds.jl")
export CBAFTData,
    CBAFTModel,
    CBAFTModification,
    CBAffineFunctionTransformation,
    CBBlockPSCPrimal,
    CBBoxData,
    CBBoxModel,
    CBBoxModelParameters,
    CBBoxOracle,
    CBBundleDLRTrustRegionProx,
    CBBundleDenseTrustRegionProx,
    CBBundleDiagonalTrustRegionProx,
    CBBundleHKWeight,
    CBBundleIdProx,
    CBBundleLowRankTrustRegionProx,
    CBBundleRQBWeight,
    CBBundleSolver,
    CBBundleTerminator,
    CBCFunction,
    CBCFunctionMinorantExtender,
    CBCMgramdense,
    CBCMgramsparse,
    CBCMgramsparse_withoutdiag,
    CBCMlowrankdd,
    CBCMlowranksd,
    CBCMlowrankss,
    CBCMsingleton,
    CBCMsymdense,
    CBCMsymsparse,
    CBClock,
    CBCoeffmatInfo,
    CBDensePSCPrimal,
    CBGB_rand,
    CBGramSparsePSCPrimal,
    CBGroundsetModification,
    CBIndexmatrix,
    CBLPGroundset,
    CBMatrix,
    CBMatrixCBSolver,
    CBMicroseconds,
    CBMinRes,
    CBMinorant,
    CBMinorantPointer,
    CBMinorantUseData,
    CBMode,
    CBModelUpdate,
    CBNNCBoxSupportFunction,
    CBNNCBoxSupportMinorantExtender,
    CBNNCBoxSupportModification,
    CBNNCData,
    CBNNCModel,
    CBNNCModelParameters,
    CBPCG,
    CBPSCAffineFunction,
    CBPSCAffineMinorantExtender,
    CBPSCAffineModification,
    CBPSCBundleParameters,
    CBPSCData,
    CBPSCModel,
    CBPSCModelParameters,
    CBPSCVariableMetricSelection,
    CBPrimalMatrix,
    CBPsqmr,
    CBQPConeModelBlock,
    CBQPDirectKKTSolver,
    CBQPIterativeKKTHASolver,
    CBQPIterativeKKTHAeqSolver,
    CBQPKKTSolverComparison,
    CBQPKKTSubspaceHPrecond,
    CBQPSolver,
    CBQPSolverParameters,
    CBQPSumModelBlock,
    CBSOCData,
    CBSOCModel,
    CBSOCModelParameters,
    CBSOCSupportFunction,
    CBSOCSupportMinorantExtender,
    CBSOCSupportModification,
    CBSparseCoeffmatMatrix,
    CBSparsePSCPrimal,
    CBSparsemat,
    CBSparsesym,
    CBSumBundle,
    CBSumBundleParameters,
    CBSumModel,
    CBSymmatrix,
    CBUQPConeModelBlock,
    CBUQPSolver,
    CBUQPSumModelBlock,
    CBUnconstrainedGroundset,
    CBVariableMetricSVDSelection,
    cb_Aasen_Lsolve,
    cb_Aasen_Ltsolve,
    cb_Aasen_factor!,
    cb_Aasen_solve,
    cb_Aasen_tridiagsolve,
    cb_B_times!,
    cb_Chol_Lmult,
    cb_Chol_Lsolve,
    cb_Chol_Ltmult,
    cb_Chol_Ltsolve,
    cb_Chol_factor!,
    cb_Chol_inverse,
    cb_Chol_scaleLi,
    cb_Chol_scaleLt,
    cb_Chol_solve,
    cb_Diag,
    cb_Gram_ip,
    cb_ItSys_mult!,
    cb_LDLfactor!,
    cb_LDLinverse,
    cb_LDLsolve,
    cb_QPallow_UQPSolver!,
    cb_QPapply_modification!,
    cb_QPclear!,
    cb_QPconstrained,
    cb_QPensure_feasibility!,
    cb_QPget_KKTsolver!,
    cb_QPget_blockA_norm!,
    cb_QPget_blockH_norm!,
    cb_QPget_dual_infeasibility_eps,
    cb_QPget_lower_bound,
    cb_QPget_lower_bound!,
    cb_QPget_lower_bound_gap_eps,
    cb_QPget_maxiter,
    cb_QPget_min_objective_relprec,
    cb_QPget_nbh_lb,
    cb_QPget_nbh_ub,
    cb_QPget_objective_gap_eps,
    cb_QPget_primal_infeasibility_eps,
    cb_QPget_solution!,
    cb_QPget_system_size!,
    cb_QPget_upper_bound,
    cb_QPget_upper_bound_gap_eps,
    cb_QPget_use_neighborhood,
    cb_QPget_use_predictor_corrector,
    cb_QPget_use_socqp,
    cb_QPinit_KKTdata!,
    cb_QPinit_KKTsystem!,
    cb_QPis_feasible!,
    cb_QPprefer_UQPSolver,
    cb_QPprint_statistics!,
    cb_QPresolve!,
    cb_QPset_KKTsolver!,
    cb_QPset_allow_UQPSolver!,
    cb_QPset_dual_infeasibility_eps!,
    cb_QPset_lower_and_upper_bounds!,
    cb_QPset_lower_bound_gap_eps!,
    cb_QPset_maxiter!,
    cb_QPset_min_objective_relprec!,
    cb_QPset_nbh_bounds!,
    cb_QPset_objective_gap_eps!,
    cb_QPset_parameters!,
    cb_QPset_primal_infeasibility_eps!,
    cb_QPset_upper_bound_gap_eps!,
    cb_QPset_use_neighborhood!,
    cb_QPset_use_predictor_corrector!,
    cb_QPset_use_socqp!,
    cb_QPsolve!,
    cb_QPsolve_KKTsystem!,
    cb_QPstart_modification!,
    cb_QPsupports_updates!,
    cb_QPsupports_yfixing!,
    cb_QPupdate!,
    cb_QR_concat_right!,
    cb_QR_factor,
    cb_QR_factor!,
    cb_QR_factor_relpiv!,
    cb_QR_solve!,
    cb_Q_times,
    cb_Qt_times,
    cb_abs,
    cb_abs!,
    cb_active,
    cb_add_BCSchur_diagonal!,
    cb_add_BDBt!,
    cb_add_Bs,
    cb_add_BtinvsysB!,
    cb_add_H,
    cb_add_Hx,
    cb_add_Schur_mult!,
    cb_add_Schur_rhs!,
    cb_add_append_blocks!,
    cb_add_append_rows!,
    cb_add_append_variables!,
    cb_add_append_vars!,
    cb_add_apply_factor!,
    cb_add_coeff!,
    cb_add_coeffs!,
    cb_add_contributions!,
    cb_add_delete_blocks!,
    cb_add_delete_rows!,
    cb_add_delete_vars!,
    cb_add_diagonal_scaling,
    cb_add_function!,
    cb_add_local_sys!,
    cb_add_localrhs!,
    cb_add_localsys!,
    cb_add_model!,
    cb_add_modelx_aggregate!,
    cb_add_offset!,
    cb_add_projection,
    cb_add_reassign_blocks!,
    cb_add_reassign_rows!,
    cb_add_reassign_variables!,
    cb_add_reassign_vars!,
    cb_add_reset_generating_primal!,
    cb_add_solver!,
    cb_add_variable_metric!,
    cb_add_xinv_kron_z!,
    cb_addmeto,
    cb_addprodto,
    cb_adjust_multiplier!,
    cb_adjust_trace!,
    cb_aggregate,
    cb_aggregate!,
    cb_aggregate_Gram_matrix!,
    cb_aggregate_primal_data!,
    cb_aggregated!,
    cb_analyze_modification,
    cb_append!,
    cb_append_blocks!,
    cb_append_columns!,
    cb_append_constraints!,
    cb_append_to_old,
    cb_append_variables!,
    cb_appended_blockdim,
    cb_appended_rowdim,
    cb_appended_vardim,
    cb_apply_Hinv,
    cb_apply_modification!,
    cb_apply_modified_transform,
    cb_apply_to_PSCAffine,
    cb_apply_to_bounds,
    cb_apply_to_costs,
    cb_apply_to_factor,
    cb_apply_to_offset,
    cb_apply_to_vars,
    cb_apply_variable_metric!,
    cb_argument_changes,
    cb_assign_Gram_matrix!,
    cb_block,
    cb_block_modifications,
    cb_blockdim,
    cb_bundle_size,
    cb_call_primal_extender!,
    cb_candidate!,
    cb_ceil,
    cb_ceil!,
    cb_center_modified!,
    cb_check_center_validity_by_candidate!,
    cb_check_correctness,
    cb_check_support,
    cb_clear!,
    cb_clear_aggregates!,
    cb_clear_cand_minorants!,
    cb_clear_fail_counts!,
    cb_clear_fails!,
    cb_clear_model!,
    cb_clear_terminated!,
    cb_clone,
    cb_clone!,
    cb_clone_BundleParameters,
    cb_clone_VariableMetricSelection!,
    cb_clone_minorant,
    cb_coeff,
    cb_coeff!,
    cb_col,
    cb_col_nonzeros,
    cb_coldim,
    cb_colhouse,
    cb_colip,
    cb_cols,
    cb_colsip,
    cb_compute!,
    cb_compute_QP_costs!,
    cb_compute_local_directions!,
    cb_compute_step!,
    cb_computed_Schur_step!,
    cb_computed_step!,
    cb_concat_below,
    cb_concat_below!,
    cb_concat_right,
    cb_concat_right!,
    cb_cond_number_mult!,
    cb_constrained,
    cb_constraints_cost!,
    cb_contains_nan!,
    cb_contains_support,
    cb_contribute_initial_bundle!,
    cb_contribute_new_minorants!,
    cb_copy_traforows,
    cb_delete_blocks!,
    cb_delete_cols!,
    cb_delete_columns!,
    cb_delete_principal_submatrix!,
    cb_delete_rows!,
    cb_delete_variables!,
    cb_deleted_block_indices,
    cb_deleted_row_indices,
    cb_deleted_var_indices,
    cb_deleted_variables_are_zero,
    cb_dense,
    cb_descent_update!,
    cb_destroy!,
    cb_diag,
    cb_diagonal_bounds_scaling_update!,
    cb_diagonal_scaling_heuristic_update!,
    cb_dim,
    cb_dim2,
    cb_dim_constraints!,
    cb_dim_model!,
    cb_display,
    cb_display_model_values!,
    cb_dnorm_sqr,
    cb_do_step!,
    cb_dual_norm_squared,
    cb_dualviol_2normsqr!,
    cb_eig,
    cb_elapsed_time,
    cb_empty,
    cb_enlarge_below!,
    cb_enlarge_right!,
    cb_ensure_feasibility!,
    cb_equal,
    cb_equals,
    cb_eval_function!,
    cb_eval_model,
    cb_eval_model!,
    cb_evaluate,
    cb_evaluate_projection!,
    cb_evaluate_trace,
    cb_extend!,
    cb_extend_Box!,
    cb_extend_Ritz!,
    cb_extend_SOC!,
    cb_extract_SOCvector!,
    cb_find,
    cb_find_number,
    cb_floor,
    cb_floor!,
    cb_form_bundlevecs!,
    cb_from_dim,
    cb_generate_minorant!,
    cb_genmult,
    cb_get_A,
    cb_get_Ab,
    cb_get_Bt!,
    cb_get_C!,
    cb_get_D,
    cb_get_H,
    cb_get_Hchol,
    cb_get_Q,
    cb_get_QPcoeff_time,
    cb_get_QPsolve_time,
    cb_get_Ritz_values,
    cb_get_UVlambda!,
    cb_get_X!,
    cb_get_Z!,
    cb_get_activedim,
    cb_get_add_offset,
    cb_get_additional_factor,
    cb_get_additional_offset,
    cb_get_aft,
    cb_get_aggr_dnormsqr,
    cb_get_aggregate,
    cb_get_aggregate_offset,
    cb_get_append_cols,
    cb_get_append_costs,
    cb_get_append_rhs,
    cb_get_append_rows,
    cb_get_appended_vardim,
    cb_get_approximate_primal,
    cb_get_approximate_slacks,
    cb_get_arg_offset,
    cb_get_arg_trafo,
    cb_get_augvalfails,
    cb_get_augvalfailslimit,
    cb_get_avg_reduction,
    cb_get_block_append,
    cb_get_boxx!,
    cb_get_bundle,
    cb_get_bundle_data,
    cb_get_bundle_parameters,
    cb_get_bundleweight,
    cb_get_c,
    cb_get_cand_gs_val,
    cb_get_cand_minorant,
    cb_get_cand_objval,
    cb_get_cand_ub,
    cb_get_cand_y,
    cb_get_candidate,
    cb_get_candidate_primal,
    cb_get_candidate_value,
    cb_get_center,
    cb_get_center_gs_val,
    cb_get_center_minorant!,
    cb_get_center_objval,
    cb_get_center_primal,
    cb_get_center_ub,
    cb_get_center_y,
    cb_get_cntobjeval,
    cb_get_coeff,
    cb_get_colindex,
    cb_get_colinfo,
    cb_get_colval,
    cb_get_constant_minorant,
    cb_get_contributed_model_aggregate,
    cb_get_corr!,
    cb_get_costs,
    cb_get_cutval,
    cb_get_data,
    cb_get_data!,
    cb_get_dense_cnt,
    cb_get_dense_coeff_store!,
    cb_get_descent_step,
    cb_get_descent_steps,
    cb_get_dim,
    cb_get_do_variable_metric,
    cb_get_dualval,
    cb_get_edge,
    cb_get_edge_rep,
    cb_get_err,
    cb_get_eval_time,
    cb_get_evalaugmodel_time,
    cb_get_factored,
    cb_get_fixed_active_bounds,
    cb_get_fun_coeff,
    cb_get_fun_offset,
    cb_get_function_factor,
    cb_get_function_minorant!,
    cb_get_function_status,
    cb_get_generating_primal,
    cb_get_generating_primal!,
    cb_get_grammatrix,
    cb_get_groundset,
    cb_get_groundset_id,
    cb_get_gs_aggregate,
    cb_get_gs_minorant,
    cb_get_ijval,
    cb_get_increase_factor,
    cb_get_infinity,
    cb_get_innerit,
    cb_get_iter,
    cb_get_keepsize,
    cb_get_last_alpha!,
    cb_get_last_weight,
    cb_get_latest_minorants!,
    cb_get_lbindex,
    cb_get_lbounds,
    cb_get_lby,
    cb_get_linear_cost,
    cb_get_lmin_invM1!,
    cb_get_local_dualcost,
    cb_get_local_model_aggregate,
    cb_get_local_primalcost,
    cb_get_lower_bounds!,
    cb_get_make_aggr_time,
    cb_get_map_to_old_variables,
    cb_get_maxeigval_factor,
    cb_get_maxit,
    cb_get_maxiter,
    cb_get_maxweight,
    cb_get_mineigval_factor,
    cb_get_minorant,
    cb_get_minweight,
    cb_get_mode,
    cb_get_model,
    cb_get_model_aggregate,
    cb_get_model_aggregate!,
    cb_get_model_calls_delete!,
    cb_get_model_data,
    cb_get_modeldcstr!,
    cb_get_modeldx!,
    cb_get_modeleps,
    cb_get_modelfails,
    cb_get_modelfailslimit,
    cb_get_modelval,
    cb_get_modelx!,
    cb_get_modification_id,
    cb_get_mu!,
    cb_get_mu_info,
    cb_get_mu_stats!,
    cb_get_n_contributors,
    cb_get_n_descent_steps,
    cb_get_n_functions,
    cb_get_n_inner_iterations,
    cb_get_n_inner_updates,
    cb_get_n_latest_minorants,
    cb_get_n_oracle_calls,
    cb_get_nbh_info,
    cb_get_nblocks,
    cb_get_new_index,
    cb_get_new_vardim,
    cb_get_next_weight,
    cb_get_next_weight_set,
    cb_get_nmult,
    cb_get_nncx!,
    cb_get_nncz!,
    cb_get_null_step,
    cb_get_objevallimit,
    cb_get_objval,
    cb_get_offset,
    cb_get_offset_append,
    cb_get_old_X!,
    cb_get_old_Z!,
    cb_get_old_mu!,
    cb_get_old_nncx!,
    cb_get_old_nncz!,
    cb_get_old_s!,
    cb_get_old_socx!,
    cb_get_old_socz!,
    cb_get_old_vardim,
    cb_get_old_y!,
    cb_get_oldfactor,
    cb_get_opAt!,
    cb_get_oracle_object!,
    cb_get_oraclefails,
    cb_get_oraclefailslimit,
    cb_get_positive,
    cb_get_posteval_time,
    cb_get_precond_rank!,
    cb_get_preeval_time,
    cb_get_primal,
    cb_get_primal!,
    cb_get_primaleigs,
    cb_get_primalval,
    cb_get_primalvecs,
    cb_get_prob_stats!,
    cb_get_prox,
    cb_get_pscx!,
    cb_get_qp_solver!,
    cb_get_qpfails,
    cb_get_qpfailslimit,
    cb_get_recomp,
    cb_get_recomplimit,
    cb_get_reset_primal,
    cb_get_residual_norm,
    cb_get_ret_code,
    cb_get_rhslb,
    cb_get_rhslbind,
    cb_get_rhsub,
    cb_get_rhsubind,
    cb_get_rowindex,
    cb_get_rowinfo,
    cb_get_rowval,
    cb_get_s!,
    cb_get_scalefactor,
    cb_get_selection_method,
    cb_get_sgnorm,
    cb_get_shallowcut,
    cb_get_skip_extension,
    cb_get_skippedsize,
    cb_get_socdim!,
    cb_get_socx!,
    cb_get_socz!,
    cb_get_solver,
    cb_get_starting_point,
    cb_get_status,
    cb_get_store,
    cb_get_subgradient,
    cb_get_sumaugvalfails,
    cb_get_sumbundle,
    cb_get_suminnerit,
    cb_get_summodelfails,
    cb_get_sumoraclefails,
    cb_get_sumqpfails,
    cb_get_sumrecomp,
    cb_get_sumupdatecnt,
    cb_get_suppcol,
    cb_get_suppind,
    cb_get_sysviol_constraints!,
    cb_get_sysviol_model!,
    cb_get_t_precond_mult!,
    cb_get_term_corr,
    cb_get_termeps,
    cb_get_terminate,
    cb_get_terminated,
    cb_get_terminator,
    cb_get_termprec,
    cb_get_timelimit,
    cb_get_topvecs,
    cb_get_trace!,
    cb_get_ubindex,
    cb_get_ubounds,
    cb_get_uby,
    cb_get_upper_bounds!,
    cb_get_use_linval,
    cb_get_use_variable_metric,
    cb_get_use_yfixing,
    cb_get_var_append,
    cb_get_variable_metric_selection,
    cb_get_weight,
    cb_get_weightu,
    cb_get_x!,
    cb_get_y!,
    cb_get_yfixed,
    cb_gramip,
    cb_groundset_changes_suffice_for_identity!,
    cb_guess_curvature,
    cb_handles!,
    cb_has_bundle_data,
    cb_has_bundle_for,
    cb_has_contributions,
    cb_has_roots,
    cb_has_working_roots,
    cb_hhmmss,
    cb_hhmmssdd,
    cb_house,
    cb_ignore_groundset_modification,
    cb_incorporate!,
    cb_init!,
    cb_init_data!,
    cb_init_diag!,
    cb_init_problem!,
    cb_init_size!,
    cb_init_support!,
    cb_init_svec!,
    cb_init_system!,
    cb_initialization_needed,
    cb_initialize!,
    cb_inner_line_search!,
    cb_insert_col!,
    cb_insert_row!,
    cb_install_external_aggregate!,
    cb_inv,
    cb_inv!,
    cb_ip,
    cb_ip_min_max,
    cb_is_DLR,
    cb_is_feasible!,
    cb_lb_function!,
    cb_lb_model,
    cb_left_genmult,
    cb_left_right_prod,
    cb_left_right_product!,
    cb_line_search!,
    cb_linesearch,
    cb_localsys_mult!,
    cb_ls!,
    cb_make_model_aggregate!,
    cb_make_symmatrix,
    cb_map_to_old_blocks,
    cb_map_to_old_rows,
    cb_map_to_old_variables,
    cb_mapped_variables_are_equal,
    cb_max,
    cb_maxcols,
    cb_maxrows,
    cb_mfile_data,
    cb_mfile_output,
    cb_min,
    cb_mincols,
    cb_minrows,
    cb_model,
    cb_model_aggregate_modified!,
    cb_modified_transform_argument,
    cb_multiply!,
    cb_new_block_indices,
    cb_new_blockdim,
    cb_new_initial_oraclemodification,
    cb_new_row_indices,
    cb_new_rowdim,
    cb_new_var_indices,
    cb_new_vardim,
    cb_new_variables_are_zero,
    cb_newsize!,
    cb_next!,
    cb_nnls,
    cb_no_additions_or_deletions_in_rows,
    cb_no_additions_or_deletions_in_vars,
    cb_no_modification,
    cb_nonzeros,
    cb_nonzeros!,
    cb_norm,
    cb_norm2,
    cb_normDsquared,
    cb_norm_sqr,
    cb_norm_squared,
    cb_normalize_sumbundle!,
    cb_nsubmodels,
    cb_nullstep_update!,
    cb_number_aggregated,
    cb_nzcoldim,
    cb_objective_value,
    cb_offset,
    cb_offset_gives_value_at_origin,
    cb_offset_gives_value_at_origin!,
    cb_old_blockdim,
    cb_old_rowdim,
    cb_old_vardim,
    cb_one_user,
    cb_only_scalars_change,
    cb_out,
    cb_output_aft_data,
    cb_output_bundle_data,
    cb_pivot_permute!,
    cb_pivot_permute_cols!,
    cb_pivot_permute_rows!,
    cb_pop_aft!,
    cb_postgenmult,
    cb_precondM1!,
    cb_precond_invG1!,
    cb_precond_invG1tran!,
    cb_precond_size!,
    cb_pregenmult,
    cb_prepare_BCSchur_JLprecond!,
    cb_preserves_identity,
    cb_primal_ip,
    cb_primalviol_2normsqr!,
    cb_principal_submatrix,
    cb_print_id,
    cb_print_line_summary,
    cb_print_problem_data,
    cb_print_problem_data_to_mfile,
    cb_print_statistics,
    cb_print_status,
    cb_print_termination_code,
    cb_print_time,
    cb_prodvec_flops,
    cb_project,
    cb_projected_clone!,
    cb_projection!,
    cb_propose_BCSchur_pcsubspace!,
    cb_provide_model_aggregate!,
    cb_push_aft!,
    cb_qp_cost_indices,
    cb_qp_mfile_data,
    cb_rand,
    cb_rand!,
    cb_rand_normal!,
    cb_rank2add,
    cb_rankadd,
    cb_reassign_blocks!,
    cb_reassign_coeffs!,
    cb_reassign_columns!,
    cb_reassign_variables!,
    cb_recompute_center!,
    cb_recursive_copy_data_of!,
    cb_recursive_delete_and_clear!,
    cb_reduce_length!,
    cb_reinit_function_model!,
    cb_remove_contributions!,
    cb_remove_model!,
    cb_reset_function_factor!,
    cb_reset_starting_point!,
    cb_reset_t_precond_mult!,
    cb_resolve!,
    cb_restart_x!,
    cb_restart_y!,
    cb_right_genmult,
    cb_rint,
    cb_rint!,
    cb_round,
    cb_round!,
    cb_roundhundredths,
    cb_roundsecs,
    cb_row,
    cb_row_nonzeros,
    cb_rowdim,
    cb_rowhouse,
    cb_rowip,
    cb_rows,
    cb_rowsip,
    cb_save,
    cb_scale!,
    cb_scale_cols!,
    cb_scale_minorant!,
    cb_scale_primal_data!,
    cb_scale_rows!,
    cb_scaled_index,
    cb_scaled_index_subset,
    cb_scaledrankadd,
    cb_scaling_indices,
    cb_select_model!,
    cb_set!,
    cb_set_D!,
    cb_set_active_bounds_fixing!,
    cb_set_aggr_dnormsqr!,
    cb_set_append_to_old!,
    cb_set_augvalfailslimit!,
    cb_set_bundle_parameters!,
    cb_set_bundleweight!,
    cb_set_cand_minorant!,
    cb_set_cbout!,
    cb_set_check_correctness!,
    cb_set_clock!,
    cb_set_data!,
    cb_set_defaults!,
    cb_set_do_yfixing!,
    cb_set_eval_limit!,
    cb_set_fun_coeff!,
    cb_set_fun_offset!,
    cb_set_groundset_id!,
    cb_set_infinity!,
    cb_set_inner_update_limit!,
    cb_set_lower_bound!,
    cb_set_mL!,
    cb_set_mN!,
    cb_set_max_Ritzvecs!,
    cb_set_max_bundlesize!,
    cb_set_max_modelsize!,
    cb_set_max_new!,
    cb_set_max_updates!,
    cb_set_max_weight!,
    cb_set_maxeigval_factor!,
    cb_set_maxit!,
    cb_set_maxiter!,
    cb_set_maxweight!,
    cb_set_min_weight!,
    cb_set_mineigval_factor!,
    cb_set_minweight!,
    cb_set_model!,
    cb_set_model_calls_delete!,
    cb_set_modeleps!,
    cb_set_modelfailslimit!,
    cb_set_modification_id!,
    cb_set_n_latest_minorants!,
    cb_set_new_center!,
    cb_set_new_center_point!,
    cb_set_next_weight!,
    cb_set_nullstep_updates!,
    cb_set_objevallimit!,
    cb_set_offset!,
    cb_set_oldfactor!,
    cb_set_oracle!,
    cb_set_oraclefailslimit!,
    cb_set_out!,
    cb_set_parent_information!,
    cb_set_point!,
    cb_set_primal!,
    cb_set_prox!,
    cb_set_prox_diagonal!,
    cb_set_qp_solver!,
    cb_set_qp_solver_parameters!,
    cb_set_qp_xstart!,
    cb_set_qp_ystart!,
    cb_set_qpfailslimit!,
    cb_set_qpsolver!,
    cb_set_recomplimit!,
    cb_set_scalefactor!,
    cb_set_selection_method!,
    cb_set_skip_extension!,
    cb_set_starting_point!,
    cb_set_subspace!,
    cb_set_sumbundle!,
    cb_set_sumbundle_parameters!,
    cb_set_term_relprec!,
    cb_set_termbounds!,
    cb_set_termeps!,
    cb_set_terminator!,
    cb_set_time_limit!,
    cb_set_timelimit!,
    cb_set_tol!,
    cb_set_upper_bound!,
    cb_set_use_linval!,
    cb_set_use_yfixing!,
    cb_set_variable_metric!,
    cb_set_variable_metric_selection!,
    cb_set_weight_update!,
    cb_set_weightu!,
    cb_set_yfixed!,
    cb_shift_diag!,
    cb_shuffle!,
    cb_sign,
    cb_sign!,
    cb_skron,
    cb_solve!,
    cb_solve_constrsys!,
    cb_sortindex,
    cb_sparse,
    cb_sparseDiag,
    cb_sparse_argument_changes,
    cb_sparsemult,
    cb_sparsify!,
    cb_sqr,
    cb_sqr!,
    cb_sqrt,
    cb_sqrt!,
    cb_start!,
    cb_start_augmodel!,
    cb_start_modification!,
    cb_start_sumaugmodel!,
    cb_starting_x!,
    cb_starting_y!,
    cb_store_SOCvec!,
    cb_store_svec,
    cb_subassign!,
    cb_subspace,
    cb_subtract_z,
    cb_suggest_mu!,
    cb_sum,
    cb_sumbundle_mode!,
    cb_sumcols,
    cb_sumrows,
    cb_support_in,
    cb_support_rankadd,
    cb_support_xbpeya!,
    cb_supports_dense_variable_metric,
    cb_supports_diagonal_bounds_scaling,
    cb_supports_diagonal_variable_metric,
    cb_supports_lowrank_variable_metric,
    cb_svec,
    cb_svec_projection!,
    cb_sveci,
    cb_swap,
    cb_swap_colsij!,
    cb_swap_rowsij!,
    cb_swapij!,
    cb_symscale,
    cb_synchronize_ids!,
    cb_termination_code,
    cb_time,
    cb_times_B!,
    cb_times_Q,
    cb_to_dim,
    cb_trace,
    cb_tracedual,
    cb_transform_argument,
    cb_transform_minorant,
    cb_transform_minorants,
    cb_transpose,
    cb_transpose!,
    cb_tril,
    cb_tril!,
    cb_tril_solve!,
    cb_triu,
    cb_triu!,
    cb_triu_solve!,
    cb_unif_long!,
    cb_update!,
    cb_update_QP_costs!,
    cb_update_model!,
    cb_valid,
    cb_variable_modifications,
    cb_wall_time,
    cb_weight_changed,
    cb_xbpeya,
    cb_xdim,
    cb_xetriu_yza!,
    cb_xeya!,
    cb_xeyapzb,
    cb_xpetriu_yza!,
    cb_xpeya!,
    cb_ydim,
    cb_zero,
    cbm_child,
    cbm_inactive,
    cbm_root,
    cbm_unavailable,
    cbmu_descent_step,
    cbmu_new_subgradient,
    cbmu_null_step