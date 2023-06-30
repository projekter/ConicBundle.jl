@doc raw"""
    CBPSCModelParameters(modelsize::Integer, bundlesize::Integer = 10, updaterule::Integer = 0, incr::Integer = -1)

constructor for size parameters
"""
CBPSCModelParameters(modelsize::Integer, bundlesize::Integer = 10, updaterule::Integer = 0, incr::Integer = -1) = CBPSCModelParameters(@ccall libcb.cb_pscmodelparameters_new2(modelsize::Cint, bundlesize::Cint, updaterule::Cint, incr::Cint)::Ptr{Cvoid})

@doc raw"""
    CBPSCModelParameters(bp::CBBundleParameters, incr::Integer = -1)

copy constructor for BundleParameters
"""
CBPSCModelParameters(bp::CBBundleParameters, incr::Integer = -1) = CBPSCModelParameters(@ccall libcb.cb_pscmodelparameters_new3(bp.data::Ptr{Cvoid}, incr::Cint)::Ptr{Cvoid})

@doc raw"""
    CBPSCModelParameters(sms::CBPSCModelParameters)

copy constructor
"""
CBPSCModelParameters(sms::CBPSCModelParameters) = CBPSCModelParameters(@ccall libcb.cb_pscmodelparameters_new4(sms.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_clone_BundleParameters(self::CBPSCModelParameters)

clone
"""
cb_clone_BundleParameters(self::CBPSCModelParameters) = CBBundleParameters(@ccall libcb.cb_pscmodelparameters_clone_bundleparameters(self.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_select_model!(self::CBPSCModelParameters, modelvecs::CBMatrix, model_aggregate::CBMinorantPointer, topvecs::CBMatrix, Ritz_values::CBMatrix, primal_Ritzval::Real, primaleigs::CBMatrix, primalvecs::CBMatrix, primal_aggregate::CBMinorantPointer, primal_aggregate_coeff::Real, growthrate::Real, primalgrowth::CBMatrix, dualgrowth::CBMatrix, cand_Ritzvec::CBMatrix, cand_Ritzval::CBMatrix, oracle::Union{<:CBPSCOracle,Nothing}, modification_id::Integer, function_task::CBFunctionTask, function_factor::Real, model_update::CBModelUpdate, center_id::Integer, center_y::CBMatrix, cand_id::Integer, cand_y::CBMatrix, model_maxviol::Real, diffval_center_aggregate::Real, H::CBBundleProxObject)

* @brief PSCModel calls this for selecting the next positive semidefinite model

      @param[in,out] modelvecs
          * on input: orthonormal basis of the subspace used in the last
            semidefinite model
          * on output: orthonormal subspace basis to be used in the next
            semidefinite model, it includes the subspace of the primalvecs of
            the keepsize largest primaleigs. In particular for null steps
      the remaining columns of primalvecs with nonzero primaleigs
            have to be included with the current aggregate in a new aggregate
            of the model; this will be taken care of outside this routine
            aftwards

      @param[in,out] model_aggregate (MinorantPointer)
          aggregate in use in the last and then the next model.

      @param[in,out] topvecs
          * on input: orthonormal basis of the collected subspace that is
            supposed to approximate the eigenspace to the largest eigenvalues,
      see Ritz_values
          * on output: same thing but maybe reduced in size

      @param[in,out] Ritz_values
          * on input: Ritz_values in cand_y for the vectors in topvecs
          * on output: same thing but maybe reduced in size as in topvecs

      @param[in,out] activedim
          * on input: the dimension of the subspace (first columns in topvecs)
             regarded as active in the last iterations of the bundle subproblem solution
    * on output: the dimension of the subspace (first columns of topvecs)
      regarded as (maybe weakly) active now
             (typically it will increase during null steps and may be tightened
             at descent steps)

      @param[out] keepsize
          columns 0..keepsize-1 of primalvecs (corresponding to the keepsize
    largest primaleigs) are included in modelvecs. The remaining
          columns need to be aggregated afterwards into the aggregate of the model

      @param[out] skippedsize
          columns activedim..activedim+skippedsize-1 of topvecs should be used
          for setting up the scaling matrix after descent steps.

      @param[in] primal_Ritzval
          the (common) Ritz value of the active subspace of the model
          (if not available use some guess like Ritz_values(0))

      @param[in] primaleigs
          eigenvalues of the last primal semdifinite model matrix
          (sorted nonincreasingly)

      @param[in] primalvecs
          corresponding (orthonormal) eigenvectors to primaleigs

      @param[in] primal_aggregate (MinorantPointer)
          aggregate in use in the last model.

      @param[in] primal_aggregate_coeff
          coefficient on how strongly the aggregated was used in the last
          primal solution to the model

      @param[in] growthrate (Real)
          factor <X,Z>/<X^-,Z^->, where X^- and Z^- are the last but one
          iterates of the interior point method

      @param[in] primalgrowth (Matrix)
          factor by which primaleigs changed in the last interior point iteration

      @param[in] dualgrowth (Matrix)
          factor by which the dual Ritz values to primalvecs changed
    during the last interior point iteration

      @param[in] cand_Ritzvec (const Matrix&)
          the (orthonormal) vectors returned by the evaluation call to the oracle

      @param[in] cand_Ritzval (const Matrix&)
    the Ritz values of cand_Ritzvec returned by the evaluation call to the oracle

      @param[in] oracle
          gives access to the evaluation oracle

      @param[in] modification_id
         the identifier for the current version of the function accounting for dynamic modifications

       @param[in] function_task
          see FunctionTask

      @param[in] function_factor
           interpreted according to function_task and the coefficients sum up to at most this value

      @param[in] model_update
          informs about whether cand_y is the result of a null_step or descent_step or a new set up.

      @param[in] center_id
          the identifier of the center point

      @param[in] center_y
          the center point

      @param[in] cand_id
          the identifier of the candidate point

      @param[in] cand_y
          the candidate (mostly differnt from the center), close to it the model should be good

      @param[in] model_maxviol
          a minorant violated by this would have caused a null step

      @param[in] diffval_center_aggregate
          difference of center value to aggregate value (nonnegative, without function_factor)

      @param[in] H
          the variable metric used in the proximal term (function_factor is already removed in this)

      @return
       - 0 on success
       - 1 on failure
    
"""
function cb_select_model!(self::CBPSCModelParameters, modelvecs::CBMatrix, model_aggregate::CBMinorantPointer, topvecs::CBMatrix, Ritz_values::CBMatrix, primal_Ritzval::Real, primaleigs::CBMatrix, primalvecs::CBMatrix, primal_aggregate::CBMinorantPointer, primal_aggregate_coeff::Real, growthrate::Real, primalgrowth::CBMatrix, dualgrowth::CBMatrix, cand_Ritzvec::CBMatrix, cand_Ritzval::CBMatrix, oracle::Union{<:CBPSCOracle,Nothing}, modification_id::Integer, function_task::CBFunctionTask, function_factor::Real, model_update::CBModelUpdate, center_id::Integer, center_y::CBMatrix, cand_id::Integer, cand_y::CBMatrix, model_maxviol::Real, diffval_center_aggregate::Real, H::CBBundleProxObject)
    skippedsize = Ref{Int}()
    keepsize = Ref{Int}()
    activedim = Ref{Int}()
    @ccall libcb.cb_pscmodelparameters_select_model(self.data::Ptr{Cvoid}, modelvecs.data::Ptr{Cvoid}, model_aggregate.data::Ptr{Cvoid}, topvecs.data::Ptr{Cvoid}, Ritz_values.data::Ptr{Cvoid}, activedim::Ref{Int}, keepsize::Ref{Int}, skippedsize::Ref{Int}, primal_Ritzval::Cdouble, primaleigs.data::Ptr{Cvoid}, primalvecs.data::Ptr{Cvoid}, primal_aggregate.data::Ptr{Cvoid}, primal_aggregate_coeff::Cdouble, growthrate::Cdouble, primalgrowth.data::Ptr{Cvoid}, dualgrowth.data::Ptr{Cvoid}, cand_Ritzvec.data::Ptr{Cvoid}, cand_Ritzval.data::Ptr{Cvoid}, (isnothing(oracle) ? C_NULL : oracle.data)::Ptr{Cvoid}, modification_id::Cint, function_task::Cint, function_factor::Cdouble, model_update.data::Ptr{Cvoid}, center_id::Cint, center_y.data::Ptr{Cvoid}, cand_id::Cint, cand_y.data::Ptr{Cvoid}, model_maxviol::Cdouble, diffval_center_aggregate::Cdouble, H.data::Ptr{Cvoid})::Cint
    return activedim[], keepsize[], skippedsize[]
end

