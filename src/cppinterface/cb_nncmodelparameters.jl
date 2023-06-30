@doc raw"""
    CBNNCModelParameters(modelsize::Integer, bundlesize::Integer = 10, updaterule::Integer = 0, incr::Integer = -1)

constructor for size parameters
"""
CBNNCModelParameters(modelsize::Integer, bundlesize::Integer = 10, updaterule::Integer = 0, incr::Integer = -1) = CBNNCModelParameters(@ccall libcb.cb_nncmodelparameters_new2(modelsize::Cint, bundlesize::Cint, updaterule::Cint, incr::Cint)::Ptr{Cvoid})

@doc raw"""
    CBNNCModelParameters(bp::CBBundleParameters, incr::Integer = -1)

copy constructor for BundleParameters
"""
CBNNCModelParameters(bp::CBBundleParameters, incr::Integer = -1) = CBNNCModelParameters(@ccall libcb.cb_nncmodelparameters_new3(bp.data::Ptr{Cvoid}, incr::Cint)::Ptr{Cvoid})

@doc raw"""
    CBNNCModelParameters(sms::CBNNCModelParameters)

copy constructor
"""
CBNNCModelParameters(sms::CBNNCModelParameters) = CBNNCModelParameters(@ccall libcb.cb_nncmodelparameters_new4(sms.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_clone_BundleParameters(self::CBNNCModelParameters)

clone
"""
cb_clone_BundleParameters(self::CBNNCModelParameters) = CBBundleParameters(@ccall libcb.cb_nncmodelparameters_clone_bundleparameters(self.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_select_model!(self::CBNNCModelParameters, model::CBMinorantBundle, coefficients::CBMatrix, activity_indicators::CBMatrix, aggregate::CBMinorantPointer, center_minorant::CBMinorantPointer, cand_minorants::CBMinorantBundle, old_minorants::CBMinorantBundle, oracle::Union{<:CBMatrixFunctionOracle,Nothing}, function_task::CBFunctionTask, function_factor::Real, model_update::CBModelUpdate, center_id::Integer, center_y::CBMatrix, cand_id::Integer, cand_y::CBMatrix, model_maxviol::Real, H::CBBundleProxObject)

* @brief NNCModel calls this for selecting the next minorants for a polyhedral model

      @param[in,out] model
          contains the model of the previous subproblem or a minimal starting choice

      @param[in,out]  coefficients
          contains the coefficients resulting from the last bundel subproblem (on input and on output they should generate the aggregate, if the aggregate is valid)

      @param[in,out]  activity_indicators
          contains activity_indicators resulting from the last bundel subproblem (currently 1 if considered active and 0 otherwise). Those indicated active in the bundle on input will also be indicated active

      @param[in] aggregate
          if valid, it holds the current aggregate. It arises as nonnegative combinations of the model by the coefficients (so it inculde the current function_factor). This property of model and coefficients should be maintained on output

      @param[in] center_minorant
          if valid, it holds the current aggregate. It arises as nonnegative combinations of the model by the coefficients, and that should be maintained on output

      @param[in] cand_minorants
          holds the vector of minorants returned by evaluation in the candidate point

      @param[in] old_minorants
          the vector of MinorantPointer gives additional minorants
    collected over time (some may be duplicates, some are
    most certainly already contained in the bundle on input)

      @param[in] oracle
          gives access to the evaluation oracle

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

      @param[in] H
          the variable metric used in the proximal term (function_factor is already removed in this)

      @return
       - 0 on success
       - 1 on failure
    
"""
cb_select_model!(self::CBNNCModelParameters, model::CBMinorantBundle, coefficients::CBMatrix, activity_indicators::CBMatrix, aggregate::CBMinorantPointer, center_minorant::CBMinorantPointer, cand_minorants::CBMinorantBundle, old_minorants::CBMinorantBundle, oracle::Union{<:CBMatrixFunctionOracle,Nothing}, function_task::CBFunctionTask, function_factor::Real, model_update::CBModelUpdate, center_id::Integer, center_y::CBMatrix, cand_id::Integer, cand_y::CBMatrix, model_maxviol::Real, H::CBBundleProxObject) = @ccall libcb.cb_nncmodelparameters_select_model(self.data::Ptr{Cvoid}, model.data::Ptr{Cvoid}, coefficients.data::Ptr{Cvoid}, activity_indicators.data::Ptr{Cvoid}, aggregate.data::Ptr{Cvoid}, center_minorant.data::Ptr{Cvoid}, cand_minorants.data::Ptr{Cvoid}, old_minorants.data::Ptr{Cvoid}, (isnothing(oracle) ? C_NULL : oracle.data)::Ptr{Cvoid}, function_task::Cint, function_factor::Cdouble, model_update.data::Ptr{Cvoid}, center_id::Cint, center_y.data::Ptr{Cvoid}, cand_id::Cint, cand_y.data::Ptr{Cvoid}, model_maxviol::Cdouble, H.data::Ptr{Cvoid})::Cint

