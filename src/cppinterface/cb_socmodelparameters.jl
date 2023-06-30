@doc raw"""
    CBSOCModelParameters(modelsize::Integer, bundlesize::Integer = 10, updaterule::Integer = 0, incr::Integer = -1)

constructor for size parameters
"""
CBSOCModelParameters(modelsize::Integer, bundlesize::Integer = 10, updaterule::Integer = 0, incr::Integer = -1) = CBSOCModelParameters(@ccall libcb.cb_socmodelparameters_new2(modelsize::Cint, bundlesize::Cint, updaterule::Cint, incr::Cint)::Ptr{Cvoid})

@doc raw"""
    CBSOCModelParameters(bp::CBBundleParameters, incr::Integer = -1)

copy constructor for BundleParameters
"""
CBSOCModelParameters(bp::CBBundleParameters, incr::Integer = -1) = CBSOCModelParameters(@ccall libcb.cb_socmodelparameters_new3(bp.data::Ptr{Cvoid}, incr::Cint)::Ptr{Cvoid})

@doc raw"""
    CBSOCModelParameters(sms::CBSOCModelParameters)

copy constructor
"""
CBSOCModelParameters(sms::CBSOCModelParameters) = CBSOCModelParameters(@ccall libcb.cb_socmodelparameters_new4(sms.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_clone_BundleParameters(self::CBSOCModelParameters)

clone
"""
cb_clone_BundleParameters(self::CBSOCModelParameters) = CBBundleParameters(@ccall libcb.cb_socmodelparameters_clone_bundleparameters(self.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_select_model!(self::CBSOCModelParameters, modelvecs::CBMatrix, aggrvec::CBMatrix, cand_SOCval::Real, cand_SOCvec::CBMatrix, center_SOCval::Real, center_SOCvec::CBMatrix, SOCvecs::CBMatrix, oracle::Union{<:CBSOCOracle,Nothing}, function_task::CBFunctionTask, function_factor::Real, model_update::CBModelUpdate, center_id::Integer, center_y::CBMatrix, cand_id::Integer, cand_y::CBMatrix, model_maxviol::Real, H::CBBundleProxObject)

* @brief SOCModel calls this for selecting the next minorants for a polyhedral model

    @param[in,out] modelvecs
        the columns span the barx subspace of the SOC face (the model is initialized if modelvecs.coldim()>0). On output it has to span at least the subspace spanned by aggrvec and cand_SOCvec

    @param[in,out] aggrvec
        current aggregate soc vector (includes all coordinates and @a function_factor according to @a function_task); if not valid it has coldim==0.

    @param[in] cand_SOCval
        lower bound on the candidate value

    @param[in] cand_SOCvec
        the SOCvector generating cand_SOCval (includes all coordinates but no function_factor)

    @param[in] center_SOCval
        lower bound on the center value

    @param[in] center_SOCvec
        the SOCvector generating center_SOCval (includes all coordinates but no function_factor)

    @param[in] SOCvecs
        collects the barx parts of the old SOCvecs

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
cb_select_model!(self::CBSOCModelParameters, modelvecs::CBMatrix, aggrvec::CBMatrix, cand_SOCval::Real, cand_SOCvec::CBMatrix, center_SOCval::Real, center_SOCvec::CBMatrix, SOCvecs::CBMatrix, oracle::Union{<:CBSOCOracle,Nothing}, function_task::CBFunctionTask, function_factor::Real, model_update::CBModelUpdate, center_id::Integer, center_y::CBMatrix, cand_id::Integer, cand_y::CBMatrix, model_maxviol::Real, H::CBBundleProxObject) = @ccall libcb.cb_socmodelparameters_select_model(self.data::Ptr{Cvoid}, modelvecs.data::Ptr{Cvoid}, aggrvec.data::Ptr{Cvoid}, cand_SOCval::Cdouble, cand_SOCvec.data::Ptr{Cvoid}, center_SOCval::Cdouble, center_SOCvec.data::Ptr{Cvoid}, SOCvecs.data::Ptr{Cvoid}, (isnothing(oracle) ? C_NULL : oracle.data)::Ptr{Cvoid}, function_task::Cint, function_factor::Cdouble, model_update.data::Ptr{Cvoid}, center_id::Cint, center_y.data::Ptr{Cvoid}, cand_id::Cint, cand_y.data::Ptr{Cvoid}, model_maxviol::Cdouble, H.data::Ptr{Cvoid})::Cint

