@doc raw"""
    CBSumBundleParameters(modelsize::Integer, bundlesize::Integer = 10, updaterule::Integer = 0, incr::Integer = -1)

constructor for customized parameters
"""
CBSumBundleParameters(modelsize::Integer, bundlesize::Integer = 10, updaterule::Integer = 0, incr::Integer = -1) = CBSumBundleParameters(@ccall libcb.cb_sumbundleparameters_new2(modelsize::Cint, bundlesize::Cint, updaterule::Cint, incr::Cint)::Ptr{Cvoid})

@doc raw"""
    CBSumBundleParameters(sbp::CBSumBundleParameters)

copy constructor
"""
CBSumBundleParameters(sbp::CBSumBundleParameters) = CBSumBundleParameters(@ccall libcb.cb_sumbundleparameters_new3(sbp.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_clone_BundleParameters(self::CBSumBundleParameters)

clone
"""
cb_clone_BundleParameters(self::CBSumBundleParameters) = CBBundleParameters(@ccall libcb.cb_sumbundleparameters_clone_bundleparameters(self.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_select_model!(self::CBSumBundleParameters, model_indices::CBIndexmatrix, cand_id::Integer, cand_y::CBMatrix, minorants::CBMinorantBundle, aggr_index::Integer, model_maxviol::Real, H::CBBundleProxObject, model_update::CBModelUpdate)

* @brief FunctionModel and SumBundleHandler call this for selecting the next minorants for a polyhedral model

      @param[out] model_indices
          the indices of minorants selected for the model; index 0 is always the aggregate indicated by the input index aggr_index

      @param[in] cand_id
          the identifier of the candidate point supplied next

      @param[in] cand_y
          the candidate (differnte from the center), close to it the model should be good

      @param[in] minorants
          the vector of MinorantPointer gives the minorants
    out of which the model should be selected.

      @param[in] aggr_index
          the index of the aggregate within the minorants

      @param[in] model_maxviol
          a minorant violated by this would have caused a null step

      @param[in] H
          the proximal term used for determining the given cand_y

      @param[in] model_update
          informs about whether cand_y is the result of a null_step or descent_step or aomw other a new set up.

      @return
       - 0 on success
       - 1 on failure
    
"""
cb_select_model!(self::CBSumBundleParameters, model_indices::CBIndexmatrix, cand_id::Integer, cand_y::CBMatrix, minorants::CBMinorantBundle, aggr_index::Integer, model_maxviol::Real, H::CBBundleProxObject, model_update::CBModelUpdate) = @ccall libcb.cb_sumbundleparameters_select_model(self.data::Ptr{Cvoid}, model_indices.data::Ptr{Cvoid}, cand_id::Cint, cand_y.data::Ptr{Cvoid}, minorants.data::Ptr{Cvoid}, aggr_index::Cint, model_maxviol::Cdouble, H.data::Ptr{Cvoid}, model_update.data::Ptr{Cvoid})::Cint

