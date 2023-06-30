@doc raw"""
    CBBoxModelParameters(cbinc::Integer = -1)

default constructor
"""
CBBoxModelParameters(cbinc::Integer = -1) = CBBoxModelParameters(@ccall libcb.cb_boxmodelparameters_new(cbinc::Cint)::Ptr{Cvoid})

@doc raw"""
    CBBoxModelParameters(bp::CBBundleParameters, cbinc::Integer = -1)

copy constructor for BundleParameters
"""
CBBoxModelParameters(bp::CBBundleParameters, cbinc::Integer = -1) = CBBoxModelParameters(@ccall libcb.cb_boxmodelparameters_new2(bp.data::Ptr{Cvoid}, cbinc::Cint)::Ptr{Cvoid})

@doc raw"""
    CBBoxModelParameters(sms::CBBoxModelParameters)

copy constructor
"""
CBBoxModelParameters(sms::CBBoxModelParameters) = CBBoxModelParameters(@ccall libcb.cb_boxmodelparameters_new3(sms.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_clone_BundleParameters(self::CBBoxModelParameters)

clone
"""
cb_clone_BundleParameters(self::CBBoxModelParameters) = CBBundleParameters(@ccall libcb.cb_boxmodelparameters_clone_bundleparameters(self.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_select_model!(self::CBBoxModelParameters, box_model::CBMinorantBundle, box_coeff::CBMatrix, box_indicators::CBMatrix, box_coords::CBIndexmatrix, box_complvalues::CBMatrix, nnc_model::CBMinorantBundle, nnc_coeff::CBMatrix, nnc_indicators::CBMatrix, coord_switching::CBMatrix, minorants::CBMinorantBundle, cand_minorant::CBMinorantPointer, cand_boxvec::CBPrimalMatrix, aggr_boxvec::CBPrimalMatrix, aggr_scaleval::Real, oracle::Union{<:CBBoxOracle,Nothing}, modification_id::Integer, function_task::CBFunctionTask, function_factor::Real, model_update::CBModelUpdate, center_id::Integer, center_y::CBMatrix, cand_id::Integer, cand_y::CBMatrix, model_maxviol::Real, H::CBBundleProxObject)

* @brief BoxModel calls this for selecting the next coordinates
         for a specialized polyhedral model with a box part and an
         nnc part for aggregates, see the general explanation of the
         class.

         There is little experience on how to do this.
         The current routine is minimalistic and simply uses
         the weighted number of switches in each coordinate
         between lower and upper bounds to select coordinates.

         To allow maybe better choices in future implementations
         the arguments try to pass all potentially relevant
         information items.

         On putput the coefficient values of the new model must be
         feasible and have to generate the same aggregate (the aggregate
         is maintained), and the new candidate minorant must be in the
         feasible set.

      @param[in,out] box_model
         the boxmodel holds the minorants describing the BoxBlock part of
         the model for selected coordinates and, unless exact or empty,
         in the last position the complement coordinates of a feasible
         point (e.g. aggr_boxvec)

      @param[in,out] box_coeff
         - on input coefficients of BoxBlock determined in last
           BundleMethod::eval_augmodel giving rise (together with nnc_coeff)
           to the aggrgate,
         - on output they match the new model and still give rise to
           the same aggregate

      @param[in,out] box_indicators
         indicators for activity of box minorants, indicators may but need
         not be maintained

      @param[in,out] box_coords
         the coordinates selected to have their respective interval range
         in the model

      @param[in,out] box_complvalues
         a point with 0 in the box_coords and feasible
         coordinate values in the complement (if not empty)

      @param[in,out] nnc_model
         if box_model is empty, this spans at least the aggregate (if available)
         and the candidate (always); if box_model is not empty but not
         the entire box, nnc_model typically holds one of the candidate or
         the aggregate; if box_model is the entire box, nnc_model is empty

      @param[in,out] nnc_coeff
         - on input coefficients of BoxBlock determined in last
           BundleMethod::eval_augmodel giving rise (together with box_coeff)
           to the aggrgate,
         - on output they match the new model and still give rise to
           the same aggregate

      @param[in,out] nnc_indicators
         indicators for activity of minorants, indicators may but need
         not be maintained

      @param[in] coord_switching
         keeps track of which coordinates where changing the most in the past
         by forming a weighted average in BoxModel::eval_function()

      @param[in] minorants
          the vector of MinorantPointer gives additional minorants
    collected over time (some may be duplicates also of those in
          nnc_model)

     @param[in] cand_minorant
         the (eps)sugradient linear minorant returned
         by BoxModel::eval_function for the candidate (without function factor)

     @param[in] cand_boxvec
         the maximizer over the box for the current candidate

     @param[in] aggr_boxvec
         the primal aggregate vector in the box (without function_factor)
         giving rise to the aggregate; not initialized if zerodimensional

     @param[in] aggr_scaleval
         0<= aggr_scalevale <= function_factor, ==function_factor if
         function_task==Objective_Function; the aggregate with
         function_factor is aggr_boxvec*aggr_scaleval;

      @param[in] oracle
          gives access to lower and upper bounds of the box

      @param[in] modification_id
          the identifier of the current function version to be used in generating
          specialized minorants corresponding to the coordinate vectors

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
cb_select_model!(self::CBBoxModelParameters, box_model::CBMinorantBundle, box_coeff::CBMatrix, box_indicators::CBMatrix, box_coords::CBIndexmatrix, box_complvalues::CBMatrix, nnc_model::CBMinorantBundle, nnc_coeff::CBMatrix, nnc_indicators::CBMatrix, coord_switching::CBMatrix, minorants::CBMinorantBundle, cand_minorant::CBMinorantPointer, cand_boxvec::CBPrimalMatrix, aggr_boxvec::CBPrimalMatrix, aggr_scaleval::Real, oracle::Union{<:CBBoxOracle,Nothing}, modification_id::Integer, function_task::CBFunctionTask, function_factor::Real, model_update::CBModelUpdate, center_id::Integer, center_y::CBMatrix, cand_id::Integer, cand_y::CBMatrix, model_maxviol::Real, H::CBBundleProxObject) = @ccall libcb.cb_boxmodelparameters_select_model(self.data::Ptr{Cvoid}, box_model.data::Ptr{Cvoid}, box_coeff.data::Ptr{Cvoid}, box_indicators.data::Ptr{Cvoid}, box_coords.data::Ptr{Cvoid}, box_complvalues.data::Ptr{Cvoid}, nnc_model.data::Ptr{Cvoid}, nnc_coeff.data::Ptr{Cvoid}, nnc_indicators.data::Ptr{Cvoid}, coord_switching.data::Ptr{Cvoid}, minorants.data::Ptr{Cvoid}, cand_minorant.data::Ptr{Cvoid}, cand_boxvec.data::Ptr{Cvoid}, aggr_boxvec.data::Ptr{Cvoid}, aggr_scaleval::Cdouble, (isnothing(oracle) ? C_NULL : oracle.data)::Ptr{Cvoid}, modification_id::Cint, function_task::Cint, function_factor::Cdouble, model_update.data::Ptr{Cvoid}, center_id::Cint, center_y.data::Ptr{Cvoid}, cand_id::Cint, cand_y.data::Ptr{Cvoid}, model_maxviol::Cdouble, H.data::Ptr{Cvoid})::Cint

