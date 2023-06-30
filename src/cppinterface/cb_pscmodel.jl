@doc raw"""
    cb_clear!(self::CBPSCModel)

resets all data to initial status of this class, also the bundle parameters
"""
cb_clear!(self::CBPSCModel) = @ccall libcb.cb_pscmodel_clear(self.data::Ptr{Cvoid})::Cvoid

@doc raw"""
    cb_set_variable_metric_selection!(self::CBPSCModel, vms::Union{<:CBVariableMetricSelection,Nothing} = nothing)

delete old selector and set a new one (0 is allowed resulting in no local selector); if vms casts to  PSCVariableMetricSelection*, the oracle is immediately set to this one
"""
cb_set_variable_metric_selection!(self::CBPSCModel, vms::Union{<:CBVariableMetricSelection,Nothing} = nothing) = @ccall libcb.cb_pscmodel_set_variable_metric_selection(self.data::Ptr{Cvoid}, (isnothing(vms) ? C_NULL : vms.data)::Ptr{Cvoid})::Cint

@doc raw"""
    CBPSCModel(fo::Union{<:CBPSCOracle,Nothing}, fun_factor::Real = 1., fun_task::CBFunctionTask = cbft_objective_function, cbinc::Integer = -1)

construct a model for the MatrixFunctionOracle pointed to by @a fo
"""
CBPSCModel(fo::Union{<:CBPSCOracle,Nothing}, fun_factor::Real = 1., fun_task::CBFunctionTask = cbft_objective_function, cbinc::Integer = -1) = CBPSCModel(@ccall libcb.cb_pscmodel_new((isnothing(fo) ? C_NULL : fo.data)::Ptr{Cvoid}, fun_factor::Cdouble, fun_task::Cint, cbinc::Cint)::Ptr{Cvoid})

@doc raw"""
    cb_get_oracle_object!(self::CBPSCModel)

returns the oracle
"""
cb_get_oracle_object!(self::CBPSCModel) = CBModifiableOracleObject(@ccall libcb.cb_pscmodel_get_oracle_object(self.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_lb_function!(self::CBPSCModel, y_id::Integer, y::CBMatrix)

see SumBlockModel::lb_function
"""
cb_lb_function!(self::CBPSCModel, y_id::Integer, y::CBMatrix) = @ccall libcb.cb_pscmodel_lb_function(self.data::Ptr{Cvoid}, y_id::Cint, y.data::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_get_data!(self::CBPSCModel)

see SumBlockModel::get_data
"""
cb_get_data!(self::CBPSCModel) = CBBundleData(@ccall libcb.cb_pscmodel_get_data(self.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_get_data(self::CBPSCModel)

see SumBlockModel::get_data
"""
cb_get_data(self::CBPSCModel) = CBBundleData(@ccall libcb.cb_pscmodel_get_data2(self.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_set_data!(self::CBPSCModel, bd::Union{<:CBBundleData,Nothing})

see SumBlockModel::set_data
"""
cb_set_data!(self::CBPSCModel, bd::Union{<:CBBundleData,Nothing}) = @ccall libcb.cb_pscmodel_set_data(self.data::Ptr{Cvoid}, (isnothing(bd) ? C_NULL : bd.data)::Ptr{Cvoid})::Cint

@doc raw"""
    cb_get_approximate_primal(self::CBPSCModel)

see SumBlockModel::get_approximate_primal
"""
cb_get_approximate_primal(self::CBPSCModel) = CBPrimalData(@ccall libcb.cb_pscmodel_get_approximate_primal(self.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_get_center_primal(self::CBPSCModel)

see SumBlockModel::get_center_primal
"""
cb_get_center_primal(self::CBPSCModel) = CBPrimalData(@ccall libcb.cb_pscmodel_get_center_primal(self.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_get_candidate_primal(self::CBPSCModel)

see SumBlockModel::get_candidate_primal
"""
cb_get_candidate_primal(self::CBPSCModel) = CBPrimalData(@ccall libcb.cb_pscmodel_get_candidate_primal(self.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_call_primal_extender!(self::CBPSCModel, prex::CBPrimalExtender)

see SumBlockModel::call_primal_extender
"""
cb_call_primal_extender!(self::CBPSCModel, prex::CBPrimalExtender) = @ccall libcb.cb_pscmodel_call_primal_extender(self.data::Ptr{Cvoid}, prex.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_set_bundle_parameters!(self::CBPSCModel, bp::CBBundleParameters)

see SumBlockModel::set_bundle_parameters
"""
cb_set_bundle_parameters!(self::CBPSCModel, bp::CBBundleParameters) = @ccall libcb.cb_pscmodel_set_bundle_parameters(self.data::Ptr{Cvoid}, bp.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_get_bundle_parameters(self::CBPSCModel)

see SumBlockModel::get_bundle_parameters
"""
cb_get_bundle_parameters(self::CBPSCModel) = CBBundleParameters(@ccall libcb.cb_pscmodel_get_bundle_parameters(self.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_get_ret_code(self::CBPSCModel)

see SumBlockModel::get_ret_code()
"""
cb_get_ret_code(self::CBPSCModel) = @ccall libcb.cb_pscmodel_get_ret_code(self.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_output_bundle_data(self::CBPSCModel)

this outputs the qp model data recursively for testing purposes
"""
cb_output_bundle_data(self::CBPSCModel) = @ccall libcb.cb_pscmodel_output_bundle_data(self.data::Ptr{Cvoid})::Cvoid

@doc raw"""
    cb_set_out!(self::CBPSCModel, pril::Integer = 1)

set output and outputlevel of warnings and errors recursively, see CBout
"""
cb_set_out!(self::CBPSCModel, pril::Integer = 1) = @ccall libcb.cb_pscmodel_set_out(self.data::Ptr{Cvoid}, pril::Cint)::Cvoid

