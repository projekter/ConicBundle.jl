@doc raw"""
    cb_clear!(self::CBBoxModel)

resets all data to initial status of this class, also the bundle parameters
"""
cb_clear!(self::CBBoxModel) = @ccall libcb.cb_boxmodel_clear(self.data::Ptr{Cvoid})::Cvoid

@doc raw"""
    CBBoxModel(fo::Union{<:CBBoxOracle,Nothing}, fun_factor::Real = 1., fun_task::CBFunctionTask = cbft_objective_function, cbinc::Integer = -1)

construct a model for the MatrixBoxOracle pointed to by @a fo
"""
CBBoxModel(fo::Union{<:CBBoxOracle,Nothing}, fun_factor::Real = 1., fun_task::CBFunctionTask = cbft_objective_function, cbinc::Integer = -1) = CBBoxModel(@ccall libcb.cb_boxmodel_new((isnothing(fo) ? C_NULL : fo.data)::Ptr{Cvoid}, fun_factor::Cdouble, fun_task::Cint, cbinc::Cint)::Ptr{Cvoid})

@doc raw"""
    cb_get_oracle_object!(self::CBBoxModel)

returns the oracle
"""
cb_get_oracle_object!(self::CBBoxModel) = CBModifiableOracleObject(@ccall libcb.cb_boxmodel_get_oracle_object(self.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_lb_function!(self::CBBoxModel, y_id::Integer, y::CBMatrix)

see SumBlockModel::lb_function
"""
cb_lb_function!(self::CBBoxModel, y_id::Integer, y::CBMatrix) = @ccall libcb.cb_boxmodel_lb_function(self.data::Ptr{Cvoid}, y_id::Cint, y.data::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_get_data!(self::CBBoxModel)

see SumBlockModel::get_data
"""
cb_get_data!(self::CBBoxModel) = CBBundleData(@ccall libcb.cb_boxmodel_get_data(self.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_get_data(self::CBBoxModel)

see SumBlockModel::get_data
"""
cb_get_data(self::CBBoxModel) = CBBundleData(@ccall libcb.cb_boxmodel_get_data2(self.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_set_data!(self::CBBoxModel, bd::Union{<:CBBundleData,Nothing})

see SumBlockModel::set_data
"""
cb_set_data!(self::CBBoxModel, bd::Union{<:CBBundleData,Nothing}) = @ccall libcb.cb_boxmodel_set_data(self.data::Ptr{Cvoid}, (isnothing(bd) ? C_NULL : bd.data)::Ptr{Cvoid})::Cint

@doc raw"""
    cb_get_approximate_primal(self::CBBoxModel)

see SumBlockModel::get_approximate_primal
"""
cb_get_approximate_primal(self::CBBoxModel) = CBPrimalData(@ccall libcb.cb_boxmodel_get_approximate_primal(self.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_get_center_primal(self::CBBoxModel)

see SumBlockModel::get_cneter_primal
"""
cb_get_center_primal(self::CBBoxModel) = CBPrimalData(@ccall libcb.cb_boxmodel_get_center_primal(self.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_get_candidate_primal(self::CBBoxModel)

see SumBlockModel::get_candidate_primal
"""
cb_get_candidate_primal(self::CBBoxModel) = CBPrimalData(@ccall libcb.cb_boxmodel_get_candidate_primal(self.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_call_primal_extender!(self::CBBoxModel, prex::CBPrimalExtender)

see SumBlockModel::call_primal_extender
"""
cb_call_primal_extender!(self::CBBoxModel, prex::CBPrimalExtender) = @ccall libcb.cb_boxmodel_call_primal_extender(self.data::Ptr{Cvoid}, prex.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_set_bundle_parameters!(self::CBBoxModel, bp::CBBundleParameters)

see SumBlockModel::set_bundle_parameters
"""
cb_set_bundle_parameters!(self::CBBoxModel, bp::CBBundleParameters) = @ccall libcb.cb_boxmodel_set_bundle_parameters(self.data::Ptr{Cvoid}, bp.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_get_bundle_parameters(self::CBBoxModel)

see SumBlockModel::get_bundle_parameters
"""
cb_get_bundle_parameters(self::CBBoxModel) = CBBundleParameters(@ccall libcb.cb_boxmodel_get_bundle_parameters(self.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_get_ret_code(self::CBBoxModel)

see SumBlockModel::get_ret_code()
"""
cb_get_ret_code(self::CBBoxModel) = @ccall libcb.cb_boxmodel_get_ret_code(self.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_set_out!(self::CBBoxModel, pril::Integer = 1)

set output and outputlevel of warnings and errors recursively, see CBout
"""
cb_set_out!(self::CBBoxModel, pril::Integer = 1) = @ccall libcb.cb_boxmodel_set_out(self.data::Ptr{Cvoid}, pril::Cint)::Cvoid

