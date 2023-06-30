@doc raw"""
    CBPSCVariableMetricSelection(cbincr::Integer = -1)

default constructor
"""
CBPSCVariableMetricSelection(cbincr::Integer = -1) = CBPSCVariableMetricSelection(@ccall libcb.cb_pscvariablemetricselection_new(cbincr::Cint)::Ptr{Cvoid})

@doc raw"""
    cb_get_selection_method(self::CBPSCVariableMetricSelection)

returns selection_method
"""
cb_get_selection_method(self::CBPSCVariableMetricSelection) = @ccall libcb.cb_pscvariablemetricselection_get_selection_method(self.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_set_selection_method!(self::CBPSCVariableMetricSelection, sm::Integer)

sets selection_method
"""
cb_set_selection_method!(self::CBPSCVariableMetricSelection, sm::Integer) = @ccall libcb.cb_pscvariablemetricselection_set_selection_method(self.data::Ptr{Cvoid}, sm::Cint)::Cvoid

@doc raw"""
    cb_get_oldfactor(self::CBPSCVariableMetricSelection)

returns the parameter
"""
cb_get_oldfactor(self::CBPSCVariableMetricSelection) = @ccall libcb.cb_pscvariablemetricselection_get_oldfactor(self.data::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_set_oldfactor!(self::CBPSCVariableMetricSelection, of::Real)

sets the parameter
"""
cb_set_oldfactor!(self::CBPSCVariableMetricSelection, of::Real) = @ccall libcb.cb_pscvariablemetricselection_set_oldfactor(self.data::Ptr{Cvoid}, of::Cdouble)::Cvoid

@doc raw"""
    cb_get_maxeigval_factor(self::CBPSCVariableMetricSelection)

returns the parameter
"""
cb_get_maxeigval_factor(self::CBPSCVariableMetricSelection) = @ccall libcb.cb_pscvariablemetricselection_get_maxeigval_factor(self.data::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_set_maxeigval_factor!(self::CBPSCVariableMetricSelection, ef::Real)

sets the parameter
"""
cb_set_maxeigval_factor!(self::CBPSCVariableMetricSelection, ef::Real) = @ccall libcb.cb_pscvariablemetricselection_set_maxeigval_factor(self.data::Ptr{Cvoid}, ef::Cdouble)::Cvoid

@doc raw"""
    cb_get_mineigval_factor(self::CBPSCVariableMetricSelection)

returns the parameter
"""
cb_get_mineigval_factor(self::CBPSCVariableMetricSelection) = @ccall libcb.cb_pscvariablemetricselection_get_mineigval_factor(self.data::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_set_mineigval_factor!(self::CBPSCVariableMetricSelection, ef::Real)

sets the parameter
"""
cb_set_mineigval_factor!(self::CBPSCVariableMetricSelection, ef::Real) = @ccall libcb.cb_pscvariablemetricselection_set_mineigval_factor(self.data::Ptr{Cvoid}, ef::Cdouble)::Cvoid

@doc raw"""
    cb_set_oracle!(self::CBPSCVariableMetricSelection, psco::Union{<:CBPSCOracle,Nothing})

sets the oracle pointer to this value (NULL is allowed, but calling add_variable_metric() then results in a WARNING and an error is returned); this is called by PSCModel when installing this
"""
cb_set_oracle!(self::CBPSCVariableMetricSelection, psco::Union{<:CBPSCOracle,Nothing}) = @ccall libcb.cb_pscvariablemetricselection_set_oracle(self.data::Ptr{Cvoid}, (isnothing(psco) ? C_NULL : psco.data)::Ptr{Cvoid})::Cvoid

@doc raw"""
    cb_add_variable_metric!(self::CBPSCVariableMetricSelection, H::CBVariableMetric, y_id::Integer, y::CBMatrix, descent_step::Bool, weightu::Real, model_maxviol::Real, indices::Union{<:CBIndexmatrix,Nothing}, bundle_data::CBVariableMetricBundleData)

see ConicBundle::VariableMetricSelection::add_variable_metric(); here it must be possible to cast bundle_data to PSCData&, otherwise the routine returns an error. It only adds something after a descent_step
"""
cb_add_variable_metric!(self::CBPSCVariableMetricSelection, H::CBVariableMetric, y_id::Integer, y::CBMatrix, descent_step::Bool, weightu::Real, model_maxviol::Real, indices::Union{<:CBIndexmatrix,Nothing}, bundle_data::CBVariableMetricBundleData) = @ccall libcb.cb_pscvariablemetricselection_add_variable_metric(self.data::Ptr{Cvoid}, H.data::Ptr{Cvoid}, y_id::Cint, y.data::Ptr{Cvoid}, descent_step::Cint, weightu::Cdouble, model_maxviol::Cdouble, (isnothing(indices) ? C_NULL : indices.data)::Ptr{Cvoid}, bundle_data.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_clone_VariableMetricSelection!(self::CBPSCVariableMetricSelection)

clone: the values are only preserved for those contained in the constructor: n_latest_minorants, selection_method and oldfactor
"""
cb_clone_VariableMetricSelection!(self::CBPSCVariableMetricSelection) = CBVariableMetricSelection(@ccall libcb.cb_pscvariablemetricselection_clone_variablemetricselection(self.data::Ptr{Cvoid})::Ptr{Cvoid})

