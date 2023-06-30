@doc raw"""
    CBVariableMetricSVDSelection(cbincr::Integer = -1)

default constructor
"""
CBVariableMetricSVDSelection(cbincr::Integer = -1) = CBVariableMetricSVDSelection(@ccall libcb.cb_variablemetricsvdselection_new(cbincr::Cint)::Ptr{Cvoid})

@doc raw"""
    CBVariableMetricSVDSelection(in_n_latest_minorants::Integer, in_selection_method::Integer, in_oldfactor::Real = 0., cbincr::Integer = -1)

constructor for specifying values for n_latest_minorants and selection_method
"""
CBVariableMetricSVDSelection(in_n_latest_minorants::Integer, in_selection_method::Integer, in_oldfactor::Real = 0., cbincr::Integer = -1) = CBVariableMetricSVDSelection(@ccall libcb.cb_variablemetricsvdselection_new2(in_n_latest_minorants::Cint, in_selection_method::Cint, in_oldfactor::Cdouble, cbincr::Cint)::Ptr{Cvoid})

@doc raw"""
    cb_get_UVlambda!(self::CBVariableMetricSVDSelection, U::CBMatrix, V::CBMatrix, lam::CBMatrix, cand::CBMatrix)

for current ongoing experiments with variable metric routines
"""
cb_get_UVlambda!(self::CBVariableMetricSVDSelection, U::CBMatrix, V::CBMatrix, lam::CBMatrix, cand::CBMatrix) = @ccall libcb.cb_variablemetricsvdselection_get_uvlambda(self.data::Ptr{Cvoid}, U.data::Ptr{Cvoid}, V.data::Ptr{Cvoid}, lam.data::Ptr{Cvoid}, cand.data::Ptr{Cvoid})::Cvoid

@doc raw"""
    cb_get_n_latest_minorants(self::CBVariableMetricSVDSelection)

returns n_latest_minorants
"""
cb_get_n_latest_minorants(self::CBVariableMetricSVDSelection) = @ccall libcb.cb_variablemetricsvdselection_get_n_latest_minorants(self.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_set_n_latest_minorants!(self::CBVariableMetricSVDSelection, nlm::Integer)

sets n_latest_minorants
"""
cb_set_n_latest_minorants!(self::CBVariableMetricSVDSelection, nlm::Integer) = @ccall libcb.cb_variablemetricsvdselection_set_n_latest_minorants(self.data::Ptr{Cvoid}, nlm::Cint)::Cvoid

@doc raw"""
    cb_get_selection_method(self::CBVariableMetricSVDSelection)

returns selection_method
"""
cb_get_selection_method(self::CBVariableMetricSVDSelection) = @ccall libcb.cb_variablemetricsvdselection_get_selection_method(self.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_set_selection_method!(self::CBVariableMetricSVDSelection, sm::Integer)

sets selection_method
"""
cb_set_selection_method!(self::CBVariableMetricSVDSelection, sm::Integer) = @ccall libcb.cb_variablemetricsvdselection_set_selection_method(self.data::Ptr{Cvoid}, sm::Cint)::Cvoid

@doc raw"""
    cb_add_variable_metric!(self::CBVariableMetricSVDSelection, H::CBVariableMetric, y_id::Integer, y::CBMatrix, descent_step::Bool, weightu::Real, model_maxviol::Real, indices::Union{<:CBIndexmatrix,Nothing}, bundle_data::CBVariableMetricBundleData)

see ConicBundle::VariableMetricSelection::add_variable_metric()
"""
cb_add_variable_metric!(self::CBVariableMetricSVDSelection, H::CBVariableMetric, y_id::Integer, y::CBMatrix, descent_step::Bool, weightu::Real, model_maxviol::Real, indices::Union{<:CBIndexmatrix,Nothing}, bundle_data::CBVariableMetricBundleData) = @ccall libcb.cb_variablemetricsvdselection_add_variable_metric(self.data::Ptr{Cvoid}, H.data::Ptr{Cvoid}, y_id::Cint, y.data::Ptr{Cvoid}, descent_step::Cint, weightu::Cdouble, model_maxviol::Cdouble, (isnothing(indices) ? C_NULL : indices.data)::Ptr{Cvoid}, bundle_data.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_clone_VariableMetricSelection!(self::CBVariableMetricSVDSelection)

clone: the values are only preserved for those contained in the constructor: n_latest_minorants, selection_method and oldfactor
"""
cb_clone_VariableMetricSelection!(self::CBVariableMetricSVDSelection) = CBVariableMetricSelection(@ccall libcb.cb_variablemetricsvdselection_clone_variablemetricselection(self.data::Ptr{Cvoid})::Ptr{Cvoid})

