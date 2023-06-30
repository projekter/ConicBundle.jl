/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/VariableMetric.cxx
    This file is part of ConciBundle, a C/C++ library for convex optimization.

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.

***************************************************************************** */



#include "VariableMetric.hxx"
#include "AffineFunctionTransformation.hxx"

using namespace CH_Matrix_Classes;

namespace ConicBundle {

  VariableMetricModel::VariableMetricModel(CBout* cb, int cbinc) :CBout(cb, cbinc) {
  }

  VariableMetricModel::~VariableMetricModel() {
  }

  int VariableMetricModel::add_variable_metric(VariableMetric& /* H */,
    CH_Matrix_Classes::Integer /* y_id */,
    const CH_Matrix_Classes::Matrix& /* y */,
    bool /* descent_step */,
    CH_Matrix_Classes::Real /* weightu */,
    CH_Matrix_Classes::Real /* model_maxviol */,
    const CH_Matrix_Classes::Indexmatrix* /* indices */) {
    return 0;
  }

  VariableMetricBundleData::~VariableMetricBundleData() {
  }

  VariableMetricSelection::VariableMetricSelection(CBout* cb, int cbincr) :CBout(cb, cbincr) {
  }

  VariableMetricSelection::~VariableMetricSelection() {
  }

  int VariableMetricSelection::add_variable_metric(VariableMetric& /* H */,
    CH_Matrix_Classes::Integer /* y_id */,
    const CH_Matrix_Classes::Matrix& /* y */,
    bool /* descent_step */,
    CH_Matrix_Classes::Real /* weightu */,
    CH_Matrix_Classes::Real /* model_maxviol */,
    const CH_Matrix_Classes::Indexmatrix* /* indices */,
    VariableMetricBundleData& /* bundle_data */) {
    return 0;
  }


  VariableMetric::~VariableMetric() {
    delete vm_selection;
  }


}

