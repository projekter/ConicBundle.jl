/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/SumBundleParametersObject.cxx
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



#include "SumBundleParametersObject.hxx"

using namespace CH_Matrix_Classes;

namespace ConicBundle {

  int SumBundleParametersObject::init(const BundleParameters& bp)
  {
    int err=BundleParameters::init(bp);
    const SumBundleParametersObject* p=dynamic_cast<const SumBundleParametersObject*>(&bp);
    if (p){
      acceptable_mode=p->acceptable_mode;
      delete vm_selection;
      vm_selection=0;
      if (p->vm_selection)
	vm_selection=p->vm_selection->clone_VariableMetricSelection();
      else
	vm_selection=0;
    }
    return err;
  }


  SumBundleParametersObject::~SumBundleParametersObject()
  {delete vm_selection;}


} // end namespace ConicBundle
