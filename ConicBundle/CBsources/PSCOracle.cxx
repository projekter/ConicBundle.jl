/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/PSCOracle.cxx
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



#include "PSCOracle.hxx"
#include "matop.hxx"

using namespace CH_Matrix_Classes; 

namespace ConicBundle {

  PSCBundleParameters::~PSCBundleParameters(){}
  
  int PSCBundleParameters::init(const BundleParameters& bp)
    {
      int err=BundleParameters::init(bp);
      const PSCBundleParameters* p=dynamic_cast<const PSCBundleParameters*>(&bp);
      if (p){
	psc_model_size=(p->psc_model_size<1)?1:p->psc_model_size;
	psc_bundle_size=(p->psc_bundle_size<psc_model_size)?psc_model_size:p->psc_bundle_size;
	psc_new_subgradients=(p->psc_new_subgradients<1)?1:p->psc_new_subgradients;
	psc_keep=(p->psc_keep<0)?0:p->psc_keep;
	psc_aggregates=(p->psc_aggregates<1)?1:p->psc_aggregates;
	psc_tolerance=(p->psc_tolerance<eps_Real)?eps_Real:p->psc_tolerance;
	psc_update_rule=p->psc_update_rule;
      }
      return err;
    }

}
