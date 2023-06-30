/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/SumModelParameters.hxx
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



#ifndef CONICBUNDLE_SUMMODELPARAMETERS_HXX
#define CONICBUNDLE_SUMMODELPARAMETERS_HXX

/**  @file SumModelParameters.hxx
    @brief Header declaring the class ConicBundle::SumModelParameters
    @version 1.0
    @date 2014-07-24
    @author Christoph Helmberg
*/

#include "SumModelParametersObject.hxx"
#include "SumBundleParametersObject.hxx"

namespace ConicBundle {
  /** @ingroup InternalBundleModel

  */
  //@{

  /** @brief specifies the suggestion procedure for including models in a joint SumBundle model of SumModel


  */

  class SumModelParameters : public SumModelParametersObject {
  private:
  public:

    /// default constructor
    SumModelParameters(CBout* cb = 0, int cbinc = -1) :
      SumModelParametersObject(cb, cbinc) {
    }

    /// constructor setting the SumBundleParameters to use for passing on; the BundleParameters also get this value
    SumModelParameters(const SumBundleParametersObject& sbp,
      int in_max_local_models = 0) :
      CBout(sbp), SumModelParametersObject() {
      init(sbp);
      max_local_models = in_max_local_models;
    }

    ///copy constructor
    SumModelParameters(const SumModelParameters& smp) :
      CBout(smp), SumModelParametersObject() {
      init(smp);
    }

    ///virtual destructor
    virtual ~SumModelParameters();

    /// return a new clone object of this on the heap (caller needs to delete the result)
    virtual BundleParameters* clone_BundleParameters() const {
      return new SumModelParameters(*this);
    }

    /// generate a suggestion list on which to include in the SumBundle (the othersi should keep their local models) and what kind of (Sum)BundleParameters to use. bp must be NULL on input; if NULL on ouput there is no preference; if not NULL the object pointed to stays property of *this
    virtual int select_models(SumModel::ModelMap& modelmap);

  };


  //@}

}

#endif

