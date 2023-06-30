/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/FunctionObjectModification.hxx
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




#ifndef CONICBUNDLE_FUNCTIONOBJECTMODIFICATION_HXX
#define CONICBUNDLE_FUNCTIONOBJECTMODIFICATION_HXX


/**  @file FunctionObjectModification.hxx
    @brief Header declaring the classes ConicBundle::FunctionObjectModification and ConicBundle::FunObjModMap
    @version 1.0
    @date 2014-11-12
    @author Christoph Helmberg
*/

#include <map>
#include "AFTModification.hxx"

namespace ConicBundle {

  /** @ingroup dynamic_modification_support
  */
  //@{

  /** @brief abstract interface for informing the bundle model routines
      about changes in the oracle function and, possibly, for modifying
      the affine transformation of the oracle function's arguments.

      This is just a container for the pointers to the
      (one or two) modification(s). Nothing is deleted here, but
      it may well be deleted by the function that receives the container ...

      Dynamic changes in the oracle function (during the runtime of the
      bundle solver) are a subtle issue and should only be used if you are
      aware of its mathematical implications!


  */



  class FunctionObjectModification {
  private:
    /// pointer to a description of modifications of the oracle (and this typically changes the model as well, but not necessarily its tranformation), may be NULL. Any object pointed to here is NOT deleted on destruction of this.
    OracleModification* oracle_modification;
    /// pointer to a description of the modifications of the AffineFunctionTransformation of the model (does not change the model, but the transformation), may be NULL. Any object pointed to here is NOT deleted on destruction of this
    AFTModification* trafo_aftmdf;

  public:
    /// specifies Modifications for oracle and possibly its AFT
    FunctionObjectModification(OracleModification* in_oracle_md = 0,
      AFTModification* in_trafo_aftm = 0) :
      oracle_modification(in_oracle_md), trafo_aftmdf(in_trafo_aftm) {
    }

    /// copy
    FunctionObjectModification(const FunctionObjectModification& fom) :
      oracle_modification(fom.oracle_modification), trafo_aftmdf(fom.trafo_aftmdf) {
    }

    /// no deletions are carried out here
    virtual ~FunctionObjectModification() {
    }

    /// assignment operator (copies the pointers, no new objects)
    FunctionObjectModification& operator=(const FunctionObjectModification& fom) {
      oracle_modification = fom.oracle_modification;
      trafo_aftmdf = fom.trafo_aftmdf;
      return *this;
    }

    /// set an AFTModification for the AffineFunctionTransformation of the argument 
    virtual void set_oracle_modification(OracleModification* oracle_mdf) {
      oracle_modification = oracle_mdf;
    }

    /// set an AFTModification for the AffineFunctionTransformation of the argument 
    virtual void set_aft_modification(AFTModification* in_aftmdf) {
      trafo_aftmdf = in_aftmdf;
    }

    /// returns the current AFTModification 
    virtual AFTModification* get_aft_modification() {
      return trafo_aftmdf;
    }

    /// returns the current AFTModification for the const case 
    virtual const AFTModification* get_aft_modification() const {
      return trafo_aftmdf;
    }

    /// returns the current AFTModification 
    virtual OracleModification* get_oracle_modification() {
      return oracle_modification;
    }

    /// returns the current AFTModification for the const case
    virtual const OracleModification* get_oracle_modification() const {
      return oracle_modification;
    }

  };


  /** @brief In order to pass on problem modifications, FunObjModMap allows to specify a FunctionObjectModification for each FunctionObject whose changes are not fully described by the default changes induced by those on the Groundset.

     In BundleModel::apply_modification each BundleModel will check whether the FunctionObject representing its oracle is present in the map. If so, it carries out the modifications described here, or it resorts to the defaul actions based on the GroundsetModifications of the Groundset.
  */
  typedef  std::map<const FunctionObject*, FunctionObjectModification> FunObjModMap;



  //@}

}

#endif

