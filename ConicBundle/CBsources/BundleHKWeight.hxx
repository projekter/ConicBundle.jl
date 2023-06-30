/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/BundleHKWeight.hxx
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



#ifndef CONICBUNDLE_BUNDLEHKWEIGHT_HXX
#define CONICBUNDLE_BUNDLEHKWEIGHT_HXX

#include "BundleWeight.hxx"

namespace ConicBundle {

  /** @ingroup InternalBundleSolver

  */
  //@{


  /** @brief Routine for selecting the weight of the quadratic/proximal term
       withing BundleSolver implementing BundleWeight along the paper by
       Helmberg and Kiwiel.
  */

  class BundleHKWeight : public BundleWeight {
    Groundset* groundset;  ///< the groundset, may not be zero in init and *_update
    BundleModel* model;    ///< may be zero at any time
    CH_Matrix_Classes::Integer iweight; ///< counting the iterations since last change
    CH_Matrix_Classes::Real weight; ///< the current weight
    CH_Matrix_Classes::Real epsweight; ///< value for guessing the function variation
    CH_Matrix_Classes::Real minweight; ///< lower bound on the weight
    CH_Matrix_Classes::Real maxweight; ///< upper bound on the weight
    CH_Matrix_Classes::Real modelmax; ///< largest model value encountered

    bool weightchanged;    ///< true if last choose_* call modified the weight
    bool next_weight_set;  ///< true if set_nextweight was just called, reset at *_update
    int nullstep_updates; ///< true is weight may be changed after null steps [default: true]

    CH_Matrix_Classes::Matrix valuelevels; ///< for judging null step progress
    CH_Matrix_Classes::Matrix ratiolevels; ///< for judging null step progress

    CH_Matrix_Classes::Real mR; ///< parameter for reduction criterion in descent steps

  public:
    /// the parameter mRin gets the value for accepting descent steps, bwp may be used to communicate the previous values use by another routine 
    BundleHKWeight(CH_Matrix_Classes::Real mRin = .5, BundleWeight* bwp = 0, const CBout* cbo = 0, int incr = -1);
    ///
    ~BundleHKWeight() {
    }


    ///set default values for 'constant' parameters, e.g. minweight and maxweight
    virtual void set_defaults();

    ///set nullstep update strategy (0 ... original, 1 ... none, 2 ... enlarge if subsequence of three norm increases is found 
    virtual void set_nullstep_updates(int nu = 0) {
      nullstep_updates = nu;
    }

    ///reset all adaptive variables and parameters
    virtual void clear();

    ///compute first weight and set some parameters
    int init(CH_Matrix_Classes::Real aggr_dnmormsqr, Groundset* groundset, BundleModel* model);

    /// <=0 leaves everything unchanged and does nothing
    void set_next_weight(CH_Matrix_Classes::Real u) {
      if (u <= 0.) return;
      weight = CH_Matrix_Classes::max(u, 1e-10); weightchanged = true; next_weight_set = true;
      modelmax = CB_minus_infinity; iweight = 0;
    }

    /// <=0 means no bound 
    void set_minweight(CH_Matrix_Classes::Real mw) {
      minweight = mw;
      if (minweight > 0) {
        if ((weight > 0) && (weight < minweight))
          weight = minweight;
        if ((maxweight > 0) && (maxweight < minweight))
          maxweight = minweight;
      }
    }

    /// true if the next weight was prespecified externally
    bool get_next_weight_set() const {
      return next_weight_set;
    }

    /// 
    CH_Matrix_Classes::Real get_minweight() const {
      return minweight;
    }

    /// <=0 means no bound
    virtual void set_maxweight(CH_Matrix_Classes::Real mw) {
      maxweight = mw;
      if (maxweight > 0) {
        if (weight > maxweight) {
          weight = maxweight;
        }
        if ((minweight > 0) && (minweight > maxweight))
          minweight = maxweight;
      }
    }

    ///
    CH_Matrix_Classes::Real get_maxweight() const {
      return maxweight;
    }

    /// returns current value of the weight
    CH_Matrix_Classes::Real get_weight() const;

    /// returns true if last call of *_update modified current value of tau, else 0
    bool weight_changed() const;

    /// determine next weight after a descent step
    int descent_update(CH_Matrix_Classes::Real newval,
      CH_Matrix_Classes::Real oldval,
      CH_Matrix_Classes::Real modelval,
      const CH_Matrix_Classes::Matrix& y,
      const CH_Matrix_Classes::Matrix& newy,
      CH_Matrix_Classes::Real normsubg2,
      BundleProxObject* Hp);

    /// determine next weight after a null step
    int nullstep_update(CH_Matrix_Classes::Real newval,
      CH_Matrix_Classes::Real oldval,
      CH_Matrix_Classes::Real modelval,
      const CH_Matrix_Classes::Matrix& y,
      const CH_Matrix_Classes::Matrix& newy,
      MinorantPointer& new_minorant,
      MinorantPointer& aggregate,
      CH_Matrix_Classes::Real nullstep_bound,
      CH_Matrix_Classes::Real normsubg2,
      BundleProxObject* Hp);

    /// reinitialize after modifications
    virtual int apply_modification(const GroundsetModification& gsmdf);

  };

}

#endif

