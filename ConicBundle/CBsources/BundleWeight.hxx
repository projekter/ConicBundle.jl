/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/BundleWeight.hxx
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



#ifndef CONICBUNDLE_BUNDLEWEIGHT_HXX
#define CONICBUNDLE_BUNDLEWEIGHT_HXX


/**  @file BundleWeight.hxx
    @brief Header declaring the class ConicBundle::BundleWeight
    @version 1.0
    @date 2014-08-05
    @author Christoph Helmberg
*/


#include "matrix.hxx"
#include "clock.hxx"
#include "CBSolver.hxx"
#include "MinorantPointer.hxx"
#include "Groundset.hxx"
#include "BundleModel.hxx"

namespace ConicBundle {

  /** @ingroup InternalBundleSolver

  */
  //@{


  /** @brief Abstract interface for BundleSolver providing routines that determine
       the weight of the quadratic term in the augmented model. It also allows the
       user to specify bounds on the weights by setting minweight and maxweight or
       to choose the weight to be used in the next iteration.
  */
  class BundleWeight : public CBout {
  public:
    /// default constructor with the possibility to set the output
    BundleWeight(const CBout* cb = 0, int incr = -1) :CBout(cb, incr) {
    }

    ///
    virtual ~BundleWeight(); //implemented in hk

    /// Sets default values for 'constant' parameters, e.g. minweight and maxweight
    virtual void set_defaults() = 0;

    /// Resets all adaptive variables and parameters
    virtual void clear() = 0;

    /** @brief Compute the first weight u and set some parameters,

       @param[in] aggr_dnormsqr
           the squared dual norm of the current aggregate

       @param[in] groundset
           pointer to the current groundset

       @param[in] model
           pointer to the current model

       @return
       - 0 on success
       - 1 on failure
    */
    virtual int init(CH_Matrix_Classes::Real aggr_dnormsqr,
      Groundset* groundset, BundleModel* model) = 0;

    /** @brief Allows the user to set the weight of the very next iteration of
        the bundle method

        As it is set by the user, the weight will be accepted as long as it is
        strictly positive, even if it is outside the interval
        [minweight,maxweight]. Nonpositive values, however, will leave
        everything unchangend.
    */
    virtual void set_next_weight(CH_Matrix_Classes::Real u) = 0;

    /** @brief Sets a lower bound for the weight. Nonpositive values may be used
        to indicate that the weight is allowed to get arbitrarily close to zero.
     */
    virtual void set_minweight(CH_Matrix_Classes::Real umin) = 0;

    /// Returns the current value of the minweight 
    virtual CH_Matrix_Classes::Real get_minweight() const = 0;

    /** @brief Sets an upper bound for the weight. Nonpositive values may be
        used to indicate that the weight is allowed to get arbitrarily large.
     */
    virtual void set_maxweight(CH_Matrix_Classes::Real umax) = 0;

    /// Returns the current value of the maxweight 
    virtual CH_Matrix_Classes::Real get_maxweight() const = 0;

    /// Returns the current value of the weight 
    virtual CH_Matrix_Classes::Real get_weight() const = 0;

    /// returns true if the current weight has been set externally by the user 
    virtual bool get_next_weight_set() const = 0;

    /** @brief Returns false only if the last call to descent_update() or
        nullstep_update() did not modify the weight and no other routine
        influenced the weight since then, otherwise it returns true
     */
    virtual bool weight_changed() const = 0;

    /** @brief The BundleSolver calls this for computing the next weight if
        the candidate will result in a descent step

      @param[in] newval
          objective value in the candidate (next center)

      @param[in] oldval
          objective value in the current (old) center

      @param[in] modelval
          value of the current model in the candidate

      @param[in] y
          current (old) center

      @param[in] newy
          candidate (next center)

      @param[in] aggr_dnormsqr
          squared dual norm of new_aggregate

      @param[in] Hp
          pointer to scaling information

       @return
       - 0 on success
       - 1 on failure
    */
    virtual int descent_update(CH_Matrix_Classes::Real newval,
      CH_Matrix_Classes::Real oldval,
      CH_Matrix_Classes::Real modelval,
      const CH_Matrix_Classes::Matrix& y,
      const CH_Matrix_Classes::Matrix& newy,
      CH_Matrix_Classes::Real aggr_dnormsqr,
      BundleProxObject* Hp) = 0;

    /** @brief The BundleSolver calls this for computing the next weight if if
        the candidate will result in a null step

      @param[in] newval
          objective value in the candidate

      @param[in] oldval
          objective value in the current center

      @param[in] modelval
          value of the current model in the candidate

      @param[in] y
          current center

      @param[in] newy
          candidate

      @param[in] new_minorant
          the new subgradient obtained in the candidate (only some rules need this)

      @param[in] aggregate
          the aggregate generating the candidate

      @param[in] nullstep_bound
          evaluationts giving a higer value than this result in null steps

      @param[in] aggr_dnormsqr
          squared dual norm of the aggregate

      @param[in] Hp
          pointer to scaling/variable metric information

      @return
       - 0 on success
       - 1 on failure
    */
    virtual int nullstep_update(CH_Matrix_Classes::Real newval,
      CH_Matrix_Classes::Real oldval,
      CH_Matrix_Classes::Real modelval,
      const CH_Matrix_Classes::Matrix& y,
      const CH_Matrix_Classes::Matrix& newy,
      MinorantPointer& new_minorant,
      MinorantPointer& aggregate,
      CH_Matrix_Classes::Real nullstep_bound,
      CH_Matrix_Classes::Real aggr_dnormsqr,
      BundleProxObject* Hp) = 0;


    /** @brief if the BundleSolver is called to modify the groundset it also calls this

        The task is then to update the stored information appropriately so that dimensions match and the proximal term makes sense for the modified problem
    */
    virtual int apply_modification(const GroundsetModification& gsmdf) = 0;

  };

  //@}

}

#endif

