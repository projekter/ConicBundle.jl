/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/BoxOracle.hxx
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



#ifndef CONICBUNDLE_BOXORACLE_HXX
#define CONICBUNDLE_BOXORACLE_HXX

/**  @file BoxOracle.hxx
    @brief Header declaring the classes ConicBundle::BoxOracle and ConicBundle::BoxBundleParameters  (needed for ConicBundle::BoxModel)
    @version 1.0
    @date 2017-09-29
    @author Christoph Helmberg
*/

//------------------------------------------------------------

#include "MatrixCBSolver.hxx"
#include "matrix.hxx"

//------------------------------------------------------------


namespace ConicBundle {

  /**@defgroup box_oracle box oracle with special purpose cutting model

     @brief This oracle has the only purpose of passing the box description data
     to a specialized internal implementation of a support function over a
     box. The objective is supplied by a suitable AffineFunctionTransformation .

      Internally, the SumBlockModel implementation of the the support
      function and its bundle model is BoxModel.
  */
  //@{


    /** @brief Interface for extending PrimalData, e.g., in Lagrangian
        relaxation of column generation approaches

        This object has to be created for Modifications in BoxOracle
        that should preserve primal information, if in the course of
        evaluating the oracle one notices, that additional primal
        variables are needed and the old primal variables need to be
        updated accordingly.

        The object will be deleted by ConicBundle after use.

    */

  class BoxPrimalExtender : public PrimalExtender {
  public:
    ///
    virtual ~BoxPrimalExtender() {
    };

    /// like in PrimalExtender, called by ConicBundle to update internal PrimalData objects, has to return 0 on success 
    virtual int extend(PrimalData&) = 0;

    /// called by ConicBundle to update internal Ritz_vectors, has to return 0 on success 
    virtual int extend_Box(CH_Matrix_Classes::Matrix& /* Boxvecs */) = 0;
  };



  /**@brief Oracle interface providing the lower and upper bounds for
     the internally implemented support function over this box for an
     affine function given by an AffineFucntionTransfomation or,
     equivalently, Lagrangian relaxation of linear programs over box domains.
     No modifications of the box groundset are supported so far.

     The support function is fully implemented in BoxModel.
  */

  class BoxOracle : public ModifiableOracleObject {
  private:
    CH_Matrix_Classes::Matrix lb;  ///< column vector of lower bounds, same dimension as ub
    CH_Matrix_Classes::Matrix ub;  ///< column vector of upper bounds, same dimension as lb

  public:

    /// constructor initializing lower and upper bounds (must have the same dimesnion, not checked)
    BoxOracle(const CH_Matrix_Classes::Matrix& in_lb,
      const CH_Matrix_Classes::Matrix& in_ub) :
      lb(in_lb), ub(in_ub) {
    }

    /// destructor
    virtual ~BoxOracle();

    /** @brief returns the lower bounds vector of the box
     */
    const CH_Matrix_Classes::Matrix& get_lower_bounds() {
      return lb;
    }

    /** @brief returns the upper bounds vector of the box
     */
    const CH_Matrix_Classes::Matrix& get_upper_bounds() {
      return ub;
    }

    /**@brief This routine need not be implemented unless variables
      (constraints in Lagrangean relaxation) are added or deleted on
      the fly

      The routine is only called by the solver if the variables indeed get
      modified by the solver or a modification is passed on by the user via
      the solver interface. @a oracle_modification is then used to either
      transfer user supplied instructions to the oracle on how to modify
      itself or to inform the oracle about changes the solver was asked to
      perform on the variables. If available, the solver will also show the
      effect of these changes on the center point in @a new_center and @a
      old_center; if these are not available then they hold NULL. A user
      supplied @a oracle_modification will be checked for consistency with the
      actual changes in the variables and mismatches will cause failures.

      The remaining variables are output variables by which the oracle tells
      the solver which information has a chance to be preserved in view of
      these changes.  If e.g. the deletion of some nonzero variables
      invalidates the function value in the new center, the oracle has to set
      discard_objective_in_center=true.  If the entire model cannot be
      preserved (this includes the aggregates and the function values), the
      oracle needs to set discard_model=true; If only aggregate minorants
      cannot be preserved, the oracle needs to set discard_aggregates=true; in
      the current implementation of BoxModel this removes all minorants
      generated by generate_minorant() because each of them is typically
      generated by more than one Ritz vector. Whenever new variables were
      added, the model can only be preserved if the remaining minorants (maybe
      without aggregates) can be extended for these new variables. In this
      case the oracle has to supply the appropriate MinorantExtender via @a
      minorant_extender and only those minorants will be kept for which this
      operation succeeds.

      Return value 0 indicates that these actions allow to continue without
      errors, other return values result in an overall error on these changes.
    */
    virtual
      int
      apply_modification
      (
        const OracleModification& /* oracle_modification */,
        const CH_Matrix_Classes::Matrix* /* new_center */,
        const CH_Matrix_Classes::Matrix* /* old_center */,
        bool& /* discard_objective_in_center */,
        bool& /* discard_model */,
        bool& /* discard_aggregates */,
        MinorantExtender*& /* minorant_extender */
      ) {
      return 1;
    }



    /**@brief switch on/off some correctnes checks on the oracle */
    virtual
      bool
      check_correctness() const {
      return true;
    }


  };


  //@}

}
#endif

