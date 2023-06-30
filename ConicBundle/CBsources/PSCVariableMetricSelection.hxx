/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/PSCVariableMetricSelection.hxx
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



#ifndef CONICBUNDLE_PSCVARIABLEMETRICSELECTION_HXX
#define CONICBUNDLE_PSCVARIABLEMETRICSELECTION_HXX


/**  @file PSCVariableMetricSelection.hxx
    @brief Header declaring the classes ConicBundle::PSCVariableMetricSelection
    @version 1.0
    @date 2019-03-12
    @author Christoph Helmberg
*/


#include "PSCOracle.hxx"
#include "VariableMetric.hxx"

namespace ConicBundle {


  /** @ingroup InternalVariableMetric   */


//@{


  /** @brief iplementation of a VariableMetricSelection routine for PSCModel for forming and adding variable metric information to a BundleProxObject

   The implementation builds on the method describe by Helmberg, Overton and Rendl.

   The amount of data available for forming the metric depends
   curcially on the model selection routine employed in
   PSCModelParameters::select_model().  Maybe the choice in particular
   of skippedsize there is more relevant than the parameters here.

   Currently there is only one method and the setting of @a selection_method is irrelevant.
 */

  class PSCVariableMetricSelection : public VariableMetricSelection {
  private:
    /// specifies which of the current experimental routines should be used for computing the metric
    CH_Matrix_Classes::Integer selection_method;

    /// for values of oldfactor between 0. and 1. the SVD is taken of (1-oldfactor)*vecH*vecH'+oldfactor*oldvecH*oldvcH'; for values <=0 there is no old contribution, for values >=1 oldvecH is reused (not yet implemented)
    CH_Matrix_Classes::Real oldfactor;

    /// maxeigval_factor*weightu is the maximum value allowed for low_rank eigenvalues, bigger eigen_values are cut off to avoid numerical diffculties.   
    CH_Matrix_Classes::Real maxeigval_factor;

    /// mineigval_factor* max_eigval is the smallest low rank eigenvalue to be included in the low rank approximation (relative to the biggest eigenvaleu max_eigval of the low rank approximation)
    CH_Matrix_Classes::Real mineigval_factor;

    PSCOracle* oracle; ///< needed for PSCOracle::left_right_product(), set once in the beginning but not owned by *this

    CH_Matrix_Classes::Matrix lamH; ///< local temporary storage for low rank eigenvalues
    CH_Matrix_Classes::Matrix vecH; ///< local temporary storage for low rank eigenvectors
    CH_Matrix_Classes::Matrix tmpQ1; ///< local temporary storage for active eigenspace basis 
    CH_Matrix_Classes::Matrix tmpQ2; ///< local temporary storage for inactive eigenspace basis 
    CH_Matrix_Classes::Matrix Q1AQ2; ///< local temporary storage for Q_1^T*A*Q2
    CH_Matrix_Classes::Matrix tmpvec; ///< local temporary storage
    CH_Matrix_Classes::Symmatrix S; ///< local storage for eigenvalue/SVD computations

  public:
    /// destructor
    ~PSCVariableMetricSelection() {
    }

    /// default constructor
    PSCVariableMetricSelection(CBout* cb = 0, int cbincr = -1) :
      VariableMetricSelection(cb, cbincr),
      selection_method(-1),
      oldfactor(0.),
      maxeigval_factor(1e6),
      mineigval_factor(1e-3),
      oracle(0) {
    }

    ///returns selection_method
    virtual CH_Matrix_Classes::Integer get_selection_method() const {
      return selection_method;
    }

    ///sets selection_method 
    virtual void set_selection_method(CH_Matrix_Classes::Integer sm) {
      selection_method = sm;
    }

    ///returns the parameter
    virtual CH_Matrix_Classes::Real get_oldfactor() const {
      return oldfactor;
    }

    ///sets the parameter
    virtual void set_oldfactor(CH_Matrix_Classes::Real of) {
      oldfactor = CH_Matrix_Classes::min(1., CH_Matrix_Classes::max(0., of));
    }

    ///returns the parameter
    virtual CH_Matrix_Classes::Real get_maxeigval_factor() const {
      return maxeigval_factor;
    }

    ///sets the parameter
    virtual void set_maxeigval_factor(CH_Matrix_Classes::Real ef) {
      maxeigval_factor = CH_Matrix_Classes::max(1., ef);
    }

    ///returns the parameter
    virtual CH_Matrix_Classes::Real get_mineigval_factor() const {
      return mineigval_factor;
    }

    ///sets the parameter
    virtual void set_mineigval_factor(CH_Matrix_Classes::Real ef) {
      mineigval_factor = CH_Matrix_Classes::max(CH_Matrix_Classes::eps_Real, CH_Matrix_Classes::min(1., ef));
    }

    ///sets the oracle pointer to this value (NULL is allowed, but calling add_variable_metric() then results in a WARNING and an error is returned); this is called by PSCModel when installing this
    virtual void set_oracle(PSCOracle* psco) {
      oracle = psco;
    }

    ///see ConicBundle::VariableMetricSelection::add_variable_metric(); here it must be possible to cast bundle_data to PSCData&, otherwise the routine returns an error. It only adds something after a descent_step
    int add_variable_metric(VariableMetric& H,
      CH_Matrix_Classes::Integer y_id,
      const CH_Matrix_Classes::Matrix& y,
      bool descent_step,
      CH_Matrix_Classes::Real weightu,
      CH_Matrix_Classes::Real model_maxviol,
      const CH_Matrix_Classes::Indexmatrix* indices,
      VariableMetricBundleData& bundle_data);

    /// clone: the values are only preserved for those contained in the constructor: n_latest_minorants, selection_method and oldfactor
    VariableMetricSelection* clone_VariableMetricSelection() {
      return new PSCVariableMetricSelection(this, 0);
    }
  };



  //@}

}

#endif

