/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/BundleDenseTrustRegionProx.hxx
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


#ifndef CONICBUNDLE_BUNDLEDENSETRUSTREGIONPROX_HXX
#define CONICBUNDLE_BUNDLEDENSETRUSTREGIONPROX_HXX


/**  @file BundleDenseTrustRegionProx.hxx
    @brief Header declaring the class ConicBundle::BundleDenseTrustRegionProx
    @version 1.0
    @date 2017-10-05
    @author Christoph Helmberg
*/


#include "BundleProxObject.hxx"

namespace ConicBundle {

  /**@ingroup InternalBundleProxObject
   */

   //@{

  /** @brief implements the abstract interface ConicBundle::BundleProxObject for \f$\|y-\hat{y}\|_H^2\f$
      for general symmetric H+weight*I (H is assumed to be positive semidefinite without checking)
      giving rise to an augmented model with dense variable metric
  */

  class BundleDenseTrustRegionProx : public BundleProxObject {
  private:
    /// the metric matrix without weightu added to it
    CH_Matrix_Classes::Symmatrix H;

    /// if is_facotred==ture this holds the Cholesky factor of (H+weightu*I)
    mutable CH_Matrix_Classes::Symmatrix Hchol;

    /// true iff Hchol is computed for the current H and weightu
    mutable bool is_factored;

    /// holds the Cholesky factor of (H+weightu*I) with old_fixed_ind deleted if the dimension fits
    CH_Matrix_Classes::Symmatrix Hind_chol;


    /// the current weightu added to H
    CH_Matrix_Classes::Real weightu;

    /// the correction value for correcting termination precision
    CH_Matrix_Classes::Real corr_val;

    /// the fixed indices for which the QP_costs were computed 
    CH_Matrix_Classes::Indexmatrix old_fixed_ind;
    /// for compute_QP_costs() this holds the non_fixed part of the subgradients 
    CH_Matrix_Classes::Matrix _A;
    /// for compute_QP_costs() this holds the constant subgradient
    CH_Matrix_Classes::Matrix _b;
    /// for compute_QP_costs() this holds the offset values of the subgradients
    CH_Matrix_Classes::Matrix _c;
    /// for compute_QP_costs() this holds the constant offset
    CH_Matrix_Classes::Real _delta;
    /// for compute_QP_costs() this holds the non_fixed part of center_y
    CH_Matrix_Classes::Matrix _y;
    /// for compute_QP_costs() this holds L^{-1}*_A 
    CH_Matrix_Classes::Matrix LinvA;
    /// for compute_QP_costs() this holds L^T*_y
    CH_Matrix_Classes::Matrix Lty;
    /// the old quadratic cost matrix 
    CH_Matrix_Classes::Symmatrix oldQ;
    /// the old linear cost term; 
    CH_Matrix_Classes::Matrix oldd;
    /// the old costant cost term; 
    CH_Matrix_Classes::Real oldoffset;


    /// if values should be computed for a new subset indices, this is stored here
    const CH_Matrix_Classes::Indexmatrix* new_indices;

    /// in add_variable_metric() the stack serves to transform the given data
    std::vector<const AffineFunctionTransformation*> aft_stack;

    /// computes the correction value, here min(1,dim/trace(H))
    void compute_corr()
      //{ corr_val=trace(H); if (corr_val>H.rowdim()) corr_val=H.rowdim()/corr_val; else corr_val=1.;}
    {
      corr_val = CH_Matrix_Classes::min(1., H.rowdim() / (H.rowdim() * weightu + trace(H)));
    }

  public:

    /// initialize to this Matrix and set the variable_metric option (false by default)
    BundleDenseTrustRegionProx(const CH_Matrix_Classes::Symmatrix& Hin,
      VariableMetricSelection* vp = 0,
      bool local_metric = false,
      bool bounds_scaling = false,
      CBout* cb = 0,
      int cbinc = -1) :
      BundleProxObject(vp, local_metric, bounds_scaling, cb, cbinc),
      H(Hin), is_factored(false) {
      weightu = 1.;
    }

    /// initialize H to the zero Matrix of this dimension (on the diagonal the weight will be added)  and set the variable_metric option (false by default)
    BundleDenseTrustRegionProx(CH_Matrix_Classes::Integer dim = 0,
      VariableMetricSelection* vp = 0,
      bool local_metric = false,
      bool bounds_scaling = false,
      CBout* cb = 0,
      int cbinc = -1) :
      BundleProxObject(vp, local_metric, bounds_scaling, cb, cbinc),
      H(dim, 0.), is_factored(false) {
      assert(dim >= 0); weightu = 1.;
    }

    ///
    virtual ~BundleDenseTrustRegionProx() {
    }

    /// set the weight of the proximal term
    void set_weightu(CH_Matrix_Classes::Real in_weightu);

    /// returns the current weight of the proximal term
    CH_Matrix_Classes::Real get_weightu() const {
      return weightu;
    }

    ///returns a correction factor for termination precision if the quadratic term is strong
    CH_Matrix_Classes::Real get_term_corr() const {
      return corr_val;
    }

    /// set H with the information, whether it is factored
    const CH_Matrix_Classes::Symmatrix& init(CH_Matrix_Classes::Symmatrix& in_H) {
      is_factored = false; return H = in_H;
    }

    /// returns true iff get_Hchol() returns the factord matrix of H with weightu
    bool get_factored() const {
      return is_factored;
    }

    /// returns the metric matrix without weightu
    const CH_Matrix_Classes::Symmatrix& get_H() const {
      return H;
    }

    /// returns the stored factorization of H with weightu (up to date if get_factored()==true)
    const CH_Matrix_Classes::Symmatrix& get_Hchol() const {
      return Hchol;
    }

    /// returns the order of the matrix
    CH_Matrix_Classes::Integer dim() const {
      return H.rowdim();
    }

    /// returns H(i,j) (without including weightu)
    const CH_Matrix_Classes::Real& operator()(CH_Matrix_Classes::Integer i, CH_Matrix_Classes::Integer j) {
      return H(i, j);
    }

    /// returns $\|B\|^2_H$ (with weightu included in H)
    virtual CH_Matrix_Classes::Real norm_sqr(const CH_Matrix_Classes::Matrix& B) const;

    /// returns $\|B\|^2_{H^{-1}}$ (with weightu included in H)
    virtual CH_Matrix_Classes::Real dnorm_sqr(const MinorantPointer& B) const;

    ///return true if H is of the form diagonal matrix plus Gram matrix of a low rank matrix
    virtual bool is_DLR() const {
      return false;
    }

    /// add H to the dense symmetric matrix as a principal submatrix starting at position start_index
    virtual int add_H(CH_Matrix_Classes::Symmatrix& big_sym, CH_Matrix_Classes::Integer start_index = 0) const;

    /// adds \f$alpha*Hx\f$ to outplusHx and returns this
    virtual CH_Matrix_Classes::Matrix& add_Hx(const CH_Matrix_Classes::Matrix& x,
      CH_Matrix_Classes::Matrix& outplusHx,
      CH_Matrix_Classes::Real alpha = 1.) const {
      outplusHx.xpeya(x, weightu * alpha); return genmult(H, x, outplusHx, alpha, 1.);
    }

    /// returns \f$H^{-1}x\f$
    virtual CH_Matrix_Classes::Matrix& apply_Hinv(CH_Matrix_Classes::Matrix& x) const;

    /// returns a suitable approximation for preconditioning, see BundleProxObject::get_precond 
    virtual void get_precond(CH_Matrix_Classes::Matrix& inD, const CH_Matrix_Classes::Matrix*& Vp) const {
      inD = diag(H); inD += weightu; Vp = 0;
    }


    /// computes the dual QP costs Q, d, and the constant offset to the bundle subproblem, see BundleProxObject::compute_QP_costs
    virtual int compute_QP_costs(CH_Matrix_Classes::Symmatrix& Q,
      CH_Matrix_Classes::Matrix& d,
      CH_Matrix_Classes::Real& offset,
      const MinorantPointer& constant_minorant,
      const MinorantBundle& bundle,
      const CH_Matrix_Classes::Matrix& y,
      const MinorantPointer& groundset_minorant,
      CH_Matrix_Classes::Indexmatrix* yfixed);


    /// updates the dual QP costs Q, d, and the constant offset to the bundle subproblem, see BundleProxObject::update_QP_costs
    virtual int update_QP_costs(CH_Matrix_Classes::Symmatrix& delta_Q,
      CH_Matrix_Classes::Matrix& delta_d,
      CH_Matrix_Classes::Real& delta_offset,
      const MinorantPointer& constant_minorant,
      const MinorantBundle& bundle,
      const CH_Matrix_Classes::Matrix& center_y,
      const MinorantPointer& groundset_minorant,
      const MinorantPointer& delta_groundset_minorant,
      const CH_Matrix_Classes::Indexmatrix& delta_index,
      CH_Matrix_Classes::Indexmatrix* yfixed);


    /// when BundleSolver is called to modify the groundset it also calls this 
    virtual int apply_modification(const GroundsetModification& gsmdf);

    /** @brief in order to allow for fixed variables, this generates a clone restricted to the given indices
   */
    virtual BundleProxObject* projected_clone(const CH_Matrix_Classes::Indexmatrix& indices) {
      CH_Matrix_Classes::Symmatrix S;
      H.principal_submatrix(indices, S);
      for (CH_Matrix_Classes::Integer i = 0; i < S.rowdim(); i++)
        S(i, i) -= weightu;
      BundleDenseTrustRegionProx* pp = new BundleDenseTrustRegionProx(S, 0, get_use_local_metric(), get_use_bounds_scaling(), this, 0);
      pp->set_weightu(weightu);
      pp->apply_factor(factor);
      if (get_variable_metric_selection())
        pp->set_variable_metric_selection(get_variable_metric_selection()->clone_VariableMetricSelection());
      return pp;
    }


    /** @brief this implementation does not support a diagonal scaling heuristic,
        therefore the following routine has to return true.
     */
    virtual bool supports_diagonal_bounds_scaling() const {
      return false;
    }

    /** @brief if supported, D_update has to contain nonnegative numbers that
        are permanently added to the diagonal here. It is important to keep
        track of this change only if afterwards update_QP_costs is called before
        compute_QP_costs. In this case the only nonzero enries in D_update must
        be those of delta_index
     */
    virtual int diagonal_bounds_scaling_update(const CH_Matrix_Classes::Matrix& /* D_update */) {
      return 1;
    }

    /// returns true if dynamic scaling with dense symmetric matrices is supported
    virtual bool supports_dense_variable_metric() const {
      return true;
    }

    /// returns true if dynamic scaling with low rank structure is supported
    virtual bool supports_lowrank_variable_metric() const {
      return true;
    }

    /// returns true if dynamic scaling with diagonal matrices is supported
    virtual bool supports_diagonal_variable_metric() const {
      return true;
    }

    /// see DynamicScaling
    virtual int apply_variable_metric(VariableMetricModel* groundset,
      VariableMetricModel* model,
      const CH_Matrix_Classes::Matrix& aggr,
      CH_Matrix_Classes::Integer y_id,
      const CH_Matrix_Classes::Matrix& y,
      bool descent_step,
      CH_Matrix_Classes::Real& current_weight,
      CH_Matrix_Classes::Real model_maxviol,
      const CH_Matrix_Classes::Indexmatrix* new_indices = 0);

    /// see BundleProxObject::add_variable_metric()
    virtual int add_variable_metric(CH_Matrix_Classes::Symmatrix& addH);

    /// see BundleProxObject::add_lowrank_variable_metric()
    virtual int add_variable_metric(CH_Matrix_Classes::Matrix& diagH,
      CH_Matrix_Classes::Matrix& vecH);

    /// see  BundleProxObject::push_aft();
    virtual int push_aft(const AffineFunctionTransformation* aft) {
      aft_stack.push_back(aft); return 0;
    }

    /// see  BundleProxObject::pop_aft();
    virtual int pop_aft() {
      aft_stack.pop_back();
      return 0;
    }

    /// output the description of the scaling in mfile-suitable format
    virtual int mfile_data(std::ostream& out) const;


  };


  //@}
}

#endif

