/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/BundleProxObject.hxx
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



#ifndef CONICBUNDLE_BUNDLEPROXOBJECT_HXX
#define CONICBUNDLE_BUNDLEPROXOBJECT_HXX


/**  @file BundleProxObject.hxx
    @brief Header declaring the class ConicBundle::BundleProxObject
    @version 1.0
    @date 2014-07-23
    @author Christoph Helmberg
*/


#include "AffineFunctionTransformation.hxx"
#include "QPSolverObject.hxx"
#include "VariableMetric.hxx"

namespace ConicBundle {

class BundleModel;
class Groundset;

/** @defgroup InternalBundleProxObject  Quadratic Proximal Terms

   @brief Together with a weight (see BundleWeight), the proximal term
   \f$\|y-\hat{y}\|_H^2\f$ (\f$H\f$ positive definite) of the
   augmented cutting model plays the role of a step length control and
   offers a possibility to bring in some metric information via the
   choice of \f$H\f$ (see \ref InternalVariableMetric). Two aspects
   are of importance in choosing this matrix.
  
   - First order methods are, unfortunately, extremely sensitive to
     scaling. Depending on the nature of the problem it may be
     possible to choose \f$H\f$ so as to add some information about
     the curvature (e.g. the positive semidefinite part of the
     Hessian) or estimates and reasonable relative step sizes into the
     model. This step size also plays a role in the selection of
     subgradient information wihtin the cutting model.

   - Depending on the structural properties of \f$H\f$ (diagonal, low
     rank, full, etc.) it may be computationally very costly or rather
     easy to solve the various QP suproblems arising in the bundle
     method

   The classes of this section help to set up appropriate
   specializations that allow to exploit structural properties of
   \f$H\f$ and to collect/compute the cost-coefficients for the QPs
   independent of the actual cutting models and ground sets in use.
   They may be brought to use via MatrixCBSolver::set_prox().

   The abstract base class for providing quadratic prox terms is BundleProxObject .
   ConicBundle offers several implementations of this, some of them
   support VariableMetric, some do not. Currently all use the BundleWeight
   as an additive term on the diagonal in a trust region fashion:

   - BundleIdProx 
     (Identity, no dynamic scaling, weight as a factor)

   - BundleDiagonalTrustRegionProx 
     (diagonal matrix, offers dynamic diagonal/bounds scaling, additive trust region weight)

   - BundleLowRankTrustRegionProx 
     (low rank Gram matrix, offers low rank variable metric, additive trust region weight)

   - BundleDLRTrustRegionProx (not properly tested yet)
     (low rank + diagonal, offers low rank/diagonal variable metric, additive trust region weight)

   - BundleDenseTrustRegionProx (dense symmetric, offers dense/low rank/diagonal variable metric, additive trust region weight; mainly for comparative testing)

   
*/

//@{



/** @brief abstract interface that allows to use different \f$H\f$-norms
     \f$\|y-\hat{y}\|_H^2\f$ with a positive definite matrix \f$H\f$ in the
     proximal term of the augmented model of ConicBundle::BundleSolver.

     ConicBundle::BundleIdProx is the identitiy \f$H=I\f$. There are
     further variants for diagonal scaling, low rank scaling or full scaling
     by a positive definite matrix.
 */


  class BundleProxObject: public VariableMetric, public QPSolverProxObject 
{
protected:
  /// constant for possible use if QP coefficients are to be computed in parallel (still purely experimental)
  static const CH_Matrix_Classes::Integer xdim_threshold = 100; 

  /// used to accumulate a compensation factor for function_factor; this factor is not included in H but can be applied externally via get_factor() if required
  CH_Matrix_Classes::Real factor;

  ///the QP may signal short steps that seem due to the quadratic term by setting this counter via set_short_QPsteps() 
  CH_Matrix_Classes::Integer short_QPsteps;

public:
  /// default constructor, switching on dynamic scaling only works for classes with corresponding support
  BundleProxObject(VariableMetricSelection* vp=0,bool local_scaling=false,bool bounds_scaling=false,CBout* cb=0,int cbincr=-1);
  //:DynamicScaling(dynamic_scaling,local_scaling,bounds_scaling){}
  //somehow inline did not work here, so I moved it to BundleProxObject.cxx

  ///
  virtual ~BundleProxObject(); //implemented in BundleScaling.cxx

  ///in future computations the following weight must be included in the quadratic term (replacing the previous one)
  virtual void set_weightu(CH_Matrix_Classes::Real weightu)=0;

  /// may be used to indicate seemingly conservative step sizes possibly due to the quadratic term
  virtual void set_short_QPsteps(CH_Matrix_Classes::Integer shortQPst)
  {short_QPsteps=shortQPst;}
  
  /// retrieves the number of conservative step sizes possibly due to the quadratic term passed on to this
  virtual CH_Matrix_Classes::Integer get_short_QPsteps()
  {return short_QPsteps;}
  
  
  /// allows AFTModel and SumModel to accumulate a compensation factor for tracing the effects of recursive applications of function_factor in update_model (does not affect H but can be retrieved by get_factor() for this purpose)
  int apply_factor(CH_Matrix_Classes::Real f)
  {if (f<=0.) return 1; factor*=f; return 0; }

  /// returns the current accumulated compensation factor by which H would need to be scaled in order to reflect the curvature relative to the current function without function_factor;   
  CH_Matrix_Classes::Real get_factor()
  {return factor;}
  
  ///return the current weight incorporated in the quadratic term
  virtual CH_Matrix_Classes::Real get_weightu() const=0;

  ///return a correction factor for termination precision if the quadratic term is strong
  virtual CH_Matrix_Classes::Real get_term_corr(void) const {return 1.;}
   
  ///return \f$\|B\|^2_H\f$
  virtual CH_Matrix_Classes::Real norm_sqr(const CH_Matrix_Classes::Matrix& B) const =0;
  ///return \f$\|B\|^2_{H^{-1}}\f$
  virtual CH_Matrix_Classes::Real dnorm_sqr(const MinorantPointer& B) const =0;

  ///return true if H is of the form diagonal matrix plus Gram matrix of a low rank matrix
  virtual bool is_DLR() const =0;

  /// add H to the dense symmetric matrix as a principal submatrix starting at position start_index
  virtual int add_H(CH_Matrix_Classes::Symmatrix& big_sym,CH_Matrix_Classes::Integer start_index=0) const=0;

  ///add \f$alpha*Hx\f$ to outplusHx and return this
  virtual CH_Matrix_Classes::Matrix& add_Hx(const CH_Matrix_Classes::Matrix& x,
					    CH_Matrix_Classes::Matrix& outplusHx,
					    CH_Matrix_Classes::Real alpha=1.) const =0;

  /** @brief returns \f$H^{-1}x\f$, possibly transformed by some AffineFunctionTransformation

      If a stack of instances of AffineFunctionTransformation has been added
      by push_aft() and A_1 to A_k are on the stack with A_k on top,
      then it returns \f$ A_k\cdots A_1\cdot H^{-1}\cdot A_1^\top\cdots A_k^\top\cdot x\f$
   */
  virtual CH_Matrix_Classes::Matrix& apply_Hinv(CH_Matrix_Classes::Matrix& x) const =0;

  /** @brief return \f$Diag(D)+VV^T\f$ which either equals $H$ exactly (iff is_DLR() returns true)
      or, if not, serves as an approximation hopefully suitable for
      preconditioning. If $V$ is available it may be large, so return only a
      pointer to avoid copying. If no $V$ is available, then put Vp==0 on
      exit.
  */
  virtual void get_precond(CH_Matrix_Classes::Matrix& D, 
			   const CH_Matrix_Classes::Matrix*& Vp) const=0;
  
  /** @brief Compute the dual QP costs Q, d, and the constant offset
      to the bundle subproblem
      
      Let \f$Y\f$ denote the (convex) feasible ConicBundle::Groundset
      and \f$i_Y\f$ its indicator function, let \f$\hat y\f$ be the
      center of stability given by @a center_y, let \f$ (\gamma,g) \f$
      be the aggregate ground set minorant \f$\gamma+g^\top y\le
      i_Y(y)\f$ described by @a gs_subg_offset and @a gs_subg, let
      \f$H\f$ denote the positive definite scaling matrix with weight
      \f$u\f$ given by this class, and let \f$f\f$ denote the
      objective for which ConicBundle::BundleModel holds the model
      \f$W\subseteq\{(\sigma,s)\colon \sigma+s^\top y\le f(y)\ \forall
      y\in Y\}\f$, which is assumed to be described in the form
      \f$W=\{(\delta+c^\top x,b+Ax):x\in X\}\f$ for some compact
      convex set \f$X\f$ of dimension @a xdim.  The routine computes
      the cost coefficients for the convex quadratic problem over
      \f$x\in X\f$ corresponding to the saddle problem

       \f[ \max_{x\in X}\min_{y\in\mathbf{R}^n} [(b+Ax+g)^\top y + \gamma+\delta+c^\top x +\frac{u}2\|y-\hat y\|_H^2],\f]

This is indeed a convex QP in x because for each x the inner
      minimization over y is unconstrained and can be solved explicitly. So
      this gives rise to a problem

       \f[ \max_{x\in X} -\frac12 x^\top Qx+d^\top x+\rho\f]

How to compute the cost coefficients \f$Q, d, \rho\f$ efficiently depends
     on the structural properties of \f$H\f$ and therefore the costs are
     computed in this class.

     The data \f$[\delta,b]\f$ is provided by the offset constant minorant
     (its offset is \f$\delta\f$, its linear part \f$b\f$), likewise
     the columns of \f$[c^T;A]\f$ are provided by the xdim minorants of the
     bundle and \f$x\in X\f$ describes the admissible convex combinations
     which are not needed here. \f$H\f$ and \f$u\f$
     are the quadratic term and the weight, the latter is set in
     set_weightu(). The cost coefficients are \f$Q=\frac1u A^\top H^{-1}A\f$,
     d=...

If some of the design variables are at their bounds and likely to
     only move out of bounds next they might be required to stay fixed
     by the calling routine and thus may not contribute to the cost
     coefficients of the dual QP other than by their contribution to
     the constant offset. These variables are indicated by a positive
     (nonzero) value in yfixed.
    
     @param[out] Q (CH_Matrix_Classes::Symmatrix&)
         the quadratic cost term to be computed by this routine

     @param[out] d (CH_Matrix_Classes::Matrix&)
         the linear cost term to be computed by this routine

     @param[out] offset (CH_Matrix_Classes::Real&)
         the constant cost offset term to be computed by this routine

     @param[in] constant_minorant (const MinorantPointer&)
         together with the groundset_minorant its offset 
	 will yield \f$\delta\f$, and its linear part \f$b\f$),

     @param[in] bundle (const MinorantBundle&)
         their offsets form \f$c\f$, their linear parts the columns of \f$A\f$.

     @param[in] center_y (const CH_Matrix_Classes::Matrix&)
         the current center of stability of the bundle method,
         it is needed for computing linear cost term and 
         constant coefficient

     @param[in] groundset_minorant (const MinorantPointer&)
         this is the aggregate subgradient \f$\eta\f$ of the ground set
         it acts like a shift of the vector b forming the linear cost
         term for the y variables. 

     @param[in,out] yfixed (CH_Matrix_Classes::Indexmatrix*)
         - if not NULL, on input the values of yfixed[i] are
           + 0  if y_i us to be treated as a normal/non-fixed variable 
           + 1  if y_i is to be treated as a fixed value again
           + 2  if y_i is newly fixed for the first time
         - on output the last group, those with value 2 have to
           be set to 1 so that only values 0 and 1 remain  

     @return
       - 0 on success
       - 1 on failure
  */
  virtual int compute_QP_costs(CH_Matrix_Classes::Symmatrix& Q,
			       CH_Matrix_Classes::Matrix& d,
			       CH_Matrix_Classes::Real& offset,
                               const MinorantPointer& constant_minorant,
                               const MinorantBundle& bundle,
			       const CH_Matrix_Classes::Matrix& center_y,
			       const MinorantPointer& groundset_minorant,
			       CH_Matrix_Classes::Indexmatrix* yfixed) =0;

  

  /** @brief computes the update of the dual QP cost terms d and offset returned by
      compute_QP_costs() for changes in the ground set aggregate \f$(\gamma,\eta)\f$ and in @a yfixed.
  
     The update of the groundset aggregate \f$g\f$ is given by @a delta_gs_subg and @a delta_index, more precisely, old_gs_subg[delta_index[i]]=gs_subg[delta_index[i]]-delta_gs_subg[i], where the current groundset aggregate @a gs_subg is the input paramter. @a delta_gs_subg_offset gives the change of \f$\gamma\f$.

     The entries of @a yfixed have the same meaning as in  compute_QP_costs(): if yfixed[i]==0, variable y[i] is treated normally (not fixed); for yfixed[i]==1 the value was already considered fixed at center_y[i] in the computations before so it causes no changes at all; for yfixed[i]==2 the value is considered newly fixed at center_y[i], so its influence in the quadratic cost term Q and the linear term d have to be corrected (this case is the only reason for having @a delta_Q in the list of arguments). After this correction, put yfixed[i]=1 to signal that in the QP costs coordinate i is considered as fixed.


    @param[out] delta_Q (CH_Matrix_Classes::Symmatrix&)
         the change in the quadratic cost term to be computed 
         by this routine

     @param[out] delta_d (CH_Matrix_Classes::Matrix&)
         the change in linear cost term to be computed by this routine

     @param[out] delta_offset (CH_Matrix_Classes::Real&)
         the change in the constant cost offset term to be computed 
          by this routine

     @param[in] constant_minorant (const MinorantPointer&)
         must be identical to the one used in compute_QP_costs

     @param[in] bundle (const CH_Matrix_Classes::Integer)
         must be identical to the one used in compute_QP_costs

     @param[in] center_y (const CH_Matrix_Classes::Matrix&)
         the current center of stability of the bundle method,

     @param[in] groundset_minorant (const MinorantPointer&)
         this is the aggregate subgradient \f$g\f$ of the ground set

     @param[in] delta_groundset_minorant (const MinorantPointer&)
         this is the delta added to the previous value of the 
         groundset aggregate to get the current aggregate 
	 subgradient \f$g\f$=groundset_minorant of the ground set

     @param[in] delta_index (const CH_Matrix_Classes::Indexmatrix&)
         holds the inidices where groundset_minorant or yfixed changed

     @param[in,out] yfixed (CH_Matrix_Classes::Indexmatrix*)
         - if not NULL, on input the values of yfixed[i] are
           + 0  if y_i us to be treated as a normal/non-fixed variable 
           + 1  if y_i is to be treated as a fixed value again
           + 2  if y_i is newly fixed for the first time
         - on output the last group, those with value 2 have to
           be set to 1 so that only values 0 and 1 remain  

    @return
       - 0 on success
       - 1 on failure

    @see compute_QP_costs
  */
  virtual int update_QP_costs(CH_Matrix_Classes::Symmatrix& delta_Q, 
			      CH_Matrix_Classes::Matrix& delta_d,  
			      CH_Matrix_Classes::Real& delta_offset,
                              const MinorantPointer& constant_minorant,
                              const MinorantBundle& bundle,
			      const CH_Matrix_Classes::Matrix& center_y,
			      const MinorantPointer& groundset_minorant,
			      const MinorantPointer& delta_groundset_minorant,
			      const CH_Matrix_Classes::Indexmatrix& delta_index,
			      CH_Matrix_Classes::Indexmatrix* yfixed) =0;

  /** @brief if the BundleSolver is called to modify the groundset it also calls this 

      The task is then to update the diagonal scaling information appropriately so that dimensions match and the proximal term makes sense for the modified problem.
  */
  virtual int apply_modification(const GroundsetModification& gsmdf)=0;

  /** @brief in order to allow for fixed variables, this generates a clone restricted to the given indices
   */
  virtual BundleProxObject* projected_clone(const CH_Matrix_Classes::Indexmatrix& indices) =0;

};


  //@}

}

#endif

