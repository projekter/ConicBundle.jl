/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/UnconstrainedGroundset.hxx
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



#ifndef CONICBUNDLE_UNCONSTRAINEDGROUNDSET_HXX
#define CONICBUNDLE_UNCONSTRAINEDGROUNDSET_HXX


/**  @file UnconstrainedGroundset.hxx
    @brief Header declaring the class ConicBundle::UnconstrainedGroundset
    @version 1.0
    @date 2019-08-21
    @author Christoph Helmberg
*/

#include "Groundset.hxx"
#include "LPGroundsetModification.hxx"

namespace ConicBundle {

  /** @ingroup InternalGroundset */
  //@{

  /** @brief implements an unconstrained groundset

      This class implements an unconstrained ground set.


   */

  class UnconstrainedGroundset :virtual public Groundset {
  protected:

    /// current dimension of the ground set
    CH_Matrix_Classes::Integer dim;

    /// nonnegative update counter for recognizing changes in the groundset description
    CH_Matrix_Classes::Integer groundset_id;

    /// the starting point
    CH_Matrix_Classes::Matrix starting_point;

    /// the ground set aggregate, here giving the linear coefficients of the affine cost function
    MinorantPointer gs_aggregate;

    /// if set to true in constrained versions, y-variables that are at strongly active bounds should be fixed to this value by indicating so in @a yfixed
    bool use_yfixing;

    /// in constrained versions, when @a use_yfixing is true, this indicates whether a coordinate should be considered as free (0), fixed (1), or newly fixed (2)
    CH_Matrix_Classes::Indexmatrix yfixed;

    ///the QPSolver for solving the bundle subproblem with Groundset
    QPSolverObject* qp_solver;

    //data for describing the QP
    const BundleProxObject* Hp;  ///< points to the quadratic cost matrix inside ensure_feasibility() and candidate(), otherwise ==NULL
    CH_Matrix_Classes::Matrix c; ///< linear cost term for QP-subproblems
    CH_Matrix_Classes::Real gamma; ///< constant offset for QP-subproblems


  public:
    /** @brief reset everything to initial state for an unconstrained ground set of dimension @a indim

        Note that a ground set is allowed to have dimension zero. This
        will lead to evaluating a function without arguments and is
        a realistic scenario in Lagrangean relaxation of cutting plane
        approaches if no cutting planes have been added yet.

        If the changes to the ground set are to be counted by groundset_id,
        then it makes sense to enter the appropriate value in in_groundset_id.
     */
    virtual void clear(CH_Matrix_Classes::Integer indim = 0, CH_Matrix_Classes::Integer in_groundset_id = 0);

    /// calls clear() with the same parameters
    UnconstrainedGroundset(CH_Matrix_Classes::Integer indim = 0,
      const CH_Matrix_Classes::Matrix* start_val = 0,
      const CH_Matrix_Classes::Matrix* costs = 0,
      const CH_Matrix_Classes::Real offset = 0.,
      CH_Matrix_Classes::Integer in_groundset_id = 0);

    virtual ~UnconstrainedGroundset();

    /// returns the current groundset_id, increased values indicate changes in the ground set  
    CH_Matrix_Classes::Integer get_groundset_id() const {
      return groundset_id;
    }

    /// sets the groundset_id to the desired value, increasing it is safer here because this is used to indicate changes
    void set_groundset_id(CH_Matrix_Classes::Integer gsid) {
      groundset_id = gsid;
    }

    /// returns the dimension of the ground set, i.e., the length of the variables vector y
    virtual CH_Matrix_Classes::Integer get_dim() const {
      return dim;
    }

    /** @brief returns false if the feasible set is the entire space (unconstrained optimization), true otherwise.

        The current class implements the unconstrained case and always returns false.
    */
    virtual bool constrained() const {
      return false;
    }

    /** @brief on input value in_groundset_id the input y was feasible. Return true if the id did not change, otherwise check if y is still feasible for the given precision.

     The routine is called by the internal bundle solver to check
     whether the given center is still valid (in some applications the
     groundset might change during the runtime of the bundle method),
     where validity of y was already checked at a point in time when the
     groundset had the in_groundset_id. If the groundset_id is still the
     same, then y is simply assumed to be still correct (the precision
     is not even looked at in this case). Otherwise the routine checks
     the validitiy of y with respect to the given precision but does not
     enforce validity. It returns true if y is valid and false otherwise.
    */
    virtual bool is_feasible(CH_Matrix_Classes::Integer& in_groundset_id,
      const CH_Matrix_Classes::Matrix& y,
      CH_Matrix_Classes::Real relprec = 1e-10);

    /** @brief if the groundset_id changed, it checks feasibility of y with respect to the given precision. If infeasible it replaces y by its projection with respect to the norm of Hp and sets ychanged to true.

     The routine is called by the internal bundle solver to check
     whether the given center is still valid (in some applications the
     groundset might change during the runtime of the bundle method),
     where, if @a ychanged is false on input, validity of @a y was
     already checked at a point in time when the groundset had
     the @a in_groundset_id. If @a ychanged==false and the groundset_id is
     still the same, then @a y is simply assumed to be still correct
     (the precision is not even looked at in this case).  Otherwise
     the routine checks the validitiy of @a y with respect to the given
     precision. If feasible, it returns the new
     groundset_id in @a in_groundset_id and keeps @a ychanged unaltered.
     If @a y is infeasible, the rountine computes its projection
     onto the feasible set with respect to the norm of @a Hp (if ==0 then
     the Euclidean norm is used), stores it in @a y, sets @a ychanged to
     true, sets @a in_groundset_id to the current groundset_id and returns
     0. Should anything go wrong, it returns 1.

     This concrete base class represents the unconstrained case, so
     feasiblity only checks the dimension and never requires projections.
     */

    virtual int ensure_feasibility(CH_Matrix_Classes::Integer& in_groundset_id,
      CH_Matrix_Classes::Matrix& y,
      bool& ychanged,
      BundleProxObject* Hp = 0,
      CH_Matrix_Classes::Real relprec = 1e-10);

    /// returns a pointer to an internal QPSolverObject that is able to solve bundle suproblems efficiently for this kind of groundset and scaling; if solves_model_without_gs == true the qp solver does not include the groundset and the groundset has to be dealt with by the Gauss Seidel approach 
    virtual QPSolverObject* get_qp_solver(bool& solves_model_without_gs,
      BundleProxObject* Hp);

    /// set parameters for the QP_Solver  
    int set_qp_solver_parameters(QPSolverParametersObject* /* in_qpsp */) {
      return 0;
    }

    /// returns a stored starting point, note: this need not be feasible; if generated automatically, its dimension is correct.  
    virtual const CH_Matrix_Classes::Matrix& get_starting_point() const {
      return starting_point;
    }

    /// stores the a new starting point irrespective of whether it is feasible or not and returns 0 if it feasible, 1 if it is infeasible
    virtual int set_starting_point(const CH_Matrix_Classes::Matrix& vec) {
      starting_point = vec; int dummy = -1; return (is_feasible(dummy, vec) == false);
    }

    /** @brief for a given model aggregate compute the groundset aggregate and the resulting (feasible) candidate

  Let \f$ (\sigma,s) \f$ be the aggregate minorant \f$\sigma+s^\top y\le
  f(y)\f$ of the cost function described by @a model_subg_offset and @a
  model_subg, let \f$\hat y\f$ be the center of stability given by @a
  center_y, let \f$H\f$ denote the positive definite scaling matrix with
  weight \f$u\f$ given by @a Hp, let \f$Y\f$ denote this (convex)
  feasible ground set and \f$i_Y\f$ its indicator function, then this
  computes a saddle point \f$(y,(\gamma,g))\f$ of

  \f[ \max_{(\gamma,g)\in\partial i_Y}\min_{y\in\mathbf{R}^n} [(s+g)^\top y + \gamma+\sigma +\frac{u}2\|y-\hat y\|_H^2],\f]

  The resulting \f$y\f$ is feasible and stored in @a newy.
       @a linval gets \f$\gamma+\sigma+(g+s)^\top y\f$,
       @a subgnorm2 gets \f$u\|y-\hat x\|_H^2\f$,
       @a augval gets @a linval+@a subgnorm2 /2.

  If @a delta_groundset_subg is not NULL, also
       @a elta_groundset_subg_offset and @a delta index are assumed to be
       not NULL. Then all changes from the previous groundset aggregate
       minorant to the new groundset aggregate minorant \f$(\gamma,g)\f$ are
       stored in @a delta_groundset_subg_offset and @a delta_groundset_subg
       and @a delta_index holds the indices of the nonzero changes
       (mostly the groundset aggregate is sparse, e.g. due to complementarity).
       This is used in BundleProxObject::update_QP_costs().

  On input \f$\varepsilon=\f$ @a relprec, \f$\bar f=\f$ @a
       center_value, and \f$\underline{f}=\f$ @a augval serve to form
       appropriate stopping criteria if solving the saddle point problem
       requires a nonlinear convex optimization method. The method is
       then assumed to produce a primal solution \f$y\in Y\f$ of value
       \f$\bar s\f$ and a dual solution \f$(\gamma,g)\f$ of value
       \f$\underline{g}\f$ with the properties \f$\underline g-\underline
       f\ge 0\f$ and \f$\bar f-\bar s\ge 0\f$ and \f$\bar
       s-\underline{g}\le\varepsilon(|\bar s|+1.)\f$.  In particular
       \f$y\in Y\f$ is assumed to hold to machine precision.

  If use_yfixing is true (the fixing heuristic is switched on),
       then \f$y_i=\hat y_i\f$ is required to hold for all i with
       yfixed(i)!=0, so these coordinates are not allowed to change. In
       particular, this routine may also set yfixed(i)=2 for new
       coordinates i, where 2 is used to indicate newly fixed
       variables. These will be reset to 1 in
       BundlesScaling::compute_QP_costs() and
       BundlesScaling::update_QP_costs()
       when this information has been digested.

  In some derived classes a scaling heuristic is called that may influence
  the scaling @a Hp so as to avoid going outside the feasible region too far.

    @param[out] gs_id
        the current groundset_id

    @param[out] newy
        the next candidate y (feasible)

    @param[out] cand_gs_val
        the value of the groundset minorant in the candidate y (=groundset objective)

    @param[out] linval (CH_Matrix_Classes::Real&) value of linear minorant in y

    @param[in,out] augval_lb (CH_Matrix_Classes::Real&)
        - on input: lower bound value of augmented model in previous (maybe infeasible) candidate
        - on output: lower bound value of augmented model in y

    @param[out] augval_ub (CH_Matrix_Classes::Real&)
        - on output: upper bound value of augmented model in y

    @param[out] subgnorm2 (CH_Matrix_Classes::Real&)
        squared Hp-norm (with weight) of joint groundset and model aggregate

    @param[in] center_y (const CH_Matrix_Classes::Matrix&)
        center of stability (feasible)

    @param[in] center_value (CH_Matrix_Classes::Real)
        function value in center_y

    @param[in] model_minorant (const MinorantPoiner&)
        aggregate linear minorant of the cost function

    @param[in,out] Hp (ConicBundle::BundleProxObject*)
        pointer to scaling matrix H, may be influenced by a scaling heuristic

    @param[in,out] delta_groundset_minorant (MinorantPointer*)
        if not NULL, the change in groundset aggregate will be stored here

    @param[in,out] delta_index (CH_Matrix_Classes::Indexmatrix*)
        must be not NULL iff delta_groundset_subg!=NULL or yfixed has changed,
        will store nonzero indices of delta_groundset_subg

    @param[in] relprec (CH_Matrix_Classes::Real)
        relative precision for termination in QP computations

    @return 0 on success, != 0 on failure

    */
    virtual int candidate(CH_Matrix_Classes::Integer& gs_id,
      CH_Matrix_Classes::Matrix& newy,
      CH_Matrix_Classes::Real& cand_gs_val,
      CH_Matrix_Classes::Real& linval,
      CH_Matrix_Classes::Real& augval_lb,
      CH_Matrix_Classes::Real& augval_ub,
      CH_Matrix_Classes::Real& subgnorm2,
      const CH_Matrix_Classes::Matrix& center_y,
      CH_Matrix_Classes::Real center_value,
      const MinorantPointer& model_minorant,
      BundleProxObject* Hp,
      MinorantPointer* delta_groundset_minorant = 0,
      CH_Matrix_Classes::Indexmatrix* delta_index = 0,
      CH_Matrix_Classes::Real relprec = 1e-2);

    /// returns the groundset aggregate computed in candidate()
    virtual const MinorantPointer& get_gs_aggregate() const {
      return gs_aggregate;
    }

    /// returns the linear minorant valid on the entire ground set (e.g. a linear cost funciton)
    virtual const MinorantPointer& get_gs_minorant() const {
      return gs_aggregate;
    }

    /// if not NULL (iff get_use_yfixing()==false) it returns the vector yfixed with yfixed(i)=0 if not fixed, =1 is fixed already, =2 if newly fixed
    virtual const CH_Matrix_Classes::Indexmatrix* get_yfixed() const {
      return (use_yfixing ? &yfixed : 0);
    }

    /// if not NULL (iff get_use_yfixing()==false) returns the vector yfixed with yfixed(i)=0 if not fixed, =1 is fixed already, =2 if newly fixed
    virtual CH_Matrix_Classes::Indexmatrix* set_yfixed() {
      return (use_yfixing ? &yfixed : 0);
    }

    /// true if the cooridinate fixing heuristic is switched on (only constrained cases)
    bool get_use_yfixing() const {
      return use_yfixing;
    }

    /// set to true to switch on the cooridinate fixing heuristic (only constrained cases)
    void set_use_yfixing(bool uyf) {
      use_yfixing = uyf;
    }

    /// return a new modification object on the heap that is initialized for modification of *this  
    GroundsetModification* start_modification() {
      return new LPGroundsetModification(dim, 0, this);
    }

    /// change the groundset description as specified by the argument 
    int apply_modification(const GroundsetModification& mdf);

    /// m-file output routine for debugging or testing in Matlab (not yet working)
    int mfile_data(std::ostream& out) const;

    ///output settings
    void set_cbout(const CBout* cb, int incr = -1) {
      VariableMetricModel::set_cbout(cb, incr);
      if (qp_solver)
        qp_solver->set_cbout(this);
    }
  };



  //@}

}

#endif

