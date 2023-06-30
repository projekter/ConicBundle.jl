/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/QPSolverObject.hxx
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


#ifndef CONICBUNDLE_QPSOLVEROBJECT_HXX
#define CONICBUNDLE_QPSOLVEROBJECT_HXX

/**  @file QPSolverObject.hxx
    @brief Header declaring the abstract class ConicBundle::QPSolverObject
    @version 1.0
    @date 2020-03-17
    @author Christoph Helmberg
*/


#include "QPModelDataObject.hxx"

namespace ConicBundle {


  /** @defgroup InternalQPSolverInterface Abstract interface for qp solvers

     @brief The classes define the interface for qp solvers to be compatible with
     ConicBundle::UnconstrainedGroundset and ConicBundle::LPGroundset. Solvers
     obeying these requirements may be installed by
     ConicBundle::MatrixCBSolver::set_qp_solver().

     The model description consists of a cutting model and a quadratic
     prox term.  The cutting model needs to be collected via the
     properties of a QPModelDataPointer, see \ref InternalQPModelData.
     Default quadratic proxterms are passed to the solver in the form of
     a ConicBundle::BundleProxObject, see \ref InternalBundleProxObject,
     but the unspecific QPSolverProxObject opens some room for other
     contexts. Likewise the unspecified QPSolverParameterObject only
     serves as a place holder for a class specificying the paraemers of
     the current solver.

     Implementations are ConicBundle::QPSolver for general a general
     ConicBundle::LPGroundset and ConicBundle::UQPSolver (see \ref
     ConstrainedQPSolver) for ConicBundle::UnconstrainedGroundset as well
     as special prox combinations with a box constrained
     ConicBundle::LPGroundset.
  */

  //@{

  /** @brief in order to pass parameters to a customized QPSolverObject, derive the parameters from this object; no other requirements */

  class QPSolverParametersObject {
  public:
    /// virtual destructor
    virtual ~QPSolverParametersObject();
  };

  /** @brief in order to pass a ConicBundle::BundleProxObject, see \ref InternalBundleProxObject, to a custzomized QPSolverObject, the prox classes are also derived from this object */

  class QPSolverProxObject {
  public:
    /// virtual destructor
    virtual ~QPSolverProxObject();

    ///return true if H is of the form diagonal matrix plus Gram matrix of a low rank matrix
    virtual bool is_DLR() const = 0;

    /// add H to the dense symmetric matrix as a principal submatrix starting at position start_index
    virtual int add_H(CH_Matrix_Classes::Symmatrix& big_sym, CH_Matrix_Classes::Integer start_index = 0) const = 0;

    ///add \f$Hx\f$ to outplusHx and return this
    virtual CH_Matrix_Classes::Matrix& add_Hx(const CH_Matrix_Classes::Matrix& x,
      CH_Matrix_Classes::Matrix& outplusHx,
      CH_Matrix_Classes::Real alpha = 1.) const = 0;

    /** @brief returns \f$H^{-1}x\f$, possibly transformed by some AffineFunctionTransformation

        If a stack of instances of AffineFunctionTransformation has been added
        by push_aft() and A_1 to A_k are on the stack with A_k on top,
        then it returns \f$ A_k\cdots A_1\cdot H^{-1}\cdot A_1^\top\cdots A_k^\top\cdot x\f$
     */
    virtual CH_Matrix_Classes::Matrix& apply_Hinv(CH_Matrix_Classes::Matrix& x) const = 0;

    /** @brief return \f$Diag(D)+VV^T\f$ which either equals $H$ exactly (iff is_DLR() returns true)
        or, if not, serves as an approximation hopefully suitable for
        preconditioning. If $V$ is available it may be large, so return only a
        pointer to avoid copying. If no $V$ is available, then put Vp==0 on
        exit.
    */
    virtual void get_precond(CH_Matrix_Classes::Matrix& D,
      const CH_Matrix_Classes::Matrix*& Vp) const = 0;


    /// may be used to indicate seemingly conservative step sizes possibly due to the quadratic term
    virtual void set_short_QPsteps(CH_Matrix_Classes::Integer shortQPst) = 0;

    /// retrieves the number of conservative step sizes possibly due to the quadratic term passed on to this
    virtual CH_Matrix_Classes::Integer get_short_QPsteps() = 0;

  };


  /** @brief abstract interface for a QPSolver

    In QPsolve() and QPupdate() the QPSolverProxObject gives access to the quadratic term of the cost function, maybe in several different ways, e.g. by  a routine computing the cost terms explicitly or as a matrix times vector interface or as a low rank description with diagonal; the object pointed to just provides information but is neither owned nor deleted by this

  */

  class QPSolverObject : public virtual QPModelDataPointer {
  public:
    /// clear
    virtual void QPclear() = 0;

    ///  default constructor
    QPSolverObject(CBout* cb = 0, int cbinc = -1) :QPModelDataPointer(cb, cbinc) {
    }

    /// virtual destructor
    virtual ~QPSolverObject();

    /// the parameter object passed here will be owned and deleted by *this
    virtual int QPset_parameters(QPSolverParametersObject* params) = 0;

    /** @brief in the case of box constraints it may be worth to fix some variables to their upper or lower bounds; return true if the QPsolver supports this in QPSolve

    */

    virtual bool QPsupports_yfixing() = 0;

    /** @brief return true iff the code supports QPupdate(), i.e., it supports external updates of the groundset aggregate in order to model constraints not included explicitly in the QP's model
    */
    virtual bool QPsupports_updates() = 0;

    /** @brief solve for the data described by @a Hp and
        QPModelDataPointer for the center of stability center_y
        terminating when the relative precision relprec of the objective
        value with respect to lower_bound and upper_bound is reached

        The main termination criterion (besides sufficient feasibility)
        is that the lower bound on the objective value is above
        lower_bound (most of the time this is the most recent augmented
        model value and this should make sure that there is some increase in
        the augemented value) and that it is sufficiently close to the current
        upper bound with relative precision relprec*(current_upper_bound-lower_bound)
        It should never exceed upper_bound (relprec is mostly a pretty large
        number for a relative precision); in practice it may happen
        that for high precision requirements it is not possible to
        get strictly above lower_bound, which indicates the limit of
        achievable numerical precision is reached.

        The groundset aggregated @a gs_aggr is typically only required when
        using an unconstrained QP solver, if simple groundset constraints are
        modelled in a Gauss-Seidel fashion. If QP solvers model all groundset
        constraints explicitly this can be ignored. Otherwise the groundset
        aggregate is updated each time to be the best groundset minorant
        for keeping the candidate inside the feasible set for the current
        model aggregate.

        The parameter @a yfixed holds a nonzero value in each index where
        the center_y value should no longer be changed by the qp solver.
        This simply requires some precprocessing which eliminates the
        corresponding variables but this need not be implemented if
        QPsupports_yfixing() returns false, then the parameter may simply
        be ignored. yfixed(i) has value
        - 0 if the coordinate is not fixed
        - 1 if the coordinate has already been fixed for some time
        - 2 if the coordinate is newly fixed. In this latter case
            QPsolve() should reset it to value 1 as soon as it has take note of this.


     */
    virtual int QPsolve(const CH_Matrix_Classes::Matrix& center_y,
      CH_Matrix_Classes::Real lower_bound,
      CH_Matrix_Classes::Real upper_bound,
      CH_Matrix_Classes::Real relprec,
      QPSolverProxObject* Hp,
      const MinorantPointer& gs_aggr,
      CH_Matrix_Classes::Indexmatrix* yfixed) = 0;

    /** @brief resolve after a call to QPsolve for a modified groundset minorant, whose changes are described in delta_gs_subg for the indices in delta_index

        The routine is only called if the solution of QPsolve does not
        produce a feasible candidate by itself and the candidate has to
        be made feasible by adapting the groundset aggregate @a gs_aggr in
        a suitable way. This never happens if a QP solver includes all
        groundset constraints explitictly. Currently it may happen in
        ConicBundle when the unconstrained QP solver is used together
        with box constraints. If it happens, @a delta_gs_aggr gives the
        change from the previous groundset aggregate to the new
        groundset aggregate as well as a list of indices in
        @a delta_index on which coordinates changed. This information
        should help to efficiently update the QP cost function without
        the need to fully form in particular the quadratic cost matrix.

        @a yfixed has the same meaning as for QPsolve

     */
    virtual int QPupdate(const CH_Matrix_Classes::Matrix& center_y,
      CH_Matrix_Classes::Real lower_bound,
      CH_Matrix_Classes::Real upper_bound,
      CH_Matrix_Classes::Real relprec,
      QPSolverProxObject* Hp,
      const MinorantPointer& gs_aggr,
      CH_Matrix_Classes::Indexmatrix* yfixed,
      const MinorantPointer& delta_gs_aggr,
      const CH_Matrix_Classes::Indexmatrix& delta_index) = 0;

    /** @brief if in the model description some trace/penalty values were adapted, this may require resolving without any other change in information

     */

    virtual int QPresolve(CH_Matrix_Classes::Real lower_bound,
      CH_Matrix_Classes::Real upper_bound,
      CH_Matrix_Classes::Real relprec) = 0;

    /// returns 0 and the solution to the last solving call if the data is available, otherwise it returns 1
    virtual int QPget_solution(CH_Matrix_Classes::Real& augval_lb,
      CH_Matrix_Classes::Real& augval_ub,
      CH_Matrix_Classes::Matrix& new_point,
      CH_Matrix_Classes::Real& gs_aggr_offset,
      CH_Matrix_Classes::Matrix& gs_aggr_gradient) = 0;

    /// returns the last lower bound used for termination
    virtual CH_Matrix_Classes::Real QPget_lower_bound() = 0;

    /// return a new modification object on the heap that is initialized for modification of *this, return 0 if no modifications applicable   
    virtual GroundsetModification* QPstart_modification() = 0;

    /// apply the modification, return 0 if successful and 1 if unsuccessful 
    virtual int QPapply_modification(const GroundsetModification& mdf) = 0;

    /// check whether the point y is feasible with respect to the constraints describing the groundset of the QP for y
    virtual bool QPis_feasible(const CH_Matrix_Classes::Matrix& y,
      CH_Matrix_Classes::Real relprec = 1e-10) = 0;

    /// makes y feasible if it is not feasible for the groundset of the QP, see Groundset::ensure_feasibility()
    virtual int QPensure_feasibility(CH_Matrix_Classes::Matrix& y,
      bool& ychanged,
      QPSolverProxObject* Hp,
      CH_Matrix_Classes::Real relprec = 1e-10) = 0;

    /// returns true if, for the current constraints and the requested ProxObject, it might be better to use the internal unconstrained QP solver (which can deal with box constraints by a work-around)
    virtual bool QPprefer_UQPSolver(QPSolverProxObject*) const = 0;

    /// returns false if the feasible set is the entire space (unconstrained optimization), true otherwise. 
    virtual bool QPconstrained() const = 0;

    /// returns true if and only if there exist box constraints and these are the only constraints; if there are box constraints (and maybe others), the pointers return the reference to them (full dimensional dense vectors lb und ub; lbind and ubind are index vectors giving only the indices i of entries lb(i)>-infty and ub(i)<infty respectively sorted by increasing index values); if the corresponding objects do not exist, the value returned is null 
    virtual bool QPboxconstrained(const CH_Matrix_Classes::Matrix*& lb,
      const CH_Matrix_Classes::Matrix*& ub,
      const CH_Matrix_Classes::Indexmatrix*& lbind,
      const CH_Matrix_Classes::Indexmatrix*& ubind) const = 0;

    /// allows to output some implementation dependent statistics on run time behaviour
    virtual std::ostream& QPprint_statistics(std::ostream& out, int printlevel = 0) = 0;

  };



  //@}

}

#endif

