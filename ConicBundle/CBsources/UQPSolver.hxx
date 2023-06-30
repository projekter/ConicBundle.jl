/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/UQPSolver.hxx
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


#ifndef CONICBUNDLE_UQPSOLVER_HXX
#define CONICBUNDLE_UQPSOLVER_HXX

/**  @file UQPSolver.hxx
    @brief Header declaring the class ConicBundle::UQPSolver
    @version 1.0
    @date 2020-03-18
    @author Christoph Helmberg
*/

#include <vector>
#include "QPSolverObject.hxx"
#include "UQPModelBlock.hxx"
#include "clock.hxx"
#include "Groundset.hxx"

namespace ConicBundle {

/** @defgroup UnconstrainedQPSolver Internal QP Solver for unconstrained groundsets

    If the groundset (see e.g. UnconstrainedGroundset) of the convex
    optimization problem is unconstrained or if it only involves box
    constraints (as a special case of an LPGroundset) towgether with a
    diagonal proximal term (eg BundleIdProx or
    BundleDiagonalTrustRegionProx) an unconstrained QP solver suffices
    to determine the next candidate as the solution of the bundle
    subproblem. Because ConicBundle is designed to handle also
    nonpolyhedral cutting models arising from support functions over
    the second order cone and positive semidefinite cone the general
    solver is an intererior point method. More precisely, it
    implements a primal dual predictor corrector method.

    The main visible object and managing core of the solver is
    UQPSolver.  It serves as an interface for ConicBundle as it
    implements a QPSolverObject (for finding the next candidate via
    solving the quadratic bundle subproblem) and a QPDataPointer (for
    collection the cutting model information).
      
    The cutting models are realized by two implemtented classes
    derived from UQPModelBlock which brings together the abstract
    classes QPModelDataObject (interface for the cutting models of
    ConicBundle) and UQPModelBlockObject (interface for the interior
    point solver):

    - UQPConeModelBlock : it allows to describe any of the cutting models
      employed by the classes NNCModel, SOCModel, PSCModel, BoxModel
   
    - UQPSumModelBlock : it forms the sum of several UQPModelBlock instances
      as required in SumModel

    The efficiency of the UQPSolver in comparison to the full QPSolver
    relies on first eliminating the unconstrained design variables (a
    kind of Schur complement step in the full KKT system of the saddle
    point prpboem) and solving the bundle subproblem in terms of the
    model variables only. Quite often the computation of the QP cost
    terms for this reduced problem is the most costly operation.
    How to do this efficiently depends mostly on the proximal term
    used. Therefore the routines for computing the reduced QP cost
    terms for the UQPSolver are implemented directly within each 
    ProxObject by the routine BundleProxObject::compute_QP_costs(),
    which is called inside UQPSolver::QPsolve() for getting the
    quadratic cost terms. 

    When some of the convex functions are used with a FunctionTask in
    AdaptivePenaltyFunction, a change in the penalty term only affects
    the right hand side of the trace constraint in the models. It does
    not affect the unconstrained QP cost terms and if the cutting
    model is not changed otherwise they need not be recomputed. Thus,
    resolving for adpated penalty terms is particularly efficient and
    the BundleSolver then calls UQPSolver::resolve() instead of
    UQPSolver::QPsolve().

    A major challenge is handling box constraints with an
    unconstrained solver, because the candidate may have to be pushed
    back into the feasible set thereby losing optimality with respect
    to the model representation.  How to resolve this via potentially
    several repeated unconstrained solves with updated groundset
    aggregates in a Gauss-Seidl fashion is described by Helmberg and
    Kiwiel. It works best, if the optimal groundset aggregate pushing
    the candidate into the feasible set can be determined by explicit
    computations as in the case of a diagonal prox term (it also works
    for more general groundsets but then requires solving a separate
    QP over the groundset). If the model quality of the current model
    aggregate is poor for this corrected candidate, the BundleSolver
    may decide to recompute the subproblem with the corrected groundset
    aggregate. The latter only enters the cost terms of the QP.
    Doing the update of these costs efficiently is again
    important. This happens in BundleProxObject::update_QP_costs(),
    which is called by UQPSolver::QPupdate() which itself is called
    when the BundleSolver decides to recompute the bundel
    subproblem for the same cutting model with corrected groundset
    aggregate.    
    
    A feature that is particularly easy to handle in solving 
    box constrained variants by the unconstrained approach is
    the fixing of some variables to the bounds. It simply means
    that in computing the unconstrained QP costs these variables
    may be considered to be constant values. So while the list
    of fixed values is passed through the interface of the
    UQPSolver it actually doesn't have to worry about it other
    than passing it on to the BundleProxObject which then
    takes care of this infomation in computing the QP costs.
*/
//@{

/** @brief unconstrained QP solver combining the properties of a
    QPModelDataPointer and QPSolverObject

    For solving the bundle subproblem the standard entry point is
    QPsolve() (see \ref UnconstrainedQPSolver for QPupdate() or
    QPresolve()).  The cutting model is given by an internal
    UQPModelBlock obtained via being passed to the BundleModel as a
    QPModelDataPointer. The BundleProxObject provides the routines for
    computing the quadratic cost terms for eliminated design
    variables, so this solver is blind to the structure of the
    Groundset. 
     
    After obtaining the quadratic cost function from the
    BundleProxObject, the solver calls solve() (or update() or
    resove()) which intializes the variables also of
    UQPModelBlock, determines
    the starting point and calls iterate(). This routine then
    computes for the current point the
    primal and dual violations, the opimality gap, etc.  As long as
    the termination criterion (steered by upper and lower bounds
    and termination precision) is not met, it computes predictor 
    corrector steps by calling predcorr_step(). 
    The latter sets up and solves a primal dual KKT system 
    in cooperation with the routines of the UQPModelBlock
    (with underlying implementations in UQPConeModelBlock and
    UQPSumModelBlock) and intermediate/susbsequent calls to 
    UQPModelBlock.line_search() and select_new_mu() for 
    determining the next step sizes and barrier parameter.

 */
  
  class UQPSolver: public virtual UQPModelPointer, public virtual QPSolverObject
{
 private:
  /// the current/latest prox term, mainly needed for computing the QP cost matrices and possibly for delivering the solution later in QPget_solution
  BundleProxObject* Hp;
  /// the current/latest center point, mainly needed for computint the QP cost matrices and possibly reconstructing the solution later in QPget_solution
  const CH_Matrix_Classes::Matrix* center_yp;
  /// the current/latest groundset aggeregate, mainly needed for computing the QP cost matrices and possibly reconstructing the solution later in QPget_solution
  MinorantPointer gs_aggr;  
  
  //--- problem desription
  CH_Matrix_Classes::Symmatrix Q;  ///< quadratic cost matrix (positive definite)
  CH_Matrix_Classes::Matrix c;     ///< linear cost matrix
  CH_Matrix_Classes::Real offset;  ///< constant added to objective

  CH_Matrix_Classes::Matrix A;     ///< a (full) blockmatrix mainly for the trace constraints, it is collected by calling model_block
  CH_Matrix_Classes::Matrix b;     ///< right hand side to A, it is collected by calling model_block

  //--- termination parameters
  CH_Matrix_Classes::Real lowerbound;  ///< lower bound on the final objective value
  CH_Matrix_Classes::Real upperbound; ///< upper bound on the final objective value
  CH_Matrix_Classes::Real termeps;   ///< termination precision
  CH_Matrix_Classes::Integer maxiter; ///< upper limit on the number of interior poitn iterations (none if <0)

  //--- global variables
  CH_Matrix_Classes::Matrix x;         ///< current model x as one joint vector
  CH_Matrix_Classes::Matrix y;         ///< current dual variable to constraints
  CH_Matrix_Classes::Matrix dx;        ///< step for x
  CH_Matrix_Classes::Matrix dy;        ///< step for y

  CH_Matrix_Classes::Real primalval;   ///< primal objective value
  CH_Matrix_Classes::Real dualval;     ///< dual objective value

  CH_Matrix_Classes::Real mu;          ///< barrier parameter
  CH_Matrix_Classes::Matrix Qx;        ///< holds Q*x for current x
  CH_Matrix_Classes::Matrix rhs;       ///< right hand side of the KKT system

  bool run_starting_point;         ///< in restarting go back to the starting point

  CH_Matrix_Classes::Integer iter;     ///< counts number o iterations
  CH_Matrix_Classes::Integer status;   ///< termination status of last call to solve or update

  //--- temporary variables, global only for memory managment purposes
  CH_Matrix_Classes::Symmatrix Qplus;   ///< L*L^T factorization of Q+blockdiag
  CH_Matrix_Classes::Matrix LinvAt;     ///< =L^-1*A^T
  CH_Matrix_Classes::Symmatrix sysdy;   ///< system matrix for dy
  CH_Matrix_Classes::Matrix rd;         ///< dual slack rd=c-Qx-At*y (=-z if feasible)
  CH_Matrix_Classes::Matrix xcorr;      ///< correction value for x
  mutable CH_Matrix_Classes::Matrix tmpvec;  ///< temporary vector for reducing reallocations      
  mutable CH_Matrix_Classes::Matrix tmpxvec; ///< temporary vector for reducing reallocations      
  mutable CH_Matrix_Classes::Matrix tmpyvec; ///< temporary vector for reducing reallocations      


  //--- statistics
  CH_Matrix_Classes::Integer sum_iter;    ///< sum over all interior point iterations
  CH_Matrix_Classes::Integer sum_choliter; ///< sum over Cholesky facotrizations 
  CH_Tools::Clock clock; ///< for timing
  CH_Tools::Microseconds sum_choltime;  ///< sum of time spent in Cholesky
  CH_Tools::Microseconds QPcoeff_time;  ///< time spent in computing the QP coefficients
  CH_Tools::Microseconds QPsolve_time; ///< time spent in solving the QP

  
  //--- private routines
  /// chooses the next value of the barrier parameter
  void select_new_mu(const CH_Matrix_Classes::Matrix& dx,
		     const CH_Matrix_Classes::Matrix& dy,
		     const CH_Matrix_Classes::Matrix& rhs_residual);
  /// carry out a predictor corrector step 
  int predcorr_step(CH_Matrix_Classes::Real& alpha);
  /// calls predcorr_step interatively until termination
  int iterate();

 public:
  /// reset to "empty/no" model
  void clear();

  /// reset parameters to default values
  void set_defaults(); 

  /// default constructor
  UQPSolver(CBout* cb=0,int cbinc=-1):UQPModelPointer(cb,cbinc)
  {clear();set_defaults();}
  /// destructor
  ~UQPSolver(){}

  
  //----  proper own routines
  /// reserve memory for this size
  void init_size(CH_Matrix_Classes::Integer maxdim){Q.newsize(maxdim); x.newsize(maxdim,1);}

  /// returns the quadratic cost matrix 
  const CH_Matrix_Classes::Symmatrix& get_Q(void) const {return Q;}
  /// returns the linear cost vector 
  const CH_Matrix_Classes::Matrix& get_c(void) const {return c;}
  /// returns the constant cost offset value 
  CH_Matrix_Classes::Real get_offset(void) const {return offset;}

  /// sets the termination lower and upper bounds 
  void set_termbounds(CH_Matrix_Classes::Real lb,CH_Matrix_Classes::Real ub){lowerbound=lb;upperbound=ub;}
  /// sets the termination precision
  void set_termeps(CH_Matrix_Classes::Real te){termeps=te;}
  /// sets the upper bound on the number of interior point iterations (<0 means no bound)
  void set_maxiter(CH_Matrix_Classes::Integer mi){maxiter=mi;}

  /// return the number of iterations of the last solve
  CH_Matrix_Classes::Integer get_iter() const {return iter;}
  /// return the status of the last solve
  CH_Matrix_Classes::Integer get_status() const {return status;}
  /// return the termination precision
  CH_Matrix_Classes::Real get_termeps() const {return termeps;}
  /// return the upper bound on interior point interations
  CH_Matrix_Classes::Integer get_maxiter() const {return maxiter;}

  /// return the primal objective value (lower bound) of the last solve
  CH_Matrix_Classes::Real get_primalval()const {return primalval;}
  /// return the dual objective value (upper bound) of the last solve
  CH_Matrix_Classes::Real get_dualval()const {return dualval; }

  /// return the joint model vector (primal solution) produced by the last solve
  const CH_Matrix_Classes::Matrix& get_x() {return x;}
  /// return the joint dual vector (dual solution) produced by the last solve
  const CH_Matrix_Classes::Matrix& get_y() {return y;}

  /// solve the QP for this cost function from scratch
  int solve(const CH_Matrix_Classes::Symmatrix& Q,const CH_Matrix_Classes::Matrix& c,CH_Matrix_Classes::Real offset);
  // returns status
  //   0 ... if solve terminated correctly for the standard stopping criterion
  //   1 ... if terminated by maxiter or precision tolerances
  //   2 ... if the factorization of the system failed
  //   3 ... if one of the blocks failed 
  
  /// resolve the QP for the same cost function as last time with slightly modified feasible set
  int resolve();
  // returns status
  //   0 ... if solve terminated correctly for the standard stopping criterion
  //   1 ... if terminated by maxiter or precision tolerances
  //   2 ... if the factorization of the system failed
  //   3 ... if one of the blocks failed 

  /// resolve the QP for the same feasible set but update the cost terms by the given argumetns first  
  int update(const CH_Matrix_Classes::Symmatrix& dQ,const CH_Matrix_Classes::Matrix& dc,CH_Matrix_Classes::Real doffset);
  // solves the system for c+=dc offset+=doffset and returns status
  //   0 ... if solve terminated correctly for the standard stopping criterion
  //   1 ... if terminated by maxiter or precision tolerances
  //   2 ... if the factorization of the system failed
  //   3 ... if one of the blocks failed 

  ///output some statistical information on performance
  std::ostream& print_statistics(std::ostream& out) const;

  ///save the current settings and values
  std::ostream& save(std::ostream& out) const;

  ///restore the settings and values
  std::istream& restore(std::istream& in);

  //------------  QPSolverOject routines

  /// calls clear(), i.e. reinitialize completely
  void QPclear(){clear();}

  /// does nothing here
  int QPset_parameters(QPSolverParametersObject* /* params */)
  {
    //not used here
    return 0;
  }


  /// does nothing here (unconstrained case)
  int apply_modification(const GroundsetModification& /* mdf */)
  {return 0;}

  /// no difficulty if the BundleProxObject does it
  bool QPsupports_yfixing()
  { return true; }

  /** @brief return true iff the code supports QPupdate(), i.e., it supports external updates of the groundset aggregate in order to model constraints not included explicitly in the QP's model   
  */
  virtual bool QPsupports_updates()
  { return true;}

  /// return the lower bound on the objective value of the bundle subproblem
  CH_Matrix_Classes::Real QPget_lower_bound()
  { return get_primalval(); }

  /// see QPSolverObject::QPsolve() and  \ref UnconstrainedQPSolver 
  int QPsolve(const CH_Matrix_Classes::Matrix& center_y,
	      CH_Matrix_Classes::Real lower_bound,
	      CH_Matrix_Classes::Real upper_bound,
	      CH_Matrix_Classes::Real relprec,
	      QPSolverProxObject* Hp,
	      const MinorantPointer& gs_aggr,
	      CH_Matrix_Classes::Indexmatrix* yfixed);

  /// see QPSolverObject::QPupdate() and  \ref UnconstrainedQPSolver 
  int QPupdate(const CH_Matrix_Classes::Matrix& center_y,
	       CH_Matrix_Classes::Real lower_bound,
	       CH_Matrix_Classes::Real upper_bound,
	       CH_Matrix_Classes::Real relprec,
	       QPSolverProxObject* Hp,
	       const MinorantPointer& gs_aggr,
	       CH_Matrix_Classes::Indexmatrix* yfixed,
	       const MinorantPointer& delta_gs_subg,
	       const CH_Matrix_Classes::Indexmatrix& delta_index);

  /// see QPSolverObject::QPresolve() and  \ref UnconstrainedQPSolver 
  int QPresolve(CH_Matrix_Classes::Real lower_bound,
		CH_Matrix_Classes::Real upper_bound,
		CH_Matrix_Classes::Real relprec);

  /// the unconstrained solver can only provide this information if the references to the input data of QPsolve() and QPudate() are still available and unchanged, otherwise the behavior is undefined and will hopefully return 1 if not valid
  int QPget_solution(CH_Matrix_Classes::Real&  augval_lb,
		     CH_Matrix_Classes::Real&  augval_ub,
		     CH_Matrix_Classes::Matrix& new_point,
		     CH_Matrix_Classes::Real&  gsaggr_offset,
		     CH_Matrix_Classes::Matrix& gsaggr_gradient);

  ///as this IS the internal UQPSolver, the answer doesn't matter, but the solver would certainly say yes
  bool QPprefer_UQPSolver(QPSolverProxObject*) const
  { return true; }  
  
  /// it is always unconstrained
  bool QPconstrained() const
  {return false;}

  /// returns 0 because no modifications are applicable   
  virtual GroundsetModification* QPstart_modification()
  {return 0;}

  /// no modifications need to be carried out here as there is no data to be modified, so any modification succeeds 
  int QPapply_modification(const GroundsetModification& /* mdf */)
  {return 0;}

  /// for the unconstrained solver any point is feasible, because it cannot even check the dimension of the design space
  bool QPis_feasible(const CH_Matrix_Classes::Matrix& /* y */,
		     CH_Matrix_Classes::Real /* relprec */=1e-10)
  {return true;}

  /// for the unconstrained solver any point is feasible, because it cannot even check the dimension of the design space, so any y is left unchanged
  virtual int QPensure_feasibility(CH_Matrix_Classes::Matrix& /* y */,
				   bool& ychanged,
				   QPSolverProxObject* /* Hp */,
				   CH_Matrix_Classes::Real /* relprec */=1e-10)
  { ychanged=false; return 0;}

  /// the unconstrained solver does not have this inforamion and does not know about it, so it returns NULL pointers and false
  bool QPboxconstrained(const CH_Matrix_Classes::Matrix*& lb,
			const CH_Matrix_Classes::Matrix*& ub,
			const CH_Matrix_Classes::Indexmatrix*& lbind,
			const CH_Matrix_Classes::Indexmatrix*& ubind) const
  { lb=ub=0; lbind=ubind=0; return false;}

  /// outputs some statistical data about solver performance
    std::ostream& QPprint_statistics(std::ostream& out,int /* printlevel */ =0)
  {return print_statistics(out);}
  

};

  //@}

  
}

#endif

