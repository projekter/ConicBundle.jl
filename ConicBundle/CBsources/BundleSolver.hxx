/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/BundleSolver.hxx
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



#ifndef CONICBUNDLE_BUNDLESOLVER_HXX
#define CONICBUNDLE_BUNDLESOLVER_HXX


/**  @file BundleSolver.hxx
    @brief Header declaring the class ConicBundle::BundleSolver
    @version 1.0
    @date 2014-08-05
    @author Christoph Helmberg
*/


#include "clock.hxx"
#include "QPSolverObject.hxx"
#include "Groundset.hxx"
#include "BundleModel.hxx"
#include "BundleWeight.hxx"
#include "BundleTerminator.hxx"


namespace ConicBundle {

/** @defgroup InternalBundleSolver Internal Bundle Solver of ConicBundle
   @brief Solve \f$min_{y\in\mathbf{R}^m}  f(y)\f$ 
    for a convex function f given in the form of a ConicBundle::BundleModel.
    The Problem and the Solver are set up and called from
    one of the solver interfaces, e.g., ConicBundle::CBSolver or 
    ConicBundle::MatrixCBSolver and are typically not accessed directly.

*/
//@{


/** @brief This is the internal bundle solver managing descent/null steps
    with inner updates for approximating the proximal point by alternatingly
    adapting the groundset aggregate and the model aggregate. 

    The main input paramter needed is a BundleModel representing the cutting
    model of the convex function. By default the Groundset is the entire space
    (unconstrained), but other groundsets can be used, see LPGroundset for an
    example. Each groundset needs to provide a suitable, QPSolver

    The default proximal term is the squared Euclidean norm, but convergence
    may be improved considerably if a better convex quadratic term is
    available, which can then be supplied via the class BundleProxObject. There
    is also a variable metric heuristic that manipulates the quadratic term 
    but it will only help sometimes. In addition, the influence
    of the quadratic term is steered by a BundleWeight for which there is a
    rather useful dynamic default heuristic somewhat similar to guiding a
    trust region radius. The quadratic term and the weight together determine
    the norm that is used to measure the "step length", its dual norm is used
    to measure the norm of the combined aggregate linear minorant (groundset 
    and model aggregate).

    The termination criterion is supplied by the class BundleTerminator, but
    the BundleSolver may also directly be forced to stop after a certain
    number of oracle evaluations or every descent step.

    The main work is solving the quadratic bundle subproblems for determining
    the next groundset aggregate and the next model aggregate. Because these
    may profit from exploiting knowledge about current structural properties
    these subproblems will be solved by a QPSolver of the GroundSet which
    provides a general interface for supplying model information.
    This leaves the rather trivial tasks of calling these
    routines in the right sequence to the BundleSolver. 

    As ConicBundle is designed for modifiying the functions on the fly,
    a lot of bookkeeping is needed to keep track of changes. The basic
    idea is to keep track of versions by some kind of update counters.
    The BundleSolver generates the sequence of points and therefore
    has an integer identification number for each point it generates.
    Whenever a point is modified it gets a new id. Similarly the 
    groundset as well as the function each have a separate id that
    has to be increased whenever it is modified, so that it can be
    recognized easily whether certain points or evaluations are
    still valid by checking whether the corresponding id values still 
    match. Likewise the BundleModel has to assign a new id to
    each new model aggregate minorant it generates.

    There is a fundamental difference between groundset aggregate and
    model aggregate. A valid groundset aggregate is always easily available,
    it is simply the constant zero function (potentially replaced by a linear
    cost function), so the groundset is always 
    assumed to have a valid groundset aggregate available, which is usually
    the one stemming from the last ground set quadratic subproblem.
    For the function model it is much more difficult to provide and
    keep a valid model aggregate, in particular in the presence of
    function changes, so a local copy of the last version of this model 
    aggregate is kept in the BundleSolver, which also explains the need
    of a separate id in order to recognize its validity.

    Changes to the groundset as well as the functions are always 
    assumed to be initiated by BundleSolver::apply_modification()
    which then calls Groundset::apply_modfication() and
    BundleModel::apply_modification(), which finally
    propagates the modification to the oracles.  

*/


class BundleSolver: public BundleTerminatorData, public CBout
{
private:
  //-----------------------------------------------------------------------
  /** @name Model, Groundset, QP solver and customizable parameters and parameter classes
   */
  //@{

  /// The pointer to the cutting model of the function to be optimized; the model is always external and will not be deleted
  BundleModel *model;       

  /// pointer to the groundset information, tpyically occupied by basic internal variants that are deleted on destruction
  Groundset* groundset;
  /// if true the groundset is not deleted on destructoin. The variable is set to true if a groundset is implanted from outside 
  bool external_groundset;

  /// pointer to the class providing the termination criterion, will be deleted before replacement or on destruction
  BundleTerminator *terminator; 

  /// pointer to the class providing the routine for choosing the weight u for the quadratic proximal term, will be deleted before replacement or on destruction
  BundleWeight *bundleweight;  

  
  /// pointer to an external clock for timing statistics, not deleted on destruction
  const CH_Tools::Clock* clockp;  

  //@}

  //-----------------------------------------------------------------------
  /** @name  the following parameters are all reset by set_defaults
   */

  //@{

  /// acceptance factor for descent steps (default = 0.1, no need to change this)
  CH_Matrix_Classes::Real mL;             
  /// nullstep factor for inexact evaluation >= mL (default = 0.15, no need to change this)
  CH_Matrix_Classes::Real mN;
  /// by default true, it then uses the aggregate linear minorant for the acceptance criterion, if false cutval (the model value) is used instead            
  bool use_linval;      
  /// fixed parameter for the precision of the model before attempting a step (default = 0.05, not much use to change this)
  CH_Matrix_Classes::Real modeleps;       
  /// maximum number of inner Gauss-Seidel steps before attempting a step, use a negative value for no bound, default = 10
  CH_Matrix_Classes::Integer max_updates; 

  /// pointer to the quadratic proximal term; will be deleted before replacement or on destruction 
  BundleProxObject* Hp;

  /** @brief employ a variable metric heuristic (0 none, 1 diagonal scaling, 2 diagonal scling plus bounds scaling in the GroundSet, 3 lowrank scaling)
   */
  int do_variable_metric; 
  
  /// default true checks validity of new subgradients against center value, set to false to switch this off 
  bool check_center_validity;


  //@}

  //-----------------------------------------------------------------------
  /** @name  center, candidate, function values and aggregates
   */

  //@{
  
  /// increased whenever a candidate point is formed or a center point is modified (-1 if not initialized, >=0 otherwise)
  CH_Matrix_Classes::Integer point_id; 

  /// point_id of the center of stability
  CH_Matrix_Classes::Integer center_id;  
  /// center of stability 
  CH_Matrix_Classes::Matrix center_y; 
  /// function version id the center upper bound center_ub was computed for
  CH_Matrix_Classes::Integer center_fid; 
  /// upper bound on the function value in the center
  CH_Matrix_Classes::Real center_ub;     
  /// relativer precision of the upper bound function evaluation in the center
  CH_Matrix_Classes::Real center_relprec; 
  /// groundset id for which validity of the center was checked
  CH_Matrix_Classes::Integer center_gid; 
  /// value of the groundset cost function in cand_y for center_gid
  CH_Matrix_Classes::Real center_gs_val;


  /// point_id of the most recent candidate
  CH_Matrix_Classes::Integer cand_id;     
  /// the most recent candidate
  CH_Matrix_Classes::Matrix cand_y;  
  /// function version id the candidate upper bound cand_ub was computed for
  CH_Matrix_Classes::Integer cand_fid;    
  /// upper bound on the function value in the candidate
  CH_Matrix_Classes::Real cand_ub;        
  /// relativer precision of the upper bound function evaluation in the candidate
  CH_Matrix_Classes::Real cand_relprec;   
  /// groundset id for which validity of the candidate was checked
  CH_Matrix_Classes::Integer cand_gid;
  /// value of the groundset cost functon in cand_y for cand_gid
  CH_Matrix_Classes::Real cand_gs_val;
 
  /// the most recent model aggregate linear minorant  
  MinorantPointer model_aggregate;  
  /// the aggregate id assigned by the BundleModel when this model aggregate was generated 
  CH_Matrix_Classes::Integer model_aggregate_id; 
 


  //@}

  //-----------------------------------------------------------------------
  /** @name  augmented and model values, central algorithmic values
   */

  //@{

  /// lower bound on value of the last augmented model (lower bound, has to increase throughout)
  CH_Matrix_Classes::Real augval_lb;
  /// upper bound on value of the last augmented model (lower bound, has to increase throughout)
  CH_Matrix_Classes::Real augval_ub;
  /// value of the aggregate linear minorant in the current candidate  
  CH_Matrix_Classes::Real linval;
  /// value of the cutting model in the current candidate
  CH_Matrix_Classes::Real cutval;
  /// if use_linval==true this is linval, otherwise cutval
  CH_Matrix_Classes::Real modelval;
  /// modelprec=(cutval-linval)/max(center_ub-linval,1e-16);
  CH_Matrix_Classes::Real modelprec;

  /// the current weight for the proximal term
  CH_Matrix_Classes::Real weightu;
  /// the squared dual norm of the aggregate (with respect to the quadratic proximal term)
  CH_Matrix_Classes::Real aggr_dnormsqr;    

  /// true if the last completed iteration was a descent step (if false, null_step need not be true!, e.g. termination might occur before either took place)  
  bool descent_step;        
  /// true if the last completed iteration was a null step (if false, descent_step need not be true!, e.g. termination might occur before either took place)  
  bool null_step;           

  /// value returned by the model's quadratic bundle subproblem solver
  int retcode;
  /// value resulting from the last call to BundleTerminator::check_termination() via *terminator
  int terminate;               

   
  //@}

  //-----------------------------------------------------------------------
  /** @name  QP_solver and timing aspects for it 
   */

  //@{

  /// pointer to the qp_solver provided and owned by groundset, 
  QPSolverObject *qp_solver;  

  /// if true the qp_solver does not include the groundset constraints and they have to be take care of in a Gauss Seidel fashion
  bool qp_solves_model_without_gs;

  CH_Tools::Microseconds QPcoeff_time;  ///< time spent in computing the QP coefficients
  CH_Tools::Microseconds QPsolve_time; ///< time spent in solving the QP
  CH_Tools::Microseconds make_aggr_time; ///< time spent in forming the aggregate
  CH_Tools::Microseconds evalaugmodel_time; ///< total time spent in solving the quadratic bundle subproblem
  mutable int qp_save_cnt;

  //-----------------------------------------------------------------------
  /** @name  variables for organizational and statistical purposes
   */

  //@{

  /// set to true if function, groundset, proxima term or the like might have been changed since the last iteration
  bool might_be_modified;
  /// set to true if the aggregate might have changed and a first aggregate and a lower bound for augval have to be determined before resuming
  bool initialize_model;
  /// set to true if the subgradient in the candidate should be used as first aggregate
  bool use_cand_for_aggregate;

  /// number of Gauss-Seidel updates since last null or descent step
  CH_Matrix_Classes::Integer updatecnt;
  /// sum over all Gauss-Seidel updates since last clear()
  CH_Matrix_Classes::Integer sumupdatecnt;
  /// number of bundle solver iterations (= null steps) since last descent step
  CH_Matrix_Classes::Integer innerit;    
  /// sum of all bundle solver iterations since last clear()
  CH_Matrix_Classes::Integer suminnerit; 
  /// number of calls to the evaluation oracle by the bundle solver (since last clear())
  CH_Matrix_Classes::Integer cntobjeval; 
  /// number of recomputations of the value in the center (due to errors or numerical problems, since last descent step)
  CH_Matrix_Classes::Integer recomp;     
  /// sum of all recomputatoins since last clear() or clear_fails()
  CH_Matrix_Classes::Integer sumrecomp;  //sum over all recomputations
  /// number of failures in calling a qp-subproblem since last descent/null step
  CH_Matrix_Classes::Integer qpfails;      //counts number of fails in calls to eval_aug_model
  /// sum of all qp failures since last clear() or clear_fails() 
  CH_Matrix_Classes::Integer sumqpfails;   
  /// number of failures in calls to BundleModel::eval_model since last descent/null step
  CH_Matrix_Classes::Integer modelfails; 
  /// sum of all model failures since last clear() or clear_fails()
  CH_Matrix_Classes::Integer summodelfails;//sums up number of fails
  /// number of iterations which failed to increase augval by a realtive epsilon since last descent step
  CH_Matrix_Classes::Integer augvalfails; 
  /// sum over all augvalfails since last clear() or clear_fails()
  CH_Matrix_Classes::Integer sumaugvalfails;  
  /// number of failures in calling BundleMethod::eval_function() (the oracle) since last descent step
  CH_Matrix_Classes::Integer oraclefails; 
  /// sum over all oraclefails since last clear() or clear_fails()
  CH_Matrix_Classes::Integer sumoraclefails;//sums up number of fails

  /// number of descent steps since last clear()
  CH_Matrix_Classes::Integer descent_steps;        //counts number of descent_steps
  /// number of oracle evaluations that yield a subgradient that does not improve the model much in the candidate (mostly due to the approximation to the quadratic subproblem being too rough)
  CH_Matrix_Classes::Integer shallowcut;      

  //@}

  //-----------------------------------------------------------------------
  /** @name  internal subroutines
   */

  //@{
  /// called if do_variable_metric >0, may change the BundleProxObject and exploit the model_aggregate
  int variable_metric(const CH_Matrix_Classes::Indexmatrix* new_indices=0);

  /// performs Gauss-Seidel steps for updating model and groundset aggregate by qp subproblems until a given model precision or update limit is reached 
  int solve_model(void);

  /** @brief Evaluates the augmented model with respect to the center of stability.

Let \f$Y\f$ denote the (convex) feasible ground set and \f$i_Y\f$ its
indicator function, let \f$\hat y\f$ be the center of stability given
by @a center_y, let \f$ (\gamma,g) \f$ be the aggregate ground set
minorant \f$\gamma+g^\top y\le i_Y(y)\f$ described by @a
gs_subg_offset and @a gs_subg, let \f$H\f$ denote the positive
definite scaling matrix with weight \f$u\f$ given by @a Hp, and
let $f$ denote the objective for which this class holds the model
\f$W\subseteq\{(\sigma,s)\colon \sigma+s^\top y\le f(y)\ \forall y\in Y\}\f$,
which is assumed to be convex and compact.
The routine computes the saddle value of

\f[ \max_{(\sigma,s)\in W}\min_{y\in\mathbf{R}^n} [(s+g)^\top y + \gamma+\sigma +\frac{u}2\|y-\hat y\|_H^2],\f]

Of the corresponding saddle point \f$y\f$ and \f$(\sigma,s)\f$
only the aggregate minorant \f$(\sigma,s)\f$ needs to be computed 
explicitly and will be requested later by a separate call to
get_model_aggregate().
 
On input @a augval_lb gives a lower bound on the saddle value 
that has to be reached or exceeded in the optimization process,
while @a center_ub gives an upper bound and is in fact the upper
bound on the function value in @a center_y as returned by 
eval_function() or recompute_center(). The required relative 
precision is given by \a relprec and is relative to the
gap @a center_ub-@a augval_lb. More precisely, if the code implements
a feasible primal dual method where ub is the upper bound (and 
already less or
equal to @a center_ub to machine precision) and
lb is the lower bound (already greater or equal to @a augval_lb to machine
precision), the method may stop if ub-lb<=min(.5*(ub-@a augval_lb),relprec*(@a center_ub-lb)).

The computed lower bound on the saddle value is then returned in @a augval_lb
(and should exceed the input value, otherwise numerical difficulties
are assumed by the bundle solver). 

Note that @a center_ub_fid and @a center_ub are also call be reference
values, because in some cases it may be useful to modify the
objective function during this computation! Indeed, if the objective
contains a penalty term for some constraint that cannot be included
easily in the ground set implementation, it may be necessary to
increase this penalty term if, e.g., the next candidate turns out no to
improve on the violation of this constraint. In this case it is
allowed to increase the value of @a center_ub and to inform the
caller about this general change by also increasing @a center_id.

@param[in,out] augval_lb     
    - on input lower bound on saddle value, 
    - on output increase lower bound of required precision       

@param[in,out] center_ub_fid
     always gives the modification id of the function for which the value @a ub
     was computed, on input as well as on output (may change!)

@param[in,out] center_ub  
    - on input upper bound on saddle value
    - on output possibly increased function value in center_y due 
      to penalty components
     
@param[in] relprec 
    relative precision requirements

@param[in] center_gs_value
    the value of the groundste aggregate minorant in center_y  

@return
    - 0 ... if all is ok, use get_augmodel_sol() for retrieving the solution
    - -1 ... if setting up the QP failed
    - otherwise it returns the satus returned by the internal QP_Solver

@see BundleProxObject::compute_QP_costs()

  */

  int eval_augmodel(CH_Matrix_Classes::Real& augval_lb,
		    CH_Matrix_Classes::Integer& center_ub_fid, 
		    CH_Matrix_Classes::Real& center_ub,
		    CH_Matrix_Classes::Real relprec,
		    CH_Matrix_Classes::Real center_gs_value);
	     
  /** @brief reevaluate the augmented model for updated groundset aggregate w.r.t. the previously called eval_augmodel()
  
ConicBundle determines the candidate by alternatingly evaluating the
mode for a new model aggregate and then computing an adapted groundset
aggregate by calling Groundset::candidate(). Shifting the candidate
by adapting the groundset aggregate might result in a candidate having
a poor model value and the augmented model is then reevaluated by
calling this routine, which should try to make use of the fact that
a similar problem has been solved already (either in eval_augmodel() or
in reeval_augmodel() if several corrections are needed) 
for a slightly different groundset aggregate. The only difference 
to eval_augmodel() is therefore that the change in the groundset aggregate 
is given with respect to the previous call and that these should be
used to update the cost coefficients of the convex QP efficiently
via the routine BundleProxObject::update_QP_costs().

The update of the groundset aggregate \f$g\f$ is given by @a
delta_gs_subg and @a delta_index, more precisely,
old_gs_subg[delta_index[i]]=gs_subg[delta_index[i]]-delta_gs_subg[i],
where the current groundset aggregate @a gs_subg is the input
paramter. @a delta_gs_subg_offset gives the change of \f$\gamma\f$.

@param[in,out] augval_lb 
    - on input lower bound on saddle value
    - on output increased lower bound of required precision       

@param[in,out] center_ub_fid
     always gives the modification id of the function for which the value @a ub
     was computed, on input as well as on output (may change!)

@param[in,out] center_ub  
    - on input upper bound on saddle value
    - on output possibly increased function value in center_y 
      due to penalty components

@param[in] relprec 
    relative precision requirements

@param[in] center_gs_value
    the value of the groundste aggregate minorant in center_y  

@param[in] delta_gs_minorant
    this is the delta added to the previous value of the groundset 
    minorant to get the current aggregate of the ground set

@param[in,out] delta_index 
    must be not NULL iff delta_groundset_subg!=NULL, store nonzero 
    indices of delta_groundset_subg and chang
 
@return
    - 0 ... if all is ok, use get_augmodel_sol() for retrieving the solution
    - -1 ... if setting up the QP failed
    - otherwise it returns the satus returned by the internal QP_Solver

@see eval_augmodel(), BundleProxObject::update_QP_costs()

  */

  int reeval_augmodel(CH_Matrix_Classes::Real& augval_lb,
		      CH_Matrix_Classes::Integer& center_ub_fid, 
		      CH_Matrix_Classes::Real& center_ub,
		      CH_Matrix_Classes::Real relprec,
		      CH_Matrix_Classes::Real center_gs_value,
		      const MinorantPointer& delta_gs_minorant,
		      const CH_Matrix_Classes::Indexmatrix& delta_index);

  //@}
    
public:

  //-----------------------------------------------------------------------
  /** @name  Intialization, defaults and clear 
   */

  //@{

  /// calls clear() 
  BundleSolver(const CBout* cbo=0,int incr=-1);

  /// calls initialize(CH_Matrix_Classes::Integer, BundleModel*)
  BundleSolver(CH_Matrix_Classes::Integer dim, BundleModel* bp=0,const CBout* cbo=0,int incr=-1);

  /// calls initialize(Groundset*,BundleModel*)
  BundleSolver(Groundset* gs,BundleModel* bp=0,const CBout* cbo=0,int incr=-1);

  /// deletes terminator, bundleweight and Hp and, if external_groundset==false, also groundset
  ~BundleSolver();
        
  /// resets all parameters to default values and calls BundleTerminator::set_defaults() for *terminator and BundleWeight::set_defaults() for *bundleweight
  void set_defaults();  

  /// resets all variables and pointers to classes to initial state and calls set_defaults()
  void clear(); 

  /// resets all fail counts to zero (call this to resume computations
  void clear_fails();  

  /// calls clear(), initializes an unconstrained groundset to this dimension and sets the bundle model to bp
  int initialize(CH_Matrix_Classes::Integer dim, BundleModel* bp=0);

  /// calls clear(), initializes the groundset to gs and the bundle model to bp
  int initialize(Groundset* gs,BundleModel* bp=0);

  /// set output and outputlevel of warnings and errors recursively, see CBout
  void set_out(std::ostream* o=0,int pril=1)
  {
    CBout::set_out(o,pril);
  }

  /// set output and outputlevel of warnings and errors recursively with CBout
  void set_cbout(const CBout* cb,int incr)
  {
    CBout::set_out(cb->get_out_ptr(),cb->get_print_level()+incr);
  }
 //@}

  //-----------------------------------------------------------------------
  /** @name  solve and similar main routines for user interaction
   */

  //@{

  /// set/change the model that should be optimized over (for the existing groundset and starting point) 
  int set_model(BundleModel* bp);

  /// replace the current center by *yp or, if yp==0, by the default starting point
  int set_new_center(const CH_Matrix_Classes::Matrix* yp=0);  // yp==0 uses default starting point

  /// execute at most maxsteps iterations of the bundle method stopping before if termination occurs or stop_at_descent_steps==true and a descent step occurs; maxsteps<=0 indicates no bound on the steps
  int solve(int maxsteps=0,bool stop_at_descent_steps=false);
  
  /** @brief print a one line summary about the current state of progress of the algorithm to @a out

      Assuming that an external clock has been set, the line looks as follows
      
      hh:mm:ss.hh endit dd ii uu weightu aggr_dnorm modelval center_ub (N)
      
      These are
      - hh:mm:ss.hh hours, mintues, seconds, hundreth of secondes
      - "endit" this string is short for "end of iteration" and is 
        convenient for the unix command grep; in the case the code 
        has terminated it shows "_endit" instead
      - dd gives the number of descent steps as in @a descent_steps
      - ii gives the total number of iterations of the bundle method, 
        counting descent steps and null steps
      - uu gives the total number of quadratic bundle subproblem 
        evaluations as counted by @a sumupdatecnt 
      - weight gives the value of the weight factor for the proximal 
        term used in the last quadratic bundle subproblem
      - aggr_dnorm gives the dual norm of the aggregate (dual with 
        respect to the quadratic proximal term)
      - modelval gives the model value used for deciding on null or 
        descent step; this is @a linval if @a use_linval ==true and
	@a cutval otherwise
      - center_ub gives the upper bound computed by the oracle on the 
         function value in the center
      - "N" is only shown if the code has not terminated and a null
         step just occured before returning (null_step is true)      
   */
  std::ostream& print_line_summary(std::ostream& out) const;

  /// modify the groundset and the oracle function(s) as described by GroundsetModification and FunObjModMap
  int apply_modification(const GroundsetModification& gsmdf,const FunObjModMap& funmdfmap);

  /// modify the groundset as described by GroundsetModification and inform the oracle function(s) about this change (calls the other apply_modification for an empty FunObjModMap)
  int apply_modification(const GroundsetModification& gsmdf)
  {
    FunObjModMap funmdfmap; 
    return apply_modification(gsmdf,funmdfmap);
  }

  
  //@}

  //-----------------------------------------------------------------------
  /** @name  set algorithmic parameters
   */

  //@{

  /// replace the previous BundleTerminator by bt; bt will be deleted when replaced or on destruction of this
  void set_terminator(BundleTerminator* bt)
  {delete terminator; terminator=bt;}

  /// replace the previous BundleWeight by bw; bw will be deleted when replaced or on destruction of this
  void set_bundleweight(BundleWeight* bw);

  /// set the required model precision by @a in_eps (if it is positive)
  void set_modeleps(CH_Matrix_Classes::Real in_eps)
  {if (in_eps>0.) modeleps=in_eps;}

  /// set the acceptance level for descent steps (rather don't change this!)
  void set_mL(CH_Matrix_Classes::Real in_mL)
  {if ((in_mL>0)&&(in_mL<1.)) mL=in_mL;}

  /// set the acceptance level for null steps (mL<=in_mN<1., rather don't change this!)
  void set_mN(CH_Matrix_Classes::Real in_mN)
  {if ((in_mN>=mL)&&(in_mN<1.)) mN=in_mN;}

  /// if set to true, the value of the aggregate in the candidate is used for deciding on null or descent step, otherwise the model value
  void set_use_linval(bool ul)
  {use_linval=ul;}

  /// if set to true, the groundset may use a heuristic to decide whether a variable is fixed to one of its bounds (often helps to reduce inner update iterations)
  void set_do_yfixing(bool dofix){groundset->set_use_yfixing(dofix);}

  /// 0 ... use no scaling, 1 ... use a scaling heuristic, 2 ... also allow groundset to influence the scaling so as to favor feasibility, 
  int set_variable_metric(int ds);  
  
  /// set the prox term to the given diagonal matrix
  int set_prox_diagonal(const CH_Matrix_Classes::Matrix& insc);

  /// set the proximal term by replacing BundleProxObject with bsp, the latter is deleted when replaced or on destruction of this 
  int set_prox(BundleProxObject* bsp);

  /// set the external clock to be used for output
  void set_clock(const CH_Tools::Clock& myclock)
  {clockp=&myclock;}

  /// set the maximum number of Gauss-Seidel iterations until the next evaluations for descent/null step, use negative numbers for infinite, 0 or 1 for at most 1
  void set_max_updates(CH_Matrix_Classes::Integer mu)
  {max_updates=mu;}
    
  //@}

  //-----------------------------------------------------------------------
  /** @name  get results of the computations
   */

  //@{

  /// returns the value of the last call to BundleTerminator::check_termination()
  int get_terminate() const 
  {return terminate;}

  /// returns the upper bound on the objective value in @a center_y  
  CH_Matrix_Classes::Real get_center_objval() const 
  {return center_ub+center_gs_val;}
  /// returns the upper bound on the objective in @a center_y returned by the oracle 
  CH_Matrix_Classes::Real get_center_ub() const 
  {return center_ub;}
  /// returns the groundset objective in @a center_y r 
  CH_Matrix_Classes::Real get_center_gs_val() const 
  {return center_gs_val;}
  /// returns the current center of stability @a center_y (after a descent step this is the same as the candidate)
  const CH_Matrix_Classes::Matrix& get_center_y() const 
  {return center_y;} 

  /// returns the upper bound on the objective in @a cand_y  
  CH_Matrix_Classes::Real get_cand_objval() const 
  {return cand_ub+cand_gs_val;}
  /// returns the upper bound on the objective in @a cand_y returned by the oracle 
  CH_Matrix_Classes::Real get_cand_ub() const 
  {return cand_ub;}
  /// returns the groundset objective value in @a cand_y 
  CH_Matrix_Classes::Real get_cand_gs_val() const 
  {return cand_gs_val;}
  /// returns the most recent candidate @a cand_y (after a descent step this is the same as the candidate)
  const CH_Matrix_Classes::Matrix& get_cand_y() const 
  {return cand_y;} 

  /// returns the dual norm squared of the current aggregate (dual w.r.t. the quadratic proximal term) 
  CH_Matrix_Classes::Real get_aggr_dnormsqr() const 
  {return aggr_dnormsqr;}
  /// returns the offset of the current aggregate linear minorant (should be called before any modifications, otherwise this may no longer be correct)
  CH_Matrix_Classes::Real get_aggregate_offset() const 
  {return groundset->get_gs_aggregate().offset()+model_aggregate.offset();}
  /// returns the linear term of the current aggregate linear minorant (this should be called before any modifications, otherwise this may no longer be correct or may even cause an error)
  void get_aggregate(CH_Matrix_Classes::Matrix& aggregate) const 
  {
    CH_Matrix_Classes::Real dummy;
    aggregate.init(groundset->get_dim(),1,0.); 
    model_aggregate.get_minorant(dummy,aggregate,0,1.,false);
    groundset->get_gs_aggregate().get_minorant(dummy,aggregate,0,1.,true);
  }

  /// returns the linear term of the current groundset aggregate linear minorant 
  const MinorantPointer& get_gs_aggregate() const 
  {return groundset->get_gs_aggregate();}
 
  /// returns the linear term of the latest model aggregate linear minorant 
  const MinorantPointer& get_model_aggregate() const 
  {return model_aggregate;} 
  
  /// returns the model value in the candidate that was used for deciding on null/descent step  
  CH_Matrix_Classes::Real get_modelval() const 
  {return modelval;}

  /// returns the weight for the proximal term used in the last quadratic subproblem 
  CH_Matrix_Classes::Real get_weight() const {return Hp->get_weightu();}
  /// returns the correction factor used in the termination criterion to compensate the strength of the proximal term 
  CH_Matrix_Classes::Real get_term_corr() const {return Hp->get_term_corr();}

  /// returns true if the latest iteration resulted in a descent step (note, get_descent_step() and get_null_step() may both return false e.g. if termination occurs)
  bool get_descent_step() const {return descent_step;}
  /// returns true if the latest iteration resulted in a null step (note, get_descent_step() and get_null_step() may both return false e.g. if termination occurs)
  bool get_null_step() const {return null_step;}

  /// if the groundset has constraints and set_do_yfixing was set with true enty i of the returned matrix is !=0 if the coordinate was fixed at one of its bounds and 0 if the coordinate is still free
  const CH_Matrix_Classes::Indexmatrix* get_yfixed() const {return groundset->get_yfixed();}


  //@}

  //-----------------------------------------------------------------------
  /** @name  get algorithmic and statistical parameters
   */

  //@{


  /// returns the model precision @a modeleps 
  CH_Matrix_Classes::Real get_modeleps() const {return modeleps;}
  /// returns true if the aggregate linear minorant is used for the model value
  int get_use_linval() const {return use_linval;}
  /// returns true if the proximal term is not of the type BundleIdProx, i.e. if it is not simply the squared Euclidean norm
  int get_do_variable_metric() const {return do_variable_metric;}
  /// returns the value of the variable @a do_dynamic_scaling, see there
  bool get_use_variable_metric() const {return Hp->employ_variable_metric();}
    
  /// returns the number of calls to the oracle since the last clear()
  CH_Matrix_Classes::Integer get_cntobjeval() const {return cntobjeval;}
  /// returns the number of descent steps since the last clear()
  CH_Matrix_Classes::Integer get_descent_steps() const {return descent_steps;}
  /// returns the number of bundle method iterations since the last descent step
  CH_Matrix_Classes::Integer get_innerit() const {return innerit;}
  /// returns the number of bundle method itrations since the last clear()
  CH_Matrix_Classes::Integer get_suminnerit() const {return suminnerit;}
  /// returns the number of model qp subproblems since the last clear()
  CH_Matrix_Classes::Integer get_sumupdatecnt() const {return sumupdatecnt;}
  /// returns the number of oracle reevaluations for the center due to numerical problems since the last descent step 
  CH_Matrix_Classes::Integer get_recomp() const {return recomp;}
  /// returns the number of oracle reevaluations for the center due to numerical problems since the last clear() or clear_fails() 
  CH_Matrix_Classes::Integer get_sumrecomp() const {return sumrecomp;} 
  /// returns the number of fails in qp subproblems since the last null/descent step
  CH_Matrix_Classes::Integer get_qpfails() const {return qpfails;}
  /// returns the number of fails in qp subproblems since the last clear() or clear_fails()
  CH_Matrix_Classes::Integer get_sumqpfails() const {return sumqpfails;}
  /// returns the number of fails in model evaluatoins since the last null/descent step
  CH_Matrix_Classes::Integer get_modelfails() const {return modelfails;}
  /// returns the number of fails in model evaluations since the last clear() or clear_fails()
  CH_Matrix_Classes::Integer get_summodelfails() const {return summodelfails;}
  /// returns the number of failures to increase the augmented model value since the last null/descent step
  CH_Matrix_Classes::Integer get_augvalfails() const {return augvalfails;}
  /// returns the number of failures to increase the augmented model value since the last clear() or clear_fails()
  CH_Matrix_Classes::Integer get_sumaugvalfails() const {return sumaugvalfails;}
  /// returns the number of fails in oracle evaluations since the last descent step
  CH_Matrix_Classes::Integer get_oraclefails() const {return oraclefails;}
  /// returns the number of fails in oracle evaluations since the last clear() or clear_fails()
  CH_Matrix_Classes::Integer get_sumoraclefails() const {return sumoraclefails;}

  /// returns the number of oracle evaluations that returned an epsilon subgradient that improved the model by a dangerously small amount (mostly this is due to solving the QP subproblems only approximately)
  CH_Matrix_Classes::Integer get_shallowcut() const {return shallowcut;}
  
  /// return a pointer to the groundset 
  const Groundset* get_groundset() const {return groundset;}
  /// return a pointer to the cutting model
  const BundleModel* get_model() const {return model;}
  /// return a pointer to the quadratic term of the proximal term
  BundleProxObject* get_prox() const {return Hp;}
  /// return a pointer to the termination criterion
  BundleTerminator* get_terminator() const {return terminator;}
  /// return a pointer to the class for updating the weightu of the proximal term 
  BundleWeight*  get_bundleweight() const {return bundleweight;}

  /// return time spent in computing the cost coefficients of the quadratic bundle subproblem
  CH_Tools::Microseconds get_QPcoeff_time() const {return QPcoeff_time;} 
  /// return time spent in solving the quadratic bundle subproblem
  CH_Tools::Microseconds get_QPsolve_time() const {return QPsolve_time;}
  /// return time spent in forming the model aggregate
  CH_Tools::Microseconds get_make_aggr_time() const {return make_aggr_time;}
  /// return time spent in total for the quadratic bundle subproblem
  CH_Tools::Microseconds get_evalaugmodel_time() const {return evalaugmodel_time;}
  /// output some time statistic paramters
  std::ostream& print_statistics(std::ostream& out) const;
  /// output the data of the Gauss-Seidel qp to in an m file format
  int qp_mfile_data(const CH_Matrix_Classes::Matrix& center_y,
		    const BundleProxObject* Hp,
		    const MinorantPointer& gs_subg,
		    const CH_Matrix_Classes::Symmatrix& Q,
		    const CH_Matrix_Classes::Matrix& c,
		    CH_Matrix_Classes::Real offset,
		    const CH_Matrix_Classes::Indexmatrix& yfixed
		    ) const;
  //@}


};

  //@}

}

#endif

