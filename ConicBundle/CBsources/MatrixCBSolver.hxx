/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/MatrixCBSolver.hxx
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



#ifndef CONICBUNDLE_MATRIXCBSOLVER_HXX
#define CONICBUNDLE_MATRIXCBSOLVER_HXX

/**  @file MatrixCBSolver.hxx
    @brief Header declaring the classes ConicBundle::MatrixCBSolver, ConicBundle::MatrixFunctionOracle, ConicBundle::PrimalMatrix, ConicBundle::MatrixMinorant, ConicBundle::ModifiableOracleObject
    @version 1.0
    @date 2005-03-01
    @author Christoph Helmberg
*/


#include "CBSolver.hxx"
#include "matrix.hxx"
#include "UQPSolver.hxx"
#include "BundleProxObject.hxx"
#include "BundleData.hxx"
#include "BundleWeight.hxx"
#include "FunctionObjectModification.hxx"
#include "AffineFunctionTransformation.hxx"
#include "SumBundleParametersObject.hxx"

//------------------------------------------------------------


namespace ConicBundle {

  /**@defgroup cxxmatrixinterface Interface to ConicBundle for the Language C++ using Matrix Classes
     @brief Solve \f$min_{y\in\mathbf{R}^m}  f_0(y) + f_1(y) + ... + f_k(y)\f$
      for convex functions f_i, the y-variables may be free, bounded,
      box constrained, or linearly constrained. The most important steps are the following.

     Before starting, note that this interface relies on several classes
     defined in CH_Matrix_Classes, in particular on
     CH_Matrix_Classes::Matrix and CH_Matrix_Classes::Indexmatrix. The
     full functionality of ConicBundle is only available with these
     classes. If you prefer an interface without having to use these,
     please use the
     \ref cxxinterface "Interface to ConicBundle for the Language C++".
     We now give a short overview of the most important
     steps in using the ConicBundle solver.

     <b>Setting up the Problem, the Functions, and the Main Loop</b>

     First create a new problem/solver ConicBundle::MatrixCBSolver, let
     us call it solver for brevity.

     Next call ConicBundle::MatrixCBSolver::init_problem() to set the dimension
     of the design variables/argument space, possibly together with box constraints,
     starting values and and linear cost terms. In prinicple,
     further linear constraints could be installed using
     ConicBundle::MatrixCBSolver::append_constraints() but this should
     be avoided whenever possible (additional constraints
     may entail a significant loss in efficiency in the internal quadratic
     subproblems. Also note that even simple bounds require a lot more work when
     using general variable metric approaches.

     Now set up each of your functions f_i as a
     ConicBundle::MatrixFunctionOracle (see \ref cbsolver for specialized oracle classes).  The solver will call the routine
     ConicBundle::MatrixFunctionOracle::evaluate(), and via this you have to
     supply, for a given argument, a function value and an affine
     ConicBundle::Minorant (you might want to use the constructors provided in
     MatrixMinorant instead). A Minorant is described by its function value and
     (epsilon) subgradient (=the gradient if the function is differentiable).
     In the oracle this function evaluate() is the only function that you
     definitely have to provide, see the miniature example below.

     The function oracles have to be added to the solver
     using the routine  ConicBundle::MatrixCBSolver::add_function().
     In adding a function \f$f\f$ you may (but need not) specify a multiplicative
     factor \f$\gamma\f$ (==1. by default) for \f$f\f$ and whether the function
     serves as an ObjectiveFunction (the solver uses the value \f$\gamma f(y)\f$ directly),
     as a ConstantPenaltyFunction (the solver uses the value \f$\gamma\max\{0,f(y)\}\f$)
     or as an AdaptivePenatlyFunction (the solver uses the value \f$\gamma\max\{0,f(y)\}\f$
     but may further increase \f$\gamma\f$ until the value indeed goes to zero).
     If your function works just on a subset/subspace of the variables, it may be worth
     to specify the function on its own original space and use a
     ConicBundle::AffineFunctionTransformation on the argument, so that
     the solver uses, e.g. \f$ \gamma\cdot(\delta+b^Ty+\rho\cdot f(c+Ay))\f$.

     Once all functions are added, the optimization process can be
     started. If you know a good starting point (and have not already
     supplied it  in init_problem()) then set it with
     ConicBundle::MatrixCBSolver::set_new_center_point() now, otherwise
     the method will pick the zero vector or, in the case of box
     constraints, the point closest to zero as starting point.

     Finally, call ConicBundle::MatrixCBSolver::solve() and retrieve
     ConicBundle::MatrixCBSolver::termination_code() for getting the
     reason for termination. Via parameters in solve() you may also
     tell the solver to return after each descent step or also
     after a maximum number of null steps. This then only interrupts
     computations and calling solve() again continues as if there
     was not break at all.

     <b>Retrieving Some Solution Information</b>

     Whenever the solver returns form solve() you can retrieve the
     current objective value by
     ConicBundle::MatrixCBSolver::get_objval() and the argument leading
     to this value by ConicBundle::MatrixCBSolver::get_center(). For
     some screen output, use ConicBundle::MatrixCBSolver::set_out().

     <b>Lagrangean Relaxation, Primal Approximations, and Cutting Planes</b>

     If you are optimizing the Lagrange multipliers of a Lagrangean relaxation,
     you might be interested in getting an approximation to your primal optimal
     solution. This can be done by specifying in each function for each
     minorant/(epsilon) subgradient the corresponding primal vectors that
     generate it, see PrimalData (with implementation PrimalMatrix) in
     ConicBundle::Minorant as a start. Then for each of your functions, you can
     retrieve the current primal approximation using
     ConicBundle::MatrixCBSolver::get_approximate_primal().

     If you are just relaxing linear constraints on your variables, an elegant
     alternative is to use an AffineFunctionTransformation as follows.  For
     concreteness, assume you want to minimize
     \f$f(y)=\max_{x\in\mathcal{X}}[h(x)+(b-Bx)^Ty]\f$ where
     \f$\mathcal{X}\subset\mathbf{R}^n\f$ is some compact set for which the
     inner max is easy to compute for any \f$y\f$. Now, it might actually be
     better to implement the function \f$\bar f(\bar
     c)=\max_{x\in\mathcal{X}}[h(x)+\bar c^Tx]\f$ and use it with the
     AffineFunctionTransformation \f$ b^Ty+ \bar f(-B^Ty)\f$ (so
     \f$A=-B^T\f$). Then if \f$\bar x\f$ is the maximizer for a given \f$\bar
     c=-B^Ty\f$, the corresponding minorant to be returned by the oracle of
     \f$\bar f\f$ is described by the offset \f$h(\bar x)\f$ with the
     corresponding subgradient being just \f$\bar x\f$ itself. So in this case
     the primal data is the minorant itself and there is no need to provide it
     separately. For the case that \f$\mathcal{X}\f$ is a box domain, see
     \ref box_oracle for a tuned model. A slightly different possibility
     is offered by the example \ref implemented_matrix_function_oracle
     with an NNCBoxSupportFunction,
     but this is likely a lot less efficient.

     In addition, you might plan to improve your primal relaxation via
     cutting planes if they are strongly violated by the current primal
     approximation. This requires introducing new multipliers for the
     ground set and typically also changes the dimension of your
     function \f$f\f$ on the fly, unless you are in the situation of
     \f$\bar f\f$ above, where only the AffineFunctionTransformation
     changes. In any case you will need
     ConicBundle::MatrixCBSolver::append_variables().  There you have to
     provide for your function either a corresponding derived
     implementation of the base class ConicBundle::OracleModification
     or, in the case of \f$\bar f\f$ above, just a corresponding
     ConicBundle::AFTModification for your transformation.  In the
     latter case this is all you need to provide, because \f$\bar f\f$
     is not changed at all, but in the general
     case of a general \f$f\f$ you will also need to implement the
     oracle's subroutine MatrixFunctionOracle::apply_modification()
     which will be called by the solver with your modification data so
     that the function adapts according to this data. In most cases
     adapting the class ConicBundle::Modification like in
     ConicBundle::LPGroundsetModification will allow to do all you need
     with minimal work.

     Typically, when changing the dimension, the cutting model of the
     bundle method is lost. This is no problem in the case of \f$\bar f\f$
     above if the newly added multipliers are zero initially, because
     then the minorants with transformation can and will be regenerated
     automatically from the unchanged minorants \f$\bar x\f$. In the case
     of general \f$f\f$, if you provided primal data with your
     minorants, you may also know how to extend your minorants on
     newly appended coordinates. If you know how to do this, then you
     can provide a ConicBundle::MinorantExtender when the solver calls
     MatrixFunctionOracle::apply_modification(). The MinorantExtender
     will be applied to all minorants of the model and if this succeeds
     the model is not lost. A rich example implementation
     (with a separate cutting model) of all these features is provided
     by the PSCAffineFunction implementation of
     the PSCOracle for solving the duals to primal
     semidefinite programs with nonempty bounded feasible sets,
     see \ref abstract_psc_oracle and \ref implemented_psc_oracle.

     If you want to get rid of primal constraints/dual variables, have a
     look at ConicBundle::MatrixCBSolver::set_active_bounds_fixing(),
     ConicBundle::MatrixCBSolver::get_fixed_active_bounds(), maybe also
     at ConicBundle::MatrixCBSolver::get_approximate_slacks() and then
     use ConicBundle::MatrixCBSolver::delete_variables(). This may again
     need to be accompanied by OracleModifcation implementations and
     corresponding actions in your oracle.

     @include mat_mini_ex.cxx

  */
  //@{

  /** @brief If in Lagrangean relaxation primal solutions are in the form of a real vector
or, more generally a matrix, then an approximate primal solution can be generated by supplying primal information of this form for each epsilon subgradient within ConicBundle::MatrixFunctionOracle::evaluate().

      In many applications, e.g. in Lagrangean relaxation, the convex
      minimization problem arises as the dual of a convex primal
      maximization problem. In this case one is typically interested
      in obtaining a primal approximate solution in addition to the
      dual solution. Under reasonable conditions this is possible if the
      primal solutions that give rise to the subgradients are
      aggregated along with the subgradients within the bundle algorithm.
      If the primal data can be represented as a CH_Matrix_Classes::Matrix
      then the user has to supply in the oracle for each sugradient the
      corresponding primal data in a PrimalMatrix and the algorithm will do the rest.
      Observe that a PrimalMatrix can be used exactly in the same way as
      a CH_Matrix_Classes::Matrix and that they are assignable among each other.

      The primal data has to be supplied within ConicBundle::MatrixFunctionOracle::Evaluate()
      and can be retrieved via the methods
      ConicBundle::MatrixCBSolver::get_approximate_primal() and
      ConicBundle::MatrixCBSolver::get_center_primal()
 */

  class PrimalMatrix :public PrimalData, public CH_Matrix_Classes::Matrix {

  public:
    /// empty matrix
    PrimalMatrix() {
    }
    /** @brief generate a matrix of size nr x nc but WITHOUT initializing the memory

        If initializing the memory externally and CONICBUNDLE_DEBUG is defined, please use
        set_init() via matrix.set_init(true) in order to avoid warnings concerning improper
        initialization
    */
    PrimalMatrix(CH_Matrix_Classes::Integer nr, CH_Matrix_Classes::Integer nc) :Matrix(nr, nc) {
    }
    /// generate a matrix of size nr x nc initializing all elements to the value d
    PrimalMatrix(CH_Matrix_Classes::Integer r, CH_Matrix_Classes::Integer c, CH_Matrix_Classes::Real d) :Matrix(r, c, d) {
    }
    /// copy constructor, *this=pm
    PrimalMatrix(const PrimalMatrix& pm) :PrimalData(), Matrix(pm) {
    }
    /// copy constructor, *this=pm
    PrimalMatrix(const CH_Matrix_Classes::Matrix& pm) :Matrix(pm) {
    }
    /// copy operator, *this=pm
    PrimalMatrix& operator=(const CH_Matrix_Classes::Matrix& pd) {
      Matrix::operator=(pd); return *this;
    }

    /// produces a new PrimalMatrix that is a copy of itself; the caller has to delete the returned object at some point
    PrimalData* clone_primal_data() const {
      return new PrimalMatrix(*this);
    }

    /// multiply *this Matrix with myfactor and add itsfactor*it (it must dynamic_cast to a PrimalMatrix)
    int aggregate_primal_data(const PrimalData& it, double itsfactor) {
      const PrimalMatrix* pd = dynamic_cast<const PrimalMatrix*>(&it);
      if (pd != 0) {
        xpeya(*pd, itsfactor);
        return 0;
      }
      return 1;
    }

    /// multiply/scale *this with a nonnegative myfactor
    int scale_primal_data(double myfactor) {
      *this *= myfactor; return 0;
    }

  };

  /**@brief Minorant interface supporting Matrix classes; simple constructros for subgradients given by column vectors in form of a CH_Matrix_Classes::Matrix or a CH_Matrix_Classes::Sparsemat

      Suppose in the evaluation of your oracle at the current point \f$y\f$ you determine the subgradient \f$s\f$ for the subgradient inequality

      \f[ f(z)\ge f(y)+\langle s,z-y\rangle\quad\forall z\in\mathbf{R}^m, \f]

      then use \f$f(y)\f$ for @a offset and \f$s\f$ for @a subg in the
      constructors and keep the default value @a offset_at_origin == false.

      If, on the other hand, your oracle implements a support function for some compact set \f$\mathcal{X}\subset\mathbf{R}^m\f$ like

      \f[ f(y) = \max_{x\in\mathcal{X}} x^\top y\f]

      then it is actually more efficient to return 0 for @a offset,
      a maximizing \f$x\f$ in @ subg and to put @a offset_at_origin = true.

   */

  class MatrixMinorant : public Minorant {
  public:

    /// initializaton by a Matrix column vector; if offset_at_origin==true, the offset gives the value at the origin (in contrast to the value at the current evaluation point)
    MatrixMinorant(CH_Matrix_Classes::Real offset,
      const CH_Matrix_Classes::Matrix& subg,
      PrimalData* primal = 0,
      bool offset_at_origin = false) :
      Minorant(offset_at_origin, double(offset), int(subg.rowdim()), (const double*)subg.get_store(), 0, 1., primal) {
      assert(subg.coldim() == 1);
    }

    /// initializaton by a Sparsemat column vector; if offset_at_origin==true, the offset gives the value at the origin (in contrast to the value at the current evaluation point)
    MatrixMinorant(CH_Matrix_Classes::Real offset,
      const CH_Matrix_Classes::Sparsemat& subg,
      PrimalData* primal = 0,
      bool offset_at_origin = false) :
      Minorant(offset_at_origin, double(offset), int(subg.nonzeros()), (const double*)(subg.get_colval().get_store()), (const int*)(subg.get_colindex().get_store()), 1., primal) {
      assert(subg.coldim() == 1);
    }

    ///generates a MatrixMinorat copy of a Minorant scaled by factor (deault =1.) and optionally without cloning the primal (if it exists)
    MatrixMinorant(const Minorant* mnrt,
      CH_Matrix_Classes::Real factor = 1.,
      bool with_primal = true)
      :Minorant(mnrt, double(factor), with_primal) {
    }

    ///
    ~MatrixMinorant() {
    }

    /// produces a clone of itself
    Minorant* clone_minorant(double factor = 1., bool with_primal = true) const {
      return new MatrixMinorant(this, CH_Matrix_Classes::Real(factor), with_primal);
    }

  };

  /**@brief ModifiableOracle provides all oracles with a uniform interface for a modification routine and an on/off switch for internal correctness checks
   */
  class ModifiableOracleObject : public FunctionObject {
  public:
    ///virtual destructor
    virtual ~ModifiableOracleObject();

    /**@brief This routine need not be implemented unless variables
      (constraints in Lagrangean relaxation) are added or deleted on
      the fly

      The routine is only called by the solver if the variables indeed get
      modified. @a oracle_modification is then used to either transfer
      user supplied instructions to the oracle on how to modify itself
      or to inform the oracle about changes the solver was asked to perform
      on the variables. If available, the solver will also show the effect
      of these changes on the center point in @a new_center and @a old_center;
      if these are not available then they hold NULL. A user supplied
      @a oracle_modification will be checked for consistency with the
      actual changes in the variables and mismatches will cause failures.

      The remaining variables are output variables by which the oracle
      tells the solver which information has a chance to be preserved
      in view of these changes.  If e.g. the deletion of some nonzero
      variables invalidates the function value in the new center, the
      oracle has to set discard_objective_in_center=true.  If the
      entire model cannot be preserved (this includes the aggregates
      and the function values), the oracle needs to set
      discard_model=true; If only aggregate minorants cannot be preserved,
      the oracle needs to set discard_aggregates=true. Whenever new
      variables were added, the model can only be preserved if the
      remaining minorants (maybe without aggregates) can be extended
      for these new variables. In this case the oracle has to supply
      the appropriate MinorantExtender via @a minorant_extender and
      only those minorants will be kept for which this operation succeeds.

      Return value 0 indicates that these actions allow to continue without
      errors, other return values result in an overall error on these changes.
      The default implementation just returns 1 and does not allow
      modifications.
    */

    virtual
      int
      apply_modification
      (
        const OracleModification& oracle_modification,
        const CH_Matrix_Classes::Matrix* new_center,
        const CH_Matrix_Classes::Matrix* old_center,
        bool& discard_objective_in_center,
        bool& discard_model,
        bool& discard_aggregates,
        MinorantExtender*& minorant_extender
      );

    /**@brief switch on/off some correctnes checks on the oracle */
    virtual
      bool
      check_correctness() const {
      return false;
    }

  };


  /**@brief Oracle interface (abstract class). For each of your functions, provide an instance of a derived class.

     The oracle interface is used to describe and pass convex
     objective functions to the ConicBundle::MatrixCBSolver.  The
     dimension of the argument vector of the function (or, if an
     ConicBundle::AffineFunctionTransformation is supplied for this
     function, the dimension of the transformation argument) must be
     set in ConicBundle::MatrixCBSolver::init_problem() and the
     functions are then added to the solver by
     ConicBundle::MatrixCBSolver::add_function().

     If the sum of several such functions is to be minimized it is the
     task of the user to guarantee, that all dimensions match, e.g. by
     using an appropriate AffineFunctionTransformation for each.

     If the function corresponds to Lagrangean relaxation of a primal
     maximization problem one may want to generate a primal
     approximate solution. In this case, return in the function
     evaluate() the generating primal objects within each minorant.
     If no primal objects are included, there will be no primal
     aggregation.

     If primal aggregation is used then it is possible to implement
     a primal cutting plane framework. This requires the introduction
     of new (dual) variables in the design space of the function.
     In this case the function MinorantExtender returned in
     apply_modification() serves the
     purpose of filling in the missing coordinates in existing
     subgradients. If this feature is not needed, the function
     may be used as is and need not be reimplemented.
  */

  class MatrixFunctionOracle : public ModifiableOracleObject {
  public:

    /** @brief Called by the solver. Has to Return function value and
               at least one (epsilon) subgradient and, possibly for Lagrangean
               relaxation, some primal data.

        The evaluation method is the main interface to the bundle
        solver.  The solver calls this method to obtain for the @a
        current_point (its dimension is set in
        ConicBundle::CBSolver::init_problem() or via a
        AffineFunctionTransformation specified in
        MatrixCBSolver::add_function()) the @a objective_value and
        (epsilon) subgradient information.  In any call several
        epsilon subgradients may be returned in @a minorants along
        with their offset, but at least one has to be returend. Each
        subgradient and offset describes a linear minorant to the
        convex function and is used by the bundle method to form a
        cutting model of the convex function.

        In many applications, computing the function value is an iterative
        process that approaches the true function value from below.
        On input the code offers a bound in @a objective_value for the
        function value, above which it is certain that the code will reject
        the current point. If in the iterative
        process a lower bound on the function value exceeds this bound,
        then it is sufficient to return, instead of the true function value
        and a subgradient, the current lower bound and a vector so that
        together they describe a supporting hyperplane (i.e. a linear minorant)
        to the function at this point.

        If the function corresponds to Lagrangean relaxation of a
        primal maximization problem one may want to generate a primal
        approximate solution. For this purpose the minorants may also
  hold information on the primal data. If at each call and for each
  epsilon subgradient the corresponding generating primal object
        (must be derived from ConicBundle::PrimalData, e.g., a
        ConicBundle::PrimalDVector) is stored, then the code automatically
        performs the aggregation corresponding to the aggregation
        of the subgradients on the primal data objects. The
        primal approximate solution is finally delivered by
        the methods ConicBundle::CBSolver::get_approximate_primal()
        or  ConicBundle::CBSolver::get_center_primal().

  All minorants passed to the solver must be objects allocated
        on the heap. The ownership of these objects is transferred
        to the solver and the solver will destroy them eventually,
        so DO NOT delete them yourself!

        If no primal aggregation is desired, simply do not touch
        @a primal_data or clear it.

        @param[in] current_point (const Matrix&)
           argument of the function as a column vector (e.g. the Lagrange multipliers)

        @param[in] relprec (double)
           relative precision requirement for objective values
           that may lead to descent steps (this precision is not
           required if it is already certain that the function
           value will be too poor)

        @param[in,out]  objective_value (double&)
         - on input:
           value gives the threshold for a null step; you may stop,
           if a cutting plane yields at least this;
         - on output:
           return an upper bound on the true function value within
           @a relprec *(abs(objval)+1.), if there is no linear minorant
           cutting above the threshold specified in objective_value on input.
           Otherwise the return value should be the max of cut_values.

        @param[out]  minorants (std::vector<Minorant*>)
           for returning linear minorants (subgradients and their offset
     values). At least one must be returned and the one at index 0
     must maximize the value at the current point over all further
     ones given here. In particular its value at point should be
     above the threshold or be within relprec *(abs(objval)+1.)
           of the true objective.  Within the minorants primal data
           may be supplied if this should be aggregated along.

         @param[out]  primal_extender (PrimalExtender*&)
           if primal_data of previous calls has now to be updated
           due to changes in the primal problem -- e.g., this may happen
           in column generation -- one may return a pointer to
           PrimalExtender object on the heap. This object will be used
           by ConicBundle to update all its internally stored
           primal_data objects by calling PrimalExtender::extend on
           each of these, afterwards ConicBundle deletes primal_extender
           If this is not needed, the variable holds 0.

        @return
           -  0,  all correct
           -  !=0, failure. This does not necessarily terminate the bundle method.
               Termination is forced only if no new subgradient is returned.
     */
    virtual
      int
      evaluate
      (
        const  CH_Matrix_Classes::Matrix& current_point,
        CH_Matrix_Classes::Real relprec,
        CH_Matrix_Classes::Real& objective_value,
        std::vector<Minorant*>& minorants,
        PrimalExtender*& primal_extender
      )
      = 0;


    /**@brief switch on/off some correctnes checks on the oracle */
    virtual
      bool
      check_correctness() const {
      return true;
    }
  };


  class MatrixCBSolverData;

  /**@brief  The Full Conic Bundle method solver invoked by ConicBundle::MatrixCBSolver(), it uses a separate cutting model for each function

  Minimizes the sum of convex functions that are given via
  ConicBundle::MatrixFunctionOracle interfaces, see
  \ref cxxmatrixinterface "the text explaining the C++ interface for Matrix Classes" for
  a quick overview.

  It provides special support for Lagrangean relaxation by generating
  primal approximate solutions if such information is provided in the
  function oracles.

  Based on these primal approximations it is also possible to implement
  cutting plane schemes. Routines for adding and deleting corresponding
  dual variables as well as a framework for extending subgradients in order
  not to loose the cutting model are available.

  */
  class MatrixCBSolver : public CBout {
  private:
    MatrixCBSolverData* data_;                       ///< pointer to internal solver data
    MatrixCBSolver(const MatrixCBSolver&);           ///< not available, blocked deliberately
    MatrixCBSolver& operator= (const MatrixCBSolver&);///< not available, blocked deliberately
  public:

    /// default constructor allows to set output level options from start (see also set_out())
    MatrixCBSolver(std::ostream* out = 0, int print_level = 0);
    /// constructor for setting output options depending on existing CBout objects 
    MatrixCBSolver(const CBout* cb, int incr = 0);
    ///
    ~MatrixCBSolver();

    //------------------------------------------------------------
    /**@name Initialization */
    //@{

    /** @brief Clears all data structures and problem information
        but keeps ouptut settings and algorithmic parameter settings
    */
    void clear();

    /** @brief Sets default values for algorithmic parameters that are not function specific (e.g., relative precision, weight and weight bounds for the augmentedproblem, etc.)
    */
    void set_defaults();

    /** @brief Initializes the problem by setting up the design space
        (the dimension and possibly box constraints on the variables)

    Clears all data structures and sets the dimension @ m for a new problem.
    for solving   min_{y in R^m}  f_0(y) + f_1(y) + ...
    Box constraints may be specified for y. (The functions f_i must be added
     by add_function()).

     Lower and/or upper bounds must be speicified for all variables
     or for none of them. To specify no bounds at all, give Null
     pointers. Otherwise use ConicBundle::CB_minus_infinity for
     unbounded below and ConicBundle::CB_plus_infinity for unbounded above.
     For NULL pointers, unbounded will be used as default for all
     variables. Specifying bounds selectively is also possible
     by set_lower_bound() or set_upper_bound(). For further constraints
     see append_constraints().

     @param[in] dim  (int)
         the dimension of the argument/design space/the number of Lagrange multipliers

     @param[in] lbounds  (const Matrix*)
         If NULL, all variables are considered unbounded below,
         otherwise lbounds[i] gives the minimum feasible value for variable y[i],
         use ConicBundle::CB_minus_infinity for unbounded below.

     @param[in] ubounds (const Matrix*)
         If NULL, all variables are considered unbounded above,
         otherwise ubounds[i] gives the maximum feasible value for variable y[i],
         use ConicBundle::CB_plus_infinity for unbounded above.

     @param[in] startval (const Matrix*)
        If NULL, the starting values are obtained by projecting
        zero onto the feasible set given by the lower and upper bounds
  resulting from the arguments before

     @param[in] costs (const Matrix*)
         Use this in order to specify linear costs on the variables in addition
   to the functions (may be convenient in Lagrangean relaxation for
   the right hand side of coupling contsraints); NULL is equivalent
   to costs zero.

     @param[in] offset (Real)
         Use this in order to specify linear costs on the variables in addition
   to the functions (may be convenient in Lagrangean relaxation for
   the right hand side of coupling contsraints); NULL is equivalent
   to costs zero.

     @return
        - 0 on success
        - != 0 otherwise


     */
    int init_problem(int dim,
      const CH_Matrix_Classes::Matrix* lbounds = 0,
      const CH_Matrix_Classes::Matrix* ubounds = 0,
      const CH_Matrix_Classes::Matrix* startval = 0,
      const CH_Matrix_Classes::Matrix* costs = 0,
      CH_Matrix_Classes::Real offset = 0.);



    /** @brief Adds a function, typically derived from ConicBundle::FunctionOracle. If the dimension does not match the current one, specify an affine function transformation to map the current ground set to the argument of the function.

     Besides the standard ConicBundle::MatrixFunctionOracle
     the interface only accepts a few other prespecified derivations
     of the class FunctionObject that come along with the CH_Matrix_Classes interface
     (e.g. for semidefinite and second order cones). Functions not derived from
     these will fail to be added and return a value !=0.

     The @a fun_factor allows to specify a scaling factor for the function. @a fun_factor must be a strictly positive number.

     The ConicBundle::FunctionTask @a fun_task specifies whether the
     function is to be used as a an ObjectiveFunction, a
     ConstantPenaltyFunction with @a fun_factor as maximum penalty
     factor, or as an AdaptivePenaltyFunction with @a fun_factor at
     initial penalty guess that might be increased or decreased over
     time.

     The AffineFunctionTransformation @a aft may be used to modify the
     argument and give an additional affine term (linear term plus
     offset). For adding an affine term there are several other
     possibilities, e.g. in init_problem(), so there is no need to do so
     here. If, however, an existing function implementation requires only
     some subset of the variables, it is more convenient to supply
     a corresponding @a aft instead of reimplementing the function.

     @a argument_list_may_change_dynamically sets a flag on how to
     treat the function arguments when variables are added or deleted.
     If the arguments may not change, any changes in the variables are
     mapped to an adaptation of an internal
     AffineFunctionTransformation so that the function does not notice
     the changes in the variables.  If arguments may change, the
     function oracle should be a ModifiableOracleObject and react
     accordingly to the changes in its
     ModifiableOracleObject::apply_modification() routine.


    @return
      - 0 on success
      - != 0 otherwise
    */

    int
      add_function
      (FunctionObject& function,
        CH_Matrix_Classes::Real fun_factor = 1.,
        FunctionTask fun_task = ObjectiveFunction,
        AffineFunctionTransformation* aft = 0,
        bool argument_list_may_change_dynamically = false);

    /**@brief Sets lower bound for variable i,
  use ConicBundle::CB_minus_infinity for unbounded from below.

       The algorithm may have to adapt the center point aftwards.
       In this case the old function values will be marked as outdated and
       will be recomputed at the next call to e.g. solve().

    @return
      - 0 on success
      - != 0 otherwise
    */
    int set_lower_bound(int i, double lb);

    /**@brief Sets upper bound for variable i,
  use ConicBundle::CB_plus_infinity for unbounded from below.

       The algorithm may have to adapt the center point aftwards.
       In this case the old function values will be marked as outdated and
       will be recomputed at the next call to e.g. solve().

     @return
      - 0 on success
      - != 0 otherwise
    */
    int set_upper_bound(int i, double ub);


    /** @brief Append new variables (always in last postions in this order).

     @attention Be sure to include a desription of required changes to your
       functions via @a affected_functions_with_modifications

     @param[in] n_append  (int)
       number of variables to append (always in last position in the same order)

     @param[in] lbounds  (const Matrix*)
        If NULL, all appended variables are considered unbounded below,
        otherwise lbounds[i] gives the minimum feasible value for variable y[i],
        use ConicBundle::CB_minus_infinity for unbounded below.

     @param[in] ubounds (const Matrix*)
        If NULL, all appended variables are considered unbounded above,
        otherwise ubounds[i] gives the maximum feasible value for variable y[i],
        use ConicBundle::CB_plus_infinity for unbounded above.

     @param[in] constraint_columns (const Sparsemat*)
        This must be NULL unless append_constraints() has been used before for
  specifying linear constraints on the ground set; if there are constraints,
  NULL is interpreted as appending zero columns to the constraints,
  otherwise the the number of rows of the Sparsemant has to match the
  current number of linear constraints and the number of columns the
  must equal n_append.

     @param[in] startval (const Matrix*)
        If NULL, the starting values are obtained by projecting
        zero onto the feasible set given by the lower and upper bounds
  resulting from the arguments before

     @param[in] costs (const Matrix*)
         Use this in order to specify linear costs on the variables in addition
   to the functions (may be convenient in Lagrangean relaxation for
   the right hand side of coupling contsraints); NULL is equivalent
   to costs zero.

     @param[in,out] affected_functions_with_modifications (const FunObjModMap*)
        If NULL, default actions are performed on all functions. In
        particular, those admitting dynamic argument changes will get
        the new variables appended at the end of their argument vector
        and (eventually) their apply_modification() routines will be
        called informing them about the groundset changes; for those
        not admitting changes in their arguments, their corresponding
        (possibly newly created) affine function transformation will
        be set up to ignore the new arguments.  If !=NULL, for the
        listed functions (and their parents up to the root function)
        the default appending action is performed unless their
        FunctionObjectModification entry gives explicit modification
        instructions which are then applied instead. For all functions
        NOT listed in the map and not having modified offsprings their
        corresponding aft will be set up to ignore the new variables.

     @return
        - 0 on success
        - != 0 otherwise


    */
    int append_variables(int n_append,
      const CH_Matrix_Classes::Matrix* lbounds = 0,
      const CH_Matrix_Classes::Matrix* ubounds = 0,
      const CH_Matrix_Classes::Sparsemat* constraint_columns = 0,
      const CH_Matrix_Classes::Matrix* startval = 0,
      const CH_Matrix_Classes::Matrix* costs = 0,
      const FunObjModMap* affected_functions_with_modifications = 0);

    /** @brief Deletes variables corresponding to the specified indices.

        The indices of the remaining variables are reassigned so that they
        are consecutive again, the routine returns in @a map_to_old
        a vector giving for each new index of these remaining variables
        the old coordinate.

     @attention Be sure to include a desription of required changes to your
       functions via @a affected_functions_with_modifications

     @param[in] delete_indices  (const Indexmatrix&)
        the entries delete_indices[i] specify the indices of the variables
        to be deleted

     @param[out] map_to_old  (Indexmatrix&)
        after the call, element map_to_old[i] gives the old index (before the call)
        of the variable that now has index position i.

     @param[in] affected_functions_with_modifications (const FunObjModMap*)
        If NULL, default actions are performed on all functions. In
        particular, for those admitting dynamic argument changes all
  those variables will be deleted whose row in a corresponding
  updated affine function transformation (so after deletion of
  the columns of the incoming variables) corresponds to the zero
  map (i.e., offset and matrix row are both zero), furthermore
  identity transformations will be preserved. For those
        not admitting changes in their arguments, their corresponding
        (possibly newly created) affine function transformation will
        only get the columns deleted, but there will be no row deleltions.
  If !=NULL, for the listed functions (and their parents up to the
  root function) the default deletion action is performed unless their
        FunctionObjectModification entry gives explicit modification
        instructions which are then applied instead. For all functions
        NOT listed in the map and not having modified offsprings their
        corresponding aft will be set up to keep the arguments unchanged.

     @return
        - 0 on success
        - != 0 otherwise

    */
    int delete_variables(const CH_Matrix_Classes::Indexmatrix& delete_indices,
      CH_Matrix_Classes::Indexmatrix& map_to_old,
      const FunObjModMap* affected_functions_with_modifications = 0);

    /** @brief Reassigns variables to new index positions by mapping to position @a i
        the variable that previously had index @a assign_new_from_old[i].

        Old variables, that are not mapped to any position will be deleted.
        It is not allowed to generate several copies of old variables.

     @attention Be sure to include a desription of required changes to your
       functions via @a affected_functions_with_modifications

     @param[in] assign_new_from_old  (const IVector&)
        entry assign_new_from_old[i] specifies
        the old index of the variable, that has to be copied to index position i.

     @param[in] affected_functions_with_modifications (const FunObjModMap*)
        If NULL, default actions are performed on all functions. In
        particular, for those admitting dynamic argument changes all
  those variables will be deleted whose row in a corresponding
  updated affine function transformation (so after mapping
  the columns of the incoming variables) correspond to the zero
  map (i.e., offset and matrix row are both zero); furthermore,
  if the transformation was the identity to start with, this will
  be preserved by mapping the arguments in the same way. For those
        not admitting changes in their arguments, their corresponding
        (possibly newly created) affine function transformation will
        only get the columns mapped, but there will be no row deleltions.
  If !=NULL, for the listed functions (and their parents up to the
  root function) the default deletion action is performed unless their
        FunctionObjectModification entry gives explicit modification
        instructions which are then applied instead. For all functions
        NOT listed in the map and not having modified offsprings their
        corresponding aft will be set up to keep the arguments unchanged.

     @return
        - 0 on success
        - != 0 otherwise

    */
    int reassign_variables(const CH_Matrix_Classes::Indexmatrix& assign_new_from_old,
      const FunObjModMap* affected_functions_with_modifications = 0);

    /** @brief append \f$rhslb\le Ay \le rhsub\f$ as linear constraints on the groundset variables \f$y\f$. \f$A\f$ has   @a append_n_rows new rows with coefficients given in  append_rows (if NULL, use default value),  lower bounds append_rhslb (if NULL use default) and upper bounds append_rhsub (if NULL use default)

     @param[in] append_n_rows  (Integer)
        (nonnegative) number of rows to be appended
  as linear constraints on the ground set

     @param[in] append_rows  (Sparsemat*)
        describes the coefficients of the linear constraints; the
  number of rows must match append_n_rows, the number of columns
  must match the current dimension of the groundset;
  if NULL, all coefficients are considered zero.

     @param[in] append_rhslb  (Matrix*)
        specifies lower bounds on the values of the constraints;
        the  number of rows must match append_n_rows;
  if NULL, all coefficients are considered CB_minus_infinity.

     @param[in] append_rhsub  (Matrix*)
        specifies upper bound on the values of the constraints;
        the  number of rows must match append_n_rows;
  if NULL, all coefficients are considered CB_plus_infinity.

     @return
        - 0 on success
        - != 0 otherwise

     */

    int append_constraints(CH_Matrix_Classes::Integer append_n_rows,
      const CH_Matrix_Classes::Sparsemat* append_rows = 0,
      const CH_Matrix_Classes::Matrix* append_rhslb = 0,
      const CH_Matrix_Classes::Matrix* append_rhsub = 0);

    //@}

    //------------------------------------------------------------
    /**@name Basic algorithmic routines and parameters */
    //@{

    /** @brief solves or does a prescribed number of iterations

  Bundle methods solve a problem by a sequence of so called
        descent steps that actually bring progress by moving from the
        current "center point" to a new center with better objective.
        A descent step may consist of several function evaluations
        (null steps), that lead to no immediate progress but mainly
  improve a cutting model of the objective function close to
        the current center point.  A minimizer to the model is
        accepted as descent step if the function value at this point
        satisfies a sufficient decrease criterion in comparison to the
        decrease predicted by the model. Having found a descent step,
        the next center is automatically shifted to this successful
        candidate.  Termination criteria may stop the process of
        seeking for a descent step, in which case the current center
        is kept and the routine termination_code() returns the
        termination code.

        Restarting, after each descent step, the bundle method from scratch
        with the new center as starting point does not endanger convergence.
        Therefore, a descent step is the smallest unit, after which
        user interaction can take place safely. To allow this there
  is a flag stop_at_descent_steps that will cause the code to
  return after the next descent step.

  If you know what your are doing, you may also use the input
        parameter maxsteps to force the algorithm to return after
        at most maxsteps null steps. Calling solve again
        without any intermediate problem configurations will then
        simply continue the process where it stopped and convergence
        is save. During null steps one may not decrease the weight
        or delete nonzero variables of the center or the current candidate!

        In a Lagrangean relaxation cutting plane approach one may want
        to separate and enlarge the dimension after a certain number
        of null steps. In this case the code will try to preserve the model,
        given appropriate subgradient extension routines have been
        provided. If the model cannot be extended, it has to be
        discarded (if subgradient extension is not successful
        this is done automatically), and the algorithm will be restarted
        from the current center point.

      @param[in] maxsteps (int)
          if maxsteps>0 the code returns after at most so many null steps

      @param[in] stop_at_descent_steps (int)
          if true the code also returns whenever a descent step occured

      @return
        - 0 on success
        - != 0 otherwise

    */
    int
      solve(int maxsteps = 0, bool stop_at_descent_steps = false);


    /** @brief Returns the termination code of the bundle algorithm for the latest descent step

      For resetting all counters relevant for termination see clear_fail_counts() .

      @return
      -  0  :    Not terminated.
             (Continue with the next solve())
      -  1  :    Relative precision criterion satisfied. (See set_term_relprec())
      -  2  :    Timelimit exceeded.
             (Currently the C interface does not offer a timelimit.)
      -  4  :    Maximum number of function reevaluations exceeded.
             (Indicates that there is a problem with one of the function
             oracles that seems to deliver no valid upper bounds on the true
             function value for descent steps)
      -  8  :    Maximum number of quadratic subproblem failures exceeded.
             (Indicates that the numerical limits of the inner quadratic
             programming solver are reached, no further progress expected)
      - 16  :    maximum number of model evaluation failures exceeded
             (Indicates that the numerical limits of the setup of the
             subproblem are reached, no further progress expected)
      - 32  :    maximum number of failures to increase the augmented model value exceeded
             (Indicates that the numerical limits  of the interplay between
             subproblem and quadratic programming solver are reached,
             no further progress expected)
       - 64  :   maximum number of oracle calls (function evaluations) exceeded,
                 see set_eval_limit()
       - 128  :   maximum number of oracle failures exceeded.
             This refers to function evaluations that terminate with insufficient
       precision but still provide a new approximate subgradient. A failure typically
             indicates numerical difficulties with the precision requirements.
             (Currently the interface does not allow to manipulate the limit, it is set to 10)


    */
    int
      termination_code() const;

    /** @brief Outputs a text version of termination code, see termination_code().

      @return
        - 0 on success
        - != 0 otherwise

    */
    std::ostream&
      print_termination_code(std::ostream& out);

    /** @brief Returns the objective value resulting from last descent
        step (initially undefined). If no problem modification routines
        were called since then, it is the objective value at the point
        returned by get_center().
    */
    double
      get_objval() const;

    /** @brief Returns the next center point that was produced by the latest call
  to solve (in some problem modification routines the
  center point may be updated immediately, in others the center point
  will be corrected automatically directly before starting
  the next descent step and its values may be infeasible till then).

      @return
        - 0 on success
        - != 0 otherwise
    */
    int
      get_center
      (CH_Matrix_Classes::Matrix& center) const;


    /** @brief Returns Euclidean norm of the latest aggregate subgradient.
     */
    double
      get_sgnorm() const;

    /** @brief Returns the latest aggregate subgradient (of the entire problem with groundset as provided by the solver)

    @return
      - 0 on success
      - != 0 otherwise

    */
    int
      get_subgradient
      (CH_Matrix_Classes::Matrix& subgradient) const;

    /** @brief Returns the cutting model value resulting from last call to
        solve() (initially undefined).
    */
    double
      get_cutval() const;

    /** @brief Returns the objective value computed in the last step of solve(),
        independent of whether this was a descent step or a null step (initially undefined).

        If no problem modification routines were called since then, it is the
        objective value at the point returned by get_candidate(). If this
        last evaluation led to a descent step, then it is the same value as
        in get_objval().
    */
    double
      get_candidate_value() const;

    /** @brief Returns the last point, the "candidate", at which the function
        was evaluated in solve().

        If this evaluation lead to a descent step, it is the same point as
  in get_center().

      @return
        - 0 on success
        - != 0 otherwise
    */
    int
      get_candidate
      (CH_Matrix_Classes::Matrix& center) const;


    //@}

    //------------------------------------------------------------
    /**@name Advanced algorithmic routines and parameters */
    //@{

    /** @brief Sets the relative precision requirements for successful termination
               (default 1e-5).

     @param[in] term_relprec (double)
       The algorithm stops with termination code 1, if predicted progress for
       the next step is less than term_relprec times
       absolute function value plus one.

    @return
      - 0 on success
      - != 0 otherwise

    */
    int
      set_term_relprec
      (const double term_relprec);

    /** @brief Set the starting point/center that will be used in the
        next call to  solve(). Each call
        to this routine causes an immediate evaluation of all oracles.

     @return
        - 0 on success
        - != 0 otherwise
    */
    int
      set_new_center_point
      (const CH_Matrix_Classes::Matrix& center_point);


    /** @brief Returns the return value of the latest evaluation call
        to this @a function.
    */
    int
      get_function_status
      (const FunctionObject& function) const;

    /** @brief Returns the multipliers for the box constraints on the design variables;
       in Lagrangean relaxation they may be interpreted as primal slacks
 for inequality constraints.
    @return
       - 0 on success
       - != 0 otherwise
   */
    int
      get_approximate_slacks(CH_Matrix_Classes::Matrix&) const;

    /** @brief returns the current approximate primal solution corresponding
        to the aggregate subgradient of the specified @a function.

        PrimalData solutions must have been supplied in all previous
        calls to evaluate; In this case it returns the current approximate
        primal solution aggregated alongside with the aggregate subgradient.
        A primal solution may not be available after addition of constraints,
        if extension of the aggregate subgradient to the new coordinates failed.
  If no primal data is availalbe, the function returns NULL.

     @return
        - pointer to the primal data of the aggregate of this function object
        - 0 if no primal is available
    */
    const PrimalData*
      get_approximate_primal
      (const FunctionObject& function)
      const;

    /** @brief Returns the primal solution corresponding to the best epsilon
  subgradient returned in the evaluation of the specified @a function
        at the current center point. If no primal data is availalbe,
  the function returns NULL.

     @return
        - pointer to the primal data of the minorant returned on evaluation
    of this function object at the current center
        - 0 if no primal is available
    */
    const PrimalData*
      get_center_primal
      (const FunctionObject& function)
      const;


    /** @brief Returns the primal solution corresponding to the best epsilon
  subgradient returned in the evaluation of the specified @a function
        at the point get_candidate. If no primal data is availalbe,
  the function returns NULL.

     @return
        - pointer to the primal data of the minorant returned on evaluation
    of this function object at the current candidate
        - 0 if no primal is available
    */
    const PrimalData*
      get_candidate_primal
      (const FunctionObject& function)
      const;


    /** @brief Starts/ends the use of a common SumBundle of the given bundle_size
        with a heuristic rule for selecting up to n_local_models in each bundle iteration

        If the function is the sum of many functions, having a local model for every one
        of them may result in a huge quadratic subproblem. It may then be better to
  form a common model of most of the functions, where a heuristic dynamically
  selects a few of the functions, for which a local model seems worth while.
  Whether such a common model should be used, how many subgradients it should contain,
  and how many local models are to be selected at most are the parameters
  set here.

  Setting these parameters only has an effect if bundle models of functions
  are present. If further functions are added later, the call should be repeated.

  This interface provides a simpler access to the SumBundle features by using
  some default parameter choices that could be set separately in
  set_bundle_parameters() and set_sumbundle_parameters() in a refined way.

     @param[in] use_sumbundle (bool)
        use value true to switch the sumbundle on, use value false to switch it off

     @param[in] n_local_models (int)
         upper bound on the number of local models to be used on top of the sumbundle's model,
   negative values correspond to no upper bound and all functions may have local models

     @param[in] bundle_parameters (const BundleParameters*)
         the maximum number of subgradients to be used in forming the SumBundle model,
   values <=1 are set to 2;

     @param[in] strategy (int)
         this is currently in experimental stage and allows to choose among some internal
         sumbundle strategies (currently 0,1,2,11 are available)

     @return
        - 0 on success
        - != 0 otherwise

    */
    int
      set_sumbundle
      (bool use_sumbundle,
        int n_local_models = -1,
        const BundleParameters* bundle_parameters = 0,
        int strategy = 1);

    /** @brief Sets the maximum number of subgradients used in forming the
        cutting model of the specified @a function

        Quite often a very small model, e.g., 2, yields very fast iterations
        and good progress in time (sometimes at the cost of more evaluations).
        By limited numerical experience, a significant reduction in the number of
        evaluations can  only be expected if the bundle is large enough to
        wrap the function rather tightly. Quite frequently, unfortunately,
        this entails that solving the quadratic subproblems
        is more expensive than function evaluation.

        The meaning of this routine may differ from standard for
  predefined special functions with special bundle types.

     @param[in] max_modelsize (int)
         maximum number of subgradients to be used in forming the cutting model

     @param[in] function
       if the aggregate subgradient of a particular function is desired,
       provide the pointer here, otherwise this referrs to the root function
       (if there is only one function to be optimized over, this is this single
       function, otherwise it is the sum of functions)

     @return
        - 0 on success
        - != 0 otherwise

    */
    int
      set_max_modelsize
      (int max_modelsize, const FunctionObject* function = 0);

    /** @brief Sets the maximum number of subgradients stored for use in
  forming the model or determining scaling information, it must be as
  least as large as max_modelsize (and is increased to this if not)

  The meaning of this routine may differ from standard for
  predefined special functions with special bundle types.


     @param[in] max_bundlesize (int)
       maximum number of subgradients stored for use in forming the model

     @param[in] function
       if the aggregate subgradient of a particular function is desired,
       provide the pointer here, otherwise this referrs to the root function
       (if there is only one function to be optimized over, this is this single
       function, otherwise it is the sum of functions)



     @return
       - 0 on success
       - != 0 otherwise

    */
    int
      set_max_bundlesize
      (int max_bundlesize, const FunctionObject* function = 0);

    /**@brief Sets the maximum bundlesize and the maximum number of new subgradients
        added in a bundle update of the cutting model for the specified @a function.
        The meaning of this routine may differ from standard for
  predefined special functions with special bundle types.

     @param[in] params (const BundleParameters&)
       some update parameters for the cutting model, see e.g. ConicBundle::BundleParameters

     @param[in] function
       if the aggregate subgradient of a particular function is desired,
       provide the pointer here, otherwise this referrs to the root function
       (if there is only one function to be optimized over, this is this single
       function, otherwise it is the sum of functions)


     @return
       - 0 on success
       - != 0 otherwise

    */
    int
      set_bundle_parameters
      (const BundleParameters& params, const FunctionObject* function = 0);

    /** @brief Retrieves current bundle parameters (not the actual size in use!)
       as set for the cutting model of the specified @a function.

 This may differ for predefined special
 functions with derived BundleParameter classes.

 If the code is asked to optimize over the sum of several functions,
 it usually does this with a separate model for each function. If there
 are too many function for this, it may be worth to consider using
 the SumBundle features. For this see also set_sumbundle_parameters().
 If the root function is a sum of functions, passing a
       SumModelParametersObject here allows to specify how many local models
 should be kept by SumModelParametersObject::set_max_local_models()
       and how these should be selected. A possible implementation for this
 is given in SumModelParameters.

    @param[in] function
      if the aggregate subgradient of a particular function is desired,
      provide the pointer here, otherwise this referrs to the root function
      (if there is only one function to be optimized over, this is this single
      function, otherwise it is the sum of functions)

    @return
      - 0 on success
      - != 0 otherwise

   */
    const BundleParameters*
      get_bundle_parameters
      (const FunctionObject* function = 0) const;


    /**@brief Specifies the behavior of the model (of the specified function)
        concerning requests to join or start a SumBundle that subsumes several
  models instead of providing a separate model for each funciton.

        The abstract interface for these Parameters is specified in
        SumBundleParametersObject, a concrete implementation is
        SumBundleParameters. Besides the usual BundleParameters
        the new main parameter is specified in
  SumBundleParametersObject::set_acceptable_mode(), see there.

     @param[in] params (const BundleParameters&)
       some update parameters for the cutting model, see e.g. ConicBundle::BundleParameters

     @param[in] function
       if the aggregate subgradient of a particular function is desired,
       provide the pointer here, otherwise this referrs to the root function
       (if there is only one function to be optimized over, this is this single
       function, otherwise it is the sum of functions)


     @return
       - 0 on success
       - != 0 otherwise

    */
    int
      set_sumbundle_parameters
      (const SumBundleParametersObject& params, const FunctionObject* function = 0);

    /** @brief Returns all current bundle data of the cutting
        model of the specified @a function.

  This may differ for predefined special
  functions with derived classes.

     @param[in] function
       if the aggregate subgradient of a particular function is desired,
       provide the pointer here, otherwise this referrs to the root function
       (if there is only one function to be optimized over, this is this single
       function, otherwise it is the sum of functions)

     @return
       - 0 on success
       - != 0 otherwise

    */
    const BundleData*
      get_bundle_data
      (const FunctionObject* function = 0) const;

    /** @brief Clears cutting model, subgradients and stored function values
       for the specified @a function (but only for the given one, not recursively)

       There should be no need to call this if the modification
       routines of this interface were used correctly. If, however,
       the oracle is modified by other means outside this interface,
       this has to be called whenever the specified function was
       modified so that the old subgradients and/or primal generators
       are no longer valid.

     @param[in] function
       if the aggregate subgradient of a particular function is desired,
       provide the pointer here, otherwise this referrs to the root function
       (if there is only one function to be optimized over, this is this single
       function, otherwise it is the sum of functions)

      @return
        - 0 on success
        - != 0 otherwise

    */
    int reinit_function_model
    (const FunctionObject* function = 0);

    /** @brief Clears the aggregate parts of the cutting model of this @a function
       (but only for the given one, not recursively)

       There should be no need to call this if the modification
       routines of this interface were used correctly. If, however,
       the oracle is modified by other means outside this interface,
       this has to be called whenever the specified function was
       modified so that the old aggregate subgradients and/or primal
       generators are no longer valid.

      @param[in] function (const FunctionObject&)
        the function added in add_function()

      @return
        - 0 on success
        - != 0 otherwise

    */
    int clear_aggregates
    (const FunctionObject* function = 0);

    /** @brief Asks @a function to call @a primal_extender for each of its primal objects (see
     also FunctionOracle::evaluate() )

     If the function is the Lagrangian dual of a primal problem and primal_data
     returned previous calls to the oracle has now to be updated due to changes
     in the primal problem -- e.g., this may happen in column generation -- the
     call causes updates of all internally stored primal_data objects by calling
     PrimalExtender::extend on each of these.

      @param[in] function (const FunctionObject&)
        the function added in add_function()

      @param[in] primal_extender (PrimalExtender&)
        the object holding the extension function for primal_data

      @return
        - 0 on success
        - 1 if for this function it is not possible to use a primal_extender
        - 2 if the primal_extender would be applicable but there is no primal_data
    */
    int call_primal_extender
    (const FunctionObject& function,
      PrimalExtender& primal_extender);

    /** @brief Returns the current weight for the quadratic term in the augmented subproblem
    (may be interpreted as 1./step_size or 1./trustregion-radius).
    */
    double
      get_last_weight() const;

    /** @brief Returns the next weight 	for the quadratic term in the augmented subproblem
        suggested by the internal weight updating heuristic
    */
    double
      get_next_weight() const;

    /** @brief Sets the  weight (>0) to be used in the quadratic term
        of the next augmented subproblem
        (may be interpreted as 1./step_size or 1./trustregion-radius).

        Independent of whether the weight violates current min- and max-bounds
        set in set_min_weight() and set_max_weight(), the next model will
        be computed for this value. Thereafter, however, it will be updated as
        usual; in particular, it may be truncated by min and max bounds
        immediately after the first subproblem.

        In order to guarantee a constant weight (e.g. 1 is frequently a reasonable
        choice if the automatic default heuristic performs poorly), set the min and max
         bounds to the same value, too.

      @param[in] weight (double)

      @return
        - 0 on success
        - != 0 otherwise
    */
    int
      set_next_weight
      (const double weight);

    /** @brief Sets a lower bound on the  weight for the quadratic term of the
        augmented subproblem.

        Nonpositive values indicate no bound.
        The new value shows its effect only at first dynamic change of
        the weight.

     @param[in] min_weight (double)

     @return
        - 0 on success
        - != 0 otherwise

    */
    int
      set_min_weight
      (const double min_weight);

    /** @brief Sets an upper bound on the  weight for the quadratic term of the
        augmented subproblem.

        Nonpositive values indicate no bound.
        The new value shows its effect only at first dynamic change of
        the weight.

      @param[in] max_weight (double)

      @return
        - 0 on success
        - != 0 otherwise
    */
    int
      set_max_weight
      (const double max_weight);

    /** @brief Replaces the internal update routine for choosing the weight used in the proximal term; input NULL reinstalls the default routine.

        The BundleWeight class instance pointed to will be deleted on
        construction, i.e., ownership is passe over to the solver.

      @param[in] bw
        replace internal update routine by bw, value 0 reinstalls the default routine

      @return
        - 0 on success
        - != 0 otherwise
    */
    int
      set_weight_update
      (BundleWeight* bw);

    /** @brief Adjusts on all conic functions the penalty parameter for
  conic violations to twice the trace of the primal approximation.

        This routine is only needed for conic function objects such
        as the nonnegative cone, the second order cone and
        the semidefinite cone if no good upper bound on the trace of
        feasible points is known and has to be determined automatically.

        If after some time, the trace values settle, the upper bounds
        on the trace may be way to high and can then be reset with this
        call.

      @return
        - 0 on success
        - != 0 otherwise
    */
    int
      adjust_multiplier
      (void);


    /** @brief Use a variable metric heuristic or switch off general metrics alltogether.
       (variable metric resets the quadratic term e.g. to some diagonal matrix,
       switching it off resets the quadratic term to the identity times the weight)

     @param[in] do_variable_metric (int)
        - 0 switch off the scaling heuristic
        - 1 use a diagonal scaling heuristic
        - 2 use a diagonal scaling heuristic combined with one for the bounds
        - 3 use a low rank scaling heuristic
        - 4 use a low rank scaling heuristic combined with a diagonal term
        - 5 use a dense scaling heuristic

      @return
        - 0 on success
        - != 0 otherwise
    */
    int
      set_variable_metric
      (int do_variable_metric);


    /** @brief For variable metric install the BundleProxObject pointed to; the object is passed to the solver who will delete it on termination or when replaced

      @param[in] proxp (BundleProxObject*)
          replace the current BundleProxObject by this object on the heap;
    NULL is allowed and results in the default choice;
    the object pointed to will be deleted by the solver

      @return
        - 0 on success
        - != 0 otherwise
    */
    int
      set_prox
      (BundleProxObject* proxp);

    /** @brief If set to true (the default is false), variables may be
        fixed automatically to active bounds if these are strongly
        active (i.e., the corresponding multipliers are big) and
  the center values are also right on these bounds already.

        The coordinates to be fixed are redetermined in each
        call following a descent step or a change of the function.
        An indicator vector of the variables fixed during the last call
        can be obtained via the routine get_fixed_active_bounds().

        Setting this value to true might improve the performance
        of the algorithm in some instances but there is no
        convergence theory. It might be particularly helpful
        within Lagrangian relaxation if a primal cutting plane
        approach is used and non-tight inequalities should be
        eliminated quickly (fixing then indicates large primal
        slack values as these are the dual variables to the bounds
  on the Lagrange mulitpliers). Furthermore, if the value
  of a variable is fixed to zero, the variable can typically
  be deleted without affecting the validity of the current
  cutting model and function values.

      @param[in] allow_fixing (bool)

      @return
        - 0 on success
        - != 0 otherwise
     */
    void
      set_active_bounds_fixing
      (bool allow_fixing);

    /** @brief clears all fail counts on numerical function oder model failures,
  may be useful if this caused premature termination.

      @return
        - 0 on success
        - != 0 otherwise
    */
    void
      clear_fail_counts
      (void);

    /** @brief Sets an upper bound on the number of calls to the oracle (use negative numbers for no limit).

        If this number is reached, the algorithm will terminate
  independently of whether the last step was a descent or
  a null step. A negative number will be interepreted as
  no limit.

      @param[in] eval_limit (Integer)

      @return
        - 0 on success
        - != 0 otherwise
    */
    void
      set_eval_limit
      (CH_Matrix_Classes::Integer eval_limit);

    /** @brief Set an upper bound on the number of inner updates for the
        cutting model with primal slacks within one null step (use negative numbers for no limit).

        A negative number will be interepreted as no limit, i.e.,
        the updates will be done till a certain precision of the
        cutting model is achieved.

      @param[in] update_limit (Integer)

      @return
        - 0 on success
        - != 0 otherwise
    */
    void
      set_inner_update_limit
      (CH_Matrix_Classes::Integer update_limit);

    /** @brief Set an upper bound on the number of seconds (user time, use negative numbers for no limit)

      @param[in] time_limit (Integer)

      @return
        - 0 on success
        - != 0 otherwise
    */
    void
      set_time_limit
      (CH_Matrix_Classes::Integer time_limit);

    /* * @brief Set parameters for the internal QP solver, possibly after first exchanging the solver with a new one

      The objects passed need to be heap objects; their ownership is transferred
      to the bundle code and they will be deleted there. The pointers may be null,
      in which case nothing is done for this object

      @param[in] qpparams (QPSolverParametersObject*)

      @param[in] newqpsolver (QPSolverObject*)

      @return
        - 0 on success
        - != 0 otherwise
    */

    int set_qp_solver(QPSolverParametersObject* qpparams,
      QPSolverObject* newqpsolver = 0);

    //@}

    //------------------------------------------------------------
    /**@name Look up basic paramaters (dimension, number of functions, ...)*/
    //@{
    /** @brief Returns the current dimension of the design space/argument
           or -1 if no dimension is set.
    */
    int get_dim() const;

    /** @brief Returns the current number of functions in the problem.
    */
    int get_n_functions() const;

    /** @brief Returns the number of function evaluations
    */
    int get_n_oracle_calls() const;

    /** @brief Returns the number of function descent setps
    */
    int get_n_descent_steps() const;

    /** @brief Returns the number of inner iterations of the bundle method
    */
    int get_n_inner_iterations() const;

    /** @brief Returns the number of inner multiplier updates for the box constraints
    */
    int get_n_inner_updates() const;

    /** @brief returns true if the last evaluation of the last call to solve()
  resulted in a descent step

  Mind: if there was no (succesdful) evaluation, neither get_descent_step() nor
  get_null_step() will return true;
    */
    bool get_descent_step() const;

    /** @brief returns true if the last evaluation of the last call to solve()
  resulted in a null step

  Mind: if there was no (successful) evaluation, neither get_descent_step() nor
  get_null_step() will return true;
    */
    bool get_null_step() const;

    /** @brief If a linear cost vector was specified, costs will hold these values, otherwise the vector is initialized to zero (for the current dimension)
    */
    int get_costs(CH_Matrix_Classes::Matrix& costs) const;

    /** @brief Returns a pointer to the vector of lower bounds or null if there is no such vector
    */
    const CH_Matrix_Classes::Matrix* get_lbounds() const;

    /** @brief Returns a pointer to the vector of upper bounds or null if there is no such vector
    */
    const CH_Matrix_Classes::Matrix* get_ubounds() const;

    /** @brief Returns NULL or (iff active bound fixing is turned on
  in set_active_bounds_fixing()) the indicator vector of
  variables temporarily fixed to the center value due to
  significantly positive multipliers for the box constraints.

        A variable gets fixed to the bound only if the center is
  already a the bound and in some iteration the dual variables
  to the bound constraint indicate that the bound is strongly
  active also for the candidate. Of course this migh just hold
  for one candidate and there is no guarantee that the bound
  is also strongly active in an optimal solution. Thus, this
  mainly a heuristic to eliminate less important variables
  quickly from entering the subproblem.
     */
    const CH_Matrix_Classes::Indexmatrix* get_fixed_active_bounds() const;

    /** @brief Returns the pointer to the current prox term of the bundle solver
    */
    BundleProxObject* get_prox() const;

    /// returns true if for function some modifications are pending, the arguments give information on these
    bool
      pending_oracle_modification(const FunctionObject& function,
        CH_Matrix_Classes::Integer& old_dim,
        CH_Matrix_Classes::Integer& new_dim,
        CH_Matrix_Classes::Integer& append_dim,
        const CH_Matrix_Classes::Indexmatrix*& map_to_old,
        const CH_Matrix_Classes::Indexmatrix*& deleted_indices,
        const CH_Matrix_Classes::Indexmatrix*& new_indices,
        const OracleModification*& oracle_modification
      ) const;

    //@}

    //------------------------------------------------------------
    /**@name Output */
    //@{

    /** @brief Specifies the output level (out==NULL: no output at all,
           out!=NULL and level=0: errors and warnings,
           level>0 increasingly detailed information)

     @param[in] out  (std::ostream*)
       direct all output to (*out). If out==NULL, there will be no output at all.

     @param[in] print_level (int)

     Output levels for print_level:
      -  0 ... no output except for errors and warnings
      -  1 ... line summary after each descent step
      - >1 ... undocumented and increasingly detailed log information.
             These higher levels should only be used if requested
             for debugging purposes.

      Example for level 1:

\verbatim
00:00:00.00 endit  1   1   1   563.    563.  39041.188  39043.162
00:00:00.00 endit  2   2   2   563.    559.  38488.165  38490.200
00:00:00.00 endit  3   3   3   56.3    555.  33014.533  33211.856
00:00:00.00 endit  4   4   4   5.63    517. -14306.459  2738.0343
00:00:00.00 endit  5   5   5   4.04    148. -2692.1131  2.2150883
00:00:00.00 endit  6   6   6   4.01    1.29  1.7908952  2.0000581
00:00:00.00 endit  7   7   7   3.95  0.0213  1.9999387  2.0000000
00:00:00.00 _endit  8   8   8   3.95 2.94e-05  2.0000000  2.0000000

Column 1      2     3   4   5    6       7       8          9
\endverbatim
      - Column 1: computation time in hh:mm:ss.dd,
      - Column 2: "endit" is convenient for grep and stands for "end of iteration".
         Iterations with termination_code()!=0 are marked with "_endit".
      - Column 3: number of descent steps
      - Column 4: number of descent and null steps. Up to initialization calls
         and reevaluations, this is the number of evaluation calls
         to the function oracles from within the bundle method.
         In the example all calls led to descent steps.
      - Column 5: number of innermost iterations. It differs from column 5 only in the
          case of variables with bounds in which case it gives the number of updates
          of the multipliers for the bounds (or primal slacks in Lagrangean
          relaxation). Exceedingly high numbers in this column indicate that
          some variables are constantly at their bounds and it might be
          possible to improve convergence by deleting them (i.e. set them
          as constants to this bound and remove the variable).
      - Column 6: the weight of the quadratic term in the augmented problem.
      - Column 7: the norm of the aggregate subgradient. If it is small,
          say below 0.1, then mostly this is good indication that the
          objective value is close to optimal.
      - Column 8: the value of the cutting model in the last candidate point. It
          is always a lower bound on the true function value in this point
      - Column 9: the objective value in the latest point that led to a descent
          step, i.e., the point returend by get_center(). Whenever
          termination_code() returns 0 this is also the objective
          value of the latest evaluation call to the function oracles and
          the value in the center point of the next iteration.
     */
    void
      set_out
      (std::ostream* out = 0, int print_level = 1);

    /// print a one line summary of important evaluation data
    std::ostream& print_line_summary(std::ostream& out) const;

    /// print a cryptic summary of computation times of important components 
    std::ostream& print_statistics(std::ostream& out) const;

    //@}

  };

  //@}
}

#endif



