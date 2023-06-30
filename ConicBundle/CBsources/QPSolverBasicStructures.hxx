/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/QPSolverBasicStructures.hxx
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


#ifndef CONICBUNDLE_QPSOLVERBASICSTRUCTURES_HXX
#define CONICBUNDLE_QPSOLVERBASICSTRUCTURES_HXX

/**  @file QPSolverBasicStructures.hxx
    @brief Header declaring the classes ConicBundle::QPSolverBasicInterface, ConicBundle::CentralPathPoint, ConicBundle::QPSolverBasicStructures
    @version 1.0
    @date 2020-03-18
    @author Christoph Helmberg
*/



#include <vector>
#include "clock.hxx"
#include "QPModelBlock.hxx"
#include "QPSolverParameters.hxx"
#include "SOCIPProxBlock.hxx"

namespace ConicBundle {


  /** @ingroup ConstrainedQPSolver
  */
  //@{

   /** @brief defines the abstract interface for QPSolverBasicStructures
       that the basic routines there need to access the cost and
       constraint data; it also defines the interface to
       QPSolverBasicStructures, how to call those routines and
       how to retrieve the results.

       Ignoring the model part, the routines solve
       \f$min_{x\in\{l\le x\le u:rhsl<=Ax<=rhsu\}}frac12x^TQx+c^Tx+\gamma\f$

       NOTE: as this is taken in part from an independent code, the roles
       of y and x switch here in comparison to their standard use in ConicBundle.

       Q needs to be positive definite.
       Q and A are supposed to be available by matrix-vector multiplication and
       the diagonal of Q must be extractable for use in a preconditioner.
       If the bound constraints are few, they still must be provided in full length,
       but they may be acessed by indexvectors.

       \f$rhsl<=Ax<=rhsu\f$ is implemented by introducing an additional variable \f$s\f$
       with \f$Ax+s=0\f$ and \f$rhsl<=-s<=rhsu\f$. If in certain componnets
       the lower and upper bounds are identical, \f$s\f$ is fixed to this value.

       The interpretation of the dual variables may be read of from the Lagrangian
       \f$L(x,s,y,z_l,z_u,z_{rl},z_{ru})=\frac12x^TQx+c^Tx+\gamma+(Ax+s)^Ty+(l-x)^Tz_l+(x-u)^Tz_u+(s+rhsl)^Tz_{rl}+(-s-rhsu)^Tz_{ru}\f$

       The main entry point is QPIsolve() where I is short for internal; see
       the general explanations in QPSolverBasicStructures for an outline of
       the inner workings thereof.
   */

  class QPSolverBasicInterface : public virtual CBout {
  public:
    /// virtual destructor
    virtual ~QPSolverBasicInterface();

    //------------------------------------------------------------
    /**@name interface for QPSolverBasicStructures for getting the data */
    //@{

    /// dimension of the quadratic variable x
    virtual CH_Matrix_Classes::Integer QPget_xdim() const = 0;

    /// dimension of the dual variables to the linear constraints (number of rows of A)
    virtual CH_Matrix_Classes::Integer QPget_ydim() const = 0;

    /// vecout += quadratic term * vecin
    virtual CH_Matrix_Classes::Matrix& QPadd_Qx(const CH_Matrix_Classes::Matrix& vecin, CH_Matrix_Classes::Matrix& vecout) const = 0;

    /// outplusAx += A * xin
    virtual CH_Matrix_Classes::Matrix& QPadd_Ax(const CH_Matrix_Classes::Matrix& xin, CH_Matrix_Classes::Matrix& outplusAx) const = 0;

    /// outplusAty += transpose(A) * yin
    virtual CH_Matrix_Classes::Matrix& QPadd_Aty(const CH_Matrix_Classes::Matrix& yin, CH_Matrix_Classes::Matrix& outplusAty) const = 0;

    /// returns the linear cost vector for the quadratic variable x
    virtual const CH_Matrix_Classes::Matrix& QPget_c() const = 0;

    /// returns a potential constant offset of the cost function
    virtual CH_Matrix_Classes::Real QPget_gamma() const = 0;

    /// returns rhslb of    rhslb <= Ax
    virtual const CH_Matrix_Classes::Matrix& QPget_rhslb() const = 0;

    /// returns rhsub of    Ax <= rhsub
    virtual const CH_Matrix_Classes::Matrix& QPget_rhsub() const = 0;

    /// returns the indices where rhslb is not minus infinity  (sorted increasingly)
    virtual const CH_Matrix_Classes::Indexmatrix& QPget_rhslbind() const = 0;

    /// returns the indices where rhsub is not plus infinity (sorted increasingly)
    virtual const CH_Matrix_Classes::Indexmatrix& QPget_rhsubind() const = 0;

    /// returns lb of   lb <= x
    virtual const CH_Matrix_Classes::Matrix& QPget_lb() const = 0;

    /// returns ub of   x <= ub
    virtual const CH_Matrix_Classes::Matrix& QPget_ub() const = 0;

    /// returns the indices where lb is not minus infinity (sorted increasingly)
    virtual const CH_Matrix_Classes::Indexmatrix& QPget_lbind() const = 0;

    /// returns the indices where ub is not plus infinity (sorted increasingly)
    virtual const CH_Matrix_Classes::Indexmatrix& QPget_ubind() const = 0;
    //@}

    //------------------------------------------------------------
    /**@name interface for calling the routines and querying the result of QPSolverBasicStructures
     */
     //@{

     /// reset all values of the internal basic structures and variables
    virtual void QPIclear() = 0;

    /// pass on new parameters, ownership is tranferred to this, will be deleted here
    virtual int QPset_solver_parameters(QPSolverParameters* params) = 0;

    /// set a starting value for the quadratic variables
    virtual int QPset_startx(const CH_Matrix_Classes::Matrix& startx) = 0;

    /// set a starting value for the slack variables of the constraints
    virtual int QPset_starts(const CH_Matrix_Classes::Matrix& starts) = 0;

    /// set a starting value for the dual variables of the constraints
    virtual int QPset_starty(const CH_Matrix_Classes::Matrix& starty) = 0;

    /// set a starting value for the dual variables of the lower bounds on the quadratic variables
    virtual int QPset_startzlb(const CH_Matrix_Classes::Matrix& startlb) = 0;

    /// set a starting value for the dual variables of the upper bounds on the quadratic variables
    virtual int QPset_startzub(const CH_Matrix_Classes::Matrix& startub) = 0;

    /// set a starting value for the barrier parameter mu
    virtual int QPset_startmu(CH_Matrix_Classes::Real startmu) = 0;

    /// solve the problem to the precision specified by the parameters; if the problem is resolved for slightly modified cost data (the feasible set should not change), setting reinitialize=false and a skip_factor for previous "central path" solutions might help to speed up the resolve
    virtual int QPIsolve(bool reinitialize = true, CH_Matrix_Classes::Real skip_factor = -1.) = 0;

    /// after QPIsolve, retrieve the primal objective value (of the quadratic variables)
    virtual CH_Matrix_Classes::Real QPget_primalval() const = 0;

    /// after QPIsolve, retrieve the dual objective value (quadratic and dual variables)
    virtual CH_Matrix_Classes::Real QPget_dualval() const = 0;

    /// return the current parameters for potential modification (do not delete!)
    virtual QPSolverParameters* QPget_parameters() const = 0;

    /// return the current value of the quadratic variables (after a successful QPIsolve the optimal solution)
    virtual const CH_Matrix_Classes::Matrix& QPget_x() const = 0;

    /// return the current value of the dual variables to the constraints (after a successful QPIsolve the optimal solution)
    virtual const CH_Matrix_Classes::Matrix& QPget_y() const = 0;

    /// return the current value of the dual variables to the lower bounds of the quadratic variables (after a successful QPIsolve the optimal solution)
    virtual const CH_Matrix_Classes::Matrix& QPget_zlb() const = 0;

    /// return the current value of the dual variables to the upper bounds of the quadratic variables (after a successful QPIsolve the optimal solution)
    virtual const CH_Matrix_Classes::Matrix& QPget_zub() const = 0;

    /// return the current value of the slack variables of the constraints (after a successful QPIsolve the optimal solution)
    virtual const CH_Matrix_Classes::Matrix& QPget_s() const = 0;

    /// return the current value of the dual variables to the lower bounds of the constraints (after a successful QPIsolve the optimal solution)
    virtual const CH_Matrix_Classes::Matrix& QPget_rhszlb() const = 0;

    /// return the current value of the dual variables to the upper bounds of the constraints (after a successful QPIsolve the optimal solution)
    virtual const CH_Matrix_Classes::Matrix& QPget_rhszub() const = 0;

    /// return the current value of the barrier parameter mu 
    virtual CH_Matrix_Classes::Real QPget_mu() const = 0;

    //@}

  };


  /** @brief currently not in use, storing the points of the central path might help to restart faster if only the cost terms are modified slightly; an instance of this class stores one such point
  */

  class QPCentralPathPoint {
  public:
    /// quadratic variable
    CH_Matrix_Classes::Matrix x;

    /// dual variable
    CH_Matrix_Classes::Matrix y;

    /// constraint slacks
    CH_Matrix_Classes::Matrix s;

    /// dual to quadratic variable lower bound
    CH_Matrix_Classes::Matrix zlb;

    /// dual to quadratic variable upper bound
    CH_Matrix_Classes::Matrix zub;

    /// dual to lower bound on constraints
    CH_Matrix_Classes::Matrix rhszlb;

    /// dual to upper bound on constraints
    CH_Matrix_Classes::Matrix rhszub;

    /// barrier parameter
    CH_Matrix_Classes::Real mu;

    /// assignment operator
    QPCentralPathPoint& operator=(const QPCentralPathPoint& ppp) {
      x = ppp.x;
      y = ppp.y;
      s = ppp.s;
      zlb = ppp.zlb;
      zub = ppp.zub;
      rhszlb = ppp.rhszlb;
      rhszub = ppp.rhszub;
      mu = ppp.mu;
      return *this;
    }

    /// default constructor
    QPCentralPathPoint() {
    }

    /// copy constructor
    QPCentralPathPoint(const QPCentralPathPoint& ppp) {
      *this = ppp;
    }

    /// initialization constructor
    QPCentralPathPoint(
      const CH_Matrix_Classes::Matrix& ix,
      const CH_Matrix_Classes::Matrix& iy,
      const CH_Matrix_Classes::Matrix& is,
      const CH_Matrix_Classes::Matrix& izlb,
      const CH_Matrix_Classes::Matrix& izub,
      const CH_Matrix_Classes::Matrix& irhszlb,
      const CH_Matrix_Classes::Matrix& irhszub,
      const CH_Matrix_Classes::Real imu
    ) :
      x(ix), y(iy), s(is), zlb(izlb), zub(izub), rhszlb(irhszlb), rhszub(irhszub), mu(imu) {
    }

  };



  /** @brief provides the basic variables and implements basic routines
      for the primal dual interior point solver interface of
      QPSolverBasicInterface (see there for the problem description and
      the switched roles of x and y) except for setting up and solving
      the primal dual KKT system, which is left to a specialized
      QPKKTSolverObject provided in QPSolverParameters
      to allow the exploitation of structural properties in the
      data. The basic problem description must also still be provided in
      a derived class that implements the remaining abstract data
      retrieval functions. The bundle data is made available by deriving
      this class from QPModelPointer, see \ref InternalQPModelData.

      The entry point is QPIsolve() (the I is short for internal), which
      intializes the variables, sets the starting point, initializes the
      QPKKTSolverObject specified by QPSolverParameters via
      QPKKTSolverObject::QPinit_KKTdata() and calls QPiterate(). This
      routine then computes for the current point the primal and dual
      violations, the opimality gap, etc.  As long as the termination
      criterion steered by the QPSolverParameters is not met, it
      computes predictor corrector steps by calling
      QPpred_corr_step(). The latter sets up a primal dual KKT system
      via a call to QPKKTSolverObject::QPinit_KKTsystem() and
      then solves this for predictor and corrector via calls
      to QPKKTSolverObject::QPsolve_KKTsystem() with intermediate/susbsequent
      calls to QPlinesearch() and QPselect_mu() for determining the
      next step sizes and barrier parameter.

      The most important routines of the model described in the QPModelBlockObject
      that are required in this process are

      - QPModelBlockObject::get_mu_info() for selecting the barrier parameter

      - QPModelBlockObject::linesearch() for performing the line search of the model variables

      - QPModelBlockObject::do_step() for carrying out the current step with the given step length

      - QPModelBlockObject::globalx_cost() for retrieving the cost contribution of the model for the current design variable values

      - QPModelBlockObject::constraints_cost() for retrieving the cost contribution of the model constraints

      - QPModelBlockObject::add_Bt_modelx() for adding the Bundle information times model variables (= the aggregate) to the dual feasibility constraint

      - QPModelBlockObject::dualviol_2normsqr() for retrieving the Euclidean norm of the model's dual violation

      - QPModelBlockObject::primalviol_2normsqr() for retrieving the Euclidean norm of the model's primal violation

      - QPModelBlockObject::reset_starting_point() for initializing the model variables to a feasible starting point with respect to the given design variable values



   */


  class QPSolverBasicStructures : public QPSolverBasicInterface, public QPModelPointer, public virtual QPSolverObject {
  private:
    QPSolverParameters* paramsp;   ///< parameters giving termination bounds etc

    CH_Matrix_Classes::Matrix x;          ///< quadratic (primal) variables 
    CH_Matrix_Classes::Matrix y;          ///< dual variables to the affine constraints
    CH_Matrix_Classes::Matrix s;          ///< slack variables to the affine constraints
    CH_Matrix_Classes::Matrix zlb;        ///< dual variables to the lower bounds on the primal variables
    CH_Matrix_Classes::Matrix zub;        ///< dual variables to the upper bounds on the primal variables
    CH_Matrix_Classes::Matrix rhszlb;     ///< dual variables to the lower bounds on the constraints
    CH_Matrix_Classes::Matrix rhszub;     ///< dual variables to the upper bounds on the constraints

    CH_Matrix_Classes::Real mu;           ///< barrier parameter

    CH_Matrix_Classes::Real last_mu;      ///< barrier parameter of the previous step
    CH_Matrix_Classes::Real next_mu;      ///< barrier parameter computed for the next step if use_predictor_corrector==false

    CH_Matrix_Classes::Real last_theta;      ///< neighborhood parameter of the previous step
    CH_Matrix_Classes::Real next_theta;      ///< neighborhood parameter computed for the next step if use_predictor_corrector==false

    CH_Matrix_Classes::Real last_alpha;   ///< step size of the previous step

    CH_Matrix_Classes::Real sigma;        ///< current reduction factor used for mu

    CH_Matrix_Classes::Integer iter;      ///< counts the number of iterations/steps
    CH_Matrix_Classes::Integer large_predictor_cnt;  ///< increased if predictor promises a good step but mu is only decreased by a little

    std::vector<QPCentralPathPoint> central_path;  ///< stores the sequence of points, potentially useful for restarting

    CH_Matrix_Classes::Real primalval;      ///< primal objective value (if feasible)
    CH_Matrix_Classes::Real dualval;        ///< dual objective value (if feasible)

    CH_Matrix_Classes::Matrix primalviol;     ///< primal violation vector
    CH_Matrix_Classes::Matrix dualviol;       ///< dual violation vector
    CH_Matrix_Classes::Matrix yviol;          ///< violation between dual varialbes and duals to constraint bounds
    CH_Matrix_Classes::Real n2primalviol;       ///< Euclidean norm of primal violation
    CH_Matrix_Classes::Real n2dualviol;         ///< Euclidean norm of dual violation
    CH_Matrix_Classes::Real n2modelprimalviol;  ///< Euclidean norm of primal violation of the bundle model 
    CH_Matrix_Classes::Real n2modeldualviol;    ///< Euclidean norm of dual violation of the bundel model
    CH_Matrix_Classes::Real n2yviol;             ///< Euclidean norm of yviol

    CH_Matrix_Classes::Matrix shortzlb;    ///< lb duals condensed on indices with finite variable bounds
    CH_Matrix_Classes::Matrix shortzub;    ///< ub duals condensed on indices with finite variable bounds
    CH_Matrix_Classes::Matrix slacklb;     ///< slack to variable lower bounds
    CH_Matrix_Classes::Matrix slackub;     ///< slack to variable upper bounds

    CH_Matrix_Classes::Matrix shortrhszlb; ///< lb duals condensed on indices with finite constraint bounds
    CH_Matrix_Classes::Matrix shortrhszub; ///< ub duals condensed on indices with finite constraint bounds
    CH_Matrix_Classes::Matrix rhsslacklb;  ///< slack to constraint lower bounds 
    CH_Matrix_Classes::Matrix rhsslackub;  ///< slack to constraint upper bounds 

    CH_Matrix_Classes::Matrix dshortzlb;    ///< step for shortslb
    CH_Matrix_Classes::Matrix dshortzub;    ///< step for shortzub
    CH_Matrix_Classes::Matrix dslacklb;     ///< step for dslacklb
    CH_Matrix_Classes::Matrix dslackub;     ///< step for dslackub

    CH_Matrix_Classes::Matrix dshortrhszlb; ///< step for shortrhszlb
    CH_Matrix_Classes::Matrix dshortrhszub; ///< step for shortrhszub
    CH_Matrix_Classes::Matrix drhsslacklb;  ///< step for rhsslacklb
    CH_Matrix_Classes::Matrix drhsslackub;  ///< step for rhsslackub

    CH_Matrix_Classes::Matrix solx;  ///< last solution to the KKT system (primal Newton step)
    CH_Matrix_Classes::Matrix soly;  ///< last solution to the KKT system (dual Newton step)

    CH_Matrix_Classes::Matrix old_x;          ///< old quadratic (primal) variables 
    CH_Matrix_Classes::Matrix old_y;          ///< old dual variables to the affine constraints
    CH_Matrix_Classes::Matrix old_s;          ///< old slack variables to the affine constraints
    CH_Matrix_Classes::Matrix old_zlb;        ///< old dual variables to the lower bounds on the primal variables
    CH_Matrix_Classes::Matrix old_zub;        ///< old dual variables to the upper bounds on the primal variables
    CH_Matrix_Classes::Matrix old_rhszlb;     ///< old dual variables to the lower bounds on the constraints
    CH_Matrix_Classes::Matrix old_rhszub;     ///< old dual variables to the upper bounds on the constraints
    CH_Matrix_Classes::Real old_mu; ///< old barrier parameter



    SOCIPProxBlock socqp; ///< holds the second order cone model of the quadratic term if this is to be used according to the parameter settings (still with numerical difficulties)

    /// compute function values and violation; if init==true, this is the first time
    int QPcompute_values(bool init);

    /// line search for keeping a vector nonnegative
    void QPvec_linesearch(CH_Matrix_Classes::Real& alpha,
      const CH_Matrix_Classes::Matrix& vec,
      const CH_Matrix_Classes::Matrix& dvec) const;

    /// line search on all vectors, also calling the model line search if present
    void QPlinesearch(CH_Matrix_Classes::Real& alpha);

    /// choose the next value of the barrier parameter
    virtual int QPselect_mu(CH_Matrix_Classes::Real& mu,
      CH_Matrix_Classes::Real alpha);

    /// choose the next step size and value of the barrier parameter within the neighborhood
    virtual int QPselect_step_mu(CH_Matrix_Classes::Real& mu,
      CH_Matrix_Classes::Real& alpha,
      bool centering);

    /// choose the next step size by the local neighborhoods and then the next barrier paramater  
    virtual int QPselect_localstep_mu(CH_Matrix_Classes::Real& mu,
      CH_Matrix_Classes::Real& alpha,
      bool centering);

    /// predictor corrector step
    int QPpredcorr_step(CH_Matrix_Classes::Real& alpha,   ///< output value step size  
      CH_Matrix_Classes::Real& prec,    ///< in/out  precision 
      CH_Matrix_Classes::Real& next_mu, ///< in/out  barrier parameter
      bool use_predcorr,                ///< if true, use predictor to find next_mu
      bool use_nbh,                     ///< if true, use neighborhood for step and next_mu
      bool centering                    ///< if true, go for a feasible centered point for the current (or maybe even larger) mu 
    );

    /// call QPpredcorr_step repeatedly until termination
    int QPiterate();

  protected:

    //Statistics
    mutable CH_Tools::Clock clock;  ///< for timing purposes
    mutable CH_Tools::Microseconds QPcoeff_time;  ///< time spent in computing the QP coefficients
    mutable CH_Tools::Microseconds QPsolve_time; ///< time spent in solving the QP
    mutable CH_Tools::Microseconds QPinitKKT_time; ///< time spent initializing the KKT system
    mutable CH_Tools::Microseconds QPpredictor_time; ///< time spent in computing predictor 
    mutable CH_Tools::Microseconds QPcorrector_time; ///< time spent in computing predictor 
    mutable CH_Tools::Microseconds QPlinesearch_time; ///< time spent in the line searches
    mutable CH_Tools::Microseconds QPcurrentprecprep_time; ///< time spent in preparing the preconditioner
    mutable CH_Tools::Microseconds QPprecprep_time; ///< time spent in preparing the preconditioner
    mutable CH_Tools::Microseconds QPprecprepmodel_time; ///< time spent in preparing the preconditioner
    mutable CH_Tools::Microseconds QPprecsolve_time; ///< time spent in preconditioning solves
    mutable CH_Tools::Microseconds QPmatmult_time; ///< time spent in matrix vector multiplications
    mutable CH_Tools::Microseconds QP_time; ///< time spent in preparing the preconditioner

  public:
    /// default constructor
    QPSolverBasicStructures(QPSolverParameters* params = 0, CBout* cb = 0);
    /// virtual destructor, deletes the parameters
    virtual ~QPSolverBasicStructures() {
      delete paramsp; paramsp = 0;
    }


    /// reset all values of the internal basic structures and variables
    virtual void QPIclear();

    /// pass on new parameters, ownership is tranferred to this, will be deleted here
    virtual int QPset_solver_parameters(QPSolverParameters* params) {
      assert(params); if (params != 0) {
        delete paramsp; paramsp = params;
      } return 0;
    }

    /// set a starting value for the quadratic variables
    virtual int QPset_startx(const CH_Matrix_Classes::Matrix& startx) {
      assert(startx.dim() == QPget_xdim()); x = startx; return 0;
    }

    /// set a starting value for the slack variables of the constraints
    virtual int QPset_starts(const CH_Matrix_Classes::Matrix& starts) {
      assert(starts.dim() == QPget_ydim()); s = starts; return 0;
    }

    /// set a starting value for the dual variables of the constraints
    virtual int QPset_starty(const CH_Matrix_Classes::Matrix& starty) {
      assert(starty.dim() == QPget_ydim()); y = starty; return 0;
    }

    /// set a starting value for the dual variables of the lower bounds on the quadratic variables
    virtual int QPset_startzlb(const CH_Matrix_Classes::Matrix& startzlb) {
      assert(startzlb.dim() == QPget_xdim()); zlb = startzlb; return 0;
    }

    /// set a starting value for the dual variables of the upper bounds on the quadratic variables
    virtual int QPset_startzub(const CH_Matrix_Classes::Matrix& startzub) {
      assert(startzub.dim() == QPget_xdim()); zub = startzub; return 0;
    }

    /// set a starting value for the barrier parameter mu
    virtual int QPset_startmu(CH_Matrix_Classes::Real startmu) {
      mu = startmu; return 0;
    }

    /// solve the problem to the precision specified by the parameters; if the problem is resolved for slightly modified cost data (the feasible set should not change), setting reinitialize=false and a skip_factor for previous "central path" solutions might help to speed up the resolve
    virtual int QPIsolve(bool reinitialize = true, CH_Matrix_Classes::Real skip_factor = -1.);

    /// after QPIsolve, retrieve the primal objective value (of the quadratic variables)
    virtual CH_Matrix_Classes::Real QPget_primalval() const {
      return primalval;
    }

    /// after QPIsolve, retrieve the dual objective value (quadratic and dual variables)
    virtual CH_Matrix_Classes::Real QPget_dualval() const {
      return dualval;
    }

    /// return the current parameters for potential modification (do not delete!)
    virtual QPSolverParameters* QPget_parameters() const {
      return paramsp;
    }

    /// return the current value of the quadratic variables (after a successful QPIsolve the optimal solution)
    virtual const CH_Matrix_Classes::Matrix& QPget_x() const {
      return x;
    }

    /// return the current value of the quadratic variables (after a successful QPIsolve the optimal solution) together with an activity estimate (1 ... active, 0 ... at bounds); in xout inactive components are set to the respective bound
    virtual int QPget_x(CH_Matrix_Classes::Matrix& xout,
      CH_Matrix_Classes::Indexmatrix& x_activity) const;

    /// return the current value of the dual variables to the constraints (after a successful QPIsolve the optimal solution)
    virtual const CH_Matrix_Classes::Matrix& QPget_y() const {
      return y;
    }

    /// return the current value of the dual variables to the lower bounds of the quadratic variables (after a successful QPIsolve the optimal solution)
    virtual const CH_Matrix_Classes::Matrix& QPget_zlb() const {
      return zlb;
    }

    /// return the current value of the dual variables to the upper bounds of the quadratic variables (after a successful QPIsolve the optimal solution)
    virtual const CH_Matrix_Classes::Matrix& QPget_zub() const {
      return zub;
    }

    /// return the current value of the slack variables of the constraints (after a successful QPIsolve the optimal solution)
    virtual const CH_Matrix_Classes::Matrix& QPget_s() const {
      return s;
    }

    /// return the current value of the slack variables of the constraints (after a successful QPIsolve the optimal solution) together with an activity estimate (1 ... active, 0 ... at bounds); in sout inactive components are set to the respective bound
    virtual int QPget_s(CH_Matrix_Classes::Matrix& sout,
      CH_Matrix_Classes::Indexmatrix& s_activity) const;

    /// return the current value of the dual variables to the lower bounds of the constraints (after a successful QPIsolve the optimal solution)
    virtual const CH_Matrix_Classes::Matrix& QPget_rhszlb() const {
      return rhszlb;
    }

    /// return the current value of the dual variables to the upper bounds of the constraints (after a successful QPIsolve the optimal solution)
    virtual const CH_Matrix_Classes::Matrix& QPget_rhszub() const {
      return rhszub;
    }

    /// return the second order cone data which models the quadratic part in case the paramaeters hold use_socqp==true
    virtual const SOCIPProxBlock& QPget_socqp() const {
      return socqp;
    }

    /// give access to the second order cone data which models the quadratic part in case the paramaeters hold use_socqp==true
    virtual SOCIPProxBlock& QPset_socqp() {
      return socqp;
    }

    /// return the current value of the barrier parameter mu 
    virtual CH_Matrix_Classes::Real QPget_mu() const {
      return mu;
    }

    /// return the number of iterations where predictor promises large progress but the barrier parameter is reduced only by a little
    virtual CH_Matrix_Classes::Integer QPget_large_predictor_cnt() const {
      return large_predictor_cnt;
    }

    /// return the number of iterations until termination 
    virtual CH_Matrix_Classes::Integer QPget_iter() const {
      return iter;
    }
  };


  //@}

}

#endif

