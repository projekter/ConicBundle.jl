/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/QPSolverParameters.hxx
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


#ifndef CONICBUNDLE_QPSOLVERPARAMETERS_HXX
#define CONICBUNDLE_QPSOLVERPARAMETERS_HXX

/**  @file QPSolverParameters.hxx
    @brief Header declaring the classes ConicBundle::QPSolverParameters
    @version 1.0
    @date 2020-03-18
    @author Christoph Helmberg
*/



#include "QPKKTSolverObject.hxx"

namespace ConicBundle {


  /** @ingroup ConstrainedQPSolver
   */
   //@{

   /** @brief parameters for steering the termination criteria and solution method of the solver
    */

  class QPSolverParameters : public virtual QPSolverParametersObject, public virtual CBout {
  private:
    //--- preset termination parameters
    CH_Matrix_Classes::Real min_objective_relprec; ///< minimum relative precision in objective
    CH_Matrix_Classes::Integer maxiter; ///< maxium iteration number, negative means no bound

    //--- selecting system and system solver
    bool allow_unconstrained; ///< allow switching to the internal UQPSolver if applicable

    //--- adaptive termination parameters set by the code
    CH_Matrix_Classes::Real objective_gap_eps;   ///< relative precision for gap between lower and upper objective values
    CH_Matrix_Classes::Real primal_infeasibility_eps;  ///< absolute precision for primal feasibility (groundset)
    CH_Matrix_Classes::Real dual_infeasibility_eps; ///< absolute precision for dual feasibility 
    CH_Matrix_Classes::Real lower_bound;  ///< optimal dual objective value should exceed that
    CH_Matrix_Classes::Real upper_bound;  ///< optimal primal objective should not exceed that
    CH_Matrix_Classes::Real lower_bound_gap_eps; ///< relative size requirement of objective gap to primal value minus lower_bound
    CH_Matrix_Classes::Real upper_bound_gap_eps; ///< relative size requirement of objective gap to upper_bound - dual value 

    QPKKTSolverObject* KKTsolver; ///< provides the routine described in QPKKTSolverObject, see also \ref ConstrainedQPSolver

    bool use_predictor_corrector; ///< default true, set to false if just one solve per KKT system is desired
    bool use_neighborhood; ///< default false, set to true if the line search should ensure staying in the neighborhood
    CH_Matrix_Classes::Real nbh_ub; ///< curve searches try to stay inside the neighborhood to this value
    CH_Matrix_Classes::Real nbh_lb; ///< the barrier parameter reduction aims for this value and below this value a predictor step is allowed

    bool use_socqp; ///< default false, set to true if the quadratic part should be modelled via a second order cone

    /// blocked copy constructor 
    QPSolverParameters(const QPSolverParameters& /*params*/);

    /// blocked assignment operator
    QPSolverParameters& operator=(const QPSolverParameters&);

  public:
    /// default constructor
    QPSolverParameters(CBout* cb = 0, int incr = -1);    ///> defaults are set in QPSolver.cxx
    /// virtual destructor
    virtual ~QPSolverParameters() {
      delete KKTsolver;
    }

    /// get this variable value
    CH_Matrix_Classes::Real QPget_min_objective_relprec() const {
      return min_objective_relprec;
    }
    /// get this variable value
    CH_Matrix_Classes::Real QPget_objective_gap_eps() const {
      return objective_gap_eps;
    }
    /// get this variable value
    CH_Matrix_Classes::Real QPget_primal_infeasibility_eps() const {
      return primal_infeasibility_eps;
    }
    /// get this variable value
    CH_Matrix_Classes::Real QPget_dual_infeasibility_eps() const {
      return dual_infeasibility_eps;
    }
    /// get this variable value
    CH_Matrix_Classes::Real QPget_lower_bound() const {
      return lower_bound;
    }
    /// get this variable value
    CH_Matrix_Classes::Real QPget_upper_bound() const {
      return upper_bound;
    }
    /// get this variable value
    CH_Matrix_Classes::Real QPget_lower_bound_gap_eps() const {
      return lower_bound_gap_eps;
    }
    /// get this variable value
    CH_Matrix_Classes::Real QPget_upper_bound_gap_eps() const {
      return upper_bound_gap_eps;
    }
    /// get this variable value
    CH_Matrix_Classes::Integer QPget_maxiter() const {
      return maxiter;
    }
    ///  get this variable value
    QPKKTSolverObject* QPget_KKTsolver() {
      return KKTsolver;
    };
    /// get this variable value
    bool QPget_use_predictor_corrector() const {
      return use_predictor_corrector;
    }
    /// get this variable value
    bool QPget_use_neighborhood() const {
      return use_neighborhood;
    }
    /// get this variable value
    CH_Matrix_Classes::Real QPget_nbh_ub() const {
      return nbh_ub;
    }
    /// get this variable value
    CH_Matrix_Classes::Real QPget_nbh_lb() const {
      return nbh_lb;
    }
    /// get this variable value
    bool QPget_use_socqp() const {
      return use_socqp;
    }


    /// set this variable value
    int QPset_min_objective_relprec(CH_Matrix_Classes::Real eps) {
      min_objective_relprec = (eps > 0) ? eps : min_objective_relprec; return (eps <= 0);
    }
    /// set this variable value
    int QPset_objective_gap_eps(CH_Matrix_Classes::Real eps) {
      objective_gap_eps = (eps > 0) ? eps : objective_gap_eps; return (eps <= 0);
    }
    /// set this variable value
    int QPset_primal_infeasibility_eps(CH_Matrix_Classes::Real eps) {
      primal_infeasibility_eps = (eps > 0) ? eps : primal_infeasibility_eps; return (eps <= 0);
    }
    /// set this variable value
    int QPset_dual_infeasibility_eps(CH_Matrix_Classes::Real eps) {
      dual_infeasibility_eps = (eps > 0) ? eps : dual_infeasibility_eps; return (eps <= 0);
    }
    /// set this variable value
    int QPset_lower_and_upper_bounds(CH_Matrix_Classes::Real lb, CH_Matrix_Classes::Real ub) {
      if (lb > ub) return 1; lower_bound = lb; upper_bound = ub; return 0;
    }
    /// set this variable value
    int QPset_lower_bound_gap_eps(CH_Matrix_Classes::Real eps) {
      lower_bound_gap_eps = (eps > 0) ? eps : lower_bound_gap_eps; return (eps <= 0);
    }
    /// set this variable value
    int QPset_upper_bound_gap_eps(CH_Matrix_Classes::Real eps) {
      upper_bound_gap_eps = (eps > 0) ? eps : upper_bound_gap_eps; return (eps <= 0);
    }
    /// set this variable value
    int QPset_maxiter(CH_Matrix_Classes::Integer mi) {
      maxiter = mi; return 0;
    }

    /// delete previous solver and replace by the new one (should not be zero when calling the solver)
    int QPset_KKTsolver(QPKKTSolverObject* in_KKTsolver) {
      delete KKTsolver; KKTsolver = in_KKTsolver; return 0;
    }

    /// if set to true (=default), a predictor corrector approach is used for solving the KKT system (solve for barrier parameter mu=0, guess mu, solve again for this mu including some bilinear perturbation) otherwise the barrier parameter is set apriori and the step computed in one solve 
    int QPset_use_predictor_corrector(bool upc) {
      use_predictor_corrector = upc; return 0;
    }

    /// if set to true (default: false), the lines search is carried out with respect to the neighborhood polynomial ensuring the afte this step the point is again inside the neighborhood of the central path
    int QPset_use_neighborhood(bool nbh) {
      use_neighborhood = nbh; return 0;
    }

    /// set the upper bound on the neighborhood that should be ensured in curve searches; ensures eps_Real<=nbhlb<=nbhub (nbhub should be < 1. and <=.35 is safe)
    int QPset_nbh_bounds(CH_Matrix_Classes::Real nbhlb, CH_Matrix_Classes::Real nbhub) {
      nbh_lb = CH_Matrix_Classes::max(CH_Matrix_Classes::eps_Real, nbhlb);
      nbh_ub = CH_Matrix_Classes::max(nbh_lb, nbhub);
      return 0;
    }

    /// if set to true (default: false), the quadratic term is modelled via a second order cone approach
    int QPset_use_socqp(bool s) {
      use_socqp = s; return 0;
    }

    /// set to true/false if switching to the unconstrained solver is allowed or not
    int QPset_allow_UQPSolver(bool allow) {
      allow_unconstrained = allow; return 0;
    }
    /// set to true/false if switching to the unconstrained solver is allowed or not
    bool QPallow_UQPSolver() {
      return allow_unconstrained;
    }
  };





  //@}

}

#endif

