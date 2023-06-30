/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/QPKKTSolverObject.hxx
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


#ifndef CONICBUNDLE_QPKKTSOLVEROBJECT_HXX
#define CONICBUNDLE_QPKKTSOLVEROBJECT_HXX

/**  @file QPKKTSolverObject.hxx
    @brief Header declaring the class ConicBundle::QPKKTSolverObject
    @version 1.0
    @date 2020-03-18
    @author Christoph Helmberg
*/

#include "QPSolverObject.hxx"
#include "QPModelBlockObject.hxx"

namespace ConicBundle {

  class QPSolverParameters;

  /** @ingroup ConstrainedQPSolver
   */
   //@{

   /** @brief abstract class for setting up and solving the primal dual KKT System within QPSolverBasicStructures

       For initializing such a solver, first call QPinit_KKTdata(), whose
       arguments must give access to the data objects describing the basic
       problem data of all KKT Systems in this call to the QP solver.
       QPinit_KKTdata() is only called once in the beginning of solving a
       QP probem and gives the following data:

       - \f$H\f$ describes the prox term given by the QPSolverProxObject,

       - \f$A\f$ is the matrix giving the linear constraints of the groundset

       - \f$B\f$ describe the bundle with subgradient information given as
         rows in the QPModelBlockObject,

       - \f$C\f$ describes the collected trace/penalty constraints of the models
         in QPModelBlockObject.

       For forming each particular primal dual KKT system to be used in
       one predictor and mostly followed by one corrector solve, the
       routine QPinit_KKTsystem() and the model will provide positive
       (definite) barrier contributions. QPinit_KKTsystem() will be
       called at the beginning of each interior point iteration (before
       computing predictor and corrector step). The barrier contributions
       are represented here by

       - \f$D_x\f$ (box constraints on variables, KKTdiagx in QPinit_KKTsystem()),

       - \f$D_A\f$ (maybe box constraints on linear constraints \f$A\f$ unless 0 for
         equations, KKTdiagy in QPinit_KKTsystem()),

       - \f$D_B\f$ (model variables obtained implicitly from the QPModelBlockObject),

       - \f$D_C\f$ (maybe slack variables on trace constraints of the
         model unless 0 for equations, obtained implicitly from the
         QPModelBlockObject).

       Finally one predictor and one corrector call to QPsolve_KKTsystem() per
       interior point iteration asks to solve for the right hand side values

       - \f$ r_H\f$ "dual" right hand side (dualrhs in QPsolve_KKTsystem())

       - \f$r_A\f$ "primal" right hand side of the ground set (primalrhs in QPsolve_KKTsystem())

       - \f$r_B\f$ "bundle" right hand side of the model (provided by QPModelBlockObject)

       - \f$r_C\f$  and a model constraint right hand side (provided by QPModelBlockObject).

       With these the primal dual KKT System of a QPKKTSolverObject reads

       \f[
        \left[\begin{array}{cccc}
        H + D_x & A^\top & B^\top & 0 \\
        A & -D_A & 0 & 0 \\
        B & 0 & -D_B & C \\
        0 & 0 & C & D_C
        \end{array} \right]
        \left(\begin{array}{c}x\\y\\ \xi \\ \zeta\end{array}\right)
        =
        \left(\begin{array}{c}r_H\\r_A\\ r_B \\r_C\end{array}\right)
       \f]

       Correspondingly we will speak of the blocks H, A, B, C of the system.

       Internally the QPKKTSolverObject may of course form variants via partial
       Schur complements and the like in order to best match the structural
       requirements of the model or application at hand.

       One choice involves modelling the quadratic term H internally by
       a second order cone. Eliminating this part just results in using
       a scaled version of H. This is the @a Hfactor used in QPinit_KKTsystem().
       Without the second order cone has value 1.

       Some implemented examples are indicated in \ref ConstrainedQPSolver.


    */

  class QPKKTSolverObject : public virtual CBout {
  protected:

    //--- data describing the KKT system
    QPSolverProxObject* Hp;  ///< points to the quadratic cost representation, may NOT be NULL afer init
    QPModelBlockObject* model; ///< points to the cutting model information, may be NULL
    const CH_Matrix_Classes::Sparsemat* A;  ///< points to a possibly present constraint matrix, may be NULL
    const CH_Matrix_Classes::Indexmatrix* eq_indices; ///< if not NULL, these rows of A correspond to equations; needed for checking applicability of this Object

    CH_Matrix_Classes::Real blockH_norm;  ///< for judging violation: (an estimate of) the norm of the H-row in the system
    CH_Matrix_Classes::Real blockA_norm;  ///< for judging violation: (an estimate of) the norm of the A-row in the system

  public:
    /// reset data to empty
    virtual void clear() {
      Hp = 0; model = 0; A = 0; eq_indices = 0; blockH_norm = 1.; blockA_norm = 1.;
    }

    /// default constructor
    QPKKTSolverObject(CBout* cb = 0, int cbinc = -1) :
      CBout(cb, cbinc), Hp(0), model(0), A(0), eq_indices(0), blockH_norm(1.), blockA_norm(1.) {
    }

    /// virtual destructor
    virtual ~QPKKTSolverObject();


    /// returns 1 if this class is not applicable in the current data situation, otherwise it stores the data pointers and these need to stay valid throught the use of the other routines but are not deleted here
    virtual int QPinit_KKTdata(QPSolverProxObject* Hp, ///< may not be be NULL 
      QPModelBlockObject* model, ///< may be NULL
      const CH_Matrix_Classes::Sparsemat* A, ///< may be NULL
      const CH_Matrix_Classes::Indexmatrix* eq_indices ///< if not NULL these rows of A correspond to equations
    ) = 0;

    /// set up the primal dual KKT system for being solved for predictor and corrector rhs in QPsolve_KKTsystem
    virtual int QPinit_KKTsystem(const CH_Matrix_Classes::Matrix& KKTdiagx,
      const CH_Matrix_Classes::Matrix& KKTdiagy,
      CH_Matrix_Classes::Real Hfactor,
      CH_Matrix_Classes::Real prec,
      QPSolverParameters* params) = 0;

    /// solve the KKTsystem to precision prec for the given right hand sides that have been computed for the value rhsmu of the barrier parameter and in which a rhscorr fraction (out of [0,1] of the corrector term have been included; in iterative solvers solx and soly may be used as starting points
    virtual int QPsolve_KKTsystem(CH_Matrix_Classes::Matrix& solx,
      CH_Matrix_Classes::Matrix& soly,
      const CH_Matrix_Classes::Matrix& primalrhs,
      const CH_Matrix_Classes::Matrix& dualrhs,
      CH_Matrix_Classes::Real rhsmu,
      CH_Matrix_Classes::Real rhscorr,
      CH_Matrix_Classes::Real prec,
      QPSolverParameters* params) = 0;

    /// for judging violation this returns (an estimate of) the norm of the H-row in the latest system
    virtual CH_Matrix_Classes::Real QPget_blockH_norm() {
      return blockH_norm;
    }

    /// for judging violation this returns (an estimate of) the norm of the A-row in the latest system
    virtual CH_Matrix_Classes::Real QPget_blockA_norm() {
      return blockA_norm;
    }

    /// for evaluation purposes with iterative solvers, return the number of matrix vector multiplications since the latest call to QPinit_KKTsystem (returns 0 for direct solvers)
    virtual CH_Matrix_Classes::Integer QPget_nmult() const {
      return 0;
    }

    /// for evaluation purposes with iterative solvers, return an estimate of the condition number (returns -1. for direct solvers or if not supplied)
    virtual CH_Matrix_Classes::Real QPget_condition_number() {
      return -1.;
    }

    /// for evaluation purposes with iterative solvers, return the rank of the precondiontioner used (or the number of n-vector multiplications per call)
    virtual CH_Matrix_Classes::Integer QPget_precond_rank() {
      return -1;
    }

    /// for evaluation purposes with iterative solvers, return the size of the system matrix
    virtual CH_Matrix_Classes::Integer QPget_system_size() {
      return 0;
    }
  };





  //@}

}

#endif

