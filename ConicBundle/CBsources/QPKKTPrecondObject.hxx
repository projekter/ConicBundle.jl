/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/QPKKTPrecondObject.hxx
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



#ifndef CONICBUNDLE_QPKKTPRECONDOBJECT_HXX
#define CONICBUNDLE_QPKKTPRECONDOBJECT_HXX

/**  @file QPKKTPrecondObject.hxx
    @brief Header declaring the abstract class ConicBundle::QPKKTPrecondObject
    @version 1.0
    @date 2020-03-20
    @author Christoph Helmberg

*/

#include "QPSolverObject.hxx"
#include "QPModelBlockObject.hxx"
#include "clock.hxx"

namespace ConicBundle {

class QPSolverParameters;

  
/** @ingroup ConstrainedQPSolver
 */
  //@{

  /** @brief Abstract Interface for preconditioners to be used with a
      QPIterativeKKTSolver and a CH_Matrix_Classes::IterativeSolverObject,
      see \ref IterativeSolverInterfaces.  It will depend on the
      system setup and the solver method which preconditioning
      routines are called and what requirements the preconditioners
      have to fulfill. Feasible combinations lie in the responsibility
      of the caller and are not checked for correctness.

      As describe for QPIterativeKKTSolver the main structure is
      as follows:
 
      - init_data() is called from QPIterativeKKTSolver::QPinit_KKTdata()
        once at the beginning of a new QP problem to communicate the
	basic problem data

      - init_system() is called from QPIterativeKKTSolver::QPinit_KKTsystem()
        once at each iteration of the interior point solver when setting
        when setting up the system before solving. Here the preconditioners
	are precomputed so that their use in precondM1() and precondM2() is
	efficient

      - precondM1() and precondM2() are called by the actual iterative solver
        CH_Matrix_Classes::IterativeSolverObject 
  */
  
class QPKKTPrecondObject: public virtual CBout
{
protected:

  //--- data describing the KKT system
  QPSolverProxObject* Hp;  ///< points to the quadratic cost representation, may NOT be NULL afer init
  QPModelBlockObject* model; ///< points to the cutting model information, may be NULL
  const CH_Matrix_Classes::Sparsemat* A;  ///< points to a possibly present constraint matrix, may be NULL
  const CH_Matrix_Classes::Indexmatrix* eq_indices; ///< if not NULL, these rows of A correspond to equations; needed for checking applicability of this Object
  bool SchurComplAineq; ///< if true, the inequalities of A are Schur complemented into the H block

  CH_Matrix_Classes::Real Hfactor; ///< the prox term of Hp is multiplied by this
  
public:
  /// reset data to empty
  virtual void clear()
  {Hp=0;model=0;A=0;eq_indices=0;}
  
  /// default constructor
  QPKKTPrecondObject(CBout* cb=0,int cbinc=-1):
    CBout(cb,cbinc),Hp(0),model(0),A(0),eq_indices(0),SchurComplAineq(false)
    {}

  /// virtual destructor
  virtual ~QPKKTPrecondObject();
  
  /// returns 1 if this class is not applicable in the current data situation, otherwise it stores the data pointers and these need to stay valid throught the use of the other routines but are not deleted here
  virtual int init_data(QPSolverProxObject* Hp, ///< may not be be NULL 
			QPModelBlockObject* model, ///< may be NULL
			const CH_Matrix_Classes::Sparsemat* A, ///< may be NULL
			const CH_Matrix_Classes::Indexmatrix* eq_indices, ///< if not NULL these rows of A correspond to equations
			bool SchurComplAineq ///< if true, the inequalities of A are Schur complemented into the H block
			)=0;
  
  /// set up the primal dual KKT system for being solved for predictor and corrector rhs; the input objects KKTdiagx, KKTdiagy and Hfactor will not change during use of the preconditioner, so it suffices to store the address if they are need during application of the preconditioner 
  virtual int init_system(const CH_Matrix_Classes::Matrix& KKTdiagx,
			  const CH_Matrix_Classes::Matrix& KKTdiagy,
			  CH_Matrix_Classes::Real Hfactor,
			  CH_Matrix_Classes::Real prec,
			  QPSolverParameters* params) =0;

  ///return (an estimate of) the minimum eigenvalue of the preconditioner M1^{-1}; this is used, e.g., to correct the precission in MINRES
  virtual CH_Matrix_Classes::Real get_lmin_invM1() {return 1.;}
  
  ///returns M1^{-1}*vec; default: M1=I
  virtual int precondM1(CH_Matrix_Classes::Matrix& /* vec */) {return 0;}
  
  ///returns M2^{-1}vec; default: M2=I
  virtual int precondM2(CH_Matrix_Classes::Matrix& /* vec */) {return 0;}

  ///for estimating the condition number with M1=G*G^T this returns G^{-1}*vec; default: G=I
  virtual int precond_invG1(CH_Matrix_Classes::Matrix& /* vec */) {return 0;}

  ///for estimating the condition number with M1=G*G^T this returns G^{-T}*vec; default: G=I
  virtual int precond_invG1tran(CH_Matrix_Classes::Matrix& /* vec */) {return 0;}

  ///for estimating the condition number directly for the preconditioned part only; negative numbers indicate that the routine is not implemented
  virtual CH_Matrix_Classes::Integer precond_size() {return -1;}
  
  ///for estimating the condition number directly for the preconditioned part only
  virtual int cond_number_mult(CH_Matrix_Classes::Matrix& /* vec */,
			       const CH_Matrix_Classes::Matrix& /* KKTdiagx */,
			       const CH_Matrix_Classes::Matrix& /* KKTdiagy */)
  {return 1;}

  ///for estimating the condition number directly for the preconditioned part only
  virtual  CH_Matrix_Classes::Real get_condition_number(const CH_Matrix_Classes::Matrix& KKTdiagx,
							const CH_Matrix_Classes::Matrix& KKTdiagy);
  
  /// for evaluation purposes with iterative solvers, return the rank of the precondiontioner used (or the number of n-vector multiplications per call)
  virtual CH_Matrix_Classes::Integer get_precond_rank()
  {return -1;}

  /// for evaluation purposes with iterative solvers, return the time spent in the multiplication with the preconditioner
  virtual CH_Tools::Microseconds get_t_precond_mult()
  {return 0;}

  /// for evaluation purposes with iterative solvers, reset the time spent in the multiplication with the preconditioner to zero
  virtual void reset_t_precond_mult()
  {}
};




  //@}
  
}

#endif

