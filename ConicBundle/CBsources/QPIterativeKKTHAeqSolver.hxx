/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/QPIterativeKKTHAeqSolver.hxx
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


#ifndef CONICBUNDLE_QPITERATIVEKKTHAEQSOLVER_HXX
#define CONICBUNDLE_QPITERATIVEKKTHAEQSOLVER_HXX

/**  @file QPIterativeKKTSolver.hxx
    @brief Header declaring the class ConicBundle::QPIterativeKKTHASolver
    @version 1.0
    @date 2020-03-20
    @author Christoph Helmberg
*/



#include "QPIterativeKKTSolver.hxx"

namespace ConicBundle {

/** @ingroup ConstrainedQPSolver 
 */
//@{

/** @brief Iterative solver for the reduced symmetric, in general
    indefinite primal dual KKT System within QPSolverBasicStructures,
    where only the block H and the equality rows of A remain, 
    the inequalities of A as well as B and C are removed by Schur
    complements. If there are no equalities in A, PCG may used, but MinRes
    seems to be more stable.

    See the text to QPKKTSolverObject for the terminology of the
    primal dual KKT System and the general outline. 

    For this class the preconditioner offered by QPKKTSubspaceHAeqPrecond 
    may be worth to try.
    
    The most important routines of the model described in the QPModelBlockObject 
    that are required here are
    (besides sizes and multiplications with the bundle matrix B)

    - QPModelBlockObject::add_Schur_rhs() for adding the right hand side contribution to the H block as induced by the Schur complement

    - QPModelBlockObject::compute_step() for computing the model solution step given the solution for the H block

    - QPModelBlockObject::add_Schur_mult() for adding the multiplication with the Schur complemented B and C blocks 
    
    

 */

  class QPIterativeKKTHAeqSolver: public QPIterativeKKTSolver
{
protected:

  //--- data describing the KKT system
  //QPSolverProxObject* Hp;  ///< points to the quadratic cost representation, may NOT be NULL afer init
  //QPModelBlockObject* model; ///< points to the cutting model information, may be NULL
  //const CH_Matrix_Classes::Sparsemat* A;  ///< points to a possibly present constraint matrix, may be NULL
  //const Indexmatrix* eq_indices; ///< if not NULL, these rows of A correspond to equations; needed for checking applicability of this Object

  // CH_Matrix_Classes::IterativeSolverObject* solver; ///< the iterative solution method to be used, may not be NULL, deleted on destruction or replacemnt
  // QPKKTPrecondObject* precond; ///< a preconditioning routine compatible with this system AND the solver, may be NULL (no preconditioning); consistency cannot be checked here and must be guaranteed externally; deleted on destruction or replacement

  // CH_Matrix_Classes::Matrix KKTdiagx;  ///< diagonal term to be added to the H-block due to bounds on the qudratic variables
  // CH_Matrix_Classes::Matrix KKTdiagy;  ///< diagonal term to be added to the A-block due to bounds on the constraints
  
  // CH_Matrix_Classes::Matrix sysrhs;  ///< the right hand side of the system
  // CH_Matrix_Classes::Matrix sol;     ///< solution vector of the entire system
  // CH_Matrix_Classes::Matrix solmod;  ///< temporary solution part of the model variables 
  // CH_Matrix_Classes::Matrix solcstr; ///< temporary solution part of the model constraint variables

  // CH_Matrix_Classes::Matrix in_vecx;  ///< temporary storage for x part of in_vec in ItSys_mult
  // CH_Matrix_Classes::Matrix in_vecy;  ///< temporary storage for y part of in_vec in ItSys_mult
  // CH_Matrix_Classes::Matrix in_model; ///< temporary storage for model part of in_vec in ItSys_mult

  CH_Matrix_Classes::Matrix xrecord;      ///< testing phase attempt to use iterates of the iterative solver for preconditiong, stored here
  CH_Matrix_Classes::Integer recordstep;  ///< each recordstep iteration the iterate is stored in xrecord

  CH_Tools::Clock clock; ///< for taking the time spent in ItSys_mult
  CH_Tools::Microseconds t_itsys_mult; ///< time spent in ItSys_mult
  
public:
  // reset data to empty but keep solver and preconditioner
  //virtual void clear();

  /// default constructor
  QPIterativeKKTHAeqSolver(CH_Matrix_Classes::IterativeSolverObject* insolver,QPKKTPrecondObject* inprecond=0,CBout* cb=0,int cbinc=-1):
    CBout(cb,cbinc),QPIterativeKKTSolver(insolver,inprecond,cb,cbinc)
  {}

  /// virtual destructor
  virtual ~QPIterativeKKTHAeqSolver();
  
  
  /// returns 1 if this class is not applicable in the current data situation, otherwise it stores the data pointers and these need to stay valid throught the use of the other routines but are not deleted here
  virtual int QPinit_KKTdata(QPSolverProxObject* Hp, ///< may not be be NULL 
  			     QPModelBlockObject* model, ///< may be NULL
  			     const CH_Matrix_Classes::Sparsemat* A, ///< may be NULL
  			     const CH_Matrix_Classes::Indexmatrix* eq_indices ///< if not NULL these rows of A correspond to equations
  			     );

  // set up the primal dual KKT system for being solved for predictor and corrector rhs in QPsolve_KKTsystem
  // virtual int QPinit_KKTsystem(const CH_Matrix_Classes::Matrix& KKTdiagx,
  // 			       const CH_Matrix_Classes::Matrix& KKTdiagy,
  // 			       CH_Matrix_Classes::Real prec,
  // 			       QPSolverParameters* params);

  /// solve the KKTsystem to precision prec for the given right hand sides that have been computed for the value rhsmu of the barrier parameter and in which a rhscorr fraction (out of [0,1] of the corrector term have been included; in iterative solvers solx and soly may be used as starting points
  virtual int QPsolve_KKTsystem(CH_Matrix_Classes::Matrix& solx,
				 CH_Matrix_Classes::Matrix& soly,
				 const CH_Matrix_Classes::Matrix& primalrhs,
				 const CH_Matrix_Classes::Matrix& dualrhs,
				 CH_Matrix_Classes::Real rhsmu,
				 CH_Matrix_Classes::Real rhscorr,
				 CH_Matrix_Classes::Real prec,
				 QPSolverParameters* params);

  //returns the right hand side vector (dense)
  //virtual const CH_Matrix_Classes::Matrix& ItSys_rhs() 
  //{return sysrhs;}
  
  ///returns out_vec=(system matrix)*in_vec
  virtual int ItSys_mult(const CH_Matrix_Classes::Matrix& in_vec,CH_Matrix_Classes::Matrix& out_vec);
  
  //returns M1^{-1}*vec; default: M1=I
  //virtual int ItSys_precondM1(CH_Matrix_Classes::Matrix& vec) 
  //{return (precond?precond->precondM1(vec):0);}
  
  //returns M2^{-1}vec; default: M2=I
  //virtual int ItSys_precondM2(CH_Matrix_Classes::Matrix& vec) 
  //{return (precond?precond->precondM2(vec):0);}

  /// for evaluation purposes with iterative solvers, return the size of the system matrix
  virtual CH_Matrix_Classes::Integer QPget_system_size() 
  {return KKTdiagx.rowdim()+(eq_indices?eq_indices->rowdim():0);}
  
  };





  //@}

}

#endif

