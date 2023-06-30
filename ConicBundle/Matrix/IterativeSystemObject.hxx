/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  Matrix/IterativeSystemObject.hxx
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



#ifndef CH_MATRIX_CLASSES__ITERATIVESYSTEMOBJECT_HXX
#define CH_MATRIX_CLASSES__ITERATIVESYSTEMOBJECT_HXX

/**  @file IterativeSystemObject.hxx
    @brief Header declaring the abstract classes CH_Matrix_Classes::IterativeSystemObject and CH_Matrix_Classes::IterativeSolverObject 
    @version 1.0
    @date 2011-12-29
    @author Christoph Helmberg

*/

#ifndef CH_MATRIX_CLASSES__MATRIX_HXX
#include "matrix.hxx"
#endif

namespace CH_Matrix_Classes {
  
/**@defgroup IterativeSolverInterfaces Interfaces and Classes for Iterative Solvers like PCG, MINRES and PSQMR 

   @brief Uniform interface for routines for solving positive definite
   or indefinite symmetric or unsymmetric systems Ax=b by iterative
   methods with or without preconditioning. The meaning, the
   requirements and the use of the preconditioner routines depends on
   the methods.
*/
  //@{

  /** @brief Abstract base class for supplying the system for an iterative solver
 
      It will depend on the solution method which preconditioning
      routines are called and what requirements the system and the
      preconditioners have to fulfill.
  */
  
  class IterativeSystemObject
  {
  public:
    virtual ~IterativeSystemObject();

    ///returns the right hand side vector (dense)
    virtual const Matrix& ItSys_rhs() =0;
  
    ///returns out_vec=(system matrix)*in_vec
    virtual int ItSys_mult(const Matrix& in_vec,Matrix& out_vec) =0;
  
    ///returns M1^{-1}*vec; default: M1=I
    virtual int ItSys_precondM1(Matrix& /* vec */) {return 0;}

    ///returns M2^{-1}vec; default: M2=I
    virtual int ItSys_precondM2(Matrix& /* vec */) {return 0;}

  };


  /** @brief Abstract interface to iterative methods for solving Ax=b given by an IterativeSystemObject

  */

  class IterativeSolverObject
  {
  public:
    ///
    virtual ~IterativeSolverObject();
    
   /** @name Set and Get Parameters
       
       There should be no need to set any parameters, default values should be
       available and reasonable.
   */
    //@{
  
    /// set maximum number of iterations 
    virtual void set_maxit(Integer in_maxit)=0;
    /// get maximum number of iterations 
    virtual Integer get_maxit() const =0;

    /// returns the error code of the last call
    virtual int get_err() const =0;
    /// returns the number of matrix-vector multiplications of the last call
    virtual Integer get_nmult() const =0;
    /// returns the residual norm of last call
    virtual Real get_residual_norm() const =0;
    /// returns the average of the achieved reduction factor per iteration 
    virtual Real get_avg_reduction() const =0;
    /// returns the (absolute) precision requirement for termination used in the last call 
    virtual Real get_termprec() const =0;

    //@}
    
    /// compute the solution for system into x with (absolute) residual precision termprec
    virtual int compute(IterativeSystemObject& system, ///< the system information with precond
                        Matrix& x,        ///< on input: starting point (x.dim()==0 uses 0-vector), on output: approx. sol.,
			Real termprec,  ///< !absolute! termination precision, stop if residual norm<= termprec,
			Matrix* storex=0, ///< if not null, store initial x and the x of i*recordstep 
			Integer storestep=0 ///< if 0 and xrecord!=0 store only the initial x
			)=0;

   /** @name Input/Output
       
   */
    //@{
  
    /// set output stream and level of detail of log output (for debugging) 
    virtual void set_out(std::ostream* out=0,int print_level=1)=0;

    //@}
  };


  //@}
  
}

#endif

