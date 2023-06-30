/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  Matrix/pcg.hxx
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



#ifndef CH_MATRIX_CLASSES__PCG_HXX
#define CH_MATRIX_CLASSES__PCG_HXX

/**  @file pcg.hxx
    @brief Header declaring CH_Matrix_Classes::PCG for preconditioned conjugate gradient for positive definite systems with positive definite preconditioner M1
    @version 1.0
    @date 2011-12-29
    @author Christoph Helmberg

*/

#include "IterativeSystemObject.hxx"

namespace CH_Matrix_Classes {

  /**@ingroup IterativeSolverInterfaces
  */
  //@{

  /** @brief Preconditioned Conjugate Gradient method for solving Ax=b with (symmetric) positive definite matrix A
      and positive definite preconditioner M1 where the
      preconditioned system is A'x'=b' with A'=M1^{-.5}AM1^{-.5},
      x'=M1^{.5}*x and b'=M1^{-.5}b.  System matrix and preconditioner are
      provided by a CH_Matrix_Classes::IterativeSystemObject

  */

  class PCG : public IterativeSolverObject {
  private:
    Integer maxit; ///< maximum number of matrix vector multiplications 
    Real resnorm;   ///< residual norm
    Real old_resnorm; ///< residual norm in the previous step
    Real avg_reduction; ///<average over the reduction factors
    Real termprec;     ///< absolute precision required in the last call
    Integer nmult;  ///< number of matrix vector multiplications/iterations
    Integer err;    ///< error code

    Matrix Ap;   ///< result of matrix multiplication
    Matrix r;    ///< residual
    Matrix z;    ///< temporary residual for step computation
    Matrix p;    ///< step

    Integer stall; ///< counts number of time reduction was less than .999


    std::ostream* myout;   ///< everything is output to *myout, default: 0 for no output
    int print_level;   ///<  level of iteration information that should be displayed

  public:
    /// default constructor
    PCG(std::ostream* out = 0, int pril = -1);
    ///
    ~PCG() {
    }

    /** @name Set and Get Parameters

        There should be no need to set any parameters, default values should be
        available and reasonable.
    */
    //@{

    /// set maximum number of iterations 
    void set_maxit(Integer in_maxit) {
      maxit = in_maxit;
    }
    /// get maximum number of iterations 
    Integer get_maxit() const {
      return maxit;
    }

    /// returns the error code of the last call
    int get_err() const {
      return err;
    }
    /// returns the number of matrix-vector multiplications of the last call
    Integer get_nmult() const {
      return nmult;
    }
    /// returns the residual norm of last call
    Real get_residual_norm() const {
      return resnorm;
    }
    /// returns the average of the achieved reduction factor per iteration 
    virtual Real get_avg_reduction() const {
      return avg_reduction;
    }
    /// return the (absolute) precision requirement for termination used in the last call 
    virtual Real get_termprec() const {
      return termprec;
    }

    //@}

    /// compute the solution for system into x with (absolute) residual precision termprec 
    int compute(IterativeSystemObject& system, ///< the system informatin with precond
      Matrix& x,        ///< on input: starting point (x.dim()==0 uses 0-vector), on output: approx. sol.,
      Real termprec,  ///< !absolute! termination precission, stop if residual norm<= termprec,
      Matrix* storex = 0, ///< if not null, store initial x and the x of i*recordstep 
      Integer storestep = 0 ///< if 0 and xrecord!=0 store only the initial x
    );

    /** @name Input/Output

    */
    //@{

    /// set output stream and level of detail of log output (for debugging) 
    virtual void set_out(std::ostream* in_out = 0, int in_print_level = 1) {
      myout = in_out; print_level = in_print_level;
    }

    //@}
  };

  //@}

}

#endif

