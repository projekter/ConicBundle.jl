/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/QPDirectKKTSolver.hxx
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


#ifndef CONICBUNDLE_QPDIRECTKKTSOLVER_HXX
#define CONICBUNDLE_QPDIRECTKKTSOLVER_HXX

/**  @file QPDirectKKTSolver.hxx
    @brief Header declaring the class ConicBundle::QPDirectKKTSolver
    @version 1.0
    @date 2020-03-18
    @author Christoph Helmberg
*/



#include "QPKKTSolverObject.hxx"

namespace ConicBundle {

  /** @ingroup ConstrainedQPSolver
   */
   //@{

   /** @brief implements a direct KKT Solver variant of QPKKTSolverObject

       See the text to QPKKTSolverObject for the terminology of the
       primal dual KKT System and the general outline.

       For solving the system the current class first forms the Schur
       complement with respect to \f$H+D_x\f$; this is done in dependence
       of whether the prox term pointed to by Hp is of the form of
       Diagonal plus low rank (-> low rank inversion) or dense or other
       (dense Cholesky). The computed Schur complement matrix
       [A; B; 0]H^{-1}[A^T,B^T,0] is stored in Schur_complement before
       it is added to ABC. The code then continues in dependence of the
       value of @a factorize_ABC

       - if factorize_ABC==false, it next factorizes the AB blocks by Cholesky,
         then the final C block by Cholesky and proceeds in the other direction
         for solving the system. If in this process a Cholesky factorization
         fails, the code sets factorize_ABC=true and re

       - if factorize_ABC==true, it uses Aasen on the indefinite remaining system
         of blocks A, B, and C.

       The most important routines of the model described in the QPModelBlockObject
       that are required here are
       (besides sizes and multiplications with the bundle matrix B)

       - QPModelBlockObject::add_BDBt() in forming the low rank Schur complement of H

       - QPModelBlockObject::add_localsys() in adding the barrier diagonal of the
         B block to the the Schur_complement

       - QPModelBlockObject::add_localrhs() for constructing the right hand side

       - QPModelBlockObject::computed_step() for communicating the solution step to the model


    */

  class QPDirectKKTSolver : public QPKKTSolverObject {
  private:

    //--- data describing the KKT system
    //BundleProxObject* Hp;  ///< points to the quadratic cost representation, may NOT be NULL afer init
    //QPModelBlockObject* model; ///< points to the cutting model information, may be NULL
    //const CH_Matrix_Classes::Sparsemat* A;  ///< points to a possibly present constraint matrix, may be NULL
    //const Indexmatrix* eq_indices; ///< if not NULL, these rows of A correspond to equations; needed for checking applicability of this Object

    CH_Matrix_Classes::Integer dim;   ///< dimension of the quadratic term
    CH_Matrix_Classes::Integer Anr;   ///< number of rows in the constrain matrix =(A==0)?0:A->rowdim();
    CH_Matrix_Classes::Integer bsz;   ///< size of the bundle information in the cutting model
    CH_Matrix_Classes::Integer csz;   ///< number of constraints within the cutting model
    const CH_Matrix_Classes::Matrix* Vp; ///< used for low rank precondititioning

    //--- data for computing and storing precomputed parts of the system
    CH_Matrix_Classes::Real Hfactor; ///< scalar factor for H term
    CH_Matrix_Classes::Matrix Diag_inv; ///< inverse of the KKT diagonal in low rank versions
    CH_Matrix_Classes::Symmatrix Qchol; ///< Cholesky facotr of KKT Q or of its low rank part
    CH_Matrix_Classes::Symmatrix Schur_complement; ///< Schur complement of ABC with Q if no bounds
    CH_Matrix_Classes::Symmatrix AQiAt_inv; ///<  Cholesky or Aasen factor of ABC with Schur complement
    CH_Matrix_Classes::Symmatrix CABinvCt_inv; ///< for forming partial Schur complements
    CH_Matrix_Classes::Matrix LinvC; ///< for partial Schur complement w.r.t. C block
    CH_Matrix_Classes::Indexmatrix piv; ///< pivot sequence for Aasen

    //--- parameters for the solution method
    bool factorize_ABC; ///< if true, use Aasen to factor the ABC block

    /// computes the Schur complement if Hp->is_DLR()==TRUE, i.e., the quadratic term is representable as diagonal plus low rank
    int compute_DLR_Schur_complement(CH_Matrix_Classes::Symmatrix& sc);

    /// computes the Schur compelement if the quadratic term is only retrievable in dense form
    int compute_dense_Schur_complement(CH_Matrix_Classes::Symmatrix& sc);

  public:
    /// reset data to empty
    virtual void clear();

    /// default constructor
    QPDirectKKTSolver(bool in_factorize_ABC = false, CBout* cb = 0, int cbinc = -1) :
      QPKKTSolverObject(cb, cbinc), factorize_ABC(in_factorize_ABC) {
      clear();
    }

    /// virtual destructor
    virtual ~QPDirectKKTSolver();


    /// returns 1 if this class is not applicable in the current data situation, otherwise it stores the data pointers and these need to stay valid throught the use of the other routines but are not deleted here
    virtual int QPinit_KKTdata(QPSolverProxObject* Hp, ///< may not be be NULL 
      QPModelBlockObject* model, ///< may be NULL
      const CH_Matrix_Classes::Sparsemat* A, ///< may be NULL
      const CH_Matrix_Classes::Indexmatrix* eq_indices ///< if not NULL these rows of A correspond to equations
    );

    /// set up the primal dual KKT system for being solved for predictor and corrector rhs in QPsolve_KKTsystem
    virtual int QPinit_KKTsystem(const CH_Matrix_Classes::Matrix& KKTdiagx,
      const CH_Matrix_Classes::Matrix& KKTdiagy,
      CH_Matrix_Classes::Real Hfactor,
      CH_Matrix_Classes::Real prec,
      QPSolverParameters* params);

    /// solve the KKTsystem to precision prec for the given right hand sides that have been computed for the value rhsmu of the barrier parameter and in which a rhscorr fraction (out of [0,1] of the corrector term have been included; in iterative solvers solx and soly may be used as starting points
    virtual int QPsolve_KKTsystem(CH_Matrix_Classes::Matrix& solx,
      CH_Matrix_Classes::Matrix& soly,
      const CH_Matrix_Classes::Matrix& primalrhs,
      const CH_Matrix_Classes::Matrix& dualrhs,
      CH_Matrix_Classes::Real rhsmu,
      CH_Matrix_Classes::Real rhscorr,
      CH_Matrix_Classes::Real prec,
      QPSolverParameters* params);

  };

  //@}

}

#endif

