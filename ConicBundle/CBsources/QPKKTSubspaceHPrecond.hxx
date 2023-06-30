/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/QPKKTSubspaceHPrecond.hxx
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



#ifndef CONICBUNDLE_QPKKTSUBSPACEHPRECOND_HXX
#define CONICBUNDLE_QPKKTSUBSPACEHPRECOND_HXX

/**  @file QPKKTSubspaceHPrecond.hxx
    @brief Header declaring the abstract class ConicBundle::QPKKTSubspaceHPrecond
    @version 1.0
    @date 2020-03-20
    @author Christoph Helmberg

*/

#include "QPKKTPrecondObject.hxx"

namespace ConicBundle {

  class QPSolverParameters;


  /** @ingroup ConstrainedQPSolver
   */
   //@{

   /** @brief Subspace projection preconditioner for the H-block of the
       KKT-System assuming that B and C have been Schur complemented
       into the H-block. For the A-Block, if it is there, the identity
       is used. If there is no A-block, PCG may be used.

       See the text to QPKKTSolverObject for the terminology of the
       primal dual KKT System and the general outline.

       This implementation of a QPKKTPrecondObject is suitable
       e.g. for QPIterativeKKTHASolver and QPIterativeKKTHAeqSolver.

       Depending on @a method the subspace is chosen deterministically
       or randomly generated. The code is still rather experimental and
       therefore currently a multitude of variants are implemented for
       testing.

       The current variants of @a method are (recommendation: use 31)

       -  0 ... no preconditioning

       -  1 ... diagonal preconditioning

       - 10 ... exact preconditioning (full low rank subspace, way to expensive in general)

       - 11 ... Johnson-Lindenstrauss preconditioning, Achlioptas {-1,0,1} version

       - 12 ... JL preconditioning, normal distribution N(0,1/cols)

       - 20 ... randomized orthogonal subspace, updated from one system to the next

       - 21 ... randomized orthogonal subspace, regenerated for every system

       - 3[0-9] ... deterministic subspace selection, where the i in
          "3i" uses the minmum value 10^i for including a column.
          31 is the recommended version.

       The most important routines of the model described in the
       QPModelBlockObject that are required here are (besides sizes and
       multiplications with the bundle matrix B)

       - QPModelBlockObject::prepare_BCSchur_JLprecond() for multiplying the Schur complemented model part with a subspace

       - QPModelBlockObject::propose_BCSchur_pcsubspace() for deterministically appending subspace generated columns of the Schur complemented model part

   */

  class QPKKTSubspaceHPrecond : public QPKKTPrecondObject {
  protected:

    //--- data describing the KKT system
    // QPSolverProxObject* Hp;  ///< points to the quadratic cost representation, may NOT be NULL afer init
    // QPModelBlockObject* model; ///< points to the cutting model information, may be NULL
    // const CH_Matrix_Classes::Sparsemat* A;  ///< points to a possibly present constraint matrix, may be NULL
    // const CH_Matrix_Classes::Indexmatrix* eq_indices; ///< if not NULL, these rows of A correspond to equations; needed for checking applicability of this Object

    CH_Matrix_Classes::Integer method;     ///< selects generation method for preconditioner

    CH_Matrix_Classes::Matrix diagH;       ///< diagonal of the prox term
    const CH_Matrix_Classes::Matrix* Vp;   ///< lowrank part of the prox term

    CH_Matrix_Classes::Integer last_nmult;  ///< last number of multplications in iterative solver
    CH_Matrix_Classes::Real  max_sigma; ///< if >0 it gives the last maximum singular value found

    CH_Matrix_Classes::Matrix Diag_inv;        ///< diagonal with KKT diagonal part, inverted, possibly the sqrt of it (if there is a low rank part)
    CH_Matrix_Classes::Real diaginvval;        ///< if the value is > 0 then Diag_inv must be this values times the all ones vector

    CH_Matrix_Classes::Matrix subspace;        ///< the subspace used for projection
    CH_Matrix_Classes::Matrix lowrank;         ///< the selected lowrank representation
    CH_Matrix_Classes::Indexmatrix pivlowrank; ///< the pivot sequence used
    CH_Matrix_Classes::Matrix eigvals;         ///< the eigenvalues of the low rank approximation on the subspace
    CH_Matrix_Classes::Matrix eigvecs;         ///< the eigenvectors within the subspace

    CH_Matrix_Classes::Matrix tmpmat;          ///< temporary matrix
    CH_Matrix_Classes::Matrix tmpvec;          ///< temporary matrix
    CH_Matrix_Classes::Symmatrix tmpsym;       ///< temporary symmetric matrix
    CH_Matrix_Classes::Matrix keepeigs;        ///< eigenvalues used to update the subspace
    CH_Matrix_Classes::Matrix keepvecs;        ///< eignevectors used to update the subspace
    CH_Matrix_Classes::Matrix rotmat;          ///< rotation to update the subspace

    CH_Matrix_Classes::Matrix Q; ///< [only used in testing] QR factorization of D^(-.5)*lowrank*eigvecs

    CH_Tools::Clock clock; ///< for taking the time spent in various parts
    CH_Tools::Microseconds t_gen_subspace;  ///< time spent in generating the subspace
    CH_Tools::Microseconds t_comp_lowrank;  ///< time spent in generating the lowrank matrix by multiplying with subspace
    CH_Tools::Microseconds t_comp_svd;   ///< time spent in computing the singular value decomposition
    CH_Tools::Microseconds t_precond_mult;   ///< time spent in multiplying with the preconditioner


  public:
    // reset data to empty
    //virtual void clear()
    //{Hp=0;model=0;A=0;eq_indices=0;}

    /// default constructor
    QPKKTSubspaceHPrecond(CH_Matrix_Classes::Integer inmethod = 0, CBout* cb = 0, int cbinc = -1);

    /// virtual destructor
    virtual ~QPKKTSubspaceHPrecond();

    /// returns 1 if this class is not applicable in the current data situation, otherwise it stores the data pointers and these need to stay valid throught the use of the other routines but are not deleted here
    virtual int init_data(QPSolverProxObject* Hp, ///< may not be be NULL 
      QPModelBlockObject* model, ///< may be NULL
      const CH_Matrix_Classes::Sparsemat* A, ///< may be NULL
      const CH_Matrix_Classes::Indexmatrix* eq_indices, ///< if not NULL these rows of A correspond to equations
      bool SchurComplAineq ///< if true, the inequalities of A are Schur complemented into the H block
    );

    /// set up the primal dual KKT system for being solved for predictor and corrector rhs; the input objects KKTdiagx and KKTdiagy will not change during use of the preconditioner, so it suffices to store the address if they are need during application of the preconditioner 
    virtual int init_system(const CH_Matrix_Classes::Matrix& KKTdiagx,
      const CH_Matrix_Classes::Matrix& KKTdiagy,
      CH_Matrix_Classes::Real Hfactor,
      CH_Matrix_Classes::Real prec,
      QPSolverParameters* params);

    ///returns M1^{-1}*vec; default: M1=I
    virtual int precondM1(CH_Matrix_Classes::Matrix& vec);

    //returns M2^{-1}vec; default: M2=I
    //virtual int precondM2(CH_Matrix_Classes::Matrix& /* vec */) {return 0;}

    ///if the method admits this, let the subspace be chosen externally
    virtual int set_subspace(const CH_Matrix_Classes::Matrix& insubspace) {
      if (method / 10 == 3) subspace = insubspace; return 0;
    }

    ///return (an estimate of) the minimum eigenvalue of the preconditioner M1^{-1}; this is used, e.g., to correct the precission in MINRES
    virtual CH_Matrix_Classes::Real get_lmin_invM1();

    ///for estimating the condition number with M1=G*G^T this returns G^{-1}*vec; default: G=I
    virtual int precond_invG1(CH_Matrix_Classes::Matrix& vec);

    ///for estimating the condition number with M1=G*G^T this returns G^{-T}*vec; default: G=I
    virtual int precond_invG1tran(CH_Matrix_Classes::Matrix& vec);

    ///for estimating the condition number directly for the preconditioned part only; negative numbers indicate that the routine is not implemented
    virtual CH_Matrix_Classes::Integer precond_size() {
      return diagH.rowdim();
    }

    ///for estimating the condition number directly for the preconditioned part only
    virtual int cond_number_mult(CH_Matrix_Classes::Matrix& vec,
      const CH_Matrix_Classes::Matrix& KKTdiagx,
      const CH_Matrix_Classes::Matrix& KKTdiagy);
    /// for evaluation purposes with iterative solvers, return the rank of the precondiontioner used (or the number of n-vector multiplications per call)
    virtual CH_Matrix_Classes::Integer get_precond_rank() {
      return eigvals.rowdim();
    }

    /// for evaluation purposes with iterative solvers, return the time spent in the multiplication with
    virtual CH_Tools::Microseconds get_t_precond_mult() {
      return t_precond_mult;
    }

    /// for evaluation purposes with iterative solvers, reset the time spent in the multiplication with the preconditioner to zero
    virtual void reset_t_precond_mult() {
      t_precond_mult = 0;
    }
  };




  //@}

}

#endif

