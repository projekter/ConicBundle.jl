/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/PSCAffineFunction.hxx
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



#ifndef CONICBUNDLE_PSCAFFINEFUNCTION_HXX
#define CONICBUNDLE_PSCAFFINEFUNCTION_HXX

/**  @file PSCAffineFunction.hxx
    @brief Header declaring the classes ConicBundle::PSCAffineFunction, ConicBundle::PSCAffineMinorantExtender (part of an implementation of ConicBundle::PSCModel)
    @version 1.0
    @date 2017-02-04
    @author Christoph Helmberg
*/

//------------------------------------------------------------

#include <map>
#include "PSCOracle.hxx"
#include "PSCPrimal.hxx"
#include "Bigmatrix.hxx"
#include "PSCAffineModification.hxx"

//------------------------------------------------------------

namespace ConicBundle {

  /**@defgroup implemented_psc_oracle implemention of a PSCOracle (PSCAffineFunction)
    @brief PSCAffineFunction is an implementation of ConicBundle::PSCOracle for the minimization of the maximum eigenvalue of an affine matrix function or, equivalently, Lagrangian relaxation of Linear Programs over the cone of positive semidefinite matrices.

     The class PSCAffineFunction implements a general purpose version of
     PSCOracle for minimizing the (\f$\gamma\f$-weighted) maximum eigenvalue
     \f$\gamma\lambda_{\max} F(y)\f$ of an affine matrix function
     \f$F(y)=C+\sum_{i=1}^my_iA_i\f$ , see \ref abstract_psc_oracle for
     an explanation of the general setting and the connections
     to Lagrangian relaxation of SDPs of the form

     \f[
     \mbox{maximize }\langle C,X\rangle \mbox{ subject to } \langle
     -A_i,X\rangle \begin{array}{c}\le\\=\\\ge\end{array} b_i,~~ \langle
     I,X\rangle \begin{array}{c}\le\\=\end{array}\gamma,~~ X\succeq 0
     \f]

     PSCAffineFunction supports a rich variety of special choices for the
     matrices \f$C\f$ and \f$A_i\f$ and uses (for large matrices) an iterative
     Lanczos method to compute eigenvalues and eigenvectors of \f$F(y)\f$.  In
     particular, the matrices \f$C\f$ and \f$A_i\f$ may consist of a single or
     several diagonal blocks (the dimensions of the blocks have to be consistent
     for \f$C\f$ and the \f$A_i\f$), each block being a coefficient matrix of
     the abstract class ConicBundle::Coeffmat . The following versions of this
     abstract class are implemented so far and may be used directly as blocks

     - CMsymdense, a general purpose dense symmetric matrix building upon
       CH_Matrix_Classes::Symmatrix

     - CMsymsparse, a sparse symmetric matrix building upon
       CH_Matrix_Classes::Sparsesym

     - CMgramdense, a matrix \f$BB^T\f$ with \f$B\f$ a dense rectangular
       CH_Matrix_Classes::Matrix

     - CMgramsparse, a matrix \f$BB^T\f$ with \f$B\f$ a sparse rectangular
       CH_Matrix_Classes::Sparsemat

     - CMgramsparse_withoutdiag, a matrix \f$BB^T-\mbox{Diag}(BB^T)\f$ for a
       rectangular CH_Matrix_Classes::Sparsemat (useful in some quadratic 0-1
       constraint representations)

     - CMsingleton, a symmetric matrix with just one nonzero entry (by symmetry
       two if not on the diagonal)

     - CMlowrankdd, a matrix \f$AB+BA^T\f$ with each of \f$A,B\f$ a dense
       rectangular CH_Matrix_Classes::Matrix

     - CMlowrankss, a matrix \f$AB+BA^T\f$ with each of \f$A,B\f$ a sparse
       rectangular CH_Matrix_Classes::Sparsemat

     - CMlowranksd, a matrix \f$AB+BA^T\f$ with \f$A\f$ a sparse rectangular
       CH_Matrix_Classes::Sparsemat and \f$B\f$ a dense rectangular
       CH_Matrix_Classes::Matrix

     The block diagonal coefficient matrices \f$C\f$ and \f$A_i\f$ are
     each organized as instances of a ConicBundle::SparseCoeffmatMatrix
     which has the blocks as rows and each of its columns corresponds
     to one symmetric block matrix consisting of the corresponding
     diagonal blocks. In fact, it is more useful to
     think of the diagonal blocks \f$j=1,\dots,k\f$ as corresponding to
     separate semidefinite variables \f$X_j\f$, so that for each block
     \f$j\f$ there are coefficient matrices \f$C_j\f$ and \f$A_{ji}\f$,
     \f$i=1,\dots,m\f$ (many possibly of value zero)

     \f[
     X=\left[\begin{array}{@{}c@{}c@{}c@{}} X_1 &  & 0 \\[-1ex]  & \ddots & \\[-.5ex] 0 & &X_k\end{array}\right],\quad
     C=\left[\begin{array}{@{}c@{}c@{}c@{}} C_1 &  & 0 \\[-1ex]  & \ddots & \\[-.5ex] 0 & &C_k\end{array}\right],\quad
     A_i=\left[\begin{array}{@{}c@{}c@{}c@{}} A_{1i} &  & 0 \\[-1ex]  & \ddots & \\[-.5ex] 0 & &A_{ki}\end{array}\right]~ (i=1,\dots,m),
     \f]

     resulting in an SDP problem

     \f[ \mbox{maximize }\sum_{j=1}^k\langle C_j,X_j\rangle
     \mbox{ subject to } \sum_{j=1}^k\langle -A_{ji},X_j\rangle
     \begin{array}{c}\le\\=\\\ge\end{array} b_i~ (i=1,\dots,m),~~ \sum_{j=1}^k\langle I,X_j\rangle
     \begin{array}{c}\le\\=\end{array}\gamma,~~ X_j\succeq 0~ (j=1,\dots,k)\f]

     In this sense  ConicBundle::SparseCoeffmatMatrix is simply a representation
     of a sparse matrix having a (pointer to a) matrix \f$A_{ji}\f$ as element \f$(j,i)\f$. In PSCAffineFunction this matrix is called opAt.

     The standard use of an PSCAffineFunction is to fully initialize it on
     construction with matrices C and opAt of type SparseCoeffmatMatrix (and
     maybe a generating primal, see below) and then to add this function to the
     MatrixCBSolver by MatrixCBSolver::add_function(). If dynamic changes to this PSCAffineFunction are
     required afterwards, use the class ConicBundle::PSCAffineModification within the
     corresponding problem modification routines MatrixCBSolver::append_variables(),
     MatrixCBSolver::reassign_variables(), MatrixCBSolver::delete_variables()
     of the MatrixCBSolver interface.

     In order to allow the solver (or rather PSCModel) to aggregate the
     eigenvector information to primal approximations of the primal semidefinite
     \f$X\f$ (or the semidefinite blocks \f$X_j\f$), one may install
     a pointer to a generating PSCPrimal on construction or via applying
     PSCAffineModification::add_reset_generating_primal(). Several versions of
     PSCPrimal are implemented and ready to use:

     - DensePSCPrimal, general purpose dense \f$X\f$ based on
       CH_Matrix_Classes::Symmatrix, this is only useful for small matrices or
       small blocks, but it also keeps all information.

     - SparsePSCPrimal, it is initialized with a nonzero support structure and
       then the matrix entries are only collected on this support. The
       implementation is based on CH_Matrix_Classes::Sparsesym and is
       particularly useful if all matrices \f$C\f$ and \f$A_i\f$ have a common
       small support, because then all inner products are formed correctly with
       this reduced support of the true \f$X\f$.

     - GramSparsePSCPrimal, this is tuned for best representing the current
       optimal solution of the bundle subproblem of the solver (or rather
       PSCModel) by a matrix of the form \f$X\approx PP^T+S\f$ where \f$P\f$
       is the part due to the subspace description of the model and \f$S\f$ is the
       sparse part collected for the aggregate matrices of the model on a
       selected support as in SparsePSCPrimal.

     - BlockPSCPrimal allows to represent the primal information of a diagonal
       block partition of the full \f$X\f$ by any choice of an PSCPrimal, in
       particular the ones given before. if \f$X\f$ is itself a block matrix,
       typically the block partition of BlockPSCPrimal will be chosen to match
       the one of \f$X\f$ but it is not required to match. One main purpose in
       this choice is to allow inner products to work correctly
       for any new contraints \f$A_i\f$ that are added dynamically.

     When the solver (or rather PSCModel) calls the implementation of
     PSCOracle::evaluate() for some \f$y\f$, PSCAffineFunction collects
     the affine matrix function \f$F_j(y)=C_j+\sum y_iA_{ji}\f$ of each block
     \f$j\f$ in a separate Bigmatrix, the latter is derived from a
     CH_Matrix_Classes::Lanczosmatrix so as to serve as input for the iterative
     Lanczos eigenvalue solver CH_Matrix_Classes::Lanczpol. Then
     PSCAffineFunction starts a separate eigensolver (but so far not in
     parallel) for the bigmatrix of each block. If the block is small, a
     standard eigenvalue method is used, otherwise the Lanczos method is
     employed. In any case the method will not only yield the maximum
     eigenvector (or a Ritz vector with sufficiently large Ritz value to ensure
     a null step), but several other Ritz vectors as well, that might help to
     improve the quality of the model in PSCModel. The routine
     PSCAffineFunction::set_max_Ritzvecs() allows to specify how many of the
     Ritz vectors should be passed on to PSCModel as a result of the call to
     PSCAffineFunction::evaluate().

     For facilitating input and output PSCAffineFunction offers
     - PSCAffineFunction::print_problem_data() that outputs the full function description
       so that it can be read again by read_problem_data
     - PSCAffineFunction::read_problem_data() reads the problem data as output by print_problem_data()
     - PSCAffineFunction::set_out() and PSCAffineFunction::set_cbout() work as described in ConicBundle::CBout

  */

  //@{

  class PSCAffineFunction;

  /** @brief Implementation of MinorantExtender for PSCAffineFunction

      This object will be returned as an object on the heap
      by PSCAffineFunction::apply_modification
      and will be deleted by ConicBundle after use.
      Its purpose is to fill up further coordinates of the minorant
      on basis of the primal information stored in the minorants (if
      this is possible).

  */

  class PSCAffineMinorantExtender : public CBout, public MinorantExtender {
  private:
    /// the oracle this MinorantExtender was generated by, needed for retrieving problem data
    PSCAffineFunction* amf;

  public:
    ///the PSCAffineFunction pointed to has to be valid for the livetime of this object
    PSCAffineMinorantExtender(PSCAffineFunction* amf);

    ///
    ~PSCAffineMinorantExtender();

    /*@brief called by ConicBundle to update internal Minorant objects, has to return 0 on success

        @param[in,out] minorant  (Minorant&)
            it holds a (possibly aggregated) minorant that was generated
            from minorants returned by oracle calls, e.g. as in
      FunctionOracle::evaluate() If PrimalData was provided in these
      minorants, this will be aggregated along and will also be
      available in this minorant.

        @param[in] n_coords (int)
            the number of coordinate positions that have to be filled in

        @param[out] new_subgradient_values  (DVector &)
      the indices of these coordinate positions (sorted in
      strictly increasing order)

  @return
           -  0 on success,
           -  1 if extension/update is impossible

     */
    int extend(Minorant& minorant, int n_coords, const int* indices);

  };



  /// forward definition, helps to avoid fixing the eigenvalue solver implementation in the header 
  class AMFMaxEigSolver;

  /**@brief general purpose implementation of PSCOracle as explained in \ref implemented_psc_oracle
   */

  class PSCAffineFunction : public PSCOracle, public CBout {
  private:
    SparseCoeffmatMatrix C; ///< block diagonal representation of \f$C\f$ as in \ref implemented_psc_oracle

    /// holds the coefficient matrices for the variables
    SparseCoeffmatMatrix opAt;

    /// This PSCPrimal can be set from outside and serves for generating further primals by cloning 
    PSCPrimal* generating_primal;

    std::vector<AMFMaxEigSolver*> maxeigsolver; ///< the eigenvalue solver of each block

    CH_Matrix_Classes::Matrix last_bigmat_y;  ///< the last y used in costructing the matrix of the affine matrix function

    CH_Matrix_Classes::Integer maxvecs; ///< the maximum number of Ritz vectors to be retrurned

    bool check_correctness_flag; ///< if true, ConicBundle employs some additional consistency checks 

    /// compute the Bigmatrix representation for the given point 
    int form_bigmatrix(const CH_Matrix_Classes::Matrix& current_point);

    /// applies the PSCAffineModfication amfmod to the current function 
    int apply_modification(const PSCAffineModification& amfmod);

    /// resets all to the initial empty state 
    void clear();

  public:
    /**@name Initialization and setting parameters*/
    //@{
    ///
    PSCAffineFunction(const CBout* cb = 0, int incr = -1);

    /** @brief initialize the PSCAffineFunction with its matrices and possible a generating_primal

  C and opAt define the constant (block-)offset and the linear
  (block-)matrix function as described in the general text of
  PSCAffineFunction

  generating_primal defines in what form primal matrices should be
  aggregated. If the argument is NULL then no primal aggregation will
  take place.  The control over the generating primal is passed over to
  this. This will delete an existing generating primal whenever a new
  generating primal is set or upon destruction.

  The final two arguments allow to set the output, see CBout.
    */
    PSCAffineFunction(const SparseCoeffmatMatrix& C,
      const SparseCoeffmatMatrix& opAt,
      PSCPrimal* generating_primal = 0,
      const CBout* cb = 0, int incr = -1);

    ///
    ~PSCAffineFunction();

    /// set the maximum number of Ritz vectors returned in evaluate(), see \ref implemented_psc_oracle

    /// if set to true, ConicBundle employs some additional consistency checks
    void set_check_correctness(bool chk) {
      check_correctness_flag = chk;
    }

    /// set the maximum number of new Ritzvectors returned by evaluate(); values<1 default to 5
    void set_max_Ritzvecs(CH_Matrix_Classes::Integer maxv) {
      maxvecs = (maxv > 1) ? maxv : 5;
    }

    //@}

    //----------- Oracle Implementation of PSCOracle ----------
    /**@name Implementations of PSCOracle routines */
    //@{

    /// see PSCOracle::generate_minorant()
    Minorant* generate_minorant(const CH_Matrix_Classes::Matrix& P);

    /// see PSCOracle::svec_projection() 
    int svec_projection(CH_Matrix_Classes::Matrix& svec_offset,
      CH_Matrix_Classes::Matrix& svec_coeffs,
      const CH_Matrix_Classes::Matrix& P,
      const CH_Matrix_Classes::Indexmatrix* index_subset = 0);



    /// see PSCOracle::evaluate()
    int evaluate(const CH_Matrix_Classes::Matrix& current_point,
      const CH_Matrix_Classes::Matrix& bundlevecs,
      const double relprec,
      const double Ritz_bound,
      CH_Matrix_Classes::Matrix& Ritz_vectors,
      CH_Matrix_Classes::Matrix& Ritz_values,
      PSCPrimalExtender*& primal_extender);

    /// see PSCOracle::evaluate_projection()
    int evaluate_projection(const CH_Matrix_Classes::Matrix& current_point,
      const CH_Matrix_Classes::Matrix& P,
      const double relprec,
      CH_Matrix_Classes::Matrix& projected_Ritz_vectors,
      CH_Matrix_Classes::Matrix& projected_Ritz_values);

    /// see PSCOracle::left_right_product()
    int left_right_product(int i,
      const CH_Matrix_Classes::Matrix& E,
      const CH_Matrix_Classes::Matrix& F,
      CH_Matrix_Classes::Matrix& G);


    /** @brief see PSCOracle::apply_modification() for the general use, here oracle_modification has a special role if it can be cast to an PSCAffineModification

  if oracle_modification cannot be cast to an PSCAffineModification it is assumed that all append modifications amount to have already been carried out on *this seperately before this routine is called. In particular, it is only checked whether the new dimension matches the one given by oracle_modification, the old dimension is ignored. If this does not hold, the routine stops with an error. Otherwise it checks the other stuff as if a suitable PSCAffineModification has just been executed.
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


    /// see PSCOracle::check_correctness()
    virtual
      bool
      check_correctness() const {
      return check_correctness_flag;
    }
    //@}


    //------------------  routines for querying the problem ----------------
    /**@name routines for querying data of the problem */
    //@{

    /** returns the row representation of the coefficient matrices
     (each entry of the map represents a row by a SparseCoeffmatVector).
    */
    const SparseCoeffmatMatrix& get_opAt() {
      return opAt;
    }

    /** returns the block representation of the coefficient matrices
     (each entry of the map represents a block by a SparseCoeffmatVector).
    */
    const SparseCoeffmatMatrix& get_C() {
      return C;
    }

    /// returns the current setting concerning the generation of an PSCPrimal (0 for none)
    const PSCPrimal* get_generating_primal(void) {
      return generating_primal;
    }

    //@}


   //------------------  routines for IO ----------------
    /**@name routines for supporting input and output */
    //@{

    /// see ConicBundle::CBout
    void  set_out(std::ostream* o = 0, int pril = 1);
    /// see ConicBundle::CBout
    void  set_cbout(const CBout* cb, int incr = -1);

    /// write the problem description to out so that it can be read again by read_problem_data()
    std::ostream& print_problem_data(std::ostream& out) const;

    /// clear() and read the problem from in in the format written by print_problem_data()
    std::istream& read_problem_data(std::istream& in);

    /// undocumented highly volatile variant for external testing 
    std::ostream& print_problem_data_to_mfile(std::ostream& out, CH_Matrix_Classes::Integer blocknr) const;
    //std::istream& read_problem_data_from_mfile(std::istream& in);


  };


}

//@}

#endif

