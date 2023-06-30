/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/Bigmatrix.hxx
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


#ifndef CONICBUNDLE_BIGMATRIX_HXX
#define CONICBUNDLE_BIGMATRIX_HXX

/**  @file Bigmatrix.hxx
    @brief Header declaring the classes ConicBundle::Bigmatrix  (needed for ConicBundle::MatrixSDPfunction)
    @version 1.0
    @date 2015-03-24
    @author Christoph Helmberg
*/

#include "lanczos.hxx"
#include "SparseCoeffmatMatrix.hxx"

namespace ConicBundle {
  /** @ingroup implemented_psc_oracle
   */
   //@{

   /** @brief AffineMatrixFunction needs to compute the maximum eigenvalue of an affine matrix function \f$F(y)=C+\sum y_iA_i\f$. This class prepares \f$F(y)\f$ in useful form for iterative eigenvalue solvers.

      In AffineMatrixFunction the matrices \f$C\f$ and \f$A_i\f$ are of type Coeffmat, the
      information provided by these matrices is used to decide which parts should
      be represented in a sparse form and which parts in a dense form. The intention
      is that in the end matrix vector products can be computed efficiently, so
      that e.g. a Lanczcos algorithm can be used. For this purpose Bigmatrix is a
      publically derived CH_Matrix_Classes::Lanczosmatrix.

    */

  class Bigmatrix :
    public CH_Matrix_Classes::Lanczosmatrix,
    protected CH_Matrix_Classes::Memarrayuser {
  private:
    /// for collecting the indices of sparse data this seems useful 
    typedef std::map<CH_Matrix_Classes::Integer, CH_Matrix_Classes::Real> Indexhash;

    CH_Matrix_Classes::Real    tol; ///< tolerance for considering an element to be zero

    CH_Matrix_Classes::Integer dim; ///< size of the matrix

    CH_Matrix_Classes::Integer nz;  ///< total number of nonzeros (except diagonal)
    CH_Matrix_Classes::Indexmatrix colnz; ///< nonzeros of each column
    std::vector<CH_Matrix_Classes::Indexmatrix> rowind; ///< nonzero row indices in each column
    std::vector<CH_Matrix_Classes::Matrix> rowval; ///< corresponding nonzero row values

    CH_Matrix_Classes::Matrix di;   ///< diagonal elements 
    std::vector<Indexhash> rowhash;  ///< for finding nonzero indices 

    std::vector<CoeffmatPointer> mcp; ///< pointers to dense matrices 

    CH_Matrix_Classes::Matrix mcv; ///< multiplier values for the dense matrices

    mutable CH_Matrix_Classes::Integer nmult; ///< counts the number of Mat*vec operations

    //---
    bool use_dense; ///< usually false, set to true if matrix is sufficiently dense
    mutable CH_Matrix_Classes::Symmatrix symrep;  ///< if use_dense is true, the matrix is stored here
    bool symrep_init; ///< flag whether symrep is computed already

    //--- temporary variables for lanczosmult
    mutable CH_Matrix_Classes::Matrix At; ///< auxilliary variable used in matrix vector products
    mutable CH_Matrix_Classes::Matrix Bt; ///< auxilliary variable used in matrix vector products
    CH_Matrix_Classes::Indexmatrix collecti; ///< for collecting the sparse part
    CH_Matrix_Classes::Indexmatrix collectj; ///< for collecting the sparse part
    CH_Matrix_Classes::Matrix collectval; ///< for collecting the sparse part

  public:
    ///
    Bigmatrix();
    ///
    ~Bigmatrix();
    /// clear all data
    void clear();
    /// return number of matrix vector products
    CH_Matrix_Classes::Integer get_nmult() const {
      return nmult;
    }
    /// reset the counter of matrix vector products to zero
    void reset_nmult() {
      nmult = 0;
    }
    /// compute a good representation of \f$F(y)=C+\sum y_iA_i\f$ (if dense==true directly a dense one) 
    int init(const CH_Matrix_Classes::Matrix& yi, CH_Matrix_Classes::Integer indim, const CoeffmatPointer C,
      const SparseCoeffmatVector* Ai, const bool dense = false);

    /// sets the tolerance for considering computed values as zeros
    void set_tol(CH_Matrix_Classes::Real t) {
      tol = t;
    }

    /// return true if a dense representation was computed
    int get_dense() const {
      return use_dense;
    }

    /// the order of the matrix as required by CH_Matrix_Classes::Lanczosmatrix
    virtual CH_Matrix_Classes::Integer lanczosdim() const;

    /// the estimated flop count of a matrix vector multiplication as required by CH_Matrix_Classes::Lanczosmatrix
    virtual CH_Matrix_Classes::Integer lanczosflops() const;

    /// computes  B = bigMatrix * A (A and B must not refer to the same memory)  as required by CH_Matrix_Classes::Lanczosmatrix
    virtual int lanczosmult(const CH_Matrix_Classes::Matrix& A, CH_Matrix_Classes::Matrix& B) const;

    /// write a dense version of the current matrix to S
    int make_symmatrix(CH_Matrix_Classes::Symmatrix& S) const;

    /// compute the dense version locally and return a reference to it
    const CH_Matrix_Classes::Symmatrix& get_symrep() const {
      make_symmatrix(symrep); return symrep;
    }

    /// output information
    friend std::ostream& operator<<(std::ostream& out, const Bigmatrix&);
  };

  //@}

} //end namespace ConicBundle

#endif

