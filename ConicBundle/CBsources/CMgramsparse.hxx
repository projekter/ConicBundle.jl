/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/CMgramsparse.hxx
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



#ifndef CONICBUNDLE_CMGRAMSPARSE_HXX
#define CONICBUNDLE_CMGRAMSPARSE_HXX

/**  @file CMgramsparse.hxx
    @brief Header declaring the class ConicBundle::CMgramsparse (needed for ConicBundle::MatrixSDPfunction)
    @version 1.0
    @date 2015-03-24
    @author Christoph Helmberg
*/

#include "Coeffmat.hxx"
#include "CMgramdense.hxx"

namespace ConicBundle {


  /** @ingroup implemented_psc_oracle
   */
   //@{

    /**@brief implements a Gram matrix \f$\pm AA^T\f$ as Coeffmat for a sparse rectangular CH_Matrix_Classes::Sparsemat \f$A\f$ (for use with MatrixSDPfunction, see \ref implemented_psc_oracle).

        While a Gram matrix by itself is always positive semidefinite, there is an extra
        flag for using its negative version.
   */


  class CMgramsparse : public Coeffmat {
  private:
    CH_Matrix_Classes::Sparsemat A;  ///< the Coeffmat acts like A*A^T or its negative
    bool positive; ///< if true use A*A^T, if false use -A*A^T
  public:
    /// copy Ain, the flag for positive/negative and store the user information
    CMgramsparse(const CH_Matrix_Classes::Sparsemat& Ain, bool pos = true, CoeffmatInfo* cip = 0) {
      A = Ain; positive = pos; CM_type = CM_gramsparse; infop = cip;
    }
    ///
    virtual ~CMgramsparse() {
    }

    ///makes an explicit copy of itself and returns a pointer to it 
    virtual Coeffmat* clone() const {
      return new CMgramsparse(A, positive, ConicBundle::clone(infop));
    }

    ///returns the order of the represented symmetric matrix
    virtual CH_Matrix_Classes::Integer dim() const {
      return A.rowdim();
    }

    ///returns the value of the matrix element (i,j)
    virtual CH_Matrix_Classes::Real operator()(CH_Matrix_Classes::Integer i, CH_Matrix_Classes::Integer j) const {
      return  (positive ? 1. : -1.) * CH_Matrix_Classes::ip(A.row(i), A.row(j));
    }

    ///returns a dense symmetric constraint matrix
    virtual void make_symmatrix(CH_Matrix_Classes::Symmatrix& S) const {
      if (positive) {
        CH_Matrix_Classes::rankadd(A, S);
      } else {
        CH_Matrix_Classes::rankadd(A, S, -1.);
      }
    }

    ///returns the Frobenius norm of the matrix
    virtual CH_Matrix_Classes::Real norm(void) const {
      CH_Matrix_Classes::Symmatrix B; return CH_Matrix_Classes::norm2(CH_Matrix_Classes::rankadd(A, B, 1., 0., 1));
    }

    ///delivers a new object on the heap corresponding to the matrix P^T(*this)P, the caller is responsible for deleting the object
    virtual Coeffmat* subspace(const CH_Matrix_Classes::Matrix& P) const {
      CH_Matrix_Classes::Matrix B; return new CMgramdense(CH_Matrix_Classes::genmult(P, A, B, 1., 0., 1), positive, ConicBundle::clone(infop));
    }

    ///multiply constraint permanentely by d; this is to allow scaling or sign changes in the constraints
    virtual void multiply(CH_Matrix_Classes::Real d) {
      if (d < 0.) {
        A *= sqrt(-d); positive = !positive;
      } else {
        A *= sqrt(d);
      }
      if (infop) infop->multiply(d);
    }

    ///returns ip(*this,S)=trace(*this*S), the trace inner product
    virtual CH_Matrix_Classes::Real ip(const CH_Matrix_Classes::Symmatrix& S) const {
      CH_Matrix_Classes::Matrix B; if (positive) return CH_Matrix_Classes::ip(CH_Matrix_Classes::genmult(S, A, B), A);
      else return -CH_Matrix_Classes::ip(CH_Matrix_Classes::genmult(S, A, B), A);
    }

    ///returns ip(*this,PP^T)=trace P^T(*this)P
    virtual CH_Matrix_Classes::Real gramip(const CH_Matrix_Classes::Matrix& P) const {
      CH_Matrix_Classes::Matrix B; CH_Matrix_Classes::genmult(P, A, B, 1., 0., 1);
      if (positive) return CH_Matrix_Classes::ip(B, B); else return -CH_Matrix_Classes::ip(B, B);
    }

    ///returns ip(*this,QQ^T)=trace Q^T(*this)Q for Q=P.rows(start_row,start_row+dim-1) 
    virtual CH_Matrix_Classes::Real gramip(const CH_Matrix_Classes::Matrix& P, CH_Matrix_Classes::Integer start_row, const CH_Matrix_Classes::Matrix* Lam = 0) const {
      CH_Matrix_Classes::Matrix B(A.coldim(), P.coldim(), 0.);
      CH_Matrix_Classes::Real* bp = B.get_store();
      const CH_Matrix_Classes::Integer brd = B.rowdim();
      const CH_Matrix_Classes::Integer pcd = P.coldim();
      const CH_Matrix_Classes::Integer prd = P.rowdim();
      const CH_Matrix_Classes::Real* pp = P.get_store() + start_row;
      const CH_Matrix_Classes::Indexmatrix& colinfo = A.get_colinfo();
      const CH_Matrix_Classes::Real* vp = (A.get_colval()).get_store();
      const CH_Matrix_Classes::Integer* ip = (A.get_colindex()).get_store();
      for (CH_Matrix_Classes::Integer j = 0; j < colinfo.rowdim(); j++) {
        CH_Matrix_Classes::Integer indj = colinfo(j, 0);
        for (CH_Matrix_Classes::Integer i = colinfo(j, 1); --i >= 0;) {
          CH_Matrix_Classes::mat_xpeya(pcd, bp + indj, brd, pp + (*ip++), prd, (*vp++));
        }
      }
      CH_Matrix_Classes::Real trval = 0;
      if (Lam == 0) {
        trval = CH_Matrix_Classes::ip(B, B);
      } else {
        assert(Lam->dim() == P.coldim());
        const CH_Matrix_Classes::Real* bp = B.get_store();
        const CH_Matrix_Classes::Real* lp = Lam->get_store();
        for (CH_Matrix_Classes::Integer i = 0; i < B.coldim(); i++, bp += B.rowdim())
          trval += (*lp++) * CH_Matrix_Classes::mat_ip(B.rowdim(), bp);
      }
      return (positive ? trval : -trval);
    }

    ///computes S+=d*(*this);
    virtual void addmeto(CH_Matrix_Classes::Symmatrix& S, CH_Matrix_Classes::Real d = 1.) const {
      if (positive) CH_Matrix_Classes::rankadd(A, S, d, 1.); else CH_Matrix_Classes::rankadd(A, S, -d, 1.);
    }

    ///computes B+=d*(*this)*C
    virtual void addprodto(CH_Matrix_Classes::Matrix& B, const CH_Matrix_Classes::Matrix& C, CH_Matrix_Classes::Real d = 1.) const {
      CH_Matrix_Classes::Matrix D;
      if (positive) CH_Matrix_Classes::genmult(A, CH_Matrix_Classes::genmult(A, C, D, 1., 0., 1), B, d, 1.);
      else CH_Matrix_Classes::genmult(A, CH_Matrix_Classes::genmult(A, C, D, 1., 0., 1), B, -d, 1.);
    }

    ///computes B+=d*(*this)*C
    virtual void addprodto(CH_Matrix_Classes::Matrix& B, const CH_Matrix_Classes::Sparsemat& C, CH_Matrix_Classes::Real d = 1.) const {
      CH_Matrix_Classes::Matrix D;
      if (positive) CH_Matrix_Classes::genmult(A, CH_Matrix_Classes::genmult(A, C, D, 1., 0., 1), B, d, 1.);
      else CH_Matrix_Classes::genmult(A, CH_Matrix_Classes::genmult(A, C, D, 1., 0., 1), B, -d, 1.);
    }

    /// computes R=P^T*(*this)*Q
    virtual void left_right_prod(const CH_Matrix_Classes::Matrix& P, const CH_Matrix_Classes::Matrix& Q, CH_Matrix_Classes::Matrix& R) const {
      CH_Matrix_Classes::Matrix tmp1; CH_Matrix_Classes::genmult(P, A, tmp1, 1., 0., 1, 0);
      CH_Matrix_Classes::Matrix tmp2; CH_Matrix_Classes::genmult(A, Q, tmp2, 1., 0., 1, 0);
      if (positive) CH_Matrix_Classes::genmult(tmp1, tmp2, R, 1., 0., 0, 0);
      else CH_Matrix_Classes::genmult(tmp1, tmp2, R, -1., 0., 0, 0);
    }

    ///returns an estimate of number of flops to compute addprodto for a vector
    virtual CH_Matrix_Classes::Integer prodvec_flops() const {
      return 4 * A.nonzeros();
    }

    ///returns 1 if its structure is as bad as its dense symmetric representation, otherwise 0
    virtual int dense() const {
      return 0;
    }

    ///returns 0 if not sparse, otherwise 1
    virtual int sparse() const {
      return 0;
    }

    /// returns 0 if not sparse. If it is sparse it returns 1 and the nonzero structure in I,J and v, where v is multiplied by d. Only the upper triangle (including diagonal) is delivered
    virtual int sparse(CH_Matrix_Classes::Indexmatrix& /* I */,
      CH_Matrix_Classes::Indexmatrix& /* J */,
      CH_Matrix_Classes::Matrix& /* val */,
      CH_Matrix_Classes::Real /* d=1. */)const {
      return 0;
    }

    /// returns 0 if the support of the costraint matrix is not contained in the support of the sparse symmetric matrix S, 1 if it is contained.
    virtual int support_in(const CH_Matrix_Classes::Sparsesym& /* S */) const {
      return 0;
    }

    ///returns the inner product of the constraint matrix with S
    virtual CH_Matrix_Classes::Real ip(const CH_Matrix_Classes::Sparsesym& S) const {
      if (positive) return CH_Matrix_Classes::ip(A, S * A); else return -CH_Matrix_Classes::ip(A, S * A);
    }

    ///computes S=P^T*(*this)*P
    virtual void project(CH_Matrix_Classes::Symmatrix& S, const CH_Matrix_Classes::Matrix& P) const {
      CH_Matrix_Classes::Matrix B; CH_Matrix_Classes::genmult(P, A, B, 1., 0., 1);
      if (positive) CH_Matrix_Classes::rankadd(B, S); else CH_Matrix_Classes::rankadd(B, S, -1.);
    }

    ///computes S+=Q^T(*this)Q for Q=P.rows(start_row,start_row+dim-1) 
    virtual void add_projection(CH_Matrix_Classes::Symmatrix& S, const CH_Matrix_Classes::Matrix& P, CH_Matrix_Classes::Real alpha = 1., CH_Matrix_Classes::Integer start_row = 0) const {
      CH_Matrix_Classes::Matrix B(A.coldim(), P.coldim(), 0.);
      CH_Matrix_Classes::Real* bp = B.get_store();
      const CH_Matrix_Classes::Integer brd = B.rowdim();
      const CH_Matrix_Classes::Integer pcd = P.coldim();
      const CH_Matrix_Classes::Integer prd = P.rowdim();
      const CH_Matrix_Classes::Real* pp = P.get_store() + start_row;
      const CH_Matrix_Classes::Indexmatrix& colinfo = A.get_colinfo();
      const CH_Matrix_Classes::Real* vp = (A.get_colval()).get_store();
      const CH_Matrix_Classes::Integer* ip = (A.get_colindex()).get_store();
      for (CH_Matrix_Classes::Integer j = 0; j < colinfo.rowdim(); j++) {
        CH_Matrix_Classes::Integer indj = colinfo(j, 0);
        for (CH_Matrix_Classes::Integer i = colinfo(j, 1); --i >= 0;) {
          CH_Matrix_Classes::mat_xpeya(pcd, bp + indj, brd, pp + (*ip++), prd, (*vp++));
        }
      }
      if (positive) CH_Matrix_Classes::rankadd(B, S, alpha, 1., 1); else CH_Matrix_Classes::rankadd(B, S, -alpha, 1., 1);
    }

    ///computes C= alpha*(*this)*B^(T if btrans) + beta*C, C is also returned
    virtual const CH_Matrix_Classes::Matrix& postgenmult(const CH_Matrix_Classes::Matrix& B, CH_Matrix_Classes::Matrix& C,
      CH_Matrix_Classes::Real alpha = 1., CH_Matrix_Classes::Real beta = 0., int btrans = 0) const {
      CH_Matrix_Classes::Matrix D;
      return CH_Matrix_Classes::genmult(A, CH_Matrix_Classes::genmult(A, B, D, 1., 0., 1, btrans), C, (positive ? 1. : -1.) * alpha, beta);
    }

    ///computes C= alpha*B^(T if btrans)*(*this) + beta*C, C is also returned
    virtual const CH_Matrix_Classes::Matrix& pregenmult(const CH_Matrix_Classes::Matrix& B, CH_Matrix_Classes::Matrix& C,
      CH_Matrix_Classes::Real alpha = 1., CH_Matrix_Classes::Real beta = 0., int btrans = 0) const {
      CH_Matrix_Classes::Matrix D;
      return CH_Matrix_Classes::genmult(CH_Matrix_Classes::genmult(B, A, D, 1., 0., btrans), A, C, (positive ? 1. : -1.) * alpha, beta, 0, 1);
    }

    ///returns 1, if p is the same derived class and entries differ by less than tol, otherwise zero
    virtual int equal(const Coeffmat* p, double tol = 1e-6) const {
      const CMgramsparse* pp = dynamic_cast<const CMgramsparse*>(p);
      if (pp == 0)
        return 0;
      return (positive == pp->positive) && CH_Matrix_Classes::equal(A, pp->A, tol);
    }

    ///display constraint information
    virtual std::ostream& display(std::ostream& o) const {
      o << "CMgramsparse\n"; A.display(o); return o;
    }

    ///put entire contents onto ostream with the class type in the beginning so that the derived class can be recognized by in().
    virtual std::ostream& out(std::ostream& o) const {
      return o << "GRAM_SPARSE\n" << positive << "\n" << A;
    }

    ///counterpart to out(), does not read the class type, though. This is assumed to have been read in order to generate the correct class
    virtual std::istream& in(std::istream& i) {
      return i >> positive >> A;
    }

    /// constructor with istream and possibly additional user information
    CMgramsparse(std::istream& is, CoeffmatInfo* cip = 0) {
      CM_type = CM_gramsparse; infop = cip; in(is);
    }

    //--- specific routines
    ///returns the const reference to the internal matrix A forming the Gram matrix
    const CH_Matrix_Classes::Sparsemat& get_A() const {
      return A;
    }
    ///returns the flag on whether the Gram matrix is used in positive or negative form
    bool get_positive() const {
      return positive;
    }
  };

  //@}

}
#endif

