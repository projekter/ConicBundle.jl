/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/CMgramsparse_withoutdiag.hxx
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



#ifndef CONICBUNDLE_CMGRAMSPARSE_WITHOUTDIAG_HXX
#define CONICBUNDLE_CMGRAMSPARSE_WITHOUTDIAG_HXX

/**  @file CMgramsparse_withoutdiag.hxx
    @brief Header declaring the class ConicBundle::CMgramsparse_withoutdiag (needed for ConicBundle::AffineMatrixFunction)
    @version 1.0
    @date 2015-03-24
    @author Christoph Helmberg
*/

#include "Coeffmat.hxx"
#include "CMsymdense.hxx"

namespace ConicBundle {

  /** @ingroup implemented_psc_oracle
   */
   //@{

    /**@brief implements a Gram matrix \f$\pm (AA^T-\mbox{Diag}(AA^T))\f$ with zero diagonal as Coeffmat for a sparse rectangular CH_Matrix_Classes::Sparsemat \f$A\f$ (for use with MatrixSDPfunction, see \ref implemented_psc_oracle).

        As for CMgramsparse there is an extra
        flag for using the positive or the  negative version of the matrix.
   */

  class CMgramsparse_withoutdiag : public Coeffmat {
  private:
    CH_Matrix_Classes::Sparsemat A;  ///< the Coeffmat acts like A*A^T-Diag(A*A^T) or its negative 
    CH_Matrix_Classes::Sparsesym di; ///< this is the precomputed diagonal di=diag(A*A^T)
    bool positive; ///< if true use A*A^T-Diag, if false use -A*A^T+Diag
  public:
    /// copy Ain, the flag for positive/negative and store the user information
    CMgramsparse_withoutdiag(const CH_Matrix_Classes::Sparsemat& Ain, bool pos = true, CoeffmatInfo* cip = 0) {
      A = Ain; di = sparseDiag(rowsip(A)); positive = pos; CM_type = CM_gramsparsewd; infop = cip;
    }
    ///
    virtual ~CMgramsparse_withoutdiag() {
    }


    ///makes an explicit copy of itself and returns a pointer to it 
    virtual Coeffmat* clone() const {
      return new CMgramsparse_withoutdiag(A, positive, ConicBundle::clone(infop));
    }

    ///returns the order of the represented symmetric matrix
    virtual CH_Matrix_Classes::Integer dim() const {
      return A.rowdim();
    }

    ///returns the value of the matrix element (i,j)
    virtual CH_Matrix_Classes::Real operator()(CH_Matrix_Classes::Integer i, CH_Matrix_Classes::Integer j) const {
      if (i == j) return 0.; return  (positive ? 1. : -1.) * CH_Matrix_Classes::ip(A.row(i), A.row(j));
    }

    ///returns a dense symmetric constraint matrix
    virtual void make_symmatrix(CH_Matrix_Classes::Symmatrix& S) const {
      if (positive) {
        CH_Matrix_Classes::rankadd(A, S);
        for (CH_Matrix_Classes::Integer i = 0; i < di.get_colindex().dim(); i++) {
          CH_Matrix_Classes::Integer j = di.get_colindex()(i);
          S(j, j) -= di.get_colval()(i);
        }
      } else {
        CH_Matrix_Classes::rankadd(A, S, -1.);
        for (CH_Matrix_Classes::Integer i = 0; i < di.get_colindex().dim(); i++) {
          CH_Matrix_Classes::Integer j = di.get_colindex()(i);
          S(j, j) += di.get_colval()(i);
        }
      }
    }

    ///returns the Frobenius norm of the matrix
    virtual CH_Matrix_Classes::Real norm(void) const {
      CH_Matrix_Classes::Real sum = 0.; CH_Matrix_Classes::Matrix vec;
      for (CH_Matrix_Classes::Integer i = 0; i < A.get_rowinfo().rowdim(); i++) {
        CH_Matrix_Classes::genmult(A, A.row(A.get_rowinfo()(i, 0)), vec, 1., 0., 0, 1);
        sum += CH_Matrix_Classes::ip(vec, vec) - di(A.get_rowinfo()(i, 0), A.get_rowinfo()(i, 0));
      }
      return std::sqrt(sum);
    }

    ///delivers a new object on the heap corresponding to the matrix P^T(*this)P, the caller is responsible for deleting the object
    virtual Coeffmat* subspace(const CH_Matrix_Classes::Matrix& P) const {
      CH_Matrix_Classes::Matrix B;
      CH_Matrix_Classes::genmult(P, A, B, 1., 0., 1);
      CH_Matrix_Classes::Symmatrix S;
      CH_Matrix_Classes::rankadd(B, S);
      CH_Matrix_Classes::genmult(P, di, B, 1., 0., 1);
      S -= CH_Matrix_Classes::Symmatrix(B * P);
      if (!positive) S *= -1;
      return new CMsymdense(S, ConicBundle::clone(infop));
    }

    ///multiply constraint permanentely by d; this is to allow scaling or sign changes in the constraints
    virtual void multiply(CH_Matrix_Classes::Real d) {
      di *= std::abs(d);
      if (d < 0.) {
        A *= sqrt(-d); positive = !positive;
      } else {
        A *= sqrt(d);
      }
      if (infop) infop->multiply(d);
    }

    ///returns ip(*this,S)=trace(*this*S), the trace inner product
    virtual CH_Matrix_Classes::Real ip(const CH_Matrix_Classes::Symmatrix& S) const {
      CH_Matrix_Classes::Matrix B; if (positive) return CH_Matrix_Classes::ip(CH_Matrix_Classes::genmult(S, A, B), A) - CH_Matrix_Classes::ip(S, di);
      else return -CH_Matrix_Classes::ip(CH_Matrix_Classes::genmult(S, A, B), A) + CH_Matrix_Classes::ip(S, di);
    }

    ///returns ip(*this,PP^T)=trace P^T(*this)P
    virtual CH_Matrix_Classes::Real gramip(const CH_Matrix_Classes::Matrix& P) const {
      CH_Matrix_Classes::Real sum = 0.;
      for (CH_Matrix_Classes::Integer i = 0; i < di.get_colindex().dim(); i++) {
        CH_Matrix_Classes::Matrix vec = P.row(di.get_colindex()(i));
        sum += CH_Matrix_Classes::ip(vec, vec) * di.get_colval()(i);
      }
      CH_Matrix_Classes::Matrix B; CH_Matrix_Classes::genmult(P, A, B, 1., 0., 1);
      if (positive) return CH_Matrix_Classes::ip(B, B) - sum;
      else return -CH_Matrix_Classes::ip(B, B) + sum;
    }

    ///returns ip(*this,QQ^T)=trace Q^T(*this)Q for Q=P.rows(start_row,start_row+dim-1) 
    virtual CH_Matrix_Classes::Real gramip(const CH_Matrix_Classes::Matrix& P, CH_Matrix_Classes::Integer start_row, const CH_Matrix_Classes::Matrix* Lam = 0) const {
      CH_Matrix_Classes::Real trval = 0.;
      for (CH_Matrix_Classes::Integer i = 0; i < di.get_colindex().dim(); i++) {
        CH_Matrix_Classes::Matrix vec = P.row(di.get_colindex()(i) + start_row);
        trval -= CH_Matrix_Classes::rowip(P, start_row + di.get_colindex()(i), Lam) * di.get_colval()(i);
      }
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
      if (Lam == 0) {
        trval += CH_Matrix_Classes::ip(B, B);
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
      if (!positive) d *= -1;
      CH_Matrix_Classes::rankadd(A, S, d, 1.); S -= di * d;
    }

    ///computes B+=d*(*this)*C
    virtual void addprodto(CH_Matrix_Classes::Matrix& B, const CH_Matrix_Classes::Matrix& C, CH_Matrix_Classes::Real d = 1.) const {
      CH_Matrix_Classes::Matrix D; if (!positive) d *= -1;
      CH_Matrix_Classes::genmult(A, CH_Matrix_Classes::genmult(A, C, D, 1., 0., 1), B, d, 1.);
      CH_Matrix_Classes::genmult(di, C, B, -d, 1., 0);
    }

    ///computes B+=d*(*this)*C
    virtual void addprodto(CH_Matrix_Classes::Matrix& B, const CH_Matrix_Classes::Sparsemat& C, CH_Matrix_Classes::Real d = 1.) const {
      if (!positive) d *= -1;
      CH_Matrix_Classes::Matrix D;
      CH_Matrix_Classes::genmult(A, CH_Matrix_Classes::genmult(A, C, D, 1., 0., 1), B, d, 1.);
      CH_Matrix_Classes::genmult(di, C, B, -d, 1., 0);
    }

    /// computes R=P^T*(*this)*Q
    virtual void left_right_prod(const CH_Matrix_Classes::Matrix& P, const CH_Matrix_Classes::Matrix& Q, CH_Matrix_Classes::Matrix& R) const {
      R.init(P.coldim(), Q.coldim(), 0.);
      for (CH_Matrix_Classes::Integer i = 0; i < di.get_colindex().dim(); i++) {
        genmult(P.row(di.get_colindex()(i)), Q.row(di.get_colindex()(i)), R, di.get_colval()(i), 1., 1);
      }
      CH_Matrix_Classes::Matrix tmp1; CH_Matrix_Classes::genmult(P, A, tmp1, 1., 0., 1, 0);
      CH_Matrix_Classes::Matrix tmp2; CH_Matrix_Classes::genmult(A, Q, tmp2, 1., 0., 1, 0);
      if (positive) CH_Matrix_Classes::genmult(tmp1, tmp2, R, 1., -1., 0, 0);
      else CH_Matrix_Classes::genmult(tmp1, tmp2, R, -1., 1., 0, 0);
    }

    ///returns an estimate of number of flops to compute addprodto for a vector
    virtual CH_Matrix_Classes::Integer prodvec_flops() const {
      return 4 * A.nonzeros() + 2 * di.get_colindex().dim();
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
      if (positive) return CH_Matrix_Classes::ip(A, S * A) - CH_Matrix_Classes::ip(di, S); else return -CH_Matrix_Classes::ip(A, S * A) + CH_Matrix_Classes::ip(di, S);
    }

    ///computes S=P^T*(*this)*P
    virtual void project(CH_Matrix_Classes::Symmatrix& S, const CH_Matrix_Classes::Matrix& P) const {
      CH_Matrix_Classes::Matrix B; CH_Matrix_Classes::genmult(P, A, B, 1., 0., 1);
      if (positive) CH_Matrix_Classes::rankadd(B, S); else CH_Matrix_Classes::rankadd(B, S, -1.);
      for (CH_Matrix_Classes::Integer i = 0; i < di.get_colindex().dim(); i++) {
        rankadd(P.row(di.get_colindex()(i)), S, (positive ? -1. : 1.) * di.get_colval()(i), 1., 1);
      }
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
      for (CH_Matrix_Classes::Integer i = 0; i < di.get_colindex().dim(); i++) {
        rankadd(P.row(di.get_colindex()(i) + start_row), S, (positive ? -alpha : alpha) * di.get_colval()(i), 1., 1);
      }
    }

    ///computes C= alpha*(*this)*B^(T if btrans) + beta*C, C is also returned
    virtual const CH_Matrix_Classes::Matrix& postgenmult(const CH_Matrix_Classes::Matrix& B, CH_Matrix_Classes::Matrix& C,
      CH_Matrix_Classes::Real alpha = 1., CH_Matrix_Classes::Real beta = 0., int btrans = 0) const {
      CH_Matrix_Classes::Matrix D;
      CH_Matrix_Classes::genmult(A, CH_Matrix_Classes::genmult(A, B, D, 1., 0., 1, btrans), C, (positive ? 1. : -1.) * alpha, beta);
      return CH_Matrix_Classes::genmult(di, B, C, (positive ? -1. : 1.) * alpha, 1., btrans);
    }

    ///computes C= alpha*B^(T if btrans)*(*this) + beta*C, C is also returned
    virtual const CH_Matrix_Classes::Matrix& pregenmult(const CH_Matrix_Classes::Matrix& B, CH_Matrix_Classes::Matrix& C,
      CH_Matrix_Classes::Real alpha = 1., CH_Matrix_Classes::Real beta = 0., int btrans = 0) const {
      CH_Matrix_Classes::Matrix D;
      return CH_Matrix_Classes::genmult(CH_Matrix_Classes::genmult(B, A, D, 1., 0., btrans), A, C, (positive ? 1. : -1.) * alpha, beta, 0, 1);
      return CH_Matrix_Classes::genmult(B, di, C, (positive ? -1. : 1.) * alpha, 1., btrans);
    }

    ///returns 1, if p is the same derived class and entries differ by less than tol, otherwise zero
    virtual int equal(const Coeffmat* p, double tol = 1e-6) const {
      const CMgramsparse_withoutdiag* pp = dynamic_cast<const CMgramsparse_withoutdiag*>(p);
      if (pp == 0)
        return 0;
      return (positive == pp->positive) && CH_Matrix_Classes::equal(A, pp->A, tol);
    }

    ///display constraint information
    virtual std::ostream& display(std::ostream& o) const {
      o << "CMgramsparse_withoutdiag\n"; A.display(o); return o;
    }

    ///put entire contents onto ostream with the class type in the beginning so that the derived class can be recognized by in().
    virtual std::ostream& out(std::ostream& o) const {
      return o << "GRAM_SPARSE_WITHOUTDIAG\n" << positive << "\n" << A;
    }

    ///counterpart to out(), does not read the class type, though. This is assumed to have been read in order to generate the correct class
    virtual std::istream& in(std::istream& i) {
      i >> positive >> A; di = sparseDiag(rowsip(A)); return i;
    }

    /// constructor with istream and possibly additional user information
    CMgramsparse_withoutdiag(std::istream& is, CoeffmatInfo* cip = 0) {
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

