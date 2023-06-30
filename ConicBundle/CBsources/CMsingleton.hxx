/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/CMsingleton.hxx
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



#ifndef CONICBUNDLE_CMSINGLETON_HXX
#define CONICBUNDLE_CMSINGLETON_HXX

/**  @file CMsingleton.hxx
    @brief Header declaring the class ConicBundle::CMsingleton (needed for ConicBundle::AffineMatrixFunction)
    @version 1.0
    @date 2015-03-24
    @author Christoph Helmberg
*/


#include "CMgramdense.hxx"
#include "CMlowrankdd.hxx"
#include "sparssym.hxx"

namespace ConicBundle {

  /** @ingroup implemented_psc_oracle
   */
   //@{

    /**@brief implements a Coeffmat having just one nonzero element (or two by symmetry) for use with MatrixSDPfunction, see \ref implemented_psc_oracle.

   */


  class CMsingleton : public Coeffmat {
  private:
    CH_Matrix_Classes::Integer nr; ///< order/size of the matrix
    CH_Matrix_Classes::Integer ii; ///< row index i of the nonzero element
    CH_Matrix_Classes::Integer jj; ///< column index j of the nonzero element
    CH_Matrix_Classes::Real    val; ///< value of the element
  public:
    /// the order is innr and the nonzero element (ini,inj) has values inval
    CMsingleton(CH_Matrix_Classes::Integer innr, CH_Matrix_Classes::Integer ini, CH_Matrix_Classes::Integer inj, CH_Matrix_Classes::Real inval, CoeffmatInfo* cip = 0) {
      nr = innr; ii = ini; jj = inj; val = inval; CM_type = CM_singleton; infop = cip;
    }
    ///
    virtual ~CMsingleton() {
    }

    ///makes an explicit copy of itself and returns a pointer to it 
    virtual Coeffmat* clone() const {
      return new CMsingleton(nr, ii, jj, val, ConicBundle::clone(infop));
    }

    ///returns the order of the represented symmetric matrix
    virtual CH_Matrix_Classes::Integer dim() const {
      return nr;
    }

    ///returns the value of the matrix element (i,j)
    virtual CH_Matrix_Classes::Real operator()(CH_Matrix_Classes::Integer i, CH_Matrix_Classes::Integer j) const {
      if (((i == ii) && (j == jj)) || ((i == jj) && (j == ii))) return val; return 0.;
    }

    ///returns a dense symmetric constraint matrix (useful for testing)
    virtual void make_symmatrix(CH_Matrix_Classes::Symmatrix& S) const {
      S.init(nr, 0.); S(ii, jj) = val;
    }

    ///returns the Frobenius norm of the matrix
    virtual CH_Matrix_Classes::Real norm(void) const {
      if (ii == jj) return fabs(val); return sqrt(2.) * fabs(val);
    }

    ///delivers a new object on the heap corresponding to the matrix P^T(*this)P, the caller is responsible for deleting the object
    virtual Coeffmat* subspace(const CH_Matrix_Classes::Matrix& P) const {
      if (ii == jj) {
        return new CMgramdense(sqrt(fabs(val)) * (P.row(ii)).transpose(), val > 0, ConicBundle::clone(infop));
      }
      return new CMlowrankdd((P.row(ii)).transpose() * val, (P.row(jj)).transpose(), ConicBundle::clone(infop));
    }

    ///multiply constraint permanentely by d; this is to allow scaling or sign changes in the constraints
    virtual void multiply(CH_Matrix_Classes::Real d) {
      val *= d; if (infop) infop->multiply(d);
    }

    ///returns ip(*this,S)=trace(*this*S), the trace inner product
    virtual CH_Matrix_Classes::Real ip(const CH_Matrix_Classes::Symmatrix& S) const {
      if (ii == jj) return val * S(ii, ii); return 2. * val * S(ii, jj);
    }

    ///returns ip(*this,PP^T)=trace P^T(*this)P
    virtual CH_Matrix_Classes::Real gramip(const CH_Matrix_Classes::Matrix& P) const {
      if (ii == jj)
        return val * CH_Matrix_Classes::mat_ip(P.coldim(), P.get_store() + ii, P.rowdim(), P.get_store() + ii, P.rowdim());
      return 2. * val * CH_Matrix_Classes::mat_ip(P.coldim(), P.get_store() + ii, P.rowdim(), P.get_store() + jj, P.rowdim());
    }

    ///returns ip(*this,QQ^T)=trace Q^T(*this)Q for Q=P.rows(start_row,start_row+dim-1) 
    virtual CH_Matrix_Classes::Real gramip(const CH_Matrix_Classes::Matrix& P, CH_Matrix_Classes::Integer start_row, const CH_Matrix_Classes::Matrix* Lam = 0) const {
      if (Lam == 0) {
        if (ii == jj)
          return val * CH_Matrix_Classes::mat_ip(P.coldim(), P.get_store() + ii + start_row, P.rowdim(), P.get_store() + ii + start_row, P.rowdim());
        return 2. * val * CH_Matrix_Classes::mat_ip(P.coldim(), P.get_store() + ii + start_row, P.rowdim(), P.get_store() + jj + start_row, P.rowdim());
      } else {
        assert(Lam->dim() == P.coldim());
        if (ii == jj)
          return val * CH_Matrix_Classes::mat_ip(P.coldim(), P.get_store() + ii + start_row, P.rowdim(), P.get_store() + ii + start_row, P.rowdim(), Lam->get_store(), 1);
        return 2. * val * CH_Matrix_Classes::mat_ip(P.coldim(), P.get_store() + ii + start_row, P.rowdim(), P.get_store() + jj + start_row, P.rowdim(), Lam->get_store(), 1);
      }
    }

    ///computes S+=d*(*this);
    virtual void addmeto(CH_Matrix_Classes::Symmatrix& S, CH_Matrix_Classes::Real d = 1.) const {
      S(ii, jj) += d * val;
    }

    ///comutes B+=d*(*this)*C
    virtual void addprodto(CH_Matrix_Classes::Matrix& B, const CH_Matrix_Classes::Matrix& C, CH_Matrix_Classes::Real d = 1.) const {
      CH_Matrix_Classes::mat_xpeya(C.coldim(), B.get_store() + ii, B.rowdim(), C.get_store() + jj, C.rowdim(), d * val);
      if (ii == jj) return;
      CH_Matrix_Classes::mat_xpeya(C.coldim(), B.get_store() + jj, B.rowdim(), C.get_store() + ii, C.rowdim(), d * val);
    }

    ///computes B+=d*(*this)*C
    virtual void addprodto(CH_Matrix_Classes::Matrix& B, const CH_Matrix_Classes::Sparsemat& C, CH_Matrix_Classes::Real d = 1.) const {
      if (C.coldim() == 1) {
        B(ii) += d * val * C(jj);
        if (ii != jj) B(jj) += d * val * C(ii);
        return;
      }
      CH_Matrix_Classes::Sparsemat tmp(C.row(jj));
      for (CH_Matrix_Classes::Integer i = 0; i < tmp.nonzeros(); i++) {
        CH_Matrix_Classes::Integer indi, indj; CH_Matrix_Classes::Real v;
        tmp.get_edge(i, indi, indj, v);
        B(ii, indj) += d * val * v;
      }
      if (ii == jj) return;
      tmp.init(C.row(ii));
      {
        for (CH_Matrix_Classes::Integer i = 0; i < tmp.nonzeros(); i++) {
          CH_Matrix_Classes::Integer indi, indj; CH_Matrix_Classes::Real v;
          tmp.get_edge(i, indi, indj, v);
          B(jj, indj) += d * val * v;
        }
      }
    }

    /// computes R=P^T*(*this)*Q
    virtual void left_right_prod(const CH_Matrix_Classes::Matrix& P, const CH_Matrix_Classes::Matrix& Q, CH_Matrix_Classes::Matrix& R) const {
      if (ii == jj) CH_Matrix_Classes::genmult(P.row(ii), Q.row(ii), R, val, 0., 1, 0);
      else {
        CH_Matrix_Classes::genmult(P.row(ii), Q.row(jj), R, val, 0., 1, 0);
        CH_Matrix_Classes::genmult(P.row(jj), Q.row(ii), R, val, 1., 1, 0);
      }
    }

    ///returns an estimate of number of flops to compute addprodto for a vector
    virtual CH_Matrix_Classes::Integer prodvec_flops() const {
      return (ii == jj) ? 2 : 4;
    }

    ///returns 1 if its structure is as bad as its dense symmetric representation, otherwise 0
    virtual int dense() const {
      return 0;
    }

    ///returns 0 if not sparse, otherwise 1
    virtual int sparse() const {
      return 1;
    }

    /// returns 0 if not sparse. If it is sparse it returns 1 and the nonzero structure in I,J and v, where v is multiplied by d. Only the upper triangle (including diagonal) is delivered
    virtual int sparse(CH_Matrix_Classes::Indexmatrix& I, CH_Matrix_Classes::Indexmatrix& J, CH_Matrix_Classes::Matrix& v, CH_Matrix_Classes::Real d = 1.)const {
      I.init(1, 1, ii); J.init(1, 1, jj); v.init(1, 1, val * d); return 1;
    }

    /// returns 0 if the support of the costraint matrix is not contained in the support of the sparse symmetric matrix S, 1 if it is contained.
    virtual int support_in(const CH_Matrix_Classes::Sparsesym& S) const {
      return S.check_support(ii, jj);
    }

    ///returns the inner product of the constraint matrix with S
    virtual CH_Matrix_Classes::Real ip(const CH_Matrix_Classes::Sparsesym& S) const {
      if (ii == jj) return val * S(ii, jj); return 2. * val * S(ii, jj);
    }

    ///computes S=P^T*(*this)*P
    virtual void project(CH_Matrix_Classes::Symmatrix& S, const CH_Matrix_Classes::Matrix& P) const {
      S.newsize(P.coldim()); chk_set_init(S, 1);
      const CH_Matrix_Classes::Real* ap = P.get_store() + ii - nr;
      CH_Matrix_Classes::Real* sp = S.get_store();
      if (ii == jj) {
        for (CH_Matrix_Classes::Integer i = P.coldim(); --i >= 0;) {
          const CH_Matrix_Classes::Real* aap = ap; CH_Matrix_Classes::Real a = *(aap = (ap += nr));
          *sp++ = a * a * val; a *= val;
          for (CH_Matrix_Classes::Integer j = i; --j >= 0;) *sp++ = a * (*(aap += nr));
        }
        return;
      }
      const CH_Matrix_Classes::Real* bp = P.get_store() + jj - nr;
      for (CH_Matrix_Classes::Integer i = P.coldim(); --i >= 0;) {
        const CH_Matrix_Classes::Real* aap; CH_Matrix_Classes::Real a = *(aap = (ap += nr)) * val;
        const CH_Matrix_Classes::Real* bbp; CH_Matrix_Classes::Real b = *(bbp = (bp += nr));
        *sp++ = 2. * a * b; b *= val;
        for (CH_Matrix_Classes::Integer j = i; --j >= 0;) *sp++ = a * (*(bbp += nr)) + b * (*(aap += nr));
      }
    }

    ///computes S+=Q^T(*this)Q for Q=P.rows(start_row,start_row+dim-1) 
    virtual void add_projection(CH_Matrix_Classes::Symmatrix& S, const CH_Matrix_Classes::Matrix& P, CH_Matrix_Classes::Real alpha = 1., CH_Matrix_Classes::Integer start_row = 0) const {
      CH_Matrix_Classes::Integer pnr = P.rowdim();
      const CH_Matrix_Classes::Real* ap = P.get_store() + start_row + ii - pnr;
      CH_Matrix_Classes::Real* sp = S.get_store();
      const CH_Matrix_Classes::Real alval = alpha * val;
      if (ii == jj) {
        for (CH_Matrix_Classes::Integer i = P.coldim(); --i >= 0;) {
          const CH_Matrix_Classes::Real* aap = ap;
          CH_Matrix_Classes::Real a = *(aap = (ap += pnr));
          (*sp++) += a * a * alval; a *= alval;
          for (CH_Matrix_Classes::Integer j = i; --j >= 0;) (*sp++) += a * (*(aap += pnr));
        }
        return;
      }
      const CH_Matrix_Classes::Real* bp = P.get_store() + start_row + jj - pnr;
      for (CH_Matrix_Classes::Integer i = P.coldim(); --i >= 0;) {
        const CH_Matrix_Classes::Real* aap;
        CH_Matrix_Classes::Real a = *(aap = (ap += pnr)) * alval;
        const CH_Matrix_Classes::Real* bbp;
        CH_Matrix_Classes::Real b = *(bbp = (bp += pnr));
        (*sp++) += 2. * a * b; b *= alval;
        for (CH_Matrix_Classes::Integer j = i; --j >= 0;) (*sp++) += a * (*(bbp += pnr)) + b * (*(aap += pnr));
      }
    }

    ///computes C= alpha*(*this)*B^(T if btrans) + beta*C, C is also returned
    virtual const CH_Matrix_Classes::Matrix& postgenmult(const CH_Matrix_Classes::Matrix& B, CH_Matrix_Classes::Matrix& C,
      CH_Matrix_Classes::Real alpha = 1., CH_Matrix_Classes::Real beta = 0., int btrans = 0) const {
      CH_Matrix_Classes::Sparsesym S(nr, 1, &ii, &jj, &val);
      return CH_Matrix_Classes::genmult(S, B, C, alpha, beta, btrans);
    }

    ///computes C= alpha*B^(T if btrans)*(*this) + beta*C, C is also returned
    virtual const CH_Matrix_Classes::Matrix& pregenmult(const CH_Matrix_Classes::Matrix& B, CH_Matrix_Classes::Matrix& C,
      CH_Matrix_Classes::Real alpha = 1., CH_Matrix_Classes::Real beta = 0., int btrans = 0) const {
      CH_Matrix_Classes::Sparsesym S(nr, 1, &ii, &jj, &val);
      CH_Matrix_Classes::Matrix D;
      return CH_Matrix_Classes::genmult(B, S, C, alpha, beta, btrans);
    }

    ///returns 1, if p is the same derived class and entries differ by less than tol, otherwise zero
    virtual int equal(const Coeffmat* p, double tol = 1e-6) const {
      const CMsingleton* pp = dynamic_cast<const CMsingleton*>(p);
      if (pp == 0)
        return 0;
      return ((nr == pp->nr) && (ii == pp->ii) && (jj == pp->jj)
        && (fabs(val - pp->val) < tol));
    }

    ///display constraint information
    virtual std::ostream& display(std::ostream& o) const {
      o << "CMsingleton\n"; o << nr << " " << ii << " " << jj << " " << val << "\n"; return o;
    }


    ///put entire contents onto outstream with the class type in the beginning so that the derived class can be recognized by in().
    virtual std::ostream& out(std::ostream& o) const {
      return o << "SINGLETON\n" << nr << " " << ii << " " << jj << " " << val << "\n";
    }
    //put entire contents onto outstream with the class type in the beginning so
    //that the derived class can be recognized.

  ///counterpart to out(), does not read the class type, though. This is assumed to have been read in order to generate the correct class
    virtual std::istream& in(std::istream& is) {
      if (!(is >> nr >> ii >> jj >> val)) {
        if (CH_Matrix_Classes::materrout) (*CH_Matrix_Classes::materrout) << "*** ERROR: CMsingleton::in(): reading from input failed";
        is.clear(std::ios::failbit);
        return is;
      }
      if (nr < 0) {
        if (CH_Matrix_Classes::materrout) (*CH_Matrix_Classes::materrout) << "*** ERROR: CMsingleton::in(): dimension of matrix must positive";
        if (CH_Matrix_Classes::materrout) (*CH_Matrix_Classes::materrout) << "but is " << nr << std::endl;
        is.clear(std::ios::failbit);
        return is;
      }
      if ((ii < 0) || (ii > nr)) {
        if (CH_Matrix_Classes::materrout) (*CH_Matrix_Classes::materrout) << "*** ERROR: CMsingleton::in(): row index outside range, ";
        if (CH_Matrix_Classes::materrout) (*CH_Matrix_Classes::materrout) << 0 << "<=" << ii << "<" << nr << std::endl;
        is.clear(std::ios::failbit);
        return is;
      }
      if ((jj < 0) || (jj > nr)) {
        if (CH_Matrix_Classes::materrout) (*CH_Matrix_Classes::materrout) << "*** ERROR: CMsingleton::in(): column index outside range, ";
        if (CH_Matrix_Classes::materrout) (*CH_Matrix_Classes::materrout) << 0 << "<=" << jj << "<" << nr << std::endl;
        is.clear(std::ios::failbit);
        return is;
      }
      return is;
    }

    /// constructor with instream and possibly additional user information
    CMsingleton(std::istream& is, CoeffmatInfo* cip = 0) {
      CM_type = CM_singleton; infop = cip; in(is);
    }

    //--- specific routines
    /// return the nonzero entry information
    int get_ijval(CH_Matrix_Classes::Integer& i, CH_Matrix_Classes::Integer& j, CH_Matrix_Classes::Real& v) const {
      i = ii; j = jj; v = val; return 0;
    }

  };

  //@}
}

#endif

