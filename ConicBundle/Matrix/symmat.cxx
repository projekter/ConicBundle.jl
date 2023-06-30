/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  Matrix/symmat.cxx
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



#include <stdlib.h>
#include "symmat.hxx"


using namespace CH_Tools;

namespace CH_Matrix_Classes {

  const Mtype Symmatrix::mtype = MTsymmetric;

  // **************************************************************************
  //                                Constructors
  // **************************************************************************

  Symmatrix& Symmatrix::xeya(const Symmatrix& A, Real d) {
    chk_init(A);
    newsize(A.nr);
    chk_set_init(*this, 1);
    if (d == 1.) {
      mat_xey((nr * (nr + 1)) / 2, m, A.get_store()); return *this;
    }
    if (d == 0.) {
      mat_xea((nr * (nr + 1)) / 2, m, 0.); return *this;
    }
    if (d == -1.) {
      mat_xemy((nr * (nr + 1)) / 2, m, A.get_store()); return *this;
    }
    mat_xeya((nr * (nr + 1)) / 2, m, A.get_store(), d);
    return *this;
  }

  Symmatrix& Symmatrix::xpeya(const Symmatrix& A, Real d) {
    chk_add(*this, A);
    if (d == 1.) {
      mat_xpey((nr * (nr + 1)) / 2, m, A.get_store()); return *this;
    }
    if (d == 0.) {
      return *this;
    }
    if (d == -1.) {
      mat_xmey((nr * (nr + 1)) / 2, m, A.get_store()); return *this;
    }
    mat_xpeya((nr * (nr + 1)) / 2, m, A.get_store(), d);
    return *this;
  }


  // returns C=beta*C+alpha* A*A^T, where A may be transposed
  Symmatrix& rankadd(const Matrix& A, Symmatrix& C, Real alpha, Real beta, int trans)
#ifdef WITH_BLAS
  {
    chk_init(A);
    Integer nr = (trans ? A.coldim() : A.rowdim());
    Integer nk = (trans ? A.rowdim() : A.coldim());
#if (CONICBUNDLE_DEBUG>=1)
    if (beta != 0.) {
      chk_init(C);
      if (C.nr != nr) {
        MEmessage(MatrixError(ME_dim, "rankadd: dimensions don't match", MTmatrix));;
      }
    }
#endif
    if ((nr == 0) || (nk == 0)) {
      if (beta == 0.) return C.init(nr, 0.);
      return C *= beta;
    }
    if (beta == 0.) {
      C.newsize(nr);
      chk_set_init(C, 1);
    }
    Matrix tmp(nr, nr); chk_set_init(tmp, 1);
    Real* tp = tmp.get_store();
    Real* cp = C.get_store();
    for (Integer i = nr; i > 0; --i) {
      mat_xey(i, tp, cp);
      cp += i;
      tp += nr + 1;
    }
    cblas_dsyrk(CblasColMajor, CblasLower, (trans ? CblasTrans : CblasNoTrans),
      nr, nk, alpha, A.get_store(), A.rowdim(), beta, tmp.get_store(), tmp.rowdim());
    tp = tmp.get_store();
    cp = C.get_store();
    for (Integer i = nr; i > 0; --i) {
      mat_xey(i, cp, tp);
      cp += i;
      tp += nr + 1;
    }
    return C;
  }
#else
  {
    chk_init(A);
    if (trans) { //A^T*A
      if (beta == 0.) {
        C.init(A.coldim(), 0.);
      } else {
        chk_colseq(C, A);
        if (beta != 1.) C *= beta;
      }
      if (alpha == 0.) return C;
      const Integer nr = C.rowdim();
      const Integer n = A.rowdim();
#ifdef WITH_OMP
      const Integer blocksz = 300;
      if (nr < blocksz) {
#endif
        Real* mp = C.get_store();
        const Real* ap = A.get_store() - n;
        for (Integer i = nr; --i >= 0;) {
          const Real* aap = (ap += n);
          (*mp++) += alpha * mat_ip(n, ap, aap);
          for (Integer j = i; --j >= 0;) {
            (*mp++) += alpha * mat_ip(n, ap, aap += n);
          }
        }
#ifdef WITH_OMP
      } else {
        Real* const cp = C.get_store();
        const Real* const ap = A.get_store();
#pragma omp parallel for schedule(dynamic,1)
        for (Integer i = 0; i < nr; i++) {
          Real* mp = cp + i * nr - (i * (i - 1)) / 2;
          const Real* const aip = ap + i * n;
          const Real* ajp = aip;
          (*mp++) += alpha * mat_ip(n, aip, ajp);
          for (Integer j = nr - i; --j > 0;) {
            (*mp++) += alpha * mat_ip(n, aip, ajp += n);
          }
        }
      }
#endif       
      return C;
    } else { //A*A^T
      if (beta == 0.) {
        C.init(A.rowdim(), 0.);
      } else {
        chk_rowseq(C, A);
        if (beta != 1.) C *= beta;
      }
      if (alpha == 0.) return C;
      const Integer nr = C.rowdim();
      Integer nc2 = (nr * (nr + 1)) / 2;
      Integer n = A.coldim();

      const Real* ap = A.get_store() - 1;
      Real* mp = C.get_store();
      for (Integer i = n; --i >= 0;) {
        for (Integer j = nr; --j >= 0; ) {
          const Real* aap;
          const Real d1 = *(aap = ++ap);
          if (d1 != 0.) {
            const Real d = d1 * alpha;
            (*mp++) += d * d1;
            for (const Real* endp = aap + j; ++aap <= endp;)
              (*mp++) += d * (*aap);
          } else {
            mp += j + 1;
          }
        }
        mp -= nc2;
      }
      /*
      for(Integer i=0;i<n;i++){
        const Integer r=i*nr;
 #pragma omp parallel for schedule(static,10)
        for(Integer j=0;j<nr;j++){
    const Real d=A.get_store()[r+j]*alpha;
    mat_xpeya(nr-j,C.get_store()+j*nr-(j*(j-1))/2,A.get_store()+r+j,d);
        }
      }
      */
      return C;
    }
  }
#endif

  // returns C=beta*C+alpha* A*A^T, where A may be transposed
  Symmatrix& scaledrankadd(const Matrix& A, const Matrix& D, Symmatrix& C, Real alpha, Real beta, int trans) {
    chk_init(A);
    if (trans) { //A^T*A
      assert(D.dim() == A.rowdim());
      if (beta == 0.) {
        C.init(A.coldim(), 0.);
      } else {
        chk_colseq(C, A);
        if (beta != 1.) C *= beta;
      }
      if (alpha == 0.) return C;
      const Integer nr = C.rowdim();
      const Integer n = A.rowdim();
      Real* mp = C.get_store();
      Matrix vec(n, 1); chk_set_init(vec, 1);
      const Real* const vend = vec.get_store() + n;
      const Real* ap = A.get_store() - n;
      for (Integer i = nr; --i >= 0;) {
        const Real* aap = (ap += n);
        //scale the next "row"-vector
        const Real* ddp = D.get_store();
        const Real* aip = ap;
        for (Real* vp = vec.get_store(); vp != vend;)
          *vp++ = (*aip++) * (*ddp++);
        //take the inner products of the remaing vectors with this
        for (Integer j = i + 1; --j >= 0;) {
          Real sumval = 0.;
          for (const Real* vp = vec.get_store(); vp != vend;)
            sumval += *aap++ * (*vp++);
          (*mp++) += alpha * sumval;
        }
      }
    } else { //A*A^T
      assert(A.coldim() == D.dim());
      if (beta == 0.) {
        C.init(A.rowdim(), 0.);
      } else {
        chk_rowseq(C, A);
        if (beta != 1.) C *= beta;
      }
      if (alpha == 0.) return C;
      const Integer nr = C.rowdim();
      Integer nc2 = (nr * (nr + 1)) / 2;
      Integer n = A.coldim();

      const Real* ap = A.get_store();
      Real* mp = C.get_store();
      const Real* dp = D.get_store();
      for (Integer i = n; --i >= 0;) {
        Real da = (*dp++) * alpha;
        if (da != 0.) {
          const Real* endp = ap + nr;
          for (Integer j = nr; --j >= 0; ) {
            const Real* aap = ap++;
            const Real d1 = *aap;
            if (d1 != 0.) {
              const Real d = d1 * da;
              (*mp++) += d * d1;
              for (; ++aap != endp;)
                (*mp++) += d * (*aap);
            } else {
              mp += j + 1;
            }
          }
          mp -= nc2;
        } else {
          ap += nr;
        }
      }
    }
    return C;
  }

  // returns C=beta*C+alpha* sym(A*B^T) [or sym(A^T*B)]
  Symmatrix& rank2add(const Matrix& A, const Matrix& B, Symmatrix& C,
    Real alpha, Real beta, int trans)
#ifdef WITH_BLAS
  {
    chk_add(A, B);
    Integer nr = (trans ? A.coldim() : A.rowdim());
    Integer nk = (trans ? A.rowdim() : A.coldim());
#if (CONICBUNDLE_DEBUG>=1)
    if (beta != 0.) {
      chk_init(C);
      if (C.nr != nr) {
        MEmessage(MatrixError(ME_dim, "rank2add: dimensions don't match", MTmatrix));;
      }
    }
#endif
    if ((nr == 0) || (nk == 0)) {
      if (beta == 0.) return C.init(nr, 0.);
      return C *= beta;
    }
    if (beta == 0.) {
      C.newsize(nr);
      chk_set_init(C, 1);
    }
    Matrix tmp(nr, nr); chk_set_init(tmp, 1);
    Real* tp = tmp.get_store();
    Real* cp = C.get_store();
    for (Integer i = nr; i > 0; --i) {
      mat_xey(i, tp, cp);
      cp += i;
      tp += nr + 1;
    }
    cblas_dsyr2k(CblasColMajor, CblasLower, (trans ? CblasTrans : CblasNoTrans),
      nr, nk, alpha / 2., A.get_store(), A.rowdim(), B.get_store(), B.rowdim(), beta,
      tmp.get_store(), tmp.rowdim());
    tp = tmp.get_store();
    cp = C.get_store();
    for (Integer i = nr; i > 0; --i) {
      mat_xey(i, cp, tp);
      cp += i;
      tp += nr + 1;
    }
    return C;
  }
#else
  {
    chk_add(A, B);
    if (trans) { //A^T*B
      if (beta == 0.) {
        C.init(A.coldim(), 0.);
      } else {
        chk_colseq(C, A);
        if (beta != 1.) C *= beta;
      }
      if (alpha == 0.) return C;
      Integer nr = C.rowdim();
      Integer n = A.rowdim();
      Real* mp = C.get_store();
      const Real* abp = A.get_store() - n;
      const Real* bbp = B.get_store() - n;
      Real d = alpha / 2.;
      for (Integer i = nr; --i >= 0;) {
        const Real* bp;
        const Real* ap;
        *mp++ += alpha * mat_ip(n, ap = (abp += n), bp = (bbp += n));
        for (Integer j = i; --j >= 0;) {
          (*mp++) += d * (mat_ip(n, abp, bp += n) + mat_ip(n, ap += n, bbp));
        }
      }
      return C;
    } else {     //A*B^T
      if (beta == 0.) {
        C.init(A.rowdim(), 0.);
      } else {
        chk_rowseq(C, A);
        if (beta != 1.) C *= beta;
      }
      if (alpha == 0.) return C;
      Integer nr = C.rowdim();
      Integer nc2 = (nr * (nr + 1)) / 2;
      Integer n = A.coldim();
      const Real* abp = A.get_store() - 1;
      const Real* bbp = B.get_store() - 1;
      Real* mp = C.get_store();
      for (Integer i = n; --i >= 0;) {
        for (Integer j = nr; --j >= 0; ) {
          const Real* ap;
          const Real* bp;
          Real ad = *(ap = ++abp);
          Real bd = alpha * *(bp = ++bbp);
          (*mp++) += ad * bd;
          ad *= alpha / 2.;
          bd /= 2.;
          for (Integer k = j; --k >= 0;)
            (*mp++) += bd * (*(++ap)) + ad * (*(++bp));
        }
        mp -= nc2;
      }
      return C;
    }
  }
#endif


  Symmatrix& Symmatrix::xetriu_yza(const Matrix& A, const Matrix& B, Real d) {
    chk_add(A, B);
    newsize(A.coldim());
    chk_set_init(*this, 1);
    if (d == 0.) {
      mat_xea((nr * (nr + 1)) / 2, m, 0.); return *this;
    }
    if (d == 1.) {
      Real* mp = m;
      const Real* ap = A.get_store();
      Integer n = A.rowdim();
      const Real* bbp = B.get_store() - n;
      for (Integer i = nr; i > 0; i--) {
        const Real* bp = (bbp += n);
        for (Integer j = i; --j >= 0;) {
          *mp++ = mat_ip(n, ap, bp);
          bp += n;
        }
        ap += n;
      }
      return *this;
    }
    if (d == -1.) {
      Real* mp = m;
      const Real* ap = A.get_store();
      Integer n = A.rowdim();
      const Real* bbp = B.get_store() - n;
      for (Integer i = nr; i > 0; i--) {
        const Real* bp = (bbp += n);
        for (Integer j = i; --j >= 0;) {
          *mp++ = -mat_ip(n, ap, bp);
          bp += n;
        }
        ap += n;
      }
      return *this;
    }
    Real* mp = m;
    const Real* ap = A.get_store();
    Integer n = A.rowdim();
    const Real* bbp = B.get_store() - n;
    for (Integer i = nr; i > 0; i--) {
      const Real* bp = (bbp += n);
      for (Integer j = i; --j >= 0;) {
        *mp++ = d * mat_ip(n, ap, bp);
        bp += n;
      }
      ap += n;
    }
    return *this;
  }

  Symmatrix& Symmatrix::xpetriu_yza(const Matrix& A, const Matrix& B, Real d) {
    chk_add(A, B);
    chk_colseq(*this, A);
    if (d == 0.) {
      return *this;
    }
    if (d == 1.) {
      Real* mp = m;
      const Real* ap = A.get_store();
      Integer n = A.rowdim();
      for (Integer i = 0; i < nr; i++) {
        const Real* bp = B.get_store() + i * n;
        for (Integer j = i; j < nr; j++) {
          (*mp++) += mat_ip(n, ap, bp);
          bp += n;
        }
        ap += n;
      }
      return *this;
    }
    if (d == -1.) {
      Real* mp = m;
      const Real* ap = A.get_store();
      Integer n = A.rowdim();
      for (Integer i = 0; i < nr; i++) {
        const Real* bp = B.get_store() + i * n;
        for (Integer j = i; j < nr; j++) {
          (*mp++) -= mat_ip(n, ap, bp);
          bp += n;
        }
        ap += n;
      }
      return *this;
    }
    Real* mp = m;
    const Real* ap = A.get_store();
    Integer n = A.rowdim();
    for (Integer i = 0; i < nr; i++) {
      const Real* bp = B.get_store() + i * n;
      for (Integer j = i; j < nr; j++) {
        (*mp++) += d * mat_ip(n, ap, bp);
        bp += n;
      }
      ap += n;
    }
    return *this;
  }

  Matrix& genmult(const Symmatrix& A, const Matrix& B, Matrix& C,
    Real alpha, Real beta, int btrans)
    //returns C=beta*C+alpha*A*B, where A and B may be transposed
#ifdef WITH_BLAS
  {
    chk_init(A);
    chk_init(B);
    Integer nr = A.rowdim();
    Integer nc = (btrans ? B.rowdim() : B.coldim());
#if (CONICBUNDLE_DEBUG>=1)
    if (beta != 0.) chk_init(C);
    if ((nr != (btrans ? B.coldim() : B.rowdim())) ||
      ((beta != 0.) && ((C.rowdim() != nr) || (C.coldim() != nc)))
      ) {
      MEmessage(MatrixError(ME_dim, "genmult(Symmatrix,Matrix,Matrix): dimensions don't match", MTmatrix));;
    }
#endif
    if ((nr == 0) || (nc == 0)) {
      if (beta == 0.) return C.init(nr, nc, 0.);
      return C *= beta;
    }
    if (beta == 0.) {
      C.newsize(nr, nc);
      chk_set_init(C, 1);
    }
    Matrix tmp(nr, nr); chk_set_init(tmp, 1);
    Real* tp = tmp.get_store();
    const Real* ap = A.get_store();
    for (Integer i = nr; i > 0; --i) {
      mat_xey(i, tp, ap);
      ap += i;
      tp += nr + 1;
    }
    if (btrans == 0) {
      cblas_dsymm(CblasColMajor, CblasLeft, CblasLower, nr, nc, alpha,
        tmp.get_store(), tmp.rowdim(), B.get_store(), B.rowdim(), beta,
        C.get_store(), C.rowdim());
    } else {
      const Integer br = B.rowdim();
      const Integer bc = B.coldim();
      Matrix tmpB(bc, br); chk_set_init(tmpB, 1);
      Real* tp = tmpB.get_store();
      const Real* bp = B.get_store();
      for (Integer i = 0; i < br; i++) {
        mat_xey(br, tp, bc, bp, 1);
        bp += br;
        tp++;
      }
      cblas_dsymm(CblasColMajor, CblasLeft, CblasLower, nr, nc, alpha,
        tmp.get_store(), tmp.rowdim(), tmpB.get_store(), tmpB.rowdim(), beta,
        C.get_store(), C.rowdim());
    }
    return C;
  }
#else
  {
    chk_init(A);
    chk_init(B);
    Integer nr, nc;
    nr = A.nr;
#if (CONICBUNDLE_DEBUG>=1)
    Integer nm = A.nr;
#endif
    if (btrans) {
      nc = B.nr;
#if (CONICBUNDLE_DEBUG>=1)
      if (nm != B.nc) {
        MEmessage(MatrixError(ME_dim, "genmult(Symmatrix,Matrix,Matrix): dimensions don't match", MTsymmetric));;
      }
#endif
    } else {
      nc = B.nc;
#if (CONICBUNDLE_DEBUG>=1)
      if (nm != B.nr) {
        MEmessage(MatrixError(ME_dim, "genmult(Symmatrix,Matrix,Matrix): dimensions don't match", MTsymmetric));;
      }
#endif
    }
    if (beta != 0.) {
      chk_init(C);
#if (CONICBUNDLE_DEBUG>=1)
      if ((nr != C.nr) || (nc != C.nc)) {
        MEmessage(MatrixError(ME_dim, "genmult(Symmatrix,Matrix,Matrix): dimensions don't match", MTsymmetric));;
      }
#endif
      if (beta != 1.) C *= beta;
    } else {
      C.init(nr, nc, 0.);
    }
    if (alpha == 0.) return C;
    if ((btrans) && (nc > 1)) {
      /*
        const Real *ap=A.get_store();
        for(Integer j=0;j<nm;j++){
            mat_xpeya(nc,C.get_store()+j,nr,B.get_store()+j*B.nc,1,alpha*(*ap++));
            for(Integer i=j+1;i<nr;i++){
                Real d=alpha*(*ap++);
                mat_xpeya(nc,C.get_store()+i,nr,B.get_store()+j*B.nc,1,d);
                mat_xpeya(nc,C.get_store()+j,nr,B.get_store()+i*B.nc,1,d);
            }
        }
      */
      Real* cp = C.m;
      for (Integer col = 0; col < nc; col++) {
        const Real* bp = B.m + col;
        const Real* ap = A.m; //the symmetric store
        for (Integer i = nr; --i >= 0;) {
          Real bi = (*bp);
          Real dsum = alpha * bi * (*ap++);
          Real* cjp = cp + 1;
          for (Integer j = i; --j >= 0;) {
            Real aij = alpha * (*ap++);
            (*cjp++) += bi * aij;
            bp += nc;
            dsum += aij * (*bp);
          }
          (*cp++) += dsum;
          bp -= (i - 1) * nc;
        }
      }
    } else {
      /*
        const Real *ap=A.get_store();
        for(Integer j=0;j<nm;j++){
            mat_xpeya(nc,C.get_store()+j,nr,B.get_store()+j,B.nr,alpha*(*ap++));
            for(Integer i=j+1;i<nr;i++){
                Real d=alpha*(*ap++);
                mat_xpeya(nc,C.get_store()+i,nr,B.get_store()+j,B.nr,d);
                mat_xpeya(nc,C.get_store()+j,nr,B.get_store()+i,B.nr,d);
            }
        }
      */
      const Real* bp = B.m;
      Real* cp = C.m;
      for (Integer col = nc; --col >= 0;) {
        const Real* ap = A.m; //the symmetric store
        for (Integer i = nr; i > 0;) {
          Real d = alpha * (*bp);
          (*cp++) += alpha * mat_ip(i--, bp++, ap++);
          if (d != 0.) {
            mat_xpeya(i, cp, ap, d);
          }
          ap += i;
        }
      }

    }
    return C;
  }
#endif

  Matrix& genmult(const Matrix& A, const Symmatrix& B, Matrix& C,
    Real alpha, Real beta, int atrans)
    //returns C=beta*C+alpha*A*B, where A and B may be transposed
#ifdef WITH_BLAS
  {
    chk_init(B);
    chk_init(A);
    Integer nr = (atrans ? A.coldim() : A.rowdim());
    Integer nc = B.rowdim();
#if (CONICBUNDLE_DEBUG>=1)
    if (beta != 0.) chk_init(C);
    if ((nc != (atrans ? A.rowdim() : A.coldim())) ||
      ((beta != 0.) && ((C.rowdim() != nr) || (C.coldim() != nc)))
      ) {
      MEmessage(MatrixError(ME_dim, "genmult(Symmatrix,Matrix,Matrix): dimensions don't match", MTmatrix));;
    }
#endif
    if ((nr == 0) || (nc == 0)) {
      if (beta == 0.) return C.init(nr, nc, 0.);
      return C *= beta;
    }
    if (beta == 0.) {
      C.newsize(nr, nc);
      chk_set_init(C, 1);
    }
    Matrix tmp(nc, nc); chk_set_init(tmp, 1);
    Real* tp = tmp.get_store();
    const Real* bp = B.get_store();
    for (Integer i = nc; i > 0; --i) {
      mat_xey(i, tp, bp);
      bp += i;
      tp += nc + 1;
    }
    if (atrans == 0) {
      cblas_dsymm(CblasColMajor, CblasRight, CblasLower, nr, nc, alpha,
        tmp.get_store(), tmp.rowdim(), A.get_store(), A.rowdim(), beta,
        C.get_store(), C.rowdim());
    } else {
      const Integer ar = A.rowdim();
      const Integer ac = A.coldim();
      Matrix tmpA(ac, ar); chk_set_init(tmpA, 1);
      Real* tp = tmpA.get_store();
      const Real* ap = A.get_store();
      for (Integer i = 0; i < ar; i++) {
        mat_xey(ar, tp, ac, ap, 1);
        ap += ar;
        tp++;
      }
      cblas_dsymm(CblasColMajor, CblasRight, CblasLower, nr, nc, alpha,
        tmp.get_store(), tmp.rowdim(), tmpA.get_store(), tmpA.rowdim(), beta,
        C.get_store(), C.rowdim());
    }
    return C;
  }
#else
  {
    chk_init(A);
    chk_init(B);
    Integer nr, nm, nc;
    if (atrans) {
      nr = A.nc;
      nm = A.nr;
    } else {
      nr = A.nr;
      nm = A.nc;
    }
    nc = B.nr;
#if (CONICBUNDLE_DEBUG>=1)
    if (nm != B.nr) {
      MEmessage(MatrixError(ME_dim, "genmult(Symmatrix,Matrix,Matrix): dimensions don't match", MTsymmetric));;
    }
#endif
    if (beta != 0.) {
      chk_init(C);
#if (CONICBUNDLE_DEBUG>=1)
      if ((nr != C.nr) || (nc != C.nc)) {
        MEmessage(MatrixError(ME_dim, "genmult(Symmatrix,Matrix,Matrix): dimensions don't match", MTsymmetric));;
      }
#endif
      if (beta != 1.) C *= beta;
    } else {
      C.init(nr, nc, 0.);
    }
    if (alpha == 0.) return C;
    if (atrans) {
      const Real* bp = B.get_store();
      for (Integer j = 0; j < nc; j++) {
        mat_xpeya(nr, C.get_store() + j * nr, 1, A.get_store() + j, nm, alpha * (*bp++));
        for (Integer i = j + 1; i < nm; i++) {
          Real d = alpha * (*bp++);
          mat_xpeya(nr, C.get_store() + j * nr, 1, A.get_store() + i, nm, d);
          mat_xpeya(nr, C.get_store() + i * nr, 1, A.get_store() + j, nm, d);
        }
      }
    } else {
      const Real* bp = B.get_store();
      for (Integer j = 0; j < nc; j++) {
        mat_xpeya(nr, C.get_store() + j * nr, A.get_store() + j * nr, alpha * (*bp++));
        for (Integer i = j + 1; i < nm; i++) {
          Real d = alpha * (*bp++);
          mat_xpeya(nr, C.get_store() + j * nr, A.get_store() + i * nr, d);
          mat_xpeya(nr, C.get_store() + i * nr, A.get_store() + j * nr, d);
        }
      }
    }
    return C;
  }
#endif

  Symmatrix& Symmatrix::xeya(const Matrix& A, double d) {
    chk_init(A);
#if (CONICBUNDLE_DEBUG>=1)
    if (A.nr != A.nc)
      MEmessage(MEdim(A.nr, A.nc, 0, 0, "Symmatrix::Symmatrix(const Matrix& A) A not square", MTsymmetric));
    is_init = 0;
#endif
    newsize(A.nr);
    Real* mp = m;
    Real f = d / 2.;
    for (Integer i = 0; i < nr; i++) {
      Real* matcp = A.m + i * nr + i;
      Real* matrp = matcp + nr;
      (*mp++) = (*matcp++) * d;
      for (Integer j = i + 1; j < nr; j++, matrp += nr)
        (*mp++) = ((*matcp++) + (*matrp)) * f;
    }
    chk_set_init(*this, 1);
    return *this;
  }

  Symmatrix& Symmatrix::xpeya(const Matrix& A, double d) {
    chk_add(*this, A);
    Real* mp = m;
    Real f = d / 2.;
    for (Integer i = 0; i < nr; i++) {
      Real* matcp = A.m + i * nr + i;
      Real* matrp = matcp + nr;
      (*mp++) += (*matcp++) * d;
      for (Integer j = i + 1; j < nr; j++, matrp += nr)
        (*mp++) += ((*matcp++) + (*matrp)) * f;
    }
    return *this;
  }

  Symmatrix& Symmatrix::xeya(const Indexmatrix& A, double d) {
    chk_init(A);
#if (CONICBUNDLE_DEBUG>=1)
    if (A.nr != A.nc)
      MEmessage(MEdim(A.nr, A.nc, 0, 0, "Symmatrix::Symmatrix(const Matrix& A) A not square", MTsymmetric));
    is_init = 0;
#endif
    newsize(A.nr);
    Real* mp = m;
    Real f = d / 2.;
    for (Integer i = 0; i < nr; i++) {
      Integer* matcp = A.m + i * nr + i;
      Integer* matrp = matcp + nr;
      (*mp++) = (*matcp++) * d;
      for (Integer j = i + 1; j < nr; j++, matrp += nr)
        (*mp++) = ((*matcp++) + (*matrp)) * f;
    }
    chk_set_init(*this, 1);
    return *this;
  }

  Symmatrix& Symmatrix::xpeya(const Indexmatrix& A, double d) {
    chk_add(*this, A);
    Real* mp = m;
    Real f = d / 2.;
    for (Integer i = 0; i < nr; i++) {
      Integer* matcp = A.m + i * nr + i;
      Integer* matrp = matcp + nr;
      (*mp++) += (*matcp++) * d;
      for (Integer j = i + 1; j < nr; j++, matrp += nr)
        (*mp++) += ((*matcp++) + (*matrp)) * f;
    }
    return *this;
  }

  void Symmatrix::newsize(Integer inr) {
#if (CONICBUNDLE_DEBUG>=1)
    if (inr < 0)
      MEmessage(MEdim(inr, inr, 0, 0, "Symmatrix::newsize(Integer,Integer) dim<0", MTsymmetric));
    is_init = 0;
#endif
    if (inr != nr) {
      nr = inr;
      Integer n2 = (nr * (nr + 1)) / 2;
      if (n2 > mem_dim) {
        memarray->free(m); m = 0;
        mem_dim = Integer(memarray->get(n2, m));
        if (mem_dim < n2)
          MEmessage(MEmem(n2,
            "Symmatrix::Symmatrix(Integer,Integer,Real) not enough memory", MTsymmetric));
      }
    }
  }

  Symmatrix& Symmatrix::shift_diag(Real s) {
    chk_init(*this);
    Real* mp = m;
    for (Integer i = 0; i < nr; i++) {
      *mp += s;
      mp += nr - i;
    }
    return *this;
  }

  void Symmatrix::display(std::ostream& out, int precision, int width,
    int screenwidth) const {
    chk_init(*this);
    out << "Symmatrix(" << nr << ")" << std::endl;
    if (nr == 0) return;
    if (precision == 0) precision = 4;
    out.precision(precision);
    if (width == 0) width = precision + 6;
    if (screenwidth == 0) screenwidth = 80;
    Integer colnr = screenwidth / (width + 1);
    Integer k, i, j;
    Integer maxk = nr / colnr + ((nr % colnr) > 0);
    Integer maxj;
    for (k = 0; k < maxk; k++) {
      out << "columns " << k * colnr << " to " << min(nr, (k + 1) * colnr) - 1 << std::endl;
      for (i = 0; i < nr; i++) {
        maxj = min((k + 1) * colnr, nr);
        for (j = k * colnr; j < maxj; j++) {
          out << ' '; out.width(width); out << (*this)(i, j);
        }
        out << std::endl;
      }
    }
  }

  void Symmatrix::mfile_output(std::ostream& out, ///< output stream
    int precision,   ///< number of most significant digits, default=16
    int width       ///< field width, default = precision+6
  ) const {
    chk_init(*this);
    if (precision <= 0)
      precision = 16;
    std::streamsize old_prec = out.precision();
    if (width <= 0)
      width = precision + 6;
    out << "[";
    out.precision(precision);
    for (Integer i = 0; i < nr; i++) {
      for (Integer j = 0; j < nr; j++) {
        out << " ";
        out.width(width);
        out << (*this)(i, j);
      }
      if (i < nr - 1)
        out << "\n";
    }
    out << "];\n";
    out.precision(old_prec);
  }



  Matrix Symmatrix::col(Integer c) const {
    chk_init(*this);
#if (CONICBUNDLE_DEBUG>=1)
    if ((c < 0) || (c >= nr))
      MEmessage(MErange(nr, 0, nr, c, "Symmatrix::col(Integer) index out of range", MTsymmetric));
#endif
    Matrix v(nr, 1);
    Real* mp = m + c;
    Real* vp = v.m;
    Integer i = 0;
    for (; i < c; i++) {
      (*vp++) = *mp;
      mp += nr - i - 1;
    }
    for (; i < nr; i++) (*vp++) = *mp++;
    chk_set_init(v, 1);
    return v;
  }

  Matrix Symmatrix::row(Integer r) const {
    chk_init(*this);
#if (CONICBUNDLE_DEBUG>=1)
    if ((r < 0) || (r >= nr))
      MEmessage(MErange(nr, r, nr, 0, "Symmatrix::row(Integer) index out of range", MTsymmetric));
#endif
    Matrix v(1, nr);
    Real* mp = m + r;
    Real* vp = v.m;
    Integer i = 0;
    for (; i < r; i++) {
      (*vp++) = *mp;
      mp += nr - i - 1;
    }
    for (; i < nr; i++) (*vp++) = *mp++;
    chk_set_init(v, 1);
    return v;
  }

  Matrix Symmatrix::cols(const Indexmatrix& vec) const  //returns cols as indexed by vec
  {
    chk_init(*this);
    Matrix A(nr, vec.dim());
    for (Integer j = 0; j < vec.dim(); j++) {
      Integer c = vec(j);
      assert((0 <= c) && (c < nr));
      Real* mp = m + c;
      Real* vp = A.m + j * nr;
      Integer i = 0;
      for (; i < c; i++) {
        (*vp++) = *mp;
        mp += nr - i - 1;
      }
      for (; i < nr; i++) (*vp++) = *mp++;
    }
    chk_set_init(A, 1);
    return A;
  }

  Matrix Symmatrix::rows(const Indexmatrix& vec) const  //returns rows as indexed by vec
  {
    chk_init(*this);
    Integer rdim = vec.dim();
    Matrix A(rdim, nr);
    for (Integer j = 0; j < rdim; j++) {
      Integer r = vec(j);
      assert((0 <= r) && (r < nr));
      Real* mp = m + r;
      Real* vp = A.m + j;
      Integer i = 0;
      for (; i < r; i++, vp += rdim) {
        *vp = *mp;
        mp += nr - i - 1;
      }
      for (; i < nr; i++, vp += rdim)
        *vp = *mp++;
    }
    chk_set_init(A, 1);
    return A;
  }

  Symmatrix& Symmatrix::swapij(Integer i, Integer j) {
    assert((0 <= i) && (i < nr) && (0 <= j) && (j < nr));
    if (j == i)
      return *this;
    if (j < i) {
      Integer h = i; i = j; j = h;
    }
    //swap the diagonal elements first
    Real* mp1 = m + i * nr - (i * (i - 1)) / 2;
    Real* mp2 = m + j * nr - (j * (j - 1)) / 2;
    Real d = *mp1; *mp1 = *mp2; *mp2 = d;
    mp1 = m + i;
    mp2 = m + j;
    //rows 0 to i-1
    for (Integer k = 0; k < i;) {
      d = *mp1; *mp1 = *mp2; *mp2 = d;
      k++;
      mp1 += nr - k;
      mp2 += nr - k;
    }
    //skip exchange of (i,i) and (i,j)
    mp1++; mp2 += nr - i - 1;
    //rows i+1 to j-1
    for (Integer k = i + 1; k < j;) {
      d = *mp1; *mp1++ = *mp2; *mp2 = d;
      k++;
      mp2 += nr - k;
    }
    //skip exchange of (i,j) and (j,j)
    mp1++; mp2++;
    //rows j+1 to nr-1
    for (Integer k = j + 1; k < nr; k++) {
      d = *mp1; *mp1++ = *mp2; *mp2++ = d;
    }
    return *this;
  }


  Symmatrix& Symmatrix::pivot_permute(const Indexmatrix& piv, bool inverse) {
    chk_init(piv);
    chk_init(*this);
    assert(nr == piv.dim());

    if (!inverse) {
      for (Integer i = 0; i < nr; i++) {
        Integer j = piv(i);
        if (j == i)
          continue;
        this->swapij(i, j);
      }
    } else {
      for (Integer i = nr; --i >= 0;) {
        Integer j = piv(i);
        if (j == i)
          continue;
        this->swapij(i, j);
      }
    }

    return *this;
  }

  Symmatrix& Symmatrix::principal_submatrix(const Indexmatrix& ind, Symmatrix& S) const {
    assert((ind.dim() == 0) || ((min(ind) >= 0) && (max(ind) < nr)));
    S.newsize(ind.dim()); chk_set_init(S, 1);
    // for(Integer i=0;i<ind.dim();i++){
    //   Integer indi=ind(i);
    //   S(i,i)=(*this)(indi,indi);
    //   for(Integer j=i+1;j<ind.dim();j++)
    //     S(i,j)=(*this)(indi,ind(j));
    // }
    Real* Smp = S.m;
    const Integer* const indend = ind.get_store() + ind.dim();
    for (const Integer* indip = ind.get_store(); indip < indend;) {
      Integer indi = (*indip++);
      const Real* const mp = m + indi * nr - (indi * (indi - 1)) / 2 - indi;
      (*Smp++) = *(mp + indi);
      for (const Integer* indjp = indip; indjp < indend;) {
        Integer indj = (*indjp++);
        if (indj >= indi)
          (*Smp++) = *(mp + indj);
        else
          (*Smp++) = *(S.m + (indj * nr - (indj * (indj - 1)) / 2) + indi - indj);
      }
    }
    return S;
  }

  Symmatrix& Symmatrix::delete_principal_submatrix(const Indexmatrix& ind, bool sorted_increasingly) {
    Integer idim = ind.dim();
    if (idim == 0)
      return *this;
    if (idim == nr) {
      init(0, 0.);
      return *this;
    }
    const Integer* indp;
    Indexmatrix sind;
    if (!sorted_increasingly) {
      sortindex(ind, sind);
      sind = ind(sind);
      indp = sind.get_store();
    } else {
      indp = ind.get_store();
    }
    const Integer* const indpend = indp + idim;
    assert((0 <= indp[0]) && (indp[idim - 1] < nr));
    const Real* mp = m;
    Real* np = m;
    for (Integer i = 0; i < nr; i++) {
      if ((indp != indpend) && (i == *indp)) {
        //skip this row
        mp += nr - i;
        indp++;
        assert((indp == indpend) || (*(indp - 1) < *indp));
        continue;
      }
      *np++ = *mp++;
      Integer j = i + 1;
      for (const Integer* jindp = indp; jindp != indpend; jindp++) {
        while (j++ < *jindp) {
          *np++ = *mp++;
        }
        mp++;
      }
      while (j++ < nr)
        *np++ = *mp++;
    }
    nr = nr - idim;
    return *this;
  }

  Symmatrix& Symmatrix::enlarge_below(Integer addn) {
    if (addn <= 0)
      return *this;
    chk_set_init(*this, 0);
    Integer nn = nr + addn;
    Integer n2 = (nn * (nn + 1)) / 2;
    if (n2 < mem_dim) {
      Real* mp = m + (nr * (nr + 1)) / 2;
      Real* np = m + n2 - (addn * (addn + 1)) / 2;
      for (Integer i = 1; i < nr; i++) {
        np -= addn;
        for (Integer j = i; --j >= 0;)
          *(--np) = *(--mp);
      }
      nr += addn;
      return *this;
    }
    // get new storage
    Real* nm = 0;
    mem_dim = Integer(memarray->get(max(n2, 2 * mem_dim), nm));
    if (mem_dim < n2)
      MEmessage(MEmem(n2,
        "Symmatrix::enlarge_below() not enough memory", MTsymmetric));
    Real* mp = m;
    Real* np = nm;
    for (Integer i = 0; i < nr; i++) {
      for (Integer j = i; j < nr; j++)
        *np++ = *mp++;
      np += addn;
    }
    memarray->free(m);
    m = nm;
    nr += addn;
    return *this;
  }

  Symmatrix& Symmatrix::enlarge_below(Integer addn, Real d) {
    if (addn <= 0)
      return *this;
    Integer nn = nr + addn;
    Integer n2 = (nn * (nn + 1)) / 2;
    if (n2 < mem_dim) {
      Real* mp = m + (nr * (nr + 1)) / 2;
      Real* np = m + n2;
      for (Integer j = 0; j < (addn * (addn + 1)) / 2; j++)
        *(--np) = d;
      for (Integer i = 1; i < nr; i++) {
        for (Integer j = 0; j < addn; j++)
          *(--np) = d;
        for (Integer j = i; --j >= 0;)
          *(--np) = *(--mp);
      }
      for (Integer j = 0; j < addn; j++)
        *(--np) = d;
      nr += addn;
      return *this;
    }
    // get new storage
    Real* nm = 0;
    mem_dim = Integer(memarray->get(max(n2, 2 * mem_dim), nm));
    if (mem_dim < n2)
      MEmessage(MEmem(n2,
        "Symmatrix::enlarge_below() not enough memory", MTsymmetric));
    Real* mp = m;
    Real* np = nm;
    for (Integer i = 0; i < nr; i++) {
      for (Integer j = i; j < nr; j++)
        *np++ = *mp++;
      for (Integer j = 0; j < addn; j++)
        *np++ = d;
    }
    for (Integer i = 0; i < addn; i++) {
      for (Integer j = i; j < addn; j++)
        *np++ = d;
    }
    memarray->free(m);
    m = nm;
    nr += addn;
    return *this;
  }


  // *****************************************************************************
  //                  Interaction with other classes
  // *****************************************************************************


  // **************************************************************************
  //                               friends
  // **************************************************************************

  void svec(const Symmatrix& A, Matrix& sv,
    Real a, bool add,
    Integer startindex_vec, Integer startindex_A, Integer blockdim) {
    chk_init(A);
    assert(startindex_A >= 0);
    if (blockdim < 0) {
      blockdim = A.nr - startindex_A;
    }
    const Integer step_to_next_row = A.nr - (startindex_A + blockdim);
    assert(step_to_next_row >= 0);
    Integer vdim = (blockdim * (blockdim + 1)) / 2;
    const Real asqrt2 = a * ::sqrt(2.); //1.41421356237310;
    const Real* mp = A.m + A.nr * startindex_A - (startindex_A * (startindex_A - 1)) / 2;
    if (add == false) {
      if (startindex_vec < 0) {
        sv.newsize(vdim, 1);
        startindex_vec = 0;
      }
      assert(sv.dim() >= vdim + startindex_vec);
      Real* vp = sv.get_store() + startindex_vec;
      for (Integer i = blockdim; --i >= 0;) {
        (*vp++) = a * (*mp++);
        for (Integer j = i; --j >= 0;) {
          (*vp++) = asqrt2 * (*mp++);
        }
        mp += step_to_next_row;
      }
    } else {
      if (startindex_vec < 0)
        startindex_vec = 0;
      assert(sv.dim() >= vdim + startindex_vec);
      Real* vp = sv.get_store() + startindex_vec;
      for (Integer i = blockdim; --i >= 0;) {
        (*vp++) += a * (*mp++);
        for (Integer j = i; --j >= 0;) {
          (*vp++) += asqrt2 * (*mp++);
        }
        mp += step_to_next_row;
      }
    }
    chk_set_init(sv, 1);
  }

  void sveci(const Matrix& sv, Symmatrix& A,
    Real a, bool add,
    Integer startindex_vec, Integer startindex_A, Integer blockdim) {
    chk_init(sv);
    assert(startindex_vec >= 0);
    if (blockdim < 0) {
      if ((add == false) && (startindex_A < 0)) {
        // reinitialize and require exact fit
        blockdim = Integer(::sqrt(Real(8 * (sv.dim() - startindex_vec) + 1)) - 1 + .1) / 2;
#if (CONICBUNDLE_DEBUG>=1)
        if ((blockdim * (blockdim + 1)) / 2 != sv.dim() - startindex_vec)
          MEmessage(MatrixError(ME_unspec,
            "sveci(): dimension of svec does not permit conversion",
            MTsymmetric));
#endif
        A.newsize(blockdim);
        startindex_A = 0.;
      } else {
        if (startindex_A < 0)
          startindex_A = 0;
        blockdim = A.nr - startindex_A;
      }
    } else if (startindex_A < 0) {
      if (add == false) {
        // reinitialize 
        A.newsize(blockdim);
      }
      startindex_A = 0.;
    }

    const Integer step_to_next_row = A.nr - (startindex_A + blockdim);
    assert(sv.dim() >= (blockdim * (blockdim + 1)) / 2 + startindex_vec);
    const Real asqrt2inv = a / ::sqrt(2.); //1./1.41421356237310;
    const Real* vp = sv.get_store() + startindex_vec;
    Real* mp = A.m + A.nr * startindex_A - (startindex_A * (startindex_A - 1)) / 2;
    if (add == false) {
      for (Integer i = blockdim; --i >= 0;) {
        (*mp++) = a * (*vp++);
        for (Integer j = i; --j >= 0;) {
          (*mp++) = asqrt2inv * (*vp++);
        }
        mp += step_to_next_row;
      }
    } else {
      for (Integer i = blockdim; --i >= 0;) {
        (*mp++) += a * (*vp++);
        for (Integer j = i; --j >= 0;) {
          (*mp++) += asqrt2inv * (*vp++);
        }
        mp += step_to_next_row;
      }
    }
    chk_set_init(A, 1);
  }

  // initialize from an svec stored in a real array (or matrix)
  void Symmatrix::init_svec(Integer inr, const Real* dp, Integer incr, Real d) {
    assert(inr >= 0);
    if (nr != inr)
      newsize(inr);
    Real* mp = m;
    const Real invsqrt2 = 1. / ::sqrt(2);
    if (d == 1.) {
      if (incr == 1) {
        for (Integer i = inr; --i >= 0;) {
          (*mp++) = (*dp++);
          for (Integer j = i; --j >= 0;)
            (*mp++) = (*dp++) * invsqrt2;
        }
      } else {
        for (Integer i = inr; --i >= 0;) {
          (*mp++) = (*dp);
          dp += incr;
          for (Integer j = i; --j >= 0;) {
            (*mp++) = (*dp) * invsqrt2;
            dp += incr;
          }
        }
      }
    } else {
      const Real dinv = d * invsqrt2;
      if (incr == 1) {
        for (Integer i = inr; --i >= 0;) {
          (*mp++) = (*dp++) * d;
          for (Integer j = i; --j >= 0;)
            (*mp++) = (*dp++) * dinv;
        }
      } else {
        const Real dinv = d * invsqrt2;
        for (Integer i = inr; --i >= 0;) {
          (*mp++) = (*dp) * d;
          dp += incr;
          for (Integer j = i; --j >= 0;) {
            (*mp++) = (*dp) * dinv;
            dp += incr;
          }
        }
      }
    }
    chk_set_init(*this, 1);
  }

  // initialize from an svec stored in a real array (or matrix)
  void Symmatrix::store_svec(Real* dp, Integer incr, Real d) const {
    const Real* mp = m;
    const Real sqrt2 = ::sqrt(2);
    if (d == 1.) {
      if (incr == 1) {
        for (Integer i = nr; --i >= 0;) {
          (*dp++) = (*mp++);
          for (Integer j = i; --j >= 0;)
            (*dp++) = (*mp++) * sqrt2;
        }
      } else {
        for (Integer i = nr; --i >= 0;) {
          (*dp) = (*mp++);
          dp += incr;
          for (Integer j = i; --j >= 0;) {
            (*dp) = (*mp++) * sqrt2;
            dp += incr;
          }
        }
      }
    } else {
      const Real dsqrt2 = d * sqrt2;
      if (incr == 1) {
        for (Integer i = nr; --i >= 0;) {
          (*dp++) = (*mp++);
          for (Integer j = i; --j >= 0;)
            (*dp++) = (*mp++) * dsqrt2;
        }
      } else {
        for (Integer i = nr; --i >= 0;) {
          (*dp) = (*mp++);
          dp += incr;
          for (Integer j = i; --j >= 0;) {
            (*dp) = (*mp++) * dsqrt2;
            dp += incr;
          }
        }
      }
    }
  }


  void skron(const Symmatrix& A, const Symmatrix& B, Symmatrix& S, Real alpha, bool add, Integer startindex_S) {
    chk_add(A, B);
    const Real sqrt2 = std::sqrt(2.); // 1.41421356237310;
    Integer n = A.rowdim();
    Integer sdim = (n * (n + 1)) / 2;
    if (startindex_S < 0) {
      if (add == false) {
        S.newsize(sdim); chk_set_init(S, 1);
      }
      startindex_S = 0;
    }
    assert(S.nr >= startindex_S + sdim);
    const Integer step_to_next_row = S.rowdim() - (startindex_S + sdim);
    assert(step_to_next_row >= 0);
    Matrix sv(sdim, 1);
    chk_set_init(sv, 1);
    Real* Sp = S.get_store() + startindex_S * S.nr - (startindex_S * (startindex_S - 1)) / 2;
    Matrix Ah(A);
    Matrix Bh(B);
    for (Integer i = 0; i < n; i++) {

      //row corresponding to diagonal
      {
        const Real* ai = Ah.get_store() + i * n + i; //points now to i-th elem of i-th column of A
        const Real* bi = Bh.get_store() + i * n + i; //points now to i-th elem of i-thcolumn of B
        Real* svp = sv.get_store();
        for (Integer l = i; l < n; l++) {
          const Real ail = *ai++; //ai points afterwards to the element A(l+1,i)
          const Real bil = *bi++; //bi points afterwards to the element B(l+1,i)
          const Real* aik = ai;
          const Real* bik = bi;
          (*svp++) = ail * bil;
          for (Integer k = l + 1; k < n; k++) {
            (*svp++) = sqrt2 * (ail * (*bik++) + (*aik++) * bil) / 2.;
          }
        }
      }
      if (add)
        mat_xpeya(sdim, Sp, sv.get_store(), alpha);
      else
        mat_xeya(sdim, Sp, sv.get_store(), alpha);
      Sp += sdim + step_to_next_row;
      sdim--;

      for (Integer j = i + 1; j < n; j++) {
        const Real* ai = Ah.get_store() + i * n + i; //points now to i-th elem in i-th column of A
        const Real* bi = Bh.get_store() + i * n + i; //points now to i-th elem in i-th column of B
        const Real* aj = Ah.get_store() + j * n + i; //points now to i-th elem in j-th column of A
        const Real* bj = Bh.get_store() + j * n + i; //points now to i-th elem in j-th column of B

        //row correpsonding to offdiagonal
        Real* svp = sv.get_store();
        {
          const Real* aik = ai + j - i; //points now to A(j,i)
          const Real* bik = bi + j - i; //points now to B(j,i)
          const Real* ajk = aj + j - i; //points now to A(j,j)
          const Real* bjk = bj + j - i; //points now to B(j,j)
          const Real aii = *ai++; //afterwards ai points to A(i+1,i)
          const Real bii = *bi++; //afterwards bi points to B(i+1,i)
          const Real aji = *aj++; //afterwards aj points now to A(i+1,j)
          const Real bji = *bj++; //afterwards bj points now to B(i+1,j)
          for (Integer k = j; k < n; k++) {
            (*svp++) = (aii * (*bjk++) + (*aik++) * bji + aji * (*bik++) + (*ajk++) * bii) / 2.;
          }
        }
        for (Integer l = i + 1; l < n; l++) {
          const Real ail = *ai++; //afterwards ai points to A(l+1,i)
          const Real bil = *bi++; //afterwards bi points to B(l+1,i)
          const Real ajl = *aj++; //afterwards aj points to A(l+1,j)
          const Real bjl = *bj++; //afterwards bj points to B(l+1,j)
          (*svp++) = sqrt2 * (ail * bjl + ajl * bil) / 2.;
          const Real* aik = ai; //aik points now to A(l+1,i)
          const Real* bik = bi; //bik points now to B(l+1,i)
          const Real* ajk = aj; //ajk points now to A(l+1,j)
          const Real* bjk = bj; //bjk points now to B(l+1,j)
          for (Integer k = l + 1; k < n; k++) {
            (*svp++) = (ail * (*bjk++) + (*aik++) * bjl + ajl * (*bik++) + (*ajk++) * bil) / 2.;
          }
        }
        if (add)
          mat_xpeya(sdim, Sp, sv.get_store(), alpha);
        else
          mat_xeya(sdim, Sp, sv.get_store(), alpha);
        Sp += sdim + step_to_next_row;
        sdim--;
      }
    }

  }

  void symscale(const Symmatrix& A, const Matrix& B, Symmatrix& S, Real alpha, Real beta, int Btrans) {
    assert(A.get_store() != S.get_store());
    if (Btrans == 0) {
      chk_mult(A, B);
      Integer n = B.coldim();
      Integer m = B.rowdim();
      Matrix tmp;
      if (alpha != 0.)
        genmult(A, B, tmp);
      const Real* bmp = B.get_store();
      if (beta == 0.) {
        if (alpha == 0.)
          S.init(n, 0.);
        else {
          S.newsize(n); chk_set_init(S, 1);
          Real* smp = S.get_store();
          for (Integer i = 0; i < n; i++) {
            const Real* amp = tmp.get_store() + i * m;
            for (Integer j = i; j < n; j++) {
              (*smp++) = alpha * mat_ip(m, bmp, amp);
              amp += m;
            }
            bmp += m;
          }
        }
      } else if (beta == 1.) {
        if (S.nr != n) {
          MEmessage(MatrixError(ME_dim, "symscale: dimensions do not match", MTmatrix));;
        }
        if (alpha != 0) {
          Real* smp = S.get_store();
          for (Integer i = 0; i < n; i++) {
            const Real* amp = tmp.get_store() + i * m;
            for (Integer j = i; j < n; j++) {
              (*smp++) += alpha * mat_ip(m, bmp, amp);
              amp += m;
            }
            bmp += m;
          }
        }
      } else {
        if (alpha == 0) {
          if (S.nr != n) {
            MEmessage(MatrixError(ME_dim, "symscale: dimensions do not match", MTmatrix));;
          }
          S *= beta;
        } else {
          Real* smp = S.get_store();
          for (Integer i = 0; i < n; i++) {
            const Real* amp = tmp.get_store() + i * m;
            for (Integer j = i; j < n; j++, smp++) {
              (*smp) = beta * (*smp) + alpha * mat_ip(m, bmp, amp);
              amp += m;
            }
            bmp += m;
          }
        }
      }
      //assert(norm2(S-transpose(B)*A*B)<1e-10);
    } else {
      chk_mult(B, A);
      Integer m = B.coldim();
      Integer n = B.rowdim();
      Matrix tmp;
      if (alpha != 0.)
        genmult(A, B, tmp, 1., 0., 1);
      const Real* amp = tmp.get_store();
      if (beta == 0.) {
        if (alpha == 0.)
          S.init(n, 0.);
        else {
          S.newsize(n); chk_set_init(S, 1);
          Real* smp = S.get_store();
          for (Integer i = n; i > 0; --i) {
            const Real* bmp = B.get_store() + n - i;
            mat_xeya(i, smp, bmp, alpha * (*amp++));
            for (Integer j = 1; j < m; j++) {
              bmp += n;
              mat_xpeya(i, smp, bmp, alpha * (*amp++));
            }
            smp += i;
          }
        }
      } else if (beta == 1.) {
        if (S.nr != n) {
          MEmessage(MatrixError(ME_dim, "symscale: dimensions do not match", MTmatrix));;
        }
        if (alpha != 0) {
          Real* smp = S.get_store();
          for (Integer i = n; i > 0; --i) {
            const Real* bmp = B.get_store() + n - i;
            mat_xpeya(i, smp, bmp, alpha * (*amp++));
            for (Integer j = 1; j < m; j++) {
              bmp += n;
              mat_xpeya(i, smp, bmp, alpha * (*amp++));
            }
            smp += i;
          }
        }
      } else {
        if (alpha == 0) {
          if (S.nr != n) {
            MEmessage(MatrixError(ME_dim, "symscale: dimensions do not match", MTmatrix));;
          }
          S *= beta;
        } else {
          Real* smp = S.get_store();
          for (Integer i = n; i > 0; --i) {
            const Real* bmp = B.get_store() + n - i;
            mat_xbpeya(i, smp, bmp, alpha * (*amp++), beta);
            for (Integer j = 1; j < m; j++) {
              bmp += n;
              mat_xpeya(i, smp, bmp, alpha * (*amp++));
            }
            smp += i;
          }
        }
      }
      //assert(norm2(S-B*A*transpose(B))<1e-10);
    }
  }


  Matrix diag(const Symmatrix& A) {
    chk_init(A);
    Matrix v(A.nr, 1);
    for (Integer i = 0; i < A.nr; i++)
      v(i) = A(i, i);
    chk_set_init(v, 1);
    return v;
  }

  Symmatrix Diag(const Matrix& A) {
    chk_init(A);
    Symmatrix S(A.dim(), 0.);
    for (Integer i = 0; i < A.dim(); i++)
      S(i, i) = A(i);
    chk_set_init(S, 1);
    return S;
  }

  Matrix sumrows(const Symmatrix& A) {
#if (CONICBUNDLE_DEBUG>=1)
    if (A.nr == 0)
      MEmessage(MEdim(A.nr, A.nr, 0, 0, "sumrows(const Symmatrix&) dimension zero", MTsymmetric));
    if (!A.is_init)
      MEmessage(MatrixError(ME_unspec, "sumrows(const Symmatrix& A) A is not initialized", MTsymmetric));
#endif
    Matrix v(1, A.nr);
    Real sum;
    Integer i, j;
    for (j = 0; j < A.nr; j++) {
      sum = 0.;
      for (i = 0; i < A.nr; i++)
        sum += A(i, j);
      v(j) = sum;
    }
#if (CONICBUNDLE_DEBUG>=1)
    v.set_init(1);
#endif 
    return v;
  }

  Matrix sumcols(const Symmatrix& A) {
#if (CONICBUNDLE_DEBUG>=1)
    if (A.nr == 0)
      MEmessage(MEdim(A.nr, A.nr, 0, 0, "sumcols(const Symmatrix&) dimension zero", MTsymmetric));
    if (!A.is_init)
      MEmessage(MatrixError(ME_unspec, "sumcols(const Symmatrix& A) A is not initialized", MTsymmetric));
#endif
    Matrix v(A.nr, 1);
    Real sum;
    Integer i, j;
    for (i = 0; i < A.nr; i++) {
      sum = 0.;
      for (j = 0; j < A.nr; j++)
        sum += A(i, j);
      v(i) = sum;
    }
#if (CONICBUNDLE_DEBUG>=1)
    v.set_init(1);
#endif 
    return v;
  }

  Real sum(const Symmatrix& A) {
#if (CONICBUNDLE_DEBUG>=1)
    if (A.nr == 0)
      MEmessage(MEdim(A.nr, A.nr, 0, 0, "sum(const Symmatrix&) dimension zero", MTsymmetric));
    if (!A.is_init)
      MEmessage(MatrixError(ME_unspec, "sum(const Symmatrix& A) A is not initialized", MTsymmetric));
#endif
    Real s = 0.;
    Real s2 = 0;
    Integer i, j;
    const Real* ap = A.m;
    for (i = A.nr; --i >= 0;) {
      s += (*ap++);
      s2 = 0.;
      for (j = i; --j >= 0;) {
        s2 += (*ap++);
      }
      s += 2. * s2;
    }
    return s;
  }

  Matrix maxrows(const Symmatrix& A) {
#if (CONICBUNDLE_DEBUG>=1)
    if (A.nr == 0)
      MEmessage(MEdim(A.nr, A.nr, 0, 0, "maxrows(const Symmatrix&) dimension zero", MTsymmetric));
    if (!A.is_init)
      MEmessage(MatrixError(ME_unspec, "maxrows(const Symmatrix& A) A is not initialized", MTsymmetric));
#endif
    Matrix v(1, A.nr);
    Real maxd;
    Integer i, j;
    for (j = 0; j < A.nr; j++) {
      maxd = A(0, j);
      for (i = 1; i < A.nr; i++)
        maxd = max(maxd, A(i, j));
      v(j) = maxd;
    }
#if (CONICBUNDLE_DEBUG>=1)
    v.set_init(1);
#endif 
    return v;
  }

  Matrix maxcols(const Symmatrix& A) {
#if (CONICBUNDLE_DEBUG>=1)
    if (A.nr == 0)
      MEmessage(MEdim(A.nr, A.nr, 0, 0, "maxcols(const Symmatrix&) dimension zero", MTsymmetric));
    if (!A.is_init)
      MEmessage(MatrixError(ME_unspec, "maxcols(const Symmatrix& A) A is not initialized", MTsymmetric));
#endif
    Matrix v(A.nr, 1);
    Real maxd;
    Integer i, j;
    for (i = 0; i < A.nr; i++) {
      maxd = A(i, 0);
      for (j = 1; j < A.nr; j++)
        maxd = max(maxd, A(i, j));
      v(i) = maxd;
    }
#if (CONICBUNDLE_DEBUG>=1)
    v.set_init(1);
#endif 
    return v;
  }

  Real max(const Symmatrix& A) {
#if (CONICBUNDLE_DEBUG>=1)
    if (A.nr == 0)
      MEmessage(MEdim(A.nr, A.nr, 0, 0, "max(const Symmatrix&) dimension zero", MTsymmetric));
    if (!A.is_init)
      MEmessage(MatrixError(ME_unspec, "max(const Symmatrix& A) A is not initialized", MTsymmetric));
#endif
    Real maxd = min_Real;
    const Real* const mpend = A.m + (A.nr * (A.nr + 1)) / 2;
    for (const Real* mp = A.m; mp < mpend; mp++) {
      if (maxd < *mp)
        maxd = *mp;
    }
    return maxd;
  }

  Matrix minrows(const Symmatrix& A) {
#if (CONICBUNDLE_DEBUG>=1)
    if (A.nr == 0)
      MEmessage(MEdim(A.nr, A.nr, 0, 0, "minrows(const Symmatrix&) dimension zero", MTsymmetric));
    if (!A.is_init)
      MEmessage(MatrixError(ME_unspec, "minrows(const Symmatrix& A) A is not initialized", MTsymmetric));
#endif
    Matrix v(1, A.nr);
    Real mind;
    Integer i, j;
    for (j = 0; j < A.nr; j++) {
      mind = A(0, j);
      for (i = 1; i < A.nr; i++)
        mind = min(mind, A(i, j));
      v(j) = mind;
    }
#if (CONICBUNDLE_DEBUG>=1)
    v.set_init(1);
#endif 
    return v;
  }

  Matrix mincols(const Symmatrix& A) {
#if (CONICBUNDLE_DEBUG>=1)
    if (A.nr == 0)
      MEmessage(MEdim(A.nr, A.nr, 0, 0, "mincols(const Symmatrix&) dimension zero", MTsymmetric));
    if (!A.is_init)
      MEmessage(MatrixError(ME_unspec, "mincols(const Symmatrix& A) A is not initialized", MTsymmetric));
#endif
    Matrix v(A.nr, 1);
    Real mind;
    Integer i, j;
    for (i = 0; i < A.nr; i++) {
      mind = A(i, 0);
      for (j = 1; j < A.nr; j++)
        mind = min(mind, A(i, j));
      v(i) = mind;
    }
#if (CONICBUNDLE_DEBUG>=1)
    v.set_init(1);
#endif 
    return v;
  }

  Real min(const Symmatrix& A) {
#if (CONICBUNDLE_DEBUG>=1)
    if (A.nr == 0)
      MEmessage(MEdim(A.nr, A.nr, 0, 0, "min(const Symmatrix&) dimension zero", MTsymmetric));
    if (!A.is_init)
      MEmessage(MatrixError(ME_unspec, "min(const Symmatrix& A) A is not initialized", MTsymmetric));
#endif
    Real mind = max_Real;
    const Real* const mpend = A.m + (A.nr * (A.nr + 1)) / 2;
    for (const Real* mp = A.m; mp < mpend; mp++) {
      if (mind > *mp)
        mind = *mp;
    }
    return mind;
  }

  Real trace(const Symmatrix& A) {
#if (CONICBUNDLE_DEBUG>=1)
    if (!A.is_init)
      MEmessage(MatrixError(ME_unspec, "trace(const Symmatrix& A) A is not initialized", MTsymmetric));
#endif
    Real sum = 0.;
    for (Integer i = 0; i < A.nr; i++) sum += A(i, i);
    return sum;
  }

  Real ip(const Symmatrix& A, const Symmatrix& B) {
#if (CONICBUNDLE_DEBUG>=1)
    if (B.nr != A.nr)
      MEmessage(MEdim(A.nr, A.nr, B.nr, B.nr, "ip(const Symmatrix&,const Symmatrix&) wrong dimensions", MTsymmetric));
    if ((!A.is_init) || (!B.is_init))
      MEmessage(MatrixError(ME_unspec, "ip(const Symmatrix&,const Symmatrix&) not initialized", MTsymmetric));
#endif
    Real sum = 0.;
    Real s2;
    Integer i, j;
    const Real* ap = A.m;
    const Real* bp = B.m;
    for (i = A.nr; --i >= 0;) {
      sum += (*ap++) * (*bp++);
      s2 = 0.;
      for (j = i; --j >= 0;) {
        s2 += (*ap++) * (*bp++);
      }
      sum += 2. * s2;
    }
    return sum;
  }

  Real ip(const Symmatrix& A, const Matrix& B) {
    chk_add(A, B);
    Matrix C(A);
    return ip(C, B);
  }

  Real ip(const Matrix& A, const Symmatrix& B) {
    chk_add(A, B);
    Matrix C(B);
    return ip(A, C);
  }

  Symmatrix abs(const Symmatrix& A) {
    chk_init(A);
    Symmatrix B; B.newsize(A.nr);
    Real* ap = A.m;
    Real* bp = B.m;
    for (Integer i = (A.nr * (A.nr + 1)) / 2; --i >= 0;)
      (*bp++) = fabs(*ap++);
    chk_set_init(B, 1);
    return B;
  }

  std::ostream& operator<<(std::ostream& o, const Symmatrix& A) {
    chk_init(A);
    o << A.nr << '\n';
    Integer i, j;
    for (i = 0; i < A.nr; i++) {
      for (j = i; j < A.nr; j++) o << ' ' << A(i, j);
      o << '\n';
    }
    return o;
  }

  std::istream& operator>>(std::istream& in, Symmatrix& A) {
    Real d;
    in >> d;
    Integer nr = Integer(d + .5);
    if (nr < 0)
      MEmessage(MEdim(nr, nr, 0, 0, "operator>>(std::istream&,Symmatrix&) dimension negative", MTsymmetric));
    A.newsize(nr);
    for (Integer i = 0; i < nr; i++)
      for (Integer j = i; j < nr; j++)
        in >> A(i, j);
    chk_set_init(A, 1);
    return in;
  }

  // **************************************************************************
  //                       Matrix:: Symmatrix specific implementations
  // **************************************************************************

  Matrix& Matrix::xeya(const Symmatrix& A, Real d) {
    chk_init(A);
    newsize(A.nr, A.nr);
    const Real* matp = A.m;
    if (d == 1.) {
      for (Integer i = 0; i < nr; i++) {
        Real* mcp = m + i * nr + i;
        Real* mrp = mcp + nr;
        (*mcp++) = (*matp++);
        for (Integer j = i + 1; j < nr; j++, mrp += nr) {
          (*mcp++) = (*mrp) = (*matp++);
        }
      }
    } else if (d == -1.) {
      for (Integer i = 0; i < nr; i++) {
        Real* mcp = m + i * nr + i;
        Real* mrp = mcp + nr;
        (*mcp++) = -(*matp++);
        for (Integer j = i + 1; j < nr; j++, mrp += nr) {
          (*mcp++) = (*mrp) = -(*matp++);
        }
      }
    } else {
      for (Integer i = 0; i < nr; i++) {
        Real* mcp = m + i * nr + i;
        Real* mrp = mcp + nr;
        (*mcp++) = (*matp++) * d;
        for (Integer j = i + 1; j < nr; j++, mrp += nr) {
          (*mcp++) = (*mrp) = (*matp++) * d;
        }
      }
    }
    chk_set_init(*this, 1);
    return *this;
  }

  Matrix& Matrix::xpeya(const Symmatrix& A, Real d) {
    chk_add(*this, A);
    const Real* matp = A.m;
    if (d == 1.) {
      for (Integer i = 0; i < nr; i++) {
        Real* mcp = m + i * nr + i;
        Real* mrp = mcp + nr;
        (*mcp++) += (*matp++);
        for (Integer j = i + 1; j < nr; j++, mrp += nr) {
          (*mcp++) += (*matp);
          (*mrp) += (*matp++);
        }
      }
    } else if (d == -1.) {
      for (Integer i = 0; i < nr; i++) {
        Real* mcp = m + i * nr + i;
        Real* mrp = mcp + nr;
        (*mcp++) -= (*matp++);
        for (Integer j = i + 1; j < nr; j++, mrp += nr) {
          (*mcp++) -= (*matp);
          (*mrp) -= (*matp++);
        }
      }
    } else {
      for (Integer i = 0; i < nr; i++) {
        Real* mcp = m + i * nr + i;
        Real* mrp = mcp + nr;
        (*mcp++) += (*matp++) * d;
        for (Integer j = i + 1; j < nr; j++, mrp += nr) {
          Real f = d * (*matp++);
          (*mcp++) += f;
          (*mrp) += f;
        }
      }
    }
    chk_set_init(*this, 1);
    return *this;
  }

}

