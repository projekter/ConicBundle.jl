/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  Matrix/chol.cxx
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
#include "mymath.hxx"
#include "symmat.hxx"



namespace CH_Matrix_Classes {

  int Symmatrix::Chol_factor(Real tol) {
    for (Integer k = 0; k < nr; k++) {

      //---- compute factorization
      /*
        (*this)(k,k)=sqrt((*this)(k,k));
        for(Integer i=k+1;i<nr;i++) (*this)(k,i)/=(*this)(k,k);
        for(Integer i=k+1;i<nr;i++) {
          for(j=i;j<nr;j++){
            (*this)(i,j)-=(*this)(k,j)*(*this)(k,i);
          }
        }
      */
      Real* mp2 = m + nr * k - ((k - 1) * k) / 2;
      if (*mp2 < tol)
        return k + 1;
      Real d = (*mp2 = ::sqrt(*mp2));
      mp2++;
      for (Integer i = nr - k; --i > 0;) {
        (*mp2++) /= d;
      }
      // #ifdef WITH_OMP
      //    const Integer ompthresh=10000;
      //    if (nr-k<2*ompthresh){
      // #endif
      const Real* mp = mp2;
      {
        for (Integer i = nr - k; --i > 0;) {
          Integer j;
          mp -= (j = i);
          const Real d = *mp;
          for (; --j >= 0;) {
            (*mp2++) -= (*mp++) * d;
          }
        }
      }
      // #ifdef WITH_OMP
      //    }
      //    else {
      //      mp2-=nr;
      //      #pragma omp parallel for schedule(static,1)
      //      for(Integer i=k+1;i<nr-ompthresh;i++) {
      //        const Real *rip=mp2+i;       
      //        mat_xpeya(nr-i,m+nr*i-((i-1)*i)/2,rip,-*rip);
      //      } 
      //      const Real *mp=mp2+nr;
      //      Real *cmp=m+nr*(nr-ompthresh)-((nr-ompthresh-1)*(nr-ompthresh))/2;
      //      {for(Integer i=ompthresh+1;--i>0;){
      // 	 Integer j;
      // 	 mp-=(j=i);
      // 	 const Real d=*mp;
      // 	 for(;--j>=0;){
      // 	   (*cmp++)-=(*mp++)*d;
      // 	 }
      //        }}
      //    }
      // #endif
    }
    return 0;
  }

  int Symmatrix::Chol_solve(Matrix& x) const {
    chk_mult(*this, x);

    for (Integer k = 0; k < x.coldim(); k++) { //solve for and overwrite column k of x
      Real* xbase = x.m + k * nr;
      //---- solve Lr=xbase
      for (Integer i = 0; i < nr; i++) {
        Real f = xbase[i];
        Real* rp = xbase;
        const Real* mp = m + i;
        const Integer h = nr - i;
        for (Integer j = nr; --j >= h;) {
          f -= (*rp++) * (*mp);
          mp += j;
        }
        *rp = f / (*mp);
      }

      //---- solve L'r=r
      {
        for (Integer i = nr; --i >= 0;) {
          Real* rp = xbase + i;
          Real f = (*rp++);
          const Real* mp = m + i * nr - (i * (i - 1)) / 2 + 1;
          for (Integer j = nr - i; --j > 0;)
            f -= (*rp++) * (*mp++);
          xbase[i] = f / (*(mp - (nr - i)));
        }
      }
    }

    return 0;
  }

  int Symmatrix::Chol_Lsolve(Matrix& x) const {
    chk_mult(*this, x);


    for (Integer k = 0; k < x.coldim(); k++) { //solve for and overwrite column k of x
      Real* xbase = x.m + k * nr;
      //---- solve Lr=xbase
      for (Integer i = 0; i < nr; i++) {
        Real f = xbase[i];
        Real* rp = xbase;
        const Real* mp = m + i;
        const Integer h = nr - i;
        for (Integer j = nr; --j >= h;) {
          f -= (*rp++) * (*mp);
          mp += j;
        }
        *rp = f / (*mp);
      }
    }

    return 0;
  }

  int Symmatrix::Chol_Ltsolve(Matrix& x) const {
    chk_mult(*this, x);


    for (Integer k = 0; k < x.coldim(); k++) { //solve for and overwrite column k of x
      Real* xbase = x.m + k * nr;
      //---- solve L'r=r
      {
        for (Integer i = nr; --i >= 0;) {
          Real* rp = xbase + i;
          Real f = (*rp++);
          const Real* mp = m + i * nr - (i * (i - 1)) / 2 + 1;
          for (Integer j = nr - i; --j > 0;)
            f -= (*rp++) * (*mp++);
          xbase[i] = f / (*(mp - (nr - i)));
        }
      }
    }

    return 0;
  }

  int Symmatrix::Chol_scaleLi(Symmatrix& S) const {
    chk_mult(*this, S);

    //---- compute tmp=L^{-1}S
    Matrix tmp(S);
    int retval;
    if ((retval = Chol_Lsolve(tmp))) {
      return retval;
    }
    //---- compute tmp=SL^{-T}
    tmp.transpose();

    //---- compute tmp=L^{-1}SL^{-T}
    if ((retval = Chol_Lsolve(tmp))) {
      return retval;
    }
    S.init(tmp);

    return 0;
  }

  int Symmatrix::Chol_scaleLt(Symmatrix& S) const {
    chk_mult(*this, S);

    Matrix tmpmat(S);
    Integer n = tmpmat.rowdim();
    Matrix tmpvec(n, Integer(1));
    Real* sp = S.get_store();
    for (Integer i = 0; i < n; i++) {
      //form the first product for row i
      //tmpvec=(*this)(i,Range(i,n-1))*tmpmat(Range(i,n-1),Range(0,n-1));
      const Real* tp = m + i * n - (i * (i - 1)) / 2;
      const Real* mp = tmpmat.get_store() + i * n + i;
      Real* vp = tmpvec.get_store();
      for (Integer j = i; j < n; j++) {
        *vp++ = mat_ip(n - i, tp, mp);
        mp += n;
      }
      //form the second product for row i
      vp = tmpvec.get_store();
      for (Integer j = i; j < n; j++) {
        //S(i,j)=ip(tmpvec(Range(j:n-1)),(*this)(Range(j,n-1),j);
        *sp++ = mat_ip(n - j, vp, tp);
        vp++;
        tp += n - j;
      }
    }
    return 0;
  }


  int Symmatrix::Chol_Lmult(Matrix& x) const {
    chk_mult(*this, x);

    Matrix tmpvec(nr, 1);
    Integer xnc = x.coldim();

    Real* xp = x.get_store();
    for (Integer j = 0; j < xnc; j++) {
      mat_xey(nr, tmpvec.get_store(), xp);
      const Real* Lp = m;
      const Real* tp = tmpvec.get_store();
      mat_xeya(nr, xp, Lp, *tp);
      xp++;
      Lp += nr;
      tp++;
      for (Integer i = nr; --i > 0;) {
        mat_xpeya(i, xp, Lp, *tp);
        xp++;
        Lp += i;
        tp++;
      }
    }
    return 0;
  }


  int Symmatrix::Chol_Ltmult(Matrix& x) const {
    chk_mult(*this, x);

    Integer xnc = x.coldim();
    Real* xp = x.get_store();
    const Real* Lp = m;

    for (Integer i = nr; i > 0; --i) {
      Real* rowxp = xp++;
      for (Integer j = 0; j < xnc; j++, rowxp += nr) {
        *rowxp = mat_ip(i, Lp, rowxp);
      }
      Lp += i;
    }
    return 0;
  }


  int Symmatrix::Chol_inverse(Symmatrix& S) const {
    chk_init(*this);
    S.newsize(nr);
    Real* sp = S.m;
    Real f;
    Real* mp;
    Real* sp1;
    Integer e, i, j, indi;
    for (e = 0; e < nr; e++) {

      indi = e * nr - ((e - 1) * e) / 2 - e;
      sp = S.m + indi;

      //---- solve Lr=v
      indi += e;
      sp[e] = 1. / (*(m + indi));
      indi++;
      mat_xeya(nr - e - 1, sp + e + 1, m + indi, -sp[e]);
      for (i = e + 1; i < nr; i++) {
        indi += nr - i;
        sp[i] /= (*(m + indi));
        indi++;
        mat_xpeya(nr - i - 1, sp + i + 1, m + indi, -sp[i]);
      }

      //---- solve L'r=r
      indi = (nr * (nr + 1)) / 2 - 1;
      for (i = nr; --i >= e;) {
        mp = m + indi;
        sp1 = sp + i;
        *sp1 /= *mp;
        mp -= nr - i;
        f = *sp1;
        for (j = i; --j >= e;) {
          (*(--sp1)) -= f * (*mp);
          mp -= nr - j;
        }
        indi -= nr - i + 1;
      }
    }
    chk_set_init(S, 1);

    return 0;
  }

  int Symmatrix::Chol_factor(Indexmatrix& piv, Real tol) {
    piv.init(nr, 1, Integer(0));
    Integer rank = 0;
    Integer i, j, k;
    Real d;
    Real* mp;
    Real* mp2;
    for (k = 0; k < nr; k++) {

      //---- find maximal diagonal element or stop if a diagonal element is negative
      Integer maxind = k;
      mp = m + nr * k - ((k - 1) * k) / 2;
      Real maxval = (*mp);
      if (maxval < 0.) {
        piv.nr = rank;
        return 1;
      }
      for (i = k + 1; i < nr; i++) {
        mp += nr - i + 1;
        d = (*mp);
        if (d < 0.) {
          piv.nr = rank;
          return 1;
        }
        if (d > maxval) {
          maxind = i;
          maxval = d;
        }
      }

      //---- check whether maxval is still large enough
      if (maxval < tol) break;
      piv(rank) = maxind;
      rank++;

      //---- exchange row k and row maxind
      if (maxind != k) {
        mp = m + k;
        mp2 = m + maxind;
        for (i = 0; i < k; i++) {
          d = *mp; *mp = *mp2; *mp2 = d;
          mp += nr - i - 1; mp2 += nr - i - 1;
        }
        d = *mp; *mp = *(mp2 = (m + nr * maxind - ((maxind - 1) * maxind) / 2)); *mp2 = d;
        mp2 = mp + (maxind - k) + nr - k - 1; mp++;
        for (i = k + 1; i < maxind; i++) {
          d = *mp; *mp = *mp2; *mp2 = d;
          mp++; mp2 += nr - i - 1;
        }
        for (i = maxind + 1; i < nr; i++) {
          mp++; mp2++;
          d = *mp; *mp = *mp2; *mp2 = d;
        }
      }

      //---- compute factorization
      /*
      (*this)(k,k)=sqrt((*this)(k,k));
      for(i=k+1;i<nr;i++) (*this)(k,i)/=(*this)(k,k);
      for(i=k+1;i<nr;i++) {
          for(j=i;j<nr;j++){
              (*this)(i,j)-=(*this)(k,j)*(*this)(k,i);
          }
      }
      */
      mp = m + nr * k - ((k - 1) * k) / 2;
      d = (*mp = ::sqrt(*mp));
      mp++;
      for (i = nr - k - 1; --i >= 0;) {
        (*mp++) /= d;
      }
      mp2 = mp;
      for (i = nr - k - 1; --i >= 0;) {
        mp -= (j = i + 1);
        d = *mp;
        for (; --j >= 0;) {
          (*mp2++) -= (*mp++) * d;
        }
      }
    }
    piv.nr = rank;
    return 0;
  }

  int Symmatrix::Chol_solve(Matrix& x, const Indexmatrix& piv) const {
    chk_mult(*this, x);

    Real f;
    Integer i, j, h;
    Real* mp;
    Real* rp;

    //permute x
    Integer pivdim = piv.dim();
    for (i = 0; i < pivdim; i++) {
      h = piv(i);
      if (h == i) continue;
      mp = x.m + i;
      rp = x.m + h;
      mat_swap(x.coldim(), mp, nr, rp, nr);
    }

    Integer k;
    Real* xbase;

    for (k = 0; k < x.coldim(); k++) { //solve for and overwrite column k of x
      xbase = x.m + k * nr;
      //---- solve Lr=xbase
      for (i = 0; i < pivdim; i++) {
        f = xbase[i];
        rp = xbase;
        mp = m + i;
        h = nr - i;
        for (j = nr; --j >= h;) {
          f -= (*rp++) * (*mp);
          mp += j;
        }
        *rp = f / (*mp);
      }

      //---- solve L'r=r
      for (i = pivdim; --i >= 0;) {
        rp = xbase + i;
        f = (*rp++);
        mp = m + i * nr - (i * (i - 1)) / 2 + 1;
        for (j = nr; --j > i;) f -= (*rp++) * (*mp++);
        xbase[i] = f / (*(mp - (nr - i)));
      }
    }

    //repermute x
    for (i = pivdim; --i >= 0;) {
      h = piv(i);
      if (h == i) continue;
      mp = x.m + i;
      rp = x.m + h;
      mat_swap(x.coldim(), mp, nr, rp, nr);
    }

    return 0;
  }

  int Symmatrix::Chol_inverse(Symmatrix& S, const Indexmatrix& piv) const {
    chk_init(*this);
    S.newsize(nr);
    Real* sp = S.m;
    Real f;
    Real* mp;
    Real* sp1;
    Integer e, i, j, indi;
    for (e = 0; e < nr; e++) {

      indi = e * nr - ((e - 1) * e) / 2 - e;
      sp = S.m + indi;

      //---- solve Lr=v
      indi += e;
      sp[e] = 1. / (*(m + indi));
      indi++;
      mat_xeya(nr - e - 1, sp + e + 1, m + indi, -sp[e]);
      for (i = e + 1; i < nr; i++) {
        indi += nr - i;
        sp[i] /= (*(m + indi));
        indi++;
        mat_xpeya(nr - i - 1, sp + i + 1, m + indi, -sp[i]);
      }

      //---- solve L'r=r
      indi = (nr * (nr + 1)) / 2 - 1;
      for (i = nr; --i >= e;) {
        mp = m + indi;
        sp1 = sp + i;
        *sp1 /= *mp;
        mp -= nr - i;
        f = *sp1;
        for (j = i; --j >= e;) {
          (*(--sp1)) -= f * (*mp);
          mp -= nr - j;
        }
        indi -= nr - i + 1;
      }
    }
    chk_set_init(S, 1);

    //repermute S
    for (i = nr; --i >= 0;) {
      e = piv(i);
      if (e == i) continue;
      mp = S.m + i;
      sp1 = S.m + e;
      for (j = 0; j < i; j++) {
        f = *mp; *mp = *sp1; *sp1 = f;
        mp += nr - j - 1; sp1 += nr - j - 1;
      }
      f = *mp; *mp = *(sp1 = (S.m + nr * e - ((e - 1) * e) / 2)); *sp1 = f;
      sp1 = mp + (e - i) + nr - i - 1; mp++;
      for (j = i + 1; j < e; j++) {
        f = *mp; *mp = *sp1; *sp1 = f;
        mp++; sp1 += nr - j - 1;
      }
      for (j = e + 1; j < nr; j++) {
        mp++; sp1++;
        f = *mp; *mp = *sp1; *sp1 = f;
      }
    }

    return 0;
  }

}

