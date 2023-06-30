/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  Matrix/qr.cxx
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
#include "heapsort.hxx"
#include "gb_rand.hxx"


using namespace CH_Tools;

namespace CH_Matrix_Classes {


  int Matrix::QR_factor(Real tol) {
    chk_init(*this);
#if (CONICBUNDLE_DEBUG>=1)
    if (nr < nc)
      MEmessage(MatrixError(ME_unspec, "QR_factor(): #rows < #cols ", MTmatrix));
#endif
    if ((aux_m == 0) || (aux_dim < nr)) {
      memarray->free(aux_m); aux_m = 0;
      aux_dim = Integer(memarray->get(nr, aux_m));
      aux_task = qr;
    }

    Integer j;
    for (j = 0; j < nc; j++) {
      // Matrix v(house(*this,j,j,tol));
      // rowhouse(*this,v,j,j);
      // mat_xey(nr-j-1,m+j*nr+j+1,v.m+j+1);

      //--- compute Housholder for pivot column
      //--- Matrix v(house(*this,r,r,tol));
      Real* base = m + j * nr + j;
      Real tau = mat_ip(nr - j, base);
      if (tau < tol) {
        return j + 1;
      }
      Real gamma = *(base);           //beta=x(0)
      Real lam = gamma;
      if (gamma < 0.)
        lam -= std::sqrt(tau);
      else
        lam += std::sqrt(tau);
      mat_xmultea(nr - j - 1, base + 1, 1. / lam);
      Real beta = -2. / (1 + (tau - gamma * gamma) / (lam * lam));
      *(aux_m + j) = beta;
      //--- rowhouse(*this,v,j,j);
      (*base) = (lam > 0) ? -std::sqrt(tau) : std::sqrt(tau);
      for (Integer jj = j + 1; jj < nc; jj++) {
        Real d = beta * (*(m + jj * nr + j));
        d += beta * mat_ip(nr - j - 1, m + jj * nr + j + 1, base + 1);
        *(m + jj * nr + j) += d;
        mat_xpeya(nr - j - 1, m + jj * nr + j + 1, base + 1, d);
      }
    }
    return 0;
  }

  int Matrix::QR_factor(Matrix& Q, Real tol) {
    QR_factor(tol);
    Matrix v(nr, 1, 1.);
    Q = diag(v);
    Integer j;
    mat_xea(nr, v.m, 0.);
    for (j = nc - 1; j >= 0; j--) {
      v(j) = 1.;
      mat_xey(nr - j - 1, v.m + j + 1, m + j * nr + j + 1);
      rowhouse(Q, v, j, j);
      mat_xea(nr - j - 1, m + j * nr + j + 1, 0.);
    }
    return 0;
  }

  int Matrix::QR_factor(Indexmatrix& piv, Real tol) {
    chk_init(*this);
    piv.init(Range(0, nc - 1));

    if ((aux_m == 0) || (aux_dim < nr)) {
      memarray->free(aux_m); aux_m = 0;
      aux_dim = Integer(memarray->get(nr, aux_m));
      aux_task = qr;
    }

    Real tau = min_Real; //squared norm of pivot column   
    Integer k = -1;      //index for pivot element
    Matrix c(nc, 1);    //norm squared of remaining columns for selecting pivot
    chk_set_init(c, 1);
    Real nrmbnd = 0.;
    for (Integer j = 0; j < nc; j++) {
      Real d = mat_ip(nr, m + j * nr);
      nrmbnd += d;
      c(j) = d;
      if (d > tau) {
        tau = d;
        k = j;
      }
    }
    nrmbnd *= tol * tol;
    Integer r = 0;       //rank
    while ((tau >= tol) && (r < nc) && (r < nr)) {
      //--- swap pivot if required
      if (r != k) {
        Integer j = piv(k); piv(k) = piv(r); piv(r) = j;
        mat_swap(nr, m + r * nr, m + k * nr);
        c(k) = c(r);
      }
      //--- compute Housholder for pivot column
      //--- Matrix v(house(*this,r,r,tol));
      Real* base = m + r * nr + r;
      Real gamma = *(base);           //beta=x(0)
      Real lam = gamma;
      if (gamma < 0.)
        lam -= std::sqrt(tau);
      else
        lam += std::sqrt(tau);
      mat_xmultea(nr - r - 1, base + 1, 1. / lam);
      Real beta = -2. / (1 + (tau - gamma * gamma) / (lam * lam));
      *(aux_m + r) = beta;
      //--- rowhouse(*this,v,r,r);
      (*base) = (lam > 0) ? -std::sqrt(tau) : std::sqrt(tau);
      for (Integer jj = r + 1; jj < nc; jj++) {
        Real d = beta * (*(m + jj * nr + r));
        d += beta * mat_ip(nr - r - 1, m + jj * nr + r + 1, base + 1);
        *(m + jj * nr + r) += d;
        mat_xpeya(nr - r - 1, m + jj * nr + r + 1, base + 1, d);
      }
      //---
      r++;
      //--- update c
      const Real* mp = m + r * nr + r - 1;
      Real* cp = c.m + r;
      Real nrm = 0.;
      tau = min_Real;
      for (Integer j = r; j < nc; j++, mp += nr) {
        Real d = (*cp++) -= (*mp) * (*mp);
        nrm += d;
        if (d > tau) {
          k = j;
          tau = d;
        }
      }
      if (nrm < nrmbnd)
        break;
    }
    return r;
  }

  int Matrix::QR_factor_relpiv(Indexmatrix& piv, Real tol) {
    assert(tol > 0.);
    chk_init(*this);
    aux_task = none;
    piv.init(Range(0, nc - 1));
    Integer j;
    Integer r = 0;       //rank
    Matrix c(nc, 1);    //norm squared of remaining columns for selecting pivot
    chk_set_init(c, 1);
    Matrix d(nc, 1);    //original norm squared of remaining columns for selecting pivot
    chk_set_init(d, 1);
    for (j = 0; j < nc; j++) {
      d(j) = c(j) = mat_ip(nr, m + j * nr);
    }
    Integer k;         //index for pivot element
    Real* mp;          //pointer to m (this) 
    Real* cp;          //pointer to c
    Real* dp;          //pointer to d
    while ((r < nc) && (r < nr)) {
      //--- find correct pivot and swap
      cp = c.m + r;
      dp = d.m + r;
      Real red = (*dp >= tol) ? *cp / *dp : 0.;
      dp++;
      Real tau = *cp++;
      k = r;
      for (j = r + 1; j < nc; j++, cp++, dp++) {
        Real redj = (*dp >= tol) ? *cp / *dp : 0.;
        if ((redj < tol) || (tau >= *cp)) continue;
        tau = *cp;
        red = redj;
        k = j;
      }
      if ((red < tol) || (tau < tol)) break;
      if (r != k) {
        j = piv(k); piv(k) = piv(r); piv(r) = j;
        mat_swap(nr, m + r * nr, m + k * nr);
        c(k) = c(r);
        d(k) = d(r);
      }
      //--- compute Housholder for pivot column
      Matrix v(house(*this, r, r, tol));
      rowhouse(*this, v, r, r);
      mat_xey(nr - r - 1, m + r * nr + r + 1, v.m + r + 1);
      r++;
      //--- update c
      mp = m + r * nr + r - 1;
      cp = c.m + r;
      for (j = nc - r; --j >= 0; mp += nr) (*cp++) -= (*mp) * (*mp);
    }
    return r;
  }

  int Matrix::QR_factor(Matrix& Q, Indexmatrix& piv, Real tol) {
    Integer r = QR_factor(piv, tol);
    Matrix v(nr, 1, 1.);
    Q = diag(v);
    Integer j;
    mat_xea(nr, v.m, 0.);
    for (j = r - 1; j >= 0; j--) {
      v(j) = 1.;
      mat_xey(nr - j - 1, v.m + j + 1, m + j * nr + j + 1);
      rowhouse(Q, v, j, j);
      mat_xea(nr - j - 1, m + j * nr + j + 1, 0.);
    }
    return r;
  }

  int Matrix::Qt_times(Matrix& A, Integer r) const {
    chk_init(A);
    assert(A.nr == nr);
    if (aux_task == qr) {
      for (Integer j = 0; j < r; j++) {
        Real beta = aux_m[j];
        Real* vp = m + j * nr + j + 1;
        for (Integer k = 0; k < A.nc; k++) {
          Real d = beta * (*(A.m + k * nr + j));
          d += beta * mat_ip(nr - j - 1, A.m + k * nr + j + 1, vp);
          *(A.m + k * nr + j) += d;
          mat_xpeya(nr - j - 1, A.m + k * nr + j + 1, vp, d);
        }
      }
    } else {
      Matrix v(nr, 1);
      chk_set_init(v, 1);
      for (Integer j = 0; j < r; j++) {
        v(j) = 1.;
        mat_xey(nr - j - 1, v.m + j + 1, m + j * nr + j + 1);
        rowhouse(A, v, j, 0);
      }
    }
    return 0;
  }

  int Matrix::Q_times(Matrix& A, Integer r) const {
    chk_init(A);
    assert(A.nr == nr);
    if (aux_task == qr) {
      for (Integer j = r; --j >= 0;) {
        Real beta = aux_m[j];
        Real* vp = m + j * nr + j + 1;
        for (Integer k = 0; k < A.nc; k++) {
          Real d = beta * (*(A.m + k * nr + j));
          d += beta * mat_ip(nr - j - 1, A.m + k * nr + j + 1, vp);
          *(A.m + k * nr + j) += d;
          mat_xpeya(nr - j - 1, A.m + k * nr + j + 1, vp, d);
        }
      }
    } else {
      Matrix v(nr, 1);
      chk_set_init(v, 1);
      for (Integer j = r; --j >= 0;) {
        v(j) = 1.;
        mat_xey(nr - j - 1, v.m + j + 1, m + j * nr + j + 1);
        rowhouse(A, v, j, 0);  //(A,v,i,j)
      }
    }
    return 0;
  }

  int Matrix::times_Q(Matrix& A, Integer r) const {
    chk_init(A);
    Matrix v(nr, 1);
    chk_set_init(v, 1);
    for (Integer j = 0; j < r; j++) {
      v(j) = 1.;
      mat_xey(nr - j - 1, v.m + j + 1, m + j * nr + j + 1);
      colhouse(A, v, j, 0);
    }
    return 0;
  }

  int Matrix::QR_solve(Matrix& rhs, Real tol) {
    chk_init(*this);
    chk_init(rhs);
    //Matrix orig(*this);
    Indexmatrix piv;
    Integer r = QR_factor(piv);
    //orig=orig.cols(piv);
    Qt_times(rhs, r);
    Matrix X(nc, rhs.nc);
    for (Integer j = 0; j < rhs.nc; j++) {
      mat_xey(r, X.m + j * X.nr, rhs.m + j * X.nr);
      mat_xea(nc - r, X.m + j * X.nr + r, 0.);
    }
    chk_set_init(X, 1);
    if (r == nc) {
      int err;
      if ((err = triu_solve(X, tol))) return err;
    } else {
      //--- construct transposed upper triangle
      Matrix A(nc, r);
      for (Integer j = 0; j < r; j++) {
        mat_xea(j, A.m + j * A.nr, 0.);
        mat_xey(A.nr - j, A.m + j * A.nr + j, 1, m + j * nr + j, nr);
      }
      chk_set_init(A, 1);
      //--- full rank QR factorization for this
      A.QR_factor();

      //--- upper triangle is now transposed of intended system
      for (Integer i = 0; i < r; i++) {  //solve for variable row i
        Real d = A(i, i);
        if (::fabs(d) < tol) return i + 1;
        Real* Xp = X.m + i;
        for (Integer j = 0; j < X.nc; j++, Xp += X.nr) {
          (*Xp) /= d;
          mat_xpeya(r - i - 1, Xp + 1, 1, A.m + (i + 1) * A.nr + i, A.nr, -(*Xp));
        }
      }
      A.Q_times(X, r);

    }
    rhs.newsize(X.nr, X.nc);
    for (Integer i = 0; i < X.nr; i++) {
      mat_xey(X.nc, rhs.m + piv(i), X.nc, X.m + i, X.nc);
    }
    chk_set_init(rhs, 1);
    return 0;
  }

  int Matrix::QR_concat_right(const Matrix& A, Indexmatrix& piv, Integer r, Real tol) {
    chk_mult(*this, A, 1, 0);
    chk_init(piv);
#if (CONICBUNDLE_DEBUG>=1)
    if ((piv.rowdim() != nc) || (piv.coldim() != 1) || ((nc > 0) && ((max(piv) != nc - 1) || (min(piv) != 0))))
      MEmessage(MatrixError(ME_unspec, "Matrix::QR_concat_right(...): piv is not a column permutation vector of *this", MTmatrix));
    if ((r < 0) || (r > min(nc, nr)))
      MEmessage(MatrixError(ME_unspec, "Matrix::QR_concat_right(...): r is not the rank of the current QR factorization of *this", MTmatrix));
#endif
    //concatenate the new columns
    piv.concat_below(Range(nc, nc + A.nc - 1));
    Integer oldnc = nc;
    this->concat_right(A);

    if (aux_task == qr) {

      //transform the new columns by the Householder vectors of the r first columns and collect the normbnd
      Real nrmbnd = mat_ip(A.nr * A.nc, A.m);
      for (Integer j = 0; j < r; j++) {
        nrmbnd += mat_ip(j + 1, m + nr * j);
        const Real* base = m + nr * j + j + 1;
        Real beta = aux_m[j];
        for (Integer k = oldnc; k < nc; k++) {
          Real d = beta * (*(m + k * nr + j));
          d += beta * mat_ip(nr - j - 1, m + k * nr + j + 1, base);
          *(m + k * nr + j) += d;
          mat_xpeya(nr - j - 1, m + k * nr + j + 1, base, d);
        }
      }

      //if the old rank r is smaller than oldnc, move the old "zero" columns to the back
      for (Integer i = r; i < oldnc; i++) {
        Integer j = piv(i); piv(i) = piv(i + A.nc); piv(i + A.nc) = j;
        mat_swap(nr, m + nr * i, m + nr * (i + A.nc));
      }
      oldnc = r;

      //now start pivoting with respect to the new columns
      Real tau = min_Real; //squared norm of pivot column   
      Integer k = -1;         //index for pivot element
      Matrix c(A.nc, 1);    //norm squared of remaining columns for selecting pivot
      chk_set_init(c, 1);
      for (Integer j = 0; j < A.nc; j++) {
        Real d = mat_ip(nr - r, m + (r + j) * nr + r);
        c(j) = d;
        if (d > tau) {
          tau = d;
          k = r + j;
        }
      }
      nrmbnd *= tol * tol;
      while ((tau >= tol) && (r < oldnc + A.nc) && (r < nr)) {
        //--- swap pivot if required
        if (r != k) {
          Integer j = piv(k); piv(k) = piv(r); piv(r) = j;
          mat_swap(nr, m + r * nr, m + k * nr);
          c(k - oldnc) = c(r - oldnc);
        }
        //--- compute Housholder for pivot column
        //--- Matrix v(house(*this,r,r,tol));
        Real* base = m + r * nr + r;
        Real gamma = *(base);           //beta=x(0)
        Real lam = gamma;
        if (lam < 0.)
          lam -= std::sqrt(tau);
        else
          lam += std::sqrt(tau);
        mat_xmultea(nr - r - 1, base + 1, 1. / lam);
        Real beta = -2. / (1 + (tau - gamma * gamma) / (lam * lam));
        *(aux_m + r) = beta;
        //--- rowhouse(*this,v,r,r);
        (*base) = (lam > 0) ? -std::sqrt(tau) : std::sqrt(tau);
        for (Integer k = r + 1; k < nc; k++) {
          Real d = beta * (*(m + k * nr + r));
          d += beta * mat_ip(nr - r - 1, m + k * nr + r + 1, base + 1);
          *(m + k * nr + r) += d;
          mat_xpeya(nr - r - 1, m + k * nr + r + 1, base + 1, d);
        }
        //---
        r++;
        //--- update c
        const Real* mp = m + r * nr + r - 1;
        Real* cp = c.m + r - oldnc;
        Real nrm = 0.;
        tau = min_Real;
        for (Integer j = r; j < oldnc + A.nc; j++, mp += nr) {
          Real d = (*cp++) -= (*mp) * (*mp);
          nrm += d;
          if (d > tau) {
            k = j;
            tau = d;
          }
        }
        if (nrm < nrmbnd)
          break;
      }
    }

    else {
      //transform the new columns by the Householder vectors of the r first columns
      Matrix v(nr, 1); chk_set_init(v, 1);
      Real* vp = v.m;
      for (Integer j = 0; j < r; j++, vp++) {
        *vp = 1.;
        mat_xey(nr - j - 1, vp + 1, m + j * nr + j + 1);
        rowhouse(*this, v, j, oldnc);
        *vp = 0.;
      }

      //if the old rank r is smaller than oldnc, move the old "zero" columns to the back
      for (Integer i = r; i < oldnc; i++) {
        Integer j = piv(i); piv(i) = piv(i + A.nc); piv(i + A.nc) = j;
        mat_swap(nr, m + nr * i, m + nr * (i + A.nc));
      }
      oldnc = r;

      //now start pivoting with respect to the new columns
      Matrix c(A.nc, 1);    //norm squared of remaining columns for selecting pivot
      for (Integer j = 0; j < A.nc; j++) {
        c(j) = mat_ip(nr - r, m + (j + r) * nr + r);
      }
      Integer k;         //index for pivot element
      Real* mp;          //pointer to m (this) 
      Real* cp;          //pointer to c
      Integer pivotnc = oldnc + A.nc;
      while ((r < pivotnc) && (r < nr)) {
        //--- find correct pivot and swap
        cp = c.m + r - oldnc;
        Real tau = *cp++;
        k = r;
        for (Integer j = r + 1; j < pivotnc; j++, cp++) {
          if (tau >= *cp) continue;
          tau = *cp;
          k = j;
        }
        if (tau < tol) break;

        if (r != k) {
          Integer j = piv(k); piv(k) = piv(r); piv(r) = j;
          mat_swap(nr, m + r * nr, m + k * nr);
          c(k - oldnc) = c(r - oldnc);
        }
        //--- compute Housholder for pivot column
        Matrix v(house(*this, r, r, tol));
        rowhouse(*this, v, r, r);
        mat_xey(nr - r - 1, m + r * nr + r + 1, v.m + r + 1);
        r++;
        //--- update c
        mp = m + r * nr + r - 1;
        cp = c.m + r - oldnc;
        for (Integer j = pivotnc - r; --j >= 0; mp += nr)
          (*cp++) -= (*mp) * (*mp);
      }
    }

    return r;

  }


  // *************************************************************************
  //                              house
  // *************************************************************************

  //GvL 2ndEd.

  Matrix house(const Matrix& x, Integer i, Integer j, Real tol) {
    chk_init(x);
#if (CONICBUNDLE_DEBUG>=1)
    if ((i < 0) || (j < 0) || (x.nr <= i) || (x.nc <= j))
      MEmessage(MEdim(x.nr, x.nc, j, j, "house: x and j do not match", MTmatrix));
#endif
    Real* xbase = x.m + j * x.nr + i;
    Real mu = ::sqrt(mat_ip(x.nr - i, xbase));
    Real beta = 1.;
    if (mu > tol) {
      beta = *(xbase);           //beta=x(0)
      if (beta < 0) beta -= mu;
      else beta += mu;
    }
    Matrix v(x.nr, 1);
    mat_xea(i, v.m, 0.);
    *(v.m + i) = 1.;
    mat_xeya(x.nr - i - 1, v.m + i + 1, xbase + 1, 1. / beta);
    chk_set_init(v, 1);
    return v;
  }

  // *************************************************************************
  //                              rowhouse
  // *************************************************************************

  //Housholder premultiplication of A with Householder-vector v (GvL 2nd p197)
  //Integer i: gives the first non zero entry of v
  //Integer j: transformtaion applied to columns j to A.coldim()-1 of A 

  int rowhouse(Matrix& A, const Matrix& v, Integer i, Integer j) {
    chk_init(v);
    chk_init(A);
#if (CONICBUNDLE_DEBUG>=1)
    if ((A.nr != v.nr) || (v.nc != 1))
      MEmessage(MEdim(A.nr, A.nc, v.nr, v.nc, "rowhouse: A and v do not match", MTmatrix));
    if ((i < 0) || (j < 0) || (i >= A.nr))
      MEmessage(MEdim(A.nr, A.nc, i, j, "rowhouse: A and (i,j) do not match", MTmatrix));
#endif
    Integer m, n;
    A.dim(m, n);
    Real beta = -2. / mat_ip(m - i, v.m + i);
    Real d;
    Integer k;
    //w=beta*transpose(A(i:m-1,j:n-1))*v(i:m-1);
    //A(i:m-1,j:n-1)+=v(i:m-1)*transpose(w);
    for (k = j; k < n; k++) {
      d = beta * mat_ip(m - i, A.m + k * m + i, v.m + i);
      mat_xpeya(m - i, A.m + k * m + i, v.m + i, d);
    }
    return 0;
  }

  // *************************************************************************
  //                              colhouse
  // *************************************************************************

  int colhouse(Matrix& A, const Matrix& v, Integer i, Integer j) {
    chk_init(v);
    chk_init(A);
#if (CONICBUNDLE_DEBUG>=1)
    if ((A.nc != v.nr) || (v.nc != 1))
      MEmessage(MEdim(A.nr, A.nc, v.nr, v.nc, "colhouse: A and v do not match", MTmatrix));
    if ((i < 0) || (j < 0) || (i >= A.nc) || (j >= A.nr))
      MEmessage(MEdim(A.nr, A.nc, i, j, "colhouse: A and (i,j) do not match", MTmatrix));
#endif
    Integer m, n;
    A.dim(m, n);
    Real beta = -2. / mat_ip(n - i, v.m + i);
    Matrix w(m, 1);
    chk_set_init(w, 1);
    Integer k;
    //w=beta*A(j:m-1,i:n-1)*v(i:n-1);
    for (k = j; k < m; k++)
      w(k) = beta * mat_ip(n - i, A.m + i * m + k, m, v.m + i, 1);
    //A(j:m-1,i:n-1)+=w*transpose(v);
    for (k = i; k < n; k++)
      mat_xpeya(m - j, A.m + k * m + j, w.m + j, v(k));
    return 0;
  }

}

