/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  Matrix/psqmr.cxx
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



#include <cmath>
#include <limits>
#include "mymath.hxx"
#include "psqmr.hxx"


using namespace CH_Tools;

namespace CH_Matrix_Classes {

  // *************************************************************************
  //                              constructor(s)
  // *************************************************************************

  Psqmr::Psqmr(std::ostream* out, int pril) :
    myout(out), print_level(pril) {
    maxit = -1;
    resnorm = 0.;
    avg_reduction = -1.;
    nmult = 0;
    err = 0;
  }

  // *************************************************************************
  //                              compute
  // *************************************************************************

  //this closely follows the paper 
  //R. W. Freund and N. M. Nachtigal "A new Krylov-subspace method for symmetric 
  //indefinite linear system"
  //and an implementation of it in QSDP by K.C.Toh

  int Psqmr::compute(IterativeSystemObject& sys,
    Matrix& x,
    Real in_termprec,
    Matrix* storex,
    Integer storestep) {
    if ((x.dim() > 0) && (x.dim() != (sys.ItSys_rhs()).dim())) {
      if (myout)
        (*myout) << "**** ERROR Psqmr::compute(...): dimensions of input vector and system do not match" << std::endl;
      err = 1;
      return err;
    }
    avg_reduction = 0.1;     //never mind the initial value, it will be ignored
    nmult = 0;
    termprec = in_termprec;
    Matrix Sq;
    Matrix r(sys.ItSys_rhs());
    if (x.dim() == 0) {
      x.init((sys.ItSys_rhs()).dim(), 1, 0.);
      Sq.init((sys.ItSys_rhs()).dim(), 1, 0.);
    } else {
      if (sys.ItSys_mult(x, Sq)) {
        if (myout)
          (*myout) << "****ERROR Psqmr::compute(...): matrix-vector multiplication subroutine returned an error" << std::endl;
        err = 2;
        if ((myout) && (print_level >= 1))
          (*myout) << "PSQMR " << std::setw(2) << nmult << "(" << err << "): " << resnorm << std::endl;
        return err;
      }
      nmult++;
      r -= Sq;
    }
    Matrix res(r);
    old_resnorm = resnorm = norm2(res);
    avg_reduction = 0.1;     //never mind the initial value, it will be ignored
    if (resnorm <= termprec) {
      err = 0;
      return err;
    }

    Matrix q(r);
    if (sys.ItSys_precondM1(q)) {
      if (myout)
        (*myout) << "****ERROR Psqmr::compute(...): preconditioner subroutine returned an error" << std::endl;
      err = 2;
      if ((myout) && (print_level >= 1))
        (*myout) << "PSQMR " << std::setw(2) << nmult << "(" << err << "): " << resnorm << std::endl;
      return err;
    }
    Real old_tau = norm2(q);
    if (sys.ItSys_precondM2(q)) {
      if (myout)
        (*myout) << "****ERROR Psqmr::compute(...): preconditioner subroutine returned an error" << std::endl;
      err = 2;
      if ((myout) && (print_level >= 1))
        (*myout) << "PSQMR " << std::setw(2) << nmult << "(" << err << "): " << resnorm << std::endl;
      return err;
    }
    Real old_rho = ip(r, q);
    Real old_theta = 0;
    Matrix d(x.dim(), 1, 0.);
    Matrix Sd(x.dim(), 1, 0.);
    Matrix u;

    //============ main loop
    Integer maxmult = max(30, x.dim());
    if (maxit > 0) {
      maxmult = maxit;
    }
    while (nmult < maxmult) {

      //--- step 1) matrix-vector multiplication

      if (sys.ItSys_mult(q, Sq)) {
        if (myout)
          (*myout) << "****ERROR Psqmr::compute(...): matrix-vector multiplication subroutine returned an error" << std::endl;
        err = 2;
        if ((myout) && (print_level >= 1))
          (*myout) << "_PSQMR " << std::setw(2) << nmult << "(" << err << "): " << resnorm << std::endl;
        return err;
      }
      nmult++;
      Real sigma = ip(q, Sq);
      if (abs(sigma) < 1e-10 * abs(old_rho)) {
        err = 0;
        if ((myout) && (print_level >= 1))
          (*myout) << "_PSQMR " << std::setw(2) << nmult << "(" << err << "): " << resnorm << std::endl;
        return err;
      }
      Real alpha = old_rho / sigma;
      r.xpeya(Sq, -alpha);

      //--- step 2) precondition by M1 and compute new x

      u.init(r);
      if (sys.ItSys_precondM1(u)) {
        if (myout)
          (*myout) << "****ERROR Psqmr::compute(...): preconditioner subroutine returned an error" << std::endl;
        err = 2;
        if ((myout) && (print_level >= 1))
          (*myout) << "PSQMR " << std::setw(2) << nmult << "(" << err << "): " << resnorm << std::endl;
        return err;
      }

      Real theta = norm2(u) / old_tau;
      Real c = 1. / std::sqrt(1 + theta * theta);
      Real tau = old_tau * theta * c;
      Real gamma = c * c * old_theta * old_theta;
      Real eta = c * c * alpha;

      xbpeya(d, q, eta, gamma);  //d=d*gamma+q*beta
      x += d;
      if (storex) {
        if (storex->coldim() == 0)
          storex->concat_right(x);
        else if ((storestep > 0) && (nmult % storestep == 0))
          storex->concat_right(x);
      }
      // Matrix tmpvec;sys.ItSys_mult(x,tmpvec);   //TEST
      // std::cout<<" psqmrviolx="<<norm2(sys.ItSys_rhs()-tmpvec)<<" normd="<<norm2(d)<<std::endl ; //TEST
      xbpeya(Sd, Sq, eta, gamma);  //Sd=Sd*gamma+Sq*beta
      res -= Sd;
      old_resnorm = resnorm;
      resnorm = norm2(res);
      //for the avg_reduction we start once nmult==2, so its initial value is ignored
      if ((nmult >= 2) && (nmult < 11)) {
        avg_reduction = (avg_reduction * (nmult - 2) + resnorm / old_resnorm) / (nmult - 1);
      }
      if (nmult >= 11) {
        avg_reduction = avg_reduction * 9 / 10. + resnorm / old_resnorm / 10.;
      }

      if ((myout) && (print_level >= 2))
        (*myout) << "PSQMR " << std::setw(2) << nmult << ": " << resnorm << " red=" << resnorm / old_resnorm << " avg_red : " << avg_reduction << std::endl;
      if ((myout) && (print_level >= 3))
        (*myout) << "res : " << res << std::endl;
      if (resnorm <= termprec) {
        err = 0;
        if ((myout) && (print_level >= 1))
          (*myout) << "_PSQMR " << std::setw(2) << nmult << "(" << err << "): " << resnorm << std::endl;
        return err;
      }

      //--- step 2) precondition by M2 

      if (sys.ItSys_precondM2(u)) {
        if (myout)
          (*myout) << "****ERROR Psqmr::compute(...): preconditioner subroutine returned an error" << std::endl;
        err = 2;
        if ((myout) && (print_level >= 1))
          (*myout) << "_PSQMR " << std::setw(2) << nmult << "(" << err << "): " << resnorm << std::endl;
        return err;
      }
      Real rho = ip(r, u);
      if (abs(old_rho) < 1e-10 * abs(rho)) {
        err = 0;
        if ((myout) && (print_level >= 1))
          (*myout) << "_PSQMR " << std::setw(2) << nmult << "(" << err << "): " << resnorm << std::endl;
        return err;
      }
      Real beta = rho / old_rho;
      q *= beta; q += u;
      old_tau = tau;
      old_rho = rho;
      old_theta = theta;
    }

    err = 3;
    if ((myout) && (print_level >= 1))
      (*myout) << "_PSQMR " << std::setw(2) << nmult << "(" << err << "): " << resnorm << std::endl;
    return err;
  }

}
