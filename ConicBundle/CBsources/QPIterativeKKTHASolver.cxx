/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/QPIterativeKKTHASolver.cxx
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


#include "pcg.hxx"
#include "QPIterativeKKTHASolver.hxx"

using namespace CH_Matrix_Classes;

namespace ConicBundle {

  QPIterativeKKTHASolver::~QPIterativeKKTHASolver() {
  }

  // *************************************************************************
  //                             QPsolve_KKTsystem
  // *************************************************************************

  // solve the KKT System

  /// on input
  /// dualrhs = -(Qx+A'y+G'modelx+c) +mu (...)
  /// primalrhs = -(Ax+s) + mu ()...
  int QPIterativeKKTHASolver::QPsolve_KKTsystem(Matrix& solx, Matrix& soly,
    const Matrix& primalrhs,
    const Matrix& dualrhs,
    Real rhsmu,
    Real rhscorr,
    Real prec,
    QPSolverParameters* /* params */) {
    maxit_bnd = max(Integer(100 * std::log(primalrhs.rowdim() + dualrhs.rowdim())), min(maxit_bnd, 2 * (primalrhs.rowdim() + dualrhs.rowdim())));

    //---- prepare the right hand side
    sysrhs.init(dualrhs);
    int status1 = 0;
    if (model) {
      status1 = model->add_Schur_rhs(sysrhs, 0, rhsmu, rhscorr);
      if (status1) {
        if (cb_out()) {
          get_out() << "**** WARNING: QPIterativeKKTHASolver::QPsolve_KKTsystem(....): model->add_Schur_rhs failed and returned " << status1 << std::endl;
        }
      }
    }
    sysrhs.concat_below(primalrhs);

    //---- call the solver
    assert(solver);
    Real termprec = prec * min(1., norm2(sysrhs));
    Integer ntries = 5;
    int status2;
    do {
      solver->set_maxit(maxit_bnd);
      status2 = solver->compute(*this, sol, termprec);
      maxit_bnd = max(maxit_bnd, 2 * solver->get_maxit());
      if ((solver->get_residual_norm() > 1.1 * solver->get_termprec()) && (solver->get_nmult() >= solver->get_maxit())) {
        if (termprec < 1e-7 * min(1., norm2(sysrhs)))
          termprec *= 10.;
        Matrix tmpvec;
        tmpvec.rand(sol.rowdim(), 1);
        tmpvec *= (0.1 * solver->get_residual_norm() / norm2(tmpvec));
        sol += tmpvec;
      }
    } while (
      (solver->get_residual_norm() > 1.1 * solver->get_termprec())
      && (solver->get_nmult() >= solver->get_maxit())
      && (--ntries > 0)
      );
    if (status2) {
      if (cb_out()) {
        get_out() << "**** WARNING: QPIterativeKKTHASolver::QPsolve_KKTsystem(....): solver->compute failed and returned " << status2 << std::endl;
      }
    }

    //---- split up the solution
    solx.init(dualrhs.rowdim(), 1, sol.get_store());
    soly.init(primalrhs.rowdim(), 1, sol.get_store() + dualrhs.rowdim());

    int status3 = 0;
    if (model) {
      status3 = model->compute_step(solx);
      if (status3) {
        if (cb_out()) {
          get_out() << "**** WARNING: QPIterativeKKTHASolver::QPsolve_KKTsystem(....): model->compute_step failed and returned " << status3 << std::endl;
        }
      }
    }

    return std::abs(status1) + std::abs(status2) + std::abs(status3);
  }

  // *************************************************************************
  //                             ItSys_mult
  // *************************************************************************

  int QPIterativeKKTHASolver::ItSys_mult(const Matrix& in_vec, Matrix& out_vec) {
    nmult++;

    out_vec.newsize(in_vec.rowdim(), 1);
    in_vecx.init(KKTdiagx.rowdim(), 1, in_vec.get_store());

    out_vec = in_vecx;
    out_vec %= KKTdiagx;
    Hp->add_Hx(in_vecx, out_vec, Hfactor);

    if (model) {
      model->add_Schur_mult(in_vecx, out_vec);
    }

    if (KKTdiagy.rowdim() > 0) {
      in_vecy.init(KKTdiagy.rowdim(), 1, in_vec.get_store() + KKTdiagx.rowdim());
      genmult(*A, in_vecy, out_vec, 1., 1., 1);
      in_vecy %= KKTdiagy;
      genmult(*A, in_vecx, in_vecy, 1., -1.);
      out_vec.concat_below(in_vecy);
    }

    Real xnormsqr = mat_ip(KKTdiagx.rowdim(), in_vec.get_store());
    if (xnormsqr > 1e-10) {
      Real Ritz = mat_ip(KKTdiagx.rowdim(), in_vec.get_store(), out_vec.get_store());
      Ritz /= xnormsqr;
      if (Ritz > blockH_norm)
        blockH_norm = Ritz;
    }

    return 0;
  }


}

