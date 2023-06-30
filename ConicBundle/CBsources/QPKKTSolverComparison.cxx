/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/QPKKTSolverComparison.cxx
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



#include <iomanip>
#include <sstream>
#include <fstream>
#include "QPKKTSolverComparison.hxx"

using namespace CH_Matrix_Classes;
using namespace CH_Tools;

namespace ConicBundle {

  QPKKTSolverComparison::QPKKTSolverComparison(CBout* cb, int cbinc) :
    QPKKTSolverObject() {
    //testmodel=0;
    testHp = 0;
    clear();
    set_cbout(cb, cbinc);
  }

  QPKKTSolverComparison::~QPKKTSolverComparison() {
    clear();
  }

  void QPKKTSolverComparison::clear() {
    for (unsigned int i = 0; i < solver.size(); i++) {
      delete solver[i];
      if ((i > 0) && (model[i] != 0)) {
        model[i]->recursive_delete_and_clear();
        delete model[i];
      }
    }
    //if (testmodel)
    //  testmodel->recursive_delete_and_clear();
    //delete testmodel;
    //testmodel=0;
    solver.clear();
    model.clear();
    probdata.clear();
  }

  int QPKKTSolverComparison::violation(Matrix& violvec, const Matrix& solx, const Matrix& soly, QPModelBlockObject* inmodel, QPKKT_SolverStats& stats) {


    //first block of rows belonging to H
    Matrix primalviol = solx;
    primalviol %= testKKTdiagx;
    testHp->add_Hx(solx, primalviol, testHfactor);

    if (testA) {
      genmult(*testA, soly, primalviol, 1., 1., 1);
    }

    if (inmodel) {
      inmodel->B_times(inmodel->get_dx(), primalviol, 1., 1., 1, 0);
    }
    primalviol -= testdualrhs;
    stats.Hviol = norm2(primalviol);
    violvec = primalviol;

    //second block of rows belonging to the constraints on the design variables
    Matrix dualviol = soly;
    if (testA) {
      dualviol %= testKKTdiagy;
      genmult(*testA, solx, dualviol, 1., -1.);
    }
    dualviol -= testprimalrhs;
    stats.Aviol = norm2(dualviol);
    violvec.concat_below(dualviol);

    //third block of rows belonging to the model block
    Matrix modelviol(0, 1, 0.);
    Matrix modelcviol(0, 1, 0.);
    if (inmodel) {
      modelviol = inmodel->get_sysviol_model(solx);
      if ((cb_out(2)) && (norm2(modelviol) > 1e-8))
        get_out() << " modelviol=" << modelviol;
      modelcviol = inmodel->get_sysviol_constraints();
      if ((cb_out(2)) && (norm2(modelcviol) > 1e-8))
        get_out() << " modelcviol=" << modelcviol;
    }
    stats.Bviol = norm2(modelviol);
    stats.Cviol = norm2(modelcviol);
    violvec.concat_below(modelviol);
    violvec.concat_below(modelcviol);
    stats.sysviol = norm2(violvec);

    if (cb_out(3))
      get_out() << " n2(viol)=" << norm2(violvec) << "[" << norm2(primalviol) << "," << norm2(dualviol) << "," << norm2(modelviol) << "," << norm2(modelcviol) << "]\n";

    return 0;
  }

  int QPKKTSolverComparison::add_solver(QPKKTSolverObject* solverp, const char* name) {
    if (solverp == 0)
      return 0;
    solver.push_back(solverp);
    solvername.push_back(std::string(name));
    model.push_back(0);
    return 0;
  }

  int QPKKTSolverComparison::QPinit_KKTdata(QPSolverProxObject* Hp,
    QPModelBlockObject* inmodel,
    const Sparsemat* A,
    const Indexmatrix* eq_indices
  ) {
    int err = 0;
    {
      Matrix D;
      const Matrix* Vp;
      Hp->get_precond(D, Vp);
      QPKKT_ProbStats prob;
      prob.Qdim = D.rowdim();
      prob.Vdim = (Vp ? Vp->coldim() : 0);
      prob.Arowdim = (A ? A->rowdim() : 0);
      prob.Aeqdim = (eq_indices ? eq_indices->rowdim() : 0);
      prob.Bdim = (inmodel ? inmodel->dim_model() : 0);
      prob.Cdim = (inmodel ? inmodel->dim_constraints() : 0);
      probdata.push_back(prob);
    }

    QPKKT_ProbStats& prob = probdata[probdata.size() - 1];
    for (unsigned int i = 0; i < solver.size(); i++) {
      if (i == 0) {
        model[0] = inmodel;
      } else {
        if (model[i]) {
          model[i]->recursive_delete_and_clear();
          delete model[i];
        }
        model[i] = (inmodel ? inmodel->clone() : 0);
        assert((inmodel == 0) || (model[i]));
      }
      clock.start();
      int erri = solver[i]->QPinit_KKTdata(Hp, model[i], A, eq_indices);
      prob.inittime.push_back(clock.time());
      if ((erri) && (cb_out())) {
        get_out() << "**** ERROR in QPKKTSolverComparison::QPinit_KKTdata(): solver " << i << " returned " << erri << std::endl;
      }
      if (i == 0)
        err = erri;

    }

    // if (testmodel){
    //   testmodel->recursive_delete_and_clear();
    //   delete testmodel;
    // }
    // testmodel=(inmodel?inmodel->clone():0);
    testHp = Hp;
    testA = A;

    return err;
  }

  int QPKKTSolverComparison::QPinit_KKTsystem(const Matrix& KKTdiagx,
    const Matrix& KKTdiagy,
    Real Hfactor,
    Real prec,
    QPSolverParameters* params) {
    int err = 0;
    for (unsigned int i = 1; i < solver.size(); i++) {
      if (model[i])
        model[i]->recursive_copy_data_of(model[0]);
    }
    // if (testmodel)
    //   testmodel->recursive_copy_data_of(model[0]);

    testHfactor = Hfactor;
    testKKTdiagx = KKTdiagx;
    testKKTdiagy = KKTdiagy;

    QPKKT_ProbStats& prob = probdata[probdata.size() - 1];
    {
      QPKKT_KKTStats kkt;
      kkt.prec = prec;
      kkt.mu = max_Real;
      prob.kktdata.push_back(kkt);
    }

    QPKKT_KKTStats& kkt = prob.kktdata[prob.kktdata.size() - 1];

    for (unsigned int i = 0; i < solver.size(); i++) {
      QPKKT_SolverStats sol;
      clock.start();
      int erri = solver[i]->QPinit_KKTsystem(KKTdiagx, KKTdiagy, Hfactor, prec, params);
      sol.preptime = clock.time();
      sol.prepnmult = solver[i]->QPget_nmult();
      sol.cond = solver[i]->QPget_condition_number();
      sol.rank = solver[i]->QPget_precond_rank();
      kkt.sdata.push_back(sol);
      if ((erri) && (cb_out())) {
        get_out() << "**** ERROR in QPKKTSolverComparison::QPinit_KKTsystem(): solver " << i << " returned " << erri << std::endl;
      }
      if (i == 0)
        err = erri;
    }
    return err;
  }

  int QPKKTSolverComparison::QPsolve_KKTsystem(Matrix& solx,
    Matrix& soly,
    const Matrix& primalrhs,
    const Matrix& dualrhs,
    Real rhsmu,
    Real rhscorr,
    Real prec,
    QPSolverParameters* params) {
    testprimalrhs = primalrhs;
    testdualrhs = dualrhs;

    QPKKT_ProbStats& prob = probdata[probdata.size() - 1];
    QPKKT_KKTStats& kkt = prob.kktdata[prob.kktdata.size() - 1];
    kkt.prec = min(kkt.prec, prec);
    if (rhsmu > 0.)
      kkt.mu = min(kkt.mu, rhsmu);

    int err = 0;
    insolx = solx;
    insoly = soly;
    Real relprec = max(1e-10 * rhsmu, prec) * max(1., max(norm2(primalrhs), norm2(dualrhs)));
    for (unsigned int i = 0; i < solver.size(); i++) {
      QPKKT_SolverStats& sol = kkt.sdata[i];
      solx = insolx;
      soly = insoly;
      if (rhsmu == 0.) {
        sol.prednmult -= solver[i]->QPget_nmult();
      } else {
        sol.corrnmult -= solver[i]->QPget_nmult();
      }
      clock.start();
      int erri = solver[i]->QPsolve_KKTsystem(solx, soly, primalrhs, dualrhs,
        rhsmu, rhscorr, prec, params);
      if (rhsmu == 0.) {
        sol.predtime += clock.time();
        sol.prednmult += solver[i]->QPget_nmult();
      } else {
        sol.corrtime += clock.time();
        sol.corrnmult += solver[i]->QPget_nmult();
      }
      if ((erri) && (cb_out())) {
        get_out() << "**** ERROR in QPKKTSolverComparison::QPsolve_KKTsystem(): solver " << i << " returned " << erri << std::endl;
      }
      if (i == 0) {
        err = erri;
        outsolx = solx;
        outsoly = soly;
        //if (testmodel)
        //  testmodel->compute_step(solx);
      } else {
        if (cb_out() && (norm2(solx - outsolx) > relprec)) {
          get_out() << "**** WARNING in QPKKTSolverComparison::QPsolve_KKTsystem(): solution of solver " << i << " differs a lot, norm2(solx-outsolx)=" << norm2(solx - outsolx) << ">=" << relprec << std::endl;
        }
        if (cb_out(2)) {
          get_out() << " reldiff[" << i << "]=" << norm2(solx - outsolx) / (norm2(outsolx) + 1.);
        }
      }
      Matrix violvec;
      violation(violvec, solx, soly, model[i], sol);
    }
    solx = outsolx;
    soly = outsoly;
    return err;
  }


  Real QPKKTSolverComparison::QPget_blockH_norm() {
    Real retval = 0.;
    for (unsigned int i = 0; i < solver.size(); i++) {
      Real val = solver[i]->QPget_blockH_norm();
      if (i == 0)
        retval = val;
    }
    return retval;
  }


  Real QPKKTSolverComparison::QPget_blockA_norm() {
    Real retval = 0.;
    for (unsigned int i = 0; i < solver.size(); i++) {
      Real val = solver[i]->QPget_blockA_norm();
      if (i == 0)
        retval = val;
    }
    return retval;
  }


  int QPKKTSolverComparison::get_mu_stats(Real lbmu,
    Real ubmu,
    Indexmatrix& dims,
    Matrix& mu,
    Matrix& prepsecs,
    Matrix& predsecs,
    Matrix& corrsecs,
    Indexmatrix& predcalls,
    Indexmatrix& corrcalls,
    Matrix& cond,
    Indexmatrix& pccols,
    Matrix& sysviol) {
    Integer cnt = 0;
    for (unsigned i = 0; i < probdata.size(); i++)
      cnt += probdata[i].cnt_mu_cols(lbmu, ubmu);
    Integer nsolvers = Integer(solvername.size());
    dims.newsize(6, cnt); dims.init(6, 0, Integer(0));
    mu.newsize(1, cnt); mu.init(1, 0, 0.);
    prepsecs.newsize(nsolvers, cnt); prepsecs.init(nsolvers, 0, 0.);
    predsecs.newsize(nsolvers, cnt); predsecs.init(nsolvers, 0, 0.);
    corrsecs.newsize(nsolvers, cnt); corrsecs.init(nsolvers, 0, 0.);
    predcalls.newsize(nsolvers, cnt); predcalls.init(nsolvers, 0, Integer(0));
    corrcalls.newsize(nsolvers, cnt); corrcalls.init(nsolvers, 0, Integer(0));
    cond.newsize(nsolvers, cnt); cond.init(nsolvers, 0, 0.);
    pccols.newsize(nsolvers, cnt); pccols.init(nsolvers, 0, Integer(0));
    sysviol.newsize(nsolvers, cnt); sysviol.init(nsolvers, 0, 0.);

    for (unsigned i = 0; i < probdata.size(); i++)
      probdata[i].get_mu_stats(lbmu, ubmu, dims, mu,
        prepsecs, predsecs, corrsecs,
        predcalls, corrcalls, cond, pccols, sysviol);
    return 0;
  }

  /// return one data column per subproblem (more efficient to append than lines) with the sum of the time/calls/etc. 
  int QPKKTSolverComparison::get_prob_stats(Indexmatrix& dims,
    Indexmatrix& iterations,
    Matrix& lastmu,
    Matrix& prepsecs,
    Matrix& predsecs,
    Matrix& corrsecs,
    Indexmatrix& predcalls,
    Indexmatrix& corrcalls) {
    Integer cnt = Integer(probdata.size());
    Integer nsolvers = Integer(solvername.size());
    dims.newsize(6, cnt); dims.init(6, 0, Integer(0));
    iterations.newsize(1, cnt); iterations.init(1, 0, Integer(0));
    lastmu.newsize(1, cnt); lastmu.init(1, 0, 0.);
    prepsecs.newsize(nsolvers, cnt); prepsecs.init(nsolvers, 0, 0.);
    predsecs.newsize(nsolvers, cnt); predsecs.init(nsolvers, 0, 0.);
    corrsecs.newsize(nsolvers, cnt); corrsecs.init(nsolvers, 0, 0.);
    predcalls.newsize(nsolvers, cnt); predcalls.init(nsolvers, 0, Integer(0));
    corrcalls.newsize(nsolvers, cnt); corrcalls.init(nsolvers, 0, Integer(0));

    for (unsigned i = 0; i < probdata.size(); i++)
      probdata[i].get_prob_stats(dims, iterations, lastmu,
        prepsecs, predsecs, corrsecs,
        predcalls, corrcalls);
    return 0;
  }

  std::ostream& operator<<(std::ostream& out, const QPKKTSolverComparison& q) {
    out << " " << q.solvername.size() << "\n";
    for (unsigned i = 0; i < q.solvername.size(); i++) {
      out << q.solvername[i] << "\n";
    }
    out << q.probdata.size() << "\n";
    for (unsigned i = 0; i < q.probdata.size(); i++) {
      out << q.probdata[i] << "\n";
    }
    return out;
  }

  std::istream& operator>>(std::istream& in, QPKKTSolverComparison& q) {
    int sz;
    in >> sz;
    q.solvername.resize(unsigned(sz));
    for (unsigned i = 0; i < q.solvername.size(); i++) {
      in >> q.solvername[i];
    }
    in >> sz;
    q.probdata.resize(unsigned(sz));
    for (unsigned i = 0; i < q.probdata.size(); i++) {
      in >> q.probdata[i];
    }
    return in;
  }

  std::ostream& operator<<(std::ostream& out, const QPKKT_SolverStats& s) {
    out.precision(12);
    out << " " << std::setw(6) << s.preptime;
    out << " " << std::setw(6) << s.predtime;
    out << " " << std::setw(6) << s.corrtime;
    out << " " << std::setw(5) << s.prepnmult;
    out << " " << std::setw(5) << s.prednmult;
    out << " " << std::setw(5) << s.corrnmult;
    out << " " << std::setw(14) << s.cond;
    out << " " << std::setw(4) << s.rank;
    out << " " << std::setw(15) << s.Hviol;
    out << " " << std::setw(15) << s.Aviol;
    out << " " << std::setw(15) << s.Bviol;
    out << " " << std::setw(15) << s.Cviol;
    out << " " << std::setw(15) << s.sysviol;
    return out;
  }

  std::istream& operator>>(std::istream& in, QPKKT_SolverStats& s) {
    in >> s.preptime;
    in >> s.predtime;
    in >> s.corrtime;
    in >> s.prepnmult;
    in >> s.prednmult;
    in >> s.corrnmult;
    in >> s.cond;
    in >> s.rank;
    in >> s.Hviol;
    in >> s.Aviol;
    in >> s.Bviol;
    in >> s.Cviol;
    in >> s.sysviol;
    return in;
  }

  std::ostream& operator<<(std::ostream& out, const QPKKT_KKTStats& k) {
    out.precision(12);
    out << "\n" << k.prec << " " << k.mu << " " << k.sdata.size();
    for (unsigned i = 0; i < k.sdata.size(); i++)
      out << "\n" << k.sdata[i];
    return out;
  }

  std::istream& operator>>(std::istream& in, QPKKT_KKTStats& k) {
    in >> k.prec >> k.mu;
    int sz;
    in >> sz;
    k.sdata.resize(unsigned(sz));
    for (unsigned i = 0; i < k.sdata.size(); i++)
      in >> k.sdata[i];
    return in;
  }

  std::ostream& operator<<(std::ostream& out, const QPKKT_ProbStats& p) {
    out << " " << p.Qdim << " " << p.Vdim << " " << p.Arowdim << " " << p.Aeqdim;
    out << " " << p.Bdim << " " << p.Cdim << "\n" << p.inittime.size();
    for (unsigned i = 0; i < p.inittime.size(); i++)
      out << " " << p.inittime[i];
    out << "\n" << p.kktdata.size();
    for (unsigned i = 0; i < p.kktdata.size(); i++)
      out << " " << p.kktdata[i];
    return out;
  }

  std::istream& operator>>(std::istream& in, QPKKT_ProbStats& p) {
    in >> p.Qdim >> p.Vdim >> p.Arowdim >> p.Aeqdim;
    in >> p.Bdim >> p.Cdim;
    int sz;
    in >> sz;
    p.inittime.resize(unsigned(sz));
    for (unsigned i = 0; i < p.inittime.size(); i++)
      in >> p.inittime[i];
    in >> sz;
    p.kktdata.resize(unsigned(sz));
    for (unsigned i = 0; i < p.kktdata.size(); i++)
      in >> p.kktdata[i];
    return in;
  }

}

