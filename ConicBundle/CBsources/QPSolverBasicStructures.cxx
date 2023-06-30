/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/QPSolverBasicStructures.cxx
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
#include "QPSolverBasicStructures.hxx"

using namespace CH_Matrix_Classes;

namespace ConicBundle {



  QPSolverBasicInterface::~QPSolverBasicInterface() {
  }


  QPSolverBasicStructures::QPSolverBasicStructures(QPSolverParameters* params, CBout* cb) :CBout(cb) {
    if (params == 0) {
      paramsp = new QPSolverParameters;
    } else {
      paramsp = params;
    }
    QPIclear();
  }

  void QPSolverBasicStructures::QPIclear() {
    mu = -1.;
    last_mu = -1.;
    next_mu = -1.;
    last_theta = -1.;
    next_theta = -1.;
    last_alpha = -1.;
    sigma = .1;
    primalval = min_Real;
    dualval = max_Real;
    n2primalviol = max_Real;
    n2dualviol = max_Real;

    large_predictor_cnt = 0;
    iter = 0;

    x.init(0, 1, 0.);
    y.init(0, 1, 0.);
    s.init(0, 1, 0.);
    zlb.init(0, 1, 0.);
    zub.init(0, 1, 0.);
    rhszlb.init(0, 1, 0.);
    rhszub.init(0, 1, 0.);

    old_x.init(0, 1, 0.);
    old_y.init(0, 1, 0.);
    old_s.init(0, 1, 0.);
    old_zlb.init(0, 1, 0.);
    old_zub.init(0, 1, 0.);
    old_rhszlb.init(0, 1, 0.);
    old_rhszub.init(0, 1, 0.);

    primalviol.init(0, 1, 0.);
    dualviol.init(0, 1, 0.);
    yviol.init(0, 1, 0.);
    shortzlb.init(0, 1, 0.);
    shortzub.init(0, 1, 0.);
    slacklb.init(0, 1, 0.);
    slackub.init(0, 1, 0.);

    shortrhszlb.init(0, 1, 0.);
    shortrhszub.init(0, 1, 0.);
    rhsslacklb.init(0, 1, 0.);
    rhsslackub.init(0, 1, 0.);

    central_path.clear();
  }


  // *************************************************************************
  //                            QPget_x
  // *************************************************************************

  int QPSolverBasicStructures::QPget_x(Matrix& xout,
    Indexmatrix& x_activity) const {
    xout = x;
    x_activity.init(xout.rowdim(), 1, Integer(1));

    if (old_x.rowdim() != x.rowdim())
      return 0;

    /*
    if (central_path.size()<=1){
      return 0;
    }

    unsigned oldind=unsigned(central_path.size())-2;
    Real old_mu=central_path[oldind].mu;
    const Matrix& old_x=central_path[oldind].x;
    const Matrix& old_zlb=central_path[oldind].zlb;
    const Matrix& old_zub=central_path[oldind].zub;
    */
    const Matrix& lb = QPget_lb();
    const Matrix& ub = QPget_ub();
    const Indexmatrix& lbind = QPget_lbind();
    const Indexmatrix& ubind = QPget_ubind();

    Real tapia_factor = mu / old_mu;
    //std::cout<<" lb[";  //TEST

    for (Integer i = 0; i < lbind.dim(); i++) {
      Integer ind = lbind(i);
      Real slack = x(ind) - lb(ind);
      Real dual = zlb(ind);

      if (slack > 1e-2 * max(1., std::fabs(lb(ind))))
        continue;

      if (dual < 1e-2 * slack)
        continue;

      if (slack < 1e-12 * max(1., std::fabs(lb(ind)))) {
        xout(ind) = lb(ind);
        x_activity(ind) = 0;
        //std::cout<<" "<<ind;  //TEST
        continue;
      }

      if ((mu < 0.01) && (slack < 0.01 * std::sqrt(mu)) && (dual > std::sqrt(mu))) {
        xout(ind) = lb(ind);
        x_activity(ind) = 0;
        //std::cout<<" "<<ind;  //TEST
        continue;
      }

      //check whether reduction is faster than for the dual
      if (tapia_factor < .8) {
        Real primal_tapia = slack / (old_x(ind) - lb(ind));
        Real dual_tapia = zlb(ind) / old_zlb(ind);
        if (
          (primal_tapia < 0.8) &&
          (primal_tapia < .5 * dual_tapia)
          ) {
          //goes to zero faster than the dual value
          xout(ind) = lb(ind);
          x_activity(ind) = 0;
          //std::cout<<" "<<ind;  //TEST
        }
      }
    }
    //std::cout<<"]"<<std::endl;  //TEST

    for (Integer i = 0; i < ubind.dim(); i++) {
      Integer ind = ubind(i);
      Real slack = ub(ind) - x(ind);
      Real dual = zub(ind);

      if (slack > 1e-2 * max(1., std::fabs(ub(ind))))
        continue;

      if (dual < 1e-2 * slack)
        continue;

      if (slack < 1e-12 * max(1., std::fabs(ub(ind)))) {
        xout(ind) = ub(ind);
        x_activity(ind) = 0;
        continue;
      }

      if ((mu < 0.01) && (slack < 0.01 * std::sqrt(mu)) && (dual > std::sqrt(mu))) {
        xout(ind) = ub(ind);
        x_activity(ind) = 0;
        continue;
      }

      //check whether reduction is faster than for the dual
      if (tapia_factor < .8) {
        Real primal_tapia = slack / (ub(ind) - old_x(ind));
        Real dual_tapia = zub(ind) / old_zub(ind);
        if (
          (primal_tapia < 0.8) &&
          (primal_tapia < .5 * dual_tapia)
          ) {
          //goes to zero faster than the dual value
          xout(ind) = ub(ind);
          x_activity(ind) = 0;
        }
      }
    }

    return 0;
  }

  // *************************************************************************
  //                            QPget_s
  // *************************************************************************

  int QPSolverBasicStructures::QPget_s(Matrix& sout,
    Indexmatrix& s_activity) const {
    sout = s;
    s_activity.init(sout.rowdim(), 1, Integer(1));

    const Matrix& lb = QPget_rhslb();
    const Matrix& ub = QPget_rhsub();
    const Indexmatrix& lbind = QPget_rhslbind();
    const Indexmatrix& ubind = QPget_rhsubind();

    //find the equations
    Integer j = 0;
    for (Integer i = 0; (i < lbind.rowdim()); i++) {
      Integer ind = lbind(i);
      while ((j < ubind.rowdim()) && (ubind(j) < ind))
        j++;
      if (j >= ubind.rowdim())
        break;
      if ((ind == ubind(j)) && (ub(ind) - lb(ind) <= 1e-12 * max(1., std::fabs(ub(ind))))) {
        sout(ind) = lb(ind);
        s_activity(ind) = 0;
      }
    }

    if (old_s.rowdim() != s.rowdim())
      return 0;

    /*
    if (central_path.size()<=1){
      return 0;
    }

    unsigned oldind=unsigned(central_path.size())-2;
    Real old_mu=central_path[oldind].mu;
    const Matrix& old_s=central_path[oldind].s;
    const Matrix& old_rhszlb=central_path[oldind].rhszlb;
    const Matrix& old_rhszub=central_path[oldind].rhszub;
    */

    Real tapia_factor = mu / old_mu;

    for (Integer i = 0; i < lbind.dim(); i++) {
      Integer ind = lbind(i);
      Real slack = -s(ind) - lb(ind);
      Real dual = zlb(ind);

      if (slack > 1e-2 * max(1., std::fabs(lb(ind))))
        continue;

      if (dual < 1e-2 * slack)
        continue;

      if (slack < 1e-12 * max(1., std::fabs(lb(ind)))) {
        sout(ind) = -lb(ind);
        s_activity(ind) = 0;
        continue;
      }

      if ((mu < 0.01) && (slack < 0.01 * std::sqrt(mu)) && (dual > std::sqrt(mu))) {
        sout(ind) = -lb(ind);
        s_activity(ind) = 0;
        continue;
      }

      //check whether reduction is faster than for the dual
      if (tapia_factor < .8) {
        Real primal_tapia = slack / (-old_s(ind) - lb(ind));
        Real dual_tapia = zlb(ind) / old_rhszlb(ind);
        if (
          (primal_tapia < 0.8) &&
          (primal_tapia < .5 * dual_tapia)
          ) {
          //goes to zero faster than the dual value
          sout(ind) = -lb(ind);
          s_activity(ind) = 0;
        }
      }
    }

    for (Integer i = 0; i < ubind.dim(); i++) {
      Integer ind = ubind(i);
      Real slack = ub(ind) + s(ind);
      Real dual = zub(ind);

      if (slack > 1e-2 * max(1., std::fabs(ub(ind))))
        continue;

      if (dual < 1e-2 * slack)
        continue;

      if (slack < 1e-12 * max(1., std::fabs(ub(ind)))) {
        sout(ind) = -ub(ind);
        s_activity(ind) = 0;
        continue;
      }

      if ((mu < 0.01) && (slack < 0.01 * std::sqrt(mu)) && (dual > std::sqrt(mu))) {
        sout(ind) = -ub(ind);
        s_activity(ind) = 0;
        continue;
      }

      //check whether reduction is faster than for the dual
      if (tapia_factor < .8) {
        Real primal_tapia = slack / (ub(ind) + old_s(ind));
        Real dual_tapia = zub(ind) / old_rhszub(ind);
        if (
          (primal_tapia < 0.8) &&
          (primal_tapia < .5 * dual_tapia)
          ) {
          //goes to zero faster than the dual value
          sout(ind) = -ub(ind);
          s_activity(ind) = 0;
        }
      }
    }

    return 0;
  }

  // *************************************************************************
  //                            QPcompute_values
  // *************************************************************************

  int QPSolverBasicStructures::QPcompute_values(bool init) {
    int err = 0;
    // compute primal and dual objective value, initialize "short" versions

    primalviol.init(s);
    dualval = ip(y, primalviol) + QPget_gamma();
    QPadd_Ax(x, primalviol);
    if (paramsp->QPget_use_socqp()) {
      dualviol.init(QPget_c());
      primalval = 0.;
    } else {
      dualviol.init(x.rowdim(), 1, 0.);
      QPadd_Qx(x, dualviol);
      primalval = ip(dualviol, x) / 2.;
      dualval -= primalval;
      dualviol += QPget_c();
    }
    primalval += ip(QPget_c(), x) + QPget_gamma();
    QPadd_Aty(y, dualviol);

    shortzlb.newsize(QPget_lbind().dim(), 1); chk_set_init(shortzlb, 1);
    slacklb.newsize(QPget_lbind().dim(), 1); chk_set_init(slacklb, 1);
    for (Integer i = 0; i < shortzlb.dim(); i++) {
      Integer ind = QPget_lbind()(i);
      Real d = zlb(ind);
      shortzlb(i) = d;
      dualviol(ind) -= d;
      Real lbi = QPget_lb()(ind);
      dualval += lbi * d;
      slacklb(i) = x(ind) - lbi;
    }
    shortzub.newsize(QPget_ubind().dim(), Integer(1)); chk_set_init(shortzub, 1);
    slackub.newsize(QPget_ubind().dim(), Integer(1)); chk_set_init(slackub, 1);
    for (Integer i = 0; i < shortzub.dim(); i++) {
      Integer ind = QPget_ubind()(i);
      Real d = zub(ind);
      shortzub(i) = d;
      dualviol(ind) += d;
      Real ubi = QPget_ub()(ind);
      dualval -= ubi * d;
      slackub(i) = ubi - x(ind);
    }
    yviol.init(QPget_ydim(), 1, 0.);
    shortrhszlb.newsize(QPget_rhsubind().dim(), 1); chk_set_init(shortrhszlb, 1);
    rhsslacklb.newsize(QPget_rhsubind().dim(), 1); chk_set_init(rhsslacklb, 1);
    shortrhszub.newsize(QPget_rhslbind().dim(), 1); chk_set_init(shortrhszub, 1);
    rhsslackub.newsize(QPget_rhslbind().dim(), 1); chk_set_init(rhsslackub, 1);
    Integer lbi = 0;
    Integer lbind;
    if (lbi < QPget_rhsubind().dim()) {
      lbind = QPget_rhsubind()(lbi);
    } else {
      lbind = QPget_ydim();
    }
    Integer ubi = 0;
    Integer ubind;
    if (ubi < QPget_rhslbind().dim()) {
      ubind = QPget_rhslbind()(ubi);
    } else {
      ubind = QPget_ydim();
    }
    while ((lbind < QPget_ydim()) || (ubind < QPget_ydim())) {
      if (lbind == ubind) {
        Real zl = rhszlb(lbind);
        shortrhszlb(lbi) = zl;
        Real rl = -QPget_rhsub()(lbind);
        dualval += (rl - s(lbind)) * zl;
        rhsslacklb(lbi) = s(lbind) - rl;

        Real zu = rhszub(ubind);
        shortrhszub(ubi) = zu;
        Real ru = -QPget_rhslb()(ubind);
        dualval += (s(ubind) - ru) * zu;
        rhsslackub(ubi) = ru - s(ubind);

        yviol(lbind) = y(lbind) - zl + zu;

        lbi++;
        if (lbi < QPget_rhsubind().dim()) {
          lbind = QPget_rhsubind()(lbi);
        } else {
          lbind = QPget_ydim();
        }
        ubi++;
        if (ubi < QPget_rhslbind().dim()) {
          ubind = QPget_rhslbind()(ubi);
        } else {
          ubind = QPget_ydim();
        }
      } else if (lbind < ubind) {
        Real zl = rhszlb(lbind);
        shortrhszlb(lbi) = zl;
        Real rl = -QPget_rhsub()(lbind);
        dualval += (rl - s(lbind)) * zl;
        rhsslacklb(lbi) = s(lbind) - rl;

        yviol(lbind) = y(lbind) - zl;

        lbi++;
        if (lbi < QPget_rhsubind().dim()) {
          lbind = QPget_rhsubind()(lbi);
        } else {
          lbind = QPget_ydim();
        }
      } else {
        Real zu = rhszub(ubind);
        shortrhszub(ubi) = zu;
        Real ru = -QPget_rhslb()(ubind);
        dualval += (s(ubind) - ru) * zu;
        rhsslackub(ubi) = ru - s(ubind);

        yviol(ubind) = y(ubind) + zu;

        ubi++;
        if (ubi < QPget_rhslbind().dim()) {
          ubind = QPget_rhslbind()(ubi);
        } else {
          ubind = QPget_ydim();
        }
      }
    }

    if (model_block) {
      primalval += model_block->globalx_cost(x) + model_block->constraints_cost();
      model_block->add_Bt_modelx(dualval, dualviol);
      n2modeldualviol = std::sqrt(model_block->dualviol_2normsqr());
      n2modelprimalviol = std::sqrt(model_block->primalviol_2normsqr());
    } else {
      n2modeldualviol = 0.;
      n2modelprimalviol = 0.;
    }

    if (paramsp->QPget_use_socqp()) {
      if ((init) && (socqp.reset_starting_point(x, dualviol, mu))) {
        if (cb_out()) {
          get_out() << "**** ERROR in QPSolverBasicStructures::QPpredcorr_step(...): socqp.reset_starting_point(...) failed" << std::endl;
        }
        err++;
      } else if (socqp.add_prox_contrib(primalval, dualval, dualviol)) {
        if (cb_out()) {
          get_out() << "**** ERROR in QPSolverBasicStructures::QPpredcorr_step(...): socqp.add_prox_contrib(...) failed" << std::endl;
        }
        err++;
      }
    }


    n2primalviol = norm2(primalviol);
    n2dualviol = norm2(dualviol);
    n2yviol = norm2(yviol);
    return err;
  }

  // *************************************************************************
  //                            QPselect_localstep_mu
  // *************************************************************************

  int QPSolverBasicStructures::QPselect_localstep_mu(Real& mu,
    Real& stepsize,
    bool centering) {

    Integer mudim = 0;
    Real globminval = max_Real;
    Real globmaxval = min_Real;
    Real minv, maxv;
    Real tr_xz = 0.;
    Real tr_xdzpdxz = 0.;
    Real tr_dxdz = 0.;

    Real nbh_ubnd = paramsp->QPget_nbh_ub();
    Real nbh_lbnd = paramsp->QPget_nbh_lb();
    stepsize = 1.;

    if (slacklb.rowdim() > 0) {
      mudim += slacklb.rowdim();
      tr_xz += ip_min_max(slacklb, shortzlb, minv, maxv);
      tr_xdzpdxz += ip(slacklb, dshortzlb) + ip(dslacklb, shortzlb);
      tr_dxdz += ip(dslacklb, dshortzlb);

      if (minv < globminval)
        globminval = minv;
      if (maxv > globmaxval)
        globmaxval = maxv;

      if (cb_out(3))
        get_out() << " slb[" << minv << "," << maxv << "]";
    }

    if (slackub.rowdim() > 0) {
      mudim += slackub.rowdim();
      tr_xz += ip_min_max(slackub, shortzub, minv, maxv);
      tr_xdzpdxz += ip(slackub, dshortzub) + ip(dslackub, shortzub);
      tr_dxdz += ip(dslackub, dshortzub);

      if (minv < globminval)
        globminval = minv;
      if (maxv > globmaxval)
        globmaxval = maxv;

      if (cb_out(3))
        get_out() << " sub[" << minv << "," << maxv << "]";
    }

    if (rhsslacklb.rowdim() > 0) {
      mudim += rhsslacklb.rowdim();
      tr_xz += ip_min_max(rhsslacklb, shortrhszlb, minv, maxv);
      tr_xdzpdxz += ip(rhsslacklb, dshortrhszlb) + ip(drhsslacklb, shortrhszlb);
      tr_dxdz += ip(drhsslacklb, dshortrhszlb);

      if (minv < globminval)
        globminval = minv;
      if (maxv > globmaxval)
        globmaxval = maxv;

      if (cb_out(3))
        get_out() << " rslb[" << minv << "," << maxv << "]";
    }

    if (rhsslackub.rowdim() > 0) {
      mudim += rhsslackub.rowdim();
      tr_xz += ip_min_max(rhsslackub, shortrhszub, minv, maxv);
      tr_xdzpdxz += ip(rhsslackub, dshortrhszub) + ip(drhsslackub, shortrhszub);
      tr_dxdz += ip(drhsslackub, dshortrhszub);

      if (minv < globminval)
        globminval = minv;
      if (maxv > globmaxval)
        globmaxval = maxv;

      if (cb_out(3))
        get_out() << " rsub[" << minv << "," << maxv << "]";
    }


    if (model_block) {
      if (model_block->get_mu_info(mudim, tr_xz, tr_xdzpdxz, tr_dxdz, globminval, globmaxval)) {
        if (cb_out()) {
          get_out() << "**** ERROR in QPSolverBasicStructures::QPselect_localstep_mu(...): model_block->get_mu_info(...) failed" << std::endl;
        }
      }
    }

    if (paramsp->QPget_use_socqp()) {
      if (socqp.get_mu_info(mudim, tr_xz, tr_xdzpdxz, tr_dxdz, globminval, globmaxval)) {
        if (cb_out()) {
          get_out() << "**** ERROR in QPSolverBasicStructures::QPselect_localstep_mu(...): socqp.get_mu_info(...) failed" << std::endl;
        }
      }
    }

    if (cb_out(3))
      get_out() << std::endl;


    //first do centering, only by linesearch 
    if (centering) {
      QPlinesearch(stepsize);
      Real smu = (tr_xz + stepsize * (tr_xdzpdxz + stepsize * tr_dxdz)) / mudim;
      mu = smu;
      //mu=.5*mu+.5*smu;
      //mu=.9*mu+.1*smu;
      sigma = mu / smu;
      if (cb_out(3)) {
        get_out() << " centering smu=" << smu << " sigma=" << sigma << " nextmu=" << mu << std::endl;
      }
      return 0;
    }

    const Real mu_xz = tr_xz / mudim;
    const Real mu_xdzpdxz = tr_xdzpdxz / mudim;
    const Real mu_dxdz = tr_dxdz / mudim;
    const Real mu_at_one = mu_xz + mu_xdzpdxz + mu_dxdz;

    Real nrmsqr_xz = 0.;
    Real nrmsqr_xdzpdxz = 0.;
    Real nrmsqr_dxdz = 0.;
    Real ip_xz_xdzpdxz = 0.;
    Real ip_xz_dxdz = 0.;
    Real ip_dxdz_xdzpdxz = 0.;
    Real max_nbh = 0.;

    if (slacklb.rowdim() > 0) {
      for (Integer i = 0; i < slacklb.rowdim(); i++) {
        NNC_nbh_stepsize(slacklb(i), shortzlb(i), dslacklb(i), dshortzlb(i),
          mu_xz, mu_xdzpdxz, mu_dxdz, mu_at_one,
          nbh_ubnd, stepsize, max_nbh,
          nrmsqr_xz, nrmsqr_xdzpdxz, nrmsqr_dxdz,
          ip_xz_xdzpdxz, ip_xz_dxdz, ip_dxdz_xdzpdxz);
      }
    }
    if (slackub.rowdim() > 0) {
      for (Integer i = 0; i < slackub.rowdim(); i++) {
        NNC_nbh_stepsize(slackub(i), shortzub(i), dslackub(i), dshortzub(i),
          mu_xz, mu_xdzpdxz, mu_dxdz, mu_at_one,
          nbh_ubnd, stepsize, max_nbh,
          nrmsqr_xz, nrmsqr_xdzpdxz, nrmsqr_dxdz,
          ip_xz_xdzpdxz, ip_xz_dxdz, ip_dxdz_xdzpdxz);
      }
    }
    if (rhsslacklb.rowdim() > 0) {
      for (Integer i = 0; i < rhsslacklb.rowdim(); i++) {
        NNC_nbh_stepsize(rhsslacklb(i), shortrhszlb(i), drhsslacklb(i), dshortrhszlb(i),
          mu_xz, mu_xdzpdxz, mu_dxdz, mu_at_one,
          nbh_ubnd, stepsize, max_nbh,
          nrmsqr_xz, nrmsqr_xdzpdxz, nrmsqr_dxdz,
          ip_xz_xdzpdxz, ip_xz_dxdz, ip_dxdz_xdzpdxz);
      }
    }
    if (rhsslackub.rowdim() > 0) {
      for (Integer i = 0; i < rhsslackub.rowdim(); i++) {
        NNC_nbh_stepsize(rhsslackub(i), shortrhszub(i), drhsslackub(i), dshortrhszub(i),
          mu_xz, mu_xdzpdxz, mu_dxdz, mu_at_one,
          nbh_ubnd, stepsize, max_nbh,
          nrmsqr_xz, nrmsqr_xdzpdxz, nrmsqr_dxdz,
          ip_xz_xdzpdxz, ip_xz_dxdz, ip_dxdz_xdzpdxz);
      }
    }
    if (model_block) {
      if (model_block->get_nbh_info(mudim, tr_xz, tr_xdzpdxz, tr_dxdz, nbh_ubnd,
        stepsize, max_nbh,
        nrmsqr_xz, nrmsqr_xdzpdxz, nrmsqr_dxdz,
        ip_xz_xdzpdxz, ip_xz_dxdz, ip_dxdz_xdzpdxz)) {
        if (cb_out()) {
          get_out() << "**** ERROR in QPSolverBasicStructures::QPselect_localstep_mu(...): model_block->get_nbh_info(...) failed" << std::endl;
        }
      }
    }

    if (paramsp->QPget_use_socqp()) {
      if (socqp.get_nbh_info(mudim, tr_xz, tr_xdzpdxz, tr_dxdz, nbh_ubnd,
        stepsize, max_nbh,
        nrmsqr_xz, nrmsqr_xdzpdxz, nrmsqr_dxdz,
        ip_xz_xdzpdxz, ip_xz_dxdz, ip_dxdz_xdzpdxz)) {
        if (cb_out()) {
          get_out() << "**** ERROR in QPSolverBasicStructures::QPselect_localstep_mu(...): socqp.get_nbh_info(...) failed" << std::endl;
        }
      }
    }


    Real omu = tr_xz / mudim;
    assert(omu > 0);
    Real gmu = max(eps_Real * tr_xz, (tr_xz + tr_xdzpdxz)) / mudim;
    Real nmu = max(eps_Real * tr_xz, (tr_xz + tr_xdzpdxz + tr_dxdz)) / mudim;

    Real otheta = std::sqrt(nrmsqr_xz) / omu;
    last_theta = otheta;


    Real theta = max(nbh_ubnd, otheta);
    if (cb_out(3)) {
      get_out().precision(6);
      get_out() << "\n lnbh: alpha=" << stepsize << " tr_xz=" << tr_xz << " tr_dxdz=" << tr_dxdz << " omu=" << omu << " gmu=" << gmu << " nmu=" << nmu << " otheta=" << otheta << " theta=" << theta << std::endl;
    }

    //form the polynomial of the squared distance to the central path minus neighborhood threshold
    Real q0 = nrmsqr_xz - sqr(theta / mudim * tr_xz);
    Real q1 = 2. * ip_xz_xdzpdxz - 2 * sqr(theta / mudim) * tr_xz * tr_xdzpdxz;
    Real q2 = nrmsqr_xdzpdxz + 2. * ip_xz_dxdz - sqr(theta / mudim) * (sqr(tr_xdzpdxz) + 2 * tr_xz * tr_dxdz);
    Real q3 = 2. * ip_dxdz_xdzpdxz - 2 * sqr(theta / mudim) * tr_xdzpdxz * tr_dxdz;
    Real q4 = nrmsqr_dxdz - sqr(theta * tr_dxdz / mudim);

    //Real maxnrmq=max(std::fabs(q0),max(std::fabs(q1),max(std::fabs(q2),max(std::fabs(q3),std::fabs(q4)))));

    if (cb_out(3)) {
      get_out().precision(12);
      get_out() << " q0=" << q0 << " q1=" << q1 << " q2=" << q2 << " q3=" << q3 << " q4=" << q4 << std::endl;
    }


    // if (centering){
    //   Real smu=(tr_xz+stepsize*(tr_xdzpdxz+stepsize*tr_dxdz))/mudim;
    //   //Real sqrnrm=nrmsqr_xz+stepsize*(2*ip_xz_xdzpdxz+stepsize*(nrmsqr_xdzpdxz+2*ip_xz_dxdz+stepsize*(2*ip_dxdz_xdzpdxz+stepsize*nrmsqr_dxdz)));
    //   //next_theta=std::sqrt(sqrnrm)/smu;
    //   next_theta=max_nbh;
    //   if (q1>0.){
    //      mu=10*mu;
    //   }
    //   else {
    //     mu=.9*min(gmu,smu)+.1*nmu;
    //   }
    //   sigma=mu/smu;
    //   if (cb_out(3)){
    //     get_out()<<" centering smu="<<smu<<" next_theta="<<next_theta<<" sigma="<<sigma<<" nextmu="<<mu<<std::endl;
    //   }
    //   return 0;
    // }  

    //ploynomial coefficients of the derivative dq
    Real dq0 = q1;
    Real dq1 = 2 * q2;
    Real dq2 = 3 * q3;
    Real dq3 = 4 * q4;

    Real qsv = (((q4 * stepsize + q3) * stepsize + q2) * stepsize + q1) * stepsize + q0;
    Real dqsv = ((dq3 * stepsize + dq2) * stepsize + dq1) * stepsize + dq0;
    Real smu = (tr_xz + stepsize * (tr_xdzpdxz + stepsize * tr_dxdz)) / mudim;
    Real sqrnrm = nrmsqr_xz + stepsize * (2 * ip_xz_xdzpdxz + stepsize * (nrmsqr_xdzpdxz + 2 * ip_xz_dxdz + stepsize * (2 * ip_dxdz_xdzpdxz + stepsize * nrmsqr_dxdz)));
    if (((smu < 0.) || (sqrnrm < 0.)) && (cb_out())) {
      get_out() << "**** WARNING in QPSolverBasicStructures::QPselect_localstep_mu(...): squared distance to central path at step computes to " << sqrnrm << " with mu =" << smu << " (both need to be positive)" << std::endl;
    }
    if ((std::fabs(qsv - (sqrnrm - sqr(theta * smu))) > 1e-6 * omu) && (cb_out())) {
      get_out().precision(12);
      get_out() << "\n trxz=" << tr_xz << " trdxpz=" << tr_xdzpdxz << " trdxdz=" << tr_dxdz << std::endl;
      get_out() << " n2xz=" << nrmsqr_xz << " n2dxpz=" << nrmsqr_xdzpdxz << " n2dxdz=" << nrmsqr_dxdz << std::endl;
      get_out() << " ixzdxpz=" << ip_xz_xdzpdxz << " ixzdxdz=" << ip_xz_dxdz << " idxpzdxdz=" << ip_dxdz_xdzpdxz << std::endl;
      get_out() << " s-q0=" << q0 << " q1=" << q1 << " q2=" << q2 << " q3=" << q3 << " q4=" << q4 << std::endl;
      get_out() << " st=" << stepsize << " sqrnrm=" << sqrnrm << " smu=" << smu << " theta=" << theta << " qsv=" << qsv << " diff=" << qsv - (sqrnrm - sqr(theta * smu)) << std::endl;
    }
    //if ((q0<0.)&&(qsv>1e-6*max(1.,omu))){
    //  assert((q0>0)||(qsv<1e-6*max(1.,omu)));
    //}
    smu = max(smu, eps_Real * tr_xz);
    sqrnrm = max(sqrnrm, eps_Real * tr_xz);
    Real stheta = std::sqrt(sqrnrm) / smu;

    next_theta = max_nbh;
    if (cb_out(3)) {
      get_out() << " qsv=" << qsv << " dqsv=" << dqsv;
      get_out() << " smu=" << smu << " stheta=" << stheta << " max_nbh=" << max_nbh << std::endl;
    }

    if ((max_nbh > nbh_ubnd + 1e-6) && (stepsize < 0.0001)) {
      sigma = 1.1;
    } else {
      if (gmu < 1e-6 * omu) {
        sigma = 1. - 0.35 / std::sqrt(mudim);
      } else {
        sigma = min(gmu / omu, 1. - 0.35 / std::sqrt(mudim));
      }
      //sigma=std::pow(sigma,min(4.,stepsize*nbh_lbnd/stheta));
      sigma = std::pow(sigma, min(4., stepsize * nbh_lbnd / max_nbh));
    }
    if (max_nbh < nbh_ubnd)
      sigma = min(sigma, 1. - 0.35 / std::sqrt(mudim));
    mu = sigma * smu;

    if (cb_out(3))
      get_out() << " sigma=" << sigma << " nextmu=" << mu << std::endl;

    return 0;
  }

  // *************************************************************************
  //                            QPselect_step_mu
  // *************************************************************************

  int QPSolverBasicStructures::QPselect_step_mu(Real& mu,
    Real& stepsize,
    bool centering) {

    Integer mudim = 0;
    Real globminval = max_Real;
    Real globmaxval = min_Real;
    Real minv, maxv;
    Real tr_xz = 0.;
    Real tr_xdzpdxz = 0.;
    Real tr_dxdz = 0.;

    Real nbh_ubnd = paramsp->QPget_nbh_ub();
    Real nbh_lbnd = paramsp->QPget_nbh_lb();
    stepsize = 1.;

    if (slacklb.rowdim() > 0) {
      mudim += slacklb.rowdim();
      tr_xz += ip_min_max(slacklb, shortzlb, minv, maxv);
      tr_xdzpdxz += ip(slacklb, dshortzlb) + ip(dslacklb, shortzlb);
      tr_dxdz += ip(dslacklb, dshortzlb);

      if (minv < globminval)
        globminval = minv;
      if (maxv > globmaxval)
        globmaxval = maxv;

      if (cb_out(3))
        get_out() << " slb[" << minv << "," << maxv << "]";
    }

    if (slackub.rowdim() > 0) {
      mudim += slackub.rowdim();
      tr_xz += ip_min_max(slackub, shortzub, minv, maxv);
      tr_xdzpdxz += ip(slackub, dshortzub) + ip(dslackub, shortzub);
      tr_dxdz += ip(dslackub, dshortzub);

      if (minv < globminval)
        globminval = minv;
      if (maxv > globmaxval)
        globmaxval = maxv;

      if (cb_out(3))
        get_out() << " sub[" << minv << "," << maxv << "]";
    }

    if (rhsslacklb.rowdim() > 0) {
      mudim += rhsslacklb.rowdim();
      tr_xz += ip_min_max(rhsslacklb, shortrhszlb, minv, maxv);
      tr_xdzpdxz += ip(rhsslacklb, dshortrhszlb) + ip(drhsslacklb, shortrhszlb);
      tr_dxdz += ip(drhsslacklb, dshortrhszlb);

      if (minv < globminval)
        globminval = minv;
      if (maxv > globmaxval)
        globmaxval = maxv;

      if (cb_out(3))
        get_out() << " rslb[" << minv << "," << maxv << "]";
    }

    if (rhsslackub.rowdim() > 0) {
      mudim += rhsslackub.rowdim();
      tr_xz += ip_min_max(rhsslackub, shortrhszub, minv, maxv);
      tr_xdzpdxz += ip(rhsslackub, dshortrhszub) + ip(drhsslackub, shortrhszub);
      tr_dxdz += ip(drhsslackub, dshortrhszub);

      if (minv < globminval)
        globminval = minv;
      if (maxv > globmaxval)
        globmaxval = maxv;

      if (cb_out(3))
        get_out() << " rsub[" << minv << "," << maxv << "]";
    }


    if (model_block) {
      if (model_block->get_mu_info(mudim, tr_xz, tr_xdzpdxz, tr_dxdz, globminval, globmaxval)) {
        if (cb_out()) {
          get_out() << "**** ERROR in QPSolverBasicStructures::QPselect_step_mu(...): model_block->get_mu_info(...) failed" << std::endl;
        }
      }
    }

    if (paramsp->QPget_use_socqp()) {
      if (socqp.get_mu_info(mudim, tr_xz, tr_xdzpdxz, tr_dxdz, globminval, globmaxval)) {
        if (cb_out()) {
          get_out() << "**** ERROR in QPSolverBasicStructures::QPselect_step_mu(...): socqp.get_mu_info(...) failed" << std::endl;
        }
      }
    }

    const Real mu_xz = tr_xz / mudim;
    const Real mu_xdzpdxz = tr_xdzpdxz / mudim;
    const Real mu_dxdz = tr_dxdz / mudim;
    const Real mu_at_one = mu_xz + mu_xdzpdxz + mu_dxdz;

    Real nrmsqr_xz = 0.;
    Real nrmsqr_xdzpdxz = 0.;
    Real nrmsqr_dxdz = 0.;
    Real ip_xz_xdzpdxz = 0.;
    Real ip_xz_dxdz = 0.;
    Real ip_dxdz_xdzpdxz = 0.;
    Real max_nbh = 0.;

    if (slacklb.rowdim() > 0) {
      for (Integer i = 0; i < slacklb.rowdim(); i++) {
        NNC_nbh_stepsize(slacklb(i), shortzlb(i), dslacklb(i), dshortzlb(i),
          mu_xz, mu_xdzpdxz, mu_dxdz, mu_at_one,
          nbh_ubnd, stepsize, max_nbh,
          nrmsqr_xz, nrmsqr_xdzpdxz, nrmsqr_dxdz,
          ip_xz_xdzpdxz, ip_xz_dxdz, ip_dxdz_xdzpdxz);
      }
    }
    if (slackub.rowdim() > 0) {
      for (Integer i = 0; i < slackub.rowdim(); i++) {
        NNC_nbh_stepsize(slackub(i), shortzub(i), dslackub(i), dshortzub(i),
          mu_xz, mu_xdzpdxz, mu_dxdz, mu_at_one,
          nbh_ubnd, stepsize, max_nbh,
          nrmsqr_xz, nrmsqr_xdzpdxz, nrmsqr_dxdz,
          ip_xz_xdzpdxz, ip_xz_dxdz, ip_dxdz_xdzpdxz);
      }
    }
    if (rhsslacklb.rowdim() > 0) {
      for (Integer i = 0; i < rhsslacklb.rowdim(); i++) {
        NNC_nbh_stepsize(rhsslacklb(i), shortrhszlb(i), drhsslacklb(i), dshortrhszlb(i),
          mu_xz, mu_xdzpdxz, mu_dxdz, mu_at_one,
          nbh_ubnd, stepsize, max_nbh,
          nrmsqr_xz, nrmsqr_xdzpdxz, nrmsqr_dxdz,
          ip_xz_xdzpdxz, ip_xz_dxdz, ip_dxdz_xdzpdxz);
      }
    }
    if (rhsslackub.rowdim() > 0) {
      for (Integer i = 0; i < rhsslackub.rowdim(); i++) {
        NNC_nbh_stepsize(rhsslackub(i), shortrhszub(i), drhsslackub(i), dshortrhszub(i),
          mu_xz, mu_xdzpdxz, mu_dxdz, mu_at_one,
          nbh_ubnd, stepsize, max_nbh,
          nrmsqr_xz, nrmsqr_xdzpdxz, nrmsqr_dxdz,
          ip_xz_xdzpdxz, ip_xz_dxdz, ip_dxdz_xdzpdxz);
      }
    }
    if (model_block) {
      if (model_block->get_nbh_info(mudim, tr_xz, tr_xdzpdxz, tr_dxdz, nbh_ubnd,
        stepsize, max_nbh,
        nrmsqr_xz, nrmsqr_xdzpdxz, nrmsqr_dxdz,
        ip_xz_xdzpdxz, ip_xz_dxdz, ip_dxdz_xdzpdxz)) {
        if (cb_out()) {
          get_out() << "**** ERROR in QPSolverBasicStructures::QPselect_step_mu(...): model_block->get_nbh_info(...) failed" << std::endl;
        }
      }
    }

    if (paramsp->QPget_use_socqp()) {
      if (socqp.get_nbh_info(mudim, tr_xz, tr_xdzpdxz, tr_dxdz, nbh_ubnd,
        stepsize, max_nbh,
        nrmsqr_xz, nrmsqr_xdzpdxz, nrmsqr_dxdz,
        ip_xz_xdzpdxz, ip_xz_dxdz, ip_dxdz_xdzpdxz)) {
        if (cb_out()) {
          get_out() << "**** ERROR in QPSolverBasicStructures::QPselect_step_mu(...): socqp.get_nbh_info(...) failed" << std::endl;
        }
      }
    }


    Real omu = tr_xz / mudim;
    assert(omu > 0);
    Real gmu = max(eps_Real * tr_xz, (tr_xz + tr_xdzpdxz)) / mudim;
    Real nmu = max(eps_Real * tr_xz, (tr_xz + tr_xdzpdxz + tr_dxdz)) / mudim;

    Real otheta = std::sqrt(nrmsqr_xz) / omu;
    last_theta = otheta;


    Real theta = max(nbh_ubnd, otheta);
    if (cb_out(3)) {
      get_out().precision(6);
      get_out() << "\n nbh: tr_xz=" << tr_xz << " tr_dxdz=" << tr_dxdz << " omu=" << omu << " gmu=" << gmu << " nmu=" << nmu << " otheta=" << otheta << " theta=" << theta << std::endl;
    }

    //form the polynomial of the squared distance to the central path minus neighborhood threshold
    Real q0 = nrmsqr_xz - sqr(theta / mudim * tr_xz);
    Real q1 = 2. * ip_xz_xdzpdxz - 2 * sqr(theta / mudim) * tr_xz * tr_xdzpdxz;
    Real q2 = nrmsqr_xdzpdxz + 2. * ip_xz_dxdz - sqr(theta / mudim) * (sqr(tr_xdzpdxz) + 2 * tr_xz * tr_dxdz);
    Real q3 = 2. * ip_dxdz_xdzpdxz - 2 * sqr(theta / mudim) * tr_xdzpdxz * tr_dxdz;
    Real q4 = nrmsqr_dxdz - sqr(theta * tr_dxdz / mudim);

    Real maxnrmq = max(std::fabs(q0), max(std::fabs(q1), max(std::fabs(q2), max(std::fabs(q3), std::fabs(q4)))));

    if (cb_out(3)) {
      get_out().precision(12);
      get_out() << " q0=" << q0 << " q1=" << q1 << " q2=" << q2 << " q3=" << q3 << " q4=" << q4 << std::endl;
    }

    //if centering is requested or otheta is still outside
    if ((centering) || (otheta > 1 - 1e-6)) {
      QPlinesearch(stepsize);
      Real smu = (tr_xz + stepsize * (tr_xdzpdxz + stepsize * tr_dxdz)) / mudim;
      Real sqrnrm = nrmsqr_xz + stepsize * (2 * ip_xz_xdzpdxz + stepsize * (nrmsqr_xdzpdxz + 2 * ip_xz_dxdz + stepsize * (2 * ip_dxdz_xdzpdxz + stepsize * nrmsqr_dxdz)));
      next_theta = std::sqrt(sqrnrm) / smu;
      // if ((next_theta>.99)||(stepsize<0.01))
      //   mu=max(10.*mu,std::sqrt(sqrnrm)/nbh_ubnd);
      // else if (next_theta<nbh_ubnd) {
      //   mu=(1.-0.35/std::sqrt(mudim))*mu;
      // }
      //mu=max(mu,min(10*mu,std::sqrt(sqrnrm)/nbh_ubnd));
      if (q1 > 0.) {
        mu = 10 * mu;
      }
      // if (q1>0.){
      //   mu=max(10*mu,omu);
      // }
      // else {
      //   mu=max(mu,omu);
      // }
      sigma = mu / smu;
      if (cb_out(3)) {
        get_out() << " centering step=" << stepsize << " smu=" << smu << " stheta=" << next_theta << " sigma=" << sigma << " nextmu=" << mu << std::endl;
      }
      return 0;
    }


    //if theta is too large, try to require a relative decrease to this
    bool origq = true;
    if (theta > nbh_ubnd + 1e-6) {
      origq = false;
      Real delta = nbh_ubnd - theta;
      q1 -= omu * omu * 2. * delta * theta;
      //q2-=omu*omu*delta*delta+2*delta*theta*sqr(theta/mudim)*(sqr(tr_xdzpdxz)+2*tr_xz*tr_dxdz);
      //q3-=delta*delta*sqr(theta/mudim)*(sqr(tr_xdzpdxz)+2*tr_xz*tr_dxdz);
      if (cb_out(3)) {
        get_out() << " q0=" << q0 << " q1=" << q1 << " q2=" << q2 << " q3=" << q3 << " q4=" << q4 << std::endl;
      }
    }

    // find the maximum stepval in [0,1] so that q(stepval)<=0.
    // note: q(0)=q0<=0.
    Real ubstep = 1.;
    Real stepval = 1.;
    bool rootsearch = true;


    //ploynomial coefficients of the derivative dq
    Real dq0 = q1;
    Real dq1 = 2 * q2;
    Real dq2 = 3 * q3;
    Real dq3 = 4 * q4;

    //ploynomial coefficients of the second derivative d2q
    Real d2q0 = dq1;
    Real d2q1 = 2 * dq2;
    Real d2q2 = 3 * dq3;

    //discriminant for finding the roots of the second derivative

    Real discr = sqr(d2q1) - 4. * d2q2 * d2q0;

    if (std::fabs(d2q2) < eps_Real * maxnrmq) {
      //only a cubic
      if (std::fabs(dq3) < eps_Real * maxnrmq) {
        //only a quadratic
        if (std::fabs(q2) < eps_Real * maxnrmq) {
          //only linear
          //std::cout<<" linear "<<std::endl;
          rootsearch = false;
          if (std::fabs(q1) < eps_Real * maxnrmq) {
            //constant;
            stepval = 1.;
          } else {
            if (q1 < 0.)
              stepval = 1.;
            else {
              stepval = min(-q0 / q1, 1.);
            }
          }
        } else {
          //full quadratic
          //std::cout<<" quadratic "<<std::endl;
          if (q2 < 0) {  //concave case
            if (q1 <= 0.) {
              ubstep = 1.;
              stepval = 1.;
              rootsearch = false;
            } else {
              Real bx = -q1 / q2 / 2.;  //>0.
              Real bv = (q2 * bx + q1) * bx + q0;
              if (bv < 0) {
                ubstep = 1.;
                stepval = 1;
                rootsearch = false;
              } else {
                ubstep = min(1., bx);
                stepval = ubstep / 2.;
                rootsearch = true;
              }
            }
          } else {  //convex case
            Real bx = -q1 / q2 / 2.;  //minimum
            if (bx >= 1.) {
              ubstep = 1.;
              stepval = 1.;
              rootsearch = false;
            } else {
              ubstep = 1.;
              stepval = 1.;
              rootsearch = true;
            }
          }
        }
      } else { //full cubic,
        Real infl = -dq1 / dq2 / 2.; //inflection point
        Real qinfl = ((q3 * infl + q2) * infl + q1) * infl + q0;
        //find zeros of derivative
        discr = sqr(dq1) - 4. * dq2 * dq0;
        //std::cout<<" cubic discr="<<discr<<std::endl;
        if (discr <= eps_Real * omu) { //no zeros or a single zero of the derviative, function is monotone
          if (q3 < 0) {
            ubstep = 1.;
            stepval = 1.;
            rootsearch = false;
          } else {
            if ((infl >= 1.) || ((infl > 0.) && (qinfl >= 0.))) {
              ubstep = min(1., infl);
              stepval = 0.;
              rootsearch = true;
            } else {
              ubstep = min(1., infl);
              stepval = 1.;
              rootsearch = true;
            }
          }
        } //end no extrema
        else {  //minimum and maximum of full cubic exist
          Real x0 = (-dq1 - discr) / 2. / dq2;
          Real x1 = (-dq1 + discr) / 2. / dq2;
          Real qx0 = ((q3 * x0 + q2) * x0 + q1) * x0 + q0;
          Real qx1 = ((q3 * x1 + q2) * x1 + q1) * x1 + q0;
          //std::cout<<" fcubic x0="<<x0<<" qx0="<<qx0<<" x1="<<x1<<" qx1="<<qx1<<std::endl;
          if (q3 > 0) { //grows towards infinity, x0 is maximum
            if ((x0 >= 1.) || ((x0 > 0.) && (qx0 >= 0.))) {
              ubstep = min(1., infl);
              stepval = 0.;
              rootsearch = true;
            } else if (x1 >= 1.) {
              ubstep = 1.;
              stepval = 1.;
              rootsearch = false;
            } else {
              ubstep = 1.;
              stepval = 1.;
              rootsearch = true;
            }
          } else { //q3<0., grows towards minus infinity, x1 is maxium
            if (x0 >= 1.) {
              ubstep = 1.;
              stepval = 1.;
              rootsearch = false;
            } else if ((x1 >= 1.) || ((x1 > 0.) && (qx1 >= 0.))) {
              ubstep = min(x1, 1.);
              stepval = max(0., min(ubstep, infl));
              rootsearch = true;
            } else {
              ubstep = 1.;
              stepval = 1.;
              rootsearch = false;
            }
          }
        }
      } //end full cubic
    } //end d2q2 zero
    else if (discr <= eps_Real * maxnrmq) {
      //full 4th degree with one extremum
      //convex or concave case
      if (q4 > 0) {  //convex
        stepval = 1.;
        ubstep = 1.;
        Real qone = q4 + q3 + q2 + q1 + q0;
        Real dqone = dq3 + dq2 + dq1 + dq0;
        if (qone < 0.) {
          //value at step size 1 is feasible
          rootsearch = false;
        } else {
          //value at step size 1 is insfeasible
          if (dqone > 0) {
            //positive slope, root search will find zero 
            rootsearch = true;
          } else {
            //negative slope up to step size 1, so 1 is the least infeasible choice
            rootsearch = false;
          }
        }
        //std::cout<<" cvx discr="<<discr<<" q(1)="<<qone<<" q'(1)="<<dqone<<" search="<<rootsearch<<std::endl;
      } else {
        //std::cout<<" ccv discr="<<discr<<", should not happen"<<std::endl;
        // q4<0, concave
        if (dq0 <= 0.) {
          //negative slope in dq0 will stay so to step size 1, 1 is the best choice
          stepval = 1.;
          ubstep = 1.;
          rootsearch = false;
        } else {
          //positive solpe in 0 find the maximum, i.e. the zero of the derivative
          Real x = 0;
          Real der = dq0;
          //Newton should not overshoot by strict concavity
          while ((x < 1. - 1e-6) && (der > 1e-6 * maxnrmq)) {
            x += -der / ((d2q2 * x + d2q1) * x + d2q0);
            x = min(x, 1.);
            der = ((dq3 * x + dq2) * x + dq1) * x + dq0;
          }
          Real qx = (((q4 * x + q3) * x + q2) * x + q1) * x + q0;
          if (qx >= 0.) {
            //the max is positive, start searching for value 0 from step size 0
            ubstep = min(1., qx);
            stepval = 0;
            rootsearch = true;
          } else {
            //the max is negative [0,1], so take the maxium step 
            stepval = 1.;
            ubstep = 1.;
            rootsearch = false;
          }
        }
      }
    } else {
      // two inflection points
      discr = std::sqrt(discr);
      Real infl0 = (-d2q1 - discr) / 2. / d2q2;
      Real infl1 = (-d2q1 + discr) / 2. / d2q2;
      if (d2q2 < 0.)
        swap(infl0, infl1);
      Real qinfl0 = (((q4 * infl0 + q3) * infl0 + q2) * infl0 + q1) * infl0 + q0;
      Real qinfl1 = (((q4 * infl1 + q3) * infl1 + q2) * infl1 + q1) * infl1 + q0;
      Real d1infl0 = ((dq3 * infl0 + dq2) * infl0 + dq1) * infl0 + dq0;
      Real d1infl1 = ((dq3 * infl1 + dq2) * infl1 + dq1) * infl1 + dq0;

      // Real d2infl0=(d2q2*infl0+d2q1)*infl0+d2q0;
      // Real d2infl1=(d2q2*infl1+d2q1)*infl1+d2q0;

      //find the roots of the derivative which is a cubic 
      Real d0 = dq2 * dq2 - 3 * dq3 * dq1;
      Real d1 = 2 * dq2 * dq2 * dq2 - 9 * dq3 * dq2 * dq1 + 27 * dq3 * dq3 * dq0;
      Real dprime = d1 * d1 - 4 * d0 * d0 * d0;
      Real r_dp = 0.;
      Real phi_dp = 0.;
      if (dprime > 0) {
        r_dp = sqrt(dprime);
        if (d1 < 0) {
          r_dp -= d1;
          phi_dp = 3.1415926535897932;
        } else {
          r_dp += d1;
        }
      } else {
        r_dp = std::sqrt(-dprime + d1 * d1);
        phi_dp = std::acos(d1 / r_dp);
      }
      r_dp /= 2.;
      r_dp = pow(r_dp, 1 / 3.);
      phi_dp /= 3.;
      Real x0 = (dq2 + std::cos(phi_dp) * r_dp + (d0 / r_dp) * std::cos(-phi_dp)) / (-3. * dq3);
      Real Im_x0 = std::sin(phi_dp) * r_dp + (d0 / r_dp) * std::sin(-phi_dp);
      phi_dp += 2. * 3.1415926535897932 / 3.;
      Real x1 = (dq2 + std::cos(phi_dp) * r_dp + (d0 / r_dp) * std::cos(-phi_dp)) / (-3. * dq3);
      Real Im_x1 = std::sin(phi_dp) * r_dp + (d0 / r_dp) * std::sin(-phi_dp);
      phi_dp += 2. * 3.1415926535897932 / 3.;
      Real x2 = (dq2 + std::cos(phi_dp) * r_dp + (d0 / r_dp) * std::cos(-phi_dp)) / (-3. * dq3);
      Real Im_x2 = std::sin(phi_dp) * r_dp + (d0 / r_dp) * std::sin(-phi_dp);
      if (x0 > x1) {
        swap(x0, x1);
        swap(Im_x0, Im_x1);
      }
      if (x2 < x0) {
        Real d = x2; x2 = x1; x1 = x0; x0 = d;
        d = Im_x2; Im_x2 = Im_x1; Im_x1 = Im_x0; Im_x0 = d;
      } else if (x2 < x1) {
        swap(x1, x2);
        swap(Im_x1, Im_x2);
      }
      Real qx0 = (((q4 * x0 + q3) * x0 + q2) * x0 + q1) * x0 + q0;
      Real qx1 = (((q4 * x1 + q3) * x1 + q2) * x1 + q1) * x1 + q0;
      Real qx2 = (((q4 * x2 + q3) * x2 + q2) * x2 + q1) * x2 + q0;
      // Real d1x0=((dq3*x0+dq2)*x0+dq1)*x0+dq0;
      // Real d1x1=((dq3*x1+dq2)*x1+dq1)*x1+dq0;
      // Real d1x2=((dq3*x2+dq2)*x2+dq1)*x2+dq0;
      // std::cout.precision(12);
      // std::cout<<" x0="<<x0<<" ("<<qx0<<",i"<<Im_x0<<","<<d1x0<<")";
      // std::cout<<" infl0="<<infl0<<" ("<<qinfl0<<" "<<d1infl0<<","<<d2infl0<<")";
      // std::cout<<" x1="<<x1<<" ("<<qx1<<",i"<<Im_x1<<","<<d1x1<<")";
      // std::cout<<" infl1="<<infl1<<" ("<<qinfl1<<" "<<d1infl1<<","<<d2infl1<<")";
      // std::cout<<" x2="<<x2<<" ("<<qx2<<",i"<<Im_x2<<","<<d1x2<<")";
      // std::cout<<std::endl;

      // assert((std::fabs(Im_x0)>1e-8)||(x0<=infl0+1e-10*(1.+std::fabs(x0))));
      // assert((std::fabs(Im_x1)>1e-8)||(infl0<=x1+1e-10*(1.+std::fabs(infl0))));
      // assert((std::fabs(Im_x1)>1e-8)||(x1<=infl1+1e-10*(1.+std::fabs(x1))));
      // assert((std::fabs(Im_x2)>1e-8)||(infl1<=x2+1e-10*(1.+std::fabs(infl1))));


      if (q4 > 0) {
        if (d1infl0 <= 0) { //both infl have negative deriv, min is in x2
          if (x2 >= 1.) {
            stepval = 1.;
            ubstep = 1.;
            rootsearch = false;
          } else {
            stepval = 1.;
            ubstep = 1.;
            rootsearch = true;
          }
        } else if (d1infl1 >= 0) { //both infl have positive deriv, min is in x0
          if (x0 >= 1.) {
            stepval = 1.;
            ubstep = 1.;
            rootsearch = false;
          } else if ((infl1 >= 1.) || ((infl1 > 0.) && (qinfl1 >= 0.))) {
            ubstep = min(1., infl1);
            stepval = max(0., min(ubstep, infl0));
            rootsearch = true;
          } else {
            stepval = 1.;
            ubstep = 1.;
            rootsearch = true;
          }
        } else { //there is a max in the middle
          assert(std::fabs(Im_x0) < 1e-8 * maxnrmq);
          assert(std::fabs(Im_x1) < 1e-8 * maxnrmq);
          assert(std::fabs(Im_x2) < 1e-8 * maxnrmq);
          if (x0 >= 1.) {
            stepval = 1.;
            ubstep = 1.;
            rootsearch = false;
          } else if ((x1 >= 1.) || ((x1 > 0.) && (qx1 >= 0.))) {
            ubstep = min(x1, 1.);
            stepval = max(0., min(ubstep, infl0));
            rootsearch = true;
          } else if (x2 >= 1.) {
            stepval = 1.;
            ubstep = 1.;
            rootsearch = false;
          } else {
            stepval = 1.;
            ubstep = 1.;
            rootsearch = true;
          }
        }
      } else {  //q4<0
        if (d1infl1 <= 0) { //both infl have negative derivatives, max in x0
          if ((x0 >= 1.) || ((x0 > 0.) && (qx0 > 0.))) {
            ubstep = min(1., x0);
            stepval = 0.;
            rootsearch = true;
          } else {
            ubstep = 1.;
            stepval = 1.;
            rootsearch = false;
          }
        } else if (d1infl0 >= 0) {//both infl have positive derivatives, max in x2
          if ((infl0 >= 1.) || ((infl0 > 0.) && (qinfl0 >= 0.))) {
            ubstep = min(1., infl0);
            stepval = 0.;
            rootsearch = true;
          } else if ((x2 >= 1.) || ((x2 > 0.) && (qx2 >= 0.))) {
            ubstep = min(1., x2);
            stepval = max(0., min(ubstep, infl1));
            rootsearch = true;
          } else {
            ubstep = 1.;
            stepval = 1.;
            rootsearch = false;
          }
        } else { //there is a minimum between the inflection points
          assert(std::fabs(Im_x0) < 1e-8 * maxnrmq);
          assert(std::fabs(Im_x1) < 1e-8 * maxnrmq);
          assert(std::fabs(Im_x2) < 1e-8 * maxnrmq);
          if ((x0 >= 1.) || ((x0 > 0.) && (qx0 >= 0.))) {
            ubstep = min(1., x0);
            stepval = 0.;
            rootsearch = true;
          } else if (x1 >= 1.) {
            ubstep = 1.;
            stepval = 1.;
            rootsearch = false;
          } else if ((x2 >= 1.) || ((x2 > 0.) && (qx2 >= 0.))) {
            ubstep = min(1., x2);
            stepval = max(0., min(ubstep, infl1));
            rootsearch = true;
          } else {
            ubstep = 1.;
            stepval = 1.;
            rootsearch = false;
          }
        } //end there is min in the middle
      } // end q4<0.
    } // end two inflection points

    Real qsv = (((q4 * stepval + q3) * stepval + q2) * stepval + q1) * stepval + q0;
    Real dqsv = ((dq3 * stepval + dq2) * stepval + dq1) * stepval + dq0;
    if (rootsearch) {
      Real lbstep = 0.;
      do {
        //std::cout<<" rs["<<lbstep<<","<<ubstep<<"]("<<qsv<<","<<dqsv<<")";
        assert(lbstep <= stepval);
        assert(stepval <= ubstep + 1.e-10);
        assert(ubstep <= 1. + 1.e-10);
        if (qsv <= 0.) {
          lbstep = stepval;
          if (dqsv > 0.) {
            stepval += min(-qsv / dqsv, 0.999 * (ubstep - stepval));
          } else {
            stepval += 0.999 * (ubstep - stepval);
          }
        } else {
          if (dqsv > 0.) {
            ubstep = stepval;
            stepval += max(-qsv / dqsv, 0.999 * (lbstep - stepval));
          } else {
            if (cb_out()) {
              get_out() << "**** WARNING in QPSolverBasicStructures::QPselect_mu(): root search should never get here,";
              get_out() << " lbstep=" << lbstep << " stepval=" << stepval << " ubstep=" << ubstep << " qsv=" << qsv << " dqsv=" << dqsv << std::endl;
            }
            break;
          }
        }
        qsv = (((q4 * stepval + q3) * stepval + q2) * stepval + q1) * stepval + q0;
        dqsv = ((dq3 * stepval + dq2) * stepval + dq1) * stepval + dq0;
      } while ((std::fabs(qsv) > 1e-8 * omu) && (ubstep - lbstep > 1e-8));
      //std::cout<<" frs["<<lbstep<<","<<ubstep<<"]("<<qsv<<","<<dqsv<<")";
    }

    Real smu = (tr_xz + stepval * (tr_xdzpdxz + stepval * tr_dxdz)) / mudim;
    Real sqrnrm = nrmsqr_xz + stepval * (2 * ip_xz_xdzpdxz + stepval * (nrmsqr_xdzpdxz + 2 * ip_xz_dxdz + stepval * (2 * ip_dxdz_xdzpdxz + stepval * nrmsqr_dxdz)));
    if (((smu < 0.) || (sqrnrm < 0.)) && (cb_out())) {
      get_out() << "**** WARNING in QPSolverBasicStructures::QPselect_step_mu(...): squared distance to central path at step computes to " << sqrnrm << " with mu =" << smu << " (both need to be positive)" << std::endl;
    }
    if ((origq) && (std::fabs(qsv - (sqrnrm - sqr(theta * smu))) > 1e-6 * omu) && (cb_out())) {
      get_out().precision(12);
      get_out() << "\n trxz=" << tr_xz << " trdxpz=" << tr_xdzpdxz << " trdxdz=" << tr_dxdz << std::endl;
      get_out() << " n2xz=" << nrmsqr_xz << " n2dxpz=" << nrmsqr_xdzpdxz << " n2dxdz=" << nrmsqr_dxdz << std::endl;
      get_out() << " ixzdxpz=" << ip_xz_xdzpdxz << " ixzdxdz=" << ip_xz_dxdz << " idxpzdxdz=" << ip_dxdz_xdzpdxz << std::endl;
      get_out() << " s-q0=" << q0 << " q1=" << q1 << " q2=" << q2 << " q3=" << q3 << " q4=" << q4 << std::endl;
      get_out() << " st=" << stepval << " sqrnrm=" << sqrnrm << " smu=" << smu << " theta=" << theta << " qsv=" << qsv << " diff=" << qsv - (sqrnrm - sqr(theta * smu)) << std::endl;
    }
    if ((q0 < 0.) && (qsv > 1e-6 * max(1., omu))) {
      assert((q0 > 0) || (qsv < 1e-6 * max(1., omu)));
    }
    sqrnrm = max(sqrnrm, eps_Real * tr_xz);
    if ((stepval > 1 - 1e-3) && ((smu < 1e-6 * omu) || (std::sqrt(sqrnrm) > 1.01 * nbh_ubnd * smu))) {
      //try to avoid numerical difficulties
      if (std::fabs(tr_dxdz) < 1e-4 * tr_xz) {
        if (tr_xdzpdxz < -.9 * tr_xz)
          stepval = min(1 - 1e-2, -(1 - 1e-2) * tr_xz / tr_xdzpdxz);
        else
          stepval = 1 - 1e-2;
      } else {
        Real q_p = tr_xdzpdxz / tr_dxdz / 2.;
        Real q_q = (1 - 1e-2) * tr_xz / tr_dxdz;
        Real d = sqr(q_p) - q_q;
        if (d < 0.) {
          stepval = 1 - 1e-2;
        } else {
          if (tr_dxdz > 0.)
            stepval = min(1 - 1e-2, -q_p - std::sqrt(d));
          else
            stepval = min(1 - 1e-2, -q_p + std::sqrt(d));
        }
      }
      assert((stepval > 0.) && (stepval <= 1.));
      smu = (tr_xz + stepval * (tr_xdzpdxz + stepval * tr_dxdz)) / mudim;
      sqrnrm = nrmsqr_xz + stepval * (2 * ip_xz_xdzpdxz + stepval * (nrmsqr_xdzpdxz + 2 * ip_xz_dxdz + stepval * (2 * ip_dxdz_xdzpdxz + stepval * nrmsqr_dxdz)));
    }

    smu = max(smu, eps_Real * tr_xz);
    sqrnrm = max(sqrnrm, eps_Real * tr_xz);
    Real stheta = std::sqrt(sqrnrm) / smu;

    //--- if theta is too large, search for a step size that minimizes theta in [0,stepval]
    if ((otheta > nbh_ubnd + 1e-6) && (stepval > 1e-2) && (otheta - stheta < 0.01 * stepval * (otheta - nbh_ubnd))) {
      q0 = nrmsqr_xz - sqr(theta * omu);
      assert(q0 <= 1e-8 * omu);
      q1 = 2. * ip_xz_xdzpdxz - 2 * sqr(theta / mudim) * tr_xz * tr_xdzpdxz;
      q2 = nrmsqr_xdzpdxz + 2. * ip_xz_dxdz - sqr(theta / mudim) * (sqr(tr_xdzpdxz) + 2 * tr_xz * tr_dxdz);
      q3 = 2. * ip_dxdz_xdzpdxz - 2 * sqr(theta / mudim) * tr_xdzpdxz * tr_dxdz;
      q4 = nrmsqr_dxdz - sqr(theta * tr_dxdz / mudim);
      if (cb_out(3)) {
        get_out() << " stheta(" << stepval << ")=" << stheta << ">" << nbh_ubnd;
        get_out() << " th-q0=" << q0 << " q1=" << q1 << " q2=" << q2 << " q3=" << q3 << " q4=" << q4 << std::endl;
      }
      //ploynomial coefficients of the derivative dq, need to find zeros of this in [0, stepval]
      dq0 = q1;
      dq1 = 2 * q2;
      dq2 = 3 * q3;
      dq3 = 4 * q4;
      //ploynomial coefficients of the second derivative d2q for testing minima and maxima
      d2q0 = dq1;
      d2q1 = 2 * dq2;
      d2q2 = 3 * dq3;

      Real bestval = q0;
      Real bestx = 0.;
      Real valstep = q0 + stepval * (q1 + stepval * (q2 + stepval * (q3 + stepval * q4)));
      if (valstep < bestval) {
        bestval = valstep;
        bestx = stepval;
      }

      //find the roots of the derivative which is a cubic 
      Real d0 = dq2 * dq2 - 3 * dq3 * dq1;
      Real d1 = 2 * dq2 * dq2 * dq2 - 9 * dq3 * dq2 * dq1 + 27 * dq3 * dq3 * dq0;
      Real dprime = d1 * d1 - 4 * d0 * d0 * d0;
      Real r_dp = 0.;
      Real phi_dp = 0.;
      if (dprime > 0) {
        r_dp = sqrt(dprime);
        if (d1 < 0) {
          r_dp -= d1;
          phi_dp = 3.1415926535897932;
        } else {
          r_dp += d1;
        }
      } else {
        r_dp = std::sqrt(-dprime + d1 * d1);
        phi_dp = std::acos(d1 / r_dp);
      }
      r_dp /= 2.;
      r_dp = pow(r_dp, 1 / 3.);
      phi_dp /= 3.;
      Real x0 = (dq2 + std::cos(phi_dp) * r_dp + (d0 / r_dp) * std::cos(-phi_dp)) / (-3. * dq3);
      Real Im_x0 = std::sin(phi_dp) * r_dp + (d0 / r_dp) * std::sin(-phi_dp);
      phi_dp += 2. * 3.1415926535897932 / 3.;
      Real x1 = (dq2 + std::cos(phi_dp) * r_dp + (d0 / r_dp) * std::cos(-phi_dp)) / (-3. * dq3);
      Real Im_x1 = std::sin(phi_dp) * r_dp + (d0 / r_dp) * std::sin(-phi_dp);
      phi_dp += 2. * 3.1415926535897932 / 3.;
      Real x2 = (dq2 + std::cos(phi_dp) * r_dp + (d0 / r_dp) * std::cos(-phi_dp)) / (-3. * dq3);
      Real Im_x2 = std::sin(phi_dp) * r_dp + (d0 / r_dp) * std::sin(-phi_dp);
      if (x0 > x1) {
        swap(x0, x1);
        swap(Im_x0, Im_x1);
      }
      if (x2 < x0) {
        Real d = x2; x2 = x1; x1 = x0; x0 = d;
        d = Im_x2; Im_x2 = Im_x1; Im_x1 = Im_x0; Im_x0 = d;
      } else if (x2 < x1) {
        swap(x1, x2);
        swap(Im_x1, Im_x2);
      }
      Real qx0 = (((q4 * x0 + q3) * x0 + q2) * x0 + q1) * x0 + q0;
      if ((0 < x0) && (x0 < stepval) && (qx0 < bestval)) {
        bestval = qx0;
        bestx = x0;
      }
      Real qx1 = (((q4 * x1 + q3) * x1 + q2) * x1 + q1) * x1 + q0;
      if ((0 < x1) && (x1 < stepval) && (qx1 < bestval)) {
        bestval = qx1;
        bestx = x1;
      }
      Real qx2 = (((q4 * x2 + q3) * x2 + q2) * x2 + q1) * x2 + q0;
      if ((0 < x2) && (x2 < stepval) && (qx2 < bestval)) {
        bestval = qx2;
        bestx = x2;
      }
      stepval = bestx;
      qsv = bestval;
      dqsv = dq0 + stepval * (dq1 + stepval * (dq2 + stepval * dq3));
      smu = (tr_xz + stepval * (tr_xdzpdxz + stepval * tr_dxdz)) / mudim;
      sqrnrm = nrmsqr_xz + stepval * (2 * ip_xz_xdzpdxz + stepval * (nrmsqr_xdzpdxz + 2 * ip_xz_dxdz + stepval * (2 * ip_dxdz_xdzpdxz + stepval * nrmsqr_dxdz)));
      if (((smu < 0.) || (sqrnrm < 0.)) && (cb_out())) {
        get_out() << "**** WARNING in QPSolverBasicStructures::QPselect_step_mu(...): squared distance to central path at step computes to " << sqrnrm << " with mu =" << smu << " (both need to be positive)" << std::endl;
      }
      smu = max(smu, eps_Real * tr_xz);
      sqrnrm = max(sqrnrm, eps_Real * tr_xz);
      stheta = std::sqrt(sqrnrm) / smu;
    }


    next_theta = stheta;
    if (cb_out(3)) {
      get_out() << " stepval=" << stepval << " qsv=" << qsv << " dqsv=" << dqsv;
      get_out() << " smu=" << smu << " stheta=" << stheta << std::endl;
    }

    stepsize = stepval;
    if ((otheta > nbh_ubnd + 1e-6) && (stepsize < 0.01)) {
      if ((gmu > 0.9 * omu) && (q1 > 0.))
        sigma = 1.1;
      else
        sigma = 1.;
    } else {
      if (gmu < 1e-6 * omu) {
        sigma = 1. - 0.35 / std::sqrt(mudim);
      } else {
        sigma = gmu / omu;
      }
      sigma = std::pow(sigma, min(4., stepval * nbh_lbnd / stheta));
    }
    if (stheta < nbh_ubnd)
      sigma = min(sigma, 1. - 0.35 / std::sqrt(mudim));
    mu = sigma * smu;

    if (cb_out(3))
      get_out() << " sigma=" << sigma << " nextmu=" << mu << std::endl;

    return 0;
  }

  // *************************************************************************
  //                            QPselect_mu
  // *************************************************************************

  int QPSolverBasicStructures::QPselect_mu(Real& mu, Real stepsize) {

    Integer mudim = 0;
    Real globminval = max_Real;
    Real globmaxval = min_Real;
    Real minv, maxv;
    Real tr_xz = 0.;
    Real tr_xdzpdxz = 0.;
    Real tr_dxdz = 0.;

    if (slacklb.rowdim() > 0) {
      mudim += slacklb.rowdim();
      tr_xz += ip_min_max(slacklb, shortzlb, minv, maxv);
      tr_xdzpdxz += ip(slacklb, dshortzlb) + ip(dslacklb, shortzlb);
      tr_dxdz += ip(dslacklb, dshortzlb);

      if (minv < globminval)
        globminval = minv;
      if (maxv > globmaxval)
        globmaxval = maxv;

      if (cb_out(3))
        get_out() << " slb[" << minv << "," << maxv << "]";
    }

    if (slackub.rowdim() > 0) {
      mudim += slackub.rowdim();
      tr_xz += ip_min_max(slackub, shortzub, minv, maxv);
      tr_xdzpdxz += ip(slackub, dshortzub) + ip(dslackub, shortzub);
      tr_dxdz += ip(dslackub, dshortzub);

      if (minv < globminval)
        globminval = minv;
      if (maxv > globmaxval)
        globmaxval = maxv;

      if (cb_out(3))
        get_out() << " sub[" << minv << "," << maxv << "]";
    }

    if (rhsslacklb.rowdim() > 0) {
      mudim += rhsslacklb.rowdim();
      tr_xz += ip_min_max(rhsslacklb, shortrhszlb, minv, maxv);
      tr_xdzpdxz += ip(rhsslacklb, dshortrhszlb) + ip(drhsslacklb, shortrhszlb);
      tr_dxdz += ip(drhsslacklb, dshortrhszlb);

      if (minv < globminval)
        globminval = minv;
      if (maxv > globmaxval)
        globmaxval = maxv;

      if (cb_out(3))
        get_out() << " rslb[" << minv << "," << maxv << "]";
    }

    if (rhsslackub.rowdim() > 0) {
      mudim += rhsslackub.rowdim();
      tr_xz += ip_min_max(rhsslackub, shortrhszub, minv, maxv);
      tr_xdzpdxz += ip(rhsslackub, dshortrhszub) + ip(drhsslackub, shortrhszub);
      tr_dxdz += ip(drhsslackub, dshortrhszub);

      if (minv < globminval)
        globminval = minv;
      if (maxv > globmaxval)
        globmaxval = maxv;

      if (cb_out(3))
        get_out() << " rsub[" << minv << "," << maxv << "]";
    }


    if (model_block) {
      if (model_block->get_mu_info(mudim, tr_xz, tr_xdzpdxz, tr_dxdz, globminval, globmaxval)) {
        if (cb_out()) {
          get_out() << "**** ERROR in QPSolverBasicStructures::QPselect_mu(...): model_block->get_mu_info(...) failed" << std::endl;
        }
      }
    }

    if (paramsp->QPget_use_socqp()) {
      if (socqp.get_mu_info(mudim, tr_xz, tr_xdzpdxz, tr_dxdz, globminval, globmaxval)) {
        if (cb_out()) {
          get_out() << "**** ERROR in QPSolverBasicStructures::QPselect_mu(...): socqp.get_mu_info(...) failed" << std::endl;
        }
      }
    }

    if (last_mu < 0) {
      assert(mu > 0.);
      last_mu = mu;
    }

    Real thresholdmu = min(mu, last_mu);

    Real omu = tr_xz / mudim;
    Real nmu = (tr_xz + stepsize * (tr_xdzpdxz + stepsize * tr_dxdz)) / mudim;

    if (paramsp->QPget_use_predictor_corrector())
      mu = omu;
    else
      mu = max(eps_Real, nmu);

    if (cb_out(3)) {
      get_out() << " gmu=" << mu;
    }



    //if (n2dualviol>max(mu,paramsp->QPget_dual_infeasibility_eps())){
    //  sigma=1.;
    //  thresholdmu*=max(.99,1.-1./std::sqrt(mudim));
    //}
    //else {

      //sigma=min(.99,pow(ipdsz/ipsz,max(2.,3*sqr(stepsize))));
      //sigma=min(.99,pow(ipdsz/ipsz,max(2.,3*sqr(min(stepsize,globminval/globmaxval)))));
    if ((stepsize < 0.25) || (globminval / globmaxval < 0.3)) {
      sigma = max(0.99, 1 - 1. / std::sqrt(mudim));
      //mu=last_mu;
    } else
      sigma = max(sigma / 2., min(.99, pow(nmu / omu, max(2., 4. * sqr(min(stepsize, globminval / globmaxval))))));
    if ((stepsize > .95) && (globminval / globmaxval > .5))
      sigma = min(sigma, max(0.1, 1. - globminval / globmaxval));
    if ((paramsp->QPget_use_predictor_corrector()) &&
      (last_alpha > .9) &&
      (last_mu < 1.1 * mu) &&
      (nmu > .5 * omu) &&
      (globminval / globmaxval > 0.5)) {
      thresholdmu *= 0.1;
      large_predictor_cnt++;
      if (cb_out(3)) {
        get_out() << " large_pred(" << large_predictor_cnt << ")";
      }
    } else {
      thresholdmu *= max(0.99, 1 - 1. / std::sqrt(mudim));
    }
    //}

    if (cb_out(3)) {
      get_out() << " omu=" << omu << " (" << globminval / globmaxval;
      if (globminval < 0.1 * globmaxval) {
        get_out() << "[" << globminval << "," << globmaxval << "]";
      }
      get_out() << ") nmu=" << nmu << " (" << stepsize << "," << large_predictor_cnt << ") mu=" << mu << " sigma=" << sigma << " tmu=" << thresholdmu;
    }

    mu = sigma * mu;
    // if ((sigma<1.)||(sigma*mu<10.*thresholdmu)) 
    //  mu=min(thresholdmu,mu);

    if (cb_out(3))
      get_out() << " nextmu=" << mu << std::endl;

    return 0;
  }

  // *************************************************************************
  //                            QPlinesearch
  // *************************************************************************

  void QPSolverBasicStructures::QPvec_linesearch(Real& alpha,
    const Matrix& vec,
    const Matrix& dvec) const {
    const Real* vp;
    const Real* dvp;

    vp = vec.get_store();
    dvp = dvec.get_store();
    for (Integer i = vec.dim(); --i >= 0;) {
      Real d = *dvp++;
      if (d < 0.) {
        alpha = min(alpha, -(*vp++) / d);
      } else vp++;
    }
  }

  // *************************************************************************
  //                            QPlinesearch
  // *************************************************************************

  void QPSolverBasicStructures::QPlinesearch(Real& alpha) {
    alpha = 2.;

    QPvec_linesearch(alpha, slacklb, dslacklb);
    QPvec_linesearch(alpha, slackub, dslackub);
    QPvec_linesearch(alpha, shortzlb, dshortzlb);
    QPvec_linesearch(alpha, shortzub, dshortzub);
    QPvec_linesearch(alpha, rhsslacklb, drhsslacklb);
    QPvec_linesearch(alpha, rhsslackub, drhsslackub);
    QPvec_linesearch(alpha, shortrhszlb, dshortrhszlb);
    QPvec_linesearch(alpha, shortrhszub, dshortrhszub);

    if (model_block) {
      if (model_block->linesearch(alpha)) {
        if (cb_out()) {
          get_out() << "**** ERROR in QPSolverBasicStructures::QPlinesearch(...): model_block->linesearch(.) failed" << std::endl;
        }
      }
    }

    if (paramsp->QPget_use_socqp()) {
      if (socqp.linesearch(alpha)) {
        if (cb_out()) {
          get_out() << "**** ERROR in QPSolverBasicStructures::QPlinesearch(...): socqp.lineserach(.) failed" << std::endl;
        }
      }
    }

    if (alpha > 1.)
      alpha = 1.;
    else
      alpha *= 0.95;


    //TEST OUTPUT
    /*
    if (alpha<0.1){
      Matrix tmp;
      tmp.init(slacklb);
      tmp.concat_right(dslacklb);
      tmp.concat_right(shortzlb);
      tmp.concat_right(dshortzlb);
      std::cout<<" ls_lb "<<tmp;
      tmp.init(slackub);
      tmp.concat_right(dslackub);
      tmp.concat_right(shortzub);
      tmp.concat_right(dshortzub);
      std::cout<<" ls_ub "<<tmp;
      tmp.init(rhsslacklb);
      tmp.concat_right(drhsslacklb);
      tmp.concat_right(shortrhszlb);
      tmp.concat_right(dshortrhszlb);
      std::cout<<" ls_rhslb "<<tmp;
      tmp.init(rhsslackub);
      tmp.concat_right(drhsslackub);
      tmp.concat_right(shortrhszub);
      tmp.concat_right(dshortrhszub);
      std::cout<<" ls_rhsub "<<tmp;
    }
    */
  }

  // *************************************************************************
  //                           predcorr_step
  // *************************************************************************

  // do a predictor corrector step

  int QPSolverBasicStructures::QPpredcorr_step(Real& alpha,
    Real& prec,
    Real& next_mu,
    bool use_predcorr,
    bool use_nbh,
    bool centering) {
    assert((paramsp) && (paramsp->QPget_KKTsolver()));

    int status = 0;

    //--- prepare system for solving
    CH_Tools::Microseconds time_initKKT = clock.time();

    Matrix KKTdiagx(QPget_xdim(), 1, 0.);
    for (Integer i = 0; i < shortzub.dim(); i++) {
      Integer ind = QPget_ubind()(i);
      KKTdiagx(ind) += shortzub(i) / slackub(i);
      dualviol(ind) -= shortzub(i);
    }
    for (Integer i = 0; i < shortzlb.dim(); i++) {
      Integer ind = QPget_lbind()(i);
      KKTdiagx(ind) += shortzlb(i) / slacklb(i);
      dualviol(ind) += shortzlb(i);
    }
    Matrix KKTdiagy(QPget_ydim(), 1, 0.);
    Matrix dsrhs(QPget_ydim(), 1, 0.);
    for (Integer i = 0; i < shortrhszub.dim(); i++) {
      Integer ind = QPget_rhslbind()(i);
      KKTdiagy(ind) += shortrhszub(i) / rhsslackub(i);
    }
    for (Integer i = 0; i < shortrhszlb.dim(); i++) {
      Integer ind = QPget_rhsubind()(i);
      KKTdiagy(ind) += shortrhszlb(i) / rhsslacklb(i);
    }
    for (Integer i = 0; i < KKTdiagy.rowdim(); i++) {
      const Real d = KKTdiagy(i);
      assert(d >= 0.);
      if (d > 0.) {
        KKTdiagy(i) = 1. / d;
        primalviol(i) -= y(i) / d;
        dsrhs(i) = y(i);
      }
    }

    //--- compute predictor rhs
    primalviol *= -1.;
    dualviol *= -1.;

    if (!use_predcorr) {
      for (Integer i = 0; i < slacklb.dim(); i++) {
        dualviol(QPget_lbind()(i)) += next_mu / slacklb(i);
      }
      for (Integer i = 0; i < slackub.dim(); i++) {
        dualviol(QPget_ubind()(i)) -= next_mu / slackub(i);
      }
      for (Integer i = 0; i < rhsslacklb.dim(); i++) {
        Integer ind = QPget_rhsubind()(i);
        primalviol(ind) -= next_mu / rhsslacklb(i) * KKTdiagy(ind);
        dsrhs(ind) -= next_mu / rhsslacklb(i);
      }
      for (Integer i = 0; i < rhsslackub.dim(); i++) {
        Integer ind = QPget_rhslbind()(i);
        primalviol(ind) += next_mu / rhsslackub(i) * KKTdiagy(ind);
        dsrhs(ind) += next_mu / rhsslackub(i);
      }
    }

    Matrix primal_rhs(primalviol);
    Matrix dual_rhs(dualviol);
    Real Hfactor = 1.;

    if (paramsp->QPget_use_socqp()) {
      if (socqp.add_prox_sysrhs(dual_rhs, Hfactor, use_predcorr ? 0. : next_mu, 0.)) {
        if (cb_out()) {
          get_out() << "**** ERROR in QPSolverBasicStructures::QPpredcorr_step(...): socqp.add_prox_sysrhs(...) failed" << std::endl;
        }
      }
    }

    if (paramsp->QPget_KKTsolver()->QPinit_KKTsystem(KKTdiagx, KKTdiagy, Hfactor, prec, paramsp)) {
      if (cb_out(1)) {
        get_out() << "*** WARNING: QPSolverBasicStructures::QPpredcorr_step(): QPinit_system() failed " << std::endl;
      }
      status = 3;
      QPinitKKT_time += clock.time() - time_initKKT;
      return status;
    }
    QPinitKKT_time += clock.time() - time_initKKT;

    CH_Tools::Microseconds time_predictor = clock.time();


    //Matrix solx;
    //Matrix soly;

    //std::cout<<" pc="<<use_predcorr<<" mu="<<(use_predcorr?0.:next_mu)<<std::endl;
    if (paramsp->QPget_KKTsolver()->QPsolve_KKTsystem(solx, soly, primal_rhs, dual_rhs, use_predcorr ? 0. : next_mu, 0., prec, paramsp)) {
      if (cb_out(1)) {
        get_out() << "*** WARNING: QPSolverBasicStructures::QPpredcorr_step(): QPsolve_system() failed for predictor" << std::endl;
      }
      status = 2;
      if (use_predcorr)
        QPpredictor_time += clock.time() - time_predictor;
      else
        QPcorrector_time += clock.time() - time_predictor;
      return status;
    }

    //--- compute predictor step

    if (paramsp->QPget_use_socqp()) {
      if (socqp.compute_step(solx)) {
        if (cb_out()) {
          get_out() << "**** ERROR in QPSolverBasicStructures::QPpredcorr_step(...): socqp.compute_step(.) failed for predictor" << std::endl;
        }
      }
    }

    Matrix tmpslackinv;
    dslacklb.init(solx(QPget_lbind()));
    dslackub.init(solx(QPget_ubind()), -1.);
    dshortzlb.init(slacklb, -1.);
    dshortzlb.inv();
    if (!use_predcorr)
      tmpslackinv.init(dshortzlb, -next_mu);
    dshortzlb %= dslacklb;
    dshortzlb -= 1;
    dshortzlb %= shortzlb;
    if (!use_predcorr)
      dshortzlb += tmpslackinv;
    dshortzub.init(slackub, -1.);
    dshortzub.inv();
    if (!use_predcorr)
      tmpslackinv.init(dshortzub, -next_mu);
    dshortzub %= dslackub;
    dshortzub -= 1;
    dshortzub %= shortzub;
    if (!use_predcorr)
      dshortzub += tmpslackinv;
    Matrix ds(soly, -1.);
    ds -= dsrhs;
    ds %= KKTdiagy;
    drhsslacklb.init(ds(QPget_rhsubind()));
    drhsslackub.init(ds(QPget_rhslbind()), -1.);
    dshortrhszlb.init(rhsslacklb, -1.);
    dshortrhszlb.inv();
    if (!use_predcorr)
      tmpslackinv.init(dshortrhszlb, -next_mu);
    dshortrhszlb %= drhsslacklb;
    dshortrhszlb -= 1.;
    dshortrhszlb %= shortrhszlb;
    if (!use_predcorr)
      dshortrhszlb += tmpslackinv;
    dshortrhszub.init(rhsslackub, -1.);
    dshortrhszub.inv();
    if (!use_predcorr)
      tmpslackinv.init(dshortrhszub, -next_mu);
    dshortrhszub %= drhsslackub;
    dshortrhszub -= 1.;
    dshortrhszub %= shortrhszub;
    if (!use_predcorr)
      dshortrhszub += tmpslackinv;

    if (use_predcorr)
      QPpredictor_time += clock.time() - time_predictor;
    else
      QPcorrector_time += clock.time() - time_predictor;


    // //TESTING
    // std::cout<<" predictor-step, use_pred="<<use_predcorr<<" prec="<<prec<<std::endl;
    // {
    //   Matrix tmpvec(dual_rhs,-1.);
    //   QPadd_Aty(soly,tmpvec);
    //   if ((model_block)&&(model_block->dim_model()>0)){
    //     model_block->B_times(model_block->get_dx(),tmpvec,1.,1.,1,0);
    //   }
    //   QPadd_Qx(solx*Hfactor,tmpvec);
    //   std::cout<<" Hblock= "<<norm2(tmpvec)<<std::endl;
    // }
    // Matrix testdualinf(QPget_c()); 
    // QPadd_Aty(y,testdualinf);
    // QPadd_Aty(soly,testdualinf);
    // for(Integer i=0;i<shortzub.dim();i++)
    //   testdualinf(QPget_ubind()(i))+=shortzub(i)+dshortzub(i);
    // for(Integer i=0;i<shortzlb.dim();i++)
    //   testdualinf(QPget_lbind()(i))-=(shortzlb(i)+dshortzlb(i));
    // if ((model_block)&&(model_block->dim_model()>0)){
    //   //std::cout<<"c+Q(x+dx)="<<testdualinf;
    //   Matrix modelx=model_block->get_x()+model_block->get_dx();
    //   //std::cout<<"modelx+modeldx="<<modelx;
    //   model_block->B_times(modelx,testdualinf,1.,1.,1,0);

    //   //Real dummy;
    //   //for (Integer i=0;i<Integer(model_block->get_bundle().size());i++){
    //   //  model_block->get_bundle()[unsigned(i)].get_minorant(dummy,testdualinf,0,modelx(i),true);;
    //   //}

    // }
    // if (paramsp->QPget_use_socqp()){
    //   socqp.test_sysviol(testdualinf,solx);
    // }
    // else {
    //   QPadd_Qx(x,testdualinf);
    //   QPadd_Qx(solx,testdualinf);
    // }
    // Integer minind,maxind;
    // std::cout<<" norm2(Q*(x+dx)+c+transpose(A)*(y+dy)+transpose(B)*(modx+moddx)-zlb-dzlb+zub+dzub)="<<norm2(testdualinf)<<" ["<<min(testdualinf,&minind)<<"("<<minind<<"),"<<max(testdualinf,&maxind)<<"("<<maxind<<")]"<<std::endl;
    // Matrix testyinf(y.rowdim(),1,0.);
    // for(Integer i=0;i<QPget_rhslbind().dim();i++){
    //   Integer ind=QPget_rhslbind()(i);
    //   testyinf(ind)=y(ind)+soly(ind);
    // }
    // for(Integer i=0;i<QPget_rhsubind().dim();i++){
    //   Integer ind=QPget_rhsubind()(i);
    //   testyinf(ind)=y(ind)+soly(ind);
    // }
    // for(Integer i=0;i<shortrhszub.dim();i++)
    //   testyinf(QPget_rhslbind()(i))+=shortrhszub(i)+dshortrhszub(i);
    // for(Integer i=0;i<shortrhszlb.dim();i++)
    //   testyinf(QPget_rhsubind()(i))-=(shortrhszlb(i)+dshortrhszlb(i));
    // std::cout<<" norm2(y+dy-rzlb-drzlb+rzub+drzub)="<<norm2(testyinf)<<std::endl; 
    // Matrix testprimalinf(s);
    // testprimalinf+=ds;
    // QPadd_Ax(x,testprimalinf);
    // QPadd_Ax(solx,testprimalinf);
    // std::cout<<" norm2(A*(x+dx)+s+ds)="<<norm2(testprimalinf)<<std::endl; 
    // testdualinf.init(slacklb%shortzlb+slacklb%dshortzlb+dslacklb%shortzlb);
    // if (!use_predcorr)
    //   testdualinf-=next_mu;
    // std::cout<<" norm2(sl%zl+sl%dzl+dsl%dzl)="<<norm2(testdualinf)<<std::endl;
    // testdualinf.init(slackub%shortzub+slackub%dshortzub+dslackub%shortzub);
    // if (!use_predcorr)
    //   testdualinf-=next_mu;
    // std::cout<<" norm2(su%zu+su%dzu+dsl%dzu)="<<norm2(testdualinf)<<std::endl;
    // testprimalinf.init(rhsslacklb%shortrhszlb+rhsslacklb%dshortrhszlb+drhsslacklb%shortrhszlb);
    // if (!use_predcorr)
    //   testprimalinf-=next_mu;
    // std::cout<<" norm2(rsl%rzl+rsl%drzl+drsl%drzl)="<<norm2(testprimalinf)<<std::endl;
    // testprimalinf.init(rhsslackub%shortrhszub+rhsslackub%dshortrhszub+drhsslackub%shortrhszub);
    // if (!use_predcorr)
    //   testprimalinf-=next_mu;
    // std::cout<<" norm2(rsu%rzu+rsu%drzu+drsu%drzu)="<<norm2(testprimalinf)<<std::endl;
    // if ((model_block)&&(model_block->dim_model()>0)){
    //   std::cout<<" norm2(modelviol)="<<norm2(model_block->get_sysviol_model(solx))<<std::endl;
    //   std::cout<<" norm2(constrviol)="<<norm2(model_block->get_sysviol_constraints())<<std::endl;
    // }


    //--- compute corrector
    if (use_predcorr) {

      bool no_dxdz = false;

      //--- select new mu
      CH_Tools::Microseconds time_linesearch = clock.time();
      if (use_nbh) {
        QPselect_localstep_mu(next_mu, alpha, centering);
        if (alpha < 0.01)
          no_dxdz = true;
      } else {
        QPlinesearch(alpha);
        if (alpha < 0.25) {
          no_dxdz = true;
        } else {
          QPselect_mu(next_mu, alpha);
        }
      }

      QPlinesearch_time += clock.time() - time_linesearch;


      if (cb_out(3)) {
        get_out() << " nbh=" << use_nbh << " dxdz=" << !no_dxdz << std::endl;
      }

      CH_Tools::Microseconds time_corrector = clock.time();

      // //BEGIN TEST OUTPUT
      // Matrix old_dzub(dshortzub);
      // Matrix old_dzlb(dshortzlb);
      // Matrix old_dsub(dslackub);
      // Matrix old_dslb(dslacklb);
      // Matrix old_drzub(dshortrhszub);
      // Matrix old_drzlb(dshortrhszlb);
      // Matrix old_drsub(drhsslackub);
      // Matrix old_drslb(drhsslacklb);
      // //END TEST OUTPUT

      time_corrector = clock.time();

      for (Integer i = 0; i < slacklb.dim(); i++) {
        dualviol(QPget_lbind()(i)) += (next_mu - (no_dxdz ? 0. : dshortzlb(i) * dslacklb(i))) / slacklb(i);
      }
      for (Integer i = 0; i < slackub.dim(); i++) {
        dualviol(QPget_ubind()(i)) -= (next_mu - (no_dxdz ? 0. : dshortzub(i) * dslackub(i))) / slackub(i);
      }
      for (Integer i = 0; i < rhsslacklb.dim(); i++) {
        Integer ind = QPget_rhsubind()(i);
        Real d = (-next_mu + (no_dxdz ? 0. : dshortrhszlb(i) * drhsslacklb(i))) / rhsslacklb(i);
        primalviol(ind) += d * KKTdiagy(ind);
        dsrhs(ind) += d;
      }
      for (Integer i = 0; i < rhsslackub.dim(); i++) {
        Integer ind = QPget_rhslbind()(i);
        Real d = (next_mu - (no_dxdz ? 0 : dshortrhszub(i) * drhsslackub(i))) / rhsslackub(i);
        primalviol(ind) += d * KKTdiagy(ind);
        dsrhs(ind) += d;
      }

      primal_rhs.init(primalviol);
      dual_rhs.init(dualviol);

      if (paramsp->QPget_use_socqp()) {
        Real dummy;
        if (socqp.add_prox_sysrhs(dual_rhs, dummy, next_mu, no_dxdz ? 0. : 1.)) {
          if (cb_out()) {
            get_out() << "**** ERROR in QPSolverBasicStructures::QPpredcorr_step(...): socqp.add_prox_sysrhs(...) failed" << std::endl;
          }
        }
        assert(dummy == Hfactor);
      }

      //keep old solution as starting point
      //std::cout<<" dxdz="<<(no_dxdz?0.:1.)<<" mu="<<next_mu<<std::endl;
      if (paramsp->QPget_KKTsolver()->QPsolve_KKTsystem(solx, soly, primal_rhs, dual_rhs, next_mu, no_dxdz ? 0. : 1., prec, paramsp)) {
        if (cb_out(1)) {
          get_out() << "*** WARNING: QPSolverBasicStructures::QPpredcorr_step(): QPsolve_system() failed for corrector" << std::endl;
        }
        status = 2;
        return status;
      }

      //--- compute step

      if (paramsp->QPget_use_socqp()) {
        if (socqp.compute_step(solx)) {
          if (cb_out()) {
            get_out() << "**** ERROR in QPSolverBasicStructures::QPpredcorr_step(...): socqp.compute_step(.) failed for corrector" << std::endl;
          }
        }
      }


      for (Integer i = 0; i < slacklb.dim(); i++) {
        Real dslb = solx(QPget_lbind()(i));
        Real slb = slacklb(i);
        dshortzlb(i) = (next_mu - (no_dxdz ? 0 : dshortzlb(i) * dslacklb(i))) / slb - (1. + dslb / slb) * shortzlb(i);
        dslacklb(i) = dslb;
      }
      for (Integer i = 0; i < slackub.dim(); i++) {
        Real dsub = -solx(QPget_ubind()(i));
        Real sub = slackub(i);
        dshortzub(i) = (next_mu - (no_dxdz ? 0. : dshortzub(i) * dslackub(i))) / sub - (1. + dsub / sub) * shortzub(i);
        dslackub(i) = dsub;
      }
      ds.init(soly, -1.);
      ds -= dsrhs;
      ds %= KKTdiagy;
      for (Integer i = 0; i < rhsslacklb.dim(); i++) {
        Integer ind = QPget_rhsubind()(i);
        Real dslb = ds(ind);
        Real slb = rhsslacklb(i);
        dshortrhszlb(i) = (next_mu - (no_dxdz ? 0. : dshortrhszlb(i) * drhsslacklb(i))) / slb - (1. + dslb / slb) * shortrhszlb(i);
        drhsslacklb(i) = dslb;
      }
      for (Integer i = 0; i < rhsslackub.dim(); i++) {
        Integer ind = QPget_rhslbind()(i);
        Real dsub = -ds(ind);
        Real sub = rhsslackub(i);
        dshortrhszub(i) = (next_mu - (no_dxdz ? 0. : dshortrhszub(i) * drhsslackub(i))) / sub - (1. + dsub / sub) * shortrhszub(i);
        drhsslackub(i) = dsub;
      }

      QPcorrector_time += clock.time() - time_corrector;

      // //TESTING

      // std::cout<<" corrector-step"<<std::endl;
      // std::cout<<"lcompl="<<norm2(slacklb%shortzlb+dslacklb%shortzlb+slacklb%dshortzlb+old_dslb%old_dzlb-next_mu)<<std::endl;
      // std::cout<<"ucompl="<<norm2(slackub%shortzub+dslackub%shortzub+slackub%dshortzub+old_dsub%old_dzub-next_mu)<<std::endl;
      // std::cout<<"rlcompl="<<norm2(rhsslacklb%shortrhszlb+drhsslacklb%shortrhszlb+rhsslacklb%dshortrhszlb+old_drslb%old_drzlb-next_mu)<<std::endl;
      // std::cout<<"rucompl="<<norm2(rhsslackub%shortrhszub+drhsslackub%shortrhszub+rhsslackub%dshortrhszub+old_drsub%old_drzub-next_mu)<<std::endl;
      // testdualinf.init(QPget_c());
      // //std::cout<<" x="<<x;
      // //std::cout<<" dx="<<solx;
      // //std::cout<<"c="<<testdualinf;
      // QPadd_Aty(y,testdualinf);
      // QPadd_Aty(soly,testdualinf);
      // for(Integer i=0;i<shortzub.dim();i++)
      //   testdualinf(QPget_ubind()(i))+=shortzub(i)+dshortzub(i);
      // for(Integer i=0;i<shortzlb.dim();i++)
      //   testdualinf(QPget_lbind()(i))-=(shortzlb(i)+dshortzlb(i));
      // if ((model_block)&&(model_block->dim_model()>0)){
      //   //std::cout<<"c+Q(x+dx)="<<testdualinf;
      //   Matrix modelx=model_block->get_x()+model_block->get_dx();
      //   //std::cout<<"modelx+modeldx="<<modelx;
      //   model_block->B_times(modelx,testdualinf,1.,1.,1,0);
      //   //Real dummy;
      //   //for (Integer i=0;i<Integer(model_block->get_bundle().size());i++){
      //   //model_block->get_bundle()[unsigned(i)].get_minorant(dummy,testdualinf,0,modelx(i),true);;
      //   //}
      // }
      // if (paramsp->QPget_use_socqp()){
      //   socqp.test_sysviol(testdualinf,solx);
      // }
      // else {
      //   QPadd_Qx(x,testdualinf);
      //   QPadd_Qx(solx,testdualinf);
      // }
      // Integer minind,maxind;
      // std::cout<<" norm2(Q*(x+dx)+c+transpose(A)*(y+dy)-zlb-dzlb+zub+dzub)="<<norm2(testdualinf)<<" ["<<min(testdualinf,&minind)<<"("<<minind<<"),"<<max(testdualinf,&maxind)<<"("<<maxind<<")]"<<std::endl; 
      // if ((model_block)&&(model_block->dim_model()>0)){
      //   std::cout<<" norm2(modelviol)="<<norm2(model_block->get_sysviol_model(solx))<<std::endl;
      //   std::cout<<" norm2(constrviol)="<<norm2(model_block->get_sysviol_constraints())<<std::endl;
      // }
      // std::cout.precision(10);


      // //  std::cout<<"y="<<transpose(y)<<std::endl;
      // //  std::cout<<"soly="<<transpose(soly)<<std::endl;
      // //  std::cout<<"izub="<<transpose(QPget_rhslbind());
      // //  std::cout<<"rzub="<<transpose(shortrhszub);
      // //  std::cout<<"drzub="<<transpose(dshortrhszub);
      // //  std::cout<<"rzub+drzub="<<transpose(shortrhszub+dshortrhszub);
      // //  std::cout<<"izlb="<<transpose(QPget_rhsubind());
      // //  std::cout<<"rzlb="<<transpose(shortrhszlb);
      // //  std::cout<<"drzlb="<<transpose(dshortrhszlb);
      // //  std::cout<<"rzlb+drzlb="<<transpose(shortrhszlb+dshortrhszlb);

      // testyinf.init(y.rowdim(),1,0.);
      // testprimalinf.init(y);
      // for(Integer i=0;i<QPget_rhslbind().dim();i++){
      //   Integer ind=QPget_rhslbind()(i);
      //   testyinf(ind)=y(ind)+soly(ind);
      // }
      // for(Integer i=0;i<QPget_rhsubind().dim();i++){
      //   Integer ind=QPget_rhsubind()(i);
      //   testyinf(ind)=y(ind)+soly(ind);
      // }
      // for(Integer i=0;i<shortrhszub.dim();i++){
      //   testyinf(QPget_rhslbind()(i))+=shortrhszub(i)+dshortrhszub(i);
      //   testprimalinf(QPget_rhslbind()(i))+=(next_mu-old_drsub(i)*old_drzub(i))/rhsslackub(i);
      // }
      // for(Integer i=0;i<shortrhszlb.dim();i++){
      //   testyinf(QPget_rhsubind()(i))-=(shortrhszlb(i)+dshortrhszlb(i));
      //   testprimalinf(QPget_rhsubind()(i))+=(-next_mu+old_drslb(i)*old_drzlb(i))/rhsslacklb(i);
      // }

      // //  std::cout<<"testyinf="<<transpose(testyinf)<<std::endl;
      // //  std::cout<<"dsrhs="<<transpose(dsrhs)<<std::endl;
      // //  std::cout<<"tpi="<<transpose(testprimalinf)<<std::endl;

      // std::cout<<"dsrhsinf="<<norm2(KKTdiagy%(dsrhs-testprimalinf))<<std::endl;
      // std::cout<<"ds-eq="<<norm2(ds+KKTdiagy%(soly+testprimalinf))<<std::endl;
      // testdualinf.init(testprimalinf,-1.);
      // testdualinf%=KKTdiagy;
      // QPadd_Ax(x,testdualinf);
      // testdualinf+=s;
      // std::cout<<"drhsinf="<<norm2(testdualinf+primalviol)<<std::endl;
      // QPadd_Ax(solx,testdualinf);
      // testprimalinf.init(soly,-1.);
      // testprimalinf%=KKTdiagy;
      // testdualinf+=testprimalinf;
      // std::cout<<"primalinf="<<norm2(testdualinf)<<std::endl;

      // std::cout<<" norm2(y+dy-rzlb-drzlb+rzub+drzub)="<<norm2(testyinf)<<std::endl; 
      // testprimalinf.init(s);
      // testprimalinf+=ds;
      // QPadd_Ax(x,testprimalinf);
      // QPadd_Ax(solx,testprimalinf);
      // std::cout<<" norm2(A*(x+dx)+s+ds)="<<norm2(testprimalinf)<<std::endl; 



    } //endif use_predcorr

    //--- linesearch and step to next point
    CH_Tools::Microseconds time_linesearch = clock.time();

    if (use_nbh) {
      QPselect_localstep_mu(next_mu, alpha, centering);
    } else {
      QPlinesearch(alpha);
      QPselect_mu(next_mu, alpha);
    }
    QPlinesearch_time += clock.time() - time_linesearch;

    //--- step to the next point  
    x.xpeya(solx, alpha);
    y.xpeya(soly, alpha);
    slackub.init(QPget_ub()(QPget_ubind()));
    slackub -= x(QPget_ubind());
    slacklb.init(x(QPget_lbind()));
    slacklb -= QPget_lb()(QPget_lbind());
    shortzlb.xpeya(dshortzlb, alpha);
    shortzub.xpeya(dshortzub, alpha);
    s.xpeya(ds, alpha);
    rhsslacklb.xpeya(drhsslacklb, alpha);
    rhsslackub.xpeya(drhsslackub, alpha);
    shortrhszlb.xpeya(dshortrhszlb, alpha);
    shortrhszub.xpeya(dshortrhszub, alpha);

    if (model_block) {
      if (model_block->do_step(alpha, x)) {
        if (cb_out()) {
          get_out() << "**** ERROR in QPSolverBasicStructures::QPpredcorr_step(...): model_block->do_step(.) failed" << std::endl;
        }
      }
    }

    if (paramsp->QPget_use_socqp()) {
      if (socqp.do_step(alpha, x)) {
        if (cb_out()) {
          get_out() << "**** ERROR in QPSolverBasicStructures::QPpredcorr_step(...): socqp.do_step(..) failed" << std::endl;
        }
      }
    }

    last_alpha = alpha;

    return status;
  }


  // *************************************************************************
  //                            iterate
  // *************************************************************************

  // loop till convergence to optimal solution

  int QPSolverBasicStructures::QPiterate() {
    Real objective_gap_eps = paramsp->QPget_objective_gap_eps();
    objective_gap_eps = min(objective_gap_eps, paramsp->QPget_min_objective_relprec());
    Real lower_bound_gap_eps = paramsp->QPget_lower_bound_gap_eps();
    Real upper_bound_gap_eps = paramsp->QPget_upper_bound_gap_eps();
    Real primeps = paramsp->QPget_primal_infeasibility_eps();
    Real dualeps = paramsp->QPget_dual_infeasibility_eps();
    Real ub = paramsp->QPget_upper_bound();
    Real lb = paramsp->QPget_lower_bound();
    Integer maxiter = paramsp->QPget_maxiter();


    // display parameters
    if (cb_out(1)) {
      get_out() << "    mi=" << maxiter;
      get_out() << " obge="; get_out().precision(5); get_out() << objective_gap_eps;
      get_out() << " lbge="; get_out().precision(5); get_out() << lower_bound_gap_eps;
      get_out() << " ubge="; get_out().precision(5); get_out() << upper_bound_gap_eps;
      get_out() << " pvie="; get_out().precision(5); get_out() << primeps;
      get_out() << " dvie="; get_out().precision(5); get_out() << dualeps;
      get_out() << " lb="; get_out().precision(10); get_out().width(11); get_out() << lb;
      get_out() << " ub="; get_out().precision(10); get_out().width(11); get_out() << ub; get_out().precision(2);
      get_out() << " xdim=" << QPget_xdim();
      get_out() << std::endl;
    }

    // compute primal and dual objective value, initialize "short" versions
    int status = 0;
    if ((status = QPcompute_values(true))) {
      if (cb_out(1)) {
        get_out() << "*** ERROR: QPSolverBasicStructures::QPiterate: QPcompute_values() returned " << status << std::endl;
      }
      return status;
    }

    Real alpha = 1.;
    Real prec = 1e-3;
    iter = 0;

    //central_path.push_back(QPCentralPathPoint(x,y,s,zlb,zub,rhszlb,rhszub,mu));

    //output
    if (cb_out(1)) {
      get_out().precision(2);
      get_out() << "    "; get_out().width(2); get_out() << iter;
      //get_out()<<"("<<central_path.size()-1<<")";
      get_out() << ": gappd="; get_out().width(7); get_out() << primalval - dualval;
      get_out() << " pv="; get_out().precision(8); get_out().width(10); get_out() << primalval;
      get_out() << " dv="; get_out().precision(8); get_out().width(10); get_out() << dualval;
      get_out().precision(2);
      get_out() << " epsob="; get_out().width(7); get_out() << objective_gap_eps * (fabs(primalval) + 1.);
      get_out() << " epspl="; get_out().width(7); get_out() << lower_bound_gap_eps * (primalval - lb);
      get_out() << " epsdu="; get_out().width(7); get_out() << upper_bound_gap_eps * (ub - dualval);
      get_out() << " pviol="; get_out().width(10); get_out() << n2primalviol;
      get_out() << " dviol="; get_out().width(10); get_out() << n2dualviol;
      get_out() << " yviol="; get_out().width(10); get_out() << n2yviol;
      if (model_block) {
        get_out() << " mpviol="; get_out().width(10); get_out() << n2modelprimalviol;
        get_out() << " mdviol="; get_out().width(10); get_out() << n2modeldualviol;
      }
      get_out() << " mu=";
      get_out().width(7); get_out() << mu;
      get_out() << std::endl;
    }


    //start with at  least one centering step
    do {

      iter++;

      old_x = x;
      old_y = y;
      old_s = s;
      old_zlb = zlb;
      old_zub = zub;
      old_rhszlb = rhszlb;
      old_rhszub = rhszub;
      old_mu = mu;

      //prec=min(1e-6,min(0.01*max(n2dualviol,n2primalviol),min(0.5*prec,1e-3*mu)));
      prec = min(1e-6, 1e-3 * mu);

      last_mu = mu;
      status = QPpredcorr_step(alpha, prec, mu,
        false, paramsp->QPget_use_neighborhood(), true);
      //write the solution back to the long vectors (before calling compute_values!)
      for (Integer i = 0; i < shortzub.dim(); i++) {
        zub(QPget_ubind()(i)) = shortzub(i);
      }
      for (Integer i = 0; i < shortzlb.dim(); i++) {
        zlb(QPget_lbind()(i)) = shortzlb(i);
      }
      for (Integer i = 0; i < shortrhszub.dim(); i++) {
        rhszub(QPget_rhslbind()(i)) = shortrhszub(i);
      }
      for (Integer i = 0; i < shortrhszlb.dim(); i++) {
        rhszlb(QPget_rhsubind()(i)) = shortrhszlb(i);
      }

      if (status) {
        if (cb_out(1)) {
          get_out() << "*** WARNING: QPSolverBasicStructures::QPiterate: QPpredcorr_step() returned " << status << std::endl;
        }
        return status;
      }

      //--- compute new values

      if ((status = QPcompute_values(false))) {
        if (cb_out(1)) {
          get_out() << "*** ERROR: QPSolverBasicStructures::QPiterate: QPcompute_values() returned " << status << std::endl;
        }
        return status;
      }

      if ((n2dualviol < 1e-8 * paramsp->QPget_KKTsolver()->QPget_blockH_norm()) && (dualval - ub > 1e-10 * (ub - lb)) && (cb_out())) {
        get_out().precision(12);
        get_out() << "**** WARNING QPSolverBasicStructures::QPiterate: dualval =" << dualval << " > " << ub << " = upper bound " << std::endl;
        Matrix tvec(x.rowdim(), 1, 0.);
        QPadd_Qx(x, tvec);
        get_out() << " quad_cost=" << ip(tvec, x) / 2.;
        get_out() << " lin_cost=" << ip(QPget_c(), x);
        get_out() << " offset=" << QPget_gamma();
        tvec.init(x.rowdim(), 1, 0.);
        Real tval = 0.;
        model_block->add_Bt_modelx(tval, tvec);
        get_out() << " dual_modcost=" << tval;
        get_out() << " dual_modval=" << tval + ip(tvec, x);
        if (model_block) {
          model_block->display_model_values(x, get_out());
        }
      }

      // // TEST begin
      // Matrix Qx(QPget_xdim(),1,0.);
      // QPadd_Qx(x,Qx);
      // std::cout<<" xQx/2="<<.5*ip(Qx,x);
      // std::cout<<" c'x="<<ip(QPget_c(),x);
      // std::cout<<" gamma="<<QPget_gamma();
      // if (model_block){
      //   std::cout<<" mxcost="<<model_block->globalx_cost(x);
      //   std::cout<<" mccost="<<model_block->constraints_cost();
      //   Matrix tmpvec(x.rowdim(),1,0.);
      //   Real tmpval=0.;
      //   model_block->add_Bt_modelx(tmpval,tmpvec);
      //   std::cout<<" mdcost="<<tmpval;
      // }
      // std::cout<<std::endl;
      // // TEST end

      //central_path.push_back(QPCentralPathPoint(x,y,s,zlb,zub,rhszlb,rhszub,mu));

      //output
      if (cb_out(1)) {
        get_out().precision(2);
        get_out() << "   c"; get_out().width(2); get_out() << iter;
        //get_out()<<"("<<central_path.size()-1<<")";
        get_out() << ": gappd="; get_out().width(7); get_out() << primalval - dualval;
        get_out() << " pv="; get_out().precision(8); get_out().width(10); get_out() << primalval;
        get_out() << " dv="; get_out().precision(8); get_out().width(10); get_out() << dualval;
        get_out().precision(2);
        get_out() << " epsob="; get_out().width(7); get_out() << objective_gap_eps * (fabs(primalval) + 1.);
        get_out() << " epspl="; get_out().width(7); get_out() << lower_bound_gap_eps * (primalval - lb);
        get_out() << " epsdu="; get_out().width(7); get_out() << upper_bound_gap_eps * (ub - dualval);
        get_out() << " pviol="; get_out().width(10); get_out() << n2primalviol / paramsp->QPget_KKTsolver()->QPget_blockA_norm();
        get_out() << " dviol="; get_out().width(10); get_out() << n2dualviol / paramsp->QPget_KKTsolver()->QPget_blockH_norm();
        get_out() << " yviol="; get_out().width(10); get_out() << n2yviol;
        if (model_block) {
          get_out() << " mpviol="; get_out().width(10); get_out() << n2modelprimalviol;
          get_out() << " mdviol="; get_out().width(10); get_out() << n2modeldualviol;
        }
        get_out() << " mu=";
        get_out().width(7); get_out() << mu;
        get_out() << " alpha="; get_out().width(4); get_out() << alpha;
        get_out() << std::endl;
        if (cb_out(2)) {
          get_out() << "    An=" << paramsp->QPget_KKTsolver()->QPget_blockA_norm();
          get_out() << " Hn=" << paramsp->QPget_KKTsolver()->QPget_blockH_norm();
          get_out() << std::endl;
          get_out() << "    tco=" << QPcoeff_time;
          get_out() << " tso=" << QPsolve_time;
          get_out() << " tiK=" << QPinitKKT_time;
          get_out() << " tpred=" << QPpredictor_time;
          get_out() << " tcorr=" << QPcorrector_time;
          get_out() << " tls=" << QPlinesearch_time;
          get_out() << " tpp=" << QPprecprep_time;
          get_out() << " tpm=" << QPprecprepmodel_time;
          get_out() << " tps=" << QPprecsolve_time;
          get_out() << " tmm=" << QPmatmult_time;
          get_out() << std::endl;
        }

        // get_out()<<"    iter="<<iter;
        // if((maxiter>0)&&(iter>=maxiter))
        // 	get_out()<<" T";
        // else
        // 	get_out()<<" F";
        // get_out()<<" pv-lb=";
        // if (primalval-lb <=1e-12*(fabs(primalval)+1.))
        // 	get_out()<<"T";
        // else
        // get_out()<<"F";
        // get_out()<<" ub-dv=";
        // if(ub-dualval<=1e-12*(fabs(dualval)+1.))
        // 	get_out()<<"T";
        // else
        // 	get_out()<<"F";
        // get_out()<<" alpha=";
        // if (alpha<=1e-8)
        // 	get_out()<<"T";
        // else
        // 	get_out()<<"F";
        // get_out()<<" pv-dv=";
        // if (primalval-dualval<=min(objective_gap_eps*(fabs(primalval)+1.),min(lower_bound_gap_eps*(primalval-lb),upper_bound_gap_eps*(ub-dualval))))
        // 	get_out()<<"T";
        // else
        // 	get_out()<<"F";
        // get_out()<<" pviol=";
        // if ((n2primalviol<=primeps)&&(n2modeldualviol<=primeps))
        // 	get_out()<<"T";
        // else
        // 	get_out()<<"F";
        // get_out()<<" dviol=";
        // if ((n2dualviol<=dualeps)&&(n2modelprimalviol<=primeps))
        // 	get_out()<<"T";
        // else
        // 	get_out()<<"F";
        // get_out()<<std::endl;
      }

    } while (
      ((maxiter < 0) || (iter < maxiter))
      &&
      (paramsp->QPget_use_neighborhood())
      &&
      (
        //(next_theta<0.)
        //||
        //(next_theta>paramsp->QPget_nbh_ub())
        //||
        (n2primalviol > primeps * paramsp->QPget_KKTsolver()->QPget_blockA_norm())
        ||
        (n2modeldualviol > primeps)
        ||
        (n2dualviol > dualeps * paramsp->QPget_KKTsolver()->QPget_blockH_norm())
        ||
        (n2yviol > dualeps)
        ||
        (n2modelprimalviol > dualeps)
        )
      );

    if ((last_alpha > .99) && (sigma > .95))
      mu = .9 * mu;
    sigma = .9; //in order to avoid initial misjudgement in the precision selection

    while (
      ((maxiter < 0) || (iter < maxiter))
      &&
      ((alpha > 1e-8) || (std::fabs(mu - last_mu) > 1e-8 * mu))
      &&
      (
        (//not sufficiently feasible
          (!paramsp->QPget_use_neighborhood())
          && (
            (n2primalviol > primeps * paramsp->QPget_KKTsolver()->QPget_blockA_norm())
            ||
            (n2modeldualviol > primeps)
            ||
            (n2dualviol > dualeps * paramsp->QPget_KKTsolver()->QPget_blockH_norm())
            ||
            (n2yviol > dualeps)
            ||
            (n2modelprimalviol > dualeps)
            )
          )
        ||
        (//not optimal
          (primalval - lb > 1e-12 * (fabs(primalval) + 1.))
          &&
          (ub - dualval > 1e-12 * (fabs(dualval) + 1.))
          &&
          (primalval - dualval > min(objective_gap_eps * (fabs(primalval) + 1.), min(lower_bound_gap_eps * (primalval - lb), upper_bound_gap_eps * (ub - dualval))))
          )
        // (
        //  ((dualval<lb)&&(primalval-dualval>(100*eps_Real)*(fabs(primalval)+1.)))
        //   ||
        //  (primalval>0.9*ub+0.1*dualval)
        //  //(primalval>0.1*ub+0.9*dualval)
        // )
        //(primalval>0.5*ub+0.5*dualval)
        )
      ) {

      iter++;

      old_x = x;
      old_y = y;
      old_s = s;
      old_zlb = zlb;
      old_zub = zub;
      old_rhszlb = rhszlb;
      old_rhszub = rhszub;
      old_mu = mu;

      //--- find the step direction

      // //Real old_gap=primalval-dualval;
      // Real old_pviol=max(n2primalviol,eps_Real);
      // Real old_dviol=max(n2dualviol,eps_Real);
      // //Real old_yviol=n2yviol;
      // Real old_mpviol=max(n2modelprimalviol,eps_Real);
      // Real old_mdviol=max(n2modeldualviol,eps_Real);

      // if (primalval-dualval<=min(objective_gap_eps*(fabs(primalval)+1.),min(lower_bound_gap_eps*(primalval-lb),upper_bound_gap_eps*(ub-dualval)))){
      //   //objective precision is reached, need to increase feasibility precision
      //   /*
      //   if ((n2primalviol>primeps)||(n2modeldualviol>primeps))
      // 	prec=min(prec,0.1*primeps);
      //   if ((n2dualviol>dualeps)||(n2yviol>dualeps)||(n2modelprimalviol>dualeps))
      // 	prec=min(prec,0.1*dualeps);
      //   */
      //   if ((n2primalviol>primeps*paramsp->QPget_KKTsolver()->QPget_blockA_norm())||
      // 	  (n2modeldualviol>primeps)||
      // 	  (n2dualviol>dualeps*paramsp->QPget_KKTsolver()->QPget_blockH_norm())||
      // 	  (n2yviol>dualeps)||
      // 	  (n2modelprimalviol>dualeps)){
      // 	fixed_mu=true;
      //   }
      // }

      if ((sigma > .95) || (alpha < 1.e-2)) {
        //prec=max(eps_Real*mu,min(1e-6,min(0.01*max(n2dualviol,n2primalviol),min(0.5*prec,1e-3*mu))));
        prec = max(eps_Real * mu, min(1e-6, min(max(1e-8, 0.01 * max(n2dualviol, n2primalviol)), min(0.9 * prec, 1e-3 * mu))));
      } else {
        prec = min(1e-2 * mu, 1e-6);
      }

      bool use_predcorr = (next_theta <= paramsp->QPget_nbh_lb()) && (paramsp->QPget_use_predictor_corrector());

      last_mu = mu;
      status = QPpredcorr_step(alpha,
        prec, mu,
        use_predcorr,
        paramsp->QPget_use_neighborhood(),
        false);
      //write the solution back to the long vectors (before calling compute_values!)
      for (Integer i = 0; i < shortzub.dim(); i++) {
        zub(QPget_ubind()(i)) = shortzub(i);
      }
      for (Integer i = 0; i < shortzlb.dim(); i++) {
        zlb(QPget_lbind()(i)) = shortzlb(i);
      }
      for (Integer i = 0; i < shortrhszub.dim(); i++) {
        rhszub(QPget_rhslbind()(i)) = shortrhszub(i);
      }
      for (Integer i = 0; i < shortrhszlb.dim(); i++) {
        rhszlb(QPget_rhsubind()(i)) = shortrhszlb(i);
      }
      if (status) {
        if (cb_out(1)) {
          get_out() << "*** WARNING: QPSolverBasicStructures::QPiterate: QPpredcorr_step() returned " << status << std::endl;
        }
        return status;
      }

      //--- compute new values
      if ((status = QPcompute_values(false))) {
        if (cb_out(1)) {
          get_out() << "*** ERROR: QPSolverBasicStructures::QPiterate: QPcompute_values() returned " << status << std::endl;
        }
        return status;
      }

      if ((n2dualviol < 1e-8 * paramsp->QPget_KKTsolver()->QPget_blockH_norm()) && (dualval - ub > 1e-10 * (ub - lb)) && (cb_out())) {
        get_out().precision(12);
        get_out() << "**** WARNING QPSolverBasicStructures::QPiterate: dualval =" << dualval << " > " << ub << " = upper bound " << std::endl;
        Matrix tvec(x.rowdim(), 1, 0.);
        QPadd_Qx(x, tvec);
        get_out() << " quad_cost=" << ip(tvec, x) / 2.;
        get_out() << " lin_cost=" << ip(QPget_c(), x);
        get_out() << " offset=" << QPget_gamma();
        tvec.init(x.rowdim(), 1, 0.);
        Real tval = 0.;
        model_block->add_Bt_modelx(tval, tvec);
        get_out() << " dual_modcost=" << tval;
        get_out() << " dual_modval=" << tval + ip(tvec, x);
        if (model_block) {
          model_block->display_model_values(x, get_out());
        }
      }

      // // TEST begin
      // Matrix Qx(QPget_xdim(),1,0.);
      // QPadd_Qx(x,Qx);
      // std::cout<<" xQx/2="<<.5*ip(Qx,x);
      // std::cout<<" c'x="<<ip(QPget_c(),x);
      // std::cout<<" gamma="<<QPget_gamma();
      // if (model_block){
      //   std::cout<<" mxcost="<<model_block->globalx_cost(x);
      //   std::cout<<" mccost="<<model_block->constraints_cost();
      //   Matrix tmpvec(x.rowdim(),1,0.);
      //   Real tmpval=0.;
      //   model_block->add_Bt_modelx(tmpval,tmpvec);
      //   std::cout<<" mdcost="<<tmpval;
      // }
      // std::cout<<std::endl;
      // // TEST end

      // if (
      // 	(primalval-dualval<=min(objective_gap_eps*(fabs(primalval)+1.),min(lower_bound_gap_eps*(primalval-lb),upper_bound_gap_eps*(ub-dualval))))||
      // 	((alpha<0.01)&&(!paramsp->QPget_use_neighborhood()))||
      // 	((fixed_mu)&&(alpha<0.99))||
      // 	((!full_step)&&(!paramsp->QPget_use_neighborhood())&&
      // 	 (
      // 	  ((n2primalviol>primeps*paramsp->QPget_KKTsolver()->QPget_blockA_norm())&&(n2primalviol>old_pviol*max(0.9,1.-alpha)))||
      // 	  ((n2dualviol>dualeps*paramsp->QPget_KKTsolver()->QPget_blockH_norm())&&(n2dualviol>old_dviol*max(0.9,1.-alpha)))||
      // 	  ((n2modelprimalviol>dualeps)&&(n2modelprimalviol>old_mpviol*max(0.9,1.-alpha)))||
      // 	  ((n2modeldualviol>dualeps)&&(n2modeldualviol>old_mdviol*max(0.9,1.-alpha)))
      // 	  //||((n2yviol>dualeps)&&(n2yviol>old_yviol*max(0.9,1.-alpha)))
      // 	  //||(primalval-dualval>old_gap*max(0.9,1.-alpha))
      // 	  )
      // 	 )){
      //   /*
      //   if (fixed_mu&&(alpha<0.1)){
      // 	mu*=10.;
      // 	std::cout<<" mu increased to "<<mu<<std::endl;
      //   }
      //   */
      //   if (cb_out(3)){
      // 	get_out()<<" fixed_mu:";
      // 	if (primalval-dualval<=min(objective_gap_eps*(fabs(primalval)+1.),min(lower_bound_gap_eps*(primalval-lb),upper_bound_gap_eps*(ub-dualval))))
      // 	  get_out()<<" p-d";
      // 	if ((alpha<0.01)&&(!paramsp->QPget_use_neighborhood()))
      // 	  get_out()<<" alpha="<<alpha;
      // 	if ((fixed_mu)&&(alpha<0.99))
      // 	  get_out()<<" fmu&al<.99 ";
      // 	if ((!full_step)&&(!paramsp->QPget_use_neighborhood())&&(n2primalviol>primeps*paramsp->QPget_KKTsolver()->QPget_blockA_norm())&&(n2primalviol>old_pviol*max(0.9,1.-alpha)))
      // 	  get_out()<<" pviol ";
      // 	if ((!full_step)&&(!paramsp->QPget_use_neighborhood())&&(n2dualviol>dualeps*paramsp->QPget_KKTsolver()->QPget_blockH_norm())&&(n2dualviol>old_dviol*max(0.9,1.-alpha)))
      // 	  get_out()<<" dviol ";
      // 	if ((!full_step)&&(!paramsp->QPget_use_neighborhood())&&(n2modelprimalviol>dualeps)&&(n2modelprimalviol>old_mpviol*max(0.9,1.-alpha)))
      // 	  get_out()<<" mpviol ";
      // 	if ((!full_step)&&(!paramsp->QPget_use_neighborhood())&&(n2modeldualviol>dualeps)&&(n2modeldualviol>old_mdviol*max(0.9,1.-alpha)))
      // 	  get_out()<<" mdviol ";
      // 	get_out()<<std::endl;
      //   }
      //   fixed_mu=true;	    
      // }
      // else {
      //   fixed_mu=false;
      // }


      //central_path.push_back(QPCentralPathPoint(x,y,s,zlb,zub,rhszlb,rhszub,mu));

      //output
      if (cb_out(1)) {
        get_out().precision(2);
        get_out() << "    "; get_out().width(2); get_out() << iter;
        //get_out()<<"("<<central_path.size()-1<<")";
        get_out() << ": gappd="; get_out().width(7); get_out() << primalval - dualval;
        get_out() << " pv="; get_out().precision(8); get_out().width(10); get_out() << primalval;
        get_out() << " dv="; get_out().precision(8); get_out().width(10); get_out() << dualval;
        get_out().precision(2);
        get_out() << " epsob="; get_out().width(7); get_out() << objective_gap_eps * (fabs(primalval) + 1.);
        get_out() << " epspl="; get_out().width(7); get_out() << lower_bound_gap_eps * (primalval - lb);
        get_out() << " epsdu="; get_out().width(7); get_out() << upper_bound_gap_eps * (ub - dualval);
        get_out() << " pviol="; get_out().width(10); get_out() << n2primalviol / paramsp->QPget_KKTsolver()->QPget_blockA_norm();
        get_out() << " dviol="; get_out().width(10); get_out() << n2dualviol / paramsp->QPget_KKTsolver()->QPget_blockH_norm();
        get_out() << " yviol="; get_out().width(10); get_out() << n2yviol;
        if (model_block) {
          get_out() << " mpviol="; get_out().width(10); get_out() << n2modelprimalviol;
          get_out() << " mdviol="; get_out().width(10); get_out() << n2modeldualviol;
        }
        get_out() << " mu=";
        get_out().width(7); get_out() << mu;
        get_out() << " alpha="; get_out().width(4); get_out() << alpha;
        get_out() << std::endl;
        if (cb_out(2)) {
          get_out() << "    An=" << paramsp->QPget_KKTsolver()->QPget_blockA_norm();
          get_out() << " Hn=" << paramsp->QPget_KKTsolver()->QPget_blockH_norm();
          get_out() << std::endl;
          get_out() << "    tco=" << QPcoeff_time;
          get_out() << " tso=" << QPsolve_time;
          get_out() << " tiK=" << QPinitKKT_time;
          get_out() << " tpred=" << QPpredictor_time;
          get_out() << " tcorr=" << QPcorrector_time;
          get_out() << " tls=" << QPlinesearch_time;
          get_out() << " tpp=" << QPprecprep_time;
          get_out() << " tpm=" << QPprecprepmodel_time;
          get_out() << " tps=" << QPprecsolve_time;
          get_out() << " tmm=" << QPmatmult_time;
          get_out() << std::endl;
        }

        // get_out()<<"    iter="<<iter;
        // if((maxiter>0)&&(iter>=maxiter))
        // 	get_out()<<" T";
        // else
        // 	get_out()<<" F";
        // get_out()<<" pv-lb=";
        // if (primalval-lb <=1e-12*(fabs(primalval)+1.))
        // 	get_out()<<"T";
        // else
        // get_out()<<"F";
        // get_out()<<" ub-dv=";
        // if(ub-dualval<=1e-12*(fabs(dualval)+1.))
        // 	get_out()<<"T";
        // else
        // 	get_out()<<"F";
        // get_out()<<" alpha=";
        // if (alpha<=1e-8)
        // 	get_out()<<"T";
        // else
        // 	get_out()<<"F";
        // get_out()<<" pv-dv=";
        // if (primalval-dualval<=min(objective_gap_eps*(fabs(primalval)+1.),min(lower_bound_gap_eps*(primalval-lb),upper_bound_gap_eps*(ub-dualval))))
        // 	get_out()<<"T";
        // else
        // 	get_out()<<"F";
        // get_out()<<" pviol=";
        // if ((n2primalviol<=primeps)&&(n2modeldualviol<=primeps))
        // 	get_out()<<"T";
        // else
        // 	get_out()<<"F";
        // get_out()<<" dviol=";
        // if ((n2dualviol<=dualeps)&&(n2modelprimalviol<=primeps))
        // 	get_out()<<"T";
        // else
        // 	get_out()<<"F";
        // get_out()<<std::endl;
      }

    } //end while


    if (cb_out() && (
      (n2primalviol > primeps * paramsp->QPget_KKTsolver()->QPget_blockA_norm())
      ||
      (n2modeldualviol > primeps)
      ||
      (n2dualviol > dualeps * paramsp->QPget_KKTsolver()->QPget_blockH_norm())
      ||
      (n2yviol > dualeps)
      ||
      (n2modelprimalviol > dualeps)
      )) {
      get_out() << "**** WARNING: QPSolverBasicStructures::QPiterate(): numerical difficulties might be ahead, terminating with poor feasibility; primalviol=" << n2primalviol << " dualviol=" << n2dualviol << " modeldualviol=" << n2modeldualviol << " modelprimalviol=" << n2modelprimalviol << " yviol=" << n2yviol << std::endl;
    }

    if (cb_out(1)) {
      get_out() << "    term: iter=" << iter;
      if ((maxiter > 0) && (iter >= maxiter))
        get_out() << " T";
      else
        get_out() << " F";
      get_out() << " pv-lb=";
      if (primalval - lb <= 1e-12 * (fabs(primalval) + 1.))
        get_out() << "T";
      else
        get_out() << "F";
      get_out() << " ub-dv=";
      if (ub - dualval <= 1e-12 * (fabs(dualval) + 1.))
        get_out() << "T";
      else
        get_out() << "F";
      get_out() << " alpha=";
      if (alpha <= 1e-8)
        get_out() << "T";
      else
        get_out() << "F";
      get_out() << " pv-dv=";
      if (primalval - dualval <= min(objective_gap_eps * (fabs(primalval) + 1.), min(lower_bound_gap_eps * (primalval - lb), upper_bound_gap_eps * (ub - dualval))))
        // if(((dualval>=lb)||(primalval-dualval<=objective_gap_eps*(fabs(primalval)+1.)))
        //    &&
        //    (primalval<=0.9*ub+0.1*dualval)
        //    )
        get_out() << "T";
      else
        get_out() << "F";
      get_out() << " pviol=";
      if ((n2primalviol <= primeps * paramsp->QPget_KKTsolver()->QPget_blockA_norm()) && (n2modeldualviol <= primeps))
        get_out() << "T";
      else
        get_out() << "F";
      get_out() << " dviol=";
      if ((n2dualviol <= dualeps * paramsp->QPget_KKTsolver()->QPget_blockH_norm()) && (n2modelprimalviol <= primeps))
        get_out() << "T";
      else
        get_out() << "F";
      get_out() << std::endl;
    }

    if ((maxiter < 0) || (iter < maxiter)) {
      status = 0;
    } else {
      status = 1;
    }

    //write the solution back to the long vectors
    for (Integer i = 0; i < shortzub.dim(); i++) {
      zub(QPget_ubind()(i)) = shortzub(i);
    }
    for (Integer i = 0; i < shortzlb.dim(); i++) {
      zlb(QPget_lbind()(i)) = shortzlb(i);
    }
    for (Integer i = 0; i < shortrhszub.dim(); i++) {
      rhszub(QPget_rhslbind()(i)) = shortrhszub(i);
    }
    for (Integer i = 0; i < shortrhszlb.dim(); i++) {
      rhszlb(QPget_rhsubind()(i)) = shortrhszlb(i);
    }

    return  status;
  }

  // *************************************************************************
  //                             QPIsolve
  // *************************************************************************

  // call starting_point and loop till convergence to optimal solution

  int QPSolverBasicStructures::QPIsolve(bool reinit, Real /* skip_factor */) {
    assert(paramsp->QPget_KKTsolver());
    Integer xdim = QPget_xdim();
    Integer ydim = QPget_ydim();
    //Integer moddim=0;
    //Integer modcdim=0;
    //if (model_block){
    //  moddim=model_block->dim_model();
    //  modcdim=model_block->dim_constraints();
    //}

    //if ((reinit)||(mu<=0.)||(skip_factor<0)||(central_path.size()==0)) {
    reinit = true;
    mu = 100. * std::log(xdim);
    central_path.clear();


    /*
       }
       else {
       int ind=min(3,int(central_path.size()));
       //int ind=int(Real(central_path.size())*skip_factor);
       //if (ind<0)
       //  ind=0;
       x=central_path[unsigned(ind)].x;
       y=central_path[unsigned(ind)].y;
       s=central_path[unsigned(ind)].s;
       zlb=central_path[unsigned(ind)].zlb;
       zub=central_path[unsigned(ind)].zub;
       rhszlb=central_path[unsigned(ind)].rhszlb;
       rhszub=central_path[unsigned(ind)].rhszub;
       mu=max(central_path[unsigned(ind)].mu,1e-3);
       central_path.resize(unsigned(ind));
       }
    */
    old_x.init(0, 1, 0.);
    old_y.init(0, 1, 0.);
    old_s.init(0, 1, 0.);
    old_zlb.init(0, 1, 0.);
    old_zub.init(0, 1, 0.);
    old_rhszlb.init(0, 1, 0.);
    old_rhszub.init(0, 1, 0.);

    last_alpha = -1.;
    last_mu = -1.;
    next_mu = -1.;
    last_theta = -1.;
    next_theta = -1.;
    large_predictor_cnt = 0;


    //----  initialize primal and dual starting point if necessary
    bool initx = false;
    bool initzlb = false;
    bool initzub = false;
    bool inits = false;
    bool initrhszlb = false;
    bool initrhszub = false;
    if ((reinit) || (x.dim() != xdim)) {
      x.init(xdim, 1, 0.);
      initx = true;
    }
    if ((reinit) || (zlb.dim() != xdim)) {
      zlb.init(xdim, 1, 0.);
      initzlb = true;
    }
    if ((reinit) || (zub.dim() != xdim)) {
      zub.init(xdim, 1, 0.);
      initzub = true;
    }
    if ((reinit) || (s.dim() != ydim)) {
      s.init(-QPget_rhslb());
      inits = true;
    }
    if ((reinit) || (rhszlb.dim() != ydim)) {
      rhszlb.init(ydim, 1, 0.);
      initrhszlb = true;
    }
    if ((reinit) || (rhszub.dim() != ydim)) {
      rhszub.init(ydim, 1, 0.);
      initrhszub = true;
    }
    if ((reinit) || (y.dim() != ydim)) {
      y.init(ydim, 1, 0.);
    }

    //ensure strictly interior starting point x
    Real sqrtmu = std::sqrt(mu);
    Integer lbi = 0;
    Integer lbind;
    if (lbi < QPget_lbind().dim()) {
      lbind = QPget_lbind()[lbi++];
    } else {
      lbind = xdim;
    }
    Integer ubi = 0;
    Integer ubind;
    if (ubi < QPget_ubind().dim()) {
      ubind = QPget_ubind()[ubi++];
    } else {
      ubind = xdim;
    }
    while ((lbind < xdim) || (ubind < xdim)) {
      if (lbind == ubind) {
        Real lb = QPget_lb()[lbind];
        Real ub = QPget_ub()[ubind];
        assert(1e-8 * std::fabs(ub) < ub - lb);
        if ((initx) || (x[lbind] <= lb) || (x[lbind] >= ub))
          x[lbind] = (lb + ub) / 2.;
        if (initzlb) {
          zlb[lbind] = mu / (x[lbind] - lb);
        }
        if (initzub) {
          zub[ubind] = mu / (ub - x[ubind]);
        }
        if (lbi < QPget_lbind().dim()) {
          lbind = QPget_lbind()[lbi++];
        } else {
          lbind = xdim;
        }
        if (ubi < QPget_ubind().dim()) {
          ubind = QPget_ubind()[ubi++];
        } else {
          ubind = xdim;
        }
      } else if (lbind < ubind) {
        Real lb = QPget_lb()[lbind];
        if ((initx) || (x[lbind] <= lb))
          x[lbind] = lb + sqrtmu;
        //x[lbind]=lb+max(1.,0.1*std::fabs(lb));
        if (initzlb) {
          zlb[lbind] = mu / (x[lbind] - lb);
        }
        if (lbi < QPget_lbind().dim()) {
          lbind = QPget_lbind()[lbi++];
        } else {
          lbind = xdim;
        }
      } else {
        Real ub = QPget_ub()[ubind];
        if ((initx) || (x[ubind] >= ub))
          x[ubind] = ub - sqrtmu;
        //x[ubind]=ub-max(1.,0.1*std::fabs(ub));
        if (initzub) {
          zub[ubind] = mu / (ub - x[ubind]);
        }
        if (ubi < QPget_ubind().dim()) {
          ubind = QPget_ubind()[ubi++];
        } else {
          ubind = xdim;
        }
      }
    }


    //ensure strictly interior starting point s
    //minus the lower bound on the rhs is the upper bound on s
    //minus the upper bound on the rhs is the lower bound on s
    lbi = 0;
    if (lbi < QPget_rhsubind().dim()) {
      lbind = QPget_rhsubind()[lbi++];
    } else {
      lbind = ydim;
    }
    ubi = 0;
    if (ubi < QPget_rhslbind().dim()) {
      ubind = QPget_rhslbind()[ubi++];
    } else {
      ubind = ydim;
    }
    while ((lbind < ydim) || (ubind < ydim)) {
      if (lbind == ubind) {
        Real lb = -QPget_rhsub()[lbind];
        Real ub = -QPget_rhslb()[ubind];
        assert(1e-8 * std::fabs(ub) < ub - lb);
        if ((inits) || (s[lbind] <= lb) || (s[lbind] >= ub))
          s[lbind] = (lb + ub) / 2.;
        if (initrhszlb) {
          rhszlb[lbind] = mu / (s[lbind] - lb);
        }
        if (initrhszub) {
          rhszub[ubind] = mu / (ub - s[ubind]);
        }
        y[lbind] = rhszlb[lbind] - rhszub[lbind];
        if (lbi < QPget_rhsubind().dim()) {
          lbind = QPget_rhsubind()[lbi++];
        } else {
          lbind = ydim;
        }
        if (ubi < QPget_rhslbind().dim()) {
          ubind = QPget_rhslbind()[ubi++];
        } else {
          ubind = ydim;
        }
      } else if (lbind < ubind) {
        Real lb = -QPget_rhsub()[lbind];
        if ((inits) || (s[lbind] <= lb))
          s[lbind] = lb + sqrtmu;
        //s[lbind]=lb+max(1.,0.1*std::fabs(lb));
        if (initrhszlb) {
          rhszlb[lbind] = mu / (s[lbind] - lb);
        }
        y[lbind] = rhszlb[lbind];
        if (lbi < QPget_rhsubind().dim()) {
          lbind = QPget_rhsubind()[lbi++];
        } else {
          lbind = ydim;
        }
      } else {
        Real ub = -QPget_rhslb()[ubind];
        if ((inits) || (s[ubind] >= ub))
          s[ubind] = ub - sqrtmu;
        //s[ubind]=ub-max(1.,0.1*std::fabs(ub));
        if (initrhszub) {
          rhszub[ubind] = mu / (ub - s[ubind]);
        }
        y[ubind] = -rhszub[ubind];
        if (ubi < QPget_rhslbind().dim()) {
          ubind = QPget_rhslbind()[ubi++];
        } else {
          ubind = ydim;
        }
      }
    }

    //std::cout<<"x="<<x;

    if (model_block) {
      if (model_block->reset_starting_point(x, mu)) {
        if (cb_out()) {
          get_out() << "**** ERROR in QPSolverBasicStructures::QPIsolve(..): model_block->reset_starting_point() failed" << std::endl;
        }
      }
    }


    //start solving
    int status = QPiterate();

    return status;
  }


}

