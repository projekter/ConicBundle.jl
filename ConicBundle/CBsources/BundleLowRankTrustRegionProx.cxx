/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/BundleLowRankTrustRegionProx.cxx
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



#include <math.h>
#include <stdlib.h>
#include "BundleLowRankTrustRegionProx.hxx"


using namespace CH_Matrix_Classes;

namespace ConicBundle {

  // *****************************************************************************
  //                       BundleLowRankTrustRegionProx::clean
  // *****************************************************************************

  void BundleLowRankTrustRegionProx::clean() {
    if (lamH.dim() == 0)
      return;

    //compute the SVD
    Indexmatrix piv;
    scaled_vecH = vecH;
    lamH.sqrt();
    scaled_vecH.scale_cols(lamH);
    Symmatrix S;
    rankadd(scaled_vecH, S, 1., 0., 1);
    Matrix P;
    S.eig(P, lamH, false);
    Real lmax = max(1e-10, lamH(0));
    Integer cnt = 0;
    Integer ubcnt = min(max_columns, lamH.dim());
    while ((cnt < ubcnt) && (lamH(cnt) > 1e-6 * lmax))
      cnt++;
    lamH.reduce_length(cnt);
    P.delete_cols(Range(cnt, P.coldim() - 1));
    lamHi = lamH;
    lamHi.sqrt();
    lamHi.inv();
    P.scale_cols(lamHi);
    genmult(scaled_vecH, P, vecH);
    needs_cleaning = false;

    if (cb_out(2)) {
      get_out() << " BLRclean(" << lamH.dim();
      if (lamH.dim() > 0)
        get_out() << "," << lamH(0) << "," << lamH(lamH.dim() - 1);
      get_out() << ")" << std::endl;
    }
  }

  // *****************************************************************************
  //                       BundleLowRankTrustRegionProx::init
  // *****************************************************************************

  void BundleLowRankTrustRegionProx::init(const CH_Matrix_Classes::Matrix& in_vecH,
    const CH_Matrix_Classes::Matrix& in_lamH) {
    //assert column vector
    assert(lamH.coldim() == 1);
    //assert matching dimsions
    assert(vecH.coldim() == lamH.rowdim());
    //assert "orthogonality indicator"
    assert(norm2(transpose(vecH) * vecH - Diag(Matrix(vecH.coldim(), 1, 1.))) < 1e-6);
    assert(std::fabs(ip(vecH, vecH) - vecH.coldim()) < 1e-10 * max(1, vecH.coldim()));

    vecH = in_vecH;
    lamH = in_lamH;
    lamHi.init(Integer(0), Integer(1), 0.);
    sqrtlamHi.init(Integer(0), Integer(1), 0.);
    scaled_vecH.init(0, 0, 0.);
    LinvindHt.init(0, 0, 0.);
    old_fixed_ind.init(0, 0, Integer(0));
  }


  // *****************************************************************************
  //                       BundleLowRankTrustRegionProx::set_weightu
  // *****************************************************************************

  void BundleLowRankTrustRegionProx::set_weightu(CH_Matrix_Classes::Real in_weightu) {
    if (std::fabs(weightu - in_weightu) > 1e-10 * fabs(weightu)) {
      weightu = in_weightu;
      lamHi.init(Integer(0), Integer(1), 0.);
      sqrtlamHi.init(Integer(0), Integer(1), 0.);
      lamHi.init(0, 0, 0.);
      sqrtlamHi.init(0, 0, 0.);
      LinvindHt.init(0, 0, 0.);
      old_fixed_ind.init(0, 0, Integer(0));
      compute_corr();
    }
  }

  // *****************************************************************************
  //                       BundleLowRankTrustRegionProx::norm_sqr
  // *****************************************************************************

  Real BundleLowRankTrustRegionProx::norm_sqr(const Matrix& B) const {
    if (lamH.dim() == 0) {
      return weightu * ip(B, B);
    }
    Matrix tmpmat;
    genmult(vecH, B, tmpmat, 1., 0., 1);
    return weightu * ip(B, B) + normDsquared(tmpmat, lamH);
  }


  // *****************************************************************************
  //                       BundleLowRankTrustRegionProx::dnorm_sqr
  // *****************************************************************************

  Real BundleLowRankTrustRegionProx::dnorm_sqr(const MinorantPointer& B) const {
    if (lamH.dim() == 0) {
      return B.dual_norm_squared() / weightu;
    }
    if (lamHi.dim() == 0) {
      lamHi = lamH;
      lamHi += weightu;
      lamHi.inv();
      lamHi %= lamH;
    }
    Real dummy;
    Matrix tmpvec(vecH.rowdim(), 1); chk_set_init(tmpvec, 1);
    B.get_minorant(dummy, tmpvec, 0);
    Matrix tmpmat;
    genmult(vecH, tmpvec, tmpmat, 1., 0., 1);
    return (B.dual_norm_squared() - normDsquared(tmpmat, lamHi)) / weightu;
  }

  // *****************************************************************************
  //                                add_H
  // *****************************************************************************

  int BundleLowRankTrustRegionProx::add_H(Symmatrix& big_sym,
    Integer start_index) const {
    assert((start_index >= 0) && (big_sym.rowdim() - start_index >= vecH.rowdim()));
    for (Integer i = 0; i < vecH.rowdim(); i++)
      big_sym(i + start_index, i + start_index) += weightu;
    if (lamH.rowdim() > 0) {
      if (big_sym.rowdim() == vecH.rowdim()) {
        scaledrankadd(vecH, lamH, big_sym, 1., 1.);
      } else {
        Symmatrix tmpsym;
        scaledrankadd(vecH, lamH, tmpsym);
        for (Integer i = 0; i < tmpsym.rowdim(); i++) {
          for (Integer j = i; j < tmpsym.rowdim(); j++)
            big_sym(i + start_index, j + start_index) += tmpsym(i, j);
        }
      }
    }
    return 0;
  }

  // *****************************************************************************
  //                       BundleLowRankTrustRegionProx::add_Hx
  // *****************************************************************************

  Matrix& BundleLowRankTrustRegionProx::add_Hx(const Matrix& x,
    Matrix& outplusHx,
    Real alpha) const {
    outplusHx.xpeya(x, weightu * alpha);
    if (lamH.dim() == 0) {
      return outplusHx;
    }
    Matrix tmpmat;
    genmult(vecH, x, tmpmat, alpha, 0., 1);
    tmpmat.scale_rows(lamH);
    return genmult(vecH, tmpmat, outplusHx, 1., 1.);
  }

  // *****************************************************************************
  //                       BundleLowRankTrustRegionProx::get_precond
  // *****************************************************************************

  void BundleLowRankTrustRegionProx::get_precond(Matrix& inD,
    const Matrix*& Vp) const {
    inD.init(vecH.rowdim(), 1, weightu);
    if (lamH.rowdim() > 0) {
      if (scaled_vecH.dim() == 0) {
        Matrix tmpvec(lamH);
        tmpvec.sqrt();
        scaled_vecH.init(vecH);
        scaled_vecH.scale_cols(tmpvec);
      }
      Vp = &scaled_vecH;
    } else {
      Vp = 0;
    }
  }

  // *****************************************************************************
  //                       BundleLowRankTrustRegionProx::compute_QP_costs
  // *****************************************************************************


  int BundleLowRankTrustRegionProx::compute_QP_costs(Symmatrix& Q,
    Matrix& d,
    Real& offset,
    const MinorantPointer& constant_minorant,
    const MinorantBundle& bundle,
    const Matrix& y,
    const MinorantPointer& groundset_minorant,
    Indexmatrix* yfixed) {
    //--- determine the fixed indices and values
    Indexmatrix ind(y.dim(), Integer(1));
    ind.init(0, 0, 0.);
    Matrix val(y.dim(), Integer(1));
    val.init(0, 0, 0.);
    _y.newsize(y.dim(), Integer(1)); chk_set_init(_y, 1);
    Integer ydim = 0;
    if (yfixed) {
      for (Integer i = 0; i < yfixed->dim(); i++) {
        if ((*yfixed)(i) > 0) {
          (*yfixed)(i) = 1;
          ind.concat_below(i);
          val.concat_below(y(i));
        } else {
          _y(ydim++) = y(i);
        }
      }
      _y.reduce_length(ydim);
    } else {
      _y = y;
      ydim = y.dim();
    }
    assert(ydim == y.dim() - ind.dim());
    bool fixed_values = (ind.dim() > 0);
    bool fixed_changed = !equal(old_fixed_ind, ind);
    if (fixed_changed) {
      old_fixed_ind = ind;
    }

    //--- compute the factorization on the current subset of indices
    Matrix tmpmat;
    if (lamHi.dim() == 0) {
      lamHi = lamH;
      lamHi += weightu;
      lamHi.inv();
      lamHi %= lamH;
    }
    if (!fixed_values) {
      if (LinvindHt.dim() != vecH.rowdim()) {
        LinvindHt.init(vecH, 1. / sqrt(weightu), 1);
        tmpmat = lamHi;
        tmpmat.sqrt();
        LinvindHt.scale_rows(tmpmat);
      }
    } else if ((fixed_changed) ||
      ((lamH.rowdim() > 0) && (LinvindHt.coldim() + old_fixed_ind.dim() != vecH.rowdim()))
      ) {
      LinvindHt.init(vecH, 1. / sqrt(weightu), 1);
      LinvindHt.delete_cols(old_fixed_ind);
      Symmatrix tmpsym;
      rankadd(LinvindHt, tmpsym, weightu);
      for (Integer i = 0; i < lamH.rowdim(); i++)
        tmpsym(i, i) += 1. / lamH(i);
      if (tmpsym.Chol_factor()) {
        if (cb_out())
          get_out() << "*** ERROR in BundleLowRankTrustRegionProx::compute_QP_costs(...): tmpsym.Chol_factor failed" << std::endl;
        return 1;
      }
      tmpsym.Chol_Lsolve(LinvindHt);
    }

    //--- get the bundle data into matrix form
    Integer xdim = Integer(bundle.size());
    _A.newsize(ydim, xdim); chk_set_init(_A, 1);
    _b.newsize(ydim, Integer(1)); chk_set_init(_b, 1);
    _c.newsize(xdim, Integer(1)); chk_set_init(_c, 1);
    _delta = 0.;

    if (groundset_minorant.get_minorant(_delta, _b, 0, 1., false, fixed_values ? &ind : 0, fixed_values ? &val : 0)) {
      if (cb_out())
        get_out() << "*** ERROR in BundleDiagonalProx::compute_QP_costs(...): groundset_minorant.get_minorant failed" << std::endl;
      return 1;
    }
    if ((!constant_minorant.empty()) && (constant_minorant.get_minorant(_delta, _b, 0, 1., true, fixed_values ? &ind : 0, fixed_values ? &val : 0))) {
      if (cb_out())
        get_out() << "*** ERROR in BundleDiagonalProx::compute_QP_costs(...): groundset_minorant.get_minorant failed" << std::endl;
      return 1;
    }

    for (Integer i = 0; i < xdim; i++) {
      if (bundle[unsigned(i)].get_minorant(_c(i), _A, i, 1, false, fixed_values ? &ind : 0, fixed_values ? &val : 0)) {
        if (cb_out())
          get_out() << "*** ERROR in BundleDiagonalProx::compute_QP_costs(...): groundset_minorant.get_minorant failed" << std::endl;
        return 1;
      }
    }

    //compute the linear and part of the constant coefficient
    genmult(LinvindHt, _A, tmpmat);
    rankadd(tmpmat, Q, -1., 0., 1);
    old_LinvQ = Q;
    rankadd(_A, Q, 1. / weightu, 1., 1);

    genmult(LinvindHt, _b, oldd);
    offset = _delta + ip(_b, _y) - (ip(_b, _b) / weightu - ip(oldd, oldd)) / 2.;

    tmpmat.init(_y);
    tmpmat.xpeya(_b, -1. / weightu);
    genmult(LinvindHt, oldd, tmpmat, 1., 1., 1);

    d = _c;
    genmult(_A, tmpmat, d, 1., 1., 1);

    oldd = d;
    oldoffset = offset;


    return 0;
  }


  // *****************************************************************************
  //                       BundleLowRankTrustRegionProx::update_QP_costs
  // *****************************************************************************

  int BundleLowRankTrustRegionProx::update_QP_costs(Symmatrix& delta_Q,
    Matrix& delta_d,
    Real& delta_offset,
    const MinorantPointer& /* constant_minorant */,
    const MinorantBundle& /* bundle */,
    const Matrix& /* y */,
    const MinorantPointer& /* subg */,
    const MinorantPointer& delta_subg,
    const Indexmatrix& delta_index,
    Indexmatrix* yfixed) {

    //--- determine the changes in fixed indices
    Integer xdim = _A.coldim();
    Integer change_dim = delta_index.dim();
    Indexmatrix new_fixed_ind(change_dim, Integer(1));
    new_fixed_ind.init(0, 0, 0.);
    Indexmatrix new_fixed_newind(change_dim, Integer(1));
    new_fixed_newind.init(0, 0, 0.);

    _delta += delta_subg.offset();

    Integer corr_cnt = 0;
    for (Integer j = 0; j < delta_index.dim(); j++) {
      Integer ind = delta_index[j];
      if (yfixed) {
        switch ((*yfixed)(ind)) {
        case 0:
        {
          //was and stayed not fixed, but the subgradient changed
          //treat this later, if no case 2 occured
          Real sgval = delta_subg.coeff(ind);
          while ((corr_cnt < old_fixed_ind.dim()) && (old_fixed_ind(corr_cnt) < ind))
            corr_cnt++;
          ind -= corr_cnt;
          _b[ind] += sgval;
          break;
        }
        case 2:
        {
          new_fixed_ind.concat_below(ind);
          while ((corr_cnt < old_fixed_ind.dim()) && (old_fixed_ind(corr_cnt) < ind))
            corr_cnt++;
          ind -= corr_cnt;
          new_fixed_newind.concat_below(ind);
          (*yfixed)(ind) = 1;
          break;
        }
        default:
          if (cb_out())
            get_out() << "*** ERROR in BundleLowRankTrustRegionProx::update_QP_costs(...):  internal error, yfixed(" << ind << ")=" << (*yfixed)(ind) << " should not occur here" << std::endl;
          std::abort();
        }
      } else {
        Real sgval = delta_subg.coeff(ind);
        while ((corr_cnt < old_fixed_ind.dim()) && (old_fixed_ind(corr_cnt) < ind))
          corr_cnt++;
        ind -= corr_cnt;
        _b[ind] += sgval;
      }
    }

    //--- check whether Hind_chol changed and recompute all if so
    Matrix tmpmat;

    if (new_fixed_ind.dim() > 0) {
      //prepare LinvindHt
      old_fixed_ind.concat_below(new_fixed_ind);
      Indexmatrix sind;
      sortindex(old_fixed_ind, sind);
      old_fixed_ind = old_fixed_ind(sind);
      LinvindHt.init(vecH, 1. / weightu, 1);
      LinvindHt.delete_cols(old_fixed_ind);
      Symmatrix tmpsym;
      rankadd(LinvindHt, tmpsym, weightu);
      for (Integer i = 0; i < lamH.rowdim(); i++)
        tmpsym(i, i) += 1. / lamH(i);
      if (tmpsym.Chol_factor()) {
        if (cb_out())
          get_out() << "*** ERROR in BundleLowRankTrustRegionProx::update_QP_costs(...): tmpsym.Chol_factor failed" << std::endl;
        return 1;
      }
      tmpsym.Chol_Lsolve(LinvindHt);

      //update the values of _delta and _c and partly Q by the newly fixed values
      _delta += ip(_b(new_fixed_ind), _y(new_fixed_ind));
      tmpmat = _A.rows(new_fixed_newind);
      genmult(tmpmat, _y(new_fixed_newind), _c, 1., 1., 1);
      rankadd(tmpmat, delta_Q, -1 / weightu, 0., 1);

      //eliminate the newly fixed rows
      _A.delete_rows(new_fixed_ind);
      _b.delete_rows(new_fixed_ind);
      _y.delete_rows(new_fixed_ind);

      //recompute the low rank part of Q
      delta_Q -= old_LinvQ;
      genmult(LinvindHt, _A, tmpmat);
      rankadd(tmpmat, old_LinvQ, -1., 0., 1);
      delta_Q += old_LinvQ;


    }//endif (new_fixed_ind)
    else {
      delta_Q.init(xdim, 0.);
    }


    delta_offset = -oldoffset;
    delta_d = -oldd;

    genmult(LinvindHt, _b, oldd);
    oldoffset = _delta + ip(_b, _y) - (ip(_b, _b) / weightu - ip(oldd, oldd)) / 2.;
    delta_offset += oldoffset;

    tmpmat.init(_y);
    tmpmat.xpeya(_b, -1. / weightu);
    genmult(LinvindHt, oldd, tmpmat, 1., 1., 1);

    oldd = _c;
    genmult(_A, tmpmat, oldd, 1., 1., 1);
    delta_d += oldd;


    return 0;
  }


  // *****************************************************************************
  //                       BundleLowRankTrustRegionProx::apply_modification
  // *****************************************************************************

  int BundleLowRankTrustRegionProx::apply_modification(const GroundsetModification& gsmdf) {
    Integer olddim = vecH.rowdim();
    if (gsmdf.old_vardim() != olddim) {
      if (cb_out())
        get_out() << "**** ERROR BundleLowRankTrustRegionProx::apply_modification: dim=" << olddim << " but modification assumes " << gsmdf.old_vardim() << std::endl;
      return 1;
    }

    vecH.enlarge_below(gsmdf.appended_vardim(), 0.);

    if (gsmdf.map_to_old_variables()) {
      vecH = vecH.rows(*(gsmdf.map_to_old_variables()));
    }
    compute_corr();

    scaled_vecH.init(0, 0, 0.);
    lamHi.init(0, 0, 0.);
    sqrtlamHi.init(0, 0, 0.);
    LinvindHt.init(0, 0, 0.);
    old_fixed_ind.init(0, 0, Integer(0));

    return 0;
  }

  // *****************************************************************************
  //                       BundleLowRankTrustRegionProx::projected_clone
  // *****************************************************************************

        /** @brief in order to allow for fixed variables, this generates a clone restricted to the given indices
     */
  BundleProxObject* BundleLowRankTrustRegionProx::projected_clone(const CH_Matrix_Classes::Indexmatrix& indices) {
    Matrix tmpvecH(vecH);
    Matrix tmplamH(lamH);
    tmplamH.sqrt();
    tmpvecH.scale_cols(tmplamH);
    Matrix tmpmat = tmpvecH.rows(indices);
    Indexmatrix piv;
    Integer r = tmpmat.QR_factor(piv);
    tmpvecH = tmpmat.rows(Range(0, r));
    Symmatrix S;
    rankadd(tmpvecH, S, 1., 0., 1);
    Real maxval = 1.;
    for (Integer i = 0; i < S.rowdim(); i++)
      maxval = max(maxval, S(i, i));
    S /= maxval;
    if (S.eig(tmpvecH, tmplamH)) {
      if (cb_out())
        get_out() << "**** WARNING BundleLowRankTrustRegionProx::projected_clone(): eig failed" << std::endl;
    }
    tmplamH *= maxval;
    tmpvecH.enlarge_below(tmpmat.rowdim(), 0.);
    tmpmat.Q_times(tmpvecH, r);

    BundleLowRankTrustRegionProx* pp = new BundleLowRankTrustRegionProx(tmpvecH, tmplamH, 0, get_use_local_metric(), this, 0);
    pp->set_weightu(weightu);
    pp->apply_factor(factor);
    if (get_variable_metric_selection())
      pp->set_variable_metric_selection(get_variable_metric_selection()->clone_VariableMetricSelection());
    return pp;
  }



  // *****************************************************************************
  //         BundleLowRankTrustRegionProx::apply_variable_metric
  // *****************************************************************************

  int BundleLowRankTrustRegionProx::apply_variable_metric(VariableMetricModel* groundset,
    VariableMetricModel* model,
    const Matrix& /* aggr */,
    Integer y_id,
    const Matrix& y,
    bool descent_step,
    Real& current_weight,
    Real model_maxviol,
    const Indexmatrix* in_new_indices) {
    assert(aft_stack.size() == 0);

    if (in_new_indices) {
      if (cb_out())
        get_out() << "**** WARNING BundleLowRankTrustRegionProx::apply_variable_metric(): this implementation does not support updates on subsets of indices -> updating everything" << std::endl;
    }


    new_indices = in_new_indices;

    if (y.dim() == 0)
      return 0;
    aft_stack.clear();

    if ((weightu <= 0.) || (current_weight <= 0.)) {
      //not initialized, initialize
      if (current_weight > 0.) {
        weightu = current_weight;
      } else {
        if (weightu <= 0.)
          weightu = 1.;
        current_weight = weightu;
      }
      vecH.init(y.dim(), 0, 0.);
      lamH.init(0, 1, 0.);
    }

    if ((descent_step) || (vecH.rowdim() != y.dim())) {
      weightu = max(current_weight, 1e-10);
      vecH.init(y.dim(), 0, 0.);
      lamH.init(0, 1, 0.);
      new_indices = 0;
    }


    max_columns = min(max(y.dim() / 5, 30), y.dim());
    needs_cleaning = false;

    int err = 0;
    if (groundset->variable_metric_transform()->add_variable_metric(*this, y_id, y,
      descent_step,
      weightu,
      model_maxviol,
      in_new_indices)) {
      if (cb_out())
        get_out() << "**** WARNING BundleLowRankTrustRegionProx::apply_variable_metric(): groundset->add_variable_metric(...) failed " << std::endl;
      err++;
    }
    if ((model) && (model->variable_metric_transform()->add_variable_metric(*this, y_id, y,
      descent_step,
      weightu,
      model_maxviol,
      in_new_indices))) {
      if (cb_out())
        get_out() << "**** WARNING BundleLowRankTrustRegionProx::apply_variable_metric(): model->transform()->add_variable_metric(...) failed " << std::endl;
      err++;
    }

    if ((needs_cleaning) || (lamH.dim() > max_columns))
      clean();

    if (cb_out(2)) {
      get_out() << " BLRTRS(" << lamH.dim();
      if (lamH.dim() > 0)
        get_out() << "," << lamH(0) << "," << lamH(lamH.dim() - 1);
      get_out() << ")";
    }

    current_weight = weightu;

    assert(aft_stack.size() == 0);
    new_indices = 0;
    aft_stack.clear();

    compute_corr();

    scaled_vecH.init(0, 0, 0.);
    lamHi.init(0, 0, 0.);
    sqrtlamHi.init(0, 0, 0.);
    LinvindHt.init(0, 0, 0.);
    old_fixed_ind.init(0, 0, Integer(0));

    return err;
  }


  // *****************************************************************************
  //         BundleLowRankTrustRegionProx::add_lowrank_variable_metric
  // *****************************************************************************

  int BundleLowRankTrustRegionProx::add_variable_metric(Matrix& /* diagH */,
    Matrix& in_vecH) {
    if (in_vecH.coldim() == 0)
      return 0;

    assert((aft_stack.size() != 0) || (in_vecH.rowdim() == vecH.rowdim()));

    //apply the transposed trafos down the stack
    Real fun_factor = 1.;
    Matrix tmpmat;
    for (Integer i = Integer(aft_stack.size()); --i >= 0;) {
      fun_factor *= aft_stack[unsigned(i)]->get_fun_coeff();
      if (aft_stack[unsigned(i)]->get_arg_trafo()) {
        genmult(*aft_stack[unsigned(i)]->get_arg_trafo(), in_vecH, tmpmat, 1., 0., 1);
        swap(tmpmat, in_vecH);
      }
    }

    //add vecH*diag(lam_vecH)*transpose(vecH) to H
    needs_cleaning = true; //(vecH.coldim()>0);
    vecH.concat_right(in_vecH);
    lamH.enlarge_below(in_vecH.coldim(), fun_factor);

    if ((needs_cleaning) || (lamH.dim() > 2 * max_columns))
      clean();

    return 0;
  }


  // *****************************************************************************
  //         BundleLowRankTrustRegionProx::apply_Hinv
  // *****************************************************************************

  Matrix& BundleLowRankTrustRegionProx::apply_Hinv(Matrix& x) const {
    //first transform x if necessary
    Real fun_coeff = 1.;
    Matrix tmpmat;
    for (Integer i = Integer(aft_stack.size()); --i >= 0;) {
      fun_coeff *= aft_stack[unsigned(i)]->get_fun_coeff();
      if (aft_stack[unsigned(i)]->get_arg_trafo()) {
        genmult(*aft_stack[unsigned(i)]->get_arg_trafo(), x, tmpmat, 1., 0., 1);
        swap(tmpmat, x);
      }
    }

    //apply the inverse
    x /= weightu;
    if (lamH.dim() != 0) {
      if (lamHi.dim() == 0) {
        lamHi = lamH;
        lamHi += weightu;
        lamHi.inv();
        lamHi %= lamH;
      }
      Matrix tmpmat;
      genmult(vecH, x, tmpmat, 1., 0., 1);
      tmpmat.scale_rows(lamHi);
      genmult(vecH, tmpmat, x, -1., 1.);
    }

    //apply the trafos up the stack  
    for (unsigned i = 0; i < aft_stack.size(); i++) {
      fun_coeff *= aft_stack[unsigned(i)]->get_fun_coeff();
      if (aft_stack[unsigned(i)]->get_arg_trafo()) {
        genmult(*aft_stack[unsigned(i)]->get_arg_trafo(), x, tmpmat);
        swap(tmpmat, x);
      }
    }

    if (fun_coeff != 1.)
      x *= fun_coeff;

    return x;
  }

  // *****************************************************************************
  //                       BundleLowRankTrustRegionProx::mfile_data
  // *****************************************************************************


  int BundleLowRankTrustRegionProx::mfile_data(std::ostream& out) const {
    out << "clear weightu qp_vecHscale qp_lamHscale;\n";
    out << "weightu=";
    out.precision(16);
    out.width(18);
    out << weightu << "\n";
    out << "qp_vecHscale=[";
    for (Integer i = 0; i < vecH.rowdim(); i++) {
      for (Integer j = 0; j < vecH.coldim(); j++) {
        out << " ";
        out.precision(16);
        out.width(18);
        out << vecH(i, j);
      }
      if (i < vecH.rowdim() - 1)
        out << "\n";
    }
    out << "];\n";
    out << "qp_lamHscale=[";
    for (Integer i = 0; i < lamH.rowdim(); i++) {
      for (Integer j = 0; j < lamH.coldim(); j++) {
        out << " ";
        out.precision(16);
        out.width(18);
        out << lamH(i, j);
      }
      if (i < lamH.rowdim() - 1)
        out << "\n";
    }
    out << "];\n";
    return 0;
  }





}  //namespace ConicBundle
