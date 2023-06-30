/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/BundleDenseTrustRegionProx.cxx
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
#include "BundleDenseTrustRegionProx.hxx"
#include "mymath.hxx"
#include "sparssym.hxx"
#include "Groundset.hxx"
#include "BundleModel.hxx"


using namespace CH_Matrix_Classes;

namespace ConicBundle {

  // *****************************************************************************
  //                       BundleDenseTrustRegionProx::set_weightu
  // *****************************************************************************

  void BundleDenseTrustRegionProx::set_weightu(CH_Matrix_Classes::Real in_weightu) {
    if (fabs(weightu - in_weightu) > 1e-10 * fabs(weightu)) {
      is_factored = false;
      Hind_chol.init(0, 0.);
      weightu = in_weightu;
      old_fixed_ind.init(0, 0, Integer(0));
      compute_corr();
    }
  }

  // *****************************************************************************
  //                       BundleDenseTrustRegionProx::norm_sqr
  // *****************************************************************************

  Real BundleDenseTrustRegionProx::norm_sqr(const Matrix& B) const {
    return ip(B, H * B) + weightu * ip(B, B);
  }


  // *****************************************************************************
  //                       BundleDenseTrustRegionProx::dnorm_sqr
  // *****************************************************************************

  Real BundleDenseTrustRegionProx::dnorm_sqr(const MinorantPointer& B) const {
    if (!is_factored) {
      is_factored = true;
      Hchol = H;
      for (Integer i = 0; i < Hchol.rowdim(); i++)
        Hchol(i, i) += weightu;
      if (Hchol.Chol_factor()) {
        if (cb_out())
          get_out() << "ERROR in BundleDenseTrustRegionProx::update_eta_step(...): H.Chol_factor() failed" << std::endl;
        return 1;
      }
    }
    Real dummy;
    Matrix tmpmat(H.rowdim(), 1); chk_set_init(tmpmat, 1);
    B.get_minorant(dummy, tmpmat, 0);
    Hchol.Chol_Lsolve(tmpmat);
    return ip(tmpmat, tmpmat);
  }


  // *****************************************************************************
  //                                add_H
  // *****************************************************************************

  int BundleDenseTrustRegionProx::add_H(Symmatrix& big_sym,
    Integer start_index) const {
    assert((start_index >= 0) && (big_sym.rowdim() - start_index >= H.rowdim()));
    if (big_sym.rowdim() == H.rowdim()) {
      big_sym += H;
      for (Integer i = 0; i < H.rowdim(); i++)
        big_sym(i, i) += weightu;
    } else {
      for (Integer i = 0; i < H.rowdim(); i++) {
        big_sym(i + start_index, i + start_index) += weightu + H(i, i);
        for (Integer j = i + 1; j < H.rowdim(); j++)
          big_sym(i + start_index, j + start_index) += H(i, j);
      }
    }
    return 0;
  }



  // *****************************************************************************
  //                       BundleDenseTrustRegionProx::compute_QP_costs
  // *****************************************************************************


  int BundleDenseTrustRegionProx::compute_QP_costs(Symmatrix& Q,
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
    if (!fixed_values) {
      if (!is_factored) {
        is_factored = true;
        Hchol = H;
        for (Integer i = 0; i < Hchol.rowdim(); i++)
          Hchol(i, i) += weightu;
        if (Hchol.Chol_factor()) {
          if (cb_out())
            get_out() << "*** ERROR in BundleDenseTrustRegionProx::compute_QP_costs(...): H.Chol_factor failed" << std::endl;
          return 1;
        }
      }
      Hind_chol = Hchol;
    } else if ((fixed_changed) || (Hind_chol.rowdim() + old_fixed_ind.dim() != H.rowdim())) {
      Hind_chol = H;
      Hind_chol.delete_principal_submatrix(old_fixed_ind);
      for (Integer i = 0; i < Hind_chol.rowdim(); i++)
        Hind_chol(i, i) += weightu;
      if (Hind_chol.Chol_factor()) {
        if (cb_out())
          get_out() << "*** ERROR in BundleDenseTrustRegionProx::compute_QP_costs(...): Hind_chol.Chol_factor failed" << std::endl;
        return 1;
      }
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
    if (constant_minorant.get_minorant(_delta, _b, 0, 1., true, fixed_values ? &ind : 0, fixed_values ? &val : 0)) {
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
    LinvA.init(_A);
    Hind_chol.Chol_Lsolve(LinvA);
    rankadd(LinvA, Q, 1., 0., 1);
    Lty.init(_y);
    Hind_chol.Chol_Ltmult(Lty);
    Matrix Linvb(_b);
    Hind_chol.Chol_Lsolve(Linvb);
    offset = _delta - ip(Linvb, Linvb) / 2 + ip(_b, _y);
    Linvb -= Lty;
    d = _c;
    genmult(LinvA, Linvb, d, -1., 1., 1);

    oldQ = Q;
    oldd = d;
    oldoffset = offset;

    return 0;
  }


  // *****************************************************************************
  //                       BundleDenseTrustRegionProx::update_QP_costs
  // *****************************************************************************

  int BundleDenseTrustRegionProx::update_QP_costs(Symmatrix& delta_Q,
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
            get_out() << "*** ERROR in BundleDenseTrustRegionProx::update_QP_costs(...):  internal error, yfixed(" << ind << ")=" << (*yfixed)(ind) << " should not occur here" << std::endl;
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
    bool Hchanged = (new_fixed_ind.dim() > 0) || (Hind_chol.rowdim() + old_fixed_ind.dim() != H.rowdim());
    if (Hchanged) {
      //prepare Hind_chol
      if (old_fixed_ind.dim() + new_fixed_ind.dim() == 0) {
        if (!is_factored) {
          is_factored = true;
          Hchol = H;
          for (Integer i = 0; i < Hchol.rowdim(); i++)
            Hchol(i, i) += weightu;
          if (Hchol.Chol_factor()) {
            if (cb_out())
              get_out() << "*** ERROR in BundleDenseTrustRegionProx::update_QP_costs(...): Hchol.Chol_factor failed" << std::endl;
            return 1;
          }
        }
        Hind_chol = Hchol;
      } else {
        if (new_fixed_ind.dim() > 0) {
          old_fixed_ind.concat_below(new_fixed_ind);
          Indexmatrix sind;
          sortindex(old_fixed_ind, sind);
          old_fixed_ind = old_fixed_ind(sind);
          //update the values of _delta and _c by the newly fixed values
          _delta += ip(_b(new_fixed_newind), _y(new_fixed_newind));
          genmult(_A.rows(new_fixed_newind), _y(new_fixed_newind), _c, 1., 1., 1);
          //eliminate the newly fixed rows
          _A.delete_rows(new_fixed_newind);
          _b.delete_rows(new_fixed_newind);
          _y.delete_rows(new_fixed_newind);
        }
        Hind_chol = H;
        Hind_chol.delete_principal_submatrix(old_fixed_ind);
        for (Integer i = 0; i < Hind_chol.rowdim(); i++)
          Hind_chol(i, i) += weightu;
        if (Hind_chol.Chol_factor()) {
          if (cb_out())
            get_out() << "*** ERROR in BundleDenseTrustRegionProx::update_QP_costs(...): Hind_chol.Chol_factor failed" << std::endl;
          return 1;
        }
      }

      //recompute also Q 
      LinvA.init(_A);
      Hind_chol.Chol_Lsolve(LinvA);
      delta_Q.init(oldQ, -1.);
      rankadd(LinvA, oldQ, 1., 0., 1);
      delta_Q += oldQ;
      Lty.init(_y);
      Hind_chol.Chol_Ltmult(Lty);

    }//endif (Hchanged)
    else {
      delta_Q.init(xdim, 0.);
    }

    Matrix Linvb(_b);
    Hind_chol.Chol_Lsolve(Linvb);

    delta_offset = -oldoffset;
    oldoffset = _delta - ip(Linvb, Linvb) / 2 + ip(_b, _y);
    delta_offset += oldoffset;

    delta_d.init(oldd, -1.);
    Linvb -= Lty;
    oldd = _c;
    genmult(LinvA, Linvb, oldd, -1., 1., 1);
    delta_d += oldd;

    return 0;
  }


  // *****************************************************************************
  //                       BundleDenseTrustRegionProx::apply_modification
  // *****************************************************************************

  int BundleDenseTrustRegionProx::apply_modification(const GroundsetModification& gsmdf) {
    Integer olddim = H.rowdim();
    if (gsmdf.old_vardim() != olddim) {
      if (cb_out())
        get_out() << "**** ERROR BundleDenseTrustRegionProx::apply_modification: dim=" << olddim << " but modification assumes " << gsmdf.old_vardim() << std::endl;
      return 1;
    }

    H.enlarge_below(gsmdf.appended_vardim(), 0.);
    for (Integer i = olddim; i < H.rowdim(); i++)
      H(i, i) = weightu;

    if (gsmdf.map_to_old_variables()) {
      H = H.principal_submatrix(*(gsmdf.map_to_old_variables()));
    }
    compute_corr();

    Hind_chol.init(0, 0.);
    old_fixed_ind.init(0, 0, Integer(0));

    return 0;
  }


  // *****************************************************************************
  //         BundleDenseTrustRegionProx::apply_variable_metric
  // *****************************************************************************

  int BundleDenseTrustRegionProx::apply_variable_metric(VariableMetricModel* groundset,
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
        get_out() << "**** WARNING BundleDenseTrustRegionProx::apply_variable_metric(): this implementation does not support updates on subsets of indices -> updating everything" << std::endl;
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
      H.init(0, 0.);
    }
    if ((descent_step) || (H.rowdim() != y.dim())) {
      weightu = max(current_weight, 1e-10);
      H.init(y.dim(), 0.);
      new_indices = 0;
    }

    int err = 0;
    if (groundset->variable_metric_transform()->add_variable_metric(*this, y_id, y,
      descent_step,
      weightu,
      model_maxviol,
      in_new_indices)) {
      if (cb_out())
        get_out() << "**** WARNING BundleDenseTrustRegionProx::apply_variable_metric(): groundset->add_variable_metric(...) failed " << std::endl;
      err++;
    }
    if ((model) && (model->variable_metric_transform()->add_variable_metric(*this, y_id, y,
      descent_step,
      weightu,
      model_maxviol,
      in_new_indices))) {
      if (cb_out())
        get_out() << "**** WARNING BundleDenseTrustRegionProx::apply_variable_metric(): model->transform()->add_variable_metric(...) failed " << std::endl;
      err++;
    }

    current_weight = weightu;

    assert(aft_stack.size() == 0);
    new_indices = 0;
    aft_stack.clear();

    compute_corr();


    is_factored = false;
    Hind_chol.init(0, 0.);
    old_fixed_ind.init(0, 0, Integer(0));

    return err;
  }


  // *****************************************************************************
  //         BundleDenseTrustRegionProx::add_lowrank_variable_metric
  // *****************************************************************************

  int BundleDenseTrustRegionProx::add_variable_metric(Matrix& diagH,
    Matrix& vecH) {
    if ((vecH.coldim() == 0) && (diagH.dim() == 0))
      return 0;


    //apply the transposed trafos down the stack
    Real fun_factor = 1.;
    bool no_trafo = true;
    for (Integer i = Integer(aft_stack.size()); --i >= 0;) {
      fun_factor *= aft_stack[unsigned(i)]->get_fun_coeff();
      if (aft_stack[unsigned(i)]->get_arg_trafo()) {
        if (vecH.coldim() > 0) {
          genmult(*aft_stack[unsigned(i)]->get_arg_trafo(), vecH, LinvA, 1., 0., 1);
          swap(LinvA, vecH);
        }
        if (diagH.rowdim() > 0) {
          if (no_trafo) {
            no_trafo = false;
            LinvA.init(*aft_stack[unsigned(i)]->get_arg_trafo(), 1., 1);
            diagH.sqrt();
            LinvA.scale_cols(diagH);
          } else {
            genmult(*aft_stack[unsigned(i)]->get_arg_trafo(), diagH, LinvA, 1., 0., 1);
          }
          swap(LinvA, diagH);
        }
      }
    }

    //add vecH*diag(lam_vecH)*transpose(vecH) to H
    if (vecH.coldim() > 0)
      rankadd(vecH, H, fun_factor, 1.);
    if (diagH.rowdim() > 0) {
      if (!no_trafo) {
        rankadd(diagH, H, fun_factor, 1.);
      } else {
        for (Integer i = 0; i < H.rowdim(); i++)
          H(i, i) += fun_factor * diagH(i);
      }
    }

    return 0;
  }

  // *****************************************************************************
  //         BundleDenseTrustRegionProx::add_variable_metric
  // *****************************************************************************

  int BundleDenseTrustRegionProx::add_variable_metric(Symmatrix& addH) {
    if (addH.rowdim() == 0)
      return 0;

    assert((aft_stack.size() != 0) || (addH.rowdim() == H.rowdim()));


    //apply the transposed trafos down the stack
    Real fun_coeff = 1.;
    bool no_trafo = true;
    for (unsigned s = 0; s < aft_stack.size(); s++) {
      fun_coeff *= aft_stack[s]->get_fun_coeff();
      if (aft_stack[s]->get_arg_trafo()) {
        if (no_trafo) {
          addH.eig(LinvA, Lty);
          Lty.sqrt();
          LinvA.scale_cols(Lty);
          no_trafo = false;
        }
        genmult(*aft_stack[s]->get_arg_trafo(), LinvA, Lty, 1., 0., 1);
        swap(LinvA, Lty);
      }
    }
    if (!no_trafo) {
      rankadd(LinvA, H, fun_coeff, 1., 1);
    } else {
      H.xpeya(addH, fun_coeff);
    }

    return 0;
  }


  // *****************************************************************************
  //         BundleDenseTrustRegionProx::apply_Hinv
  // *****************************************************************************

  Matrix& BundleDenseTrustRegionProx::apply_Hinv(Matrix& x) const {
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

    //appply the inverse
    if (!is_factored) {
      is_factored = true;
      Hchol = H;
      for (Integer i = 0; i < Hchol.rowdim(); i++)
        Hchol(i, i) += weightu;
      if (Hchol.Chol_factor()) {
        if (cb_out())
          get_out() << "ERROR in BundleDenseTrustRegionProx::apply_Hinv(...): H.Chol_factor() failed" << std::endl;
        return x;
      }
    }
    Hchol.Chol_solve(x);

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
  //                       BundleDenseTrustRegionProx::mfile_data
  // *****************************************************************************


  int BundleDenseTrustRegionProx::mfile_data(std::ostream& out) const {
    out << "clear weightu qp_Hscale;\n";
    out << "weightu=";
    out.precision(16);
    out.width(18);
    out << weightu << "\n";
    out << "qp_Hscale=[";
    for (Integer i = 0; i < H.rowdim(); i++) {
      for (Integer j = 0; j < H.coldim(); j++) {
        out << " ";
        out.precision(16);
        out.width(18);
        out << H(i, j);
      }
      if (i < H.rowdim() - 1)
        out << "\n";
    }
    out << "];\n";
    return 0;
  }

}  //namespace ConicBundle
