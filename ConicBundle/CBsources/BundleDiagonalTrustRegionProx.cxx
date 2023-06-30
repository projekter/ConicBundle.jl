/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/BundleDiagonalTrustRegionProx.cxx
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
#include <iomanip>
#include "BundleDiagonalTrustRegionProx.hxx"
#include "mymath.hxx"
#include "sparssym.hxx"
#include "Groundset.hxx"
#include "BundleModel.hxx"


using namespace CH_Matrix_Classes;

namespace ConicBundle {


  // *****************************************************************************
  //                       BundleDiagonalTrustRegionProx::set_weightu
  // *****************************************************************************

  void BundleDiagonalTrustRegionProx::set_weightu(CH_Matrix_Classes::Real in_weightu) {
    if (fabs(weightu - in_weightu) < 1e-10 * fabs(weightu))
      return;

    D += (in_weightu - weightu);
    weightu = in_weightu;
    compute_corr();
    oldmap.clear();
    oldQ.init(0, 0.);
  }

  // *****************************************************************************
  //                       BundleDiagonalTrustRegionProx::norm_sqr
  // *****************************************************************************

  Real BundleDiagonalTrustRegionProx::norm_sqr(const Matrix& B) const {
    Real ipsum = normDsquared(B, D);
    assert(fabs(ipsum - ip(B, sparseDiag(D) * B)) < 1e-6);
    return ipsum;
  }


  // *****************************************************************************
  //                                add_H
  // *****************************************************************************

  int BundleDiagonalTrustRegionProx::add_H(Symmatrix& big_sym,
    Integer start_index) const {
    assert((start_index >= 0) && (big_sym.rowdim() - start_index >= D.rowdim()));
    for (Integer i = 0; i < D.rowdim(); i++)
      big_sym(i + start_index, i + start_index) += D(i);
    return 0;
  }

  // *****************************************************************************
  //                       BundleDiagonalTrustRegionProx::compute_QP_costs
  // *****************************************************************************


  int BundleDiagonalTrustRegionProx::compute_QP_costs(Symmatrix& Q,
    Matrix& d,
    Real& offset,
    const MinorantPointer& constant_minorant,
    const MinorantBundle& bundle,
    const Matrix& y,
    const MinorantPointer& groundset_minorant,
    Indexmatrix* yfixed) {
    //--- reset update_Dvalue
    update_Dvalue.init(0, 0, 0.);

    //--- get the bundle data into matrix form
    Indexmatrix ind(y.dim(), 1);
    ind.init(0, 0, 0.);
    Matrix val(y.dim(), 1);
    val.init(0, 0, 0.);
    Matrix _y(y.dim(), 1); chk_set_init(_y, 1);
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
    // bool fixed_changed=!equal(old_fixed_ind,ind);
    // if (fixed_changed){
    //   old_fixed_ind=ind;
    //   oldmap.clear();
    // }


    Integer xdim = Integer(bundle.size());
    Matrix _A(ydim, xdim); chk_set_init(_A, 1);
    Matrix _b(ydim, 1); chk_set_init(_b, 1);
    Matrix _c(xdim, 1); chk_set_init(_c, 1);
    Real _delta = 0.;

    if (groundset_minorant.get_minorant(_delta, _b, 0, 1., false, fixed_values ? &ind : 0, fixed_values ? &val : 0)) {
      if (cb_out())
        get_out() << "*** ERROR in BundleDiagonalTrustRegionProx::compute_QP_costs(...): groundset_minorant.get_minorant failed" << std::endl;
      return 1;
    }
    if ((!constant_minorant.empty()) && (constant_minorant.get_minorant(_delta, _b, 0, 1., true, fixed_values ? &ind : 0, fixed_values ? &val : 0))) {
      if (cb_out())
        get_out() << "*** ERROR in BundleDiagonalTrustRegionProx::compute_QP_costs(...): constant_minorant.get_minorant failed" << std::endl;
      return 1;
    }

    // Indexmatrix oldind(xdim,1,Integer(-1));
    // MinorantPointerMap::const_iterator it;
    for (Integer i = 0; i < xdim; i++) {
      // it=oldmap.find(bundle[unsigned(i)]);
      // if (it!=oldmap.end()){
      //  oldind(i)=it->second;
      // }
      if (bundle[unsigned(i)].get_minorant(_c(i), _A, i, 1, false, fixed_values ? &ind : 0, fixed_values ? &val : 0)) {
        if (cb_out())
          get_out() << "*** ERROR in BundleDiagonalTrustRegionProx::compute_QP_costs(...): bundle[" << i << "].get_minorant failed" << std::endl;
        return 1;
      }
    }

    //compute the linear and part of the constant coefficient
    Matrix tmpvec(D, -1.);
    if (fixed_values)
      tmpvec.delete_rows(ind);
    tmpvec.inv();
    tmpvec %= _b;
    offset = _delta + ip(_b, _y) + .5 * ip(tmpvec, _b);
    tmpvec += _y;
    d = _c;
    genmult(_A, tmpvec, d, 1., 1., 1);

    //compute the quadratic cost term
    tmpvec.init(D);
    if (fixed_values)
      tmpvec.delete_rows(ind);
    tmpvec.inv();
    tmpvec.sqrt();
    _A.scale_rows(tmpvec);

    Q.newsize(xdim); chk_set_init(Q, 1);
    for (Integer i = 0; i < xdim; i++) {
      // Integer ii=oldind(i);    //this concept never paid off
      // if (ii>=0){
      //   Q(i,i)=oldQ(ii,ii);
      //   for(Integer j=i+1;j<xdim;j++){
      // 	Integer jj=oldind(j);
      // 	if (jj>=0){
      // 	  Q(i,j)=oldQ(ii,jj);
      // 	}
      // 	else {
      // 	  Q(i,j)=mat_ip(ydim,_A.get_store()+i*ydim,_A.get_store()+j*ydim);
      // 	}
      //   }
      // }
      // else {
      for (Integer j = i; j < xdim; j++) {
        Q(i, j) = mat_ip(ydim, _A.get_store() + i * ydim, _A.get_store() + j * ydim);
      }
      // }
    }

    //rankadd(_A,oldQ,1.,0.,1);
    //assert(norm2(Q-oldQ)<1e-8*(1.+max(diag(oldQ))));

    oldQ = Q;

    // //---------- for testing purposes
    // Integer testdim=y.dim();
    // Real test_offset=0.;
    // Matrix testg(testdim,1,0.); 
    // groundset_minorant.get_minorant(test_offset,testg,0,1.);
    // if (!constant_minorant.empty())
    //   constant_minorant.get_minorant(test_offset,testg,0,1.,true);
    // Matrix testB(testdim,xdim,0.); chk_set_init(testB,1);
    // Matrix test_d(1,xdim,0.);  chk_set_init(test_d,1);
    // for (Integer i=0;i<xdim;i++){
    //   bundle[unsigned(i)].get_minorant(test_d(i),testB,i,1.);
    // }
    // test_offset+=ip(y,testg);
    // genmult(y,testB,test_d,1.,1.,1,0);
    // Matrix test_D(D);
    // test_D.sqrt();
    // test_D.inv();
    // testB.scale_rows(test_D);
    // testg.scale_rows(test_D);
    // Indexmatrix delind;
    // for (Integer i=0;i<testdim;i++){
    //   if ((yfixed)&&((*yfixed)(i)!=0)){
    //     delind.concat_below(i);
    //   }
    // }
    // if (delind.rowdim()>0){
    //   testB.delete_rows(delind);
    //   testg.delete_rows(delind);
    // }
    // Symmatrix test_Q;
    // rankadd(testB,test_Q,1./weightu,0.,1);
    // genmult(testg,testB,test_d,-1./weightu,1.,1,0);
    // test_offset-=ip(testg,testg)/weightu/2.;
    // test_d.transpose();
    // std::cout<<"\n norm2(test_Q-Q)="<<norm2(test_Q-Q);
    // std::cout<<" norm2(test_d-d)="<<norm2(test_d-d);
    // std::cout<<" fabs(test_offset-offset)="<<std::fabs(test_offset-offset);
    // std::cout<<" norm2(testB-_A)="<<norm2(testB-_A);
    // std::cout<<std::endl;

    // last_Q=Q;
    // last_d=d;
    // last_offset=offset;

    return 0;

  }


  // *****************************************************************************
  //                       BundleDiagonalTrustRegionProx::update_QP_costs
  // *****************************************************************************

  int BundleDiagonalTrustRegionProx::update_QP_costs(Symmatrix& delta_Q,
    Matrix& delta_d,
    Real& delta_offset,
    const MinorantPointer& constant_minorant,
    const MinorantBundle& bundle,
    const Matrix& y,
    const MinorantPointer& groundset_minorant,
    const MinorantPointer& delta_subg,
    const Indexmatrix& delta_index,
    Indexmatrix* yfixed) {
    Integer xdim = Integer(bundle.size());
    delta_d.init(xdim, 1, 0.);
    delta_Q.init(xdim, 0.);
    Matrix tmpvec(xdim, 1, 0.);
    delta_offset = delta_subg.offset();

    //--- for each y_index that needs updating 
    for (Integer j = 0; j < delta_index.dim(); j++) {

      Integer ind = delta_index[j];
      Real uv = delta_subg.coeff(ind);
      Real oldbb = groundset_minorant.coeff(ind) - uv;
      if (!constant_minorant.empty())
        oldbb += constant_minorant.coeff(ind);
      for (Integer i = 0; i < xdim; i++)
        tmpvec(i) = bundle[unsigned(i)].coeff(ind);
      Real yy = y(ind);
      Real Dnew = D(ind);

      if (yfixed) {
        switch ((*yfixed)(ind)) {
        case 0:
        {
          //was and stayed not fixed, but groundset_minorant and the quadratic term may have changed
          if ((update_Dvalue.rowdim() > 0) && (update_Dvalue(ind) != 0.)) {
            Real Dold = Dnew - update_Dvalue(ind);
            Real updD = (1. / Dnew - 1. / Dold);
            delta_offset -= updD * oldbb * oldbb / 2.;
            delta_d.xpeya(tmpvec, -oldbb * updD);
            rankadd(tmpvec, delta_Q, updD, 1.);
          }
          if (uv != 0.) {
            delta_offset += uv * (yy - (oldbb + uv / 2.) / Dnew);
            delta_d.xpeya(tmpvec, -uv / Dnew);
          }
          break;
        }
        case 2:
        {
          //newly fixed, remove old eta (in this case -uv) and quadratic term
          assert((update_Dvalue.rowdim() == 0) || (update_Dvalue(ind) == 0.));
          delta_offset += (oldbb * oldbb) / 2. / Dnew + uv * yy;
          delta_d.xpeya(tmpvec, oldbb / Dnew);
          rankadd(tmpvec, delta_Q, -1. / Dnew, 1.);
          (*yfixed)(ind) = 1;
          break;
        }
        default:
          if (cb_out())
            get_out() << "*** ERROR in BundleDiagonalTrustRegionProx::update_QP_cosgts(...):  internal error, yfixed(" << ind << ")=" << (*yfixed)(ind) << " should not occur here" << std::endl;
          std::abort();
        }
      } else {
        if ((update_Dvalue.rowdim() > 0) && (update_Dvalue(ind) != 0.)) {
          Real Dold = Dnew - update_Dvalue(ind);
          Real updD = (1. / Dnew - 1. / Dold);
          delta_offset -= updD * oldbb * oldbb / 2.;
          delta_d.xpeya(tmpvec, -oldbb * updD);
          rankadd(tmpvec, delta_Q, updD, 1.);
        }
        if (uv != 0.) {
          delta_offset += uv * (yy - (oldbb + uv / 2.) / Dnew);
          delta_d.xpeya(tmpvec, -uv / Dnew);
        }
      }
    }

    // //---------- for testing purposes
    // Integer testdim=y.dim();
    // Real test_offset=0.;
    // Matrix testg(testdim,1,0.); 
    // groundset_minorant.get_minorant(test_offset,testg,0,1.);
    // if (!constant_minorant.empty())
    //   constant_minorant.get_minorant(test_offset,testg,0,1.,true);
    // Matrix testB(testdim,xdim,0.); chk_set_init(testB,1);
    // Matrix test_d(1,xdim,0.);  chk_set_init(test_d,1);
    // for (Integer i=0;i<xdim;i++){
    //   bundle[unsigned(i)].get_minorant(test_d(i),testB,i,1.);
    // }
    // test_offset+=ip(y,testg);
    // genmult(y,testB,test_d,1.,1.,1,0);
    // Matrix test_D(D);
    // test_D.sqrt();
    // test_D.inv();
    // testB.scale_rows(test_D);
    // testg.scale_rows(test_D);
    // Indexmatrix delind;
    // for (Integer i=0;i<testdim;i++){
    //   if ((yfixed)&&((*yfixed)(i)!=0){
    //     delind.concat_below(i);
    //   }
    // }
    // if (delind.rowdim()>0){
    //   testB.delete_rows(delind);
    //   testg.delete_rows(delind);
    // }
    // Symmatrix test_Q;
    // rankadd(testB,test_Q,1./weightu,0.,1);
    // genmult(testg,testB,test_d,-1./weightu,1.,1,0);
    // test_offset-=ip(testg,testg)/weightu/2.;
    // test_d.transpose();
    // std::cout<<"\n norm2(test_Q-new_Q)="<<norm2(test_Q-last_Q-delta_Q);
    // std::cout<<" norm2(test_d-new_d)="<<norm2(test_d-last_d-delta_d);
    // std::cout<<" fabs(test_offset-new_offset)="<<std::fabs(test_offset-last_offset-delta_offset);
    // std::cout<<std::endl;

    // last_Q+=delta_Q;
    // last_d+=delta_d;
    // last_offset+=delta_offset;

    return 0;
  }


  // *****************************************************************************
  //                       BundleDiagonalTrustRegionProx::apply_modification
  // *****************************************************************************

  int BundleDiagonalTrustRegionProx::apply_modification(const GroundsetModification& gsmdf) {
    if (gsmdf.old_vardim() != D.rowdim()) {
      if (cb_out())
        get_out() << "**** ERROR BundleDiagonalTrustRegionProx::apply_modification: dim=" << D.rowdim() << " but modification assumes " << gsmdf.old_vardim() << std::endl;
      return 1;
    }
    oldmap.clear();
    oldQ.init(0, 0.);
    old_fixed_ind.init(0, 0, Integer(0));
    if (gsmdf.no_modification())
      return 0;
    D.concat_below(Matrix(gsmdf.appended_vardim(), 1, weightu));
    if (gsmdf.map_to_old_variables()) {
      D = D.rows(*(gsmdf.map_to_old_variables()));
    }
    compute_corr();
    return 0;
  }


  // *****************************************************************************
  //         BundleDiagonalTrustRegionProx::apply_variable_metric
  // *****************************************************************************

  int BundleDiagonalTrustRegionProx::apply_variable_metric(VariableMetricModel* /* groundset */,
    VariableMetricModel* model,
    const Matrix& in_aggr,
    Integer y_id,
    const Matrix& y,
    bool descent_step,
    Real& current_weight,
    Real model_maxviol,
    const Indexmatrix* in_new_indices) {
    assert(aft_stack.size() == 0);

    new_indices = in_new_indices;
    aggr = &in_aggr;

    if (y.dim() == 0)
      return 0;

    if ((weightu <= 0.) || (current_weight <= 0.)) {
      //not initialized, initialize
      if (current_weight > 0.) {
        weightu = current_weight;
      } else {
        if (weightu <= 0.)
          weightu = 1.;
        current_weight = weightu;
      }
      old_damping = 1.;
    }
    weightu = max(current_weight, 1e-10);
    if (D.dim() != y.dim()) {
      new_indices = 0;
    }

    D.init(y.dim(), 1, 0.);
    aft_stack.clear();


    int err = 0;
    if ((model) && (model->variable_metric_transform()->add_variable_metric(*this, y_id, y,
      descent_step,
      weightu,
      model_maxviol,
      in_new_indices))) {
      if (cb_out())
        get_out() << "**** WARNING BundleDiagonalTrustRegionProx::apply_variable_metric(): model->transform()->add_variable_metric(...) failed " << std::endl;
      err++;
    }


    if (old_D.dim() == D.dim()) {
      D *= .1;
      D.xpeya(old_D, .9);
      if (!descent_step) {
        for (Integer i = 0; i < D.rowdim(); i++)
          if (D(i) < old_D(i))
            D(i) = old_D(i);
      }
    }
    old_D = D;

    if (cb_out(2)) {
      get_out() << "\nBDTR(" << sum(D) / D.dim() << "," << max(D) << "," << min(D) << ")";
    }

    //damping
    bool damping = false;
    if (weightu < 1.) {
      if (descent_step) {
        D += weightu;
        Real daggrnorm = normDsquared(*aggr, D, 0, 1);
        Real aggrnorm = sqr(norm2(*aggr));
        D -= weightu;
        if (daggrnorm < .5 * aggrnorm) {
          damping = true;
          old_damping = daggrnorm / aggrnorm;
          D *= old_damping;
          if ((cb_out(2)) && (damping))
            get_out() << " damping=" << old_damping;
        } else
          old_damping = 1.;
      } else if (old_damping != 1.) {
        damping = true;
        D *= old_damping;
      }
    }

    D += weightu;


    if (cb_out(2)) {
      get_out() << std::endl;
    }

    current_weight = weightu;

    assert(aft_stack.size() == 0);
    new_indices = 0;
    aft_stack.clear();

    compute_corr();
    oldmap.clear();
    oldQ.init(0, 0.);
    return err;
  }



  // *****************************************************************************
  //         BundleDiagonalTrustRegionProx::add_variable_metric
  // *****************************************************************************

  int BundleDiagonalTrustRegionProx::add_variable_metric(Matrix& diagH,
    Matrix& vecH) {
    if ((diagH.dim() == 0) && (vecH.dim() == 0))
      return 0;

    //vecH is ignored but the diagonal values could be added if so desired
    if (vecH.dim() > 0) {
      //apply the transformations down the stack
      Real fun_factor = 1.;
      Matrix tmpmat;
      for (Integer i = Integer(aft_stack.size()); --i >= 0;) {
        fun_factor *= aft_stack[unsigned(i)]->get_fun_coeff();
        if (aft_stack[unsigned(i)]->get_arg_trafo()) {
          genmult(*aft_stack[unsigned(i)]->get_arg_trafo(), vecH, tmpmat, 1., 0., 1);
          swap(tmpmat, vecH);
        }
      }
      //take the diagonal of vecH*vecH^T  
      if (new_indices == 0) {
        Real aggrnorm = norm2(*aggr);
        Real dnorm = aggrnorm * aggrnorm / weightu;
        //Real dnormbound=(weightu>.5)?dnorm/2.:aggrsqu/2.;
        if (cb_out(2))
          get_out() << " BDTRadd " << dnorm << " ";
        Integer dim = D.dim();
        const Real* xp = vecH.get_store();
        const Real* const dpend = D.get_store() + dim;
        const Real* const agend = aggr->get_store() + dim;
        for (Integer j = 0; j < vecH.coldim(); j++) {

          const Real* ap = aggr->get_store();
          Real normvec = mat_ip(dim, xp, xp);
          Real ipvecval = sqr(mat_ip(dim, xp, ap));
          if (cb_out(2)) {
            get_out() << "(" << ipvecval * weightu / (weightu + normvec) << ")";
          }
          Real ipdiagval = 0.;
          while (ap != agend)
            ipdiagval += sqr((*xp++) * (*ap++));
          if (10 * ipvecval < ipdiagval)
            continue;
          Real factor = 1.;
          if (ipvecval < ipdiagval)
            factor = ipvecval / ipdiagval;
          if (cb_out(2)) {
            get_out() << factor << " ";
          }
          factor *= fun_factor;
          xp -= dim;
          Real* dp = D.get_store();
          while (dp != dpend) {
            Real s = *xp++;
            *dp++ += s * s * factor;
          }
        }
      } else {
        for (Integer j = 0; j < vecH.coldim(); j++) {
          for (Integer i = 0; i < new_indices->dim(); i++) {
            Integer ind = (*new_indices)(i);
            D(ind) += sqr(vecH(ind, j)) * fun_factor;
          }
        }
      }
    }

    //apply the transposed trafos down the stack
    if (diagH.dim() > 0) {
      Real fun_coeff = 1.;
      Integer trafos = -1;
      for (unsigned s = 0; s < aft_stack.size(); s++) {
        fun_coeff *= aft_stack[s]->get_fun_coeff();
        if (aft_stack[s]->get_arg_trafo()) {
          trafos = Integer(s);
        }
      }
      if (trafos >= 0) {
        Sparsemat sm = *(aft_stack[unsigned(trafos)]->get_arg_trafo());
        Indexmatrix rowinfo = sm.get_rowinfo();
        Matrix vec1, vec2;
        for (Integer j = 0; j < rowinfo.rowdim(); j++) {
          Integer ind = rowinfo(j, 0);
          if (diagH(ind) <= 0.)
            continue;
          vec1 = sm.row(ind);
          vec1.transpose();
          for (Integer s = trafos; --s >= 0;) {
            if (aft_stack[unsigned(s)]->get_arg_trafo()) {
              genmult(*aft_stack[unsigned(s)]->get_arg_trafo(), vec1, vec2, 1., 0., 1);
              swap(vec1, vec2);
            }
          }
          vec1 %= vec1;
          if (new_indices == 0)
            D.xpeya(vec1, fun_coeff * diagH(ind));
          else {
            Real d = fun_coeff * diagH(ind);
            for (Integer i = 0; i < new_indices->dim(); i++) {
              Integer k = (*new_indices)(i);
              D(k) += vec1(k) * d;
            }
          }
        }
      } else {
        if (new_indices == 0)
          D.xpeya(diagH, fun_coeff);
        else {
          for (Integer i = 0; i < new_indices->dim(); i++) {
            Integer k = (*new_indices)(i);
            D(k) += diagH(k) * fun_coeff;
          }
        }
      }
    }

    return 0;
  }

  // *****************************************************************************
  //         BundleDiagonalTrustRegionProx::apply_Hinv
  // *****************************************************************************

  Matrix& BundleDiagonalTrustRegionProx::apply_Hinv(Matrix& x) const {
    if ((aft_stack.size() == 0) && (x.coldim() == 1)) {
      return x /= D;
    }

    //apply the transposed trafos down the stack
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
    Real* xp = x.get_store();
    const Real* const dpend = D.get_store() + D.dim();
    for (Integer i = 0; i < x.coldim(); i++) {
      const Real* dp = D.get_store();
      while (dp != dpend)
        *xp++ /= *dp++;
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
  //                       BundleDiagonalTrustRegionProx::mfile_data
  // *****************************************************************************


  int BundleDiagonalTrustRegionProx::mfile_data(std::ostream& out) const {
    out << "clear eta yfixed qp_weight qp_diagscalr;\n";
    out << "qp_diagscale=[";
    for (Integer i = 0; i < D.dim(); i++) {
      out.precision(16);
      out.width(18);
      out << D(i);
      if (i < D.dim() - 1)
        out << "\n";
    }
    out << "];\n";
    return 0;
  }


}  //namespace ConicBundle
