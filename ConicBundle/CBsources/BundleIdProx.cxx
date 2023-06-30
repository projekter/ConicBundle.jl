/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/BundleIdProx.cxx
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
#include "BundleIdProx.hxx"
#include "mymath.hxx"

using namespace CH_Matrix_Classes;

namespace ConicBundle {


  // *****************************************************************************
  //                                add_H
  // *****************************************************************************

  int BundleIdProx::add_H(Symmatrix& big_sym,
    Integer start_index) const {
    assert((start_index >= 0) && (big_sym.rowdim() - start_index >= dim));
    for (Integer i = 0; i < dim; i++) {
      big_sym(i + start_index, i + start_index) += weightu;
    }
    return 0;
  }

  // *****************************************************************************
  //                       BundleIdProx::compute_QP_costs
  // *****************************************************************************


  int BundleIdProx::compute_QP_costs(Symmatrix& Q,
    Matrix& d,
    Real& offset,
    const MinorantPointer& constant_minorant,
    const MinorantBundle& bundle,
    const Matrix& y,
    const MinorantPointer& groundset_minorant,
    Indexmatrix* yfixed) {


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
    _A.newsize(ydim, xdim); chk_set_init(_A, 1);
    _b.newsize(ydim, 1); chk_set_init(_b, 1);
    _c.newsize(xdim, 1); chk_set_init(_c, 1);
    _delta = 0.;

    if (groundset_minorant.get_minorant(_delta, _b, 0, 1., false, fixed_values ? &ind : 0, fixed_values ? &val : 0)) {
      if (cb_out())
        get_out() << "*** ERROR in BundleIdProx::compute_QP_costs(...): groundset_minorant.get_minorant failed" << std::endl;
      return 1;
    }
    if ((!constant_minorant.empty()) && (constant_minorant.get_minorant(_delta, _b, 0, 1., true, fixed_values ? &ind : 0, fixed_values ? &val : 0))) {
      if (cb_out())
        get_out() << "*** ERROR in BundleIdProx::compute_QP_costs(...): groundset_minorant.get_minorant failed" << std::endl;
      return 1;
    }


    //Indexmatrix oldind(xdim,1,Integer(-1));
    //MinorantPointerMap::const_iterator it;
    for (Integer i = 0; i < xdim; i++) {
      //it=oldmap.find(bundle[unsigned(i)]);
      //if (it!=oldmap.end()){
      //  oldind(i)=it->second;
      //}
      if (bundle[unsigned(i)].get_minorant(_c(i), _A, i, 1, false, fixed_values ? &ind : 0, fixed_values ? &val : 0)) {
        if (cb_out())
          get_out() << "*** ERROR in BundleIdProx::compute_QP_costs(...): groundset_minorant.get_minorant failed" << std::endl;
        return 1;
      }
    }

    //Real weightfactor=oldweightu/weightu;
    Q.newsize(xdim); chk_set_init(Q, 1);
    for (Integer i = 0; i < xdim; i++) {
      // Integer ii=oldind(i);    //this concept never paid off
      // if (ii>=0){
      //   Q(i,i)=weightfactor*oldQ(ii,ii);
      //   std::cout<<" old("<<ii<<","<<ii<<")";
      //   for(Integer j=i+1;j<xdim;j++){
      // 	Integer jj=oldind(j);
      // 	if (jj>=0){
      // 	  Q(i,j)=weightfactor*oldQ(ii,jj);
      // 	  std::cout<<" old("<<ii<<","<<jj<<")";
      // 	}
      // 	else {
      // 	  Q(i,j)=bundle[unsigned(i)].ip(bundle[unsigned(j)],fixed_values?&ind:0)/weightu;
      // 	}
      //   }
      // }
      // else {
      for (Integer j = i; j < xdim; j++) {
        Q(i, j) = bundle[unsigned(i)].ip(bundle[unsigned(j)], fixed_values ? &ind : 0) / weightu;
      }
      //}
    }


    //rankadd(_A,oldQ,1./weightu,0.,1);
    //assert(norm2(Q-oldQ)<1e-8*(1.+max(diag(oldQ))));

    oldQ = Q;
    oldweightu = weightu;

    //Matrix _
    d = _c;
    Matrix tmpvec;
    tmpvec.init(_b, -1 / weightu);
    tmpvec += _y;
    genmult(_A, tmpvec, d, 1., 1., 1);
    //Real _
    offset = _delta - ip(_b, _b) / 2. / weightu + ip(_b, _y);


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

    // //Q=test_Q;
    // //d=test_d;
    // //offset=test_offset;

    // last_Q=Q;
    // last_d=d;
    // last_offset=offset;


    return 0;
  }


  // *****************************************************************************
  //                       BundleIdProx::update_QP_costs
  // *****************************************************************************

  int BundleIdProx::update_QP_costs(
    Symmatrix& delta_Q,
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

      if (yfixed) {
        switch ((*yfixed)(ind)) {
        case 0:
        {
          //was and stayed not fixed, but groundset_minorant may have changed
          if (uv != 0.) {
            delta_offset += uv * (yy - (oldbb + uv / 2.) / weightu);
            delta_d.xpeya(tmpvec, -uv / weightu);
          }
          break;
        }
        case 2:
        {
          //newly fixed, remove old eta (in this case -uv) and quadratic term
          delta_offset += (oldbb * oldbb) / 2. / weightu + uv * yy;
          delta_d.xpeya(tmpvec, oldbb / weightu);
          rankadd(tmpvec, delta_Q, -1. / weightu, 1.);
          (*yfixed)(ind) = 1;
          break;
        }
        default:
          if (cb_out())
            get_out() << "*** ERROR in BundleIdProx::update_QP_costs(...):  internal error, yfixed(" << ind << ")=" << (*yfixed)(ind) << " should not occur here" << std::endl;
          std::abort();
        }
      } else {
        if (uv != 0.) {
          delta_offset += uv * (yy - (oldbb + uv / 2.) / weightu);
          delta_d.xpeya(tmpvec, -uv / weightu);
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
    // std::cout<<"\n norm2(new_Q-test_Q)="<<norm2(last_Q+delta_Q-test_Q);
    // std::cout<<" norm2(new_d-test_d)="<<norm2(last_d+delta_d-test_d);
    // std::cout<<" fabs(new_offset-test_offset)="<<std::fabs(last_offset+delta_offset-test_offset);
    // std::cout<<std::endl;

    // delta_Q=test_Q-last_Q;
    // delta_d=test_d-last_d;
    // delta_offset=test_offset-last_offset;

    // last_Q=test_Q;
    // last_d=test_d;
    // last_offset=test_offset;

    return 0;
  }

  // *****************************************************************************
  //                       BundleIdProx::apply_modification
  // *****************************************************************************

  int BundleIdProx::apply_modification(const GroundsetModification& gsmdf) {
    if (gsmdf.old_vardim() != dim) {
      if (cb_out())
        get_out() << "**** ERROR BundleIdProx::apply_modification: dim=" << dim << " but modification assumes " << gsmdf.old_vardim() << std::endl;
      return 1;
    }
    dim = gsmdf.new_vardim();
    if (!gsmdf.no_modification()) {
      oldmap.clear();
      oldQ.init(0, 0.);
      old_fixed_ind.init(0, 0, Integer(0));
    }
    return 0;
  }


  // *****************************************************************************
  //                       BundleIdProx::mfile_data
  // *****************************************************************************


  int BundleIdProx::mfile_data(std::ostream& out) const {
    out << "clear qp_weight qp_diagscalr;\n";
    out << "qp_weight=" << weightu << ";\n";
    return 0;
  }



}  //namespace ConicBundle
