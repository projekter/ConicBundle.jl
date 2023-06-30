/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/SOCSupportFunction.cxx
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




#include <string.h>
#include "SOCSupportFunction.hxx"

using namespace CH_Matrix_Classes;

namespace ConicBundle {

  //*****************************************************************************
  //                            SOCSupportMinorantExtender
  //*****************************************************************************

  SOCSupportMinorantExtender::SOCSupportMinorantExtender(SOCSupportFunction* in_fun) :
    CBout() {
    assert(in_fun);
    fun = in_fun;
    set_cbout(fun, 0);
  }

  SOCSupportMinorantExtender::~SOCSupportMinorantExtender() {
  }

  int SOCSupportMinorantExtender::extend(Minorant& /* mnrt */,
    int /* nc */,
    const int* /* indices */) {
    //all new elements are zero automatically, nothing to do
    return 0;
  }

  //*****************************************************************************
  //                            SOCSupportFunction
  //*****************************************************************************

  SOCSupportFunction::SOCSupportFunction(Integer in_socdim,
    const CBout* cb,
    int incr) :CBout(cb, incr), socdim(in_socdim) {
    if (cb_out()) {
      if (socdim < 1)
        get_out() << "**** ERROR: SOCSupportFunction::SOCSupportFunction(....): second order cond dimension must be at least 1 but is " << socdim << std::endl;
    }
    assert(socdim >= 1);
  }

  Minorant* SOCSupportFunction::generate_minorant(const Matrix& SOCvec) {
    assert(socdim == SOCvec.dim());
    return new Minorant(true, 0., SOCvec.dim(), SOCvec.get_store());
  }

  int SOCSupportFunction::extract_SOCvector(Matrix& vec, const Minorant* mnrt) {
    assert(socdim == vec.dim());
    assert(mnrt);
    int n;
    const double* val;
    const int* ind;
    if (mnrt->get_coeffs(n, val, ind)) {
      if (cb_out())
        get_out() << "**** ERROR: SOCSupportFunction::extract_SOCvector(....): failed to get coefficients from the minorant" << std::endl;
      return 1;
    }
    if (ind == 0) {
      if (n > socdim) {
        if (cb_out())
          get_out() << "**** ERROR: SOCSupportFunction::extract_SOCvector(....): nonzero indices of the minorant exceed the second order cone dimension" << std::endl;
        return 1;
      }
      if (n == socdim) {
        vec.init(socdim, 1, val);
      } else {
        //fill up remainder;
        vec.newsize(socdim, 1);
        vec.init(n, 1, val);
        vec.enlarge_below(socdim - n, 0.);
      }
    } else {
      vec.init(socdim, 1, 0.);
      for (Integer i = 0; i < n; i++) {
        Integer indi = *ind++;
        if (indi >= socdim) {
          if (cb_out())
            get_out() << "**** ERROR: SOCSupportFunction::extract_SOCvector(....): nonzero indices of the minorant exceed the second order cone dimension" << std::endl;
          return 1;
        }
        vec(indi) = *val++;
      }
    }
    //check validity of the vector
    Real d = vec(0);
    vec(0) = 0.;
    Real n2 = norm2(vec);
    vec(0) = d;
    if ((d < 0) || (n2 > d * (1 + 1e-10))) {
      if (cb_out())
        get_out() << "**** ERROR: SOCSupportFunction::extract_SOCvector(....): no valid second order cone vector: zero component has value " << d << " and other compnonents have norm " << n2 << std::endl;
      return 1;
    }
    return 0;
  }

  int SOCSupportFunction::projection(Matrix& offset,
    Matrix& coeffs,
    const Matrix& bar_P,
    const Indexmatrix* index_subset) {
    assert(bar_P.rowdim() == socdim - 1);
    // c is zero
    offset.init(1, bar_P.coldim() + 1, 0.);
    // A is the identity
    if (index_subset == 0) {
      coeffs.newsize(socdim, bar_P.coldim() + 1); chk_set_init(coeffs, 1);
      Real* cop = coeffs.get_store();
      *cop++ = 1;
      Integer n = socdim - 1;
      mat_xea(n, cop, 0.);
      cop += n;
      const Real* pp = bar_P.get_store();
      for (Integer j = 0; j < bar_P.coldim(); j++) {
        *cop++ = 0.;
        mat_xey(n, cop, pp);
        cop += n;
        pp += n;
      }
    } else {
      Integer m = index_subset->dim();
      coeffs.newsize(m, bar_P.coldim() + 1); chk_set_init(coeffs, 1);
      for (Integer i = 0; i < m; i++) {
        Integer ind = (*index_subset)(i);
        if (ind == 0) {
          Real* cop = coeffs.get_store() + i;
          *cop = 1.;
          cop += m;
          mat_xea(bar_P.coldim(), cop, m, 0.);
        } else {
          Real* cop = coeffs.get_store() + i;
          *cop = 0.;
          cop += m;
          mat_xey(bar_P.coldim(), cop, m, bar_P.get_store() + ind - 1, socdim - 1);
        }
      }
    }

    // //--- TEST
    // Matrix P(bar_P.rowdim()+1,bar_P.coldim()+1,0.);
    // P(0,0)=1.;
    // P.subassign(Range(1,bar_P.rowdim()),Range(1,bar_P.coldim()),bar_P);
    // if (index_subset)
    //   P=P.rows(*index_subset);
    // Real testn2=norm2(coeffs-P);
    // std::cout<<" (TEST SOCSupportFunction::projection norm="<<testn2<<")"<<std::endl;
    // assert(testn2<1e-8);


    return 0;
  }

  int SOCSupportFunction::evaluate(const Matrix& y,
    Real /* relprec */,
    Real& val,
    Matrix& x,
    SOCPrimalExtender*& primal_extender) {
    if (y.dim() != socdim) {
      if (cb_out())
        get_out() << "**** ERROR: SOCSupportFunction::evaluate(.....): mismatch in argument dimension = " << y.dim() << " != " << socdim << " = second order cone dimension" << std::endl;
      return 1;
    }
    primal_extender = 0;
    x = y;
    x(0) = 0.;
    Real n2 = norm2(x);
    val = y(0) + n2;
    if (n2 > 1e-10)
      x /= n2;
    x(0) = 1.;
    return 0;
  }

  int SOCSupportFunction::evaluate_projection(const Matrix& y,
    const Matrix& P,
    Real /* relprec */,
    Real& val) {
    if (y.dim() != socdim) {
      if (cb_out())
        get_out() << "**** ERROR: SOCSupportFunction::evaluate_projection(....): mismatch in argument dimension = " << y.dim() << " != " << socdim << " = second order cone dimension" << std::endl;
      return 1;
    }
    if (P.rowdim() != socdim - 1) {
      if (cb_out())
        get_out() << "**** ERROR: SOCSupportFunction::evaluate_projection(....): row size of projection = " << P.rowdim() << " != " << socdim - 1 << " = second order cone dimension-1" << std::endl;
      return 1;
    }
    Matrix bary(socdim - 1, 1, y.get_store() + 1);
    Matrix barx;
    genmult(P, bary, barx, 1., 0., 1);
    Real n2 = norm2(barx);
    val = y(0) + n2;
    return 0;
  }

  int SOCSupportFunction::apply_modification(const SOCSupportModification& mod) {
    int err = 0;

    if (mod.no_modification()) {
      return err;
    }

    socdim = mod.new_vardim();

    if (socdim < 1) {
      err++;
      if (cb_out()) {
        get_out() << "**** ERROR in SOCSupportFunction::apply_modification(.): new second order cone dimension should be at least 1 but is " << socdim << std::endl;
      }
    }

    return err;
  }

  int SOCSupportFunction::apply_modification(const OracleModification& omod,
    const Matrix* new_center,
    const Matrix* old_center,
    bool& discard_obj_in_center,
    bool& discard_model,
    bool& discard_aggregates,
    MinorantExtender*& extender) {
    int err = 0;

    SOCSupportModification mod(0, this, 0);
    const SOCSupportModification* modp = dynamic_cast<const SOCSupportModification*>(&omod);
    if (modp == 0) {
      mod.clear(socdim);
      modp = &mod;
      if (mod.incorporate(omod)) {
        err++;
        if (cb_out()) {
          get_out() << "**** ERROR in SOCSupportFunction::apply_modification(...): inocorporating a general oraclemodification failed" << std::endl;
        }
      }
    }

    discard_obj_in_center = false;
    discard_model = false;
    discard_aggregates = false;
    extender = 0;

    if (modp->no_modification()) {
      return err;
    }

    //--- check whether deletions concern only zero variables/zero matrices
    if ((old_center == 0) || (new_center == 0)) {
      discard_obj_in_center = true;
    } else {
      if ((!modp->deleted_variables_are_zero(*old_center)) ||
        (!modp->mapped_variables_are_equal(*new_center, *old_center))) {
        discard_obj_in_center = true;
      }
    }

    //--- apply the modification
    if (apply_modification(*modp)) {
      err++;
      if (cb_out()) {
        get_out() << "**** ERROR in SOCSupportFunction::apply_modification(...): modifcation failed for the affine matrix function" << std::endl;
      }
    }

    //--- check whether additions have zero variables and require/allow extensions
    if (!discard_obj_in_center) {
      if (!modp->new_variables_are_zero(*new_center)) {
        discard_obj_in_center = true;
      }
    }

    //provide an extender if needed
    if ((!discard_model) && (!discard_aggregates)) {
      extender = new SOCSupportMinorantExtender(this);
    }

    return err;
  }


  void  SOCSupportFunction::set_out(std::ostream* o, int pril) {
    CBout::set_out(o, pril);
  }

  void  SOCSupportFunction::set_cbout(const CBout* cb, int incr) {
    CBout::set_cbout(cb, incr);
  }

  std::ostream& SOCSupportFunction::print_problem_data(std::ostream& o) const {
    o << "\nBEGIN_SOCSUPPORTFUNCTION\n";
    o << "\nDIMENSION\n";
    o.precision(12);
    o << " " << socdim;
    o << "\nEND_SOCSUPPORTFUNCTION" << std::endl;
    return o;
  }

  std::istream& SOCSupportFunction::read_problem_data(std::istream& in) {
    char next_word[80];
    if (!in.good()) {
      if (cb_out()) {
        get_out() << "*** ERROR in SOCSupportFunction::read_problem_data(): ";
        get_out() << " instream is not good";
      }
      return in;
    }
    in >> next_word;
    if (strcmp(next_word, "BEGIN_SOCSUPPORTFUNCTION") != 0) {
      if (cb_out()) {
        get_out() << "*** ERROR in SOCSupportFunction::read_problem_data(): ";
        get_out() << "expected BEGIN_SOCSUPPORTFUNCTION but got " << next_word << std::endl;
      }
      in.clear(in.rdstate() | std::ios::failbit);
      return in;
    }
    in >> next_word;
    if (strcmp(next_word, "DIMENSION") != 0) {
      if (cb_out()) {
        get_out() << "*** ERROR in SOCSupportFunction::read_problem_data(): ";
        get_out() << "expected DIMENSION but got " << next_word << std::endl;
      }
      in.clear(in.rdstate() | std::ios::failbit);
      return in;
    }
    in >> socdim;
    if (socdim < 1) {
      if (cb_out()) {
        get_out() << "*** ERROR in SOCSupportFunction::read_problem_data(): ";
        get_out() << " dimension should be at least 1 but is " << socdim << std::endl;
      }
    }
    in >> next_word;
    if (strcmp(next_word, "END_SOCSUPPORTFUNCTION") != 0) {
      if (cb_out()) {
        get_out() << "*** ERROR in SOCSupportFunction::read_problem_data(): ";
        get_out() << "expected END_SOCSUPPORTFUNCTION but got " << next_word << std::endl;
      }
      in.clear(in.rdstate() | std::ios::failbit);
      return in;
    }
    return in;
  }


  std::ostream& SOCSupportFunction::print_problem_data_to_mfile(std::ostream& o, Integer block_nr) const {
    o << "\n% BEGIN_SOCSUPPORTFUNCTION " << block_nr << "\n";
    o << "\n% DIMENSION within this block\n";
    o << "socdim{" << block_nr << "} = " << socdim;
    o << "\n% END_SOCSUPPORTFUNCTION " << block_nr << std::endl;
    return o;
  }




}
