/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/NNCBoxSupportFunction.cxx
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
#include "NNCBoxSupportFunction.hxx"

using namespace CH_Matrix_Classes;

namespace ConicBundle {

  //*****************************************************************************
  //                            NNCBoxSupportMinorantExtender
  //*****************************************************************************

  NNCBoxSupportMinorantExtender::NNCBoxSupportMinorantExtender(NNCBoxSupportFunction* in_fun) :
    CBout() {
    assert(in_fun);
    fun = in_fun;
    set_cbout(fun, 0);
  }

  NNCBoxSupportMinorantExtender::~NNCBoxSupportMinorantExtender() {
  }

  int NNCBoxSupportMinorantExtender::extend(Minorant& mnrt, int nc, const int* indices) {
    int err = 0;
    if (nc > 0) {
      for (int i = 0; i < nc; i++) {
        Integer ind = Integer(indices[i]);
        Real val = fun->get_lower_bounds()(ind);
        if (val <= 0.) {
          val = fun->get_upper_bounds()(ind);
          if (val >= 0.)
            val = 0.;
        }
        if (val != 0.) {
          if (mnrt.add_coeff(ind, val)) {
            err++;
            if (cb_out()) {
              get_out() << "**** WARNING: NNCBoxSupportMinorantExtender::extend(...): adding a new coefficint to the minorant failed for index=" << ind << std::endl;
            }
          }
        }
      }
    }

    return err;
  }

  //*****************************************************************************
  //                            NNCBoxSupportFunction
  //*****************************************************************************

  NNCBoxSupportFunction::NNCBoxSupportFunction(const Matrix& in_lb,
    const Matrix& in_ub,
    const CBout* cb,
    int incr) :CBout(cb, incr), lb(in_lb), ub(in_ub) {
    if (cb_out()) {
      if (lb.coldim() != 1)
        get_out() << "**** ERROR: NNCBoxSupportFunction::NNCBoxSupportFunction(....): lower bound vector has column dimension = " << lb.coldim() << "!= 1" << std::endl;
      if (ub.coldim() != 1)
        get_out() << "**** ERROR: NNCBoxSupportFunction::NNCBoxSupportFunction(....): upper bound vector has column dimension = " << ub.coldim() << "!= 1" << std::endl;
      if (lb.rowdim() != ub.rowdim())
        get_out() << "**** ERROR: NNCBoxSupportFunction::NNCBoxSupportFunction(....): number of rows of lower bound vector = " << lb.rowdim() << " != " << ub.rowdim() << " = number of rows of upper bound vector" << std::endl;
      if (max(lb - ub) > 0.)
        get_out() << "**** ERROR: NNCBoxSupportFunction::NNCBoxSupportFunction(....): some lower bound exceeds an upper bound by " << max(lb - ub) << std::endl;
    }
    assert(ub.coldim() == 1);
    assert(ub.coldim() == 1);
    assert(lb.rowdim() == ub.rowdim());
    assert((lb.rowdim() == 0) || (max(lb - ub) <= 0.));
  }

  int NNCBoxSupportFunction::evaluate(const Matrix& y,
    Real /* relprec */,
    Real& val,
    std::vector<Minorant*>& minorants,
    PrimalExtender*& primal_extender) {
    if (y.dim() != lb.rowdim()) {
      if (cb_out())
        get_out() << "**** ERROR: NNCBoxSupportFunction::evaluate(....): mismatch in argument dimension=" << y.dim() << " != " << lb.dim() << "= box dimension" << std::endl;
      return 1;
    }
    primal_extender = 0;
    Integer n = lb.rowdim();
    val = 0;
    const Real* yp = y.get_store();
    const Real* lp = lb.get_store();
    const Real* up = ub.get_store();
    Matrix x(n, 1); chk_set_init(x, 1);
    Real* xp = x.get_store();
    for (Integer i = 0; i < n; i++, lp++, up++) {
      Real d = *yp++;
      if (d < 0.) {
        val += d * (*lp);
        (*xp++) = *lp;
      } else if (d > 0.) {
        val += d * (*up);
        (*xp++) = *up;
      } else {
        if (*lp <= 0.) {
          if (*up >= 0.)
            (*xp++) = 0.;
          else
            (*xp++) = *up;
        } else
          (*xp++) = *lp;
      }
    }
    assert(minorants.size() == 0);
    minorants.push_back(new Minorant(true, 0., n, x.get_store()));
    return 0;
  }

  int NNCBoxSupportFunction::apply_modification(const NNCBoxSupportModification& mod) {
    int err = 0;

    if (mod.no_modification()) {
      return err;
    }

    if (mod.apply_to_bounds(lb, ub)) {
      err++;
      if (cb_out()) {
        get_out() << "**** ERROR in NNCBoxSupportFunction::apply_modification(.): modification failed" << std::endl;
      }
    }

    return err;
  }

  int NNCBoxSupportFunction::apply_modification(const OracleModification& omod,
    const Matrix* new_center,
    const Matrix* old_center,
    bool& discard_obj_in_center,
    bool& discard_model,
    bool& discard_aggregates,
    MinorantExtender*& extender) {
    int err = 0;

    NNCBoxSupportModification mod(0, this, 0);
    const NNCBoxSupportModification* modp = dynamic_cast<const NNCBoxSupportModification*>(&omod);
    if (modp == 0) {
      mod.clear(lb.rowdim());
      modp = &mod;
      if (mod.incorporate(omod)) {
        err++;
        if (cb_out()) {
          get_out() << "**** ERROR in NNCBoxSupportFunction::apply_modification(...): inocorporating a general oraclemodification failed" << std::endl;
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
        get_out() << "**** ERROR in NNCBoxSupportFunction::apply_modification(...): modifcation failed for the affine matrix function" << std::endl;
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
      extender = new NNCBoxSupportMinorantExtender(this);
    }

    return err;
  }


  void  NNCBoxSupportFunction::set_out(std::ostream* o, int pril) {
    CBout::set_out(o, pril);
  }

  void  NNCBoxSupportFunction::set_cbout(const CBout* cb, int incr) {
    CBout::set_cbout(cb, incr);
  }

  std::ostream& NNCBoxSupportFunction::print_problem_data(std::ostream& o) const {
    o << "\nBEGIN_NNCBOXSUPPORTFUNCTION\n";
    o << "\nBOUNDS\n";
    o.precision(12);
    o << " " << lb;
    o << " " << ub;
    o << "\nEND_NNCBOXSUPPORTFUNCTION" << std::endl;
    return o;
  }

  std::istream& NNCBoxSupportFunction::read_problem_data(std::istream& in) {
    char next_word[80];
    if (!in.good()) {
      if (cb_out()) {
        get_out() << "*** ERROR in NNCBoxSupportFunction::read_problem_data(): ";
        get_out() << " instream is not good";
      }
      return in;
    }
    in >> next_word;
    if (strcmp(next_word, "BEGIN_NNCBOXSUPPORTFUNCTION") != 0) {
      if (cb_out()) {
        get_out() << "*** ERROR in NNCBoxSupportFunction::read_problem_data(): ";
        get_out() << "expected BEGIN_NNCBOXSUPPORTFUNCTION but got " << next_word << std::endl;
      }
      in.clear(in.rdstate() | std::ios::failbit);
      return in;
    }
    in >> next_word;
    if (strcmp(next_word, "BOUNDS") != 0) {
      if (cb_out()) {
        get_out() << "*** ERROR in NNCBoxSupportFunction::read_problem_data(): ";
        get_out() << "expected BOUNDS but got " << next_word << std::endl;
      }
      in.clear(in.rdstate() | std::ios::failbit);
      return in;
    }
    in >> lb;
    if (!in.good()) {
      if (cb_out()) {
        get_out() << "*** ERROR in NNCBoxSupportFunction::read_problem_data(): ";
        get_out() << " instream is not good after reading lower bounds" << std::endl;
      }
      return in;
    }
    in >> ub;
    if (!in.good()) {
      if (cb_out()) {
        get_out() << "*** ERROR in NNCBoxSupportFunction::read_problem_data(): ";
        get_out() << " instream is not good after reading upper bounds" << std::endl;
      }
      return in;
    }
    if (max(lb - ub) > 0.) {
      if (cb_out()) {
        get_out() << "*** ERROR in SOCSupportFunction::read_problem_data(): ";
        get_out() << " a lower bound exceeds an upper bound" << std::endl;
      }
    }
    in >> next_word;
    if (strcmp(next_word, "END_NNCBOXSUPPORTFUNCTION") != 0) {
      if (cb_out()) {
        get_out() << "*** ERROR in NNCBoxSupportFunction::read_problem_data(): ";
        get_out() << "expected END_NNCBOXSUPPORTFUNCTION but got " << next_word << std::endl;
      }
      in.clear(in.rdstate() | std::ios::failbit);
      return in;
    }
    return in;
  }


  std::ostream& NNCBoxSupportFunction::print_problem_data_to_mfile(std::ostream& o, Integer block_nr) const {
    o << "\n% BEGIN_NNCBOXSUPPORTFUNCTION " << block_nr << "\n";
    o << "\n% BOUNDS within this block\n";
    o << "bounds{" << block_nr << ",1} = ";
    lb.mfile_output(o);
    o << "bounds{" << block_nr << ",2} = ";
    ub.mfile_output(o);
    o << "\n% END_NNCBOXSUPPORTFUNCTION " << block_nr << std::endl;
    return o;
  }




}
