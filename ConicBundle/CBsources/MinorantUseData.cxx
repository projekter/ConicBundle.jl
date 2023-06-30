/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/MinorantUseData.cxx
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



#include "MinorantUseData.hxx"

using namespace CH_Matrix_Classes;

namespace ConicBundle {


  // *****************************************************************************
  //                                 ~MinorantUseData
  // *****************************************************************************

  MinorantUseData::~MinorantUseData() {
    delete minorant;
    assert((md == 0) || (md->use_cnt > 0));
    if ((md) && (--(md->use_cnt) == 0)) {
      delete md;
    }
    md = 0;
  }

  // *****************************************************************************
  //                               get_scaleval_and_minorant
  // *****************************************************************************

  int MinorantUseData::get_scaleval_and_minorant(Real& sv, Minorant*& m) {
    int err = 0;
    if (md) {
      err = md->get_scaleval_and_minorant(sv, m);
      sv *= scaleval;
    } else {
      sv = scaleval;
      m = minorant;
    }
    return err;
  }

  // *****************************************************************************
  //                               valid
  // *****************************************************************************

  bool MinorantUseData::valid() const {
    return ((minorant != 0) && (modification_id >= 0)) ||
      ((md != 0) && (md->valid()));
  }


  // *****************************************************************************
  //                            get_minorant
  // *****************************************************************************

  Minorant* MinorantUseData::get_minorant() const {
    if (minorant)
      return minorant;
    if (md)
      return md->get_minorant();
    return 0;
  }

  // *****************************************************************************
  //                           set_modification_id
  // *****************************************************************************

  CH_Matrix_Classes::Integer& MinorantUseData::set_modification_id() {
    if (minorant)
      return modification_id;
    if (md)
      return md->set_modification_id();
    return modification_id;
  }

  // *****************************************************************************
  //                             get_modification_id
  // *****************************************************************************

  CH_Matrix_Classes::Integer MinorantUseData::get_modification_id() const {
    if (minorant)
      return modification_id;
    if (md)
      return md->get_modification_id();
    return -1;
  }

  // *****************************************************************************
  //                             get_modification_id
  // *****************************************************************************

  int MinorantUseData::synchronize_ids(Integer new_modification_id,
    Integer new_center_id,
    Integer old_center_id,
    Integer new_cand_id,
    Integer old_cand_id,
    Integer new_prex_id) {
    if (minorant) {
      modification_id = new_modification_id;
      prex_id = new_prex_id;
      std::map<Integer, Real>::iterator it;
      std::map<Integer, Real> tmp_evals;
      it = evals.find(old_center_id);
      if (it != evals.end())
        tmp_evals[new_center_id] = it->second;
      it = evals.find(old_cand_id);
      if (it != evals.end())
        tmp_evals[new_cand_id] = it->second;
      evals = tmp_evals;
      return 0;
    }
    if (md)
      return md->synchronize_ids(new_modification_id, new_center_id, old_center_id, new_cand_id, new_cand_id, new_prex_id);
    return 1;
  }


  // *****************************************************************************
  //                               aggregate
  // *****************************************************************************

  bool MinorantUseData::aggregate() const {
    if ((aggr_stat > 0) || ((minorant) && (minorant->number_aggregated() > 1)))
      return true;
    if (md)
      return md->aggregate();
    return false;
  }

  // *****************************************************************************
  //                              aggregated
  // *****************************************************************************

  int MinorantUseData::aggregated(CH_Matrix_Classes::Integer n) {
    if (md)
      return md->aggregated(n);
    if (minorant) {
      aggr_stat += n;
      return 0;
    }
    return 1;
  }

  // *****************************************************************************
  //                              offset
  // *****************************************************************************

  Real MinorantUseData::offset() const {
    if (minorant)
      return scaleval * minorant->offset();
    if (md)
      return scaleval * md->offset();
    return CB_minus_infinity;
  }

  // *****************************************************************************
  //                               coeff
  // *****************************************************************************

  Real MinorantUseData::coeff(CH_Matrix_Classes::Integer i) const {
    if (minorant)
      return scaleval * minorant->coeff(i);
    if (md)
      return scaleval * md->coeff(i);
    return CB_minus_infinity;
  }

  // *****************************************************************************
  //                               coeff
  // *****************************************************************************

  int MinorantUseData::scale(Real factor) {
    assert(use_cnt == 1);
    if (use_cnt != 1)
      return 1;
    Real d = factor * scaleval;
    scaleval = 1.;
    if (minorant) {
      if (d == 1.)
        return 0;
      else
        return minorant->scale_minorant(d);
    }
    if (md)
      return md->scale(d);
    return 1;
  }

  // *****************************************************************************
  //                            call_primal_extender
  // *****************************************************************************

  int MinorantUseData::call_primal_extender(PrimalExtender& prex, CH_Matrix_Classes::Integer in_prex_id) {
    if (minorant) {
      PrimalData* pd = minorant->get_primal();
      if ((pd == 0) || (prex_id < 0))
        return 1;
      if (in_prex_id > prex_id) {
        prex_id = in_prex_id;
        int status = prex.extend(*pd);
        if (status) {
          if (cb_out())
            get_out() << "**** WARNING: MinorantUseData::call_primal_extender(..): extend() failed" << std::endl;
          prex_id = -1;
          return status;
        }
      }
      return 0;
    }
    if (md)
      return md->call_primal_extender(prex, in_prex_id);
    return 1;
  }

  // *****************************************************************************
  //                               evaluate
  // *****************************************************************************

  Real MinorantUseData::evaluate(CH_Matrix_Classes::Integer yid, const Matrix& y, bool with_constant) const {
    if (minorant) {
      Real val;
      std::map<Integer, Real>::iterator it;
      if ((yid >= 0) && ((it = evals.find(yid)) != evals.end()))
        val = it->second;
      else {
        val = 0.;
        Integer n;
        const Real* cp;
        const Integer* ip;
        if (minorant->get_coeffs(n, cp, ip)) {
          return CB_minus_infinity;
        }
        if (ip == 0)
          val += mat_ip(n, cp, y.get_store());
        else {
          const Real* yp = y.get_store();
          for (; --n >= 0;) {
            val += (*cp++) * (*(yp + (*ip++)));
          }
        }
        if (yid >= 0)
          evals[yid] = val;
      }
      if (with_constant)
        val += minorant->offset();
      return val * scaleval;
    }
    if (md) {
      Real val = md->evaluate(yid, y, with_constant);
      if (val > CB_minus_infinity)
        return scaleval * val;
    }
    return CB_minus_infinity;
  }

  // *****************************************************************************
  //                               one_user
  // *****************************************************************************

  bool MinorantUseData::one_user() const {
    if (use_cnt > 1)
      return false;
    if (md)
      return md->one_user();
    return true;
  }



}

