/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/Minorant.cxx
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




#include "CBSolver.hxx"
#include "matrix.hxx"

#include <algorithm>
#include <map>



using namespace CH_Matrix_Classes;

namespace ConicBundle {


  /** @brief serves as the default tolerance for considering minorant entries as zero
   */
   //const double CB_minorant_zero_tolerance=eps_Real;
  const double CB_minorant_zero_tolerance = 1e-100;

  /** @brief serves as the default ratio of nonzeros to dimension for using a sparse
      representatio of a minorant
   */
  const double CB_minorant_sparsity_ratio = 0.3;
  //define CB_minorant_sparsity_ratio 0.2
  //define CB_minorant_zero_tolerance  1e-12

  //------------------------------------------------------------
  // Data class 
  //------------------------------------------------------------
  /**@brief in order to keep the interface 'clean' the interface's data is separated into an extra class

     If sparse==false, all is stored in vecval from index 0 up to
     and including maxindex (vecval may be larger).

     If sparse==false and in addition clean==true, then nz_cnt has already
     been computed. It might be that the indices vecind are not computed if
     the vector is clearly dense, then vecind.dim()==0. If vecind.dim()==nz_cnt
     then the indices in vecind are correct.

     if sparse==true all relevant data is stored in vecval and vecind from
     index 0 to nz_cnt-1 where vecval[i] belongs to index vec_ind[i]
     and vecval.dim()=vecind.dim()=nz_cnt,
     but if clean==false the entries need not be sorted and may have
     duplicates. If clean==true then the entries of vecind are in strictly
     increasing order.

   */
  class Minorant::MinorantData {
  protected:
    virtual ~MinorantData();


    //----------------------------------------
    // data
    Real offset;
    Matrix vecval;
    Indexmatrix vecind;
    PrimalData* primal;
    Real normsqu; // if <0 not computed, if >0 the current Euclidean norm squared

    Integer n_aggregated;

    bool offset_at_origin;
    Integer maxindex;
    Integer nz_cnt;
    bool sparse;
    bool clean;


    MinorantData(const Minorant::MinorantData& md, Real factor = 1., bool with_primal = false) :
      offset(md.offset* factor), vecval(md.vecval, factor), vecind(md.vecind),
      primal(0), normsqu(md.normsqu* factor* factor),
      n_aggregated(md.n_aggregated),
      offset_at_origin(md.offset_at_origin), maxindex(md.maxindex),
      nz_cnt(md.nz_cnt), sparse(md.sparse), clean(md.clean) {
      if ((with_primal) && (md.primal)) {
        primal = md.primal->clone_primal_data();
        if ((primal) && (factor != 1.))
          primal->scale_primal_data(factor);
      }
    }


    MinorantData(bool oao, Real ofs, Integer rne) {
      offset = ofs;
      offset_at_origin = oao;
      maxindex = -1;
      rne = max(rne, 0);
      vecval.newsize(rne, 1);
      vecval.init(0, 1, 0.);
      vecind.init(0, 1, Integer(0));
      nz_cnt = 0;
      n_aggregated = 1;
      sparse = false;
      clean = true;
      primal = 0;
      normsqu = 0.;
    }

  private:
    friend class Minorant;

    double coeff(Integer ind);
    int make_clean(Real tol = CB_minorant_zero_tolerance, Real spratio = CB_minorant_sparsity_ratio);
    int make_dense(Integer maxindex);
    int prepare_for_changes(Integer maxindex, Integer new_nz);
    int reassign_coeffs(Integer nel, const Integer* map_to_old);
    int scale(Real factor);
    Real norm_squared();

    MinorantData& operator= (const Minorant::MinorantData& md) {
      offset = md.offset;
      vecval = md.vecval;
      vecind = md.vecind;
      offset_at_origin = md.offset_at_origin;
      maxindex = md.maxindex;
      nz_cnt = md.nz_cnt;
      sparse = md.sparse;
      clean = md.clean;
      delete primal;
      primal = 0;
      normsqu = md.normsqu;
      n_aggregated = md.n_aggregated;
      return *this;
    }

    MinorantData& init(const Minorant::MinorantData& md, double factor) {
      if (factor == 1.) {
        *this = md;
      } else if (factor == 0.) {
        offset = 0.;
        offset_at_origin = true;
        maxindex = -1;
        vecval.newsize(md.vecval.dim(), 1);
        vecval.init(0, 1, 0.);
        vecind.init(0, 1, Integer(0));
        nz_cnt = 0;
        sparse = false;
        clean = true;
        delete primal;
        primal = 0;
        normsqu = 0.;
        n_aggregated = md.n_aggregated;
      } else {
        offset = md.offset * factor;
        vecval.init(md.vecval, factor);
        vecind = md.vecind;
        offset_at_origin = md.offset_at_origin;
        maxindex = md.maxindex;
        nz_cnt = md.nz_cnt;
        sparse = md.sparse;
        clean = md.clean;
        delete primal;
        primal = 0;
        normsqu = md.normsqu * factor * factor;
        n_aggregated = md.n_aggregated;
      }
      return *this;
    }

  };


  Minorant::MinorantData::~MinorantData() {
    delete primal;
  }

  // *****************************************************************************
  //                               Minorant::MinorantData::make_clean
  // *****************************************************************************

  int Minorant::MinorantData::make_clean(Real tol, Real spratio) {
    if (clean) {
      if (!sparse) {
        //nz_cnt is already determined but 
        if (nz_cnt < (maxindex + 1) * spratio) {
          //switch to sparse
          if (vecind.rowdim() == nz_cnt) {
            //all data available, shrink vecval;
            assert(vecval.rowdim() >= maxindex + 1);
            Real* vval = vecval.get_store();
            Real* const vstart = vval;
            const Integer* vind = vecind.get_store();
            for (Integer i = nz_cnt; --i >= 0; vind++) {
              Real d = *(vstart + *vind);
              *(vstart + *vind) = 0.;
              *vval++ = d;
            }
            vecval.reduce_length(nz_cnt);
            sparse = true;
            return 0;
          }
          //here we need to compute vecind, this is done below, continue there
        } else {
          return 0;
        }
      } else {
        return 0;
      }
    }

    if (sparse) {
      assert(!clean);
      //not clean, so some indices may be out of order or may have duplicates
      const Real abstol = tol * (1. + std::fabs(offset));
      vecind.reduce_length(nz_cnt);
      Indexmatrix sind;
      sortindex(vecind, sind);
      Matrix val(maxindex + 1, 1); chk_set_init(val, 1);
      Real* valp = val.get_store();
      Integer* indp = sind.get_store();
      Integer cnt = 0;
      const Integer* ip = sind.get_store();
      const Integer* const ipend = ip + nz_cnt;
      const Integer* const vecindp = vecind.get_store();
      const Real* const vecvalp = vecval.get_store();
      Integer oldi = -1;
      while (ip != ipend) {
        Integer newi = *(vecindp + *ip);
        assert(newi >= 0);
        if (newi == oldi) {
          *valp += vecval(*ip);
        } else {
          if ((oldi >= 0) && (fabs(*valp) > abstol)) {
            cnt++;
            valp++;
            indp++;
          }
          *valp = *(vecvalp + *ip);
          *indp = newi;
          oldi = newi;
        }
        ip++;
      }
      if ((oldi >= 0) && (fabs(*valp) > abstol)) {
        cnt++;
      }
      nz_cnt = cnt;
      swap(vecval, val);
      vecval.reduce_length(nz_cnt);
      swap(vecind, sind);
      vecind.reduce_length(nz_cnt);
      if (nz_cnt > 0) {
        maxindex = vecind(nz_cnt - 1);
      } else {
        maxindex = -1;
      }
      if (nz_cnt >= (maxindex + 1) * spratio) {
        //switch to dense representation
        vecval.enlarge_below(maxindex + 1 - nz_cnt, 0.);
        const Integer* const vecindp = vecind.get_store();
        Real* const vecvalp = vecval.get_store();
        for (Integer i = nz_cnt; --i >= 0;) {
          Integer ind = *(vecindp + i);
          if (ind == i)
            break;
          *(vecvalp + ind) = *(vecvalp + i);
          *(vecvalp + i) = 0.;
        }
        sparse = false;
      }
      clean = true;
      return 0;
    }

    // ---  check whether it is worth to switch from dense to sparse -----
    if (!clean) {
      //determine the number of nonzeros  nz_cnt 
      Integer cnt = 0;
      const Real* vval = vecval.get_store();
      const Real* const vend = vval + maxindex + 1;
      const Real abstol = tol * (1. + std::fabs(offset));
      while (vval != vend) {
        if (std::fabs(*vval++) > abstol) {
          cnt++;
        }
      }
      nz_cnt = cnt;
      vecind.init(0, 1, Integer(0));
    }
    if (nz_cnt < (maxindex + 1) * spratio) {
      //switch from dense to sparse
      //for this determine the indices first
      if (nz_cnt != vecind.rowdim()) {
        vecind.newsize(nz_cnt, 1); chk_set_init(vecind, 1);
      }
      const Real* vval = vecval.get_store();
      const Real* const vend = vval + maxindex + 1;
      const Real abstol = tol * (1. + std::fabs(offset));
      Integer* vind = vecind.get_store();
      const Integer* const viend = vind + nz_cnt;
      Integer i = 0;
      //if the first few elements are all nonzeros skip them
      while ((vval != vend) && (std::fabs(*vval) > abstol)) {
        *vind++ = i++;
        vval++;
      }
      //then skip nonzeros and shift the next ones to initial positions
      Real* vv = const_cast<Real*>(vval);
      Real* vold = vv;
      while (vind != viend) {
        while (!(std::fabs(*vold) > abstol)) {
          vold++; i++;
        }
        *vind++ = i++;
        *vv++ = *vold;
        *vold++ = 0.;
      }
      if (nz_cnt > 0) {
        maxindex = *(--vind);
      } else {
        maxindex = -1;
      }
      vecval.reduce_length(nz_cnt);
      sparse = true;
    }
    clean = true;
    return 0;
  }

  // *****************************************************************************
  //                               Minorant::MinorantData::make_dense
  // *****************************************************************************

  int Minorant::MinorantData::make_dense(Integer new_maxind) {
    Integer old_rowdim = vecval.rowdim();
    if ((old_rowdim <= new_maxind) && ((!sparse) || (clean))) {
      vecval.enlarge_below(new_maxind + 1 - old_rowdim, 0.);
    }
    if (!sparse)
      return 0;
    //--- sparse case ---
    sparse = false;
    if (clean) {
      if (old_rowdim <= maxindex)
        mat_xea(old_rowdim - nz_cnt, vecval.get_store() + nz_cnt, 0.);
      Integer* ind = vecind.get_store() + nz_cnt;
      Real* val = vecval.get_store() + nz_cnt;
      Real* const vstart = vecval.get_store();
      const Integer* const iend = vecind.get_store();
      while (ind != iend) {
        Real d = *(--val);
        *val = 0.;
        *(vstart + (*(--ind))) = d;
      }
      return 0;
    }
    //not clean and will stay not clean
    Matrix val(new_maxind + 1, 1, 0.);
    for (Integer i = 0; i < nz_cnt; i++) {
      val(vecind(i)) += vecval(i);
    }
    swap(val, vecval);
    return 0;
  }

  // *****************************************************************************
  //                               Minorant::MinorantData::prepare_for_changes
  // *****************************************************************************

  int Minorant::MinorantData::prepare_for_changes(Integer new_maxind, Integer new_nz) {
    if (!sparse) {
      //
      if ((new_maxind + 1) * CB_minorant_sparsity_ratio > max(maxindex, new_nz)) {
        //switch to sparse
        if ((!clean) || (nz_cnt != vecind.dim())) {
          vecind.newsize(max(new_maxind, maxindex), 1); chk_set_init(vecind, 1);
          const Real abstol = CB_minorant_zero_tolerance * (1. + std::fabs(offset));
          const Real* vval = vecval.get_store();
          Integer cnt = 0;
          for (Integer i = 0; i <= maxindex; i++, vval++) {
            if (std::fabs(*vval) > abstol) {
              if (i > cnt) {
                vecval[cnt] = *vval;
              }
              vecind[cnt] = i;
              cnt++;
            }
          }
          nz_cnt = cnt;
          vecind.reduce_length(nz_cnt);
          vecval.reduce_length(nz_cnt);
        } else {
          for (Integer i = 0; i < nz_cnt; i++) {
            vecval[i] = vecval[vecind[i]];
          }
          vecval.reduce_length(nz_cnt);
        }
        sparse = true;
      } else {
        if (new_maxind >= vecval.dim()) {
          vecval.enlarge_below(new_maxind + 1 - vecval.dim(), 0.);
        }
      }
      clean = false;
      return 0;
    }
    //--- sparse case ---
    if (new_maxind < maxindex) {
      new_maxind = maxindex;
    }
    int err = 0;
    if ((new_maxind + 1) * CB_minorant_sparsity_ratio < new_nz + nz_cnt) {
      err = make_dense(new_maxind);
    }
    clean = false;
    return err;
  }


  // *****************************************************************************
  //                               Minorant::MinorantData::coeff
  // *****************************************************************************

  double Minorant::MinorantData::coeff(Integer ind) {
    if (!clean)
      make_clean();
    if (ind > maxindex)
      return 0.;
    if (sparse) {
      Integer ub = nz_cnt;
      Integer lb = 0;
      const Integer* const vind = vecind.get_store();
      Integer mid;
      Integer k;
      while (lb < ub) {
        mid = (ub + lb) / 2;
        k = *(vind + mid);
        if (k < ind)
          lb = mid + 1;
        else if (k > ind)
          ub = mid;
        else break;
      }
      if (lb < ub)
        return vecval(mid);
      return 0.;
    } else {
      return vecval(ind);
    }
    return 0.;
  }

  // *****************************************************************************
  //                               Minorant::MinorantData::reassign_coeffs
  // *****************************************************************************

  int Minorant::MinorantData::reassign_coeffs(Integer nel, const Integer* map_to_old) {
    if ((!clean) && (make_clean()))
      return 1;
    normsqu = -1.;
    Matrix tmpval(nel, 1); chk_set_init(tmpval, 1);
    Real* tval = tmpval.get_store();
    Integer mi = -1;
    if (sparse) {
      for (Integer i = 0; i < nel; i++) {
        Real d = coeff(*map_to_old++);
        *tval++ = d;
        if (d != 0.)
          mi = i;
      }
    } else {
      Real* tval = tmpval.get_store();
      const Real* const vval = vecval.get_store();
      for (Integer i = 0; i < nel; i++) {
        Integer ind = *map_to_old++;
        if (ind > maxindex)
          *tval++ = 0.;
        else {
          *tval++ = *(vval + ind);
          mi = i;
        }
      }
    }
    swap(tmpval, vecval);
    maxindex = mi;
    clean = false;
    sparse = false;

    return 0;
  }

  // *****************************************************************************
  //                               Minorant::MinorantData::reassign_coeffs
  // *****************************************************************************

  int Minorant::MinorantData::scale(Real factor) {
    offset *= factor;
    if (normsqu > 0.)
      normsqu *= factor * factor;
    if (primal)
      primal->scale_primal_data(factor);
    if (sparse) {
      mat_xmultea(nz_cnt, vecval.get_store(), factor);
    } else {
      mat_xmultea(maxindex + 1, vecval.get_store(), factor);
    }
    return 0;
  }

  // *****************************************************************************
  //                               Minorant::MinorantData::norm_squared
  // *****************************************************************************

  Real Minorant::MinorantData::norm_squared() {
    if (normsqu < 0.) {
      if (!clean)
        make_clean();
      if (sparse)
        normsqu = mat_ip(nz_cnt, vecval.get_store());
      else
        normsqu = mat_ip(maxindex + 1, vecval.get_store());
    }
    return normsqu;
  }





  //------------------------------------------------------------
  // Minorant implementation  
  //------------------------------------------------------------

// *****************************************************************************
//                                 Minorant
// *****************************************************************************

  Minorant::Minorant(bool oao, double offset,
    int ne, const double* coeffs, const int* indices,
    double scale_val, PrimalData* primal) {
    data = new Minorant::MinorantData(oao, scale_val * offset, ne);
    if ((ne > 0) && (coeffs != 0))
      add_coeffs(ne, coeffs, indices, scale_val);
    set_primal(primal);
  }

  // *****************************************************************************
  //                                 Minorant
  // *****************************************************************************

  Minorant::Minorant(double offset, const DVector& subg, PrimalData* primal, bool oao) {
    data = new Minorant::MinorantData(oao, offset, int(subg.size()));
    if (subg.size() > 0)
      add_coeffs(int(subg.size()), &(subg[0]));
    set_primal(primal);
  }

  // *****************************************************************************
  //                                 Minorant
  // *****************************************************************************

  Minorant::Minorant(double offset, const DVector& subg, const IVector& ind, PrimalData* primal, bool oao) {
    assert(subg.size() == ind.size());
    data = new Minorant::MinorantData(oao, offset, int(subg.size()));
    if (subg.size() > 0)
      add_coeffs(int(subg.size()), &(subg[0]), &(ind[0]));
    set_primal(primal);
  }

  // *****************************************************************************
  //                                Minorant
  // *****************************************************************************

  Minorant::Minorant(const Minorant* mnrt, double factor, bool with_primal) {
    data = new Minorant::MinorantData(*(mnrt->data), Real(factor), with_primal);
  }

  // *****************************************************************************
  //                                 ~Minorant
  // *****************************************************************************

  Minorant::~Minorant() {
    delete data;
  }

  // *****************************************************************************
  //                               offset
  // *****************************************************************************

  double Minorant::offset() const {
    return data->offset;
  }

  // *****************************************************************************
  //                               add_offset
  // *****************************************************************************

  int Minorant::add_offset(double offset) {
    data->offset += offset;
    return 0;
  }

  // *****************************************************************************
  //                           offset_gives_value_at_origin
  // *****************************************************************************

  bool& Minorant::offset_gives_value_at_origin() {
    return data->offset_at_origin;
  }

  // *****************************************************************************
  //                        offset_gives_value_at_origin
  // *****************************************************************************

  bool Minorant::offset_gives_value_at_origin() const {
    return data->offset_at_origin;
  }

  // *****************************************************************************
  //                               coeff
  // *****************************************************************************

  double Minorant::coeff(int ind) {
    assert(ind >= 0);
    return data->coeff(ind);
  }

  // *****************************************************************************
  //                               coeff
  // *****************************************************************************

  int Minorant::add_coeff(int ind, double val) {
    assert(ind >= 0);
    data->normsqu = -1.;
    Integer new_maxind = max(data->maxindex, ind);
    int err = data->prepare_for_changes(new_maxind, 1);
    assert(!data->clean);
    if (data->sparse) {
      assert((data->nz_cnt == data->vecind.dim()) && (data->nz_cnt == data->vecval.dim()));
      data->vecind.concat_below(ind);
      data->vecval.concat_below(val);
      data->nz_cnt++;
    } else {
      data->vecval(ind) += val;
    }
    data->maxindex = new_maxind;
    return err;
  }

  // *****************************************************************************
  //                               add_coeffs
  // *****************************************************************************

  int Minorant::add_coeffs(int nel,
    const double* coeff,
    double factor,
    int start_pos) {
    if (nel <= 0)
      return 0;
    data->normsqu = -1.;
    if (coeff == 0) {
      // use this for reserving space
      if (data->vecval.dim() < nel) {
        Integer n = data->vecval.dim();
        data->vecval.enlarge_below(nel - n, 0.);
        data->vecval.reduce_length(n);
      }
      return 0;
    }
    assert(coeff);
    assert(start_pos >= 0);
    Integer new_maxind = max(nel + start_pos - 1, data->maxindex);
    data->make_dense(new_maxind);
    data->vecind.init(0, 1, Integer(0));
    data->clean = false;
    Real* vval = data->vecval.get_store() + start_pos;
    const Real* cval = coeff;
    if (start_pos <= data->maxindex) {
      const Real* const vend = vval + min(nel, data->maxindex + 1 - start_pos);
      while (vval != vend)
        (*vval++) += factor * (*cval++);
    }
    const Real* const vend = data->vecval.get_store() + start_pos + nel;
    while (vval != vend)
      (*vval++) = factor * (*cval++);
    data->maxindex = new_maxind;
    return 0;
  }

  // *****************************************************************************
  //                               add_coeffs
  // *****************************************************************************

  int Minorant::add_coeffs(int nel,
    const double* coeff,
    const int* indices,
    double factor) {
    if ((nel <= 0) || (factor == 0.))
      return 0;
    if (indices == 0) {
      return add_coeffs(nel, coeff, factor, 0);
    }
    assert(coeff);
    data->normsqu = -1.;

    //--- check whether a clean state can be kept and find new maxind
    Integer new_maxind = data->maxindex;
    const Integer* const ipend = indices + nel;
    const Integer* ip = indices;
    if ((data->clean == true) && ((data->sparse) || (data->vecval.dim() == 0))) {
      while ((ip != ipend) && (new_maxind < *ip)) {
        new_maxind = *ip++;
        //we do not check coeff for zero because we regard it as the users choice
      }
      if (ip == ipend) {
        //data may be appended
        data->vecval.enlarge_below(nel, coeff, factor);
        data->vecind.enlarge_below(nel, indices);
        data->nz_cnt += nel;
        data->maxindex = new_maxind;
        data->sparse = true;
        return 0;
      }
    }
    for (; ip != ipend; ip++) {
      if (new_maxind < *ip)
        new_maxind = *ip;
    }

    //--- data will not stay clean, add the information and clean up later
    Integer new_nz = nel;

    data->prepare_for_changes(new_maxind, new_nz);
    if (data->sparse) {
      assert((data->nz_cnt == data->vecind.dim()) && (data->nz_cnt == data->vecval.dim()));
      data->vecind.enlarge_below(Integer(nel), (const Integer*)indices);
      data->vecval.enlarge_below(Integer(nel), (const Real*)coeff, Real(factor));
      data->nz_cnt += nel;
    } else {
      Real* const valstart = data->vecval.get_store();
      while (indices != ipend) {
        valstart[*indices++] += factor * *coeff++;
      }
    }

    data->clean = false;
    data->maxindex = new_maxind;
    return 0;
  }

  // *****************************************************************************
  //                               sparsify
  // *****************************************************************************

  int Minorant::sparsify(double tol, double sparsity_ratio) {
    return data->make_clean(tol, sparsity_ratio);
  }

  // *****************************************************************************
  //                               nonzeros
  // *****************************************************************************

  int Minorant::nonzeros() {
    data->make_clean();
    return data->nz_cnt;
  }

  // *****************************************************************************
  //                             get_coeffs
  // *****************************************************************************

  int Minorant::get_coeffs(int& nel, const double*& coeffs, const int*& indices) const {
    int err = 0;
    if (!data->clean)
      err = data->make_clean();
    if (data->sparse) {
      nel = data->nz_cnt;
      coeffs = data->vecval.get_store();
      indices = data->vecind.get_store();
    } else {
      nel = data->maxindex + 1;
      coeffs = data->vecval.get_store();
      indices = 0;
    }
    return err;
  }

  // *****************************************************************************
  //                             get_coeffs
  // *****************************************************************************

  double* Minorant::get_dense_coeff_store(int nel) {
    data->normsqu = -1.;
    Integer new_maxind = max(nel - 1, data->maxindex);
    data->make_dense(new_maxind);
    data->vecind.init(0, 1, Integer(0));
    data->clean = false;
    data->maxindex = new_maxind;
    return data->vecval.get_store();
  }

  // *****************************************************************************
  //                             reassign_coeffs
  // *****************************************************************************

  int Minorant::reassign_coeffs(int nel, const int* map_to_old) {
    return data->reassign_coeffs(nel, map_to_old);
  }

  // *****************************************************************************
  //                             scale
  // *****************************************************************************

  int Minorant::scale_minorant(Real val) {
    return data->scale(val);
  }

  // *****************************************************************************
  //                             clone
  // *****************************************************************************

  Minorant* Minorant::clone_minorant(double factor, bool with_primal) const {
    return new Minorant(this, factor, with_primal);
  }

  // *****************************************************************************
  //                            get_primal()
  // *****************************************************************************

  const PrimalData* Minorant::get_primal() const {
    return data->primal;
  }

  // *****************************************************************************
  //                            get_primal()
  // *****************************************************************************

  PrimalData* Minorant::get_primal() {
    return data->primal;
  }

  // *****************************************************************************
  //                            set_primal()
  // *****************************************************************************

  void Minorant::set_primal(PrimalData* pd) {
    delete data->primal;
    data->primal = pd;
  }

  // *****************************************************************************
  //                            aggregate()
  // *****************************************************************************

  int Minorant::aggregate(const Minorant& mnrt, double factor) {
    int err0 = 0;
    data->normsqu = -1.;
    if (mnrt.get_primal()) {
      if (data->primal == 0) {
        if (data->nz_cnt != 0)
          err0 = 1;
        data->primal = mnrt.get_primal()->clone_primal_data();
        if ((factor != 1.) && (data->primal))
          data->primal->scale_primal_data(factor);
      } else {
        err0 = data->primal->aggregate_primal_data(*(mnrt.get_primal()), factor);
      }
    }
    data->n_aggregated += mnrt.number_aggregated();
    int err1 = add_offset(mnrt.offset() * factor);
    Integer n;
    const Real* cp;
    const Integer* ip;
    if (mnrt.get_coeffs(n, cp, ip))
      return 1;
    int err2 = add_coeffs(n, cp, ip, factor);
    return err0 + err1 + err2;
  }

  // *****************************************************************************
  //                            aggregate()
  // *****************************************************************************

  int Minorant::number_aggregated() const {
    return data->n_aggregated;
  }

  // *****************************************************************************
  //                            norm_squared()
  // *****************************************************************************

  double Minorant::norm_squared() const {
    return data->norm_squared();
  }


}
