/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/BoxIPBundleBlock.cxx
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


#include "BoxIPBundleBlock.hxx"

using namespace CH_Matrix_Classes;

namespace ConicBundle {

  int BoxIPBundleBlock::find_inactive_indices(Indexmatrix& inactive,
    Real trace_rhs,
    bool cautious) const {
    inactive.newsize(vecdim, 1);
    inactive.init(0, 1, Integer(0));

    Real tapia_factor = (old_mu > 0.) ? mu / old_mu : 1.;

    for (Integer i = 0; i < vecdim; i++) {
      Real lslack = x(i);
      Real uslack = -lslack;
      if (use_scaling) {
        lslack -= lb(i) * s;
        uslack += ub(i) * s;
      } else {
        lslack -= lb(i);
        uslack += ub(i);
      }
      if ((lslack > trace_rhs * (cautious ? 1e-6 : 1e-6)) && (uslack > trace_rhs * (cautious ? 1e-6 : 1e-6))) {
        //too big
        continue;
      }
      //check for relation to the dual
      if ((lz(i) < (cautious ? 1e6 : 1e3) * lslack) && (uz(i) < (cautious ? 1e6 : 1e3) * uslack)) {
        continue;
      }
      if (min(lslack, uslack) < trace_rhs * (cautious ? 1e-10 : 1e-10)) {
        //this is certainly not active
        inactive.concat_below(i);
        continue;
      }
      if ((mu < 0.01)
        && (
          ((lslack < 0.01 * std::sqrt(mu)) && (lz(i) > std::sqrt(mu))) ||
          ((uslack < 0.01 * std::sqrt(mu)) && (uz(i) > std::sqrt(mu)))
          )
        ) {
        inactive.concat_below(i);
        continue;
      }
      //check whether reduction is faster than for the dual
      if (tapia_factor < .8) {
        Real primal_tapia = lslack / (oldx(i) - (use_scaling ? s * lb(i) : lb(i)));
        Real dual_tapia = lz(i) / oldlz(i);
        if (
          (primal_tapia < 0.8) &&
          //((1.1*(primal_tapia-tapia_factor)>(dual_tapia-tapia_factor)))
          (primal_tapia < .1 * dual_tapia)
          ) {
          //goes to zero faster than the dual value
          inactive.concat_below(i);
          continue;
        }
        primal_tapia = uslack / ((use_scaling ? s * ub(i) : ub(i)) - oldx(i));
        dual_tapia = uz(i) / olduz(i);
        if (
          (primal_tapia < 0.8) &&
          //((1.1*(primal_tapia-tapia_factor)>(dual_tapia-tapia_factor)))
          (primal_tapia < .1 * dual_tapia)
          ) {
          //goes to zero faster than the dual value
          inactive.concat_below(i);
          continue;
        }
      }
    }

    if (use_scaling) {
      bool clearly_active = (
        ((s > trace_rhs * (cautious ? 1e-6 : 1e-6)) &&
          ((scalub <= 0.) || (s < scalub * (1 - 1e-6))))
        ||
        ((lt < (cautious ? 1e6 : 1e3) * s) &&
          ((scalub <= 0.) || (ut < (cautious ? 1e6 : 1e3) * s)))
        );

      if (!clearly_active) {
        if ((s < trace_rhs * (cautious ? 1e-10 : 1e-10))
          || ((scalub > 0.) && (scalub - s < scalub * (cautious ? 1e-10 : 1e-10)))
          || ((mu < 0.01) && (s < 0.01 * std::sqrt(mu)) && (lt > std::sqrt(mu)))
          || ((scalub > 0.) && (mu < 0.01) && (scalub - s < 0.01 * std::sqrt(mu)) && (ut > std::sqrt(mu)))
          ) {
          //clearly inactive
          inactive.concat_below(vecdim);
        } else {
          //check whether reduction is faster than for the dual
          Real primal_tapia = s / olds;
          Real dual_ltapia = lt / oldlt;
          if (tapia_factor < .8) {
            if ((
              (primal_tapia < 0.8) &&
              (primal_tapia < .1 * dual_ltapia)
              ) ||
              (
                (scalub > 0.) &&
                ((scalub - s) / (scalub - olds) < 0.8) &&
                ((scalub - s) / (scalub - olds) < .1 * ut / oldut)
                )
              )
              //goes to zero faster than the dual value
              inactive.concat_below(vecdim);
          }

        }
      }
    }

    if (cb_out(1)) {
      get_out() << " BoxIPBB: last_alpha=" << last_alpha;
      get_out() << " mu=" << mu << " old_mu=" << old_mu;
      get_out() << " tapia_factor=" << tapia_factor << " cautious=" << cautious << " inactive=" << transpose(inactive);
    }


    return 0;
  }

  int BoxIPBundleBlock::compute_NTscaling(void) {
    if (use_scaling) {

      lxinv.init(lb, -s);
      lxinv += x;
      lxinv.inv();

      uxinv.init(ub, s);
      uxinv -= x;
      uxinv.inv();

      hatb.init(lxinv);
      hatb %= lz;
      xiz = hatb;
      hatb %= lb;
      Cholinvb.init(uxinv);
      Cholinvb %= uz;
      xiz += Cholinvb;
      Cholinvb %= ub;
      hats = lt / s + ip(lb, hatb) + ip(ub, Cholinvb);
      if (scalub > 0.)
        hats += ut / (scalub - s);
      hatb += Cholinvb;

      sqrt_xiz.init(xiz);
      sqrt_xiz.inv();
      Cholinvb.init(hatb);
      Cholinvb %= sqrt_xiz;
      Chols = std::sqrt(hats - ip(hatb, Cholinvb));
      Cholinvb /= Chols;
      sqrt_xiz.inv();
      sqrt_xiz.sqrt();

    } else {
      lxinv.init(x);
      lxinv -= lb;
      lxinv.inv();

      uxinv.init(ub);
      uxinv -= x;
      uxinv.inv();

      xiz = lxinv;
      xiz %= lz;
      hatb = uxinv;
      hatb %= uz;
      xiz += hatb;

      sqrt_xiz.init(xiz);
      sqrt_xiz.sqrt();
    }

    return 0;
  }

  void BoxIPBundleBlock::point_changed() {
    dx.init(0, 1, 0.);
    dlz.init(0, 1, 0.);
    duz.init(0, 1, 0.);
    xiz.init(0, 1, 0.);
    compl_lrhs.init(0, 1, 0.);
    compl_urhs.init(0, 1, 0.);

    ds = 0.;
    dlt = 0.;
    dut = 0.;
    compl_ltrhs = 0.;
    compl_utrhs = 0.;
  }

  int BoxIPBundleBlock::given_dx_compute_dz(void) {
    if (use_scaling) {
      dlt = compl_ltrhs - (1. + ds / s) * lt;

      //std::cout<<" slcompl="<<((s+ds)*lt+s*dlt-s*compl_ltrhs);
      assert((s + ds) * lt + s * dlt - s * compl_ltrhs < 1e-10 * max(s, lt));

      if (scalub > 0.) {
        dut = compl_utrhs - (1. + ds / (scalub - s)) * ut;

        //std::cout<<" sucompl="<<((scalub-s-ds)*ut+(scalub-s)*dut-(scalub-s)*compl_utrhs);
        assert((scalub - s - ds) * ut + (scalub - s) * dut - (scalub - s) * compl_utrhs < 1e-10 * max(scalub - s, ut));
      }

      dlz.init(lb, ds);
      dlz -= dx;
      dlz %= lxinv;
      dlz -= 1.;
      dlz %= lz;
      dlz += compl_lrhs;
      duz.init(dx);
      duz.xpeya(ub, -ds);
      duz %= uxinv;
      duz -= 1.;
      duz %= uz;
      duz += compl_urhs;

      //std::cout<<"xlcompl="<<transpose((x+dx-(s+ds)*lb)%lz+(x-s*lb)%dlz-(x-s*lb)%compl_lrhs);
      //std::cout<<"xucompl="<<transpose(((s+ds)*ub-x-dx)%uz+(s*ub-x)%duz-(s*ub-x)%compl_urhs);

      assert(norm2((x + dx - (s + ds) * lb) % lz + (x - s * lb) % dlz - (x - s * lb) % compl_lrhs) < 1e-10 * max(norm2(x - s * lb), norm2(lz)));
      assert(norm2(((s + ds) * ub - x - dx) % uz + (s * ub - x) % duz - (s * ub - x) % compl_urhs) < 1e-10 * max(norm2(s * ub - x), norm2(uz)));
    } else {
      dlz.init(dx, -1.);
      dlz %= lxinv;
      dlz -= 1.;
      dlz %= lz;
      dlz += compl_lrhs;
      duz.init(dx);
      duz %= uxinv;
      duz -= 1.;
      duz %= uz;
      duz += compl_urhs;

      //std::cout<<"xlcompl="<<transpose((x+dx-lb)%lz+(x-lb)%dlz-(x-lb)%compl_lrhs);
      //std::cout<<"xucompl="<<transpose((ub-x-dx)%uz+(ub-x)%duz-(ub-x)%compl_urhs);

      //assert(norm2((x+dx-lb)%lz+(x-lb)%dlz-(x-lb)%compl_lrhs)<=1e-10*max(norm2(x),norm2(lz)));
      //assert(norm2((ub-x-dx)%uz+(ub-x)%duz-(ub-x)%compl_urhs)<=1e-10*max(norm2(x),norm2(uz)));
    }

    return 0;
  }

  void BoxIPBundleBlock::clear(const Matrix& inlb, const Matrix& inub, bool inscal, Real inscalub) {
    assert(inlb.rowdim() == inub.rowdim());
    assert(inlb.coldim() == 1);
    assert(inub.coldim() == 1);
    assert(min(inub - inlb) > 1e-10);


    //oracle_data=0;

    use_scaling = inscal;
    scalub = inscalub;
    if ((use_scaling == false) && (scalub > 0.)) {
      lb.init(inlb, scalub);
      ub.init(inub, scalub);
      scalub = -1.;
    } else {
      lb = inlb;
      ub = inub;
    }

    vecdim = lb.rowdim();
    scalvecdim = vecdim + (use_scaling ? 1 : 0);

    x.init(vecdim, 1, 0.);
    lz.init(vecdim, 1, 0.);
    uz.init(vecdim, 1, 0.);

    s = lt = ut = 0.;

    last_rhs_mu = 0.;
    mu = 0.;
    old_mu = 0.;
    last_alpha = 0.;

    oldx.init(vecdim, 1, 0.);
    oldlz.init(vecdim, 1, 0.);
    olduz.init(vecdim, 1, 0.);
    olds = 0;
    oldlt = 0;
    oldut = 0;

    tmpvec.newsize(vecdim, 1);
    tmpvec.init(0, 0, 0.);
    tmpmat.init(0, 0, 0.);

    //reserve space but then set back to empty
    dx.newsize(vecdim, Integer(1));
    dlz.newsize(vecdim, Integer(1));
    duz.newsize(vecdim, Integer(1));
    xiz.newsize(vecdim, Integer(1));
    compl_lrhs.newsize(vecdim, Integer(1));
    compl_urhs.newsize(vecdim, Integer(1));

    point_changed();

    bundle_dim = scalvecdim;
    diff_model.init(0, 1, 0.);
    map_to_old.init(Range(0, scalvecdim - 1));

    test_myself_call = false;
  }

  BoxIPBundleBlock::BoxIPBundleBlock(const Matrix& inlb, const Matrix& inub, bool inscal, Real inscalub, CBout* cb, int cbinc) :CBout(cb, cbinc), InteriorPointBundleBlock(cb) {
    clear(inlb, inub, inscal, inscalub);
  }

  BoxIPBundleBlock::~BoxIPBundleBlock() {
  }

  InteriorPointBundleBlock* BoxIPBundleBlock::clone() {
    BoxIPBundleBlock* p = new BoxIPBundleBlock(lb, ub, use_scaling, scalub, this, 0);
    p->copy_from(this);

    return p;
  }

  int BoxIPBundleBlock::copy_from(InteriorPointBundleBlock* inp) {
    BoxIPBundleBlock* p = dynamic_cast<BoxIPBundleBlock*>(inp);
    if (p == 0)
      return 1;

    bundle_dim = p->bundle_dim;
    diff_model = p->diff_model;
    Bt = p->Bt;
    Boffset = p->Boffset;
    sqrBnorms = p->sqrBnorms;

    //oracle_data=p->oracle_data;

    vecdim = p->vecdim;
    scalvecdim = p->scalvecdim;
    lb = p->lb;
    ub = p->ub;
    use_scaling = p->use_scaling;
    scalub = p->scalub;
    x = p->x;
    lz = p->lz;
    uz = p->uz;
    dx = p->dx;
    dlz = p->dlz;
    duz = p->duz;
    s = p->s;
    lt = p->lt;
    ut = p->ut;
    ds = p->ds;
    dlt = p->dlt;
    dut = p->dut;
    xiz = p->xiz;
    lxinv = p->lxinv;
    uxinv = p->uxinv;
    hats = p->hats;
    hatb = p->hatb;
    sqrt_xiz = p->sqrt_xiz;
    Chols = p->Chols;
    Cholinvb = p->Cholinvb;
    compl_lrhs = p->compl_lrhs;
    compl_urhs = p->compl_urhs;
    compl_ltrhs = p->compl_ltrhs;
    compl_utrhs = p->compl_utrhs;
    sys_srhs = p->sys_srhs;
    last_rhs_mu = p->last_rhs_mu;
    mu = p->mu;
    old_mu = p->old_mu;
    last_alpha = p->last_alpha;
    oldx = p->oldx;
    olds = p->olds;
    oldlz = p->oldlz;
    olduz = p->olduz;
    oldlt = p->oldlt;
    oldut = p->oldut;
    tmpvec = p->tmpvec;
    tmpmat = p->tmpmat;

    map_to_old = p->map_to_old;
    test_myself_call = p->test_myself_call;

    return 0;
  }

  Integer BoxIPBundleBlock::get_vecdim() const {
    return scalvecdim;
  }

  /// set x to value*"one" to x, or if add==true, add value*"one" to x
  int BoxIPBundleBlock::center_x(Real val, bool add) {
    point_changed();
    if (use_scaling) {
      if (add) {
        if (scalub > 0) {
          s = scalub * .5;   //never mind what others want, no interaction
        } else {
          s += (2 * vecdim + 1) * val;   //contribution to trace is now mu_dim*val
        }
      } else {
        if (scalub > 0) {
          s = scalub * .5;  //never mind what others want, no interaction
        } else {
          s = (2 * vecdim + 1) * val; //contribution to trace is now mu_dim*val
        }
      }
      x.init(lb, s * .5);
      x.xpeya(ub, s * .5);
    } else {
      x.init(lb, .5);
      x.xpeya(ub, .5);
    }

    // if (use_scaling)
    //   std::cout<<" boxs="<<s;
    // std::cout<<" boxx="<<transpose(x);

    return 0;
  }

  /// set z to value*"one" to z, or if add==true, add value*"one" to z
  int BoxIPBundleBlock::center_z(Real val, bool add) {
    point_changed();
    if (add) {
      if (use_scaling) {
        lt += val;
        ut += val;
      } else {
        lz += val;
        uz += val;
      }
    } else {
      lz.init(vecdim, 1, 1.);
      uz.init(vecdim, 1, 1.);
      if (use_scaling) {
        Real lip = ip(lb, lz);
        Real uip = ip(ub, uz);
        //choose value "a" so that   t=val/2. and t-val + a*(uip-lip) =0
        // assert(uip>lip);
        Real a = .5 * val / (uip - lip);
        lz *= a;
        uz *= a;
      }
    }
    return 0;
  }

  //this routine does not fit here really, because there need not be a shift that makes x feasible
  int BoxIPBundleBlock::set_x(const Matrix& vec, Integer startindex, Real& add_center_value) {
    //missing: scaling, adding center_values: what is the center here
    assert(vec.dim() >= startindex + scalvecdim);
    point_changed();
    Real maxlbviol = max_Real;
    Real maxubviol = max_Real;
    const Real* vp = vec.get_store() + startindex;
    Real* xp = x.get_store();
    Real* lbp = lb.get_store();
    Real* ubp = ub.get_store();
    const Real* const xend = xp + vecdim;
    while (xp != xend) {
      Real d = *vp++;
      (*xp++) = d;
      Real dlb = d - (*lbp++);
      if (dlb < maxlbviol)
        maxlbviol = dlb;
      Real dub = (*ubp++) - d;
      if (dub < maxubviol)
        maxubviol = dub;
    }
    s = *vp;
    if (s < maxlbviol)
      maxlbviol = s;
    if ((scalub > 0.) && (scalub - s < maxubviol))
      maxubviol = scalub - s;
    if (maxlbviol + maxubviol <= 2e-6)
      return 1;
    if ((maxlbviol <= 0) || (maxubviol <= 0))
      add_center_value = max(-maxlbviol, -maxubviol);
    else
      add_center_value = 0;

    return 0;
  }

  int BoxIPBundleBlock::set_z(const Matrix& vec, Integer startindex, Real& add_center_value) {
    assert(vec.dim() >= startindex + vecdim);
    point_changed();
    const Real* vp = vec.get_store() + startindex;
    Real* lzp = lz.get_store();
    Real* uzp = uz.get_store();
    const Real* const zend = lzp + vecdim;
    while (lzp != zend) {
      Real d = *vp++;
      if (d < 0) {
        *uzp++ = -d + 1.;
        *lzp++ = 1.;
      } else {
        *uzp++ = 1.;
        *lzp++ = d + 1.;
      }
    }
    lt = 0.;
    if (use_scaling) {
      lt = ip(lb, lz) - ip(ub, uz) + *vp;
      if (scalub > 0.) {
        ut = max(-lt + 1., 1.);
        lt += ut;
      }
    }
    add_center_value = max(0., -lt);

    // // TEST output begin
    // if (use_scaling){
    //   std::cout<<" lt="<<lt;
    //   if (scalub>0)
    // 	std::cout<<" ut="<<ut;
    // }
    // std::cout<<" lz="<<transpose(lz)<<" uz="<<transpose(uz);
    // // TEST output end

    return 0;
  }

  /// set z to the slack of the bundle and return a value>=0 that needs to be added to make it feasible
  int BoxIPBundleBlock::set_bundle_z(const Matrix& y,
    MinorantBundle& globalbundle,
    Integer startindex_bundle,
    Real& add_center_value) {
    diff_model.newsize(scalvecdim, 1); chk_set_init(diff_model, 1);
    for (Integer i = 0; i < scalvecdim; i++) {
      diff_model(i) = -globalbundle[unsigned(startindex_bundle + map_to_old(i))].evaluate(-1, y);
    }
    // std::cout<<" BoxIPdiffmodel="<<diff_model;
    return set_z(diff_model, 0, add_center_value);
  }


  /// on vec[startindex+0,+1 ...,+(vecdim-1)] put or add  a * x into vec for a real number a   
  int BoxIPBundleBlock::vecgetsax(Matrix& vec, Integer startindex, Real a, bool add) {
    assert(vec.dim() >= startindex + vecdim);
    if (add == false) {
      mat_xeya(vecdim, vec.get_store() + startindex, x.get_store(), a);
    } else {
      mat_xpeya(vecdim, vec.get_store() + startindex, x.get_store(), a);
    }
    if (use_scaling) {
      if (add == false)
        vec(startindex + vecdim) = a * s;
      else
        vec(startindex + vecdim) += a * s;
    }
    return 0;
  }

  /// on vec[startindex+0,+1 ...,+(vecdim-1)]   put or add a * z into vec for a real number a   
  int BoxIPBundleBlock::vecgetsaz(Matrix& vec, Integer startindex, Real a, bool add) {
    assert(vec.dim() >= startindex + vecdim);
    if (add == false) {
      mat_xeya(vecdim, vec.get_store() + startindex, lz.get_store(), a);
      mat_xpeya(vecdim, vec.get_store() + startindex, uz.get_store(), -a);
    } else {
      mat_xpeya(vecdim, vec.get_store() + startindex, lz.get_store(), a);
      mat_xpeya(vecdim, vec.get_store() + startindex, uz.get_store(), -a);
    }
    if (use_scaling) {
      if (add == false)
        vec(startindex + vecdim) = a * (lt - (scalub > 0 ? ut : 0.));
      else
        vec(startindex + vecdim) += a * (lt - (scalub > 0 ? ut : 0.));
    }
    return 0;
  }

  /// add dimensions of the primal-dual pairs to mudim and add the inner products of the primal-dual pairs of the current point to current_ip and those of the next point obtained by the given stepsize with the latest computed step to step_ip
  int BoxIPBundleBlock::get_mu_info(Integer& inmudim,
    Real& tr_xz,
    Real& tr_xdzpdxz,
    Real& tr_dxdz,
    Real& min_xz,
    Real& max_xz) const {
    assert((dx.dim() == vecdim) && (dx.coldim() == 1));

    inmudim += vecdim + scalvecdim + (scalub > 0. ? 1 : 0);
    Real minv, maxv;

    if (use_scaling) {
      //ip(x-s*lb,lz)
      tmpvec.init(lb, -s);
      tmpvec += x;
      tr_xz += ip_min_max(tmpvec, lz, minv, maxv);
      tr_xdzpdxz += ip(tmpvec, dlz);
      tmpvec.init(lb, -ds);
      tmpvec += dx;
      tr_xdzpdxz += ip(tmpvec, lz);
      tr_dxdz += ip(tmpvec, dlz);

      if (minv < min_xz)
        min_xz = minv;
      if (maxv > max_xz)
        max_xz = maxv;

      if (cb_out(2))
        get_out() << " xlz[" << minv << "," << maxv << "]";

      //ip(s*ub-x,uz)
      tmpvec.init(ub, s);
      tmpvec -= x;
      tr_xz += ip_min_max(tmpvec, uz, minv, maxv);
      tr_xdzpdxz += ip(tmpvec, duz);
      tmpvec.init(ub, ds);
      tmpvec -= dx;
      tr_xdzpdxz += ip(tmpvec, uz);
      tr_dxdz += ip(tmpvec, duz);

      if (minv < min_xz)
        min_xz = minv;
      if (maxv > max_xz)
        max_xz = maxv;

      if (cb_out(2))
        get_out() << " xuz[" << minv << "," << maxv << "]";

      //ip(s,lt)
      Real d = s * lt;
      tr_xz += d;
      tr_xdzpdxz += s * dlt + ds * lt;
      tr_dxdz += ds * dlt;

      if (d < min_xz)
        min_xz = d;
      else if (d > max_xz)
        max_xz = d;

      if (cb_out(2))
        get_out() << " slt[" << d << "]";

      if (scalub > 0.) {
        Real d = (scalub - s) * ut;
        tr_xz += d;
        tr_xdzpdxz += (scalub - s) * dut - ds * ut;
        tr_dxdz += -ds * dut;

        if (d < min_xz)
          min_xz = d;
        else if (d > max_xz)
          max_xz = d;

        if (cb_out(2))
          get_out() << " sut[" << d << "]";
      }

    } else {
      //ip(x-lb,lz)
      tmpvec.init(x);
      tmpvec -= lb;
      tr_xz += ip_min_max(tmpvec, lz, minv, maxv);
      tr_xdzpdxz += ip(tmpvec, dlz) + ip(dx, dlz);
      tr_dxdz += ip(dx, dlz);

      if (minv < min_xz)
        min_xz = minv;
      if (maxv > max_xz)
        max_xz = maxv;

      if (cb_out(2))
        get_out() << " xlz[" << minv << "," << maxv << "]";

      //ip(ub-x,uz)
      tmpvec.init(ub);
      tmpvec -= x;
      tr_xz += ip_min_max(tmpvec, uz, minv, maxv);
      tr_xdzpdxz += ip(tmpvec, duz) - ip(dx, duz);
      tr_dxdz += -ip(dx, duz);

      if (minv < min_xz)
        min_xz = minv;
      if (maxv > max_xz)
        max_xz = maxv;

      if (cb_out(2))
        get_out() << " xuz[" << minv << "," << maxv << "]";

    }

    return 0;
  }

  /// for limiting the stepsize with respect to the neighborhood this information about x_i*z_i and dx_i*dz_i is required, each block *adds* its contribution to the numbers
  int BoxIPBundleBlock::get_nbh_info(Integer inmudim,
    Real tr_xz,
    Real tr_xdzpdxz,
    Real tr_dxdz,
    Real nbh_ubnd,
    Real& alpha,
    Real& max_nbh,
    Real& nrmsqr_xz,
    Real& nrmsqr_xdzpdxz,
    Real& nrmsqr_dxdz,
    Real& ip_xz_xdzpdxz,
    Real& ip_xz_dxdz,
    Real& ip_dxdz_xdzpdxz) const {
    assert(use_scaling || (ds == 0.));
    const Real mu_xz = tr_xz / inmudim;
    const Real mu_xdzpdxz = tr_xdzpdxz / inmudim;
    const Real mu_dxdz = tr_dxdz / inmudim;
    const Real mu_at_one = mu_xz + mu_xdzpdxz + mu_dxdz;

    if (use_scaling) {
      NNC_nbh_stepsize(s, lt, ds, dlt,
        mu_xz, mu_xdzpdxz, mu_dxdz, mu_at_one,
        nbh_ubnd, alpha, max_nbh,
        nrmsqr_xz, nrmsqr_xdzpdxz, nrmsqr_dxdz,
        ip_xz_xdzpdxz, ip_xz_dxdz, ip_dxdz_xdzpdxz);

      if (scalub > 0.) {
        NNC_nbh_stepsize(scalub - s, ut, -ds, dut,
          mu_xz, mu_xdzpdxz, mu_dxdz, mu_at_one,
          nbh_ubnd, alpha, max_nbh,
          nrmsqr_xz, nrmsqr_xdzpdxz, nrmsqr_dxdz,
          ip_xz_xdzpdxz, ip_xz_dxdz, ip_dxdz_xdzpdxz);
      }
    }

    const Real* xinvp = lxinv.get_store();
    const Real* zp = lz.get_store();
    const Real* dxp = dx.get_store();
    const Real* dzp = dlz.get_store();
    for (Integer i = 0; i < lb.rowdim(); i++) {
      NNC_nbh_stepsize(1. / (*xinvp++), *zp++, *dxp++ - ds * lb(i), *dzp++,
        mu_xz, mu_xdzpdxz, mu_dxdz, mu_at_one,
        nbh_ubnd, alpha, max_nbh,
        nrmsqr_xz, nrmsqr_xdzpdxz, nrmsqr_dxdz,
        ip_xz_xdzpdxz, ip_xz_dxdz, ip_dxdz_xdzpdxz);
    }
    xinvp = uxinv.get_store();
    zp = uz.get_store();
    dxp = dx.get_store();
    dzp = duz.get_store();
    for (Integer i = 0; i < ub.rowdim(); i++) {
      NNC_nbh_stepsize(1. / (*xinvp++), *zp++, ds * ub(i) - *dxp++, *dzp++,
        mu_xz, mu_xdzpdxz, mu_dxdz, mu_at_one,
        nbh_ubnd, alpha, max_nbh,
        nrmsqr_xz, nrmsqr_xdzpdxz, nrmsqr_dxdz,
        ip_xz_xdzpdxz, ip_xz_dxdz, ip_dxdz_xdzpdxz);
    }

    return 0;
  }


  /// if necessary, reduce alpha to the biggest value so that feasibility is maintained with this step size
  int BoxIPBundleBlock::linesearch(Real& alpha) const {
    assert((dx.dim() == vecdim) && (dx.coldim() == 1));

    const Real* vp = x.get_store();
    const Real* lp = lb.get_store();
    const Real* up = ub.get_store();
    const Real* dvp = dx.get_store();
    const Real* const vendx = vp + vecdim;
    Real a = alpha;
    if (use_scaling) {
      while (vp != vendx) {
        Real d = *dvp++;
        Real l = *lp++;
        Real u = *up++;
        Real valx = *vp++;
        if (d < ds * l) {
          Real ax = (valx - s * l) / (ds * l - d);
          if (ax < a)
            a = ax;
        }
        if (d > ds * u) {
          Real ax = (s * u - valx) / (d - ds * u);
          if (ax < a)
            a = ax;
        }
      }
    } else {
      while (vp != vendx) {
        Real d = *dvp++;
        Real l = *lp++;
        Real u = *up++;
        Real valx = *vp++;
        if (d < 0.) {
          d = -(valx - l) / d;
          if (d < a)
            a = d;
        } else if (d > 0.) {
          d = (u - valx) / d;
          if (d < a)
            a = d;
        }
      }
    }
    vp = lz.get_store();
    dvp = dlz.get_store();
    const Real* const vendlz = vp + vecdim;
    while (vp != vendlz) {
      Real d = *dvp++;
      if (d < 0.) {
        d = -(*vp++) / d;
        if (d < a)
          a = d;
      } else {
        vp++;
      }
    }
    vp = uz.get_store();
    dvp = duz.get_store();
    const Real* const venduz = vp + vecdim;
    while (vp != venduz) {
      Real d = *dvp++;
      if (d < 0.) {
        d = -(*vp++) / d;
        if (d < a)
          a = d;
      } else {
        vp++;
      }
    }
    if (use_scaling) {
      if (ds < 0.) {
        Real d = -s / ds;
        if (d < a)
          a = d;
      }
      if (dlt < 0.) {
        Real d = -lt / dlt;
        if (d < a)
          a = d;
      }
      if (scalub > 0.) {
        if (ds > 0) {
          Real d = (scalub - s) / ds;
          if (d < a)
            a = d;
        }
        if (dut < 0.) {
          Real d = -ut / dut;
          if (d < a)
            a = d;
        }
      }
    }
    if (a < alpha)
      alpha = a;
    return 0;
  }

  /// compute the complementarity_rhs=rhsmu*xi-rhscorr*xi*dx*dz (wihtout "-z") for mu=rhsmu and for corrector for factor rhscorr>0., store this and add it to rhs 
  int BoxIPBundleBlock::add_muxinv(Matrix& rhs,
    Integer startindex,
    Real rhsmu,
    Real rhscorr,
    bool minus) {
    assert(rhs.dim() >= startindex + vecdim);
    assert(rhscorr >= 0.);
    assert(rhsmu >= 0.);

    //make sure, the system is available
    if (xiz.dim() != vecdim) {
      compute_NTscaling();
    }

    last_rhs_mu = rhsmu;

    if (rhsmu > 0.) {
      compl_lrhs.init(lxinv, rhsmu);
      compl_urhs.init(uxinv, rhsmu);
      if (use_scaling) {
        compl_ltrhs = rhsmu / s;
        sys_srhs = compl_ltrhs + rhsmu * (-ip(lb, lxinv) + ip(ub, uxinv));
        if (scalub > 0.) {
          compl_utrhs = rhsmu / (scalub - s);
          sys_srhs -= compl_utrhs;
        }
      }
    } else {
      compl_lrhs.init(vecdim, 1, 0.);
      compl_urhs.init(vecdim, 1, 0.);
      compl_ltrhs = 0;
      compl_utrhs = 0;
      sys_srhs = 0;
    }

    if (rhscorr > 0.) {
      assert((dx.dim() == vecdim) && (dx.coldim() == 1));
      if (use_scaling) {
        tmpvec.init(lb, -ds);
        tmpvec += dx;
        tmpvec %= dlz;
        tmpvec %= lxinv;
        compl_lrhs.xpeya(tmpvec, -rhscorr);
        sys_srhs += rhscorr * ip(tmpvec, lb);

        tmpvec.init(ub, ds);
        tmpvec -= dx;
        tmpvec %= duz;
        tmpvec %= uxinv;
        compl_urhs.xpeya(tmpvec, -rhscorr);
        sys_srhs -= rhscorr * ip(tmpvec, ub);

        compl_ltrhs -= rhscorr * ds * dlt / s;
        sys_srhs -= rhscorr * ds * dlt / s;
        if (scalub > 0.) {
          compl_utrhs += rhscorr * ds * dut / (scalub - s);
          sys_srhs -= rhscorr * ds * dut / (scalub - s);
        }
      } else {
        tmpvec.init(dx);
        tmpvec %= dlz;
        tmpvec %= lxinv;
        compl_lrhs.xpeya(tmpvec, -rhscorr);

        tmpvec.init(dx);
        tmpvec %= duz;
        tmpvec %= uxinv;
        compl_urhs.xpeya(tmpvec, rhscorr);
      }
    }

    // std::cout<<"startind="<<startindex<<std::endl;
    // std::cout<<" mu="<<rhsmu<<std::endl;
    // std::cout<<" rhscorr="<<rhscorr<<std::endl;
    // std::cout<<" compl_rhs="<<compl_rhs;
    // std::cout<<" rhsbefore="<<rhs;

    if ((rhsmu > 0) || (rhscorr > 0)) {
      mat_xpeya(vecdim, rhs.get_store() + startindex, compl_lrhs.get_store(), minus ? -1 : +1.);
      mat_xpeya(vecdim, rhs.get_store() + startindex, compl_urhs.get_store(), minus ? +1 : -1.);
      if (use_scaling) {
        rhs(startindex + vecdim) += (minus ? -sys_srhs : sys_srhs);
      }
    }

    // std::cout<<" rhsafter="<<rhs;

    return 0;
  }

  /// extract dx from rhs at startindex and compute at the same time dz (=-sys dx -z +complentarity_rhs)
  int BoxIPBundleBlock::set_dx(const CH_Matrix_Classes::Matrix& rhs, CH_Matrix_Classes::Integer startindex) {
    assert((xiz.dim() == vecdim) && (xiz.coldim() == 1));
    assert(rhs.dim() >= startindex + vecdim);
    assert((compl_lrhs.dim() == vecdim) && (compl_lrhs.coldim() == 1));

    dx.init(vecdim, 1, rhs.get_store() + startindex);
    if (use_scaling)
      ds = rhs(startindex + vecdim);

    given_dx_compute_dz();

    // dz.init(dx,-1.);
    // dz%=xiz;
    // dz+=compl_rhs;
    // dz-=z;

    // std::cout<<" setdx="<<dx;
    // std::cout<<" setdz="<<dz;
    // std::cout<<" setx="<<x;
    // std::cout<<" setz="<<z;
    // std::cout<<" setcompl_rhs="<<compl_rhs;
    // std::cout.flush();

    return 0;
  }

  /// compute dx=sysinv*rhs and at the same time dz (=-rhs-z +complentarity_rhs); 
  int BoxIPBundleBlock::set_dx_xizsolverhs(const Matrix& rhs, Integer startindex) {
    assert((xiz.dim() == vecdim) && (xiz.coldim() == 1));
    assert(rhs.dim() >= startindex + vecdim);
    assert((compl_lrhs.dim() == vecdim) && (compl_lrhs.coldim() == 1));
    assert((compl_urhs.dim() == vecdim) && (compl_urhs.coldim() == 1));

    dx.newsize(vecdim, 1); chk_set_init(dx, 1);
    Real* dxp = dx.get_store();
    const Real* const xpend = dxp + vecdim;
    const Real* vp = rhs.get_store() + startindex;
    const Real* xizp = xiz.get_store();
    while (dxp != xpend) {
      // assert(*xizp>0);
      (*dxp++) = (*vp++) / (*xizp++);
    }

    if (use_scaling) {
      //finish computation of ds and dx
      ds = ((*vp) + ip(dx, hatb)) / Chols / Chols;
      dx.xpeya(Cholinvb, ds * Chols);

      // //TEST begin
      // if (!test_myself_call){
      // 	test_myself_call=true;
      // 	Matrix testvec(dx);
      // 	testvec.concat_below(ds);
      // 	Real cmpnorm=mat_ip(scalvecdim,rhs.get_store());
      // 	apply_xiz(testvec,0);
      // 	mat_xmey(scalvecdim,testvec.get_store(),rhs.get_store()+startindex);
      // 	assert(norm2(testvec)<1e-10*(std::sqrt(cmpnorm)+1.));
      // 	test_myself_call=false;
      // }
      // //TEST end
    }

    given_dx_compute_dz();

    return 0;
  }

  /// compute sysinv*rhs into rhs 
  int BoxIPBundleBlock::apply_xizinv(Matrix& rhs, Integer startindex, bool minus) {
    assert(rhs.dim() >= startindex + vecdim);
    //make sure, the system is available
    if (xiz.dim() != vecdim) {
      compute_NTscaling();
    }

    // Matrix testvec(scalvecdim,1,rhs.get_store()+startindex); // TEST

    Real* vp = rhs.get_store() + startindex;
    const Real* const vend = vp + vecdim;
    const Real* xizp = xiz.get_store();
    if (minus) {
      while (vp != vend) {
        // assert(*xizp>0);
        (*vp++) /= -(*xizp++);
      }
    } else {
      while (vp != vend) {
        // assert(*xizp>0);
        (*vp++) /= (*xizp++);
      }
    }

    if (use_scaling) {
      if (minus) {
        (*vp) = -(*vp - mat_ip(vecdim, hatb.get_store(), rhs.get_store() + startindex)) / Chols / Chols;
      } else {
        (*vp) = (*vp + mat_ip(vecdim, hatb.get_store(), rhs.get_store() + startindex)) / Chols / Chols;
      }
      mat_xpeya(vecdim, rhs.get_store() + startindex, Cholinvb.get_store(), (*vp) * Chols);
    }

    // //TEST begin
    // if (!test_myself_call){
    //   test_myself_call=true;
    //   //std::cout<<" minus="<<minus<<" xiz="<<transpose(xiz);
    //   //std::cout<<" hats="<<hats<<" hatb="<<transpose(hatb);
    //   //std::cout<<" xiz="<<transpose(xiz);
    //   //std::cout<<" Chols="<<Chols<<" Cholinvb="<<transpose(Cholinvb);
    //   //std::cout<<" tv="<<transpose(testvec);
    //   Matrix t2vec(scalvecdim,1,rhs.get_store()+startindex);      
    //   //std::cout<<" xizinv*tv=tv2="<<transpose(t2vec);
    //   apply_xiz(t2vec,0,minus);
    //   //std::cout<<" xiz*tv2="<<transpose(t2vec);
    //   assert(norm2(t2vec-testvec)<=1e-10*(norm2(testvec)+1.));
    //   test_myself_call=false;
    // }
    // //TEST end


    return 0;
  }

  /// compute sys*rhs into rhs 
  int BoxIPBundleBlock::apply_xiz(Matrix& rhs, Integer startindex, bool minus) {
    assert(rhs.dim() >= startindex + vecdim);

    //make sure, the system is available
    if (xiz.dim() != vecdim) {
      compute_NTscaling();
    }

    // Matrix testvec(scalvecdim,1,rhs.get_store()+startindex); // TEST

    Real hatbip = 0;
    if (use_scaling) {
      hatbip = mat_ip(vecdim, rhs.get_store() + startindex, hatb.get_store());
    }

    Real* vp = rhs.get_store() + startindex;
    const Real* const vend = vp + vecdim;
    const Real* xizp = xiz.get_store();
    if (minus) {
      while (vp != vend) {
        (*vp++) *= -(*xizp++);
      }
    } else {
      while (vp != vend) {
        (*vp++) *= (*xizp++);
      }
    }

    if (use_scaling) {
      mat_xpeya(vecdim, rhs.get_store() + startindex, hatb.get_store(), minus ? *vp : -*vp);
      *vp = (*vp) * hats - hatbip;
      if (minus)
        *vp *= -1.;
    }

    // //TEST begin
    // if (!test_myself_call){
    //   test_myself_call=true;
    //   Matrix t2vec(scalvecdim,1,rhs.get_store()+startindex);
    //   apply_xizinv(t2vec,0,minus);
    //   assert(norm2(t2vec-testvec)<=1e-10*(norm2(testvec)+1.));
    //   test_myself_call=false;
    // }
    // //TEST end

    return 0;
  }

  /// move to (x+alpha*dx, z+alpha*dz)
  int BoxIPBundleBlock::do_step(Real alpha) {
    assert((dx.dim() == vecdim) && (dx.coldim() == 1));

    if ((old_mu == 0.) || (last_rhs_mu < mu)) {
      oldx = x;
      oldlz = lz;
      olduz = uz;
      old_mu = mu;
      last_alpha = alpha;

      if (use_scaling) {
        olds = s;
        oldlt = lt;
        if (scalub > 0.) {
          oldut = ut;
        }
      }
    }

    mu = last_rhs_mu;

    x.xpeya(dx, alpha);
    lz.xpeya(dlz, alpha);
    uz.xpeya(duz, alpha);

    if (use_scaling) {
      s += alpha * ds;
      lt += alpha * dlt;
      if (scalub > 0.) {
        ut += alpha * ut;
      }
    }

    // //TEST output begin
    // std::cout<<" step: alpha="<<alpha;				
    // if (use_scaling){
    //   std::cout<<" ds="<<ds<<" s="<<s;
    //   std::cout<<" dlt="<<dlt<<" lt="<<lt;
    //   if (scalub>0.){
    // 	std::cout<<" dut="<<dut<<" ut="<<ut;
    //   }
    // } 
    // std::cout<<std::endl;
    // std::cout<<" dx="<<transpose(dx);
    // std::cout<<" x="<<transpose(x);
    // std::cout<<" dlz="<<transpose(dlz);
    // std::cout<<" lz="<<transpose(lz);
    // std::cout<<" duz="<<transpose(duz);
    // std::cout<<" uz="<<transpose(uz);
    // //TEST output end

    point_changed();

    return 0;
  }

  /// move to (x+alpha*dx, z+alpha*dz), update diff_model and possibly reduce the model size if some part is below zero_threshold
  int BoxIPBundleBlock::do_bundle_step(Real alpha,
    const Matrix& y,
    MinorantBundle& globalbundle,
    Integer startindex_bundle,
    Real tracedual,
    Real /* trace_rhs */) {
    do_step(alpha);

    // // check for size reductions due to inactivity (this never worked well)
    // Indexmatrix delind;
    // find_inactive_indices(delind,trace_rhs,true);
    // if (delind.rowdim()>0){
    //   //would need to be adapted to switch to another mode if the scaling variable is deleted
    //   if (cb_out(2))
    //     get_out()<<" delete "<<delind<<std::endl;
    //   //std::cout<<" delete "<<delind<<std::endl;
    //   lb.delete_rows(delind,true);
    //   ub.delete_rows(delind,true);
    //   x.delete_rows(delind,true);
    //   lz.delete_rows(delind,true);
    //   uz.delete_rows(delind,true);
    //   oldx.delete_rows(delind,true);
    //   oldlz.delete_rows(delind,true);
    //   olduz.delete_rows(delind,true);
    //   map_to_old.delete_rows(delind,true);
    //   vecdim=x.rowdim();
    // }

    diff_model.newsize(scalvecdim, 1); chk_set_init(diff_model, 1);
    for (Integer i = 0; i < scalvecdim; i++) {
      diff_model(i) = -globalbundle[unsigned(startindex_bundle + map_to_old(i))].evaluate(-1, y);
    }
    if ((use_scaling) && (scalub <= 0.))
      diff_model(vecdim) += tracedual;
    return 0;
  }



  /// add the Schur complement to a big system matrix
  int BoxIPBundleBlock::add_AxizinvAt(const Matrix& A,
    Symmatrix& globalsys,
    bool minus,
    bool Atrans) {
    assert((Atrans == false ? A.coldim() : A.rowdim()) == vecdim);
    assert((Atrans == false ? A.rowdim() : A.coldim()) == globalsys.rowdim());

    if (xiz.dim() != vecdim) {
      compute_NTscaling();

    }

    //take the columns of A * cholL^{-T} and use rankadd for each

    for (Integer i = 0; i < vecdim; i++) {
      if (Atrans == false)
        tmpvec.init(A.rowdim(), 1, A.get_store() + i * A.rowdim(), 1, 1. / sqrt_xiz(i));
      else
        tmpvec.init(A.rowdim(), 1, A.get_store() + i, A.rowdim(), 1. / sqrt_xiz(i));
      rankadd(tmpvec, globalsys, minus ? -1. : 1., 1.);
    }

    if (use_scaling) {
      tmpvec.newsize(scalvecdim, 1); chk_set_init(tmpvec, 1);
      mat_xey(vecdim, tmpvec.get_store(), Cholinvb.get_store());
      tmpvec(vecdim) = 1. / Chols;
      genmult(A, tmpvec, tmpmat, 1., 0., Atrans);
      rankadd(tmpvec, globalsys, minus ? -1. : 1., 1.);
    }

    return 0;
  }

  ///add the "scaled" minorant outer products to globalsys, where the correct minroants start at the given index
  int BoxIPBundleBlock::add_BtinvsysB(Symmatrix& globalsys,
    const MinorantBundle& globalbundle,
    Integer startindex_bundle) {
    if (Bt.coldim() != scalvecdim) {
      Bt.newsize(globalsys.rowdim(), scalvecdim); chk_set_init(tmpmat, 1);
      for (Integer i = 0; i < scalvecdim; i++) {
        Real dummy;
        globalbundle[unsigned(startindex_bundle + map_to_old(i))].get_minorant(dummy, Bt, i);
      }
    }
    return add_AxizinvAt(Bt, globalsys);
  }


  /// add the system matrix*factor into a big system matrix starting at startindex
  int BoxIPBundleBlock::add_xiz(Symmatrix& globalsys, Integer startindex, bool minus) {
    if (xiz.dim() != vecdim) {
      compute_NTscaling();
    }

    // // TEST output start
    // std::cout<<"x="<<x;
    // std::cout<<"z="<<z;
    // std::cout<<"xiz="<<xiz;
    // // TEST output end

    const Real* xizp = xiz.get_store();
    if (minus) {
      for (Integer i = startindex; i < startindex + vecdim; i++) {
        globalsys(i, i) -= *xizp++;
      }
    } else {
      for (Integer i = startindex; i < startindex + vecdim; i++) {
        globalsys(i, i) += *xizp++;
      }
    }

    if (use_scaling) {
      Integer colind = startindex + vecdim;
      const Real* bp = hatb.get_store();
      if (minus) {
        for (Integer i = startindex; i < startindex + vecdim; i++) {
          globalsys(i, colind) += *bp++;
        }
        globalsys(colind, colind) -= hats;
      } else {
        for (Integer i = startindex; i < startindex + vecdim; i++) {
          globalsys(i, colind) -= *bp++;
        }
        globalsys(colind, colind) += hats;
      }
    }

    return 0;
  }


  Real BoxIPBundleBlock::evaluate_trace_x() {
    if ((!use_scaling) || (scalub > 0.))
      return 0.;
    return s;
  }

  Real BoxIPBundleBlock::evaluate_trace_z() {
    return sum(lz) + sum(uz) + lt + ut;
  }

  Real BoxIPBundleBlock::evaluate_trace_dx() {
    if ((!use_scaling) || (scalub > 0.))
      return 0.;

    assert(xiz.dim() == vecdim);
    assert(dx.dim() == vecdim);
    return ds;
  }

  Matrix& BoxIPBundleBlock::B_times(const Matrix& A,
    Matrix& C,
    Real alpha,
    Real beta,
    int Btrans,
    int Atrans,
    Integer startindex_model,
    MinorantBundle& globalbundle,
    Integer startindex_bundle) {
    if ((startindex_bundle == 0) && (startindex_model == 0) && (beta != 1.)) {
      if (beta == 0.)
        C.init(C.rowdim(), C.coldim(), 0.);
      else
        C *= beta;
    }
    for (Integer i = 0; i < scalvecdim; i++) {
      globalbundle[unsigned(startindex_bundle + map_to_old(i))].left_genmult(A, C, alpha, 1., (Btrans == 0), Atrans, startindex_model++);
    }
    return C;
  }

  /// C=beta*C+alpha*A*B where B and B may be transposed; carry out the model part of this beginning at startindex_model 
  Matrix& BoxIPBundleBlock::times_B(const Matrix& A,
    Matrix& C,
    Real alpha,
    Real beta,
    int Atrans,
    int Btrans,
    Integer startindex_model,
    MinorantBundle& globalbundle,
    Integer startindex_bundle) {
    if ((startindex_bundle == 0) && (startindex_model == 0) && (beta != 1.)) {
      if (beta == 0.)
        C.init(C.rowdim(), C.coldim(), 0.);
      else
        C *= beta;
    }
    for (Integer i = 0; i < scalvecdim; i++) {
      globalbundle[unsigned(startindex_bundle + map_to_old(i))].right_genmult(A, C, alpha, 1., Atrans, (Btrans == 0), startindex_model++);
    }
    return C;
  }

  ///add the main diagonal block tranpose(projection)*diagvec*projection to bigS starting at startindex
  Symmatrix& BoxIPBundleBlock::add_BDBt(const Matrix& diagvec,
    Symmatrix& bigS,
    bool minus,
    Integer startindex,
    Matrix& Bt,
    Integer startindex_model,
    MinorantBundle& globalbundle,
    Integer startindex_bundle) {
    if (minus) {
      for (Integer i = 0; i < scalvecdim; i++) {
        MinorantPointer& mp = globalbundle[unsigned(startindex_bundle + map_to_old(i))];
        for (Integer j = i + startindex_model; j < Bt.coldim(); j++) {
          bigS(startindex + startindex_model + i, startindex + j) -= mp.ip(Bt, &diagvec, j * Bt.rowdim());
        }
      }
    } else {
      for (Integer i = 0; i < scalvecdim; i++) {
        MinorantPointer& mp = globalbundle[unsigned(startindex_bundle + map_to_old(i))];
        for (Integer j = i + startindex_model; j < Bt.coldim(); j++) {
          bigS(startindex + startindex_model + i, startindex + j) += mp.ip(Bt, &diagvec, j * Bt.rowdim());
        }
      }
    }
    return bigS;
  }


  /// get the current matrix for the coupling matrix Bt in the first row of blocks
  Matrix& BoxIPBundleBlock::get_Bt(Matrix& globBt,
    Integer startindex_model,
    MinorantBundle& globalbundle,
    Integer startindex_bundle) {
    assert((&Bt != &globBt) || (startindex_model == 0));
    if (Bt.coldim() != scalvecdim) {
      Bt.newsize(globBt.rowdim(), scalvecdim); chk_set_init(Bt, 1);
      Boffset.newsize(scalvecdim, 1); chk_set_init(Boffset, 1);
      for (Integer i = 0; i < scalvecdim; i++) {
        globalbundle[unsigned(startindex_bundle + map_to_old(i))].get_minorant(Boffset(i), Bt, i);
      }
    }
    if (&Bt != &globBt) {
      assert(Bt.rowdim() == globBt.rowdim());
      assert(globBt.coldim() >= startindex_model + scalvecdim);
      mat_xey(scalvecdim * Bt.rowdim(), globBt.get_store() + startindex_model * globBt.rowdim(), Bt.get_store());
    }
    return globBt;
  }

  /// adds opB transposed times modelx (without constant affine term) to the arguments
  int BoxIPBundleBlock::add_modelx_aggregate(Real& val,
    Matrix& vec,
    MinorantBundle& globalbundle,
    Integer startindex_bundle) {
    for (Integer i = 0; i < vecdim; i++) {
      globalbundle[unsigned(startindex_bundle + map_to_old(i))].get_minorant(val, vec, 0, x(i), true);
    }
    if (use_scaling)
      globalbundle[unsigned(startindex_bundle + map_to_old(vecdim))].get_minorant(val, vec, 0, s, true);
    return 0;
  }

  /// set the model violation for the current system solution 
  int BoxIPBundleBlock::get_sysviol_model(Matrix& sysviol_model,
    Integer startindex_model,
    const Matrix& dy,
    const Real deltatrdual,
    MinorantBundle& globalbundle,
    Integer startindex_bundle) {
    for (Integer i = 0; i < vecdim; i++) {
      sysviol_model(startindex_model + i) = globalbundle[unsigned(startindex_bundle + map_to_old(i))].evaluate(-1, dy, false) - diff_model(i) + lz(i) + dlz(i) - uz(i) - duz(i);
    }
    if (use_scaling) {
      sysviol_model(startindex_model + vecdim) = globalbundle[unsigned(startindex_bundle + map_to_old(vecdim))].evaluate(-1, dy, false) - diff_model(vecdim) + lt + dlt - ip(lb, lz) - ip(lb, dlz) + ip(ub, uz) + ip(ub, duz) - (scalub > 0. ? ut + dut : deltatrdual);
    }
    return 0;
  }

  int BoxIPBundleBlock::add_trace(Matrix& vec, Real trace_dual, Integer startindex) {
    if ((use_scaling) && (scalub <= 0.))
      vec(startindex + vecdim) += trace_dual;
    return 0;
  }

  /// set the trace premultiplied by sqrt(inv(xiz)) in vec[startindex+0,...,startindex+dim_bundle()-1]
  int BoxIPBundleBlock::set_xizinvsqrt_trace(Matrix& vec,
    Integer startindex) {
    assert(scalvecdim == bundle_dim);
    if (xiz.dim() != vecdim) {
      compute_NTscaling();
    }

    for (Integer i = 0; i < vecdim; i++) {
      vec(startindex + i) = 0.;
    }
    if (use_scaling)
      vec(startindex + vecdim) = 1. / Chols;

    return 0;
  }


  ///add trace_dual*trace to diff_model for the right hand side (negative of current model violation)
  int BoxIPBundleBlock::add_trace_to_diff_model(Real trace_dual) {
    if ((use_scaling) && (scalub <= 0.))
      diff_model(vecdim) += trace_dual;
    return 0;
  }

  ///return the squared Euclidean norm of the dual model violation  
  Real BoxIPBundleBlock::dualviol_2normsqr() {
    tmpvec.init(diff_model);
    if (use_scaling) {
      mat_xmey(vecdim, tmpvec.get_store(), lz.get_store());
      mat_xpey(vecdim, tmpvec.get_store(), uz.get_store());
      tmpvec(vecdim) += -lt + ip(lz, lb) - ip(uz, ub) + (scalub > 0. ? ut : 0);
    } else {
      tmpvec -= lz;
      tmpvec += uz;
    }
    return mat_ip(scalvecdim, tmpvec.get_store());
  }


  /// If mu is not zero, always add the centering term for this mu as well;
  int BoxIPBundleBlock::set_modelrhs(Matrix& globalrhs,
    Real rhsmu,
    Real rhscorr,
    Integer startindex_model) {
    mat_xey(scalvecdim, globalrhs.get_store() + startindex_model, diff_model.get_store());
    return add_muxinv(globalrhs, startindex_model, rhsmu, rhscorr, true);
  }



  /// glob_lowrank(:,...)=(Btsyssqrtinv), trafotrace(...)=syssqrtinv*trace
  int BoxIPBundleBlock::Schur_transform_bundle(Matrix& glob_lowrank,
    MinorantBundle& globalbundle,
    Integer startindex_bundle,
    Matrix& trafotrace,
    Integer startindex_trace) {
    if (xiz.dim() != vecdim) {
      compute_NTscaling();
    }

    assert(trafotrace.rowdim() >= startindex_trace + scalvecdim);
    mat_xea(vecdim, trafotrace.get_store() + startindex_trace, 0.);
    if (use_scaling) {
      if (scalub < 0)
        trafotrace(startindex_trace + vecdim) = 1 / Chols;
      else
        trafotrace(startindex_trace + vecdim) = 0.;
    }

    assert(glob_lowrank.coldim() >= startindex_bundle + scalvecdim);

    for (Integer i = 0; i < vecdim; i++) {
      Real xizinv = 1. / xiz(i);
      Real sqrtxizinv = std::sqrt(xizinv);
      Real dummy;
      globalbundle[unsigned(startindex_bundle + map_to_old(i))].get_minorant(dummy, glob_lowrank, startindex_bundle + i, sqrtxizinv);
    }
    if (use_scaling) {
      Real dummy;
      globalbundle[unsigned(startindex_bundle + map_to_old(vecdim))].get_minorant(dummy, glob_lowrank, startindex_bundle + vecdim, 1. / Chols);
      for (Integer i = 0; i < vecdim; i++) {
        globalbundle[unsigned(startindex_bundle + map_to_old(i))].get_minorant(dummy, glob_lowrank, startindex_bundle + vecdim, Cholinvb(i), true);
      }
    }

    return 0;
  }

  /** @brief add diag(Bt*sqrt(invsys)*(I-lambda*trvec*trvec')*sqrt(invsys)*B) to diagonal

      @param[out] diagonal
         add the entries here

      @param[in] globalbundle
         the bundle vectors are [startindex_bundle ... startindex_bundle+vecdim()-1]

      @param[in] startindex_bundle
         see globalbundle

      @param[in] lambda
         may ==0., then use I in the middle
   otherwise use (I-lambda*trafortrace*trafotrace') in the middle

      @param[in] trafotrace
         holds precomputed invsys^(.5)*trace in coordinates
         [startindex_trace ... startindex_trace+vecdim()-1]

      @param[in] startindex_trace
         beginning of local trace within trafotracevec


   */
  int BoxIPBundleBlock::add_bundle_xizinv_diagonal(Matrix& diagonal,
    Matrix& ipBtrvec,
    MinorantBundle& globalbundle,
    Integer startindex_bundle,
    const Matrix& trafotrace,
    Integer startindex_trace) {
    if (xiz.dim() != vecdim) {
      compute_NTscaling();
    }

    if (!use_scaling) {

      assert(trafotrace.rowdim() >= startindex_trace + vecdim);
      Real scaleval;
      Minorant* mp;
      Integer nz;
      const Real* coeffs;
      const Integer* ind;
      for (Integer i = 0; i < bundle_dim; i++) {
        Real xizinv = 1. / sqrt_xiz(i);
        globalbundle[unsigned(startindex_bundle + map_to_old(i))].get_scaleval_and_minorant(scaleval, mp);
        if ((scaleval == 0.) || (mp == 0))
          continue;
        mp->get_coeffs(nz, coeffs, ind);
        if ((nz == 0) || (coeffs == 0))
          continue;
        xizinv *= scaleval;
        if (ind == 0) {
          for (Integer j = 0; j < nz; j++) {
            diagonal(j) += sqr(xizinv * coeffs[j]);
          }
        } else {
          for (Integer j = 0; j < nz; j++) {
            diagonal(ind[j]) += sqr(xizinv * coeffs[j]);
          }
        }
        //accumulate lambda part in tmpvec
        if (ipBtrvec.rowdim()) {
          xizinv *= trafotrace(startindex_trace + i);
          if (ind == 0) {
            for (Integer j = 0; j < nz; j++) {
              ipBtrvec(j) += xizinv * coeffs[j];
            }
          } else {
            for (Integer j = 0; j < nz; j++) {
              ipBtrvec(ind[j]) += xizinv * coeffs[j];
            }
          }
        }
      } //end for

    }

    else { // use_scaling == true

      tmpmat.init(diagonal.rowdim(), 1, 0.);  //accumulate scaling part
      assert(trafotrace.rowdim() >= startindex_trace + vecdim);
      Real scaleval;
      Minorant* mp;
      Integer nz;
      const Real* coeffs;
      const Integer* ind;
      for (Integer i = 0; i < bundle_dim - 1; i++) {
        Real xizinv = 1. / sqrt_xiz(i);
        globalbundle[unsigned(startindex_bundle + map_to_old(i))].get_scaleval_and_minorant(scaleval, mp);
        if ((scaleval == 0.) || (mp == 0))
          continue;
        mp->get_coeffs(nz, coeffs, ind);
        if ((nz == 0) || (coeffs == 0))
          continue;
        xizinv *= scaleval;
        Real cib = Cholinvb(i) * scaleval;
        if (ind == 0) {
          for (Integer j = 0; j < nz; j++) {
            diagonal(j) += sqr(xizinv * coeffs[j]);
            tmpmat(j) += cib * coeffs[j];
          }
        } else {
          for (Integer j = 0; j < nz; j++) {
            diagonal(ind[j]) += sqr(xizinv * coeffs[j]);
            tmpmat(ind[j]) += cib * coeffs[j];
          }
        }
        //accumulate lambda part in tmpvec
        if (ipBtrvec.rowdim() > 0) {
          xizinv *= trafotrace(startindex_trace + i);
          if (ind == 0) {
            for (Integer j = 0; j < nz; j++) {
              ipBtrvec(j) += xizinv * coeffs[j];
            }
          } else {
            for (Integer j = 0; j < nz; j++) {
              ipBtrvec(ind[j]) += xizinv * coeffs[j];
            }
          }
        }
      } //end for
      //add scaling column
      globalbundle[unsigned(startindex_bundle + bundle_dim - 1)].get_scaleval_and_minorant(scaleval, mp);
      if ((scaleval != 0.) && (mp != 0)) {
        mp->get_coeffs(nz, coeffs, ind);
        if ((nz != 0) && (coeffs != 0)) {
          Real xizinv = 1. / Chols;
          xizinv *= scaleval;
          if (ind == 0) {
            for (Integer j = 0; j < nz; j++) {
              tmpmat(j) += xizinv * coeffs[j];
            }
          } else {
            for (Integer j = 0; j < nz; j++) {
              tmpmat(ind[j]) += xizinv * coeffs[j];
            }
          }
        }
      }
      for (Integer i = 0; i < diagonal.rowdim(); i++) {
        diagonal(i) += sqr(tmpmat(i));
      }

      if (ipBtrvec.rowdim() > 0) {
        Real d = trafotrace(startindex_trace + bundle_dim - 1);
        for (Integer i = 0; i < diagonal.rowdim(); i++) {
          ipBtrvec(i) += d * tmpmat(i);
        }
      }

    }

    return 0;
  }


  int BoxIPBundleBlock::add_pcsubspace(Matrix& lowrank,
    Matrix& sigma_guess,
    const Matrix& Diag_inv,
    Real minval,
    Real diaginvval,
    Matrix& minus_trmult,
    Real schur_trace,
    MinorantBundle& globalbundle,
    Integer startindex_bundle) {
    if (vecdim == 0)
      return 0;

    if (xiz.dim() != vecdim) {
      compute_NTscaling();
    }

    assert(minval > 0.);

    Matrix localsqrBnorms;
    if (diaginvval > 0.) {
      if (sqrBnorms.rowdim() == 0) {
        sqrBnorms.init(scalvecdim, 1, 0.);
        for (Integer i = 0; i < scalvecdim; i++) {
          sqrBnorms(i) = globalbundle[unsigned(startindex_bundle + i)].norm_squared();
        }
      }
      localsqrBnorms.init(sqrBnorms, diaginvval);
    } else {
      localsqrBnorms.init(scalvecdim, 1, 0.);
      for (Integer i = 0; i < scalvecdim; i++) {
        localsqrBnorms(i) = globalbundle[unsigned(startindex_bundle + i)].norm_squared(&Diag_inv);
      }
    }


    //count the number of values exceeding the threshold
    Real threshold = minval * minval;
    Integer cnt = 0;
    for (Integer i = 0; i < vecdim; i++) {
      if (localsqrBnorms(i) / xiz(i) > threshold)
        cnt++;
    }

    if (cnt == 0)
      return 0;

    if (use_scaling) {
      //also include the last row
      tmpmat = Cholinvb;
      cnt++;
    }


    //reserve the memory
    Integer nextcol = lowrank.coldim();
    lowrank.enlarge_right(cnt, 0.);
    sigma_guess.enlarge_below(cnt, 0.);

    //append large enough vectors
    for (Integer i = 0; i < vecdim; i++) {
      Real d = localsqrBnorms(i) / sqrt_xiz(i);
      if (d > threshold) {
        if (cb_out(2)) {
          get_out() << " B(" << i << ";" << d << ")";
        }
        if (use_scaling) {
          //keep last row subspace close to orthogonal
          tmpmat(i) = 0;
        }
        sigma_guess(nextcol) = std::sqrt(d);
        Real dummy;
        globalbundle[unsigned(startindex_bundle + i)].get_minorant(dummy, lowrank, nextcol, 1. / sqrt_xiz(i), true);
        minus_trmult.concat_below(0.);
        nextcol++;
      }
    }

    //if vectors where appended, also append orthonormalized last row corrected for the trace term
    if (use_scaling) {

      //compute tmpmat = transpose(Cholinv)*(orthonormalized last row)
      const Real trcorr = sqr(sqr(1. / Chols)) / schur_trace;
      Real tmplast = 1. / sqr(Chols) - trcorr;
      Real nrm2 = std::sqrt(ip(tmpmat, tmpmat) + tmplast);
      tmplast = std::sqrt(tmplast) / nrm2;
      tmpmat /= nrm2;
      tmpmat /= sqrt_xiz;
      tmpmat.xpeya(Cholinvb, tmplast);
      tmpmat.concat_below(tmplast / Chols);

      //determine Bt*tmpmat
      for (Integer i = 0; i < tmpmat.rowdim(); i++) {
        Real dummy;
        globalbundle[unsigned(startindex_bundle + i)].get_minorant(dummy, lowrank, nextcol, tmpmat(i), true);
      }
      sigma_guess(nextcol) = minval * 1.1;
      if (cb_out(2)) {
        get_out() << " B(-;" << std::sqrt(colip(lowrank, nextcol)) << ")";
        //get_out()<<" B(-;"<<std::sqrt(colip(lowrank,nextcol,&Diag_inv))<<")";
      }
      minus_trmult.concat_below(tmplast / Chols);  //inner product for trace correction
      nextcol++;
    }

    assert(nextcol == lowrank.coldim());

    return 0;
  }



  /// add bundle*sqrt(inv(xiz))*subspace to glob_lowrank with bundle(:,si_bundle+1:si_bundle+dim_bundle()-1) and subspace(si_subsp:si_subsp+dim_bundle,:); sqrt(inv(xiz)) has to match that used in set_xizinvsqrt_trace()
  int BoxIPBundleBlock::add_bundle_xizinvsqrt_projection(Matrix& glob_lowrank,
    Matrix& subspace,
    Integer startindex_subspace,
    MinorantBundle& globalbundle,
    Integer startindex_bundle) {
    assert(scalvecdim == dim_bundle());
    if (xiz.dim() != vecdim) {
      compute_NTscaling();
    }

    if (!use_scaling) {

      if (startindex_subspace >= 0) {
        for (Integer i = 0; i < bundle_dim; i++) {

          MinorantPointer& mp = globalbundle[unsigned(startindex_bundle + i)];
          Real scalval = 1. / sqrt_xiz(i);
          for (Integer j = 0; j < subspace.coldim(); j++) {
            Real d = subspace(startindex_subspace + i, j) * scalval;
            if (d != 0.) {
              Real dummy;
              mp.get_minorant(dummy, glob_lowrank, j, d, true);
            }
          } //endfor j

        } //endfor i
      }

      else {
        Matrix tmpmat(glob_lowrank.coldim(), bundle_dim, 0.);
        times_B(glob_lowrank, tmpmat, 1., 0., 1, 1, 0, globalbundle, startindex_bundle);
        tmpvec = sqrt_xiz;
        tmpvec.inv();
        tmpmat.scale_cols(tmpvec);
        subspace.concat_right(tmpmat);
      }

    } //endif !use_scaling
    else {
      //use scaling

      if (startindex_subspace >= 0) {
        for (Integer i = 0; i < bundle_dim - 1; i++) {

          MinorantPointer& mp = globalbundle[unsigned(startindex_bundle + i)];
          Real scalval1 = 1. / sqrt_xiz(i);
          Real scalval2 = Cholinvb(i);
          for (Integer j = 0; j < subspace.coldim(); j++) {
            Real d = subspace(startindex_subspace + i, j) * scalval1 + subspace(startindex_subspace + vecdim, j) * scalval2;
            if (d != 0.) {
              Real dummy;
              mp.get_minorant(dummy, glob_lowrank, j, d, true);
            }
          } //endfor j

        }

        MinorantPointer& mp = globalbundle[unsigned(startindex_bundle + vecdim - 1)];
        for (Integer j = 0; j < subspace.coldim(); j++) {
          Real d = subspace(startindex_subspace + vecdim, j) / Chols;
          if (d != 0.) {
            Real dummy;
            mp.get_minorant(dummy, glob_lowrank, j, d, true);
          }
        } //endfor j
      }

      else {
        Matrix tmpmat(glob_lowrank.coldim(), bundle_dim, 0.);
        times_B(glob_lowrank, tmpmat, 1., 0., 1, 1, 0, globalbundle, startindex_bundle);
        Matrix tmpvec(glob_lowrank.coldim(), 1, tmpmat.get_store() + (bundle_dim - 1) * glob_lowrank.coldim(), 1, 1. / Chols);
        tmpmat.delete_cols(Range(bundle_dim - 1, bundle_dim - 1));
        genmult(tmpmat, Cholinvb, tmpvec, 1., 1.);
        Matrix scalvec(sqrt_xiz);
        scalvec.inv();
        tmpmat.scale_cols(scalvec);
        subspace.concat_right(tmpmat);
        subspace.concat_right(tmpvec);
      }

    } //endelse use scaling

    return 0;
  }


  /// out_vec+=BtinvsysB*in_vec
  int BoxIPBundleBlock::add_BtinvsysB_times(const Matrix& in_vec,
    Matrix& out_vec,
    Real zeta_inval,
    Real* zeta_outval,
    MinorantBundle& globalbundle,
    Integer startindex_bundle) {
    if (xiz.dim() != vecdim) {
      compute_NTscaling();
    }

    tmpvec.init(scalvecdim, 1, 0.);
    B_times(in_vec, tmpvec, 1., 0., 0, 0, 0, globalbundle, startindex_bundle);
    if (use_scaling)
      tmpvec(vecdim) -= zeta_inval;
    apply_xizinv(tmpvec, 0);
    B_times(tmpvec, out_vec, 1., 1., 1, 0, 0, globalbundle, startindex_bundle);
    if (use_scaling && zeta_outval)
      (*zeta_outval) -= tmpvec(vecdim);

    return 0;
  }

  /// compute dx (and then dz) given step_y and step_trdual on basis of the last rhs computed for the model block
  int BoxIPBundleBlock::set_dx_xizsolvestep(const Matrix& step_y,
    const Real step_trdual,
    MinorantBundle& globalbundle,
    Integer startindex_bundle) {
    //use of tmpvec is ok, because not employed in B_times() nor in set_dx_xizsolverhs();
    tmpvec.init(diff_model);
    if (use_scaling)
      tmpvec(vecdim) += step_trdual;
    B_times(step_y, tmpvec, -1., 1., 0, 0, 0, globalbundle, startindex_bundle);

    return set_dx_xizsolverhs(tmpvec, 0);
  }



  /// after the bundle subproblem is solved, this retrieves the local linear solution vector; if linx_activity is set, the values between zero and one indicate the guess on the coefficients activity level 
  int BoxIPBundleBlock::get_boxx(Matrix& boxx,
    Matrix* boxx_activity,
    Real trace_rhs,
    bool cautious) const {
    if (scalvecdim == bundle_dim) {
      boxx.init(vecdim, 1, x.get_store());
      if (use_scaling)
        boxx.concat_below(s);
      if (boxx_activity) {
        Indexmatrix inactive;
        find_inactive_indices(inactive, trace_rhs, cautious);
        boxx_activity->init(scalvecdim, 1, 1.);
        for (Integer i = 0; i < inactive.rowdim(); i++) {
          (*boxx_activity)(inactive(i)) = 0.;
        }
      }
    } else {
      boxx.init(bundle_dim, 1, 0.);
      for (Integer i = 0; i < scalvecdim; i++) {
        Integer ind = map_to_old(i);
        if ((use_scaling) && (ind == bundle_dim - 1))
          boxx(ind) = s;
        else
          boxx(ind) = x(i);
      }
      if (boxx_activity) {
        boxx_activity->init(bundle_dim, 1, 0.);
        for (Integer i = 0; i < scalvecdim; i++) {
          (*boxx_activity)(map_to_old(i)) = 1.;
        }
      }
    }

    return 0;
  }



}
