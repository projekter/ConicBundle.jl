/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/SOCIPProxBlock.cxx
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


#include "SOCIPProxBlock.hxx"

using namespace CH_Matrix_Classes;

namespace ConicBundle {

  int SOCIPProxBlock::clear_prox(const QPSolverProxObject* po) {
    Hp = po;

    rho = -1;
    drho = 0.;
    rho_rhs = 0.;
    rhoconst = 1.;

    sqrtHy_rhs.init(0, 0, 0.);

    ydim = 0;
    lrdim = 0;

    if (Hp) {
      Hp->get_precond(sqrtD, Vp);
      ydim = sqrtD.rowdim();
      sqrtD.sqrt();
      lrdim = (Vp ? Vp->coldim() : 0);
      SOCIPBlock::clear(2 + ydim + lrdim);
      if (!Hp->is_DLR()) {
        if (cb_out()) {
          get_out() << "**** WARNING SOCIPProxBlock::clear(.): proximal term is not of the form diagonal+low rank, not supported" << std::endl;
        }
        return 1;
      }
    } else {
      SOCIPBlock::clear();
    }

    return 0;
  }


  SOCIPProxBlock::SOCIPProxBlock(const QPSolverProxObject* po, CBout* cb, int cbinc) :
    CBout(cb, cbinc), SOCIPBlock(0, cb, cbinc) {
    clear_prox(po);
  }

  SOCIPProxBlock::~SOCIPProxBlock() {
  }

  int SOCIPProxBlock::reset_starting_point(const Matrix& y,
    const Matrix& sys_lhs,
    Real& mu) {
    assert(Hp);
    assert(y.rowdim() == ydim);
    assert(sys_lhs.rowdim() == ydim);
    assert(vecdim == ydim + lrdim + 2);

    //--- set x body
    x.newsize(vecdim, 1); chk_set_init(x, 1);
    Real sum_xsqrs = 0.;
    Real* xp = x.get_store() + 2;
    const Real* sp = sys_lhs.get_store();
    const Real* dp = sqrtD.get_store();
    for (Integer i = 0; i < ydim; i++) {
      Real d = (*sp++) / (*dp++);
      sum_xsqrs += sqr(d);
      *xp++ = d;
    }
    for (Integer i = 0; i < lrdim; i++) {
      *xp++ = 0.;
    }

    //--- compute sqrtHy_rhs withou rho
    sqrtHy_rhs.newsize(vecdim, 1); chk_set_init(sqrtHy_rhs, 1);
    Real* Hyp = sqrtHy_rhs.get_store();
    (*Hyp++) = 0.;
    (*Hyp++) = 1.;
    Real sum_zsqrs = 0.;
    const Real* yp = y.get_store();
    dp = sqrtD.get_store();
    for (Integer i = 0; i < ydim; i++) {
      Real d = (*yp++) * (*dp++);
      sum_zsqrs += sqr(d);
      *Hyp++ = d;
    }
    if (lrdim > 0) {
      genmult(*Vp, y, tmpvec, 1., 0., 1);
      dp = tmpvec.get_store();
      for (Integer i = 0; i < lrdim; i++) {
        Real d = (*dp++);
        *Hyp++ = d;
        sum_zsqrs += sqr(d);
      }
    }

    //--- set z body
    z.newsize(vecdim, 1); chk_set_init(z, 1);
    mat_xey(vecdim - 2, z.get_store() + 2, sqrtHy_rhs.get_store() + 2);

    //--- set rho, x0, x1, z0, z1
    rhoconst = max(1., max(std::sqrt(sum_xsqrs), std::sqrt(sum_zsqrs)) / 100.);
    x(0) = 1.001 * (sum_xsqrs / rhoconst + rhoconst + 2.) / 2.;
    rho = 1.001 * (sum_zsqrs / rhoconst + rhoconst + 2.) / 2.;
    Real d = (x(0) - rhoconst) * (rhoconst - rho) + mat_ip(ydim + lrdim, x.get_store() + 2, z.get_store() + 2);
    Real gmu = max(mu, x(0) * rho + d);
    d = std::fabs(gmu - d);
    x(0) = max(x(0), d);
    x(1) = x(0) - rhoconst;
    rho = max(rho, d / x(0));
    z(0) = rho;
    z(1) = rhoconst - rho;
    sqrtHy_rhs(0) = rho;
    sqrtHy_rhs(1) = rhoconst - rho;

    //--- set rho, x0, x1, z0, z1  //only care about centrality and not feasibility of x
    /*
    rhoconst=max(1.,std::sqrt(sum_zsqrs)/1000.);
    rho=1.001*(sum_zsqrs/rhoconst+rhoconst+2.)/2.;
    z(0)=rho;
    z(1)=rhoconst-rho;
    sqrtHy_rhs(0)=rho;
    sqrtHy_rhs(1)=rhoconst-rho;
    x.init(z,-1.);
    x(0)=rho;
    */
    return 0;
  }


  int SOCIPProxBlock::add_prox_contrib(Real& min_objective,
    Real& max_objective,
    Matrix& sys_lhs) {
    assert(sys_lhs.rowdim() == ydim);
    //std::cout<<"socsys_lhsin="<<sys_lhs;
    //std::cout<<" minin="<<min_objective<<std::endl;
    //std::cout<<" maxin="<<max_objective<<std::endl;

    min_objective += rhoconst * (rho - .5 * rhoconst);
    max_objective += rhoconst * (-x(1) - .5 * rhoconst);

    //-- add diagonal part to sys_lhs
    Real* sysp = sys_lhs.get_store();
    const Real* dp = sqrtD.get_store();
    const Real* xp = x.get_store() + 2;
    for (Integer i = 0; i < ydim; i++) {
      *sysp++ -= (*dp++) * (*xp++);
    }

    //-- add low rank part to sys_lhs
    if (Vp) {
      tmpvec.init(lrdim, 1, x.get_store() + ydim + 2);
      genmult(*Vp, tmpvec, sys_lhs, -1., 1.);
    }

    // std::cout.precision(12);
    // std::cout<<"norm2(sqrtHy_rhs-z)="<<norm2(sqrtHy_rhs-z)<<std::endl;
    // std::cout<<"n2(socsys_lhs)="<<norm2(sys_lhs)<<std::endl;
    // std::cout<<"abs(x(0)-x(1)-1.)="<<std::fabs(x(0)-x(1)-1.)<<std::endl;
    // std::cout<<" minout="<<min_objective<<std::endl;
    // std::cout<<" maxout="<<max_objective<<std::endl;

    return 0;
  }

  Real SOCIPProxBlock::primalviol_2normsqr() {
    tmpvec = sqrtHy_rhs;
    tmpvec -= z;
    return ip(tmpvec, tmpvec);
  }

  int SOCIPProxBlock::add_prox_sysrhs(Matrix& rhs,
    Real& Hfactor,
    Real rhsmu,
    Real rhscorr) {

    if (f.dim() != vecdim) {
      compute_NTscaling();
    }

    last_rhs_mu = rhsmu;

    if (rhscorr > 0.) {
      assert((dx.dim() == vecdim) && (dx.coldim() == 1));

      tmpvec.init(dx);
      apply_F(tmpvec);
      compl_rhs.init(dz, -rhscorr);
      apply_Arw(tmpvec, apply_Finv(compl_rhs));
    } else
      compl_rhs.init(vecdim, 1, 0.);

    if (rhsmu > 0.) {
      compl_rhs(0) += rhsmu;
    }

    if ((rhsmu > 0) || (rhscorr > 0)) {
      apply_Arwinv(scaled_point, compl_rhs);
      apply_Finv(compl_rhs);
    }

    /* //alternatively (does not seem to work so well?!?)
      apply_Arwinv(scaled_point,compl_rhs);
      apply_Finv(compl_rhs);
    }
    else
      compl_rhs.init(vecdim,1,0.);

    if (rhsmu>0.){
      compl_rhs(0)+=mu/gammazsqr*z(0);
      mat_xpeya(vecdim-1,compl_rhs.get_store()+1,z.get_store()+1,-mu/gammazsqr);
    }
    */

    //tmpvec.init(sHrhsmcompl);
    //Real testrho_rhs=(tmpvec(0)+tmpvec(1)-2*(f(0)+f(1))*(2*f(0)*tmpvec(0)-ip(tmpvec,f)))/sqr(omega)-1.+x(0)-x(1);

    //std::cout<<"sqrtHy_rhs= "<<sqrtHy_rhs<<" tmpvec="<<tmpvec;
    //apply_Finvsqr(tmpvec);

    compl_rhs -= x;
    tmpvec.init(compl_rhs, -1.);

    rho_rhs = -rhoconst + x(0) - x(1) - tmpvec(0) + tmpvec(1);
    //std::cout<<" errrho_rhs= "<<rho_rhs-testrho_rhs;
    //std::cout<<"rho_rhs="<<rho_rhs<<" x0="<<x(0)<<" x1="<<x(1)<<" Finvsqr(sqrtHy)="<<tmpvec<<std::endl;
    tmpvec(0) += rho_rhs * f(0) / (f(0) + f(1));
    mat_xpeya(vecdim - 1, tmpvec.get_store() + 1, f.get_store() + 1, -rho_rhs / (f(0) + f(1)));

    Real* rhsp = rhs.get_store();
    const Real* vp = tmpvec.get_store() + 2;
    const Real* dp = sqrtD.get_store();
    for (Integer i = 0; i < ydim; i++) {
      *rhsp++ -= (*dp++) * (*vp++);
    }

    if (Vp) {
      mat_xey(lrdim, tmpvec.get_store(), tmpvec.get_store() + 2 + ydim);
      tmpvec.reduce_length(lrdim);
      genmult(*Vp, tmpvec, rhs, -1., 1.);
    }

    Hfactor = 1. / sqr(omega);

    return 0;
  }

  int SOCIPProxBlock::compute_step(const Matrix& step) {
    dz = sqrtHy_rhs;
    dz -= z;
    assert(norm2(dz) < 1e-10 * z(0));
    const Real f0pf1 = f(0) + f(1);
    drho = rho_rhs * omega * omega / f0pf1 / 2.;

    const Real* sp = step.get_store();
    const Real* dp = sqrtD.get_store();
    const Real* fp = f.get_store() + 2;
    Real* dzp = dz.get_store() + 2;
    for (Integer i = 0; i < ydim; i++) {
      Real d = (*dp++) * (*sp++);
      (*dzp++) += d;
      drho += d * (*fp++);
    }

    if (Vp) {
      genmult(*Vp, step, tmpvec, 1., 0., 1);
      sp = tmpvec.get_store();
      for (Integer i = 0; i < lrdim; i++) {
        Real d = (*sp++);
        (*dzp++) += d;
        drho += d * (*fp++);
      }
    }
    drho /= f0pf1;

    dz(0) += drho;
    dz(1) -= drho;

    dx.init(dz, -1.);
    apply_Finvsqr(dx);
    dx += compl_rhs;

    // dx.init(dz,-1);
    // apply_Finvsqr(dx);
    // tmpvec.init(compl_rhs);
    // dx+=tmpvec;




    // //-- BEGIN Testing
    // Real n2compl=norm2(apply_F(tmpmat.init(dx-compl_rhs))+apply_Finv(tmpvec.init(dz)));
    // if (n2compl>1e-8*(z(0)+x(0))){
    //   std::cout.precision(4);
    //   std::cout<<" SOCcomplviol="<<n2compl<<"("<<norm2(apply_F(tmpmat.init(dx-compl_rhs))+apply_Finv(tmpvec.init(dz)))<<","<<norm2(dz-apply_Fsqr(apply_Finvsqr(tmpvec.init(dz))))<<")"<<" omega="<<omega<<" gammazsqr="<<gammazsqr<<" gammaxsqr="<<gammaxsqr<<" x(0)="<<x(0)<<" z(0)="<<z(0)<<" diff=";
    //   std::cout<<transpose(apply_F(tmpmat.init(dx-compl_rhs))+apply_Finv(tmpvec.init(dz)));
    //   std::cout<<" x="<<transpose(x)<<" z="<<transpose(z)<<" dx="<<transpose(dx)<<" dz="<<transpose(dz)<<" norm2(f)^2="<<mat_ip(vecdim,f.get_store())<<"("<<2*f(0)*f(0)-1<<") f="<<transpose(f);
    // }
    // tmpvec.init(vecdim,1,0.);
    // tmpvec(0)=-drho;
    // tmpvec(1)=drho;
    // sp=step.get_store();
    // dp=sqrtD.get_store();
    // Real *vp=tmpvec.get_store()+2;
    // for(Integer i=0;i<ydim;i++){
    //   (*vp++)=-(*dp++)*(*sp++);
    // }
    // if (Vp){
    //   genmult(*Vp,step,tmpmat,-1.,0.,1);
    //   mat_xey(lrdim,vp,tmpmat.get_store());
    // }
    // std::cout<<"dz_test= "<<norm2(tmpvec+z+dz-sqrtHy_rhs)<<std::endl;

    // tmpmat.init(dx,-1.);
    // apply_Fsqr(tmpmat);
    // Matrix tmpFcompl(compl_rhs+x);
    // apply_Fsqr(tmpFcompl);
    // tmpmat+=tmpvec;
    // tmpmat-=sqrtHy_rhs;
    // tmpmat+=tmpFcompl;
    // std::cout<<" dx-test= "<<norm2(tmpmat)<<std::endl;
    // std::cout<<" [-1 1 0_m](x+dx)+1= "<<1.+x(1)+dx(1)-x(0)-dx(0)<<std::endl;
    // //tmpmat=sqrtHy_rhs;
    // //tmpmat-=tmpFcompl;
    // //tmpmat-=tmpvec;
    // //apply_Finvsqr(tmpmat);
    // //std::cout<<" dx-test2= "<<norm2(tmpmat+dx)<<std::endl;
    // //tmpmat=sqrtHy_rhs;
    // //tmpmat-=tmpFcompl;
    // //std::cout<<"sqrtHy_rhs= "<<sqrtHy_rhs<<" tmpmat="<<tmpmat;
    // //apply_Finvsqr(tmpmat);    
    // //std::cout<<" rhorhs_test= "<<(rho_rhs+1-x(0)+x(1)+tmpmat(0)-tmpmat(1))<<std::endl;
    // //std::cout<<"rho_rhs="<<rho_rhs<<" x0="<<x(0)<<" x1="<<x(1)<<" Finvsqr(sqrtHy)="<<tmpmat<<std::endl;
    // apply_Finvsqr(tmpvec);
    // std::cout<<" drho-test= "<<-tmpvec(0)+tmpvec(1)-rho_rhs<<std::endl;
    // //--- END Testing

    return 0;
  }

  int SOCIPProxBlock::do_step(Real alpha,
    const Matrix& y) {
    rho += alpha * drho;
    drho = 0.;
    rho_rhs = 0.;

    int err = SOCIPBlock::do_step(alpha);

    //--- compute sqrtHy_rhs
    sqrtHy_rhs.newsize(vecdim, 1); chk_set_init(sqrtHy_rhs, 1);
    Real* Hyp = sqrtHy_rhs.get_store();
    (*Hyp++) = rho;
    (*Hyp++) = rhoconst - rho;
    const Real* yp = y.get_store();
    const Real* dp = sqrtD.get_store();
    for (Integer i = 0; i < ydim; i++) {
      Real d = (*yp++) * (*dp++);
      *Hyp++ = d;
    }
    if (lrdim > 0) {
      genmult(*Vp, y, tmpvec, 1., 0., 1);
      dp = tmpvec.get_store();
      for (Integer i = 0; i < lrdim; i++) {
        Real d = (*dp++);
        *Hyp++ = d;
      }
    }

    // std::cout<<" sqrtHinfeas= "<<norm2(sqrtHy_rhs-z)<<std::endl;
    // std::cout<<" pinfeas= "<<x(0)-1-x(1)<<std::endl;

    return err;
  }

  int SOCIPProxBlock::test_sysviol(Matrix& sys_lhs, const Matrix& step) {
    assert(sys_lhs.rowdim() == ydim);
    //std::cout<<"socsys_lhsin="<<sys_lhs;

    //-- add diagonal part to sys_lhs
    Real* sysp = sys_lhs.get_store();
    const Real* dp = sqrtD.get_store();
    const Real* xp = x.get_store() + 2;
    const Real* dxp = dx.get_store() + 2;
    for (Integer i = 0; i < ydim; i++) {
      *sysp++ -= (*dp++) * ((*xp++) + (*dxp++));
    }

    //-- add low rank part to sys_lhs
    if (Vp) {
      tmpvec.init(lrdim, 1, x.get_store() + ydim + 2);
      mat_xpey(lrdim, tmpvec.get_store(), dx.get_store() + ydim + 2);
      genmult(*Vp, tmpvec, sys_lhs, -1., 1.);
    }

    //output the system violation
    //std::cout<<"socsys_lhs="<<sys_lhs;
    std::cout << " [-1 1 0_m](x+dx)+1= " << rhoconst + x(1) + dx(1) - x(0) - dx(0) << std::endl;
    std::cout << " rho0-infeas= " << -rho - drho + z(0) + dz(0) << std::endl;
    std::cout << " rho1-infeas= " << rho + drho + z(1) + dz(1) - rhoconst << std::endl;

    tmpvec = sqrtHy_rhs;
    tmpvec -= z;
    tmpvec -= dz;
    tmpvec(0) += drho;
    tmpvec(1) -= drho;
    Real* vp = tmpvec.get_store() + 2;
    const Real* sp = step.get_store();
    dp = sqrtD.get_store();
    for (Integer i = 0; i < ydim; i++) {
      *vp++ += (*dp++) * (*sp++);
    }
    if (Vp) {
      genmult(*Vp, step, tmpmat, 1., 0., 1);
      mat_xpey(lrdim, vp, tmpmat.get_store());
    }
    std::cout << " sqrtHinfeas= " << norm2(tmpvec) << std::endl;
    std::cout << " socqpcompl= " << norm2(apply_F(tmpmat.init(dx - compl_rhs)) + apply_Finv(tmpvec.init(dz))) << std::endl;

    return 0;
  }


}
