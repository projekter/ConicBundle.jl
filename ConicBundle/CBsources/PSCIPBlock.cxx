/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/PSCIPBlock.cxx
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


#include "PSCIPBlock.hxx"

using namespace CH_Matrix_Classes;

namespace ConicBundle {

  void PSCIPBlock::point_changed() {
    dX.init(0, 0.);
    dZ.init(0, 0.);
    W.init(0, 0.);
    Winv.init(0, 0.);
    G.init(0, 0, 0.);
    Ginv.init(0, 0, 0.);
    D.init(0, 1, 0.);
    compl_rhs.init(0, 0.);
    lamX.init(0, 1, 0.);
    PX.init(rowdim, 0, 0.);
    Weig.init(0, 1, 0.);
    Wvec.init(0, 0, 0.);
  }

  int PSCIPBlock::compute_NTscaling() {
    bool use_Zchol = false;
    if (W.rowdim() != rowdim) {
      if (use_Zchol) {
        tmpsym = Z;
        if (tmpsym.Chol_factor(1e-20)) {
          if (cb_out()) get_out() << "*** WARNING: PSCIPBlock::compute_NTscaling(): factorizing Z failed" << std::endl;
          return 1;
        }
        tmpsym2.init(X);
        tmpsym.Chol_scaleLt(tmpsym2);
        tmpsym2.eig(tmpmat, tmpvec);
        D.newsize(rowdim, Integer(1)); chk_set_init(D, 1);
        for (Integer k = 0; k < rowdim; k++) {
          D(k) = std::sqrt(max(tmpvec(k), 1e-20));
          tmpvec(k) = std::sqrt(D(k));
        }
        Ginv.init(tmpmat);

        //G= Lam*P^T*Zchol^{-1}
        tmpsym.Chol_Ltsolve(tmpmat);
        tmpmat.scale_cols(tmpvec);
        G.init(tmpmat, 1., 1);

        //Ginv=Zchol*P*Lam^{-1}
        tmpsym.Chol_Lmult(Ginv);
        tmpvec.inv();
        Ginv.scale_cols(tmpvec);

        //assert(norm2(G*Ginv-tmpvec.init_diag(rowdim))<1e-6);
        //assert(norm2(G*Z*transpose(G)-tmpvec.init_diag(D))<1e-6);
        //assert(norm2(transpose(Ginv)*X*Ginv-tmpvec.init_diag(D))<1e-6);
      } else {
        //TEST begin
        //X.eig(tmpmat,tmpvec);
        //std::cout<<"eig(X)="<<tmpvec;
        //TEST end

        tmpsym = X;
        if (tmpsym.Chol_factor(1e-20)) {
          if (cb_out()) get_out() << "*** WARNING: PSCIPBlock::compute_NTscaling(): factorizing X failed" << std::endl;
          return 1;
        }
        tmpsym2.init(Z);
        tmpsym.Chol_scaleLt(tmpsym2);
        tmpsym2.eig(tmpmat, tmpvec);
        D.newsize(rowdim, Integer(1)); chk_set_init(D, 1);
        for (Integer k = 0; k < rowdim; k++) {
          D(k) = std::sqrt(max(tmpvec(k), 1e-20));
          tmpvec(k) = std::sqrt(D(k));
        }
        Ginv.init(tmpmat);

        //Ginv=Xcholt^{-1}*P*Lam
        tmpsym.Chol_Ltsolve(Ginv);
        Ginv.scale_cols(tmpvec);

        //G= Lam^{-1}*P^T*Xcholt
        tmpsym.Chol_Lmult(tmpmat);
        tmpvec.inv();
        tmpmat.scale_cols(tmpvec);
        G.init(tmpmat, 1., 1);

        //assert(norm2(G*Ginv-tmpvec.init_diag(rowdim))<1e-6);
        //assert(norm2(G*Z*transpose(G)-tmpvec.init_diag(D))<1e-6);
        //assert(norm2(transpose(Ginv)*X*Ginv-tmpvec.init_diag(D))<1e-6);

      }

      //W
      rankadd(G, W, 1., 0., 1);

      //Winv
      rankadd(Ginv, Winv);

      //TEST begin
      compute_Weig_Wvec();
      // std::cout<<" devdiag(W)="<<norm2(W-Diag(diag(W)))/trace(W);
      // std::cout<<" maxWeig="<<max(Weig)<<" minWeig="<<min(Weig)<<" Weig="<<transpose(Weig)<<std::endl;
      // std::cout<<" diag(W)="<<transpose(diag(W));

      //TEST end

      //assert(norm2(W*Winv-tmpvec.init_diag(rowdim))<1e-6);
      //assert(norm2(W*Z*W-X)<1e-6);
      //assert(norm2(Winv*X*Winv-Z)<1e-6);
    }
    return 0;
  }

  int PSCIPBlock::compute_Weig_Wvec() {
    int status = 0;
    if (Weig.rowdim() != rowdim) {
      status = W.eig(Wvec, Weig, false);
      if (status) {
        if (cb_out())
          get_out() << "\n**** WARNING PSCIPBlock::compute_Weig_Wvec(): W.eig failed and returned " << status << std::endl;
      }
      if (cb_out(5)) {
        get_out().precision(4);
        get_out() << " maxWeig=" << Weig(0) << " minWeig=" << Weig(rowdim - 1) << " Weig" << transpose(Weig);
        // //TEST begin
        // Matrix tmpvec,tmpmat;
        // X.eig(tmpmat,tmpvec,false);
        // std::cout<<" maxXeig="<<tmpvec(0)<<" minXeig="<<tmpvec(rowdim-1)<<" Xeig"<<transpose(tmpvec);
        // tmpvec=diag(transpose(tmpmat)*W*tmpmat);
        // std::cout<<" maxPeig="<<tmpvec(0)<<" minPeig="<<tmpvec(rowdim-1)<<" Peig"<<transpose(tmpvec);
        // //TEST end
      }

      // //TEST begin
      // Symmatrix skronW=skron(W,W);
      // Symmatrix W2(W*W);
      // Real trW2=trace(W2);
      // Matrix svecW2(svec(W2),std::sqrt(1./trW2));
      // rankadd(svecW2,skronW,-1.,1.);
      // Matrix Seig,Svec;
      // skronW.eig(Svec,Seig,false);
      // Symmatrix St(Weig.rowdim(),0.);
      // Matrix PSt;
      // Symmatrix prep;
      // std::cout<<" Wspec="<<norm2(W-Wvec*Diag(Weig)*transpose(Wvec))<<std::endl;
      // for(Integer i=0;i<Weig.rowdim();i++){
      // 	St(i,i)=sqr(Weig(i))*(1-sqr(Weig(i))/trW2);
      // 	rankadd(Wvec.col(i),prep);
      // 	PSt.concat_right(svec(prep));
      // 	for(Integer j=i+1;j<Weig.rowdim();j++)
      // 	  St(i,j)=-Weig(i)*Weig(j)*Weig(i)*Weig(j)/trW2;
      // }
      // //std::cout<<" transpose(PSt)*PSt="<<transpose(PSt)*PSt<<std::endl;
      // Matrix Steig,Stvec;
      // St.eig(Stvec,Steig,false);
      // Matrix Ptmp;
      // genmult(PSt,Stvec,Ptmp);
      // //std::cout<<" transpose(Ptmp)*Ptmp="<<transpose(Ptmp)*Ptmp<<std::endl;
      // std::cout<<" Steig="<<transpose(Steig);
      // std::cout<<" Weig="<<transpose(Weig);
      // std::cout<<" Seig="<<transpose(Seig);      
      // std::cout<<" skrWeig=\n";
      // for(Integer i=0;i<Weig.rowdim();i++){
      // 	std::cout<<" "<<Steig(i);
      // 	for(Integer j=i+1;j<Weig.rowdim();j++){
      // 	  std::cout<<" "<<Weig(i)*Weig(j);
      // 	  Steig.concat_below(Weig(i)*Weig(j));
      // 	  rank2add(Wvec.col(i),Wvec.col(j),prep,std::sqrt(2));
      // 	  Ptmp.concat_right(svec(prep));
      // 	}
      // }
      // std::cout<<std::endl;
      // std::cout<<" norm2(skronW-Ptmp*Diag(Steig)*transpose(Ptmp))="<<norm2(skronW-Ptmp*Diag(Steig)*transpose(Ptmp))<<std::endl;
      // std::cout<<" norm2(Diag(Matrix(Ptmp.rowdim(),1,1.))-Ptmp*transpose(Ptmp))="<<norm2(Diag(Matrix(Ptmp.rowdim(),1,1.))-Ptmp*transpose(Ptmp))<<std::endl;
      // //std::cout<<" Ptmp*transpose(Ptmp)="<<Ptmp*transpose(Ptmp)<<std::endl;
      // //TEST end


      Weig.sqrt();
      Wvec.scale_cols(Weig);
      Weig %= Weig;
    }
    return status;
  }


  void PSCIPBlock::clear(Integer dim) {
    rowdim = max(dim, 0);
    vecdim = (rowdim * (rowdim + 1)) / 2;
    X.init(rowdim, 0.);
    Z.init(rowdim, 0.);

    last_rhs_mu = 0.;
    mu = 0.;
    old_mu = 0.;
    last_alpha = 0.;
    oldX.init(rowdim, 0.);
    oldZ.init(rowdim, 0.);

    tmpsym.newsize(rowdim);
    tmpsym.init(0, 0.);
    tmpsym2.newsize(rowdim);
    tmpsym2.init(0, 0.);
    tmpvec.newsize(vecdim, 1);
    tmpvec.init(0, 0, 0.);
    tmpmat.init(0, 0, 0.);


    //reserve space but then set back to empty
    dX.newsize(rowdim);
    dZ.newsize(rowdim);
    W.newsize(rowdim);
    Winv.newsize(rowdim);
    G.newsize(rowdim, rowdim);
    Ginv.newsize(rowdim, rowdim);
    D.newsize(rowdim, Integer(1));
    compl_rhs.newsize(rowdim);
    point_changed();
  }

  PSCIPBlock::PSCIPBlock(Integer dim, CBout* cb, int cbinc) :CBout(cb, cbinc), InteriorPointBlock(cb, cbinc) {
    clear(dim);
  }

  PSCIPBlock::~PSCIPBlock() {
  }

  Integer PSCIPBlock::get_vecdim() const {
    return vecdim;
  }

  /// set x to value*"one" to x, or if add==true, add value*"one" to x
  int PSCIPBlock::center_x(Real val, bool add) {
    point_changed();
    if (!add) {
      X.init(rowdim, 0.);
    }
    X.shift_diag(val);
    return 0;
  }

  /// set z to value*"one" to z, or if add==true, add value*"one" to z
  int PSCIPBlock::center_z(Real val, bool add) {
    point_changed();
    if (!add) {
      Z.init(rowdim, 0.);
    }
    Z.shift_diag(val);
    return 0;
  }

  int PSCIPBlock::set_x(const Matrix& vec, Integer startindex, Real& add_center_value) {
    assert(vec.dim() >= startindex + vecdim);
    point_changed();

    const Real sqrt2 = std::sqrt(2);
    tmpvec.init(rowdim, 1, 0.);    // sum abs values of offdiagonals - value of diagonal 
    const Real* vp = vec.get_store() + startindex;
    Real* xp = X.get_store();
    Real* tip = tmpvec.get_store();
    for (Integer i = 0; i < rowdim; i++, tip++) {
      Real* tjp = tip;
      (*tjp++) -= ((*xp++) = (*vp++));
      for (Integer j = i + 1; j < rowdim; j++) {
        Real d = std::fabs(((*xp++) = (*vp++) / sqrt2));
        (*tip) += d;
        (*tjp++) += d;
      }
    }
    Real maxval = max(tmpvec); //this needs to be added to make it diagonally dominant
    if (maxval > 0.)
      add_center_value = maxval;
    else
      add_center_value = 0.;
    return 0;
  }

  int PSCIPBlock::set_z(const Matrix& vec, Integer startindex, Real& add_center_value) {
    assert(vec.dim() >= startindex + vecdim);
    point_changed();

    const Real sqrt2 = std::sqrt(2);
    tmpvec.init(rowdim, 1, 0.);    // sum abs values of offdiagonals - value of diagonal 
    const Real* vp = vec.get_store() + startindex;
    Real* zp = Z.get_store();
    Real* tip = tmpvec.get_store();
    for (Integer i = 0; i < rowdim; i++, tip++) {
      Real* tjp = tip;
      (*tjp++) -= ((*zp++) = (*vp++));
      for (Integer j = i + 1; j < rowdim; j++) {
        Real d = std::fabs(((*zp++) = (*vp++) / sqrt2));
        (*tip) += d;
        (*tjp++) += d;
      }
    }
    Real maxval = max(tmpvec); //this needs to be added to make it diagonally dominant
    if (maxval > 0.)
      add_center_value = maxval;
    else
      add_center_value = 0.;
    return 0;
  }

  /// on vec[startindex+0,+1 ...,+(vecdim-1)]   vec = b*vec + a * x for a real numbers a and b  
  int PSCIPBlock::vecgetsax(Matrix& vec, Integer startindex, Real a, bool add) {
    assert(vec.dim() >= startindex + vecdim);
    svec(X, vec, a, add, startindex);
    return 0;
  }

  /// on vec[startindex+0,+1 ...,+(vecdim-1)]   vec = b*vec + a * z for a real numbers a and b  
  int PSCIPBlock::vecgetsaz(Matrix& vec, Integer startindex, Real a, bool add) {
    assert(vec.dim() >= startindex + vecdim);
    svec(Z, vec, a, add, startindex);
    return 0;
  }


  int PSCIPBlock::get_mu_info(Integer& inmudim,
    Real& tr_xz,
    Real& tr_xdzpdxz,
    Real& tr_dxdz,
    Real& min_xz,
    Real& max_xz) const {
    assert(dX.rowdim() == rowdim);

    inmudim += rowdim;

    Real minv, maxv;
    tr_xz += ip_min_max(D, D, minv, maxv);
    tr_xdzpdxz += ip(X, dZ) + ip(dX, Z);
    tr_dxdz += ip(dX, dZ);

    if (minv < min_xz)
      min_xz = minv;
    if (maxv > max_xz)
      max_xz = maxv;

    if (cb_out(2))
      get_out() << " diagXZ[" << minv << "," << maxv << "]";

    return 0;
  }

  int PSCIPBlock::get_nbh_info(Integer inmudim,
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
    const Real mu_xz = tr_xz / inmudim;
    const Real mu_xdzpdxz = tr_xdzpdxz / inmudim;
    const Real mu_dxdz = tr_dxdz / inmudim;

    genmult(dZ, tmpvec.init(X), tmpmat);
    genmult(Z, tmpvec.init(dX), tmpmat, 1., 1.);
    genmult(G, tmpmat, tmpvec);
    genmult(tmpvec, Ginv, tmpmat);
    tmpsym2.init(tmpmat);
    for (Integer i = 0; i < rowdim; i++) {
      tmpsym2(i, i) -= mu_xdzpdxz;
    }
    //tmpsym2 now holds H_{G^{-T}}(XdZ+dXZ)-tr_xdzpdxz/mudim*I
    genmult(dX, Ginv, tmpmat);
    genmult(dZ, tmpmat, tmpvec);
    genmult(G, tmpvec, tmpmat);
    tmpsym.init(tmpmat);
    for (Integer i = 0; i < rowdim; i++)
      tmpsym(i, i) -= mu_dxdz;
    //tmpsym now holds H_{G^{-T}}(dXdZ)-tr_dxdz/mudim*I

    //Diag(D.*D) == H_{G^{-T}}(XZ))
    Real psc_nrmsqr_xz = 0;
    Real psc_nrmsqr_xdzpdxz = 0;
    Real psc_nrmsqr_dxdz = 0;
    Real psc_ip_xz_xdzpdxz = 0;
    Real psc_ip_xz_dxdz = 0;
    Real psc_ip_dxdz_xdzpdxz = 0;

    psc_nrmsqr_xdzpdxz += ip(tmpsym2, tmpsym2);
    for (Integer i = 0; i < rowdim; i++) {
      const Real d = sqr(D(i)) - mu_xz;
      psc_nrmsqr_xz += d * d;
      psc_ip_xz_xdzpdxz += d * tmpsym2(i, i);
    }
    psc_ip_dxdz_xdzpdxz += ip(tmpsym, tmpsym2);
    psc_nrmsqr_dxdz += ip(tmpsym, tmpsym);
    for (Integer i = 0; i < rowdim; i++) {
      const Real d = sqr(D(i)) - mu_xz;
      psc_ip_xz_dxdz += d * tmpsym(i, i);
    }

    nrmsqr_xz += psc_nrmsqr_xz;
    nrmsqr_xdzpdxz += psc_nrmsqr_xdzpdxz;
    nrmsqr_dxdz += psc_nrmsqr_dxdz;
    ip_xz_xdzpdxz += psc_ip_xz_xdzpdxz;
    ip_xz_dxdz += psc_ip_xz_dxdz;
    ip_dxdz_xdzpdxz += psc_ip_dxdz_xdzpdxz;

    int err = 0;
    if (alpha > 1000 * eps_Real) {
      err = control_nbh_step(alpha, max_nbh, nbh_ubnd,
        mu_xz, mu_xdzpdxz, mu_dxdz,
        psc_nrmsqr_xz, psc_nrmsqr_xdzpdxz, psc_nrmsqr_dxdz,
        psc_ip_xz_xdzpdxz, psc_ip_xz_dxdz, psc_ip_dxdz_xdzpdxz);
      if (err) {
        if (cb_out()) {
          get_out() << "*** ERROR PSCIPBlock::get_nbh_info(): control_nbh_step(.) returned " << err << std::endl;
        }
      }
    }

    return err;
  }


  /// if necessary, reduce alpha to the biggest value so that feasibility is maintained with this step size
  int PSCIPBlock::linesearch(Real& alpha) const {
    assert(dX.rowdim() == rowdim);

    Real a = min(alpha, 1.001);

    tmpsym.init(X);
    tmpsym.xpeya(dX, a);
    while ((tmpsym.Chol_factor(eps_Real)) && (a > eps_Real)) {
      a *= .8;
      tmpsym.init(X);
      tmpsym.xpeya(dX, a);
    }
    tmpsym.init(Z);
    tmpsym.xpeya(dZ, a);
    while ((tmpsym.Chol_factor(eps_Real)) && (a > eps_Real)) {
      a *= .8;
      tmpsym.init(Z);
      tmpsym.xpeya(dZ, a);
    }

    if (a < alpha)
      alpha = a;
    return 0;
  }

  /// compute the complementarity_rhs=rhsmu*xi-rhscorr*xi*dx*dz (wihtout "-z") for mu=rhsmu and for corrector for factor rhscorr>0., store this and add it to rhs 
  int PSCIPBlock::add_muxinv(Matrix& rhs,
    Integer startindex,
    Real rhsmu,
    Real rhscorr,
    bool minus) {
    assert(rhs.dim() >= startindex + vecdim);
    assert(rhscorr >= 0.);
    assert(rhsmu >= 0.);

    compute_NTscaling();

    //no -Z, because it cancels out

    if (rhscorr > 0.) {
      //form  bar_dZ=G*dZ*G^T
      Symmatrix bardX;
      symscale(dX, Ginv, bardX);
      //collect the corrector term
      tmpsym = Diag(D);
      /*
      if (last_rhs_mu>0.){
  for (Integer i=0;i<rowdim;i++){
    tmpsym(i,i)-=last_rhs_mu/D(i);
  }
      }
      */
      tmpsym += bardX;
      tmpsym.init(bardX * tmpsym);
      for (Integer i = 0; i < rowdim; i++) {
        for (Integer j = i; j < rowdim; j++) {
          tmpsym(i, j) /= (D(i) + D(j)) / 2.;
        }
      }
      if (rhsmu > 0.) {
        for (Integer i = 0; i < rowdim; i++) {
          tmpsym(i, i) += rhsmu / D(i);
        }
      }
    } else if (rhsmu > 0.) {
      //muval*Xinv
      tmpsym.init(rowdim, 0.);
      for (Integer i = 0; i < rowdim; i++) {
        tmpsym(i, i) += rhsmu / D(i);
      }
    }

    //complete and add the corrector term
    if ((rhsmu > 0.) || (rhscorr > 0.)) {
      symscale(tmpsym, Ginv, compl_rhs, 1., 0., 1);
      if (minus)
        svec(compl_rhs, rhs, -1., true, startindex);
      else
        svec(compl_rhs, rhs, 1., true, startindex);
    } else {
      compl_rhs.init(rowdim, 0.);
    }

    last_rhs_mu = rhsmu;

    return 0;
  }

  /// extract dx from rhs at startindex and compute at the same time dz (=-sys dx -z +complentarity_rhs); 
  int PSCIPBlock::set_dx(const CH_Matrix_Classes::Matrix& rhs, CH_Matrix_Classes::Integer startindex) {
    assert(W.rowdim() == rowdim);
    assert(rhs.dim() >= startindex + vecdim);
    assert(compl_rhs.rowdim() == rowdim);

    //get_sveci(dX,rhs,startindex);
    sveci(rhs, dX, 1., false, startindex, -1, rowdim);

    //TEST
    //Matrix tmpvec(rhs);
    //mat_xea(vecdim,tmpvec.get_store()+startindex,-10000.);
    //svec(dX,tmpvec,1.,false,startindex);
    //std::cout<<" svecnorm="<<norm2(tmpvec-rhs);

    tmpsym.init(dX, -1);
    tmpmat.init(Winv);
    symscale(tmpsym, tmpmat, dZ);
    dZ -= Z;
    dZ += compl_rhs;

    // std::cout<<" setdX="<<dX;
    // std::cout<<" dZ="<<dZ;


    //assert(norm2(Symmatrix(G*transpose(X*Z+dX*Z+X*dZ)*Ginv)-Symmatrix(G*compl_rhs*X*Ginv))<1e-8*(trace(X)+trace(Z)));
    return 0;
  }




  /// compute dx=sysinv*rhs and at the same time dz (=-rhs -z +complentarity_rhs); 
  int PSCIPBlock::set_dx_xizsolverhs(const Matrix& rhs, Integer startindex) {
    assert(W.rowdim() == rowdim);
    assert(rhs.dim() >= startindex + vecdim);
    assert(compl_rhs.rowdim() == rowdim);

    sveci(rhs, dZ, -1., false, startindex, -1, rowdim);
    tmpmat.init(W);
    symscale(dZ, tmpmat, dX, -1.);
    dZ += compl_rhs;
    dZ -= Z;

    //assert(norm2(Symmatrix(G*transpose(X*Z+dX*Z+X*dZ)*Ginv)-Symmatrix(G*compl_rhs*X*Ginv))<1e-8*(trace(X)+trace(Z)));
    return 0;
  }

  /// compute sysinv*rhs into rhs 
  int PSCIPBlock::apply_xizinv(Matrix& rhs, Integer startindex, bool minus) {

    compute_NTscaling();

    sveci(rhs, tmpsym, 1., false, startindex, 0, rowdim);
    tmpmat.init(W);
    symscale(tmpsym, tmpmat, tmpsym2);

    // //TEST begin
    // Symmatrix testsym1,testsym2;
    // symscale(tmpsym,G,testsym1,1.,0.,1);
    // symscale(testsym1,G,testsym2);
    // std::cout<<" GWdiff="<<norm2(tmpmat-transpose(G)*G)<<std::endl;
    // std::cout<<" Wdiff="<<norm2(tmpsym2-W*tmpsym*W)<<std::endl;
    // std::cout<<" tmpmatdiff="<<norm2(tmpsym2-transpose(tmpmat)*tmpsym*tmpmat)<<std::endl;
    // std::cout<<" testsym1diff="<<norm2(testsym1-G*tmpsym*transpose(G))<<std::endl;
    // std::cout<<" testsym2diff="<<norm2(testsym2-transpose(G)*testsym1*G)<<std::endl;
    // std::cout<<" xizinvdiff="<<norm2(tmpsym2-testsym2)<<std::endl;
    // //TEST end

    svec(tmpsym2, rhs, (minus ? -1. : 1.), false, startindex);

    return 0;
  }

  /// compute sysinv*rhs into rhs 
  int PSCIPBlock::apply_xiz(Matrix& rhs, Integer startindex, bool minus) {

    compute_NTscaling();

    sveci(rhs, tmpsym, 1., false, startindex, 0, rowdim);
    tmpmat.init(Winv);
    symscale(tmpsym, tmpmat, tmpsym2);
    svec(tmpsym2, rhs, (minus ? -1. : 1.), false, startindex);

    return 0;
  }

  /// move to (x+alpha*dx, z+alpha*dz)
  int PSCIPBlock::do_step(Real alpha) {
    assert(dX.rowdim() == rowdim);

    if ((old_mu == 0.) || (last_rhs_mu < mu)) {
      oldX = X;
      oldZ = Z;
      old_mu = mu;
      last_alpha = alpha;
    }
    mu = last_rhs_mu;

    X.xpeya(dX, alpha);
    Z.xpeya(dZ, alpha);

    if (cb_out(4)) {
      X.eig(PX, lamX, false);
      genmult(Z, PX, tmpmat);
      tmpvec.newsize(X.rowdim(), 1); chk_set_init(tmpvec, 1);
      for (Integer i = 0; i < X.rowdim(); i++) {
        tmpvec(i) = mat_ip(X.rowdim(), PX.get_store() + i * X.rowdim(), tmpmat.get_store() + i * X.rowdim());
      }
      Indexmatrix sind;
      sortindex(tmpvec, sind);
      get_out().precision(10);
      Real avgsl = ip(lamX, tmpvec);
      get_out() << " avgsl=" << avgsl << std::endl;
      for (Integer i = 0; i < X.rowdim(); i++) {
        Integer ind = sind(i);
        get_out() << "[" << std::setw(2) << ind << "," << std::setw(12) << tmpvec(ind) << "," << std::setw(12) << lamX(ind) << "]" << std::endl;
      }
      Real val0 = tmpvec(sind(0));
      tmpvec -= val0;
      lamX /= sum(lamX);
      Real devval = ip(lamX, tmpvec);
      get_out() << " devval=" << devval << "(" << sum(tmpvec) / X.rowdim() << ")" << std::endl;
    }

    point_changed();

    return 0;
  }

  /// add the Schur complement to a big system matrix
  int PSCIPBlock::add_AxizinvAt(const Matrix& A,
    Symmatrix& globalsys,
    bool minus,
    bool Atrans) {
    assert((Atrans == false ? A.coldim() : A.rowdim()) == vecdim);
    assert((Atrans == false ? A.rowdim() : A.coldim()) == globalsys.rowdim());

    compute_NTscaling();

    Integer ncols = Atrans ? A.coldim() : A.rowdim();
    tmpmat.newsize(vecdim, ncols);
    /*
    tmpmat.init(vecdim,0,0.);

    for(Integer i=0;i<ncols;i++){
      if (Atrans)
  tmpvec.init(vecdim,1,A.get_store()+i*A.rowdim());
      else
  tmpvec.init(vecdim,1,A.get_store()+i,A.rowdim());

      sveci(tmpvec,tmpsym);
      symscale(tmpsym,G,tmpsym2,1.,0.,1);
      svec(tmpsym2,tmpvec);
      tmpmat.concat_right(tmpvec);
    }
    */

    chk_set_init(tmpmat, 1);

    for (Integer i = 0; i < ncols; i++) {
      if (Atrans)
        tmpsym.init_svec(rowdim, A.get_store() + i * A.rowdim());
      else
        tmpsym.init_svec(rowdim, A.get_store() + i, A.rowdim());

      symscale(tmpsym, G, tmpsym2, 1., 0., 1);
      tmpsym2.store_svec(tmpmat.get_store() + i * vecdim);
    }

    rankadd(tmpmat, globalsys, minus ? -1. : 1., 1., 1);

    return 0;
  }

  /// add the system matrix*factor into a big system matrix starting at startindex
  int PSCIPBlock::add_xiz(Symmatrix& globalsys, Integer startindex, bool minus) {

    compute_NTscaling();

    if (minus)
      skron(Winv, Winv, globalsys, -1., true, startindex);
    else
      skron(Winv, Winv, globalsys, 1., true, startindex);


    return 0;
  }


}
