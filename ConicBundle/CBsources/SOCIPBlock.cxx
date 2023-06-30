/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/SOCIPBlock.cxx
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


#include "SOCIPBlock.hxx"

using namespace CH_Matrix_Classes;

namespace ConicBundle {


  void SOCIPBlock::point_changed()
  {
    dx.init(0,1,0.);
    dz.init(0,1,0.);
    gammaxsqr=0.;   
    gammazsqr=0.;  
    omega=0.; 
    f.init(0,1,0.);
    scaled_point.init(0,1,0.);
    compl_rhs.init(0,1,0.);
  }

// v gets Arw(u)*v

  Matrix& SOCIPBlock::apply_Arw(const Matrix& u,Matrix& v) const
  {
    assert(u.dim()==vecdim);
    assert(v.rowdim()==vecdim);
    assert(u.get_store()!=v.get_store());
    
    Real *vp=v.get_store();
    const Real* const up=u.get_store();
    Real u0=(*up);
    
    for(Integer i=0;i<v.coldim();i++,vp+=vecdim){
      Real v0=(*vp);
      Real ip_baruv=mat_ip(vecdim-1,up+1,vp+1);
      (*vp)=u0*v0+ip_baruv;
      mat_xbpeya(vecdim-1,vp+1,up+1,v0,u0);
    }
    
    return v;
  }
    
// v gets Arwinv(u)*v    
 
  Matrix& SOCIPBlock::apply_Arwinv(const Matrix& u, Matrix& v) const
  {
    assert(u.dim()==vecdim);
    assert(v.rowdim()==vecdim);
    assert(u.get_store()!=v.get_store());
    
    Real *vp=v.get_store();
    const Real* const up=u.get_store();
    Real u0=(*up);
    Real gamma2=u0*u0-mat_ip(vecdim-1,up+1);
    
    for(Integer i=0;i<v.coldim();i++,vp+=vecdim){
      Real v0=(*vp);
      Real ip_baruv=mat_ip(vecdim-1,up+1,vp+1);
      (*vp)=(u0*v0-ip_baruv)/gamma2;
      mat_xbpeya(vecdim-1,vp+1,up+1,(-v0+ip_baruv/u0)/gamma2,1./u0);
    }
    
    return v;
  }

// v gets F(f)*v    
 
  Matrix& SOCIPBlock::apply_F(Matrix& v) const
  {
    assert((f.rowdim()==vecdim)&&(f.coldim()==1));
    assert(v.rowdim()==vecdim);
  
    Real* vp=v.get_store();
    const Real* const fp=f.get_store();
    Real f0=(*fp);

    for (Integer i=0;i<v.coldim();i++,vp+=vecdim){
      Real v0=(*vp);
      Real ip_barfv=mat_ip(vecdim-1,fp+1,vp+1);
      (*vp)=(v0*f0+ip_barfv)*omega;
      mat_xbpeya(vecdim-1,vp+1,fp+1,(v0+ip_barfv/(1.+f0))*omega,omega);
    }
    return v;
  }

// v gets Finv(f)*v   
 
  Matrix& SOCIPBlock::apply_Finv(Matrix& v) const
  {
    assert((f.rowdim()==vecdim)&&(f.coldim()==1));
    assert(v.rowdim()==vecdim);

    Real* vp=v.get_store();
    const Real* const fp=f.get_store();
    Real f0=(*fp);

    for (Integer i=0;i<v.coldim();i++,vp+=vecdim){
      Real v0=(*vp);
      Real ip_barfv=mat_ip(vecdim-1,fp+1,vp+1);
      (*vp)=(v0*f0-ip_barfv)/omega;
      mat_xbpeya(vecdim-1,vp+1,fp+1,(-v0+ip_barfv/(1.+f0))/omega,1./omega);
    }
    return v;
  }

  
// v gets F(f)^2*v,    this is xiz
 
  Matrix& SOCIPBlock::apply_Fsqr(Matrix& v,bool minus) const
  {
    assert((f.rowdim()==vecdim)&&(f.coldim()==1));
    assert(v.rowdim()==vecdim);

    Real* vp=v.get_store();
    const Real* const fp=f.get_store();
    Real f0=(*fp);
    Real pos00=2*f0*f0-1; 
    Real o2=minus?-omega*omega:omega*omega;
    
    for (Integer i=0;i<v.coldim();i++,vp+=vecdim){
      Real v0=(*vp);
      Real ip_barfv=mat_ip(f.dim()-1,fp+1,vp+1);
      (*vp)=(v0*pos00+2*f0*ip_barfv)*o2;
      mat_xbpeya(vecdim-1,vp+1,fp+1,2*(f0*v0+ip_barfv)*o2,o2);
    }
    
    return v;
  }

// v gets F(f)^2*v,    this is xiz
 
  int SOCIPBlock::apply_Fsqr(Real* vp,bool minus) const
  {
    assert((f.rowdim()==vecdim)&&(f.coldim()==1));

    const Real* const fp=f.get_store();
    Real f0=(*fp);
    Real pos00=2*f0*f0-1.; 
    Real o2=minus?-omega*omega:omega*omega;
    
    Real v0=(*vp);
    Real ip_barfv=mat_ip(f.dim()-1,fp+1,vp+1);
    (*vp)=(v0*pos00+2*f0*ip_barfv)*o2;
    mat_xbpeya(vecdim-1,vp+1,fp+1,2*(f0*v0+ip_barfv)*o2,o2);
    
    return 0;
  }

// v gets Finv(f)^2*v,    this is xizinv 
 
  Matrix& SOCIPBlock::apply_Finvsqr(Matrix& v,bool minus) const
  {
    assert((f.rowdim()==vecdim)&&(f.coldim()==1));
    assert(v.rowdim()==vecdim);

    Real* vp=v.get_store();
    const Real* const fp=f.get_store();
    Real f0=(*fp);
    Real pos00=2*f0*f0-1.; 
    Real o2inv=(minus?-1.:1.)/(omega*omega);
    
    for (Integer i=0;i<v.coldim();i++,vp+=vecdim){
      Real v0=(*vp);
      Real ip_barfv=mat_ip(f.dim()-1,fp+1,vp+1);
      (*vp)=(v0*pos00-2*f0*ip_barfv)*o2inv;
      mat_xbpeya(vecdim-1,vp+1,fp+1,2.*(-f0*v0+ip_barfv)*o2inv,o2inv);
    }
    
    return v;
  }

// v gets Finv(f)^2*v,    this is xizinv 
 
  int SOCIPBlock::apply_Finvsqr(Real* vp,bool minus) const
  {
    assert((f.rowdim()==vecdim)&&(f.coldim()==1));

    const Real* const fp=f.get_store();
    Real f0=(*fp);
    Real pos00=2*f0*f0-1; 
    Real o2inv=(minus?-1.:1.)/(omega*omega);
    
    Real v0=(*vp);
    Real ip_barfv=mat_ip(f.dim()-1,fp+1,vp+1);
    (*vp)=(v0*pos00-2*f0*ip_barfv)*o2inv;
    mat_xbpeya(vecdim-1,vp+1,fp+1,2*(-f0*v0+ip_barfv)*o2inv,o2inv);
    
    return 0;
  }

  
  int SOCIPBlock::compute_NTscaling()
  {
    const Real* xp=x.get_store();
    const Real* zp=z.get_store();
    gammaxsqr=sqr(*xp);
    gammaxsqr-=mat_ip(vecdim-1,xp+1);
    if (gammaxsqr<=0.){
      if (cb_out()){
	get_out().precision(12);
	get_out()<<"**** ERROR SOCIPBlock::compute_NTscaling(): gammaxsqr="<<gammaxsqr<<"< =0 but has to be >0. ("<<sqr(x(0))<<","<<mat_ip(vecdim-1,x.get_store()+1)<<")"<<std::endl;
      }
      return 1;
    }
    gammazsqr=sqr(*zp);
    gammazsqr-=mat_ip(vecdim-1,zp+1);
    if (gammazsqr<=0.){
      if(cb_out()){
	get_out().precision(12);
	get_out()<<"**** ERROR SOCIPBlock::compute_NTscaling(): gzsqr="<<gammazsqr<<"<=0. but has to be >0. ("<<sqr(z(0))<<","<<mat_ip(vecdim-1,z.get_store()+1)<<")"<<std::endl;
      }
      return 1;
    }
    omega=std::sqrt(std::sqrt(gammazsqr/gammaxsqr));
    f.init(x,-omega);
    f(0)*=-1.;
    f.xpeya(z,1./omega);
    
    f*=1./std::sqrt(2*(std::sqrt(gammaxsqr*gammazsqr)+ip(x,z)));

    scaled_point.init(x);
    apply_F(scaled_point);

    // //TEST begin
    // //std::cout<<"f_0^2-bf^Tbf="<<2*sqr(f(0))-ip(f,f)<<std::endl;
    // apply_Finv(tmpvec.init(scaled_point));
    // //std::cout<<" x="<<x;
    // //std::cout<<" Finv(F(x))="<<tmpvec;
    // std::cout<<" norm2(x-Finv(F(x)))="<<norm2(x-tmpvec)<<std::endl;
    // apply_Finv(tmpvec);
    // apply_Fsqr(tmpvec);    
    // std::cout<<" norm2(sp-Fsqr(Finv(Finv(sp))))="<<norm2(scaled_point-tmpvec)<<std::endl;
    
    // apply_Finv(tmpvec.init(z));
    // //std::cout<<" z="<<z;
    // //std::cout<<" F(x)=scaled_point="<<scaled_point;
    // //std::cout<<" Finv(z)="<<tmpvec;
    // std::cout<<" norm2(F(x)-Finv(z))="<<norm2(scaled_point-tmpvec)<<std::endl;
    // apply_F(tmpvec);
    // //std::cout<<" F(Finv(z))="<<tmpvec;
    // std::cout<<" norm2(z-F(Finv(z)))="<<norm2(z-tmpvec)<<std::endl;
    // apply_F(tmpvec);
    // apply_Finvsqr(tmpvec);
    // std::cout<<" norm2(sp-Finvsqr(F(F(sp))))="<<norm2(scaled_point-tmpvec)<<std::endl;

    // tmpvec.rand(vecdim,1);
    // Matrix tmpvec2(tmpvec);
    // apply_F(tmpvec);
    // apply_F(tmpvec);
    // apply_Fsqr(tmpvec2);
    // //std::cout<<"F(F(tv))="<<tmpvec;
    // //std::cout<<"Fsqr(tv2)="<<tmpvec2;
    // std::cout<<"norm2(tmpvec-tmpvec2)="<<norm2(tmpvec-tmpvec2)<<std::endl;
    // tmpvec.rand(vecdim,1);
    // tmpvec2.init(tmpvec);
    // apply_Finv(tmpvec);
    // apply_Finv(tmpvec);
    // apply_Finvsqr(tmpvec2);
    // //std::cout<<"Finv(Finv(tv))="<<tmpvec;
    // //std::cout<<"Finvsqr(tv2)="<<tmpvec2;
    // std::cout<<"norm2(tmpvec-tmpvec2)="<<norm2(tmpvec-tmpvec2)<<std::endl;

    
    // tmpvec2.rand(vecdim,1)+=.1;
    // tmpvec.rand(vecdim,1);
    // Matrix tmpvec3(tmpvec);
    // apply_Arw(tmpvec2,tmpvec);
    // apply_Arwinv(tmpvec2,tmpvec);
    // std::cout<<"norm2(tv-Arwinv(Arw(tv)))="<<norm2(tmpvec-tmpvec3)<<std::endl;
    // //TEST end
    
    // assert(norm2(scaled_point-apply_Finv(tmpvec.init(z)))<1e-8*norm2(f));
    
    return 0;
  }
  
  void SOCIPBlock::clear(Integer dim)
  {
    vecdim=max(dim,0);
    x.init(vecdim,1,0.);
    z.init(vecdim,1,0.);
    
    last_rhs_mu=0.;
    mu=0.;
    old_mu=0.;
    last_alpha=0.;
    oldx.init(vecdim,1,0.);
    oldz.init(vecdim,1,0.);

    tmpvec.newsize(vecdim,1);
    tmpvec.init(0,0,0.);
    tmpmat.init(0,0,0.);
	      
//reserve space but then set back to empty
    dx.newsize(vecdim,Integer(1));
    dz.newsize(vecdim,Integer(1));
    gammaxsqr=0;
    gammazsqr=0;
    omega=0.;
    f.newsize(vecdim,Integer(1));
    scaled_point.newsize(vecdim,Integer(1));
    compl_rhs.newsize(vecdim,Integer(1));
    point_changed();
  }
  
  SOCIPBlock::SOCIPBlock(Integer dim,CBout* cb,int cbinc):InteriorPointBlock(cb,cbinc)
  {clear(dim);}

  SOCIPBlock::~SOCIPBlock()
  {}
  
  Integer SOCIPBlock::get_vecdim() const
  {return vecdim;}
  
  /// set x to value*"one" to x, or if add==true, add value*"one" to x
  int SOCIPBlock::center_x(Real val,bool add)
  {
    point_changed();
    if (!add) {
      x.init(vecdim,1,0.);
    }
    x(0)+=val;
    return 0;
  }

  /// set z to value*"one" to z, or if add==true, add value*"one" to z
  int SOCIPBlock::center_z(Real val,bool add)
  {
    point_changed();
    if (!add) {
      z.init(vecdim,1,0.);
    }
    z(0)+=val;
    return 0;
  }

  int SOCIPBlock::set_x(const Matrix& vec,Integer startindex,Real& add_center_value)
  {
    assert(vec.dim()>=startindex+vecdim);
    point_changed();
    const Real *vp=vec.get_store()+startindex;
    Real *xp=x.get_store();
    const Real * const xend = xp+vecdim;
    Real x0sqr=(*xp++)=(*vp++);
    x0sqr*=x0sqr;
    Real normsqr=0.;
    while(xp!=xend){
      Real d=*vp++;
      (*xp++)=d;
      normsqr+=d*d;
    }
    if (normsqr>x0sqr)
      add_center_value=std::sqrt(normsqr)-x(0);
    else 
      add_center_value=0;
    return 0;
  }

  int SOCIPBlock::set_z(const Matrix& vec,Integer startindex,Real& add_center_value)
  {
    assert(vec.dim()>=startindex+vecdim);
    point_changed();
    const Real *vp=vec.get_store()+startindex;
    Real *zp=z.get_store();
    const Real * const zend = zp+vecdim;
    Real z0sqr=(*zp++)=(*vp++);
    z0sqr*=z0sqr;
    Real normsqr=0.;
    while(zp!=zend){
      Real d=*vp++;
      (*zp++)=d;
      normsqr+=d*d;
    }
    if (normsqr>z0sqr)
      add_center_value=std::sqrt(normsqr)-z(0);
    else 
      add_center_value=0;
    return 0;
  }

  /// on vec[startindex+0,+1 ...,+(vecdim-1)]   vec = b*vec + a * x for a real numbers a and b  
  int SOCIPBlock::vecgetsax(Matrix& vec,Integer startindex,Real a,bool add)
  {
    assert(vec.dim()>=startindex+vecdim);
    if (add==false){
      mat_xeya(vecdim,vec.get_store()+startindex,x.get_store(),a);
    }
    else {
      mat_xpeya(vecdim,vec.get_store()+startindex,x.get_store(),a);
    }
    return 0;
  }

  /// on vec[startindex+0,+1 ...,+(vecdim-1)]   vec = b*vec + a * z for a real numbers a and b  
  int SOCIPBlock::vecgetsaz(Matrix& vec,Integer startindex,Real a,bool add)
  {
    assert(vec.dim()>=startindex+vecdim);
    if (add==false){
      mat_xeya(vecdim,vec.get_store()+startindex,z.get_store(),a);
    }
    else {
      mat_xpeya(vecdim,vec.get_store()+startindex,z.get_store(),a);
    }
    return 0;
  }


  int SOCIPBlock::get_mu_info(Integer& inmudim,
			      Real& tr_xz,
			      Real& tr_xdzpdxz,
			      Real& tr_dxdz,
			      Real& min_xz,
			      Real& max_xz) const
  {
    assert((dx.dim()==vecdim)&&(dx.coldim()==1));

    inmudim+=1;
    Real cip=ip(x,z);
    if (cb_out(2))
      get_out()<<" socxz["<<cip<<"]";
    if (cip<min_xz)
      min_xz=cip;
    if (cip>max_xz)
      max_xz=cip;

    tr_xz+=cip;
    tr_xdzpdxz+=ip(x,dz)+ip(z,dx);
    tr_dxdz+=ip(dx,dz);
    

    return 0;
  }
  
  Real SOCIPBlocktestsoc(Matrix tmp)
  {
    Real d= sqr(tmp(0))-mat_ip(tmp.dim()-1,tmp.get_store()+1);
    assert(d>eps_Real);
    return d;
  }


  int SOCIPBlock::get_nbh_info(Integer inmudim,
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
			       Real& ip_dxdz_xdzpdxz) const
  {
    assert(SOCIPBlocktestsoc(x)>eps_Real);
    assert(SOCIPBlocktestsoc(z)>eps_Real);

    const Real mu_xz=tr_xz/inmudim;
    const Real mu_xdzpdxz=tr_xdzpdxz/inmudim;
    const Real mu_dxdz=tr_dxdz/inmudim;

    tmpvec.init(dx);
    apply_F(tmpvec);
    tmpmat.init(dz);
    apply_Finv(tmpmat);
    tmp_xdzpdxz=tmpmat;
    apply_Arw(tmpvec,tmpmat);
    tmpmat(0)-=mu_dxdz;
    //tmpmat now holds dxdz-tr_dxdz/mudim*e_0
    apply_Arw(scaled_point,tmp_xdzpdxz);
    apply_Arw(scaled_point,tmpvec);
    tmp_xdzpdxz+=tmpvec;
    tmp_xdzpdxz(0)-=mu_xdzpdxz;
    //tmp_xdzpdxz now holds xdzpdxz -tr_xdzpdxz/mudim*e_0 

    tmpvec=scaled_point;
    apply_Arw(scaled_point,tmpvec);
    tmpvec(0)-=mu_xz;
    //tmpvec now holds xz-tr_xz/mudim*e_0

    const Real* xzp=tmpvec.get_store();
    const Real* xdzpdxzp=tmp_xdzpdxz.get_store();
    const Real* dxdzp=tmpmat.get_store();
    Real soc_nrmsqr_xz=0;
    Real soc_nrmsqr_xdzpdxz=0;
    Real soc_nrmsqr_dxdz=0;
    Real soc_ip_xz_xdzpdxz=0;
    Real soc_ip_xz_dxdz=0;
    Real soc_ip_dxdz_xdzpdxz=0;
    for(Integer i=0;i<vecdim;i++){
      Real xz=(*xzp++);
      Real xdzpdxz=(*xdzpdxzp++);
      Real dxdz=(*dxdzp++);
      soc_nrmsqr_xz+=sqr(xz);
      soc_nrmsqr_xdzpdxz+=sqr(xdzpdxz);
      soc_nrmsqr_dxdz+=sqr(dxdz);
      soc_ip_xz_xdzpdxz+=xz*xdzpdxz;
      soc_ip_xz_dxdz+=xz*dxdz;
      soc_ip_dxdz_xdzpdxz+=dxdz*xdzpdxz;
    }

    nrmsqr_xz+=soc_nrmsqr_xz;
    nrmsqr_xdzpdxz+=soc_nrmsqr_xdzpdxz;
    nrmsqr_dxdz+=soc_nrmsqr_dxdz;
    ip_xz_xdzpdxz+=soc_ip_xz_xdzpdxz;
    ip_xz_dxdz+=soc_ip_xz_dxdz;
    ip_dxdz_xdzpdxz+=soc_ip_dxdz_xdzpdxz;

    int err=0;
    err=linesearch(alpha);
    if (err){
      if (cb_out()){
	get_out()<<"*** ERROR PSCIPBlock::get_nbh_info(): linesearch(.) returned "<<err<<std::endl;
      }
    }
    
    if (alpha>1000*eps_Real) {
      err=control_nbh_step(alpha,max_nbh,nbh_ubnd,
			   mu_xz,mu_xdzpdxz,mu_dxdz,
			   soc_nrmsqr_xz,soc_nrmsqr_xdzpdxz,soc_nrmsqr_dxdz,
			   soc_ip_xz_xdzpdxz,soc_ip_xz_dxdz,soc_ip_dxdz_xdzpdxz);
      if (err){
	if (cb_out()){
	  get_out()<<"*** ERROR PSCIPBlock::get_nbh_info(): control_nbh_step() returned "<<err<<std::endl;
	}
      }
      
      // TEST BEGIN
      // Matrix e1(vecdim,1,0.);
      // for(Integer i=0;i<=50;i++){
      // 	Real d=i/50.;
      // 	Matrix tmpx(x+d*dx);
      // 	Matrix tmpz(z+d*dz);
      // 	e1(0)=mu_xz+d*mu_xdzpdxz+d*d*mu_dxdz;
      // 	Real q=std::sqrt(soc_nrmsqr_xz+d*(2*soc_ip_xz_xdzpdxz+d*(soc_nrmsqr_xdzpdxz+2*soc_ip_xz_dxdz+d*(2*soc_ip_dxdz_xdzpdxz+d*soc_nrmsqr_dxdz))));
      // 	//if (std::fabs(norm2(apply_Arw(apply_F(tmpx),apply_Finv(tmpz))-e1)-q)>1e-8*mu_xz){
      // 	std::cout<<" SOCnbh: al="<<alpha<<" d="<<d<<" nrm="<<norm2(apply_Arw(apply_F(tmpx),apply_Finv(tmpz))-e1)<<" q="<<q<<" th*mu="<<e1(0)*nbh_ubnd<< std::endl;
      // 	  //}
      // }
      assert(SOCIPBlocktestsoc(x+alpha*dx)>eps_Real);
      assert(SOCIPBlocktestsoc(z+alpha*dz)>eps_Real);
      // TEST END
    }
        
    return err;
  }


  /// if necessary, reduce alpha to the biggest value so that feasibility is maintained with this step size
  int SOCIPBlock::linesearch(Real& alpha) const
  {
    assert((dx.dim()==vecdim)&&(dx.coldim()==1));

    const Real *xp=x.get_store();
    const Real *dxp=dx.get_store();
    Real a=0;
    Real b=0;
    Real c=0;
    xp++;
    dxp++;
    for(Integer j=vecdim-1;--j>=0;xp++,dxp++){
      a-=(*dxp)*(*dxp);
      b-=(*dxp)*(*xp);
      c-=(*xp)*(*xp);
    }
    xp-=vecdim;
    dxp-=vecdim;
    a+=(*dxp)*(*dxp);
    b+=(*dxp)*(*xp);
    c+=(*xp)*(*xp);
    b*=2.;
    assert(c>0.);
    c*=(1.-1.e-8);

    if (alpha*(std::fabs(a)+std::fabs(b))>c){
      if (a<-eps_Real*std::fabs(b)){
	b/=a;c/=a;
	Real d=b*b/4-c;
	assert(d>0.);
	Real f=.95*(-b/2.+std::sqrt(d));
	assert(f>0);
	if (f<alpha){
	  alpha=f;
	}
      }
      else if (a<eps_Real*std::fabs(b)){
	if (b<0.){
	  Real f=.95*(-c/b);
	  assert(f>0.);
	  if (f<alpha){
	    alpha=f;
	  }
	}
      } else {
	b/=a;c/=a;
	Real d=b*b/4.-c;
	if ((b<0.)&&(d>0.)){
	  Real f=.95*(-b/2.-std::sqrt(d));
	  assert(f>0.);
	  if (f<alpha){
	    alpha=f;
	  }
	}
      }
    }

    // //BEGIN TESTING
    // {
    //   xp=x.get_store();
    //   dxp=dx.get_store();
    //   assert((*xp)+alpha*(*dxp)>0.);
    //   Real sum=sqr((*xp++)+alpha*(*dxp++));
    //   for(Integer j=vecdim-1;--j>=0;){
    // 	sum-=sqr((*xp++)+alpha*(*dxp++));
    //   }
    //   if (sum<=0.){
    // 	std::cout<<" socxsum="<<sum<<" alpha="<<alpha<<" a="<<a<<" b="<<b*a<<" c="<<c*a<<std::endl;
    //   }
    //   assert(sum>0.);
    //   std::cout<<" socx("<<alpha<<")ok";
    // }
    // //END TESTING
    
       
    //dual variables
    xp=z.get_store();
    dxp=dz.get_store();
    a=0;
    b=0;
    c=0;
    xp++;
    dxp++;
    for(Integer j=vecdim-1;--j>=0;xp++,dxp++){
      a-=(*dxp)*(*dxp);
      b-=(*dxp)*(*xp);
      c-=(*xp)*(*xp);
    }
    xp-=vecdim;
    dxp-=vecdim;
    a+=(*dxp)*(*dxp);
    b+=(*dxp)*(*xp);
    c+=(*xp)*(*xp);
    b*=2.;
    assert(c>0.);
    c*=(1.-1.e-8);

    if (alpha*(std::fabs(a)+std::fabs(b))>c){
      if (a<-eps_Real*std::fabs(b)){
	b/=a;c/=a;
	Real d=b*b/4-c;
	assert(d>0.);
	Real f=.95*(-b/2.+std::sqrt(d));
	assert(f>0);
	if (f<alpha){
	  alpha=f;
	}
      }
      else if (a<eps_Real*std::fabs(b)){
	if (b<0.){
	  Real f=.95*(-c/b);
	  assert(f>0.);
	  if (f<alpha){
	    alpha=f;
	  }
	}
      } else {
	b/=a;c/=a;
	Real d=b*b/4.-c;
	if ((b<0.)&&(d>0.)){
	  Real f=.95*(-b/2.-std::sqrt(d));
	  assert(f>0.);
	  if (f<alpha){
	    alpha=f;
	  }
	}
      }
    }
         
    // //BEGIN TESTING
    // {
    //   xp=z.get_store();
    //   dxp=dz.get_store();
    //   assert((*xp)+alpha*(*dxp)>0.);
    //   Real sum=sqr((*xp++)+alpha*(*dxp++));
    //   for(Integer j=vecdim-1;--j>=0;){
    // 	sum-=sqr((*xp++)+alpha*(*dxp++));
    //   }
    //   if (sum<=0.){
    // 	std::cout<<" socxsum="<<sum<<" alpha="<<alpha<<" a="<<a<<" b="<<b<<" c="<<c<<std::endl;
    //   }
    //   assert(sum>0.);
    //   std::cout<<" socz("<<alpha<<")ok";
    // }
    // //END TESTING

    return 0;
  }

  /// compute the complementarity_rhs=rhsmu*xi-rhscorr*xi*dx*dz (wihtout "-z") for mu=rhsmu and for corrector for factor rhscorr>0., store this and add it to rhs 
   int SOCIPBlock::add_muxinv(Matrix& rhs,
			    Integer startindex,
			    Real rhsmu,
			    Real rhscorr,
			    bool minus)
  {
    assert(rhs.dim()>=startindex+vecdim);
    assert(rhscorr>=0.);
    assert(rhsmu>=0.);

    if (f.dim()!=vecdim){
      compute_NTscaling();
    }
    
    last_rhs_mu=rhsmu;    
    
    if (rhscorr>0.){
      assert((dx.dim()==vecdim)&&(dx.coldim()==1));
      
      tmpvec.init(dx);
      apply_F(tmpvec);
      compl_rhs.init(dz,-rhscorr);
      apply_Arw(tmpvec,apply_Finv(compl_rhs)); 
    }
    else 
      compl_rhs.init(vecdim,1,0.);

    if (rhsmu>0.){
      compl_rhs(0)+=rhsmu;
    }
    
    if ((rhsmu>0)||(rhscorr>0)){
      apply_Arwinv(scaled_point,compl_rhs);
      apply_F(compl_rhs);
      const Real* cp=tmpvec.get_store();
      const Real* const cend=cp+vecdim;
      Real* vp=rhs.get_store()+startindex;
      if (minus){
	while(cp!=cend)
	  (*vp++)-=(*cp++);
      }
      else {
	while(cp!=cend)
	  (*vp++)+=(*cp++);
      }
    }

    return 0;
  }
  
  /// extract dx from rhs at startindex and compute at the same time dz (=-sys dx -z +complentarity_rhs); 
  int SOCIPBlock::set_dx(const CH_Matrix_Classes::Matrix& rhs,CH_Matrix_Classes::Integer startindex)
  {
    assert((f.dim()==vecdim)&&(f.coldim()==1));
    assert(rhs.dim()>=startindex+vecdim);
    assert((compl_rhs.dim()==vecdim)&&(compl_rhs.coldim()==1));

    dx.init(vecdim,1,rhs.get_store()+startindex);

    dz.init(dx,-1);
    apply_Fsqr(dz);
    dz-=z;
    dz+=compl_rhs;

    //assert(norm2(apply_F(tmpmat.init(dx))+apply_Finv(tmpvec.init(z+dz-compl_rhs)))<1e-8);
    return 0;
  }


  /// compute dx=sysinv*rhs and at the same time dz (=-rhs -z+complentarity_rhs); 
  int SOCIPBlock::set_dx_xizsolverhs(const Matrix& rhs,Integer startindex)
  {
    assert((f.dim()==vecdim)&&(f.coldim()==1));
    assert(rhs.dim()>=startindex+vecdim);
    assert((compl_rhs.dim()==vecdim)&&(compl_rhs.coldim()==1));
    
    dz.init(compl_rhs);
    mat_xmey(vecdim,dz.get_store(),rhs.get_store()+startindex);
    dz-=z;

    dx.init(vecdim,1,rhs.get_store()+startindex);
    apply_Finvsqr(dx);

    //assert(norm2(apply_F(tmpmat.init(dx))+apply_Finv(tmpvec.init(z+dz-compl_rhs)))<1e-8);
    return 0;
  }

  /// compute sysinv*rhs into rhs 
  int SOCIPBlock::apply_xizinv(Matrix& rhs,Integer startindex,bool minus)
  {
    if (f.dim()!=vecdim){
      compute_NTscaling();
    }
    assert((f.dim()==vecdim)&&(f.coldim()==1));
    assert(rhs.dim()>=startindex+vecdim);

    apply_Finvsqr(rhs.get_store()+startindex,minus);

    return 0;
  }

  /// compute sys*rhs into rhs 
  int SOCIPBlock::apply_xiz(Matrix& rhs,Integer startindex,bool minus)
  {
    assert(rhs.dim()>=startindex+vecdim);

    if (f.dim()!=vecdim){
      compute_NTscaling();
    }

    apply_Fsqr(rhs.get_store()+startindex,minus);

    return 0;
  }

  /// move to (x+alpha*dx, z+alpha*dz)
   int SOCIPBlock::do_step(Real alpha)
  {
    assert((dx.dim()==vecdim)&&(dx.coldim()==1));

    if ((old_mu==0.)||(last_rhs_mu<mu)){
      oldx=x;
      oldz=z;
      old_mu=mu;
      last_alpha=alpha;
    }

    mu=last_rhs_mu;
    
    x.xpeya(dx,alpha);
    z.xpeya(dz,alpha);

    // //BEGIN TESTING
    // {
    //   Real *xp=x.get_store();
    //   assert((*xp)>0.);
    //   Real sum=sqr(*xp++);
    //   for(Integer j=vecdim-1;--j>=0;){
    // 	sum-=sqr(*xp++);
    //   }
    //   assert(sum>0.);
    //   std::cout<<" stepsocx("<<alpha<<")ok";
    // }
    // {
    //   Real *zp=z.get_store();
    //   assert((*zp)>0.);
    //   Real sum=sqr(*zp++);
    //   for(Integer j=vecdim-1;--j>=0;){
    // 	sum-=sqr(*zp++);
    //   }
    //   assert(sum>0.);
    //   std::cout<<" stepsocz("<<alpha<<")ok";
    // }
    // //END TESTING


    point_changed();

    return 0;
  }

  /// add the Schur complement to a big system matrix
  int SOCIPBlock::add_AxizinvAt(const Matrix& A,
				Symmatrix& globalsys,
				bool minus,
				bool Atrans)
  {
    assert((Atrans==false?A.coldim():A.rowdim())==vecdim);
    assert((Atrans==false?A.rowdim():A.coldim())==globalsys.rowdim());
    assert((min(x)>0.)&&(min(z)>0.));

    if (f.dim()!=vecdim){
      compute_NTscaling();
    }

    tmpmat.init(A,1.,Atrans==false); //transpose
    apply_Finv(tmpmat);
    
    rankadd(tmpmat,globalsys,minus?-1.:1.,1.,1);

    return 0;
  }

  int SOCIPBlock::add_xiz(Symmatrix& globalsys,Integer startindex,bool minus)
  {

    if (f.dim()!=vecdim){
      compute_NTscaling();
    }

    //add F^2=omega^2*[2*f0^2-1, 2*f0*barf'; 2*f0*barf, I+2*barf*barf']
    const Real* fpi=f.get_store();
    Real o2=omega*omega;
    Real *sysm=globalsys.get_store()+startindex*globalsys.rowdim()-(startindex*(startindex-1))/2;
    
    if (minus) {
      //globalsys(startindex,startindex)+=2*o2;
      (*sysm)+=2*o2;
      Integer endindex=startindex+vecdim;
      for (Integer i=startindex;i<endindex;i++,fpi++){
	Real d=2.*o2*(*fpi);
	const Real* fpj=fpi;
	//globalsys(i,i)-=o2+d*(*fpj++);
	(*sysm++)-=o2+d*(*fpj++);
	for(Integer j=i+1;j<endindex;j++){
	  //globalsys(i,j)-=d*(*fpj++);
	  (*sysm++)-=d*(*fpj++);
	}
	sysm+=globalsys.rowdim()-endindex;
      }
    }
    else {
      //globalsys(startindex,startindex)-=2*o2;
      (*sysm)-=2*o2;
      Integer endindex=startindex+vecdim;
      for (Integer i=startindex;i<endindex;i++,fpi++){
	Real d=2.*o2*(*fpi);
	const Real* fpj=fpi;
	//globalsys(i,i)+=o2+d*(*fpj++);
	(*sysm++)+=o2+d*(*fpj++);
	for(Integer j=i+1;j<endindex;j++){
	  //globalsys(i,j)+=d*(*fpj++);
	  (*sysm++)+=d*(*fpj++);
	}
	sysm+=globalsys.rowdim()-endindex;
      }
    }

    return 0;
  }


  
}
