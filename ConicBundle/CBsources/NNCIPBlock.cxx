/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/NNCIPBlock.cxx
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


#include "NNCIPBlock.hxx"

using namespace CH_Matrix_Classes;

namespace ConicBundle {

  void NNCIPBlock::point_changed()
  {
    dx.init(0,1,0.);
    dz.init(0,1,0.);
    xiz.init(0,1,0.);
    compl_rhs.init(0,1,0.);
  }
  
  void NNCIPBlock::clear(Integer dim)
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
    xiz.newsize(vecdim,Integer(1));
    compl_rhs.newsize(vecdim,Integer(1));
    point_changed();
  }
  
  NNCIPBlock::NNCIPBlock(Integer dim,CBout* cb,int cbinc):InteriorPointBlock(cb,cbinc)
  {clear(dim);}

  NNCIPBlock::~NNCIPBlock()
  {}
  
  Integer NNCIPBlock::get_vecdim() const
  {return vecdim;}
  
  /// set x to value*"one" to x, or if add==true, add value*"one" to x
  int NNCIPBlock::center_x(Real val,bool add)
  {
    point_changed();
    if (add) {
      x+=val;
    }
    else {
      x.init(vecdim,1,val);
    }
    return 0;
  }

  /// set z to value*"one" to z, or if add==true, add value*"one" to z
  int NNCIPBlock::center_z(Real val,bool add)
  {
    point_changed();
    if (add) {
      z+=val;
    }
    else {
      z.init(vecdim,1,val);
    }
    return 0;
  }

  int NNCIPBlock::set_x(const Matrix& vec,Integer startindex,Real& add_center_value)
  {
    assert(vec.dim()>=startindex+vecdim);
    point_changed();
    Real minval=0;
    const Real *vp=vec.get_store()+startindex;
    Real *xp=x.get_store();
    const Real * const xend = xp+vecdim;
    while(xp!=xend){
      Real d=*vp++;
      (*xp++)=d;
      if (d<minval)
	minval=d;
    }
    if (minval<0.)
      minval*=-1;
    add_center_value=minval;
    return 0;
  }

  int NNCIPBlock::set_z(const Matrix& vec,Integer startindex,Real& add_center_value)
  {
    assert(vec.dim()>=startindex+vecdim);
    point_changed();
    Real minval=0;
    const Real *vp=vec.get_store()+startindex;
    Real *zp=z.get_store();
    const Real * const zend = zp+vecdim;
    while(zp!=zend){
      Real d=*vp++;
      (*zp++)=d;
      if (d<minval)
	minval=d;
    }
    if (minval<0.)
      minval*=-1;
    add_center_value=minval;
    return 0;
  }

  /// on vec[startindex+0,+1 ...,+(vecdim-1)]   vec = b*vec + a * x for a real numbers a and b  
  int NNCIPBlock::vecgetsax(Matrix& vec,Integer startindex,Real a,bool add)
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
  int NNCIPBlock::vecgetsaz(Matrix& vec,Integer startindex,Real a,bool add)
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

  /// add dimensions of the primal-dual pairs to mudim and add the inner products of the primal-dual pairs of the current point to current_ip and those of the next point obtained by the given stepsize with the latest computed step to step_ip
  int NNCIPBlock::get_mu_info(Integer& inmudim,
			      Real& tr_xz,
			      Real& tr_xdzpdxz,
			      Real& tr_dxdz,
			      Real& min_xz,
			      Real& max_xz) const
  {
    assert((dx.dim()==vecdim)&&(dx.coldim()==1));

    inmudim+=vecdim;
    Real minv,maxv;
    tr_xz+=ip_min_max(x,z,minv,maxv);
    tr_xdzpdxz+=ip(x,dz)+ip(z,dx);
    tr_dxdz+=ip(dx,dz);

    if (minv<min_xz)
      min_xz=minv;
    if (maxv>max_xz)
      max_xz=maxv;

    if (cb_out(2))
      get_out()<<" linxz["<<minv<<","<<maxv<<"]";
    
    return 0;
  }
  
  /// for limiting the stepsize with respect to the neighborhood this information about x_i*z_i and dx_i*dz_i is required, each block *adds* its contribution to the numbers
  int NNCIPBlock::get_nbh_info(Integer inmudim,
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
    assert(alpha>0.);
    //std::cout<<" NNCcompl="<<norm2(x%z+x%dz+dx%z-last_rhs_mu);
    //std::cout<<" "<<std::fabs(ip(x,z)+ip(dx,z)+ip(x,dz)-last_rhs_mu*vecdim);
    //std::cout<<" lmu="<<last_rhs_mu<<std::endl;
    
    const Real* xp=x.get_store();
    const Real* zp=z.get_store();
    const Real* dxp=dx.get_store();
    const Real* dzp=dz.get_store();
    const Real mu_xz=tr_xz/inmudim;
    const Real mu_xdzpdxz=tr_xdzpdxz/inmudim;
    const Real mu_dxdz=tr_dxdz/inmudim;
    const Real mu_at_one=mu_xz+mu_xdzpdxz+mu_dxdz;
    for(Integer i=0;i<vecdim;i++){
      NNC_nbh_stepsize(*xp++,*zp++,*dxp++,*dzp++,
		       mu_xz,mu_xdzpdxz,mu_dxdz,mu_at_one,
		       nbh_ubnd,alpha,max_nbh,
		       nrmsqr_xz,nrmsqr_xdzpdxz,nrmsqr_dxdz,
		       ip_xz_xdzpdxz,ip_xz_dxdz,ip_dxdz_xdzpdxz);
    }
    return 0;
  }


  
  /// if necessary, reduce alpha to the biggest value so that feasibility is maintained with this step size
  int NNCIPBlock::linesearch(Real& alpha) const
  {
    assert((dx.dim()==vecdim)&&(dx.coldim()==1));

    const Real* vp=x.get_store();
    const Real* dvp=dx.get_store();
    const Real* const vendx=vp+vecdim;
    Real a=alpha;
    while(vp!=vendx){
      Real d=*dvp++;
      if (d<0.){
	d=-(*vp++)/d;
	if (d<a)
	  a=d;
      }
      else {
	vp++;
      }
    }
    vp=z.get_store();
    dvp=dz.get_store();
    const Real* const vendz=vp+vecdim;
    while(vp!=vendz){
      Real d=*dvp++;
      if (d<0.){
	d=-(*vp++)/d;
	if (d<a)
	  a=d;
      }
      else {
	vp++;
      }
    }
    if (a<alpha)
      alpha=a;
    return 0;
  }

  /// compute the complementarity_rhs=rhsmu*xi-rhscorr*xi*dx*dz (wihtout "-z") for mu=rhsmu and for corrector for factor rhscorr>0., store this and add it to rhs 
   int NNCIPBlock::add_muxinv(Matrix& rhs,
			    Integer startindex,
			    Real rhsmu,
			    Real rhscorr,
			    bool minus)
  {
    assert(rhs.dim()>=startindex+vecdim);
    assert(rhscorr>=0.);
    assert(rhsmu>=0.);

    //make sure, the system is available
    if (xiz.dim()!=vecdim){
      compute_NTscaling();
    }

    last_rhs_mu=rhsmu;
    
    if (rhsmu>0.){
      compl_rhs.newsize(vecdim,1); chk_set_init(compl_rhs,1);
      const Real* xp=x.get_store();
      const Real* const xend=xp+vecdim;
      Real* vp=compl_rhs.get_store();
      while(xp!=xend){
	//assert(*xp>0.);
	*vp++ = rhsmu/(*xp++);
      }
    }
    else 
      compl_rhs.init(vecdim,1,0.);
    
    if (rhscorr>0.){
      assert((dx.dim()==vecdim)&&(dx.coldim()==1));
      const Real* xp=x.get_store();
      const Real* const xend=xp+vecdim;
      const Real* dxp=dx.get_store();
      const Real* dzp=dz.get_store();
      Real* vp=compl_rhs.get_store();
      while(xp!=xend){
	//assert(*xp>0.);
	*vp++ -= rhscorr*(*dxp++)*(*dzp++)/(*xp++);
      }
    }

    // std::cout<<"startind="<<startindex<<std::endl;
    // std::cout<<" mu="<<rhsmu<<std::endl;
    // std::cout<<" rhscorr="<<rhscorr<<std::endl;
    // std::cout<<" compl_rhs="<<compl_rhs;
    // std::cout<<" rhsbefore="<<rhs;

    if ((rhsmu>0)||(rhscorr>0)){
      const Real* cp=compl_rhs.get_store();
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

    // std::cout<<" rhsafter="<<rhs;

    return 0;
  }

  /// extract dx from rhs at startindex and compute at the same time dz (=-sys dx -z +complentarity_rhs); 
  int NNCIPBlock::set_dx(const Matrix& rhs,Integer startindex)
  {
    assert((xiz.dim()==vecdim)&&(xiz.coldim()==1));
    assert(rhs.dim()>=startindex+vecdim);
    assert((compl_rhs.dim()==vecdim)&&(compl_rhs.coldim()==1));

    dx.init(vecdim,1,rhs.get_store()+startindex);
    dz.init(dx,-1.);
    dz%=xiz;
    dz+=compl_rhs;
    dz-=z;

    // std::cout<<" setdx="<<dx;
    // std::cout<<" setdz="<<dz;
    // std::cout<<" setx="<<x;
    // std::cout<<" setz="<<z;
    // std::cout<<" setcompl_rhs="<<compl_rhs;
    // std::cout.flush();

    //assert(norm2(dx%z+dz%x+x%z-x%compl_rhs)<=1e-10*max(norm2(x),norm2(z)));
    return 0;
  }

  
  /// compute dx=sysinv*rhs and at the same time dz (=-rhs-z +complentarity_rhs); 
  int NNCIPBlock::set_dx_xizsolverhs(const Matrix& rhs,Integer startindex)
  {
    assert((xiz.dim()==vecdim)&&(xiz.coldim()==1));
    assert(rhs.dim()>=startindex+vecdim);
    assert((compl_rhs.dim()==vecdim)&&(compl_rhs.coldim()==1));
    
    dz.init(compl_rhs);
    mat_xmey(vecdim,dz.get_store(),rhs.get_store()+startindex);
    dz-=z;

    dx.newsize(vecdim,1); chk_set_init(dx,1);
    Real* dxp=dx.get_store();
    const Real* const xpend=dxp+vecdim;
    const Real* vp=rhs.get_store()+startindex;
    const Real* xizp=xiz.get_store();
    while (dxp!=xpend) {
      //assert(*xizp>0);
      (*dxp++) = (*vp++)/(*xizp++);
    }

    // std::cout<<" sovdx="<<dx;
    // std::cout<<" sovdz="<<dz;
    // std::cout<<" sovx="<<x;
    // std::cout<<" sovz="<<z;
    // std::cout<<" sovcompl_rhs="<<compl_rhs;
    // std::cout.flush();

    //assert(norm2(dx%z+dz%x+x%z-x%compl_rhs)<=1e-8*max(norm2(dx),norm2(dz)));
   return 0;
  }

  /// compute sysinv*rhs into rhs 
  int NNCIPBlock::apply_xizinv(Matrix& rhs,Integer startindex,bool minus)
  {
    assert(rhs.dim()>=startindex+vecdim);
    //make sure, the system is available
    if (xiz.dim()!=vecdim){
      compute_NTscaling();
    }

    Real* vp=rhs.get_store()+startindex;
    const Real* const vend=vp+vecdim;
    const Real* xizp=xiz.get_store();
    if (minus){
      while (vp!=vend) {
	//assert(*xizp>0);
	(*vp++) /= -(*xizp++);
      }
    }
    else{
      while (vp!=vend) {
	//assert(*xizp>0);
	(*vp++) /= (*xizp++);
      }
    }

    return 0;
  }

  /// compute sys*rhs into rhs 
  int NNCIPBlock::apply_xiz(Matrix& rhs,Integer startindex,bool minus)
  {
    assert(rhs.dim()>=startindex+vecdim);

    //make sure, the system is available
    if (xiz.dim()!=vecdim){
      compute_NTscaling();
    }

    Real* vp=rhs.get_store()+startindex;
    const Real* const vend=vp+vecdim;
    const Real* xizp=xiz.get_store();
    if (minus){
      while (vp!=vend) {
	(*vp++) *= -(*xizp++);
      }
    }
    else{
      while (vp!=vend) {
	(*vp++) *= (*xizp++);
      }
    }

    return 0;
  }

  /// move to (x+alpha*dx, z+alpha*dz)
   int NNCIPBlock::do_step(Real alpha)
  {
    assert((dx.dim()==vecdim)&&(dx.coldim()==1));

    if ((old_mu==0.)||(last_rhs_mu<mu)){
      oldx=x;
      oldz=z;
      old_mu=mu;
      last_alpha=alpha;
    }

    old_mu=mu;
    mu=last_rhs_mu;
    
    x.xpeya(dx,alpha);
    z.xpeya(dz,alpha);

    point_changed();

    return 0;
  }

  /// add the Schur complement to a big system matrix
  int NNCIPBlock::add_AxizinvAt(const Matrix& A,
				Symmatrix& globalsys,
				bool minus,
				bool Atrans)
  {
    assert((Atrans==false?A.coldim():A.rowdim())==vecdim);
    assert((Atrans==false?A.rowdim():A.coldim())==globalsys.rowdim());
    //assert((min(x)>0.)&&(min(z)>0.));

    if (xiz.dim()!=vecdim){
      compute_NTscaling();
	      
    }

    tmpvec(xiz);
    tmpvec.sqrt();
    tmpvec.inv();

    scaledrankadd(A,tmpvec,globalsys,minus?-1:1.,1.,Atrans);

    return 0;
  }

  /// add the system matrix*factor into a big system matrix starting at startindex
  int NNCIPBlock::add_xiz(Symmatrix& globalsys,Integer startindex,bool minus)
  {
    //assert((min(x)>0.)&&(min(z)>0.));

    if (xiz.dim()!=vecdim){
      compute_NTscaling();
    }

    //TEST output start
    // std::cout<<"x="<<x;
    // std::cout<<"z="<<z;
    // std::cout<<"xiz="<<xiz;

    const Real* xizp=xiz.get_store();
    if (minus) {
      for (Integer i=startindex;i<startindex+vecdim;i++){
	globalsys(i,i)-=*xizp++;
      }
    }
    else {
      for (Integer i=startindex;i<startindex+vecdim;i++){
	globalsys(i,i)+=*xizp++;
      }
    }

    return 0;
  }



}
