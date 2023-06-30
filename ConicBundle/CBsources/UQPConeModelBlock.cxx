/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/UQPConeModelBlock.cxx
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



#include "mymath.hxx"
#include "UQPConeModelBlock.hxx"

 
using namespace CH_Matrix_Classes;

namespace ConicBundle {
  
// ****************************************************************************
//                 routines for computing SOC directions
// ****************************************************************************

// x0   = x(0);        (must be >0.)
// barx = x(1..n-1);
//
// Arw(x)=[ x0, barx'; barx, x0*I ]
//
// gamma = sqrt(x0*x0 - barx'*barx);  (must be >0.)
//
// Arwinv(x) = (1/gamma^2)*[ x0, -barx'; -barx', (barx*barx'+gamma^2*I)/x0]:
// 
// G(x)= [x0, barx'; barx, gamma*I+barx*barx'/(gamma+x0)];
//
// Ginv(x)= [-1, 0;0, I]/gamma  -  1/( gamma*( barx'*barx / (gamma*x0)-x0)) *
//          [gamma+x0, -barx'; barx, barx*barx'/(gamma+x0)];
//
// omega=sqrt(gamma(z)/gamma(x))
//
// F(f)= omega*[f0, barf'; barf, I+barf*barf'/(1+f0)]; (gamma(f)=1)
//
// Finv(f)= (1/omega)*[f0, -barf'; -barf, I+barf*barf'/(1+f0)];
//


// Arw(x)*v

  Matrix UQPConeModelBlock::mult_Arw(const Matrix& x,const Matrix& v) const
  {
    chk_add(x,v);
    Matrix w(x,v(0));
    mat_xpeya(w.dim()-1,w.get_store()+1,v.get_store()+1,x(0));
    w(0)+=mat_ip(x.rowdim()-1,x.get_store()+1,v.get_store()+1);
    return w;
  }
  
  
// Arwinv(x)*v    
 
  Matrix UQPConeModelBlock::mult_Arwinv(const Matrix& x,const Matrix& v) const
  {
    chk_add(x,v);
    Real x0=x(0);
    Real gamma2=x0*x0-mat_ip(x.rowdim()-1,x.get_store()+1);

    //assert(x0>0.);
    //assert(gamma>0.);
    //assert(fabs(gamma-sqrt(x0*x0-mat_ip(x.rowdim()-1,x.get_store()+1,x.get_store()+1)))<1e-10*gamma);
    
    Real ip_barxv=mat_ip(x.rowdim()-1,x.get_store()+1,v.get_store()+1);
    Matrix w(x.rowdim(),1); chk_set_init(w,1);
    w(0)=(v(0)*x(0)-ip_barxv)/gamma2;
    mat_xeya(w.dim()-1,w.get_store()+1,x.get_store()+1,(-v(0)+ip_barxv/x0)/gamma2);
    mat_xpeya(w.dim()-1,w.get_store()+1,v.get_store()+1,1./x0);
    return w;
  }

//w=G(x)*v

  Matrix UQPConeModelBlock::mult_G(const Matrix& x,Real gamma,const Matrix& v) const
  {
    chk_add(x,v);
    Real x0=x(0);

    //assert(x0>0.);
    //assert(gamma>0.);
    //assert(fabs(gamma-sqrt(x0*x0-mat_ip(x.rowdim()-1,x.get_store()+1,x.get_store()+1)))<1e-10*gamma);
    
    
    Real ip_barxv=mat_ip(x.rowdim()-1,x.get_store()+1,v.get_store()+1);
    Matrix w(x,v(0));
    w(0)+=ip_barxv;
    mat_xpeya(w.dim()-1,w.get_store()+1,x.get_store()+1,ip_barxv/(x0+gamma));
    mat_xpeya(w.dim()-1,w.get_store()+1,v.get_store()+1,gamma);
    return w;
  }
   
//w=Ginv(x)*v
   
  Matrix UQPConeModelBlock::mult_Ginv(const Matrix& x,Real gamma,const Matrix& v) const
  {
    chk_add(x,v);
    Real x0=x(0);

    //assert(x0>0.);
    //assert(gamma>0.);
    //assert(fabs(gamma-sqrt(x0*x0-mat_ip(x.rowdim()-1,x.get_store()+1,x.get_store()+1)))<1e-10*gamma);
    
    Real ip_barxv=mat_ip(x.rowdim()-1,x.get_store()+1,v.get_store()+1);
    Real d=-1./(gamma*((x0*x0-gamma*gamma)/(gamma+x0)-x0));
    Matrix w(x.rowdim(),1); chk_set_init(w,1);
    w(0)=v(0)*(-1./gamma+(gamma+x0)*d)-ip_barxv*d;
    mat_xeya(w.dim()-1,w.get_store()+1,v.get_store()+1,1./gamma);
    mat_xpeya(w.dim()-1,w.get_store()+1,x.get_store()+1,(ip_barxv/(gamma+x0)-v(0))*d);
    return w;
  }
    
// w=F(f)*v    
 

// ****************************************************************************
//                           clear
// ****************************************************************************

void UQPConeModelBlock::clear()
{
  UQPModelBlock::clear();
  x.init(0,1,0.);         
  z.init(0,1,0.);
  y=0.;
  A.init(0,0,0.);
  b=0.;              
  less_or_equal=false;
  s=0.;
  scost=0.;
  nnc_dim=0;             
  soc_dim.init(0,1,Integer(0));
  soc_start.init(0,1,Integer(0));       
  Xp.clear();
  Zp.clear();
  Xinv.clear();
  lb.init(0,1,0.);
  ub.init(0,1,0.);
  uz.init(0,1,0.);
  box_uz=0;
  box_dim=0;
  box_start=0;

  qp_xstart=0;
  qp_ystart=0;
  mu_dim=0;
  restart_factor=0.;


  current_mu=0.;
  old_mu=0.;
  old_x.init(0,1,0.);                   
  old_z.init(0,1,0.);
  old_uz.init(0,1,0.);
  old_box_uz=0;
  old_y=0.;                     
  old_s=0.;   
  old_Xp.clear(); 
  old_Zp.clear(); 
  last_alpha=0.;

  zcorr.init(0,1,0.);
  uzcorr.init(0,1,0.);
  box_uzcorr=0;

  dx.init(0,1,0.);                   
  dy=0.;                   
  dz.init(0,1,0.);
  duz.init(0,1,0.);
  box_duz=0.;
  ds=0.;
    
  tmpvec.init(0,0,0.);
  tmpsvec.init(0,0,0.);
  tmpsymmat.init(0,0.);
  dX.init(0,0.);
}

// ****************************************************************************
//                           init
// ****************************************************************************

  int UQPConeModelBlock::init(const MinorantPointer& in_constant_minorant,
			  const MinorantBundle& in_bundle,
			  Integer in_nnc_dim,
			  const Indexmatrix& in_soc_dim,
			  const Indexmatrix& sdp_dim,
			  const Matrix& box_lb,
			  const Matrix& box_ub,
			  Real in_b,
			  FunctionTask ft,
			      QPModelOracleDataObject* /* oracle_data */,
			  bool scale_box )
{
  int err=0;
  
  nnc_dim=in_nnc_dim;
  if (nnc_dim<0){
    if (cb_out())
      get_out()<<"**** ERROR in UQPConeModelBlock::init(.......): linear model dimension negative: nnc_dim="<<nnc_dim<<std::endl;
    nnc_dim=0;
    err++;
  }
  soc_dim=in_soc_dim;
  b=in_b;
  if (b<0){
    if (cb_out())
      get_out()<<"**** ERROR in UQPConeModelBlock::init(.......): trace rhs value negative: "<<b<<std::endl;
    b=1.;
    err++;
  }
  box_dim=0;
  if (box_lb.dim()>0){
    assert(box_lb.coldim()==1);
    assert(box_ub.coldim()==1);
    if ((box_lb.rowdim()!=box_ub.rowdim())||(min(box_ub-box_lb)<1e-10)){
      if (cb_out())
	get_out()<<"**** ERROR in UQPConeModelBlock::init(.......): box is incocsistent or has no interior; box_lb="<<transpose(box_lb)<<" box_ub="<<transpose(box_ub)<<std::endl;
      err++;
    }
    if ((scale_box==false)&&(ft!=ObjectiveFunction)){
      if (cb_out())
	get_out()<<"**** WARNING: QPModelBlock::init(): scale_box is false but the function task is not ObjectionFunction; setting scale_box to true"<<std::endl;
      scale_box=true;
    }
    box_scaleub=false;
    box_scaling=scale_box;
    if (scale_box){
      if ((nnc_dim==0)&&(sum(soc_dim)==0)&&(sum(sdp_dim)==0)){
	box_scaleub=true;
      }
      lb=box_lb;
      ub=box_ub;
    }
    else {
      if ((nnc_dim==0)&&(sum(soc_dim)==0)&&(sum(sdp_dim)==0)){
	lb.init(box_lb,b);
	ub.init(box_ub,b);
      }
      else {
	lb=box_lb;
	ub=box_ub;
      }
    }
    box_dim=box_lb.rowdim()+(box_scaling?1:0);
  }
    
  scost=0.;
  less_or_equal=((ft!=ObjectiveFunction)&&(!box_scaleub));
  Integer xdim=nnc_dim;
  if (less_or_equal) {
    mu_dim=nnc_dim+1;
  }
  else {
    mu_dim=nnc_dim;
  }
  assert((soc_dim.dim()==0)||(min(soc_dim)>0));
  mu_dim+=soc_dim.dim();
  xdim+=sum(soc_dim);

  Xp.resize((unsigned long)(sdp_dim.dim()));
  Zp.resize((unsigned long)(sdp_dim.dim()));
  old_Xp.resize((unsigned long)(sdp_dim.dim()));
  old_Zp.resize((unsigned long)(sdp_dim.dim()));
  Xinv.resize((unsigned long)(sdp_dim.dim()));

  for(unsigned int i=0;i<Xp.size();i++){
    Integer sdim=sdp_dim(Integer(i));
    if (sdim<0){
      if (cb_out())
	get_out()<<"**** ERROR in UQPConeModelBlock::init(.......): psc cone "<<i<<" has negative dimension "<<sdim<<std::endl;
      sdim=1;
      err++;
    }
    mu_dim+=sdim;
    Xp[i].init(sdim,0.);
    Xp[i].shift_diag(1.);
    xdim+=(sdim*(sdim+1))/2;
  }
  if (box_dim>0){
    mu_dim+=2*lb.rowdim()+(box_scaling?1:0)+(box_scaleub?1:0);
    xdim+=box_dim;
    uz.init(lb.rowdim(),1,0.);
    box_uz=0.;
    duz.init(lb.rowdim(),1,0.);
    box_duz=0.;
    old_uz.init(lb.rowdim(),1,0.);
    old_box_uz=0.;
    uzcorr.init(lb.rowdim(),1,0.);
    box_uzcorr=0.;
  }
  
  x.init(xdim,1,0.);
  z.init(xdim,1,0.);
  dx.init(xdim,1,0.);
  dz.init(xdim,1,0.);
  old_x.init(xdim,1,0.);
  old_z.init(xdim,1,0.);
  zcorr.init(xdim,1,0.);
  last_alpha=0.;
  y=0.;
  dy=0.;
  old_y=0;
  s=0.;
  ds=0;
  old_s=0.;
  scorr=0;

  constant_minorant.clear();
  constant_minorant.push_back(in_constant_minorant);
  bundle.clear();
  bundle.push_back(in_bundle);
  //bundle.clear();
  //bundle.insert(bundle.end(),in_bundle().begin(),in_bundle().end());
  
  if (xdim!=Integer(in_bundle.size())){
    if (cb_out())
      get_out()<<"**** ERROR in UQPConeModelBlock::init(.......): model dimension does not match number of bundle vectors: xdim="<<xdim<<" bundle.size()="<<in_bundle.size()<<std::endl;
    err++;
  }
  
  //--- construct A
  if ((nnc_dim==0)&&(sum(soc_dim)==0)&&(sum(sdp_dim)==0)){
    A.init(0,xdim,1.);
  }
  else {
    A.init(1,xdim,1.);
    Integer tmpdim=nnc_dim;
    soc_start.newsize(soc_dim.dim(),1); chk_set_init(soc_start,1);
    {for(Integer i=0;i<soc_dim.dim();i++){
	soc_start(i)=tmpdim;
	mat_xea(soc_dim(Integer(i))-1,A.get_store()+tmpdim+1,0.);
	tmpdim+=soc_dim(Integer(i));
      }}
    {for(unsigned int i=0;i<Xp.size();i++){
	svec(Xp[i],tmpsvec);
	mat_xey(tmpsvec.dim(),A.get_store()+tmpdim,tmpsvec.get_store());
    tmpdim+=(sdp_dim(Integer(i))*(sdp_dim(Integer(i))+1))/2;
      }}
    if (lb.rowdim()>0){
      box_start=tmpdim;
      mat_xea(lb.rowdim(),A.get_store()+tmpdim,0.);
      tmpdim+=box_dim;
    }
    assert(xdim==tmpdim);
  }


  return err;
}

// ****************************************************************************
//                           adjust_trace
// ****************************************************************************

int UQPConeModelBlock::adjust_trace(Real in_b)
{
  b=in_b;
  return 0;
}

// ****************************************************************************
//                           evaluate_trace
// ****************************************************************************

Real UQPConeModelBlock::evaluate_trace() const
{
  if (A.rowdim()==0)
    return 0.;
  return mat_ip(x.rowdim(),A.get_store(),x.get_store());
}

// ****************************************************************************
//                           get_nncx
// ****************************************************************************

int UQPConeModelBlock::get_nncx(Matrix& nncx,
			    Matrix* nncx_activity,
			    bool /* cautious */)
{
  nncx.init(nnc_dim,1,x.get_store());

  if (nncx_activity){
    if (nnc_dim==0){
      nncx_activity->init(0,1,0.);
      return 0;
    }

    //determine activity
    Real tapia_factor=(old_mu>0.)?current_mu/old_mu:1.;
 
    nncx_activity->newsize(nnc_dim,1); chk_set_init(*nncx_activity,1);
    Matrix dual_tapia(nnc_dim,1); chk_set_init(dual_tapia,1);
    for(Integer i=0;i<nnc_dim;i++){
      (*nncx_activity)(i)=nncx(i)/old_x(i);
      dual_tapia(i)=z(i)/old_z(i);
    }
    Indexmatrix sind;
    sortindex(nncx,sind,false);
    Integer i0=sind(0);

    if (cb_out(2)){
      get_out()<<" last_alpha="<<last_alpha;
      get_out()<<" current_mu="<<current_mu<<" old_mu="<<old_mu<<" tapia_factor="<<tapia_factor<<std::endl;
      for(Integer i=0;i<nnc_dim;i++){
	Integer ii=sind(i);
	get_out()<<" "<<i<<"("<<ii<<"):"<<nncx(ii)<<","<<z(ii);
	get_out()<<" "<<Real(z(i0)*nncx(ii)>min(1e-3,current_mu)*z(ii)*nncx(i0));
	get_out()<<" : "<<old_x(ii)<<","<<old_z(ii);
	get_out()<<"   ("<<(*nncx_activity)(ii)<<","<<dual_tapia(ii)<<")";
	get_out()<<" "<<Real((*nncx_activity)(ii)>dual_tapia(ii))<<" ";
	if (tapia_factor<0.99){
	  get_out()<<"   ("<<((*nncx_activity)(ii)-tapia_factor)/(1.-tapia_factor);
	  get_out()<<","<<(dual_tapia(ii)-tapia_factor)/(1.-tapia_factor)<<") ";
	  get_out()<<" "<<Real(((*nncx_activity)(ii)>0.8)&&(1.1*(((*nncx_activity)(ii)-tapia_factor)/(1.-tapia_factor))>((dual_tapia(ii)-tapia_factor)/(1.-tapia_factor))));
	}
	get_out()<<std::endl;						      
      }
    }

    Integer maxi=0;
    for (Integer i=0;i<nnc_dim;i++){
      Integer ii=sind(i);
      //if (z(i0)/nncx(i0)>min(1e-3,current_mu)*z(ii)/nncx(ii))
      //	   maxi=i;
	
      if (tapia_factor<1-1e-6){
	  if((((*nncx_activity)(ii)>0.8)&&((1.1*(((*nncx_activity)(ii)-tapia_factor)/(1.-tapia_factor))>((dual_tapia(ii)-tapia_factor)/(1.-tapia_factor))))))
	    maxi=i;
      }
      else {
	if (z(i0)/nncx(i0)>current_mu*z(ii)/nncx(ii))
	  maxi=i;
      }
    }
    for(Integer i=0;i<=maxi;i++)
      (*nncx_activity)(sind(i))=1.;
    for(Integer i=maxi+1;i<nnc_dim;i++)
      (*nncx_activity)(sind(i))=0.;

  }

  return 0;
}

// ****************************************************************************
//                           get_socx
// ****************************************************************************

int UQPConeModelBlock::get_socx(Integer soci,
			     Matrix& socx,
			    Real* socx_activity,
			    bool /* cautious */) 
{
  if ((soci<0)||(soci>soc_dim.rowdim())){
    if (cb_out())
      get_out()<<"**** ERROR in UQPConeModelBlock::get_socx(...): SOC variable index "<<soci<<" needs to be in [0,"<<soc_dim.rowdim()<<")"<<std::endl;
    return 1;
  }
  Integer sdim=soc_dim(soci);
  socx.init(sdim,1,x.get_store()+soc_start(soci));

  if (socx_activity){
    if (sdim==0) {
      *socx_activity=0.;
      return 0;
    }
    
    //determine activity
    Real x0=socx(0);
    Real z0=z(soc_start(soci));
    Real oldx0=old_x(soc_start(soci));
    Real oldz0=old_z(soc_start(soci));
   
    Real tapia_factor=(old_mu>0.)?current_mu/old_mu:1.;


    *socx_activity = x0/oldx0;  
    Real dual_tapia= z0/oldz0;

    if (cb_out(2)){
      get_out()<<" last_alpha="<<last_alpha;
      get_out()<<" current_mu="<<current_mu<<" old_mu="<<old_mu<<" tapia_factor="<<tapia_factor<<std::endl;
      get_out()<<" "<<x0<<","<<z0;
      get_out()<<" : "<<oldx0<<","<<oldz0;
      get_out()<<"   ("<<(*socx_activity)<<","<<dual_tapia<<")";
      get_out()<<" "<<Real((*socx_activity)>dual_tapia)<<" ";
      if (tapia_factor<0.99){
	get_out()<<"   ("<<((*socx_activity)-tapia_factor)/(1.-tapia_factor);
	get_out()<<","<<(dual_tapia-tapia_factor)/(1.-tapia_factor)<<") ";
	get_out()<<" "<<Real(((*socx_activity)>0.8)&&(1.1*(((*socx_activity)-tapia_factor)/(1.-tapia_factor))>((dual_tapia-tapia_factor)/(1.-tapia_factor))));
      }
      get_out()<<std::endl;						      
    }

    if (tapia_factor<1-1e-6){
      if(
	 (*socx_activity>0.8)&&
	 ((1.1*((*socx_activity-tapia_factor)/(1.-tapia_factor))>((dual_tapia-tapia_factor)/(1.-tapia_factor))))
	 )
	*socx_activity=1.;
      else
	*socx_activity=0.;
    }
    else {
      if (x0>0.1*std::sqrt(current_mu)*z0)
	*socx_activity=1.;
      else
	*socx_activity=0.;
    }
  }

  return 0;
}

// ****************************************************************************
//                           get_pscx
// ****************************************************************************

int UQPConeModelBlock::get_pscx(Integer psci,
				Matrix& pscx_eigs,
				Matrix& pscx_vecs,
				Real& growthrate,
				Matrix& primalgrowth,
				Matrix& dualgrowth) 
{
  if ((psci<0)||(psci>Integer(Xp.size()))){
    if (cb_out())
      get_out()<<"**** ERROR in UQPConeModelBlock::get_socx(...): SOC variable index "<<psci<<" needs to be in [0,"<<Xp.size()<<")"<<std::endl;
    return 1;
  }

  if (cb_out(2)){
    get_out().precision(12);
    get_out()<<" UQPdualpsc="<<y;  
  }
  
  Integer sdim=Xp[unsigned(psci)].rowdim();
  if (sdim==0){
    pscx_eigs.init(0,1,0.);
    pscx_vecs.init(0,0,0.);
    growthrate=0.;
    primalgrowth.init(0,1,0.);
    dualgrowth.init(0,1,0.);
    //if (pscx_activity_rank)
    //  *pscx_activity_rank=0;
    return 0;
  }

  
  int eigstat=Xp[unsigned(psci)].eig(pscx_vecs,pscx_eigs,false);
  if (eigstat){
    if (cb_out())
      get_out()<<"**** WARNING UQPConeModelBlock::get_pscx(....): eigenvalue factorization for semidefinite Xp["<<psci<<"] failed and returned "<<eigstat<<std::endl;
  }

  growthrate=ip(Xp[unsigned(psci)],Zp[unsigned(psci)])/ip(old_Xp[unsigned(psci)],old_Zp[unsigned(psci)]);
  primalgrowth.newsize(pscx_eigs.rowdim(),1); chk_set_init(primalgrowth,1);
  dualgrowth.newsize(pscx_eigs.rowdim(),1); chk_set_init(dualgrowth,1);
    
  for(Integer i=0;i<pscx_eigs.rowdim();i++){
    if (pscx_eigs(i)<0.){
      if (cb_out()) {
	get_out()<<"*** WARNING: PUQPConeModelBlock::get_pscx(....): pscx_eigs("<<i<<")="<<pscx_eigs(i)<<"<0, setting to 0"<<std::endl;
      }
      pscx_eigs(i)=0.;
    }
    tmpvec=pscx_vecs.col(i);
    Real Xval=pscx_eigs(i);
    Real Zval=ip(tmpvec,Zp[unsigned(psci)]*tmpvec);
    Real oldXval=ip(tmpvec,old_Xp[unsigned(psci)]*tmpvec);
    Real oldZval=ip(tmpvec,old_Zp[unsigned(psci)]*tmpvec);
    primalgrowth(i)=Xval/oldXval;
    dualgrowth(i)=Zval/oldZval;
  }
  
  return 0;
  
  /*
  ///   ------------  old code ---------------  
  Integer i=sdim-1;
  while((i>=0)&&(pscx_eigs(i)<0.)){
    if (cb_out()) {
      get_out()<<"*** WARNING: PUQPConeModelBlock::get_pscx(....): pscx_eigs("<<i<<")="<<pscx_eigs(i)<<"<0, setting to 0"<<std::endl;
    }
    pscx_eigs(i)=0.;
    --i;
  }

  if (pscx_activity_rank){

    Integer rank=sdim;
    Real tapia_factor=(old_mu>0.)?current_mu/old_mu:1.;

    if (cb_out(2)){
      get_out()<<" last_alpha="<<last_alpha;
      get_out()<<" mu="<<current_mu<<" old_mu="<<old_mu;
      get_out()<<" tapia_factor="<<tapia_factor<<std::endl;
      Real X0val=pscx_eigs(0);
      Matrix tmpvec=pscx_vecs.col(0);
      Real Z0val=ip(tmpvec,Zp[unsigned(psci)]*tmpvec);
      for(Integer i=0;i<pscx_eigs.dim();i++){
	tmpvec=pscx_vecs.col(i);
	Real Xval=pscx_eigs(i);
	Real Zval=ip(tmpvec,Zp[unsigned(psci)]*tmpvec);
	Real oldXval=ip(tmpvec,old_Xp[unsigned(psci)]*tmpvec);
	Real oldZval=ip(tmpvec,old_Zp[unsigned(psci)]*tmpvec);
	Real primal_tapia=Xval/oldXval;
	Real dual_tapia=Zval/oldZval;	
	get_out()<<" "<<i<<":"<<Xval<<","<<Zval;
	get_out()<<" "<<Real(Z0val/X0val>min(1e-3,current_mu)*Zval/Xval);
	get_out()<<" : "<<oldXval<<","<<oldZval;
	get_out()<<"   ("<<primal_tapia<<","<<dual_tapia<<")";
	get_out()<<" "<<(primal_tapia>dual_tapia)<<" ";
	get_out()<<"[b"<<(Xval>0.1*b);
	get_out()<<",c"<<(Xval<1e-10*b);
	//get_out()<<",b"<<(Zval<Xval);
	get_out()<<",b"<<((Zval<1e3*Xval)&&(primal_tapia>2.*dual_tapia));
	get_out()<<",c"<<((current_mu<0.01)&&(Xval<1e-2*std::sqrt(current_mu))&&(Zval>std::sqrt(current_mu)));
	get_out()<<",b"<<(((tapia_factor<.8)&&((primal_tapia>0.8)||(primal_tapia>.9*dual_tapia)))||((tapia_factor>=.8)&&((Z0val/X0val>std::sqrt(current_mu)*Zval/Xval))));
	get_out()<<"]";
	// if (tapia_factor<0.99){
	//   get_out()<<"   ("<<(primal_tapia-tapia_factor)/(1.-tapia_factor);
	//   get_out()<<","<<(dual_tapia-tapia_factor)/(1.-tapia_factor)<<") ";
	//   get_out()<<" "<<Real((primal_tapia>0.8)&&(1.1*((primal_tapia-tapia_factor)/(1.-tapia_factor))>((dual_tapia-tapia_factor)/(1.-tapia_factor))));
	// }
	get_out()<<std::endl;						      
      }
    }

    
    Real X0val=pscx_eigs(0);
    if (X0val<1e-10*b)
      return 0;
    Matrix tmpvec=pscx_vecs.col(0);
    Real Z0val=ip(tmpvec,Zp[unsigned(psci)]*tmpvec);
    
    while (--rank>=0) {
      Real Xval=pscx_eigs(rank);
	
      //a rather large eigenvalue might not allow to exclude its eigenspace yet 
      if (Xval>0.1*b)
	break;
     
      //if the value is really small, its eigenspace is irrelevant
      if (Xval<1e-10*b)
	continue;
	
      //check whether the dual signals irrelevance of this eigenspace
      tmpvec=pscx_vecs.col(rank);
      Real Zval=ip(tmpvec,Zp[unsigned(psci)]*tmpvec);
      Real oldXval=ip(tmpvec,old_Xp[unsigned(psci)]*tmpvec);
      Real oldZval=ip(tmpvec,old_Zp[unsigned(psci)]*tmpvec);
      Real primal_tapia=Xval/oldXval;
      Real dual_tapia=Zval/oldZval;
      //if (Zval<Xval)
      if ((Zval<1e-3*Xval)&&(primal_tapia>2.*dual_tapia))
	break;
	
      if ((current_mu<0.01)&&(Xval<1e-2*std::sqrt(current_mu))&&(Zval>std::sqrt(current_mu)))
	  continue;
	
      // check whether the eigenvalue is going to zero
      // if there was a reduction from old_mu to mu, did this affect the primal more than the dual? 
	
      if (tapia_factor<.8){
	if(
	   (primal_tapia>0.8)||
	   //((1.1*(primal_tapia-tapia_factor)>(dual_tapia-tapia_factor)))
	   (primal_tapia>.9*dual_tapia)
	   )
	  //does not go to zero much faster than the dual value
	  break;
      }
      else {
	//no sufficient reduction in mu, judge by increase off the ratio of dual to primal
	if (Z0val/X0val>std::sqrt(current_mu)*Zval/Xval)
	  break;
      }
    
    }
    
    rank++;


    if (cb_out(2)){
      get_out()<<" rank="<<rank<<std::endl;
    }
    
    *pscx_activity_rank=rank;

    ///   ------------  even older code ---------------  
    //determine active rank
    Matrix Zval(sdim,1);  chk_set_init(Zval,1);
    Matrix oldXval(sdim,1); chk_set_init(oldXval,1);
    Matrix oldZval(sdim,1);  chk_set_init(oldZval,1);
    Matrix primal_tapia(sdim,1); chk_set_init(primal_tapia,1);
    Matrix dual_tapia(sdim,1); chk_set_init(dual_tapia,1);
    Matrix tmpvec;
    for (Integer i=0;i<sdim;i++){
      tmpvec=pscx_vecs.col(i);
      Zval(i)=ip(tmpvec,Zp[unsigned(psci)]*tmpvec);
      oldXval(i)=ip(tmpvec,old_Xp[unsigned(psci)]*tmpvec);
      oldZval(i)=ip(tmpvec,old_Zp[unsigned(psci)]*tmpvec);
      primal_tapia(i)=pscx_eigs(i)/oldXval(i);
      dual_tapia(i)=Zval(i)/oldZval(i);
    }

    Real tapia_factor=(old_mu>0.)?current_mu/old_mu:1.;
    
    if (cb_out(2)){
      get_out()<<" last_alpha="<<last_alpha;
      get_out()<<" mu="<<current_mu<<" old_mu="<<old_mu<<" tapia_factor="<<tapia_factor<<std::endl;
      for(Integer i=0;i<pscx_eigs.dim();i++){
	get_out()<<" "<<i<<":"<<pscx_eigs(i)<<","<<Zval(i);
	get_out()<<" "<<Real(Zval(0)/pscx_eigs(0)>min(1e-3,current_mu)*Zval(i)/pscx_eigs(i));
	get_out()<<" : "<<oldXval(i)<<","<<oldZval(i);
	get_out()<<"   ("<<primal_tapia(i)<<","<<dual_tapia(i)<<")";
	get_out()<<" "<<Real(primal_tapia(i)>dual_tapia(i))<<" ";
	if (tapia_factor<0.99){
	  get_out()<<"   ("<<(primal_tapia(i)-tapia_factor)/(1.-tapia_factor);
	  get_out()<<","<<(dual_tapia(i)-tapia_factor)/(1.-tapia_factor)<<") ";
	  get_out()<<" "<<Real((primal_tapia(i)>0.8)&&(1.1*((primal_tapia(i)-tapia_factor)/(1.-tapia_factor))>((dual_tapia(i)-tapia_factor)/(1.-tapia_factor))));
	}
	get_out()<<std::endl;						      
      }
    }
    
    Integer maxi=0;
    for (Integer i=0;i<pscx_eigs.dim();i++){
      // is the value too small in absolut terms?
      if (pscx_eigs(i)<1e-10*b)
	break;
      // is the value roughly the same as the dual or smaller?
      if (pscx_eigs(i)<std::sqrt(current_mu))
	break;
      // if there was a reduction from old_mu to mu, did this affect the primal more than the dual? 
      if (tapia_factor<.8){
	  if(((primal_tapia(i)>0.8)&&((1.1*((primal_tapia(i)-tapia_factor)/(1.-tapia_factor))>((dual_tapia(i)-tapia_factor)/(1.-tapia_factor))))))
	    maxi=i;
      }
      else {
	//no sufficient reduction in mu, judge by increase off the ratio of dual to primal
	if (Zval(0)/pscx_eigs(0)>std::sqrt(current_mu)*Zval(i)/pscx_eigs(i))
	  maxi=i;
	else 
	  break;
      }
    }

    *pscx_activity_rank=maxi+1;

    if (cb_out(2)){
      get_out()<<" actrank="<<maxi+1<<std::endl;
    }

  }
  */

  return 0;

}

// ****************************************************************************
//                           get_boxx
// ****************************************************************************


  /// after the bundle subproblem is solved, this retrieves the local linear solution vector; if nncx_activity is set, the values between zero and one indicate the guess on the coefficients activity level 
   int UQPConeModelBlock::get_boxx(Matrix& boxx,
			       Matrix* boxx_activity,
			       bool /* cautious */) 
  {
    boxx.init(box_dim,1,x.get_store()+box_start);
    
    if (boxx_activity){
      Indexmatrix inactive;
      
      //find_inactive_indices(inactive,b,false);
      boxx_activity->init(box_dim,1,1.);
      bool cautious=false;
      
      Real tapia_factor=(old_mu>0.)?current_mu/old_mu:1.;

      Real s=0.;
      if (box_scaling){
	s=x(box_start+lb.rowdim());
      }
      
      for (Integer i=0;i<lb.rowdim();i++){
	Real lslack=x(box_start+i);
	Real uslack=-lslack;
	if (box_scaling){
	  lslack-=lb(i)*s;
	  uslack+=ub(i)*s;
	}
	else {
	  s=0.;
	  lslack-=lb(i);
	  uslack+=ub(i);
	}
	if ((lslack>b*(cautious?1e-6:1e-6))&&(uslack>b*(cautious?1e-6:1e-6))){
	  //too big
	  continue;
	}
	//check for relation to the dual
	if ((z(box_start+i)<(cautious?1e6:1e3)*lslack)&&(uz(i)<(cautious?1e6:1e3)*uslack)){
	  continue;
	}
	if (min(lslack,uslack)<b*(cautious?1e-10:1e-10)){
	  //this is certainly not active
	  (*boxx_activity)(i)=0;
	  continue;
	}
	if ((current_mu<0.01)
	    &&(
	       ((lslack<0.01*std::sqrt(current_mu))&&(z(box_start+i)>std::sqrt(current_mu)))||
	       ((uslack<0.01*std::sqrt(current_mu))&&(uz(i)>std::sqrt(current_mu)))
	       )
	    ){
	  (*boxx_activity)(i)=0;
	  continue;
	}
	//check whether reduction is faster than for the dual
	if (tapia_factor<.8){
	  Real primal_tapia=lslack/(old_x(box_start+i)-(box_scaling?s*lb(i):lb(i)));
	  Real dual_tapia=z(box_start+i)/old_z(box_start+i);
	  if(
	     (primal_tapia<0.8)&&
	     //((1.1*(primal_tapia-tapia_factor)>(dual_tapia-tapia_factor)))
	     (primal_tapia<.1*dual_tapia)
	     ){
	    //goes to zero faster than the dual value
	    (*boxx_activity)(i)=0;
	    continue;
	  }
	  primal_tapia=uslack/((box_scaling?s*ub(i):ub(i))-old_x(box_start+i));
	  dual_tapia=uz(i)/old_uz(i);
	  if(
	     (primal_tapia<0.8)&&
	     //((1.1*(primal_tapia-tapia_factor)>(dual_tapia-tapia_factor)))
	     (primal_tapia<.1*dual_tapia)
	     ){
	    //goes to zero faster than the dual value
	    (*boxx_activity)(i)=0;
	    continue;
	  }
	}
      }
      
      if (box_scaling) {
	Real boxs=x(lb.rowdim()+box_start);
	Real box_lz=z(lb.rowdim()+box_start);
	bool clearly_active=(
			     ((boxs>b*(cautious?1e-6:1e-6))&&
			      ((!box_scaleub)||(boxs<b*(1-1e-6))))
			     ||
			     ((box_lz<(cautious?1e6:1e3)*boxs)&&
			      ((!box_scaleub)||(box_uz<(cautious?1e6:1e3)*boxs)))
			     );
	
	if (!clearly_active){
	  if  ((boxs<b*(cautious?1e-10:1e-10))
	       ||((box_scaleub)&&(b-boxs<b*(cautious?1e-10:1e-10)))
	       ||((current_mu<0.01)&&(boxs<0.01*std::sqrt(current_mu))&&(box_lz>std::sqrt(current_mu)))
	       ||((box_scaleub)&&(current_mu<0.01)&&(b-boxs<0.01*std::sqrt(current_mu))&&(box_uz>std::sqrt(current_mu)))
	       ){
	    //clearly inactive
	    (*boxx_activity)(lb.rowdim())=0;
	  }
	  else { 
	    //check whether reduction is faster than for the dual
	    Real primal_tapia=boxs/old_x(box_start+lb.rowdim());
	    Real dual_ltapia=box_lz/old_z(box_start+lb.rowdim());
	    if (tapia_factor<.8){
	      if((
		  (primal_tapia<0.8)&&
		  (primal_tapia<.1*dual_ltapia)
		  )||
		 (
		  (box_scaleub)&&
		  ((b-boxs)/(b-old_x(box_start+lb.rowdim()))<0.8)&&
		  ((b-boxs)/(b-old_x(box_start+lb.rowdim()))<.1*box_uz/old_box_uz)
		  )
		 )
		//goes to zero faster than the dual value
		(*boxx_activity)(lb.rowdim())=0;
	    }
	    
	  }
	}
      }
      
      if (cb_out(1)){
	get_out()<<" UQPBox: last_alpha="<<last_alpha;
	get_out()<<" mu="<<current_mu<<" old_mu="<<old_mu;
	get_out()<<" tapia_factor="<<tapia_factor<<" cautious="<<cautious<<" activity="<<transpose(*boxx_activity);
      }
       
    }
    
    return 0;
  }
    




// ****************************************************************************
//                           add_modelx_aggregate
// ****************************************************************************

int UQPConeModelBlock::add_modelx_aggregate(CH_Matrix_Classes::Real& offset,
			   CH_Matrix_Classes::Matrix& gradient)
{
  Matrix val(1,1,0.);
  if (x.rowdim()!=Integer(get_bundle().size())){
    if (cb_out())
      get_out()<<"**** ERROR in UQPConeModelBlock::add_modelx_aggregate(): x not yet computed"<<std::endl;
    return 1;
  }

  genmult(get_bundle(),x,gradient,1.,1.,0,0,&val);
  offset+=val(0);
  return 0;
}

// ****************************************************************************
//                           trace_dual
// ****************************************************************************

Real UQPConeModelBlock::tracedual(Real* out_prec) const
{
  if (out_prec)
    *out_prec=current_mu;
  return y;
}


// ****************************************************************************
//                          compute_local_directions
// ****************************************************************************

// the psd-line-search below should maybe better be done by eigenvalue computation

int UQPConeModelBlock::compute_local_directions(const Matrix& qp_dx,
					    const Matrix& qp_dy,
					    const Matrix& rhs_resid,
					    Matrix& dz,
					    Matrix& duz,
					    Real& box_duz,
					    Real& ds)
{
  dx.init(x.rowdim(),1,qp_dx.get_store()+qp_xstart);
  if (A.rowdim()>0){
    dy=qp_dy(qp_ystart);
    if (less_or_equal){
      if (y>scost) 
	ds=scorr-s*(1.+dy/(y-scost));
    }
  }

  assert((dz.rowdim()==x.rowdim())&&(dz.coldim()==1));

  //the cones are all slacks for the residual, this fixes their dz
  Integer cone_dim=x.rowdim()-box_dim;
  Integer xind=qp_xstart;
  mat_xey(cone_dim,dz.get_store(),rhs_resid.get_store()+xind);
  mat_xmey(cone_dim,dz.get_store(),z.get_store());
  xind+=cone_dim;

  //the box uses the corrector terms
  if (box_dim>0){
    if (box_scaling){

      Real xval=x(box_start+lb.rowdim());
      Real dxval=qp_dx(xind+lb.rowdim());

      // tmpvec=(xval*ub-x)
      tmpvec.init(ub,xval);
      mat_xmey(ub.rowdim(),tmpvec.get_store(),x.get_store()+box_start);
      // duz= uzcorr + ((dx-dxval*ub)/(xval*ub-x)-1.).*z
      duz.init(ub.rowdim(),1,dx.get_store()+box_start);
      duz.xpeya(ub,-dxval);
      mat_xinvhadey(ub.rowdim(),duz.get_store(),tmpvec.get_store());
      duz-=1;
      duz%=uz;
      duz+=uzcorr;
		    
      

      // tmpvec=(x-xval*lb)
      tmpvec.init(lb,-xval);
      mat_xpey(lb.rowdim(),tmpvec.get_store(),x.get_store()+box_start);
      // dz= zcorr + ((dxval*lb-dx.)/(x-xval*lb)-1.).*z
      mat_xeya(lb.rowdim(),dz.get_store()+box_start,lb.get_store(),dxval);
      mat_xmey(lb.rowdim(),dz.get_store()+box_start,dx.get_store()+box_start);
      mat_xinvhadey(ub.rowdim(),dz.get_store()+box_start,tmpvec.get_store());
      mat_xpea(lb.rowdim(),dz.get_store()+box_start,-1.);
      mat_xhadey(ub.rowdim(),dz.get_store()+box_start,z.get_store()+box_start);
      mat_xpey(ub.rowdim(),dz.get_store()+box_start,zcorr.get_store()+box_start);
     		
      dz(box_start+lb.rowdim())=zcorr(box_start+lb.rowdim())-(1.+dxval/xval)*z(box_start+lb.rowdim());

			      
      if (box_scaleub){
	box_duz= box_uzcorr-(1.+dxval/(b-xval))*box_uz;

      }

      // //TEST begin
      // Matrix tx(lb.rowdim(),1,x.get_store()+box_start);
      // Matrix tdx(lb.rowdim(),1,dx.get_store()+box_start);

      // Matrix tlz(lb.rowdim(),1,z.get_store()+box_start);
      // Matrix tlcorr(lb.rowdim(),1,zcorr.get_store());
      // Matrix lbs=tx-lb*dxval;
      // Matrix lbsi=lbs; lbsi.inv();
      // Matrix tdlz(lb.rowdim(),1,dz.get_store()+box_start);

      // Real tblz=z(box_start+lb.rowdim());
      // Real tbdlz=dz(box_start+lb.rowdim());
      // Real tbcorr=zcorr(box_start+lb.rowdim());

      // std::cout<<" tlcorr="<<transpose(tlcorr);

      // std::cout<<" res-(z+dz-uz-duz)="<<transpose(rhs_resid.rows(Range(xind,xind+lb.dim()-1))-tlz-tdlz+uz+duz);
		       
      // std::cout<<"xlcompl="<<transpose((tx+tdx-lb*(xval+dxval))%tlz+(tx-xval*lb)%tdlz-(tx-xval*lb)%tlcorr);
      // std::cout<<"xucompl="<<transpose((ub*(xval+dxval)-tx-tdx)%uz+(ub*xval-tx)%duz-(ub*xval-tx)%uzcorr);
      // std::cout<<" blcompl="<<((xval+dxval)*tblz+xval*tbdlz-xval*tbcorr)<<std::endl;

      // Real bres=rhs_resid(xind+lb.dim())-tblz-tbdlz+ip(lb,tlz+tdlz)-ip(ub,uz+duz);
      // if (box_scaleub){
      // 	std::cout<<" bucompl="<<((b-xval-dxval)*box_uz+(b-xval)*box_duz-(b-xval)*box_uzcorr)<<std::endl;
      // 	bres+=box_uz+box_duz;
      // }
      // std::cout<<" boxres="<<bres<<std::endl;
      // //Test end
      

    }
    else {
      // tmpvec=(ub-x)
      tmpvec.init(ub);
      mat_xmey(ub.rowdim(),tmpvec.get_store(),x.get_store()+box_start);
      // duz= uzcorr + (dx/(ub-x)-1.).*z
      duz.init(ub.rowdim(),1,dx.get_store()+box_start);
      mat_xinvhadey(ub.rowdim(),duz.get_store(),tmpvec.get_store());
      duz-=1;
      duz%=uz;
      duz+=uzcorr;
      
      // tmpvec=(x-lb)
      mat_xey(lb.rowdim(),tmpvec.get_store(),x.get_store()+box_start);
      tmpvec-=lb;
      // dz= zcorr + (-dx./(x-lb)-1.).*z
      mat_xeya(lb.rowdim(),dz.get_store()+box_start,dx.get_store()+box_start,-1.);
      mat_xinvhadey(ub.rowdim(),dz.get_store()+box_start,tmpvec.get_store());
      mat_xpea(lb.rowdim(),dz.get_store()+box_start,-1.);
      mat_xhadey(ub.rowdim(),dz.get_store()+box_start,z.get_store()+box_start);
      mat_xpey(ub.rowdim(),dz.get_store()+box_start,zcorr.get_store()+box_start);


      // //TEST begin
      // Matrix tx(lb.rowdim(),1,x.get_store()+box_start);
      // Matrix tdx(lb.rowdim(),1,dx.get_store()+box_start);

      // Matrix tlz(lb.rowdim(),1,z.get_store()+box_start);
      // Matrix tdlz(lb.rowdim(),1,dz.get_store()+box_start);
      // Matrix tlcorr(lb.rowdim(),1,zcorr.get_store());

      // std::cout<<"xlcompl="<<transpose((tx+tdx-lb)%tlz+(tx-lb)%tdlz-(tx-lb)%tlcorr);
      // std::cout<<"xucompl="<<transpose((ub-tx-tdx)%uz+(ub-tx)%duz-(ub-tx)%uzcorr);

      // //assert(norm2((tx+tdx-lb)%tlz+(tx-lb)%tdlz-(tx-lb)%tlcorr)<=1e-10*max(norm2(tx),norm2(tlz)));
      // //assert(norm2((ub-tx-tdx)%uz+(ub-tx)%duz-(ub-tx)%uzcorr)<=1e-10*max(norm2(tx),norm2(uz)));
      // //Test end

    }
  }
     
  return 0;
}
  
// ****************************************************************************
//                          inner_line_search
// ****************************************************************************

// the psd-line-search below should maybe better be done by eigenvalue computation

int UQPConeModelBlock::inner_line_search(Real& alpha,
				     const Matrix& qp_dx,
				     const Matrix& qp_dy,
				     const Matrix& dz,
				     const Matrix& duz,
				     Real box_duz,
				     Real ds)
{
  
  Real al=alpha;
  Real bdfac=.9;  //boundary factor
  if (less_or_equal){
    if (ds<0.)
      al=min(al,-bdfac*s/ds);
    if ((y>scost)&&(qp_dy(qp_ystart)<0.)) {
      Real d=-(y-scost)/qp_dy(qp_ystart);
      if (d<=al){
	al=bdfac*d;
      }
    }
  }

  //---linear part
  Integer xind=qp_xstart;
  {for(Integer i=0;i<nnc_dim;i++){
    if (qp_dx(xind+i)<0){
      Real d=-x(i)/qp_dx(xind+i);
      if (d<=al){
	al=bdfac*d;
      }
    } 
    if (dz(i)<0){
      Real d=-z(i)/dz(i);
      if (d<=al){
	al=bdfac*d;
      }
    } 
  }}
  xind+=nnc_dim;

  //---second order part
  for(Integer i=0;i<soc_dim.dim();i++){
    //primal variables
    const Real *xp=x.get_store()+xind-qp_xstart;
    const Real *dxp=qp_dx.get_store()+xind;
    Real a=(*dxp)*(*dxp);
    Real b=(*dxp)*(*xp);
    Real c=(*xp)*(*xp);
    xp++;
    dxp++;
    for(Integer j=soc_dim(i)-1;--j>=0;xp++,dxp++){
      a-=(*dxp)*(*dxp);
      b-=(*dxp)*(*xp);
      c-=(*xp)*(*xp);
    }
    if((fabs(a)<=eps_Real)&&(b<-eps_Real)){
      Real d=-c/b/2.;
      if (d<=al){
	al=bdfac*d;
      }
    } 
    else {
      Real d=b*b-a*c; 
      if ((a<-eps_Real)||((b<0.)&&(d>=0.))){
	Real f=(-b-sqrt(d))/a;
	if (f<=al){
	  al=bdfac*f;
	}
      } 
    }
    //dual variables
    xp=z.get_store()+xind-qp_xstart;
    dxp=dz.get_store()+xind-qp_xstart;
    a=(*dxp)*(*dxp);
    b=(*dxp)*(*xp);
    c=(*xp)*(*xp);
    xp++;
    dxp++;
    {for(Integer j=soc_dim(i)-1;--j>=0;xp++,dxp++){
      a-=(*dxp)*(*dxp);
      b-=(*dxp)*(*xp);
      c-=(*xp)*(*xp);
    }}
    if((fabs(a)<=eps_Real)&&(b<-eps_Real)){
      Real d=-c/b/2.;
      if (d<=al){
	al=bdfac*d;
      }
    } 
    else {
      Real d=b*b-a*c; 
      if ((a<-eps_Real)||((b<0.)&&(d>=0.))){
	Real f=(-b-sqrt(d))/a;
	if (f<=al){
	  al=bdfac*f;
	}
      } 
    }
    xind+=soc_dim(i);
  }

  //---semidefinite part
  {for(unsigned int i=0;i<Xp.size();i++){
    Integer Xrdim=Xp[i].rowdim();
    Integer xvdim=(Xrdim*(Xrdim+1))/2;
    tmpsvec.init(xvdim,1,qp_dx.get_store()+xind);
    sveci(tmpsvec,dX);
    Real dXtol=norm2(tmpsvec)*1e-10;
    tmpsymmat.init(dX,al); tmpsymmat+=Xp[i];
    while((tmpsymmat.Chol_factor(dXtol*al))&&(al>eps_Real)) {
      al*=.8;
      tmpsymmat.init(dX,al); tmpsymmat+=Xp[i];
    }
    tmpsvec.init(xvdim,1,dz.get_store()+xind-qp_xstart);
    sveci(tmpsvec,dX);
    Real dZtol=norm2(tmpsvec)*1e-10;
    tmpsymmat.init(dX,al); tmpsymmat+=Zp[i];
    while((tmpsymmat.Chol_factor(dZtol*al))&&(al>eps_Real)){
      al*=.8;
      tmpsymmat.init(dX,al); tmpsymmat+=Zp[i];
    }
    xind+=xvdim;  
  }}


  //box part
  if (box_dim>0){
    const Real* vp=x.get_store()+box_start;
    const Real* lp=lb.get_store();
    const Real* up=ub.get_store();
    const Real* dvp=qp_dx.get_store()+xind;
    const Real* const vendx=vp+lb.rowdim();
    if (box_scaling){
      Real box_s=x(box_start+lb.rowdim());
      Real box_ds=qp_dx(xind+lb.rowdim());
      while(vp!=vendx){
	Real d=*dvp++;
	Real l=*lp++;
	Real u=*up++;
	Real valx=*vp++;
	if (d<box_ds*l){
	  Real ax=(valx-box_s*l)/(box_ds*l-d);
	  if (ax<=al)
	    al=bdfac*ax;
	}
	if (d>box_ds*u){
	  Real ax=(box_s*u-valx)/(d-box_ds*u);
	  if (ax<=al)
	    al=bdfac*ax;
	}
      }
    }
    else {
      while(vp!=vendx){
	Real d=*dvp++;
	Real l=*lp++;
	Real u=*up++;
	Real valx=*vp++;
	if (d<0.){
	  d=-(valx-l)/d;
	  if (d<=al)
	    al=bdfac*d;
	}
	else if (d>0.){
	  d=(u-valx)/d;
	  if (d<=al)
	    al=bdfac*d;
	}
      }
    }
    vp=z.get_store()+box_start;
    dvp=dz.get_store()+box_start;
    const Real* const vendlz=vp+lb.rowdim();
    while(vp!=vendlz){
      Real d=*dvp++;
      if (d<0.){
	d=-(*vp++)/d;
	if (d<=al)
	  al=bdfac*d;
      }
      else {
	vp++;
      }
    }
    vp=uz.get_store();
    dvp=duz.get_store();
    const Real* const venduz=vp+lb.rowdim();
    while(vp!=venduz){
      Real d=*dvp++;
      if (d<0.){
	d=-(*vp++)/d;
	if (d<=al)
	  al=bdfac*d;
      }
      else {
	vp++;
      }
    }
    if (box_scaling){
      Real boxs=x(box_start+lb.rowdim());
      Real box_ds=qp_dx(xind+lb.rowdim());
      if (box_ds<0.){
	Real d=-boxs/box_ds;
	if (d<=al)
	  al=bdfac*d;
      }
      Real box_dlz=dz(box_start+lb.rowdim());
      if (box_dlz<0.){
	Real d=-z(box_start+lb.rowdim())/box_dlz;
	if (d<=al)
	  al=bdfac*d;
      }
      if (box_scaleub){
	if (box_ds>0){
	  Real d=(b-boxs)/box_ds;
	  if (d<=al)
	    al=bdfac*d;
	}
	if (box_duz<0.){
	  Real d=-box_uz/box_duz;
	  if (d<=al)
	    al=bdfac*d;
	}
      }
    }
    xind += box_dim;
  }

  
  if (al<alpha)
    alpha=al;

  return 0;
}




// ****************************************************************************
//                           starting_x
// ****************************************************************************

//generate a strictly feasible primal starting point
//store it in the qpx_range of x
//returns 0 on success, 1 on failure

int UQPConeModelBlock::starting_x(Matrix& qp_x)
{
  Real initval=b/mu_dim;
  s=0;
  if (A.rowdim()>0){
    if ((box_dim>0)&&(!box_scaling))
      initval=b/(mu_dim-2*lb.rowdim());
    x.init(A,initval);
    x.transpose();
    if (less_or_equal) 
      s=initval;
  }

  // PSC part
  for(unsigned int i=0;i<Xp.size();i++){
    Xp[i].init(Xp[i].rowdim(),0.);
    Xp[i].shift_diag(initval);
  }
  
  // box part
  if (box_dim>0){
    if (box_scaling){
      Real s=b*.5;
      if (!box_scaleub){
	s=(2*lb.rowdim()+1)*initval; //contribution to trace is now mu_dim*val
      }
      mat_xeya(lb.rowdim(),x.get_store()+box_start,lb.get_store(),.5*s);
      mat_xpeya(ub.rowdim(),x.get_store()+box_start,ub.get_store(),.5*s);
      x(box_start+lb.rowdim())=s;
    }
    else {
      mat_xeya(lb.rowdim(),x.get_store()+box_start,lb.get_store(),.5);
      mat_xpeya(ub.rowdim(),x.get_store()+box_start,ub.get_store(),.5);
    }
  }

  //copy into qp_x
  mat_xey(x.rowdim(),qp_x.get_store()+qp_xstart,x.get_store());
      
  return 0;
}


// ****************************************************************************
//                           starting_y
// ****************************************************************************

// generate a strictly feasible dual starting point
// store it in the qpy_range of y and in the qpx_range of z
// x is fixed already by a previous call to starting_x and Qx=Q*x
// returns 0 on success, 1 on failure

int UQPConeModelBlock::starting_y(Matrix& qp_y,const Matrix& qp_Qx, const Matrix& qp_c)
{
  //z = (qp_Qx - qp_c)[qp_xstart .. qp_xstart+x.rowdim()-1]
  z.init(x.rowdim(),1,qp_Qx.get_store()+qp_xstart);
  zcorr.init(x.rowdim(),1,0.);
  mat_xmey(z.dim(),z.get_store(),qp_c.get_store()+qp_xstart);

  y=1.+max(scost,0.);
  scorr=0.;
    
  //box upper bound set of dual variables
  if (box_dim>0){
    uzcorr.init(lb.rowdim(),1,0.);
    Real* lzp=z.get_store()+box_start;
    Real* uzp=uz.get_store();
    const Real* const endlzp=lzp+lb.rowdim();
    while(lzp<endlzp){
      Real val=(*lzp);
      if (val<0.){
	(*uzp++)=-val+1.;
	(*lzp++)=1.;
      }
      else {
	(*lzp++)=val+1.;
	(*uzp++)=1.;
      }
    }
    if (box_scaling){
      box_uzcorr=0.;
      Real val=(*lzp);
      val+=mat_ip(lb.rowdim(),lb.get_store(),z.get_store()+box_start);
      val-=ip(ub,uz);
      if (box_scaleub){
	box_uz=max(-val+1.,1.);
	(*lzp)=val+box_uz;
      }
      else {
	(*lzp)=val;
	y=max(y,1.-val);
      }
    }
  }

  //nnc part
  {for(Integer i=0;i<nnc_dim;i++){
    y=max(y,1.-z(i));
  }}
  Integer zind=nnc_dim;

  //soc part
  {for(Integer i=0;i<soc_dim.dim();i++){
    const Real *zp=z.get_store()+zind;
    Real a=(*zp++);
    Real b=0.;
    for(Integer j=soc_dim(i)-1;--j>=0;zp++){
      b+=(*zp)*(*zp);
    }
    y=max(y,sqrt(b)-a+1.);
    zind+=soc_dim(i);
  }}

  //psc part
  {for (unsigned int i=0;i<Zp.size();i++){
    Integer Zrdim=Xp[i].rowdim();
    Integer zvdim=(Zrdim*(Zrdim+1))/2;
    tmpsvec.init(zvdim,1,z.get_store()+zind);
    zind+=zvdim;
    sveci(tmpsvec,Zp[i]);
    y=max(y,max(sumrows(abs(Zp[i])))+1.);
  }}

  if (A.rowdim()>0){
    A.transpose();
    z.xpeya(A,y);
    A.transpose();
    {for (unsigned int i=0;i<Zp.size();i++){
	Zp[i].shift_diag(y);
      }}
    qp_y(qp_ystart)=y;
  }
  

  return 0;
}


// ****************************************************************************
//                       get_local_dualcost()
// ****************************************************************************

  CH_Matrix_Classes::Real UQPConeModelBlock::get_local_dualcost() const 
  {
    if (box_dim==0)
      return 0.;
    if (box_scaling){
      if (box_scaleub)
	return b*box_uz;
      return 0.;
    }
    return ip(ub,uz)-mat_ip(lb.rowdim(),lb.get_store(),z.get_store()+box_start);
  }
 

// ****************************************************************************
//                            get_Ab
// ****************************************************************************
  
//store the local coefficients of matrices A and b in the positions
//corresponding to qpy_range (rows) and qpx_range (columns) 
//returns 0 on success, 1 on failure

int UQPConeModelBlock::get_Ab(Matrix& qp_A,Matrix &qp_b) const
{
  if (A.rowdim()>0){
    mat_xey(A.dim(),qp_A.get_store()+qp_A.rowdim()*qp_xstart+qp_ystart,qp_A.rowdim(),A.get_store(),1);
    qp_b(qp_ystart)=b;
  }

  return 0;
}


// ****************************************************************************
//                           restart_x
// ****************************************************************************
  
//it is assumed that the problem was solved already once and is now
//resolved for a new linear cost term qp_c that resulted from the old
//one by adding qp_dc.
//on input qp_x holds the old optimal solution and on output
//the coorespoind qpx_range should be replaced by a reasonable 
//strictly feasible solution for x suitable for restarting
//(see also restart_yz)
//returns 0 on success, 1 on failure

int UQPConeModelBlock::restart_x(Matrix& qp_x,const Matrix& qp_c,const Matrix& qp_dc)
{
  Real normc=mat_ip(x.rowdim(),qp_c.get_store()+qp_xstart,qp_c.get_store()+qp_xstart);
  Real normdc=mat_ip(x.rowdim(),qp_dc.get_store()+qp_xstart,qp_dc.get_store()+qp_xstart);
  restart_factor=min(max(0.9,sqrt(1-min(1.,sqrt(normdc/normc)))),0.99999);
  if (A.rowdim()>0){
    A.transpose();
    xbpeya(x,A,(b-b*restart_factor)/mu_dim,restart_factor);
    A.transpose();
  }

  if (box_dim>0){
    mat_xpeya(ub.rowdim(),x.get_store()+box_start,ub.get_store(),.5*restart_factor);
    mat_xpeya(lb.rowdim(),x.get_store()+box_start,lb.get_store(),.5*restart_factor);
  }

  if (less_or_equal){
    s=s*restart_factor+(b-b*restart_factor)/mu_dim;
  }

  mat_xey(x.rowdim(),qp_x.get_store()+qp_xstart,x.get_store());

  Integer xind=nnc_dim;
  for (unsigned int i=0;i<Xp.size();i++){
    Integer Xrdim=Xp[i].rowdim();
    Integer xvdim=(Xrdim*(Xrdim+1))/2;
    tmpsvec.init(xvdim,1,x.get_store()+xind);
    xind+=xvdim;
    sveci(tmpsvec,Xp[i]);
  }

  
  return 0;
}


// ****************************************************************************
//                           restart_y
// ****************************************************************************
  
//this is called after restart_x (see there)
//on input qp_y holds the old optimal solution and on output
//the cooresponding qpy range should be replaced by a reasonable 
//strictly feasible solution for y suitable for restarting
//returns 0 on success, 1 on failure

int UQPConeModelBlock::restart_y(Matrix& qp_y,
			     const Matrix& qp_Qx,
			     const Matrix& qp_c,
			     const Matrix& /* qp_dc */)
{
  // tmpvec= A*y-z + qp_Qx - c
  if (A.rowdim()>0){
    A.transpose();
    xeyapzb(tmpvec,A,z,y,-1.);
    A.transpose();
  }
  else {
    tmpvec.init(z,-1.);
  }
  mat_xpey(tmpvec.dim(),tmpvec.get_store(),qp_Qx.get_store()+qp_xstart);
  mat_xmey(tmpvec.dim(),tmpvec.get_store(),qp_c.get_store()+qp_xstart);


  // make tmpvec positive definit  
  assert((!less_or_equal)||(y>scost));
  Real dy=max((y-scost)/10.,1e-3);
  {for(Integer i=0;i<nnc_dim;i++){
    dy=max(dy,fabs(tmpvec(i)));
  }}
  Integer zind=nnc_dim;
  {for(Integer i=0;i<soc_dim.dim();i++){
    const Real *zp=tmpvec.get_store()+zind;
    Real a=(*zp++);
    Real b=0.;
    for(Integer j=soc_dim(i)-1;--j>=0;zp++){
      b+=(*zp)*(*zp);
    }
    dy=max(dy,sqrt(b)+fabs(a));
    zind+=soc_dim(i);
  }}
  {for (unsigned int i=0;i<Zp.size();i++){
    Integer Zrdim=Zp[i].rowdim();
    Integer zvdim=(Zrdim*(Zrdim+1))/2;
    tmpsvec.init(zvdim,1,tmpvec.get_store()+zind);
    zind+=zvdim;
    sveci(tmpsvec,Xinv[i]);
    dy=max(dy,max(sumrows(abs(Xinv[i]))));
  }}
  if ((box_dim>0)&&(box_scaling)&&(!box_scaleub)) {
    Real val=tmpvec(box_start+lb.rowdim());
    val+=mat_ip(lb.rowdim(),lb.get_store(),z.get_store()+box_start)-ip(ub,uz);
    dy=max(dy,-val);
  }

  dy+=(1.-restart_factor)*.1;
  dy=max(dy,1e-5);
  z+=tmpvec;
  if (A.rowdim()>0.){
    A.transpose();
    z.xpeya(A,dy);
    A.transpose();
  }
  y+=dy;
  {for (unsigned int i=0;i<Zp.size();i++){
    Zp[i].shift_diag(dy);
    Zp[i]+=Xinv[i];
  }}
  
  qp_y(qp_ystart)=y;
  
  return 0;
}
    
// ****************************************************************************
//                        add_xinv_kron_z
// ****************************************************************************
  
//add the system term corresponding to (xinv kron z)
//(that arises from solving the linearized perturbed complementarity system
// x*z =0 or =mu*I for dx in the preferred search direction)
//to the diagonal block corresponding to qpx_range x qpx_range

int UQPConeModelBlock::add_xinv_kron_z(Symmatrix& barQ)
{
  Integer dind=qp_xstart;  //index into the diagonal of barQ
  {for(Integer i=0;i<nnc_dim;i++){
    barQ(dind,dind)+=z(i)/x(i);
    dind++;
  }}
  for (Integer i=0;i<soc_dim.dim();i++){
    const Real* xp=x.get_store()+dind-qp_xstart;
    const Real* zp=z.get_store()+dind-qp_xstart;
    Real xs=xp[0]*xp[0];
    for(Integer j=1;j<soc_dim(i);j++){
      xs-=xp[j]*xp[j];
    }
    Real xzip=mat_ip(soc_dim(i),xp,zp);
  
    barQ(dind,dind)+=(-xzip+2*xp[0]*zp[0])/xs;
    {for(Integer j=1;j<soc_dim(i);j++){
      barQ(dind,dind+j)+=(-zp[0]*xp[j]+xp[0]*zp[j])/xs;
      barQ(dind+j,dind+j)+=(xzip-2*zp[j]*xp[j])/xs;
      for(Integer k=j+1;k<soc_dim(i);k++){
	barQ(dind+j,dind+k)+=(-xp[j]*zp[k]-xp[k]*zp[j])/xs;
      }
    }}
    dind+=soc_dim(i);
  }
  {for (unsigned int i=0;i<Xp.size();i++){
    tmpsymmat=Xp[i];
    if (tmpsymmat.Chol_factor(1e-20)){
      if (cb_out()) get_out()<<"*** WARNING: UQPConeModelBlock::add_xinv_kron_z: factorizing Xp["<<i<<"] failed"<<std::endl;
      return 1;
    }
    tmpsymmat.Chol_inverse(Xinv[i]);
    skron(Xinv[i],Zp[i],tmpsymmat);
    for (Integer j=0;j<tmpsymmat.rowdim();j++){
      mat_xpey(tmpsymmat.rowdim()-j,
	       barQ.get_store()+barQ.rowdim()*dind-(dind*(dind-1))/2,
	       tmpsymmat.get_store()+tmpsymmat.rowdim()*j-(j*(j-1))/2);
      dind++;
    }
  }}
  //box part
  if (box_dim>0){
    if (!box_scaling){
      for (Integer i=0;i<lb.rowdim();i++,dind++){
	barQ(dind,dind)+=z(dind)/(x(dind)-lb(i))+uz(i)/(ub(i)-x(dind));
      }
    }
    else {
      //Symmatrix tQ(barQ); //TEST
      
      Real xval=x(box_start+lb.rowdim());
      Real lastval=0.;
      Integer lastind=dind+lb.rowdim();
      for (Integer i=0;i<lb.rowdim();i++,dind++){
	Real lval=z(dind)/(x(dind)-xval*lb(i));
	Real uval=uz(i)/(ub(i)*xval-x(dind));
	barQ(dind,dind)+=lval+uval;
	lval*=lb(i);
	lastval+=lb(i)*lval;
	uval*=ub(i);
	lastval+=ub(i)*uval;
	barQ(dind,lastind)-=lval+uval;
      }
      lastval+=z(dind)/x(dind);
      if (box_scaleub)
	lastval+=box_uz/(b-x(dind));
      barQ(dind,dind)+=lastval;

      // //TEST begin
      // Matrix tx(lb.rowdim(),1,x.get_store()+box_start);
      // Real zval=z(box_start+lb.rowdim());
      // Matrix tlz(lb.rowdim(),1,z.get_store()+box_start);
      // Matrix lbs=tx-xval*lb;
      // Matrix ubs=xval*ub-tx;
      // Matrix lbsi(lbs); lbsi.inv();
      // Matrix ubsi(ubs); ubsi.inv();
      // Matrix D=lbsi%tlz+ubsi%uz;
      // Real hatx=zval/xval+ip(lb,lbsi%tlz%lb)+ip(ub,ubsi%uz%ub);
      // Matrix hatb=lbsi%tlz%lb+ubsi%uz%ub;
      // if (box_scaleub)
      // 	hatx+=box_uz/(b-xval);
      // Symmatrix S(tQ.rowdim(),0.);
      // for(Integer i=0;i<lb.dim();i++){
      // 	S(dind-lb.dim()+i,dind-lb.dim()+i)=D(i);
      // 	S(dind-lb.dim()+i,dind)=-hatb(i);
      // }
      // S(dind,dind)=hatx;
      // std::cout<<" barQdiff="<<barQ-tQ-S;
      // std::cout<<" barQ="<<barQ;
      // std::cout<<" tQ="<<tQ;
      // std::cout<<" S="<<S;

      // //TEST end
    }
  }
	
  return 0;  
}

// ****************************************************************************
//                           add_local_sys
// ****************************************************************************
  
//on input: sysdy= A*barQ^{-1}*A^T    (barQ as returned in add_xinv_kron_z)
//          rhs= A*barQ^{-1}*(c-Q*x-A^T*y)-(b-A*x)
//if the block uses additional internal variables 
//(like an additional term + B*s with s>=0 in the primal feasibility constr)
//then the corresponding block terms have now to be added to sysdy and rhs,
//eg,
//   sysdy +=  B*(t^{-1} kron s)*B^T     (if t is the dual variable to s)
//   rhs   +=  B*s - B*(t^{-1} kron s)*B^T*y

int UQPConeModelBlock::add_local_sys(Symmatrix& sysdy,Matrix& /* rhs */)
{
  if (less_or_equal){
    if (s>=1e10*(y-scost)) {
      if (cb_out(1)){
	get_out()<<"*** WARNING: UQPConeModelBlock::add_local_sys adds large diagonal term: s="<<s<<" y="<<y<<std::endl;
      }
      /*
      y=0;
      for(Integer i=0;i<sysdy.rowdim();i++){
	sysdy(i,qp_ystart)=0.;
      }
      sysdy(qp_ystart,qp_ystart)=1.;
      rhs(qp_ystart)=0;
      return 0;
      */
    }
    sysdy(qp_ystart,qp_ystart)+=s/(y-scost);
  }
  return 0;
}

// ****************************************************************************
//                           suggest_mu
// ****************************************************************************
  
  //dx, dy is the predictor direction giving rise to the rhs_residual -(c-At(y+dy)-Q(x+dx)).
  //Compute the direction dz and local step and based on the predictor 
  //(x+dx,y+dy,z+dz) suggest a value for mu by specifying the
  //inner product of the dual cone variables ip_xz=ip(x,z)+ip(s,t),
  //the dimension of the conic variable space mu_dim= cone_x.rowdim+cone_s.dim
  //a value for the factor on mu to obtain the new target

int UQPConeModelBlock::suggest_mu(Real& ip_xz,
			      Integer& out_mu_dim,
			      Real& sigma,
			      const Matrix& qp_dx,
			      const Matrix& qp_dy,
			      const Matrix& rhs_residual)
{
  out_mu_dim=mu_dim;

  compute_local_directions(qp_dx,qp_dy,rhs_residual,dz,duz,box_duz,ds);
      
  Real al=1.;
  inner_line_search(al,qp_dx,qp_dy,dz,duz,box_duz,ds);

  ip_xz=mat_ip(x.rowdim()-box_dim,x.get_store(),z.get_store());
  Real ipdxdz=ip_xz;
  if ((less_or_equal)&&(y>scost)) {
    ip_xz+=s*(y-scost);
    ipdxdz+=al*al*ds*qp_dy(qp_ystart)+s*al*qp_dy(qp_ystart)+al*ds*(y-scost);
  }
  ipdxdz+=al*mat_ip(x.rowdim()-box_dim,x.get_store(),dz.get_store());
  ipdxdz+=al*mat_ip(x.rowdim()-box_dim,z.get_store(),qp_dx.get_store()+qp_xstart);
  ipdxdz+=al*al*mat_ip(x.rowdim()-box_dim,qp_dx.get_store()+qp_xstart,dz.get_store());

  //box part
  if (box_dim>0){
    Real cip=0.; //current inner product
    Real sip=0.; //step*al innner product
    Matrix tmpvec;

    if (box_scaling){
      
      Real xval=x(box_start+lb.rowdim());
      Real dxval=qp_dx(qp_xstart+box_start+lb.rowdim());

      //ip(x-s*lb,lz)
      tmpvec.init(lb,-xval);
      mat_xpey(lb.rowdim(),tmpvec.get_store(),x.get_store()+box_start);
      cip+=mat_ip(lb.rowdim(),tmpvec.get_store(),z.get_store()+box_start);
      mat_xpeya(lb.rowdim(),tmpvec.get_store(),qp_dx.get_store()+qp_xstart+box_start,al);
      tmpvec.xpeya(lb,-dxval*al);
      sip+=mat_ip(lb.rowdim(),tmpvec.get_store(),z.get_store()+box_start);
      sip+=al*mat_ip(lb.rowdim(),tmpvec.get_store(),dz.get_store()+box_start);
    
      //ip(s*ub-x,uz)
      tmpvec.init(ub,xval);
      mat_xmey(lb.rowdim(),tmpvec.get_store(),x.get_store()+box_start);
      cip+=ip(tmpvec,uz);
      mat_xpeya(lb.rowdim(),tmpvec.get_store(),qp_dx.get_store()+qp_xstart+box_start,-al);
      tmpvec.xpeya(ub,dxval*al);
      sip+=ip(tmpvec,uz)+al*ip(tmpvec,duz);

      //ip(s,lt)
      Real zval=z(box_start+lb.rowdim());
      Real dzval=dz(box_start+lb.rowdim());
      cip+=xval*zval;
      sip+=(xval+al*dxval)*(zval+al*dzval);

      if (box_scaleub){
	cip+=(b-xval)*box_uz;
	sip+=(b-xval-al*dxval)*(box_uz+al*box_duz);
      }

    }
    else {
      //ip(x-lb,lz)
      tmpvec.init(lb.rowdim(),1,x.get_store()+box_start);
      tmpvec-=lb;
      cip+=mat_ip(lb.rowdim(),tmpvec.get_store(),z.get_store()+box_start);
      mat_xpeya(lb.rowdim(),tmpvec.get_store(),qp_dx.get_store()+qp_xstart+box_start,al);
      sip+=mat_ip(lb.rowdim(),tmpvec.get_store(),z.get_store()+box_start);
      sip+=al*mat_ip(lb.rowdim(),tmpvec.get_store(),dz.get_store()+box_start);
      
      //ip(ub-x,uz)
      tmpvec.init(ub.rowdim(),1,x.get_store()+box_start,1,-1.);
      tmpvec+=ub;
      cip+=ip(tmpvec,uz);
      mat_xpeya(ub.rowdim(),tmpvec.get_store(),qp_dx.get_store()+qp_xstart+box_start,-al);
      sip+=ip(tmpvec,uz)+al*ip(tmpvec,duz);
      
    }
    ip_xz+=cip;
    ipdxdz+=sip;
  }

  
  sigma=min(1.,pow(ipdxdz/ip_xz,max(2.,3*sqr(al))));
  return 0;
}

// ****************************************************************************
//                           get_corr
// ****************************************************************************
  
//on input (w.r.t. corresponding positions)
//    xcorr = 0
//    rhs as on output of add_local_sys
//on output the corresponding positions of xcorr should hold the corrector
//term of the search direction, eg,  xcorr = mu*x^{-1} - x^{-1}*dx*dz,
//and if the block holds additional local variables as in add_local_sys then
//   rhs += B*(mu * t^{-1}- t^{-1}*dt*ds) 

int UQPConeModelBlock::get_corr(Matrix& xcorr,
			    Matrix& rhs,
			    Real mu)
{
  old_mu=current_mu;
  current_mu=mu;
  Integer xind=qp_xstart;  
  {for(Integer i=0;i<nnc_dim;i++){
      xcorr(i+xind)=zcorr(i)=(mu-dx(i)*dz(i))/x(i);
  }}
  xind+=nnc_dim;

  //------- compute corrector rhs for SOC

  for (Integer i=0;i<soc_dim.dim();i++){
    const Real* xp=x.get_store()+xind-qp_xstart;
    Real x0=xp[0];
    Real xbarsq=mat_ip(soc_dim(i)-1,xp+1,xp+1);
    Real gamma=sqrt(x0*x0-xbarsq);
    Real d=1./(gamma*(xbarsq/(gamma+x0)-x0));
    xcorr(xind)=mu*(-1/gamma-(gamma+x0)*d);
    mat_xeya(soc_dim(i)-1,xcorr.get_store()+xind+1,xp+1,mu*d);
    Matrix socx(soc_dim(i),1,xp);
    Matrix socdx(soc_dim(i),1,dx.get_store()+xind-qp_xstart);
    Matrix socdz(soc_dim(i),1,dz.get_store()+xind-qp_xstart);
    Matrix res=mult_Ginv(socx,gamma,mult_Arw(mult_G(socx,gamma,socdz),mult_Ginv(socx,gamma,socdx)));
    mat_xpeya(soc_dim(i),xcorr.get_store()+xind,res.get_store(),-1.);
    mat_xey(soc_dim(i),zcorr.get_store()+xind-qp_xstart,xcorr.get_store()+xind);
    xind+=soc_dim(i);
  }

  //------- compute corrector rhs for SDP

  {for (unsigned int i=0;i<Xp.size();i++){
    Integer Xrdim=Xp[i].rowdim();
    Integer xvdim=(Xrdim*(Xrdim+1))/2;
    tmpsvec.init(xvdim,1,dx.get_store()+xind-qp_xstart);  
    sveci(tmpsvec,tmpsymmat);
    genmult(Xinv[i],Matrix(tmpsymmat),tmpvec,-1.);
    tmpsvec.init(xvdim,1,dz.get_store()+xind-qp_xstart); 
    sveci(tmpsvec,tmpsymmat);
    genmult(tmpvec,tmpsymmat,tmpsvec);
    tmpsymmat.init(tmpsvec);
    tmpsymmat.xpeya(Xinv[i],mu);
    svec(tmpsymmat,tmpsvec);
    mat_xey(xvdim,xcorr.get_store()+xind,tmpsvec.get_store());
    mat_xey(xvdim,zcorr.get_store()+xind-qp_xstart,tmpsvec.get_store());
    xind+=xvdim;
  }}
  
  if (less_or_equal){
    //ds=-s(1+dy/(y-scost)); (mu-dy*ds)/(y-scost)=(mu+dy*s*(1+dy/(y-scost)))/(y-scost)
    if (y>scost){
      scorr=(mu-dy*ds)/(y-scost);
      rhs(qp_ystart)+= scorr;
    }
    else {
      scorr=0.;
      rhs(qp_ystart)=0.;
    }
  }

  //---------- compute corrector rhs for Box
  if (box_dim>0){
    if (box_scaling){
      tmpvec.init(lb.rowdim(),1,x.get_store()+box_start);
      Real xval=x(box_start+lb.rowdim());
      Real dxval=dx(box_start+lb.rowdim());
      tmpvec.xpeya(lb,-xval);
      tmpvec.inv();
      tmpsvec.init(lb.rowdim(),1,dx.get_store()+box_start,1,-1.);
      tmpsvec.xpeya(lb,dxval);
      mat_xhadey(lb.rowdim(),tmpsvec.get_store(),dz.get_store()+box_start);
      tmpsvec+=mu;
      tmpsvec%=tmpvec;
      mat_xey(lb.rowdim(),zcorr.get_store()+box_start,tmpsvec.get_store());
      mat_xey(lb.rowdim(),xcorr.get_store()+xind,tmpsvec.get_store());
      zcorr(box_start+lb.rowdim())=(mu-dxval*dz(box_start+lb.rowdim()))/xval;
      xcorr(xind+lb.rowdim())=zcorr(box_start+lb.rowdim());
      xcorr(xind+lb.rowdim())-=ip(lb,tmpsvec);
      
      tmpvec.init(ub,xval);
      mat_xmey(ub.rowdim(),tmpvec.get_store(),x.get_store()+box_start);
      tmpvec.inv();
      tmpsvec.init(ub.rowdim(),1,dx.get_store()+box_start);
      tmpsvec.xpeya(ub,-dxval);
      tmpsvec%=duz;
      tmpsvec+=mu;
      tmpsvec%=tmpvec;
      uzcorr=tmpsvec;
      mat_xmey(ub.rowdim(),xcorr.get_store()+xind,tmpsvec.get_store());
      xcorr(xind+lb.rowdim())+=ip(ub,tmpsvec);
      
      if (box_scaleub){
	box_uzcorr=(mu+dxval*box_duz)/(b-xval);
	xcorr(xind+lb.rowdim())-=box_uzcorr;
      }

      // //TEST begin
      // Matrix tx(lb.rowdim(),1,x.get_store()+box_start);
      // Matrix tdx(lb.rowdim(),1,dx.get_store()+box_start);

      // Matrix tlz(lb.rowdim(),1,z.get_store()+box_start);
      // Matrix tdlz(lb.rowdim(),1,dz.get_store()+box_start);
      // Matrix tlcorr(lb.rowdim(),1,zcorr.get_store()+box_start);
      // Matrix txcorr(lb.rowdim(),1,xcorr.get_store()+xind);
      // Real txbcorr=xcorr(xind+lb.rowdim());
      // Matrix lbs=tx-xval*lb;
      // Matrix ubs=xval*ub-tx;
      // Matrix lbsi(lbs); lbsi.inv();
      // Matrix ubsi(ubs); ubsi.inv();
      

      // Real tblz=z(box_start+lb.rowdim());
      // Real tbdlz=dz(box_start+lb.rowdim());
      // Real tbcorr=zcorr(box_start+lb.rowdim());


      // std::cout<<"lzcorr="<<transpose(tlcorr+((tdx-dxval*lb)%tdlz-mu)%lbsi);
      // std::cout<<"uzcorr="<<transpose(uzcorr+((dxval*ub-tdx)%duz-mu)%ubsi);
      // std::cout<<"blzcorr="<<tbcorr+(dxval*tbdlz-mu)/xval<<std::endl;
      // std::cout<<"xcorr="<<transpose(txcorr-tlcorr+uzcorr);
      // txbcorr+= +ip(lb,tlcorr) - ip(ub,uzcorr) - tbcorr;
      // if (box_scaleub){
      // 	std::cout<<" buzcorr="<<box_uzcorr-(dxval*box_duz+mu)/(b-xval)<<std::endl;
      // 	txbcorr+=box_uzcorr;
      // }
      // std::cout<<"txbcorr="<<txbcorr<<std::endl;
      // //Test end

    }
    else {
      tmpvec.init(lb.rowdim(),1,x.get_store()+box_start);
      tmpvec-=lb;
      tmpsvec.init(lb.rowdim(),1,dx.get_store()+box_start,1,-1.);
      mat_xhadey(lb.rowdim(),tmpsvec.get_store(),dz.get_store()+box_start);
      tmpsvec+=mu;
      mat_xinvhadey(lb.rowdim(),tmpsvec.get_store(),tmpvec.get_store());
      mat_xey(lb.rowdim(),zcorr.get_store()+box_start,tmpsvec.get_store());
      mat_xey(lb.rowdim(),xcorr.get_store()+xind,tmpsvec.get_store());

      tmpvec.init(ub.rowdim(),1,x.get_store()+box_start,1,-1.);
      tmpvec.xpeya(ub);
      tmpsvec.init(ub.rowdim(),1,dx.get_store()+box_start,1);
      tmpsvec%=duz;
      tmpsvec+=mu;
      mat_xinvhadey(ub.rowdim(),tmpsvec.get_store(),tmpvec.get_store());
      uzcorr=tmpsvec;
      mat_xmey(ub.rowdim(),xcorr.get_store()+xind,tmpsvec.get_store());   
    }
  }


  return 0;  
  
}

// ****************************************************************************
//                            line_search
// ****************************************************************************
  
//dx,dy,dz gives the final step direction, alphap and alphad are on input
//upper bounds on the primal and dual step size. 
//On output alpahp and alphad have to be no larger than on input and
//have to guarantee strict feasibility of the primal/dual step on
//the local variables
//(if the block has additional internal variables it has to compute the step
// direction for these now and to choose alphap and alphad so that
// strict feasibility is guaranteed for the internal variables as well)

int UQPConeModelBlock::line_search(Real& alpha,
			       const Matrix& qp_dx,
			       const Matrix& qp_dy,
			       const Matrix& rhs_residual)
{
  compute_local_directions(qp_dx,qp_dy,rhs_residual,dz,duz,box_duz,ds);
      
  inner_line_search(alpha,qp_dx,qp_dy,dz,duz,box_duz,ds);

  //in case of doubt or for testing uncomment this  (but corret it before, not yet updated
  /*
  for (Integer i=0;i<nnc_dim;i++){
    if ((x(i)+alpha*qp_dx(qp_xstart+i)<=0.)||(z(i)+alpha*qp_dz(qp_xstart+i)<=0.)){
      if (cb_out()) get_out()<<"*** ERROR: UQPConeModelBlock::line_search(): liner variables not >=0 for i="<<i<<std::endl;
    }
  }

  Integer xind=nnc_dim;
  for (Integer i=0;i<soc_dim.dim();i++){
    Matrix xp(soc_dim(i),1,x.get_store()+xind);
    mat_xpeya(soc_dim(i),xp.get_store(),qp_dx.get_store()+qp_xstart+xind,alpha);
    Real xs=xp[0]*xp[0];
    for(Integer j=1;j<soc_dim(i);j++){
      xs-=xp[j]*xp[j];
    }
    if (xs<0.){
      if (cb_out()) get_out()<<"*** ERROR: UQPConeModelBlock::line_search(): soc_x not soc for i="<<i<<std::endl;
    }
    Matrix zp(soc_dim(i),1,z.get_store()+xind);
    mat_xpeya(soc_dim(i),zp.get_store(),qp_dz.get_store()+qp_xstart+xind,alpha);
    Real zs=zp[0]*zp[0];
    for(Integer j=1;j<soc_dim(i);j++){
      zs-=zp[j]*zp[j];
    }
    if (zs<0.){
      if (cb_out()) get_out()<<"*** ERROR: UQPConeModelBlock::line_search(): soc_z not soc for i="<<i<<std::endl;
    }
    xind+=soc_dim(i);
  }
  for (Integer i=0;i<Xp.size();i++){
    Integer Xrdim=Xp[i].rowdim();
    Integer xvdim=(Xrdim*(Xrdim+1))/2;
    tmpsvec.init(xvdim,1,x.get_store()+xind);
    mat_xpeya(xvdim,tmpsvec.get_store(),qp_dx.get_store()+qp_xstart+xind,alpha);
    sveci(tmpsvec,tmpsymmat);
    if (tmpsymmat.Chol_factor()){
      if (cb_out()) get_out()<<"*** ERROR: UQPConeModelBlock::line_search(): X not psd for i="<<i<<std::endl;
    }
    tmpsvec.init(xvdim,1,z.get_store()+xind);
    mat_xpeya(xvdim,tmpsvec.get_store(),qp_dz.get_store()+qp_xstart+xind,alpha);
    sveci(tmpsvec,tmpsymmat);
    if (tmpsymmat.Chol_factor()){
      if (cb_out()) get_out()<<"*** ERROR: UQPConeModelBlock::line_search(): X not psd for i="<<i<<std::endl;
    }
    xind+=xvdim;
  }
  */

  return 0;
}

// ****************************************************************************
//                           set_point
// ****************************************************************************
  
  //x,y,z is the new point and has to be stored
  //alphap and alphad are the step size used in the step, they are passed
  //so that the block can take the same step for internal variables if
  //needed.

int UQPConeModelBlock::set_point(const Matrix& qp_x,
			     const Matrix& qp_y,
			     Real alpha )
{
  last_alpha=alpha;
  
  Integer xind=qp_xstart;
  swap(old_x,x);
  x.init(old_x.rowdim(),1,qp_x.get_store()+xind);


  old_z=z;
  z.xpeya(dz,alpha);
  zcorr.init(z.rowdim(),1,0.);

  if (A.rowdim()>0){
    old_y=y;
    y=qp_y(qp_ystart);
    if (less_or_equal){
      old_s=s;
      scorr=0;
      assert((A.coldim()==A.dim())&&(A.dim()==x.rowdim()));
      s=b-mat_ip(A.coldim(),A.get_store(),x.get_store());
    }
  }
  
  
  xind=nnc_dim+sum(soc_dim);

  for (unsigned int i=0;i<Xp.size();i++){
    swap(old_Xp[i],Xp[i]);
    swap(old_Zp[i],Zp[i]);
    
    Integer Xrdim=old_Xp[i].rowdim();
    Integer xvdim=(Xrdim*(Xrdim+1))/2;
    tmpsvec.init(xvdim,1,x.get_store()+xind);
    sveci(tmpsvec,Xp[i]);
    tmpsvec.init(xvdim,1,z.get_store()+xind);
    sveci(tmpsvec,Zp[i]);
    xind+=xvdim;
  }
  
  if (box_dim>0){
    old_uz=uz;
    uz.xpeya(duz,alpha);
    uzcorr.init(uz.rowdim(),1,0.);

    if (box_scaling && box_scaleub){
      old_box_uz=box_uz;
      box_uz+=alpha*box_duz;
      box_uzcorr=0.;
    }
  }

  //in case of doubt or for testing uncomment this
  /*
  for (Integer i=0;i<nnc_dim;i++){
    if ((x(i)<=0.)||(z(i)<=0.)){
      if (cb_out()) get_out()<<"*** ERROR: UQPConeModelBlock::set_point(): liner variables not >=0 for i="<<i<<std::endl;
    }
  }
  xind=nnc_dim;
  for (Integer i=0;i<soc_dim.dim();i++){
    const Real* xp=x.get_store()+xind;
    Real xs=xp[0]*xp[0];
    for(Integer j=1;j<soc_dim(i);j++){
      xs-=xp[j]*xp[j];
    }
    if (xs<0.){
      if (cb_out()) get_out()<<"*** ERROR: UQPConeModelBlock::set_point(): soc_x not soc for i="<<i<<std::endl;
    }
    const Real* zp=z.get_store()+xind;
    Real zs=zp[0]*zp[0];
    for(Integer j=1;j<soc_dim(i);j++){
      zs-=zp[j]*zp[j];
    }
    if (zs<0.){
      if (cb_out()) get_out()<<"*** ERROR: UQPConeModelBlock::set_point(): soc_z not soc for i="<<i<<std::endl;
    }
    xind+=soc_dim(i);
  }
  for (Integer i=0;i<Xp.size();i++){
    tmpsymmat=Xp[i];
    if (tmpsymmat.Chol_factor()){
      if (cb_out()) get_out()<<"*** ERROR: UQPConeModelBlock::set_point(): X not psd for i="<<i<<std::endl;
    }
    tmpsymmat=Zp[i];
    if (tmpsymmat.Chol_factor()){
      if (cb_out()) get_out()<<"*** ERROR: UQPConeModelBlock::set_point(): X not psd for i="<<i<<std::endl;
    }
  }
  */

  return 0;
}


    
// ****************************************************************************
//                         add_Bs
// ****************************************************************************

Matrix& UQPConeModelBlock::add_Bs(Matrix& qp_vec) const
{
  if (less_or_equal){
    qp_vec(qp_ystart)+=s;
  }
  return qp_vec;
}

// ****************************************************************************
//                         subtract_z
// ****************************************************************************

Matrix& UQPConeModelBlock::subtract_z(Matrix& qp_vec,bool with_step) const
{
  mat_xmey(x.dim(),qp_vec.get_store()+qp_xstart,z.get_store());
  if (box_dim>0){
    mat_xpey(x.dim(),qp_vec.get_store()+qp_xstart+box_start,uz.get_store());
    if (box_scaling){
      qp_vec(qp_xstart+box_start+lb.rowdim())+=mat_ip(lb.rowdim(),lb.get_store(),z.get_store()+box_start)-ip(ub,uz);
      if (box_scaleub)
	qp_vec(qp_xstart+box_start+lb.rowdim())+=box_uz;
    }
  }
  if (with_step){
    mat_xmey(x.dim(),qp_vec.get_store()+qp_xstart,dz.get_store());
    if (box_dim>0){
      mat_xpey(x.dim(),qp_vec.get_store()+qp_xstart+box_start,duz.get_store());
      if (box_scaling){
	qp_vec(qp_xstart+box_start+lb.rowdim())+=mat_ip(lb.rowdim(),lb.get_store(),dz.get_store()+box_start)-ip(ub,duz);
	if (box_scaleub)
	  qp_vec(qp_xstart+box_start+lb.rowdim())+=box_duz;
      }
    }
  }
  return qp_vec;
}


}

