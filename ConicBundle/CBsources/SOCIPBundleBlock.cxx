/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/SOCIPBundleBlock.cxx
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


#include "SOCIPBundleBlock.hxx"

using namespace CH_Matrix_Classes;

namespace ConicBundle {

  void SOCIPBundleBlock::clear(Integer dim)
  {
    SOCIPBlock::clear(dim);
    bundle_dim=dim;
    diff_model.init(0,0,0.);
    
    //oracle_data=0;
  }
  
  SOCIPBundleBlock::SOCIPBundleBlock(Integer dim,CBout* cb,int cbinc):
    CBout(cb,cbinc),InteriorPointBundleBlock(cb,cbinc),SOCIPBlock(dim,cb,cbinc)
  {clear(dim);}

  SOCIPBundleBlock::~SOCIPBundleBlock()
  {}
  

  InteriorPointBundleBlock* SOCIPBundleBlock::clone()
  {
    SOCIPBundleBlock* p=new SOCIPBundleBlock(0,this,0);
    p->copy_from(this);

    return p;
  }
  
  int SOCIPBundleBlock::copy_from(InteriorPointBundleBlock* inp)
  {
    SOCIPBundleBlock* p=dynamic_cast<SOCIPBundleBlock*>(inp);
    if (p==0)
      return 1;
   
    bundle_dim=p->bundle_dim;
    diff_model=p->diff_model;
    Bt=p->Bt;
    Boffset=p->Boffset;
    sqrBnorms=p->sqrBnorms;
    
    //oracle_data=p->oracle_data;
    
    vecdim=p->vecdim;
    x=p->x;
    z=p->z;
    dx=p->dx;
    dz=p->dz;
    gammaxsqr=p->gammaxsqr;
    gammazsqr=p->gammazsqr;
    omega=p->omega;
    f=p->f;
    scaled_point=p->scaled_point;
    compl_rhs=p->compl_rhs;
    last_rhs_mu=p->last_rhs_mu;
    mu=p->mu;
    old_mu=p->old_mu;
    last_alpha=p->last_alpha;
    oldx=p->oldx;
    oldz=p->oldz;
    tmpvec=p->tmpvec;
    tmpmat=p->tmpmat;
    tmp_xdzpdxz=p->tmp_xdzpdxz;
    
    return 0;
  }

  Real SOCIPBundleBlock::evaluate_trace_x()
  {return x(0);}

  Real SOCIPBundleBlock::evaluate_trace_z()
  {return z(0);}

  Real SOCIPBundleBlock::evaluate_trace_dx()
  {return dx(0);}


    /// C=beta*C+alpha*B*B where B and B may be transposed; carry out the model part of this beginning at startindex_model and beta for the part, that is added to (the calling routine has to make sure beta is not executed repeatedly if the same part is affected by other models as well)
Matrix& SOCIPBundleBlock::B_times(const Matrix& A,
				Matrix& C,
				Real alpha,
				Real beta,
				int Btrans,
				int Atrans,
				Integer startindex_model,
				MinorantBundle& globalbundle,
				Integer startindex_bundle)
{
  if ((startindex_bundle==0)&&(startindex_model==0)&&(beta!=1.)){
    if (beta==0.)
      C.init(C.rowdim(),C.coldim(),0.);
    else 
      C*=beta;
  }
  for (Integer i=0;i<vecdim;i++){
    globalbundle[unsigned(startindex_bundle+i)].left_genmult(A,C,alpha,1.,(Btrans==0),Atrans,startindex_model++);
  }
  return C;
}
						   
  /// C=beta*C+alpha*A*B where B and B may be transposed; carry out the model part of this beginning at startindex_model 
Matrix& SOCIPBundleBlock::times_B(const Matrix& A,
				Matrix& C,
				Real alpha,
				Real beta,
				int Atrans,
				int Btrans,
				Integer startindex_model,
				MinorantBundle& globalbundle,
				Integer startindex_bundle)
{
  if ((startindex_bundle==0)&&(startindex_model==0)&&(beta!=1.)){
    if (beta==0.)
      C.init(C.rowdim(),C.coldim(),0.);
    else 
      C*=beta;
  }
  for (Integer i=0;i<vecdim;i++){
    globalbundle[unsigned(startindex_bundle+i)].right_genmult(A,C,alpha,1.,Atrans,(Btrans==0),startindex_model++);
  }
  return C;
}

  ///add the main diagonal block tranpose(projection)*diagvec*projection to bigS starting at startindex
Symmatrix& SOCIPBundleBlock::add_BDBt(const Matrix& diagvec,
				    Symmatrix& bigS,
				    bool minus,
				    Integer startindex,
				    Matrix& Bt,
				    Integer startindex_model,
				    MinorantBundle& globalbundle,
				    Integer startindex_bundle)
{
  if (minus) {
    for (Integer i=0;i<vecdim;i++){
      MinorantPointer& mp=globalbundle[unsigned(startindex_bundle+i)];
      for(Integer j=i+startindex_model;j<Bt.coldim();j++){
	bigS(startindex+startindex_model+i,startindex+j)-=mp.ip(Bt,&diagvec,j*Bt.rowdim());
      }
    }
  }
  else {
    for (Integer i=0;i<vecdim;i++){
      MinorantPointer& mp=globalbundle[unsigned(startindex_bundle+i)];
      for(Integer j=i+startindex_model;j<Bt.coldim();j++){
	bigS(startindex+startindex_model+i,startindex+j)+=mp.ip(Bt,&diagvec,j*Bt.rowdim());
      }
    }
  }
  return bigS;
}

  /// get the current matrix for the coupling matrix Bt in the first row of blocks
Matrix& SOCIPBundleBlock::get_Bt(Matrix& globBt,
			       Integer startindex_model,
			       MinorantBundle& globalbundle,
			       Integer startindex_bundle)
{
  assert((&Bt!=&globBt)||(startindex_model==0));
  if (Bt.coldim()!=vecdim){
    Bt.newsize(globBt.rowdim(),vecdim); chk_set_init(Bt,1);
    Boffset.newsize(vecdim,1); chk_set_init(Boffset,1);
    for (Integer i=0;i<vecdim;i++){
      globalbundle[unsigned(startindex_bundle+i)].get_minorant(Boffset(i),Bt,i);
    }
  }
  if (&Bt != &globBt){
    assert(Bt.rowdim()==globBt.rowdim());
    assert(globBt.coldim()>=startindex_model+vecdim);
    mat_xey(Bt.dim(),globBt.get_store()+startindex_model*Bt.rowdim(),Bt.get_store());
  }
  return globBt;
}
  
  /// adds opB transposed times modelx (without constant affine term) to the arguments
int SOCIPBundleBlock::add_modelx_aggregate(Real& val,
				   Matrix& vec,
				   MinorantBundle& globalbundle,
				   Integer startindex_bundle)
{
  for (Integer i=0;i<vecdim;i++){
    globalbundle[unsigned(startindex_bundle+i)].get_minorant(val,vec,0,x(i),true);
  }
  return 0;
}

  /// set the model violation for the current system solution 
int SOCIPBundleBlock::get_sysviol_model(Matrix& sysviol_model,
				      Integer startindex_model,
				      const Matrix& dy,
				      const Real deltatrdual,
				      MinorantBundle& globalbundle,
				      Integer startindex_bundle)
{
  for(Integer i=0;i<vecdim;i++){
    sysviol_model(startindex_model+i)=globalbundle[unsigned(startindex_bundle+i)].evaluate(-1,dy,false)-diff_model(i)+z(i)+dz(i);
  }
  sysviol_model(startindex_model)-=deltatrdual;
  return 0;
}

  /// set z to the slack of the bundle and return a value>=0 that needs to be added to make it feasible
int SOCIPBundleBlock::set_bundle_z(const Matrix& y,
		    MinorantBundle& globalbundle,
		    Integer startindex_bundle,
		    Real& add_center_value)
{
  diff_model.newsize(vecdim,1);chk_set_init(diff_model,1); 
  for(Integer i=0;i<vecdim;i++){
    diff_model(i)=-globalbundle[unsigned(startindex_bundle+i)].evaluate(-1,y);
  }
  return SOCIPBlock::set_z(diff_model,0,add_center_value);
}

int SOCIPBundleBlock::add_trace(Matrix& vec,Real trace_dual,Integer startindex)
{
  vec(startindex)+=trace_dual;
  return 0;
}


    /// set the trace premultiplied by sqrt(inv(xiz)) in vec[startindex+0,...,startindex+dim_bundle()-1]
  int SOCIPBundleBlock::set_xizinvsqrt_trace(Matrix& vec,
					   Integer startindex)
  {
    assert(vecdim==bundle_dim);    
    if (f.dim()!=vecdim){
      compute_NTscaling();
    }

    vec(startindex)=f(0)/omega;
    for(Integer i=1;i<vecdim;i++){
      vec(startindex+i)=-f(i)/omega;
    }
    return 0;
  }
  


  ///add trace_dual*trace to diff_model for the right hand side (negative of current model violation)
int SOCIPBundleBlock::add_trace_to_diff_model(Real trace_dual)
{
  diff_model(0)+=trace_dual;
  return 0;
}

///return the squared Euclidean norm of the dual model violation  
Real SOCIPBundleBlock::dualviol_2normsqr()
{
  tmpvec.init(diff_model);
  tmpvec-=z;
  return ip(tmpvec,tmpvec);
}

  /// move to (x+alpha*dx, z+alpha*dz), update diff_model and possibly reduce the model size if some part is below zero_threshold
int SOCIPBundleBlock::do_bundle_step(Real alpha,
				     const Matrix& y,
				     MinorantBundle& globalbundle,
				     Integer startindex_bundle,
				     Real trace_dual,
				     Real /* trace_rhs */)
{
  SOCIPBlock::do_step(alpha);
  diff_model.newsize(vecdim,1);chk_set_init(diff_model,1); 
  for(Integer i=0;i<vecdim;i++){
    diff_model(i)=-globalbundle[unsigned(startindex_bundle+i)].evaluate(-1,y);
  }
  diff_model(0)+=trace_dual;
  return 0;
}

  /// If mu is not zero, always add the centering term for this mu as well;
int SOCIPBundleBlock::set_modelrhs(Matrix& globalrhs, 
				 Real rhsmu,
				 Real rhscorr,
				 Integer startindex_model)
{
  mat_xey(vecdim,globalrhs.get_store()+startindex_model,diff_model.get_store());
  return add_muxinv(globalrhs,startindex_model,rhsmu,rhscorr,true);
}

  ///add the "scaled" minorant outer products to globalsys, where the correct minroants start at the given index
int SOCIPBundleBlock::add_BtinvsysB(Symmatrix& globalsys,
			    const MinorantBundle& globalbundle,
			    Integer startindex_bundle)
{
  tmpmat.newsize(globalsys.rowdim(),vecdim); chk_set_init(tmpmat,1);
  for (Integer i=0;i<vecdim;i++){
    Real dummy;
    globalbundle[unsigned(startindex_bundle+i)].get_minorant(dummy,tmpmat,i);
  }
  return add_AxizinvAt(tmpmat,globalsys);
}



   /// glob_lowrank(:,...)=(Btsyssqrtinv), trafotrace(...)=syssqrtinv*trace
  int SOCIPBundleBlock::Schur_transform_bundle(Matrix& glob_lowrank,
					     MinorantBundle& globalbundle,
					     Integer startindex_bundle,
					     Matrix& trafotrace,
					     Integer startindex_trace)
  {
    if (f.dim()!=vecdim){
      compute_NTscaling();
    }

    assert(trafotrace.rowdim()>=startindex_trace+vecdim);
    assert(glob_lowrank.coldim()>=startindex_bundle+vecdim);

    //Finv applied to e_0
    mat_xeya(vecdim,trafotrace.get_store()+startindex_trace,f.get_store(),-1./omega);
    trafotrace(startindex_trace)*=-1;

    if (Bt.coldim()!=vecdim){
      Bt.newsize(glob_lowrank.rowdim(),0);
      get_Bt(Bt,0,globalbundle,startindex_bundle);
    }
    tmpmat.init(Bt,1.,1);   //transposed!
    apply_Finv(tmpmat);

    for (Integer i=0;i<vecdim;i++){  
      mat_xey(glob_lowrank.rowdim(),glob_lowrank.get_store()+(startindex_bundle+i)*glob_lowrank.rowdim(),1,tmpmat.get_store()+i,tmpmat.rowdim());
    }
    
    return 0;
  }

  /** @brief add diag(Bt*sqrt(invsys)*(I-lambda*trvec*trvec')*sqrt(invsys)*B) to diagonal 
   */
  int SOCIPBundleBlock::add_bundle_xizinv_diagonal(Matrix& diagonal,
						   Matrix& ipBtrvec,
						   MinorantBundle& globalbundle,
						   Integer startindex_bundle,
						   const Matrix& trafotrace,
						   Integer startindex_trace)
  {
    assert(vecdim==dim_bundle());
    if (f.dim()!=vecdim){
      compute_NTscaling();
    }

    if (Bt.coldim()!=vecdim){
      Bt.newsize(diagonal.rowdim(),0); //0 columns signals initialization
      get_Bt(Bt,0,globalbundle,startindex_bundle);
    }   

    assert(trafotrace.rowdim()>=startindex_trace+vecdim);
    for(Integer i=0;i<diagonal.rowdim();i++){
      tmpvec=Bt.row(i);
      tmpvec.transpose();
      apply_Finv(tmpvec);
      diagonal(i)+=ip(tmpvec,tmpvec);
      if (ipBtrvec.rowdim()>0)
	ipBtrvec(i)+=mat_ip(vecdim,trafotrace.get_store()+startindex_trace,tmpvec.get_store());
    }
    
    return 0;
  }
					 

int SOCIPBundleBlock::add_pcsubspace(Matrix& lowrank,
			     Matrix& sigma_guess,
			     const Matrix& Diag_inv,
			     Real minval,
			     Real diaginvval,
			     Matrix & minus_trmult,
			     Real schur_trace,
			     MinorantBundle& globalbundle,
			     Integer startindex_bundle)
  {
    //project system matrix onto the two relevant vectors f and e1
    //and compute the eigenvectors for this projection
    if (vecdim==0)
      return 0;

    if (cb_out())
      get_out()<<"\n****  WARNING:  SOCIPBundleBlock::add_pcsubspace() has not yet been tested in any way"<<std::endl;

    if (f.dim()!=vecdim){
      compute_NTscaling();
    }

    Matrix localsqrBnorms;
    if (diaginvval>0.) {
      if (sqrBnorms.rowdim()==0){
	sqrBnorms.init(vecdim,1,0.);
	for (Integer i=0;i<vecdim;i++){
	  sqrBnorms(i)=globalbundle[unsigned(startindex_bundle+i)].norm_squared();
	}
      }
      localsqrBnorms.init(sqrBnorms,diaginvval);
    }
    else {
      localsqrBnorms.init(vecdim,1,0.);
      for (Integer i=0;i<vecdim;i++){
	localsqrBnorms(i)=globalbundle[unsigned(startindex_bundle+i)].norm_squared(&Diag_inv);
      }
    }

    Real f0=f(0);
    Real favgnrm=0.;
    for(Integer i=0;i<vecdim;i++)
      favgnrm+=sqr(f(i))*localsqrBnorms(i);
    favgnrm=std::sqrt(favgnrm)/(sqr(omega)*(2*f0*f0-1.));

    Real corrfac=(1-std::sqrt(1-(2*sqr(f0)-1.)/schur_trace))/(2*sqr(f0)-1);
    //Real trFinv2tr=(8*sqr(f0)*(sqr(f0)-1.)+1.)/(sqr(sqr(omega)));
    
    // first guess value for  Finv*(I-trvec*trevc/trfac)^.5 times e_0 and barf via averaged norm
    Real d0=favgnrm*std::sqrt(sqr(f0*(1-corrfac*(2*f0*f0-1)/sqr(omega)))+(f0*f0)*sqr(2*corrfac*sqr(f0/omega)-1.))/omega;
    Real d1=favgnrm*std::fabs(1-f0)*std::sqrt(sqr(1+f0-corrfac*(1+f0)*(2*f0*f0-1)/sqr(omega))+(f0*f0-1.)*(sqr(1.-corrfac*(1+f0)*2*f0/sqr(omega)))/(f0*f0-1))/omega;

    // reserve the memory
    Integer cnt=Integer(d0>minval)+Integer(d1>minval);
    if (cnt==0)
      return 0;
    Integer nextcol=lowrank.coldim();
    lowrank.enlarge_right(cnt,0.);
    sigma_guess.enlarge_below(cnt,0.);
      
    // append the corresponding vectors
    
    if (d0>minval){
      if (cb_out(2)){
	get_out()<<" S0("<<d0<<")";
      }
      sigma_guess.concat_below(d0);
      Real dummy;
      globalbundle[unsigned(startindex_bundle)].get_minorant(dummy,lowrank,nextcol,f0/omega,true);
      for (Integer i=1;i<vecdim;i++){
	if (f(i)!=0.)
	  globalbundle[unsigned(startindex_bundle+i)].get_minorant(dummy,lowrank,nextcol,-f(i)/omega,true);
      }
      minus_trmult.concat_below(f0/omega);   //inner product with tracevector
      nextcol++;
    }

    if (d1>minval){
      if (cb_out(2)){
	get_out()<<" Sf("<<d1<<")";
      }
      sigma_guess.concat_below(d1);
      Real d=(1.-f0)/omega/std::sqrt(f0*f0-1);
      if (d0>minval) {
	//use previous column
	Real dummy;
	globalbundle[unsigned(startindex_bundle)].get_minorant(dummy,lowrank,nextcol,d*(1.+f0)-f0/d,true);
	mat_xpeya(lowrank.rowdim(),lowrank.get_store()+nextcol*lowrank.rowdim(),lowrank.get_store()+(nextcol-1)*lowrank.rowdim(),omega/d);
      }
      else {
	//compute from scratch
	Real dummy;
	globalbundle[unsigned(startindex_bundle)].get_minorant(dummy,lowrank,nextcol,d*(1.+f0),true);
	for (Integer i=1;i<vecdim;i++){
	  if (f(i)!=0.)
	    globalbundle[unsigned(startindex_bundle+i)].get_minorant(dummy,lowrank,nextcol,-d*f(i),true);
	}
      }
      minus_trmult.concat_below(d*(1.+f0));   //inner product with tracevector
      nextcol++;
    }

    assert(nextcol==lowrank.coldim());
    		 
    return 0;
  }

    /// add bundle*sqrt(inv(xiz))*subspace to glob_lowrank with bundle(:,si_bundle+1:si_bundle+dim_bundle()-1) and subspace(si_subsp:si_subsp+dim_bundle,:); sqrt(inv(xiz)) has to match that used in set_xizinvsqrt_trace()
  int SOCIPBundleBlock::add_bundle_xizinvsqrt_projection(Matrix& glob_lowrank,
						       Matrix& subspace,
						       Integer startindex_subspace,
						       MinorantBundle& globalbundle,
						       Integer startindex_bundle)
  {
    assert(vecdim==dim_bundle());
    if (f.dim()!=vecdim){
      compute_NTscaling();
    }

    if (startindex_subspace>=0){
      tmpmat.newsize(vecdim,subspace.coldim()); chk_set_init(tmpmat,1);
      for(Integer i=0;i<subspace.coldim();i++){
	mat_xey(vecdim,tmpmat.get_store()+i*vecdim,subspace.get_store()+i*subspace.rowdim()+startindex_subspace+i);
      }
      apply_Finv(tmpmat);
      B_times(tmpmat,glob_lowrank,1.,1.,1,0,0,globalbundle,startindex_bundle);
    }
    else {
      Matrix tmpmat(bundle_dim,glob_lowrank.coldim(),0.);
      B_times(glob_lowrank,tmpmat,1.,0.,0,0,0,globalbundle,startindex_bundle);
      apply_Finv(tmpmat);
      subspace.concat_right(tmpmat,1);
    }

	
    return 0;
  }


  /// out_vec+=BtinvsysB*in_vec
  int SOCIPBundleBlock::add_BtinvsysB_times(const Matrix& in_vec,
					  Matrix& out_vec,
					  Real zeta_inval,
					  Real* zeta_outval,
					  MinorantBundle& globalbundle,
					  Integer startindex_bundle)
  {
    tmpvec.init(vecdim,1,0.);
    B_times(in_vec,tmpvec,1.,0.,0,0,0,globalbundle,startindex_bundle);
    tmpvec(0)-=zeta_inval;
    SOCIPBlock::apply_xizinv(tmpvec,0);
    B_times(tmpvec,out_vec,1.,1.,1,0,0,globalbundle,startindex_bundle);
    if (zeta_outval)
      (*zeta_outval)-=tmpvec(0);
    
    return 0;
  }

  /// compute dx (and then dz) given step_y and step_trdual on basis of the last rhs computed for the model block
  int SOCIPBundleBlock::set_dx_xizsolvestep(const Matrix& step_y,
				 const Real step_trdual,
				 MinorantBundle& globalbundle,
				 Integer startindex_bundle)
  {
    tmpvec.init(diff_model);
    tmpvec(0)+=step_trdual;
    B_times(step_y,tmpvec,-1.,1.,0,0,0,globalbundle,startindex_bundle);

    //use of tmpvec is ok, because not employed in LinBlock::set_dx_xizsolverhs;
    return SOCIPBlock::set_dx_xizsolverhs(tmpvec,0);
  }
				  

  /// after the bundle subproblem is solved, this retrieves the local linear solution vector; if linx_activity is set, the values between zero and one indicate the guess on the coefficients activity level 
  int SOCIPBundleBlock::get_socx(Matrix& socx,
			       Real* socx_activity,
			       Real /* trace_rhs */,
			       bool /* cautious */) const
  {
    socx.init(x);

    if (socx_activity){
      if (vecdim==0) {
	*socx_activity=0.;
	return 0;
      }
      
      //determine activity
      Real x0=x(0);
      Real z0=z(0);
      Real oldx0=oldx(0);
      Real oldz0=oldz(0);
      
      Real tapia_factor=(old_mu>0.)?mu/old_mu:1.;
      
      *socx_activity = x0/oldx0;  
      Real dual_tapia= z0/oldz0;
      
      if (cb_out(1)){
	get_out()<<" SOCIPBB: last_alpha="<<last_alpha;
	get_out()<<" mu="<<mu<<" old_mu="<<old_mu<<" tapia_factor="<<tapia_factor<<std::endl;
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
	if (x0>0.1*std::sqrt(mu)*z0)
	  *socx_activity=1.;
	else
	  *socx_activity=0.;
      }
    }
    
    return 0;
  }




}
