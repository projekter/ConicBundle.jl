/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/NNCIPBundleBlock.cxx
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


#include "NNCIPBundleBlock.hxx"

using namespace CH_Matrix_Classes;

namespace ConicBundle {

  int NNCIPBundleBlock::find_inactive_indices(Indexmatrix& inactive,
					    Real trace_rhs,
					    bool cautious) const
  {
    inactive.newsize(vecdim,1);
    inactive.init(0,0,Integer(0));
    
    Real tapia_factor=(old_mu>0.)?mu/old_mu:1.;
    
    for (Integer i=0;i<vecdim;i++){
      if (x(i)>trace_rhs*(cautious?1e-6:1e-6)){
	//too big
	continue;
      }
      //check for relation to the dual
      if (z(i)<(cautious?1e6:1e3)*x(i)){
	continue;
      }
      if (x(i)<trace_rhs*(cautious?1e-10:1e-10)){
	//this is certainly not active
	inactive.concat_below(i);
	continue;
      }
      if ((mu<0.01)&&(x(i)<0.01*std::sqrt(mu))&&(z(i)>std::sqrt(mu))){
	inactive.concat_below(i);
	continue;
      }
      //check whether reduction is faster than for the dual
      Real primal_tapia=x(i)/oldx(i);
      Real dual_tapia=z(i)/oldz(i);
      if (tapia_factor<.8){
	if(
	   (primal_tapia<0.8)&&
	   //((1.1*(primal_tapia-tapia_factor)>(dual_tapia-tapia_factor)))
	   (primal_tapia<.1*dual_tapia)
	   )
	  //goes to zero faster than the dual value
	  inactive.concat_below(i);
      }
    }

    if (cb_out(1)){
      get_out()<<" NNCIPBB: last_alpha="<<last_alpha;
      get_out()<<" mu="<<mu<<" old_mu="<<old_mu;
      get_out()<<" tapia_factor="<<tapia_factor<<" cautious="<<cautious<<" inactive="<<transpose(inactive);
    }


    return 0;
  }
  
  void NNCIPBundleBlock::clear(Integer dim)
  {
    NNCIPBlock::clear(dim);
    bundle_dim=dim;
    diff_model.init(0,0,0.);
    map_to_old.init(Range(0,dim-1));
    
    //oracle_data=0;
  }
  
  NNCIPBundleBlock::NNCIPBundleBlock(Integer dim,CBout* cb,int cbinc):
    CBout(cb,cbinc),InteriorPointBundleBlock(cb,cbinc),NNCIPBlock(dim,cb,cbinc)
  {clear(dim);}

  NNCIPBundleBlock::~NNCIPBundleBlock()
  {}


  InteriorPointBundleBlock* NNCIPBundleBlock::clone()
  {
    NNCIPBundleBlock* p=new NNCIPBundleBlock(0,this,0);
    p->copy_from(this);

    return p;
  }
  
  int NNCIPBundleBlock::copy_from(InteriorPointBundleBlock* inp)
  {
    NNCIPBundleBlock* p=dynamic_cast<NNCIPBundleBlock*>(inp);
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
    xiz=p->xiz;
    compl_rhs=p->compl_rhs;
    last_rhs_mu=p->last_rhs_mu;
    mu=p->mu;
    old_mu=p->old_mu;
    last_alpha=p->last_alpha;
    oldx=p->oldx;
    oldz=p->oldz;
    tmpvec=p->tmpvec;
    tmpmat=p->tmpmat;
    
    map_to_old=p->map_to_old;

    return 0;
  }
  

  Real NNCIPBundleBlock::evaluate_trace_x()
  {return sum(x);}

  Real NNCIPBundleBlock::evaluate_trace_z()
  {return sum(z);}

  Real NNCIPBundleBlock::evaluate_trace_dx() 
  {return sum(dx);}

  

    /// C=beta*C+alpha*B*B where B and B may be transposed; carry out the model part of this beginning at startindex_model and beta for the part, that is added to (the calling routine has to make sure beta is not executed repeatedly if the same part is affected by other models as well)
Matrix& NNCIPBundleBlock::B_times(const Matrix& A,
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
    globalbundle[unsigned(startindex_bundle+map_to_old(i))].left_genmult(A,C,alpha,1.,(Btrans==0),Atrans,startindex_model++);
  }
  return C;
}
						   
  /// C=beta*C+alpha*A*B where B and B may be transposed; carry out the model part of this beginning at startindex_model 
Matrix& NNCIPBundleBlock::times_B(const Matrix& A,
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
    globalbundle[unsigned(startindex_bundle+map_to_old(i))].right_genmult(A,C,alpha,1.,Atrans,(Btrans==0),startindex_model++);
  }
  return C;
}

  ///add the main diagonal block tranpose(projection)*diagvec*projection to bigS starting at startindex
Symmatrix& NNCIPBundleBlock::add_BDBt(const Matrix& diagvec,
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
      MinorantPointer& mp=globalbundle[unsigned(startindex_bundle+map_to_old(i))];
      for(Integer j=i+startindex_model;j<Bt.coldim();j++){
	bigS(startindex+startindex_model+i,startindex+j)-=mp.ip(Bt,&diagvec,j*Bt.rowdim());
      }
    }
  }
  else {
    for (Integer i=0;i<vecdim;i++){
      MinorantPointer& mp=globalbundle[unsigned(startindex_bundle+map_to_old(i))];
      for(Integer j=i+startindex_model;j<Bt.coldim();j++){
	bigS(startindex+startindex_model+i,startindex+j)+=mp.ip(Bt,&diagvec,j*Bt.rowdim());
      }
    }
  }
  return bigS;
}

  /// get the current matrix for the coupling matrix Bt in the first row of blocks
Matrix& NNCIPBundleBlock::get_Bt(Matrix& globBt,
			       Integer startindex_model,
			       MinorantBundle& globalbundle,
			       Integer startindex_bundle)
{
  assert((&Bt!=&globBt)||(startindex_model==0));
  if (Bt.coldim()!=vecdim){
    Bt.newsize(globBt.rowdim(),vecdim); chk_set_init(Bt,1);
    Boffset.newsize(vecdim,1); chk_set_init(Boffset,1);
    for (Integer i=0;i<vecdim;i++){
      globalbundle[unsigned(startindex_bundle+map_to_old(i))].get_minorant(Boffset(i),Bt,i);
    }
  }
  if (&Bt != &globBt){
    assert(Bt.rowdim()==globBt.rowdim());
    assert(globBt.coldim()>=startindex_model+vecdim);
    mat_xey(vecdim*Bt.rowdim(),globBt.get_store()+startindex_model*globBt.rowdim(),Bt.get_store());
  }
  return globBt;
}
  
  /// adds opB transposed times modelx (without constant affine term) to the arguments
int NNCIPBundleBlock::add_modelx_aggregate(Real& val,
				   Matrix& vec,
				   MinorantBundle& globalbundle,
				   Integer startindex_bundle)
{
  for (Integer i=0;i<vecdim;i++){
    globalbundle[unsigned(startindex_bundle+map_to_old(i))].get_minorant(val,vec,0,x(i),true);
  }
  return 0;
}

  /// set the model violation for the current system solution 
int NNCIPBundleBlock::get_sysviol_model(Matrix& sysviol_model,
				      Integer startindex_model,
				      const Matrix& dy,
				      const Real deltatrdual,
				      MinorantBundle& globalbundle,
				      Integer startindex_bundle)
{
  for(Integer i=0;i<vecdim;i++){
    sysviol_model(startindex_model+i)=globalbundle[unsigned(startindex_bundle+map_to_old(i))].evaluate(-1,dy,false)-diff_model(i)-deltatrdual+z(i)+dz(i);
  }
  return 0;
}

  /// set z to the slack of the bundle and return a value>=0 that needs to be added to make it feasible
int NNCIPBundleBlock::set_bundle_z(const Matrix& y,
		    MinorantBundle& globalbundle,
		    Integer startindex_bundle,
		    Real& add_center_value)
{
  diff_model.newsize(vecdim,1);chk_set_init(diff_model,1); 
  for(Integer i=0;i<vecdim;i++){
    diff_model(i)=-globalbundle[unsigned(startindex_bundle+map_to_old(i))].evaluate(-1,y);
  }
  return NNCIPBlock::set_z(diff_model,0,add_center_value);
}

int NNCIPBundleBlock::add_trace(Matrix& vec,Real trace_dual,Integer startindex) 
{
  mat_xpea(vecdim,vec.get_store()+startindex,trace_dual);
  return 0;
}

    /// set the trace premultiplied by sqrt(inv(xiz)) in vec[startindex+0,...,startindex+dim_bundle()-1]
  int NNCIPBundleBlock::set_xizinvsqrt_trace(Matrix& vec,
					   Integer startindex)
  {
    assert(vecdim==bundle_dim);    
    if (xiz.dim()!=vecdim){
      compute_NTscaling();
    }

    for(Integer i=0;i<vecdim;i++){
      vec(startindex+i)=std::sqrt(1./xiz(i));
    }
    return 0;
  }
  

///add trace_dual*trace to diff_model for the right hand side (negative of current model violation)
int NNCIPBundleBlock::add_trace_to_diff_model(Real trace_dual)
{
  diff_model+=trace_dual;
  return 0;
}

///return the squared Euclidean norm of the dual model violation  
Real NNCIPBundleBlock::dualviol_2normsqr()
{
  tmpvec.init(diff_model);
  tmpvec-=z;
  return ip(tmpvec,tmpvec);
}

  /// move to (x+alpha*dx, z+alpha*dz), update diff_model and possibly reduce the model size if some part is below zero_threshold
int NNCIPBundleBlock::do_bundle_step(Real alpha,
			    const Matrix& y,
			    MinorantBundle& globalbundle,
			    Integer startindex_bundle,
			    Real tracedual,
			    Real /* trace_rhs */)
{
  NNCIPBlock::do_step(alpha);

  // check for size reductions due to inactivity   //seemed useless
  // Indexmatrix delind;
  // find_inactive_indices(delind,trace_rhs,true);
  // //if (delind.rowdim()>0){
  // if (false){
  //   if (cb_out(2))
  //     get_out()<<" delete "<<delind<<std::endl;
  //   //std::cout<<" delete "<<delind<<std::endl;
  //   x.delete_rows(delind,true);
  //   z.delete_rows(delind,true);
  //   oldx.delete_rows(delind,true);
  //   oldz.delete_rows(delind,true);
  //   map_to_old.delete_rows(delind,true);
  //   vecdim=x.rowdim();
  // }
  
  diff_model.init(vecdim,1,tracedual);
  for(Integer i=0;i<vecdim;i++){
    diff_model(i)-=globalbundle[unsigned(startindex_bundle+map_to_old(i))].evaluate(-1,y);
  }
  return 0;
}

  /// If mu is not zero, always add the centering term for this mu as well;
int NNCIPBundleBlock::set_modelrhs(Matrix& globalrhs, 
				 Real rhsmu,
				 Real rhscorr,
				 Integer startindex_model)
{
  mat_xey(vecdim,globalrhs.get_store()+startindex_model,diff_model.get_store());
  return add_muxinv(globalrhs,startindex_model,rhsmu,rhscorr,true);
}

  ///add the "scaled" minorant outer products to globalsys, where the correct minroants start at the given index
int NNCIPBundleBlock::add_BtinvsysB(Symmatrix& globalsys,
			    const MinorantBundle& globalbundle,
			    Integer startindex_bundle)
{
  tmpmat.newsize(globalsys.rowdim(),vecdim); chk_set_init(tmpmat,1);
  for (Integer i=0;i<vecdim;i++){
    Real dummy;
    globalbundle[unsigned(startindex_bundle+map_to_old(i))].get_minorant(dummy,tmpmat,i);
  }
  return add_AxizinvAt(tmpmat,globalsys);
}


   /// glob_lowrank(:,...)=(Btsyssqrtinv), trafotrace(...)=syssqrtinv*trace
  int NNCIPBundleBlock::Schur_transform_bundle(Matrix& glob_lowrank,
					     MinorantBundle& globalbundle,
					     Integer startindex_bundle,
					     Matrix& trafotrace,
					     Integer startindex_trace)
  {
    if (xiz.dim()!=vecdim){
      compute_NTscaling();
    }

    assert(trafotrace.rowdim()>=startindex_trace+vecdim);
    assert(glob_lowrank.coldim()>=startindex_bundle+vecdim);

    tmpmat.newsize(glob_lowrank.rowdim(),1);chk_set_init(tmpmat,1);
    for(Integer i=0;i<vecdim;i++){
      Real sqrtxizinv=1./std::sqrt(xiz(i));
      trafotrace(startindex_trace+i)=sqrtxizinv;
      Real dummy;
      globalbundle[unsigned(startindex_bundle+map_to_old(i))].get_minorant(dummy,glob_lowrank,startindex_bundle+i,sqrtxizinv);
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
  int NNCIPBundleBlock::add_bundle_xizinv_diagonal(CH_Matrix_Classes::Matrix& diagonal,
						   CH_Matrix_Classes::Matrix& ipBtrvec,
						   MinorantBundle& globalbundle,
						   CH_Matrix_Classes::Integer startindex_bundle,
						   const CH_Matrix_Classes::Matrix& trafotrace,
						   CH_Matrix_Classes::Integer startindex_trace)
  {
    if (xiz.dim()!=vecdim){
      compute_NTscaling();
    }

    assert(trafotrace.rowdim()>=startindex_trace+vecdim);
    Real scaleval;
    Minorant* mp;
    Integer nz;
    const Real* coeffs;
    const Integer *ind;
    for(Integer i=0;i<vecdim;i++){
      Real xizinv=1./xiz(i);
      globalbundle[unsigned(startindex_bundle+map_to_old(i))].get_scaleval_and_minorant(scaleval,mp);
      if ((scaleval==0.)||(mp==0))
	continue;
      mp->get_coeffs(nz,coeffs,ind);
      if ((nz==0)||(coeffs==0))
	continue;
      xizinv*=scaleval*scaleval;
      if (ind==0){
	for(Integer j=0;j<nz;j++){
	  diagonal(j)+=xizinv*sqr(coeffs[j]);
	}
      }
      else {
	for(Integer j=0;j<nz;j++){
	  diagonal(ind[j])+=xizinv*sqr(coeffs[j]);
	}
      }
      //accumulate lambda part in tmpvec
      if (ipBtrvec.rowdim()>0){
	xizinv=std::sqrt(xizinv)*trafotrace(startindex_trace+i);
	if (ind==0){
	  for(Integer j=0;j<nz;j++){
	    ipBtrvec(j)+=xizinv*(coeffs[j]);
	  }
	}
	else {
	  for(Integer j=0;j<nz;j++){
	    ipBtrvec(ind[j])+=xizinv*coeffs[j];
	  }
	}
      }
    } //end for

    
    return 0;
  }
					 

  int NNCIPBundleBlock::add_pcsubspace(Matrix& lowrank,
				       Matrix& sigma_guess,
				       const Matrix& Diag_inv,
				       Real minval,
				       Real diaginvval,
				       Matrix & minus_trmult,
				       Real schur_trace,
				       MinorantBundle& globalbundle,
				       Integer startindex_bundle)
  {
    if (xiz.dim()!=vecdim){
      compute_NTscaling();
    }
    assert(minval>0.);

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

      
    //count the number of values exceeding the threshold
    Real threshold=minval*minval;
    Integer cnt=0;
    tmpvec.newsize(vecdim,1); chk_set_init(tmpvec,1);
    for (Integer i=0;i<vecdim;i++){
      Real d=localsqrBnorms(i)*(1/xiz(i)*(1.-1./xiz(i)/schur_trace));
      tmpvec(i)=d;
      if (d>=threshold)
	cnt++;
    }

    if (cnt==0)
      return 0;

    //reserve the memory
    Integer nextcol=lowrank.coldim();
    lowrank.enlarge_right(cnt,0.);
    sigma_guess.enlarge_below(cnt,0.);

    //append large enough vectors
    for(Integer i=0;i<vecdim;i++){
      if (tmpvec(i)>=threshold){
	if (cb_out(2)){
	  get_out()<<" N("<<i<<";"<<tmpvec(i)<<")";
	}
	sigma_guess(nextcol)=std::sqrt(tmpvec(i));
	Real d=1./std::sqrt(xiz(i));
	Real dummy;
	globalbundle[unsigned(startindex_bundle+i)].get_minorant(dummy,lowrank,nextcol,d,true);
	minus_trmult.concat_below(d);
	nextcol++;
      }
    }
    assert(nextcol==lowrank.coldim());
    
    return 0;
  }

    /// add bundle*sqrt(inv(xiz))*subspace to glob_lowrank with bundle(:,si_bundle+1:si_bundle+dim_bundle()-1) and subspace(si_subsp:si_subsp+dim_bundle,:); sqrt(inv(xiz)) has to match that used in set_xizinvsqrt_trace()
  int NNCIPBundleBlock::add_bundle_xizinvsqrt_projection(Matrix& glob_lowrank,
						       Matrix& subspace,
						       Integer startindex_subspace,
						       MinorantBundle& globalbundle,
						       Integer startindex_bundle)
  {
    assert(vecdim==dim_bundle());
    if (xiz.dim()!=vecdim){
      compute_NTscaling();
    }

    if (startindex_subspace>=0) {
      for(Integer i=0;i<vecdim;i++){
	MinorantPointer& mp=globalbundle[unsigned(startindex_bundle+i)];
	Real scalval=std::sqrt(1./xiz(i));
	for(Integer j=0;j<subspace.coldim();j++){
	  Real d=subspace(startindex_subspace+i,j)*scalval;
	  if (d!=0.){
	    Real dummy;
	    mp.get_minorant(dummy,glob_lowrank,j,d,true);
	  }
	}
      }
    }
    else {
      Matrix tmpmat(glob_lowrank.coldim(),bundle_dim,0.);
      times_B(glob_lowrank,tmpmat,1.,0.,1,1,0,globalbundle,startindex_bundle);
      tmpvec=xiz;
      tmpvec.inv();
      tmpvec.sqrt();
      tmpmat.scale_cols(tmpvec);
      subspace.concat_right(tmpmat);
    }

    return 0;
  }

  /// out_vec+=BtinvsysB*in_vec
  int NNCIPBundleBlock::add_BtinvsysB_times(const Matrix& in_vec,
					  Matrix& out_vec,
					  Real zeta_inval,
					  Real* zeta_outval,
					  MinorantBundle& globalbundle,
					  Integer startindex_bundle)
  {
    if (xiz.dim()!=vecdim){
      compute_NTscaling();
    }

    tmpvec.init(vecdim,1,0.);
    B_times(in_vec,tmpvec,1.,0.,0,0,0,globalbundle,startindex_bundle);
    tmpvec-=zeta_inval;
    tmpvec/=xiz;
    B_times(tmpvec,out_vec,1.,1.,1,0,0,globalbundle,startindex_bundle);
    if (zeta_outval)
      (*zeta_outval)-=sum(tmpvec);
    
    return 0;
  }

  /// compute dx (and then dz) given step_y and step_trdual on basis of the last rhs computed for the model block
  int NNCIPBundleBlock::set_dx_xizsolvestep(const Matrix& step_y,
				 const Real step_trdual,
				 MinorantBundle& globalbundle,
				 Integer startindex_bundle)
  {
    tmpvec.init(diff_model);
    tmpvec+=step_trdual;
    B_times(step_y,tmpvec,-1.,1.,0,0,0,globalbundle,startindex_bundle);

    //use of tmpvec is ok, because not employed in NNCIPBlock::set_dx_xizsolverhs;
    return NNCIPBlock::set_dx_xizsolverhs(tmpvec,0);
  }
				  


  /// after the bundle subproblem is solved, this retrieves the local linear solution vector; if nncx_activity is set, the values between zero and one indicate the guess on the coefficients activity level 
   int NNCIPBundleBlock::get_nncx(Matrix& nncx,
				Matrix* nncx_activity,
				Real trace_rhs,
				bool cautious) const
  {
    if (vecdim==bundle_dim){
      nncx.init(vecdim,1,x.get_store());
      if (nncx_activity){
	Indexmatrix inactive;
	find_inactive_indices(inactive,trace_rhs,cautious);
	nncx_activity->init(vecdim,1,1.);
	for (Integer i=0;i<inactive.rowdim();i++){
	  (*nncx_activity)(inactive(i))=0.;
	}
      }
    }
    else {
      nncx.init(bundle_dim,1,0.);
      for (Integer i=0;i<vecdim;i++){
	nncx(map_to_old(i))=x(i);
      }
      if (nncx_activity){
	nncx_activity->init(bundle_dim,1,0.);
	for(Integer i=0;i<vecdim;i++){
	  (*nncx_activity)(map_to_old(i))=1.;
	}
      }			    
    }
    
    return 0;
  }
    



}
