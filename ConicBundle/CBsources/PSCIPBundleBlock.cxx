/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/PSCIPBundleBlock.cxx
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


#include "PSCIPBundleBlock.hxx"

using namespace CH_Matrix_Classes;

namespace ConicBundle {

  QPPSCOracleDataObject::~QPPSCOracleDataObject()
  {}
  
  
  int PSCIPBundleBlock::get_growth(const CH_Matrix_Classes::Matrix& lamX,
				   const CH_Matrix_Classes::Matrix& vecX,
				   CH_Matrix_Classes::Real& growthrate,
				   CH_Matrix_Classes::Matrix& primalgrowth,
				   CH_Matrix_Classes::Matrix& dualgrowth) const
  {
    growthrate=ip(X,Z)/ip(oldX,oldZ);
    primalgrowth.newsize(lamX.rowdim(),1); chk_set_init(primalgrowth,1);
    dualgrowth.newsize(lamX.rowdim(),1); chk_set_init(dualgrowth,1);
    
    for(Integer i=0;i<lamX.rowdim();i++){
      tmpvec=vecX.col(i);
      Real Xval=lamX(i);
      Real Zval=ip(tmpvec,Z*tmpvec);
      Real oldXval=ip(tmpvec,oldX*tmpvec);
      Real oldZval=ip(tmpvec,oldZ*tmpvec);
      primalgrowth(i)=Xval/oldXval;
      dualgrowth(i)=Zval/oldZval;
    }
    return 0;
  }
  

  Integer PSCIPBundleBlock::active_rank(const Matrix& eigX,
				      const Matrix& vecX,
				      Real trace_rhs,
				      bool cautious) const
  {
    Real tapia_factor=(old_mu>0.)?mu/old_mu:1.;
    Real tapia_factor_guess=ip(X,Z)/ip(oldX,oldZ);
    Real X0val,Z0val;

    Real tapiaratio_threshold=1.;  //=.95;
    //Real tapia_threshold=min(1.,std::pow(tapia_factor_guess,1./3.));
    Real tapia_threshold=std::sqrt(tapia_factor_guess);

    Matrix tapiaratio(eigX.rowdim(),1); chk_set_init(tapiaratio,1);
    for(Integer i=0;i<eigX.rowdim();i++){
      tmpvec=vecX.col(i);
      Real Xval=eigX(i);
      Real Zval=ip(tmpvec,Z*tmpvec);
      Real oldXval=ip(tmpvec,oldX*tmpvec);
      Real oldZval=ip(tmpvec,oldZ*tmpvec);
      Real primal_tapia=Xval/oldXval;
      Real dual_tapia=Zval/oldZval;
      tapiaratio(i)=primal_tapia/dual_tapia;
    }
    Indexmatrix sind;
    sortindex(tapiaratio,sind);

    Integer rank=rowdim;
    if (cb_out(1)){
      get_out()<<" PSCIPBB: last_alpha="<<last_alpha;
      get_out()<<" mu="<<mu<<" old_mu="<<old_mu;
      get_out()<<" tapia_factor="<<tapia_factor;
      get_out()<<"("<<tapia_factor_guess<<";thr="<<tapia_threshold<<")";
      get_out()<<" cautious="<<cautious<<std::endl;
      X0val=eigX(0);
      tmpvec=vecX.col(0);
      Z0val=ip(tmpvec,Z*tmpvec);
      for(Integer i=0;i<eigX.dim();i++){
	tmpvec=vecX.col(i);
	Real Xval=eigX(i);
	Real Zval=ip(tmpvec,Z*tmpvec);
	Real oldXval=ip(tmpvec,oldX*tmpvec);
	Real oldZval=ip(tmpvec,oldZ*tmpvec);
	Real primal_tapia=Xval/oldXval;
	Real dual_tapia=Zval/oldZval;	
	get_out()<<" "<<i<<":"<<Xval<<","<<Zval;
	get_out()<<" "<<Real(Z0val/X0val>min(1e-3,mu)*Zval/Xval);
	get_out()<<" : "<<oldXval<<","<<oldZval;
	get_out()<<"   ("<<primal_tapia<<","<<dual_tapia<<")";
	get_out()<<" "<<(primal_tapia/dual_tapia)<<"("<<sind(eigX.rowdim()-i-1)<<") ";
	get_out()<<"[b"<<(Xval>(cautious?1e-6:.1)*trace_rhs);
	get_out()<<",c"<<(Xval<1e-10*trace_rhs);
	get_out()<<",b"<<((Zval<(cautious?1e3:1e-3)*Xval)&&(primal_tapia>2.*dual_tapia));
	get_out()<<",c"<<((mu<0.01)&&(Xval<1e-2*std::sqrt(mu))&&(Zval>std::sqrt(mu)));
	//	get_out()<<",b"<<(((tapia_factor<.8)&&((primal_tapia>0.8)||(primal_tapia>(cautious?.1:.9)*dual_tapia)))||((tapia_factor>=.8)&&((Z0val/X0val>std::sqrt(mu)*Zval/Xval))));
	get_out()<<",b"<<(((tapia_factor<.8)&&((primal_tapia>.9)||((primal_tapia>(cautious?0.3:0.7))&&(primal_tapia>(cautious?.9:1.5)*dual_tapia))))||((Z0val/X0val>std::sqrt(mu)*Zval/Xval)));
	get_out()<<",b"<<(primal_tapia/dual_tapia>tapiaratio_threshold);
	get_out()<<"]";
	get_out()<<"{"<<(primal_tapia>tapia_threshold)<<"}";
	// if (tapia_factor<0.99){
	//   get_out()<<"   ("<<(primal_tapia-tapia_factor)/(1.-tapia_factor);
	//   get_out()<<","<<(dual_tapia-tapia_factor)/(1.-tapia_factor)<<") ";
	//   get_out()<<" "<<Real((primal_tapia>0.8)&&(1.1*((primal_tapia-tapia_factor)/(1.-tapia_factor))>((dual_tapia-tapia_factor)/(1.-tapia_factor))));
	// }
	get_out()<<std::endl;						      
      }
    }

    if (rank>0){
      X0val=eigX(0);
      if (X0val<1e-10*trace_rhs)
	return 0;
      tmpvec=vecX.col(0);
      Z0val=ip(tmpvec,Z*tmpvec);
      
      while (--rank>=0) {
	Real Xval=eigX(rank);
	
	//a rather large eigenvalue might not allow to exclude its eigenspace yet 
	if (Xval>(cautious?1e-6:0.1)*trace_rhs){
	  if (cb_out(2)){
	    get_out()<<" "<<rank<<": "<<Xval;
	    get_out()<<"[b"<<(Xval>(cautious?1e-6:.1)*trace_rhs);
	    get_out()<<"]"<<std::endl;
	  }
	  break;
	}
	
	//if the value is really small, its eigenspace is irrelevant
	if (Xval<1e-10*trace_rhs)
	  continue;
	
	//check whether the dual signals irrelevance of this eigenspace
	
	tmpvec=vecX.col(rank);
	Real Zval=ip(tmpvec,Z*tmpvec);
	Real oldXval=ip(tmpvec,oldX*tmpvec);
	Real oldZval=ip(tmpvec,oldZ*tmpvec);
	Real primal_tapia=Xval/oldXval;
	Real dual_tapia=Zval/oldZval;
	if ((Zval<(cautious?1e3:1e-3)*Xval)&&(primal_tapia>2.*dual_tapia)){
	  if (cb_out(2)){
	    get_out()<<" "<<rank<<": "<<Xval<<","<<Zval;
	    get_out()<<" : "<<oldXval<<","<<oldZval;
	    get_out()<<"   ("<<primal_tapia<<","<<dual_tapia<<")";
	    get_out()<<"[b"<<(Xval>(cautious?1e-6:.1)*trace_rhs);
	    get_out()<<",c"<<(Xval<1e-10*trace_rhs);
	    get_out()<<",b"<<((Zval<(cautious?1e3:1e-3)*Xval)&&(primal_tapia>2.*dual_tapia))<<"]"<<std::endl;
	  }
	  break;
	}
	
	if ((mu<0.01)&&(Xval<1e-2*std::sqrt(mu))&&(Zval>std::sqrt(mu)))
	  continue;
	
	// check whether the eigenvalue is going to zero
	// if there was a reduction from old_mu to mu, did this affect the primal more than the dual? 
	
	if (tapia_factor<.8){
	  if(
	     (primal_tapia>0.9)||
	     //(primal_tapia>0.8)||
	     //((1.1*(primal_tapia-tapia_factor)>(dual_tapia-tapia_factor)))
	     //(primal_tapia>(cautious?.1:.9)*dual_tapia)
	    //does not go to zero fast and clearly slower than the dual value
	     ((primal_tapia>(cautious?0.3:0.7))&&
	      (primal_tapia>(cautious?.9:1.5)*dual_tapia))
	     ) {
	    //does not go to zero much faster than the dual value
	    if (cb_out(2)){
	      get_out()<<" "<<rank<<": "<<Xval<<","<<Zval;
	      get_out()<<" : "<<oldXval<<","<<oldZval;
	      get_out()<<"   ("<<primal_tapia<<","<<dual_tapia<<")";
	      get_out()<<"[b"<<(Xval>(cautious?1e-6:.1)*trace_rhs);
	      get_out()<<",c"<<(Xval<1e-10*trace_rhs);
	      get_out()<<",b"<<((Zval<(cautious?1e3:1e-3)*Xval)&&(primal_tapia>2.*dual_tapia));
	      get_out()<<",c"<<((mu<0.01)&&(Xval<1e-2*std::sqrt(mu))&&(Zval>std::sqrt(mu)));
	      get_out()<<",b"<<((tapia_factor<.8)&&((primal_tapia>.9)||((primal_tapia>(cautious?0.3:0.7))&&(primal_tapia>(cautious?.9:1.5)*dual_tapia))));
	      get_out()<<"]"<<std::endl;
	    }
	    break;
	  }
	}
	else {
	  //no sufficient reduction in mu, judge by increase off the ratio of dual to primal
	  if (Z0val/X0val>std::sqrt(mu)*Zval/Xval){
	    if (cb_out(2)){
	      get_out()<<" "<<rank<<": "<<Xval<<","<<Zval;
	      get_out()<<" : "<<oldXval<<","<<oldZval;
	      get_out()<<"   ("<<primal_tapia<<","<<dual_tapia<<")";
	      get_out()<<"[b"<<(Xval>(cautious?1e-6:.1)*trace_rhs);
	      get_out()<<",c"<<(Xval<1e-10*trace_rhs);
	      get_out()<<",b"<<((Zval<(cautious?1e3:1e-3)*Xval)&&(primal_tapia>2.*dual_tapia));
	      get_out()<<",c"<<((mu<0.01)&&(Xval<1e-2*std::sqrt(mu))&&(Zval>std::sqrt(mu)));
	      get_out()<<",wb"<<((tapia_factor<.8)&&((Z0val/X0val>std::sqrt(mu)*Zval/Xval)));
	      get_out()<<"]"<<std::endl;
	    }
	    break;
	  }
	}
	if (primal_tapia>tapia_threshold){
	  //if (primal_tapia/dual_tapia>tapiaratio_threshold){
	  if (cb_out(2)){
	    get_out()<<" "<<rank<<":"<<Xval<<","<<Zval;
	    get_out()<<" "<<Real(Z0val/X0val>min(1e-3,mu)*Zval/Xval);
	    get_out()<<" : "<<oldXval<<","<<oldZval;
	    get_out()<<"   ("<<primal_tapia<<","<<dual_tapia<<")";
	    get_out()<<"[b"<<(Xval>(cautious?1e-6:.1)*trace_rhs);
	    get_out()<<",c"<<(Xval<1e-10*trace_rhs);
	    get_out()<<",b"<<((Zval<(cautious?1e3:1e-3)*Xval)&&(primal_tapia>2.*dual_tapia));
	    get_out()<<",c"<<((mu<0.01)&&(Xval<1e-2*std::sqrt(mu))&&(Zval>std::sqrt(mu)));
	    get_out()<<",b"<<(((tapia_factor<.8)&&((primal_tapia>.9)||((primal_tapia>(cautious?0.3:0.7))&&(primal_tapia>(cautious?.9:1.5)*dual_tapia))))||((Z0val/X0val>std::sqrt(mu)*Zval/Xval)));
	    get_out()<<",b"<<(primal_tapia/dual_tapia>1.1);
	    get_out()<<"]";
	    get_out()<<std::endl;						      
	  }
	  break;
	}
      } // end while (--rank>=0);
      
      rank++;
    }

    if (cb_out(2)){
      get_out()<<" rank="<<rank<<std::endl;
    }

    return rank;
  }
  
  int PSCIPBundleBlock::form_B(Integer dim,
			       const MinorantBundle& globalbundle,
			       Integer startindex_bundle)
  {
    if (B.rowdim()!=vecdim){
      B.newsize(vecdim,dim); chk_set_init(B,1);
      Boffset.newsize(vecdim,1); chk_set_init(Boffset,1);
      tmpmat.newsize(dim,1); chk_set_init(tmpmat,1);
      for (Integer i=0;i<vecdim;i++){
	globalbundle[unsigned(startindex_bundle+i)].get_minorant(Boffset(i),tmpmat,0);
	mat_xey(dim,B.get_store()+i,vecdim,tmpmat.get_store(),1);
      }
      //std::cout<<" B="<<B;
      Bnzcolind.newsize(dim,1);
      Bnzcolind.init(0,0,Integer(0));
      Real* bp=B.get_store();
      for(Integer j=0;j<dim;j++){
	Integer i=vecdim;
	while ((--i>=0)&&(*bp++==0.));
	if (i>=0){
	  Bnzcolind.concat_below(j);
	  bp+=i;
	}
      }
    }
    assert(B.coldim()==dim);
    return 0;
  }

  int PSCIPBundleBlock::sveciB_times_vec(const Matrix& vec,Matrix& Bvec)
  {
    assert(vec.rowdim()==rowdim);
    assert(B.rowdim()==vecdim);
    Bvec.newsize(vec.rowdim(),B.coldim()); chk_set_init(Bvec,1); 
    Integer lastci=-1;
    Real *vBp=Bvec.get_store();
    const Real* Bp=B.get_store(); 
    const Real invsqrt2= 1./std::sqrt(2.);
    const Real* vp=vec.get_store();
    const Real* const endvp=vp+rowdim;
    for (Integer ci=0;ci<Bnzcolind.rowdim();ci++){
      //skip unused columns
      Integer nextci=Bnzcolind(ci);
      Integer diff=nextci-lastci;
      mat_xea(rowdim*diff,vBp,0.);
      vBp+=rowdim*(diff-1);
      Bp+=vecdim*(diff-1);
      lastci=nextci;
      for(Integer i=rowdim;i>0;){
	Real vd=(*vp++);
	Real dip=(*Bp++)*vd;
	vBp++;
	for(;vp<endvp;){
	  Real bd=(*Bp++)*invsqrt2;
	  dip+=bd*(*vp++);
	  (*vBp++)+=bd*vd;
	}
	vBp-=(i--);
	(*vBp++)+=dip;
	vp-=i;
      }
      vp-=rowdim;
      
      // // BEGIN TEST
      // Symmatrix tsym;
      // sveci(B.col(nextci),tsym);
      // if (norm2(tsym*vec-Bvec.col(nextci))>=1e-10*(1.+norm2(vec))){
      // 	std::cout.precision(8);
      // 	std::cout<<" sveciB_times_vec: nextci="<<nextci<<" B.col(nextci)="<<transpose(B.col(nextci));
      // 	std::cout<<" tsym="<<tsym;
      // 	std::cout<<" vec="<<vec;
      // 	std::cout<<" tsym*vec="<<transpose(tsym*vec);
      // 	std::cout<<" Bvec.col(nextci)="<<transpose(Bvec.col(nextci))<<std::endl;
      // 	std::cout<<" diff="<<transpose(tsym*vec-Bvec.col(nextci))<<std::endl;
      // }
      // assert(norm2(tsym*vec-Bvec.col(nextci))<1e-10*(1.+norm2(vec)));
      // // END TEST
    }
    Integer diff=Bvec.coldim()-1-lastci;
    if (diff>0)
      mat_xea(rowdim*diff,vBp,0.);
    
    return 0;
  }
      
  
  
  void PSCIPBundleBlock::clear(Integer dim)
  {
    PSCIPBlock::clear(dim);
    bundle_dim=vecdim;
    diff_model.init(0,0,0.);
    B.init(0,0,0.);
    lamX.init(0,1,0.);
    oracle_data=0;
  }
  
  PSCIPBundleBlock::PSCIPBundleBlock(Integer dim,CBout* cb,int cbinc):
    CBout(cb,cbinc),InteriorPointBundleBlock(cb),PSCIPBlock(dim,cb)
  {clear(dim);}

  PSCIPBundleBlock::~PSCIPBundleBlock()
  {}
  
  InteriorPointBundleBlock* PSCIPBundleBlock::clone()
  {
    PSCIPBundleBlock* p=new PSCIPBundleBlock(0,this,0);
    p->copy_from(this);

    return p;
  }
  
  int PSCIPBundleBlock::copy_from(InteriorPointBundleBlock* inp)
  {
    PSCIPBundleBlock* p=dynamic_cast<PSCIPBundleBlock*>(inp);
    if (p==0)
      return 1;

   
    bundle_dim=p->bundle_dim;
    diff_model=p->diff_model;
    B=p->B;
    Boffset=p->Boffset;
    sqrBnorms=p->sqrBnorms;
    Bnzcolind=p->Bnzcolind;
    oracle_data=p->oracle_data;
    
    rowdim=p->rowdim;
    vecdim=p->vecdim;
    X=p->X;
    Z=p->Z;
    dX=p->dX;
    dZ=p->dZ;
    W=p->W;
    Winv=p->Winv;
    G=p->G;
    Ginv=p->Ginv;
    D=p->D;
    Weig=p->Weig;
    Wvec=p->Wvec;
    last_rhs_mu=p->last_rhs_mu;
    mu=p->mu;
    old_mu=p->old_mu;
    last_alpha=p->last_alpha;
    oldX=p->oldX;
    oldZ=p->oldZ;
    lamX=p->lamX;
    PX=p->PX;
    tmpsym=p->tmpsym;
    tmpsym2=p->tmpsym2;
    tmpvec=p->tmpvec;
    tmpmat=p->tmpmat;

    testmat=p->testmat;
    
    return 0;
  }

  Real PSCIPBundleBlock::evaluate_trace_x()
  {return trace(X);}

  Real PSCIPBundleBlock::evaluate_trace_z()
  {return trace(Z);}

  Real PSCIPBundleBlock::evaluate_trace_dx()
  {return trace(dX);}


    /// C=beta*C+alpha*B*A where B and A may be transposed; carry out the model part of this beginning at startindex_model and beta for the part, that is added to (the calling routine has to make sure beta is not executed repeatedly if the same part is affected by other models as well)
Matrix& PSCIPBundleBlock::B_times(const Matrix& A,
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
  if (alpha==0.)
    return C;

  if (
      (oracle_data)
      &&(Btrans)
      &&(C.rowdim()==oracle_data->get_opAt().coldim())
      &&(oracle_data->get_opAt().nzcoldim()>2*oracle_data->get_bundlevecs().rowdim())
      ){
    if (Atrans==0){
      assert(A.coldim()==C.coldim());
      Symmatrix tmpsym;
      Matrix eigvec,eigval;
      Matrix bvec;
      for(Integer j=0;j<A.coldim();j++){
	tmpsym.init_svec(oracle_data->get_bundlevecs().coldim(),A.get_store()+j*A.rowdim()+startindex_model);
	if (tmpsym.eig(eigvec,eigval)){
	  if (cb_out())
	    get_out()<<"**** WARNING PSCIPBundleBlock::B_times(): tmpsym.eig() failed"<<std::endl;
	}
	genmult(oracle_data->get_bundlevecs(),eigvec,bvec);
	eigvec.init(oracle_data->get_opAt().coldim(),1,0.);
	if (oracle_data->get_opAt().Gram_ip(eigvec,bvec,&eigval)){
	  if (cb_out())
	    get_out()<<"**** WARNING PSCIPBundleBlock::B_times(): Gram_ip() failed"<<std::endl;
	}
	assert(eigvec.rowdim()==C.rowdim());
	mat_xpeya(C.rowdim(),C.get_store()+j*C.rowdim(),eigvec.get_store(),alpha);
      }
    }
    else {
      assert(A.rowdim()==C.coldim());
      Symmatrix tmpsym;
      Matrix eigvec,eigval;
      Matrix bvec;
      for(Integer j=0;j<A.rowdim();j++){
	tmpsym.init_svec(oracle_data->get_bundlevecs().coldim(),A.get_store()+startindex_model*A.rowdim()+j,A.rowdim());
	if (tmpsym.eig(eigvec,eigval)){
	  if (cb_out())
	    get_out()<<"**** WARNING PSCIPBundleBlock::B_times(): tmpsym.eig() failed"<<std::endl;
	}
	genmult(oracle_data->get_bundlevecs(),eigvec,bvec);
	eigvec.init(oracle_data->get_opAt().coldim(),1,0.);
	if (oracle_data->get_opAt().Gram_ip(eigvec,bvec,&eigval)){
	  if (cb_out())
	    get_out()<<"**** WARNING PSCIPBundleBlock::B_times(): Gram_ip() failed"<<std::endl;
	}
	assert(eigvec.rowdim()==C.rowdim());
	mat_xpeya(C.rowdim(),C.get_store()+j*C.rowdim(),eigvec.get_store(),alpha);
      }
    }
    // // BEGIN TEST
    // Integer simod=startindex_model;
    // for (Integer i=0;i<vecdim;i++){
    //   globalbundle[unsigned(startindex_bundle+i)].left_genmult(A,testmat,-alpha,1.,(Btrans==0),Atrans,simod++);
    // }
    // if (norm2(testmat)<1e-10*(1.+norm2(A))){
    //   std::cout<<" Gok"<<std::flush;
    // }
    // else {
    //   std::cout<<" Gfail"<<std::flush;
    //   std::cout<<" diff="<<testmat;
    // }
    // // END TEST
    return C;
  }
  
  if ((vecdim==bundle_dim)&&(B.rowdim()!=vecdim)) {
    for (Integer i=0;i<vecdim;i++){
      globalbundle[unsigned(startindex_bundle+i)].left_genmult(A,C,alpha,1.,(Btrans==0),Atrans,startindex_model++);
    }
  }
  else {  
    assert(B.rowdim()==vecdim);
    if (Btrans==0){
      if (Atrans==0){
	assert(A.rowdim()==B.coldim());
	assert(A.coldim()==C.coldim());
	const Real*ap=A.get_store();
	for (Integer Acol=0;Acol<A.coldim();Acol++,ap+=A.rowdim()){
	  Real *cp=C.get_store()+Acol*C.rowdim()+startindex_model;
	  const Real* const cpend=cp+vecdim;
	  const Integer* indp=Bnzcolind.get_store();
	  for (Integer i=Bnzcolind.rowdim();--i>=0;indp++){
	    const Real aa=alpha*(*(ap+*indp));
	    if (aa!=0.){
	      const Real* bp=B.get_store()+(*indp)*vecdim;
	      for(;cp!=cpend;){
		(*cp++)+=aa*(*bp++);
	      }
	      cp-=vecdim;
	    }
	  }
	}
      } //endif ATRANS==0
      else {
	genmult(B,A,tmpmat,alpha,0.,0,Atrans);
	for (Integer i=0;i<tmpmat.coldim();i++){
	  mat_xpey(tmpmat.rowdim(),C.get_store()+startindex_model+i*C.rowdim(),tmpmat.get_store()+i*tmpmat.rowdim());
	}
      }
    }
    else { //Btrans=1, startindex_model gives starting row(column) of B
      assert(B.coldim()==C.rowdim());
      if (Atrans==0){
	assert(A.coldim()==C.coldim());
	assert(startindex_model+B.rowdim()<=A.rowdim());
	Real* cp=C.get_store();
	for(Integer Acol=0;Acol<A.coldim();Acol++){
	  const Real* ap=A.get_store()+Acol*A.rowdim()+startindex_model;
	  const Real* const apend=ap+vecdim;
	  const Integer* indp=Bnzcolind.get_store();
	  for (Integer i=Bnzcolind.rowdim();--i>=0;indp++){
	    const Real *bp=B.get_store()+(*indp)*vecdim;
	    Real val=0.;
	    for(;ap<apend;){
	      val+=(*ap++)*(*bp++);
	    }
	    ap-=vecdim;
	    (*(cp+*indp))+=alpha*val;
	  }
	  cp+=C.rowdim();
	}
      } //end if Atrans==0
      else { //Atrans==1
	tmpmat.newsize(vecdim,C.coldim()); chk_set_init(tmpmat,1);
	// if (Atrans==0){
	//   assert(A.coldim()==C.coldim());
	//   assert(startindex_model+vecdim<=A.rowdim());
	//   for(Integer i=0;i<tmpmat.coldim();i++){
	//     mat_xey(vecdim,tmpmat.get_store()+i*vecdim,A.get_store()+startindex_model+i*A.rowdim());
	//   }
	// }
	// else {
	assert(A.rowdim()==C.coldim());
	assert(startindex_model+vecdim<=A.coldim());
	for(Integer i=0;i<vecdim;i++){
	  mat_xey(A.rowdim(),tmpmat.get_store()+i,vecdim,A.get_store()+(startindex_model+i)*A.rowdim(),1);
	}
	//} 
	genmult(B,tmpmat,C,alpha,1.,1,0);
      } //end else Atrans==1
    } //end else Btrans=1
  }
    
  return C;
}
						   
  /// C=beta*C+alpha*A*B where B and B may be transposed; carry out the model part of this beginning at startindex_model 
Matrix& PSCIPBundleBlock::times_B(const Matrix& A,
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
  if (alpha==0.)
    return C;
  
  if ((vecdim==bundle_dim)&&(B.coldim()!=vecdim)){
    for (Integer i=0;i<vecdim;i++){
      globalbundle[unsigned(startindex_bundle+i)].right_genmult(A,C,alpha,1.,Atrans,(Btrans==0),startindex_model++);
    }
  }
  else {
    assert(B.rowdim()==vecdim);
    if (Btrans!=0){
      genmult(A,B,tmpmat,alpha,0.,Atrans,1);
      for (Integer i=0;i<vecdim;i++){
	mat_xpey(C.rowdim(),C.get_store()+(startindex_model+i)*C.rowdim(),1,tmpmat.get_store()+i*tmpmat.rowdim(),1);
      }
    }
    else { //Btrans==0, startindex_model refers to A
      assert(B.coldim()==C.coldim());
      tmpmat.newsize(vecdim,C.rowdim()); chk_set_init(tmpmat,1);
      if (Atrans==0){
	assert(A.rowdim()==C.rowdim());
	assert(A.coldim()>=startindex_model+vecdim);	
	for(Integer i=0;i<vecdim;i++){
	  mat_xey(A.rowdim(),tmpmat.get_store()+i,vecdim,A.get_store()+(startindex_model+i)*A.rowdim(),1);
	}
      }
      else { //Atrans==1
	assert(A.coldim()==C.rowdim());
	assert(A.rowdim()>=startindex_model+vecdim);	
	for(Integer i=0;i<A.coldim();i++){
	  mat_xey(vecdim,tmpmat.get_store()+i*vecdim,A.get_store()+startindex_model+i*A.rowdim());
	}
      }	
      genmult(tmpmat,B,C,alpha,1.,1,0);
    }
  }
  return C;
}

  ///add the main diagonal block tranpose(projection)*diagvec*projection to bigS starting at startindex
Symmatrix& PSCIPBundleBlock::add_BDBt(const Matrix& diagvec,
				    Symmatrix& bigS,
				    bool minus,
				    Integer startindex,
				    Matrix& globBt,
				    Integer startindex_model,
				    MinorantBundle& globalbundle,
				    Integer startindex_bundle)
{
  if (vecdim==bundle_dim) {
    if (minus) {
      for (Integer i=0;i<vecdim;i++){
	MinorantPointer& mp=globalbundle[unsigned(startindex_bundle+i)];
	for(Integer j=i+startindex_model;j<globBt.coldim();j++){
	  bigS(startindex+startindex_model+i,startindex+j)-=mp.ip(globBt,&diagvec,j*globBt.rowdim());
	}
      }
    }
    else {
      for (Integer i=0;i<vecdim;i++){
	MinorantPointer& mp=globalbundle[unsigned(startindex_bundle+i)];
	for(Integer j=i+startindex_model;j<globBt.coldim();j++){
	bigS(startindex+startindex_model+i,startindex+j)+=mp.ip(globBt,&diagvec,j*globBt.rowdim());
	}
      }
    }
  }
  else {
    assert(B.rowdim()==vecdim);
    assert(B.coldim()==globBt.rowdim());
    assert(diagvec.rowdim()==globBt.rowdim());
    assert(bigS.rowdim()>=startindex+globBt.coldim());
    Integer rd=globBt.rowdim();
    const Real* dp=diagvec.get_store();
    const Real* vip=globBt.get_store()+startindex_model*rd;
    if (minus) {
      for (Integer i=0;i<vecdim;i++,vip+=rd){
	const Real* vjp=vip;
	for (Integer j=i+startindex_model;j<globBt.coldim();j++,vjp+=rd){
	  bigS(startindex+startindex_model+i,startindex+startindex_model+j)-=mat_ip(B.coldim(),vip,vjp,dp);
	}
      }
    }
    else {
      for (Integer i=0;i<vecdim;i++,vip+=rd){
	const Real* vjp=vip;
	for (Integer j=i+startindex_model;j<globBt.coldim();j++,vjp+=rd){
	  bigS(startindex+startindex_model+i,startindex+startindex_model+j)+=mat_ip(B.coldim(),vip,vjp,dp);
	}
      }
    }
  }
  return bigS;
}

  /// get the current matrix for the coupling matrix Bt in the first row of blocks
Matrix& PSCIPBundleBlock::get_Bt(Matrix& globBt,
			       Integer startindex_model,
			       MinorantBundle& globalbundle,
			       Integer startindex_bundle)
{
  assert(&B!=&globBt);
  
  for (Integer i=0;i<vecdim;i++){
    Real dummy;
    globalbundle[unsigned(startindex_bundle+i)].get_minorant(dummy,globBt,startindex_model+i);
    //std::cout<<" globBt="<<globBt;
  }
  return globBt;
}
  
  /// adds opB transposed times modelx (without constant affine term) to the arguments
int PSCIPBundleBlock::add_modelx_aggregate(Real& val,
				   Matrix& vec,
				   MinorantBundle& globalbundle,
				   Integer startindex_bundle)
{
  svec(X,tmpvec);
  //std::cout<<" dv_aggrv="<<tmpvec;
  for (Integer i=0;i<vecdim;i++){
    globalbundle[unsigned(startindex_bundle+i)].get_minorant(val,vec,0,tmpvec(i),true);
  }
  return 0;
}

  /// set the model violation for the current system solution 
int PSCIPBundleBlock::get_sysviol_model(Matrix& sysviol_model,
				      Integer startindex_model,
				      const Matrix& dy,
				      const Real deltatrdual,
				      MinorantBundle& globalbundle,
				      Integer startindex_bundle)
{
  for(Integer i=0;i<vecdim;i++){
    sysviol_model(startindex_model+i)=globalbundle[unsigned(startindex_bundle+i)].evaluate(-1,dy,false)-diff_model(i);
  }
  svec(Z,sysviol_model,1.,true,startindex_model);
  svec(dZ,sysviol_model,1.,true,startindex_model);
  Integer cnt=startindex_model;
  for (Integer i=rowdim;i>0;--i){
    sysviol_model(cnt)-=deltatrdual;
    cnt+=i;
  }
  return 0;
}

  /// set z to the slack of the bundle and return a value>=0 that needs to be added to make it feasible
int PSCIPBundleBlock::set_bundle_z(const Matrix& y,
		    MinorantBundle& globalbundle,
		    Integer startindex_bundle,
		    Real& add_center_value)
{
  diff_model.newsize(vecdim,1);chk_set_init(diff_model,1); 
  for(Integer i=0;i<vecdim;i++){
    diff_model(i)=-globalbundle[unsigned(startindex_bundle+i)].evaluate(-1,y);
  }
  return PSCIPBlock::set_z(diff_model,0,add_center_value);
}

int PSCIPBundleBlock::add_trace(Matrix& vec,Real trace_dual,Integer startindex) 
{
  Integer cnt=startindex;
  for (Integer i=rowdim;i>0;--i){
    vec(cnt)+=trace_dual;
    cnt+=i;
  }
  return 0;
}

  /// set the trace premultiplied by sqrt(inv(xiz)) in vec[startindex+0,...,startindex+dim_bundle()-1]
  int PSCIPBundleBlock::set_xizinvsqrt_trace(Matrix& vec,
					   Integer startindex)
  {
    compute_NTscaling();
    //compute_Weig_Wvec();

    rankadd(G,tmpsym);
    tmpsym.store_svec(vec.get_store()+startindex);

    /*
     // W=G'*G=(skron(G',G')*svec(I))=Wvec*Wvec'
     // (W skron W) =  (Wvec skron Wvec)*(Wvec' skron Wvec')
     // id*(W skron W)*id =  (Wvec' id Wvec)' * (Wvec' id Wvec) = Diag(Weig)*Diag(Weig)
    for(Integer i=0;i<rowdim;i++){
      vec(startindex++)=Weig(i);
      for(Integer j=i+1;j<rowdim;j++){
	vec(startindex++)=0.;
      }
    }
    */


    // //TEST begin
    // {
    //   Matrix trace_vec(vecdim,1,0.);
    //   add_trace(trace_vec,1.,0);
    //   Matrix systr=trace_vec;
    //   PSCIPBlock::apply_xizinv(systr,0);
    //   Real schurtr=ip(systr,trace_vec);
      
    //   Matrix ivec(vecdim,1,0.);
    //   Matrix sys0(vecdim,vecdim,0.);
    //   for(Integer i=0;i<vecdim;i++){
    // 	ivec.init(vecdim,1,0.);ivec(i)=1.;
    // 	Real zeta_inval=ip(systr,ivec)/schurtr;
    // 	add_trace(ivec,-zeta_inval,0);	
    // 	PSCIPBlock::apply_xizinv(ivec,0);
    // 	for(Integer j=0;j<vecdim;j++)
    // 	  sys0(j,i)=ivec(j);
    //   }
    //   std::cout<<" innersys0="<<sys0;
    //   Symmatrix sysskr;
    //   skron(W,W,sysskr);
    //   rankadd(systr,sysskr,-1./ip(systr,trace_vec),1.);
    //   std::cout<<" syskr="<<sysskr;
    //   if (norm2(sys0-sysskr)>1e-10*trace(sys0)){
    // 	std::cout<<" innerskr differ ="<<sys0-sysskr;
    //   }

    //   Matrix subspace1(Diag(Matrix(vecdim,1,1.)));
    //   Matrix subspace2(Diag(Matrix(vecdim,1,1.)));
    //   Symmatrix tsym2;
    //   Symmatrix tsym;
    //   Symmatrix tsymt;
    //   Matrix gramsys1(vecdim,vecdim,0.);
    //   Matrix gramsys2(vecdim,vecdim,0.);
    //   Matrix sqrttr1(vecdim,1,0.);
    //   Matrix sqrttr2(vecdim,1,0.);
    //   //W.store_svec(sqrttr.get_store());
    //   rankadd(Wvec,tsym,1.,0.,1);
    //   tsym.store_svec(sqrttr1.get_store());
    //   Real n2sqrttr1=norm2(sqrttr1);
    //   rankadd(G,tsym);
    //   tsym.store_svec(sqrttr2.get_store());
    //   Real n2sqrttr2=norm2(sqrttr2);
    //   Real lambda=1;
    //   Matrix tmat;
    //   genmult(sqrttr1,subspace1,tmat,1./n2sqrttr1,0.,1,0);
    //   genmult(sqrttr1,tmat,subspace1,-lambda/n2sqrttr1,1.);
    //   genmult(sqrttr2,subspace2,tmat,1./n2sqrttr2,0.,1,0);
    //   genmult(sqrttr2,tmat,subspace2,-lambda/n2sqrttr2,1.);
	
    //   for(Integer i=0;i<subspace1.coldim();i++){
    // 	tsym2.init_svec(rowdim,subspace1.get_store()+i*subspace1.rowdim());
    // 	symscale(tsym2,Wvec,tsym,1.,0.,1);
    // 	tsym.store_svec(gramsys1.get_store()+i*vecdim);
    // 	tsym2.init_svec(rowdim,subspace2.get_store()+i*subspace2.rowdim());
    // 	symscale(tsym2,G,tsymt);
    // 	tsymt.store_svec(gramsys2.get_store()+i*vecdim);
    //   }
    //   Matrix sys1=gramsys1*transpose(gramsys1);
    //   Matrix sys2=gramsys2*transpose(gramsys2);
    //   std::cout<<" innersys1="<<sys1;
    //   if (norm2(sys0-sys1)>1e-10*trace(sys0)){
    // 	std::cout<<" innersys1 differ ="<<sys0-sys1;
    //   }
    //   std::cout<<std::endl;
    //   std::cout<<" innersys2="<<sys2;
    //   if (norm2(sys0-sys2)>1e-10*trace(sys0)){
    // 	std::cout<<" innersys2 differ ="<<sys0-sys2;
    //   }
    //   std::cout<<std::endl;
		      
    // }
    // //TEST end

    
    return 0;
  }
  
  ///add trace_dual*trace to diff_model for the right hand side (negative of current model violation)
int PSCIPBundleBlock::add_trace_to_diff_model(Real trace_dual)
{
  Integer cnt=0;
  for (Integer i=rowdim;i>0;--i){
    diff_model(cnt)+=trace_dual;
    cnt+=i;
  }
  return 0;
}

///return the squared Euclidean norm of the dual model violation  
Real PSCIPBundleBlock::dualviol_2normsqr()
{
  tmpvec.init(diff_model);
  svec(Z,tmpvec,-1.,true);
  return ip(tmpvec,tmpvec);
}

  /// move to (x+alpha*dx, z+alpha*dz), update diff_model and possibly reduce the model size if some part is below zero_threshold
int PSCIPBundleBlock::do_bundle_step(Real alpha,
			    const Matrix& y,
			    MinorantBundle& globalbundle,
			    Integer startindex_bundle,
			    Real trace_dual,
			    Real /* trace_rhs */)
{
  PSCIPBlock::do_step(alpha);

  int eigstat=X.eig(PX,lamX,false);
  if (eigstat){
    if (cb_out())
      get_out()<<"**** WARNING PSCIPBlock::do_step(....): eigenvalue factorization failed and returned "<<eigstat<<std::endl;
  }

  diff_model.newsize(vecdim,1);chk_set_init(diff_model,1); 
  for(Integer i=0;i<vecdim;i++){
    diff_model(i)=-globalbundle[unsigned(startindex_bundle+i)].evaluate(-1,y);
  }
  add_trace_to_diff_model(trace_dual);
  
  return 0;
}

  /// If mu is not zero, always add the centering term for this mu as well;
int PSCIPBundleBlock::set_modelrhs(Matrix& globalrhs, 
				 Real rhsmu,
				 Real rhscorr,
				 Integer startindex_model)
{
  mat_xey(vecdim,globalrhs.get_store()+startindex_model,diff_model.get_store());
  return add_muxinv(globalrhs,startindex_model,rhsmu,rhscorr,true);
}

  ///add the "scaled" minorant outer products to globalsys, where the correct minroants start at the given index
int PSCIPBundleBlock::add_BtinvsysB(Symmatrix& globalsys,
			    const MinorantBundle& globalbundle,
			    Integer startindex_bundle)
{
  form_B(globalsys.rowdim(),globalbundle,startindex_bundle);
  //return add_AxizinvAt(B,globalsys,false,true);
  compute_NTscaling();

  Integer ncols=B.coldim();
  tmpmat.init(ncols,vecdim,0.);

  for(Integer i=0;i<Bnzcolind.rowdim();i++){
    Integer ind=Bnzcolind(i);
    tmpsym.init_svec(rowdim,B.get_store()+ind*vecdim);
    symscale(tmpsym,G,tmpsym2,1.,0.,1);
    tmpsym2.store_svec(tmpmat.get_store()+ind,ncols);
  }
      
  rankadd(tmpmat,globalsys,1.,1.);
  return 0;
}


   /// glob_lowrank(:,...)=(Btsyssqrtinv), trafotrace(...)=syssqrtinv*trace
  int PSCIPBundleBlock::Schur_transform_bundle(Matrix& glob_lowrank,
					     MinorantBundle& globalbundle,
					     Integer startindex_bundle,
					     Matrix& trafotrace,
					     Integer startindex_trace)
  {
   assert(trafotrace.rowdim()>=startindex_trace+vecdim);
   assert(glob_lowrank.coldim()>=startindex_bundle+vecdim);

   if (vecdim==0)
     return 0;
   
   compute_NTscaling();

   Integer rdim=glob_lowrank.rowdim();
   form_B(rdim,globalbundle,startindex_bundle);

   /// Integer method=((norm2(W-Diag(diag(W)))/trace(W))<0.001);
   Real *wp=W.get_store();
   Real offdnorm=0.;
   Real trW=0.;
   for (Integer i=0;i<tmpsym.rowdim();i++){
     trW+=(*wp++);
     for (Integer j=1;j<tmpsym.rowdim();j++){
       offdnorm+=sqr(*wp++);
     }					       
   }

   //Integer method=(offdnorm<0.01*trW);
   Integer method=0;
    
   
   if (cb_out(4)){
     get_out()<<" PSCIPlrmethod="<<method<<std::endl;
   }
   
   switch(method){
   case 0: {
     compute_Weig_Wvec();
       
     // W=G'*G=(skron(G',G')*svec(I))=Wvec*Wvec'
     // (W skron W) =  (Wvec skron Wvec)*(Wvec' skron Wvec')
     // id*(W skron W)*id =  (Wvec' id Wvec)' * (Wvec' id Wvec) = Diag(Weig)*Diag(Weig)

     Integer ind=startindex_trace;
     for (Integer i=0;i<rowdim;i++){
       trafotrace(ind++)=Weig(i);   // !!!  -id (negative trace!)
       for(Integer j=i+1;j<rowdim;j++){
	 trafotrace(ind++)=0.;
       }
     }
     
     Real *glrp=glob_lowrank.get_store()+startindex_bundle*rdim;
     for(Integer i=0;i<rdim;i++){
       /*
       tmpvec.init(vecdim,1,B.get_store()+i*vecdim);
       sveci(tmpvec,tmpsym2);
       symscale(tmpsym2,Wvec,tmpsym,1.,0.);
       svec(tmpsym,tmpvec);
       mat_xey(vecdim,glrp+i,rdim,tmpvec.get_store(),1);
       */
       tmpsym2.init_svec(rowdim,B.get_store()+i*vecdim);
       symscale(tmpsym2,Wvec,tmpsym,1.,0.);
       tmpsym.store_svec(glrp+i,rdim);
     }
     break;
   }

   case 1:{
     //use only the diagonal of W
     Matrix diagw(diag(W));
     diagw.sqrt();
     
     Integer ind=startindex_trace;
     for (Integer i=0;i<rowdim;i++){
       Real wisqrt=diagw(i);
       trafotrace(ind++)=wisqrt*wisqrt;   // !!!  -id (negative trace!)
       for(Integer j=i+1;j<rowdim;j++){
	 trafotrace(ind++)=wisqrt*diagw(j);
       }
     }
     
     
     ind=0;
     for (Integer i=0;i<rowdim;i++){
       for (Integer j=i;j<rowdim;j++){
	 Real wij=trafotrace(startindex_trace+ind);
	 mat_xeya(glob_lowrank.rowdim(),
		  glob_lowrank.get_store()+(startindex_bundle+ind)*rdim,1,
		  B.get_store()+ind,vecdim,wij);
	 ind++;
       }
     }
     break;
   }
   }
   
   return 0;
  }


  int PSCIPBundleBlock::add_pcsubspace(Matrix& lowrank,
				       Matrix& sigma_guess,
				       const Matrix& Diag_inv,
				       Real minval,
				       Real diaginvval,
				       Matrix & minus_trmult,
				       Real schur_trace,
				       MinorantBundle& globalbundle,
				       Integer startindex_bundle)
  {
    assert(minval>0.);
    compute_NTscaling();
    compute_Weig_Wvec();

    if (Weig.rowdim()==0)
      return 0;

    Real nrmbnd=0.;
    Matrix localsqrBnorms;

    if (diaginvval>0.){
      if (sqrBnorms.rowdim()==0){
	sqrBnorms.init(vecdim,1,0.);
	for (Integer i=0;i<vecdim;i++){
	  sqrBnorms(i)=globalbundle[unsigned(startindex_bundle+i)].norm_squared();
	}
      }
      localsqrBnorms.init(sqrBnorms,diaginvval);
      nrmbnd=std::sqrt(max(localsqrBnorms));
    }    
    else {
      localsqrBnorms.init(vecdim,1,0.);
      for (Integer i=0;i<vecdim;i++){
	Real d=globalbundle[unsigned(startindex_bundle+i)].norm_squared(&Diag_inv);
	localsqrBnorms(i)=d;
	if (d>nrmbnd)
	  nrmbnd=d;
      }
      nrmbnd=std::sqrt(nrmbnd);
    }
          
    if (max(Weig)*nrmbnd<minval)
      return 0;
    
    tmpsym.init(Weig.rowdim(),0.);
    Real dmax=1.;
    Real eigijmax=0.;
    for(Integer i=0;i<Weig.rowdim();i++){
      const Real eigi=Weig(i);
      Real d=eigi*eigi*(1.-eigi*eigi/schur_trace);
      tmpsym(i,i)=d;
      if (d>dmax)
	dmax=d;
      for(Integer j=i+1;j<Weig.rowdim();j++){
	const Real eigij=eigi*Weig(j);
	if (eigij>eigijmax)
	  eigijmax=eigij;
	tmpsym(i,j)=-eigij*eigij/schur_trace;
      }
    }
    tmpsym/=dmax;
    Matrix tmpeig;
    Matrix eigvec;
    tmpsym.eig(eigvec,tmpeig,false);
    tmpeig*=dmax;
    eigijmax=std::sqrt(max(eigijmax,tmpeig(0)));

    //std::cout<<" Weig="<<transpose(Weig)<<" tmpeig="<<transpose(tmpeig);

    //count candidates and reserve space
    Real threshold=minval*minval/nrmbnd/nrmbnd;
    Indexmatrix candpairs(vecdim,2);
    candpairs.init(0,2,0);
    Indexmatrix candcnt(rowdim,1,Integer(0));
    Indexmatrix delind(rowdim,1);
    delind.init(0,1,Integer(0));
    for(Integer i=0;i<Weig.rowdim();i++){
      const Real eigi=Weig(i);
      if(tmpeig(i)<=threshold){
	delind.concat_below(i);
      }
      else {
	tmpvec.init(eigvec.col(i));
	tmpvec/=Weig;
	scaledrankadd(Wvec,tmpvec,tmpsym);
	tmpvec=svec(tmpsym);
	assert(std::fabs(norm2(tmpvec)-1.)<1e-8);
	assert(vecdim==tmpvec.rowdim());
	Real d=0;
	for(Integer k=0;k<vecdim;k++)
	  d+=localsqrBnorms(k)*sqr(tmpvec(k));
	d=std::sqrt(d*tmpeig(i));
	if (d>=minval){
	  if (cb_out(2)){
	    get_out().precision(4);
	    get_out()<<" P{"<<i<<";"<<d<<"}";
	  }
	  candcnt(i)+=1;
	  candpairs.enlarge_below(1,Integer(i));
	}
	else {
	  delind.concat_below(i);
	}
      }
      Integer j=i+1;
      while ((j<Weig.rowdim())&&(eigi*Weig(j)>threshold)){
	const Real eigij=eigi*Weig(j);
	rank2add(Wvec.col(i),Wvec.col(j),tmpsym,std::sqrt(2./eigij));
	tmpvec=svec(tmpsym);
	assert(std::fabs(norm2(tmpvec)-1.)<1e-8);
	assert(vecdim==tmpvec.rowdim());
	Real d=0;
	for(Integer k=0;k<vecdim;k++)
	  d+=localsqrBnorms(k)*sqr(tmpvec(k));
	d=std::sqrt(d*eigij);
	if (d>=minval){
	  if (cb_out(2)){
	    get_out().precision(4);
	    get_out()<<" P{"<<i<<","<<j<<";"<<d<<"}";
	  }
	  candcnt(i)++;
	  candpairs.enlarge_below(1,Integer(i));
	  candpairs(candpairs.rowdim()-1,1)=j;
	}
	j++;
      }
    }
    Integer cnt=candpairs.rowdim();

    if (cnt==0)
      return 0;

    form_B(lowrank.rowdim(),globalbundle,startindex_bundle);
    
    Integer nextcol=lowrank.coldim();
    lowrank.enlarge_right(cnt,0.);
    sigma_guess.enlarge_below(cnt,0.);
    
    //first append the diagonal elements that are large enough
    if (delind.rowdim()>0){
      tmpeig.delete_rows(delind);
      eigvec.delete_cols(delind);
    }
    Real Weignrmsqr=ip(Weig,Weig);
    assert(Weignrmsqr<(1.+100.*eps_Real)*schur_trace);
    Real Weigfactor=(1.-std::sqrt(max(0.,1-Weignrmsqr/schur_trace)))/Weignrmsqr;
    for(Integer i=0;i<tmpeig.rowdim();i++){
      //transform vector
      tmpvec.init(eigvec.col(i),1./std::sqrt(tmpeig(i)));
      tmpvec%=Weig;
      tmpvec.xpeya(Weig,-Weigfactor*ip(Weig,tmpvec));

      Real trmult=ip(Weig,tmpvec);

      //mutliply with V_W \skron V_W
      scaledrankadd(Wvec,tmpvec,tmpsym);
      tmpvec=svec(tmpsym);

      //multiply with B
      //genmult(B,tmpvec,tmpmat,1.,0.,1);
      Real *lrp=lowrank.get_store()+nextcol*lowrank.rowdim();
      for(Integer ci=0;ci<Bnzcolind.rowdim();ci++){
	Integer cind=Bnzcolind(ci);
	Real d=mat_ip(vecdim,B.get_store()+cind*B.rowdim(),tmpvec.get_store());
	*(lrp+cind)=d;
      }
      if (cb_out(2)){
      	get_out()<<" Pd("<<i<<")";
      }
      sigma_guess(nextcol)=0.;
      minus_trmult.concat_below(trmult);
      nextcol++;
    }  
    
    //now append the mixed eigenvalues if large enough
    Matrix sgvec(lowrank.rowdim(),1);
    Integer tmpmatind=-1;
    for(Integer i=0;i<candpairs.rowdim();i++){
      Integer ind0=candpairs(i,0);
      Integer ind1=candpairs(i,1);
      if (ind0==ind1)
	continue;
      const Real eigi=Weig(ind0);
      const Real eigij=eigi*Weig(ind1);
      if ((candcnt(ind0)>3)&&(ind0!=tmpmatind)){
	sveciB_times_vec(Wvec.col(ind0),tmpmat);
	tmpmatind=ind0;
      }
      Real sg=0.;
      Real *lrp=lowrank.get_store()+nextcol*lowrank.rowdim();
      if (tmpmatind==ind0){
	if (ind0==ind1){
	  //genmult(tmpmat,Wvec.col(ind0),sgvec,1.,0.,1);
	  for(Integer ci=0;ci<Bnzcolind.rowdim();ci++){
	    Integer cind=Bnzcolind(ci);
	    Real d=mat_ip(rowdim,tmpmat.get_store()+cind*tmpmat.rowdim(),Wvec.get_store()+ind1*Wvec.rowdim());
	    *(lrp+cind)=d;
	    if (diaginvval<=0.){
	      sg+=d*d*Diag_inv(cind);
	    }
	  }
	}
	else {
	  //genmult(tmpmat,Wvec.col(ind1),sgvec,std::sqrt(2.),0.,1);
	  for(Integer ci=0;ci<Bnzcolind.rowdim();ci++){
	    Integer cind=Bnzcolind(ci);
	    Real d=std::sqrt(2.)*mat_ip(rowdim,tmpmat.get_store()+cind*tmpmat.rowdim(),Wvec.get_store()+ind1*Wvec.rowdim());
	    *(lrp+cind)=d;
	    if (diaginvval<=0.){
	      sg+=d*d*Diag_inv(cind);
	    }
	    else
	      sg+=d*d;
	  }
	}
      }
      else {
	if (ind0==ind1)
	  rankadd(Wvec.col(ind0),tmpsym,1./eigi);
	else
	  rank2add(Wvec.col(ind0),Wvec.col(ind1),tmpsym,std::sqrt(2./eigij));
	tmpvec=svec(tmpsym);
	assert(std::fabs(norm2(tmpvec)-1.)<1e-10*eigij);
	const Real ld=std::sqrt(eigij);	    
	//genmult(B,tmpvec,sgvec,ld,0.,1);
	for(Integer ci=0;ci<Bnzcolind.rowdim();ci++){
	  Integer cind=Bnzcolind(ci);
	  Real d=ld*mat_ip(vecdim,B.get_store()+cind*B.rowdim(),tmpvec.get_store());
	  *(lrp+cind)=d;
	  if (diaginvval<=0.){
	    sg+=d*d*Diag_inv(cind);
	  }
	  else
	    sg+=d*d;
	}
      }
      if (diaginvval>0.)
	sg=std::sqrt(sg*diaginvval);
      else
	sg=std::sqrt(sg);
      if (cb_out(2)){
      	get_out()<<" Pm("<<ind0<<","<<ind1<<";"<<sg<<","<<(sg>=minval)<<")";
      }
      if (sg>=minval){
	sigma_guess(nextcol)=sg;
	nextcol++;
	if (ind1==ind0)
	  minus_trmult.concat_below(eigi);
	else 
	  minus_trmult.concat_below(0.);
      }
    }
    
    if (nextcol<lowrank.coldim()){
      sigma_guess.reduce_length(nextcol);
      lowrank.delete_cols(Range(nextcol,lowrank.coldim()-1));
    }

    return 0;
  }



  /** @brief add diag(Bt*sqrt(invsys)*(I-lambda*trvec*trvec')*sqrt(invsys)*B) to diagonal 
   */
  int PSCIPBundleBlock::add_bundle_xizinv_diagonal(Matrix& diagonal,
						   Matrix& ipBtrvec,
						   MinorantBundle& globalbundle,
						   Integer startindex_bundle,
						   const Matrix& trafotrace,
						   Integer startindex_trace)
  {
    assert(vecdim==dim_bundle());
    
    compute_NTscaling();
    //compute_Weig_Wvec();

    form_B(diagonal.rowdim(),globalbundle,startindex_bundle);
    
    assert(trafotrace.rowdim()>=startindex_trace+vecdim);
    for(Integer i=0;i<diagonal.rowdim();i++){
      /*
      tmpvec=B.col(i);
      sveci(tmpvec,tmpsym2);
      */
      tmpsym2.init_svec(rowdim,B.get_store()+i*vecdim);
      symscale(tmpsym2,G,tmpsym,1.,0.,1);
      svec(tmpsym,tmpvec);
      diagonal(i)+=ip(tmpvec,tmpvec);
      if (ipBtrvec.rowdim()>0)
	ipBtrvec(i)+=mat_ip(vecdim,trafotrace.get_store()+startindex_trace,tmpvec.get_store());
    }
    
    return 0;
  }
					 
  /// add bundle*sqrt(inv(xiz))*subspace to glob_lowrank with bundle(:,si_bundle+1:si_bundle+dim_bundle()-1) and subspace(si_subsp:si_subsp+dim_bundle,:); sqrt(inv(xiz)) has to match that used in set_xizinvsqrt_trace()
  int PSCIPBundleBlock::add_bundle_xizinvsqrt_projection(Matrix& glob_lowrank,
						       Matrix& subspace,
						       Integer startindex_subspace,
						       MinorantBundle& globalbundle,
						       Integer startindex_bundle)
  {
    assert(vecdim==dim_bundle());

    compute_NTscaling();

    form_B(glob_lowrank.rowdim(),globalbundle,startindex_bundle);
    
    if (startindex_subspace>=0){
      tmpmat.newsize(vecdim,subspace.coldim()); chk_set_init(tmpmat,1);
      //tmpmat.init(vecdim,0,0.);
      for(Integer i=0;i<subspace.coldim();i++){
       tmpsym2.init_svec(rowdim,subspace.get_store()+startindex_subspace+i*subspace.rowdim());
       symscale(tmpsym2,G,tmpsym);
       tmpsym.store_svec(tmpmat.get_store()+i*vecdim);
      }
      
      // //begin TEST
      // std::cout<<" trafosubsp="<<tmpmat;
      // //end TEST
      
      B_times(tmpmat,glob_lowrank,1.,1.,1,0,0,globalbundle,startindex_bundle);
      
      // //begin TEST
      // std::cout<<" globlowrank="<<glob_lowrank;
      // //end TEST
    }
    else {
      tmpmat.init(vecdim,glob_lowrank.coldim(),0.);
      B_times(glob_lowrank,tmpmat,1.,0.,0,0,0,globalbundle,startindex_bundle);
      for(Integer i=0;i<tmpmat.coldim();i++){
	tmpsym2.init_svec(rowdim,tmpmat.get_store()+i*tmpmat.rowdim());
	symscale(tmpsym2,G,tmpsym,1.,0.,1);
	tmpsym.store_svec(tmpmat.get_store()+i*tmpmat.rowdim());
      }     
      subspace.concat_right(tmpmat,1);
    }
    
    return 0;
  }




  /// out_vec+=BtinvsysB*in_vec
  int PSCIPBundleBlock::add_BtinvsysB_times(const Matrix& in_vec,
					  Matrix& out_vec,
					  Real zeta_inval,
					  Real *zeta_outval,
					  MinorantBundle& globalbundle,
					  Integer startindex_bundle)
  {
    //Matrix testout(out_vec);
    form_B(in_vec.rowdim(),globalbundle,startindex_bundle);

    tmpvec.init(vecdim,1,0.);
    B_times(in_vec,tmpvec,1.,0.,0,0,0,globalbundle,startindex_bundle);
    add_trace(tmpvec,-zeta_inval,0);
    PSCIPBlock::apply_xizinv(tmpvec,0);
    B_times(tmpvec,out_vec,1.,1.,1,0,0,globalbundle,startindex_bundle);
    if (zeta_outval){
      Integer cnt=0;
      for(Integer i=0;i<rowdim;i++){
	(*zeta_outval)-=tmpvec(cnt);
	cnt+=rowdim-i;
      }
    }

    // std::cout<<" PSCIPdiff="<<norm2(out_vec-testout-transpose(testmat)*testmat*in_vec)<<std::endl;
    
    return 0;
  }

  /// compute dx (and then dz) given step_y and step_trdual on basis of the last rhs computed for the model block
  int PSCIPBundleBlock::set_dx_xizsolvestep(const Matrix& step_y,
				 const Real step_trdual,
				 MinorantBundle& globalbundle,
				 Integer startindex_bundle)
  {
    tmpvec.init(diff_model);
    tmpvec(0)+=step_trdual;
    B_times(step_y,tmpvec,-1.,1.,0,0,0,globalbundle,startindex_bundle);

    //use of tmpvec is ok, because not employed in LinBlock::set_dx_xizsolverhs;
    return PSCIPBlock::set_dx_xizsolverhs(tmpvec,0);
  }
				  

  int PSCIPBundleBlock::get_pscx(Matrix& pscx_eigs,
				 Matrix& pscx_vecs,
				 Real& pscx_growthrate,
				 Matrix& pscx_primalgrowth,
				 Matrix& pscx_dualgrowth) 
  {
    if (lamX.rowdim()!=X.rowdim()){
      int eigstat=X.eig(PX,lamX,false);
      if (eigstat){
	if (cb_out())
	  get_out()<<"**** WARNING PSCIPBlock::get_pscx(....): eigenvalue factorization failed and returned "<<eigstat<<std::endl;
      }
    }
      
    //if (pscx_activity_rank)
    //  *pscx_activity_rank=active_rank(lamX,PX,trace_rhs,cautious);
    get_growth(lamX,PX,pscx_growthrate,pscx_primalgrowth,pscx_dualgrowth);

    pscx_eigs=lamX;
    pscx_vecs=PX;
    
    return 0;  
  }




}
