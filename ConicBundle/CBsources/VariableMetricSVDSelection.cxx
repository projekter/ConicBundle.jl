/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/VariableMetricSVDSelection.cxx
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



#include "VariableMetricSVDSelection.hxx"

using namespace CH_Matrix_Classes;

namespace ConicBundle {


//*****************************************************************************
//                    VariableMetricSVDSelection::vecH_weighted_SVDs
//*****************************************************************************

int VariableMetricSVDSelection::vecH_weighted_SVDs(Matrix& vecH,
				      Real& relprec,
				      const MinorantBundle& bundle,
				      const MinorantBundle& model,
				      const Matrix& modelcoeff,
				      const MinorantPointer& aggregate,
				      Real function_factor,
				      Integer y_id,
				      const Matrix& y,
				      Real weightu,
				      Real model_maxviol)
{
  assert(function_factor>0.);
  assert(model_maxviol>0.);
  assert(weightu>0.);

  if (cb_out(2)){
    get_out()<<"\n wSVDsdscal";
  }

      Integer dim=y.dim();
  Integer msz=Integer(model.size());   // number of minorants in the model (includes aggr)
  Integer bsz=Integer(bundle.size());

  vecH.newsize(dim,max(msz,bsz));
  vecH.init(dim,0,0.);
  Uspace.init(dim,0,0.);

  relprec=min(1e-2,max(relprec,1e-4));

  Real aggrval=aggregate.evaluate(y_id,y)/function_factor;
  Real aggrsqr=aggregate.dual_norm_squared()/function_factor/function_factor;
  Real aggrnorm=std::sqrt(aggrsqr);
  Real coeffzero=max(1e-8*function_factor,1e-3*max(modelcoeff));
  Real curveref=max(weightu,1./std::sqrt(model_maxviol));
  curveref=max(curveref,1000.);

  //---- collect the aggregate and the active model vectors and their data
  Indexmatrix sind;
  sortindex(modelcoeff,sind,false);
  MinorantBundle active(unsigned(msz+1)); active.clear();
  Matrix activecoeff(msz+1,1); activecoeff.init(0,1,0.);
  Matrix activevalue(msz+1,1); activevalue.init(0,1,0.);
  Matrix activenormsqr(msz+1,1); activenormsqr.init(0,1,0.);
  //append the aggregate
  active.push_back(aggregate);
  active[0].scale(1./function_factor);
  activecoeff.concat_below(0.);
  activevalue.concat_below(aggrval);
  activenormsqr.concat_below(aggrsqr);
  //append the active model vectors
  for (Integer i=0;i<msz;i++){
    if (modelcoeff(i)<coeffzero)
      break;
    activecoeff.concat_below(modelcoeff(i));
    const MinorantPointer& act=model[unsigned(i)];
    active.push_back(act);
    activevalue.concat_below(act.evaluate(y_id,y));
    activenormsqr.concat_below(act.dual_norm_squared());
  }
  activecoeff(0)=max(sum(modelcoeff)-sum(activecoeff),1e-3);
  
  //---- assign the minorants to active model elements
  std::vector<std::vector<unsigned> > assignedindices(active.size());
  std::vector<std::vector<Real> > assigneddenoms(active.size());
  for(unsigned i=0;i<bundle.size();i++){
    const MinorantPointer& buni=bundle[i];
    Real vali=buni.evaluate(y_id,y);
    if (vali<aggrval-100.*model_maxviol)
      continue;
    Real normsqri=buni.dual_norm_squared();
    Real minscalval=curveref;  //if above this, ignore it
    Integer minscalind=-1;
    Real mindenom;
    for (Integer j=1;j<activecoeff.rowdim();j++){
      Real ipij=buni.ip(active[unsigned(j)]);
      Real num=std::fabs(.5*(activenormsqr(j)+normsqri)-ipij);
      Real denom=model_maxviol+std::fabs(activevalue(j)-vali);
      if (minscalval*denom>num){
	minscalval=num/denom;
	mindenom=denom;
	minscalind=j;
      }
    }
    if (minscalind>0){
      if (minscalval>1e-8){
	assignedindices[unsigned(minscalind)].push_back(i);
	assigneddenoms[unsigned(minscalind)].push_back(mindenom);
      }
    }
    /*
    else {
      //use the common aggregate as reference for scaling
      Real ipij=buni.ip(active[0]);
      Real num=std::fabs(.5*(aggrsqr+normsqri)-ipij);
      if (num>1e-8*aggrsqr){
	assignedindices[0].push_back(i);
	assigneddenoms[0].push_back(max(model_maxviol+aggrval-vali,num/curveref));
      }
    }
    */
  }
    
  //---- initialize the variables for collecting the scaling information
  Real maxeig=0.;
  Matrix tmpvec(dim,1);
  Matrix tmpmat;
  Matrix modelbasis;
  Symmatrix tmpsym;
  Matrix lambda;
  Matrix lamvecs;


  //----  for each active element generate and append a scaling matrix
  for(Integer i=0;i<activecoeff.rowdim();i++){
    
    std::vector<unsigned>& mdi=assignedindices[unsigned(i)];


    if (mdi.size()==0)
      continue;

    if (cb_out(2))
      get_out()<<" "<<i<<"("<<mdi.size()<<"):";

    // collect the assigned difference vectors
    Real coeff=activecoeff(i);
    MinorantPointer amp=active[unsigned(i)];
    std::vector<Real>& mdd=assigneddenoms[unsigned(i)];
    tmpmat.newsize(dim,Integer(mdi.size())); tmpmat.init(dim,0,0.);
    tmpvec.newsize(dim,1); chk_set_init(tmpvec,1);
    
    for(unsigned j=0;j<mdi.size();j++){
      Real dummy;
      Real factor=std::sqrt(coeff/mdd[j]);
      amp.get_minorant(dummy,tmpvec,0,factor);
      bundle[mdi[j]].get_minorant(dummy,tmpvec,0,-factor,true);
      tmpmat.concat_right(tmpvec);
      if (cb_out(2)){
	get_out()<<" "<<mdi[j]<<"["<<factor<<"]";
	Real ldummy;
	Matrix ltmpvec(dim,1,0.);
	bundle[mdi[j]].get_minorant(ldummy,ltmpvec,0);
	get_out().precision(8);
	get_out()<<ldummy<<", "<<transpose(ltmpvec);
      }
    }
    // std::cout<<"   "<<transpose(tmpmat);

    //compute the SVD and the model basis
    rankadd(tmpmat,tmpsym,1.,0.,1);

    Real scale=1.;
    for(Integer j=1;j<tmpsym.rowdim();j++)
      scale=max(scale,tmpsym(j,j));
    tmpsym/=scale;

    int retval=tmpsym.eig(lamvecs,lambda,false);
    if ((min(lambda)<-1e-6)||(retval)){
      if (cb_out())
	get_out()<<"**** WARNING: VariableMetricSVDSelection::vecH_weighted_SVDs(............): S.eig failed for i="<<i<<" and returned "<<retval<<" (order="<<tmpsym.rowdim()<<") and min eigenvalue="<<min(lambda)<<std::endl;
      
    }

    lambda*=scale;
    if (norm2(lambda)<max(1e-10,1e-6*std::sqrt(activenormsqr(i)))){  
      //"nothing left" after orthogonalization, no contribution here
      modelbasis.init(dim,0,0.);
    }
    else {
      Real lmaxbound=1e-3*lambda(0);
      Integer cnt=1;
      Integer maxdir=lambda.dim();
      while ((cnt<maxdir)&&(lambda(cnt)>lmaxbound))
	cnt++;
      if (cnt<lambda.dim()){
	lambda.reduce_length(cnt);
	lamvecs.delete_cols(Range(cnt,lamvecs.coldim()-1));
      }
      lambda.sqrt();         // now lambda == Sigma
      lambda.inv();
      lamvecs.scale_cols(lambda);  //now modelbasis <- Q*Sigma^{-1} 
   
      genmult(tmpmat,lamvecs,modelbasis);  // vecH <- P = V * Q*Sigma^{-1}
    }

    //find approximate diagonal weights
    genmult(modelbasis,tmpmat,lamvecs,1.,0.,1);
    lamvecs %= lamvecs;
    lambda = maxcols(lamvecs);
    // std::cout<<lamvecs;
          
    //delete small diagonal weights 
    maxeig=max(maxeig,max(lambda));
    Real lmin=max(1e-6,relprec*maxeig);
    Indexmatrix delind(lambda.dim(),1);
    delind.init(0,0,Integer(0));
    for (Integer j=0;j<lambda.dim();j++){
      if (lambda(j)<lmin)
	delind.concat_below(j);
    }
    if (delind.dim()){
      lambda.delete_rows(delind,true);
      modelbasis.delete_cols(delind,true);
    }    
    if (lambda.rowdim()==0)
      continue;

    if (cb_out(2)){
      get_out()<<" activei="<<i<<" sz="<<mdi.size()<<" cols="<<modelbasis.coldim()<<" lmax="<<max(lambda)/coeff<<" coeff="<<coeff<<std::endl;
    }
   
    lambda.sqrt();
    modelbasis.scale_cols(lambda);

    //append to vecH
    vecH.concat_right(modelbasis);
  }
    

  //--- form SVD of vecH
  Indexmatrix piv;
  Integer r=vecH.QR_factor(piv);
  tmpmat=vecH.rows(Range(0,r-1));
  tmpmat.triu();
  rankadd(tmpmat,tmpsym,1.,0.);
  Real scale=1.;
  for(Integer i=1;i<tmpsym.rowdim();i++)
    scale=max(scale,tmpsym(i,i));
  tmpsym/=scale;
  int retval=tmpsym.eig(lamvecs,lambda,false);
  if ((min(lambda)<-1e-6)||(retval)){
    if (cb_out())
      get_out()<<"**** WARNING: VariableMetricSVDSelection::vecH_weighted_SVDs(............): S.eig failed and returned "<<retval<<" (order="<<tmpsym.rowdim()<<") and min eigenvalue="<<min(lambda)<<std::endl;
    
  }
  lambda*=scale;
  if (norm2(lambda)<max(1e-10,1e-6*aggrnorm)){  
    //"nothing left" after orthogonalization, no contribution here
    vecH.init(dim,0,0.);
    
    // Uspace=vecH;
    // lambda.init(0,1,0.);
    // candNewton*=-1./weightu;
    // candNewton+=y;
    
    return 0;    
  }
  //lambda-=weightu;
  Real lmaxbound=max(eps_Real*weightu,relprec*lambda(0));
  //Real lmaxbound=1e-6*lam(0);
  Integer cnt=0;
  Integer maxdir=lambda.dim();
  while ((cnt<maxdir)&&(lambda(cnt)>lmaxbound))
      cnt++;
  if (cnt<lambda.dim()){
    lambda.reduce_length(cnt);
    lamvecs.delete_cols(Range(cnt,lamvecs.coldim()-1),true);
  }
  if (lamvecs.coldim()>0){
    lamvecs.enlarge_below(dim-lamvecs.rowdim(),0.);
    vecH.Q_times(lamvecs,r);
    swap(lamvecs,vecH);
  }
  else {
    vecH.init(dim,0,0.);
  }
  if (norm2(lambda)<max(1e-10,1e-6*aggrnorm)){  
    //"nothing left" after orthogonalization, no contribution here
    vecH.init(dim,0,0.);
    
    // Uspace=vecH;
    // lambda.init(0,1,0.);
    // candNewton*=-1./weightu;
    // candNewton+=y;
    
    return 0;    
  }

  
  // Uspace=vecH;
  // lambda=lam;
  // tmpvec=candNewton;
  // rankadd(vecH,tmpsym,1.,0.,1);
  // for(Integer i=0;i<tmpsym.rowdim();i++)
  //   tmpsym(i,i)+=weightu/lambda(i);
  // if (tmpsym.Chol_factor())
  //   std::cout<<"factorizing tmpsym failed"<<std::endl;
  // genmult(vecH,candNewton,tmpmat,1.,0.,1);
  // tmpsym.Chol_solve(tmpmat);
  // genmult(vecH,tmpmat,candNewton,-1.,1.);
  // candNewton*=-1./weightu;
  // tmpvec+=candNewton*weightu+vecH*Diag(lambda)*transpose(vecH)*candNewton;
  // std::cout<<" cand_error="<<norm2(tmpvec);
  // candNewton+=y;
  // lambda+=weightu;
  

  //tmpmat.init_diag(vecH.coldim(),1.);
  //genmult(vecH,vecH,tmpmat,1.,-1.,1);
  //std::cout<<" norm2(eye(n)-Uspace^T*Uspace)="<<norm2(tmpmat);
  

  if (cb_out(2)){
    get_out()<<"\n wSVDs.dim="<<lambda.dim();
    if (lambda.dim()>0)
      get_out()<<" max(lam0)="<<max(lambda)<<" sum(lam0)="<<sum(lambda)<<" avg(vecHSq)="<<sum(lambda)/lambda.dim();
    //genmult(Vspace,vecH,tmpmat,1.,0.,1);
    //get_out()<<" orthogonality="<<norm2(tmpmat)<<std::endl;
    
    get_out()<<std::endl;
  }

  lambda.sqrt();        
  vecH.scale_cols(lambda);      

  //vecH.concat_right(Vspace*100.);

  return 0;
}
    

//*****************************************************************************
//                    VariableMetricSVDSelection::vecH_weighted_SVD
//*****************************************************************************

int VariableMetricSVDSelection::vecH_weighted_SVD(Matrix& vecH,
				      Real& relprec,
				      const MinorantBundle& bundle,
				      const MinorantBundle& model,
				      const Matrix& modelcoeff,
				      const MinorantPointer& aggregate,
				      Real function_factor,
				      Integer y_id,
				      const Matrix& y,
				      Real weightu,
				      Real model_maxviol)
{
  assert(function_factor>0.);
  assert(model_maxviol>0.);
  assert(weightu>0.);

  if (cb_out(2)){
    get_out()<<"\n wSVDdscal";
  }

    Integer dim=y.dim();
  Integer msz=Integer(model.size());   // number of minorants in the model (includes aggr)
  Integer bsz=Integer(bundle.size());

  vecH.init(dim,max(msz,bsz),0.);
  Uspace.init(dim,0,0.);

  relprec=1.;

  Real aggrval=aggregate.evaluate(y_id,y)/function_factor;
  Real aggrsqu=aggregate.dual_norm_squared()/function_factor/function_factor;
  Real aggrnorm=std::sqrt(aggrsqu);
  Real coeffzero=max(1e-8*function_factor,1e-4*max(modelcoeff));
  
  //---- collect the subspace spanned by the difference vectors of the active model minorants 
  Real maxval=aggrval;                // will be set to the model value in y
  vecH.init(0,0,0.);
  Matrix modelbasis(dim,msz);
  Matrix modelnorm(msz,1,0.);
  modelbasis.init(dim,0,0.);
  Matrix tmpmat,tmpvec;
  Symmatrix tmpsym;
  Matrix lam; 


  if (model.size()>1){

    //determine the matrix of difference vectors to the aggregate
    Real dummy;
    tmpmat.newsize(dim,1); chk_set_init(tmpmat,1);
    aggregate.get_minorant(dummy,tmpmat,0,1./function_factor);
    for (Integer i=0;i<msz;i++){ 
      assert(!model[unsigned(i)].empty());
      Real v=model[unsigned(i)].evaluate(y_id,y);
      if (v>maxval)
	maxval=v;
      if (modelcoeff(i)<coeffzero)
        continue;
      modelnorm(i)=std::sqrt(model[unsigned(i)].dual_norm_squared());
      //subtract vector from column i-1 of modelbasis
      modelbasis.concat_right(tmpmat);
      model[unsigned(i)].get_minorant(dummy,modelbasis,modelbasis.coldim()-1,-1.,true);
    }     

    //determine main directions via SVD

    rankadd(modelbasis,tmpsym,1.,0.,1);

    Real scale=1.;
    for(Integer i=1;i<tmpsym.rowdim();i++)
      scale=max(scale,tmpsym(i,i));
    tmpsym/=scale;

    int retval=tmpsym.eig(tmpmat,lam,false);
    if ((min(lam)<-1e-6)||(retval)){
      if (cb_out())
	get_out()<<"**** WARNING: VariableMetricSVDSelection::vecH_weighted_SVD(............): S.eig failed and returned "<<retval<<" (order="<<tmpsym.rowdim()<<") and min eigenvalue="<<min(lam)<<std::endl;
      
    }

    lam*=scale;
    if (norm2(lam)<max(1e-10,1e-6*std::sqrt(aggregate.dual_norm_squared())/function_factor)){  
      //"nothing left" after orthogonalization, no contribution here
      modelbasis.init(dim,0,0.);
    }
    else {
      Real lmaxbound=1e-3*lam(0);
      Integer cnt=1;
      Integer maxdir=lam.dim();
      while ((cnt<maxdir)&&(lam(cnt)>lmaxbound))
	cnt++;
      if (cnt<lam.dim()){
	lam.reduce_length(cnt);
	tmpmat.delete_cols(Range(cnt,modelbasis.coldim()-1));
      }
      lam.sqrt();         // now lam == Sigma
      lam.inv();
      tmpmat.scale_cols(lam);  //now modelbasis <- Q*Sigma^{-1} 
   
      genmult(modelbasis,tmpmat,vecH);  // vecH <- P = V * Q*Sigma^{-1}
      swap(modelbasis,vecH);
      vecH.init(dim,0,0.);
    }

    

  }//endif(model.size>1)

  
  Vspace=modelbasis;
  candNewton.newsize(dim,1); chk_set_init(candNewton,1);
  Real dummy;
  aggregate.get_minorant(dummy,candNewton,0,1./function_factor);
  genmult(Vspace,candNewton,tmpmat,1.,0.,1);
  genmult(Vspace,tmpmat,candNewton,-1.,1.);
  

  if (cb_out(2)){
    get_out()<<" modelsize="<<model.size()<<" SVDmodeldim="<<modelbasis.coldim()<<std::endl;
  }

  assert(maxval>=aggrval-1e10*(1+std::fabs(aggrval)));
  Real violation_eps=.9*model_maxviol/function_factor;
  Real mindiff=.01*violation_eps;
  relprec=min(1e-3,max(1e-12,mindiff)/(1.+CH_Matrix_Classes::abs(maxval)));

  if (cb_out(2)){
    get_out().precision(12);
    get_out()<<" violeps="<<violation_eps<<" maxval-aggrval="<<maxval-aggrval<<" mindiff="<<mindiff<<" relprec="<<relprec<<std::endl;
  }

  //---- select active subgradients in model
  MinorantBundle active;
  Matrix activenorm(msz,1);activenorm.init(0,1,0.);
  Matrix activecoeff(msz,1);activecoeff.init(0,1,0.);
  if (modelbasis.coldim()==0){
    //use the minorant having maximum norm*coeff
    Real maxval=-1.;
    Integer maxi=-1;
    for(Integer i=0;i<msz;i++){
      Real d=modelnorm(i)*modelcoeff(i);
      if (d>maxval){
	maxval=d;
	maxi=i;
      }
    }
    assert(maxi>=0);
    active.push_back(model[unsigned(maxi)]);
    activenorm.concat_below(modelnorm(maxi));
    activecoeff.concat_below(modelcoeff(maxi));
  }
  else {
    //select the modelbasis.coldim()+1 minorants having largest projection 
    tmpmat.init(dim,1,0.);
    activecoeff=modelcoeff;
    for (Integer i=0;i<msz;i++){ 
      model[unsigned(i)].get_minorant(dummy,tmpmat,0,modelcoeff(i));
      genmult(modelbasis,tmpmat,tmpvec,1.,0.,1);
      activecoeff(i)*=norm2(tmpvec);
    }
    Indexmatrix sind;
    sortindex(activecoeff,sind,false);
    activecoeff.init(0,1,0.);
    //make sure that they are sufficiently dissimilar 
    Integer cnt=0;
    for(Integer i=0;i<msz;i++){
      //find the best matching vector so far
      Real maxval=-2.;
      Integer maxj=-1;
      for(Integer j=0;j<cnt;j++){
	Real d=model[unsigned(i)].ip(active[unsigned(j)])/(activenorm(j)*modelnorm(i));
	if (d>maxval){
	  maxval=d;
	  maxj=j;
	}
      }
      if (maxval>.7){
	activecoeff(maxj)+=modelcoeff(i);
      }
      else {
	cnt++;
	active.push_back(model[unsigned(i)]);
	activecoeff.concat_below(modelcoeff(i));
	activenorm.concat_below(modelnorm(i));
	if (cnt>=modelbasis.coldim()+1)
	  break;
      }
    }
  }     


  /*
  Indexmatrix maxipind(Integer(model.size()),1,-1);
  Real maxcoeff=-1.;
  Integer nexti=-1;
  maxcoeff=max(modelcoeff,&nexti);
  assert((0<=nexti)&&(unsigned(nexti)<model.size()));
  Matrix maxip(modelcoeff.rowdim(),1,-1.);
  while (nexti>=0){
    std::cout<<" modeli="<<nexti<<" maxip="<<maxip(nexti)<<" coeff="<<maxcoeff<<std::endl;
    maxipind(nexti)=-1;
    active.push_back(model[unsigned(nexti)]);
    maxip(nexti)=2.;
    Real mnrtnrm=modelnorm(nexti);
    activenorm.concat_below(mnrtnrm);
    activecoeff.concat_below(maxcoeff);
    const MinorantPointer& mnrt=model[unsigned(nexti)];
    nexti=-1;
    maxcoeff=-1;
    for(Integer i=0;i<modelcoeff.rowdim();i++){
      if ((maxip(i)>=1.)||(modelcoeff(i)<coeffzero))
	continue;
      Real ipval=model[unsigned(i)].ip(mnrt)/mnrtnrm/modelnorm(i);
      if (ipval>maxip(i)){
	maxip(i)=ipval;
	maxipind(i)=activecoeff.rowdim()-1;
      }
      if ((maxip(i)<.95)&&(modelcoeff(i)>maxcoeff)){
	maxcoeff=modelcoeff(i);
	nexti=i;
      }
    }
  }
  for (Integer i=0;i<modelcoeff.rowdim();i++){
    if (maxipind(i)>0){
      std::cout<<" modeli="<<i<<" maxip(i)="<<maxip(i)<<" to "<<maxipind(i)<<std::endl;
      activecoeff(maxipind(i))+=modelcoeff(i);
    }
  }
  assert(sum(activecoeff)-sum(modelcoeff)<1e-12*(1.+sum(modelcoeff)));
  */
  
  //---- assign minorants to active subgradients
  std::vector<std::vector<Integer> > assignedmnrts(active.size());
  for(unsigned int i=0;i<bundle.size();i++){
    Real maxip=min_Real;
    unsigned int bestj((unsigned int)active.size());
    const MinorantPointer& mnrt=bundle[i];
    if (mnrt.empty()||mnrt.aggregate())
      continue;
    Real v=mnrt.evaluate(y_id,y);
    if ((v>=aggrval-1e-6*(maxval-aggrval))||(v<aggrval-1000.*(maxval-aggrval))) 
      continue;
    for(unsigned int j=0;j<active.size();j++){
      if (activecoeff(Integer(j))<coeffzero)
        continue;
      Real ip=active[j].ip(mnrt)/activenorm(Integer(j));
      if (ip>maxip){
	maxip=ip;
	bestj=j;
      }
    }
    if(bestj<active.size())
      assignedmnrts[bestj].push_back(Integer(i));
  }

  //---- for each active model vector collect the information of assigned minorants
  vecH.init(dim,0,0.);

  for(unsigned int activei=0;activei<active.size();activei++){
    Real coeff=activecoeff(Integer(activei));
    if (coeff<coeffzero)
      continue;
    std::vector<Integer>& mnrts=assignedmnrts[activei];
    if (mnrts.size()==0)
      continue;
    Matrix bundlediff(dim,Integer(mnrts.size()),0.);
    Integer bcnt=0;
    for(unsigned int mnrti=0;mnrti<mnrts.size();mnrti++){
      unsigned long bi=unsigned(mnrts[mnrti]);
      Real v=bundle[bi].evaluate(y_id,y);
      Real delta=std::sqrt(2.*(max(mindiff,aggrval-v)+violation_eps));
      Real dummy;
      active[activei].get_minorant(dummy,bundlediff,bcnt,1./delta,true);
      bundle[bi].get_minorant(dummy,bundlediff,bcnt,-1./delta,true);
      bcnt++;
    }
    bundlediff.delete_cols(Range(bcnt,bundlediff.coldim()-1));    
    genmult(modelbasis,bundlediff,tmpmat,1.,0.,1);
    genmult(modelbasis,tmpmat,bundlediff,-1.,1.);
    rankadd(bundlediff,tmpsym,1.,0.,1);
    Real scale=1.;
    for(Integer i=1;i<tmpsym.rowdim();i++)
      scale=max(scale,tmpsym(i,i));
    tmpsym/=scale;
    int retval=tmpsym.eig(tmpmat,lam,false);
    if ((min(lam)<-1e-6)||(retval)){
      if (cb_out())
	get_out()<<"**** WARNING: VariableMetricSVDSelection::vecH_weighted_SVD(............): S.eig failed and returned "<<retval<<" (order="<<tmpsym.rowdim()<<") and min eigenvalue="<<min(lam)<<std::endl;
      
    }
    lam*=scale;
    if (norm2(lam)<max(1e-10,1e-6*aggrnorm)){  
      //"nothing left" after orthogonalization, no contribution here
      continue;
    }
    Real lmaxbound=relprec*lam(0);
    //Real lmaxbound=1e-6*lam(0);
    Integer cnt=1;
    Integer maxdir=lam.dim();
    while ((cnt<maxdir)&&(lam(cnt)>lmaxbound))
      cnt++;
    if (cnt<lam.dim()){
      lam.reduce_length(cnt);
      tmpmat.delete_cols(Range(cnt,tmpmat.coldim()-1));
    }
    if (lam.rowdim()==0)
      continue;
    lam.sqrt();         // now lam == Sigma
    lam.inv();
    tmpmat.scale_cols(lam);  //now modelbasis <- Q*Sigma^{-1} 
    
    Matrix local_vecH;
    genmult(bundlediff,tmpmat,local_vecH);  // vecH <- P = V * Q*Sigma^{-1}

    genmult(local_vecH,bundlediff,tmpmat,1.,0.,1);
    tmpmat %= tmpmat;
    lam = maxcols(tmpmat);
    lam*=coeff;
          
    Real lmin=max(1e-6,relprec*max(lam));
    Indexmatrix delind(lam.dim(),1);
    delind.init(0,0,Integer(0));
    for (Integer i=0;i<lam.dim();i++){
      if (lam(i)<lmin)
	delind.concat_below(i);
    }
    if (delind.dim()){
      lam.delete_rows(delind,true);
      local_vecH.delete_cols(delind,true);
    }    
    if (lam.rowdim()==0)
      continue;

    if (cb_out(2))
      get_out()<<" activei="<<activei<<" sz="<<mnrts.size()<<" cols="<<local_vecH.coldim()<<" lmax="<<max(lam)/coeff<<" coeff="<<coeff<<std::endl;
   
    lam.sqrt();
    local_vecH.scale_cols(lam);

    vecH.concat_right(local_vecH);
  }

  //--- form SVD of vecH
  Indexmatrix piv;
  Integer r=vecH.QR_factor(piv);
  tmpmat=vecH.rows(Range(0,r-1));
  tmpmat.triu();
  rankadd(tmpmat,tmpsym,1.,0.);
  Real scale=1.;
  for(Integer i=1;i<tmpsym.rowdim();i++)
    scale=max(scale,tmpsym(i,i));
  tmpsym/=scale;
  int retval=tmpsym.eig(tmpmat,lam,false);
  if ((min(lam)<-1e-6)||(retval)){
    if (cb_out())
      get_out()<<"**** WARNING: VariableMetricSVDSelection::vecH_orthogonal_to_model(............): S.eig failed and returned "<<retval<<" (order="<<tmpsym.rowdim()<<") and min eigenvalue="<<min(lam)<<std::endl;
    
  }
  lam*=scale;
  if (norm2(lam)<max(1e-10,1e-6*aggrnorm)){  
    //"nothing left" after orthogonalization, no contribution here
    vecH.init(dim,0,0.);
    
    Uspace=vecH;
    lambda.init(0,1,0.);
    candNewton*=-1./weightu;
    candNewton+=y;
    
    return 0;    
  }
  lam-=weightu;
  Real lmaxbound=max(eps_Real*weightu,relprec*lam(0));
  //Real lmaxbound=1e-6*lam(0);
  Integer cnt=0;
  Integer maxdir=lam.dim();
  while ((cnt<maxdir)&&(lam(cnt)>lmaxbound))
      cnt++;
  if (cnt<lam.dim()){
    lam.reduce_length(cnt);
    tmpmat.delete_cols(Range(cnt,tmpmat.coldim()-1));
  }
  if (tmpmat.coldim()>0){
    tmpmat.enlarge_below(dim-tmpmat.rowdim(),0.);
    vecH.Q_times(tmpmat,r);
    swap(tmpmat,vecH);
  }
  else {
    vecH.init(dim,0,0.);
  }
  if (norm2(lam)<max(1e-10,1e-6*aggrnorm)){  
    //"nothing left" after orthogonalization, no contribution here
    vecH.init(dim,0,0.);
    
    Uspace=vecH;
    lambda.init(0,1,0.);
    candNewton*=-1./weightu;
    candNewton+=y;
    
    return 0;    
  }

  
  Uspace=vecH;
  lambda=lam;
  tmpvec=candNewton;
  rankadd(vecH,tmpsym,1.,0.,1);
  for(Integer i=0;i<tmpsym.rowdim();i++)
    tmpsym(i,i)+=weightu/lambda(i);
  if (tmpsym.Chol_factor())
    if (cb_out())
      get_out()<<"**** WARNING: VariableMetricSVDSelection::vecH_weighted_SVD(............): Chol_factor falied "<<std::endl;
  genmult(vecH,candNewton,tmpmat,1.,0.,1);
  tmpsym.Chol_solve(tmpmat);
  genmult(vecH,tmpmat,candNewton,-1.,1.);
  candNewton*=-1./weightu;
  tmpvec+=candNewton*weightu+vecH*Diag(lambda)*transpose(vecH)*candNewton;
  if (cb_out(2))
    get_out()<<" cand_error="<<norm2(tmpvec);
  candNewton+=y;
  lambda+=weightu;
  

  //tmpmat.init_diag(vecH.coldim(),1.);
  //genmult(vecH,vecH,tmpmat,1.,-1.,1);
  //std::cout<<" norm2(eye(n)-Uspace^T*Uspace)="<<norm2(tmpmat);
  
  lam.sqrt();         // now lam == Sigma
  vecH.scale_cols(lam);      


  if (cb_out(2)){
    get_out()<<"\n orthlam0.dim="<<lam.dim();
    if (lam.dim()>0)
      get_out()<<" max(lam0)="<<max(lam)<<" sum(lam0)="<<sum(lam)<<" avg(vecHSq)="<<sum(lam)/lam.dim();
    //genmult(Vspace,vecH,tmpmat,1.,0.,1);
    //get_out()<<" orthogonality="<<norm2(tmpmat)<<std::endl;
    
    get_out()<<std::endl;
  }

  //vecH.concat_right(Vspace*100.);

  return 0;
}
    

//*****************************************************************************
//                    VariableMetricSVDSelection::add_variable_metric
//*****************************************************************************

int VariableMetricSVDSelection::vecH_orthogonal_to_modelSVD(Matrix& vecH,
						Real& relprec,
						const MinorantBundle& bundle,
						const MinorantBundle& model,
						const Matrix& /* modelcoeff */,
						const MinorantPointer& aggregate,
						Real function_factor,
						Integer y_id,
						const Matrix& y,
						Real weightu,
						Real model_maxviol)
{
  assert(function_factor>0.);
  assert(model_maxviol>0.);
  assert(weightu>0.);

  if (cb_out(2)){
    get_out()<<"\n orthmodelSVDdscal"; 
  }
  
  Integer dim=y.dim();
  Integer msz=Integer(model.size());   // number of minorants in the model (includes aggr)
  Integer bsz=Integer(bundle.size());

  vecH.init(dim,max(msz,bsz),0.);
  Uspace.init(dim,0,0.);

  relprec=1.;

  Real aggrval=aggregate.evaluate(y_id,y)/function_factor;
  //Real aggrsqu=aggregate.dual_norm_squared()/function_factor/function_factor;
  //Real aggrnorm=std::sqrt(aggrsqu);
  
  //---- collect the subspace spanned by the difference vectors of the active model minorants 
  Real maxval=aggrval;                // will be set to the model value in y
  vecH.init(0,0,0.);
  Matrix modelbasis(dim,msz);
  modelbasis.init(dim,0,0.);
  Matrix tmpmat;
  Symmatrix tmpsym;
  Matrix lam; 

  if (model.size()>1){

    //determine the matrix of difference vectors to the aggregate
    Real dummy;
    tmpmat.newsize(dim,1); chk_set_init(tmpmat,1);
    aggregate.get_minorant(dummy,tmpmat,0,1./function_factor);
    for (Integer i=0;i<msz;i++){ 
      assert(!model[unsigned(i)].empty());
      Real v=model[unsigned(i)].evaluate(y_id,y);
      if (v>maxval)
	maxval=v;
      //if (modelcoeff(i)<1e-10)
      //  continue;
      //subtract vector from column i-1 of modelbasis
      modelbasis.concat_right(tmpmat);
      model[unsigned(i)].get_minorant(dummy,modelbasis,modelbasis.coldim()-1,-1.,true);
    }     

    //determine main directions via SVD

    rankadd(modelbasis,tmpsym,1.,0.,1);

    Real scale=1.;
    for(Integer i=1;i<tmpsym.rowdim();i++)
      scale=max(scale,tmpsym(i,i));
    tmpsym/=scale;

    int retval=tmpsym.eig(tmpmat,lam,false);
    if ((min(lam)<-1e-6)||(retval)){
      if (cb_out())
	get_out()<<"**** WARNING: VariableMetricSVDSelection::vecH_orthogonal_to_modelSVD(............): S.eig failed and returned "<<retval<<" (order="<<tmpsym.rowdim()<<") and min eigenvalue="<<min(lam)<<std::endl;
      
    }

    lam*=scale;
    if (norm2(lam)<max(1e-10,1e-6*std::sqrt(aggregate.dual_norm_squared())/function_factor)){  
      //"nothing left" after orthogonalization, no contribution here
      modelbasis.init(dim,0,0.);
    }
    else {
      Real lmaxbound=1e-6*lam(0);
      Integer cnt=1;
      Integer maxdir=lam.dim();
      while ((cnt<maxdir)&&(lam(cnt)>lmaxbound))
	cnt++;
      if (cnt<lam.dim()){
	lam.reduce_length(cnt);
	tmpmat.delete_cols(Range(cnt,modelbasis.coldim()-1));
      }
      lam.sqrt();         // now lam == Sigma
      lam.inv();
      tmpmat.scale_cols(lam);  //now modelbasis <- Q*Sigma^{-1} 
   
      genmult(modelbasis,tmpmat,vecH);  // vecH <- P = V * Q*Sigma^{-1}
      swap(modelbasis,vecH);
      vecH.init(dim,0,0.);
    }

    

  }//endif(model.size>1)

  Vspace=modelbasis;
  candNewton.newsize(dim,1); chk_set_init(candNewton,1);
  Real dummy;
  aggregate.get_minorant(dummy,candNewton,0,1./function_factor);
 

  if (cb_out(2)){
    get_out()<<" modelsize="<<model.size()<<" SVDmodeldim="<<modelbasis.coldim()<<std::endl;
  }

  assert(maxval>=aggrval-1e10*(1+std::fabs(aggrval)));
  Real violation_eps=.9*model_maxviol/function_factor;
  Real mindiff=.01*violation_eps;
  relprec=min(1e-3,max(1e-12,mindiff)/(1.+CH_Matrix_Classes::abs(maxval)));

  if (cb_out(2)){
    get_out().precision(12);
    get_out()<<" violeps="<<violation_eps<<" maxval-aggrval="<<maxval-aggrval<<" mindiff="<<mindiff<<" relprec="<<relprec<<std::endl;
  }

  //---- collect the difference vectors of bundle vectors and the aggregate  
  Matrix bundlediff(dim,bsz,0.);
  Integer bcnt=0;
  for (Integer i=0;i<bsz;i++){ 
    if (bundle[unsigned(i)].empty()||bundle[unsigned(i)].aggregate())
      continue;
    Real v=bundle[unsigned(i)].evaluate(y_id,y);
    if ((v>=aggrval-1e-6*(maxval-aggrval))||(v<aggrval-100.*(maxval-aggrval))) 
      continue;
    Real delta=std::sqrt(2.*(max(mindiff,aggrval-v)+violation_eps));
    Real dummy;
    aggregate.get_minorant(dummy,bundlediff,bcnt,1./delta/function_factor,true);
    bundle[unsigned(i)].get_minorant(dummy,bundlediff,bcnt,-1./delta,true);
    bcnt++;
  }
  bundlediff.delete_cols(Range(bcnt,bundlediff.coldim()-1));
  if (cb_out(2))
    get_out()<<" consider "<<bundlediff.coldim()<<" of "<<bsz<<std::endl;
  //std::cout<<" maxabs="<<max(abs(bundlediff));

  //---- compute the SVD for the vectors projected onto the orthogonal complement of modelbasis
  genmult(modelbasis,bundlediff,tmpmat,1.,0.,1);
  genmult(modelbasis,tmpmat,bundlediff,-1.,1.);

  //std::cout<<" maxabsQ="<<max(abs(bundlediff));

  rankadd(bundlediff,tmpsym,1.,0.,1);


  //std::cout<<" norm2(tmpsym-transpose(tmpmat)*tmpmat)="<<norm2(tmpsym-transpose(tmpmat)*tmpmat)<<std::endl;

  Real scale=1.;
  for(Integer i=1;i<tmpsym.rowdim();i++)
    scale=max(scale,tmpsym(i,i));
  tmpsym/=scale;

  int retval=tmpsym.eig(tmpmat,lam,false);

  //std::cout<<" norm2(tmpsym-tmpmat*Diag(lam)*transpose(tmpmat))="<<norm2(tmpsym-tmpmat*Diag(lam)*transpose(tmpmat))<<std::endl;

  if ((min(lam)<-1e-6)||(retval)){
    if (cb_out())
      get_out()<<"**** WARNING: VariableMetricSVDSelection::vecH_orthogonal_to_modelSVD(............): S.eig failed and returned "<<retval<<" (order="<<tmpsym.rowdim()<<") and min eigenvalue="<<min(lam)<<std::endl;

  }
  lam*=scale;
  if (norm2(lam)<max(1e-10,1e-6*std::sqrt(aggregate.dual_norm_squared())/function_factor)){  
    //"nothing left" after orthogonalization, no contribution here
    vecH.init(dim,0,0.);
    Uspace=vecH;
    lambda.init(0,1,0.);
    return 0;
  }
  Real lmaxbound=relprec*lam(0);
  //Real lmaxbound=1e-6*lam(0);
  Integer cnt=1;
  Integer maxdir=lam.dim();
  while ((cnt<maxdir)&&(lam(cnt)>lmaxbound))
    cnt++;
  if (cnt<lam.dim()){
    lam.reduce_length(cnt);
    tmpmat.delete_cols(Range(cnt,tmpmat.coldim()-1));
  }
  lam.sqrt();         // now lam == Sigma
  lam.inv();
  tmpmat.scale_cols(lam);  //now modelbasis <- Q*Sigma^{-1} 
   
  genmult(bundlediff,tmpmat,vecH);  // vecH <- P = V * Q*Sigma^{-1}


  //---- compute the eigenvalue approximations
  if (vecH.coldim()>0){

    genmult(vecH,bundlediff,tmpmat,1.,0.,1);
    tmpmat %= tmpmat;
    lam = maxcols(tmpmat);
    lam-=weightu;
          
    Real lmin=max(1e-6,relprec*max(lam));
    Indexmatrix delind(lam.dim(),1);
    delind.init(0,0,Integer(0));
    for (Integer i=0;i<lam.dim();i++){
      if (lam(i)<lmin)
	delind.concat_below(i);
    }
    if (delind.dim()){
      lam.delete_rows(delind,true);
      vecH.delete_cols(delind,true);
    }
   
    Uspace=vecH;
    lambda=lam;
    genmult(Vspace,candNewton,tmpmat,1.,0.,1);
    genmult(Vspace,tmpmat,candNewton,-1.,1.);
    Matrix tmpvec=candNewton;
    rankadd(vecH,tmpsym,1.,0.,1);
    for(Integer i=0;i<tmpsym.rowdim();i++)
      tmpsym(i,i)+=weightu/lambda(i);
    if (tmpsym.Chol_factor())
      if (cb_out())
	get_out()<<"**** WARNING: VariableMetricSVDSelection::vecH_orthogonal_to_modelSVD(............): Chol_factor falied "<<std::endl;
    genmult(vecH,candNewton,tmpmat,1.,0.,1);
    tmpsym.Chol_solve(tmpmat);
    genmult(vecH,tmpmat,candNewton,-1.,1.);
    candNewton*=-1./weightu;
    tmpvec+=candNewton*weightu+vecH*Diag(lambda)*transpose(vecH)*candNewton;
    if (cb_out(2))
      get_out()<<" cand_error="<<norm2(tmpvec);
    candNewton+=y;
    lambda+=weightu;
    
    tmpmat.init_diag(vecH.coldim(),1.);
    genmult(vecH,vecH,tmpmat,1.,-1.,1);
    if (cb_out(2))
      get_out()<<" norm2(eye(n)-Uspace^T*Uspace)="<<norm2(tmpmat);


    lam.sqrt();
    vecH.scale_cols(lam);

    if (cb_out(2)){
      get_out()<<"\n orthlam0.dim="<<lam.dim();
      if (lam.dim()>0)
	get_out()<<" max(lam0)="<<max(lam)<<" sum(lam0)="<<sum(lam)<<" avg(vecHSq)="<<sum(lam)/lam.dim();
      genmult(Vspace,vecH,tmpmat,1.,0.,1);
      get_out()<<" orthogonality="<<norm2(tmpmat)<<std::endl;

      get_out()<<std::endl;
    }

  }
  else {
    Uspace=vecH;
    lambda.init(0,1,0.);
  }
  //vecH.concat_right(Vspace*100.);

  return 0;
}
    

//*****************************************************************************
//                    VariableMetricSVDSelection::add_variable_metric
//*****************************************************************************

int VariableMetricSVDSelection::vecH_orthogonal_to_model(Matrix& vecH,
					     Real& relprec,
					     const MinorantBundle& bundle,
					     const MinorantBundle& model,
					     const Matrix& modelcoeff,
					     const MinorantPointer& aggregate,
					     Real function_factor,
					     Integer y_id,
					     const Matrix& y,
					     Real /* weightu */,
					     Real model_maxviol)
{
  assert(function_factor>0.);
  assert(model_maxviol>0.);

  if (cb_out(2)){
    get_out()<<"\n orthmodeldscal"; 
  }
  
  Integer dim=y.dim();
  Integer msz=Integer(model.size());   // number of minorants in the model (includes aggr)
  Integer bsz=Integer(bundle.size());

  vecH.init(dim,max(msz,bsz),0.);
  Uspace.init(dim,0,0.);

  relprec=1.;

  Real aggrval=aggregate.evaluate(y_id,y)/function_factor;
  //Real aggrsqu=aggregate.dual_norm_squared()/function_factor/function_factor;
  //Real aggrnorm=std::sqrt(aggrsqu);
  
  //---- collect the subspace spanned by the difference vectors of the active model minorants 
  Real maxval=aggrval;                // will be set to the model value in y
  Matrix modelval(Integer(model.size()),1,0.);
  vecH.init(0,0,0.);
  Matrix tmpmat;
  Matrix modelbasis(dim,msz,0.);
  modelbasis.init(dim,0,0.);
  for (Integer i=0;i<msz;i++){ 
    if (model[unsigned(i)].empty()){
      if (cb_out())
	get_out()<<"**** WARNING: VariableMetricSVDSelection::vecH_orthogonal_to_model(............): model["<<i<<"] is empty"<<std::endl;
      continue;
    }
    Real v=model[unsigned(i)].evaluate(y_id,y);
    modelval(i)=v;
    if (v>maxval)
      maxval=v;
    if (modelcoeff(i)<1e-10)
      continue;
    if (tmpmat.coldim()==0){
      Real dummy;
      tmpmat.init(dim,1,0.);
      model[unsigned(i)].get_minorant(dummy,tmpmat,0,1.,true);
    }
    else {
      //subtract vector from column i-1 of modelbasis
      Real dummy;
      modelbasis.concat_right(tmpmat);
      model[unsigned(i)].get_minorant(dummy,modelbasis,modelbasis.coldim()-1,-1.,true);
    }     
  }

  Indexmatrix tmppiv;
  //Integer modeldim=modelbasis.QR_factor(tmppiv,1e-6);
  Integer modeldim=modelbasis.QR_factor_relpiv(tmppiv,1e-3);//*min(1.,aggrsqu));
  if (cb_out(2)){
    get_out()<<"\n coldim="<<modelbasis.coldim()<<" modeldim="<<modeldim<<std::endl;
  }
  
  Vspace.init(dim,modeldim,0.);
  for(Integer i=0;i<modeldim;i++)
    Vspace(i,i)=1;
  modelbasis.Q_times(Vspace,modeldim);

  assert(maxval>=aggrval-1e10*(1+std::fabs(aggrval)));
  //assert(weightu>0.);
  Real violation_eps=.9*model_maxviol/function_factor;
  //Real maxmult=1e6*weightu/function_factor;
  Real mindiff=.01*violation_eps;
  //Real normbnd=1e-4*max(1.,aggrsqu);
  relprec=min(1e-3,max(1e-12,mindiff)/(1.+CH_Matrix_Classes::abs(maxval)));

  

  //---- collect the difference vectors of bundle vectors and the aggregate  
  Matrix bundlediff(dim,bsz,0.);
  Integer bcnt=0;
  for (Integer i=0;i<bsz;i++){ 
    if (bundle[unsigned(i)].empty()||bundle[unsigned(i)].aggregate())
      continue;
    Real v=bundle[unsigned(i)].evaluate(y_id,y);
    if ((v>aggrval)||(v<aggrval-100.*(maxval-aggrval))) 
      continue;
    Real delta=std::sqrt(2.*(max(0.,aggrval-v)+violation_eps));
    Real dummy;
    aggregate.get_minorant(dummy,bundlediff,bcnt,1./delta/function_factor,true);
    bundle[unsigned(i)].get_minorant(dummy,bundlediff,bcnt,-1./delta,true);
    bcnt++;
  }
  bundlediff.delete_cols(Range(bcnt,bundlediff.coldim()-1));
  if (cb_out(2))
    get_out()<<" consider "<<bundlediff.coldim()<<" of "<<bsz<<std::endl;
  //std::cout<<" maxabs="<<max(abs(bundlediff));

  //---- compute the SVD for the vectors projected onto the orthogonal complement of modelbasis
  tmpmat=bundlediff;
  modelbasis.Qt_times(tmpmat,modeldim);
  for(Integer j=0;j<bundlediff.coldim();j++){
    for(Integer i=0;i<modeldim;i++){
      tmpmat(i,j)=0.;
    }
  }
  modelbasis.Q_times(tmpmat,modeldim);

  //std::cout<<" maxabsQ="<<max(abs(bundlediff));

  Symmatrix tmpsym;
  rankadd(tmpmat,tmpsym,1.,0.,1);


  //std::cout<<" norm2(tmpsym-transpose(tmpmat)*tmpmat)="<<norm2(tmpsym-transpose(tmpmat)*tmpmat)<<std::endl;

  Real scale=1.;
  for(Integer i=1;i<tmpsym.rowdim();i++)
    scale=max(scale,tmpsym(i,i));
  tmpsym/=scale;

  Matrix lam; 
  int retval=tmpsym.eig(modelbasis,lam,false);

  //std::cout<<" norm2(tmpsym-tmpmat*Diag(lam)*transpose(tmpmat))="<<norm2(tmpsym-tmpmat*Diag(lam)*transpose(tmpmat))<<std::endl;

  if ((min(lam)<-1e-6)||(retval)){
    if (cb_out())
      get_out()<<"**** WARNING: VariableMetricSVDSelection::vecH_orthogonal_to_model(............): S.eig failed and returned "<<retval<<" (order="<<tmpsym.rowdim()<<") and min eigenvalue="<<min(lam)<<std::endl;

  }
  lam*=scale;
  if (norm2(lam)<max(1e-10,1e-6*std::sqrt(aggregate.dual_norm_squared())/function_factor)){  
    //"nothing left" after orthogonalization, no contribution here
    vecH.init(dim,0,0.);
    return 0;
  }
  Real lmaxbound=relprec*lam(0);
  //Real lmaxbound=1e-6*lam(0);
  Integer cnt=1;
  Integer maxdir=lam.dim();
  while ((cnt<maxdir)&&(lam(cnt)>lmaxbound))
    cnt++;
  if (cnt<lam.dim()){
    lam.reduce_length(cnt);
    modelbasis.delete_cols(Range(cnt,modelbasis.coldim()-1));
  }
  lam.sqrt();         // now lam == Sigma
  lam.inv();
  modelbasis.scale_cols(lam);  //now modelbasis <- Q*Sigma^{-1} 
   
  genmult(tmpmat,modelbasis,vecH);  // vecH <- P = V * Q*Sigma^{-1}


  //---- compute the eigenvalue approximations
  if (vecH.coldim()>0){

    genmult(vecH,bundlediff,modelbasis,1.,0.,1);
    modelbasis %= modelbasis;
    lam = maxcols(modelbasis);
    lam.sqrt();
         
    Real lmaxbound=relprec*max(lam);
    Indexmatrix delind(lam.dim(),1);
    delind.init(0,0,Integer(0));
    for (Integer i=0;i<lam.dim();i++){
      if (lam(i)<lmaxbound)
	delind.concat_below(i);
    }
    if (delind.dim()){
      lam.delete_rows(delind,true);
      vecH.delete_cols(delind,true);
    }
   
    Uspace=vecH;

    tmpmat.init_diag(vecH.coldim(),1.);
    genmult(vecH,vecH,tmpmat,1.,-1.,1);
    if (cb_out(2))
      get_out()<<" norm2(eye(n)-Uspace^T*Uspace)="<<norm2(tmpmat);


    lambda=lam;
    lambda%=lambda;

    vecH.scale_cols(lam);

    if (cb_out(2)){
      get_out()<<"\n orthlam0.dim="<<lam.dim();
      if (lam.dim()>0)
	get_out()<<" max(lam0)="<<max(lam)<<" sum(lam0)="<<sum(lam)<<" avg(vecHSq)="<<sum(lam)/lam.dim();
      genmult(Vspace,vecH,tmpmat,1.,0.,1);
      get_out()<<" orthogonality="<<norm2(tmpmat)<<std::endl;

      get_out()<<std::endl;
    }
  }


  return 0;
}
    

//*****************************************************************************
//                    VariableMetricSVDSelection::add_variable_metric
//*****************************************************************************

int VariableMetricSVDSelection::vecH_by_model(Matrix& vecH,
				  Real& relprec,
				  const MinorantBundle& bundle,
				  const MinorantBundle& model,
				  const MinorantPointer& aggregate,
				  Real function_factor,
				  Integer y_id,
				  const Matrix& y,
				  Real weightu,
				  Real model_maxviol)
{
  assert(function_factor>0.);
  assert(model_maxviol>0.);

  Integer dim=y.dim();
  vecH.init(dim,0,0.);
  relprec=1.;

  Real aggrval=aggregate.evaluate(y_id,y)/function_factor;
  Real aggrsqu=aggregate.dual_norm_squared()/function_factor/function_factor;
  Real aggrnorm=std::sqrt(aggrsqu);
  
  Real maxval=min_Real;
  Integer msz=Integer(model.size());
  Matrix modelval(Integer(model.size()),1,0.);
  Matrix modelsqu(Integer(model.size()),1,0.);
  for (Integer i=0;i<msz;i++){ 
    if (model[unsigned(i)].empty()){
      if (cb_out())
	get_out()<<"**** WARNING: VariableMetricSVDSelection::add_variable_metric(............): model["<<i<<"] is empty"<<std::endl;
      continue;
    }
    modelsqu(i)=model[unsigned(i)].dual_norm_squared();
    Real v=model[unsigned(i)].evaluate(y_id,y);
    modelval(i)=v;
    if (v>maxval)
      maxval=v;
  }
  
  Integer bsz=Integer(bundle.size());
  Matrix val(bsz,1); chk_set_init(val,1);
  for (Integer i=0;i<bsz;i++){ 
    if (bundle[unsigned(i)].empty()||bundle[unsigned(i)].aggregate())
      continue;
    Real v=bundle[unsigned(i)].evaluate(y_id,y);
    val(i)=v;
    if (v>maxval)
      maxval=v;
  }
  if (maxval==min_Real){
    return 0;
  }

  assert(model_maxviol>0.);  
  assert(maxval>=aggrval-1e10*(1+std::fabs(aggrval)));
  assert(weightu>0.);
  Real violation_eps=.9*model_maxviol/function_factor;
  Real maxmult=1e6*weightu/function_factor;
  Real mindiff=.01*violation_eps;
  Real normbnd=1e-4*max(1.,aggrsqu);
  relprec=min(1e-3,max(1e-12,mindiff)/(1.+CH_Matrix_Classes::abs(maxval)));

  //---- for each minorant compute the multiplier; for the diagonal take 
  //---- for each coordinate direction the maximum of each resulting diagonal 

  Matrix diffvecs(dim,bsz);  //stores the difference vectors used in scaling
  diffvecs.init(dim,0,0.); 

  Integer too_good=0;
  Integer too_bad=0;
  Integer too_close=0;
  Real sumnorm=0.;
  
  for (Integer i=0;i<bsz;i++){
    if (bundle[unsigned(i)].empty()||bundle[unsigned(i)].aggregate())
      continue;

    Real gnormsqu=bundle[unsigned(i)].dual_norm_squared();    
    if (val(i)+aggrnorm*std::sqrt(gnormsqu)/weightu<aggrval+violation_eps){
      too_bad++;
      continue;
    }

    //determine values for the aggregate
    Real violoffset=violation_eps-val(i);
    Real bestdv=violoffset+aggrval;
    if (bestdv<mindiff){
      too_good++;
      bestdv=mindiff;
    }  
    Real ipval=bundle[unsigned(i)].ip(aggregate)/function_factor;   
    Real diffnormsqu=aggrsqu-2*ipval+gnormsqu;
    Real bestnormsqu=.5*diffnormsqu/bestdv;
    if (min(diffnormsqu,bestnormsqu)<normbnd){
      too_close++;
      continue;
    }
    Integer bestind=-1;
    //check for smaller vectors in the bundle
    for (Integer j=0;j<msz;j++){
      Real dvj=violoffset+modelval(j);
      if (dvj<mindiff){
	//current point already violates by more than local_violation_eps
        //continue;
	dvj=mindiff;
      }
      Real ipvalj=bundle[unsigned(i)].ip(model[unsigned(j)]);
      Real dnsj=modelsqu(j)-2*ipvalj+gnormsqu;
      Real normsqu=.5*dnsj/dvj;
      if (min(dnsj,normsqu)<normbnd){
	//too close, do not include at all
	bestind=-2;
	break;
      }
      if (normsqu<bestnormsqu){
	bestdv=dvj;
	bestnormsqu=normsqu;
	bestind=j;
      }
    }
    if (bestind<-1) //signal for "almost identical minorant in model", skip this
      continue;
    //append the difference vector
    diffvecs.enlarge_right(1,0.);
    Real dummy;
    bestdv=sqrt(2.*bestdv);
    sumnorm+=std::sqrt(bestnormsqu);
    bundle[unsigned(i)].get_minorant(dummy,diffvecs,diffvecs.coldim()-1,-1./bestdv,true);
    if (bestind<0){
      aggregate.get_minorant(dummy,diffvecs,diffvecs.coldim()-1,1./function_factor/bestdv,true);
    }
    else {
      model[unsigned(bestind)].get_minorant(dummy,diffvecs,diffvecs.coldim()-1,-1./bestdv,true);
    }
  }
 
  if (cb_out(2)){
    get_out()<<"\n bsz="<<bsz<<" relprec="<<relprec<<" mindiff="<<mindiff<<" maxmult="<<maxmult<<" tg="<<too_good<<" tb="<<too_bad<<" tc="<<too_close;
    if (diffvecs.coldim()>0)
      get_out()<<" avg(norm)="<<sumnorm/diffvecs.coldim();
  }

  //compute (approximate) SVD basis for the difference of minorants and aggr
  //vecH is P in the singular value decomposition of diffvecs = P*Sigma*Q'
  Matrix barQ;  // mostly Q in singular value decomp. diffvecs = P*Sigma*Q'
  Matrix vecHsq;  //first Sigma, later eigenvalues of H= vecH*vecHsq*vecH'
  if (diffvecs.coldim()>0){ 
    Symmatrix S;
    rankadd(diffvecs,S,1.,0.,1);
    
    int retval=S.eig(barQ,vecHsq,false);
    if ((min(vecHsq)<-1e-6)||(retval)){
      if (cb_out())
	get_out()<<"**** WARNING: VariableMetricSVDSelection::add_variable_metric(...): S.eig failed and returned "<<retval<<" (order="<<S.rowdim()<<")"<<std::endl;
    }
    Real lmaxbound=relprec*vecHsq(0);
    Integer cnt=1;
    Integer maxdir=vecHsq.dim();//   min(vecHsq.dim(),20);
    while ((cnt<maxdir)&&(vecHsq(cnt)>lmaxbound))
      cnt++;
    if (cnt<vecHsq.dim()){
      vecHsq.reduce_length(cnt);
      barQ.delete_cols(Range(cnt,barQ.coldim()-1));
    }
    vecHsq.sqrt();         // now VecHsq == Sigma
    vecHsq.inv();
    barQ.scale_cols(vecHsq);  //now barQ <- Q*Sigma^{-1} 
    
    genmult(diffvecs,barQ,vecH);  // vecH <- P = V * Q*Sigma^{-1}
    
    //orthogonalize vecH, maybe not that important ...
    /*
      swap(vecHsq,vecH);
      vecHsq.QR_factor();
      vecH.init_diag(cnt,1.);
      vecH.enlarge_below(dim-cnt,0.);
      vecHsq.Q_times(vecH,cnt);
    */
    
  }

  //========= low rank entries ===========
  if (vecH.coldim()>0){

    //compute (the square root of) the diagonal for this basis
    vecHsq.inv();
    vecHsq%=vecHsq;          //  Sigma^2
    barQ.scale_cols(vecHsq); // barQ <- Q*Sigma = diffvecs'*P
    barQ.abs();
    vecHsq=maxrows(barQ);
     
    Real lmaxbound=relprec*max(vecHsq);
    Indexmatrix delind(vecHsq.dim(),1);
    delind.init(0,0,Integer(0));
    for (Integer i=0;i<vecHsq.dim();i++){
      if (vecHsq(i)<lmaxbound)
	delind.concat_below(i);
    }
    if (delind.dim()){
      vecHsq.delete_cols(delind,true);
      vecH.delete_cols(delind,true);
    }
   
    vecH.scale_cols(vecHsq);

    if (cb_out(2)){
      get_out()<<"\n modellam0.dim="<<vecHsq.dim();
      if (vecHsq.dim()>0)
	get_out()<<" max(lam0)="<<max(vecHsq)<<" sum(lam0)="<<sum(vecHsq)<<" avg(vecHSq)="<<sum(vecHsq)/vecHsq.dim();
      get_out()<<std::endl;
    }
  }

  return 0;
}

//*****************************************************************************
//                    VariableMetricSVDSelection::add_variable_metric
//*****************************************************************************

int VariableMetricSVDSelection::vecH_by_aggregate(Matrix& vecH,
				      Real& relprec,
				      const MinorantBundle& bundle,
				      const MinorantPointer& aggregate,
				      Real function_factor,
				      Integer y_id,
				      const Matrix& y,
				      Real weightu,
				      Real model_maxviol)
{
  assert(function_factor>0.);
  assert(model_maxviol>0.);

  Integer dim=y.dim();
  vecH.init(dim,0,0.);
  relprec=1.;

  Real aggrval=aggregate.evaluate(y_id,y)/function_factor;
  Real aggrsqu=aggregate.dual_norm_squared()/function_factor/function_factor;
  Real aggrnorm=std::sqrt(aggrsqu);

  Integer bsz=Integer(bundle.size());
  Matrix val(bsz,1); chk_set_init(val,1);
  Real maxval=min_Real;
  for (Integer i=0;i<bsz;i++){ 
    if (bundle[unsigned(i)].empty()||bundle[unsigned(i)].aggregate())
      continue;
    Real v=bundle[unsigned(i)].evaluate(y_id,y);
    val(i)=v;
    if (v>maxval)
      maxval=v;
  }
  if (maxval==min_Real){
    return 0;
  }

  assert(model_maxviol>0.);  
  assert(maxval>=aggrval-1e10*(1+std::fabs(aggrval)));
  assert(weightu>0.);
  Real violation_eps=.9*model_maxviol/function_factor;
  Real maxmult=1e6*weightu/function_factor;
  Real mindiff=.1*violation_eps;
  Real normbnd=1e-4*max(1.,aggrsqu);
  relprec=min(1e-3,max(1e-12,mindiff)/(1.+CH_Matrix_Classes::abs(maxval)));

  //---- for each minorant compute the multiplier; for the diagonal take 
  //---- for each coordinate direction the maximum of each resulting diagonal 
  Matrix multval(bsz,1);
  multval.init(0,0,0.);
  Matrix aggrip(bsz,1);
  aggrip.init(0,0,0.);
  Indexmatrix bundleind(bsz,1);  //store the index within old_minorants
  bundleind.init(0,1,0.);

  Integer too_good=0;
  Integer too_bad=0;
  Integer too_close=0;

  for (Integer i=0;i<bsz;i++){ //start from the most recent one
    if (bundle[unsigned(i)].empty()||bundle[unsigned(i)].aggregate())
      continue;
    Real gnormsqu=bundle[unsigned(i)].dual_norm_squared();    
    if (val(i)+aggrnorm*std::sqrt(gnormsqu)/weightu<aggrval+violation_eps){
      too_bad++;
      continue;
    }

    //Real dv=maxval-val(i);
    Real dv=violation_eps+aggrval-val(i);
    if (dv<mindiff){
      too_good++;
      dv=mindiff;
    }
    Real ipval=bundle[unsigned(i)].ip(aggregate)/function_factor;   
    Real diffnormsqu=aggrsqu-2*ipval+gnormsqu;
    if (diffnormsqu<normbnd){
      too_close++;
      continue;
    }
    bundleind.concat_below(i);
    multval.concat_below(std::sqrt(.5/dv));
    aggrip.concat_below(ipval);
  }
  assert(multval.dim()==bundleind.dim());

 
  if (cb_out(2)){
    get_out()<<"\n bsz="<<bsz<<" relprec="<<relprec<<" mindiff="<<mindiff<<" maxmult="<<maxmult<<" tg="<<too_good<<" tb="<<too_bad<<" tc="<<too_close;
    if (multval.dim()>0)
      get_out()<<" avg(sqrt(multval))="<<sum(multval)/multval.dim();
  }

  //compute (approximate) SVD basis for the difference of minorants and aggr
  Symmatrix S;

  Matrix vecHtdiff;
  Matrix vecHsq;  
  if (multval.dim()>0){ 

    S.newsize(multval.dim()); chk_set_init(S,1);
    for (Integer i=0;i<S.rowdim();i++){
      Real ipgi=aggrip(i);
      Real mvi=multval(i);
      const MinorantPointer& mpi=bundle[unsigned(bundleind(i))];
      for(Integer j=i;j<S.rowdim();j++){
	Real ipij=mpi.ip(bundle[unsigned(bundleind(j))]);
	S(i,j)=(aggrsqu-ipgi-aggrip(j)+ipij)*mvi*multval(j);
      }
    }
    
    int retval=S.eig(vecHtdiff,vecHsq,false);
    if ((min(vecHsq)<-1e-6)||(retval)){
      if (cb_out())
	get_out()<<"**** WARNING: VariableMetricSVDSelection::add_variable_metric(...): S.eig failed and returned "<<retval<<" (order="<<S.rowdim()<<")"<<std::endl;
    }
    Real lmaxbound=relprec*vecHsq(0);
    Integer cnt=1;
    Integer maxdir=vecHsq.dim();//   min(vecHsq.dim(),20);
    while ((cnt<maxdir)&&(vecHsq(cnt)>lmaxbound))
      cnt++;
    if (cnt<vecHsq.dim()){
      vecHsq.reduce_length(cnt);
      vecHtdiff.delete_cols(Range(cnt,vecHtdiff.coldim()-1));
    }
    vecHsq.sqrt();
    vecHsq.inv();
    vecHtdiff.scale_cols(vecHsq);

    vecH.newsize(dim,vecHtdiff.coldim()); chk_set_init(vecH,1);
    vecHtdiff.scale_rows(multval);
    vecHsq.init(sumrows(vecHtdiff),-1./function_factor);
    aggregate.left_genmult(vecHsq,vecH);
    for(Integer i=0;i<bundleind.dim();i++)
      bundle[unsigned(bundleind(i))].left_genmult(vecHtdiff,vecH,1.,1.,0,0,i);

    //orthogonalize vecH, maybe not that important ...
    /*
      swap(vecHsq,vecH);
      vecHsq.QR_factor();
      vecH.init_diag(cnt,1.);
      vecH.enlarge_below(dim-cnt,0.);
      vecHsq.Q_times(vecH,cnt);
    */
    
  }

  //========= low rank entries ===========
  if (vecH.dim()>0){

    //compute (the square root of) the diagonal for this basis
    vecHsq.newsize(vecH.coldim(),1); chk_set_init(vecHsq,1);
    Matrix aggrip(vecH.coldim(),1); chk_set_init(aggrip,1);
    for (Integer j=0;j<vecH.coldim();j++){
      vecHtdiff.init(dim,1,vecH.get_store()+j*dim);
      Real aggrval=aggregate.ip(vecHtdiff)/function_factor;
      aggrip(j)=aggrval;
      Real maxval=multval(0)*CH_Matrix_Classes::abs(bundle[unsigned(bundleind(0))].ip(vecHtdiff)-aggrval);
      for (Integer i=1;i<bundleind.dim();i++){
	maxval=max(maxval,multval(i)*CH_Matrix_Classes::abs(bundle[unsigned(bundleind(i))].ip(vecHtdiff)-aggrval));
      }
      vecHsq(j)=maxval;
    }    
 
    //remove columns with small diagonal entries
    /*
    Indexmatrix sind;
    sortindex(vecHsq,sind,false);
    Real dnorm=aggrsqu/weightu;
    Real dnormbound=(weightu>.5)?dnorm/2.:aggrsqu/2.;
    Real lmaxbound=relprec*vecHsq(sind(0));
    Integer cnt=0;
    while(cnt<vecHsq.dim()){
      Integer ind=sind(cnt);
      Real d=vecHsq(ind);
      if (d<lmaxbound)
        break;
      Real a=sqr(aggrip(ind))*d*weightu/(d+weightu);
      if (dnorm-a<dnormbound){
	if (cb_out(2)){
	  get_out()<<"\n aggrsqu/u="<<aggrsqu/weightu<<" dnorm("<<cnt-1<<")="<<dnorm<<" dnormbound="<<dnormbound<<" val("<<cnt<<")="<<a<<std::endl;
	}
	break;
      }
      dnorm-=a;
      cnt++;
    }
    if (cnt<vecHsq.dim()){
      vecH.delete_cols(sind(Range(cnt,sind.dim()-1)));
      vecHsq.delete_rows(sind(Range(cnt,sind.dim()-1)));
    }
    */  
    
    Real lmaxbound=relprec*max(vecHsq);
    Indexmatrix delind(vecHsq.dim(),1);
    delind.init(0,0,Integer(0));
    for (Integer i=0;i<vecHsq.dim();i++){
      if (vecHsq(i)<lmaxbound)
	delind.concat_below(i);
    }
    if (delind.dim()){
      vecHsq.delete_rows(delind,true);
      vecH.delete_cols(delind,true);
    }
   
    vecH.scale_cols(vecHsq);

    if (cb_out(2)){
      get_out()<<"\n aggrlam0.dim="<<vecHsq.dim();
      if (vecHsq.dim()>0)
	get_out()<<" max(lam0)="<<max(vecHsq)<<" sum(lam0)="<<sum(vecHsq)<<" avg(vecHSq)="<<sum(vecHsq)/vecHsq.dim();
      get_out()<<std::endl;
    }
  }
  
  return 0;

}

//*****************************************************************************
//                    VariableMetricSVDSelection::add_variable_metric
//*****************************************************************************

  int VariableMetricSVDSelection::add_variable_metric(VariableMetric& H,
						      Integer y_id,
						      const Matrix& y,
						      bool /* descent_step */,
						      Real weightu,
						      Real model_maxviol,
						      const Indexmatrix* /* indices */,
						      VariableMetricBundleData& bdata
						    )
{
  assert(weightu>0.);
  assert(model_maxviol>0.);

  if (n_latest_minorants<=0)
    return 0;


  int retval=0;

  Matrix& old_diagH=bdata.set_diagH();
  Matrix& old_lowrankH=bdata.set_lowrankH();
  Symmatrix& old_symH=bdata.set_denseH();
  const MinorantPointer& aggregate=bdata.get_aggregate();
  Real function_factor=bdata.get_function_factor();
  
  MinorantBundle minorants;
  if(bdata.get_latest_minorants(minorants,n_latest_minorants)){
    if (cb_out())
      get_out()<<"**** WARNING: VariableMetricSVDSelection::add_variable_metric(...): bdata.get_latest_minorants(..) failed"<<std::endl;
    retval++;
  }   
  
  if (aggregate.empty()||(minorants.size()==0)){
    old_diagH.init(0,0,0.);
    old_lowrankH.init(0,0,0.);
    old_symH.init(0,0.);
    return retval;
  }

  old_symH.init(0,0.); //currently not supported by this routine

  bool use_dense=H.supports_dense_variable_metric();
  bool use_lowrank=(use_dense|| H.supports_lowrank_variable_metric());
  bool use_diagonal=((!use_lowrank) && H.supports_diagonal_variable_metric());

  if (oldfactor>=1.){
    //do not modify the current metric but only use the old one
    if ((use_diagonal)||(use_lowrank)){
      assert((old_lowrankH.coldim()==0)||(old_lowrankH.rowdim()==y.rowdim()));
      assert((old_diagH.rowdim()==0)||(old_diagH.rowdim()==y.rowdim()));
      if (function_factor!=1.){
	Matrix diagH(old_diagH,function_factor);
	Matrix vecH(old_lowrankH,function_factor);
	if (H.add_variable_metric(diagH,vecH)){
	  if (cb_out())
	    get_out()<<"**** WARNING: VariableMetricSVDSelection::add_variable_metric(...): H.add_variable_metric(..) failed"<<std::endl;
	  retval++;
	}
      }
      else {
	if (H.add_variable_metric(old_diagH,old_lowrankH)){
	  if (cb_out())
	    get_out()<<"**** WARNING: VariableMetricSVDSelection::add_variable_metric(...): H.add_variable_metric(..) failed"<<std::endl;
	  retval++;
	}
      }
    }
    return retval;
  }
    
  //========= collect the required data out of bdata

  if (Integer(minorants.size())>n_latest_minorants)
    minorants.resize(unsigned(n_latest_minorants));
  
  MinorantBundle model;
  Matrix modelcoeff;
  if(bdata.get_model_data(model,modelcoeff)){
    if (cb_out())
      get_out()<<"**** WARNING: VariableMetricSVDSelection::add_variable_metric(...): bdata.get_model_minorants(..) failed"<<std::endl;
    retval++;
  }   
 
		         
  //=========== prepare for and call the selected method
  
  Matrix vecH;    //will contain the resulting low rank metric 
  Real relprec=0.;  //determined within and returned by the methods
  if (model.size()==0){
    retval=vecH_by_aggregate(vecH,relprec,minorants,aggregate,function_factor,
			     y_id,y,weightu,model_maxviol);
    if ((retval)&&(cb_out())){
      get_out()<<"**** ERROR: VariableMetricSVDSelection::add_variable_metric(): vecH_by_aggregate() returned "<<retval<<std::endl;
    }
  }
  else {
    switch(selection_method){
    case 0: default: {
      retval=vecH_by_aggregate(vecH,relprec,minorants,aggregate,function_factor,
			       y_id,y,weightu,model_maxviol);
      if ((retval)&&(cb_out())){
	get_out()<<"**** ERROR: VariableMetricSVDSelection::add_variable_metric(): vecH_by_aggregate() returned "<<retval<<std::endl;
      }
    } break;
    case 1: {
      retval=vecH_by_model(vecH,relprec,minorants,model,aggregate,function_factor,y_id,y,weightu,model_maxviol);
      if ((retval)&&(cb_out())){
	get_out()<<"**** ERROR: VariableMetricSVDSelection::add_variable_metric(): vecH_by_model() returned "<<retval<<std::endl;
      }
    } break;
    case 2: {
      retval=vecH_orthogonal_to_model(vecH,relprec,minorants,model,modelcoeff,aggregate,function_factor,y_id,y,weightu,model_maxviol);
      if ((retval)&&(cb_out())){
	get_out()<<"**** ERROR: VariableMetricSVDSelection::add_variable_metric(): vecH_orthongonal_to_model() returned "<<retval<<std::endl;
      }
    } break;
    case 3: {
      retval=vecH_orthogonal_to_modelSVD(vecH,relprec,minorants,model,modelcoeff,aggregate,function_factor,y_id,y,weightu,model_maxviol);
      if ((retval)&&(cb_out())){
	get_out()<<"**** ERROR: VariableMetricSVDSelection::add_variable_metric(): vecH_orthongonal_to_modelSVD() returned "<<retval<<std::endl;
      }
    } break;
    case 4: {
      retval=vecH_weighted_SVD(vecH,relprec,minorants,model,modelcoeff,aggregate,function_factor,y_id,y,weightu,model_maxviol);
      if ((retval)&&(cb_out())){
	get_out()<<"**** ERROR: VariableMetricSVDSelection::add_variable_metric(): vecH_weighted_SVD() returned "<<retval<<std::endl;
      }
    } break;
    case 5: {
      retval=vecH_weighted_SVDs(vecH,relprec,minorants,model,modelcoeff,aggregate,function_factor,y_id,y,weightu,model_maxviol);
      if ((retval)&&(cb_out())){
	get_out()<<"**** ERROR: VariableMetricSVDSelection::add_variable_metric(): vecH_weighted_SVDs() returned "<<retval<<std::endl;
      }
    } break;
    }
  }
      
  Integer dim=y.dim();
  
  //when precision requirements get high, form a kind of average with the previous low rank representations
  if ((oldfactor>0.)&&(relprec<1e-4)&&(old_lowrankH.rowdim()==vecH.rowdim())&&(old_lowrankH.coldim()>0)){
    if (vecH.coldim()==0){
      vecH=old_lowrankH;
    }
    else {
      assert(oldfactor<1.);
      vecH*=std::sqrt(1.-oldfactor);
      old_lowrankH *= std::sqrt(oldfactor);
      old_lowrankH.concat_right(vecH);
      Symmatrix S;
      rankadd(old_lowrankH,S,1.,0.,1);
      Matrix P,lam;
      int retval=S.eig(P,lam,false);
      if ((min(lam)<-1e-6)||(retval)){
	if (cb_out())
	  get_out()<<"**** WARNING: VariableMetricSVDSelection::add_variable_metric(...): S.eig failed for average with previous lowrank representation and returned "<<retval<<" (order="<<S.rowdim()<<") and minimum eigenvalue = "<<min(lam)<<std::endl;
      }
      Real lmaxbound=relprec*lam(0);
      Integer cnt=1;
      while ((cnt<lam.dim())&&(lam(cnt)>lmaxbound))
	cnt++;
      if (cnt<lam.dim()){
	lam.reduce_length(cnt);
	P.delete_cols(Range(cnt,P.coldim()-1));
      }
      genmult(old_lowrankH,P,vecH);
      
      if (cb_out(2)){
	get_out()<<"\n lam1.dim="<<lam.dim()<<" max(lam1)="<<max(lam)<<" sum(lam1)="<<sum(lam)<<" avg(lam1)="<<sum(lam)/lam.dim()<<std::endl;
      }
    }
    
    assert(vecH.rowdim()==dim);
  }

  //======== store determined low rank metric in BundleData ===========
  old_lowrankH=vecH;
  old_diagH.init(0,0,0.);  //only lowrank information is generated here

  //======== form diagonal entries if only these are accepted by H ===========
  Matrix diagH;
  if ((!use_lowrank)&&(use_diagonal)&&(vecH.coldim()>0)){
    const Real* rowp=vecH.get_store();
    Integer coldim=vecH.coldim();
    diagH.newsize(dim,1); chk_set_init(diagH,1);
    for(Integer i=0;i<dim;i++){
      diagH(i)=mat_ip(coldim,rowp,dim,rowp,dim);
      rowp++;
    }
    vecH.init(0,0,0.);

    if (cb_out(2)){
      get_out()<<" max(diagH0)="<<max(diagH)<<" sum(diagH0)="<<sum(diagH);
      get_out()<<std::endl;
    }

  }

  //=========== add the new metric contribution

  if ((retval==0)&&((use_diagonal)||(use_lowrank))){
    if (function_factor!=1.){
      diagH*=function_factor;
      vecH*=function_factor;
    }
    if (H.add_variable_metric(diagH,vecH)){
      if (cb_out())
	get_out()<<"**** WARNING: VariableMetricSVDSelection::add_variable_metric(...): H.add_variable_metric(..) failed"<<std::endl;
      retval++;
    }
  }

  return retval;

}


}

