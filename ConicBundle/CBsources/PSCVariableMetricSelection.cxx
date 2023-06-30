/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/PSCVariableMetricSelection.cxx
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



#include "PSCVariableMetricSelection.hxx"
#include "PSCData.hxx"

using namespace CH_Matrix_Classes;

namespace ConicBundle {


//*****************************************************************************
//                    PSCVariableMetricSelection::add_variable_metric
//*****************************************************************************

  int PSCVariableMetricSelection::add_variable_metric(VariableMetric& H,
						      Integer /* y_id */,
						      const Matrix& y,
						      bool  descent_step,
						      Real weightu,
						      Real /* model_maxviol */,
						      const Indexmatrix* /* indices */,
						      VariableMetricBundleData& bdata
						    )
{
  assert(weightu>0.);
  // assert(model_maxviol>0.);

  if (oracle==0){
    if (cb_out())
      get_out()<<"**** WARNING: PSCVariableMetricSelection::add_variable_metric(...): oracle pointer is NULL"<<std::endl;
    return 1;
  }

  PSCData* pscdata=dynamic_cast<PSCData*>(&bdata);
  if (pscdata==0){
    if (cb_out())
      get_out()<<"**** WARNING: PSCVariableMetricSelection::add_variable_metric(...): bdata could not be cast to PSCData"<<std::endl;
    return 1;
  }

  if (!descent_step)
    return 0;

  Integer dim=y.dim();
  Integer primdim=pscdata->get_keepsize();
  Integer skipstart=pscdata->get_activedim();

  Integer Ritzdim=pscdata->get_topvecs().coldim();
  assert((skipstart>0)&&(skipstart<=Ritzdim));

  if ((primdim==0)||(skipstart>=Ritzdim)){
    return 0;
  }
  
  //primalvecs is Q1 and next Q2 is stored in tmpat 
  tmpQ1.init(pscdata->get_primalvecs().rowdim(),primdim,pscdata->get_primalvecs().get_store());
  tmpQ2.init(pscdata->get_topvecs().rowdim(),Ritzdim-skipstart,pscdata->get_topvecs().get_store()+skipstart*pscdata->get_topvecs().rowdim());   //Q2
  Integer skipdim=tmpQ2.coldim();
  Real lmax=pscdata->get_Ritz_values()(0);

  Matrix lamH(primdim*skipdim,1,0.);
  for(Integer j=0;j<skipdim;j++){
    for(Integer i=0;i<primdim;i++){
      lamH(j*primdim+i)=2*pscdata->get_primaleigs()(i)/max(lmax-pscdata->get_Ritz_values()(skipstart+j),1e-6);     
    }
  }

  Matrix vecH(dim,primdim*skipdim,0.);
  lamH.sqrt();
  Matrix Q1AQ2;   //=Q_1^T*A*Q2
  for(Integer k=0;k<dim;k++){
    oracle->left_right_product(k,tmpQ1,tmpQ2,Q1AQ2);
    const Real* lamp=lamH.get_store();
    for(Integer j=0;j<skipdim;j++){
      for(Integer i=0;i<primdim;i++){
        vecH(k,j*primdim+i)=Q1AQ2(i,j)*(*lamp++);     
      }
    }
  }

  //compute diagonal representation
  Indexmatrix piv;
  Integer r=vecH.QR_factor(piv,min(0.1,1e-3*weightu));
  int retval=0;
  if (r==0) {
    vecH.init(dim,0,0.);
  }
  else {
    Real dmax=std::fabs(vecH(0,0));
    while (std::fabs(vecH(r-1,r-1))<1e-4*dmax){
      //std::cout<<" skip("<<r<<","<<vecH(r-1,r-1)<<")";
      r--;
    }
    //std::cout<<" dmax="<<dmax<<" dmin="<<vecH(r-1,r-1)<<" r="<<r<<std::endl;
    tmpQ2=vecH.rows(Range(0,r-1));
    tmpQ2.triu();
    //Matrix tmpvec;// =lamH(piv);
    //tmpvec.sqrt();
    //tmpQ2.scale_cols(tmpvec);
    Symmatrix S;
    rankadd(tmpQ2,S);
    dmax=1.;
    for (Integer i=0;i<S.rowdim();i++){
      dmax=max(dmax,S(i,i));
    }    
    S/=dmax;
    if ((retval=S.eig(tmpQ2,lamH,false))){
      if (cb_out())
	get_out()<<"**** ERROR: PSCVariableMetricSelection::add_variable_metric(): S.eig failed and returned "<<retval<<std::endl;
    }
    lamH*=dmax;
    Real eigmax=lamH(0);
    Integer neigs=0;
    while((neigs<min(lamH.dim(),max(primdim,dim/10)))&&
	  (lamH(neigs)>eigmax*mineigval_factor)){
      if (eigmax>maxeigval_factor*weightu)
	lamH(neigs)=std::sin(lamH(neigs)/eigmax*(3.1415926535/2.))*(maxeigval_factor*weightu);
      neigs++;
    }
    lamH.reduce_length(neigs);
    tmpQ2.delete_cols(Range(neigs,tmpQ2.coldim()-1));
    
    if (cb_out(2))
      get_out()<<" PSCVMS: neigs="<<neigs<<" max(lamH)="<<max(lamH)<<" min(lamH)="<<min(lamH)<<std::endl;
    //std::cout<<" PSCVMS: neigs="<<neigs<<" max(lamH)="<<max(lamH)<<" min(lamH)="<<min(lamH)<<std::endl;
    
    lamH.sqrt();
    tmpQ2.scale_cols(lamH);
    tmpQ2.enlarge_below(vecH.rowdim()-tmpQ2.rowdim(),0.);
    vecH.Q_times(tmpQ2,r);
    swap(vecH,tmpQ2);
    //tmpQ2 is free again for other use
  }
    
  //----  when precision requirements get high, form a kind of average with the previous low rank representations
  
  Real relprec=max((lmax-pscdata->get_Ritz_values()(skipstart))/lmax,0.);
  bool use_dense=H.supports_dense_variable_metric();
  bool use_lowrank=(use_dense|| H.supports_lowrank_variable_metric());
  bool use_diagonal=((!use_lowrank) && H.supports_diagonal_variable_metric());
  Matrix& old_diagH=bdata.set_diagH();
  Matrix& old_lowrankH=bdata.set_lowrankH();
  Symmatrix& old_symH=bdata.set_denseH();
  old_symH.init(0,0.); //currently not supported by this routine
  old_diagH.init(0,0,0.);//currently not supported by this routine
  
  if ((oldfactor>0.)&&(relprec<1e-4)&&(old_lowrankH.rowdim()==vecH.rowdim())&&(old_lowrankH.coldim()>0)){
    if (vecH.coldim()==0){
      vecH=old_lowrankH;
    }
    else {
      assert(oldfactor<1.);
      vecH*=std::sqrt(1.-oldfactor);
      old_lowrankH *= std::sqrt(oldfactor);
      old_lowrankH.concat_right(vecH);
      rankadd(old_lowrankH,S,1.,0.,1);
      Matrix P,lam;
      int retval=S.eig(tmpQ2,tmpvec,false);
      if ((min(tmpvec)<-1e-6)||(retval)){
	if (cb_out())
	  get_out()<<"**** WARNING: PSCVariableMetricSelection::add_variable_metric(...): S.eig failed for average with previous lowrank representation and returned "<<retval<<" (order="<<S.rowdim()<<") and minimum eigenvalue = "<<min(tmpvec)<<std::endl;
      }
      Real lmaxbound=relprec*tmpvec(0);
      Integer cnt=1;
      while ((cnt<tmpvec.dim())&&(tmpvec(cnt)>lmaxbound))
	cnt++;
      if (cnt<tmpvec.dim()){
	tmpvec.reduce_length(cnt);
	tmpQ2.delete_cols(Range(cnt,tmpQ2.coldim()-1));
      }
      genmult(old_lowrankH,tmpQ2,vecH);
      
      if (cb_out(2)){
	get_out()<<"\n lam1.dim="<<tmpvec.dim()<<" max(lam1)="<<tmpvec(lam)<<" sum(lam1)="<<sum(tmpvec)<<" avg(lam1)="<<sum(tmpvec)/tmpvec.dim()<<std::endl;
      }
    }
    
    assert(vecH.rowdim()==dim);
  }
  
  //======== store determined low rank metric in BundleData ===========
  old_lowrankH=vecH;
  old_diagH.init(0,0,0.);  //only lowrank information is generated here

  //======== form diagonal entries if only these are accepted by H ===========
  tmpvec.init(0,0,0.);
  if ((!use_lowrank)&&(use_diagonal)&&(vecH.coldim()>0)){
    const Real* rowp=vecH.get_store();
    Integer coldim=vecH.coldim();
    tmpvec.newsize(dim,1); chk_set_init(tmpvec,1);
    for(Integer i=0;i<dim;i++){
      tmpvec(i)=mat_ip(coldim,rowp,dim,rowp,dim);
      rowp++;
    }
    vecH.init(0,0,0.);

    if (cb_out(2)){
      get_out()<<" max(diagH0)="<<max(tmpvec)<<" sum(diagH0)="<<sum(tmpvec);
      get_out()<<std::endl;
    }

  }
  
  //=========== add the new metric contribution

  if ((retval==0)&&((use_diagonal)||(use_lowrank))){
    //--- function_factor is already included in the aggregate primal matrix
    if (H.add_variable_metric(tmpvec,vecH)){
      if (cb_out())
	get_out()<<"**** WARNING: PSCVariableMetricSelection::add_variable_metric(...): H.add_variable_metric(..) failed"<<std::endl;
      retval++;
    }
  }

  return retval;

}


}

