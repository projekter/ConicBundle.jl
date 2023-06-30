/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/BundleRQBWeight.cxx
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



#include "BundleRQBWeight.hxx"
#include "mymath.hxx"

 
using namespace CH_Matrix_Classes;

namespace ConicBundle {

  BundleRQBWeight::BundleRQBWeight(BundleWeight* bwp,
				   const CBout* cbo,int incr):
    BundleWeight(cbo,incr)
{
  clear();
  if (bwp){
    maxweight=bwp->get_maxweight();
    minweight=bwp->get_minweight();
    weight=bwp->get_weight();
    next_weight_set=bwp->get_next_weight_set();
    weightchanged=bwp->weight_changed();
  }
}

  BundleRQBWeight::BundleRQBWeight(Real im1,
				   Real im2,
				   Real im3,
				   Real ieta,
				   BundleWeight* bwp,
				   const CBout* cbo,int incr):
    BundleWeight(cbo,incr)
{
  clear();
  m1=min(max(im1,10.*eps_Real),1-10.*eps_Real);
  m2=min(max(im2,m1+(1-m1)*1e-6),1.);
  m3=max(10.*eps_Real,im3);
  eta=max(10.*eps_Real,ieta);
  if (bwp){
    maxweight=bwp->get_maxweight();
    minweight=bwp->get_minweight();
    weight=bwp->get_weight();
    next_weight_set=bwp->get_next_weight_set();
    weightchanged=bwp->weight_changed();
  }
}

void BundleRQBWeight::set_defaults()
{
  tL=0.;
  tR=0.;
  t=1.;
  m1=.1;
  m2=.2;
  m3=1.;
  eta=1e-6;
  frel=1e-5;
  minweight=-1;
  maxweight=-1;
} 

void BundleRQBWeight::clear()
{
  groundset=0;
  model=0;
  center_aggr.clear();
  weight=-1.;
  weightchanged=false; 
  next_weight_set=false;
  set_defaults();
}
  
  int BundleRQBWeight::init(Real norm2subg,Groundset* gs,BundleModel* mo)
{
  tL=0.;
  tR=0.;
  t=1.;
  groundset=gs;
  model=mo;
  assert(gs);
  center_aggr=groundset->get_gs_aggregate();
  if (model){
    Integer dummy;
    model->transform()->get_model_aggregate(dummy,center_aggr);
  }
  if (!next_weight_set){
    if (weight<0.) {
      Integer dim=groundset->get_dim();
      if (dim==0)
	weight=1.;
      else {
	Real d=sqrt(norm2subg);
	if (d<dim*1e-10)
	  weight=1.;
	else
	  weight=max(d,1e-4);
      }
    }
  }
  if (minweight<=0) 
    minweight=max(1e-10*weight,10.*eps_Real);
  weight=max(minweight,weight);
  if (maxweight>0) 
    weight=min(weight,maxweight);  
  weightchanged=true;
 return 0;
}

Real BundleRQBWeight::get_weight() const
{ return weight; }
    
bool BundleRQBWeight::weight_changed() const
{ return weightchanged; }

int BundleRQBWeight::delta_subg(Matrix& delta_subg,MinorantPointer& center,MinorantPointer& cand)
{
  if ((center.empty())||(cand.empty()))
    return 1;
  assert(groundset);
  delta_subg.init(groundset->get_dim(),1,0.);
  Real dummy;
  cand.get_minorant(dummy,delta_subg,0);
  center.get_minorant(dummy,delta_subg,0,-1.,true);
  return 0;
}

int BundleRQBWeight::descent_update(Real /* newval */,
				   Real oldval,
				   Real modelval,
				   const Matrix& y, 
				   const Matrix& newy,
				    Real norm2subg,
				    BundleProxObject* /* Hp */)
{
 if (weight<0) {
   if (cb_out())
     get_out()<<"**** ERROR BundleRQBWeight::descent_update(.......): negative weight "<<weight<<std::flush;
   return 1;
  }

 if (cb_out(1))
   get_out()<<"  descent step"<<std::flush;

 Real oldweight=weight;
 weightchanged=false;
 next_weight_set=false;
 tL=t;

 assert(groundset);
 Matrix dy=newy;
 dy-=y;
 Real delhat=oldval-modelval;
 if (t==1.){
   delhat1=delhat;
 }
 MinorantPointer cand_mnrt=groundset->get_gs_aggregate();
 if (model){
   Integer dummy;
   model->transform()->get_function_minorant(dummy,cand_mnrt);
 }
 if (cand_mnrt.empty()){
   if (cb_out(1))
     get_out()<<", no candidate minorant "<<std::endl;
   tL=0.;
   tR=0.;
   t=1.;
   return 1;
 }
 MinorantPointer cand_aggr=groundset->get_gs_aggregate();
 if (model){
   Integer dummy;
   model->transform()->get_model_aggregate(dummy,cand_aggr);
 }
 if (cand_aggr.empty()){
   if (cb_out(1))
     get_out()<<", no candidate aggregate "<<std::flush;
 }

 if (cand_mnrt.ip(dy)>=-m2*delhat) {  // quasi Newton udpate 

   MinorantPointer center_mnrt=groundset->get_gs_aggregate();
   if (model){
     Integer dummy;
     model->transform()->get_center_minorant(dummy,center_mnrt);
   }

   Integer best_index=0;
   Matrix dsubg;
   Real new_weight=weight;
   if (delta_subg(dsubg,center_aggr,cand_aggr)==0){
     Real ipsy=ip(dsubg,dy);
     if (ipsy>0){
       Real w=1./(ipsy/ip(dsubg,dsubg)+1./weight);
       if (w<new_weight){
	 new_weight=w;
	 best_index=1;
       }
     }
   }
   if (delta_subg(dsubg,center_aggr,cand_mnrt)==0){
     Real ipsy=ip(dsubg,dy);
     if (ipsy>0){
       Real w=1./(ipsy/ip(dsubg,dsubg)+1./weight);
       if (w<new_weight){
	 new_weight=w;
	 best_index=2;
       }
     }
   }
   if (delta_subg(dsubg,center_mnrt,cand_mnrt)==0){
     Real ipsy=ip(dsubg,dy);
     if (ipsy>0){
       Real w=1./(ipsy/ip(dsubg,dsubg)+1./weight);
       if (w<new_weight){
	 new_weight=w;
	 best_index=3;
       }
     }
   }
   if (delta_subg(dsubg,center_mnrt,cand_aggr)==0){
     Real ipsy=ip(dsubg,dy);
     if (ipsy>0){
       Real w=1./(ipsy/ip(dsubg,dsubg)+1./weight);
       if (w<new_weight){
	 new_weight=w;
	 best_index=4;
       }
     }
   }

   weight=new_weight;
   tL=0.;
   tR=0.;
   t=1.;
   if (cb_out(1))
     get_out()<<", rqb "<<weight<<"("<<best_index<<")";
 }

 else {  //no quasi Newton update

   Real linerr=oldval-cand_aggr.evaluate(-1,y);

   if ((tR>0.)||
       ((std::sqrt(norm2subg)>eta)&&
	(delhat-linerr>1e-8*linerr)&&
	(norm2subg*t*t>=frel*delhat))){ 
     if (tR==0.) {
       t*=8.;
       weight/=8.;
     
       if (cb_out(1))
	 get_out()<<", extrapolate ["<<tL<<","<<tR<<"] "<<t;
     }
     else {
       weight*=t;
       t=.5*(tL+tR);
       weight/=t;

       if (cb_out(1))
	 get_out()<<", interpolate ["<<tL<<","<<tR<<"] "<<t;
     }
   }
   else {
     if (cb_out(1))
       get_out()<<", cuttingplane "<<weight;
     tL=0.;
     tR=0.;
     t=1.;
   }
 }

 center_aggr=cand_aggr;
 
 assert((minweight<0)||(maxweight<0)||(minweight<=maxweight));

 if (minweight>0.)
   weight=max(weight,minweight);
 if (maxweight>0.)
   weight=min(weight,maxweight);
 
 if (std::fabs(weight-oldweight)>1e-10*weight) {
     weightchanged=true;
 }

 if (cb_out(1))
    get_out()<<", unew="<<weight<<std::endl;

 return 0;
}


int BundleRQBWeight::nullstep_update(
				     Real /* newval */,
				    Real oldval,
				    Real modelval,
				    const Matrix&  y , 
				     const Matrix&  /* newy */, 
				     MinorantPointer& /* new_minorant */,
				     MinorantPointer& /* aggregate */,
				     Real /* nullstep_bound */,
				     Real /* normsubg2 */,
				     BundleProxObject* /* Hp */)
{
  if (weight<0) {
   if (cb_out())
     get_out()<<"**** ERROR BundleRQBWeight::descent_update(.......): negative weight "<<weight<<std::flush;
   return 1;
  }

 if (cb_out(1))
   get_out()<<"  null step"<<std::flush;
 assert(groundset);

 Real oldweight=weight;
 weightchanged=false;
 next_weight_set=false;
 tR=t;

 Real delhat=oldval-modelval;
 if (t==1.)
   delhat1=oldval-modelval;
 MinorantPointer cand_mnrt=groundset->get_gs_aggregate();
 if (model){
   Integer dummy;
   model->transform()->get_function_minorant(dummy,cand_mnrt);
 }
 if (cand_mnrt.empty()){
   if (cb_out(1))
     get_out()<<", no candidate minorant"<<std::endl;
   tL=0.;
   tR=0.;
   t=1.;
   return 1;
 }
 Real linerr=oldval-cand_mnrt.evaluate(-1,y);

 if ((tL>0.)||((linerr > m3*delhat1)&&(delhat>frel*std::abs(oldval)))) {
   weight*=t;
   t=.5*(tL+tR);
   weight/=t;
   if (cb_out(1))
     get_out()<<", interpolate ["<<tL<<","<<tR<<"] "<<t<<std::flush;
 }
 else {
   tL=0.;
   tR=0.;
   t=1.;
   if (cb_out(1))
     get_out()<<", unchanged "<<weight<<std::flush;
 }

 assert((minweight<0)||(maxweight<0)||(minweight<=maxweight));

 if (minweight>0.)
   weight=max(weight,minweight);
 if (maxweight>0.)
   weight=min(weight,maxweight);

 if (::fabs(weight-oldweight)>1e-10*weight) {
     weightchanged=true;
 }

 if (cb_out(1))
    get_out()<<", unew="<<weight<<std::endl;

 return 0;
}

  int BundleRQBWeight::apply_modification(const GroundsetModification& /* gsmdf */)
{
 weightchanged=true;
 tL=0.;
 tR=0.;
 t=1.;
 center_aggr.clear();
 return 0;  
}

}

