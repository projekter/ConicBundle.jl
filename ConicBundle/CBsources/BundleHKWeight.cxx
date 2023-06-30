/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/BundleHKWeight.cxx
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



#include "BundleHKWeight.hxx"
#include "mymath.hxx"

 
using namespace CH_Matrix_Classes;

namespace ConicBundle {

  BundleHKWeight::BundleHKWeight(Real mRin,BundleWeight* bwp,
				 const CBout* cbo,int incr):
    BundleWeight(cbo,incr)
{
  nullstep_updates=0;
  mR=mRin;
  assert((.5<=mR)&&(mR<1.));
  clear();
  if (bwp){
    maxweight=bwp->get_maxweight();
    minweight=bwp->get_minweight();
    weight=bwp->get_weight();
    next_weight_set=bwp->get_next_weight_set();
    weightchanged=bwp->weight_changed();
  }
}

void BundleHKWeight::set_defaults()
{
  mR=0.5;
} 

void BundleHKWeight::clear()
{
  groundset=0;
  model=0;
  modelmax=CB_minus_infinity;
  minweight=maxweight=weight=-1.;
  weightchanged=false; 
  iweight=0;
  epsweight=1e30;
  next_weight_set=false;
}
  
int BundleHKWeight::init(Real norm2subg,Groundset* gs,BundleModel* mo)
{
 groundset=gs;
 model=mo;
 if (!next_weight_set){
   if (weight<0.) {
     assert(gs);
     Integer dim=gs->get_dim();
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
   minweight=max(1e-10*weight,1e-10);
 weight=max(minweight,weight);
 if (maxweight>0) 
   weight=min(weight,maxweight);  
 weightchanged=true;
 modelmax=CB_minus_infinity;
 epsweight=1e30;
 iweight=0;
 return 0;
}

Real BundleHKWeight::get_weight() const
{ return weight; }
    
bool BundleHKWeight::weight_changed() const
{ return weightchanged; }


int BundleHKWeight::descent_update(Real newval,
				   Real oldval,
				   Real modelval,
				   const Matrix& /* y */, 
				   const Matrix& /* newy */,
				   Real /* norm2subg */,
				   BundleProxObject* Hp)
{
  if (weight<0) {
   if (cb_out())
     get_out()<<"**** ERROR BundleHKWeight::descent_update(.......): negative weight "<<weight<<std::flush;
   return 1;
  }
 next_weight_set=false;
 valuelevels.init(0,1,0.);
 ratiolevels.init(0,1,0.);
     
 modelmax=max(modelmax,modelval);
 Real oldweight=weight;
 Real weightint=2*weight*(1.-(oldval-newval)/(oldval-modelval));

 //TEST
 /*
 assert(groundset);
 MinorantPointer newmnrt(groundset->get_gs_minorant());
 assert(model);
 Integer fid;
 if (model->get_function_minorant(fid,newmnrt)){
   if (cb_out(1))
     get_out()<<"**** ERROR in BundleHKWeight::nullstep_update(): model->get_fucntion_minorant failed"<<std::flush;
 }
 Real lin_approx=oldval-newmnrt.evaluate(-1,y);
 MinorantPointer aggr(groundset->get_gs_aggregate());
 Integer aid;
 if (model->get_model_aggregate(aid,aggr)){
   if (cb_out(1))
     get_out()<<"**** ERROR in BundleHKWeight::nullstep_update(): model->get_model_aggregate failed"<<std::flush;
 }
 Real kappa=.2;
 Real approx=lin_approx-kappa*(oldval-aggr.evaluate(-1,y));
 Matrix dy(y);
 dy-=newy;
 Real progr=kappa*aggr.ip(dy)-newmnrt.ip(dy);
 Real guess_factor=1.;
 if (approx>1e-6){
   guess_factor=progr/approx;
 }
 */

 if (cb_out(1))
   get_out()<<"  descent step, i_u="<<iweight<<std::flush;
 //<<" uint="<<weightint<<" ufactor="<<guess_factor<<"("<<approx<<","<<progr<<")"<<std::flush;
 if (((oldval-newval)>mR*(oldval-modelval))&&(iweight>0)){
   if (cb_out(1))
     get_out()<<" uint="<<weightint<<std::flush;
   weight=weightint;
 }
 else if (iweight>3){ //there were four, now is the 5th consecutive serious 
   if (cb_out(1))
     get_out()<<" i_u>3 u/2 "<<std::flush;
   weight/=2.;
 }
 else if (newval<modelmax){
   if (cb_out(1))
     get_out()<<" nv<linmax u/2 "<<std::flush;
   weight/=2.;
 }
 else if (Hp->get_short_QPsteps()>2) {
   if (cb_out(1))
     get_out()<<" shortQP("<<Hp->get_short_QPsteps()<<") u/2 "<<std::flush;
   weight/=2.;
 } 
 weight=max(oldweight/10.,weight);
 if (minweight>0) 
   weight=max(weight,minweight);
 if (cb_out(1))
    get_out()<<" unew="<<weight<<std::endl;
 epsweight=max(epsweight,2*(oldval-modelval));
 iweight=max(iweight+1,Integer(1));
 if (weight<oldweight) {
     weightchanged=true;
     iweight=1;
     modelmax=CB_minus_infinity;
 }
 else {
     weightchanged=false;
 }
 return 0;
}


int BundleHKWeight::nullstep_update(
				    Real newval,
				    Real oldval,
				    Real modelval,
				    const Matrix& y, 
				    const Matrix& /* newy */,
				    MinorantPointer& new_minorant,
				    MinorantPointer& aggregate,
				    Real nullstep_bound, 
				    Real normsubg2,
				    BundleProxObject* Hp)
{
  if (weight<0) {
   if (cb_out())
     get_out()<<"**** ERROR BundleHKWeight::nullstep_update(.......): negative weight "<<weight<<std::flush;
   return 1;
  }
 next_weight_set=false;
 
 Real oldweight=weight;
 switch(nullstep_updates){
 case 0: default: {
   Real weightint=2*weight*(1.-(oldval-newval)/(oldval-modelval));
   assert(groundset);
   MinorantPointer newmnrt(groundset->get_gs_minorant());
   assert(model);
   Integer fid;
   if (model->get_function_minorant(fid,newmnrt)){
     if (cb_out(1))
       get_out()<<"**** ERROR in BundleHKWeight::nullstep_update(): model->get_fucntion_minorant failed"<<std::flush;
   }
   Real lin_approx=oldval-newmnrt.evaluate(-1,y);
   //TEST
   /*
     MinorantPointer aggr(groundset->get_gs_aggregate());
     Integer aid;
     if (model->get_model_aggregate(aid,aggr)){
     if (cb_out(1))
     get_out()<<"**** ERROR in BundleHKWeight::nullstep_update(): model->get_model_aggregate failed"<<std::flush;
     }
     Real kappa=.2;
     Real approx=lin_approx-kappa*(oldval-aggr.evaluate(-1,y));
     Matrix dy(y);
     dy-=newy;
     Real progr=kappa*aggr.ip(dy)-newmnrt.ip(dy);
     Real guess_factor=1.;
     if (approx>1e-6){
     guess_factor=progr/approx;
     }
   */
   
   epsweight=min(epsweight,sqrt(normsubg2*weight)+oldval-modelval-normsubg2); 
   if (cb_out(1))
     get_out()<<"  null step, i_u="<<iweight<<" eps_v="<<epsweight<<" lapprox="<<lin_approx<<std::flush;
   //<<" uint="<<weightint<<" ufactor="<<guess_factor<<"("<<approx<<","<<progr<<")"<<std::flush;
   
   if ((lin_approx>max(epsweight,10.*(oldval-modelval)))&&(iweight<-3)){
     if (cb_out(1))
       get_out()<<" uint="<<weightint<<std::flush;
     weight=weightint;
   }
   weight=min(weight,10.*oldweight);   
   if (maxweight>0) 
     weight=min(weight,maxweight);
   break;
 }
 case 1: break;
 case 2: {
   Integer maxincr=2;
   Integer nincreases=0;
   do {
     if (nincreases>=valuelevels.rowdim()){
       valuelevels.concat_below(normsubg2);
       break;
     }
     else if (normsubg2<=valuelevels(nincreases)){
       valuelevels(nincreases)=normsubg2;
       break;
     }
   } while (++nincreases<maxincr);
   if (cb_out(1)){
     get_out()<<"  null step, i_u="<<iweight<<" ninc="<<nincreases<<std::flush;
   }
   if (nincreases>=maxincr){
     weight*=1.5;
     if (maxweight>0.)
       weight=min(weight,maxweight);
     valuelevels.init(0,1,0.);
   }
   break;
 }
 case 3: {
   Integer maxincr=2;
   Integer nincreases=0;
   Real ratio=-(oldval-newval)/(oldval-modelval);
   do {
     if (nincreases>=ratiolevels.rowdim()){
       ratiolevels.concat_below(ratio);
       break;
     }
     else if (ratio<=ratiolevels(nincreases)){
       ratiolevels(nincreases)=ratio;
       break;
     }
   } while (++nincreases<maxincr);
   if (cb_out(1)){
     get_out()<<"  null step, i_u="<<iweight<<" ratio="<<-ratio<<" ninc="<<nincreases<<std::flush;
   }
   if (nincreases>=maxincr){
     weight*=1.5;
     if (maxweight>0.)
       weight=min(weight,maxweight);
     ratiolevels.init(0,1,0.);
   }
   break;
 }
 case 4: {
   if (cb_out(1)){
     get_out()<<"  null step, i_u="<<iweight<<std::flush;
   }
   Integer maxincr=2;
   Integer nincreases=0;
   do {
     if (nincreases>=valuelevels.rowdim()){
       valuelevels.concat_below(normsubg2);
       break;
     }
     else if (normsubg2<=valuelevels(nincreases)){
       valuelevels(nincreases)=normsubg2;
       break;
     }
   } while (++nincreases<maxincr);
   if (cb_out(1)){
     get_out()<<" nsinc="<<nincreases<<std::flush;
   }
   if (nincreases>=maxincr){
     weight*=1.5;
     if (maxweight>0.)
       weight=min(weight,maxweight);
     valuelevels.init(0,1,0.);
     ratiolevels.init(0,1,0.);
     break;
   }
   Real ratio=-(oldval-newval)/(oldval-modelval);
   nincreases=0;
   do {
     if (nincreases>=ratiolevels.rowdim()){
       ratiolevels.concat_below(ratio);
       break;
     }
     else if (ratio<=ratiolevels(nincreases)){
       ratiolevels(nincreases)=ratio;
       break;
     }
   } while (++nincreases<maxincr);
   if (cb_out(1)){
     get_out()<<" ratio="<<-ratio<<" nrinc="<<nincreases<<std::flush;
   }
   if (nincreases>=maxincr){
     weight*=1.5;
     if (maxweight>0.)
       weight=min(weight,maxweight);
     ratiolevels.init(0,1,0.);
     valuelevels.init(0,1,0.);
     if (cb_out(1)){
       get_out()<<" unew="<<weight<<std::endl;
     }
   }
   break;
 }
 case 5: {
   Real maxviol=.9*(nullstep_bound-modelval);
   assert(maxviol>0.);
   Matrix diffvec(y.rowdim(),1,0.);
   Real diffoffset=0.;
   aggregate.get_minorant(diffoffset,diffvec,0);
   new_minorant.get_minorant(diffoffset,diffvec,0,-1.,true);
   Real diffval=diffoffset+ip(diffvec,y);
   if (cb_out(1)){
     get_out()<<"  null step, i_u="<<iweight<<" diffval="<<diffval;
   }
   if (diffval<0.){
     diffval=max(eps_Real,0.1*(oldval-new_minorant.evaluate(-1,y)));
     if (cb_out(1)){
       get_out()<<"("<<diffval<<")"<<std::flush;
     }
   }
   Real diffnorm=norm2(diffvec);
   diffvec/=diffnorm;
   Real weightinc=sqr(diffnorm)/(maxviol+diffval)/2.-Hp->norm_sqr(diffvec);
   if (cb_out(1)){
     get_out()<<" weightinc="<<weightinc<<std::flush;
   }
   if (weightinc>.5*weight){
     valuelevels.concat_below(weightinc);
   }
   if (valuelevels.rowdim()>=3){
     weight=min(10.*weight,weight+sum(valuelevels)/valuelevels.rowdim());
     if (maxweight>0.)
       weight=min(weight,maxweight);
     valuelevels.init(0,1,0.);
     if (cb_out(1)){
       get_out()<<" unew="<<weight<<std::endl;
     }
   }
   break;
 }
 } //end switch

 iweight=min(iweight-1,Integer(-1));
 if (weight>(1.+eps_Real)*oldweight) {
   if (cb_out(1)){
     get_out()<<" unew="<<weight<<std::endl;
   }
   iweight=-1;
   //modelmax=CB_minus_infinity;
   weightchanged=true;
 }
 else {
   weightchanged=false;
 }
 return 0;
}

  int BundleHKWeight::apply_modification(const GroundsetModification& /* gsmdf */)
{
 weightchanged=true;
 modelmax=CB_minus_infinity;
 epsweight=1e30;
 iweight=0;
 return 0;  
}


}

