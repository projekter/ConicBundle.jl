/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  Matrix/minres.cxx
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



#include <cmath>
#include <limits>
#include "mymath.hxx"
#include "minres.hxx"

 
using namespace CH_Tools;

namespace CH_Matrix_Classes {

// *************************************************************************
//                              constructor(s)
// *************************************************************************

MinRes::MinRes(std::ostream* out,int pril):
  myout(out),print_level(pril)
{
  maxit=-1;
  resnorm=0.;
  old_resnorm=0.;
  avg_reduction=-1.;
  nmult=0;
  err=0;

  stall=0;
}

// *************************************************************************
//                              compute
// *************************************************************************

//this closely follows the book by Elman, Silvester and Wathen, reprinted 2006

  int MinRes::compute(IterativeSystemObject& sys,
		      Matrix& x,
		      Real in_termprec,
		      Matrix* storex,
		      Integer storestep)
{
  Integer dim=sys.ItSys_rhs().dim();
  if ((x.dim()>0)&&(x.dim()!=dim)){
    if (myout)
      (*myout)<<"**** ERROR MinRes::compute(...): dimensions of input vector and system do not match"<<std::endl;
    err=1;
    return err;
  }
  if (x.dim()==0){
    x.init(dim,1,0.);
  }
  err=0;
  nmult=0;
  stall=0;
  termprec=in_termprec;
  
  vecvprev.init(dim,1,0.);  
  vecwprev.init(dim,1,0.);
  vecwcurr.init(dim,1,0.);
  
  vecvcurr.init(sys.ItSys_rhs());  
  vecAv.newsize(dim,1); chk_set_init(vecAv,1);
  nmult++;
  if ((err=sys.ItSys_mult(x,vecAv))){
    if (myout)
      (*myout)<<"****ERROR MinRes::compute(...): first call to system multiplication routine returned error "<<err<<std::endl;
    if ((myout)&&(print_level>=1))
      (*myout)<<"MinRes "<<std::setw(2)<<nmult<<"("<<err<<")"<<std::endl;
    return err;
  }
  vecvcurr-=vecAv;
  veczcurr.init(vecvcurr);
  if ((err=sys.ItSys_precondM1(veczcurr))){
    if (myout)
      (*myout)<<"****ERROR MinRes::compute(...): first call to preconditioning routine M1 returned error "<<err<<std::endl;
    if ((myout)&&(print_level>=1))
      (*myout)<<"MinRes "<<std::setw(2)<<nmult<<"("<<err<<")"<<std::endl;
    return err;
  }

  Real gammaprev=1.;
  Real gammacurr=std::sqrt(ip(veczcurr,vecvcurr));
  Real eta=gammacurr;
  resnorm=std::fabs(eta);
  old_resnorm=resnorm;
  avg_reduction=.1;
  Real sprev=0.;
  Real scurr=0.;
  Real cprev=1.;
  Real ccurr=1.;

  while((resnorm>termprec)&&(stall<=5*dim)&&((maxit<0)||(nmult<maxit))){
    if ((myout)&&((print_level>=2)||((print_level>=1)&&(nmult%100==0))))
      (*myout)<<" MinRes("<<nmult<<","<<termprec<<"): "<<resnorm<<std::endl;
    nmult++;
    veczprev.init(veczcurr,1./gammacurr);
    if ((err=sys.ItSys_mult(veczprev,vecAv))){
      if (myout)
	(*myout)<<"****ERROR MinRes::compute(...): system multiplication routine returned error "<<err<<std::endl;
      if ((myout)&&(print_level>=1))
	(*myout)<<"MinRes "<<std::setw(2)<<nmult<<"("<<err<<"): "<<resnorm<<std::endl;
      return err;
    }
    Real delta=ip(vecAv,veczprev);
    vecvprev*=-gammacurr/gammaprev;
    vecvprev.xpeya(vecvcurr,-delta/gammacurr);
    vecvprev+=vecAv;
    swap(vecvprev,vecvcurr);
    veczcurr.init(vecvcurr);
    if ((err=sys.ItSys_precondM1(veczcurr))){
      if (myout)
	(*myout)<<"****ERROR MinRes::compute(...): preconditioning routine M1 returned error "<<err<<std::endl;
      if ((myout)&&(print_level>=1))
	(*myout)<<"MinRes "<<std::setw(2)<<nmult<<"("<<err<<"): "<<resnorm<<std::endl;
      return err;
    }
    gammaprev=gammacurr;
    gammacurr=std::sqrt(ip(veczcurr,vecvcurr));
    Real alpha0=ccurr*delta-cprev*scurr*gammaprev;
    Real alpha1=std::sqrt(alpha0*alpha0+gammacurr*gammacurr);
    Real alpha2=scurr*delta+cprev*ccurr*gammaprev;
    Real alpha3=sprev*gammaprev;
    cprev=ccurr;
    ccurr=alpha0/alpha1;
    sprev=scurr;
    scurr=gammacurr/alpha1;
    vecwprev*=-alpha3/alpha1;
    vecwprev.xpeya(vecwcurr,-alpha2/alpha1);
    vecwprev.xpeya(veczprev,1./alpha1);
    swap(vecwprev,vecwcurr);
    x.xpeya(vecwcurr,ccurr*eta);
    if (storex){
      if (storex->coldim()==0)
	storex->concat_right(x);
      else if((storestep>0)&&(nmult%storestep==0))
	storex->concat_right(x);
    }
    
    if (std::fabs(scurr)>.999)
      stall++;
    eta*=-scurr;
    
    //for the avg_reduction we start once nmult==2, so its initial value is ignored
    old_resnorm=resnorm;
    resnorm=std::fabs(eta);
    if ((nmult>=2)&&(nmult<11)){
      avg_reduction=(avg_reduction*(nmult-2)+resnorm/old_resnorm)/(nmult-1);
    }
    if (nmult>=11){
      avg_reduction=avg_reduction*9/10.+resnorm/old_resnorm/10.;
    }
  }

  if ((myout)&&(print_level>=1))
    (*myout)<<" _MinRes("<<nmult<<"): "<<eta<<std::endl;
  
  if ((std::fabs(eta)>termprec)&&(myout)){
    (*myout)<<"**** WARNING: MinRes::compute(): MINRES did not reach required precision; eta="<<eta<<" termprec="<<termprec<<" nmult="<<nmult<<"(<"<<maxit<<") stallcnt="<<stall<<std::endl;
  }
  
  return err;
}

}
