/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  Matrix/pcg.cxx
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
#include "pcg.hxx"

using namespace CH_Tools;

namespace CH_Matrix_Classes {

// *************************************************************************
//                              constructor(s)
// *************************************************************************

PCG::PCG(std::ostream* out,int pril):
  myout(out),print_level(pril)
{
  maxit=-1;
  resnorm=0.;
  avg_reduction=-1.;
  nmult=0;
  err=0;
}

// *************************************************************************
//                              compute
// *************************************************************************

//this closely follows the book by Elman, Silvester and Wathen, reprinted 2006

  int PCG::compute(IterativeSystemObject& sys,
		   Matrix& x,
		   Real in_termprec,
		   Matrix* storex,
		   Integer storestep
		   )
{
  Integer dim=sys.ItSys_rhs().dim();
  if ((x.dim()>0)&&(x.dim()!=dim)){
    if (myout)
      (*myout)<<"**** ERROR PCG::compute(...): dimensions of input vector and system do not match"<<std::endl;
    err=1;
    return err;
  }
  err=0;
  nmult=0;
  stall=0;
  termprec=in_termprec;
  
  if (x.dim()==0){
    x.init(dim,1,0.);
    Ap.init(dim,1,0.);
  }
  else { 
    nmult++;
    if ((err=sys.ItSys_mult(x,Ap))){
      if (myout)
	(*myout)<<"****ERROR PCG::compute(...): first call to system multiplication routine returned error "<<err<<std::endl;
      if ((myout)&&(print_level>=1))
	(*myout)<<"PCG "<<std::setw(2)<<nmult<<"("<<err<<")"<<std::endl;
      return err;
    }
  }
  r.init(sys.ItSys_rhs());
  r-=Ap;
  Real n2r=norm2(r);
  z.init(r);
  if ((err=sys.ItSys_precondM1(z))){
    if (myout)
      (*myout)<<"****ERROR PCG::compute(...): first call to preconditioning routine M1 returned error "<<err<<std::endl;
    if ((myout)&&(print_level>=1))
      (*myout)<<"PCG "<<std::setw(2)<<nmult<<"("<<err<<")"<<std::endl;
    return err;
  }
  Real rho=ip(r,z);
  Integer pc_notpd=0;  //counts if preconditioner seems to not be positive definite
  if (rho<=1e-10*n2r){
    pc_notpd++;
  }
  p.init(z);

  while((rho>eps_Real*n2r)&&(stall<5*dim)&&((maxit<0)||(nmult<maxit))){
    if ((myout)&&(print_level>=2))
      (*myout)<<" PCG("<<nmult<<"): "<<resnorm<<std::endl;
    nmult++;
    if ((err=sys.ItSys_mult(p,Ap))){
      if (myout)
	(*myout)<<"****ERROR PCG::compute(...): system multiplication routine returned error "<<err<<std::endl;
      if ((myout)&&(print_level>=1))
	(*myout)<<"PCG "<<std::setw(2)<<nmult<<"("<<err<<"): "<<resnorm<<std::endl;
      return err;
    }
    Real alpha=rho/ip(Ap,p);
    x.xpeya(p,alpha);
    if (storex){
      if (storex->coldim()==0)
	storex->concat_right(x);
      else if((storestep>0)&&(nmult%storestep==0))
	storex->concat_right(x);
    }
    r.xpeya(Ap,-alpha);
    Real nextn2r=norm2(r);
    if (nextn2r>.999*n2r)
      stall++; 
    n2r=nextn2r;

    //for the avg_reduction we start once nmult==2, so its initial value is ignored
    old_resnorm=resnorm;
    resnorm=n2r;
    if ((nmult>=2)&&(nmult<11)){
      avg_reduction=(avg_reduction*(nmult-2)+resnorm/old_resnorm)/(nmult-1);
    }
    if (nmult>=11){
      avg_reduction=avg_reduction*9/10.+resnorm/old_resnorm/10.;
    }
    
    if (n2r<termprec)
      break;
    z=r;
    if ((err=sys.ItSys_precondM1(z))){
      if (myout)
	(*myout)<<"****ERROR PCG::compute(...): preconditioning routine M1 returned error "<<err<<std::endl;
      if ((myout)&&(print_level>=1))
	(*myout)<<"PCG "<<std::setw(2)<<nmult<<"("<<err<<"): "<<resnorm<<std::endl;
      return err;
    }
    Real nextrho=ip(r,z);
    if (nextrho<1e-10*n2r){
      pc_notpd++;
      if (pc_notpd>10){
	err=pc_notpd;
	if (myout)
	  (*myout)<<"****ERROR PCG::compute(...): preconditioning seems not to be positive definite, number of times this problem was encounterd ="<<err<<std::endl;
	if ((myout)&&(print_level>=1))
	  (*myout)<<"PCG "<<std::setw(2)<<nmult<<"("<<err<<"): "<<resnorm<<std::endl;
	return err;
      }
    }
    p*=(nextrho/rho);
    p+=z;
    rho=nextrho;
  }

  if ((myout)&&(print_level>=1))
    (*myout)<<" _PCG("<<nmult<<"): "<<resnorm<<std::endl;
  
  if ((resnorm>termprec)&&(myout)){
    (*myout)<<"**** WARNING: PCG::compute(): did not reach required precision; resnorm="<<resnorm<<" termprec="<<termprec<<" stallcnt="<<stall<<std::endl;
  }
  
  return err;
}

}
