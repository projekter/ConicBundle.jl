/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/Coeffmat.cxx
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



#include <string.h>
#include <stdlib.h>
#include <cctype>

#include "Coeffmat.hxx"
#include "CMsymdense.hxx"
#include "CMsymsparse.hxx"
#include "CMgramdense.hxx"
#include "CMgramsparse.hxx"
#include "CMlowrankdd.hxx"
#include "CMlowranksd.hxx"
#include "CMlowrankss.hxx"
#include "CMsingleton.hxx"
#include "CMgramsparse_withoutdiag.hxx"

using namespace CH_Matrix_Classes;

namespace ConicBundle {

Coeffmat* coeffmat_read(std::istream& in)
{
 char name[80];
 in>>std::ws;
 if (!in) {
     if (materrout) 
       (*materrout)<<"*** ERROR: coeffmat_read(): input stream broken"<<std::endl;
     return 0;
 }
 int cnt=0;
 while((in)&&(!isspace(in.peek()))&&(cnt<80)){
   in.get(name[cnt++]);
 }
 if ((!in)||(cnt==80)){
     if (materrout) 
       (*materrout)<<"*** ERROR: coeffmat_read(): failed in reading name of constraint"<<std::endl;
     return 0;
 }
 name[cnt]=0;
 Coeffmat* p=0;
 if (!strcmp(name,"SYMMETRIC_DENSE")) p=new CMsymdense(in);
 else if (!strcmp(name,"SYMMETRIC_SPARSE")) p=new CMsymsparse(in);
 else if (!strcmp(name,"GRAM_DENSE")) p=new CMgramdense(in);
 else if (!strcmp(name,"GRAM_SPARSE")) p=new CMgramsparse(in);
 else if (!strcmp(name,"LOWRANK_DENSE_DENSE")) p=new CMlowrankdd(in);
 else if (!strcmp(name,"LOWRANK_SPARSE_DENSE")) p=new CMlowranksd(in);
 else if (!strcmp(name,"LOWRANK_SPARSE_SPARSE")) p=new CMlowrankss(in);
 else if (!strcmp(name,"SINGLETON")) p=new CMsingleton(in);
 else if (!strcmp(name,"GRAM_SPARSE_WITHOUTDIAG")) p=new CMgramsparse_withoutdiag(in);
 else {
     if (CH_Matrix_Classes::materrout) (*CH_Matrix_Classes::materrout)<<"*** ERROR: coeffmat_read(): unknown constraint name :"<<name<<std::endl;
     in.clear(std::ios::failbit);
     return 0;
 }
 if (p==0) {
     if (CH_Matrix_Classes::materrout) (*CH_Matrix_Classes::materrout)<<"*** ERROR: coeffmat_read():";
     if (CH_Matrix_Classes::materrout) (*CH_Matrix_Classes::materrout)<<" failed in reading a constraint of type "<<name<<std::endl;
     in.clear(std::ios::failbit);
     return 0;
 }
 return p;
}




}

