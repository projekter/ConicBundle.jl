/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  Matrix/aasen.cxx
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



#include <stdlib.h>
#include "mymath.hxx"
#include "symmat.hxx"

 

namespace CH_Matrix_Classes {

//Golub/van Loan, Chapter 4.4.3, Matrix computations, 2nd ed., 4th printing
//slight differences: indexing starts at 0 and for i=1:n piv(i) gives the
//index to be exchanged with i in this sequence (GvL exchange with i+1)

int Symmatrix::Aasen_factor(Indexmatrix& piv)
{
  chk_init(*this);

 //TESTING
  //Symmatrix A(*this);   //original matrix
  //Matrix L(nr,nr,0.);   //unit lower triangular, first column e_1
  //L(0,0)=1.;
  //L(nr-1,nr-1)=1.;
  //Symmatrix T(nr,0.);   //tridiagonal 
  //Matrix H(nr,nr,0.);   //H=T*L^T, Hessenberg

 //initialize pivoting permutaion
 piv.init(nr,1,Integer(0));

 if (nr==0)
   return 0;

 for(Integer i=0;i<nr;i++)
   piv[i]=i;

 //provide temporary storage for h
 Real *tmpp;   
 Integer tmp_dim=Integer(memarray->get(nr,tmpp));
 if (tmp_dim<nr)
   MEmessage(MEmem(nr,
		   "Symmatrix::Aasen_factor(Indexmatrix&,Real) not enough memory",MTsymmetric));
 
 //main loop
 for(Integer j=0;j<nr;j++){
   

   //up to j-1 alpha and beta are already overwritten 
   //and L(:,j) is stored in this(:,j-1) (first column of L is e_1)

   //------ compute h(0:j) into tmpp and store alpha(j)
   //std::cout<<"========      j="<<j<<std::endl;
   if (j==0){
     *tmpp=*m;         //=alpha(0)

     //TESTING
     //H(0,0)=T(0,0)=*m;
   }
   else if (j==1){
     *tmpp=*(m+1);  //=this(1,0)=beta(0)
     *(tmpp+1)=*(m+nr);  //=this(1,1)=alpha(1)
     //TESTING
     //H(1,1)=*(m+nr);
     //H(0,1)=*(m+1);
   }
   else {
     //std::cout<<"S="<<*this<<"\n  [l beta a h]"<<std::endl;
     Real hj=0.;  
     Real lprev=0.;    
     Real lk=0.;
     Real betaprev=0;  //=beta(-1)
     //std::cout<<"[ "<<lprev<<" "<<betaprev<<" --  -- ]"<<std::endl;
     Real *hp=tmpp;
     const Real* mp=m;
     for(Integer k=0;k<j;k++){
       Real alphak=*mp++; //this(k,k)
       Real betak=*mp;
       mp+=j-k-1;
       Real lnext=(k<j-1)?*mp:1.; //=L(j,k+1)=this(j,k);
       mp+=nr-j;
       Real hk=betaprev*lprev+alphak*lk+betak*lnext;
       *hp++=hk;

       //TESTING
       //H(k,j)=hk;
       //std::cout<<"[ "<<lk<<" "<<betak<<" "<<alphak<<" "<<hk<<" ]"<<std::endl;

       hj+=lk*hk;
       lprev=lk;
       lk=lnext;
       betaprev=betak;
     }
     *hp=*mp-hj; //+=this(j,j)
     *(m+j*nr-(j*(j-1))/2)=*hp-betaprev*lprev;

     //TESTING
     //H(j,j)=*hp;
     //std::cout<<"[ "<<lk<<" --- "<<*hp-betaprev*lprev<<" "<<*hp<<" ]"<<std::endl;
   }
   //TESTING
   //T(j,j)=*(m+j*nr-(j*(j-1))/2);
   //  std::cout<<"norm2(A(0:j,j)-L(0:j,0:j)*H(0:j,j))=";
   //  std::cout<<norm2((A.col(j))(Range(0,j))-L(Range(0,j),Range(0,j))*H(Range(0,j),Range(j,j)));
   //  std::cout<<std::endl;

   //compute next column and pivot 
   if (j<nr-1){
     //overwrite this(j+1:n-1,j) with v(j+1:n-1)
     Real* mp=m+nr*j-(j*(j-1))/2+1;


     for(Integer k=0;k<j;k++)
       mat_xpeya(nr-j-1,mp,m+nr*k-(k*(k-1))/2+j+1-k,-tmpp[k+1]);


     //TESTING
     //Matrix v(nr-j-1,1,mp);
     //std::cout<<"norm2(v-(A(j+1:nr-1,j)-L(j+1:nr-1,0:j)*H(0:j,j)))=";
     //std::cout<<norm2(v-((A.col(j))(Range(j+1,nr-1))-L(Range(j+1,nr-1),Range(0,j))*H(Range(0,j),Range(j,j))));
     //std::cout<<std::endl;

     //find the maximum element in v
     Integer maxind=j+1;
     Real maxval=abs(*mp++);
     for(Integer k=j+2;k<nr;k++){
       Real val=abs(*mp++);
       if (val>maxval){
	 maxval=val;
	 maxind=k;
       }
     }
     //swap if necessary
     if (maxind>j+1){
       piv(j+1)=maxind;
       swapij(j+1,maxind);

       //TESTING
       //std::cout<<"swap "<<j+1<<" and "<<maxind<<std::endl;
       //A.swapij(j+1,maxind);
       //L.swap_rowsij(j+1,maxind);
     }
     //beta(j) is at its correct position already

     //TESTING
     //T(j+1,j)=H(j+1,j)=*(m+nr*j-(j*(j-1))/2+1);      
     
   }
   //compute L(j+2:n-1,j+1) into this(j+2:n-1,j) (v is there already)
   if (j<nr-2){
     Real *mp=m+nr*j-(j*(j-1))/2+1;
     Real betaj=*mp++;
     if (betaj!=0.){
       mat_xmultea(nr-j-2,mp,1./betaj);
     }

     //TESTING
     //L(j+1,j+1)=1;
     //for(Integer i=j+2;i<nr;i++)
     //  L(i,j+1)=(*this)(i,j);
   }
   //TESTING
   //if (j>0){
   //  std::cout<<"norm2(H(:,0:j-1)-T(:,0:j)*transpose(L(0:j-1,0:j)))=";
   //  std::cout<<norm2(H.cols(Range(0,j-1))-T.cols(Range(0,j))*CH_Matrix_Classes::transpose(L(Range(0,j-1),Range(0,j))))<<std::endl;
   //}
   //std::cout<<"norm2(A(:,0:j)-L(:,0:j+1)*H(0:j+1,0:j))=";
   //std::cout<<norm2(A.cols(Range(0,j))-L.cols(Range(0,min(j+1,nr-1)))*H(Range(0,min(j+1,nr-1)),Range(0,j)))<<std::endl;
              
 }

 memarray->free(tmpp); tmpp=0;
 
 //TESTING
 //Matrix L1(nr,nr,0.);
 //Symmatrix T1(nr,0.);
 //L1(0,0)=1.;
 //T1(0,0)=(*this)(0,0);
 //for(Integer i=1;i<nr;i++){
 //  L1(i,i)=1.;
 //  T1(i-1,i)=(*this)(i-1,i);
 //  T1(i,i)=(*this)(i,i);
 //  for (Integer j=i+1;j<nr;j++){
 //    L1(j,i)=(*this)(j,i-1);
 //  }
 //} 
 //L1.display(std::cout);
 //L.display(std::cout);
 //std::cout<<"norm2(T1-T)="<<norm2(T1-T)<<std::endl;
 //std::cout<<"norm2(L1-L)="<<norm2(L1-L)<<std::endl;
 //std::cout<<"norm2(A-L1*T1*tranpose(L1))="<<norm2(A-L1*T1*CH_Matrix_Classes::transpose(L1))<<std::endl;

 return 0;
}


int Symmatrix::Aasen_Lsolve(Matrix& x) const
{
 chk_mult(*this,x);

 //solve for and overwrite column k of x
 for(Integer k=0;k<x.coldim();k++){ 
   
   Real* xbase=x.m+k*nr;
   
   //---- solve Lr=xbase (first column is e_1, next column is in i-1)
   for(Integer i=1;i<nr;i++){
         Real f=xbase[i];
         Real *rp=xbase+1;
         const Real *mp=m+i;
         const Integer h=nr-i+1;
         for(Integer j=nr;--j>=h;) {
             f-=(*rp++)*(*mp);
             mp+=j;
         }
         *rp=f;
     }
 }

 return 0;
}

  int Symmatrix::Aasen_Ltsolve(Matrix& x) const
{
 chk_mult(*this,x);

 //solve for and overwrite column k of x (first column is e_1, next column is in i-1)
 for(Integer k=0;k<x.coldim();k++){
 
     Real* xbase=x.m+k*nr;
     //---- solve L'r=r
     {for(Integer i=nr;--i>0;){
         Real *rp=xbase+i;
         Real f=(*rp++);
         const Real *mp=m+(i-1)*nr-((i-1)*(i-2))/2+2;
         for(Integer j=nr-i;--j>0;) 
	   f-=(*rp++)*(*mp++);
         xbase[i]=f;
     }}
 }

 return 0;
}

int Symmatrix::Aasen_tridiagsolve(Matrix& x) const
{
 chk_mult(*this,x);

 if (nr==0)
   return 0;
 
 if (nr==1){
   if(*m==0.){
     if (max(abs(x))>0.)
       return -1;
     else 
       return 0;
   }
   x/=*m;
   if (x.contains_nan()){
     return -1;
   }
   return 0;
 }

 //solve for and overwrite column k of x
 //Matrix tmp(3,nr,0.);
 //provide temporary storage for h
 Real *tmpp;   
 Integer tmp_dim=Integer(memarray->get(3*nr,tmpp));
 if (tmp_dim<3*nr)
   MEmessage(MEmem(3*nr,
		   "Symmatrix::Aasen_factor(Indexmatrix&,Real) not enough memory",MTsymmetric));

 for (Integer k=0;k<x.coldim();k++){
   //forward factorization
   Real* xbase=x.m+k*nr;
   Real* xp=xbase;
   const Real* mp=m;
   Real a1=*mp;
   Real b1=*(mp+1);
   Real x1=*xp;
   for (Integer i=0;i<nr-1;i++,xp++){
     Real x2=*(xp+1);
     Real c2=*(mp+1);
     mp+=nr-i;
     Real a2=*(mp);
     Real b2;
     if (i<nr-1) 
       b2=*(mp+1);
     else
       b2=0;
     if (abs(a1)>=abs(c2)){
       *tmpp++=a1; //tmp(0,i)=a1;
       *tmpp++=b1; //tmp(1,i)=b1;
       *tmpp++=0.; //tmp(2,i)=0.;
       *xp=x1;
       if (a1!=0.){
	 Real f=-c2/a1;
	 if (isnan(f)){
	   memarray->free(tmpp); tmpp=0;
	   return -(i+1);
	 }
	 a1=a2+b1*f;
	 b1=b2;
	 x1=x2+x1*f;
       }
       else {
	 a1=a2;
	 b1=b2;
	 x1=x2;
       }
     }
     else { //abs(c2)>abs(a1), swap the two rows
       *tmpp++=c2; //tmp(0,i)=c2;
       *tmpp++=a2; //tmp(1,i)=a2;
       *tmpp++=b2; //tmp(2,i)=b2;
       *xp=x2;
       Real f=-a1/c2;
       if (isnan(f)){
	 memarray->free(tmpp); tmpp=0;
	 return -(i+1);
       }
       a1=b1+f*a2;
       b1=f*b2;
       x1+=x2*f;
     }
   }
   //backward solve
   //tmp(0,nr-1)=a1;
   //tmp(1,nr-1)=0.;
   //tmp(2,nr-1)=0.;
   *xp=x1;
   if (a1==0.){
     if (x1!=0.){
       memarray->free(tmpp); tmpp=0;
       return -nr;
     }
     else
       *xp=0.;
   }
   else {
     *xp/=a1;
     if (isnan(*xp)){
       memarray->free(tmpp); tmpp=0;
       return -nr;
     }
   }
   Real d2=*xp;
   xp--;
   --tmpp; //remove tmp(2,nr-2), which is not needed
   *xp-=*(--tmpp)*d2; // *xp-=tmp(1,nr-2)*d2;
   Real r=*(--tmpp);
   if (r==0.){
     if (*xp!=0.){
       memarray->free(tmpp); tmpp=0;
       return -(nr-1);
     }
     else
       *xp=0.;
   }
   else {
     *xp/=r;
     if (isnan(*xp)){
       memarray->free(tmpp); tmpp=0;
       return -(nr-1);
     }
   }
   Real d1=*xp;
   xp--;
   for(Integer i=nr-2;--i>=0;){
     Real d0=*xp;
     d0-=d2* *(--tmpp);  //tmp(2,i);
     d0-=d1* *(--tmpp);  //tmp(1,i);
     Real r=*(--tmpp);
     if (r==0.){
       if (d0!=0.){
	 memarray->free(tmpp); tmpp=0;
	 return -(i+1);
       }
       else
	 d0=0.;
     }
     else {
       d0/=r;
       if (isnan(d0)){
	 memarray->free(tmpp); tmpp=0;
	 return -(i+1);
       }
     }
     *xp=d0;
     xp--;
     d2=d1;
     d1=d0;
   }
 }

 memarray->free(tmpp); tmpp=0;
 
 return 0;
}  

int Symmatrix::Aasen_solve(Matrix& x,const Indexmatrix& piv) const
{
 chk_mult(*this,x);

 if (nr==0)
   return 0;
 
 if (nr==1){
   if(*m==0.){
     if (max(abs(x))>0.)
       MEmessage(MatrixError(ME_num,"aasen_solve: division by zero",MTsymmetric));
     else 
       return 0;
   }
   x/=*m;
   if (x.contains_nan()){
     return -1;
   }
   return 0;
 }

 x.pivot_permute_rows(piv);
 
 int retval=Aasen_Lsolve(x);
 if (retval!=0){
   return retval;
 }
 if (x.contains_nan()){
   return 1;
 }
 retval=Aasen_tridiagsolve(x);
 if (retval!=0){
   return retval;
 }
 if (x.contains_nan()){
   return 2;
 }
 retval=Aasen_Ltsolve(x);
 if (x.contains_nan()){
   return 3;
 }

 x.pivot_permute_rows(piv,true);
   
 return 0;
}   
   

}

