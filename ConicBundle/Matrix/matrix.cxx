/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  Matrix/matrix.cxx
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
#include <algorithm>
#include <random>
#include "matrix.hxx"
#include "heapsort.hxx"
#include "gb_rand.hxx"

#ifdef WITH_OMP
#include <omp.h>
#endif

 
using namespace CH_Tools;

namespace CH_Matrix_Classes {

const Mtype Matrix::mtype = MTmatrix;

// **************************************************************************
//                                Constructors
// **************************************************************************

  Matrix& Matrix::xeya(const Matrix& A,Real d,int atrans)
{
 chk_init(A);
 if ((atrans==0)||(A.nr<=1)||(A.nc<=1)){
   if (atrans==0){
     newsize(A.nr,A.nc);
   }
   else {
     newsize(A.nc,A.nr);
   }
   chk_set_init(*this,1);
   mat_xeya(nr*nc,m,A.m,d);
 }
 else {
   newsize(A.nc,A.nr);
   chk_set_init(*this,1);
   if (d==0.) {
     mat_xea(nr*nc,m,0.);
   }
   else {
     if (nr<nc){
       if (d==1.) {
	 for(Integer i=0;i<nr;i++)
	   mat_xey(nc,m+i,nr,A.m+i*nc,1);
       }
       else if (d==-1.) {
	 for(Integer i=0;i<nr;i++)
	   mat_xemy(nc,m+i,nr,A.m+i*nc,1);
       }
       else {
	 for(Integer i=0;i<nr;i++)
	   mat_xeya(nc,m+i,nr,A.m+i*nc,1,d);
       }
     }
     else {
       if (d==1.) {
	 for(Integer i=0;i<nc;i++)
	   mat_xey(nr,m+i*nr,1,A.m+i,nc);
       }
       else if (d==-1.) {
	 for(Integer i=0;i<nr;i++)
	   mat_xemy(nr,m+i*nr,1,A.m+i,nc);
       }
       else {
	 for(Integer i=0;i<nr;i++)
	   mat_xeya(nr,m+i*nr,1,A.m+i,nc,d);
       }
     }
   } 
 }  
 return *this;
}

Matrix& Matrix::xpeya(const Matrix& A,Real d)
{
 chk_add(*this,A);
 mat_xpeya(nr*nc,m,A.m,d);
 return *this;
}

Matrix& Matrix::xeya(const Indexmatrix& A,Real d)
{
 chk_init(A);
 newsize(A.nr,A.nc);
 chk_set_init(*this,1);
 if (d==1.) { for(Integer i=0;i<nr*nc;i++) m[i]=Real(A.m[i]); return *this;}
 if (d==0.) { mat_xea(nr*nc,m,0.); return *this;}
 if (d==-1.) { for(Integer i=0;i<nr*nc;i++) m[i]=Real(-A.m[i]); return *this;}
 for(Integer i=0;i<nr*nc;i++) m[i]=d*A.m[i];
 return *this;
}

Matrix& Matrix::xpeya(const Indexmatrix& A,Real d)
{
 chk_add(*this,A);
 if (d==1.) { for(Integer i=0;i<nr*nc;i++) m[i]+=Real(A.m[i]); return *this;}
 if (d==0.) { return *this;}
 if (d==-1.) { for(Integer i=0;i<nr*nc;i++) m[i]-=A.m[i]; return *this;}
 for(Integer i=0;i<nr*nc;i++) m[i]+=d*A.m[i];
 return *this;
}

Matrix& xbpeya(Matrix& x,const Matrix& y,Real alpha,Real beta,int ytrans)
  //returns x= alpha*y+beta*x, where y may be transposed (ytrans=1)
  //if beta==0. then x is initialized to the correct size
{
  chk_init(y);
  if (beta==0.){
    if (!ytrans){ //y is not transposed
      x.newsize(y.nr,y.nc);
      if (alpha==0.){
	mat_xea(x.dim(),x.m,0.);
      }
      else if (alpha==1.){
	mat_xey(y.dim(),x.m,y.m);
      }
      else {
	mat_xeya(y.dim(),x.m,y.m,alpha);
      }
      chk_set_init(x,1);
      return x;
    }
    else{ //y is transposed
      x.newsize(y.nc,y.nr);
      if (alpha==0.){
	mat_xea(x.dim(),x.m,0.);
      }
      else if (alpha==1.){
	if (x.nr>x.nc){
	  for (Integer i=0;i<x.nc;i++){
	    mat_xey(x.nr,x.m+i*x.nr,1,y.m+i,x.nc);
	  }
	} 
	else {
	  for (Integer i=0;i<x.nr;i++){
	    mat_xey(x.nc,x.m+i,x.nr,y.m+i*x.nc,1);
	  }
	}
      }
      else {
	if (x.nr>x.nc){
	  for (Integer i=0;i<x.nc;i++){
	    mat_xeya(x.nr,x.m+i*x.nr,1,y.m+i,x.nc,alpha);
	  }
	} 
	else {
	  for (Integer i=0;i<x.nr;i++){
	    mat_xeya(x.nc,x.m+i,x.nr,y.m+i*x.nc,1,alpha);
	  }
	}
      }
      chk_set_init(x,1);
      return x;
    }
  }
  //Now beta!=0
  if (!ytrans) { //y is not transposed
    chk_add(x,y);
    if (alpha==0.){
      if (beta==1.) return x;
      else if (beta==-1.) {
	mat_xemx(x.dim(),x.m);
        return x;
      }
      mat_xmultea(x.dim(),x.m,beta);
      return x;   
    }
    else if (beta==1.){
      if (alpha==1.){
	mat_xpey(x.dim(),x.m,y.m);
        return x;
      }
      else if (alpha==-1.){
	mat_xmey(x.dim(),x.m,y.m);
        return x;
      }
      mat_xpeya(x.dim(),x.m,y.m,alpha);
      return x; 
    }
    mat_xbpeya(x.dim(),x.m,y.m,alpha,beta);
    return x; 
  } 
  //Now y transposed
  chk_init(x);
#if (CONICBUNDLE_DEBUG>=1)
  if ((x.nr!=y.nc)||(x.nc!=y.nr)){
    MEmessage(MatrixError(ME_dim,"xbpeya: dimensions don't match",MTmatrix));;
  }
#endif 
  if (alpha==0.){
    if (beta==1.) return x;
    if (beta==-1.) {
      mat_xemx(x.dim(),x.m);
      return x;
    }
    mat_xmultea(x.dim(),x.m,beta);
    return x;   
  }
  if (beta==1.){
    if (alpha==1.){
      if (x.nr>x.nc){
	for (Integer i=0;i<x.nc;i++){
	  mat_xpey(x.nr,x.m+i*x.nr,1,y.m+i,x.nc);
	}
      } 
      else {
	for (Integer i=0;i<x.nr;i++){
	  mat_xpey(x.nc,x.m+i,x.nr,y.m+i*x.nc,1);
	}
      }
      return x;
    }
    if (alpha==-1.){
      if (x.nr>x.nc){
	for (Integer i=0;i<x.nc;i++){
	  mat_xmey(x.nr,x.m+i*x.nr,1,y.m+i,x.nc);
	}
      } 
      else {
	for (Integer i=0;i<x.nr;i++){
	  mat_xmey(x.nc,x.m+i,x.nr,y.m+i*x.nc,1);
	}
      }
      return x;
    }
      
    if (x.nr>x.nc){
      for (Integer i=0;i<x.nc;i++){
	mat_xpeya(x.nr,x.m+i*x.nr,1,y.m+i,x.nc,alpha);
      }
    } 
    else {
      for (Integer i=0;i<x.nr;i++){
	mat_xpeya(x.nc,x.m+i,x.nr,y.m+i*x.nc,1,alpha);
      }
    }
    return x;
  }
  //now beta!=1. 
  if (x.nr>x.nc){
    for (Integer i=0;i<x.nc;i++){
      mat_xbpeya(x.nr,x.m+i*x.nr,1,y.m+i,x.nc,alpha,beta);
    }
  } 
  else {
    for (Integer i=0;i<x.nr;i++){
      mat_xbpeya(x.nc,x.m+i,x.nr,y.m+i*x.nc,1,alpha,beta);
    }
  }
  return x;
}
  
Matrix& xeyapzb(Matrix& x,const Matrix& y,const Matrix& z,Real alpha,Real beta)
  //returns x= alpha*y+beta*z,
  //x is initialized to the correct size
{
  chk_add(y,z);
  x.newsize(y.nr,y.nc);
  chk_set_init(x,1);
  if (alpha==1.){
    if(beta==1.){
      mat_xeypz(x.dim(),x.m,y.m,z.m);
      return x;
    }
    if (beta==-1.){
      mat_xeymz(x.dim(),x.m,y.m,z.m);
      return x;
    }
  }
  if ((beta==1.)&&(alpha==-1.)){
    mat_xeymz(x.dim(),x.m,z.m,y.m);
    return x;
  }
  mat_xeyapzb(x.dim(),x.m,y.m,z.m,alpha,beta);
  return x;
}


Matrix& genmult(const Matrix& A,const Matrix& B,Matrix& C,
                    Real alpha,Real beta,int atrans,int btrans)
            //returns C=beta*C+alpha*A*B, where A and B may be transposed
            //C may neither be equal to A nor B
#if defined(WITH_BLAS)|| defined(BLAS_GENMULT)
{
 chk_init(A);
 chk_init(B);
 Integer nr=(atrans?A.nc:A.nr);
 Integer nc=(btrans?B.nr:B.nc);
 Integer nk=(atrans?A.nr:A.nc);
#if (CONICBUNDLE_DEBUG>=1)
 if (beta!=0.) chk_init(C);
 if ((nk!=(btrans?B.nc:B.nr))||
     ((beta!=0.)&&((C.nr!=nr)||(C.nc!=nc)))
     ){
   MEmessage(MatrixError(ME_dim,"genmult: dimensions don't match",MTmatrix));;
 }
#endif
 if ((nr==0)||(nc==0)||(nk==0)) {
   if (beta==0.) return C.init(nr,nc,0.);
   return C*=beta;
 }
 if (beta==0.) {
   C.newsize(nr,nc);
   chk_set_init(C,1);
 }
 cblas_dgemm(CblasColMajor,(atrans?CblasTrans:CblasNoTrans),
	     (btrans?CblasTrans:CblasNoTrans),nr,nc,nk,alpha,
	     A.get_store(),A.nr,B.get_store(),B.nr,beta,
	     C.get_store(),C.nr);
 return C;
}
#else
{
 chk_init(A);
 chk_init(B);
 Integer nr,nm,nc;
 if (atrans) {nr=A.nc;nm=A.nr;}
 else { nr=A.nr; nm=A.nc;}
 if (btrans) {
     nc=B.nr;
#if (CONICBUNDLE_DEBUG>=1)
     if (nm!=B.nc) {
         MEmessage(MatrixError(ME_dim,"genmult: dimensions don't match",MTmatrix));;
     }
#endif
 }
 else {
     nc=B.nc;
#if (CONICBUNDLE_DEBUG>=1)
     if (nm!=B.nr) {
         MEmessage(MatrixError(ME_dim,"genmult: dimensions don't match",MTmatrix));;
     }
#endif
 }
 if ((nr==0)||(nc==0)||(nm==0)) {
   if (beta==0.) return C.init(nr,nc,0.);
   return C*=beta;
 }
 bool initialize;
 if (beta!=0.){
   initialize=false;
   chk_init(C);
#if (CONICBUNDLE_DEBUG>=1)
   if ((nr!=C.nr)||(nc!=C.nc)) {
     MEmessage(MatrixError(ME_dim,"genmult: dimensions don't match",MTmatrix));;
   }
#endif
   if (beta!=1.) 
     C*=beta;
 }
 else {
   if (nm==0){
     C.init(nr,nc,0.);
     return C;
   }
   initialize=true;
   C.newsize(nr,nc); chk_set_init(C,1);
 }
 if (alpha==0.) return C;
 //deal with the cases where at least one of the matrices is a vector
 if (nr==1) {
   if (nc==1) {
     if (initialize)
       *(C.m)=alpha*mat_ip(nm,A.m,B.m);
     else
       *(C.m)+=alpha*mat_ip(nm,A.m,B.m);
     return C;
   }
   if (btrans) {
     Real *ap=A.m;
     Real *bp=B.m;
     Integer i=0;
     if (initialize){
       mat_xeya(nc,C.m,bp,alpha* *ap++);
       i++,bp+=nc;
     }
     for(;i<nm;i++,bp+=nc)
       mat_xpeya(nc,C.m,bp,alpha* *ap++);
     return C;
   }
   Real *cp=C.m;
   Real *bp=B.m;
   if (initialize){
     for(Integer i=0;i<nc;i++,bp+=nm)
       *cp++ = alpha*mat_ip(nm,A.m,bp);
   }
   else {
     for(Integer i=0;i<nc;i++,bp+=nm)
       *cp++ += alpha*mat_ip(nm,A.m,bp);
   }
   return C;
 }
 if (nc==1) {
   if (atrans){
     Real *cp=C.m;
     Real *ap=A.m;
     if (initialize){
       for(Integer i=0;i<nr;i++,ap+=nm)
	 *cp++ = alpha*mat_ip(nm,B.m,ap);
     }
     else {
       for(Integer i=0;i<nr;i++,ap+=nm)
	 *cp++ += alpha*mat_ip(nm,B.m,ap);
     }
     return C;
   }
   Real *ap=A.m;
   Real *bp=B.m;
   Integer i=0;
   if (initialize){
     mat_xeya(nr,C.m,ap,alpha* *bp++);
     i++,ap+=nr;
   }
   for(;i<nm;i++,ap+=nr)
     mat_xpeya(nr,C.m,ap,alpha* *bp++);
   return C;
 }     
 //now neither A nor B are vectors
 if (atrans){
   if (btrans){
     //Matrix tmpB(B,1.,1);
     Real *cp=C.m;
     if (initialize) {
       for(Integer j=0;j<nc;j++){
	 for(Integer i=0;i<nr;i++){
	   *cp++ =alpha*mat_ip(nm,A.m+i*nm,1,B.m+j,nc);
	   //*cp++ =alpha*mat_ip(nm,A.m+i*nm,tmpB.m+j*nm);
	 }
       }
     }
     else {
       for(Integer j=0;j<nc;j++){
	 for(Integer i=0;i<nr;i++){
	   *cp++ +=alpha*mat_ip(nm,A.m+i*nm,1,B.m+j,nc);
	   //*cp++ +=alpha*mat_ip(nm,A.m+i*nm,tmpB.m+j*nm);
	 }
       }
     }
   }
   else {
#ifdef WITH_OMP
     //if (true){
     if ((nr<100)||(nm<30)){
#endif
       Real *cp=C.m;
       if (initialize){
	 for(Integer j=0;j<nc;j++){
	   for(Integer i=0;i<nr;i++){
	     *cp++ =alpha*mat_ip(nm,A.m+i*nm,B.m+j*nm);
	   }
	 }
       }
       else {
	 for(Integer j=0;j<nc;j++){
	   for(Integer i=0;i<nr;i++){
	     *cp++ +=alpha*mat_ip(nm,A.m+i*nm,B.m+j*nm);
	   }
	 }
       }
#ifdef WITH_OMP
     }
     else {
       for(Integer j=0;j<nc;j++){
	 Real *const cp=C.m+j*nr;
	 const Real* const bp=B.m+j*nm;
	 const Real* const ap=A.m;
	 if (initialize){
#pragma omp parallel for schedule(static,1)
	   for(Integer i=0;i<nr;i++){
	     cp[i]=alpha*mat_ip(nm,ap+i*nm,bp);
	   }
	 }
	 else {
#pragma omp parallel for schedule(static,1)
	   for(Integer i=0;i<nr;i++){
	     cp[i]+=alpha*mat_ip(nm,ap+i*nm,bp);
	   }
	 }
       }
     }
#endif
   }
 }
 else {
   //Matrix tmpA(A,1.,1);
   if (btrans){
     //Matrix tmpB(B,1.,1);
     Real *cp=C.m;
     if (initialize){
       for(Integer j=0;j<nc;j++){
	 for(Integer i=0;i<nr;i++){
	   *cp++ =alpha*mat_ip(nm,A.m+i,nr,B.m+j,nc);
	   //*cp++ =alpha*mat_ip(nm,tmpA.m+i*nm,tmpB.m+j*nm);
	 }
       }
     }
     else {
       for(Integer j=0;j<nc;j++){
	 for(Integer i=0;i<nr;i++){
	   *cp++ +=alpha*mat_ip(nm,A.m+i,nr,B.m+j,nc);
	   //*cp++ +=alpha*mat_ip(nm,tmpA.m+i*nm,tmpB.m+j*nm);
	 }
       }
     }
   }
   else {
     Real *cp=C.m;
     if (initialize){
       for(Integer j=0;j<nc;j++){
	 for(Integer i=0;i<nr;i++){
	   *cp++ =alpha*mat_ip(nm,A.m+i,nr,B.m+j*nm,1);
	   //*cp++ =alpha*mat_ip(nm,tmpA.m+i*nm,B.m+j*nm);
	 }
       }
     }
     else {
       for(Integer j=0;j<nc;j++){
	 for(Integer i=0;i<nr;i++){
	   *cp++ +=alpha*mat_ip(nm,A.m+i,nr,B.m+j*nm,1);
	   //*cp++ +=alpha*mat_ip(nm,tmpA.m+i*nm,B.m+j*nm);
	 }
       }
     }
   }
 }
 return C;
}
#endif
         

Matrix& Matrix::init(const Realrange& r)
{
 Integer i;
 if (r.step>=0) {
     if (r.from>r.to+r.tol) i=0;
     else i=Integer((r.to+r.tol-r.from)/r.step)+1;
 }
 else {
     if (r.from<r.to-r.tol) i=0;
     else i=Integer((r.to-r.tol-r.from)/r.step)+1;
 }
 newsize(i,1);
 Real d=r.from;
 Real *mp=m;
 for(i=nr;--i>=0;){
     (*mp++)=d;
     d+=r.step;
 }
 chk_set_init(*this,1);
 return *this;
}

void Matrix::newsize(Integer inr,Integer inc)
{
 chk_range(inr,inc,-1,-1);
 chk_set_init(*this,0);
 aux_task=none;
 if ((inr==0)||(inc==0)){
     nr=inr;
     nc=inc;
     return;
 }
 if ((inr!=nr)||(inc!=nc)) {
     nr=inr;
     nc=inc;
     if (nr*nc>mem_dim){
         memarray->free(m); m=0;
         mem_dim=Integer(memarray->get(nr*nc,m));
         if (mem_dim<nr*nc)
             MEmessage(MEmem(nr*nc,
                         "Matrix::Matrix(Integer,Integer,Real) not enough memory",MTmatrix));
     }
 }
}

Matrix Matrix::operator()(const Indexmatrix &rv,const Indexmatrix &cv) const
{
 chk_init(rv);
 chk_init(cv);
 chk_init(*this);
 if ((rv.dim()==0)||(cv.dim()==0)) return Matrix(0,0,0.);
 chk_range(min(rv),min(cv),nr,nc);
 chk_range(max(rv),max(cv),nr,nc);
 Matrix A(rv.dim(),cv.dim());
 Integer i,j;
 Real *ap=A.m;
 Integer *rp;
 Integer *cp=cv.m;
 Real *mcp;           //points to current column indexed by cv
 for(j=A.nc;--j>=0;){
     mcp=m+nr*(*cp++);
     rp=rv.m;
     for(i=A.nr;--i>=0;)
         (*ap++)=mcp[*rp++];
 }
 chk_set_init(A,1);
 return A;
}

Matrix Matrix::operator()(const Indexmatrix &v) const
{
 chk_init(v);
 chk_init(*this);
 if (v.dim()==0) return Matrix(v.nr,v.nc,0.);
 chk_range(min(v),max(v),nr*nc,nr*nc);
 Matrix A(v.nr,v.nc);
 Integer i;
 Real *ap=A.m;
 Integer *vp=v.m;
 for(i=A.nr*A.nc;--i>=0;){
         (*ap++)=m[*vp++];
 }
 chk_set_init(A,1);
 return A;
}

Matrix Matrix::col(Integer c) const
{
 chk_init(*this);
 chk_range(c,0,nc,1);
 return Matrix(nr,1,m+c*nr);
}

Matrix Matrix::row(Integer r) const
{
 chk_init(*this);
 chk_range(r,0,nr,1);
 return Matrix(1,nc,m+r,nr);
}

Matrix Matrix::cols(const Indexmatrix& v) const
{
 chk_init(v);
 chk_init(*this);
 if (v.dim()==0) return Matrix(nr,v.dim(),0.);
 chk_range(min(v),max(v),nc,nc);
 Matrix A(nr,v.dim());
 Real *mp;
 Real *ap=A.m;
 Integer *vp=v.m;
 Integer i,j;
 for(i=A.nc;--i>=0;){
     mp=m+nr*(*vp++);
     for(j=nr;--j>=0;)
         *ap++=*mp++;
 }
 chk_set_init(A,1);
 return A;
}

Matrix Matrix::rows(const Indexmatrix& v) const
{
 chk_init(v);
 chk_init(*this);
 if (v.dim()==0) return Matrix(v.dim(),nc,0.);
 chk_range(min(v),max(v),nr,nr);
 Matrix A(v.dim(),nc);
 Integer *vp=v.m;
 for(Integer i=0;i<A.nr;i++){
     mat_xey(nc,A.m+i,A.nr,m+(*vp++),nr);
 }
 chk_set_init(A,1);
 return A;
}


Matrix& Matrix::swap_rowsij(Integer i,Integer j)
{
  chk_init(*this);
  assert((0<=i)&&(i<nr)&&(0<=j)&&(j<nr));
  if (i!=j)
    mat_swap(nc,m+i,nr,m+j,nr);
  return *this;
}

Matrix& Matrix::swap_colsij(Integer i,Integer j)
{
  chk_init(*this);
  assert((0<=i)&&(i<nc)&&(0<=j)&&(j<nc));
  if (i!=j)
    mat_swap(nr,m+i*nr,m+j*nr);
  return *this;
}


Matrix& Matrix::pivot_permute_rows(const Indexmatrix & piv, bool inverse)
{ 
  chk_init(*this);
  chk_init(piv);
  assert(piv.dim()==nr);

  if(!inverse){
    for(Integer i=0;i<nr;i++){
      Integer h=piv(i);
      if (h==i) continue;
      assert((0<=h)&&(h<nr));
      mat_swap(nc,m+i,nr,m+h,nr);
    }
  }
  else {
    for(Integer i=nr;--i>=0;){
      Integer h=piv(i);
      if (h==i) continue;
      assert((0<=h)&&(h<nr));
      mat_swap(nc,m+i,nr,m+h,nr);
    }
  }    

  return *this;
}

Matrix& Matrix::pivot_permute_cols(const Indexmatrix & piv, bool inverse)
{ 
  chk_init(*this);
  chk_init(piv);
  assert(piv.dim()==nc);

  if(!inverse){
    for(Integer i=0;i<nc;i++){
      Integer h=piv(i);
      if (h==i) continue;
      assert((0<=h)&&(h<nc));
      mat_swap(nr,m+i*nr,m+h*nr);
    }
  }
  else {
    for(Integer i=nc;--i>=0;){
      Integer h=piv(i);
      if (h==i) continue;
      assert((0<=h)&&(h<nc));
      mat_swap(nr,m+i*nr,m+h*nr);
    }
  }

  return *this;
}


Matrix& Matrix::transpose()
{
 chk_init(*this);
 if ((nr<=1)||(nc<=1)){
     Integer h=nr;nr=nc;nc=h;
     return *this;
 }
 Real *nm;
 Integer nmem_dim=Integer(memarray->get(nr*nc,nm));
 if (nmem_dim<nr*nc)
     MEmessage(MEmem(nr*nc,"Matrix::transpose() not enough memory",MTmatrix));
 Integer i,j;
 Real *mp=m;
 Real *nmp;
 for(j=0;j<nc;j++){
     nmp=nm+j;
     for(i=0;i<nr;i++,nmp+=nc)
         *nmp=*mp++;
 }
 i=nr;nr=nc;nc=i;
 memarray->free(m);
 m=nm;
 mem_dim=nmem_dim;
 return *this;
}

Matrix& Matrix::subassign(const Indexmatrix &rv,const Indexmatrix &cv,
                          const Matrix& A)
{
 chk_init(rv);
 chk_init(cv);
 chk_init(A);
 if ((rv.dim()==0)||(cv.dim()==0)) return *this;
 chk_range(min(rv),min(cv),nr,nc);
 chk_range(max(rv),max(cv),nr,nc);
#if (CONICBUNDLE_DEBUG>=1)
 if ((rv.dim()!=A.nr)||(cv.dim()!=A.nc))
     MEmessage(MEdim(rv.dim(),cv.dim(),A.nr,A.nc,
                 "Matrix::subassign(const Indexmatrix&,const Indexmatrix&,const Matrix&) dimensions do not match",
                 MTmatrix));
#endif
 Integer i,j;
 Real *ap=A.m;
 Integer *rp;
 Integer *cp=cv.m;
 Real *mcp;           //points to current column indexed by cv
 for(j=A.nc;--j>=0;){
     mcp=m+nr*(*cp++);
     rp=rv.m;
     for(i=A.nr;--i>=0;)
         mcp[*rp++]=(*ap++);
 }
 return *this;
}

Matrix& Matrix::subassign(const Indexmatrix &v,const Matrix& A)
{
 chk_init(v);
 chk_init(A);
 if (v.dim()==0) return *this;
 chk_range(min(v),max(v),nr*nc,nr*nc);
#if (CONICBUNDLE_DEBUG>=1)
 if (v.dim()!=A.dim())
     MEmessage(MEdim(v.dim(),1,A.dim(),1,
                 "Matrix::subassign(const Indexmatrix&,const Matrix&) dimensions do not match",
                 MTmatrix));
#endif
 Integer i;
 Real *ap=A.m;
 Integer *vp=v.m;
 for(i=A.nr*A.nc;--i>=0;){
         m[*vp++]=(*ap++);
 }
 return *this;
}

Matrix operator<(const Matrix& A,const Matrix& B) 
{
 chk_add(A,B);
 Matrix C(A.nr,A.nc);
 Real *ap=A.m;
 Real *bp=B.m;
 Real *cp=C.m;
 Integer i=A.nr*A.nc;
 while(--i>=0)
     (*cp++)=Real((*ap++)<(*bp++));
 chk_set_init(C,1);
 return C;
}

Matrix operator<=(const Matrix& A,const Matrix& B) 
{
 chk_add(A,B);
 Matrix C(A.nr,A.nc);
 Real *ap=A.m;
 Real *bp=B.m;
 Real *cp=C.m;
 Integer i=A.nr*A.nc;
 while(--i>=0)
     (*cp++)=Real((*ap++)<=(*bp++));
 chk_set_init(C,1);
 return C;
}

Matrix operator==(const Matrix& A,const Matrix& B) 
{
 chk_add(A,B);
 Matrix C(A.nr,A.nc);
 Real *ap=A.m;
 Real *bp=B.m;
 Real *cp=C.m;
 Integer i=A.nr*A.nc;
 while(--i>=0)
     (*cp++)=Real((*ap++)==(*bp++));
 chk_set_init(C,1);
 return C;
}

Matrix operator!=(const Matrix& A,const Matrix& B) 
{
 chk_add(A,B);
 Matrix C(A.nr,A.nc);
 Real *ap=A.m;
 Real *bp=B.m;
 Real *cp=C.m;
 Integer i=A.nr*A.nc;
 while(--i>=0)
     (*cp++)=Real((*ap++)!=(*bp++));
 chk_set_init(C,1);
 return C;
}


Matrix operator<(const Matrix& A,Real d)
{
 chk_init(A);
 Matrix B(A.nr,A.nc);
 Real *bp=B.m;
 Real *ap=A.m;
 Integer i=A.nr*A.nc;
 while(--i>=0)
     (*bp++)=Real((*ap++)<d);
 chk_set_init(B,1);
 return B;
}

Matrix operator>(const Matrix& A,Real d)
{
 chk_init(A);
 Matrix B(A.nr,A.nc);
 Real *bp=B.m;
 Real *ap=A.m;
 Integer i=A.nr*A.nc;
 while(--i>=0)
     (*bp++)=Real((*ap++)>d);
 chk_set_init(B,1);
 return B;
}

Matrix operator<=(const Matrix& A,Real d)
{
 chk_init(A);
 Matrix B(A.nr,A.nc);
 Real *bp=B.m;
 Real *ap=A.m;
 Integer i=A.nr*A.nc;
 while(--i>=0)
     (*bp++)=Real((*ap++)<=d);
 chk_set_init(B,1);
 return B;
}

Matrix operator>=(const Matrix& A,Real d)
{
 chk_init(A);
 Matrix B(A.nr,A.nc);
 Real *bp=B.m;
 Real *ap=A.m;
 Integer i=A.nr*A.nc;
 while(--i>=0)
     (*bp++)=Real((*ap++)>=d);
 chk_set_init(B,1);
 return B;
}

Matrix operator==(const Matrix& A,Real d)
{
 chk_init(A);
 Matrix B(A.nr,A.nc);
 Real *bp=B.m;
 Real *ap=A.m;
 Integer i=A.nr*A.nc;
 while(--i>=0)
     (*bp++)=Real((*ap++)==d);
 chk_set_init(B,1);
 return B;
}

Matrix operator!=(const Matrix& A,Real d)
{
 chk_init(A);
 Matrix B(A.nr,A.nc);
 Real *bp=B.m;
 Real *ap=A.m;
 Integer i=A.nr*A.nc;
 while(--i>=0)
     (*bp++)=Real((*ap++)!=d);
 chk_set_init(B,1);
 return B;
}


// **************************************************************************
//                                 triu
// **************************************************************************

Matrix& Matrix::triu(Integer i)
{
 chk_init(*this);
 Integer j;
 Integer n=min(nc,nr+i-1);
 for(j=0;j<n;j++){
   mat_xea(nr-max(Integer(0),j+1-i),m+j*nr+max(Integer(0),j+1-i),0.);
 }
 return *this;
}

Matrix& Matrix::tril(Integer i)
{
 chk_init(*this);
 Integer j;
 for(j=max(0,i+1);j<nc;j++){
     mat_xea(min(nr,j-i),m+j*nr,0.);
 }
 return *this;
}


Matrix& Matrix::concat_right(const Matrix& A,int Atrans)
{
 chk_init(*this);
 chk_init(A);
 if ((A.nr==0)&&(A.nc==0))
   return *this;
 if ((nr==0)&&(nc==0)){
   this->init(A,1.,Atrans);
   return *this;
 } 
#if (CONICBUNDLE_DEBUG>=1)
 if (nr!=(Atrans?A.nc:A.nr))
     MEmessage(MEdim(nr,nc,Atrans?A.nc:A.nr,Atrans?A.nr:A.nc,
		     "Matrix::concat_right(const Matrix&) dimensions do not match",MTmatrix));
#endif
 if (nr==0){
   nc+=Atrans?A.nr:A.nc;
   return *this;
 }
 if (mem_dim<nr*nc+A.nr*A.nc){
     Real *mnew;
     mem_dim=Integer(memarray->get(max(2*mem_dim,nr*nc+A.nr*A.nc),mnew));
     if (mem_dim<nr*nc+A.nr*A.nc){
         MEmessage(MEmem(nr*nc+A.nr*A.nc,
                     "Matrix::concat_right(const Matrix&) not enough memory",
                     MTmatrix));
     }
     mat_xey(nr*nc,mnew,m);
     memarray->free(m);
     m=mnew;
 }
 if (Atrans){
   assert(nr==A.nc);
   for (Integer i=0;i<A.nr;i++)
     mat_xey(nr,m+(nc+i)*nr,1,A.m+i,A.nr);
 }
 else
   mat_xey(A.nr*A.nc,m+nr*nc,A.m);
 nc+=Atrans?A.nr:A.nc;
 chk_set_init(*this,1);
 return (*this);
}

Matrix& Matrix::concat_below(const Matrix& A)
{
 chk_init(*this);
 chk_init(A);
 if ((A.nr==0)&&(A.nc==0))
   return *this;
 if ((nr==0)&&(nc==0)){
   *this=A;
   return *this;
 } 
#if (CONICBUNDLE_DEBUG>=1)
 if (nc!=A.nc)
     MEmessage(MEdim(nr,nc,A.nr,A.nc,
                 "Matrix::concat_below(const Matrix&) dimensions do not match",MTmatrix));
#endif
 if (nc==0){
   nr+=A.nr;
   return *this;
 }
 Integer i,j;
 Real* mp;
 Real* ap;
 Real* oldm=m;
 Real* op;
 int free_oldm=0;
 if (mem_dim<nr*nc+A.nr*A.nc){
   mem_dim=Integer(memarray->get(max(2*mem_dim,nr*nc+A.nr*A.nc),m));
     if (mem_dim<nr*nc+A.nr*A.nc){
         MEmessage(MEmem(nr*nc+A.nr*A.nc,
                     "Matrix::concat_right(const Matrix&) not enough memory",
                     MTmatrix));
     }
     free_oldm=1;
 }
 mp=m+nr*nc+A.nr*A.nc-1;
 ap=A.m+A.nr*A.nc-1;
 op=oldm+nr*nc-1;
 for(j=nc;--j>=0;){
     for(i=A.nr;--i>=0;)
         (*mp--)=(*ap--);
     if (mp==op) break;   //only possible for last column on same memory
     for(i=nr;--i>=0;)
         (*mp--)=(*op--);
 }
 nr+=A.nr;
 if (free_oldm) memarray->free(oldm);
 chk_set_init(*this,1);
 return (*this);
}

Matrix& Matrix::concat_below(Real d)
{
 chk_init(*this);
#if (CONICBUNDLE_DEBUG>=1)
 if ((nc!=1)&&(!((nc==0)&&(nr==0))))
     MEmessage(MEdim(nr,nc,1,1,
                 "Matrix::concat_below(Real d) matrix is not a column vector",MTmatrix));
#endif
 Integer n=nr*nc;
 if (mem_dim<n+1){
   Real* mp=m;
   mem_dim=Integer(memarray->get(max(2*mem_dim,n+1),m));
   if (mem_dim<n+1){
     MEmessage(MEmem(n+1,
                     "Matrix::concat_below(Real d) not enough memory",
                     MTmatrix));
   }
   mat_xey(n,m,mp); 
   memarray->free(mp);
 }
 m[n]=d;
 nr=n+1;nc=1;
 return (*this);
}

Matrix& Matrix::concat_right(Real d)
{
 chk_init(*this);
#if (CONICBUNDLE_DEBUG>=1)
 if ((nr!=1)&&(!((nc==0)&&(nr==0))))
     MEmessage(MEdim(nr,nc,1,1,
                 "Matrix::concat_right(Real d) matrix is not a row vector",MTmatrix));
#endif
 Integer n=nr*nc;
 if (mem_dim<n+1){
   Real* mp=m;
   mem_dim=Integer(memarray->get(max(2*mem_dim,n+1),m));
   if (mem_dim<n+1){
     MEmessage(MEmem(n+1,
                     "Matrix::concat_right(Real d) not enough memory",
                     MTmatrix));
   }
   mat_xey(n,m,mp); 
   memarray->free(mp);
 }
 m[n]=d;
 nr=1;
 nc=n+1;
 return (*this);
}

Matrix& Matrix::enlarge_right(Integer addnc)
{
  chk_init(*this);
#if (CONICBUNDLE_DEBUG>=1)
  if (addnc<0)
    MEmessage(MEdim(nr,nc,0,addnc,
		    "Matrix::enlarge_right(Integer addnc) addnc <0",MTmatrix));
#endif
  if (addnc<=0) 
    return *this;

  if (nr==0){
    nc+=addnc;
    return *this;
  } 
  if (mem_dim<nr*(nc+addnc)){
     Real *mnew;
     mem_dim=Integer(memarray->get(max(2*mem_dim,nr*(nc+addnc)),mnew));
     if (mem_dim<nr*(nc+addnc)){
       MEmessage(MEmem(nr*(nc+addnc),
                     "Matrix::enlarge_right(Integer addnc) not enough memory",
                     MTmatrix));
     }
     mat_xey(nr*nc,mnew,m);
     memarray->free(m);
     m=mnew;
  }
  //mat_xey(A.nr*A.nc,m+nr*nc,A.m);
  nc+=addnc;
  chk_set_init(*this,0);
  return (*this);
}

Matrix& Matrix::enlarge_below(Integer addnr)
{
  chk_init(*this);
#if (CONICBUNDLE_DEBUG>=1)
  if (addnr<0)
    MEmessage(MEdim(nr,nc,addnr,0,
                 "Matrix::enlarge_below(Integer addnr) addnr <0",MTmatrix));
#endif
  if (addnr<=0)
    return *this;
  if (nc==0){
    nr+=addnr;
    return *this;
  }
  Integer i,j;
  Real* mp;
  Real* oldm=m;
  Real* op;
  int free_oldm=0;
  if (mem_dim<(nr+addnr)*nc){
    mem_dim=Integer(memarray->get(max(2*mem_dim,(nr+addnr)*nc),m));
    if (mem_dim<(nr+addnr)*nc){
      MEmessage(MEmem((nr+addnr)*nc,
                     "Matrix::enlarge_below(Integer addnr) not enough memory",
                     MTmatrix));
    }
    free_oldm=1;
  }
  mp=m+(nr+addnr)*nc-1;
  op=oldm+nr*nc-1;
  for(j=nc;--j>=0;){
    mp-=addnr;
    if (mp==op) 
      break;   //only possible for last column on same memory
    for(i=nr;--i>=0;)
      (*mp--)=(*op--);
  }
  nr+=addnr;
  if (free_oldm) memarray->free(oldm);
  chk_set_init(*this,0);
  return (*this);
}

Matrix& Matrix::enlarge_right(Integer addnc,Real d)
{
  chk_init(*this);
#if (CONICBUNDLE_DEBUG>=1)
  if (addnc<0)
    MEmessage(MEdim(nr,nc,0,addnc,
		    "Matrix::enlarge_right(Integer addnc,Real d) addnc <0",MTmatrix));
#endif
  if (addnc<=0) 
    return *this;

  if (nr==0){
    nc+=addnc;
    return *this;
  } 
  if (mem_dim<nr*(nc+addnc)){
     Real *mnew;
     mem_dim=Integer(memarray->get(max(2*mem_dim,nr*(nc+addnc)),mnew));
     if (mem_dim<nr*(nc+addnc)){
       MEmessage(MEmem(nr*(nc+addnc),
                     "Matrix::enlarge_right(Integer addnc,Real d) not enough memory",
                     MTmatrix));
     }
     mat_xey(nr*nc,mnew,m);
     memarray->free(m);
     m=mnew;
  }
  mat_xea(nr*addnc,m+nr*nc,d);
  nc+=addnc;
  return (*this);
}

Matrix& Matrix::enlarge_below(Integer addnr,Real d)
{
  chk_init(*this);
#if (CONICBUNDLE_DEBUG>=1)
  if (addnr<0)
    MEmessage(MEdim(nr,nc,addnr,0,
                 "Matrix::enlarge_below(Integer addnr,Real d) addnr <0",MTmatrix));
#endif
  if (addnr<=0)
    return *this;
  if (nc==0){
    nr+=addnr;
    return *this;
  }
  Integer i,j;
  Real* mp;
  Real* oldm=m;
  Real* op;
  int free_oldm=0;
  if (mem_dim<(nr+addnr)*nc){
    mem_dim=Integer(memarray->get(max(2*mem_dim,(nr+addnr)*nc),m));
    if (mem_dim<(nr+addnr)*nc){
      MEmessage(MEmem((nr+addnr)*nc,
                     "Matrix::enlarge_below(Integer addnr,Real d) not enough memory",
                     MTmatrix));
    }
    free_oldm=1;
  }
  mp=m+(nr+addnr)*nc-1;
  op=oldm+nr*nc-1;
  for(j=nc;--j>=0;){
    for(i=addnr;--i>=0;)
          (*mp--)=d;
    if (mp==op) 
      break;   //only possible for last column on same memory
    for(i=nr;--i>=0;)
      (*mp--)=(*op--);
  }
  nr+=addnr;
  if (free_oldm) memarray->free(oldm);
  return (*this);
}

  Matrix& Matrix::enlarge_right(Integer addnc,const Real* dp,Real d)
{
  chk_init(*this);
#if (CONICBUNDLE_DEBUG>=1)
  if (addnc<0)
    MEmessage(MEdim(nr,nc,0,addnc,
		    "Matrix::enlarge_right(Integer addnc,const Real* dp,Real d) addnc <0",MTmatrix));
#endif
  if (addnc<=0) 
    return *this;

  if (nr==0){
    nc+=addnc;
    return *this;
  } 
  if (mem_dim<nr*(nc+addnc)){
     Real *mnew;
     mem_dim=Integer(memarray->get(max(2*mem_dim,nr*(nc+addnc)),mnew));
     if (mem_dim<nr*(nc+addnc)){
       MEmessage(MEmem(nr*(nc+addnc),
                     "Matrix::enlarge_right(Integer addnc,const Real* dp,Real d) not enough memory",
                     MTmatrix));
     }
     mat_xey(nr*nc,mnew,m);
     memarray->free(m);
     m=mnew;
  }
  mat_xeya(nr*addnc,m+nr*nc,dp,d);
  nc+=addnc;
  return (*this);
}

  Matrix& Matrix::enlarge_below(Integer addnr,const Real* dp,Real d)
{
  chk_init(*this);
#if (CONICBUNDLE_DEBUG>=1)
  if (addnr<0)
    MEmessage(MEdim(nr,nc,addnr,0,
                 "Matrix::enlarge_below(Integer addnr,const Real* dp,Real d) addnr <0",MTmatrix));
#endif
  if (addnr<=0)
    return *this;
  if (nc==0){
    nr+=addnr;
    return *this;
  }
  Integer i,j;
  Real* mp;
  Real* oldm=m;
  Real* op;
  int free_oldm=0;
  if (mem_dim<(nr+addnr)*nc){
    mem_dim=Integer(memarray->get(max(2*mem_dim,(nr+addnr)*nc),m));
    if (mem_dim<(nr+addnr)*nc){
      MEmessage(MEmem((nr+addnr)*nc,
                     "Matrix::enlarge_below(Integer addnr,const Real* dp,Real d) not enough memory",
                     MTmatrix));
    }
    free_oldm=1;
  }
  mp=m+(nr+addnr)*nc-1;
  dp+=addnr*nc-1;
  op=oldm+nr*nc-1;
  for(j=nc;--j>=0;){
    for(i=addnr;--i>=0;)
      (*mp--)=(*dp--)*d;
    if (mp==op) 
      break;   //only possible for last column on same memory
    for(i=nr;--i>=0;)
      (*mp--)=(*op--);
  }
  nr+=addnr;
  if (free_oldm) memarray->free(oldm);
  return (*this);
}

Matrix& Matrix::rand(Integer inr,Integer inc,GB_rand* rgp)
{
  if (rgp==0) 
    rgp=&mat_randgen;
 newsize(inr,inc);
 Real *mp=m;
 for(int i=inr*inc;--i>=0;)
     (*mp++)=rgp->next();
 chk_set_init(*this,1);
 return *this;
}

Matrix& Matrix::rand_normal(Integer inr,Integer inc,Real mean,Real variance,int generator_type)
{
 newsize(inr,inc);
 Real *mp=m;
 std::normal_distribution<double> distrib(mean,variance);
 switch(generator_type){
 case 0: default: {
   for(int i=inr*inc;--i>=0;)
     (*mp++)=distrib(mat_std_randgen);
   break;
 }
 case 1: {
   for(int i=inr*inc;--i>=0;)
     (*mp++)=distrib(mat_mt_randgen);
   break;
 }
 case 2: {
   for(int i=inr*inc;--i>=0;)
     (*mp++)=distrib(mat_mt64_randgen);
   break;
 }
 }
 chk_set_init(*this,1);
 return *this;
}

Matrix& Matrix::shuffle(CH_Tools::GB_rand* rgp)
{
  if (rgp==0)
    rgp=&mat_randgen;
  Integer dim=nr*nc;
  for(Integer i=0;i<dim;i++){
    Integer h=Integer(rgp->unif_long(dim-i));
    Real d=m[i+h];
    m[i+h]=m[i];
    m[i]=d;
  }
  return *this;
} 

Matrix& Matrix::inv(void)
{
 chk_init(*this);
 Integer i;
 for(i=nr*nc;--i>=0;) m[i]=1./m[i];
 return *this;
}

Matrix& Matrix::sqrt(void)
{
 chk_init(*this);
 Integer i;
 for(i=nr*nc;--i>=0;)
   m[i]=::sqrt(m[i]);
 return *this;
}

Matrix& Matrix::sqr(void)
{
 chk_init(*this);
 Integer i;
 for(i=nr*nc;--i>=0;)
   m[i]*=m[i];
 return *this;
}

Matrix& Matrix::sign(Real tol)
{
 chk_init(*this);
 Integer i;
 Real *mp=m;
 for(i=nr*nc;--i>=0;) {
     if (*mp>tol) *mp++=1.;
     else if (*mp<-tol) *mp++=-1.;
     else *mp++=0.;
 }
 return *this;
}

Matrix& Matrix::floor(void)
{
 chk_init(*this);
 for(Integer i=nr*nc;--i>=0;) m[i]=::floor(m[i]);
 return *this;
}

Matrix& Matrix::ceil(void)
{
 chk_init(*this);
 for(Integer i=nr*nc;--i>=0;) m[i]=::ceil(m[i]);
 return *this;
}

Matrix& Matrix::rint(void)
{
 chk_init(*this);
#ifndef __unix
 for(Integer i=nr*nc;--i>=0;) m[i]=::floor(m[i]+.5);
#else
 for(Integer i=nr*nc;--i>=0;) m[i]=::rint(m[i]);
#endif
 return *this;
}

Matrix& Matrix::round(void)
{
 chk_init(*this);
#ifndef __unix
 for(Integer i=nr*nc;--i>=0;) m[i]=::floor(m[i]+.5);
#else
 for(Integer i=nr*nc;--i>=0;) m[i]=::round(m[i]);
#endif
 return *this;
}

Matrix& Matrix::abs(void)
{
 chk_init(*this);
 for(Integer i=0;i<nr*nc;i++) m[i]=::fabs(m[i]);
 return *this;
}


bool Matrix::contains_nan()
{
  Integer i=nr*nc;
  const Real* mp=m;
  while ((--i>=0)&&(!isnan(*mp++)));
  return (i>=0);
}


Matrix& Matrix::scale_rows(const Matrix& vec)
{
 chk_init(*this);
 chk_init(vec);
#if (CONICBUNDLE_DEBUG>=1)
 if (vec.dim()!=nr) MEmessage(MEdim(nr,nc,vec.dim(),1,
                    "Matrix::scale_rows(const Matrix& vec) vec.dim() is not number of  rows",
                    MTmatrix));
#endif
 const Real *vp=vec.m;
 Integer i;
 for(i=0;i<nr;i++,vp++){
     mat_xmultea(nc,m+i,nr,*vp);
 }
 return *this;
}

Matrix& Matrix::scale_cols(const Matrix& vec)
{
 chk_init(*this);
 chk_init(vec);
#if (CONICBUNDLE_DEBUG>=1)
 if (vec.dim()!=nc) MEmessage(MEdim(nr,nc,vec.dim(),1,
                    "Matrix::scale_cols(const Matrix& vec) vec.dim() is not number of  columns",
                    MTmatrix));
#endif
 const Real *vp=vec.m;
 Integer i;
 for(i=0;i<nc;i++,vp++){
     mat_xmultea(nr,m+i*nr,*vp);
 }
 return *this;
}

Indexmatrix Matrix::find(Real tol) const
{
 chk_init(*this);
 Indexmatrix ind(nr*nc,1);
 chk_set_init(ind,1);
 Integer i,k=0;
 Real *mp=m;
 for(i=0;i<nr*nc;i++){
     if (::fabs(*mp++)>tol) ind(k++)=i;
 }
 ind.nr=k;
 return ind;
}

Indexmatrix Matrix::find_number(Real num,Real tol) const
{
 chk_init(*this);
 Indexmatrix ind(nr*nc,1);
 chk_set_init(ind,1);
 Integer i,k=0;
 Real *mp=m;
 for(i=0;i<nr*nc;i++){
     if (::fabs(*mp++-num)<tol) ind(k++)=i;
 }
 ind.nr=k;
 return ind;
}

  Matrix& Matrix::delete_rows(const Indexmatrix& ind,bool sorted_increasingly)
{
 chk_init(*this);
 chk_init(ind);
 if (ind.dim()==0) return *this;
 Indexmatrix sind;
 const Integer *ip;
 if (sorted_increasingly){
   ip=ind.get_store();
 }
 else {
   sortindex(ind,sind);
   sind=ind(sind);
   ip=sind.get_store();
 }

 chk_range(ip[0],ip[ind.dim()-1],nr,nr);
 Integer j,k;
 Real* mp=m;
 for(j=0;(j<nc);j++){
     mat_xey(ip[0],mp,m+j*nr);
     mp+=ip[0];
     for(k=1;k<ind.dim();k++){
#if (CONICBUNDLE_DEBUG>=1)
         if (ip[k]<=ip[k-1]) MEmessage(MatrixError(ME_unspec,"Matrix::delete_rows(): a row is to be deleted twice or sorted_increasingly is incorrectly set to true",MTmatrix));
#endif
         mat_xey(ip[k]-ip[k-1]-1,mp,m+j*nr+ip[k-1]+1);
         mp+=ip[k]-ip[k-1]-1;
     }
     mat_xey(nr-ip[k-1]-1,mp,m+j*nr+ip[k-1]+1);
     mp+=nr-ip[k-1]-1;
 }
 nr-=ind.dim();
 return *this;
}


  Matrix& Matrix::delete_cols(const Indexmatrix& ind,bool sorted_increasingly)
{
 chk_init(*this);
 chk_init(ind);
 if (ind.dim()==0) return *this;
 Indexmatrix sind;
 const Integer* ip;
 if (sorted_increasingly){
   ip=ind.get_store();
 }
 else {
   sortindex(ind,sind);
   sind=ind(sind);
   ip=sind.get_store();
 }
 chk_range(ip[0],ip[ind.dim()-1],nc,nc);
 Integer j,oldj,k;
 oldj=ip[0];
 k=1;
 for(j=ip[0]+1;(j<nc)&&(k<ind.dim());j++){

#if (CONICBUNDLE_DEBUG>=1)
     if (ip[k]<j) MEmessage(MatrixError(ME_unspec,"Matrix::delete_cols(): a column is to be deleted twice",MTmatrix));
#endif

     if (j==ip[k]){
         k++;
         continue;
     }

     mat_xey(nr,m+oldj*nr,m+j*nr);
     oldj++;
 }

 mat_xey(nr*(nc-j),m+oldj*nr,m+j*nr);
 nc-=ind.dim();
 return *this;
}

Matrix&  Matrix::insert_row(Integer ind,const Matrix& v)
{
  chk_init(*this);
  chk_init(v);
#if (CONICBUNDLE_DEBUG>=1)
  if ((nc!=v.dim())&&(!((nc==0)&&(nr==0)))) 
    MEmessage(MatrixError(ME_unspec,"Matrix::insert_row(): column dimensions do not match",MTmatrix));
  if ((ind<0)||(ind>nr)) MEmessage(MatrixError(ME_unspec,"Matrix::insert_row(): index out of range",MTmatrix));
#endif
  nc=v.dim();
  Real* mp;
  Real* oldm=m;
  Real* op;
  int free_oldm=0;
  if (mem_dim<(nr+1)*nc){
    mem_dim=Integer(memarray->get((nr+1)*nc,m));
    if (mem_dim<(nr+1)*nc){
      MEmessage(MEmem((nr+1)*nc,
		      "Matrix::insert_row(): not enough memory",
		      MTmatrix));
    }
    free_oldm=1;
  }
  mp=m+(nr+1)*nc-1;
  op=oldm+nr*nc-1;
  mat_xey(nr-ind,mp,-1,op,-1);
  mp-=nr-ind+1;
  op-=nr-ind;
  for(Integer j=1;j<nc;j++){
    mat_xey(nr,mp,-1,op,-1);
    mp-=nr+1;
    op-=nr;
  }
  if (free_oldm){
    mat_xey(ind,m,oldm);
    memarray->free(oldm);
  }
  mat_xey(nc,m+ind,nr+1,v.get_store(),1);
  nr++;
  return *this;
}

Matrix&  Matrix::insert_col(Integer ind,const Matrix& v)
{
  chk_init(*this);
  chk_init(v);
#if (CONICBUNDLE_DEBUG>=1)
  if ((nr!=v.dim())&&(!((nc==0)&&(nr==0)))) 
    MEmessage(MatrixError(ME_unspec,"Matrix::insert_col(): row dimensions do not match",MTmatrix));
  if ((ind<0)||(ind>nc)) MEmessage(MatrixError(ME_unspec,"Matrix::insert_col(): index out of range",MTmatrix));
#endif
  nr=v.dim();
  Real* oldm=m;
  int free_oldm=0;
  if (mem_dim<nr*(nc+1)){
    mem_dim=Integer(memarray->get(nr*(nc+1),m));
    if (mem_dim<nr*(nc+1)){
      MEmessage(MEmem(nr*(nc+1),
		      "Matrix::insert_col(): not enough memory",
		      MTmatrix));
    }
    free_oldm=1;
    mat_xey(nr*ind,m,oldm);
  }
  mat_xey(nr*(nc-ind),m+nr*(nc+1)-1,-1,oldm+nr*nc-1,-1);
  mat_xey(nr,m+ind*nr,v.get_store());
  if (free_oldm){
    memarray->free(oldm);
  }
  nc++;
  return *this;
}
   

void Matrix::display(std::ostream& out,int precision,int width,
                        int screenwidth) const
{
 out<<"Matrix("<<nr<<","<<nc<<")"<<std::endl;
 if ((nr==0)||(nc==0)) return;
 chk_init(*this);
 if (precision==0) precision=4;
 out.precision(precision);
 if (width==0) width=precision+6;
 if (screenwidth==0) screenwidth=80;
 Integer colnr=screenwidth/(width+1);
 Integer k,i,j;
 Integer maxk=nc/colnr+((nc%colnr)>0);
 Integer maxj;
 for(k=0;k<maxk;k++){
     out<<"columns "<<k*colnr<<" to "<<min(nc,(k+1)*colnr)-1<<std::endl;
     for(i=0;i<nr;i++){
         maxj=min((k+1)*colnr,nc);
         for(j=k*colnr;j<maxj;j++){
             out<<' ';out.width(width);out<<(*this)(i,j);
         }
         out<<std::endl;
     }     
 }
}

void Matrix::mfile_output(std::ostream& out, ///< output stream
			  int precision,   ///< number of most significant digits, default=16
			  int width      ///< field width, default = precision+6
			  ) const
{
  chk_init(*this);
  if (precision<=0) 
    precision=16;
  std::streamsize old_prec=out.precision();
  if (width<=0) 
    width=precision+6;
  out<<"[";
  out.precision(precision);
  for(Integer i=0;i<nr;i++){
    for(Integer j=0;j<nc;j++){
      out<<" ";
      out.width(width);
      out<<(*this)(i,j);
    }
    if (i<nr-1)
      out<<"\n";
  }
  out<<"];\n";
  out.precision(old_prec);
}


Matrix diag(const Matrix& A)
{
 chk_init(A);
 if (min(A.nr,A.nc)==1){ //make a diagonal matrix of this vector
     Integer k=max(A.nr,A.nc);
     Matrix B(k,k,0.);
     for(Integer i=0;i<k;i++)
         B(i,i)=A(i);
     return B;
 }
 return Matrix(min(A.nr,A.nc),1,A.m,A.nr+1);
}

Matrix sumrows(const Matrix& A)
{
 chk_init(A);
 Matrix v(1,A.nc,0.);
 for(Integer i=0;i<A.nr;i++)
     mat_xpey(A.nc,v.m,1,A.m+i,A.nr);
 chk_set_init(v,1);
 return v;
}

Matrix sumcols(const Matrix& A)
{
 chk_init(A);
 Matrix v(A.nr,1,0.);
 for(Integer i=0;i<A.nc;i++)
     mat_xpey(A.nr,v.m,A.m+i*A.nr);
 chk_set_init(v,1);
 return v;
}

Real sum(const Matrix& A)
{
 chk_init(A);
 Real s=0.;
 Integer i;
 Real *mp=A.m;
 for(i=A.nr*A.nc;--i>=0;){
     s+=(*mp++);
 }
 return s;
}

Real ip_min_max(const Matrix& A, const Matrix& B,Real &minval,Real& maxval)
{
  chk_add(A,B);
  Real sumval=0.;
  minval=max_Real;
  maxval=min_Real;
  const Real* ap=A.m;
  const Real* const apend=ap+A.nr*A.nc;
  const Real* bp=B.m;
  while(ap!=apend){
    Real d= *ap++ * *bp++;
    sumval+=d;
    if (d<minval)
      minval=d;
    if (d>maxval)
      maxval=d;
  }
  return sumval;
 }

Matrix maxrows(const Matrix& A)
{
 chk_init(A);
 if (A.dim()==0) return Matrix(0,0,0.);
 Matrix v(1,A.nc);
 Real maxd;
 Integer i,j;
 for(j=0;j<A.nc;j++){
     maxd=A(0,j);
     for(i=1;i<A.nr;i++)
         maxd=max(maxd,A(i,j));
     v(j)=maxd;
 }
 chk_set_init(v,1);
 return v;
}

Matrix maxcols(const Matrix& A)
{
 chk_init(A);
 if (A.dim()==0) return Matrix(0,0,0.);
 Matrix v(A.nr,1);
 Real maxd;
 Integer i,j;
 for(i=0;i<A.nr;i++){
     maxd=A(i,0);
     for(j=1;j<A.nc;j++)
         maxd=max(maxd,A(i,j));
     v(i)=maxd;
 }
 chk_set_init(v,1);
 return v;
}

Real max(const Matrix& A,Integer *iindex,Integer *jindex)
{
 chk_init(A);
 if (A.dim()==0) return min_Real;
 Real maxd;
 if (iindex==0){
     Integer i=A.nr*A.nc-1;
     const Real* ap=A.get_store();
     maxd=*ap++;
     for(;--i>=0;){
         maxd=max(maxd,*ap++);
     }
 }
 else{
     Integer dim=A.nr*A.nc-1;
     const Real* ap=A.get_store();
     maxd=*ap;
     Integer besti=0;
     for(Integer i=dim;--i>=0;){
         if (*(++ap)>maxd){
             maxd=*ap;
             besti=dim-i;
         }
     }
     if (jindex!=0){
         *jindex=besti/A.nr;
         *iindex=besti%A.nr;
     }
     else *iindex=besti;
 }
     
 return maxd;
}

Matrix minrows(const Matrix& A)
{
 chk_init(A);
 if (A.dim()==0) return Matrix(0,0,0.);
#if (CONICBUNDLE_DEBUG>=1)
 if ((A.nc==0)||(A.nr==0))
     MEmessage(MEdim(A.nc,A.nr,0,0,"minrows(const Matrix&) dimension zero",MTmatrix));
#endif
 Matrix v(1,A.nc);
 Real mind;
 Integer i,j;
 for(j=0;j<A.nc;j++){
     mind=A(0,j);
     for(i=1;i<A.nr;i++)
         mind=min(mind,A(i,j));
     v(j)=mind;
 }
 chk_set_init(v,1);
 return v;
}

Matrix mincols(const Matrix& A)
{
 chk_init(A);
 if (A.dim()==0) return Matrix(0,0,0.);
 Matrix v(A.nr,1);
 Real mind;
 Integer i,j;
 for(i=0;i<A.nr;i++){
     mind=A(i,0);
     for(j=1;j<A.nc;j++)
         mind=min(mind,A(i,j));
     v(i)=mind;
 }
 chk_set_init(v,1);
 return v;
}

Real min(const Matrix& A,Integer *iindex,Integer *jindex)
{
 chk_init(A);
 if (A.dim()==0) return max_Real;
 Real mind;
 if (iindex==0){
     Integer i=A.nr*A.nc-1;
     const Real* ap=A.get_store();
     mind=*ap++;
     for(;--i>=0;){
         mind=min(mind,*ap++);
     }
 }
 else{
     Integer dim=A.nr*A.nc-1;
     const Real* ap=A.get_store();
     mind=*ap;
     Integer besti=0;
     for(Integer i=dim;--i>=0;){
         if (*(++ap)<mind){
             mind=*ap;
             besti=dim-i;
         }
     }
     if (jindex!=0){
         *jindex=besti/A.nr;
         *iindex=besti%A.nr;
     }
     else *iindex=besti;
 }
     
 return mind;
}

void sortindex(const Matrix& vec,Indexmatrix& ind,bool nondecreasing)
{
  ind.init(Range(0,vec.dim()-1));
  Integer *ip=ind.get_store();
  const Real *vp=vec.get_store();
  //heapsort<Integer,Real>(vec.dim(),ip,vp,nondecreasing);
  if (nondecreasing)
    std::sort( ip, ip + vec.dim(), mat_less_index<Real>(vp) );
  else
    std::sort( ip, ip + vec.dim(), mat_greater_index<Real>(vp) );
}


Real trace(const Matrix& A)
{
 chk_init(A);
 Real sum=0.;
 Integer k=min(A.nr,A.nc);
 for(Integer i=0;i<k;i++) sum+=A(i,i);
 return sum;
}

///returns trace(A^TDA)=\|A\|^2_D with D=Diag(d). A may be transposed, D may be inverted 
Real normDsquared(const Matrix& A,const Matrix& d,int atrans,int dinv)
{
  chk_init(A);
  chk_init(d);
  const Real *ap=A.m;
  const Real *dp=d.m;
  Real ipsum=0.;

  if (atrans==0){
#if (CONICBUNDLE_DEBUG>=1)
    if (d.dim()!=A.nr){
    MEmessage(MatrixError(ME_dim,"normDsquared: dimensions don't match",MTmatrix));
  }
#endif 
    if (dinv==0){
      for(Integer j=A.nc;--j>=0;){
	const Real* ddp=dp;
	for (Integer i=A.nr;--i>=0;ap++){
	  ipsum+=(*ddp++)*(*ap)*(*ap);
	}
      }
    }
    else {
      for(Integer j=A.nc;--j>=0;){
	const Real* ddp=dp;
	for (Integer i=A.nr;--i>=0;ap++){
	  ipsum+=(*ap)*(*ap)/(*ddp++);
	}
      }
    }
  }
  else {
#if (CONICBUNDLE_DEBUG>=1)
    if (d.dim()!=A.nc){
      MEmessage(MatrixError(ME_dim,"normDsquared: dimensions don't match",MTmatrix));
    }
#endif 
    if (dinv==0){
      for(Integer j=A.nc;--j>=0;){
	Real colsum=0.;
	for (Integer i=A.nr;--i>=0;ap++){
	  colsum+=(*ap)*(*ap);
	}
	ipsum+=colsum*(*dp++);
      }
    }
    else {
      for(Integer j=A.nc;--j>=0;){
	Real colsum=0.;
	for (Integer i=A.nr;--i>=0;ap++){
	  colsum+=(*ap)*(*ap);
	}
	ipsum+=colsum/(*dp++);
      }
    }
  }
  return ipsum;
}



Matrix abs(const Matrix& A)
{
 chk_init(A);
 Matrix B(A.nr,A.nc);
 const Real *ap=A.get_store();
 Real *bp=B.get_store();
 for(Integer i=A.nr*A.nc;--i>=0;) *bp++=fabs(*ap++);
 chk_set_init(B,1);
 return B;
}

Matrix transpose(const Matrix& A)
{
 chk_init(A);
 if ((A.nr<=1)||(A.nc<=1)) return Matrix(A.nc,A.nr,A.m);
 Matrix B(A.nc,A.nr);
 Integer i;
 for(i=0;i<B.nc;i++){
     mat_xey(B.nr,B.m+i*B.nr,1,A.m+i,A.nr);
 }
 chk_set_init(B,1);
 return B;
}

std::ostream& operator<<(std::ostream& o,const Matrix &A)
{
 chk_init(A);
 o<<A.nr<<" "<<A.nc<<'\n';
 Integer i,j;
 for(i=0;i<A.nr;i++){
     for(j=0;j<A.nc;j++) o<<' '<<A(i,j);
     o<<'\n';
 }
 return o;
}

std::istream& operator>>(std::istream& in,Matrix &A)
{
 Real d;
 Integer nr,nc;
 in>>d;
 nr=Integer(d+.5);
 in>>d;
 nc=Integer(d+.5);
 if ((nr<0)||(nc<0))
     MEmessage(MEdim(nc,nr,0,0,"operator>>(std::istream&,Matrix&) dimension negative",MTmatrix));
 A.newsize(nr,nc);
 Integer i,j;
 for(i=0;i<nr;i++)
     for(j=0;j<nc;j++)
         in>>A(i,j);
 chk_set_init(A,1);
 return in;
}

// **************************************************************************
//                     Matrix specific Indexmatrix 
// **************************************************************************

Indexmatrix::Indexmatrix(const Matrix &A)
{
 chk_init(A);
 init_to_zero();
 newsize(A.nr,A.nc);
 Integer i;
 Integer *mp=m;
 Real *matp=A.m;
 for(i=nr*nc;--i>=0;matp++){
     *mp++=Integer((*matp>0)? *matp+.5 : *matp-.5);
 }
 chk_set_init(*this,1);
}

}

