/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  Matrix/matop.hxx
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



#ifndef CH_MATRIX_CLASSES__MATOP_HXX
#define CH_MATRIX_CLASSES__MATOP_HXX

/**  @file matop.hxx
    @brief Basic types, constants and templates for simple linear algebra routines like BLAS level 1.

    @version 1.0
    @date 2005-03-01
    @author Christoph Helmberg

*/


#include <limits.h>
#include <float.h>
#include <assert.h>
#if defined(WITH_BLAS) || defined(BLAS_GENMULT)
extern "C"{
#include "cblas.h"
}
#endif
#include <random>
#include "gb_rand.hxx"
#include "CBconfig.hxx"

/**@brief  Matrix Classes and Linear Algebra. See \ref matrixclasses for a quick introduction.
*/	
namespace CH_Matrix_Classes {


  //---------------------------------------------------------------------------

  /**@defgroup matop_types Basic Types and Constants 
   */
  //@{
  
  /// all integer numbers in calculations and indexing are of this type 
  typedef int Integer; 

  /// maximal attainable value by an #Integer
  const Integer max_Integer= INT_MAX; 

  /// minimal attainable value by an #Integer  
  const Integer min_Integer= INT_MIN;

  
  /// all real numbers in calculations are of this type 
  typedef double Real;
  
  /// maximal attainable value by a #Real
  const Real max_Real= DBL_MAX;

  /// minimal attainable value by a #Real
  const Real min_Real= -DBL_MAX;
  
  /// machine epsilon for type #Real
  const Real eps_Real= DBL_EPSILON;

  //@}

  //---------------------------------------------------------------------------

  /**@defgroup matop_templates Basic Templates for Linear Algebra
     @brief templates for simple linear algebra routines like BLAS level 1.
   */
  //@{


  /**@brief Set x[i]=a for len elements of the array x.
     @param len gives the number of elements to be processed
     @param x points to the starting address 
     @param a is the scalar value
   */
template<class Val>
inline void mat_xea( Integer len, Val *x, const Val a)
{
  const Val* const xend=x+len;
  for(;x!=xend;) 
    (*x++)=a;
}

  /**@brief Set x[i]=a for len elements of the array x incremented by incx.
     @param len gives the number of elements to be processed
     @param x points to the starting address 
     @param incx gives the increment for x in each step
     @param a is the scalar value
   */
template<class Val>
inline void mat_xea( Integer len, Val *x, const Integer incx,
                     const Val a)
{
 const Val* const xend=x+len*incx;
 for(;x!=xend;x+=incx) 
   (*x)=a;
}

#ifdef WITH_BLAS
  /**@brief Copy a Real array of length len to destination x from source y.
     @param len gives the number of elements to be processed
     @param x points to the starting address of the destination
     @param y points to the starting address of the source
   */
  inline void mat_xey( Integer len, double *x, const double* y)
  {
    cblas_dcopy(len,y,1,x,1);
  }
#endif

  /**@brief Copy an array of length len to destination x from source y.
     @param len gives the number of elements to be processed
     @param x points to the starting address of the destination
     @param y points to the starting address of the source
   */
template<class Val>
inline void mat_xey( Integer len, Val *x, const Val* y)
{
  const Val* const xend=x+len;
  for(;x!=xend;) 
    (*x++)=(*y++);
}

#ifdef WITH_BLAS
  /**@brief Copy a Real array of length len to destination x (increment incx) from source y (increment incy).
     @param len gives the number of elements to be processed
     @param x points to the starting address of the destination
     @param incx gives the increment for x in each step
     @param y points to the starting address of the source
     @param incy gives the increment for y in each step
   */
 inline void mat_xey( Integer len, double *x, const Integer incx, 
		      const double* y, const Integer incy)
 {
   cblas_dcopy(len,y,incy,x,incx);
 }
#endif

  /**@brief Copy an array of length len to destination x (increment incx) from source y (increment incy).
     @param len gives the number of elements to be processed
     @param x points to the starting address of the destination
     @param incx gives the increment for x in each step
     @param y points to the starting address of the source
     @param incy gives the increment for y in each step
   */
template<class Val>
inline void mat_xey( Integer len, Val *x, const Integer incx,
                     const Val *y, const Integer incy)
{
  const Val* const xend=x+len*incx;
  for(;x!=xend;x+=incx,y+=incy) 
    (*x)=(*y);
}

#ifdef WITH_BLAS
  /**@brief Set x[i]=-y[i] for len elements of the arrays x and y.
     @param len gives the number of elements to be processed
     @param x points to the starting address of the destination
     @param y points to the starting address of the source
   */
inline void mat_xmey( Integer len, double *x,
		      const double *y)
 {
   cblas_daxpy(len,-1.,y,1,x,1);
 }
#endif

  /**@brief Set x[i]=-y[i] for len elements of the arrays x and y.
     @param len gives the number of elements to be processed
     @param x points to the starting address of the destination
     @param y points to the starting address of the source
   */
template<class Val>
inline void mat_xmey( Integer len, Val *x, const Val *y)
{
  const Val* const xend=x+len;
  for(;x!=xend;) 
    (*x++)-=(*y++);
}

#ifdef WITH_BLAS
  /**@brief Set x[i]=-y[i] for len elements of the arrays x and y incremented by incx and incy, respectively.
     @param len gives the number of elements to be processed
     @param x points to the starting address of the destination
     @param incx gives the increment for x in each step
     @param y points to the starting address of the source
     @param incy gives the increment for y in each step
   */
inline void mat_xmey( Integer len, double *x, const Integer incx,
		      const double *y, const Integer incy)
 {
   cblas_daxpy(len,-1.,y,incy,x,incx);
 }
#endif

  /**@brief Set x[i]=-y[i] for len elements of the arrays x and y incremented by incx and incy, respectively.
     @param len gives the number of elements to be processed
     @param x points to the starting address of the destination
     @param incx gives the increment for x in each step
     @param y points to the starting address of the source
     @param incy gives the increment for y in each step
   */
template<class Val>
inline void mat_xmey( Integer len, Val *x, const Integer incx,
                          const Val *y, const Integer incy)
{
  const Val* const xend=x+len*incx;
  for(;x!=xend;x+=incx,y+=incy) 
   (*x)-=(*y);
}

#ifdef WITH_BLAS
  /**@brief Set x[i]=-x[i] for an array of length len.
     @param len gives the number of elements
     @param x points to the starting address
   */
inline void mat_xemx( Integer len, double *x)
 {
   cblas_dscal(len,-1.,x,1);
 }
#endif

  /**@brief Set x[i]=-x[i] for an array of length len.
     @param len gives the number of elements
     @param x points to the starting address
   */
template<class Val>
inline void mat_xemx( Integer len, Val *x)
{
  const Val* const xend=x+len;
  for(;x!=xend;x++) 
     (*x)=-(*x);
}

#ifdef WITH_BLAS
  /**@brief Set x[i]=-x[i] for len elements of an array incremented by incx.
     @param len gives the number of elements
     @param x points to the starting address
     @param incx gives the increment for x in each step
   */
inline void mat_xemx( Integer len, double *x, const Integer incx)
 {
   cblas_dscal(len,-1.,x,incx);
 }
#endif

  /**@brief Set x[i]=-x[i] for len elements of an array incremented by incx.
     @param len gives the number of elements
     @param x points to the starting address
     @param incx gives the increment for x in each step
   */
template<class Val>
inline void mat_xemx( Integer len, Val *x, const Integer incx)
{
  const Val* const xend=x+len*incx;
  for(;x!=xend;x+=incx) 
     (*x)=-(*x);
}

#ifdef WITH_BLAS
  /**@brief Set x[i]=-y[i] for len elements of the arrays x and y.
     @param len gives the number of elements to be processed
     @param x points to the starting address of the destination
     @param y points to the starting address of the source
   */
inline void mat_xemy( Integer len, double *x, const double *y)
 {
   cblas_dcopy(len,y,1,x,1);
   cblas_dscal(len,-1.,x,1);
 }
#endif

  /**@brief Set x[i]=-y[i] for len elements of the arrays x and y.
     @param len gives the number of elements to be processed
     @param x points to the starting address of the destination
     @param y points to the starting address of the source
   */
template<class Val>
inline void mat_xemy( Integer len, Val *x, const Val *y)
{
  const Val* const xend=x+len;
  for(;x!=xend;)
    (*x++)=-(*y++);
}

#ifdef WITH_BLAS
  /**@brief Set x[i]=-y[i] for len elements of the arrays x and y incremented by incx and incy, respectively.
     @param len gives the number of elements to be processed
     @param x points to the starting address of the destination
     @param incx gives the increment for x in each step
     @param y points to the starting address of the source
     @param incy gives the increment for y in each step
   */
inline void mat_xemy( Integer len, double *x, const Integer incx,
		      const double *y, const Integer incy)
 {
   cblas_dcopy(len,y,incy,x,incx);
   cblas_dscal(len,-1.,x,incx);
 }
#endif

  /**@brief Set x[i]=-y[i] for len elements of the arrays x and y incremented by incx and incy, respectively.
     @param len gives the number of elements to be processed
     @param x points to the starting address of the destination
     @param incx gives the increment for x in each step
     @param y points to the starting address of the source
     @param incy gives the increment for y in each step
   */
template<class Val>
inline void mat_xemy( Integer len, Val *x, const Integer incx,
                          const Val *y, const Integer incy)
{
  const Val* const xend=x+len*incx;
  for(;x!=xend;x+=incx,y+=incy)
    (*x)=-(*y);
}

#ifdef WITH_BLAS
  /**@brief Set x[i]=a*y[i] for len elements of the arrays x and y.
     @param len gives the number of elements to be processed
     @param x points to the starting address of the destination
     @param y points to the starting address of the source
     @param a is the multiplicative factor
   */
inline void mat_xeya( Integer len, double *x,
                     const double *y, double a)
 {
   cblas_dcopy(len,y,1,x,1);
   cblas_dscal(len,a,x,1);
 }
#endif

  /**@brief Set x[i]=a*y[i] for len elements of the arrays x and y.
     @param len gives the number of elements to be processed
     @param x points to the starting address of the destination
     @param y points to the starting address of the source
     @param a is the multiplicative factor
   */
template<class Val>
inline void mat_xeya( Integer len, Val *x,
                     const Val *y, const Val a)
{
  if (a==Val(0)){
    mat_xea(len,x,Val(0));
    return;
  }
  if (a==Val(1)){
    mat_xey(len,x,y);
    return;
  }
  if (a==Val(-1)){
    mat_xemy(len,x,y);
    return;
  }
  const Val* const xend=x+len;
  for(;x!=xend;) 
    (*x++)=a*(*y++);
}

#ifdef WITH_BLAS
  /**@brief Set x[i]=a*y[i] for len elements of the arrays x and y incremented by incx and incy, respectively.
     @param len gives the number of elements to be processed
     @param x points to the starting address of the destination
     @param incx gives the increment for x in each step
     @param y points to the starting address of the source
     @param incy gives the increment for y in each step
     @param a is the multiplicative factor
   */
inline void mat_xeya( Integer len, double *x, const Integer incx,
		      const double *y, const Integer incy,double a)
 {
   cblas_dcopy(len,y,incy,x,incx);
   cblas_dscal(len,a,x,incx);
 }
#endif

  /**@brief Set x[i]=a*y[i] for len elements of the arrays x and y incremented by incx and incy, respectively.
     @param len gives the number of elements to be processed
     @param x points to the starting address of the destination
     @param incx gives the increment for x in each step
     @param y points to the starting address of the source
     @param incy gives the increment for y in each step
     @param a is the multiplicative factor
   */
template<class Val>
inline void mat_xeya( Integer len, Val *x, const Integer incx,
                     const Val *y, const Integer incy, const Val a)
{
  if (a==Val(0)){
    mat_xea(len,x,incx,Val(0));
    return;
  }
  if (a==Val(1)){
    mat_xey(len,x,incx,y,incy);
    return;
  }
  if (a==Val(-1)){
    mat_xemy(len,x,incx,y,incy);
    return;
  }
  const Val* const xend=x+len*incx;
  for(;x!=xend;x+=incx,y+=incy) 
     (*x)=a*(*y);
}

#ifdef WITH_BLAS
  /**@brief Set x[i]+=y[i] for len elements of the arrays x and y.
     @param len gives the number of elements to be processed
     @param x points to the starting address of the destination
     @param y points to the starting address of the source
   */
inline void mat_xpey( Integer len, double *x,
		      const double *y)
 {
   cblas_daxpy(len,1.,y,1,x,1);
 }
#endif

  /**@brief Set x[i]+=y[i] for len elements of the arrays x and y.
     @param len gives the number of elements to be processed
     @param x points to the starting address of the destination
     @param y points to the starting address of the source
   */
template<class Val>
inline void mat_xpey( Integer len, Val *x, const Val *y)
{
  const Val* const xend=x+len;
  for(;x!=xend;) 
   (*x++)+=(*y++);
}

#ifdef WITH_BLAS
  /**@brief Set x[i]+=y[i] for len elements of the arrays x and y incremented by incx and incy, respectively.
     @param len gives the number of elements to be processed
     @param x points to the starting address of the destination
     @param incx gives the increment for x in each step
     @param y points to the starting address of the source
     @param incy gives the increment for y in each step
   */
inline void mat_xpey( Integer len, double *x, const Integer incx,
		      const double *y, const Integer incy)
 {
   cblas_daxpy(len,1.,y,incy,x,incx);
 }
#endif

  /**@brief Set x[i]+=y[i] for len elements of the arrays x and y incremented by incx and incy, respectively.
     @param len gives the number of elements to be processed
     @param x points to the starting address of the destination
     @param incx gives the increment for x in each step
     @param y points to the starting address of the source
     @param incy gives the increment for y in each step
   */
template<class Val>
inline void mat_xpey( Integer len, Val *x, const Integer incx,
                         const Val *y, const Integer incy)
{
  const Val* const xend=x+len*incx;
  for(;x!=xend;x+=incx,y+=incy) 
    (*x)+=(*y);
}


  /**@brief Set x[i]*=y[i] for len elements of the arrays x and y.
     @param len gives the number of elements to be processed
     @param x points to the starting address of the destination
     @param y points to the starting address of the source
   */
template<class Val>
inline void mat_xhadey( Integer len, Val *x, const Val *y)
{
  const Val* const xend=x+len;
  for(;x!=xend;) 
    (*x++)*=(*y++);
}

  /**@brief Set x[i]/=y[i] for len elements of the arrays x and y, no zero checking!
     @param len gives the number of elements to be processed
     @param x points to the starting address of the destination
     @param y points to the starting address of the source
   */
template<class Val>
inline void mat_xinvhadey( Integer len, Val *x, const Val *y)
{
  const Val* const xend=x+len;
  for(;x!=xend;) 
    (*x++)/=(*y++);
}

  /**@brief Set x[i]*=y[i] for len elements of the arrays x and y incremented by incx and incy, respectively.
     @param len gives the number of elements to be processed
     @param x points to the starting address of the destination
     @param incx gives the increment for x in each step
     @param y points to the starting address of the source
     @param incy gives the increment for y in each step
   */
template<class Val>
inline void mat_xhadey( Integer len, Val *x, const Integer incx,
                          const Val *y, const Integer incy)
{
  const Val* const xend=x+len*incx;
  for(;x!=xend;x+=incx,y+=incy) 
    (*x)*=(*y);
}

  /**@brief Set x[i]/=y[i] for len elements of the arrays x and y incremented by incx and incy, respectively, no zero checking!
     @param len gives the number of elements to be processed
     @param x points to the starting address of the destination
     @param incx gives the increment for x in each step
     @param y points to the starting address of the source
     @param incy gives the increment for y in each step
   */
template<class Val>
inline void mat_xinvhadey( Integer len, Val *x, const Integer incx,
                          const Val *y, const Integer incy)
{
  const Val* const xend=x+len*incx;
  for(;x!=xend;x+=incx,y+=incy) 
    (*x)/=(*y);
}

#ifdef WITH_BLAS
  /**@brief Set x[i]+=a*y[i] for len elements of the arrays x and y.
     @param len gives the number of elements to be processed
     @param x points to the starting address of the destination
     @param y points to the starting address of the source
     @param a is the multiplicative factor
   */
inline void mat_xpeya( Integer len, Real *x,
		       const Real *y, const Real a)
{
   cblas_daxpy(len,a,y,1,x,1);
 }
#endif

  /**@brief Set x[i]+=a*y[i] for len elements of the arrays x and y.
     @param len gives the number of elements to be processed
     @param x points to the starting address of the destination
     @param y points to the starting address of the source
     @param a is the multiplicative factor
   */
template<class Val>
inline void mat_xpeya( Integer len, Val *x,
                         const Val *y, const Val a)
{
  if (a==Val(0))
    return;
  if (a==Val(1)){
    mat_xpey(len,x,y);
    return;
  }
  if (a==Val(-1)){
    mat_xmey(len,x,y);
    return;
  }
  const Val* const xend=x+len;
  for(;x!=xend;) 
    (*x++)+=a*(*y++);
}

#ifdef WITH_BLAS
  /**@brief Set x[i]+=a*y[i] for len elements of the arrays x and y incremented by incx and incy, respectively.
     @param len gives the number of elements to be processed
     @param x points to the starting address of the destination
     @param incx gives the increment for x in each step
     @param y points to the starting address of the source
     @param incy gives the increment for y in each step
     @param a is the multiplicative factor
   */
inline void mat_xpeya( Integer len, Real *x, const Integer incx,
                         const Real *y, const Integer incy, const Real a)
{
  cblas_daxpy(len,a,y,incy,x,incx);
 }
#endif

  /**@brief Set x[i]+=a*y[i] for len elements of the arrays x and y incremented by incx and incy, respectively.
     @param len gives the number of elements to be processed
     @param x points to the starting address of the destination
     @param incx gives the increment for x in each step
     @param y points to the starting address of the source
     @param incy gives the increment for y in each step
     @param a is the multiplicative factor
   */
template<class Val>
inline void mat_xpeya( Integer len, Val *x, const Integer incx,
                         const Val *y, const Integer incy, const Val a)
{
  if (a==Val(0))
    return;
  if (a==Val(1)){
    mat_xpey(len,x,incx,y,incy);
    return;
  }
  if (a==Val(-1)){
    mat_xmey(len,x,incx,y,incy);
    return;
  }
  const Val* const xend=x+len*incx;
  for(;x!=xend;x+=incx,y+=incy) 
    (*x)+=a*(*y);
}

#ifdef WITH_BLAS
  /**@brief Set x[i]=a*y[i]+b*x[i] for len elements of the arrays x and y.
     @param len gives the number of elements to be processed
     @param x points to the starting address of the destination
     @param y points to the starting address of the source
     @param a is the multiplicative factor for y
     @param b is the multiplicative factor for x
   */
inline void mat_xbpeya( Integer len, Real *x,
                         const Real *y, const Real a, const Real b)
{
  cblas_dscal(len,b,x,1);
  cblas_daxpy(len,a,y,1,x,1);
}
#endif

  /**@brief Set x[i]=a*y[i]+b*x[i] for len elements of the arrays x and y.
     @param len gives the number of elements to be processed
     @param x points to the starting address of the destination
     @param y points to the starting address of the source
     @param a is the multiplicative factor for y
     @param b is the multiplicative factor for x
   */
template<class Val>
inline void mat_xbpeya( Integer len, Val *x,
                         const Val *y, const Val a,const Val b)
{
  const Val* const xend=x+len;
  if (b!=Val(1)){
    for(;x!=xend;x++) 
      (*x)=b*(*x)+a*(*y++);
  }
  else {
    for(;x!=xend;) 
      (*x++)+=a*(*y++);
  }
}

#ifdef WITH_BLAS
  /**@brief Set x[i]=a*y[i]+b*x[i] for len elements of the arrays x and y incremented by incx and incy, respectively.
     @param len gives the number of elements to be processed
     @param x points to the starting address of the destination
     @param incx gives the increment for x in each step
     @param y points to the starting address of the source
     @param incy gives the increment for y in each step
     @param a is the multiplicative factor for y
     @param b is the multiplicative factor for x
   */
inline void mat_xbpeya( Integer len, Real *x, const Integer incx,
                         const Real *y, const Integer incy, const Real a, const Real b)
{
  cblas_dscal(len,b,x,incx);
  cblas_daxpy(len,a,y,incy,x,incx);
}
#endif

  /**@brief Set x[i]=a*y[i]+b*x[i] for len elements of the arrays x and y incremented by incx and incy, respectively.
     @param len gives the number of elements to be processed
     @param x points to the starting address of the destination
     @param incx gives the increment for x in each step
     @param y points to the starting address of the source
     @param incy gives the increment for y in each step
     @param a is the multiplicative factor for y
     @param b is the multiplicative factor for x
   */
template<class Val>
inline void mat_xbpeya( Integer len, Val *x, const Integer incx,
                         const Val *y, const Integer incy, const Val a,const Val b)
{
  const Val* const xend=x+len*incx;
  if (b!=Val(1)){
    for(;x!=xend;x+=incx,y+=incy) 
      (*x)=b*(*x)+a*(*y);
  }
  else {
    for(;x!=xend;x+=incx,y+=incy) 
      (*x)+=a*(*y);
  }
}

  /**@brief Set x[i]+=a for len elements of the array x.
     @param len gives the number of elements to be processed
     @param x points to the starting address 
     @param a is the scalar value
   */
template<class Val>
inline void mat_xpea( Integer len, Val *x, const Val a)
{
  if (a==Val(0))
    return;
  const Val* const xend=x+len;
  for(;x!=xend;) 
    (*x++)+=a;
}

  /**@brief Set x[i]+=a for len elements of the array x incremented by incx.
     @param len gives the number of elements to be processed
     @param x points to the starting address 
     @param incx gives the increment for x in each step
     @param a is the scalar value
   */
template<class Val>
inline void mat_xpea( Integer len, Val *x, const Integer incx,
                         const Val a)
{
  if (a==Val(0))
    return;
  const Val* const xend=x+len*incx;
  for(;x!=xend;x+=incx) 
    (*x)+=a;
}

#ifdef WITH_BLAS
  /**@brief Set x[i]*=a for len elements of the array x.
     @param len gives the number of elements to be processed
     @param x points to the starting address 
     @param a is the scalar value
   */
inline void mat_xmultea( Integer len, double *x,double a)
 {
   cblas_dscal(len,a,x,1);
 }
#endif

  /**@brief Set x[i]*=a for len elements of the array x.
     @param len gives the number of elements to be processed
     @param x points to the starting address 
     @param a is the scalar value
   */
template<class Val>
inline void mat_xmultea( Integer len, Val *x, const Val a)
{
  if (a==Val(1))
    return;
  const Val* const xend=x+len;
  for(;x!=xend;) 
    (*x++)*=a;
}

#ifdef WITH_BLAS
  /**@brief Set x[i]*=a for len elements of the array x incremented by incx.
     @param len gives the number of elements to be processed
     @param x points to the starting address 
     @param incx gives the increment for x in each step
     @param a is the scalar value
   */
inline void mat_xmultea( Integer len, double *x, const Integer incx,double a)
 {
   cblas_dscal(len,a,x,incx);
 }
#endif

  /**@brief Set x[i]*=a for len elements of the array x incremented by incx.
     @param len gives the number of elements to be processed
     @param x points to the starting address 
     @param incx gives the increment for x in each step
     @param a is the scalar value
   */
template<class Val>
inline void mat_xmultea( Integer len, Val *x, const Integer incx,
                         const Val a)
{
  if (a==Val(1))
    return;
  const Val* const xend=x+len*incx;
  for(;x!=xend;x+=incx) 
    (*x)*=a;
}

#ifdef WITH_BLAS
  /**@brief Set x[i]/=a for len elements of the array x.
     @param len gives the number of elements to be processed
     @param x points to the starting address 
     @param a is the scalar value (it is not checked whether a is zero)
   */
inline void mat_xdivea( Integer len, double *x,double a)
 {
   cblas_dscal(len,1./a,x,1);
 }
#endif

  /**@brief Set x[i]/=a for len elements of the array x.
     @param len gives the number of elements to be processed
     @param x points to the starting address 
     @param a is the scalar value (it is not checked whether a is zero)
   */
template<class Val>
inline void mat_xdivea( Integer len, Val *x, const Val a)
{
  if (a==Val(1))
    return;
  const Val* const xend=x+len;
  for(;x!=xend;) 
    (*x++)/=a;
}

#ifdef WITH_BLAS
  /**@brief Set x[i]/=a for len elements of the array x incremented by incx.
     @param len gives the number of elements to be processed
     @param x points to the starting address 
     @param incx gives the increment for x in each step
     @param a is the scalar value (it is not checked whether a is zero)
   */
inline void mat_xdivea( Integer len, double *x, const Integer incx,double a)
 {
   cblas_dscal(len,1./a,x,incx);
 }
#endif

  /**@brief Set x[i]/=a for len elements of the array x incremented by incx.
     @param len gives the number of elements to be processed
     @param x points to the starting address 
     @param incx gives the increment for x in each step
     @param a is the scalar value (it is not checked whether a is zero)
   */
template<class Val>
inline void mat_xdivea( Integer len, Val *x, const Integer incx,
                         const Val a)
{
  if (a==Val(1))
    return;
  const Val* const xend=x+len*incx;
  for(;x!=xend;x+=incx) 
    (*x)/=a;
}

  /**@brief Set x[i]%=a for len elements of the array x.
     @param len gives the number of elements to be processed
     @param x points to the starting address 
     @param a is the scalar value (it is not checked whether a is zero)
   */
template<class Val>
inline void mat_xmodea( Integer len, Val *x, const Val a)
{
  const Val* const xend=x+len;
  for(;x!=xend;) 
    (*x++)%=a;
}

  /**@brief Set x[i]%=a for len elements of the array x incremented by incx.
     @param len gives the number of elements to be processed
     @param x points to the starting address 
     @param incx gives the increment for x in each step
     @param a is the scalar value (it is not checked whether a is zero)
   */
template<class Val>
inline void mat_xmodea( Integer len, Val *x, const Integer incx,
                         const Val a)
{
  const Val* const xend=x+len*incx;
  for(;x!=xend;x+=incx) 
    (*x)%=a;
}

#ifdef WITH_BLAS
  /**@brief Set x[i]=y[i]+z[i] for len elements of the arrays x, y and z.
     @param len gives the number of elements to be processed
     @param x points to the starting address of the destination
     @param y points to the starting address of the first source
     @param z points to the starting address of the second source
   */
inline void mat_xeypz( Integer len, Real *x,
                         const Real *y,
		         const Real *z)
{
  cblas_dcopy(len,y,1,x,1);
  cblas_daxpy(len,1.,z,1,x,1);
}
#endif

  /**@brief Set x[i]=y[i]+z[i] for len elements of the arrays x, y and z.
     @param len gives the number of elements to be processed
     @param x points to the starting address of the destination
     @param y points to the starting address of the first source
     @param z points to the starting address of the second source
   */
template<class Val>
inline void mat_xeypz( Integer len, Val *x,
                      const Val *y, const Val *z)
{
  const Val* const xend=x+len;
  for(;x!=xend;) 
    (*x++)=(*y++)+(*z++);
}

#ifdef WITH_BLAS
  /**@brief Set x[i]=y[i]+z[i] for len elements of the arrays x, y and z incremented by incx, incy, and incz, respectively.
     @param len gives the number of elements to be processed
     @param x points to the starting address of the destination
     @param incx gives the increment for x in each step
     @param y points to the starting address of the first source
     @param incy gives the increment for y in each step
     @param z points to the starting address of the second source
     @param incz gives the increment for z in each step
   */
inline void mat_xeypz( Integer len, Real *x, const Integer incx,
                         const Real *y, const Integer incy,
                         const Real *z, const Integer incz)
{
  cblas_dcopy(len,y,incy,x,incx);
  cblas_daxpy(len,1.,z,incz,x,incx);
}
#endif

  /**@brief Set x[i]=y[i]+z[i] for len elements of the arrays x, y and z incremented by incx, incy, and incz, respectively.
     @param len gives the number of elements to be processed
     @param x points to the starting address of the destination
     @param incx gives the increment for x in each step
     @param y points to the starting address of the first source
     @param incy gives the increment for y in each step
     @param z points to the starting address of the second source
     @param incz gives the increment for z in each step
   */
template<class Val>
inline void mat_xeypz( Integer len, Val *x, const Integer incx,
                       const Val *y, const Integer incy,
                       const Val *z, const Integer incz)
{
  const Val* const xend=x+len*incx;
  for(;x!=xend;x+=incx,y+=incy,z+=incz) 
    (*x)=(*y)+(*z);
}

#ifdef WITH_BLAS
  /**@brief Set x[i]=y[i]-z[i] for len elements of the arrays x, y and z.
     @param len gives the number of elements to be processed
     @param x points to the starting address of the destination
     @param y points to the starting address of the first source
     @param z points to the starting address of the second source
   */
inline void mat_xeymz( Integer len, Real *x,
                         const Real *y,
                         const Real *z)
{
  cblas_dcopy(len,y,1,x,1);
  cblas_daxpy(len,-1.,z,1,x,1);
}
#endif

  /**@brief Set x[i]=y[i]-z[i] for len elements of the arrays x, y and z.
     @param len gives the number of elements to be processed
     @param x points to the starting address of the destination
     @param y points to the starting address of the first source
     @param z points to the starting address of the second source
   */
template<class Val>
inline void mat_xeymz( Integer len, Val *x,
                       const Val *y, const Val *z)
{
  const Val* const xend=x+len;
  for(;x!=xend;) 
    (*x++)=(*y++)-(*z++);
}

#ifdef WITH_BLAS
  /**@brief Set x[i]=y[i]-z[i] for len elements of the arrays x, y and z incremented by incx, incy, and incz, respectively.
     @param len gives the number of elements to be processed
     @param x points to the starting address of the destination
     @param incx gives the increment for x in each step
     @param y points to the starting address of the first source
     @param incy gives the increment for y in each step
     @param z points to the starting address of the second source
     @param incz gives the increment for z in each step
   */
inline void mat_xeymz( Integer len, Real *x, const Integer incx,
                         const Real *y, const Integer incy,
                         const Real *z, const Integer incz)
{
  cblas_dcopy(len,y,incy,x,incx);
  cblas_daxpy(len,-1.,z,incz,x,incx);
}
#endif

  /**@brief Set x[i]=y[i]-z[i] for len elements of the arrays x, y and z incremented by incx, incy, and incz, respectively.
     @param len gives the number of elements to be processed
     @param x points to the starting address of the destination
     @param incx gives the increment for x in each step
     @param y points to the starting address of the first source
     @param incy gives the increment for y in each step
     @param z points to the starting address of the second source
     @param incz gives the increment for z in each step
   */
template<class Val>
inline void mat_xeymz( Integer len, Val *x, const Integer incx,
                       const Val *y, const Integer incy,
                       const Val *z, const Integer incz)
{
  const Val* const xend=x+len*incx;
  for(;x!=xend;x+=incx,y+=incy,z+=incz) 
    (*x)=(*y)-(*z);
}

#ifdef WITH_BLAS
  /**@brief Set x[i]=a*y[i]+b*z[i] for len elements of the arrays x, y and z.
     @param len gives the number of elements to be processed
     @param x points to the starting address of the destination
     @param y points to the starting address of the first source
     @param z points to the starting address of the second source
     @param a is the multiplicative factor for y
     @param b is the multiplicative factor for z
   */
inline void mat_xeyapzb( Integer len, Real *x,
                         const Real *y,
                         const Real *z, const Real a, const Real b)
{
  cblas_dcopy(len,y,1,x,1);
  cblas_dscal(len,a,x,1);
  cblas_daxpy(len,b,z,1,x,1);
}
#endif

  /**@brief Set x[i]=a*y[i]+b*z[i] for len elements of the arrays x, y and z.
     @param len gives the number of elements to be processed
     @param x points to the starting address of the destination
     @param y points to the starting address of the first source
     @param z points to the starting address of the second source
     @param a is the multiplicative factor for y
     @param b is the multiplicative factor for z
   */
template<class Val>
inline void mat_xeyapzb( Integer len, Val *x,
                      const Val *y, const Val *z,const Val a,const Val b)
{
  const Val* const xend=x+len;
  for(;x!=xend;) 
    (*x++)=a*(*y++)+b*(*z++);
}

#ifdef WITH_BLAS
  /**@brief Set x[i]=a*y[i]+b*z[i] for len elements of the arrays x, y and z incremented by incx, incy, and incz, respectively.
     @param len gives the number of elements to be processed
     @param x points to the starting address of the destination
     @param incx gives the increment for x in each step
     @param y points to the starting address of the first source
     @param incy gives the increment for y in each step
     @param z points to the starting address of the second source
     @param incz gives the increment for z in each step
     @param a is the multiplicative factor for y
     @param b is the multiplicative factor for z
   */
inline void mat_xeyapzb( Integer len, Real *x, const Integer incx,
                         const Real *y, const Integer incy,
                         const Real *z, const Integer incz, const Real a, const Real b)
{
  cblas_dcopy(len,y,incy,x,incx);
  cblas_dscal(len,a,x,incx);
  cblas_daxpy(len,b,z,incz,x,incx);
}
#endif

  /**@brief Set x[i]=a*y[i]+b*z[i] for len elements of the arrays x, y and z incremented by incx, incy, and incz, respectively.
     @param len gives the number of elements to be processed
     @param x points to the starting address of the destination
     @param incx gives the increment for x in each step
     @param y points to the starting address of the first source
     @param incy gives the increment for y in each step
     @param z points to the starting address of the second source
     @param incz gives the increment for z in each step
     @param a is the multiplicative factor for y
     @param b is the multiplicative factor for z
   */
template<class Val>
inline void mat_xeyapzb( Integer len, Val *x, const Integer incx,
                       const Val *y, const Integer incy,
                       const Val *z, const Integer incz,const Val a,const Val b)
{
  const Val* const xend=x+len*incx;
  for(;x!=xend;x+=incx,y+=incy,z+=incz) 
    (*x)=a*(*y)+b*(*z);
}

#ifdef WITH_BLAS
  /**@brief return sum(x[i]*y[i]) summing over len elements of the arrays x and y.
     @param len gives the number of elements to be processed
     @param x points to the starting address of the x array
     @param y points to the starting address of the y array
     @return sum_{i=0}^{len} x[i] * y[i] 
   */
inline double mat_ip( Integer len, const double *x,const double *y)
{
  return cblas_ddot(len,x,1,y,1);
}
#endif

  /**@brief return sum(x[i]*y[i]) summing over len elements of the arrays x and y.
     @param len gives the number of elements to be processed
     @param x points to the starting address of the x array
     @param y points to the starting address of the y array
     @param d if !=0 points to the starting address of the d array for diangonal weights
     @return sum_{i=0}^{len} x[i] * y[i]  (* d[i])
   */
template<class Val>
inline Val mat_ip( Integer len, const Val *x, const Val *y,const Val* d=0)
{
  Val sum=0;
  const Val* const xend=x+len;
  if (d==0)
    for(;x!=xend;) 
      sum+=(*x++)*(*y++);
  else {
    for(;x!=xend;) 
      sum+=(*x++)*(*y++)*(*d++);
  }
 return sum;
}


#ifdef WITH_BLAS
  /**@brief return sum(x[i]*y[i]) summing over len elements of the arrays x and y incremented by incx and incy, respectively.
     @param len gives the number of elements to be processed
     @param x points to the starting address of the x array
     @param incx gives the increment for x in each step
     @param y points to the starting address of the y array
     @param incy gives the increment for y in each step
     @return sum_{i=0}^{len} x[i*incx] * y[i*incy] 
   */
inline double mat_ip( Integer len, const double *x, const Integer incx,
                   const double *y, const Integer incy)
{
  return cblas_ddot(len,x,incx,y,incy);
}
#endif

  /**@brief return sum(x[i]*y[i]) summing over len elements of the arrays x and y incremented by incx and incy, respectively.
     @param len gives the number of elements to be processed
     @param x points to the starting address of the x array
     @param incx gives the increment for x in each step
     @param y points to the starting address of the y array
     @param incy gives the increment for y in each step
     @param d if !=0 points to the starting address of the d array for diangonal weights
     @param incd gives the increment for d in each step
     @return sum_{i=0}^{len} x[i*incx] * y[i*incy] (* d(i))
   */
template<class Val>
inline Val mat_ip( Integer len, const Val *x, const Integer incx,
                   const Val *y, const Integer incy,const Val* d=0,const Integer incd=1)
{
  Val sum=0;
  const Val* const xend=x+len*incx;
  if (d==0)
    for(;x!=xend;x+=incx,y+=incy) 
      sum+=(*x)*(*y);
  else 
    for(;x!=xend;x+=incx,y+=incy,d+=incd) 
      sum+=(*x)*(*y)*(*d);
 return sum;
}

  /**@brief return sum(x[yind[j]]*yval[j]) summing over elements of the dense array x and a sparse array representation of y.
 
     Any elements outside the range of the other are considered 0.
  
     @param lenx gives the number of elements of x
     @param x points to the starting address of the x array
     @param leny gives the number of nonzero elements of y
     @param yval points to the starting address of the nonzero value array of y 
     @param yind points to the strating address of the nonzero index array of y
            (the indices are assumed to be in strictly increasing order)
     @param d if !=0 points to the starting address of the d array for diangonal weights
     @return \f$sum_{i=0}^{leny} x[i[indy[j]] * y[i] (* d[i[indy[j]])\f$ 
   */
template<class Val>
inline Val mat_ip_dense_sparse(const Integer lenx, const Val *x, Integer leny, const Val *yval, const Integer* yind,const Val* d=0)
{
  Val sum=0;
  const Integer* const yend=yind+leny;
  if (d==0)
    while((yind!=yend)&&(*yind<lenx))
      sum+=(*(x+*yind++))*(*yval++);
  else {
    while((yind!=yend)&&(*yind<lenx)){
      sum+=(*(x+*yind))*(*yval++)*(*(d+*yind));
      yind++;
    }
  }
  return sum;
}

  /**@brief return sum(xval[i]*yval[j] for i,j with xind[i]==yind[j]) summing over elements of the sparse array representations of x and y.
 
     Any elements outside the range of the other are considered 0.
  
     @param lenx gives the number of nonzero elements of x
     @param xval points to the starting address of the nonzero value arra of x 
     @param xind points to the strating address of the nonzero index array of x
            (the indices are assumed to be in strictly increasing order)
     @param leny gives the number of nonzero elements of y
     @param yval points to the starting address of the nonzero value array of y
     @param yind points to the starting address of the nonzero index array of y
            (the indices are assumed to be in strictly increasing order)
     @param d if !=0 points to the starting address of the d array for diangonal weights
     @return \f$ sum_{(i,j)\in J}xval[i] * yval[j] (* d[xind[i]])\f$ where \f$J=\{(i,j)\colon xind[i]=yind[j],i=1,\dots,lenx,j=1,\dots,leny\}\f$  
   */
template<class Val>
inline Val mat_ip_sparse_sparse(const Integer lenx, const Val *xval, const Integer* xind,const Integer leny, const Val *yval, const Integer* yind,const Val* d)
{
  if ((lenx<=0)||(leny<=0))
    return 0.;
  Val sum=0;
  const Integer* const xend=xind+lenx;
  const Integer* const yend=yind+leny;
  if (d==0) {
    do {
      if (*xind==*yind){
	sum+=(*xval++)*(*yval++);
	xind++;
	if (xind==xend)
	  break;
	yind++;      
      }
      while ((yind!=yend)&&(*yind<*xind)){
	yval++;
	yind++;
      }
      if (yind==yend) 
	break;
      while ((xind!=xend)&&(*xind<*yind)){
	xval++;
	xind++;
      }
    } while(xind!=xend);
  }
  else {
    do {
      if (*xind==*yind){
	sum+=(*xval++)*(*yval++)*(*(d+*xind));
	xind++;
	if (xind==xend)
	  break;
	yind++;      
      }
      while ((yind!=yend)&&(*yind<*xind)){
	yval++;
	yind++;
      }
      if (yind==yend) 
	break;
      while ((xind!=xend)&&(*xind<*yind)){
	xval++;
	xind++;
      }
    } while(xind!=xend);
  }
  return sum;
}
/*
{
  Val sum=0.;
  if((lenx>0)&&(leny>0)){
    const Integer* const xend=xind+lenx;
    const Integer* const yend=yind+leny;
    while(true){
      Integer diff=*xind-*yind;
      if (diff==0){
	sum+=(*xval++)*(*yval++);
	if(++xind==xend) 
	  break;
	if(++yind==yend) 
	  break;
      }
      else {
	if (diff>0){
	  yval++;
	  if (++yind==yend)
	    break;
	}
	else {
	  xval++;
	  if (++xind==xend)
	    break;
	}
      }
    }
  }
  return sum;
}
*/
  /*
{
  Val sum=0.;
  if((lenx>0)&&(leny>0)){
    while(true){
      Integer diff=*xind-*yind;
      if (diff==0){
	sum+=(*xval++)*(*yval++);
	xind++;yind++;
	if ((--lenx==0)||(--leny==0))
	  break;
      }
      else {
	if (diff>0){
	  yval++;
	  yind++;
	  if (--leny==0)
	    break;
	}
	else {
	  xval++;
	  xind++;
	  if(--lenx==0) 
	    break;
	}
      }
    }
  }
  return sum;
}
  */
 

#ifdef WITH_BLAS
  /**@brief return sum(x[i]*x[i]) summing over len elements of the array x.
     @param len gives the number of elements to be processed
     @param x points to the starting address of the x array
     @return sum_{i=0}^{len} x[i] * x[i] 
   */
inline double mat_ip( Integer len, const double *x)
{
  return cblas_ddot(len,x,1,x,1);
}
#endif

  /**@brief return sum(x[i]*x[i]) summing over len elements of the array x.
     @param len gives the number of elements to be processed
     @param x points to the starting address of the array
     @return sum_{i=0}^{len} x[i] * x[i] 
   */
template<class Val>
inline Val mat_ip( Integer len, const Val *x)
{
  Val sum=0;
  const Val* const xend=x+len;
  while(x!=xend) {
    const Val d=(*x++);
    sum+=d*d;
  }
  return sum;
}


#ifdef WITH_BLAS
  /**@brief return sum(x[i]*x[i]) summing over len elements of the array x incremented by incx.
     @param len gives the number of elements to be processed
     @param x points to the starting address of the x array
     @param incx gives the increment for x in each step
     @return sum_{i=0}^{len} x[i*incx] * x[i*incx] 
   */
inline double mat_ip( Integer len, const double *x, const Integer incx)
{
  return cblas_ddot(len,x,incx,x,incx);
}
#endif

  /**@brief return sum(x[i]*x[i]) summing over len elements of the array x incremented by incx.
     @param len gives the number of elements to be processed
     @param x points to the starting address of the destination
     @param incx gives the increment for x in each step
     @return sum_{i=0}^{len} x[i*incx] * x[i*incx] 
   */
template<class Val>
inline Val mat_ip( Integer len, const Val *x, const Integer incx)
{
  Val sum=0;
  const Val* const xend=x+len*incx;
  for(;x!=xend;x+=incx) {
    const Val d=(*x);
    sum+=d*d;
  }
 return sum;
}



  /**@brief returns sum(x[i]) over len elements of the array x.
     @param len gives the number of elements to be processed
     @param x points to the starting address of the first array
     @return sum_{i=0}^{len} x[i]
   */
template<class Val>
inline Val mat_sum( Integer len, const Val *x)
{
  Val sum=0;
  const Val* const xend=x+len;
  for(;x!=xend;) 
    sum+=(*x++);
 return sum;
}

  /**@brief returns sum(x[i]) over len elements of the array x incremented by incx.
     @param len gives the number of elements to be processed
     @param x points to the starting address of the first array
     @param incx gives the increment for x in each step
     @return sum_{i=0}^{len} x[i*incx]
   */
template<class Val>
inline Val mat_sum( Integer len, const Val *x,const Integer incx)
{
  Val sum=0;
  const Val* const xend=x+len*incx;
  for(;x!=xend;x+=incx) 
    sum+=(*x);
 return sum;
}

  /**@brief returns true if the elements of the arrays x and y are exactly equal.
     @param len gives the number of elements to be processed
     @param x points to the starting address of the first array
     @param y points to the starting address of the secibd array
     @return true if x[i]=y[i], 0<=i<len
   */
template<class Val>
inline bool mat_equal( Integer len, const Val *x,const Val *y)
{
  const Val* const xend=x+len;
  while(x!=xend){
    if (*x!=*y)
      break;
    x++,y++;
  }
  return (x==xend);
}


#ifdef WITH_BLAS
  /**@brief swap values x[i] and y[i] for len elements of the arrays x and y.
     @param len gives the number of elements to be processed
     @param x points to the starting address of the first array
     @param y points to the starting address of the second array
   */
inline void mat_swap( Integer len, double *x, double *y)
{
  cblas_dswap(len,x,1,y,1);
}
#endif

  /**@brief swap values x[i] and y[i] for len elements of the arrays x and y.
     @param len gives the number of elements to be processed
     @param x points to the starting address of the first array
     @param y points to the starting address of the second array
   */
template<class Val>
inline void mat_swap( Integer len, Val *x, Val *y)
{
  const Val* const xend=x+len;
  for(;x!=xend;) 
    {Val h=*x;(*x++)=*y;(*y++)=h;}
}

#ifdef WITH_BLAS
  /**@brief swap values x[i] and y[i] for len elements of the arrays x and y incremented by incx and incy, respectively.
     @param len gives the number of elements to be processed
     @param x points to the starting address of the first array
     @param incx gives the increment for x in each step
     @param y points to the starting address of the second array
     @param incy gives the increment for y in each step
   */
inline void mat_swap( Integer len, double *x, const Integer incx, double *y, const Integer incy)
{
  cblas_dswap(len,x,incx,y,incy);
}
#endif

  /**@brief swap values x[i] and y[i] for len elements of the arrays x and y incremented by incx and incy, respectively.
     @param len gives the number of elements to be processed
     @param x points to the starting address of the first array
     @param incx gives the increment for x in each step
     @param y points to the starting address of the second array
     @param incy gives the increment for y in each step
   */
template<class Val>
inline void mat_swap( Integer len, Val *x, const Integer incx,
                      Val *y, const Integer incy)
{
  const Val* const xend=x+len*incx;
  for(;x!=xend;x+=incx,y+=incy) 
    {Val h=*x;(*x)=*y;(*y)=h;}
}

//@}


  /**@brief "less"-routine for sorting indices of value arrays by std::sort
   */
template<class Val>
struct mat_less_index {
  /// constructor initializing _values
  mat_less_index( const Val* values ):
    _values(values)
  {}
  
  /// returns (values[i] < _values[j])
  bool operator()( Integer i, Integer j) const {
    return _values[i] < _values[j];
  }
  
private:
  /// const pointer to values
  const Val* _values;
  
};


  /**@brief "greater"-routine for sorting indices of value arrays by std::sort
   */
template<class Val>
struct mat_greater_index {
  /// constructor initializing _values
  mat_greater_index( const Val* values ):
    _values(values)
  {}

  /// returns (values[i] > _values[j])
  bool operator()( Integer i, Integer j) const {
    return _values[i] > _values[j];
  }
  
private:
  /// const pointer to values
  const Val* _values;
  
};



//-----------------------------------------------------------------------------

/**@defgroup matop_globals Global Variables and Objects
   @brief The variables and objects are instantiated in memarray.cxx
*/
//@{ 


// **************************************************************************
//                                mat_randgen
// **************************************************************************

     
/** @brief common random number generator for use when a random matrix is generated.

    At any time it is possible to reset the random generator by a seed, see #CH_Tools::GB_rand 
 */
  extern CH_Tools::GB_rand mat_randgen; 
  extern std::default_random_engine mat_std_randgen; ///< optional fast alternative random number generator from std
extern std::mt19937 mat_mt_randgen; ///< optional high quality alternative random number generator from std
extern std::mt19937_64 mat_mt64_randgen; ///< optional high quality alternative random number generator from std
  

// **************************************************************************
//                                materrout
// **************************************************************************


/** @brief if not zero, this is the output stream for runtime error messages, by default it is set to &std::cout 
 */
extern std::ostream* materrout; 

//@}





//-----------------------------------------------------------------------------

/**@defgroup matop_matrixerror Classes and Functions used for Debugging 
   @brief To catch MatrixErrors in CONICBUNDLE_DEBUG mode stop/break in #CH_Matrix_Classes::MEmessage() (first routine in memarray.cxx)

   Implementation started before exceptions existed, so these routines where 
   developed for error checking during development. They are only fully 
   functioning if the code is compiled with CONICBUNDLE_DEBUG being defined. If DEBUG
   is undefined, no range checks or checks on consistency of dimensions
   are performed. 
*/

//@{

// **************************************************************************
//                               Matrix Types
// **************************************************************************

/// serves for specifying the source (matrix class or function) of the error
enum Mtype {
  MTglobalfun,   ///< error arises in a global function
  MTindexmatrix, ///< error arises in a message of #CH_Matrix_Classes::Indexmatrix
  MTmatrix,      ///< error arises in a message of #CH_Matrix_Classes::Matrix
  MTsymmetric, ///< error arises in a message of #CH_Matrix_Classes::Symmatrix
  MTsparse,    ///< error arises in a message of #CH_Matrix_Classes::Sparsemat
  MTsparsesym  ///< error arises in a message of #CH_Matrix_Classes::Sparsesym
};

// **************************************************************************
//                               "exceptions"
// **************************************************************************

/// serves for specifying the error type.
enum MEcode {
  ME_unspec, ///< unspecified error type, not in the list below
  ME_range,  ///< error due to range check
  ME_mem,    ///< error due to insufficient memory
  ME_dim,    ///< error due to inconsistent dimesions
  ME_num,    ///< error due to numerical reasons (e.g. division by zero)
  ME_warning ///< no error, but a possible source of difficulties
};

/// Such an object is generated and passed to #MEmessage(), whenever an error occurs. It holds some output information on the error.
class MatrixError
{
public:
  MEcode code;            ///< see #MEcode for allowed error types
    const char *message;  ///< the error mussage pointer must point to existing object and does not get freed
    Mtype mtype;          ///< see #Mtype for allowed error source types

    ///there is only this one constructor and no other messages
    MatrixError(MEcode c=ME_unspec,const char *mes=0,Mtype mt=MTmatrix):
        code(c),message(mes),mtype(mt){}
};

/// Such an object is generated and passed to #MEmessage() whenever some index is out of range.
class MErange:public MatrixError
{
public:
  Integer nr; ///< number of rows
  Integer r;  ///< row index
  Integer nc; ///< number of columsn
  Integer c;  ///< column index

  /// constructor 
    MErange(Integer inr,Integer ir,Integer inc,Integer ic,const char *mes=0,Mtype mt=MTmatrix):
        MatrixError(ME_range,mes,mt),nr(inr),r(ir),nc(inc),c(ic){}
};

/// Such an object is generated and passed to #MEmessage() whenever a memory allocation fails.
class MEmem:public MatrixError
{
public:
  Integer size; ///< memory size requested

  /// constructor 
    MEmem(Integer s,const char *mes=0,Mtype mt=MTmatrix):
        MatrixError(ME_mem,mes,mt),size(s){}
};

/// Such an object is generated and passed to #MEmessage() whenever matrix dimensions do not agree for a desired operation.
class MEdim:public MatrixError
{
public:
  Integer nr1;  ///< number of rows of object 1
  Integer nc1;  ///< number of columns of object 1
  Integer nr2;  ///< number of rows of object 2
  Integer nc2;  ///< number of columns of object 2

  /// constructor
    MEdim(Integer r1,Integer c1,Integer r2,Integer c2,const char *mes=0,Mtype mt=MTmatrix):
        MatrixError(ME_dim,mes,mt),nr1(r1),nc1(c1),nr2(r2),nc2(c2){}
};

///displays an error message and terminates via abort() or returns 1 in case of warnings.
int MEmessage(const MatrixError&); 


// **************************************************************************
//                         CONICBUNDLE_DEBUG error checking templates
// **************************************************************************

#if (CONICBUNDLE_DEBUG>=1)

/// checks whether A is initiliazed and raises a MatrixError if not
template<class Mat>
inline void chk_init(const Mat& A)
{
 if (!A.get_init())
     MEmessage(MatrixError(ME_unspec,"Matrix not initialized",A.get_mtype()));
}

/// checks whether matrices A and B can be added and raises the MatrixError MEdim if not
template<class Mat1,class Mat2>
inline void chk_add(const Mat1& A,const Mat2& B)
{
 chk_init(A); chk_init(B);
 if ((A.dim()==Integer(0))&&(B.dim()==Integer(0))) return;
 if ((A.rowdim()!=B.rowdim())||(A.coldim()!=B.coldim()))
   MEmessage(MEdim(A.rowdim(),A.coldim(),B.rowdim(),B.coldim(),
		   "dimensions do not match in additive expression",
		   A.get_mtype()));
}

/// checks whether matrices A and B can be multiplied and raises the MatrixError MEdim if not
template<class Mat1,class Mat2>
inline void chk_mult(const Mat1& A,const Mat2& B,int atrans=0,int btrans=0)
{
 chk_init(A); chk_init(B);
 if ((A.dim()==Integer(0))&&(B.dim()==Integer(0))) return;
 if (
     ((atrans==0)&&(btrans==0)&&(A.coldim()!=B.rowdim())) 
     ||
     ((atrans==1)&&(btrans==0)&&(A.rowdim()!=B.rowdim())) 
     ||
     ((atrans==0)&&(btrans==1)&&(A.coldim()!=B.coldim())) 
     ||
     ((atrans==1)&&(btrans==1)&&(A.rowdim()!=B.coldim())) 
     )
     MEmessage(MEdim(A.rowdim(),A.coldim(),B.rowdim(),B.coldim(),
                     "dimensions do not match in multiplicative expression",
                     A.get_mtype()));
}

/// checks whether matrices A and B have the same row dimension and raises the MatrixError MEdim if not
template<class Mat1,class Mat2>
inline void chk_rowseq(const Mat1& A,const Mat2& B)
{
 chk_init(A); chk_init(B);
 if ((A.dim()==Integer(0))&&(B.dim()==Integer(0))) return;

 if (A.rowdim()!=B.rowdim())
     MEmessage(MEdim(A.rowdim(),A.coldim(),B.rowdim(),B.coldim(),
                     "number of rows do not agree",
                     A.get_mtype()));
}

/// checks whether matrices A and B have the same column dimension and raises the MatrixError MEdim if not
template<class Mat1,class Mat2>
inline void chk_colseq(const Mat1& A,const Mat2& B)
{
 chk_init(A); chk_init(B);
 if ((A.dim()==Integer(0))&&(B.dim()==Integer(0))) return;
 if (A.coldim()!=B.coldim())
     MEmessage(MEdim(A.rowdim(),A.coldim(),B.rowdim(),B.coldim(),
                     "number of columns do not agree",
                     A.get_mtype()));
}

/// checks whether indices r and c fit into [0,ubr-1] and [0,ubc-1], respectively, and raises the MatrixError MErange if not
inline void chk_range(Integer r,Integer c,Integer ubr,Integer ubc)
{
 if ((r<0)||(c<0)||
     ((ubr>=0)&&(r>=ubr))||
     ((ubc>=0)&&(c>=ubc)))
     MEmessage(MErange(r,ubr,c,ubc,"index out of range or negative dimension"));
}

/// checks whether index i fits into [0,ub-1] and raises the MatrixError MErange if not
inline void chk_single_range(Integer i,Integer ub)
{
 if ((i<0)||
     ((ub>=0)&&(i>=ub)))
     MEmessage(MErange(i,ub,0,0,"index out of range or negative dimension"));
}

/// allows to set the initialized flag of matrix A to val, needed e.g. for Matrix::newsize. The code is not generated if CONICBUNDLE_DEBUG is not defined. 
template<class Mat>
inline void chk_set_init(Mat &A,bool val){A.set_init(val);}

#else

// if CONICBUNDLE_DEBUG is not definied all chk_functions are mapped to the null string

/// CONICBUNDLE_DEBUG being undefined, the template function is removed. Otherwise it would check, whether x is initialized 
#define chk_init(x)
/// CONICBUNDLE_DEBUG being undefined, the template function is removed. Otherwise it would check, whether matrices x and y can be added 
#define chk_add(x,y)
/// CONICBUNDLE_DEBUG being undefined, the template function is removed. Otherwise it would check, whether matrices x and y can be multiplied
template<class Mat1,class Mat2>
inline void chk_mult(const Mat1&,const Mat2&,int =0,int =0){}
/// CONICBUNDLE_DEBUG being undefined, the template function is removed. Otherwise it would check, whether matrices x and y have the same row dimension 
#define chk_rowseq(x,y)
/// CONICBUNDLE_DEBUG being undefined, the template function is removed. Otherwise it would check, whether matrices x and y have the same column dimension
#define chk_colseq(x,y)
/// CONICBUNDLE_DEBUG being undefined, the template function is removed. Otherwise it would check, whether i and j fit into [0,ubi-1] amd [0,ubj-1] respectively.  
#define chk_range(i,j,ubi,ubj)
/// CONICBUNDLE_DEBUG being undefined, the template function is removed. Otherwise it would allow to set the initialized flag of matrix x to value y. 
#define chk_set_init(x,y)
/// CONICBUNDLE_DEBUG being undefined, the template function is removed. Otherwise it would whether index i fits into [0,ub-1] and raise the MatrixError MErange if not. 
#define chk_single_range(i,ub)

#endif 

//@}

}

#endif
