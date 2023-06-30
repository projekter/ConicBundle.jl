/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/InteriorPointBlock.hxx
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



#ifndef CONICBUNDLE_INTERIORPOINTBLOCK_HXX
#define CONICBUNDLE_INTERIORPOINTBLOCK_HXX

/**  @file InteriorPointBlock.hxx
    @brief Header declaring the class ConicBundle::InteriorPointBlock
    @version 1.0
    @date 2020-05-02
    @author Christoph Helmberg
*/

#include "symmat.hxx"
#include "CBout.hxx"

namespace ConicBundle {

  /** @ingroup ConstrainedQPSolver
   */

   //@{

   /** @brief  abstract interface for interior point vector/matrix variables and routines specific to primal dual complementarity conditions of symmetric cones

       This is the general part independent of the bundle framework.

   */


  class InteriorPointBlock : public virtual CBout {
  protected:
    /// for a polynomial q0+alpha*(q1+alpha*(q2+alpha*(q3+alpha*q4))) with q0<=0 find the maximum alpha <= stepsize, that keeps its value <=0. and put stepsize=min(stepsize,alpha)
    int pol_le_zero_step(CH_Matrix_Classes::Real& stepsize,
      CH_Matrix_Classes::Real q0,
      CH_Matrix_Classes::Real q1,
      CH_Matrix_Classes::Real q2,
      CH_Matrix_Classes::Real q3,
      CH_Matrix_Classes::Real q4,
      CH_Matrix_Classes::Real abseps = 1e-10) const;

    /// find the minimizing alpha of  a polynomial q0+alpha*(q1+alpha*(q2+alpha*(q3+alpha*q4))) within 0<=alpha<=stepsize and put stepsize=alpha, where q1 is assumed to be <0  or be ==0  becoming negative for small alpa>=0.
    int minimize_pol_step(CH_Matrix_Classes::Real& stepsize,
      CH_Matrix_Classes::Real q0,
      CH_Matrix_Classes::Real q1,
      CH_Matrix_Classes::Real q2,
      CH_Matrix_Classes::Real q3,
      CH_Matrix_Classes::Real q4,
      CH_Matrix_Classes::Real abseps = 1e-10) const;

    /// find a stepsize so that outside the nbh_ubnd some progress is made towards it and inside the upper bound the step size is chosen as large as possible while staying within
    int control_nbh_step(CH_Matrix_Classes::Real& stepsize,
      CH_Matrix_Classes::Real& max_nbh,
      CH_Matrix_Classes::Real nbh_ubnd,
      CH_Matrix_Classes::Real mu_xz,
      CH_Matrix_Classes::Real mu_xdzpdxz,
      CH_Matrix_Classes::Real mu_dxdz,
      CH_Matrix_Classes::Real nrmsqr_xz,
      CH_Matrix_Classes::Real nrmsqr_xdzpdxz,
      CH_Matrix_Classes::Real nrmsqr_dxdz,
      CH_Matrix_Classes::Real ip_xz_xdzpdxz,
      CH_Matrix_Classes::Real ip_xz_dxdz,
      CH_Matrix_Classes::Real ip_dxdz_xdzpdxz) const;

  public:
    ///default constructor
    InteriorPointBlock(CBout* cb = 0, int cbinc = -1) :CBout(cb, cbinc) {
    }

    ///virtual destructor (implemented in InteriorPointBundleBlock.cxx)
    virtual ~InteriorPointBlock();

    /// the dimension of the variable
    virtual CH_Matrix_Classes::Integer get_vecdim() const = 0;

    /// set x to value*"one" to x (spend a total of mu_dim*val), or if add==true, add value*"one" to x
    virtual int center_x(CH_Matrix_Classes::Real val, bool add = false) = 0;

    /// set z to value*"one" to z, or if add==true, add value*"one" to z
    virtual int center_z(CH_Matrix_Classes::Real val, bool add = false) = 0;

    /// set x to the values of vec[startindex+0,+1 ...,+(vecdim-1)] and return in add_center_value a value>=0 that needs to be added via center_x to make it feasible
    virtual int set_x(const CH_Matrix_Classes::Matrix& vec,
      CH_Matrix_Classes::Integer startindex,
      CH_Matrix_Classes::Real& add_center_value) = 0;

    /// set z to the values of vec[startindex+0,+1 ...,+(vecdim-1)] and return a value>=0 that needs to be added via center_z to make it feasible
    virtual int set_z(const CH_Matrix_Classes::Matrix& vec,
      CH_Matrix_Classes::Integer startindex,
      CH_Matrix_Classes::Real& add_center_value) = 0;

    /// on vec[startindex+0,+1 ...,+(vecdim-1)] put or add  a * x into vec for a real number a  
    virtual int vecgetsax(CH_Matrix_Classes::Matrix& vec,
      CH_Matrix_Classes::Integer startindex,
      CH_Matrix_Classes::Real a = 1.,
      bool add = false) = 0;

    /// on vec[startindex+0,+1 ...,+(vecdim-1)] put or add a * z into vec for a real number a   
    virtual int vecgetsaz(CH_Matrix_Classes::Matrix& vec,
      CH_Matrix_Classes::Integer startindex,
      CH_Matrix_Classes::Real a = 1.,
      bool add = false) = 0;

    /// add dimensions of the primal-dual pairs to mudim and add the "trace" (the inner product with center) of the respective primal-dual pair products for the current step; update the min and max values of x_i*z_i
    virtual int get_mu_info(CH_Matrix_Classes::Integer& mudim,
      CH_Matrix_Classes::Real& tr_xz,
      CH_Matrix_Classes::Real& tr_xdzpdxz,
      CH_Matrix_Classes::Real& tr_dxdz,
      CH_Matrix_Classes::Real& min_xz,
      CH_Matrix_Classes::Real& max_xz) const = 0;

    /// for limiting the stepsize with respect to the neighborhood this information about norms and inner products of x(.)*z-tr_xz-tr_xz/mudim(.*)1, x.()*dz+dx(.)*z-tr_xdzpdxz/mudim(.*)1, and dx(.)*dz-tr_dxdz/mudim(.)*1 is required, each block *adds* its contribution to the numbers
    virtual int get_nbh_info(CH_Matrix_Classes::Integer mudim,
      CH_Matrix_Classes::Real tr_xz,
      CH_Matrix_Classes::Real tr_xdzpdxz,
      CH_Matrix_Classes::Real tr_dxdz,
      CH_Matrix_Classes::Real nbh_ubnd,
      CH_Matrix_Classes::Real& alpha,
      CH_Matrix_Classes::Real& max_nbh,
      CH_Matrix_Classes::Real& nrmsqr_xz,
      CH_Matrix_Classes::Real& nrmsqr_xdzpdxz,
      CH_Matrix_Classes::Real& nrmsqr_dxdz,
      CH_Matrix_Classes::Real& ip_xz_xdzpdxz,
      CH_Matrix_Classes::Real& ip_xz_dxdz,
      CH_Matrix_Classes::Real& ip_dxdz_xdzpdxz) const = 0;


    /// if necessary, reduce alpha to the biggest value so that feasibility is maintained with this step size
    virtual int linesearch(CH_Matrix_Classes::Real& alpha) const = 0;

    /// compute the complementarity_rhs=rhsmu*xi-rhscorr*xi*dx*dz (wihtout "-z") for mu=rhsmu and for corrector for factor rhscorr>0., store this and add it to rhs;
    virtual int add_muxinv(CH_Matrix_Classes::Matrix& rhs,
      CH_Matrix_Classes::Integer startindex,
      CH_Matrix_Classes::Real rhsmu,
      CH_Matrix_Classes::Real rhscorr,
      bool minus = false) = 0;

    /// extract dx from rhs at startindex and compute at the same time dz (=-sys dx-z  +complentarity_rhs); this may only be called after add_muxinv() was called for this point
    virtual int set_dx(const CH_Matrix_Classes::Matrix& rhs, CH_Matrix_Classes::Integer startindex) = 0;

    /// compute dx=sysinv*rhs and at the same time dz (=-rhs-z +complentarity_rhs); may only be called after add_muxinv was called for the most recent point; this may only be called after add_muxinv() was called for this point
    virtual int set_dx_xizsolverhs(const CH_Matrix_Classes::Matrix& rhs,
      CH_Matrix_Classes::Integer startindex) = 0;

    /// compute sysinv*rhs into rhs, possibly with a negative sign 
    virtual int apply_xizinv(CH_Matrix_Classes::Matrix& rhs,
      CH_Matrix_Classes::Integer startindex,
      bool minus = false) = 0;

    /// compute sys*rhs into rhs, possibly with a negative sign
    virtual int apply_xiz(CH_Matrix_Classes::Matrix& rhs,
      CH_Matrix_Classes::Integer startindex,
      bool minus = false) = 0;

    /// move to (x+alpha*dx, z+alpha*dz)
    virtual int do_step(CH_Matrix_Classes::Real alpha) = 0;


    /// add the Schur complement to a big system matrix
    virtual int add_AxizinvAt(const CH_Matrix_Classes::Matrix& A,
      CH_Matrix_Classes::Symmatrix& globalsys,
      bool minus = false,
      bool Atrans = false) = 0;

    /// add (or subract if minus==true) the system matrix to a big system matrix starting at startindex
    virtual int add_xiz(CH_Matrix_Classes::Symmatrix& globalsys,
      CH_Matrix_Classes::Integer startindex,
      bool minus = false) = 0;

    //---------------------- mainly for testing

    /// return the vector form of x 
    virtual int get_vecx(CH_Matrix_Classes::Matrix& vecx, CH_Matrix_Classes::Integer startindex) = 0;

    /// return the vector form of z
    virtual int get_vecz(CH_Matrix_Classes::Matrix& vecz, CH_Matrix_Classes::Integer startindex) = 0;

    /// return the vector form of dx, 1 if not available 
    virtual int get_vecdx(CH_Matrix_Classes::Matrix& vecdx, CH_Matrix_Classes::Integer startindex) = 0;

    /// return the vector form of dz, 1 if not available
    virtual int get_vecdz(CH_Matrix_Classes::Matrix& vecdz, CH_Matrix_Classes::Integer startindex) = 0;


  };

  //--------------------------------------------------------
  /// computes values for a neighborhood line search for a primal nonnegative cone pair
  inline void NNC_nbh_stepsize(CH_Matrix_Classes::Real x,
    CH_Matrix_Classes::Real z,
    CH_Matrix_Classes::Real dx,
    CH_Matrix_Classes::Real dz,
    CH_Matrix_Classes::Real mu_xz,
    CH_Matrix_Classes::Real mu_xdzpdxz,
    CH_Matrix_Classes::Real mu_dxdz,
    CH_Matrix_Classes::Real mu_at_one,
    CH_Matrix_Classes::Real nbh_ubnd,
    CH_Matrix_Classes::Real& alpha,
    CH_Matrix_Classes::Real& max_nbh,
    CH_Matrix_Classes::Real& nrmsqr_xz,
    CH_Matrix_Classes::Real& nrmsqr_xdzpdxz,
    CH_Matrix_Classes::Real& nrmsqr_dxdz,
    CH_Matrix_Classes::Real& ip_xz_xdzpdxz,
    CH_Matrix_Classes::Real& ip_xz_dxdz,
    CH_Matrix_Classes::Real& ip_dxdz_xdzpdxz) {
    CH_Matrix_Classes::Real xz = x * z - mu_xz;
    CH_Matrix_Classes::Real xdzpdxz = x * dz + dx * z - mu_xdzpdxz;
    CH_Matrix_Classes::Real dxdz = dx * dz - mu_dxdz;
    nrmsqr_xz += CH_Matrix_Classes::sqr(xz);
    nrmsqr_xdzpdxz += CH_Matrix_Classes::sqr(xdzpdxz);
    nrmsqr_dxdz += CH_Matrix_Classes::sqr(dxdz);
    ip_xz_xdzpdxz += xz * xdzpdxz;
    ip_xz_dxdz += xz * dxdz;
    ip_dxdz_xdzpdxz += dxdz * xdzpdxz;

    if (x + alpha * dx < CH_Matrix_Classes::eps_Real * mu_xz)
      alpha = CH_Matrix_Classes::min(alpha, CH_Matrix_Classes::max(0., (CH_Matrix_Classes::eps_Real * mu_xz - x) / dx));
    if (z + alpha * dz < CH_Matrix_Classes::eps_Real * mu_xz)
      alpha = CH_Matrix_Classes::min(alpha, CH_Matrix_Classes::max(0., (CH_Matrix_Classes::eps_Real * mu_xz - z) / dz));
    if (alpha > CH_Matrix_Classes::eps_Real) {
      CH_Matrix_Classes::Real otheta = std::fabs(xz) / mu_xz;
      CH_Matrix_Classes::Real gtheta = nbh_ubnd;
      CH_Matrix_Classes::Real theta = CH_Matrix_Classes::max(otheta, nbh_ubnd);
      if (otheta > (1 + 1e-6) * nbh_ubnd) {
        gtheta = std::sqrt(CH_Matrix_Classes::max(0., otheta * otheta + 2 * (xz * xdzpdxz - otheta * otheta * mu_xz * mu_xdzpdxz) / CH_Matrix_Classes::max(mu_at_one, 1e-6 * mu_xz)));
        if (gtheta < .1 * nbh_ubnd + 0.9 * otheta) {
          gtheta = CH_Matrix_Classes::max(.1 * gtheta + .9 * otheta, nbh_ubnd);
        } else {
          gtheta = (1 - 1e-3) * otheta + 1e-3 * nbh_ubnd;
        }
      }

      //find the longest step for the positive branch of the absolute value
      CH_Matrix_Classes::Real q0 = CH_Matrix_Classes::min(xz - theta * mu_xz, 0.);
      CH_Matrix_Classes::Real q1 = xdzpdxz - theta * mu_xdzpdxz - (gtheta - theta) * (mu_at_one);
      CH_Matrix_Classes::Real q2 = dxdz - theta * mu_dxdz;

      if (q0 + 1e-8 * q1 > CH_Matrix_Classes::eps_Real * mu_xz) {
        alpha = 0.;
      } else if (q2 > CH_Matrix_Classes::eps_Real * mu_xz) {
        //strictly convex case with min >= zero, max is upper root
        CH_Matrix_Classes::Real p = q1 / q2 / 2.;
        CH_Matrix_Classes::Real q = q0 / q2;
        alpha = CH_Matrix_Classes::min(alpha, -p + std::sqrt(p * p - q));
      } else if (q2 < -CH_Matrix_Classes::eps_Real * mu_xz) {
        //strictly concave case, restrictions only for positive q1
        if (q1 > 0.) {
          CH_Matrix_Classes::Real p = q1 / q2 / 2.;
          CH_Matrix_Classes::Real q = q0 / q2;
          CH_Matrix_Classes::Real discr = p * p - q;
          if (discr > 0) {
            alpha = CH_Matrix_Classes::min(alpha, -p - std::sqrt(discr));
          }
        }
      } else {
        //q is affine
        if (q0 + alpha * q1 > 0.) {
          alpha = CH_Matrix_Classes::min(alpha, CH_Matrix_Classes::max(0., -q0 / q1));
        }
      }

      //find the longest step for the negative branch of the absolute value
      q0 = -xz - theta * mu_xz;
      q1 = -xdzpdxz - theta * mu_xdzpdxz - (gtheta - theta) * (mu_at_one);
      q2 = -dxdz - theta * mu_dxdz;
      if (q0 + 1e-8 * q1 > CH_Matrix_Classes::eps_Real * mu_xz) {
        alpha = 0;
      } else if (q2 > CH_Matrix_Classes::eps_Real * mu_xz) {
        //strictly convex case with min >= zero, max is upper root
        CH_Matrix_Classes::Real p = q1 / q2 / 2.;
        CH_Matrix_Classes::Real q = q0 / q2;
        alpha = CH_Matrix_Classes::min(alpha, -p + std::sqrt(p * p - q));
      } else if (q2 < -CH_Matrix_Classes::eps_Real * mu_xz) {
        //strictly concave case, restrictions only for positive q1
        if (q1 > 0.) {
          CH_Matrix_Classes::Real p = q1 / q2 / 2.;
          CH_Matrix_Classes::Real q = q0 / q2;
          CH_Matrix_Classes::Real discr = p * p - q;
          if (discr > 0) {
            alpha = CH_Matrix_Classes::min(alpha, -p - std::sqrt(discr));
          }
        }
      } else {
        //q is affine
        if (q0 + alpha * q1 > 0.) {
          alpha = CH_Matrix_Classes::min(alpha, CH_Matrix_Classes::max(0., -q0 / q1));
        }
      }
    }

    CH_Matrix_Classes::Real step_nbh = std::fabs(xz + alpha * (xdzpdxz + alpha * dxdz)) / CH_Matrix_Classes::max(mu_xz + alpha * (mu_xdzpdxz + alpha * (mu_dxdz)), 1e-6 * mu_xz);
    max_nbh = CH_Matrix_Classes::max(max_nbh, step_nbh);
  }


  //@}

}

#endif

