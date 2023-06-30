/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/InteriorPointBlock.cxx
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


#include "NNCIPBlock.hxx"

using namespace CH_Matrix_Classes;

namespace ConicBundle {

  InteriorPointBlock::~InteriorPointBlock() {
  }

  int InteriorPointBlock::pol_le_zero_step(Real& stepsize,
    Real q0,
    Real q1,
    Real q2,
    Real q3,
    Real q4,
    Real abseps) const {
    if (q0 + 1e-8 * q1 > 0.) {
      stepsize = 0;
      return 0;
    }
    Real maxnrmq = max(std::fabs(q0), max(std::fabs(q1), max(std::fabs(q2), max(std::fabs(q3), std::fabs(q4)))));

    // find the maximum stepval in [0,1] so that q(stepval)<=0.
    // note: q(0)=q0<=0.
    Real ubstep = stepsize;
    Real stepval = stepsize;
    bool rootsearch = true;

    //ploynomial coefficients of the derivative dq
    Real dq0 = q1;
    Real dq1 = 2 * q2;
    Real dq2 = 3 * q3;
    Real dq3 = 4 * q4;

    //ploynomial coefficients of the second derivative d2q
    Real d2q0 = dq1;
    Real d2q1 = 2 * dq2;
    Real d2q2 = 3 * dq3;

    //discriminant for finding the roots of the second derivative
    Real discr = sqr(d2q1) - 4. * d2q2 * d2q0;

    if (std::fabs(d2q2) < eps_Real * maxnrmq) {
      //only a cubic
      if (std::fabs(dq3) < eps_Real * maxnrmq) {
        //only a quadratic
        if (std::fabs(q2) < eps_Real * maxnrmq) {
          //only linear
          //std::cout<<" linear "<<std::endl;
          rootsearch = false;
          if (std::fabs(q1) < eps_Real * maxnrmq) {
            //constant;
            stepval = stepsize;
          } else {
            if (q1 < 0.)
              stepval = stepsize;
            else {
              stepval = min(-q0 / q1, stepsize);
            }
          }
        } else {
          //full quadratic
          //std::cout<<" quadratic "<<std::endl;
          if (q2 < 0) {  //concave case
            if (q1 <= 0.) {
              ubstep = stepsize;
              stepval = stepsize;
              rootsearch = false;
            } else {
              Real bx = -q1 / q2 / 2.;  //>0.
              Real bv = (q2 * bx + q1) * bx + q0;
              if (bv < 0) {
                ubstep = stepsize;
                stepval = stepsize;
                rootsearch = false;
              } else {
                ubstep = min(stepsize, bx);
                stepval = ubstep / 2.;
                rootsearch = true;
              }
            }
          } else {  //convex case
            Real bx = -q1 / q2 / 2.;  //minimum
            if (bx >= stepsize) {
              ubstep = stepsize;
              stepval = stepsize;
              rootsearch = false;
            } else {
              ubstep = stepsize;
              stepval = stepsize;
              rootsearch = true;
            }
          }
        }
      } else { //full cubic,
        Real infl = -dq1 / dq2 / 2.; //inflection point
        Real qinfl = ((q3 * infl + q2) * infl + q1) * infl + q0;
        //find zeros of derivative
        discr = sqr(dq1) - 4. * dq2 * dq0;
        //std::cout<<" cubic discr="<<discr<<std::endl;
        if (discr <= abseps) { //no zeros or a single zero of the derviative, function is monotone
          if (q3 < 0) {
            ubstep = stepsize;
            stepval = stepsize;
            rootsearch = false;
          } else {
            if ((infl >= stepsize) || ((infl > 0.) && (qinfl >= 0.))) {
              ubstep = min(stepsize, infl);
              stepval = 0.;
              rootsearch = true;
            } else {
              ubstep = min(stepsize, infl);
              stepval = stepsize;
              rootsearch = true;
            }
          }
        } //end no extrema
        else {  //minimum and maximum of full cubic exist
          Real x0 = (-dq1 - discr) / 2. / dq2;
          Real x1 = (-dq1 + discr) / 2. / dq2;
          Real qx0 = ((q3 * x0 + q2) * x0 + q1) * x0 + q0;
          Real qx1 = ((q3 * x1 + q2) * x1 + q1) * x1 + q0;
          //std::cout<<" fcubic x0="<<x0<<" qx0="<<qx0<<" x1="<<x1<<" qx1="<<qx1<<std::endl;
          if (q3 > 0) { //grows towards infinity, x0 is maximum
            if ((x0 >= stepsize) || ((x0 > 0.) && (qx0 >= 0.))) {
              ubstep = min(stepsize, infl);
              stepval = 0.;
              rootsearch = true;
            } else if (x1 >= stepsize) {
              ubstep = stepsize;
              stepval = stepsize;
              rootsearch = false;
            } else {
              ubstep = stepsize;
              stepval = stepsize;
              rootsearch = true;
            }
          } else { //q3<0., grows towards minus infinity, x1 is maxium
            if (x0 >= stepsize) {
              ubstep = stepsize;
              stepval = stepsize;
              rootsearch = false;
            } else if ((x1 >= stepsize) || ((x1 > 0.) && (qx1 >= 0.))) {
              ubstep = min(x1, stepsize);
              stepval = max(0., min(ubstep, infl));
              rootsearch = true;
            } else {
              ubstep = stepsize;
              stepval = stepsize;
              rootsearch = false;
            }
          }
        }
      } //end full cubic
    } //end d2q2 zero
    else if (discr <= eps_Real * maxnrmq) {
      //full 4th degree with one extremum
      //convex or concave case
      if (q4 > 0) {  //convex
        stepval = stepsize;
        ubstep = stepsize;
        Real qone = (((q4 * stepval + q3) * stepval + q2) * stepval + q1) * stepval + q0;
        Real dqone = ((dq3 * stepval + dq2) * stepval + dq1) * stepval + dq0;
        if (qone < 0.) {
          //value at full step size  is feasible
          rootsearch = false;
        } else {
          //value at fall step size is insfeasible
          if (dqone > 0) {
            //positive slope, root search will find zero 
            rootsearch = true;
          } else {
            //negative slope up to step size 1, so 1 is the least infeasible choice
            rootsearch = false;
          }
        }
        //std::cout<<" cvx discr="<<discr<<" q(1)="<<qone<<" q'(1)="<<dqone<<" search="<<rootsearch<<std::endl;
      } else {
        //std::cout<<" ccv discr="<<discr<<", should not happen"<<std::endl;
        // q4<0, concave
        if (dq0 <= 0.) {
          //negative slope in dq0 will stay so to step size 1, 1 is the best choice
          ubstep = stepsize;
          stepval = stepsize;
          rootsearch = false;
        } else {
          //positive solpe in 0 find the maximum, i.e. the zero of the derivative
          Real x = 0;
          Real der = dq0;
          //Newton should not overshoot by strict concavity
          while ((x < 1. - 1e-6) && (der > 1e-6 * maxnrmq)) {
            x += -der / ((d2q2 * x + d2q1) * x + d2q0);
            x = min(x, stepsize);
            der = ((dq3 * x + dq2) * x + dq1) * x + dq0;
          }
          Real qx = (((q4 * x + q3) * x + q2) * x + q1) * x + q0;
          if (qx >= 0.) {
            //the max is positive, start searching for value 0 from step size 0
            ubstep = min(stepsize, qx);
            stepval = 0;
            rootsearch = true;
          } else {
            //the max is negative [0,1], so take the maxium step 
            ubstep = stepsize;
            stepval = stepsize;
            rootsearch = false;
          }
        }
      }
    } else {
      // two inflection points
      discr = std::sqrt(discr);
      Real infl0 = (-d2q1 - discr) / 2. / d2q2;
      Real infl1 = (-d2q1 + discr) / 2. / d2q2;
      if (d2q2 < 0.)
        swap(infl0, infl1);
      Real qinfl0 = (((q4 * infl0 + q3) * infl0 + q2) * infl0 + q1) * infl0 + q0;
      Real qinfl1 = (((q4 * infl1 + q3) * infl1 + q2) * infl1 + q1) * infl1 + q0;
      Real d1infl0 = ((dq3 * infl0 + dq2) * infl0 + dq1) * infl0 + dq0;
      Real d1infl1 = ((dq3 * infl1 + dq2) * infl1 + dq1) * infl1 + dq0;

      // Real d2infl0=(d2q2*infl0+d2q1)*infl0+d2q0;
      // Real d2infl1=(d2q2*infl1+d2q1)*infl1+d2q0;

      //find the roots of the derivative which is a cubic 
      Real d0 = dq2 * dq2 - 3 * dq3 * dq1;
      Real d1 = 2 * dq2 * dq2 * dq2 - 9 * dq3 * dq2 * dq1 + 27 * dq3 * dq3 * dq0;
      Real dprime = d1 * d1 - 4 * d0 * d0 * d0;
      Real r_dp = 0.;
      Real phi_dp = 0.;
      if (dprime > 0) {
        r_dp = sqrt(dprime);
        if (d1 < 0) {
          r_dp -= d1;
          phi_dp = 3.1415926535897932;
        } else {
          r_dp += d1;
        }
      } else {
        r_dp = std::sqrt(-dprime + d1 * d1);
        phi_dp = std::acos(d1 / r_dp);
      }
      r_dp /= 2.;
      r_dp = pow(r_dp, 1 / 3.);
      phi_dp /= 3.;
      Real x0 = (dq2 + std::cos(phi_dp) * r_dp + (d0 / r_dp) * std::cos(-phi_dp)) / (-3. * dq3);
      Real Im_x0 = std::sin(phi_dp) * r_dp + (d0 / r_dp) * std::sin(-phi_dp);
      phi_dp += 2. * 3.1415926535897932 / 3.;
      Real x1 = (dq2 + std::cos(phi_dp) * r_dp + (d0 / r_dp) * std::cos(-phi_dp)) / (-3. * dq3);
      Real Im_x1 = std::sin(phi_dp) * r_dp + (d0 / r_dp) * std::sin(-phi_dp);
      phi_dp += 2. * 3.1415926535897932 / 3.;
      Real x2 = (dq2 + std::cos(phi_dp) * r_dp + (d0 / r_dp) * std::cos(-phi_dp)) / (-3. * dq3);
      Real Im_x2 = std::sin(phi_dp) * r_dp + (d0 / r_dp) * std::sin(-phi_dp);
      if (x0 > x1) {
        swap(x0, x1);
        swap(Im_x0, Im_x1);
      }
      if (x2 < x0) {
        Real d = x2; x2 = x1; x1 = x0; x0 = d;
        d = Im_x2; Im_x2 = Im_x1; Im_x1 = Im_x0; Im_x0 = d;
      } else if (x2 < x1) {
        swap(x1, x2);
        swap(Im_x1, Im_x2);
      }
      Real qx0 = (((q4 * x0 + q3) * x0 + q2) * x0 + q1) * x0 + q0;
      Real qx1 = (((q4 * x1 + q3) * x1 + q2) * x1 + q1) * x1 + q0;
      Real qx2 = (((q4 * x2 + q3) * x2 + q2) * x2 + q1) * x2 + q0;
      // Real d1x0=((dq3*x0+dq2)*x0+dq1)*x0+dq0;
      // Real d1x1=((dq3*x1+dq2)*x1+dq1)*x1+dq0;
      // Real d1x2=((dq3*x2+dq2)*x2+dq1)*x2+dq0;
      // std::cout.precision(12);
      // std::cout<<" x0="<<x0<<" ("<<qx0<<",i"<<Im_x0<<","<<d1x0<<")";
      // std::cout<<" infl0="<<infl0<<" ("<<qinfl0<<" "<<d1infl0<<","<<d2infl0<<")";
      // std::cout<<" x1="<<x1<<" ("<<qx1<<",i"<<Im_x1<<","<<d1x1<<")";
      // std::cout<<" infl1="<<infl1<<" ("<<qinfl1<<" "<<d1infl1<<","<<d2infl1<<")";
      // std::cout<<" x2="<<x2<<" ("<<qx2<<",i"<<Im_x2<<","<<d1x2<<")";
      // std::cout<<std::endl;

      // assert((std::fabs(Im_x0)>1e-8)||(x0<=infl0+1e-10*(1.+std::fabs(x0))));
      // assert((std::fabs(Im_x1)>1e-8)||(infl0<=x1+1e-10*(1.+std::fabs(infl0))));
      // assert((std::fabs(Im_x1)>1e-8)||(x1<=infl1+1e-10*(1.+std::fabs(x1))));
      // assert((std::fabs(Im_x2)>1e-8)||(infl1<=x2+1e-10*(1.+std::fabs(infl1))));


      if (q4 > 0) {
        if (d1infl0 <= 0) { //both infl have negative deriv, min is in x2
          if (x2 >= 1.) {
            ubstep = stepsize;
            stepval = stepsize;
            rootsearch = false;
          } else {
            ubstep = stepsize;
            stepval = stepsize;
            rootsearch = true;
          }
        } else if (d1infl1 >= 0) { //both infl have positive deriv, min is in x0
          if (x0 >= 1.) {
            ubstep = stepsize;
            stepval = stepsize;
            rootsearch = false;
          } else if ((infl1 >= stepsize) || ((infl1 > 0.) && (qinfl1 >= 0.))) {
            ubstep = min(stepsize, infl1);
            stepval = max(0., min(ubstep, infl0));
            rootsearch = true;
          } else {
            ubstep = stepsize;
            stepval = stepsize;
            rootsearch = true;
          }
        } else { //there is a max in the middle
          assert(std::fabs(Im_x0) < 1e-8 * maxnrmq);
          assert(std::fabs(Im_x1) < 1e-8 * maxnrmq);
          assert(std::fabs(Im_x2) < 1e-8 * maxnrmq);
          if (x0 >= stepsize) {
            ubstep = stepsize;
            stepval = stepsize;
            rootsearch = false;
          } else if ((x1 >= stepsize) || ((x1 > 0.) && (qx1 >= 0.))) {
            ubstep = min(x1, stepsize);
            stepval = max(0., min(ubstep, infl0));
            rootsearch = true;
          } else if (x2 >= stepsize) {
            ubstep = stepsize;
            stepval = stepsize;
            rootsearch = false;
          } else {
            ubstep = stepsize;
            stepval = stepsize;
            rootsearch = true;
          }
        }
      } else {  //q4<0
        if (d1infl1 <= 0) { //both infl have negative derivatives, max in x0
          if ((x0 >= stepsize) || ((x0 > 0.) && (qx0 > 0.))) {
            ubstep = min(stepsize, x0);
            stepval = 0.;
            rootsearch = true;
          } else {
            ubstep = stepsize;
            stepval = stepsize;
            rootsearch = false;
          }
        } else if (d1infl0 >= 0) {//both infl have positive derivatives, max in x2
          if ((infl0 >= stepsize) || ((infl0 > 0.) && (qinfl0 >= 0.))) {
            ubstep = min(stepsize, infl0);
            stepval = 0.;
            rootsearch = true;
          } else if ((x2 >= stepsize) || ((x2 > 0.) && (qx2 >= 0.))) {
            ubstep = min(stepsize, x2);
            stepval = max(0., min(ubstep, infl1));
            rootsearch = true;
          } else {
            ubstep = stepsize;
            stepval = stepsize;
            rootsearch = false;
          }
        } else { //there is a minimum between the inflection points
          assert(std::fabs(Im_x0) < 1e-8 * maxnrmq);
          assert(std::fabs(Im_x1) < 1e-8 * maxnrmq);
          assert(std::fabs(Im_x2) < 1e-8 * maxnrmq);
          if ((x0 >= stepsize) || ((x0 > 0.) && (qx0 >= 0.))) {
            ubstep = min(stepsize, x0);
            stepval = 0.;
            rootsearch = true;
          } else if (x1 >= stepsize) {
            ubstep = stepsize;
            stepval = stepsize;
            rootsearch = false;
          } else if ((x2 >= stepsize) || ((x2 > 0.) && (qx2 >= 0.))) {
            ubstep = min(stepsize, x2);
            stepval = max(0., min(ubstep, infl1));
            rootsearch = true;
          } else {
            ubstep = stepsize;
            stepval = stepsize;
            rootsearch = false;
          }
        } //end there is min in the middle
      } // end q4<0.
    } // end two inflection points

    Real qsv = (((q4 * stepval + q3) * stepval + q2) * stepval + q1) * stepval + q0;
    Real dqsv = ((dq3 * stepval + dq2) * stepval + dq1) * stepval + dq0;
    if (rootsearch) {
      Real lbstep = 0.;
      do {
        //std::cout<<" rs["<<lbstep<<","<<ubstep<<"]("<<qsv<<","<<dqsv<<")";
        assert(ubstep <= stepsize + 1.e-10);
        assert(lbstep <= stepval);
        assert(stepval <= ubstep + 1.e-10);
        if (qsv <= 0.) {
          lbstep = stepval;
          if (dqsv > 0.) {
            stepval += min(-qsv / dqsv, 0.999 * (ubstep - stepval));
          } else {
            stepval += 0.999 * (ubstep - stepval);
          }
        } else {
          if (dqsv > 0.) {
            ubstep = stepval;
            stepval += max(-qsv / dqsv, 0.999 * (lbstep - stepval));
          } else {
            if (cb_out()) {
              get_out() << "**** WARNING in InteriorPointBlock::pol_le_zero_step(): root search should never get here,";
              get_out() << " lbstep=" << lbstep << " stepval=" << stepval << " ubstep=" << ubstep << " qsv=" << qsv << " dqsv=" << dqsv << std::endl;
            }
            break;
          }
        }
        qsv = (((q4 * stepval + q3) * stepval + q2) * stepval + q1) * stepval + q0;
        dqsv = ((dq3 * stepval + dq2) * stepval + dq1) * stepval + dq0;
      } while ((std::fabs(qsv) > abseps) && (ubstep - lbstep > 1e-8));
      //std::cout<<" frs["<<lbstep<<","<<ubstep<<"]("<<qsv<<","<<dqsv<<")";
    }

    stepsize = min(stepsize, max(0., stepval));
    return 0;
  }

  int InteriorPointBlock::minimize_pol_step(Real& stepsize,
    Real q0,
    Real q1,
    Real q2,
    Real q3,
    Real q4,
    Real /* abseps */) const {
    //ploynomial coefficients of the derivative dq, need to find zeros of this in [0, stepval]
    Real dq0 = q1;
    Real dq1 = 2 * q2;
    Real dq2 = 3 * q3;
    Real dq3 = 4 * q4;

    Real bestval = q0;
    Real bestx = 0.;
    Real valstep = q0 + stepsize * (q1 + stepsize * (q2 + stepsize * (q3 + stepsize * q4)));
    if (valstep < bestval) {
      bestval = valstep;
      bestx = stepsize;
    }

    //find the roots of the derivative which is a cubic 
    Real d0 = dq2 * dq2 - 3 * dq3 * dq1;
    Real d1 = 2 * dq2 * dq2 * dq2 - 9 * dq3 * dq2 * dq1 + 27 * dq3 * dq3 * dq0;
    Real dprime = d1 * d1 - 4 * d0 * d0 * d0;
    Real r_dp = 0.;
    Real phi_dp = 0.;
    if (dprime > 0) {
      r_dp = sqrt(dprime);
      if (d1 < 0) {
        r_dp -= d1;
        phi_dp = 3.1415926535897932;
      } else {
        r_dp += d1;
      }
    } else {
      r_dp = std::sqrt(-dprime + d1 * d1);
      phi_dp = std::acos(d1 / r_dp);
    }
    r_dp /= 2.;
    r_dp = pow(r_dp, 1 / 3.);
    phi_dp /= 3.;
    Real x0 = (dq2 + std::cos(phi_dp) * r_dp + (d0 / r_dp) * std::cos(-phi_dp)) / (-3. * dq3);
    Real Im_x0 = std::sin(phi_dp) * r_dp + (d0 / r_dp) * std::sin(-phi_dp);
    phi_dp += 2. * 3.1415926535897932 / 3.;
    Real x1 = (dq2 + std::cos(phi_dp) * r_dp + (d0 / r_dp) * std::cos(-phi_dp)) / (-3. * dq3);
    Real Im_x1 = std::sin(phi_dp) * r_dp + (d0 / r_dp) * std::sin(-phi_dp);
    phi_dp += 2. * 3.1415926535897932 / 3.;
    Real x2 = (dq2 + std::cos(phi_dp) * r_dp + (d0 / r_dp) * std::cos(-phi_dp)) / (-3. * dq3);
    Real Im_x2 = std::sin(phi_dp) * r_dp + (d0 / r_dp) * std::sin(-phi_dp);
    if (x0 > x1) {
      swap(x0, x1);
      swap(Im_x0, Im_x1);
    }
    if (x2 < x0) {
      Real d = x2; x2 = x1; x1 = x0; x0 = d;
      d = Im_x2; Im_x2 = Im_x1; Im_x1 = Im_x0; Im_x0 = d;
    } else if (x2 < x1) {
      swap(x1, x2);
      swap(Im_x1, Im_x2);
    }
    Real qx0 = (((q4 * x0 + q3) * x0 + q2) * x0 + q1) * x0 + q0;
    if ((0 < x0) && (x0 < stepsize) && (qx0 < bestval)) {
      bestval = qx0;
      bestx = x0;
    }
    Real qx1 = (((q4 * x1 + q3) * x1 + q2) * x1 + q1) * x1 + q0;
    if ((0 < x1) && (x1 < stepsize) && (qx1 < bestval)) {
      bestval = qx1;
      bestx = x1;
    }
    Real qx2 = (((q4 * x2 + q3) * x2 + q2) * x2 + q1) * x2 + q0;
    if ((0 < x2) && (x2 < stepsize) && (qx2 < bestval)) {
      bestval = qx2;
      bestx = x2;
    }
    stepsize = bestx;

    return 0;
  }

  int InteriorPointBlock::control_nbh_step(Real& stepsize,
    Real& max_nbh,
    Real nbh_ubnd,
    Real mu_xz,
    Real mu_xdzpdxz,
    Real mu_dxdz,
    Real nrmsqr_xz,
    Real nrmsqr_xdzpdxz,
    Real nrmsqr_dxdz,
    Real ip_xz_xdzpdxz,
    Real ip_xz_dxdz,
    Real ip_dxdz_xdzpdxz) const {
    const Real minstep = 1000. * eps_Real;
    if (stepsize < minstep)
      return 0;

    int err = 0;
    Real otheta = std::sqrt(nrmsqr_xz) / mu_xz;
    Real theta = max(nbh_ubnd, otheta);
    Real q0 = nrmsqr_xz - sqr(theta * mu_xz);
    Real q1 = 2. * ip_xz_xdzpdxz - 2 * theta * theta * mu_xz * mu_xdzpdxz;
    Real q2 = nrmsqr_xdzpdxz + 2. * ip_xz_dxdz - theta * theta * (mu_xdzpdxz * mu_xdzpdxz + 2 * mu_xz * mu_dxdz);
    Real q3 = 2. * ip_dxdz_xdzpdxz - 2 * theta * theta * mu_xdzpdxz * mu_dxdz;
    Real q4 = nrmsqr_dxdz - sqr(theta * mu_dxdz);
    Real q_infnrm = max(std::fabs(q0), max(std::fabs(q1), max(std::fabs(q2), max(std::fabs(q3), std::fabs(q4)))));

    if (cb_out(2)) {
      get_out().precision(4);
      get_out() << " nbh ina=" << stepsize;
    }


    if ((otheta > 1 - 1e-8) || (q_infnrm < 1e-4 * max(mu_xz, otheta))) {
      int retval = linesearch(stepsize);
      if (retval) {
        if (cb_out()) {
          get_out() << "*** ERROR InteriorPointBlock::control_nbh_step(): linesearch(.) returned " << retval << std::endl;
        }
        err++;
      }
      if (cb_out(2)) {
        get_out() << " lsalpha=" << stepsize;
      }
    }

    if ((stepsize > minstep) && (q_infnrm > .99e-4 * max(mu_xz, otheta))) {
      const Real abseps = 1e-8 * mu_xz;
      bool use_minimize_theta = true;
      bool skip_search = false;
      if (otheta < (1 + 1e-6) * nbh_ubnd) {
        //the current point is safely inside, make the step as long as possible
        use_minimize_theta = false;
        if ((std::fabs(q1) < abseps)
          && (std::fabs(q2) < abseps)
          && (std::fabs(q3) < abseps)
          && (std::fabs(q4) < abseps)
          ) {
          //polynomial is constant
          skip_search = true;
        }
      } else {
        //theta sill too large
        //by just using the linear approximation, the following theta could be achieved
        Real gtheta = std::sqrt(max(0., otheta * otheta + 2 * (ip_xz_xdzpdxz - otheta * otheta * mu_xz * mu_xdzpdxz) / max((mu_xz + mu_xdzpdxz + mu_dxdz), 1e-6 * mu_xz)));
        if (gtheta < .1 * nbh_ubnd + 0.9 * otheta) {
          use_minimize_theta = false;
          gtheta = max(.1 * gtheta + .9 * otheta, nbh_ubnd);
          Real mu_c1 = gtheta * mu_xdzpdxz + (gtheta - theta) * (mu_xz + mu_dxdz);
          q1 = 2. * ip_xz_xdzpdxz - 2 * theta * mu_xz * mu_c1;
          assert(q1 < 0.);
          q2 = nrmsqr_xdzpdxz + 2. * ip_xz_dxdz - (mu_c1 * mu_c1 + 2 * theta * theta * mu_xz * mu_dxdz);
          q3 = 2. * ip_dxdz_xdzpdxz - 2 * theta * mu_c1 * mu_dxdz;
          q4 = nrmsqr_dxdz - sqr(theta * mu_dxdz);
          if (cb_out(2)) {
            get_out() << " gnbh=" << gtheta;
          }
        } else if (q1 > 1e-8 * mu_xz) {
          stepsize = 0.;
          skip_search = true;
        } else {
          //q1 "=" 0
          if (q2 > 1e-8 * mu_xz) {
            //curves upwards, no step
            stepsize = 0;
            skip_search = true;
          } else if (q2 > -1e-8 * mu_xz) {
            //q2 "=" 0
            if (q3 > 1e-8 * mu_xz) {
              //curves upwards, no step
              stepsize = 0;
              skip_search = true;
            } else if (q3 > -1e-8 * mu_xz) {
              //q3 "=" 0
              if (q4 > 1e-8 * mu_xz) {
                //curves upwards, no step
                stepsize = 0;
                skip_search = true;
              } else if (q4 > -1e-8 * mu_xz) {
                //polynomial is zero, just ignore it and do not restrict the step further
                skip_search = true;
              }
            }
          }
        }
      } // end theta too large

      if (cb_out(2)) {
        get_out() << " skips=" << skip_search;
        get_out() << " min=" << use_minimize_theta;
        get_out() << " q0=" << q0;
        get_out() << " q1=" << q1;
        get_out() << " q2=" << q2;
        get_out() << " q3=" << q3;
        get_out() << " q4=" << q4;
      }

      if (!skip_search) {
        //determine the largest step possible
        int retval = pol_le_zero_step(stepsize, q0, q1, q2, q3, q4, 1e-8 * mu_xz);
        if (retval) {
          if (cb_out()) {
            get_out() << "*** ERROR InteriorPointBlock::get_nbh_info(): pol_le_zero_step() returned " << retval << std::endl;
          }
          err++;
        }
        if (use_minimize_theta) {
          //within this stepsize find the minimum
          retval = minimize_pol_step(stepsize, q0, q1, q2, q3, q4, 1e-8 * mu_xz);
          if (retval) {
            if (cb_out()) {
              get_out() << "*** ERROR InteriorPointBlock::get_nbh_info(): pol_le_zero_step() returned " << retval << std::endl;
            }
            err++;
          }
        }
      } //endif (!skip_search) after building the neighborhood polynomial
    } //endif (stepsize>minstep) after line search

    Real step_nbh = std::sqrt(std::fabs(nrmsqr_xz + stepsize * (2 * ip_xz_xdzpdxz + stepsize * (nrmsqr_xdzpdxz + 2 * ip_xz_dxdz + stepsize * (2 * ip_dxdz_xdzpdxz + stepsize * nrmsqr_dxdz))))) / max(mu_xz + stepsize * (mu_xdzpdxz + stepsize * (mu_dxdz)), 1e-6 * mu_xz);
    max_nbh = max(max_nbh, step_nbh);

    if (cb_out(2)) {
      get_out() << " outa=" << stepsize << "(" << step_nbh << ")";
    }

    return err;
  }

}
