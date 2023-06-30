/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/UQPSolver.cxx
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



#include <iomanip>
#include <sstream>
#include <fstream>
#include "UQPSolver.hxx"


using namespace CH_Matrix_Classes;

namespace ConicBundle {


  void UQPSolver::set_defaults() {
    termeps = 1e-7;
    maxiter = 100;
  }

  void UQPSolver::clear() {
    clear_model_data_ptr();
    Q.init(0, 0.);
    c.init(0, 1, 0.);
    offset = 0.;
    A.init(0, 0, 0.);
    b.init(0, 1, 0.);
    lowerbound = 0;
    upperbound = 0;
    x.init(0, 1, 0.);
    y.init(0, 1, 0.);
    primalval = 0.;
    dualval = 0.;
    mu = 0.;
    Qx.init(0, 0, 0.);
    iter = 0;
    status = 0;
    dx.init(0, 1, 0.);
    dy.init(0, 1, 0.);
    Qplus.init(0, 0.);
    LinvAt.init(0, 0, 0.);
    sysdy.init(0, 0.);
    rd.init(0, 1, 0.);
    rhs.init(0, 1, 0.);
    xcorr.init(0, 1, 0.);
    tmpvec.init(0, 1, 0.);

    sum_iter = 0;
    sum_choliter = 0;
    clock.start();
    sum_choltime = 0;
  }



  // *************************************************************************
  //                            select_new_mu
  // *************************************************************************

  void UQPSolver::select_new_mu(const Matrix& dx, const Matrix& dy, const Matrix& rhs_residual) {
    Real oldmu = mu;
    Real ip_xz;
    Real sigma;
    Integer mu_dim;
    model_block->suggest_mu(ip_xz, mu_dim, sigma, dx, dy, rhs_residual);
    sigma = max(0.01, sigma);
    sigma = min(1., sigma);
    mu = ip_xz / mu_dim;
    mu = min(oldmu, sigma * mu);
  }


  // *************************************************************************
  //                            predcorr_step
  // *************************************************************************

  // do a predictor corrector step 

  int UQPSolver::predcorr_step(Real& alpha) {
    status = 0;
    Indexmatrix piv;

    //--- build Q-part of system matrix
    Qplus = Q;
    if (model_block->add_xinv_kron_z(Qplus)) {
      if (cb_out(1)) {
        get_out() << "*** WARNING: UQPSolver::predcorr_step(): collecting system matrix failed " << std::endl;
      }
      status = 3;
      return status;
    }

    Symmatrix barQ(Qplus);                   //TEST

    //--- factorize Q-part of system matrix and compute system for y
    if (Qplus.Chol_factor(eps_Real)) {
      if (cb_out(1)) {
        get_out() << "*** WARNING: UQPSolver::predcorr_step(): factorizing system matrix failed" << std::endl;
      }
      //cout<<barQ<<std::endl;
      //cout<<Q<<std::endl;
      status = 2;
      return status;
    }

    //Symmatrix barQinv;                        //TEST
    //Matrix barQinvmat(Diag(Matrix(Qplus.dim(),1,1.)));//TEST
    //Qplus.Chol_solve(barQinvmat);//TEST
    //barQinv.init(barQinvmat);//TEST
    //std::cout<<"  norm2(barQinv*barQ-Diag(Matrix(Qplus.dim(),1,1.)))="<<norm2(barQinv*barQ-Diag(Matrix(Qplus.dim(),1,1.)))<<std::endl;//TEST

    if (y.dim() > 0) {
      //LinvAt=L^-1*transpose(A)
      xbpeya(LinvAt, A, 1., 0., 1);
      Qplus.Chol_Lsolve(LinvAt);

      //sysdy=LiAt'*LiAt
      rankadd(LinvAt, sysdy, 1., 0., 1);
    }
    //std::cout<<"  norm2(sysdy-A*barQinv*transpose(A))="<<norm2(sysdy-A*barQinv*transpose(A))<<std::endl;//TEST

    //rd=c-Qx-At*y;  
    xeyapzb(rd, c, Qx, 1., -1.);

    if (y.dim() > 0) {
      genmult(A, y, rd, -1., 1., 1);

      //std::cout<<"  norm2(rd-c+Qx+transpose(A)*y)="<<norm2(rd-c+Qx+transpose(A)*y)<<std::endl;//TEST

      //rhs= A*barQinv*rd-(b-Ax)
      tmpvec = rd;
      Qplus.Chol_Lsolve(tmpvec);
      genmult(LinvAt, tmpvec, rhs, 1., 0., 1);
      rhs -= b;
      genmult(A, x, rhs, 1., 1.);

      //std::cout<<"  norm2(rhs-A*barQinv*rd+b-A*x)="<<norm2(rhs-A*barQinv*rd+b-A*x)<<std::endl;//TEST


      //--- build dy-part of system matrix
      if (model_block->add_local_sys(sysdy, rhs)) {
        if (cb_out(1)) {
          get_out() << "*** WARNING: UQPSolver::predcorr_step(): adding local system for predictor failed" << std::endl;
        }
        status = 3;
        return status;
      }

      //Symmatrix testsysdy=sysdy; //TEST


      //--- solve predictor direction
      dy = rhs;
      sysdy.Chol_factor(piv, 1000. * eps_Real);
      if (piv.dim() < sysdy.rowdim()) {
        if (cb_out(1)) {
          get_out() << "*** WARNING: UQPSolver::predcorr_step(): factorizing for predictor failed" << std::endl;
        }
        status = 2;
        return status;
      }
      sysdy.Chol_solve(dy, piv);

      //std::cout<<"  norm2(testsysdy*dy-rhs)="<<norm2(testsysdy*dy-rhs)<<std::endl; //TEST

      //dx=barQinv*(c-At*y-Qx-At*dy)
      //dz=-(c-At*y-Qx-At*dy)-z+Q*dx;
      genmult(A, dy, dx, -1., 0., 1);
      dx += rd;
    } else {
      dy.init(0, 0, 0.);
      rhs.init(0, 0, 0.);
      dx = rd;
    }
    tmpvec.init(dx, -1.);
    Qplus.Chol_solve(dx);
    /// tmpvec = -(c-At*(y+dy)-Q(x+dx));
    genmult(Q, dx, tmpvec, 1., 1.);

    //--- select new mu (and compute dz)
    select_new_mu(dx, dy, tmpvec);

    // std::cout<<" norm2(barQ*dx-c+Qx+transpose(A)*(y+dy))="<<norm2(barQ*dx-c+Qx+transpose(A)*(y+dy))<<std::endl; //TEST
    // Matrix dual_respred=Q*(x+dx)-c+transpose(A)*(y+dy);
    // std::cout<<" norm2(Q*(x+dx)-c+transpose(A)*(y+dy)-z-dz)="<<norm2(model_block->subtract_z(dual_respred,true))<<" dual_respred="<<transpose(dual_respred)<<std::endl; //TEST
    //Symmatrix Xtest,Ztest,dX1test,dZ1test;//TEST 
    //Integer dim1=Integer(::sqrt(Real(8*x.dim()+1))-1+.1)/2;
    //if (x.dim()!=(dim1*(dim1+1))/2){//TEST
    //  sveci(x(Range(1,x.dim()-1)),Xtest);//TEST
    //  sveci(z(Range(1,z.dim()-1)),Ztest);//TEST
    //  sveci(dx(Range(1,x.dim()-1)),dX1test);//TEST
    //  sveci(dz(Range(1,x.dim()-1)),dZ1test);//TEST
    //}                //TEST
    //else {           //TEST
    //  sveci(x,Xtest);//TEST
    //  sveci(z,Ztest);//TEST
    //  sveci(dx,dX1test);//TEST
    //  sveci(dz,dZ1test);//TEST
    //}                  //TEST
    //std::cout<<" norm2(Symmatrix(Xtest*Ztest*dX1test)+Xtest*(Ztest+dZ1test)*Xtest)="<<norm2(Symmatrix(Xtest*Ztest*dX1test)+Xtest*(Ztest+dZ1test)*Xtest)<<std::endl; //TEST


    //--- compute corrector terms (and modify rhs for iternal correctors)
    xcorr.init(x.dim(), 1, 0.);
    if (model_block->get_corr(xcorr, rhs, mu)) {
      if (cb_out(1)) {
        get_out() << "*** WARNING: UQPSolver::predcorr_step(): getting local correctors failed" << std::endl;
      }
      status = 3;
      return status;
    }

    //--- solve corrector direction
    // dy = sysinv*( barQinv*xcorr+rhs) 
    if (y.dim() > 0) {
      tmpvec = xcorr;
      Qplus.Chol_Lsolve(tmpvec);
      genmult(LinvAt, tmpvec, dy, 1., 0., 1);
      dy += rhs;
      sysdy.Chol_solve(dy, piv);
      //dx=barQinv*(c-At*y-Qx-At*dy+xcorr)
      //tmpvec=-(c-At*y-Qx-At*dy)-z+Q*dx;
      genmult(A, dy, dx, -1., 0., 1);
      dx += rd;
    } else {
      dx = rd;
    }
    tmpvec.init(dx, -1.);
    dx += xcorr;
    Qplus.Chol_solve(dx);
    //tmpvec=-(c-At*y-Qx-At*dy)-z+Q*dx;
    genmult(Q, dx, tmpvec, 1., 1.);

    //--- compute dz and do the line search 
    alpha = 1.;
    if (model_block->line_search(alpha, dx, dy, tmpvec)) {
      if (cb_out(1)) {
        get_out() << "*** WARNING: UQPSolver::pred_corr_step(): line search failed" << std::endl;
      }
      status = 3;
      return status;
    }

    // std::cout<<" norm2(barQ*dx-c+Qx+transpose(A)*(y+dy)-xcorr)="<<norm2(barQ*dx-c+Qx+transpose(A)*(y+dy)-xcorr)<<std::endl; //TEST
    // Matrix dual_rescorr=Q*(x+dx)-c+transpose(A)*(y+dy);
    // std::cout<<" norm2(Q*(x+dx)-c+transpose(A)*(y+dy)-z-dz)="<<norm2(model_block->subtract_z(dual_rescorr,true))<<" drcorr="<<transpose(dual_rescorr)<<std::endl; //TEST
    //Symmatrix dX2test,dZ2test;//TEST 
    //Integer dim2=Integer(::sqrt(Real(8*dx.dim()+1))-1+.1)/2;
    //if (x.dim()!=(dim2*(dim2+1))/2){//TEST
    // sveci(dx(Range(1,x.dim()-1)),dX2test);//TEST
    //  sveci(dz(Range(1,x.dim()-1)),dZ2test);//TEST
    //}                //TEST
    //else {           //TEST
    //  sveci(dx,dX2test);//TEST
    //  sveci(dz,dZ2test);//TEST
    //}                  //TEST
    //std::cout<<"  norm2(Symmatrix(Xtest*Ztest*dX2test+dX1test*dZ1test*Xtest)+Xtest*(Ztest+dZ2test)*Xtest-mu*Xtest)="<<norm2(Symmatrix(Xtest*Ztest*dX2test+dX1test*dZ1test*Xtest)+Xtest*(Ztest+dZ2test)*Xtest-mu*Xtest)<<std::endl; //TEST


    return status;
  }

  // *************************************************************************
  //                            iterate
  // *************************************************************************

  // loop till convergence to optimal solution

  int UQPSolver::iterate() {
    //std::cout<<" norm2(A*x-b)="<<norm2(A*x-b);          //TEST
    //std::cout<<" norm2(Q*x-c+transpose(A)*y-z)="<<norm2(Q*x-c+transpose(A)*y-z)<<std::endl; //TEST
    //cout<<"Q=";Q.display(std::cout);
    //cout<<"c=";c.display(std::cout);
    //cout<<"A=";A.display(std::cout);
    //cout<<"b=";b.display(std::cout);
    //cout<<"x=";x.display(std::cout);
    //cout<<"y=";y.display(std::cout);
    //cout<<"z=";z.display(std::cout);

   // display parameters
    if (cb_out(1)) {
      get_out() << "\n    mi=" << maxiter;
      get_out() << " te="; get_out().precision(5); get_out() << termeps;
      get_out() << " lb="; get_out().precision(16); get_out().width(11); get_out() << lowerbound;
      get_out() << " ub="; get_out().precision(16); get_out().width(11); get_out() << upperbound; get_out().precision(2);
      get_out() << " Qn=" << Q.rowdim();
      get_out() << std::endl;
    }

    // compute primal and dual objective value
    // (x,y,z,and Qx are already computed)

    Real xtQxo2 = ip(x, Qx) / 2.;

    // x.transpose(); std::cout<<"x= ";x.display(std::cout); x.transpose(); //TEST
    //std::cout<<"A= ";A.display(std::cout); //TEST
    //std::cout<<"b= ";b.display(std::cout); //TEST
    // std::cout<<" xtQxo2="<<xtQxo2; //TEST
    // std::cout<<" ip(c,x)="<<ip(c,x); //TEST
    // std::cout<<" offset="<<offset<<std::endl; //TEST
    //Matrix AxBs=A*x; //TEST
    //model_block->add_Bs(AxBs); //TEST
    //std::cout<<" Fx=";get_out().width(7);get_out()<<norm2(b-AxBs); //TEST
    //std::cout<<" Fz=";get_out().width(7);get_out()<<norm2(c+z-Qx-transpose(A)*y); //TEST

    dualval = xtQxo2 + ip(b, y) + offset;
    dualval += model_block->get_local_dualcost();
    primalval = -xtQxo2 + ip(c, x) + offset;
    primalval += model_block->get_local_primalcost();
    Real ub = upperbound;
    Real lb = lowerbound;
    if (lb <= primalval) lb = primalval;

    //output
    if (cb_out(1)) {
      get_out().precision(2);
      get_out() << "    "; get_out().width(2); get_out() << iter << ":";
      get_out() << " gappd="; get_out().width(7); get_out() << dualval - primalval;
      get_out() << " pv="; get_out().precision(16); get_out().width(10); get_out() << primalval;
      get_out() << " dv="; get_out().precision(16); get_out().width(10); get_out() << dualval;
      get_out().precision(2);
      get_out() << " epsdl="; get_out().width(7); get_out() << .5 * (dualval - lb);
      get_out() << " epsup="; get_out().width(7); get_out() << termeps * (ub - primalval);
      get_out() << std::endl;
    }

    Real alpha = 1.;

    while (
      ((maxiter < 0) || (iter < maxiter))
      && (
        ((dualval > lb + eps_Real * (fabs(lb) + 1.)) && (primalval < lb + .5 * (dualval - lb)))
        ||
        (
          (dualval - lb > 1e-12 * (fabs(dualval) + 1.)) &&
          (ub - primalval > 1e-12 * (fabs(primalval) + 1.)) &&
          (alpha > 1e-8) &&
          (dualval - primalval > min(.5 * (dualval - lb), termeps * (ub - primalval)))
          )
        )
      ) {

      iter++;
      sum_iter++;

      //--- find the step direction
      sum_choliter++;
      CH_Tools::Microseconds st = clock.time();
      status = predcorr_step(alpha);
      sum_choltime += clock.time() - st;
      if (status) {
        if (cb_out(1)) {
          get_out() << "*** WARNING: UQPSolver::predcorr_step() returned " << status << std::endl;
        }
        return status;
      }


      //--- set new point with safeguard when restarting if first stepsizes get small 
      bool restart = false;
      if ((!run_starting_point) && (alpha < 1e-3) && (iter <= 3)) {
        //restart from scratch
        //initialize starting point (dimensions are ok already)
        restart = true;
        run_starting_point = true;
        model_block->starting_x(x);
        genmult(Q, x, Qx);   //Qx=Q*x
        model_block->starting_y(y, Qx, c);

        alpha = 1.;
      } else {
        //--- x+=alpha*dx, y+=alpha*dy, z+=alpha*dz, set new point
        x.xpeya(dx, alpha);
        y.xpeya(dy, alpha);

        if (model_block->set_point(x, y, alpha)) {
          if (cb_out(1)) {
            get_out() << "*** WARNING: UQPSolver::iterate(): setting new point failed" << std::endl;
          }
          status = 3;
          return status;
        }

        genmult(Q, x, Qx);
      }

      Real xtQxo2 = ip(x, Qx) / 2.;

      dualval = xtQxo2 + ip(b, y) + offset;
      dualval += model_block->get_local_dualcost();
      primalval = -xtQxo2 + ip(c, x) + offset;
      primalval += model_block->get_local_primalcost();

      //output
      if (cb_out(1)) {
        // if (cb_out(2)){
        //   get_out().precision(12);
        //   get_out()<<"     xtQxt/2="<<xtQxo2;
        //   get_out()<<"  ip(c,x)="<<ip(c,x);
        //   get_out()<<"  offset="<<offset;
        //   get_out()<<"  loccsot="<<model_block->get_local_primalcost();
        //   get_out()<<std::endl;
        // }

        get_out().precision(2);

        // Matrix AxBs=A*x; //TEST
        // model_block->add_Bs(AxBs); //TEST
        // std::cout<<" Fx=";get_out().width(7);get_out()<<norm2(b-AxBs); //TEST
        // Matrix dual_res=Qx+transpose(A)*y-c;
        // std::cout<<" Fz=";get_out().width(7);get_out()<<norm2(model_block->subtract_z(dual_res)); //TEST

        if (restart)
          get_out() << " [r]";
        else
          get_out() << "    ";
        get_out().width(2); get_out() << iter << ":";
        get_out() << " gappd="; get_out().width(7); get_out() << dualval - primalval;
        get_out() << " pv="; get_out().precision(16); get_out().width(10); get_out() << primalval;
        get_out() << " dv="; get_out().precision(16); get_out().width(10); get_out() << dualval;
        get_out().precision(2);
        get_out() << " epsdl="; get_out().width(7); get_out() << .5 * (dualval - lb);
        get_out() << " epsup="; get_out().width(7); get_out() << termeps * (ub - primalval);
        get_out() << " mu="; get_out().width(7); get_out() << mu;
        get_out() << " alpha="; get_out().width(4); get_out() << alpha;
        get_out() << std::endl;
      }

    } //end while

    if (cb_out(1)) {
      get_out() << "    term: iter=";
      if ((maxiter > 0) && (iter > maxiter))
        get_out() << "T";
      else
        get_out() << "F";
      get_out() << " dv-lb=";
      if (dualval - lb <= 1e-12 * (fabs(dualval) + 1.))
        get_out() << "T";
      else
        get_out() << "F";
      get_out() << " ub-pv=";
      if (ub - primalval <= 1e-12 * (fabs(primalval) + 1.))
        get_out() << "T";
      else
        get_out() << "F";
      get_out() << " alpha=";
      if (alpha <= 1e-8)
        get_out() << "T";
      else
        get_out() << "F";
      get_out() << " dv-pv=";
      if (dualval - primalval <= min(.5 * (dualval - lb), termeps * (ub - primalval)))
        get_out() << "T";
      else
        get_out() << "F";
      get_out() << std::endl;
    }

    if ((maxiter < 0) || (iter < maxiter)) {
      status = 0;
    } else {
      status = 1;
    }

    return  status;
  }

  // *************************************************************************
  //                             solve
  // *************************************************************************

  // call starting_point and loop till convergence to optimal solution

  int UQPSolver::solve(const Symmatrix& Qin, const Matrix& cin, Real offsetin) {
    Q = Qin;
    c = cin;
    offset = offsetin;

    run_starting_point = true;

    //initialize dimensions and primal starting point
    x.init(Q.rowdim(), 1, 0.);
    Integer ydim = 0;
    Integer xdim = 0;
    model_block->set_qp_xstart(xdim);
    model_block->set_qp_ystart(ydim);
    xdim += model_block->xdim();
    ydim += model_block->ydim();
    model_block->starting_x(x);

    //initialize dual starting point and constraints
    y.init(ydim, 1, 0.);
    A.init(ydim, xdim, 0.);
    b.init(ydim, 1, 0.);
    genmult(Q, x, Qx);   //Qx=Q*x
    model_block->starting_y(y, Qx, c);
    model_block->get_Ab(A, b);

    //std::cout<<"QP_start: -x^TQx/2+c^Tx+d="<<ip(c,x)-ip(x,Qx)/2+offset<<std::endl;

    //start solving
    mu = 1e100;
    iter = 0;

    status = iterate();

    return status;
  }

  // *************************************************************************
  //                             resolve
  // *************************************************************************

  //this routine is called if the right hand sides of the primal constraints 
  // were changed
  // resolve for the same cost function
  // call starting_point and loop till convergence to optimal solution

  int UQPSolver::resolve() {
    run_starting_point = true;

    //initialize dimensions and primal starting point
    x.init(Q.rowdim(), 1, 0.);
    Integer ydim = 0;
    Integer xdim = 0;
    model_block->set_qp_xstart(xdim);
    model_block->set_qp_ystart(ydim);
    xdim += model_block->xdim();
    ydim += model_block->ydim();
    model_block->starting_x(x);

    //initialize dual starting point and constraints
    y.init(ydim, 1, 0.);
    A.init(ydim, xdim, 0.);
    b.init(ydim, 1, 0.);
    genmult(Q, x, Qx);   //Qx=Q*x
    model_block->starting_y(y, Qx, c);
    model_block->get_Ab(A, b);

    //start solving
    mu = 1e100;
    iter = 0;

    status = iterate();


    return status;
  }

  // *************************************************************************
  //                             update
  // *************************************************************************

  // call starting_point and loop till convergence to optimal solution

  int UQPSolver::update(const Symmatrix& dQ, const Matrix& dc, Real doffset) {
    run_starting_point = false;

    Q += dQ;
    offset += doffset;
    c += dc;
    model_block->restart_x(x, c, dc);
    genmult(Q, x, Qx);   //Qx=Q*x
    model_block->restart_y(y, Qx, c, dc);

    mu = 1e100;
    iter = 0;

    status = iterate();

    return status;
  }

  // *************************************************************************
  //                             save
  // *************************************************************************

  // write all variables to out in order to enable resuming interrupted sessions

  std::ostream& UQPSolver::save(std::ostream& o) const {
    o.precision(20);
    o << Q << c << offset << "\n";
    o << primalval << "\n";
    o << dualval << "\n";
    o << termeps << "\n";
    o << maxiter << "\n";
    o << lowerbound << "\n";
    o << upperbound << "\n";
    o << iter << "\n";
    o << status << "\n";
    return o;
  }

  // *************************************************************************
  //                             restore
  // *************************************************************************

  // read all variables from "in" in order to enable resuming interrupted sessions

  std::istream& UQPSolver::restore(std::istream& in) {
    in >> Q >> c >> offset;
    in >> primalval;
    in >> dualval;
    in >> termeps;
    in >> maxiter;
    in >> lowerbound;
    in >> upperbound;
    in >> iter;
    in >> status;
    return in;
  }

  // *************************************************************************
  //                             print_statistics
  // *************************************************************************

  std::ostream& UQPSolver::print_statistics(std::ostream& out) const {
    out << " qpit " << sum_iter;
    out << " qpcit " << sum_choliter << " qpctime " << sum_choltime;
    out << " QPScoeff " << QPcoeff_time;
    out << " QPSsolve " << QPsolve_time << "\n";

    return out;
  }

  // *************************************************************************
  //                             QPsolve
  // *************************************************************************

  int UQPSolver::QPsolve(const Matrix& center_y,
    Real lower_bound,
    Real upper_bound,
    Real relprec,
    QPSolverProxObject* inHp,
    const MinorantPointer& in_gs_aggr,
    Indexmatrix* yfixed) {
    assert(get_model_data_ptr() != 0);

    Hp = dynamic_cast<BundleProxObject*>(inHp);
    if (Hp == 0) {
      if (cb_out())
        get_out() << "**** ERROR in UQPSolver::QPsolve(.......): dynamic_cast of QPSolverProxObject* to BundleProxObject* failed, cannot continue" << std::endl;
      return 1;
    }
    center_yp = &center_y;
    gs_aggr = in_gs_aggr;

    set_maxiter(100);

    //--- get cost matrices  
    CH_Tools::Microseconds coeff_start = clock.time();

    if (Hp->compute_QP_costs(Q, c, offset,
      get_model_data_ptr()->get_constant_minorant(),
      get_model_data_ptr()->get_bundle(),
      center_y,
      gs_aggr,
      yfixed)) {
      QPcoeff_time += clock.time() - coeff_start;
      if (cb_out())
        get_out() << "**** ERROR UQPSolver::QPSolve(): Hp->compute_QP_costs failed" << std::endl;
      return 1;
    }

    CH_Tools::Microseconds solve_start = clock.time();
    QPcoeff_time += solve_start - coeff_start;

    //--- call the proper quadratic program solver

    set_termbounds(lower_bound, upper_bound);
    set_termeps(relprec);

    int status = solve(Q, c, offset);
    QPsolve_time += clock.time() - solve_start;
    if (status) {
      if (cb_out()) get_out() << "**** WARNING: UQPSolver::QPSolve(): solve() failed and returned " << status << std::endl;
    }


    // if (cb_out(2)){
    //   get_out().precision(16);
    //   Matrix newy(center_y.dim(),1); chk_set_init(newy,1);
    //   Real dummy;
    //   //qp_bundle[0].get_minorant(dummy,newy,0,1,false);
    //   //get_out()<<"\noffset0="<<dummy<<" vec0="<<newy;
    //   MinorantPointer mp;
    //   assert(get_x().dim()==Integer(get_model_data_ptr()->get_bundle().size()));
    //   for (unsigned i=1;i<get_model_data_ptr()->get_bundle().size();i++)
    //     get_model_data_ptr()->get_bundle()[i].get_minorant(mp,get_x()(Integer(i)));
    //   //mp.get_minorant(dummy,newy,0,1./center_y.dim(),false);
    //   //get_out()<<"\noffset115="<<dummy<<" vec115="<<newy;
    //   get_model_data_ptr()->get_bundle()[0].get_minorant(mp,get_x()(Integer(0)));
    //   //mp.get_minorant(dummy,newy,0,1.,false);
    //   //get_out()<<"\noffsetall="<<dummy<<" vecall="<<newy;
    //   get_model_data_ptr()->get_constant_minorant().get_minorant(mp);
    //   //mp.get_minorant(dummy,newy,0,1.,false);
    //   //get_out()<<"\noffset+const="<<dummy<<" vecall+const="<<newy;
    //   MinorantPointer modelmp=mp;
    //   gs_aggr.get_minorant(mp);
    //   mp.get_minorant(dummy,newy,0,-1,false);
    //   //get_out()<<"\noffset+gs="<<-dummy<<" vecall+gs="<<-newy;
    //   Hp->apply_Hinv(newy);
    //   Real normsubg2=-mp.ip(newy);
    //   newy+=center_y;
    //   Real modelval=modelmp.evaluate(-1,newy);
    //   Real linval=mp.evaluate(-1,newy);
    //   Real gs_val=gs_aggr.evaluate(-1,newy);
    //   //Real cutval;
    //   //sbm_transform()->eval_model(cutval,-1,newy,0.1*(upper_bound-linval)/(fabs(upper_bound)+1.));
    //   //cutval+=gs_val;
    //   //Real aggrval=old_model_aggregate.evaluate(-1,newy);
    //   //aggrval+=gs_val;
    //   get_out()<<" primalval="<<-.5*ip(get_x(),Q*get_x())+ip(get_x(),c)+offset;
    //   get_out()<<" nsg2="<<normsubg2<<" modelval="<<modelval<<" gsval="<<gs_val;
    //   get_out()<<" linval="<<linval<<" augval="<<linval+normsubg2/2.<<std::endl;
    //   //get_out()<<" cutval="<<cutval<<" aggrval="<<aggrval<<" aggraug="<<aggrval+normsubg2/2.<<std::endl;
    //   //mp.get_minorant(dummy,newy,0,1.,false);
    //   //get_out()<<"\n mp_offset="<<dummy<<" mp_vec="<<newy;
    //   get_out()<<"x="<<transpose(get_x());
    //   get_out()<<"newy="<<transpose(newy);

    // }


    //qp_mfile_data(center_y,Hp,gs_subg,gs_subg_offset,Q,c,offset,yfixed);

    return status;
  }


  // *************************************************************************
  //                             QPupdate
  // *************************************************************************

  int UQPSolver::QPupdate(const Matrix& center_y,
    Real lower_bound,
    Real upper_bound,
    Real relprec,
    QPSolverProxObject* inHp,
    const MinorantPointer& in_gs_aggr,
    Indexmatrix* yfixed,
    const MinorantPointer& delta_gs_subg,
    const Indexmatrix& delta_index) {
    assert(get_model_data_ptr() != 0);
    assert(Hp);

    BundleProxObject* Hp2 = dynamic_cast<BundleProxObject*>(inHp);
    if (Hp2 != Hp) {
      if (cb_out())
        get_out() << "**** ERROR in QPSolver::QPupdate(.........): dynamic_cast of QPSolverProxObject* to BundleProxObject* produced a different ProxObject than in QPsolve(), cannot continue" << std::endl;
      return 1;
    }
    if (center_yp != &center_y) {
      get_out() << "**** WARNING in QPSolver::QPupdate(.........): center_y object differs form the one in QPsolve() but should be the same" << std::endl;
    }
    gs_aggr = in_gs_aggr;

    set_maxiter(100);

    CH_Tools::Microseconds coeff_start = clock.time();

    ///Indexmatrix yfixedold(yfixed);

    Symmatrix dQ;
    Matrix dc;
    Real doffset;

    if (Hp->update_QP_costs(dQ, dc, doffset,
      get_model_data_ptr()->get_constant_minorant(),
      get_model_data_ptr()->get_bundle(),
      center_y,
      gs_aggr,
      delta_gs_subg, delta_index, yfixed)) {
      if (cb_out())
        get_out() << "**** ERROR UQPSolver::QPupdate(): Hp->update_QP_costs failed" << std::endl;
      QPcoeff_time += clock.time() - coeff_start;
      return 1;
    }

    //TEST
    /*
    Symmatrix Q;
    Matrix c;
    Real offset;

    if (Hp->compute_QP_costs(Q,c,offset,qp_constant_minorant,qp_bundle,center_y,gs_subg,yfixed)){
      if (cb_out())
        get_out()<<"**** ERROR SumBlockModel::reeval_augmodel(): Hp->update_QP_costs failed"<<std::endl;
      return 1;
    }

    Integer xdim(Integer(qp_bundle.size()));
    Matrix t_A(center_y.dim(),xdim,0.);
    Matrix t_c(xdim,1,0.);
    Matrix t_b(center_y.dim(),1,0.);
    Real t_d=0.;
    qp_constant_minorant.get_minorant(t_d,t_b,0);
    gs_subg.get_minorant(t_d,t_b,0,1.,true);
    for (Integer i=0;i<xdim;i++){
      qp_bundle[unsigned(i)].get_minorant(t_c(i),t_A,i);
    }
    const BundleDiagonalTrustRegionScaling* dHp=dynamic_cast<BundleDiagonalTrustRegionScaling*>(Hp);
    if (dHp) {
      Matrix diagH(dHp->get_D());
      Symmatrix tQ(xdim,0.);
      Matrix tc=t_c+transpose(t_A)*center_y;
      Real td=t_d+ip(t_b,center_y);
      for(Integer i=0;i<center_y.dim();i++){
        if (yfixed(i)!=0)
    continue;
        Real di=diagH(i);
        Real bi=t_b(i)/std::sqrt(di);
        Matrix tmpvec(t_A.row(i),1./std::sqrt(di));
        tmpvec.transpose();
        rankadd(tmpvec,tQ,1.,1.);
        tc.xpeya(tmpvec,-bi);
        td-=bi*bi/2.;
      }
      std::cout<<"\n direct: norm2(Qdiff)"<<norm2(Q-tQ);
      std::cout<<" norm2(cdiff)"<<norm2(c-tc);
      std::cout<<" norm2(odiff)"<<std::fabs(offset-td);
      std::cout<<"\n update: norm2(Qdiff)"<<norm2(tQ-solver.get_Q()-dQ);
      std::cout<<" norm2(cdiff)"<<norm2(tc-solver.get_c()-dc);
      std::cout<<" norm2(odiff)"<<std::fabs(td-solver.get_offset()-doffset);
      std::cout<<std::endl;
    }

    std::cout<<" norm2(Qdiff)"<<norm2(Q-solver.get_Q()-dQ);
    std::cout<<" norm2(cdiff)"<<norm2(c-solver.get_c()-dc);
    std::cout<<" norm2(odiff)"<<std::fabs(offset-solver.get_offset()-doffset);
    std::cout<<std::endl;
    if (norm2(c-solver.get_c()-dc)>1e-8*(1+norm2(c))){
      Matrix tmp(c-solver.get_c()-dc);
      std::cout<<"errmat="<<tmp;
      tmp.init(center_y.dim(),1,0.);
      Real dummy;
      delta_gs_subg.get_minorant(dummy,tmp,0);
      tmp.enlarge_right(3,0.);
      for (Integer j=0;j<delta_index.dim();j++){
        Integer ind=delta_index(j);
        tmp(ind,1)=-1.;
        tmp(ind,2)=yfixedold(ind);
        tmp(ind,3)=yfixed(ind);
      }
      std::cout<<"delta="<<tmp;
    }
    */

    CH_Tools::Microseconds solve_start = clock.time();
    QPcoeff_time += solve_start - coeff_start;

    //--- call the quadratic program solver

    set_termbounds(lower_bound, upper_bound);
    set_termeps(relprec);
    int status = update(dQ, dc, doffset);
    //status=solver.solve(Q,c,offset);
    QPsolve_time += clock.time() - solve_start;
    if (status) {
      if (cb_out()) get_out() << "UQPSolver::QPupdate(): update() failed ..." << std::endl;
    }

    /*
    if (cb_out(2)){
      get_out().precision(16);
      Matrix newy(center_y.dim(),1); chk_set_init(newy,1);
      Real dummy;
      //qp_bundle[0].get_minorant(dummy,newy,0,1,false);
      //get_out()<<"\noffset0="<<dummy<<" vec0="<<newy;
      MinorantPointer mp;
      assert(solver.get_x().dim()==Integer(qp_bundle.size()));
      for (unsigned i=1;i<qp_bundle.size();i++)
        qp_bundle[i].get_minorant(mp,solver.get_x()(Integer(i)));
      //mp.get_minorant(dummy,newy,0,1./center_y.dim(),false);
      //get_out()<<"\noffset115="<<dummy<<" vec115="<<newy;
      qp_bundle[0].get_minorant(mp,solver.get_x()(Integer(0)));
      //mp.get_minorant(dummy,newy,0,1.,false);
      //get_out()<<"\noffsetall="<<dummy<<" vecall="<<newy;
      qp_constant_minorant.get_minorant(mp);
      //mp.get_minorant(dummy,newy,0,1.,false);
      //get_out()<<"\noffset+const="<<dummy<<" vecall+const="<<newy;
      gs_subg.get_minorant(mp);
      mp.get_minorant(dummy,newy,0,-1,false);
      //get_out()<<"\noffset+gs="<<-dummy<<" vecall+gs="<<-newy;
      Hp->apply_Hinv(newy);
      Real normsubg2=-mp.ip(newy);
      newy+=center_y;
      Real linval=mp.evaluate(-1,newy);
      Real gs_val=gs_subg.evaluate(-1,newy);
      Real cutval;
      sbm_transform()->eval_model(cutval,-1,newy,0.1*(center_ub+center_gs_val-linval)/(fabs(center_ub+center_gs_val)+1.));
      cutval+=gs_val;
      Real aggrval=old_model_aggregate.evaluate(-1,newy);
      aggrval+=gs_val;
      get_out()<<" primalval="<<-.5*ip(solver.get_x(),Q*solver.get_x())+ip(solver.get_x(),c)+offset;
      get_out()<<" linval="<<linval<<" cutval="<<cutval<<" aggrval="<<aggrval<<" augval="<<linval+normsubg2/2.<<" aggraug="<<aggrval+normsubg2/2.<<std::endl;
      //mp.get_minorant(dummy,newy,0,1.,false);
      //get_out()<<"\n mp_offset="<<dummy<<" mp_vec="<<newy;
      //get_out()<<"x="<<solver.get_x();
    }
    */

    //qp_mfile_data(y,Hp);

    return status;
  }


  // *************************************************************************
  //                             QPresolve
  // *************************************************************************


  int UQPSolver::QPresolve(Real lower_bound,
    Real upper_bound,
    Real relprec) {
    set_termbounds(lower_bound, upper_bound);
    set_termeps(relprec);
    CH_Tools::Microseconds solve_start = clock.time();
    int status = resolve();
    QPsolve_time += clock.time() - solve_start;
    return status;
  }

  // *************************************************************************
  //                             QPget_solution
  // *************************************************************************

  int UQPSolver::QPget_solution(Real& augval_lb,
    Real& augval_ub,
    Matrix& newy,
    Real& gsaggr_offset,
    Matrix& gsaggr_gradient) {
    if (Hp == 0)
      return 1;
    augval_lb = get_primalval();
    augval_ub = get_dualval();
    //augval_ub=augval_lb;
    gsaggr_gradient.newsize(center_yp->dim(), 1); chk_set_init(gsaggr_gradient, 1);
    gs_aggr.get_minorant(gsaggr_offset, gsaggr_gradient, 0, 1., false);
    newy.init(gsaggr_gradient, -1.);
    Real dummy = -gsaggr_offset;
    get_model_data_ptr()->get_constant_minorant().get_minorant(dummy, newy, 0, -1., true);
    assert(get_x().dim() == Integer(get_model_data_ptr()->get_bundle().size()));
    for (unsigned i = 0; i < get_model_data_ptr()->get_bundle().size(); i++)
      get_model_data_ptr()->get_bundle()[i].get_minorant(dummy, newy, 0, -get_x()(Integer(i)), true);
    Hp->apply_Hinv(newy);
    newy += *center_yp;
    return 0;
  }




}

