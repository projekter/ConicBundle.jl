/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  mat_main.cxx
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



#include <iostream>
#include <iomanip>
#include "MatrixCBSolver.hxx"

using namespace std;
using namespace ConicBundle;
using namespace CH_Matrix_Classes;

class QuadraticFunction: public MatrixFunctionOracle
{
private:
  Matrix qcenter;
  double step;

  //f(x)=(x(0)-qcenter(0))^2+(x(1)-qcenter(1))^2
  void quadratic(const Matrix& x,double& val,Matrix& subg) const
  {
    Matrix s=x-qcenter;
    val=ip(s,s);
    subg=2.*s;
  }

public:
  QuadraticFunction(double x,double y,double st=0.1):qcenter(2,1,0.),step(st)
  { qcenter(0)=x; qcenter(1)=y; }

  virtual int evaluate( const Matrix& dual, /* argument/Lagrange multipliers */
                        Real relprec,
		        Real&  objective_value,
			vector<Minorant*>& minorants,
			PrimalExtender*&
		      );
 
};

int QuadraticFunction::evaluate( 
			      const Matrix& x, 
                              Real /* relprec */,
		              Real& objective_value,
			      vector<Minorant*>& minorants, 
			      PrimalExtender*&
			      )
  {
    //construct first subgradient and objective value
    Real val;
    PrimalMatrix h;
    quadratic(x,val,h);
    objective_value=val;
    minorants.push_back(new MatrixMinorant(val,h,h.clone_primal_data()));

    //std::cout<<"val0="<<val<<" h0="<<h;
  
    //add a few further subgradient (not needed)
    Matrix dx(2,1,step);
    quadratic(x+dx,val,h);
    val-=ip(h,dx);
    minorants.push_back(new MatrixMinorant(val,h,h.clone_primal_data()));

    //std::cout<<"val1="<<val<<" h1="<<h;

    dx(1)*=-1.; 
    quadratic(x+dx,val,h);
    val-=ip(h,dx);
    minorants.push_back(new MatrixMinorant(val,h,h.clone_primal_data()));

    //std::cout<<"val2="<<val<<" h2="<<h;

    dx(0)*=-1.; 
    quadratic(x+dx,val,h);
    val-=ip(h,dx);
    minorants.push_back(new MatrixMinorant(val,h,h.clone_primal_data()));

    //std::cout<<"val3="<<val<<" h3="<<h;

    dx(1)*=-1.; 
    quadratic(x+dx,val,h);
    val-=ip(h,dx);
    minorants.push_back(new MatrixMinorant(val,h,h.clone_primal_data()));

    //std::cout<<"val4="<<val<<" h4="<<h<<std::endl;

    return 0;
  }
 

int mat_main()
{
  MatrixCBSolver solver(&cout,1);

  QuadraticFunction fun0(100,100,0.);
  QuadraticFunction fun1(100,98,0.1);

  solver.init_problem(2);

  solver.add_function(fun0); 
  solver.add_function(fun1);

  solver.set_term_relprec(1e-9);           //relative precision for termination

  BundleParameters bp;
  bp.set_max_model_size(3);
  bp.set_max_bundle_size(10);
  solver.set_bundle_parameters(bp,&fun0);  

  bp.set_max_model_size(5);            
  bp.set_max_bundle_size(10);            
  solver.set_bundle_parameters(bp,&fun1);


  //solver.set_next_weight(5.);            //set the initial weight
  
  int cnt=0;

  Matrix center(2,1,0.);
  center(0)=17.;
  center(1)=23.;
  if (solver.set_new_center_point(center)){
    cout<<"**** ERROR: main(): solver.set_new_center_point() failed"<<endl;
  }

  do {

    int retval;

    
    /* call solve with option to stop after the next descent step */
    if ((retval=solver.solve(0,true))){
      cout<<"solve() returned "<<retval<<endl;
      return 1;
    }

    
    /* get solution information */
    double obj=solver.get_objval();
    Matrix y;
    if ((retval=solver.get_center(y))){
      cout<<"get_center() returned "<<retval<<endl;
      return 1;
    }
    double u=solver.get_last_weight();

    cout<<cnt<<": "<<obj<<" ["<<y(0)<<","<<y(1)<<"] "<<u<<endl; 


    /* get primal information */
    const PrimalMatrix* x0=dynamic_cast<const PrimalMatrix*>(solver.get_approximate_primal(fun0));
    if (x0==0){
      cout<<"get_approximate_primal() failed for function 0"<<endl;
      return 1;
    }

    cout<<" xf0: "<<" ["<<(*x0)(0)<<","<<(*x0)(1)<<"] "<<endl;

    const PrimalMatrix* x1=dynamic_cast<const PrimalMatrix*>(solver.get_approximate_primal(fun1));
    if (x1==0){
      cout<<"get_approximate_primal() failed for function 1"<<endl;
      return 1;
    }

    cout<<" xf1: "<<" ["<<(*x1)(0)<<","<<(*x1)(1)<<"] "<<endl;


    /* set some paramters */
    solver.set_min_weight(0.01);
    solver.set_max_weight(100.);

    cnt++;

  } while (!solver.termination_code());

  solver.print_termination_code(cout);

  return 0;
}

int main()
{
  return mat_main();
}
