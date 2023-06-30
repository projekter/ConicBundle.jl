/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  cxx_main.cxx
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
#include "CBSolver.hxx"

using namespace ConicBundle;
using namespace std;

class QuadraticFunction: public FunctionOracle
{
private:
  double x_center;
  double y_center;
  double step;

  //f(x,y)=(x-x_center)^2+(y-y_center)^2
  void quadratic(double x,double y,double& val,DVector& subg) const
  {
    double s0=x-x_center;
    double s1=y-y_center;
    
    val=s0*s0+s1*s1;
    subg.clear();
    subg.push_back(2.*s0);
    subg.push_back(2.*s1);
  }

public:
  QuadraticFunction(double x,double y,double st=0.1):
    x_center(x),y_center(y),step(st){}

  int evaluate(  const double* dual, /* argument/Lagrange multipliers */
		 double relprec,
		 double&   objective_value,
		 vector<Minorant*>&  minorants,
		 PrimalExtender*&
		 );
 
};

int QuadraticFunction::evaluate( 
			      const double* dual,
                              double /* relprec */,
		              double&   objective_value,
			      vector<Minorant*>& minorants,
			      PrimalExtender*&
			      )
  {
    //construct first subgradient and objective value
    double val;
    PrimalDVector h;
    quadratic(dual[0],dual[1],val,h);
    objective_value=val;
    minorants.push_back(new Minorant(val,h,h.clone_primal_data()));
  
    //add a few further subgradient (not needed)
    quadratic(dual[0]+step,dual[1]+step,val,h);
    val-=h[0]*step+h[1]*step;
    minorants.push_back(new Minorant(val,h,h.clone_primal_data()));

    quadratic(dual[0]+step,dual[1]-step,val,h);
    val-=h[0]*step-h[1]*step;
    minorants.push_back(new Minorant(val,h,h.clone_primal_data()));

    quadratic(dual[0]-step,dual[1]-step,val,h);
    val-=-h[0]*step-h[1]*step;
    minorants.push_back(new Minorant(val,h,h.clone_primal_data()));

    quadratic(dual[0]-step,dual[1]+step,val,h);
    val-=-h[0]*step+h[1]*step;
    minorants.push_back(new Minorant(val,h,h.clone_primal_data()));

    return 0;
  }
 

int cxx_main()
{
  CBSolver solver(&cout,1);

  QuadraticFunction fun0(100,100,0.);
  QuadraticFunction fun1(100,98,0.1);

  solver.init_problem(2);

  solver.add_function(fun0); 
  solver.add_function(fun1);

  solver.set_term_relprec(1e-9);           //relative precision for termination

  BundleParameters bp;
  bp.set_max_model_size(3);
  bp.set_max_bundle_size(10);
  solver.set_bundle_parameters(fun0,bp);  

  bp.set_max_model_size(5);            
  bp.set_max_bundle_size(10);            
  solver.set_bundle_parameters(fun1,bp);

  //solver.set_next_weight(5.);            //set the initial weight
  
  int cnt=0;

  DVector center;
  center.push_back(17.);
  center.push_back(23.);
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
    DVector y;
    if ((retval=solver.get_center(y))){
      cout<<"get_center() returned "<<retval<<endl;
      return 1;
    }
    double u=solver.get_last_weight();

    cout<<cnt<<": "<<obj<<" ["<<y[0]<<","<<y[1]<<"] "<<u<<endl; 


    /* get primal information */
    const PrimalDVector* x0=dynamic_cast<const PrimalDVector*>(solver.get_approximate_primal(fun0));
    if(x0==0){
      cout<<"get_approximate_primal() failed for function 0"<<endl;
      return 1;
    }


    cout<<" xf0: "<<" ["<<(*x0)[0]<<","<<(*x0)[1]<<"] "<<endl;

    const PrimalDVector* x1=dynamic_cast<const PrimalDVector*>(solver.get_approximate_primal(fun1));
    if (x1==0){
      cout<<"get_approximate_primal() failed for function 1"<<endl;
      return 1;
    }

    cout<<" xf1: "<<" ["<<(*x1)[0]<<","<<(*x1)[1]<<"] "<<endl;


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
  return cxx_main();
}
