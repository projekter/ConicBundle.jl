/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/QPSolverParameters.cxx
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
#include "QPSolverParameters.hxx"
#include "QPDirectKKTSolver.hxx"

using namespace CH_Matrix_Classes;

namespace ConicBundle {


  QPSolverParameters::QPSolverParameters(CBout* cb,int incr):CBout(cb,incr)
{
  min_objective_relprec=0.1;
  objective_gap_eps=1e-8;
  lower_bound_gap_eps=1e-8;
  upper_bound_gap_eps=1e-8;
  primal_infeasibility_eps=1e-8;
  dual_infeasibility_eps=1e-8;
  lower_bound=min_Real;
  upper_bound=max_Real;
  maxiter=100;

  allow_unconstrained=true;

  KKTsolver=new QPDirectKKTSolver(false,cb,incr);
  use_predictor_corrector=true;
  use_neighborhood=false;
  //nbh_ub=.99;  
  //nbh_lb=.9;  
  nbh_ub=.9; 
  nbh_lb=.6; 
  use_socqp=false;
}




}

