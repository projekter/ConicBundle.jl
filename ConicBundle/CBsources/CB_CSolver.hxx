/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/CB_CSolver.hxx
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



#ifndef CONICBUNDLE_CB_CSOLVER_HXX
#define CONICBUNDLE_CB_CSOLVER_HXX

/**  @file CB_CSolver.hxx
    @brief Header declaring the class CB_CSolver
    @version 1.0
    @date 2014-07-23
    @author Christoph Helmberg
*/

#include <map>
#include "cb_cinterface.h"
#include "MatrixCBSolver.hxx"
#include "CFunction.hxx"

/** @defgroup internal_cinterface Internal implementation of the "C" interface
 
   @brief Internally the "C" interface is implemented as follows. 
   cb_construct_problem() generates a CB_CSolver that contains a 
   MatrixCBSolver and feeds each c-evaluation function
   to this solver in the form of a CFunction. All calls are
   then passed on to this solver.
*/

//@{

/** @brief Interface class for implementing the language "C" interface

   This interface is provided to maintain compatibility with former
   versions. The option no_bundle is no longer available, because
   it is not compatabile with the use of penalty modes. It is 
   now replaced by the sumbundle feature.

 */

class CB_CSolver
{
private: 
public:
  bool no_bundle;  ///< if true, swith on the minimal sumbundle version
  std::map<void*, ConicBundle::CFunction*> funmap; ///< maps to the c functions
  ConicBundle::MatrixCBSolver* solver; ///< the actual solver

  /// constructor;
  CB_CSolver(bool no_bundle);

  /// destructor
  ~CB_CSolver();
	  
};

//@}

#endif

