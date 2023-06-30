/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/CFunction.hxx
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



#ifndef CONICBUNDLE_CFUNCTION_HXX
#define CONICBUNDLE_CFUNCTION_HXX

/**  @file CFunction.hxx
    @brief Header declaring the classes ConicBundle::CFunction and ConicBundle::CFunctionMinorantExtender
    @version 1.0
    @date 2014-07-23
    @author Christoph Helmberg
*/

//------------------------------------------------------------

#include "MatrixCBSolver.hxx"
#include "cb_cinterface.h"

//------------------------------------------------------------


namespace ConicBundle {
/** @ingroup internal_cinterface
*/
//@{
  
/**@brief for the "C" interface this maps c oracles to the standard function oracle with matrix classes
       
*/	

  class CFunction: public CBout, public MatrixFunctionOracle 
  {
  private:
    void* function_key;  ///< identifier for c-code
    cb_functionp oracle; ///< c-function for evaluate
    cb_subgextp subgext; ///< c-function for subgradient extension
    CH_Matrix_Classes::Integer primaldim;   ///< length of primal vectors; uses PrimalMatrix
    CH_Matrix_Classes::Integer max_new;     ///< maximum number of new vectors per call

  public:
    ///constructor
    CFunction(void* fk,cb_functionp fp,cb_subgextp se=0,int prdim=0);
    ///destructor
    ~CFunction(){};

    ///set the maximum number of new subgardients per evaluations
    void set_max_new(CH_Matrix_Classes::Integer mn)
    {max_new=CH_Matrix_Classes::max(CH_Matrix_Classes::Integer(1),mn);}

    /// see MatrixFunctionOracle::evaluate() for explanations 
    int evaluate(
		 const  CH_Matrix_Classes::Matrix& current_point,
		 double relprec,
		 double&  objective_value,
		 std::vector<Minorant*>& minorants,
		 PrimalExtender*& 
     );

    /// see MatrixFunctionOracle::apply_modfication() for explanations 
    int apply_modification
    (
     const OracleModification& oracle_modification ,
     const CH_Matrix_Classes::Matrix* new_center ,
     const CH_Matrix_Classes::Matrix* old_center ,
     bool& discard_objective_in_center ,
     bool& discard_model , 
     bool& discard_aggregates ,
     MinorantExtender*& minorant_extender 
     );

       
  };

  //@}
  
/** @ingroup internal_cinterface
*/
//@{
  
  /**@brief MinorantExtender for CFunction

   */

  class CFunctionMinorantExtender: public CBout, public MinorantExtender
  {
  private:
    void* function_key;  ///< identifier for c-code
    cb_subgextp subgext; ///< c-function for subgradient extension

    public:
    /// constructor
    CFunctionMinorantExtender(void* fk,cb_subgextp se):function_key(fk),subgext(se)
    { assert(se); }

    /// destructor
    ~CFunctionMinorantExtender(){}

    /// see MinorantExtender::extend() for explanations
    int extend(Minorant& minorant,int n_coords,const int* indices);
  };

  //@}

}
#endif

