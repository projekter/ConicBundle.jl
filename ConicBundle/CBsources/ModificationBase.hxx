/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/ModificationBase.hxx
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




#ifndef CONICBUNDLE_MODIFICATIONBASE_HXX
#define CONICBUNDLE_MODIFICATIONBASE_HXX


/**  @file Modification.hxx
    @brief Header declaring the class ConicBundle::ModificationBase
    @version 1.0
    @date 2016-12-30
    @author Christoph Helmberg
*/

#include <map>
#include "CBout.hxx"
#include "indexmat.hxx"

namespace ConicBundle {

/** @ingroup dynamic_modification_support
*/
//@{

/** @brief basic routines for the base classes for reorganizing maps 
 */

class ModificationBase: public CBout
{
protected:

  //-----------------------------------------------------------------------
  /** @name Internal routines for updating maps in reassign and delete steps
   */
  //@{

  /** @brief internal routine for changing the maps in either of 
      add_reassign_vars and add_reassign_rows 

     @param[in,out] adapt_this_map
            holds either var_map_to_old or row_map_to_old

     @param[in,out] del_ind
            holds either var_del_ind or row_del_ind

     @param[in,out] new_ind
            holds either var_new_ind or row_new_ind

     @param[out] append_del_ind
            holds the indices to be deleted in the variables/rows
            that up till now should have been appended

     @param[in] input_map
            holds the input reassignment map_to_old for the indices

     @param[in,out] append_dim
            holds either var_append_dim or row_append_dim

     @param[in] olddim
            gives either var_olddim or row_olddim 

     @param[in] newdim
            gives either var_newdim or row_newdim at input time 

     @return number of errors detected in the input map,
            in case of errors nothing is changed at all
  */
  int adapt_map_to_old(CH_Matrix_Classes::Indexmatrix*& adapt_this_map,
		       CH_Matrix_Classes::Indexmatrix*& del_ind,
		       CH_Matrix_Classes::Indexmatrix*& new_ind,
		       CH_Matrix_Classes::Indexmatrix& append_del_ind,
		       const CH_Matrix_Classes::Indexmatrix& input_map,
		       CH_Matrix_Classes::Integer& append_dim,
		       CH_Matrix_Classes::Integer olddim,
		       CH_Matrix_Classes::Integer newdim) const;

  /** @brief internal routine that converts in either of 
      add_delete_vars and add_delete_rows the list of indices to be
      deleted into a map_to_old that may the be treated in
      add_reassign_vars or add_reassign_rows respectively.
      Returns the number of errors detected in the input map.
      In case of errors, the returned map will not be valid.
  */
  int form_map_to_old(CH_Matrix_Classes::Indexmatrix& map_to_old,
		      const CH_Matrix_Classes::Indexmatrix& del_ind,
		      CH_Matrix_Classes::Integer dim) const;

  //@}

public:
  //-----------------------------------------------------------------------
  /** @name Constructors and initialization 
   */
  //@{

  ///
  virtual ~ModificationBase();

  ///
  ModificationBase(const CBout* cbout=0,int incr=0):CBout(cbout,incr){}

  //@}
  

};





  //@}

}

#endif

