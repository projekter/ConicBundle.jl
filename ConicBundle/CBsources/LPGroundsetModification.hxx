/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/LPGroundsetModification.hxx
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




#ifndef CONICBUNDLE_LPGROUNDSETMODIFICATION_HXX
#define CONICBUNDLE_LPGROUNDSETMODIFICATION_HXX


/**  @file LPGroundsetModification.hxx
    @brief Header declaring the class ConicBundle::LPGroundsetModification
    @version 1.0
    @date 2014-11-10
    @author Christoph Helmberg
*/

#include "GroundsetModification.hxx"

namespace ConicBundle {

  /** @ingroup dynamic_modification_support
  */
  //@{

  /** @brief Collects modifications for the linearly constrained LPGroundset for appending, deleting or reassigning variables (with bounds and starting values) and constraints (with bounds)

      This class is implemented by restricting the possibilities of ConicBundle::Modification to the operations involving variables, so all routines simply call those of Modification. See the description there for usage and default values.

      The class is derived from GroundsetModification which only describes changes of variables. If this class is viewed as a GroundsetModification, it also only reveals the effects on the variables.

   */


  class LPGroundsetModification : public GroundsetModification {
  protected:
    //Modification mdf;

  public:
    //-----------------------------------------------------------------------
    /** @name Constructors and initialization
     */
     //@{

     /// reset modifications to an unmodified object currently having var_olddim variables and row_olddim constraints,  calls Modification::clear
    void clear(CH_Matrix_Classes::Integer var_olddim, CH_Matrix_Classes::Integer row_olddim) {
      mdf.clear(var_olddim, row_olddim); offset = 0.;
    }

    /// calls clear(var_olddim,row_olddim)
    LPGroundsetModification(CH_Matrix_Classes::Integer var_olddim, CH_Matrix_Classes::Integer row_olddim, const CBout* cbo = 0, int incr = -1) :GroundsetModification(0, cbo, incr) {
      clear(var_olddim, row_olddim);
    }

    virtual ~LPGroundsetModification();


    //@}

    //-----------------------------------------------------------------------
    /** @name  Routines for adding modifications
     */
     //@{


     /** @brief set the lower bound on variable with index @a ind to value @a lb

          If value @a lb exceeds CB_plus_infinity or the index is out of
          range, the return value is 1 and no changes are carried out,
          otherwise it returns 0.
     */
    int add_set_lb(CH_Matrix_Classes::Integer ind, CH_Matrix_Classes::Real lb) {
      return mdf.add_set_lb(ind, lb);
    }

    /** @brief set the upper bound on variable with index @a ind to value @a ub

         If value @a ub is below CB_minus_infinity or the index is out of
         range, the return value is 1 and no changes are carried out,
         otherwise it returns 0.
    */
    int add_set_ub(CH_Matrix_Classes::Integer ind, CH_Matrix_Classes::Real ub) {
      return mdf.add_set_ub(ind, ub);
    }

    //int add_append_vars(CH_Matrix_Classes::Integer append_dim,
    //		      const CH_Matrix_Classes::Matrix* start_val=0)

    /// append @a append_dim new variables with lower bounds @a append_lb (NULL: default value) and upper bounds append_ub (NULL: default value), corresponding columns @a append_cols (NULL: default value) and @a start_val as initial value (NULL: default value) by calling Modification::add_append_vars
    int add_append_vars(CH_Matrix_Classes::Integer append_dim,
      const CH_Matrix_Classes::Matrix* append_lb,
      const CH_Matrix_Classes::Matrix* append_ub,
      const CH_Matrix_Classes::Sparsemat* append_cols,
      const CH_Matrix_Classes::Matrix* start_val,
      const CH_Matrix_Classes::Matrix* append_costs) {
      return mdf.add_append_vars(append_dim, append_lb, append_ub, append_cols, start_val, append_costs);
    }

    //int add_reassign_vars(const CH_Matrix_Classes::Indexmatrix& map_to_old);

    //int add_delete_vars(const CH_Matrix_Classes::Indexmatrix& del_ind,
    //		      CH_Matrix_Classes::Indexmatrix& map_to_old);

    /** @brief set the lower bound on row right hand side with index @a ind to value @a rhslb

        If value @a rhslb exceeds CB_plus_infinity or the index is out of
        range, the return value is 1 and no changes are carried out,
        otherwise it returns 0.
    */
    int add_set_rhslb(CH_Matrix_Classes::Integer ind, CH_Matrix_Classes::Real rhslb) {
      return mdf.add_set_rhslb(ind, rhslb);
    }

    /** @brief set the upper bound on row right hand side with index @a ind to value @a rhsub

        If value @a rhsub is below CB_minus_infinity or the index is out of
        range, the return value is 1 and no changes are carried out,
        otherwise it returns 0.
    */
    int add_set_rhsub(CH_Matrix_Classes::Integer ind, CH_Matrix_Classes::Real rhsub) {
      return mdf.add_set_rhslb(ind, rhsub);
    }

    /// append @a append_dim new rows as in append_rows (if NULL, use default value) with lower bound append_rhslb (if NULL use default) and upper bound append_rhsub (if NULL use default), calls Modification::add_append_rows
    int add_append_rows(CH_Matrix_Classes::Integer append_dim,
      const CH_Matrix_Classes::Sparsemat* append_rows,
      const CH_Matrix_Classes::Matrix* append_rhslb,
      const CH_Matrix_Classes::Matrix* append_rhsub) {
      return mdf.add_append_rows(append_dim, append_rows, append_rhslb, append_rhsub);
    }

    //int add_reassign_rows(const CH_Matrix_Classes::Indexmatrix& map_to_old);

    //int add_delete_rows(const CH_Matrix_Classes::Indexmatrix& rows_del_ind,
    //		      CH_Matrix_Classes::Indexmatrix& rows_map_to_old);

    //int incorporate(const GroundsetModification& m);

    /// incorporate the LPGroundsetModification @a m into this one; calls Modification::incorporate
    /// incorporate the OracleModification @a m (it should only contain variable changes, but this is not checked!) into this one; calls Modification::incorporate
    int incorporate(const OracleModification& m) {
      const LPGroundsetModification* lpgm = dynamic_cast<const LPGroundsetModification*>(&m);
      if (lpgm) {
        offset += lpgm->offset; return mdf.incorporate(lpgm->mdf);
      }
      return GroundsetModification::incorporate(m);
    }

    /// returns a new object on the heap, that allows to incorporate this but starts off from a function whose input argument dimension is old_var_dim 
    OracleModification* new_initial_oraclemodification(int old_var_dim) const {
      return new LPGroundsetModification(old_var_dim, 0);
    }

    //@}

    //-----------------------------------------------------------------------
    /** @name  Routines for applying the collected modifications
     */

     //@{

     //int apply_to_vars(CH_Matrix_Classes::Matrix& vars) const

     /// carry out the collected modifications on the given vectors of lower and upper bounds on the variables by calling Modification::apply_to_vars
    int apply_to_bounds(CH_Matrix_Classes::Matrix& lb,
      CH_Matrix_Classes::Matrix& ub) const {
      return mdf.apply_to_vars(0, &lb, &ub, 0);
    }

    //int apply_to_costs(CH_Matrix_Classes::Matrix& costs,CH_Matrix_Classes::Real& offset) const

    /// carry out the collected modification on the constraint matrix and the corresponding lower and upper bounds by calling Modification::apply_to_rows
    int apply_to_rows(CH_Matrix_Classes::Sparsemat& rows,
      CH_Matrix_Classes::Matrix& rhslb,
      CH_Matrix_Classes::Matrix& rhsub) const {
      return mdf.apply_to_rows(&rows, &rhslb, &rhsub);
    }

    //@}

    //-----------------------------------------------------------------------
    /** @name  Routines for querying properties of the collected modifications
     */

     //@{

     //bool no_modification() const;

     //bool no_additions_or_deletions_in_vars() const

     //bool deleted_variables_are_zero(const CH_Matrix_Classes::Matrix& oldpoint) const

     //bool new_variables_are_zero(const CH_Matrix_Classes::Matrix& newpoint) const

     //bool mapped_variables_are_equal(const CH_Matrix_Classes::Matrix& newpoint,
     //				  const CH_Matrix_Classes::Matrix& oldpoint) const

     //@}

     //-----------------------------------------------------------------------
     /** @name  Routines for retrieving the detailed collected modifications
      */

      //@{

      //CH_Matrix_Classes::Integer old_vardim() const;
      //CH_Matrix_Classes::Integer new_vardim() const;
      //CH_Matrix_Classes::Integer appended_vardim() const;

      /// returns the number of rows before modification (given on initialization)
    CH_Matrix_Classes::Integer old_rowdim() const {
      return mdf.old_rowdim();
    }
    /// returns the number of rows once all stored modifications have been performed
    CH_Matrix_Classes::Integer new_rowdim() const {
      return mdf.new_rowdim();
    }
    /// returns the number of rows that are appended (due to later reassignments they may no longer be located at the end)
    CH_Matrix_Classes::Integer appended_rowdim() const {
      return mdf.appended_rowdim();
    }

    //const CH_Matrix_Classes::Indexmatrix* map_to_old_variables() const;

    //const CH_Matrix_Classes::Indexmatrix* deleted_var_indices() const;

    //const CH_Matrix_Classes::Indexmatrix* new_var_indices() const;

    //const CH_Matrix_Classes::Indexmatrix* map_to_old_rows() const;

    //const CH_Matrix_Classes::Indexmatrix* deleted_row_indices() const;

    //const CH_Matrix_Classes::Indexmatrix* new_row_indices() const;

  };



  //@}

}

#endif

