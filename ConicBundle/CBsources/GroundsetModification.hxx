/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/GroundsetModification.hxx
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




#ifndef CONICBUNDLE_GROUNDSETMODIFICATION_HXX
#define CONICBUNDLE_GROUNDSETMODIFICATION_HXX


/**  @file GroundsetModification.hxx
    @brief Header declaring the class ConicBundle::GroundsetModification
    @version 1.0
    @date 2014-11-10
    @author Christoph Helmberg
*/

#include "Modification.hxx"

namespace ConicBundle {

  class Minorant;

  /** @ingroup dynamic_modification_support
  */
  //@{

  /** @brief Collects modifications for the unconstrained Groundset for appending, deleting or reassigning variables

      This class is implemented by restricting the possibilities of ConicBundle::Modification to the operations involving variables, so all routines simply call those of Modification. See the description there for usage and default values.

   */


  class GroundsetModification : public CBout, public OracleModification {
  protected:
    /// this class provides a restricted interface to this Modification instance where all modifications are organized and stored
    Modification mdf;
    /// specifies which value should be added to the offset of the groundset objective
    CH_Matrix_Classes::Real offset;

  public:
    //-----------------------------------------------------------------------
    /** @name Constructors and initialization
     */
     //@{

     /// constructor, calls modification constructor
    GroundsetModification(CH_Matrix_Classes::Integer var_olddim = 0, const CBout* cbo = 0, int incr = -1) :
      CBout(cbo, incr), mdf(var_olddim, 0), offset(0.) {
    }

    /// reset modifications to an unmodified object currently having var_olddim variables, calls Modification::clear
    void clear(CH_Matrix_Classes::Integer var_olddim) {
      mdf.clear(var_olddim, 0); offset = 0.;
    }


    virtual ~GroundsetModification();


    //@}

    //-----------------------------------------------------------------------
    /** @name  Routines for adding modifications
     */
     //@{

     /// append @a append_dim new variables with start_val as initial values (if NULL, use default value), calls Modification::add_append_vars
    int add_append_vars(CH_Matrix_Classes::Integer append_dim,
      const CH_Matrix_Classes::Matrix* start_val = 0,
      const CH_Matrix_Classes::Matrix* costs = 0) {
      return mdf.add_append_vars(append_dim, 0, 0, 0, start_val, costs);
    }


    /// reassign the variables as given in @a map_to_old, calls Modification::add_reassign_vars
    int add_reassign_vars(const CH_Matrix_Classes::Indexmatrix& map_to_old) {
      return mdf.add_reassign_vars(map_to_old);
    }

    /// delete the variables indexed by @a del_ind; for each new index @a map_to_old returns the old one; calls Modification::add_delete_vars
    int add_delete_vars(const CH_Matrix_Classes::Indexmatrix& del_ind,
      CH_Matrix_Classes::Indexmatrix& map_to_old) {
      return mdf.add_delete_vars(del_ind, map_to_old);
    }

    ///add to the current offset the value @a delta
    int add_offset(CH_Matrix_Classes::Real delta) {
      offset += delta; return 0;
    }

    //@}

    //-----------------------------------------------------------------------
    /** @name  Routines for applying the collected modifications
     */

     //@{

     /// carry out the collected modifications on the given vector by calling Modification::apply_to_vars
    int apply_to_vars(CH_Matrix_Classes::Matrix& vars) const {
      return mdf.apply_to_vars(&vars, 0, 0, 0);
    }

    /// carry out the collected modifications on the given vector by calling Modification::apply_to_vars
    int apply_to_costs(CH_Matrix_Classes::Matrix& costs, CH_Matrix_Classes::Real& in_offset) const {
      in_offset += offset; return mdf.apply_to_vars(0, 0, 0, &costs);
    }

    //@}

    //-----------------------------------------------------------------------
    /** @name  Routines for querying properties of the collected modifications
     */

     //@{


     /// returns true if no modifications need to be executed except possibly an offset change for the ground set minorant
    bool no_modification() const {
      return mdf.no_modification();
    }

    /// if set to true, no deletions/reassignments may be present or specified in the future, only appensions are allowed 
    int set_append_to_old(bool append_only) {
      return mdf.set_append_to_old(append_only);
    }
    /// returns true if this only contains appending operations and incorporating this is done with respect to the old dimension 
    bool append_to_old() const {
      return mdf.append_to_old();
    }

    /// returns true if no variables were added or deleted (allows permutations), false otherwise
    bool no_additions_or_deletions_in_vars() const {
      return ((mdf.appended_vardim() == 0) && ((mdf.map_to_old_variables() == 0) || (mdf.map_to_old_variables()->dim() == mdf.old_vardim())));
    }

    /// returns true if all entries deleted in @a oldpoint (must be a vector of length old_vardim()) are 0,  false otherwise
    bool deleted_variables_are_zero(const CH_Matrix_Classes::Matrix& oldpoint) const {
      return mdf.deleted_variables_are_zero(oldpoint);
    }

    /// returns true if all entries in newpoint (must be a vector of length new_vardim()) that correspond to new variables have value 0 and false otherwise
    bool new_variables_are_zero(const CH_Matrix_Classes::Matrix& newpoint) const {
      return mdf.new_variables_are_zero(newpoint);
    }

    /// returns true if the values in newpoint (must be a vector of length new_vardim()) that correspond to old variables match the old values stored in oldpoint (must be a vector of length old_vardim()) and false otherwise
    bool mapped_variables_are_equal(const CH_Matrix_Classes::Matrix& newpoint,
      const CH_Matrix_Classes::Matrix& oldpoint) const {
      return mdf.mapped_variables_are_equal(newpoint, oldpoint);
    }

    //@}

    //-----------------------------------------------------------------------
    /** @name  Routines implementing the abstract OracleModifcation messages
     */

     //@{

     /// returns the number of variables before modification
    int get_old_vardim() const {
      return mdf.old_vardim();
    }
    /// returns the number of variables once all stored modifications have been performed
    int get_new_vardim() const {
      return mdf.new_vardim();
    }
    /// returns the number of variables that are appended (due to later reassignmentds they may no longer be located at the end)
    int get_appended_vardim() const {
      return mdf.appended_vardim();
    }

    /// returns null if there are no index changes, otherwise the Indexmatrix pointed to is a vector whose i-th entry holds the old index of the new i-th variable (injective!), index values exceeding old_vardim() refer to newly appended variables  
    const int* get_map_to_old_variables() const {
      return (mdf.map_to_old_variables()) == 0 ? 0 : (const int*)(mdf.map_to_old_variables()->get_store());
    }

    /// returns the change in the offste value of the groundset minorant
    CH_Matrix_Classes::Real get_add_offset() const {
      return offset;
    }

    /// returns the change in the offste value of the groundset minorant
    const CH_Matrix_Classes::Matrix* get_append_costs() const {
      return mdf.get_var_append_costs();
    }

    /// incorporate the OracleModification @a m (it should only contain variable changes, but this is not checked!) into this one; calls Modification::incorporate
    int incorporate(const OracleModification& m) {
      const GroundsetModification* gm = dynamic_cast<const GroundsetModification*>(&m);
      if (gm) {
        offset += gm->offset; return mdf.incorporate(gm->mdf);
      }
      return 1;
    }

    /// returns a new object on the heap, that allows to incorporate this but starts off from a function whose input argument dimension is old_var_dim 
    OracleModification* new_initial_oraclemodification(int old_var_dim) const {
      return new GroundsetModification(old_var_dim);
    }

    /// append append_dim further variables at the end of the argument vector (without specifying their effect, they may e.g. be ignored by the function)
    int add_append_variables(int append_dim) {
      return add_append_vars(append_dim);
    }

    /// reorder and resize the variables as given by the first new_dim entries of map_to_old_indices; each former index may appear at most once in this list, so new_dim < get_new_vardim()
    int add_reassign_variables(int new_dim, const int* map_to_old_indices) {
      CH_Matrix_Classes::Indexmatrix map_to_old(new_dim, 1, map_to_old_indices);
      return add_reassign_vars(map_to_old);
    }

    //@}

    //-----------------------------------------------------------------------
    /** @name  Routines for retrieving the detailed collected modifications
     */

     //@{

     /// returns the number of variables before modification (given on initialization)
    CH_Matrix_Classes::Integer old_vardim() const {
      return mdf.old_vardim();
    }
    /// returns the number of variables once all stored modifications have been performed
    CH_Matrix_Classes::Integer new_vardim() const {
      return mdf.new_vardim();
    }
    /// returns the number of variables that are appended (due to later reassignmentds they may no longer be located at the end)
    CH_Matrix_Classes::Integer appended_vardim() const {
      return mdf.appended_vardim();
    }

    /// returns null if there are no index changes, otherwise the Indexmatrix pointed to is a vector whose i-th entry holds the old index of the new i-th variable (injective!), index values exceeding old_vardim() refer to newly appended variables  
    const CH_Matrix_Classes::Indexmatrix* map_to_old_variables() const {
      return mdf.map_to_old_variables();
    }

    /// returns null if there were no deletions, otherwise the Indexmatrix pointed to is a vector holding the deleted old variable indices in increasing order   
    const CH_Matrix_Classes::Indexmatrix* deleted_var_indices() const {
      return mdf.deleted_var_indices();
    }

    /// returns null if no variables were added, otherwise the Indexmatrix pointed ato is a vector holding the new indices of the new variables in increasing order
    const CH_Matrix_Classes::Indexmatrix* new_var_indices() const {
      return mdf.new_var_indices();
    }

    //@}

    //-----------------------------------------------------------------------
    /** @name  Output settings
     */

     //@{

     /// see CBout::set_out
    void set_out(std::ostream* out = 0, int print_level = 1) {
      CBout::set_out(out, print_level); mdf.set_out(out, print_level);
    }

    /// see CBout::set_out
    void set_cbout(const CBout* cb = 0, int incr = -1) {
      CBout::set_cbout(cb, incr); mdf.set_cbout(this, 0);
    }

    //@}

  };



  //@}

}

#endif

