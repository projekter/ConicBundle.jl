/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/AFTModification.hxx
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




#ifndef CONICBUNDLE_AFTMODIFICATION_HXX
#define CONICBUNDLE_AFTMODIFICATION_HXX


/**  @file AFTModification.hxx
    @brief Header declaring the class ConicBundle::AFTModification
    @version 1.0
    @date 2014-11-10
    @author Christoph Helmberg
*/

#include "Modification.hxx"

namespace ConicBundle {

/** @ingroup dynamic_modification_support
*/
//@{

/** @brief collects modifications for an AffineFunctionTransformation for the scaling and offset constants as well as for appending, deleting or reassigning columns and rows of the transformation data

   The matrix modifications are implemented by restricting the possiblities of ConicBundle::Modification, so most routines simply call those of Modification. See the description there for usage and default values (one exception: for the right hand side offset the default value is 0. here).

   A somewhat peculiar property of AffineFunctionTransformations is that a NULL matrix is interpreted as the identity. If modifications are consistent with this interpretation (number of rows and columns stays the same and NULL is used for the transformation matrix), the NULL matrix is preserved. 
 */


  class AFTModification:public CBout, public OracleModification
{
protected:
  /// this class provides a restricted interface to this Modification instance where all modifications are organized and stored  
  Modification mdf;

  /// multiply the factor by this value
  CH_Matrix_Classes::Real factor;
  /// add this value to the offset
  CH_Matrix_Classes::Real offset;

  /// if *this does not explicitly change the transformation, usually the groundset modifications will be applied in AffineFunctionTransformation::apply_modification(); if this is to be avoided, set ignore_gs_mdf=true.
  bool ignore_gs_mdf;
  
  mutable int preserves_id_flag; ///< for lazy evaulation, -1 if not yet determined, 0 if false, 1  if true
  
public:
  //-----------------------------------------------------------------------
  /** @name Constructors and initialization 
   */
  //@{


  /// initialize and reset to an unmodified object currently having var_olddim columns and row_olddim rows
  AFTModification(CH_Matrix_Classes::Integer var_olddim=0,CH_Matrix_Classes::Integer row_olddim=0,bool ignore_groundset_modification=false);
  virtual ~AFTModification();

  /// reset modifications to an unmodified object currently having var_olddim columns and row_olddim rows,  calls Modification::clear
  void clear(CH_Matrix_Classes::Integer var_olddim,CH_Matrix_Classes::Integer row_olddim);

  //@}
  
  //-----------------------------------------------------------------------
  /** @name  Routines for adding modifications
   */
  //@{

  
  /// append @a append_dim new variables/columns with values specified by @a append_cols (NULL: default value) and @a linear_costs coefficients (NULL: default value) by calling Modification::add_append_vars
  int add_append_vars(CH_Matrix_Classes::Integer append_dim,
		      const CH_Matrix_Classes::Sparsemat* append_cols,
		      const CH_Matrix_Classes::Matrix* linear_costs)
  {preserves_id_flag=-1;return mdf.add_append_vars(append_dim,0,0,append_cols,0,linear_costs);}

  ///reassign the variables/columns as given in @a map_to_old, calls Modification::add_reassign_vars
  int add_reassign_vars(const CH_Matrix_Classes::Indexmatrix& map_to_old)
  {preserves_id_flag=-1;return mdf.add_reassign_vars(map_to_old);}

  /// delete the variables indexed by @a del_ind; for each new index @a map_to_old returns the old one; calls Modification::add_delete_vars
  int add_delete_vars(const CH_Matrix_Classes::Indexmatrix& del_ind,
		      CH_Matrix_Classes::Indexmatrix& map_to_old)
  {preserves_id_flag=-1;return mdf.add_delete_vars(del_ind,map_to_old);}

  /// append @a append_dim new rows as in append_rows (if NULL, use default value) with affine offset append_rhs (if NULL use default value 0.), calls Modification::add_append_rows
  int add_append_rows(CH_Matrix_Classes::Integer append_dim,
		      const CH_Matrix_Classes::Sparsemat* append_rows,
		      const CH_Matrix_Classes::Matrix* append_rhs)
  {preserves_id_flag=-1;return mdf.add_append_rows(append_dim,append_rows,append_rhs,0);}

  ///reassign the rows as given in @a map_to_old, calls Modification::add_reassign_vars
  int add_reassign_rows(const CH_Matrix_Classes::Indexmatrix& map_to_old)
  {preserves_id_flag=-1;return mdf.add_reassign_rows(map_to_old);}

  ///delete the rows indexed by @a rows_del_ind; for each new index @a rows_map_to_old returns the old one; calls Modification::add_delete_rows
  int add_delete_rows(const CH_Matrix_Classes::Indexmatrix& rows_del_ind,
		      CH_Matrix_Classes::Indexmatrix& rows_map_to_old)
  {preserves_id_flag=-1;return mdf.add_delete_rows(rows_del_ind,rows_map_to_old);}
 
  ///multiply the current factor by the value @a times_factor (>=0, returns 1 if <0.) 
  virtual int add_apply_factor(CH_Matrix_Classes::Real times_factor)
  {if (times_factor<0.) return 1; factor*=times_factor; return 0;}

  ///add to the currecnt offset the value @a delta
  int add_offset(CH_Matrix_Classes::Real delta)
  {offset+=delta; return 0;}

 
  /// incorporate the AFTModification @a m into this one; after factor and offset are dealt with it calls Modification::incorporate
  int incorporate(const AFTModification& m);

  //@}

  //-----------------------------------------------------------------------
  /** @name  Routines for applying the collected modifications
   */

  //@{

  /// carry out the collected modifications on the given vector by calling Modification::apply_to_vars
  int apply_to_costs(CH_Matrix_Classes::Matrix*& linear_cost) const;

  /// carry out the collected modifications on the given data by calling Modification::apply_to_rows (with default value 0. for appended default entries for rhs)
  int apply_to_rows(CH_Matrix_Classes::Sparsemat*& rows,
		    CH_Matrix_Classes::Matrix*& rhs) const;

  /// multiply @a f by the factor
  int apply_to_factor(CH_Matrix_Classes::Real& f) const 
  {f*=factor;return 0;}

  /// add to @a o the offset
  int apply_to_offset(CH_Matrix_Classes::Real& o) const 
  {o+=offset;return 0;}

  /// transform the vector as if the Modification had been carried out
  const CH_Matrix_Classes::Matrix& apply_modified_transform(CH_Matrix_Classes::Matrix& out_y,
			       const CH_Matrix_Classes::Matrix &in_y,
			       const CH_Matrix_Classes::Sparsemat* arg_trafo,
			       const CH_Matrix_Classes::Matrix* arg_offset) const;
  //@}

  //-----------------------------------------------------------------------
  /** @name  Routines for querying properties of the collected modifications
   */

  //@{
  
  /// returns true if no modifications need to be executed
  virtual bool no_modification() const
  { return ((offset==0)&&(factor==1.)&&(mdf.no_modification())); }

  /// if set to true, no deletions/reassignments may be present or specified in the future, only appensions are allowed 
  int set_append_to_old(bool append_only)
  {return mdf.set_append_to_old(append_only);}
  /// returns true if this only contains appending operations and incorporating this is done with respect to the old dimension 
  bool append_to_old() const
  {return mdf.append_to_old();}

  /// returns true if linear_cost, matrix and affine rhs offset are not changed
  virtual bool only_scalars_change() const
  {return mdf.no_modification();}

  /// returns true if groundset_modifications should be ignored
  virtual bool ignore_groundset_modification() const
  {return ignore_gs_mdf;}

  /// returns true if for an AFT with argtrafo==0 the changes in the ground set reflect all modifications
  virtual bool groundset_changes_suffice_for_identity();

  /// returns true if the modifications are consistent with the AffineFunctionTransformation matrix staying the identity 
  bool preserves_identity() const;

  
 /// returns true if no columns/variables were added or deleted (allows permutations), false otherwise
  bool no_additions_or_deletions_in_vars() const;

 /// returns true if no rows were added or deleted (allows permutations), false otherwise
  bool no_additions_or_deletions_in_rows() const;

 /// returns true if all entries deleted in @a oldpoint (must be a vector of length old_vardim()) are 0,  false otherwise
  bool deleted_variables_are_zero(const CH_Matrix_Classes::Matrix& oldpoint) const
  {return mdf.deleted_variables_are_zero(oldpoint);}

 /// returns true if all entries in newpoint (must be a vector of length new_vardim()) that correspond to new variables have value 0 and false otherwise
  bool new_variables_are_zero(const CH_Matrix_Classes::Matrix& newpoint) const
  {return mdf.new_variables_are_zero(newpoint);}

  /// returns true if the values in newpoint (must be a vector of length new_vardim()) that correspond to old variables match the old values stored in oldpoint (must be a vector of length old_vardim()) and false otherwise
  bool mapped_variables_are_equal(const CH_Matrix_Classes::Matrix& newpoint,
				  const CH_Matrix_Classes::Matrix& oldpoint) const
  {return mdf.mapped_variables_are_equal(newpoint,oldpoint);}

  //@}

  //-----------------------------------------------------------------------
  /** @name  Routines implementing the abstract OracleModifcation messages
   */

  //@{

  /// returns the number of variables before modification
  int get_old_vardim() const {return mdf.old_vardim();}
  /// returns the number of variables once all stored modifications have been performed
  int get_new_vardim() const {return mdf.new_vardim();}
  /// returns the number of variables that are appended (due to later reassignmentds they may no longer be located at the end)
  int get_appended_vardim() const {return mdf.appended_vardim();}

  /// returns null if there are no index changes, otherwise the Indexmatrix pointed to is a vector whose i-th entry holds the old index of the new i-th variable (injective!), index values exceeding old_vardim() refer to newly appended variables  
  const int* get_map_to_old_variables() const
  {return (mdf.map_to_old_variables())==0?0:(const int*)(mdf.map_to_old_variables()->get_store());}

  /// incorporate the OracleModification @a m (it should only contain variable changes, but this is not checked!) into this one; calls Modification::incorporate
  int incorporate(const OracleModification& m)
  { 
    const AFTModification* aftm=dynamic_cast<const AFTModification*>(&m);
    if (aftm) return incorporate(*aftm); 
    return 1;
  }

    /// returns a new object on the heap, that allows to incorporate this but starts off from a function whose input argument dimension is old_var_dim 
  OracleModification* new_initial_oraclemodification(int old_var_dim) const 
  {return new AFTModification(old_var_dim,old_rowdim());}
       
    /// append append_dim further variables at the end of the argument vector (without specifying their effect, they may e.g. be ignored by the function)
  int add_append_variables(int append_dim)
  {return add_append_vars(CH_Matrix_Classes::Integer(append_dim),0,0);}

    /// reorder and resize the variables as given by the first new_dim entries of map_to_old_indices; each former index may appear at most once in this list, so new_dim < get_new_vardim()
    virtual int add_reassign_variables(int new_dim,const int* map_to_old_indices)
  { CH_Matrix_Classes::Indexmatrix maptoold(CH_Matrix_Classes::Integer(new_dim),1,map_to_old_indices); 
    return add_reassign_vars(maptoold);}

  //@}

  //-----------------------------------------------------------------------
  /** @name  Routines for retrieving the detailed collected modifications
   */

  //@{

  /// returns the number of variables before modification (given on initialization)
  CH_Matrix_Classes::Integer old_vardim() const 
  {return mdf.old_vardim();}
  /// returns the number of variables once all stored modifications have been performed
  CH_Matrix_Classes::Integer new_vardim() const
  {return mdf.new_vardim();}
  /// returns the number of variables that are appended (due to later reassignmentds they may no longer be located at the end)
   CH_Matrix_Classes::Integer appended_vardim() const
  {return mdf.appended_vardim();}

  /// returns the number of rows before modification (given on initialization)
  CH_Matrix_Classes::Integer old_rowdim() const
  {return mdf.old_rowdim();}
  /// returns the number of rows once all stored modifications have been performed
   CH_Matrix_Classes::Integer new_rowdim() const
  {return mdf.new_rowdim();}
  /// returns the number of rows that are appended (due to later reassignments they may no longer be located at the end)
  CH_Matrix_Classes::Integer appended_rowdim() const
  {return mdf.appended_rowdim();}

  /// returns null if there are no index changes, otherwise the Indexmatrix pointed to is a vector whose i-th entry holds the old index of the new i-th variable (injective!), index values exceeding old_vardim() refer to newly appended variables  
  const CH_Matrix_Classes::Indexmatrix* map_to_old_variables() const
  {return mdf.map_to_old_variables();}
 
  /// returns null if there were no deletions, otherwise the Indexmatrix pointed to is a vector holding the deleted old variable indices in increasing order   
   const CH_Matrix_Classes::Indexmatrix* deleted_var_indices() const
  {return mdf.deleted_var_indices();}

  /// returns null if no variables were added, otherwise the Indexmatrix pointed to is a vector holding the new indices of the new variables in increasing order
  const CH_Matrix_Classes::Indexmatrix* new_var_indices() const
  {return mdf.new_var_indices();}

  /// returns null if there are no index changes, otherwise the Indexmatrix pointed to is a vector whose i-th entry holds the old index of the new i-th row (injective!), index values exceeding old_rowdim() refer to newly appended rows  
  const CH_Matrix_Classes::Indexmatrix* map_to_old_rows() const
  {return mdf.map_to_old_rows();}

  /// returns null if there were no deletions, otherwise the Indexmatrix pointed to is a vector holding the deleted old row indices in increasing order   
  const CH_Matrix_Classes::Indexmatrix* deleted_row_indices() const
  {return mdf.deleted_row_indices();}
  
  /// returns null if no variables were added, otherwise the Indexmatrix pointed to is a vector holding the new indices of the new rows in increasing order
  const CH_Matrix_Classes::Indexmatrix* new_row_indices() const
  {return mdf.new_row_indices();}

  /// returns the value to be added to the offset
  CH_Matrix_Classes::Real get_additional_offset() const
  {return offset;}

  /// returns null if nothing or default values have to be appended, otherwise it points to a sparse matrix whose columns need to be appended to the matrix
  CH_Matrix_Classes::Real get_additional_factor() const
  {return factor;}

  /// returns null if nothing or default values have to be appended, otherwise it points to a sparse matrix whose columns need to be appended to the matrix
  const CH_Matrix_Classes::Sparsemat* get_append_cols() const
  {return mdf.get_var_append_cols();}

  /// returns null if nothing or default values have to be appended, otherwise it points to a matrix whose entries need to be appended to the cost vector
  const CH_Matrix_Classes::Matrix* get_append_costs() const
  {return mdf.get_var_append_costs();}

  /// returns null if nothing or default values have to be appended, otherwise it points to a sparse matrix whose rows need to be appended to the matrix
  const CH_Matrix_Classes::Sparsemat* get_append_rows() const
  {return mdf.get_row_append_mat();}

  /// returns null if nothing or default values have to be appended, otherwise it points to a matrix whose rows need to be appended to the argument offset
  const CH_Matrix_Classes::Matrix* get_append_rhs() const
  {return mdf.get_row_append_rhslb();}

  //@}

  //-----------------------------------------------------------------------
  /** @name  Output settings
   */

  //@{

  /// see CBout::set_out
   void set_out(std::ostream* out=0,int print_level=1)
  {CBout::set_out(out,print_level);mdf.set_out(out,print_level);}

  //@}

};



  //@}

}

#endif

