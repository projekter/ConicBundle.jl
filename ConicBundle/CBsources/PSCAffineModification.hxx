/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/PSCAffineModification.hxx
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




#ifndef CONICBUNDLE_PSCAffineMODIFICATION_HXX
#define CONICBUNDLE_PSCAffineMODIFICATION_HXX


/**  @file PSCAffineModification.hxx
    @brief Header declaring the class ConicBundle::PSCAffineModification
    @version 1.0
    @date 2016-09-29
    @author Christoph Helmberg
*/

#include <map>
#include "CBout.hxx"
#include "CBSolver.hxx"
#include "SparseCoeffmatMatrix.hxx"
#include "ModificationBase.hxx"

namespace ConicBundle {

  class PSCPrimal;

  /** @ingroup implemented_psc_oracle

  */
  //@{



  /** @brief class for collecting and organizing a sequence of changes to block diagonal symmetric affine matrix functions so that it can be carried out in one step later on;

     The general setting assumed here is to support changes for a function of the form

     \f[ \left[\begin{array}{ccc}C_1& & 0\\&\ddots&\\0& &C_k\end{array}\right]+\sum_{j=1}^m y_j\left[\begin{array}{ccc}A_{1,j}& & 0\\&\ddots&\\0& &A_{k,j}\end{array}\right] \f]

     The data is thought as given by a (column) vector of pointers to matrices C_i and
     a k times m matrix calA where row i holds the pointers to matrices A_{i,j} for variable y_j, j=1,...,m. Note, because matrices are considered potentially hug objects
     this class never copies or deletes any matrices but only deals with the
     pointers to the matrices!

     Supported are an arbitrary sequence of

     - appending variables (add_append_vars()), i.e. appending further
       "columns" of blockdiagonal matrices  to \f$A\f$.

     - reassigning the variables (add_reassign_vars()) via an index vector whose
       entries specify which variable will be moved to this position. Each
       existing variable index may appear at most once in this index vector, but
       not all indices need to appear. If they do not appear they will be
       deleted.

     - deleting variables (add_delete_vars()); the indices to be deleted are
       specfied by an index vector, which is then used to generate an index
       vector as a map to the old indices and then uses this with the previous
       reassign variables.

     - appending diagonal blocks (add_append_blocks()), i.e. appending further rows
       to A and further elements to C.

     - reassigning the blocks (add_reassign_blocks()) in the same style as
       reassigning variables

     - deleting diagonal blocks (add_delete_blocks()) in the same style as deleting
       variables

     Each next modification has to be given/added with respect to the
     current virtual number of variables and blocks as if the
     compiled modifications had been carried out already. If certain
     elements are not specified, they are assumed to be zero

     The modifications are compiled within this class (and not yet executed)
     to one bulk of data of appending information for variables
     and blocks and reassignment information, so that the entire sequence of
     transformations is then achieved in the following sequence of four steps
     (the sequence is important here!)

     1. appending the compiled variable (and column) data

     3. reassigning the variables (and columns) by the compiled variable reassignment map

     4. appending the compiled block (row) data

     5. reassigning the blocks (rows) by the compiled block reassignment map

     Further modifications can be incorporated into this one if the
     dimensions match the current values.

     In the end the routines apply_to_vars() and apply_to_blocks() can
     be used to carry out the modifications on a SparseCoeffmatMatrix.

   */

  class PSCAffineModification : public ModificationBase, public OracleModification {
  private:

    /// if true, only appending operations are allowed and incorporating this has a different effect
    bool append_only;


    //first execute the changes in the variables then in the constrains

    //------ changes in variable space
    //first carry out the append step, afterwards the delete/reassign step

    //-----------------------------------------------------------------------
    /** @name Modification information for variables/columns
     */
     //@{

     ///initial number (or dimension) of variables
    CH_Matrix_Classes::Integer var_olddim;
    ///number (or dimension) of variables after all listed modifications have been applied
    CH_Matrix_Classes::Integer var_newdim;

    ///number of variables appended
    CH_Matrix_Classes::Integer var_append_dim;

    ///columns to be appended for these new variables
    SparseCoeffmatMatrix var_append;

    ///indices of variables that will be deleted by the reassignment of indices by map_to_old after all additions have been carried out
    CH_Matrix_Classes::Indexmatrix* var_del_ind;
    ///the variables are rearranged so that the new index i had previously (after additions and before deletion) the index map_to_old(i) 
    CH_Matrix_Classes::Indexmatrix* var_map_to_old;
    ///in the end the appended new variables have these positions 
    CH_Matrix_Classes::Indexmatrix* var_new_ind;

    //@}

    //------ changes in the block (the rows)
    //first carry out the append step, afterwards the reassign/delete step

    //-----------------------------------------------------------------------
    /** @name Modification information for blocks (rows)
     */
     //@{

     ///initial number of blocks and their dimension
    CH_Matrix_Classes::Indexmatrix block_olddim;
    /// number of blocks and their dimension after all listed modifications have been applied
    CH_Matrix_Classes::Indexmatrix block_newdim;

    /// number of blocks appended with their dimension
    CH_Matrix_Classes::Indexmatrix block_append_dim;
    ///offset blocks to be appended;
    SparseCoeffmatMatrix offset_append;
    ///rows to be appended;
    SparseCoeffmatMatrix block_append;
    ///indices of rows that will be deleted by the reassignment of indices by map_to_old after all additions have been carried out
    CH_Matrix_Classes::Indexmatrix* block_del_ind;
    ///the rows are rearranged so that the new index i had previously (after additions and before deletions) the index map_to_old(i) 
    CH_Matrix_Classes::Indexmatrix* block_map_to_old;
    ///in the end the appended new rows have these positions 
    CH_Matrix_Classes::Indexmatrix* block_new_ind;

    //@}

    //-----------------------------------------------------------------------
    /** @name Modification information for primal aggregation
     */
     //@{
     ///
    bool reset_primal;
    /// primal to reset to
    PSCPrimal* generating_primal;
    /// if set to true, no extension is performed on newly added variables
    bool skip_extension;

    //@}

  public:

    //-----------------------------------------------------------------------
    /** @name Constructors and initialization
     */
     //@{

     ///
    virtual ~PSCAffineModification();


    /// calls clear() with these parameters
    PSCAffineModification(CH_Matrix_Classes::Integer var_olddim,
      const CH_Matrix_Classes::Indexmatrix& block_olddim,
      const CBout* cb = 0, int incr = 0);

    /** @brief resets all variables so that the object to be modified has
        starting size var_olddim (number of variables) and block_olddim
        (number of rows) and no modifications

        The actual old data is not needed at this point,
        the changes on it will be collected and excuted in the
        routines apply_to_vars and apply_to_blocks

        Setting the parameter ensure_start_val_box_feasibility to true
        will cause the algorithm to check in add_append_vars() whether
        the input values are within the given bounds and in
        apply_to_vars() it will project all start_values onto the bounds
        for all, old and new, indices (which might have been changed by
        then).  If it is false (default), all values will be accepted as
        given.

        Setting the parameter ensure_bounds_consistency to true
        (default) will raise errors in add_append_vars() and in
        apply_to_vars() whenever there are lower bounds greater than the
        respective upper bounds so as to avoid trivial
        infeasibilities. This check is omitted if set to false.

        The remaining values give the values of plus and minus infinity
        the no bounds should exceed. These are the default values at the
        same time (maybe it might be good to have a separate default value,
        but this is not implemented here).
    */
    int clear(CH_Matrix_Classes::Integer var_olddim,
      const CH_Matrix_Classes::Indexmatrix& block_olddim);

    //@}

    //-----------------------------------------------------------------------
    /** @name  Routines for adding modifications
     */
     //@{


     /** @brief append information on new variables at the respective ends

        @param append_dim
               number of variables (or columns ) to
         be appended

        @param append_cols
               if NULL, append zero matrices, otherwise it must point to
         a sparse matrix of size new_rowdim() times @a append_dim
               that is to be appended to the constraint matrix on the right

        @return number of errors; if errors occured, none of the new changes are performed
     */
    int add_append_vars(CH_Matrix_Classes::Integer append_dim,
      const SparseCoeffmatMatrix* append_cols);

    /** @brief reassign the current variable indices (with modifications) as specified by @a map_to_old

        @a map_to_old must specify an injective map (no two values match)
        into indices 0 up to new_vardim()-1 (not all need to appear). The variable getting index
        i (for i=0 to map_to_old.dim()-1) is the variable with current
        index map_to_old(i) (current refers to considering all previous
        modifications as having been carried out already). The return
        value is the number of errors in @a map_to_old. If such occured, this
        reassign is not performed.
     */
    int add_reassign_vars(const CH_Matrix_Classes::Indexmatrix& map_to_old);

    /** @brief delete the variables indexed by the vector del_ind and
        return the index changes of the others in a vector map_to_old

        @a del_ind need not be ordered in any way, but any index in 0 to
        new_vardim()-1 may appear at most once. On output the dimension
        of the column vector @a map_to_old gives the number of remaining
        variables and its entry i holds the index the variable with new
        index i had before the deletion. The return value is the number
        of errors in @a del_ind. If such occured, this deletion is
        not performed and @a map_to_old may contain garbage.
    */
    int add_delete_vars(const CH_Matrix_Classes::Indexmatrix& del_ind,
      CH_Matrix_Classes::Indexmatrix& map_to_old);




    /** @brief append information on new rows at the respective ends

       @param append_dim
              the dimension gives the number of (diagonal) blocks and each
        entry the order of the block to be appended

       @param append_blocks
              if NULL, append zero matrics, otherwise it must point
              to a sparse matrix of size @a append_dim times new_vardim()
              that is to be appended to the constraint matrix below.

       @param append_offsets
              if NULL, append zero matrices, otherwise it must point to
        a column vector of size @a append_dim that is to be appended
              to the vector of offsets


       @return number of errors; if errors occured, none of the new changes are performed
    */
    int add_append_blocks(const CH_Matrix_Classes::Indexmatrix& append_dim,
      const SparseCoeffmatMatrix* append_offsets,
      const SparseCoeffmatMatrix* append_blocks);

    /** @brief reassign the current row indices (with modifications) as specified by @a map_to_old

        @a map_to_old must specify an injective map (no two values
        match) into indices 0 up to new_rowdim()-1 (not all need to
        appear).  The row getting index i (for i=0 to
        map_to_old.dim()-1) is the row with current index map_to_old(i)
        (current refers to considering all previous modifications as
        having been carried out already). The return value is the number
        of errors in @a map_to_old. If such occured, this reassign is
        not performed.
     */
    int add_reassign_blocks(const CH_Matrix_Classes::Indexmatrix& map_to_old);

    /** @brief delete the rows indexed by the vector del_ind and
        return the index changes of the others in a vector map_to_old

        @a del_ind need not be ordered in any way, but any index in 0 to
        new_rowdim()-1 may appear at most once. On output the dimension
        of the column vector @a map_to_old gives the number of remaining
        rows and its entry i holds the index the row with new index i
        had before the deletion. The return value is the number of
        errors in @a del_ind. If such occured, this deletion is not
        performed and @a map_to_old may contain garbage.
    */
    int add_delete_blocks(const CH_Matrix_Classes::Indexmatrix& del_ind,
      CH_Matrix_Classes::Indexmatrix& map_to_old);


    /** @brief replace the current generating primal by new_generating_primal (on the heap, will be deleted here). It may be NULL to switch of generating primals.

        Any changes here cause the deletion of all current aggregate minorants.
        If not NULL, the object pointed to is on the heap and control over it
        is passed over to *this, so *this will make sure it is deleted at due
        time.

        PSCPrimal gives no information on the order of the matrices involved, so no consistency checks are done here.
    */
    int add_reset_generating_primal(PSCPrimal* new_generating_primal);

    /** @brief if this time no extension is possible for newly added variables
        with the availabel generating primal, set this to true;
     */
    int set_skip_extension(bool skip);

    //@}

    //-----------------------------------------------------------------------
    /** @name  Routines for applying the collected modifications
     */

     //@{

     /** @brief carry out the collected modifications on the data describing the PSCAffineFunction

         If a specific parameter is NULL, no changes are performed on it,
         if it is not null, it must point to a SparseCoeffmatMatrix
         whose sizes correspond to the "old" data. Then the following
         operations will be performed on it in this sequence:

         1. new variables (columns) are appended to the matrix

         2. if reassignment information is given, the columns are
            mapped/reorderd as given by *map_to_old_variables()

         3. new blocks (rows) are appended to the offset and the matrix

         4. if reassignment information is given, the rows are
            mappen/reorderd as given by *map_to_old_blocks()

         @param[in,out] offset
            if not NULL, this points to the offset vector.

         @param[in,out] matrix
            if not NULL, this points to the old matrix.

         @return the number of dimension errors of non NULL inputs,
            if any, no modifications are made to any inputs.
     */
    int apply_to_PSCAffine(SparseCoeffmatMatrix* offset, SparseCoeffmatMatrix* matrix) const;

    //@}

    //-----------------------------------------------------------------------
    /** @name  Routines for querying properties of the collected modifications
     */

     //@{

     /// returns true if no modifications need to be executed
    bool no_modification() const;
    /// if set to true, no deletions/reassignments may be present or specified in the future, only appensions are allowed 
    int set_append_to_old(bool append_only);
    /// returns true if this only contains appending operations and incorporating this is done with respect to the old dimension 
    bool append_to_old() const {
      return append_only;
    }
    /// returns true if all entries deleted in @a oldpoint (must be a vector of length old_vardim()) or the corresponding entries in the old matrix are 0 and false otherwise
    bool deleted_variables_are_zero(const CH_Matrix_Classes::Matrix& oldpoint, const SparseCoeffmatMatrix& oldmat) const;
    /// returns true if for all indices of new variables the entries in newpoint (must be a vector of length new_vardim()) or the matrices in newmat are 0 and false otherwise
    bool new_variables_are_zero(const CH_Matrix_Classes::Matrix& newpoint, const SparseCoeffmatMatrix& newmat) const;
    /// returns true if the values in newpoint (must be a vector of length new_vardim()) that correspond to old variables match the old values stored in oldpoint (must be a vector of length old_vardim()) and false otherwise
    bool mapped_variables_are_equal(const CH_Matrix_Classes::Matrix& newpoint,
      const CH_Matrix_Classes::Matrix& oldpoint) const;
    /// returns true if some modifications are performed on the block structure
    bool variable_modifications() const {
      return ((var_append_dim > 0) || (var_map_to_old));
    }
    /// returns true if some modifications are performed on the block structure
    bool block_modifications() const {
      return ((block_append_dim.dim() > 0) || (block_map_to_old));
    }

    //@}

    //-----------------------------------------------------------------------
    /** @name  Routines for retrieving the detailed collected modifications
     */

     //@{

     /// returns the number of variables before modification (given on initialization)
    CH_Matrix_Classes::Integer old_vardim() const {
      return var_olddim;
    }
    /// returns the number of variables once all stored modifications have been performed
    CH_Matrix_Classes::Integer new_vardim() const {
      return var_newdim;
    }
    /// returns the number of variables that are appended (due to later reassignmentds they may no longer be located at the end)
    CH_Matrix_Classes::Integer appended_vardim() const {
      return var_append_dim;
    }

    /// returns the number of rows before modification (given on initialization)
    const CH_Matrix_Classes::Indexmatrix& old_blockdim() const {
      return block_olddim;
    }
    /// returns the number of rows once all stored modifications have been performed
    const CH_Matrix_Classes::Indexmatrix& new_blockdim() const {
      return block_newdim;
    }
    /// returns the number of rows that are appended (due to later reassignments they may no longer be located at the end)
    const CH_Matrix_Classes::Indexmatrix& appended_blockdim() const {
      return block_append_dim;
    }


    /// returns null if there are no index changes, otherwise the Indexmatrix pointed to is a vector whose i-th entry holds the old index of the new i-th variable (injective!), index values exceeding old_vardim() refer to newly appended variables  
    const CH_Matrix_Classes::Indexmatrix* map_to_old_variables() const {
      return var_map_to_old;
    }

    /// returns null if there were no deletions, otherwise the Indexmatrix pointed to is a vector holding the deleted old variable indices in increasing order   
    const CH_Matrix_Classes::Indexmatrix* deleted_var_indices() const {
      return var_del_ind;
    }

    /// returns null if no variables were added, otherwise the Indexmatrix pointed to is a vector holding the new indices of the new variables in increasing order
    const CH_Matrix_Classes::Indexmatrix* new_var_indices() const {
      return var_new_ind;
    }

    /// returns null if there are index changes, otherwise the Indexmatrix pointed to is a vector whose i-th entry holds the old index of the new i-th row (injective!), index values exceeding old_rowdim() refer to newly appended rows  
    const CH_Matrix_Classes::Indexmatrix* map_to_old_blocks() const {
      return block_map_to_old;
    }

    /// returns null if there were no deletions, otherwise the Indexmatrix pointed to is a vector holding the deleted old row indices in increasing order   
    const CH_Matrix_Classes::Indexmatrix* deleted_block_indices() const {
      return block_del_ind;
    }

    /// returns null if no rows were added, otherwise the Indexmatrix pointed ato is a vector holding the new indices of the new rows in increasing order
    const CH_Matrix_Classes::Indexmatrix* new_block_indices() const {
      return block_new_ind;
    }

    /// returns null if nothing or default values have to be appended, otherwise it points to a sparse matrix whose columns need to be appended to the matrix
    const SparseCoeffmatMatrix& get_var_append() const {
      return var_append;
    }
    /// returns null if nothing or default values have to be appended, otherwise it points to a sparse matrix whose rows need to be appended to the matrix
    const SparseCoeffmatMatrix& get_block_append() const {
      return block_append;
    }
    /// returns null if nothing or default values have to be appended, otherwise it points to a matrix whose entries need to be appended to the right hand side lower bounds vector
    const SparseCoeffmatMatrix& get_offset_append() const {
      return offset_append;
    }
    /// returns true if the generating primal is to be replaced by the one stored here
    bool get_reset_primal() const {
      return reset_primal;
    }
    /// returns the generating primal pointer stored here (may be NULL); if get_reset_primal() is true, the PSCAffineFunction should either delete its generating primal (NULL) or replace its generating primal by a clone of this one
    const PSCPrimal* get_generating_primal() const {
      return generating_primal;
    }

    /// returns true if the generating primal is to be replaced by the one stored here
    bool get_skip_extension() const {
      return skip_extension;
    }

    //@}


    //-----------------------------------------------------------------------
    /** @name  Routines implementing the abstract OracleModifcation messages
     */

     //@{

     /// returns the number of variables before modification
    int get_old_vardim() const {
      return old_vardim();
    }
    /// returns the number of variables once all stored modifications have been performed
    int get_new_vardim() const {
      return new_vardim();
    }
    /// returns the number of variables that are appended (due to later reassignmentds they may no longer be located at the end)
    int get_appended_vardim() const {
      return appended_vardim();
    }

    /// returns null if there are no index changes, otherwise the Indexmatrix pointed to is a vector whose i-th entry holds the old index of the new i-th variable (injective!), index values exceeding old_vardim() refer to newly appended variables  
    const int* get_map_to_old_variables() const {
      return (map_to_old_variables()) == 0 ? 0 : (const int*)(map_to_old_variables()->get_store());
    }

    /** @brief add the modification specified in @a m on top of
        the modifications collected so far

        If m is in fact an PSCAffineModification,
        the old_vardim() of modification @a m must be
        identical to new_vardim() of this and
        old_rowdim() of modification @a m must be identical to
        new_rodim() of this.  The return value is the number of
        errors in this respect. If such occured, this incorporation
        is not performed.

        A general OracleModification @a m should only contain variable
        changes (this is not checked) and appending variables appends
        zero blocks.
    */
    int incorporate(const OracleModification& m);

    /// returns a new object on the heap, that allows to incorporate this but starts off from a function whose input argument dimension is old_var_dim 
    OracleModification* new_initial_oraclemodification(int old_var_dim) const {
      return new PSCAffineModification(CH_Matrix_Classes::Integer(old_var_dim), block_olddim, this, 0);
    }

    /// append append_dim further variables at the end of the argument vector (without specifying their effect, they may e.g. be ignored by the function)
    int add_append_variables(int append_dim) {
      return add_append_vars(CH_Matrix_Classes::Integer(append_dim), 0);
    }

    /// reorder and resize the variables as given by the first new_dim entries of map_to_old_indices; each former index may appear at most once in this list, so new_dim < get_new_vardim()
    virtual int add_reassign_variables(int new_dim, const int* map_to_old_indices) {
      CH_Matrix_Classes::Indexmatrix maptoold(CH_Matrix_Classes::Integer(new_dim), 1, map_to_old_indices);
      return add_reassign_vars(maptoold);
    }

    //@}

    //-----------------------------------------------------------------------
    /** @name  Output settings
     */

     //@{

     /// see CBout::set_cbout
    void set_cbout(const CBout* out, int incr = -1) {
      CBout::set_cbout(out, incr); var_append.set_cbout(this); block_append.set_cbout(this);
    }

    /// see CBout::set_out
    void set_out(std::ostream* out = 0, int print_level = 1) {
      CBout::set_out(out, print_level); set_cbout(this, 0);
    }

    //@}

  };





  //@}

}

#endif

