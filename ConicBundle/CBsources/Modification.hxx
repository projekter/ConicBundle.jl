/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/Modification.hxx
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




#ifndef CONICBUNDLE_MODIFICATION_HXX
#define CONICBUNDLE_MODIFICATION_HXX


/**  @file Modification.hxx
    @brief Header declaring the class ConicBundle::Modification
    @version 1.0
    @date 2014-08-29
    @author Christoph Helmberg
*/

#include <map>
#include "ModificationBase.hxx"
#include "CBSolver.hxx"
#include "sparsmat.hxx"

namespace ConicBundle {

  /** @defgroup dynamic_modification_support dynamic modification support for adding/deleting/modifying variables, oracles and constraints
     @brief Most classes are specializations/reinterpretations of the class Modification

    Dynamic additions and deletions of variables during the optimization
    process occur frequently in Lagrangian relaxation of primal linear
    cutting planes that are separated on basis of a current approximate
    solution. In the context of sums of convex functions they frequently
    represent coupling ressource constraints between function subgroups
    that compete for this ressource. The typical features are thus

    - variables have to be added/deleted that are common to some of the functions
      but do not appear in all of them

    - cutting models of unaffected functions can and should be preserved

    - if the new variable is a multiplier of a new linear constraint, it
      typically starts with value zero, so even for the changed
      functions the function values are not affected. Furthermore, most
      subgradients arise as primal violation \f$b-Ax\f$ of the relaxed
      constraints \f$Ax=b\f$ for some specific PrimalData \f$x\f$. If
      this PrimalData is still available for an old subgradient, a new
      constraints \f$a_i^\top x=b_i\f$ in \f$Ax=b\f$ just gives rise to
      a new coordinate \f$b_i-a_i^\top x\f$ in this subgradient and the
      new entries can be filled in. If this can be done for all new
      coordinates of the subgradients in use, then also the cutting
      model of the modified functions can be preserved in spite of the
      changes.

    - in bundle methods it does not hurt too much to discard a model after
      a descent step, but if models cannot be preserved throughout sequences
      of null steps convergence is seriously endangered and the quality
      of the primal approximations is poor as well.

    If the Lagrangian relaxation of a changing set of linear constraints
    is modelled via an AffineFunctionTransformation on the support
    function of the primal groundset, then all aspects can be handled
    and checked by ConicBundle directly. Indeed, in this case, the
    function representing the ground set will not be affected by the
    changes in the argument function. In several applications, however,
    the changes might also affect the function itself (e.g. if it uses a
    column generation approach on the primal side) and ConicBundle then
    needs the cooperation of the function in order to keep a model alive
    as long as possible.  The attempt to support this in any thinkable
    form entails rather cumbersome details.

    - All changes, those affecting any AffineFunctionTransformation or
      any Oracle, have to be communicated and excuted through the solver
      (and from this through its respective SumBlockModel) via so called
      modifications. The solver then adapts the groundset and possibly
      modifies some AffineFunctionTransformation automatically so that
      only those functions are affected that need to be. The models of
      affected functions communicate the local changes via modifications
      to their respective oracles by calling
      e.g. FunctionOracle::apply_modification().  Only if in this the
      oracle asserts that the modifcations may allow to preserve the
      cutting model by possibly providing extension routines for
      previous subgradients, will the model try to save the cutting
      model.

    - All modifications are accumulated as long as possible before they
      are actually excuted (via a routine called e.g. Modification::incorporate())
      This should help to keep the previously computed data valid as long
      as possible and may help to continue even if one modification turns out
      to fail the consistency checks in the accumulation process. Thus
      the function may not yet change even if a modification has been entered
      for it. Still, further modifications always have to be formulated as
      if the previous ones are already executed for the accumulation process
      to work independent of the actual point in time, when part or all of
      the modifcations are carried out.

    - The modifications induced on the Groundset by
      MatrixCBSolver::append_variables(), MatrixCBSolver::delete_variables()
      and MatrixCBSolver::reassign_variables() are generated by
      MatrixCBSolver automatically (via GroundsetModification and
      LPGroundsetModification), but a function specific
      FunctionObjectModification to an AffineFunctionTransformation (via
      AFTModification) or to the oracles themselves (via an appropriate
      implementation of an OracleModification and e.g.
      MatrixFunctionOracle::apply_modification()) have to be given
      explicitly and consistently! Read the instructions there
      carefully.

    - Changing an oracle during a call to e.g. MatrixFunctionOracle::evaluate()
      because of a primal column generation approach will typically also allow
      to preserve the model if a PrimalExtender is returned in the last
      argument to the calling SumBlockModel, so that it can adapt the PrimalData
      of past subgradients that are still in use there.

    Examples modifications for the various oracle classes are
    NNCBoxSupportModification, SOCSupportModification, PSCAffineModification and
    AFTModification. The classes Modification and GroundsetModification come
    in handy in this.

  */
  //@{

  /** @brief base class for collecting and organizing a sequence of changes to linear data so that it can be carried out in one step later on; this class comprises all  features, derived ones are specializations with partial reinterpretations.

     The general setting assumed here is to support changes for a problem of the form

     \f[ \min\{c^\top x: a \le Ax\le b, l\le x \le u\} \f]

     where for some given @a old dimensions \f$\bar n,\bar m\in\mathbf{N}_0\f$
     the data are \f$c,l,u\in\mathbf{R}^{\bar n}\f$ (linear cost,lower and upper
     bounds on the variables), \f$a,b\in\mathbf{R}^{\bar m}\f$ (lower and upper
     bounds on the constraint values), \f$A\in\mathbf{R}^{\bar m\times\bar n}\f$
     (the constraint matrix). In addition the class supports manipulating a
     vector of starting values \f$\bar x\in\mathbf{R}^{\bar n}\f$. These
     not be feasible for the data (for some infeasible methods this might even
     be desirable), so their values are not checked. If desired, checking
     feasibility with respect to \f$l,u\f$ can be turned on.

     Supported are an arbitrary sequence of

     - setting new values for lower/upper bounds \f$a,b,l,u\f$ on selected indices
       (add_set_lb(), add_set_ub(), add_set_rhslb(), add_set_rhsub()), by the
       default consistency of the bounds is checked when applying the changes
       (typically it cannot be done before in lack of the necessary information
       regarding the opposite bound or would require simultaneous setting of
       lower and upper bounds which most applications do not want to do).

     - appending variables (add_append_vars()), i.e. appending further elements to
       \f$c,l,u,\bar x\f$ as well as further columns to \f$A\f$.

     - reassigning the variables (add_reassign_vars()) via an index vector whose
       entries specify which variable will be moved to this position. Each
       existing variable index may appear at most once in this index vector, but
       not all indices need to appear. If they do not appear they will be
       deleted. This affects \f$c,l,u,\bar x\f$ and the columns of \f$A\f$.

     - deleting variables (add_delete_vars()); the indices to be deleted are
       specfied by an index vector, which is then used to generates an index
       vector as a map to the old indices and then uses this with the previous
       reassign variables.

     - appending constraints rows (add_append_rows()), i.e. appending further rows
       to \f$A\f$ and further elements to \f$a,b\f$.

     - reassigning the constraint rows (add_reassign_rows()) in the same style as
       reassigning variables

     - deleting constraint rows (add_delete_rows()) in the same style as deleting
       variables

     Each next modification has to be given/added with respect to the
     current virtual number of variables and constraint rows as if the
     compiled modifications had been carried out already. If certain
     elements are not specified, they are assumed to have the following
     default values (if used at all; they can also be set to other
     values on construction or by clear()):

     - 0 for the linear cost \f$c,\bar x\f$ and matrix elements of \f$A\f$

     - CB_minus_infinity for \f$a,l\f$

     - CB_plus_infinity for \f$b,u\f$

     The modifications are compiled within this class (and not yet executed)
     to one bulk of data of bound inromations, appending information for variables
     and rows and reassignment information, so that the entire sequence of
     transformations is then achieved in the following sequence of four steps
     (the sequence is important here!)

     1. setting the requested bounds to their new values

     2. appending the compiled variable (and column) data

     3. reassigning the variables (and columns) by the compiled variable reassignment map

     4. appending the compiled constraint row data

     5. reassigning the constraint rows by the compiled row reassignment map

     Another modification can be incorporated into this one if either its
     starting dimensions match the current values or if the modification starts
     from the same old dimensions and only appends information. In the latter
     case any "missed" appensions of *this are executed with the neutral
     element, likewise the reassignments of *this are executed and then the
     resultig modified parts of the modification to incorporate are appended.

     In the end the routines apply_to_vars() and apply_to_rows() can
     be used to carry out the modifications.

   */

  class Modification : public ModificationBase {
  private:
    /// for storing updates of bounds
    typedef std::map<CH_Matrix_Classes::Integer, CH_Matrix_Classes::Real> Realmap;

    /// if true, only appending operations are allowed and incorporating this has a different effect
    bool append_only;

    //first execute the changes in the variables then in the constrains

    //------ changes in variable space
    //first carry out the append step, afterwards the delete/reassign step, finally the bound changes

    //-----------------------------------------------------------------------
    /** @name Modification information for variables/columns
     */
     //@{

     ///initial number of variables
    CH_Matrix_Classes::Integer var_olddim;
    ///number of variables after all listed modifications have been applied
    CH_Matrix_Classes::Integer var_newdim;

    ///first set selected new values in variable lower bounds; if NULL there are no changes 
    Realmap* var_set_lb;
    ///and set selected new values in variable upper bounds; if NULL there are no changes 
    Realmap* var_set_ub;
    ///number of variables appended
    CH_Matrix_Classes::Integer var_append_dim;
    ///lower bounds on appended variables; if NULL, use CB_minus_infinity
    CH_Matrix_Classes::Matrix* var_append_lb;
    ///upper bounds on appended variables; if NULL, use CB_plus_infinity
    CH_Matrix_Classes::Matrix* var_append_ub;
    ///columns to be appended for these new variables; if NULL, use zero
    CH_Matrix_Classes::Sparsemat* var_append_cols;
    ///starting values for appended variables; if NULL, use the projection of 0 on [lb,ub]
    CH_Matrix_Classes::Matrix* var_start_val;
    ///linear costs for appended variables; if NULL, use zero
    CH_Matrix_Classes::Matrix* var_append_costs;
    ///indices of variables that will be deleted by the reassignment of indices by map_to_old after all additions have been carried out
    CH_Matrix_Classes::Indexmatrix* var_del_ind;
    ///the variables are rearranged so that the new index i had previously (after additions and before deletion) the index map_to_old(i) 
    CH_Matrix_Classes::Indexmatrix* var_map_to_old;
    ///in the end the appended new variables have these positions 
    CH_Matrix_Classes::Indexmatrix* var_new_ind;

    //@}

    //------ changes in row space
    //first carry out the append step, afterwards the reassign/delete step

    //-----------------------------------------------------------------------
    /** @name Modification information for rows
     */
     //@{

     ///initial number of rows
    CH_Matrix_Classes::Integer row_olddim;
    /// number of rows after all listed modifications have been applied
    CH_Matrix_Classes::Integer row_newdim;

    ///first set selected new values in rhs lower bounds; if NULL there are no changes 
    Realmap* row_set_rhslb;
    ///and set selected new values in rhs upper bounds; if NULL there are no changes 
    Realmap* row_set_rhsub;
    ///number of rows appended
    CH_Matrix_Classes::Integer row_append_dim;
    ///rows to be appended as these new rows; if NULL, use zero
    CH_Matrix_Classes::Sparsemat* row_append_mat;
    ///new lower bound on rhs values for these new rows; if NULL, use CB_minus_inifinity
    CH_Matrix_Classes::Matrix* row_append_rhslb;
    ///new upper bound on rhs values for these new rows; if NULL, use CB_plus_infintiy
    CH_Matrix_Classes::Matrix* row_append_rhsub;
    ///indices of rows that will be deleted by the reassignment of indices by map_to_old after all additions have been carried out
    CH_Matrix_Classes::Indexmatrix* row_del_ind;
    ///the rows are rearranged so that the new index i had previously (after additions and before deletions) the index map_to_old(i) 
    CH_Matrix_Classes::Indexmatrix* row_map_to_old;
    ///in the end the appended new rows have these positions 
    CH_Matrix_Classes::Indexmatrix* row_new_ind;

    //@}

    //-----------------------------------------------------------------------
    /** @name parameters for consistency checking and default values
     */
     //@{
    bool enforce_bounds_consistency; ///< produce errors if bounds cause infeasibilities when used in apply_to_vars(), default true
    bool enforce_start_val_box_feasibility; ///< produce errors if given starting values violate bounds; in apply_to_vars() project default values onto the bound interval, default false
    CH_Matrix_Classes::Real start_val_default; ///< default value for starting point, default 0.
    CH_Matrix_Classes::Real bounds_minus_infinity;    ///< default value for variable lower bounds, default CB_minus_infinity, must be not greater than bounds_plus_infinity
    CH_Matrix_Classes::Real bounds_plus_infinity;    ///< default value for variable upper bounds, default CB_plus_infinity, must be not smaller than bounds_minus_infinity
    CH_Matrix_Classes::Real rhs_minus_infinity; ///< default value for right hand side lower bound, default CB_minus_infinity, must be not greater than rhs_plus_infinity
    CH_Matrix_Classes::Real rhs_plus_infinity; ///< default value for right hand side upper bound, default CB_plus_infinity, must be not smaller than rhs_minus_infinity
    CH_Matrix_Classes::Real cost_default; ///< default value for costs, default 0.
    //@}


  public:

    //-----------------------------------------------------------------------
    /** @name Constructors and initialization
     */
     //@{

     ///
    virtual ~Modification();


    /// calls clear() with these parameters and initializes CBout
    Modification(CH_Matrix_Classes::Integer var_olddim,
      CH_Matrix_Classes::Integer row_olddim,
      const CBout* cb = 0,
      int incr = -1,
      bool ensure_start_val_box_feasibility = false,
      bool ensure_bounds_consistency = true,
      CH_Matrix_Classes::Real start_val_def = 0.,
      CH_Matrix_Classes::Real bounds_minus_infty = CB_minus_infinity,
      CH_Matrix_Classes::Real bounds_plus_infty = CB_plus_infinity,
      CH_Matrix_Classes::Real rhs_minus_infty = CB_minus_infinity,
      CH_Matrix_Classes::Real rhs_plus_infty = CB_plus_infinity,
      CH_Matrix_Classes::Real cost_def = 0.);

    /** @brief resets all variables so that the object to be modified has
        starting size var_olddim (number of variables) and row_olddim
        (number of rows) and no modifications

        The actual old data is not needed at this point,
        the changes on it will be collected and excuted in the
        routines apply_to_vars and apply_to_rows

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
      CH_Matrix_Classes::Integer row_olddim,
      bool ensure_start_val_box_feasibility = false,
      bool ensure_bounds_consistency = true,
      CH_Matrix_Classes::Real start_val_def = 0.,
      CH_Matrix_Classes::Real bounds_minus_infty = CB_minus_infinity,
      CH_Matrix_Classes::Real bounds_plus_infty = CB_plus_infinity,
      CH_Matrix_Classes::Real rhs_minus_infty = CB_minus_infinity,
      CH_Matrix_Classes::Real rhs_plus_infty = CB_plus_infinity,
      CH_Matrix_Classes::Real cost_def = 0.);

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
    int add_set_lb(CH_Matrix_Classes::Integer ind, CH_Matrix_Classes::Real lb);

    /** @brief set the upper bound on variable with index @a ind to value @a ub

         If value @a ub is below CB_minus_infinity or the index is out of
         range, the return value is 1 and no changes are carried out,
         otherwise it returns 0.
    */
    int add_set_ub(CH_Matrix_Classes::Integer ind, CH_Matrix_Classes::Real ub);

    /** @brief append information on new variables at the respective ends

       @param append_dim
              number of variables (or columns of the constraint matrix) to
        be appended

       @param append_lb
              if NULL, append default values, otherwise it must point to
        a column vector of size @a append_dim that is to be appended
              to the vector of lower bounds on the variable values.
              Errors occur whenever a value is greater than CB_plus_infinity
        or the corresponding upper bound (if specified).
              Warnings are issued whenever a value equals CB_plus_infinity,
              is smaller than CB_minus_infinity or is greater than
              the corresponding start_val (if given).

       @param append_ub
              if NULL, append default values, otherwise it must point to
        a column vector of size @a append_dim that is to be appended
              to the vector of upper bounds on the variable values
              Errors occur whenever a value is smaller than CB_minus_infinity
        or the corresponding lower bound (if specified).
              Warnings are issued whenever a value equals CB_minus_infinity,
              is greater than CB_plus_infinity or is smaller than
              the corresponding start_val (if given).

       @param append_cols
              if NULL, append default values, otherwise it must point to
        a sparse matrix of size new_rowdim() times @a append_dim
              that is to be appended to the constraint matrix on the right

       @param append_start_val
              if NULL, append default values, otherwise it must point to
        a column vector of size @a append_dim that is to be appended
              to the vector of starting values for the variables. Starting
              values may violate the bounds, but this causes a warning.

       @param append_linear_costs
              if NULL, append default values, otherwise it must point to
        a column vector of size @a append_dim that is to be appended
              to the cost vector for the variables

       @return number of errors; if errors occured, none of the new changes are performed
    */
    int add_append_vars(CH_Matrix_Classes::Integer append_dim,
      const CH_Matrix_Classes::Matrix* append_lb,
      const CH_Matrix_Classes::Matrix* append_ub,
      const CH_Matrix_Classes::Sparsemat* append_cols,
      const CH_Matrix_Classes::Matrix* append_start_val,
      const CH_Matrix_Classes::Matrix* append_linear_costs);

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



    /** @brief set the lower bound on row right hand side with index @a ind to value @a rhslb

        If value @a rhslb exceeds CB_plus_infinity or the index is out of
        range, the return value is 1 and no changes are carried out,
        otherwise it returns 0.
    */
    int add_set_rhslb(CH_Matrix_Classes::Integer ind, CH_Matrix_Classes::Real rhslb);

    /** @brief set the upper bound on row right hand side with index @a ind to value @a rhsub

        If value @a rhsub is below CB_minus_infinity or the index is out of
        range, the return value is 1 and no changes are carried out,
        otherwise it returns 0.
    */
    int add_set_rhsub(CH_Matrix_Classes::Integer ind, CH_Matrix_Classes::Real rhsub);

    /** @brief append information on new rows at the respective ends

       @param append_dim
              number of rows (constraints) to be appended

       @param append_rows
              if NULL, append default values, otherwise it must point
              to a sparse matrix of size @a append_dim times new_vardim()
              that is to be appended to the constraint matrix below.

       @param append_rhslb
              if NULL, append default values, otherwise it must point to
        a column vector of size @a append_dim that is to be appended
              to the vector of right hand side lower bounds

       @param append_rhsub
              if NULL, append default values, otherwise it must point to
        a column vector of size @a append_dim that is to be appended
              to the vector of right hand side upper bounds

       @return number of errors; if errors occured, none of the new changes are performed
    */
    int add_append_rows(CH_Matrix_Classes::Integer append_dim,
      const CH_Matrix_Classes::Sparsemat* append_rows,
      const CH_Matrix_Classes::Matrix* append_rhslb,
      const CH_Matrix_Classes::Matrix* append_rhsub);

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
    int add_reassign_rows(const CH_Matrix_Classes::Indexmatrix& map_to_old);

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
    int add_delete_rows(const CH_Matrix_Classes::Indexmatrix& del_ind,
      CH_Matrix_Classes::Indexmatrix& map_to_old);

    /** @brief add the modification specified in @a m on top of
        the modifications collected so far

        For this, the old_vardim() of modification @a m must be
        identical to new_vardim() of this and
        old_rowdim() of modification @a m must be identical to
        new_rodim() of this.  The return value is the number of
        errors in this respect. If such occured, this incorporation
        is not performed.
     */
    virtual int incorporate(const Modification& m);

    //@}

    //-----------------------------------------------------------------------
    /** @name  Routines for applying the collected modifications
     */

     //@{

     /** @brief carry out the collected modifications on the given vectors

         If a specific parameter is NULL, no changes are performed on it,
         if it is not null, it must point to a column vector of length
         old_vardim(). Then the following operations will be performed
         on it in this sequence:

         1. new information is appended (or default values, if the
            corresponding append information pointer is NULL)

         2. if reassignment information is given, the values are
            mapped as given by *map_to_old_variables()

         @param[in,out] vars
            if not NULL, this points to the old vector of variable
      values. The start values will be appended to it before
      a possible index reassignment.

         @param[in,out] lb
            if not NULL, this points to the old vector of lower bounds
      on the variables.

         @param[in,out] ub
            if not NULL, this points to the old vector of upper bounds
      on the variables

         @param[in,out] cost
            if not NULL, this points to the old cost vector

         @return the number of dimension errors of non NULL inputs,
            if any, no modifications are made to any inputs.
     */
    int apply_to_vars(CH_Matrix_Classes::Matrix* vars,
      CH_Matrix_Classes::Matrix* lb,
      CH_Matrix_Classes::Matrix* ub,
      CH_Matrix_Classes::Matrix* cost) const;

    /** @brief carry out the collected modifications on the given vectors

        If a specific parameter is NULL, no changes are performed on it.
        If it is not null, @a rows must point to a sparse matrix of
        size old_rowdim() times old_vardim(), @a rhslb or @a rhsub
        must point to a column vector of length old_rowdim().
        Then the following operations will be performed on it in
        this sequence (variable/column operations only apply to @a rows) :

        1. new column information is appended to @a rows (or default values,
           if the corresponding append information pointer is NULL)

        2. if column reassignment information is given, the columns of @a rows
           are mapped as given by *map_to_old_variables()

        3. new row information is appended (or default values,
           if the corresponding append information pointer is NULL)

        4. if row reassignment information is given, the rows are
           mapped as given by *map_to_old_rows()

        @param[in,out] rows
           if not NULL, this points to the old sparse matrix of
     constraint rows.

        @param[in,out] rhslb
           if not NULL, this points to the old vector of right hand
     side lower bounds

        @param[in,out] rhsub
           if not NULL, this points to the old vector of right hand
     side upper bounds

        @return the number of dimension errors of non NULL inputs,
           if any, no modifications are made to any inputs.
    */
    int apply_to_rows(CH_Matrix_Classes::Sparsemat* rows,
      CH_Matrix_Classes::Matrix* rhslb,
      CH_Matrix_Classes::Matrix* rhsub) const;

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
    /// returns true if all entries deleted in @a oldpoint (must be a vector of length old_vardim()) are 0 and false otherwise
    bool deleted_variables_are_zero(const CH_Matrix_Classes::Matrix& oldpoint) const;
    /// returns true if all entries in newpoint (must be a vector of length new_vardim()) that correspond to new variables have value 0 and false otherwise
    bool new_variables_are_zero(const CH_Matrix_Classes::Matrix& newpoint) const;
    /// returns true if the values in newpoint (must be a vector of length new_vardim()) that correspond to old variables match the old values stored in oldpoint (must be a vector of length old_vardim()) and false otherwise
    bool mapped_variables_are_equal(const CH_Matrix_Classes::Matrix& newpoint,
      const CH_Matrix_Classes::Matrix& oldpoint) const;

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
    CH_Matrix_Classes::Integer old_rowdim() const {
      return row_olddim;
    }
    /// returns the number of rows once all stored modifications have been performed
    CH_Matrix_Classes::Integer new_rowdim() const {
      return row_newdim;
    }
    /// returns the number of rows that are appended (due to later reassignments they may no longer be located at the end)
    CH_Matrix_Classes::Integer appended_rowdim() const {
      return row_append_dim;
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
    const CH_Matrix_Classes::Indexmatrix* map_to_old_rows() const {
      return row_map_to_old;
    }

    /// returns null if there were no deletions, otherwise the Indexmatrix pointed to is a vector holding the deleted old row indices in increasing order   
    const CH_Matrix_Classes::Indexmatrix* deleted_row_indices() const {
      return row_del_ind;
    }

    /// returns null if no rows were added, otherwise the Indexmatrix pointed ato is a vector holding the new indices of the new rows in increasing order
    const CH_Matrix_Classes::Indexmatrix* new_row_indices() const {
      return row_new_ind;
    }

    /// returns null if no variable lower bounds change, otherwise it points to the map holding the (index,value) pairs with the new values 
    const std::map<CH_Matrix_Classes::Integer, CH_Matrix_Classes::Real>*
      get_var_set_lb() const {
      return var_set_lb;
    }
    /// returns null if no variable upper bounds change, otherwise it points to the map holding the (index,value) pairs with the new values 
    const std::map<CH_Matrix_Classes::Integer, CH_Matrix_Classes::Real>*
      get_var_set_ub() const {
      return var_set_ub;
    }
    /// returns null if nothing or default values have to be appended, otherwise it points to a sparse matrix whose columns need to be appended to the matrix
    const CH_Matrix_Classes::Sparsemat* get_var_append_cols() const {
      return var_append_cols;
    }
    /// returns null if nothing or default values have to be appended, otherwise it points to a matrix whose entries need to be appended to the lower bounds vector
    const CH_Matrix_Classes::Matrix* get_var_append_lb() const {
      return var_append_lb;
    }
    /// returns null if nothing or default values have to be appended, otherwise it points to a matrix whose entries need to be appended to the upper bounds vector
    const CH_Matrix_Classes::Matrix* get_var_append_ub() const {
      return var_append_ub;
    }
    /// returns null if nothing or default values have to be appended, otherwise it points to a matrix whose entries need to be appended to the starting values vector
    const CH_Matrix_Classes::Matrix* get_var_start_val() const {
      return var_start_val;
    }
    /// returns null if no right hand side lower bounds change, otherwise it points to the map holding the (index,value) pairs with the new values 
    const std::map<CH_Matrix_Classes::Integer, CH_Matrix_Classes::Real>*
      get_row_set_rhslb() const {
      return row_set_rhslb;
    }
    /// returns null if no right hand side upper bounds change, otherwise it points to the map holding the (index,value) pairs with the new values 
    const std::map<CH_Matrix_Classes::Integer, CH_Matrix_Classes::Real>*
      get_row_set_rhsub() const {
      return row_set_rhsub;
    }
    /// returns null if nothing or default values have to be appended, otherwise it points to a matrix whose entries need to be appended to the cost vector
    const CH_Matrix_Classes::Matrix* get_var_append_costs() const {
      return var_append_costs;
    }
    /// returns null if nothing or default values have to be appended, otherwise it points to a sparse matrix whose rows need to be appended to the matrix
    const CH_Matrix_Classes::Sparsemat* get_row_append_mat() const {
      return row_append_mat;
    }
    /// returns null if nothing or default values have to be appended, otherwise it points to a matrix whose entries need to be appended to the right hand side lower bounds vector
    const CH_Matrix_Classes::Matrix* get_row_append_rhslb() const {
      return row_append_rhslb;
    }
    /// returns null if nothing or default values have to be appended, otherwise it points to a matrix whose entries need to be appended to the right hand side upper bounds vector
    const CH_Matrix_Classes::Matrix* get_row_append_rhsub() const {
      return row_append_rhsub;
    }

    //@}

  };





  //@}

}

#endif

