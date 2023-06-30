@doc raw"""
    CBMatrixCBSolver(print_level::Integer = 0)

default constructor allows to set output level options from start (see also set_out())
"""
CBMatrixCBSolver(print_level::Integer = 0) = CBMatrixCBSolver(@ccall libcb.cb_matrixcbsolver_new(print_level::Cint)::Ptr{Cvoid})

@doc raw"""
    cb_clear!(self::CBMatrixCBSolver)

* @brief Clears all data structures and problem information
        but keeps ouptut settings and algorithmic parameter settings
    
"""
cb_clear!(self::CBMatrixCBSolver) = @ccall libcb.cb_matrixcbsolver_clear(self.data::Ptr{Cvoid})::Cvoid

@doc raw"""
    cb_set_defaults!(self::CBMatrixCBSolver)

* @brief Sets default values for algorithmic parameters that are not function specific (e.g., relative precision, weight and weight bounds for the augmentedproblem, etc.)
    
"""
cb_set_defaults!(self::CBMatrixCBSolver) = @ccall libcb.cb_matrixcbsolver_set_defaults(self.data::Ptr{Cvoid})::Cvoid

@doc raw"""
    cb_init_problem!(self::CBMatrixCBSolver, dim::Integer, lbounds::Union{<:CBMatrix,Nothing} = nothing, ubounds::Union{<:CBMatrix,Nothing} = nothing, startval::Union{<:CBMatrix,Nothing} = nothing, costs::Union{<:CBMatrix,Nothing} = nothing, offset::Real = 0.)

* @brief Initializes the problem by setting up the design space
        (the dimension and possibly box constraints on the variables)

    Clears all data structures and sets the dimension @ m for a new problem.
    for solving   min_{y in R^m}  f_0(y) + f_1(y) + ...
    Box constraints may be specified for y. (The functions f_i must be added
     by add_function()).

     Lower and/or upper bounds must be speicified for all variables
     or for none of them. To specify no bounds at all, give Null
     pointers. Otherwise use ConicBundle::CB_minus_infinity for
     unbounded below and ConicBundle::CB_plus_infinity for unbounded above.
     For NULL pointers, unbounded will be used as default for all
     variables. Specifying bounds selectively is also possible
     by set_lower_bound() or set_upper_bound(). For further constraints
     see append_constraints().

     @param[in] dim  (int)
         the dimension of the argument/design space/the number of Lagrange multipliers

     @param[in] lbounds  (const Matrix*)
         If NULL, all variables are considered unbounded below,
         otherwise lbounds[i] gives the minimum feasible value for variable y[i],
         use ConicBundle::CB_minus_infinity for unbounded below.

     @param[in] ubounds (const Matrix*)
         If NULL, all variables are considered unbounded above,
         otherwise ubounds[i] gives the maximum feasible value for variable y[i],
         use ConicBundle::CB_plus_infinity for unbounded above.

     @param[in] startval (const Matrix*)
        If NULL, the starting values are obtained by projecting
        zero onto the feasible set given by the lower and upper bounds
  resulting from the arguments before

     @param[in] costs (const Matrix*)
         Use this in order to specify linear costs on the variables in addition
   to the functions (may be convenient in Lagrangean relaxation for
   the right hand side of coupling contsraints); NULL is equivalent
   to costs zero.

     @param[in] offset (Real)
         Use this in order to specify linear costs on the variables in addition
   to the functions (may be convenient in Lagrangean relaxation for
   the right hand side of coupling contsraints); NULL is equivalent
   to costs zero.

     @return
        - 0 on success
        - != 0 otherwise


     
"""
cb_init_problem!(self::CBMatrixCBSolver, dim::Integer, lbounds::Union{<:CBMatrix,Nothing} = nothing, ubounds::Union{<:CBMatrix,Nothing} = nothing, startval::Union{<:CBMatrix,Nothing} = nothing, costs::Union{<:CBMatrix,Nothing} = nothing, offset::Real = 0.) = @ccall libcb.cb_matrixcbsolver_init_problem(self.data::Ptr{Cvoid}, dim::Cint, (isnothing(lbounds) ? C_NULL : lbounds.data)::Ptr{Cvoid}, (isnothing(ubounds) ? C_NULL : ubounds.data)::Ptr{Cvoid}, (isnothing(startval) ? C_NULL : startval.data)::Ptr{Cvoid}, (isnothing(costs) ? C_NULL : costs.data)::Ptr{Cvoid}, offset::Cdouble)::Cint

@doc raw"""
    cb_add_function!(self::CBMatrixCBSolver, function_::CBFunctionObject, fun_factor::Real = 1., fun_task::CBFunctionTask = cbft_objective_function, aft::Union{<:CBAffineFunctionTransformation,Nothing} = nothing, argument_list_may_change_dynamically::Bool = false)

* @brief Adds a function, typically derived from ConicBundle::FunctionOracle. If the dimension does not match the current one, specify an affine function transformation to map the current ground set to the argument of the function.

     Besides the standard ConicBundle::MatrixFunctionOracle
     the interface only accepts a few other prespecified derivations
     of the class FunctionObject that come along with the CH_Matrix_Classes interface
     (e.g. for semidefinite and second order cones). Functions not derived from
     these will fail to be added and return a value !=0.

     The @a fun_factor allows to specify a scaling factor for the function. @a fun_factor must be a strictly positive number.

     The ConicBundle::FunctionTask @a fun_task specifies whether the
     function is to be used as a an ObjectiveFunction, a
     ConstantPenaltyFunction with @a fun_factor as maximum penalty
     factor, or as an AdaptivePenaltyFunction with @a fun_factor at
     initial penalty guess that might be increased or decreased over
     time.

     The AffineFunctionTransformation @a aft may be used to modify the
     argument and give an additional affine term (linear term plus
     offset). For adding an affine term there are several other
     possibilities, e.g. in init_problem(), so there is no need to do so
     here. If, however, an existing function implementation requires only
     some subset of the variables, it is more convenient to supply
     a corresponding @a aft instead of reimplementing the function.

     @a argument_list_may_change_dynamically sets a flag on how to
     treat the function arguments when variables are added or deleted.
     If the arguments may not change, any changes in the variables are
     mapped to an adaptation of an internal
     AffineFunctionTransformation so that the function does not notice
     the changes in the variables.  If arguments may change, the
     function oracle should be a ModifiableOracleObject and react
     accordingly to the changes in its
     ModifiableOracleObject::apply_modification() routine.


    @return
      - 0 on success
      - != 0 otherwise
    
"""
cb_add_function!(self::CBMatrixCBSolver, function_::CBFunctionObject, fun_factor::Real = 1., fun_task::CBFunctionTask = cbft_objective_function, aft::Union{<:CBAffineFunctionTransformation,Nothing} = nothing, argument_list_may_change_dynamically::Bool = false) = @ccall libcb.cb_matrixcbsolver_add_function(self.data::Ptr{Cvoid}, function_.data::Ptr{Cvoid}, fun_factor::Cdouble, fun_task::Cint, (isnothing(aft) ? C_NULL : aft.data)::Ptr{Cvoid}, argument_list_may_change_dynamically::Cint)::Cint

@doc raw"""
    cb_set_lower_bound!(self::CBMatrixCBSolver, i::Integer, lb::Real)

*@brief Sets lower bound for variable i,
  use ConicBundle::CB_minus_infinity for unbounded from below.

       The algorithm may have to adapt the center point aftwards.
       In this case the old function values will be marked as outdated and
       will be recomputed at the next call to e.g. solve().

    @return
      - 0 on success
      - != 0 otherwise
    
"""
cb_set_lower_bound!(self::CBMatrixCBSolver, i::Integer, lb::Real) = @ccall libcb.cb_matrixcbsolver_set_lower_bound(self.data::Ptr{Cvoid}, i::Cint, lb::Cdouble)::Cint

@doc raw"""
    cb_set_upper_bound!(self::CBMatrixCBSolver, i::Integer, ub::Real)

*@brief Sets upper bound for variable i,
  use ConicBundle::CB_plus_infinity for unbounded from below.

       The algorithm may have to adapt the center point aftwards.
       In this case the old function values will be marked as outdated and
       will be recomputed at the next call to e.g. solve().

     @return
      - 0 on success
      - != 0 otherwise
    
"""
cb_set_upper_bound!(self::CBMatrixCBSolver, i::Integer, ub::Real) = @ccall libcb.cb_matrixcbsolver_set_upper_bound(self.data::Ptr{Cvoid}, i::Cint, ub::Cdouble)::Cint

@doc raw"""
    cb_append_variables!(self::CBMatrixCBSolver, n_append::Integer, lbounds::Union{<:CBMatrix,Nothing} = nothing, ubounds::Union{<:CBMatrix,Nothing} = nothing, constraint_columns::Union{<:CBSparsemat,Nothing} = nothing, startval::Union{<:CBMatrix,Nothing} = nothing, costs::Union{<:CBMatrix,Nothing} = nothing)

* @brief Append new variables (always in last postions in this order).

     @attention Be sure to include a desription of required changes to your
       functions via @a affected_functions_with_modifications

     @param[in] n_append  (int)
       number of variables to append (always in last position in the same order)

     @param[in] lbounds  (const Matrix*)
        If NULL, all appended variables are considered unbounded below,
        otherwise lbounds[i] gives the minimum feasible value for variable y[i],
        use ConicBundle::CB_minus_infinity for unbounded below.

     @param[in] ubounds (const Matrix*)
        If NULL, all appended variables are considered unbounded above,
        otherwise ubounds[i] gives the maximum feasible value for variable y[i],
        use ConicBundle::CB_plus_infinity for unbounded above.

     @param[in] constraint_columns (const Sparsemat*)
        This must be NULL unless append_constraints() has been used before for
  specifying linear constraints on the ground set; if there are constraints,
  NULL is interpreted as appending zero columns to the constraints,
  otherwise the the number of rows of the Sparsemant has to match the
  current number of linear constraints and the number of columns the
  must equal n_append.

     @param[in] startval (const Matrix*)
        If NULL, the starting values are obtained by projecting
        zero onto the feasible set given by the lower and upper bounds
  resulting from the arguments before

     @param[in] costs (const Matrix*)
         Use this in order to specify linear costs on the variables in addition
   to the functions (may be convenient in Lagrangean relaxation for
   the right hand side of coupling contsraints); NULL is equivalent
   to costs zero.

     @param[in,out] affected_functions_with_modifications (const FunObjModMap*)
        If NULL, default actions are performed on all functions. In
        particular, those admitting dynamic argument changes will get
        the new variables appended at the end of their argument vector
        and (eventually) their apply_modification() routines will be
        called informing them about the groundset changes; for those
        not admitting changes in their arguments, their corresponding
        (possibly newly created) affine function transformation will
        be set up to ignore the new arguments.  If !=NULL, for the
        listed functions (and their parents up to the root function)
        the default appending action is performed unless their
        FunctionObjectModification entry gives explicit modification
        instructions which are then applied instead. For all functions
        NOT listed in the map and not having modified offsprings their
        corresponding aft will be set up to ignore the new variables.

     @return
        - 0 on success
        - != 0 otherwise


    
"""
cb_append_variables!(self::CBMatrixCBSolver, n_append::Integer, lbounds::Union{<:CBMatrix,Nothing} = nothing, ubounds::Union{<:CBMatrix,Nothing} = nothing, constraint_columns::Union{<:CBSparsemat,Nothing} = nothing, startval::Union{<:CBMatrix,Nothing} = nothing, costs::Union{<:CBMatrix,Nothing} = nothing) = @ccall libcb.cb_matrixcbsolver_append_variables(self.data::Ptr{Cvoid}, n_append::Cint, (isnothing(lbounds) ? C_NULL : lbounds.data)::Ptr{Cvoid}, (isnothing(ubounds) ? C_NULL : ubounds.data)::Ptr{Cvoid}, (isnothing(constraint_columns) ? C_NULL : constraint_columns.data)::Ptr{Cvoid}, (isnothing(startval) ? C_NULL : startval.data)::Ptr{Cvoid}, (isnothing(costs) ? C_NULL : costs.data)::Ptr{Cvoid})::Cint

@doc raw"""
    cb_delete_variables!(self::CBMatrixCBSolver, delete_indices::CBIndexmatrix, map_to_old::CBIndexmatrix)

* @brief Deletes variables corresponding to the specified indices.

        The indices of the remaining variables are reassigned so that they
        are consecutive again, the routine returns in @a map_to_old
        a vector giving for each new index of these remaining variables
        the old coordinate.

     @attention Be sure to include a desription of required changes to your
       functions via @a affected_functions_with_modifications

     @param[in] delete_indices  (const Indexmatrix&)
        the entries delete_indices[i] specify the indices of the variables
        to be deleted

     @param[out] map_to_old  (Indexmatrix&)
        after the call, element map_to_old[i] gives the old index (before the call)
        of the variable that now has index position i.

     @param[in] affected_functions_with_modifications (const FunObjModMap*)
        If NULL, default actions are performed on all functions. In
        particular, for those admitting dynamic argument changes all
  those variables will be deleted whose row in a corresponding
  updated affine function transformation (so after deletion of
  the columns of the incoming variables) corresponds to the zero
  map (i.e., offset and matrix row are both zero), furthermore
  identity transformations will be preserved. For those
        not admitting changes in their arguments, their corresponding
        (possibly newly created) affine function transformation will
        only get the columns deleted, but there will be no row deleltions.
  If !=NULL, for the listed functions (and their parents up to the
  root function) the default deletion action is performed unless their
        FunctionObjectModification entry gives explicit modification
        instructions which are then applied instead. For all functions
        NOT listed in the map and not having modified offsprings their
        corresponding aft will be set up to keep the arguments unchanged.

     @return
        - 0 on success
        - != 0 otherwise

    
"""
cb_delete_variables!(self::CBMatrixCBSolver, delete_indices::CBIndexmatrix, map_to_old::CBIndexmatrix) = @ccall libcb.cb_matrixcbsolver_delete_variables(self.data::Ptr{Cvoid}, delete_indices.data::Ptr{Cvoid}, map_to_old.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_reassign_variables!(self::CBMatrixCBSolver, assign_new_from_old::CBIndexmatrix)

* @brief Reassigns variables to new index positions by mapping to position @a i
        the variable that previously had index @a assign_new_from_old[i].

        Old variables, that are not mapped to any position will be deleted.
        It is not allowed to generate several copies of old variables.

     @attention Be sure to include a desription of required changes to your
       functions via @a affected_functions_with_modifications

     @param[in] assign_new_from_old  (const IVector&)
        entry assign_new_from_old[i] specifies
        the old index of the variable, that has to be copied to index position i.

     @param[in] affected_functions_with_modifications (const FunObjModMap*)
        If NULL, default actions are performed on all functions. In
        particular, for those admitting dynamic argument changes all
  those variables will be deleted whose row in a corresponding
  updated affine function transformation (so after mapping
  the columns of the incoming variables) correspond to the zero
  map (i.e., offset and matrix row are both zero); furthermore,
  if the transformation was the identity to start with, this will
  be preserved by mapping the arguments in the same way. For those
        not admitting changes in their arguments, their corresponding
        (possibly newly created) affine function transformation will
        only get the columns mapped, but there will be no row deleltions.
  If !=NULL, for the listed functions (and their parents up to the
  root function) the default deletion action is performed unless their
        FunctionObjectModification entry gives explicit modification
        instructions which are then applied instead. For all functions
        NOT listed in the map and not having modified offsprings their
        corresponding aft will be set up to keep the arguments unchanged.

     @return
        - 0 on success
        - != 0 otherwise

    
"""
cb_reassign_variables!(self::CBMatrixCBSolver, assign_new_from_old::CBIndexmatrix) = @ccall libcb.cb_matrixcbsolver_reassign_variables(self.data::Ptr{Cvoid}, assign_new_from_old.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_append_constraints!(self::CBMatrixCBSolver, append_n_rows::Integer, append_rows::Union{<:CBSparsemat,Nothing} = nothing, append_rhslb::Union{<:CBMatrix,Nothing} = nothing, append_rhsub::Union{<:CBMatrix,Nothing} = nothing)

* @brief append \f$rhslb\le Ay \le rhsub\f$ as linear constraints on the groundset variables \f$y\f$. \f$A\f$ has   @a append_n_rows new rows with coefficients given in  append_rows (if NULL, use default value),  lower bounds append_rhslb (if NULL use default) and upper bounds append_rhsub (if NULL use default)

     @param[in] append_n_rows  (Integer)
        (nonnegative) number of rows to be appended
  as linear constraints on the ground set

     @param[in] append_rows  (Sparsemat*)
        describes the coefficients of the linear constraints; the
  number of rows must match append_n_rows, the number of columns
  must match the current dimension of the groundset;
  if NULL, all coefficients are considered zero.

     @param[in] append_rhslb  (Matrix*)
        specifies lower bounds on the values of the constraints;
        the  number of rows must match append_n_rows;
  if NULL, all coefficients are considered CB_minus_infinity.

     @param[in] append_rhsub  (Matrix*)
        specifies upper bound on the values of the constraints;
        the  number of rows must match append_n_rows;
  if NULL, all coefficients are considered CB_plus_infinity.

     @return
        - 0 on success
        - != 0 otherwise

     
"""
cb_append_constraints!(self::CBMatrixCBSolver, append_n_rows::Integer, append_rows::Union{<:CBSparsemat,Nothing} = nothing, append_rhslb::Union{<:CBMatrix,Nothing} = nothing, append_rhsub::Union{<:CBMatrix,Nothing} = nothing) = @ccall libcb.cb_matrixcbsolver_append_constraints(self.data::Ptr{Cvoid}, append_n_rows::Cint, (isnothing(append_rows) ? C_NULL : append_rows.data)::Ptr{Cvoid}, (isnothing(append_rhslb) ? C_NULL : append_rhslb.data)::Ptr{Cvoid}, (isnothing(append_rhsub) ? C_NULL : append_rhsub.data)::Ptr{Cvoid})::Cint

@doc raw"""
    cb_solve!(self::CBMatrixCBSolver, maxsteps::Integer = 0, stop_at_descent_steps::Bool = false)

* @brief solves or does a prescribed number of iterations

  Bundle methods solve a problem by a sequence of so called
        descent steps that actually bring progress by moving from the
        current "center point" to a new center with better objective.
        A descent step may consist of several function evaluations
        (null steps), that lead to no immediate progress but mainly
  improve a cutting model of the objective function close to
        the current center point.  A minimizer to the model is
        accepted as descent step if the function value at this point
        satisfies a sufficient decrease criterion in comparison to the
        decrease predicted by the model. Having found a descent step,
        the next center is automatically shifted to this successful
        candidate.  Termination criteria may stop the process of
        seeking for a descent step, in which case the current center
        is kept and the routine termination_code() returns the
        termination code.

        Restarting, after each descent step, the bundle method from scratch
        with the new center as starting point does not endanger convergence.
        Therefore, a descent step is the smallest unit, after which
        user interaction can take place safely. To allow this there
  is a flag stop_at_descent_steps that will cause the code to
  return after the next descent step.

  If you know what your are doing, you may also use the input
        parameter maxsteps to force the algorithm to return after
        at most maxsteps null steps. Calling solve again
        without any intermediate problem configurations will then
        simply continue the process where it stopped and convergence
        is save. During null steps one may not decrease the weight
        or delete nonzero variables of the center or the current candidate!

        In a Lagrangean relaxation cutting plane approach one may want
        to separate and enlarge the dimension after a certain number
        of null steps. In this case the code will try to preserve the model,
        given appropriate subgradient extension routines have been
        provided. If the model cannot be extended, it has to be
        discarded (if subgradient extension is not successful
        this is done automatically), and the algorithm will be restarted
        from the current center point.

      @param[in] maxsteps (int)
          if maxsteps>0 the code returns after at most so many null steps

      @param[in] stop_at_descent_steps (int)
          if true the code also returns whenever a descent step occured

      @return
        - 0 on success
        - != 0 otherwise

    
"""
cb_solve!(self::CBMatrixCBSolver, maxsteps::Integer = 0, stop_at_descent_steps::Bool = false) = @ccall libcb.cb_matrixcbsolver_solve(self.data::Ptr{Cvoid}, maxsteps::Cint, stop_at_descent_steps::Cint)::Cint

@doc raw"""
    cb_termination_code(self::CBMatrixCBSolver)

* @brief Returns the termination code of the bundle algorithm for the latest descent step

      For resetting all counters relevant for termination see clear_fail_counts() .

      @return
      -  0  :    Not terminated.
             (Continue with the next solve())
      -  1  :    Relative precision criterion satisfied. (See set_term_relprec())
      -  2  :    Timelimit exceeded.
             (Currently the C interface does not offer a timelimit.)
      -  4  :    Maximum number of function reevaluations exceeded.
             (Indicates that there is a problem with one of the function
             oracles that seems to deliver no valid upper bounds on the true
             function value for descent steps)
      -  8  :    Maximum number of quadratic subproblem failures exceeded.
             (Indicates that the numerical limits of the inner quadratic
             programming solver are reached, no further progress expected)
      - 16  :    maximum number of model evaluation failures exceeded
             (Indicates that the numerical limits of the setup of the
             subproblem are reached, no further progress expected)
      - 32  :    maximum number of failures to increase the augmented model value exceeded
             (Indicates that the numerical limits  of the interplay between
             subproblem and quadratic programming solver are reached,
             no further progress expected)
       - 64  :   maximum number of oracle calls (function evaluations) exceeded,
                 see set_eval_limit()
       - 128  :   maximum number of oracle failures exceeded.
             This refers to function evaluations that terminate with insufficient
       precision but still provide a new approximate subgradient. A failure typically
             indicates numerical difficulties with the precision requirements.
             (Currently the interface does not allow to manipulate the limit, it is set to 10)


    
"""
cb_termination_code(self::CBMatrixCBSolver) = @ccall libcb.cb_matrixcbsolver_termination_code(self.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_print_termination_code(self::CBMatrixCBSolver)

* @brief Outputs a text version of termination code, see termination_code().

      @return
        - 0 on success
        - != 0 otherwise

    
"""
cb_print_termination_code(self::CBMatrixCBSolver) = @ccall libcb.cb_matrixcbsolver_print_termination_code(self.data::Ptr{Cvoid})::Cvoid

@doc raw"""
    cb_get_objval(self::CBMatrixCBSolver)

* @brief Returns the objective value resulting from last descent
        step (initially undefined). If no problem modification routines
        were called since then, it is the objective value at the point
        returned by get_center().
    
"""
cb_get_objval(self::CBMatrixCBSolver) = @ccall libcb.cb_matrixcbsolver_get_objval(self.data::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_get_center(self::CBMatrixCBSolver, center::CBMatrix)

* @brief Returns the next center point that was produced by the latest call
  to solve (in some problem modification routines the
  center point may be updated immediately, in others the center point
  will be corrected automatically directly before starting
  the next descent step and its values may be infeasible till then).

      @return
        - 0 on success
        - != 0 otherwise
    
"""
cb_get_center(self::CBMatrixCBSolver, center::CBMatrix) = @ccall libcb.cb_matrixcbsolver_get_center(self.data::Ptr{Cvoid}, center.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_get_sgnorm(self::CBMatrixCBSolver)

* @brief Returns Euclidean norm of the latest aggregate subgradient.
     
"""
cb_get_sgnorm(self::CBMatrixCBSolver) = @ccall libcb.cb_matrixcbsolver_get_sgnorm(self.data::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_get_subgradient(self::CBMatrixCBSolver, subgradient::CBMatrix)

* @brief Returns the latest aggregate subgradient (of the entire problem with groundset as provided by the solver)

    @return
      - 0 on success
      - != 0 otherwise

    
"""
cb_get_subgradient(self::CBMatrixCBSolver, subgradient::CBMatrix) = @ccall libcb.cb_matrixcbsolver_get_subgradient(self.data::Ptr{Cvoid}, subgradient.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_get_cutval(self::CBMatrixCBSolver)

* @brief Returns the cutting model value resulting from last call to
        solve() (initially undefined).
    
"""
cb_get_cutval(self::CBMatrixCBSolver) = @ccall libcb.cb_matrixcbsolver_get_cutval(self.data::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_get_candidate_value(self::CBMatrixCBSolver)

* @brief Returns the objective value computed in the last step of solve(),
        independent of whether this was a descent step or a null step (initially undefined).

        If no problem modification routines were called since then, it is the
        objective value at the point returned by get_candidate(). If this
        last evaluation led to a descent step, then it is the same value as
        in get_objval().
    
"""
cb_get_candidate_value(self::CBMatrixCBSolver) = @ccall libcb.cb_matrixcbsolver_get_candidate_value(self.data::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_get_candidate(self::CBMatrixCBSolver, center::CBMatrix)

* @brief Returns the last point, the "candidate", at which the function
        was evaluated in solve().

        If this evaluation lead to a descent step, it is the same point as
  in get_center().

      @return
        - 0 on success
        - != 0 otherwise
    
"""
cb_get_candidate(self::CBMatrixCBSolver, center::CBMatrix) = @ccall libcb.cb_matrixcbsolver_get_candidate(self.data::Ptr{Cvoid}, center.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_set_term_relprec!(self::CBMatrixCBSolver, term_relprec::Real)

* @brief Sets the relative precision requirements for successful termination
               (default 1e-5).

     @param[in] term_relprec (double)
       The algorithm stops with termination code 1, if predicted progress for
       the next step is less than term_relprec times
       absolute function value plus one.

    @return
      - 0 on success
      - != 0 otherwise

    
"""
cb_set_term_relprec!(self::CBMatrixCBSolver, term_relprec::Real) = @ccall libcb.cb_matrixcbsolver_set_term_relprec(self.data::Ptr{Cvoid}, term_relprec::Cdouble)::Cint

@doc raw"""
    cb_set_new_center_point!(self::CBMatrixCBSolver, center_point::CBMatrix)

* @brief Set the starting point/center that will be used in the
        next call to  solve(). Each call
        to this routine causes an immediate evaluation of all oracles.

     @return
        - 0 on success
        - != 0 otherwise
    
"""
cb_set_new_center_point!(self::CBMatrixCBSolver, center_point::CBMatrix) = @ccall libcb.cb_matrixcbsolver_set_new_center_point(self.data::Ptr{Cvoid}, center_point.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_get_function_status(self::CBMatrixCBSolver, function_::CBFunctionObject)

* @brief Returns the return value of the latest evaluation call
        to this @a function.
    
"""
cb_get_function_status(self::CBMatrixCBSolver, function_::CBFunctionObject) = @ccall libcb.cb_matrixcbsolver_get_function_status(self.data::Ptr{Cvoid}, function_.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_get_approximate_slacks(self::CBMatrixCBSolver, param0::CBMatrix)

* @brief Returns the multipliers for the box constraints on the design variables;
       in Lagrangean relaxation they may be interpreted as primal slacks
 for inequality constraints.
    @return
       - 0 on success
       - != 0 otherwise
   
"""
cb_get_approximate_slacks(self::CBMatrixCBSolver, param0::CBMatrix) = @ccall libcb.cb_matrixcbsolver_get_approximate_slacks(self.data::Ptr{Cvoid}, param0.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_get_approximate_primal(self::CBMatrixCBSolver, function_::CBFunctionObject)

* @brief returns the current approximate primal solution corresponding
        to the aggregate subgradient of the specified @a function.

        PrimalData solutions must have been supplied in all previous
        calls to evaluate; In this case it returns the current approximate
        primal solution aggregated alongside with the aggregate subgradient.
        A primal solution may not be available after addition of constraints,
        if extension of the aggregate subgradient to the new coordinates failed.
  If no primal data is availalbe, the function returns NULL.

     @return
        - pointer to the primal data of the aggregate of this function object
        - 0 if no primal is available
    
"""
cb_get_approximate_primal(self::CBMatrixCBSolver, function_::CBFunctionObject) = CBPrimalData(@ccall libcb.cb_matrixcbsolver_get_approximate_primal(self.data::Ptr{Cvoid}, function_.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_get_center_primal(self::CBMatrixCBSolver, function_::CBFunctionObject)

* @brief Returns the primal solution corresponding to the best epsilon
  subgradient returned in the evaluation of the specified @a function
        at the current center point. If no primal data is availalbe,
  the function returns NULL.

     @return
        - pointer to the primal data of the minorant returned on evaluation
    of this function object at the current center
        - 0 if no primal is available
    
"""
cb_get_center_primal(self::CBMatrixCBSolver, function_::CBFunctionObject) = CBPrimalData(@ccall libcb.cb_matrixcbsolver_get_center_primal(self.data::Ptr{Cvoid}, function_.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_get_candidate_primal(self::CBMatrixCBSolver, function_::CBFunctionObject)

* @brief Returns the primal solution corresponding to the best epsilon
  subgradient returned in the evaluation of the specified @a function
        at the point get_candidate. If no primal data is availalbe,
  the function returns NULL.

     @return
        - pointer to the primal data of the minorant returned on evaluation
    of this function object at the current candidate
        - 0 if no primal is available
    
"""
cb_get_candidate_primal(self::CBMatrixCBSolver, function_::CBFunctionObject) = CBPrimalData(@ccall libcb.cb_matrixcbsolver_get_candidate_primal(self.data::Ptr{Cvoid}, function_.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_set_sumbundle!(self::CBMatrixCBSolver, use_sumbundle::Bool, n_local_models::Integer = -1, bundle_parameters::Union{<:CBBundleParameters,Nothing} = nothing, strategy::Integer = 1)

* @brief Starts/ends the use of a common SumBundle of the given bundle_size
        with a heuristic rule for selecting up to n_local_models in each bundle iteration

        If the function is the sum of many functions, having a local model for every one
        of them may result in a huge quadratic subproblem. It may then be better to
  form a common model of most of the functions, where a heuristic dynamically
  selects a few of the functions, for which a local model seems worth while.
  Whether such a common model should be used, how many subgradients it should contain,
  and how many local models are to be selected at most are the parameters
  set here.

  Setting these parameters only has an effect if bundle models of functions
  are present. If further functions are added later, the call should be repeated.

  This interface provides a simpler access to the SumBundle features by using
  some default parameter choices that could be set separately in
  set_bundle_parameters() and set_sumbundle_parameters() in a refined way.

     @param[in] use_sumbundle (bool)
        use value true to switch the sumbundle on, use value false to switch it off

     @param[in] n_local_models (int)
         upper bound on the number of local models to be used on top of the sumbundle's model,
   negative values correspond to no upper bound and all functions may have local models

     @param[in] bundle_parameters (const BundleParameters*)
         the maximum number of subgradients to be used in forming the SumBundle model,
   values <=1 are set to 2;

     @param[in] strategy (int)
         this is currently in experimental stage and allows to choose among some internal
         sumbundle strategies (currently 0,1,2,11 are available)

     @return
        - 0 on success
        - != 0 otherwise

    
"""
cb_set_sumbundle!(self::CBMatrixCBSolver, use_sumbundle::Bool, n_local_models::Integer = -1, bundle_parameters::Union{<:CBBundleParameters,Nothing} = nothing, strategy::Integer = 1) = @ccall libcb.cb_matrixcbsolver_set_sumbundle(self.data::Ptr{Cvoid}, use_sumbundle::Cint, n_local_models::Cint, (isnothing(bundle_parameters) ? C_NULL : bundle_parameters.data)::Ptr{Cvoid}, strategy::Cint)::Cint

@doc raw"""
    cb_set_max_modelsize!(self::CBMatrixCBSolver, max_modelsize::Integer, function_::Union{<:CBFunctionObject,Nothing} = nothing)

* @brief Sets the maximum number of subgradients used in forming the
        cutting model of the specified @a function

        Quite often a very small model, e.g., 2, yields very fast iterations
        and good progress in time (sometimes at the cost of more evaluations).
        By limited numerical experience, a significant reduction in the number of
        evaluations can  only be expected if the bundle is large enough to
        wrap the function rather tightly. Quite frequently, unfortunately,
        this entails that solving the quadratic subproblems
        is more expensive than function evaluation.

        The meaning of this routine may differ from standard for
  predefined special functions with special bundle types.

     @param[in] max_modelsize (int)
         maximum number of subgradients to be used in forming the cutting model

     @param[in] function
       if the aggregate subgradient of a particular function is desired,
       provide the pointer here, otherwise this referrs to the root function
       (if there is only one function to be optimized over, this is this single
       function, otherwise it is the sum of functions)

     @return
        - 0 on success
        - != 0 otherwise

    
"""
cb_set_max_modelsize!(self::CBMatrixCBSolver, max_modelsize::Integer, function_::Union{<:CBFunctionObject,Nothing} = nothing) = @ccall libcb.cb_matrixcbsolver_set_max_modelsize(self.data::Ptr{Cvoid}, max_modelsize::Cint, (isnothing(function_) ? C_NULL : function_.data)::Ptr{Cvoid})::Cint

@doc raw"""
    cb_set_max_bundlesize!(self::CBMatrixCBSolver, max_bundlesize::Integer, function_::Union{<:CBFunctionObject,Nothing} = nothing)

* @brief Sets the maximum number of subgradients stored for use in
  forming the model or determining scaling information, it must be as
  least as large as max_modelsize (and is increased to this if not)

  The meaning of this routine may differ from standard for
  predefined special functions with special bundle types.


     @param[in] max_bundlesize (int)
       maximum number of subgradients stored for use in forming the model

     @param[in] function
       if the aggregate subgradient of a particular function is desired,
       provide the pointer here, otherwise this referrs to the root function
       (if there is only one function to be optimized over, this is this single
       function, otherwise it is the sum of functions)



     @return
       - 0 on success
       - != 0 otherwise

    
"""
cb_set_max_bundlesize!(self::CBMatrixCBSolver, max_bundlesize::Integer, function_::Union{<:CBFunctionObject,Nothing} = nothing) = @ccall libcb.cb_matrixcbsolver_set_max_bundlesize(self.data::Ptr{Cvoid}, max_bundlesize::Cint, (isnothing(function_) ? C_NULL : function_.data)::Ptr{Cvoid})::Cint

@doc raw"""
    cb_set_bundle_parameters!(self::CBMatrixCBSolver, params::CBBundleParameters, function_::Union{<:CBFunctionObject,Nothing} = nothing)

*@brief Sets the maximum bundlesize and the maximum number of new subgradients
        added in a bundle update of the cutting model for the specified @a function.
        The meaning of this routine may differ from standard for
  predefined special functions with special bundle types.

     @param[in] params (const BundleParameters&)
       some update parameters for the cutting model, see e.g. ConicBundle::BundleParameters

     @param[in] function
       if the aggregate subgradient of a particular function is desired,
       provide the pointer here, otherwise this referrs to the root function
       (if there is only one function to be optimized over, this is this single
       function, otherwise it is the sum of functions)


     @return
       - 0 on success
       - != 0 otherwise

    
"""
cb_set_bundle_parameters!(self::CBMatrixCBSolver, params::CBBundleParameters, function_::Union{<:CBFunctionObject,Nothing} = nothing) = @ccall libcb.cb_matrixcbsolver_set_bundle_parameters(self.data::Ptr{Cvoid}, params.data::Ptr{Cvoid}, (isnothing(function_) ? C_NULL : function_.data)::Ptr{Cvoid})::Cint

@doc raw"""
    cb_get_bundle_parameters(self::CBMatrixCBSolver, function_::Union{<:CBFunctionObject,Nothing} = nothing)

* @brief Retrieves current bundle parameters (not the actual size in use!)
       as set for the cutting model of the specified @a function.

 This may differ for predefined special
 functions with derived BundleParameter classes.

 If the code is asked to optimize over the sum of several functions,
 it usually does this with a separate model for each function. If there
 are too many function for this, it may be worth to consider using
 the SumBundle features. For this see also set_sumbundle_parameters().
 If the root function is a sum of functions, passing a
       SumModelParametersObject here allows to specify how many local models
 should be kept by SumModelParametersObject::set_max_local_models()
       and how these should be selected. A possible implementation for this
 is given in SumModelParameters.

    @param[in] function
      if the aggregate subgradient of a particular function is desired,
      provide the pointer here, otherwise this referrs to the root function
      (if there is only one function to be optimized over, this is this single
      function, otherwise it is the sum of functions)

    @return
      - 0 on success
      - != 0 otherwise

   
"""
cb_get_bundle_parameters(self::CBMatrixCBSolver, function_::Union{<:CBFunctionObject,Nothing} = nothing) = CBBundleParameters(@ccall libcb.cb_matrixcbsolver_get_bundle_parameters(self.data::Ptr{Cvoid}, (isnothing(function_) ? C_NULL : function_.data)::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_set_sumbundle_parameters!(self::CBMatrixCBSolver, params::CBSumBundleParametersObject, function_::Union{<:CBFunctionObject,Nothing} = nothing)

*@brief Specifies the behavior of the model (of the specified function)
        concerning requests to join or start a SumBundle that subsumes several
  models instead of providing a separate model for each funciton.

        The abstract interface for these Parameters is specified in
        SumBundleParametersObject, a concrete implementation is
        SumBundleParameters. Besides the usual BundleParameters
        the new main parameter is specified in
  SumBundleParametersObject::set_acceptable_mode(), see there.

     @param[in] params (const BundleParameters&)
       some update parameters for the cutting model, see e.g. ConicBundle::BundleParameters

     @param[in] function
       if the aggregate subgradient of a particular function is desired,
       provide the pointer here, otherwise this referrs to the root function
       (if there is only one function to be optimized over, this is this single
       function, otherwise it is the sum of functions)


     @return
       - 0 on success
       - != 0 otherwise

    
"""
cb_set_sumbundle_parameters!(self::CBMatrixCBSolver, params::CBSumBundleParametersObject, function_::Union{<:CBFunctionObject,Nothing} = nothing) = @ccall libcb.cb_matrixcbsolver_set_sumbundle_parameters(self.data::Ptr{Cvoid}, params.data::Ptr{Cvoid}, (isnothing(function_) ? C_NULL : function_.data)::Ptr{Cvoid})::Cint

@doc raw"""
    cb_get_bundle_data(self::CBMatrixCBSolver, function_::Union{<:CBFunctionObject,Nothing} = nothing)

* @brief Returns all current bundle data of the cutting
        model of the specified @a function.

  This may differ for predefined special
  functions with derived classes.

     @param[in] function
       if the aggregate subgradient of a particular function is desired,
       provide the pointer here, otherwise this referrs to the root function
       (if there is only one function to be optimized over, this is this single
       function, otherwise it is the sum of functions)

     @return
       - 0 on success
       - != 0 otherwise

    
"""
cb_get_bundle_data(self::CBMatrixCBSolver, function_::Union{<:CBFunctionObject,Nothing} = nothing) = CBBundleData(@ccall libcb.cb_matrixcbsolver_get_bundle_data(self.data::Ptr{Cvoid}, (isnothing(function_) ? C_NULL : function_.data)::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_reinit_function_model!(self::CBMatrixCBSolver, function_::Union{<:CBFunctionObject,Nothing} = nothing)

* @brief Clears cutting model, subgradients and stored function values
       for the specified @a function (but only for the given one, not recursively)

       There should be no need to call this if the modification
       routines of this interface were used correctly. If, however,
       the oracle is modified by other means outside this interface,
       this has to be called whenever the specified function was
       modified so that the old subgradients and/or primal generators
       are no longer valid.

     @param[in] function
       if the aggregate subgradient of a particular function is desired,
       provide the pointer here, otherwise this referrs to the root function
       (if there is only one function to be optimized over, this is this single
       function, otherwise it is the sum of functions)

      @return
        - 0 on success
        - != 0 otherwise

    
"""
cb_reinit_function_model!(self::CBMatrixCBSolver, function_::Union{<:CBFunctionObject,Nothing} = nothing) = @ccall libcb.cb_matrixcbsolver_reinit_function_model(self.data::Ptr{Cvoid}, (isnothing(function_) ? C_NULL : function_.data)::Ptr{Cvoid})::Cint

@doc raw"""
    cb_clear_aggregates!(self::CBMatrixCBSolver, function_::Union{<:CBFunctionObject,Nothing} = nothing)

* @brief Clears the aggregate parts of the cutting model of this @a function
       (but only for the given one, not recursively)

       There should be no need to call this if the modification
       routines of this interface were used correctly. If, however,
       the oracle is modified by other means outside this interface,
       this has to be called whenever the specified function was
       modified so that the old aggregate subgradients and/or primal
       generators are no longer valid.

      @param[in] function (const FunctionObject&)
        the function added in add_function()

      @return
        - 0 on success
        - != 0 otherwise

    
"""
cb_clear_aggregates!(self::CBMatrixCBSolver, function_::Union{<:CBFunctionObject,Nothing} = nothing) = @ccall libcb.cb_matrixcbsolver_clear_aggregates(self.data::Ptr{Cvoid}, (isnothing(function_) ? C_NULL : function_.data)::Ptr{Cvoid})::Cint

@doc raw"""
    cb_call_primal_extender!(self::CBMatrixCBSolver, function_::CBFunctionObject, primal_extender::CBPrimalExtender)

* @brief Asks @a function to call @a primal_extender for each of its primal objects (see
     also FunctionOracle::evaluate() )

     If the function is the Lagrangian dual of a primal problem and primal_data
     returned previous calls to the oracle has now to be updated due to changes
     in the primal problem -- e.g., this may happen in column generation -- the
     call causes updates of all internally stored primal_data objects by calling
     PrimalExtender::extend on each of these.

      @param[in] function (const FunctionObject&)
        the function added in add_function()

      @param[in] primal_extender (PrimalExtender&)
        the object holding the extension function for primal_data

      @return
        - 0 on success
        - 1 if for this function it is not possible to use a primal_extender
        - 2 if the primal_extender would be applicable but there is no primal_data
    
"""
cb_call_primal_extender!(self::CBMatrixCBSolver, function_::CBFunctionObject, primal_extender::CBPrimalExtender) = @ccall libcb.cb_matrixcbsolver_call_primal_extender(self.data::Ptr{Cvoid}, function_.data::Ptr{Cvoid}, primal_extender.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_get_last_weight(self::CBMatrixCBSolver)

* @brief Returns the current weight for the quadratic term in the augmented subproblem
    (may be interpreted as 1./step_size or 1./trustregion-radius).
    
"""
cb_get_last_weight(self::CBMatrixCBSolver) = @ccall libcb.cb_matrixcbsolver_get_last_weight(self.data::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_get_next_weight(self::CBMatrixCBSolver)

* @brief Returns the next weight 	for the quadratic term in the augmented subproblem
        suggested by the internal weight updating heuristic
    
"""
cb_get_next_weight(self::CBMatrixCBSolver) = @ccall libcb.cb_matrixcbsolver_get_next_weight(self.data::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_set_next_weight!(self::CBMatrixCBSolver, weight::Real)

* @brief Sets the  weight (>0) to be used in the quadratic term
        of the next augmented subproblem
        (may be interpreted as 1./step_size or 1./trustregion-radius).

        Independent of whether the weight violates current min- and max-bounds
        set in set_min_weight() and set_max_weight(), the next model will
        be computed for this value. Thereafter, however, it will be updated as
        usual; in particular, it may be truncated by min and max bounds
        immediately after the first subproblem.

        In order to guarantee a constant weight (e.g. 1 is frequently a reasonable
        choice if the automatic default heuristic performs poorly), set the min and max
         bounds to the same value, too.

      @param[in] weight (double)

      @return
        - 0 on success
        - != 0 otherwise
    
"""
cb_set_next_weight!(self::CBMatrixCBSolver, weight::Real) = @ccall libcb.cb_matrixcbsolver_set_next_weight(self.data::Ptr{Cvoid}, weight::Cdouble)::Cint

@doc raw"""
    cb_set_min_weight!(self::CBMatrixCBSolver, min_weight::Real)

* @brief Sets a lower bound on the  weight for the quadratic term of the
        augmented subproblem.

        Nonpositive values indicate no bound.
        The new value shows its effect only at first dynamic change of
        the weight.

     @param[in] min_weight (double)

     @return
        - 0 on success
        - != 0 otherwise

    
"""
cb_set_min_weight!(self::CBMatrixCBSolver, min_weight::Real) = @ccall libcb.cb_matrixcbsolver_set_min_weight(self.data::Ptr{Cvoid}, min_weight::Cdouble)::Cint

@doc raw"""
    cb_set_max_weight!(self::CBMatrixCBSolver, max_weight::Real)

* @brief Sets an upper bound on the  weight for the quadratic term of the
        augmented subproblem.

        Nonpositive values indicate no bound.
        The new value shows its effect only at first dynamic change of
        the weight.

      @param[in] max_weight (double)

      @return
        - 0 on success
        - != 0 otherwise
    
"""
cb_set_max_weight!(self::CBMatrixCBSolver, max_weight::Real) = @ccall libcb.cb_matrixcbsolver_set_max_weight(self.data::Ptr{Cvoid}, max_weight::Cdouble)::Cint

@doc raw"""
    cb_set_weight_update!(self::CBMatrixCBSolver, bw::Union{<:CBBundleWeight,Nothing})

* @brief Replaces the internal update routine for choosing the weight used in the proximal term; input NULL reinstalls the default routine.

        The BundleWeight class instance pointed to will be deleted on
        construction, i.e., ownership is passe over to the solver.

      @param[in] bw
        replace internal update routine by bw, value 0 reinstalls the default routine

      @return
        - 0 on success
        - != 0 otherwise
    
"""
cb_set_weight_update!(self::CBMatrixCBSolver, bw::Union{<:CBBundleWeight,Nothing}) = @ccall libcb.cb_matrixcbsolver_set_weight_update(self.data::Ptr{Cvoid}, (isnothing(bw) ? C_NULL : bw.data)::Ptr{Cvoid})::Cint

@doc raw"""
    cb_adjust_multiplier!(self::CBMatrixCBSolver)

* @brief Adjusts on all conic functions the penalty parameter for
  conic violations to twice the trace of the primal approximation.

        This routine is only needed for conic function objects such
        as the nonnegative cone, the second order cone and
        the semidefinite cone if no good upper bound on the trace of
        feasible points is known and has to be determined automatically.

        If after some time, the trace values settle, the upper bounds
        on the trace may be way to high and can then be reset with this
        call.

      @return
        - 0 on success
        - != 0 otherwise
    
"""
cb_adjust_multiplier!(self::CBMatrixCBSolver) = @ccall libcb.cb_matrixcbsolver_adjust_multiplier(self.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_set_variable_metric!(self::CBMatrixCBSolver, do_variable_metric::Integer)

* @brief Use a variable metric heuristic or switch off general metrics alltogether.
       (variable metric resets the quadratic term e.g. to some diagonal matrix,
       switching it off resets the quadratic term to the identity times the weight)

     @param[in] do_variable_metric (int)
        - 0 switch off the scaling heuristic
        - 1 use a diagonal scaling heuristic
        - 2 use a diagonal scaling heuristic combined with one for the bounds
        - 3 use a low rank scaling heuristic
        - 4 use a low rank scaling heuristic combined with a diagonal term
        - 5 use a dense scaling heuristic

      @return
        - 0 on success
        - != 0 otherwise
    
"""
cb_set_variable_metric!(self::CBMatrixCBSolver, do_variable_metric::Integer) = @ccall libcb.cb_matrixcbsolver_set_variable_metric(self.data::Ptr{Cvoid}, do_variable_metric::Cint)::Cint

@doc raw"""
    cb_set_prox!(self::CBMatrixCBSolver, proxp::Union{<:CBBundleProxObject,Nothing})

* @brief For variable metric install the BundleProxObject pointed to; the object is passed to the solver who will delete it on termination or when replaced

      @param[in] proxp (BundleProxObject*)
          replace the current BundleProxObject by this object on the heap;
    NULL is allowed and results in the default choice;
    the object pointed to will be deleted by the solver

      @return
        - 0 on success
        - != 0 otherwise
    
"""
cb_set_prox!(self::CBMatrixCBSolver, proxp::Union{<:CBBundleProxObject,Nothing}) = @ccall libcb.cb_matrixcbsolver_set_prox(self.data::Ptr{Cvoid}, (isnothing(proxp) ? C_NULL : proxp.data)::Ptr{Cvoid})::Cint

@doc raw"""
    cb_set_active_bounds_fixing!(self::CBMatrixCBSolver, allow_fixing::Bool)

* @brief If set to true (the default is false), variables may be
        fixed automatically to active bounds if these are strongly
        active (i.e., the corresponding multipliers are big) and
  the center values are also right on these bounds already.

        The coordinates to be fixed are redetermined in each
        call following a descent step or a change of the function.
        An indicator vector of the variables fixed during the last call
        can be obtained via the routine get_fixed_active_bounds().

        Setting this value to true might improve the performance
        of the algorithm in some instances but there is no
        convergence theory. It might be particularly helpful
        within Lagrangian relaxation if a primal cutting plane
        approach is used and non-tight inequalities should be
        eliminated quickly (fixing then indicates large primal
        slack values as these are the dual variables to the bounds
  on the Lagrange mulitpliers). Furthermore, if the value
  of a variable is fixed to zero, the variable can typically
  be deleted without affecting the validity of the current
  cutting model and function values.

      @param[in] allow_fixing (bool)

      @return
        - 0 on success
        - != 0 otherwise
     
"""
cb_set_active_bounds_fixing!(self::CBMatrixCBSolver, allow_fixing::Bool) = @ccall libcb.cb_matrixcbsolver_set_active_bounds_fixing(self.data::Ptr{Cvoid}, allow_fixing::Cint)::Cvoid

@doc raw"""
    cb_clear_fail_counts!(self::CBMatrixCBSolver)

* @brief clears all fail counts on numerical function oder model failures,
  may be useful if this caused premature termination.

      @return
        - 0 on success
        - != 0 otherwise
    
"""
cb_clear_fail_counts!(self::CBMatrixCBSolver) = @ccall libcb.cb_matrixcbsolver_clear_fail_counts(self.data::Ptr{Cvoid})::Cvoid

@doc raw"""
    cb_set_eval_limit!(self::CBMatrixCBSolver, eval_limit::Integer)

* @brief Sets an upper bound on the number of calls to the oracle (use negative numbers for no limit).

        If this number is reached, the algorithm will terminate
  independently of whether the last step was a descent or
  a null step. A negative number will be interepreted as
  no limit.

      @param[in] eval_limit (Integer)

      @return
        - 0 on success
        - != 0 otherwise
    
"""
cb_set_eval_limit!(self::CBMatrixCBSolver, eval_limit::Integer) = @ccall libcb.cb_matrixcbsolver_set_eval_limit(self.data::Ptr{Cvoid}, eval_limit::Cint)::Cvoid

@doc raw"""
    cb_set_inner_update_limit!(self::CBMatrixCBSolver, update_limit::Integer)

* @brief Set an upper bound on the number of inner updates for the
        cutting model with primal slacks within one null step (use negative numbers for no limit).

        A negative number will be interepreted as no limit, i.e.,
        the updates will be done till a certain precision of the
        cutting model is achieved.

      @param[in] update_limit (Integer)

      @return
        - 0 on success
        - != 0 otherwise
    
"""
cb_set_inner_update_limit!(self::CBMatrixCBSolver, update_limit::Integer) = @ccall libcb.cb_matrixcbsolver_set_inner_update_limit(self.data::Ptr{Cvoid}, update_limit::Cint)::Cvoid

@doc raw"""
    cb_set_time_limit!(self::CBMatrixCBSolver, time_limit::Integer)

* @brief Set an upper bound on the number of seconds (user time, use negative numbers for no limit)

      @param[in] time_limit (Integer)

      @return
        - 0 on success
        - != 0 otherwise
    
"""
cb_set_time_limit!(self::CBMatrixCBSolver, time_limit::Integer) = @ccall libcb.cb_matrixcbsolver_set_time_limit(self.data::Ptr{Cvoid}, time_limit::Cint)::Cvoid

@doc raw"""
    cb_set_qp_solver!(self::CBMatrixCBSolver, qpparams::Union{<:CBQPSolverParametersObject,Nothing}, newqpsolver::Union{<:CBQPSolverObject,Nothing} = nothing)

 * @brief Set parameters for the internal QP solver, possibly after first exchanging the solver with a new one

      The objects passed need to be heap objects; their ownership is transferred
      to the bundle code and they will be deleted there. The pointers may be null,
      in which case nothing is done for this object

      @param[in] qpparams (QPSolverParametersObject*)

      @param[in] newqpsolver (QPSolverObject*)

      @return
        - 0 on success
        - != 0 otherwise
    
"""
cb_set_qp_solver!(self::CBMatrixCBSolver, qpparams::Union{<:CBQPSolverParametersObject,Nothing}, newqpsolver::Union{<:CBQPSolverObject,Nothing} = nothing) = @ccall libcb.cb_matrixcbsolver_set_qp_solver(self.data::Ptr{Cvoid}, (isnothing(qpparams) ? C_NULL : qpparams.data)::Ptr{Cvoid}, (isnothing(newqpsolver) ? C_NULL : newqpsolver.data)::Ptr{Cvoid})::Cint

@doc raw"""
    cb_get_dim(self::CBMatrixCBSolver)

* @brief Returns the current dimension of the design space/argument
           or -1 if no dimension is set.
    
"""
cb_get_dim(self::CBMatrixCBSolver) = @ccall libcb.cb_matrixcbsolver_get_dim(self.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_get_n_functions(self::CBMatrixCBSolver)

* @brief Returns the current number of functions in the problem.
    
"""
cb_get_n_functions(self::CBMatrixCBSolver) = @ccall libcb.cb_matrixcbsolver_get_n_functions(self.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_get_n_oracle_calls(self::CBMatrixCBSolver)

* @brief Returns the number of function evaluations
    
"""
cb_get_n_oracle_calls(self::CBMatrixCBSolver) = @ccall libcb.cb_matrixcbsolver_get_n_oracle_calls(self.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_get_n_descent_steps(self::CBMatrixCBSolver)

* @brief Returns the number of function descent setps
    
"""
cb_get_n_descent_steps(self::CBMatrixCBSolver) = @ccall libcb.cb_matrixcbsolver_get_n_descent_steps(self.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_get_n_inner_iterations(self::CBMatrixCBSolver)

* @brief Returns the number of inner iterations of the bundle method
    
"""
cb_get_n_inner_iterations(self::CBMatrixCBSolver) = @ccall libcb.cb_matrixcbsolver_get_n_inner_iterations(self.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_get_n_inner_updates(self::CBMatrixCBSolver)

* @brief Returns the number of inner multiplier updates for the box constraints
    
"""
cb_get_n_inner_updates(self::CBMatrixCBSolver) = @ccall libcb.cb_matrixcbsolver_get_n_inner_updates(self.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_get_descent_step(self::CBMatrixCBSolver)

* @brief returns true if the last evaluation of the last call to solve()
  resulted in a descent step

  Mind: if there was no (succesdful) evaluation, neither get_descent_step() nor
  get_null_step() will return true;
    
"""
cb_get_descent_step(self::CBMatrixCBSolver) = Bool(@ccall libcb.cb_matrixcbsolver_get_descent_step(self.data::Ptr{Cvoid})::Cint)

@doc raw"""
    cb_get_null_step(self::CBMatrixCBSolver)

* @brief returns true if the last evaluation of the last call to solve()
  resulted in a null step

  Mind: if there was no (successful) evaluation, neither get_descent_step() nor
  get_null_step() will return true;
    
"""
cb_get_null_step(self::CBMatrixCBSolver) = Bool(@ccall libcb.cb_matrixcbsolver_get_null_step(self.data::Ptr{Cvoid})::Cint)

@doc raw"""
    cb_get_costs(self::CBMatrixCBSolver, costs::CBMatrix)

* @brief If a linear cost vector was specified, costs will hold these values, otherwise the vector is initialized to zero (for the current dimension)
    
"""
cb_get_costs(self::CBMatrixCBSolver, costs::CBMatrix) = @ccall libcb.cb_matrixcbsolver_get_costs(self.data::Ptr{Cvoid}, costs.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_get_lbounds(self::CBMatrixCBSolver)

* @brief Returns a pointer to the vector of lower bounds or null if there is no such vector
    
"""
cb_get_lbounds(self::CBMatrixCBSolver) = CBMatrix(@ccall libcb.cb_matrixcbsolver_get_lbounds(self.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_get_ubounds(self::CBMatrixCBSolver)

* @brief Returns a pointer to the vector of upper bounds or null if there is no such vector
    
"""
cb_get_ubounds(self::CBMatrixCBSolver) = CBMatrix(@ccall libcb.cb_matrixcbsolver_get_ubounds(self.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_get_fixed_active_bounds(self::CBMatrixCBSolver)

* @brief Returns NULL or (iff active bound fixing is turned on
  in set_active_bounds_fixing()) the indicator vector of
  variables temporarily fixed to the center value due to
  significantly positive multipliers for the box constraints.

        A variable gets fixed to the bound only if the center is
  already a the bound and in some iteration the dual variables
  to the bound constraint indicate that the bound is strongly
  active also for the candidate. Of course this migh just hold
  for one candidate and there is no guarantee that the bound
  is also strongly active in an optimal solution. Thus, this
  mainly a heuristic to eliminate less important variables
  quickly from entering the subproblem.
     
"""
cb_get_fixed_active_bounds(self::CBMatrixCBSolver) = CBIndexmatrix(@ccall libcb.cb_matrixcbsolver_get_fixed_active_bounds(self.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_get_prox(self::CBMatrixCBSolver)

* @brief Returns the pointer to the current prox term of the bundle solver
    
"""
cb_get_prox(self::CBMatrixCBSolver) = CBBundleProxObject(@ccall libcb.cb_matrixcbsolver_get_prox(self.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_set_out!(self::CBMatrixCBSolver, print_level::Integer = 1)

* @brief Specifies the output level (out==NULL: no output at all,
           out!=NULL and level=0: errors and warnings,
           level>0 increasingly detailed information)

     @param[in] out  (std::ostream*)
       direct all output to (*out). If out==NULL, there will be no output at all.

     @param[in] print_level (int)

     Output levels for print_level:
      -  0 ... no output except for errors and warnings
      -  1 ... line summary after each descent step
      - >1 ... undocumented and increasingly detailed log information.
             These higher levels should only be used if requested
             for debugging purposes.

      Example for level 1:

\verbatim
00:00:00.00 endit  1   1   1   563.    563.  39041.188  39043.162
00:00:00.00 endit  2   2   2   563.    559.  38488.165  38490.200
00:00:00.00 endit  3   3   3   56.3    555.  33014.533  33211.856
00:00:00.00 endit  4   4   4   5.63    517. -14306.459  2738.0343
00:00:00.00 endit  5   5   5   4.04    148. -2692.1131  2.2150883
00:00:00.00 endit  6   6   6   4.01    1.29  1.7908952  2.0000581
00:00:00.00 endit  7   7   7   3.95  0.0213  1.9999387  2.0000000
00:00:00.00 _endit  8   8   8   3.95 2.94e-05  2.0000000  2.0000000

Column 1      2     3   4   5    6       7       8          9
\endverbatim
      - Column 1: computation time in hh:mm:ss.dd,
      - Column 2: "endit" is convenient for grep and stands for "end of iteration".
         Iterations with termination_code()!=0 are marked with "_endit".
      - Column 3: number of descent steps
      - Column 4: number of descent and null steps. Up to initialization calls
         and reevaluations, this is the number of evaluation calls
         to the function oracles from within the bundle method.
         In the example all calls led to descent steps.
      - Column 5: number of innermost iterations. It differs from column 5 only in the
          case of variables with bounds in which case it gives the number of updates
          of the multipliers for the bounds (or primal slacks in Lagrangean
          relaxation). Exceedingly high numbers in this column indicate that
          some variables are constantly at their bounds and it might be
          possible to improve convergence by deleting them (i.e. set them
          as constants to this bound and remove the variable).
      - Column 6: the weight of the quadratic term in the augmented problem.
      - Column 7: the norm of the aggregate subgradient. If it is small,
          say below 0.1, then mostly this is good indication that the
          objective value is close to optimal.
      - Column 8: the value of the cutting model in the last candidate point. It
          is always a lower bound on the true function value in this point
      - Column 9: the objective value in the latest point that led to a descent
          step, i.e., the point returend by get_center(). Whenever
          termination_code() returns 0 this is also the objective
          value of the latest evaluation call to the function oracles and
          the value in the center point of the next iteration.
     
"""
cb_set_out!(self::CBMatrixCBSolver, print_level::Integer = 1) = @ccall libcb.cb_matrixcbsolver_set_out(self.data::Ptr{Cvoid}, print_level::Cint)::Cvoid

@doc raw"""
    cb_print_line_summary(self::CBMatrixCBSolver)

print a one line summary of important evaluation data
"""
cb_print_line_summary(self::CBMatrixCBSolver) = @ccall libcb.cb_matrixcbsolver_print_line_summary(self.data::Ptr{Cvoid})::Cvoid

@doc raw"""
    cb_print_statistics(self::CBMatrixCBSolver)

print a cryptic summary of computation times of important components
"""
cb_print_statistics(self::CBMatrixCBSolver) = @ccall libcb.cb_matrixcbsolver_print_statistics(self.data::Ptr{Cvoid})::Cvoid

cb_get_solver(self::CBMatrixCBSolver) = CBBundleSolver(@ccall libcb.cb_matrixcbsolver_get_solver(self.data::Ptr{Cvoid})::Ptr{Cvoid})

