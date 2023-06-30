@doc raw"""
    cb_clear!(self::CBUnconstrainedGroundset, indim::Integer = 0, in_groundset_id::Integer = 0)

* @brief reset everything to initial state for an unconstrained ground set of dimension @a indim

        Note that a ground set is allowed to have dimension zero. This
        will lead to evaluating a function without arguments and is
        a realistic scenario in Lagrangean relaxation of cutting plane
        approaches if no cutting planes have been added yet.

        If the changes to the ground set are to be counted by groundset_id,
        then it makes sense to enter the appropriate value in in_groundset_id.
     
"""
cb_clear!(self::CBUnconstrainedGroundset, indim::Integer = 0, in_groundset_id::Integer = 0) = @ccall libcb.cb_unconstrainedgroundset_clear(self.data::Ptr{Cvoid}, indim::Cint, in_groundset_id::Cint)::Cvoid

@doc raw"""
    CBUnconstrainedGroundset(indim::Integer = 0, start_val::Union{<:CBMatrix,Nothing} = nothing, costs::Union{<:CBMatrix,Nothing} = nothing, offset::Real = 0., in_groundset_id::Integer = 0)

calls clear() with the same parameters
"""
CBUnconstrainedGroundset(indim::Integer = 0, start_val::Union{<:CBMatrix,Nothing} = nothing, costs::Union{<:CBMatrix,Nothing} = nothing, offset::Real = 0., in_groundset_id::Integer = 0) = CBUnconstrainedGroundset(@ccall libcb.cb_unconstrainedgroundset_new(indim::Cint, (isnothing(start_val) ? C_NULL : start_val.data)::Ptr{Cvoid}, (isnothing(costs) ? C_NULL : costs.data)::Ptr{Cvoid}, offset::Cdouble, in_groundset_id::Cint)::Ptr{Cvoid})

@doc raw"""
    cb_get_groundset_id(self::CBUnconstrainedGroundset)

returns the current groundset_id, increased values indicate changes in the ground set
"""
cb_get_groundset_id(self::CBUnconstrainedGroundset) = @ccall libcb.cb_unconstrainedgroundset_get_groundset_id(self.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_set_groundset_id!(self::CBUnconstrainedGroundset, gsid::Integer)

sets the groundset_id to the desired value, increasing it is safer here because this is used to indicate changes
"""
cb_set_groundset_id!(self::CBUnconstrainedGroundset, gsid::Integer) = @ccall libcb.cb_unconstrainedgroundset_set_groundset_id(self.data::Ptr{Cvoid}, gsid::Cint)::Cvoid

@doc raw"""
    cb_get_dim(self::CBUnconstrainedGroundset)

returns the dimension of the ground set, i.e., the length of the variables vector y
"""
cb_get_dim(self::CBUnconstrainedGroundset) = @ccall libcb.cb_unconstrainedgroundset_get_dim(self.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_constrained(self::CBUnconstrainedGroundset)

* @brief returns false if the feasible set is the entire space (unconstrained optimization), true otherwise.

        The current class implements the unconstrained case and always returns false.
    
"""
cb_constrained(self::CBUnconstrainedGroundset) = Bool(@ccall libcb.cb_unconstrainedgroundset_constrained(self.data::Ptr{Cvoid})::Cint)

@doc raw"""
    cb_is_feasible!(self::CBUnconstrainedGroundset, y::CBMatrix, relprec::Real = 1e-10)

* @brief on input value in_groundset_id the input y was feasible. Return true if the id did not change, otherwise check if y is still feasible for the given precision.

     The routine is called by the internal bundle solver to check
     whether the given center is still valid (in some applications the
     groundset might change during the runtime of the bundle method),
     where validity of y was already checked at a point in time when the
     groundset had the in_groundset_id. If the groundset_id is still the
     same, then y is simply assumed to be still correct (the precision
     is not even looked at in this case). Otherwise the routine checks
     the validitiy of y with respect to the given precision but does not
     enforce validity. It returns true if y is valid and false otherwise.
    
"""
function cb_is_feasible!(self::CBUnconstrainedGroundset, y::CBMatrix, relprec::Real = 1e-10)
    in_groundset_id = Ref{Int}()
    Bool(@ccall libcb.cb_unconstrainedgroundset_is_feasible(self.data::Ptr{Cvoid}, in_groundset_id::Ref{Int}, y.data::Ptr{Cvoid}, relprec::Cdouble)::Cint)
    return in_groundset_id[]
end

@doc raw"""
    cb_ensure_feasibility!(self::CBUnconstrainedGroundset, y::CBMatrix, ychanged::Bool, Hp::Union{<:CBBundleProxObject,Nothing} = nothing, relprec::Real = 1e-10)

* @brief if the groundset_id changed, it checks feasibility of y with respect to the given precision. If infeasible it replaces y by its projection with respect to the norm of Hp and sets ychanged to true.

     The routine is called by the internal bundle solver to check
     whether the given center is still valid (in some applications the
     groundset might change during the runtime of the bundle method),
     where, if @a ychanged is false on input, validity of @a y was
     already checked at a point in time when the groundset had
     the @a in_groundset_id. If @a ychanged==false and the groundset_id is
     still the same, then @a y is simply assumed to be still correct
     (the precision is not even looked at in this case).  Otherwise
     the routine checks the validitiy of @a y with respect to the given
     precision. If feasible, it returns the new
     groundset_id in @a in_groundset_id and keeps @a ychanged unaltered.
     If @a y is infeasible, the rountine computes its projection
     onto the feasible set with respect to the norm of @a Hp (if ==0 then
     the Euclidean norm is used), stores it in @a y, sets @a ychanged to
     true, sets @a in_groundset_id to the current groundset_id and returns
     0. Should anything go wrong, it returns 1.

     This concrete base class represents the unconstrained case, so
     feasiblity only checks the dimension and never requires projections.
     
"""
function cb_ensure_feasibility!(self::CBUnconstrainedGroundset, y::CBMatrix, ychanged::Bool, Hp::Union{<:CBBundleProxObject,Nothing} = nothing, relprec::Real = 1e-10)
    in_groundset_id = Ref{Int}()
    @ccall libcb.cb_unconstrainedgroundset_ensure_feasibility(self.data::Ptr{Cvoid}, in_groundset_id::Ref{Int}, y.data::Ptr{Cvoid}, ychanged::Ref{Cint}, (isnothing(Hp) ? C_NULL : Hp.data)::Ptr{Cvoid}, relprec::Cdouble)::Cint
    return in_groundset_id[]
end

@doc raw"""
    cb_get_qp_solver!(self::CBUnconstrainedGroundset, solves_model_without_gs::Bool, Hp::Union{<:CBBundleProxObject,Nothing})

returns a pointer to an internal QPSolverObject that is able to solve bundle suproblems efficiently for this kind of groundset and scaling; if solves_model_without_gs == true the qp solver does not include the groundset and the groundset has to be dealt with by the Gauss Seidel approach
"""
cb_get_qp_solver!(self::CBUnconstrainedGroundset, solves_model_without_gs::Bool, Hp::Union{<:CBBundleProxObject,Nothing}) = CBQPSolverObject(@ccall libcb.cb_unconstrainedgroundset_get_qp_solver(self.data::Ptr{Cvoid}, solves_model_without_gs::Ref{Cint}, (isnothing(Hp) ? C_NULL : Hp.data)::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_set_qp_solver_parameters!(self::CBUnconstrainedGroundset, param0::Union{<:CBQPSolverParametersObject,Nothing})

set parameters for the QP_Solver
"""
cb_set_qp_solver_parameters!(self::CBUnconstrainedGroundset, param0::Union{<:CBQPSolverParametersObject,Nothing}) = @ccall libcb.cb_unconstrainedgroundset_set_qp_solver_parameters(self.data::Ptr{Cvoid}, (isnothing(param0) ? C_NULL : param0.data)::Ptr{Cvoid})::Cint

@doc raw"""
    cb_get_starting_point(self::CBUnconstrainedGroundset)

returns a stored starting point, note: this need not be feasible; if generated automatically, its dimension is correct.
"""
cb_get_starting_point(self::CBUnconstrainedGroundset) = (@ccall libcb.cb_unconstrainedgroundset_get_starting_point(self.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

@doc raw"""
    cb_set_starting_point!(self::CBUnconstrainedGroundset, vec::CBMatrix)

stores the a new starting point irrespective of whether it is feasible or not and returns 0 if it feasible, 1 if it is infeasible
"""
cb_set_starting_point!(self::CBUnconstrainedGroundset, vec::CBMatrix) = @ccall libcb.cb_unconstrainedgroundset_set_starting_point(self.data::Ptr{Cvoid}, vec.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_candidate!(self::CBUnconstrainedGroundset, newy::CBMatrix, center_y::CBMatrix, center_value::Real, model_minorant::CBMinorantPointer, Hp::Union{<:CBBundleProxObject,Nothing}, delta_groundset_minorant::Union{<:CBMinorantPointer,Nothing} = nothing, delta_index::Union{<:CBIndexmatrix,Nothing} = nothing, relprec::Real = 1e-2)

* @brief for a given model aggregate compute the groundset aggregate and the resulting (feasible) candidate

  Let \f$ (\sigma,s) \f$ be the aggregate minorant \f$\sigma+s^\top y\le
  f(y)\f$ of the cost function described by @a model_subg_offset and @a
  model_subg, let \f$\hat y\f$ be the center of stability given by @a
  center_y, let \f$H\f$ denote the positive definite scaling matrix with
  weight \f$u\f$ given by @a Hp, let \f$Y\f$ denote this (convex)
  feasible ground set and \f$i_Y\f$ its indicator function, then this
  computes a saddle point \f$(y,(\gamma,g))\f$ of

  \f[ \max_{(\gamma,g)\in\partial i_Y}\min_{y\in\mathbf{R}^n} [(s+g)^\top y + \gamma+\sigma +\frac{u}2\|y-\hat y\|_H^2],\f]

  The resulting \f$y\f$ is feasible and stored in @a newy.
       @a linval gets \f$\gamma+\sigma+(g+s)^\top y\f$,
       @a subgnorm2 gets \f$u\|y-\hat x\|_H^2\f$,
       @a augval gets @a linval+@a subgnorm2 /2.

  If @a delta_groundset_subg is not NULL, also
       @a elta_groundset_subg_offset and @a delta index are assumed to be
       not NULL. Then all changes from the previous groundset aggregate
       minorant to the new groundset aggregate minorant \f$(\gamma,g)\f$ are
       stored in @a delta_groundset_subg_offset and @a delta_groundset_subg
       and @a delta_index holds the indices of the nonzero changes
       (mostly the groundset aggregate is sparse, e.g. due to complementarity).
       This is used in BundleProxObject::update_QP_costs().

  On input \f$\varepsilon=\f$ @a relprec, \f$\bar f=\f$ @a
       center_value, and \f$\underline{f}=\f$ @a augval serve to form
       appropriate stopping criteria if solving the saddle point problem
       requires a nonlinear convex optimization method. The method is
       then assumed to produce a primal solution \f$y\in Y\f$ of value
       \f$\bar s\f$ and a dual solution \f$(\gamma,g)\f$ of value
       \f$\underline{g}\f$ with the properties \f$\underline g-\underline
       f\ge 0\f$ and \f$\bar f-\bar s\ge 0\f$ and \f$\bar
       s-\underline{g}\le\varepsilon(|\bar s|+1.)\f$.  In particular
       \f$y\in Y\f$ is assumed to hold to machine precision.

  If use_yfixing is true (the fixing heuristic is switched on),
       then \f$y_i=\hat y_i\f$ is required to hold for all i with
       yfixed(i)!=0, so these coordinates are not allowed to change. In
       particular, this routine may also set yfixed(i)=2 for new
       coordinates i, where 2 is used to indicate newly fixed
       variables. These will be reset to 1 in
       BundlesScaling::compute_QP_costs() and
       BundlesScaling::update_QP_costs()
       when this information has been digested.

  In some derived classes a scaling heuristic is called that may influence
  the scaling @a Hp so as to avoid going outside the feasible region too far.

    @param[out] gs_id
        the current groundset_id

    @param[out] newy
        the next candidate y (feasible)

    @param[out] cand_gs_val
        the value of the groundset minorant in the candidate y (=groundset objective)

    @param[out] linval (CH_Matrix_Classes::Real&) value of linear minorant in y

    @param[in,out] augval_lb (CH_Matrix_Classes::Real&)
        - on input: lower bound value of augmented model in previous (maybe infeasible) candidate
        - on output: lower bound value of augmented model in y

    @param[out] augval_ub (CH_Matrix_Classes::Real&)
        - on output: upper bound value of augmented model in y

    @param[out] subgnorm2 (CH_Matrix_Classes::Real&)
        squared Hp-norm (with weight) of joint groundset and model aggregate

    @param[in] center_y (const CH_Matrix_Classes::Matrix&)
        center of stability (feasible)

    @param[in] center_value (CH_Matrix_Classes::Real)
        function value in center_y

    @param[in] model_minorant (const MinorantPoiner&)
        aggregate linear minorant of the cost function

    @param[in,out] Hp (ConicBundle::BundleProxObject*)
        pointer to scaling matrix H, may be influenced by a scaling heuristic

    @param[in,out] delta_groundset_minorant (MinorantPointer*)
        if not NULL, the change in groundset aggregate will be stored here

    @param[in,out] delta_index (CH_Matrix_Classes::Indexmatrix*)
        must be not NULL iff delta_groundset_subg!=NULL or yfixed has changed,
        will store nonzero indices of delta_groundset_subg

    @param[in] relprec (CH_Matrix_Classes::Real)
        relative precision for termination in QP computations

    @return 0 on success, != 0 on failure

    
"""
function cb_candidate!(self::CBUnconstrainedGroundset, newy::CBMatrix, center_y::CBMatrix, center_value::Real, model_minorant::CBMinorantPointer, Hp::Union{<:CBBundleProxObject,Nothing}, delta_groundset_minorant::Union{<:CBMinorantPointer,Nothing} = nothing, delta_index::Union{<:CBIndexmatrix,Nothing} = nothing, relprec::Real = 1e-2)
    subgnorm2 = Ref{Float64}()
    augval_ub = Ref{Float64}()
    augval_lb = Ref{Float64}()
    linval = Ref{Float64}()
    cand_gs_val = Ref{Float64}()
    gs_id = Ref{Int}()
    @ccall libcb.cb_unconstrainedgroundset_candidate(self.data::Ptr{Cvoid}, gs_id::Ref{Int}, newy.data::Ptr{Cvoid}, cand_gs_val::Ref{Float64}, linval::Ref{Float64}, augval_lb::Ref{Float64}, augval_ub::Ref{Float64}, subgnorm2::Ref{Float64}, center_y.data::Ptr{Cvoid}, center_value::Cdouble, model_minorant.data::Ptr{Cvoid}, (isnothing(Hp) ? C_NULL : Hp.data)::Ptr{Cvoid}, (isnothing(delta_groundset_minorant) ? C_NULL : delta_groundset_minorant.data)::Ptr{Cvoid}, (isnothing(delta_index) ? C_NULL : delta_index.data)::Ptr{Cvoid}, relprec::Cdouble)::Cint
    return gs_id[], cand_gs_val[], linval[], augval_lb[], augval_ub[], subgnorm2[]
end

@doc raw"""
    cb_get_gs_aggregate(self::CBUnconstrainedGroundset)

returns the groundset aggregate computed in candidate()
"""
cb_get_gs_aggregate(self::CBUnconstrainedGroundset) = (@ccall libcb.cb_unconstrainedgroundset_get_gs_aggregate(self.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

@doc raw"""
    cb_get_gs_minorant(self::CBUnconstrainedGroundset)

returns the linear minorant valid on the entire ground set (e.g. a linear cost funciton)
"""
cb_get_gs_minorant(self::CBUnconstrainedGroundset) = (@ccall libcb.cb_unconstrainedgroundset_get_gs_minorant(self.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

@doc raw"""
    cb_get_yfixed(self::CBUnconstrainedGroundset)

if not NULL (iff get_use_yfixing()==false) it returns the vector yfixed with yfixed(i)=0 if not fixed, =1 is fixed already, =2 if newly fixed
"""
cb_get_yfixed(self::CBUnconstrainedGroundset) = CBIndexmatrix(@ccall libcb.cb_unconstrainedgroundset_get_yfixed(self.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_set_yfixed!(self::CBUnconstrainedGroundset)

if not NULL (iff get_use_yfixing()==false) returns the vector yfixed with yfixed(i)=0 if not fixed, =1 is fixed already, =2 if newly fixed
"""
cb_set_yfixed!(self::CBUnconstrainedGroundset) = CBIndexmatrix(@ccall libcb.cb_unconstrainedgroundset_set_yfixed(self.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_get_use_yfixing(self::CBUnconstrainedGroundset)

true if the cooridinate fixing heuristic is switched on (only constrained cases)
"""
cb_get_use_yfixing(self::CBUnconstrainedGroundset) = Bool(@ccall libcb.cb_unconstrainedgroundset_get_use_yfixing(self.data::Ptr{Cvoid})::Cint)

@doc raw"""
    cb_set_use_yfixing!(self::CBUnconstrainedGroundset, uyf::Bool)

set to true to switch on the cooridinate fixing heuristic (only constrained cases)
"""
cb_set_use_yfixing!(self::CBUnconstrainedGroundset, uyf::Bool) = @ccall libcb.cb_unconstrainedgroundset_set_use_yfixing(self.data::Ptr{Cvoid}, uyf::Cint)::Cvoid

@doc raw"""
    cb_start_modification!(self::CBUnconstrainedGroundset)

return a new modification object on the heap that is initialized for modification of *this
"""
cb_start_modification!(self::CBUnconstrainedGroundset) = CBGroundsetModification(@ccall libcb.cb_unconstrainedgroundset_start_modification(self.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_apply_modification!(self::CBUnconstrainedGroundset, mdf::CBGroundsetModification)

change the groundset description as specified by the argument
"""
cb_apply_modification!(self::CBUnconstrainedGroundset, mdf::CBGroundsetModification) = @ccall libcb.cb_unconstrainedgroundset_apply_modification(self.data::Ptr{Cvoid}, mdf.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_mfile_data(self::CBUnconstrainedGroundset)

m-file output routine for debugging or testing in Matlab (not yet working)
"""
cb_mfile_data(self::CBUnconstrainedGroundset) = @ccall libcb.cb_unconstrainedgroundset_mfile_data(self.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_set_cbout!(self::CBUnconstrainedGroundset, incr::Integer = -1)

output settings
"""
cb_set_cbout!(self::CBUnconstrainedGroundset, incr::Integer = -1) = @ccall libcb.cb_unconstrainedgroundset_set_cbout(self.data::Ptr{Cvoid}, incr::Cint)::Cvoid

