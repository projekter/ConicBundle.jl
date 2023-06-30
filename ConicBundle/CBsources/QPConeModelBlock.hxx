/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/QPConeModelBlock.hxx
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


#ifndef CONICBUNDLE_QPCONEMODELBLOCK_HXX
#define CONICBUNDLE_QPCONEMODELBLOCK_HXX

/**  @file QPConeModelBlock.hxx
    @brief Header declaring the class ConicBundle::QPConeModelBlock
    @version 1.0
    @date 2020-03-23
    @author Christoph Helmberg
*/


#include "QPModelBlock.hxx"
#include "NNCIPBundleBlock.hxx"
#include "SOCIPBundleBlock.hxx"
#include "PSCIPBundleBlock.hxx"
#include "BoxIPBundleBlock.hxx"

namespace ConicBundle {


  /** @ingroup ConstrainedQPSolver
   */

   //@{

   /** @brief  implements a QPModelBlock for conic cutting models in QPSolver

      note, the bundle at the final level might not be the same as the
      one entered in the intial call, because it might have undergone
      affine function transformations.

      The starting point of the model is the call to init() as described
      in QPConeModelDataObject::init(), which sets
      the bundle (subgradient or minorant information)
      with the corresponding description of the model variables and
      constraints specifying on how the bundle inforamtion is to be
      combined. The description possibilities here allow to formulate
      the default ConicBundle models NNCModel, SOCModel, PSCModel and
      BoxModel. The sequence of the cones with their dimensions must
      match the sequence given in the bundle. All cones are combined
      by at most one single row trace constraint, which might even be
      missing in the case of a full box model. The right hand side
      of the trace constraint serves as function factor in the equality
      case or, in the inequalitiy case, as penalty parameter. This
      trace constraint (and if present, its slack variable)
      is handeld directly by this QPModleBlock; the interior point
      contributions of the cone parts are each handled by a corresponding
      InteriorPointBundleBlock implementation NNCIPBundleBlock,
      SOCIPBundleBlock, PSCIPBundleBlock or BoxIPBundleBlock.

      Each call to a QPModelBlockObject is then implemented as
      handling the coupling constraint combined with calling each
      InterioPointBundleBlock for its contribution.
   */

  class QPConeModelBlock : public QPModelBlock, public QPConeModelDataObject {
  private:
    //data that stays constant after init
    CH_Matrix_Classes::Integer modeldim; ///< combined size of the model variables
    CH_Matrix_Classes::Integer mu_dim; ///< combined size of the interioer point cone dimension relevant for computing the barrier parameter
    FunctionTask ft; ///< role of the coupling constraint (equality constraint for ObjectiveFunction, otherwise inequality constraint)
    bool use_trace;   ///< true iff the trace constraint is used at all

    //values related to the barrier parameter
    CH_Matrix_Classes::Real last_rhs_mu;  ///< the last mu used in rhs computations
    CH_Matrix_Classes::Real mu;           ///< in a step mu gets the value of last_rhs_mu
    CH_Matrix_Classes::Real old_mu;       ///< in a step old_mu gets the value of mu before this gets last_rhs_mu

    // primal-dual variables for the trace constraint
    CH_Matrix_Classes::Matrix trace_vec;  ///< coefficnet vector of the trace constraint, stays constant unless box takes part in scaling;
    CH_Matrix_Classes::Real trace_rhs; ///< right hand side value of the trace constraint
    CH_Matrix_Classes::Real trace_slack;  ///< value of the slack variable to the trace constraint if this is an inequality 
    CH_Matrix_Classes::Real trace_dual;  ///< dual variable to the trace constraint

    CH_Matrix_Classes::Real trace_delta_slack; ///< step of the trace slack
    CH_Matrix_Classes::Real trace_delta_dual; ///< step of the trace dual variable

    CH_Matrix_Classes::Real sys_trace;   ///< trace_slack/trace_dual + possibly a box contribution, initialized if >0.
    CH_Matrix_Classes::Real complrhs_trace;  ///< stores the complementarity value determined in computing the right hand side value

    // data that stays fixed as long as the working point is not changed
    CH_Matrix_Classes::Matrix diff_model;  ///< trace_dual*trace_vec - bundlevalues (should be == z)
    CH_Matrix_Classes::Real diff_trace;    ///< sum of traces - trace_rhs  (should be == -trace_slack)

    // data that records the rhs used depending on predictor/corrector settings 
    CH_Matrix_Classes::Matrix sysrhs_model; ///< KKT system right hand side computed for the model 
    CH_Matrix_Classes::Real sysrhs_trace; ///< KKT system right hand side computed for the trace constraint
    CH_Matrix_Classes::Matrix sysinv_trace;   ///< gets apply_xizinv(trace_vec), only formed when needed 
    CH_Matrix_Classes::Matrix Btsysinv_trace;   ///< gets transpose(B)*sysinv_trce, only formed when needed 
    CH_Matrix_Classes::Real schur_trace;   ///< gets ip(trace_vec,sysinv_trace)[+slack/dual]>0!, only formed when needed, otherwise -1.

    // basic interior point blocks for the cones
    NNCIPBundleBlock* nncblock; ///< if not NULL it points to the NNC cone of the model 
    std::vector<SOCIPBundleBlock*> socblock; ///< the SOC blocks of the model
    std::vector<PSCIPBundleBlock*> pscblock; ///< the PSC blocks of the model
    BoxIPBundleBlock* boxblock; ///< the box block of the model

    std::vector<InteriorPointBundleBlock*> block; ///< all active blocks lined up for easier loop design

    QPModelOracleDataObject* oracle_data; ///< the pointer to the oracle_data, if given on init


    /// determine the value of the trace for the given model variable
    CH_Matrix_Classes::Real evaluate_trace(CH_Matrix_Classes::Matrix& vec);

    ///< if the model variable changes, this resets all values computed for the oblivious model values
    void modelx_changed();

  public:

    /// C=beta*C+alpha*B*A where B and A may be transposed; carry out the model part of this beginning at startindex_model and beta for the part, that is added to (the calling routine has to make sure beta is not executed repeatedly if the same part is affected by other models as well)
    virtual CH_Matrix_Classes::Matrix&
      B_times(const CH_Matrix_Classes::Matrix& A,
        CH_Matrix_Classes::Matrix& C,
        CH_Matrix_Classes::Real alpha,
        CH_Matrix_Classes::Real beta,
        int Btrans,
        int Atrans,
        CH_Matrix_Classes::Integer startindex_model,
        MinorantBundle& globalbundle,
        CH_Matrix_Classes::Integer startindex_bundle);

    /// C=beta*C+alpha*A*B where A and B may be transposed; carry out the model part of this beginning at startindex_model 
    virtual CH_Matrix_Classes::Matrix&
      times_B(const CH_Matrix_Classes::Matrix& A,
        CH_Matrix_Classes::Matrix& C,
        CH_Matrix_Classes::Real alpha,
        CH_Matrix_Classes::Real beta,
        int Atrans,
        int Btrans,
        CH_Matrix_Classes::Integer startindex_model,
        MinorantBundle& globalbundle,
        CH_Matrix_Classes::Integer startindex_bundle);

    ///add the main diagonal block tranpose(projection)*diagvec*projection to bigS starting at startindex
    virtual CH_Matrix_Classes::Symmatrix&
      add_BDBt(const CH_Matrix_Classes::Matrix& diagvec,
        CH_Matrix_Classes::Symmatrix& bigS,
        bool minus,
        CH_Matrix_Classes::Integer startindex,
        CH_Matrix_Classes::Matrix& Bt,
        CH_Matrix_Classes::Integer startindex_model,
        MinorantBundle& globalbundle,
        CH_Matrix_Classes::Integer startindex_bundle);

    /// get the current matrix for the coupling matrix Bt in the first row of blocks
    virtual CH_Matrix_Classes::Matrix&
      get_Bt(CH_Matrix_Classes::Matrix& Bt,
        CH_Matrix_Classes::Integer startindex_model,
        MinorantBundle& global_bundle,
        CH_Matrix_Classes::Integer startindex_bundle);

    /// set the local modelx value in modelx beginning with startindex (initialize it, do not add)
    virtual int get_modelx(CH_Matrix_Classes::Matrix& modelx,
      CH_Matrix_Classes::Integer startindex_model);

    /// set the local modeldx value in modeldx beginning with startindex (initialize it, do not add)
    virtual int get_modeldx(CH_Matrix_Classes::Matrix& modeldx,
      CH_Matrix_Classes::Integer startindex_model);

    /// set the local modeldcstr value in modeldcstr beginning with startindex (initialize it, do not add)
    virtual int get_modeldcstr(CH_Matrix_Classes::Matrix& modeldcstr,
      CH_Matrix_Classes::Integer startindex_constraints);

    /// adds opB transposed times modelx (with offsets but without constant affine term) to the arguments
    virtual int add_modelx_aggregate(CH_Matrix_Classes::Real& val,
      CH_Matrix_Classes::Matrix& vec,
      MinorantBundle& global_bundle,
      CH_Matrix_Classes::Integer startindex_bundle);

    /// set the model violation for the current system solution 
    virtual int get_sysviol_model(CH_Matrix_Classes::Matrix& modelvec,
      CH_Matrix_Classes::Integer startindex_model,
      const CH_Matrix_Classes::Matrix& y_plus_dy,
      MinorantBundle& global_bundle,
      CH_Matrix_Classes::Integer startindex_bundle);

    /// set the constraint violation for the current system solution starting at this index
    virtual int get_sysviol_constraints(CH_Matrix_Classes::Matrix& constrvec,
      CH_Matrix_Classes::Integer startindex_constr);

    virtual void display_model_values(const CH_Matrix_Classes::Matrix& y,
      MinorantBundle& global_bundle,
      CH_Matrix_Classes::Integer startindex_bundle,
      std::ostream& out);

    // bundlevalues holds the negative evaluation of the bundle for the current y 
    virtual int reset_starting_point(const CH_Matrix_Classes::Matrix& y,
      CH_Matrix_Classes::Real mu,
      MinorantBundle& global_bundle,
      CH_Matrix_Classes::Integer startindex_bundle);

    /// compute the step in the model space given the step in the design space
    virtual int compute_step(const CH_Matrix_Classes::Matrix& ystep,
      MinorantBundle& global_bundle,
      CH_Matrix_Classes::Integer startindex_bundle);

    /// store this computed step locally and compute the missing local dual step information
    virtual int computed_step(const CH_Matrix_Classes::Matrix& modelxstep,
      CH_Matrix_Classes::Integer startindex_model,
      const CH_Matrix_Classes::Matrix& modelconstrstep,
      CH_Matrix_Classes::Integer startindex_constr);

    /// move in the last computed step direction by a step of length alpha and compute and store the violation in this point for later use in 
    virtual int do_step(CH_Matrix_Classes::Real alpha,
      const CH_Matrix_Classes::Matrix& y,
      MinorantBundle& global_bundle,
      CH_Matrix_Classes::Integer startindex_bundle);

    /// If mu is not zero, always add the centering term for this mu as well, if append is false, add the Schur complement rhs for add_BtinvsysB, if append is true, fill in the rhs of the local system starting at startindex for the model and at startindex_constraints for the modelconstraints
    virtual int add_localrhs(CH_Matrix_Classes::Matrix& globalrhs,
      CH_Matrix_Classes::Real rhsmu,
      CH_Matrix_Classes::Real rhscorr,
      CH_Matrix_Classes::Integer startindex_model,
      CH_Matrix_Classes::Integer startindex_constraints,
      bool append,
      MinorantBundle& bundle,
      CH_Matrix_Classes::Integer startindex_bundel);

    ///add the "scaled" minorant outer products to globalsys, where the correct minorants start at the given index
    virtual int add_BtinvsysB(CH_Matrix_Classes::Symmatrix& globalsys,
      MinorantBundle& bundle,
      CH_Matrix_Classes::Integer startindex_bundle);

    /** @brief  given the Cholesky factorization LL' of *minus* the blocks A and B (contraints on design variables and Bundle-modelx) and LinvABrhs, solve for the local constraints C and add the new contribution of tracedual*LinvTrace to LinvABsol; store the tracedual in Crhs_and_sol but not yet locally (this will be done by computed_step() ).
     */
    virtual int solve_constrsys(const CH_Matrix_Classes::Symmatrix& ABchol,
      const CH_Matrix_Classes::Matrix& LinvABrhs,
      CH_Matrix_Classes::Matrix& LinvABsol,
      CH_Matrix_Classes::Integer startindex_model,
      CH_Matrix_Classes::Matrix& Crhs_and_sol,
      CH_Matrix_Classes::Integer startindex_constraints);



  public:
    /// reset to "empty/no" model
    void clear();

    /// default constructor
    QPConeModelBlock(CBout* cb = 0, int cbinc = -1);

    /// virtual destructor
    virtual ~QPConeModelBlock();

    /// return a cloned object on the heap
    virtual QPModelBlockObject* clone();

    /// usually the objects of the recursive block structure and not deleted in a clear. If needed, this can be invoked explicitly here, e.g., in order to clean up clones
    virtual void recursive_delete_and_clear();

    /// sofar this is only needed for some comparative evaluations; to work *this must be a clone of the the argument and the content of the objects is copied recursively
    virtual int recursive_copy_data_of(QPModelBlockObject*);

    //------------ QPModelBlockObject routines

    /// sets up the model with bundle information and how to combine it, see QPConeModelDataObject::init() for a detailed description
    int init(const MinorantPointer& constant_minorant,
      const MinorantBundle& bundle,
      CH_Matrix_Classes::Integer nnc_dim,
      const CH_Matrix_Classes::Indexmatrix& soc_dim,
      const CH_Matrix_Classes::Indexmatrix& psc_dim,
      const CH_Matrix_Classes::Matrix& box_lb,
      const CH_Matrix_Classes::Matrix& box_ub,
      CH_Matrix_Classes::Real b,
      FunctionTask ft,
      QPModelOracleDataObject* oracle_data = 0,
      bool scale_box = true);

    /// change the right hand side of the trace constraint to b
    int adjust_trace(CH_Matrix_Classes::Real b);


    /// evaluate the left hand side of the trace constraint for modelx
    CH_Matrix_Classes::Real evaluate_trace() const;

    /// get the right hand side of the trace constraint
    CH_Matrix_Classes::Real get_trace();

    /// get the linear part of modelx (and a guess, which of them are active, in {0.,1.})
    int get_nncx(CH_Matrix_Classes::Matrix& nncx,
      CH_Matrix_Classes::Matrix* nncx_activity = 0,
      bool cautious = false);

    /// get the SOC part of modelx (and a guess whether the entire cone is active
    int get_socx(CH_Matrix_Classes::Integer i,
      CH_Matrix_Classes::Matrix& socx,
      CH_Matrix_Classes::Real* socx_activity,
      bool cautious = false);

    /// get the PSC part of modelx (and a guess on the rank of the active part)
    int get_pscx(CH_Matrix_Classes::Integer i,
      CH_Matrix_Classes::Matrix& pscx_eigs,
      CH_Matrix_Classes::Matrix& pscx_vecs,
      CH_Matrix_Classes::Real& pscx_growthrate,
      CH_Matrix_Classes::Matrix& pscx_primalgrowth,
      CH_Matrix_Classes::Matrix& pscx_dualgrowth);

    /// get the box part of modelx (and a guess, which of the bounds are active, in {0.,1.})
    int get_boxx(CH_Matrix_Classes::Matrix& boxx,
      CH_Matrix_Classes::Matrix* linx_activity = 0,
      bool cautious = false);

    /// return the value of the dual variable to the trace consrat == support function value
    CH_Matrix_Classes::Real tracedual(CH_Matrix_Classes::Real* prec = 0) const;


    //-----------  QPBlock routines

    /// returns the dimension of the model set (here the same as the bundle size)
    virtual CH_Matrix_Classes::Integer dim_model();

    /// returns the dimension of the system describing the model set (may contain further constraints)
    virtual CH_Matrix_Classes::Integer dim_constraints();

    ///returns the dual upper bound to the model value (the trace weighted sum of the dual trace variables); it returns 0. if no model is contained
    virtual CH_Matrix_Classes::Real constraints_cost();

    ///return squared Euclidean norm of constraint violation of modelx
    virtual CH_Matrix_Classes::Real primalviol_2normsqr();

    ///return the squared Euclidean norm of the model violation  
    virtual CH_Matrix_Classes::Real dualviol_2normsqr();

    /// add dimensions of the primal-dual pairs to mudim and add the "trace" (the inner product with center) of the respective primal-dual pair products for the current step; update the min and max values of x_i*z_i
    virtual int get_mu_info(CH_Matrix_Classes::Integer& mudim,
      CH_Matrix_Classes::Real& tr_xz,
      CH_Matrix_Classes::Real& tr_xdzpdxz,
      CH_Matrix_Classes::Real& tr_dxdz,
      CH_Matrix_Classes::Real& min_xz,
      CH_Matrix_Classes::Real& max_xz) const;

    /// for limiting the stepsize with respect to the neighborhood this information about norms and inner products of x(.)*z-tr_xz-tr_xz/mudim(.*)1, x.()*dz+dx(.)*z-tr_xdzpdxz/mudim(.*)1, and dx(.)*dz-tr_dxdz/mudim(.)*1 is required, each block *adds* its contribution to the numbers
    virtual int get_nbh_info(CH_Matrix_Classes::Integer mudim,
      CH_Matrix_Classes::Real tr_xz,
      CH_Matrix_Classes::Real tr_xdzpdxz,
      CH_Matrix_Classes::Real tr_dxdz,
      CH_Matrix_Classes::Real nbh_ubnd,
      CH_Matrix_Classes::Real& alpha,
      CH_Matrix_Classes::Real& max_nbh,
      CH_Matrix_Classes::Real& nrmsqr_xz,
      CH_Matrix_Classes::Real& nrmsqr_xdzpdxz,
      CH_Matrix_Classes::Real& nrmsqr_dxdz,
      CH_Matrix_Classes::Real& ip_xz_xdzpdxz,
      CH_Matrix_Classes::Real& ip_xz_dxdz,
      CH_Matrix_Classes::Real& ip_dxdz_xdzpdxz) const;

    /// if necessary, reduce alpha to the biggest value so that feasibility is maintained with this step size
    virtual int linesearch(CH_Matrix_Classes::Real& alpha) const;

    ///add the local system for the model at startindex, the local for the constraints at startindex_constraints
    virtual int add_localsys(CH_Matrix_Classes::Symmatrix& globalsys,
      CH_Matrix_Classes::Integer startindex_model,
      CH_Matrix_Classes::Integer startindex_constraints);

    /** @brief multiply the local system diagonal block consisting of the model and local contraints rows and columns by in_vec[startindex_model+0,...,+dim_model(),startindex_constraints+0,...,+dim_constraints]  into the same coordinates of out_vec. */
    virtual int localsys_mult(const CH_Matrix_Classes::Matrix& in_vec,
      CH_Matrix_Classes::Matrix& out_vec,
      CH_Matrix_Classes::Integer startindex_model,
      CH_Matrix_Classes::Integer startindex_constraints);

    /** @brief add the diagonal of the Schur complemented blocks belonging to bundle and local constraints (used for diagonal preconditioning)
    */
    virtual int add_BCSchur_diagonal(CH_Matrix_Classes::Matrix& diagonal,
      MinorantBundle& globalbundle,
      CH_Matrix_Classes::Integer startindex_bundle);

    /** @brief append to lowrank "large" columns that should serve well for generating a low rank projection of the Schur complemented model part. For each column i the coordinate sigma_guess(i) gives the Diag_inv-norm for this column. The parameter minval asks to ignore columns whose norms are smaller than minval. If diaginvval is positive, the vector Diag_inv is this value times the all ones vector.

      On input lowrank must have the correct number of rows already but may
      have 0 columns.
    */
    virtual int propose_BCSchur_pcsubspace(CH_Matrix_Classes::Matrix& lowrank,
      CH_Matrix_Classes::Matrix& sigma_guess,
      const CH_Matrix_Classes::Matrix& Diag_inv,
      CH_Matrix_Classes::Real minval,
      CH_Matrix_Classes::Real diaginvval,
      MinorantBundle& globalbundle,
      CH_Matrix_Classes::Integer startindex_bundle);

    /** @brief compute the preconditioning low-rank representation of the Schur complementd blocks belonging to bundle and local constraints by adding a Johnson-Lindenstrauss projection onto the given subspace to glob_lowrank

    */
    virtual int prepare_BCSchur_JLprecond(CH_Matrix_Classes::Matrix& glob_lowrank,
      CH_Matrix_Classes::Matrix& subspace,
      bool append_globtransp_times_mat_to_subspace,
      MinorantBundle& globalbundle,
      CH_Matrix_Classes::Integer startindex_bundle);


    /// add the contributions to glob_rhs of the Schur complemented model block, and return local_rhs of the non complemented constraint block in the rows/columns/diagonal block starting at startindex_constraints
    virtual int add_Schur_rhs(CH_Matrix_Classes::Matrix& glob_rhs,
      CH_Matrix_Classes::Matrix* local_rhs,
      CH_Matrix_Classes::Real rhsmu,
      CH_Matrix_Classes::Real rhscorr,
      CH_Matrix_Classes::Integer startindex_constraints,
      MinorantBundle& globalbundle,
      CH_Matrix_Classes::Integer startindex_bundle);

    /// multiply in_vec with the local contribution to the global main block and add it to out_vec; the other local multiplications are carried out externally with the information provide in prepare_Schur_precond and are not done here.
    virtual int add_Schur_mult(const CH_Matrix_Classes::Matrix& in_vec,
      CH_Matrix_Classes::Matrix& out_vec,
      const CH_Matrix_Classes::Matrix* in_cvec,
      CH_Matrix_Classes::Matrix* out_cvec,
      CH_Matrix_Classes::Integer startindex_constraints,
      MinorantBundle& globalbundle,
      CH_Matrix_Classes::Integer startindex_bundle);

    /// use the computed step information to also compute the steps of the complemented parts
    virtual int computed_Schur_step(const CH_Matrix_Classes::Matrix& xstep,
      const CH_Matrix_Classes::Matrix& local_step,
      CH_Matrix_Classes::Integer startindex_model,
      MinorantBundle& globalbundle,
      CH_Matrix_Classes::Integer startindex_bundle);

    ///output settings
    void set_cbout(const CBout* cb, int incr = -1) {
      CBout::set_cbout(cb, incr);
      for (unsigned int i = 0; i < block.size(); i++)
        block[i]->set_cbout(this);
    }



  };


  //@}

}

#endif

