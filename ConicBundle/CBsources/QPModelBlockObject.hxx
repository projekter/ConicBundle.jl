/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/QPModelBlockObject.hxx
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


#ifndef CONICBUNDLE_QPMODELBLOCKOBJECT_HXX
#define CONICBUNDLE_QPMODELBLOCKOBJECT_HXX

/**  @file QPModelBlockObject.hxx
    @brief Header declaring the classes ConicBundle::QPModelBlockObject
    @version 1.0
    @date 2019-05-22
    @author Christoph Helmberg
*/


#include "symmat.hxx"

namespace ConicBundle {


  /** @ingroup ConstrainedQPSolver
   */

   //@{

   /** @brief abstract interface for model blocks in the constrained QPSolver

       Besides basic routines for retrieving sizes and multiplying with
       the bundle matrix several components of the constrained QPSolver
       call special routines, whose use within the corresponding setting
       is sketched there.

       Routines needed within
       - QPSolverBasicStructures

         + QPModelBlockObject::get_mu_info() for selecting the barrier parameter

         + QPModelBlockObject::linesearch() for performing the line search of the model variables

         + QPModelBlockObject::do_step() for carrying out the current step with the given step length

         + QPModelBlockObject::globalx_cost() for retrieving the cost contribution of the model for the current design variable values

         + QPModelBlockObject::constraints_cost() for retrieving the cost contribution of the model constraints

         + QPModelBlockObject::add_Bt_modelx() for adding the Bundle information times model variables (= the aggregate) to the dual feasibility constraint

         + QPModelBlockObject::dualviol_2normsqr() for retrieving the Euclidean norm of the model's dual violation

         + QPModelBlockObject::primalviol_2normsqr() for retrieving the Euclidean norm of the model's primal violation

         + QPModelBlockObject::reset_starting_point() for initializing the model variables to a feasible starting point with respect to the given design variable values

       - QPDirectKKTSolver

         + QPModelBlockObject::add_BDBt() in forming the low rank Schur complement of H

         + QPModelBlockObject::add_localsys() in adding the barrier diagonal of the
           B block to the the Schur_complement

         + QPModelBlockObject::add_localrhs() for constructing the right hand side

         + QPModelBlockObject::computed_step() for communicating the solution step to the model

       - QPIterativeKKTSolver

         + QPModelBlockObject::add_localrhs() for constructing the right hand side

         + QPModelBlockObject::computed_step() for communicating the solution step to the model

         + QPModelBlockObject::localsys_mult() for the multiplication with the B and C blocks
           on the main diagonal

       - QPIterativeKKTHASolver

         + QPModelBlockObject::add_Schur_rhs() for adding the right hand side contribution to the H block as induced by the Schur complement

         + QPModelBlockObject::compute_step() for computing the model solution step given the solution for the H block

         + QPModelBlockObject::add_Schur_mult() for adding the multiplication with the Schur complemented B and C blocks

       - QPKKTSubspaceHPrecond

         + QPModelBlockObject::prepare_BCSchur_JLprecond() for adding the right hand side contribution to the H block as induced by the Schur complement

     QPModelBlock combines the QPModelBlockObject with the properties of
     QPModelDataObject for obtaining the model data from a
     BundleModel. Actual implementations are QPSumModelBlock (for
     managing the sum of several Model Blocks) and QPConeModelBlock.
   */

  class QPModelBlockObject {
  public:

    virtual ~QPModelBlockObject();

    /// return a cloned object on the heap; sofar this is only needed for some comparative evaluations
    virtual QPModelBlockObject* clone() = 0;

    /// sofar this is only needed for some comparative evaluations; usually the objects of the recursive block structure and not deleted in a clear. If needed, this can be invoked explicitly here, e.g., in order to clean up clones
    virtual void recursive_delete_and_clear() = 0;

    /// sofar this is only needed for some comparative evaluations; to work *this must be a clone of the the argument and the content of the objects is copied recursively
    virtual int recursive_copy_data_of(QPModelBlockObject*) = 0;

    /// returns the dimension of the model set
    virtual CH_Matrix_Classes::Integer dim_model() = 0;

    /// returns the dimension of the system describing the model set (may contain further constraints)
    virtual CH_Matrix_Classes::Integer dim_constraints() = 0;

    /// returns the value of constant offset plus global linear cost term for the current globalx
    virtual CH_Matrix_Classes::Real globalx_cost(const CH_Matrix_Classes::Matrix& globalx) = 0;

    ///returns the dual upper bound to the model value (the trace weighted sum of the dual trace variables)
    virtual CH_Matrix_Classes::Real constraints_cost() = 0;

    /// adds opB transposed times modelx and constant affine term to the arguments
    virtual int add_Bt_modelx(CH_Matrix_Classes::Real& val,
      CH_Matrix_Classes::Matrix& vec) = 0;



    /// computes and returns C=alpha*B*A+beta*C where B and A may be transposed; C needs to have the correct size on input but will be initialized if beta==0.
    virtual CH_Matrix_Classes::Matrix& B_times(const CH_Matrix_Classes::Matrix& A,
      CH_Matrix_Classes::Matrix& C,
      CH_Matrix_Classes::Real alpha = 1.,
      CH_Matrix_Classes::Real beta = 0.,
      int Btrans = 0,
      int Atrans = 0) = 0;

    /// computes and returns C=alpha*A*B+beta*C where A and B may be transposed; C needs to have the correct size on input but will be initialized if beta==0.
    virtual CH_Matrix_Classes::Matrix& times_B(const CH_Matrix_Classes::Matrix& A,
      CH_Matrix_Classes::Matrix& C,
      CH_Matrix_Classes::Real alpha = 1.,
      CH_Matrix_Classes::Real beta = 0.,
      int Atrans = 0,
      int Btrans = 0) = 0;

    /// add B*Diag(diagvec)*Bt to C  in the principal block starting at startindex
    virtual CH_Matrix_Classes::Symmatrix& add_BDBt(const CH_Matrix_Classes::Matrix& diagvec,
      CH_Matrix_Classes::Symmatrix& S,
      bool minus = false,
      CH_Matrix_Classes::Integer startindex = 0) = 0;

    /// get the current matrix for the coupling matrix Bt in the first row; Bt needs to have the correct size on input (the size of globalx may not be known in the block) but will be initialized
    virtual CH_Matrix_Classes::Matrix& get_Bt(CH_Matrix_Classes::Matrix& Bt,
      CH_Matrix_Classes::Integer start_col = 0) = 0;


    /// initialize the model variables to a strictly feasible "central" starting point; this is the first call when the next QP problem is solved, so other initialization steps may be appropriate as well here. 
    virtual int reset_starting_point(const CH_Matrix_Classes::Matrix& firsty, CH_Matrix_Classes::Real mu) = 0;

    ///Euclidean norm of constraint violation of modelx
    virtual CH_Matrix_Classes::Real primalviol_2normsqr() = 0;

    ///Euclidean norm of constraint violation on the dual model side for the current point transmitted in the latest call to starting_point() or  do_step() respectively  
    virtual CH_Matrix_Classes::Real dualviol_2normsqr() = 0;

    /// add dimensions of the primal-dual pairs to mudim and add the "trace" (the inner product with center) of the respective primal-dual pair products for the current step; update the min and max values of x_i*z_i
    virtual int get_mu_info(CH_Matrix_Classes::Integer& mudim,
      CH_Matrix_Classes::Real& tr_xz,
      CH_Matrix_Classes::Real& tr_xdzpdxz,
      CH_Matrix_Classes::Real& tr_dxdz,
      CH_Matrix_Classes::Real& min_xz,
      CH_Matrix_Classes::Real& max_xz) const = 0;

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
      CH_Matrix_Classes::Real& ip_dxdz_xdzpdxz) const = 0;

    /// if necessary, reduce alpha to the biggest value so that feasibility is maintained with this step size
    virtual int linesearch(CH_Matrix_Classes::Real& alpha) const = 0;

    /// compute the step in the model space given the step in the design space
    virtual int compute_step(const CH_Matrix_Classes::Matrix& ystep) = 0;

    /// store the computed step and compute the missing dual step information
    virtual int computed_step(const CH_Matrix_Classes::Matrix& modelxstep,
      const CH_Matrix_Classes::Matrix& modelconstrstep) = 0;

    /// move in the last computed step direction by a step of length alpha
    virtual int do_step(CH_Matrix_Classes::Real alpha,
      const CH_Matrix_Classes::Matrix& nexty) = 0;


    /** @brief  compute the right hand sides for the rows correpsonding to the local model and contraints and add them to the global system accoriding to appennd. If mu is not zero, always add the centering term for this mu as well, if append is false, add the Schur complement rhs for add_BtinvsysB, if append is true, append the rhs of the local system

        On input the model positions of globalrhs may, e,g. be  filled with  -B ysysinv rhsy if the Q block is eliminated in a Schur complement approach
     */
    virtual int add_localrhs(CH_Matrix_Classes::Matrix& globalrhs,
      CH_Matrix_Classes::Real rhsmu,
      CH_Matrix_Classes::Real rhscorr,
      CH_Matrix_Classes::Integer startindex_model,
      CH_Matrix_Classes::Integer startindex_constraints,
      bool append) = 0;


    /** @brief add the Schur complement of the local system Matrix to the global one so that, if minus == false  the positive definite block of model is added or, if minus == true, subtracted in the leading block starting at index 0
     */
    virtual int add_BtinvsysB(CH_Matrix_Classes::Symmatrix& globalsys) = 0;

    /** @brief add the local system Matrix to the global one so that the positive definite complementarity block of model is added in the block starting at startindex_model and with its constraints blocks starting in startindex_constraints; if any startindex_value is negative, do not add the respecective part
     */
    virtual int add_localsys(CH_Matrix_Classes::Symmatrix& globalsys,
      CH_Matrix_Classes::Integer startindex_model,
      CH_Matrix_Classes::Integer startindex_constraints) = 0;

    /** @brief  given the Cholesky factorization LL' of *minus* the blocks A and B (contraints on design variables and Bundle-modelx), solve for the local constraints C and then for AB and return the solution to the ABC-variables without storing them yet (this will be done by computed_step() ).
     */
    virtual int solve_constrsys(const CH_Matrix_Classes::Symmatrix& ABchol,
      CH_Matrix_Classes::Matrix& ABCrhs_and_sol,
      CH_Matrix_Classes::Integer startindex_model,
      CH_Matrix_Classes::Integer startindex_constraints) = 0;

    /** @brief multiply the local system diagonal block consisting of the model and local contraints rows and columns by in_vec[startindex_model+0,...,+dim_model(),startindex_constraints+0,...,+dim_constraints]  into the same coordinates of out_vec. */
    virtual int localsys_mult(const CH_Matrix_Classes::Matrix& in_vec,
      CH_Matrix_Classes::Matrix& out_vec,
      CH_Matrix_Classes::Integer startindex_model,
      CH_Matrix_Classes::Integer startindex_constraints) = 0;


    /** @brief add the diagonal of the Schur complemented blocks belonging to bundle and local constraints (used for diagonal preconditioning)
    */
    virtual int add_BCSchur_diagonal(CH_Matrix_Classes::Matrix& diagonal) = 0;


    /** @brief append to lowrank "large" columns that should serve well for generating a low rank projection of the Schur complemented model part. For each column i the coordinate sigma_guess(i) gives the Diag_inv-norm for this column. The parameter minval asks to ignore columns whose norms are smaller than minval. If diaginvval is positive, the vector Diag_inv is this value times the all ones vector.

      On input lowrank must have the correct number of rows already but may
      have 0 columns.
    */
    virtual int propose_BCSchur_pcsubspace(CH_Matrix_Classes::Matrix& lowrank,
      CH_Matrix_Classes::Matrix& sigma_guess,
      const CH_Matrix_Classes::Matrix& Diag_inv,
      CH_Matrix_Classes::Real minval,
      CH_Matrix_Classes::Real diaginvval = -1.) = 0;

    /** @brief compute the preconditioning low-rank representation of the Schur complemented blocks belonging to bundle and local constraints by adding a Johnson-Lindenstrauss projection onto the given subspace to glob_lowrank; if append_glob_times_mat_to_subspace==true go the opposite direction with the transpose of glob_lowrank and append the result as new columns to subspace.
    */
    virtual int prepare_BCSchur_JLprecond(CH_Matrix_Classes::Matrix& glob_lowrank,
      CH_Matrix_Classes::Matrix& subspace,
      bool append_globtransp_times_mat_to_subspace = false) = 0;

    /// add the contributions to glob_diagonal and glob_rhs of the Schur complemented parts, and return local_rhs, local_globblock, local_diagblock of the non complemented parts 
    virtual int add_Schur_rhs(CH_Matrix_Classes::Matrix& glob_rhs,
      CH_Matrix_Classes::Matrix* local_rhs,
      CH_Matrix_Classes::Real rhsmu,
      CH_Matrix_Classes::Real rhscorr) = 0;

    /// multiply in_Qvec and in_Cvec with the local Schur complement contribution to the main block and the local constraints block, add it to out_Qvec and outCvec.
    virtual int add_Schur_mult(const CH_Matrix_Classes::Matrix& in_Qvec,
      CH_Matrix_Classes::Matrix& out_Qvec,
      const CH_Matrix_Classes::Matrix* in_Cvec = 0,
      CH_Matrix_Classes::Matrix* out_Cvec = 0) = 0;

    /// use the computed step information to also compute the steps of the complemented parts
    virtual int computed_Schur_step(const CH_Matrix_Classes::Matrix& xstep,
      const CH_Matrix_Classes::Matrix& local_step) = 0;

    //-----------  mainly for testing the correctness of the system solution

    /// get the vector formed by all model x variables
    virtual CH_Matrix_Classes::Matrix& get_x() = 0;
    /// get the vector formed by all delta model x variables
    virtual CH_Matrix_Classes::Matrix& get_dx() = 0;
    /// get the vector formed by all delta dual constraint variables
    virtual CH_Matrix_Classes::Matrix& get_dcstr() = 0;
    /// get the model violation for the current system solution with this dy
    virtual CH_Matrix_Classes::Matrix& get_sysviol_model(const CH_Matrix_Classes::Matrix& dy) = 0;
    /// get the constraint violation for the current system solution
    virtual CH_Matrix_Classes::Matrix& get_sysviol_constraints() = 0;


    /// output the model values in a readable format for testing 
    virtual void display_model_values(const CH_Matrix_Classes::Matrix& y, std::ostream& out) = 0;

  };



  //@}

}

#endif

