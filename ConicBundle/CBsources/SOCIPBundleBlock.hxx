/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/SOCIPBundleBlock.hxx
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



#ifndef CONICBUNDLE_SOCIPBUNDLEBLOCK_HXX
#define CONICBUNDLE_SOCIPBUNDLEBLOCK_HXX

/**  @file SOCIPBundleBlock.hxx
    @brief Header declaring the class ConicBundle::SOCIPBundleBlock
    @version 1.0
    @date 2020-05-02
    @author Christoph Helmberg
*/

#include "InteriorPointBundleBlock.hxx"
#include "SOCIPBlock.hxx"

namespace ConicBundle {

  /** @ingroup ConstrainedQPSolver
   */

   //@{

   /** @brief  interior point variables and routines specific to primal dual complementarity conditions of a second order cone with special routines for handling the bundle and the trace constraint

       adds bundle functionality to SOCIPBlock
   */

  class SOCIPBundleBlock : public InteriorPointBundleBlock, public SOCIPBlock {
  private:
    CH_Matrix_Classes::Matrix Bt;       ///< matrix representation of bundle
    CH_Matrix_Classes::Matrix Boffset; ///< columnvector of bundle offsets
    CH_Matrix_Classes::Matrix sqrBnorms; ///< if rowdim>0 it holds the squared norm of each subgradient in the same sequence as in B

  public:
    /// reset all point information to zero for dimension dim, the rest to zero
    void clear(CH_Matrix_Classes::Integer dim = 0);

    /// default constructor, allows to intialize the dimension
    SOCIPBundleBlock(CH_Matrix_Classes::Integer dim = 0, CBout* cb = 0, int cbinc = -1);

    /// destructor
    ~SOCIPBundleBlock();

    /// returns a clone
    virtual InteriorPointBundleBlock* clone();

    /// copies to content of the argument to this; to work *this must be a clone of the argument; sofar this is only needed for comparative testing
    virtual int copy_from(InteriorPointBundleBlock*);

    /// returns the number of consecutive bundle elements this cone makes use of
    virtual CH_Matrix_Classes::Integer dim_bundle() const {
      return bundle_dim;
    }

    /// return the "trace" value of the current point
    virtual CH_Matrix_Classes::Real evaluate_trace_x();

    /// return the "trace" value of the current point
    virtual CH_Matrix_Classes::Real evaluate_trace_z();

    /// return the change in "trace" value caused by the current step
    virtual CH_Matrix_Classes::Real evaluate_trace_dx();

    /// add alpha*trace_vec to vec starting at startindex
    virtual int add_trace(CH_Matrix_Classes::Matrix& vec,
      CH_Matrix_Classes::Real alpha,
      CH_Matrix_Classes::Integer startindex);

    /// set the trace premultiplied by sqrt(inv(xiz)) in vec[startindex+0,...,startindex+dim_bundle()-1]
    virtual int set_xizinvsqrt_trace(CH_Matrix_Classes::Matrix& vec,
      CH_Matrix_Classes::Integer startindex);

    // after the bundle subproblem is solved, this retrieves the local linear solution vector; if linx_activity is set, the values between zero and one indicate the guess on the coefficients activity level 
    // virtual int get_linx(CH_Matrix_Classes::Matrix& linx,
    // 		       CH_Matrix_Classes::Matrix* linx_activity=0) const;

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

    /// adds opB transposed times modelx (without constant affine term) to the arguments
    virtual int add_modelx_aggregate(CH_Matrix_Classes::Real& val,
      CH_Matrix_Classes::Matrix& vec,
      MinorantBundle& global_bundle,
      CH_Matrix_Classes::Integer startindex_bundle);

    /// set the model violation for the current system solution for the precomputed rhs on basis of the y and tracedual set in connection with set_bundle_z()/add_trace_to_diff_model or do_bundle_step()
    virtual int get_sysviol_model(CH_Matrix_Classes::Matrix& sysviol_model,
      CH_Matrix_Classes::Integer startindex_model,
      const CH_Matrix_Classes::Matrix& dy,
      const CH_Matrix_Classes::Real deltatrdual,
      MinorantBundle& global_bundle,
      CH_Matrix_Classes::Integer startindex_bundle);

    /// set z to the slack of the bundle and return a value>=0 that needs to be added to make it feasible
    virtual int set_bundle_z(const CH_Matrix_Classes::Matrix& y,
      MinorantBundle& global_bundle,
      CH_Matrix_Classes::Integer startindex_bundle,
      CH_Matrix_Classes::Real& add_center_value);

    ///add trace_dual*trace to diff_model for the right hand side (negative of current model violation)
    virtual int add_trace_to_diff_model(CH_Matrix_Classes::Real trace_dual);

    ///return the squared Euclidean norm of the dual model violation  
    virtual CH_Matrix_Classes::Real dualviol_2normsqr();

    /// move to (x+alpha*dx, z+alpha*dz), update diff_model and possibly reduce the model size if some part is too small relative to trace_rhs
    virtual int do_bundle_step(CH_Matrix_Classes::Real alpha,
      const CH_Matrix_Classes::Matrix& y,
      MinorantBundle& globalbundle,
      CH_Matrix_Classes::Integer startindex_bundle,
      CH_Matrix_Classes::Real tracedual,
      CH_Matrix_Classes::Real trace_rhs);

    /// If mu is not zero, always add the centering term for this mu as well;
    virtual int set_modelrhs(CH_Matrix_Classes::Matrix& globalrhs,
      CH_Matrix_Classes::Real rhsmu,
      CH_Matrix_Classes::Real rhscorr,
      CH_Matrix_Classes::Integer startindex_model);

    ///add the "scaled" minorant outer products to globalsys, where the correct minroants start at the given index
    virtual int add_BtinvsysB(CH_Matrix_Classes::Symmatrix& globalsys,
      const MinorantBundle& bundle,
      CH_Matrix_Classes::Integer startindex_bundle);

    /** @brief provides basic information for low rank preconditioning (in the extreme case solving) for the case of Schur complementing the model and the trace into the quadratic block

        @param[in] glob_lowrank
           set the low rank contribution Bt*invsys^(.5)
           without considering the correction for the trace vector
     for the columns startindex_bundle ... startindex_bundle+vecdim()-1

        @param[in] globalbundle
           only input for forming Bt if needed

        @param[in] startindex_bundle
           input for forming Bt together with globalbundle and column
           positions within glob_lowrank

        @param[in] trafotrace
           put invsys^(.5)*trace in coordinates
           startindex_trace ... startindex_trace+vecdim()-1

        @param[in] startindex_trace
           beginning of local trace within trafotracevec

    */
    int Schur_transform_bundle(CH_Matrix_Classes::Matrix& glob_lowrank,
      MinorantBundle& globalbundle,
      CH_Matrix_Classes::Integer startindex_bundle,
      CH_Matrix_Classes::Matrix& trafotrace,
      CH_Matrix_Classes::Integer startindex_trace);


    /** @brief add diag(Bt*sqrt(invsys)*(I-lambda*trvec*trvec')*sqrt(invsys)*B) to diagonal

        @param[out] diagonal
           add the diagonal entries diag(Bt*invsys*B) here

        @param[out] ipBtrvec
           add the vector Bt*sqrt(invsys)*trvec here

        @param[in] globalbundle
           the bundle vectors are [startindex_bundle ... startindex_bundle+vecdim()-1]

        @param[in] startindex_bundle
           see globalbundle

        @param[in] trafotrace
           holds precomputed invsys^(.5)*trace in coordinates
           [startindex_trace ... startindex_trace+vecdim()-1]

        @param[in] startindex_trace
           beginning of local trace within trafotracevec


     */
    virtual int add_bundle_xizinv_diagonal(CH_Matrix_Classes::Matrix& diagonal,
      CH_Matrix_Classes::Matrix& ipBtrvec,
      MinorantBundle& globalbundle,
      CH_Matrix_Classes::Integer startindex_bundle,
      const CH_Matrix_Classes::Matrix& trafotrace,
      CH_Matrix_Classes::Integer startindex_trace);


    /** @brief append to lowrank "large" columns that should serve well
        for generating a low rank projection of the Schur complemented
        model part. For each column i the coordinate sigma_guess(i)
        gives the Diag_inv-norm for this column. The parameter minval
        asks to ignore columns whose norms are smaller than minval. If
        diaginvval is positive, the vector Diag_inv is this value times
        the all ones vector.

        On input lowrank must have the correct number of rows already but may
        have 0 columns.

        For most new columns of lowrank it will be necessary to subtract
        the same tracevector-column according to the inner product of
        the column stored in minus_trmult, so for each new column the
        inner product of the generating vector with tracevector has to
        be appended to minus_trmult.
    */
    virtual int add_pcsubspace(CH_Matrix_Classes::Matrix& lowrank,
      CH_Matrix_Classes::Matrix& sigma_guess,
      const CH_Matrix_Classes::Matrix& Diag_inv,
      CH_Matrix_Classes::Real minval,
      CH_Matrix_Classes::Real diaginvval,
      CH_Matrix_Classes::Matrix& minus_trmult,
      CH_Matrix_Classes::Real schur_trace,
      MinorantBundle& globalbundle,
      CH_Matrix_Classes::Integer startindex_bundle);


    /// add bundle*sqrt(inv(xiz))*subspace to glob_lowrank with bundle(:,si_bundle+1:si_bundle+dim_bundle()-1) and subspace(si_subsp:si_subsp+dim_bundle,:); sqrt(inv(xiz)) has to match that used in set_xizinvsqrt_trace()
    int add_bundle_xizinvsqrt_projection(CH_Matrix_Classes::Matrix& glob_lowrank,
      CH_Matrix_Classes::Matrix& subspace,
      CH_Matrix_Classes::Integer startindex_subsspace,
      MinorantBundle& globalbundle,
      CH_Matrix_Classes::Integer startindex_bundle);

    /// out_vec+=BtinvsysB*in_vec
    virtual int add_BtinvsysB_times(const CH_Matrix_Classes::Matrix& in_vec,
      CH_Matrix_Classes::Matrix& out_vec,
      CH_Matrix_Classes::Real zeta_inval,
      CH_Matrix_Classes::Real* zeta_outval,
      MinorantBundle& globalbundle,
      CH_Matrix_Classes::Integer startindex_bundle);

    /// compute dx (and then dz) given step_y and step_trdual on basis of the last rhs computed for the model block
    virtual int set_dx_xizsolvestep(const CH_Matrix_Classes::Matrix& step_y,
      const CH_Matrix_Classes::Real step_trdual,
      MinorantBundle& globalbundle,
      CH_Matrix_Classes::Integer startindex_bundle);

    /// after the bundle subproblem is solved, this retrieves the local SOCIP solution vector; if socx_activity is set, the values between zero and one indicate the guess on the coefficients activity level 
    virtual int get_socx(CH_Matrix_Classes::Matrix& socx,
      CH_Matrix_Classes::Real* socx_activity,
      CH_Matrix_Classes::Real trace_rhs,
      bool cautious) const;



  };


  //@}

}

#endif

