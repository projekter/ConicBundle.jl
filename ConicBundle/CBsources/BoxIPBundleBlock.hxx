/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/BoxIPBundleBlock.hxx
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



#ifndef CONICBUNDLE_BOXIPBUNDLEBLOCK_HXX
#define CONICBUNDLE_BOXIPBUNDLEBLOCK_HXX

/**  @file BoxIPBundleBlock.hxx
    @brief Header declaring the class ConicBundle::BoxIPBundleBlock
    @version 1.0
    @date 2020-05-02
    @author Christoph Helmberg
*/

#include "InteriorPointBundleBlock.hxx"

namespace ConicBundle {

  /** @ingroup InternalQPSolver
   */

   //@{

   /** @brief  QPBundleBlock interface for the interior point routines specific to the primal dual complementarity conditions of a scaled, fulldimensional box given by vectors lb < ub. It resolves the scaling internally. There is no pure QPBlock version of this.

   */

  class BoxIPBundleBlock : public virtual InteriorPointBundleBlock {
  protected:
    CH_Matrix_Classes::Integer vecdim;  ///< size of the model vector 
    CH_Matrix_Classes::Integer scalvecdim; ///< size of the model vector with scaling
    CH_Matrix_Classes::Matrix lb; ///< lower bounds of the box, lb < ub
    CH_Matrix_Classes::Matrix ub; ///< upper bounds of the box, lb < ub
    bool use_scaling;                ///< whether scaling is to used or s=1. throughout
    CH_Matrix_Classes::Real scalub; ///<  >0 only if use_scaling==true, and if so then enforce s <= scalub  

    CH_Matrix_Classes::Matrix x;     ///< "primal" point x,   lb*s<= x <= ub*s
    CH_Matrix_Classes::Matrix lz;     ///< "dual" point lz for lower bound >=0
    CH_Matrix_Classes::Matrix uz;     ///< "dual" point uz for upper bound >=0
    CH_Matrix_Classes::Matrix dx;    ///< current step for x
    CH_Matrix_Classes::Matrix dlz;    ///< current step for lz
    CH_Matrix_Classes::Matrix duz;    ///< current step for uz

    CH_Matrix_Classes::Real s;       ///< "primal" scaling size >=0
    CH_Matrix_Classes::Real lt;       ///< "dual" to scaling size >= 0
    CH_Matrix_Classes::Real ut;       ///< "dual" to scaling size <= scalval;
    CH_Matrix_Classes::Real ds;      ///< current step for s
    CH_Matrix_Classes::Real dlt;      ///< current step for lt (for s>=0 
    CH_Matrix_Classes::Real dut;      ///< current step for ut (for s<=scaval) 

    CH_Matrix_Classes::Matrix xiz;        ///< lz(i)/(x(i)-s*lb) + uz(i)/(s*ub-x(i))
    CH_Matrix_Classes::Matrix lxinv;      ///< 1./(x(i)-s*lb)
    CH_Matrix_Classes::Matrix uxinv;      ///< 1./(s*ub-x(i))
    CH_Matrix_Classes::Real hats;       ///< [ut/(scalub-s)]+lt/s+sum(lb(i)^2*lz(i)/(x(i)-s*lb) + ub(i)^2*uz(i)/(s*ub-x(i)))
    CH_Matrix_Classes::Matrix hatb;       ///< (lb(i)*lz(i)/(x(i)-s*lb) + ub(i)*uz(i)/(s*ub-x(i))
    CH_Matrix_Classes::Matrix sqrt_xiz; ///< xiz.^{.5}
    CH_Matrix_Classes::Real Chols;      ///< sqrt(hats-Cholb'*Cholb)
    CH_Matrix_Classes::Matrix Cholinvb; ///< xiz.^{-1}.*hatb/Chols

    CH_Matrix_Classes::Matrix compl_lrhs;  ///< rhs used in solving the complementarity line for lb
    CH_Matrix_Classes::Matrix compl_urhs;  ///< rhs used in solving the complementarity line for ub
    CH_Matrix_Classes::Real compl_ltrhs; ///< rhs usd in solving the complementartiy line for t
    CH_Matrix_Classes::Real compl_utrhs; ///< rhs usd in solving the complementartiy line for t
    CH_Matrix_Classes::Real sys_srhs; ///< rhs of s in the system with uz,lz,t eliminated

    CH_Matrix_Classes::Real last_rhs_mu;  ///< the last mu used in rhs computations
    CH_Matrix_Classes::Real mu;           ///< in a step mu gets the value of last_rhs_mu
    CH_Matrix_Classes::Real old_mu;       ///< in a step old_mu gets the value of mu before this gets last_rhs_mu

    CH_Matrix_Classes::Real last_alpha;   ///< last alpha used in do_step()

    CH_Matrix_Classes::Matrix oldx;       ///< point before x
    CH_Matrix_Classes::Real olds;          ///< point before s  
    CH_Matrix_Classes::Matrix oldlz;       ///< point before lz
    CH_Matrix_Classes::Matrix olduz;       ///< point before uz
    CH_Matrix_Classes::Real oldlt;         ///< point before lt 
    CH_Matrix_Classes::Real oldut;         ///< point before ut

    CH_Matrix_Classes::Matrix Bt;       ///< matrix representation of bundle
    CH_Matrix_Classes::Matrix Boffset; ///< columnvector of bundle offsets
    CH_Matrix_Classes::Matrix sqrBnorms; ///< if rowdim>0 it holds the squared norm of each subgradient in the same sequence as in B

    mutable CH_Matrix_Classes::Matrix tmpvec;  ///< temporary vector to avoid reallocations
    mutable CH_Matrix_Classes::Matrix tmpmat;  ///< temporary matrix to avoid reallocatoins

    CH_Matrix_Classes::Indexmatrix map_to_old; ///< currently not in use (used in failed attempt to eliminate some inactive variables on the fly)

    bool test_myself_call; ///< if true, this is already called by a routine for checking a result, so don't produce infinite test loops by calling this routine again 

    /// returns in inactive, in increasing order, the indices considered inactive (small in comparison to trace_bound or going to zero faster than their dual values)
    int find_inactive_indices(CH_Matrix_Classes::Indexmatrix& inactive,
      CH_Matrix_Classes::Real trace_bound,
      bool cautious) const;

    /// clear variables that are no longer valid for the current point
    void point_changed();

    /// computes the step information for the complementarity parts 
    int compute_NTscaling(void);

    /// does what is says
    int given_dx_compute_dz(void);

  public:
    /// reset all point information to zero for dimension dim, the rest to zero
    void clear(const CH_Matrix_Classes::Matrix& lb,
      const CH_Matrix_Classes::Matrix& ub,
      bool use_scaling,
      CH_Matrix_Classes::Real scalub);

    /// if use_scaling is true, then if scalub>0 it gives the maximum scaling value allowed, while a negative value just requires nonnegative scalings. If use_scaling is false and scalub>0. then work with (lb*scalub) and (ub*scalub);
    BoxIPBundleBlock(const CH_Matrix_Classes::Matrix& lb,
      const CH_Matrix_Classes::Matrix& ub,
      bool use_scaling = false,
      CH_Matrix_Classes::Real scalub = -1.,
      CBout* cb = 0, int cbinc = -1);

    ~BoxIPBundleBlock();

    /// returns a clone
    virtual InteriorPointBundleBlock* clone();

    /// copies to content of the argument to this; to work *this must be a clone of the argument; sofar this is only needed for comparative testing
    virtual int copy_from(InteriorPointBundleBlock*);

    /// returns the length of the varialbe vector
    virtual CH_Matrix_Classes::Integer get_vecdim() const;

    /// set x to value*"one" to x, or if add==true, add value*"one" to x
    virtual int center_x(CH_Matrix_Classes::Real val, bool add = false);

    /// set z to value*"one" to z, or if add==true, add value*"one" to z
    virtual int center_z(CH_Matrix_Classes::Real val, bool add = false);

    /// the base description "set x to the values of vec[startindex+0,+1 ...,+(vecdim-1)] and return in add_center_value a value>=0 that needs to be added to make it feasible" does not fit here, because there need not be a shift that makes x feasible for the box. In this case it returns an error.
    virtual int set_x(const CH_Matrix_Classes::Matrix& vec, CH_Matrix_Classes::Integer startindex, CH_Matrix_Classes::Real& add_center_value);

    /// set z to the values of vec[startindex+0,+1 ...,+(vecdim-1)] and add sufficient center to make z feasible, return this value>=0 in added_center_value
    virtual int set_z(const CH_Matrix_Classes::Matrix& vec, CH_Matrix_Classes::Integer startindex, CH_Matrix_Classes::Real& add_center_value);

    /// on vec[startindex+0,+1 ...,+(vecdim-1)] put or add  a * x into vec for a real number a  
    virtual int vecgetsax(CH_Matrix_Classes::Matrix& vec,
      CH_Matrix_Classes::Integer startindex,
      CH_Matrix_Classes::Real a = 1.,
      bool add = false);

    /// on vec[startindex+0,+1 ...,+(vecdim-1)] put or add a * z into vec for a real number a   
    virtual int vecgetsaz(CH_Matrix_Classes::Matrix& vec,
      CH_Matrix_Classes::Integer startindex,
      CH_Matrix_Classes::Real a = 1.,
      bool add = false);

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

    /// compute the complementarity_rhs=rhsmu*xi-rhscorr*xi*dx*dz (wihtout "-z") for mu=rhsmu and for corrector for factor rhscorr>0., store this and add it to rhs 
    virtual int add_muxinv(CH_Matrix_Classes::Matrix& rhs,
      CH_Matrix_Classes::Integer startindex,
      CH_Matrix_Classes::Real rhsmu,
      CH_Matrix_Classes::Real rhscorr,
      bool minus = false);


    /// extract dx from rhs at startindex and compute at the same time dz (=-sys dx-z  +complentarity_rhs); 
    virtual int set_dx(const CH_Matrix_Classes::Matrix& rhs, CH_Matrix_Classes::Integer startindex);

    /// compute dx=sysinv*rhs and at the same time dz (=-rhs -z +complentarity_rhs); 
    virtual int set_dx_xizsolverhs(const CH_Matrix_Classes::Matrix& rhs, CH_Matrix_Classes::Integer startindex);

    /// compute sysinv*rhs into rhs, possibly with a negative sign 
    virtual int apply_xizinv(CH_Matrix_Classes::Matrix& rhs,
      CH_Matrix_Classes::Integer startindex,
      bool minus = false);

    /// compute sys*rhs into rhs, possibly with a negative sign
    virtual int apply_xiz(CH_Matrix_Classes::Matrix& rhs,
      CH_Matrix_Classes::Integer startindex,
      bool minus = false);

    /// move to (x+alpha*dx, z+alpha*dz)
    virtual int do_step(CH_Matrix_Classes::Real alpha);

    /// add the Schur complement to a big system matrix
    virtual int add_AxizinvAt(const CH_Matrix_Classes::Matrix& A,
      CH_Matrix_Classes::Symmatrix& globalsys,
      bool minus = false,
      bool Atrans = false);

    /// add (or subract if minus==true) the system matrix to a big system matrix starting at startindex
    virtual int add_xiz(CH_Matrix_Classes::Symmatrix& globalsys, CH_Matrix_Classes::Integer startindex, bool minus = false);



    //---------------------- mainly for testing

    /// return the vector form of x 
    virtual int get_vecx(CH_Matrix_Classes::Matrix& vecx, CH_Matrix_Classes::Integer startindex) {
      CH_Matrix_Classes::mat_xey(vecdim, vecx.get_store() + startindex, x.get_store());
      if (use_scaling)
        vecx(startindex + vecdim) = s;
      return 0;
    }

    /// return the vector form of z
    virtual int get_vecz(CH_Matrix_Classes::Matrix& vecz, CH_Matrix_Classes::Integer startindex) {
      CH_Matrix_Classes::mat_xey(vecdim, vecz.get_store() + startindex, lz.get_store());
      CH_Matrix_Classes::mat_xmey(vecdim, vecz.get_store() + startindex, uz.get_store());
      if (use_scaling)
        vecz(startindex + vecdim) = lt - (scalub > 0. ? ut : 0.);
      return 0;
    }

    /// return the vector form of dx, 1 if not available 
    virtual int get_vecdx(CH_Matrix_Classes::Matrix& vecdx, CH_Matrix_Classes::Integer startindex) {
      if (dx.dim() != vecdim)
        return 1;
      CH_Matrix_Classes::mat_xey(vecdim, vecdx.get_store() + startindex, dx.get_store());
      if (use_scaling)
        vecdx(startindex + vecdim) = ds;
      return 0;
    }


    /// return the vector form of dz, 1 if not available
    virtual int get_vecdz(CH_Matrix_Classes::Matrix& vecdz, CH_Matrix_Classes::Integer startindex) {
      if (dlz.dim() != vecdim)
        return 1;
      CH_Matrix_Classes::mat_xey(vecdim, vecdz.get_store() + startindex, dlz.get_store());
      CH_Matrix_Classes::mat_xmey(vecdim, vecdz.get_store() + startindex, duz.get_store());
      if (use_scaling)
        vecdz(startindex + vecdim) = dlt - (scalub > 0. ? dut : 0.);
      return 0;
    }

    virtual CH_Matrix_Classes::Integer dim_bundle() const {
      return bundle_dim;
    }

    virtual CH_Matrix_Classes::Real evaluate_trace_x();

    /// return the "trace" value of the current point
    virtual CH_Matrix_Classes::Real evaluate_trace_z();

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

        - glob_lowrank: set the low rank contribution Bt*invsys^(.5)
          without considering the correction for the trace vector
    for the columns startindex_bundle ... startindex_bundle+vecdim()-1

        - globalbundle: only input for forming Bt if needed

        - startindex_bundle: input for forming Bt together with globalbundle
          and column positions within glob_lowrank

        - trafotracevec: put invsys^(.5)*trace in coordinates
          startindex_trace ... startindex_trace+vecdim()-1

        - startindex_trace: beginning of local trace within trafotracevec

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


    /// after the bundle subproblem is solved, this retrieves the local linear solution vector; if linx_activity is set, the values between zero and one indicate the guess on the coefficients activity level 
    int get_boxx(CH_Matrix_Classes::Matrix& linx,
      CH_Matrix_Classes::Matrix* linx_activity,
      CH_Matrix_Classes::Real trace_rhs,
      bool cautious) const;

    /// returns the dual cost of the box constraints
    CH_Matrix_Classes::Real constraints_cost() {
      if (use_scaling) {
        if (scalub > 0.) {
          return scalub * ut;
        } else {
          return 0.;
        }
      }
      CH_Matrix_Classes::Real val = CH_Matrix_Classes::ip(ub, uz) - CH_Matrix_Classes::ip(lb, lz);
      return val;
    }

  };


  //@}

}

#endif

