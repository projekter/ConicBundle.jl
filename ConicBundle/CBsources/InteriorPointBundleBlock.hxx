/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/InteriorPointBundleBlock.hxx
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



#ifndef CONICBUNDLE_INTERIORPOINTBUNDLEBLOCK_HXX
#define CONICBUNDLE_INTERIORPOINTBUNDLEBLOCK_HXX

/**  @file InteriorPointBundleBlock.hxx
    @brief Header declaring the class ConicBundle::InteriorPointBundleBlock
    @version 1.0
    @date 2020-05-02
    @author Christoph Helmberg
*/

#include "InteriorPointBlock.hxx"
#include "MinorantPointer.hxx"
#include "QPModelDataObject.hxx"

namespace ConicBundle {
  
/** @ingroup ConstrainedQPSolver 
 */

//@{

/** @brief  abstract interface for interior point routines specific to primal dual complementarity conditions of a symmetric cone  with special routines for handling the bundle and the trace constraint 

    this builds upon InteriorPointBlock which comprises all non bundle
    specific routines
*/


class InteriorPointBundleBlock: public virtual InteriorPointBlock
{
protected:
  CH_Matrix_Classes::Integer bundle_dim;   ///< dimension of the bundle (in case of subspace projections this may differ from vecdim
  CH_Matrix_Classes::Matrix diff_model;  ///< negative evaluation of minorants in current point

public:
  /// default constructor
  InteriorPointBundleBlock(CBout* cb=0,int cbinc=-1):CBout(cb,cbinc){bundle_dim=0;}

  /// virtual destructor
  virtual ~InteriorPointBundleBlock(); 
  
  /// returns a clone; sofar this is only needed for comparative testing
  virtual InteriorPointBundleBlock* clone()=0;
  
  /// copies to content of the argument to this; to work *this must be a clone of the argument; sofar this is only needed for comparative testing
  virtual int copy_from(InteriorPointBundleBlock*)=0;
  
  //virtual CH_Matrix_Classes::Integer get_vecdim() const =0;

  ///allows to pass on additional information about the oracle if required
  virtual void set_oracle_data(QPModelOracleDataObject* /*oracle_data*/)
  {}

  /// returns the number of consecutive bundle elements this cone makes use of
  virtual CH_Matrix_Classes::Integer dim_bundle() const=0;

  /// return the "trace" value of the current point
  virtual CH_Matrix_Classes::Real evaluate_trace_x() =0;

  /// return the "trace" value of the current point
  virtual CH_Matrix_Classes::Real evaluate_trace_z() =0;

  /// return the change in "trace" value caused by the current step
  virtual CH_Matrix_Classes::Real evaluate_trace_dx() =0;

  /// add alpha*trace_vec to vec starting at startindex
  virtual int add_trace(CH_Matrix_Classes::Matrix& vec,
			CH_Matrix_Classes::Real alpha,
			CH_Matrix_Classes::Integer startindex) =0;

  /// set the trace premultiplied by sqrt(inv(xiz)) in vec[startindex+0,...,startindex+dim_bundle()-1]
  virtual int set_xizinvsqrt_trace(CH_Matrix_Classes::Matrix& vec,
				   CH_Matrix_Classes::Integer startindex) =0;
  
  // set x to value*"one" to x, or if add==true, add value*"one" to x
  //virtual int center_x(CH_Matrix_Classes::Real val,bool add=false)=0;

  // set z to value*"one" to z, or if add==true, add value*"one" to z
  //virtual int center_z(CH_Matrix_Classes::Real val,bool add=false)=0;

  // set x to the values of vec[startindex+0,+1 ...,+(vecdim-1)] and return in add_center_value a value>=0 that needs to be added to make it feasible
  //virtual int set_x(const CH_Matrix_Classes::Matrix& vec,
  //		    CH_Matrix_Classes::Integer startindex,
  //		    CH_Matrix_Classes::Real& add_center_value)=0;

  // set z to the values of vec[startindex+0,+1 ...,+(vecdim-1)] and return a value>=0 that needs to be added to make it feasible
  //virtual int set_z(const CH_Matrix_Classes::Matrix& vec,
  //		    CH_Matrix_Classes::Integer startindex,
  //		    CH_Matrix_Classes::Real& add_center_value)=0;

  // on vec[startindex+0,+1 ...,+(vecdim-1)] put or add  a * x into vec for a real number a  
  //virtual int vecgetsax(CH_Matrix_Classes::Matrix& vec,
  //			CH_Matrix_Classes::Integer startindex,
  //			CH_Matrix_Classes::Real a=1.,
  //			bool add=false)=0;

  // on vec[startindex+0,+1 ...,+(vecdim-1)] put or add a * z into vec for a real number a   
  //virtual int vecgetsaz(CH_Matrix_Classes::Matrix& vec,
  //			CH_Matrix_Classes::Integer startindex,
  //			CH_Matrix_Classes::Real a=1.,
  //			bool add=false)=0;

  // add dimensions of the primal-dual pairs to mudim and add the inner products of the primal-dual pairs of the current point to current_ip and those of the next point obtained by the given stepsize with the latest computed step to step_ip
  //virtual int get_mu_info(CH_Matrix_Classes::Integer& inmudim,
  //		  CH_Matrix_Classes::Real& current_ip,
  //		  CH_Matrix_Classes::Real& step_ip,
  //		  CH_Matrix_Classes::Real stepsize) const =0;

  // if necessary, reduce alpha to the biggest value so that feasibility is maintained with this step size
  //virtual int linesearch(CH_Matrix_Classes::Real& alpha) const=0;

  // compute the complementarity_rhs=rhsmu*xi-rhscorr*xi*dx*dz (wihtout "-z") for mu=rhsmu and for corrector for factor rhscorr>0., store this and add it to rhs;
  //virtual int add_muxinv(CH_Matrix_Classes::Matrix& rhs,
  //			 CH_Matrix_Classes::Integer startindex,
  //			 CH_Matrix_Classes::Real rhsmu,
  //			 CH_Matrix_Classes::Real rhscorr,
  //			 bool minus=false)=0;

    /// extract dx from rhs at startindex and compute at the same time dz (=-sys dx-z  +complentarity_rhs); this may only be called after add_muxinv() was called for this point
  //virtual int set_dx(const CH_Matrix_Classes::Matrix& rhs,CH_Matrix_Classes::Integer startindex)=0;

  // compute dx=sysinv*rhs and at the same time dz (=-rhs-z +complentarity_rhs); may only be called after add_muxinv was called for the most recent point; this may only be called after add_muxinv() was called for this point
  //virtual int set_dx_xizsolverhs(const CH_Matrix_Classes::Matrix& rhs,
  //				CH_Matrix_Classes::Integer startindex)=0;

  // compute sysinv*rhs into rhs, possibly with a negative sign 
  //virtual int apply_xizinv(CH_Matrix_Classes::Matrix& rhs,
  //			   CH_Matrix_Classes::Integer startindex,
  //			   bool minus=false)=0;

  // compute sys*rhs into rhs, possibly with a negative sign
  //virtual int apply_xiz(CH_Matrix_Classes::Matrix& rhs,
  //			CH_Matrix_Classes::Integer startindex,
  //			bool minus=false)=0;

  // move to (x+alpha*dx, z+alpha*dz)
  //virtual int do_step(CH_Matrix_Classes::Real alpha)=0;
  //

  // add the Schur complement to a big system matrix
  //virtual int add_AxizinvAt(const CH_Matrix_Classes::Matrix& A,
  //			    CH_Matrix_Classes::Symmatrix& globalsys,
  //			    bool minus=false)=0;

  //add (or subract if minus==true) the system matrix to a big system matrix starting at startindex
  //virtual int add_xiz(CH_Matrix_Classes::Symmatrix& globalsys,
  //		      CH_Matrix_Classes::Integer startindex,
  //		      bool minus=false)=0;

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
	  CH_Matrix_Classes::Integer startindex_bundle)=0;
						   
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
	  CH_Matrix_Classes::Integer startindex_bundle)=0;

  ///add the main diagonal block tranpose(projection)*diagvec*projection to bigS starting at startindex
  virtual CH_Matrix_Classes::Symmatrix&
  add_BDBt(const CH_Matrix_Classes::Matrix& diagvec,
	   CH_Matrix_Classes::Symmatrix& bigS,
	   bool minus,
	   CH_Matrix_Classes::Integer startindex,
	   CH_Matrix_Classes::Matrix& Bt,
	   CH_Matrix_Classes::Integer startindex_model,
	   MinorantBundle& globalbundle,
	   CH_Matrix_Classes::Integer startindex_bundle)=0;

  /// get the current matrix for the coupling matrix Bt in the first row of blocks
  virtual CH_Matrix_Classes::Matrix&
  get_Bt(CH_Matrix_Classes::Matrix& Bt,
	 CH_Matrix_Classes::Integer startindex_model,
	 MinorantBundle& global_bundle,
	 CH_Matrix_Classes::Integer startindex_bundle)=0;
  
  /// adds opB transposed times modelx (without constant affine term) to the arguments
  virtual int add_modelx_aggregate(CH_Matrix_Classes::Real& val,
				   CH_Matrix_Classes::Matrix& vec,
				   MinorantBundle& global_bundle,
				   CH_Matrix_Classes::Integer startindex_bundle)=0;

  /// set the model violation for the current system solution for the precomputed rhs on basis of the y and tracedual set in connection with set_bundle_z()/add_trace_to_diff_model or do_bundle_step()
  virtual int get_sysviol_model(CH_Matrix_Classes::Matrix& sysviol_model,
				CH_Matrix_Classes::Integer startindex_model,
				const CH_Matrix_Classes::Matrix& dy,
				const CH_Matrix_Classes::Real deltatrdual,
				MinorantBundle& global_bundle,
				CH_Matrix_Classes::Integer startindex_bundle)=0;

  /// set z to the slack of the bundle and return a value>=0 that needs to be added to make it feasible
  virtual int set_bundle_z(const CH_Matrix_Classes::Matrix& y,
		    MinorantBundle& global_bundle,
		    CH_Matrix_Classes::Integer startindex_bundle,
		    CH_Matrix_Classes::Real& add_center_value)=0;

  ///add trace_dual*trace to diff_model for the right hand side (negative of current model violation)
  virtual int add_trace_to_diff_model(CH_Matrix_Classes::Real trace_dual)=0;

  ///return the squared Euclidean norm of the dual model violation  
  virtual CH_Matrix_Classes::Real dualviol_2normsqr()=0;

  /// move to (x+alpha*dx, z+alpha*dz), update diff_model and possibly reduce the model size if some part is too small relative to trace_rhs
  virtual int do_bundle_step(CH_Matrix_Classes::Real alpha,
		      const CH_Matrix_Classes::Matrix& y,
		      MinorantBundle& globalbundle,
		      CH_Matrix_Classes::Integer startindex_bundle,
		      CH_Matrix_Classes::Real tracedual,
		      CH_Matrix_Classes::Real trace_rhs)=0;

  /// If mu is not zero, always add the centering term for this mu as well;
  virtual int set_modelrhs(CH_Matrix_Classes::Matrix& globalrhs, 
			   CH_Matrix_Classes::Real rhsmu,
			   CH_Matrix_Classes::Real rhscorr,
			   CH_Matrix_Classes::Integer startindex_model)=0;

  ///add the "scaled" minorant outer products to globalsys, where the correct minroants start at the given index
  virtual int add_BtinvsysB(CH_Matrix_Classes::Symmatrix& globalsys,
			    const MinorantBundle& bundle,
			    CH_Matrix_Classes::Integer startindex_bundle) =0;

  /** @brief provides basic information for low rank preconditioning (in the extreme case solving) for the case of Schur complementing the model and the trace into the quadratic block
   
      @param[out] glob_lowrank
         set the low rank contribution Bt*invsys^(.5) 
         without considering the correction for the trace vector  
	 for the columns startindex_bundle ... startindex_bundle+vecdim()-1

      @param[in] globalbundle
         only input for forming Bt if needed
       
      @param[in] startindex_bundle
         input for forming Bt together with globalbundle and column 
         positions within glob_lowrank

      @param[out] trafotrace
         put invsys^(.5)*trace in coordinates 
         startindex_trace ... startindex_trace+vecdim()-1

      @param[in] startindex_trace
         beginning of local trace within trafotracevec
        
  */
  virtual int Schur_transform_bundle(CH_Matrix_Classes::Matrix& glob_lowrank,
				     MinorantBundle& globalbundle,
				     CH_Matrix_Classes::Integer startindex_bundle,
				     CH_Matrix_Classes::Matrix& trafotrace,
				     CH_Matrix_Classes::Integer startindex_trace)=0;


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
					 CH_Matrix_Classes::Integer startindex_trace)=0;
					 


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
			     CH_Matrix_Classes::Matrix & minus_trmult,
			     CH_Matrix_Classes::Real schur_trace,
			     MinorantBundle& globalbundle,
			     CH_Matrix_Classes::Integer startindex_bundle) =0;

  /// add bundle*sqrt(inv(xiz))*subspace to glob_lowrank with bundle(:,si_bundle+1:si_bundle+dim_bundle()-1) and subspace(si_subsp:si_subsp+dim_bundle,:); sqrt(inv(xiz)) has to match that used in set_xizinvsqrt_trace(); if startindex_subspace is negative, append transpose(glob_lowrank)*bundle*sqrt(inv(xiz)) as new columns to subspace
  virtual int add_bundle_xizinvsqrt_projection(CH_Matrix_Classes::Matrix& glob_lowrank,
					       CH_Matrix_Classes::Matrix& subspace,
					       CH_Matrix_Classes::Integer startindex_subsspace,
					       MinorantBundle& globalbundle,
					       CH_Matrix_Classes::Integer startindex_bundle) =0;


  /// out_vec+=BtinvsysB*in_vec
  virtual int add_BtinvsysB_times(const CH_Matrix_Classes::Matrix& in_vec,
				  CH_Matrix_Classes::Matrix& out_vec,
				  CH_Matrix_Classes::Real zeta_inval,
				  CH_Matrix_Classes::Real* zeta_outval,
				  MinorantBundle& globalbundle,
				  CH_Matrix_Classes::Integer startindex_bundle)=0;

  /// compute dx (and then dz) given step_y and step_trdual on basis of the last rhs computed for the model block
  virtual int set_dx_xizsolvestep(const CH_Matrix_Classes::Matrix& step_y,
				  const CH_Matrix_Classes::Real step_trdual,
				  MinorantBundle& globalbundle,
				  CH_Matrix_Classes::Integer startindex_bundle)=0;

  
  //---------------------- mainly for testing

  // return the vector form of x 
  //virtual int get_vecx(CH_Matrix_Classes::Matrix& vecx,CH_Matrix_Classes::Integer startindex)=0;
  
  // return the vector form of z
  //virtual int get_vecz(CH_Matrix_Classes::Matrix& vecz,CH_Matrix_Classes::Integer startindex)=0;
  
  // return the vector form of dx, 1 if not available 
  //virtual int get_vecdx(CH_Matrix_Classes::Matrix& vecdx,CH_Matrix_Classes::Integer startindex)=0;
  
  // return the vector form of dz, 1 if not available
  // virtual int get_vecdz(CH_Matrix_Classes::Matrix& vecdz,CH_Matrix_Classes::Integer startindex)=0;
  

};

  


  //@}

}

#endif

