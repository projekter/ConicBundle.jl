/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/QPSumModelBlock.hxx
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


#ifndef CONICBUNDLE_QPSUMMODELBLOCK_HXX
#define CONICBUNDLE_QPSUMMODELBLOCK_HXX

/**  @file QPSumModelBlock.hxx
    @brief Header declaring the class ConicBundle::QPSumModelBlock
    @version 1.0
    @date 2019-05-22
    @author Christoph Helmberg
*/


#include "QPModelBlock.hxx"

namespace ConicBundle {


/** @ingroup ConstrainedQPSolver 
 */

//@{

/** @brief  implements a (virtual) cutting model being built of a (possibly recursive) sum of QPModelBlock cutting model instances for QPSolver 

   Note, the bundle at the final level might not be the same as the
   one in the recursive calls, because it might have undergone an
   AffineFunctionTransformation.

   Not much is happening here besides passing on the calls
   to the various submodels with information on where to find
   its own (possibly transformed) bundle information and where
   to store the requested information in respective global objects
*/
  
class QPSumModelBlock: public virtual QPModelBlock, public virtual QPSumModelDataObject
{
private:
  /// the container pointing to the QPModelBlock instances comprised in the sum 
  std::vector<QPModelBlock*> blocks;
  
  CH_Matrix_Classes::Integer modeldim;    ///< the joint dimension of the modle,  may change after any change of the point
  CH_Matrix_Classes::Integer constrdim;   ///< the joint number of constraints within the models, this may change after any change of the point


public:
  /// reset to "empty/no" model
  void clear();

  /// default constructor
  QPSumModelBlock(CBout* cb=0,int cbincr=-1);

  /// virtual destructor
  virtual ~QPSumModelBlock();

  /// return a cloned object on the heap
  virtual QPModelBlockObject* clone();

  /// usually the objects of the recursive block structure and not deleted in a clear. If needed, this can be invoked explicitly here, e.g., in order to clean up clones
  virtual void recursive_delete_and_clear();

  /// sofar this is only needed for some comparative evaluations
  virtual int recursive_copy_data_of(QPModelBlockObject*);

  /// add another model at the end of the list
  int append(QPModelDataObject* inblock);

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

  /// adds opB transposed times modelx (without constant affine term) to the arguments
  virtual int add_modelx_aggregate(CH_Matrix_Classes::Real& val,
				   CH_Matrix_Classes::Matrix& vec,
				   MinorantBundle& global_bundle,
				   CH_Matrix_Classes::Integer startindex_bundle);

  /// set the model violation for the current system solution 
  virtual int get_sysviol_model(CH_Matrix_Classes::Matrix& modelvec,
				CH_Matrix_Classes::Integer startindex_model,
				const CH_Matrix_Classes::Matrix& dy,
				MinorantBundle& global_bundle,
				CH_Matrix_Classes::Integer startindex_bundle);

  /// set the constraint violation for the current system solution starting at this index
  virtual int get_sysviol_constraints(CH_Matrix_Classes::Matrix& constrvec,
				      CH_Matrix_Classes::Integer startindex_constr);

  /// output for debbuging purposes 
  virtual void display_model_values(const CH_Matrix_Classes::Matrix& y,
				    MinorantBundle& global_bundle,
				    CH_Matrix_Classes::Integer startindex_bundle,
				    std::ostream& out);

  /// reset the starting point for this value of the design variables y 
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

  ///add the "scaled" minorant outer products to globalsys, where the correct minroants start at the given index
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

  /// for passing on the solution information after solvin the system
  virtual int computed_Schur_step(const CH_Matrix_Classes::Matrix& xstep,
				  const CH_Matrix_Classes::Matrix& local_step,
				  CH_Matrix_Classes::Integer startindex_model,
				  MinorantBundle& globalbundle,
				  CH_Matrix_Classes::Integer startindex_bundle);
  


  //-----------  QPBlock routines

  /// returns the dimension of the model set (here the same as the bundle size)
  virtual CH_Matrix_Classes::Integer dim_model();

  /// returns the dimension of the system describing the model set (may contain further constraints)
  virtual CH_Matrix_Classes::Integer dim_constraints();

  ///returns the dual upper bound to the model value (the trace weighted sum of the dual trace variables); it returns 0. if no model is contained
  virtual CH_Matrix_Classes::Real constraints_cost(void);
      
  ///return the squared Euclidean norm of constraint violation of modelx
  virtual CH_Matrix_Classes::Real primalviol_2normsqr();
   
  ///return the squared Euclidean norm of the dual model violation  
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

    
  ///add the "scaled" minorant outer products to globalsys, 
  //  virtual int add_BtinvsysB(CH_Matrix_Classes::Symmatrix& globalsys);

  virtual int add_localsys(CH_Matrix_Classes::Symmatrix& globalsys,
			   CH_Matrix_Classes::Integer startindex_model,
			   CH_Matrix_Classes::Integer startindex_constraints);

 /** @brief multiply the local system diagonal block consisting of the model and local contraints rows and columns by in_vec[startindex_model+0,...,+dim_model(),startindex_constraints+0,...,+dim_constraints]  into the same coordinates of out_vec. */ 
  virtual int localsys_mult(const CH_Matrix_Classes::Matrix& in_vec,
			    CH_Matrix_Classes::Matrix& out_vec,
			    CH_Matrix_Classes::Integer startindex_model,
			    CH_Matrix_Classes::Integer startindex_constraints);
  
    
};


  //@}

}

#endif

