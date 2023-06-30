/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/QPModelBlock.hxx
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


#ifndef CONICBUNDLE_QPMODELBLOCK_HXX
#define CONICBUNDLE_QPMODELBLOCK_HXX

/**  @file QPModelBlock.hxx
    @brief Header declaring the classes ConicBundle::QPModelBlock, ConicBundle::QPModelPointer
    @version 1.0
    @date 2019-05-22
    @author Christoph Helmberg
*/


#include "QPModelDataObject.hxx"
#include "QPModelBlockObject.hxx"

namespace ConicBundle {


/** @ingroup ConstrainedQPSolver 
 */

//@{

/** @brief combines and provides basic functionalities of QPModelDataObject and QPModelBlockObject, but is still abstract

    This class serves as a base class for actual implementations of the
    models and provides some basic variables and functionalities required
    for uniform use by BundleModel and QPSolver.
    
    In particular it provides storage for and access to the bundle
    (and possibly an additional constant minorant) of the underlying
    cutting model(s). It also ensures consistent handling of each 
    AffineFunctionTransformation of this data. It might be worth
    to consider to collect the transformations and not to execute
    them at once, but direct execution of these transformations
    is the current approach.

    Most routines of QPModelBlockObject are transcribed to variants
    working with indices on where to find the bundle information of
    a particular submodel, so that it makes use of the information
    with the requested AffineFunctionTransformation applied to it. 
    
    Actual implementations are QPSumModelBlock (for managing the sum
    of several model blocks) and QPConeModelBlock.
 */
  
  class QPModelBlock: public virtual QPModelDataObject, public QPModelBlockObject
{
protected:
    ///constant offset minorant (fixed affine function to be added to the model); each aft creates a new one with push_back
  std::vector<MinorantPointer> constant_minorant;
  ///the minorants forming the cutting model(s); how to combine them within the model(s) comprised in *this is described in derived classes; each aft creates a new one with push_back
  std::vector<MinorantBundle> bundle;

  CH_Matrix_Classes::Matrix modelx;               ///< the current vector of model variables of all models comprised in *this
  CH_Matrix_Classes::Matrix Bt;                   ///< if the matrix of the bundle information has to be formed at least once, it is then stored here for later use 
  CH_Matrix_Classes::Matrix modeldx;              ///< only for testing
  CH_Matrix_Classes::Matrix modeldcstr;           ///< only for testing
  CH_Matrix_Classes::Matrix sysviol_model;        ///< only for testing
  CH_Matrix_Classes::Matrix sysviol_constraints;  ///< only for testing
  
  MinorantPointer modelx_aggregate; ///< if asked to form the aggregate for the current modlex, it is also stored here for later use

  /// whenever modelx is (going to be) changed, the information collected for the oblivious modelx is reset to "empty" here
  virtual void modelx_changed()
  {
    modelx.init(0,1,0.);

    modeldx.init(0,0,0.);
    modeldcstr.init(0,0,0.);
    sysviol_model.init(0,0,0.);
    sysviol_constraints.init(0,0,0.);
  }

public:

  /// reset to uninitialized state (no model)
  void clear()
  { constant_minorant.clear(); bundle.clear();modelx_aggregate.clear();
    modelx_changed();Bt.init(0,0,0.);}

  /// default constructor
  QPModelBlock(CBout* cb=0,int cbinc=-1):QPModelDataObject(cb,cbinc),QPModelBlockObject()
  {}

  /// virtual destructor
  virtual ~QPModelBlock();

  /// usually the objects of the recursive block structure and not deleted in a clear. If needed, this can be invoked explicitly here, e.g., in order to clean up clones
  virtual void recursive_delete_and_clear()=0;

  /// return a cloned object on the heap
  virtual QPModelBlockObject* clone()=0;

    ///gives reading access to a constant offset minorant
  virtual const MinorantPointer& get_constant_minorant() const
  { return constant_minorant.back(); }
  
  ///gives reading access to the bundle minorants of the cutting model
  virtual const MinorantBundle& get_bundle() const
  { return bundle.back();}  

  ///gives access to a constant offset minorant
  virtual MinorantPointer& get_constant_minorant() 
  { return constant_minorant.back(); }
  
  ///gives access to the bundle minorants of the cutting model
  virtual MinorantBundle& get_bundle() 
  { return bundle.back();}  

  ///applies the AffineFunctionTransformation to constant_minorant and bundle, where (if given) only the global_indices of the transformed subgradients are required which need the local_indices only. If precomputed is given, it may contain some or contains afterwards a map from original minorant to transformed minorant; retunrs 0 on success
  virtual int push_aft(const AffineFunctionTransformation* inaft,
			const CH_Matrix_Classes::Indexmatrix* global_indices,
			const CH_Matrix_Classes::Indexmatrix* local_indices,
			std::map<MinorantPointer,MinorantPointer>* precomputed=0);

  /// undo the last push_aft
  virtual int pop_aft();

  /// returns bundle.size() of the model (may differ from dim_model in case of fixing some subspace to zero)
  virtual CH_Matrix_Classes::Integer dim_bundle()
  {return CH_Matrix_Classes::Integer(get_bundle().size());}

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
						   
  /// C=beta*C+alpha*A*B where B and B may be transposed; carry out the model part of this beginning at startindex_model 
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
  
  /// set the local modelx value in modelx beginning with startindex (initialize it, do not add)
  virtual int get_modelx(CH_Matrix_Classes::Matrix& modelx,
			 CH_Matrix_Classes::Integer startindex_model)=0;

  /// set the local modeldx value in modeldx beginning with startindex (initialize it, do not add)
  virtual int get_modeldx(CH_Matrix_Classes::Matrix& modeldx,
			  CH_Matrix_Classes::Integer startindex_model)=0;

  /// set the local modeldx value in modeldx beginning with startindex (initialize it, do not add)
  virtual int get_modeldcstr(CH_Matrix_Classes::Matrix& modeldcstr,
			     CH_Matrix_Classes::Integer startindex_constraints)=0;

  /// adds opB transposed times modelx (without constant affine term) to the arguments
  virtual int add_modelx_aggregate(CH_Matrix_Classes::Real& val,
				   CH_Matrix_Classes::Matrix& vec,
				   MinorantBundle& global_bundle,
				   CH_Matrix_Classes::Integer startindex_bundle)=0;

  /// set the model violation for the current system solution 
  virtual int get_sysviol_model(CH_Matrix_Classes::Matrix& modelvec,
				CH_Matrix_Classes::Integer startindex_model,
				const CH_Matrix_Classes::Matrix& dy,
				MinorantBundle& global_bundle,
				CH_Matrix_Classes::Integer startindex_bundle)=0;

  /// set the constraint violation for the current system solution starting at this index
  virtual int get_sysviol_constraints(CH_Matrix_Classes::Matrix& constrvec,
				      CH_Matrix_Classes::Integer startindex_constr)=0;

  /// reset the starting point for the current y 
  virtual int reset_starting_point(const CH_Matrix_Classes::Matrix& y,
				   CH_Matrix_Classes::Real mu,
				   MinorantBundle& global_bundle,
				   CH_Matrix_Classes::Integer startindex_bundle)=0;

  /// compute the step in the model space given the step in the design space
  virtual int compute_step(const CH_Matrix_Classes::Matrix& ystep,
			   MinorantBundle& global_bundle,
			   CH_Matrix_Classes::Integer startindex_bundle)=0;

    /// store this computed step locally and compute the missing local dual step information
  virtual int computed_step(const CH_Matrix_Classes::Matrix& modelxstep,
			    CH_Matrix_Classes::Integer startindex_model,
			    const CH_Matrix_Classes::Matrix& modelconstrstep,
			    CH_Matrix_Classes::Integer startindex_constr)=0;

  /// move in the last computed step direction by a step of length alpha and compute and store the violation in this point for later use in 
  virtual int do_step(CH_Matrix_Classes::Real alpha,
		      const CH_Matrix_Classes::Matrix& y,
		      MinorantBundle& global_bundle,
		      CH_Matrix_Classes::Integer startindex_bundle)=0;

  /// If mu is not zero, always add the centering term for this mu as well, if append is false, add the Schur complement rhs for add_BtinvsysB, if append is true, fill in the rhs of the local system starting at startindex for the model and at startindex_constraints for the modelconstraints
  virtual int add_localrhs(CH_Matrix_Classes::Matrix& globalrhs, 
			   CH_Matrix_Classes::Real rhsmu,
			   CH_Matrix_Classes::Real rhscorr,
			   CH_Matrix_Classes::Integer startindex_model,
			   CH_Matrix_Classes::Integer startindex_constraints,
			   bool append,
			   MinorantBundle& bundle,
			   CH_Matrix_Classes::Integer startindex_bundle) =0;

  ///add the "scaled" minorant outer products to globalsys, where the correct minorants start at the given index
  virtual int add_BtinvsysB(CH_Matrix_Classes::Symmatrix& globalsys,
			    MinorantBundle& bundle,
			    CH_Matrix_Classes::Integer startindex_bundle) =0;

  /** @brief  given the Cholesky factorization LL' of *minus* the blocks A and B (contraints on design variables and Bundle-modelx) and LinvABrhs, solve for the local constraints C and add the new contribution of tracedual*LinvTrace to LinvABsol; store the tracedual in Crhs_and_sol but not yet locally (this will be done by computed_step() ). 
   */
  virtual int solve_constrsys(const CH_Matrix_Classes::Symmatrix& ABchol,
			      const CH_Matrix_Classes::Matrix& LinvABrhs,
			      CH_Matrix_Classes::Matrix& LinvABsol,
			      CH_Matrix_Classes::Integer startindex_model,
			      CH_Matrix_Classes::Matrix& Crhs_and_sol,
			      CH_Matrix_Classes::Integer startindex_constraints) =0;


  /** @brief add the diagonal of the Schur complemented blocks belonging to bundle and local constraints (used for diagonal preconditioning)
  */ 
  virtual int add_BCSchur_diagonal(CH_Matrix_Classes::Matrix& diagonal,
				   MinorantBundle& globalbundle,
				   CH_Matrix_Classes::Integer startindex_bundle)=0;


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
					 CH_Matrix_Classes::Integer startindex_bundle)=0;


  /** @brief compute the preconditioning low-rank representation of the Schur complementd blocks belonging to bundle and local constraints by adding a Johnson-Lindenstrauss projection onto the given subspace to glob_lowrank   
  
  */ 
  virtual int prepare_BCSchur_JLprecond(CH_Matrix_Classes::Matrix& glob_lowrank,
					CH_Matrix_Classes::Matrix& subspace,
					bool append_globtransp_times_mat_to_subspace,
					MinorantBundle& globalbundle,
					CH_Matrix_Classes::Integer startindex_bundle)=0;

  /// add the contributions to glob_rhs of the Schur complemented model block, and return local_rhs of the non complemented constraint block in the rows/columns/diagonal block starting at startindex_constraints
  virtual int add_Schur_rhs(CH_Matrix_Classes::Matrix& glob_rhs,
			    CH_Matrix_Classes::Matrix* local_rhs,
			    CH_Matrix_Classes::Real rhsmu,
			    CH_Matrix_Classes::Real rhscorr,
			    CH_Matrix_Classes::Integer startindex_constraints,
			    MinorantBundle& globalbundle,
			    CH_Matrix_Classes::Integer startindex_bundle)=0;

  /// multiply in_vec with the local contribution to the global main block and add it to out_vec; the other local multiplications are carried out externally with the information provide in prepare_Schur_precond and are not done here.
  virtual int add_Schur_mult(const CH_Matrix_Classes::Matrix& in_vec,
			     CH_Matrix_Classes::Matrix& out_vec,
			     const CH_Matrix_Classes::Matrix* in_cvec,
			     CH_Matrix_Classes::Matrix* out_cvec,
			     CH_Matrix_Classes::Integer startindex_constraints,
			     MinorantBundle& globalbundle,
			     CH_Matrix_Classes::Integer startindex_bundle)=0;

  /// informs the model about the step computed
  virtual int computed_Schur_step(const CH_Matrix_Classes::Matrix& xstep,
				  const CH_Matrix_Classes::Matrix& local_step,
				  CH_Matrix_Classes::Integer startindex_model,
				  MinorantBundle& globalbundle,
				  CH_Matrix_Classes::Integer startindex_bundle)=0;
  
  
  /// returns the value of constant offset plus global linear cost term for the current globalx
  virtual CH_Matrix_Classes::Real globalx_cost(const CH_Matrix_Classes::Matrix& globalx);

  /// adds opB transposed times modelx and constant affine term to the arguments
  virtual int add_modelx_aggregate(CH_Matrix_Classes::Real& val,
				   CH_Matrix_Classes::Matrix& vec);

  /// adds opB transposed times modelx and constant affine term to the arguments
  virtual int add_Bt_modelx(CH_Matrix_Classes::Real& val,
			    CH_Matrix_Classes::Matrix& vec)
  { return add_modelx_aggregate(val,vec); }
  
  /// computes and returns C=alpha*B*A+beta*C where B and A may be transposed; C needs to have the correct size on input but will be initialized if beta==0.
  virtual CH_Matrix_Classes::Matrix& B_times(const CH_Matrix_Classes::Matrix& A,
					  CH_Matrix_Classes::Matrix& C,
					  CH_Matrix_Classes::Real alpha=1.,
					  CH_Matrix_Classes::Real beta=0.,
					  int Btrans=0,
					     int Atrans=0);

  /// computes and returns C=alpha*A*(B)+beta*C where A and B may be transposed; C needs to have the correct size on input but will be initialized if beta==0.
   virtual CH_Matrix_Classes::Matrix& times_B(const CH_Matrix_Classes::Matrix& A,
					      CH_Matrix_Classes::Matrix& C,
					      CH_Matrix_Classes::Real alpha=1.,
					      CH_Matrix_Classes::Real beta=0.,
					      int Atrans=0,
					      int Btrans=0);   

  /// add B*Diag(diagvec)*Bt to S  in the principal block starting at startindex
  virtual CH_Matrix_Classes::Symmatrix& add_BDBt(const CH_Matrix_Classes::Matrix& diagvec,
						 CH_Matrix_Classes::Symmatrix& S,
						 bool minus=false,
						 CH_Matrix_Classes::Integer startindex=0);


  /// store the coupling matrix Bt (first block row in the system) starting at column start_col (enlarge Bt first if start_col==Bt.coldim(), but the row dimension must be correct)
  virtual CH_Matrix_Classes::Matrix& get_Bt(CH_Matrix_Classes::Matrix& Bt,
					    CH_Matrix_Classes::Integer start_col=0);
  
  /// get the vector formed by all model x variables
  virtual CH_Matrix_Classes::Matrix& get_x();
 
  /// get the vector formed by all delta model x variables
  virtual CH_Matrix_Classes::Matrix& get_dx();

  /// get the vector formed by all delta model x variables
  virtual CH_Matrix_Classes::Matrix& get_dcstr();

  /// get the model violation for the current system solution
  virtual CH_Matrix_Classes::Matrix& get_sysviol_model(const CH_Matrix_Classes::Matrix& dy);
  
  /// get the constraint violation for the current system solution
  virtual CH_Matrix_Classes::Matrix& get_sysviol_constraints();  

  /// for test outputs
  virtual void display_model_values(const CH_Matrix_Classes::Matrix& y,
				    MinorantBundle& global_bundle,
				    CH_Matrix_Classes::Integer startindex_bundle,
				    std::ostream& out)=0;

  /// for test outputs
  virtual void display_model_values(const CH_Matrix_Classes::Matrix& y,
				    std::ostream& out);

  /// initialize the model variables to a strictly feasible "central" starting point; this is the first call when the next QP problem is solved, so other initialization steps may be appropriate as well here. 
  int reset_starting_point(const CH_Matrix_Classes::Matrix& y,CH_Matrix_Classes::Real mu);
   
  /// compute the step in the model space given the step in the design space
  int compute_step(const CH_Matrix_Classes::Matrix& ystep);

  /// store the computed step and compute the missing dual step information
  int computed_step(const CH_Matrix_Classes::Matrix& modelxstep,
		    const CH_Matrix_Classes::Matrix& modelconstrstep);

  /// move in the last computed step direction by a step of length alpha
  int do_step(CH_Matrix_Classes::Real alpha,const CH_Matrix_Classes::Matrix& nexty);

    /// if mu is not zero, always add the centering term for this mu as well, if append is false, add the Schur complement rhs for add_BtinvsysB, if append is true, append the rhs of the local system
  int add_localrhs(CH_Matrix_Classes::Matrix& globalrhs,
		   CH_Matrix_Classes::Real rhsmu,
		   CH_Matrix_Classes::Real rhscorr,
		   CH_Matrix_Classes::Integer startindex_model,
		   CH_Matrix_Classes::Integer startindex_constraints,
		   bool append);

  ///add the "scaled" minorant outer products to globalsys
  int add_BtinvsysB(CH_Matrix_Classes::Symmatrix& globalsys);

  /** @brief  given the Cholesky factorization LL' of *minus* the blocks A and B (contraints on design variables and Bundle-modelx), solve for the local constraints C and then for AB and return the solution to the ABC-variables without storing them yet (this will be done by computed_step() ). 
   */
  virtual int solve_constrsys(const CH_Matrix_Classes::Symmatrix& ABchol,
			      CH_Matrix_Classes::Matrix& ABCrhs_and_sol,
			      CH_Matrix_Classes::Integer startindex_model,
			      CH_Matrix_Classes::Integer startindex_constraints);


  /** @brief add the diagonal of the Schur complemented blocks belonging to bundle and local constraints (used for diagonal preconditioning)
  */ 
  virtual int add_BCSchur_diagonal(CH_Matrix_Classes::Matrix& diagonal);


  /** @brief append to lowrank "large" columns that should serve well for generating a low rank projection of the Schur complemented model part. For each column i the coordinate sigma_guess(i) gives the Diag_inv-norm for this column. The parameter minval asks to ignore columns whose norms are smaller than minval. If diaginvval is positive, the vector Diag_inv is this value times the all ones vector. 

    On input lowrank must have the correct number of rows already but may
    have 0 columns.  
  */ 
  virtual int propose_BCSchur_pcsubspace(CH_Matrix_Classes::Matrix& lowrank,
					 CH_Matrix_Classes::Matrix& sigma_guess,
					 const CH_Matrix_Classes::Matrix& Diag_inv,
					 CH_Matrix_Classes::Real minval,
					 CH_Matrix_Classes::Real diaginvval=-1.);

  /** @brief compute the preconditioning low-rank representation of the Schur complementd blocks belonging to bundle and local constraints by adding a Johnson-Lindenstrauss projection onto the given subspace to glob_lowrank   
  
  */ 
  virtual int prepare_BCSchur_JLprecond(CH_Matrix_Classes::Matrix& glob_lowrank,
					CH_Matrix_Classes::Matrix& subspace,
					bool append_globtransp_times_mat_to_subspace=false);

  /// add the contributions to glob_diagonal and glob_rhs of the Schur complemented parts, and return local_rhs, local_globblock, local_diagblock of the non complemented parts 
  virtual int add_Schur_rhs(CH_Matrix_Classes::Matrix& glob_rhs,
			    CH_Matrix_Classes::Matrix* local_rhs,
			    CH_Matrix_Classes::Real rhsmu,
			    CH_Matrix_Classes::Real rhscorr);

  /// multiply in_vec with the local contribution to the global main block and add it to out_vec; the other local multiplications are carried out externally with the information provide in prepare_Schur_precond and are not done here.
  virtual int add_Schur_mult(const CH_Matrix_Classes::Matrix& in_Qvec,
			     CH_Matrix_Classes::Matrix& out_Qvec,
			     const CH_Matrix_Classes::Matrix* in_Cvec=0,
			     CH_Matrix_Classes::Matrix* out_Cvec=0);

  /// use the computed step information to also compute the steps of the complemented parts
  virtual int computed_Schur_step(const CH_Matrix_Classes::Matrix& xstep,
				  const CH_Matrix_Classes::Matrix& local_step);
   

};

 /** @brief Implementation of the QPModelDataPointer Interface for BundelModel for generating the correct type of blocks for QPSolver and for setting the final block in the solver 

     The objects generated here are the implementations
     QPSumModelBlock and QPConeModelBlock of the base class QPModelBlock. 
     The ownership of the generated objects is passed over to the
     calling routine. Also the pointer to the initial block stored
     here only gives access to the block, but the object is not
     owned and may not be deleted here. 
*/

  
class QPModelPointer: public virtual QPModelDataPointer
{
protected:
  QPModelBlock* model_block;  ///< stores a pointer to the current starting block giving access to the cutting model(s) [it does not own or delete this object]
public:
  /// default constructor
  QPModelPointer(CBout* cb=0,int cbinc=-1):QPModelDataPointer(cb,cbinc),model_block(0){}
  /// virtual destructor
  virtual ~QPModelPointer();

  /// set the pointer to NULL
  void clear_model_data_ptr()
  {model_block=0;}

  
  ///store the pointer to the object if it matches the required type for the QP solver, otherwise return a nonzero value as error; this is used in the models to return the local qp model data
  int set_model_data(QPModelDataObject* inbp)
  { model_block=dynamic_cast<QPModelBlock*>(inbp); return (model_block==0);}

  ///returns a new QPSumModelDataObject, that has to be deleted by the caller. The argument is optional and allows to potentially generate different blocks for different derived BundleModel objects; this is used in SumModel to collect the models of the various oracles that are summed over 
  QPSumModelDataObject* generate_summodel_data(BundleModel* bmp=0);
 
  ///returns a new QPConeModelDataObject suitable for the default conic BundleModel implementations; it has to be deleted by the caller. The argument is optional and allows to potentially generate specialized objects for special BundleModel objects 
  QPConeModelDataObject* generate_conemodel_data(BundleModel* bmp=0);
  
  /// returns the pointer value
  QPModelDataObject* get_model_data_ptr() const {return model_block;}
};

  
  //@}

}

#endif

