/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/QPSolver.hxx
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


#ifndef CONICBUNDLE_QPSOLVER_HXX
#define CONICBUNDLE_QPSOLVER_HXX

/**  @file QPSolver.hxx
    @brief Header declaring the class ConicBundle::QPSolver
    @version 1.0
    @date 2020-03-18
    @author Christoph Helmberg
*/

#include "QPSolverBasicStructures.hxx"
#include "Groundset.hxx"
#include "LPGroundsetModification.hxx"

namespace ConicBundle {


/** @defgroup ConstrainedQPSolver Internal QP Solver for linearly constrained groundsets

    If the groundset (see e.g. LPGroundset) of the convex
    optimization problem involves box constraints with a general
    quadratic proximal term (i.e. not just a diagonal) or general
    linear constraints even with diagonal proximal term, a constrained
    QP solver is required to determine the next candidate as the
    solution of the bundle subproblem. Because ConicBundle is designed
    to handle also nonpolyhedral cutting models arising from support
    functions over the second order cone and positive semidefinite
    cone the general solver is an intererior point method. More
    precisely, it implements a primal dual predictor corrector method.

    The main visible object is QPSolver. This mainly serves as an
    interface for ConicBundle as it implements a QPSolverObject (for
    finding the next candidate via solving the quadratic bundle
    subproblem) and a QPDataPointer (for collection the cutting model
    information). The internal solver routines are mostly defined in
    QPSolverBasicStructures. The basic structures serve to describe
    all parts in connection with the linearly constraint groundset and
    the general iterative procedure of when to collect which data,
    when to solver the primal dual KKT-system and when to perform line
    searches. The central part of setting up and solving the primal
    dual KKT-system is supplied externally via the QPSolverParameters,
    which must specify a QPKKTSolverObject.  Currently the following
    options exist for such an QPKKTSolverObject:

    - QPDirectKKTSolver : by a Schur complement the quadratic term 
      is computed into the lines of the linear constraints and the
      constraints describing the models. This system is then solved
      either with an aditional Schur complement by Cholesky or
      directly in indefinite form by Aasen's method. 
      
    - QPIterativeKKTSolver : It mainly offers a general interface for
      iterative methods, but also provides the matrix times vector
      multiplication for the entire KKT-System.  This iterative solver
      class always just provides the matrix vector product for a
      special variant of the KKT-System in the style of a
      CH_Matrix_Classes::IterativeSystemObject, the actual iterative
      solver with a suitable preconditioner must always be specified
      during construction by passing a
      CH_Matrix_Classes::IterativeSolverObject and possibly a
      QPKKTPrecondObject that fits the iterative solver and the system
      variant. So far no specific preconditioner is implemented for
      this system configuration, but CH_Matrix_Classes::Psqmr and
      CH_Matrix_Classes::MinRes are possible iterative solver options. 
     
    - QPIterativeKKTHASolver : is dervied from QPIterativeKKTSolver and
      differs in that in the KKT-system the model block is Schur
      complemented into the lines of the prox-termin (the H block),
      but the linear constraints of the ground set (the A block) is
      left untouched.  So the KKT-System is reduced to the H and A
      blocks and the class implements the matrix vector product for
      this system variant. A preconditioning class for this setting
      is provided in QPKKTSubspaceHPrecond, which provides various
      options for obtaining a preconditioner for the H-block (and
      so far only for the H-block) via subspace projection and
      low rank inversion. If the linear constraint block A is friendly,
      this works fine with the iterative solvers CH_Matrix_Classes::MinRes
      CH_Matrix_Classes::Psqmr. If there is no A at all (but maybe 
      box constraints), the best method is likely CH_Matrix_Classes::PCG.
      
    - QPIterativeKKTHAeqSolver : is dervied from QPIterativeKKTSolver and
      differs to QPIterativeKKTHASolver in that in the KKT-system also
      the inequality constraints of the A block are Schur complemented
      into the H block. Thus, in the presence of equality constraints
      the iterative solvers CH_Matrix_Classes::MinRes is needed,
      if there are inequality constraints in A only then 
      CH_Matrix_Classes::PCG can also be used. 
      The combination with MinRes seems to be a bit more stable and
      together with preconditioner 31 of  QPKKTSubspaceHPrecond
      this is currently the recommended combination for the iterative case. 
      
    The cutting models are realized by two implemtented classes
    derived from QPModelBlock which brings together the abstract
    classes QPModelDataObject (interface for the cutting modles of
    ConicBundle) and QPModelBlockObject (interface for
    QPSolverBasicInterface and the QPKKTSolverObject)

    - QPConeModelBlock : it allows to describe any of the cutting models
      employed by the classes NNCModel, SOCModel, PSCModel, BoxModel
      by tying together the following implementations of 
      the abstract class InteriorPointBundleBlock
      * NNCIPBundleBlock (building on NNCIPBlock for a nonnegative cone)
      * SOCIPBundleBlock (building on SOCIPBlock for a second order cone)
      * PSCIPBundleBlock (building on PSCIPBlock for a positive semidefnite cone)
      * BoxIPBundleBlock (building also on NNCIPBlock for the bounds) 
   
    - QPSumModelBlock : it forms the sum of several QPModelBlock instances
      as required in SumModel

 
*/
//@{


class GroundsetModification;

    /** @brief QPSolver is the access point for ConicBundle to the internal constrained QP Solver, see \ref ConstrainedQPSolver
   */

  class QPSolver:public virtual QPSolverBasicStructures
{
private:
  /** @brief this class facilitates switching between original data and preprocessed data that removes fixed variables 
   */
  class QPProblemData {
  public:
    /** @name data describing the feasible set (usually this does not change over several calls
     */
    //@{
    CH_Matrix_Classes::Integer dim;  ///< dimension of the y ground set
    CH_Matrix_Classes::Matrix lby;    ///< lower bounds on y, same dim as y
    CH_Matrix_Classes::Matrix uby;    ///< upper bounds on y, same dim as y
    CH_Matrix_Classes::Indexmatrix lbindex;  ///< indices i with y(i) lower bounded (sorted increasingly)
    CH_Matrix_Classes::Indexmatrix ubindex;  ///< indices i with y(i) upper bounded (sorted increasingly)
    CH_Matrix_Classes::Sparsemat A;  ///< feasible set is lby<=y<=uby with rhslb<=Ay<=rhsub
    CH_Matrix_Classes::Matrix rhslb; ///< right hand side lower bounds 
    CH_Matrix_Classes::Matrix rhsub; ///< right hand side upper bounds
    CH_Matrix_Classes::Indexmatrix rhslbindex;  ///< indices i with (Gy)(i) lower bounded (sorted increasingly; used for slacks, equations are not listed)
    CH_Matrix_Classes::Indexmatrix rhsubindex;  ///< indices i with (Gy)(i) upper bounded (sorted increasingly, used for slacks, equations are not listed)
    CH_Matrix_Classes::Indexmatrix rhseqindex;
    //@}
    
    
    /** @name data describing the quadratic cost function,  contributors are proximal term, groundset minorant and cutting model
     */
  //@{
  //data for describing the quadratic cost matrix Q
    BundleProxObject* Hp;  ///< points to the quadratic cost matrix inside ensure_feasibility() and candidate(), otherwise ==NULL
    CH_Matrix_Classes::Matrix c; ///< linear cost term for QP-subproblems
    CH_Matrix_Classes::Real gamma; ///< constant offset for QP-subproblems
    
    CH_Matrix_Classes::Matrix groundset_c; ///< linear cost term of the groundset
    CH_Matrix_Classes::Real groundset_gamma; ///< constant offset of the groundset
    //@}
  };
  
  /** @name data describing the ground set for the internal qp solver
   */
  //@{
  /// collects the original problem data and the data describing the quadratic costs
  QPProblemData original_data;

  /// collects the preprocessed problem data (some fixed variables may be eliminated) and the data describing the quadratic costs
  QPProblemData preproc_data;
  
  CH_Matrix_Classes::Indexmatrix preproc_fixed; ///< 1 if coordinate has fixed value, 0 otherwise; dimension is (0,0) if not initialized
  CH_Matrix_Classes::Matrix preproc_fixedval; ///< if fixed it gives the value the coordinate is fixed to, otherwise the value is zero 
  CH_Matrix_Classes::Indexmatrix preproc_indices; ///< the sequence of original non fixed indices used in preproc_data; dimension is (0,0) if not initialized
  AffineFunctionTransformation preproc_aft; ///< the transformation of the bundle
  std::map<MinorantPointer,MinorantPointer> preproc_bundle_projection; ///< stores the transformed minorants of the bundle
  
  QPProblemData* qp_data; ///< either points to original_data or to fixing_data
  //@}

  /** @name the solution data of the last solve in the original data space
   */
  //@{
  CH_Matrix_Classes::Real sol_val_lb;
  CH_Matrix_Classes::Real sol_val_ub;
  CH_Matrix_Classes::Matrix sol_point;
  CH_Matrix_Classes::Real gs_aggr_offset;
  CH_Matrix_Classes::Matrix gs_aggr_gradient;
  //@}
  
  
  
  /** @name termination parameters for the QP
   */
  //@{
  CH_Matrix_Classes::Real primal_infeasibility_eps; ///< primal infeasibility termination parameter 
  CH_Matrix_Classes::Real dual_infeasibility_eps; ///< dual infeasibility termination parameter
  CH_Matrix_Classes::Real lower_bound_gap_eps; ///< distance to lower bound termination parameter
  CH_Matrix_Classes::Real upper_bound_gap_eps; ///< distance to upper bound termination parameter
  CH_Matrix_Classes::Real objective_gap_eps; ///< relative precision for QP objective termination
  //@}

  /** @name internal support routines
   */
  //@{
  
  /// given lby,uby,rhslb,rhsub compute lbindex,ubindex,rhslbindex,rhsubindex
  int determine_indices(QPProblemData& qpd);

  /// initialize preproc_data for 
  int preprocess_data(const CH_Matrix_Classes::Matrix& center_y,
		      CH_Matrix_Classes::Indexmatrix* yfixed,
		      bool& no_changes);

  int postprocess_data(bool round_to_active_bounds);

  //@}
  
protected:

  /** @name implementations of QPBasicInterface routines not needed externally
   */
  //@{
  /// QP-solver-interface routine, returns primal dimension (length of y)
  CH_Matrix_Classes::Integer QPget_xdim() const {return qp_data->dim;}
  /// QP-solver-interface routine, returns primal dimension (length of y)
  CH_Matrix_Classes::Integer QPget_ydim() const {return qp_data->A.rowdim();}
 
  /// QP-solver-interface routine for adding quadratic matrix times vector
  CH_Matrix_Classes::Matrix&  QPadd_Qx(const CH_Matrix_Classes::Matrix& xin,CH_Matrix_Classes::Matrix& outplusQx) const 
  {return qp_data->Hp->add_Hx(xin,outplusQx);}
 
  /// QP-sover-interface routine for linear cost term
  const CH_Matrix_Classes::Matrix& QPget_c() const {return qp_data->c;}

  /// QP-sover-interface routine for constant cost term
  CH_Matrix_Classes::Real QPget_gamma() const {return qp_data->gamma;}

  /// QP-sover-interface routine for constraint lower bounds
  const CH_Matrix_Classes::Matrix& QPget_rhslb() const {return qp_data->rhslb;}

  /// QP-sover-interface routine for constraint upper bounds
  const CH_Matrix_Classes::Matrix& QPget_rhsub() const {return qp_data->rhsub;}

  /// QP-sover-interface routine for indices of constraint lower bound slacks 
  const CH_Matrix_Classes::Indexmatrix& QPget_rhslbind() const {return qp_data->rhslbindex;}

  /// QP-sover-interface routine for indices of constraint upper bound slacks 
  const CH_Matrix_Classes::Indexmatrix& QPget_rhsubind() const {return qp_data->rhsubindex;}

  /// QP-solver-interface routine for variable lower bounds
  const CH_Matrix_Classes::Matrix& QPget_lb() const {return qp_data->lby;}

  /// QP-solver-interface routine for variable upper bounds
  const CH_Matrix_Classes::Matrix& QPget_ub() const {return qp_data->uby;}

  /// QP-solver-interface routine for indices of variable lower bounds > ConicBundle::CB_minus_infinity
  const CH_Matrix_Classes::Indexmatrix& QPget_lbind() const {return qp_data->lbindex;}

  /// QP-solver-interface routine for indices of variable upper < ConicBundle::CB_plus_infinity
  const CH_Matrix_Classes::Indexmatrix& QPget_ubind() const {return qp_data->ubindex;}

  /// QP-solver-interface routine for adding constraint matrix times vector
  CH_Matrix_Classes::Matrix& QPadd_Ax(const CH_Matrix_Classes::Matrix& xin,CH_Matrix_Classes::Matrix& outplusAx) const 
  {return CH_Matrix_Classes::genmult(qp_data->A,xin,outplusAx,1.,1.);}
  
  /// QP-solver-interface routine for adding transposed constraint matrix times vector
  CH_Matrix_Classes::Matrix& QPadd_Aty(const CH_Matrix_Classes::Matrix& yin,CH_Matrix_Classes::Matrix& outplusAty) const 
  {return CH_Matrix_Classes::genmult(qp_data->A,yin,outplusAty,1.,1.,1);}
  //@}
  

public:
  
  /// (re)initialize to empty
  void QPclear();

  /// default constructor
  QPSolver(CBout* cb=0,int cbinc=-1):
    CBout(cb,cbinc),QPSolverBasicStructures(0,cb)
  { original_data.Hp=0;preproc_data.Hp=0; QPclear();}
  	   
  ~QPSolver(){}

  //------------  QPSolverOject routines

  /// check whether the parameters are QPSolverParameters and set them if so
  int QPset_parameters(QPSolverParametersObject* params)
  {
    QPSolverParameters* qpp=dynamic_cast<QPSolverParameters*>(params);
    if (qpp==0)
      return 1;
    return QPset_solver_parameters(qpp);
  }

  /// yfixing is currently not supported, this returns false.
  bool QPsupports_yfixing()
  { return false;}

  /** @brief return true iff the code supports QPupdate(), i.e., it supports external updates of the groundset aggregate in order to model constraints not included explicitly in the QP's model   
  */
  virtual bool QPsupports_updates()
  { return false;}
  
  /// returns the current lower bound on the optimal value (if feasibility is good enough)
  CH_Matrix_Classes::Real QPget_lower_bound()
  { return QPget_dualval(); }

  /// solve the current bundle subproblem so that precision requirements are met (see \ref InternalQPSolverInterface)
  int QPsolve(const CH_Matrix_Classes::Matrix& center_y,
	      CH_Matrix_Classes::Real lower_bound,
	      CH_Matrix_Classes::Real upper_bound,
	      CH_Matrix_Classes::Real relprec,
	      QPSolverProxObject* Hp,
	      const MinorantPointer& gs_aggr,
	      CH_Matrix_Classes::Indexmatrix* yfixed);

  /// solve the bundle subproblem for updated box multipliers so that precision requirements are met (see \ref InternalQPSolverInterface). This routine is typically not called for this solver, because box constraints are included explicitly. 
  int QPupdate(const CH_Matrix_Classes::Matrix& center_y,
	       CH_Matrix_Classes::Real lower_bound,
	       CH_Matrix_Classes::Real upper_bound,
	       CH_Matrix_Classes::Real relprec,
	       QPSolverProxObject* Hp,
	       const MinorantPointer& gs_aggr,
	       CH_Matrix_Classes::Indexmatrix* yfixed,
	       const MinorantPointer& delta_gs_aggr,
	       const CH_Matrix_Classes::Indexmatrix& delta_index);

  /// resolve the bundle subproblem (usually because of modified penalty parameters) so that precision requirements are met (see \ref InternalQPSolverInterface)
  int QPresolve(CH_Matrix_Classes::Real lower_bound,
		CH_Matrix_Classes::Real upper_bound,
		CH_Matrix_Classes::Real relprec);

  /// retrieve the solution produced (see \ref InternalQPSolverInterface)
  int QPget_solution(CH_Matrix_Classes::Real& augval_lb,
		     CH_Matrix_Classes::Real& augval_ub,
		     CH_Matrix_Classes::Matrix& new_point,
		     CH_Matrix_Classes::Real& gsaggr_offset,
		     CH_Matrix_Classes::Matrix& gsaggr_gradient);

  /// currently it does nothing
  std::ostream& QPprint_statistics(std::ostream& out,int /* printlevel*/ =0)
  {return out;}
  
  /// return a new modification object on the heap that is initialized for modification of *this  
  GroundsetModification* QPstart_modification()
  {return new LPGroundsetModification(original_data.dim,original_data.A.rowdim(),this);}

  /// groundset changes are communicated to the solver here
  int QPapply_modification(const GroundsetModification& mdf);

  /// check feasiblity of y for the current groundset constraints
  bool QPis_feasible(const CH_Matrix_Classes::Matrix& y,
		   CH_Matrix_Classes::Real relprec=1e-10);

  /// makes y feasible if not so, see Groundset::ensure_feasibility()
  int QPensure_feasibility(CH_Matrix_Classes::Matrix& y,
			   bool& ychanged,
			   QPSolverProxObject* inHp,
			   CH_Matrix_Classes::Real relprec=1e-10);

  /// solve the quadratic problem for the given cost function and precision (without cutting model, usually for finding feasible starting points)
  int solve(BundleProxObject* Hp,
	    const CH_Matrix_Classes::Matrix& c,
	    CH_Matrix_Classes::Real gamma,
	    CH_Matrix_Classes::Real lowerbound,
	    CH_Matrix_Classes::Real upperbound,
	    CH_Matrix_Classes::Real relprec,
	    CH_Matrix_Classes::Real skip_factor);
	   
  /// returns true if, for the current constraints and the requested ProxObject, it might be better to use the internal unconstrained QP solver (which can deal with box constraints by a work-around)
  bool QPprefer_UQPSolver(QPSolverProxObject*) const;

  /// returns false if the feasible set is the entire space (unconstrained optimization), true otherwise. 
  bool QPconstrained() const
  {return (original_data.lbindex.rowdim()+original_data.ubindex.rowdim()+original_data.A.rowdim()==0)?false:true;}

  /// number of linear constraints
  CH_Matrix_Classes::Integer rowdim() const {return original_data.A.rowdim();}

  /// returns the lower bounds on y
  const CH_Matrix_Classes::Matrix& get_lby() const { return original_data.lby;}

  /// returns the upper bounds on y
  const CH_Matrix_Classes::Matrix& get_uby() const { return original_data.uby;}

  /// returns the indices of variable lower bounds > ConicBundle::CB_minus_infinity
  const CH_Matrix_Classes::Indexmatrix& get_lbindex() const { return original_data.lbindex;}

  /// returns the indices of variable lower bounds < ConicBundle::CB_plus_infinity
  const CH_Matrix_Classes::Indexmatrix& get_ubindex() const { return original_data.ubindex;}

  /// returns the constraint matrix of the feasible set
  const CH_Matrix_Classes::Sparsemat& get_A() const { return original_data.A;}

  /// returns the constraint lower bounds
  const CH_Matrix_Classes::Matrix& get_rhslb() const {return original_data.rhslb;}

  /// returns the constraint upper bounds
  const CH_Matrix_Classes::Matrix& get_rhsub() const {return original_data.rhsub;}

  /// returns the indices with constraint lower bound slacks
  const CH_Matrix_Classes::Indexmatrix& get_rhslbind() const {return original_data.rhslbindex;}

  /// returns the indices with constraint upper bound slacks
  const CH_Matrix_Classes::Indexmatrix& get_rhsubind() const {return original_data.rhsubindex;}

  /// return true if only box constraints are present and return the present box constraints
  bool QPboxconstrained(const CH_Matrix_Classes::Matrix*& lb, ///< lower bounds
			const CH_Matrix_Classes::Matrix*& ub, ///< upper bounds
			const CH_Matrix_Classes::Indexmatrix*& lbind, ///< indices with "finite lower bounds"
			const CH_Matrix_Classes::Indexmatrix*& ubind ///< indices with "finite upper bounds" 
			) const
  {
    lb=&get_lby();
    ub=&get_uby();
    lbind=&get_lbindex();
    ubind=&get_ubindex();
    if ((!QPconstrained())||(original_data.A.rowdim()!=0))
      return false;
    return true;
  }

  /// output the data describing the QP in m-file style 
  int mfile_data(std::ostream& out) const;

  /// set output settings
  void set_cbout(const CBout* cb, int incr=-1)
  {
    CBout::set_cbout(cb,incr);
    if (model_block)
      model_block->set_cbout(cb,incr);
  }

};




  //@}

}

#endif

