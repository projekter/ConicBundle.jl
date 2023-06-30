/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/LPGroundset.hxx
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




#ifndef CONICBUNDLE_LPGROUNDSET_HXX
#define CONICBUNDLE_LPGROUNDSET_HXX


/**  @file LPGroundset.hxx
    @brief Header declaring the class ConicBundle::Groundset
    @version 1.0
    @date 2014-07-14
    @author Christoph Helmberg
*/

#include "Groundset.hxx"
#include "QPSolver.hxx"
#include "UQPSolver.hxx"

namespace ConicBundle {

/** @ingroup InternalGroundset
 * @{
*/

/** @brief represents a polyhedral feasible ground set by affine constraints as in Linear Programming  with special support for simple box constraints
 
The polyhedral feasible set is described in the form
   \f[Y=\{y\in\mathbf{R}^m\colon l\le y\le u, r_l\le Ay \le r_u\} \f]
   where \f$l\le u\f$ are \f$m\f$-vectors, \f$A\f$ is an \f$m\times n\f$
   matrix and \f$r_l\le r_u\f$ are \f$n\f$-vectors. \f$m\f$ and \f$n\f$
   may be zero (empty), lower and upper bounds may take values
   ConicBundle::CB_minus_infinity and ConicBundle::CB_plus_infinity to
   implement appropriate inequalties. 
   
For general polyhedral ground sets the two routines ensure_feasibility()
   and candidate() need to solve a general convex QP; this is currently
   done by interior point methods.

If there are only box constraints and a diagonal scaling matrix for the
   proximal term, then all computations are done directly witout the
   need to solve a full convex QP. In this case the routine also offers
   a scaling heuristic and a heuristic for fixing y variables at their
   bounds.

*/


class LPGroundset:public Groundset
{
private:
  /// current dimension of the ground set
  CH_Matrix_Classes::Integer dim;  

  /// nonnegative update counter for recognizing changes in the groundset description
  CH_Matrix_Classes::Integer groundset_id;

  /// the starting point
  CH_Matrix_Classes::Matrix starting_point;

  /// the ground set aggregate, here giving the linear coefficients of the affine cost function
  MinorantPointer gs_aggregate;

  /// if set to true in constrained versions, y-variables that are at strongly active bounds should be fixed to this value by indicating so in @a yfixed
  bool use_yfixing;                       

  /// in constrained versions, when @a use_yfixing is true, this indicates whether a coordinate should be considered as free (0), fixed (1), or newly fixed (2)
  CH_Matrix_Classes::Indexmatrix yfixed;

  //the QPSolver for solving the bundle subproblem with Groundset
  QPSolverObject* qpsolver;   ///< general solver for all groundsets
  UQPSolver uqpsolver;  ///< special solver for easy groundsets
  QPSolverObject* qp_solver; ///< the selected solver
  bool solve_model_without_gs; ///< selected varaint of solving the subproblem


  //data for describing the QP
  BundleProxObject* Hp;  ///< points to the quadratic cost matrix inside ensure_feasibility() and candidate(), otherwise ==NULL
  CH_Matrix_Classes::Matrix c; ///< linear cost term for QP-subproblems
  CH_Matrix_Classes::Real gamma; ///< constant offset for QP-subproblems

  MinorantPointer gs_minorant; ///<  affine cost function for the ground set

  //variable metric information

  VariableMetricSelection* vm_selection; ///< specifieable routine for computing a metric contribution of the groundset, see \ref InternalVariableMetric
  
  bool store_old_minorants; ///< false by default, set to true if variable metric is used
  MinorantBundle old_minorants; ///< minorants collected for variable metric purpose
  CH_Matrix_Classes::Integer max_minorants; ///< maximum number of minorants
  CH_Matrix_Classes::Integer minorant_nexti; ///< used by the variable metric low rank heuristic 
  CH_Matrix_Classes::Matrix old_lowrank; ///< used by the variable metric low rank heuristic 
  CH_Matrix_Classes::Matrix old_diagonal; ///< used by the variable metric low rank heuristic 
  CH_Matrix_Classes::Symmatrix old_sym; ///< used by the variable metric heuristic 
  
protected:

public:
  /// resets all values as described in Groundset::clear()
  virtual void clear( CH_Matrix_Classes::Integer indim=0, 
		      CH_Matrix_Classes::Integer in_groundset_id=0);

  /// calls clear() with the same parameters
  LPGroundset(CBout* cb=0):CBout(cb),vm_selection(0){qpsolver=new QPSolver(cb); uqpsolver.set_cbout(cb); clear(0,0);}

  /// allows to specify the groundset in the constructor, zero is allowed everywhere
  LPGroundset(CH_Matrix_Classes::Integer dim,
	      const CH_Matrix_Classes::Matrix* lbyp=0,
	      const CH_Matrix_Classes::Matrix* ubyp=0,
	      const CH_Matrix_Classes::Sparsemat* Gp=0,
	      const CH_Matrix_Classes::Matrix* rhslbp=0,
	      const CH_Matrix_Classes::Matrix* rhsubp=0,
	      const CH_Matrix_Classes::Matrix* start_val=0,
	      const CH_Matrix_Classes::Matrix* costs=0,
	      const CH_Matrix_Classes::Real offset=0., 
	      CH_Matrix_Classes::Integer in_groundset_id=0,
	      CBout* cb=0);
  
  ~LPGroundset(){delete qpsolver; delete vm_selection;}

  /// returns the current groundset_id, increased values indicate changes in the ground set  
  CH_Matrix_Classes::Integer get_groundset_id() const {return groundset_id;}

  /// sets the groundset_id to the desired value, increasing it is safer here because this is used to indicate changes
  void set_groundset_id(CH_Matrix_Classes::Integer gsid){groundset_id=gsid;}

  /// returns the dimension of the ground set, i.e., the length of the variables vector y
  virtual CH_Matrix_Classes::Integer get_dim() const 
  {return dim;}

  ///Set the qp solver's parameters to qpparams (if not null); if the second argument qpsolver is also given, the old solver is first discarded and replaced by this new solver and then the parameters are set (if given).  Any object passed here will be owned and deleted by *this. For correct continuaton a new qpsolver needs to have the same feasible set as the current solver but this must be ensured by the caller.  
  int set_qpsolver(QPSolverParametersObject* qpparams,QPSolverObject* qpsolver=0);

  /// returns false if the feasible set is the entire space (unconstrained optimization), true otherwise. 
  bool constrained() const
  {return qpsolver->QPconstrained();}

  /// number of linear constraints
  //CH_Matrix_Classes::Integer rowdim() const {return qpsolver.rowdim();}

  /// returns the lower bounds vector on y if it exists
  const CH_Matrix_Classes::Matrix* get_lby() const;

  /// returns the upper bounds vector on y if it exists
  const CH_Matrix_Classes::Matrix* get_uby() const;

  /// returns the indices of variable lower bounds > ConicBundle::CB_minus_infinity
  //const CH_Matrix_Classes::Indexmatrix& get_lbindex() const { return qpsolver.get_lbindex();}

  /// returns the indices of variable lower bounds < ConicBundle::CB_plus_infinity
  //const CH_Matrix_Classes::Indexmatrix& get_ubindex() const { return qpsolver.get_ubindex();}

  /// returns the constraint matrix of the feasible set
  //const CH_Matrix_Classes::Sparsemat& get_G() const { return qpsolver.get_A();}

  /// returns the constraint lower bounds
  //const CH_Matrix_Classes::Matrix& get_rhslb() const {return qpsolver.get_rhslb();}

  /// returns the constraint upper bounds
  //const CH_Matrix_Classes::Matrix& get_rhsub() const {return qpsolver.get_rhsub();}

  /// returns the indices with constraint lower bound slacks
  //const CH_Matrix_Classes::Indexmatrix& get_rhslbind() const {return qpsolver.get_rhslbind();}

  /// returns the indices with constraint upper bound slacks
  //const CH_Matrix_Classes::Indexmatrix& get_rhsubind() const {return qpsolver.get_rhsubind();}

  /// returns true if still feasible, see Groundset::is_feasible() 
  virtual bool is_feasible(CH_Matrix_Classes::Integer& in_groundset_id,
			   const CH_Matrix_Classes::Matrix& y,
			   CH_Matrix_Classes::Real relprec=1e-10);

  /// makes y feasible if not so, see Groundset::ensure_feasibility()
  int ensure_feasibility(CH_Matrix_Classes::Integer& in_groundset_id,
			 CH_Matrix_Classes::Matrix& y,
			 bool& ychanged,
			 BundleProxObject* Hp,
			 CH_Matrix_Classes::Real relprec=1e-10);

  /// returns a pointer to an internal QPSolverObject that is able to solve bundle suproblems efficiently for this kind of groundset and scaling; if solves_model_without_gs == true the qp solver does not include the groundset and the groundset has to be dealt with by the Gauss Seidel approach 
  virtual QPSolverObject* get_qp_solver(bool& solves_model_without_gs,
					BundleProxObject* Hp); 

  /// set parameters for the QP_Solver  
  int set_qp_solver_parameters(QPSolverParametersObject* /* in_qpsp */)
  {
    return 0;
  }

  /// returns a stored starting point, note: this need not be feasible; if generated automatically, its dimension is correct.  
  virtual const CH_Matrix_Classes::Matrix& get_starting_point() const 
  {return starting_point;}

  /// stores the a new starting point irrespective of whether it is feasible or not and returns 0 if it feasible, 1 if it is infeasible
  virtual int set_starting_point(const CH_Matrix_Classes::Matrix& vec)
  {starting_point=vec; int dummy=-1; return (is_feasible(dummy,vec)==false);}

  /// computes the next ground set minorant and candidate, see Groundset::candidate()
  virtual int candidate(CH_Matrix_Classes::Integer& gs_id,
                        CH_Matrix_Classes::Matrix& newy,
			CH_Matrix_Classes::Real& cand_gs_val,
			CH_Matrix_Classes::Real& linval,
			CH_Matrix_Classes::Real& augval_lb,
			CH_Matrix_Classes::Real& augval_ub,
			CH_Matrix_Classes::Real& subgnorm2,
			const CH_Matrix_Classes::Matrix& center_y,
			CH_Matrix_Classes::Real center_value,
			const MinorantPointer& model_minorant,
			BundleProxObject* Hp,
			MinorantPointer* delta_groundset_minorant=0,
			CH_Matrix_Classes::Indexmatrix* delta_index=0,
			CH_Matrix_Classes::Real relprec=1e-2);
  
  /// returns the groundset aggregate computed in candidate()
  virtual const MinorantPointer& get_gs_aggregate() const { return gs_aggregate;}

  /// returns the linear minorant valid on the entire ground set (e.g. a linear cost funciton)
  const MinorantPointer& get_gs_minorant() const { return gs_minorant;}

  /// if not NULL (iff get_use_yfixing()==false) it returns the vector yfixed with yfixed(i)=0 if not fixed, =1 is fixed already, =2 if newly fixed
  virtual const CH_Matrix_Classes::Indexmatrix* get_yfixed() const
  {return (use_yfixing? &yfixed : 0);}

  /// if not NULL (iff get_use_yfixing()==false) returns the vector yfixed with yfixed(i)=0 if not fixed, =1 is fixed already, =2 if newly fixed
  virtual CH_Matrix_Classes::Indexmatrix* set_yfixed()
  {return (use_yfixing? &yfixed : 0);}

  /// true if the cooridinate fixing heuristic is switched on (only constrained cases)
  bool get_use_yfixing() const {return use_yfixing;}

  /// set to true to switch on the cooridinate fixing heuristic (only constrained cases)
  void set_use_yfixing(bool uyf) {use_yfixing=uyf;}

  /// delete old selector and set a new one (0 is allowed resulting in no local selector)
  int set_variable_metric_selection(VariableMetricSelection* vms=0)
  { delete vm_selection; vm_selection=vms; return 0;}

  /// delete old selector and set a new one (0 is allowed resulting in no local selector)
  VariableMetricSelection* get_variable_metric_selection() const
  { return vm_selection;}

  /// see VariableMetric
  int add_variable_metric(VariableMetric& H,
			  CH_Matrix_Classes::Integer y_id,
			  const CH_Matrix_Classes::Matrix& y,
			  bool descent_step,
			  CH_Matrix_Classes::Real weightu,
			  CH_Matrix_Classes::Real model_maxviol,
			  const CH_Matrix_Classes::Indexmatrix* indices=0);

  /// propagates the call to QPSolverObject::QPstart_modification() of the current qpsolver
  GroundsetModification* start_modification()
  { return qpsolver->QPstart_modification();}

  /// change the groundset description as specified by the argument 
  int apply_modification(const GroundsetModification& mdf);


  /// m-file output routine for debugging or testing in Matlab (not yet working)
  int mfile_data(std::ostream& out) const;

  ///output settings
  void set_cbout(const CBout* cb, int incr=-1)
  {
    VariableMetricModel::set_cbout(cb,incr);
    qpsolver->set_cbout(this);
    uqpsolver.set_cbout(this);
  }

};

/** 
 * @}
*/



}

#endif

