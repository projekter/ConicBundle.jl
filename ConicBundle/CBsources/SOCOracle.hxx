/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/SOCOracle.hxx
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



#ifndef CONICBUNDLE_SOCORACLE_HXX
#define CONICBUNDLE_SOCORACLE_HXX

/**  @file SOCOracle.hxx
    @brief Header declaring the classes ConicBundle::SOCOracle and ConicBundle::SOCBundleParameters  (needed for ConicBundle::SOCModel)
    @version 1.0
    @date 2017-09-29
    @author Christoph Helmberg
*/

//------------------------------------------------------------

#include "MatrixCBSolver.hxx"
#include "matrix.hxx"

//------------------------------------------------------------


namespace ConicBundle {

/**@defgroup abstract_soc_oracle abstract second order cone oracle
  @brief for minimizing the support function over the second order cone with \f$x_0=1\f$ for an affine cost function or, equivalently, Lagrangian relaxation of Linear Programs over the second order cone 

   For a vector \f$c\in\mathbf{R}^n\f$ and a matrix \f$A\in\mathbf{E}^{n\times n}\f$ (explicitly or implicitly given), the abstract SOCOracle represents one of the
  following convex functions \f$\mathbf{R}^m\to\mathbf{R}\f$ depending on the
  selected parameter ConicBundle::FunctionTask and the scalar factor
  \f$\gamma>0\f$ when adding the function to the problem by
  MatrixCBsolver::add_function(),

   \f[ f_\gamma(y)=\gamma \max\{{1 \choose \bar x}^\top(c+Ay)\colon 1\ge\|\bar x\|\}\f]
   \f[ f^+_\gamma(y)=\gamma \max\{0,{1 \choose \bar x}^\top(c+Ay)\colon 1\ge\|\bar x\|\}\f]
   \f[ f^+_\infty(y)=\gamma \max\{0,{1 \choose \bar x}^\top(c+Ay)\colon 1\ge\|\bar x\|\} \mbox{ with }\gamma\to \infty. \f]

   While \f$\gamma\f$ is the fixed value given explicitly in
   MatrixCBsolver::add_function() in the first two cases, in the last case only
   the initial value is supplied and then solver and model successively
   determine new values so as to ensure (if feasible) that in the end
   solutions \f$y\f$ (almost) satisfy \f$\max\{{1 \choose \bar x}^\top(c+Ay)\colon 1\ge\|\bar x\|\}\le 0\f$.

   The solver exploits a SOCOracle by using a specialized SOCModel which
   implements a cutting model similar to the spectral bundle method, i.e., the
   cutting surface is determined by optimizing over a succsively updated
   subspace for \f$\bar x\f$. The approach is explicitly intended for large scale
   applications, therefore SOCOracle never supplies the coefficient matrices
   explicitly but provides access to computing their action on given vectors
   or matrices and to the evaluation of the function for given \f$y\f$. In
   this sense it is a matrix free oracle. For an implementation that supports
   sparse matrices see \ref implemented_soc_oracle. In special situations it
   may also be worth to implement a corresponding specialized class derived
   from SOCOracle.
  
   If in setting up groundset within MatrixCBSolver the \f$y\f$ variables are
   introduced with a linear cost vector \f$b\f$ and appropriate sign
   constraints, the three functions above may be interpreted as the Lagrangian
   dual to a primal semidefinite program of the form

   \f[ \mbox{maximize }c^\top x\rangle 
       \mbox{ subject to } 
       A^\top x\begin{array}{c}\le\\=\\\ge\end{array}b,
       x_0\begin{array}{c}\le\\=\end{array}\gamma, 
    ~~ x_0\ge \|\bar x\| \mbox{ with }x={x_0\choose\bar x}\f]

   The constraint on \f$x_0\f$ is always
   added implicitly and its type is guided by the choice of
   ConicBundle::FunctionTask in MatrixCBSolver::add_function(). In particular
 
   - ObjectiveFunction represents the constraint \f$x_0=\gamma\f$.
     Lagrangian relaxation of the constraints \f$A^\top x\f$ by
     multipliers \f$y\f$ then results in the dual function
     \f$f_\gamma\f$. Note, if the other constraints already imply
     \f$x_0 =\gamma\f$, the dual solution \f$y\f$ will not be unique,
     because adding a constant to the first component of \f$ c+Ay\f$
     will give rise to the same objective value.

   - ConstantPenaltyFunction represents the constraint
     \f$x_0\le\gamma\f$.  Lagrangian relaxation of the constraints
     \f$A^\top x\f$ by multipliers \f$y\f$ results in the dual
     function \f$f_\gamma^+\f$. It corresponds to a primal feasible
     set having $x_0$ bounded by \f$\gamma\f$.  If a bound on the size
     of the primal optimal solution is known in advance, this is a
     good choice for \f$\gamma\f$, and this variant is then more
     suitable than the next and last one.

   - AdaptivePenaltyFunction represents the constraint
     \f$x_0\le\gamma\f$ for \f$\gamma \to \infty\f$, so in fact this only
     leads to \f$x_0\ge 0\f$.  Lagrangian relaxation of the constraints 
     \f$A^\top x\f$ by multipliers \f$y\f$ results in
     the dual function \f$f_\infty^+\f$.  The given initial value of
     \f$\gamma\f$ will be increased by SOCModel whenever there is no
     progress in forcing the support function value of the affine cost
     function to a value of at most 0. If the primal problem has a
     bounded optimal solution, this will end at some point, but
     \f$\gamma\f$ might get too large for keeping numerical
     stability. In general this option is still a bit hazardous, the
     other two variants are clearly to be preferred. In particular, if
     a good bound on reasonable sizes of optimal solutions are known,
     it is suggested to choose ConstantPenaltyFunction instead.

   In terms of the primal SOC view the sense of the constraints
   \f$A_{\bullet,i}^\top x \begin{array}{c}\le\\=\\\ge\end{array} b_i\f$ is
   determined by the sign constraints on the dual variables \f$y_i\f$ and
   these must be encoded in the ConicBundle::Groundset of the y variables,
   these bounds and the right hand side are _not_ provided by this oracle.
   The right hand side can be added e.g. as a linear cost term for the
   variables \f$y\f$ or alternatively via an AffineFunctionTransformation.

   The bundle method computes the cutting model by successively
   computing and collecting extreme rays of the second order cone with
   \f$x_0=1\f$ in addition to the unit vector having \f$x_0=1\f$ and
   \f$\bar x=0\f$. The proper convex/conic combination of these
   (computed by the method) successively approximates an optimal solution
   \f$x\f$ to the primal SOC. As this primal approximation is often
   the central object of interest (e.g. it is the basis for cutting
   plane approaches) the method can be instructed to collect this
   primal information explicitly by supplying some MatrixPrimal in the
   implementation of SOCOracle::generate_minorant().  

    For a concrete implementation of SOCOracle using \ref matrixclasses see
    \ref implemented_soc_oracle. 

    Internally, the SumBlockModel implementation of the bundle model
    is SOCModel.
*/
//@{


  /** @brief Interface for extending PrimalData, e.g., in Lagrangian
      relaxation of column generation approaches

      This object has to be created and returned in
      SOCOracle::evaluate if in the course of evaluating the oracle
      one notices, that additional primal variables are needed and the
      old primal variables need to be updated accordingly.

      The object will be deleted by ConicBundle after use.

  */

  class SOCPrimalExtender: public PrimalExtender
  {
  public:
    ///
    virtual ~SOCPrimalExtender(){}

    /// like in PrimalExtender, called by ConicBundle to update internal PrimalData objects, has to return 0 on success 
    virtual int extend(PrimalData&)=0;

    /// called by ConicBundle to update internal SOC vectors, has to return 0 on success 
    virtual int extend_SOC(CH_Matrix_Classes::Matrix& /* vectors */ )=0;
  };


  /**@brief Oracle interface for minimization of the support function over the seoncd order cone with \f$x_0=1\f$ for an affine cost function or, equivalently, Lagrangian relaxation of second order cone programs

   Within the setting explained in \ref abstract_soc_oracle the
   abstract class SOCOracle defines a matrix free interface to the 
   support function and the action of the affine cost function
   \f$c+Ay\f$ for the second order bundle cutting model
   implemented in SOCModel. In particular, it provides in
   
   - generate_minorant() the minorant corresponding to a SOC element 
     (with primal data for primal aggregation)

   - extract_SOCvector() returns for a minorant the SOC vector that gives
     rise to this minorant (for this all minorants must include this primal data)
   
   - projection() the projection of the affine cost function data 
     onto a subspace of \f$\bar x\f$ spanned by \f$P\f$ 

   - evaluate() for given \f$y\f$ the function value for the 
     affine cost vector or a SOC vector with sufficiently high SOC_value 
     to exceed the threshold for null steps in the bundle method 

   - evaluate_projection() for given \f$y\f$ the same as in evaluate
     but now for the projected affine matrix function where \f$\bar x\f$
     is restricted to the span of \f$P\f$.
 
   - and the function check_correctness() with the same meaning as in
     ConicBundle::FunctionOracle
   
   For a concrete implementation of SOCOracle see ConicBundle::AffineMatrixFunction .  
  */

  class SOCOracle: public ModifiableOracleObject 
  {
  public:

    /** @brief generates the minorant that arises from the given second order
	cone vector. The minorant must allow some primal aggregation.
	
	The returned minorant must satisfy
	offset_gives_value_at_origin()==true.  The minorant is passed
	over to the caller and will be deleted there.  From a convex
	combination of such minorants the routine extract_SOCvec()
	must be able to reproduce a suitable SOCvec, see there.
    */
    virtual Minorant* generate_minorant(const CH_Matrix_Classes::Matrix& SOCvec)=0;

    /** @brief given a minorant arising from convex combinations of minorants 
	produced by generate_minorant() it determines a \f$(x_0,\bar x)\f$
	generating it and returns SOCvec\f$=\bar x/\|\bar x\|\f$ or
	some unit vector if \f$\bar x=0\f$.
	
	This routine is needed when returning from a polyhedral model
	used by a SumBundleHandler in order to reconstruct an aggregate 
        direction for the conic model. Note that generate_minorant has
        to provide suitable primal information for this, so that the
        convex combinations of the minorants allow this reconstruction. 
    */
    virtual int extract_SOCvector(CH_Matrix_Classes::Matrix& SOCvec, const Minorant* SOCminorant)=0;
    
    /** @brief compute projected coefficient matrices

        The projection matrix is \f$P=\left[\begin{array}{cc}1 & 0 \\ 0 & \bar
	P \end{array}\right]\f$. offset\f$=c^\top P\f$ always has to be
	computed.  If index_subset==0, row i of coeffs is filled with row i of
	\f$A^\top P\f$ for i=1,dots,m If index_subset!=0, row i of coeffs holds the row
	(*index_subset)(i)
    */
    virtual int projection(
			   CH_Matrix_Classes::Matrix& offset,
			   CH_Matrix_Classes::Matrix& coeffs,
			   const CH_Matrix_Classes::Matrix& bar_P,
			   const CH_Matrix_Classes::Indexmatrix* index_subset=0) =0;
    
    /** computes \f$\max\{{1 \choose \bar x}^\top(c+Ay)\colon 1\ge\|\bar x\|\}\f$ and
        returns the value and a maximizing vector \f${1 \choose \bar x}\f$.
	@return 0 on success */
    virtual
    int
    evaluate
    (
     /// argument = current variables = position where to evaluate the function
     const CH_Matrix_Classes::Matrix& current_point,
     /** The columns of the matrix are orthonormal bundle vectors and 
         may help to construct good startig vectors for the eigenvalue 
         computation by iterative methods but have no other use. 
         bundlevecs may have dimension 0. */
     const CH_Matrix_Classes::Real relprec,
     /** on input: gives the threshold for a null step; a vector is good enough to
         yield sufficient improvement if its value exceeds 
         this SOC_value. 

         on output: the value obtained for the returned SOC_vector; if it is smaller than the input SOC_value, then ist must be guaranteed to lie within 
         relprec*(max(abs(SOC_value),1.)) of the true optimum  */
     CH_Matrix_Classes::Real& SOC_value,
     /**  @param[out] SOC_vector return a maximizing vector \f${1\choose \bar x}\f$ */
     CH_Matrix_Classes::Matrix& SOC_vector,
     /** @param[out] primal_extender (PrimalExtender*&) if primal_data or
           vectors provided in (minonrants of) previous calls now have to
           be updated due to changes in the primal problem -- e.g., this may
           happen in column generation -- one may return a pointer to a
           SOCPrimalExtender object on the heap.  This object will be used
           by ConicBundle to update all its internally stored vectors and
           primal_data objects in its minorants by calling
           PrimalExtender::extend on each of these.  Afterwards ConicBundle
           deletes primal_extender.  If this is not needed, set the variable
           to NULL. */
     SOCPrimalExtender*& primal_extender 
     )
     = 0;

    /** compute SOC maximum  for the restricted \f$\bar x\f$ subspace spanned by P
	@return 0 on success */
    virtual
    int
    evaluate_projection
    (
     /// argument = current variable values = position where to evaluate the function
     const CH_Matrix_Classes::Matrix& current_point,
     /** orthogonal matrix defining the projection; */
     const CH_Matrix_Classes::Matrix& P, 
     /** relative precision requirement */
     const CH_Matrix_Classes::Real relprec,
     /** output: ortho */
     CH_Matrix_Classes::Real& projected_SOC_value
     )
     = 0;


    /**@brief This routine need not be implemented unless variables
      (constraints in Lagrangean relaxation) are added or deleted on
      the fly

      The routine is only called by the solver if the variables indeed get
      modified by the solver or a modification is passed on by the user via
      the solver interface. @a oracle_modification is then used to either
      transfer user supplied instructions to the oracle on how to modify
      itself or to inform the oracle about changes the solver was asked to
      perform on the variables. If available, the solver will also show the
      effect of these changes on the center point in @a new_center and @a
      old_center; if these are not available then they hold NULL. A user
      supplied @a oracle_modification will be checked for consistency with the
      actual changes in the variables and mismatches will cause failures.

      The remaining variables are output variables by which the oracle tells
      the solver which information has a chance to be preserved in view of
      these changes.  If e.g. the deletion of some nonzero variables
      invalidates the function value in the new center, the oracle has to set
      discard_objective_in_center=true.  If the entire model cannot be
      preserved (this includes the aggregates and the function values), the
      oracle needs to set discard_model=true; If only aggregate minorants
      cannot be preserved, the oracle needs to set discard_aggregates=true; in
      the current implementation of SOCModel this removes all minorants
      generated by generate_minorant() because each of them is typically 
      generated by more than one Ritz vector. Whenever new variables were
      added, the model can only be preserved if the remaining minorants (maybe
      without aggregates) can be extended for these new variables. In this
      case the oracle has to supply the appropriate MinorantExtender via @a
      minorant_extender and only those minorants will be kept for which this
      operation succeeds.
	      
      Return value 0 indicates that these actions allow to continue without
      errors, other return values result in an overall error on these changes.
    */
    virtual 
    int 
    apply_modification
    (
     const OracleModification& /* oracle_modification */,
     const CH_Matrix_Classes::Matrix* /* new_center */,
     const CH_Matrix_Classes::Matrix* /* old_center */,
     bool& /* discard_objective_in_center */,
     bool& /* discard_model */, 
     bool& /* discard_aggregates */,
     MinorantExtender*& /* minorant_extender */
     )
    {return 1;}



    /**@brief switch on/off some correctnes checks on the oracle */
    virtual
    bool
    check_correctness() const
    {return true;}
      

  };



  /**@brief Bundle parameters for SOCModel

   The current bundle update routine implemented in SOCModel 
   may be controlled to some extent by passing SOCBundleParameters,
   but the update is well tuned already and the parameters are
   used in a quite different meaning than for the usual polyhedral
   model,so it seems better not to meddle with them.

   The old parameters inherited from BundleParameters have the following
   meaning.
    
   BundleParameters::max_model_size gives an upper bound on the subspace dimension,
   which is at least three here 
   
   BundleParameters::max_bundle_size gives an upper bound on the number of 
   second order cone vectors that may stored for selecting 
   the subspace 
   
  */ 

  class SOCBundleParameters: public BundleParameters
  {
  public:
   //int max_model_size;       //inherited   
   //int max_bundle_size;       //inherited

       
    //currently no further items needed

  };

//@}

}
#endif

