/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/PSCOracle.hxx
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



#ifndef CONICBUNDLE_PSCORACLE_HXX
#define CONICBUNDLE_PSCORACLE_HXX

/**  @file PSCOracle.hxx
    @brief Header declaring the classes ConicBundle::PSCOracle and ConicBundle::PSCBundleParameters  (needed for ConicBundle::PSCModel)
    @version 1.0
    @date 2017-09-29
    @author Christoph Helmberg
*/

//------------------------------------------------------------

#include "MatrixCBSolver.hxx"

//------------------------------------------------------------


namespace CH_Matrix_Classes {
  class Indexmatrix;
  class Matrix;
}


namespace ConicBundle {

/**@defgroup abstract_psc_oracle abstract positive semidefinite cone oracle
  @brief for the minimization of the maximum eigenvalue of affine matrix functions (see also \ref implemented_psc_oracle) or, equivalently, Lagrangian relaxation of Linear Programs over the cone of positive semidefinite matrices 

   The oracle for the positive semidefinite cone is based on the
   observation, that the support function over the positive 
   semidefinite matrices with trace 1 is the maximum eigenvalue 
   function,
   
   \f[ \lambda_{\max}(C)=\max\{\langle C,X\rangle\colon \langle I,X\rangle=1,X\succeq 0\}\f]

   For (explicitly or implicitly given) symmetric matrices \f$C\f$ and
  \f$A_i\f$ of order \f$n\f$ for \f$i=1,\dots,m\f$, the abstract
  PSCOracle represents one of the following convex functions
  \f$\mathbf{R}^m\to\mathbf{R}\f$ depending on the selected parameter
  ConicBundle::FunctionTask and the scalar factor \f$\gamma>0\f$ when
  adding the function to the problem by
  MatrixCBsolver::add_function(),

   \f[ f_\gamma(y)=\gamma \lambda_{\max}(C+\sum y_iA_i), \f]
   \f[ f^+_\gamma(y)=\gamma \max\{0,\lambda_{\max}(C+\sum y_iA_i)\}. \f]
   \f[ f^+_\infty(y)=\gamma \max\{0,\lambda_{\max}(C+\sum y_iA_i)\} \mbox{ with }\gamma\to \infty. \f]

   While \f$\gamma\f$ is the fixed value given explicitly in
   MatrixCBsolver::add_function() in the first two cases, in the last case only
   the initial value is supplied and then solver and model successively
   determine new values so as to ensure (if feasible) that in the end
   solutions \f$y\f$ (almost) satisfy \f$\lambda_{\max}(C+\sum y_iA_i)\le
   0\f$.

   The solver exploits a PSCOracle by using a specialized PSCModel which
   implements the cutting model of the spectral bundle method, i.e., the
   cutting surface is determined by an approximation of the active eigenspace
   of the maximum eigenvalue successively constructed by requesting maximum
   eigenvalue/eigenvector evaluations via PSCOracle::evaluate(). The
   approach is explicitly intended for large scale matrices, therefore
   PSCOracle never supplies the coefficient matrices explicitly but
   provides access to computing their action on given vectors or matrices and
   to the evaluation of the mamximum eigenvalue/eigenvector for given
   \f$y\f$. In this sense it is a matrix free oracle. For an implementation
   that supports various special kinds of matrices (e.g. sparse and low rank
   matrices) see \ref implemented_psc_oracle. Of course special situations may
   allow for significantly more efficient eigenvalue/eigenvector computations
   than the general purpose iterative Lanczos method included there. In this
   case it is recommended to implement a corresponding specialized class
   derived from PSCOracle.
  
   If in setting up groundset within MatrixCBSolver the \f$y\f$ variables are
   introduced with a linear cost vector \f$b\f$ and appropriate sign
   constraints, the three functions above may be interpreted as the Lagrangian
   dual to a primal semidefinite program of the form

   \f[ \mbox{maximize }\langle C,X\rangle 
       \mbox{ subject to } 
       \langle -A_i,X\rangle \begin{array}{c}\le\\=\\\ge\end{array} b_i,
    ~~ \langle I,X\rangle \begin{array}{c}\le\\=\end{array}\gamma,
    ~~ X\succeq 0 \f]

   The constraint \f$\langle I,X\rangle\f$ on the trace is always
   added implicitly and its type is guided by the choice of
   ConicBundle::FunctionTask in MatrixCBSolver::add_function(). In particular
 
   - FunctionTask::ObjectiveFunction represents the constraint \f$X\succeq 0,
     \langle I,X\rangle=\gamma\f$.  Lagrangian relaxation of the constraints
     on \f$\langle -A_i,X\rangle\f$ by multipliers \f$y_i\f$ then results in
     the dual function \f$f_\gamma\f$. It corresponds to a primal feasible set
     having constant trace \f$\gamma\f$ Note, if the other constraints already
     imply \f$\langle I,X\rangle =\gamma\f$, the dual solution \f$y\f$ will
     not be unique, because adding a constant to the diagonal of \f$ C+\sum
     y_iA_i\f$ will give rise to the same objective value.

   - FunctionTask::ConstantPenaltyFunction represents the constraint
     \f$X\succeq 0, \langle I,X\rangle\le\gamma\f$.  Lagrangian relaxation of
     the constraints on \f$\langle -A_i,X\rangle\f$ by multipliers \f$y_i\f$
     results in the dual function \f$f_\gamma^+\f$. It corresponds to a primal
     feasible set having bounded trace \f$\gamma\f$.  If a bound on the trace
     of the primal optimal solution is known in advance, this is a good choice
     for \f$\gamma\f$, and this variant is then more suitable than the next
     and last one.

   - FunctionTask::AdaptivePenaltyFunction represents the constraint
     \f$X\succeq 0, \langle I,X\rangle\le\gamma\f$ for \f$\gamma \to
     \infty\f$, so in fact only \f$X\succeq 0\f$.  Lagrangian relaxation of
     the constraints on \f$\langle -A_i,X\rangle\f$ by multipliers \f$y_i\f$
     results in the dual function \f$f_\infty^+\f$.  The given initial value
     of \f$\gamma\f$ will be increased by PSCModel whenever there is no
     progress in forcing the maximum eigenvalue of the affine matrix function
     to a value of at most 0. If the primal problem has a bounded optimal
     solution, this will end at some point, but \f$\gamma\f$ might get too
     large for keeping numerical stability. In general this option is still a
     bit hazardous, the other two variants are clearly to be preferred. In
     particular, if a good bound on reasonable sizes of optimal solutions are
     known, it is suggested to choose ConstantPenaltyFunction instead.

   In terms of the primal SDP view the sense of the constraints \f$\langle
   -A_i,X\rangle \begin{array}{c}\le\\=\\\ge\end{array} b_i\f$ is determined by
   the sign constraints on the dual variables \f$y_i\f$ and these must be
   encoded in the ConicBundle::Groundset of the y variables, these bounds and
   the right hand side are _not_ provided by this oracle.  The right hand side
   can be added e.g. as a linear cost term for the variables \f$y\f$ or
   alternatively via an AffineFunctionTransformation.

   It is convenient to summarize the effect of the \f$A_i\f$ by introducing the mapping 
   (\f$\mathbf{S}^n\f$ refers to the set of real symmetric matrices of order \f$n\f$)

   \f[
      \mbox{opA}\colon \mathbf{S}^n\to\mathbf{R}^m, X\mapsto 
      \left[\begin{array}{@{}c@{}}
      \langle A_1,X\rangle\\[-1ex]
      \vdots\\[-.5ex]
      \langle A_m,X\rangle\end{array}\right]
   \f]

   The name \f$\mbox{opA}\f$ or of its adjoint 

   \f[
   \mbox{opAt}\colon \mathbf{R}^m\to\mathbf{S}^n, y\mapsto \sum_{i=1}^m y_iA_i 
   \f]

   will appear frequently in the names of 
   the members of PSCOracle and its derived classes.     

   The bundle method computes the cutting model by successively
   computing and collecting eingenvectors to the maximum eigenvalue of
   the maximum eigenvalue function. The proper convex/conic
   combination of these, not needed explicitly yet implicitly computed
   by the method, successively approximates an optimal solution
   \f$X\f$ to the primal SDP. As this primal approximation is often
   the central object of interest (e.g. it is the basis for cutting
   plane approaches) the method can be instructed to collect this
   primal information explicitly by supplying some PSCPrimal (similar
   to MatrixPrimal in MatrixFunctionOracle) in the
   implementation of PSCOracle::generate_minorant().  As the bundle
   method is designed for large scale matrices, collecting the full
   matrices (implemented in DensePSCPrimal) might be computationally
   too expensive, but other variants like collecting only the entries
   on some sparse support (e.g. in SparsePSCPrimal), sparse combined
   with a low rank dense part (e.g. in GramSparsePSCPrimal), or
   collecting each of these just on blocks (BlockPSCPrimal) can be
   realized by proper implementations of PSCPrimal.

    For a concrete implementation of PSCOracle using \ref matrixclasses see
    \ref implemented_psc_oracle. 

    Internally, the SumBlockModel implementation of the spectral bundle model
    is PSCModel.
*/
//@{


  /** @brief Interface for extending PrimalData, e.g., in Lagrangian relaxation of column generation approaches

      This object has to be created and returned in PSCOracle::evaluate
      if in the course of evaluating the oracle one notices, that additional primal
      variables are needed and the old primal variables need to be updated accordingly.

      The object will be deleted by ConicBundle after use.

  */

  class PSCPrimalExtender: public PrimalExtender
  {
  public:
    ///
    virtual ~PSCPrimalExtender(){};

    /// like in PrimalExtender, called by ConicBundle to update internal PrimalData objects, has to return 0 on success 
    virtual int extend(PrimalData&)=0;

    /// called by ConicBundle to update internal Ritz_vectors, has to return 0 on success 
    virtual int extend_Ritz(CH_Matrix_Classes::Matrix& /* Ritzvecs */ )=0;
  };


  /**@brief Oracle interface for minimization of the maximum eigenvalue of an affine matrix function or, equivalently, Lagrangian relaxation of semidefinite programs

   Within the setting explained in \ref abstract_psc_oracle the
   abstract class PSCOracle defines a matrix free interface to maximum
   eigenvalue/eigenvector and the action of the affine matrix function
   \f$C+\sum_{i=1}^my_iA_i\f$ for the spectral bundle cutting model
   implemented in PSCModel. In particular, it provides in
   
   - generate_minorant() the minorant corresponding to a gram matrix PP^T
     (maybe with additional primal data for primal aggregation)

   - svec_projection() the projection \f$P^T*A*P\f$ onto a subspace spanned by
     \f$P\f$ for matrices \f$C\f$ and \f$A_i\f$

   - evaluate() for given \f$y\f$ the maximum eigenvalue/eigenvector of the 
     affine matrix function or a Ritz_vector with sufficiently high Ritz_value 
     to exceed the threshold for null steps in the bundle method 

   - evaluate_projection() for given \f$y\f$ the same as in evaluate but now 
     for the projected affine matrix function \f$P^T(C+\sum y_iA_i)P\f$ where 
     the columns in $P$ span the desired subspace 
 
   - left_right_prod() the product \f$G=P^T*(A_i)*Q\f$ for matrices \f$P\f$
     and \f$Q\f$ (required for the bundle update and scaling heuristics in
     PSCModel)

   - and the function check_correctness() with the same meaning as in
     ConicBundle::FunctionOracle
   
   For a concrete implementation of PSCOracle see ConicBundle::AffineMatrixFunction .  
  */

  class PSCOracle: public ModifiableOracleObject 
  {
  public:

    /** @brief generates the minorant that arises from the gram matrix PP^T (maybe including some primal information for primal aggregation). 
	
	The returned minorant must satisfy offset_gives_value_at_origin()==true.
	The minorant is passed over to the caller and will be deleted there. 
    */
    virtual Minorant* generate_minorant(const CH_Matrix_Classes::Matrix& P)=0;
    
    /** @brief compute svec representation of the projected coefficient matrices

	This always includes the constant cofficient matrix, svec_offset=svec(P^T*C*P).
	If index_subset==0, row i of svec_coeffs is filled with svec(P^T*A_i*P)^T
	for i=1,dots,m
	If index_subset!=0, row i of svec_coeffs holds this for A_{(*index_subset)(i)}
    */
    virtual int svec_projection(
				CH_Matrix_Classes::Matrix& svec_offset,
				CH_Matrix_Classes::Matrix& svec_coeffs,
				const CH_Matrix_Classes::Matrix& P,
				const CH_Matrix_Classes::Indexmatrix* index_subset=0) =0;
    
    /** compute maximum eigenvalue of the matrix C+opAt(y)
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
     const CH_Matrix_Classes::Matrix& bundlevecs, 
     /** relative precision requirement for objective values
	 leading to descent steps */
     const double relprec,
     /** gives the threshold for a null step; a vector is good enough to
         yield sufficient improvement if its Ritz value exceeds 
         this Ritz_bound. Otherwise the largest of the returned
         Ritz_values must be guaranteed to lie within 
         relprec*(max(abs(Ritz_bound),1.)) of the maximum eigenvalue */
     const double Ritz_bound,
     /** on input: if not zero dimensional, the columns of the matrix contains
         the orthonormal family of Ritz_vectors returned in the previous call;
         these may help to construct good starting vectors for the eigenvalue
         computation by iterative methods.

         on output: the columns of the matrix form an orthonormal family
         that has to contain at least one vector \f$v\f$ with the following
         property: \f$v\f$ is an eigenvector to the maximum eigenvalue or it
         has a Ritz value \f$v^T(C-\sum y_iA_i)v\f$ that exceeds the 
         null step bound @a Ritz_bound. */
     CH_Matrix_Classes::Matrix& Ritz_vectors,
     /** Ritz_values corresponding to the Ritz_vectors */
     CH_Matrix_Classes::Matrix&  Ritz_values,
     /** @param[out] primal_extender (PrimalExtender*&) if primal_data or
           Ritz_vectors provided in minonrants of previous calls now have to
           be updated due to changes in the primal problem -- e.g., this may
           happen in column generation -- one may return a pointer to a
           PSCPrimalExtender object on the heap.  This object will be used
           by ConicBundle to update all its internally stored Ritz_vectors and
           primal_data objects in its minorants by calling
           PrimalExtender::extend on each of these.  Afterwards ConicBundle
           deletes primal_extender.  If this is not needed, set the variable
           to NULL. */
     PSCPrimalExtender*& primal_extender 
     )
     = 0;

    /** compute maximum eigenvalue of P^T(C-opAt(y))P
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
     const double relprec,
     /** output: orthognal matrix that contains at least one column 
         vector that is an eigenvector to the maximum eigenvalue 
         of the projected Matrix */
     CH_Matrix_Classes::Matrix& projected_Ritz_vectors,
     /** output: Ritz_values corresponding to the Ritz_vectors */
     CH_Matrix_Classes::Matrix& projected_Ritz_values  
     )
     = 0;



    /** fills row i with vec(P^T*(A_i)*Q)^T  (needed for the bundle update and the scaling heuristic in PSCModel)
	@return 0 on success
    */
    virtual
    int 
    left_right_product
    (
     /// index of the coefficient matrix 
     int i,
     /// left multiplication matrix 
     const CH_Matrix_Classes::Matrix& P,
     /// right multiplication matrix
     const CH_Matrix_Classes::Matrix& Q,
     /// the result is stored in this matrix. 
     CH_Matrix_Classes::Matrix& G
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
      the current implementation of PSCModel this removes all minorants
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



  /**@brief Bundle parameters for PSCModel, we recommend no to modify them

   The current bundle update routine implemented in PSCModel 
   may be controlled to some extent by passing PSCBundleParameters,
   but the update is well tuned already and the parameters are
   used in a quite different meaning than for the usual polyhedral
   model, so it seems better not to meddle with them.

   The inherited BundleParameters are potentially used as 
   parameters for a SumBundle representation (if employed)

  */ 

  class PSCBundleParameters: virtual public BundleParameters
  {
  public:

    /** gives a lower bound on the subspace dimension
	that is added to the active subspace identified by Tapia indicators
	for the bundle subspace used in the model
    */
    int psc_model_size;       
   
    /** gives an upper bound on the number of old Ritz-vectors 
	that may stored for generating second order information
    */
    int psc_bundle_size;      
    
    /// gives an upper bound on the number of new Ritz-vectors added in each bundle update
    
   int psc_new_subgradients;  

    /** minimum size of the "next best subspace" in addition to the bundle
	that is kept for generating good scaling information 
    */
    int psc_keep;

    /// maximum number of aggregate matrices allowed (currently the choice is at most one)
    int psc_aggregates;

    /** relative precision tolerance for including Ritz vectors in the active dimension if their value does not differ by more than this relative precision*/
    double psc_tolerance;

    /// selection parameter in case several update rules are available for the semidefinite model
    int psc_update_rule;

    /**@brief initialize to given values */
    virtual int init(const BundleParameters& bp);

    /// default constructor
    PSCBundleParameters():
      BundleParameters(),psc_model_size(3),psc_bundle_size(100),
      psc_new_subgradients(5),psc_keep(3),
      psc_aggregates(1),psc_tolerance(1e-6),psc_update_rule(0)
    {}
    
   /// "copy" constructor
   PSCBundleParameters(const BundleParameters& bp):
      BundleParameters(),psc_model_size(3),psc_bundle_size(100),
      psc_new_subgradients(5),psc_keep(3),
      psc_aggregates(1),psc_tolerance(1e-6),psc_update_rule(0)
    {init(bp);}
    
    ///
    ~PSCBundleParameters();

    };

//@}

}
#endif

