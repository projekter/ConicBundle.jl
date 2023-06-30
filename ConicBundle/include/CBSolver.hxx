/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  include/CBSolver.hxx
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



#ifndef CONICBUNDLE_CBSOLVER_HXX
#define CONICBUNDLE_CBSOLVER_HXX

/**  @file CBSolver.hxx
    @brief Header declaring the classes ConicBundle::CBSolver, ConicBundle::FunctionOracle and ConicBundle::PrimalData
    @version 1.0
    @date 2005-03-01
    @author Christoph Helmberg
*/

#include <assert.h>
#include <iostream>
#include <vector>

//------------------------------------------------------------

/**@brief   conic bundle method solver for sum of convex functions. See the \ref ConicBundle_Manual for a quick introduction.
	@author  C. Helmberg
*/	
namespace ConicBundle {

/**@defgroup cxxinterface Interface to ConicBundle for the Language C++
   @brief Solve \f$min_{y\in\mathbf{R}^m}  f_0(y) + f_1(y) + ... + f_k(y)\f$ 
   for convex functions f_i, the y-variables may be bounded or 
   box constrained. The most important steps are explained here

   <b>Setting up the Problem, the Functions, and the Main Loop</b>

   First create a new problem/solver ConicBundle::CBSolver, let us
   call it solver for brevity.  In invoking the constructor a boolean
   flag may be set to enforce the use of a minimal bundle model that
   employs only one common aggregate and one common subgradient for
   all functions, so basically "no bundle", which may be favorable if
   fast iterations and/or little memory consumption are essential.
 

   Next, set the dimension of the design
   variables/argument as well as possible box constraints on these by
   ConicBundle::CBSolver::init_problem().

   Now set up each of your functions f_i as a
   ConicBundle::FunctionOracle.  Via the routine
   ConicBundle::FunctionOracle::evaluate() you will supply, for a
   given argument, the function value and a subgradient (=the gradient
   if the function is differentiable) to the solver. The function
   evaluate() is the only function that you definitely have to
   provide, see the miniature example below.

   The function oracles have to be added to the solver 
   using the routine  ConicBundle::CBSolver::add_function().

   Once all functions are added, the optimization process can be
   started. If you know a good starting point then set it with
   ConicBundle::CBSolver::set_new_center_point() now, otherwise the
   method will pick the zero vector or, in the case of box
   constraints, the point closest to zero as starting point.

   Finally, call ConicBundle::MatrixCBSolver::solve() and retrieve 
   ConicBundle::MatrixCBSolver::termination_code() for getting the
   reason for termination. Via parameters in solve() you may also
   tell the solver to return after each descent step or also
   after a maximum number of null steps. This then only interrupts
   computations and calling solve() again continues as if there
   was not break at all. 

   <b>Setting up the Problem, the Functions, and the Main Loop</b>


   After the first call to ConicBundle::CBSolver::solve()
   you can retrieve, at any time, the current objective value by
   ConicBundle::CBSolver::get_objval() and the argument leading to
   this value by ConicBundle::CBSolver::get_center(). For some screen
   output, use ConicBundle::CBSolver::set_out().

   <b>Lagrangean Relaxation, Primal Approximations, and Cutting Planes</b>


   If you are optimizing the Lagrange multipliers of a Lagrangean relaxation, 
   you might be interested in getting an approximation to your primal optimal
   solution. This can be done by specifying in each function
   for each (epsilon) subgradient the corresponding primal
   vectors that generate it, see the parameter primal_solutions in 
   ConicBundle::FunctionOracle::evaluate() as a
   start. Then for each of your functions, you can retrieve
   the current primal approximation using 
   ConicBundle::CBSolver::get_approximate_primal().

   If, in addition, you plan to improve your primal relaxation
   via cutting planes, that are strongly violated by the current
   primal approximation, you should have a look at 
   ConicBundle::CBSolver::append_variables(),
   ConicBundle::FunctionOracle::subgradient_extension() and
   ConicBundle::CBSolver::reinit_function_model(). If you want to
   get rid of primal constraints/dual variables, use
   ConicBundle::CBSolver::get_approximate_slacks()
   and ConicBundle::CBSolver::delete_variables().

   @include cxx_mini_ex.cxx

*/

//@{

  /** @brief serves as the value  "minus infinity", i.e., all bounds <= 
     this value are set to this value and are regarded as minus infinity
   */
  extern const double CB_plus_infinity;
  /** @brief serves as the value  "plus infinity", i.e., all bounds >= 
     this value are set to this value and are regarded as plus infinity
   */
  extern const double CB_minus_infinity;
  
  //------------------------------------------------------------
  
 
  /** @brief A dense vector of double, arguments and subgradients are 
     specified like this 
  */
  typedef std::vector<double> DVector;

  /** @brief A dense vector of int, index vectors for deleting/reorganizing 
     variables are specified like this  
  */
  typedef std::vector<int>    IVector;

  /** @brief In Lagrangean relaxation an approximate primal solution can 
       be generated by supplying primal information derived from this 
       abstract class for each epsilon subgradient within 
       ConicBundle::FunctionOracle::evaluate(). 

      In many applications, e.g. in Lagrangean relaxation, the convex
      minimization problem arises as the dual of a convex primal  
      maximization problem. In this case one is typically interested
      in obtaining a primal approximate solution in addition to the
      dual solution. Under reasonable conditions this is possible if the
      primal solutions that give rise to the subgradients are aggregated 
      along with the subgradients maintained by the bundle algorithm.
      Within the interface the user can achieve this by specifying and 
      returning in the oracle an appropriate derivation of the class 
      ConicBundle::PrimalData. The classical case of a primal 
      ConicBundle::DVector is supplied as an example implementation below. 

      The primal data has to be supplied within ConicBundle::FunctionOracle::evaluate()
      and can be retrieved via the methods 
      ConicBundle::CBSolver::get_approximate_primal() and
      ConicBundle::CBSolver::get_center_primal()
  */

  class PrimalData
  {

  public:

    ///
    virtual ~PrimalData();

    /// Returns a newly generated identical Object
    virtual PrimalData* clone_primal_data() const=0;  

    /// add factor*it to this (types will need to be compatible to this)
    virtual int aggregate_primal_data(const PrimalData& it,double factor=1.)=0;
	
    /// multiply/scale *this with a nonnegative myfactor
    virtual int scale_primal_data(double myfactor)=0;

  };

  /** @brief Interface for extending PrimalData, e.g., in Lagrangian relaxation of column generation approaches

      This object has to be created and returned in FunctionOracle::evaluate
      if in the course of evaluating the oracle one notices, that addional primal
      variables are needed and the old primal variables need to be updated accordingly.

      The object will be deleted by ConicBundle after use.

  */

  class PrimalExtender
  {
  public:
    ///
    virtual ~PrimalExtender();

    /// called by ConicBundle to update internal PrimalData objects, has to return 0 on success 
    virtual int extend(PrimalData&)=0;
  };



  /** @brief serves as the default tolerance for considering minorant entries as zero
   */
  extern const double CB_minorant_zero_tolerance;

  /** @brief serves as the default ratio of nonzeros to dimension for using a sparse
      representatio of a minorant
   */
  extern const double CB_minorant_sparsity_ratio;


/** @brief Each function \f$f\f$ represented by a FunctionModel is
    equipped with a @a function_factor \f$\gamma>0\f$ (it defaults to
    1.) and may be declared as a usual objective function (default) or
    as a penalty function \f$\gamma \max\{f,0\}\f$ with either a
    constant penalty factor or an adaptive penalty factor 
    \f$\gamma\in(0,\infty)\f$.

    @attention AdaptivePenaltyFunction is still very experimental and
    not yet reliable

    - ObjectiveFunction   adds \f$\gamma f\f$

    - ConstantPenaltyFunction  adds \f$\gamma \max\{f,0\}\f$

    - AdaptivePenaltyFunction  adds \f$\gamma \max\{f,0\}\f$ for \f$\gamma\to\infty\f$ 
      so that the result should resemble the indicator of the level set to level 0
      (in exact penalty cases there will be no need to go to infinity) 
 */

  enum FunctionTask {
    ObjectiveFunction=0,              
    ConstantPenaltyFunction=1,
    AdaptivePenaltyFunction=2
  };


  /** @brief basic function object (abstract class). It serves for using the same interface
      on distinct oracle types, but is not yet needed in the standard C++ interface.
 
     Besides the standard ConicBundle::FunctionOracle 
     the interface only accepts a number of prespecified derivations 
     of this class that come along with the CH_Matrix_Classes interface
     (e.g. for semidefinite and second order cones). 
 */
  class FunctionObject
  {
  public:
    ///
    virtual ~FunctionObject();
  };

/** @brief this is used to describe affine minorants of convex functions that will be used for generating cutting models of these functions.
  
    They have to be returned in calls to evaluation oracles like in 
    FunctionOracle::evaluate() or MatrixFunctionOracle::evaluate().

    The Minorant is of the form 
        f(point)=offset+ inner_product(coeff,point),
    offset and coeff must be specified by the user by using the
    constructors or add_offset() and add_coeffs(). The dimension need
    not be specified (it is clear from point) and it is always assumed
    without actual checking that all nonzero coefficients are restricted
    to the dimension. 

    Sometimes it is more natural for the user to specify the offset at
    the point where the evaluation took place, in this case the
    minorant would have the form 
        f(point)=offset+ inner_prodcut(coeff,point-evalpoint).
    If the user wants to avoid computing the offset with respect to
    the origin, it suffices to set offset_at_origin=false at 
    initialization berfore returning the minorant. 

    Once a minorant is returned in an oracle evaluation, control over
    this object is handed over to ConicBundle. ConicBundle may want to 
    and is allowed to clone or modify the minorants and will eventually 
    delete this object. Thus a minorant (and all its associated date)
    must be an object on the heap.
    
    The implementation of this class is given in Minorant.cxx
    
 */

  class Minorant
  {    
  public:


    /** @brief Minorant constructor for a dense subgradient specified via a DVector 

      Suppose in the evaluation of your oracle at the current point \f$y\f$ you determine a subgradient \f$s\f$ for the subgradient inequality

      \f[ f(z)\ge f(y)+\langle s,z-y\rangle\quad\forall z\in\mathbf{R}^m, \f]
      
      then use \f$f(y)\f$ for @a offset and \f$s\f$ in the form of a DVector in @a subg and keep the default value @a offset_at_origin == false.

      If, on the other hand, your oracle implements a support function for some compact set \f$\mathcal{X}\subset\mathbf{R}^m\f$ like
    
      \f[ f(y) = \max_{x\in\mathcal{X}} x^\top y\f]

      then it is actually more efficient to return 0 for @a offset,
      a maximizing \f$x\f$ in @ subg and to put @a offset_at_origin = true.

      You may also pass over a PrimalData object in primal that will be then be
      owned and later deleted by the minorant and aggregated along with the minornat
     */
    Minorant(double offset,const DVector& subg,PrimalData* primal=0,bool offset_at_origin=false);

    /** @brief Minorant constructor for a sparse subgradient where the nonzero values are specified by DVector and the corresponding indices by an IVector 

      Suppose in the evaluation of your oracle at the current point \f$y\f$ you determine a subgradient \f$s\f$ with few nonzeros for the subgradient inequality

      \f[ f(z)\ge f(y)+\langle s,z-y\rangle\quad\forall z\in\mathbf{R}^m, \f]
      
      then use \f$f(y)\f$ for @a offset and pass \f$s\f$ by giving the nonzeros in the DVector @a subg_val, the corresponding indices in a IVector of the same length in @a subg_ind and keep the default value @a offset_at_origin == false.

      If, on the other hand, your oracle implements a support function for some compact set \f$\mathcal{X}\subset\mathbf{R}^m\f$ like
    
      \f[ f(y) = \max_{x\in\mathcal{X}} x^\top y\f]

      then it is actually more efficient to return 0 for @a offset,
      a sparse maximizing \f$x\f$ in @a subg_val and @a subg_ind  and to put @a offset_at_origin = true.

      You may also pass over a PrimalData object in primal that will be then be
      owned and later deleted by the minorant and aggregated along with the minornat
     */
    Minorant(double offset,const DVector& subg_val,const IVector& subg_ind,PrimalData* primal=0,bool offset_at_origin=false);
    
    /** @brief default initializes a zero minorant, see full explanation otherwise; NOTE: if the offset supplied by the minorant refers to the evaluation point (i.e., if it is the function value at the point of evaluation), set @a offset_at_origin = false !

    In many applications, in particular for Lagrangian relaxation, 
    the offset is more naturally given for the origin, and giving it this way
    also reduces computational cost a bit, so @a offset_at_origin = true  
    is the suggested default. If, however, the minorant arises from a subgradient 
    inequality by evaluation in the current point, you may as well give the function 
    value as offset directly and set @a offset_at_origin=false;  
 
    The data specifying the minorant may be set here or (part of) it may be
    entered/added later. The meaning of the other parameters is as follows.

    @a offset gives the constant value (if offset_at_origin = false then relative to
    the evaluation point)

    @a n_elements gives the number of elements specified by @a coeffs
    (and possibly indices), but if @a coeffs == NULL it just asks to reserve space
    for that many coefficients

    If @a coeffs is not NULL, it points to an array of doubles of size
    at least @a n_elements. If @a indices == NULL then coeff[i] belongs to position
    i for i=0 to n_elements-1, otherwise coeff[i] belongs to position indices[i].
    All data is copied, the arrays are not modified, not used later and not
    deleted here.

    The entire input data (offset and coefficients) is multplied by @a scale_val
    to give the final minorant (scale_val is not memorized internally but executed
    immediately).

    If the minorant arises from @a PrimalData and the primal data should be aggregated
    along, it may be entered here or in the routine Minorant::set_primal(). The
    object pointed to is then owned by Minorant and will be deleted by Minorant on its 
    destruction. 

    */
    Minorant(bool offset_at_origin=true,
	     double offset=0.,
	     int n_elementes=0,
	     const double* coeffs=0,
	     const int* indices=0,
             double scale_val=1.,
	     PrimalData* primal=0
	     );

    /// they main purpose of this constructor is to allow easy cloning for derived classes
    Minorant(const Minorant* mnrt,double factor=1.,bool with_primal=false);

    ///
    virtual ~Minorant();

    /// returns the current offset value
    virtual double offset() const;
    /// adds this value to the current offset value
    virtual int add_offset(double value);

    /// allows to specify/modify whether the offset refers to origin or to the point of evaluation at which this minorant was supplied
    virtual bool& offset_gives_value_at_origin();  
    /// true if the offset refers to origin and false, if the offset refers to the value of the  minorant in the point of the oracle evaluation at which this minorant was supplied  
    virtual bool offset_gives_value_at_origin() const;  

    /// returns the value of the coefficient in coordinate i
    virtual double coeff(int i);

    /// adds the value to the coefficient in coordinate i
    virtual int add_coeff(int i,double value);
    
    /// adds the n values (possibly multiplied by factor) to consecutive coefficients starting at start_pos (by default 0); this always converts the minroant into a dense vector first and adds then
    virtual int add_coeffs(int n_elements,const double* values,double factor=1.,int start_pos=0);

    /// adds the n values (possibly multiplied by factor) to the coefficients indicated by indices (if zero, this calls the other add_coeffs for the consecutive version)
    virtual int add_coeffs(int n_elements,const double* values,const int* indices,double factor=1.);

    /// converts to sparse format with zeros of absolut value at most tol*(fabs(offset)+1) if the given ratio of elements is zero in relation to the maximum nonzero index
    virtual int sparsify(double tol=CB_minorant_zero_tolerance,double sparsity_ratio=CB_minorant_sparsity_ratio); 

    /// returns the number of nonzero coefficients
    virtual int nonzeros(); 

    /**@brief returns the list of all nonzero coefficients by returning pointers to the arrays and their common length n_elements 
     
       If on return coeffs==NULL this is interpreted as the zero vector. 
       If on return indices==NULL, the indices are 0 to n_elements-1.
       If on return indices!=NULL, the entries of this array will be sorted in strictly increasing order
    */
    virtual int get_coeffs(int& n_elements,const double*& coeffs,const int*& indices) const;

    /**@brief If the returned pointer is not NULL it gives direct access to the array of current coefficient values with indices 0 up to n_elements-1;

      If the return value is NULL, the representation may not be 
      available or access to the store may not be granted; in this case
      other routines like get_coeffs and add_coeffs have to be used.
      
      This routine is mainly intended for increasing efficiency in
      some internal computations; the validity of the pointer returned
      may get lost with any call to any other routine of this Minorant,
      so during manipulations of the stored values no other routines
      should be called. Needless to say, this routine should only be
      used by experts. 
    */
    virtual double* get_dense_coeff_store(int n_elements);

    /// resorts (and deletes) coefficients so that afterwards it has n_elements and the new coeff(i) has the previous value of coeff(map_to_old_coeff(i))   
      virtual int reassign_coeffs(int n_elements,const int* map_to_old_coeffs);

    /// if the minorant is generated by PrimalData and this should be aggregated along, insert a heap object of it here (see PrimalData)
    virtual void set_primal(PrimalData*);
  
    /// returns NULL if there is no primal data and otherwise a pointer to it
    virtual PrimalData* get_primal();
    
    /// returns NULL if there is no primal data and otherwise a pointer to it (const version)
    virtual const PrimalData* get_primal() const;

    ///generates a full copy (multiplied by factor) on the heap and returns a pointer to it (it also includes a clone of the primal data if with_primal==true, otherwise the copy will have no primal data)
    virtual Minorant* clone_minorant(double factor=1.,bool with_primal=true) const;

    ///adds factor*minorant to this and does this also for the primal if it is availabe
    virtual int aggregate(const Minorant& minorant,double factor=1.);

    ///returns the number of minorants aggregated in this one, value 1 thus means not aggregated
    virtual int number_aggregated() const;

    ///mutliply offset and coefficients (and PrimalData, if given) by scale_val
    virtual int scale_minorant(double scale_val);

    ///return the squared Euclidean norm
    virtual double norm_squared() const;

  protected:
    /// implementation details are hidden on purpose
    class MinorantData; 
    /// implementation details are hidden on purpose
    mutable MinorantData* data;

    Minorant(const Minorant&){}             ///< forbidden, blocked deliberately
    Minorant& operator=(const Minorant&);    ///< forbidden, blocked deliberately
  };

  /** @brief Interface for extending a Minorant, e.g., in Lagrangian Relaxation of cutting plane approaches
 
      Such an object may be returned in a call to apply_modification 
      of an oracle, whenever additional variables are introduced 
      (e.g. as Lagrange multipliers to new constraints) on the fly by some 
      problem modification. It will be this objects task to fill in the new 
      coordinates in minorants, typically by means of PrimalData supplied and 
      updated along with the minorants that were returned by oracle calls,
      e.g. as in FunctionOracle::evaluate().

      The object will be deleted by ConicBundle after use.

         The solver calls this routine whenever new variables have been added on
         the fly in order to extend old subgradients to the new coordinates.
         If primal data was supplied for the subgradients then  
         @a generating_primal holds a pointer to this (possibly aggregated) data,
         otherwise it is NULL. 
 
         In the presence of primal data, the new coordinates correspond to 
         the violation of the new primal constraints. These have to be 
         returned in  @a new_subgradient_values; more precisely, 
         for i=0 to @a varialb_indices.size()-1 the element 
         @a new_subgradient_values[i] has to hold the subgradient information  
         of constraint @a variable_indices[i]; 
 
         If generating_primal is NULL, then the routine can only successfully
         extend the subgradients if the new coordinates have no influence
         on the function; then the new subgradient coordinates are all zero
         and the components of @a new_subgradient_values have to be initialized 
         to zero.

  */

  class MinorantExtender
  {
  public:
    ///
    virtual ~MinorantExtender();

    /** @brief called by ConicBundle to update internal Minorant objects, has to return 0 on success

        @param[in,out] minorant  (Minorant&)
            it holds a (possibly aggregated) minorant that was generated
            from minorants returned by oracle calls, e.g. as in 
	    FunctionOracle::evaluate() If PrimalData was provided in these
	    minorants, this will be aggregated along and will also be 
	    available in this minorant.

        @param[in] n_coords (int)
            the number of coordinate positions that have to be filled in

        @param[out] indices  (const int*)
	    the indices of these coordinate positions (sorted in
	    strictly increasing order)
      
	@return
           -  0 on success,
           -  1 if extension is impossible 
      
     */ 
    virtual int extend(Minorant& minorant,int n_coords,const int* indices)=0;
  };

  /** @brief Base class for informing oracles (or the solver) about dynamic changes in the number and sorting of the variables, if such changes occur at all

      This base class may either be extended by the user to inform
      her/his oracles via the solver on how to adapt dynamically, or
      it is used by the solver to inform an oracle about changes that
      it is was asked to perform on the variables (and only if such
      occured) in order to ensure that the oracle works consistently
      with this.  E.g. for a FunctionOracle the solver uses the
      routine ConicBundle::FunctionOracle::apply_modification() for
      this. Via its return values there the oracle should then tell
      the model/solver how to react to these changes.

      Any changes can always be regrouped so that they result in two steps
      - first appending some new variables (if there are new ones at all)
      - then mapping the (maybe extended) positions to their new places (note,
        this allows to delete variables while clarifying how this is done)

      So the base class only returns this information:
      - the dimension before the changes in get_old_vardim()
      - the dimension after all changes are carried out in get_new_vardim()
      - the number of variables first appended to the old variables in get_appended_vardim()
      - the map from the new indices 0..get_new_vardim()-1 to the old plus appended indices in get_map_to_old_variables()
      
   */
  class OracleModification
  {
  public:
    ///
    virtual ~OracleModification();

    /// returns the old dimension (before the changes occured), a value >=0
    virtual int get_old_vardim() const =0;
    /// returns the current dimension (after the changes occured), a value >=0
    virtual int get_new_vardim() const =0;
    /// returns the number of variables appended, a value >=0
    virtual int get_appended_vardim() const =0;
    /// returns get_new_vardim() many distinct indices out of 0 ... get_old_vardim()+get_appended_vardim()-1, so "old" here refers to old plus appended variables; if this returns 0, there are no rearrangements 
    virtual const int* get_map_to_old_variables() const =0;

    /** @brief incorporate the next_modification in this one so that this one now performs both

       In order to allow collecting several modifications in one, this routine
       must change the current modification so that its application results in
       the equivalent of its current and then the next_modification; this may
       fail if sequence of modification or the implemented types are not
       compatible; 

       Typically the next_modification can only be incorporated if its old
       variable dimension matches the new variable dimension after the
       application of *this. However, if in new_modifications append_to_old()
       returns true (it then only consist of appending operation set up for the 
       old dimensio of this), then both need to have the same old variable 
       dimension and the appending operations of next_modification are carried
       out so that the stored reassignments in *this only affect to old
       coordinates of the appended information of next_modifcation, the new 
       coordinates added by *this are inserted as zeros in the appensions
       stored in next_modification and the new coordinates of next_modification 
       are appended at the end. 
       
    */
    virtual int incorporate(const OracleModification& next_modification)=0;

    /// returns a new object on the heap, that allows to incorporate this but starts off from a function whose input argument dimension is old_var_dim 
    virtual OracleModification* new_initial_oraclemodification(int old_var_dim) const =0;
       
    /// append append_dim further variables at the end of the argument vector (without specifying their effect, they may e.g. be ignored by the function)
    virtual int add_append_variables(int append_dim)=0;

    /// reorder and resize the variables as given by the first new_dim entries of map_to_old_indices; each former index may appear at most once in this list, so new_dim < get_new_vardim()
    virtual int add_reassign_variables(int new_dim,const int* map_to_old_indices)=0; 

    /// returns true if no modifications need to be executed
    virtual bool no_modification() const=0;

    /// if set to true, no deletions/reassignments may be present or specified in the future, only appensions are allowed 
    virtual int set_append_to_old(bool ) =0;
 
    /// returns true if this only contains appending operations 
    virtual bool append_to_old() const=0;

  };


  /** @brief If in Lagrangean relaxation primal solutions are in the form of a ConicBundle::DVector, then an approximate primal solution can be generated by supplying primal information of this form for each epsilon subgradient within ConicBundle::FunctionOracle::evaluate(). 

      In many applications, e.g. in Lagrangean relaxation, the convex
      minimization problem arises as the dual of a convex primal  
      maximization problem. In this case one is typically interested
      in obtaining a primal approximate solution in addition to the
      dual solution. Under reasonable conditions this is possible if the
      primal solutions that give rise to the subgradients are 
      aggregated along with the subgradients within the bundle algorithm.
      If the primal data can be represented as a DVector then the user
      has to supply in the oracle for each sugradient the corresponding primal  
      data in a PrimalDvector and the algorithm will do the rest. 
      Observe that a PrimalDVector can be used exactly in the same way as
      a DVector and that they are assignable among each other.

      The primal data has to be supplied within ConicBundle::FunctionOracle::Evaluate()
      and can be retrieved via the methods 
      ConicBundle::CBSolver::get_approximate_primal() and
      ConicBundle::CBSolver::get_center_primal()
 */
 
  class PrimalDVector: public PrimalData, public DVector
  {

  public:
    ///
    PrimalDVector(){}
    /// construct a primal double vector with n elements
    PrimalDVector(int n)                   :DVector((unsigned long)(n))   {}
    /// construct a primal double vector with n elements initialized to d
    PrimalDVector(int n ,double d)         :DVector((unsigned long)(n),d) {}
    /// copy constructor
    PrimalDVector(const PrimalDVector& pd) :PrimalData(),DVector(pd)  {}
    /// conversion from a DVector
    PrimalDVector(const DVector& pd)       :DVector(pd)  {}
    /// assignment of a DVector
    PrimalDVector& operator=(const DVector& pd) { DVector::operator= (pd); return *this; }
    
    /// produces a new PrimalDVector that is a copy of itself
    PrimalData* clone_primal_data() const{return new PrimalDVector(*this);}

    /// copy its information to *this
    int assign_primal_data(const PrimalData& it,double factor=1.)
    {
      const PrimalDVector* pd=dynamic_cast<const PrimalDVector*>(&it);
      assert(pd!=0);
      DVector::operator=(*pd);
      if (factor!=1)	
	for(unsigned int i=0;i<size();i++){
	  (*this)[i]*=factor;
	}
      return 0;
    }

    /// if it is a PrimalDVector of the same dimension, add factor*it to *this
    int aggregate_primal_data(const PrimalData& it,double factor=1.)
    {
      const PrimalDVector* pd=dynamic_cast<const PrimalDVector*>(&it);
      assert(pd!=0);
      if (size()!=pd->size()) return 1;
      for(unsigned int i=0;i<size();i++){
	(*this)[i]+=(*pd)[i]*factor;
      }
      return 0;
    }

    /// multiply/scale *this with a nonnegative myfactor
    virtual int scale_primal_data(double myfactor)
    { 
      if (myfactor!=1.){
	for(unsigned int i=0;i<size();i++){
	  (*this)[i]*=myfactor;
	}
      }
      return 0;
      
    }

  };


  /**@brief oracle interface (abstract class). For each of your functions, provide a derived class.  

     The oracle interface is used to describe and pass convex objective 
     functions to the ConicBundle::CBSolver. 
     The dimension of the argument vector of the function must be
     set in ConicBundle::CBSolver::init_problem() and the functions
     are then added to the solver by ConicBundle::CBSolver::add_function(). 

     If the sum of several such functions is to be minimized it
     is the task of the user to guarantee, that all dimensions
     match.

     If the function corresponds to Lagrangean relaxation of a 
     primal maximization problem one may want to generate a primal 
     approximate solution. In this case, return in the function
     FunctionOracle::evaluate() within the minorant
     the generating primal objects for each subgradient/minorant.
     If no primal objects are included, 
     there will be no primal aggregation.

     If primal aggregation is used then it is possible to implement
     a primal cutting plane framework. This requires the introduction
     of new (dual) variables in the design space of the function. 
     In this case a heap object MinorantExtender must be returend
     in the call of  FunctionOracle::get_minorant_extender(); this
     must be able to fill in the missing coordinates in existing
     minorants/subgradients maybe on basis of the associated primal
     data stored in the minorants. If this feature is not needed, 
     the function may be used as is and need not be reimplemented.
  */

  class FunctionOracle: public FunctionObject 
  {
  public:
        
    /** @brief Called by the solver. Has to Return function value and 
               at least one (epsilon) subgradient and, possibly for Lagrangean 
               relaxation, some primal data.

        The evaluation method is the main interface to the bundle
        solver.  The solver calls this method to obtain for the @a
        current_point (its dimension is set in
        ConicBundle::CBSolver::init_problem()) the @a objective_value
        and (epsilon) subgradient information.  In any call several
        epsilon subgradients may be returned in @a minorants along
        with their offset, but at least one has to be returend. Each
        subgradient and offset describes a linear minorant to the
        convex function and is used by the bundle method to form a
        cutting model of the convex function.
       
        In many applications, computing the function value is an iterative
        process that approaches the true function value from below. 
        On input the code offers a bound in @a objective_value for the 
        function value, above which it is certain that the code will reject 
        the current point. If in the iterative
        process a lower bound on the function value exceeds this bound, 
        then it is sufficient to return, instead of the true function value
        and a subgradient, the current lower bound and a vector so that
        together they describe a supporting hyperplane (i.e. a linear minorant) 
        to the function at this point. 

        If the function corresponds to Lagrangean relaxation of a 
        primal maximization problem one may want to generate a primal 
        approximate solution. For this purpose the minorants may also 
	hold information on the primal data. If at each call and for each 
	epsilon subgradient the corresponding generating primal object 
        (must be derived from ConicBundle::PrimalData, e.g., a 
        ConicBundle::PrimalDVector) is stored, then the code automatically 
        performs the aggregation corresponding to the aggregation
        of the subgradients on the primal data objects. The
        primal approximate solution is finally delivered by
        the methods ConicBundle::CBSolver::get_approximate_primal()
        or  ConicBundle::CBSolver::get_center_primal(). 

	All minorants passed to the solver must be objects allocated
        on the heap. The ownership of these objects is transferred 
        to the solver and the solver will destroy them eventually,
        so DO NOT delete them yourself!

        If no primal aggregation is desired, simply do not touch 
        @a primal_data or clear it.  

        @param[in] current_point (const double*)
           argument of the function (e.g. the Lagrange multipliers)

        @param[in] relprec (double)
           relative precision requirement for objective values 
           that may lead to descent steps (this precision is not
           required if it is already certain that the function 
           value will be too poor)

        @param[in,out]  objective_value (double&) 
         - on input: 
           value gives the threshold for a null step; you may stop, 
           if a cutting plane yields at least this; 
         - on output: 
           return an upper bound on the true function value within
           @a relprec *(abs(objval)+1.), if there is no linear minorant 
           cutting above the threshold specified in objective_value on input. 
           Otherwise the return value should be the max of cut_values.

        @param[out]  minorants (std::vector<Minorant*>)
           for returning linear minorants (subgradients and their offset
	   values). At least one must be returned and the one at index 0
	   must maximize the value at the current point over all further
	   ones given here. In particular its value at point should be
	   above the threshold or be within relprec *(abs(objval)+1.)
           of the true objective. Within the minorants primal data
           may be supplied if this should be aggregated along.

        @param[out]  primal_extender (PrimalExtender*&)
           if primal_data provided in minonrants of previous calls has 
	   now to be updated due to changes in the primal problem 
	   -- e.g., this may happen in column generation -- one may 
           return a pointer to PrimalExtender object on the heap. 
           This object will be used by ConicBundle to update all its 
           internally stored primal_data objects in its minorants
           by calling PrimalExtender::extend on each of these 
           (but not on those supplied by the new minorants).
	   Afterwards ConicBundle deletes primal_extender.
           If this is not needed, the variable holds 0.
  
        @return  
           -  0,  all correct
           -  !=0, failure. This does not necessarily terminate the bundle method.
               Termination is forced only if no new subgradient is returned. 
     */
    virtual
    int
    evaluate
    (
     const double* current_point,
     double relprec,
     double& objective_value,
     std::vector<Minorant*>&  minorants,
     PrimalExtender*& primal_extender 
     )= 0;


    /**@brief This routine need not be implemented unless variables
      (constraints in Lagrangean relaxation) are added or deleted on
      the fly

      The routine is only called by the solver if the variables indeed get
      modified. @a oracle_modification is then used to either transfer
      user supplied instructions to the oracle on how to modify itself
      or to inform the oracle about changes the solver was asked to perform
      on the variables. If available, the solver will also show the effect
      of these changes on the center point in @a new_center and @a old_center;
      if these are not available then they hold NULL. A user supplied
      @a oracle_modification will be checked for consistency with the
      actual changes in the variables and mismatches will cause failures.

      The remaining variables are output variables by which the oracle
      tells the solver which information has a chance to be preserved
      in view of these changes.  If e.g. the deletion of some nonzero
      variables invalidates the function value in the new center, the
      oracle has to set discard_objective_in_center=true.  If the
      entire model cannot be preserved (this includes the aggregates
      and the function values), the oracle needs to set
      discard_model=true; If only aggregate minorants cannot be
      preserved, the oracle needs to set
      discard_aggregates=true. Whenever new variables were added, the
      model can only be preserved if the remaining minorants (maybe
      without aggregates) can be extended for these new variables. In
      this case the oracle has to supply the appropriate
      MinorantExtender via @a minorant_extender. A given
      minorant_extender will only be applied if new variables were
      added and indices passed to it then refer to the situation
      *after* the changes (append and map-operation of
      OracleModification) were executed on them. If the operation
      fails for any of the minorants, the entire model will be
      discarded.
	      
      Return value 0 indicates that these actions allow to continue without
      errors, other return values result in an overall error on these changes.
    */

    virtual 
    int 
    apply_modification
    (
     const OracleModification& oracle_modification ,
     const double* new_center,
     const double* old_center,
     bool& discard_objective_in_center,
     bool& discard_model, 
     bool& discard_aggregates,
     MinorantExtender*& minorant_extender
     );
    /*{return 1;}*/
  
 
    /**@brief switch on/off some correctness checks on the oracle */
    virtual
    bool
    check_correctness() const
    {return true;}
      
       
  };
  

  /**@brief Serves for specifying parameters regarding the construction of cutting models.

  */

  class BundleParameters
  {
  protected:
    int max_model_size;       ///< maximum number of minorants to be selected for the cutting model (numbers<=1 for no limit, numbers >=2 impose a strict limit)
    int max_bundle_size;      ///< suggested maximum number of latest minorants stored for use in a model, for constructing variable metric information etc. (negative numbers give no preference; the size may be increased internally in case of confliciting requirements, eg. in n_model_size or by variable metric routines) 
    int update_rule;      ///< in case several update rules are available

  public:
    /// initialize to given values
    virtual int init(const BundleParameters& bp)
    {
      max_model_size=bp.max_model_size;
      max_bundle_size=bp.max_bundle_size;
      update_rule=bp.update_rule;
      return 0;
    }

    /// returns the value of the variable 
    virtual int get_max_model_size() const 
    { return max_model_size; }
    /// returns the value of the variable 
    virtual int get_max_bundle_size() const 
    { return max_bundle_size; }
    /// returns the value of the variable 
    virtual int get_update_rule() const 
    { return update_rule; }

    /// sets the value of the variable
    virtual int set_max_model_size(int mms)
    { max_model_size=mms; return 0; }
    /// sets the value of the variable
    virtual int set_max_bundle_size(int mbs)
    { max_bundle_size=mbs; return 0; }
    /// sets the value of the variable
    virtual int set_update_rule(int ur)
    { update_rule=ur; return 0; }

    /// often works well: small model of size 2 and some history in bundle size for use in scaling
    BundleParameters(const BundleParameters& bp)
    {BundleParameters::init(bp);}

    /// often works well for fast initial progress: small model of size 2 and some history in bundle size for use in scaling; default values give no preference at all
    BundleParameters(int modelsize=-1,int bundlesize=-1,int updaterule=-1)
    {
      max_model_size=modelsize;
      max_bundle_size=bundlesize;
      update_rule=updaterule;
    }

    /// virtual destructor
    virtual ~BundleParameters(); 

    /// return a new clone object of this on the heap (caller needs to delete the result)
    virtual BundleParameters* clone_BundleParameters() const
    { return new BundleParameters(*this); }

  };


  class MatrixCBSolver;
  
  /**@brief   Bundle method solver.
     
  Minimizes the sum of convex functions that are given via 
  ConicBundle::FunctionOracle interfaces, see 
  \ref cxxinterface "the text explaining the C++ interface" for
  a quick overview.  

  It provides special support for Lagrangean relaxation by generating
  primal approximate solutions if such information is provided in the
  function oracles. 

  Based on these primal approximations it is also possible to implement 
  cutting plane schemes. Routines for adding and deleting corresponding
  dual variables as well as a framework for extending subgradients in order
  not to loose the cutting model are available. 
  
  */	
  class CBSolver
  {
   private:
    MatrixCBSolver* solver;                     ///< pointer to internal solver

    CBSolver ( const CBSolver& );               ///< not available, blocked deliberately
    CBSolver& operator= ( const CBSolver& );    ///< not available, blocked deliberately
  public:    
    
    /// default constructor allows to set output level options from start (see also set_out())  
    CBSolver(std::ostream* out=0,int print_level=0);
    ///
    virtual ~CBSolver();
    
    //------------------------------------------------------------
    /**@name Initialization */
    //@{
       
    /** @brief Clears all data structures and problem information 
        but keeps ouptut settings and algorithmic parameter settings 
    */
    void clear();

    /** @brief Sets default values for algorithmic parameters that are not function specific (e.g., relative precision, weight and weight bounds for the augmentedproblem, etc.) 
    */
    void set_defaults();

    /** @brief Initializes the problem by setting up the design space 
        (the dimension and possible box constraints of the variables)
 
    Clears all data structures and sets the dimension @ m for a new problem.
    for solving   min_{y in R^m}  f_0(y) + f_1(y) + ...
    Box constraints may be specified for y. (The functions f_i must be added 
     by add_function()). 
  
     Lower and/or upper bounds must be speicified for all variables
     or for none of them. To specify no bounds at all, give Null
     pointers. Otherwise use ConicBundle::CB_minus_infinity for 
     unbounded below and ConicBundle::CB_plus_infinity for unbounded above.
     For NULL pointers, unbounded will be used as default for all
     variables. Specifying bounds selectively is also possible
     by set_lower_bound() or set_upper_bound().

     @param[in] dim  (int)
         the dimension of the argument/design space/the number of Lagrange multipliers
   
     @param[in] lbounds  (const DVector*)
         If NULL, all variables are considered unbounded below,
         otherwise lowerb[i] gives the minimum feasible value for variable y[i],
         use ConicBundle::CB_minus_infinity for unbounded below.

     @param[in] ubounds (const DVector*)
         If NULL, all variables are considered unbounded above,
         otherwise upperb[i] gives the maximum feasible value for variable y[i],
         use ConicBundle::CB_plus_infinity for unbounded above.

     @return 
        - 0 on success
        - != 0 otherwise

 
     */
    int init_problem(int dim, 
		     const DVector* lbounds=0,
		     const DVector* ubounds=0);
    
    /** @brief Adds a function, typically derived from ConicBundle::FunctionOracle; all functions added must have the same argument dimension.

     Besides the standard ConicBundle::FunctionOracle 
     the interface only accepts a few other prespecified derivations 
     of the class FunctionObject that come along with the CH_Matrix_Classes interface
     (e.g. for semidefinite and second order cones). Functions not derived from
     these will fail to be added and return a value !=0.

    @return 
      - 0 on success
      - != 0 otherwise
    */

    int
    add_function
    ( FunctionObject& function );
    
    /**@brief Sets lower bound for variable i,
	use ConicBundle::CB_minus_infinity for unbounded from below.

       The algorithm may have to adapt the center point aftwards.
       In this case the old function values will be marked as outdated and 
       will be recomputed at the next call to e.g. solve().
 
    @return 
      - 0 on success
      - != 0 otherwise
    */
    int set_lower_bound(int i, double lb);

    /**@brief Sets upper bound for variable i,
	use ConicBundle::CB_plus_infinity for unbounded from below.

       The algorithm may have to adapt the center point aftwards.
       In this case the old function values will be marked as outdated and 
       will be recomputed at the next call to e.g. solve().

     @return 
      - 0 on success
      - != 0 otherwise
    */
    int set_upper_bound(int i, double ub); 


    /** @brief Append new variables (always in last postions in this order).

       If 0 is feasible for the new coordinates then this is selected as
       starting value for the new coordinates; otherwise, the number
       closest to zero is used. If all new coordinates can be set to zero 
       then it is assumed that for an existing center point the 
       function values need not be recomputed (this is e.g. the case in 
       Lagrangean relaxation; if this is not correct call 
       reinit_function_model() below). 
       Otherwise the old function values will be marked as outdated and 
       will be recomputed at the next call to e.g. solve().

     @attention Be sure to update your objective functions so that they can 
       handle the new variables before you call this and any further ConicBundle 
       routines that require function evaluations. 
       Also, these operations may lead to inavailability of certain 
       other data such as subgradients and primal approximations.
 
     @param[in] n_append  (int)
       number of variables to append (always in last position in the same order)
  
     @param[in] lbounds  (DVector*)
        If NULL, all appended variables are considered unbounded below,
        otherwise lowerb[i] gives the minimum feasible value for variable y[i],
        use ConicBundle::CB_minus_infinity for unbounded below.

     @param[in] ubounds (DVector*)
        If NULL, all appended variables are considered unbounded above,
        otherwise upperb[i] gives the maximum feasible value for variable y[i],
        use ConicBundle::CB_plus_infinity for unbounded above.

     @return 
        - 0 on success
        - != 0 otherwise


    */
    int append_variables(int n_append, 
			 const DVector* lbounds=0,
			 const DVector* ubounds=0);

    /** @brief Deletes variables corresponding to the specified indices. 

        The indices of the remaining variables are reassigned so that they
        are consecutive again, the routine returns in @a map_to_old 
        a vector giving for each new index of these remaining variables 
        the old coordinate.

        If all of the deleted variables are zero, function values are assumed 
        to remain correct (if this is not so, call reinit_function_model() below)
        Otherwise the old function values will be marked as outdated and 
        will be recomputed at the next call to e.g. solve().

     @attention Be sure to update your objective functions so that they can 
        handle the new variables before you call any further ConicBundle 
        routines that require function evaluations.
        Also, these operations may lead to inavailability of certain 
        other data such as subgradients and primal approximations.

     @param[in] delete_indices  (const IVector&)
        the entries delete_indices[i] specify the indices of the variables 
        to be deleted
  
     @param[out] map_to_old  (IVector&) 
        after the call, element map_to_old[i] gives the old index (before the call)
        of the variable that now has index position i.
  
     @return 
        - 0 on success
        - != 0 otherwise

    */
    int delete_variables(const IVector& delete_indices,IVector& map_to_old);

    /** @brief Reassigns variables to new index positions by mapping to position @a i 
        the variable that previously had index @a assign_new_from_old[i].

        Old variables, that are not mapped to any position will be deleted.
        It is allowed to generate several copies of old variables. 

        If all of the deleted variables as well as new multiple copies are zero, 
        function values are assumed to remain correct (if this is not so, 
        call reinit_function_model() below).
        Otherwise the old function values will be marked as outdated and 
        will be recomputed at the next call to e.g. solve().

     @attention Be sure to update your objective functions so that they can 
        handle the new variables before you call any further ConicBundle 
        routines that require function evaluations.
        Also, these operations may lead to inavailability of certain 
        other data such as subgradients and primal approximations.

     @param[in] assign_new_from_old  (const IVector&)
        entry assign_new_from_old[i] specifies
        the old index of the variable, that has to be copied to index position i.
    
     @return 
        - 0 on success
        - != 0 otherwise

    */
    int reassign_variables(const IVector& assign_new_from_old);

    //@}

    //------------------------------------------------------------
    /**@name Basic Algorithmic Routines and Parameters */
    //@{
    
    /** @brief solves or does a prescribed number of iterations

	Bundle methods solve a problem by a sequence of so called
        descent steps that actually bring progress by moving from the
        current "center point" to a new center with better objective.
        A descent step may consist of several function evaluations
        (null steps), that lead to no immediate progress but mainly
	improve a cutting model of the objective function close to
        the current center point.  A minimizer to the model is
        accepted as descent step if the function value at this point
        satisfies a sufficient decrease criterion in comparison to the
        decrease predicted by the model. Having found a descent step,
        the next center is automatically shifted to this successful
        candidate.  Termination criteria may stop the process of
        seeking for a descent step, in which case the current center
        is kept and the routine termination_code() returns the
        termination code.

        Restarting, after each descent step, the bundle method from scratch 
        with the new center as starting point does not endanger convergence.
        Therefore, a descent step is the smallest unit, after which
        user interaction can take place safely. To allow this there
	is a flag stop_at_descent_steps that will cause the code to
	return after the next descent step.

	If you know what your are doing, you may also use the input
        parameter maxsteps to force the algorithm to return after
        at most maxsteps null steps. Calling solve() again
        without any intermediate problem configurations will then
        simply continue the process where it stopped and convergence
        is save. During null steps one may not decrease the weight
        or delete nonzero variables of the center or the current candidate!
     
        In a Lagrangean relaxation cutting plane approach one may want
        to separate and enlarge the dimension after a certain number
        of null steps. In this case the code will try to preserve the model,
        given appropriate subgradient extension routines have been 
        provided. If the model cannot be extended, it has to be
        discarded (if subgradient extension is not successful
        this is done automatically), and the algorithm will be restarted
        from the current center point. 
	
 
      @param[in] maxsteps (int)
          if maxsteps>0 the code returns after at most so many null steps 

      @param[in] stop_at_descent_steps (int)
          if true the code also returns whenever a descent step occured 

      @return 
        - 0 on success
        - != 0 otherwise
    
    */
    int 
    solve(int maxsteps=0,bool stop_at_descent_steps=false);
    
    /** @brief Returns the termination code of the bundle algorithm for the latest descent step

      For resetting all counters relevant for termination see clear_fail_counts() .
 
      @return
      -  0  :    Not terminated. 
             (Continue with the next solve())
      -  1  :    Relative precision criterion satisfied. (See set_term_relprec())
      -  2  :    Timelimit exceeded. 
             (Currently the C interface does not offer a timelimit.)
      -  4  :    Maximum number of function reevaluations exceeded. 
             (Indicates that there is a problem with one of the function 
             oracles that seems to deliver no valid upper bounds on the true 
             function value for descent steps)
      -  8  :    Maximum number of quadratic subproblem failures exceeded.
             (Indicates that the numerical limits of the inner quadratic 
             programming solver are reached, no further progress expected) 
      - 16  :    maximum number of model evaluation failures exceeded
             (Indicates that the numerical limits of the setup of the 
             subproblem are reached, no further progress expected) 
      - 32  :    maximum number of failures to increase the augmented model value exceeded
             (Indicates that the numerical limits  of the interplay between 
             subproblem and quadratic programming solver are reached, 
             no further progress expected) 
       - 64  :   maximum number of oracle calls (function evaluations) exceeded,
                 see set_eval_limit()
       - 128  :   maximum number of oracle failures exceeded. 
             This refers to function evaluations that terminate with insufficient 
	     precision but still provide a new approximate subgradient. A failure typically
             indicates numerical difficulties with the precision requirements.  
             (Currently the interface does not allow to manipulate the limit, it is set to 10)

    */
    int
    termination_code() 
      const; 
    
    /** @brief Outputs a text version of termination code, see termination_code(). 

      @return 
        - 0 on success
        - != 0 otherwise

    */
    std::ostream&
    print_termination_code(std::ostream& out); 
    
    /** @brief Returns the objective value resulting from last descent 
        step (initially undefined). If no problem modification routines 
        were called since then, it is the objective value at the point 
        returned by get_center().
    */
    double 
    get_objval()
      const;
    
    /** @brief Returns the next center point that was produced by the latest call  
	to solve() (in some problem modification routines the
	center point may be updated immediately, in others the center point 
	will be corrected automatically directly before starting 
	the next descent step and its values may be infeasible till then).

      @return 
        - 0 on success
        - != 0 otherwise
    */
    int 
    get_center
    ( DVector& center )
      const;
    
    
    /** @brief Returns Euclidean norm of the latest aggregate subgradient.
     */
    double 
    get_sgnorm()
      const; 

    /** @brief Returns the latest aggregate subgradient.

    @return 
      - 0 on success
      - != 0 otherwise

    */
    int 
    get_subgradient
    ( DVector& subgradient )
      const;
 
    /** @brief Returns the objective value computed in the last step of solve(), 
        independent of whether this was a descent step or a null step (initially undefined). 

        If no problem modification routines were called since then, it is the 
        objective value at the point returned by get_candidate(). If this 
        last evaluation led to a descent step, then it is the same value as
        in get_objval().
    */
    virtual double 
    get_candidate_value()
      const;
    
    /** @brief Returns the last point, the "candidate", at which the function 
        was evaluated in solve(). 

        If this evaluation lead to a descent step, it is the same point as 
	in get_center().

      @return 
        - 0 on success
        - != 0 otherwise
    */
    virtual int 
    get_candidate
    ( DVector& candidate )
      const;
    
    
    //@}
    
    //------------------------------------------------------------
    /**@name Advanced Algorithmic Routines and Parameters */
    //@{

    /** @brief Sets the relative precision requirements for successful termination
               (default 1e-5).

     @param[in] term_relprec (double)
       The algorithm stops with termination code 1, if predicted progress for 
       the next step is less than term_relprec times
       absolute function value plus one. 

    @return 
      - 0 on success
      - != 0 otherwise

    */
    int 
    set_term_relprec
    ( const double term_relprec );
    
    /** @brief Set the starting point/center that will be used in the 
        next call to  solve(). Each call
        to this routine causes an immediate evaluation of all oracles. 

     @return  
        - 0 on success 
        - != 0 otherwise
    */
    int
    set_new_center_point
    ( const DVector& center_point );

    
    /** @brief Returns the return value of the latest evaluation call
        to this @a function.
    */
    int 
    get_function_status
    ( const FunctionObject& function )
      const;

    /** @brief Returns the multipliers for the box constraints on the design variables;
        in Lagrangean relaxation they may be interpreted as primal slacks
	for inequality constraints.
     @return 
        - 0 on success
        - != 0 otherwise
    */
    int 
    get_approximate_slacks
    ( DVector& center )
      const;
    
    /** @brief returns the current approximate primal solution corresponding
        to the aggregate subgradient of the specified @a function. 
        
        PrimalData solutions must have been supplied in all previous 
        calls to evaluate; In this case it returns the current approximate 
        primal solution aggregated alongside with the aggregate subgradient.  
        A primal solution may not be available after addition of constraints, 
        if extension of the aggregate subgradient to the new coordinates failed.
	If no primal data is availalbe, the function returns NULL.

     @return 
        - pointer to the primal data of the aggregate of this function object
        - 0 if no primal is available
    */
    const PrimalData*
    get_approximate_primal
    ( const FunctionObject& function)
      const; 
    
    /** @brief Returns the primal solution corresponding to the best epsilon
	subgradient returned in the evaluation of the specified @a function 
        at the current center point. If no primal data is availalbe, 
	the function returns NULL.
  
     @return 
        - pointer to the primal data of the minorant returned on evaluation 
	  of this function object at the current center
        - 0 if no primal is available
    */
    const PrimalData*
    get_center_primal
    ( const FunctionObject& function)
      const; 


    /** @brief Returns the primal solution corresponding to the best epsilon
	subgradient returned in the evaluation of the specified @a function 
        at the point get_candidate. If no primal data is availalbe, 
	the function returns NULL.
       
     @return 
        - pointer to the primal data of the minorant returned on evaluation 
	  of this function object at the current candidate
        - 0 if no primal is available
    */
    const PrimalData*
    get_candidate_primal
    (const FunctionObject& function)
      const; 


    /** @brief Sets the maximum number of subgradients used in forming the
        cutting model of the specified @a function  

        Quite often a very small model, e.g., 2, yields very fast iterations
        and good progress in time (sometimes at the cost of more evaluations).
        By limited numerical experience, a significant reduction in the number of 
        evaluations can  only be expected if the bundle is large enough to 
        wrap the function rather tightly. Quite frequently, unfortunately,
        this entails that solving the quadratic subproblems
        is more expensive than function evaluation. 

        The meaning of this routine may differ from standard for 
	predefined special functions with special bundle types.

     @param[in] function (const FunctionObject&)
        the function added in add_function()

     @param[in] max_modelsize (int)
         maximum number of subgradients to be used in forming the cutting model 

     @return 
        - 0 on success
        - != 0 otherwise

    */
    int 
    set_max_modelsize
    (  const FunctionObject& function, int max_modelsize );

    /** @brief Sets the maximum number of subgradients stored for use in
	forming the model or determining scaling information, it must be as
	least as large as max_modelsize (and is increased to this if not)

        The meaning of this routine may differ from standard for 
	predefined special functions with special bundle types.

     @param[in] function (const FunctionObject&)
        the function added in add_function()

     @param[in] max_bundlesize (int)
       maximum number of new epsilon subgradients to be used in bundle updates 

     @return 
       - 0 on success
       - != 0 otherwise

    */
    int 
    set_max_bundlesize
    (  const FunctionObject& function, int max_bundlesize );

    /**@brief Sets the maximum bundlesize and the maximum number of new subgradients
        added in a bundle update of the cutting model for the specified @a function.
        The meaning of this routine may differ from standard for 
	predefined special functions with special bundle types.

     @param[in] function (const FunctionObject&)
        the function added in add_function()

     @param[in] params (const BundleParameters&)
       some update parameters for the cutting model, see e.g. ConicBundle::BundleParameters

     @return 
       - 0 on success
       - != 0 otherwise

    */
    int 
    set_bundle_parameters
    (  const FunctionObject& function,
       const BundleParameters& params );

    /** @brief Retrieves current bundle parameters (not the actual size in use!)
        as set for the cutting model of the specified @a function.

	This may differ for predefined special 
	functions with derived BundleParameter classes. 

     @param[in] function (const FunctionObject&)
        the function added in add_function()

     @return 
       - 0 if no such function or such data is available
       - otherwise a pointer to the BundleParameters 

    */
    const BundleParameters* 
    get_bundle_parameters
    ( const FunctionObject& function )
      const;


    /** @brief Clears cutting model, subgradients and stored function values 
       for the specified @a function
 
        This has to be called whenever the specified function was modified
        so that the old subgradients and/or primal generators are no longer
        valid. 

      @param[in] function (const FunctionObject&)
        the function added in add_function()

      @return 
        - 0 on success
        - != 0 otherwise

    */
    int reinit_function_model
    ( const FunctionObject& function );
    
    /** @brief Asks @a function to call @a primal_extender for each of its primal objects (see
     also FunctionOracle::evaluate() )

     If the function is the Lagrangian dual of a primal problem and primal_data 
     returned previous calls to the oracle has now to be updated due to changes 
     in the primal problem -- e.g., this may happen in column generation -- the 
     call causes updates of all internally stored primal_data objects by calling 
     PrimalExtender::extend on each of these.

      @param[in] function (const FunctionObject&)
        the function added in add_function()

      @param[in] primal_extender (PrimalExtender&)
        the object holding the extension function for primal_data

      @return 
        - 0 on success
        - 1 if for this function it is not possible to use a primal_extender
        - 2 if the primal_extender would be applicable but there is no primal_data
    */
    int call_primal_extender
    (const FunctionObject& function,
     PrimalExtender& primal_extender);

    /** @brief Returns the current weight for the quadratic term in the augmented subproblem
    (may be interpreted as 1./step_size or 1./trustregion-radius).
    */
    double 
    get_last_weight() const; 

    /** @brief Sets the  weight (>0) to be used in the quadratic term 
        of the next augmented subproblem 
        (may be interpreted as 1./step_size or 1./trustregion-radius).

        Independent of whether the weight violates current min- and max-bounds 
        set in set_min_weight() and set_max_weight(), the next model will 
        be computed for this value. Thereafter, however, it will be updated as 
        usual; in particular, it may be truncated by min and max bounds 
        immediately after the first subproblem. 

        In order to guarantee a constant weight (e.g. 1 is frequently a reasonable 
        choice if the automatic default heuristic performs poorly), set the min and max
         bounds to the same value, too.

      @param[in] weight (double)

      @return 
        - 0 on success
        - != 0 otherwise
    */
    int 
    set_next_weight
    ( const double weight );
    
    /** @brief Sets a lower bound on the  weight for the quadratic term of the 
        augmented subproblem.

        Nonpositive values indicate no bound. 
        The new value shows its effect only at first dynamic change of 
        the weight.

     @param[in] min_weight (double)

     @return 
        - 0 on success
        - != 0 otherwise

    */
    int 
    set_min_weight
    ( const double min_weight ); 
    
    /** @brief Sets an upper bound on the  weight for the quadratic term of the 
        augmented subproblem.

        Nonpositive values indicate no bound. 
        The new value shows its effect only at first dynamic change of 
        the weight.

      @param[in] max_weight (double)

      @return 
        - 0 on success
        - != 0 otherwise
    */
    int 
    set_max_weight
    ( const double max_weight );
    
    
    /** @brief Use a scaling heuristic or switch off scaling alltogether.
       (the scaling heuristic resets the quadratic term to some diagonal matrix,
       switching it off resets the diagonal term to the identity)

     @param[in] do_scaling (int) 
        - 0 switch off the scaling heuristic
        - 1 use a diagonal scaling heuristic
        - 2 use a diagonal scaling heuristic combined with one for the bounds
	
      @return 
        - 0 on success
        - != 0 otherwise
    */
    virtual int
    set_variable_metric
    ( int do_scaling );

    /** @brief If set to true (the default is false),  
        some variables will be fixed automatically to the center 
        value if their bounds are strongly active (i.e., the
        corresponding multipliers are big). 

        The coordinates to be fixed are redetermined in each
        call following a descent step or a change of the function. 
        An indicator vector of the variables fixed in the last call 
        can be obtained via the routine get_fixed_active_bounds().

        Setting this value to true might improve the performance
        of the algorithm in some instances but there is no 
        convergence theory. It might be particularly helpful
        within Lagrangian relaxation if a primal cutting plane
        approach is used and non-tight inequalities should be
        eliminated quickly (fixing then indicates large primal 
        slack values). 

      @param[in] allow_fixing (bool)

     */ 
    virtual void
    set_active_bounds_fixing
    ( bool allow_fixing );

    /** @brief clears all fail counts on numerical function oder model failures,
	may be useful if this caused premature termination. 

    */
    virtual void 
    clear_fail_counts
    (void);
    
    /** @brief Sets an upper bound on the number of calls to the oracle (use negative numbers for no limit). 

        If this number is reached, the algorithm will terminate 
	independently of whether the last step was a descent or 
	a null step. A negative number will be interepreted as 
	no limit.
 
      @param[in] eval_limit (int)
    */
    virtual void 
    set_eval_limit
    (int eval_limit);

    /** @brief Set an upper bound on the number of inner updates for the
        cutting model with primal slacks within one null step (use negative numbers for no limit). 

        A negative number will be interepreted as no limit, i.e., 
        the updates will be done till a certain precision of the 
        cutting model is achieved.
 
      @param[in] update_limit (int)
      
    */
    virtual void 
    set_inner_update_limit
    (int update_limit);
 
    //@}
    
    //------------------------------------------------------------
    /**@name Look Up Basic Paramaters (dimension, number of functions, ...)*/
    //@{

    /** @brief Returns the current dimension of the design space/argument
           or -1 if no dimension is set. 
    */
    int get_dim();

    /** @brief Returns the current number of functions in the problem. 
    */
    int get_n_functions();

    /** @brief Returns the indicator vector of variables temporarily fixed to 
	the center value due to significantly positive multipliers
        for the box constraints.

        Such a fixing indicates that the corresponding 
        variables would like to stay at their bounds.
        There will be nonzero entries only if set_active_bounds_fixing() 
	is switched on.
     */
    int get_fixed_active_bounds(IVector& fixed_active_bounds) const;

    //@}

    //------------------------------------------------------------
    /**@name Output */
    //@{

    /** @brief Specifies the output level (out==NULL: no output at all, 
           out!=NULL and level=0: errors and warnings, 
           level>0 increasingly detailed information)

     @param[in] out  (ostream*) 
       direct all output to (*out). If out==NULL, there will be no output at all.

     @param[in] print_level (int)

     Output levels for print_level: 
      -  0 ... no output except for errors and warnings 
      -  1 ... line summary after each descent step  
      - >1 ... undocumented and increasingly detailed log information. 
             These higher levels should only be used if requested 
             for debugging purposes. 

      Example for level 1:

\verbatim
00:00:00.00 endit  1   1   1   563.    563.  39041.188  39043.162
00:00:00.00 endit  2   2   2   563.    559.  38488.165  38490.200
00:00:00.00 endit  3   3   3   56.3    555.  33014.533  33211.856
00:00:00.00 endit  4   4   4   5.63    517. -14306.459  2738.0343
00:00:00.00 endit  5   5   5   4.04    148. -2692.1131  2.2150883
00:00:00.00 endit  6   6   6   4.01    1.29  1.7908952  2.0000581
00:00:00.00 endit  7   7   7   3.95  0.0213  1.9999387  2.0000000
00:00:00.00 _endit  8   8   8   3.95 2.94e-05  2.0000000  2.0000000

Column 1      2     3   4   5    6       7       8          9
\endverbatim
      - Column 1: computation time in hh:mm:ss.dd,
      - Column 2: "endit" is convenient for grep and stands for "end of iteration".
         Iterations with termination_code()!=0 are marked with "_endit".
      - Column 3: number of descent steps
      - Column 4: number of descent and null steps. Up to initialization calls
         and reevaluations, this is the number of evaluation calls 
         to the function oracles from within the bundle method. 
         In the example all calls led to descent steps.
      - Column 5: number of innermost iterations. It differs from column 5 only in the
          case of variables with bounds in which case it gives the number of updates
          of the multipliers for the bounds (or primal slacks in Lagrangean 
          relaxation). Exceedingly high numbers in this column indicate that
          some variables are constantly at their bounds and it might be 
          possible to improve convergence by deleting them (i.e. set them
          as constants to this bound and remove the variable). 
      - Column 6: the weight of the quadratic term in the augmented problem.
      - Column 7: the norm of the aggregate subgradient. If it is small,
          say below 0.1, then mostly this is good indication that the
          objective value is close to optimal. 
      - Column 8: the value of the cutting model in the last candidate point. It
          is always a lower bound on the true function value in this point
      - Column 9: the objective value in the latest point that led to a descent 
          step, i.e., the point returend by get_center(). Whenever 
          termination_code() returns 0 this is also the objective
          value of the latest evaluation call to the function oracles and
          the value in the center point of the next iteration. 
     */
    void 
    set_out
    (std::ostream* out=0,int print_level=1);

    //@}
    
  };

  //@}
}



#endif



