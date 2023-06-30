/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/QPModelDataObject.hxx
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


#ifndef CONICBUNDLE_QPMODELDATAOBJECT_HXX
#define CONICBUNDLE_QPMODELDATAOBJECT_HXX


/**  @file QPModelDataObject.hxx
    @brief Header declaring the class ConicBundle::QPModelDataObject, ConicBundle::QPSumModelDataObject, ConicBundle::QPConeModelDataObject, ConicBundle::QPModelDataPointer
    @version 1.0
    @date 2019-05-18
    @author Christoph Helmberg
*/

#include "MinorantPointer.hxx"
#include "AffineFunctionTransformation.hxx"

namespace ConicBundle {




  /** @defgroup InternalQPModelData Interface for qp solver dependent model descriptions

     @brief The classes define the interface for obtaining a solver dependent description
     of the bundle subproblem. The solver has to provide (or be) a QPModelDataPointer,
     which generates on request (issued in ConicBundle::BundleModel::start_augmodel() )
     a possibly model dependent QPModelDataObject that the model then fills with its
     information. Derived from this base class QPModelDataObject there is
     QPSumModelDataObject used by ConciBundle::SumModel and QPConeModelDataObject
     that is general enough to fit the models of ConicBundle::NNCModel,
     ConicBundle::SOCModel, ConicBundle::PSCModel, ConicBundle::BoxModel, and
     ConicBundle::AFTModel. The cutting model description typically consists of the
     collected subgradient information and the support set describing the feasible
     combinations of the subgradient information. A potentially present sequence of
     AffineFunctionTransformation instances modifies this subgradient information
     and the routine QPModelDataPointer::apply_aft() serves this purpose.
     Several classes derived from BundleProxObject also use this interface to
     compute the quadratic costs for the internal ConicBundle::UQPSolver.
  */
  //@{

    /** @brief abstract interface declaring the uniform outer appearance of general cutting models produced by the various models for use in generating the model data so that QP Solvers can provide their own tuned data blocks
     */

  class QPModelDataObject : public virtual CBout {
  public:
    /// clear
    virtual void clear() = 0;
    /// default constructor
    QPModelDataObject(CBout* cb = 0, int cbinc = -1) :CBout(cb, cbinc) {
    }
    /// virtual destructor
    virtual ~QPModelDataObject();

    ///gives reading access to a constant offset minorant
    virtual const MinorantPointer& get_constant_minorant() const = 0;

    ///gives reading access to the (possibly transformed) bundle minorants of the cutting model
    virtual const MinorantBundle& get_bundle() const = 0;

    ///gives access to a constant offset minorant
    virtual MinorantPointer& get_constant_minorant() = 0;

    ///gives access to the bundle minorants of the cutting model so that an AffineFunctionsTransformation and similar can be applied
    virtual MinorantBundle& get_bundle() = 0;

    /** @brief applies the AffineFunctionTransformation to
        constant_minorant and bundle, where (if given) only the
        global_indices of the transformed subgradients are required
        which need the local_indices only.

        If precomputed is given, it may contain some or contains
        afterwards a map from original minorant to transformed minorant;
        returns 0 on success.

        The given aft and indices will be stored
        and may be used during the solution of the QP subproblem,
        therefore the objects pointed to have to remain valid throughout
        their use in the QP.

        The AFT has to pushed onto the model data after the subgradient
        data has been initialized. Multiple AFTs have to pushed in the
        correct sequence so that the first AFT pushed is the first to
        be applied to the subgradient data, etc.
    */
    virtual int push_aft(const AffineFunctionTransformation* inaft,
      const CH_Matrix_Classes::Indexmatrix* global_indices,
      const CH_Matrix_Classes::Indexmatrix* local_indices,
      std::map<MinorantPointer, MinorantPointer>* precomputed = 0) = 0;

    /// undo the top push_aft this 
    virtual int pop_aft() = 0;

    /// this routine is added here mainly to allow for comperative testing; for the current value of x it adds the aggregate minorant (offset and gradient) to the arguments. Caution: the bundle used might have been modified by AffineFuncitonTransformations and need not be the same as initialized with.
    virtual int add_modelx_aggregate(CH_Matrix_Classes::Real& offset,
      CH_Matrix_Classes::Matrix& gradient) = 0;



  };



  /** @brief abstract interface extension of QPModelDataObject for use in SumModel where several blocks are concatenated for a common bundle
   */
  class QPSumModelDataObject : public virtual QPModelDataObject {
  public:
    /// default constructor
    QPSumModelDataObject(CBout* cb = 0, int cbinc = -1) :QPModelDataObject(cb, cbinc) {
    }

    /// virtual destructor
    virtual ~QPSumModelDataObject();

    /// do nothing if inblock==0, otherwise add the constant_minorant of inblock and append the bundle of inblock; returns 0 on success
    virtual int append(QPModelDataObject* inblock) = 0;
  };

  /** @brief abstract base class for passing additional oracle information
      to the QP

     This class is only intended for passing pointers to information
     via QPConeModelDataObject::init(), but not for storing information.

     If passed to the QPSolver, the information is assumed to exist
     throught the time the QPSolver works on solving the bundle
     subproblem for the current model.  The Solver is allowed to use any
     of the objects pointed to in different routines maybe at the same
     time, so in providing the information routines special care must be
     taken to avoid side effects.

  */

  class QPModelOracleDataObject {
  public:
    ///virtual destructor
    virtual ~QPModelOracleDataObject();

  };

  /** @brief abstract interface extension of QPModelDataObject to allow uniform generation of tuned quadratic solver model blocks in model forming routines for specializd oracles
     */
  class QPConeModelDataObject : public virtual QPModelDataObject {
  public:
    /// default constructor
    QPConeModelDataObject(CBout* cb = 0, int cbinc = -1) :QPModelDataObject(cb, cbinc) {
    }

    /// virtual destructor
    virtual ~QPConeModelDataObject();

    /** @brief set the model data of the current local model, only positive trace values are allowed

        @param[in] constant_minorant
           serves as a common offset for all minorants, may be empty

        @param[in] bundle
           the minorant information (not necessarily all feasible
     minorants by themselves); the sequence must match the
     sequence of the variables in the nonempty cones and box
     bounds in the argument list. The linear combination of
     this minorant information by feasible
     conic variable settings observing the trace constraint
     must give rise to a valid minorant (-> the "aggregate")
     of the underlying convex function.

        @param[in] nncdim
           the dimension of the nonnegative cone; together with
     the trace constraint restricting the value of the sum
     of the single variables this typically models convex
     combinations of minorants.

        @param[in] socdim
           a vector giving the dimensions of the second order
     cone variables to be used in this sequence.

        @param[in] pscdim
           a vector giving the order of the positive semidefinite
     matrix variables. In this sequence each matrix variable
     X of order n has to have a matching sequence of n*(n+1)/2
     items in bundle, the orderinig within being as in
     CH_Matrix_Classes::svec().

        @param[in] box_lb
           the vector specifies the lower bounds (with matching
     upper bounds in box_ub) of a corresponding sequence
     of box constrained variables

        @param[in] box_ub
           the vector specifies the upper bounds (with matching
     lower bounds in box_lb) of a corresponding sequence
     of box constrained variables

        @param[in] rhs_trace
           a positive value specifiying the right hand side of
     the trace constraint (if present), i.e., the sum
     of the NNC-variables, the 0-indexed entries of the
     SOC vectors, the diagonal entries of the PSC matrix
     variables, and one common mulitplicative factor variable
     for jointly scaling the values of box_lb and box_ub
           (this implicit variable scales the size of the box if
     there is one).

        @param[in] ft
           if the value is ObjectiveFunction, the trace constraint
           is used as equality constraint and the rhs_value serves
     as a function factor. Otherwise the trace constraint
     is used as in inequality constraint with an additional
     slack variable; this models a variable function factor
     between 0 and rhs_trace, so that the function is used
     as penalty function with penalty value rhs_trace.

        @param[in] oracle_data
           Some QP solvers may be able to exploit direct access to
     the oracle or some other additional data for improving
           their performance. If this is the case, this allows
     to pass such information. The object pointed to has
     to stay valid throughout the solution process of this
     QP subproblem and will not be destroyed here.

        @param[in] scale_box
           this variable is only used if the box dimension is not zero
     and all other cones have dimension zero.
     In this case, scale_box==true allows the box to be scaled
     whereas scale_box==false fixes the common scaling variable
     of box_lb and box_ub to rhs_trace, in this case
     there is no trace constraint.

        @return
           0 if all is ok, all other values refer to some error

     */
    virtual int init(const MinorantPointer& constant_minorant,
      const MinorantBundle& bundle,
      CH_Matrix_Classes::Integer nncdim,
      const CH_Matrix_Classes::Indexmatrix& socdim,
      const CH_Matrix_Classes::Indexmatrix& pscdim,
      const CH_Matrix_Classes::Matrix& box_lb,
      const CH_Matrix_Classes::Matrix& box_ub,
      CH_Matrix_Classes::Real rhs_trace,
      FunctionTask ft,
      QPModelOracleDataObject* oracle_data = 0,
      bool scale_box = true) = 0;

    ///  only positive trace values are allowed, this 
    virtual int adjust_trace(CH_Matrix_Classes::Real next_rhs_trace) = 0;

    /// after the bundle subpolbem is sovled returns the sum of the trace values of x of the local model
    virtual CH_Matrix_Classes::Real evaluate_trace() const = 0;

    /// after the bundle subproblem is solved, this retrieves the local nonnegative cone  solution vector; if nncx_activity is set, the return values between zero and one indicate the guess on the coefficients activity level 
    virtual int get_nncx(CH_Matrix_Classes::Matrix& nncx,
      CH_Matrix_Classes::Matrix* nncx_activity = 0,
      bool cautious = false) = 0;

    /// after the bundle subproblem is solved, this retrieves the local soc solution vector i; if socx_activity is set, the return value between zero and one indicates the guess on the activity level  
    virtual int get_socx(CH_Matrix_Classes::Integer i,
      CH_Matrix_Classes::Matrix& socx,
      CH_Matrix_Classes::Real* socx_activity = 0,
      bool cautious = false) = 0;

    /// after the bundle subproblem is solved, this retrieves the local psc solution i as eigenvalue factorization; the growthrate is <x,z>/<x^⁻,z^->, where x^- and z^- refer to the previous and x, z to the current. 
    virtual int get_pscx(CH_Matrix_Classes::Integer i,
      CH_Matrix_Classes::Matrix& pscx_eigs,
      CH_Matrix_Classes::Matrix& pscx_vecs,
      CH_Matrix_Classes::Real& pscx_growthrate,
      CH_Matrix_Classes::Matrix& pscx_primalgrowth,
      CH_Matrix_Classes::Matrix& pscx_dualgrowth) = 0;

    /// after the bundle subproblem is solved, this retrieves the local box solution vector; if boxx_activity is set, the values between zero and one indicate the guess on the coefficients activity level 
    virtual int get_boxx(CH_Matrix_Classes::Matrix& x,
      CH_Matrix_Classes::Matrix* x_activity = 0,
      bool cautious = false) = 0;

    /// after the bundle subproblem is solved, this returns the dual variable to the trace constraint which yields an upper bound to the model value for the local model without the constant minorant; if prec is given, the termination precision used in its computation (here typically the current value of the barrier parameter) is returned here
    virtual CH_Matrix_Classes::Real tracedual(CH_Matrix_Classes::Real* prec = 0) const = 0;

  };


  class BundleModel;

  /** @brief abstract interface for automatically generating and retrieving the qp blocks with the model data that are needed for specific implementations of the qp solver.

      In BundleModel::start_augmodel() the model may already have a previous
      QPModelDataObject. It first test compatibility of this object by using
      set_model_data(); if this fails, it ask for generating an returning
      a new one by calling either generate_summodel_data() or generate_conemodel_data();
      The model then enters the cutting model information into the
      available or freshly generated object und uses set_model_data to
      store the pointer to final model information in the QPModelDataPointer.
      This information can then be retrieved by teh caller via get_model_data_ptr().

   */

  class QPModelDataPointer : public virtual CBout {
  public:
    /// default constructor
    QPModelDataPointer(CBout* cb = 0, int cbinc = -1) :CBout(cb, cbinc) {
    }

    /// virtual destructor
    virtual ~QPModelDataPointer();

    /// set the pointer to NULL
    virtual void clear_model_data_ptr() = 0;

    ///store the pointer to the object if it matches the required type for the QP solver, otherwise return a nonzero value as error; this is used in the models to return the local qp model data
    virtual int set_model_data(QPModelDataObject*) = 0;

    ///returns a new QPSumModelDataObject, that has to be deleted by the caller. The argument is optional and allows to potentially generate different blocks for different derived BundleModel objects; this is used in SumModel to collect the models of the various oracles that are summed over 
    virtual QPSumModelDataObject* generate_summodel_data(BundleModel* bmp = 0) = 0;

    ///returns a new QPConeModelDataObject suitable for the default conic BundleModel implementations; it has to be deleted by the caller. The argument is optional and allows to potentially generate specialized objects for special BundleModel objects 
    virtual QPConeModelDataObject* generate_conemodel_data(BundleModel* bmp = 0) = 0;

    /// returns the pointer value
    virtual QPModelDataObject* get_model_data_ptr() const = 0;
  };


  //@}

}

#endif

