/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/AffineFunctionTransformation.hxx
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



#ifndef CONICBUNDLE_AFFINEFUNCTIONTRANSFORMATION_HXX
#define CONICBUNDLE_AFFINEFUNCTIONTRANSFORMATION_HXX

/**  @file AffineFunctionTransformation.hxx
    @brief Header declaring the classes ConicBundle::AffineFunctionTransformation
    @version 1.0
    @date 2014-08-04
    @author Christoph Helmberg
*/

#include "CBout.hxx"
#include "CBSolver.hxx"
#include "MinorantPointer.hxx"
#include "AFTModification.hxx"
#include "GroundsetModification.hxx"


namespace ConicBundle {

  /** @ingroup cxxmatrixinterface
  */
  //@{


  //-----------------------------------------------------------------------------
  //                                AffineFunctionTransformation
  //-----------------------------------------------------------------------------

  /** @brief transform a function f(z) to
      fun_coeff*f(arg_offset+arg_trafo*y)+linear_cost*y+fun_offset (scales the
      value, substitutes argument z=c+Ay, and adds an affine term b'y+d)

      The AffineFunctionTransformation (AFT for short) allows to use an implemented function
      oracle for \f$f\colon\mathbf{R}^n\to\mathbf{R},\quad z\mapsto f(z)\f$ as a
      function oracle of the function

      \f[ g:\mathbf{R}^m\to\mathbf{R},\quad y\mapsto g(y)=\sigma f(c+Ay)+b^\top y+\delta \f]

      where
      - \f$\sigma\ge 0\f$ (Real @a fun_coeff, default value 1) scales the function value (0 is allowed to cancel the function for cases where removing it might cause difficulties)
      - \f$c+Ay\f$ specifies an affine transformation of the argument with
           + \f$c\in\mathbf{R}^n\f$ (Matrix pointed to by @a arg_offset,
             NULL pointer represents \f$c=0\f$)
     + \f$A\in\mathbf{R}^{n\times m}\f$ (Sparsemat pointed to by @a arg_trafo, NULL pointer represents \f$A=I_n\f$, so indeed the IDENTITY and NOT the zero matrix)
      - \f$b\in\mathbf{R}^m\f$ (Matrix pointed to by @a linear_cost, NULL pointer represents \f$b=0\f$) adds a linear cost term \f$b^\top y\f$
      - \f$\delta\in\mathbf{R}\f$ (Real @a fun_offset, default value=0) adds a constant offset

      If (some of) the pointers to matrices are specified, ownership of the
      matrix objects pointed to will be taken over by the
      AFT and the objects will be deleted by this class
      at the end of their use.

      Note, if @a arg_offset and @a arg_trafo are both zero, then the argument
      does not change at all. If all three matrix pointers are zero, there is no
      information on the dimension of the transformation, so it will be
      determined by the arguments supplied.

      Using value 0. for @a fun_coeff is allowed and removes the influence of
      the function completely, leaving only the affine function in y specified by
      @a fun_offset and @a linear_cost.

      Internally, the SumBlockModel realization of the AFT is achieved by AFTModel.

  */

  class AffineFunctionTransformation : public CBout, public FunctionObject {
  private:
    CH_Matrix_Classes::Real fun_coeff; ///< function value is multiplied by this 
    CH_Matrix_Classes::Real fun_offset;  ///< constant added in the end 
    CH_Matrix_Classes::Matrix* linear_cost;   ///< NULL means no linear cost 
    CH_Matrix_Classes::Matrix* arg_offset;    ///< NULL means no offset
    CH_Matrix_Classes::Sparsemat* arg_trafo;  ///< NULL means identity (!!!)

    MinorantPointer constant_minorant; ///<  corresponds to fun_offset + <*linear_cost,.>

    mutable GroundsetModification local_gsmdf; ///< this describes the effect of the last modification on the image space (if there was one at all)

    bool model_calls_delete; ///< tells the model whether it should delete this at the end of its use or just leave it alone (i.e., this AFT is then owned by someone else)

  public:
    /// sets the parameters of the transformation. The ownership of objects pointed to is passed to *this (they will be deleted here). If *this is entered into an AFTModel, model_calls_delete==true tells the AFTModel to delete this AffineFunctionTransformation at the end.  
    virtual int init(CH_Matrix_Classes::Real fun_coeff = 1.,
      CH_Matrix_Classes::Real fun_offset = 0.,
      CH_Matrix_Classes::Matrix* linear_cost = 0,
      CH_Matrix_Classes::Matrix* arg_offset = 0,
      CH_Matrix_Classes::Sparsemat* arg_trafo = 0,
      bool model_calls_delete = true);

    /// calls init()
    AffineFunctionTransformation(CH_Matrix_Classes::Real in_fun_coeff = 1.,
      CH_Matrix_Classes::Real in_fun_offset = 0.,
      CH_Matrix_Classes::Matrix* in_linear_cost = 0,
      CH_Matrix_Classes::Matrix* in_arg_offset = 0,
      CH_Matrix_Classes::Sparsemat* in_arg_trafo = 0,
      bool in_model_calls_delete = true,
      CBout* cbo = 0, int incr = -1) :
      CBout(cbo, incr), linear_cost(0), arg_offset(0), arg_trafo(0) {
      init(in_fun_coeff, in_fun_offset, in_linear_cost, in_arg_offset, in_arg_trafo, in_model_calls_delete);
    }

    /// deletes @a linear_cost, @a arg_offset and @a arg_trafo
    virtual ~AffineFunctionTransformation();

    /// retruns true if the model has to delete this
    bool get_model_calls_delete() {
      return model_calls_delete;
    }

    /// set to true if the model has to delete this, to false if it is destructed elsewhere
    void set_model_calls_delete(bool mcd) {
      model_calls_delete = mcd;
    }

    /// returns true if not the identity
    bool argument_changes() const {
      if ((arg_offset != 0) || (arg_trafo != 0)) return true; return false;
    }

    /** @brief returns true if arg_trafo influences at most two thirds of the
  entries of local_argument
     */
    bool sparse_argument_changes() const {
      return ((arg_trafo != 0) && (3 * arg_trafo->get_rowinfo().rowdim() <= arg_trafo->rowdim()));
    }

    /** @brief returns true if the transformation maps each index (in col_ind
  if !=0) onto at most one index and vice versa (out of row_ind if !=0)
     */
    bool scaled_index_subset(const CH_Matrix_Classes::Indexmatrix* col_ind,
      const CH_Matrix_Classes::Indexmatrix* row_ind) const;

    /** @brief returns false if index is mapped to more than one index,
  otherwise true with mapped_index==-1 if mapped to zero, else >=0 and
  @a coeff gives the coefficient
     */
    bool scaled_index(CH_Matrix_Classes::Integer& mapped_index,
      CH_Matrix_Classes::Real& coeff,
      CH_Matrix_Classes::Integer index) const;

    /// returns the factor for the function
    CH_Matrix_Classes::Real get_fun_coeff() const {
      return fun_coeff;
    }

    /// allows to set the factor for the function
    CH_Matrix_Classes::Real& set_fun_coeff() {
      return fun_coeff;
    }

    /// returns the constant offset for the funciton
    CH_Matrix_Classes::Real get_fun_offset() const {
      return fun_offset;
    }

    /// allows to set the constant offset for the funciton
    CH_Matrix_Classes::Real& set_fun_offset() {
      return fun_offset;
    }

    /// returns the pointer to the linear term added to the funciton
    const CH_Matrix_Classes::Matrix* get_linear_cost() const {
      return linear_cost;
    }

    /// returns the pointer to the constant offset added to the argument (not neeeded in the code)
    const CH_Matrix_Classes::Matrix* get_arg_offset() const {
      return arg_offset;
    }

    /// returns the pointer to the linear transformation of the argument (not neeeded in the code)
    const CH_Matrix_Classes::Sparsemat* get_arg_trafo() const {
      return arg_trafo;
    }

    /// returns the value of the linear cost coefficient for @a i>=0 and for i==-1 the constant offset    
    CH_Matrix_Classes::Real get_linear_cost(CH_Matrix_Classes::Integer i) const {
      if (i < 0) return fun_offset;
      if (linear_cost) return (*linear_cost)(i);
      return 0.;
    }

    /// return the constant linear minorant corresponding to fun_offset+<*linear_cost,.>
    const MinorantPointer& get_constant_minorant() const {
      return constant_minorant;
    }

    /// returns the dimension of the input argument or -1 if it is unknown  
    CH_Matrix_Classes::Integer from_dim() const {
      if (arg_trafo)
        return arg_trafo->coldim();
      if (linear_cost)
        return linear_cost->rowdim();
      if (arg_offset) //arg_trafo==NULL is the identity, so this is ok
        return arg_offset->rowdim();
      return -1; //no information on the dimension
    }

    /// returns the dimension of the output argument or -1 if it is unknown  
    CH_Matrix_Classes::Integer to_dim() const {
      if (arg_trafo)
        return arg_trafo->rowdim();
      if (arg_offset)
        return arg_offset->rowdim();
      if (linear_cost) //arg_trafo==NULL is the identity, so this is ok 
        return linear_cost->rowdim();
      return -1; //no information on the dimension
    }

    /// only copies the elements that are effected by the image of arg_trafo  
    int copy_traforows(CH_Matrix_Classes::Matrix& copy_to,
      const CH_Matrix_Classes::Matrix& copy_from) const;

    /** @brief computes @a transformed_offset and, if @a this is not the identity,@a transformed_y and returns either @a input_y (id) or @a transformed_y

  On input @a transformed_y is supposed to be of dimension zero, unless
  it was already the result of a previous transformation for exactly
  this AFT with this data (no intermediate modifications). This is
  important, because if @a transformed_y has the correct size, only the
  data changed by the transformation is overwritten and the constant
  parts are assumed to be already initialized. On output it is again of
  dimension 0 if the transformation is the identity.

        The @a transformed_offset = fun_offset+ip(linear_cost,@a input_y) is
        the value needed in objective_value() together with the result of the
        evaluation of the function for @a transformed_y in order to compute
        the transformed objective value.
    */
    const CH_Matrix_Classes::Matrix&
      transform_argument(CH_Matrix_Classes::Matrix& transformed_y,
        CH_Matrix_Classes::Real& transformed_offset,
        const CH_Matrix_Classes::Matrix& input_y) const;


    /** @brief given the modification aftmdf or if 0, gsmdf, compute the
        transformed argument that would arise after this modification as in
  transform_argument()

  On input @a transformed_y is supposed to be of dimension zero.
  The matrix returned is transformed_y unless the modification
  preserves the identity transformation. In this case, the
        returned matrix is input_y.
    */
    const CH_Matrix_Classes::Matrix&
      modified_transform_argument(CH_Matrix_Classes::Matrix& transformed_y,
        const CH_Matrix_Classes::Matrix& input_y,
        const AFTModification* aftmdf,
        const GroundsetModification& gsmdf) const;


    /** @brief  if @a offset is the value computed in transform_argument and
       @a function_value results form an evaluation in the respective point,
       the routine returns the objective value obtained by the affine transformation
     */
    CH_Matrix_Classes::Real objective_value(CH_Matrix_Classes::Real offset,
      CH_Matrix_Classes::Real function_value) const {
      return offset + fun_coeff * function_value;
    }

    /** @brief transforms the in linear minorant (scaled by alpha) and
        initializes or adds it to the out linear minorant

        if out_minorant.empty()==true, the out_minorant is initialized,
        otherwise the information is added.

        the constant_minorant of the transformation is only added if
        add_trafo_minorant is set to true.

  If given, provided_row_indices and needed_column_indices must be as
        specified in qp_cost_indices, i.e. of in_minorant only the
        provided_row_indices (all if NULL) will be used to compute only the
        needed_column_indices (as rows, all if NULL) of out_minorant.
     */
    int transform_minorant(MinorantPointer& out_minorant,
      const MinorantPointer& in_minorant,
      CH_Matrix_Classes::Real alpha = 1.,
      bool add_trafo_minorant = false,
      const CH_Matrix_Classes::Indexmatrix* provided_row_indices = 0,
      const CH_Matrix_Classes::Indexmatrix* needed_column_indices = 0) const;

    /** @brief transforms several linear minorants (scaled by alpha) and
        initializes or adds them to the out linear minorants

        for out_minorant[i].empty()==true the out_minorant is initialized,
        otherwise the information is added.
     */
    int transform_minorants(MinorantBundle& out_minorants,
      const MinorantBundle& in_minorants,
      CH_Matrix_Classes::Real alpha = 1.) const;



    /** @brief if the algorithm only requires the indices given in @a
       needed_col_indices (NULL means all indices) then it suffices to supply
       the @a provide_row_indices as @a indices in transform_minorant (a 0x0
       return matrix again means all indices, while a 0x1 matrix means no
       indices). Any input or output indices must be in strictly increasing
       order.
     */
    int qp_cost_indices(CH_Matrix_Classes::Indexmatrix& provide_row_indices,
      const CH_Matrix_Classes::Indexmatrix* needed_col_indices) const;

    /** @brief if AFTModel::add_diagonal_scaling() is called with @a indices
  specified by @a col_indices, then AFTmodel has to provide a
  diagonal scaling matrix of dimension to_dim() as @a in_diagscale
  in add_diagonal_scaling() with the entries in the output vector
  @a row_indices computed correctly.
     */
    int scaling_indices(CH_Matrix_Classes::Indexmatrix& row_indices,
      const CH_Matrix_Classes::Indexmatrix& col_indices) const;

    /** @brief transform the scaling matrix in_diagscale of the untransformed
  function model and add it as described in
  BundleMethod::add_diagonal_scaling.

        If indices are given, they must be sorted in striclty increasing order.
    */
    int add_diagonal_scaling(CH_Matrix_Classes::Matrix& diagscale,
      const CH_Matrix_Classes::Indexmatrix* indices,
      CH_Matrix_Classes::Real alpha,
      const CH_Matrix_Classes::Matrix& in_diagscale) const;

    /** @brief returns information on the changes in the ground set of the
        transformed arguments and checks whether after the modification the
  nonzero image of the transformed minorants will not differ from before;
        returns NULL on errors.

  If arg_trafo==NULL (act as identity) and if aftmdf!=NULL adds
  any explicit matrix parts or has modifications not consistent with
  maintaining the identity, arg_trafo is virtually first set to an explicit
  identity before executing the modification. If arg_trafo==NULL and
  aftmdf==NULL, the transformations of gsmdf are applied to linear_cost
  and arg_offset.

  If the input @a gsmdf and the AFT modifications in @a aftmdf result in
  a different GroundsetModification of the transformed space, this is
  computed into the member @a local_gsmdf and the returned pointer
  points to this @a local_gsmdf. If all transformations are passed on
  exactly as in @a gsmdf (because the trafo still acts as the identity),
  the returned pointer points to @a gsmdf.

        If the changes affect the transformation at all in its effect
  on the nonzero image of transformed minorants,
        @a minorant_trafo_differs will be set to true otherwise to false.

     */
    const GroundsetModification*
      analyze_modification(bool& minorant_trafo_differs,
        const AFTModification* aftmdf,
        const GroundsetModification& gsmdf) const;

    /** @brief if arg_trafo==NULL (act as identity) and if aftmdf!=NULL adds
  any explicit matrix parts or has modifications not consistent with
  maintaining the identity, arg_trafo is first set to an explicit
  identity before executing the modification. If arg_trafo==NULL and
  aftmdf==NULL, the transformations of gsmdf are applied to linear_cost
  and arg_offset.
     */
    int
      apply_modification(const AFTModification* aftmdf,
        const GroundsetModification& gsmdf);


    /// for testing purposes this outputs the data in mfile readable form
    virtual std::ostream& output_aft_data(std::ostream& out) const;


  };



  //@} 

}

#endif

