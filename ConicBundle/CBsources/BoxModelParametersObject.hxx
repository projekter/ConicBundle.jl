/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/BoxModelParametersObject.hxx
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



#ifndef CONICBUNDLE_BOXMODELPARAMETERSOBJECT_HXX
#define CONICBUNDLE_BOXMODELPARAMETERSOBJECT_HXX

/**  @file BoxModelParametersObject.hxx
    @brief Header declaring the class ConicBundle::BoxModelParametersObject
    @version 1.0
    @date 2020-03-28
    @author Christoph Helmberg
*/

#include "SumBlockModel.hxx"
#include "BoxData.hxx"

namespace ConicBundle {
  /** @ingroup InternalBundleModel

  */
  //@{

  /** @brief abstract interface for BoxModel for the model selection routine
    select_model()

    Because BoxModel implements a model for an AffineMatrixTransformation
    given support function over the box, all subgradients are actually
    points inside the box (convex combinations of its extreme points).

    The possibilities supported by BoxModel allow to specify a (possibly
    empty) subset of coordinates (@a box_coords) that are then bounded
    by their respective intervalls, together with one point (@a
    box_complvalues) contained in the box of the remaining coordinates,
    plus an arbitrary number of standard subgradients (@a nnc_model)
    with non negative cone coefficients (@a nnc_coeff). The trace
    constraint makes sure that the resulting point is (a nonnegative
    multiple @a box_scaleval of) a convex combination of the box
    representation @a box_model and the @a nnc_model.

    If the coordinate subset comprises all coordinates -- the full or
    exact model of the box -- no other points are needed and no other
    information may be specified. In this case @a box_model and @a
    box_coeff have the same size as @a box_coords and @a box_coeff
    directly gives the coordinates of the point.

    If the coordinates are a nonempty proper subset of all coordinates,
    the size of @a box_model and @a box_coeff is the size of @a
    box_coords plus one, where the last element of @a box_model holds
    the information of @a box_complvalues and the last element of @a
    box_coeff gives the scaling of the entire box, possibly as convex
    combination factor to the elements of @a nnc_model. In particular,
    assume gamma=box_coeff(box_coeff.rowdim()-1)>0 and let 0<=i<@a
    box_coeff.rowdim()-1 be an index into @a box_coords with orignial
    coordinate j=@a box_coords(i), then lb(j)<= @a box_coeff(i)/gamma <=
    ub(j) is the value of coordinate j.

    It does not matter, which point is put into @a box_complvalues
    as long as the entire model together with the points in @a nnc_model
    covers the aggregate and the new candidate subgradient.
  */

  class BoxModelParametersObject : public virtual CBout, public BundleParameters {
  protected:
    //int max_model_size;       // maximum number of minorants to be selected for the cutting model (numbers<=1 for no limit, numbers >=2 impose a strict limit)
    //int max_bundle_size;      // suggested maximum number of latest minorants stored for use in a model, for constructing variable metric information etc. (negative numbers give no preference; the size may be increased internally in case of confliciting requirements, eg. in n_model_size or by variable metric routines) 
    //int update_rule;      // in case several update rules are available

    //bool enforce_separate_model; // when forming combined models for a sum of functions, enforce a separate model for this function

  public:
    /// initialize BundleParameters to given values, if bp is a SumModleParametersObject alos set n_local_models to its values, otherwise leave it unchanged
    int init(const BundleParameters& bp) {
      BundleParameters::init(bp);
      const BoxModelParametersObject* mpo = dynamic_cast<const BoxModelParametersObject*>(&bp);
      if (mpo) {
        set_cbout(mpo, 0);
      }
      return 0;
    }

    /// default constructor
    BoxModelParametersObject(const CBout* cb = 0, int cbinc = -1) :
      CBout(cb, cbinc), BundleParameters() {
    }

    /// copy constructor for BundleParameters
    BoxModelParametersObject(const BundleParameters& bp, const CBout* cb = 0, int cbinc = -1) :
      CBout(cb, cbinc), BundleParameters(bp) {
    }

    ///copy constructor
    BoxModelParametersObject(const BoxModelParametersObject& sms) :
      CBout(sms), BundleParameters(sms) {
    }

    ///virtual destructor, implemented in BoxModelParameters.cxx
    virtual ~BoxModelParametersObject();

    /** @brief BoxModel calls this for selecting the next coordinates
         for a specialized polyhedral model with a box part and an
         nnc part for aggregates, see the general explanation of the
         class.

         There is little experience on how to do this.
         The current routine is minimalistic and simply uses
         the weighted number of switches in each coordinate
         between lower and upper bounds to select coordinates.

         To allow maybe better choices in future implementations
         the arguments try to pass all potentially relevant
         information items.

         On putput the coefficient values of the new model must be
         feasible and have to generate the same aggregate (the aggregate
         is maintained), and the new candidate minorant must be in the
         feasible set.

      @param[in,out] box_model
         the boxmodel holds the minorants describing the BoxBlock part of
         the model for selected coordinates and, unless exact or empty,
         in the last position the complement coordinates of a feasible
         point (e.g. aggr_boxvec)

      @param[in,out] box_coeff
         - on input coefficients of BoxBlock determined in last
           BundleMethod::eval_augmodel giving rise (together with nnc_coeff)
           to the aggrgate,
         - on output they match the new model and still give rise to
           the same aggregate

      @param[in,out] box_indicators
         indicators for activity of box minorants, indicators may but need
         not be maintained

      @param[in,out] box_coords
         the coordinates selected to have their respective interval range
         in the model

      @param[in,out] box_complvalues
         a point with 0 in the box_coords and feasible
         coordinate values in the complement (if not empty)

      @param[in,out] nnc_model
         if box_model is empty, this spans at least the aggregate (if available)
         and the candidate (always); if box_model is not empty but not
         the entire box, nnc_model typically holds one of the candidate or
         the aggregate; if box_model is the entire box, nnc_model is empty

      @param[in,out] nnc_coeff
         - on input coefficients of BoxBlock determined in last
           BundleMethod::eval_augmodel giving rise (together with box_coeff)
           to the aggrgate,
         - on output they match the new model and still give rise to
           the same aggregate

      @param[in,out] nnc_indicators
         indicators for activity of minorants, indicators may but need
         not be maintained

      @param[in] coord_switching
         keeps track of which coordinates where changing the most in the past
         by forming a weighted average in BoxModel::eval_function()

      @param[in] minorants
          the vector of MinorantPointer gives additional minorants
    collected over time (some may be duplicates also of those in
          nnc_model)

     @param[in] cand_minorant
         the (eps)sugradient linear minorant returned
         by BoxModel::eval_function for the candidate (without function factor)

     @param[in] cand_boxvec
         the maximizer over the box for the current candidate

     @param[in] aggr_boxvec
         the primal aggregate vector in the box (without function_factor)
         giving rise to the aggregate; not initialized if zerodimensional

     @param[in] aggr_scaleval
         0<= aggr_scalevale <= function_factor, ==function_factor if
         function_task==Objective_Function; the aggregate with
         function_factor is aggr_boxvec*aggr_scaleval;

      @param[in] oracle
          gives access to lower and upper bounds of the box

      @param[in] modification_id
          the identifier of the current function version to be used in generating
          specialized minorants corresponding to the coordinate vectors

      @param[in] function_task
          see FunctionTask

      @param[in] function_factor
           interpreted according to function_task and the coefficients sum up to at most this value

      @param[in] model_update
          informs about whether cand_y is the result of a null_step or descent_step or a new set up.

      @param[in] center_id
          the identifier of the center point

      @param[in] center_y
          the center point

      @param[in] cand_id
          the identifier of the candidate point

      @param[in] cand_y
          the candidate (mostly differnt from the center), close to it the model should be good

      @param[in] model_maxviol
          a minorant violated by this would have caused a null step

      @param[in] H
          the variable metric used in the proximal term (function_factor is already removed in this)

      @return
       - 0 on success
       - 1 on failure
    */
    virtual int select_model(MinorantBundle& box_model,
      CH_Matrix_Classes::Matrix& box_coeff,
      CH_Matrix_Classes::Matrix& box_indicators,
      CH_Matrix_Classes::Indexmatrix& box_coords,
      CH_Matrix_Classes::Matrix& box_complvalues,
      MinorantBundle& nnc_model,
      CH_Matrix_Classes::Matrix& nnc_coeff,
      CH_Matrix_Classes::Matrix& nnc_indicators,
      const CH_Matrix_Classes::Matrix& coord_switching,
      const MinorantBundle& minorants,
      const MinorantPointer cand_minorant,
      const PrimalMatrix cand_boxvec,
      const PrimalMatrix aggr_boxvec,
      CH_Matrix_Classes::Real aggr_scaleval,
      BoxOracle* oracle,
      CH_Matrix_Classes::Integer modification_id,
      FunctionTask function_task,
      CH_Matrix_Classes::Real function_factor,
      BundleModel::ModelUpdate model_update,
      CH_Matrix_Classes::Integer center_id,
      const CH_Matrix_Classes::Matrix& center_y,
      CH_Matrix_Classes::Integer cand_id,
      const CH_Matrix_Classes::Matrix& cand_y,
      CH_Matrix_Classes::Real model_maxviol,
      BundleProxObject& H) = 0;


  };


  //@}

}

#endif

