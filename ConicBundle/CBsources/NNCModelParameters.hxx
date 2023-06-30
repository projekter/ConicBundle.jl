/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/NNCModelParameters.hxx
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



#ifndef CONICBUNDLE_NNCMODELPARAMETERS_HXX
#define CONICBUNDLE_NNCMODELPARAMETERS_HXX


/**  @file NNCModelParameters.hxx
    @brief Header declaring the class ConicBundle::NNCModelParameters
    @version 1.0
    @date 2020-03-28
    @author Christoph Helmberg
*/

#include "NNCModelParametersObject.hxx"

namespace ConicBundle {

  /** @ingroup InternalBundleModel

  */
  //@{


  /** @brief default model selection routines for NNCModel

      It is always made sure that the aggregate and the minorant with maximum
      value are contained in the model. If max_model_size==2 this is all there
      is. If max_model_size<2, there is no limit on the model size, otherwise
      max_model_size is strictly observed.

      update_rule in {0,1} are implemented by select_active_model() and
      are minimal fuss methods for keeping those non aggregate minorants
      that indicators show as active and aggregating the rest in either
      all steps (udate_rule 0) or just in descent steps or when the size
      limit is reached (update_rule 1).

      uldate_rule 2 is an older heuristic that is still present
      if someone happens to miss it, but it won't ever be
      documented.

  */
  class NNCModelParameters : public NNCModelParametersObject {
  private:
    /** @brief this implements update_rule 0 and 1: Add maximizing
       minorant, keep positive indicators if possible. Udate rule 1
       grows up to maximum size during null steps before aggregating
       again.

        - update_rule==0 (grow=false): add the maximizing minorant, keep
          all non aggregated minorants that are selected by indicators
          and aggregate the rest; if these are too many, include them by
          largest coefficient values. Note: here aggregates are always
          combined to a single aggregate, even when active.

        - update_rule==1: identical to update_rule==0 during descent steps
          and whenever the maximum size is exceeded with the next minorant
    (put grow=false). Otherwise (put grow=true) the minorant with
    maximum value is added without eliminating other minorants.
    Thus the model grows during null steps, maybe collapsing again
    when the size limit is reached.
    */
    int select_active_model(MinorantBundle& bundle,
      CH_Matrix_Classes::Matrix& coefficients,
      CH_Matrix_Classes::Matrix& indicators,
      CH_Matrix_Classes::Integer cand_id,
      const CH_Matrix_Classes::Matrix& cand_y,
      const MinorantBundle& minorants,
      bool grow);


    /// heuristic for updating the model grown over some time (update_rule==2); do not expect any documentation
    int select_model_heuristic(MinorantBundle& bundle,
      CH_Matrix_Classes::Matrix& coefficients,
      CH_Matrix_Classes::Matrix& activity_indicators,
      CH_Matrix_Classes::Integer cand_id,
      const CH_Matrix_Classes::Matrix& cand_y,
      const MinorantBundle& cand_minorants,
      const MinorantPointer& aggregate,
      FunctionTask function_task,
      CH_Matrix_Classes::Real function_factor);


  public:
    /// default constructor with the possibility to set the output
    NNCModelParameters(const CBout* cb = 0, int incr = -1) :
      NNCModelParametersObject(cb, incr) {
    }

    /// constructor for size parameters
    NNCModelParameters(int modelsize, int bundlesize = 10, int updaterule = 0, const CBout* cb = 0, int incr = -1) :
      NNCModelParametersObject(BundleParameters(modelsize, bundlesize, updaterule), cb, incr) {
    }

    /// copy constructor for BundleParameters
    NNCModelParameters(const BundleParameters& bp, const CBout* cb = 0, int incr = -1) :
      NNCModelParametersObject(bp, cb, incr) {
    }

    ///copy constructor
    NNCModelParameters(const NNCModelParameters& sms) :
      CBout(sms), NNCModelParametersObject(sms) {
    }

    /// destructor
    virtual ~NNCModelParameters();

    /// clone
    BundleParameters* clone_BundleParameters() const {
      return new NNCModelParameters(*this);
    }


    /** @brief NNCModel calls this for selecting the next minorants for a polyhedral model

      @param[in,out] model
          contains the model of the previous subproblem or a minimal starting choice

      @param[in,out]  coefficients
          contains the coefficients resulting from the last bundel subproblem (on input and on output they should generate the aggregate, if the aggregate is valid)

      @param[in,out]  activity_indicators
          contains activity_indicators resulting from the last bundel subproblem (currently 1 if considered active and 0 otherwise). Those indicated active in the bundle on input will also be indicated active

      @param[in] aggregate
          if valid, it holds the current aggregate. It arises as nonnegative combinations of the model by the coefficients (so it inculde the current function_factor). This property of model and coefficients should be maintained on output

      @param[in] center_minorant
          if valid, it holds the current aggregate. It arises as nonnegative combinations of the model by the coefficients, and that should be maintained on output

      @param[in] cand_minorants
          holds the vector of minorants returned by evaluation in the candidate point

      @param[in] old_minorants
          the vector of MinorantPointer gives additional minorants
    collected over time (some may be duplicates, some are
    most certainly already contained in the bundle on input)

      @param[in] oracle
          gives access to the evaluation oracle

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
    virtual int select_model(MinorantBundle& model,
      CH_Matrix_Classes::Matrix& coefficients,
      CH_Matrix_Classes::Matrix& activity_indicators,
      const MinorantPointer& aggregate,
      const MinorantPointer& center_minorant,
      const MinorantBundle& cand_minorants,
      const MinorantBundle& old_minorants,
      MatrixFunctionOracle* oracle,
      FunctionTask function_task,
      CH_Matrix_Classes::Real function_factor,
      BundleModel::ModelUpdate model_update,
      CH_Matrix_Classes::Integer center_id,
      const CH_Matrix_Classes::Matrix& center_y,
      CH_Matrix_Classes::Integer cand_id,
      const CH_Matrix_Classes::Matrix& cand_y,
      CH_Matrix_Classes::Real model_maxviol,
      BundleProxObject& H);


  };

  //@}

}

#endif

