/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/SumBundleParametersObject.hxx
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



#ifndef CONICBUNDLE_SUMBUNDLEPARAMETERSOBJECT_HXX
#define CONICBUNDLE_SUMBUNDLEPARAMETERSOBJECT_HXX


/**  @file SumBundleParametersObject.hxx
    @brief Header declaring the class ConicBundle::SumBundleParametersObject
    @version 1.0
    @date 2020-03-28
    @author Christoph Helmberg
*/


#include "BundleModel.hxx"
#include "SumBundle.hxx"

namespace ConicBundle {

  /** @ingroup InternalBundleModel

  */
  //@{


  /** @brief abstract interface for SumBundleHandler for the model selection and variable metric seletcion routines

      In designing model selection routines for the SumBundleHandler
      note that model selection on basis of the previous model does not
      make too much sense here. SumBundle and SumBundleHandler are
      mainly designed for frequent switching of contributing functions,
      so the model is likely formed for a different function each
      time. It is therefore rather pointless to rely on the coefficients
      of the previous model. Simple geometric properties might provide
      better guidelines.

      In contrast to the typical *ModelParametersObject this interface
      not only provides the abstract model selection routine
      select_model() but also gives access to a potentially available
      variable metric selection routine
      get_variable_metric_selection(). The reason is, that there is no
      explicit access to the various SumBundleHandler instances other
      than via the recusive use of
      SumBlockModel::sumbundle_contribution() that must pass on all
      settings to all SumBundleHandler instacnes. Therefore all
      relevant parameters are passed on within this one type of object.


  */
  class SumBundleParametersObject : public virtual CBout, public BundleParameters {
  private:
    /// see set_acceptable_mode() for explanations
    SumBundle::Mode acceptable_mode;

  protected:

    /// points to a VariableMetricSelection object owned by *this, may be NULL
    VariableMetricSelection* vm_selection;


  public:

    /// copies the values of BundleParameters, leaving vm_selection unchanged, unless bp is a SumBundleParametersObject, then it clones the VariableMetricSelection pointed to by bp vm_selection (or also sets it to 0). 
    virtual int init(const BundleParameters& bp);

    /// default constructor with the possibility to set the output
    SumBundleParametersObject(const CBout* cb = 0, int incr = -1) :
      CBout(cb, incr), BundleParameters(), acceptable_mode(SumBundle::inactive), vm_selection(0) {
    }

    /// copy constructor for BundleParameters
    SumBundleParametersObject(const BundleParameters& bp, const CBout* cb = 0, int incr = -1) :
      CBout(cb, incr), BundleParameters(), acceptable_mode(SumBundle::inactive), vm_selection(0) {
      init(bp);
    }

    /// copy constructor
    SumBundleParametersObject(const SumBundleParametersObject& bp) :
      CBout(&bp, 0), BundleParameters(), acceptable_mode(SumBundle::inactive), vm_selection(0) {
      init(bp);
    }

    ///
    virtual ~SumBundleParametersObject();

    /// see set_acceptable_mode() for possible settings and their meaning
    virtual SumBundle::Mode get_acceptable_mode() const {
      return acceptable_mode;
    }

    /** @brief specifies the sumbundle modes the model accepts

        The possible settings have the following meaning
        - SumBundle::unavailable ... do not use SumBundle here and remove an existing bundlehandler
        - SumBundle::inactive ... play along with any of the choices that are suggested when called
        - SumBundle::root ... always use a SumBundle with *this as root, do not accept another role suggested by a parent; start a SumBundle here if not yet available
        - SumBundle::child ... only accept the role SumBundle::child and SumBundle::inactive if suggested by a parent but switch to local model when suggested to be SumBundle::root
     */
    virtual int set_acceptable_mode(SumBundle::Mode mode) {
      acceptable_mode = mode; return 0;
    }

    /// delete old selector and set a new one (0 is allowed resulting in no local selector)
    int set_variable_metric_selection(VariableMetricSelection* vms = 0) {
      delete vm_selection; vm_selection = vms; return 0;
    }

    /// return pointer to an installed VariableMetricSelection object or NULL
    VariableMetricSelection* get_variable_metric_selection() const {
      return vm_selection;
    }


    /** @brief the SumBundleHandler calls this for selecting the next minorants for a polyhedral model

      @param[out] model_indices
          the indices of minorants selected for the model; index 0 is always the aggregate indicated by the input index aggr_index

      @param[in] cand_id
          the identifier of the candidate point supplied next

      @param[in] cand_y
          the candidate (differnte from the center), close to it the model should be good

      @param[in] minorants
          the vector of MinorantPointer gives the minorants
    out of which the model should be selected.

      @param[in] aggr_index
          the index of the aggregate within the minorants

      @param[in] model_maxviol
          a minorant violated by this would have caused a null step

      @param[in] H
          the proximal term used for determining the given cand_y

      @param[in] model_update
          informs about whether cand_y is the result of a null_step or descent_step or aomw other a new set up.

      @return
       - 0 on success
       - 1 on failure
    */
    virtual int select_model(CH_Matrix_Classes::Indexmatrix& model_indices,
      CH_Matrix_Classes::Integer cand_id,
      const CH_Matrix_Classes::Matrix& cand_y,
      const MinorantBundle& minorants,
      CH_Matrix_Classes::Integer aggr_index,
      CH_Matrix_Classes::Real model_maxviol,
      BundleProxObject& H,
      BundleModel::ModelUpdate model_update) = 0;


  };

  //@}

}

#endif

