/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/ConeModel.hxx
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



#ifndef CONICBUNDLE_CONEMODEL_HXX
#define CONICBUNDLE_CONEMODEL_HXX

/**  @file ConeModel.hxx
    @brief Header declaring the class ConicBundle::ConeModel
    @version 1.0
    @date 2017-09-28
    @author Christoph Helmberg
*/

#include "SumBlockModel.hxx"

namespace ConicBundle {

  /** @ingroup InternalBundleModel
  */
  //@{

  /**@brief provides a general abstract base class for conic cutting models with SumBundle capabilities  (as a derived class of SumBlockModel)

     the rountines deal with all administrative aspects of
     SumBundle and SumBundleHandler and call local versions
     for the oracle and cutting model specific parts.

     Note, this base class only allows for either having or contributing
     to an active SumBundle model or for having an active local model,
     but not both can be active at the same time. Before initialization
     it may happen that both are inactive, but once the model is
     initialized, exactly one of data.sumbundle.active() and @a
     data.local_model_active will be true (the correct value will be
     ensured by ConeModel). This still allows that both models are
     updated simultaneously, sumbundle by ConeModel and the full cutting
     model by the derived class, if this is so desired. Usually this is
     difficult because the QP solution is missing for one of them.

     The local model is under full control of the derived class; the
     derived class also has to manage its own QPModelDataObject. In
     terms of the aggregate the derived class has to use and manipulate
     @a data.local_model_aggregate exclusively and only if it is active.
     If the sumbundle is active, also the @a data.local_model_aggregate
     will be handled ConeModel.  The full aggregate data.model_aggregate
     is always and exclusively handled by ConeModel.

     For the data.sumbundle the only point of contact for the derived class
     is sumbundle_mode(); see the explanations given there.

   */


  class ConeModel : public SumBlockModel {
  protected:


    //----------------------------------------------------------------------
   /** @name abstract methods to be filled with oracle and cutting model specific content  */
   //@{ 


   /// see SumBlockModel::eval_function() but ignore ub_fid and ub with function_factor and function_task, this is set outside. Set the data.cand* information and in particular data.cand_ub_mid to data.get_modification_id() if the evaluation results make some sense even in the presence of errors. Please measure the oracle's SumBlockModel::eval_time as well as SumBlockModel::preeval_time and SumBlockModel::posteval_time here if this information is of interest
    virtual int evaluate_oracle(CH_Matrix_Classes::Integer y_id,
      const CH_Matrix_Classes::Matrix& y,
      CH_Matrix_Classes::Real nullstep_bound,
      CH_Matrix_Classes::Real relprec) = 0;

    /// see SumBlockModel::eval_model() but ignore function_factor and function_task, this is set by ConeModel::eval_model(); note: evalmodel_time is already kept by ConeModel
    virtual int evaluate_local_model(CH_Matrix_Classes::Real& lb,
      CH_Matrix_Classes::Integer y_id,
      const CH_Matrix_Classes::Matrix& y,
      CH_Matrix_Classes::Real relprec) = 0;

    /** @brief see BundleModel::update_model() for the general setting;

        The routine is called irrespective of whether the local model is
        active or not; in particular, data.do_step() has to be called
        here if model_update==descent_step (this has the purpose to
        allow the local model to use the data of the old center before
        it is overwritten by the current candidate);

        model_maxviol and H are as passed to ConeModel by BundleModel::update_model();

        note: updatemodel_time is already kept by ConeModel
    */
    virtual int update_local_model(ModelUpdate model_update,
      CH_Matrix_Classes::Integer center_id,
      const CH_Matrix_Classes::Matrix& center_y,
      CH_Matrix_Classes::Integer y_id,
      const CH_Matrix_Classes::Matrix& y,
      CH_Matrix_Classes::Real model_maxviol,
      BundleProxObject& H) = 0;


    /// the local model has to take over after sumbundle was active; current_model holds the final model and coefficients so that the local model can take over whatever it needs, in particualr the aggregate
    virtual int switch_to_local(const MinorantBundle& sumbundle_model,
      const CH_Matrix_Classes::Matrix& sumbundle_coeff) = 0;

    /// the sumbundle is initialized (not necessarily switched on) and needs to know the normalized aggregate of the local model without the influence of the function_factor as well as the coefficient reflecting the influnce of the function_factor; this is called in alternative to start_local_augmodel, so some information has to be available
    virtual int get_aggregate_and_coefficient(MinorantPointer& sumbundle_aggr,
      CH_Matrix_Classes::Real& aggr_coeff) = 0;

    /// see SumBlockModel::start_augmodel(); the local QP block is exclusively organized here; the routine is only called, if the local modle is active
    virtual int start_local_augmodel(QPModelDataPointer& qpblockp,
      CH_Matrix_Classes::Integer cand_id,
      const CH_Matrix_Classes::Matrix& cand_y,
      const CH_Matrix_Classes::Indexmatrix* needed_indexsubset) = 0;

    /// see SumBlockModel::provide_model_aggregate(), only called if local model is active
    virtual int provide_local_aggregate(CH_Matrix_Classes::Integer y_id,
      const CH_Matrix_Classes::Matrix& iny) = 0;

    /// see SumBlockModel::make_model_aggregate for the general rationale; the routine is only called, if the local modle is active; if the local model decides to increase the function factor, it reports back this value in new_function_factor (only allowed if keep_penalty_fixed==false); the new data.function_factor and the modification ids are the set in ConeModel
    virtual int make_local_aggregate(CH_Matrix_Classes::Real& new_function_factor,
      bool keep_penalty_fixed) = 0;

    /// see SumBlockModel::adjust_multiplier();  the routine is only called for adaptive penalty functions and if the local model is active; if the local model decides to adapt the function factor (here it is typically decreased), it reports back this value in new_function_factor; the new data.function_factor and the modification ids are the set in ConeModel
    virtual int local_adjust_multiplier(CH_Matrix_Classes::Real& new_function_factor) = 0;

    /// returns true if candidate information (always the results of the most recent function evaulation) is available and up to date (no model changes in the meantime) 
    virtual bool candidate_available() = 0;
    //e.g. { return (data.cand_ub_mid==data.modification_id)&&(data.cand_Ritzvec.coldim()>0); }

    /// if it cannot be provided it returns a minorant with valid()==false
    virtual MinorantPointer& provide_cand_minorant() = 0;

    /// if it cannot be provided it returns a minorant with valid()==false
    virtual MinorantPointer& provide_center_minorant() = 0;


    //@}




  public:

    /// resets all data to initial status of this class, also the bundle parameters
    virtual void clear();

    // set_variable_metric_selection     //as in SumBlockData

    /// construct a model for the MatrixFunctionOracle pointed to by @a fo
    ConeModel(CBout* cb = 0, int incr = -1);

    /// calls clear() and deletes the bundle parameters
    virtual ~ConeModel();


    //----------------------------------------------------------------------
    /** @name implementations of abstract class BundleModel (maybe overloading some of SumBlockModel) */
    //@{


    /// see BundleModel::eval_function
    int eval_function(CH_Matrix_Classes::Integer& ub_fid,
      CH_Matrix_Classes::Real& ub,
      CH_Matrix_Classes::Integer y_id,
      const CH_Matrix_Classes::Matrix& y,
      CH_Matrix_Classes::Real nullstep_bound,
      CH_Matrix_Classes::Real relprec);

    /// see BundleModel::eval_model
    int eval_model(CH_Matrix_Classes::Real& lb,
      CH_Matrix_Classes::Integer y_id,
      const CH_Matrix_Classes::Matrix& y,
      CH_Matrix_Classes::Real relprec);

    //get_model_aggregate(three parameters) //as in SumBlockModel

    /// see BundleModel::update_model
    int update_model(ModelUpdate model_update,
      CH_Matrix_Classes::Integer center_id,
      const CH_Matrix_Classes::Matrix& center_y,
      CH_Matrix_Classes::Integer cand_id,
      const CH_Matrix_Classes::Matrix& cand_y,
      CH_Matrix_Classes::Real model_maxviol,
      BundleProxObject& H);

    //synchronize_ids       //as in SumBlockModel

    //center_modified        //as in SumBlockModel

    /// see BundleModel::recompute_center
    int recompute_center(CH_Matrix_Classes::Integer& new_center_fid,
      CH_Matrix_Classes::Real& new_center_ub,
      CH_Matrix_Classes::Integer center_id,
      const CH_Matrix_Classes::Matrix& y,
      bool accept_only_higher_values = false,
      CH_Matrix_Classes::Real relprec = -1.);



    /// see BundleModel::model_aggregate_modified
    bool model_aggregate_modified(CH_Matrix_Classes::Integer old_model_aggregate_id) {
      return get_data()->model_aggregate_modified(old_model_aggregate_id);
    }

    /// see BundleModel::provide_model_aggregate
    int provide_model_aggregate(CH_Matrix_Classes::Integer y_id,
      const CH_Matrix_Classes::Matrix& y);

    // apply_modification     // as in SumBlockModel

    /// see BundleModel::check_center_validity_by_candidate
    int check_center_validity_by_candidate(bool& cand_minorant_is_below,
      CH_Matrix_Classes::Integer center_id,
      const CH_Matrix_Classes::Matrix& center_y);

    //@}

    //----------------------------------------------------------------------
    /** @name implemented/overloaded routines of SumBlockModel */
    //@{

    // get_oracle_object() still abstract

    // start_augmodel see SumBlockModel::start_augmodel

    /// see SumBlockModel::make_model_aggregate
    int make_model_aggregate(bool& penalty_parameter_increased,
      bool keep_penalty_fixed);

    //get_model_aggregate(... aft=0)   //as in SumBlockModel

    /// see SumBlockModel::lb_function
    CH_Matrix_Classes::Real lb_function(CH_Matrix_Classes::Integer y_id,
      const CH_Matrix_Classes::Matrix& y);

    /// see SumBlockModel::get_function_minorant()
    int get_function_minorant(MinorantPointer& minorant,
      const AffineFunctionTransformation* aft = 0);

    /// see SumBlockModel::get_center_minorant()
    int get_center_minorant(MinorantPointer& minorant,
      const AffineFunctionTransformation* aft = 0);

    /// see SumBlockModel::adjust_multiplier()
    int adjust_multiplier(bool& values_may_have_changed);

    /// see SumBlockModel::sumbundle_mode
    int sumbundle_mode(SumBundle::Mode& mode,
      SumBundleHandler* bh = 0,
      AffineFunctionTransformation* aft = 0);

    ///see BundleModel::start_augmodel() for the first four parameters; for the others see sumbundle_contribution()
    virtual int start_sumaugmodel(QPModelDataPointer& blockp,
      CH_Matrix_Classes::Integer cand_id,
      const CH_Matrix_Classes::Matrix& cand_y,
      const CH_Matrix_Classes::Indexmatrix* indices = 0,
      SumBundleHandler* bh = 0,
      SumBundle::Mode mode = SumBundle::inactive,
      AffineFunctionTransformation* aft = 0);


    /// see SumBlockModle::update_model(ModelUpdate,CH_Matrix_Classes::Integer,const CH_Matrix_Classes::Matrix&,CH_Matrix_Classes::Real,CH_Matrix_Classes::Real,CH_Matrix_Classes::Real&,CH_Matrix_Classes::Real&)
    virtual int update_model(ModelUpdate model_update,
      CH_Matrix_Classes::Integer center_id,
      const CH_Matrix_Classes::Matrix& center_y,
      CH_Matrix_Classes::Integer cand_id,
      const CH_Matrix_Classes::Matrix& cand_y,
      CH_Matrix_Classes::Real model_maxviol,
      BundleProxObject& H,
      CH_Matrix_Classes::Real& model_deviation,
      CH_Matrix_Classes::Real& model_curvature);


    // SumBlockModel::get_data is still abstract

     //@}

     //----------------------------------------------------------------------
     /** @name messages for direct get/set requestss */
     //@{


     // SumBlockModel::get_approximate_primal is still abstract

     // SumBlockModel::get_center_primal is still abstract

     // SumBlockModel::get_candidate_primal is still abstract

     // SumBlockModel::call_primal_extender is still asbtract

     // SumBlockModel::set_bundle_parameters is still abstract

     // SumBlockModel::get_bundle_parameters is still abstract

     // see SumBlockModel::clear_model 
    void clear_model(bool discard_minorants_only);

    /// see SumBlockModel::clear_aggregates
    void clear_aggregates();

    /// SumBlockModel::get_ret_code() is still abstract


      //@}

    //get_QPcoeff_time()      // as in SumBlockModel

    //get_QPsolve_time()      // as in SumBlockModel

    //get_evalaugmodel_time() // as in SumBlockModel

    //get_evalmodel_time()    // as in SumBlockModel

    //get_updatemodel_time()  // as in SumBlockModel

    //get_eval_time()         // as in SumBlockModel

    //set_out()               // as in SumBlockModel

    //set_cbout()  // calls set_out in SumBlockModel



  };

  //@}

}

#endif

