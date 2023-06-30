/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/SumModel.hxx
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



#ifndef CONICBUNDLE_SUMMODEL_HXX
#define CONICBUNDLE_SUMMODEL_HXX

/**  @file SumModel.hxx
    @brief Header declaring the class ConicBundle::SumModel
    @version 1.0
    @date 2014-07-24
    @author Christoph Helmberg
*/

#include <list>
#include "SumBlockModel.hxx"
#include "AFTModel.hxx"

namespace ConicBundle {
  /** @ingroup InternalBundleModel

  */
  //@{

  class SumModelParametersObject;

  /**@brief gives a SumBlockModel interface to a sum of (possibly
     AffineFunctionTransformation transformed) SumBlockModels. This
     allows e.g. to optimize over a nonnegative combination of convex
     functions, possibly each with its own cutting model.

     Here the @a modification_id of BundleData is increased whenever
     a model is added or deleted from the sum.

     If SumModel contains many separate functions, having a separate
     model for each of the functions may be unnecessary and inefficient.
     Thus SumModel contains a heuristic (may be disabled) that decides
     which of the functions may keep its separate model and which is
     simply providing its subgradient information to feed a common
     polyhedral model for the remaining functions. The heuristic decides
     on the fly, whether it is worth to switch for a certain function to
     a local model or to the common polyhedral model. Thus each function
     needs to keep track of the subgradients and the aggregate in the
     global polyhedral model once it is told do so, so that switching
     back and forth does not endanger convergence. This requires some
     communication and SumBlockModel provides the interface for this,
     see the routines
     - SumBlockModel::sumbundle_contribution()
     - SumBlockModel::remove_contributions()

     Due to the different types of functions (cost and penalty) SumModle
     has to use of a different global model for each type appearing. For
     keeping consistent with the recursive structure it also needs to
     assume that each of its functions has these three different types
     of models. This makes the whole business a bit clumsy ...

   */


  class SumModel : public SumBlockModel {
  public:
    //==================  models description ====================
    /** @brief For each model added to the sum, there is an object ModelData that stores required information

        In addition to the pointer @a oracle_model to the model itself, it holds
        a pointer @a aft_model to an optinal AFTModel that allows to apply an
        AffineFunctionTransformation (AFT) to the function of this model. If no
        AFT is used, the @a aft_model pointer will NULL. An AFT may be specified
        when setting up the problem or it will be generated automatically if a
        FunctionObjectModification for this model's function specifies a
        modification for an AFT.  The AFTModel will be deleted in the end and,
        by default, also the AFT it contains. If initially an AFT is provided at
        construction and this AFT should not be deleted, set @a
        delete_this_aft_on_exit to false.

        The only other variable is @a aggregate_id which is -1 if not
        initialized or gives the identifier returned by the last call
        to get_model_aggregate.

     */
    class ModelData {
      friend SumModel;   ///< SumModel may access any of this directly
    private:
      /// this oracle model is part of the sum, it will NOT be deleted on exit
      SumBlockModel* oracle_model;
      /// -1 if not initialized, otherwise the maximum of the identifiers returned by the various calls to evaluation routines
      CH_Matrix_Classes::Integer function_id;
      /// -1 if not initialized, otherwise the identifier returned by the call giving the center_ub_value
      CH_Matrix_Classes::Integer center_ub_fid;
      /// -1 if not initialized, otherwise the identifier returned by the calls giving the cand_ub_value
      CH_Matrix_Classes::Integer cand_ub_fid;
      /// -1 if not initialized, otherwise the identifier returned by the last call to get_model_aggregate
      CH_Matrix_Classes::Integer aggregate_id;

      ///for judging whether to use a local model or to suggest contribution to summodel, function_value(y)-aggregate(y) as returned by the last call to the extended SumBlockModel::update_model()
      CH_Matrix_Classes::Real model_deviation;

      ///for judging whether to use a local model or to suggest contribution to summodel, this is changed/updated by  SumModelSelection::select_model()
      CH_Matrix_Classes::Real deviation_contribution;

      ///for judging whether to use a local model or to suggest contribution to summodel, an estimate of returned by the last call to the extended SumBlockModel::update_model()
      CH_Matrix_Classes::Real model_curvature;

      ///for judging whether to use a local model or to suggest contribution to summodel, this is changed/updated by  SumModelSelection::select_model()
      CH_Matrix_Classes::Real weighted_curvature;

      ///suggested mode for this model set by the last call to SumModelSelection::select_model()
      SumBundle::Mode suggested_mode;

      /// returns @a aft_model if this is not NULL, otherwise @a oracle_model
    public:
      SumBlockModel* model() const {
        return oracle_model->sbm_transform();
      }
      ///
      ~ModelData();
      /// sets the oracle_model
      ModelData(SumBlockModel* oracle_model);

      /// get recent difference between model and function value in candidate
      CH_Matrix_Classes::Real get_model_deviation() const {
        return model_deviation;
      }

      /// get recent estimate of the function's "curvature" in candidate
      CH_Matrix_Classes::Real get_model_curvature() const {
        return model_curvature;
      }

      /// get and set the accumulated contribution of this function to the total deviation
      CH_Matrix_Classes::Real& set_deviation_contribution() {
        return deviation_contribution;
      }

      /// get and set the accumulated curvature information of this function
      CH_Matrix_Classes::Real& set_weighted_curvature() {
        return weighted_curvature;
      }

      /// get and set the sumbundle mode suggested to the model of the function
      SumBundle::Mode& set_suggested_mode() {
        return suggested_mode;
      }
    };

    /// a ModelMap maps FunctionObject pointers (these must be unique) to their ModelData
    typedef std::map<const FunctionObject*, ModelData*> ModelMap;

  private:
    //==================  dummy oracle  =====================
    /// the oracle object of this model
    ModifiableOracleObject dummy_oracle;

    //--- data describing the model
    /// the relevant data for running BundleSolver consistently
    BundleData data;

    /// this holds the ModelData indexed by the (unique) FunctionObject pointer  
    ModelMap modelmap;

    /// if the sumbundle is active, this collects all models that have a separate augmented model part and need to be called there 
    std::list<const ModelData*> local_augmodel;

    //=================== data of function evaluations ==================
    ///number of (recursive) calls to eval_function or recompute_center
    CH_Matrix_Classes::Integer ncalls;

    /// the minorant in the candidate (if collected explicitly)
    MinorantPointer cand_minorant;

    /// the minorant in the center (if collected explicitly)
    MinorantPointer center_minorant;


    //===================  quadratic subproblems ==================
    /// collects the blocks of the submodels 
    QPSumModelDataObject* block;

    /// parameters and routines for choosing the models for SumBundle
    SumModelParametersObject* model_selection;


  public:
    /// resets the SumModel to its initial state, in particular it removes but does not delete any Models added in add_model() (their AFTs are deleted if not explicitly denied) 
    void clear();

    /// constructor
    SumModel(CBout* cb = 0);

    /// calls clear()
    ~SumModel();

    ///returns the submodel for FunctionObject fo if it in this model, otherwise 0
    const SumBlockModel* model(const FunctionObject* fo) const;

    ///returns the number of submodels in this model
    CH_Matrix_Classes::Integer nsubmodels() const {
      return CH_Matrix_Classes::Integer(modelmap.size());
    }

    /// adds the @a model as submodel to this model (if this model may have submodels); if any error occurs the model is not added
    int add_model(SumBlockModel* model);

    /// remove the submodel identified by @a fo from this model, this does NOT destruct the model. It returns the pointer to the model if there is one, otherwise 0 
    SumBlockModel* remove_model(const FunctionObject* fo);


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

    //eval_augmodel           //as in SumBlockModel

    //reeval_augmodel         //as in SumBlockModel

    //get_model_aggregate(three parameters) //as in SumBlockModel

    /// see BundleModel::update_model
    int update_model(ModelUpdate model_update,
      CH_Matrix_Classes::Integer center_id,
      const CH_Matrix_Classes::Matrix& center_y,
      CH_Matrix_Classes::Integer cand_id,
      const CH_Matrix_Classes::Matrix& cand_y,
      CH_Matrix_Classes::Real model_maxviol,
      BundleProxObject& H);

    /// see BundleModel::synchronize_ids
    int synchronize_ids(CH_Matrix_Classes::Integer& new_center_ub_fid,
      CH_Matrix_Classes::Integer new_center_id,
      CH_Matrix_Classes::Integer old_center_id,
      CH_Matrix_Classes::Integer& new_cand_ub_fid,
      CH_Matrix_Classes::Integer new_cand_id,
      CH_Matrix_Classes::Integer old_cand_id,
      CH_Matrix_Classes::Integer& new_aggregate_id);

    /// see BundleModel::center_modified
    bool center_modified(CH_Matrix_Classes::Integer& function_modification_id,
      CH_Matrix_Classes::Integer center_id);

    /// see BundleModel::recompute_center
    int recompute_center(CH_Matrix_Classes::Integer& new_center_ub_fid,
      CH_Matrix_Classes::Real& new_center_ub,
      CH_Matrix_Classes::Integer center_id,
      const CH_Matrix_Classes::Matrix& y,
      bool accept_only_higher_values = false,
      CH_Matrix_Classes::Real relprec = -1.);

    /// see BundleModel::model_aggregate_modified
    bool model_aggregate_modified(CH_Matrix_Classes::Integer old_model_aggregate_id);

    /// see BundleModel::provide_model_aggregate
    int provide_model_aggregate(CH_Matrix_Classes::Integer y_id,
      const CH_Matrix_Classes::Matrix& y);

    /// see BundleModel::apply_modfication()
    int apply_modification(bool& no_changes,
      const GroundsetModification& gsmdf,
      const FunObjModMap& funmdfmap,
      CH_Matrix_Classes::Integer new_center_id,
      const CH_Matrix_Classes::Matrix& new_center,
      CH_Matrix_Classes::Integer old_center_id,
      const CH_Matrix_Classes::Matrix& old_center);


    /// see BundleModel::add_dynamic_scaling
    int add_variable_metric(VariableMetric& H,
      CH_Matrix_Classes::Integer y_id,
      const CH_Matrix_Classes::Matrix& y,
      bool descent_step,
      CH_Matrix_Classes::Real weightu,
      CH_Matrix_Classes::Real model_maxviol,
      const CH_Matrix_Classes::Indexmatrix* indices = 0);

    /// see BundleModel::check_center_validity_by_candidate
    int check_center_validity_by_candidate(bool& cand_minorant_is_below,
      CH_Matrix_Classes::Integer center_id,
      const CH_Matrix_Classes::Matrix& center_y);

    //@}

    //----------------------------------------------------------------------
    /** @name implemented/overloaded routines of SumBlockModel */
    //@{

    /// as AFTModel has no oracle of its own, this returns the dummy oracle
    ModifiableOracleObject* get_oracle_object() {
      return &dummy_oracle;
    }

    // start_augmodel as in SumBlockModel::start_augmode()l

    /// see SumBlockModel::make_model_aggregate
    int make_model_aggregate(bool& penalty_parameter_increased,
      bool keep_penalty_fixed);

    // get_model_aggregate(six parameters)     //as in SumBlockModel 

    /// see SumBlockModel::lb_function()
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


    /// see SumBlockModle::update_model(ModelUpdate,CH_Matrix_Classes::Integer,const CH_Matrix_Classes::Matrix&,CH_Matrix_Classes::Integer,const CH_Matrix_Classes::Matrix&,BundleProxObject&,CH_Matrix_Classes::Real&,CH_Matrix_Classes::Real&)
    int update_model(ModelUpdate model_update,
      CH_Matrix_Classes::Integer center_id,
      const CH_Matrix_Classes::Matrix& center_y,
      CH_Matrix_Classes::Integer cand_id,
      const CH_Matrix_Classes::Matrix& cand_y,
      CH_Matrix_Classes::Real null_step_factor,
      BundleProxObject& H,
      CH_Matrix_Classes::Real& model_deviation,
      CH_Matrix_Classes::Real& model_curvature);

    /// see SumBlockModel::get_data()
    BundleData* get_data() {
      return &data;
    }

    /// see SumBlockModel::get_data()
    const BundleData* get_data() const {
      return &data;
    }

    /// see SumBlockModel::set_data()
    int set_data(BundleData* bd) {
      return data.init(bd);
    }

    //@}


    //----------------------------------------------------------------------
    /** @name overloaded SumBlockModel messages for direct get/set requests */
    //@{

    /// this has no primals, so it returns 1, see SumBlockModel::get_approximate_primal
    const PrimalData* get_approximate_primal() const {
      return 0;
    }

    /// this has no primals, so it returns 1, see SumBlockModel::get_center_primal
    const PrimalData* get_center_primal() const {
      return 0;
    }

    /// this has no primals, so it returns 1, see SumBlockModel::get_candidate_primal 
    const PrimalData* get_candidate_primal() const {
      return 0;
    }

    /// this has no primals, so it returns 1, see SumBlockModel::call_primal_extender 
    virtual int call_primal_extender(PrimalExtender& /* pext */) {
      return 1;
    }

    /// if bp is a SumBundleParametersObject, it is cloned, otherwise the paraemters are copied
    int set_bundle_parameters(const BundleParameters& bp);

    ///  see SumBlockModel::get_bundle_parameters()  
    const BundleParameters* get_bundle_parameters() const;

    /// see SumBlockModel::clear_model
    void clear_model(bool discard_minorants_only = false);

    /// see SumBlockModel::clear_aggregates
    void clear_aggregates();

    //get_ret_code()      // as in SumBlockModel

    //@}

    //---------------------------------------------------------------------- 
    /**@name overloaded routines of SumBlockModel for output and statistics on solution times */
    //@{

    //get_QPcoeff_time()      // as in SumBlockModel

    //get_QPsolve_time()      // as in SumBlockModel

    //get_evalaugmodel_time() // as in SumBlockModel

    //get_evalmodel_time()    // as in SumBlockModel

    //get_updatemodel_time()  // as in SumBlockModel

    //get_eval_time()         // as in SumBlockModel

    /// set output and outputlevel of warnings and errors recursively, see CBout
    void set_out(std::ostream* o = 0, int pril = 1);

    //set_cbout()  // calls set_out in SumBlockModel

    /// see SumBlockModel::get_preeval_time()
    virtual CH_Tools::Microseconds get_preeval_time() const {
      CH_Tools::Microseconds ms(preeval_time);
      for (ModelMap::const_iterator it = modelmap.begin(); it != modelmap.end(); it++) {
        ms += it->second->model()->get_preeval_time();
      }
      return ms;
    }
    /// see SumBlockModel::get_eval_time()
    virtual CH_Tools::Microseconds get_eval_time() const {
      CH_Tools::Microseconds ms(0);
      for (ModelMap::const_iterator it = modelmap.begin(); it != modelmap.end(); it++) {
        ms += it->second->model()->get_eval_time();
      }
      return ms;
    }

    /// see SumBlockModel::get_posteval_time()
    virtual CH_Tools::Microseconds get_posteval_time() const {
      CH_Tools::Microseconds ms(posteval_time);
      for (ModelMap::const_iterator it = modelmap.begin(); it != modelmap.end(); it++) {
        ms += it->second->model()->get_posteval_time();
      }
      return ms;
    }


    //@}


  };


  //@}

}

#endif

