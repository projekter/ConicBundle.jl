/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/SumBlockModel.hxx
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



#ifndef CONICBUNDLE_SUMBLOCKMODEL_HXX
#define CONICBUNDLE_SUMBLOCKMODEL_HXX


/**  @file SumBlockModel.hxx
    @brief Header declaring the class ConicBundle::SumBlockModel
    @version 1.0
    @date 2014-07-14
    @author Christoph Helmberg
*/

#include "clock.hxx"
#include "MatrixCBSolver.hxx"
#include "BundleModel.hxx"
#include "AffineFunctionTransformation.hxx"
#include "BundleData.hxx"
#include "SumBundleHandler.hxx"

namespace ConicBundle {


  /** @ingroup InternalBundleModel

  */
  //@{

  class AFTModel;

  /** @brief abstract interface extending BundleModel so that any such
  model can be used alone or within SumModel and so that it supports
  AffineFunctionTransformation as well as switching to a SumBundle with
  some SumBundleHandler

  The BundleSolver will alway assume that there is one cutting model and it does
  not need to know how this looks like. A SumBlockModel like the SumModel is
  used to present this single model view to the BundleSolver, yet it allows
  e.g. in the case of SumModel for a sum of convex functions to use a separate
  cutting model for each of the functions. Of course it is also possible to set
  up a common single model via an appropriate oracle or to recursively build a
  sum of models of sums of models, the latter with some loss of efficiency. The
  joint quadratic bundle subproblem will consist of a quadratic cost matrix
  coupling all models (start_augmodel() serves to collect the bundle model,
  BundleProxObject forms the quadratic data from this), but the generating
  descriptions for the linear minorants of each model will each be a separate
  block of constraints (obtained in start_augmodel()). After the
  quadratic subproblem is solved the local aggregates are formed via
  make_model_aggregate(), afterwards the joint aggregate of the entire sum will
  need to be collected by recursive calls to get_model_aggregate().

  If an objective is a sum of convex functions, most functions will just work on
  part of the coordinates or on an affine transformation of the
  coordinates. Here we implement this by supplying each model with the
  possibility to be combined with an AFTModel, which first performs an
  AffineFunctionTransformation of the input arguments.

  Depending on its ConicBundle::FunctionTask each function may be used as a
  simple cost function or it may be of penalty type for its level set to the
  value zero. Functions of penalty type may have a prespecified constant penalty
  parameter or an adaptive penalty parameter that may be increased or decreased
  upon need (dynamic updates of the penalty parameter take place in
  make_model_aggregate() and adjust_multiplier()).  Intuitively speeking, the
  penalty parameter allows to increase the length of subgradients, thereby
  increasing their slope so that the iterates get forced to lie within the level
  set in the limit. Whenever the feasible set is compact and the functions are
  finite valued, a finite penalty value is sufficient (exact penalty) if on the
  boundary of the level set of the penalty function the slope is bounded away
  from zero.

  If SumModel contains many separate functions, having a separate model for each
  of the functions may be unnecessary and inefficient.  Making use of additional
  information provided in the extended routine update_model(ModelUpdate,CH_Matrix_Classes::Integer,const CH_Matrix_Classes::Matrix&,CH_Matrix_Classes::Real&,CH_Matrix_Classes::Real&,CH_Matrix_Classes::Real&,CH_Matrix_Classes::Real&,CH_Matrix_Classes::Real&)
  the class SumModel offers a heuristic that suggests which of the functions
  should keep a separate model and which should only contribute its subgradient
  information to a common polyhedral model for the remaining functions. The
  heuristic gives a new suggestion on the fly whenever the next model is
  formed. Thus each function needs to keep track of the subgradients and the
  aggregate in the global polyhedral model should it want to participate in the
  common global model, so that switching back and forth does not endanger
  convergence. Note that no model is forced to contribute and that models may
  even decide to contribute one part but keep another part as a local model. Due
  to the different types of functions (fixed, bounded, dynamic) SumModel has to
  use a different common model for each type appearing. For keeping consistent
  with the recursive structure it also needs to assume that each of its
  functions has these three different types of models. This makes the whole
  business a bit clumsy. For uniform treatment, each Model is equipped with a
  SumBundle, and this is employed only if the Model is given a SumBundleHandler
  and informed to use it via the routine sumbundle_contribution(). If a
  model/function wants to stop contributing to a SumBundle of some parent model
  without destroying the bundle information of the others e.g. because it is
  modified, it can do so by calling sumbundle_remove_contributions().

  The class also offers a number of information collection routines
  that are not needed for running BundleSolver but are
  requested in many applications, e.g., for providing the common
  subgradient in the center, information about the PrimalData
  generating the aggregate, or some computation time measurements.

  */

  class SumBlockModel : public BundleModel {
  protected:
    CH_Tools::Clock clock;      ///< for collecting time statistics
    CH_Tools::Microseconds evalmodel_time; ///< total time spent in evaluating the model in eval_model()
    CH_Tools::Microseconds updatemodel_time; ///< total time spent in updating the bundle in update_model()
    CH_Tools::Microseconds eval_time; ///< total time spent in the oracle in eval_function() or calls to eval_function of children
    CH_Tools::Microseconds preeval_time; ///< total time spent in eval_function() before the oracle call
    CH_Tools::Microseconds posteval_time; ///< total time spent in eval_function() after the oracle call
    CH_Tools::Microseconds metric_time; ///< total time spent in add_variable_metric()

    //--- every SumBlockModel can be equipped with its own affine function transformation
    AFTModel* aftmodel; ///< if not NULL this points to an AFTModel describing an AffineFunctionTransformation

    //--- every SumBlockModel can be equipped with its own VariableMetricSelection
    VariableMetricSelection* vm_selection; ///< if not NULL this points to the selection routine for computing the local contribution to the metric

    //--- the bundle handler for contributing to a parent sumbundle
    /// if the sumbundle is used, this points to the handler that operates it
    SumBundleHandler* bundlehandler;
    /// this holds the default settings whether and how to start or contribute to a SumBundle
    SumBundleParametersObject* sumbundle_parameters;

    //--- for testing purposes
    /// for testing purposes
    mutable MinorantPointer old_model_aggregate;


  public:
    /// destructor
    virtual ~SumBlockModel();

    /// resets all data to the initial status of this class, also @a aftmodel and @a vm_selection are deleted 
    virtual void clear();

    /// first it discards an old affine function transformation if there is one, then it sets the new one, if there ins one. 
    int initialize_aft(AffineFunctionTransformation* aft = 0);

    /// delete old selector and set a new one (0 is allowed resulting in no local selector)
    virtual int set_variable_metric_selection(VariableMetricSelection* vms = 0) {
      delete vm_selection; vm_selection = vms; return 0;
    }

    /// delete old selector and set a new one (0 is allowed resulting in no local selector)
    virtual VariableMetricSelection* get_variable_metric_selection() const {
      return vm_selection;
    }

    /** @brief in apply_modification this routine is needed to check whether the aftmodel is modified already

       If an aftmodel is present or fundmdfmap contains no data for the aftmodel
       of this, the routine returns false. If, however, no aftmodel is present
       so far and in spite of this funmdfmap contains modifications for
       it, these were no yet carried out but need to be carried out before the
       data is passed on to *this. So in this case the routine initializes
       an aftmodel and returns true (otherwise it always returns false).
       Whenever this routine returns true, the apply_modification routine
       has to pass the call on to the new aftmodell first so that it is
       called itself again with the correct parameters.

     */
    bool call_aftmodel_first(const FunObjModMap& funmdfmap);

    /// calls clear
    SumBlockModel(CBout* cb = 0, int cbinc = -1);

    /// return the function oracle this provides a model for, or some dummy oracle
    virtual ModifiableOracleObject* get_oracle_object() = 0;

    ///returns the submodel for FunctionObject fo if it in this model, otherwise 0
    virtual const SumBlockModel* model(const FunctionObject* /* fo */) const {
      return 0;
    }

    ///returns the number of submodels in this model (direct ones, not all descendants; most have none, but see e.g. SumModel)
    virtual CH_Matrix_Classes::Integer nsubmodels() const {
      return 0;
    }

    /// adds the @a model as submodel to this model (if this model may have submodels) 
    virtual int add_model(SumBlockModel* /* model */) {
      return 1;
    }

    /// remove the submodel identified by @a fo from this model, this does NOT destruct the model. It returns the pointer to the model if there is one, otherwise 0 
    virtual SumBlockModel* remove_model(const FunctionObject* /* fo */) {
      return 0;
    }

    /// remove the submodel identified by the given models FunctionObject from this model, this does NOT destruct the model. It returns the pointer to the model if there is one, otherwise 0 
    virtual SumBlockModel* remove_model(SumBlockModel* model) {
      if (model != 0)
        return remove_model(dynamic_cast<FunctionObject*>(model->get_oracle_object()));
      return 0;
    }



    //----------------------------------------------------------------------
    /** @name some default implementations of abstract class BundleModel */
    //@{

    //virtual int eval_function(...)=0; is still abstract 

    //virtual int eval_model(...) =0;  is still abstract

    ///see BundleModel::start_augmodel, here it just moves on to start_sumaugmodel
    virtual int start_augmodel(QPModelDataPointer& blockp,
      CH_Matrix_Classes::Integer cand_id,
      const CH_Matrix_Classes::Matrix& cand_y,
      const CH_Matrix_Classes::Indexmatrix* indices = 0);


    // virtual int start_augmodel(......) =0; is still abstract

    // virtual int make_model_aggregate(..)=0; is still abstract 

    ///see BundleModel::get_model_aggregate
    virtual int get_model_aggregate(CH_Matrix_Classes::Integer& model_aggregate_id,
      MinorantPointer& model_aggregate);

    //see BundleModel::update_model() it is still abstract
    virtual int update_model(ModelUpdate model_update,
      CH_Matrix_Classes::Integer center_id,
      const CH_Matrix_Classes::Matrix& center_y,
      CH_Matrix_Classes::Integer cand_id,
      const CH_Matrix_Classes::Matrix& cand_y,
      CH_Matrix_Classes::Real model_maxviol,
      BundleProxObject& H) = 0;

    ///see BundleModel::synchronize_ids
    virtual int synchronize_ids(CH_Matrix_Classes::Integer& new_center_ub_fid,
      CH_Matrix_Classes::Integer new_center_id,
      CH_Matrix_Classes::Integer old_center_id,
      CH_Matrix_Classes::Integer& new_cand_ub_fid,
      CH_Matrix_Classes::Integer new_cand_id,
      CH_Matrix_Classes::Integer old_cand_id,
      CH_Matrix_Classes::Integer& new_aggregate_id);

    ///see BundleModel::center_modified
    virtual bool center_modified(CH_Matrix_Classes::Integer& center_fid,
      CH_Matrix_Classes::Integer center_id);

    //virtual bool recompute_center(...) =0;  is still abstract

    ///see BundleModel::model_aggregate_modified
    virtual bool model_aggregate_modified(CH_Matrix_Classes::Integer old_model_aggregate_id);

    //virtual int provide_model_aggregate(...) =0; is still abstract

    /** @brief passes on modification information about function and ground set changes

        In addition to the explanation given in BundleModel::apply_modification()
        the variable no_changes also acts as an input variable.

        If on input no_change==false and the model contributes as a child
        to sumbundle, it has to remove these contributions before excuting
        any changes, even if there are no changes within.
    */

    virtual int apply_modification(bool& no_changes,
      const GroundsetModification& gsmdf,
      const FunObjModMap& funmdfmap,
      CH_Matrix_Classes::Integer new_center_id,
      const CH_Matrix_Classes::Matrix& new_center,
      CH_Matrix_Classes::Integer old_center_id,
      const CH_Matrix_Classes::Matrix& old_center);

    /// see VariableMetricModel::add_variable_metric; this implements the non local option
    int add_variable_metric(VariableMetric& H,
      CH_Matrix_Classes::Integer y_id,
      const CH_Matrix_Classes::Matrix& y,
      bool descent_step,
      CH_Matrix_Classes::Real weightu,
      CH_Matrix_Classes::Real model_maxviol,
      const CH_Matrix_Classes::Indexmatrix* indices = 0);

    // virtual int check_center_validity_by_candidate(...)=0; is still abstract

    /// see BundleModel::get_function_minorant()
    int get_function_minorant(CH_Matrix_Classes::Integer& function_modification_id, MinorantPointer& minorant);

    /// see BundleModel::get_center_minorant()
    int get_center_minorant(CH_Matrix_Classes::Integer& function_modification_id, MinorantPointer& minorant);

    /// if an affine function transformation is defined for this model, return it, otherwise use this
    BundleModel* transform();

    //@}

    //----------------------------------------------------------------------
    /**@name new virtual functions needed for the classes AFTModel and SumModel */
    //@{

    /// like BundleModel::transform() this allows SumBlockModel routines to use a given transformation automatically
    SumBlockModel* sbm_transform();



    /** @brief returns the model aggregate if available.

      The current routine differs from BundelModle::get_model_aggregate
      in the two additional parameters @a all_parts, and
      @a aft.  If @a add==false, the aggregate has to be initiliazed. If
      == true, the aggregate has to be added. If @a all_parts==true, the
      entire aggregate is used. If @a all_parts==false, then the local
      parts are used that correspnd to model parts handled by *this
      rather than a sumbundle of a parent. If @a aft!=0, the
      AffineFunctionTransformation pointed to has first to be applied
      and the result has to be stored/added in/to the aggregate.

    @param[out] model_aggregate_id
        nonnegative counter for quickly checking validity of the aggregate

    @param[in,out] model_aggregate
        the aggregate linear model minorant

    @param[in] all_parts
        - if == true, use the entire aggregate
        - if == false, use only the local parts that are not taken care of by a parent sumbundle

    @param[in] aft
        - if == 0 store/add the local aggregate directly
        - if != 0 apply this transformation first to the local aggregate before storing/adding it

    @return
      - 0 if the information is available
          (if eval_augmodel() did not return 0, the information will not
          satisfy the precision requirements but may still be available)
      - 1 if the desired information is not available
    */

    virtual int get_model_aggregate(CH_Matrix_Classes::Integer& model_aggregate_id,
      MinorantPointer& model_aggregate,
      bool all_parts,
      const AffineFunctionTransformation* aft = 0);

    /** @brief returns a __quick__ lower bound for the function value at y (eg by a previous subgradient)
    */
    virtual CH_Matrix_Classes::Real lb_function(CH_Matrix_Classes::Integer y_id,
      const CH_Matrix_Classes::Matrix& y);

    /** @brief returns the minorant corresponding to the subgradient inequality returned by the last function evaluation (if the function has not changed since then). If the input minorant is not empty, the local minorant is added.
     */
    virtual int get_function_minorant(MinorantPointer& function_minorant,
      const AffineFunctionTransformation* aft = 0) = 0;

    /** @brief returns the minorant corresponding to the subgradient inequality returned by the function evaluation for the current center (if the function has not changed since then). If the input minorant is not empty, the local minorant is added.
     */
    virtual int get_center_minorant(MinorantPointer& center_minorant,
      const AffineFunctionTransformation* aft = 0) = 0;


    //-------- messages needed, if online problem modifications are desired
    //         these messages are passed on recursively by SumModel and AFTModel


    /** @brief for conic subproblems with adjustable multipliers, reset the
    multiplier to twice the current trace and set values_may_have_changed to
    true in this case. Otherwise leave values_may_have_changed unmodified.

    This message is recursively passed on. Afterwards a recompute_center has to
    be called if values_may_have_changed is true.

    @param[in,out] values_may_have_changed
        - On input: if use_local_model is false and this value is true,
          the adaptive penalty parameter was changed and needs to be adapted
          here as well.
        - On output: if the function value is affected by an adjustment it
          must be set to true.

    @return 0 on success, 1 on failure

    */

    virtual int adjust_multiplier(bool& values_may_have_changed) = 0;


    /** @brief called by start_sumaugmodel in order to suggest (or to force) *this to opt in or out of a common SumBundle handled by the SumBundleHandler

     The routine decides on contributing to the parents sumbundle or
     not (unless force==true). Its outcome may be, that
     *this does not contribute or that it contributes just some part of
     its subgradient information. As the parent can only discern
     between contributing and not contributing, it is the task of *this
     to feed the parent's Bundlehandler *bh with the correct
     information.  In particular, if the model is willing to contribute
     at all, it should always update its own SumBundle by the same
     rules in order to allow for switching back and forth between
     common and local model without siginificant loss of information.

     The parent's bundlehandler and its corresponding aft is only
     passed here and has to be memorized by the local handler if it
     intends to contribute or communicate with the parent about the
     sumbundle at some point in time.

     The basic assumption is that when the calling parent passes the
     next (the same or another) handler bh (with corresponding aft)
     while *this is already contributing, this next handler and aft is
     consistent with the previous information in that the contribution
     of *this is not changed and *this will be able to continue by only
     updating its contribution.  If this is not the case, the caller
     has to first force this to not contribute to the old handler
     (i.e. to remove its contribution) and then call the routine again
     with its new handler.

     If the SumBundleHandler bh==NULL, the parent does not maintain a
     sumbundle and the model will be strictly local.

       On input:
       - bh is the SumBundleHandler for updating the parent's sumbundle
         (NULL means there is no parent)
       - If mode==child, the parent wants *this to contribute to its
         sumbundle via the bundle handler bh but it never forces *this to
         contribute. In particular *this may decide to contribute only
         part of its model. Whatever decision *this takes, it has to
   make sure that in all other routines it passes on the corresponding
         information consistently.
       - If mode==inactive and force==false, the parent suggests not to
         contribute (maybe not any more, in this case use the bundlehandler
         to remove the contribution) and *this should use its own
         local model.
       - If mode==inactive and force==true, *this has to remove its
         contribution to summodel and switch back to its own model
         which still may be a SumBundle and
       - If mode==root the caller asks (or forces) *this to use
         its own sumbundle in root mode and neither to use its local
         model nor to contribute to the parent.
       - if no bundle parameters are passed (bp==0) default choices are used.
       - aft tells which affine function transformation to use when
         contributing to the parent model. If it is 0, no transformation
         is required

       On output:
       - if call_my_local_model==true the submodel keeps some
         local model part that is not contributed (it may have a contribution
   as well), so it must be called when generating the quadratic
         subproblem; This information will only be needed in start_augmodel()
         and get_row() and will be updated in the beginning of
   the parent's start_augmodel()
       - if call_my_local_model==false all information is incorporated
         in the summodel; this may only occur if bh!=0 and mode and force
         allow this.

       It returns 0 if things worked out as they should, 1 otherwise.

   */
    virtual int sumbundle_mode(SumBundle::Mode& mode,
      SumBundleHandler* bh = 0,
      AffineFunctionTransformation* aft = 0) = 0;

    ///see BundleModel::start_augmodel() for the first four parameters; for the others see sumbundle_mode()
    virtual int start_sumaugmodel(QPModelDataPointer& blockp,
      CH_Matrix_Classes::Integer cand_id,
      const CH_Matrix_Classes::Matrix& cand_y,
      const CH_Matrix_Classes::Indexmatrix* indices = 0,
      SumBundleHandler* bh = 0,
      SumBundle::Mode mode = SumBundle::inactive,
      AffineFunctionTransformation* aft = 0) = 0;


    /** @brief generate the next cutting model and store the center information in the case of a descent step. This version also allows for updating a common model owned by parent SumModel.

      If model_update==null_step , the next model has to contain at least
      all convex combinations of the current subgradient of eval_function()
      and the aggregate minorant of (re)eval_augmodel(). update_model may
      only be called if at least one of both  is available.

      The point passed in the argument is the candidate where the last
      function evaluation took place, its id is passed on so that
      some consistency check is possible without having to store the
      entire vector. The coordinates of the point may help in setting
      up a model of good quality close to this point

      If model_upate is descent_step and the last function
      evaluation took place for this candidate, the function evaluation
      data and point identifier (not nec. the point itself) has to be
      stored as the function value of the center.

      If model_update is new_subgradient the model may wish to use new
      subgradient information available in the candidate for inclusion
      in the model; if the modle is empty it has to include it so that
      the modle is initialized. The center is not moved in this case.

      Whenever the model is already initialized before the call to this
      routine, the output variable model_deviation returns the
      difference between the value of the local aggregate and the true
      function in cand_y, otherwise it is set to zero.

      If there is a parent SumBundle, it has already been updated
      except that -- if *this is indeed contributing -- it needs the
      information about the new subgradient from *this and the
      bundlehandler passed in sumbundle_contribution() knows where
      to add this.
      The update of the common and local sumbundle now requires the
      following steps:

      If model_update is new_subgradient the model may wish to use new
      subgradient information available in the candidate for inclusion
      in the model; if the modle is empty it has to include it so that
      the modle is initialized. The center is not moved in this case.

      The default implementation sets default choices for
      model_deviation and model_curvature and then calls the routine
      update_model(ModelUpdate,CH_Matrix_Classes::Integer,const CH_Matrix_Classes::Matrix&,CH_Matrix_Classes::Real,CH_Matrix_Classes::Real).

      @param[in] model_update
          - if modle_update==new_subgradient, then the latest function evaluation is
            due to a separate funtion evaluation; include the new subgradient
            but do not expect any informtion on a new aggregate (this should always
            be called after setting a new center, in particular on initialization
            so that the model contains at least one subgradient).
          - if model_update==descent step, then the next step is a descent step;
      in this case the model may be restarted from scratch and there is no
      need not include the last aggregate; move the center to the candidate
    - if mode_update==null_step, then the next tep is a null step;
            in this case the update has to preserve the aggregate.


      @param[in] center_id
          the point id of the current center

      @param[in] center_y
          the coordinates of the current center,
          its availability may help to compute the relevance of minorants

      @param[in] cand_id
          the point id of the candidate, so that the model can check
    consistency of the evaluation data in moving the center for
          a descent step

      @param[in] cand_y
          the next model should be good close to this point

      @param[in] model_maxviol
          a minorant violated by this would have resulted in a null step

      @param[in] H
          the weight intended for the proximal term of the next model

      @param[out] model_deviation
          the difference of function value and aggregate in cand_y

      @param[out] model_curvature
          an estimate of the curvature when violating minorants is allowed by model_maxviol

      @return
        - 0 ... if the requirements on the model could be met and no errors occured
        - 1 ... otherwise
    */

    virtual int update_model(ModelUpdate model_update,
      CH_Matrix_Classes::Integer center_id,
      const CH_Matrix_Classes::Matrix& center_y,
      CH_Matrix_Classes::Integer cand_id,
      const CH_Matrix_Classes::Matrix& cand_y,
      CH_Matrix_Classes::Real model_maxviol,
      BundleProxObject& H,
      CH_Matrix_Classes::Real& model_deviation,
      CH_Matrix_Classes::Real& model_curvature) = 0;


    //@}

    //----------------------------------------------------------------------
    /** @name messages for direct get/set requestss */
    //@{

    /// returns @a aftmodel
    AFTModel* get_aftmodel() {
      return aftmodel;
    }

    /// returns @a aftmodel (as a const variant)
    const AFTModel* get_aftmodel() const {
      return aftmodel;
    }

    /** @brief returns the essential information on the bundle date for continuing from the same point lateron, see BundleData */
    virtual BundleData* get_data() = 0;

    /** @brief const version returning the essential information on the bundle date for continuing from the same point lateron, see BundleData */
    virtual const BundleData* get_data() const = 0;

    /** @brief reinstalls the bundle date from get_data() for continuing from this previous point, see BundleData */
    virtual int set_data(BundleData*) {
      return 1;
    }

    /** @brief if primal data is provided by the oracle then the primal
        corresponding to the current aggregate is formed in @a primal
        (this may differ if @a primal is a recognized derived class)
    */
    virtual const PrimalData* get_approximate_primal() const {
      return get_data()->get_approximate_primal();
    }

    /** @brief if primal data is provided by the oracle then the primal
        corresponding to the best eps-subgradient of the evaluation in
        the current center is returned (this may differ if @a primal is a
        recognized derived class)
     */
    virtual const PrimalData* get_center_primal() const {
      return get_data()->get_center_primal();
    }

    /** @brief if primal data is provided by the oracle then the primal
        corresponding to the best eps-subgradient of the evaluation in
        the current center is returned (this may differ if @a primal is a
        recognized derived class)
     */
    virtual const PrimalData* get_candidate_primal() const {
      return get_data()->get_candidate_primal();
    }

    /** @brief if the function is the Lagrangian dual and primal_data of
       previous calls has now to be updated due to changes in the primal
       problem -- e.g., this may happen in column generation -- the
       problem updates all its internally stored primal_data objects by
       calling PrimalExtender::extend on each of these. returns 0 on
       success, 1 if not applicable to this function, 2 if it would be
       applicable but there is no primal data.
     */
    virtual int call_primal_extender(PrimalExtender& pext) {
      return get_data()->call_primal_extender(pext);
    }

    /** @brief set max_bundle_size and max_model_size (this may
        differ if the parameter is a recognized derived class); model
        blocks without bundle return 1.
    */
    virtual int set_bundle_parameters(const BundleParameters&) {
      return 1;
    }

    /** @brief returns the current parameter settings; model blocks
        without bundle parameters return NULL.
    */
    virtual const BundleParameters* get_bundle_parameters() const {
      return 0;
    }

    /** @brief set max_bundle_size and max_model_size (this may
        differ if the parameter is a derived class, in particular
        a SumBunldeParametersObject); model
        blocks without bundle return 1.
    */
    virtual int set_sumbundle_parameters(const BundleParameters&);

    /** @brief returns the current parameter settings; model blocks
        without bundle parameters return NULL.
    */
    virtual const SumBundleParametersObject* get_sumbundle_parameters() const {
      return sumbundle_parameters;
    }

    /** @brief modifications of this specific problem were such that old subgradient data and function values have to be removed completely; do so. In the special case of some particular support functions it may be possile to regenerate some of the minorants and keep the core of the model of the support set; in this case set discard_minorants_only=true.
     */
    virtual void clear_model(bool discard_minorants_only = false) {
      get_data()->clear_model(discard_minorants_only);
    }

    /** @brief remove all current aggregate cutting planes, but not the other parts of the model; returns 0 on success, 1 on failure */
    virtual void clear_aggregates() {
      get_data()->clear_aggregates();
    }

    /** @brief for functions given by an oracle: return the return value that
        was produced by the last call to the function evaluation routine;
        model blocks without oracle return 0.
    */
    virtual int get_ret_code() const {
      return 0;
    }


    //@}

    //---------------------------------------------------------------------- 
    /**@name output and statistics on solution times */
    //@{

    /// return time spent in total for evaluating the model in eval_model()
    CH_Tools::Microseconds get_evalmodel_time() const {
      return evalmodel_time;
    }
    /// return time spent in total for updating the model in update_model()
    CH_Tools::Microseconds get_updatemodel_time() const {
      return updatemodel_time;
    }
    /// return time spent in total for the oracle in eval_function()
    virtual CH_Tools::Microseconds get_preeval_time() const {
      return preeval_time;
    }
    /// return time spent in total for the oracle in eval_function()
    virtual CH_Tools::Microseconds get_eval_time() const {
      return eval_time;
    }
    /// return time spent in total for the oracle in eval_function()
    virtual CH_Tools::Microseconds get_posteval_time() const {
      return posteval_time;
    }
    /// return time spent in total for the add_variable_metric routine
    CH_Tools::Microseconds get_metric_time() const {
      return metric_time;
    }

    ///output the timing statistics 
    std::ostream& print_statistics(std::ostream& out) const;

    /// set output and outputlevel of warnings and errors recursively, see CBout
    void set_out(std::ostream* o = 0, int pril = 1) {
      CBout::set_out(o, pril);
      get_data()->set_cbout(this);
      if (vm_selection)
        vm_selection->set_cbout(this, 0);
    }

    /// set output and outputlevel of warnings and errors recursively with CBout
    void set_cbout(const CBout* cb, int incr) {
      set_out(cb->get_out_ptr(), cb->get_print_level() + incr);
    }

    //@}
  };

  //@}

}

#endif

