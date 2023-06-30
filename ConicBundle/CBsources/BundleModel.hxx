/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/BundleModel.hxx
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



#ifndef CONICBUNDLE_BUNDLEMODEL_HXX
#define CONICBUNDLE_BUNDLEMODEL_HXX


/**  @file BundleModel.hxx
    @brief Header declaring the class ConicBundle::BundleModel
    @version 1.0
    @date 2014-07-14
    @author Christoph Helmberg
*/


#include "QPModelDataObject.hxx"
#include "BundleProxObject.hxx"
#include "GroundsetModification.hxx"
#include "FunctionObjectModification.hxx"

namespace ConicBundle {

  /** @defgroup InternalBundleModel Classes for General and Specialized Cutting Models for Various Objective Functions

      The minimal requirements on a cutting model for BundleSolver are
      laid down in the abstract base class BundleModel. Derived from
      this is the abstract base class SumBlockModel. The latter provides
      the actual foundation for all implemented cutting models, it also
      feeds the common quadratic solver with model data for the quadratic
      bundle subproblem. Any implemented SumBlockModel can either be used
      directly as a separate single model for BundleSolver or it can be
      added to a SumModel. A SumModel represents and manages the sum
      over several models of convex functions (each must be an
      implemented SumBlockModel).  For each implemented SumBlockModel
      the actual model data is kept in a class (derived from)
      BundleData. For the following types of convex functions there
      exist the following pairs of model and data

      - NNCModel and NNCData (a ConeModel for non negative cones; it implements a polyhedral model
        for a general MatrixFunctionOracle / FunctionOracle)

      - SOCModel and SOCData (a ConeModel for the second order cone; it implements a second order cone model  for a SOCOracle)

      - PSCModel and PSCData (a ConeModel for the positve semidefinite cone; it implements a spectral bundle model
        for a maximum eigenvalue oracle PSCOracle)

      - BoxModel and BoxData (a ConeModel implementing a special purpose polyhedral model
        for a support function over a box given via a BoxOracle)

      - SumModel and BundleData (implementing the sum over instances of SumBlockModel)

      - AFTModel and AFTData (implementing an Affine Function or an
        Affine Function Transformation of some associated SumBlockModel)



      Each SumBlockModel may be thought of as offering the choice between two models.
      In the first the function is represented by its actual local model implemented
      in this particulay model class. In the second the local minorant information is
      used to build a general polyhedral SumBundle that may itself either be used
      as a local model or may contribute to a parent SumBundle, so that
      the minorant information of several functions is collected and added up in
      one common polyhedral model. The SumBundle is organized by a SumBundleHandler
      that is set in SumBlockModel::start_sumaugmodel() after calling SumBlockModel::sumbundle_mode().

  */
  //@{

    /** @brief abstract interface for BundleSolver giving access to all objective function specific bundle routines and model descriptions. In particular it hides the cutting model and the oracle.

  Viewed from the BundleSolver, the basic data maintained by the BundleModel are (see also ConicBundle::BundleData)
     - the function value in the center with the relative
       precision required at its computation (but it need
       NOT store the center or its subgradient)
     - the function value in the last candidate together with
       its subgradient and the required relative precision
     - the linear aggregate minorant to the function as collected
       over time.
     - a nonnegative counter modification_id that is increased by the model
       whenever  changes in the function might invalidate the computed
       function value in the center
     - a nonnegative counter aggregate_id that is increased by
       the model whenever changes in the function or model might
       invalidate the aggregate minorant

  Note, it is not assumed that the model keeps track
  of the current center point, candidate point or the ground set;
  if such data is relevant, it is given by arguments.

  The most important functionality required is
     - to evaluate the function at a given point returning an
       upper bound on it in required precision (eval_function())
     - to evaluate the model in a given point (eval_model())
     - to pass the model description for the bundle subproblem
       as a QPModelDataObject via QPModelDataPointer (start_augmodel())
     - to construct the next modelaggregate from the solution information
       in the QPModelDataObject after the bundle subproblem was solved
       (make_model_aggregate())
     - to return the most recent model aggregate
       (get_model_aggregate())
     - to update the model and center data according to the last candidate
       computations and the null/descent step decision (update_model())

  Some convex functions offer the possibility of providing
  reasonable variable metric information for the proximal term. Limited
  support for adapting the proximal term by a variable metric heuristic
  is provided by the routine add_variable_metric().


  ConicBundle is designed to handle problem changes on the fly, so a
  number of routines are devoted to recognizing and synchronizing
  external changes in the groundset data, the model or the function
  itself.  In particular the bundle solver keeps track of the counter
  value of @a center_id that it receives along with data associated with
  the center and of the counter value of @a aggregate_id that it
  receives along with data associated with the model aggregate. When
  these values turn out to be smaller than the ones stored in the model,
  then the data needs to be synchronized before going on. The following
  routines serve this purpose:

      - center_modified() checks whether center_id was increased
      - recompute_center() reevaluates function in the center if this needed
      - model_aggregate_modified() checks wether aggregate_id was increased
      - provide_model_aggregate() provides a valid model aggregate, if the
        current one is no longer valid
      - apply_modification() passes on information on problem changes to the model
      - synchronize_ids() allows to reset all ids to a common value from
        outside (might be needed for branch and bound at some point in time)

  The most typical implementation mistake when using ConicBundle is
  that the user returns incorrect function values and subgradients.
  In order to support recognizing this at least a little bit, it is
  possible to turn on a sanity check based on testing for consistency
  namely  whether the new subgradient inequality computed in the candidate
  holds for the upper bound computed in the current center.
  This is done in the routine check_center_validity_by_candidate().
  If this fails, the center value might be small or the new
  subgradient might be wrong.
   */

  class BundleModel :public VariableMetricModel {

  public:

    BundleModel(CBout* cb = 0, int cbinc = -1); ///< constructor (cb allows to set output options)
    virtual ~BundleModel(); ///< virtual destructor

    /// for informing update_model() at what stage it is called to update the bundle so that the amount of information available is clear
    enum ModelUpdate {
      new_subgradient, ///< the latest function evaluation and its subgradient arise from an extra evaluation of the function
      descent_step, ///< the latest function evaluation and its subgradient give rise to a descent step, preserving the aggregate is not that important
      null_step ///< the latest function evaluation and its subgradient resulted in a null step, preserve the aggregate!
    };


    /** @brief evaluates the objective function in @a y and returns an upper bound in @a ub within relative precision @a relprec.

      If evaluated by an iterative method that provides upper (@a ub) and lower bounds (lb),
      the method may stop when the lower bound is above the @a nullstep_bound.
      Evalutation may also stop, if (ub-lb)<@a relprec*(|ub|+1) is satisfied.

      @param[out] ub_fid
          gives the modification id of the function for which the value @a ub
    was computed

      @param[out] ub
          (relprec-tight) upper bound on the function value in y

      @param[in] y_id
          the identification number of the point to evaluate at
          (to ensure consistency in routines like recompute_center())

      @param[in] y
          the point to evaluate at.

      @param[in] nullstep_bound (CH_Matrix_Classes::Real)
          if a lower bound on the function value by a minorant
          is above this value, a null step will be made
          and any minorant above this bound suffices to ensure
          convergence

      @param[in] relprec (CH_Matrix_Classes::Real)
          if the nullstep_bound is not reached, then evalutation may
          stop, if (ub-lb)<relprec*(|ub|+1) is satisfied for an upper bound
          ub and a lower bound lb on the objective

      @return
          - 0 ... if all is ok
          - >0 ... if solution could not be computed due to fatal errors
          - <0 ... if solution could not be computed to desired precision
    */
    virtual int eval_function(CH_Matrix_Classes::Integer& ub_fid,
      CH_Matrix_Classes::Real& ub,
      CH_Matrix_Classes::Integer y_id,
      const CH_Matrix_Classes::Matrix& y,
      CH_Matrix_Classes::Real nullstep_bound,
      CH_Matrix_Classes::Real relprec) = 0;


    /** @brief evaluate the current cutting model in the given point

      If evaluated by an iterative method that provides upper (ub) and
      lower bounds (@a lb) it may stop when  (ub-lb)<@a relprec*(|ub|+1) is satisfied.

      The result @a lb is only used immediately and not used again later,
      therefore there is no need for passing the modification id of the
      function.

      @param[out] lb
          (relprec-tight) upper bound on the model value in y

      @param[in] y_id
          the point id of the point to evaluate at

      @param[in] y
          the point to evaluate at

      @param[in] relprec
          in an iterative method evalutation may stop, if
          (ub-lb)<relprec*(|ub|+1) is satisfied for an upper bound ub
          and a lower bound lb on the model value

      @return
          - 0 ... if all is ok
          - >0 ... if solution could not be computed due to fatal errors
          - <0 ... if solution could not be computed to desired precision
    */

    virtual int eval_model(CH_Matrix_Classes::Real& lb,
      CH_Matrix_Classes::Integer y_id,
      const CH_Matrix_Classes::Matrix& y,
      CH_Matrix_Classes::Real relprec) = 0;




    /** @brief the returned @a blockp points to a description of the
        variables and constraints generating the cutting model as
        required by the current QP Sollver that provides@a bolckp, see
        QPModelDataObject (the object stays property of this), the
        constant_minorant and the bundle hold the corresponding model
        data


        If indices != NULL, only the coefficients specified by indices will
        be retrieved from constant_minorant and bundle.


      @param[inout] blockp
          pointer for generating a bundle description object suitable for the current QP solver

      @param[in] cand_id
          the point id of the latest candidate

      @param[in] cand_y
          the coordinates of the latest candidate

      @param[in] indices
          if not NULL the subgradient coordinates use in the bundle subproblem
          will just consider these indices, the other coordinates will not be
    used

      @return
        - 0 on succes,
        - !=0 if the necessary information for forming the aggregate is not avilable

    */

    virtual int start_augmodel(QPModelDataPointer& blockp,
      CH_Matrix_Classes::Integer cand_id,
      const CH_Matrix_Classes::Matrix& cand_y,
      const CH_Matrix_Classes::Indexmatrix* indices = 0) = 0;


    /** @brief after the common QP is solved, this call asks to form the new
        aggregate from the solution. If keep_penalty_fixed==false the model may
        decide to increase some internal penalty parameter, has then to
        report this in penalty_parameter_increased but need not form the aggregate.

        This model has access to the solution information directly via the
        QPModelDataObject, of which a pointer has been generated/provided
        in start_augmodel().

        If the model's function contains a penalty component for some constraint
        not represented in Groundset, the model may find out by inspecting the
        solution that this penalty would need to be increased. Because
        increasing a penalty parameter changes the cost function on the fly,
        this requires a lot of attention in forming the bundle subproblem and
        the BundleSolver. Indeed, it might lead to infinite loops in the case of
        infeasible problems.  Therefore it is only allowed to do so if @a
        keep_penalty_fixed == false and it has then to set @a
        penalty_parameter_increased to true. Increasing the parameter will
        always lead to a call to recompute_center(), where the model has to ensure
        that the value with the new penalty is computed correctly. Furthermore,
        the quadratic bundle subproblem has to be resolved in this case.

        @param[out] penalty_parameter_increased
            true if increased (only allowed if keep_penalty_fixed==false), false otherwise

        @param[in] keep_penalty_fixed
            if true, a possibly present penalty parameter may not be changed

        @return
            - 0 on succes,
            - !=0 if the necessary information for forming the aggregate is not avilable

    */
    virtual int make_model_aggregate(bool& penalty_parameter_increased, bool keep_penalty_fixed) = 0;


    /** @brief returns the model aggregate if available.

      The model aggregate is usually the best linear minorant
      \f$(\sigma,s)\f$=(@a model_aggregate_offset,@a model_aggregate)
      of the model of the objective function \f$f\f$
      computed in successful calls to BundleSolver::eval_augmodel(),
      BundleSolver::reeval_augmodel(), see there for its
      precise definition as one part of the saddle point solution.
      It is typically not availabel initially or after problem
      modifications. The nonnegative counter @a model_aggregate_id
      serves to identify the validity of the aggregate via
      the routine model_aggregate_modified(). If no aggregate
      is available, calling the routine provide_model_aggregate()
      will generate some reasonable replacement aggregate with
      a larger aggregate_id that can then be retrieved with the
      present routine.

    @param[out] model_aggregate_id
        nonnegative counter for quickly checking validity of the aggregate

    @param[in,out] model_aggregate
        If empty on input, model_aggregate is initialized to the aggregate,
        otherwise the aggregate is added to model_aggregates.

    @return
      - 0 if the information is available
          (if eval_augmodel() did not return 0, the information will not
          satisfy the precision requirements but may still be available)
      - 1 if the desired information is not available
    */

    virtual int get_model_aggregate(CH_Matrix_Classes::Integer& model_aggregate_id,
      MinorantPointer& model_aggregate) = 0;


    /** @brief generate the next cutting model and store the center information in the case of a descent step

      If model_update==null_step, the next model has to contain at least
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
      data and point identifier (not nec. the point itself) needs to be
      stored as the new function value of the center.

      If model_update is new_subgradient the model may wish to use new
      subgradient information available in the candidate for inclusion
      in the model; if the modle is empty it has to include it so that
      the modle is initialized. The center is not moved in this case.

      @param[in] model_update
          - if model_update==new_subgradient, then the latest function evaluation is
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
          the point id of the center, so that the model can check
    consistency of the evaluation data in moving the center for
          a descent step

      @param[in] center_y
          this may help to estimate the model changes relative to this point

      @param[in] cand_id
          the point id of the candidate, so that the model can check
    consistency of the evaluation data in moving the center for
          a descent step

      @param[in] cand_y
          the next model should be good close to this point

      @param[in] model_maxviol
          a minorant violated by this ammount would have resulted in a null step

      @param[in] H
          the weight intended for the proximal term of the next model.
          it should not be changed here, but may use computation routines
          of H

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
      BundleProxObject& H) = 0;


    /** @brief reset all id counters and references to zero and
        center_id and cand_id to the given values if consistent data
        is available

        After the completion of a successful bundle iteration (null or
        descent) all relevant evaluation and bundle data must be
        available and up to date, even if typically the candidate data is
        not needed any more. This function is intended for  being
        called at exactly this point in the algorithm. If data is indeed
        consistent with this, it resets all center and candidate ids to
        the predifined ids, all other current ids are passed on and if
        the corresponding results are still valid, the respective
        ids should also be set to zero (the other ids reflect "the"
        aggregate, the current version of a specific function or
        transformation, so at most one valid object of the
        other ids should be around and that will get the new version
        number 0). If not valid any more, -1 should be passed back
        instead. The purpose is to make generating consistent
        restarting information easier (this might be useful e.g. in
        conjunction with branch and bound or when linking in new
        functions at some point in time).

        In particular, if nonnegative ids are given for center or candidate ids
        and their respective function value on input, a new nonnegative function
        id for this value is returned if all relevant data is available and
        consistent, otherwise the returned function id will be -1 to indicate
        that the old value needs to be recomputed.  A negative function id on
        input will always result in a -1 being returned. If the target id of the
        center or candidate point is nonnegative, the stored ids must match the
        old versions or the stored versions MUST be discarded and this will
        automatically result in also returning -1 for the repsective function
        ids.  If candidate and center are different they may not be assigned the
        same new id (checked by assert). Likewise, if they are identical they
        may not be assigned distinct nonnegative numbers (one negative is
        allowed, again checked by assert).

        @param[in,out] new_center_ub_fid
            - on input: the function id the value ub was computed for
                 (in this case all point and evaluation data is expected
           to be consistent for the center)
           or -1 if it is known to be out of date
      - on output: the sychronized new function id (expected to be 0) if
           the value is still valid, -1 otherwise

        @param[in] new_center_id
            the newly assigned id for the center point or -1 if the
      model may discard its center information. new_center_id
      may only be >=0 if also old_center_id>=0.

        @param[in] old_center_id
            the previously assigned id for the center point or -1 if
      none is available. If old_center_id>=0 and new_center_id>=0
      the old_center_id must match the center id stored in the model.

        @param[in,out] new_cand_ub_fid
            same as for center

        @param[in] new_cand_id
            same as for center

        @param[in] old_cand_id
            same as for center

        @param[in,out] new_aggregate_id
            - on input: the aggregate id that was returned by get_model_aggregate
           or -1 if it is known to be invalid
      - on output: the sychronized aggregate id (should be 0) if
           the aggregate is still valid, -1 otherwise

    */

    virtual int synchronize_ids(CH_Matrix_Classes::Integer& new_center_ub_fid,
      CH_Matrix_Classes::Integer new_center_id,
      CH_Matrix_Classes::Integer old_center_id,
      CH_Matrix_Classes::Integer& new_cand_ub_fid,
      CH_Matrix_Classes::Integer new_cand_id,
      CH_Matrix_Classes::Integer old_cand_id,
      CH_Matrix_Classes::Integer& new_aggregate_id) = 0;

    /** @brief returns true if the evaluation data for the known
        @a function_modification_id and for the identifier @a center_id
        for the center point is no longer valid or available at the model.

        This usually means that some function or model modifications made it
        impossible to maintain or ensure the validity of the previously computed
        data (the subgradient in the center has no relevance here but would
        then also need to be removed).

        If the model contains further models and relies on their evaluations,
        this call has to be propagated recursively to these and all
        intermediate results that rely on wrong computations have to be marked
        as no longer valid.

        If there were changes but the model implementation is aware that
        these changes do not invalidate the center data, then it may
        also return a new @a function_modification_id and still return false.

        If it returned true, call recompute_center to update the center data.

        @param[in,out] function_modification_id
            - on input: the function id that was returned by routine
                  computing the value in the current center
                        (may be -1 if not initialized)
      - on output: the current function id (always >=0)

        @param[in] center_id
            the point id of the center point (always >=0)

        @return
            - true if the center value should be recomputed
                   (maybe due to recent changes in the function)
            - false if the old evaluation is still valid
    */

    virtual bool center_modified(CH_Matrix_Classes::Integer& function_modification_id,
      CH_Matrix_Classes::Integer center_id) = 0;


    /** @brief after modifications of the problem the center information may have  to be recomputed partially or completely

      The BundleSolver calls this routine for initialization, after problem
      modifications and if numerical problems indicate that the value
      center_ub is actually too small. Some models might call it themselves
      if, e.g., updates of a penalty parameter during (re)eval_augmodle()
      require the recomputation of the function value in the center.

      In the cases where the center was already available before (@a
      center_id matches the one stored) an actual recomputation might be
      avoidable if the problem modifications did not affect the function
      and relative precision requirements did not increase (@a relprec
      negative is used to indicate this). Yet even in some of these
      cases higher precision requirements (for this @a relprec needs to
      be compared to the relative precision of the last center
      computaton, which must be stored by the class) or the requirement
      to produce a value at least as high as the old one (@a
      accept_only_higher_values is set to true in this case and new)
      make a reevalution necessary anyways.

      @param[out] new_center_ub_fid
    the function id for which the value of
          @a new_center_ub was now computed

      @param[out] new_center_ub
          upper bound on the function value in center_y, hopefully in desired precision

      @param[in] center_id
          the point id of the center

      @param[in] center_y
          the coordinate vector of the center

      @param[in] accept_only_higher_values
          if false any result of the recomputation is stored as new center
          information, if true (in this case center_id must match the
          stored number) only the higher result of the previous (the
    value must still be available there) and the current
    computation is stored. It will only be recomputed, however,
    if the relative precision requirements have increased so
    that there is hope for a different result in deterministic
    computations.

      @param[in] relprec
          negative values indicate that the same relative precision is
          used as in the stored center computation (which must include
          center_ub and its "personal" relprec) and may only occur if
          the center is still the same (center_id matches the stored
          value). If relprec has smaller value than in the stored
          computational results, then the function is reevaluated even
          if the center did not change.

      @return
          - 0 ... if all is ok
          - >0 ... if solution could not be computed due to fatal errors
          - <0 ... if solution could not be computed to desired precision

    */
    virtual int recompute_center(CH_Matrix_Classes::Integer& new_center_ub_fid,
      CH_Matrix_Classes::Real& new_center_ub,
      CH_Matrix_Classes::Integer center_id,
      const CH_Matrix_Classes::Matrix& center_y,
      bool accept_only_higher_values = false,
      CH_Matrix_Classes::Real relprec = -1.) = 0;


    /** @brief returns true if the data about the aggregate minorant has changed w.r.t. old_model_aggregate_id. In this case call provide_model_aggregate and then add_model_aggregate to update the model_aggregate */

    virtual bool model_aggregate_modified(CH_Matrix_Classes::Integer old_model_aggregate_id) = 0;



    /** @brief makes sure that the model_aggregate returned by add_model_aggregate
        is actually a minorant contained in the next cutting model.

        This should only be called if model_aggregate_modified()
        returned true (e.g. on intialization). The point @a y is either
        the last candidate or the center and should hint at
        where the model should be good, this should help to choose a good
        initial aggregate if several possibilities exist.
    */

    virtual int provide_model_aggregate(CH_Matrix_Classes::Integer y_id,
      const CH_Matrix_Classes::Matrix& y) = 0;



    /** @brief passes on modification information about function and ground set changes

        If there is no center information or the current groundset modifications
        modify the center too much so that there is only hope to preserve some
        model information (the aggregate is of particular importance) but there
        is no hope to keep other previous evaluation information, then the
        routine is called with new_center_id=old_center_id=-1 and
        empty matrices in new_center and old_center will do.

        If the current groundset modifications modify the center just a
        little or in a very canocial way so that there is a chance to
        preserve previous modle information and function evaluation
        results, then provide the old @a old_center_id, the @a old_center
        coordinates, and the current (new) @a center_id and @a center_y
        coordinates as they arise from the ConicBundle::GroundsetModification
        @a gsmdf (if the ids match, they must be identical). For this class
        the information how to modify itself is given by the
        ConicBundle::FunctionObjectModification @a fundmdfmap[@a
        oracle], where @a oracle is the pointer to the function object
        (function oracle) that gives rise to this model.  If no such
        Modification is contained in the map, then there is no need to
        adapt to anything but the groundset changes.  The class has to
        decide by itself which model components and function evaluation
        results can be preserved under these modificiations.

        Depending on the modification, all function evaluation results may well
        be deleted or marked as outdated. The routine will try to preserve the
        model but in general this will not be possible. The output variable
        no_change is set to true if the modifications had no effect on the
        validitiy of any data returned so far; if no_change==false, the old
        data is no longer valid and needs recomputation.
    */

    virtual int apply_modification(bool& no_changes,
      const GroundsetModification& gsmdf,
      const FunObjModMap& funmdfmap,
      CH_Matrix_Classes::Integer new_center_id,
      const CH_Matrix_Classes::Matrix& new_center,
      CH_Matrix_Classes::Integer old_center_id,
      const CH_Matrix_Classes::Matrix& old_center) = 0;


    /** @brief consistency check for oracle computations: test if the subgradient inequality arising out of the last eval_function holds for @a center_y.

    @param[out] cand_minorant_is_below is set to true if subgradient inequality holds, otherwise false.

    @param[in] center_id
        point id of the center, it should match the one stored here

    @param[in] center_y
        the coordinates of the current center (need not be available inside this class)

    @return 0 if all required data is available, !=0 if something failed

    */

    virtual int check_center_validity_by_candidate(bool& cand_minorant_is_below,
      CH_Matrix_Classes::Integer center_id,
      const CH_Matrix_Classes::Matrix& center_y) = 0;

    /** @brief returns the minorant corresponding to the subgradient inequality returned by the last function evaluation. If the minorant is not empty on input, the local minorant is added to it. If no subgradient is available or modifications invalidated it, the minorant will be set to empty.
     */
    virtual int get_function_minorant(CH_Matrix_Classes::Integer& function_modification_id, MinorantPointer& minorant) = 0;

    /** @brief returns the minorant corresponding to the subgradient inequality returned by the function evaluation for the current center. If the minorant is not empty on input, the local minorant is added to it. If no subgradient is available or modifications invalidated it, the minorant will be set to empty.
     */
    virtual int get_center_minorant(CH_Matrix_Classes::Integer& function_modification_id, MinorantPointer& minorant) = 0;


    /** @brief Overload this in order apply transformations in between.

    The bundle solver never calls the above routines directly but always
    via a call to transform. This allows to insert affine function
    transformations, if so desired. */
    virtual BundleModel* transform() {
      return this;
    }

    /// replaces variable_metric_transform by transform
    VariableMetricModel* variable_metric_transform() {
      return transform();
    }

  };



  //@}

}

#endif

