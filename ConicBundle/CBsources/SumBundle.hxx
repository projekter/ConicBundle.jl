/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/SumBundle.hxx
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



#ifndef CONICBUNDLE_SUMBUNDLE_HXX
#define CONICBUNDLE_SUMBUNDLE_HXX


/**  @file SumBundle.hxx
    @brief Header declaring the class ConicBundle::SumBundle (see ConicBundle::SumBlockModel)
    @version 1.0
    @date 2015-11-27
    @author Christoph Helmberg
*/


#include <map>
#include "VariableMetric.hxx"
#include "AffineFunctionTransformation.hxx"

namespace ConicBundle {


  /** @defgroup dynamic_submodel_selection dynamic submodel selection

     Dynamic submodel selection is switched on or off via MatrixCBSolver::set_sumbundle().

     Dynamic submodel selection works for sums of convex functions by
     maintaining a common model part for a dynamic selection of the
     oracles. For each oracle in the sum it can be decided whether it
     has a specialized separate model (called local model) or whether it
     contriubtes to the common model. A possible contributor to a sum
     of convex functions may itself be a sum of convex functions
     resulting in a recursive tree structure with the sums being the
     inner nodes and the leaves being the final oracles representing
     "elementary" functions. A common model may be initiated as the local
     model of any node (inner or leaf) which then acts as the root of this
     particular common model.

     Internally dynamic sumbundle selection builds upon the following
     common structures in SumBlockModel, which is the abstract contributing
     block to the inner SumModel nodes (these organize sums of models of
     separate oracles, they are the inner nodes and typical roots of common
     models):

     - SumBundle: Every implementation of a SumBlockModel
       maintains a SumBundle in its associcated BundleData,
       that allows to organize the minorants and aggregates
       in a unified way

     - SumBundleHandler: Each model acitvely participating in a common
       model  has its own SumBundleHandler that either acts as root
       or communicates with its parent on what has to be done. Its tasks comprise
       + the recursively coordinated update of the models including
         adding and removing contributions to the parent's common model
       + as the root, the forming of the contributing QPBlock of the
         quadratic bundle subproblem
       + the forming of the aggregate on its own level
       + in case of dynamic oracle modifications, the corresponding
         modifications of the minorants

     - SumBlockModel::sumbundle_contribution(): a call to this starts
       the decision process on which models take part in a common model
       and which don't. This routine may be used to ask any model to
       act as a root but it is also used repeatedly by a SumModel root
       in the following way. Whenever SumModel::start_augmodel() is
       called, it uses the local model precisions of its descendants
       to work out suggestions on who should contribute to the common
       model and who might profit from a better local model. It then
       asks all its descendants via sumbundle_contribution() to follow
       this suggestion (which they might refuse to do).

  */
  //@{



  /** @brief class for use with SumBlockModel and BundlData for storing and
      managing a common bundle describing (part of) the model

      Organizing a common bundle for a selection of functions in a sum of
      functions may be split into three mayor parts

      - the data itself (SumBundle) that represents the collection of minorants
        and their primals and may contribute to a parent's SumBundle

      - a SumBundleHandler for the data that updates the minorants according to
        some updating scheme in a consistent way over all participating models.

      - a model deciding part (the base class SumBlockModel), that determines
        which parts go into the common part (SumModel), that provides an
        interface to the quadratic bundle subproblem, and that deals with
        function dependent issues (FunctionModel and others).

      The SumModel plays the role of the initiator to the model deciding
      part if its SumBundle has a part set to a socalled root mode (see
      below).  Right before setting up the next quadratic subproblem, it
      tries to determine dynamically for which of its submodels (each a
      SumBlockModel) it might be more efficient to include them in its
      common SumBundle and suggests this to the submodels. Each submodel
      may still decide on its own whether it is willing to contribute to
      the common SumBundle (maybe also only in part) or not.  Because a
      submodel may choose to contribute or not, each submodel has to
      supply the information that it wants to be taken care of in the
      common model by itself.  Thus the minorants to be used in forming
      the common SumBundle must be fed to the SumBundleHandler by the
      SumBundleHandler of the submodel. Furthermore, if the submodel
      decides to change its state of contributing, it has to communicate
      this information accordingly via its SumBundleHandler and it has
      to be able to extract the necessary information from its global
      part in order to continue without loosing the aggregate. The
      potential splitting of the model into a local and a global part
      also requires some care in computing the model value itself so
      that no parts are counted twice.

      In order to be able to switch between local and global model
      frequently without much loss in quality of the model, each
      submodel keeps updating its own potential contribution to the
      sumbundle by its SumBundleHandler even if it is currently not
      included in the SumBundle of the caller. Whenever it switches to
      the global model, it has to communicate its global model
      contribution's aggregate to its SumBundleHandler and the bundle
      handler has to be able to incorporate that in the common
      SumBundle's aggregate, so the local SumBundleHandler has to know
      the SumBundleHandler of the parent. In order to allow proper
      updating of the aggregate of the global model, the parent's bundle
      handler always has to keep the aggregate explicitly as a separate
      minorant in the SumBundle. The submodel itself requires no
      information on who uses its global model or how it is used, but it
      must always be aware of whether it is used or not.  The submodel
      is informed whether it may choose its global model part to be used
      or not by a caller via calling
      SumBlockModel::sumbundle_contribution() with a SumBundleHandler or
      with a NULL pointer respectively. Between two calls of this
      routine the local SumBundleHandler will always interact with this
      as the parent's SumBundleHandler, so the parent's
      SumBundleHandler may not change or be deleted during two calls.
      If a radical change is desired, this always requires a new call of
      SumBlockModel::sumbundle_contribution().  Whenever the external
      use of the global model is switched off, the submodel is
      responsible for including the global part of the model in its
      local model again. In particular, the submodel may do so by
      extracting the aggregate or by setting up its own handler and the
      interface to the global model part ensures that it will always be
      able to do so. If the function itself changes, the model first
      removes its contribution by calling the parent's bundlhandler and
      changes its data only afterwards, so that the parent's model stays
      correct.

      Because of the three different kinds of FunctionTask and their
      different uses in the quadratic subproblem, the SumBundle of
      SumModel needs to handle three different parts, one for each
      FunctionTask.  In contrast, in FunctionModel only one FunctionTask
      is active each time and so FunctionModel will only make use of one
      of them at each time.  As the contributions of the SumBlockModel
      classes to the SumBundle may change dynamically, each part will
      only be fully active if there are contributions to it by some
      function. Each part can be in three modes:

      - if the SumBundle::Mode is SumBundle::root the SumBlockModel
        owning the SumBundle must include all contributions in its
        own model and add it to the quadratic bundle suproblem.
        The mode may be root even without any contributors; then the
        respective bundle data will not be updated and it only serves
        to keep the SumBundle of its potential but currently inactive
        children alive in that the SumBundleHandler uses it to tell
        these potential children how to update the common bundle. A root
        part can have a parent at the same time, but in this case it is
        not contributing to the parent. The SumBundleHandler should,
        however, align its updating rules with those of the parent,
        so that it can switch from root to child at any time.

      - if the SumBundle::Mode is SumBundle::child the SumBlockModel
        contributes its SumBundle to a parent's SumBundle and has
        to leave this part out of its own model. A part of a
        SumBundle can be a child only if it has contributors,
        otherwise it is inactive or a root. For a child part
        the bundle always has to be updated according to the parent's
        rules.

      - if the SumBundle::Mode is SumBundle::inactive, it should only
        have contributors at the leaves (if at all). The
        SumBundleHandler has to update the bundle at leaves with
        contributors because they will potentially contribute at some
        later time, but the SumBundleHandler has to use inactive parts
        in between without contributors for transporting the information
        about the common update rules. inactive parts are currently not
        used in any model, these parts have to be taken care of elsewhere.

  */

  class SumBundle : public CBout {
  public:
    /// specifies for different parts of the sumbundle whether it is active and who is responsible for handling it
    enum Mode {
      root, ///< the SumBundle/part is active and starts here, it does not contribute to parent
      child, ///< the SumBundle/part is active and contributes to the parent, no control here
      inactive, ///< the SumBundle/part is maybe maintained but currently not in use anywhere
      unavailable ///< the SumBundle/part is not even maintained
    };

    ///
    SumBundle();
    ///copy constructor
    SumBundle(const SumBundle& sb);
    ///
    ~SumBundle();

    /// initialize 
    void init(const SumBundle& sb);

    /// sets the modification_id to id
    void synchronize_ids(CH_Matrix_Classes::Integer new_modification_id,
      CH_Matrix_Classes::Integer new_center_id,
      CH_Matrix_Classes::Integer old_center_id,
      CH_Matrix_Classes::Integer new_cand_id,
      CH_Matrix_Classes::Integer old_cand_id,
      CH_Matrix_Classes::Integer new_prex_id = 0);

    /// returns true if BData exists for this mode
    bool has_bundle_for(FunctionTask ft) const;

    /// returns true if BData exists 
    bool has_bundle_data() const;

    /// if BData exists for this mode it returns its bundle_size (possibly 0), otherwise 0
    int bundle_size(FunctionTask ft) const;

    /// returns true if one of its parts is a root
    bool has_roots() const;

    /// returns true if one of its parts is a root with n_contributors>0
    bool has_working_roots() const;

    /// returns true if one of its parts is not inactive
    bool active() const;

    /// returns true if one of its parts is a child
    bool has_contributions() const;


    // gets the corresponding valid flag (call only if a has_bundle_for(ft)==true)
    //bool get_valid(FunctionTask ft) const;

    /// gets the corresponding mode (call only if a has_bundle_for(ft)==true)
    Mode get_mode(FunctionTask ft) const;

    /// gets the corresponding function factor (call only if a has_bundle_for(ft)==true)
    CH_Matrix_Classes::Real get_function_factor(FunctionTask ft) const;

    /// gets the corresponding n_contributors (call only if a has_bundle_for(ft)==true)
    CH_Matrix_Classes::Integer get_n_contributors(FunctionTask ft) const;

    /// gets the corresponding minorants (call only if a has_bundle_for(ft)==true)
    const MinorantBundle& get_bundle(FunctionTask ft) const;

    /// gets the corresponding aggregation coefficients (call only if a has_bundle_for(ft)==true)
    const CH_Matrix_Classes::Matrix& get_coeff(FunctionTask ft) const;

    /// gets the corresponding aggregate (call only if a has_bundle_for(ft)==true)
    const MinorantPointer& get_aggregate(FunctionTask ft) const;

    /// gets the corresponding candidate minorant (call only if a has_bundle_for(ft)==true)
    const MinorantPointer& get_cand_minorant(FunctionTask ft) const;

    /// evaluate the model values for all root parts, it may only be called if there are such because otherwise lb would need to be -infty 
    int eval_model(CH_Matrix_Classes::Real& lb,
      CH_Matrix_Classes::Integer yid,
      const CH_Matrix_Classes::Matrix& y,
      CH_Matrix_Classes::Real nullstep_bound,
      CH_Matrix_Classes::Real relprec) const;

    /// returns a *quick* lower bound for the root parts, it may only be called if there are such because otherwise lb would need to be -infty 
    CH_Matrix_Classes::Real lb_model(CH_Matrix_Classes::Integer yid,
      const CH_Matrix_Classes::Matrix& y) const;

    /// get the aggregate that is due to root sumbundle parts handled here
    int get_local_model_aggregate(MinorantPointer& aggregate,
      CH_Matrix_Classes::Real factor = 1.,
      const AffineFunctionTransformation* aft = 0) const;

    /// get the aggregate that is due to child parts of the sumbundle, which are contributed to parents
    int get_contributed_model_aggregate(MinorantPointer& aggregate,
      CH_Matrix_Classes::Real factor = 1.,
      const AffineFunctionTransformation* aft = 0) const;

    /// call this primal extender for the primals; if there minorants but no primals or if it fails, return 1
    int call_primal_extender(PrimalExtender& prex,
      CH_Matrix_Classes::Integer prex_id,
      FunctionTask ft);

    /// rearrange/extend the minorants according to the given groundset modifications 
    int apply_modification(const GroundsetModification& gsmdf,
      CH_Matrix_Classes::Integer mod_id,
      MinorantExtender* mex,
      FunctionTask ft);

  private:
    friend class SumBundleHandler;

    /** For every SumBundle part it stores the current activity mode and the validity status of its aggregate.

      The mode may only be changed by the local bundlehandler. The model
      has to tell the bundlehandler on initialization, which parts
      it is allowed to use.

      If the mode is root, the bundle was or will be used directly in the
      quadratic bundle subproblem. If the local handler knows about a
      parent handler, it keeps the parent's bundlesize part arranged in
      the same way. It may use more columns, however, if so desired,
      and forms the aggregate in its own way, even if it has to be stored
      in the same position in update_model with the coefficients being
      function_factor for the aggregate and 0 otherwise in order to allow
      contributing to the parent at later points in time.

      If the mode is child, the bundle is incorporated in the parent's
      summodel, so it may not be used to form a quadratic bundle subproblem
      and the local handler has to follow the update rules of the parent
      in detail.

      If the mode is inactive, it is currently not in use anywhere, but it is
      kept alive and up to date in order to allow for switching the mode quickly
      (this is only necessary if there is a local contributing part). In
      particular the local handler keeps the minorants in the same order
      and the model feeds the new minorant in every iteration. Again the
      purpose is to allow contributing at later points in time with
      consistent information.

      If a sumbundle part is active (root or child), the corresponding valid
      flag is false, if some problem modification may have caused the bundle
      information of this to have changed (the BundleHandler needs to make sure
      that parts contributed to parents are removed before that) to the effect
      that the bundle with its coefficients might no longer be in sync with the
      stored aggregate.

      The sumbundle parts usable by the handler must be initialized by the
      model. If the model has not initialized it, the Bundlehandler is not
      allowed to use it. A part that has no contributors cannot be active.
      The only way to establish that a sumbundle is a

      The coefficients of a part are meaningful only if the bundle is active.

    */

    class BData : public CBout {
    public:
      /// ==root is handled here, ==child is handled at a parent, ==inactive if currently not in use but kept in sync with the parent
      Mode mode;

      ///the number of different subbundles added to this sumbundle part
      CH_Matrix_Classes::Integer n_contributors;

      /// for external use all offsets and minorants have to be mulitplied by this factor
      CH_Matrix_Classes::Real function_factor;

      /// the coefficients for forming the next aggregate generated by the last quadratic bundle subproblem, they sum up to at most the value function_factor
      CH_Matrix_Classes::Matrix coeff;

      /// the minorants collected and not discarded over time, they do NOT include function_factor yet
      MinorantBundle bundle;

      /// This is empty or stores last aggregate computed. If not empty, it includes the function_factor 
      MinorantPointer aggregate;

      /// the minorant in the current candidate is stored and collected here and does NOT include the function_facotr
      MinorantPointer cand_minorant;

      /// default constructor
      BData();
      /// copy constructor (uses clear and init)
      BData(const BData& bd);
      /// destructor
      ~BData();
      /// resets all to initial state
      void clear(CH_Matrix_Classes::Real fun_factor = 1.);
      /// make a copy
      void init(const BData& bd);
      /// sets the modification_id to id
      void synchronize_ids(CH_Matrix_Classes::Integer new_modification_id,
        CH_Matrix_Classes::Integer new_center_id,
        CH_Matrix_Classes::Integer old_center_id,
        CH_Matrix_Classes::Integer new_cand_id,
        CH_Matrix_Classes::Integer old_cand_id,
        CH_Matrix_Classes::Integer new_prex_id = 0);
      /// make a copy (uses init)
      BData& operator=(const BData& bd);

      /// evaluate the model value 
      int eval_model(CH_Matrix_Classes::Real& lb, CH_Matrix_Classes::Integer yid, const CH_Matrix_Classes::Matrix& y) const;

      /// returns a *quick* lower bound for the model value
      CH_Matrix_Classes::Real lb_model(CH_Matrix_Classes::Integer yid, const CH_Matrix_Classes::Matrix& y) const;

      /// get the aggregate multiplied by factor, add it if add==true, and with applying the AffineFunctionTransformation if aft!=0 
      int get_model_aggregate(MinorantPointer& aggregate,
        CH_Matrix_Classes::Real factor = 1.,
        const AffineFunctionTransformation* aft = 0) const;

      /// call this primal extender for the primals; if there are minorants but no primals or if it fails, return the count of fails
      int call_primal_extender(PrimalExtender& prex,
        CH_Matrix_Classes::Integer prex_id);

      /// rearrange/extend the minorants according to the given groundset modifications 
      int apply_modification(const GroundsetModification& gsmdf,
        CH_Matrix_Classes::Integer mod_id,
        MinorantExtender* mex);
    };

    /// for storing separate data for ObjectiveFunciton, ConstantPenaltyFunction, AdaptivePenaltyFunction
    typedef std::vector<BData> BDataVector;

    /// stores separate bundle data for ObjectiveFunciton, ConstantPenaltyFunction, AdaptivePenaltyFunction
    BDataVector bdata;

    /// clear all data and set the mode of all parts to inactive
    void clear();

    /// clear all data and remove this part
    void clear(FunctionTask ft);

    /// copy by calling initialize(const SumBundle&)
    SumBundle& operator=(const SumBundle& sb);

    ///reset an existing or create a new BData entry for ft
    void init(FunctionTask ft, CH_Matrix_Classes::Real fun_factor = 1.);

    /** @brief if mode==root or mode==child and *this has a part marked as child, switch the mode to root; if mode==inactive, deactivate child and root parts
     */
    void take_control(Mode mode);

    // allows to set the corresponding valid flag (call only if a has_bundle_for(ft)==true)
    //bool& set_valid(FunctionTask ft);

    /// allows to set the corresponding mode (call only if a has_bundle_for(ft)==true)
    Mode& set_mode(FunctionTask ft);

    /// allows to set the corresponding function factor (call only if a has_bundle_for(ft)==true)
    CH_Matrix_Classes::Real& set_function_factor(FunctionTask ft);

    /// allows to set the corresponding n_contributors (call only if a has_bundle_for(ft)==true)
    CH_Matrix_Classes::Integer& set_n_contributors(FunctionTask ft);

    /// allows to set the corresponding minorants (call only if a has_bundle_for(ft)==true)
    MinorantBundle& set_bundle(FunctionTask ft);

    /// allows to set the corresponding aggregation coefficients (call only if a has_bundle_for(ft)==true)
    CH_Matrix_Classes::Matrix& set_coeff(FunctionTask ft);

    /// allows to set the corresponding aggregate (call only if a has_bundle_for(ft)==true)
    MinorantPointer& set_aggregate(FunctionTask ft);

    /// allows to set the corresponding aggregate (call only if a has_bundle_for(ft)==true)
    MinorantPointer& set_cand_minorant(FunctionTask ft);

    /// if active, get the approximate primal (if add==false, initialize if possible and  set add=true in this case)
    const PrimalData* get_approximate_primal(FunctionTask ft) const;


  };


  //@}

}

#endif

