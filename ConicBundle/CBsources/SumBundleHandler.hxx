/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/SumBundleHandler.hxx
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



#ifndef CONICBUNDLE_SUMBUNDLEHANDLER_HXX
#define CONICBUNDLE_SUMBUNDLEHANDLER_HXX


/**  @file SumBundleHandler.hxx
    @brief Header declaring the class ConicBundle::SumBundleHandler (see ConicBundle::SumBlockModel)
    @version 1.0
    @date 2014-07-22
    @author Christoph Helmberg
*/


#include "SumBundle.hxx"
//#include "BundleModel.hxx"
#include "SumBundleParametersObject.hxx"

namespace ConicBundle {

/** @ingroup dynamic_submodel_selection
*/
//@{

/** @brief routines for updating and handling SumBundle components, possibly  by cooperating over several recursive levels

 On construction the handler is initialized with the sumbundle it works on 
 with a map that informs the handler on which parts of the sumbundle it has to
 work on, and if so with which function factor. The bundle handler will only
 treat these parts and ignore the ones not listed. If it also gets a set of 
 bundle parameters it will use this same set for all relevant parts.

 The mode stored in the SumBundle (and lateron the information stored in the
 bundle handler) governs what the handler has to do with the sumbundle.  ...
 Whether the bundle handler actually does something with the parts depends on
 the SumBundle mode root, child, inactive, and also on whether there is a
 SumBundleHandler of the parent model that it has to cooperate with.
    
 */


class SumBundleHandler: public CBout
{ 
private:
  /// this is the bundle this handler works on
  SumBundle* sumbundle;

  /// if there is a parent sumbundle (that this might contribute to or not), this is its handler
  SumBundleHandler* parent_handler;
  ///use this transformation when writing back to the parent bundle
  const AffineFunctionTransformation* aft;

  /// allows to store organizational information for each FunctionTask type of bundle
  class BundleInformation
  {
  public:
    SumBundleParametersObject* sbp; ///< holds bundle and model sizes and the model selection routine 
    CH_Matrix_Classes::Integer bundle_size;      ///< current size of the bundle
    CH_Matrix_Classes::Integer aggr_index;       ///< index of the aggregate
    CH_Matrix_Classes::Integer new_index;        ///< index where to write the new subgradient to
    CH_Matrix_Classes::Indexmatrix map_to_old;   ///< if not empty, map_to_old(i) holds the minorant index of i before the most recent update (like in reassign_minorants)
    CH_Matrix_Classes::Indexmatrix model_indices; ///< the indices of the bundle minorants employed in the model, only needed at a root, not propagated recursively
    CH_Matrix_Classes::Matrix old_lowrank;  ///< in dynamic scaling the previous low rank scaling matrix
    CH_Matrix_Classes::Matrix old_diagonal; ///< in dynamic scaling the previous diagonal scaling matrix
    CH_Matrix_Classes::Symmatrix old_sym; ///< in dynamic scaling the previous symmetric scaling matrix

    QPConeModelDataObject* block;                           ///< the feasible set of the local quadratic model
    CH_Matrix_Classes::Real increase_factor;     ///< the factor by which function_factor was increased

    //the following avoid introducing further function argumetns in start_augmodel
    CH_Matrix_Classes::Real tmodel_maxviol; ///< the violation parameter used in the latest call to update_model
    CH_Matrix_Classes::Real tweightu; ///< the weight used in the latest call to update_model 
    BundleModel::ModelUpdate tmodel_update; ///< the ModelUpdate parameter in the latest call to updated_model

    /// default constructor
    BundleInformation();

    /// destructor
    ~BundleInformation();

    ///copies/clones only the information regarding the bundle selection; initializes other members
    BundleInformation& operator=(const BundleInformation& bi);
  };

  
  /// for storing separate data for ObjectiveFunciton, ConstantPenaltyFunction, AdaptivePenaltyFunction
  typedef std::vector<BundleInformation*> BundleInformationVector;
  /// stores separate BundleInformation for ObjectiveFunciton, ConstantPenaltyFunction, AdaptivePenaltyFunction
  BundleInformationVector bundleinfo;

  /// called by the constructor; initializes the bundleinfo depending on whether as sumbundle part exists already or  not 
  int init(FunctionTask ft,CH_Matrix_Classes::Real funfactor,const BundleParameters* bp);

  /// sets max_bundle_size and max_model_size for this FunctionTask; this may be increased internally if the mode and/or the parent handler require this 
  int set_bundle_parameters(FunctionTask ft,const BundleParameters& bp);

  /// store the (normalized) aggregate for FunctionTask ft in its position in the bundle, potentially with forming its primal and update coeff to reflect this (the aggregate is assumed to be available in sumbundle->get_aggregate(ft), but not its primal)
  int store_aggregate(FunctionTask ft); 

  /** @brief afterwars position i is held by the minorant the had position map_to_old(i) before; the current bundle size is reduced to map_to_old.dim();

  This also updates coeff. If a minorant appears multiple times, the
  value of coeff is only copied for the first occurence, later ones get
  value 0.; If a minorant with a positive coeff is not copied, this
  results in an error. Without errors sum(coeff) will keep its value.
  The values of aggr_index and new_index are not changed.
  */
  int reassign_minorants(const CH_Matrix_Classes::Indexmatrix& map_to_old,FunctionTask ft);


  /// the only change to SumBlockModel::make_model_aggregate() is that FunctionTask refers to the part for which the aggregate should be formed
  int make_model_aggregate(FunctionTask ft);

  /// see SumBlockModel::provide_model_aggregate)(, do this for this FunctionTask
  int provide_model_aggregate(FunctionTask ft);

  /** @brief updates the sumbundle even if not active for that FunctionTask; in this it tries to do the same as the parent handler bh. It stores the aggregate at the same place and provides new room at the same place; the new subgradient is *not* entered here but in a subsequent call to set_cand_minorant()

      Whenever the sumbundle has no contributors, update_model just makes
      sure that the update rules of an existing parent are stored so that they
      can be communicated to subbundles if needed. If there is no parent, it
      uses the default rules to update an imaginary empty bundle and provides
      an imaginary location for the new information and by doing so it may
      enlarge the imaginary bundle.

      The only situation in which a sumbundle may have a contributor but the
      bundlesize is zero is when a model is started. In this case the update
      just makes room for a first minorant which must then be added by
      set_cand_minorant();
   */
  int update_model(BundleModel::ModelUpdate model_update,
		   CH_Matrix_Classes::Integer center_id,
		   const CH_Matrix_Classes::Matrix& center_y,
		   CH_Matrix_Classes::Integer cand_id,
		   const CH_Matrix_Classes::Matrix& cand_y,
		   CH_Matrix_Classes::Real model_maxviol,
		   BundleProxObject& H,
		   FunctionTask ft);

  /// evaluate the model value  for this FunctionTask
  int eval_model(CH_Matrix_Classes::Real& lb,CH_Matrix_Classes::Integer yid,const CH_Matrix_Classes::Matrix& y,FunctionTask ft) const;

  /// returns a *quick* lower bound for the model value for this FunctionTask
  CH_Matrix_Classes::Real lb_model(CH_Matrix_Classes::Integer yid,const CH_Matrix_Classes::Matrix& y,FunctionTask ft) const;
  
  /// once all new minorants have been collected at this level, they are passed to the next if required for the FunctionTask
  int contribute_new_minorants(FunctionTask ft);

  ///remove valid contributions and set the states correspondingly for the FunctionTasks
  int remove_contributions(FunctionTask ft);

  ///add valid contributions and set the states correspondingly for the FunctionTask
  int add_contributions(FunctionTask ft);

/// align the bundle to the parent handler bh for this FunctionTask so that contributing to it is possible
  int align_bundle(bool use_parent_update,FunctionTask ft);

 /// align the bundle to the parent handler bh so that conributing to it is possible
  int align_bundle(bool use_parent_update);

  /// first calls remove_contributions(), then discards all current models
  void clear_model(FunctionTask ft);

  /// first calls remove_contributions(), then discard all aggregate minorants in the current model
  void clear_aggregates(FunctionTask ft);


public:

  ///initialize the handler to handle sb and therein the parts for which factor_map sets the (initial) function factors; all of them get the same bundle parameters if specified; otherwise default values are used
  SumBundleHandler(SumBundle& sb,const std::map<FunctionTask,CH_Matrix_Classes::Real>& factor_map,const BundleParameters* bp=0);
  ///
  ~SumBundleHandler();

  /// returns the sumbundle 
  const SumBundle* get_sumbundle() const {return sumbundle;}

  /// returns true if this FunctionTask is handled
  bool handles(FunctionTask ft)
  {return (sumbundle->has_bundle_for(ft));}

  /// returns the index of the newest subgradient in the bundle 
  CH_Matrix_Classes::Integer get_new_index(FunctionTask ft) const
  {const BundleInformation* bi=bundleinfo[ft]; return (bi==0)?-1:bi->new_index;}

  /// returns true if the corresponding part is root with contributions but has bundle_size 0
  bool initialization_needed(FunctionTask ft) const;

  /// returns true if one of the parts is root with contributions but has bundle_size 0
  bool initialization_needed() const;

  /** @brief sets a parent handler, an aft and prepares for the next in_mode. If the respective pointers are null and in_mode is active, a part in child mode is changed to root.

  If the parent pointer is not null, the routine calls align_bundle in
  order to match the existing bundle to the parent's.  Therefore, if
  the parent pointer is not null and the parent has a positive bundle
  size, set_parent_information may only be called if all parts that
  have contributors also have a bundle containing an aggregate,
  i.e. they need a bundlesize of at least one.

  If in_mode==root, the handler has to perpare the sumbundle to serve as root.

  If in_mode==inactive, the handler has to remove any relevant
  contributions and to deactivate the sumbundle

  If in_mode==child, the action depends on the the current mode. If
  the mode is child or root already, no changes are required. In
  particular, if the mode is root and the in_mode child is desired,
  this needs to be achieved later by a call to add_contributions()
  (once all bundle components have been collected). If the current
  mode is inactive, the sumbundle's mode is changed to root and again
  an add_contributions() will be needed to change that to child.

  */

  int set_parent_information(SumBundleHandler *parent_sbh,const AffineFunctionTransformation* aft,SumBundle::Mode in_mode); 

  /// resets the value of the function factor for this part of sumbundle
  int reset_function_factor(FunctionTask ft,CH_Matrix_Classes::Real factor);

  /// sets max_bundle_size and max_model_size for all parts; this may be increased internally if the mode and/or the parent handler require this 
  int set_bundle_parameters(const BundleParameters& bp);


  /// returns the increase factor for the unbounded part of sumbundle, otherwise 1.
  CH_Matrix_Classes::Real get_increase_factor() const;

  /// updates the sumbundle maybe even if not active; in this it tries to do the same as bh, storing the aggregate at the same place and providing new room at the same place; the new subgradient is *not* entered here but in set_cand_minorant()
  int update_model(BundleModel::ModelUpdate model_update,
		   CH_Matrix_Classes::Integer center_id,
		   const CH_Matrix_Classes::Matrix& center_y,
		   CH_Matrix_Classes::Integer cand_id,
		   const CH_Matrix_Classes::Matrix& cand_y,
		   CH_Matrix_Classes::Real model_maxviol,
		   BundleProxObject& H);


  /// evaluate the model value 
  int eval_model(CH_Matrix_Classes::Real& lb,CH_Matrix_Classes::Integer yid,const CH_Matrix_Classes::Matrix& y) const;

  /// returns a *quick* lower bound for the model value
  CH_Matrix_Classes::Real lb_model(CH_Matrix_Classes::Integer yid,const CH_Matrix_Classes::Matrix& y) const;
  

  /** brief bring the bundle into a normal form so that contributions may be added and subtracted consistently

  The normal form has the aggregate w.r.t. to the last bundle cofficients
  stored in column aggr->index in scaled form so that in coeff all weight is
  set to it. No other columns are modified.

  When add_contributions is called, it is assumed that parent bundle and
  contributing bundle are in this form. Some care has to be taken that this
  normalization is happening for parent and children in coordinated form so
  that remove_contribution() does not cause havoc. For this, the sumbundle
  should always be normalized in (or right before) calling
  SumBlockModel::sumbundle_contribution() before starting any interaction with
  the parent or the children.
  */
  int normalize_sumbundle();

  /// once all new minorants have been collected at this level, they are passed to the next if required
  int contribute_new_minorants();

  ///remove own contributions to the parent and set the states correspondingly
  int remove_contributions();

  /** @brief add root parts of this sumbundle to the parent's sumbundle and set the states correspondingly 

      If the current mode is child, nothing is done, because it is
      assumed that the sumbundle is already part of the parent. Thus,
      only root parts are considerd new and are added. This is only
      done, if the parenhandler accepts these kind of parts.

      When a contribution is added to a parent that is currently
      inactive or child, any respective parent's contributions to the
      parent's parent are first removed and then the mode of the
      parent is set to root. Thus, the parent has to call
      add_contributions afterwards, if it still wants to contribute to
      its own parent.

      add_contribution should only be called from
      SumBlockModel::sumbundle_contribution() with a normalized bundle,
      i.e., normalize_sumbundle() should have been called before. 
   */
  int add_contributions();

  /// start the augmented model block for FunctionTask ft if to be handled here and increase xdim correspondingly; if there are several, use start_augmodel(QP_SumBlock&,Integer&) instead
  int start_augmodel(QPModelDataPointer& bp,
		     CH_Matrix_Classes::Integer cand_id,
		     const CH_Matrix_Classes::Matrix& cand_y,
		     const CH_Matrix_Classes::Indexmatrix* indices,
		     FunctionTask ft);

  /// add augmented model blocks to the sumblock for parts to be handled here and increase xdim correspondingly
  int start_augmodel(QPModelDataPointer& bp,
		     QPSumModelDataObject& sumblock,
		     CH_Matrix_Classes::Integer cand_id,
		     const CH_Matrix_Classes::Matrix& cand_y,
		     const CH_Matrix_Classes::Indexmatrix* indices=0);


  /// see SumBlockModel::make_model_aggregate
  int make_model_aggregate(bool& increased, bool fixed);
  
  /// see SumBlockModel::provide_model_aggregate
  int provide_model_aggregate();
  
  /// see SumBlockModel::adjust_multiplier
  int adjust_multiplier(bool& values_may_have_changed);
  
  /** @brief (re)initialize the bundle of the respective part
 
     This is only possible, if the handler handles a bundle for this.
     If not this causes an error. The main purpose is really 
     to start the bundle on the fly if a parent starts a sumbundle
     in the middle of the computation. The information in coeff
     is assumed to yield the current aggregate.
  
     If this bundle part has been contributed berfore (its mode is child), the
     contribution is first removed from the parent. Then any existing
     information except for the function factor is discarded and replaced by
     the new bundle information with the number of contributions set to 1.

     If the size of primals is not zero, it must have one nonzero entry per 
     column of minorants and all future calls via set_cand_minorant() also have
     to provide exactly one primal for each update.
 
  */
  int contribute_initial_bundle(FunctionTask ft,
				const MinorantBundle& bundle_minorants,
				const CH_Matrix_Classes::Matrix& coeff);
   
  /** @brief replace the aggregate by one from outside
 
     This is only possible, if the handler handles a bundle for this and 
     the bundle is not active. The main purpose is to switch from an
     external model to a newly contributing sumbundle that has been
     updated all along but not been in use. Installing the external
     aggregate then ensures convergence.  
  
     If primal is not zero, it must already have one nonzero entry per 
     column of existing minorants and all future calls via set_cand_minorant() 
     also have to provide exactly one primal for each update.
 
  */
  int install_external_aggregate(FunctionTask ft,
				 const MinorantPointer& aggr,
				 CH_Matrix_Classes::Real aggr_coeff);
   

  /** @brief set the new minorant information of the candidate
  
     Only one minorant can take this position. 
     This is mainly due to that only one primal information can be associated
     with each minorant. Thus adding must involve compatible primals and
     this is easiest to guarantee if there is only one inital type at the
     lowest level. Contributions to parents will not care about the primals.
   
  */ 
  int set_cand_minorant(FunctionTask ft,
		       const MinorantPointer& minorant);

  /// first calls remove_contributions(), then discards all current models
  void clear_model();

  /// first calls remove_contributions(), then discards all aggregate minorants in the current model
  void clear_aggregates();

  /// before contributing the new evaluation results, the candidates have to be cleared
  void clear_cand_minorants();

  /// first calls remove_contributions(), then appends the data to the bundle of minorants
  int append_vars(FunctionTask ft,const CH_Matrix_Classes::Matrix& append_mat);

  /// first calls remove_contributions(), then reassigns the rows of the bundle of minorants according to map_to_old
  int reassign_vars(FunctionTask ft,const CH_Matrix_Classes::Indexmatrix& map_to_old);

  /// see DynamicScaling   
  int add_variable_metric(FunctionTask ft,
			  VariableMetric& H,
			  CH_Matrix_Classes::Integer yid,
			  const CH_Matrix_Classes::Matrix& y,
			  bool descent_step,
			  CH_Matrix_Classes::Real weightu,
			  CH_Matrix_Classes::Real model_maxviol,
			  const CH_Matrix_Classes::Indexmatrix* indices=0);

  /// see DynamicScaling   
  int add_variable_metric(VariableMetric& H,
			  CH_Matrix_Classes::Integer yid,
			  const CH_Matrix_Classes::Matrix& center_y,
			  bool descent_step,
			  CH_Matrix_Classes::Real weightu,
			  CH_Matrix_Classes::Real model_maxviol,
			  const CH_Matrix_Classes::Indexmatrix* indices=0);


  /// computes an estimate of the current curvature to support SumModel in the selection of submodles 
  CH_Matrix_Classes::Real guess_curvature(const MinorantBundle& mnrts,
					  const CH_Matrix_Classes::Indexmatrix& selected_indices,
					  CH_Matrix_Classes::Integer cand_id,
					  const CH_Matrix_Classes::Matrix& cand_y,
					  CH_Matrix_Classes::Real model_maxviol) const;

};




  //@}

}

#endif

