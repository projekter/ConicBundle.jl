/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/VariableMetric.hxx
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



#ifndef CONICBUNDLE_VARIABLEMETRIC_HXX
#define CONICBUNDLE_VARIABLEMETRIC_HXX


/**  @file VariableMetric.hxx
    @brief Header declaring the classes ConicBundle::VariableMetricModel, ConicBundle::VariableMetricBundleData, ConicBundle::VariableMetricSelection, and  ConicBundle::VariableMetric
    @version 1.0
    @date 2019-03-12
    @author Christoph Helmberg
*/


#include "CBout.hxx"
#include "MinorantPointer.hxx"

namespace ConicBundle {

  class AffineFunctionTransformation;

/** @defgroup InternalVariableMetric   Variable Metric 

   @brief The classes and routines here help to dynamically specify the part
   of the quadratic term within the proximal term (see \ref
   InternalBundleProxObject) that is responsible for the shape of the
   level sets of the prox term in the bundle subproblem. In the smooth
   case this quadratic contribution should be the Hessian of the
   function. Generally the term should make it easy to move into
   directions that are still considered promising and where
   "curvature" also seems to allow for larger steps.

   The class VariableMetric gives the abstract interface to classes
   describing the quadratic term. VariableMetric implementations may
   offer one or several possible ways to add second order information
   (diagonal, low rank, or dense) for cutting models. To provide this,
   the latter tpyically use collected subgradient or curvature
   information depending on the oracle/model in use. An istance of the 
   class VariableMetricSelection may be passed to VariableMetric via
   the constructor or VariableMetric::set_variable_metric_selection(),
   that provides a default routine for constructing and adding such 
   information from default data provided by a call of the model to 
   VariableMetricSelection::add_variable_metric() with the
   data requested in VariableMetricBundleData . VariableMetric 
   specifies several further interface routines that are needed
   for the BundleSolver and the BundleModel to ensure treatment consistent
   with the use of AffineFunctionTransformation, SumBundle etc. or
   options on which models should contribute at what stage to the 
   VariableMetric.

   For a model to provide customized "second order" information it
   must be a VariableMetricModel and the model has to implement the
   virtual routine VariableMetricModel::add_variable_metric() (all
   SumBlockModel implementations in ConicBundle have this).  This
   routine may query VariableMetric on which matrix structures are
   supported (diagonal, low rank, or dense) and the model may use a
   corresponding VariableMetric::add_variable_metric() routine to add
   its second order data. Potentional transformations due to a
   possible cascade of AffineFunctionTransformation applications will
   be performed automatically when used within the standard framework;
   for this the VariableMetric class is informed about such
   transformations recursively by VariableMetric::push_aft() and
   VariableMetric::pop_aft(). Alternatively the model may be equipped
   with a separate implemented version of a VariableMetricSelection class 
   that does so. The VariableMetricBundleData passed on should then give
   access to the respective BundleData of the model via a suitable
   dynamic cast.


   A VariableMetric implementation that supports a dynamic metric 
   allows to switch on or off various scaling components:
   
   - Whenever a VariableMetricSelection object is supllied,  this 
     switches on the variable metirc routines. How many of the latest
     gradients are used in the construction of the metric information
     by the selection routine is passed on in the second argument when 
     it calls the routine VariableMetricBundleData::get_latest_minorants().
     This number may be regarded as a suggestion to store and provide
     that many of the most recent minorants. 
     If the number is <=0, no information is requested from the models 
     and the initial scaling matrix is used as the basis.

   - local_scaling only has an affect if a VariableMetricSelection
     object is specified in VariableMetric. If local_scaling==true,
     every model that contributes its own local model to the bundle
     subproblem should add its own local scaling contribution; in
     contrast, local_scaling==false asks to only use the global
     minorants at the top level without any local resolution to be
     used for forming the scaling matrix, so the first model called
     should provide all the relevant information without recursive
     calls to subsequent models.

   - bounds_scaling is independent of the other two and refers to
     a heuristic within LPGroundset for adding diagonal weights if
     the variables are box constrained and close to their bounds.
     If the current step would cause the variable to exceed the
     bound, the diagonal elements of the scaling matrix are increased 
     so that when considering the diagonal of the scaling matrix only 
     the variables just end up right on the boundary. This heuristic
     can only be expected to work well with diagonal variants of
     VariableMetric.
   
*/

//@{

  class VariableMetric;

/// declares the interface that a BundelModel needs to provide for contributing to VariableMetric information

class VariableMetricModel: public virtual CBout
{
public:
  /// constructor for passing on ouptut information
  VariableMetricModel(CBout* cb=0,int cbincr=-1);
  /// virtual destructor
  virtual ~VariableMetricModel();

  /** @brief  add to the variable metric information H some model dependent "second order" information of the function modelled

    The routine may add nothing if there is not sufficient information.
    Otherwise see the functions in VariableMetric
    for how to enter information depending on how H supports dynamic scaling.
   
    @param[in,out] H
      BundleScaling H holds the current scaling matrix, provides
      the information needed to add to it.

    @param[in] y_id
      the point id of the point supplied

    @param[in] y
      the scaling should be relevant around this point, 
      e.g. in the smooth case the (nonnegative) diagonal or low rank representation of its Hessian 

    @param[in] descent_step
      true if and only if the solver just performed a descent step; if not true, the dynamic scaling routine may only add further semidefintie matrices but may not make the term smaller in the semidefinite sense.
     
    @param[in] weightu
      the positive scalar value currently used for weighting the quadratic (trust region) term 
     
    @param[in] model_maxviol
      for nonsmooth functions the aggregate together with the quadratic term should aim at keeping the value above the value of the minorant minus this model_maxviol (where all minorants are assumed or shifted to cut below the aggregate by at least .1 this value)

    @param[in] indices
       if not NULL, the effect of the changes of the quadratic term must be
       restricted to these indices
     
    @return
      - 0 ... on success or if not supported (= adds nothing), 
      - 1 ... on failure,

  */

  virtual int add_variable_metric(VariableMetric& H,
				  CH_Matrix_Classes::Integer y_id ,
				  const CH_Matrix_Classes::Matrix& y,
				  bool descent_step ,
				  CH_Matrix_Classes::Real weightu ,
				  CH_Matrix_Classes::Real model_maxviol,
				  const CH_Matrix_Classes::Indexmatrix* indices =0);

  /** @brief Overload this in order apply transformations in between. 

  VariableMetric never calls the above routines directly but always 
  via a call to transform. This allows to insert affine function 
  transformations, if so desired. */

  virtual VariableMetricModel* variable_metric_transform() {return this;} 
  
};


/** @brief abstract interface providing the bundle data that is typically needed in 
VariableMetricSelection classes.

  The standard models all employ derived classes of this class so as to
  allow customized implementations of VariableMetricSelection classes access
  to the full bundle information via a dynamic_cast.

*/

class VariableMetricBundleData
{
private:
  CH_Matrix_Classes::Matrix vmbd_lowrankH; ///< used by the variable metric low rank heuristic (without function_factor)
  CH_Matrix_Classes::Matrix vmbd_diagH; ///< used by the variable metric low rank heuristic (without function_factor)
  CH_Matrix_Classes::Symmatrix vmbd_denseH; ///< used by the variable metric heuristic 

public:
  /// virtual destructor
  virtual ~VariableMetricBundleData();

  /// the factor by which the minorants (except for the aggregate) need to be multiplied in order to match the current function scaling
  virtual CH_Matrix_Classes::Real get_function_factor() const=0;

  /// latest_minorants should return up to the max_number latest minorants whose information should be incorporated in the metric model; the number may fall short of max_number; if the number returned is too large, only the first max_number will usually be considered; the minorants still need to be mutliplied by function_factor
  virtual int get_latest_minorants(MinorantBundle& latest_minorants,
				   CH_Matrix_Classes::Integer max_number) =0;

  /// model_mionrants should hold the minorants currently used in the model; the list may be empty, it may and typically will contain other minorants than returned in get_latest_minorants(); the minorants still need to be mutliplied by function_factor; modell_coeff needs to return a column vector with the same number of elements as in model_minorants and in the same sorting; mostly it will hold the convex combination of the model_minorants that gives rist to the aggregate, but it may also contain other nonnegativer numbers indicating the importance of the respective minorants 
  virtual int get_model_data(MinorantBundle& model_minorants,
			     CH_Matrix_Classes::Matrix& model_coeff) const=0;

  /// the aggregate minorant as currently in use in the bundle method; it already includes the function_factor
  virtual const MinorantPointer& get_aggregate() const=0;

  /// allows to retrieve dense variable metric information stored here
  virtual const CH_Matrix_Classes::Symmatrix& get_denseH() const
  {return vmbd_denseH;}

  /// allows to retrieve the dense variable metric information generated in the previous call and  allows to store the new one in the end 
  virtual CH_Matrix_Classes::Symmatrix& set_denseH()
  {return vmbd_denseH;}

  /// allows to retrieve low rank variable metric information stored here
  virtual const CH_Matrix_Classes::Matrix& get_lowrankH() const
  {return vmbd_lowrankH;}

  /// allows to retrieve the low rank variable metric information generated in the previous call and  allows to store the new one in the end 
  virtual CH_Matrix_Classes::Matrix& set_lowrankH()
  {return vmbd_lowrankH;}

  /// allows to retrieve the diagonal variable metric information generated in the previous call and  allows to store the new one in the end 
  virtual const CH_Matrix_Classes::Matrix& get_diagH() const
  {return vmbd_diagH;}

  /// allows to retrieve the diagonal variable metric information generated in the previous call and  allows to store the new one in the end 
  virtual CH_Matrix_Classes::Matrix& set_diagH()
  {return vmbd_diagH;}
};

  
/** @brief abstract interface, that allows to specify a routine for providing or computing a suitable variable metric for specific models. Per current default, no such information is generated.
*/

class VariableMetricSelection: public virtual CBout
{
public:
  /// constructor for passing on ouptut information
  VariableMetricSelection(CBout* cb=0,int cbincr=-1);
  /// virtual destructor
  virtual ~VariableMetricSelection(); 

  /** @brief this routine allows to supply a customized scaling heuristic to a VariableMetricModel

    The routine may add nothing if there is not sufficient information.
    Otherwise see the functions in VariableMetric
    for how to enter information depending on how H supports dynamic scaling.
   
    @param[in,out] H
      VariableMetric H holds the current matrix and provides
      the information needed to add to it.

    @param[in] y_id
      the point id of the point supplied

    @param[in] y
      the scaling should be relevant around this point, 
      e.g. in the smooth case the (nonnegative) diagonal or low rank representation of its Hessian 

    @param[in] descent_step
      true if and only if the solver just performed a descent step; if not true, the dynamic scaling routine may only add further semidefintie matrices but may not make the term smaller in the semidefinite sense.
     
    @param[in] weightu
      the positive scalar value currently used for weighting the quadratic (trust region) term 
     
    @param[in] model_maxviol
      for nonsmooth functions the aggregate together with the quadratic term should aim at keeping the value above the value of the minorant minus this model_maxviol (where all minorants are assumed or shifted to cut below the aggregate by at least .1 this value)

    @param[in] indices
       if not NULL, the effect of the changes of the quadratic term must be
       restricted to these indices
     
    @param[in] bundle_data
       the bundle information on the model for which the variable metric should
       be added. The 
     
    @return
      - 0 ... on success or if not supported (= adds nothing), 
      - 1 ... on failure,

   */
  virtual int add_variable_metric(VariableMetric& H,
				  CH_Matrix_Classes::Integer y_id ,
				  const CH_Matrix_Classes::Matrix& y,
				  bool descent_step,
				  CH_Matrix_Classes::Real weightu,
				  CH_Matrix_Classes::Real model_maxviol,
				  const CH_Matrix_Classes::Indexmatrix* indices,
				  VariableMetricBundleData& bundle_data)=0;

  ///clone
  virtual VariableMetricSelection* clone_VariableMetricSelection() =0;
				  
};








  /** @brief interface class that allows a VariableMetricModel to contribute information to a VariableMetric object for use in a BundleProxObject

 */

class VariableMetric: public virtual CBout
{
private: 
  /// NULL also signals not to call dynamic VariableMetric routines; if a corresponding object has been transferred to here, variable metric should be used. The object serves as a default routine for contributing variable metric information derived from the latest minorants
  VariableMetricSelection* vm_selection;
  
  /// flag for whether the scaling information should be collected from active local models; if false, only the root model may provide the scaling information for all
  bool use_local_metric;
  
  /// flag for whether LPGroundset may use its diagonal update heuristic for bounds
  bool use_bounds_scaling;

public:
  /// virtual destructor
  virtual ~VariableMetric();

  /// default constructor; if vp is not zero, ownership of *vp is passed over to *this and *this will delete vp 
  VariableMetric(VariableMetricSelection* vp=0,
		 bool use_loc_metric=false,
		 bool use_bnds_scaling=false,
		 CBout* cbo=0,int cbinc=-1):
    CBout(cbo,cbinc),
    vm_selection(vp),
    use_local_metric(use_loc_metric),
    use_bounds_scaling(use_bnds_scaling)
  {
    vm_selection=vp;
  }
  
  ///returns 0 or an available VariableMetricSelection object, that may be employed for computing a variable metric term in accordance with the value of use_local_metric
  VariableMetricSelection* get_variable_metric_selection() const {return vm_selection;}

  ///sets use_variable_metric; the object passed is then owned by this
  void set_variable_metric_selection(VariableMetricSelection* vp=0)
  {delete vm_selection;vm_selection=vp;}

  /// returns true if add_dense_variable_metric() is supported 
  virtual bool supports_dense_variable_metric() const
  {return false;}

  /// returns true if add_variable_metric() does not ignore the low rank argument vecH
  virtual bool supports_lowrank_variable_metric() const
  {return false;}

  /// returns true if add_variable_metric() does not ignore the diagonal argument diagH
  virtual bool supports_diagonal_variable_metric() const
  {return false;}

  /// returns true if some dynamic scaling is supported and switched on
  bool employ_variable_metric() const
  {return  (vm_selection)&&(supports_diagonal_variable_metric()||
	    supports_lowrank_variable_metric()||
	    supports_dense_variable_metric());}

  /// the BundleSolver starts an update of this by dynamic scaling by calling this in every step; negative parameters give no preferences, but the global aggregate (groundset+model) has to be provided and if descent_step==false the scaling may only increase
  virtual int apply_variable_metric(VariableMetricModel* /* groundset */,
				    VariableMetricModel* /* model */,
				    const CH_Matrix_Classes::Matrix& /* aggregate */,
				    CH_Matrix_Classes::Integer /* y_id */,
				    const CH_Matrix_Classes::Matrix& /* y */,
				    bool /* descent_step */,
				    CH_Matrix_Classes::Real& /* current_weight */,
				    CH_Matrix_Classes::Real /* model_maxviol */,
				    const CH_Matrix_Classes::Indexmatrix* /* new_indices */=0)
  {return 1;}
  
  /** @brief adds a suitable modification of symH (symH may be modified in this) to the scaling matrix H

      The suitable modification depends on the on the affine transformations 
      applied to the model or to the form supported by this class.
      The latter may, e.g., be only a diagonal matrix, a low rank matrix or
      a dense one. 

   */
  virtual int add_variable_metric(CH_Matrix_Classes::Symmatrix& /* symH */)
  {return 1;}

  /** @brief adds (a suitable modification of) Diag(diagH)+vecH*transpose(vecH) to the scaling matrix H (either matrix may be modified in this)

      The suitable modification depends on the on the affine transformations 
      applied to the model or to the form supported by this class.
      The latter may, e.g., be only a diagonal matrix, a low rank matrix or
      a dense one. 

      Any of the two arguments may be an empty matrix (0x0)and this part is
      then ignored.

   */
  virtual int add_variable_metric(CH_Matrix_Classes::Matrix& /* diagH */,
				  CH_Matrix_Classes::Matrix& /* vecH */)
  {return 1;}
  

  /// this AffineFunctionTransformation has to be used before applying Hinv in 
  virtual int push_aft(const AffineFunctionTransformation* /* aft */)
  {return 1;}
  
  /// removes the top most aft (without deleting it!)
  virtual int pop_aft()
  {return 1;}


  ///returns use_bounds_scaling
  bool get_use_bounds_scaling() const {return use_bounds_scaling;}

  ///sets use_bounds_scaling
  void set_use_bounds_scaling(bool bounds_scaling)
  {use_bounds_scaling=bounds_scaling;}

  /** @brief if the respective implementation supports a diagonal bounds scaling heuristic, the following routine has to return true; see also diagonal_scaling_heuristic_update()  
  */
  virtual bool supports_diagonal_bounds_scaling() const {return false;}
  
  /** @brief if the respective implementation supports a diagonal bounds scaling heuristic, the following routine has to return true; see also diagonal_bounds_scaling_update()  
  */
  bool employ_diagonal_bounds_scaling() const 
  {return (use_bounds_scaling && supports_diagonal_bounds_scaling());}
  
  /** @brief if supported, D_update has to contain nonnegative numbers that are permanently added to the diagonal here. 

   It is important to keep track of this change only if afterwards update_QP_costs is called before compute_QP_costs. In this case the nonzero entries in D_update must be a subset of the indices in delta_index
  */
  virtual int diagonal_bounds_scaling_update(const CH_Matrix_Classes::Matrix& /* D_update */)
  {return 1;}

  ///returns use_local_metric
  bool get_use_local_metric() const {return use_local_metric;}

  ///sets use_local_metric
  void set_use_local_metric(bool local_metric)
  {use_local_metric=local_metric;}


};

    

  //@}

}

#endif

