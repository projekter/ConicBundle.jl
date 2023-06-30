/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/VariableMetricSVDSelection.hxx
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



#ifndef CONICBUNDLE_VARIABLEMETRICSVDSELECTION_HXX
#define CONICBUNDLE_VARIABLEMETRICSVDSELECTION_HXX


/**  @file VariableMetricSVDSelection.hxx
    @brief Header declaring the classes ConicBundle::VariableMetricSVDSelection
    @version 1.0
    @date 2019-03-12
    @author Christoph Helmberg
*/


#include "VariableMetric.hxx"

namespace ConicBundle {


  /** @ingroup InternalVariableMetric   */


//@{


  /** @brief general implementation of a VariableMetricSelection routine to form and add variable metric information to a BundleProxObject generically mainly from the aggregate and a collection of minorants given by a MinorantBundle

   the class currently also includes a number of experimental routines 
   for computing suitable low rank metric information from the latest minorant data
 */

class VariableMetricSVDSelection: public virtual VariableMetricSelection
{
private: 
  /// a positive number means dynamic scaling should be used and suggests to store and use this number of the most recent minorants for this; if <=0, no metric information will be added
  CH_Matrix_Classes::Integer n_latest_minorants;
  
  /// specifies which of the current experimental routines should be used for computing the metric
  CH_Matrix_Classes::Integer selection_method;

  /// for values of oldfactor between 0. and 1. the SVD is taken of (1-oldfactor)*vecH*vecH'+oldfactor*oldvecH*oldvcH'; for values <=0 there is no old contribution, for values >=1 oldvecH is reused
  CH_Matrix_Classes::Real oldfactor;
  
  /// experimental data that tries to collect the current Vspace by an orthonormal basis (only some routines attempt this)
  CH_Matrix_Classes::Matrix Vspace;

  /// experimental data that tries to collect the current Uspace information by an orthonormal basis  (only some routines attempt this)
  CH_Matrix_Classes::Matrix Uspace;

  /// experimental data that tries to identify the curvature (Hessian eigenvalues) in Uspace directions (only some routines attempt this)
  CH_Matrix_Classes::Matrix lambda;

  /// experimental candidate that results from an approximate Newton step computaiton based on Vspace, Uspace and lambda (only some routines compute this)
  CH_Matrix_Classes::Matrix candNewton;
    

  /// scaling heuristic based on collected minorants relative to the aggregate alone 
  int vecH_by_aggregate(CH_Matrix_Classes::Matrix& vecH,
			CH_Matrix_Classes::Real& relprec,
			const MinorantBundle& bundle,
			const MinorantPointer& aggregate,
			CH_Matrix_Classes::Real function_factor,
			CH_Matrix_Classes::Integer y_id,
			const CH_Matrix_Classes::Matrix& y,
			CH_Matrix_Classes::Real weightu,
			CH_Matrix_Classes::Real violation_eps);

  /// scaling heuristic based on collected minorants relative to aggregate and model
  int vecH_by_model(CH_Matrix_Classes::Matrix& vecH,
		    CH_Matrix_Classes::Real& relprec,
		    const MinorantBundle& bundle,
		    const MinorantBundle& model,
		    const MinorantPointer& aggregate,
		    CH_Matrix_Classes::Real function_factor,
		    CH_Matrix_Classes::Integer y_id,
		    const CH_Matrix_Classes::Matrix& y,
		    CH_Matrix_Classes::Real weightu,
		    CH_Matrix_Classes::Real violation_eps);

  /// scaling heuristic based on collected minorants relative to aggregate and orthogonalized to model
  int vecH_orthogonal_to_model(CH_Matrix_Classes::Matrix& vecH,
			       CH_Matrix_Classes::Real& relprec,
			       const MinorantBundle& bundle,
			       const MinorantBundle& model,
			       const CH_Matrix_Classes::Matrix& modelcoeff,
			       const MinorantPointer& aggregate,
			       CH_Matrix_Classes::Real function_factor,
			       CH_Matrix_Classes::Integer y_id,
			       const CH_Matrix_Classes::Matrix& y,
			       CH_Matrix_Classes::Real weightu,
			       CH_Matrix_Classes::Real violation_eps);

  /// scaling heuristic based on collected minorants relative to aggregate and orthogonalized to model
  int vecH_orthogonal_to_modelSVD(CH_Matrix_Classes::Matrix& vecH,
				  CH_Matrix_Classes::Real& relprec,
				  const MinorantBundle& bundle,
				  const MinorantBundle& model,
				  const CH_Matrix_Classes::Matrix& modelcoeff,
				  const MinorantPointer& aggregate,
				  CH_Matrix_Classes::Real function_factor,
				  CH_Matrix_Classes::Integer y_id,
				  const CH_Matrix_Classes::Matrix& y,
				  CH_Matrix_Classes::Real weightu,
				  CH_Matrix_Classes::Real violation_eps);

  /// scaling heuristic based weighted sum of collected minorants relative to active subgradients and orthogonalized to model
  int vecH_weighted_SVD(CH_Matrix_Classes::Matrix& vecH,
				  CH_Matrix_Classes::Real& relprec,
				  const MinorantBundle& bundle,
				  const MinorantBundle& model,
				  const CH_Matrix_Classes::Matrix& modelcoeff,
				  const MinorantPointer& aggregate,
				  CH_Matrix_Classes::Real function_factor,
				  CH_Matrix_Classes::Integer y_id,
				  const CH_Matrix_Classes::Matrix& y,
				  CH_Matrix_Classes::Real weightu,
				  CH_Matrix_Classes::Real violation_eps);

  /// scaling heuristic based weighted sum of SVD Hessians collected via minorants relative to active subgradients 
  int vecH_weighted_SVDs(CH_Matrix_Classes::Matrix& vecH,
				  CH_Matrix_Classes::Real& relprec,
				  const MinorantBundle& bundle,
				  const MinorantBundle& model,
				  const CH_Matrix_Classes::Matrix& modelcoeff,
				  const MinorantPointer& aggregate,
				  CH_Matrix_Classes::Real function_factor,
				  CH_Matrix_Classes::Integer y_id,
				  const CH_Matrix_Classes::Matrix& y,
				  CH_Matrix_Classes::Real weightu,
				  CH_Matrix_Classes::Real violation_eps);

public:
  /// destructor
  ~VariableMetricSVDSelection()
  {}

  /// default constructor
  VariableMetricSVDSelection(CBout* cb=0,int cbincr=-1):
    VariableMetricSelection(cb,cbincr),
    n_latest_minorants(50),
    selection_method(-1),
    oldfactor(0.)
  {}

  ///  constructor for specifying values for n_latest_minorants and selection_method
  VariableMetricSVDSelection(CH_Matrix_Classes::Integer in_n_latest_minorants,
			     CH_Matrix_Classes::Integer in_selection_method,
			     CH_Matrix_Classes::Real in_oldfactor=0.,
			     CBout* cb=0,int cbincr=-1):
    VariableMetricSelection(cb,cbincr),
    n_latest_minorants(in_n_latest_minorants),
    selection_method(in_selection_method),
    oldfactor(in_oldfactor)
  {}
  
  ///for current ongoing experiments with variable metric routines
  void get_UVlambda(CH_Matrix_Classes::Matrix& U, 
		    CH_Matrix_Classes::Matrix& V,
		    CH_Matrix_Classes::Matrix& lam,
		    CH_Matrix_Classes::Matrix& cand)
  {U=Uspace;V=Vspace;lam=lambda;cand=candNewton;}

  ///returns n_latest_minorants
  CH_Matrix_Classes::Integer get_n_latest_minorants() const
  {return n_latest_minorants;}

  ///sets n_latest_minorants
  void set_n_latest_minorants(CH_Matrix_Classes::Integer nlm)
  {n_latest_minorants=nlm;}

  ///returns selection_method
  CH_Matrix_Classes::Integer get_selection_method() const
  {return selection_method;}

  ///sets selection_method
  void set_selection_method(CH_Matrix_Classes::Integer sm) 
  {selection_method=sm;}

  ///see ConicBundle::VariableMetricSelection::add_variable_metric()
  int add_variable_metric(VariableMetric& H,
			  CH_Matrix_Classes::Integer y_id ,
			  const CH_Matrix_Classes::Matrix& y,
			  bool descent_step,
			  CH_Matrix_Classes::Real weightu,
			  CH_Matrix_Classes::Real model_maxviol,
			  const CH_Matrix_Classes::Indexmatrix* indices,
			  VariableMetricBundleData& bundle_data);

  /// clone: the values are only preserved for those contained in the constructor: n_latest_minorants, selection_method and oldfactor
  VariableMetricSelection* clone_VariableMetricSelection() {
    return new VariableMetricSVDSelection(n_latest_minorants,selection_method,oldfactor,this,0);
  }
};



  //@}

}

#endif

