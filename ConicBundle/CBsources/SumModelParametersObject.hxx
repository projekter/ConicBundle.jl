/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/SumModelParametersObject.hxx
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



#ifndef CONICBUNDLE_SUMMODELPARAMETERSOBJECT_HXX
#define CONICBUNDLE_SUMMODELPARAMETERSOBJECT_HXX

/**  @file SumModelParametersObject.hxx
    @brief Header declaring the class ConicBundle::SumModelParametersObject
    @version 1.0
    @date 2014-07-24
    @author Christoph Helmberg
*/

#include "SumModel.hxx"
#include "SumBundleParametersObject.hxx"

namespace ConicBundle {
/** @ingroup InternalBundleModel

*/
//@{

/** @brief abstract interface for SumModel for the models selection routine when using SumBundle and SumBundleHandler

    Whenever a joint model should be formed dynamically for some of
    the functions within SumModel in order to achieve a reduced model
    size for those that do not need a detailed description, the 
    routine select_models() provides a suggestion which models should
    be asked to contribute to this joint model. 

    The parameter @a max_local_models gives a preference on the maximum
    number of models that may keep their own separate model, all others
    should be asked to join.

    Except for possibly @a update_rule the other parameters of
    BundleParameters are likely not needed here, because the separate
    SumBundleParametersObject is relevant for the SumBundleHandler,
    but the class helps to provide a uniform interface.
*/

  class SumModelParametersObject: public virtual CBout, public BundleParameters
{
protected:
  //int max_model_size;       // maximum number of minorants to be selected for the cutting model (numbers<=1 for no limit, numbers >=2 impose a strict limit)
  //int max_bundle_size;      // suggested maximum number of latest minorants stored for use in a model, for constructing variable metric information etc. (negative numbers give no preference; the size may be increased internally in case of confliciting requirements, eg. in n_model_size or by variable metric routines) 
  //int update_rule;      // in case several update rules are available

  int max_local_models; ///< when combining models, aim for this number of local models
    
public:
  /// initialize BundleParameters to the given values; if bp is a SumModelParametersObject alos set max_local_models and possibly clone sum_bundle_handler correspondigly; if bp is a SumBundleParametersObject, clone it into sum_bundle_handerl;
  int init(const BundleParameters& bp){
    BundleParameters::init(bp);
    const SumModelParametersObject* smpo=dynamic_cast<const SumModelParametersObject*>(&bp);
    if (smpo){
      max_local_models=smpo->max_local_models;
      return 0;
    }
    return 0;
  }

  /// default constructor
  SumModelParametersObject(CBout* cb=0,int cbinc=-1):
    CBout(cb,cbinc),BundleParameters(),max_local_models(0)
  {}

  ///virtual destructor, implemented in SumModelParameters.cxx
  virtual ~SumModelParametersObject();
  
  /// returns the value of the variable 
  virtual int get_max_local_models() const
  { return max_local_models;}
 
  /// returns the value of the variable 
  virtual int set_max_local_models(int mlm)
  { max_local_models=mlm; return 0;}
 
  /// generate a suggestion list on which to include in the SumBundle (the othersi should keep their local models) 
  virtual int select_models(SumModel::ModelMap& modelmap)=0;

};


  //@}

}

#endif

