/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/BundleRQBWeight.hxx
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



#ifndef CONICBUNDLE_BUNDLERQBWEIGHT_HXX
#define CONICBUNDLE_BUNDLERQBWEIGHT_HXX

#include "BundleWeight.hxx"

namespace ConicBundle {

/** @ingroup InternalBundleSolver 

*/
//@{


/** @brief Routine for selecting the weight of the quadratic/proximal term
     within BundleSolver implementing BundleWeight roughly along the paper
     C. Lemarechal, C. Sagastizabal, Variable metric bundle methods: From conceptual to implementable forms, Mathematical Programming 76 (1997) 393-410.
     Implemented obeserving details of the email communication with 
     Claudia Sagastizabal.
      
*/

class BundleRQBWeight: public BundleWeight
{
  Groundset* groundset;  ///< the groundset, may not be zero in init and *_update
  BundleModel* model;    ///< may be zero at any time
  MinorantPointer center_aggr;  ///< the aggregate stored for the current center, may be empty if not available
  CH_Matrix_Classes::Real weight; ///< the current weight

  CH_Matrix_Classes::Real minweight; ///< lower bound on the weight, no bound if negative
  CH_Matrix_Classes::Real maxweight; ///< upper bound on the weight, no bound of negative
  
  bool weightchanged;    ///< true if last choose_* call modified the weight
  bool next_weight_set;  ///< true if set_nextweight was just called, reset at *_update
  
  CH_Matrix_Classes::Real m1; ///< parameter for descent step
  CH_Matrix_Classes::Real m2; ///< parameter in (m1,1) for using quasi udpate
  CH_Matrix_Classes::Real m3; ///< parameter for "accepting null steps"
  CH_Matrix_Classes::Real eta; ///< parameter for "zero subgradients"
  CH_Matrix_Classes::Real tL; ///< lower bound for recording extrapolation
  CH_Matrix_Classes::Real tR; ///< upper bound for recording interpolation
  CH_Matrix_Classes::Real t; ///< recorded inter/extrapolation
  CH_Matrix_Classes::Real frel; ///< realtive precision w.r.t. function value  
  CH_Matrix_Classes::Real delhat1; ///< oldval-modelval at t==1 unchanged

  /// computes cand-center into delta_subg and retunrs 0 if both are not empty; otherwise it returns 1
  int delta_subg(CH_Matrix_Classes::Matrix& delta_subg,MinorantPointer& center,MinorantPointer& cand);
    
public:
  /// bwp may be used to communicate the previous values used by another routine 
  BundleRQBWeight(BundleWeight* bwp=0,
		  const CBout* cbo=0,
		  int incr=-1);

  /// bwp may be used to communicate the previous values used by another routine 
  BundleRQBWeight(CH_Matrix_Classes::Real m1,
		  CH_Matrix_Classes::Real m2=.2,
		  CH_Matrix_Classes::Real m3=1.,
		  CH_Matrix_Classes::Real eta=1e-6,
		  BundleWeight* bwp=0,
		  const CBout* cbo=0,
		  int incr=-1);
  ///
  ~BundleRQBWeight(){}

 
  ///set default values for 'constant' parameters, e.g. minweight and maxweight
  virtual void set_defaults(); 

  ///reset all adaptive variables and parameters
  virtual void clear();        
  

  ///compute first weight and set some parameters
  int init(CH_Matrix_Classes::Real aggr_dnmormsqr,Groundset* groundset,BundleModel* model);
  
  /// <=0 leaves everything unchanged and does nothing
  void set_next_weight(CH_Matrix_Classes::Real u)
  { if (u<=0.) return; 
  weight=CH_Matrix_Classes::max(u,1e-10); weightchanged=true;next_weight_set=true;}
  
  /// <=0 means no bound 
  void set_minweight(CH_Matrix_Classes::Real mw)
  { 
    minweight=mw; 
    if (minweight>0){
      if ((weight>0)&&(weight<minweight))
	weight=minweight;
      if ((maxweight>0)&&(maxweight<minweight)) 
	maxweight=minweight;
    }
  }

  /// true if the next weight was prespecified externally
  bool get_next_weight_set() const {return next_weight_set;}

  /// 
  CH_Matrix_Classes::Real get_minweight() const {return minweight;}
  
  /// <=0 means no bound
  virtual void set_maxweight(CH_Matrix_Classes::Real mw)
  { 
    maxweight=mw; 
    if (maxweight>0){
      if (weight>maxweight){
	weight=maxweight;
      }
      if ((minweight>0)&&(minweight>maxweight)) 
	minweight=maxweight;
    }
  }
  
  ///
  CH_Matrix_Classes::Real get_maxweight() const {return maxweight;}

  /// returns current value of the weight
  CH_Matrix_Classes::Real get_weight() const;
  
  /// returns true if last call of *_update modified current value of tau, else 0
  bool weight_changed() const;
  
  /// determine next weight after a descent step
  int descent_update(CH_Matrix_Classes::Real newval,
		     CH_Matrix_Classes::Real oldval,
		     CH_Matrix_Classes::Real modelval,
		     const CH_Matrix_Classes::Matrix& y, 
		     const CH_Matrix_Classes::Matrix& newy,
		     CH_Matrix_Classes::Real normsubg2,
		     BundleProxObject* Hp);
  
  /// determine next weight after a null step
  int nullstep_update(CH_Matrix_Classes::Real newval,
		      CH_Matrix_Classes::Real oldval,
		      CH_Matrix_Classes::Real modelval,
		      const CH_Matrix_Classes::Matrix& y, 
		      const CH_Matrix_Classes::Matrix& newy,
		      MinorantPointer& new_minorant,
		      MinorantPointer& aggregate,
		      CH_Matrix_Classes::Real nullstep_bound,
		      CH_Matrix_Classes::Real normsubg2,
		      BundleProxObject *Hp);
      
  /// reinitialize after modifications
  virtual int apply_modification(const GroundsetModification& gsmdf);

};

}

#endif

