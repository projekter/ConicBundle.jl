/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/BundleIdProx.hxx
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



#ifndef CONICBUNDLE_BUNDLEIDPROX_HXX
#define CONICBUNDLE_BUNDLEIDPROX_HXX


/**  @file BundleIdProx.hxx
    @brief Header declaring the class ConicBundle::BundleIdProx
    @version 1.0
    @date 2014-08-06
    @author Christoph Helmberg
*/


#include "BundleProxObject.hxx"

namespace ConicBundle {

/**@ingroup InternalBundleProxObject
 */

  //@{

 /** @brief implements the abstract interface ConicBundle::BundleProxObject for \f$\|y-\hat{y}\|_H^2\f$ with H=weight*I,
     giving rise to a pure augmented model without scaling
 */

  class BundleIdProx: public BundleProxObject
{
private:
  /// the weight for the proximal term
  CH_Matrix_Classes::Real weightu;
  /// the dimension of the identity
  CH_Matrix_Classes::Integer dim;

  /// The MinorantPointerMap serves to locate an identical MinorantPointer in a previous bundle in order to reduce the amount of computations
  typedef std::map<MinorantPointer,CH_Matrix_Classes::Integer> MinorantPointerMap;
  /// identifies which MinorantPointer was used last time in which position
  MinorantPointerMap oldmap;
  /// the old quadratic cost matrix; this is where oldmap points into 
  CH_Matrix_Classes::Symmatrix oldQ;
  /// the weight value used for computing the old quadratic cost matrix 
  CH_Matrix_Classes::Real oldweightu;
  /// the old fixed indices for which oldmap and oldQ were computed 
  CH_Matrix_Classes::Indexmatrix old_fixed_ind;
  /// for storing the relevant linear parts of the minorants
  CH_Matrix_Classes::Matrix _A;
  /// for storing the relevant linear parts of the constant minorant
  CH_Matrix_Classes::Matrix _b;
  /// for storing the offsets of the minorants
  CH_Matrix_Classes::Matrix _c;
  /// for storing the offset of the constant minorant
  CH_Matrix_Classes::Real _delta;

  // CH_Matrix_Classes::Symmatrix last_Q;  ///< for testing updates 
  //CH_Matrix_Classes::Matrix last_d;     ///< for testing updates
  // CH_Matrix_Classes::Real last_offset;  ///< for testing updates

public:
  /// initialize with dimension and weight
  BundleIdProx(CH_Matrix_Classes::Integer d=0,CH_Matrix_Classes::Real w=1.,CBout* cb=0, int cbinc=-1):
    BundleProxObject(0,false,false,cb,cbinc)
  {dim=d;weightu=CH_Matrix_Classes::max(CH_Matrix_Classes::eps_Real,w);oldweightu=1.;}
  ///
  virtual ~BundleIdProx(){}
   
  /// set the weight of the proximal term
  virtual void set_weightu(CH_Matrix_Classes::Real in_weightu)
  {weightu=CH_Matrix_Classes::max(CH_Matrix_Classes::eps_Real,in_weightu);}

  /// returns the current weight of the proximal term
  virtual CH_Matrix_Classes::Real get_weightu() const
  {return weightu;}

  

  /// returns the correction factor for the termination criterion, here min(1,1/weight)
  virtual CH_Matrix_Classes::Real get_term_corr(void) const {return CH_Matrix_Classes::min(1.,1/weightu);}

  /// returns \f$\|B\|^2_H\f$
  virtual CH_Matrix_Classes::Real norm_sqr(const CH_Matrix_Classes::Matrix& B) const
  {return weightu*CH_Matrix_Classes::ip(B,B);}

  /// returns \f$\|B\|^2_{H^{-1}}\f$
  virtual CH_Matrix_Classes::Real dnorm_sqr(const MinorantPointer& B) const
  {return B.dual_norm_squared()/weightu;}
    
  ///return true if H is of the form diagonal matrix plus Gram matrix of a low rank matrix
  virtual bool is_DLR() const
  {return true;}

  /// add H to the dense symmetric matrix as a principal submatrix starting at position start_index
  virtual int add_H(CH_Matrix_Classes::Symmatrix& big_sym,CH_Matrix_Classes::Integer start_index=0) const;

  /// adds \f$alpha*Hx\f$ to outplusHx and returns this
  virtual CH_Matrix_Classes::Matrix& add_Hx(const CH_Matrix_Classes::Matrix& x,
					    CH_Matrix_Classes::Matrix& outplusHx,
					    CH_Matrix_Classes::Real alpha=1.) const
  {return outplusHx.xpeya(x,weightu*alpha);}

  /// returns \f$H^{-1}x\f$
  virtual CH_Matrix_Classes::Matrix& apply_Hinv(CH_Matrix_Classes::Matrix& x) const 
  {return x/=weightu;}

  /// returns a suitable approximation for preconditioning, see BundleProxObject::get_precond 
  virtual void get_precond(CH_Matrix_Classes::Matrix& D,const CH_Matrix_Classes::Matrix*& Vp) const
  {D.init(dim,1,weightu);Vp=0;}
  
   
  /// computes the dual QP costs Q, d, and the constant offset to the bundle subproblem, see BundleProxObject::compute_QP_costs
  virtual int compute_QP_costs(CH_Matrix_Classes::Symmatrix& Q,
			       CH_Matrix_Classes::Matrix& d,
			       CH_Matrix_Classes::Real& offset,
			       const MinorantPointer& constant_minorant,
			       const MinorantBundle& bundle, 
			       const CH_Matrix_Classes::Matrix& y,
			       const MinorantPointer& groundset_minorant,
			       CH_Matrix_Classes::Indexmatrix* yfixed);
  

  /// updates the dual QP costs Q, d, and the constant offset to the bundle subproblem, see BundleProxObject::update_QP_costs
  virtual int update_QP_costs(CH_Matrix_Classes::Symmatrix& delta_Q, 
			      CH_Matrix_Classes::Matrix& delta_d,  
			      CH_Matrix_Classes::Real& delta_offset,
                              const MinorantPointer& constant_minorant,
                              const MinorantBundle& bundle,
			      const CH_Matrix_Classes::Matrix& center_y,
			      const MinorantPointer& groundset_minorant,
			      const MinorantPointer& delta_groundset_minorant,
			      const CH_Matrix_Classes::Indexmatrix& delta_index,
			      CH_Matrix_Classes::Indexmatrix* yfixed);

  /// when BundleSolver is called to modify the groundset it also calls this 
  virtual int apply_modification(const GroundsetModification& gsmdf);

    /** @brief in order to allow for fixed variables, this generates a clone restricted to the given indices
   */
  virtual BundleProxObject* projected_clone(const CH_Matrix_Classes::Indexmatrix& indices) 
  { return new BundleIdProx(indices.dim(),weightu,this,0); }

  /// output the description of the prox term in mfile-suitable format
  virtual int mfile_data(std::ostream& out) const;

};


  //@}
}

#endif

