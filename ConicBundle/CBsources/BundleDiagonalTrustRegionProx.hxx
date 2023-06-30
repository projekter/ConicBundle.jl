/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/BundleDiagonalTrustRegionProx.hxx
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



#ifndef CONICBUNDLE_BUNDLEDIAGONALTRUSTREGIONPROX_HXX
#define CONICBUNDLE_BUNDLEDIAGONALTRUSTREGIONPROX_HXX


/**  @file BundleDiagonalTrustRegionProx.hxx
    @brief Header declaring the class ConicBundle::BundleDiagonalTrustRegionProx
    @version 1.0
    @date 2016-08-06
    @author Christoph Helmberg
*/


#include "BundleProxObject.hxx"

namespace ConicBundle {

/**@ingroup InternalBundleProxObject
 */

  //@{

 /** @brief implements the abstract interface ConicBundle::BundleProxObject for \f$\|y-\hat{y}\|_H^2\f$ with H=D+weight*I,
     where D is a diagonal matrix, giving rise to an augmented model with diagonal scaling
 */

class BundleDiagonalTrustRegionProx: public BundleProxObject
{
private:
  /// the current weightu added into the diagonal D
  CH_Matrix_Classes::Real weightu;

  /// the diaongal of the diagonal matrix with the weightu added into it
  CH_Matrix_Classes::Matrix D;

  /// the correction value for correcting termination precision
  CH_Matrix_Classes::Real corr_val;

  /// used for the diagonal scaling heuristic
  CH_Matrix_Classes::Matrix update_Dvalue;
  
  /// The MinorantPointerMap serves to locate an identical MinorantPointer in a previous bundle in order to reduce the amount of computations
  typedef std::map<MinorantPointer,CH_Matrix_Classes::Integer> MinorantPointerMap;
  /// identifies which MinorantPointer was used last time in which position
  MinorantPointerMap oldmap;

  /// the old quadratic cost matrix; this is where oldmap points into 
  CH_Matrix_Classes::Symmatrix oldQ;

  /// the old fixed indices for which oldmap and oldQ were computed 
  CH_Matrix_Classes::Indexmatrix old_fixed_ind;

  /// if values should be computed for a new subset indices, this is stored here
  const CH_Matrix_Classes::Indexmatrix* new_indices;

  /// for judging whether a dynamic diagonal contribution does to aggr what it should ...
  const CH_Matrix_Classes::Matrix* aggr;

  /// for inserting some conservatism in dynamic updating of the diagonal
  CH_Matrix_Classes::Matrix old_D;

  ///corrects excessive step reductions by dynamic scaling, recomputed in descent steps
  CH_Matrix_Classes::Real old_damping;

  /// in add_dynamic_scaling() the stack serves to transform the given data
  std::vector<const AffineFunctionTransformation*> aft_stack;

  // CH_Matrix_Classes::Symmatrix last_Q;  ///< for testing updates 
  // CH_Matrix_Classes::Matrix last_d;     ///< for testing updates
  // CH_Matrix_Classes::Real last_offset;  ///< for testing updates

  
  /// computes the correction value, here min(1,dim/trace(D)))
  void compute_corr()
  //{ corr_val=sum(D)-weightu*D.dim(); if (corr_val>D.dim()) corr_val=D.dim()/corr_val; else corr_val=1.;}
  {corr_val=CH_Matrix_Classes::min(1.,D.rowdim()/sum(D));}

public:
  /// initialize to this diagonal matrix
  BundleDiagonalTrustRegionProx(const CH_Matrix_Classes::Matrix& Din,
				VariableMetricSelection* vp=0,
				bool local_scaling=false,
				bool bounds_scaling=false,
				CBout* cb=0,
				int cbinc=-1):
    BundleProxObject(vp,local_scaling,bounds_scaling,cb,cbinc),
    weightu(1.),D(Din)
  {D+=weightu;compute_corr();}
  
  /// initialize to a diaognal matrix d*identity of dimesion dim
  BundleDiagonalTrustRegionProx(CH_Matrix_Classes::Integer dim,
				CH_Matrix_Classes::Real d,
				VariableMetricSelection* vp=0,
				bool local_scaling=false,
				bool bounds_scaling=false,
				CBout* cb=0,
				int cbinc=-1):
    BundleProxObject(vp,local_scaling,bounds_scaling,cb,cbinc),
    weightu(1.),D(dim,1,weightu+d)
  {assert(dim>=0);assert(d>=0.);compute_corr();}

  /// initialize to a zero diaognal matrix of dimesion dim
  BundleDiagonalTrustRegionProx(CH_Matrix_Classes::Integer dim=0,
				VariableMetricSelection* vp=0,
				bool local_scaling=false,
				bool bounds_scaling=false,
				CBout* cb=0,
				int cbinc=-1):
    BundleProxObject(vp,local_scaling,bounds_scaling,cb,cbinc),
    weightu(1.),D(dim,1,weightu)
  {assert(dim>=0);compute_corr();}

  /// destructor
  virtual ~BundleDiagonalTrustRegionProx(){}

  /// set the weight of the proximal term
  void set_weightu(CH_Matrix_Classes::Real in_weightu);

  /// returns the current weight of the proximal term
  CH_Matrix_Classes::Real get_weightu() const
  {return weightu*factor;}
  
  ///returns a correction factor for termination precision if the quadratic term is strong
  CH_Matrix_Classes::Real get_term_corr() const
  {return corr_val*factor;}

  /// returns the diagonal D of the diagonal scaling matrix
  const CH_Matrix_Classes::Matrix& get_D() const {return D;}

  /// set the diagonal (it needs to be >=0 but this is not checked)
  void set_D(CH_Matrix_Classes::Matrix& in_D)
  {D=in_D+weightu;compute_corr();oldmap.clear();oldQ.init(0,0.);}

  /// returns the dimension of the diagonal
  CH_Matrix_Classes::Integer dim() const {return D.dim();}

  /// return the i-th element of the diagonal matrix D
  CH_Matrix_Classes::Real operator()(CH_Matrix_Classes::Integer i) const {return D(i);}

  /// return the i-th element of the diagonal matrix D
  CH_Matrix_Classes::Real operator[](CH_Matrix_Classes::Integer i) const {return D[i];}
   
  /// returns \f$\|B\|^2_H\f$
  virtual CH_Matrix_Classes::Real norm_sqr(const CH_Matrix_Classes::Matrix& B) const;

  /// returns \f$\|B\|^2_{H^{-1}}\f$
  virtual CH_Matrix_Classes::Real dnorm_sqr(const MinorantPointer& B) const
  {return B.dual_norm_squared(&D)/factor;}
  
  ///return true if H is of the form diagonal matrix plus Gram matrix of a low rank matrix
  virtual bool is_DLR() const
  {return true;}

  /// add H to the dense symmetric matrix as a principal submatrix starting at position start_index
  virtual int add_H(CH_Matrix_Classes::Symmatrix& big_sym,CH_Matrix_Classes::Integer start_index=0) const;

  /// adds \f$alpha*Hx\f$ to outplusHx and returns this
  virtual CH_Matrix_Classes::Matrix& add_Hx(const CH_Matrix_Classes::Matrix& x,
					    CH_Matrix_Classes::Matrix& outplusHx,
					    CH_Matrix_Classes::Real alpha=1.) const
  { return outplusHx.xpeya(D%x,alpha); }

  /// returns \f$H^{-1}x\f$
  virtual CH_Matrix_Classes::Matrix& apply_Hinv(CH_Matrix_Classes::Matrix& x) const ;

  /// returns a suitable approximation for preconditioning, see BundleProxObject::get_precond 
  virtual void get_precond(CH_Matrix_Classes::Matrix& inD,const CH_Matrix_Classes::Matrix*& Vp) const
  { inD=D; Vp=0; }


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
  virtual BundleProxObject* projected_clone(const CH_Matrix_Classes::Indexmatrix& indices ) 
  {
    CH_Matrix_Classes::Matrix tmpD=D(indices);
    tmpD-=weightu;
    BundleDiagonalTrustRegionProx* pp=new BundleDiagonalTrustRegionProx(tmpD,0,get_use_local_metric(),get_use_bounds_scaling(),this,0);
    pp->set_weightu(weightu);
    pp->apply_factor(factor);
    if (get_variable_metric_selection())
      pp->set_variable_metric_selection(get_variable_metric_selection()->clone_VariableMetricSelection());
    return pp;
  }

  /** @brief this implementation supports a diagonal scaling heuristic for 
      bounds in the groundset, therefore the following routine has to return true.
   */
  virtual bool supports_diagonal_bounds_scaling() const {return true;}
  
  /** @brief if supported, D_update has to contain nonnegative numbers that
      are permanently added to the diagonal here. It is important to keep
      track of this change only if afterwards update_QP_costs is called before
      compute_QP_costs. In this case the only nonzero enries in D_update must
      be those of delta_index
   */
  virtual int diagonal_bounds_scaling_update(const CH_Matrix_Classes::Matrix& D_update)
  {
    assert(min(D_update)>=0.);
    update_Dvalue=D_update;D+=D_update;compute_corr();
    oldmap.clear();oldQ.init(0,0.);
    return 0;
  }
  
  /// returns true if dynamic scaling with low rank matrices is supported
  virtual bool supports_lowrank_variable_metric() const {return false;}

  /// returns true if dynamic scaling with diagonal matrices is supported
  virtual bool supports_diagonal_variable_metric() const {return true;}

  /// see VariableMetric
  virtual int apply_variable_metric(VariableMetricModel* groundset,
				    VariableMetricModel* model,
				    const CH_Matrix_Classes::Matrix& aggr,
				    CH_Matrix_Classes::Integer y_id,
				    const CH_Matrix_Classes::Matrix& y,
				    bool descent_step,
				    CH_Matrix_Classes::Real& current_weight,
				    CH_Matrix_Classes::Real model_maxviol,
				    const CH_Matrix_Classes::Indexmatrix* new_indices =0);
  
  /// see BundleProxObject::add_dynamic_scaling()
  virtual int add_variable_metric(CH_Matrix_Classes::Matrix& diagH,
				  CH_Matrix_Classes::Matrix& vecH);

  /// see  BundleProxObject::push_aft();
  virtual int push_aft(const AffineFunctionTransformation* aft)
  { aft_stack.push_back(aft); return 0;}
  
  /// see  BundleProxObject::pop_aft();
  virtual int pop_aft()
  { 
    aft_stack.pop_back(); 
    return 0;
  }

  /// output the description of the scaling in mfile-suitable format
  virtual int mfile_data(std::ostream& out) const;
};


  //@}
}

#endif

