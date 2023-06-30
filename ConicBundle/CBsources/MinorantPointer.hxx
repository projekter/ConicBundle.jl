/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/MinorantPointer.hxx
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



#ifndef CONICBUNDLE_MINORANTPOINTER_HXX
#define CONICBUNDLE_MINORANTPOINTER_HXX


/**  @file MinorantPointer.hxx
    @brief Header declaring the class ConicBundle::MinorantPointer
    @version 1.0
    @date 2016-02-17
    @author Christoph Helmberg
*/

#include "MinorantUseData.hxx"
#include "sparsmat.hxx"

namespace ConicBundle {

  class GroundsetModification;

/**@ingroup InternalBundleSolver
 */

  class MinorantPointer;

  /// a bundle is a vector with MinorantPointer entries
  typedef std::vector<MinorantPointer> MinorantBundle;

/** @brief points to MinorantUseData that may be shared by many and allows computations with Minorants

    A minorant pointer is _empty_ if it does not point to any MinorantUseData
 */

class MinorantPointer: public CBout
{
private:
  /// if null, it is regarded as not initialized or _empty_
  MinorantUseData* md;
  /// if -1 it is invalid or _empty_, otherwise it gives the modification id of the function that it was created for, at that time identical to the one in the MinorantUseData

  /// reduces the use_cnt of the MinorantUseData it points to and deletes it if this reaches 0; afterwards it is _empty_
  void delete_data();

  /// first deletes the old data and then fills the information with the new
  void new_data(MinorantUseData* in_md);

  /// check whether the minorant is owend by this pointer (use_cnt==1) and if not, clone it first
  int prepare_for_changes(CH_Matrix_Classes::Real factor=1.,bool with_primal=false);

public:  
  /// calls new_data with the data of mp and copies the modification_id
  MinorantPointer& operator=(const MinorantPointer& mp);

  /// if factor!=1 it generates another MinorantUseData referring to the one of mp, otherwise it simply uses the same one (empty stays empty)
  void init(const MinorantPointer& mp,CH_Matrix_Classes::Real factor=1.,bool enforce_copy=false);

  /// if mp==0 it becomes empty, otherwise it creates and then points to a MinorantUseData for holding mp with this modification_id and factor
  void init(Minorant* mp,CH_Matrix_Classes::Integer modification_id=-1,CH_Matrix_Classes::Real factor=1.);
  
  /// declares the pointer _empty_
  MinorantPointer():md(0)
  {}

  /// initialize this to point to in_md
  MinorantPointer(MinorantUseData* in_md):CBout(),md(0)
  {new_data(in_md);}

  /// calls new_data
  MinorantPointer(const MinorantPointer& mp):CBout(mp),md(0)
  {new_data(mp.md);}

  /// calls init(const MinorantPointer&,CH_Matrix_Classes::Real)
  MinorantPointer(const MinorantPointer& mp,CH_Matrix_Classes::Real factor):
    CBout(mp),md(0)
  {init(mp,factor);}
  
  /// calls init(Minorant*,CH_Matrix_Classes::Integer,CH_Matrix_Classes::Real)
  MinorantPointer(Minorant* mnrt,CH_Matrix_Classes::Integer modification_id,CH_Matrix_Classes::Real factor=1.):
    CBout(),md(0)
  {init(mnrt,modification_id,factor);}
  
  /// calls delete_data
  virtual ~MinorantPointer(){ delete_data();}

  /// afterwards the pointer is _empty_ (calls delete_data())
  void clear() { delete_data();}

  /// if not empty it sets the modification_id to its new id and reinitializes the evaluation map
  int synchronize_ids(CH_Matrix_Classes::Integer new_modification_id,
		       CH_Matrix_Classes::Integer new_center_id,
		       CH_Matrix_Classes::Integer old_center_id,
		       CH_Matrix_Classes::Integer new_cand_id,
		       CH_Matrix_Classes::Integer old_cand_id,
		       CH_Matrix_Classes::Integer new_prex_id=0);

  /// returns true if the pointer is _empty_
  bool empty() const {return (md==0);}

  /// returns true if the pointer is not empty and the data is valid
  bool valid() const {return (md!=0)&&(md->valid());}
 
  /// returns true if the pointer is not empty but all entrys (also the offset) are zero 
  bool zero() const;
 
  /**@brief returns ture if it points to a combination of minorants from the same function;
     
     A minorant is considered an aggregate if it is obtained from an aggregate 
     or it was computed by a call to aggregate(MinorantBundle&,Real); 
     It is not an aggregate if it is the sum of non aggregate minorants from 
     different functions.

     An _empty_ pointer is not an aggregate.
  */
  bool aggregate() const;

  /// returns true if not valid or this is the only active pointer to the minorant
  bool one_user() const;

  /// returns the Minorant *this points to or 0 if _empty_
  const Minorant* get_minorant() const
  {if (empty()) return 0; return md->get_minorant();}

  /// returns the Minorant *this points to with its scaling value or 0 if _empty_
  int get_scaleval_and_minorant(CH_Matrix_Classes::Real& sv,Minorant*& m)
  {if (empty()) { m=0;sv=0.;return 0;} return md->get_scaleval_and_minorant(sv,m);}

  /// returns the primal of the Minorant *this points to or 0 if _empty_ (first carrying out any pending scalings on the minorant)
  const PrimalData* get_primal(); 

  /// returns the offset of the minorant (including the internal scalings) or CB_minus_infinity if _empty_
  CH_Matrix_Classes::Real offset() const;

  /// returns coefficient i of the minorant (including the internal scalings) or 0. if _empty_
  CH_Matrix_Classes::Real coeff(CH_Matrix_Classes::Integer i) const;

  /// returns the number of nonzero coefficients
  CH_Matrix_Classes::Integer nonzeros() const
  {if (empty()) return 0; return md->get_minorant()->nonzeros();}

  /// multiply the MinorantPointer by val (this is an external factor for the minorant and possibly its primal information, but it will be used in aggreagation and when retrieving the approximate primal)
  int scale(CH_Matrix_Classes::Real val);

  /// executes prex on the primal if its id is smaller then prex_id; returns != if this fails or not primal data is available 
  int call_primal_extender(PrimalExtender& prex,CH_Matrix_Classes::Integer prex_id);

  /// if valid and the local modification_id is smaller than mod_id, it modifies the coefficients as described by GroundsetModification; if successful, the mod_id is assigend and the function returns 0; otherwise it is invalidated and returns !=0; the costs of gsmdf are used only if apply_gsmdf_costs=false in whic case mex must be NULL
  int apply_modification(const GroundsetModification& gsmdf,
			 CH_Matrix_Classes::Integer mod_id,
			 MinorantExtender* mex,
			 bool apply_gsmdf_costs=false);

  /// negative ids are allowed and indicate there is no need to memorize this result, returns CB_minus_infinity if empty, otherwise offset+ip(minorant,y)
  CH_Matrix_Classes::Real evaluate(CH_Matrix_Classes::Integer yid,
				   const CH_Matrix_Classes::Matrix& y,bool with_constant=true) const;
  
  /// add the offset; if empty, initialize to offset with zero linear part, if not the only user (use_cnt>1), clone it first
  int add_offset(CH_Matrix_Classes::Real offset);

  /**@brief store/add the minorant in offset and a matrix column, possibly skipping indices skip_fixed. The coordinats of the latter ones are multiplied by the given values or 0 and added to offset

     The dimensions of the matrix must already fit the requirements on input.
     If add is false (as by default), the full length of the column 
     is initialized to the gradient and filled up with zeros where needed.
     skip_fixed and fixed_vals must both be given or both not be given.
     If they are both given, both must have the same length.
     If skip_fixed is given, it must have strictly increasing indices.
     The skip part of the routine is used in implementations of BundleScaling::get_QP_costs.
   */
  int get_minorant(CH_Matrix_Classes::Real& offset,
		   CH_Matrix_Classes::Matrix& mat,
		   CH_Matrix_Classes::Integer column,
		   CH_Matrix_Classes::Real alpha=1.,
		   bool add=false,
		   const CH_Matrix_Classes::Indexmatrix* skip_fixed=0,
		   const CH_Matrix_Classes::Matrix* fixed_vals=0) const; 


  /**@brief store/add the minorant in/to mp, possibly scaled by alpha

     If mp is empty, store it there, if mp is not empty, add it.

     if *this is empty, it causes an error.
   */
  int get_minorant(MinorantPointer& mp,
		   CH_Matrix_Classes::Real alpha=1.) const; 

  /**@brief store/add the minorant in/to mp, possibly scaled by alpha and transformed by sp, which possibly requires only the indices of provided_row_indices to compute possibly only the indices in needed_col_indices 

     If mp is empty, store it there, if mp is not empty, add it.

     if sp==0, it is treated as the identity

     If provided_row_indices or needed_col_indices is given, its indices must
     be in strictly increasing order.  For sp==0 both must coincide and will
     probably be ignored by just returning a scaled reference to *this in mp.

     if *this is empty, it causes an error.
   */
  int get_minorant(MinorantPointer& mp,
		   CH_Matrix_Classes::Real alpha,
		   const CH_Matrix_Classes::Sparsemat* sp,
		   const CH_Matrix_Classes::Indexmatrix* provided_row_indices=0,
		   const CH_Matrix_Classes::Indexmatrix* needed_col_indices=0,
		   bool enforce_copy=false) const; 

  /**@brief collect the aggregate as the nonnegative linear combination of the bundle minorants (here *this may not be part of the bundle!). If *this is not empty, add this aggregate to *this. Only aggregation carries along primal data!

     All minorants in the bundle need to be valid. Aggregating over no minorants 
     is treated as a zero minorant without any primal information.
   */
  int aggregate(const MinorantBundle& minorants,const CH_Matrix_Classes::Matrix& coeff,CH_Matrix_Classes::Real factor=1.);

  /**@brief aggregate the minorant*itsfactor to *this, both must be valid for this. Only aggregation carries along primal data!
  */
  int aggregate(const MinorantPointer& minorant,double itsfactor=1.);

  /**@brief computes the inner product of the two minorants; if skip_fixed!=NULL the corrsponding indices are not considered, if ipdiag!=0 the inner product is taken with respect to this diagonal matrix, i.e.  sum_i mp(i)*(*this)(i)*(*ipdiag)(i)
   */
  CH_Matrix_Classes::Real ip(const MinorantPointer& mp,const CH_Matrix_Classes::Indexmatrix* skip_fixed=0,const CH_Matrix_Classes::Matrix* ipdiag=0) const;

  /**@brief computes the inner product with m; if ipdiag!=0 the inner product is taken with respect to this diagonal matrix, i.e.  sum_i mp(startindex_m+i)*(*this)(i)*(*ipdiag)(i) */
  CH_Matrix_Classes::Real ip(const CH_Matrix_Classes::Matrix& m,const CH_Matrix_Classes::Matrix* ipdiag=0,CH_Matrix_Classes::Integer startindex_m=0) const;

  /// they are equal if they point to the same object, compares the addresses of this objects
  bool operator<(const MinorantPointer& mp) const
  {return (md<mp.md);}  
  
  /// they are equal if they point to the same object, compares the addresses of this objects
  bool operator==(const MinorantPointer& mp) const
  {return (md==mp.md);}
 
  /// they are equal if they point to the same object, compares the addresses of this objects
  bool operator>(const MinorantPointer& mp) const
  {return (md>mp.md);}

  /// they are equal if they point to the same object or are both 0. If not, they differ if their matrix representations differ; if not, they differ if the entries differ by at least tol*(1.+fabs(this->offset())
  bool equals(const MinorantPointer& mp,CH_Matrix_Classes::Real tol=1e-10) const;

  ///Compute the norm squared of this for the given diagonal matrix D (identity if not given), i.e. \f$\|(*this)\|^2_{D}\f$
  CH_Matrix_Classes::Real norm_squared(const CH_Matrix_Classes::Matrix *D=0) const;

  ///Compute the dual norm squared of this for the given diagonal matrix D (identity if not given), i.e. \f$\|(*this)\|^2_{D^{-1}}\f$
  CH_Matrix_Classes::Real dual_norm_squared(const CH_Matrix_Classes::Matrix *D=0) const;

  /// computes and returns C=alpha*(*this)*B+beta*C where B and *this may be transposed and *this is considered to be a column with thisindex in a bigger matrix
  CH_Matrix_Classes::Matrix& left_genmult(const CH_Matrix_Classes::Matrix& B,
					  CH_Matrix_Classes::Matrix& C,
					  CH_Matrix_Classes::Real alpha=1.,
					  CH_Matrix_Classes::Real beta=0.,
					  int thistrans=0,
					  int btrans=0,
					  CH_Matrix_Classes::Integer thisindex=0) const;

  /// computes and returns C=alpha*A*(*this)+beta*C where A and *this may be transposed and *this is considered to be a column with thisindex in a bigger matrix
  CH_Matrix_Classes::Matrix& right_genmult(const CH_Matrix_Classes::Matrix& A,
					   CH_Matrix_Classes::Matrix& C,
					   CH_Matrix_Classes::Real alpha=1.,
					   CH_Matrix_Classes::Real beta=0.,
					   int atrans=0,
					   int thistrans=0,
					   CH_Matrix_Classes::Integer thisindex=0) const;

  ///output the Minorant in a nice format
  std::ostream& display(std::ostream& out,int precision=8) const;
};

  /// computes and returns C=alpha*A*B+beta*C where A and B  may be transposed and A is considered to have the gradients of the minorants of the bundle as columns. Because the bundle has no fixed row dimension, the dimension of C has to be compatible at input to serve as size in the untranposed case. If Coffset is given and A is not transposed, the offsets are treated as an extra row to be computed into Coffset, if A is transposed, the vector of offsets is added to Coffset with the same alpha and beta interpretation but without multiplication
  CH_Matrix_Classes::Matrix& genmult(const MinorantBundle& A,
				     const CH_Matrix_Classes::Matrix& B,
				     CH_Matrix_Classes::Matrix& C,
				     CH_Matrix_Classes::Real alpha=1.,
				     CH_Matrix_Classes::Real beta=0.,
				     int Atrans=0,
				     int Btrans=0,
				     CH_Matrix_Classes::Matrix* Coffset=0);

  /// computes and returns C=alpha*A*B+beta*C where A and B  may be transposed and B is considered to have the gradients of the minorants of the bundle as columns. Because the bundle has no fixed row dimension, the dimension of C has to be compatible at input to serve as size in the tranposed case. If Coffset is given and A is transposed, the offsets are treated as an extra column to be computed into Coffset, if A is not transposed, the row vector of offsets is added to the presized row vector Coffset with the same alpha and beta interpretation but without multiplication
  CH_Matrix_Classes::Matrix& genmult(const CH_Matrix_Classes::Matrix& A,
				     const MinorantBundle& B,
				     CH_Matrix_Classes::Matrix& C,
				     CH_Matrix_Classes::Real alpha=1.,
				     CH_Matrix_Classes::Real beta=0.,
				     int Atrans=0,
				     int Btrans=0,
				     CH_Matrix_Classes::Matrix* Coffset=0);


  //@}


}


#endif

