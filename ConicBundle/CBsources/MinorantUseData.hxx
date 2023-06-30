/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/MinorantUseData.hxx
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



#ifndef CONICBUNDLE_MINORANTUSEDATA_HXX
#define CONICBUNDLE_MINORANTUSEDATA_HXX


/**  @file MinorantUseData.hxx
    @brief Header declaring the class ConicBundle::MinorantUseData
    @version 1.0
    @date 2016-02-17
    @author Christoph Helmberg
*/

#include <map>
#include "CBSolver.hxx"
#include "CBout.hxx"
#include "matrix.hxx"


namespace ConicBundle {

  /**@ingroup  InternalBundleSolver
   */


   /** @brief stores use information for the model that the minorant is generated for or to the data the model refers to

       The data always relates to the function/model the MinorantPointer was generated for;
       If a MinorantPointer from another function/model makes use of this, it generates
       its own MinorantUseData pointing to this MinorantUseData. This make it possible to
       recognize/enable dynamic function changes at any level.
    */

  class MinorantUseData : public CBout {
  private:
    friend class MinorantPointer;

    /// the number of MinorantPointer instance pointing to this; once it reaches 0 the corresponding Minorant Pointer will delete it
    CH_Matrix_Classes::Integer use_cnt;
    /// the value of the modification id of the function at the time of generation or when modified; if minorant!=0 and modification_id<0, the minorant is no longer valid
    CH_Matrix_Classes::Integer modification_id;
    /// the number of minorants aggregated in this (only caused by MinorantPointer::aggregate
    CH_Matrix_Classes::Integer aggr_stat;
    /// gives the id of the last primal extender call, -1 for failed, 0 initially
    CH_Matrix_Classes::Integer prex_id;
    /// value by which the minorant or the MinorantUseData has to be scaled
    CH_Matrix_Classes::Real scaleval;

    /// if minorant!=0, this records evaluations for given point ids
    mutable std::map<CH_Matrix_Classes::Integer, CH_Matrix_Classes::Real> evals;

    /// this points to the minorant or maybe to nothing if not valid
    Minorant* minorant;
    /// if there is no minorant, it should be in here (recursively)
    MinorantUseData* md;

    /// forbidden, blocked deriberately
    MinorantUseData(const MinorantUseData&) :CBout() {
      assert(false);
    }
    /// forbidden, blocked deriberately
    MinorantUseData& operator=(const MinorantUseData&) {
      assert(false); return *this;
    }

  public:
    /// constructor for a new minorant with scaling @a sval and modification id @a modfi_id
    MinorantUseData(Minorant* mp,
      CH_Matrix_Classes::Real sval,
      CH_Matrix_Classes::Integer modif_id) :
      CBout(), modification_id(modif_id), scaleval(sval) {
      assert(mp); use_cnt = 0; minorant = mp; md = 0; aggr_stat = 0; prex_id = 0;
    }
    /// constructor for a recursively containing further MinorantUseData with an additional scaling factor @a sval  
    MinorantUseData(MinorantUseData* mdp, CH_Matrix_Classes::Real sval) :
      CBout(), scaleval(sval) {
      assert((mdp == 0) || (mdp->use_cnt >= 1));
      use_cnt = 0; minorant = 0; md = mdp; md->use_cnt++;
      modification_id = -1; aggr_stat = 0; prex_id = 0;
    }
    ///
    ~MinorantUseData();

    /// by a recursive call return the collected scaling value and the final minorant
    int get_scaleval_and_minorant(CH_Matrix_Classes::Real& sv, Minorant*& m);

    /// check validity recursively
    bool valid() const;

    /// return the final minorant (by a recursive call) or 0 if there is none
    Minorant* get_minorant() const;

    /// returns the modification id also for overwriting
    CH_Matrix_Classes::Integer& set_modification_id();

    /// returns the modification id
    CH_Matrix_Classes::Integer get_modification_id() const;

    /// sets the modification_id to its new id and reinitializes the evaluation map
    int synchronize_ids(CH_Matrix_Classes::Integer new_modification_id,
      CH_Matrix_Classes::Integer new_center_id,
      CH_Matrix_Classes::Integer old_center_id,
      CH_Matrix_Classes::Integer new_cand_id,
      CH_Matrix_Classes::Integer old_cand_id,
      CH_Matrix_Classes::Integer new_prex_id);

    /// returns true if the minorant is an aggregate of several minorants
    bool aggregate() const;

    /// add @a n to the number couting the aggregations
    int aggregated(CH_Matrix_Classes::Integer n);

    /// return the offset of the minorant
    CH_Matrix_Classes::Real offset() const;

    /// return coefficient @a i of the minorant
    CH_Matrix_Classes::Real coeff(CH_Matrix_Classes::Integer i) const;

    /// carries through the scaling for the underlying minorant, afterwards scaleval==1., may only be carried out if use_cnt==1
    int scale(CH_Matrix_Classes::Real factor);

    /// apply the primal extender to the minorant if @a in_prex_id indicates a new extension
    int call_primal_extender(PrimalExtender& prex, CH_Matrix_Classes::Integer in_prex_id);

    /// evaluate the minorant for @a y unluess @a yid allows to retrieve a previous evaluation
    CH_Matrix_Classes::Real evaluate(CH_Matrix_Classes::Integer yid, const CH_Matrix_Classes::Matrix& y, bool with_constant = true) const;

    /// returns true if not valid or use_cnt==1 recursively 
    bool one_user() const;


  };


  //@}


}


#endif

