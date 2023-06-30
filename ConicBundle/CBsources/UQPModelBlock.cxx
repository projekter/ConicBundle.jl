/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/UQPModelBlock.cxx
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



#include "UQPSumModelBlock.hxx"
#include "UQPConeModelBlock.hxx"


namespace ConicBundle {

  UQPModelBlock::~UQPModelBlock()
  {}

  // *****************************************************************************
//                                push_aft
// *****************************************************************************

  int UQPModelBlock::push_aft(const AffineFunctionTransformation* aft,
			      const CH_Matrix_Classes::Indexmatrix* global_indices,
			      const CH_Matrix_Classes::Indexmatrix* local_indices,
			      std::map<MinorantPointer,MinorantPointer>* precomputed)
  {
    assert(constant_minorant.size()>0);
    assert(bundle.size()>0);
    constant_minorant.push_back(constant_minorant.back());
    bundle.push_back(bundle.back());
    if (aft==0)
      return 0;
    
    MinorantPointer& cm=get_constant_minorant();
    MinorantBundle& bun=get_bundle();

    MinorantPointer tmpm;
    int err=0;
    if (precomputed==0){
      if (!cm.empty()){
	if (aft->transform_minorant(tmpm,cm,1.,true,local_indices,global_indices)){
	  err++;
	  if (cb_out(0)){
	    get_out()<<"\n**** ERROR: QPModelBlock::apply_aft(..): transform_minorant failed for constant_minorant"<<std::endl;
	  }
	}
	cm=tmpm;
      }
      for(unsigned int i=0;i<bun.size();i++){
	tmpm.clear();
	if (aft->transform_minorant(tmpm,bun[i],1.,false,local_indices,global_indices)){
	  err++;
	  if (cb_out(0)){
	    get_out()<<"\n**** ERROR: QPModelBlock::apply_aft(..): transform_minorant failed for minorant "<<i<<std::endl;
	  }
	}
	bun[i]=tmpm;
      }
      return err;
    }

    // precomputed is available, deal with it
    std::map<MinorantPointer,MinorantPointer>::iterator mapit;
    if (!cm.empty()){
      mapit=precomputed->find(cm);
      if ((mapit==precomputed->end())||(! mapit->second.valid())){
	if (aft->transform_minorant(tmpm,cm,1.,false,local_indices,global_indices)){
	  if (cb_out(0)){
	    get_out()<<"\n**** ERROR: QPModelBlock::apply_aft(..): transform_minorant failed for constant_minorant"<<std::endl;
	  }
	  err++;
	}
	(*precomputed)[cm]=tmpm;
	cm=tmpm;
      }
      else {
	cm=mapit->second;
      }
    }
    for(unsigned int i=0;i<bun.size();i++){
      mapit=precomputed->find(bun[i]);
      if ((mapit==precomputed->end())||(!mapit->second.valid())){
	tmpm.clear();
	if (aft->transform_minorant(tmpm,bun[i],1.,false,local_indices,global_indices)){
	  if (cb_out(0)){
	    get_out()<<"\n**** ERROR: QPModelBlock::apply_aft(..): transform_minorant failed for minorant "<<i<<std::endl;
	  }
	  err++;
	}
	(*precomputed)[bun[i]]=tmpm;
	bun[i]=tmpm;
      }
      else {
	bun[i]=mapit->second;
      }
    }

    return err;
    
  }

  // *****************************************************************************
//                                pop_aft
// *****************************************************************************

  int UQPModelBlock::pop_aft()
  {
    assert(constant_minorant.size()==bundle.size());
    if (bundle.size()<=1)
      return 1;
    bundle.pop_back();
    constant_minorant.pop_back();
    return 0;
  }
    
UQPModelPointer::~UQPModelPointer()
  {}

  QPSumModelDataObject* UQPModelPointer::generate_summodel_data(BundleModel* )
  { return new UQPSumModelBlock(this); }
  
  QPConeModelDataObject*  UQPModelPointer::generate_conemodel_data(BundleModel*)
  { return new UQPConeModelBlock(this); }
 
}
