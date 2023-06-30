/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/CFunction.cxx
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



#include "CFunction.hxx"
#include "GroundsetModification.hxx"
 
using namespace CH_Matrix_Classes;

namespace ConicBundle {

  CFunction::CFunction(void* fk,cb_functionp fp,cb_subgextp se,int prdim)
  {
    function_key=fk;
    assert(fp!=0);
    oracle=fp;
    subgext=se;
    primaldim=prdim;
    max_new=1;
  }
       
  int CFunction::evaluate(
			  const  Matrix& current_point,
			  double relprec,
			  double&  obval,
			  std::vector<Minorant*>&  mnrt,
			  PrimalExtender*& 
			)
  {
    assert(mnrt.size()==0);
    Integer dim=current_point.dim();
    Matrix subg(dim,Integer(max_new)); chk_set_init(subg,1); 
    Matrix val(max_new,1);             chk_set_init(val,1);
    
    int n_new;
    int ret_code;
    if (primaldim>0){
      Matrix x(primaldim,max_new);  chk_set_init(x,1);
      ret_code =(*oracle)(
			  function_key,
			  const_cast<double *>(current_point.get_store()),
			  relprec,
			  max_new,
			  &obval,
			  &n_new,
			  val.get_store(),
			  subg.get_store(),
			  x.get_store()
			  );
      for (Integer i=0;i<n_new;i++){
	mnrt.push_back(new Minorant(false,val(i),dim,subg.get_store()+i*dim,0,1.,new PrimalMatrix(x.col(i))));
      }
    }
    else {
      ret_code =(*oracle)(
			  function_key,
			  const_cast<double *>(current_point.get_store()),
			  relprec,
			  max_new,
			  &obval,
			  &n_new,
			  val.get_store(),
			  subg.get_store(),
			  0
			  );
      for (Integer i=0;i<n_new;i++){
	mnrt.push_back(new Minorant(false,val(i),dim,subg.get_store()+i*dim));
      }
    }
    
    return ret_code;

  }

  int CFunction::apply_modification
  (
   const OracleModification& omod,
   const CH_Matrix_Classes::Matrix* new_center ,
   const CH_Matrix_Classes::Matrix* old_center ,
   bool& discard_obj_in_center ,
   bool& discard_model , 
   bool& discard_aggregates ,
   MinorantExtender*& extender 
   )
  {
    int err=0;

    const GroundsetModification* modp=dynamic_cast<const GroundsetModification*>(&omod);
    if (modp==0) {
      err++;
      if (cb_out()){
	get_out()<<"**** ERROR in CFunction::apply_modification(...): the oraclemodification is no GroundsetModification"<<std::endl;
	}

      discard_obj_in_center=true;
      discard_model=true;
      discard_aggregates=true;
      extender=0;
    }
    
    discard_obj_in_center=false;
    discard_model=false;
    discard_aggregates=false;
    extender=0;

    if (modp->no_modification()){
      return err;
    }

    //--- check whether deletions concern only zero variables/zero matrices
    if ((old_center==0)||(new_center==0)) {
      discard_obj_in_center=true;
    }
    else {
      if((!modp->deleted_variables_are_zero(*old_center))||
	 (!modp->mapped_variables_are_equal(*new_center,*old_center))){
	discard_obj_in_center=true;
      }
    }
    
    //--- check whether additions have zero variables and require/allow extensions
    if (!discard_obj_in_center){
      if(!modp->new_variables_are_zero(*new_center)){
	discard_obj_in_center=true;
      }
    }
  
    //provide an extender if needed
    if ((!discard_model)&&(!discard_aggregates)&&(subgext)){
      extender=new CFunctionMinorantExtender(function_key,subgext);
    }

    return err;
  }


    int CFunctionMinorantExtender::extend(Minorant& minorant,int n_coords,const int* indices)
  { 
    if (subgext==0) return 1;
    int ret_code;
    Matrix new_subgradient_values(n_coords,1,0.);
    PrimalData* generating_primal=minorant.get_primal();
    if (generating_primal==0){
      ret_code = subgext(
			 function_key,
			 0,
			 n_coords,
			 const_cast<Integer *>(indices),
			 new_subgradient_values.get_store()
			 );
    }
    else {
      const PrimalMatrix* xp=dynamic_cast<const PrimalMatrix *>(generating_primal);
      assert(xp!=0);
      ret_code = subgext(
			 function_key,
			 const_cast<double *>(xp->get_store()),
			 n_coords,
			 const_cast<Integer *>(indices),
			 new_subgradient_values.get_store()
			 );
    }
    for (int i=0;i<n_coords;i++){
      if (new_subgradient_values(i)!=0.)
	minorant.add_coeff(indices[i],new_subgradient_values(i));
    }
     
    return ret_code;
       
  }

}

