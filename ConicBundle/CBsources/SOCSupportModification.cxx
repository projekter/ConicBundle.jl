/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/SOCSupportModification.cxx
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



#include "SOCSupportModification.hxx"

using namespace CH_Matrix_Classes;

namespace ConicBundle {


  SOCSupportModification::~SOCSupportModification()
  {}

  /// reassign the variables as given in @a map_to_old, calls Modification::add_reassign_vars
  int SOCSupportModification::add_reassign_vars(const CH_Matrix_Classes::Indexmatrix& map_to_old)
  {
    if (map_to_old.dim()==0) {
      if (cb_out())
	get_out()<<"**** ERROR in SOCSupportModification::add_reassign_vars(.): map_to_old is empty but must preserve element 0 because of its special role in the second order cone"<<std::endl;
      return 1;
    }
    if (map_to_old.dim()==0) {
      if (cb_out())
	get_out()<<"**** ERROR in SOCSupportModification::add_reassign_vars(.): map_to_old(0)="<<map_to_old(0)<<" != 0, but it has preserve element 0 due to its special role in the second order cone"<<std::endl;
      return 1;
    }
    return mdf.add_reassign_vars(map_to_old);
  }

  /// delete the variables indexed by @a del_ind; for each new index @a map_to_old returns the old one; calls Modification::add_delete_vars
  int SOCSupportModification::add_delete_vars(const CH_Matrix_Classes::Indexmatrix& del_ind,
		      CH_Matrix_Classes::Indexmatrix& map_to_old)
  {
    if (min(del_ind)==0){
      if (cb_out())
	get_out()<<"**** ERROR in SOCSupportModification::add_delete_vars(..): del_ind contains 0,  but element 0 may not be deleted due to its special role in the second order cone"<<std::endl;
      return 1;
    }
    return mdf.add_delete_vars(del_ind,map_to_old);
  }

  /// incorporate the OracleModification @a m (it should only contain variable changes, but this is not checked!) into this one; calls Modification::incorporate
  int SOCSupportModification::incorporate(const OracleModification& m)
  {
    const SOCSupportModification* gm=dynamic_cast<const SOCSupportModification*>(&m);
    if (gm==0)
      return 1;
    if ((gm->map_to_old_variables()!=0)&&
	((gm->map_to_old_variables()->dim()==0)||
	 ((*(gm->map_to_old_variables()))(0)!=0))
	){
      if (cb_out())
	get_out()<<"**** ERROR in SOCSupportModification::incorporate(.): map_to_old has to preserve element 0 due to its special role in the second order cone but does not"<<std::endl;
      return 1;
    }
				      
    return mdf.incorporate(gm->mdf);
  }


  
}

