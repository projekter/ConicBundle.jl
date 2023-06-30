/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/Modification.cxx
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



#include <typeinfo>
#include "Modification.hxx"

using namespace CH_Matrix_Classes;

namespace ConicBundle {

  // *****************************************************************************
  //                              Modification::Modification
  // *****************************************************************************

  Modification::Modification(Integer volddim, Integer rolddim,
    const CBout* cb,
    int incr,
    bool ensure_start_val_box_feasibility,
    bool ensure_bounds_consistency,
    Real start_val_def,
    Real lb_def,
    Real ub_def,
    Real rhslb_def,
    Real rhsub_def,
    Real cost_def) :
    ModificationBase(cb, incr) {
    append_only = false;

    var_set_lb = 0;
    var_set_ub = 0;
    var_newdim = var_olddim = 0;
    var_append_dim = 0;
    var_append_cols = 0;
    var_append_lb = 0;
    var_append_ub = 0;
    var_start_val = 0;
    var_append_costs = 0;
    var_del_ind = 0;
    var_map_to_old = 0;
    var_new_ind = 0;

    row_set_rhslb = 0;
    row_set_rhsub = 0;
    row_newdim = row_olddim;
    row_append_dim = 0;
    row_append_mat = 0;
    row_append_rhslb = 0;
    row_append_rhsub = 0;
    row_del_ind = 0;
    row_map_to_old = 0;
    row_new_ind = 0;

    clear(volddim, rolddim, ensure_start_val_box_feasibility, ensure_bounds_consistency, start_val_def,
      lb_def, ub_def, rhslb_def, rhsub_def, cost_def);
  }

  // *****************************************************************************
  //                              Modification::~Modification
  // *****************************************************************************

  Modification::~Modification() {
    clear(0, 0);
  }

  // *****************************************************************************
  //                              Modification::Modification
  // *****************************************************************************

  int Modification::clear(Integer volddim, Integer rolddim,
    bool ensure_start_val_box_feasibility,
    bool ensure_bounds_consistency,
    Real start_val_def,
    Real lb_def,
    Real ub_def,
    Real rhslb_def,
    Real rhsub_def,
    Real cost_def) {
    append_only = false;
    var_olddim = max(0, volddim);
    row_olddim = max(0, rolddim);

    enforce_bounds_consistency = ensure_bounds_consistency;
    enforce_start_val_box_feasibility = ensure_start_val_box_feasibility;
    bounds_minus_infinity = lb_def;
    bounds_plus_infinity = ub_def;
    int err = 0;
    if (bounds_minus_infinity > bounds_plus_infinity) {
      if (cb_out()) {
        get_out() << "**** ERROR: Modification::clear(): default lower bound=" << bounds_minus_infinity << " greater than default upper bound=" << bounds_plus_infinity << std::endl;
      }
      err++;
    }
    rhs_minus_infinity = rhslb_def;
    rhs_plus_infinity = rhsub_def;
    if (rhs_minus_infinity > rhs_plus_infinity) {
      if (cb_out()) {
        get_out() << "**** ERROR: Modification::clear(): default right hand side lower bound=" << rhs_minus_infinity << " greater than default right hand side upper bound=" << rhs_plus_infinity << std::endl;
      }
      err++;
    }
    start_val_default = start_val_def;
    if ((enforce_start_val_box_feasibility) &&
      ((start_val_default < bounds_minus_infinity) || (start_val_default > bounds_plus_infinity))
      ) {
      if (cb_out()) {
        get_out() << "**** ERROR: Modification::clear(): default start value=" << start_val_default << " outside default bounds [" << bounds_minus_infinity << "," << bounds_plus_infinity << "]" << std::endl;
      }
      err++;
      if (start_val_default < bounds_minus_infinity)
        start_val_default = bounds_minus_infinity;
      if (start_val_default > bounds_plus_infinity)
        start_val_default = bounds_plus_infinity;
    }
    cost_default = cost_def;

    var_newdim = var_olddim;
    delete var_set_lb;
    var_set_lb = 0;
    delete var_set_ub;
    var_set_ub = 0;
    var_append_dim = 0;
    delete var_append_cols;
    var_append_cols = 0;
    delete var_append_lb;
    var_append_lb = 0;
    delete var_append_ub;
    var_append_ub = 0;
    delete var_start_val;
    var_start_val = 0;
    delete var_append_costs;
    var_append_costs = 0;
    delete var_del_ind;
    var_del_ind = 0;
    delete var_map_to_old;
    var_map_to_old = 0;
    delete var_new_ind;
    var_new_ind = 0;

    row_newdim = row_olddim;
    delete row_set_rhslb;
    row_set_rhslb = 0;
    delete row_set_rhsub;
    row_set_rhsub = 0;
    row_append_dim = 0;
    delete row_append_mat;
    row_append_mat = 0;
    delete row_append_rhslb;
    row_append_rhslb = 0;
    delete row_append_rhsub;
    row_append_rhsub = 0;
    delete row_del_ind;
    row_del_ind = 0;
    delete row_map_to_old;
    row_map_to_old = 0;
    delete row_new_ind;
    row_new_ind = 0;

    return err;
  }



  // *****************************************************************************
  //                           Modification::add_set_lb
  // *****************************************************************************

  int Modification::add_set_lb(Integer ind, Real val) {
    int err = 0;
    if (append_only) {
      if (cb_out())
        get_out() << "**** ERROR: Modification::add_set_lb(..): append_to_old is set to true, so this operations is not allowed" << std::endl;
      err++;
    }
    if ((ind < 0) || (ind >= var_newdim)) {
      if (cb_out())
        get_out() << "**** ERROR: Modification::add_set_lb(..): index=" << ind << " exceeds the range, it must lie between 0 and " << var_newdim - 1 << std::endl;
      err++;
    }
    if (val > bounds_plus_infinity) {
      if (cb_out()) get_out() << "**** ERROR: Modification::add_set_lb(..): lower bound exceeds plus_infinity: " << val << std::endl;
      err++;
    }
    if (val == bounds_plus_infinity) {
      if (cb_out()) get_out() << "**** WARNING: Modification::add_set_lb(...): lower bound equals plus_infinity: " << val << std::endl;
    }
    if (val < bounds_minus_infinity) {
      if (cb_out()) get_out() << "**** WARNING: Modification::add_set_lb(...): lower bound is smaller than minus_infinity: " << val << std::endl;
    }
    if (err) {
      return 1;
    }

    if (var_map_to_old)
      ind = (*var_map_to_old)(ind);

    if (ind < var_olddim) {
      if (var_set_lb == 0) {
        var_set_lb = new Realmap;
      }
      (*var_set_lb)[ind] = val;
      assert((var_set_lb->size() > 0) && (var_set_lb->size() <= (unsigned long)(var_olddim)));
    } else {
      assert(var_append_dim > ind - var_olddim);
      if (var_append_lb == 0) {
        var_append_lb = new Matrix(var_append_dim, 1, bounds_minus_infinity);
      }
      (*var_append_lb)(ind - var_olddim) = val;
    }

    return 0;
  }


  // *****************************************************************************
  //                           Modification::add_set_ub
  // *****************************************************************************

  int Modification::add_set_ub(Integer ind, Real val) {
    int err = 0;
    if (append_only) {
      if (cb_out())
        get_out() << "**** ERROR: Modification::add_set_ub(..): append_to_old is set to true, so this operations is not allowed" << std::endl;
      err++;
    }
    if ((ind < 0) || (ind >= var_newdim)) {
      if (cb_out())
        get_out() << "**** ERROR: Modification::add_set_ub(..): index=" << ind << " exceeds the range, it must lie between 0 and " << var_newdim - 1 << std::endl;
      err++;
    }
    if (val < bounds_minus_infinity) {
      if (cb_out()) get_out() << "**** ERROR: Modification::add_set_ub(..): upper bound is smaller than minus_infinity: " << val << std::endl;
      err++;
    }
    if (val == bounds_minus_infinity) {
      if (cb_out()) get_out() << "**** WARNING: Modification::add_set_ub(...): upper bound equals minus_infinity: " << val << std::endl;
    }
    if (val > bounds_plus_infinity) {
      if (cb_out()) get_out() << "**** WARNING: Modification::add_set_ub(...): upper bound exceeds plus_infinity: " << val << std::endl;
    }
    if (err) {
      return 1;
    }

    if (var_map_to_old)
      ind = (*var_map_to_old)(ind);

    if (ind < var_olddim) {
      if (var_set_ub == 0) {
        var_set_ub = new Realmap;
      }
      (*var_set_ub)[ind] = val;
      assert((var_set_ub->size() > 0) && (var_set_ub->size() <= (unsigned long)(var_olddim)));
    } else {
      assert(var_append_dim > ind - var_olddim);
      if (var_append_ub == 0) {
        var_append_ub = new Matrix(var_append_dim, 1, bounds_plus_infinity);
      }
      (*var_append_ub)(ind - var_olddim) = val;
    }

    return 0;
  }



  // *****************************************************************************
  //                           Modification::add_append_vars
  // *****************************************************************************

  int Modification::add_append_vars(Integer append_dim,
    const Matrix* append_lb,
    const Matrix* append_ub,
    const Sparsemat* append_cols,
    const Matrix* start_val,
    const Matrix* append_costs) {
    int err = 0;
    if (append_dim < 0) {
      err++;
      if (cb_out())
        get_out() << "**** ERROR: Modification::add_append_variables(...): append_dim=" << append_dim << ", cannot append negative number of variables" << std::endl;
    }
    if ((append_lb != 0) && ((append_lb->coldim() != 1) || (append_lb->rowdim() != append_dim))) {
      err++;
      if (cb_out())
        get_out() << "**** ERROR: Modification::add_append_variables(...): vector of lower bounds has size " << append_lb->rowdim() << " x " << append_lb->coldim() << " but should have size " << append_dim << " x 1" << std::endl;
    }
    if ((append_ub != 0) && ((append_ub->coldim() != 1) || (append_ub->rowdim() != append_dim))) {
      err++;
      if (cb_out())
        get_out() << "**** ERROR: Modification::add_append_variables(...): vector of upper bounds has size " << append_ub->rowdim() << " x " << append_ub->coldim() << " but should have size " << append_dim << " x 1" << std::endl;
    }
    if (
      (append_cols != 0) &&
      (
        (append_cols->coldim() != append_dim) ||
        (append_cols->rowdim() != row_newdim)
        )
      ) {
      err++;
      if (cb_out())
        get_out() << "**** ERROR: Modification::add_append_variables(...): matrix to be appended to constraints has size " << append_cols->rowdim() << " x " << append_cols->coldim() << " but should have size " << row_newdim << " x " << append_dim << std::endl;
    }
    /*
    if (
        (add_to_rhs!=0)&&
        (
         (add_to_rhs->coldim()!=1)||
         ((dim>0)&&(add_to_rhs->rowdim()!=rhs.rowdim()))
    )
        ){
      err++;
      if (cb_out())
        get_out()<<"**** ERROR: Modification::add_append_variables(...): vector to be added to rhs has size "<<add_to_rhs->rowdim()<<" x "<<add_to_rhs->coldim()<<" but should have size "<<rhs.rowdim()<<" x 1"<<std::endl;
    }
    */
    if ((start_val != 0) && ((start_val->coldim() != 1) || (start_val->rowdim() != append_dim))) {
      err++;
      if (cb_out())
        get_out() << "**** ERROR: Modification::add_append_variables(...): vector of starting values has size " << start_val->rowdim() << " x " << start_val->coldim() << " but should have size " << append_dim << " x 1" << std::endl;
    }
    if ((append_costs != 0) && ((append_costs->coldim() != 1) || (append_costs->rowdim() != append_dim))) {
      err++;
      if (cb_out())
        get_out() << "**** ERROR: Modification::add_append_variables(...): vector of linear costs has size " << append_costs->rowdim() << " x " << append_costs->coldim() << " but should have size " << append_dim << " x 1" << std::endl;
    }

    //check correctnes of values of append_lb, append_ub, and start_val
    if ((append_lb != 0) || (append_ub != 0)) {
      for (Integer i = 0; i < append_dim; i++) {
        if (append_lb) {
          if ((*append_lb)[i] > bounds_plus_infinity) {
            if (cb_out()) get_out() << "**** ERROR: Modification::add_append_variables(...): lower bound of coordinate " << i << " exceeds plus_infinity: " << (*append_lb)[i] << std::endl;
            err++;
          }
          if ((*append_lb)[i] == bounds_plus_infinity) {
            if (cb_out()) get_out() << "**** WARNING: Modification::add_append_variables(...): lower bound of coordinate " << i << " equals plus_infinity: " << (*append_lb)[i] << std::endl;
          }
          if ((*append_lb)[i] < bounds_minus_infinity) {
            if (cb_out()) get_out() << "**** WARNING: Modification::add_append_variables(...): lower bound of coordinate " << i << " is smaller than minus_infinity: " << (*append_lb)[i] << std::endl;
          }
          if ((start_val) && ((*start_val)[i] < (*append_lb)[i])) {
            if (enforce_start_val_box_feasibility) {
              if (cb_out()) get_out() << "**** ERROR: Modification::add_append_variables(...): lower_bound[" << i << "]=" << (*append_lb)[i] << " is greater than start_val[" << i << "]=" << (*start_val)[i] << std::endl;
              err++;
            } else {
              if (cb_out()) get_out() << "**** WARNING: Modification::add_append_variables(...): lower_bound[" << i << "]=" << (*append_lb)[i] << " is greater than start_val[" << i << "]=" << (*start_val)[i] << std::endl;
            }
          }
        }
        if (append_ub) {
          if ((*append_ub)[i] < bounds_minus_infinity) {
            if (cb_out()) get_out() << "**** ERROR: Modification::add_append_variables(...): upper bound of coordinate " << i << " is smaller than minus_infinity: " << (*append_ub)[i] << std::endl;
            err++;
          }
          if ((*append_ub)[i] == bounds_minus_infinity) {
            if (cb_out()) get_out() << "**** WARNING: Modification::add_append_variables(...): upper bound of coordinate " << i << " equals minus_infinity: " << (*append_ub)[i] << std::endl;
          }
          if ((*append_ub)[i] > bounds_plus_infinity) {
            if (cb_out()) get_out() << "**** WARNING: Modification::add_append_variables(...): upper bound of coordinate " << i << " exceeds plus_infinity: " << (*append_ub)[i] << std::endl;
          }
          if ((append_lb) && ((*append_ub)[i] < (*append_lb)[i])) {
            if (cb_out()) get_out() << "**** ERROR: Modification::add_append_variables(...): lower bound " << (*append_lb)[i] << " of coordinate " << i << " exceeds upper bound " << (*append_ub)[i] << std::endl;
            err++;
          }
          if ((start_val) && ((*start_val)[i] > (*append_ub)[i])) {
            if (enforce_start_val_box_feasibility) {
              if (cb_out()) get_out() << "**** ERROR: Modification::add_append_variables(...): upper_bound[" << i << "]=" << (*append_ub)[i] << " is greater than start_val[" << i << "]=" << (*start_val)[i] << std::endl;
              err++;
            } else {
              if (cb_out()) get_out() << "**** WARNING: Modification::add_append_variables(...): upper_bound[" << i << "]=" << (*append_ub)[i] << " is greater than start_val[" << i << "]=" << (*start_val)[i] << std::endl;
            }
          }
        } //endif append_ub
      } //endfor
    } //endif append_lb or append_ub

    if (err) {
      return err;
    }

    if (append_dim == 0) {
      return 0;
    }

    var_append_dim += append_dim;
    var_newdim += append_dim;

    if (var_new_ind == 0)
      var_new_ind = new Indexmatrix(0, 1, Integer(0));
    var_new_ind->concat_below(Range(var_newdim - append_dim, var_newdim - 1));

    if (append_lb) {
      if (var_append_lb == 0)
        var_append_lb = new Matrix(var_append_dim - append_dim, 1, bounds_minus_infinity);
      var_append_lb->concat_below(*append_lb);
    } else {
      if (var_append_lb) {
        var_append_lb->concat_below(Matrix(append_dim, 1, bounds_minus_infinity));
      }
    }

    if (append_ub) {
      if (var_append_ub == 0)
        var_append_ub = new Matrix(var_append_dim - append_dim, 1, bounds_plus_infinity);
      var_append_ub->concat_below(*append_ub);
    } else {
      if (var_append_ub) {
        var_append_ub->concat_below(Matrix(append_dim, 1, bounds_plus_infinity));
      }
    }

    if (append_costs) {
      if (var_append_costs == 0)
        var_append_costs = new Matrix(var_append_dim - append_dim, 1, cost_default);
      var_append_costs->concat_below(*append_costs);
    } else {
      if (var_append_costs) {
        var_append_costs->concat_below(Matrix(append_dim, 1, cost_default));
      }
    }

    if (start_val) {
      if (var_start_val == 0) {
        var_start_val = new Matrix(var_append_dim - append_dim, 1, start_val_default);
        // for (Integer i=0;i<var_start_val->rowdim();i++){
        // 	if ((var_append_lb)&&((*var_append_lb)(i)>start_val_default)){
        // 	  (*var_start_val)(i)=(*var_append_lb)(i);
        // 	  continue;
        // 	}
        // 	if ((var_append_ub)&&((*var_append_ub)(i)<start_val_default)){
        // 	  (*var_start_val)(i)=(*var_append_ub)(i);
        // 	  continue;
        // 	}
        // }
      }
      var_start_val->concat_below(*start_val);
    } else {
      if (var_start_val) {
        var_start_val->concat_below(Matrix(append_dim, 1, start_val_default));
        // for (Integer i=var_append_dim-append_dim;i<var_start_val->rowdim();i++){
        // 	if ((var_append_lb)&&((*var_append_lb)(i)>start_val_default)){
        // 	  (*var_start_val)(i)=(*var_append_lb)(i);
        // 	  continue;
        // 	}
        // 	if ((var_append_ub)&&((*var_append_ub)(i)<start_val_default)){
        // 	  (*var_start_val)(i)=(*var_append_ub)(i);
        // 	  continue;
        // 	}
        // }
      }
    }

    if (var_map_to_old)
      var_map_to_old->concat_below(Range(var_olddim + var_append_dim - append_dim, var_olddim + var_append_dim - 1));

    if (append_cols) {
      if (var_append_cols == 0)
        var_append_cols = new Sparsemat(row_olddim, var_append_dim - append_dim);

      if (row_map_to_old == 0) {
        if (row_append_dim == 0) {
          var_append_cols->concat_right(*append_cols);
          if (row_append_mat)
            row_append_mat->concat_right(Sparsemat(0, append_dim));
        } else {
          var_append_cols->concat_right(append_cols->rows(Range(0, row_olddim - 1)));
          if (row_append_mat == 0)
            row_append_mat = new Sparsemat(row_append_dim, var_newdim - append_dim);
          row_append_mat->concat_right(append_cols->rows(Range(row_olddim, row_newdim - 1)));
        }
      } else { //row_map_to_old!=0; //need to transform indices back
        Indexmatrix indi1;
        Indexmatrix indj1;
        Matrix val1;
        append_cols->get_edge_rep(indi1, indj1, val1);
        Integer nz1 = indi1.dim();
        Integer nz2 = 0;
        for (Integer i = 0; i < nz1; i++) {
          Integer ind = (*row_map_to_old)(indi1(i));
          indi1(i) = ind;
          if (ind >= row_olddim)
            nz2++;
        }
        if (nz2 == 0) {
          var_append_cols->concat_right(Sparsemat(row_olddim, append_dim, nz1, indi1, indj1, val1));
          if (row_append_mat)
            row_append_mat->concat_right(Sparsemat(row_append_dim, append_dim));
        } else { //split the matrix in two parts
          Indexmatrix indi2(nz2, 1, Integer(0));
          Indexmatrix indj2(nz2, 1, Integer(0));
          Matrix val2(nz2, 1, 0.);
          Integer cnt1 = 0;
          Integer cnt2 = 0;
          for (Integer i = 0; i < nz1; i++) {
            Integer ind = indi1(i);
            if (ind < row_olddim) {
              indi1(cnt1) = indi1(i);
              indj1(cnt1) = indj1(i);
              val1(cnt1) = val1(i);
              cnt1++;
            } else {
              indi2(cnt2) = indi1(i) - row_olddim;
              indj2(cnt2) = indj1(i);
              val2(cnt2) = val1(i);
              cnt2++;
            }
          }
          assert(cnt1 + cnt2 == nz1);
          var_append_cols->concat_right(Sparsemat(row_olddim, append_dim, cnt1, indi1, indj1, val1));
          if (row_append_mat == 0) {
            indj2 += var_newdim - append_dim;
            row_append_mat = new Sparsemat(row_append_dim, var_newdim, cnt2, indi2, indj2, val2);
          } else {
            row_append_mat->concat_right(Sparsemat(row_append_dim, append_dim, cnt2, indi2, indj2, val2));
          }
        } //end of split the matrix in two
      } //end of else //row_map_to_old!=0
    } //end of if append_cols != 0	     
    else {
      if (var_append_cols) {
        var_append_cols->concat_right(Sparsemat(row_olddim, append_dim));
      }
      if (row_append_mat)
        row_append_mat->concat_right(Sparsemat(row_append_dim, append_dim));
    }


    assert((var_append_lb == 0) || (var_append_lb->rowdim() == var_append_dim));
    assert((var_append_ub == 0) || (var_append_ub->rowdim() == var_append_dim));
    assert((var_append_costs == 0) || (var_append_costs->rowdim() == var_append_dim));
    assert((var_start_val == 0) || (var_start_val->rowdim() == var_append_dim));
    assert((var_append_cols == 0) || ((var_append_cols->rowdim() == row_olddim) && (var_append_cols->coldim() == var_append_dim)));
    assert((var_map_to_old == 0) || (var_map_to_old->coldim() == 1));
    assert((var_map_to_old == 0) || (var_map_to_old->rowdim() + var_del_ind->rowdim() == var_olddim + var_append_dim));
    assert((row_append_rhslb == 0) || (row_append_rhslb->rowdim() == row_append_dim));
    assert((row_append_rhsub == 0) || (row_append_rhsub->rowdim() == row_append_dim));
    assert((row_append_mat == 0) || ((row_append_mat->coldim() == var_newdim) && (row_append_mat->rowdim() == row_append_dim)));
    assert((row_map_to_old == 0) || (row_map_to_old->coldim() == 1));
    assert((row_map_to_old == 0) || (row_map_to_old->rowdim() + row_del_ind->rowdim() == row_olddim + row_append_dim));


    return 0;
  }

  // *****************************************************************************
  //                              Modification::add_reassign_variables
  // *****************************************************************************

  int Modification::add_reassign_vars(const Indexmatrix& map_to_old) {
    int err = 0;
    if (append_only) {
      if (cb_out())
        get_out() << "**** ERROR: Modification::add_reassign_vars(..): append_to_old is set to true, so this operations is not allowed" << std::endl;
      err++;
      return err;
    }

    if (map_to_old.dim() == 0) {
      var_newdim = 0;
      var_append_dim = 0;

      if (var_map_to_old == 0)
        var_map_to_old = new Indexmatrix(0, 1, Integer(0));
      else
        var_map_to_old->init(0, 1, Integer(0));
      if (var_del_ind == 0)
        var_del_ind = new Indexmatrix(Range(0, var_olddim - 1));
      else
        var_del_ind->init(Range(0, var_olddim - 1));

      delete var_new_ind;
      var_new_ind = 0;

      delete var_set_lb;
      var_set_lb = 0;
      delete var_set_ub;
      var_set_ub = 0;
      delete var_append_cols;
      var_append_cols = 0;
      delete var_append_lb;
      var_append_lb = 0;
      delete var_append_ub;
      var_append_ub = 0;
      delete var_start_val;
      var_start_val = 0;
      delete var_append_costs;
      var_append_costs = 0;

      //row_newdim does not change
      //row_append_dim does not change
      delete row_append_mat;
      row_append_mat = 0;
      //row_append_rhs does not change
      //row_del_ind does not change
      //row_map_to_old does not change;

      assert((var_append_lb == 0) || (var_append_lb->rowdim() == var_append_dim));
      assert((var_append_ub == 0) || (var_append_ub->rowdim() == var_append_dim));
      assert((var_append_costs == 0) || (var_append_costs->rowdim() == var_append_dim));
      assert((var_start_val == 0) || (var_start_val->rowdim() == var_append_dim));
      assert((var_append_cols == 0) || ((var_append_cols->rowdim() == row_olddim) && (var_append_cols->coldim() == var_append_dim)));
      assert((var_map_to_old == 0) || (var_map_to_old->coldim() == 1));
      assert((var_map_to_old == 0) || (var_map_to_old->rowdim() + var_del_ind->rowdim() == var_olddim + var_append_dim));
      assert((row_append_rhslb == 0) || (row_append_rhslb->rowdim() == row_append_dim));
      assert((row_append_rhsub == 0) || (row_append_rhsub->rowdim() == row_append_dim));
      assert((row_append_mat == 0) || ((row_append_mat->coldim() == var_newdim) && (row_append_mat->rowdim() == row_append_dim)));
      assert((row_map_to_old == 0) || (row_map_to_old->coldim() == 1));
      assert((row_map_to_old == 0) || (row_map_to_old->rowdim() + row_del_ind->rowdim() == row_olddim + row_append_dim));

      return 0;
    }

    if (row_append_mat) {
      (*row_append_mat) = row_append_mat->cols(map_to_old);
    }
    Indexmatrix append_del_ind;

    err = adapt_map_to_old(var_map_to_old, var_del_ind, var_new_ind, append_del_ind, map_to_old, var_append_dim, var_olddim, var_newdim);
    if (err) {
      if (cb_out())
        get_out() << "**** ERROR: Modification::add_reassign_variables(...): adapt_map_to_old(...) failed and returned " << err << std::endl;
      return err;
    }

    var_newdim = var_map_to_old->dim();
    if (var_append_lb)
      var_append_lb->delete_rows(append_del_ind);
    if (var_append_ub)
      var_append_ub->delete_rows(append_del_ind);
    if (var_append_cols)
      var_append_cols->delete_cols(append_del_ind);
    if (var_start_val)
      var_start_val->delete_rows(append_del_ind);
    if (var_append_costs)
      var_append_costs->delete_rows(append_del_ind);
    //never mind about deleting newly useless entries in var_set_lb and var_set_ub

    assert((var_append_lb == 0) || (var_append_lb->rowdim() == var_append_dim));
    assert((var_append_ub == 0) || (var_append_ub->rowdim() == var_append_dim));
    assert((var_append_costs == 0) || (var_append_costs->rowdim() == var_append_dim));
    assert((var_start_val == 0) || (var_start_val->rowdim() == var_append_dim));
    assert((var_append_cols == 0) || ((var_append_cols->rowdim() == row_olddim) && (var_append_cols->coldim() == var_append_dim)));
    assert((var_map_to_old == 0) || (var_map_to_old->coldim() == 1));
    assert((var_map_to_old == 0) || (var_map_to_old->rowdim() + var_del_ind->rowdim() == var_olddim + var_append_dim));
    assert((row_append_rhslb == 0) || (row_append_rhslb->rowdim() == row_append_dim));
    assert((row_append_rhsub == 0) || (row_append_rhsub->rowdim() == row_append_dim));
    assert((row_append_mat == 0) || ((row_append_mat->coldim() == var_newdim) && (row_append_mat->rowdim() == row_append_dim)));
    assert((row_map_to_old == 0) || (row_map_to_old->coldim() == 1));
    assert((row_map_to_old == 0) || (row_map_to_old->rowdim() + row_del_ind->rowdim() == row_olddim + row_append_dim));

    return 0;
  }

  // *****************************************************************************
  //                              Modification::add_delete_variables
  // *****************************************************************************

  int Modification::add_delete_vars(const Indexmatrix& delete_variables,
    Indexmatrix& map_to_old) {
    int err = 0;
    if (append_only) {
      if (cb_out())
        get_out() << "**** ERROR: Modification::add_delete_vars(..): append_to_old is set to true, so this operations is not allowed" << std::endl;
      err++;
      return err;
    }

    err = form_map_to_old(map_to_old, delete_variables, var_newdim);
    if (err) {
      if (cb_out())
        get_out() << "**** ERROR: Modification::add_delete_vars(...): form_map_to_old(...) failed and returned " << err << std::endl;
      return err;
    }

    err = add_reassign_vars(map_to_old);
    if (err) {
      if (cb_out())
        get_out() << "**** ERROR: Modification::add_delete_vars(...): add_reassign_vars(...) failed and returned " << err << std::endl;
    }

    assert((var_append_lb == 0) || (var_append_lb->rowdim() == var_append_dim));
    assert((var_append_ub == 0) || (var_append_ub->rowdim() == var_append_dim));
    assert((var_append_costs == 0) || (var_append_costs->rowdim() == var_append_dim));
    assert((var_start_val == 0) || (var_start_val->rowdim() == var_append_dim));
    assert((var_append_cols == 0) || ((var_append_cols->rowdim() == row_olddim) && (var_append_cols->coldim() == var_append_dim)));
    assert((var_map_to_old == 0) || (var_map_to_old->coldim() == 1));
    assert((var_map_to_old == 0) || (var_map_to_old->rowdim() + var_del_ind->rowdim() == var_olddim + var_append_dim));
    assert((row_append_rhslb == 0) || (row_append_rhslb->rowdim() == row_append_dim));
    assert((row_append_rhsub == 0) || (row_append_rhsub->rowdim() == row_append_dim));
    assert((row_append_mat == 0) || ((row_append_mat->coldim() == var_newdim) && (row_append_mat->rowdim() == row_append_dim)));
    assert((row_map_to_old == 0) || (row_map_to_old->coldim() == 1));
    assert((row_map_to_old == 0) || (row_map_to_old->rowdim() + row_del_ind->rowdim() == row_olddim + row_append_dim));



    return err;
  }

  // *****************************************************************************
  //                           Modification::add_set_rhslb
  // *****************************************************************************

  int Modification::add_set_rhslb(Integer ind, Real val) {
    int err = 0;
    if (append_only) {
      if (cb_out())
        get_out() << "**** ERROR: Modification::add_set_rhslb(..): append_to_old is set to true, so this operations is not allowed" << std::endl;
      err++;
    }
    if ((ind < 0) || (ind >= row_newdim)) {
      if (cb_out())
        get_out() << "**** ERROR: Modification::add_set_rhslb(..): index=" << ind << " exceeds the range, it must lie between 0 and " << row_newdim - 1 << std::endl;
      err++;
    }
    if (val > rhs_plus_infinity) {
      if (cb_out()) get_out() << "**** ERROR: Modification::add_set_rhslb(..): lower bound exceeds plus_infinity: " << val << std::endl;
      err++;
    }
    if (val == rhs_plus_infinity) {
      if (cb_out()) get_out() << "**** WARNING: Modification::add_set_rhslb(...): lower bound equals plus_infinity: " << val << std::endl;
    }
    if (val < rhs_minus_infinity) {
      if (cb_out()) get_out() << "**** WARNING: Modification::add_set_rhslb(...): lower bound is smaller than minus_infinity: " << val << std::endl;
    }
    if (err) {
      return 1;
    }

    if (row_map_to_old)
      ind = (*row_map_to_old)(ind);

    if (ind < row_olddim) {
      if (row_set_rhslb == 0) {
        row_set_rhslb = new Realmap;
      }
      (*row_set_rhslb)[ind] = val;
      assert((row_set_rhslb->size() > 0) && (row_set_rhslb->size() <= (unsigned long)(row_olddim)));
    } else {
      assert(row_append_dim > ind - row_olddim);
      if (row_append_rhslb == 0) {
        row_append_rhslb = new Matrix(row_append_dim, 1, rhs_minus_infinity);
      }
      (*row_append_rhslb)(ind - row_olddim) = val;
    }

    return 0;
  }


  // *****************************************************************************
  //                           Modification::add_set_rhsub
  // *****************************************************************************

  int Modification::add_set_rhsub(Integer ind, Real val) {
    int err = 0;
    if (append_only) {
      if (cb_out())
        get_out() << "**** ERROR: Modification::add_set_rhsub(..): append_to_old is set to true, so this operations is not allowed" << std::endl;
      err++;
    }
    if ((ind < 0) || (ind >= row_newdim)) {
      if (cb_out())
        get_out() << "**** ERROR: Modification::add_set_rhsub(..): index=" << ind << " exceeds the range, it must lie between 0 and " << row_newdim - 1 << std::endl;
      err++;
    }
    if (val < rhs_minus_infinity) {
      if (cb_out()) get_out() << "**** ERROR: Modification::add_set_rhsub(..): upper bound is smaller than minus_infinity: " << val << std::endl;
      err++;
    }
    if (val == rhs_minus_infinity) {
      if (cb_out()) get_out() << "**** WARNING: Modification::add_set_rhsub(...): upper bound equals minus_infinity: " << val << std::endl;
    }
    if (val > rhs_plus_infinity) {
      if (cb_out()) get_out() << "**** WARNING: Modification::add_set_rhsub(...): upper bound exceeds plus_infinity: " << val << std::endl;
    }
    if (err) {
      return 1;
    }

    if (row_map_to_old)
      ind = (*row_map_to_old)(ind);

    if (ind < row_olddim) {
      if (row_set_rhsub == 0) {
        row_set_rhsub = new Realmap;
      }
      (*row_set_rhsub)[ind] = val;
      assert((row_set_rhsub->size() > 0) && (row_set_rhsub->size() <= (unsigned long)(row_olddim)));
    } else {
      assert(row_append_dim > ind - row_olddim);
      if (row_append_rhsub == 0) {
        row_append_rhsub = new Matrix(row_append_dim, 1, rhs_plus_infinity);
      }
      (*row_append_rhsub)(ind - row_olddim) = val;
    }

    return 0;
  }




  // *****************************************************************************
  //                        Modification::add_append_rows
  // *****************************************************************************

  int Modification::add_append_rows(Integer append_dim,
    const Sparsemat* append_rows,
    const Matrix* append_rhslb,
    const Matrix* append_rhsub) {
    int err = 0;
    if ((append_rows) && (
      (append_rows->rowdim() != append_dim) ||
      (append_rows->coldim() != var_newdim)
      )) {
      err++;
      if (cb_out())
        get_out() << "**** ERROR: Modification::add_append_rows(...): matrix to be appended has size " << append_rows->rowdim() << " x " << append_rows->coldim() << " but should have size " << append_dim << " x " << var_newdim << std::endl;
    }
    if ((append_rhslb) && ((append_rhslb->coldim() != 1) || (append_rhslb->rowdim() != append_dim))) {
      err++;
      if (cb_out())
        get_out() << "**** ERROR: Modification::add_append_rows(...): vector to be appended to rhslb has size " << append_rhslb->rowdim() << " x " << append_rhslb->coldim() << " but should have size " << append_dim << " x 1" << std::endl;
    }
    if ((append_rhsub) && ((append_rhsub->coldim() != 1) || (append_rhsub->rowdim() != append_dim))) {
      err++;
      if (cb_out())
        get_out() << "**** ERROR: Modification::add_append_rows(...): vector to be appended to rhsub has size " << append_rhsub->rowdim() << " x " << append_rhsub->coldim() << " but should have size " << append_dim << " x 1" << std::endl;
    }

    if ((append_rhslb != 0) || (append_rhsub != 0)) {
      for (Integer i = 0; i < append_dim; i++) {
        if (append_rhslb) {
          if ((*append_rhslb)[i] > rhs_plus_infinity) {
            if (cb_out()) get_out() << "**** ERROR: Modification::add_append_rows(...): lower bound of row " << i << " exceeds plus_infinity: " << (*append_rhslb)[i] << std::endl;
            err++;
          }
          if ((*append_rhslb)[i] == rhs_plus_infinity) {
            if (cb_out()) get_out() << "**** WARNING: Modification::add_append_rows(...): lower bound of row " << i << " equals plus_infinity: " << (*append_rhslb)[i] << std::endl;
          }
          if ((*append_rhslb)[i] < rhs_minus_infinity) {
            if (cb_out()) get_out() << "**** WARNING: Modification::add_append_rows(...): lower bound of row " << i << " is smaller than minus_infinity: " << (*append_rhslb)[i] << std::endl;
          }
        }
        if (append_rhsub) {
          if ((*append_rhsub)[i] < rhs_minus_infinity) {
            if (cb_out()) get_out() << "**** ERROR: Modification::add_append_rows(...): upper bound of row " << i << " is smaller than minus_infinity: " << (*append_rhsub)[i] << std::endl;
            err++;
          }
          if ((*append_rhsub)[i] == rhs_minus_infinity) {
            if (cb_out()) get_out() << "**** WARNING: Modification::add_append_rows(...): upper bound of row " << i << " equals minus_infinity: " << (*append_rhsub)[i] << std::endl;
          }
          if ((*append_rhsub)[i] > rhs_plus_infinity) {
            if (cb_out()) get_out() << "**** WARNING: Modification::add_append_rows(...): upper bound of row " << i << " exceeds plus_infinity: " << (*append_rhsub)[i] << std::endl;
          }
          if ((append_rhslb) && ((*append_rhsub)[i] < (*append_rhslb)[i])) {
            if (cb_out()) get_out() << "**** ERROR: Modification::add_append_rows(...): lower bound " << (*append_rhslb)[i] << " of row " << i << " exceeds upper bound " << (*append_rhsub)[i] << std::endl;
            err++;
          }
        } //endif append_rhsub
      } //endfor
    } //endif append_rhslb or append_rhsub

    if (err)
      return 1;

    if (append_dim == 0) {
      return 0;
    }

    row_newdim += append_dim;
    row_append_dim += append_dim;

    if (row_new_ind == 0)
      row_new_ind = new Indexmatrix;
    row_new_ind->concat_below(Range(row_newdim - append_dim, row_newdim - 1));


    if (append_rows) {
      if (row_append_mat == 0) {
        row_append_mat = new Sparsemat(row_append_dim - append_dim, var_newdim);
      }
      row_append_mat->concat_below(*append_rows);
    } else {
      if (row_append_mat)
        row_append_mat->concat_below(Sparsemat(append_dim, var_newdim));
    }

    if (append_rhslb) {
      if (row_append_rhslb == 0) {
        row_append_rhslb = new Matrix(row_append_dim - append_dim, 1, rhs_minus_infinity);
      }
      row_append_rhslb->concat_below(*append_rhslb);
    } else {
      if (row_append_rhslb)
        row_append_rhslb->concat_below(Matrix(append_dim, 1, rhs_minus_infinity));
    }

    if (append_rhsub) {
      if (row_append_rhsub == 0) {
        row_append_rhsub = new Matrix(row_append_dim - append_dim, 1, rhs_plus_infinity);
      }
      row_append_rhsub->concat_below(*append_rhsub);
    } else {
      if (row_append_rhsub)
        row_append_rhsub->concat_below(Matrix(append_dim, 1, rhs_plus_infinity));
    }

    if (row_map_to_old)
      row_map_to_old->concat_below(Range(row_olddim + row_append_dim - append_dim, row_olddim + row_append_dim - 1));

    assert((var_append_lb == 0) || (var_append_lb->rowdim() == var_append_dim));
    assert((var_append_ub == 0) || (var_append_ub->rowdim() == var_append_dim));
    assert((var_append_costs == 0) || (var_append_costs->rowdim() == var_append_dim));
    assert((var_start_val == 0) || (var_start_val->rowdim() == var_append_dim));
    assert((var_append_cols == 0) || ((var_append_cols->rowdim() == row_olddim) && (var_append_cols->coldim() == var_append_dim)));
    assert((var_map_to_old == 0) || (var_map_to_old->coldim() == 1));
    assert((var_map_to_old == 0) || (var_map_to_old->rowdim() + var_del_ind->rowdim() == var_olddim + var_append_dim));
    assert((row_append_rhslb == 0) || (row_append_rhslb->rowdim() == row_append_dim));
    assert((row_append_rhsub == 0) || (row_append_rhsub->rowdim() == row_append_dim));
    assert((row_append_mat == 0) || ((row_append_mat->coldim() == var_newdim) && (row_append_mat->rowdim() == row_append_dim)));
    assert((row_map_to_old == 0) || (row_map_to_old->coldim() == 1));
    assert((row_map_to_old == 0) || (row_map_to_old->rowdim() + row_del_ind->rowdim() == row_olddim + row_append_dim));


    return 0;
  }


  // *****************************************************************************
  //                              Modification::add_reassign_rows
  // *****************************************************************************

  int Modification::add_reassign_rows(const Indexmatrix& map_to_old) {
    int err = 0;
    if (append_only) {
      if (cb_out())
        get_out() << "**** ERROR: Modification::add_reassign_rows(..): append_to_old is set to true, so this operations is not allowed" << std::endl;
      err++;
      return err;
    }
    if (map_to_old.dim() == 0) {
      row_newdim = 0;
      row_append_dim = 0;

      if (row_map_to_old == 0)
        row_map_to_old = new Indexmatrix(0, 1, Integer(0));
      else
        row_map_to_old->init(0, 1, Integer(0));
      if (row_del_ind == 0)
        row_del_ind = new Indexmatrix(Range(0, row_olddim - 1));
      else
        row_del_ind->init(Range(0, row_olddim - 1));

      delete row_set_rhslb;
      row_set_rhslb = 0;
      delete row_set_rhsub;
      row_set_rhsub = 0;
      delete var_append_cols;
      var_append_cols = 0;
      delete row_append_mat;
      row_append_mat = 0;
      delete row_append_rhslb;
      row_append_rhslb = 0;
      delete row_append_rhsub;
      row_append_rhsub = 0;
      delete row_new_ind;
      row_new_ind = 0;

      assert((var_append_lb == 0) || (var_append_lb->rowdim() == var_append_dim));
      assert((var_append_ub == 0) || (var_append_ub->rowdim() == var_append_dim));
      assert((var_append_costs == 0) || (var_append_costs->rowdim() == var_append_dim));
      assert((var_start_val == 0) || (var_start_val->rowdim() == var_append_dim));
      assert((var_append_cols == 0) || ((var_append_cols->rowdim() == row_olddim) && (var_append_cols->coldim() == var_append_dim)));
      assert((var_map_to_old == 0) || (var_map_to_old->coldim() == 1));
      assert((var_map_to_old == 0) || (var_map_to_old->rowdim() + var_del_ind->rowdim() == var_olddim + var_append_dim));
      assert((row_append_rhslb == 0) || (row_append_rhslb->rowdim() == row_append_dim));
      assert((row_append_rhsub == 0) || (row_append_rhsub->rowdim() == row_append_dim));
      assert((row_append_mat == 0) || ((row_append_mat->coldim() == var_newdim) && (row_append_mat->rowdim() == row_append_dim)));
      assert((row_map_to_old == 0) || (row_map_to_old->coldim() == 1));
      assert((row_map_to_old == 0) || (row_map_to_old->rowdim() + row_del_ind->rowdim() == row_olddim + row_append_dim));

      return 0;
    }

    Indexmatrix append_del_ind;

    err = adapt_map_to_old(row_map_to_old, row_del_ind, row_new_ind, append_del_ind, map_to_old, row_append_dim, row_olddim, row_newdim);
    if (err) {
      if (cb_out())
        get_out() << "**** ERROR: Modification::add_reassign_rows(...): adapt_map_to_old(...) failed and returned " << err << std::endl;
      return err;
    }


    row_newdim = row_map_to_old->dim();
    if (row_append_mat) {
      row_append_mat->delete_rows(append_del_ind);
    }
    if (row_append_rhslb) {
      row_append_rhslb->delete_rows(append_del_ind);
    }
    if (row_append_rhsub) {
      row_append_rhsub->delete_rows(append_del_ind);
    }

    assert((var_append_lb == 0) || (var_append_lb->rowdim() == var_append_dim));
    assert((var_append_ub == 0) || (var_append_ub->rowdim() == var_append_dim));
    assert((var_append_costs == 0) || (var_append_costs->rowdim() == var_append_dim));
    assert((var_start_val == 0) || (var_start_val->rowdim() == var_append_dim));
    assert((var_append_cols == 0) || ((var_append_cols->rowdim() == row_olddim) && (var_append_cols->coldim() == var_append_dim)));
    assert((var_map_to_old == 0) || (var_map_to_old->coldim() == 1));
    assert((var_map_to_old == 0) || (var_map_to_old->rowdim() + var_del_ind->rowdim() == var_olddim + var_append_dim));
    assert((row_append_rhslb == 0) || (row_append_rhslb->rowdim() == row_append_dim));
    assert((row_append_rhsub == 0) || (row_append_rhsub->rowdim() == row_append_dim));
    assert((row_append_mat == 0) || ((row_append_mat->coldim() == var_newdim) && (row_append_mat->rowdim() == row_append_dim)));
    assert((row_map_to_old == 0) || (row_map_to_old->coldim() == 1));
    assert((row_map_to_old == 0) || (row_map_to_old->rowdim() + row_del_ind->rowdim() == row_olddim + row_append_dim));


    return 0;
  }

  // *****************************************************************************
  //                              Modification::add_delete_rows
  // *****************************************************************************

  int Modification::add_delete_rows(const Indexmatrix& delete_rows,
    Indexmatrix& map_to_old) {
    int err = 0;
    if (append_only) {
      if (cb_out())
        get_out() << "**** ERROR: Modification::add_delete_rows(..): append_to_old is set to true, so this operations is not allowed" << std::endl;
      err++;
      return err;
    }
    err = form_map_to_old(map_to_old, delete_rows, row_newdim);
    if (err) {
      if (cb_out())
        get_out() << "**** ERROR: Modification::add_delete_rows(...): form_map_to_old(...) failed and returned " << err << std::endl;
      return err;
    }

    err = add_reassign_rows(map_to_old);
    if (err) {
      if (cb_out())
        get_out() << "**** ERROR: Modification::add_delete_rows(...): add_reassign_rows(...) failed and returned " << err << std::endl;
    }

    assert((var_append_lb == 0) || (var_append_lb->rowdim() == var_append_dim));
    assert((var_append_ub == 0) || (var_append_ub->rowdim() == var_append_dim));
    assert((var_append_costs == 0) || (var_append_costs->rowdim() == var_append_dim));
    assert((var_start_val == 0) || (var_start_val->rowdim() == var_append_dim));
    assert((var_append_cols == 0) || ((var_append_cols->rowdim() == row_olddim) && (var_append_cols->coldim() == var_append_dim)));
    assert((var_map_to_old == 0) || (var_map_to_old->coldim() == 1));
    assert((var_map_to_old == 0) || (var_map_to_old->rowdim() + var_del_ind->rowdim() == var_olddim + var_append_dim));
    assert((row_append_rhslb == 0) || (row_append_rhslb->rowdim() == row_append_dim));
    assert((row_append_rhsub == 0) || (row_append_rhsub->rowdim() == row_append_dim));
    assert((row_append_mat == 0) || ((row_append_mat->coldim() == var_newdim) && (row_append_mat->rowdim() == row_append_dim)));
    assert((row_map_to_old == 0) || (row_map_to_old->coldim() == 1));
    assert((row_map_to_old == 0) || (row_map_to_old->rowdim() + row_del_ind->rowdim() == row_olddim + row_append_dim));


    return err;
  }

  // *****************************************************************************
  //                              Modification::incorporate
  // *****************************************************************************

  int Modification::incorporate(const Modification& m) {
    int err = 0;
    if (m.append_to_old()) {
      if ((m.var_olddim != var_olddim) || (m.row_olddim != row_olddim)) {
        err++;
        if (cb_out())
          get_out() << "**** ERROR: Modification::incorporate(const Modfication& m): append_to_old is set to true; old variable and row dimension=(" << m.var_olddim << "," << m.row_olddim << ") of m must match the values (" << var_newdim << "," << var_olddim << ") of this but does not" << std::endl;
        return err;
      }
      //first append the columns
      Sparsemat* rowapp = 0;
      if (m.var_append_dim > 0) {

        if (var_append_lb) {
          if (m.var_append_lb)
            var_append_lb->concat_below(*(m.var_append_lb));
          else
            var_append_lb->enlarge_below(m.var_append_dim, m.bounds_minus_infinity);
        } else {
          if (m.var_append_lb) {
            var_append_lb = new Matrix(var_append_dim, 1, bounds_minus_infinity);
            var_append_lb->concat_below(*(m.var_append_lb));
          } else {
            if (bounds_minus_infinity != m.bounds_minus_infinity) {
              var_append_lb = new Matrix(var_append_dim, 1, bounds_minus_infinity);
              var_append_lb->enlarge_below(m.var_append_dim, m.bounds_minus_infinity);
            }
          }
        }

        if (var_append_ub) {
          if (m.var_append_ub)
            var_append_ub->concat_below(*(m.var_append_ub));
          else
            var_append_ub->enlarge_below(m.var_append_dim, m.bounds_plus_infinity);
        } else {
          if (m.var_append_ub) {
            var_append_ub = new Matrix(var_append_dim, 1, bounds_plus_infinity);
            var_append_ub->concat_below(*(m.var_append_ub));
          } else {
            if (bounds_plus_infinity != m.bounds_plus_infinity) {
              var_append_ub = new Matrix(var_append_dim, 1, bounds_plus_infinity);
              var_append_ub->enlarge_below(m.var_append_dim, m.bounds_plus_infinity);
            }
          }
        }

        if (var_append_cols) {
          if (m.var_append_cols)
            var_append_cols->concat_right(*(m.var_append_cols));
          else
            var_append_cols->concat_right(Sparsemat(var_append_cols->rowdim(), m.var_append_dim));
        } else {
          if (m.var_append_cols) {
            var_append_cols = new Sparsemat(row_olddim, var_append_dim);
            var_append_cols->concat_right(*(m.var_append_cols));
          }
        }

        if (var_start_val) {
          if (m.var_start_val)
            var_start_val->concat_below(*(m.var_start_val));
          else
            var_start_val->enlarge_below(m.var_append_dim, m.start_val_default);
        } else {
          if (m.var_start_val) {
            var_start_val = new Matrix(var_append_dim, 1, start_val_default);
            var_start_val->concat_below(*(m.var_start_val));
          } else {
            if (start_val_default != m.start_val_default) {
              var_start_val = new Matrix(var_append_dim, 1, start_val_default);
              var_start_val->enlarge_below(m.var_append_dim, m.start_val_default);
            }
          }
        }

        if (var_append_costs) {
          if (m.var_append_costs)
            var_append_costs->concat_below(*(m.var_append_costs));
          else
            var_append_costs->enlarge_below(m.var_append_dim, m.cost_default);
        } else {
          if (m.var_append_costs) {
            var_append_costs = new Matrix(var_append_dim, 1, cost_default);
            var_append_costs->concat_below(*(m.var_append_costs));
          } else {
            if (cost_default != m.cost_default) {
              var_append_costs = new Matrix(var_append_dim, 1, cost_default);
              var_append_costs->enlarge_below(m.var_append_dim, m.cost_default);
            }
          }
        }

        if (m.row_append_mat) {
          rowapp = new Sparsemat(m.row_append_mat->cols(Range(0, var_olddim - 1)));
          rowapp->concat_right(Sparsemat(rowapp->rowdim(), var_append_dim));
          rowapp->concat_right(m.row_append_mat->cols(Range(var_olddim, m.var_newdim - 1)));
        }


        if (var_map_to_old) {
          var_map_to_old->concat_below(Range(var_olddim + var_append_dim, var_olddim + var_append_dim + m.var_append_dim - 1));
          if (rowapp)
            (*rowapp) = rowapp->cols(*var_map_to_old);
        }

        if (var_new_ind) {
          var_new_ind->concat_below(Range(var_newdim, var_newdim + m.var_append_dim - 1));
        } else
          var_new_ind = new Indexmatrix(Range(var_newdim, var_newdim + m.var_append_dim - 1));

        if (row_append_mat)
          row_append_mat->concat_right(Sparsemat(row_append_mat->rowdim(), m.var_append_dim));

        var_append_dim += m.var_append_dim;
        var_newdim += m.var_append_dim;

      } // endif (m.var_append_dim>0)

      if (m.row_append_dim > 0) {

        if (row_append_mat) {
          if (rowapp) {
            row_append_mat->concat_below(*rowapp);
            delete rowapp;
            rowapp = 0;
          } else if (m.row_append_mat)
            row_append_mat->concat_below(*(m.row_append_mat));
          else
            row_append_mat->concat_below(Sparsemat(m.row_append_dim, row_append_mat->coldim()));
        } else {
          if (rowapp) {
            row_append_mat = new Sparsemat(row_append_dim, rowapp->coldim());
            row_append_mat->concat_below(*rowapp);
            delete rowapp;
            rowapp = 0;
          } else if (m.row_append_mat) {
            row_append_mat = new Sparsemat(row_append_dim, m.row_append_mat->coldim());
            row_append_mat->concat_below(*(m.row_append_mat));
          }
        }
        assert(rowapp == 0);

        if (row_append_rhslb) {
          if (m.row_append_rhslb)
            row_append_rhslb->concat_below(*(m.row_append_rhslb));
          else
            row_append_rhslb->enlarge_below(m.row_append_dim, m.rhs_minus_infinity);
        } else {
          if (m.row_append_rhslb) {
            row_append_rhslb = new Matrix(row_append_dim, 1, rhs_minus_infinity);
            row_append_rhslb->concat_below(*(m.row_append_rhslb));
          } else {
            if (rhs_minus_infinity != m.rhs_minus_infinity) {
              row_append_rhslb = new Matrix(row_append_dim, 1, rhs_minus_infinity);
              row_append_rhslb->enlarge_below(m.row_append_dim, m.rhs_minus_infinity);
            }
          }
        }

        if (row_append_rhsub) {
          if (m.row_append_rhsub)
            row_append_rhsub->concat_below(*(m.row_append_rhsub));
          else
            row_append_rhsub->enlarge_below(m.row_append_dim, m.rhs_plus_infinity);
        } else {
          if (m.row_append_rhsub) {
            row_append_rhsub = new Matrix(row_append_dim, 1, rhs_plus_infinity);
            row_append_rhsub->concat_below(*(m.row_append_rhsub));
          } else {
            if (rhs_plus_infinity != m.rhs_plus_infinity) {
              row_append_rhsub = new Matrix(row_append_dim, 1, rhs_plus_infinity);
              row_append_rhsub->enlarge_below(m.row_append_dim, m.rhs_plus_infinity);
            }
          }
        }

        if (row_map_to_old)
          row_map_to_old->concat_below(Range(row_olddim + row_append_dim, row_olddim + row_append_dim + m.row_append_dim - 1));

        if (row_new_ind)
          row_new_ind->concat_below(Range(row_newdim, row_newdim + m.row_append_dim - 1));
        else
          row_new_ind = new Indexmatrix(Range(row_newdim, row_newdim + m.row_append_dim - 1));

        row_append_dim += m.row_append_dim;
        row_newdim += m.row_append_dim;
      } // endif (m.row_append_dim>0)

      assert((var_append_lb == 0) || (var_append_lb->rowdim() == var_append_dim));
      assert((var_append_ub == 0) || (var_append_ub->rowdim() == var_append_dim));
      assert((var_append_costs == 0) || (var_append_costs->rowdim() == var_append_dim));
      assert((var_start_val == 0) || (var_start_val->rowdim() == var_append_dim));
      assert((var_append_cols == 0) || ((var_append_cols->rowdim() == row_olddim) && (var_append_cols->coldim() == var_append_dim)));
      assert((var_map_to_old == 0) || (var_map_to_old->coldim() == 1));
      assert((var_map_to_old == 0) || (var_map_to_old->rowdim() + var_del_ind->rowdim() == var_olddim + var_append_dim));
      assert((row_append_rhslb == 0) || (row_append_rhslb->rowdim() == row_append_dim));
      assert((row_append_rhsub == 0) || (row_append_rhsub->rowdim() == row_append_dim));
      assert((row_append_mat == 0) || ((row_append_mat->coldim() == var_newdim) && (row_append_mat->rowdim() == row_append_dim)));
      assert((row_map_to_old == 0) || (row_map_to_old->coldim() == 1));
      assert((row_map_to_old == 0) || (row_map_to_old->rowdim() + row_del_ind->rowdim() == row_olddim + row_append_dim));

      return err;
    }


    if (m.var_olddim != var_newdim) {
      if (cb_out())
        get_out() << "**** ERROR: Modification::incorporate(const Modfication& m): old variable dimension=" << m.var_olddim << " of m must match new variable dimension=" << var_newdim << " of this" << std::endl;
      err++;
    }
    if (m.row_olddim != row_newdim) {
      if (cb_out())
        get_out() << "**** ERROR: Modification::incorporate(const Modfication& m): old row dimsion=" << m.row_olddim << " of m must match new row dimension=" << row_newdim << " of this" << std::endl;
      err++;
    }
    if (err)
      return err;

    if (m.var_set_lb) {
      for (Realmap::const_iterator it = m.var_set_lb->begin(); it != m.var_set_lb->end(); it++) {
        if (add_set_lb(it->first, it->second)) {
          if (cb_out())
            get_out() << "**** ERROR: Modification::incorporate(const Modfication& m): incorporating the new lower bound value=" << it->second << " for index=" << it->first << " failed " << std::endl;
          err++;
        }
        assert(err == 0);
      }
    }

    if (m.var_set_ub) {
      for (Realmap::const_iterator it = m.var_set_ub->begin(); it != m.var_set_ub->end(); it++) {
        if (add_set_ub(it->first, it->second)) {
          if (cb_out())
            get_out() << "**** ERROR: Modification::incorporate(const Modfication& m): incorporating the new upper bound value=" << it->second << " for index=" << it->first << " failed " << std::endl;
          err++;
        }
        assert(err == 0);
      }
    }

    if (m.row_set_rhslb) {
      for (Realmap::const_iterator it = m.row_set_rhslb->begin(); it != m.row_set_rhslb->end(); it++) {
        if (add_set_rhslb(it->first, it->second)) {
          if (cb_out())
            get_out() << "**** ERROR: Modification::incorporate(const Modfication& m): incorporating the new right hand side lower bound value=" << it->second << " for index=" << it->first << " failed " << std::endl;
          err++;
        }
        assert(err == 0);
      }
    }

    if (m.row_set_rhsub) {
      for (Realmap::const_iterator it = m.row_set_rhsub->begin(); it != m.row_set_rhsub->end(); it++) {
        if (add_set_rhsub(it->first, it->second)) {
          if (cb_out())
            get_out() << "**** ERROR: Modification::incorporate(const Modfication& m): incorporating the new right hand side upper bound value=" << it->second << " for index=" << it->first << " failed " << std::endl;
          err++;
        }
        assert(err == 0);
      }
    }


    assert((var_append_lb == 0) || (var_append_lb->rowdim() == var_append_dim));
    assert((var_append_ub == 0) || (var_append_ub->rowdim() == var_append_dim));
    assert((var_append_costs == 0) || (var_append_costs->rowdim() == var_append_dim));
    assert((var_start_val == 0) || (var_start_val->rowdim() == var_append_dim));
    assert((var_append_cols == 0) || ((var_append_cols->rowdim() == row_olddim) && (var_append_cols->coldim() == var_append_dim)));
    assert((var_map_to_old == 0) || (var_map_to_old->coldim() == 1));
    assert((var_map_to_old == 0) || (var_map_to_old->rowdim() + var_del_ind->rowdim() == var_olddim + var_append_dim));
    assert((row_append_rhslb == 0) || (row_append_rhslb->rowdim() == row_append_dim));
    assert((row_append_rhsub == 0) || (row_append_rhsub->rowdim() == row_append_dim));
    assert((row_append_mat == 0) || ((row_append_mat->coldim() == var_newdim) && (row_append_mat->rowdim() == row_append_dim)));
    assert((row_map_to_old == 0) || (row_map_to_old->coldim() == 1));
    assert((row_map_to_old == 0) || (row_map_to_old->rowdim() + row_del_ind->rowdim() == row_olddim + row_append_dim));

    if (m.var_append_dim > 0) {
      Matrix* var_lb = m.var_append_lb;
      if ((var_lb == 0) && (bounds_minus_infinity != m.bounds_minus_infinity))
        var_lb = new Matrix(m.var_append_dim, 1, m.bounds_minus_infinity);
      Matrix* var_ub = m.var_append_ub;
      if ((var_ub == 0) && (bounds_plus_infinity != m.bounds_plus_infinity))
        var_ub = new Matrix(m.var_append_dim, 1, m.bounds_plus_infinity);
      Matrix* start = m.var_start_val;
      if ((start == 0) && (start_val_default != m.start_val_default))
        start = new Matrix(m.var_append_dim, 1, m.start_val_default);
      Matrix* costs = m.var_append_costs;
      if ((costs == 0) && (cost_default != m.cost_default))
        costs = new Matrix(m.var_append_dim, 1, cost_default);

      err += add_append_vars(m.var_append_dim, var_lb, var_ub, m.var_append_cols, start, costs);

      if (var_lb != m.var_append_lb)
        delete var_lb;
      if (var_ub != m.var_append_ub)
        delete var_ub;
      if (start != m.var_start_val)
        delete start;
      if (costs != m.var_append_costs)
        delete costs;

    }
    assert(err == 0);

    assert((var_append_lb == 0) || (var_append_lb->rowdim() == var_append_dim));
    assert((var_append_ub == 0) || (var_append_ub->rowdim() == var_append_dim));
    assert((var_append_costs == 0) || (var_append_costs->rowdim() == var_append_dim));
    assert((var_start_val == 0) || (var_start_val->rowdim() == var_append_dim));
    assert((var_append_cols == 0) || ((var_append_cols->rowdim() == row_olddim) && (var_append_cols->coldim() == var_append_dim)));
    assert((var_map_to_old == 0) || (var_map_to_old->coldim() == 1));
    assert((var_map_to_old == 0) || (var_map_to_old->rowdim() + var_del_ind->rowdim() == var_olddim + var_append_dim));
    assert((row_append_rhslb == 0) || (row_append_rhslb->rowdim() == row_append_dim));
    assert((row_append_rhsub == 0) || (row_append_rhsub->rowdim() == row_append_dim));
    assert((row_append_mat == 0) || ((row_append_mat->coldim() == var_newdim) && (row_append_mat->rowdim() == row_append_dim)));
    assert((row_map_to_old == 0) || (row_map_to_old->coldim() == 1));
    assert((row_map_to_old == 0) || (row_map_to_old->rowdim() + row_del_ind->rowdim() == row_olddim + row_append_dim));

    if (m.var_map_to_old)
      err += add_reassign_vars(*(m.var_map_to_old));
    assert(err == 0);

    assert((var_append_lb == 0) || (var_append_lb->rowdim() == var_append_dim));
    assert((var_append_ub == 0) || (var_append_ub->rowdim() == var_append_dim));
    assert((var_append_costs == 0) || (var_append_costs->rowdim() == var_append_dim));
    assert((var_start_val == 0) || (var_start_val->rowdim() == var_append_dim));
    assert((var_append_cols == 0) || ((var_append_cols->rowdim() == row_olddim) && (var_append_cols->coldim() == var_append_dim)));
    assert((var_map_to_old == 0) || (var_map_to_old->coldim() == 1));
    assert((var_map_to_old == 0) || (var_map_to_old->rowdim() + var_del_ind->rowdim() == var_olddim + var_append_dim));
    assert((row_append_rhslb == 0) || (row_append_rhslb->rowdim() == row_append_dim));
    assert((row_append_rhsub == 0) || (row_append_rhsub->rowdim() == row_append_dim));
    assert((row_append_mat == 0) || ((row_append_mat->coldim() == var_newdim) && (row_append_mat->rowdim() == row_append_dim)));
    assert((row_map_to_old == 0) || (row_map_to_old->coldim() == 1));
    assert((row_map_to_old == 0) || (row_map_to_old->rowdim() + row_del_ind->rowdim() == row_olddim + row_append_dim));

    if (m.row_append_dim > 0) {
      Matrix* row_rhslb = m.row_append_rhslb;
      if ((row_rhslb == 0) && (rhs_minus_infinity != m.rhs_minus_infinity))
        row_rhslb = new Matrix(m.row_append_dim, 1, m.rhs_minus_infinity);
      Matrix* row_rhsub = m.row_append_rhsub;
      if ((row_rhsub == 0) && (rhs_plus_infinity != m.rhs_plus_infinity))
        row_rhsub = new Matrix(m.row_append_dim, 1, m.rhs_plus_infinity);

      err += add_append_rows(m.row_append_dim, m.row_append_mat, row_rhslb, row_rhsub);

      if (row_rhslb != m.row_append_rhslb)
        delete row_rhslb;
      if (row_rhsub != m.row_append_rhsub)
        delete row_rhsub;

    }
    assert(err == 0);

    assert((var_append_lb == 0) || (var_append_lb->rowdim() == var_append_dim));
    assert((var_append_ub == 0) || (var_append_ub->rowdim() == var_append_dim));
    assert((var_append_costs == 0) || (var_append_costs->rowdim() == var_append_dim));
    assert((var_start_val == 0) || (var_start_val->rowdim() == var_append_dim));
    assert((var_append_cols == 0) || ((var_append_cols->rowdim() == row_olddim) && (var_append_cols->coldim() == var_append_dim)));
    assert((var_map_to_old == 0) || (var_map_to_old->coldim() == 1));
    assert((var_map_to_old == 0) || (var_map_to_old->rowdim() + var_del_ind->rowdim() == var_olddim + var_append_dim));
    assert((row_append_rhslb == 0) || (row_append_rhslb->rowdim() == row_append_dim));
    assert((row_append_rhsub == 0) || (row_append_rhsub->rowdim() == row_append_dim));
    assert((row_append_mat == 0) || ((row_append_mat->coldim() == var_newdim) && (row_append_mat->rowdim() == row_append_dim)));
    assert((row_map_to_old == 0) || (row_map_to_old->coldim() == 1));
    assert((row_map_to_old == 0) || (row_map_to_old->rowdim() + row_del_ind->rowdim() == row_olddim + row_append_dim));

    if (m.row_map_to_old)
      err += add_reassign_rows(*(m.row_map_to_old));
    assert(err == 0);

    assert((var_append_lb == 0) || (var_append_lb->rowdim() == var_append_dim));
    assert((var_append_ub == 0) || (var_append_ub->rowdim() == var_append_dim));
    assert((var_append_costs == 0) || (var_append_costs->rowdim() == var_append_dim));
    assert((var_start_val == 0) || (var_start_val->rowdim() == var_append_dim));
    assert((var_append_cols == 0) || ((var_append_cols->rowdim() == row_olddim) && (var_append_cols->coldim() == var_append_dim)));
    assert((var_map_to_old == 0) || (var_map_to_old->coldim() == 1));
    assert((var_map_to_old == 0) || (var_map_to_old->rowdim() + var_del_ind->rowdim() == var_olddim + var_append_dim));
    assert((row_append_rhslb == 0) || (row_append_rhslb->rowdim() == row_append_dim));
    assert((row_append_rhsub == 0) || (row_append_rhsub->rowdim() == row_append_dim));
    assert((row_append_mat == 0) || ((row_append_mat->coldim() == var_newdim) && (row_append_mat->rowdim() == row_append_dim)));
    assert((row_map_to_old == 0) || (row_map_to_old->coldim() == 1));
    assert((row_map_to_old == 0) || (row_map_to_old->rowdim() + row_del_ind->rowdim() == row_olddim + row_append_dim));



    return err;
  }

  // *****************************************************************************
  //                              Modification::apply_to_vars
  // *****************************************************************************

  int Modification::apply_to_vars(Matrix* vars,
    Matrix* lb,
    Matrix* ub,
    Matrix* costs) const {
    int err = 0;
    if ((vars) && ((vars->rowdim() != var_olddim) || (vars->coldim() != 1))) {
      err++;
      if (cb_out())
        get_out() << "**** ERROR: Modification::apply_to_vars(...): vector of variables has size " << vars->rowdim() << " x " << vars->coldim() << " but should have size " << var_olddim << " x 1" << std::endl;
    }
    if ((lb) && ((lb->rowdim() != var_olddim) || (lb->coldim() != 1))) {
      err++;
      if (cb_out())
        get_out() << "**** ERROR: Modification::apply_to_vars(...): vector of lower bounds has size " << lb->rowdim() << " x " << lb->coldim() << " but should have size " << var_olddim << " x 1" << std::endl;
    }
    if ((ub) && ((ub->rowdim() != var_olddim) || (ub->coldim() != 1))) {
      err++;
      if (cb_out())
        get_out() << "**** ERROR: Modification::apply_to_vars(...): vector of upper bounds has size " << ub->rowdim() << " x " << ub->coldim() << " but should have size " << var_olddim << " x 1" << std::endl;
    }
    if ((costs) && ((costs->rowdim() != var_olddim) || (costs->coldim() != 1))) {
      err++;
      if (cb_out())
        get_out() << "**** ERROR: Modification::apply_to_vars(...): vector of costs has size " << costs->rowdim() << " x " << costs->coldim() << " but should have size " << var_olddim << " x 1" << std::endl;
    }
    if (err)
      return err;

    //first change the bounds
    if ((var_set_lb) && (lb)) {
      for (Realmap::const_iterator it = var_set_lb->begin(); it != var_set_lb->end(); it++)
        (*lb)(it->first) = it->second;
    }
    if ((var_set_ub) && (ub)) {
      for (Realmap::const_iterator it = var_set_ub->begin(); it != var_set_ub->end(); it++)
        (*ub)(it->first) = it->second;
    }
    if ((enforce_bounds_consistency) && (lb) && (ub)) {
      if (var_set_lb) {
        for (Realmap::const_iterator it = var_set_lb->begin(); it != var_set_lb->end(); it++) {
          if ((*lb)(it->first) > (*ub)(it->first)) {
            if (cb_out())
              get_out() << "**** ERROR: Modification::apply_to_vars(...): new variable lower bound for mapped index " << it->first << " resulted in a lower bound=" << (*lb)(it->first) << " greater than the upper bound=" << (*ub)(it->first) << std::endl;
            err++;
          }
        }
      }
      if (var_set_ub) {
        for (Realmap::const_iterator it = var_set_ub->begin(); it != var_set_ub->end(); it++) {
          if ((*lb)(it->first) > (*ub)(it->first)) {
            if (cb_out())
              get_out() << "**** ERROR: Modification::apply_to_vars(...): new variable upper bound for mapped index " << it->first << " resulted in a lower bound=" << (*lb)(it->first) << " greater than the upper bound=" << (*ub)(it->first) << std::endl;
            err++;
          }
        }
      }
    }

    //next append 
    if (var_append_dim > 0) {
      if (costs) {
        if (var_append_costs) {
          costs->concat_below(*var_append_costs);
        } else {
          costs->concat_below(Matrix(var_append_dim, 1, cost_default));
        }
      }
      if ((enforce_bounds_consistency) && (lb) && (ub) && (var_append_lb) && (var_append_ub)) {
        for (Integer i = 0; i < var_append_dim; i++) {
          if ((*var_append_lb)(i) > (*var_append_ub)(i)) {
            if (cb_out())
              get_out() << "**** ERROR: Modification::apply_to_vars(...): the bounds appended variable index" << i << " resulted in a lower bound=" << (*var_append_lb)(i) << " greater than the upper bound=" << (*var_append_ub)(i) << std::endl;
            err++;
          }
        }
      }
      if (lb) {
        if (var_append_lb) {
          lb->concat_below(*var_append_lb);
        } else {
          lb->concat_below(Matrix(var_append_dim, 1, bounds_minus_infinity));
        }
      }
      if (ub) {
        if (var_append_ub) {
          ub->concat_below(*var_append_ub);
        } else {
          ub->concat_below(Matrix(var_append_dim, 1, bounds_plus_infinity));
        }
      }

      if (vars) {
        if (var_start_val) {
          vars->concat_below(*var_start_val);
        } else {
          vars->concat_below(Matrix(var_append_dim, 1, start_val_default));
        }
        if ((enforce_start_val_box_feasibility) && ((lb) || (ub))) {
          for (Integer i = 0; i < vars->dim(); i++) {
            if ((lb) && ((*vars)(i) < (*lb)(i))) {
              (*vars)(i) = (*lb)(i);
              continue;
            }
            if ((ub) && ((*vars)(i) > (*ub)(i))) {
              (*vars)(i) = (*ub)(i);
            }
          }
        }
      } //endif (vars)			   
    } //endif (var_append_dim>0)

    //now apply the map if needed
    if (var_map_to_old) {
      if (vars)
        (*vars) = vars->rows(*var_map_to_old);
      if (lb)
        (*lb) = lb->rows(*var_map_to_old);
      if (ub)
        (*ub) = ub->rows(*var_map_to_old);
      if (costs)
        (*costs) = costs->rows(*var_map_to_old);
    }

    return err;
  }

  // *****************************************************************************
  //                              Modification::apply_to_rows
  // *****************************************************************************

  int Modification::apply_to_rows(Sparsemat* rows,
    Matrix* rhslb, Matrix* rhsub) const {
    int err = 0;
    if ((rows) && ((rows->coldim() != var_olddim) || (rows->rowdim() != row_olddim))) {
      err++;
      if (cb_out())
        get_out() << "**** ERROR: Modification::apply_to_rows(...): matrix has size " << rows->rowdim() << " x " << rows->coldim() << " but should have size " << row_olddim << " x " << var_olddim << std::endl;
    }
    if ((rhslb) && ((rhslb->rowdim() != row_olddim) || (rhslb->coldim() != 1))) {
      err++;
      if (cb_out())
        get_out() << "**** ERROR: Modification::apply_to_rows(...): rhslb vector has size " << rhslb->rowdim() << " x " << rhslb->coldim() << " but should have size " << row_olddim << " x 1" << std::endl;
    }
    if ((rhsub) && ((rhsub->rowdim() != row_olddim) || (rhsub->coldim() != 1))) {
      err++;
      if (cb_out())
        get_out() << "**** ERROR: Modification::apply_to_rows(...): rhsub vector has size " << rhsub->rowdim() << " x " << rhsub->coldim() << " but should have size " << row_olddim << " x 1" << std::endl;
    }

    if (err)
      return err;

    //first change the bounds
    if ((row_set_rhslb) && (rhslb)) {
      for (Realmap::const_iterator it = row_set_rhslb->begin(); it != row_set_rhslb->end(); it++)
        (*rhslb)(it->first) = it->second;
    }
    if ((row_set_rhsub) && (rhsub)) {
      for (Realmap::const_iterator it = row_set_rhsub->begin(); it != row_set_rhsub->end(); it++)
        (*rhsub)(it->first) = it->second;
    }
    if ((enforce_bounds_consistency) && (rhslb) && (rhsub)) {
      if (row_set_rhslb) {
        for (Realmap::const_iterator it = row_set_rhslb->begin(); it != row_set_rhslb->end(); it++) {
          if ((*rhslb)(it->first) > (*rhsub)(it->first)) {
            if (cb_out())
              get_out() << "**** ERROR: Modification::apply_to_rows(...): new right hand side lower bound for mapped index " << it->first << " resulted in a lower bound=" << (*rhslb)(it->first) << " greater than the upper bound=" << (*rhsub)(it->first) << std::endl;
            err++;
          }
        }
      }
      if (row_set_rhsub) {
        for (Realmap::const_iterator it = row_set_rhsub->begin(); it != row_set_rhsub->end(); it++) {
          if ((*rhslb)(it->first) > (*rhsub)(it->first)) {
            if (cb_out())
              get_out() << "**** ERROR: Modification::apply_to_rows(...): new right hand side upper bound for mapped index " << it->first << " resulted in a lower bound=" << (*rhslb)(it->first) << " greater than the upper bound=" << (*rhsub)(it->first) << std::endl;
            err++;
          }
        }
      }
    }

    //next append the columns
    if ((var_append_dim > 0) && (rows)) {
      if (var_append_cols)
        rows->concat_right(*var_append_cols);
      else
        rows->concat_right(Sparsemat(row_olddim, var_append_dim));
    }

    //now reorder the colums if necessary
    if ((var_map_to_old) && (rows)) {
      *rows = rows->cols(*var_map_to_old);
    }

    //next append the new rows
    if (row_append_dim > 0) {
      if (rows) {
        if (row_append_mat) {
          rows->concat_below(*row_append_mat);
        } else {
          rows->concat_below(Sparsemat(row_append_dim, var_newdim));
        }
      }
      if (rhslb) {
        if (row_append_rhslb) {
          rhslb->concat_below(*row_append_rhslb);
        } else {
          rhslb->concat_below(Matrix(row_append_dim, 1, rhs_minus_infinity));
        }
      }
      if (rhsub) {
        if (row_append_rhsub) {
          rhsub->concat_below(*row_append_rhsub);
        } else {
          rhsub->concat_below(Matrix(row_append_dim, 1, rhs_plus_infinity));
        }
      }
      if ((enforce_bounds_consistency) && (rhslb) && (rhsub) && (row_append_rhslb) && (row_append_rhsub)) {
        for (Integer i = 0; i < row_append_dim; i++) {
          if ((*row_append_rhslb)(i) > (*row_append_rhsub)(i)) {
            if (cb_out())
              get_out() << "**** ERROR: Modification::apply_to_rows(...): the bounds appended with row index" << i << " resulted in a lower bound=" << (*row_append_rhslb)(i) << " greater than the upper bound=" << (*row_append_rhsub)(i) << std::endl;
            err++;
          }
        }
      }
    }

    //now reorder the row if necessary
    if (row_map_to_old) {
      if (rows)
        *rows = rows->rows(*row_map_to_old);
      if (rhslb)
        *rhslb = rhslb->rows(*row_map_to_old);
      if (rhsub)
        *rhsub = rhsub->rows(*row_map_to_old);
    }

    return err;
  }

  // *****************************************************************************
  //                              Modification::no_modification
  // *****************************************************************************

  bool Modification::no_modification() const {
    if ((var_set_lb == 0) &&
      (var_set_ub == 0) &&
      (var_append_dim == 0) &&
      (var_map_to_old == 0) &&
      (row_set_rhslb == 0) &&
      (row_set_rhsub == 0) &&
      (row_append_dim == 0) &&
      (row_map_to_old == 0)
      )
      return true;
    return false;
  }

  // *****************************************************************************
  //                              Modification::set_append_to_old
  // *****************************************************************************

  int Modification::set_append_to_old(bool ao) {
    if (ao) {
      if (!(
        (var_set_lb == 0) &&
        (var_set_ub == 0) &&
        (var_map_to_old == 0) &&
        (row_set_rhslb == 0) &&
        (row_set_rhsub == 0) &&
        (row_map_to_old == 0)
        )) {
        if (cb_out()) {
          get_out() << "**** ERROR Modification::set_append_to_old(.): failed to set append_to_old to true because some operations stored here are not only of append type" << std::endl;
        }
        return 1;
      }
    }
    append_only = ao;
    return 0;
  }

  // *****************************************************************************
  //                    Modification::deleted_variables_are_zero
  // *****************************************************************************

  bool Modification::deleted_variables_are_zero(const Matrix& oldpoint) const {
    if ((oldpoint.rowdim() != var_olddim) || (oldpoint.coldim() != 1)) {
      if (cb_out()) {
        get_out() << "**** ERROR Modification::deleted_variables_are_zero(.): point has dim=" << oldpoint.rowdim() << " x " << oldpoint.coldim() << " but should have dim=" << var_olddim << " x 1" << std::endl;
      }
      std::abort();
    }
    if (var_del_ind) {
      Integer n = var_del_ind->rowdim();
      for (Integer i = 0; i < n; i++) {
        Integer ind = (*var_del_ind)(i);
        if (oldpoint(ind) != 0.) {
          return false;
        }
      }
    }
    return true;
  }

  // *****************************************************************************
  //                    Modification::new_variables_are_zero
  // *****************************************************************************

  bool Modification::new_variables_are_zero(const Matrix& newpoint) const {
    if ((newpoint.rowdim() != var_newdim) || (newpoint.coldim() != 1)) {
      if (cb_out()) {
        get_out() << "**** ERROR Modification::new_variables_are_zero(.): point has dim=" << newpoint.rowdim() << " x " << newpoint.coldim() << " but should have dim=" << var_newdim << " x 1" << std::endl;
      }
      std::abort();
    }
    if (var_new_ind) {
      Integer n = var_new_ind->rowdim();
      for (Integer i = 0; i < n; i++) {
        Integer ind = (*var_new_ind)(i);
        if (newpoint(ind) != 0.) {
          return false;
        }
      }
    }
    return true;
  }

  // *****************************************************************************
  //                    Modification::mapped_variables_are_equal
  // *****************************************************************************

  bool Modification::mapped_variables_are_equal(const Matrix& newpoint, const Matrix& oldpoint) const {
    int err = 0;
    if ((oldpoint.rowdim() != var_olddim) || (oldpoint.coldim() != 1)) {
      if (cb_out()) {
        get_out() << "**** ERROR Modification::deleted_variables_are_equal(.): point has dim=" << oldpoint.rowdim() << " x " << oldpoint.coldim() << " but should have dim=" << var_olddim << " x 1" << std::endl;
      }
      err++;
    }
    if ((newpoint.rowdim() != var_newdim) || (newpoint.coldim() != 1)) {
      if (cb_out()) {
        get_out() << "**** ERROR Modification::new_variables_are_equal(.): point has dim=" << oldpoint.rowdim() << " x " << newpoint.coldim() << " but should have dim=" << var_newdim << " x 1" << std::endl;
      }
      err++;
    }
    if (err)
      std::abort();

    if (var_map_to_old) {
      Integer n = var_map_to_old->rowdim();
      for (Integer i = 0; i < n; i++) {
        Integer ind = (*var_map_to_old)(i);
        if ((ind < var_olddim) && (std::fabs(newpoint(i) - oldpoint(ind)) > 1e2 * eps_Real * (1. + std::fabs(newpoint(i))))) {
          return false;
        }
      }
    } else {
      for (Integer i = 0; i < var_olddim; i++) {
        if (std::fabs(newpoint(i) - oldpoint(i)) > 1e2 * eps_Real * (1. + std::fabs(newpoint(i)))) {
          return false;
        }
      }
    }
    return true;
  }


}

