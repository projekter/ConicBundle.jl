/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/PSCAffineModification.cxx
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
#include "PSCAffineFunction.hxx"
#include "GroundsetModification.hxx"

using namespace CH_Matrix_Classes;

namespace ConicBundle {

  // *****************************************************************************
  //                              PSCAffineModification::PSCAffineModification
  // *****************************************************************************

  PSCAffineModification::PSCAffineModification(Integer volddim,
    const Indexmatrix& bolddim,
    const CBout* cb,
    int incr) :
    ModificationBase(cb, incr) {
    append_only = false;

    var_newdim = var_olddim = 0;
    var_append_dim = 0;
    block_olddim.init(0, 1, Integer(0));
    block_newdim = block_olddim;
    var_append.init(block_olddim, var_append_dim);
    var_del_ind = 0;
    var_map_to_old = 0;
    var_new_ind = 0;

    block_append_dim.init(0, 1, Integer(0));
    block_append.init(block_append_dim, var_newdim);
    offset_append.init(block_append_dim, 1);
    block_del_ind = 0;
    block_map_to_old = 0;
    block_new_ind = 0;

    reset_primal = false;
    generating_primal = 0;
    skip_extension = false;

    clear(volddim, bolddim);
  }

  // *****************************************************************************
  //                              PSCAffineModification::~PSCAffineModification
  // *****************************************************************************

  PSCAffineModification::~PSCAffineModification() {
    clear(0, Indexmatrix(0, 1, Integer(0)));
  }

  // *****************************************************************************
  //                              PSCAffineModification::PSCAffineModification
  // *****************************************************************************

  int PSCAffineModification::clear(Integer volddim, const Indexmatrix& bolddim) {
    int err = 0;
    if (bolddim.dim() > 0) {
      if (min(bolddim) < 0) {
        if (cb_out()) {
          get_out() << "**** ERROR: PSCAffineModification::clear(..): input block dimension vector has negative entries, block_olddim=" << block_olddim;
        }
        err++;
      }
      if (bolddim.coldim() != 1) {
        if (cb_out()) {
          get_out() << "**** WARNING: PSCAffineModification::clear(..): input block dimension is treated as a column vector but has coldim=" << bolddim.coldim() << std::endl;
        }
      }
    }
    if (volddim < 0) {
      if (cb_out()) {
        get_out() << "**** WARNING: PSCAffineModification::clear(..): input variable dimension =" << volddim << " < 0, it is treated as 0" << std::endl;
      }
    }
    if (err) {
      if (cb_out()) {
        get_out() << "**** ERROR: PSCAffineModification::clear(..): fatal errors occured, so nothing is cleared; return value =" << err << std::endl;
      }
      return err;
    }

    append_only = false;

    var_olddim = max(0, volddim);
    block_olddim = bolddim;

    var_newdim = var_olddim;
    var_append_dim = 0;
    block_newdim = block_olddim;
    var_append.init(block_olddim, var_append_dim);
    delete var_del_ind;
    var_del_ind = 0;
    delete var_map_to_old;
    var_map_to_old = 0;
    delete var_new_ind;
    var_new_ind = 0;

    block_append_dim.init(0, 0, Integer(0));
    block_append.init(block_append_dim, var_newdim);
    offset_append.init(block_append_dim, 1);
    delete block_del_ind;
    block_del_ind = 0;
    delete block_map_to_old;
    block_map_to_old = 0;
    delete block_new_ind;
    block_new_ind = 0;

    reset_primal = false;
    delete generating_primal;
    generating_primal = 0;
    skip_extension = false;

    return err;
  }



  // *****************************************************************************
  //                           PSCAffineModification::add_append_vars
  // *****************************************************************************

  int PSCAffineModification::add_append_vars(Integer append_dim,
    const SparseCoeffmatMatrix* append_cols) {
    int err = 0;
    if (append_dim < 0) {
      err++;
      if (cb_out())
        get_out() << "**** ERROR: PSCAffineModification::add_append_variables(..): append_dim=" << append_dim << ", cannot append negative number of variables" << std::endl;
    }
    if (append_cols != 0) {
      if (((var_newdim > 0) || (block_newdim.dim() > 0)) && (block_newdim.dim() != append_cols->rowdim())) {
        err++;
        if (cb_out())
          get_out() << "**** ERROR: PSCAffineModification::add_append_variables(..): current number of blocks =" << block_newdim.dim() << " != " << append_cols->rowdim() << "= the number of blocks in the appended columns" << std::endl;
      } else {
        for (Integer i = 0; i < block_newdim.dim(); i++) {
          if (block_newdim(i) != append_cols->blockdim(i)) {
            err++;
            if (cb_out())
              get_out() << "**** ERROR: PSCAffineModification::add_append_variables(..): order of block(" << i << ") =" << block_newdim(i) << " != " << append_cols->blockdim(i) << "= order of block(" << i << ") in the appended columns" << std::endl;
          }
        }
      }
    }

    if (err) {
      if (cb_out())
        get_out() << "**** ERROR: PSCAffineModification::add_append_variables(..): fatal input errors occured, leaving without appending anything; return value=" << err << std::endl;
      return err;
    }

    if ((block_newdim.dim() == 0) && (append_cols) && (append_cols->rowdim() > 0)) {
      assert(var_newdim == 0);
      if (add_append_blocks(append_cols->blockdim(), 0, 0)) {
        if (cb_out())
          get_out() << "**** ERROR: PSCAffineModification::add_append_vars(...): add_append_blocks(.) failed" << std::endl;
        return err;
      }
    }

    if (append_dim == 0) {
      return err;
    }

    var_append_dim += append_dim;
    var_newdim += append_dim;

    if (var_new_ind == 0)
      var_new_ind = new Indexmatrix(0, 1, Integer(0));
    var_new_ind->concat_below(Range(var_newdim - append_dim, var_newdim - 1));

    if (var_map_to_old)
      var_map_to_old->concat_below(Range(var_olddim + var_append_dim - append_dim, var_olddim + var_append_dim - 1));

    if (append_cols) {

      if (block_map_to_old == 0) {
        if (block_append_dim.dim() == 0) {
          var_append.append_columns(*append_cols);
          block_append.append_columns(SparseCoeffmatMatrix(Indexmatrix(0, 0, Integer(0)), append_cols->coldim()));
        } else {
          Indexmatrix ind(Range(0, block_olddim.dim() - 1));
          var_append.append_columns(*append_cols, &ind);
          ind.init(Range(block_olddim.dim(), block_newdim.dim() - 1));
          block_append.append_columns(*append_cols, &ind);
        }
      } else { //block_map_to_old!=0; //need to transform indices back
        Indexmatrix old_blocks = -block_olddim;
        Indexmatrix new_blocks(block_append_dim.dim(), 1, Integer(0));
        for (Integer i = 0; i < block_map_to_old->dim(); i++) {
          Integer j = (*block_map_to_old)(i);
          if (j >= old_blocks.dim()) {
            new_blocks(j - old_blocks.dim()) = i;
          } else {
            old_blocks(j) = i;
          }
        }
        var_append.append_columns(*append_cols, &old_blocks);
        block_append.append_columns(*append_cols, &new_blocks);
      } //end of else //row_map_to_old!=0
    } //end of if append_cols != 0	     
    else {
      var_append.append_columns(SparseCoeffmatMatrix(block_olddim, append_dim));
      block_append.append_columns(SparseCoeffmatMatrix(block_append_dim, append_dim));
    }


    assert((var_append.rowdim() == block_olddim.dim()) && (var_append.coldim() == var_append_dim));
    assert((block_append.rowdim() == block_append_dim.dim()) && (block_append.coldim() == var_newdim));
    assert((var_map_to_old == 0) || (var_map_to_old->coldim() == 1));
    assert((var_map_to_old == 0) || (var_map_to_old->rowdim() + var_del_ind->rowdim() == var_olddim + var_append_dim));
    assert((var_map_to_old == 0) || (var_map_to_old->rowdim() == var_newdim));
    assert((block_map_to_old == 0) || (block_map_to_old->coldim() == 1));
    assert((block_map_to_old == 0) || (block_map_to_old->rowdim() + block_del_ind->rowdim() == block_olddim.dim() + block_append_dim.dim()));
    assert((block_map_to_old == 0) || (block_map_to_old->rowdim() == block_newdim.rowdim()));

    return 0;
  }

  // *****************************************************************************
  //                              PSCAffineModification::add_reassign_variables
  // *****************************************************************************

  int PSCAffineModification::add_reassign_vars(const Indexmatrix& map_to_old) {
    int err = 0;
    if (append_only) {
      if (cb_out())
        get_out() << "**** ERROR: PSCAffineModification::add_reassign_vars(.): append_to_old is set to true, so this operations is not allowed" << std::endl;
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

      var_append.init(block_olddim, 0);

      //block_newdim does not change
      //block_append_dim does not change
      block_append.init(block_append_dim, 0);
      //offset_append does not change
      //block_del_ind does not change
      //block_map_to_old does not change;

      assert((var_append.rowdim() == block_olddim.dim()) && (var_append.coldim() == var_append_dim));
      assert((block_append.rowdim() == block_append_dim.dim()) && (block_append.coldim() == var_newdim));
      assert((var_map_to_old == 0) || (var_map_to_old->coldim() == 1));
      assert((var_map_to_old == 0) || (var_map_to_old->rowdim() + var_del_ind->rowdim() == var_olddim + var_append_dim));
      assert((var_map_to_old == 0) || (var_map_to_old->rowdim() == var_newdim));
      assert((block_map_to_old == 0) || (block_map_to_old->coldim() == 1));
      assert((block_map_to_old == 0) || (block_map_to_old->rowdim() + block_del_ind->rowdim() == block_olddim.dim() + block_append_dim.dim()));
      assert((block_map_to_old == 0) || (block_map_to_old->rowdim() == block_newdim.rowdim()));

      return 0;
    }

    block_append.reassign_columns(map_to_old);

    Indexmatrix append_del_ind;

    err = adapt_map_to_old(var_map_to_old, var_del_ind, var_new_ind, append_del_ind, map_to_old, var_append_dim, var_olddim, var_newdim);
    if (err) {
      if (cb_out())
        get_out() << "**** ERROR: PSCAffineModification::add_reassign_variables(...): adapt_map_to_old(...) failed and returned " << err << std::endl;
      return err;
    }

    var_newdim = var_map_to_old->dim();
    var_append.delete_columns(append_del_ind);

    assert((var_append.rowdim() == block_olddim.dim()) && (var_append.coldim() == var_append_dim));
    assert((block_append.rowdim() == block_append_dim.dim()) && (block_append.coldim() == var_newdim));
    assert((var_map_to_old == 0) || (var_map_to_old->coldim() == 1));
    assert((var_map_to_old == 0) || (var_map_to_old->rowdim() + var_del_ind->rowdim() == var_olddim + var_append_dim));
    assert((var_map_to_old == 0) || (var_map_to_old->rowdim() == var_newdim));
    assert((block_map_to_old == 0) || (block_map_to_old->coldim() == 1));
    assert((block_map_to_old == 0) || (block_map_to_old->rowdim() + block_del_ind->rowdim() == block_olddim.dim() + block_append_dim.dim()));
    assert((block_map_to_old == 0) || (block_map_to_old->rowdim() == block_newdim.rowdim()));

    return 0;
  }

  // *****************************************************************************
  //                              PSCAffineModification::add_delete_variables
  // *****************************************************************************

  int PSCAffineModification::add_delete_vars(const Indexmatrix& delete_variables,
    Indexmatrix& map_to_old) {
    int err = 0;
    if (append_only) {
      if (cb_out())
        get_out() << "**** ERROR: PSCAffineModification::add_delete_vars(..): append_to_old is set to true, so this operations is not allowed" << std::endl;
      err++;
      return err;
    }
    err = form_map_to_old(map_to_old, delete_variables, var_newdim);
    if (err) {
      if (cb_out())
        get_out() << "**** ERROR: PSCAffineModification::add_delete_vars(...): form_map_to_old(...) failed and returned " << err << std::endl;
      return err;
    }

    err = add_reassign_vars(map_to_old);
    if (err) {
      if (cb_out())
        get_out() << "**** ERROR: PSCAffineModification::add_delete_vars(...): add_reassign_vars(...) failed and returned " << err << std::endl;
    }

    assert((var_append.rowdim() == block_olddim.dim()) && (var_append.coldim() == var_append_dim));
    assert((block_append.rowdim() == block_append_dim.dim()) && (block_append.coldim() == var_newdim));
    assert((var_map_to_old == 0) || (var_map_to_old->coldim() == 1));
    assert((var_map_to_old == 0) || (var_map_to_old->rowdim() + var_del_ind->rowdim() == var_olddim + var_append_dim));
    assert((var_map_to_old == 0) || (var_map_to_old->rowdim() == var_newdim));
    assert((block_map_to_old == 0) || (block_map_to_old->coldim() == 1));
    assert((block_map_to_old == 0) || (block_map_to_old->rowdim() + block_del_ind->rowdim() == block_olddim.dim() + block_append_dim.dim()));
    assert((block_map_to_old == 0) || (block_map_to_old->rowdim() == block_newdim.rowdim()));

    return err;
  }



  // *****************************************************************************
  //                        PSCAffineModification::add_append_blocks
  // *****************************************************************************

  int PSCAffineModification::add_append_blocks(const Indexmatrix& append_dim,
    const SparseCoeffmatMatrix* append_offsets,
    const SparseCoeffmatMatrix* append_blocks) {
    int err = 0;
    if (append_dim.dim() > 0) {
      if (append_dim.coldim() != 1) {
        err++;
        if (cb_out())
          get_out() << "**** ERROR: PSCAffineModification::add_append_blocks(...): append_dim.coldim()=" << append_dim.coldim() << " but should be 1 (a column vector)" << std::endl;
      }
      if (min(append_dim) < 0) {
        err++;
        if (cb_out())
          get_out() << "**** ERROR: PSCAffineModification::add_append_blocks(...): min(append_dim)=" << min(append_dim) << " but no blocks of negative order are allowed" << std::endl;
      }
    }
    if (append_offsets) {
      if ((append_offsets->rowdim() != append_dim.dim()) ||
        (append_offsets->coldim() != 1)) {
        err++;
        if (cb_out())
          get_out() << "**** ERROR: PSCAffineModification::add_append_blocks(...): matrix of offsets to be appended has size " << append_offsets->rowdim() << " x " << append_offsets->coldim() << " but should have size " << append_dim.dim() << " x 1" << std::endl;
      } else {
        for (Integer i = 0; i < append_dim.dim(); i++) {
          if (append_offsets->blockdim(i) != append_dim(i)) {
            err++;
            if (cb_out())
              get_out() << "**** ERROR: PSCAffineModification::add_append_blocks(...): block " << i << " of matrix to be appended has size " << append_offsets->blockdim(i) << " != " << append_dim(i) << ", which is the the input block size given by append_dim(" << i << ")" << std::endl;
          }
        }
      }
    }
    if (append_blocks) {
      if ((append_blocks->rowdim() != append_dim.dim()) ||
        (((var_newdim > 0) || (block_newdim.dim() > 0)) && (append_blocks->coldim() != var_newdim))) {
        err++;
        if (cb_out())
          get_out() << "**** ERROR: PSCAffineModification::add_append_blocks(...): matrix to be appended has size " << append_blocks->rowdim() << " x " << append_blocks->coldim() << " but should have size " << append_dim.dim() << " x " << var_newdim << std::endl;
      } else {
        for (Integer i = 0; i < append_dim.dim(); i++) {
          if (append_blocks->blockdim(i) != append_dim(i)) {
            err++;
            if (cb_out())
              get_out() << "**** ERROR: PSCAffineModification::add_append_blocks(...): block " << i << " of matrix to be appended has size " << append_blocks->blockdim(i) << " != " << append_dim(i) << ", which is the the input block size given by append_dim(" << i << ")" << std::endl;
          }
        }
      }
    }


    if (err) {
      if (cb_out())
        get_out() << "**** ERROR: PSCAffineModification::add_append_blocks(...): fatal input errors occured, leaving without appending anything; return value=" << err << std::endl;
      return err;
    }

    if ((var_newdim == 0) && (append_blocks) && (append_blocks->coldim() > 0)) {
      assert(block_newdim.dim() == 0);
      if (add_append_vars(append_blocks->coldim(), 0)) {
        if (cb_out())
          get_out() << "**** ERROR: PSCAffineModification::add_append_blocks(...): add_append_vars(.) failed" << std::endl;
        return err;
      }
    }

    if (append_dim.dim() == 0) {
      return err;
    }

    block_newdim.concat_below(append_dim);
    block_append_dim.concat_below(append_dim);

    if (block_new_ind == 0)
      block_new_ind = new Indexmatrix;
    block_new_ind->concat_below(Range(block_newdim.dim() - append_dim.dim(), block_newdim.dim() - 1));

    if (append_blocks) {
      if (block_append.append_blocks(*append_blocks)) {
        if (cb_out())
          get_out() << "**** ERROR: PSCAffineModification::add_append_blocks(...); block_append.append_blocks failed (error should have been detected on input)" << std::endl;
        err++;
      }
    } else {
      if (block_append.append_blocks(SparseCoeffmatMatrix(append_dim, var_newdim))) {
        if (cb_out())
          get_out() << "**** ERROR: PSCAffineModification::add_append_blocks(...); empty block_append.append_blocks failed (error should have been detected on input)" << std::endl;
        err++;
      }
    }


    if (append_offsets) {
      if (offset_append.append_blocks(*append_offsets)) {
        if (cb_out())
          get_out() << "**** ERROR: PSCAffineModification::add_append_blocks(...); offset_append.append_blocks failed (error should have been detected on input)" << std::endl;
        err++;
      }
    } else {
      if (offset_append.append_blocks(SparseCoeffmatMatrix(append_dim, 1))) {
        if (cb_out())
          get_out() << "**** ERROR: PSCAffineModification::add_append_blocks(...); empty offset_append.append_blocks failed (error should have been detected on input)" << std::endl;
        err++;
      }
    }

    if (block_map_to_old)
      block_map_to_old->concat_below(Range(block_olddim.dim() + block_append_dim.dim() - append_dim.dim(), block_olddim.dim() + block_append_dim.dim() - 1));


    assert((var_append.rowdim() == block_olddim.dim()) && (var_append.coldim() == var_append_dim));
    assert((block_append.rowdim() == block_append_dim.dim()) && (block_append.coldim() == var_newdim));
    assert((var_map_to_old == 0) || (var_map_to_old->coldim() == 1));
    assert((var_map_to_old == 0) || (var_map_to_old->rowdim() + var_del_ind->rowdim() == var_olddim + var_append_dim));
    assert((var_map_to_old == 0) || (var_map_to_old->rowdim() == var_newdim));
    assert((block_map_to_old == 0) || (block_map_to_old->coldim() == 1));
    assert((block_map_to_old == 0) || (block_map_to_old->rowdim() + block_del_ind->rowdim() == block_olddim.dim() + block_append_dim.dim()));
    assert((block_map_to_old == 0) || (block_map_to_old->rowdim() == block_newdim.rowdim()));

    return 0;
  }


  // *****************************************************************************
  //                              PSCAffineModification::add_reassign_blocks
  // *****************************************************************************

  int PSCAffineModification::add_reassign_blocks(const Indexmatrix& map_to_old) {
    int err = 0;
    if (append_only) {
      if (cb_out())
        get_out() << "**** ERROR: PSCAffineModification::add_reassign_blocks(.): append_to_old is set to true, so this operations is not allowed" << std::endl;
      err++;
      return err;
    }
    if (map_to_old.dim() == 0) {
      block_newdim.init(0, 1, Integer(0));
      block_append_dim.init(0, 1, Integer(0));

      if (block_map_to_old == 0)
        block_map_to_old = new Indexmatrix(0, 1, Integer(0));
      else
        block_map_to_old->init(0, 1, Integer(0));
      if (block_del_ind == 0)
        block_del_ind = new Indexmatrix(Range(0, block_olddim.dim() - 1));
      else
        block_del_ind->init(Range(0, block_olddim.dim() - 1));

      var_append.init(block_olddim, var_append.coldim());
      block_append.init(Indexmatrix(0, 1, Integer(0)), var_newdim);
      offset_append.init(Indexmatrix(0, 1, Integer(0)), 1);

      assert((var_append.rowdim() == block_olddim.dim()) && (var_append.coldim() == var_append_dim));
      assert((block_append.rowdim() == block_append_dim.dim()) && (block_append.coldim() == var_newdim));
      assert((var_map_to_old == 0) || (var_map_to_old->coldim() == 1));
      assert((var_map_to_old == 0) || (var_map_to_old->rowdim() + var_del_ind->rowdim() == var_olddim + var_append_dim));
      assert((var_map_to_old == 0) || (var_map_to_old->rowdim() == var_newdim));
      assert((block_map_to_old == 0) || (block_map_to_old->coldim() == 1));
      assert((block_map_to_old == 0) || (block_map_to_old->rowdim() + block_del_ind->rowdim() == block_olddim.dim() + block_append_dim.dim()));
      assert((block_map_to_old == 0) || (block_map_to_old->rowdim() == block_newdim.rowdim()));


      return 0;
    }

    Indexmatrix append_del_ind;

    Integer bappdim = block_append_dim.dim();

    err = adapt_map_to_old(block_map_to_old, block_del_ind, block_new_ind, append_del_ind, map_to_old, bappdim, block_olddim.dim(), block_newdim.dim());
    if (err) {
      if (cb_out())
        get_out() << "**** ERROR: PSCAffineModification::add_reassign_blocks(...): adapt_map_to_old(...) failed and returned " << err << std::endl;
      return err;
    }

    block_newdim = block_newdim(map_to_old);
    block_append_dim.delete_rows(append_del_ind);

    block_append.delete_blocks(append_del_ind);
    offset_append.delete_blocks(append_del_ind);

    assert((var_append.rowdim() == block_olddim.dim()) && (var_append.coldim() == var_append_dim));
    assert((block_append.rowdim() == block_append_dim.dim()) && (block_append.coldim() == var_newdim));
    assert((var_map_to_old == 0) || (var_map_to_old->coldim() == 1));
    assert((var_map_to_old == 0) || (var_map_to_old->rowdim() + var_del_ind->rowdim() == var_olddim + var_append_dim));
    assert((var_map_to_old == 0) || (var_map_to_old->rowdim() == var_newdim));
    assert((block_map_to_old == 0) || (block_map_to_old->coldim() == 1));
    assert((block_map_to_old == 0) || (block_map_to_old->rowdim() + block_del_ind->rowdim() == block_olddim.dim() + block_append_dim.dim()));
    assert((block_map_to_old == 0) || (block_map_to_old->rowdim() == block_newdim.rowdim()));

    return 0;
  }

  // *****************************************************************************
  //                              PSCAffineModification::add_delete_blocks
  // *****************************************************************************

  int PSCAffineModification::add_delete_blocks(const Indexmatrix& delete_blocks,
    Indexmatrix& map_to_old) {
    int err = form_map_to_old(map_to_old, delete_blocks, block_newdim.dim());
    if (err) {
      if (cb_out())
        get_out() << "**** ERROR: PSCAffineModification::add_delete_blocks(...): form_map_to_old(...) failed and returned " << err << std::endl;
      return err;
    }

    err = add_reassign_blocks(map_to_old);
    if (err) {
      if (cb_out())
        get_out() << "**** ERROR: PSCAffineModification::add_delete_blocks(...): add_reassign_blocks(...) failed and returned " << err << std::endl;
    }

    assert((var_append.rowdim() == block_olddim.dim()) && (var_append.coldim() == var_append_dim));
    assert((block_append.rowdim() == block_append_dim.dim()) && (block_append.coldim() == var_newdim));
    assert((var_map_to_old == 0) || (var_map_to_old->coldim() == 1));
    assert((var_map_to_old == 0) || (var_map_to_old->rowdim() + var_del_ind->rowdim() == var_olddim + var_append_dim));
    assert((var_map_to_old == 0) || (var_map_to_old->rowdim() == var_newdim));
    assert((block_map_to_old == 0) || (block_map_to_old->coldim() == 1));
    assert((block_map_to_old == 0) || (block_map_to_old->rowdim() + block_del_ind->rowdim() == block_olddim.dim() + block_append_dim.dim()));
    assert((block_map_to_old == 0) || (block_map_to_old->rowdim() == block_newdim.rowdim()));

    return err;
  }

  // *****************************************************************************
  //                 PSCAffineModification::add_reset_generating_primal
  // *****************************************************************************

  int PSCAffineModification::add_reset_generating_primal(PSCPrimal* new_generating_primal) {
    int err = 0;
    if (append_only) {
      if (cb_out())
        get_out() << "**** ERROR: PSCAffineModification::add_reset_generating_primal(..): append_to_old is set to true, so this operation is not allowed" << std::endl;
      err++;
      return err;
    }
    delete generating_primal;
    reset_primal = true;
    generating_primal = new_generating_primal;
    return 0;
  }

  // *****************************************************************************
  //                 PSCAffineModification::set_skip_extension
  // *****************************************************************************

  int PSCAffineModification::set_skip_extension(bool se) {
    skip_extension = se;
    return 0;
  }

  // *****************************************************************************
  //                        PSCAffineModification::incorporate
  // *****************************************************************************

  int PSCAffineModification::incorporate(const OracleModification& om) {
    int err = 0;

    const PSCAffineModification* m = dynamic_cast<const PSCAffineModification*>(&om);

    if (append_only) {
      if (om.get_old_vardim() != var_olddim) {
        if (cb_out())
          get_out() << "**** ERROR: PSCAffineModification::incorporate(const OracleModfication& om): append_to_old is true; old variable dimension=" << om.get_old_vardim() << " of input om must match old variable dimension=" << var_olddim << " of this" << std::endl;
        err++;
      }
      if (m) {
        if (m->block_olddim.dim() != block_olddim.dim()) {
          if (cb_out())
            get_out() << "**** ERROR: PSCAffineModification::incorporate(.): append_to_old is true; input is PSCAffineModification and old row dimension=" << m->block_olddim.dim() << " of input must match old row dimension=" << block_olddim.dim() << " of this" << std::endl;
          err++;
        } else {
          for (Integer i = 0; i < block_newdim.dim(); i++) {
            if (m->block_olddim(i) != block_newdim(i)) {
              if (cb_out())
                get_out() << "**** ERROR: PSCAffineModification::incorporate(.): append_to_old is true; input is PSCAffineModification and its block " << i << " has old order =" << m->block_olddim(i) << " != " << block_olddim(i) << " = old order of this, both must match" << std::endl;
              err++;
            }
          }
        }
      }

      if (err) {
        if (cb_out())
          get_out() << "**** ERROR: PSCAffineModification::incorporate(.): fatal input errors occured for append_to_old == true, leaving without incorporating anything; return value=" << err << std::endl;
        return err;
      }

      // if ((m)&&(m->reset_primal)){
      //   if (reset_primal) {
      // 	if (cb_out())
      // 	  get_out()<<"**** WARNING: PSCAffineModification::incorporate(.): the current incorporation of a Modification in append mode changes the once more"<<std::endl;
      //   }
      //   reset_primal=true;
      //   delete generating_primal;
      //   if (m->generating_primal){
      // 	generating_primal=dynamic_cast<PSCPrimal*>(m->generating_primal->clone_primal_data());
      // 	assert(generating_primal);
      //   }
      //   else 
      // 	generating_primal=0;
      // }


      SparseCoeffmatMatrix blockapp;
      if (om.get_appended_vardim() == 0) {
        if (m)
          blockapp = m->block_append;
      } else { //om.get_appended_vardim()>0
        if (m)
          var_append.append_columns(m->var_append);
        else
          var_append.append_columns(SparseCoeffmatMatrix(var_append.blockdim(), om.get_appended_vardim()));

        if (block_append_dim.dim() > 0) {
          block_append.append_columns(SparseCoeffmatMatrix(block_append.blockdim(), om.get_appended_vardim()));
        }

        if (m) {
          if (var_append_dim > 0) {
            Indexmatrix ind(Range(0, var_olddim - 1));
            ind.enlarge_below(var_append_dim, -1);
            ind.concat_below(Range(var_olddim, var_olddim + m->var_append_dim - 1));
            blockapp.init(m->block_append_dim, 0);
            blockapp.append_columns(m->block_append, 0, &ind);
          } else {
            blockapp = m->block_append;
          }
        }

        if (var_map_to_old) {
          var_map_to_old->concat_below(Range(var_olddim + var_append_dim, var_olddim + var_append_dim + om.get_appended_vardim() - 1));
          if (m) {
            SparseCoeffmatMatrix tmp(blockapp.blockdim(), 0);
            tmp.append_columns(blockapp, 0, var_map_to_old);
            blockapp = tmp;
          }
        }

        if (var_new_ind)
          var_new_ind->concat_below(Range(var_newdim, var_newdim + om.get_appended_vardim() - 1));
        else
          var_new_ind = new Indexmatrix(Range(var_newdim, var_newdim + om.get_appended_vardim() - 1));

        block_append.append_columns(SparseCoeffmatMatrix(block_append_dim, om.get_appended_vardim()));

        var_append_dim += om.get_appended_vardim();
        var_newdim += om.get_appended_vardim();

      } // endif om.get_appended_vardim()>0

      if ((m) && (m->block_append_dim.dim() > 0)) {

        block_append.append_blocks(blockapp);
        offset_append.append_blocks(m->offset_append);

        if (block_map_to_old)
          block_map_to_old->concat_below(Range(block_olddim.dim() + block_append_dim.dim(), block_olddim.dim() + block_append_dim.dim() + m->block_append_dim.dim() - 1));

        if (block_new_ind)
          block_new_ind->concat_below(Range(block_newdim.dim(), block_newdim.dim() + m->block_append_dim.dim() - 1));
        else
          block_new_ind = new Indexmatrix(Range(block_newdim.dim() + block_append_dim.dim(), block_newdim.dim() + m->block_append_dim.dim() - 1));

        block_append_dim.concat_below(m->block_append_dim);
        block_newdim.concat_below(m->block_append_dim);

      } // endif (m.row_append_dim>0)

      assert((var_append.rowdim() == block_olddim.dim()) && (var_append.coldim() == var_append_dim));
      assert((block_append.rowdim() == block_append_dim.dim()) && (block_append.coldim() == var_newdim));
      assert((var_map_to_old == 0) || (var_map_to_old->coldim() == 1));
      assert((var_map_to_old == 0) || (var_map_to_old->rowdim() + var_del_ind->rowdim() == var_olddim + var_append_dim));
      assert((var_map_to_old == 0) || (var_map_to_old->rowdim() == var_newdim));
      assert((block_map_to_old == 0) || (block_map_to_old->coldim() == 1));
      assert((block_map_to_old == 0) || (block_map_to_old->rowdim() + block_del_ind->rowdim() == block_olddim.dim() + block_append_dim.dim()));
      assert((block_map_to_old == 0) || (block_map_to_old->rowdim() == block_newdim.rowdim()));

      return err;
    } // endif append

    if (om.get_old_vardim() != var_newdim) {
      if (cb_out())
        get_out() << "**** ERROR: PSCAffineModification::incorporate(const OracleModfication& om): old variable dimension=" << om.get_old_vardim() << " of input om must match new variable dimension=" << var_newdim << " of this" << std::endl;
      err++;
    }
    if (m) {
      if (m->block_olddim.dim() != block_newdim.dim()) {
        if (cb_out())
          get_out() << "**** ERROR: PSCAffineModification::incorporate(.): input is PSCAffineModification and old row dimension=" << m->block_olddim.dim() << " of input must match new row dimension=" << block_newdim.dim() << " of this" << std::endl;
        err++;
      } else {
        for (Integer i = 0; i < block_newdim.dim(); i++) {
          if (m->block_olddim(i) != block_newdim(i)) {
            if (cb_out())
              get_out() << "**** ERROR: PSCAffineModification::incorporate(.): input is PSCAffineModification and its block " << i << " has old order =" << m->block_olddim(i) << " != " << block_newdim(i) << " = order after modifiction of this, both must match" << std::endl;
            err++;
          }
        }
      }
    }

    if (err) {
      if (cb_out())
        get_out() << "**** ERROR: PSCAffineModification::incorporate(.): fatal input errors occured, leaving without incorporating anything; return value=" << err << std::endl;
      return err;
    }

    if (m) {
      if (m->var_append_dim > 0)
        err += add_append_vars(m->var_append_dim, &(m->var_append));
      assert(err == 0);

      if (m->var_map_to_old)
        err += add_reassign_vars(*(m->var_map_to_old));
      assert(err == 0);

      if (m->block_append_dim.dim() > 0)
        err += add_append_blocks(m->block_append_dim, &m->offset_append, &m->block_append);
      assert(err == 0);

      if (m->block_map_to_old)
        err += add_reassign_blocks(*(m->block_map_to_old));
      assert(err == 0);

      if (m->reset_primal) {
        reset_primal = true;
        delete generating_primal;
        if (m->generating_primal) {
          generating_primal = dynamic_cast<PSCPrimal*>(m->generating_primal->clone_primal_data());
          assert(generating_primal);
        } else
          generating_primal = 0;
      }
      if (m->skip_extension == true)
        skip_extension = true;
    } else {
      if (om.get_appended_vardim() > 0)
        err += add_append_variables(om.get_appended_vardim());
      assert(err == 0);

      if (om.get_map_to_old_variables()) {
        const GroundsetModification* gm = dynamic_cast<const GroundsetModification*>(&om);
        if (gm) {
          assert(gm->map_to_old_variables());
          err += add_reassign_vars(*(gm->map_to_old_variables()));
        } else {
          err += add_reassign_vars(Indexmatrix(om.get_new_vardim(), 1, om.get_map_to_old_variables()));
        }
      }
      assert(err == 0);
    }

    assert((var_append.rowdim() == block_olddim.dim()) && (var_append.coldim() == var_append_dim));
    assert((block_append.rowdim() == block_append_dim.dim()) && (block_append.coldim() == var_newdim));
    assert((var_map_to_old == 0) || (var_map_to_old->coldim() == 1));
    assert((var_map_to_old == 0) || (var_map_to_old->rowdim() + var_del_ind->rowdim() == var_olddim + var_append_dim));
    assert((var_map_to_old == 0) || (var_map_to_old->rowdim() == var_newdim));
    assert((block_map_to_old == 0) || (block_map_to_old->coldim() == 1));
    assert((block_map_to_old == 0) || (block_map_to_old->rowdim() + block_del_ind->rowdim() == block_olddim.dim() + block_append_dim.dim()));
    assert((block_map_to_old == 0) || (block_map_to_old->rowdim() == block_newdim.rowdim()));

    return err;
  }

  // *****************************************************************************
  //                              PSCAffineModification::apply_to_PSCAffine
  // *****************************************************************************

  int PSCAffineModification::apply_to_PSCAffine(SparseCoeffmatMatrix* offset, SparseCoeffmatMatrix* matrix) const {
    int err = 0;
    if (offset) {
      if ((offset->rowdim() != block_olddim.dim()) ||
        (offset->coldim() > 1) || ((offset->rowdim() != 0) && (offset->coldim() != 1))) {
        err++;
        if (cb_out())
          get_out() << "**** ERROR: PSCAffineModification::apply_to_PSCAffine(..): input offset matrix has size " << offset->rowdim() << " x " << offset->coldim() << " but should have size " << block_olddim.dim() << " x 1" << std::endl;
      } else {
        for (Integer i = 0; i < block_olddim.dim(); i++) {
          if (offset->blockdim(i) != block_olddim(i)) {
            err++;
            if (cb_out())
              get_out() << "**** ERROR: PSCAffineModification::apply_to_PSCAffine(..): block " << i << " of input offset has order " << offset->blockdim(i) << " != " << block_olddim(i) << ", which is the block size given by block_olddim(" << i << ")" << std::endl;
          }
        }
      }
    }
    if (matrix) {
      if ((matrix->rowdim() != block_olddim.dim()) ||
        (((matrix->rowdim() != 0) || (matrix->coldim() != 0)) && (matrix->coldim() != var_olddim))) {
        err++;
        if (cb_out())
          get_out() << "**** ERROR: PSCAffineModification::apply_to_PSCAffine(..): input matrix has size " << matrix->rowdim() << " x " << matrix->coldim() << " but should have size " << block_olddim.dim() << " x " << var_olddim << std::endl;
      } else {
        for (Integer i = 0; i < block_olddim.dim(); i++) {
          if (matrix->blockdim(i) != block_olddim(i)) {
            err++;
            if (cb_out())
              get_out() << "**** ERROR: PSCAffineModification::apply_to_PSCAffine(..): block " << i << " of input matrix has order " << matrix->blockdim(i) << " != " << block_olddim(i) << ", which is the block size given by block_olddim(" << i << ")" << std::endl;
          }
        }
      }
    }

    if (err) {
      if (cb_out())
        get_out() << "**** ERROR: PSCAffineModification::apply_to_PSCAffine(..): fatal input errors occured, leaving without applying anything; return value=" << err << std::endl;
      return err;
    }

    if (var_append_dim > 0) {
      if ((matrix) && (matrix->append_columns(var_append))) {
        err++;
      }
    }
    if (var_map_to_old) {
      if ((matrix) && (matrix->reassign_columns(*var_map_to_old))) {
        err++;
      }
    }
    if (block_append_dim.dim() > 0) {
      if ((matrix) && (matrix->append_blocks(block_append))) {
        err++;
      }
      if ((offset) && (offset->append_blocks(offset_append))) {
        err++;
      }
    }
    if (block_map_to_old) {
      if ((matrix) && (matrix->reassign_blocks(*block_map_to_old))) {
        err++;
      }
      if ((offset) && (offset->reassign_blocks(*block_map_to_old))) {
        err++;
      }
    }

    return err;
  }



  // *****************************************************************************
  //                              PSCAffineModification::no_modification
  // *****************************************************************************

  bool PSCAffineModification::no_modification() const {
    if (
      (var_append_dim == 0) &&
      (var_map_to_old == 0) &&
      (block_append_dim.dim() == 0) &&
      (block_map_to_old == 0) &&
      (reset_primal == false)
      )
      return true;
    return false;
  }

  // *****************************************************************************
  //                            PSCAffineModification::set_append_to_old
  // *****************************************************************************

  int PSCAffineModification::set_append_to_old(bool ao) {
    if (ao) {
      if (!(
        (var_map_to_old == 0) &&
        (block_map_to_old == 0)
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
  //                    PSCAffineModification::deleted_variables_are_zero
  // *****************************************************************************

  bool PSCAffineModification::deleted_variables_are_zero(const Matrix& oldpoint, const SparseCoeffmatMatrix& oldmat) const {
    if ((oldpoint.rowdim() != var_olddim) || (oldpoint.coldim() != 1)) {
      if (cb_out()) {
        get_out() << "**** ERROR PSCAffineModification::deleted_variables_are_zero(.): point has dim=" << oldpoint.rowdim() << " x " << oldpoint.coldim() << " but should have dim=" << var_olddim << " x 1" << std::endl;
      }
      std::abort();
    }
    if ((oldmat.coldim() != var_olddim) || (oldmat.rowdim() != block_olddim.dim())) {
      if (cb_out()) {
        get_out() << "**** ERROR PSCAffineModification::deleted_variables_are_zero(.): matrix has dim=" << oldmat.rowdim() << " x " << oldmat.coldim() << " but should have dim=" << var_olddim << " x " << block_olddim.dim() << std::endl;
      }
      std::abort();
    }
    if (var_del_ind) {
      Integer n = var_del_ind->rowdim();
      for (Integer i = 0; i < n; i++) {
        Integer ind = (*var_del_ind)(i);
        if ((oldpoint(ind) != 0.) && (oldmat.column(ind))) {
          return false;
        }
      }
    }
    return true;
  }

  // *****************************************************************************
  //                    PSCAffineModification::new_variables_are_zero
  // *****************************************************************************

  bool PSCAffineModification::new_variables_are_zero(const Matrix& newpoint, const SparseCoeffmatMatrix& newmat) const {
    if ((newpoint.rowdim() != var_newdim) || (newpoint.coldim() != 1)) {
      if (cb_out()) {
        get_out() << "**** ERROR PSCAffineModification::new_variables_are_zero(.): point has dim=" << newpoint.rowdim() << " x " << newpoint.coldim() << " but should have dim=" << var_newdim << " x 1" << std::endl;
      }
      std::abort();
      if ((newmat.coldim() != var_newdim) || (newmat.rowdim() != block_newdim.dim())) {
        if (cb_out()) {
          get_out() << "**** ERROR PSCAffineModification::deleted_variables_are_zero(.): matrix has dim=" << newmat.rowdim() << " x " << newmat.coldim() << " but should have dim=" << var_newdim << " x " << block_newdim.dim() << std::endl;
        }
        std::abort();
      }
    }
    if (var_new_ind) {
      Integer n = var_new_ind->rowdim();
      for (Integer i = 0; i < n; i++) {
        Integer ind = (*var_new_ind)(i);
        assert((var_map_to_old == 0) || ((*var_map_to_old)(ind) == i + var_olddim));
        if ((newpoint(ind) != 0.) && (newmat.column(ind))) {
          return false;
        }
      }
    }
    return true;
  }

  // *****************************************************************************
  //                    PSCAffineModification::mapped_variables_are_equal
  // *****************************************************************************

  bool PSCAffineModification::mapped_variables_are_equal(const Matrix& newpoint, const Matrix& oldpoint) const {
    int err = 0;
    if ((oldpoint.rowdim() != var_olddim) || (oldpoint.coldim() != 1)) {
      if (cb_out()) {
        get_out() << "**** ERROR PSCAffineModification::deleted_variables_are_equal(.): point has dim=" << oldpoint.rowdim() << " x " << oldpoint.coldim() << " but should have dim=" << var_olddim << " x 1" << std::endl;
      }
      err++;
    }
    if ((newpoint.rowdim() != var_newdim) || (newpoint.coldim() != 1)) {
      if (cb_out()) {
        get_out() << "**** ERROR PSCAffineModification::new_variables_are_equal(.): point has dim=" << oldpoint.rowdim() << " x " << newpoint.coldim() << " but should have dim=" << var_newdim << " x 1" << std::endl;
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

