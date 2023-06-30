/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/MatrixCBSolver.cxx
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



#include "MatrixCBSolver.hxx"
#include "ModificationTreeData.hxx"
#include "BundleSolver.hxx"
#include "LPGroundset.hxx"
#include "LPGroundsetModification.hxx"
#include "FunctionObjectModification.hxx"
#include "SumBundle.hxx"
#include "SumModel.hxx"
#include "NNCModel.hxx"
#include "PSCModel.hxx"
#include "SOCModel.hxx"
#include "BoxModel.hxx"
#include "BundleIdProx.hxx"
#include "BundleDiagonalTrustRegionProx.hxx"
#include "SumBundleParameters.hxx"
#include "SumModelParameters.hxx"

#include <algorithm>
#include <map>

//------------------------------------------------------------


using namespace CH_Tools;
using namespace CH_Matrix_Classes;

namespace ConicBundle {

  const double CB_plus_infinity = max_Real / 10000.;
  const double CB_minus_infinity = -CB_plus_infinity;

  ModifiableOracleObject::~ModifiableOracleObject() {
  }

  int ModifiableOracleObject::apply_modification
  (
    const OracleModification& /* oracle_modification */,
    const CH_Matrix_Classes::Matrix* /* new_center */,
    const CH_Matrix_Classes::Matrix* /* old_center */,
    bool& /* discard_objective_in_center */,
    bool& /* discard_model */,
    bool& /* discard_aggregates */,
    MinorantExtender*& /* minorant_extender */
  ) {
    return 1;
  }

  //------------------------------------------------------------
  // Wrapper for transforming FunctionOracle to MatrixFunctionOracle
  //------------------------------------------------------------


  class FunctionOracleWrapper : public MatrixFunctionOracle, public CBout {
  private:
    FunctionOracle& oracle;
  public:
    FunctionOracleWrapper(FunctionOracle& o) : oracle(o) {
    }
    ~FunctionOracleWrapper() {
    }


    /** see MatrixFunctionOracle for explanations */
    int evaluate(const  Matrix& y, double relprec, double& objective_value,
      std::vector<Minorant*>& minorants,
      PrimalExtender*& primal_extender) {
      //---- evaluate
      return oracle.evaluate(y.get_store(), relprec, objective_value,
        minorants, primal_extender);

    }

    int apply_modification(
      const OracleModification& oracle_modification,
      const CH_Matrix_Classes::Matrix* new_center,
      const CH_Matrix_Classes::Matrix* old_center,
      bool& discard_objective_in_center,
      bool& discard_model,
      bool& discard_aggregates,
      MinorantExtender*& minorant_extender
    ) {
      return oracle.apply_modification(oracle_modification,
        new_center ? new_center->get_store() : 0,
        old_center ? old_center->get_store() : 0,
        discard_objective_in_center,
        discard_model,
        discard_aggregates,
        minorant_extender);
    }

    bool check_correctness() const {
      return oracle.check_correctness();
    }

  };


  //------------------------------------------------------------
  // Problem structure and modification class
  //------------------------------------------------------------

  //------------------------------------------------------------
  // Data class
  //------------------------------------------------------------
  /** in order to keep the interface 'clean' the interface's data is separated into an extra class */
  class MatrixCBSolverData : public CBout {
    friend class MatrixCBSolver;

    MatrixCBSolverData(const MatrixCBSolverData&); // blocked
    MatrixCBSolverData& operator= (const MatrixCBSolverData&); // blocked

  public:

    //----------------------------------------
    // data
    BundleSolver solver;
    LPGroundset groundset;
    LPGroundsetModification* gs_modif;
    ModificationTreeData* root;
    FunctionMap  fun_model;
    Clock        myclock;
    std::vector<FunctionOracleWrapper*> wrappers;

    void set_cbout(const CBout* cb, int incr = -1) {
      CBout::set_cbout(cb, incr);
      solver.set_cbout(this, 0);
      groundset.set_cbout(this, 0);
      gs_modif->set_cbout(this, 0);
      if (root)
        root->set_cbout(this, 0);
      for (FunctionMap::iterator it = fun_model.begin();
        it != fun_model.end();
        ++it) {
        it->second->set_cbout(this, 0);
      }
    }

    void set_out(std::ostream* o = 0, int pril = 1) {
      CBout::set_out(o, pril);
      set_cbout(this, 0);
    }


    int clear_modifications() {
      delete gs_modif;
      gs_modif = dynamic_cast<LPGroundsetModification*>(groundset.start_modification());
      assert(gs_modif);
      if (root) {
        if (root->clear_subtree_modification()) {
          if (cb_out()) {
            get_out() << "**** ERROR MatrixCBSolverData::clear_modifications(): root->clear_subtree_modification() failed" << std::endl;
          }
          return 1;
        }
      }
      return 0;
    }

    int update_tree_modification(Integer add_dim, const Indexmatrix* map_to_old, const FunObjModMap* funmodmap) {
      if (root == 0)
        return 0;

      int retval = 0;
      FunObjModMap extfunmodmap;
      if (funmodmap) {
        extfunmodmap = *funmodmap;
        for (FunObjModMap::const_iterator it = funmodmap->begin(); it != funmodmap->end(); it++) {
          FunctionMap::iterator fit = fun_model.find(it->first);
          if (fit == fun_model.end()) {
            if (cb_out())
              get_out() << "**** WARNING: MatrixCBSolverData::update_tree_modifications(...): the map giving affected functions lists functions that are not part of the current setting" << std::endl;
            retval++;
          } else if (((add_dim > 0) || (map_to_old)) && (fit->second->add_parents_to_map(extfunmodmap))) {
            if (cb_out())
              get_out() << "**** ERROR: MatrixCBSolverData::update_tree_modifications(...): adding parent functions for a function of the list of affected functions to this list failed" << std::endl;
            retval++;
          }
        }
        funmodmap = &extfunmodmap;
      }

      if (retval == 0) {
        retval = root->update_subtree_modification(add_dim, map_to_old, funmodmap);
        if (retval) {
          if (cb_out())
            get_out() << "**** ERROR: MatrixCBSolverData::update_tree_modifications(...): updating the modifications of the functions failed and returned " << retval << std::endl;
        }
      }

      if (retval) {
        if (cb_out())
          get_out() << "****      MatrixCBSolverData::update_tree_modifications(...): errors occured with partial modification updates, clearing all pending modifications" << std::endl;
        clear_modifications();
      }

      return retval;
    }


    int modifications_performed() {
      delete gs_modif;
      gs_modif = dynamic_cast<LPGroundsetModification*>(groundset.start_modification());
      assert(gs_modif);
      if (root) {
        if (root->subtree_modification_performed()) {
          if (cb_out()) {
            get_out() << "**** ERROR MatrixCBSolverData::clear_modifications(): root->subtree_modification_performed() failed" << std::endl;
          }
          return 1;
        }
      }
      return 0;
    }

    int apply_modifications() {
      FunObjModMap fomm;
      if ((root) && (root->collect_subtree_modification(fomm))) {
        if (cb_out()) {
          get_out() << "**** ERROR MatrixCBSolverData::apply_modifications(): root->collect_subtree_modification() failed" << std::endl;
        }
        return 1;
      }
      if ((gs_modif->no_modification()) && (fomm.size() == 0))
        return 0;
      if (solver.apply_modification(*gs_modif, fomm)) {
        if (cb_out()) {
          get_out() << "**** ERROR MatrixCBSolverData::apply_modifications(): solver.apply_modification failed" << std::endl;
        }
        return 1;
      }
      modifications_performed();
      return 0;
    }

    void clear(void) {
      if ((gs_modif) && (apply_modifications())) {
        if (cb_out()) {
          get_out() << "**** ERROR MatrixCBSolverData::clear(): executing pending modifications before clear() failed" << std::endl;
        }
      }
      groundset.clear(0, 0);
      delete gs_modif;
      gs_modif = dynamic_cast<LPGroundsetModification*>(groundset.start_modification());
      assert(gs_modif);
      if (root) {
        root->delete_descendants(fun_model);
        fun_model.clear();
        delete root;
        root = 0;
      }
      for (unsigned int i = 0; i < wrappers.size(); i++)
        delete wrappers[i];
      wrappers.clear();

      solver.initialize(&groundset, 0);
      solver.set_clock(myclock);
      myclock.start();
    };

    void set_defaults(void) {
      solver.set_defaults();
    }


    ///
    MatrixCBSolverData(const CBout* cb, int incr = -1) :CBout(cb, incr), gs_modif(0), root(0) {
      solver.set_cbout(this, 0);
      groundset.set_cbout(this, 0);
      clear();
      set_defaults();
    }

    ///
    ~MatrixCBSolverData() {
      clear();
      delete gs_modif;
      gs_modif = 0;
    };


  };

  //------------------------------------------------------------
  // CBmethod implementation - mostly just wrapped to CBmethodData
  //------------------------------------------------------------

  //--------------------
  MatrixCBSolver::MatrixCBSolver(std::ostream* out, int pril) :CBout(out, pril) {
    data_ = new MatrixCBSolverData(this, 0);
    assert(data_);
    data_->set_cbout(this, 0);
  }

  //--------------------
  MatrixCBSolver::MatrixCBSolver(const CBout* cb, int incr) :CBout(cb, incr) {
    data_ = new MatrixCBSolverData(this, 0);
    assert(data_);
    data_->set_cbout(this, 0);
  }

  //--------------------
  MatrixCBSolver::~MatrixCBSolver() {
    assert(data_);
    delete data_;
  }

  //--------------------
  void MatrixCBSolver::clear() {
    assert(data_);
    data_->clear();
  }

  //--------------------
  void MatrixCBSolver::set_defaults() {
    assert(data_);
    data_->set_defaults();
  }

  //--------------------
  int MatrixCBSolver::init_problem(int dim,
    const Matrix* lbounds,
    const Matrix* ubounds,
    const Matrix* startval,
    const Matrix* costs,
    Real offset) {
    assert(data_);
    clear();
    int err = 0;
    if (data_->gs_modif->add_append_vars(dim, lbounds, ubounds, 0, startval, costs)) {
      if (data_->cb_out())
        data_->get_out() << "**** ERROR: MatrixCBSolver::init_problem(...): "
        << " setting initial dimension failed" << std::endl;
      err++;
    }
    if (data_->gs_modif->add_offset(offset)) {
      if (data_->cb_out())
        data_->get_out() << "**** ERROR: MatrixCBSolver::init_problem(...): "
        << " setting initial offset failed" << std::endl;
      err++;
    }
    return 0;
  }

  //--------------------
  int MatrixCBSolver::add_function(FunctionObject& function,
    Real fun_factor,
    FunctionTask fun_task,
    AffineFunctionTransformation* aft,
    bool dyn_arg) {
    assert(data_);

    if (data_->cb_out(10)) {
      data_->get_out() << "\n  entering  MatrixCBSolver::add_function" << std::endl;
    }

    if (data_->fun_model.find(&function) != data_->fun_model.end()) {
      if (data_->cb_out())
        data_->get_out() << "**** ERROR: MatrixCBSolver::add_function(...): "
        << "function already added and cannot be added twice; add a copy instead" << std::endl;
      return 1;
    }

    int retval = 0;
    if (fun_factor <= 0.) {
      if (data_->cb_out())
        data_->get_out() << "**** ERROR: MatrixCBSolver::add_function(...): "
        << "function factor=" << fun_factor << " but has to be strictly positive" << std::endl;
      retval++;
    }

    if ((dynamic_cast<AffineFunctionTransformation*>(&function))
      && ((fun_factor != 1.) || (fun_task != ObjectiveFunction))) {
      if (data_->cb_out())
        data_->get_out() << "**** ERROR: MatrixCBSolver::add_function(...): the call adds an AffineFunctionTransformation as function, but for these only function factor == 1. (it is =" << fun_factor << ") with function task ObjectiveFunction (it is =" << fun_task << ") are allowed" << std::endl;
      retval++;
    }

    if ((aft) && (data_->gs_modif->new_vardim() != aft->from_dim())) {
      if (data_->cb_out())
        data_->get_out() << "**** ERROR: MatrixCBSolver::add_function(...): "
        << "the input dimension of the affine function transformation =" << aft->from_dim() << " doesn not match the current problem dimension = " << data_->gs_modif->new_vardim() << std::endl;
      retval++;
    }


    //execute accmulated changes
    if ((retval == 0) && (data_->apply_modifications())) {
      if (data_->cb_out()) {
        data_->get_out() << "**** ERROR MatrixCBSolver::add_function(.): apply_modifications failed" << std::endl;
      }
      retval++;
    }

    SumBlockModel* sbm = 0;
    FunctionOracleWrapper* ow = 0;

    //determine the right sumblockmodel, equip it with the aft if necessary and add it
    if (retval == 0) {
      bool unknown_derivation = true;
      if (dynamic_cast<FunctionOracle*>(&function)) {
        unknown_derivation = false;
        FunctionOracleWrapper* ow = new FunctionOracleWrapper(*dynamic_cast<FunctionOracle*>(&function));
        if (ow == 0) {
          if (data_->cb_out())
            data_->get_out() << "**** ERROR: MatrixCBSolver::add_function(...): construction of wrapper interface for function failed" << std::endl;
          retval++;
        } else {
          data_->wrappers.push_back(ow);
          sbm = new NNCModel(ow, fun_factor, fun_task, this);
        }
      } else if (dynamic_cast<MatrixFunctionOracle*>(&function)) {
        unknown_derivation = false;
        sbm = new NNCModel(dynamic_cast<MatrixFunctionOracle*>(&function), fun_factor, fun_task, this);
      } else if (dynamic_cast<PSCOracle*>(&function)) {
        unknown_derivation = false;
        sbm = new PSCModel(dynamic_cast<PSCOracle*>(&function), fun_factor, fun_task, this);
      } else if (dynamic_cast<SOCOracle*>(&function)) {
        unknown_derivation = false;
        sbm = new SOCModel(dynamic_cast<SOCOracle*>(&function), fun_factor, fun_task, this);
      } else if (dynamic_cast<BoxOracle*>(&function)) {
        unknown_derivation = false;
        sbm = new BoxModel(dynamic_cast<BoxOracle*>(&function), fun_factor, fun_task, this);
      } else if (dynamic_cast<AffineFunctionTransformation*>(&function)) {
        unknown_derivation = false;
        sbm = new AFTModel(0, dynamic_cast<AffineFunctionTransformation*>(&function), 0, false, this);
      }

      if (unknown_derivation) {
        if (data_->cb_out())
          data_->get_out() << "**** ERROR: MatrixCBSolver::add_function(...): unknown derivation of FunctionObject" << std::endl;
        retval++;
      } else if (sbm == 0) {
        if ((retval == 0) && (data_->cb_out()))
          data_->get_out() << "**** ERROR: MatrixCBSolver::add_function(...): construction of submodel for function failed" << std::endl;
        retval++;
      } else if ((aft) && (sbm->initialize_aft(aft))) {
        assert(retval == 0);
        if (data_->cb_out())
          data_->get_out() << "**** ERROR: MatrixCBSolver::add_function(...): could not initialize the affine function transformation" << std::endl;
        retval++;
      } else {
        if (data_->fun_model.size() == 1) {
          ModificationTreeData* old_root = data_->root;
          SumBlockModel* sumbl = new SumModel;
          data_->root = new ModificationTreeData(sumbl->get_oracle_object(), 0, sumbl, data_->gs_modif->new_vardim(), -1, 0, this);
          if (data_->root->add_child(old_root)) {
            if (data_->cb_out())
              data_->get_out() << "**** ERROR: MatrixCBSolver::add_function(...): in forming the sum of the first two functions adding the first failed" << std::endl;
            retval++;
          }
          if (retval) {
            delete data_->root;
            data_->root = old_root;
          } else {
            data_->fun_model[sumbl->get_oracle_object()] = data_->root;
            data_->solver.set_model(data_->root->get_model());
          }
        }
      }

      if (retval) {
        delete sbm;
        if (ow)
          delete ow;
      } else {
        ModificationTreeData* mtd = new ModificationTreeData(&function, ow, sbm, data_->gs_modif->new_vardim(), dyn_arg ? -1 : (((aft == 0) || (aft->to_dim() < 0)) ? data_->gs_modif->new_vardim() : aft->to_dim()), aft, this);

        if (data_->root == 0) {
          assert(data_->fun_model.size() == 0);
          data_->root = mtd;
          data_->solver.set_model(data_->root->get_model());
        } else {
          if (data_->root->add_child(mtd)) {
            if (data_->cb_out())
              data_->get_out() << "**** ERROR: MatrixCBSolver::add_function(...): in forming the sum of the first two functions adding the first failed" << std::endl;
            retval++;
          }
        }

        if (retval)
          delete mtd;
        else
          data_->fun_model[&function] = mtd;
      }
    }


    if (data_->cb_out(10)) {
      data_->get_out() << "\n  leaving  MatrixCBSolver::append_variables with return value" << retval << std::endl;
    }

    return retval;


    /*
    BaseSOCOracle* soc=dynamic_cast<BaseSOCOracle*>(&function);
    if (soc) {
      SocModel* socp=new SocModel(soc);
      if (socp==0) {
  if (data_->cb_out())
          data_->get_out()<< "**** ERROR: MatrixCBSolver::add_function(...): construction of subproblen for function failed" <<std::endl;
  return 1;
      }

      data_->init=0;
      if (data_->model.add_function(socp)){
  if (data_->cb_out())
          data_->get_out()<< "**** ERROR: MatrixCBSolver::add_function(...): could not add function to SumModel " <<std::endl;
  delete socp;
  return 1;
      }

      data_->fun_model[ &function ] = socp;
      return 0;
    }

    BaseConeOracle* con=dynamic_cast<BaseConeOracle*>(&function);
    if (con) {
      ConeModel* conp=new ConeModel(con);
      if (conp==0) {
  if (data_->cb_out())
          data_->get_out()<< "**** ERROR: MatrixCBSolver::add_function(...): construction of subproblen for function failed" <<std::endl;
  return 1;
      }

      data_->init=0;
      if (data_->model.add_function(conp)){
  if (data_->cb_out())
          data_->get_out()<< "**** ERROR: MatrixCBSolver::add_function(...): could not add function to SumModel " <<std::endl;
  delete conp;
  return 1;
      }

      data_->fun_model[ &function ] = conp;
      return 0;
    }
    */

  }


  //----------------------------------------
  // append new variables (always in last postions in this order)
  int MatrixCBSolver::append_variables(int n_append,
    const Matrix* lbounds,
    const Matrix* ubounds,
    const Sparsemat* columns,
    const Matrix* startval,
    const Matrix* costs,
    const FunObjModMap* funmodmap) {
    assert(data_);

    if (data_->cb_out(10)) {
      data_->get_out() << "\n  entering  MatrixCBSolver::append_variables" << std::endl;
    }

    int retval = 0;
    if ((n_append == 0) && ((funmodmap == 0) || (funmodmap->size() == 0)))
      return 0;
    if (n_append < 0) {
      if (data_->cb_out())
        data_->get_out() << "**** ERROR: MatrixCBSolver::append_variables(...): n_append <0" << std::endl;
      retval++;
    } else {
      if ((lbounds != 0) && (lbounds->dim() != n_append)) {
        if (data_->cb_out())
          data_->get_out() << "**** ERROR: MatrixCBSolver::append_variables(...): lower bounds vector does not match added dimension" << std::endl;
        retval++;
      }
      if ((ubounds != 0) && (ubounds->dim() != n_append)) {
        if (data_->cb_out())
          data_->get_out() << "**** ERROR: MatrixCBSolver::append_variables(...): upper bounds vector does not match added dimension" << std::endl;
        retval++;
      }
      if ((columns != 0) && (columns->coldim() != n_append)) {
        if (data_->cb_out())
          data_->get_out() << "**** ERROR: MatrixCBSolver::append_variables(...): number of constraint columns does not match added dimension" << std::endl;
        retval++;
      }
      if ((columns != 0) && (columns->rowdim() != data_->gs_modif->new_rowdim())) {
        if (data_->cb_out())
          data_->get_out() << "**** ERROR: MatrixCBSolver::append_variables(...): number of rows of constraint columns =" << columns->rowdim() << " does not match number of constraints =" << data_->gs_modif->new_rowdim() << std::endl;
        retval++;
      }
      if ((startval != 0) && (startval->dim() != n_append)) {
        if (data_->cb_out())
          data_->get_out() << "**** ERROR: MatrixCBSolver::append_variables(...): vector of starting values does not match added dimension" << std::endl;
        retval++;
      }
      if ((costs != 0) && (costs->dim() != n_append)) {
        if (data_->cb_out())
          data_->get_out() << "**** ERROR: MatrixCBSolver::append_variables(...): costs vector does not match added dimension" << std::endl;
        retval++;
      }
    }

    //check correctnes of values of lbounds and ubounds and warn about startval violations
    if ((retval == 0) && ((lbounds != 0) || (ubounds != 0) || (startval != 0))) {
      for (Integer i = 0; i < n_append; i++) {
        if (lbounds) {
          if ((*lbounds)[i] > CB_plus_infinity) {
            if (data_->cb_out())
              data_->get_out() << "**** ERROR: MatrixCBSolver::append_variables(...): lower bound of coordinate " << i << " exceeds plus_infinity: " << (*lbounds)[i] << std::endl;
            retval++;
          }
          if ((*lbounds)[i] == CB_plus_infinity) {
            if (data_->cb_out())
              data_->get_out() << "**** WARNING: MatrixCBSolver::append_variables(...): lower bound of coordinate " << i << " equals plus_infinity: " << (*lbounds)[i] << std::endl;
          }
          if ((*lbounds)[i] < CB_minus_infinity) {
            if (data_->cb_out())
              data_->get_out() << "**** WARNING: MatrixCBSolver::append_variables(...): lower bound of coordinate " << i << " is smaller than minus_infinity: " << (*lbounds)[i] << std::endl;
          }
        } //endif Lbounds

        if (ubounds) {
          if ((*ubounds)[i] < CB_minus_infinity) {
            if (data_->cb_out())
              data_->get_out() << "**** ERROR: MatrixCBSolver::append_variables(...): upper bound of coordinate " << i << " exceeds minus_infinity: " << (*ubounds)[i] << std::endl;
            retval++;
          }
          if ((*ubounds)[i] == CB_minus_infinity) {
            if (data_->cb_out())
              data_->get_out() << "**** WARNING: MatrixCBSolver::append_variables(...): upper bound of coordinate " << i << " equals minus_infinity: " << (*ubounds)[i] << std::endl;
          }
          if ((*ubounds)[i] > CB_plus_infinity) {
            if (data_->cb_out())
              data_->get_out() << "**** WARNING: MatrixCBSolver::append_variables(...): upper bound of coordinate " << i << " exceeds plus_infinity: " << (*ubounds)[i] << std::endl;
          }
          if ((lbounds) && ((*ubounds)[i] < (*lbounds)[i])) {
            if (data_->cb_out())
              data_->get_out() << "**** ERROR: MatrixCBSolver::append_variables(...): lower bound " << (*lbounds)[i] << " of coordinate " << i << " exceeds upper bound " << (*ubounds)[i] << std::endl;
            retval++;
          }
        } //endif ubounds

        if (startval) {
          if (lbounds) {
            if ((*startval)[i] < (*lbounds)[i]) {
              if (data_->cb_out())
                data_->get_out() << "**** WARNING: MatrixCBSolver::append_variables(...): starting value[" << i << "]=" << (*startval)[i] << " is below lower bound [" << i << "]= " << (*lbounds)[i] << std::endl;
            }
          } else if ((*startval)[i] < CB_minus_infinity) {
            if (data_->cb_out())
              data_->get_out() << "**** WARNING: MatrixCBSolver::append_variables(...): starting value[" << i << "]=" << (*startval)[i] << " is below minus_infinity" << std::endl;
          }
          if (ubounds) {
            if ((*startval)[i] > (*ubounds)[i]) {
              if (data_->cb_out())
                data_->get_out() << "**** WARNING: MatrixCBSolver::append_variables(...): starting value[" << i << "]=" << (*startval)[i] << " exceeds upper bound [" << i << "]= " << (*ubounds)[i] << std::endl;
            }
          } else if ((*startval)[i] > CB_plus_infinity) {
            if (data_->cb_out())
              data_->get_out() << "**** WARNING: MatrixCBSolver::append_variables(...): starting value[" << i << "]=" << (*startval)[i] << " exceeds plus_infinity" << std::endl;
          }
        }
      } //endfor
    } //endif lbounds or ubounds or startval

    //set the modification of the ground set
    if (retval == 0) {

      retval = data_->gs_modif->add_append_vars(n_append, lbounds, ubounds, columns, startval, costs);
      if ((retval) && (data_->cb_out())) {
        data_->get_out() << "**** ERROR: MatrixCBSolver::append_variables(...): setting ground set changes failed and returned " << retval << std::endl;
      }
    }

    //set the modifications of the functions
    if ((retval == 0) && (data_->update_tree_modification(n_append, 0, funmodmap))) {
      if (data_->cb_out())
        data_->get_out() << "**** ERROR: MatrixCBSolver::append_variables(...): update_tree_modification failed" << std::endl;
      retval++;
    }

    //
    if (data_->cb_out(10)) {
      data_->get_out() << "\n  leaving  MatrixCBSolver::append_variables with return value " << retval << std::endl;
    }

    return retval;
  }

  //----------------------------------------
  // delete variables
  int MatrixCBSolver::delete_variables(const Indexmatrix& del_indices,
    Indexmatrix& map_to_old,
    const FunObjModMap* funmodmap) {
    assert(data_);

    if (data_->cb_out(10)) {
      data_->get_out() << "\n  entering  MatrixCBSolver::delete_variables" << std::endl;
    }

    int retval = 0;
    if (del_indices.dim() == 0) {
      map_to_old.init(Range(0, data_->gs_modif->new_vardim() - 1));
    } else {
      Indexmatrix del_vec;
      sortindex(del_indices, del_vec);
      del_vec = del_indices(del_vec);
      if ((del_vec(0) < 0) || (del_vec(del_vec.dim() - 1) >= data_->gs_modif->new_vardim())) {
        if (data_->cb_out())
          data_->get_out() << "**** ERROR: MatrixCBSolver::delete_function(...): variable indices for deletion out of range" << std::endl;
        retval++;
      }
      for (Integer i = 1; i < del_vec.dim(); i++) {
        if (del_vec(i - 1) == del_vec(i)) {
          if (data_->cb_out())
            data_->get_out() << "**** ERROR: MatrixCBSolver::delete_function(...): multiple indices for deletion, e.g. " << del_vec(i) << std::endl;
          retval++;
          break;
        }
      }
      if (retval == 0) {
        map_to_old.init(Range(0, data_->gs_modif->new_vardim() - 1));
        map_to_old.delete_rows(del_vec);
      }
    }

    if ((retval == 0) && ((del_indices.dim() != 0) || ((funmodmap != 0) && (funmodmap->size() > 0)))) {
      retval = reassign_variables(map_to_old);
    }

    if (data_->cb_out(10)) {
      data_->get_out() << "\n  leaving  MatrixCBSolver::delete_variables" << std::endl;
    }

    return retval;
  }


  //----------------------------------------
  // reassign variables
  int MatrixCBSolver::reassign_variables(const Indexmatrix& avec,
    const FunObjModMap* funmodmap) {
    assert(data_);

    if (data_->cb_out(10)) {
      data_->get_out() << "\n  entering  MatrixCBSolver::reassign_variables" << std::endl;
    }

    int retval = 0;
    if ((avec.dim() > 0) && ((min(avec) < 0) || (max(avec) >= data_->gs_modif->new_vardim()))) {
      if (data_->cb_out())
        data_->get_out() << "**** ERROR: MatrixCBSolver::reassign_variables(...): variable indices out of range" << std::endl;
      retval++;
    }

    //set the modifications of the ground set
    if ((retval == 0) && (data_->gs_modif->add_reassign_vars(avec))) {
      if (data_->cb_out())
        data_->get_out() << "**** ERROR: MatrixCBSolver::reassign_variables(...): reassigning groundset variables failed" << std::endl;
      retval++;
    }

    //set the modifications of the functions
    if ((retval == 0) && (data_->update_tree_modification(0, &avec, funmodmap))) {
      if (data_->cb_out())
        data_->get_out() << "**** ERROR: MatrixCBSolver::reassign_variables(...): update_tree_modification failed" << std::endl;
      retval++;
    }


    if (data_->cb_out(10)) {
      data_->get_out() << "\n  leaving  MatrixCBSolver::reassign_variables with return value" << retval << std::endl;
    }

    return retval;
  }

  int MatrixCBSolver::append_constraints(Integer n_append,
    const Sparsemat* append_rows,
    const Matrix* append_rhslb,
    const Matrix* append_rhsub) {
    assert(data_);

    if (data_->cb_out(10)) {
      data_->get_out() << "\n  entering  MatrixCBSolver::append_constraints" << std::endl;
    }

    int retval = 0;
    if (n_append == 0)
      return 0;
    if (n_append < 0) {
      if (data_->cb_out())
        data_->get_out() << "**** ERROR: MatrixCBSolver::append_constraints(...): n_append <0" << std::endl;
      retval++;
    } else if (data_->gs_modif->add_append_rows(n_append, append_rows, append_rhslb, append_rhsub)) {
      if (data_->cb_out())
        data_->get_out() << "**** ERROR: MatrixCBSolver::append_constraints(...): add_append_rows failed for the groundset modification" << std::endl;
      retval++;
    }

    if (data_->cb_out(10)) {
      data_->get_out() << "\n  leaving  MatrixCBSolver::append_constraints with return value" << retval << std::endl;
    }

    return retval;
  }


  //----------------------------------------
  int MatrixCBSolver::set_lower_bound(int i, double lb) {
    assert(data_);
    if ((i < 0) || (i >= data_->gs_modif->new_vardim()))
      return 1;

    if (lb < CB_minus_infinity)
      return data_->gs_modif->add_set_lb(i, CB_minus_infinity);
    return data_->gs_modif->add_set_lb(i, lb);
  }

  //----------------------------------------
  int MatrixCBSolver::set_upper_bound(int i, double ub) {
    assert(data_);
    if ((i < 0) || (i >= data_->gs_modif->new_vardim()))
      return 1;

    if (ub > CB_plus_infinity)
      return data_->gs_modif->add_set_ub(i, CB_plus_infinity);
    return data_->gs_modif->add_set_ub(i, ub);
  }


  //--------------------
  int MatrixCBSolver::solve(int maxsteps, bool stop_at_descent_steps) {
    assert(data_);

    if (data_->cb_out(10)) {
      data_->get_out() << "\n  entering  MatrixCBSolver::solve" << std::endl;
    }


    if (data_->apply_modifications()) {
      if (data_->cb_out()) {
        data_->get_out() << "**** ERROR MatrixCBSolver::solve(.): apply_modifications failed" << std::endl;
      }
      return 1;
    }

    int status = data_->solver.solve(maxsteps, stop_at_descent_steps);
    if (status) {
      if (data_->cb_out()) {
        data_->get_out() << "**** ERROR MatrixCBSolver::solve(.): solve failed and returned " << status << std::endl;
      }
    }

    return status;
  }

  //--------------------
  int MatrixCBSolver::termination_code() const {
    assert(data_);
    return data_->solver.get_terminate();
  }

  //--------------------
  std::ostream& MatrixCBSolver::print_termination_code(std::ostream& out) const {
    assert(data_);
    data_->solver.get_terminator()->print_status(out);
    return out;
  }

  //--------------------
  double MatrixCBSolver::get_objval() const {
    assert(data_);

    return data_->solver.get_center_objval();
  }

  //--------------------
  int MatrixCBSolver::get_center(Matrix& y) const {
    assert(data_);

    y = data_->solver.get_center_y();
    return 0;
  }

  //--------------------
  double MatrixCBSolver::get_candidate_value() const {
    assert(data_);

    return data_->solver.get_cand_objval();
  }

  //--------------------
  int MatrixCBSolver::get_candidate(Matrix& y) const {
    assert(data_);

    y = data_->solver.get_cand_y();
    return 0;
  }

  //--------------------
  int MatrixCBSolver::get_approximate_slacks(Matrix& eta) const {
    assert(data_);
    eta.newsize(data_->groundset.get_dim(), 1); chk_set_init(eta, 1);
    assert(data_->solver.get_gs_aggregate().valid());
    Real dummy;
    return data_->solver.get_gs_aggregate().get_minorant(dummy, eta, 0, -1., false);
  }

  //--------------------
  const PrimalData* MatrixCBSolver::get_approximate_primal(const FunctionObject& function) const {
    assert(data_);

    if (data_->fun_model.find(&function) == data_->fun_model.end())
      return 0;
    return data_->fun_model[&function]->get_model()->get_approximate_primal();
  }

  //--------------------
  const PrimalData* MatrixCBSolver::get_center_primal(const FunctionObject& function) const {
    assert(data_);

    if (data_->fun_model.find(&function) == data_->fun_model.end())
      return 0;
    return data_->fun_model[&function]->get_model()->get_center_primal();
  }

  //--------------------
  const PrimalData* MatrixCBSolver::get_candidate_primal(const FunctionObject& function) const {
    assert(data_);

    if (data_->fun_model.find(&function) == data_->fun_model.end())
      return 0;
    return data_->fun_model[&function]->get_model()->get_candidate_primal();
  }

  //--------------------
  double MatrixCBSolver::get_sgnorm() const {
    assert(data_);
    return sqrt(data_->solver.get_aggr_dnormsqr());
  }

  //--------------------
  int MatrixCBSolver::get_subgradient(Matrix& subg) const {
    assert(data_);
    data_->solver.get_aggregate(subg);
    return 0;
  }

  //--------------------
  double MatrixCBSolver::get_cutval() const {
    assert(data_);

    return data_->solver.get_modelval();
  }

  //--------------------
  int MatrixCBSolver::get_function_status(const FunctionObject& function) const {
    assert(data_);
    if (data_->fun_model.find(&function) == data_->fun_model.end()) {
      if (data_->cb_out())
        data_->get_out() << "**** ERROR: MatrixCBSolver::get_function_status(.): cannot find this function object in this problem" << std::endl;
      return 1;
    }
    return data_->fun_model[&function]->get_model()->get_ret_code();
  }


  //--------------------
  int MatrixCBSolver::set_sumbundle(bool use_sumbundle, int n_local_models, const BundleParameters* bp, int strategy) {
    assert(data_);
    if (data_->root) {

      SumBundleParameters sbp;
      if (bp)
        sbp.init(*bp);

      if (use_sumbundle) {
        if (dynamic_cast<SumModel*>(data_->root->get_model())) {
          SumModelParameters smp(this);
          if ((bp) && (smp.init(*bp))) {
            if (cb_out()) {
              get_out() << "**** WARNING: MatrixCBSolver::set_sumbundle(....): initializing SumModelParameters failed" << std::endl;
            }
          }
          smp.set_max_local_models(n_local_models);
          smp.set_update_rule(strategy);
          data_->root->get_model()->set_bundle_parameters(smp);
        }
        sbp.set_acceptable_mode(SumBundle::root);
      } else {
        sbp.set_acceptable_mode(SumBundle::unavailable);
      }
      return data_->root->get_model()->set_sumbundle_parameters(sbp);
    }
    if (cb_out()) {
      get_out() << "**** WARNING: MatrixCBSolver::set_sumbundle(....): currently no bundle model available, so this call has no effect" << std::endl;
    }
    return 0;
  }

  //--------------------
  int MatrixCBSolver::set_max_bundlesize(int mb, const FunctionObject* function) {
    assert(data_);
    SumBlockModel* model = data_->root->get_model();
    if (function) {
      if (data_->fun_model.find(function) == data_->fun_model.end()) {
        if (data_->cb_out())
          data_->get_out() << "**** ERROR: MatrixCBSolver::set_max_bundlesize(..): cannot find this function object in this problem" << std::endl;
        return 1;
      }
      model = data_->fun_model[function]->get_model();
    }

    if (model->get_bundle_parameters()) {
      BundleParameters* bp = model->get_bundle_parameters()->clone_BundleParameters();
      if (bp == 0) {
        if (data_->cb_out())
          data_->get_out() << "**** ERROR: MatrixCBSolver::set_max_bundlesize(..): cloning the bundle parameters of the model to the function object failed" << std::endl;
        return 1;
      }
      int msz = model->get_bundle_parameters()->get_max_model_size();
      bp->set_max_bundle_size(max(mb, msz));
      if (model->set_bundle_parameters(*bp)) {
        if (data_->cb_out())
          data_->get_out() << "**** ERROR: MatrixCBSolver::set_max_bundlesize(..): setting the new bundle parameters failed" << std::endl;
        delete bp;
        return 1;
      }
      delete bp;
    } else {
      if (data_->cb_out())
        data_->get_out() << "**** ERROR: MatrixCBSolver::set_max_bundlesize(..): the model to function object offers no bundle parameters for setting this size" << std::endl;
      return 1;
    }
    return 0;
  }

  //--------------------
  int MatrixCBSolver::set_max_modelsize(int ms, const FunctionObject* function) {
    assert(data_);
    SumBlockModel* model = data_->root->get_model();
    if (function) {
      if (data_->fun_model.find(function) == data_->fun_model.end()) {
        if (data_->cb_out())
          data_->get_out() << "**** ERROR: MatrixCBSolver::set_max_modelsize(..): cannot find this function object in this problem" << std::endl;
        return 1;
      }
      model = data_->fun_model[function]->get_model();
    }
    if (model->get_bundle_parameters()) {
      BundleParameters* bp = model->get_bundle_parameters()->clone_BundleParameters();
      if (bp == 0) {
        if (data_->cb_out())
          data_->get_out() << "**** ERROR: MatrixCBSolver::set_max_modelsize(..): cloning the bundle parameters of the model to the function object failed" << std::endl;
        return 1;
      }
      int bsz = model->get_bundle_parameters()->get_max_bundle_size();
      bp->set_max_model_size(ms);
      bp->set_max_bundle_size(max(ms, bsz));
      if (model->set_bundle_parameters(*bp)) {
        if (data_->cb_out())
          data_->get_out() << "**** ERROR: MatrixCBSolver::set_max_modelsize(..): setting the new bundle parameters failed" << std::endl;
        delete bp;
        return 1;
      }
      delete bp;
    } else {
      if (data_->cb_out())
        data_->get_out() << "**** ERROR: MatrixCBSolver::set_max_modelsize(..): the model to function object offers no bundle parameters for setting this size" << std::endl;
      return 1;
    }
    return 0;
  }

  //--------------------
  int MatrixCBSolver::set_bundle_parameters(const BundleParameters& bp,
    const FunctionObject* function) {
    assert(data_);
    SumBlockModel* model = data_->root->get_model();
    if (function) {
      if (data_->fun_model.find(function) == data_->fun_model.end()) {
        if (data_->cb_out())
          data_->get_out() << "**** ERROR: MatrixCBSolver::set_bundle_parameters(..): cannot find this function object in this problem" << std::endl;
        return 1;
      }
      model = data_->fun_model[function]->get_model();
    }
    return model->set_bundle_parameters(bp);
  }

  //--------------------
  const BundleParameters* MatrixCBSolver::get_bundle_parameters(const FunctionObject* function) const {
    assert(data_);
    SumBlockModel* model = data_->root->get_model();
    if (function) {
      if (data_->fun_model.find(function) == data_->fun_model.end()) {
        if (data_->cb_out())
          data_->get_out() << "**** ERROR: MatrixCBSolver::get_bundle_parameters(..): cannot find this function object in this problem" << std::endl;
        return 0;
      }
      model = data_->fun_model[function]->get_model();
    }
    return model->get_bundle_parameters();
  }

  //--------------------
  int MatrixCBSolver::set_sumbundle_parameters(const SumBundleParametersObject& bp,
    const FunctionObject* function) {
    assert(data_);
    SumBlockModel* model = data_->root->get_model();
    if (function) {
      if (data_->fun_model.find(function) == data_->fun_model.end()) {
        if (data_->cb_out())
          data_->get_out() << "**** ERROR: MatrixCBSolver::set_sumbundle_parameters(..): cannot find this function object in this problem" << std::endl;
        return 1;
      }
      model = data_->fun_model[function]->get_model();
    }
    return model->set_sumbundle_parameters(bp);
  }

  //--------------------
  const BundleData* MatrixCBSolver::get_bundle_data(const FunctionObject* function) const {
    assert(data_);
    SumBlockModel* model = data_->root->get_model();
    if (function) {
      if (data_->fun_model.find(function) == data_->fun_model.end()) {
        if (data_->cb_out())
          data_->get_out() << "**** ERROR: MatrixCBSolver::get_bundle_data(..): cannot find this function object in this problem" << std::endl;
        return 0;
      }
      model = data_->fun_model[function]->get_model();
    }
    return model->get_data();
  }

  //--------------------
  int MatrixCBSolver::reinit_function_model(const FunctionObject* function) {
    assert(data_);
    SumBlockModel* model = data_->root->get_model();
    if (function) {
      if (data_->fun_model.find(function) == data_->fun_model.end()) {
        if (data_->cb_out())
          data_->get_out() << "**** ERROR: MatrixCBSolver::reinit_function_model(..): cannot find this function object in this problem" << std::endl;
        return 1;
      }
      model = data_->fun_model[function]->get_model();
    }
    if (data_->apply_modifications()) {
      if (data_->cb_out()) {
        data_->get_out() << "**** ERROR MatrixCBSolver::reinit_function_model(.): apply_modifications failed" << std::endl;
      }
      return 1;
    }

    model->clear_model();
    return 0;
  }

  //--------------------
  int MatrixCBSolver::clear_aggregates(const FunctionObject* function) {
    assert(data_);
    SumBlockModel* model = data_->root->get_model();
    if (function) {
      if (data_->fun_model.find(function) == data_->fun_model.end()) {
        if (data_->cb_out())
          data_->get_out() << "**** ERROR: MatrixCBSolver::clear_aggregatex(..): cannot find this function object in this problem" << std::endl;
        return 1;
      }
      model = data_->fun_model[function]->get_model();
    }
    if (data_->apply_modifications()) {
      if (data_->cb_out()) {
        data_->get_out() << "**** ERROR MatrixCBSolver::clear_aggregates(.): apply_modifications failed" << std::endl;
      }
      return 1;
    }

    model->clear_aggregates();
    return 0;
  }

  //--------------------
  int MatrixCBSolver::call_primal_extender(const FunctionObject& function, PrimalExtender& primal_extender) {
    assert(data_);
    if (data_->fun_model.find(&function) == data_->fun_model.end()) {
      if (data_->cb_out())
        data_->get_out() << "**** ERROR: MatrixCBSolver::call_primal_extender(.): cannot find this function object in this problem" << std::endl;
      return 1;
    }

    if (data_->apply_modifications()) {
      if (data_->cb_out()) {
        data_->get_out() << "**** ERROR MatrixCBSolver::do_descent_step(.): apply_modifications failed" << std::endl;
      }
      return 1;
    }

    return data_->fun_model[&function]->get_model()->call_primal_extender(primal_extender);
  }

  int MatrixCBSolver::set_term_relprec(const double term_relprec) {
    assert(data_);
    data_->solver.get_terminator()->set_termeps(term_relprec);
    return 0;
  }

  double MatrixCBSolver::get_last_weight() const {
    assert(data_);
    return data_->solver.get_weight();
  }

  double MatrixCBSolver::get_next_weight() const {
    assert(data_);
    return data_->solver.get_bundleweight()->get_weight();
  }

  int MatrixCBSolver::set_next_weight(const double weight) {
    assert(data_);
    if (weight <= 0) return 1;
    data_->solver.get_bundleweight()->set_next_weight(weight);
    return 0;
  }

  int MatrixCBSolver::set_min_weight(const double weight) {
    assert(data_);
    data_->solver.get_bundleweight()->set_minweight(weight);
    return 0;
  }

  int MatrixCBSolver::set_max_weight(const double weight) {
    assert(data_);
    data_->solver.get_bundleweight()->set_maxweight(weight);
    return 0;
  }

  int MatrixCBSolver::set_weight_update(BundleWeight* bw) {
    assert(data_);
    data_->solver.set_bundleweight(bw);
    return 0;
  }

  int MatrixCBSolver::set_variable_metric(int do_var_metric) {
    assert(data_);
    data_->solver.set_variable_metric(do_var_metric);
    return 0;
  }

  int MatrixCBSolver::set_prox(BundleProxObject* bsp) {
    assert(data_);
    data_->solver.set_prox(bsp);
    return 0;
  }

  void MatrixCBSolver::set_active_bounds_fixing(bool allow_fixing) {
    assert(data_);
    data_->solver.set_do_yfixing(allow_fixing);
  }

  void MatrixCBSolver::clear_fail_counts(void) {
    assert(data_);
    data_->solver.clear_fails();
  }

  void MatrixCBSolver::set_eval_limit(Integer eval_limit) {
    assert(data_);
    data_->solver.get_terminator()->set_objevallimit(eval_limit);
  }

  void MatrixCBSolver::set_inner_update_limit(Integer update_limit) {
    assert(data_);
    data_->solver.set_max_updates(update_limit);
  }

  void MatrixCBSolver::set_time_limit(Integer tl) {
    assert(data_);
    if (tl > 0)
      data_->solver.get_terminator()->set_timelimit(&(data_->myclock), Microseconds(tl));
    else
      data_->solver.get_terminator()->set_timelimit(0, Microseconds(0));
  }

  int MatrixCBSolver::set_qp_solver(QPSolverParametersObject* qpparams,
    QPSolverObject* newqpsolver) {
    assert(data_);
    return data_->groundset.set_qpsolver(qpparams, newqpsolver);
  }

  int MatrixCBSolver::set_new_center_point(const Matrix& center_point) {
    assert(data_);
    if (data_->gs_modif->new_vardim() != center_point.dim()) {
      if (data_->cb_out()) {
        data_->get_out() << "**** ERROR MatrixCBSolver::set_new_center_point(.): dimension=" << center_point.dim() << " of input vector does not match current groundset dimension=" << data_->gs_modif->new_vardim() << std::endl;
      }
      return 1;
    }

    if (data_->apply_modifications()) {
      if (data_->cb_out()) {
        data_->get_out() << "**** ERROR MatrixCBSolver::set_new_center_point(.): apply_modifications required before setting new center, but it failed" << std::endl;
      }
      return 1;
    }

    return data_->solver.set_new_center(&center_point);
  }

  int MatrixCBSolver::adjust_multiplier(void) {
    assert(data_);

    if (data_->apply_modifications()) {
      if (data_->cb_out()) {
        data_->get_out() << "**** ERROR MatrixCBSolver::adjust_multiplier(): apply_modifications was required before adjusting (adjusting should be done before modifications), but it failed" << std::endl;
      }
      return 1;
    }

    bool dummy;
    if (data_->root)
      return data_->root->get_model()->adjust_multiplier(dummy);
    return 0;
  }

  int MatrixCBSolver::get_dim() const {
    assert(data_);
    if (data_->apply_modifications()) {
      if (data_->cb_out()) {
        data_->get_out() << "**** ERROR MatrixCBSolver::get_dim(): apply_modifications was required but failed" << std::endl;
      }
      return 1;
    }
    return data_->groundset.get_dim();
  }

  int MatrixCBSolver::get_n_functions() const {
    assert(data_);
    return int(data_->fun_model.size());
  }

  int MatrixCBSolver::get_n_oracle_calls() const {
    assert(data_);
    return data_->solver.get_cntobjeval();
  }

  int MatrixCBSolver::get_n_descent_steps() const {
    assert(data_);
    return data_->solver.get_descent_steps();
  }

  int MatrixCBSolver::get_n_inner_iterations() const {
    assert(data_);
    return data_->solver.get_suminnerit();
  }

  int MatrixCBSolver::get_n_inner_updates() const {
    assert(data_);
    return data_->solver.get_sumupdatecnt();
  }

  bool MatrixCBSolver::get_descent_step() const {
    assert(data_);
    return data_->solver.get_descent_step();
  }

  bool MatrixCBSolver::get_null_step() const {
    assert(data_);
    return data_->solver.get_null_step();
  }

  int  MatrixCBSolver::get_costs(Matrix& costs) const {
    assert(data_);

    if (data_->apply_modifications()) {
      if (data_->cb_out()) {
        data_->get_out() << "**** ERROR MatrixCBSolver::get_costs(.): apply_modifications required before getting new costs, but it failed" << std::endl;
      }
      return 1;
    }
    costs.newsize(data_->groundset.get_dim(), 1); chk_set_init(costs, 1);
    Real dummy;
    return (data_->groundset.get_gs_minorant()).get_minorant(dummy, costs, 0);
  }

  const Matrix* MatrixCBSolver::get_lbounds() const {
    assert(data_);

    if (data_->apply_modifications()) {
      if (data_->cb_out()) {
        data_->get_out() << "**** ERROR MatrixCBSolver::get_lbounds(.): apply_modifications required before getting new lower bounds, but it failed" << std::endl;
      }
    }

    return data_->groundset.get_lby();
  }

  const Matrix* MatrixCBSolver::get_ubounds() const {
    assert(data_);

    if (data_->apply_modifications()) {
      if (data_->cb_out()) {
        data_->get_out() << "**** ERROR MatrixCBSolver::get_ubounds(.): apply_modifications required before getting new upper bounds, but it failed failed" << std::endl;
      }
    }

    return data_->groundset.get_uby();
  }

  const Indexmatrix* MatrixCBSolver::get_fixed_active_bounds() const {
    assert(data_);
    return data_->solver.get_yfixed();
  }

  BundleProxObject* MatrixCBSolver::get_prox() const {
    assert(data_);
    return data_->solver.get_prox();
  }

  bool MatrixCBSolver::pending_oracle_modification(const FunctionObject& function,
    Integer& old_dim,
    Integer& new_dim,
    Integer& append_dim,
    const Indexmatrix*& map_to_old,
    const Indexmatrix*& deleted_indices,
    const Indexmatrix*& new_indices,
    const OracleModification*& oracle_modification
  ) const {
    assert(data_);

    if (data_->fun_model.find(&function) == data_->fun_model.end()) {
      if (data_->cb_out())
        data_->get_out() << "**** Warning: MatrixCBSolver::pending_oracle_modification(........): cannot find this function object in this problem" << std::endl;
      old_dim = new_dim = append_dim = 0;
      map_to_old = deleted_indices = new_indices = 0;
      oracle_modification = 0;
      return false;
    }

    return data_->fun_model[&function]->pending_oracle_modification(old_dim, new_dim, append_dim, map_to_old, deleted_indices, new_indices, oracle_modification);
  }

  void MatrixCBSolver::set_out(std::ostream* o, int pril) {
    assert(data_);
    data_->set_out(o, pril);
  }


  std::ostream& MatrixCBSolver::print_line_summary(std::ostream& out) const {
    assert(data_);
    data_->solver.print_line_summary(out);
    if (data_->root)
      data_->root->get_model()->sbm_transform()->print_statistics(out);
    return out;
  }

  std::ostream& MatrixCBSolver::print_statistics(std::ostream& out) const {
    assert(data_);
    data_->solver.print_statistics(out);
    if (data_->root)
      data_->root->get_model()->sbm_transform()->print_statistics(out);
    return out;
  }

  const BundleSolver* MatrixCBSolver::get_solver(void) const {
    assert(data_);
    return &data_->solver;
  }


} //end namespace ConicBundle

