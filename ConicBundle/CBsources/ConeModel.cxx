/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/ConeModel.cxx
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



#include "mymath.hxx"
#include "ConeModel.hxx"
#include "BundleIdProx.hxx"

using namespace CH_Tools;
using namespace CH_Matrix_Classes;

namespace ConicBundle {

  // *****************************************************************************
  //                            clear()
  // *****************************************************************************

  void ConeModel::clear() {
    if (bundlehandler) {
      if (bundlehandler->remove_contributions()) {
        if (cb_out())
          get_out() << "**** ERROR: ConeModel::clear(...): bundlehandler->remove_contributions failed" << std::endl;
      }
    }

    SumBlockModel::clear();
  }

  // *****************************************************************************
  //                            ConeModel
  // *****************************************************************************

  ConeModel::ConeModel(CBout* cb, int incr) :
    CBout(cb, incr), SumBlockModel(cb, incr) {
  }

  // *****************************************************************************
  //                            ConeModel
  // *****************************************************************************

  ConeModel::~ConeModel() {
  }

  // *****************************************************************************
  //                            eval_function
  // *****************************************************************************


  int ConeModel::eval_function(Integer& ub_fid,
    Real& ub,
    Integer y_id,
    const Matrix& y,
    Real nullstep_bound,
    Real relprec) {
    assert(relprec > 0.);

    if (get_data()->dim < 0)
      get_data()->dim = y.dim();
    assert(get_data()->dim == y.dim());

    //delete old info
    ub_fid = get_data()->cand_ub_mid = -1;
    get_data()->cand_relprec = relprec;
    get_data()->cand_id = y_id;

    if (bundlehandler) {
      bundlehandler->clear_cand_minorants();
      get_data()->max_old_minorants = max(get_data()->max_old_minorants, sumbundle_parameters->get_max_bundle_size());
    }

    int err = evaluate_oracle(y_id, y, nullstep_bound, relprec);
    if ((err) && (cb_out())) {
      get_out() << "**** ERROR: ConeModel::eval_function(): evaluate_oracle returned " << err << std::endl;
      err++;
    }

    if (get_data()->cand_ub_mid == get_data()->modification_id) {
      if ((get_data()->function_task != ObjectiveFunction) && (get_data()->cand_ub < 0))
        ub = 0.;
      else
        ub = get_data()->function_factor * get_data()->cand_ub;
      ub_fid = get_data()->cand_ub_mid;

      if (cb_out(1)) {
        get_out().precision(8);
        get_out() << " cand_ub=" << get_data()->cand_ub << " ub=" << ub;
        if (cb_out(3)) {
          Real lb = provide_cand_minorant().evaluate(y_id, y) * get_data()->function_factor;
          get_out() << " lb=" << lb;
        }
        get_out() << std::endl;
      }

      //update old minorants if necessary
      get_data()->max_old_minorants = max(get_data()->max_old_minorants, get_bundle_parameters()->get_max_bundle_size());
      if (get_data()->max_old_minorants > 0) {
        MinorantPointer& cand_minorant = provide_cand_minorant();
        if ((!cand_minorant.valid()) || (get_data()->store_old_minorant(cand_minorant))) {
          if (cb_out())
            get_out() << "**** WARNING ConeModel::eval_function(): get_data()->store_old_minorant() failed" << std::endl;
        }
      }

      //update contribution to the sumbundel
      if ((bundlehandler) && (get_data()->sumbundle.has_bundle_for(get_data()->function_task))) {
        MinorantPointer& cand_minorant = provide_cand_minorant();
        if (cand_minorant.valid()) {
          bundlehandler->set_cand_minorant(get_data()->function_task, cand_minorant);
          if (get_data()->sumbundle.has_contributions()) {
            if (bundlehandler->contribute_new_minorants()) {
              if (cb_out())
                get_out() << "**** ERROR ConeModel::eval_function(): bundlehandler->contribute_new_minorants failed" << std::endl;
              err++;
            }
          }
        }
      }
    }

    return err;
  }

  // //==================================== cut begin
  // //---- evaluate  	
  // CH_Tools::Microseconds start_eval=clock.time();
  // int err=0;
  // Real cmpval=((get_data()->cand_ub_mid>=0)?0.1*CH_Matrix_Classes::abs(get_data()->cand_ub):1.);

  // Real eigval_bound=nullstep_bound/get_data()->function_factor;
  // //Real eigval_relprec=0.1*relprec/get_data()->function_factor*(fabs(nullstep_bound)+1.)/(fabs(eigval_bound)+1.);
  // //Real eigval_relprec=0.1*relprec*(cmpval+1./get_data()->function_factor)/(cmpval+1.);
  // Real eigval_relprec=max(10*eps_Real,0.5*relprec/min(1.,(1.-relprec))*(cmpval+1./get_data()->function_factor)/(cmpval+1.));
  // nr_eval++;

  // Matrix tmpmat;

  // if (get_data()->topvecs.coldim()>get_data()->bundlevecs.coldim()){
  //   if (get_data()->activedim<=0)
  //     tmpmat=get_data()->topvecs;
  //   else
  //     tmpmat.init(get_data()->topvecs.rowdim(),min(get_data()->topvecs.coldim(),get_data()->activedim),get_data()->topvecs.get_store());
  // }
  // else {
  //   if (get_data()->activedim<=0)
  //     tmpmat=get_data()->bundlevecs;
  //   else
  //     tmpmat.init(get_data()->bundlevecs.rowdim(),min(get_data()->bundlevecs.coldim(),get_data()->activedim),get_data()->bundlevecs.get_store());
  // }

  // //if (get_data()->bundlevecs.coldim()>0)
  // //  tmpmat=get_data()->bundlevecs;
  // //if (get_data()->topvecs.coldim()>0)
  // //  tmpmat=get_data()->topvecs;
  // ConePrimalExtender* pep=0;	

  // preeval_time+=clock.time()-start_eval;  
  // start_eval=clock.time();

  // ret_code = oracle->evaluate(y,tmpmat,eigval_relprec,eigval_bound,
  // 			      get_data()->cand_Ritzvec,get_data()->cand_Ritzval,pep);
  // eval_time+=clock.time()-start_eval;
  // start_eval=clock.time();

  // if (get_data()->cand_Ritzval.dim()==0) {
  //   if (cb_out()) 
  //     get_out()<<"**** ERROR: ConeModel::eval_function(): function returned no Ritz_values"<<std::endl;
  //   err++;
  // }

  // if (get_data()->cand_Ritzval.dim()!=get_data()->cand_Ritzvec.coldim()){
  //   if (cb_out()) 
  //     get_out()<<"**** ERROR: ConeModel::eval_function(): function returned different number of Ritz_values and Ritz_vectors"<<std::endl;
  //   err++;
  // }

  // //--- if needed, apply the modifications to the old data
  // if (pep){
  //   if (call_primal_extender(*pep)){
  //     if (cb_out()) 
  // 	get_out()<<"**** WARNING: ConeModel::eval_function(): call_primal_extender failed"<<std::endl;
  //     err++;
  //   }
  //   delete pep;
  // }

  // if ((get_data()->bundlevecs.coldim()>0)&&(get_data()->bundlevecs.rowdim()!=get_data()->cand_Ritzvec.rowdim())){
  //   if (cb_out()) 
  //     get_out()<<"**** ERROR: ConeModel::eval_function(): order of the matrix changed without having called ConeModel::clear_model first"<<std::endl;
  //   err++;
  // }

  // if (err){
  //   posteval_time+=clock.time()-start_eval;
  //   return err;
  // }


  // //--- update the most important subspace (thereby maybe correcting the function value)
  // Real lmax=max(get_data()->cand_Ritzval);

  // //compute a basis of the subspace spanned by topvecs/bundlevecs and cand_Ritzvec into cand_Ritzvec
  // if (get_data()->topvecs.coldim()>0)
  //   tmpmat=get_data()->topvecs;
  // else
  //   tmpmat=get_data()->bundlevecs;
  // tmpmat.concat_right(get_data()->cand_Ritzvec);
  // Indexmatrix piv;
  // Integer r1=tmpmat.QR_factor(piv);
  // tmpmat.delete_cols(Range(r1,tmpmat.coldim()-1));
  // piv.init(Range(0,r1-1));
  // get_data()->topvecs.init(tmpmat.rowdim(),tmpmat.coldim(),0.);
  // for(Integer i=0;i<get_data()->topvecs.coldim();i++) 
  //   get_data()->topvecs(i,i)=1.;
  // tmpmat.Q_times(get_data()->topvecs,r1);

  // //compute the eigenvalue decomposition of the projected matrix
  // Matrix Ritz_vecs;
  // int ret_val=oracle->evaluate_projection(y,get_data()->topvecs,1e-10,Ritz_vecs,get_data()->Ritz_values);
  // if (ret_val){
  //   if (cb_out()){
  //     get_out()<<"**** WARNING: ConeModel::eval_function(): evaluate_projection returned code = "<<ret_val<<std::endl;
  //   }
  // }      
  // genmult(get_data()->topvecs,Ritz_vecs,tmpmat);
  // swap(get_data()->topvecs,tmpmat);

  // if (cb_out(1))
  //   get_out()<<"  Cone:";

  // if (cb_out(2))
  //   get_out()<<" lmax="<<lmax<<" Ritz_values="<<transpose(get_data()->Ritz_values);

  // lmax=get_data()->Ritz_values(0);

  // ub_fid=get_data()->cand_ub_mid=get_data()->modification_id;
  // get_data()->cand_ub=lmax+eigval_relprec*(fabs(lmax)+1.);



  // if (ret_code){
  //   if (cb_out()) get_out()<<"**** WARNING: ConeModel::eval_function(): function returned code = "<<ret_code<<std::endl;
  //   posteval_time+=clock.time()-start_eval;
  //   return ret_code;
  // }

  // if (ub>=CB_plus_infinity){
  //   if (cb_out()) get_out()<<"**** WARNING: ConeModel::eval_function(): >=CB_plus_infinity"<<std::endl;
  // }    

  // posteval_time+=clock.time()-start_eval;

  // //==================================== cut end




// *****************************************************************************
//                              eval_model
// *****************************************************************************

//evaluate the current cutting model in $y$ 

  int ConeModel::eval_model(Real& lb,
    Integer y_id,
    const Matrix& y,
    Real relprec) {
    CH_Tools::Microseconds evalmodel_start = clock.time();


    int err = 0;
    if ((((bundlehandler == 0) || (!get_data()->sumbundle.active())) && (!get_data()->get_local_model_active())) ||
      ((bundlehandler != 0) && (bundlehandler->initialization_needed()))) {
      if (candidate_available()) {
        BundleIdProx H;
        int retval = update_model(new_subgradient, y_id, y, y_id, y, max(1e-6, 1e-3 * std::fabs(get_data()->cand_ub)), H);
        if (retval) {
          if (cb_out())
            get_out() << "**** ERROR ConeModel::eval_model: update_model failed and returned " << retval << std::endl;
          err++;
        }
      } else
        err++;
    }

    lb = min_Real;
    if (err == 0) {
      if ((bundlehandler) && (get_data()->sumbundle.active())) {
        if (get_data()->sumbundle.has_working_roots()) {
          if (bundlehandler->eval_model(lb, y_id, y)) {
            if (cb_out())
              get_out() << "**** ERROR ConeModel::eval_model: get_data()->sumbundle.eval_model failed" << std::endl;
            err++;
          }
        }
      } else {
        if (evaluate_local_model(lb, y_id, y, relprec)) {
          if (cb_out())
            get_out() << "**** ERROR ConeModel::eval_model: get_data()->sumbundle.evaluate_local_model failed" << std::endl;
          err++;
        } else {
          if ((get_data()->function_task == ObjectiveFunction) || (lb > 0))
            lb *= get_data()->function_factor;
          else
            lb = 0.;
        }
      }
    }

    if (cb_out(3)) {
      get_out().precision(12);
      get_out() << " modval=" << lb;
    }

    evalmodel_time += clock.time() - evalmodel_start;

    return err;
  }


  // *****************************************************************************
  //                                update_model
  // *****************************************************************************

  int ConeModel::update_model(ModelUpdate model_update,
    Integer center_id,
    const Matrix& center_y,
    Integer y_id,
    const Matrix& y,
    Real model_maxviol,
    BundleProxObject& H) {
    CH_Tools::Microseconds update_start = clock.time();
    int err = 0;

    //call first the local part, it hast to do at least the step
    if (update_local_model(model_update,
      center_id, center_y,
      y_id, y,
      model_maxviol,
      H)
      ) {
      if (cb_out())
        get_out() << "**** ERROR ConeModel::update_model: update_local_model failed" << std::endl;
      err++;
    }
    get_data()->local_model_active = true;

    if (bundlehandler) {
      if (get_data()->sumbundle.active())
        get_data()->local_model_active = false;

      //rearrange the existing bundle and include the new element
      if (bundlehandler->update_model(model_update, center_id, center_y, y_id, y, model_maxviol, H)) {
        if (cb_out())
          get_out() << "**** ERROR ConeModel::update_model: bundlehandler->update_model failed" << std::endl;
        err++;
      }
    }

    updatemodel_time += clock.time() - update_start;

    //topvecs holds the Ritz-vectors sorted nonincreasingly by Ritz-values

    return err;
  }


  // *****************************************************************************
  //                                update_model
  // *****************************************************************************


  int ConeModel::update_model(ModelUpdate model_update,
    Integer center_id,
    const Matrix& center_y,
    Integer y_id,
    const Matrix& y,
    Real model_maxviol,
    BundleProxObject& H,
    Real& model_deviation,
    Real& model_curvature) {

    CH_Tools::Microseconds update_start = clock.time();
    assert(get_data());
    assert((y_id == get_data()->center_id) || (y_id == get_data()->cand_id));

    if ((!get_data()->aggregate.valid()) || ((y_id != get_data()->center_id) && (y_id != get_data()->cand_id))) {
      model_deviation = 0.;
    } else {
      if (get_data()->cand_id == y_id) {
        model_deviation = get_data()->cand_ub * get_data()->function_factor;
      } else {
        model_deviation = get_data()->center_ub * get_data()->function_factor;
      }
      model_deviation -= get_data()->aggregate.evaluate(y_id, y);
    }

    if ((bundlehandler) && (model_curvature >= 0.)) {
      Integer n = min(Integer(get_data()->old_minorants.size()), Integer(10));
      //Integer n= Integer(get_data()->old_minorants.size());
      Indexmatrix ind(n + 1, 1);
      ind.init(0, 1, Integer(0));
      Integer indi = get_data()->minorant_nexti;
      for (Integer i = 0; i < n; i++) {
        if (--indi < 0)
          indi = Integer(get_data()->old_minorants.size()) - 1;
        ind.concat_below(indi);
      }
      if (get_data()->aggregate.valid()) {
        ind.concat_below(Integer(get_data()->old_minorants.size()));
        get_data()->old_minorants.push_back(MinorantPointer(get_data()->aggregate, 1 / get_data()->function_factor));
      }
      model_curvature = bundlehandler->guess_curvature(get_data()->old_minorants, ind, y_id, y, model_maxviol / get_data()->function_factor);
      if (get_data()->aggregate.valid())
        get_data()->old_minorants.pop_back();
    }

    updatemodel_time += clock.time() - update_start;

    return update_model(model_update, center_id, center_y, y_id, y, model_maxviol, H);
  }



  // *****************************************************************************
  //                                recompute_center
  // *****************************************************************************

  int ConeModel::recompute_center(Integer& new_center_fid,
    Real& new_ub,
    Integer center_id,
    const Matrix& center_y,
    bool accept_only_higher,
    Real relprec) {
    if (get_data()->dim < 0)
      get_data()->dim = center_y.dim();
    assert(get_data()->dim == center_y.dim());

    if ((relprec <= 0) && ((center_id != get_data()->center_id) || (get_data()->center_relprec <= 0))) {
      Real rp = 1e-3;
      if ((get_data()->center_relprec > 0) && (rp > get_data()->center_relprec))
        rp = get_data()->center_relprec;
      if (cb_out()) {
        get_out() << "\n**** WARNING: ConeModel::recompute_center: called with relprec=" << relprec << "<=0 but no previous precision available for center_id=" << center_id << ", using default choice " << rp << std::endl;
      }
      relprec = rp;
    }
    if (relprec <= 0)
      relprec = get_data()->center_relprec;

    if ((accept_only_higher) &&
      ((get_data()->center_id != center_id) || (get_data()->center_ub_mid != get_data()->modification_id))
      ) {
      if (cb_out()) {
        get_out() << "\n**** WARNING: ConeModel::recompute_center: called with accept_only_higher=true but there is no valid previous value or get_data()->center_id=" << get_data()->center_id << " differs from center_id=" << center_id << "; this should never happen, setting accept_only_higher to false" << std::endl;
      }
      accept_only_higher = false;
    }

    if ((get_data()->center_id == center_id) &&
      (get_data()->center_ub_mid == get_data()->modification_id) &&
      (relprec >= get_data()->center_relprec)
      ) {
      new_center_fid = get_data()->center_ub_mid;
      new_ub = 0;

      if ((get_data()->function_task == ObjectiveFunction) || (get_data()->center_ub > 0.))
        new_ub = get_data()->center_ub * get_data()->function_factor;

      return 0;
    }

    int status = eval_function(new_center_fid, new_ub, center_id, center_y, CB_plus_infinity, relprec);
    if (status > 0) {
      if (cb_out()) {
        get_out() << "\n**** ERROR: ConeModel::recompute_center: eval_function failed and retured " << status << std::endl;
      }
      return status;
    }
    if ((!accept_only_higher) || (new_ub > get_data()->center_ub)) {
      status = get_data()->do_step(center_id);
      if (status) {
        if (cb_out()) {
          get_out() << "\n**** ERROR: ConeModel::recompute_center: get_data()->do_step failed and retured " << status << std::endl;
        }
        return status;
      }
    } else {
      get_data()->center_ub_mid = get_data()->modification_id;
      new_ub = 0;
      if ((get_data()->function_task == ObjectiveFunction) || (get_data()->center_ub > 0.))
        new_ub = get_data()->center_ub * get_data()->function_factor;

      assert(new_center_fid == get_data()->modification_id);
    }

    return status;
  }


  // *****************************************************************************
  //                        provide_model_aggregate
  // *****************************************************************************


  int ConeModel::provide_model_aggregate(Integer y_id, const Matrix& iny) {
    if (get_data()->aggregate.valid()) {
      return 0;
    }

    int err = 0;
    //--- first ensure that a bundle exists, initializing it if possible
    if ((((bundlehandler == 0) || (!get_data()->sumbundle.active())) && (!get_data()->get_local_model_active())) ||
      ((bundlehandler != 0) && (bundlehandler->initialization_needed()))) {
      BundleIdProx H;
      int retval = update_model(new_subgradient, y_id, iny, y_id, iny, max(1e-6, std::fabs(get_data()->cand_ub)), H);
      if (retval) {
        if (cb_out())
          get_out() << "**** ERROR ConeModel::provide_model_aggregate: update_model failed and returned " << retval << std::endl;
        err++;
      }
    }

    if (err == 0) {
      //--- form an aggregate 
      get_data()->local_aggregate.clear();
      get_data()->aggregate.clear();


      if ((bundlehandler) && (get_data()->sumbundle.active())) {
        //----  sumbundle is active, get it there
        if (bundlehandler->provide_model_aggregate()) {
          if (cb_out())
            get_out() << "**** ERROR ConeModel::provide_model_aggregate: bundlehandler->provide_model_aggregate() failed" << std::endl;
          err++;
        } else if (get_data()->sumbundle.get_mode(get_data()->function_task) == SumBundle::root) {
          if (get_data()->sumbundle.get_local_model_aggregate(get_data()->local_aggregate)) {
            if (cb_out())
              get_out() << "**** ERROR ConeModel::provide_model_aggregate: get_local_model_aggregate failed" << std::endl;
            err++;
          } else
            get_data()->aggregate = get_data()->local_aggregate;
        } else {
          if (get_data()->sumbundle.get_contributed_model_aggregate(get_data()->aggregate)) {
            if (cb_out())
              get_out() << "**** ERROR ConeModel::provide_model_aggregate: get_contributed_model_aggregate failed" << std::endl;
            err++;
          }
        }
      } else {
        //----  local model is active, get it there
        assert(get_data()->get_local_model_active());
        if (provide_local_aggregate(y_id, iny)) {
          if (cb_out())
            get_out() << "**** ERROR ConeModel::provide_model_aggregate: get_contributed_model_aggregate failed" << std::endl;
          err++;
        } else
          get_data()->aggregate = get_data()->local_aggregate;
      }

      if (get_data()->aggregate.valid()) {
        get_data()->aggregate_id++;
      }
    }

    return err;
  }

  // *****************************************************************************
  //                         check_center_validity_by_candidate
  // *****************************************************************************


  int ConeModel::check_center_validity_by_candidate(bool& cand_minorant_is_below,
    Integer center_id,
    const Matrix& center_y) {
    int err = 0;
    if ((get_data()->center_id != center_id) ||
      (get_data()->center_ub_mid != get_data()->modification_id)
      ) {
      if (cb_out())
        get_out() << "**** WARNING: ConeModel::check_center_validity_by_candidate: no center available for checking validity" << std::endl;
      err++;
    }

    MinorantPointer& cand_minorant = provide_cand_minorant();
    if ((err) || (!cand_minorant.valid())) {
      if (cb_out())
        get_out() << "**** WARNING: ConeModel::check_center_validity_by_candidate: providing candidate minorant failed" << std::endl;
      err++;
    } else {
      Real lin_approx = get_data()->function_factor * cand_minorant.evaluate(center_id, center_y);
      Real center_val = 0.;
      if ((get_data()->function_task == ObjectiveFunction) || (get_data()->center_ub > 0.))
        center_val = get_data()->center_ub * get_data()->function_factor;

      if (center_val - lin_approx < -max(eps_Real, min(get_data()->cand_relprec, get_data()->center_relprec)) * (fabs(center_val) + 1.)) {  //minorant cuts above!
        cand_minorant_is_below = false;
        if (cb_out()) {
          get_out().precision(12);
          get_out() << "**** WARNING: ConeModel::check_center_validity_by_candidate: center_val=" << center_val << " is below lin_approx=" << lin_approx << std::endl;
        }
      } else {
        cand_minorant_is_below = true;
      }
    }

    return err;
  }

  // *****************************************************************************
  //                              sumbundle_mode
  // *****************************************************************************

  int ConeModel::sumbundle_mode(SumBundle::Mode& mode,
    SumBundleHandler* bh,
    AffineFunctionTransformation* aft) {
    assert((mode != SumBundle::child) || (bh != 0));

    SumBundle::Mode next_mode = sumbundle_parameters->get_acceptable_mode();
    if ((next_mode == SumBundle::root) ||
      ((next_mode == SumBundle::inactive) && (mode == SumBundle::root))
      ) {
      next_mode = SumBundle::root;
    } else if ((next_mode == SumBundle::unavailable) ||
      ((bh == 0) && ((mode != SumBundle::root) || (next_mode == SumBundle::child)))
      ) {
      next_mode = SumBundle::unavailable;
    } else {   //note: bh!=0 (all cases==0 are covered)
      assert(bh);
      //next_mode is either inactive or child
      if ((next_mode == SumBundle::child) && (mode == SumBundle::root)) {
        next_mode = SumBundle::inactive;
      } else if ((mode != SumBundle::unavailable) && (!bh->handles(get_data()->function_task))) {
        next_mode = SumBundle::unavailable;
      } else { //mode = inactive or child or unavailable
        next_mode = mode;
      }
    }

    mode = next_mode;

    SumBundle::Mode current_mode = SumBundle::unavailable;
    if (get_data()->sumbundle.has_bundle_for(get_data()->function_task)) {
      current_mode = get_data()->sumbundle.get_mode(get_data()->function_task);
    }
    assert(!(((current_mode == SumBundle::root) || (current_mode == SumBundle::child)) && (get_data()->get_local_model_active()))); //note: initially both might be inactive

    //--- if no longer a child, remove contributions
    if ((current_mode == SumBundle::child) && (mode != SumBundle::child)) {
      assert(bundlehandler);
      assert(bundlehandler->handles(get_data()->function_task));
      bundlehandler->remove_contributions();
    }

    //switch to local model if necessary
    int err = 0;
    if (((current_mode == SumBundle::root) || (current_mode == SumBundle::child)) &&
      ((mode == SumBundle::unavailable) || (mode == SumBundle::inactive))) {
      assert(bundlehandler);
      assert(bundlehandler->handles(get_data()->function_task));
      assert(get_data()->sumbundle.get_n_contributors(get_data()->function_task) > 0);
      if (switch_to_local(get_data()->sumbundle.get_bundle(get_data()->function_task),
        get_data()->sumbundle.get_coeff(get_data()->function_task))) {
        if (cb_out()) {
          get_out() << "*** WARNING: ConeModel::sumbundle_mode(): switch_to_local() failed" << std::endl;
        }
        err++;
      }
      get_data()->local_model_active = true;
    }

    //initialize sumbundle if necessary
    MinorantPointer init_aggr;
    Real init_coeff;
    if (((mode == SumBundle::root) || (mode == SumBundle::child)) &&
      ((current_mode == SumBundle::unavailable) || (current_mode == SumBundle::inactive))) {
      if (get_aggregate_and_coefficient(init_aggr, init_coeff)) {
        if (cb_out()) {
          get_out() << "*** WARNING: ConeModel::sumbundle_mode(): initialize_sumbundle() failed" << std::endl;
        }
        err++;
      }
    }

    //carry out the decision:
    //--- if mode==SumBundle::unavailable or no parent handler, remove an existing bundlehandler;
    //    if mode==SumBundle::inactive with parent handler, set it to inactive but keep following 
    if ((mode == SumBundle::unavailable) || (mode == SumBundle::inactive)) {
      get_data()->local_model_active = true;
      if (mode == SumBundle::unavailable) {
        bh = 0;
      }
      if (bundlehandler != 0) {
        if (bundlehandler->set_parent_information(bh, aft, SumBundle::inactive)) {
          if (cb_out()) {
            get_out() << "*** ERROR: ConeModel::sumbundle_mode(): SumBundleHandler::set_parent_information() failed" << std::endl;
          }
          err++;
        }
        if (bh == 0) {
          bundlehandler->clear_model();
          delete bundlehandler;
          bundlehandler = 0;
          mode = SumBundle::unavailable;
        }
      }
    } else {
      //--- mode is either root or child
      get_data()->local_model_active = false;

      //--- first make sure a bundlehandler exists 
      if (bundlehandler == 0) {
        assert(get_data()->sumbundle.active() == false);
        std::map<FunctionTask, Real> factor_map;
        factor_map[get_data()->function_task] = get_data()->function_factor;
        bundlehandler = new SumBundleHandler(get_data()->sumbundle, factor_map, sumbundle_parameters);
        bundlehandler->set_cbout(this, 0);
      }

      //--- make sure it has a bundle
      if ((bundlehandler->handles(get_data()->function_task)) &&
        (!get_data()->sumbundle.active())
        ) {
        if (bundlehandler->reset_function_factor(get_data()->function_task, get_data()->function_factor)) {
          if (cb_out()) {
            get_out() << "*** ERROR: ConeModel::sumbundle_mode(): bundlehandler->reset_function_factor failed" << std::endl;
          }
          err++;
        }
        if (get_data()->sumbundle.get_n_contributors(get_data()->function_task) == 0) {
          MinorantBundle new_model;
          new_model.push_back(init_aggr);
          Matrix new_coeff(1, 1, init_coeff);
          if (bundlehandler->contribute_initial_bundle(get_data()->function_task,
            new_model,
            new_coeff)) {
            if (cb_out()) {
              get_out() << "*** ERROR: ConeModel::sumbundle_mode(): bundlehandler->contribute_initial_bundle() failed" << std::endl;
            }
            err++;
          }
          if (bundlehandler->normalize_sumbundle()) {
            if (cb_out()) {
              get_out() << "*** ERROR: ConeModel::sumbundle_mode(): bundlehandler->normalize_sumbundle() failed" << std::endl;
            }
            err++;
          }
        } else {
          // the sumbundle has been updated along, only the aggregate needs to be installed
          if ((!init_aggr.valid()) ||
            (bundlehandler->install_external_aggregate(get_data()->function_task,
              init_aggr,
              init_coeff))
            ) {
            if (cb_out()) {
              get_out() << "*** ERROR: ConeModel::sumbundle_mode(): bundlehandler->install_external_aggregate() failed" << std::endl;
            }
            err++;
          }
        }
      }

      //--- now it has a bundle (maybe an empty one), set the parent and mode information
      if (bundlehandler->set_parent_information(bh, aft, mode)) {
        if (cb_out()) {
          get_out() << "*** ERROR: ConeModel::sumbundle_mode(): SumBundleHandler::set_parent_information() failed" << std::endl;
        }
        err++;
      }

      //--- if it is a new root, set the local sumbundle_parameters
      if ((mode == SumBundle::root) && (current_mode == SumBundle::child)) {
        if (bundlehandler->set_bundle_parameters(*sumbundle_parameters)) {
          if (cb_out()) {
            get_out() << "*** ERROR: ConeModel::sumbundle_mode(): set_bundle_parameters() failed" << std::endl;
          }
          err++;
        }
      }

      //--- if it is a child, contribute to the parent
      if (mode == SumBundle::child) {
        //add the contribution to the parent if it exists
        if (bundlehandler->add_contributions()) {
          if (cb_out()) {
            get_out() << "*** ERROR: ConeModel::sumbundle_mode(): add_contributions() failed" << std::endl;
          }
          err++;
        }
      }
    }

    return err;
  }


  // *****************************************************************************
  //                                start_sumaugmodel
  // *****************************************************************************


  int ConeModel::start_sumaugmodel(QPModelDataPointer& qpblockp,
    CH_Matrix_Classes::Integer cand_id,
    const CH_Matrix_Classes::Matrix& cand_y,
    const CH_Matrix_Classes::Indexmatrix* indices,
    SumBundleHandler* bh,
    SumBundle::Mode mode,
    AffineFunctionTransformation* aft) {
    //first check what to do in terms of the sumbundle
    int err = 0;
    if (sumbundle_mode(mode, bh, aft)) {
      if (cb_out()) {
        get_out() << "**** ERROR ConeModel::start_augmodel(): sumbundle_mode() failed" << std::endl;
      }
      err++;
    } else {
      if (mode == SumBundle::child) {
        qpblockp.clear_model_data_ptr();
      } else if (mode == SumBundle::root) {
        //use the local summodel
        assert(bundlehandler);
        assert(get_data()->sumbundle.active());
        if (bundlehandler->start_augmodel(qpblockp, cand_id, cand_y, indices, get_data()->function_task)) {
          if (cb_out()) {
            get_out() << "*** ERROR: ConeModel::start_augmodel(): bundlehandler->start_augmodel() failed" << std::endl;
          }
          err++;
        } else {
          get_data()->local_aggregate.clear();
          get_data()->aggregate.clear();
        }
      } else {
        if (start_local_augmodel(qpblockp, cand_id, cand_y, indices)) {
          if (cb_out()) {
            get_out() << "*** ERROR: ConeModel::start_augmodel(): bundlehandler->start_augmodel() failed" << std::endl;
          }
          err++;
        } else {
          get_data()->local_aggregate.clear();
          get_data()->aggregate.clear();
        }
      }
    }

    return err;
  }


  // *****************************************************************************
  //                                make_model_aggregate
  // *****************************************************************************

  int ConeModel::make_model_aggregate(bool& penalty_changed,
    bool keep_penalty_fixed) {

    penalty_changed = false;
    Real old_function_factor = get_data()->function_factor;
    Real new_function_factor = old_function_factor;
    get_data()->local_aggregate.clear();
    get_data()->aggregate.clear();

    int err = 0;
    if (bundlehandler) {
      //even if not active this is needed to maybe update the coefficients
      bool local_penalty_changed = false;
      if (bundlehandler->make_model_aggregate(local_penalty_changed, keep_penalty_fixed)) {
        if (cb_out(0)) {
          get_out() << "\n**** ERROR:  ConeModel::make_model_aggregate(...): make_model_aggregate failed for bundlehandler" << std::endl;
        }
        err++;
      } else if ((get_data()->function_task == AdaptivePenaltyFunction) &&
        (get_data()->sumbundle.active()) &&
        (local_penalty_changed)) {
        penalty_changed = true;
        new_function_factor *= bundlehandler->get_increase_factor();
      } else if (get_data()->sumbundle.get_mode(get_data()->function_task) == SumBundle::root) {
        if (get_data()->sumbundle.get_local_model_aggregate(get_data()->local_aggregate)) {
          if (cb_out())
            get_out() << "**** ERROR ConeModel::make_model_aggregate(): get_local_model_aggregate failed" << std::endl;
          err++;
        } else {
          get_data()->aggregate = get_data()->local_aggregate;
          get_data()->aggregate_id++;
        }
      } else if (get_data()->sumbundle.get_mode(get_data()->function_task) == SumBundle::child) {
        if (get_data()->sumbundle.get_contributed_model_aggregate(get_data()->aggregate)) {
          if (cb_out())
            get_out() << "**** ERROR ConeModel::make_model_aggregate(): get_contributed_model_aggregate failed" << std::endl;
          err++;
          assert(!get_data()->aggregate.valid());
        } else {
          get_data()->aggregate_id++;
        }
      }
    }

    if (get_data()->get_local_model_active()) {
      if (make_local_aggregate(new_function_factor, keep_penalty_fixed)) {
        if (cb_out()) {
          get_out() << "**** ERROR ConeModel::make_model_aggregate(): make_local_aggregate() failed" << std::endl;
        }
        err++;
      } else {
        if (std::fabs(new_function_factor - old_function_factor) < 1e-10 * old_function_factor) {
          get_data()->aggregate = get_data()->local_aggregate;
          get_data()->aggregate_id++;
        }
      }
    }

    if (std::fabs(new_function_factor - old_function_factor) > 1e-12 * old_function_factor) {
      penalty_changed = true;
      if (new_function_factor <= old_function_factor) {
        if (cb_out()) {
          get_out() << "**** ERROR ConeModel::make_model_aggregate(): the penalty is reported to have changed but the new function_factor=" << new_function_factor << " is not bigger than the old value =" << old_function_factor << std::endl;
        }
        err++;
      }
      assert(new_function_factor > (1 - 1e-12) * old_function_factor);
      get_data()->function_factor = new_function_factor;
      if (get_data()->cand_ub_mid == get_data()->modification_id)
        get_data()->cand_ub_mid++;
      if (get_data()->center_ub_mid == get_data()->modification_id)
        get_data()->center_ub_mid++;
      get_data()->modification_id++;
      if (bundlehandler) {
        if (bundlehandler->reset_function_factor(get_data()->function_task, get_data()->function_factor)) {
          if (cb_out()) {
            get_out() << "**** ERROR ConeModel::adjust_multiplier(): bundlehandler->reset_function_factor failed" << std::endl;
          }
          err++;
        }
      }
    }

    return err;
  }

  // *****************************************************************************
  //                                 lb_function
  // *****************************************************************************

  //returns a *quick* lower bound for the function value at y
  //(eg by a previous subgradient)

  Real ConeModel::lb_function(Integer y_id, const Matrix& y) {
    assert((get_data()->dim < 0) || (get_data()->dim == y.dim()));
    if (get_data()->dim < 0)
      get_data()->dim = y.dim();

    Real lb = CB_minus_infinity;

    if (get_data()->aggregate.valid()) {
      lb = get_data()->aggregate.evaluate(y_id, y);
    } else if (get_data()->center_ub_mid == get_data()->modification_id) {
      MinorantPointer& mnrt = provide_center_minorant();
      if (mnrt.valid()) {
        lb = get_data()->function_factor * mnrt.evaluate(y_id, y);
      }
    } else if (get_data()->cand_ub_mid == get_data()->modification_id) {
      MinorantPointer& mnrt = provide_cand_minorant();
      if (mnrt.valid()) {
        lb = get_data()->function_factor * mnrt.evaluate(y_id, y);
      }
    }

    if (get_data()->function_task != ObjectiveFunction)
      lb = max(0., lb);

    return lb;
  }

  // *****************************************************************************
  //                            get_function_minorant
  // *****************************************************************************


  int ConeModel::get_function_minorant(MinorantPointer& minorant,
    const AffineFunctionTransformation* aft) {
    if (get_data()->cand_ub_mid != get_data()->modification_id) {
      if (cb_out())
        get_out() << "**** WARNING: ConeModel::get_function_minorant(..): no candidate available" << std::endl;
      return 1;
    }
    MinorantPointer& mnrt = provide_cand_minorant();
    if (!mnrt.valid()) {
      if (cb_out())
        get_out() << "**** WARNING: ConeModel::get_function_minorant(..): provide_cand_minorant failed" << std::endl;
      return 1;
    }

    Real alpha = get_data()->function_factor;
    if ((get_data()->function_task != ObjectiveFunction) && (get_data()->cand_ub <= 0.))
      alpha = 0.;

    if (aft) {
      if (aft->transform_minorant(minorant, mnrt, alpha, true)) {
        if (cb_out())
          get_out() << "**** WARNING: ConeModel::get_function_minorant(..): collecting minorant information failed" << std::endl;
        return 1;
      }
      return 0;
    }

    if (mnrt.get_minorant(minorant, alpha)) {
      if (cb_out())
        get_out() << "**** WARNING: ConeModel::get_function_minorant(..): get_minorant failed" << std::endl;
      return 1;
    }

    return 0;
  }

  // *****************************************************************************
  //                            get_center_minorant
  // *****************************************************************************


  int ConeModel::get_center_minorant(MinorantPointer& minorant,
    const AffineFunctionTransformation* aft) {
    if (get_data()->center_ub_mid != get_data()->modification_id) {
      if (cb_out())
        get_out() << "**** WARNING: ConeModel::get_center_minorant(..): no center available" << std::endl;
      return 1;
    }

    MinorantPointer& mnrt = provide_center_minorant();
    if (!mnrt.valid()) {
      if (cb_out())
        get_out() << "**** WARNING: ConeModel::get_center_minorant(..): provide_center_minorant failed" << std::endl;
      return 1;
    }

    Real alpha = get_data()->function_factor;
    if ((get_data()->function_task != ObjectiveFunction) && (get_data()->center_ub <= 0.))
      alpha = 0.;

    if (aft) {
      if (aft->transform_minorant(minorant, mnrt, alpha, true)) {
        if (cb_out())
          get_out() << "**** WARNING: ConeModel::get_center_minorant(..): aft->transform_minorant information failed" << std::endl;
        return 1;
      }
      return 0;
    }

    if (mnrt.get_minorant(minorant, alpha)) {
      if (cb_out())
        get_out() << "**** WARNING: ConeModel::get_center_minorant(..): get_minorant failed" << std::endl;
      return 1;
    }

    return 0;
  }

  // *****************************************************************************
  //                            adjust_multiplier
  // *****************************************************************************


  int ConeModel::adjust_multiplier(bool& values_may_have_changed) {
    values_may_have_changed = false;
    if (get_data()->function_task != AdaptivePenaltyFunction) {
      //no multipliers of relevance in this class
      return 0;
    }

    int err = 0;
    Real old_function_factor = get_data()->function_factor;
    Real new_function_factor = old_function_factor;
    if (get_data()->sumbundle.active()) {
      //all is taken care of by the sumbundlehandler, but check for changes in function_factor
      assert(bundlehandler);
      bool maybe_changes = false;
      if (bundlehandler->adjust_multiplier(maybe_changes)) {
        if (cb_out())
          get_out() << "**** WARNING: ConeModel::adjust_multiplier(...): bundlehandler->adjust_multiplier() failed" << std::endl;
        err++;
      }
      values_may_have_changed |= maybe_changes;
      if (bundlehandler->get_increase_factor() != 1.) {
        new_function_factor *= bundlehandler->get_increase_factor();
      }
    } else {
      if (local_adjust_multiplier(new_function_factor)) {
        if (cb_out()) {
          get_out() << "**** ERROR ConeModel::adjust_multiplier(.): local_adjuct_multiplier() failed" << std::endl;
        }
        err++;
      }
    }

    if (std::fabs(new_function_factor - old_function_factor) > 1e-10 * old_function_factor) {
      values_may_have_changed = true;
      assert(new_function_factor > 0.);
      get_data()->function_factor = new_function_factor;
      if (get_data()->cand_ub_mid == get_data()->modification_id)
        get_data()->cand_ub_mid++;
      if (get_data()->center_ub_mid == get_data()->modification_id)
        get_data()->center_ub_mid++;
      get_data()->modification_id++;
      get_data()->local_aggregate.clear();
      get_data()->aggregate.clear();
      if (bundlehandler) {
        if (bundlehandler->reset_function_factor(get_data()->function_task, get_data()->function_factor)) {
          if (cb_out()) {
            get_out() << "**** ERROR ConeModel::adjust_multiplier(): bundlehandler->reset_function_factor failed" << std::endl;
          }
          err++;
        }
      }

      if (cb_out(1))
        get_out() << " newmult=" << get_data()->function_factor << std::endl;
    }
    return 0;
  }

  // *****************************************************************************
  //                                clear_model
  // *****************************************************************************

  void ConeModel::clear_model(bool discard_minorants_only) {
    if (get_data()->sumbundle.bundle_size(get_data()->function_task) > 0) {
      if (bundlehandler) {
        bundlehandler->clear_model();
      } else {
        std::map<FunctionTask, Real> fmap;
        fmap[get_data()->function_task] = get_data()->function_factor;
        bundlehandler = new SumBundleHandler(get_data()->sumbundle, fmap);
        bundlehandler->set_cbout(this, 0);
        bundlehandler->clear_model();
        delete bundlehandler;
        bundlehandler = 0;
      }
    }

    get_data()->clear_model(discard_minorants_only);
  }

  // *****************************************************************************
  //                                clear_aggregates
  // *****************************************************************************

  void ConeModel::clear_aggregates() {
    if (get_data()->sumbundle.bundle_size(get_data()->function_task) > 0) {
      if (bundlehandler) {
        bundlehandler->clear_aggregates();
      } else {
        std::map<FunctionTask, Real> fmap;
        fmap[get_data()->function_task] = get_data()->function_factor;
        bundlehandler = new SumBundleHandler(get_data()->sumbundle, fmap);
        bundlehandler->set_cbout(this, 0);
        bundlehandler->clear_aggregates();
        delete bundlehandler;
        bundlehandler = 0;
      }
    }

    get_data()->clear_aggregates();

  }


}
