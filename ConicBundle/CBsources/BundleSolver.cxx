/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/BundleSolver.cxx
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


#include <fstream>
#include <cmath>
#include <cstdlib>
#include <typeinfo>
#include <sstream>
#include "mymath.hxx"
#include "BundleSolver.hxx"
#include "UnconstrainedGroundset.hxx"
#include "LPGroundset.hxx"
#include "LPGroundsetModification.hxx"
#include "BundleHKWeight.hxx"
#include "BundleIdProx.hxx"
#include "BundleDiagonalTrustRegionProx.hxx"
#include "BundleLowRankTrustRegionProx.hxx"
#include "BundleDLRTrustRegionProx.hxx"
#include "BundleDenseTrustRegionProx.hxx"
#include "UQPSolver.hxx"
#include "VariableMetricSVDSelection.hxx"

using namespace CH_Matrix_Classes;

namespace ConicBundle {


// *****************************************************************************
//                                 BundleSolver
// *****************************************************************************

  BundleSolver::BundleSolver(const CBout* cbo,int incr):CBout(cbo,incr)
{
  external_groundset=false;
  groundset=0;
  model=0; 
  terminator=0;
  bundleweight=0;
  Hp=0;
  clockp=0;

  qp_save_cnt=0;
 
  clear();
  set_defaults();
}

// *****************************************************************************
//                                 BundleSolver
// *****************************************************************************

  BundleSolver::BundleSolver(Integer in_dim,BundleModel* bp,const CBout* cbo,int incr):CBout(cbo,incr)
{
  external_groundset=false;
  groundset=0;
  model=0; 
  terminator=0;
  bundleweight=0;
  Hp=0;
  clockp=0; 

  initialize(in_dim,bp);
}

// *****************************************************************************
//                                 BundleSolver
// *****************************************************************************

  BundleSolver::BundleSolver(Groundset* gs,BundleModel* bp,const CBout* cbo,int incr):CBout(cbo,incr)
{
  external_groundset=false;
  groundset=0;
  model=0; 
  terminator=0;
  bundleweight=0;
  Hp=0;
  clockp=0; 
 
  initialize(gs,bp);
}

// *****************************************************************************
//                                 ~BundleSolver
// *****************************************************************************

BundleSolver::~BundleSolver()
{
  if (!external_groundset)
    delete groundset;

  delete terminator; 
  delete bundleweight;
  delete Hp;
  
}

// *****************************************************************************
//                                 set_defaults
// *****************************************************************************

// usually called on construction

void BundleSolver::set_defaults()
{
  //modeleps=0.6; 
  modeleps=0.05; 
  mL=0.1;   
  mN=0.15; //must be > mL!


  check_center_validity=true;

  delete Hp;
  Hp=new BundleIdProx;
  Hp->set_cbout(this);
  
  max_updates=10; 
  use_linval=true;
  if (terminator) terminator->set_defaults();
  if (bundleweight) bundleweight->set_defaults();
}

// *****************************************************************************
//                                   clear
// *****************************************************************************

// usually called on construction, resets all matrices and calls set_defaults

void BundleSolver::clear()
{
  if (!external_groundset)
    delete groundset;
  groundset=0;
  external_groundset=false;

  model=0;
  qp_solver=0;

  delete terminator;
  terminator=new BundleTerminator(this,0);
  terminate=0;

  delete bundleweight;
  bundleweight=new BundleHKWeight(.5,0,this,0);

  
  //HP is initialized in set_defaults()

  clockp=0; 

  point_id=-1;

  center_id=-1;
  center_y.init(0,1,0.);
  center_fid=-1;
  center_ub=0.;
  center_relprec=1e-3;
  center_gid=-1;
  center_gs_val=0.;

  cand_id=-1;
  cand_y.init(0,1,0.);
  cand_fid=-1;
  cand_ub=0.;
  cand_relprec=1e-3;
  cand_gid=-1;
  cand_gs_val=0.;

  model_aggregate.init(new Minorant,cand_fid);
  model_aggregate_id=-1;
  modelval=min_Real;

  weightu=-1;      
  aggr_dnormsqr=-1;

  updatecnt=0;
  sumupdatecnt=0;
  retcode=0;

  initialize_model=true;
  use_cand_for_aggregate=false;

  cand_ub=0.;
  linval=0.;
  cutval=0.;
  augval_lb=min_Real;
  augval_ub=max_Real;

  descent_step=false;
  null_step=false;

  innerit=0;
  suminnerit=0; 
  cntobjeval=0; 

  recomp=0;     
  sumrecomp=0;  
  qpfails=0;    
  sumqpfails=0; 
  modelfails=0;   
  summodelfails=0;
  augvalfails=0;                    
  sumaugvalfails=0;
  oraclefails=0;
  sumoraclefails=0;

  descent_steps=0;        

  shallowcut=0; 
  modelprec=0;

  set_defaults();
  might_be_modified=true;
}

// *****************************************************************************
//                                   clear_fails
// *****************************************************************************

// resets all fail values

void BundleSolver::clear_fails()
{
  recomp=0;     
  sumrecomp=0;  
  qpfails=0;    
  sumqpfails=0; 
  modelfails=0;   
  summodelfails=0;
  augvalfails=0;                    
  sumaugvalfails=0;
  oraclefails=0;
  sumoraclefails=0;
}

// *****************************************************************************
//                                  initialize
// *****************************************************************************

int BundleSolver::initialize(Integer dim,BundleModel* bp)
{
  int err=0;
  if (dim<0) {
    if (cb_out()){
      get_out()<<"**** ERROR BundleSolver::initialize(): dim="<<dim<<" is negative"<<std::endl;
    }
    err=1;
  }
  if (err)
    return err;

  clear();  
  groundset=new UnconstrainedGroundset(dim);
  groundset->set_cbout(this);

  GroundsetModification gsmdf(0);
  gsmdf.add_append_vars(groundset->get_dim());
  Hp->apply_modification(gsmdf);

  model=bp;
  return 0;
}

// *****************************************************************************
//                                  initialize
// *****************************************************************************

int BundleSolver::initialize(Groundset* gp,BundleModel* bp)
{
  int err=0;
  if (gp==0) {
    if (cb_out()){
      get_out()<<"**** ERROR BundleSolver::initialize(): NULL pointer to groundset "<<std::endl;
    }
    err=1;
  }
  if (err)
    return err;

  clear();  
  external_groundset=true;
  groundset=gp;

  GroundsetModification gsmdf(0);
  gsmdf.add_append_vars(groundset->get_dim());
  Hp->apply_modification(gsmdf);

  return set_model(bp);
}

// *****************************************************************************
//                                  initialize
// *****************************************************************************

int BundleSolver::set_model(BundleModel* bp)
{
  if (bp!=model){
    model=bp;
    center_fid=-1;
    cand_fid=-1;
    center_relprec=1e-3;  
    initialize_model=true;
    if (model){
      center_ub=CB_plus_infinity;
      cand_ub=CB_plus_infinity;
      model_aggregate.clear();
    }
    else { 
      center_ub=0.;
      cand_ub=0.;
      model_aggregate.init(new Minorant,cand_fid);
    }
    model_aggregate_id=-1;    
    augval_lb=linval=cutval=modelval=CB_minus_infinity;
    augval_ub=CB_plus_infinity;
  }

  return 0;
}

// *****************************************************************************
//                                   set_bundlweight
// *****************************************************************************

void BundleSolver::set_bundleweight(BundleWeight* bw)
{
  if (bw) {
    delete bundleweight; 
    bundleweight=bw;
  }
  else if (dynamic_cast<BundleHKWeight*>(bundleweight)==0){
    BundleWeight* bw=new BundleHKWeight(.5,bundleweight,this,0);
    delete bundleweight; 
    bundleweight=bw;
  }  
}

// *****************************************************************************
//                                   set_variable_metric
// *****************************************************************************

int BundleSolver::set_variable_metric(int ds)
{
  assert(groundset);
  switch(ds){
  case 0: default: {
    delete Hp; 
    Hp=new BundleIdProx(groundset->get_dim(),1.,this);
    Hp->set_cbout(this);
    initialize_model=true;
    break;
  }
  case 1: {
    delete Hp;
    VariableMetricSVDSelection* vms=new VariableMetricSVDSelection(50,5,0.,this);
    assert(vms);
    Hp=new BundleDiagonalTrustRegionProx(groundset->get_dim(),1.,vms,false,false,this);
    break;
  }
  case 2: {    
    delete Hp;
    VariableMetricSVDSelection* vms=new VariableMetricSVDSelection(50,5,0.,this);
    assert(vms);
    Hp=new BundleDiagonalTrustRegionProx(groundset->get_dim(),1.,vms,true,true,this);
    initialize_model=true;
    break;
  }
  case 3: {    
    delete Hp;
    VariableMetricSVDSelection* vms=new VariableMetricSVDSelection(50,5,0.,this);
    assert(vms);
    Hp=new BundleLowRankTrustRegionProx(groundset->get_dim(),vms,true,this);
    initialize_model=true;
    break;
  }
  case 4: {    
    delete Hp;
    VariableMetricSVDSelection* vms=new VariableMetricSVDSelection(50,5,0.,this);
    assert(vms);
    Hp=new BundleDLRTrustRegionProx(groundset->get_dim(),vms,true,this);
    initialize_model=true;
    break;
  }
  case 5: {    
    delete Hp;
    VariableMetricSVDSelection* vms=new VariableMetricSVDSelection(50,5,0.,this);
    assert(vms);
    Hp=new BundleDenseTrustRegionProx(groundset->get_dim(),vms,true,false,this);
    initialize_model=true;
    break;
  }
  }
  return 0;
}

// *****************************************************************************
//                                   set_prox_diagonal
// *****************************************************************************

int BundleSolver::set_prox_diagonal(const CH_Matrix_Classes::Matrix& insc)
{
  delete Hp;
  Matrix tmp(insc);
  tmp.inv();
  BundleDiagonalTrustRegionProx* dp=new BundleDiagonalTrustRegionProx(tmp);
  Hp=dp;
  Hp->set_cbout(this);
  initialize_model=true;
  return 0;
}

// *****************************************************************************
//                                   set_prox
// *****************************************************************************

int BundleSolver::set_prox(BundleProxObject* Sp)
{
  if (Sp!=0){
    delete Hp;
    Hp=Sp;
    initialize_model=true;
    return 0;
  }
  if (cb_out())
    get_out()<<"**** WARNING BundleSolver::set_prox(BundleProxObject*) was called with 0 pointer"<<std::endl;
  return 1;
}


// *****************************************************************************
//                                set_new_center
// *****************************************************************************

//set a new center point. If the input is NULL then a default
//starting point is constructed.
//returns 0 on success, 1 on failure 

int BundleSolver::set_new_center(const Matrix* yp) 
{ 
  if (cb_out(10)){
    get_out()<<"\n  entering  BundleSolver::set_new_center"<<std::endl;
  }

  bool ychanged=false;
  if ((yp==0)||(yp->dim()!=groundset->get_dim())) {
    //try to use the old center as starting point instead
    if ((yp)&&(yp->dim()!=groundset->get_dim())){
      if (cb_out()){
	get_out()<<"**** WARNING: BundleSolver::set_new_center(...): dim(starting point)="<<yp->dim()<<" != groundset->get_dim()="<<groundset->get_dim()<<", try to use old center intstead "<<std::endl;
      }
    }
    if (center_y.dim()!=groundset->get_dim()){
      ychanged=true;
      center_y.init(groundset->get_starting_point());
      if (center_y.dim()!=groundset->get_dim()){
	if (cb_out()){
	  get_out()<<"**** WARNING: BundleSolver::set_new_center(...): dim(groundset->get_starting point()="<<center_y.dim()<<" != groundset->get_dim()="<<groundset->get_dim()<<", use zero vector instead"<<std::endl;
	}
	center_y.init(groundset->get_dim(),1,0.);
      }
    }
  }
  else {
    if (yp!=&center_y){
      center_y=*yp;
      ychanged=true;
    }
    else // initialize if any of the ids is negative
      ychanged=(center_gid<0)||(center_id<0);
  }

  if (ychanged)
    center_gid=-1;
  int old_center_gid=center_gid;
  int status=groundset->ensure_feasibility(center_gid,center_y,ychanged,Hp,1e-6);
  if (status){
    if (cb_out()){
      get_out()<<"**** ERROR: BundleSolver::set_new_center(...): groundset->ensure_feasibility failed in determinig the center and returned "<<status<<std::endl;
    }
    if (cb_out(10)){
      get_out()<<"\n  leaving  BundleSolver::set_new_center with status="<<status<<std::endl;
    }
    return status;
  }
  if ((ychanged)||(center_gid>old_center_gid)){
    center_id=++point_id;
    center_gs_val=get_gs_aggregate().evaluate(center_id,center_y);
  }
  
  //center_relprec=1e-6;  //simply use the previous value
  if (model){
    status=model->transform()->recompute_center(center_fid,center_ub,center_id,center_y,false,center_relprec);
    cntobjeval++;
  }
  else {
    center_ub=0.;
    center_fid=-1;
  }

  if (cb_out(1)){
    get_out().precision(12);
    get_out()<<"\n recomp center_ub="<<center_ub<<" center_gs_val="<<center_gs_val<<" center_objval="<<center_ub+center_gs_val<<std::endl;
  }

  use_cand_for_aggregate=false;
  if (status){
    if (cb_out()){
      get_out()<<"**** ERROR: BundleSolver::recompute_center(...): recompute_center failed in the new center and returned "<<status<<std::endl;
    }
    if (cb_out(10)){
      get_out()<<"\n  leaving  BundleSolver::set_new_center with status="<<status<<std::endl;
    }
    return status;
  }

  
  if ((model)&&(model->transform()->update_model(BundleModel::new_subgradient,center_id,center_y,center_id,center_y,max(1e-6,1e-3*std::fabs(center_ub)),*Hp))){
    if (cb_out()) {
      get_out()<<"**** WARNING: BundleSolver::solve(): update_model for new_center returned an error\n";  
    }
  }

  if (cb_out(10)){
    get_out()<<"\n  leaving  BundleSolver::set_new_center with status="<<status<<std::endl;
  }

  return status;
}

// *****************************************************************************
//                               variable_metric
// *****************************************************************************


  int BundleSolver::variable_metric(const CH_Matrix_Classes::Indexmatrix* new_indices )
{
  if (!Hp->employ_variable_metric())
    return 0;

  MinorantPointer gaggr=model_aggregate;
  groundset->get_gs_aggregate().get_minorant(gaggr);

  Real dummy_offset;
  Matrix global_aggregate(center_y.dim(),1); chk_set_init(global_aggregate,1);
  gaggr.get_minorant(dummy_offset,global_aggregate,0);

  Real local_weightu=weightu<0?1.:weightu;
  Real aggrval=gaggr.evaluate(center_id,center_y);
  Real model_maxviol = center_ub+center_gs_val-aggrval;
  if (model_maxviol<1e-12*(1.+std::fabs(center_ub))){
    Real aggrdnormsqr=gaggr.dual_norm_squared()/local_weightu;
    model_maxviol = max(1e-12*(1.+std::fabs(center_ub)),(1-mN)*(center_ub+center_gs_val-(aggrval-aggrdnormsqr))-.5*aggrdnormsqr);
  }
  
  if (cb_out(2))
    get_out()<<" variable_metric: wu="<<local_weightu<<" mv="<<model_maxviol<<std::endl;

  if (Hp->apply_variable_metric(groundset,
				model,
				global_aggregate,
				center_id,
				center_y,
				descent_step,
				local_weightu,
				model_maxviol,
				new_indices)){
    if (cb_out()) 
      get_out()<<"**** WARNING BundleSolver::variable_metric(): Hp->apply_variable_metric(...) failed "<<std::endl;
    return 1;
  }
  //if (weightu!=bundleweight->get_weight()){
  //  bundleweight->set_next_weight(weightu);
  //}
  if (cb_out(2)){
    Real aggrnorm=std::sqrt(gaggr.dual_norm_squared());
    Real aggrHnorm=std::sqrt(Hp->dnorm_sqr(gaggr));
    get_out()<<" norm(aggr)="<<aggrnorm<<" Hnorm(aggr)="<<aggrHnorm;
    if (aggrnorm>eps_Real) 
      get_out()<<" red1="<<aggrHnorm/aggrnorm;
    Matrix tmpvec(global_aggregate);
    Hp->apply_Hinv(tmpvec);
    Real aggrH2norm=norm2(tmpvec);
    get_out()<<" H2norm(aggr)"<<aggrH2norm;
    if ((aggrnorm>eps_Real)&&(aggrH2norm>eps_Real)) 
      get_out()<<" red2="<<aggrH2norm/aggrnorm<<" ip="<<ip(tmpvec,global_aggregate)/aggrH2norm/aggrnorm;
    get_out()<<std::endl;
  }


  return 0;
}

// *****************************************************************************
//                                 eval_augmodel
// *****************************************************************************

int BundleSolver::eval_augmodel(Real& augbound,
				Integer& center_ub_fid,
				Real& center_ub, 
				Real relprec,
				Real center_gs_val)
{
  //for testing
  //Integer dummy;
  //old_model_aggregate.clear();
  //sbm_transform()->get_model_aggregate(dummy,old_model_aggregate);

  //initialize blocks and variable dimensions
  CH_Tools::Microseconds coeff_start;
  if (clockp)
    coeff_start=clockp->time();

  //qp_solver->set_cbout(this,6);
  //qp_solver->QPset_proxterm(Hp); 

  qp_solver->clear_model_data_ptr();
  if (model->transform()->start_augmodel(*qp_solver,cand_id,cand_y)){
    if (cb_out()){
      get_out()<<"\n**** WARNING BundleSolver::eval_augmodel(...): start_augmodel failed"<<std::endl;
    }
    return -1;
  }
  assert(qp_solver->get_model_data_ptr()!=0);

  CH_Tools::Microseconds solve_start;
  if (clockp){
    solve_start=clockp->time();
    QPcoeff_time+=solve_start-coeff_start;
  }

  //--- call the quadratic program solver
  int status=qp_solver->QPsolve(center_y,augbound,max(augbound+eps_Real*std::fabs(augbound),center_ub+center_gs_val),relprec,Hp,groundset->get_gs_aggregate(),groundset->set_yfixed());
  if (status){
    if (cb_out()) 
      get_out()<<"**** ERROR BundleSolver::eval_augmodel(): qp_solver->QPSolver failed and returned "<<status<<std::endl;
    if (clockp){
      QPsolve_time+=clockp->time()-solve_start;
      evalaugmodel_time+=clockp->time()-coeff_start;
    }
  }
  if (clockp)
    QPsolve_time+=clockp->time()-solve_start;


  //--- get the objective values and the aggregate subgradient
  //    if the cone multiplier was increased recompute the solution

  Integer it=0;
  do {
    Real oldaugbound=augbound;
    augbound=qp_solver->QPget_lower_bound();
    if (augbound<oldaugbound-1e-12*(1.+CH_Matrix_Classes::abs(oldaugbound))){
      if (cb_out()){
	get_out().precision(16);
	get_out()<<"**** WARNING: BundleSolver::eval_augmodel(): qp_solver->solve() returned "<<augbound<<" below old augval_lb="<<oldaugbound<<std::endl;
      }
    }

    bool penalty_changed=false;
    bool keep_penalty_fixed=true;
    if (it++<1) keep_penalty_fixed=false;
    CH_Tools::Microseconds make_aggr_start;
    if (clockp)
      make_aggr_start=clockp->time();
    int linstat=model->transform()->make_model_aggregate(penalty_changed,keep_penalty_fixed);
    if (clockp)
      make_aggr_time+=clockp->time()-make_aggr_start;
    if (linstat){
      if (cb_out()) get_out()<<"**** WARNING: BundleSolver::eval_augmodel(): make_model_aggregate failed and returned "<<linstat<<std::endl;
      status|=linstat;
    }

    if (!penalty_changed) break;
    int restat=model->transform()->recompute_center(center_ub_fid,center_ub,center_id,center_y);
    if  (restat<0){
      if (cb_out()) get_out()<<"**** WARNING: BundleSolver::eval_augmodel(): recompute_center failed to obtain required precision and returned "<<restat<<std::endl;
    }
    if  (restat>0){
      if (cb_out()) get_out()<<"**** WARNING: BundleSolver::eval_augmodel(): recompute_center failed and returned "<<restat<<std::endl;
      status|=restat;
    }
    if (linstat||(restat>0))
      break;
    if (clockp)
      solve_start=clockp->time();
    int status2=qp_solver->QPresolve(augbound,center_ub+center_gs_val,relprec);
    if (clockp)
      QPsolve_time+=clockp->time()-solve_start;
    if (status){
      if (cb_out()) get_out()<<"**** WARNING: BundleSolver::eval_augmodel(): qp_solver->QPresolve() failed and returned "<<status<<std::endl;
      status|=status2;
    }

    //qp_mfile_data(y,Hp);
  }while (true);
  if (clockp)
    evalaugmodel_time+=clockp->time()-coeff_start;

 
  if (cb_out(10)){
    get_out()<<"\n  leaving  BundleSolver::eval_augmodel with status="<<status<<std::endl;
  }
  return status;
}

// *****************************************************************************
//                                 reeval_augmodel
// *****************************************************************************

int BundleSolver::reeval_augmodel(Real& augbound,
				  Integer& center_ub_fid,
				  Real& center_ub, 
				  Real relprec,
				  Real center_gs_val,
				  const MinorantPointer& delta_gs_subg,
				  const Indexmatrix& delta_index)
{

  CH_Tools::Microseconds solve_start;
  if (clockp){
    solve_start=clockp->time();
  }
  //--- call the quadratic program solver

  int status=qp_solver->QPupdate(center_y,augbound,center_ub+center_gs_val,relprec,
				 Hp,groundset->get_gs_aggregate(),groundset->set_yfixed(),
				 delta_gs_subg,delta_index);

  //status=solver.solve(Q,c,offset);
  if (clockp)
    QPsolve_time+=clockp->time()-solve_start;

  if (status){
    if (cb_out()) get_out()<<"BundleSolver::reeval_augmodel(): solver->QPupdate() failed and returned "<<status<<std::endl;
  }


  //--- get the objective values and the aggregate subgradient
  //    if the cone the cone multiplier was increased recompute the solution
  Integer it=0;
  do {
    Real oldaugbound=augbound;
    augbound=qp_solver->QPget_lower_bound();
    if (oldaugbound>augbound){
      if (cb_out()){
	get_out().precision(16);
	get_out()<<"**** WARNING: BundleSolver::reeval_augmodel(): solver.solve() returned primalval="<<augbound<<" below old augval_lb="<<oldaugbound<<std::endl;
      }
    }
    bool penalty_changed=false;
    bool keep_penalty_fixed=true;
    if (it++<1) keep_penalty_fixed=false;
    CH_Tools::Microseconds make_aggr_start;
    if (clockp)
      make_aggr_start=clockp->time();
    int linstat=model->transform()->make_model_aggregate(penalty_changed,keep_penalty_fixed);
    if (clockp)
        make_aggr_time+=clockp->time()-make_aggr_start;
    if (linstat){
      if (cb_out()) get_out()<<"**** WARNING: BundleSolver::reeval_augmodel(): make_model_aggreagte failed and returned "<<linstat<<std::endl;
      status|=linstat;
    }
    if (!penalty_changed) break;
    int restat=model->transform()->recompute_center(center_ub_fid,center_ub,center_id,center_y);
    if  (restat<0){
      if (cb_out()) get_out()<<"**** WARNING: BundleSolver::reeval_augmodel(): recompute_center failed to obtain required precision and returned "<<restat<<std::endl;
    }
    if  (restat>0){
      if (cb_out()) get_out()<<"**** WARNING: BundleSolver::reeval_augmodel(): recompute_center failed and returned "<<restat<<std::endl;
      status|=restat;
    }
    if (linstat||(restat>0))
      break;
    if (clockp)
      solve_start=clockp->time();
    int status=qp_solver->QPresolve(augbound,center_ub+center_gs_val,relprec);
    if (clockp)
      QPsolve_time+=clockp->time()-solve_start;
    if (status){
      if (cb_out()) get_out()<<"**** WARNING: BundleSolver::reeval_augmodel(): solver.QPresolve() failed and returned "<<status<<std::endl;
    }
      
    //qp_mfile_data(y,Hp);
  } while (true);
  

  if (cb_out(10)){
    get_out()<<"\n  leaving  BundleSolver::reeval_augmodel with status="<<status<<std::endl;
  }
  if (clockp)
    evalaugmodel_time+=clockp->time()-solve_start;

  return status;
}


// *****************************************************************************
//                                 solve_model
// *****************************************************************************

// loops over model optimization and lagrange updates till model optimizer
// is sufficiently close to a feasible model solution.

int BundleSolver::solve_model()
{     
  assert(model);
  updatecnt=0;
  retcode=0;
  qpfails=0;         //is increased whenever the quadratic bundle subproblem
                     //could not be solved to sufficient precision
  modelfails=0;

  MinorantPointer delta_groundset_aggregate;
  Indexmatrix delta_index;

  Real oldaugval= CB_minus_infinity;
 
  //--- iterate till model is sufficiently precise
  do {
   
    updatecnt++;
    sumupdatecnt++;
    if (cb_out(1)){
      get_out()<<"  upd"<<updatecnt<<":";
    }
    



    //--- solve the quadratic model or update it
    //this will result in an increase augval and even center_ub
    //(center_ub may be increased due to an exact penalty approach)
    oldaugval=max(oldaugval,augval_lb);   
    int errcode;
    if (updatecnt==1){
      Real localprec=0.9;
      //if (center_ub+center_gs_val-augval_lb>1e-2*(fabs(center_ub+center_gs_val)+1.))
      //	localprec=1e-3;
      errcode=eval_augmodel(augval_lb,
			    center_fid,center_ub,		
			    localprec,
			    center_gs_val);
    }
    else {
      errcode=reeval_augmodel(augval_lb,
			      center_fid,center_ub,
			      0.1,
			      center_gs_val,
			      delta_groundset_aggregate,
			      delta_index);
    }
    
    retcode=retcode || errcode;
    
    if (errcode){
      if (cb_out()) get_out()<<"**** WARNING BundleSolver::solve_model(): solving quadratic model failed "<<errcode<<std::endl;
      qpfails++;
      sumqpfails++;
    }
    
    //--- get the aggregate minorant of the model
    model_aggregate.clear();
    if (model->transform()->get_model_aggregate(model_aggregate_id,model_aggregate)){
      if (cb_out()) 
	get_out()<<"**** ERROR BundleSolver::solve_model(): no augmented model solution available"<<std::endl;
      terminate=terminator->check_termination(this);
      return 1;
    }

    // //TESTING BEGIN
    // {
    //   Real tmpdummy;
    //   Matrix tmpvec(center_y.rowdim(),1,0.);
    //   model_aggregate.get_minorant(tmpdummy,tmpvec,0,-1.,false);
    //   qp_solver->get_model_data_ptr()->add_modelx_aggregate(tmpdummy,tmpvec);
    //   if (std::fabs(tmpdummy)+norm2(tmpvec)>1e-8){
    // 	std::cout<<"**** ERROR BundleSolver::solve_model(): model aggregates differ: diffoffset="<<tmpdummy<<" norm2(diffvec)="<<norm2(tmpvec);
    // 	std::cout<<" diffinds=";
    // 	for(Integer i=0;i<tmpvec.dim();i++){
    // 	  if (std::fabs(tmpvec(i))>1e-8)
    // 	    std::cout<<" ("<<i<<","<<tmpvec(i)<<")";
    // 	}
    // 	std::cout<<std::endl;
    //   }
    // }
    // //TESTING END

    //--- determine step and compute the value of the linearized and the augmented model in cand_y 
    //Real relprec=max(eps_Real,min(1e-3,max(0.01,mN-mL)*(center_ub+center_gs_val-linval)/(fabs(center_ub+center_gs_val)+1.)));
    Real relprec=max(eps_Real,min(0.1,max(0.01,mN-mL)*(center_ub+center_gs_val-linval)/(fabs(center_ub+center_gs_val)+1.)));

    if(groundset->candidate(cand_gid,cand_y,cand_gs_val,linval,augval_lb,augval_ub,aggr_dnormsqr,center_y,center_ub+center_gs_val,model_aggregate,Hp,&delta_groundset_aggregate,&delta_index,relprec)){
      if (cb_out())
	get_out()<<"**** ERROR BundleSolver::solve_model(): candidate(...) failed"<<std::endl;
      qpfails++;
      sumqpfails++;
      terminate=terminator->check_termination(this);
      return 1;
    }
    cand_id=++point_id;

#ifdef CONICBUNDLE_DEBUG
    {
      Integer dummy;
      MinorantPointer aggr=groundset->get_gs_aggregate();
      if (model)
	model->transform()->get_model_aggregate(dummy,aggr);
      if((std::fabs(aggr.evaluate(cand_id,cand_y)-linval)>relprec*(1+std::fabs(linval)))&&(cb_out())){
	get_out().precision(12);
	get_out()<<"**** WARNING BundleSolver::solve_model(): linval="<<linval<<" is too far from aggregate value="<<aggr.evaluate(cand_id,cand_y)<<" for relprec="<<relprec<<std::endl;
      }
    }
#endif

    //--- if there were changes in cand_y compute the true model value
    if (((delta_index.rowdim()>0)&&(modeleps>0))||(augval_lb>center_ub+center_gs_val)){
      int modelerr=model->transform()->eval_model(cutval,cand_id,cand_y,0.1*(max(eps_Real*(1+std::fabs(linval)),center_ub+center_gs_val-linval))/(fabs(center_ub+center_gs_val)+1.));
      if (modelerr){
	if (cb_out()) get_out()<<"**** WARNING BundleSolver::solve_model(): evaluation of cutting model failed: "<<modelerr<<std::endl;
	modelfails++;
	summodelfails++;
      }
      if (cb_out(1)){
	get_out().precision(8);
	get_out()<<" modelval="<<cutval;
	get_out()<<" gsval="<<cand_gs_val;
      }
      cutval+=cand_gs_val;
      if (cutval<linval) {
	if (cb_out()&&(cutval<linval-1e-10*(std::fabs(linval)+1.))){
	  get_out()<<"**** WARNING BundleSolver::solve_model(): something is strange, agregate value linval="<<linval<<">"<<cutval<<"=cutval which should be the max over linval (value of the aggregate) and further linear minorants..."<<std::endl;
	}
	cutval=linval;
      }
    }
    else {
      cutval=linval;
    }
    modelprec=(cutval-linval)/max(center_ub+center_gs_val-linval,1e-16);
    
    //--- output some information on current iteration
    if (cb_out(1)){
      get_out().precision(10);
      get_out()<<" augval_lb="<<augval_lb;
      get_out()<<" augval_ub="<<augval_ub;
      get_out().precision(8);
      get_out()<<" flin="<<linval;
      if (delta_index.rowdim()>0) {
	get_out().precision(8);get_out()<<" fhat="<<cutval;
	get_out()<<" (";get_out().precision(2);
	get_out()<<modelprec<<")";
      }
      get_out().precision(2);
      get_out()<<" n2="<<aggr_dnormsqr;
      get_out()<<" d2="<<norm2(cand_y-center_y);
      get_out().precision(6);get_out()<<std::endl;
    }
    
    //--- check for termination
    if (use_linval){
      modelval=linval;
    }
    else {
      modelval=cutval;
    }

    terminate=terminator->check_termination(this);
    
  } while (
	   (terminate==0) &&
	   (qp_solver->QPsupports_updates()) &&
	   (modeleps>0) &&
	   (augval_lb>oldaugval+eps_Real*(fabs(oldaugval)+1.)) &&
	   (cutval-linval>modeleps*(center_ub+center_gs_val-linval)) &&  //modelprecision
	   ((max_updates<0)||(updatecnt<max_updates))
	   );
  
  return retcode;
}


// *****************************************************************************
//                                 solve
// *****************************************************************************

// performs null steps till a serious step is encountered or maxsteps is reached

int BundleSolver::solve(int maxsteps,bool stop_at_descent_steps)
{
  if (cb_out(10)){
    get_out()<<"\n  entering  BundleSolver::solve"<<std::endl;
  }

  if (model)
    model->transform()->set_cbout(this);

  if (!groundset) {
    if (cb_out()) 
      get_out()<<"**** ERROR: BundleSolver::solve(...): the groundset specification is missing"<<std::endl;
    return 1;
  }
  else {
    groundset->set_cbout(this);
  }

  if (!Hp) {
    if (cb_out()) get_out()<<"**** ERROR: BundleSolver::solve(...): the specification of the quadratic term is missing"<<std::endl;
    return 1;
  }
  else {
    Hp->set_cbout(this);
  }
   
  if (!bundleweight) { 
    if (cb_out()) get_out()<<"**** ERROR: BundleSolver::solve(...): no routine for choosing bundleweight specified"<<std::endl;
    return 1;
  }
  else {
    bundleweight->set_cbout(this);
  }

  if (!terminator) { 
    if (cb_out()) get_out()<<"**** ERROR: BundleSolver::solve(...): no routine for checking termination specified"<<std::endl;
    return 1;
  }
  else {
    terminator->set_cbout(this);
  }
  
  terminate=0;
  might_be_modified=true;

  //--- compute null steps till a serious step is achieved
  do{

    if (might_be_modified){

      might_be_modified=false;
      
      //--- check whether the center is still valid and udate it if not
      if (!groundset->is_feasible(center_gid,center_y)){
	descent_step=false;
	null_step=false;
	int retval=set_new_center(0);
	if (retval){
	  if (cb_out()){
	    get_out()<<"**** ERROR: BundleSolver::solve(...): set_new_center failed and returned "<<retval<<std::endl;
	  }
	  return retval;
	}
	initialize_model=true;
      }
      else if (model) {
	if (model->transform()->center_modified(center_fid,center_id)){
	  descent_step=false;
	  null_step=false;
	  int retval=set_new_center(&center_y);
	  if (retval){
	    if (cb_out()){
	      get_out()<<"**** ERROR: BundleSolver::solve(...): set_new_center failed and returned "<<retval<<std::endl;
	    }
	    return retval;
	  }
	}
      }
      else {
	descent_step=false;
	null_step=false;
	assert(center_ub==0.);
	assert(center_fid==-1);
      }
      
      //--- check whether the model_aggregate is still valid and update if not
      if ((model)&&(model->transform()->model_aggregate_modified(model_aggregate_id))){
	descent_step=false;
	null_step=false;
	int retval;
	if ((use_cand_for_aggregate)&&(groundset->is_feasible(cand_gid,cand_y))){
	  retval=model->transform()->provide_model_aggregate(cand_id,cand_y);
	}
	else {
	  retval=model->transform()->provide_model_aggregate(center_id,center_y);
	}
	if (retval){
	  if (cb_out()){
	    get_out()<<"**** ERROR: BundleSolver::solve(...): provide_model_aggregate() failed and returned "<<retval<<std::endl;
	  }
	  return retval;
	}
	model_aggregate.clear();
	retval=model->transform()->get_model_aggregate(model_aggregate_id,model_aggregate);
	if (retval){
	  if (cb_out()){
	    get_out()<<"**** ERROR: BundleSolver::solve(...): get_model_aggregate() failed and returned "<<retval<<std::endl;
	  }
	  return retval;
	}
	initialize_model=true;
      }    
    } //endif (might_be_modified)
    
    if (initialize_model){
      innerit=0;
    }

    //--- initialize or update the metric
    if (((initialize_model)||(descent_step))&&(Hp->employ_variable_metric()))
      variable_metric();
      
    descent_step=false;
    null_step=false;
      
    //--- (re)initialize bundleweight 
    if (weightu<0) {
      Hp->set_weightu(1.);
      MinorantPointer aggregate=groundset->get_gs_aggregate();
      if (model)
	model_aggregate.get_minorant(aggregate);
      bundleweight->init(Hp->dnorm_sqr(aggregate),groundset,model);
      initialize_model=true;
    }
    if (fabs(weightu-bundleweight->get_weight())>1e-10*fabs(weightu)){
      initialize_model=true;
    }
    
    //--- complete initialization of the initial model
    if (initialize_model){
      
      //--- g(.)=linval+<subg-eta,.-y> is a minorant of f contained in the model 
      //    compute a lower bound on the model and the augmented model value
      //    for linval and augval_lb, respectively by means of g
      //    (minimizer of g(.)+||.-y||_H^2/2. is  cand_y = y- H^{-1}subg)
      recomp=0;
      terminate=0;
      augvalfails=0;
      oraclefails=0;  
    }
    
    innerit++;
    suminnerit++;
    
    weightu=bundleweight->get_weight();
    Hp->set_weightu(weightu);
    
    if (cb_out(1)){
      get_out()<<" ir"<<innerit;
      if (clockp) {get_out()<<"("<<*clockp<<")";}
      get_out()<<": u=";get_out().precision(4);get_out()<<weightu<<std::endl;
    }

    qp_solver=groundset->get_qp_solver(qp_solves_model_without_gs,Hp);
    if (qp_solver==0){
      if (cb_out())
	get_out()<<"**** ERROR BundleSolver::solve(): groundset->get_qp_solver(..) failed"<<std::endl;
      if (cb_out(10)){
	get_out()<<"\n leaving  BundleSolver::solve: return value 1"<<std::endl;
      }    
      return 1;
    }
    
    //--- after a serious step or a change of the weight
    //    provide a lower bound for the augmented value
    //    and a guess for the groundset aggregate
    Real lastaugval=CB_minus_infinity;
    if (
        (model==0)||
	(
	 (qp_solves_model_without_gs)&&
	 ((initialize_model)||(bundleweight->weight_changed())
	  )
	 )
	){
      augval_lb=CB_minus_infinity;
      if(groundset->candidate(cand_gid,cand_y,cand_gs_val,linval,augval_lb,augval_ub,aggr_dnormsqr,center_y,center_ub+center_gs_val,model_aggregate,Hp)){
	if (cb_out())
	  get_out()<<"**** ERROR BundleSolver::solve(): candidate(...) failed"<<std::endl;
	if (cb_out(10)){
	  get_out()<<"\n leaving  BundleSolver::solve: return value 1"<<std::endl;
	}    
	return 1;
      }
      cand_id=++point_id;
      modelval=cutval=linval; 
      augval_lb=linval+0.999*(aggr_dnormsqr/2.); 
      augval_ub=center_ub; 
      initialize_model=false;
    
      //check whether this candidate solves the model to sufficient precision already
      if (model){
	int modelerr=model->transform()->eval_model(cutval,cand_id,cand_y,max(eps_Real,0.1*(center_ub+center_gs_val-linval)/(fabs(center_ub+center_gs_val)+1.)));
	cutval+=cand_gs_val;
	if (modelerr){
	  if (cb_out()) get_out()<<"**** WARNING BundleSolver::solve(): evaluation of cutting model failed: "<<modelerr<<std::endl;
	  modelfails++;
	  summodelfails++;
	}
      }
      else {
	cand_ub=cutval=linval;
      }
      modelprec=(cutval-linval)/max(center_ub+center_gs_val-linval,eps_Real);
      if ((modelprec<-1e-10*(std::fabs(linval)+1))&&(cb_out())){
	get_out()<<"**** WARNING BundleSolver::solve(): something is strange, agregate value linval="<<linval<<">"<<cutval<<"=cutval which should be the max over linval (value of the aggregate) and further linear minorants..."<<std::endl;
																												    cutval=linval;
      }
      if (use_linval){
	modelval=linval;
      }
      else {
	modelval=cutval;
      }

      //check for termination    
    
      if ((model)||(cutval-linval>0.001*(center_ub+center_gs_val-linval))){  //modelprecision

	//--- solve quadratic model
	lastaugval=augval_lb;
	if (model)
	  solve_model();
	else {
	  if (cb_out()) 
	    get_out()<<"**** WARNING BundleSolver::solve(): internal error, the solver should never get here "<<std::endl;
	}
      }
      else {
	terminate=terminator->check_termination(this);
	//this is only to avoid an error message below
	lastaugval=linval-10*eps_Real*(fabs(linval)+1);
      }

    }
    else {
      //--- solve quadratic model
      assert(model!=0);
      if ((initialize_model)||(bundleweight->weight_changed())){
	//provide a guess for augval_lb
	MinorantPointer mp(model_aggregate);
	groundset->get_gs_minorant().get_minorant(mp,1.);
	//the next +.1 avoids a warning if augval is too exact to be increased
	augval_lb=mp.evaluate(center_id,center_y)-(.5+.1)*Hp->dnorm_sqr(mp);
	//don't care too much about the candidate, but there needs to be one
	cand_id=++point_id;
	cand_y=center_y;
	initialize_model=false;
      }
      lastaugval=augval_lb;
      solve_model();
    }
    
    if (augval_lb<=lastaugval+eps_Real*(fabs(lastaugval)+1)){
      if (cb_out(0)){
    	get_out().precision(12);
    	get_out()<<"**** WARNING: BundleSolver::solve(): could not increase augval_lb="<<augval_lb<<" lastaugval="<<lastaugval<<"\n";
    	get_out()<<"             maybe precision requirements are too high or\n";
    	get_out()<<"             the weight was decreased instead of increased during consecutive null steps"<<std::endl;
    	get_out().precision(4);
      }
      augvalfails++;
      sumaugvalfails++;
    }
        
    if (terminate) {  //termination checked in solve_model
      if (cb_out(0))
	print_line_summary(get_out());
      break;
    }

    //---- determine nullstep_bound and descent_bound
    Real center_objval = center_ub+center_gs_val;
    Real nullstep_bound = center_objval-mN*(center_objval-modelval);
    Real descent_bound  = center_objval-mL*(center_objval-modelval);
    assert(mN<.5);
    Real model_maxviol = (1-mN)*(center_objval-modelval)-.5*aggr_dnormsqr;

    if ((augval_ub<.99*center_objval+.01*augval_lb)&&(augval_ub>0.01*center_objval+.99*augval_lb)) {
      nullstep_bound=mN*augval_ub+(1.-mN)*center_objval;
      descent_bound=mL*augval_ub+(1.-mL)*center_objval;
    }

    if (descent_bound<cutval){
      if (cb_out(0))
	get_out()<<"\n    imprecise model evaluation enforces null step (cure: increase update iterations or model precision)"<<std::endl;
    }
  
    //--- evaluate function at the candidate cand_y
    
    cand_relprec=max(eps_Real,min(1e-3,max(0.01,mN-mL)*(center_objval-nullstep_bound)/(fabs(center_objval)+1.)));
    int status1=0;
    if (model) {
      status1=model->transform()->eval_function(cand_fid,cand_ub,cand_id,cand_y,nullstep_bound-cand_gs_val,max(eps_Real,cand_relprec*(1.+std::fabs(center_objval))/(1.+std::fabs(center_ub))));
      cntobjeval++;
      if (status1<0) {//no accurate solution, yet it may suffice to continue
	if (cb_out()) get_out()<<"**** WARNING: BundleSolver::solve(): function evaluation failed to provide a sufficiently accurate solution and returned: "<<status1<<std::endl;
	oraclefails++;
	sumoraclefails++;
      }
      else if (status1>0) { //no new solution information at all, return with error
	if (cb_out()) get_out()<<"**** ERROR: BundleSolver::solve(): function evaluation failed beyond numerical problems and returned "<<status1<<std::endl;
	if (cb_out(10)){
	  get_out()<<"\n leaving  BundleSolver::solve: return value 1"<<std::endl;
	}    
	return 1;
      }
    }

    //--- if the function id increased check validity of center and aggregate
    if (cand_fid>center_fid){
      assert(model);
      if ((model->transform()->center_modified(center_fid,center_id))||
	  (model->transform()->model_aggregate_modified(model_aggregate_id))){
	might_be_modified=true;
	if (cb_out()){
	  get_out()<<"**** ERROR: BundleSolver::solve():   center or aggregate invalidated by dynamic function changes, candidate is discarded and recomputed"<<std::endl;
	}
	continue;
      } 
    }
    
    //--- check for potential problems with relative precision of all kinds
    Real cand_objval=cand_ub+cand_gs_val;
#ifdef CONICBUNDLE_DEBUG
    {
      Integer dummy;
      MinorantPointer cand_mnrt=groundset->get_gs_minorant();
      if (model)
	model->transform()->get_function_minorant(dummy,cand_mnrt);
      if ((std::fabs(cand_mnrt.evaluate(cand_id,cand_y)-cand_objval)>cand_relprec*(1+std::fabs(cand_objval)))&&(cb_out())){
	get_out().precision(12);
	get_out()<<"**** WARNING BundleSolver::solve(): cand_objval="<<cand_objval<<" is too far from minorant value="<<cand_mnrt.evaluate(cand_id,cand_y)<<" for relprec="<<cand_relprec<<std::endl;
      }
    }
#endif
    if (status1==0){
      if(cand_objval>descent_bound){  //upper bound will produce a null step
	if (cand_objval<=center_objval-mL*(center_objval-cutval)){
	  if (cb_out(0)){
	    get_out()<<"    null step might be due to insuffcient model precision"<<std::endl;
	  }
	}

	//check if new cut is good enough to improve the model ...

	if (center_objval-cand_objval>.8*(center_objval-cutval)){ 
	  //subgradient won't yield much improvement
	  if (cb_out(0)){
	    get_out()<<"  shallow cut (subgradient won't yield much improvement)";get_out().precision(3);
	  }
	  shallowcut++;
	}
      }
    }
    
    //--- output solution info
    if (cb_out(1)){
      get_out()<<"  nreval="<<cntobjeval;
      get_out().precision(8);get_out()<<" center_objval="<<center_objval;
      get_out().precision(8);get_out()<<" cand_objval="<<cand_objval;
      if (use_linval) {
	get_out().precision(8);get_out()<<" lin_model=";
      }
      else {
	get_out().precision(8);get_out()<<" full_model=";
      }
      get_out()<<modelval;
      get_out().precision(2);get_out()<<" model_prec="<<(center_objval-cand_objval)/(center_objval-modelval);
      get_out().precision(4);get_out()<<std::endl;
    }
    
    if (cutval>cand_objval+max(eps_Real,1e-3*cand_relprec)*(fabs(cand_objval)+1.)){  //minorant cuts above new function value! 
      if (cb_out()) {
	get_out().precision(12);
	get_out()<<"**** ERROR: BundleSolver::solve(): new objective value "<<cand_objval<<" is smaller than minorizing model value ="<<cutval<<" \n";
	get_out()<<"****        maybe the objective value was not computed to sufficient precision\n";
	get_out()<<"****        or the sugbradient information provided by a previous oracle call was wrong.\n";
	get_out().precision(4);
	return -1;
      }
    }

    //--- check for descent step; if so, move to this point and possibly exit the loop
    if ((!status1)&&(cand_objval>=nullstep_bound)){
      assert(model);
      //null step, update the model
      null_step=true;
      descent_step=false;
      Integer dummy;
      MinorantPointer cand_mnrt=groundset->get_gs_minorant();
      if (model)
	model->transform()->get_function_minorant(dummy,cand_mnrt);
      MinorantPointer aggr=groundset->get_gs_aggregate();
      if (model)
	model->transform()->get_model_aggregate(dummy,aggr);
      //--- choose a new weight (greater or equal the current one)
      bundleweight->nullstep_update(cand_objval,center_objval,
				    modelval,center_y,cand_y,
				    cand_mnrt,aggr,nullstep_bound,aggr_dnormsqr,Hp);
      if ((model)&&(model->transform()->update_model(BundleModel::null_step,center_id,center_y,cand_id,cand_y,model_maxviol,*Hp))){
	if (cb_out()) {
	  get_out()<<"**** WARNING: BundleSolver::solve(): update_model for null_step returned an error\n";  
	}
      }
      use_cand_for_aggregate=true;

      if (cb_out(1))
	print_line_summary(get_out());
    }
    else {
      //descent step (serious step)
      descent_step=true;
      null_step=false;
      bundleweight->descent_update(cand_objval,center_objval,modelval,center_y,cand_y,aggr_dnormsqr,Hp);
      if ((model)&&(model->transform()->update_model(BundleModel::descent_step,center_id,center_y,cand_id,cand_y,model_maxviol,*Hp))){
	if (cb_out()) {
	  get_out()<<"**** WARNING: BundleSolver::solve(): update_model for descent_step returned an error\n";  
	}
      }
      center_id=cand_id;
      center_y=cand_y;
      center_fid=cand_fid;
      center_ub=cand_ub;
      center_relprec=cand_relprec;
      center_gid=cand_gid;
      center_gs_val=cand_gs_val;
      center_objval=center_ub+center_gs_val;
      descent_steps++;
      initialize_model=true;
      use_cand_for_aggregate=false;

      if (cb_out(0))
	print_line_summary(get_out());

      /*
      if (cb_out(2)){
	Indexmatrix sind=sortindex(center_y,false);
	get_out()<<" sorted_center_y="<<transpose(center_y(sind));
	get_out()<<" sind="<<transpose(sind);
      }
      */
      
      if (stop_at_descent_steps)
	break;
      else {
	continue;
      }
    }

      
    //--- check for validity of the old objective value
    //    by means of the new subgradient = new minorant

    if ((model)&&((check_center_validity)||(null_step && (augval_lb>center_ub+center_gs_val)))){
      bool cand_minorant_is_below=true;
      int status2=0;
      if (check_center_validity){
	status2=model->transform()->check_center_validity_by_candidate(cand_minorant_is_below,center_id,center_y);
	if (status2) {  //no new solution information at all, return with error
	  if (cb_out()) 
	    get_out()<<"**** WARNING: BundleSolver::solve(): checking center value by new subgradient failed"<<std::endl;
	}
	else if (!cand_minorant_is_below){
	  if (cb_out()) {
	    get_out().precision(12);
	    get_out()<<"**** WARNING: BundleSolver::solve(): new subgradient yields greater objective value in center (center_ub="<<center_ub<<"),\n";
	    get_out()<<"****          maybe the old evaluation in the center was not computed to sufficient precision\n";
	    get_out()<<"****          or the sugbradient information provided by the last oracle call was wrong.\n";
	    get_out()<<"****          Trying a function reevaluation in the center with higher precision..."<<std::endl;
	  }
	}
      }
      
      if ((!cand_minorant_is_below)||(null_step && (augval_lb>center_ub+center_gs_val))){  
	recomp++;
	sumrecomp++;
	  
	Real try_ub=center_ub;
	Real try_relprec=.001*center_relprec;
	status1=model->transform()->recompute_center(center_fid,try_ub,center_id,center_y,true,try_relprec);
	cntobjeval++;
	if (status1) {
	  if (cb_out()) get_out()<<"**** WARNING: BundleSolver::solve(): reevaluation returned error code: "<<status1<<std::endl;
	}
	
	if (try_ub>center_ub){
	  center_ub=try_ub;
	  center_objval=center_ub+cand_gs_val;
	  center_relprec=try_relprec;
	  if (cb_out()) get_out()<<"****          reevaluation in center yields higher center_objval="<<center_objval<<std::endl;         
	}
	else {
	  if (cb_out()) get_out()<<"****          reevaluation in center yields no better value"<<std::endl;
	}
      }
    }
    
        
  }while(--maxsteps);

  /*
  if (descent_step){
    int err=0;
    if (model)
      err=model->transform()->synchronize_ids(center_fid,0,center_id,cand_fid,0,cand_id,model_aggregate_id);
    center_id=0;
    cand_id=0;
    if (err) {
      if (cb_out()) 
	get_out()<<"**** ERROR: BundleSolver::solve(): synchronize_ids returned error code "<<err<<std::endl;
      return err;
    }
  }
  if (null_step) {
    int err=0;
    if (model)
      err=model->transform()->synchronize_ids(center_fid,0,center_id,cand_fid,1,cand_id,model_aggregate_id);
    center_id=0;
    cand_id=1;
    if (err) {
      if (cb_out()) 
	get_out()<<"**** ERROR: BundleSolver::solve(): synchronize_ids returned error code "<<err<<std::endl;
      return err;
    }
  }
  */

  if (cb_out(10)){
    get_out()<<"\n leaving  BundleSolver::solve():"<<std::endl;
  }    

 return 0;
} 


// *****************************************************************************
//                                 print_line_summary
// *****************************************************************************

// outputs most important data within one line

std::ostream& BundleSolver::print_line_summary(std::ostream& o) const
{
  o.setf(std::ios::showpoint);
  o.fill(' ');
  if (clockp) o<<*clockp;
  if (!terminate){
    if (null_step)
      o<<" enditn ";
    else
      o<<" enditd ";
  }
  else 
    o<<" endit_ ";
  o.width(3);o<<descent_steps;o.precision(3);
  o<<" ";o.width(4);o<<suminnerit;
  o<<" ";o.width(4);o<<sumupdatecnt;
  o<<" ";o.width(6);o<<weightu;
  o<<" ";o.width(7);o<<sqrt(aggr_dnormsqr);
  o<<" ";o.precision(10);o.width(12);o<<modelval;
  o<<" ";o.width(12);o<<center_ub+center_gs_val;
  if ((!terminate)&& null_step){
    o<<" ";o.precision(12);o.width(14);o<<augval_lb;
    o<<" ";o.precision(2);o.width(6);o<<(center_ub+center_gs_val-cand_ub-cand_gs_val)/(center_ub+center_gs_val-modelval);
  }
  o<<std::endl;
  return o;
}


// *****************************************************************************
//                                 apply_modification
// *****************************************************************************

int BundleSolver::apply_modification(const GroundsetModification& gsmdf,
				     const FunObjModMap& funmdfmap)
{
  might_be_modified=true;
  if ((typeid(gsmdf)==typeid(LPGroundsetModification))&&
      (typeid(*groundset)==typeid(Groundset))&&
      (!external_groundset)){
    LPGroundset* lpgs=new LPGroundset;
    if (lpgs==0){
      if (cb_out()){
	get_out()<<"**** ERROR: BundleSolver::apply_modification: new LPGroundset failed, no changes excecuted"<<std::endl;
      }
      return 1;
    }
    lpgs->set_cbout(this);
    lpgs->clear(groundset->get_dim(),groundset->get_groundset_id());
    lpgs->set_use_yfixing(groundset->get_use_yfixing());
    delete groundset;
    groundset=lpgs;
  }

  int retval=groundset->apply_modification(gsmdf);
  if (retval){
    if (cb_out())
      get_out()<<"**** ERROR: BundleSolver::apply_modification: groundset.apply_modification() failed and returned "<<retval<<std::endl;
    return retval;
  }

  bool no_changes=true;
  bool no_center=true;
  Integer old_center_id=center_id;
  Matrix old_center=center_y;
  if (center_id>=0){
    if (!gsmdf.no_modification()){
      gsmdf.apply_to_vars(center_y);
      center_gid=-1;
      if (groundset->is_feasible(center_gid,center_y)){
	center_id=++point_id;
	no_center=false;
      }
    }
    else {
      no_center=false;
    }
    if ((!no_center)&&(model)){
      retval=model->transform()->apply_modification(no_changes,gsmdf,funmdfmap,center_id,center_y,old_center_id,old_center);
    }
  }
  if ((no_center)&&(model)){
    retval=model->transform()->apply_modification(no_changes,gsmdf,funmdfmap,-1,center_y,-1,old_center);
    if (retval){
      if (cb_out())
	get_out()<<"**** ERROR: BundleSolver::apply_modification: model->transform()->apply_modification() failed and returned "<<retval<<std::endl;
      return retval;
    }
  }

  retval=Hp->apply_modification(gsmdf);
  if (retval){
    if (cb_out())
      get_out()<<"**** ERROR: BundleSolver::apply_modification: Hp->apply_modification() failed and returned "<<retval<<std::endl;
    return retval;
  }

  retval=bundleweight->apply_modification(gsmdf);
  if (retval){
    if (cb_out())
      get_out()<<"**** ERROR: BundleSolver::apply_modification: bundleweight->apply_modification() failed and returned "<<retval<<std::endl;
    return retval;
  }

  //separate scaling for new indices currently not supported
  if (false&&Hp->employ_variable_metric()&&(gsmdf.new_var_indices())&& null_step) {
    //for descent steps the diagonal is fully reinitialized in the loop
    //so this only deals with new entries added after null steps
    retval= variable_metric(gsmdf.new_var_indices());
    if (retval){
      if (cb_out())
	get_out()<<"**** ERROR: BundleSolver::apply_modification: variable_metric failed and returned "<<retval<<std::endl;
      return 1;
      
    }
  }

  if ((!gsmdf.no_modification())||(!no_changes)||(no_center)){
    initialize_model=true;
  }

  return 0;
}
// *****************************************************************************
//                              print_statistics
// *****************************************************************************

std::ostream& BundleSolver::print_statistics(std::ostream& out) const
{
  out<<" QPcoeff "<<QPcoeff_time;
  out<<" QPsolve "<<QPsolve_time;
  out<<" make_aggr "<<make_aggr_time;
  qp_solver->QPprint_statistics(out);
  return out;
}
  
// *****************************************************************************
//                              qp_mfile_data
// *****************************************************************************

int BundleSolver::qp_mfile_data(const CH_Matrix_Classes::Matrix& /* center_y */,
				const BundleProxObject* /* Hp */ ,
				const MinorantPointer& /* gs_subg */,
				const CH_Matrix_Classes::Symmatrix& /* Q */,
				const CH_Matrix_Classes::Matrix& /* c */,
				CH_Matrix_Classes::Real /* offset */,
				const CH_Matrix_Classes::Indexmatrix& /* yfixed */
				) const
{
  /*   //currently not in use
  UQPSolver* uqp_solver=dynamic_cast<UQPSolver*>(qp_solver);
  if ((uqp_solver)&&(uqp_solver->get_qp_save_cnt()>qp_save_cnt)){
    qp_save_cnt++;
    std::stringstream probname;
    probname<<"qp_qsdp_prob_"<<qp_save_cnt<<".m";
    std::ofstream fout((probname.str()).c_str());
    fout<<"% file "<<probname.str()<<"\n";
    
    fout<<"\n% computed cost coefficients of the QP\n";
    
    fout.precision(16);
    fout<<"\n Q=[";
    for (Integer i=0;i<Q.rowdim();i++){
      for (Integer j=0;j<Q.coldim();j++){
	fout<<" "<<Q(i,j);
      }
      if (i<Q.rowdim()-1)
	fout<<"\n";
    }
    fout<<"];";
    fout<<"\n c=[";
    for (Integer i=0;i<c.dim();i++){
	fout<<" "<<c(i);
    }
    fout<<"]';";
    fout<<"\n offset="<<offset<<";\n";

    fout<<"\n% this results from the following model information\n";

    fout<<"\n% current center yhat and its fixed indices\n\n";   
    fout<<"yhat=[";
    for (Integer i=0;i<center_y.dim();i++){
	fout<<" "<<center_y(i);
    }
    fout<<"]';\n";
    fout<<"yfixed=[";
    for (Integer i=0;i<yfixed.dim();i++){
	fout<<" "<<yfixed(i);
    }
    fout<<"]';\n";

				 
    fout<<"\n% scaling information (sofar either a weight or diagonal scaling)\n\n";
    Hp->mfile_data(fout);

    fout<<"\n% ground set\n";

    fout<<"\n gs_offset="<<gs_subg.offset()<<";";
    fout<<"\n gs_subg=[";
    for (Integer i=0;i<center_y.dim();i++){
	fout<<" "<<gs_subg.coeff(i);
    }
    fout<<"]';\n";

    //fout<<"\n% model\n";
    //output_bundle_data(fout); 

    fout<<"\n% some consistency checks\n";
    fout<<" yused=-yfixed+1;\n";
    fout<<" smat=subgmat.*(yused*ones(1,xdim));";
    fout<<" norm(c-(costs+smat'*(yhat-rhs-gs_subg)))\n";
    fout<<" norm(Q-smat'*smat)\n";

    fout.close(); 
    qp_save_cnt=uqp_solver->get_qp_save_cnt();
  }
  */
  return 0;
}

}
