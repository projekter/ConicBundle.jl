/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/QPSolver.cxx
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



#include <iomanip>
#include <sstream>
#include <fstream>
#include "QPSolver.hxx"
#include "LPGroundsetModification.hxx"
#include "BundleIdProx.hxx"
#include "QPDirectKKTSolver.hxx"

using namespace CH_Matrix_Classes;

namespace ConicBundle {
  
// *************************************************************************
//                       QPSolver::determine_indices
// *************************************************************************
  
int QPSolver::determine_indices(QPProblemData& qpd)
{
  assert(qpd.lby.dim()==qpd.dim);
  assert(qpd.uby.dim()==qpd.dim);
  qpd.lbindex.newsize(qpd.dim,1);
  qpd.lbindex.init(0,1,Integer(0));
  qpd.ubindex.newsize(qpd.dim,1);
  qpd.ubindex.init(0,1,Integer(0));
  int err=0;
  for (Integer i=0;i<qpd.dim;i++){
    if (qpd.lby(i)>CB_minus_infinity)
      qpd.lbindex.concat_below(i);
    if (qpd.uby(i)<CB_plus_infinity)
      qpd.ubindex.concat_below(i);
    if (qpd.lby(i)>qpd.uby(i)){
       if (cb_out())
	 get_out()<<"**** ERROR: QPSolver::determine_indices(): for index="<<i<<" lower bound="<<qpd.lby(i)<<" exceeds upper bound="<<qpd.uby(i)<<std::endl;
       err++;
    }     
  }

  assert(qpd.rhslb.dim()==qpd.A.rowdim());
  assert(qpd.rhsub.dim()==qpd.A.rowdim());
  qpd.rhslbindex.newsize(qpd.A.rowdim(),1);
  qpd.rhslbindex.init(0,1,Integer(0));
  qpd.rhsubindex.newsize(qpd.A.rowdim(),1);
  qpd.rhsubindex.init(0,1,Integer(0));
  qpd.rhseqindex.newsize(qpd.A.rowdim(),1);
  qpd.rhseqindex.init(0,1,Integer(0));
  for (Integer i=0;i<qpd.A.rowdim();i++){
    if ((qpd.rhsub(i)<CB_plus_infinity)&&(fabs(qpd.rhsub(i)-qpd.rhslb(i))<=1e-8*max(1.,std::fabs(qpd.rhsub(i))))){          //equation
      qpd.rhslb(i)=qpd.rhsub(i);
      qpd.rhseqindex.concat_below(i);
      continue;
    }
    if (qpd.rhslb(i)>CB_minus_infinity)  // greater or equal lb(i)
      qpd.rhslbindex.concat_below(i);
    if (qpd.rhsub(i)<CB_plus_infinity)   // lower or equal ub(i)
      qpd.rhsubindex.concat_below(i);
    if (qpd.rhslb(i)>qpd.rhsub(i)){
       if (cb_out())
	 get_out()<<"**** ERROR: QPSolver::determine_indices(): for index="<<i<<" right hand side lower bound="<<qpd.lby(i)<<" exceeds right hand sid upper bound="<<qpd.uby(i)<<std::endl;
       err++;
    }     
  }
  assert(err==0);
  return err;
}
 
// *************************************************************************
//                       QPSolver::preprocess_data
// *************************************************************************
  
int QPSolver::preprocess_data(const Matrix& center_y,Indexmatrix* yfixed,bool& no_changes)
{
  int err=0;
  no_changes=true;

  if (preproc_fixed.coldim()==0){
    preproc_fixed.init(original_data.dim,1,Integer(0));
    preproc_fixedval.init(original_data.dim,1,0.);
  }
  
  if ((preproc_indices.coldim()==0)||(yfixed)){
    preproc_indices.newsize(original_data.dim,1);chk_set_init(preproc_indices,1);
    //find the fixed indices
    //std::cout<<" fix["; //TEST
    Integer cnt=0;
    for (Integer i=0;i<original_data.dim;i++){
      if ((preproc_fixed(i))||((yfixed)&&((*yfixed)(i)>0))||
	  (original_data.uby(i)-original_data.lby(i)<1e-10*max(1.,std::fabs(original_data.uby(i))))){
	if ((preproc_fixed(i)==0)||(std::fabs(preproc_fixedval(i)-center_y(i))>=eps_Real*max(std::fabs(center_y(i)),1.)))
	  no_changes=false;
	preproc_fixedval(i)=center_y(i);
	preproc_fixed(i)=1;
	//std::cout<<" "<<i<<"("<<center_y(i)<<")"; //TEST
	if (yfixed)
	  (*yfixed)(i)=1;
      }
      else {
	preproc_indices(cnt++)=i;
      }
    }
    //std::cout<<"]"<<std::endl; //TEST
    preproc_indices.reduce_length(cnt);
  }

  if (!no_changes){
    //store the original data and generate the data with fixed parts eliminated
    preproc_data.dim=preproc_indices.rowdim();
    preproc_data.lby=original_data.lby(preproc_indices);
    preproc_data.uby=original_data.uby(preproc_indices);
    preproc_data.rhslb=original_data.rhslb;
    preproc_data.rhsub=original_data.rhsub;
    Matrix tmprhs;
    genmult(original_data.A,preproc_fixedval,tmprhs);
    for(Integer i=0;i<preproc_data.rhslb.rowdim();i++){
      if (preproc_data.rhslb(i)>CB_minus_infinity){
	preproc_data.rhslb(i)=max(CB_minus_infinity,preproc_data.rhslb(i)-tmprhs(i));
      }
      if (preproc_data.rhsub(i)<CB_plus_infinity){
	preproc_data.rhsub(i)=min(CB_plus_infinity,preproc_data.rhsub(i)-tmprhs(i));
      }
    }
    preproc_data.A=original_data.A.cols(preproc_indices);
    err=determine_indices(preproc_data);
    if ((err)&&(cb_out())){
      get_out()<<"**** ERROR in QPSolver::preprocess_data(...): determine_indices failed and returned "<<err<<std::endl;
    }
    preproc_data.groundset_gamma=original_data.groundset_gamma+ip(original_data.groundset_c,preproc_fixedval);
    preproc_data.groundset_c=original_data.groundset_c(preproc_indices);

    // add an AFT to the block
    Indexmatrix indj(Range(0,preproc_data.dim-1));
    Matrix val(preproc_data.dim,1,1.);
    preproc_aft.init(1.,0.,0,new Matrix(preproc_fixedval),new Sparsemat(original_data.dim,preproc_data.dim,preproc_data.dim,preproc_indices,indj,val));
    preproc_bundle_projection.clear();

    // //TEST begin
    // {
    //   const MinorantBundle& bun=get_model_data_ptr()->get_bundle();
    //   std::cout<<" bundle_before=\n";
    //   for(unsigned i=0;i<bun.size();i++)
    // 	bun[i].display(std::cout,6);
    // }
    // //TEST end
    
  }

  if (preproc_indices.rowdim()==original_data.dim)
    qp_data=&original_data;
  else {
    qp_data=&preproc_data;
    get_model_data_ptr()->push_aft(&preproc_aft,0,0,&preproc_bundle_projection);
  }
  
  // //TEST begin
  // {
  //   std::cout<<" bundle_after=\n";
  //   const MinorantBundle& bun=get_model_data_ptr()->get_bundle();
  //   for(unsigned i=0;i<bun.size();i++)
  //     bun[i].display(std::cout,6);
  // }
  // //TEST end

  return err;
}


// *************************************************************************
//                       QPSolver::postprocess_data
// *************************************************************************

int QPSolver::postprocess_data(bool round_to_active_bounds)
{
  //-----get the solution
  sol_val_lb=QPget_dualval();
  sol_val_ub=QPget_primalval();
  gs_aggr_gradient.init(qp_data->groundset_c);
  QPadd_Aty(QPget_y(),gs_aggr_gradient);
  gs_aggr_gradient+=QPget_zub();
  gs_aggr_gradient-=QPget_zlb();
  gs_aggr_offset=qp_data->groundset_gamma;
  gs_aggr_offset+=ip(qp_data->lby,QPget_zlb())-ip(qp_data->uby,QPget_zub());
  gs_aggr_offset+=-ip(qp_data->rhsub,QPget_rhszlb())+ip(qp_data->rhslb,QPget_rhszub());
  for (Integer i=0;i<qp_data->rhslb.rowdim();i++){
    if (qp_data->rhslb(i)==qp_data->rhsub(i)){
      gs_aggr_offset-=QPget_y()(i)*qp_data->rhslb(i);
    }
  }
  
  if (round_to_active_bounds){
    Indexmatrix x_act;
    QPget_x(sol_point,x_act);
    //check whether this needs corrections on the constraint side
    if ((qp_data->A.rowdim()>0)&&(sum(x_act)<x_act.rowdim())){
      Matrix sol_s;
      Indexmatrix s_act;
      QPget_s(sol_s,s_act);
      Matrix ds;
      genmult(qp_data->A,sol_point,ds,1.,0.);
      //check violation of equality contsraints
      Real viol=0.;
      bool ineqviol=false;
      for (Integer i=0;i<ds.rowdim();i++){
	if (s_act(i)==0){ //should hold as equation
	  viol+=sqr(sol_s(i)+ds(i));
	}
	else {// should be within bounds
	  if((qp_data->rhslb(i)>ds(i))||(qp_data->rhsub(i)<ds(i))){
	    ineqviol=true;
	    break;
	  }
	}
      }

      if ((!ineqviol)&&(std::sqrt(viol)>QPget_parameters()->QPget_primal_infeasibility_eps())){
	//violation of equalitiy constraints too large, try least squares correction
	//form right hand side
	ds+=sol_s;
	//form matrix
	Indexmatrix dxind(x_act.rowdim(),1);
	dxind.init(0,0,Integer(0));
	for(Integer i=0;i<x_act.rowdim();i++){
	  if (x_act(i)){
	    dxind.concat_below(i);
	  }
	}
	Sparsemat SparseAt(qp_data->A.cols(dxind));
	SparseAt.transpose();
	Matrix At(dxind.rowdim()+sum(s_act),ds.rowdim());
        At.init(SparseAt);
	At.enlarge_below(sum(s_act),0.);
	Indexmatrix dsind(s_act.rowdim(),1);
	dsind.init(0,0,Integer(0));
	for(Integer i=0;i<s_act.rowdim();i++){
	  if (s_act(i)){
	    At(dxind.rowdim()+dsind.rowdim(),i)=1.;
	    dsind.concat_below(i);
	  }
	}
	//solve by QR transposed
	Indexmatrix piv;
	Integer r=At.QR_factor(piv);
	if (r==At.coldim()){
	  Matrix dx=ds(piv);
	  Matrix L=At.rows(Range(0,r-1));
	  L.triu();
	  L.transpose();
	  L.tril_solve(ds,eps_Real);
	  dx.enlarge_below(At.rowdim()-r);
	  At.Q_times(dx,r);
	  for(Integer i=0;i<dxind.rowdim();i++){
	    Integer ii=dxind(i);
	    Real& sx=sol_point(ii);
	    sx+=dx(i);
	    if ((sx<qp_data->lby(ii))||(sx>qp_data->uby(ii)))
	      ineqviol=true;
	  }
	  for(Integer i=0;i<dsind.rowdim();i++){
	    Integer ii=dsind(i);
	    Real& sx=sol_s(ii);
	    sx+=dx(dxind.rowdim()+i);
	    if ((sx<qp_data->rhslb(ii))||(sx>qp_data->rhsub(ii)))
	      ineqviol=true;
	  }
	  assert(norm2(qp_data->A*sol_point+sol_s)<1e-6*norm2(ds));
	}
	else { //rank difficulties, give up correcting
	  ineqviol=true;
	}
      } //endif constraints need corrections
      if (ineqviol){
	//rounding did not produce a feasible solution
	round_to_active_bounds=false;
      }
    } //endif there are contsraints 
  }
      
  if (!round_to_active_bounds){
    sol_point=QPget_x();
  }

  //if necessary, insert fixed coordinates and correct offset
  if (qp_data==&preproc_data){
    Matrix tmpvec=sol_point;
    sol_point=preproc_fixedval;
    sol_point.subassign(preproc_indices,tmpvec);
    tmpvec=gs_aggr_gradient;
    gs_aggr_gradient=original_data.groundset_c;
    gs_aggr_offset-=ip(gs_aggr_gradient,preproc_fixedval);
    gs_aggr_gradient.subassign(preproc_indices,tmpvec);

    get_model_data_ptr()->pop_aft();
    qp_data=&original_data;
  }
  
  return 0;
}

  // *************************************************************************
//                             QPSolver::clear
// *************************************************************************

void QPSolver::QPclear()
{  
  original_data.dim=0;
  original_data.lby.init(original_data.dim,1,CB_minus_infinity);              
  original_data.uby.init(original_data.dim,1,CB_plus_infinity);             
  original_data.lbindex.init(0,1,Integer(0));  
  original_data.ubindex.init(0,1,Integer(0));
  original_data.groundset_c.init(original_data.dim,1,0.);
  original_data.groundset_gamma=0.;
  original_data.A.init(0,original_data.dim);
  original_data.rhslb.init(0,1,0.);
  original_data.rhsub.init(0,1,0.);
  original_data.rhslbindex.init(0,1,Integer(0));  
  original_data.rhsubindex.init(0,1,Integer(0));
  original_data.rhseqindex.init(0,1,Integer(0));
  original_data.Hp=0;
  
  delete preproc_data.Hp;
  preproc_data.Hp=0;

  preproc_fixed.init(0,0,Integer(0));
  preproc_fixedval.init(0,0,0.);
  preproc_indices.init(0,0,Integer(0));
  preproc_aft.init();
  preproc_bundle_projection.clear();

  sol_val_lb=min_Real;
  sol_val_ub=max_Real;
  sol_point.init(original_data.dim,1,0.);
  gs_aggr_offset=0.;
  gs_aggr_gradient.init(original_data.dim,1,0.);

  primal_infeasibility_eps=1e-8;  
  dual_infeasibility_eps=1e-7;
  lower_bound_gap_eps=.99;    // ensures that the lower bound increases a little bit
  upper_bound_gap_eps=.8;     //
  objective_gap_eps=1e-10;


  QPIclear();
  clock.start();
}


// *************************************************************************
//                             QPSolver::QPsolve
// *************************************************************************

int QPSolver::QPsolve(const Matrix& center_y,
		      Real lower_bound,
		      Real upper_bound,
		      Real relprec,
		      QPSolverProxObject* inHp,
		      const MinorantPointer& /* gs_aggr */,
		      Indexmatrix* yfixed)
{
  assert(get_model_data_ptr()!=0);
  assert(upper_bound>lower_bound);
  
  original_data.Hp=dynamic_cast<BundleProxObject*>(inHp);
  if (original_data.Hp==0){
    if (cb_out())
      get_out()<<"**** ERROR in QPSolver::QPsolve(.......): dynamic_cast of QPSolverProxObject* to BundleProxObject* failed, cannot continue"<<std::endl;
    return 1;
  }

  CH_Tools::Microseconds coeff_start=clock.time();

  //--- get cost matrices
  bool no_changes;
  if (preprocess_data(center_y,yfixed,no_changes)){
    if (cb_out())
      get_out()<<"**** ERROR in QPSolver::QPsolve(.......): preprocess_data failed"<<std::endl;
    return 1;
  }

  const Matrix* cyp=&center_y;
  Matrix tmpy;
  if (qp_data==&preproc_data){
    delete preproc_data.Hp;
    preproc_data.Hp=dynamic_cast<BundleProxObject*>(original_data.Hp->projected_clone(preproc_indices));
    if (preproc_data.Hp==0){
      if (cb_out())
	get_out()<<"**** ERROR in QPSolver::QPsolve(.......): dynamic_cast of prox clone to BundleProxObject* failed, cannot continue "<<std::endl;
      return 1;
    }
    tmpy=center_y(preproc_indices);
    cyp=&tmpy;
  }

  qp_data->c.init(qp_data->dim,1,0.);
  qp_data->Hp->add_Hx(*cyp,qp_data->c);
  qp_data->gamma=qp_data->groundset_gamma+.5*ip(qp_data->c,*cyp);   //H-norm squared of center_y
  qp_data->c*=-1;
  qp_data->c+=qp_data->groundset_c;
  if (!get_model_data_ptr()->get_constant_minorant().empty())
    get_model_data_ptr()->get_constant_minorant().get_minorant(qp_data->gamma,qp_data->c,0,1.,true);


  //set the stopping criteria
  QPSolverParameters* paramsp=QPget_parameters();
  paramsp->QPset_objective_gap_eps(relprec); 
  paramsp->QPset_primal_infeasibility_eps(min(1e-3,primal_infeasibility_eps*(max(1.,norm2(qp_data->A)))));
  paramsp->QPset_dual_infeasibility_eps(min(1e-3,dual_infeasibility_eps*max(1.,upper_bound-lower_bound)));
  paramsp->QPset_lower_and_upper_bounds(lower_bound,upper_bound);
  paramsp->QPset_lower_bound_gap_eps(lower_bound_gap_eps);
  paramsp->QPset_upper_bound_gap_eps(upper_bound_gap_eps);
  paramsp->QPset_maxiter(min(5*qp_data->dim,Integer(std::log(qp_data->dim)*100)));

  if (paramsp->QPget_KKTsolver()->QPinit_KKTdata(qp_data->Hp,model_block,&(qp_data->A),&(qp_data->rhseqindex))){
    if (cb_out())
      get_out()<<"**** WARNING in QPSolver::QPsolve(.......): choice of KKTsolver in QPSolverParameters is not compatible with the data of the BundleSubproblem; trying to reset to a QPDirectKKTsolver"<<std::endl;
    QPDirectKKTSolver* KKTsolver=new QPDirectKKTSolver(false,this,0);
    if (KKTsolver->QPinit_KKTdata(qp_data->Hp,model_block,&(qp_data->A),&(qp_data->rhseqindex))){
      if (cb_out())
	get_out()<<"**** ERROR in QPSolver::QPsolve(.......): attempt to initialize QPDirectKKTsolver failed as well, aborting"<<std::endl;
      return 1;
    }
    paramsp->QPset_KKTsolver(KKTsolver);
  }

  if (paramsp->QPget_use_socqp()){
    if (QPset_socqp().clear_prox(qp_data->Hp)){
      if (cb_out())
	get_out()<<"**** WARNING in QPSolver::QPsolve(.......): request for second order model for the proximal term, but the proximal term and solver are incompatible for this, using standard settings instead"<<std::endl;
      paramsp->QPset_use_socqp(false);
    }
    else {
      QPset_socqp().set_cbout(this);
    }
  }

  CH_Tools::Microseconds solve_start=clock.time();
  QPcoeff_time+=solve_start-coeff_start;

  //solve the QP
  Real skip_factor=-1.;
  int status=QPIsolve((skip_factor<0.),skip_factor);
  if (status){
    if (cb_out()){
      get_out()<<"**** WARNING in QPSolver::QPsolve(): QPIsolve returned "<<status<<std::endl;
    } 
  }

  if (postprocess_data(yfixed!=0)){
    if (cb_out()){
      get_out()<<"**** WARNING in QPSolver::QPsolve(): postprocess_data() failed "<<std::endl;
    } 
  }


  if ((QPget_dualval()-upper_bound>1e-10*(upper_bound-lower_bound))&&(cb_out())){
    get_out().precision(12);
    get_out()<<"**** WARNING in QPSolver::QPsolve(): dualval = "<<QPget_dualval();
    get_out()<<" > "<<upper_bound<<" = upper bound"<<std::endl;
    Real augval_lb;
    Real augval_ub;
    Matrix newy;
    Real gsaggrval;
    Matrix gsaggr;
    QPget_solution(augval_lb,augval_ub,newy,gsaggrval,gsaggr);
    get_out()<<" augval_lb="<<augval_lb;
    get_out()<<" augval_ub="<<augval_ub;
    get_out()<<" gsval(newy)="<<gsaggrval+ip(gsaggr,newy);
    Matrix dy(center_y);
    get_out()<<" gsval(ceny)="<<gsaggrval+ip(gsaggr,dy);
    dy-=newy;
    get_out()<<" Hnormsqr(dy)/2="<<original_data.Hp->norm_sqr(dy)/2.<<std::endl;
    get_out()<<" groundset_gamma="<<original_data.groundset_gamma;
    get_out()<<" Hnormsqr(cy)/2="<<original_data.gamma-original_data.groundset_gamma<<std::endl;
    dy.init(center_y.rowdim(),1,0.);
    original_data.Hp->add_Hx(center_y,dy);
    get_out()<<" addnorm(cy)/2="<<ip(dy,center_y)/2.;
    get_out()<<" Hnormsqr(cy)/2="<<original_data.Hp->norm_sqr(center_y)/2.<<std::endl;
  }

  original_data.Hp->set_short_QPsteps(QPget_large_predictor_cnt());

  QPsolve_time+=clock.time()-solve_start;


  // //TEST begin
  // std::cout.precision(12);
  // for (Integer i=0;i<original_data.dim;i++){
  //   if ((sol_point(i)<original_data.lby(i)-eps_Real)||
  // 	(sol_point(i)>original_data.uby(i)+eps_Real))
  //     std::cout<<"**** TEST ERROR: sol_point("<<i<<")="<<sol_point(i)<<" violates bounds ["<<original_data.lby(i)<<","<<original_data.uby(i)<<std::endl;
  // }
  // Matrix tmpvec;
  // if (original_data.A.rowdim()>0){
  //   genmult(original_data.A,sol_point,tmpvec);
  //   for (Integer i=0;i<tmpvec.rowdim();i++){
  //     if ((tmpvec(i)<original_data.rhslb(i)-1e-8*(1.+std::fabs(tmpvec(i))))||
  // 	  (tmpvec(i)>original_data.rhsub(i)+1e-8*(1.+std::fabs(tmpvec(i)))))
  // 	std::cout<<"**** TEST ERROR: (A*sol_point)("<<i<<")="<<tmpvec(i)<<" violates bounds ["<<original_data.rhslb(i)<<","<<original_data.rhsub(i)<<std::endl;
  //   }
  // }
  // Real tmpval=gs_aggr_offset;
  // tmpvec=gs_aggr_gradient;
  // if (model_block)
  //   model_block->add_modelx_aggregate(tmpval,tmpvec);
  // Real eval=tmpval+ip(sol_point,tmpvec)+original_data.Hp->norm_sqr(center_y-sol_point)/2;
  // if (std::fabs(eval-sol_val)>1e-6*(1.+std::fabs(eval))){
  //   std::cout<<"**** TEST WARNING: testeval="<<eval<<" sol_val="<<sol_val<<std::endl;
  // }
  // //TEST end
  
  return status;
}


// *************************************************************************
//                             QPSolver::QPupdate
// *************************************************************************

int QPSolver::QPupdate(const Matrix& center_y,
		       Real lower_bound,
		       Real upper_bound,
		       Real relprec,
		       QPSolverProxObject* inHp,
		       const MinorantPointer& /* gs_aggr */,
		       Indexmatrix* yfixed,
		       const MinorantPointer& delta_gs_aggr,
		       const Indexmatrix& /* delta_index */)
{
  assert(false); ///currently this does not work and does not make sesne, because the groundset_aggregate is fully included and cannot be added or changed from outside
  
  assert(get_model_data_ptr()!=0);
  assert(original_data.Hp);

  BundleProxObject* Hp2=dynamic_cast<BundleProxObject*>(inHp);
  if (Hp2!=original_data.Hp){
    if (cb_out())
      get_out()<<"**** ERROR in QPSolver::QPupdate(.........): dynamic_cast of QPSolverProxObject* to BundleProxObject* produced a different ProxObject than in QPsolve(), cannot continue"<<std::endl;
    return 1;
  }

  CH_Tools::Microseconds coeff_start=clock.time();

  //--- get cost matrices

  bool no_changes;
  if (preprocess_data(center_y,yfixed,no_changes)){
    if (cb_out())
      get_out()<<"**** ERROR in QPSolver::QPupdate(.......): preprocess_data failed"<<std::endl;
    return 1;
  }

  delta_gs_aggr.get_minorant(original_data.gamma,original_data.c,0,1.,true);
  if (qp_data==&preproc_data){
    Real val;
    Matrix grad;
    delta_gs_aggr.get_minorant(val,grad,0,1.,true);
    preproc_data.gamma+=ip(grad,preproc_fixedval);
    preproc_data.c+=grad(preproc_indices);

    if ((!no_changes)&&(norm2(grad)>eps_Real)){
      if (cb_out()){
	get_out()<<"**** WARNING in QPSolve::QPupdated(): updates with changes in the fixed indices will not work, because the external gs_aggregate is not accumulated explicitly "<<std::endl;
      } 
    }
  }

  if (!no_changes){
    //newly fixed indices, rebuild cost matrix from scratch 
    const Matrix* cyp=&center_y;
    Matrix tmpy;
    if (qp_data==&preproc_data){
      delete preproc_data.Hp;
      preproc_data.Hp=dynamic_cast<BundleProxObject*>(original_data.Hp->projected_clone(preproc_indices));
      if (preproc_data.Hp==0){
	if (cb_out())
	  get_out()<<"**** ERROR in QPSolver::QPsolve(.......): dynamic_cast of prox clone to BundleProxObject* failed, cannot continue "<<std::endl;
	return 1;
      }
      tmpy=center_y(preproc_indices);
      cyp=&tmpy;
    }
    
    qp_data->c.init(qp_data->dim,1,0.);
    qp_data->Hp->add_Hx(*cyp,qp_data->c);
    qp_data->gamma=qp_data->groundset_gamma+.5*ip(qp_data->c,*cyp);   //H-norm squared of center_y
    qp_data->c*=-1;
    qp_data->c+=qp_data->groundset_c;
    if (!get_model_data_ptr()->get_constant_minorant().empty())
      get_model_data_ptr()->get_constant_minorant().get_minorant(qp_data->gamma,qp_data->c,0,1.,true);
  }
        
  //set the stopping criteria
  QPSolverParameters* paramsp=QPget_parameters();
  paramsp->QPset_objective_gap_eps(relprec); 
  paramsp->QPset_primal_infeasibility_eps(min(1e-3,primal_infeasibility_eps*(max(1.,norm2(qp_data->A)))));
  paramsp->QPset_dual_infeasibility_eps(dual_infeasibility_eps*max(1.,upper_bound-lower_bound));
  paramsp->QPset_lower_and_upper_bounds(lower_bound,upper_bound);
  paramsp->QPset_lower_bound_gap_eps(lower_bound_gap_eps);
  paramsp->QPset_upper_bound_gap_eps(upper_bound_gap_eps);
  paramsp->QPset_maxiter(min(5*qp_data->dim,Integer(std::log(qp_data->dim)*100)));

  if ((!no_changes)&&(paramsp->QPget_KKTsolver()->QPinit_KKTdata(qp_data->Hp,model_block,&(qp_data->A),&(qp_data->rhseqindex)))){
    if (cb_out())
      get_out()<<"**** WARNING in QPSolver::QPupdate(.......): choice of KKTsolver in QPSolverParameters is not compatible with the data of the BundleSubproblem; trying to reset to a QPDirectKKTsolver"<<std::endl;
    QPDirectKKTSolver* KKTsolver=new QPDirectKKTSolver(false,this,0);
    if (KKTsolver->QPinit_KKTdata(qp_data->Hp,model_block,&(qp_data->A),&(qp_data->rhseqindex))){
      if (cb_out())
	get_out()<<"**** ERROR in QPSolver::QPsolve(.......): attempt to initialize QPDirectKKTsolver failed as well, aborting"<<std::endl;
      return 1;
    }
    paramsp->QPset_KKTsolver(KKTsolver);
  }

  CH_Tools::Microseconds solve_start=clock.time();
  QPcoeff_time+=solve_start-coeff_start;

  //solve the QP
  Real skip_factor=-1.;
  int status=QPIsolve((skip_factor<0.),skip_factor);
  if (status){
    if (cb_out()){
      get_out()<<"**** WARNING in QPSolve::QPupdate(): QPIsolve returned "<<status<<std::endl;
    } 
  }

  if (postprocess_data(yfixed!=0)){
    if (cb_out()){
      get_out()<<"**** WARNING in QPSolver::QPupdate(): postprocess_data() failed "<<std::endl;
    } 
  }

  original_data.Hp->set_short_QPsteps(QPget_large_predictor_cnt());


  QPsolve_time+=clock.time()-solve_start;

  return status;
}

 
// *************************************************************************
//                             QPSolver::QPresolve
// *************************************************************************


int QPSolver::QPresolve(Real lower_bound,
			Real upper_bound,
			Real relprec)
{
  assert(get_model_data_ptr()!=0);
  assert(original_data.Hp);

  CH_Tools::Microseconds coeff_start=clock.time();
 
  //--- get cost matrices
  if (preproc_indices.rowdim()==original_data.dim)
    qp_data=&original_data;
  else {
    qp_data=&preproc_data;
    get_model_data_ptr()->push_aft(&preproc_aft,0,0,&preproc_bundle_projection);
  }
  
  //set the stopping criteria
  QPSolverParameters* paramsp=QPget_parameters();
  paramsp->QPset_objective_gap_eps(relprec); 
  paramsp->QPset_primal_infeasibility_eps(min(1e-3,primal_infeasibility_eps*(max(1.,norm2(qp_data->A)))));
  paramsp->QPset_dual_infeasibility_eps(min(1e-3,dual_infeasibility_eps*max(1.,upper_bound-lower_bound)));
  paramsp->QPset_lower_and_upper_bounds(lower_bound,upper_bound);
  paramsp->QPset_lower_bound_gap_eps(lower_bound_gap_eps);
  paramsp->QPset_upper_bound_gap_eps(upper_bound_gap_eps);
  paramsp->QPset_maxiter(min(5*qp_data->dim,Integer(std::log(qp_data->dim)*100)));

  CH_Tools::Microseconds solve_start=clock.time();
  QPcoeff_time+=solve_start-coeff_start;

  //solve the QP
  Real skip_factor=-1.;
  int status=QPIsolve((skip_factor<0.),skip_factor);
  if (status){
    if (cb_out()){
      get_out()<<"**** WARNING in QPSolve::QPsolver(): QPIsolve returned "<<status<<std::endl;
    } 
  }

  if (postprocess_data(false)){
    if (cb_out()){
      get_out()<<"**** WARNING in QPSolver::QPupdate(): postprocess_data() failed "<<std::endl;
    } 
  }

  original_data.Hp->set_short_QPsteps(QPget_large_predictor_cnt());

  QPsolve_time+=clock.time()-solve_start;

  return status;
}

// *************************************************************************
//                             QPSolver::QPget_solution
// *************************************************************************

int QPSolver::QPget_solution(CH_Matrix_Classes::Real& augval_lb,
			     CH_Matrix_Classes::Real& augval_ub,
			     CH_Matrix_Classes::Matrix& new_point,
			     CH_Matrix_Classes::Real& gsaggr_offset,
			     CH_Matrix_Classes::Matrix& gsaggr_gradient)
{
  new_point=sol_point;
  augval_lb=sol_val_lb;
  augval_ub=sol_val_ub;
  gsaggr_offset=gs_aggr_offset;
  gsaggr_gradient=gs_aggr_gradient;
  
  return 0;
}


// *************************************************************************
//                             QPSolver::is_feasible
// *************************************************************************

bool QPSolver::QPis_feasible(const Matrix &y,Real relprec)
{
  //check dimension
  if (y.dim()!=original_data.dim){
    return false;
  }
  
  //check bounds
  if (original_data.lbindex.dim()>0){
    for(Integer i=0;i<original_data.lbindex.dim();i++){
      Integer ind=original_data.lbindex(i);
      if (y(ind)<original_data.lby(ind)){
	return false;
      }
    }
  }
  if (original_data.ubindex.dim()>0){
    for(Integer i=0;i<original_data.ubindex.dim();i++){
      Integer ind=original_data.ubindex(i);
      if (y(ind)>original_data.uby(ind)){
	return false;
      }
    }
  }
  
  //check constraints
  if (original_data.A.rowdim()!=0){
    Matrix tmp;
    genmult(original_data.A,y,tmp);
    for (Integer i=0;i<tmp.rowdim();i++){
      if (original_data.rhslb(i)==original_data.rhsub(i)){
	if (fabs(tmp(i)-original_data.rhslb(i))>100.*relprec*(fabs(original_data.rhslb(i))+1.)){
	  return false;
	}
      }
      else {
	if ((original_data.rhslb(i)!=CB_minus_infinity)&&(tmp(i)<original_data.rhslb(i)-relprec*(fabs(original_data.rhslb(i))+1.)))
	  return false;
	if ((original_data.rhsub(i)!=CB_plus_infinity)&&(tmp(i)>original_data.rhsub(i)+relprec*(fabs(original_data.rhsub(i))+1.)))
	  return false;
      }
    }
  }

  return true;
}
 
// *****************************************************************************
//                             ensure_feasibility
// *****************************************************************************

int QPSolver::QPensure_feasibility(Matrix& y,
				   bool& ychanged,
				   QPSolverProxObject* inHp,
				   Real relprec)
{

  original_data.Hp=dynamic_cast<BundleProxObject*>(inHp);
  if (original_data.Hp==0){
    if (cb_out())
      get_out()<<"**** ERROR in QPSolver::QPensure_feasibility(....): dynamic_cast of QPSolverProxObject* to BundleProxObject* failed, cannot continue"<<std::endl;
    return 1;
  }

  //check dimension
  if (y.dim()!=original_data.dim){
    y.init(original_data.lby.rowdim(),1,0.);
    ychanged=true;
  }
  
  //check and enforce lower bounds
  if (original_data.lbindex.dim()>0){
    for(Integer i=0;i<original_data.lbindex.dim();i++){
      Integer ind=original_data.lbindex(i);
      if (y(ind)<original_data.lby(ind)){
	ychanged=true;
	y(ind)=original_data.lby(ind);
      }
    }
  }
  if (original_data.ubindex.dim()>0){
    for(Integer i=0;i<original_data.ubindex.dim();i++){
      Integer ind=original_data.ubindex(i);
      if (y(ind)>original_data.uby(ind)){
	ychanged=true;
	y(ind)=original_data.uby(ind);
      }
    }
  }
  
  //check equality constraints
  bool violated=false;
  if (original_data.A.rowdim()!=0){
    Matrix tmp;
    genmult(original_data.A,y,tmp);
    for (Integer i=0;i<tmp.rowdim();i++){
      if (original_data.rhslb(i)==original_data.rhsub(i)){
	if (fabs(tmp(i)-original_data.rhslb(i))>relprec*(fabs(original_data.rhslb(i))+1.)){
	  violated=true;
	  break;
	}
      }
      else {
	if (((original_data.rhslb(i)!=CB_minus_infinity)&&(tmp(i)<original_data.rhslb(i)-relprec*(fabs(original_data.rhslb(i))+1.)))||
	    ((original_data.rhsub(i)!=CB_plus_infinity)&&(tmp(i)>original_data.rhsub(i)+relprec*(fabs(original_data.rhsub(i))+1.)))){
	  violated=true;
	  break;
	}
      }
    }
  }
  if (!violated){
    return 0;
  }
    
  //set the cost coefficients so as to minimize the distance to y
  original_data.c.init(original_data.dim,1,0.);
  original_data.Hp->add_Hx(y,original_data.c);
  original_data.gamma=ip(y,original_data.c)/2.;
  original_data.c*=-1;

  qp_data=&original_data;

  //set the stopping criteria
  QPSolverParameters* paramsp=QPget_parameters();
  paramsp->QPset_objective_gap_eps(1.e10);  //no high precision required
  paramsp->QPset_primal_infeasibility_eps(1e-8*(1.+norm2(original_data.A)));
  paramsp->QPset_dual_infeasibility_eps(1e10);
  paramsp->QPset_lower_and_upper_bounds(-1e10,max_Real);
  paramsp->QPset_lower_bound_gap_eps(1e10);
  paramsp->QPset_upper_bound_gap_eps(1e10);

  if (paramsp->QPget_KKTsolver()->QPinit_KKTdata(original_data.Hp,model_block,&original_data.A,&original_data.rhseqindex)){
    if (cb_out())
      get_out()<<"**** WARNING in QPSolver::QPsolve(.......): choice of KKTsolver in QPSolverParameters is not compatible with the data of the BundleSubproblem; trying to reset to a QPDirectKKTsolver"<<std::endl;
    QPDirectKKTSolver* KKTsolver=new QPDirectKKTSolver(false,this,0);
    if (KKTsolver->QPinit_KKTdata(original_data.Hp,model_block,&original_data.A,&original_data.rhseqindex)){
      if (cb_out())
	get_out()<<"**** ERROR in QPSolver::QPsolve(.......): attempt to initialize QPDirectKKTsolver failed as well, aborting"<<std::endl;
      return 1;
    }
    paramsp->QPset_KKTsolver(new QPDirectKKTSolver);
  }
  
  if (paramsp->QPget_use_socqp()){
    if (QPset_socqp().clear_prox(qp_data->Hp)){
      if (cb_out())
	get_out()<<"**** WARNING in QPSolver::QPsolve(.......): request for second order model for the proximal term, but the proximal term and solver are incompatible for this, using standard settings instead"<<std::endl;
      paramsp->QPset_use_socqp(false);
    }
    else {
      QPset_socqp().set_cbout(this);
    }
  }


  //solve the QP
  int status=QPIsolve();
  if (status){
    if (cb_out()){
      get_out()<<"**** WARNING in QPSolver::ensure_feasibility: QPIsolve() returned "<<status<<std::endl;
    } 
  }
  
  //retrieve the solution
  y=QPget_x();
  ychanged=true;
  
  return status;
}

// *************************************************************************
//                             QPSolver::QPprefer_UQPSolver
// *************************************************************************

bool QPSolver::QPprefer_UQPSolver(QPSolverProxObject* prox) const
{
  BundleProxObject* inHp=dynamic_cast<BundleProxObject*>(prox);
  if ((inHp!=0)&&(QPget_parameters()->QPallow_UQPSolver())){
    Matrix diagQ;
    const Matrix* Vp;
    inHp->get_precond(diagQ,Vp);
    if (
	((!QPconstrained())&&(inHp->is_DLR()))
	||
	((original_data.A.rowdim()==0)&&(inHp->is_DLR())&&((Vp==0)||(Vp->coldim()==0)))
	){
      return true;
    }
  }
  return false;
}

// *************************************************************************
//                             QPSolver::solve
// *************************************************************************

//solve the quadratic problem for the given cost function and precision

int QPSolver::solve(BundleProxObject* inHp,
		    const Matrix& inc,
		    Real ingamma,
		    Real lowerbound,
		    Real upperbound,
		    Real relprec,
		    Real skip_factor)
{
  assert(inHp);
  assert(original_data.dim==inc.rowdim());
  assert(1==inc.coldim());

  original_data.Hp=inHp;
  original_data.c=inc;
  original_data.gamma=ingamma;
  
  qp_data=&original_data;

  //set the stopping criteria
  QPSolverParameters* paramsp=QPget_parameters();
  paramsp->QPset_objective_gap_eps(relprec); 
  paramsp->QPset_primal_infeasibility_eps(primal_infeasibility_eps*(max(1.,norm2(original_data.A))));
  paramsp->QPset_dual_infeasibility_eps(dual_infeasibility_eps*max(1.,upperbound-lowerbound));
  paramsp->QPset_lower_and_upper_bounds(lowerbound,upperbound);
  paramsp->QPset_lower_bound_gap_eps(lower_bound_gap_eps);
  paramsp->QPset_upper_bound_gap_eps(upper_bound_gap_eps);

  if (paramsp->QPget_KKTsolver()->QPinit_KKTdata(original_data.Hp,model_block,&original_data.A,&original_data.rhseqindex)){
    if (cb_out())
      get_out()<<"**** WARNING in QPSolver::solve(.......): choice of KKTsolver in QPSolverParameters is not compatible with the data of the BundleSubproblem; trying to reset to a QPDirectKKTsolver"<<std::endl;
    QPDirectKKTSolver* KKTsolver=new QPDirectKKTSolver(false,this,0);
    if (KKTsolver->QPinit_KKTdata(original_data.Hp,model_block,&original_data.A,&original_data.rhseqindex)){
      if (cb_out())
	get_out()<<"**** ERROR in QPSolver::solve(.......): attempt to initialize QPDirectKKTsolver failed as well, aborting"<<std::endl;
      return 1;
    }
    paramsp->QPset_KKTsolver(new QPDirectKKTSolver);
  }

						 
  //solve the QP
  int status=QPIsolve((skip_factor<0.),skip_factor);
  if (status){
    if (cb_out()){
      get_out()<<"**** WARNING in QPSolver::solve: QPIsolve returned "<<status<<std::endl;
    } 
  }

  if (postprocess_data(false)){
    if (cb_out()){
      get_out()<<"**** WARNING in QPSolver::solve: postprocess_data() failed "<<std::endl;
    } 
  }

  return status;
}


// *************************************************************************
//                             QPSolver::apply_modification
// *************************************************************************

int QPSolver::QPapply_modification(const GroundsetModification& mdf)
{
  const LPGroundsetModification* lpmdf=dynamic_cast<const LPGroundsetModification*>(&mdf);
  if (lpmdf){
    if ((original_data.A.rowdim()!=lpmdf->old_rowdim())||(original_data.dim!=lpmdf->old_vardim())){
      if (cb_out())
	get_out()<<"**** ERROR: QPSolver::apply_modification: there are "<<original_data.dim<<" variables and "<<original_data.A.rowdim()<<" constraints but modification assumes "<<lpmdf->old_vardim()<<" variables and "<<lpmdf->old_rowdim()<<" constraints"<<std::endl;
      return 1;
    }
    if (lpmdf->no_modification())
      return 0;

    int err=0;
    preproc_fixed.init(0,0,Integer(0));
    preproc_fixedval.init(0,0,0.);
    preproc_indices.init(0,0,Integer(0));
    preproc_aft.init();
    preproc_bundle_projection.clear();

    sol_val_lb=min_Real;
    sol_val_ub=max_Real;
    sol_point.init(0,0,0.);
    gs_aggr_offset=0.;
    gs_aggr_gradient.init(0,0,0.);
    
    original_data.dim=lpmdf->new_vardim();
    int retval=lpmdf->apply_to_bounds(original_data.lby,original_data.uby);
    if (retval){
      if (cb_out())
	get_out()<<"**** ERROR: QPSolver::apply_modification: apply to bounds failed and returned "<<retval<<std::endl;
      return retval;
    }
    retval=mdf.apply_to_costs(original_data.groundset_c,original_data.groundset_gamma);
    if  (retval){
      if (cb_out())
	get_out()<<"**** ERROR: QPSolver::apply_modification: apply to costs failed and returned "<<retval<<std::endl;
      return retval;
    }
    retval=lpmdf->apply_to_rows(original_data.A,original_data.rhslb,original_data.rhsub);
    if (retval){
      if (cb_out())
	get_out()<<"**** ERROR: QPSolver::apply_modification: apply to rows failed and returned "<<retval<<std::endl;
      err++;
    }
    if (err)
      return err;
  }
  else { //lpmdf==0, so this is some general Groundset, cannot rely on apply_modification
    if (original_data.dim!=mdf.old_vardim()){
      if (cb_out())
	get_out()<<"**** ERROR: QPSolver::apply_modification: there are "<<original_data.dim<<" variables, but modification assumes "<<mdf.old_vardim()<<" variables"<<std::endl;
      return 1;
    }
    if (mdf.no_modification())
      return 0;
    
    int err=0;
    original_data.dim=mdf.new_vardim();
    int retval=mdf.apply_to_costs(original_data.groundset_c,original_data.groundset_gamma);
    if  (retval){
      if (cb_out())
	get_out()<<"**** ERROR: QPSolver::apply_modification: apply to costs failed and returned "<<retval<<std::endl;
      return retval;
    }
    if (mdf.appended_vardim()>0){
      original_data.lby.concat_below(Matrix(mdf.appended_vardim(),1,CB_minus_infinity));
      original_data.uby.concat_below(Matrix(mdf.appended_vardim(),1,CB_plus_infinity));
      original_data.A.concat_right(Sparsemat(original_data.A.rowdim(),mdf.appended_vardim()));
    }
    if (mdf.map_to_old_variables()){
      original_data.lby=original_data.lby.rows(*(mdf.map_to_old_variables()));
      original_data.uby=original_data.uby.rows(*(mdf.map_to_old_variables()));
      original_data.A=original_data.A.cols(*(mdf.map_to_old_variables()));
    }
    if (err)
      return err;
  }
  int retval=determine_indices(original_data);
  if (retval){
    if (cb_out())
      get_out()<<"**** ERROR: QPSolver::apply_modification(): determine_indices() failed and retured"<<retval<<std::endl;
  }
  return retval;
}


// *************************************************************************
//                       QPSolver::mfile_data
// *************************************************************************

int QPSolver::mfile_data(std::ostream& out) const
{
  out<<"clear Aindi Aindj Aval A rhs lby uby\n";
  out<<"lby=[";
  for (Integer i=0;i<original_data.lby.dim();i++){
    out.precision(16);
    out.width(18);
    out<<original_data.lby(i);
    if (i<original_data.lby.dim()-1)
      out<<"\n";
  }
  out<<"];\n";
  out<<"uby=[";
  for (Integer i=0;i<original_data.lby.dim();i++){
    out.precision(16);
    out.width(18);
    out<<original_data.uby(i);
    if (i<original_data.uby.dim()-1)
      out<<"\n";
  }
  out<<"];\n";
  out<<"rhslb=[";
  for (Integer i=0;i<original_data.rhslb.dim();i++){
    out.precision(16);
    out.width(18);
    out<<original_data.rhslb(i);
    if (i<original_data.rhslb.dim()-1)
      out<<"\n";
  }
  out<<"];\n";
  out<<"rhsub=[";
  for (Integer i=0;i<original_data.rhsub.dim();i++){
    out.precision(16);
    out.width(18);
    out<<original_data.rhsub(i);
    if (i<original_data.rhsub.dim()-1)
      out<<"\n";
  }
  out<<"];\n";
  Indexmatrix indi,indj;
  Matrix val;
  original_data.A.get_edge_rep(indi,indj,val);
  out<<"Aindi=[";
  for (Integer i=0;i<indi.dim();i++){
    out.precision(16);
    out.width(18);
    out<<indi(i);
    if (i<indi.dim()-1)
      out<<"\n";
  }
  out<<"];\n";
  out<<"Aindj=[";
  for (Integer i=0;i<indj.dim();i++){
    out.precision(16);
    out.width(18);
    out<<indj(i);
    if (i<indj.dim()-1)
      out<<"\n";
  }
  out<<"];\n";
  out<<"Aval=[";
  for (Integer i=0;i<val.dim();i++){
    out.precision(16);
    out.width(18);
    out<<val(i);
    if (i<val.dim()-1)
      out<<"\n";
  }
  out<<"];\n";
  out<<"A=sparse(Aindi,Aindj,Aval,"<<original_data.A.rowdim()<<","<<original_data.A.coldim()<<");\n";
  return 0;
}


}

