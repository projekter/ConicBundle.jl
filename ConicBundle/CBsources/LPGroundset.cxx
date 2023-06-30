/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/LPGroundset.cxx
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
#include "mymath.hxx"
#include "LPGroundset.hxx"
#include "LPGroundsetModification.hxx"
#include "BundleIdProx.hxx"
#include "MatrixCBSolver.hxx"

using namespace CH_Matrix_Classes;

namespace ConicBundle {




// *****************************************************************************
//                              LPGroundset::LPGroundset
// *****************************************************************************

LPGroundset::LPGroundset(Integer dim,
			 const Matrix* lbyp,
			 const Matrix* ubyp,
			 const Sparsemat* Gp,
			 const Matrix* rhslbp,
			 const Matrix* rhsubp,
			 const Matrix* start_val,
			 const Matrix* costs,
			 const Real offset, 
			 Integer in_groundset_id,
			 CBout* cb):
  CBout(cb),vm_selection(0)
{
  qpsolver=new QPSolver(this);
  LPGroundset::clear(0,in_groundset_id);

  LPGroundsetModification mdf(0,0);
  mdf.add_append_vars(dim,lbyp,ubyp,0,start_val,costs);
  mdf.add_offset(offset);
  Integer nrows=0;
  if (Gp)
    nrows=Gp->rowdim();
  else if (rhslbp)
    nrows=rhslbp->dim();
  else if (rhsubp)
    nrows=rhsubp->dim();
  mdf.add_append_rows(nrows,Gp,rhslbp,rhsubp);

  apply_modification(mdf);
}

// *****************************************************************************
//                              LPGroundset::clear
// *****************************************************************************

  void LPGroundset::clear(Integer indim, Integer gs_id)
{
  Hp=0;
  dim=0;
  indim=(indim>0)?indim:0;
  assert(gs_id>=0);
  groundset_id=gs_id;
  
  gs_minorant.init(new Minorant,groundset_id);
  gs_aggregate=gs_minorant;
  use_yfixing=false;
  yfixed.init(dim,1,Integer(0));
  starting_point.init(dim,1,0.);
  c.init(0,1,0.);
  gamma=0;

  delete vm_selection;
  vm_selection=0;
  store_old_minorants=false;
  old_minorants.clear();
  max_minorants=100;
  minorant_nexti=0;
  old_lowrank.init(0,0,0.);
  old_diagonal.init(0,1,0.);
  old_sym.init(0,0.);

  qpsolver->QPclear();
  LPGroundsetModification mdf(0,0);
  mdf.add_append_vars(indim,0,0,0,0,0);
  groundset_id--;  //next operation increments the counter
  apply_modification(mdf);
}

// *****************************************************************************
//                           LPGroundset::get_lby()
// *****************************************************************************
  
int LPGroundset::set_qpsolver(QPSolverParametersObject* inqpparams,
			      QPSolverObject* inqpsolver)
{
  if (inqpsolver){
    delete qpsolver;
    qpsolver=inqpsolver;
  }
  if (inqpparams)
    qpsolver->QPset_parameters(inqpparams);
  return 0;
}

// *****************************************************************************
//                           LPGroundset::get_lby()
// *****************************************************************************
  
const Matrix* LPGroundset::get_lby() const
{
  const Matrix* p_lby; 
  const Matrix* p_uby;  
  const Indexmatrix* p_lbindex; 
  const Indexmatrix* p_ubindex; 

  qpsolver->QPboxconstrained(p_lby,p_uby,p_lbindex,p_ubindex);
  return p_lby;
}

// *****************************************************************************
//                           LPGroundset::get_uby()
// *****************************************************************************
  
const Matrix* LPGroundset::get_uby() const
{
  const Matrix* p_lby; 
  const Matrix* p_uby;  
  const Indexmatrix* p_lbindex; 
  const Indexmatrix* p_ubindex; 

  qpsolver->QPboxconstrained(p_lby,p_uby,p_lbindex,p_ubindex);
  return p_uby;
}

// *****************************************************************************
//                           LPGroundset::is_feasible()
// *****************************************************************************

bool LPGroundset::is_feasible(Integer& gs_id,const Matrix &y,Real relprec)
{
  if (gs_id==groundset_id)
    return true;

  if (!qpsolver->QPis_feasible(y,relprec))
    return false;
  
  gs_id=groundset_id;
  return true;
}

// *****************************************************************************
//                             ensure_feasibility
// *****************************************************************************

int LPGroundset::ensure_feasibility(Integer& gs_id, Matrix& y,bool& ychanged,BundleProxObject* inHp,Real relprec)
{
  if ((!ychanged)&&(gs_id==groundset_id)){
    return 0;
  }

  int status=qpsolver->QPensure_feasibility(y,ychanged,inHp,relprec);
  if ((status)&&(cb_out())){
      get_out()<<"**** ERROR in LPGroundset::ensure_feasibility(....): qpsolver->QPensure_fesibility returned "<<status<<std::endl;
  }

  gs_id=groundset_id;
  return status;
}

// *****************************************************************************
//                             get_qp_solver
// *****************************************************************************

QPSolverObject* LPGroundset::get_qp_solver(bool& qpwogs,
					   BundleProxObject* inHp)
{
  if (qpsolver->QPprefer_UQPSolver(inHp)) {  
    qpwogs=solve_model_without_gs= qpsolver->QPconstrained();
    qp_solver=&uqpsolver;
  }
  else {
    qpwogs=solve_model_without_gs=false; 
    qp_solver=qpsolver;
  }
  return qp_solver;
}
  

// *************************************************************************
//                              candidate
// *************************************************************************

//augval==CB_minus_infinity means that no good lower bound was available

int LPGroundset::candidate(Integer& gs_id,
			   Matrix& newy,
			   Real& cand_gs_val,
			   Real& linval,
			   Real& augval_lb,
			   Real& augval_ub,
			   Real& normsubg2,
			   const Matrix& center_y,
			   Real center_value,
			   const MinorantPointer& model_subg,
			   BundleProxObject* inHp,
			   MinorantPointer* delta_groundset_aggregate,
			   Indexmatrix* delta_index,
			   Real relprec)
{
  

  assert(inHp);
  //assert(augval_lb<center_value+max(eps_Real,relprec)*(fabs(center_value)+1.));
  assert((delta_groundset_aggregate==0)==(delta_index==0));
  bool initialize=(delta_groundset_aggregate==0);
  Hp=inHp;
  gs_id=groundset_id;

  if ((!solve_model_without_gs)&&(qp_solver==qpsolver)){
    //get candidate directly as the solution
    Real gsval;
    Matrix gsagg;
    qpsolver->QPget_solution(augval_lb,augval_ub,newy,gsval,gsagg);
    cand_gs_val=gs_minorant.evaluate(-1,newy);
    gs_aggregate=MinorantPointer(new MatrixMinorant(gsval,gsagg,0,true),groundset_id);
    MinorantPointer aggr(model_subg);
    aggr.aggregate(gs_aggregate,1.);
    normsubg2=Hp->dnorm_sqr(aggr);
    linval=aggr.evaluate(-1,newy);
    //assert(std::fabs(2.*(augval_lb-linval)-normsubg2)<=max(1e-6,relprec)*(std::fabs(augval_lb)+1.));

    if (delta_groundset_aggregate){
      delta_groundset_aggregate->init(new Minorant,groundset_id);
      assert(delta_index);
      delta_index->init(0,1,Integer(0));
    }

    // //TESTING BEGIN
    //  std::cout.precision(12);
    //  std::cout<<"TEST linval="<<linval<<"("<<model_subg.evaluate(-1,newy)<<"+"<<gs_aggregate.evaluate(-1,newy)<<"="<<aggr.evaluate(-1,newy)<<")";
    //  std::cout<<" augval_lb="<<augval_lb<<"(+="<<normsubg2/2.<<",<="<<aggr.evaluate(-1,center_y)-Hp->dnorm_sqr(aggr)/2.<<")"<<std::endl;
    //  std::cout<<" nsg2="<<normsubg2<<" HPdn(fa)="<<Hp->dnorm_sqr(aggr)<<" HPdn(ma)="<<Hp->dnorm_sqr(model_subg)<<" HPdn(ga)="<<Hp->dnorm_sqr(gs_aggregate)<<std::endl; 
    //   assert(std::abs(linval-model_subg.evaluate(-1,newy)-gs_aggregate.evaluate(-1,newy))<1e-6*(std::fabs(linval)+1.));
    //  assert(-1e-2*relprec*(std::fabs(augval_lb)+1.)<=aggr.evaluate(-1,center_y)-Hp->dnorm_sqr(aggr)/2.-augval_lb);
    // //TESTING END

    const Indexmatrix* lbind;
    const Matrix* lby;
    const Indexmatrix* ubind;
    const Matrix* uby;
    qp_solver->QPboxconstrained(lby,uby,lbind,ubind);
    if (use_yfixing){
      if (lbind){ 
	for(Integer i=0;i<lbind->dim();i++){
	  Integer ind=(*lbind)(i);
	  Real lb=(*lby)(ind);
	  if ((std::fabs(center_y(ind)-lb)<1e-10*max(1.,std::fabs(lb)))
	      &&(std::fabs(newy(ind)-lb)<1e-10*max(1.,std::fabs(lb)))
	      ){
	    if (yfixed(ind)==0) {
	      yfixed(ind)=2;
	    }
	  }
	}
      }
      if (ubind){
	for(Integer i=0;i<ubind->dim();i++){
	  Integer ind=(*ubind)(i);
	  Real ub=(*uby)(ind);
	  if ((std::fabs(center_y(ind)-ub)<1e-10*max(1.,std::fabs(ub)))
	      &&(std::fabs(newy(ind)-ub)<1e-10*max(1.,std::fabs(ub)))
	      ){
	    if (yfixed(ind)==0) {
	      yfixed(ind)=2;
	    }
	  }
	}
      }
    }

    Hp=0;
    return 0;
  }
  
  //first deal with the trivial case of no bounds and no constraints
  if (! qpsolver->QPconstrained()){
    Real dummy;
    newy.newsize(dim,1); chk_set_init(newy,1);
    model_subg.get_minorant(dummy,newy,0,-1.,false);
    gs_minorant.get_minorant(dummy,newy,0,-1.,true);
    Hp->apply_Hinv(newy);
    normsubg2=-model_subg.ip(newy)-gs_minorant.ip(newy);
    newy+=center_y;
    cand_gs_val=gs_aggregate.evaluate(-1,newy);
    linval=model_subg.evaluate(-1,newy)+cand_gs_val;
    //    assert(augval_lb<linval+normsubg2/2.+1e-12*(1.+std::fabs(linval+normsubg2/2.)));
    augval_lb=augval_ub=linval+normsubg2/2.;

    if (delta_groundset_aggregate){
      delta_groundset_aggregate->init(new Minorant,groundset_id);
      assert(delta_index);
      delta_index->init(0,1,Integer(0));
    }

    Hp=0;
    return 0;
  }
  
  //determine the rough structure of the scaling matrix
  Matrix diagQ;
  const Matrix* Vp;
  Hp->get_precond(diagQ,Vp);
  
  Real old_subg_offset=gs_aggregate.offset();
  Real gs_aggregate_offset=gs_minorant.offset();

  const Indexmatrix* p_lbindex; 
  const Matrix* p_lby; 
  const Indexmatrix* p_ubindex; 
  const Matrix* p_uby;  

  //next deal with the case that we have diagonal scaling and only box constraints
  if ((qpsolver->QPboxconstrained(p_lby,p_uby,p_lbindex,p_ubindex))&&(Hp->is_DLR())&&((Vp==0)||(Vp->coldim()==0))){
    const Indexmatrix& lbindex=*p_lbindex; //=qpsolver->get_lbindex();
    const Matrix& lby=*p_lby; //=qpsolver->get_lby();
    const Indexmatrix& ubindex=*p_ubindex; //=qpsolver->get_ubindex();
    const Matrix& uby=*p_uby;  //=qpsolver->get_uby();
    Real dummy;
    Matrix subg(dim,1); chk_set_init(subg,1); 
    model_subg.get_minorant(dummy,subg,0,1.,false);
    gs_minorant.get_minorant(dummy,subg,0,1.,true);
    newy.init(subg,-1.);
    newy/=diagQ;
    newy+=center_y;

    // //TESTING BEGIN
    // Matrix tmpnewy;
    // Matrix tmpvec;
    // Real tmpdummy;
    // Real tmpaugval;
    // if (qp_solver->QPget_solution(tmpaugval,tmpnewy,tmpdummy,tmpvec)==0){
    //   qp_solver->get_model_data_ptr()->add_modelx_aggregate(tmpdummy,tmpvec);
    //   std::cout<<" offsetdiff="<<std::fabs(dummy-tmpdummy);
    //   std::cout<<" subgdiff="<<norm2(subg-tmpvec);
    //   std::cout<<" newydiff="<<norm2(newy-tmpnewy);
    //   tmpnewy-=center_y;
    //   std::cout<<" augvaldiff="<<std::fabs(tmpaugval-(dummy+ip(subg,newy)+.5*normDsquared(tmpnewy,diagQ)));
    // }
    // //TESTING END
    
    //--- NOTE: This newy is not the optimal solution of the QP,
    //but the one that would result from using no bound constraints at all.
    //If you want to get the QP solution:
    //use gs_aggr instead of gs_minorant and set the
    //fixed values to center_y to get the optimal solution of the QP
    
    int cnt_fixed=0;
    int cnt_newfixed=0;
    int cnt_free=0;
    int cnt_changed=0;
    int cnt_increased=0;
    
    bool D_heuristic=Hp->employ_diagonal_bounds_scaling();
    Matrix D_update;
    if (D_heuristic)
      D_update.init(dim,1,0.);
    
    //--- update subg, if necessary 
    //    (so that the step satisfies the box constraints)
    Matrix update_value(dim,1); chk_set_init(update_value,1);
    Indexmatrix update_index(dim,1); chk_set_init(update_index,1);
    Matrix gs_aggr(dim,1); chk_set_init(gs_aggr,1);
    Real aggrdummy;
    gs_aggregate.get_minorant(aggrdummy,gs_aggr,0,1.,false);
    gs_minorant.get_minorant(aggrdummy,gs_aggr,0,-1.,true);
    
    Integer lbi=0;
    Integer lbind;
    if (lbi< lbindex.dim()){
      lbind=lbindex(lbi++);
    }
    else {
      lbind=dim;
    }
    Integer ubi=0;
    Integer ubind;
    if (ubi< ubindex.dim()){
      ubind=ubindex(ubi++);
    }
    else {
      ubind=dim;
    }
    Integer nz=0;
    while((lbind<dim)||(ubind<dim)){
      Integer bind=min(lbind,ubind);
      if (bind==lbind){
	if (lbi< lbindex.dim()){
	  lbind=lbindex[lbi++];
	}
	else {
	  lbind=dim;
	}
      }
      if (bind==ubind){
	if (ubi< ubindex.dim()){
	  ubind=ubindex[ubi++];
	}
	else {
	  ubind=dim;
	}
      } 
      Real old_subg=gs_aggr(bind);
      Real d=newy(bind);
      if ((use_yfixing)&&(!initialize)&&(yfixed(bind)>0)){
	yfixed(bind)=1; //flag for was fixed and stays fixed
	cnt_fixed++;
	//gs_aggregate(bind)=(d-center_y(bind))*diagQ(bind);  //if gs_aggregate needs to reflect the corresponding "slack" this would be the formula
	gs_aggr(bind)=0.;
	newy(bind)=center_y(bind);
	continue;
      }
      if (d<lby(bind)-eps_Real*(fabs(lby(bind))+1.)){
	//--- new value would be smaller than lower bound
	gs_aggr(bind)=(d-lby(bind))*diagQ(bind);   // <0
	//gs_aggregate_offset will be modified afterwards
	newy(bind)=lby(bind);
	if ((use_yfixing)&&(center_y(bind)<=lby(bind)+eps_Real*(fabs(lby(bind))+1.))){
	  //y is at lower bound, fix this value for the time being
	  gs_aggr(bind)=0;
	  if ((initialize)||(yfixed(bind)>0)){
	    yfixed(bind)=1;  //flag for was fixed and stays fixed
	    cnt_fixed++;
	  }
	  else {
	    yfixed(bind)=2;  //flag for newly fixed
	    cnt_newfixed++;
	    cnt_fixed++;
	    update_value(nz)=-old_subg;
	    update_index(nz)=bind;
	    nz++;
	  }
	}
	else {
	  //y is above lower bound
	  Real diff=gs_aggr(bind)-old_subg;
	  if (fabs(diff)>eps_Real*(fabs(old_subg)+1.)){
            if ((D_heuristic)&&(diff<-eps_Real*(fabs(old_subg)+1.))){
	      Real oldD=diagQ(bind);
	      Real newD=-subg(bind)/(0.9*lby(bind)+.1*d-center_y(bind));
	      newD=min(newD,max(oldD,(fabs(lby(bind))+1.)*1e6));
	      Real Dupd=newD-oldD;
	      if (Dupd>eps_Real*(oldD+1.)){
		diagQ(bind)=newD;
		D_update(bind)=Dupd;
		d=center_y(bind)-subg(bind)/diagQ(bind);
		assert(d<=lby(bind));
		gs_aggr(bind)=(d-lby(bind))*diagQ(bind);
		diff=gs_aggr(bind)-old_subg;
		cnt_increased++;
	      }
	    }
	    yfixed(bind)=0;  //flag for not fixed
	    cnt_free++;
	    update_value(nz)=diff;
	    update_index(nz)=bind;
	    cnt_changed++;
	    nz++;
	  }
	} 
	gs_aggregate_offset-=lby(bind)*gs_aggr(bind);
      }
      else if (d>uby(bind)+eps_Real*(fabs(uby(bind))+1.)){
	//--- new value would be greater than upper bound
	gs_aggr(bind)=(d-uby(bind))*diagQ(bind);   // >0
	newy(bind)=uby(bind);
	if ((use_yfixing)&&(center_y(bind)>=uby(bind)-eps_Real*(fabs(uby(bind))+1.))){
	  gs_aggr(bind)=0;
	  //y is at upper bound, fix this value for the time being
	  if ((initialize)||(yfixed(bind)>0)){
	    yfixed(bind)=1;  //flag for was fixed and stays fixed
	    cnt_fixed++;
	  }
	  else {
	    yfixed(bind)=2;  //flag for newly fixed
	    cnt_newfixed++;
	    cnt_fixed++;
	    update_value(nz)=-old_subg;
	    update_index(nz)=bind;
	    nz++;
	  }
	}
	else {
	  // y is below upper bound
	  Real diff=gs_aggr(bind)-old_subg;
	  if (fabs(diff)>eps_Real*(fabs(old_subg)+1.)){
            if ((D_heuristic)&&(diff>eps_Real*(fabs(old_subg)+1.))){
	      Real oldD=diagQ(bind);
	      Real newD=-subg(bind)/(0.9*uby(bind)+.1*d-center_y(bind));
	      newD=min(newD,max(oldD,(fabs(uby(bind))+1.)*1e6));
	      Real Dupd=newD-oldD;
	      if (Dupd>eps_Real*(oldD+1.)){
		diagQ(bind)=newD;
		D_update(bind)=Dupd;
		d=center_y(bind)-subg(bind)/diagQ(bind);
		assert(d>uby(bind));
		gs_aggr(bind)=(d-uby(bind))*diagQ(bind);
		diff=gs_aggr(bind)-old_subg;
		cnt_increased++;
	      }
	    }
	    yfixed(bind)=0;  //flag for not fixed 
	    cnt_free++;
	    update_value(nz)=gs_aggr(bind)-old_subg;
	    update_index(nz)=bind;
	    cnt_changed++;
	    nz++;
	  }
	} 
	gs_aggregate_offset-=uby(bind)*gs_aggr(bind);
      }
      else {  
	//--- new value within bounds
	yfixed(bind)=0; // flag for not fixed 
	cnt_free++;
	gs_aggr(bind)=0.;
	if (old_subg!=0.) {
	  update_value(nz)=-old_subg;
	  update_index(nz)=bind;
	  cnt_changed++;
	  nz++;
	}
      }
    }
    update_value.reduce_length(nz);
    update_index.reduce_length(nz);
    if (cb_out(0)){
      get_out()<<" init="<<initialize;
      get_out()<<" free="<<cnt_free;
      get_out()<<" fixed="<<cnt_fixed;
      get_out()<<" (+"<<cnt_newfixed<<") ";
      if(D_heuristic){
	get_out()<<" incr="<<cnt_increased<<"(<="<<max(D_update(update_index))<<")";
      }
      get_out()<<" changed="<<cnt_changed<<"("<<norm2(update_value)<<")\n";
    }

    
    newy-=center_y;
    normsubg2=normDsquared(newy,diagQ);
    newy+=center_y;

    //by complementarity we may drop the new subgradients contribution to linval
    //cand_gs_val=gs_minorant.evaluate(-1,newy);

    gs_minorant.get_minorant(dummy,gs_aggr,0,1.,true);
    Minorant* minorant=new Minorant(true,gs_aggregate_offset);
    minorant->add_coeffs(gs_aggr.dim(),gs_aggr.get_store());
    //minorant->sparsify(1e-12*std::sqrt(normsubg2));
    gs_aggregate.init(minorant,groundset_id);
    cand_gs_val=gs_aggregate.evaluate(-1,newy);

    linval=model_subg.evaluate(-1,newy)+cand_gs_val;
    if ((cnt_changed>0)&&(augval_lb-relprec*(fabs(augval_lb)+1.)>linval+normsubg2/2.)){
      if (cb_out(2)){
	get_out().precision(16);
	get_out()<<"**** WARNING LPGroundset::candidate(): maybe numerical problems ahead: current value of augval_lb="<<augval_lb<<" decreases to "<<linval+normsubg2/2.<<std::endl;
      }
    }
    augval_lb=augval_ub=linval+normsubg2/2.;

    if(store_old_minorants){
      if (minorant_nexti==Integer(old_minorants.size())){
	old_minorants.push_back(gs_aggregate);
      }
      else {
	old_minorants[unsigned(minorant_nexti)]=gs_aggregate;
      }
      minorant_nexti++;
      minorant_nexti%=max_minorants;
    }
    
    if (delta_groundset_aggregate){
      Minorant* minorant=new Minorant(true,gs_aggregate_offset-old_subg_offset);
      minorant->add_coeffs(update_value.dim(),update_value.get_store(),update_index.get_store());
      delta_groundset_aggregate->init(minorant,groundset_id);
      delta_index->init(update_index);
    }

    if (D_heuristic) {
      if (inHp->diagonal_bounds_scaling_update(D_update)){
	if (cb_out()){
	  get_out()<<"**** ERROR LPGroundset::candidate(): inHp->diagonal_bounds_scaling_update failed even though inHp->supports_diagonal_scaling_heuristic returned true"<<std::endl;
	}
        Hp=0;
	return 1;
      }	  
    }

    Hp=0;
    return 0;
  }
    
  //for the remaining cases we solve a QP

  //if augval==CB_minus_infinity compute a useful lower bound to
  //allow for a reasonable stopping criterion in th QP code
  Real dummy;
  Matrix tmpvec;
  //bool tmpvecinit=false;
  if (augval_lb==CB_minus_infinity){
    //tmpvecinit=true;
    tmpvec.newsize(dim,1); chk_set_init(tmpvec,1);
    model_subg.get_minorant(dummy,tmpvec,0,1.,false);
    gs_aggregate.get_minorant(dummy,tmpvec,0,1.,true);
    c.init(tmpvec);
    Hp->apply_Hinv(c);
    augval_lb=augval_ub=model_subg.offset()+gs_aggregate.offset()+ip(tmpvec,center_y)-.5*ip(tmpvec,c);
  }

  //compute the cost coefficients for the bundle subproblem

  /*
  c.init(dim,1,0.);
  Hp->add_Hx(center_y,c);
  gamma=gs_minorant.offset()+.5*ip(c,center_y);
  c*=-1;
  gs_minorant.get_minorant(dummy,c,0,1.,true);
  if (qpsolver->get_model_block_ptr()==0){
    model_subg.get_minorant(gamma,c,0,1.,true);
  }
  else {
    qpsolver->get_model_block_ptr()->get_constant_minorant().get_minorant(gamma,c,0,1.,true);
  }
  */

    /*
  if (!tmpvecinit){
    tmpvec.newsize(dim,1); chk_set_init(tmpvec,1);
    model_subg.get_minorant(dummy,tmpvec,0,1.,false);
    gs_aggregate.get_minorant(dummy,tmpvec,0,1.,true);
  }
  paramsp->QPset_dual_infeasibility_eps(dual_infeasibility_eps*(center_value-augval_lb+norm2(tmpvec)/Hp->get_term_corr()));
  std::cout<<" die="<<dual_infeasibility_eps<<" cv-av="<<center_value-augval_lb<<" n2(tv)="<<norm2(tmpvec)<<" tcorr="<<Hp->get_term_corr()<<" ntv/tc="<<norm2(tmpvec)/Hp->get_term_corr()<<std::endl;
  */


  //determine the required precision and call the solver
  relprec=min(max(eps_Real,relprec),0.05*(max(eps_Real,center_value-augval_lb))/(fabs(center_value)+1.));
  //Real skip_factor=(delta_groundset_aggregate==0)?.4:.7;
  //int status=qpsolver->solve(Hp,c,gamma,augval_lb,max(augval_lb+eps_Real,center_value),relprec,skip_factor);
  int status=qpsolver->QPsolve(center_y,augval_lb,max(augval_lb+eps_Real,center_value),relprec,
			       Hp,gs_aggregate,&yfixed);
  if ((status)&&(cb_out())){
    get_out()<<"**** WARNING in LPGroundset::candidate(): qpsolver->QPsolve() returned "<<status<<std::endl;
  }
  
  
  //retrieve the solution and the new groundset aggregate
  if (delta_groundset_aggregate)
    *delta_groundset_aggregate=gs_aggregate;
  Real dummy_lb,dummy_ub;
  qpsolver->QPget_solution(dummy_lb,dummy_ub,newy,gs_aggregate_offset,tmpvec);
  Minorant* minorant=new MatrixMinorant(gs_aggregate_offset,tmpvec,0,true);
  gs_aggregate.init(minorant,groundset_id);

  newy-=center_y;
  normsubg2=Hp->norm_sqr(newy);
  

  //TEST
  /*
  Matrix tmpmat(gs_aggregate);
  tmpmat+=model_subg;
  Hp->apply_Hinv(tmpmat);
  if (norm2(tmpmat+newy)>1e-6*norm2(newy)){
    if (cb_out()){
      get_out()<<"**** WARNING in LPGroundset::candidate: Hinv*aggregate and step differ by norm="<<norm2(tmpmat+newy)<<std::endl;
    } 
  }
  */
    
  newy+=center_y;

  if(store_old_minorants){
    if (minorant_nexti==Integer(old_minorants.size())){
      old_minorants.push_back(gs_aggregate);
    }
    else {
      old_minorants[unsigned(minorant_nexti)]=gs_aggregate;
    }
    minorant_nexti++;
    minorant_nexti%=max_minorants;
  }
    

  if (delta_groundset_aggregate){
     delta_groundset_aggregate->scale(-1.);
     gs_aggregate.get_minorant(*delta_groundset_aggregate,1.);
    Integer nz;
    const Real* valp;
    const Integer *indp;
    delta_groundset_aggregate->get_minorant()->get_coeffs(nz,valp,indp);
    if (indp!=0){
      delta_index->init(nz,1,indp);
    }
    else {
      delta_index->init(Range(0,nz-1));
    }
   /*
    delta_groundset_aggregate->get_minorant(gs_aggregate_offset,tmpvec,0,-1,true);
    minorant=new Minorant(true,gs_aggregate_offset);
    minorant->add_coeffs(tmpvec.dim(),tmpvec.get_store());
    minorant->sparsify(1e-12*normsubg2);
    delta_groundset_aggregate->init(minorant,0);
    Integer nz;
    const Real* valp;
    const Integer *indp;
    minorant->get_coeffs(nz,valp,indp);
    if (indp!=0){
      delta_index->init(nz,1,indp);
    }
    else {
      delta_index->init(Range(0,nz-1));
    }
    */
  }

  Real oldaugval=augval_lb;
  cand_gs_val=gs_aggregate.evaluate(-1,newy);
  linval=cand_gs_val+model_subg.evaluate(-1,newy);
  augval_lb=augval_ub=linval+normsubg2/2.;
  if ((normsubg2>relprec)&&(augval_lb<oldaugval)){
    if (cb_out()){
      get_out().precision(10);
      get_out()<<"**** WARNING in LPGroundset::candidate: Hp->candidate failed to increase augval="<<oldaugval<<" by producing value="<<augval_lb<<status<<std::endl;
    } 
  }
  Hp=0;
  return status;
}

// *************************************************************************
//                       LPGroundset::mfile_data
// *************************************************************************

int LPGroundset::mfile_data(std::ostream& out) const
{
  out<<"clear gs_subg gs_sugb_offset yfixed G rhs lby uby\n";
  out<<"gs_subg=[";
  for (Integer i=0;i<dim;i++){
    out.precision(16);
    out.width(18);
    out<<gs_aggregate.coeff(i);
    if (i<dim-1)
      out<<"\n";
  }
  out<<"];\n";
  out<<"gs_subg_offset="<<gs_aggregate.offset()<<";\n";
  out<<"yfixed=[";
  for (Integer i=0;i<yfixed.dim();i++){
    out.precision(16);
    out.width(18);
    out<<yfixed(i);
    if (i<yfixed.dim()-1)
      out<<"\n";
  }
  out<<"];\n"; 
  return qpsolver->mfile_data(out);
}

// *****************************************************************************
//                                add_variable_metric
// *****************************************************************************

int LPGroundset::add_variable_metric(VariableMetric& H,
				     Integer y_id,
				     const Matrix& y,
				     bool descent_step,
				     Real weightu,
				     Real model_maxviol,
				     const Indexmatrix* indices)
{
  if (!H.employ_variable_metric()){
    old_minorants.clear();
    minorant_nexti=0;
    old_lowrank.init(0,0,0.);
    old_diagonal.init(0,1,0.);
    old_sym.init(0,0.);
    return 0;
  }
  
  store_old_minorants=true;
      
  if (gs_aggregate.empty()||(old_minorants.size()<=1)){
    old_lowrank.init(0,0,0.);
    old_diagonal.init(0,1,0.);
    old_sym.init(0,0.);

    return 0;
  }

  class VariableMetricLPGroundsetData: public VariableMetricBundleData
  {
  public:
    Integer* max_minorants;
    Integer minorant_nexti;
    MinorantBundle* old_minorants;
    Symmatrix* old_sym;
    Matrix* old_diagonal;
    Matrix* old_lowrank;
    MinorantPointer* aggregate;    

    Real get_function_factor() const
    { return 1.;}

    int get_latest_minorants(MinorantBundle& latest_minorants,
			     Integer max_number)
    {
      *max_minorants=max(*max_minorants,max_number);
      if (Integer(old_minorants->size())<=max_number)
	latest_minorants=*old_minorants;
      else if (max_number>0){
	latest_minorants.resize(unsigned(max_number));
	long unsigned int ind=unsigned(minorant_nexti);
	for(unsigned int i=0;i<latest_minorants.size();i++){
	  if (ind==0)
	    ind=old_minorants->size();
	  --ind;
	  latest_minorants[i]=(*old_minorants)[ind];
	}
      }
      else {
	latest_minorants.clear();
      }
      return 0;
    }

    int get_model_data(MinorantBundle& model_minorants,Matrix& model_coeff) const
    {
      model_minorants.clear();
      model_coeff.init(0,0,0.);
      return 0;
    }

    const MinorantPointer& get_aggregate() const
    { return *aggregate;}

    Symmatrix& set_denseH()
    { return *old_sym; }

    Matrix& set_lowrankH()
    { return *old_lowrank; }

    Matrix& set_diagH()
    { return *old_diagonal; }
  };


  VariableMetricSelection* vms=vm_selection;
  if (vm_selection==0)
    vm_selection=H.get_variable_metric_selection();
  if (vm_selection!=0){
    VariableMetricLPGroundsetData vmsbd;
    vmsbd.max_minorants=&max_minorants;
    vmsbd.minorant_nexti=minorant_nexti;
    vmsbd.old_sym=&old_sym;
    vmsbd.old_diagonal=&old_diagonal;
    vmsbd.old_lowrank=&old_lowrank;
    vmsbd.aggregate=&gs_aggregate;

    if (vms->add_variable_metric(H,y_id,y,descent_step,
				 weightu,model_maxviol,
				 indices,vmsbd)){
      if (cb_out())
      get_out()<<"**** WARNING: LPGroundset::add_dynamic_scaling(): vm_selection->add_variable_metric(........) failed"<<std::endl;

      return 1;
    }
  }

  return 0;
}    

// *****************************************************************************
//                           apply_modification
// *****************************************************************************

int LPGroundset::apply_modification(const GroundsetModification& mdf)
{
  if (dim!=mdf.old_vardim()){
    if (cb_out())
      get_out()<<"**** ERROR: LPGroundset::apply_modification: there are "<<dim<<" variables but modification assumes "<<mdf.old_vardim()<<" variables"<<std::endl;
    return 1;
  }

  int err=0;
  dim=mdf.new_vardim();
  yfixed.init(dim,1,0);
  groundset_id++; 
  if (gs_minorant.apply_modification(mdf,groundset_id,0,true)){
    if (cb_out())
      get_out()<<"**** ERROR: LPGroundset::apply_modification(.): modification of the groundset minorant failed"<<std::endl;
    err++;
  }
  gs_aggregate=gs_minorant;
  
  if (mdf.apply_to_vars(starting_point)){
    if (cb_out())
      get_out()<<"**** ERROR: Groundset::apply_modification(.): modification of the starting point failed"<<std::endl;
    err++;
  }
  
  int retval = qpsolver->QPapply_modification(mdf);
  if (retval){
    if (cb_out())
      get_out()<<"**** ERROR: LPGroundset::apply_modification: apply to qpsolver failed and returned "<<retval<<std::endl;
    err++;
  }
  if (err)
    return err;
  
  int dummy=-1;
  if (!is_feasible(dummy,starting_point)){
    if (cb_out())
      get_out()<<"**** WARNING: Groundset::apply_modification(.): starting point is not feasible (but this is allowed)"<<std::endl;
  }
  return retval;
}



}

