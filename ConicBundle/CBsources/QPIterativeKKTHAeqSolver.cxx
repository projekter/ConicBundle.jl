/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/QPIterativeKKTHAeqSolver.cxx
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


#include "pcg.hxx"
#include "QPIterativeKKTHAeqSolver.hxx"
#include "QPKKTSubspaceHPrecond.hxx"
#include "clock.hxx"

using namespace CH_Matrix_Classes;

namespace ConicBundle {

  QPIterativeKKTHAeqSolver::~QPIterativeKKTHAeqSolver()
  {}
 
// *************************************************************************
//                            QPinit_KKTdata
// *************************************************************************

int QPIterativeKKTHAeqSolver::QPinit_KKTdata(QPSolverProxObject* in_Hp,
					     QPModelBlockObject* in_model, 
					     const Sparsemat* in_A, 
					     const Indexmatrix* in_eq_indices)
{
  assert(solver);
  assert(in_Hp);
  clear();
  
  Hp=in_Hp;
  model=in_model;
  A=in_A;
  eq_indices=in_eq_indices;

  if (A)
    blockA_norm=max(1.,norm2(*A));

  // if ((dynamic_cast<CH_Matrix_Classes::PCG*>(solver))&&(A||model)&&cb_out())
  //   get_out()<<"**** WARNING in QPIterativeKKTSolver::QPinit_KKTdata(): a PCG solver might not be suitable for this likely indefinite system"<<std::endl;

  int err=0;
  if (precond){
    err=precond->init_data(Hp,model,A,eq_indices,true);
    if ((err)&&(cb_out()))
      get_out()<<"**** WARNING in QPIterativeKKTSolver::QPinit_KKTdata(): precond->init_data() returned"<<err<<std::endl;
  }
    

  sol.init(0,0,0.);
  xrecord.init(0,0,0.);
  //if (dynamic_cast<QPKKTSubspaceHPrecond*>(precond))
  //  recordstep=10;
  //else 
    recordstep=-1;

  maxit_bnd=-1;
  nmult=0;

  return err;
}


// *************************************************************************
//                             QPsolve_KKTsystem
// *************************************************************************

// solve the KKT System

/// on input
/// dualrhs = -(Qx+A'y+G'modelx+c) +mu (...)
/// primalrhs = -(Ax+s) + mu ()...
int QPIterativeKKTHAeqSolver::QPsolve_KKTsystem(Matrix& solx,Matrix& soly,
					 const Matrix& primalrhs,
					 const Matrix& dualrhs,
					 Real rhsmu,
					 Real rhscorr,
					 Real prec,
					 QPSolverParameters* /* params */)
{
  clock.start();
  t_itsys_mult=0;
  if (precond)
    precond->reset_t_precond_mult();
      
  maxit_bnd=max(100*Integer(std::log(primalrhs.rowdim()+dualrhs.rowdim())),min(maxit_bnd,2*(primalrhs.rowdim()+dualrhs.rowdim())));

  //---- prepare the right hand side
  sysrhs.init(dualrhs);
  int status1=0;
  if (model){
    status1=model->add_Schur_rhs(sysrhs,0,rhsmu,rhscorr);
    if (status1){
      if (cb_out()){
	get_out()<<"**** WARNING: QPIterativeKKTHAeqSolver::QPsolve_KKTsystem(....): model->add_Schur_rhs failed and returned "<<status1<<std::endl;
      }
    }
  }
  Integer eqnr=(eq_indices)?eq_indices->rowdim():0;
  if (primalrhs.rowdim()>eqnr){
    Matrix tmpvec(primalrhs);
    Integer j=0;
    for(Integer i=0;i<eqnr;i++){
      Integer ii=(*eq_indices)(i);
      for(;j<ii;j++){
	tmpvec(j)/=KKTdiagy(j);
      }
      tmpvec(j++)=0;
    }
    for(;j<tmpvec.rowdim();j++){
      tmpvec(j)/=KKTdiagy(j);
    }
    genmult(*A,tmpvec,sysrhs,1.,1.,1);
    if (eqnr)
      sysrhs.concat_below(primalrhs.rows(*eq_indices));
  }
  else {
    sysrhs.concat_below(primalrhs);
  }

  //--- TEST on diagonal
  /*
  Matrix diagH;
  const Matrix *Vp;
  Hp->get_precond(diagH,Vp);
  diagH+=KKTdiagx;
  std::cout.precision(10);
  if (Vp){
    const Real* vp=Vp->get_store();
    for(Integer j=0;j<Vp->coldim();j++){
      Real* dp=diagH.get_store();
      const Real* const dpend=dp+diagH.rowdim();
      for(;dp!=dpend;vp++,dp++)
	(*dp)+=(*vp)*(*vp);
    }
  } //endif (Vp)
  //std::cout<<"\n dH="<<transpose(diagH);
  
  if (model){
    model->add_BCSchur_diagonal(diagH);
  }
  //std::cout<<" dM="<<transpose(diagH);
	
  if(A){	  
    Integer eqi=0;
    Integer eqind=(eqi<eqnr)?(*eq_indices)(eqi++):A->rowdim();
    for(Integer rinfoi=0;rinfoi<A->get_rowinfo().rowdim();rinfoi++){
      Integer rind=A->get_rowinfo()(rinfoi,0);
      while(eqind<rind)
	eqind=(eqi<eqnr)?(*eq_indices)(eqi++):A->rowdim();
      if (rind<eqind){
	Integer rnz=A->get_rowinfo()(rinfoi,1);
	Integer rbase=A->get_rowinfo()(rinfoi,2);
	Real d=1./KKTdiagy(rind);
	for(Integer i=0;i<rnz;i++){
	  diagH(A->get_rowindex()(rbase+i))+=d*sqr(A->get_rowval()(rbase+i));
	}		
      }
    }
  } //endif
  //std::cout<<" dA="<<transpose(diagH);
	

  Real derr=0.;
  for(Integer i=0;i<diagH.rowdim();i++){
    Matrix testinvec(diagH.rowdim()+eqnr,1,0.);
    testinvec(i)=1.;
    Matrix testoutvec;
    ItSys_mult(testinvec,testoutvec);
    if (std::fabs(testoutvec(i)-diagH(i))>1e-6*(std::fabs(testoutvec(i))+1.)){
      std::cout<<" ["<<i<<","<<testoutvec(i)<<","<<diagH(i)<<"] "<<std::flush;
    }
    derr+=sqr(testoutvec(i)-diagH(i));
  }
  if (std::sqrt(derr)>1e-8){
    std::cout<<" derr="<<std::sqrt(derr)<<std::endl;
  }
  */
  
  //---- call the solver
  assert(solver);
  Real termprec=prec*min(1.,norm2(sysrhs));
  if (precond) {
    termprec*=std::sqrt(precond->get_lmin_invM1());
  }
  Matrix *xsp=0;
  Integer ntries=5;
  int status2;
  if (sol.dim()==sysrhs.dim()){
    //perturb the old vector a little bit
    Matrix tmpvec;
    tmpvec.rand(sol.rowdim(),1);
    tmpvec-=.5;
    tmpvec*=(1e-6*norm2(sol)/sol.dim());
    sol+=tmpvec;
  }
  do {
    if (recordstep>=0){
      xrecord.newsize(sysrhs.rowdim(),100);
      xrecord.init(sysrhs.rowdim(),0,0.);
      xsp=&xrecord;
    }
    solver->set_maxit(maxit_bnd);
    status2=solver->compute(*this,sol,termprec,xsp,recordstep);
    maxit_bnd=max(maxit_bnd,Integer(1.2*solver->get_nmult()));
    
    if (xsp){
      xrecord.concat_right(sol);
      xrecord.delete_rows(Range(KKTdiagx.rowdim(),xrecord.rowdim()-1));
      dynamic_cast<QPKKTSubspaceHPrecond*>(precond)->set_subspace(xrecord);
    }
    if ((solver->get_residual_norm()>1.1*solver->get_termprec())&&(solver->get_nmult()>=solver->get_maxit())){
      if (termprec<1e-7*min(1.,norm2(sysrhs)))
	termprec*=10.;
      Matrix tmpvec;
      tmpvec.rand(sol.rowdim(),1);
      tmpvec-=.5;
      tmpvec*=(0.1*solver->get_residual_norm()/norm2(tmpvec));
      sol+=tmpvec;
    }
  } while (
	   (solver->get_residual_norm()>1.1*solver->get_termprec())
	   &&(solver->get_nmult()>=solver->get_maxit())
	   &&(--ntries>0)
	   );
	   
  if (status2){
    if (cb_out()){
      get_out()<<"**** WARNING: QPIterativeKKTHAeqSolver::QPsolve_KKTsystem(....): solver->compute failed and returned "<<status2<<std::endl;
    }
  }

  //---- split up the solution
  solx.init(dualrhs.rowdim(),1,sol.get_store());
  if (primalrhs.rowdim()>eqnr){
    soly.init(primalrhs,-1.);
    genmult(*A,solx,soly,1.,1.);
    Integer j=0;
    for(Integer i=0;i<eqnr;i++){
      Integer ii=(*eq_indices)(i);
      for(;j<ii;j++){
	soly(j)/=KKTdiagy(j);
      }
      soly(j++)=sol(dualrhs.rowdim()+i);
    }
    for(;j<soly.rowdim();j++){
      soly(j)/=KKTdiagy(j);
    }
  }
  else {
    soly.init(primalrhs.rowdim(),1,sol.get_store()+dualrhs.rowdim());
  }

  int status3=0;
  if (model){
    status3=model->compute_step(solx);
    if (status3){
      if (cb_out()){
	get_out()<<"**** WARNING: QPIterativeKKTHAeqSolver::QPsolve_KKTsystem(....): model->compute_step failed and returned "<<status3<<std::endl;
      }
    }
  }

  if (cb_out(3)){
    get_out()<<" solt "<< clock.time()<<"["<<t_itsys_mult;
    if (precond)
      get_out()<<","<<precond->get_t_precond_mult();
    get_out()<<"]";
  }

  return std::abs(status1)+std::abs(status2)+std::abs(status3);
}

// *************************************************************************
//                             ItSys_mult
// *************************************************************************

int QPIterativeKKTHAeqSolver::ItSys_mult(const Matrix& in_vec,Matrix& out_vec)
{
  CH_Tools::Microseconds t_start=clock.time();
  
  nmult++;
  
  out_vec.newsize(in_vec.rowdim(),1);
  in_vecx.init(KKTdiagx.rowdim(),1,in_vec.get_store());

  out_vec=in_vecx;
  out_vec%=KKTdiagx;
  Hp->add_Hx(in_vecx,out_vec,Hfactor);
  //std::cout<<"\n out_vecH="<<transpose(out_vec);
    
  if (model){
    model->add_Schur_mult(in_vecx,out_vec);
  }
  //std::cout<<"out_vecM="<<transpose(out_vec);
    
  if (KKTdiagy.rowdim()>0) {
    Integer eqnr=(eq_indices)?eq_indices->rowdim():0;
    if (KKTdiagy.rowdim()==eqnr){
      assert(norm2(KKTdiagy)<1e-6*KKTdiagy.rowdim());
      in_vecy.init(KKTdiagy.rowdim(),1,in_vec.get_store()+KKTdiagx.rowdim());
      genmult(*A,in_vecy,out_vec,1.,1.,1);
      //in_vecy%=KKTdiagy;
      genmult(*A,in_vecx,in_vecy,1.,0.);
      out_vec.concat_below(in_vecy);
    }
    else if (eqnr==0) {
      genmult(*A,in_vecx,in_vecy,1.,0.);
      in_vecy/=KKTdiagy;
      genmult(*A,in_vecy,out_vec,1.,1.,1);
    }
    else {
      genmult(*A,in_vecx,in_vecy,1.,0.);
      Matrix tmpvec(eqnr,1); chk_set_init(tmpvec,1);
      Integer j=0;
      for(Integer i=0;i<eqnr;i++){
	Integer ii=(*eq_indices)(i);
	for(;j<ii;j++){
	  in_vecy(j)/=KKTdiagy(j);
	}
	tmpvec(i)=in_vecy(j);
	in_vecy(j++)=in_vec(KKTdiagx.rowdim()+i);
      }
      for(;j<in_vecy.rowdim();j++){
	in_vecy(j)/=KKTdiagy(j);
      }
      genmult(*A,in_vecy,out_vec,1.,1.,1);
      out_vec.concat_below(tmpvec);
    }
  }
  //std::cout<<"out_vecA="<<transpose(out_vec);

  Real xnormsqr=mat_ip(KKTdiagx.rowdim(),in_vec.get_store());
  if (xnormsqr>1e-10){
    Real Ritz=mat_ip(KKTdiagx.rowdim(),in_vec.get_store(),out_vec.get_store());
    Ritz/=xnormsqr;
    if (Ritz>blockH_norm)
      blockH_norm=Ritz;
  }

  t_itsys_mult+=clock.time()-t_start;

  return 0;
}


}

