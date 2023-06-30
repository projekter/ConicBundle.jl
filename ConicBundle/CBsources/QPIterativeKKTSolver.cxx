/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/QPIterativeKKTSolver.cxx
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


#include "QPIterativeKKTSolver.hxx"
#include "lanczpol.hxx"

using namespace CH_Matrix_Classes;

namespace ConicBundle {

  QPIterativeKKTSolver::~QPIterativeKKTSolver()
  {delete solver; delete precond;}
 
// *************************************************************************
//                             clear
// *************************************************************************

void QPIterativeKKTSolver::clear()
{
  if (precond)
    precond->clear();
  
  Hfactor=1.;
  QPKKTSolverObject::clear();

}


// *************************************************************************
//                            QPinit_KKTdata
// *************************************************************************

int QPIterativeKKTSolver::QPinit_KKTdata(QPSolverProxObject* in_Hp,
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
    err=precond->init_data(Hp,model,A,eq_indices,false);
    if ((err)&&(cb_out()))
      get_out()<<"**** WARNING in QPIterativeKKTSolver::QPinit_KKTdata(): precond->init_data() returned"<<err<<std::endl;
  }
    

  sol.init(0,0,0.);
  
  maxit_bnd=-1;
  nmult=0;

  return err;
}


// *************************************************************************
//                             QPinit_KKTsystem
// *************************************************************************

// set up the KKT System

// KKTdiagx and KTTdiagy are both the positive compelementary systems,
// but the one of KKTdiagy is on the dual side to KKTdiagx and taken with minus 
  
  
int QPIterativeKKTSolver::QPinit_KKTsystem(const Matrix& inKKTdiagx,
					   const Matrix& inKKTdiagy,
					   Real Hfac,
					   Real prec,
					   QPSolverParameters* params)
{
  nmult=0;

  Hfactor=Hfac;
  KKTdiagx=inKKTdiagx;
  KKTdiagy=inKKTdiagy;
  
  int status=0;
  if (precond){
    status=precond->init_system(KKTdiagx,KKTdiagy,Hfac,prec,params);
    if ((status)&&(cb_out()))
      get_out()<<"**** WARNING in QPIterativeKKTSolver::QPinit_KKTdata(): precond->init_data() returned"<<status<<std::endl;
  }

  
  //sol.init(0,0,0.); //start each predictor correcto step anew

  return status;
}

// *************************************************************************
//                             QPsolve_KKTsystem
// *************************************************************************

// solve the KKT System

/// on input
/// dualrhs = -(Qx+A'y+G'modelx+c) +mu (...)
/// primalrhs = -(Ax+s) + mu ()...
int QPIterativeKKTSolver::QPsolve_KKTsystem(Matrix& solx,Matrix& soly,
					 const Matrix& primalrhs,
					 const Matrix& dualrhs,
					 Real rhsmu,
					 Real rhscorr,
					 Real prec,
					 QPSolverParameters* /* params */)
{
  maxit_bnd=max(Integer(100*std::log(primalrhs.rowdim()+dualrhs.rowdim())),min(maxit_bnd,2*(primalrhs.rowdim()+dualrhs.rowdim())));

  //---- prepare the right hand side
  sysrhs.init(dualrhs);
  sysrhs.concat_below(primalrhs);
  Integer start_mod=dualrhs.rowdim()+primalrhs.rowdim();
  int status1=0;
  if (model) {
    sysrhs.enlarge_below(model->dim_model()+model->dim_constraints(),0.);
    Integer start_cstr=start_mod+model->dim_model();
    status1=model->add_localrhs(sysrhs,rhsmu,rhscorr,start_mod,start_cstr,true);
    if (status1){
      if (cb_out()){
	get_out()<<"**** WARNING: QPIterativeKKTSolver::QPsolve_KKTsystem(....): model->add_local_rhs failed and returned "<<status1<<std::endl;
      }
    }
  }

  //---- call the solver
  assert(solver);
  Real termprec=prec*min(1.,norm2(sysrhs));
  Integer ntries=5;
  int status2;
  do {
    solver->set_maxit(maxit_bnd);
    status2=solver->compute(*this,sol,termprec);
    maxit_bnd=max(maxit_bnd,2*solver->get_maxit());
    if ((solver->get_residual_norm()>1.1*solver->get_termprec())&&(solver->get_nmult()>=solver->get_maxit())){
      if (termprec<1e-7*min(1.,norm2(sysrhs)))
	termprec*=10.;
      Matrix tmpvec;
      tmpvec.rand(sol.rowdim(),1);
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
      get_out()<<"**** WARNING: QPIterativeKKTSolver::QPsolve_KKTsystem(....): solver->compute failed and returned "<<status2<<std::endl;
    }
  }

  //---- split up the solution
  solx.init(dualrhs.rowdim(),1,sol.get_store());
  soly.init(primalrhs.rowdim(),1,sol.get_store()+dualrhs.rowdim());

  int status3=0;
  if (model){
    solmod.init(model->dim_model(),1,sol.get_store()+start_mod);
    Integer start_cstr=start_mod+model->dim_model();
    solcstr.init(model->dim_constraints(),1,sol.get_store()+start_cstr);
    status3=model->computed_step(solmod,solcstr);
    if (status3){
      if (cb_out()){
	get_out()<<"**** WARNING: QPIterativeKKTSolver::QPsolve_KKTsystem(....): model->computed_step failed and returned "<<status3<<std::endl;
      }
    }
  }

  return std::abs(status1)+std::abs(status2)+std::abs(status3);
}

// *************************************************************************
//                             ItSys_mult
// *************************************************************************

int QPIterativeKKTSolver::ItSys_mult(const Matrix& in_vec,Matrix& out_vec)
{
  nmult++;
  out_vec.newsize(in_vec.rowdim(),1);
  in_vecx.init(KKTdiagx.rowdim(),1,in_vec.get_store());
    
  //first block of rows belonging to Q
  out_vec=in_vecx;
  out_vec%=KKTdiagx;
  Hp->add_Hx(in_vecx,out_vec,Hfactor);
  
  if (A){
    in_vecy.init(KKTdiagy.rowdim(),1,in_vec.get_store()+KKTdiagx.rowdim());
    genmult(*A,in_vecy,out_vec,1.,1.,1);
  }
    
  if (model){
    in_model.init(model->dim_model(),1,in_vec.get_store()+KKTdiagx.rowdim()+KKTdiagy.rowdim());
    model->B_times(in_model,out_vec,1.,1.,1,0);
  }
    
  //second block of rows belonging to the constraints on the design variables
  if (A) {
    in_vecy%=KKTdiagy;
    genmult(*A,in_vecx,in_vecy,1.,-1.);
    out_vec.concat_below(in_vecy);
  }
    
  //third block of rows belonging to the model block
  if (model){
    out_vec.enlarge_below(model->dim_model()+model->dim_constraints());
    chk_set_init(out_vec,1);
    model->localsys_mult(in_vec,out_vec,
			       KKTdiagx.rowdim()+KKTdiagy.rowdim(),
			       KKTdiagx.rowdim()+KKTdiagy.rowdim()+model->dim_model());
    in_vecy.newsize(model->dim_model(),1); chk_set_init(in_vecy,1);
    model->B_times(in_vecx,in_vecy,1.,0.,0,0);
    mat_xpey(in_vecy.rowdim(),out_vec.get_store()+KKTdiagx.rowdim()+KKTdiagy.rowdim(),in_vecy.get_store());
  }

  Real xnormsqr=mat_ip(KKTdiagx.rowdim(),in_vec.get_store());
  if (xnormsqr>1e-10){
    Real Ritz=mat_ip(KKTdiagx.rowdim(),in_vec.get_store(),out_vec.get_store());
    Ritz/=xnormsqr;
    if (Ritz>blockH_norm)
      blockH_norm=Ritz;
  }
  
  return 0;
}

// *************************************************************************
//                             QPget_condition_number
// *************************************************************************

  Real QPIterativeKKTSolver::QPget_condition_number()
  {
    if (precond)
      return precond->get_condition_number(KKTdiagx,KKTdiagy);
	
    class KKTmat: public CH_Matrix_Classes::Lanczosmatrix{
    private:
      QPIterativeKKTSolver* itsol;
      bool positive;
    public:
      KKTmat(QPIterativeKKTSolver* initsol):itsol(initsol),positive(true)
      {assert(initsol);}

      void set_positive(bool pos)
      {positive=pos;}
      
      ///returns the order of the (virtual) symmetric matrix
      Integer lanczosdim() const
      {return itsol->QPget_system_size();}
  
      ///returns a rough estimate on the number of flops needed by lanczosmult() for a  vector
      Integer lanczosflops() const
      {return lanczosdim()*2*(1+(itsol->A?itsol->A->rowdim():0)+(itsol->model?itsol->model->dim_model():0));}

      ///computes  B = (*this) * A; A and B must not be the same object!
      int lanczosmult(const Matrix& A,Matrix& B) const
      {
	Matrix aa(A.rowdim(),1); chk_set_init(aa,1);
	Matrix bb(A.rowdim(),1); chk_set_init(bb,1);
	B.newsize(A.rowdim(),A.coldim()); chk_set_init(B,1);
	for(Integer j=0;j<A.coldim();j++){
	  mat_xeya(A.rowdim(),aa.get_store(),A.get_store()+j*A.rowdim(),positive?1:-1.);
	  if (itsol->precond){
	    itsol->precond->precond_invG1tran(aa);
	    itsol->ItSys_mult(aa,bb);
	    itsol->precond->precond_invG1(bb);
	  }
	  else {
	    itsol->ItSys_mult(aa,bb);
	  }
	  mat_xey(A.rowdim(),B.get_store()+j*A.rowdim(),bb.get_store());
	}
	return 0;
      }
    };

    KKTmat kkt(this);
    Matrix eigval;
    Matrix eigvec;
    Lanczpol lancz;
    //lancz.set_out(get_out_ptr(),get_print_level());
    //lancz.set_out(&std::cout,10);
    lancz.set_nchebit(10);
    lancz.set_maxiter(10);
    //lancz.set_maxiter(50);
    lancz.set_maxmult(500);
    //lancz.set_maxmult(1000);
    lancz.set_relprec(1.e-4);
    int err=lancz.compute(&kkt,eigval,eigvec,1);
    if ((err)&&(cb_out(10))){
      get_out()<<"**** ERROR in QPIterativeKKTSolver::QPget_condition_number(): lancz.compute returned "<<err<<std::endl;
    }
    Real maxeig=max(eigval);
    if (eigval.dim()==0){
      lancz.get_lanczosvecs(eigval,eigvec);
      maxeig=max(eigval);
      if (cb_out(10)){
	get_out()<<" **** WARNING: eigval.dim==0";
	get_out()<<" Ritzval.dim="<<eigval.dim()<<" maxeig=max(Ritzval)="<<maxeig<<std::endl;
      }
    }
    kkt.set_positive(false);
    eigval.init(0,0,0.);
    eigvec.init(0,0,0.);
    err=lancz.compute(&kkt,eigval,eigvec,1);
    if ((err)&&(cb_out(10))){
      get_out()<<"**** ERROR in QPIterativeKKTSolver::QPget_condition_number(): lancz.compute returned "<<err<<std::endl;
    }
    Real mineig=-max(eigval);
    if (eigval.dim()==0){
      lancz.get_lanczosvecs(eigval,eigvec);
      mineig=-max(eigval);
      if (cb_out(10)){
	get_out()<<" **** WARNING: eigval.dim==0";
	get_out()<<" Ritzval.dim="<<eigval.dim()<<" mineig=-max(Ritzval)="<<mineig<<std::endl;
      }
    }
    Real cond=maxeig/mineig;
    if (cb_out(0)){
      get_out()<<" lanczcond="<<maxeig<<"/"<<mineig<<"="<<cond<<std::endl;
    }

    // //BEGINT test
    // if(kkt.lanczosdim()<500){
    //   kkt.set_positive(true);
    //   Matrix A(kkt.lanczosdim(),kkt.lanczosdim(),0.);
    //   for (Integer i=0;i<A.rowdim();i++)
    // 	A(i,i)=1.;
    //   Matrix B;
    //   kkt.lanczosmult(A,B);
    //   std::cout<<" lanczcond="<<maxeig<<"/"<<mineig<<"="<<cond<<std::endl;
    //   std::cout<<" norm2(B-transpose(B))="<<norm2(B-transpose(B));
    //   Symmatrix S(B);
    //   S.eig(A,B);
    //   std::cout<<" maxeig="<<max(B)<<" mineig="<<min(B)<<" cond="<<max(B)/min(B);
    //   std::cout<<" eigs="<<transpose(B);
    // }
    // //END test
    
    return cond;
  }

}

