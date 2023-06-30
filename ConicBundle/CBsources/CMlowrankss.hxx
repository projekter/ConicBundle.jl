/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/CMlowrankss.hxx
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


#ifndef CONICBUNDLE_CMLOWRANKSS_HXX
#define CONICBUNDLE_CMLOWRANKSS_HXX

/**  @file CMlowrankss.hxx
    @brief Header declaring the class ConicBundle::CMlowrankss (needed for ConicBundle::AffineMatrixFunction)
    @version 1.0
    @date 2015-03-24
    @author Christoph Helmberg
*/

#include "Coeffmat.hxx"
#include "CMlowrankdd.hxx"

namespace ConicBundle {

/** @ingroup implemented_psc_oracle
 */
//@{

 /**@brief implements a low rank matrix \f$AB^T+BA^T\f$ as Coeffmat with \f$A,B\f$ each a sparse rectangular CH_Matrix_Classes::Sparsemat (for use with MatrixSDPfunction, see \ref implemented_psc_oracle).
  */

  class CMlowrankss: public Coeffmat
  {
  private:
    CH_Matrix_Classes::Sparsemat A; ///< this is A in A*B^T+B*A^T
    CH_Matrix_Classes::Sparsemat B; ///< this is B in A*B^T+B*A^T
  public:
    ///copy Ain, Bin and store the user information
    CMlowrankss(const CH_Matrix_Classes::Sparsemat& Ain,const CH_Matrix_Classes::Sparsemat& Bin,CoeffmatInfo* cip=0)
    {A=Ain;B=Bin;CM_type=CM_lowrankss;infop=cip;}
    ///
    virtual ~CMlowrankss(){}
    
    ///makes an explicit copy of itself and returns a pointer to it 
    virtual Coeffmat* clone() const 
    {return new CMlowrankss(A,B,ConicBundle::clone(infop));}

    ///returns the order of the represented symmetric matrix
    virtual CH_Matrix_Classes::Integer dim() const { return A.rowdim(); }
    
    ///returns the value of the matrix element (i,j)
    virtual CH_Matrix_Classes::Real operator()(CH_Matrix_Classes::Integer i,CH_Matrix_Classes::Integer j) const 
    { return CH_Matrix_Classes::ip(A.row(i),B.row(j))+ CH_Matrix_Classes::ip(B.row(i),A.row(j));}

    ///returns a dense symmetric constraint matrix
    virtual void make_symmatrix(CH_Matrix_Classes::Symmatrix& S) const
    { CH_Matrix_Classes::rank2add(A,CH_Matrix_Classes::Matrix(B),S,2.); }

    ///returns the Frobenius norm of the matrix
    virtual CH_Matrix_Classes::Real norm(void) const
    { CH_Matrix_Classes::Matrix C,D; CH_Matrix_Classes::genmult(A,CH_Matrix_Classes::Matrix(B),C,1.,0.,1); CH_Matrix_Classes::genmult(C,C,D);
      CH_Matrix_Classes::Real d=2.*CH_Matrix_Classes::trace(D); CH_Matrix_Classes::genmult(A,A,C,1.,0.,1); CH_Matrix_Classes::genmult(B,B,D,1.,0.,1);
      return sqrt(2.*CH_Matrix_Classes::ip(C,D)+d);}

     ///delivers a new object on the heap corresponding to the matrix P^T(*this)P, the caller is responsible for deleting the object
    virtual Coeffmat* subspace(const CH_Matrix_Classes::Matrix& P) const
    {CH_Matrix_Classes::Matrix C,D; CH_Matrix_Classes::genmult(P,A,C,1.,0.,1); CH_Matrix_Classes::genmult(P,B,D,1.,0.,1);
      return new CMlowrankdd(C,D,ConicBundle::clone(infop)); }
    
    ///multiply constraint permanentely by d; this is to allow scaling or sign changes in the constraints
    virtual void multiply(CH_Matrix_Classes::Real d)
    { A*=d; if (infop) infop->multiply(d);}

    ///returns ip(*this,S)=trace(*this*S), the trace inner product
    virtual CH_Matrix_Classes::Real ip(const CH_Matrix_Classes::Symmatrix& S) const
    { CH_Matrix_Classes::Matrix C; return 2.*CH_Matrix_Classes::ip(CH_Matrix_Classes::genmult(S,A,C),B); }
    
    ///returns ip(*this,PP^T)=trace P^T(*this)P
    virtual CH_Matrix_Classes::Real gramip(const CH_Matrix_Classes::Matrix& P) const
    { CH_Matrix_Classes::Matrix C,D; CH_Matrix_Classes::genmult(P,A,C,1.,0.,1); CH_Matrix_Classes::genmult(P,B,D,1.,0.,1);
      return 2.*CH_Matrix_Classes::ip(C,D); }

    ///returns ip(*this,QQ^T)=trace Q^T(*this)Q for Q=P.rows(start_row,start_row+dim-1) 
    virtual CH_Matrix_Classes::Real gramip(const CH_Matrix_Classes::Matrix& P,CH_Matrix_Classes::Integer start_row,const CH_Matrix_Classes::Matrix* Lam=0) const
    { 
      CH_Matrix_Classes::Matrix C(A.coldim(),P.coldim(),0.);
      CH_Matrix_Classes::Real *bp=C.get_store();
      const CH_Matrix_Classes::Integer brd=C.rowdim();
      const CH_Matrix_Classes::Integer pcd=P.coldim();
      const CH_Matrix_Classes::Integer prd=P.rowdim();
      const CH_Matrix_Classes::Real *pp=P.get_store()+start_row;
      const CH_Matrix_Classes::Indexmatrix &colinfo=A.get_colinfo();
      const CH_Matrix_Classes::Real *vp=(A.get_colval()).get_store();
      const CH_Matrix_Classes::Integer *ip=(A.get_colindex()).get_store();
      for(CH_Matrix_Classes::Integer j=0;j<colinfo.rowdim();j++){
	CH_Matrix_Classes::Integer indj=colinfo(j,0);
	for(CH_Matrix_Classes::Integer i=colinfo(j,1);--i>=0;){
	  CH_Matrix_Classes::mat_xpeya(pcd,bp+indj,brd,pp+(*ip++),prd,(*vp++));
	}
      }
      CH_Matrix_Classes::Matrix D(B.coldim(),P.coldim(),0.);
      bp=D.get_store();
      const CH_Matrix_Classes::Indexmatrix &colinfo2=B.get_colinfo();
      vp=(B.get_colval()).get_store();
      ip=(B.get_colindex()).get_store();
      {for(CH_Matrix_Classes::Integer j=0;j<colinfo2.rowdim();j++){
	CH_Matrix_Classes::Integer indj=colinfo2(j,0);
	for(CH_Matrix_Classes::Integer i=colinfo2(j,1);--i>=0;){
	  CH_Matrix_Classes::mat_xpeya(pcd,bp+indj,brd,pp+(*ip++),prd,(*vp++));
	}
      }}
      CH_Matrix_Classes::Real trval=0;
      if (Lam==0) {
	trval=CH_Matrix_Classes::ip(C,D);
      }
      else {
	assert(Lam->dim()==P.coldim());
	const CH_Matrix_Classes::Real *cp=C.get_store();
	const CH_Matrix_Classes::Real *dp=D.get_store();
	const CH_Matrix_Classes::Real *lp=Lam->get_store();
	for (CH_Matrix_Classes::Integer i=0;i<C.coldim();i++,cp+=C.rowdim(),dp+=D.rowdim())
	  trval+=(*lp++)*CH_Matrix_Classes::mat_ip(C.rowdim(),cp,dp);
      }
      return 2.*trval; 
    }

    ///computes S+=d*(*this);
    virtual void addmeto(CH_Matrix_Classes::Symmatrix& S,CH_Matrix_Classes::Real d=1.) const
    { CH_Matrix_Classes::rank2add(A,CH_Matrix_Classes::Matrix(B),S,2.*d,1.);}

    ///computes D+=d*(*this)*C
    virtual void addprodto(CH_Matrix_Classes::Matrix& D,const CH_Matrix_Classes::Matrix&C ,CH_Matrix_Classes::Real d=1.) const
    {CH_Matrix_Classes::Matrix E; CH_Matrix_Classes::genmult(A,CH_Matrix_Classes::genmult(B,C,E,1.,0.,1),D,d,1.);
     CH_Matrix_Classes::genmult(B,CH_Matrix_Classes::genmult(A,C,E,1.,0.,1),D,d,1.);}

    ///computes D+=d*(*this)*C
    virtual void addprodto(CH_Matrix_Classes::Matrix& D,const CH_Matrix_Classes::Sparsemat&C ,CH_Matrix_Classes::Real d=1.) const
    {CH_Matrix_Classes::Matrix E; CH_Matrix_Classes::genmult(A,CH_Matrix_Classes::genmult(B,C,E,1.,0.,1),D,d,1.);
     CH_Matrix_Classes::genmult(B,CH_Matrix_Classes::genmult(A,C,E,1.,0.,1),D,d,1.);}

    /// computes R=P^T*(*this)*Q
    virtual void left_right_prod(const CH_Matrix_Classes::Matrix& P,const CH_Matrix_Classes::Matrix& Q,CH_Matrix_Classes::Matrix& R) const
    {
      CH_Matrix_Classes::Matrix tmp1; CH_Matrix_Classes::genmult(P,A,tmp1,1.,0.,1,0);
      CH_Matrix_Classes::Matrix tmp2; CH_Matrix_Classes::genmult(B,Q,tmp2,1.,0.,1,0);
      CH_Matrix_Classes::genmult(tmp1,tmp2,R,1.,0.,0,0);
      CH_Matrix_Classes::genmult(P,B,tmp1,1.,0.,1,0);
      CH_Matrix_Classes::genmult(A,Q,tmp2,1.,0.,1,0);
      CH_Matrix_Classes::genmult(tmp1,tmp2,R,1.,1.,0,0);
    }

    ///returns an estimate of number of flops to compute addprodto for a vector
    virtual CH_Matrix_Classes::Integer prodvec_flops() const 
    { return 4*A.nonzeros()+4*B.nonzeros(); }

    ///returns 1 if its structure is as bad as its dense symmetric representation, otherwise 0
    virtual int dense() const
    {return 0;}
    
    ///returns 0 if not sparse, otherwise 1
    virtual int sparse() const
    { return 0;}
    
    /// returns 0 if not sparse. If it is sparse it returns 1 and the nonzero structure in I,J and val, where val is multiplied by d. Only the upper triangle (including diagonal) is delivered
    virtual int sparse(CH_Matrix_Classes::Indexmatrix& /* I */,
		       CH_Matrix_Classes::Indexmatrix& /* J */,
		       CH_Matrix_Classes::Matrix& /* val */,
		       CH_Matrix_Classes::Real /* d=1. */)const
    {return 0;}
    
    /// returns 0 if the support of the costraint matrix is not contained in the support of the sparse symmetric matrix S, 1 if it is contained.
    virtual int support_in(const CH_Matrix_Classes::Sparsesym& /* S */) const
    {return 0;}

    ///returns the inner product of the constraint matrix with S
    virtual CH_Matrix_Classes::Real ip(const CH_Matrix_Classes::Sparsesym& S) const
    {return 2.*CH_Matrix_Classes::ip(A,S*B);}
    
    ///computes S=P^T*(*this)*P
    virtual void project(CH_Matrix_Classes::Symmatrix& S,const CH_Matrix_Classes::Matrix& P) const
    {CH_Matrix_Classes::Matrix C,D; CH_Matrix_Classes::genmult(P,A,C,1.,0.,1); CH_Matrix_Classes::genmult(P,B,D,1.,0.,1);
     CH_Matrix_Classes::rank2add(C,D,S,2.);}
    // S=P^t*(*this)*P

    ///computes S+=Q^T(*this)Q for Q=P.rows(start_row,start_row+dim-1) 
    virtual void add_projection(CH_Matrix_Classes::Symmatrix& S,const CH_Matrix_Classes::Matrix& P,CH_Matrix_Classes::Real alpha=1.,CH_Matrix_Classes::Integer start_row=0) const
    {
      CH_Matrix_Classes::Matrix C(A.coldim(),P.coldim(),0.);
      CH_Matrix_Classes::Real *bp=C.get_store();
      const CH_Matrix_Classes::Integer brd=C.rowdim();
      const CH_Matrix_Classes::Integer pcd=P.coldim();
      const CH_Matrix_Classes::Integer prd=P.rowdim();
      const CH_Matrix_Classes::Real *pp=P.get_store()+start_row;
      const CH_Matrix_Classes::Indexmatrix &colinfo=A.get_colinfo();
      const CH_Matrix_Classes::Real *vp=(A.get_colval()).get_store();
      const CH_Matrix_Classes::Integer *ip=(A.get_colindex()).get_store();
      for(CH_Matrix_Classes::Integer j=0;j<colinfo.rowdim();j++){
	CH_Matrix_Classes::Integer indj=colinfo(j,0);
	for(CH_Matrix_Classes::Integer i=colinfo(j,1);--i>=0;){
	  CH_Matrix_Classes::mat_xpeya(pcd,bp+indj,brd,pp+(*ip++),prd,(*vp++));
	}
      }
      CH_Matrix_Classes::Matrix D(B.coldim(),P.coldim(),0.);
      bp=D.get_store();
      const CH_Matrix_Classes::Indexmatrix &colinfo2=B.get_colinfo();
      vp=(B.get_colval()).get_store();
      ip=(B.get_colindex()).get_store();
      {for(CH_Matrix_Classes::Integer j=0;j<colinfo2.rowdim();j++){
	CH_Matrix_Classes::Integer indj=colinfo2(j,0);
	for(CH_Matrix_Classes::Integer i=colinfo2(j,1);--i>=0;){
	  CH_Matrix_Classes::mat_xpeya(pcd,bp+indj,brd,pp+(*ip++),prd,(*vp++));
	}
      }}
      CH_Matrix_Classes::rank2add(C,D,S,2.*alpha,1.,1); 
    }

    ///computes C= alpha*(*this)*D^(T if dtrans) + beta*C, C is also returned
    virtual const CH_Matrix_Classes::Matrix& postgenmult(const CH_Matrix_Classes::Matrix& D,CH_Matrix_Classes::Matrix& C,
			     CH_Matrix_Classes::Real alpha=1.,CH_Matrix_Classes::Real beta=0.,int dtrans=0) const
    { 
      CH_Matrix_Classes::Matrix E; 
      CH_Matrix_Classes::genmult(A,CH_Matrix_Classes::genmult(B,D,E,1.,0.,1,dtrans),C,alpha,beta);
      return CH_Matrix_Classes::genmult(B,CH_Matrix_Classes::genmult(A,D,E,1.,0.,1,dtrans),C,alpha,1.);
    }

    ///computes C= alpha*D^(T if dtrans)*(*this) + beta*C, C is also returned
    virtual const CH_Matrix_Classes::Matrix& pregenmult(const CH_Matrix_Classes::Matrix& D,CH_Matrix_Classes::Matrix& C,
			     CH_Matrix_Classes::Real alpha=1.,CH_Matrix_Classes::Real beta=0.,int dtrans=0) const
    { 
      CH_Matrix_Classes::Matrix E; 
      CH_Matrix_Classes::genmult(CH_Matrix_Classes::genmult(D,A,E,1.,0.,dtrans),B,C,alpha,beta,0,1);
      return CH_Matrix_Classes::genmult(CH_Matrix_Classes::genmult(D,B,E,1.,0.,dtrans),A,C,alpha,1.,0,1);
    }

    ///returns 1, if p is the same derived class and entries differ by less than tol, otherwise zero
    virtual int equal(const Coeffmat* p,double tol=1e-6) const
    {
      const CMlowrankss *pp=dynamic_cast<const CMlowrankss *>(p);
      if (pp==0) 
	return 0;
      if (CH_Matrix_Classes::equal(A,pp->A,tol) && CH_Matrix_Classes::equal(B,pp->B,tol))
	return 1;
      if (CH_Matrix_Classes::equal(A,pp->B,tol) && CH_Matrix_Classes::equal(B,pp->A,tol))
	return 1;
      return 0;
    }

    ///display constraint information
    virtual std::ostream& display(std::ostream& o) const 
    {o<<"CMlowrankss\n";A.display(o);B.display(o);return o;}
    
    ///put entire contents onto ostream with the class type in the beginning so that the derived class can be recognized by in().
    virtual std::ostream& out(std::ostream& o) const
    {return o<<"LOWRANK_SPARSE_SPARSE\n"<<A<<B;}

    ///counterpart to out(), does not read the class type, though. This is assumed to have been read in order to generate the correct class
    virtual std::istream& in(std::istream& i)
    {return i>>A>>B;     
      if((A.rowdim()!=B.rowdim())||(A.coldim()!=B.coldim())){
         i.clear(std::ios::failbit);
         if (CH_Matrix_Classes::materrout) (*CH_Matrix_Classes::materrout)<<"*** ERROR: CMlowrankss::in(): dimensions of A and B do not match"<<std::endl;
      }
      return i;
    }

    /// constructor with istream and possibly additional user information
    CMlowrankss(std::istream& is,CoeffmatInfo* cip=0)
    {CM_type=CM_lowrankss;infop=cip;in(is);}

  };

  //@}

}

#endif

