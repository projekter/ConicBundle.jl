/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/MinorantPointer.cxx
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




#include "MinorantPointer.hxx"
#include "GroundsetModification.hxx"

using namespace CH_Matrix_Classes;

namespace ConicBundle {

// *****************************************************************************
//                                 delete_data()
// *****************************************************************************

void MinorantPointer::delete_data()
{
  assert((md==0)||(md->use_cnt>0));
  if ((md)&&(--(md->use_cnt)==0)){
    delete md;
  }
  md=0;
}

// *****************************************************************************
//                                 new_data()
// *****************************************************************************

  void MinorantPointer::new_data(MinorantUseData* in_md)
{
  delete_data();
  md=in_md;
  if (md){
    md->use_cnt++;
  }
  assert((md==0)||(md->use_cnt>0));
}

// *****************************************************************************
//                              perpare_for_changes
// *****************************************************************************

  int MinorantPointer::prepare_for_changes(double factor,bool with_primal)
{
  assert(md);
  if (md->one_user())
    return md->scale(factor);
  Real sv;
  Minorant* mp;
  if (md->get_scaleval_and_minorant(sv,mp))
    return 1;
  if (mp==0)
    return 1;
  init(mp->clone_minorant(factor*sv,with_primal),md->get_modification_id(),1.);
  return 0;
}


// *****************************************************************************
//                                oeprator=
// *****************************************************************************
  
MinorantPointer& MinorantPointer::operator=(const MinorantPointer& mp)
{
  new_data(mp.md);
  return *this;
}

// *****************************************************************************
//                                  init(MinorantPointer)
// *****************************************************************************
  
  void MinorantPointer::init(const MinorantPointer& mp,Real factor,bool enforce_copy)
{
  if (enforce_copy){
    Real sv;
    Minorant* mnrtp;
    if ((mp.md==0)||(mp.md->get_scaleval_and_minorant(sv,mnrtp)))
      new_data(0);
    else 
      init(mnrtp->clone_minorant(sv*factor,false),mp.md->get_modification_id(),1.);
  }
  else if ((factor==1.)||(mp.md==0)) 
    new_data(mp.md);
  else
    new_data(new MinorantUseData(mp.md,factor));
}


// *****************************************************************************
//                                  init(Minorant*)
// *****************************************************************************

void MinorantPointer::init(Minorant* mp,Integer modification_id,Real factor)
{
  if (mp) {
    new_data(new MinorantUseData(mp,factor,modification_id));
  }
  else{
    delete_data(); 
  }
}
  
// *****************************************************************************
//                    synchronize_ids
// *****************************************************************************

int MinorantPointer::synchronize_ids(CH_Matrix_Classes::Integer new_modification_id,
		       CH_Matrix_Classes::Integer new_center_id,
		       CH_Matrix_Classes::Integer old_center_id,
		       CH_Matrix_Classes::Integer new_cand_id,
		       CH_Matrix_Classes::Integer old_cand_id,
		       CH_Matrix_Classes::Integer new_prex_id)
{
  if (md)
    return md->synchronize_ids(new_modification_id,
			new_center_id,old_center_id,
			new_cand_id,old_cand_id,
			new_prex_id);
  return 0;
}

// *****************************************************************************
//                                 zero
// *****************************************************************************

  bool MinorantPointer::zero() const
  {
    if (empty()) 
      return false;
    Real sv;
    Minorant *mp;
    if (md->get_scaleval_and_minorant(sv,mp))
      return false;
    if (mp==0)
      return false;
    if (sv==0.) 
      return true;
    if (mp->offset()!=0.)
      return false;
    Integer n;
    const Real* cp;
    const Integer* ip;
    if (mp->get_coeffs(n,cp,ip))
      return false;
    if (n==0)
      return true;
    while((--n>=0)&&(*cp==0.));
    return n<0;
  }

// *****************************************************************************
//                                aggregate
// *****************************************************************************

bool MinorantPointer::aggregate() const
{
  if (md)
    return md->aggregate();
  return false;
}

// *****************************************************************************
//                               one_user
// *****************************************************************************

  bool MinorantPointer::one_user() const
  {
    if (md)
      return md->one_user();
    return true;
  }

// *****************************************************************************
//                               get_primal
// *****************************************************************************

  const PrimalData* MinorantPointer::get_primal() 
  {
    if (md==0)
      return 0;
    Real sv;
    Minorant *mp;
    if (md->get_scaleval_and_minorant(sv,mp))
      return 0;
    
    if (sv==1.)
      return mp->get_primal();

    if (mp->get_primal()){
      if (prepare_for_changes(1.,true))
	return 0;
      return md->get_minorant()->get_primal();
    }

    return 0;
  }

// *****************************************************************************
//                                offset
// *****************************************************************************

Real MinorantPointer::offset() const
{
  if (md) 
    return md->offset(); 
  return CB_minus_infinity;
}
  
// *****************************************************************************
//                                coeff
// *****************************************************************************

Real MinorantPointer::coeff(Integer i) const
{
  if (md) 
    return md->coeff(i); 
  return 0.;
}

// *****************************************************************************
//                                coeff
// *****************************************************************************

int MinorantPointer::scale(Real val)
{
  if (empty())
    return 1;
  if (val==1.)
    return 0;
  if (md->one_user()){
    md->scaleval*=val;
    return 0;
  }
  return prepare_for_changes(val,true);
}

// *****************************************************************************
//                             call_primal_extender
// *****************************************************************************

 int MinorantPointer::call_primal_extender(PrimalExtender& prex,Integer prex_id)
 {
   if (md)
     return md->call_primal_extender(prex,prex_id);
   return 1;
 }

// *****************************************************************************
//                              apply_modification
// *****************************************************************************

  int MinorantPointer::apply_modification(const GroundsetModification& gsmdf,Integer mod_id,MinorantExtender* mex,bool apply_costs)
{
  assert((md==0)||(mod_id>=md->get_modification_id()));
  if ((md)&&(md->get_modification_id()==mod_id)){
    //the data has already been updated
    return 0;
  }
  if (!valid()) {
    return 1;
  }
  if (gsmdf.no_modification()&&((!apply_costs)||(gsmdf.get_add_offset()==0.))){
    md->set_modification_id()=mod_id;
    return 0;
  }
  
  md->evals.clear();
  Minorant *mp=md->get_minorant();
  
  //unless apply_costs==true appended values are assumed to be zero and are filled in afterwards
  if (apply_costs){
    mp->add_offset(gsmdf.get_add_offset());
    
    const Matrix* costs=gsmdf.get_append_costs();
    if (costs){
      if (mp->add_coeffs(costs->dim(),costs->get_store(),1.,gsmdf.old_vardim())){
	md->set_modification_id()=-1;
	return 1;
      }
    }
  }
  
  if (gsmdf.map_to_old_variables()){
    if (mp->reassign_coeffs(gsmdf.map_to_old_variables()->dim(),gsmdf.map_to_old_variables()->get_store())){
      md->set_modification_id()=-1;
      return 1;
    }
  }
  
  if ((!apply_costs)&&(gsmdf.appended_vardim()>0)){
    assert((gsmdf.new_var_indices())&&(gsmdf.new_var_indices()->dim()==gsmdf.appended_vardim()));
    if (mex==0){
      md->set_modification_id()=-1;
      return 1;
    }
    assert(md->minorant==mp);
    if (md->scaleval!=1.){
      mp->scale_minorant(md->scaleval);
      md->scaleval=1.;
    }
    int status=mex->extend(*mp,int(gsmdf.appended_vardim()),(const int*)(gsmdf.new_var_indices()->get_store()));
    if (status) {
      md->set_modification_id()=-1;
      return status;
    }
  }
  
  md->set_modification_id()=mod_id;
  return 0;
}
  

// *****************************************************************************
//                               evaluate
// *****************************************************************************

Real MinorantPointer::evaluate(Integer yid,const Matrix& y,bool with_constant) const
{
  if (md)
    return md->evaluate(yid,y,with_constant);
  return CB_minus_infinity;
}
  

// *****************************************************************************
//                              add_offset
// *****************************************************************************

int MinorantPointer::add_offset(Real offset)
{
  if (empty()){
    init(new Minorant(true,offset));
    return 0;
  }
  if (offset==0.)
    return 0;
  if (md->one_user()){
    Real sv;
    Minorant *mp;
    if (md->get_scaleval_and_minorant(sv,mp))
      return 1;
    if (std::fabs(sv)<1e-10*std::fabs(offset)){
      md->scale(1.);
      mp->add_offset(offset);
    }
    else if (mp->add_offset(offset/sv))
      return 1;
    return 0;
  }

  if (prepare_for_changes()){
    return 1;
  }
  Real sv;
  Minorant *mp;
  if (md->get_scaleval_and_minorant(sv,mp))
    return 1;
  if (sv==0.)
    return 1;
  if (mp->add_offset(offset/sv))
    return 1;
  return 0;
}

// *****************************************************************************
//                                coeff
// *****************************************************************************

int MinorantPointer::get_minorant(Real& offset,
				  Matrix& mat,
				  Integer column,
				  Real alpha,
				  bool add,
				  const Indexmatrix* skip_fixed,
				  const Matrix* fixed_vals) const
{
  assert(((skip_fixed==0)&&(fixed_vals==0))||((skip_fixed)&&(fixed_vals)));
  assert((skip_fixed==0)||(skip_fixed->dim()==fixed_vals->dim()));
  
  if (empty())
    return 1;
  Real sv;
  Minorant *mp;
  if (md->get_scaleval_and_minorant(sv,mp))
    return 1;
  if (mp==0)
    return 1;
  Integer n;
  const Real* cp;
  const Integer* ip;
  if (mp->get_coeffs(n,cp,ip))
    return 1;
  Integer rdim=mat.rowdim();
  Real* matp=mat.get_store()+column*rdim;
  alpha*=sv;
  if (alpha==0.){
    if (!add){
      offset=0.;
      mat_xea(rdim,matp,0.);
    }
    return 0;
  } 

  /*
  Integer mnrtdim;
  assert(((mnrtdim=rdim+(skip_fixed?skip_fixed->dim():0)))||(mnrtdim==0));
  assert(((ip!=0)&&(n<=mnrtdim)&&(max(Indexmatrix(n,1,ip)<mnrtdim)))||((ip==0)&&(n<=mnrtdim)));
  Matrix mnrt;
  assert(&(mnrt=ip?Matrix(Sparsemat(mnrtdim,Integer(1),n,ip,Indexmatrix(n,1,Integer(0)).get_store(),cp)):Matrix(Sparsemat(mnrtdim,Integer(1),n,Indexmatrix(Range(0,n-1)).get_store(),Indexmatrix(n,1,Integer(0)).get_store(),cp))));
  Matrix test;
  assert(&(test=((add?mat.col(column):Matrix(rdim,1,0.))+alpha*(skip_fixed?mnrt.delete_rows(*skip_fixed):mnrt)))||(true));
  */

  if (!add){

    //---  initialize
    offset=alpha*mp->offset();
    if ((skip_fixed==0)||(skip_fixed->dim()==0)){
      if (ip==0){
	n=min(n,rdim);
	mat_xeya(n,matp,cp,alpha);
	mat_xea(rdim-n,matp+n,0.);
      }
      else {
	Integer nexti=0;
	for(Integer h=n;--h>=0;){
	  if (*ip>=rdim)
	    break;
	  mat_xea(*ip-nexti,matp+nexti,0.);
	  nexti=*ip+1;
	  *(matp+(*ip++))=(*cp++)*alpha;
	}
	mat_xea(rdim-nexti,matp+nexti,0.);
      }
    }
    else {
      //skip_fixed has elements
      assert(fixed_vals);
      Integer sn=skip_fixed->dim();
      const Integer* sfp=skip_fixed->get_store();
      const Real* fvp=fixed_vals->get_store();
      if (ip==0){
	n=min(n,rdim+sn);
	Integer skipi=0;
	Integer k=0;
	while((k<sn)&&(*sfp<n)){
	  mat_xeya(*sfp-skipi,matp+skipi-k,cp+skipi,alpha);
	  skipi=*sfp+1;
	  offset+=alpha*(*(cp+*sfp++))*(*fvp++);
	  k++;
	}
	mat_xeya(n-skipi,matp+skipi-k,cp+skipi,alpha);
	mat_xea(rdim-(n-k),matp+(n-k),0.);
      }
      else {
	Integer nexti=0;
	Integer h=0;
	Integer k=0;
	while((h<n)&&(k<sn)){
	  if (*sfp<*ip) {
	    k++;
	    sfp++;
	    fvp++;
	    continue;
	  }
	  if (*sfp==*ip) {
	    offset+=alpha*(*cp++)*(*fvp++);
	    sfp++;
	    k++;
	    ip++;
	    h++;
	    continue;
	  }
	  if (*ip-k>=rdim)
	    break;
	  mat_xea(*ip-k-nexti,matp+nexti,0.);
	  nexti=*ip-k+1;
	  *(matp+(*ip++)-k)=(*cp++)*alpha;
	  h++;
	}
	assert((k==sn)||(h==n)||((*sfp>=*ip)&&(*ip-k>=rdim)));
	while((h<n)&&(*ip-k<rdim)){
	  mat_xea(*ip-k-nexti,matp+nexti,0.);
	  nexti=*ip-k+1;
	  *(matp+(*ip++)-k)=(*cp++)*alpha;
	  h++;
	}
	mat_xea(rdim-nexti+k,matp+nexti,0.);
      }
    }
  }

  else {

    //---  add
    offset+=alpha*mp->offset();
    if ((skip_fixed==0)||(skip_fixed->dim()==0)){
      if (ip==0){
	n=min(n,rdim);
	mat_xpeya(n,matp,cp,alpha);
      }
      else {
	for(Integer h=n;--h>=0;){
	  if (*ip>=rdim)
	    break;
	  *(matp+(*ip++))+=(*cp++)*alpha;
	}
      }
    }
    else {
      //skip_fixed has elements
      assert(fixed_vals);
      Integer sn=skip_fixed->dim();
      const Integer* sfp=skip_fixed->get_store();
      const Real* fvp=fixed_vals->get_store();
      if (ip==0){
	n=min(n,rdim+sn);
	Integer skipi=0;
	Integer k=0;
	while((k<sn)&&(*sfp<n)){
	  mat_xpeya(*sfp-skipi,matp+skipi-k,cp+skipi,alpha);
	  skipi=*sfp+1;
	  offset+=alpha*(*(cp+*sfp++))*(*fvp++);
	  k++;
	}
	mat_xpeya(n-skipi,matp+skipi-k,cp+skipi,alpha);
      }
      else {
	Integer h=0;
	Integer k=0;
	while((h<n)&&(k<sn)){
	  if (*sfp<*ip) {
	    k++;
	    sfp++;
	    fvp++;
	    continue;
	  }
	  if (*sfp==*ip) {
	    offset+=alpha*(*cp++)*(*fvp++);
	    sfp++;
	    k++;
	    ip++;
	    h++;
	    continue;
	  }
	  if (*ip>=rdim)
	    break;
	  *(matp+(*ip++)-k) += (*cp++)*alpha;
	  h++;
	}
	assert((k==sn)||(h==n)||((*sfp>=*ip)&&(*ip>=rdim)));
	while((h<n)&&(*ip<rdim)){
	  *(matp+(*ip++)-k) += (*cp++)*alpha;
	  h++;
	}
      }
    }
  }

  /*
  assert(norm2(test-mat.col(column))<1e-10*(1.+norm2(test)));
  */

  return 0;
}
  
  
// *****************************************************************************
//                               get_minorant
// *****************************************************************************

/// if mp is empty, it constructs new MinorantUseData pointing to this minorant 

int MinorantPointer::get_minorant(MinorantPointer& mp,
				  Real alpha) const
{
  if (empty())
    return 1;
  assert(valid());
  if (mp.empty()){
    mp.init(*this,alpha);
    return 0;
  }
  if (mp.empty()){
    mp.init(new Minorant,0,1.);
  }
  else {
    if (mp.prepare_for_changes())
      return 1;
  }
  Real out_sv;
  Minorant *out_mnrt;
  if (mp.md->get_scaleval_and_minorant(out_sv,out_mnrt))
    return 1;
  if (out_sv==0.){
    mp.init(new Minorant,0,1.);
    mp.md->get_scaleval_and_minorant(out_sv,out_mnrt);
    assert(out_sv!=0.);
  }
  assert(mp.valid());
  Real in_sv;
  Minorant *in_mnrt;
  if (md->get_scaleval_and_minorant(in_sv,in_mnrt))
    return 1;
  if (in_mnrt==0)
    return 1;
  alpha*= in_sv/out_sv;
  if (alpha==0.){
    return 0;
  }

  if (aggregate())
    mp.md->aggregated(1);

  out_mnrt->add_offset(in_mnrt->offset()*alpha);
   
  Integer n;
  const Real* cp;
  const Integer* ip;
  if (in_mnrt->get_coeffs(n,cp,ip))
    return 1;

  if (n<=0)
    return 0;  //no coefficients need to be added

  return out_mnrt->add_coeffs(n,cp,ip,alpha);
}


// *****************************************************************************
//                               get_minorant
// *****************************************************************************

/// if mp is empty, it constructs new MinorantUseData pointing to this minorant (cloned if transformed more than by scaling); if it is not empty, and its use_cnt==1 it adds the Minorant, if its use_cnt>1 it clones it first 

int MinorantPointer::get_minorant(MinorantPointer& mp,
				  Real alpha,
				  const Sparsemat* sp,
				  const Indexmatrix* provided_row_indices,
				  const Indexmatrix* needed_col_indices,
				  bool enforce_copy) const
{
  if (empty())
    return 1;
  assert(valid());
  if ((mp.empty())&&(sp==0)){
    mp.init(*this,alpha,enforce_copy);
    return 0;
  }
  if (mp.empty()){
    mp.init(new Minorant,0,1.);
  }
  else {
    if (mp.prepare_for_changes())
      return 1;
  }
  Real out_sv;
  Minorant *out_mnrt;
  if (mp.md->get_scaleval_and_minorant(out_sv,out_mnrt))
    return 1;
  if (out_sv==0.){
    mp.init(new Minorant,0,1.);
    mp.md->get_scaleval_and_minorant(out_sv,out_mnrt);
    assert(out_sv!=0.);
  }
  assert(mp.valid());
  Real in_sv;
  Minorant *in_mnrt;
  if (md->get_scaleval_and_minorant(in_sv,in_mnrt))
    return 1;
  if (in_mnrt==0)
    return 1;
  alpha*= in_sv/out_sv;
  if (alpha==0.){
    return 0;
  }

  if (aggregate())
    mp.md->aggregated(1);

  out_mnrt->add_offset(in_mnrt->offset()*alpha);
   
  Integer n;
  const Real* cp;
  const Integer* ip;
  if (in_mnrt->get_coeffs(n,cp,ip))
    return 1;

  if ((n<=0)
      ||(provided_row_indices && (provided_row_indices->dim()==0))
      ||(needed_col_indices && (needed_col_indices->dim()==0))
      )
    return 0;  //no coefficients need to be added

  //if provided_row_indices are given, reduce the elements to those needed
  Matrix tmpcp;
  Indexmatrix tmpip;
  if (provided_row_indices){
    Integer pdim=provided_row_indices->dim();
    if (ip==0){
      //dense case
      Integer cnt=0;
      tmpcp.newsize(min(pdim,n),1); chk_set_init(tmpcp,1);
      const Integer *ind=provided_row_indices->get_store();
      Real* val=tmpcp.get_store();
      for(;(cnt<pdim)&&(*ind<n);cnt++){
	*val++=*(cp+*ind++);
      }
      if (cnt==0)
	return 0;
      n=cnt;
      tmpcp.reduce_length(n);
      cp=tmpcp.get_store();
      ip=provided_row_indices->get_store();
    }
    else {
      //sparse case
      Integer cnt=0;
      tmpcp.newsize(min(pdim,n),1); chk_set_init(tmpcp,1);
      tmpip.newsize(min(pdim,n),1); chk_set_init(tmpip,1);
      const Integer *pind=provided_row_indices->get_store();
      const Integer* const pind_end=pind+pdim;
      const Integer* const ip_end=ip+n;
      Real* val=tmpcp.get_store();
      Integer* ind=tmpip.get_store();
      while(ip!=ip_end){
	if (*ip==*pind){
	  cnt++;
	  *val++=*cp++;
	  *ind++=*ip++;
	  if (ip==ip_end)
	    break;
	  pind++;
	}
	while((pind!=pind_end)&&(*pind<*ip)){
	  pind++;
	}
	if (pind==pind_end)
	  break;
        while((ip!=ip_end)&&(*ip<*pind)){
	  ip++;
	  cp++;
	}
      }
      if (cnt==0)
	return 0;
      n=cnt;
      tmpcp.reduce_length(n);
      tmpip.reduce_length(n);
      cp=tmpcp.get_store();
      ip=tmpip.get_store();
    }
  }
    
  //transform the remaining coefficients
  if (sp==0){
    return out_mnrt->add_coeffs(n,cp,ip,alpha);
  }


  Integer sp_nc=sp->coldim();
  Real* tp=out_mnrt->get_dense_coeff_store(sp_nc);
  Matrix tmpval(0,1,0.);
  if (tp==0){
    tmpval.init(sp_nc,1,0.);
    tp=tmpval.get_store();
  }
  //Matrix test;
  //Real testn2;
  if (ip==0){
    //assert((testn2=norm2(test=Matrix(Sparsemat(sp->rowdim(),1,n,Indexmatrix(Range(0,n-1)),Indexmatrix(n,1,Integer(0)),Matrix(n,1,cp)))))>=0.);
    //dense case
    const Integer *aip=sp->get_rowindex().get_store();
    const Real *avp=sp->get_rowval().get_store();
    Integer i=sp->get_rowinfo().rowdim();
    const Integer *iip=sp->get_rowinfo().get_store();
    const Integer *nzp=iip+i;
    for(;--i>=0;){
      if (*iip>=n) break;
      Real val=*(cp+(*iip++));
      if (val==0.){
	aip+=*nzp;
	avp+=*nzp++;
	continue;
      }
      val*=alpha;
      for(Integer j=(*nzp++);--j>=0;){
	*(tp+*aip++)+=val*(*avp++);
      }
    }
    //Matrix tmp(tmpval);
    //assert(norm2(genmult(*sp,test,tmp,-alpha,1.,1))>1e-8*testn2);
  }
  else {
    //assert((testn2=norm2(test=Matrix(Sparsemat(sp->rowdim(),1,n,Indexmatrix(n,1,ip),Indexmatrix(n,1,Integer(0)),Matrix(n,1,cp)))))>=0.);
    //sparse case
    const Integer *aip=sp->get_rowindex().get_store();
    const Real *avp=sp->get_rowval().get_store();

    Integer i=sp->get_rowinfo().rowdim();
    const Integer *iip=sp->get_rowinfo().get_store();
    const Integer *nzp=iip+i;

    while ((n>0)&&(i>0)){
      if (*ip<*iip){
	ip++;
	cp++;
	--n;
	continue;
      }
      if (*iip<*ip){
	iip++;
	aip+=*nzp;
	avp+=*nzp++;
	--i;
	continue;
      }
      Real val=alpha* *cp++;
      assert(val!=0.);
      ip++;
      --n;
      for(Integer j=(*nzp++);--j>=0;){
	*(tp+*aip++)+=val*(*avp++);
      }
      iip++;
      --i;
    }
    //Matrix tmp(tmpval);
    //assert(norm2(genmult(*sp,test,tmp,-alpha,1.,1))>1e-8*testn2);
  }

  if ((needed_col_indices)&&(tp==tmpval.get_store())){
    Integer ni_i=needed_col_indices->dim();
    const Integer* ip=needed_col_indices->get_store();
    for (Integer i=0;i<ni_i;i++){
      tmpval(i)=tmpval(*ip++);
    }
    return out_mnrt->add_coeffs(ni_i,(const double*)(tmpval.get_store()),
				(const int *)needed_col_indices->get_store());  
  }

  if (tp==tmpval.get_store())
    return out_mnrt->add_coeffs(sp_nc,(const double*)(tmpval.get_store()));
  
  return 0;
}

// *****************************************************************************
//                              aggregate
// *****************************************************************************

  int MinorantPointer::aggregate(const MinorantBundle& minorants,const Matrix& coeff,Real factor)
{
  assert(Integer(minorants.size())==coeff.dim());
  //assert(min(coeff)>=0.);  //not possible for a PSCBlock in QPModelBlock
  assert(factor>=0.);
  if (minorants.size()==0){
    if (cb_out())
      get_out()<<"**** WARNING: MinorantPointer::aggregate(..): aggregating over zero minorants leads to a zero minorant without primal information ... "<<std::endl;
  }
  int err=0;
  if (factor>0.){
    for (unsigned int i=0;i<minorants.size();i++){
      if (coeff(Integer(i))==0.) 
	continue;
      if (!minorants[i].valid()){
	if (cb_out())
	  get_out()<<"**** ERROR: MinorantPointer::aggregate(..): minorant "<<i<<" is not valid, skipping it"<<std::endl;
	err++;
	continue;
      }
      if (empty()) {
	init(minorants[i],coeff(Integer(i))*factor);
	continue;
      }
      if (aggregate(minorants[i],coeff(Integer(i))*factor)){
	if (cb_out())
	  get_out()<<"**** ERROR: MinorantPointer::aggregate(..): aggregating minorant "<<i<<" failed"<<std::endl;
      }
    }
  }
  if (empty()){
    init(new Minorant,0,1.);
  }
 
  return err;
}

   

// *****************************************************************************
//                              aggregate
// *****************************************************************************

  int MinorantPointer::aggregate(const MinorantPointer& minorant, double itsfactor)
{
  assert((valid())&&(minorant.valid()));
  if (prepare_for_changes(1.,true))
    return 1;

  Real in_sv;
  Minorant *in_mnrt;
  if (minorant.md->get_scaleval_and_minorant(in_sv,in_mnrt))
    return 1;
  itsfactor*=in_sv;
  if (itsfactor==0.){
    return 0;
  }
  
  Real my_sv;
  Minorant *my_mnrt;
  if (md->get_scaleval_and_minorant(my_sv,my_mnrt))
    return 1;
  if (my_sv==0.){
    init(minorant,itsfactor);
    return 0;
  }
  itsfactor/=my_sv;

  int err= my_mnrt->aggregate((const Minorant&)(*in_mnrt),double(itsfactor));
  if (err){
    if (cb_out())
      get_out()<<"**** ERROR: MinorantPointer::aggregate(...): minorant->aggregate failed"<<std::endl;
    err++;
  }

  md->aggregated(1);

  return err;
}


// *****************************************************************************
//                               ip
// *****************************************************************************

Real MinorantPointer::ip(const MinorantPointer& mp,
			 const Indexmatrix* skip_fixed,
			 const Matrix* ipdiag) const
{
  assert(valid());
  
  Real xmult;
  Minorant *xmnrt;
  md->get_scaleval_and_minorant(xmult,xmnrt);
  if (xmult==0.){
    return 0.;
  }
  Real ymult;
  Minorant *ymnrt;
  mp.md->get_scaleval_and_minorant(ymult,ymnrt);
  if (ymult==0.){
    return 0.;
  }
  Integer lenx;
  const Real* xval;
  const Integer* xind;
  xmnrt->get_coeffs(lenx,xval,xind);
  if (lenx==0)
    return 0.;
  
  Integer leny;
  const Real* yval;
  const Integer* yind;
  ymnrt->get_coeffs(leny,yval,yind);
  if (leny==0)
    return 0.;

  const Real* dval=0;
  if (ipdiag)
    dval=ipdiag->get_store();
  
  xmult *= ymult;

  //TEST begin
  bool test=false;
  Real testval=0.;
  Real testeps=1e-10;
  if (test) {
    Integer dim=max(lenx,leny);
    if (yind)
      dim=max(dim,yind[leny-1]+1);
    if (xind)
      dim=max(dim,xind[lenx-1]+1);
    if (skip_fixed)
      dim=max(dim,(*skip_fixed)(skip_fixed->dim()-1)+1);
    if (ipdiag){
      assert(dim<ipdiag->dim());
      dim=ipdiag->dim();
    }
    Real dummy;
    Matrix xvec(dim,1,0.);
    this->get_minorant(dummy,xvec,0,1.);
    Matrix yvec(dim,1,0.);
    mp.get_minorant(dummy,yvec,0,1.);
    if (skip_fixed){
      for(Integer i=0;i<skip_fixed->dim();i++)
	xvec((*skip_fixed)(i))=0.;
    }
    if (ipdiag){
      for(Integer i=0;i<dim;i++)
	testval+=xvec(i)*yvec(i)*(*ipdiag)(i);
    }
    else
      testval=CH_Matrix_Classes::ip(xvec,yvec);
    testeps=testeps*(std::fabs(testval)+1.);
  }
  //TEST end

  
  
  if ((skip_fixed==0)||(skip_fixed->dim()==0)){
    if (xind==0) {
      if (yind==0) {
	Integer len=min(lenx,leny);
	Real sum=mat_ip(len,xval,yval,dval);
	assert((!test)||(std::fabs(sum-testval)<testeps));
	return xmult*sum;
      }
      Real sum=mat_ip_dense_sparse(lenx,xval,leny,yval,yind,dval);
      assert((!test)||(std::fabs(sum-testval)<testeps));
      return xmult*sum;
    }
    if (yind==0){
      Real sum=mat_ip_dense_sparse(leny,yval,lenx,xval,xind,dval);
      assert((!test)||(std::fabs(sum-testval)<testeps));
      return xmult*sum;
    }
    Real sum=mat_ip_sparse_sparse(lenx,xval,xind,leny,yval,yind,dval);
    assert((!test)||(std::fabs(sum-testval)<testeps));
    return xmult*sum;
  }

  Real sum=0;
  if (dval==0) {
    if ((xind==0)&&(yind==0)){
      Integer len=min(lenx,leny);
      const Integer *sind=skip_fixed->get_store();
      const Integer* const send=sind+skip_fixed->dim();
      const Real* const xstart=xval;
      while(sind!=send){
	Integer nexts=*sind++;
	if (len<=nexts) 
	  break;
	const Real* const xend=xstart+nexts;
	while(xval!=xend)
	  sum+=(*xval++)*(*yval++);
	xval++,yval++;
      }
      const Real* const xend=xstart+len;
      while(xval!=xend)
	sum+=(*xval++)*(*yval++);
      assert((!test)||(std::fabs(sum-testval)<testeps));
      return xmult*sum;
    }
    if((xind!=0)&&(yind!=0)){
      const Integer* const xend=xind+lenx;
      const Integer* const yend=yind+leny;
      const Integer* const send=skip_fixed->get_store()+skip_fixed->dim();
      const Integer *sind=skip_fixed->get_store();
      while(xind!=xend){
	while ((yind!=yend)&&(*yind<*xind)){
	  yind++;yval++;
	}
	if (yind==yend) break;
	while ((xind!=xend)&&(*xind<*yind)){
	  xind++;xval++;
	}
	if (xind==xend) break;
	if (*xind==*yind){
	  while((sind!=send)&&(*sind<*xind))
	    sind++;
	  if ((sind!=send)&&(*sind==*xind)){
	    xval++;
	    yval++;
	    xind++;
	    yind++;
	    sind++;
	    continue;
	  }
	  sum+=(*xval++)*(*yval++);
	  xind++;yind++;
	}
      }
      assert((!test)||(std::fabs(sum-testval)<testeps));
      return xmult*sum;
    }
    
    if (yind==0){
      const Real* h=xval;xval=yval;yval=h;
      const Integer* k=xind;xind=yind;yind=k;
      swap(lenx,leny);
    }
    
    const Integer* const yend=yind+leny;
    const Integer *sind=skip_fixed->get_store();
    const Integer* const send=sind+skip_fixed->dim();
    while(sind!=send){
      if (*sind==*yind){
	yval++;
	yind++;
	sind++;
      }
      while ((yind!=yend)&&(*yind<*sind)&&(*yind<lenx)){
	sum+=(*(xval+(*yind++)))*(*yval++);
      }
      if ((yind==yend)||(*yind>=lenx)) 
	break;
      while ((sind!=send)&&(*sind<*yind)){
	sind++;
      }
    }
    assert((!test)||(std::fabs(sum-testval)<testeps));
    return xmult*sum;
  }
  else { //dval!=0
    if ((xind==0)&&(yind==0)){
      Integer len=min(lenx,leny);
      const Integer *sind=skip_fixed->get_store();
      Integer nexts;
      const Real* const xstart=xval;
      Integer sn=skip_fixed->dim();
      while(--sn>=0){
	nexts=*sind++;
	if (len<nexts) 
	  break;
	const Real* const xend=xstart+nexts;
	while(xval!=xend)
	  sum+=(*xval++)*(*yval++)*(*dval++);
	xval++,yval++,dval++;
      }
      if ((sn>=0)&&(xval<xstart+len)){
	const Real* const xend=xstart+len;
	while(xval!=xend)
	  sum+=(*xval++)*(*yval++)*(*dval++);
      }
      assert((!test)||(std::fabs(sum-testval)<testeps));
      return xmult*sum;
    }
    if((xind!=0)&&(yind!=0)){
      const Integer* const xend=xind+lenx;
      const Integer* const yend=yind+leny;
      const Integer* const send=skip_fixed->get_store()+skip_fixed->dim();
      const Integer *sind=skip_fixed->get_store();
      while(xind!=xend){
	while ((yind!=yend)&&(*yind<*xind)){
	  yind++;yval++;
	}
	if (yind==yend) break;
	while ((xind!=xend)&&(*xind<*yind)){
	  xind++;xval++;
	}
	if (xind==xend) break;
	if (*xind==*yind){
	  while((sind!=send)&&(*sind<*xind))
	    sind++;
	  if ((sind!=send)&&(*sind==*xind)){
	    xval++;
	    yval++;
	    xind++;
	    yind++;
	    sind++;
	    continue;
	  }
	  sum+=(*xval++)*(*yval++)*(*(dval+*xind));
	  xind++;yind++;
	}
      }
      assert((!test)||(std::fabs(sum-testval)<testeps));
      return xmult*sum;
    }
    
    if (yind==0){
      const Real* h=xval;xval=yval;yval=h;
      const Integer* k=xind;xind=yind;yind=k;
      swap(lenx,leny);
    }
    
    const Integer* const yend=yind+leny;
    const Integer *sind=skip_fixed->get_store();
    const Integer* const send=sind+skip_fixed->dim();
    while(sind!=send){
      if (*sind==*yind){
	yval++;
	yind++;
	sind++;
      }
      while ((yind!=yend)&&(*yind<*sind)&&(*yind<lenx)){
	sum+=(*(xval+(*yind)))*(*yval++)*(*(dval+*yind));
	yind++;
      }
      if ((yind==yend)||(*yind>=lenx)) 
	break;
      while ((sind!=send)&&(*sind<*yind)){
	sind++;
      }
    }
  }

  assert((!test)||(std::fabs(sum-testval)<testeps));
  return xmult*sum;
}



// *****************************************************************************
//                                coeff
// *****************************************************************************

Real MinorantPointer::ip(const Matrix& vec,
			 const Matrix* ipdiag,
			 Integer m_start) const
{
  assert(valid());
  chk_init(vec);
  
  Real xmult;
  Minorant *xmnrt;
  md->get_scaleval_and_minorant(xmult,xmnrt);
  if (xmult==0.){
    return 0.;
  }
  Integer lenx;
  const Real* xval;
  const Integer* xind;
  xmnrt->get_coeffs(lenx,xval,xind);
  if (lenx<=0)
    return 0.;

  const Real* dval=0;
  if (ipdiag)
    dval=ipdiag->get_store();

  if (xind)
    return xmult*mat_ip_dense_sparse(vec.dim(),vec.get_store()+m_start,lenx,xval,xind,dval);
  lenx=min(lenx,vec.dim());
  return xmult*mat_ip(lenx,xval,vec.get_store()+m_start,dval);
}

// *****************************************************************************
//                                equals
// *****************************************************************************

bool MinorantPointer::equals(const MinorantPointer& mp,Real tol) const
{
  assert(valid()&&mp.valid());
  if (*this==mp)
    return true;

  Real xmult;
  Minorant *xmnrt;
  md->get_scaleval_and_minorant(xmult,xmnrt);
  Real ymult;
  Minorant *ymnrt;
  mp.md->get_scaleval_and_minorant(ymult,ymnrt);
  if ((xmult==0.)&&(ymult==0.))
    return true;

  Real xf=xmult*xmnrt->offset();

  const Real abstol=tol*(1.+std::fabs(xf));

  if (std::fabs(xf-ymult*ymnrt->offset())>abstol)
    return false;

  Integer lenx;
  const Real* xval;
  const Integer* xind;
  xmnrt->get_coeffs(lenx,xval,xind);

  Integer leny;
  const Real* yval;
  const Integer* yind;
  ymnrt->get_coeffs(leny,yval,yind);

  if ((lenx!=leny)||((xind==0)!=(yind==0)))
    return false;

  if (ymult==0.) {
    const Real* const xend=xval+lenx;
    while ((xval!=xend)&&(std::fabs(*xval)<abstol)) 
      xval++;
    if (xval==xend)
      return true;
    return false;
  }

  xmult/=ymult;
  
  assert(lenx==leny);
  if (xind==0){
    assert(yind==0);
    const Real* const xend=xval+lenx;
    while ((xval!=xend)&&(std::fabs(xmult* *xval-*yval++)<abstol)) 
      xval++;
    if (xval==xend)
      return true;
    return false;
  }
  assert(yind!=0);
  const Integer* const yend=yind+leny;
  while ((yind!=yend)&&(*yind==*xind++))
    yind++;
  if (yind!=yend)
    return false;
  
  const Real* const xend=xval+lenx;
  while ((xval!=xend)&&(std::fabs(xmult* *xval-*yval++)<abstol)) 
    xval++;
  if (xval==xend)
    return true;
  return false;
}
  
// *****************************************************************************
//                                norm_squared
// *****************************************************************************

  Real MinorantPointer::norm_squared(const Matrix *D) const
{
  assert(valid());

  Real xmult;
  Minorant *xmnrt;
  md->get_scaleval_and_minorant(xmult,xmnrt);
  if (xmult==0.)
    return 0;
  if (D==0)
    return xmult*xmult*xmnrt->norm_squared();
  Integer lenx;
  const Real* xval;
  const Integer* xind;
  xmnrt->get_coeffs(lenx,xval,xind);

  Real sum=0.;
  assert(lenx<=D->dim());
  if (xind==0){
    const Real* const xend=xval+lenx;
    const Real* dval=D->get_store(); 
    while(xval!=xend){
      assert(*dval>0);
      const Real d=*xval++;
      sum+=d*d*(*dval++);
    }
  }
  else {
    const Integer* const indend=xind+lenx;
    const Real* const dval=D->get_store();
    while (xind!=indend){
      assert(*xind<D->dim());
      const Real d=*xval++;
      sum+=d*d*(*(dval+*xind++));
    }
  }

  return sum*xmult*xmult;
}

// *****************************************************************************
//                                dual_norm_squared
// *****************************************************************************

  Real MinorantPointer::dual_norm_squared(const Matrix *D) const
{
  assert(valid());

  Real xmult;
  Minorant *xmnrt;
  md->get_scaleval_and_minorant(xmult,xmnrt);
  if (xmult==0.)
    return 0;
  if (D==0)
    return xmult*xmult*xmnrt->norm_squared();
  Integer lenx;
  const Real* xval;
  const Integer* xind;
  xmnrt->get_coeffs(lenx,xval,xind);

  Real sum=0.;
  assert(lenx<=D->dim());
  if (xind==0){
    const Real* const xend=xval+lenx;
    const Real* dval=D->get_store(); 
    while(xval!=xend){
      assert(*dval>0);
      const Real d=*xval++;
      sum+=d*d/(*dval++);
    }
  }
  else {
    const Integer* const indend=xind+lenx;
    const Real* const dval=D->get_store();
    while (xind!=indend){
      assert(*xind<D->dim());
      const Real d=*xval++;
      sum+=d*d/(*(dval+*xind++));
    }
  }

  return sum*xmult*xmult;
}

// *****************************************************************************
//                                left_genmult
// *****************************************************************************

  Matrix& MinorantPointer::left_genmult(const Matrix& B,
					Matrix &C,
					Real alpha,
					Real beta,
					int thistrans,
					int btrans,
					Integer thisindex) const
{
  assert(valid());
  assert(thisindex>=0);

  const Integer Bnr=B.rowdim();
  const Integer Bnc=B.coldim();
  const Integer Cnr=C.rowdim();
  const Integer Cnc=C.coldim();
  
  Real xmult;
  Minorant* xmnrt;
  md->get_scaleval_and_minorant(xmult,xmnrt);
 
  alpha*=xmult;

  Integer lenx;
  const Real* xval;
  const Integer* xind;
  xmnrt->get_coeffs(lenx,xval,xind);

  if (thistrans){

    assert(Cnr>thisindex);
    const Real* bp=B.get_store();
    Real* cp=C.get_store()+thisindex;
    
    if (btrans==0){ //B is not transposed

      assert(Cnc==Bnc);
            
      if (beta==0.) { //initialize C  
	//ensure size requirements for
	
	if ((alpha==0.)||(lenx==0)){
	  mat_xea(Cnc,cp,Cnr,0.);
	}
	else{
	  if (xind){
	    mat_xeya(Cnc,cp,Cnr,bp+(*xind++),Bnr,alpha*(*xval++));
	    for (Integer i=1;i<lenx;i++)
	      mat_xpeya(Cnc,cp,Cnr,bp+(*xind++),Bnr,alpha*(*xval++));
	  }
	  else {
	    for(Integer i=0;i<Cnc;i++,cp+=Cnr,bp+=Bnr)
	      *cp = alpha*mat_ip(lenx,bp,xval);
	  }
	}
      }	
      else { //add to C

	if ((alpha==0.)||(lenx==0)){
	  if (beta!=1.)
	    mat_xmultea(Cnc,cp,Cnr,beta);
	}
	else{
	  if (xind){
	    mat_xbpeya(Bnc,cp,Cnr,bp+(*xind++),Bnr,alpha*(*xval++),beta);
	    for (Integer i=1;i<lenx;i++)
	      mat_xpeya(Bnc,cp,Cnr,bp+(*xind++),Bnr,alpha*(*xval++));
	  }
	  else {
	    for(Integer i=0;i<Cnc;i++,cp+=Cnr,bp+=Bnr){
	      *cp *= beta;
	      *cp += alpha*mat_ip(lenx,bp,xval);
	    }
	  }
	}
      }

    }
    else {  //B is transposed
      assert(Cnc==Bnr);

      if (beta==0.) { //initialize C  
	
	if ((alpha==0.)||(lenx==0)){
	  mat_xea(Cnc,cp,Cnr,0.);
	}
	else{
	  if (xind){
	    mat_xeya(Cnc,cp,Cnr,bp+(*xind++)*Bnr,1,alpha*(*xval++));
	    for (Integer i=1;i<lenx;i++)
	      mat_xpeya(Cnc,cp,Cnr,bp+(*xind++)*Bnr,1,alpha*(*xval++));
	  }
	  else {
	    mat_xeya(Cnc,cp,Cnr,bp,1,alpha*(*xval++));
	    bp+=Bnr;
	    for(Integer i=1;i<lenx;i++,bp+=Bnr)
	      mat_xpeya(Cnc,cp,Cnr,bp,1,alpha*(*xval++));
	  }
	}
      }	
      else { //add to C

	if ((alpha==0.)||(lenx==0)){
	  if (beta!=1.)
	    mat_xmultea(Cnc,cp,Cnr,beta);
	}
	else{
	  if (xind){
	    mat_xbpeya(Cnc,cp,Cnr,bp+(*xind++)*Bnr,1,alpha*(*xval++),beta);
	    for (Integer i=1;i<lenx;i++)
	      mat_xpeya(Cnc,cp,Cnr,bp+(*xind++)*Bnr,1,alpha*(*xval++));
	  }
	  else {
	    mat_xbpeya(Cnc,cp,Cnr,bp,1,alpha*(*xval++),beta);
	    bp+=Bnr;
	    for(Integer i=1;i<lenx;i++,bp+=Bnr)
	      mat_xpeya(Cnc,cp,Cnr,bp,1,alpha*(*xval++));
	  }
	}
      }

    }
      
  }

  else { //this is not transposed, outer product mode

    assert(Cnr>=lenx);
    assert((xind==0)||(Cnr>xind[lenx-1]));
 
    Real* cp=C.get_store();
    
    if (btrans==0) { // B is not transposed

      assert(Bnr>thisindex);
      assert(Cnc==Bnc);

      const Real* bp=B.get_store()+thisindex;
      
      if (beta==0.) { //initialize C
	
	if ((alpha==0)||(lenx==0)){
	  mat_xea(Cnr*Cnc,cp,0.);
	}
	else {
	  if (xind){
	    mat_xea(Cnr*Cnc,cp,0.);
	    for (Integer i=0;i<lenx;i++)
	      mat_xeya(Cnc,cp+(*xind++),Cnr,bp,Bnr,alpha*(*xval++));
	  }
	  else {
	    for (Integer i=0;i<lenx;i++,cp++)
	      mat_xeya(Cnc,cp,Cnr,bp,Bnr,alpha*(*xval++));
	    if (Cnr>lenx){
	      for (Integer i=0;i<Cnc;i++,cp+=Cnr)
		mat_xea(Cnr-lenx,cp,0.);
	    }
	  }
	}	
      }

      else { // add to C

	if ((alpha==0)||(lenx==0)){
	  if (beta!=1.)
	    mat_xmultea(Cnr*Cnc,cp,beta);
	}
	else {
	  if (xind){
	    if (beta!=1.)
	      mat_xmultea(Cnr*Cnc,cp,beta);
	    for (Integer i=0;i<lenx;i++)
	      mat_xpeya(Cnc,cp+(*xind++),Cnr,bp,Bnr,alpha*(*xval++));
	  }
	  else {
	    for (Integer i=0;i<lenx;i++,cp++)
	      mat_xbpeya(Cnc,cp,Cnr,bp,Bnr,alpha*(*xval++),beta);
	    if ((beta!=1.)&&(Cnr>lenx)){
	      for (Integer i=0;i<Cnc;i++,cp+=Cnr)
		mat_xmultea(Cnr-lenx,cp,beta);
	    }
	  }
	}	
      }
    }

    else {  // B is transposed

      assert(Bnc>thisindex);
      assert(Cnc==Bnr);

      const Real* bp=B.get_store()+thisindex*Bnr;
      
      if (beta==0.) { //initialize C
	
	if ((alpha==0)||(lenx==0)){
	  mat_xea(Cnr*Cnc,cp,0.);
	}
	else {
	  if (xind){
	    mat_xea(Cnr*Cnc,cp,0.);
	    for (Integer i=0;i<lenx;i++)
	      mat_xeya(Cnc,cp+(*xind++),Cnr,bp,1,alpha*(*xval++));
	  }
	  else {
	    for (Integer i=0;i<lenx;i++,cp++)
	      mat_xeya(Cnc,cp,Cnr,bp,1,alpha*(*xval++));
	    for (Integer i=0;i<Cnc;i++,cp+=Cnr)
	      mat_xea(Cnr-lenx,cp,0.);
	  }
	}	
      }

      else { // add to C

	if ((alpha==0)||(lenx==0)){
	  if (beta!=1.)
	    mat_xmultea(Cnr*Cnc,cp,beta);
	}
	else {
	  if (xind){
	    if (beta!=1.)
	      mat_xmultea(Cnr*Cnc,cp,beta);
	    for (Integer i=0;i<lenx;i++)
	      mat_xpeya(Cnc,cp+(*xind++),Cnr,bp,1,alpha*(*xval++));
	  }
	  else {
	    for (Integer i=0;i<lenx;i++,cp++)
	      mat_xbpeya(Cnc,cp,Cnr,bp,1,alpha*(*xval++),beta);
	    if ((beta!=1.)&&(Cnr>lenx)){
	      for (Integer i=0;i<Cnc;i++,cp+=Cnr)
		mat_xmultea(Cnr-lenx,cp,beta);
	    }
	  }
	}
	
      }
      
    }
  }
    
  return C;
}

// *****************************************************************************
//                                right_genmult
// *****************************************************************************

Matrix& MinorantPointer::right_genmult(const Matrix& A,
				      Matrix &C,
				      Real alpha,
				      Real beta,
				      int atrans,
				      int thistrans,
				      Integer thisindex) const
{
  assert(valid());
  assert(thisindex>=0);
  
  const Integer Anr=A.rowdim();
  const Integer Cnr=C.rowdim();
  const Integer Cnc=C.coldim();

  Real xmult;
  Minorant* xmnrt;
  md->get_scaleval_and_minorant(xmult,xmnrt);
 
  alpha*=xmult;

  Integer lenx;
  const Real* xval;
  const Integer* xind;
  xmnrt->get_coeffs(lenx,xval,xind);



  if (thistrans==0){

    assert(Cnc>thisindex);
    Real* cp=C.get_store()+thisindex*Cnr;

    if (atrans==0) {
      
      assert(Cnr==Anr);

      const Real* ap=A.get_store();
      
      if (beta==0.){ //initialize C
	
	if ((alpha==0)||(lenx==0)){
	  mat_xea(Cnr,cp,0.);
	}
	else {
	  if (xind){
	    mat_xeya(Cnr,cp,ap+(*xind++)*Anr,alpha*(*xval++));
	    for (Integer i=1;i<lenx;i++)
	      mat_xpeya(Cnr,cp,ap+(*xind++)*Anr,alpha*(*xval++));
	  }
	  else {
	    mat_xeya(Cnr,cp,ap,alpha*(*xval++));
	    ap+=Anr;
	    for (Integer i=1;i<lenx;i++,ap+=Anr)
	      mat_xpeya(Cnr,cp,ap,alpha*(*xval++));
	  }
	}
      }

      else { //add to C

	if ((alpha==0)||(lenx==0)){
	  if (beta!=1.)
	    mat_xmultea(Cnr,cp,beta);
	}
	else {
	  if (xind){
	    mat_xbpeya(Cnr,cp,ap+(*xind++)*Anr,alpha*(*xval++),beta);
	    for (Integer i=1;i<lenx;i++)
	      mat_xpeya(Cnr,cp,ap+(*xind++)*Anr,alpha*(*xval++));
	  }
	  else {
	    mat_xbpeya(Cnr,cp,ap,alpha*(*xval++),beta);
	    ap+=Anr;
	    for (Integer i=1;i<lenx;i++,ap+=Anr)
	      mat_xpeya(Cnr,cp,ap,alpha*(*xval++));
	  }
	}
      
      }
    }

    else {  // A is transposed
      
      assert(Cnr==A.coldim());

      const Real* ap=A.get_store();
      
      if (beta==0.){ //initialize C
	
	if ((alpha==0)||(lenx==0)){
	  mat_xea(Cnr,cp,0.);
	}
	else {
	  if (xind){
	    mat_xeya(Cnr,cp,1,ap+(*xind++),Anr,alpha*(*xval++));
	    for (Integer i=1;i<lenx;i++)
	      mat_xpeya(Cnr,cp,1,ap+(*xind++),Anr,alpha*(*xval++));
	  }
	  else {
	    for (Integer i=0;i<Cnr;i++,ap+=Anr)
	      *cp++ =  alpha*mat_ip(lenx,ap,xval);
	  }
	}
      }
      
      else { //add to C

	if ((alpha==0)||(lenx==0)){
	  if (beta!=1.)
	    mat_xmultea(Cnr,cp,beta);
	}
	else {
	  if (xind){
	    mat_xbpeya(Cnr,cp,1,ap+(*xind++),Anr,alpha*(*xval++),beta);
	    for (Integer i=1;i<lenx;i++)
	      mat_xpeya(Cnr,cp,1,ap+(*xind++),Anr,alpha*(*xval++));
	  }
	  else {
	    for (Integer i=0;i<Cnr;i++,ap+=Anr){
	      *cp *= beta;
	      *cp++ +=  alpha*mat_ip(lenx,ap,xval);
	    }
	  }
	}
	
      }
    }
    
  }
  
  else {  //this is transposed, outer product mode

    assert(Cnc>=lenx);
    assert((xind==0)||(Cnc>xind[lenx-1]));

    Real* cp=C.get_store();

    if (atrans==0) {
      assert(Cnr==Anr);
      assert(A.coldim()>thisindex);

      const Real* ap=A.get_store()+thisindex*Anr;

      if (beta==0.) { //initialize C

	if ((alpha==0.)||(lenx==0)){
	  mat_xea(Cnr*Cnc,cp,0.);
	}
	else {
	  if (xind) {
	    mat_xea(Cnr*Cnc,cp,0.);
	    for (Integer i=0;i<lenx;i++)
	      mat_xeya(Cnr,cp+(*xind++)*Cnr,ap,alpha*(*xval++));
	  }
	  else {
	    for (Integer i=0;i<lenx;i++,cp+=Cnr)
	      mat_xeya(Cnr,cp,ap,alpha*(*xval++));
	    mat_xea(Cnr*(Cnc-lenx),cp,0.);
	  }
	}
      }
      else { // add to C
	if ((alpha==0.)||(lenx==0)){
	  if (beta!=1.)
	    mat_xmultea(Cnr*Cnc,cp,beta);
	}
	else {
	  if (xind) {
	    if (beta!=1.)
	      mat_xmultea(Cnr*Cnc,cp,beta);
	    for (Integer i=0;i<lenx;i++)
	      mat_xpeya(Cnr,cp+(*xind++)*Cnr,ap,alpha*(*xval++));
	  }
	  else {
	    for (Integer i=0;i<lenx;i++,cp+=Cnr)
	      mat_xbpeya(Cnr,cp,ap,alpha*(*xval++),beta);
	    if (beta!=1.)
	      mat_xmultea(Cnr*(Cnc-lenx),cp,beta);
	  }
	}
      }

    }

    else { //A is transposed
      assert(Cnr==A.coldim());
      assert(Anr>thisindex);

      const Real* ap=A.get_store()+thisindex;

      if (beta==0.) { //initialize C

	if ((alpha==0.)||(lenx==0)){
	  mat_xea(Cnr*Cnc,cp,0.);
	}
	else {
	  if (xind) {
	    mat_xea(Cnr*Cnc,cp,0.);
	    for (Integer i=0;i<lenx;i++)
	      mat_xeya(Cnr,cp+(*xind++)*Cnr,1,ap,Anr,alpha*(*xval++));
	  }
	  else {
	    for (Integer i=0;i<lenx;i++,cp+=Cnr)
	      mat_xeya(Cnr,cp,1,ap,Anr,alpha*(*xval++));
	    mat_xea(Cnr*(Cnc-lenx),cp,0.);
	  }
	}
      }
      else { // add to C
	if ((alpha==0.)||(lenx==0)){
	  if (beta!=1.)
	    mat_xmultea(Cnr*Cnc,cp,beta);
	}
	else {
	  if (xind) {
	    if (beta!=1.)
	      mat_xmultea(Cnr*Cnc,cp,beta);
	    for (Integer i=0;i<lenx;i++)
	      mat_xpeya(Cnr,cp+(*xind++)*Cnr,1,ap,Anr,alpha*(*xval++));
	  }
	  else {
	    for (Integer i=0;i<lenx;i++,cp+=Cnr)
	      mat_xbpeya(Cnr,cp,1,ap,Anr,alpha*(*xval++),beta);
	    if (beta!=1.)
	      mat_xmultea(Cnr*(Cnc-lenx),cp,beta);
	  }
	}
      }

    }

  }

  return C;
}

// *****************************************************************************
//                                display
// *****************************************************************************

std::ostream& MinorantPointer::display(std::ostream& out,int precision) const
{
  if (empty()){ 
    out<<"[empty]";
  }
  else {
    if (precision<1)
      precision=8;
    out.precision(precision);
 
    Real xmult;
    Minorant* xmnrt;
    md->get_scaleval_and_minorant(xmult,xmnrt);
      
    Integer lenx;
    const Real* xval;
    const Integer* xind;
    xmnrt->get_coeffs(lenx,xval,xind);

    out<<"["<<xmnrt->offset()*xmult<<"; "<<lenx;
    for (Integer i=0;i<lenx;i++){
      out<<"\n";
      out.width(precision+4);
      out<<(*xval++)*xmult;
      if (xind){
	out<<" ";
	out.width(8);
	out<<*xind++;
      }
    }
    out<<" ]\n";
  }

  return out;
}


// *****************************************************************************
//                                genmult Bundle=A
// *****************************************************************************

Matrix& genmult(const MinorantBundle& bundle,
		const Matrix& B,
		Matrix& C,
		Real alpha,
		Real beta,
		int Atrans,
		int Btrans,
		Matrix* Coffset)
{
  //TEST begin
  /*
  Matrix testC(C);
  Matrix testC0;
  if (Coffset)
    testC0.init(*Coffset);
  */
  //TEST end
  
  int err=0;
  if (bundle.size()==0){
    assert(C.coldim()==(Btrans?B.rowdim():B.coldim()));
    assert((Atrans==0)||(C.rowdim()==0));
    if (Atrans==0){
      if (beta==0.) {//initialize to zero
	C.init(C.rowdim(),C.coldim(),0.);
      }
      else
	C*=beta;
      if (Coffset){
	assert(Coffset->rowdim()==1);
	assert(Coffset->coldim()==(Btrans?B.rowdim():B.coldim()));
	if (beta==0.)
	  Coffset->init(1,Coffset->coldim(),0.);
	else
	  *Coffset*=beta;
      }
    }
  }
  else {
    for (unsigned i=0;i<bundle.size();i++){
      bundle[i].left_genmult(B,C,alpha,((Atrans==1)||(i==0))?beta:1.,Atrans,Btrans,Integer(i));
    }
    if (Coffset){
      if ((Atrans)&&
	  (Coffset->rowdim()==Integer(bundle.size()))&&
	  (Coffset->coldim()==1)){
	if (beta==0.)
	  for (unsigned i=0;i<bundle.size();i++)
	    (*Coffset)(Integer(i))=alpha*bundle[i].offset();
	else{
	  for (unsigned i=0;i<bundle.size();i++){
	    (*Coffset)(Integer(i))*=beta;
	    (*Coffset)(Integer(i))+=alpha*bundle[i].offset();
	  }
	}
      }
      else if ((Atrans==0)&&
	       (Coffset->coldim()==(Btrans?B.rowdim():B.coldim()))&&
	       (Coffset->rowdim()==1)){
	Matrix tmpvec(1,Integer(bundle.size())); chk_set_init(tmpvec,1);
	for(unsigned i=0;i<bundle.size();i++)
	  tmpvec(Integer(i))=alpha*bundle[i].offset();
	genmult(tmpvec,B,*Coffset,1.,beta,0,Btrans);
      }
      else {
	err++;
      }
    }
  }

  
  if (err)
    MEmessage(MatrixError(ME_unspec,"genmult(const MinorantBundle& bundle,const Matrix& B,Matrix& C,Real alpha,Real beta,int Atrans,int Btrans,Matrix* Coffset) failed"));
  
  //TEST begin
  /*
  Integer Growdim=(Atrans?(Btrans?B.coldim():B.rowdim()):C.rowdim());
  Matrix G(Growdim,Integer(bundle.size()),0.);
  Matrix G0(Integer(bundle.size()),1,0.);
  for(Integer i=0;i<G.coldim();i++){
    bundle[unsigned(i)].get_minorant(G0(i),G,i);
  }
  genmult(G,B,testC,alpha,beta,Atrans,Btrans);
  Real Cn2=norm2(C-testC);
  if (Cn2>1e-10){
    std::cout<<" Cn2="<<Cn2;
    std::cout<<" diff="<<C-testC<<std::endl;
  }
  if (Coffset){
    if (Atrans==0){
      G0.transpose();
      genmult(G0,B,testC0,alpha,beta,Atrans,Btrans);
    }
    else {
      if (beta==0.)
	testC0.init(G0,alpha);
      else {
	testC0*=beta;
	testC0.xpeya(G0,alpha);
      }
    }
    Real C0n2=norm2(*Coffset-testC0);
    if (C0n2>1e-10){
      std::cout<<" C0n2="<<C0n2;
      std::cout<<" diff="<<*Coffset-testC0<<std::endl;
    }
  }
  */
  //TEST end
  
  return C;
}


// *****************************************************************************
//                                genmult Bundle=B
// *****************************************************************************

  Matrix& genmult(const Matrix& A,
		  const MinorantBundle& bundle,
		  Matrix& C,
		  Real alpha,
		  Real beta,
		  int Atrans,
		  int Btrans,
		  Matrix* Coffset)
  {
    int err=0;
    if (bundle.size()==0){
      assert(C.rowdim()==(Atrans?A.coldim():A.rowdim()));
      assert((Btrans==1)||(C.coldim()==0));
      if (Btrans==1){
	if (beta==0.)//initialize to zero
	  C.init(C.rowdim(),C.coldim(),0.);
	else
	  C*=beta;
	if (Coffset){
	  assert(Coffset->coldim()==1);
	  assert(Coffset->rowdim()==(Atrans?A.coldim():A.rowdim()));
	  if (beta==0.)
	    Coffset->init(Coffset->rowdim(),1,0.);
	  else
	    *Coffset *= beta;
	}
      }
    }
    else {
      for (unsigned i=0;i<bundle.size();i++){
	bundle[i].right_genmult(A,C,alpha,((Btrans==0)||(i==0))?beta:1.,Atrans,Btrans,Integer(i));
      }
      if (Coffset){
	if ((Btrans==0)&&
	    (Coffset->rowdim()==1)&&
	    (Coffset->coldim()==Integer(bundle.size()))){
	  if (beta==0.)
	    for (unsigned i=0;i<bundle.size();i++)
	      (*Coffset)(Integer(i))=alpha*bundle[i].offset();
	  else{
	    for (unsigned i=0;i<bundle.size();i++){
	      (*Coffset)(Integer(i))*=beta;
	      (*Coffset)(Integer(i))+=alpha*bundle[i].offset();
	    }
	  }
	}
	else if ((Btrans)&&
		 (Coffset->rowdim()==(Atrans?A.coldim():A.rowdim()))&&
		 (Coffset->coldim()==1)){
	  Matrix tmpvec(Integer(bundle.size()),1); chk_set_init(tmpvec,1);
	  for(unsigned i=0;i<bundle.size();i++)
	    tmpvec(Integer(i))=alpha*bundle[i].offset();
	  genmult(A,tmpvec,*Coffset,1.,beta,Atrans,0);
	}
	else {
	  err++;
	}
      }
    }
    
    if (err)
      MEmessage(MatrixError(ME_unspec,"genmult(const Matrix& A,const MinorantBundle& bundle,Matrix& C,Real alpha,Real beta,int Atrans,int Btrans,Matrix* Coffset) failed"));
    return C;
  }

}


