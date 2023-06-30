/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/SOCModelParameters.cxx
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



#include "SOCModelParameters.hxx"

using namespace CH_Matrix_Classes;

namespace ConicBundle {


SOCModelParametersObject::~SOCModelParametersObject()
{}
 
SOCModelParameters::~SOCModelParameters()
{}
  
//******************************************************************************
//                              select_model
//******************************************************************************

int SOCModelParameters::select_model(Matrix& modelvecs,
				     const Matrix& aggrvec,
				     Real /* cand_SOCval */,
				     const Matrix& cand_SOCvec,
				     Real /* center_SOCval */,
				     const Matrix& /* center_SOCvec */,
				     const Matrix& SOCvecs,
				     SOCOracle* /* oracle */,
				     FunctionTask /* function_task */,
				     Real /* function_factor */,
				     BundleModel::ModelUpdate /* model_update */,
				     Integer /* center_id */,
				     const Matrix& /* center_y */,
				     Integer /* y_id */,
				     const Matrix& /* y */,
				     Real /* model_maxviol */,
				     BundleProxObject& /* H */)
{
  assert((aggrvec.coldim()==1)||(modelvecs.coldim()==0));
  assert((SOCvecs.coldim()==0)||(cand_SOCvec.rowdim()==SOCvecs.rowdim()+1));

  Integer SOCdim=cand_SOCvec.rowdim();
  Integer max_columns=max_model_size>2?max_model_size-1:max(2,SOCdim-1);

  Matrix tmpvec;
  modelvecs.init(0,0,0.);
  Real tol=1e-10;

  //-- insert the aggregate
  if (aggrvec.coldim()==1) {
    tmpvec.init(aggrvec.rowdim()-1,1,aggrvec.get_store()+1);
    Real n2=norm2(tmpvec);
    if (n2<tol){
      tmpvec.init(aggrvec.rowdim()-1,1,0.);
      tmpvec(0)=1.;
    }
    else
      tmpvec/=n2;
    modelvecs.concat_right(tmpvec);
  }
  
  //-- insert the candidate
  tmpvec.init(cand_SOCvec.rowdim()-1,1,cand_SOCvec.get_store()+1);
  if (modelvecs.coldim()==1){
    //orthogonalize to the aggregate
    tmpvec.xpeya(modelvecs,-ip(tmpvec,modelvecs));
  } 
  Real n2=norm2(tmpvec);
  if (n2<tol){
    tmpvec.init(cand_SOCvec.rowdim()-1,1,0.);
    tmpvec(0)=1.;
  }
  else
    tmpvec/=n2;
  modelvecs.concat_right(tmpvec);
  
  assert(modelvecs.coldim()>0);

  //add further columns if bundle_parameters suggest this
  if ((max_columns>modelvecs.coldim())&&(SOCvecs.coldim()>1)){
    // number of vectors to add
    Integer k=max(max_columns-modelvecs.coldim(),SOCvecs.coldim()-1);
    //orthogonalize the vectors to the existing ones
    Matrix tmpmat=SOCvecs;
    tmpmat-=modelvecs*genmult(modelvecs,SOCvecs,tmpvec,1.,0.,1);
    //compute basis to k largest singular values
    Symmatrix S;
    if (tmpmat.rowdim()<=tmpmat.coldim()){
      rankadd(tmpmat,S);
      S.eig(tmpmat,tmpvec,false);
      Integer i=0;
      while ((i<k+1)&&(tmpvec(i)>tol))
	i++;
      tmpmat.delete_cols(Range(i,tmpmat.coldim()-1),true);
      modelvecs.concat_right(tmpmat);
    }
    else {
      rankadd(tmpmat,S,1.,0.,1);
      Matrix P;
      S.eig(P,tmpvec,false);
      Integer i=0;
      while ((i<k+1)&&(tmpvec(i)>tol))
	i++;
      P.delete_cols(Range(i,P.coldim()-1));
      genmult(tmpmat,P,tmpvec);
      modelvecs.concat_right(tmpvec);
    }
    //orthogonalize once more
    modelvecs.QR_factor();
    tmpmat.init(modelvecs.rowdim(),min(modelvecs.coldim(),max(2,max_columns)),0.);
    for (Integer i=0;i<tmpmat.coldim();i++)
      tmpmat(i,i)=1.;
    modelvecs.Q_times(tmpmat,modelvecs.coldim());
    swap(modelvecs,tmpmat);
  }

  return 0;
}




} // end namespace ConicBundle
