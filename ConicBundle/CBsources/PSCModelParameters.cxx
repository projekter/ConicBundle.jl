/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/PSCModelParameters.cxx
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



#include "PSCModelParameters.hxx"

using namespace CH_Matrix_Classes;

namespace ConicBundle {


  PSCModelParametersObject::~PSCModelParametersObject() {
  }

  PSCModelParameters::~PSCModelParameters() {
  }

  // *****************************************************************************
  //                                get_minorant
  // *****************************************************************************

  //Compute the subgradient corresponding to the matrix P*Diag(d)*P'/sum(d)
  //and store it in a new column. P has orthonormal columns and d is 
  //a nonengative vector of appropriate size.

  int PSCModelParametersObject::get_minorant(MinorantPointer& mp,
    Real& mp_coeff,
    const Matrix& P,
    const Matrix& d,
    PSCOracle* oracle,
    Integer modification_id) {
    assert(P.coldim() == d.rowdim());
    assert((!mp.empty()) || (P.coldim() > 0));
    Real sumd = sum(d);

    if ((sumd < eps_Real) && (!mp.empty())) {
      return 0; //nothing to aggregate
    }

    if ((mp_coeff < eps_Real) || (mp.empty())) {
      mp.clear();
      mp_coeff = 0.;
    } else {
      mp.scale(mp_coeff / (sumd + mp_coeff));
    }

    Matrix tmpvec;
    Matrix tmpmat;
    if (sumd > eps_Real) {
      tmpvec.init(d, 1. / (sumd + mp_coeff));
      for (Integer i = 0; i < tmpvec.rowdim(); i++) {
        Real r = tmpvec(i);
        if (r < 0) {
          if (r < -1e-10 * sumd) {
            if (cb_out())
              get_out() << "**** WARNING PSCModelParametersObject::get_minornat(......): negativ d(" << i << ")=" << r << " truncated to 0" << std::endl;
          }
          tmpvec(i) = 0.;
        }
      }
      tmpvec.sqrt();
      tmpmat = P;
      tmpmat.scale_cols(tmpvec);
    } else {
      if (mp.empty())
        tmpmat = P.col(0);
    }

    int err = 0;
    if (tmpmat.coldim() > 0) {
      Minorant* mnrt = oracle->generate_minorant(tmpmat);
      if ((mnrt == 0) || (!mnrt->offset_gives_value_at_origin())) {
        if (cb_out())
          get_out() << "**** ERROR PSCModelParametersObject::get_minorant(......): oracle->generate_minorant failed to generate a minorant or returned one with offset_gives_value_at_origin()==false" << std::endl;
        return 1;
      }
      if (mp.empty())
        mp.init(mnrt, modification_id);
      else
        err = mp.aggregate(MinorantPointer(mnrt, modification_id));
    }

    mp_coeff += sumd;
    return err;
  }




  //******************************************************************************
  //                              select_model
  //******************************************************************************

  int PSCModelParameters::select_model(Matrix& modelvecs,
    MinorantPointer& model_aggregate,
    Matrix& topvecs,
    Matrix& Ritz_values,
    Integer& activedim,
    Integer& keepsize,
    Integer& skippedsize,
    Real primal_Ritzval,
    const Matrix& primaleigs,
    const Matrix& primalvecs,
    const MinorantPointer& primal_aggregate,
    Real primal_aggregate_coeff,
    Real growthrate,
    const Matrix& primalgrowth,
    const Matrix& dualgrowth,
    const Matrix& /* cand_Ritzvec */,
    const Matrix& /* cand_Ritzval */,
    PSCOracle* oracle,
    Integer modification_id,
    FunctionTask function_task,
    Real function_factor,
    BundleModel::ModelUpdate model_update,
    Integer /* center_id */,
    const Matrix& /* center_y */,
    Integer  /* cand_id */,
    const Matrix& cand_y,
    Real model_maxviol,
    Real diffval_center_aggregate,
    BundleProxObject& /* H */) {
    //topvecs holds the Ritz-vectors sorted nonincreasingly by Ritz-value

    if (cb_out(2)) {
      get_out() << " PSCMP";
    }

    Integer vecdim = topvecs.rowdim();
    if (model_update != BundleModel::null_step) {
      n_nullsteps = 0;
    } else {
      n_nullsteps++;
    }

    if ((function_task != ObjectiveFunction) && (primal_Ritzval < 1e-4)) {
      primal_Ritzval = -1e-8 * min(1., function_factor);
    }

    Matrix tmpmat;
    Matrix tmpvec;

    //--- determine an estimate of the rank of the primal solution
    assert(primaleigs.rowdim() == primalgrowth.rowdim());
    Integer primalrank = 0;
    Real growththreshold = min(.95, std::sqrt(growthrate));
    if ((update_rule == 0) && (model_update == BundleModel::descent_step)) {
      primalrank = 0;
      while ((primalrank < primalgrowth.rowdim()) && (primalgrowth(primalrank++) > growththreshold));
    } else {
      primalrank = primalgrowth.rowdim();
      while ((--primalrank >= 0) && (primalgrowth(primalrank) < growththreshold));
      primalrank++;
    }

    if (cb_out(2)) {
      get_out() << " growthr=" << growthrate << " thrsh=" << growththreshold << std::endl;
      get_out().precision(8);
      for (Integer i = 0; i < primaleigs.rowdim(); i++)
        get_out() << std::setw(3) << i << ": " << std::setw(13) << primaleigs(i) << " (" << std::setw(13) << primalgrowth(i) << "," << std::setw(13) << dualgrowth(i) << ")\n";
      if (model_update == BundleModel::descent_step)
        get_out() << " descrank=";
      else
        get_out() << " nullrank=";
      get_out() << primalrank;
    }

    //--- determine a Ritz value above which all Ritz vectors are considered active/relevant
    Real new_cutoffval = primal_Ritzval - .1 * diffval_center_aggregate;

    assert(primalrank <= primalvecs.coldim());
    if (primalvecs.coldim() > 0) {
      oracle->evaluate_projection(cand_y, primalvecs, 1e-10, tmpmat, tmpvec);
      if (cb_out(2)) {
        get_out() << " prim_Ritz=" << transpose(tmpvec);
      }

      new_cutoffval = tmpvec(max(primalrank - 1, 0)) - .01 * diffval_center_aggregate;
    } else {
      new_cutoffval = min(new_cutoffval, Ritz_values(min(Ritz_values.rowdim() - 1, max(primalrank, 0))) - .01 * diffval_center_aggregate);
    }

    if (model_update != BundleModel::null_step) {
      cutoffval = new_cutoffval;
    } else {
      cutoffval = min(.95 * cutoffval + .05 * new_cutoffval, new_cutoffval);
    }

    //Real mingap=1e-4*(primaleigs.rowdim()>0?primaleigs(0):min(1.,function_factor));
    //Real gapsz=max(mingap,0.01*(Ritz_values(0)-primal_Ritzval));

    //--- determine which Ritz vectors should be included in forming the new model
    Integer actind = 0;
    while ((actind < Ritz_values.rowdim()) && (Ritz_values(actind) >= cutoffval))
      actind++;

    Real new_gapsz = Ritz_values(0) - Ritz_values(actind - 1);
    Real new_Ritz_skipval = cutoffval - 5 * max(new_gapsz, model_maxviol / function_factor);
    if (model_update != BundleModel::null_step) {
      gapsz = new_gapsz;
      Ritz_skipval = new_Ritz_skipval;
    } else {
      gapsz = max(new_gapsz, .95 * gapsz + .05 * new_gapsz);
      Ritz_skipval = min(new_Ritz_skipval, .95 * Ritz_skipval + .05 * new_Ritz_skipval);
    }

    if (cb_out(2)) {
      get_out().precision(8);
      get_out() << " cutoffval=" << cutoffval << "(" << new_cutoffval << ")";
      get_out() << " aggrc=" << primal_aggregate_coeff << " actind=" << actind << " primal_Ritzval=" << primal_Ritzval << " gap=" << gapsz << "(" << new_gapsz << ") skipval=" << Ritz_skipval << "(" << new_Ritz_skipval << ")";
    }

    if (model_update == BundleModel::descent_step) {
      //-- reset to active part
      //if (model_update==BundleModel::descent_step){
      activedim = min(Ritz_values.dim(), max(min(actind, activedim + 3), primalrank + 1));
      keepsize = min(primalrank + 1, primalvecs.coldim());
    } else {
      //-- ensure some increase
      if ((function_task != ObjectiveFunction) && (primal_Ritzval < 1e-12) && (Ritz_values(0) < 1e-12)) {
        keepsize = min(max(primalrank + 1, keepsize), primalvecs.coldim());
        activedim = min(Ritz_values.dim(), max(min(actind, activedim + 3), keepsize));
      } else {
        keepsize = min(max(primalrank + 1, keepsize + 1), primalvecs.coldim());
        activedim = min(Ritz_values.dim(), max(max(min(actind, activedim + 3), activedim + 1), keepsize));
      }
    }
    assert(activedim <= Ritz_values.dim());


    if (cb_out(2)) {
      get_out() << " keepsize=" << keepsize << " activedim=" << activedim;
    }

    //---- find the number of vectors to keep on top of activedim for the quadratic term

    Integer skpsz = 0;
    while ((activedim + skpsz < Ritz_values.dim()) && (Ritz_values(activedim + skpsz) > Ritz_skipval)) {
      skpsz++;
    }

    skippedsize = max(min(10, Ritz_values.dim() - activedim), min(20, skpsz));

    topvecs.delete_cols(Range(activedim + skippedsize, topvecs.coldim() - 1));
    Ritz_values.reduce_length(topvecs.coldim());

    if (cb_out(2)) {
      get_out() << " skippedsize=" << skippedsize << "(" << skpsz << ") maxgap=" << Ritz_values(0) - Ritz_values(skippedsize + activedim - 1) << " (lmax=" << Ritz_values(0) << ")";
    }

    //--- select new bundlevectors from primalvecs and topvecs
    //--- form an orthonormal basis of these 

    Integer maxaddcols = activedim;

    tmpmat.newsize(vecdim, keepsize + maxaddcols);
    tmpmat.init(vecdim, keepsize, primalvecs.get_store());

    Integer aggregation_mode = 0;
    // 0: primalvecs[keepsize ...], 1: orthogonalize against modelvecs


    switch (update_rule) {
    case 0: default:
    {
      Real acceptval = cutoffval;
      if (update_rule == 2) {
        acceptval = Ritz_skipval;
      } else if (update_rule == 3) {
        acceptval = cutoffval - 2 * gapsz;
      } else if (update_rule == 4) {
        acceptval = cutoffval - 10 * gapsz;
      } else if (update_rule == 5) {
        acceptval = cutoffval - n_nullsteps * gapsz;
      } else if (update_rule == 6) {
        acceptval = cutoffval - (5 + n_nullsteps) * gapsz;
      }
      if (cb_out(2)) {
        get_out() << " acceptval=" << acceptval;
      }

      //-- include subspace of Ritz vector to maximum Ritz value (~ new subgradient)
      Indexmatrix piv;  //pivotindices for QR of tmpmat
      Integer r = tmpmat.QR_factor(piv);
      modelvecs.init(vecdim, 1, topvecs.get_store());
      r = tmpmat.QR_concat_right(modelvecs, piv, r);
      Real residualtrace = 0.;

      if (r < tmpmat.coldim()) {
        tmpmat.delete_cols(Range(r, tmpmat.coldim() - 1));
        piv.init(Range(0, r - 1));
      } else {
        residualtrace += std::fabs(tmpmat(r - 1, r - 1));
        if (cb_out(2)) {
          get_out() << "\n " << 0 << ": " << std::fabs(tmpmat(r - 1, r - 1));
        }
      }

      //-- add the subspace of the maxaddcols largest Ritz value vectors one by one
      //  if the residual vector is not too small in norm or Ritz value
      for (Integer i = 1; i < maxaddcols; i++) {
        modelvecs.init(vecdim, 1, topvecs.get_store() + i * vecdim);
        r = tmpmat.QR_concat_right(modelvecs, piv, r);

        if (r == tmpmat.coldim()) {
          Real d = std::fabs(tmpmat(r - 1, r - 1));
          residualtrace += d;
          if (cb_out(2)) {
            get_out() << "\n " << i << ": " << d;
          }
          if (d < 1e-8) {
            r--;
            if (cb_out(2)) {
              get_out() << "s";
            }
          } else {
            tmpvec.init(vecdim, 1, 0.);
            tmpvec(r - 1) = 1.;

            tmpmat.Q_times(tmpvec, r);

            Matrix eigvec;
            Matrix eigval;
            oracle->evaluate_projection(cand_y, tmpvec, 1e-10, eigvec, eigval);
            if (cb_out(2)) {
              get_out() << " " << eigval(0);
            }
            if (eigval(0) < acceptval) {
              r--;
              if (cb_out(2)) {
                get_out() << "s";
              }
            }
          }
        } else {
          if (cb_out(2)) {
            get_out() << " [" << i << "] ";
          }
        }
        if (r < tmpmat.coldim()) {
          tmpmat.delete_cols(Range(r, tmpmat.coldim() - 1));
          piv.reduce_length(r);
        }
      }
      if (r == keepsize)
        modelvecs.init(primalvecs.rowdim(), keepsize, primalvecs.get_store());
      else {
        modelvecs.init(tmpmat.rowdim(), r - keepsize, 0.);
        for (Integer i = 0; i < modelvecs.coldim(); i++)
          modelvecs(keepsize + i, i) = 1.;
        tmpmat.Q_times(modelvecs, r);
        swap(tmpmat, modelvecs);
        modelvecs.init(primalvecs.rowdim(), keepsize, primalvecs.get_store());
        modelvecs.concat_right(tmpmat);
      }
      if (cb_out(2)) {
        get_out() << " residualtrace=" << residualtrace;
      }
      break;
    }

    } //end switch(update_rule)


    //--- collect the information of primalvecs not covered by modelvecs
    //--- and aggregate it into the model modelaggregate

    model_aggregate = primal_aggregate;
    Real  model_aggregate_coeff = primal_aggregate_coeff;
    switch (aggregation_mode) {
    case 0: default:
    {
      if (keepsize == primaleigs.dim()) {
        if (model_aggregate_coeff < eps_Real * function_factor) {
          model_aggregate.clear();
          model_aggregate_coeff = 0.;
        }
      } else {
        Indexmatrix ind(Range(keepsize, primaleigs.dim() - 1));
        tmpmat.init(primalvecs.rowdim(), primaleigs.rowdim() - keepsize, primalvecs.get_store() + keepsize * primalvecs.rowdim());
        tmpvec.init(primaleigs.rowdim() - keepsize, 1, primaleigs.get_store() + keepsize);
        if (get_minorant(model_aggregate, model_aggregate_coeff, tmpmat, tmpvec, oracle, modification_id)) {
          if (cb_out())
            get_out() << "**** ERROR PSCModelParameters::select_model(): get_minorant failed for model_aggregate" << std::endl;
        }
      }
      break;
    }
    case 1:
    {
      tmpmat = modelvecs;
      Indexmatrix piv;
      Integer r = tmpmat.QR_factor(piv);
      assert(r == modelvecs.coldim());
      Indexmatrix piv2(piv);
      Integer r2 = tmpmat.QR_concat_right(primalvecs, piv2, r);
      if (r2 == r) {
        if (model_aggregate_coeff < eps_Real * function_factor) {
          model_aggregate.clear();
          model_aggregate_coeff = 0.;
        }
      } else {
        piv2.delete_rows(Range(0, r - 1));
        piv2 -= r;
        Matrix R;
        Symmatrix S;
        Matrix tmpmat2;

        // //TEST BEGIN
        // R=tmpmat(Range(0,r2-1),Range(r,r+primaleigs.rowdim()-1));
        // std::cout<<"r="<<r<<" r2="<<r2<<"  R="<<R;
        // R.triu(-r);
        // std::cout<<"triu(R)="<<R;
        // tmpvec=primaleigs(piv2);
        // std::cout<<" sum(tv)="<<sum(tmpvec)<<" sum(pe)="<<sum(primaleigs);
        // scaledrankadd(R,primaleigs(piv2),S);
        // S.eig(R,tmpvec);
        // tmpmat2.init(tmpmat.rowdim(),tmpvec.rowdim(),0.);
        // for (Integer j=0;j<tmpvec.rowdim();j++){
        //   for (Integer i=0;i<r2;i++){
        //     tmpmat2(i,j)=R(i,j);
        //   }
        // }
        // tmpmat.Q_times(tmpmat2,r2);
        // Real dummy1=primal_aggregate_coeff;
        // std::cout<<" orth="<<norm2(transpose(tmpmat2)*tmpmat2-Diag(Matrix(tmpvec.rowdim(),1,1.)))<<" tr="<<primal_aggregate_coeff+sum(tmpvec);
        // //std::cout<<" primag1=";model_aggregate.display(std::cout);
        // if(get_minorant(model_aggregate,dummy1,tmpmat2,tmpvec,oracle,modification_id)){
        //   if (cb_out())
        //     get_out()<<"**** ERROR PSCModelParameters::select_model(): get_minorant failed for test1 model_aggregate"<<std::endl;
        // }
        // tmpvec.init(cand_y.rowdim(),2,0.);
        // Real val=0.;
        // model_aggregate.get_minorant(val,tmpvec,0,dummy1);
        // std::cout<<" dummy1="<<dummy1<<" val="<<val;
        // model_aggregate=primal_aggregate;
        // Real dummy2=primal_aggregate_coeff;
        // //std::cout<<" primag2=";model_aggregate.display(std::cout);
        // if(get_minorant(model_aggregate,dummy2,primalvecs,primaleigs,oracle,modification_id)){
        //   if (cb_out())
        //     get_out()<<"**** ERROR PSCModelParameters::select_model(): get_minorant failed for test2 model_aggregate"<<std::endl;
        // }
        // model_aggregate.get_minorant(val,tmpvec,1,-dummy2,true);
        // std::cout<<transpose(tmpvec)<<" dummy2="<<dummy2<<" val="<<val;
        // std::cout<<" updatePSCM test1devnorm="<<norm2(tmpvec*Matrix(2,1,1.))<<" valdev="<<std::fabs(val)<<std::endl;
        // model_aggregate=primal_aggregate;
        // //TEST END

        R = tmpmat(Range(0, r2 - 1), Range(r, r + primaleigs.rowdim() - 1));
        R.triu(-r);
        tmpvec = primaleigs(piv2);
        tmpvec.sqrt();
        R.scale_cols(tmpvec);
        //std::cout<<" R="<<R;

        Matrix B(primaleigs.rowdim(), r2 - r); chk_set_init(B, 1);
        for (Integer j = 0; j < tmpvec.rowdim(); j++) {
          for (Integer i = r; i < r2; i++) {
            B(j, i - r) = R(i, j);
          }
        }
        //std::cout<<" B="<<B;

        Indexmatrix Bpiv;
        Integer Br = B.QR_factor(Bpiv);
        B.times_Q(R, Br);
        if (cb_out(2)) {
          get_out() << " r=" << r << " r2=" << r2 << " Br=" << Br << std::endl;
          get_out() << "R*U=" << R;
        }
        R.delete_cols(Range(Br, R.coldim() - 1));
        //std::cout<<" R="<<R;
        rankadd(R, S, 1., 0., 1);
        //std::cout<<" S="<<S<<std::endl;
        assert(Br == S.rowdim());
        S.eig(B, tmpvec);
        genmult(R, B, tmpmat2);
        B = tmpvec;
        assert(min(B) > eps_Real);
        B.sqrt();
        B.inv();
        tmpmat2.scale_cols(B);
        //std::cout<<" orth="<<norm2(transpose(tmpmat2)*tmpmat2-Diag(Matrix(Br,1,1.)))<<std::endl;
        tmpmat2.enlarge_below(tmpmat.rowdim() - r2, 0.);
        tmpmat.Q_times(tmpmat2, r2);
        model_aggregate_coeff = primal_aggregate_coeff;
        if (get_minorant(model_aggregate, model_aggregate_coeff, tmpmat2, tmpvec, oracle, modification_id)) {
          if (cb_out())
            get_out() << "**** ERROR PSCModelParameters::select_model(): get_minorant failed for model_aggregate" << std::endl;
        }
      }
      break;
    }
    }// end switch(aggregation_mode)


    //--- output Ritz_values of the model vectors
    if (cb_out(2)) {
      //get_out()<<" orth="<<norm2(transpose(modelvecs)*modelvecs-Diag(Matrix(modelvecs.coldim(),1,1.)));
      oracle->evaluate_projection(cand_y, modelvecs, 1e-10, tmpmat, tmpvec);
      get_out() << " modelaggrcoeff=" << model_aggregate_coeff << "(" << primal_aggregate_coeff << ")";
      get_out() << " Ritz_values=" << transpose(tmpvec);
    }

    return 0;

  }

} // end namespace ConicBundle





