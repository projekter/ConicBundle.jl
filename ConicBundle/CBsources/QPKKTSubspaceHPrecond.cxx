/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/QPKKTSubspaceHPrecond.cxx
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
#include "QPKKTSubspaceHPrecond.hxx"

using namespace CH_Matrix_Classes;

namespace ConicBundle {

  QPKKTSubspaceHPrecond::QPKKTSubspaceHPrecond(CH_Matrix_Classes::Integer inmethod, CBout* cb, int cbinc) :
    CBout(cb, cbinc) {
    method = inmethod;
  }

  QPKKTSubspaceHPrecond::~QPKKTSubspaceHPrecond() {
  }

  // *************************************************************************
  //                                init_data
  // *************************************************************************

  int QPKKTSubspaceHPrecond::init_data(QPSolverProxObject* in_Hp,
    QPModelBlockObject* in_model,
    const Sparsemat* in_A,
    const Indexmatrix* in_eq_indices,
    bool in_SchurComplAineq) {
    assert(in_Hp);
    clear();
    Hp = in_Hp;
    model = in_model;
    A = in_A;
    eq_indices = in_eq_indices;
    SchurComplAineq = in_SchurComplAineq;

    Hp->get_precond(diagH, Vp);
    subspace.init(0, 0, 0.);
    last_nmult = 0;
    max_sigma = -1.;

    diaginvval = -1.;

    return 0;
  }


  // *************************************************************************
  //                               init_system
  // *************************************************************************

  // set up the KKT System

  // KKTdiagx and KTTdiagy are both the positive compelementary systems,
  // but the one of KKTdiagy is on the dual side to KKTdiagx and taken with minus 


  int QPKKTSubspaceHPrecond::init_system(const Matrix& KKTdiagx,
    const Matrix& KKTdiagy,
    Real Hfac,
    Real /* prec */,
    QPSolverParameters* /* params */) {
    clock.start();
    t_gen_subspace = 0;
    t_comp_lowrank = 0;
    t_comp_svd = 0;
    t_precond_mult = 0;
    Q.init(0, 0, 0.);

    int status = 0;
    Hfactor = Hfac;
    assert(Hfactor > 0.);


    Diag_inv.init(diagH, Hfactor);
    Diag_inv += KKTdiagx;

    Real mindiag = max_Real;
    Real maxdiag = 0.;
    for (Integer i = 0; i < Diag_inv.dim(); i++) {
      Real d = Diag_inv(i);
      if (d < mindiag)
        mindiag = d;
      if (d > maxdiag)
        maxdiag = d;
    }
    assert(mindiag > 0.);
    assert(mindiag <= maxdiag);
    diaginvval = -1.;
    if ((maxdiag - mindiag) < 1e-10 * maxdiag)
      diaginvval = maxdiag;


    Integer nrows = 0;
    if (model)
      nrows += model->dim_model();
    if (Vp)
      nrows += Vp->coldim();
    if ((A) && (SchurComplAineq)) {
      nrows += A->rowdim() - ((eq_indices) ? eq_indices->rowdim() : 0);
    }

    Integer add_cols = 0;

    lowrank.init(Diag_inv.rowdim(), 0, 0.);


    //*************    select subspace   ****************
    if (nrows == 0) {
      subspace.init(nrows, 0, 0.);
    } else {
      switch (method / 10) {
      case 0:
      {
        //========  no preconditioning

        //method  0   ...   no preconditioning
        //method  1   ...   diagonal of H preconditioning

        subspace.init(nrows, 0, 0.);

        if (method == 1) {

          diaginvval = -1.;  //no longer valid, because Diag_inv will be modified

          if (Vp) {
            const Real* vp = Vp->get_store();
            for (Integer j = 0; j < Vp->coldim(); j++) {
              Real* dp = Diag_inv.get_store();
              const Real* const dpend = dp + Diag_inv.rowdim();
              for (; dp != dpend; vp++, dp++)
                (*dp) += Hfactor * (*vp) * (*vp);
            }
          } //endif (Vp)

          if (model) {
            model->add_BCSchur_diagonal(Diag_inv);
          }

          if ((A) && (SchurComplAineq)) {
            Integer eqnr = (eq_indices) ? eq_indices->rowdim() : 0;
            Integer eqi = 0;
            Integer eqind = (eqi < eqnr) ? (*eq_indices)(eqi++) : A->rowdim();
            for (Integer rinfoi = 0; rinfoi < A->get_rowinfo().rowdim(); rinfoi++) {
              Integer rind = A->get_rowinfo()(rinfoi, 0);
              while (eqind < rind)
                eqind = (eqi < eqnr) ? (*eq_indices)(eqi++) : A->rowdim();
              if (rind < eqind) {
                Integer rnz = A->get_rowinfo()(rinfoi, 1);
                Integer rbase = A->get_rowinfo()(rinfoi, 2);
                Real d = 1. / KKTdiagy(rind);
                for (Integer i = 0; i < rnz; i++) {
                  Diag_inv(A->get_rowindex()(rbase + i)) += d * sqr(A->get_rowval()(rbase + i));
                }
              }
            }
          } //endif (A)&&(SchurComplAineq)

        } //endif (method==1)

        break;
      } //end case 0

      case 1:
      {
        //========  exact  or  pure random Johnson-Lindenstrauss preconditioning

        //method  10   ...   exact preconditioning
        //method  11   ...   JL preconditioning, Achlioptas {-1,0,1} version
        //method  12   ...   JL preconditioning, normal distribution N(0,1/cols)

        //safeguarded Dasgupta/Gupta bound version of JL 
        Integer ncols = min(min(nrows, Diag_inv.rowdim()), min(Integer(200. * log(Real(Diag_inv.rowdim())) + .5), max(nrows / 3, 20)));
        //correct by experience from last iteration
        ncols = min(min(ncols, nrows), max(3, min(3 + 2 * eigvals.rowdim(), Integer(std::sqrt(last_nmult * (last_nmult / 4. + nrows / 4.)) - (last_nmult / 2.) + .5))));


        assert(ncols <= nrows);

        if ((ncols > .9 * nrows) || (method % 10 == 0)) {
          //==>  special case: full/exact precondtioning <==
          subspace.init_diag(nrows, 1.);
          if (cb_out(3)) {
            get_out() << " exact precond, ncols=" << nrows;
          }
        } else {
          //------ generate a random matrix
          switch (method % 10) {
          case 1:
          {
            // use Johnson Lindenstrauss scaling in the version of Achlioptas 2001
            Real scal = std::sqrt(3. / Real(max(ncols, 1)));
            subspace.init(nrows, ncols, 0.);
            for (Integer j = 0; j < ncols; j++) {
              for (Integer i = 0; i < nrows; i++) {
                Integer r = Integer(mat_randgen.unif_long(6));
                if (r == 0)
                  subspace(i, j) = -scal;
                else if (r == 1)
                  subspace(i, j) = scal;
              }
            }
            if (cb_out(3)) {
              get_out() << " JL-A precond, ev.rdim=" << eigvals.rowdim() << " ncols=" << ncols;
            }

            break;
          }

          case 2: default:
          {
            // use normal distribution with Johnson-Lindenstrauss dimension
            subspace.rand_normal(nrows, ncols, 0., 1. / Real(ncols));

            if (cb_out(3)) {
              get_out() << " JL precond, ev.rdim=" << eigvals.rowdim() << " ncols=" << ncols;
            }

            break;
          }
          } //end switch on subspace rand generator
        }

        break;
      }

      case 2:
      {
        //======== randomized projection low rank approximation

        // method 20 ... randomized subspace is updated from one system to the next
        // method 21 ... randomized subspace is regenerated for every system

        //safeguarded Dasgupta/Gupta bound version of JL 
        Integer ncols = min(min(nrows, Diag_inv.rowdim()), min(Integer(200. * log(Real(Diag_inv.rowdim())) + .5), max(nrows / 3, 20)));

        //correct by experience from last iteration
        ncols = min(min(ncols, nrows), max(3, min(3 + 2 * eigvals.rowdim(), Integer(std::sqrt(last_nmult * (last_nmult / 4. + nrows / 4.)) - (last_nmult / 2.) + .5))));

        assert(ncols < nrows);

        //------  generate random columns
        if ((subspace.coldim() == 0) || (method % 10 == 1)) {
          if (ncols > .9 * nrows) {
            subspace.init_diag(nrows, 1.);

            if (cb_out(3)) {
              get_out() << " exact precond, ncols=" << nrows;
            }
          } else {
            tmpmat.rand_normal(nrows, ncols, 0., 1.);
            Indexmatrix piv;
            ncols = tmpmat.QR_factor(piv);
          }
          if (cb_out(3)) {
            get_out() << " rand subsp-precond, ";
          }
        } else {
          tmpmat.init(subspace);
          Integer r = tmpmat.QR_factor(pivlowrank);
          add_cols = max(3, min(ncols, Integer(sqrt(last_nmult) / 2.)) - subspace.coldim());
          subspace.rand_normal(tmpmat.rowdim(), add_cols, 0., 1.);
          ncols = tmpmat.QR_concat_right(subspace, pivlowrank, r);
          if (cb_out(3)) {
            get_out() << " rand subspinc-precond lastmult=" << last_nmult << " addcols=" << add_cols;
          }
        }
        subspace.init(nrows, ncols, 0.);
        for (Integer i = 0; i < ncols; i++)
          subspace(i, i) = 1.;
        tmpmat.Q_times(subspace, ncols);

        if (cb_out(3)) {
          get_out() << " newdim=" << ncols;
        }

        break;
      }

      case 3:
      {
        //======== proposed/constructed projection for low rank approximation
        Real minval = pow(10., (method % 10) * .5);
        if (cb_out(3)) {
          get_out() << " sing_est[" << method << "](" << minval << ";";
        }
        Matrix tmpvec;
        Integer nmod = 0;
        Integer nV = 0;
        Integer nA = 0;
        Diag_inv.inv();
        diaginvval = 1. / diaginvval;
        tmpvec.init(0, 1, 0.);
        lowrank.init(Diag_inv.rowdim(), 0, 0.);
        if (model) {
          model->propose_BCSchur_pcsubspace(lowrank, tmpvec, Diag_inv, minval, diaginvval);
          nmod = lowrank.coldim();
        }
        if (Vp) {
          for (Integer i = 0; i < Vp->coldim(); i++) {
            Real d;
            if (diaginvval > 0.) {
              d = std::sqrt(diaginvval * colip(*Vp, i));
            } else {
              d = std::sqrt(colip(*Vp, i, &Diag_inv));
            }
            if (d > minval) {
              tmpvec.concat_below(d);
              lowrank.concat_right(Vp->col(i));
              nV++;
            }
          }
        }
        if (((A) && (SchurComplAineq)) &&
          ((eq_indices == 0) || (eq_indices->rowdim() < A->rowdim()))) {
          Integer eqi = 0;
          Integer eqind = (((eq_indices == 0) || (eqi >= eq_indices->rowdim())) ? -1 : (*eq_indices)(eqi));
          for (Integer i = 0; i < A->rowdim(); i++) {
            if (i == eqind) {
              eqi++;
              if (eqi < eq_indices->rowdim())
                eqind = (*eq_indices)(eqi);
              else
                eqind = -1;
            } else {
              Real d;
              if (diaginvval > 0.) {
                d = std::sqrt(diaginvval * rowip(*A, i) / KKTdiagy(i));
              } else {
                d = std::sqrt(rowip(*A, i, &Diag_inv) / KKTdiagy(i));
              }
              if (d > minval) {
                std::cout << " A" << i << "[" << d << "]";
                tmpvec.concat_below(d);
                tmpmat.init(A->row(i), std::sqrt(1. / KKTdiagy(i)));
                tmpmat.transpose();
                lowrank.concat_right(tmpmat);
                nA++;
              }
            }
          }
        }
        if (cb_out(3)) {
          get_out() << minval << ";" << nmod << "+" << nV << "+" << nA << ")";
        }
        //assert(norm2(transpose(subspace)*subspace-Diag(Matrix(subspace.coldim(),1,1.)))<1e-10*subspace.coldim());
        // // TEST begin
        // Real d=norm2(transpose(subspace)*subspace-Diag(Matrix(subspace.coldim(),1,1.)));
        // if (d>1e-10*subspace.coldim()){
        // 	std::cout<<" nrmdev="<<d<<" dev="<<transpose(subspace)*subspace-Diag(Matrix(subspace.coldim(),1,1.));
        // }
        // // TEST end
        Diag_inv.inv();
        diaginvval = 1. / diaginvval;
        break;
      }

      case 4:
      {
        //======== externally given subspace, assumed to be orthonormalized
        assert((subspace.coldim() == 0) || (subspace.rowdim() == Diag_inv.rowdim()));
        //normalize
        tmpmat.init(subspace);
        Integer r = tmpmat.QR_factor(pivlowrank);
        subspace.init(tmpmat.rowdim(), r, 0.);
        assert(r <= tmpmat.rowdim());
        for (Integer i = 0; i < r; i++) {
          subspace(i, i) = 1.;
        }
        tmpmat.Q_times(subspace, r);
        break;
      }

      }  //end switch subspace selection
    } //endif (nrows>0)  subspace selection

    t_gen_subspace = clock.time();

    //*************    multiply with subspace if needed   ****************

    Diag_inv.inv();
    diaginvval = 1. / diaginvval;

    if ((subspace.coldim() > 0) && (lowrank.coldim() == 0)) {

      //------ multiply with the low rank part
      if (cb_out(3))
        get_out() << " project " << subspace.coldim() << "(" << subspace.rowdim() << ")";

      lowrank.init(Diag_inv.rowdim(), subspace.coldim(), 0.);
      Integer start_row = 0;
      if (model) {
        tmpmat = subspace.rows(Range(start_row, start_row + model->dim_model() - 1));
        model->prepare_BCSchur_JLprecond(lowrank, tmpmat);
        start_row += model->dim_model();
      }
      if (Vp) {
        tmpmat = subspace.rows(Range(start_row, start_row + Vp->coldim() - 1));
        genmult(*Vp, tmpmat, lowrank, std::sqrt(Hfactor), 1.);
        start_row += Vp->coldim();
      }
      if ((A) && (SchurComplAineq) &&
        ((eq_indices == 0) || (eq_indices->rowdim() < A->rowdim()))) {
        tmpmat.init(A->rowdim(), subspace.coldim(), 0.);
        Integer row = 0;
        Integer eqnr = (eq_indices) ? eq_indices->rowdim() : 0;
        for (Integer i = 0; i < eqnr; i++) {
          Integer ii = (*eq_indices)(i);
          for (; row < ii; row++) {
            mat_xeya(tmpmat.coldim(), tmpmat.get_store() + row, tmpmat.rowdim(), subspace.get_store() + start_row + row - i, subspace.rowdim(), 1. / std::sqrt(KKTdiagy(row)));
          }
          row++;
        }
        for (; row < tmpmat.rowdim(); row++) {
          mat_xeya(tmpmat.coldim(), tmpmat.get_store() + row, tmpmat.rowdim(), subspace.get_store() + start_row + row - eqnr, subspace.rowdim(), 1. / std::sqrt(KKTdiagy(row)));
        }
        genmult(*A, tmpmat, lowrank, 1., 1., 1);
        start_row += A->rowdim() - ((eq_indices == 0) ? 0 : eq_indices->rowdim());
      }

    }

    t_comp_lowrank = clock.time() - t_gen_subspace;

    //*************   compute SVD (and update subspace)   ****************

    eigvals.init(0, 0, 0.);
    eigvecs.init(0, 0, 0.);

    if ((subspace.coldim() > 0) || (lowrank.coldim() > 0)) {

      //--- determine singular value decomposition
      if (lowrank.coldim() < 5) {
        if (diaginvval > 0.)
          rankadd(lowrank, tmpsym, diaginvval, 0., 1);
        else
          scaledrankadd(lowrank, Diag_inv, tmpsym, 1., 0., 1);
      } else {
        tmpmat.init(lowrank, 1., 1);
        if (diaginvval > 0.)
          rankadd(tmpmat, tmpsym, diaginvval, 0.);
        else
          scaledrankadd(tmpmat, Diag_inv, tmpsym, 1., 0.);
      }
      if (cb_out(3)) {
        get_out() << " diag=" << transpose(diag(tmpsym));
      }

      Real maxdiag = 0.;
      Real mindiag = max_Real;
      for (Integer i = 0; i < tmpsym.rowdim(); i++) {
        Real d = tmpsym(i, i);
        if (d > maxdiag)
          maxdiag = d;
        if (d < mindiag)
          mindiag = d;
      }
      if ((mindiag < 0.) && (cb_out())) {
        get_out() << "**** ERROR: QPKKTSubspaceHPrecond::init_system(): in preconditioning method " << method << " the Gram matrix has a negative diagonal element, might be numerica ldifficulties" << std::endl;
      }
      tmpsym *= 1. / maxdiag;
      int retcode;
      if ((retcode = tmpsym.eig(eigvecs, eigvals, false))) {
        if (cb_out())
          get_out() << "**** ERROR: QPKKTSubspaceHPrecond::init_system(): tmpsym.eig failed for preconditioning method " << method << " and returned " << retcode << " for matrix " << tmpsym << ", which is A*A' for A=" << lowrank << std::endl;
        status++;
      }
      eigvals *= maxdiag;
      max_sigma = eigvals(0);

      // //BEGIN TEST
      // if (eigvals.rowdim()>0){
      //   //--- remove small eigenvalues,  add one and invert
      //   Integer keep=0;
      //   Matrix leigvals=eigvals;
      //   Matrix leigvecs=eigvecs;
      //   while((keep<leigvals.rowdim())&&(leigvals(keep)>1.)){
      // 	keep++;
      //   }
      //   if (keep<leigvals.rowdim()){
      // 	leigvals.reduce_length(keep);
      // 	leigvecs.delete_cols(Range(keep,leigvecs.coldim()-1));
      //   }

      //   genmult(lowrank,leigvecs,Q);
      //   tmpvec=Diag_inv;
      //   tmpvec.sqrt();
      //   Q.scale_rows(tmpvec);
      //   rotmat=leigvals;
      //   rotmat.inv();
      //   rotmat.sqrt();
      //   Q.scale_cols(rotmat);
      //   if(norm2(transpose(Q)*Q-Diag(Matrix(Q.coldim(),1,1.)))>1e-8*max(leigvals)){
      // 	std::cout<<" eig["<<max(leigvals)<<","<< min(leigvals)<<"] bad Q'*Q="<<transpose(Q)*Q;
      //   }	
      //   int lerr=Q.QR_factor();
      //   if (lerr){
      // 	if (cb_out())
      // 	  get_out()<<"**** ERROR in QPKKTSubspaceHPrecond test: QR_factor failed and returned "<<lerr<<std::endl;
      //   }
      //   else {
      // 	if ((method/10<=9)&&(subspace.coldim()>0)){
      // 	  genmult(subspace,leigvecs,rotmat);
      // 	}
      // 	else {
      // 	  rotmat.init(0,0,0.);
      // 	}
      // 	Real maxlam=0;
      // 	Matrix tmpvec2;
      // 	std::cout<<" testmult";
      // 	for(Integer i=0;i<leigvals.rowdim();i++){
      // 	  Real lam=std::sqrt(leigvals(i));
      // 	  std::cout<<" ("<<i<<","<<lam;
      // 	  tmpvec2.init(Q.rowdim(),1,0.);
      // 	  keepvecs.init(Q.rowdim(),1,0.);
      // 	  keepvecs(i)=1.;
      // 	  Q.Q_times(keepvecs,Q.coldim());
      // 	  keepvecs%=tmpvec;
      // 	  if (model)
      // 	    model->add_Schur_mult(keepvecs,tmpvec2);
      // 	  tmpvec2%=tmpvec;
      // 	  keepvecs/=tmpvec;
      // 	  std::cout<<","<<std::sqrt(ip(tmpvec2,keepvecs));
      // 	  if (rotmat.coldim()==leigvals.rowdim()){
      // 	    keepvecs=rotmat.col(i);
      // 	    tmpvec2.init(Q.rowdim(),1,0.);
      // 	    if (model)
      // 	      model->prepare_BCSchur_JLprecond(tmpvec2,keepvecs);
      // 	    tmpvec2%=tmpvec;
      // 	    std::cout<<","<<std::sqrt(ip(tmpvec2,tmpvec2));
      // 	  }
      // 	  std::cout<<")"<<std::flush;
      // 	  if (lam>maxlam)
      // 	    maxlam=lam;
      // 	}
      //   }
      //   Q.init(0,0,0.);
      // }
      // //END TEST


      if (eigvals.rowdim() > 0) {
        //--- if subspace is used, update subspace   (??? check for enlargements ???)
        if (method == 20) {
          if (subspace.coldim() > subspace.rowdim() * .9) {
            subspace.init_diag(subspace.rowdim(), 1.);
            if (cb_out(3)) {
              get_out() << " subspupd rdim=" << subspace.rowdim() << " cdim=" << subspace.coldim() << " identity";
            }
          } else {
            Real eigthresh = max(10., exp(.1 * log(max(eigvals)) + .9 * log(max(1e-12, min(eigvals)))));
            Integer keep = min(3, eigvals.rowdim());
            //Real eigthresh=exp(log(max(eigvals))/3.+log(max(1e-12,min(eigvals)))*2./3.);
            if (method == 60) {
              eigthresh = 10;
              keep = 0;
            }
            if (method == 70) {
              eigthresh = 100;
              keep = 0;
            }
            while ((keep < eigvals.rowdim()) && (eigvals(keep) > eigthresh)) {
              keep++;
            }
            if (last_nmult < subspace.coldim())
              keep = min(keep, subspace.coldim() - add_cols);
            keepvecs.init(eigvecs.rowdim(), keep, eigvecs.get_store());
            keepeigs.init(keep, 1, eigvals.get_store(), 1.);
            //std::cout<<" eigthresh="<<maxdiag*eigthresh<<" keepeigs="<<transpose(keepeigs);
            genmult(subspace, keepvecs, tmpmat);
            subspace = tmpmat;
            if (cb_out(3)) {
              get_out() << " subspupd keep=" << keep;
            }
          }
        }
      }

      //--- remove small eigenvalues,  add one and invert
      Integer keep = 0;
      while ((keep < eigvals.rowdim()) && (eigvals(keep) > 1.)) {
        keep++;
      }
      if (keep < eigvals.rowdim()) {
        eigvals.reduce_length(keep);
        eigvecs.delete_cols(Range(keep, eigvecs.coldim() - 1));
      }

      if (keep > 0) {
        eigvals += 1.;
        eigvals.inv();
        if (cb_out(3)) {
          get_out() << " pcglrnc=" << lowrank.coldim() << " rank=" << eigvals.rowdim() << " lmax=" << eigvals(0) << " lmin=" << eigvals(eigvals.rowdim() - 1) << " eigs=" << transpose(eigvals) << std::endl;
        }

        //--- check if multiplying lowrank with eigvecs is better to do now
        if (last_nmult * (lowrank.coldim() - eigvals.rowdim()) > eigvals.rowdim()) {
          genmult(lowrank, eigvecs, tmpmat);
          lowrank = tmpmat;
          eigvecs.init(0, 0, 0.);
        }

      } //endif keep>0

    }

    last_nmult = 0;

    t_comp_svd = clock.time() - (t_comp_lowrank + t_gen_subspace);

    if (cb_out(3)) {
      get_out() << " prept " << clock.time();
      get_out() << "[" << t_gen_subspace;
      get_out() << "," << t_comp_lowrank;
      get_out() << "," << t_comp_svd;
      get_out() << "]";
    }

    return status;
  }

  // *************************************************************************
  //                             precondM1
  // *************************************************************************

    ///return (an estimate of) the minimum eigenvalue of the preconditioner M1^{-1}; this is used, e.g., to correct the precission in MINRES

  Real QPKKTSubspaceHPrecond::get_lmin_invM1() {
    if (method == 0)
      return 1.;
    if (eigvals.rowdim() > 0) {
      Real sum = 0;
      for (Integer i = 0; i < eigvals.rowdim(); i++) {
        sum += std::log(eigvals(i));
      }
      sum /= Diag_inv.rowdim();
      return min(Diag_inv) * std::exp(sum);
    }
    return min(Diag_inv);
  }

  // *************************************************************************
  //                             precondM1
  // *************************************************************************


  int QPKKTSubspaceHPrecond::precondM1(Matrix& vec) {

    last_nmult++;

    if (method != 0) {
      clock.start();

      if (diaginvval > 0.) {
        tmpmat.init(Diag_inv.rowdim(), 1, vec.get_store(), 1, diaginvval);
      } else {
        tmpmat.init(Diag_inv.rowdim(), 1, vec.get_store());
        tmpmat %= Diag_inv;
      }

      if (eigvals.rowdim() > 0) {
        genmult(lowrank, tmpmat, rotmat, 1., 0., 1);
        if (eigvecs.coldim() == eigvals.rowdim()) {
          genmult(eigvecs, rotmat, tmpvec, 1., 0., 1);
          tmpvec %= eigvals;
          genmult(eigvecs, tmpvec, rotmat);
        } else
          rotmat %= eigvals;
        genmult(lowrank, rotmat, tmpvec);
        if (diaginvval > 0.) {
          tmpmat.xpeya(tmpvec, -diaginvval);
        } else {
          tmpvec %= Diag_inv;
          tmpmat -= tmpvec;
        }
      }

      mat_xey(Diag_inv.rowdim(), vec.get_store(), tmpmat.get_store());

      // if (A){
      // 	//currently use the identity, but maybe something else would be better
      // }

      t_precond_mult += clock.time();
    }


    return 0;
  }

  // *************************************************************************
  //                            cond_number_mult
  // *************************************************************************

  ///for estimating the condition number with M1=G*G^T this returns G^{-1}*vec; default: G=I
  int QPKKTSubspaceHPrecond::cond_number_mult(Matrix& vec,
    const Matrix& KKTdiagx,
    const Matrix& KKTdiagy) {
    assert(vec.rowdim() == Diag_inv.rowdim());

    if (method != 0) {
      assert(vec.dim() == Diag_inv.dim());
      assert(min(Diag_inv) > eps_Real);

      tmpvec = Diag_inv;
      tmpvec.sqrt();

      if (eigvals.rowdim() > 0) {
        //make sure Q is available
        if (Q.coldim() != eigvals.rowdim()) {
          //if not, form D{-.5}*lowrank*eigvecs and compute Householder QR
          if (eigvecs.coldim() == eigvals.rowdim()) {
            genmult(lowrank, eigvecs, Q);
          } else {
            //lowrank is already mutlipied by eigvecs
            Q = lowrank;
          }
          Q.scale_rows(tmpvec);

          int err = Q.QR_factor();
          if (err) {
            if (cb_out())
              get_out() << "**** ERROR in QPKKTSubspaceHPrecond::cond_number_mult(): the QR_factor failed and returned " << err << std::endl;
            return err;
          }

          // //BEGIN TEST
          // Q.Qt_times(tmpmat,Q.coldim());
          // Real maxlam=0;
          // for(Integer i=0;i<eigvals.rowdim();i++){
          //   Real lam=std::sqrt(1./eigvals(i)-1.);
          //   if (lam>maxlam)
          //     maxlam=lam;
          //   lam*=(tmpmat(i,i)<0?-1.:1.);
          //   if(std::fabs(tmpmat(i,i)-lam)>1e-8*maxlam){
          //     std::cout.precision(10);
          //     std::cout<<" testfail: i="<<i<<" lam="<<lam<<" tmpmat="<<tmpmat<<std::endl;
          //     std::abort();
          //   }
          //   tmpmat(i,i)-=lam;
          // }
          // assert(norm2(tmpmat)<1e-8*maxlam);
          // //END TEST

        }

        //multiply with G^{-T}
        for (Integer i = 0; i < eigvals.rowdim(); i++) {
          Real d = std::sqrt(eigvals(i));
          vec(i) *= d;
        }
        Q.Q_times(vec, Q.coldim());

      }

      //Matrix tvec=vec; //TEST

      // multiply with (I+D^(-.5)*fatV*fatV'*D^(-.5))
      tmpmat = vec;
      tmpmat %= tmpvec;
      if (Vp) {
        genmult(*Vp, tmpmat, keepvecs, 1., 0., 1);
        genmult(*Vp, keepvecs, rotmat);
      } else {
        rotmat.init(Diag_inv.rowdim(), 1, 0.);
      }
      if (model) {
        model->add_Schur_mult(tmpmat, rotmat);
      }
      Integer eqnr = (eq_indices) ? eq_indices->rowdim() : 0;
      if ((A) && (KKTdiagy.rowdim() > eqnr)) {
        genmult(*A, tmpmat, keepvecs);
        Integer j = 0;
        for (Integer i = 0; i < eqnr; i++) {
          Integer ii = (*eq_indices)(i);
          for (; j < ii; j++) {
            keepvecs(j) /= KKTdiagy(j);
          }
          keepvecs(j++) = 0.;
        }
        for (; j < keepvecs.rowdim(); j++) {
          keepvecs(j) /= KKTdiagy(j);
        }
        genmult(*A, keepvecs, rotmat, 1., 1., 1);
      }
      rotmat %= tmpvec;
      vec += rotmat;

      // //BEGIN TEST
      // {
      //   Real HRitz=ip(tvec,vec)/ip(tvec,tvec);
      //   Matrix tv2=tvec;
      //   if (eigvals.rowdim()>0){
      // 	Q.Qt_times(tv2,Q.coldim());
      // 	for(Integer i=0;i<eigvals.rowdim();i++){
      // 	  tv2(i)*=1/eigvals(i);
      // 	}
      // 	Q.Q_times(tv2,Q.coldim());
      //   }
      //   Real CondRitz=ip(tvec,tv2)/ip(tvec,tvec);
      //   if (HRitz*(1.+1e-6)<CondRitz){
      // 	std::cout.precision(10);
      // 	std::cout<<" HRitz="<<HRitz<<" CondRitz="<<CondRitz; 
      // 	tv2.init(tvec.rowdim(),1,0.);
      // 	if (model)
      // 	  model->add_Schur_mult(tvec,tv2);
      // 	std::cout<<" Hmod="<<ip(tv2,tvec)/ip(tvec,tvec);
      // 	genmult(lowrank,tvec,tv2,1.,0.,1);
      // 	std::cout<<" CondV="<<ip(tv2,tv2)/ip(tvec,tvec);
      // 	tv2%=tv2;
      // 	tv2/=ip(tvec,tvec);
      // 	std::cout<<" elems="<<transpose(tv2)<<std::endl;
      //   }
      // }
      // //END TEST

      if (eigvals.rowdim() > 0) {
        //multiply with G^{-1}
        Q.Qt_times(vec, Q.coldim());

        for (Integer i = 0; i < eigvals.rowdim(); i++) {
          Real d = std::sqrt(eigvals(i));
          vec(i) *= d;
        }
      }
    } else {
      //multiply tmpmat by the H-block of the system into vec
      tmpmat.init(Diag_inv.rowdim(), 1, vec.get_store());
      vec %= KKTdiagx;
      Hp->add_Hx(tmpmat, vec, Hfactor);
      if (model) {
        model->add_Schur_mult(tmpmat, vec);
      }
      Integer eqnr = (eq_indices) ? eq_indices->rowdim() : 0;
      if ((A) && (KKTdiagy.rowdim() > eqnr)) {
        keepvecs.newsize(KKTdiagy.rowdim(), 1); chk_set_init(keepvecs, 1);
        genmult(*A, tmpmat, keepvecs, 1., 0.);
        Integer j = 0;
        for (Integer i = 0; i < eqnr; i++) {
          Integer ii = (*eq_indices)(i);
          for (; j < ii; j++) {
            keepvecs(j) /= KKTdiagy(j);
          }
          keepvecs(j++) = 0.;
        }
        for (; j < keepvecs.rowdim(); j++) {
          keepvecs(j) /= KKTdiagy(j);
        }
        genmult(*A, keepvecs, vec, 1., 1., 1);
      }
    }

    return 0;
  }


  // *************************************************************************
  //                             precond_invG1
  // *************************************************************************

  ///for estimating the condition number with M1=G*G^T this returns G^{-1}*vec; default: G=I
  int QPKKTSubspaceHPrecond::precond_invG1(Matrix& vec) {
    if (method != 0) {
      assert(vec.rowdim() >= Diag_inv.rowdim());

      tmpmat.init(Diag_inv.rowdim(), 1, vec.get_store());
      tmpvec = Diag_inv;
      tmpvec.sqrt();
      tmpmat %= tmpvec;

      if (eigvals.rowdim() > 0) {
        //make sure Q is available
        if (Q.coldim() != eigvals.rowdim()) {
          //if not, form D{-.5}*lowrank*eigvecs and compute Householder QR
          if (eigvecs.coldim() == eigvals.rowdim()) {
            genmult(lowrank, eigvecs, Q);
          } else {
            //lowrank is already mutlipied by eigvecs
            Q = lowrank;
          }

          int err = Q.QR_factor();
          if (err) {
            if (cb_out())
              get_out() << "**** ERROR in QPKKTSubspaceHPrecond::precond_invG1(): QR_factor failed and returned " << err << std::endl;
            return err;
          }
        }

        //multiply with G^{-1}
        Q.Qt_times(tmpmat, Q.coldim());

        for (Integer i = 0; i < eigvals.rowdim(); i++) {
          Real d = std::sqrt(eigvals(i));
          tmpmat(i) *= d;
        }
      }

      mat_xey(Diag_inv.rowdim(), vec.get_store(), tmpmat.get_store());
    }
    return 0;
  }


  // *************************************************************************
  //                             precond_invG1trans
  // *************************************************************************

  ///for estimating the condition number with M1=G*G^T this returns G^{-T}*vec; default: G=I
  int QPKKTSubspaceHPrecond::precond_invG1tran(Matrix& vec) {
    if (method != 0) {
      assert(vec.rowdim() >= Diag_inv.rowdim());
      assert(min(Diag_inv) > eps_Real);
      tmpmat.init(Diag_inv.rowdim(), 1, vec.get_store());

      if (eigvals.rowdim() > 0) {
        //make sure Q is available
        if (Q.coldim() != eigvals.rowdim()) {
          //if not, form D{-.5}*lowrank*eigvecs and compute Householder QR
          if (eigvecs.coldim() == eigvals.rowdim()) {
            genmult(lowrank, eigvecs, Q);
          } else {
            //lowrank is already mutlipied by eigvecs
            Q = lowrank;
          }
          int err = Q.QR_factor();
          if (err) {
            if (cb_out())
              get_out() << "**** ERROR in QPKKTSubspaceHPrecond::precond_invG1tran(): QR_factor failed and returned " << err << std::endl;
            return err;
          }
        }

        //multiply with G^{-T}
        for (Integer i = 0; i < eigvals.rowdim(); i++) {
          Real d = std::sqrt(eigvals(i));
          vec(i) *= d;
        }
        Q.Q_times(vec, Q.coldim());

      }

      tmpvec = Diag_inv;
      tmpvec.sqrt();
      tmpmat %= tmpvec;

      mat_xey(Diag_inv.rowdim(), vec.get_store(), tmpmat.get_store());


    }
    return 0;
  }



}

