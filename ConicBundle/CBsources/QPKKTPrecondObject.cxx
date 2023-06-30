/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/QPKKTPrecondObject.cxx
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
#include "lanczpol.hxx"
#include "QPKKTPrecondObject.hxx"

using namespace CH_Matrix_Classes;

namespace ConicBundle {

  QPKKTPrecondObject::~QPKKTPrecondObject() {
  }

  // *************************************************************************
  //                             QPget_condition_number
  // *************************************************************************

  Real QPKKTPrecondObject::get_condition_number(const Matrix& KKTdiagx,
    const Matrix& KKTdiagy) {
    if (precond_size() < 0)
      return -1.;

    class KKTmat : public CH_Matrix_Classes::Lanczosmatrix {
    private:
      QPKKTPrecondObject* itsol;
      const Matrix& KKTdiagx;
      const Matrix& KKTdiagy;
      bool positive;
    public:
      KKTmat(QPKKTPrecondObject* initsol,
        const Matrix& inKKTdiagx,
        const Matrix& inKKTdiagy) :
        itsol(initsol), KKTdiagx(inKKTdiagx), KKTdiagy(inKKTdiagy), positive(true) {
        assert(initsol);
      }

      void set_positive(bool pos) {
        positive = pos;
      }

      ///returns the order of the (virtual) symmetric matrix
      Integer lanczosdim() const {
        return itsol->precond_size();
      }

      ///returns a rough estimate on the number of flops needed by lanczosmult() for a  vector
      Integer lanczosflops() const {
        return lanczosdim() * 2 * (1 + (itsol->A ? (itsol->A->rowdim() - (itsol->eq_indices ? itsol->eq_indices->rowdim() : 0)) : 0) + (itsol->model ? itsol->model->dim_model() : 0));
      }

      ///computes  B = (*this) * A; A and B must not be the same object!
      int lanczosmult(const Matrix& A, Matrix& B) const {
        Matrix bb(A.rowdim(), 1); chk_set_init(bb, 1);
        B.newsize(A.rowdim(), A.coldim()); chk_set_init(B, 1);
        int err = 0;
        for (Integer j = 0; j < A.coldim(); j++) {
          mat_xeya(A.rowdim(), bb.get_store(), A.get_store() + j * A.rowdim(), positive ? 1 : -1.);
          err = itsol->cond_number_mult(bb, KKTdiagx, KKTdiagy);
          if (err)
            break;
          mat_xey(A.rowdim(), B.get_store() + j * A.rowdim(), bb.get_store());
        }
        return err;
      }
    };

    KKTmat kkt(this, KKTdiagx, KKTdiagy);
    Matrix eigval;
    Matrix eigvec;
    Lanczpol lancz;
    //lancz.set_out(get_out_ptr(),get_print_level());
    //lancz.set_out(&std::cout,10);
    lancz.set_nchebit(10);
    lancz.set_maxiter(10);
    lancz.set_maxmult(500);
    lancz.set_relprec(1.e-4);
    int err = lancz.compute(&kkt, eigval, eigvec, 1);
    if ((err) && (cb_out(10))) {
      get_out() << "**** ERROR in QPKKTPrecondObject::get_condition_number(): lancz.compute returned " << err << std::endl;
    }
    Real maxeig = max(eigval);
    if (eigval.dim() == 0) {
      lancz.get_lanczosvecs(eigval, eigvec);
      maxeig = max(eigval);
      if (cb_out(10)) {
        get_out() << " **** WARNING: eigval.dim==0";
        get_out() << " Ritzval.dim=" << eigval.dim() << " maxeig=max(Ritzval)=" << maxeig << std::endl;
      }
    }
    kkt.set_positive(false);
    eigval.init(0, 0, 0.);
    eigvec.init(0, 0, 0.);
    err = lancz.compute(&kkt, eigval, eigvec, 1);
    if ((err) && (cb_out(10))) {
      get_out() << "**** ERROR in QPKKTPrecondObject::get_condition_number(): lancz.compute returned " << err << std::endl;
    }
    Real mineig = -max(eigval);
    if (eigval.dim() == 0) {
      lancz.get_lanczosvecs(eigval, eigvec);
      mineig = -max(eigval);
      if (cb_out(10)) {
        get_out() << " **** WARNING: eigval.dim==0";
        get_out() << " Ritzval.dim=" << eigval.dim() << " mineig=-max(Ritzval)=" << mineig << std::endl;
      }
    }
    Real cond = maxeig / mineig;
    if (cb_out(2)) {
      get_out() << " lanczcond=" << maxeig << "/" << mineig << "=" << cond << std::endl;
    }

    // //BEGINT test
    // if(kkt.lanczosdim()<500){
    //   kkt.set_positive(true);
    //   Matrix A(kkt.lanczosdim(),kkt.lanczosdim(),0.);
    //   for (Integer i=0;i<A.rowdim();i++)
    //     A(i,i)=1.;
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

