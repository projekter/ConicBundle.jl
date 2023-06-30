/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/UQPSumModelBlock.cxx
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



#include "UQPSumModelBlock.hxx"


using namespace CH_Matrix_Classes;

namespace ConicBundle {

  UQPSumModelBlock::~UQPSumModelBlock() {
  }

  int UQPSumModelBlock::append(QPModelDataObject* inblock) {
    if (inblock == 0)
      return 0;
    UQPModelBlock* qpbp = dynamic_cast<UQPModelBlock*>(inblock);
    if (qpbp == 0) {
      if (cb_out())
        get_out() << "**** ERROR in QP_SumModelBlock::append(): not a derived UQPModelBlock*" << std::endl;
      return 1;
    }
    blocks.push_back(qpbp);


    if (bundle.size() == 0) {
      //initialize
      bundle.push_back(inblock->get_bundle());
      constant_minorant.push_back(inblock->get_constant_minorant());
    } else {
      if (!inblock->get_constant_minorant().empty())
        inblock->get_constant_minorant().get_minorant(get_constant_minorant());
      get_bundle().insert(get_bundle().end(), inblock->get_bundle().begin(), inblock->get_bundle().end());
    }
    return 0;
  }

  Integer UQPSumModelBlock::xdim() const {
    Integer dim = 0;
    for (unsigned int i = 0; i < blocks.size(); i++) dim += blocks[i]->xdim();
    return dim;
  }

  Integer UQPSumModelBlock::ydim() const {
    Integer dim = 0;
    for (unsigned int i = 0; i < blocks.size(); i++) dim += blocks[i]->ydim();
    return dim;
  }

  int UQPSumModelBlock::set_qp_xstart(Integer x_start_index) {
    int retval = 0;
    xstart = x_start_index;
    xend = xstart;
    for (unsigned int i = 0; i < blocks.size(); i++) {
      retval |= blocks[i]->set_qp_xstart(xend);
      xend += blocks[i]->xdim();
    }
    return retval;
  }

  int UQPSumModelBlock::set_qp_ystart(Integer y_start_index) {
    int retval = 0;
    ystart = y_start_index;
    yend = ystart;
    for (unsigned int i = 0; i < blocks.size(); i++) {
      retval |= blocks[i]->set_qp_ystart(yend);
      yend += blocks[i]->ydim();
    }
    return retval;
  }

  int UQPSumModelBlock::starting_x(Matrix& qp_x) {
    int retval = 0;
    for (unsigned int i = 0; i < blocks.size(); i++) {
      retval |= blocks[i]->starting_x(qp_x);
    }
    return retval;
  }

  int UQPSumModelBlock::starting_y(Matrix& qp_y,
    const Matrix& qp_Qx,
    const Matrix& qp_c) {
    int retval = 0;
    for (unsigned int i = 0; i < blocks.size(); i++) {
      retval |= blocks[i]->starting_y(qp_y, qp_Qx, qp_c);
    }
    return retval;
  }

  Real UQPSumModelBlock::get_local_primalcost() const {
    Real sum = 0.;
    for (unsigned int i = 0; i < blocks.size(); i++) {
      sum += blocks[i]->get_local_primalcost();
    }
    return sum;
  }

  Real UQPSumModelBlock::get_local_dualcost() const {
    Real sum = 0.;
    for (unsigned int i = 0; i < blocks.size(); i++) {
      sum += blocks[i]->get_local_dualcost();
    }
    return sum;
  }


  int UQPSumModelBlock::get_Ab(Matrix& qp_A, Matrix& qp_b) const {
    int retval = 0;
    for (Integer j = xstart; j < xend; j++) { //set entire block to zero
      mat_xea(yend - ystart, qp_A.get_store() + qp_A.rowdim() * j + ystart, 0.);
    }
    for (unsigned int i = 0; i < blocks.size(); i++) {
      retval |= blocks[i]->get_Ab(qp_A, qp_b);
    }
    return retval;
  }



  int UQPSumModelBlock::restart_x(Matrix& qp_x, const Matrix& qp_c, const Matrix& qp_dc) {
    int retval = 0;
    for (unsigned int i = 0; i < blocks.size(); i++) {
      retval |= blocks[i]->restart_x(qp_x, qp_c, qp_dc);
    }
    return retval;
  }


  int UQPSumModelBlock::restart_y(Matrix& qp_y,
    const Matrix& qp_Qx,
    const Matrix& qp_c,
    const Matrix& qp_dc) {
    int retval = 0;
    for (unsigned int i = 0; i < blocks.size(); i++) {
      retval |= blocks[i]->restart_y(qp_y, qp_Qx, qp_c, qp_dc);
    }
    return retval;
  }

  int UQPSumModelBlock::add_xinv_kron_z(Symmatrix& barQ) {
    int retval = 0;
    for (unsigned int i = 0; i < blocks.size(); i++) {
      retval |= blocks[i]->add_xinv_kron_z(barQ);
    }
    return retval;
  }


  int UQPSumModelBlock::add_local_sys(Symmatrix& sysdy, Matrix& rhs) {
    int retval = 0;
    for (unsigned int i = 0; i < blocks.size(); i++) {
      retval |= blocks[i]->add_local_sys(sysdy, rhs);
    }
    return retval;
  }


  int UQPSumModelBlock::suggest_mu(Real& ip_xz,
    Integer& mu_dim,
    Real& sigma,
    const Matrix& qp_dx,
    const Matrix& qp_dy,
    const Matrix& rhs_residual) {
    int retval = 0;
    ip_xz = 0.;
    mu_dim = 0;
    sigma = 0.;
    Real ipxz, s;
    Integer md;
    for (unsigned int i = 0; i < blocks.size(); i++) {
      retval |= blocks[i]->suggest_mu(ipxz, md, s, qp_dx, qp_dy, rhs_residual);
      ip_xz += ipxz;
      mu_dim += md;
      sigma = max(s, sigma);
    }
    return retval;
  }


  int UQPSumModelBlock::get_corr(Matrix& xcorr,
    Matrix& rhs,
    Real mu) {
    int retval = 0;
    for (unsigned int i = 0; i < blocks.size(); i++) {
      retval |= blocks[i]->get_corr(xcorr, rhs, mu);
    }
    return retval;
  }


  int UQPSumModelBlock::line_search(Real& alpha,
    const Matrix& qp_dx,
    const Matrix& qp_dy,
    const Matrix& rhs_residual) {
    int retval = 0;
    for (unsigned int i = 0; i < blocks.size(); i++) {
      retval |= blocks[i]->line_search(alpha, qp_dx, qp_dy, rhs_residual);
    }
    return retval;
  }


  int UQPSumModelBlock::set_point(const Matrix& qp_x,
    const Matrix& qp_y,
    Real alpha) {
    int retval = 0;
    for (unsigned int i = 0; i < blocks.size(); i++) {
      retval |= blocks[i]->set_point(qp_x, qp_y, alpha);
    }
    return retval;
  }


  int UQPSumModelBlock::add_modelx_aggregate(Real& offset,
    Matrix& gradient) {
    int retval = 0;
    for (unsigned int i = 0; i < blocks.size(); i++) {
      retval |= blocks[i]->add_modelx_aggregate(offset, gradient);
    }
    return retval;
  }

  void UQPSumModelBlock::set_out(std::ostream* o, int pril) {
    CBout::set_out(o, pril);
    for (unsigned int i = 0; i < blocks.size(); i++) {
      blocks[i]->set_out(o, pril);
    }
  }

  void UQPSumModelBlock::set_cbout(const CBout* cb, int incr) {
    CBout::set_out(cb->get_out_ptr(), cb->get_print_level() + incr);
    for (unsigned int i = 0; i < blocks.size(); i++) {
      blocks[i]->set_cbout(this, 0);
    }
  }

  //---------------- for debugging purposes

  Matrix& UQPSumModelBlock::add_Bs(Matrix& qp_vec) const {
    for (unsigned int i = 0; i < blocks.size(); i++) {
      blocks[i]->add_Bs(qp_vec);
    }
    return qp_vec;
  }

  Matrix& UQPSumModelBlock::subtract_z(Matrix& dual_residual, bool with_step) const {
    for (unsigned int i = 0; i < blocks.size(); i++) {
      blocks[i]->subtract_z(dual_residual, with_step);
    }
    return dual_residual;
  }

}

