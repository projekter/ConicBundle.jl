/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/PSCPrimal.cxx
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



#include "PSCPrimal.hxx"

using namespace CH_Matrix_Classes;


namespace ConicBundle {


  BlockPSCPrimal::BlockPSCPrimal(const BlockPSCPrimal& pr, double factor) : PSCPrimal() {
    assert(factor >= 0.);
    Xdim = pr.Xdim;
    for (std::map<Integer, PSCPrimal*>::const_iterator i = pr.primal.begin(); i != pr.primal.end(); ++i) {
      primal[i->first] = dynamic_cast<PSCPrimal*>(i->second->clone_primal_data());
      assert(primal[i->first]);
      if (factor != 1.)
        primal[i->first]->scale_primal_data(factor);
    }

  }

  BlockPSCPrimal::BlockPSCPrimal(const Indexmatrix& Xd, const std::map<Integer, PSCPrimal*>& pr, double factor) : PSCPrimal() {
    assert(factor >= 0.);
    Xdim = Xd;
    for (std::map<Integer, PSCPrimal*>::const_iterator i = pr.begin(); i != pr.end(); ++i) {
      if (i->second == 0) continue;
      assert(i->first < Xdim.dim());
      primal[i->first] = i->second;
      if (factor != 1.)
        primal[i->first]->scale_primal_data(factor);
    }
  }

  BlockPSCPrimal::~BlockPSCPrimal() {
    for (std::map<Integer, PSCPrimal*>::iterator i = primal.begin(); i != primal.end(); ++i) {
      delete i->second;
    }
    primal.clear();
  }

  /// for each element aij in the support set aij=<P.row(i),P.row(j)> 
  int BlockPSCPrimal::assign_Gram_matrix(const Matrix& P) {
    assert(sum(Xdim) == P.rowdim());
    Integer start_row = 0;
    Integer cnt = 0;
    for (std::map<Integer, PSCPrimal*>::iterator i = primal.begin(); i != primal.end(); ++i) {
      while (cnt < i->first) {
        start_row += Xdim(cnt);
        cnt++;
      }
      i->second->assign_Gram_matrix(P.rows(Range(start_row, start_row + Xdim(cnt) - 1)));
    }
    return 0;
  }

  /** multiply this with myfactor and add itsfactor*it to this
      (it must also be a SparsePSCPrimal and on the same support) */
  int BlockPSCPrimal::aggregate_primal_data(const PrimalData& it, double factor) {
    assert(factor >= 0.);
    const BlockPSCPrimal* p = dynamic_cast<const BlockPSCPrimal*>(&it);
    if (p == 0) return 1;
    assert(norm2(Xdim - p->Xdim) < 0.1);
    for (std::map<Integer, PSCPrimal*>::iterator i = primal.begin(); i != primal.end(); ++i) {
      std::map<Integer, PSCPrimal*>::const_iterator j = p->primal.find(i->first);
      if ((j == p->primal.end()) || (j->second == 0)) return 1;
      i->second->aggregate_primal_data(*j->second, factor);
    }
    return 0;
  }

  /// multiply this with myfactor and add itsfactor*P*P^T to this
  int BlockPSCPrimal::aggregate_Gram_matrix(const Matrix& P, double factor) {
    assert(sum(Xdim) == P.rowdim());
    Integer start_row = 0;
    Integer cnt = 0;
    for (std::map<Integer, PSCPrimal*>::iterator i = primal.begin(); i != primal.end(); ++i) {
      while (cnt < i->first) {
        start_row += Xdim(cnt);
        cnt++;
      }
      i->second->aggregate_Gram_matrix(P.rows(Range(start_row, start_row + Xdim(cnt) - 1)), factor);
    }
    return 0;
  }

  PSCPrimal* BlockPSCPrimal::block(Integer i) const {
    std::map<Integer, PSCPrimal*>::const_iterator j = primal.find(i);
    if (j == primal.end()) return 0;
    return j->second;
  }


  /// multiply this with myfactor and add itsfactor*P*P^T to this
  int BlockPSCPrimal::scale_primal_data(double factor) {
    for (std::map<Integer, PSCPrimal*>::iterator i = primal.begin(); i != primal.end(); ++i) {
      i->second->scale_primal_data(factor);
    }
    return 0;
  }

  /// if compatible evaluate value=ip(*this,A.column[i])
  int BlockPSCPrimal::primal_ip(Real& value,
    const SparseCoeffmatMatrix& A,
    Integer col) const {
    if ((col < 0) || (col >= A.coldim()) || (Xdim.dim() != A.blockdim().dim()))
      return 1;
    const SparseCoeffmatVector* colA = A.column(col);
    value = 0.;
    if (colA == 0)
      return 0;
    for (Integer j = 0; j < Xdim.dim(); j++) {
      if (Xdim(j) != A.blockdim(j))
        return 1;
    }
    for (SparseCoeffmatVector::const_iterator it = colA->begin(); it != colA->end(); it++) {
      std::map<Integer, PSCPrimal*>::const_iterator pit = primal.find(it->first);
      if (pit != primal.end()) {
        const DensePSCPrimal* dp = dynamic_cast<const DensePSCPrimal*>(pit->second);
        if (dp) {
          value += it->second->ip(*dp);
          continue;
        }
        const SparsePSCPrimal* sp = dynamic_cast<const SparsePSCPrimal*>(pit->second);
        if (sp) {
          if (it->second->support_in(*sp) == 0) {
            return 1;
          }
          value += it->second->ip(*sp);
          continue;
        }
        return 1;
      }
    }
    return 0;
  }

  /// returns a newly generated identical Object
  PrimalData* BlockPSCPrimal::clone_primal_data() const {
    return new BlockPSCPrimal(*this);
  }




}

