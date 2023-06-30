/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/BundleTerminator.hxx
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



#ifndef CONICBUNDLE_BUNDLETERMINATOR_HXX
#define CONICBUNDLE_BUNDLETERMINATOR_HXX


/**  @file BundleTerminator.hxx
    @brief Header declaring the classes ConicBundle::BundleTerminatorData, ConicBundle::BundleTerminator
    @version 1.0
    @date 2014-08-06
    @author Christoph Helmberg
*/


#include "matrix.hxx"
#include "clock.hxx"
#include "CBSolver.hxx"
#include "CBout.hxx"

namespace ConicBundle {

  /** @ingroup InternalBundleSolver
  */

  //@{

    /** @brief abstract interface for BundleTerminator providing the data needed for deciding on termination
     */
  class BundleTerminatorData {
  public:
    ///
    virtual ~BundleTerminatorData(); //implemented in BundleTerminator.cxx

    /// returns the number of calls to the oracle, i.e., to BundleProblem::eval_function()
    virtual CH_Matrix_Classes::Integer get_cntobjeval() const = 0;
    /// returns the number of reevalutions in center points (if their function value violates the most recent subgradient inequality)
    virtual CH_Matrix_Classes::Integer get_sumrecomp() const = 0;
    /// returns the number of descent and null steps
    virtual CH_Matrix_Classes::Integer get_suminnerit() const = 0;
    /// returns the number of times, the call to a quadratic subproblem failed
    virtual CH_Matrix_Classes::Integer get_sumqpfails() const = 0;
    /// returns the upper bound on the function value in the current center of stability
    virtual CH_Matrix_Classes::Real get_center_objval() const = 0;
    /// returns the upper bound on the function value in the current candidate
    virtual CH_Matrix_Classes::Real get_cand_objval() const = 0;
    /// returns the model value in the current candidate
    virtual CH_Matrix_Classes::Real get_modelval() const = 0;
    /// returns the norm of the aggregate subgradient squared (for the norm dual to quadratic augmented term) 
    virtual CH_Matrix_Classes::Real get_aggr_dnormsqr() const = 0;
    /// returns the number of failed callls to model evaluations 
    virtual CH_Matrix_Classes::Integer get_summodelfails() const = 0;
    /// returns the number of times, the augmented model value could not be increased since the last descent step
    virtual CH_Matrix_Classes::Integer get_augvalfails() const = 0;
    /// returns the number of times, the augmented model value could not be increased over all descent and null steps
    virtual CH_Matrix_Classes::Integer get_sumaugvalfails() const = 0;
    /// returns the number of times, the oracle returend some error code since the last descent step
    virtual CH_Matrix_Classes::Integer get_oraclefails() const = 0;
    /// returns the number of times, the oracle returend some error code over all descent and null steps
    virtual CH_Matrix_Classes::Integer get_sumoraclefails() const = 0;
    /// returns a correction factor for the termination criterion if the quadratic term is big (the proximal term is strong)
    virtual CH_Matrix_Classes::Real get_term_corr() const = 0;

  };

  //@}

  class BundleSolver;
  class NBundleSolver;

  /** @ingroup InternalBundleSolver
  */

  //@{
    /** @brief basic class implementing termination criteria for BundleSolver, may also serve as base class for other termination criteria
     */
  class BundleTerminator : public CBout {
    friend class BundleSolver;
    friend class NBundleSolver;

  protected:
    /// termination precision
    CH_Matrix_Classes::Real termeps;
    /// if the value is positive, do not terminate due to sufficient precision as long as the squared dual norm of the aggregate exceeds this 
    CH_Matrix_Classes::Real aggregate_dnormsqr;

    /// if a computation time criterion is set, this pointer is also set and points at the clock to be used
    const CH_Tools::Clock* clockp;
    /// if clockp is not Null then this gives the upper time limit in Microseconds
    CH_Tools::Microseconds timelimit;

    /// upper limit on number returned in BundleTerminatorData::get_sumrecomp, <0 if no limit
    CH_Matrix_Classes::Integer recomplimit;

    /// upper limit on number returned in BundleTerminatorData::get_sumqpfails, <0 if no limit
    CH_Matrix_Classes::Integer qpfailslimit;

    /// upper limit on number returned in BundleTerminatorData::get_summodelfails, <0 if no limit
    CH_Matrix_Classes::Integer modelfailslimit;

    /// upper limit on number returned in BundleTerminatorData::get_sumaugvalfails, <0 if no limit
    CH_Matrix_Classes::Integer augvalfailslimit;

    /// upper limit on number returned in BundleTerminatorData::get_cntobjeval, <0 if no limit
    CH_Matrix_Classes::Integer objevallimit;

    /// upper limit on number returned in BundleTerminatorData::get_sumoraclefails, <0 if no limit
    CH_Matrix_Classes::Integer oraclefailslimit;

    /// result of last call to check_termination
    int terminated;


    /** @brief computes and returns a termination code for the BundleTerminationData passed to this

        The termination code is composed by setting various bits for signaling
        various reasons for termination. The values/bits mean the following:
        -  0    not terminated,
        -  1    precision achieved,
        -  2    timelimit exceeded
        -  4    recomplimit exceeded
        -  8    qpfailslimit exceeded
        - 16    modelfailslimit exceeded
        - 32    augvalfailslimit exceeded
        - 64    objevallimit exceeded
        - 128   oraclefailslimit exceeded
    */

    virtual int check_termination(BundleTerminatorData* sb) {
      terminated = 0;
      if (clockp)
        terminated |= 2 * (clockp->time() >= timelimit);
      if (recomplimit >= 0)
        terminated |= 4 * (recomplimit <= sb->get_sumrecomp());
      if (qpfailslimit >= 0)
        terminated |= 8 * (qpfailslimit <= sb->get_sumqpfails());
      if (modelfailslimit >= 0)
        terminated |= 16 * (modelfailslimit <= sb->get_summodelfails());
      if (augvalfailslimit >= 0)
        terminated |= 32 * (augvalfailslimit <= sb->get_sumaugvalfails());
      if (objevallimit >= 0)
        terminated |= 64 * (objevallimit <= sb->get_cntobjeval());
      if (oraclefailslimit >= 0)
        terminated |= 128 * (oraclefailslimit <= sb->get_sumoraclefails());
      if ((sb->get_suminnerit() < 10) && (sb->get_aggr_dnormsqr() > 0.1)) return terminated;
      if ((aggregate_dnormsqr > 0) && (sb->get_aggr_dnormsqr() > aggregate_dnormsqr)) return terminated;
      CH_Matrix_Classes::Real center_objval = sb->get_center_objval();
      terminated |= (center_objval - sb->get_modelval() <=
        termeps * (CH_Matrix_Classes::abs(center_objval) + 1.) * sb->get_term_corr());
      return terminated;
    }

  public:
    /// sets the default parameter values
    virtual void set_defaults() {
      termeps = 1e-5; clockp = 0; timelimit.set_infinity(true);
      recomplimit = 100; qpfailslimit = 100; modelfailslimit = 100;
      augvalfailslimit = 10; oraclefailslimit = 10; objevallimit = -1;
      aggregate_dnormsqr = -1.;
    }

    /// resets @a terminated
    virtual void clear() {
      terminated = 0;
    }

    /// calls set_defaults() and clear()
    BundleTerminator(const CBout* cb = 0, int incr = -1) :CBout(cb, incr) { /*save_function=0;*/ set_defaults(); clear();
    }

    ///
    virtual ~BundleTerminator();  //implemented in bundleterminator.cxx

    /// set the termination precision (>0!) 
    virtual void set_termeps(CH_Matrix_Classes::Real teps) {
      if (teps > 0.) termeps = teps;
    }
    /// returns the current termination precision
    virtual CH_Matrix_Classes::Real get_termeps() const {
      return termeps;
    }
    /// set an upper bound for the dual norm squared of the aggregate at termination (or <=0 if no such bound is desired) 
    virtual void set_aggr_dnormsqr(CH_Matrix_Classes::Real sg) {
      aggregate_dnormsqr = sg;
    }
    /// returns the current bound for the dual norm squared of the aggregate
    virtual CH_Matrix_Classes::Real get_aggr_dnormsqr() const {
      return aggregate_dnormsqr;
    }
    /// set cp==0 for no timelimit, otherwise specify clock and microseconds
    virtual void set_timelimit(const CH_Tools::Clock* cp, CH_Tools::Microseconds tl) {
      clockp = cp; timelimit = tl;
    }
    /// returns the timelimit value
    virtual CH_Tools::Microseconds get_timelimit() const {
      return timelimit;
    }
    /// set upper bound on the value returned by BundleTerminatorData::get_sumrecomp, <0 if no limit
    virtual void set_recomplimit(CH_Matrix_Classes::Integer rl) {
      recomplimit = rl;
    }
    /// returns the current value of this parameter
    virtual CH_Matrix_Classes::Integer get_recomplimit() const {
      return recomplimit;
    }
    /// set upper bound on the value returned by BundleTerminatorData::get_sumqpfails, <0 if no limit
    virtual void set_qpfailslimit(CH_Matrix_Classes::Integer ql) {
      qpfailslimit = ql;
    }
    /// returns the current value of this parameter
    virtual CH_Matrix_Classes::Integer get_qpfailslimit() const {
      return qpfailslimit;
    }
    /// set upper bound on the value returned by BundleTerminatorData::get_summodelfails, <0 if no limit
    virtual void set_modelfailslimit(CH_Matrix_Classes::Integer ml) {
      modelfailslimit = ml;
    }
    /// returns the current value of this parameter
    virtual CH_Matrix_Classes::Integer get_modelfailslimit() const {
      return modelfailslimit;
    }
    /// set upper bound on the value returned by BundleTerminatorData::get_sumaugvalfails(), <0 if no limit
    virtual void set_augvalfailslimit(CH_Matrix_Classes::Integer al) {
      augvalfailslimit = al;
    }
    /// returns the current value of this parameter
    virtual CH_Matrix_Classes::Integer get_augvalfailslimit() const {
      return augvalfailslimit;
    }
    /// set upper bound on the value returned by BundleTerminatorData::get_cntobjeval, <0 if no limit
    virtual void set_objevallimit(CH_Matrix_Classes::Integer ol) {
      objevallimit = ol;
    }
    /// returns the current value of this parameter
    virtual CH_Matrix_Classes::Integer get_objevallimit() const {
      return objevallimit;
    }
    /// set upper bound on the value returned by BundleTerminatorData::get_sumoraclefails, <0 if no limit
    virtual void set_oraclefailslimit(CH_Matrix_Classes::Integer ol) {
      oraclefailslimit = ol;
    }
    /// returns the current value of this parameter
    virtual CH_Matrix_Classes::Integer get_oraclefailslimit() const {
      return oraclefailslimit;
    }

    /// return the termination code returned in the last call to check_termination()  
    virtual int get_terminated() const {
      return terminated;
    }
    /// reset the termination code to zero (this does not remove the reason for the termination; for this, set other bounds or clear the numbers supplied by BundleTerminationData)
    virtual void clear_terminated() {
      terminated = 0;
    }

    /// output an explanation string for the current termination code 
    virtual void print_status(std::ostream& o) const {
      o << "termination status: " << terminated;
      if (terminated == 0) {
        o << " (not terminated)" << std::endl; return;
      }
      if (terminated & 1) {
        o << ", relative precision criterion satisfied";
      }
      if (terminated & 2) {
        o << ", timelimit exceeded";
      }
      if (terminated & 4) {
        o << ", function reevaluation limit exceeded";
      }
      if (terminated & 8) {
        o << ", limit of QP failures exceeded";
      }
      if (terminated & 16) {
        o << ", limit of model failures exceeded";
      }
      if (terminated & 32) {
        o << ", limit of augmented model failures exceeded";
      }
      if (terminated & 64) {
        o << ", limit of calls to evaluation oracle exceeded";
      }
      if (terminated & 128) {
        o << ", limit of failed oracle calls exceeded";
      }
      o << std::endl;
    }

    /// output current parameter settings
    virtual std::ostream& save(std::ostream& o) const {
      o.precision(20);
      o << termeps << "\n" << timelimit << "\n" << recomplimit << "\n" << qpfailslimit << "\n" << modelfailslimit << "\n" << augvalfailslimit << "\n" << objevallimit << "\n" << oraclefailslimit << "\n" << aggregate_dnormsqr << "\n" << terminated << "\n";
      return o;
    }

    /// input new parameter settings, the clock must be set in addition!
    virtual std::istream& restore(std::istream& in) {
      in >> termeps;
      in >> timelimit;
      in >> recomplimit;
      in >> qpfailslimit;
      in >> modelfailslimit;
      in >> augvalfailslimit;
      in >> objevallimit;
      in >> oraclefailslimit;
      in >> aggregate_dnormsqr;
      in >> terminated;
      return in;
    }
  };

  //@}

}

#endif

