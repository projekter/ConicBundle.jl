/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/QPKKTSolverComparison.hxx
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


#ifndef CONICBUNDLE_QPKKTSOLVERCOMPARISON_HXX
#define CONICBUNDLE_QPKKTSOLVERCOMPARISON_HXX

/**  @file QPKKTSolverComparison.hxx
    @brief Header declaring the class ConicBundle::QPKKTSolverComparision
    @version 1.0
    @date 2021-03-30
    @author Christoph Helmberg
*/

#include <string>
#include "QPKKTSolverObject.hxx"
#include "clock.hxx"

namespace ConicBundle {


  /** @ingroup ConstrainedQPSolver
   */
   //@{


   /** @brief Used for collecting statics in QPKKTSolverComparison: For each bundle subproblem there is a block QPKKT_ProbStats, which holds for each KKT system a block of QPKKT_KKTStats, which holds for each solver a block QPKKT_SolverStats (this).

      The violation values typically refer to the solution of the corrector.
    */
  class QPKKT_SolverStats {
  public:
    CH_Tools::Microseconds preptime;  ///< time spent in setting up the system
    CH_Tools::Microseconds predtime;  ///< time spent in solving the predictor
    CH_Tools::Microseconds corrtime;  ///< time spent in solving the corrector
    CH_Matrix_Classes::Integer prepnmult; ///< matrix multiplications in setting up the system
    CH_Matrix_Classes::Integer prednmult; ///< matrix multiplications in solving the predictor 
    CH_Matrix_Classes::Integer corrnmult; ///< matrix multiplications in solving the correcotr
    CH_Matrix_Classes::Real cond;   ///< estimate of the condition number
    CH_Matrix_Classes::Integer rank;  ///< the rank used in the preconditioner
    CH_Matrix_Classes::Real Hviol;   ///< norm of the residual for the H block  
    CH_Matrix_Classes::Real Aviol;   ///< norm of the residual for the A block
    CH_Matrix_Classes::Real Bviol;   ///< norm of the residual for the B block
    CH_Matrix_Classes::Real Cviol;   ///< norm of the residual for the C block
    CH_Matrix_Classes::Real sysviol; ///< norm of the residual of the entire KKT system

    /// constructor
    QPKKT_SolverStats() {
      preptime = predtime = corrtime = 0;
      prepnmult = prednmult = corrnmult = 0;
      cond = -1.;
      rank = -1;
      Hviol = Aviol = Bviol = Cviol = sysviol = 0.;
    }

    /// destructor
    ~QPKKT_SolverStats() {
    }

    /// (file-)output
    friend std::ostream& operator<<(std::ostream& out, const QPKKT_SolverStats& s);

    /// (file-)input
    friend std::istream& operator>>(std::ostream& in, QPKKT_SolverStats& s);
  };

  /// output QPKKT_SolverStats (for files)
  std::ostream& operator<<(std::ostream& out, const QPKKT_SolverStats& s);

  /// input QPKKT_SolverStats (for files)
  std::istream& operator>>(std::ostream& in, QPKKT_SolverStats& s);

  //@}


/** @ingroup ConstrainedQPSolver
 */
 //@{

 /** @brief Used for collecting statics in QPKKTSolverComparison: For each bundle subproblem there is a block QPKKT_ProbStats, which holds for each KKT system a block of QPKKT_KKTStats (this), which holds for each solver a block QPKKT_SolverStats.

  */
  class QPKKT_KKTStats {
  public:
    CH_Matrix_Classes::Real prec;       ///< this is the (minimal) precision requirement in the calls
    CH_Matrix_Classes::Real mu;         ///< this is the (minimal) mu>0 value used in the calls 
    std::vector<QPKKT_SolverStats> sdata; ///< one entry per solver

    /// constructor
    QPKKT_KKTStats() {
      prec = mu = 0.;
    }

    /// destructor
    ~QPKKT_KKTStats() {
    }

    /// appends the statsitics data of this KKT system as the last column of the repsective matrices
    int append_col(CH_Matrix_Classes::Matrix& muvals, ///< the barrier parameter, 
      CH_Matrix_Classes::Matrix& prepsecs, ///< preparation time (one row for each solver)
      CH_Matrix_Classes::Matrix& predsecs, ///< predictor time (one row for each solver)
      CH_Matrix_Classes::Matrix& corrsecs, ///< corrector time (one row for each solver)
      CH_Matrix_Classes::Indexmatrix& predcalls, ///< predictor matrix vector multiplications (one row for each solver)
      CH_Matrix_Classes::Indexmatrix& corrcalls, ///< corrector matrix vector multiplications (one row for each solver)
      CH_Matrix_Classes::Matrix& cond, ///< condition number estimate (one row for each solver)
      CH_Matrix_Classes::Indexmatrix& pccols, ///< rank of the predictor (one row for each solver)
      CH_Matrix_Classes::Matrix& sysviol) ///< norm of the system residual for corrector (one row for each solver)
    {
      CH_Matrix_Classes::Integer ind = muvals.coldim();
      muvals.concat_right(mu);
      assert(prepsecs.rowdim() == CH_Matrix_Classes::Integer(sdata.size()));
      assert(predsecs.rowdim() == CH_Matrix_Classes::Integer(sdata.size()));
      assert(corrsecs.rowdim() == CH_Matrix_Classes::Integer(sdata.size()));
      assert(predcalls.rowdim() == CH_Matrix_Classes::Integer(sdata.size()));
      assert(corrcalls.rowdim() == CH_Matrix_Classes::Integer(sdata.size()));
      assert(cond.rowdim() == CH_Matrix_Classes::Integer(sdata.size()));
      assert(pccols.rowdim() == CH_Matrix_Classes::Integer(sdata.size()));
      assert(sysviol.rowdim() == CH_Matrix_Classes::Integer(sdata.size()));
      prepsecs.enlarge_right(1, 0.);
      predsecs.enlarge_right(1, 0.);
      corrsecs.enlarge_right(1, 0.);
      predcalls.enlarge_right(1, CH_Matrix_Classes::Integer(0));
      corrcalls.enlarge_right(1, CH_Matrix_Classes::Integer(0));
      cond.enlarge_right(1, CH_Matrix_Classes::Integer(0));
      pccols.enlarge_right(1, CH_Matrix_Classes::Integer(0));
      sysviol.enlarge_right(1, CH_Matrix_Classes::Integer(0));
      for (CH_Matrix_Classes::Integer i = 0; i < CH_Matrix_Classes::Integer(sdata.size()); i++) {
        prepsecs(i, ind) += sdata[unsigned(i)].preptime;
        predsecs(i, ind) += sdata[unsigned(i)].predtime;
        corrsecs(i, ind) = sdata[unsigned(i)].corrtime;
        predcalls(i, ind) = sdata[unsigned(i)].prednmult;
        corrcalls(i, ind) = sdata[unsigned(i)].corrnmult;
        cond(i, ind) = sdata[unsigned(i)].cond;
        pccols(i, ind) = sdata[unsigned(i)].rank;
        sysviol(i, ind) = sdata[unsigned(i)].sysviol;
      }
      return 0;
    }

    /// adds the respecitve values of the solvers to last column. It is used for collecting cummlative data for each subproblem per solver. Columns are subproblems, rows are solvers
    int add_col(CH_Matrix_Classes::Matrix& prepsecs, ///< preparation time (one row for each solver)
      CH_Matrix_Classes::Matrix& predsecs, ///< predictor time (one row for each solver)
      CH_Matrix_Classes::Matrix& corrsecs, ///< corrector time (one row for each solver)
      CH_Matrix_Classes::Indexmatrix& predcalls, ///< predictor matrix vector multiplications (one row for each solver)
      CH_Matrix_Classes::Indexmatrix& corrcalls ///< corrector matrix vector multiplications (one row for each solver)
    ) {
      CH_Matrix_Classes::Integer ind = prepsecs.coldim() - 1;
      for (CH_Matrix_Classes::Integer i = 0; i < CH_Matrix_Classes::Integer(sdata.size()); i++) {
        prepsecs(i, ind) += sdata[unsigned(i)].preptime;
        predsecs(i, ind) += sdata[unsigned(i)].predtime;
        corrsecs(i, ind) += sdata[unsigned(i)].corrtime;
        predcalls(i, ind) += sdata[unsigned(i)].prednmult;
        corrcalls(i, ind) += sdata[unsigned(i)].corrnmult;
      }
      return 0;
    }

    /// (file-)output
    friend std::ostream& operator<<(std::ostream& out, const QPKKT_KKTStats& k);

    /// (file-)input
    friend std::istream& operator>>(std::ostream& in, QPKKT_KKTStats& k);

  };

  /// output QPKKT_KKZStats (for files)
  std::ostream& operator<<(std::ostream& out, const QPKKT_KKTStats& k);

  /// input QPKKT_KKZStats (for files)
  std::istream& operator>>(std::istream& in, QPKKT_KKTStats& k);

  //@}


/** @ingroup ConstrainedQPSolver
 */
 //@{

 /** @brief Used for collecting statics in QPKKTSolverComparison: For each bundle subproblem there is a block QPKKT_ProbStats (this), which holds for each KKT system a block of QPKKT_KKTStats, which holds for each solver a block QPKKT_SolverStats.

  */
  class QPKKT_ProbStats {
  public:
    CH_Matrix_Classes::Integer Qdim;       ///< the order of the H block (the quadratic/proximal term)
    CH_Matrix_Classes::Integer Vdim;       ///< the low-rank rank of the H matrix
    CH_Matrix_Classes::Integer Arowdim;    ///< rowdim of the constraints
    CH_Matrix_Classes::Integer Aeqdim;     ///< with this number of equations
    CH_Matrix_Classes::Integer Bdim;       ///< number of vectors in the bundle
    CH_Matrix_Classes::Integer Cdim;       ///< number of constraints in the model
    std::vector<CH_Tools::Microseconds> inittime; ///< time spent in initializing the solver for a subproblem
    std::vector<QPKKT_KKTStats> kktdata; ///< one entry per KKT system

    /// constructor
    QPKKT_ProbStats() {
      Qdim = Vdim = Arowdim = Aeqdim = Bdim = Cdim = 0;
    }

    /// destructor
    ~QPKKT_ProbStats() {
    }

    /// for collecting the statistics for an interval of barrier parameter values, this first determines the number of KKT systems stored with barrier parameter in this range
    CH_Matrix_Classes::Integer cnt_mu_cols(CH_Matrix_Classes::Real lbmu,
      CH_Matrix_Classes::Real ubmu) {
      CH_Matrix_Classes::Integer sumval = 0;
      for (unsigned int i = 0; i < kktdata.size(); i++)
        if ((lbmu <= kktdata[i].mu) && (kktdata[i].mu < ubmu))
          sumval++;
      return sumval;
    }

    /// collects the statistics for an interval of barrier parameter values by appending a column for each KKT system (the rows then correspond the solvers)
    int get_mu_stats(CH_Matrix_Classes::Real lbmu,
      CH_Matrix_Classes::Real ubmu,
      CH_Matrix_Classes::Indexmatrix& dims,
      CH_Matrix_Classes::Matrix& mu,
      CH_Matrix_Classes::Matrix& prepsecs,
      CH_Matrix_Classes::Matrix& predsecs,
      CH_Matrix_Classes::Matrix& corrsecs,
      CH_Matrix_Classes::Indexmatrix& predcalls,
      CH_Matrix_Classes::Indexmatrix& corrcalls,
      CH_Matrix_Classes::Matrix& cond,
      CH_Matrix_Classes::Indexmatrix& pccols,
      CH_Matrix_Classes::Matrix& sysviol) {
      CH_Matrix_Classes::Integer cnt = 0;
      for (unsigned int i = 0; i < kktdata.size(); i++) {
        if ((lbmu <= kktdata[i].mu) && (kktdata[i].mu < ubmu)) {
          kktdata[i].append_col(mu, prepsecs, predsecs, corrsecs, predcalls, corrcalls, cond, pccols, sysviol);
          cnt++;
        }
      }
      assert(dims.coldim() + cnt == mu.coldim());
      CH_Matrix_Classes::Integer startcol = dims.coldim();
      if (dims.coldim() == 0)
        dims.init(6, cnt, CH_Matrix_Classes::Integer(0));
      else {
        assert(dims.rowdim() == 6);
        dims.enlarge_right(cnt, CH_Matrix_Classes::Integer(0));
      }
      for (CH_Matrix_Classes::Integer i = 0; i < cnt; i++) {
        dims(0, startcol + i) = Qdim;
        dims(1, startcol + i) = Vdim;
        dims(2, startcol + i) = Arowdim;
        dims(3, startcol + i) = Aeqdim;
        dims(4, startcol + i) = Bdim;
        dims(5, startcol + i) = Cdim;
      }
      return 0;
    }

    ///collects the cummulative values for each quadratic (bundle) subproblem. For each problem a new column is appended; in dims the rows give the dimensions stored in QPKKT_ProbStats; iterations gives the number of KKT systems in each subproblem; lastmu gives the value of the barrier parameter used in the last KKT system of the subproblem; for the other parameters the rows correspond to the cummulative values for each solvers 
    int get_prob_stats(CH_Matrix_Classes::Indexmatrix& dims, ///< per column/subproblem: order of H, lowrank of H, # rows of A, # equalities in A, # rows of B, # rows of C
      CH_Matrix_Classes::Indexmatrix& iterations, ///< number of KKTsystem of each subproblem
      CH_Matrix_Classes::Matrix& lastmu, ///< value of the barrier parameter in the last KKT system of the subproblem
      CH_Matrix_Classes::Matrix& prepsecs, ///< cummulative time spent in setting up the subproblem and systems (one row per solver)
      CH_Matrix_Classes::Matrix& predsecs, ///< cummulative time spent in solving the predictors (one row per solver)
      CH_Matrix_Classes::Matrix& corrsecs, ///< cummulative time spent in solving the correctors (one row per solver)
      CH_Matrix_Classes::Indexmatrix& predcalls, ///< cummulative number of matrix vector multiplications in the predictors (one row per solver)
      CH_Matrix_Classes::Indexmatrix& corrcalls ///< cummulative number of matrix vector multiplications in the correctors (one row per solver)
    ) {
      assert(dims.coldim() == iterations.coldim());
      assert(dims.coldim() == lastmu.coldim());
      assert(dims.coldim() == prepsecs.coldim());
      assert(dims.coldim() == predsecs.coldim());
      assert(dims.coldim() == corrsecs.coldim());
      assert(dims.coldim() == predcalls.coldim());
      assert(dims.coldim() == corrcalls.coldim());
      CH_Matrix_Classes::Integer ind = dims.coldim();
      if (dims.coldim() == 0)
        dims.init(6, 1, CH_Matrix_Classes::Integer(0));
      else {
        assert(dims.rowdim() == 6);
        dims.enlarge_right(1, CH_Matrix_Classes::Integer(0));
      }
      dims(0, ind) = Qdim;
      dims(1, ind) = Vdim;
      dims(2, ind) = Arowdim;
      dims(3, ind) = Aeqdim;
      dims(4, ind) = Bdim;
      dims(5, ind) = Cdim;
      iterations.concat_right(CH_Matrix_Classes::Integer(kktdata.size()));
      lastmu.concat_right(kktdata[kktdata.size() - 1].mu);
      assert(inittime.size() == unsigned(prepsecs.rowdim()));
      prepsecs.enlarge_right(1, 0.);
      for (unsigned int i = 0; i < inittime.size(); i++)
        prepsecs(CH_Matrix_Classes::Integer(i), ind) += inittime[i];
      predsecs.enlarge_right(1, 0.);
      corrsecs.enlarge_right(1, 0.);
      predcalls.enlarge_right(1, 0.);
      corrcalls.enlarge_right(1, 0.);
      for (unsigned int i = 0; i < kktdata.size(); i++) {
        kktdata[i].add_col(prepsecs, predsecs, corrsecs, predcalls, corrcalls);
      }
      return 0;
    }


    /// (file-)output
    friend std::ostream& operator<<(std::ostream& out, const QPKKT_ProbStats& p);

    /// (file-)input
    friend std::istream& operator>>(std::istream& in, QPKKT_ProbStats& p);
  };

  /// output QPKKT_ProbStats (for files)
  std::ostream& operator<<(std::ostream& out, const QPKKT_ProbStats& p);

  /// input QPKKT_ProbStats (for files)
  std::istream& operator>>(std::istream& in, QPKKT_ProbStats& p);

  //@}



/** @ingroup ConstrainedQPSolver
 */
 //@{

 /** @brief This is a pseudosolver designed for producing comparative statistics on the performance of mainly iterative solvers for the interior point KKT systems of QPSolver.

     Use add_solver() to enter the solvers that are to be compared. The first solver entered is used for actually returning the computations as if this first solver would be the true solver. The others work on exactly the same KKT system and model data via cloning.

     The statistics are stored as follows: for each bundle subproblem there is a block QPKKT_ProbStats, which holds for each KKT system a block of QPKKT_KKTStats, which holds for each solver a block QPKKT_SolverStats.

  */

  class QPKKTSolverComparison : public QPKKTSolverObject {
  private:
    std::vector<QPKKTSolverObject*> solver;  ///< solvers add via add_solver()
    std::vector<std::string> solvername;    ///< name of the solver for output
    std::vector<QPKKT_ProbStats> probdata;  ///< for each bundle subproblem there is a block QPKKT_ProbStats (stored here), which holds for each KKT system a block of QPKKT_KKTStats, which holds for each solver a block QPKKT_SolverStats.
    std::vector<QPModelBlockObject*> model; ///< a seperate clone of the model for each solver in order to avoid side effects

    QPSolverProxObject* testHp; ///< the proximal term/quadratic cost matrix; this points to external object
    const CH_Matrix_Classes::Sparsemat* testA; ///< the constraint matrix; this points to external object
    CH_Matrix_Classes::Matrix testKKTdiagx; ///< stores the input diagonal for the H block
    CH_Matrix_Classes::Matrix testKKTdiagy; ///< stores the input diagonal for the A block
    CH_Matrix_Classes::Real testHfactor; ///< stores a factor for testHp, usually ==1. unless a second order cone variant is used
    CH_Matrix_Classes::Matrix testprimalrhs; ///< stores the input rhs for the A block
    CH_Matrix_Classes::Matrix testdualrhs; ///< stores the input rhs for the H block

    CH_Tools::Clock clock; ///< for taking the time

    CH_Matrix_Classes::Matrix insolx;   ///< stores the input solution of the H block (e.g. of predictor)
    CH_Matrix_Classes::Matrix insoly;  ///< stores the input solution of the A block (e.g. of predictor)
    CH_Matrix_Classes::Matrix outsolx; ///< stores the output solution of the H block (first solver)
    CH_Matrix_Classes::Matrix outsoly; ///< stores the output solution of the A block (first solver)

    /// computes and stores the norms of the residualy of the KKT system
    int violation(CH_Matrix_Classes::Matrix& violvec, const CH_Matrix_Classes::Matrix& solx, const CH_Matrix_Classes::Matrix& soly, QPModelBlockObject* inmodel, QPKKT_SolverStats& stats);

  public:
    /// reset data to empty
    virtual void clear();


    /// default constructor
    QPKKTSolverComparison(CBout* cb = 0, int cbinc = -1);

    /// virtual destructor
    virtual ~QPKKTSolverComparison();

    /// the first solver added is the reference solver
    virtual int add_solver(QPKKTSolverObject* solver, const char* name);

    /// returns 1 if this class is not applicable in the current data situation, otherwise it stores the data pointers and these need to stay valid throught the use of the other routines but are not deleted here
    virtual int QPinit_KKTdata(QPSolverProxObject* Hp, ///< may not be be NULL 
      QPModelBlockObject* model, ///< may be NULL
      const CH_Matrix_Classes::Sparsemat* A, ///< may be NULL
      const CH_Matrix_Classes::Indexmatrix* eq_indices ///< if not NULL these rows of A correspond to equations
    );

    /// set up the primal dual KKT system for being solved for predictor and corrector rhs in QPsolve_KKTsystem
    virtual int QPinit_KKTsystem(const CH_Matrix_Classes::Matrix& KKTdiagx,
      const CH_Matrix_Classes::Matrix& KKTdiagy,
      CH_Matrix_Classes::Real Hfactor,
      CH_Matrix_Classes::Real prec,
      QPSolverParameters* params);

    /// solve the KKTsystem to precision prec for the given right hand sides that have been computed for the value rhsmu of the barrier parameter and in which a rhscorr fraction (out of [0,1] of the corrector term have been included; in iterative solvers solx and soly may be used as starting points
    virtual int QPsolve_KKTsystem(CH_Matrix_Classes::Matrix& solx,
      CH_Matrix_Classes::Matrix& soly,
      const CH_Matrix_Classes::Matrix& primalrhs,
      const CH_Matrix_Classes::Matrix& dualrhs,
      CH_Matrix_Classes::Real rhsmu,
      CH_Matrix_Classes::Real rhscorr,
      CH_Matrix_Classes::Real prec,
      QPSolverParameters* params);

    /// for judging violation this returns (an estimate of) the norm of the H-row in the latest system
    virtual CH_Matrix_Classes::Real QPget_blockH_norm();

    /// for judging violation this returns (an estimate of) the norm of the A-row in the latest system
    virtual CH_Matrix_Classes::Real QPget_blockA_norm();

    //------- for later statistical evaluation purposes

    /// returns the names of the solvers
    std::vector<std::string>& get_solvernames() {
      return solvername;
    }

    /// returns the stored statistical data
    std::vector<QPKKT_ProbStats>& get_probdata() {
      return probdata;
    }

    /// return those data columns (each a KKT system; columns are more efficient to append than lines) that fall into the given lower and upper bounds on mu
    int get_mu_stats(CH_Matrix_Classes::Real lbmu, ///< get KKT systems with last barrier paremeter at least this lower bound
      CH_Matrix_Classes::Real ubmu, ///< get KKT systems with last barrier paremeter at most this upper bound
      CH_Matrix_Classes::Indexmatrix& dims, ///< order of H, lowrank of H,  # rows of A, # equalities in A, # rows of B, # rows of C
      CH_Matrix_Classes::Matrix& mu, ///< value of the barrier parameter
      CH_Matrix_Classes::Matrix& prepsecs, ///< time spent in setting up the subproblem and systems (one row per solver)
      CH_Matrix_Classes::Matrix& predsecs, ///< time spent in solving the predictors (one row per solver)
      CH_Matrix_Classes::Matrix& corrsecs, ///< time spent in solving the correctors (one row per solver)
      CH_Matrix_Classes::Indexmatrix& predcalls, ///< number of matrix vector multiplications in the predictors (one row per solver)
      CH_Matrix_Classes::Indexmatrix& corrcalls, ///< number of matrix vector multiplications in the correctors (one row per solver)
      CH_Matrix_Classes::Matrix& cond, ///< condition number of the system (one row per solver)
      CH_Matrix_Classes::Indexmatrix& pccols, ///< number of columns in the low rank preconditioners
      CH_Matrix_Classes::Matrix& sysviol ///< norm of the (corrector) residual of the entire KKT system
    );

    /// return one data column per subproblem (more efficient to append than lines) with the sum of the time/calls/etc. 
    int get_prob_stats(CH_Matrix_Classes::Indexmatrix& dims,///< per column/subproblem: order of H, lowrank of H, # rows of A, # equalities in A, # rows of B, # rows of C
      CH_Matrix_Classes::Indexmatrix& iterations, ///< number of KKTsystem of each subproblem
      CH_Matrix_Classes::Matrix& lastmu, ///< value of the barrier parameter in the last KKT system of the subproblem
      CH_Matrix_Classes::Matrix& prepsecs, ///< cummulative time spent in setting up the subproblem and systems (one row per solver)
      CH_Matrix_Classes::Matrix& predsecs, ///< cummulative time spent in solving the predictors (one row per solver)
      CH_Matrix_Classes::Matrix& corrsecs, ///< cummulative time spent in solving the correctors (one row per solver)
      CH_Matrix_Classes::Indexmatrix& predcalls, ///< cummulative number of matrix vector multiplications in the predictors (one row per solver)
      CH_Matrix_Classes::Indexmatrix& corrcalls ///< cummulative number of matrix vector multiplications in the correctors (one row per solver)
    );

    /// (file-) output
    friend std::ostream& operator<<(std::ostream& out, const QPKKTSolverComparison& q);

    /// (file-) input
    friend std::istream& operator>>(std::istream& in, QPKKTSolverComparison& q);
  };


  /// output the collected statistical data for this  QPKKTSolverComparison (for files)
  std::ostream& operator<<(std::ostream& out, const QPKKTSolverComparison& q);

  /// input the stored statistical data of this  QPKKTSolverComparison (for files)
  std::istream& operator>>(std::istream& in, QPKKTSolverComparison& q);

  //@}



}

#endif

