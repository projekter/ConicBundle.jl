/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/SOCSupportFunction.hxx
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



#ifndef CONICBUNDLE_SOCSUPPORTFUNCTION_HXX
#define CONICBUNDLE_SOCSUPPORTFUNCTION_HXX

/**  @file SOCSupportFunction.hxx
    @brief Header declaring the classes ConicBundle::SOCSupportFunction, ConicBundle::AMFMinorantExtender (part of an implementation of ConicBundle::FunctionModel)
    @version 1.0
    @date 2017-09-30
    @author Christoph Helmberg
*/

//------------------------------------------------------------

#include <map>
#include "SOCSupportModification.hxx"
#include "SOCOracle.hxx"

//------------------------------------------------------------

namespace ConicBundle {

  /**@defgroup implemented_soc_oracle implemention of a SOCOracle (SOCSupportFunction)
    @brief SOCSupportFunction is an implementation of ConicBundle::SOCOracle for the minimization of the support function over the second order cone with x0=1, which may be used for Lagrangian relaxation of linear programs over the second order cone.

     The class SOCSupportFunction implements a general purpose version of SOCOracle for minimizing the support function

     \f[f:\mathbf{R}^n\to\mathbf{R},\quad f(\tilde c)=\max\{\tilde c^\top{1 \choose \bar x}\colon 1\ge\|\bar x\|\}\f]

     When adding this function to the solver by MatrixCBsolver::add_function(),
     one may specify a weight \f$\gamma>0\f$, a ConicBundle::FunctionTask, and an
     AffineFunctionTransformation \f$F\colon \mathbf{R}^m\to\mathbf{R}^n,
     F(y)=c+Ay\f$ for given \f$c\in\mathbf{R}^n\f$ and
     \f$A\in\mathbf{R}^{n\times m}\f$. Depending on the
     ConicBundle::FunctionTask, the solver actually minimizes the following:

     - FunctionTask::ObjectiveFunction

       \f[ f_{\gamma,F}(y)=\gamma f(F(y))=\gamma\max\{(c+Ay)^\top {1 \choose \bar x}\colon 1\ge\|\bar x\|\},\f]

       If in setting up the groundset within MatrixCBSolver the \f$y\f$ variables
       are introduced with a linear cost vector \f$b\mathbf{R}^m\f$ and
       appropriate sign constraints, this corresponds to Lagrangian relaxation
       of the linear constraints in the second order cone program

       \f[ \mbox{maximize } c^\top x
           \mbox{ subject to } A^\top x \begin{array}{c}\le\\=\\\ge\end{array} b,
           x={x_0 \choose \bar x}\mbox{ with }\|\bar x\|\le x_0=\gamma
       \f]

     - FunctionTask::ConstantPenaltyFunction

       \f[ f^+_{\gamma,F}(y)=\gamma\max\{0,f(F(y))\},\f]

       or, as above in the Lagrangean relaxation setting

       \f[ \mbox{maximize } c^\top x
           \mbox{ subject to } A^\top x \begin{array}{c}\le\\=\\\ge\end{array} b,
           x={x_0 \choose \bar x}\mbox{ with }\|\bar x\|\le x_0\le\gamma
       \f]

     - FunctionTask::AdaptivePenaltyFunction

       \f[ f^+_{\infty,F}(y)=\gamma\max\{0,f(F(y))\}\mbox{ with }\gamma\to \infty,\f]

       or, as above in the Lagrangean relaxation setting

       \f[ \mbox{maximize } c^\top x
           \mbox{ subject to } A^\top x \begin{array}{c}\le\\=\\\ge\end{array} b,
           x={x_0 \choose \bar x}\mbox{ with }\|\bar x\|\le x_0\le\gamma\mbox{ and }\gamma\to \infty
       \f]

     The standard use of a SOCSupportFunction is to initialize it on
     construction with dimension>=1 and then to add this function to the
     MatrixCBSolver. If dynamic changes to this SOCSupportFunction are required
     afterwards (only adding, deleting, rearranging coordinates with index != 0
     is allowed), use the class ConicBundle::SOCSupportModification within the
     corresponding problem modification routines of the MatrixCBSolver
     interface.

     The minorant generated for the input cost vector \f$\tilde c\f$ is a
     maximizing element of the second order cone, so it is itself the
     basic primal data and no extra primal data is needed.

     For facilitating input and output, SOCSupportFunction offers
     - SOCSupportFunction::print_problem_data() that outputs the full function description
       so that it can be read again by read_problem_data
     - SOCSupportFunction::read_problem_data() reads the problem data as output by print_problem_data()
     - SOCSupportFunction::set_out() and SOCSupportFunction::set_cbout() work as described in ConicBundle::CBout

  */
  //@{

    ///
  class SOCSupportFunction;

  /** @brief Implementation of MinorantExtender for SOCSupportFunction

      This object will be returned as an object on the heap
      by SOCSupportFunction::apply_modification
      and will be deleted by ConicBundle after use.
      Its purpose is to fill up further coordinates of the minorant
      (if this is possible).

  */

  class SOCSupportMinorantExtender : public CBout, public MinorantExtender {
  private:
    /// the oracle this MinorantExtender was generated by, needed for retrieving problem data
    SOCSupportFunction* fun;

  public:
    ///the SOCSupportFunction pointed to has to be valid for the livetime of this object
    SOCSupportMinorantExtender(SOCSupportFunction* fun);

    ///
    ~SOCSupportMinorantExtender();

    /*@brief called by ConicBundle to update internal Minorant objects, has to return 0 on success

        for each relevant index i the minorant value is set to the projection
  of 0 onto the interval [lower_bound(i),upper_bound(i)]

        @param[in,out] minorant  (Minorant&)
            it holds a (possibly aggregated) minorant that was generated
            from minorants returned by oracle calls, e.g. as in
      FunctionOracle::evaluate()

        @param[in] n_coords (int)
            the number of coordinate positions that have to be filled in

        @param[out] new_subgradient_values  (DVector &)
      the indices of these coordinate positions (sorted in
      strictly increasing order)

  @return
           -  0 on success,
           -  1 if extension/update is impossible

     */
    int extend(Minorant& minorant, int n_coords, const int* indices);

  };



  /**@brief general purpose implementation of SOCOracle as explained in \ref implemented_soc_oracle
   */

  class SOCSupportFunction : public SOCOracle, public CBout {
  private:
    /// dimension of the second order cone (for consistency checks)
    CH_Matrix_Classes::Integer socdim;

    /// applies the SOCSupportModfication mod to the current function
    int apply_modification(const SOCSupportModification& mod);

  public:
    /**@name Initialization */
    //@{

    /// initialize with dimension >= 1 (and output options)
    SOCSupportFunction(CH_Matrix_Classes::Integer socdim, const CBout* cb = 0, int incr = -1);

    ///
    ~SOCSupportFunction() {
    }

    //@}

    //----------- Oracle Implementation of SOCOracle ----------
    /**@name Implementations of SOCOracle routines */
    //@{

    /// see SOCOracle::generate_minorant()
    Minorant* generate_minorant(const CH_Matrix_Classes::Matrix& SOCvec);

    /// see SOCOracle::extract_SOCvector()
    int extract_SOCvector(CH_Matrix_Classes::Matrix& SOCvec, const Minorant* SOCminorant);

    /// see SOCOracle::projection()
    int projection(CH_Matrix_Classes::Matrix& offset,
      CH_Matrix_Classes::Matrix& coeffs,
      const CH_Matrix_Classes::Matrix& bar_P,
      const CH_Matrix_Classes::Indexmatrix* index_subset = 0);

    /// see SOCOracle::evaluate()
    virtual
      int
      evaluate
      (
        const CH_Matrix_Classes::Matrix& current_point,
        const CH_Matrix_Classes::Real relprec,
        CH_Matrix_Classes::Real& SOC_value,
        CH_Matrix_Classes::Matrix& SOC_vector,
        SOCPrimalExtender*& primal_extender
      );

    /// see SOCOracle::evaluate()
    virtual
      int
      evaluate_projection
      (
        const CH_Matrix_Classes::Matrix& current_point,
        const CH_Matrix_Classes::Matrix& P,
        const CH_Matrix_Classes::Real relprec,
        CH_Matrix_Classes::Real& projected_SOC_value
      );


    /** @brief see SOCOracle::apply_modification() for the general use, here oracle_modification has a special role if it can be cast to an AMFModification

  if oracle_modification cannot be cast to an SOCSupportModification it is assumed that all append modifications amount to have already been carried out on *this seperately before this routine is called. In particular, it is only checked whether the new dimension matches the one given by oracle_modification, the old dimension is ignored. If this does not hold, the routine stops with an error. Otherwise it checks the other stuff as if a suitable SOCSupportModification has just been executed.
     */
    virtual
      int
      apply_modification
      (
        const OracleModification& oracle_modification,
        const CH_Matrix_Classes::Matrix* new_center,
        const CH_Matrix_Classes::Matrix* old_center,
        bool& discard_objective_in_center,
        bool& discard_model,
        bool& discard_aggregates,
        MinorantExtender*& minorant_extender
      );


    /// see SOCOracle::check_correctness() (true only needed for debugging)
    virtual
      bool
      check_correctness() const {
      return true;
    }
    //@}


    //------------------  routines for querying the problem ----------------
    /**@name routines for querying data of the problem */
    //@{

    /// returns the dimension of the second order cone
    CH_Matrix_Classes::Integer get_socdim() {
      return socdim;
    }


    //@}


   //------------------  routines for IO ----------------
    /**@name routines for supporting input and output */
    //@{

    /// see ConicBundle::CBout
    void  set_out(std::ostream* o = 0, int pril = 1);
    /// see ConicBundle::CBout
    void  set_cbout(const CBout* cb, int incr = -1);

    /// write the problem description to out so that it can be read again by read_problem_data()
    std::ostream& print_problem_data(std::ostream& out) const;

    /// clear() and read the problem from in in the format written by print_problem_data()
    std::istream& read_problem_data(std::istream& in);

    /// undocumented highly volatile variant for external testing 
    std::ostream& print_problem_data_to_mfile(std::ostream& out, CH_Matrix_Classes::Integer blocknr) const;
    //std::istream& read_problem_data_from_mfile(std::istream& in);


  };


  //@}

}

#endif

