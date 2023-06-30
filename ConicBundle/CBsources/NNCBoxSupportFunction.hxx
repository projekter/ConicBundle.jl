/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/NNCBoxSupportFunction.hxx
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



#ifndef CONICBUNDLE_NNCBOXSUPPORTFUNCTION_HXX
#define CONICBUNDLE_NNCBOXSUPPORTFUNCTION_HXX

/**  @file NNCBoxSupportFunction.hxx
    @brief Header declaring the classes ConicBundle::NNCBoxSupportFunction, ConicBundle::NNCBoxSupportMinorantExtender (part of an implementation of ConicBundle::FunctionModel)
    @version 1.0
    @date 2017-09-30
    @author Christoph Helmberg
*/

//------------------------------------------------------------

#include <map>
#include "MatrixCBSolver.hxx"
#include "NNCBoxSupportModification.hxx"

//------------------------------------------------------------

namespace ConicBundle {

/**@defgroup implemented_matrix_function_oracle implemention of MatrixFunctionOracle (NNCBoxSupportFunction)

   @brief NNCBoxSupportFunction is an example implementation of
   ConicBundle::MatrixFunctionOracle for the minimization of the
   support function over a box (it is, however, recommended to use a
   BoxOracle instead, which should be more efficient). It may be used for
   Lagrangian relaxation of linear programs over boxes or, using
   FunctionTask::AdaptivePenaltyFunction, over unbounded domains with
   finite optima

   NOTE: For support functions over boxes the BoxOracle with its
   specialized BoxModel should be much better suited and more
   efficient. This is just an additional possible option and serves as
   a good example implementation of a MatricFunctionOracle.

   The class NNCBoxSupportFunction implements a general purpose version
   of MatrixFunctionOracle for minimizing the support function
   
   \f[f:\mathbf{R}^n\to\mathbf{R},\quad f(\tilde c)=\max\{\tilde c^\top x\colon x\in[l,u]\},\f] 

   where \f$l,u\in\mathbf{R}^n\f$ are given lower and upper bounds
   (\f$l\le u\f$ componentwise). When adding this function to the
   solver by MatrixCBsolver::add_function(), one may specify a weight
   \f$\gamma>0\f$, a ConicBundle::FunctionTask, and an
   AffineFunctionTransformation \f$F\colon
   \mathbf{R}^m\to\mathbf{R}^n, F(y)=c+Ay\f$ for given
   \f$c\in\mathbf{R}^n\f$ and \f$A\in\mathbf{R}^{n\times
   m}\f$. Depending on the ConicBundle::FunctionTask, the solver
   actually minimizes the following:

   - FunctionTask::ObjectiveFunction 

     \f[ f_{\gamma,F}(y)=\gamma f(F(y))=\gamma\max\{(c+Ay)^\top x\colon x\in[l,u]\},\f] 

     If in setting up the groundset within MatrixCBSolver the \f$y\f$
     variables are introduced with a linear cost vector
     \f$b\mathbf{R}^m\f$ and appropriate sign constraints, this
     corresponds to Lagrangian relaxation of the linear constraints in
     the linear program

     \f[ \mbox{maximize } c^\top x 
         \mbox{ subject to } A^\top x \begin{array}{c}\le\\=\\\ge\end{array} b, 
         x\in[\gamma l,\gamma u]
     \f]
     
   - FunctionTask::ConstantPenaltyFunction 

     \f[ f^+_{\gamma,F}(y)=\gamma\max\{0,f(F(y))\},\f]

     or, as above in the Lagrangean relaxation setting

     \f[ \mbox{maximize } c^\top x 
         \mbox{ subject to } A^\top x \begin{array}{c}\le\\=\\\ge\end{array} b, 
         x\in[\delta l,\delta u]\mbox{ with }0\le \delta\le \gamma 
     \f]
     
   - FunctionTask::AdaptivePenaltyFunction 

     \f[ f^+_{\infty,F}(y)=\gamma\max\{0,f(F(y))\}\mbox{ with }\gamma\to \infty.\f]

     or, as above in the Lagrangean relaxation setting

     \f[ \mbox{maximize } c^\top x 
         \mbox{ subject to } A^\top x \begin{array}{c}\le\\=\\\ge\end{array} b, 
         x\in[\delta l,\delta u]\mbox{ with }0\le \delta\le \gamma  \mbox{ and }\gamma\to \infty
     \f]
     
   The standard use of a NNCBoxSupportFunction is to fully initialize it
   on construction with the data l and u and then to add this function
   to the MatrixCBSolver. If dynamic changes to this
   NNCBoxSupportFunction are required afterwards (only adding and
   deleting box coordinates is allowed), use the class
   ConicBundle::NNCBoxSupportModification within the corresponding
   problem modification routines of the MatrixCBSolver interface.

   The minorant generated for the input cost vector \f$\tilde c\f$ is
   a maximizing vertex of the box, so it is itself the basic primal
   data and no extram primal data is needed.

   For facilitating input and output, NNCBoxSupportFunction offers
   - NNCBoxSupportFunction::print_problem_data() that outputs the full function description
     so that it can be read again by read_problem_data
   - NNCBoxSupportFunction::read_problem_data() reads the problem data as output by print_problem_data()
   - NNCBoxSupportFunction::set_out() and NNCBoxSupportFunction::set_cbout() work as described in ConicBundle::CBout 
 
*/

//@{

  class NNCBoxSupportFunction;

  /** @brief Implementation of MinorantExtender for NNCBoxSupportFunction

      This object will be returned as an object on the heap 
      by NNCBoxSupportFunction::apply_modification
      and will be deleted by ConicBundle after use.
      Its purpose is to fill up further coordinates of the minorant
      (if this is possible).

  */

  class NNCBoxSupportMinorantExtender: public CBout, public MinorantExtender
  {
  private:
    /// the oracle this MinorantExtender was generated by, needed for retrieving problem data
    NNCBoxSupportFunction* fun;

  public:
    ///the NNCBoxSupportFunction pointed to has to be valid for the livetime of this object
    NNCBoxSupportMinorantExtender(NNCBoxSupportFunction* fun);
			
    ///
    ~NNCBoxSupportMinorantExtender();

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
    int extend(Minorant& minorant,int n_coords,const int* indices);

  };



  /**@brief general purpose implementation of MatrixFunctionOracle as explained in \ref implemented_matrix_function_oracle (see, however, the likely better choice of a BoxOracle)
   */

  class NNCBoxSupportFunction: public MatrixFunctionOracle, public CBout
  {
  private:
    CH_Matrix_Classes::Matrix lb; ///< column vector of lower bounds
    CH_Matrix_Classes::Matrix ub; ///< column vector of upper bounds
    
    /// applies the NNCBoxSupportModfication mod to the current function
    int apply_modification(const NNCBoxSupportModification& mod);
    
  public:
    /**@name Initialization */
    //@{
    /// intialize with lower bounds vector lb and upper bounds vector ub (and output options), both must be column vectors of the same length and lb<=ub componentwise, length 0 results in objective value 0
    NNCBoxSupportFunction(const CH_Matrix_Classes::Matrix& lb,
		       const CH_Matrix_Classes::Matrix& ub,
		       const CBout* cb=0,int incr=-1);    

    ///
    ~NNCBoxSupportFunction(){}

    //@}

    //----------- Oracle Implementation of MatrixFunctionOracle ----------
    /**@name Implementations of MatrixFunctionOracle routines */
    //@{

    /// see MatrixFunctionOracle::evaluate()
    int evaluate(const  CH_Matrix_Classes::Matrix& current_point,
		 CH_Matrix_Classes::Real relprec,
		 CH_Matrix_Classes::Real& objective_value,
		 std::vector<Minorant*>&  minorants,
		 PrimalExtender*& primal_extender);


    /** @brief see MatrixFunctionOracle::apply_modification() for the general use, here oracle_modification has a special role if it can be cast to an AMFModification

	if oracle_modification cannot be cast to an NNCBoxSupportModification it is assumed that all append modifications amount to have already been carried out on *this seperately before this routine is called. In particular, it is only checked whether the new dimension matches the one given by oracle_modification, the old dimension is ignored. If this does not hold, the routine stops with an error. Otherwise it checks the other stuff as if a suitable NNCBoxSupportModification has just been executed.
     */
    virtual 
    int 
    apply_modification
    (
     const OracleModification& oracle_modification ,
     const CH_Matrix_Classes::Matrix* new_center,
     const CH_Matrix_Classes::Matrix* old_center,
     bool& discard_objective_in_center,
     bool& discard_model, 
     bool& discard_aggregates,
     MinorantExtender*& minorant_extender
    );

    
    /// see MatrixFunctionOracle::check_correctness()   (true only needed for debugging)
     virtual
    bool
    check_correctness() const
    {return true;}
    //@}
    

    //------------------  routines for querying the problem ----------------
    /**@name routines for querying data of the problem */
    //@{

    /// returns the column vector of lower bounds
    const CH_Matrix_Classes::Matrix& get_lower_bounds() { return lb;}
    
    /// retunrs the column vector of upper bounds
    const CH_Matrix_Classes::Matrix& get_upper_bounds() {return ub;}
            
    //@}


   //------------------  routines for IO ----------------
    /**@name routines for supporting input and output */
    //@{

    /// see ConicBundle::CBout
    void  set_out(std::ostream* o=0,int pril=1);
    /// see ConicBundle::CBout
    void  set_cbout(const CBout* cb,int incr=-1);

    /// write the problem description to out so that it can be read again by read_problem_data()
    std::ostream& print_problem_data(std::ostream& out) const;

    /// clear() and read the problem from in in the format written by print_problem_data()
    std::istream& read_problem_data(std::istream& in);

    /// undocumented highly volatile variant for external testing 
    std::ostream& print_problem_data_to_mfile(std::ostream& out, CH_Matrix_Classes::Integer blocknr) const;
    //std::istream& read_problem_data_from_mfile(std::istream& in);


  };


}

  //@}

#endif

