/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/SOCIPProxBlock.hxx
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



#ifndef CONICBUNDLE_SOCIPPROXBLOCK_HXX
#define CONICBUNDLE_SOCIPPROXBLOCK_HXX

/**  @file SOCIPProxBlock.hxx
    @brief Header declaring the class ConicBundle::SOCIPProxBlock
    @version 1.0
    @date 2020-05-02
    @author Christoph Helmberg
*/

#include "SOCIPBlock.hxx"
#include "QPSolverObject.hxx"

namespace ConicBundle {

  /** @ingroup ConstrainedQPSolver
   */

   //@{

   /** @brief  interior point variables and routines specific to primal dual complementarity conditions of a second order cone with special routines for handling the prox term

       adds prox term functionality to SOCIPBlock

       In the hope to get better centering properities for the quadratic part this
       models the prox term by a second order cone approach.
   */

  class SOCIPProxBlock : public SOCIPBlock {
  private:
    CH_Matrix_Classes::Real rho;    ///< its value is >=  (rhoconst^2+<Hy,y>)/2/rhoconst
    CH_Matrix_Classes::Real drho;   ///< step for rho
    CH_Matrix_Classes::Real rhoconst; ///< the constant for rho

    CH_Matrix_Classes::Real rho_rhs;  ///< 
    CH_Matrix_Classes::Matrix sqrtHy_rhs; ///< = e_2-[0;0;-H^.5]y-[-1;1;0]rho  but without -z

    const QPSolverProxObject* Hp; ///< points to the current prox object
    CH_Matrix_Classes::Matrix sqrtD; ///< diagonal of Hp
    const CH_Matrix_Classes::Matrix* Vp; ///< low rank contribution of Hp if available
    CH_Matrix_Classes::Integer ydim;    ///< length of the diagonal
    CH_Matrix_Classes::Integer lrdim;   ///< number of columns in the low rank part

    //----- for testing
    CH_Matrix_Classes::Matrix tmp_sysrhs;

  public:
    /// reset all point information to zero for dimension dim, the rest to zero
    int clear_prox(const QPSolverProxObject*);

    /// default constructor, allows to intialize the dimension
    SOCIPProxBlock(const QPSolverProxObject* po = 0, CBout* cb = 0, int cbinc = -1);

    /// destructor
    ~SOCIPProxBlock();

    /// reset the starting point for the current y and aim for mu or increase it if too small
    virtual int reset_starting_point(const CH_Matrix_Classes::Matrix& y,
      const CH_Matrix_Classes::Matrix& sys_lhs,
      CH_Matrix_Classes::Real& mu);

    /// add <Hy,y>/2=rho-.5 to min_objective, add x_2-.5 to max_objective and compute sys_lhs+=[0,0,-H^.5]*x 
    virtual int add_prox_contrib(CH_Matrix_Classes::Real& min_objective,
      CH_Matrix_Classes::Real& max_objective,
      CH_Matrix_Classes::Matrix& sys_lhs);

    ///return the squared Euclidean norm of the violation  
    virtual CH_Matrix_Classes::Real primalviol_2normsqr();

    //  /// add dimensions of the primal-dual pairs to mudim and add the "trace" (the inner product with center) of the respective primal-dual pair products for the current step; update the min and max values of x_i*z_i
    // virtual int get_mu_info(CH_Matrix_Classes::Integer& mudim,
    // 			  CH_Matrix_Classes::Real& tr_xz,
    // 			  CH_Matrix_Classes::Real& tr_xdzpdxz,
    // 			  CH_Matrix_Classes::Real& tr_dxdz,
    // 			  CH_Matrix_Classes::Real& min_xz,
    // 			  CH_Matrix_Classes::Real& max_xz) const;

    // /// for limiting the stepsize with respect to the neighborhood this information about norms and inner products of x(.)*z-tr_xz-tr_xz/mudim(.*)1, x.()*dz+dx(.)*z-tr_xdzpdxz/mudim(.*)1, and dx(.)*dz-tr_dxdz/mudim(.)*1 is required, each block *adds* its contribution to the numbers
    // virtual int get_nbh_info(CH_Matrix_Classes::Integer mudim,
    // 			   CH_Matrix_Classes::Real tr_xz,
    // 			   CH_Matrix_Classes::Real tr_xdzpdxz,
    // 			   CH_Matrix_Classes::Real tr_dxdz,
    // 			   CH_Matrix_Classes::Real& nrmsqr_xz,
    // 			   CH_Matrix_Classes::Real& nrmsqr_xdzpdxz,
    // 			   CH_Matrix_Classes::Real& nrmsqr_dxdz,
    // 			   CH_Matrix_Classes::Real& ip_xz_xdzpdxz,
    // 			   CH_Matrix_Classes::Real& ip_xz_dxdz,
    // 			   CH_Matrix_Classes::Real& ip_dxdz_xdzpdxz) const=0;

    // /// if necessary, reduce alpha to the biggest value so that feasibility is maintained with this step size
    // virtual int linesearch(CH_Matrix_Classes::Real& alpha) const;

    /// If mu is not zero, always add the centering term for this mu as well;
    virtual int add_prox_sysrhs(CH_Matrix_Classes::Matrix& rhs,
      CH_Matrix_Classes::Real& Hfactor,
      CH_Matrix_Classes::Real rhsmu,
      CH_Matrix_Classes::Real rhscorr);

    /// given the step in the design space, determine the second order cone steps 
    virtual int compute_step(const CH_Matrix_Classes::Matrix& step);

    /// move to (x+alpha*dx, z+alpha*dz), update diff_model and possibly reduce the model size if some part is too small relative to trace_rhs
    virtual int do_step(CH_Matrix_Classes::Real alpha,
      const CH_Matrix_Classes::Matrix& y);

    /// for testing purposes, output and contribute to the system violation by computing sys_lhs+=[0,0,-H^.5]*(x+dx) 
    virtual int test_sysviol(CH_Matrix_Classes::Matrix& sys_lhs,
      const CH_Matrix_Classes::Matrix& step);
  };


  //@}

}

#endif

