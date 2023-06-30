/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/UQPSumModelBlock.hxx
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



#ifndef CONICBUNDLE_UQPSUMMODELBLOCK_HXX
#define CONICBUNDLE_UQPSUMMODELBLOCK_HXX

/**  @file UQPSumModelBlock.hxx
    @brief Header declaring the class ConicBundle::UQPSumModelBlock
    @version 1.0
    @date 2020-03-22
    @author Christoph Helmberg
*/


#include "UQPModelBlock.hxx"

namespace ConicBundle {

  /** @ingroup UnconstrainedQPSolver
  */

  //@{

  /** @brief  implements a (virtual) cutting model being built of a (possibly recursive) sum of UQPModelBlock cutting model instances for UQPSolver

     Note, the bundle at the final level might not be the same as the
     one in the recursive calls, because it might have undergone an
     AffineFunctionTransformation.

     Not much is happening here besides passing on the calls
     to the various submodels with information on where to find
     its own (possibly transformed) bundle information and where
     to store the requested information in respective global objects
  */


  class UQPSumModelBlock : public QPSumModelDataObject, public UQPModelBlock {
  private:
    /// the container pointing to the UQPModelBlock instances comprised in the sum 
    std::vector<UQPModelBlock*> blocks;

    CH_Matrix_Classes::Integer xstart; ///< starting index of the x variables subsumed here within the globale vector x
    CH_Matrix_Classes::Integer xend;  ///< ending index of local x
    CH_Matrix_Classes::Integer ystart; ///< starting index of the y variables subsumed here within the globale vector y
    CH_Matrix_Classes::Integer yend; ///< ending index of local y


  public:
    /// reset to "empty/no" model 
    void clear() {
      blocks.clear(); xstart = 0; xend = 0; ystart = 0; yend = 0; UQPModelBlock::clear();
    }

    /// default constructor
    UQPSumModelBlock(CBout* cb = 0, int incr = -1) :QPSumModelDataObject(cb, incr), UQPModelBlock(cb, incr) {
      clear();
    }

    /// virtual destructor
    ~UQPSumModelBlock();

    /// add another model at the end of the list
    int append(QPModelDataObject* inblock);

    //-----------  UQPModelBlock routines

    /// sum of all xdim of the sublocks
    CH_Matrix_Classes::Integer xdim() const;

    /// sum of all ydim of the sublocks
    CH_Matrix_Classes::Integer ydim() const;

    /// set the starting index of x also for all subblocks
    int set_qp_xstart(CH_Matrix_Classes::Integer x_start_index);

    /// set the starting index of y also for all subblocks
    int set_qp_ystart(CH_Matrix_Classes::Integer y_start_index);

    /// get the starting x of all subblocks
    int starting_x(CH_Matrix_Classes::Matrix& qp_x);

    /// get the starting y information of all subblocks
    int starting_y(CH_Matrix_Classes::Matrix& qp_y,
      const CH_Matrix_Classes::Matrix& qp_Qx,
      const CH_Matrix_Classes::Matrix& qp_c);

    /// get joint primalcost of all subblocks
    CH_Matrix_Classes::Real get_local_primalcost() const;

    /// get joint dualcost of all subblocks
    CH_Matrix_Classes::Real get_local_dualcost() const;

    /// get the A matrix of all subblocks and store it consistently
    int get_Ab(CH_Matrix_Classes::Matrix& qp_A, CH_Matrix_Classes::Matrix& qp_b) const;

    /// get a good restarting x of all subblocks for this change in costs
    int restart_x(CH_Matrix_Classes::Matrix& qp_x, const CH_Matrix_Classes::Matrix& qp_c, const CH_Matrix_Classes::Matrix& qp_dc);


    /// get a good restarting y  of all subblocks for this change in costs
    int restart_y(CH_Matrix_Classes::Matrix& qp_y,
      const CH_Matrix_Classes::Matrix& qp_Qx,
      const CH_Matrix_Classes::Matrix& qp_c,
      const CH_Matrix_Classes::Matrix& qp_dc);

    /// add this for all subblocks
    int add_xinv_kron_z(CH_Matrix_Classes::Symmatrix& barQ);


    /// add this for all subblocks
    int add_local_sys(CH_Matrix_Classes::Symmatrix& sysdy, CH_Matrix_Classes::Matrix& rhs);

    /// get this from all subblocks
    int suggest_mu(CH_Matrix_Classes::Real& ip_xz,
      CH_Matrix_Classes::Integer& mu_dim,
      CH_Matrix_Classes::Real& sigma,
      const CH_Matrix_Classes::Matrix& qp_dx,
      const CH_Matrix_Classes::Matrix& qp_dy,
      const CH_Matrix_Classes::Matrix& rhs_residual);


    /// get this from all subblocks
    int get_corr(CH_Matrix_Classes::Matrix& xcorr,
      CH_Matrix_Classes::Matrix& rhs,
      CH_Matrix_Classes::Real mu);


    /// get/do this from/for all subblocks
    int line_search(CH_Matrix_Classes::Real& alpha,
      const CH_Matrix_Classes::Matrix& qp_dx,
      const CH_Matrix_Classes::Matrix& qp_dy,
      const CH_Matrix_Classes::Matrix& rhs_residual);


    /// do this for all subblocks
    int set_point(const CH_Matrix_Classes::Matrix& qp_x,
      const CH_Matrix_Classes::Matrix& qp_y,
      CH_Matrix_Classes::Real alpha);

    /// do this for all subblocks
    int add_modelx_aggregate(CH_Matrix_Classes::Real& offset,
      CH_Matrix_Classes::Matrix& gradient);

    /// do this for all subblocks
    void set_out(std::ostream* o = 0, int pril = 1);

    /// do this for all subblocks
    void set_cbout(const CBout* cb, int incr = -1);

    //---------------- for debugging purposes

    /// do this for all subblocks
    CH_Matrix_Classes::Matrix& add_Bs(CH_Matrix_Classes::Matrix& qp_vec) const;

    /// do this for all subblocks
    CH_Matrix_Classes::Matrix& subtract_z(CH_Matrix_Classes::Matrix& dual_residual, bool with_step = false) const;

  };

  //@}

}

#endif

