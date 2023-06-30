/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/UQPConeModelBlock.hxx
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



#ifndef CONICBUNDLE_UQPCONEMODELBLOCK_HXX
#define CONICBUNDLE_UQPCONEMODELBLOCK_HXX

/**  @file UQPConeModelBlock.hxx
    @brief Header declaring the class ConicBundle::UQPConeModelBlock
    @version 1.0
    @date 2020-03-25
    @author Christoph Helmberg
*/

#include <vector>
#include "UQPModelBlock.hxx"


namespace ConicBundle {

  /** @ingroup UnconstrainedQPSolver
   */

   //@{

   /** @brief  implements a UQPModelBlock for conic cutting models in UQPSolver

      note, the bundle at the final level might not be the same as the
      one entered in the intial call, because it might have undergone
      affine function transformations.

      The starting point of the model is the call to init() as described
      in QPConeModelDataObject::init(), which sets
      the bundle (subgradient or minorant information)
      with the corresponding description of the model variables and
      constraints specifying on how the bundle inforamtion is to be
      combined. The description possibilities here allow to formulate
      the default ConicBundle models NNCModel, SOCModel, PSCModel and
      BoxModel. The sequence of the cones with their dimensions must
      match the sequence given in the bundle. All cones are combined
      by at most one single row trace constraint, which might even be
      missing in the case of a full box model. The right hand side
      of the trace constraint serves as function factor in the equality
      case or, in the inequalitiy case, as penalty parameter.



   */


  class UQPConeModelBlock : public UQPModelBlock, public QPConeModelDataObject {
  private:
    ///joint local primal variables as a vector, x=(vecx',svec(X1)',svec(X2)',...,svec(Xk)',boxx)'
    CH_Matrix_Classes::Matrix x;

    ///joint local dual slack variables as a vector, z=(vecz',svec(Z1)',svec(Z2)',...,svec(Zk)',boxlz)'
    CH_Matrix_Classes::Matrix z;

    ///here, A is only a row vector (the trace), A=(ones(vecz)', svec(I1)',....,svec(Ik)')
    CH_Matrix_Classes::Matrix A;

    /// the trace right hand side
    CH_Matrix_Classes::Real b;

    ///if !=0, then Ax <= b
    int less_or_equal;
    CH_Matrix_Classes::Real s;      ///< primal slack variable if inequality
    CH_Matrix_Classes::Real scost;  ///< cost of primal slack if ineqquality
    CH_Matrix_Classes::Real y;      ///< dual variable for A

    CH_Matrix_Classes::Integer nnc_dim;    ///< dimension of leading linear part

    CH_Matrix_Classes::Indexmatrix soc_dim; ///< dimensions of second order cone variables
    CH_Matrix_Classes::Indexmatrix soc_start; ///< first element of i-th SOC variable in x

    std::vector<CH_Matrix_Classes::Symmatrix> Xp; ///< local primal psd variables X 
    std::vector<CH_Matrix_Classes::Symmatrix> Zp; ///< local dual slack psd variables Z
    std::vector<CH_Matrix_Classes::Symmatrix> Xinv; ///< inverses of local primal psd variables X 
    //std::vector<CH_Matrix_Classes::Symmatrix> Xchol; //inverses of local primal psd variables X 
    //std::vector<CH_Matrix_Classes::Symmatrix> Zchol; //inverses of local primal psd variables Z

    CH_Matrix_Classes::Matrix lb;    ///< box lower bounds
    CH_Matrix_Classes::Matrix ub;    ///< box upper bounds
    bool box_scaling; ///< true if box takes part in scaling
    bool box_scaleub; ///< if true, the scaling is upper bounded by b


    //CH_Matrix_Classes::Matrix lz;    // box dual variable for lower bounds is in z
    CH_Matrix_Classes::Matrix uz;    ///< box dual variable for upper bounds(if box_scaleub)
    //CH_Matrix_Classes::Real box_lz;  //dual variable for lower bound on box_s is in z
    CH_Matrix_Classes::Real box_uz;  ///< dual variable for upper bound on box_s


    CH_Matrix_Classes::Integer box_dim;  ///< the size of lb + scaling (if so)

    CH_Matrix_Classes::Integer box_start;  ///< the index of box_x in qp_x


    CH_Matrix_Classes::Integer qp_xstart; ///< the first index of the local x in qp_x
    CH_Matrix_Classes::Integer qp_ystart; ///< the index of the local y in qp_y

    CH_Matrix_Classes::Integer mu_dim;    ///< = vecx.dim + sum of the matrix orders [ + 1]
    CH_Matrix_Classes::Real restart_factor; ///< after a cost update, how far to go back into the interior


    //---- copies of old variables for analyzing activity
    CH_Matrix_Classes::Real current_mu;   ///< current value of the barrier parameter
    CH_Matrix_Classes::Real old_mu;  ///< previous value of the barrier parameter
    CH_Matrix_Classes::Matrix old_x;  ///< previous value of x                 
    CH_Matrix_Classes::Matrix old_z;  ///< previous value of z
    CH_Matrix_Classes::Matrix old_uz;   ///< previous value of uz
    CH_Matrix_Classes::Real old_box_uz;  ///< previous value of box_uz 
    CH_Matrix_Classes::Real old_y;  ///< previous value of y                   
    CH_Matrix_Classes::Real old_s;  ///< previous value of s
    std::vector<CH_Matrix_Classes::Symmatrix> old_Xp; ///< previous value of Xp
    std::vector<CH_Matrix_Classes::Symmatrix> old_Zp; ///< previous value of Zp
    CH_Matrix_Classes::Real last_alpha; ///< most recent value of step size alpha

    //---- for forming or testing the corrector the old step is copied here
    CH_Matrix_Classes::Matrix zcorr;   ///< value used in corrector of z
    CH_Matrix_Classes::Matrix uzcorr;  ///< value used in corrector of uz
    CH_Matrix_Classes::Real box_uzcorr; ///< value used in corrector of box_uz 
    CH_Matrix_Classes::Real scorr;    ///< value used in corrector of s

    mutable CH_Matrix_Classes::Matrix dx;    ///< step for x
    mutable CH_Matrix_Classes::Real dy;      ///< step for y             
    mutable CH_Matrix_Classes::Matrix dz;    ///< step for z                
    mutable CH_Matrix_Classes::Matrix duz;   ///< step for uz
    mutable CH_Matrix_Classes::Real box_duz; ///< step for box_uz
    mutable CH_Matrix_Classes::Real ds;      ///< step for s

    //---- temporary variables
    mutable CH_Matrix_Classes::Matrix tmpvec;  ///< temporary variable to reduce reallocations
    mutable CH_Matrix_Classes::Matrix tmpsvec; ///< temporary variable to reduce reallocations
    mutable CH_Matrix_Classes::Symmatrix tmpsymmat; ///< temporary variable to reduce reallocations
    mutable CH_Matrix_Classes::Symmatrix dX; ///< temporary variable to reduce reallocations

    //the Block has to memorize and update its local x,y,z variables 
    //(and all additional ones that are not visible to the outside)
    //in particular it is assumed to keep a copy of the current variables
    //after the calls to the following functions:
    //  starting_x, starting_yz, restart_x, restart_yz, set_point

    /// used for second order cone computations
    CH_Matrix_Classes::Matrix mult_Arw(const CH_Matrix_Classes::Matrix& x, const CH_Matrix_Classes::Matrix& v) const;

    /// used for second order cone computations
    CH_Matrix_Classes::Matrix mult_Arwinv(const CH_Matrix_Classes::Matrix& x, const CH_Matrix_Classes::Matrix& v) const;

    /// used for second order cone computations
    CH_Matrix_Classes::Matrix mult_G(const CH_Matrix_Classes::Matrix& x, CH_Matrix_Classes::Real gamma, const CH_Matrix_Classes::Matrix& v) const;

    /// used for second order cone computations
    CH_Matrix_Classes::Matrix mult_Ginv(const CH_Matrix_Classes::Matrix& x, CH_Matrix_Classes::Real gamma, const CH_Matrix_Classes::Matrix& v) const;



  public:
    /// default constructor
    UQPConeModelBlock(CBout* cb = 0, int cbinc = -1) :UQPModelBlock(cb, cbinc) {
    }

    /// destructor
    ~UQPConeModelBlock() {
    }

    /// reinitialize to "empty/no" model
    void clear();


    /// sets up the model with bundle information and how to combine it, see QPConeModelDataObject::init() for a detailed description
    int init(const MinorantPointer& constant_minorant,
      const MinorantBundle& bundle,
      CH_Matrix_Classes::Integer nnc_dim,
      const CH_Matrix_Classes::Indexmatrix& soc_dim,
      const CH_Matrix_Classes::Indexmatrix& sdp_dim,
      const CH_Matrix_Classes::Matrix& box_lb,
      const CH_Matrix_Classes::Matrix& box_ub,
      CH_Matrix_Classes::Real b,
      FunctionTask ft,
      QPModelOracleDataObject* oracle_data = 0,
      bool scale_box = true);

    /// change the right hand side of the trace constraint to b
    int adjust_trace(CH_Matrix_Classes::Real b);

    /// evaluate the left hand side of the trace constraint for modelx
    CH_Matrix_Classes::Real evaluate_trace() const;

    /// get the right hand side of the trace constraint
    CH_Matrix_Classes::Real get_trace() {
      return b;
    }

    /// get the linear part of modelx (and a guess, which of them are active, in {0.,1.})
    int get_nncx(CH_Matrix_Classes::Matrix& nncx,
      CH_Matrix_Classes::Matrix* nncx_activity = 0,
      bool cautious = false);

    /// get the SOC part of modelx (and a guess whether the entire cone is active
    int get_socx(CH_Matrix_Classes::Integer i,
      CH_Matrix_Classes::Matrix& socx,
      CH_Matrix_Classes::Real* socx_activity,
      bool cautious = false);

    /// get the PSC part of modelx (and a guess on the rank of the active part)
    int get_pscx(CH_Matrix_Classes::Integer i,
      CH_Matrix_Classes::Matrix& pscx_eigs,
      CH_Matrix_Classes::Matrix& pscx_vecs,
      CH_Matrix_Classes::Real& pscx_growthrate,
      CH_Matrix_Classes::Matrix& pscx_primalgrowth,
      CH_Matrix_Classes::Matrix& pscx_dualgrowth);

    /// get the box part of modelx (and a guess, which of the bounds are active, in {0.,1.})
    int get_boxx(CH_Matrix_Classes::Matrix& /* linx */,
      CH_Matrix_Classes::Matrix* /* linx_activity */ = 0,
      bool cautious = false);

    /// adds opB transposed times modelx (with offsets but without constant affine term) to the arguments
    int add_modelx_aggregate(CH_Matrix_Classes::Real& offset,
      CH_Matrix_Classes::Matrix& gradient);


    /// return the value of the dual variable to the trace consrat == support function value
    CH_Matrix_Classes::Real tracedual(CH_Matrix_Classes::Real* prec = 0) const;

    /// get the current dual non negative cone point (of the solution) 
    int get_nncz(CH_Matrix_Classes::Matrix& vecz) {
      vecz.init(nnc_dim, 1, z.get_store()); return 0;
    }
    /// get the current dual second order cone point to cone i (of the solution) 
    int get_socz(CH_Matrix_Classes::Integer i, CH_Matrix_Classes::Matrix& vecz) {
      vecz.init(soc_dim(i), 1, z.get_store() + soc_start(i)); return 0;
    }
    /// get the current primal positive semidefinite cone point to cone i (of the solution) 
    int get_X(CH_Matrix_Classes::Integer i, CH_Matrix_Classes::Symmatrix& X) {
      X.init(Xp[(unsigned long)(i)]); return 0;
    }
    /// get the current dual positive semidefinite cone point to cone i (of the solution) 
    int get_Z(CH_Matrix_Classes::Integer i, CH_Matrix_Classes::Symmatrix& Z) {
      Z.init(Zp[(unsigned long)(i)]); return 0;
    }
    /// get the current dual value of the trace constraint
    CH_Matrix_Classes::Real get_y(void) {
      return y;
    }
    /// get the current slack value of the trace constraint
    CH_Matrix_Classes::Real get_s(void) {
      if (less_or_equal) return s; return 0.;
    }
    /// get the current barrier parameter
    CH_Matrix_Classes::Real get_mu(void) {
      return current_mu;
    }


    /// get the previous primal non negative cone point (of the solution) 
    int get_old_nncx(CH_Matrix_Classes::Matrix& vecx) {
      vecx.init(nnc_dim, 1, old_x.get_store()); return 0;
    }
    /// get the previous dual non negative cone point (of the solution) 
    int get_old_nncz(CH_Matrix_Classes::Matrix& vecz) {
      vecz.init(nnc_dim, 1, old_z.get_store()); return 0;
    }
    /// get the previous primal second order cone point to cone i (of the solution) 
    int get_old_socx(CH_Matrix_Classes::Integer i, CH_Matrix_Classes::Matrix& vecx) {
      vecx.init(soc_dim(i), 1, old_x.get_store() + soc_start(i)); return 0;
    }
    /// get the previous dual second order cone point to cone i (of the solution) 
    int get_old_socz(CH_Matrix_Classes::Integer i, CH_Matrix_Classes::Matrix& vecz) {
      vecz.init(soc_dim(i), 1, old_z.get_store() + soc_start(i)); return 0;
    }
    /// get the previous primal positive semidefinite cone point to cone i (of the solution) 
    int get_old_X(CH_Matrix_Classes::Integer i, CH_Matrix_Classes::Symmatrix& X) {
      X.init(old_Xp[(unsigned long)(i)]); return 0;
    }
    /// get the previous dual positive semidefinite cone point to cone i (of the solution) 
    int get_old_Z(CH_Matrix_Classes::Integer i, CH_Matrix_Classes::Symmatrix& Z) {
      Z.init(old_Zp[(unsigned long)(i)]); return 0;
    }
    /// get the previous dual value of the trace constraint
    CH_Matrix_Classes::Real get_old_y(void) {
      return old_y;
    }
    /// get the previous slack value of the trace constraint
    CH_Matrix_Classes::Real get_old_s(void) {
      if (less_or_equal) return old_s; return 0.;
    }
    /// get the previous barrier parameter
    CH_Matrix_Classes::Real get_old_mu(void) {
      return old_mu;
    }
    /// get the most recent step size
    CH_Matrix_Classes::Real get_last_alpha() {
      return last_alpha;
    }

    /// given the steps of the global part, compute the step of the local part
    int compute_local_directions(const CH_Matrix_Classes::Matrix& qp_dx,
      const CH_Matrix_Classes::Matrix& qp_dy,
      const CH_Matrix_Classes::Matrix& rhs_resid,
      CH_Matrix_Classes::Matrix& dz,
      CH_Matrix_Classes::Matrix& duz,
      CH_Matrix_Classes::Real& box_ds,
      CH_Matrix_Classes::Real& ds);

    /// perform a line search for the given direction and return a feasible step length
    int inner_line_search(CH_Matrix_Classes::Real& alpha,
      const CH_Matrix_Classes::Matrix& qp_dx,
      const CH_Matrix_Classes::Matrix& qp_dy,
      const CH_Matrix_Classes::Matrix& dz,
      const CH_Matrix_Classes::Matrix& duz,
      CH_Matrix_Classes::Real box_ds,
      CH_Matrix_Classes::Real ds);

    //-----------  QP_Block routines (see there)

    /// dimension of externally visible primal variables
    CH_Matrix_Classes::Integer xdim() const {
      return x.rowdim();
    }

    ///dimension of externally visible dual variables
    CH_Matrix_Classes::Integer ydim() const {
      return A.rowdim();
    }

    /// the indices of the local variables correspond to the indices of the qp variables x and z starting with this index; returns 0 on success, 1 on failure
    int set_qp_xstart(CH_Matrix_Classes::Integer x_start_index) {
      qp_xstart = x_start_index; return 0;
    }

    /// the indices of the local variables correspond to the indices of the qp variables y starting with this index; returns 0 on success, 1 on failure
    int set_qp_ystart(CH_Matrix_Classes::Integer y_start_index) {
      qp_ystart = y_start_index; return 0;
    }

    /// generate a strictly feasible primal starting point, store it in the qpx_range of x; returns 0 on success, 1 on failure
    int starting_x(CH_Matrix_Classes::Matrix& qp_x);

    /// generate a strictly feasible dual starting point, store it in the qpy_range of y,  x is fixed already by a previous call to starting_x and Qx=Q*x; returns 0 on success, 1 on failure
    int starting_y(CH_Matrix_Classes::Matrix& qp_y,
      const CH_Matrix_Classes::Matrix& qp_Qx,
      const CH_Matrix_Classes::Matrix& qp_c);

    /// store the local coefficients of matrices A and b in the positions corresponding to qpy_range (rows) and qpx_range (columns); returns 0 on success, 1 on failure
    int get_Ab(CH_Matrix_Classes::Matrix& qp_A, CH_Matrix_Classes::Matrix& qp_b) const;

    /** @brief it is assumed that the problem was solved already once and is now
     resolved for a new linear cost term qp_c that resulted from the old
     one by adding qp_dc.

     on input qp_x holds the old optimal solution and on output
     the coorespoind qpx_range should be replaced by a reasonable
     strictly feasible solution for x suitable for restarting
     (see also restart_yz)

     returns 0 on success, 1 on failure
   */
    int restart_x(CH_Matrix_Classes::Matrix& qp_x,
      const CH_Matrix_Classes::Matrix& qp_c,
      const CH_Matrix_Classes::Matrix& qp_dc);

    /** @brief this is called after restart_x (see there)

       on input qp_y and qp_z hold the old optimal solution and on output
       the coorespoind qpy/qpx_range should be replaced by a reasonable
       strictly feasible solution for y/z suitable for restarting

       returns 0 on success, 1 on failure
    */
    int restart_y(CH_Matrix_Classes::Matrix& qp_y,
      const CH_Matrix_Classes::Matrix& qp_Qx,
      const CH_Matrix_Classes::Matrix& qp_c,
      const CH_Matrix_Classes::Matrix& qp_dc);

    /// add the system term corresponding to (xinv kron z) (that arises from solving the linearized perturbed complementarity system x*z =0 or =mu*I for dx in the preferred search direction) to the diagonal block corresponding to qpx_range x qpx_range
    int add_xinv_kron_z(CH_Matrix_Classes::Symmatrix& barQ);

    /** @brief add the local system informatoin

      on input:
               sysdy= A*barQ^{-1}*A^T    (barQ as returned in add_xinv_kron_z)
               rhs= A*barQ^{-1}*(c-Q*x-A^T*y)-(b-A*x)

      if the block uses additional internal variables
      (like an additional term + B*s with s>=0 in the primal feasibility constr)
      then the corresponding block terms have now to be added to sysdy and rhs, eg,
         sysdy +=  B*(t^{-1} kron s)*B^T     (if t is the dual variable to s)
         rhs   +=  B*s - B*(t^{-1} kron s)*B^T*y
    */
    int add_local_sys(CH_Matrix_Classes::Symmatrix& sysdy,
      CH_Matrix_Classes::Matrix& rhs);

    /// returns the current local primal cost contribution <d,s>
    CH_Matrix_Classes::Real get_local_primalcost() const {
      return (less_or_equal ? s * scost : 0.);
    }

    /// returns the current local dual cost contribution
    CH_Matrix_Classes::Real get_local_dualcost() const;

    /** @brief supply the information for the choice of the next barrier parameter value

       dx, dy is the predictor direction giving rise to the
       rhs_residual -(c-At(y+dy)-Q(x+dx)). Compute the direction dz and
       local step and based on the predictor (x+dx,y+dy,z+dz) suggest a
       value for mu by specifying the inner product of the dual cone
       variables ip_xz=ip(x,z)+ip(s,t), the dimension of the conic
       variable space mu_dim= cone_x.dim+cone_s.dim a value for the
       factor on mu to obtain the new target
   */
    int suggest_mu(CH_Matrix_Classes::Real& ip_xz,
      CH_Matrix_Classes::Integer& mu_dim,
      CH_Matrix_Classes::Real& sigma,
      const CH_Matrix_Classes::Matrix& qp_dx,
      const CH_Matrix_Classes::Matrix& qp_dy,
      const CH_Matrix_Classes::Matrix& rhs_residual);

    /** @brief supply the information for the corrector

     on input (w.r.t. corresponding positions)
          xcorr = 0
          rhs as on output of add_local_sys

     on output the corresponding positions of xcorr should hold the corrector
     term of the search direction, eg,  xcorr = mu*x^{-1} - x^{-1}*dx*dz,
     and if the block holds additional local variables as in add_local_sys then

             rhs += B*(mu * t^{-1}- t^{-1}*dt*ds)

     has to be called after suggest_mu which computes the other directions
   */
    int get_corr(CH_Matrix_Classes::Matrix& xcorr,
      CH_Matrix_Classes::Matrix& rhs,
      CH_Matrix_Classes::Real mu);

    /** @brief perform a line search for the block variables

      dx,dy give the final step direction, alpha is on input
      an upper bound on the step size.

      On output alpha has to be no larger than on input and
      has to guarantee strict feasibility of the primal/dual step on
      the local variables.

      The block has to compute the step direction dz as well as
      for additional internal variables now and to choose alpha so
      that strict feasibility is guaranteed for the internal
      variables as well
   */
    int line_search(CH_Matrix_Classes::Real& alpha,
      const CH_Matrix_Classes::Matrix& qp_dx,
      const CH_Matrix_Classes::Matrix& qp_dy,
      const CH_Matrix_Classes::Matrix& rhs_residual);

    ///x,y,z is the new point and has to be stored, alpha is the step size used in the step, it is passed so thatthe block can take the same step for internal variables if needed.
    int set_point(const CH_Matrix_Classes::Matrix& qp_x,
      const CH_Matrix_Classes::Matrix& qp_y,
      CH_Matrix_Classes::Real alpha);




    //---------------- for debugging purposes

    ///add the local product of matrices B and s in the positions corresponding to qpy_range (rows) and return qp_vec; returns 0 on success, 1 on failure
    CH_Matrix_Classes::Matrix& add_Bs(CH_Matrix_Classes::Matrix& qp_vec) const;

    /// add the contributions of the dual slacks and return dual_residual returns 0 on success, 1 on failure
    CH_Matrix_Classes::Matrix& subtract_z(CH_Matrix_Classes::Matrix& dual_residual, bool with_step = false) const;


  };


  //@}

}

#endif

