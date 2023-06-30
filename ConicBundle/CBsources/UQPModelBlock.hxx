/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/UQPModelBlock.hxx
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



#ifndef CONICBUNDLE_UQPBLOCK_HXX
#define CONICBUNDLE_UQPBLOCK_HXX

/**  @file UQPModelBlock.hxx
    @brief Header declaring the classes ConicBundle::UQPModelBlock, ConicBundle::UQPModelPointer
    @version 1.0
    @date 2020-03-22
    @author Christoph Helmberg
*/

#include "symmat.hxx"
#include "UQPModelBlockObject.hxx"
#include "QPModelDataObject.hxx"

namespace ConicBundle {

  /** @ingroup UnconstrainedQPSolver
  */

  //@{

  /** @brief combines and provides basic functionalities of QPModelDataObject and UQPModelBlockObject, but is still abstract

      This class serves as a base class for actual implementations of the
      models and provides some basic variables and functionalities required
      for uniform use by BundleModel and QPSolver.

      In particular it provides storage for and access to the bundle
      (and possibly an additional constant minorant) of the underlying
      cutting model(s). It also ensures consistent handling of each
      AffineFunctionTransformation of this data. It might be worth
      to consider to collect the transformations and not to execute
      them at once, but direct execution of these transformations
      is the current approach.

      Actual implementations are UQPSumModelBlock (for managing the sum
      of several model blocks) and UQPConeModelBlock.
   */

  class UQPModelBlock : public virtual QPModelDataObject, public UQPModelBlockObject {
  protected:

    ///constant offset minorant (fixed affine function to be added to the model)
    std::vector<MinorantPointer> constant_minorant;
    ///the minorants forming the cutting model; how to combine them is described in derived classes
    std::vector<MinorantBundle> bundle;

  public:

    /// reset to "empty/no" model information
    void clear() {
      constant_minorant.clear(); bundle.clear();
    }

    /// default constructor
    UQPModelBlock(CBout* cb = 0, int cbinc = -1) :QPModelDataObject(cb, cbinc), UQPModelBlockObject() {
    }

    /// vritual destructor
    virtual ~UQPModelBlock();

    ///gives reading access to a constant offset minorant
    virtual const MinorantPointer& get_constant_minorant() const {
      return constant_minorant.back();
    }

    ///gives reading access to the bundle minorants of the cutting model
    virtual const MinorantBundle& get_bundle() const {
      return bundle.back();
    }

    ///gives access to a constant offset minorant
    virtual MinorantPointer& get_constant_minorant() {
      return constant_minorant.back();
    }

    ///gives access to the bundle minorants of the cutting model
    virtual MinorantBundle& get_bundle() {
      return bundle.back();
    }

    ///applies the AffineFunctionTransformation to constant_minorant and bundle, where (if given) only the global_indices of the transformed subgradients are required which need the local_indices only. If precomputed is given, it may contain some or contains afterwards a map from original minorant to transformed minorant; retunrs 0 on success
    virtual int push_aft(const AffineFunctionTransformation* inaft,
      const CH_Matrix_Classes::Indexmatrix* global_indices,
      const CH_Matrix_Classes::Indexmatrix* local_indices,
      std::map<MinorantPointer, MinorantPointer>* precomputed = 0);

    /// undo the last push_aft
    virtual int pop_aft();

    //virtual CH_Matrix_Classes::Integer xdim() const =0; //dimension of externally visible primal variables
    //virtual CH_Matrix_Classes::Integer ydim() const =0; //dimension of externally visible dual variables

    //virtual int set_qp_xstart(CH_Matrix_Classes::Integer x_start_index)=0;
    //the indices of the local variables correspond to the indices 
    //of the qp variables x and z starting with this index
    //returns 0 on success, 1 on failure

    //virtual int set_qp_ystart(CH_Matrix_Classes::Integer y_start_index)=0;
    //the indices of the local variables correspond to the indices 
    //of the qp variables y starting with this index
    //returns 0 on success, 1 on failure

    //virtual int starting_x(CH_Matrix_Classes::Matrix& qp_x)=0;
    //generate a strictly feasible primal starting point
    //store it in the qpx_range of x
    //returns 0 on success, 1 on failure

    // virtual int starting_y(CH_Matrix_Classes::Matrix& qp_y,
    // 			 const CH_Matrix_Classes::Matrix& qp_Qx,
    // 			 const CH_Matrix_Classes::Matrix& qp_c)=0;
    //generate a strictly feasible dual starting point
    //store it in the qpy_range of y
    //x is fixed already by a previous call to starting_x and Qx=Q*x
    //returns 0 on success, 1 on failure

    //virtual CH_Matrix_Classes::Real get_local_primalcost() const=0;
    //returns the current local primal cost contribution <d,s>

    //virtual CH_Matrix_Classes::Real get_local_dualcost() const=0;
    //returns the current local dual cost contribution

    //virtual int get_Ab(CH_Matrix_Classes::Matrix& qp_A,CH_Matrix_Classes::Matrix &qp_b) const=0;
    //store the local coefficients of matrices A and b in the positions
    //corresponding to qpy_range (rows) and qpx_range (columns) 
    //returns 0 on success, 1 on failure

    // virtual int restart_x(CH_Matrix_Classes::Matrix& qp_x,
    // 			const CH_Matrix_Classes::Matrix& qp_c,
    // 			const CH_Matrix_Classes::Matrix& qp_dc)=0;
    //it is assumed that the problem was solved already once and is now
    //resolved for a new linear cost term qp_c that resulted from the old
    //one by adding qp_dc.
    //on input qp_x holds the old optimal solution and on output
    //the coorespoind qpx_range should be replaced by a reasonable 
    //strictly feasible solution for x suitable for restarting
    //(see also restart_yz)
    //returns 0 on success, 1 on failure

    // virtual int restart_y(CH_Matrix_Classes::Matrix& qp_y,
    // 			const CH_Matrix_Classes::Matrix& qp_Qx,
    // 			const CH_Matrix_Classes::Matrix& qp_c,
    // 			const CH_Matrix_Classes::Matrix& qp_dc)=0;
    //this is called after restart_x (see there)
    //on input qp_y and qp_z hold the old optimal solution and on output
    //the coorespoind qpy/qpx_range should be replaced by a reasonable 
    //strictly feasible solution for y/z suitable for restarting
    //returns 0 on success, 1 on failure

    //virtual int add_xinv_kron_z(CH_Matrix_Classes::Symmatrix& barQ)=0;
    //add the system term corresponding to (xinv kron z)
    //(that arises from solving the linearized perturbed complementarity system
    // x*z =0 or =mu*I for dx in the preferred search direction)
    //to the diagonal block corresponding to qpx_range x qpx_range

    //virtual int add_local_sys(CH_Matrix_Classes::Symmatrix& sysdy,CH_Matrix_Classes::Matrix& rhs)=0;
    //on input: sysdy= A*barQ^{-1}*A^T    (barQ as returned in add_xinv_kron_z)
    //          rhs= A*barQ^{-1}*(c-Q*x-A^T*y)-(b-A*x)
    //if the block uses additional internal variables 
    //(like an additional term + B*s with s>=0 in the primal feasibility constr)
    //then the corresponding block terms have now to be added to sysdy and rhs,
    //eg,
    //   sysdy +=  B*(t^{-1} kron s)*B^T     (if t is the dual variable to s)
    //   rhs   +=  B*s - B*(t^{-1} kron s)*B^T*y

    // virtual int suggest_mu(CH_Matrix_Classes::Real& ip_xz,
    // 			 CH_Matrix_Classes::Integer& mu_dim,
    // 			 CH_Matrix_Classes::Real& sigma,
    // 			 const CH_Matrix_Classes::Matrix& qp_dx,
    // 			 const CH_Matrix_Classes::Matrix& qp_dy,
    // 			 const CH_Matrix_Classes::Matrix& rhs_residual)=0;
    //dx, dy is the predictor direction giving rise to the rhs_residual -(c-At(y+dy)-Q(x+dx)).
    //Compute the direction dz and local step and based on the predictor 
    //(x+dx,y+dy,z+dz) suggest a value for mu by specifying the
    //inner product of the dual cone variables ip_xz=ip(x,z)+ip(s,t),
    //the dimension of the conic variable space mu_dim= cone_x.dim+cone_s.dim
    //a value for the factor on mu to obtain the new target

    // virtual int get_corr(CH_Matrix_Classes::Matrix& xcorr,
    // 		       CH_Matrix_Classes::Matrix& rhs,
    // 		       CH_Matrix_Classes::Real mu)=0;
    //on input (w.r.t. corresponding positions)
    //    xcorr = 0
    //    rhs as on output of add_local_sys
    //on output the corresponding positions of xcorr should hold the corrector
    //term of the search direction, eg,  xcorr = mu*x^{-1} - x^{-1}*dx*dz,
    //and if the block holds additional local variables as in add_local_sys then
    //   rhs += B*(mu * t^{-1}- t^{-1}*dt*ds)
    //has to be called after suggest_mu which computes the other directions

    // virtual int line_search(CH_Matrix_Classes::Real& alpha,
    // 			  const CH_Matrix_Classes::Matrix& qp_dx,
    //                         const CH_Matrix_Classes::Matrix& qp_dy,
    // 			  const CH_Matrix_Classes::Matrix& rhs_residual)=0;
    //dx,dy give the final step direction, alpha is on input
    //an upper bound on the step size. 
    //On output alpha has to be no larger than on input and
    //has to guarantee strict feasibility of the primal/dual step on
    //the local variables.
    //The block has to compute the step direction dz as well as
    //for additional internal variables now and to choose alpha so
    // that strict feasibility is guaranteed for the internal
    // variables as well

    // virtual int set_point(const CH_Matrix_Classes::Matrix& qp_x,
    // 			const CH_Matrix_Classes::Matrix& qp_y,
    // 			CH_Matrix_Classes::Real alpha)=0;
    //x,y,z is the new point and has to be stored
    //alpha is the step size used in the step, it is passed so that 
    //the block can take the same step for internal variables if needed.


    //---------------- for debugging purposes

    //virtual CH_Matrix_Classes::Matrix& add_Bs(CH_Matrix_Classes::Matrix& qp_vec) const=0;
    //add the local product of matrices B and s in the positions
    //corresponding to qpy_range (rows) and return qp_vec
    //returns 0 on success, 1 on failure

    //virtual CH_Matrix_Classes::Matrix& subtract_z(CH_Matrix_Classes::Matrix& dual_residual,bool with_step=false) const=0;
    //add the contributions of the dual slacks and return dual_residual
    //returns 0 on success, 1 on failure

  };


  /** @brief Interface in BundelSolver for generating the correct type of blocks for UQPSolver and for setting the final block in the solver

     The objects generated here are the implementations
     QPSumModelBlock and QPConeModelBlock of the base class QPModelBlock.
     The ownership of the generated objects is passed over to the
     calling routine. Also the pointer to the initial block stored
     here only gives access to the block, but the object is not
     owned and may not be deleted here.

*/

  class UQPModelPointer : public virtual QPModelDataPointer {
  protected:
    UQPModelBlock* model_block; ///< stores a pointer to the current starting block giving access to the cutting model(s) [it does not own or delete this object]

  public:

    /// default constructor
    UQPModelPointer(CBout* cb = 0, int cbinc = -1) :QPModelDataPointer(cb, cbinc), model_block(0) {
    }

    /// virtual destructor
    virtual ~UQPModelPointer();

    /// set the pointer to NULL
    void clear_model_data_ptr() {
      model_block = 0;
    }

    ///store the pointer to the object if it matches the required type for the QP solver, otherwise return a nonzero value as error; this is used in the models to return the local qp model data
    int set_model_data(QPModelDataObject* inbp) {
      model_block = dynamic_cast<UQPModelBlock*>(inbp); return (model_block == 0);
    }

    ///returns a new QPSumModelDataObject, that has to be deleted by the caller. The argument is optional and allows to potentially generate different blocks for different derived BundleModel objects; this is used in SumModel to collect the models of the various oracles that are summed over 
    QPSumModelDataObject* generate_summodel_data(BundleModel* bmp = 0);

    ///returns a new QPConeModelDataObject suitable for the default conic BundleModel implementations; it has to be deleted by the caller. The argument is optional and allows to potentially generate specialized objects for special BundleModel objects 
    QPConeModelDataObject* generate_conemodel_data(BundleModel* bmp = 0);

    /// returns the pointer value
    QPModelDataObject* get_model_data_ptr() const {
      return model_block;
    }
  };

  //@}

}

#endif

