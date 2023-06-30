/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/PSCPrimal.hxx
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



#ifndef CONICBUNDLE_PSCPRIMAL_HXX
#define CONICBUNDLE_PSCPRIMAL_HXX

/**  @file PSCPrimal.hxx
    @brief Header declaring the classes ConicBundle::PSCPrimal, ConicBundle::DensePSCPrimal, ConicBundle::SparsePSCPrimal, ConicBundle::GramSparsePSCPrimal, ConicBundle::BlockPSCPrimal (designed for ConicBundle::PSCOracle and ConicBundle::AffineMatrixFunction)
    @version 1.0
    @date 2017-02-04
    @author Christoph Helmberg
*/

//------------------------------------------------------------

#include <map>
#include "CBout.hxx"
#include "PSCOracle.hxx"
#include "SparseCoeffmatMatrix.hxx"

//------------------------------------------------------------

namespace ConicBundle {

  /**@ingroup implemented_psc_oracle
   */
   //@{

    /** @brief PSCPrimal is the corresponding positive semidefinite object for
        PSCOracle like PrimalMatrix for a MatrixFunctionOracle.

       It allows the system to aggregate eigenvector information obtained in
       successive evaluations of the maximum eigenvalue function in PSCOralce
       to a positive semidefinite matrix so that in the Lagrangian relaxation
       view the primal optimal semidefinite matrix can be approximated, see \ref
       abstract_psc_oracle.

       ConicBundle::PSCModel does rely on primal information. It also does not
       care about whether and in which form primal information is
       represented. This is all hidden in the Minorant information obtained via
       PSCOracle::generate_minorant(). If provided, the primal information can
       be retrieved from the model by the routines
       PSCModel::get_approximate_primal() and PSCModel::get_center_primal()
       or via the corresponding routines in the interface.
    */

  class PSCPrimal : public PrimalData {

  public:
    /// 
    virtual ~PSCPrimal() {
    }

    /// returns a newly generated identical Object
    virtual PrimalData* clone_primal_data() const = 0;

    ///intialize the primal corresponding to the (positive semidefinite) Gram matrix P*P^T 
    virtual int assign_Gram_matrix(const CH_Matrix_Classes::Matrix& P) = 0;

    /// add factor*pd to this
    virtual int aggregate_primal_data(const PrimalData& pd, double factor = 1.) = 0;

    /// add factor*P*P^T to this
    virtual int aggregate_Gram_matrix(const CH_Matrix_Classes::Matrix& P, double factor = 1.) = 0;

    /// multiply/scale *this with a nonnegative factor
    virtual int scale_primal_data(double factor) = 0;

    /// if compatible evaluate value=ip(*this,A.column[i])
    virtual int primal_ip(CH_Matrix_Classes::Real& value,
      const SparseCoeffmatMatrix& A,
      CH_Matrix_Classes::Integer column) const = 0;
  };



  /** @brief implements a general purpose dense symmetric PSCPrimal based on CH_Matrix_Classes::Symmatrix

      DensePSCPrimal is pubically derived from PSCPrimal _and_ CH_Matrix_Classes::Symmatrix, so it may be used directly like a symmetric matrix
   */

  class DensePSCPrimal : public PSCPrimal, public CH_Matrix_Classes::Symmatrix {
  private:

  public:
    ///initialize to a symmetric matrix of size 0
    DensePSCPrimal() {
    }
    /// copy constructor
    DensePSCPrimal(const DensePSCPrimal& symmat, double factor = 1.) : PSCPrimal(), CH_Matrix_Classes::Symmatrix(symmat, factor) {
    }
    /// copy constructor from a CH_Matrix_Classes::Symmatrix
    DensePSCPrimal(const CH_Matrix_Classes::Symmatrix& symmat, double factor = 1.) : PSCPrimal(), CH_Matrix_Classes::Symmatrix(symmat, factor) {
    }
    /// direct size initialization to a zero matrix
    DensePSCPrimal(CH_Matrix_Classes::Integer n) : PSCPrimal(), CH_Matrix_Classes::Symmatrix(n, 0.) {
    }
    /// assigns a symmetric matrix
    const DensePSCPrimal& operator=(const CH_Matrix_Classes::Symmatrix& symmat) {
      init(symmat); return *this;
    }

    /// returns a newly generated identical Object, see PrimalData::clone_primal_data()
    virtual PrimalData* clone_primal_data() const {
      return new DensePSCPrimal(*this);
    }


    /// assign P*P^T to this
    virtual int assign_Gram_matrix(const CH_Matrix_Classes::Matrix& P) {
      CH_Matrix_Classes::rankadd(P, *this); return 0;
    }

    /// add factor*it to this (it must also be a DensePSCPrimal) 
    virtual int aggregate_primal_data(const PrimalData& it, double factor = 1.) {
      const DensePSCPrimal* pd = dynamic_cast<const DensePSCPrimal*>(&it);
      assert(pd != 0);
      CH_Matrix_Classes::xbpeya(*this, *pd, factor, 1.);
      return 0;
    }

    /// add factor*P*P^T to this
    virtual int aggregate_Gram_matrix(const CH_Matrix_Classes::Matrix& P, double factor = 1.) {
      CH_Matrix_Classes::rankadd(P, *this, factor, 1.); return 0;
    }

    /// multiply/scale *this with a nonnegative factor
    virtual int scale_primal_data(double factor) {
      assert(factor >= 0.); *this *= factor; return 0;
    }

    /// if compatible evaluate value=ip(*this,A.column[i])
    virtual int primal_ip(CH_Matrix_Classes::Real& value,
      const SparseCoeffmatMatrix& A,
      CH_Matrix_Classes::Integer column) const {
      if ((column < 0) || (column >= A.coldim()) || (A.blockdim().dim() != 1) || (A.blockdim()(0) != rowdim()))
        return 1;
      const SparseCoeffmatVector* block = A.block(0);
      if (block == 0) {
        value = 0.;
      } else {
        SparseCoeffmatVector::const_iterator it = block->find(column);
        if (it == block->end())
          value = 0.;
        else
          value = it->second->ip(*this);
      }
      return 0;
    }

  };


  /** @brief implements a sparse symmetric PSCPrimal collecting data only on a sparse, prespecified support; it is based on CH_Matrix_Classes::Sparsesym

      SparsePSCPrimal is pubically derived from PSCPrimal _and_ CH_Matrix_Classes::Sparseym, so it may be used directly like a sparse symmetric matrix
   */

  class SparsePSCPrimal : public PSCPrimal, public CH_Matrix_Classes::Sparsesym {
  public:
    /// copy constructor from a CH_Matrix_Classes::Sparssym, only the support of this matrix  will be used in all Gram operations
    SparsePSCPrimal(const CH_Matrix_Classes::Sparsesym& sps, double factor = 1.) : PSCPrimal(), CH_Matrix_Classes::Sparsesym(sps, factor) {
    }
    /// copy constructor, only the same support will be used in all Gram operations
    SparsePSCPrimal(const SparsePSCPrimal& pr, double factor = 1.) : PSCPrimal(), CH_Matrix_Classes::Sparsesym(pr, factor) {
    }
    ~SparsePSCPrimal() {
    }
    /// assigns this Sparsesym, only the same support will be used in all Gram operations
    SparsePSCPrimal& operator=(const CH_Matrix_Classes::Sparsesym& sdp) {
      init(sdp); return *this;
    }

    /// returns a newly generated identical Object
    PrimalData* clone_primal_data() const {
      return new SparsePSCPrimal(*this);
    }

    /// for each element aij in the support set aij=<P.row(i),P.row(j)> 
    int assign_Gram_matrix(const CH_Matrix_Classes::Matrix& P) {
      CH_Matrix_Classes::support_rankadd(P, *this); return 0;
    }

    /// if it is a SparseSDPRimal, add factor*it to this on the support of this
    int aggregate_primal_data(const PrimalData& it, double factor = 1.) {
      const SparsePSCPrimal* pd = dynamic_cast<const SparsePSCPrimal*>(&it);
      assert(pd != 0);
      this->support_xbpeya(*pd, factor, 1.);
      return 0;
    }

    /// add factor*P*P^T on the support to this
    int aggregate_Gram_matrix(const CH_Matrix_Classes::Matrix& P, double factor = 1.) {
      CH_Matrix_Classes::support_rankadd(P, *this, factor, 1.); return 0;
    }

    /// multiply/scale *this with a nonnegative factor
    virtual int scale_primal_data(double factor) {
      assert(factor >= 0.); *this *= factor; return 0;
    }

    /// if compatible evaluate value=ip(*this,A.column[i])
    virtual int primal_ip(CH_Matrix_Classes::Real& value,
      const SparseCoeffmatMatrix& A,
      CH_Matrix_Classes::Integer column) const {
      if ((column < 0) || (column >= A.coldim()) || (A.blockdim().dim() != 1) || (A.blockdim()(0) != rowdim()))
        return 1;
      const SparseCoeffmatVector* block = A.block(0);
      if (block == 0)
        value = 0.;
      else {
        SparseCoeffmatVector::const_iterator it = block->find(column);
        if (it == block->end())
          value = 0.;
        else {
          if (it->second->support_in(*this) == 0) {
            return 1;
          }
          value = it->second->ip(*this);
        }
      }
      return 0;
    }
  };


  /** @brief represents an PSCPrimal as the sum \f$PP^T+S\f$ of a Gram matrix and a sparse symmetric matrix \f$S\f$

      GramSparsePSCPrimal is pubically derived from PSCPrimal _and_ CH_Matrix_Classes::Sparseym, so it may be used directly like a sparse symmetric matrix, but this does then _not_ include or involve the Gram matrix part!
   */

  class GramSparsePSCPrimal : public PSCPrimal, public CH_Matrix_Classes::Sparsesym {
  protected:
    CH_Matrix_Classes::Matrix gramblock; ///< the gram matrix part is gramblock*transpose(grampblock)

  public:
    /// initialize to the given sparse symmetric matrix, the gram part is zero
    GramSparsePSCPrimal(const CH_Matrix_Classes::Sparsesym& sps, double factor = 1.) : PSCPrimal(), CH_Matrix_Classes::Sparsesym(sps, factor) {
    }
    /// copy constructor
    GramSparsePSCPrimal(const GramSparsePSCPrimal& pr, double factor = 1.) :
      PSCPrimal(),
      CH_Matrix_Classes::Sparsesym(pr, factor),
      gramblock(pr.gramblock, std::sqrt(factor)) {
      assert(factor >= 0.);
    }
    ///
    ~GramSparsePSCPrimal() {
    }
    /// assign the sparse symmetric matrix to this and set the gram part to zero
    GramSparsePSCPrimal& operator=(const CH_Matrix_Classes::Sparsesym& sdp) {
      init(sdp); gramblock.init(0, 0, 0.); return *this;
    }
    /// copy the information
    GramSparsePSCPrimal& operator=(const GramSparsePSCPrimal& sdp) {
      init(sdp); gramblock = sdp.gramblock; return *this;
    }

    /// returns the matrix \f$P\f$ giving rise to the Gram matrix \f$PP^T\f$
    const CH_Matrix_Classes::Matrix& get_grammatrix() const {
      return  gramblock;
    }

    /// returns a newly generated identical Object
    PrimalData* clone_primal_data() const {
      return new GramSparsePSCPrimal(*this);
    }

    /// Set the grammatrix part to \f$P\f$ and set all values on the sparse support to zero (but keep this support even if it is zero now!)
    int assign_Gram_matrix(const CH_Matrix_Classes::Matrix& P) {
      gramblock = P; CH_Matrix_Classes::support_rankadd(CH_Matrix_Classes::Matrix(gramblock.rowdim(), 1, 0.), *this, 0.); return 0;
    }

    /**@brief if it is a GramSparsePSCPrimal or SparseSDPRimal, add factor*it to this on only the support of this sparsematrix

       Even if it is a GramSparsePSCPrimal and it has a nontirival Gram matrix part,
       this part is only added to the sparse part of this on the support of the
       sparse part of this. No attempt is made to enlarge the Gram part.
     */
    int aggregate_primal_data(const PrimalData& it, double factor = 1.) {
      assert(factor >= 0.);
      const GramSparsePSCPrimal* pd = dynamic_cast<const GramSparsePSCPrimal*>(&it);
      if (pd != 0) {
        this->support_xbpeya(*pd, factor, 1.);
        if (pd->gramblock.dim() != 0)
          CH_Matrix_Classes::support_rankadd(pd->gramblock, *this, factor);
        return 0;
      }
      const SparsePSCPrimal* ps = dynamic_cast<const SparsePSCPrimal*>(&it);
      if (ps != 0) {
        this->support_xbpeya(*ps, factor, 1.);
        return 0;
      }
      return 1;
    }

    /**@brief add factor*P*P^T to this, collecting
       all available information only in the sparse part, even the own Gram part.

       This operation more or less converts this to a SparsePSCPrimal. The point
       is that this routine is typically only called by PSCModel when aggregating
       the information in a single aggregate matrix over several steps and it is
       pointless to try to keep any gram information in this case.
     */
    int aggregate_Gram_matrix(const CH_Matrix_Classes::Matrix& P, double factor = 1.) {
      if (gramblock.dim() != 0) {
        CH_Matrix_Classes::support_rankadd(gramblock, *this);
        gramblock.init(0, 0, 0.);
      }
      CH_Matrix_Classes::support_rankadd(P, *this, factor, 1.);
      return 0;
    }

    /// multiply/scale *this with a nonnegative factor
    virtual int scale_primal_data(double factor) {
      assert(factor >= 0.); *this *= factor; gramblock *= std::sqrt(factor); return 0;
    }

    /// if compatible evaluate value=ip(*this,A.column[i])
    virtual int primal_ip(CH_Matrix_Classes::Real& value,
      const SparseCoeffmatMatrix& A,
      CH_Matrix_Classes::Integer column) const {
      if ((column < 0) || (column >= A.coldim()) || (A.blockdim().dim() != 1) || (A.blockdim()(0) != rowdim()))
        return 1;
      const SparseCoeffmatVector* block = A.block(0);
      if (block == 0) {
        value = 0.;
      } else {
        SparseCoeffmatVector::const_iterator it = block->find(column);
        if (it == block->end())
          value = 0.;
        else {
          if (it->second->support_in(*this) == 0) {
            return 1;
          }
          value = it->second->ip(*this);
          if (gramblock.dim() > 0)
            value += it->second->gramip(gramblock);
        }
      }
      return 0;
    }
  };

  /** @brief implements a block diagonal PSCPrimal consisting of several PSCPrimal blocks

      The partition is generated on construction and at this point the various kinds
      of PSCPrimal variants are also assigned to each block by a map. For some blocks
      the assignment may be empty so that no information is aggregated for these.
      Note, the partition need not match some natural partition of the original
      coefficient data, even if it may and will in most practical applications.
   */



  class BlockPSCPrimal : public PSCPrimal {
  private:
    CH_Matrix_Classes::Indexmatrix Xdim; ///< stores the sizes of the block
    std::map<CH_Matrix_Classes::Integer, PSCPrimal*> primal; ///< holds the PSCPrimal of blcok i if there is one
  public:
    /// copy constructor
    BlockPSCPrimal(const BlockPSCPrimal& pr, double factor = 1.);
    /** in this version control over the objects pointed to in pr
  is passed to this. this will delete them upon its destruction */
    BlockPSCPrimal(const CH_Matrix_Classes::Indexmatrix& Xd, const std::map<CH_Matrix_Classes::Integer, PSCPrimal*>& pr, double factor = 1.);
    ///
    ~BlockPSCPrimal();

    /// returns a newly generated identical Object
    virtual PrimalData* clone_primal_data() const;

    /**@brief for each element aij in the support set aij=<P.row(i),P.row(j)>

       The rows of P corresponding to each block are passed on the
       with a corresponding call to the PSCPrimal of this block
    */
    virtual int assign_Gram_matrix(const CH_Matrix_Classes::Matrix& P);

    /** add factor*it to this
        (it must also be a BlockPSCPrimal with the same block partition) */
    virtual int aggregate_primal_data(const PrimalData& it, double factor = 1.);

    /**@brief add factor*P*P^T to this

       The rows of P corresponding to each block are passed on the
       with a corresponding call to the PSCPrimal of this block
    */
    virtual int aggregate_Gram_matrix(const CH_Matrix_Classes::Matrix& P, double factor = 1.);

    /// returns the number of blocks 
    CH_Matrix_Classes::Integer get_nblocks() const {
      return Xdim.dim();
    }
    /// returns the size of block i
    CH_Matrix_Classes::Integer blockdim(CH_Matrix_Classes::Integer i) const {
      return Xdim(i);
    }
    /// returns the PSCPrimal of block i if there is one, 0 otherwise
    PSCPrimal* block(CH_Matrix_Classes::Integer i) const;

    /// multiply/scale *this with a nonnegative factor
    virtual int scale_primal_data(double factor);

    /// if compatible evaluate value=ip(*this,A.column[i])
    virtual int primal_ip(CH_Matrix_Classes::Real& value,
      const SparseCoeffmatMatrix& A,
      CH_Matrix_Classes::Integer column) const;

  };

  //@}

}


#endif

