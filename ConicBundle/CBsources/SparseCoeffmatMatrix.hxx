/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/SparseCoeffmatMatrix.hxx
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



#ifndef CONICBUNDLE_SPARSECOEFFMATMATRIX_HXX
#define CONICBUNDLE_SPARSECOEFFMATMATRIX_HXX

/**  @file SparseCoeffmatMatrix.hxx
    @brief Header declaring the classes ConicBundle::SparseCoeffmatMatrix  (needed for ConicBundle::PSCAffineFunction)
    @version 1.0
    @date 2016-12-20
    @author Christoph Helmberg
*/

#include "Coeffmat.hxx"

namespace ConicBundle {

  class PSCPrimal;

/** @ingroup implemented_psc_oracle
 */
//@{ 

  /// convenient for initializing SparseCoeffmatMatrix via the sparse (block_i,column_i,Coeffmat_i), i=1,...,nz (nonzeros) format with Indexmatrix,Indexmatrix,CoeffmatVector
  typedef std::vector<CoeffmatPointer> CoeffmatVector;
  /// this is used to extract a row/block or a column from a SparseCoeffmatMatrix
  typedef std::map<CH_Matrix_Classes::Integer,CoeffmatPointer> SparseCoeffmatVector;

  /** @brief stores/organizes the CoeffmatPointer pointers to Coeffmat matrices 
      with the purpose of describing the block diagonal symmetric matrices used in PSCAffineFunction

      In an \f$k\times m\f$ SparseCoeffmatMatrix A, each row \f$i=1,\dots,k\f$ 
      describes a linear symmetric matrix function \f$\sum_{j=1}^my_iA_i\f$ 
      of matrices \f$A_i\f$ of order blockdim(i). Alternatively each 
      column \f$j=1,\dots,m\f$ of A 
      describes a block diagonal matrix consisting of k diagonal blocks 
      with block i of size blockdim(i).

      This class offers a number of routines for modifications, like
      appending further blocks or columns, deleting, reassigning indices, 
      setting single entries etc.

      Likewise it allows to extract single rows or columns in the form
      of a SparseCoeffmatVector 

      For each Coeffmat coefficient matrix the matrix only stores a 
      CoeffmatPointer pointing to it, so if the Coeffmat pointed to is 
      changed externally, this matrix class will not notice this but
      simply keep pointing to the changed object.
   */
  class SparseCoeffmatMatrix : virtual public CBout
{
public:
  /// this will be used for a (lazy) column representation of the matrix
  typedef std::map<CH_Matrix_Classes::Integer,SparseCoeffmatVector> SCMcolrep;

private:
  CH_Matrix_Classes::Indexmatrix block_dim; ///< column vector holding for each row the dimension of the block
  CH_Matrix_Classes::Integer col_dim; ///< number of columns in the SparseCoeffmatMatrix
  CH_Matrix_Classes::Indexmatrix dense_cnt; ///< number of dense matrices per block
  
  /// each row/block is likely to have at least one nonzero entry, therefore the block representation uses a vector of SparseCoeffmatVector (instead of a map) 
  typedef std::vector<SparseCoeffmatVector> SCMblockrep;
  /// this is the basic block representation, it holds for each row/block a sparse representation of the row via a SparseCoeffmatVector; this is the real work horse
  SCMblockrep blockrep;

  
  mutable SCMcolrep* colrep; ///< this column representation is only formed on demand and deleted on changes

 
  /// rebuilds the column representation from the block representation (if needed) 
  void form_colrep() const;
  
    
public:
  /// copy
  SparseCoeffmatMatrix& operator=(const SparseCoeffmatMatrix& );

  ///clears all (empty 0 times 0 matrix)  
  void clear();

  ///first calls clear() and then it sets the new values (if one of block_ind, col_ind, or coeff_vec is !=NULL, all must be !=NULL and of the same size)
  int init(const CH_Matrix_Classes::Indexmatrix& block_dim,
	   CH_Matrix_Classes::Integer col_dim,
	   const CH_Matrix_Classes::Indexmatrix* block_ind=0,
	   const CH_Matrix_Classes::Indexmatrix* col_ind=0,
	   const CoeffmatVector* coeff_vec=0);

  ///set the output and call clear()
  SparseCoeffmatMatrix(const CBout* cb=0,
		       int incr=-1):
    CBout(cb,incr),col_dim(0),blockrep(),colrep(0)
  {clear();}

  ///set the output and call clear()
  SparseCoeffmatMatrix(const SparseCoeffmatMatrix& S,const CBout* cb=0,
		       int incr=-1):
    CBout(cb,incr),col_dim(0),blockrep(),colrep(0)
  { *this=S; }
  

  ///set the output and call init() for the given sparse information
  SparseCoeffmatMatrix(const CH_Matrix_Classes::Indexmatrix& in_block_dim,
		       CH_Matrix_Classes::Integer in_col_dim,
		       const CH_Matrix_Classes::Indexmatrix* block_ind=0,
		       const CH_Matrix_Classes::Indexmatrix* col_ind=0,
		       const CoeffmatVector* coeff_vec=0,
		       const CBout* cb=0,
		       int incr=-1):
    CBout(cb,incr),col_dim(0),blockrep(),colrep(0)
  { init(in_block_dim,in_col_dim,block_ind,col_ind,coeff_vec); }

  ///calls clear() and exits
  ~SparseCoeffmatMatrix();


  /// returns a column vector with row i giving the order of block i
  const CH_Matrix_Classes::Indexmatrix& blockdim() const {return block_dim;}
  /// returns the order of block i
  CH_Matrix_Classes::Integer blockdim(CH_Matrix_Classes::Integer i) const;
  /// returns the number of columns (i.e. of blockdiagonal matrices) 
  CH_Matrix_Classes::Integer coldim() const {return col_dim;}
  /// returns the number of blocks
  CH_Matrix_Classes::Integer rowdim() const {return block_dim.rowdim();}
  /// returns the number of columns with nonzero coefficient matrices
  CH_Matrix_Classes::Integer nzcoldim() const
  {form_colrep(); return CH_Matrix_Classes::Integer(colrep->size()); }

  /// sets the CoeffmatPointer of block i in column j (blockdiagonal matrix j) (may be empty)  
  int set(CH_Matrix_Classes::Integer i,CH_Matrix_Classes::Integer j,const CoeffmatPointer& cm);

  /// sets the CoeffmatPointer of block i in column j (blockdiagonal matrix j) to point to cm (or deletes it if cm==0)  
  int set(CH_Matrix_Classes::Integer i,CH_Matrix_Classes::Integer j,Coeffmat* cm)
  {CoeffmatPointer cmp(cm); return set(i,j,cmp);}

  /// returns the CoeffmatPointer of block i in column j (blockdiagonal matrix j) (may be empty)
  CoeffmatPointer operator()(CH_Matrix_Classes::Integer i,CH_Matrix_Classes::Integer j) const;

  /// returns NULL if the block is empty and a pointer to a Vector of CoeffmatPointers otherwise
  const SparseCoeffmatVector* block(CH_Matrix_Classes::Integer i) const;
  /// returns NULL if the column is empty and a pointer to a Vector of CoeffmatPointers otherwise
  const SparseCoeffmatVector* column(CH_Matrix_Classes::Integer i) const;

  
  /** @brief append append_mat (or its submatrix given by blocks and/or cols) below
     
      If blocks or cols have negative entries, each of these represent a row of blocks of zeros whose order is the negative of the given value or a column of zeros of appropriate size (the size of the negative value does not matter for columns)
   */
  int append_blocks(const SparseCoeffmatMatrix& append_mat,
		    const CH_Matrix_Classes::Indexmatrix* blocks=0,
		    const CH_Matrix_Classes::Indexmatrix* cols=0);

  /** @brief append append_mat (or its submatrix given by blocks and/or cols) to the right
     
      If blocks or cols have negative entries, each of these represent a row of blocks of zeros whose order is the negative of the given value (if columns exist already, these need to match the order of the existing blocks) or a column of zeros of appropriate size (the size of the negative value does not matter for columns)
   */
  int append_columns(const SparseCoeffmatMatrix& append_mat,
		     const CH_Matrix_Classes::Indexmatrix* blocks=0,
		     const CH_Matrix_Classes::Indexmatrix* cols=0);
		    
  /// afterwards the new block i is the previous block map_to_old(i); no multiple appearances are allowed, but not all have to appear (these are deleted)
  int reassign_blocks(const CH_Matrix_Classes::Indexmatrix& map_to_old);
  /// afterwards the new column i is the previous column map_to_old(i); no multiple appearances are allowed, but not all have to appear (these are deleted)
  int reassign_columns(const CH_Matrix_Classes::Indexmatrix& map_to_old);
  
  /// generates a map_to_old by deleting the specified indices and calls reassign_blocks; this map_to_old will be returned if the corresponding pointer is not NULL 
  int delete_blocks(const CH_Matrix_Classes::Indexmatrix& delete_indices, CH_Matrix_Classes::Indexmatrix* map_to_old=0);
  /// generates a map_to_old by deleting the specified indices and calls reassign_columns; this map_to_old will be returned if the corresponding pointer is not NULL 
  int delete_columns(const CH_Matrix_Classes::Indexmatrix& delete_indices, CH_Matrix_Classes::Indexmatrix* map_to_old=0);
  
  /// returns the column representation of the matrix
  const SCMcolrep* get_colrep() const
  { form_colrep(); return colrep;}

  /// returns the number of dense matrices in block i
  CH_Matrix_Classes::Integer get_dense_cnt(CH_Matrix_Classes::Integer i) const
  { return dense_cnt(i);}

  /// useful for testing purposes; true iff both have the same size and both have the same pointers in the same positions
  bool operator==(const SparseCoeffmatMatrix& mat) const;
  
  /// useful for testing purposes; false iff both have the same size and both have the same pointers in the same positions
  bool operator!=(const SparseCoeffmatMatrix& mat) const
  {return !(*this==mat);}

  /// computes the inner products of (selected) columns (which represent block diagonal symmetric matrices) with the Gram matrix P*P^T (or, if Lam is given, P*Diag(Lam)*P^T) into the column vector ipvec(j)=ip(P*P^T,A.column((*ind)(j)}) (j=0,...,ind->dim()-1); if ind==NULL, use all columns
  int Gram_ip(CH_Matrix_Classes::Matrix& ipvec,const CH_Matrix_Classes::Matrix& P,const CH_Matrix_Classes::Matrix* Lam=0,const CH_Matrix_Classes::Indexmatrix* ind=0) const;

  /// computes the inner product of the block diagonal symmetric matrix stored in colummn j with the Gram matrix P*P^T into ipval=ip(P*P^T,A.column(j))
  int Gram_ip(CH_Matrix_Classes::Real& ipval,const CH_Matrix_Classes::Matrix& P,CH_Matrix_Classes::Integer j) const;

  /// computes the inner products of (selected) columns (which represent block diagonal symmetric matrices) with the primal into the column vector ipvec(j)=ip(*primal,A.column((*ind)(j)}) (j=0,...,ind->dim()-1); if ind==NULL, use all columns
  int primal_ip(CH_Matrix_Classes::Matrix& ipvec,const PSCPrimal* primal,const CH_Matrix_Classes::Indexmatrix* ind=0) const;

  /// computes the inner product of the block diagonal symmetric matrix stored in colummn j with the primal into ipval=ip(*primal,A.column(j))
  int primal_ip(CH_Matrix_Classes::Real& value,const PSCPrimal* primal,CH_Matrix_Classes::Integer j) const;
  
  /// computes S=P^T*A.column(j)*P
  int project(CH_Matrix_Classes::Symmatrix& S, const CH_Matrix_Classes::Matrix& P,const CH_Matrix_Classes::Integer j) const;
  
  
};


  //@}
}

#endif

