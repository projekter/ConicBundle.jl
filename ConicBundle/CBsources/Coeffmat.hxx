/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/Coeffmat.hxx
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



#ifndef CONICBUNDLE_COEFFMAT_HXX
#define CONICBUNDLE_COEFFMAT_HXX

/**  @file Coeffmat.hxx
    @brief Header declaring the classes ConicBundle::Coeffmat, ConicBundle::CoeffmatPointer, ConicBundle::CoeffmatInfo, ConicBundle::CMIname, and enum tpye ConicBundle::Coeffmattype  (needed for ConicBundle::PSCAffineFunction)
    @version 1.0
    @date 2017-02-03
    @author Christoph Helmberg
*/

#include <string>
#include <map>
#include "memarray.hxx"
#include "symmat.hxx"
#include "sparssym.hxx"
#include "CBout.hxx"


namespace ConicBundle {


  /** @ingroup implemented_psc_oracle
   */
   //@{

   /// for recognizing the type when writing and reading the problem
  enum Coeffmattype {
    CM_unspec = 0,      ///<any user defined constraint may use this
    CM_symdense = 1,    ///< for CMsymdense
    CM_symsparse = 2,   ///< for CMsymsparse
    CM_lowrankdd = 3,   ///< for CMlowrankdd
    CM_lowranksd = 4,   ///< for CMlowranksd
    CM_lowrankss = 5,   ///< for CMlowrankss
    CM_gramdense = 6,   ///< for CMgramdense
    CM_gramsparse = 7,  ///< for CMgramsparse
    CM_singleton = 8,   ///< for CMsingleton
    CM_gramsparsewd = 9 ///< for CMgramsparse_withoutdiag
  };

  //@}

  /** @ingroup implemented_psc_oracle
   */
   //@{

     /**@brief allows to memorize the scalings applied to a Coeffmat and offers the basis
        for storing further user defined informations on a Coeffmat
      */

  class CoeffmatInfo {
  protected:
    CH_Matrix_Classes::Real scalefactor; ///< scaling value
  public:
    /// default value is 1 for no scaling
    CoeffmatInfo(CH_Matrix_Classes::Real sf = 1.) :scalefactor(sf) {
    }
    ///
    virtual ~CoeffmatInfo() {
    }

    /// generates a new copy of itself on the heap
    virtual CoeffmatInfo* clone() const {
      return new CoeffmatInfo(scalefactor);
    }

    /// returns the scale factor
    CH_Matrix_Classes::Real get_scalefactor() const {
      return scalefactor;
    }
    /// sets the scale factor
    void set_scalefactor(CH_Matrix_Classes::Real sf) {
      scalefactor = sf;
    }
    /// scales the scale factor
    void multiply(CH_Matrix_Classes::Real sf) {
      scalefactor *= sf;
    }
    /// output a name if there is one for recognizing the type
    virtual std::ostream& print_id(std::ostream& out) const {
      return out << "NoInfo";
    }
  };

  //@}

  /** @ingroup implemented_psc_oracle
   */
   //@{

   /**@brief extends CoeffmatInfo to store a name (e.g. of the constraint it represents)
   */

  class CMIName : public CoeffmatInfo {
  private:
    std::string name; ///< well, this is the name
  public:
    /// initalize name and scale factor
    CMIName(std::string in_name, CH_Matrix_Classes::Real sf = 1.) :CoeffmatInfo(sf), name(in_name) {
    }
    ///
    ~CMIName() {
    }
    /// generates a new copy of itself on the heap
    CoeffmatInfo* clone() const {
      return new CMIName(name, scalefactor);
    }
    /// outputs the name
    std::ostream& print_id(std::ostream& out) const {
      return out << name.c_str();
    }
  };

  //@}

  /** @ingroup implemented_psc_oracle
   */
   //@{

     /// if cip is not zero, it calls and returns cip->clone() and 0 otherwise
  inline CoeffmatInfo* clone(const CoeffmatInfo* cip) {
    return (cip == 0) ? 0 : cip->clone();
  }


  //@}

  /** @ingroup implemented_psc_oracle
   */
   //@{

     /**@brief defines a base class for coefficient matrices in semidefinite programming, in particular for use with MatrixSDPfunction, see \ref implemented_psc_oracle.

        The class supports several different needs, e.g. it should work well
        with the requirements of MaxEigOracle as well as those of the
        eigenvalue computations with Bigmatrix. Furthermore there is
        support for writing and reading the various implementations.

      */

  class Coeffmat : protected CH_Matrix_Classes::Memarrayuser {
  private:
    friend class CoeffmatPointer;
    /// Each CoeffmatPointer pointing to his increase the use_cnt by one and it reduces it by one once it stops pointing to it; if use_cnt turns zero and deletion_by_CoeffmatPointer is true, the object will be deleted by the respective CoeffmatPointer
    mutable CH_Matrix_Classes::Integer use_cnt;
    /// set deletion_by_CoeffmatPointer==true if a CoeffmatPointer reducing the use_cnt to zero should delete this 
    bool deletion_by_CoeffmatPointer;
  protected:
    Coeffmattype CM_type;  ///< in order to enable type identification
    CoeffmatInfo* infop;  ///< allows the user to specify and output additional information
  public:
    /// default constructor; set del_by_CoeffmatPointer==true if a CoeffmatPointer reducing the use_cnt to zero should delete this 
    Coeffmat(bool del_by_CoeffmatPointer = true) {
      CM_type = CM_unspec; infop = 0; use_cnt = 0; deletion_by_CoeffmatPointer = del_by_CoeffmatPointer;
    }
    ///
    virtual ~Coeffmat() {
      delete infop;
    }

    ///if set to true any CoeffmatPointer that reduces use_cnt to zero will delete this object
    void set_deletion_by_CoeffmatPointer(bool dbMP) {
      deletion_by_CoeffmatPointer = dbMP;
    }
    /// returns its Coeffmattype
    virtual Coeffmattype get_type() const {
      return CM_type;
    }
    /// returns the user information
    virtual CoeffmatInfo* get_info() {
      return  infop;
    }
    /// returns the user information in const form 
    virtual const CoeffmatInfo* get_info() const {
      return  infop;
    }
    /// deletes the old and sets new user information 
    virtual void set_info(CoeffmatInfo* cip) {
      delete infop; infop = cip;
    }


    ///makes an explicit copy of itself and returns a pointer to it 
    virtual Coeffmat* clone() const = 0;

    ///returns the order of the represented symmetric matrix
    virtual CH_Matrix_Classes::Integer dim() const = 0;

    ///returns the value of the matrix element (i,j)
    virtual CH_Matrix_Classes::Real operator()(CH_Matrix_Classes::Integer i, CH_Matrix_Classes::Integer j) const = 0;

    ///returns a dense symmetric constraint matrix (useful for testing)
    virtual void make_symmatrix(CH_Matrix_Classes::Symmatrix& S) const = 0;

    ///returns the Frobenius norm of the matrix
    virtual CH_Matrix_Classes::Real norm(void) const = 0;


    ///delivers a new object on the heap corresponding to the matrix P^TAP, the caller is responsible for deleting the object
    virtual Coeffmat* subspace(const CH_Matrix_Classes::Matrix& P) const = 0;

    ///multiply constraint permanentely by d; this is to allow scaling or sign changes in the constraints
    virtual void multiply(CH_Matrix_Classes::Real d) = 0;

    ///returns ip(*this,S)=trace(*this*S), the trace inner product
    virtual CH_Matrix_Classes::Real ip(const CH_Matrix_Classes::Symmatrix& S) const = 0;

    ///returns ip(*this,PP^T)=trace P^T(*this)P
    virtual CH_Matrix_Classes::Real gramip(const CH_Matrix_Classes::Matrix& P) const = 0;

    ///returns ip(*this,Q*Diag(Lam)*Q^T)=trace Q^T(*this)Q*Diag(Lam) for Q=P.rows(start_row,start_row+dim-1), where Diag(Lam) is the identity if not given 
    virtual CH_Matrix_Classes::Real gramip(const CH_Matrix_Classes::Matrix& P, CH_Matrix_Classes::Integer start_row, const CH_Matrix_Classes::Matrix* Lam = 0) const = 0;

    ///computes S+=d*(*this);
    virtual void addmeto(CH_Matrix_Classes::Symmatrix& S, CH_Matrix_Classes::Real d = 1.) const = 0;

    ///comutes A+=d*(*this)*B
    virtual void addprodto(CH_Matrix_Classes::Matrix& A, const CH_Matrix_Classes::Matrix& B, CH_Matrix_Classes::Real d = 1.) const = 0;

    ///computes A+=d*(*this)*B
    virtual void addprodto(CH_Matrix_Classes::Matrix& A, const CH_Matrix_Classes::Sparsemat& B, CH_Matrix_Classes::Real d = 1.) const = 0;

    ///computes R=P^T*(*this)*Q
    virtual void left_right_prod(const CH_Matrix_Classes::Matrix& P, const CH_Matrix_Classes::Matrix& Q, CH_Matrix_Classes::Matrix& R) const = 0;

    ///returns an estimate of number of flops to compute addprodto for a vector
    virtual CH_Matrix_Classes::Integer prodvec_flops() const = 0;

    ///returns 1 if its structure is as bad as its dense symmetric representation, otherwise 0
    virtual int dense() const = 0;

    ///returns 0 if not sparse, otherwise 1
    virtual int sparse() const = 0;

    /// returns 0 if not sparse. If it is sparse it returns 1 and the nonzero structure in I,J and val, where val is multiplied by d. Only the upper triangle (including diagonal) is delivered
    virtual int sparse(CH_Matrix_Classes::Indexmatrix& I, CH_Matrix_Classes::Indexmatrix& J, CH_Matrix_Classes::Matrix& val, CH_Matrix_Classes::Real d = 1.)const = 0;

    /// returns 0 if the support of the costraint matrix is not contained in the support of the sparse symmetric matrix A, 1 if it is contained.
    virtual int support_in(const CH_Matrix_Classes::Sparsesym& A) const = 0;

    ///returns the inner product of the constraint matrix with A
    virtual CH_Matrix_Classes::Real ip(const CH_Matrix_Classes::Sparsesym& A) const = 0;

    ///computes S=P^T*(*this)*P
    virtual void project(CH_Matrix_Classes::Symmatrix& S, const CH_Matrix_Classes::Matrix& P) const = 0;

    ///computes S+=alpha*Q^T(*this)Q for Q=P.rows(start_row,start_row+dim-1) 
    virtual void add_projection(CH_Matrix_Classes::Symmatrix& S, const CH_Matrix_Classes::Matrix& P, CH_Matrix_Classes::Real alpha = 1., CH_Matrix_Classes::Integer start_row = 0) const = 0;

    ///computes C= alpha*(*this)*B^(T if btrans) + beta*C, C is also returned
    virtual const CH_Matrix_Classes::Matrix& postgenmult(const CH_Matrix_Classes::Matrix& B, CH_Matrix_Classes::Matrix& C,
      CH_Matrix_Classes::Real alpha = 1., CH_Matrix_Classes::Real beta = 0., int btrans = 0) const = 0;

    ///computes C= alpha*B^(T if btrans)*(*this) + beta*C, C is also returned
    virtual const CH_Matrix_Classes::Matrix& pregenmult(const CH_Matrix_Classes::Matrix& B, CH_Matrix_Classes::Matrix& C,
      CH_Matrix_Classes::Real alpha = 1., CH_Matrix_Classes::Real beta = 0., int btrans = 0) const = 0;

    ///returns 1, if p is the same derived class and entries differ by less than tol, otherwise zero
    virtual int equal(const Coeffmat* p, double tol = 1e-6) const = 0;

    ///display constraint information
    virtual std::ostream& display(std::ostream& o) const = 0;

    ///put entire contents onto outstream with the class type in the beginning so that the derived class can be recognized by in().
    virtual std::ostream& out(std::ostream& o) const = 0;

    ///counterpart to out(), does not read the class type, though. This is assumed to have been read in order to generate the correct class
    virtual std::istream& in(std::istream& i) = 0;

  };

  //@}

  /** @ingroup implemented_psc_oracle
   */
   //@{

     /**@brief reads the next Coeffmat from in into an object on the heap and returns a pointer to it. The caller has to destruct the object.
      */
  Coeffmat* coeffmat_read(std::istream& in);

  //@}

  /** @ingroup implemented_psc_oracle
   */
   //@{

     /**@brief pointer class for Coeffmat for deleting objects on the heap if Coefmat::use_cnt is reduced to zero and deletion is allowed.
      */

  class CoeffmatPointer {
  private:
    /// holds the pointer to Coeffmat, may be NULL == empty
    Coeffmat* coeffmatp;

  public:
    /// initialize the Pointer to hold cfmp (==NULL is allowed for empty); for the coeffmat stored previously it reduces the use_cnt and deletes it if zero.
    void init(Coeffmat* cfmp) {
      if ((coeffmatp) &&
        (--(coeffmatp->use_cnt) <= 0) &&
        (coeffmatp->deletion_by_CoeffmatPointer)) {
        assert(coeffmatp->use_cnt == 0);
        delete coeffmatp;
      }
      coeffmatp = cfmp;
      if (coeffmatp)
        coeffmatp->use_cnt++;
    }

    /// default initialization to an empty pointer 
    CoeffmatPointer() :coeffmatp(0) {
    }

    /// calls init() (cfmp==NULL is allowed)
    CoeffmatPointer(Coeffmat* cfmp) :coeffmatp(0) {
      init(cfmp);
    }

    /// calls init(mp.coeffmatp)
    CoeffmatPointer(const CoeffmatPointer& mp) :coeffmatp(0) {
      init(mp.coeffmatp);
    }

    /// calls init(0) first
    ~CoeffmatPointer() {
      init(0);
    }

    /// calls init(mp.coeffmatp)
    CoeffmatPointer& operator=(const CoeffmatPointer& mp) {
      init(mp.coeffmatp); return *this;
    }

    /// calls init(cfmp), NULL is allowed
    CoeffmatPointer& operator=(Coeffmat* cfmp) {
      init(cfmp); return *this;
    }

    /// compares the pointers
    bool operator==(const CoeffmatPointer& mp) const {
      return coeffmatp == mp.coeffmatp;
    }

    /// compares the pointers
    bool operator==(const Coeffmat* cfmp) const {
      return coeffmatp == cfmp;
    }

    /// compares the pointers
    bool operator!=(const CoeffmatPointer& mp) const {
      return coeffmatp != mp.coeffmatp;
    }

    /// compares the pointers
    bool operator!=(const Coeffmat* cfmp) const {
      return coeffmatp != cfmp;
    }

    /// compares the pointers
    bool operator<(const CoeffmatPointer& mp) const {
      return coeffmatp < mp.coeffmatp;
    }

    /// compares the pointers
    bool operator<(const Coeffmat* cfmp) const {
      return coeffmatp < cfmp;
    }

    /// compares the pointers
    bool operator>(const CoeffmatPointer& mp) const {
      return coeffmatp > mp.coeffmatp;
    }

    /// compares the pointers
    bool operator>(const Coeffmat* cfmp) const {
      return coeffmatp > cfmp;
    }

    /// returns the object pointed to (may not be called if empty)
    Coeffmat& operator*() const {
      assert(coeffmatp); return *coeffmatp;
    }

    /// returns the pointers (may not be called if empty)
    Coeffmat* operator->() const {
      assert(coeffmatp); return coeffmatp;
    }

    /// returns the pointer (may return NULL)
    Coeffmat* ptr() const {
      return coeffmatp;
    }

  };

  //@}



}

#endif

