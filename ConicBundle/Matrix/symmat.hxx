/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  Matrix/symmat.hxx
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



#ifndef CH_MATRIX_CLASSES__SYMMAT_HXX
#define CH_MATRIX_CLASSES__SYMMAT_HXX

/**  @file symmat.hxx
    @brief Header declaring the class CH_Matrix_Classes::Symmatrix for symmetric matrices with Real elements
    @version 1.0
    @date 2005-03-01
    @author Christoph Helmberg

*/

#ifndef CH_MATRIX_CLASSES__MATRIX_HXX
#include "matrix.hxx"
#endif

namespace CH_Matrix_Classes {


  // **************************************************************************
  //                            class definition
  // **************************************************************************

  /**@defgroup Symmatrixgroup Symmatrix (dense, real, symmetric, n by n)
  */
  //@{

  /** @brief %Matrix class of symmetric matrices with real values of type #Real

      Internally a symmetric matrix of size nr x nr is stored in a one dimensional
      array of ::Real variables, the elements of the lower triangle
      are arranged in columnwise order (a11,a21,...,anr1,a22,a32,...).

      Any matrix element can be indexed by (i,j) or directly by the one dimensional
      index (i+j*nr). The latter view directly corresponds to the vec() operator
      often used in the linear algebra literature, i.e., the matrix is
      transformed to a vector by stacking the columns on top of each other.

      NOTE: Any change of A(i,j) also changes A(j,i) as both variables are identical!

   */
  class Symmatrix : protected Memarrayuser {
    friend class Matrix;
    friend class Sparsesym;
    friend class Sparsemat;

  private:

    static const Mtype mtype;   ///< used for MatrixError templates (runtime type information was not yet existing)
    Integer mem_dim;   ///< amount of memory currently allocated
    Integer nr;        ///< number of rows = number of columns
    Real* m;     ///< lower triangle stored columnwise (a11,a21,...,anr1,a22,.....)

    bool is_init;   ///< flag whether memory is initialized, it is only used if CONICBUNDLE_DEBUG is defined

    /// initialize the matrix to a 0x0 matrix without storage
    inline void init_to_zero();

    /// a subroutine needed internally for eigenvalue computations (eigval.cxx)
    Integer tred2(Integer nm, Integer n, Real* a, Real* d, Real* e, Real* z) const;
    /// a subroutine needed internally for eigenvalue computations (eigval.cxx)
    Integer imtql2(Integer nm, Integer n, Real* d, Real* e, Real* z) const;

  public:

    //----------------------------------------------------
    //----  constructors, destructor, and initialization
    //----------------------------------------------------


    /** @name Constructors, Destructor, and Initialization (Members)
     */
     //@{

     /// empty matrix
    inline Symmatrix();
    /// copy constructor, *this=d*A
    inline Symmatrix(const Symmatrix& A, double d = 1.);
    /** @brief generate a matrix of size nr x nr but WITHOUT initializing the memory

        If initializing the memory externally and CONICBUNDLE_DEBUG is defined, please use
        set_init() via matrix.set_init(true) in order to avoid warnings concerning improper
        initialization
    */
    inline Symmatrix(Integer nr);
    /// generate a matrix of size nr x nr initializing all elements to the value d
    inline Symmatrix(Integer nr, Real d);
    /// generate a matrix of size nr x nr initializing the elements from the (one dimensional) array dp, which must have the elements arranged consecutively in internal order
    inline Symmatrix(Integer nr, const Real* dp);
    ///
    inline ~Symmatrix();

#if (CONICBUNDLE_DEBUG>=1)
    /// after external initialization, call matrix.set_init(true) (not needed if CONICBUNDLE_DEBUG is undefined)
    void set_init(bool i) {
      is_init = i;
    }
    /// returns true if the matrix has been declared initialized (not needed if CONICBUNDLE_DEBUG is undefined)
    int get_init() const {
      return is_init;
    }
#else 
    /// after external initialization, call matrix.set_init(true) (not needed if CONICBUNDLE_DEBUG is undefined) 
    void set_init(bool /* i */) {
    }
    /// returns true if the matrix has been declared initialized (not needed if CONICBUNDLE_DEBUG is undefined)
    bool get_init() const {
      return true;
    }
#endif

    /// initialize to *this=A*d
    inline Symmatrix& init(const Symmatrix& A, double d = 1.);
    /// initialize to *this=d*(A+transpose(A))/2.
    inline Symmatrix& init(const Matrix& A, double d = 1.);
    /// initialize to *this=d*(A+transpose(A))/2.
    inline Symmatrix& init(const Indexmatrix& A, double d = 1.);
    /// initialize to *this=A*d
    inline Symmatrix& init(const Sparsesym& A, Real d = 1.);
    /// intialize *this to a matrix of size nr x nr initializing all elements to the value d
    inline Symmatrix& init(Integer nr, Real d);
    /// generate a matrix of size nr x nc initializing the elements from the (one dimensional) array dp which must have the elements arranged consecutively in internal order
    inline Symmatrix& init(Integer nr, const Real* dp);

    /** @brief resize the matrix to nr x nr elements but WITHOUT initializing the memory

        If initializing the memory externally and CONICBUNDLE_DEBUG is defined, please use
        set_init() via matrix.set_init(true) in order to avoid warnings concerning improper
        initialization
    */
    void newsize(Integer n);      //resize matrix without Initialization


    //@}

    /** @name Conversions from other Matrix Classes (Members)
     */
     //@{

     /// (*this)=d*(A+transpose(A))/2.
    inline Symmatrix(const Matrix&, double d = 1.);
    /// (*this)=d*(A+transpose(A))/2.
    inline Symmatrix(const Indexmatrix&, double d = 1.);
    /// (*this)=d*A
    inline Symmatrix(const Sparsesym& A, Real d = 1.);

    //@}

    //----------------------------------------------------
    //----  size and type information
    //----------------------------------------------------

    /** @name Size and Type Information (Members)
     */
     //@{

     /// returns the number of rows in _nr and _nc
    void dim(Integer& _nr, Integer& _nc) const {
      _nr = nr; _nc = nr;
    }

    /// returns the dimension rows * columns when the matrix is regarded as a vector
    Integer dim() const {
      return nr;
    }

    /// returns the row dimension
    Integer rowdim() const {
      return nr;
    }

    /// returns the column dimension
    Integer coldim() const {
      return nr;
    }

    /// returns the type of the matrix, MTsymmetric
    Mtype get_mtype() const {
      return mtype;
    }

    //@}


    //--------------------------------
    //----  Indexing and Submatrices
    //--------------------------------

    /** @name Indexing and Submatrices (Members)
     */
     //@{

     /// returns reference to element (i,j) of the matrix (rowindex i, columnindex j)
    inline Real& operator()(Integer i, Integer j);

    /// returns reference to element (i) of the matrix if regarded as vector of stacked columns [element (i%rowdim, i/rowdim)]
    inline Real& operator()(Integer i);

    /// returns value of element (i,j) of the matrix (rowindex i, columnindex j)
    inline Real operator()(Integer i, Integer j) const;

    /// returns value of element (i) of the matrix if regarded as vector of stacked columns [element (i%rowdim, i/rowdim)]
    inline Real operator()(Integer i) const;

    /// returns column i copied to a new Matrix 
    Matrix col(Integer i) const;
    /// returns row i copied to a new Matrix 
    Matrix row(Integer i) const;
    /// returns a matrix of size this->rowdim() x vec.dim(), with column i a copy of column vec(i) of *this
    Matrix cols(const Indexmatrix& vec) const;  //returns cols as indexed by vec
    /// returns a matrix of size vec.dim() x this->coldim(), with row i a copy of row vec(i) of *this
    Matrix rows(const Indexmatrix& vec) const;  //returns rows as indexed by vec

    /// swaps rows (and columns) i and j
    Symmatrix& swapij(Integer i, Integer j);

    /// for i=0 to rowdim row (and column) i of this matrix is swapped with row piv(j); for inverse=true the inverse permutation is generated
    Symmatrix& pivot_permute(const Indexmatrix& piv, bool inverse = false);

    ///returns S and in S the principal submatrix indexed by ind (multiple indices are allowed)
    Symmatrix& principal_submatrix(const Indexmatrix& ind, Symmatrix& S) const;

    ///returns the principal submatrix indexed by ind (multiple indices are allowed)
    inline Symmatrix principal_submatrix(const Indexmatrix& ind) const;


    ///returns this afte deleting the principal submatrix indexed by ind (no repetitions!); 
    Symmatrix& delete_principal_submatrix(const Indexmatrix& ind, bool sorted_increasingly = false);

    ///increases the order of the matrix by appending storage for further addn rows and columns (marked as not initiliazed if addn>0, no changes if addn<=0)
    Symmatrix& enlarge_below(Integer addn);

    ///increases the order of the matrix by appending storage for further addn rows and columns initialized to d (no changes if addn<=0);
    Symmatrix& enlarge_below(Integer addn, Real d);



    /// returns the current address of the internal value array; use cautiously, do not use delete!
    Real* get_store() {
      return m;
    }   //use cautiously, do not use delete!
/// returns the current address of the internal value array; use cautiously!
    const Real* get_store() const {
      return m;
    }   //use cautiously

//@}


/** @name Indexing and Submatrices (Friends)
 */
 //@{

 /// returns a column vector v consisting of the elements v(i)=(*this)(i,i), 0<=i<row dimension 
    friend Matrix diag(const Symmatrix& A);    //=(A(1,1),A(2,2),...)^t
    /// returns a symmetric diagonal matrix S of order A.dim() with vec(A) on the diagonal, i.e., S(i,i)=A(i) for all i and S(i,j)=0 for i!=j 
    friend Symmatrix Diag(const Matrix& A);    //vec(A) on the diagonal

    /// swap the content of the two matrices A and B (involves no copying)
    friend inline void swap(Symmatrix& A, Symmatrix& B);

    //@}

    //------------------------------
    //----  BLAS-like Routines
    //------------------------------

    /** @name BLAS-like Routines (Members)
     */
     //@{

     ///sets *this=d*A and returns *this
    Symmatrix& xeya(const Symmatrix& A, Real d = 1.);   //*this=d*A
    ///sets *this+=d*A and returns *this  
    Symmatrix& xpeya(const Symmatrix& A, Real d = 1.);  //*this+=d*A;


    //@}

    /** @name BLAS-like Routines (Friends)
     */
     //@{

     /// returns C=beta*C+alpha* A*A^T, where A may be transposed. If beta==0. then C is initiliazed to the correct size.
    friend Symmatrix& rankadd(const Matrix& A, Symmatrix& C,
      Real alpha, Real beta, int trans);
    /// returns C=beta*C+alpha* A*D*A^T, where D is a vector representing a diagonal matrix and A may be transposed; if beta==0. then C is initialized to the correct size
    friend Symmatrix& scaledrankadd(const Matrix& A, const Matrix& D, Symmatrix& C,
      Real alpha, Real beta, int trans);
    /// returns C=beta*C+alpha*(A*B^T+B*A^T)/2 [or for transposed (A^T*B+B^T*A)/2]. If beta==0. then C is initiliazed to the correct size.
    friend Symmatrix& rank2add(const Matrix& A, const Matrix& B, Symmatrix& C,
      Real alpha, Real beta, int trans);


    ///returns x= alpha*y+beta*x; if beta==0. then x is initialized to the correct size
    friend Symmatrix& xbpeya(Symmatrix& x, const Symmatrix& y, Real alpha, Real beta);


    ///returns x= alpha*y+beta*z; x is initialized to the correct size
    friend Symmatrix& xeyapzb(Symmatrix& x, const Symmatrix& y, const Symmatrix& z, Real alpha, Real beta);

    ///returns C=beta*C+alpha*A*B, where B may be transposed; C must not be equal to B; if beta==0. then C is initialized to the correct size
    friend Matrix& genmult(const Symmatrix& A, const Matrix& B, Matrix& C,
      Real alpha, Real beta, int btrans);

    ///returns C=beta*C+alpha*A*B, where A may be transposed; C must not be equal to A; if beta==0. then C is initialized to the correct size
    friend Matrix& genmult(const Matrix& A, const Symmatrix& B, Matrix& C,
      Real alpha, Real beta, int atrans);

    //@}

    //------------------------------
    //----  usual operators
    //------------------------------

    /** @name Usual Arithmetic Operators (Members)
     */
     //@{

     ///
    inline Symmatrix& operator=(const Symmatrix& A);
    ///
    inline Symmatrix& operator+=(const Symmatrix& A);
    ///
    inline Symmatrix& operator-=(const Symmatrix& A);
    /// ATTENTION: this is redefined as the Hadamard product, (*this)(i,j)=(*this)(i,j)*A(i,j) for all i<=j
    inline Symmatrix& operator%=(const Symmatrix& A);  //Hadamard product
    ///
    inline Symmatrix operator-() const;

    /// 
    inline Symmatrix& operator*=(Real d);
    /// ATTENTION: d is NOT checked for 0
    inline Symmatrix& operator/=(Real d);
    /// sets (*this)(i,j)+=d for all i<=j
    inline Symmatrix& operator+=(Real d);
    /// sets (*this)(i,j)-=d for all i<=j
    inline Symmatrix& operator-=(Real d);

    ///transposes itself (at almost no cost)  
    Symmatrix& transpose() {
      return *this;
    }

    //@}

    /** @name Usual Arithmetic Operators (Friends)
     */
     //@{

     /// 
    friend inline Matrix operator*(const Symmatrix& A, const Symmatrix& B);
    /// ATTENTION: this is redefined as the Hadamard product and sets (i,j)=A(i,j)*B(i,j) for all i<=j  
    friend inline Symmatrix operator%(const Symmatrix& A, const Symmatrix& B);
    ///
    friend inline Symmatrix operator+(const Symmatrix& A, const Symmatrix& B);
    ///
    friend inline Symmatrix operator-(const Symmatrix& A, const Symmatrix& B);
    ///
    friend inline Matrix operator*(const Symmatrix& A, const Matrix& B);
    ///
    friend inline Matrix operator*(const Matrix& A, const Symmatrix& B);
    ///
    friend inline Matrix operator+(const Symmatrix& A, const Matrix& B);
    ///
    friend inline Matrix operator+(const Matrix& A, const Symmatrix& B);
    ///
    friend inline Matrix operator-(const Symmatrix& A, const Matrix& B);
    ///
    friend inline Matrix operator-(const Matrix& A, const Symmatrix& B);

    ///
    friend inline Symmatrix operator*(const Symmatrix& A, Real d);
    ///
    friend inline Symmatrix operator*(Real d, const Symmatrix& A);
    ///
    friend inline Symmatrix operator/(const Symmatrix& A, Real d);
    /// returns (i,j)=A(i,j)+d for all i<=j
    friend inline Symmatrix operator+(const Symmatrix& A, Real d);
    /// returns (i,j)=A(i,j)+d for all i<=j
    friend inline Symmatrix operator+(Real d, const Symmatrix& A);
    /// returns (i,j)=A(i,j)-d for all i<=j
    friend inline Symmatrix operator-(const Symmatrix& A, Real d);
    /// returns (i,j)=d-A(i,j) for all i<=j
    friend inline Symmatrix operator-(Real d, const Symmatrix& A);

    /// (drop it or use a constructor instead)
    friend inline Symmatrix transpose(const Symmatrix& A);

    //@}

    //------------------------------------------
    //----  Connections to other Matrix Classes
    //------------------------------------------

    /** @name Connections to other Classes (Members)
     */
     //@{

     ///sets *this=d*(A+transpose(A))/2. and returns *this
    Symmatrix& xeya(const Matrix& A, Real d = 1.);   //*this=d*A
    ///sets *this+=d*(A+transpose(A))/2. and returns *this
    Symmatrix& xpeya(const Matrix& A, Real d = 1.);  //*this+=d*A;
    ///sets *this=d*(A+transpose(A))/2. and returns *this
    Symmatrix& xeya(const Indexmatrix& A, Real d = 1.);   //*this=d*A
    ///sets *this+=d*(A+transpose(A))/2. and returns *this
    Symmatrix& xpeya(const Indexmatrix& A, Real d = 1.);  //*this+=d*A;
    ///sets *this=d*A and returns *this 
    Symmatrix& xeya(const Sparsesym& A, Real d = 1.);   //*this=d*A
    ///sets *this+=d*A and returns *this 
    Symmatrix& xpeya(const Sparsesym& A, Real d = 1.);  //*this+=d*A;

    ///sets *this(i,j), i<=j to the upper triangle of the matrix product d*transpose(A)*B
    Symmatrix& xetriu_yza(const Matrix& A, const Matrix& B, Real d = 1.);
    ///adds to *this(i,j), i<=j the upper triangle of the matrix product d*transpose(A)*B
    Symmatrix& xpetriu_yza(const Matrix& A, const Matrix& B, Real d = 1.);
    ///sets *this(i,j), i<=j to the upper triangle of the matrix product d*transpose(A)*B  
    Symmatrix& xetriu_yza(const Sparsemat& A, const Matrix& B, Real d = 1.);
    ///adds to *this(i,j), i<=j the upper triangle of the matrix product d*transpose(A)*B
    Symmatrix& xpetriu_yza(const Sparsemat& A, const Matrix& B, Real d = 1.);

    ///
    inline Symmatrix& operator=(const Sparsesym& A);
    ///
    inline Symmatrix& operator+=(const Sparsesym& A);
    ///
    inline Symmatrix& operator-=(const Sparsesym& A);

    //@}

    /** @name Connections to other Classes (Friends)
     */
     //@{

     ///returns C=beta*C+alpha*A*B, where B may be transposed; if beta==0. then C is initialized to the correct size
    friend Matrix& genmult(const Symmatrix& A, const Sparsemat& B, Matrix& C,
      Real alpha, Real beta, int btrans);

    ///returns C=beta*C+alpha*A*B, where A may be transposed; if beta==0. then C is initialized to the correct size
    friend Matrix& genmult(const Sparsemat& A, const Symmatrix& B, Matrix& C,
      Real alpha, Real beta, int atrans);

    /// returns C=beta*C+alpha* A*A^T, where A may be transposed; if beta==0. then C is initialized to the correct size
    friend Symmatrix& rankadd(const Sparsemat& A, Symmatrix& C,
      Real alpha, Real beta, int trans);

    /// returns C=beta*C+alpha* A*D*A^T, where D is a vector representing a diagonal matrix and A may be transposed; if beta==0. then C is initialized to the correct size
    friend Symmatrix& scaledrankadd(const Sparsemat& A, const Matrix& D, Symmatrix& C,
      Real alpha, Real beta, int trans);

    /// returns C=beta*C+alpha*(A*B^T+B*A^T)/2 [or for transposed (A^T*B+B^T*A)/2]. If beta==0. then C is initiliazed to the correct size.
    friend Symmatrix& rank2add(const Sparsemat& A, const Matrix& B, Symmatrix& C,
      Real alpha, Real beta, int trans);

    //@}

   //------------------------------
   //----  Elementwise Operations
   //------------------------------


   /** @name Elementwise Operations (Friends)
    */
    //@{

    /// returns a Symmatrix with elements abs(A(i,j))
    friend Symmatrix abs(const Symmatrix& A);

    //@}

    //----------------------------
    //----  Numerical Methods
    //----------------------------

    /** @name Numerical Methods (Members)
     */
     //@{

     /// shifts the diagonal by s, i.e., (*this)(i,i)+=s for all i
    Symmatrix& shift_diag(Real s);             //add s to all diagonal elements

    //----- LDL Factorization (for positive definite matrices, no pivots so far!)
    /// computes LDLfactorization (implemented only for positive definite matrices so far, no pivoting), (*this) is overwritten by the factorization; returns 1 if diagonal elements go below tol
    int LDLfactor(Real tol = 1e-10);
    /// computes, after LDLfactor was executed succesfully, the solution to (*old_this)x=rhs; rhs is overwritten by the solution; always returns 0; NOTE: there is NO check against division by zero  
    int LDLsolve(Matrix& x) const;      //call LDLfactor before using LDLsolve!
    /// computes, after LDLfactor was executed succesfully, the inverse to (*old_this) and stores it in S (numerically not too wise); always returns 0; NOTE: there is NO check against division by zero  
    int LDLinverse(Symmatrix& S) const; //call LDLfactor before using LDLinverse!

    //----- Cholesky Factorization with pivoting
    /// computes the Cholesky factorization, for positive definite matrices only, (*this) is overwritten by the factorization; there is no pivoting; returns 1 if diagonal elements go below tol
    int Chol_factor(Real tol = 1e-10); //stores fact. in *this
    /// computes, after Chol_factor was executed succesfully, the solution to (*old_this)x=rhs; rhs is overwritten by the solution; always returns 0; NOTE: there is NO check against division by zero  
    int Chol_solve(Matrix& x) const; //call _factor before
    /// computes, after Chol_factor was executed succesfully, the inverse to (*old_this) and stores it in S (numerically not too wise); always returns 0; NOTE: there is NO check against division by zero  
    int Chol_inverse(Symmatrix& S) const; // --- " ---
    /// computes, after Chol_factor into LL^T was executed succesfully, the solution to Lx=rhs; rhs is overwritten by the solution; always returns 0; NOTE: there is NO check against division by zero  
    int Chol_Lsolve(Matrix& rhs) const; // --- " ---, solve only Ly=x
    /// computes, after Chol_factor into LL^T was executed succesfully, the solution to L^Tx=rhs; rhs is overwritten by the solution; always returns 0; NOTE: there is NO check against division by zero  
    int Chol_Ltsolve(Matrix& rhs) const; // --- " ---, solve only L^Ty=x
    /// computes, after Chol_factor into LL^T was executed succesfully, L^{-1}SL^{-T} overwriting S
    int Chol_scaleLi(Symmatrix& S) const;
    /// computes, after Chol_factor into LL^T was executed succesfully, L^TSL overwriting S
    int Chol_scaleLt(Symmatrix& S) const;
    /// computes, after Chol_factor into LL^T was executed succesfully,  L*rhs, overwriting rhs by the result; always returns 0;  
    int Chol_Lmult(Matrix& rhs) const;
    /// computes, after Chol_factor into LL^T was executed succesfully,  L^Trhs, overwriting rhs by the result; always returns 0;  
    int Chol_Ltmult(Matrix& rhs) const;
    /// computes the Cholesky factorization with pivoting, for positive semidefinite matrices only, (*this) is overwritten by the factorization; on termination piv.dim() is the number of positive pivots>=tol; returns 1 if negative diagonal element is encountered during computations, 0 otherwise.
    int Chol_factor(Indexmatrix& piv, Real tol = 1e-10); //stores fact. in *this
    /// computes, after Chol_factor(Indexmatrix&,Real) with pivoting was executed succesfully, the solution to (*old_this)*x=rhs(piv); rhs is overwritten by the solution arranged in original unpermuted order; always returns 0; NOTE: there is NO check against division by zero  
    int Chol_solve(Matrix& x, const Indexmatrix& piv) const; //call _factor before
    /// computes, after Chol_factor(Indexmatrix&,Real) with pivoting was executed succesfully, the inverse to (*old_this) and stores it in S (the pivoting permutation is undone in S); NOTE: there is NO check against division by zero  
    int Chol_inverse(Symmatrix& S, const Indexmatrix& piv) const; // --- " ---

    //----- Aasen factorization with pivoting
    /// computes Aasen factorization LTL^T with pivoting, where L is unit lower triangular with first colum e_1 and T is tridiagonal; (*this) is overwritten by the factorization, with column i of L being stored in column i-1 of (*this); always returns 0;
    int Aasen_factor(Indexmatrix& piv);
    /// computes, after Aasen_factor into LTL^T was executed, the solution to Lx=rhs; rhs is overwritten by the solution; always returns 0; 
    int Aasen_Lsolve(Matrix& x) const;
    /// computes, after Aasen_factor into LTL^T was executed, the solution to L^Tx=rhs; rhs is overwritten by the solution; always returns 0; 
    int Aasen_Ltsolve(Matrix& x) const;
    /// computes, after Aasen_factor into LTL^T was executed, the solution to Tx=rhs; rhs is overwritten by the solution;if the solution fails due to division by zero (=system not solvable) the return value is -(rowindex+1) where this occured in the backsolve
    int Aasen_tridiagsolve(Matrix& x) const;
    /// computes, after Aasen_factor into LTL^T was executed, the solution to (*old_this)x=rhs; rhs is overwritten by the solution; if the solution fails due to division by zero (=system not solvable) the return value is -(rowindex+1) where this occured in the backsolve
    int Aasen_solve(Matrix& x, const Indexmatrix& piv) const; //call Aasen_factor before!

    //----- Eigenvalue decomposition
    /// computes an eigenvalue decomposition P*Diag(d)*tranpose(P)=(*this) by symmetric QR; returns 0 on success, 
    Integer eig(Matrix& P, Matrix& d, bool sort_non_decreasingly = true) const;
    //if return value ==0 then P contains eigenvectors as columns
    //and d is a column vector containing the eigenvalues.

    //@}

    /** @name Numerical Methods (Friends)
     */
     //@{

     /// returns the sum of the diagonal elements A(i,i) over all i
    friend Real trace(const Symmatrix& A);                  //=sum(diag(A))
    ///returns the usual inner product of A and B, i.e., the sum of A(i,j)*B(i,j) over all i,j 
    friend Real ip(const Symmatrix& A, const Symmatrix& B); //=trace(B^t*A)
    ///returns the usual inner product of A and B, i.e., the sum of A(i,j)*B(i,j) over all i,j 
    friend Real ip(const Matrix& A, const Symmatrix& B);    //=trace(B^t*A)
    ///returns the usual inner product of A and B, i.e., the sum of A(i,j)*B(i,j) over all i,j 
    friend Real ip(const Symmatrix& A, const Matrix& B);    //=trace(B^t*A)
    ///returns the Frobenius norm of A, i.e., the square root of the sum of A(i,j)*A(i,j) over all i,j 
    friend inline Real norm2(const Symmatrix& A); //{return sqrt(ip(A,A));}

    ///returns a row vector holding the sum over all rows, i.e., (1 1 ... 1)*A
    friend Matrix sumrows(const Symmatrix& A); //=(1 1 1 ... 1)*A
    ///returns a column vector holding the sum over all columns, i.e., A*(1 1 ... 1)^T
    friend Matrix sumcols(const Symmatrix& A); //=A*(1 1 ... 1)^t
    ///returns the sum over all elements of A, i.e., (1 1 ... 1)*A*(1 1 ... 1)^T
    friend Real sum(const Symmatrix& A);       //=(1 1 ... 1)*A*(1 1 ... 1)^t

    //---- svec and skron

    /// the symmetric vec operator stacks the lower triangle of A to a n*(n+1)/2 vector with the same norm2 as A; here it sets svec(A)=[a11,sqrt(2)a12,...,sqrt(2)a1n,a22,...,sqrt(2)a(n-1,n),ann]', multiplies it by a and sets or adds (if add==true) it to v starting from startindex_vec possibly restricted to the subblock of order blockdim (whenever >=0, else blockdim is set to A.rowdim()-startindex_A) starting from startindex_A (must be >=0); if add==false and startindex_vec<0 then vec is also reinitialzed to the appropriate size
    friend void svec(const Symmatrix& A, Matrix& v, Real a, bool add, Integer startindex_vec, Integer startindex_A, Integer blockdim);
    /// the symmetric vec operator, stacks the lower triangle of A to a n*(n+1)/2 vector with the same norm2 as A; i.e., it returns svec(A)=[a11,sqrt(2)a12,...,sqrt(2)a1n,a22,...,sqrt(2)a(n-1,n),ann]'
    friend inline Matrix svec(const Symmatrix& A);
    /// the inverse operator to svec, extracts from v at startindex_vec (>=0) the symmetric matrix of blockdim adding its mutliple by a into A starting at startindex_A; if add==false and startindex_A<0 A is initialized to the size of blockdim; if the latter is also negative then v.dim()-startindex_vec must match an exact order and  matrix A is initialized to this size. In all other cases the size of the symmetric matrix determines the missing parameters and vec.dim-startindex_vec 
    friend void sveci(const Matrix& v, Symmatrix& A, Real a, bool add, Integer startindex_vec, Integer startindex_A, Integer blockdim);

    // initialize from an svec stored in a real array (or matrix)
    void init_svec(Integer nr, const Real* dp, Integer incr = 1, Real d = 1.);
    // store in the form of an svec in the real array (or matrix)
    void store_svec(Real* dp, Integer incr = 1, Real d = 1.) const;

    /// the symmetric Kronecker product, defined via (A skron B)svec(C)=(BCA'+ACB')/2; sets or adds (if add==true) the symmetric matrix a*(A skron B) into S starting at startindex_S; if add==false and startindex_S<0, S is initialzed to the correct size 
    friend inline Symmatrix skron(const Symmatrix& A, const Symmatrix& B, Real alpha, bool add, Integer startindex_S);
    /// def symmetric Kronecker product (A skron B)svec(C)=(BCA'+ACB')/2; sets S=alpha*(A skron B) or S*=... (if add==true) possibly shifted to the block starting at startindex_S;  if add==false and startindex_S<0, S is initialzed to the correct size 
    friend void skron(const Symmatrix& A, const Symmatrix& B, Symmatrix& S, Real a, bool add, Integer startindex_S);
    /// sets S=beta*S+alpha*B'*A*B for symmatrix A and matrix B
    friend void symscale(const Symmatrix& A, const Matrix& B, Symmatrix& S, Real a, Real b, int btrans);

    //@}

    //---------------------------------------------
    //----  Comparisons / Max / Min / sort / find
    //---------------------------------------------

    /** @name Comparisons, Max, Min, Sort, Find (Friends)
     */
     //@{

      /// returns a row vector holding in each column the minimum over all rows in this column
    friend Matrix minrows(const Symmatrix& A); //minimum of each column (over rows)
    /// returns a column vector holding in each row the minimum over all columns in this row
    friend Matrix mincols(const Symmatrix& A); //minimum of each row (over colums)
    /// returns the minimum value over all elements of the matrix
    friend Real min(const Symmatrix& A);       //minimum of matrix
    /// returns a row vector holding in each column the maximum over all rows in this column
    friend Matrix maxrows(const Symmatrix& A); //similar
    /// returns a column vector holding in each row the maximum over all columns in this row
    friend Matrix maxcols(const Symmatrix& A); //similar
    /// returns the maximum value over all elements of the matrix
    friend Real max(const Symmatrix& A);       //similar


    //@}

    //--------------------------------
    //----  Input / Output
    //--------------------------------


    /** @name Input, Output (Members)
     */
     //@{

     /** @brief displays a matrix in a pretty way for bounded screen widths; for variables of value zero default values are used.
      */
    void display(std::ostream& out, ///< output stream
      int precision = 0,   ///< number of most significant digits, default=4
      int width = 0,       ///< field width, default = precision+6
      int screenwidth = 0  ///< maximum number of characters in one output line, default = 80
    ) const;

    /** @brief outputs a matrix A in the format "[ A(0,1) ... A(0,nc-1)\n ... A(nr-1,nc-1)];\n" so that it can be read e.g. by octave as an m-file
     */
    void mfile_output(std::ostream& out, ///< output stream
      int precision = 16,   ///< number of most significant digits, default=16
      int width = 0       ///< field width, default = precision+6
    ) const;

    //@}


    /** @name Input, Output (Friends)
     */
     //@{

     ///output format (nr and nc are ::Integer values, all others ::Real values): \n nr nc \\n A(1,1) A(1,2) ... A(1,nc) \\n A(2,1) ... A(nr,nc) \\n
    friend std::ostream& operator<<(std::ostream& o, const Symmatrix& A);

    ///input format (nr and nc are ::Integer values, all others ::Real values): \n nr nc \\n A(1,1) A(1,2) ... A(1,nc) \\n A(2,1) ... A(nr,nc) \\n
    friend std::istream& operator>>(std::istream& i, Symmatrix& A);

    //@}  


  };


  //@}

  // **************************************************************************
  //                make non inline friends available outside
  // **************************************************************************

    /// returns a column vector v consisting of the elements v(i)=(*this)(i,i), 0<=i<row dimension 
  Matrix diag(const Symmatrix& A);    //=(A(1,1),A(2,2),...)^t

  /// returns a symmetric diagonal matrix S of order A.dim() with vec(A) on the diagonal, i.e., S(i,i)=A(i) for all i and S(i,j)=0 for i!=j 
  Symmatrix Diag(const Matrix& A);    //vec(A) on the diagonal

  /// returns C=beta*C+alpha* A*A^T, where A may be transposed. If beta==0. then C is initiliazed to the correct size.
  Symmatrix& rankadd(const Matrix& A, Symmatrix& C,
    Real alpha = 1., Real beta = 0., int trans = 0);

  /// returns C=beta*C+alpha* A*D*A^T, where D is a vector representing a diagonal matrix and A may be transposed; if beta==0. then C is initialized to the correct size
  Symmatrix& scaledrankadd(const Matrix& A, const Matrix& D, Symmatrix& C,
    Real alpha = 1., Real beta = 0., int trans = 0);

  /// returns C=beta*C+alpha*(A*B^T+B*A^T)/2 [or for transposed (A^T*B+B^T*A)/2]. If beta==0. then C is initiliazed to the correct size.
  Symmatrix& rank2add(const Matrix& A, const Matrix& B, Symmatrix& C,
    Real alpha = 1., Real beta = 0., int trans = 0);

  ///returns C=beta*C+alpha*A*B, where B may be transposed; C must not be equal to B; if beta==0. then C is initialized to the correct size
  Matrix& genmult(const Symmatrix& A, const Matrix& B, Matrix& C,
    Real alpha = 1., Real beta = 0., int btrans = 0);

  ///returns C=beta*C+alpha*A*B, where A may be transposed; C must not be equal to A; if beta==0. then C is initialized to the correct size
  Matrix& genmult(const Matrix& A, const Symmatrix& B, Matrix& C,
    Real alpha = 1., Real beta = 0., int atrans = 0);

  ///returns C=beta*C+alpha*A*B, where B may be transposed; if beta==0. then C is initialized to the correct size
  Matrix& genmult(const Symmatrix& A, const Sparsemat& B, Matrix& C,
    Real alpha = 1., Real beta = 0., int btrans = 0);

  ///returns C=beta*C+alpha*A*B, where A may be transposed; if beta==0. then C is initialized to the correct size
  Matrix& genmult(const Sparsemat& A, const Symmatrix& B, Matrix& C,
    Real alpha = 1., Real beta = 0., int atrans = 0);

  /// returns C=beta*C+alpha* A*A^T, where A may be transposed; if beta==0. then C is initialized to the correct size
  Symmatrix& rankadd(const Sparsemat& A, Symmatrix& C,
    Real alpha = 1., Real beta = 0., int trans = 0);

  /// returns C=beta*C+alpha* A*D*A^T, where D is a vector representing a diagonal matrix and A may be transposed; if beta==0. then C is initialized to the correct size
  Symmatrix& scaledrankadd(const Sparsemat& A, const Matrix& D, Symmatrix& C,
    Real alpha = 1., Real beta = 0., int trans = 0);

  /// returns C=beta*C+alpha*(A*B^T+B*A^T)/2 [or for transposed (A^T*B+B^T*A)/2]. If beta==0. then C is initiliazed to the correct size.
  Symmatrix& rank2add(const Sparsemat& A, const Matrix& B, Symmatrix& C,
    Real alpha = 1., Real beta = 0., int trans = 0);

  /// returns a Symmatrix with elements abs(A(i,j))
  Symmatrix abs(const Symmatrix& A);

  /// returns the sum of the diagonal elements A(i,i) over all i
  Real trace(const Symmatrix& A);                  //=sum(diag(A))

  ///returns the usual inner product of A and B, i.e., the sum of A(i,j)*B(i,j) over all i,j 
  Real ip(const Symmatrix& A, const Symmatrix& B); //=trace(B^t*A)

  ///returns the usual inner product of A and B, i.e., the sum of A(i,j)*B(i,j) over all i,j 
  Real ip(const Matrix& A, const Symmatrix& B);    //=trace(B^t*A)

  ///returns the usual inner product of A and B, i.e., the sum of A(i,j)*B(i,j) over all i,j 
  Real ip(const Symmatrix& A, const Matrix& B);    //=trace(B^t*A)

  ///returns a row vector holding the sum over all rows, i.e., (1 1 ... 1)*A
  Matrix sumrows(const Symmatrix& A); //=(1 1 1 ... 1)*A

  ///returns a column vector holding the sum over all columns, i.e., A*(1 1 ... 1)^T
  Matrix sumcols(const Symmatrix& A); //=A*(1 1 ... 1)^t

  ///returns the sum over all elements of A, i.e., (1 1 ... 1)*A*(1 1 ... 1)^T
  Real sum(const Symmatrix& A);       //=(1 1 ... 1)*A*(1 1 ... 1)^t


  /// the symmetric vec operator stacks the lower triangle of A to a n*(n+1)/2 vector with the same norm2 as A; here it sets svec(A)=[a11,sqrt(2)a12,...,sqrt(2)a1n,a22,...,sqrt(2)a(n-1,n),ann]', multiplies it by a and sets or adds (if add==true) it to v starting from startindex_vec possibly restricted to the subblock of order blockdim (whenever >=0, else blockdim is set to A.rowdim()-startindex_A) starting from startindex_A (must be >=0); if add==false and startindex_vec<0 then vec is also reinitialzed to the appropriate size
  void svec(const Symmatrix& A, Matrix& sv, Real a = 1., bool add = false, Integer startindex_vec = -1, Integer startindex_A = 0, Integer blockdim = -1); //stores svec in sv

  /// the inverse operator to svec, extracts from v at startindex_vec (>=0) the symmetric matrix of blockdim adding its mutliple by a into A starting at startindex_A; if add==false and startindex_A<0 A is initialized to the size of blockdim; if the latter is also negative then v.dim()-startindex_vec must match an exact order and  matrix A is initialized to this size. In all other cases the size of the symmetric matrix determines the missing parameters and vec.dim-startindex_vec 
  void sveci(const Matrix& sv, Symmatrix& A, Real a = 1., bool add = false, Integer startindex_vec = 0, Integer startindex_A = -1, Integer blockdim = -1);

  /// def symmetric Kronecker product (A skron B)svec(C)=(BCA'+ACB')/2; sets S=alpha*(A skron B) or S*=... (if add==true) possibly shifted to the block starting at startindex_S;  if add==false and startindex_S<0, S is initialzed to the correct size 
  void skron(const Symmatrix& A, const Symmatrix& B, Symmatrix& S, Real alpha = 1., bool add = false, Integer startindex_S = -1);

  /// sets S=beta*S+alpha*B'*A*B for symmatrix A and matrix B
  void symscale(const Symmatrix& A, const Matrix& B, Symmatrix& S, Real alpha = 1., Real beta = 0., int btrans = 0);

  /// returns a row vector holding in each column the minimum over all rows in this column
  Matrix minrows(const Symmatrix& A); //minimum of each column (over rows)

  /// returns a column vector holding in each row the minimum over all columns in this row
  Matrix mincols(const Symmatrix& A); //minimum of each row (over colums)

  /// returns the minimum value over all elements of the matrix
  Real min(const Symmatrix& A);       //minimum of matrix

  /// returns a row vector holding in each column the maximum over all rows in this column
  Matrix maxrows(const Symmatrix& A); //similar

  /// returns a column vector holding in each row the maximum over all columns in this row
  Matrix maxcols(const Symmatrix& A); //similar

  /// returns the maximum value over all elements of the matrix
  Real max(const Symmatrix& A);       //similar

  ///output format (nr and nc are ::Integer values, all others ::Real values): \n nr nc \\n A(1,1) A(1,2) ... A(1,nc) \\n A(2,1) ... A(nr,nc) \\n
  std::ostream& operator<<(std::ostream& o, const Symmatrix& A);

  ///input format (nr and nc are ::Integer values, all others ::Real values): \n nr nc \\n A(1,1) A(1,2) ... A(1,nc) \\n A(2,1) ... A(nr,nc) \\n
  std::istream& operator>>(std::istream& i, Symmatrix& A);


  // **************************************************************************
  //                    implemenation of inline functions
  // **************************************************************************

  inline void Symmatrix::init_to_zero() {
    nr = 0; mem_dim = 0; m = 0;
#if (CONICBUNDLE_DEBUG>=1)
    is_init = true;
#endif
  }

  inline Symmatrix& Symmatrix::init(const Symmatrix& A, double d) {
    return xeya(A, d);
  }

  inline Symmatrix& Symmatrix::init(const Matrix& A, double d) {
    return xeya(A, d);
  }

  inline Symmatrix& Symmatrix::init(const Indexmatrix& A, double d) {
    return xeya(A, d);
  }

  inline Symmatrix& Symmatrix::init(Integer inr, Real d) {
    newsize(inr);
    mat_xea((nr * (nr + 1)) / 2, m, d);
    chk_set_init(*this, 1);
    return *this;
  }

  inline Symmatrix& Symmatrix::init(Integer inr, const Real* pd) {
    newsize(inr);
    mat_xey((nr * (nr + 1)) / 2, m, pd);
    chk_set_init(*this, 1);
    return *this;
  }

  inline Symmatrix::Symmatrix() {
    init_to_zero();
  }

  inline Symmatrix::Symmatrix(const Symmatrix& A, double d) :Memarrayuser() {
    init_to_zero(); xeya(A, d);
  }

  inline Symmatrix::Symmatrix(const Matrix& M, double d) {
    init_to_zero(); xeya(M, d);
  }

  inline Symmatrix::Symmatrix(const Indexmatrix& M, double d) {
    init_to_zero(); xeya(M, d);
  }

  inline Symmatrix::Symmatrix(Integer i) {
    init_to_zero(); newsize(i);
  }

  inline Symmatrix::Symmatrix(Integer i, Real d) {
    init_to_zero(); init(i, d);
  }

  inline Symmatrix::Symmatrix(Integer i, const Real* pd) {
    init_to_zero(); init(i, pd);
  }

  inline Symmatrix::~Symmatrix() {
    memarray->free(m);
  }

  inline Real& Symmatrix::operator()(Integer i, Integer j) {
    chk_range(i, j, nr, nr);
    //if (j>=i) return m[((i*((nr<<1)-(i+1)))>>1)+j];
    //return m[((j*((nr<<1)-(j+1)))>>1)+i];
    //now unsigned is used to avoid "cc1plus: warning: assuming signed overflow does not occur when assuming that (X + c) < X is always false [-Wstrict-overflow]"
    if (unsigned(j) >= unsigned(i)) return m[((unsigned(i) * ((unsigned(nr) << 1) - (unsigned(i) + 1))) >> 1) + unsigned(j)];
    return m[((unsigned(j) * ((unsigned(nr) << 1) - (unsigned(j) + 1))) >> 1) + unsigned(i)];
  }

  inline Real Symmatrix::operator()(Integer i, Integer j) const {
    chk_range(i, j, nr, nr);
    if (unsigned(j) >= unsigned(i)) return m[((unsigned(i) * ((unsigned(nr) << 1) - (unsigned(i) + 1))) >> 1) + unsigned(j)];
    return m[((unsigned(j) * ((unsigned(nr) << 1) - (unsigned(j) + 1))) >> 1) + unsigned(i)];
  }

  inline Real& Symmatrix::operator()(Integer i) {
    chk_range(i, 0, nr * nr, 1);
    // Integer j=i/nr;
    // i=i%nr;
    // if (j>=i) return m[((i*((nr<<1)-(i+1)))>>1)+j];
    // return m[((j*((nr<<1)-(j+1)))>>1)+i];
    const unsigned long uj = unsigned(i / nr);
    const unsigned long ui = unsigned(i % nr);
    if (uj >= ui) return m[((ui * ((unsigned(nr) << 1) - (ui + 1))) >> 1) + uj];
    return m[((uj * ((unsigned(nr) << 1) - (uj + 1))) >> 1) + ui];
  }

  inline Real Symmatrix::operator()(Integer i) const {
    chk_range(i, 0, nr * nr, 1);
    const unsigned long uj = unsigned(i / nr);
    const unsigned long ui = unsigned(i % nr);
    if (uj >= ui) return m[((ui * ((unsigned(nr) << 1) - (ui + 1))) >> 1) + uj];
    return m[((uj * ((unsigned(nr) << 1) - (uj + 1))) >> 1) + ui];
  }

  inline Symmatrix Symmatrix::principal_submatrix(const Indexmatrix& ind) const {
    Symmatrix S;
    principal_submatrix(ind, S);
    return S;
  }

  /// swap the content of the two matrices A and B (involves no copying)
  inline void swap(Symmatrix& A, Symmatrix& B) {
    Real* hm = A.m; A.m = B.m; B.m = hm;
    Integer hi = A.nr; A.nr = B.nr; B.nr = hi;
    hi = A.mem_dim; A.mem_dim = B.mem_dim; B.mem_dim = hi;
#if (CONICBUNDLE_DEBUG>=1)
    bool hb = A.is_init; A.is_init = B.is_init; B.is_init = hb;
#endif
  }

  ///returns x= alpha*y+beta*x; if beta==0. then x is initialized to the correct size
  inline Symmatrix& xbpeya(Symmatrix& x, const Symmatrix& y, Real alpha = 1., Real beta = 0.)
    //returns x= alpha*y+beta*x, where y may be transposed (ytrans=1)
    //if beta==0. then x is initialized to the correct size
  {
    if (beta == 0.) {
      chk_init(y);
      x.init(y, alpha);
    } else {
      chk_add(x, y);
      mat_xbpeya((x.rowdim() * (x.rowdim() + 1)) / 2,
        x.get_store(), y.get_store(), alpha, beta);
    }
    return x;
  }

  ///returns x= alpha*y+beta*z; x is initialized to the correct size
  inline Symmatrix& xeyapzb(Symmatrix& x, const Symmatrix& y, const Symmatrix& z, Real alpha = 1., Real beta = 1.)
    //returns x= alpha*y+beta*z
    //x is initialized to the correct size
  {
    chk_add(y, z);
    x.newsize(y.rowdim()); chk_set_init(x, 1);
    mat_xeyapzb((x.rowdim() * (x.rowdim() + 1)) / 2,
      x.get_store(), y.get_store(), z.get_store(), alpha, beta);
    return x;
  }

  inline Symmatrix& Symmatrix::operator=(const Symmatrix& A) {
    return xeya(A);
  }
  inline Symmatrix& Symmatrix::operator+=(const Symmatrix& A) {
    return xpeya(A);
  }
  inline Symmatrix& Symmatrix::operator-=(const Symmatrix& A) {
    return xpeya(A, -1.);
  }
  inline Symmatrix& Symmatrix::operator%=(const Symmatrix& A) {
    chk_add(*this, A); mat_xhadey((nr * (nr + 1)) / 2, m, A.m); return *this;
  }
  inline Symmatrix  Symmatrix::operator-() const {
    return Symmatrix(*this, -1.);
  }

  inline Symmatrix& Symmatrix::operator*=(Real d) {
    chk_init(*this); mat_xmultea(nr * (nr + 1) / 2, m, d); return *this;
  }
  inline Symmatrix& Symmatrix::operator/=(Real d) {
    chk_init(*this); mat_xmultea(nr * (nr + 1) / 2, m, 1. / d); return *this;
  }
  inline Symmatrix& Symmatrix::operator+=(Real d) {
    chk_init(*this); mat_xpea(nr * (nr + 1) / 2, m, d); return *this;
  }
  inline Symmatrix& Symmatrix::operator-=(Real d) {
    chk_init(*this); mat_xpea(nr * (nr + 1) / 2, m, -d); return *this;
  }

  /// returns a Matrix that equals A*B
  inline Matrix operator*(const Symmatrix& A, const Symmatrix& B) {
    Matrix C; return genmult(Matrix(A), B, C);
  }

  /// returns a Matrix that equals A%B (where % is overloaded as elementwise multiplication)
  inline Symmatrix operator%(const Symmatrix& A, const Symmatrix& B) {
    Symmatrix C(A); C %= B; return C;
  }

  /// returns a Matrix that equals A+B
  inline Symmatrix operator+(const Symmatrix& A, const Symmatrix& B) {
    Symmatrix C(A); return C.xpeya(B);
  }

  /// returns a Matrix that equals A-B
  inline Symmatrix operator-(const Symmatrix& A, const Symmatrix& B) {
    Symmatrix C(A); return C.xpeya(B, -1.);
  }

  /// returns a Matrix that equals A*B
  inline Matrix operator*(const Symmatrix& A, const Matrix& B) {
    Matrix C; return genmult(A, B, C);
  }

  /// returns a Matrix that equals A*B
  inline Matrix operator*(const Matrix& A, const Symmatrix& B) {
    Matrix C; return genmult(A, B, C);
  }

  /// returns a Matrix that equals A+B
  inline Matrix operator+(const Symmatrix& A, const Matrix& B) {
    Matrix C(A); return C.xpeya(B);
  }

  /// returns a Matrix that equals A+B
  inline Matrix operator+(const Matrix& A, const Symmatrix& B) {
    Matrix C(A); return C.xpeya(B);
  }

  /// returns a Matrix that equals A-B
  inline Matrix operator-(const Symmatrix& A, const Matrix& B) {
    Matrix C(A); return C.xpeya(B, -1);
  }

  /// returns a Matrix that equals A-B
  inline Matrix operator-(const Matrix& A, const Symmatrix& B) {
    Matrix C(A); return C.xpeya(B, -1);
  }

  /// returns a Symmatrix that equals A*d
  inline Symmatrix operator*(const Symmatrix& A, Real d) {
    return Symmatrix(A, d);
  }

  /// returns a Symmatrix that equals A*d
  inline Symmatrix operator*(Real d, const Symmatrix& A) {
    return Symmatrix(A, d);
  }

  /// returns a Symmatrix that equals A/d; ATTENTION: no check against division by zero
  inline Symmatrix operator/(const Symmatrix& A, Real d) {
    return Symmatrix(A, 1. / d);
  }

  /// returns a Symmatrix that equals A+d (d is added to each element)
  inline Symmatrix operator+(const Symmatrix& A, Real d) {
    Symmatrix B(A); return B += d;
  }

  /// returns a Symmatrix that equals A+d (d is added to each element)
  inline Symmatrix operator+(Real d, const Symmatrix& A) {
    Symmatrix B(A); return B += d;
  }

  /// returns a Symmatrix that equals A-d (d is subtracted from each element)
  inline Symmatrix operator-(const Symmatrix& A, Real d) {
    Symmatrix B(A); return B -= d;
  }

  /// returns a Symmatrix that equals d-A (each element subtracted from d)
  inline Symmatrix operator-(Real d, const Symmatrix& A) {
    Symmatrix B(A, -1.); return B += d;
  }


  /// the symmetric vec operator, stacks the lower triangle of A to a n*(n+1)/2 vector with the same norm2 as A; i.e., it returns svec(A)=[a11,sqrt(2)a12,...,sqrt(2)a1n,a22,...,sqrt(2)a(n-1,n),ann]'
  inline Matrix svec(const Symmatrix& A) //=(A(11),sqrt(2)A(12),...A(22)..)^t
  {
    Matrix sv; svec(A, sv); return sv;
  }

  /// the symmetric Kronecker product, defined via (A skron B)svec(C)=(BCA'+ACB')/2; sets or adds (if add==true) the symmetric matrix a*(A skron B) into S starting at startindex_S; if add==false and startindex_S<0, S is initialzed to the correct size 
  inline Symmatrix skron(const Symmatrix& A, const Symmatrix& B, Real alpha = 1., bool add = false, Integer startindex_S = -1) {
    Symmatrix S; skron(A, B, S, alpha, add, startindex_S); return S;
  }

  ///returns the Frobenius norm of A, i.e., the square root of the sum of A(i,j)*A(i,j) over all i,j 
  inline Real norm2(const Symmatrix& A) {
    return ::sqrt(ip(A, A));
  }

  /// returns a copy of A (drop it or use a constructor instead)
  inline Symmatrix transpose(const Symmatrix& A) {
    return Symmatrix(A);
  }


  inline Matrix::Matrix(const Symmatrix& A, Real d)

  {
    init_to_zero(); xeya(A, d);
  }

  inline Matrix& Matrix::init(const Symmatrix& A, Real d)

  {
    return xeya(A, d);
  }

  inline Matrix& Matrix::operator=(const Symmatrix& A)

  {
    return xeya(A);
  }

  inline Matrix& Matrix::operator*=(const Symmatrix& A)

  {
    Matrix C; return *this = genmult(*this, A, C);
  }

  inline Matrix& Matrix::operator+=(const Symmatrix& A)

  {
    return xpeya(A);
  }

  inline Matrix& Matrix::operator-=(const Symmatrix& A)

  {
    return xpeya(A, -1.);
  }


}

#ifndef CH_MATRIX_CLASSES__SPARSMAT_HXX
#include "sparsmat.hxx"
#endif

#endif

