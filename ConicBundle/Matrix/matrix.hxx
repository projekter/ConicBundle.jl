/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  Matrix/matrix.hxx
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


#ifndef CH_MATRIX_CLASSES__MATRIX_HXX
#define CH_MATRIX_CLASSES__MATRIX_HXX

/**  @file matrix.hxx
    @brief Header declaring the classes CH_Matrix_Classes::Realrange and  CH_Matrix_Classes::Matrix having Real elements
    @version 1.0
    @date 2005-03-01
    @author Christoph Helmberg

*/

#include <math.h>
#ifndef CH_MATRIX_CLASSES__INDEXMAT_HXX
#include "indexmat.hxx"
#endif


namespace CH_Matrix_Classes {

  //everything involving a "Sparsesym" is implemented in "sparssym.cxx/hxx"
  //everything else involving a "Sparsemat" is implemented in "sparsmat.cxx/hxx"
  //everything else involving a "Symmatrix" is implemented in "symmat.cxx/hxx"

  // **************************************************************************
  //                               Range
  // **************************************************************************

  /**@defgroup realrangegroup  Realrange (of real numbers with step size)
  */
  //@{

/// allows to specify a range of real values via (from, to, step,tol) meaning {x=from+i*step:x in(from-tol,to+tol),i in {0,1,2,...}}
  class Realrange {
  public:
    /// value of starting element
    Real from;
    /// upper bound on element values 
    Real to;
    /// step distance between successive elements starting with from
    Real step;
    /// tolerance for chekcing equality
    Real tol;
    /// constructor for {d: d=from+k*step for some k in N_0 and d<=to+tol}
    Realrange(Real _from, Real _to, Real _step = 1, Real _tol = 1e-8) :
      from(_from), to(_to), step(_step), tol(_tol) {
    }
    /// destructor, nothing to do
    ~Realrange() {
    }
  };

  //@}

// **************************************************************************
//                               Matrix
// **************************************************************************

/**@defgroup Matrixgroup Matrix (dense, real, m by n)
*/
//@{

/** @brief %Matrix class for real values of type #Real

    Internally a matrix of size nr x nc is stored in a one dimensional array of ::Real variables,
    the elements are arranged in columnwise order
    (a11,a21,...,anr1,a12,a22,...).

    Any matrix element can be indexed by (i,j), which internally refers to m[i+j*nr],
    or directly by the one dimensional index (i+j*nr). The latter view directly
    corresponds to the vec() operator often used in the linear algebra literature,
    i.e., the matrix is transformed to a vector by stacking the columns on top
    of each other.
 */
  class Matrix : protected Memarrayuser {
    friend class Indexmatrix;
    friend class Symmatrix;
    friend class Sparsemat;
    friend class Sparsesym;

  private:

    static const Mtype mtype;   ///< used for MatrixError templates (runtime type information was not yet existing)
    Integer mem_dim;   ///< amount of memory currently allocated
    Integer nr,        ///< number of rows
      nc;        ///< number of columns
    Real* m;           ///< pointer to store, order is columnwise (a11,a21,...,anr1,a12,a22,.....)

    //--- auxilliary storage space for repeated use of information, e.g., in qr
    /// the auxilliary storage aux_m may be used for different pruposes which are listed here; this is used to indicate the latest purpose in aux_task 
    enum AuxTask {
      none, ///< aux_m currently not in use 
      qr ///< aux_m  is now / was last  used by a QR factorization routine which stores the beta multipliers to be used in QR_mult type routines 
    };
    AuxTask aux_task;  ///< the current purpose of the auxilliary information
    Integer aux_dim;   ///< amount of memory currently allocated in auxm
    Real* aux_m;      ///< pointer to supplementary store for repeated use in some routines like qr

    bool is_init;      ///< flag whether memory is initialized, it is only used if CONICBUNDLE_DEBUG is defined

    /// initialize the matrix to a 0x0 matrix without storage
    inline void init_to_zero();

  public:

    //----------------------------------------------------
    //----  constructors, destructor, and initialization
    //----------------------------------------------------


    /** @name Constructors, Destructor, and Initialization (Members)
     */
     //@{

     /// empty matrix
    inline Matrix();
    /// copy constructor, *this=d*A
    inline Matrix(const Matrix&, Real d = 1., int atrans = 0);
    /// generate a column vector holding the elements of this Realrange
    inline Matrix(const Realrange&);
    /** @brief generate a matrix of size nr x nc but WITHOUT initializing the memory

        If initializing the memory externally and CONICBUNDLE_DEBUG is defined, please use
        set_init() via matrix.set_init(true) in order to avoid warnings concerning improper
        initialization
    */
    inline Matrix(Integer nr, Integer nc);
    /// generate a matrix of size nr x nc initializing all elements to the value d
    inline Matrix(Integer nr, Integer nc, Real d);
    /// generate a matrix of size nr x nc initializing the elements from the (one dimensional) array dp with increment incr and scaled by d
    inline Matrix(Integer nr, Integer nc, const Real* dp, Integer incr = 1, Real d = 1.);
    /// copy a std::vector<#Real> to a column vector of the same size and content
    inline Matrix(const std::vector<Real>& vec);
    ///
    ~Matrix() {
      memarray->free(m); if (aux_m) memarray->free(aux_m);
    }

#if (CONICBUNDLE_DEBUG>=1)
    /// after external initialization, call matrix.set_init(true) (not needed if CONICBUNDLE_DEBUG is undefined)
    void set_init(bool i) {
      is_init = i;
    }
    /// returns true if the matrix has been declared initialized (not needed if CONICBUNDLE_DEBUG is undefined)
    bool get_init() const {
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

    /// initialize to *this=A*d where A may be transposed
    inline Matrix& init(const Matrix& A, Real d = 1., int atrans = 0);
    /// initialize to *this=A*d  
    inline Matrix& init(const Indexmatrix& A, Real d = 1.);
    /// initialize to *this=A*d  
    inline Matrix& init(const Sparsemat& A, Real d = 1.); //*this=A*d;
    /// initialize to *this=A*d  
    inline Matrix& init(const Symmatrix& S, Real d = 1.);
    /// initialize to *this=A*d  
    inline Matrix& init(const Sparsesym&, Real d = 1.);
    /// initialize *this to a column vector holding the elements of Realrange
    Matrix& init(const Realrange&);
    /// intialize *this to a matrix of size nr x nc initializing all elements to the value d
    inline Matrix& init(Integer nr, Integer nc, Real d);
    /// generate a matrix of size nr x nc initializing the elements from the (one dimensional) array dp with increment incr and scaled by d
    inline Matrix& init(Integer nr, Integer nc, const Real* dp, Integer incr = 1, Real d = 1.);
    /// use std::vector<#Real> to initialize this to a column vector of the same size and content
    inline Matrix& init(const std::vector<Real>& vec);
    ///initialize to a diagonal nr x nr matrix with constant diagonal value d
    inline Matrix& init_diag(int nr, Real d = 1.);
    ///initialize to a diagonal matrix with diagonal given by vec
    inline Matrix& init_diag(const Matrix& vec, Real d = 1.);
    ///initialize to a diagonal matrix with diagonal given by vec
    inline Matrix& init_diag(const Indexmatrix& vec, Real d = 1.);

    /** @brief resize the matrix to nr x nc elements but WITHOUT initializing the memory

        If initializing the memory externally and CONICBUNDLE_DEBUG is defined, please use
        set_init() via matrix.set_init(true) in order to avoid warnings concerning improper
        initialization
    */
    void newsize(Integer nr, Integer nc);  //resize matrix without initialization


    //@}

    /** @name Conversions from other Matrix Classes (Members)
     */
     //@{

     /// (*this)=d*A
    inline Matrix(const Indexmatrix& A, Real d = 1.);
    /// (*this)=d*A
    inline Matrix(const Sparsemat& A, Real d = 1.);
    /// (*this)=d*A
    inline Matrix(const Symmatrix& S, Real d = 1.);
    /// (*this)=d*A
    inline Matrix(const Sparsesym&, Real d = 1.);

    //@}

    //----------------------------------------------------
    //----  size and type information
    //----------------------------------------------------

    /** @name Size and Type Information (Members)
     */
     //@{

     /// returns the number of rows in _nr and the number of columns in _nc
    void dim(Integer& _nr, Integer& _nc) const {
      _nr = nr; _nc = nc;
    }

    /// returns the dimension rows * columns when the matrix is regarded as a vector
    Integer dim() const {
      return nr * nc;
    }

    /// returns the row dimension
    Integer rowdim() const {
      return nr;
    }

    /// returns the column dimension
    Integer coldim() const {
      return nc;
    }

    /// returns the type of the matrix, MTmatrix
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
    inline Real& operator()(Integer i);      //index to vector of stacked columns

    /// returns value of element (i,j) of the matrix (rowindex i, columnindex j)
    inline Real operator()(Integer i, Integer j) const;

    /// returns value of element (i) of the matrix if regarded as vector of stacked columns [element (i%rowdim, i/rowdim)]
    inline Real operator()(Integer i) const; //index to vector of stacked columns

    /// returns a new submatrix as indexed by vecrow and veccol, A(i,j)=(*this)(vecrow(i),veccol(j)) for 0<=i<vecrow.dim(), 0<=j<veccol.dim() 
    Matrix operator()(const Indexmatrix& vecrow, const Indexmatrix& veccol) const;

    /// returns a new matrix B of the same shape as A with B(i,j)=(*this)(A(i),A(j)) for 0<=i<A.rowdim(), 0<=j<A.coldim() 
    Matrix operator()(const Indexmatrix& A) const;


    /// returns reference to element (i) of the matrix if regarded as vector of stacked columns [element (i%rowdim, i/rowdim)]
    inline Real& operator[](Integer i);      //{return (*this)(i);}
    /// returns value of element (i) of the matrix if regarded as vector of stacked columns [element (i%rowdim, i/rowdim)]
    inline Real operator[](Integer i) const; //{return (*this)(i);}

    /// returns column i copied to a new matrix 
    Matrix col(Integer i) const;          //returns column i as column vector
    /// returns row i copied to a new matrix 
    Matrix row(Integer i) const;          //returns row i as row vector
    /// returns a matrix of size this->rowdim() x vec.dim(), with column i a copy of column vec(i) of *this
    Matrix cols(const Indexmatrix& vec) const;  //returns cols as indexed by vec
    /// returns a matrix of size vec.dim() x this->coldim(), with row i a copy of row vec(i) of *this
    Matrix rows(const Indexmatrix& vec) const;  //returns rows as indexed by vec
    /// row i of this matrix is swapped with row j
    Matrix& swap_rowsij(Integer i, Integer j);
    /// column i of this matrix is swapped with column j
    Matrix& swap_colsij(Integer i, Integer j);
    /// for i=0 to rowdim row i of this matrix is swapped with row piv(i); for inverse =true the inverse permutation is generated
    Matrix& pivot_permute_rows(const Indexmatrix& piv, bool inverse = false);
    /// for j=0 to coldim column  j of this matrix is swapped with column piv(j); for inverse =true the inverse permutation is generated
    Matrix& pivot_permute_cols(const Indexmatrix& piv, bool inverse = false);

    /// keeps everything above and including diagonal d, everything below is set to zero, returns *this 
    Matrix& triu(Integer d = 0);    // (*this)(i,j)=0 for j<i+d
    /// keeps everything below and including diagonal d, everything above is set to zero, returns *this   
    Matrix& tril(Integer d = 0);    // (*this)(i,j)=0 for j>i+d

    /// assigns A to a submatrix of *this,  (*this)(vecrow(i),veccol(j))=A(i,j) for 0<=i<vecrow.dim(), 0<=j<veccol.dim()
    Matrix& subassign(const Indexmatrix& vecrow, const Indexmatrix& veccol,
      const Matrix& A);
    //(*this)(vecrow(i),veccol(j))=A(i,j) for all i,j
 /// assigns vector A to a subvector of *this,  (*this)(vec(i))=A(i) for 0<=i<vec.dim(), *this, vec, and A may be rectangular matrices, their dimesions are not changed, returns *this
    Matrix& subassign(const Indexmatrix& vec, const Matrix& A);
    //(*this)(vec(i))=A(i);
    //*this, vec, and A may be rect matrices but will be used as vectors

 /// all rows indexed by vector ind are deleted, no row should appear twice in ind, remaining rows are moved up keeping their order, returns *this
    Matrix& delete_rows(const Indexmatrix& ind, bool sorted_increasingly = false);
    /// all colmuns indexed by vector ind are deleted, no column should appear twice in ind, remaining columns are moved left keeping their order, returns *this
    Matrix& delete_cols(const Indexmatrix& ind, bool sorted_increasingly = false);

    /// insert the row vector v before row i, 0<=i<= row dimension, for i==row dimension the row is appended below; appending to a 0x0 matrix is allowed, returns *this
    Matrix& insert_row(Integer i, const Matrix& v); // 0<=i<=nr, v vector
    /// insert a column before column i, 0<=i<= column dimension, for i==column dimension the column is appended at the right; appending to a 0x0 matrix is allowed, returns *this
    Matrix& insert_col(Integer i, const Matrix& v); // 0<=i<=nc, v vector

    /// (*this) is set to a column vector of length min{max{0,n},dim()}; usually used to truncate a vector, returns *this
    inline Matrix& reduce_length(Integer n);
    //interpret *this as a column vector and cut at n 

 /// concats matrix A (or its tranpose) to the right of *this, A or *this may be the 0x0 matrix initally, returns *this
    Matrix& concat_right(const Matrix& A, int Atrans = 0);
    /// concats matrix A to the bottom of *this, A or *this may be the 0x0 matrix initally, returns *this
    Matrix& concat_below(const Matrix& A);
    /// concat value d at the bottom of *this, *this must be a column vector or the 0x0 matrix, returns *this 
    Matrix& concat_right(Real d);   //only for column vectors (nr==1)
    /// concat value d at the right of *this, *this must be a row vector or the 0x0 matrix, returns *this 
    Matrix& concat_below(Real d);   //only for row vectors (nc==1)
    /// enlarge the matrix by addnc>=0 columns without intializaton of the new columns, returns *this (marked as not initialized if nr>0)
    Matrix& enlarge_right(Integer addnc);
    /// enlarge the matrix by addnr>=0 rows without intializaton of the new rows, returns *this (marked as not initialized if addnr>0)
    Matrix& enlarge_below(Integer addnr);
    /// enlarge the matrix by addnc>=0 columns intializing the new columns by value d, returns *this 
    Matrix& enlarge_right(Integer addnc, Real d);
    /// enlarge the matrix by addnr>=0 rows intializing the new rows by value d, returns *this 
    Matrix& enlarge_below(Integer addnr, Real d);
    /// enlarge the matrix by addnc>=0 columns intializing the new columns by the values pointed to by dp times d, returns *this 
    Matrix& enlarge_right(Integer addnc, const Real* dp, Real d = 1.);
    /// enlarge the matrix by addnr>=0 rows intializing the new rows by the values pointed to by dp times d, returns *this 
    Matrix& enlarge_below(Integer addnr, const Real* dp, Real d = 1.);

    /// returns the current address of the internal value array; use cautiously, do not use delete!
    Real* get_store() {
      return m;
    }       //use cautiously, do not use delete!
/// returns the current address of the internal value array; use cautiously!
    const Real* get_store() const {
      return m;
    }   //use cautiously


//@}


/** @name Indexing and Submatrices (Friends)
 */
 //@{

 /// returns a column vector v consisting of the elements v(i)=(*this)(i,i), 0<=i<min(row dimension,column dimension) 
    friend Matrix diag(const Matrix& A);      //=(A(1,1),A(2,2),...)^t

    /// retuns a matrix that keeps the upper triangle of A starting with diagonal d, i.e., (i,j)=A(i,j) for 0<=i<row dimension, max(0,i+d)<=j<column dimension, and sets (i,j)=0 otherwise
    friend Matrix triu(const Matrix& A, Integer i);
    /// retuns a matrix that keeps the lower triangle of A starting with diagonal d, i.e., (i,j)=A(i,j) for 0<=i<row dimension, 0<=j<min(i+d+1,column dimension), and sets (i,j)=0 otherwise
    friend Matrix tril(const Matrix& A, Integer i);


    /// returns a new matrix [A, B], i.e., it concats matrices A and B rowwise; A or B may be a 0x0 matrix
    friend Matrix concat_right(const Matrix& A, const Matrix& B);
    /// returns a bew matrix [A; B], i.e., it concats matrices A and B columnwise; A or B may be a 0x0 matrix
    friend Matrix concat_below(const Matrix& A, const Matrix& B);

    /// swap the content of the two matrices A and B (involves no copying)
    friend void swap(Matrix& A, Matrix& B);

    //@}


    //------------------------------
    //----  BLAS-like Routines
    //------------------------------

    /** @name BLAS-like Routines (Members)
     */
     //@{

     ///sets *this=d*A where A may be transposed and returns *this
    Matrix& xeya(const Matrix& A, Real d = 1., int atrans = 0);   //*this=d*A
    ///sets *this+=d*A and returns *this  
    Matrix& xpeya(const Matrix& A, Real d = 1.);  //*this+=d*A;
    ///sets *this=d*A and returns *this
    Matrix& xeya(const Indexmatrix& A, Real d = 1.);   //*this=d*A
    ///sets *this+=d*A and returns *this  
    Matrix& xpeya(const Indexmatrix& A, Real d = 1.);  //*this+=d*A;

    //@}

    /** @name BLAS-like Routines (Friends)
     */
     //@{

     ///returns x= alpha*y+beta*x, where y may be transposed (ytrans=1); if beta==0. then x is initialized to the correct size
    friend Matrix& xbpeya(Matrix& x, const Matrix& y, Real alpha, Real beta, int ytrans);

    ///returns x= alpha*y+beta*z; x is initialized to the correct size
    friend Matrix& xeyapzb(Matrix& x, const Matrix& y, const Matrix& z, Real alpha, Real beta);

    ///returns C=beta*C+alpha*A*B, where A and B may be transposed; C must not be equal to A and B; if beta==0. then C is initialized to the correct size
    friend Matrix& genmult(const Matrix& A, const Matrix& B, Matrix& C,
      Real alpha, Real beta, int atrans, int btrans);
    //@}


    //------------------------------
    //----  usual operators
    //------------------------------

    /** @name Usual Arithmetic Operators (Members)
     */
     //@{

     ///
    inline Matrix& operator=(const Matrix& A);
    ///
    inline Matrix& operator*=(const Matrix& s);
    ///
    inline Matrix& operator+=(const Matrix& v);
    ///
    inline Matrix& operator-=(const Matrix& v);
    /// ATTENTION: this is redefined as the Hadamard product, (*this)(i,j)=(*this)(i,j)*A(i,j) for all i,j
    inline Matrix& operator%=(const Matrix& A);
    /// ATTENTION: this is redefined to act componentwise without checking for zeros, (*this)(i,j)=(*this)(i,j)/A(i,j) for all i,j
    inline Matrix& operator/=(const Matrix& A);
    ///
    inline Matrix operator-() const;

    /// 
    inline Matrix& operator*=(Real d);
    /// ATTENTION: d is NOT checked for 0
    inline Matrix& operator/=(Real d);
    /// sets (*this)(i,j)+=d for all i,j
    inline Matrix& operator+=(Real d);
    /// sets (*this)(i,j)-=d for all i,j
    inline Matrix& operator-=(Real d);

    ///transposes itself (cheap for vectors, expensive for matrices)
    Matrix& transpose();            //transposes itself

    //@}

    /** @name Usual Arithmetic Operators (Friends)
     */
     //@{

     ///
    friend Matrix operator*(const Matrix& A, const Matrix& B);
    ///
    friend Matrix operator+(const Matrix& A, const Matrix& B);
    ///
    friend Matrix operator-(const Matrix& A, const Matrix& B);
    /// ATTENTION: this is redefined as the Hadamard product, C(i,j)=A(i,j)*B(i,j) for all i,j
    friend Matrix operator%(const Matrix& A, const Matrix& B);
    /// ATTENTION: this is redefined to act componentwise without checking for zeros, C(i,j)=A(i,j)/B(i,j) for all i,j
    friend Matrix operator/(const Matrix& A, const Matrix& B);
    /// 
    friend Matrix operator*(const Matrix& A, Real d);
    /// 
    friend Matrix operator*(Real d, const Matrix& A);
    /// ATTENTION: d is NOT checked for 0
    friend Matrix operator/(const Matrix& A, Real d);
    /// returns (i,j)=A(i,j)+d for all i,j
    friend Matrix operator+(const Matrix& A, Real d);
    /// returns (i,j)=A(i,j)+d for all i,j
    friend Matrix operator+(Real d, const Matrix& A);
    /// returns (i,j)=A(i,j)-d for all i,j
    friend Matrix operator-(const Matrix& A, Real d);
    /// returns (i,j)=d-A(i,j) for all i,j
    friend Matrix operator-(Real d, const Matrix& A);

    ///
    friend Matrix transpose(const Matrix& A);

    //@}

    //------------------------------------------
    //----  Connections to other Matrix Classes
    //------------------------------------------

    /** @name Connections to other Classes (Members)
     */
     //@{


     ///sets *this=d*A and returns *this
    Matrix& xeya(const Symmatrix& A, Real d = 1.);  //*this=A*d;
    ///sets *this+=d*A and returns *this  
    Matrix& xpeya(const Symmatrix& A, Real d = 1.); //*this+=A*d;
    ///
    inline Matrix& operator=(const Symmatrix& S);
    ///
    inline Matrix& operator*=(const Symmatrix& S);
    ///
    inline Matrix& operator+=(const Symmatrix& S);
    ///
    inline Matrix& operator-=(const Symmatrix& S);

    ///sets *this=d*A and returns *this
    Matrix& xeya(const Sparsesym& A, Real d = 1.);  //*this=A*d;
    ///sets *this+=d*A and returns *this  
    Matrix& xpeya(const Sparsesym& A, Real d = 1.); //*this+=A*d;
    ///
    inline Matrix& operator=(const Sparsesym&);
    ///
    inline Matrix& operator*=(const Sparsesym& S);
    ///
    inline Matrix& operator+=(const Sparsesym& S);
    ///
    inline Matrix& operator-=(const Sparsesym& S);

    ///sets *this=d*A and returns *this
    Matrix& xeya(const Sparsemat& A, Real d = 1.);  //*this=A*d;
    ///sets *this+=d*A and returns *this  
    Matrix& xpeya(const Sparsemat& A, Real d = 1.); //*this+=A*d;
    ///
    inline Matrix& operator=(const Sparsemat& A);
    ///
    inline Matrix& operator*=(const Sparsemat& A);
    ///
    inline Matrix& operator+=(const Sparsemat& A);
    ///
    inline Matrix& operator-=(const Sparsemat& A);

    //@}

    /** @name Connections to other Classes (Friends)
     */
     //@{

     /// interpret A as a vector and copy it to a std::vector<double> which is also returned 
    friend std::vector<double>& assign(std::vector<double>& vec, const Matrix& A);


    ///returns C=beta*C+alpha*A*B, where A and B may be transposed; C must not be equal to A and B; if beta==0. then C is initialized to the correct size
    friend Matrix& genmult(const Symmatrix& A, const Matrix& B, Matrix& C,
      Real alpha, Real beta, int btrans);
    ///returns C=beta*C+alpha*A*B, where A and B may be transposed; C must not be equal to A and B; if beta==0. then C is initialized to the correct size
    friend Matrix& genmult(const Matrix& A, const Symmatrix& B, Matrix& C,
      Real alpha, Real beta, int atrans);

    ///returns C=beta*C+alpha*A*B, where A and B may be transposed; C must not be equal to A and B; if beta==0. then C is initialized to the correct size
    friend Matrix& genmult(const Sparsesym& A, const Matrix& B, Matrix& C,
      Real alpha, Real beta, int btrans);
    ///returns C=beta*C+alpha*A*B, where A and B may be transposed; C must not be equal to A and B; if beta==0. then C is initialized to the correct size
    friend Matrix& genmult(const Matrix& A, const Sparsesym& B, Matrix& C,
      Real alpha, Real beta, int atrans);

    ///returns C=beta*C+alpha*A*B, where A and B may be transposed; C must not be equal to A and B; if beta==0. then C is initialized to the correct size
    friend Matrix& genmult(const Sparsemat& A, const Matrix& B, Matrix& C,
      Real alpha, Real beta, int atrans, int btrans);
    ///returns C.col(colC)=beta*C.col(colC)+alpha*A*B.col(colB), where A and B may be transposed first; C must not be equal to A and B; if beta==0. then C is initialized, but the size of C must be correct already
    friend Matrix& genmult(const Sparsemat& A, const Matrix& B, Integer colB, Matrix& C, Integer colC, Real alpha, Real beta, int atrans, int btrans);
    ///returns C=beta*C+alpha*A*B, where A and B may be transposed; C must not be equal to A and B; if beta==0. then C is initialized to the correct size
    friend Matrix& genmult(const Matrix& A, const Sparsemat& B, Matrix& C,
      Real alpha, Real beta, int atrans, int btrans);

    //@}

    //------------------------------
    //----  Elementwise Operations
    //------------------------------

    /** @name Elementwise Operations (Members)
     */
     //@{

     /// resize *this to an nr x nc matrix and assign to (i,j) a random number uniformly from [0,1] for all i,j
    Matrix& rand(Integer nr, Integer nc, CH_Tools::GB_rand* random_generator = 0);
    /// resize *this to an nr x nc matrix and assign to (i,j) a random number from the normal distribution with given mean and variance (generators: 0 std, 1 mt, 2 mt64)
    Matrix& rand_normal(Integer nr, Integer nc, Real mean = 0., Real variance = 1., int generator_type = 0);
    /// shuffle the elements randomly (does not change dimensions)
    Matrix& shuffle(CH_Tools::GB_rand* random_generator = 0);
    /// sets (*this)(i,j)=1./(*this)(i,j) for all i,j and returns *this
    Matrix& inv(void);                        //reciprocal value componentwise
    /// sets (*this)(i,j)=sqrt((*this)(i,j)) for all i,j and returns *this
    Matrix& sqrt(void);
    /// sets (*this)(i,j)=sqr((*this)(i,j)) for all i,j and returns *this
    Matrix& sqr(void);
    /// sets (*this)(i,j)=sign((*this)(i,j),tol) for all i,j using ::sign(double,double) and returns *this
    Matrix& sign(Real tol = 1e-12);
    /// sets (*this)(i,j)=floor((*this)(i,j)) for all i,j and returns *this
    Matrix& floor(void);
    /// sets (*this)(i,j)=ceil((*this)(i,j)) for all i,j and returns *this
    Matrix& ceil(void);
    /// sets (*this)(i,j)=rint((*this)(i,j)) for all i,j and returns *this
    Matrix& rint(void);
    /// sets (*this)(i,j)=round((*this)(i,j)) for all i,j and returns *this
    Matrix& round(void);
    /// sets (*this)(i,j)=abs((*this)(i,j)) for all i,j and returns *this
    Matrix& abs(void);

    //@}

    /** @name Elementwise Operations (Friends)
     */
     //@{

     /// return a nr x nc matrix with (i,j) assigned a random number uniformly from [0,1] for all i,j
    friend Matrix rand(Integer nr, Integer nc, CH_Tools::GB_rand* random_generator);
    /// returns a matrix with elements (i,j)=1./((*this)(i,j)) for all i,j; ATTENTION: no check for division by zero 
    friend Matrix inv(const Matrix& A);
    /// returns a matrix with elements (i,j)=sqrt((*this)(i,j)) for all i,j 
    friend Matrix sqrt(const Matrix& A);
    /// returns a matrix with elements (i,j)=sqr((*this)(i,j)) for all i,j 
    friend Matrix sqr(const Matrix& A);
    /// returns a matrix with elements (i,j)=sign((*this)(i,j)) for all i,j using ::sign(double,double)
    friend Matrix sign(const Matrix& A, Real tol);
    /// returns a matrix with elements (i,j)=floor((*this)(i,j)) for all i,j 
    friend Matrix floor(const Matrix& A);
    /// returns a matrix with elements (i,j)=ceil((*this)(i,j)) for all i,j 
    friend Matrix ceil(const Matrix& A);
    /// returns a matrix with elements (i,j)=rint((*this)(i,j)) for all i,j 
    friend Matrix rint(const Matrix& A);
    /// returns a matrix with elements (i,j)=round((*this)(i,j)) for all i,j 
    friend Matrix round(const Matrix& A);
    /// returns a matrix with elements (i,j)=abs((*this)(i,j)) for all i,j 
    friend Matrix abs(const Matrix& A);

    //@}

    //----------------------------
    //----  Numerical Methods
    //----------------------------

    /** @name Numerical Methods (Members)
     */
     //@{
     ///returns true if any matrix element returns true on std::isnan
    bool contains_nan();

    ///scales each row i of (*this) by vec(i), i.e., (*this)=diag(vec)*(*this), and returns (*this)
    Matrix& scale_rows(const Matrix& vec);    //A=diag(vec)*A
    ///scales each column i of (*this) by vec(i), i.e., (*this)=(*this)*diag(vec), and returns (*this)
    Matrix& scale_cols(const Matrix& vec);    //A=A*diag(vec)

    //----- triangular routines
    /// solves (*this)*x=rhs for x by back substitution regarding (*this) as an upper triangle matrix and stores x in rhs. Returns 0 on success, otherwise i+1 if abs(*this)(i,i)<tol and the remaining row of rhs is nonzero.
    int triu_solve(Matrix& rhs, Real tol = 1e-10);   //only elements i>=j are used
    /// solves (*this)*x=rhs for x by forward substitution regarding (*this) as an upper triangle matrix and stores x in rhs. Returns 0 on success, otherwise i+1 if abs(*this)(i,i)<tol and the reduced row of rhs is nonzero.
    int tril_solve(Matrix& rhs, Real tol = 1e-10);   //only elements i<=j are used


    //----- QR-Factorization
    /// computes a Householder QR_factorization overwriting (*this); returns 0 on success, otherwise column index +1 if the norm of the column is below tol when reaching it. 
    int QR_factor(Real tol = 1e-10);                            //factorization stored in *this
    /// computes a Householder QR_factorization computing the Q matrix explicitly and setting (*this)=R; it always returns 0 
    int QR_factor(Matrix& Q, Real tol = 1e-10);                   //afterwards *this is R
    /// computes a Householder QR_factorization computing matrices Q and R explicitly and leaving (*this) unchanged; it always returns 0 
    int QR_factor(Matrix& Q, Matrix& R, Real tol) const;    //*this is unchanged

    /// computes a Householder QR_factorization with pivoting  and overwriting (*this); the pivoting permutation is stored in piv; returns the rank 
    int QR_factor(Indexmatrix& piv, Real tol = 1e-10);
    /// computes a Householder QR_factorization with pivoting on those with tolerable norm and tolerable reduction in norm and overwriting (*this); the pivoting permutation is stored in piv; returns the rank 
    int QR_factor_relpiv(Indexmatrix& piv, Real tol = 1e-10);
    /// Computes a Householder QR_factorization with pivoting computing the Q matrix explicitly and setting (*this)=R; the pivoting permutation is stored in piv; returns the rank 
    int QR_factor(Matrix& Q, Indexmatrix& piv, Real tol = 1e-10);
    /// computes a Householder QR_factorization with pivoting computing matrices Q and R explicitly and leaving (*this) unchanged; the pivoting permutation is stored in piv; returns the rank 
    int QR_factor(Matrix& Q, Matrix& R, Indexmatrix& piv, Real tol) const;
    /// computes A=transpose(Q)*A, assuming a housholder Q is coded in the first r columns of the lower triangle of (*this); it always returns 0 
    int Qt_times(Matrix& A, Integer r) const;
    /// computes A=Q*A, assuming a housholder Q is coded in the first r columns of the lower triangle of (*this); it always returns 0 
    int Q_times(Matrix& A, Integer r) const;
    /// computes A=A*Q, assuming a housholder Q is coded in the first r columns of the lower triangle of (*this); it always returns 0 
    int times_Q(Matrix& A, Integer r) const;
    /** @brief solves (*this)*x=rhs by factorizing and overwriting (*this); rhs is overwritten with the solution.  Returns 0 on success, otherwise i+1 if in the backsolve abs(*this)(i,i)<tol and the reduced row of rhs is nonzero.

        To avoid overwriting (*this), use the appropriate version
        of #QR_factor and #triu_solve.
     */
    int QR_solve(Matrix& rhs, Real tol = 1e-10);

    /** extend the current Householder QR-factorization stored in this by appending
        and factorizing the columns of A yielding a new QR_factorization

        @param[in] A contains the addtionial columns to be factorized

        @param[in,out] piv contains on input, the permution vector returned by QR_factor
            with pivoting, and on output the new entire permutation

        @param[in] r gives the initial rank of *this as returned by QR_factor with pviaton

        @param[in] tol gives the tolerance for regarding a vector as having norm zero

        @return the rank of the new QR-facotrization
    */
    int QR_concat_right(const Matrix& A, Indexmatrix& piv, Integer r, Real tol = 1e-10);


    //----- Least squares  
    ///computes a least squares solution by #QR_solve, overwriting (*this). rhs is overwritten with the solution. In fact, the full code is return this->QR_solve(rhs,tol);
    int ls(Matrix& rhs, Real tol);

    /** @brief computes a nonnegative least squares solution; rhs is overwritten by the solution; if dual!=0, the dual variables are stored there; returns 0 on success, 1 on failure

      Computes the least squares solution of min ||Ax-b|| s.t. x >=0;\n
      The KKT system A'*A*x - A'*b - l = 0; x >=0, l>=0, x'*l=0 is solved
      solved by interior point method with QR-solution of the extended system.

      The current implementation is based on
      [P. Matsoms, "Sparse Linear Least Squares Problems in Optimization",
      Comput. Opt. and Appl., 7, 89-110 (1997)] but is only
      a quick and rather sloppy implementation of it ...
    */
    int nnls(Matrix& rhs, Matrix* dual = 0, Real tol = 1e-10) const;

    //@}

    /** @name Numerical Methods (Friends)
     */
     //@{

     /// returns the sum of the diagonal elements A(i,i) over all i
    friend Real trace(const Matrix& A);               //=sum(diag(A))
    /// returns the usual inner product of A and B, i.e., the sum of A(i,j)*B(i,j) over all i,j 
    friend Real ip(const Matrix& A, const Matrix& B); //=trace(B^t*A)
    /// returns in addition to the usual inner product of A and B the minimum and maximum value of A(i,j)*B(i,j) over all (i,j)
    friend Real ip_min_max(const Matrix& A, const Matrix& B, Real& minval, Real& maxval); //=trace(B^t*A), minval=min(B_ij*A_ij), maxval=max(B_ij*A_ij)
    ///returns the squared Frobenius norm of column j of A, i.e., the sum of A(i,j)*A(i,j) over all i with possibly (if scaling!=0) each term i multiplied by (*scaling)(i)
    friend Real colip(const Matrix& A, Integer j, const Matrix* scaling);
    ///returns the squared Frobenius norm of row i of A, i.e., the sum of A(i,j)*A(i,j) over all j  with possibly (if scaling!=0) each term j multiplied by (*scaling)(j)
    friend Real rowip(const Matrix& A, Integer i, const Matrix* scaling);
    ///returns the column vector of the squared Frobenius norm of all columns j of A, i.e., the sum of A(i,j)*A(i,j) over all i for each j
    friend Matrix colsip(const Matrix& A);
    ///returns the row vector of the squared Frobenius norm of all rowd i of A, i.e., the sum of A(i,j)*A(i,j) over all j for each i
    friend Matrix rowsip(const Matrix& A);
    ///returns the Frobenius norm of A, i.e., the square root of the sum of A(i,j)*A(i,j) over all i,j 
    friend Real norm2(const Matrix& A);               //=sqrt(ip(A,A));
    ///returns trace(A^TDA)=\|A\|^2_D with D=Diag(d). A may be transposed, D may be inverted but there is no check for division by zero
    friend Real normDsquared(const Matrix& A, const Matrix& d, int atrans, int dinv);

    ///returns a row vector holding the sum over all rows, i.e., (1 1 ... 1)*A
    friend Matrix sumrows(const Matrix& A);   //=(1 1 1 ... 1)*A
    ///returns a column vector holding the sum over all columns, i.e., A*(1 1 ... 1)^T
    friend Matrix sumcols(const Matrix& A);   //=A*(1 1 ... 1)^t
    ///returns the sum over all elements of A, i.e., (1 1 ... 1)*A*(1 1 ... 1)^T
    friend Real sum(const Matrix& A);         //=(1 1 ... 1)*A*(1 1 ... 1)^t

    //----- Householder rotations
    /// returns the Householder vector of size A.rowdim() for the subcolumn A(i:A.rowdim(),j) 
    friend Matrix house(const Matrix& A, Integer i, Integer j, Real tol);
    /// Housholder pre-multiplication of A with Householder vector v; the first nonzero of v is index i, the multplication is applied to all columns of A with index >=j; always returns 0
    friend int rowhouse(Matrix& A, const Matrix& v, Integer i, Integer j);
    /// Housholder post-multiplication of A with Householder vector v; the first nonzero of v is index i, the multplication is applied to all rows of A with index >=j; always returns 0
    friend int colhouse(Matrix& A, const Matrix& v, Integer i, Integer j);

    /// computes a Householder QR factorization of A and outputs Q and R leaving A unchanged; always returns 0  
    friend int QR_factor(const Matrix& A, Matrix& Q, Matrix& R, Real tol);
    /// computes a Householder QR factorization of A with pivating. It outputs Q, R, and the pivoting permuation in piv; returns the rank of A  
    friend int QR_factor(const Matrix& A, Matrix& Q, Matrix& R, Indexmatrix& piv, Real tol);

    //@}    

    //---------------------------------------------
    //----  Comparisons / Max / Min / sort / find
    //---------------------------------------------

    /** @name Find (Members)
     */
     //@{

     /// returns an Indexmatrix ind so that (*this)(ind(i)) 0<=i<ind.dim() runs through all nonzero elements  
    Indexmatrix find(Real tol = 1e-10) const;   //finds nonzeros
    /// returns an Indexmatrix ind so that (*this)(ind(i)) 0<=i<ind.dim() runs through all elements of value num 
    Indexmatrix find_number(Real num = 0., Real tol = 1e-10) const;

    //@}

    /** @name Comparisons, Max, Min, Sort, Find (Friends)
     */
     //@{

     /// returns a matrix having elements (i,j)=Real(A(i,j)<B(i,j)) for all i,j
    friend Matrix operator<(const Matrix& A, const Matrix& B);
    /// returns a matrix having elements (i,j)=Real(A(i,j)>B(i,j)) for all i,j
    friend Matrix operator>(const Matrix& A, const Matrix& B);
    /// returns a matrix having elements (i,j)=Real(A(i,j)<=B(i,j)) for all i,j
    friend Matrix operator<=(const Matrix& A, const Matrix& B);
    /// returns a matrix having elements (i,j)=Real(A(i,j)>=B(i,j)) for all i,j
    friend Matrix operator>=(const Matrix& A, const Matrix& B);
    /// returns a matrix having elements (i,j)=Real(A(i,j)==B(i,j)) for all i,j
    friend Matrix operator==(const Matrix& A, const Matrix& B);
    /// returns a matrix having elements (i,j)=Real(A(i,j)!=B(i,j)) for all i,j
    friend Matrix operator!=(const Matrix& A, const Matrix& B);
    /// returns a matrix having elements (i,j)=Real(A(i,j)<d) for all i,j
    friend Matrix operator<(const Matrix& A, Real d);
    /// returns a matrix having elements (i,j)=Real(A(i,j)>d) for all i,j
    friend Matrix operator>(const Matrix& A, Real d);
    /// returns a matrix having elements (i,j)=Real(A(i,j)<=d) for all i,j
    friend Matrix operator<=(const Matrix& A, Real d);
    /// returns a matrix having elements (i,j)=Real(A(i,j)>=d) for all i,j
    friend Matrix operator>=(const Matrix& A, Real d);
    /// returns a matrix having elements (i,j)=Real(A(i,j)==d) for all i,j
    friend Matrix operator==(const Matrix& A, Real d);
    /// returns a matrix having elements (i,j)=Real(A(i,j)!=d) for all i,j
    friend Matrix operator!=(const Matrix& A, Real d);
    /// returns a matrix having elements (i,j)=Real(d<A(i,j)) for all i,j
    friend Matrix operator<(Real d, const Matrix& A);
    /// returns a matrix having elements (i,j)=Real(d>A(i,j)) for all i,j
    friend Matrix operator>(Real d, const Matrix& A);
    /// returns a matrix having elements (i,j)=Real(d<=A(i,j)) for all i,j
    friend Matrix operator<=(Real d, const Matrix& A);
    /// returns a matrix having elements (i,j)=Real(d>=A(i,j)) for all i,j
    friend Matrix operator>=(Real d, const Matrix& A);
    /// returns a matrix having elements (i,j)=Real(d==A(i,j)) for all i,j
    friend Matrix operator==(Real d, const Matrix& A);
    /// returns a matrix having elements (i,j)=Real(d!=A(i,j)) for all i,j
    friend Matrix operator!=(Real d, const Matrix& A);

    /// returns true if both matrices have the same size and the same elements    
    friend bool equal(const Matrix& A, const Matrix& B);

    /// returns a row vector holding in each column the minimum over all rows in this column
    friend Matrix minrows(const Matrix& A);   //min of each column (over the rows)
    /// returns a column vector holding in each row the minimum over all columns in this row
    friend Matrix mincols(const Matrix& A);   //min of each row (over the columns)
    /// returns the minimum value over all elements of the matrix
    friend Real min(const Matrix& A, Integer* iindex, Integer* jindex);
    /// returns a row vector holding in each column the maximum over all rows in this column
    friend Matrix maxrows(const Matrix& A);   //similar
    /// returns a column vector holding in each row the maximum over all columns in this row
    friend Matrix maxcols(const Matrix& A);   //similar
    /// returns the maximum value over all elements of the matrix
    friend Real max(const Matrix& A, Integer* iindex, Integer* jindex);

    /// returns an Indexmatrix ind so that vec(ind(0))<=vec(ind(1))<=...<=vec(ind(vec.dim()-1)) (vec may be rectangular, set nondecreasing=false for opposite order)
    friend Indexmatrix sortindex(const Matrix& vec, bool nondecreasing);
    /// sets ind so that vec(ind(0))<=vec(ind(1))<=...<=vec(ind(vec.dim()-1)) (vec may be rectangular, set nondecreasing=false for opposite order)
    friend void sortindex(const Matrix& vec, Indexmatrix& ind, bool nondecreasing);

    /// returns an Indexmatrix ind so that A(ind(i)) 0<=i<ind.dim() runs through all nonzero elements with abs(A(j))>tol 
    friend Indexmatrix find(const Matrix& A, Real tol);
    /// returns an Indexmatrix ind so that A(ind(i)) 0<=i<ind.dim() runs through all elements of A having value num, i.e., abs(A(j)-num)<tol 
    friend Indexmatrix find_number(const Matrix& A, Real num, Real tol);

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
    friend std::ostream& operator<<(std::ostream& o, const Matrix& v);
    ///input format (nr and nc are ::Integer values, all others ::Real values): \n nr nc \\n A(1,1) A(1,2) ... A(1,nc) \\n A(2,1) ... A(nr,nc) \\n
    friend std::istream& operator>>(std::istream& i, Matrix& v);

    //@}  

  };

  //@}

  // **************************************************************************
  //                make non inline friends available outside
  // **************************************************************************

  /// returns a column vector v consisting of the elements v(i)=(*this)(i,i), 0<=i<min(row dimension,column dimension)   
  Matrix diag(const Matrix& A);      //=(A(1,1),A(2,2),...)^t

  ///returns x= alpha*y+beta*x, where y may be transposed (ytrans=1); if beta==0. then x is initialized to the correct size
  Matrix& xbpeya(Matrix& x, const Matrix& y, Real alpha = 1., Real beta = 0., int ytrans = 0);

  ///returns x= alpha*y+beta*z; x is initialized to the correct size
  Matrix& xeyapzb(Matrix& x, const Matrix& y, const Matrix& z, Real alpha = 1., Real beta = 1.);

  ///returns C=beta*C+alpha*A*B, where A and B may be transposed; C must not be equal to A and B; if beta==0. then C is initialized to the correct size
  Matrix& genmult(const Matrix& A, const Matrix& B, Matrix& C,
    Real alpha = 1., Real beta = 0., int atrans = 0, int btrans = 0);

  /// returns a transposed matrix of A
  Matrix transpose(const Matrix& A);

  /// interpret A as a vector and copy it to a std::vector<double> which is also returned 
  std::vector<double>& assign(std::vector<double>& vec, const Matrix& A);

  ///returns C=beta*C+alpha*A*B, where A and B may be transposed; C must not be equal to A and B; if beta==0. then C is initialized to the correct size
  Matrix& genmult(const Symmatrix& A, const Matrix& B, Matrix& C,
    Real alpha, Real beta, int btrans);

  ///returns C=beta*C+alpha*A*B, where A and B may be transposed; C must not be equal to A and B; if beta==0. then C is initialized to the correct size
  Matrix& genmult(const Matrix& A, const Symmatrix& B, Matrix& C,
    Real alpha, Real beta, int atrans);

  ///returns C=beta*C+alpha*A*B, where A and B may be transposed; C must not be equal to A and B; if beta==0. then C is initialized to the correct size
  Matrix& genmult(const Sparsesym& A, const Matrix& B, Matrix& C,
    Real alpha, Real beta, int btrans);

  ///returns C=beta*C+alpha*A*B, where A and B may be transposed; C must not be equal to A and B; if beta==0. then C is initialized to the correct size
  Matrix& genmult(const Matrix& A, const Sparsesym& B, Matrix& C,
    Real alpha, Real beta, int atrans);

  ///returns C=beta*C+alpha*A*B, where A and B may be transposed; C must not be equal to A and B; if beta==0. then C is initialized to the correct size
  Matrix& genmult(const Sparsemat& A, const Matrix& B, Matrix& C,
    Real alpha, Real beta, int atrans, int btrans);

  ///returns C=beta*C+alpha*A*B, where A and B may be transposed; C must not be equal to A and B; if beta==0. then C is initialized to the correct size
  Matrix& genmult(const Matrix& A, const Sparsemat& B, Matrix& C,
    Real alpha, Real beta, int atrans, int btrans);

  /// returns a matrix with elements (i,j)=abs((*this)(i,j)) for all i,j 
  Matrix abs(const Matrix& A);

  ///=sum(diag(A))  
  Real trace(const Matrix& A);

  ///returns trace(A^TDA)=\|A\|^2_D with D=Diag(d). A may be transposed, D may be inverted but there is no check for division by zero
  Real normDsquared(const Matrix& A, const Matrix& d, int atrans = 0, int dinv = 0);

  ///returns a row vector holding the sum over all rows, i.e., (1 1 ... 1)*A
  Matrix sumrows(const Matrix& A);   //=(1 1 1 ... 1)*A

  ///returns a column vector holding the sum over all columns, i.e., A*(1 1 ... 1)^T
  Matrix sumcols(const Matrix& A);   //=A*(1 1 ... 1)^t

  ///returns the sum over all elements of A, i.e., (1 1 ... 1)*A*(1 1 ... 1)^T
  Real sum(const Matrix& A);         //=(1 1 ... 1)*A*(1 1 ... 1)^t

  /// returns the Householder vector of size A.rowdim() for the subcolumn A(i:A.rowdim(),j) 
  Matrix house(const Matrix& x, Integer i = 0, Integer j = 0, Real tol = 1e-10);

  /// Housholder pre-multiplication of A with Householder vector v; the first nonzero of v is index i, the multplication is applied to all columns of A with index >=j; always returns 0
  int rowhouse(Matrix& A, const Matrix& v, Integer i = 0, Integer j = 0);

  /// Housholder post-multiplication of A with Householder vector v; the first nonzero of v is index i, the multplication is applied to all rows of A with index >=j; always returns 0
  int colhouse(Matrix& A, const Matrix& v, Integer i = 0, Integer j = 0);

  /// returns a matrix having elements (i,j)=Real(A(i,j)<B(i,j)) for all i,j
  Matrix operator<(const Matrix& A, const Matrix& B);

  /// returns a matrix having elements (i,j)=Real(A(i,j)<=B(i,j)) for all i,j
  Matrix operator<=(const Matrix& A, const Matrix& B);

  /// returns a matrix having elements (i,j)=Real(A(i,j)==B(i,j)) for all i,j
  Matrix operator==(const Matrix& A, const Matrix& B);

  /// returns a matrix having elements (i,j)=Real(A(i,j)!=B(i,j)) for all i,j
  Matrix operator!=(const Matrix& A, const Matrix& B);

  /// returns a matrix having elements (i,j)=Real(A(i,j)<d) for all i,j
  Matrix operator<(const Matrix& A, Real d);

  /// returns a matrix having elements (i,j)=Real(A(i,j)>d) for all i,j
  Matrix operator>(const Matrix& A, Real d);

  /// returns a matrix having elements (i,j)=Real(A(i,j)<=d) for all i,j
  Matrix operator<=(const Matrix& A, Real d);

  /// returns a matrix having elements (i,j)=Real(A(i,j)>=d) for all i,j
  Matrix operator>=(const Matrix& A, Real d);

  /// returns a matrix having elements (i,j)=Real(A(i,j)==d) for all i,j
  Matrix operator==(const Matrix& A, Real d);

  /// returns a matrix having elements (i,j)=Real(A(i,j)!=d) for all i,j
  Matrix operator!=(const Matrix& A, Real d);

  /// returns a row vector holding in each column the minimum over all rows in this column
  Matrix minrows(const Matrix& A);   //min of each column (over the rows)

  /// returns a column vector holding in each row the minimum over all columns in this row
  Matrix mincols(const Matrix& A);   //min of each row (over the columns)

  /// returns the minimum value over all elements of the matrix
  Real min(const Matrix& A, Integer* iindex = 0, Integer* jindex = 0);

  /// returns a row vector holding in each column the maximum over all rows in this column
  Matrix maxrows(const Matrix& A);   //similar

  /// returns a column vector holding in each row the maximum over all columns in this row
  Matrix maxcols(const Matrix& A);   //similar

  /// returns the maximum value over all elements of the matrix
  Real max(const Matrix& A, Integer* iindex = 0, Integer* jindex = 0);

  /// returns an Indexmatrix ind so that vec(ind(0))<=vec(ind(1))<=...<=vec(ind(vec.dim()-1)) (vec may be rectangular, set nondecreasing=false for opposite order)
  Indexmatrix sortindex(const Matrix& vec, bool nondecreasing = true);

  /// sets ind so that vec(ind(0))<=vec(ind(1))<=...<=vec(ind(vec.dim()-1)) (vec may be rectangular, set nondecreasing=false for opposite order)
  void sortindex(const Matrix& vec, Indexmatrix& ind, bool nondecreasing = true);

  ///output format (nr and nc are ::Integer values, all others ::Real values): \n nr nc \\n A(1,1) A(1,2) ... A(1,nc) \\n A(2,1) ... A(nr,nc) \\n
  std::ostream& operator<<(std::ostream& o, const Matrix& v);

  ///input format (nr and nc are ::Integer values, all others ::Real values): \n nr nc \\n A(1,1) A(1,2) ... A(1,nc) \\n A(2,1) ... A(nr,nc) \\n
  std::istream& operator>>(std::istream& i, Matrix& v);

  // **************************************************************************
  //                   implementation of inline functions
  // **************************************************************************

  inline void Matrix::init_to_zero() {
    nr = nc = 0; mem_dim = 0; m = 0;
    aux_task = none;
    aux_dim = 0;
    aux_m = 0;
    chk_set_init(*this, 1);
  }

  inline Matrix& Matrix::init(Integer inr, Integer inc, Real d) {
    newsize(inr, inc);
    mat_xea(nr * nc, m, d);
    chk_set_init(*this, 1);
    return *this;
  }

  inline Matrix& Matrix::init(Integer inr, Integer inc, const Real* p, Integer incr, Real d) {
    newsize(inr, inc);
    if (d == 1.) {
      if (incr == 1) mat_xey(nr * nc, m, p);
      else mat_xey(nr * nc, m, Integer(1), p, incr);
    } else {
      if (incr == 1) mat_xeya(nr * nc, m, p, d);
      else mat_xeya(nr * nc, m, Integer(1), p, incr, d);
    }
    chk_set_init(*this, 1);
    return *this;
  }

  inline Matrix& Matrix::init(const Matrix& A, Real d, int atrans) {
    return xeya(A, d, atrans);
  }

  inline Matrix& Matrix::init(const Indexmatrix& A, Real d) {
    return xeya(A, d);
  }

  inline Matrix& Matrix::init(const std::vector<Real>& vec) {
    newsize(Integer(vec.size()), 1); chk_set_init(*this, 1);
    for (Integer i = 0; i < nr; i++) m[i] = vec[(unsigned long)(i)];
    return *this;
  }

  inline Matrix& Matrix::init_diag(int n, Real d) {
    init(n, n, 0.);
    mat_xea(n, m, n + 1, d);
    return *this;
  }

  inline Matrix& Matrix::init_diag(const Matrix& vec, Real d) {
    Integer n = vec.dim();
    init(n, n, 0.);
    if (d == 1.)
      mat_xey(n, m, n + 1, vec.get_store(), 1);
    else
      mat_xeya(n, m, n + 1, vec.get_store(), 1, d);
    return *this;
  }

  inline Matrix& Matrix::init_diag(const Indexmatrix& vec, Real d) {
    Integer n = vec.dim();
    init(n, n, 0.);
    if (d == 1.)
      for (Integer i = 0; i < n; i++) m[i * (n + 1)] = Real(vec.m[i]);
    else
      for (Integer i = 0; i < n; i++) m[i * (n + 1)] = vec.m[i] * d;
    return *this;
  }


  inline Matrix::Matrix() :Memarrayuser() {
    init_to_zero();
  }

  inline Matrix::Matrix(const Matrix& A, Real d, int atrans) : Memarrayuser() {
    init_to_zero();
    xeya(A, d, atrans);
  }

  inline Matrix::Matrix(const Indexmatrix& A, Real d) : Memarrayuser() {
    init_to_zero();
    xeya(A, d);
  }

  inline Matrix::Matrix(Integer inr, Integer inc) : Memarrayuser() {
    init_to_zero();
    newsize(inr, inc);
  }

  inline Matrix::Matrix(const Realrange& range) : Memarrayuser() {
    init_to_zero();
    init(range);
  }

  inline Matrix::Matrix(Integer inr, Integer inc, Real d) : Memarrayuser() {
    init_to_zero();
    init(inr, inc, d);
  }

  inline Matrix::Matrix(Integer inr, Integer inc, const Real* p, Integer incr, Real d) : Memarrayuser() {
    init_to_zero();
    init(inr, inc, p, incr, d);
  }

  inline Matrix::Matrix(const std::vector<Real>& vec) : Memarrayuser() {
    init_to_zero();
    init(vec);
  }

  inline Real& Matrix::operator()(Integer i, Integer j) {
    chk_range(i, j, nr, nc);
    return m[j * nr + i];
  }

  inline Real& Matrix::operator()(Integer i) {
    chk_range(i, 0, nr * nc, 1);
    return m[i];
  }

  inline Real Matrix::operator()(Integer i, Integer j) const {
    chk_range(i, j, nr, nc);
    return m[j * nr + i];
  }

  inline Real Matrix::operator()(Integer i) const {
    chk_range(i, 0, nr * nc, 1);
    return m[i];
  }

  inline Real& Matrix::operator[](Integer i) {
    return (*this)(i);
  }

  inline Real Matrix::operator[](Integer i) const {
    return (*this)(i);
  }

  inline Matrix& Matrix::reduce_length(Integer n) {
    nr = min(nr * nc, max(Integer(0), n)); nc = 1; return *this;
  }



  /// returns the usual inner product of A and B, i.e., the sum of A(i,j)*B(i,j) over all i,j 
  inline Real ip(const Matrix& A, const Matrix& B) {
    chk_add(A, B);
    return mat_ip(A.nc * A.nr, A.m, B.m);
  }


  ///returns the squared Frobenius norm of column j of A, i.e., the sum of A(i,j)*A(i,j) over all i with possibly (if scaling!=0) each term i multiplied by (*scaling)(i)
  inline Real colip(const Matrix& A, Integer j, const Matrix* scaling = 0) {
    chk_init(A);
    chk_range(j, j, A.nc, A.nc);
    if (scaling == 0)
      return mat_ip(A.nr, A.m + j * A.nr);
    assert(scaling->dim() == A.nr);
    const Real* mp = A.m + j * A.nr;
    const Real* const mend = mp + A.nr;
    const Real* sp = scaling->m;
    Real sum = 0;
    for (; mp != mend;) {
      const Real d = (*mp++);
      sum += d * d * (*sp++);
    }
    return sum;
  }

  ///returns the squared Frobenius norm of row i of A, i.e., the sum of A(i,j)*A(i,j) over all j with possibly (if scaling!=0) each term j multiplied by (*scaling)(j)
  inline Real rowip(const Matrix& A, Integer i, const Matrix* scaling = 0) {
    chk_init(A);
    chk_range(i, i, A.nr, A.nr);
    if (scaling == 0)
      return mat_ip(A.nc, A.m + i, A.nr);
    assert(scaling->dim() == A.nc);
    const Real* mp = A.m + i;
    const Real* const mend = mp + A.nc * A.nr;
    const Real* sp = scaling->m;
    Real sum = 0;
    for (; mp != mend; mp += A.nr) {
      const Real d = (*mp);
      sum += d * d * (*sp++);
    }
    return sum;
  }

  ///returns the column vector of the squared Frobenius norm of all columns j of A, i.e., the sum of A(i,j)*A(i,j) over all i for each j
  inline Matrix colsip(const Matrix& A) {
    chk_init(A);
    Matrix tmp(1, A.nc); chk_set_init(tmp, 1);
    for (Integer j = 0; j < A.nc; j++)
      tmp[j] = colip(A, j);
    return tmp;
  }

  ///returns the row vector of the squared Frobenius norm of all rowd i of A, i.e., the sum of A(i,j)*A(i,j) over all j for each i
  inline Matrix rowsip(const Matrix& A) {
    chk_init(A);
    Matrix tmp(A.nr, 1, 0.);
    Real* mp = A.m;
    for (Integer j = 0; j < A.nc; j++) {
      Real* tp = tmp.m;
      for (Integer i = 0; i < A.nr; i++, mp++)
        (*tp++) += (*mp) * (*mp);
    }
    return tmp;
  }

  ///returns the Frobenius norm of A, i.e., the square root of the sum of A(i,j)*A(i,j) over all i,j 
  inline Real norm2(const Matrix& A) {
    chk_init(A);
    return ::sqrt(mat_ip(A.nc * A.nr, A.m));
  }

  inline Matrix& Matrix::operator=(const Matrix& A) {
    return xeya(A);
  }

  inline Matrix& Matrix::operator*=(const Matrix& A) {
    Matrix C; return xeya(genmult(*this, A, C));
  }

  inline Matrix& Matrix::operator+=(const Matrix& A) {
    return xpeya(A);
  }

  inline Matrix& Matrix::operator-=(const Matrix& A) {
    return xpeya(A, -1.);
  }

  Matrix& Matrix::operator%=(const Matrix& A) {
    chk_add(*this, A); mat_xhadey(nr * nc, m, A.m); return *this;
  }

  Matrix& Matrix::operator/=(const Matrix& A) {
    chk_add(*this, A); mat_xinvhadey(nr * nc, m, A.m); return *this;
  }

  inline Matrix Matrix::operator-() const {
    return Matrix(*this, -1.);
  }

  inline Matrix& Matrix::operator*=(Real d) {
    chk_init(*this); mat_xmultea(nr * nc, m, d); return *this;
  }

  inline Matrix& Matrix::operator/=(Real d) {
    chk_init(*this); mat_xmultea(nr * nc, m, 1. / d); return *this;
  }

  inline Matrix& Matrix::operator+=(Real d) {
    chk_init(*this); mat_xpea(nr * nc, m, d); return *this;
  }

  inline Matrix& Matrix::operator-=(Real d) {
    chk_init(*this); mat_xpea(nr * nc, m, -d); return *this;
  }

  /// returns Matrix equal to A*B
  inline Matrix operator*(const Matrix& A, const Matrix& B) {
    Matrix C; return genmult(A, B, C);
  }

  /// returns Matrix equal to A+B
  inline Matrix operator+(const Matrix& A, const Matrix& B) {
    Matrix C; return xeyapzb(C, A, B, 1., 1.);
  }

  /// returns Matrix equal to A-B
  inline Matrix operator-(const Matrix& A, const Matrix& B) {
    Matrix C; return xeyapzb(C, A, B, 1., -1.);
  }

  /// ATTENTION: this is redefined as the Hadamard product, C(i,j)=A(i,j)*B(i,j) for all i,j
  inline Matrix operator%(const Matrix& A, const Matrix& B) {
    Matrix C(A); return C %= B;
  }

  /// ATTENTION: this is redefined to act componentwise without checking for zeros, C(i,j)=A(i,j)/B(i,j) for all i,j
  inline Matrix operator/(const Matrix& A, const Matrix& B) {
    Matrix C(A); return C /= B;
  }

  /// returns Matrix equal to A*d
  inline Matrix operator*(const Matrix& A, Real d) {
    return Matrix(A, d);
  }

  /// returns Matrix equal to d*A
  inline Matrix operator*(Real d, const Matrix& A) {
    return Matrix(A, d);
  }

  /// returns Matrix equal to A/d;  ATTENTION: d is NOT checked for 0
  inline Matrix operator/(const Matrix& A, Real d) {
    return Matrix(A, 1. / d);
  }

  /// returns (i,j)=A(i,j)+d for all i,j
  inline Matrix operator+(const Matrix& A, Real d) {
    Matrix B(A); return B += d;
  }

  /// returns (i,j)=A(i,j)+d for all i,j
  inline Matrix operator+(Real d, const Matrix& A) {
    Matrix B(A); return B += d;
  }

  /// returns (i,j)=A(i,j)-d for all i,j
  inline Matrix operator-(const Matrix& A, Real d) {
    Matrix B(A); return B -= d;
  }

  /// returns (i,j)=d-A(i,j) for all i,j
  inline Matrix operator-(Real d, const Matrix& A) {
    Matrix B(A, -1.); return B += d;
  }

  inline int Matrix::QR_factor(Matrix& Q, Matrix& R, Real tol = 1e-10) const  //*this is unchanged
  {
    R = *this; return R.QR_factor(Q, tol);
  }
  inline int Matrix::QR_factor(Matrix& Q, Matrix& R, Indexmatrix& piv, Real tol = 1e-10) const {
    R = *this; return R.QR_factor(Q, piv, tol);
  }
  inline int Matrix::ls(Matrix& rhs, Real tol = 1e-10) {
    return this->QR_solve(rhs, tol);
  }

  /// computes a Householder QR factorization of A and outputs Q and R leaving A unchanged; always returns 0  
  inline int QR_factor(const Matrix& A, Matrix& Q, Matrix& R, Real tol) {
    return A.QR_factor(Q, R, tol);
  }

  /// computes a Householder QR factorization of A with pivating. It outputs Q, R, and the pivoting permuation in piv; returns the rank of A  
  inline int QR_factor(const Matrix& A, Matrix& Q, Matrix& R, Indexmatrix& piv, Real tol) {
    return A.QR_factor(Q, R, piv, tol);
  }

  /// retuns a matrix that keeps the upper triangle of A starting with diagonal d, i.e., (i,j)=A(i,j) for 0<=i<row dimension, max(0,i+d)<=j<column dimension, and sets (i,j)=0 otherwise
  inline Matrix triu(const Matrix& A, Integer i = 0) {
    Matrix B(A); B.triu(i); return B;
  }

  /// retuns a matrix that keeps the lower triangle of A starting with diagonal d, i.e., (i,j)=A(i,j) for 0<=i<row dimension, 0<=j<min(i+d+1,column dimension), and sets (i,j)=0 otherwise
  inline Matrix tril(const Matrix& A, Integer i = 0) {
    Matrix B(A); B.tril(i); return B;
  }

  /// returns a new matrix [A, B], i.e., it concats matrices A and B rowwise; A or B may be a 0x0 matrix
  inline Matrix concat_right(const Matrix& A, const Matrix& B) {
    Matrix C(A.dim() + B.dim(), 1); C = A; C.concat_right(B); return C;
  }

  /// returns a bew matrix [A; B], i.e., it concats matrices A and B columnwise; A or B may be a 0x0 matrix
  inline Matrix concat_below(const Matrix& A, const Matrix& B) {
    Matrix C(A.dim() + B.dim(), 1); C = A; C.concat_below(B); return C;
  }

  /// swap the content of the two matrices A and B (involves no copying)
  inline void swap(Matrix& A, Matrix& B) {
    Real* hm = A.m; A.m = B.m; B.m = hm;
    Integer hi = A.nr; A.nr = B.nr; B.nr = hi;
    hi = A.nc; A.nc = B.nc; B.nc = hi;
    hi = A.mem_dim; A.mem_dim = B.mem_dim; B.mem_dim = hi;
#if (CONICBUNDLE_DEBUG>=1)
    bool hb = A.is_init; A.is_init = B.is_init; B.is_init = hb;
#endif
  }

  /// return a nr x nc matrix with (i,j) assigned a random number uniformly from [0,1] for all i,j
  inline Matrix rand(Integer rows, Integer cols, CH_Tools::GB_rand* random_generator = 0) {
    Matrix A; return A.rand(rows, cols, random_generator);
  }

  /// returns a matrix with elements (i,j)=1./((*this)(i,j)) for all i,j; ATTENTION: no check for division by zero 
  inline Matrix inv(const Matrix& A) {
    Matrix B(A); return B.inv();
  }

  /// returns a matrix with elements (i,j)=sqrt((*this)(i,j)) for all i,j 
  inline Matrix sqrt(const Matrix& A) {
    Matrix B(A); return B.sqrt();
  }

  /// returns a matrix with elements (i,j)=sqr((*this)(i,j)) for all i,j 
  inline Matrix sqr(const Matrix& A) {
    Matrix B(A); return B.sqr();
  }

  /// returns a matrix with elements (i,j)=sign((*this)(i,j)) for all i,j using ::sign(double,double)
  inline Matrix sign(const Matrix& A, Real tol = 1e-12) {
    Matrix B(A); return B.sign(tol);
  }

  /// returns a matrix with elements (i,j)=floor((*this)(i,j)) for all i,j 
  inline Matrix floor(const Matrix& A) {
    Matrix B(A); return B.floor();
  }

  /// returns a matrix with elements (i,j)=ceil((*this)(i,j)) for all i,j 
  inline Matrix ceil(const Matrix& A) {
    Matrix B(A); return B.floor();
  }

  /// returns a matrix with elements (i,j)=rint((*this)(i,j)) for all i,j 
  inline Matrix rint(const Matrix& A) {
    Matrix B(A); return B.floor();
  }

  /// returns a matrix with elements (i,j)=round((*this)(i,j)) for all i,j 
  inline Matrix round(const Matrix& A) {
    Matrix B(A); return B.floor();
  }


  /// returns a matrix having elements (i,j)=Real(A(i,j)>d) for all i,j
  inline Matrix operator>(const Matrix& A, const Matrix& B) {
    return B < A;
  }

  /// returns a matrix having elements (i,j)=Real(A(i,j)>=d) for all i,j
  inline Matrix operator>=(const Matrix& A, const Matrix& B) {
    return B <= A;
  }
  /// returns a matrix having elements (i,j)=Real(d<A(i,j)) for all i,j
  inline Matrix operator<(Real d, const Matrix& A) {
    return A > d;
  }
  /// returns a matrix having elements (i,j)=Real(d>A(i,j)) for all i,j
  inline Matrix operator>(Real d, const Matrix& A) {
    return A < d;
  }

  /// returns a matrix having elements (i,j)=Real(d<=A(i,j)) for all i,j
  inline Matrix operator<=(Real d, const Matrix& A) {
    return A >= d;
  }

  /// returns a matrix having elements (i,j)=Real(d>=A(i,j)) for all i,j
  inline Matrix operator>=(Real d, const Matrix& A) {
    return A <= d;
  }

  /// returns a matrix having elements (i,j)=Real(d==A(i,j)) for all i,j
  inline Matrix operator==(Real d, const Matrix& A) {
    return A == d;
  }

  /// returns a matrix having elements (i,j)=Real(d!=A(i,j)) for all i,j
  inline Matrix operator!=(Real d, const Matrix& A) {
    return A != d;
  }

  /// returns true if both matrices have the same size and the same elements    
  inline bool equal(const Matrix& A, const Matrix& B) {
    return ((A.nc == B.nc) && (A.nr == B.nr) && (mat_equal(A.nr * A.nc, A.get_store(), B.get_store())));
  }

  inline Indexmatrix sortindex(const Matrix& vec, bool nondecreasing) {
    Indexmatrix ind; sortindex(vec, ind, nondecreasing); return ind;
  }


  /// returns an Indexmatrix ind so that A(ind(i)) 0<=i<ind.dim() runs through all nonzero elements with abs(A(j))>tol 
  inline Indexmatrix find(const Matrix& A, Real tol = 1e-10) {
    return A.find(tol);
  }
  /// returns an Indexmatrix ind so that A(ind(i)) 0<=i<ind.dim() runs through all elements of A having value num, i.e., abs(A(j)-num)<tol 
  inline Indexmatrix find_number(const Matrix& A, Real num = 0., Real tol = 1e-10) {
    return A.find_number(num, tol);
  }

  inline std::vector<double>& assign(std::vector<double>& vec, const Matrix& A) {
    chk_init(A);
    vec.resize((unsigned long)(A.dim()));
    for (Integer i = 0; i < A.dim(); i++) vec[(unsigned long)(i)] = A(i);
    return vec;
  }

}


//in order to define the inlines associated with the other matrix
//classes, the hxx-files of these are included here, as well. 

#ifndef CH_MATRIX_CLASSES__SYMMAT_HXX
#include "symmat.hxx"
#endif

#endif
