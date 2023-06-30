@doc raw"""
    CBMatrix()

empty matrix
"""
CBMatrix() = CBMatrix(@ccall libcb.cb_matrix_new()::Ptr{Cvoid})

@doc raw"""
    CBMatrix(param0::CBMatrix, d::Real = 1., atrans::Integer = 0)

copy constructor, *this=d*A
"""
CBMatrix(param0::CBMatrix, d::Real = 1., atrans::Integer = 0) = CBMatrix(@ccall libcb.cb_matrix_new2(param0.data::Ptr{Cvoid}, d::Cdouble, atrans::Cint)::Ptr{Cvoid})

@doc raw"""
    CBMatrix(param0::AbstractRange{<:Real}, param0_tol::Real = 1e-8)

generate a column vector holding the elements of this Realrange
"""
CBMatrix(param0::AbstractRange{<:Real}, param0_tol::Real = 1e-8) = CBMatrix(@ccall libcb.cb_matrix_new3(first(param0)::Cdouble, last(param0)::Cdouble, step(param0)::Cdouble, param0_tol::Cdouble)::Ptr{Cvoid})

@doc raw"""
    CBMatrix(nr::Integer, nc::Integer)

* @brief generate a matrix of size nr x nc but WITHOUT initializing the memory

        If initializing the memory externally and CONICBUNDLE_DEBUG is defined, please use
        set_init() via matrix.set_init(true) in order to avoid warnings concerning improper
        initialization
    
"""
CBMatrix(nr::Integer, nc::Integer) = CBMatrix(@ccall libcb.cb_matrix_new4(nr::Cint, nc::Cint)::Ptr{Cvoid})

@doc raw"""
    CBMatrix(nr::Integer, nc::Integer, d::Real)

generate a matrix of size nr x nc initializing all elements to the value d
"""
CBMatrix(nr::Integer, nc::Integer, d::Real) = CBMatrix(@ccall libcb.cb_matrix_new5(nr::Cint, nc::Cint, d::Cdouble)::Ptr{Cvoid})

@doc raw"""
    CBMatrix(nr::Integer, nc::Integer, dp::Union{<:AbstractVector{Cdouble},Nothing}, d::Real = 1.)

generate a matrix of size nr x nc initializing the elements from the (one dimensional) array dp with increment incr and scaled by d
"""
CBMatrix(nr::Integer, nc::Integer, dp::Union{<:AbstractVector{Cdouble},Nothing}, d::Real = 1.) = GC.@preserve dp begin
    CBMatrix(@ccall libcb.cb_matrix_new6(nr::Cint, nc::Cint, dp::Ptr{Cdouble}, stride(dp, 1)::Cint, d::Cdouble)::Ptr{Cvoid})
end

@doc raw"""
    cb_init!(self::CBMatrix, A::CBMatrix, d::Real = 1., atrans::Integer = 0)

initialize to *this=A*d where A may be transposed
"""
cb_init!(self::CBMatrix, A::CBMatrix, d::Real = 1., atrans::Integer = 0) = (@ccall libcb.cb_matrix_init(self.data::Ptr{Cvoid}, A.data::Ptr{Cvoid}, d::Cdouble, atrans::Cint)::Ptr{Cvoid}; return self)

@doc raw"""
    cb_init!(self::CBMatrix, A::CBIndexmatrix, d::Real = 1.)

initialize to *this=A*d
"""
cb_init!(self::CBMatrix, A::CBIndexmatrix, d::Real = 1.) = (@ccall libcb.cb_matrix_init2(self.data::Ptr{Cvoid}, A.data::Ptr{Cvoid}, d::Cdouble)::Ptr{Cvoid}; return self)

@doc raw"""
    cb_init!(self::CBMatrix, A::CBSparsemat, d::Real = 1.)

initialize to *this=A*d
"""
cb_init!(self::CBMatrix, A::CBSparsemat, d::Real = 1.) = (@ccall libcb.cb_matrix_init3(self.data::Ptr{Cvoid}, A.data::Ptr{Cvoid}, d::Cdouble)::Ptr{Cvoid}; return self)

@doc raw"""
    cb_init!(self::CBMatrix, S::CBSymmatrix, d::Real = 1.)

initialize to *this=A*d
"""
cb_init!(self::CBMatrix, S::CBSymmatrix, d::Real = 1.) = (@ccall libcb.cb_matrix_init4(self.data::Ptr{Cvoid}, S.data::Ptr{Cvoid}, d::Cdouble)::Ptr{Cvoid}; return self)

@doc raw"""
    cb_init!(self::CBMatrix, param0::CBSparsesym, d::Real = 1.)

initialize to *this=A*d
"""
cb_init!(self::CBMatrix, param0::CBSparsesym, d::Real = 1.) = (@ccall libcb.cb_matrix_init5(self.data::Ptr{Cvoid}, param0.data::Ptr{Cvoid}, d::Cdouble)::Ptr{Cvoid}; return self)

@doc raw"""
    cb_init!(self::CBMatrix, param0::AbstractRange{<:Real}, param0_tol::Real = 1e-8)

initialize *this to a column vector holding the elements of Realrange
"""
cb_init!(self::CBMatrix, param0::AbstractRange{<:Real}, param0_tol::Real = 1e-8) = (@ccall libcb.cb_matrix_init6(self.data::Ptr{Cvoid}, first(param0)::Cdouble, last(param0)::Cdouble, step(param0)::Cdouble, param0_tol::Cdouble)::Ptr{Cvoid}; return self)

@doc raw"""
    cb_init!(self::CBMatrix, nr::Integer, nc::Integer, d::Real)

intialize *this to a matrix of size nr x nc initializing all elements to the value d
"""
cb_init!(self::CBMatrix, nr::Integer, nc::Integer, d::Real) = (@ccall libcb.cb_matrix_init7(self.data::Ptr{Cvoid}, nr::Cint, nc::Cint, d::Cdouble)::Ptr{Cvoid}; return self)

@doc raw"""
    cb_init!(self::CBMatrix, nr::Integer, nc::Integer, dp::Union{<:AbstractVector{Cdouble},Nothing}, d::Real = 1.)

generate a matrix of size nr x nc initializing the elements from the (one dimensional) array dp with increment incr and scaled by d
"""
cb_init!(self::CBMatrix, nr::Integer, nc::Integer, dp::Union{<:AbstractVector{Cdouble},Nothing}, d::Real = 1.) = GC.@preserve dp begin
    (@ccall libcb.cb_matrix_init8(self.data::Ptr{Cvoid}, nr::Cint, nc::Cint, dp::Ptr{Cdouble}, stride(dp, 1)::Cint, d::Cdouble)::Ptr{Cvoid}; return self)
end

@doc raw"""
    cb_init_diag!(self::CBMatrix, nr::Integer, d::Real = 1.)

initialize to a diagonal nr x nr matrix with constant diagonal value d
"""
cb_init_diag!(self::CBMatrix, nr::Integer, d::Real = 1.) = (@ccall libcb.cb_matrix_init_diag(self.data::Ptr{Cvoid}, nr::Cint, d::Cdouble)::Ptr{Cvoid}; return self)

@doc raw"""
    cb_init_diag!(self::CBMatrix, vec::CBMatrix, d::Real = 1.)

initialize to a diagonal matrix with diagonal given by vec
"""
cb_init_diag!(self::CBMatrix, vec::CBMatrix, d::Real = 1.) = (@ccall libcb.cb_matrix_init_diag2(self.data::Ptr{Cvoid}, vec.data::Ptr{Cvoid}, d::Cdouble)::Ptr{Cvoid}; return self)

@doc raw"""
    cb_init_diag!(self::CBMatrix, vec::CBIndexmatrix, d::Real = 1.)

initialize to a diagonal matrix with diagonal given by vec
"""
cb_init_diag!(self::CBMatrix, vec::CBIndexmatrix, d::Real = 1.) = (@ccall libcb.cb_matrix_init_diag3(self.data::Ptr{Cvoid}, vec.data::Ptr{Cvoid}, d::Cdouble)::Ptr{Cvoid}; return self)

@doc raw"""
    cb_newsize!(self::CBMatrix, nr::Integer, nc::Integer)

* @brief resize the matrix to nr x nc elements but WITHOUT initializing the memory

        If initializing the memory externally and CONICBUNDLE_DEBUG is defined, please use
        set_init() via matrix.set_init(true) in order to avoid warnings concerning improper
        initialization
    
"""
cb_newsize!(self::CBMatrix, nr::Integer, nc::Integer) = @ccall libcb.cb_matrix_newsize(self.data::Ptr{Cvoid}, nr::Cint, nc::Cint)::Cvoid

@doc raw"""
    CBMatrix(A::CBIndexmatrix, d::Real = 1.)

(*this)=d*A
"""
CBMatrix(A::CBIndexmatrix, d::Real = 1.) = CBMatrix(@ccall libcb.cb_matrix_new7(A.data::Ptr{Cvoid}, d::Cdouble)::Ptr{Cvoid})

@doc raw"""
    CBMatrix(A::CBSparsemat, d::Real = 1.)

(*this)=d*A
"""
CBMatrix(A::CBSparsemat, d::Real = 1.) = CBMatrix(@ccall libcb.cb_matrix_new8(A.data::Ptr{Cvoid}, d::Cdouble)::Ptr{Cvoid})

@doc raw"""
    CBMatrix(S::CBSymmatrix, d::Real = 1.)

(*this)=d*A
"""
CBMatrix(S::CBSymmatrix, d::Real = 1.) = CBMatrix(@ccall libcb.cb_matrix_new9(S.data::Ptr{Cvoid}, d::Cdouble)::Ptr{Cvoid})

@doc raw"""
    CBMatrix(param0::CBSparsesym, d::Real = 1.)

(*this)=d*A
"""
CBMatrix(param0::CBSparsesym, d::Real = 1.) = CBMatrix(@ccall libcb.cb_matrix_new10(param0.data::Ptr{Cvoid}, d::Cdouble)::Ptr{Cvoid})

@doc raw"""
    cb_dim2(self::CBMatrix)

returns the number of rows in _nr and the number of columns in _nc
"""
function cb_dim2(self::CBMatrix)
    _nc = Ref{Int}()
    _nr = Ref{Int}()
    @ccall libcb.cb_matrix_dim(self.data::Ptr{Cvoid}, _nr::Ref{Int}, _nc::Ref{Int})::Cvoid
    return _nr[], _nc[]
end

@doc raw"""
    cb_dim(self::CBMatrix)

returns the dimension rows * columns when the matrix is regarded as a vector
"""
cb_dim(self::CBMatrix) = @ccall libcb.cb_matrix_dim2(self.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_rowdim(self::CBMatrix)

returns the row dimension
"""
cb_rowdim(self::CBMatrix) = @ccall libcb.cb_matrix_rowdim(self.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_coldim(self::CBMatrix)

returns the column dimension
"""
cb_coldim(self::CBMatrix) = @ccall libcb.cb_matrix_coldim(self.data::Ptr{Cvoid})::Cint

@doc raw"""
    Base.setindex!(self::CBMatrix, value::Real, i::Integer, j::Integer)

returns reference to element (i,j) of the matrix (rowindex i, columnindex j)
"""
Base.setindex!(self::CBMatrix, value::Real, i::Integer, j::Integer) = @ccall libcb.cb_matrix_set(self.data::Ptr{Cvoid}, i::Cint, j::Cint, value::Cdouble)::Cvoid

@doc raw"""
    Base.setindex!(self::CBMatrix, value::Real, i::Integer)

returns reference to element (i) of the matrix if regarded as vector of stacked columns [element (i%rowdim, i/rowdim)]
"""
Base.setindex!(self::CBMatrix, value::Real, i::Integer) = @ccall libcb.cb_matrix_set2(self.data::Ptr{Cvoid}, i::Cint, value::Cdouble)::Cvoid

@doc raw"""
    Base.getindex(self::CBMatrix, i::Integer, j::Integer)

returns value of element (i,j) of the matrix (rowindex i, columnindex j)
"""
Base.getindex(self::CBMatrix, i::Integer, j::Integer) = @ccall libcb.cb_matrix_get(self.data::Ptr{Cvoid}, i::Cint, j::Cint)::Cdouble

@doc raw"""
    Base.getindex(self::CBMatrix, i::Integer)

returns value of element (i) of the matrix if regarded as vector of stacked columns [element (i%rowdim, i/rowdim)]
"""
Base.getindex(self::CBMatrix, i::Integer) = @ccall libcb.cb_matrix_get2(self.data::Ptr{Cvoid}, i::Cint)::Cdouble

@doc raw"""
    Base.getindex(self::CBMatrix, vecrow::CBIndexmatrix, veccol::CBIndexmatrix)

returns a new submatrix as indexed by vecrow and veccol, A(i,j)=(*this)(vecrow(i),veccol(j)) for 0<=i<vecrow.dim(), 0<=j<veccol.dim()
"""
Base.getindex(self::CBMatrix, vecrow::CBIndexmatrix, veccol::CBIndexmatrix) = CBMatrix(@ccall libcb.cb_matrix_new_get(self.data::Ptr{Cvoid}, vecrow.data::Ptr{Cvoid}, veccol.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    Base.getindex(self::CBMatrix, A::CBIndexmatrix)

returns a new matrix B of the same shape as A with B(i,j)=(*this)(A(i),A(j)) for 0<=i<A.rowdim(), 0<=j<A.coldim()
"""
Base.getindex(self::CBMatrix, A::CBIndexmatrix) = CBMatrix(@ccall libcb.cb_matrix_new_get2(self.data::Ptr{Cvoid}, A.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_col(self::CBMatrix, i::Integer)

returns column i copied to a new matrix
"""
cb_col(self::CBMatrix, i::Integer) = CBMatrix(@ccall libcb.cb_matrix_new_col(self.data::Ptr{Cvoid}, i::Cint)::Ptr{Cvoid})

@doc raw"""
    cb_row(self::CBMatrix, i::Integer)

returns row i copied to a new matrix
"""
cb_row(self::CBMatrix, i::Integer) = CBMatrix(@ccall libcb.cb_matrix_new_row(self.data::Ptr{Cvoid}, i::Cint)::Ptr{Cvoid})

@doc raw"""
    cb_cols(self::CBMatrix, vec::CBIndexmatrix)

returns a matrix of size this->rowdim() x vec.dim(), with column i a copy of column vec(i) of *this
"""
cb_cols(self::CBMatrix, vec::CBIndexmatrix) = CBMatrix(@ccall libcb.cb_matrix_new_cols(self.data::Ptr{Cvoid}, vec.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_rows(self::CBMatrix, vec::CBIndexmatrix)

returns a matrix of size vec.dim() x this->coldim(), with row i a copy of row vec(i) of *this
"""
cb_rows(self::CBMatrix, vec::CBIndexmatrix) = CBMatrix(@ccall libcb.cb_matrix_new_rows(self.data::Ptr{Cvoid}, vec.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_swap_rowsij!(self::CBMatrix, i::Integer, j::Integer)

row i of this matrix is swapped with row j
"""
cb_swap_rowsij!(self::CBMatrix, i::Integer, j::Integer) = (@ccall libcb.cb_matrix_swap_rowsij(self.data::Ptr{Cvoid}, i::Cint, j::Cint)::Ptr{Cvoid}; return self)

@doc raw"""
    cb_swap_colsij!(self::CBMatrix, i::Integer, j::Integer)

column i of this matrix is swapped with column j
"""
cb_swap_colsij!(self::CBMatrix, i::Integer, j::Integer) = (@ccall libcb.cb_matrix_swap_colsij(self.data::Ptr{Cvoid}, i::Cint, j::Cint)::Ptr{Cvoid}; return self)

@doc raw"""
    cb_pivot_permute_rows!(self::CBMatrix, piv::CBIndexmatrix, inverse::Bool = false)

for i=0 to rowdim row i of this matrix is swapped with row piv(i); for inverse =true the inverse permutation is generated
"""
cb_pivot_permute_rows!(self::CBMatrix, piv::CBIndexmatrix, inverse::Bool = false) = (@ccall libcb.cb_matrix_pivot_permute_rows(self.data::Ptr{Cvoid}, piv.data::Ptr{Cvoid}, inverse::Cint)::Ptr{Cvoid}; return self)

@doc raw"""
    cb_pivot_permute_cols!(self::CBMatrix, piv::CBIndexmatrix, inverse::Bool = false)

for j=0 to coldim column  j of this matrix is swapped with column piv(j); for inverse =true the inverse permutation is generated
"""
cb_pivot_permute_cols!(self::CBMatrix, piv::CBIndexmatrix, inverse::Bool = false) = (@ccall libcb.cb_matrix_pivot_permute_cols(self.data::Ptr{Cvoid}, piv.data::Ptr{Cvoid}, inverse::Cint)::Ptr{Cvoid}; return self)

@doc raw"""
    cb_triu!(self::CBMatrix, d::Integer = 0)

keeps everything above and including diagonal d, everything below is set to zero, returns *this
"""
cb_triu!(self::CBMatrix, d::Integer = 0) = (@ccall libcb.cb_matrix_triu(self.data::Ptr{Cvoid}, d::Cint)::Ptr{Cvoid}; return self)

@doc raw"""
    cb_tril!(self::CBMatrix, d::Integer = 0)

keeps everything below and including diagonal d, everything above is set to zero, returns *this
"""
cb_tril!(self::CBMatrix, d::Integer = 0) = (@ccall libcb.cb_matrix_tril(self.data::Ptr{Cvoid}, d::Cint)::Ptr{Cvoid}; return self)

@doc raw"""
    cb_subassign!(self::CBMatrix, vecrow::CBIndexmatrix, veccol::CBIndexmatrix, A::CBMatrix)

assigns A to a submatrix of *this,  (*this)(vecrow(i),veccol(j))=A(i,j) for 0<=i<vecrow.dim(), 0<=j<veccol.dim()
"""
cb_subassign!(self::CBMatrix, vecrow::CBIndexmatrix, veccol::CBIndexmatrix, A::CBMatrix) = (@ccall libcb.cb_matrix_subassign(self.data::Ptr{Cvoid}, vecrow.data::Ptr{Cvoid}, veccol.data::Ptr{Cvoid}, A.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

@doc raw"""
    cb_subassign!(self::CBMatrix, vec::CBIndexmatrix, A::CBMatrix)

assigns vector A to a subvector of *this,  (*this)(vec(i))=A(i) for 0<=i<vec.dim(), *this, vec, and A may be rectangular matrices, their dimesions are not changed, returns *this
"""
cb_subassign!(self::CBMatrix, vec::CBIndexmatrix, A::CBMatrix) = (@ccall libcb.cb_matrix_subassign2(self.data::Ptr{Cvoid}, vec.data::Ptr{Cvoid}, A.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

@doc raw"""
    cb_delete_rows!(self::CBMatrix, ind::CBIndexmatrix, sorted_increasingly::Bool = false)

all rows indexed by vector ind are deleted, no row should appear twice in ind, remaining rows are moved up keeping their order, returns *this
"""
cb_delete_rows!(self::CBMatrix, ind::CBIndexmatrix, sorted_increasingly::Bool = false) = (@ccall libcb.cb_matrix_delete_rows(self.data::Ptr{Cvoid}, ind.data::Ptr{Cvoid}, sorted_increasingly::Cint)::Ptr{Cvoid}; return self)

@doc raw"""
    cb_delete_cols!(self::CBMatrix, ind::CBIndexmatrix, sorted_increasingly::Bool = false)

all colmuns indexed by vector ind are deleted, no column should appear twice in ind, remaining columns are moved left keeping their order, returns *this
"""
cb_delete_cols!(self::CBMatrix, ind::CBIndexmatrix, sorted_increasingly::Bool = false) = (@ccall libcb.cb_matrix_delete_cols(self.data::Ptr{Cvoid}, ind.data::Ptr{Cvoid}, sorted_increasingly::Cint)::Ptr{Cvoid}; return self)

@doc raw"""
    cb_insert_row!(self::CBMatrix, i::Integer, v::CBMatrix)

insert the row vector v before row i, 0<=i<= row dimension, for i==row dimension the row is appended below; appending to a 0x0 matrix is allowed, returns *this
"""
cb_insert_row!(self::CBMatrix, i::Integer, v::CBMatrix) = (@ccall libcb.cb_matrix_insert_row(self.data::Ptr{Cvoid}, i::Cint, v.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

@doc raw"""
    cb_insert_col!(self::CBMatrix, i::Integer, v::CBMatrix)

insert a column before column i, 0<=i<= column dimension, for i==column dimension the column is appended at the right; appending to a 0x0 matrix is allowed, returns *this
"""
cb_insert_col!(self::CBMatrix, i::Integer, v::CBMatrix) = (@ccall libcb.cb_matrix_insert_col(self.data::Ptr{Cvoid}, i::Cint, v.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

@doc raw"""
    cb_reduce_length!(self::CBMatrix, n::Integer)

(*this) is set to a column vector of length min{max{0,n},dim()}; usually used to truncate a vector, returns *this
"""
cb_reduce_length!(self::CBMatrix, n::Integer) = (@ccall libcb.cb_matrix_reduce_length(self.data::Ptr{Cvoid}, n::Cint)::Ptr{Cvoid}; return self)

@doc raw"""
    cb_concat_right!(self::CBMatrix, A::CBMatrix, Atrans::Integer = 0)

concats matrix A (or its tranpose) to the right of *this, A or *this may be the 0x0 matrix initally, returns *this
"""
cb_concat_right!(self::CBMatrix, A::CBMatrix, Atrans::Integer = 0) = (@ccall libcb.cb_matrix_concat_right(self.data::Ptr{Cvoid}, A.data::Ptr{Cvoid}, Atrans::Cint)::Ptr{Cvoid}; return self)

@doc raw"""
    cb_concat_below!(self::CBMatrix, A::CBMatrix)

concats matrix A to the bottom of *this, A or *this may be the 0x0 matrix initally, returns *this
"""
cb_concat_below!(self::CBMatrix, A::CBMatrix) = (@ccall libcb.cb_matrix_concat_below(self.data::Ptr{Cvoid}, A.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

@doc raw"""
    cb_concat_right!(self::CBMatrix, d::Real)

concat value d at the bottom of *this, *this must be a column vector or the 0x0 matrix, returns *this
"""
cb_concat_right!(self::CBMatrix, d::Real) = (@ccall libcb.cb_matrix_concat_right2(self.data::Ptr{Cvoid}, d::Cdouble)::Ptr{Cvoid}; return self)

@doc raw"""
    cb_concat_below!(self::CBMatrix, d::Real)

concat value d at the right of *this, *this must be a row vector or the 0x0 matrix, returns *this
"""
cb_concat_below!(self::CBMatrix, d::Real) = (@ccall libcb.cb_matrix_concat_below2(self.data::Ptr{Cvoid}, d::Cdouble)::Ptr{Cvoid}; return self)

@doc raw"""
    cb_enlarge_right!(self::CBMatrix, addnc::Integer)

enlarge the matrix by addnc>=0 columns without intializaton of the new columns, returns *this (marked as not initialized if nr>0)
"""
cb_enlarge_right!(self::CBMatrix, addnc::Integer) = (@ccall libcb.cb_matrix_enlarge_right(self.data::Ptr{Cvoid}, addnc::Cint)::Ptr{Cvoid}; return self)

@doc raw"""
    cb_enlarge_below!(self::CBMatrix, addnr::Integer)

enlarge the matrix by addnr>=0 rows without intializaton of the new rows, returns *this (marked as not initialized if addnr>0)
"""
cb_enlarge_below!(self::CBMatrix, addnr::Integer) = (@ccall libcb.cb_matrix_enlarge_below(self.data::Ptr{Cvoid}, addnr::Cint)::Ptr{Cvoid}; return self)

@doc raw"""
    cb_enlarge_right!(self::CBMatrix, addnc::Integer, d::Real)

enlarge the matrix by addnc>=0 columns intializing the new columns by value d, returns *this
"""
cb_enlarge_right!(self::CBMatrix, addnc::Integer, d::Real) = (@ccall libcb.cb_matrix_enlarge_right2(self.data::Ptr{Cvoid}, addnc::Cint, d::Cdouble)::Ptr{Cvoid}; return self)

@doc raw"""
    cb_enlarge_below!(self::CBMatrix, addnr::Integer, d::Real)

enlarge the matrix by addnr>=0 rows intializing the new rows by value d, returns *this
"""
cb_enlarge_below!(self::CBMatrix, addnr::Integer, d::Real) = (@ccall libcb.cb_matrix_enlarge_below2(self.data::Ptr{Cvoid}, addnr::Cint, d::Cdouble)::Ptr{Cvoid}; return self)

@doc raw"""
    cb_enlarge_right!(self::CBMatrix, addnc::Integer, dp::Union{<:AbstractVector{Cdouble},Nothing}, d::Real = 1.)

enlarge the matrix by addnc>=0 columns intializing the new columns by the values pointed to by dp times d, returns *this
"""
cb_enlarge_right!(self::CBMatrix, addnc::Integer, dp::Union{<:AbstractVector{Cdouble},Nothing}, d::Real = 1.) = GC.@preserve dp begin
    (LinearAlgebra.chkstride1(dp); @ccall libcb.cb_matrix_enlarge_right3(self.data::Ptr{Cvoid}, addnc::Cint, dp::Ptr{Cdouble}, d::Cdouble)::Ptr{Cvoid}; return self)
end

@doc raw"""
    cb_enlarge_below!(self::CBMatrix, addnr::Integer, dp::Union{<:AbstractVector{Cdouble},Nothing}, d::Real = 1.)

enlarge the matrix by addnr>=0 rows intializing the new rows by the values pointed to by dp times d, returns *this
"""
cb_enlarge_below!(self::CBMatrix, addnr::Integer, dp::Union{<:AbstractVector{Cdouble},Nothing}, d::Real = 1.) = GC.@preserve dp begin
    (LinearAlgebra.chkstride1(dp); @ccall libcb.cb_matrix_enlarge_below3(self.data::Ptr{Cvoid}, addnr::Cint, dp::Ptr{Cdouble}, d::Cdouble)::Ptr{Cvoid}; return self)
end

@doc raw"""
    cb_get_store(self::CBMatrix)

returns the current address of the internal value array; use cautiously!
"""
cb_get_store(self::CBMatrix) = @ccall libcb.cb_matrix_get_store2(self.data::Ptr{Cvoid})::Ptr{Cdouble}

@doc raw"""
    cb_diag(A::CBMatrix)

returns a column vector v consisting of the elements v(i)=(*this)(i,i), 0<=i<min(row dimension,column dimension)
"""
cb_diag(A::CBMatrix) = CBMatrix(@ccall libcb.cb_matrix_new_diag(A.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_triu(A::CBMatrix, i::Integer)

retuns a matrix that keeps the upper triangle of A starting with diagonal d, i.e., (i,j)=A(i,j) for 0<=i<row dimension, max(0,i+d)<=j<column dimension, and sets (i,j)=0 otherwise
"""
cb_triu(A::CBMatrix, i::Integer) = CBMatrix(@ccall libcb.cb_matrix_new_triu(A.data::Ptr{Cvoid}, i::Cint)::Ptr{Cvoid})

@doc raw"""
    cb_tril(A::CBMatrix, i::Integer)

retuns a matrix that keeps the lower triangle of A starting with diagonal d, i.e., (i,j)=A(i,j) for 0<=i<row dimension, 0<=j<min(i+d+1,column dimension), and sets (i,j)=0 otherwise
"""
cb_tril(A::CBMatrix, i::Integer) = CBMatrix(@ccall libcb.cb_matrix_new_tril(A.data::Ptr{Cvoid}, i::Cint)::Ptr{Cvoid})

@doc raw"""
    cb_concat_right(A::CBMatrix, B::CBMatrix)

returns a new matrix [A, B], i.e., it concats matrices A and B rowwise; A or B may be a 0x0 matrix
"""
cb_concat_right(A::CBMatrix, B::CBMatrix) = CBMatrix(@ccall libcb.cb_matrix_new_concat_right(A.data::Ptr{Cvoid}, B.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_concat_below(A::CBMatrix, B::CBMatrix)

returns a bew matrix [A; B], i.e., it concats matrices A and B columnwise; A or B may be a 0x0 matrix
"""
cb_concat_below(A::CBMatrix, B::CBMatrix) = CBMatrix(@ccall libcb.cb_matrix_new_concat_below(A.data::Ptr{Cvoid}, B.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_swap(A::CBMatrix, B::CBMatrix)

swap the content of the two matrices A and B (involves no copying)
"""
cb_swap(A::CBMatrix, B::CBMatrix) = @ccall libcb.cb_matrix_swap(A.data::Ptr{Cvoid}, B.data::Ptr{Cvoid})::Cvoid

@doc raw"""
    cb_xeya!(self::CBMatrix, A::CBMatrix, d::Real = 1., atrans::Integer = 0)

sets *this=d*A where A may be transposed and returns *this
"""
cb_xeya!(self::CBMatrix, A::CBMatrix, d::Real = 1., atrans::Integer = 0) = (@ccall libcb.cb_matrix_xeya(self.data::Ptr{Cvoid}, A.data::Ptr{Cvoid}, d::Cdouble, atrans::Cint)::Ptr{Cvoid}; return self)

@doc raw"""
    cb_xpeya!(self::CBMatrix, A::CBMatrix, d::Real = 1.)

sets *this+=d*A and returns *this
"""
cb_xpeya!(self::CBMatrix, A::CBMatrix, d::Real = 1.) = (@ccall libcb.cb_matrix_xpeya(self.data::Ptr{Cvoid}, A.data::Ptr{Cvoid}, d::Cdouble)::Ptr{Cvoid}; return self)

@doc raw"""
    cb_xeya!(self::CBMatrix, A::CBIndexmatrix, d::Real = 1.)

sets *this=d*A and returns *this
"""
cb_xeya!(self::CBMatrix, A::CBIndexmatrix, d::Real = 1.) = (@ccall libcb.cb_matrix_xeya2(self.data::Ptr{Cvoid}, A.data::Ptr{Cvoid}, d::Cdouble)::Ptr{Cvoid}; return self)

@doc raw"""
    cb_xpeya!(self::CBMatrix, A::CBIndexmatrix, d::Real = 1.)

sets *this+=d*A and returns *this
"""
cb_xpeya!(self::CBMatrix, A::CBIndexmatrix, d::Real = 1.) = (@ccall libcb.cb_matrix_xpeya2(self.data::Ptr{Cvoid}, A.data::Ptr{Cvoid}, d::Cdouble)::Ptr{Cvoid}; return self)

@doc raw"""
    cb_xbpeya(x::CBMatrix, y::CBMatrix, alpha::Real, beta::Real, ytrans::Integer)

returns x= alpha*y+beta*x, where y may be transposed (ytrans=1); if beta==0. then x is initialized to the correct size
"""
cb_xbpeya(x::CBMatrix, y::CBMatrix, alpha::Real, beta::Real, ytrans::Integer) = (@ccall libcb.cb_matrix_xbpeya(x.data::Ptr{Cvoid}, y.data::Ptr{Cvoid}, alpha::Cdouble, beta::Cdouble, ytrans::Cint)::Ptr{Cvoid}; return self)

@doc raw"""
    cb_xeyapzb(x::CBMatrix, y::CBMatrix, z::CBMatrix, alpha::Real, beta::Real)

returns x= alpha*y+beta*z; x is initialized to the correct size
"""
cb_xeyapzb(x::CBMatrix, y::CBMatrix, z::CBMatrix, alpha::Real, beta::Real) = (@ccall libcb.cb_matrix_xeyapzb(x.data::Ptr{Cvoid}, y.data::Ptr{Cvoid}, z.data::Ptr{Cvoid}, alpha::Cdouble, beta::Cdouble)::Ptr{Cvoid}; return self)

@doc raw"""
    cb_genmult(A::CBMatrix, B::CBMatrix, C::CBMatrix, alpha::Real, beta::Real, atrans::Integer, btrans::Integer)

returns C=beta*C+alpha*A*B, where A and B may be transposed; C must not be equal to A and B; if beta==0. then C is initialized to the correct size
"""
cb_genmult(A::CBMatrix, B::CBMatrix, C::CBMatrix, alpha::Real, beta::Real, atrans::Integer, btrans::Integer) = (@ccall libcb.cb_matrix_genmult(A.data::Ptr{Cvoid}, B.data::Ptr{Cvoid}, C.data::Ptr{Cvoid}, alpha::Cdouble, beta::Cdouble, atrans::Cint, btrans::Cint)::Ptr{Cvoid}; return self)

Base.copy!(self::CBMatrix, A::CBMatrix) = (@ccall libcb.cb_matrix_assign(self.data::Ptr{Cvoid}, A.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

MA.operate!(::typeof(*), self::CBMatrix, s::CBMatrix) = (@ccall libcb.cb_matrix_times(self.data::Ptr{Cvoid}, s.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

MA.promote_operation(::typeof(*), ::Type{<:CBMatrix}, ::Type{<:CBMatrix}) = CBMatrix

MA.operate!(::typeof(+), self::CBMatrix, v::CBMatrix) = (@ccall libcb.cb_matrix_plus(self.data::Ptr{Cvoid}, v.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

MA.promote_operation(::typeof(+), ::Type{<:CBMatrix}, ::Type{<:CBMatrix}) = CBMatrix

MA.operate!(::typeof(-), self::CBMatrix, v::CBMatrix) = (@ccall libcb.cb_matrix_minus(self.data::Ptr{Cvoid}, v.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

MA.promote_operation(::typeof(-), ::Type{<:CBMatrix}, ::Type{<:CBMatrix}) = CBMatrix

@doc raw"""
    MA.operate!(::typeof(%), self::CBMatrix, A::CBMatrix)

ATTENTION: this is redefined as the Hadamard product, (*this)(i,j)=(*this)(i,j)*A(i,j) for all i,j
"""
MA.operate!(::typeof(%), self::CBMatrix, A::CBMatrix) = (@ccall libcb.cb_matrix_rem(self.data::Ptr{Cvoid}, A.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

MA.promote_operation(::typeof(%), ::Type{<:CBMatrix}, ::Type{<:CBMatrix}) = CBMatrix

@doc raw"""
    MA.operate!(::typeof(/), self::CBMatrix, A::CBMatrix)

ATTENTION: this is redefined to act componentwise without checking for zeros, (*this)(i,j)=(*this)(i,j)/A(i,j) for all i,j
"""
MA.operate!(::typeof(/), self::CBMatrix, A::CBMatrix) = (@ccall libcb.cb_matrix_divide(self.data::Ptr{Cvoid}, A.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

MA.promote_operation(::typeof(/), ::Type{<:CBMatrix}, ::Type{<:CBMatrix}) = CBMatrix

Base.:-(self::CBMatrix) = CBMatrix(@ccall libcb.cb_matrix_new_minus(self.data::Ptr{Cvoid})::Ptr{Cvoid})

MA.operate!(::typeof(*), self::CBMatrix, d::Real) = (@ccall libcb.cb_matrix_times2(self.data::Ptr{Cvoid}, d::Cdouble)::Ptr{Cvoid}; return self)

MA.promote_operation(::typeof(*), ::Type{<:CBMatrix}, ::Type{<:Real}) = CBMatrix

@doc raw"""
    MA.operate!(::typeof(/), self::CBMatrix, d::Real)

ATTENTION: d is NOT checked for 0
"""
MA.operate!(::typeof(/), self::CBMatrix, d::Real) = (@ccall libcb.cb_matrix_divide2(self.data::Ptr{Cvoid}, d::Cdouble)::Ptr{Cvoid}; return self)

MA.promote_operation(::typeof(/), ::Type{<:CBMatrix}, ::Type{<:Real}) = CBMatrix

@doc raw"""
    MA.operate!(::typeof(+), self::CBMatrix, d::Real)

sets (*this)(i,j)+=d for all i,j
"""
MA.operate!(::typeof(+), self::CBMatrix, d::Real) = (@ccall libcb.cb_matrix_plus2(self.data::Ptr{Cvoid}, d::Cdouble)::Ptr{Cvoid}; return self)

MA.promote_operation(::typeof(+), ::Type{<:CBMatrix}, ::Type{<:Real}) = CBMatrix

@doc raw"""
    MA.operate!(::typeof(-), self::CBMatrix, d::Real)

sets (*this)(i,j)-=d for all i,j
"""
MA.operate!(::typeof(-), self::CBMatrix, d::Real) = (@ccall libcb.cb_matrix_minus2(self.data::Ptr{Cvoid}, d::Cdouble)::Ptr{Cvoid}; return self)

MA.promote_operation(::typeof(-), ::Type{<:CBMatrix}, ::Type{<:Real}) = CBMatrix

@doc raw"""
    cb_transpose!(self::CBMatrix)

transposes itself (cheap for vectors, expensive for matrices)
"""
cb_transpose!(self::CBMatrix) = (@ccall libcb.cb_matrix_transpose(self.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

Base.:*(A::CBMatrix, B::CBMatrix) = CBMatrix(@ccall libcb.cb_matrix_new_times(A.data::Ptr{Cvoid}, B.data::Ptr{Cvoid})::Ptr{Cvoid})

Base.:+(A::CBMatrix, B::CBMatrix) = CBMatrix(@ccall libcb.cb_matrix_new_plus(A.data::Ptr{Cvoid}, B.data::Ptr{Cvoid})::Ptr{Cvoid})

Base.:-(A::CBMatrix, B::CBMatrix) = CBMatrix(@ccall libcb.cb_matrix_new_minus2(A.data::Ptr{Cvoid}, B.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    Base.:%(A::CBMatrix, B::CBMatrix)

ATTENTION: this is redefined as the Hadamard product, C(i,j)=A(i,j)*B(i,j) for all i,j
"""
Base.:%(A::CBMatrix, B::CBMatrix) = CBMatrix(@ccall libcb.cb_matrix_new_rem(A.data::Ptr{Cvoid}, B.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    Base.:/(A::CBMatrix, B::CBMatrix)

ATTENTION: this is redefined to act componentwise without checking for zeros, C(i,j)=A(i,j)/B(i,j) for all i,j
"""
Base.:/(A::CBMatrix, B::CBMatrix) = CBMatrix(@ccall libcb.cb_matrix_new_divide(A.data::Ptr{Cvoid}, B.data::Ptr{Cvoid})::Ptr{Cvoid})

Base.:*(A::CBMatrix, d::Real) = CBMatrix(@ccall libcb.cb_matrix_new_times2(A.data::Ptr{Cvoid}, d::Cdouble)::Ptr{Cvoid})

Base.:*(d::Real, A::CBMatrix) = CBMatrix(@ccall libcb.cb_matrix_new_times3(d::Cdouble, A.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    Base.:/(A::CBMatrix, d::Real)

ATTENTION: d is NOT checked for 0
"""
Base.:/(A::CBMatrix, d::Real) = CBMatrix(@ccall libcb.cb_matrix_new_divide2(A.data::Ptr{Cvoid}, d::Cdouble)::Ptr{Cvoid})

@doc raw"""
    Base.:+(A::CBMatrix, d::Real)

returns (i,j)=A(i,j)+d for all i,j
"""
Base.:+(A::CBMatrix, d::Real) = CBMatrix(@ccall libcb.cb_matrix_new_plus2(A.data::Ptr{Cvoid}, d::Cdouble)::Ptr{Cvoid})

@doc raw"""
    Base.:+(d::Real, A::CBMatrix)

returns (i,j)=A(i,j)+d for all i,j
"""
Base.:+(d::Real, A::CBMatrix) = CBMatrix(@ccall libcb.cb_matrix_new_plus3(d::Cdouble, A.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    Base.:-(A::CBMatrix, d::Real)

returns (i,j)=A(i,j)-d for all i,j
"""
Base.:-(A::CBMatrix, d::Real) = CBMatrix(@ccall libcb.cb_matrix_new_minus3(A.data::Ptr{Cvoid}, d::Cdouble)::Ptr{Cvoid})

@doc raw"""
    Base.:-(d::Real, A::CBMatrix)

returns (i,j)=d-A(i,j) for all i,j
"""
Base.:-(d::Real, A::CBMatrix) = CBMatrix(@ccall libcb.cb_matrix_new_minus4(d::Cdouble, A.data::Ptr{Cvoid})::Ptr{Cvoid})

cb_transpose(A::CBMatrix) = CBMatrix(@ccall libcb.cb_matrix_new_transpose(A.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_xeya!(self::CBMatrix, A::CBSymmatrix, d::Real = 1.)

sets *this=d*A and returns *this
"""
cb_xeya!(self::CBMatrix, A::CBSymmatrix, d::Real = 1.) = (@ccall libcb.cb_matrix_xeya3(self.data::Ptr{Cvoid}, A.data::Ptr{Cvoid}, d::Cdouble)::Ptr{Cvoid}; return self)

@doc raw"""
    cb_xpeya!(self::CBMatrix, A::CBSymmatrix, d::Real = 1.)

sets *this+=d*A and returns *this
"""
cb_xpeya!(self::CBMatrix, A::CBSymmatrix, d::Real = 1.) = (@ccall libcb.cb_matrix_xpeya3(self.data::Ptr{Cvoid}, A.data::Ptr{Cvoid}, d::Cdouble)::Ptr{Cvoid}; return self)

Base.copy!(self::CBMatrix, S::CBSymmatrix) = (@ccall libcb.cb_matrix_assign2(self.data::Ptr{Cvoid}, S.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

MA.operate!(::typeof(*), self::CBMatrix, S::CBSymmatrix) = (@ccall libcb.cb_matrix_times3(self.data::Ptr{Cvoid}, S.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

MA.promote_operation(::typeof(*), ::Type{<:CBMatrix}, ::Type{<:CBSymmatrix}) = CBMatrix

MA.operate!(::typeof(+), self::CBMatrix, S::CBSymmatrix) = (@ccall libcb.cb_matrix_plus3(self.data::Ptr{Cvoid}, S.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

MA.promote_operation(::typeof(+), ::Type{<:CBMatrix}, ::Type{<:CBSymmatrix}) = CBMatrix

MA.operate!(::typeof(-), self::CBMatrix, S::CBSymmatrix) = (@ccall libcb.cb_matrix_minus3(self.data::Ptr{Cvoid}, S.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

MA.promote_operation(::typeof(-), ::Type{<:CBMatrix}, ::Type{<:CBSymmatrix}) = CBMatrix

@doc raw"""
    cb_xeya!(self::CBMatrix, A::CBSparsesym, d::Real = 1.)

sets *this=d*A and returns *this
"""
cb_xeya!(self::CBMatrix, A::CBSparsesym, d::Real = 1.) = (@ccall libcb.cb_matrix_xeya4(self.data::Ptr{Cvoid}, A.data::Ptr{Cvoid}, d::Cdouble)::Ptr{Cvoid}; return self)

@doc raw"""
    cb_xpeya!(self::CBMatrix, A::CBSparsesym, d::Real = 1.)

sets *this+=d*A and returns *this
"""
cb_xpeya!(self::CBMatrix, A::CBSparsesym, d::Real = 1.) = (@ccall libcb.cb_matrix_xpeya4(self.data::Ptr{Cvoid}, A.data::Ptr{Cvoid}, d::Cdouble)::Ptr{Cvoid}; return self)

Base.copy!(self::CBMatrix, param0::CBSparsesym) = (@ccall libcb.cb_matrix_assign3(self.data::Ptr{Cvoid}, param0.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

MA.operate!(::typeof(*), self::CBMatrix, S::CBSparsesym) = (@ccall libcb.cb_matrix_times4(self.data::Ptr{Cvoid}, S.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

MA.promote_operation(::typeof(*), ::Type{<:CBMatrix}, ::Type{<:CBSparsesym}) = CBMatrix

MA.operate!(::typeof(+), self::CBMatrix, S::CBSparsesym) = (@ccall libcb.cb_matrix_plus4(self.data::Ptr{Cvoid}, S.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

MA.promote_operation(::typeof(+), ::Type{<:CBMatrix}, ::Type{<:CBSparsesym}) = CBMatrix

MA.operate!(::typeof(-), self::CBMatrix, S::CBSparsesym) = (@ccall libcb.cb_matrix_minus4(self.data::Ptr{Cvoid}, S.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

MA.promote_operation(::typeof(-), ::Type{<:CBMatrix}, ::Type{<:CBSparsesym}) = CBMatrix

@doc raw"""
    cb_xeya!(self::CBMatrix, A::CBSparsemat, d::Real = 1.)

sets *this=d*A and returns *this
"""
cb_xeya!(self::CBMatrix, A::CBSparsemat, d::Real = 1.) = (@ccall libcb.cb_matrix_xeya5(self.data::Ptr{Cvoid}, A.data::Ptr{Cvoid}, d::Cdouble)::Ptr{Cvoid}; return self)

@doc raw"""
    cb_xpeya!(self::CBMatrix, A::CBSparsemat, d::Real = 1.)

sets *this+=d*A and returns *this
"""
cb_xpeya!(self::CBMatrix, A::CBSparsemat, d::Real = 1.) = (@ccall libcb.cb_matrix_xpeya5(self.data::Ptr{Cvoid}, A.data::Ptr{Cvoid}, d::Cdouble)::Ptr{Cvoid}; return self)

Base.copy!(self::CBMatrix, A::CBSparsemat) = (@ccall libcb.cb_matrix_assign4(self.data::Ptr{Cvoid}, A.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

MA.operate!(::typeof(*), self::CBMatrix, A::CBSparsemat) = (@ccall libcb.cb_matrix_times5(self.data::Ptr{Cvoid}, A.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

MA.promote_operation(::typeof(*), ::Type{<:CBMatrix}, ::Type{<:CBSparsemat}) = CBMatrix

MA.operate!(::typeof(+), self::CBMatrix, A::CBSparsemat) = (@ccall libcb.cb_matrix_plus5(self.data::Ptr{Cvoid}, A.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

MA.promote_operation(::typeof(+), ::Type{<:CBMatrix}, ::Type{<:CBSparsemat}) = CBMatrix

MA.operate!(::typeof(-), self::CBMatrix, A::CBSparsemat) = (@ccall libcb.cb_matrix_minus5(self.data::Ptr{Cvoid}, A.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

MA.promote_operation(::typeof(-), ::Type{<:CBMatrix}, ::Type{<:CBSparsemat}) = CBMatrix

@doc raw"""
    cb_genmult(A::CBSymmatrix, B::CBMatrix, C::CBMatrix, alpha::Real, beta::Real, btrans::Integer)

returns C=beta*C+alpha*A*B, where A and B may be transposed; C must not be equal to A and B; if beta==0. then C is initialized to the correct size
"""
cb_genmult(A::CBSymmatrix, B::CBMatrix, C::CBMatrix, alpha::Real, beta::Real, btrans::Integer) = (@ccall libcb.cb_matrix_genmult2(A.data::Ptr{Cvoid}, B.data::Ptr{Cvoid}, C.data::Ptr{Cvoid}, alpha::Cdouble, beta::Cdouble, btrans::Cint)::Ptr{Cvoid}; return self)

@doc raw"""
    cb_genmult(A::CBMatrix, B::CBSymmatrix, C::CBMatrix, alpha::Real, beta::Real, atrans::Integer)

returns C=beta*C+alpha*A*B, where A and B may be transposed; C must not be equal to A and B; if beta==0. then C is initialized to the correct size
"""
cb_genmult(A::CBMatrix, B::CBSymmatrix, C::CBMatrix, alpha::Real, beta::Real, atrans::Integer) = (@ccall libcb.cb_matrix_genmult3(A.data::Ptr{Cvoid}, B.data::Ptr{Cvoid}, C.data::Ptr{Cvoid}, alpha::Cdouble, beta::Cdouble, atrans::Cint)::Ptr{Cvoid}; return self)

@doc raw"""
    cb_genmult(A::CBSparsesym, B::CBMatrix, C::CBMatrix, alpha::Real, beta::Real, btrans::Integer)

returns C=beta*C+alpha*A*B, where A and B may be transposed; C must not be equal to A and B; if beta==0. then C is initialized to the correct size
"""
cb_genmult(A::CBSparsesym, B::CBMatrix, C::CBMatrix, alpha::Real, beta::Real, btrans::Integer) = (@ccall libcb.cb_matrix_genmult4(A.data::Ptr{Cvoid}, B.data::Ptr{Cvoid}, C.data::Ptr{Cvoid}, alpha::Cdouble, beta::Cdouble, btrans::Cint)::Ptr{Cvoid}; return self)

@doc raw"""
    cb_genmult(A::CBMatrix, B::CBSparsesym, C::CBMatrix, alpha::Real, beta::Real, atrans::Integer)

returns C=beta*C+alpha*A*B, where A and B may be transposed; C must not be equal to A and B; if beta==0. then C is initialized to the correct size
"""
cb_genmult(A::CBMatrix, B::CBSparsesym, C::CBMatrix, alpha::Real, beta::Real, atrans::Integer) = (@ccall libcb.cb_matrix_genmult5(A.data::Ptr{Cvoid}, B.data::Ptr{Cvoid}, C.data::Ptr{Cvoid}, alpha::Cdouble, beta::Cdouble, atrans::Cint)::Ptr{Cvoid}; return self)

@doc raw"""
    cb_genmult(A::CBSparsemat, B::CBMatrix, C::CBMatrix, alpha::Real, beta::Real, atrans::Integer, btrans::Integer)

returns C=beta*C+alpha*A*B, where A and B may be transposed; C must not be equal to A and B; if beta==0. then C is initialized to the correct size
"""
cb_genmult(A::CBSparsemat, B::CBMatrix, C::CBMatrix, alpha::Real, beta::Real, atrans::Integer, btrans::Integer) = (@ccall libcb.cb_matrix_genmult6(A.data::Ptr{Cvoid}, B.data::Ptr{Cvoid}, C.data::Ptr{Cvoid}, alpha::Cdouble, beta::Cdouble, atrans::Cint, btrans::Cint)::Ptr{Cvoid}; return self)

@doc raw"""
    cb_genmult(A::CBSparsemat, B::CBMatrix, colB::Integer, C::CBMatrix, colC::Integer, alpha::Real, beta::Real, atrans::Integer, btrans::Integer)

returns C.col(colC)=beta*C.col(colC)+alpha*A*B.col(colB), where A and B may be transposed first; C must not be equal to A and B; if beta==0. then C is initialized, but the size of C must be correct already
"""
cb_genmult(A::CBSparsemat, B::CBMatrix, colB::Integer, C::CBMatrix, colC::Integer, alpha::Real, beta::Real, atrans::Integer, btrans::Integer) = (@ccall libcb.cb_matrix_genmult7(A.data::Ptr{Cvoid}, B.data::Ptr{Cvoid}, colB::Cint, C.data::Ptr{Cvoid}, colC::Cint, alpha::Cdouble, beta::Cdouble, atrans::Cint, btrans::Cint)::Ptr{Cvoid}; return self)

@doc raw"""
    cb_genmult(A::CBMatrix, B::CBSparsemat, C::CBMatrix, alpha::Real, beta::Real, atrans::Integer, btrans::Integer)

returns C=beta*C+alpha*A*B, where A and B may be transposed; C must not be equal to A and B; if beta==0. then C is initialized to the correct size
"""
cb_genmult(A::CBMatrix, B::CBSparsemat, C::CBMatrix, alpha::Real, beta::Real, atrans::Integer, btrans::Integer) = (@ccall libcb.cb_matrix_genmult8(A.data::Ptr{Cvoid}, B.data::Ptr{Cvoid}, C.data::Ptr{Cvoid}, alpha::Cdouble, beta::Cdouble, atrans::Cint, btrans::Cint)::Ptr{Cvoid}; return self)

@doc raw"""
    cb_rand!(self::CBMatrix, nr::Integer, nc::Integer, random_generator::Union{<:CBGB_rand,Nothing} = nothing)

resize *this to an nr x nc matrix and assign to (i,j) a random number uniformly from [0,1] for all i,j
"""
cb_rand!(self::CBMatrix, nr::Integer, nc::Integer, random_generator::Union{<:CBGB_rand,Nothing} = nothing) = (@ccall libcb.cb_matrix_rand(self.data::Ptr{Cvoid}, nr::Cint, nc::Cint, (isnothing(random_generator) ? C_NULL : random_generator.data)::Ptr{Cvoid})::Ptr{Cvoid}; return self)

@doc raw"""
    cb_rand_normal!(self::CBMatrix, nr::Integer, nc::Integer, mean::Real = 0., variance::Real = 1., generator_type::Integer = 0)

resize *this to an nr x nc matrix and assign to (i,j) a random number from the normal distribution with given mean and variance (generators: 0 std, 1 mt, 2 mt64)
"""
cb_rand_normal!(self::CBMatrix, nr::Integer, nc::Integer, mean::Real = 0., variance::Real = 1., generator_type::Integer = 0) = (@ccall libcb.cb_matrix_rand_normal(self.data::Ptr{Cvoid}, nr::Cint, nc::Cint, mean::Cdouble, variance::Cdouble, generator_type::Cint)::Ptr{Cvoid}; return self)

@doc raw"""
    cb_shuffle!(self::CBMatrix, random_generator::Union{<:CBGB_rand,Nothing} = nothing)

shuffle the elements randomly (does not change dimensions)
"""
cb_shuffle!(self::CBMatrix, random_generator::Union{<:CBGB_rand,Nothing} = nothing) = (@ccall libcb.cb_matrix_shuffle(self.data::Ptr{Cvoid}, (isnothing(random_generator) ? C_NULL : random_generator.data)::Ptr{Cvoid})::Ptr{Cvoid}; return self)

@doc raw"""
    cb_inv!(self::CBMatrix)

sets (*this)(i,j)=1./(*this)(i,j) for all i,j and returns *this
"""
cb_inv!(self::CBMatrix) = (@ccall libcb.cb_matrix_inv(self.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

@doc raw"""
    cb_sqrt!(self::CBMatrix)

sets (*this)(i,j)=sqrt((*this)(i,j)) for all i,j and returns *this
"""
cb_sqrt!(self::CBMatrix) = (@ccall libcb.cb_matrix_sqrt(self.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

@doc raw"""
    cb_sqr!(self::CBMatrix)

sets (*this)(i,j)=sqr((*this)(i,j)) for all i,j and returns *this
"""
cb_sqr!(self::CBMatrix) = (@ccall libcb.cb_matrix_sqr(self.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

@doc raw"""
    cb_sign!(self::CBMatrix, tol::Real = 1e-12)

sets (*this)(i,j)=sign((*this)(i,j),tol) for all i,j using ::sign(double,double) and returns *this
"""
cb_sign!(self::CBMatrix, tol::Real = 1e-12) = (@ccall libcb.cb_matrix_sign(self.data::Ptr{Cvoid}, tol::Cdouble)::Ptr{Cvoid}; return self)

@doc raw"""
    cb_floor!(self::CBMatrix)

sets (*this)(i,j)=floor((*this)(i,j)) for all i,j and returns *this
"""
cb_floor!(self::CBMatrix) = (@ccall libcb.cb_matrix_floor(self.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

@doc raw"""
    cb_ceil!(self::CBMatrix)

sets (*this)(i,j)=ceil((*this)(i,j)) for all i,j and returns *this
"""
cb_ceil!(self::CBMatrix) = (@ccall libcb.cb_matrix_ceil(self.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

@doc raw"""
    cb_rint!(self::CBMatrix)

sets (*this)(i,j)=rint((*this)(i,j)) for all i,j and returns *this
"""
cb_rint!(self::CBMatrix) = (@ccall libcb.cb_matrix_rint(self.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

@doc raw"""
    cb_round!(self::CBMatrix)

sets (*this)(i,j)=round((*this)(i,j)) for all i,j and returns *this
"""
cb_round!(self::CBMatrix) = (@ccall libcb.cb_matrix_round(self.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

@doc raw"""
    cb_abs!(self::CBMatrix)

sets (*this)(i,j)=abs((*this)(i,j)) for all i,j and returns *this
"""
cb_abs!(self::CBMatrix) = (@ccall libcb.cb_matrix_abs(self.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

@doc raw"""
    cb_rand(nr::Integer, nc::Integer, random_generator::Union{<:CBGB_rand,Nothing})

return a nr x nc matrix with (i,j) assigned a random number uniformly from [0,1] for all i,j
"""
cb_rand(nr::Integer, nc::Integer, random_generator::Union{<:CBGB_rand,Nothing}) = CBMatrix(@ccall libcb.cb_matrix_new_rand(nr::Cint, nc::Cint, (isnothing(random_generator) ? C_NULL : random_generator.data)::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_inv(A::CBMatrix)

returns a matrix with elements (i,j)=1./((*this)(i,j)) for all i,j; ATTENTION: no check for division by zero
"""
cb_inv(A::CBMatrix) = CBMatrix(@ccall libcb.cb_matrix_new_inv(A.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_sqrt(A::CBMatrix)

returns a matrix with elements (i,j)=sqrt((*this)(i,j)) for all i,j
"""
cb_sqrt(A::CBMatrix) = CBMatrix(@ccall libcb.cb_matrix_new_sqrt(A.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_sqr(A::CBMatrix)

returns a matrix with elements (i,j)=sqr((*this)(i,j)) for all i,j
"""
cb_sqr(A::CBMatrix) = CBMatrix(@ccall libcb.cb_matrix_new_sqr(A.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_sign(A::CBMatrix, tol::Real)

returns a matrix with elements (i,j)=sign((*this)(i,j)) for all i,j using ::sign(double,double)
"""
cb_sign(A::CBMatrix, tol::Real) = CBMatrix(@ccall libcb.cb_matrix_new_sign(A.data::Ptr{Cvoid}, tol::Cdouble)::Ptr{Cvoid})

@doc raw"""
    cb_floor(A::CBMatrix)

returns a matrix with elements (i,j)=floor((*this)(i,j)) for all i,j
"""
cb_floor(A::CBMatrix) = CBMatrix(@ccall libcb.cb_matrix_new_floor(A.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_ceil(A::CBMatrix)

returns a matrix with elements (i,j)=ceil((*this)(i,j)) for all i,j
"""
cb_ceil(A::CBMatrix) = CBMatrix(@ccall libcb.cb_matrix_new_ceil(A.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_rint(A::CBMatrix)

returns a matrix with elements (i,j)=rint((*this)(i,j)) for all i,j
"""
cb_rint(A::CBMatrix) = CBMatrix(@ccall libcb.cb_matrix_new_rint(A.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_round(A::CBMatrix)

returns a matrix with elements (i,j)=round((*this)(i,j)) for all i,j
"""
cb_round(A::CBMatrix) = CBMatrix(@ccall libcb.cb_matrix_new_round(A.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_abs(A::CBMatrix)

returns a matrix with elements (i,j)=abs((*this)(i,j)) for all i,j
"""
cb_abs(A::CBMatrix) = CBMatrix(@ccall libcb.cb_matrix_new_abs(A.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_contains_nan!(self::CBMatrix)

returns true if any matrix element returns true on std::isnan
"""
cb_contains_nan!(self::CBMatrix) = Bool(@ccall libcb.cb_matrix_contains_nan(self.data::Ptr{Cvoid})::Cint)

@doc raw"""
    cb_scale_rows!(self::CBMatrix, vec::CBMatrix)

scales each row i of (*this) by vec(i), i.e., (*this)=diag(vec)*(*this), and returns (*this)
"""
cb_scale_rows!(self::CBMatrix, vec::CBMatrix) = (@ccall libcb.cb_matrix_scale_rows(self.data::Ptr{Cvoid}, vec.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

@doc raw"""
    cb_scale_cols!(self::CBMatrix, vec::CBMatrix)

scales each column i of (*this) by vec(i), i.e., (*this)=(*this)*diag(vec), and returns (*this)
"""
cb_scale_cols!(self::CBMatrix, vec::CBMatrix) = (@ccall libcb.cb_matrix_scale_cols(self.data::Ptr{Cvoid}, vec.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

@doc raw"""
    cb_triu_solve!(self::CBMatrix, rhs::CBMatrix, tol::Real = 1e-10)

solves (*this)*x=rhs for x by back substitution regarding (*this) as an upper triangle matrix and stores x in rhs. Returns 0 on success, otherwise i+1 if abs(*this)(i,i)<tol and the remaining row of rhs is nonzero.
"""
cb_triu_solve!(self::CBMatrix, rhs::CBMatrix, tol::Real = 1e-10) = @ccall libcb.cb_matrix_triu_solve(self.data::Ptr{Cvoid}, rhs.data::Ptr{Cvoid}, tol::Cdouble)::Cint

@doc raw"""
    cb_tril_solve!(self::CBMatrix, rhs::CBMatrix, tol::Real = 1e-10)

solves (*this)*x=rhs for x by forward substitution regarding (*this) as an upper triangle matrix and stores x in rhs. Returns 0 on success, otherwise i+1 if abs(*this)(i,i)<tol and the reduced row of rhs is nonzero.
"""
cb_tril_solve!(self::CBMatrix, rhs::CBMatrix, tol::Real = 1e-10) = @ccall libcb.cb_matrix_tril_solve(self.data::Ptr{Cvoid}, rhs.data::Ptr{Cvoid}, tol::Cdouble)::Cint

@doc raw"""
    cb_QR_factor!(self::CBMatrix, tol::Real = 1e-10)

computes a Householder QR_factorization overwriting (*this); returns 0 on success, otherwise column index +1 if the norm of the column is below tol when reaching it.
"""
cb_QR_factor!(self::CBMatrix, tol::Real = 1e-10) = @ccall libcb.cb_matrix_qr_factor(self.data::Ptr{Cvoid}, tol::Cdouble)::Cint

@doc raw"""
    cb_QR_factor!(self::CBMatrix, Q::CBMatrix, tol::Real = 1e-10)

computes a Householder QR_factorization computing the Q matrix explicitly and setting (*this)=R; it always returns 0
"""
cb_QR_factor!(self::CBMatrix, Q::CBMatrix, tol::Real = 1e-10) = @ccall libcb.cb_matrix_qr_factor2(self.data::Ptr{Cvoid}, Q.data::Ptr{Cvoid}, tol::Cdouble)::Cint

@doc raw"""
    cb_QR_factor(self::CBMatrix, Q::CBMatrix, R::CBMatrix, tol::Real)

computes a Householder QR_factorization computing matrices Q and R explicitly and leaving (*this) unchanged; it always returns 0
"""
cb_QR_factor(self::CBMatrix, Q::CBMatrix, R::CBMatrix, tol::Real) = @ccall libcb.cb_matrix_qr_factor3(self.data::Ptr{Cvoid}, Q.data::Ptr{Cvoid}, R.data::Ptr{Cvoid}, tol::Cdouble)::Cint

@doc raw"""
    cb_QR_factor!(self::CBMatrix, piv::CBIndexmatrix, tol::Real = 1e-10)

computes a Householder QR_factorization with pivoting  and overwriting (*this); the pivoting permutation is stored in piv; returns the rank
"""
cb_QR_factor!(self::CBMatrix, piv::CBIndexmatrix, tol::Real = 1e-10) = @ccall libcb.cb_matrix_qr_factor4(self.data::Ptr{Cvoid}, piv.data::Ptr{Cvoid}, tol::Cdouble)::Cint

@doc raw"""
    cb_QR_factor_relpiv!(self::CBMatrix, piv::CBIndexmatrix, tol::Real = 1e-10)

computes a Householder QR_factorization with pivoting on those with tolerable norm and tolerable reduction in norm and overwriting (*this); the pivoting permutation is stored in piv; returns the rank
"""
cb_QR_factor_relpiv!(self::CBMatrix, piv::CBIndexmatrix, tol::Real = 1e-10) = @ccall libcb.cb_matrix_qr_factor_relpiv(self.data::Ptr{Cvoid}, piv.data::Ptr{Cvoid}, tol::Cdouble)::Cint

@doc raw"""
    cb_QR_factor!(self::CBMatrix, Q::CBMatrix, piv::CBIndexmatrix, tol::Real = 1e-10)

Computes a Householder QR_factorization with pivoting computing the Q matrix explicitly and setting (*this)=R; the pivoting permutation is stored in piv; returns the rank
"""
cb_QR_factor!(self::CBMatrix, Q::CBMatrix, piv::CBIndexmatrix, tol::Real = 1e-10) = @ccall libcb.cb_matrix_qr_factor5(self.data::Ptr{Cvoid}, Q.data::Ptr{Cvoid}, piv.data::Ptr{Cvoid}, tol::Cdouble)::Cint

@doc raw"""
    cb_QR_factor(self::CBMatrix, Q::CBMatrix, R::CBMatrix, piv::CBIndexmatrix, tol::Real)

computes a Householder QR_factorization with pivoting computing matrices Q and R explicitly and leaving (*this) unchanged; the pivoting permutation is stored in piv; returns the rank
"""
cb_QR_factor(self::CBMatrix, Q::CBMatrix, R::CBMatrix, piv::CBIndexmatrix, tol::Real) = @ccall libcb.cb_matrix_qr_factor6(self.data::Ptr{Cvoid}, Q.data::Ptr{Cvoid}, R.data::Ptr{Cvoid}, piv.data::Ptr{Cvoid}, tol::Cdouble)::Cint

@doc raw"""
    cb_Qt_times(self::CBMatrix, A::CBMatrix, r::Integer)

computes A=transpose(Q)*A, assuming a housholder Q is coded in the first r columns of the lower triangle of (*this); it always returns 0
"""
cb_Qt_times(self::CBMatrix, A::CBMatrix, r::Integer) = @ccall libcb.cb_matrix_qt_times(self.data::Ptr{Cvoid}, A.data::Ptr{Cvoid}, r::Cint)::Cint

@doc raw"""
    cb_Q_times(self::CBMatrix, A::CBMatrix, r::Integer)

computes A=Q*A, assuming a housholder Q is coded in the first r columns of the lower triangle of (*this); it always returns 0
"""
cb_Q_times(self::CBMatrix, A::CBMatrix, r::Integer) = @ccall libcb.cb_matrix_q_times(self.data::Ptr{Cvoid}, A.data::Ptr{Cvoid}, r::Cint)::Cint

@doc raw"""
    cb_times_Q(self::CBMatrix, A::CBMatrix, r::Integer)

computes A=A*Q, assuming a housholder Q is coded in the first r columns of the lower triangle of (*this); it always returns 0
"""
cb_times_Q(self::CBMatrix, A::CBMatrix, r::Integer) = @ccall libcb.cb_matrix_times_q(self.data::Ptr{Cvoid}, A.data::Ptr{Cvoid}, r::Cint)::Cint

@doc raw"""
    cb_QR_solve!(self::CBMatrix, rhs::CBMatrix, tol::Real = 1e-10)

* @brief solves (*this)*x=rhs by factorizing and overwriting (*this); rhs is overwritten with the solution.  Returns 0 on success, otherwise i+1 if in the backsolve abs(*this)(i,i)<tol and the reduced row of rhs is nonzero.

        To avoid overwriting (*this), use the appropriate version
        of #QR_factor and #triu_solve.
     
"""
cb_QR_solve!(self::CBMatrix, rhs::CBMatrix, tol::Real = 1e-10) = @ccall libcb.cb_matrix_qr_solve(self.data::Ptr{Cvoid}, rhs.data::Ptr{Cvoid}, tol::Cdouble)::Cint

@doc raw"""
    cb_QR_concat_right!(self::CBMatrix, A::CBMatrix, piv::CBIndexmatrix, r::Integer, tol::Real = 1e-10)

* extend the current Householder QR-factorization stored in this by appending
        and factorizing the columns of A yielding a new QR_factorization

        @param[in] A contains the addtionial columns to be factorized

        @param[in,out] piv contains on input, the permution vector returned by QR_factor
            with pivoting, and on output the new entire permutation

        @param[in] r gives the initial rank of *this as returned by QR_factor with pviaton

        @param[in] tol gives the tolerance for regarding a vector as having norm zero

        @return the rank of the new QR-facotrization
    
"""
cb_QR_concat_right!(self::CBMatrix, A::CBMatrix, piv::CBIndexmatrix, r::Integer, tol::Real = 1e-10) = @ccall libcb.cb_matrix_qr_concat_right(self.data::Ptr{Cvoid}, A.data::Ptr{Cvoid}, piv.data::Ptr{Cvoid}, r::Cint, tol::Cdouble)::Cint

@doc raw"""
    cb_ls!(self::CBMatrix, rhs::CBMatrix, tol::Real)

computes a least squares solution by #QR_solve, overwriting (*this). rhs is overwritten with the solution. In fact, the full code is return this->QR_solve(rhs,tol);
"""
cb_ls!(self::CBMatrix, rhs::CBMatrix, tol::Real) = @ccall libcb.cb_matrix_ls(self.data::Ptr{Cvoid}, rhs.data::Ptr{Cvoid}, tol::Cdouble)::Cint

@doc raw"""
    cb_nnls(self::CBMatrix, rhs::CBMatrix, dual::Union{<:CBMatrix,Nothing} = nothing, tol::Real = 1e-10)

* @brief computes a nonnegative least squares solution; rhs is overwritten by the solution; if dual!=0, the dual variables are stored there; returns 0 on success, 1 on failure

      Computes the least squares solution of min ||Ax-b|| s.t. x >=0;\n
      The KKT system A'*A*x - A'*b - l = 0; x >=0, l>=0, x'*l=0 is solved
      solved by interior point method with QR-solution of the extended system.

      The current implementation is based on
      [P. Matsoms, "Sparse Linear Least Squares Problems in Optimization",
      Comput. Opt. and Appl., 7, 89-110 (1997)] but is only
      a quick and rather sloppy implementation of it ...
    
"""
cb_nnls(self::CBMatrix, rhs::CBMatrix, dual::Union{<:CBMatrix,Nothing} = nothing, tol::Real = 1e-10) = @ccall libcb.cb_matrix_nnls(self.data::Ptr{Cvoid}, rhs.data::Ptr{Cvoid}, (isnothing(dual) ? C_NULL : dual.data)::Ptr{Cvoid}, tol::Cdouble)::Cint

@doc raw"""
    cb_trace(A::CBMatrix)

returns the sum of the diagonal elements A(i,i) over all i
"""
cb_trace(A::CBMatrix) = @ccall libcb.cb_matrix_trace(A.data::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_ip(A::CBMatrix, B::CBMatrix)

returns the usual inner product of A and B, i.e., the sum of A(i,j)*B(i,j) over all i,j
"""
cb_ip(A::CBMatrix, B::CBMatrix) = @ccall libcb.cb_matrix_ip(A.data::Ptr{Cvoid}, B.data::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_ip_min_max(A::CBMatrix, B::CBMatrix)

returns in addition to the usual inner product of A and B the minimum and maximum value of A(i,j)*B(i,j) over all (i,j)
"""
function cb_ip_min_max(A::CBMatrix, B::CBMatrix)
    maxval = Ref{Float64}()
    minval = Ref{Float64}()
    @ccall libcb.cb_matrix_ip_min_max(A.data::Ptr{Cvoid}, B.data::Ptr{Cvoid}, minval::Ref{Float64}, maxval::Ref{Float64})::Cdouble
    return minval[], maxval[]
end

@doc raw"""
    cb_colip(A::CBMatrix, j::Integer, scaling::Union{<:CBMatrix,Nothing})

returns the squared Frobenius norm of column j of A, i.e., the sum of A(i,j)*A(i,j) over all i with possibly (if scaling!=0) each term i multiplied by (*scaling)(i)
"""
cb_colip(A::CBMatrix, j::Integer, scaling::Union{<:CBMatrix,Nothing}) = @ccall libcb.cb_matrix_colip(A.data::Ptr{Cvoid}, j::Cint, (isnothing(scaling) ? C_NULL : scaling.data)::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_rowip(A::CBMatrix, i::Integer, scaling::Union{<:CBMatrix,Nothing})

returns the squared Frobenius norm of row i of A, i.e., the sum of A(i,j)*A(i,j) over all j  with possibly (if scaling!=0) each term j multiplied by (*scaling)(j)
"""
cb_rowip(A::CBMatrix, i::Integer, scaling::Union{<:CBMatrix,Nothing}) = @ccall libcb.cb_matrix_rowip(A.data::Ptr{Cvoid}, i::Cint, (isnothing(scaling) ? C_NULL : scaling.data)::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_colsip(A::CBMatrix)

returns the column vector of the squared Frobenius norm of all columns j of A, i.e., the sum of A(i,j)*A(i,j) over all i for each j
"""
cb_colsip(A::CBMatrix) = CBMatrix(@ccall libcb.cb_matrix_new_colsip(A.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_rowsip(A::CBMatrix)

returns the row vector of the squared Frobenius norm of all rowd i of A, i.e., the sum of A(i,j)*A(i,j) over all j for each i
"""
cb_rowsip(A::CBMatrix) = CBMatrix(@ccall libcb.cb_matrix_new_rowsip(A.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_norm2(A::CBMatrix)

returns the Frobenius norm of A, i.e., the square root of the sum of A(i,j)*A(i,j) over all i,j
"""
cb_norm2(A::CBMatrix) = @ccall libcb.cb_matrix_norm2(A.data::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_normDsquared(A::CBMatrix, d::CBMatrix, atrans::Integer, dinv::Integer)

returns trace(A^TDA)=\|A\|^2_D with D=Diag(d). A may be transposed, D may be inverted but there is no check for division by zero
"""
cb_normDsquared(A::CBMatrix, d::CBMatrix, atrans::Integer, dinv::Integer) = @ccall libcb.cb_matrix_normdsquared(A.data::Ptr{Cvoid}, d.data::Ptr{Cvoid}, atrans::Cint, dinv::Cint)::Cdouble

@doc raw"""
    cb_sumrows(A::CBMatrix)

returns a row vector holding the sum over all rows, i.e., (1 1 ... 1)*A
"""
cb_sumrows(A::CBMatrix) = CBMatrix(@ccall libcb.cb_matrix_new_sumrows(A.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_sumcols(A::CBMatrix)

returns a column vector holding the sum over all columns, i.e., A*(1 1 ... 1)^T
"""
cb_sumcols(A::CBMatrix) = CBMatrix(@ccall libcb.cb_matrix_new_sumcols(A.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_sum(A::CBMatrix)

returns the sum over all elements of A, i.e., (1 1 ... 1)*A*(1 1 ... 1)^T
"""
cb_sum(A::CBMatrix) = @ccall libcb.cb_matrix_sum(A.data::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_house(A::CBMatrix, i::Integer, j::Integer, tol::Real)

returns the Householder vector of size A.rowdim() for the subcolumn A(i:A.rowdim(),j)
"""
cb_house(A::CBMatrix, i::Integer, j::Integer, tol::Real) = CBMatrix(@ccall libcb.cb_matrix_new_house(A.data::Ptr{Cvoid}, i::Cint, j::Cint, tol::Cdouble)::Ptr{Cvoid})

@doc raw"""
    cb_rowhouse(A::CBMatrix, v::CBMatrix, i::Integer, j::Integer)

Housholder pre-multiplication of A with Householder vector v; the first nonzero of v is index i, the multplication is applied to all columns of A with index >=j; always returns 0
"""
cb_rowhouse(A::CBMatrix, v::CBMatrix, i::Integer, j::Integer) = @ccall libcb.cb_matrix_rowhouse(A.data::Ptr{Cvoid}, v.data::Ptr{Cvoid}, i::Cint, j::Cint)::Cint

@doc raw"""
    cb_colhouse(A::CBMatrix, v::CBMatrix, i::Integer, j::Integer)

Housholder post-multiplication of A with Householder vector v; the first nonzero of v is index i, the multplication is applied to all rows of A with index >=j; always returns 0
"""
cb_colhouse(A::CBMatrix, v::CBMatrix, i::Integer, j::Integer) = @ccall libcb.cb_matrix_colhouse(A.data::Ptr{Cvoid}, v.data::Ptr{Cvoid}, i::Cint, j::Cint)::Cint

@doc raw"""
    cb_find(self::CBMatrix, tol::Real = 1e-10)

returns an Indexmatrix ind so that (*this)(ind(i)) 0<=i<ind.dim() runs through all nonzero elements
"""
cb_find(self::CBMatrix, tol::Real = 1e-10) = CBIndexmatrix(@ccall libcb.cb_matrix_new_find(self.data::Ptr{Cvoid}, tol::Cdouble)::Ptr{Cvoid})

@doc raw"""
    cb_find_number(self::CBMatrix, num::Real = 0., tol::Real = 1e-10)

returns an Indexmatrix ind so that (*this)(ind(i)) 0<=i<ind.dim() runs through all elements of value num
"""
cb_find_number(self::CBMatrix, num::Real = 0., tol::Real = 1e-10) = CBIndexmatrix(@ccall libcb.cb_matrix_new_find_number(self.data::Ptr{Cvoid}, num::Cdouble, tol::Cdouble)::Ptr{Cvoid})

@doc raw"""
    Base.:<(A::CBMatrix, B::CBMatrix)

returns a matrix having elements (i,j)=Real(A(i,j)<B(i,j)) for all i,j
"""
Base.:<(A::CBMatrix, B::CBMatrix) = CBMatrix(@ccall libcb.cb_matrix_new_less(A.data::Ptr{Cvoid}, B.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    Base.:>(A::CBMatrix, B::CBMatrix)

returns a matrix having elements (i,j)=Real(A(i,j)>B(i,j)) for all i,j
"""
Base.:>(A::CBMatrix, B::CBMatrix) = CBMatrix(@ccall libcb.cb_matrix_new_greater(A.data::Ptr{Cvoid}, B.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    Base.:<=(A::CBMatrix, B::CBMatrix)

returns a matrix having elements (i,j)=Real(A(i,j)<=B(i,j)) for all i,j
"""
Base.:<=(A::CBMatrix, B::CBMatrix) = CBMatrix(@ccall libcb.cb_matrix_new_lessequal(A.data::Ptr{Cvoid}, B.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    Base.:>=(A::CBMatrix, B::CBMatrix)

returns a matrix having elements (i,j)=Real(A(i,j)>=B(i,j)) for all i,j
"""
Base.:>=(A::CBMatrix, B::CBMatrix) = CBMatrix(@ccall libcb.cb_matrix_new_greaterequal(A.data::Ptr{Cvoid}, B.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    Base.:(==)(A::CBMatrix, B::CBMatrix)

returns a matrix having elements (i,j)=Real(A(i,j)==B(i,j)) for all i,j
"""
Base.:(==)(A::CBMatrix, B::CBMatrix) = CBMatrix(@ccall libcb.cb_matrix_new_equal(A.data::Ptr{Cvoid}, B.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    Base.:<(A::CBMatrix, d::Real)

returns a matrix having elements (i,j)=Real(A(i,j)<d) for all i,j
"""
Base.:<(A::CBMatrix, d::Real) = CBMatrix(@ccall libcb.cb_matrix_new_less2(A.data::Ptr{Cvoid}, d::Cdouble)::Ptr{Cvoid})

@doc raw"""
    Base.:>(A::CBMatrix, d::Real)

returns a matrix having elements (i,j)=Real(A(i,j)>d) for all i,j
"""
Base.:>(A::CBMatrix, d::Real) = CBMatrix(@ccall libcb.cb_matrix_new_greater2(A.data::Ptr{Cvoid}, d::Cdouble)::Ptr{Cvoid})

@doc raw"""
    Base.:<=(A::CBMatrix, d::Real)

returns a matrix having elements (i,j)=Real(A(i,j)<=d) for all i,j
"""
Base.:<=(A::CBMatrix, d::Real) = CBMatrix(@ccall libcb.cb_matrix_new_lessequal2(A.data::Ptr{Cvoid}, d::Cdouble)::Ptr{Cvoid})

@doc raw"""
    Base.:>=(A::CBMatrix, d::Real)

returns a matrix having elements (i,j)=Real(A(i,j)>=d) for all i,j
"""
Base.:>=(A::CBMatrix, d::Real) = CBMatrix(@ccall libcb.cb_matrix_new_greaterequal2(A.data::Ptr{Cvoid}, d::Cdouble)::Ptr{Cvoid})

@doc raw"""
    Base.:(==)(A::CBMatrix, d::Real)

returns a matrix having elements (i,j)=Real(A(i,j)==d) for all i,j
"""
Base.:(==)(A::CBMatrix, d::Real) = CBMatrix(@ccall libcb.cb_matrix_new_equal2(A.data::Ptr{Cvoid}, d::Cdouble)::Ptr{Cvoid})

@doc raw"""
    Base.:<(d::Real, A::CBMatrix)

returns a matrix having elements (i,j)=Real(d<A(i,j)) for all i,j
"""
Base.:<(d::Real, A::CBMatrix) = CBMatrix(@ccall libcb.cb_matrix_new_less3(d::Cdouble, A.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    Base.:>(d::Real, A::CBMatrix)

returns a matrix having elements (i,j)=Real(d>A(i,j)) for all i,j
"""
Base.:>(d::Real, A::CBMatrix) = CBMatrix(@ccall libcb.cb_matrix_new_greater3(d::Cdouble, A.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    Base.:<=(d::Real, A::CBMatrix)

returns a matrix having elements (i,j)=Real(d<=A(i,j)) for all i,j
"""
Base.:<=(d::Real, A::CBMatrix) = CBMatrix(@ccall libcb.cb_matrix_new_lessequal3(d::Cdouble, A.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    Base.:>=(d::Real, A::CBMatrix)

returns a matrix having elements (i,j)=Real(d>=A(i,j)) for all i,j
"""
Base.:>=(d::Real, A::CBMatrix) = CBMatrix(@ccall libcb.cb_matrix_new_greaterequal3(d::Cdouble, A.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    Base.:(==)(d::Real, A::CBMatrix)

returns a matrix having elements (i,j)=Real(d==A(i,j)) for all i,j
"""
Base.:(==)(d::Real, A::CBMatrix) = CBMatrix(@ccall libcb.cb_matrix_new_equal3(d::Cdouble, A.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_equal(A::CBMatrix, B::CBMatrix)

returns true if both matrices have the same size and the same elements
"""
cb_equal(A::CBMatrix, B::CBMatrix) = Bool(@ccall libcb.cb_matrix_equal(A.data::Ptr{Cvoid}, B.data::Ptr{Cvoid})::Cint)

@doc raw"""
    cb_minrows(A::CBMatrix)

returns a row vector holding in each column the minimum over all rows in this column
"""
cb_minrows(A::CBMatrix) = CBMatrix(@ccall libcb.cb_matrix_new_minrows(A.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_mincols(A::CBMatrix)

returns a column vector holding in each row the minimum over all columns in this row
"""
cb_mincols(A::CBMatrix) = CBMatrix(@ccall libcb.cb_matrix_new_mincols(A.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_min(A::CBMatrix, iindex::Union{<:AbstractVector{Cint},Nothing}, jindex::Union{<:AbstractVector{Cint},Nothing})

returns the minimum value over all elements of the matrix
"""
cb_min(A::CBMatrix, iindex::Union{<:AbstractVector{Cint},Nothing}, jindex::Union{<:AbstractVector{Cint},Nothing}) = GC.@preserve iindex jindex begin
    (LinearAlgebra.chkstride1(jindex); LinearAlgebra.chkstride1(iindex); @ccall libcb.cb_matrix_min(A.data::Ptr{Cvoid}, iindex::Ptr{Cint}, jindex::Ptr{Cint})::Cdouble)
end

@doc raw"""
    cb_maxrows(A::CBMatrix)

returns a row vector holding in each column the maximum over all rows in this column
"""
cb_maxrows(A::CBMatrix) = CBMatrix(@ccall libcb.cb_matrix_new_maxrows(A.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_maxcols(A::CBMatrix)

returns a column vector holding in each row the maximum over all columns in this row
"""
cb_maxcols(A::CBMatrix) = CBMatrix(@ccall libcb.cb_matrix_new_maxcols(A.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_max(A::CBMatrix, iindex::Union{<:AbstractVector{Cint},Nothing}, jindex::Union{<:AbstractVector{Cint},Nothing})

returns the maximum value over all elements of the matrix
"""
cb_max(A::CBMatrix, iindex::Union{<:AbstractVector{Cint},Nothing}, jindex::Union{<:AbstractVector{Cint},Nothing}) = GC.@preserve iindex jindex begin
    (LinearAlgebra.chkstride1(jindex); LinearAlgebra.chkstride1(iindex); @ccall libcb.cb_matrix_max(A.data::Ptr{Cvoid}, iindex::Ptr{Cint}, jindex::Ptr{Cint})::Cdouble)
end

@doc raw"""
    cb_sortindex(vec::CBMatrix, nondecreasing::Bool)

returns an Indexmatrix ind so that vec(ind(0))<=vec(ind(1))<=...<=vec(ind(vec.dim()-1)) (vec may be rectangular, set nondecreasing=false for opposite order)
"""
cb_sortindex(vec::CBMatrix, nondecreasing::Bool) = CBIndexmatrix(@ccall libcb.cb_matrix_new_sortindex(vec.data::Ptr{Cvoid}, nondecreasing::Cint)::Ptr{Cvoid})

@doc raw"""
    cb_sortindex(vec::CBMatrix, ind::CBIndexmatrix, nondecreasing::Bool)

sets ind so that vec(ind(0))<=vec(ind(1))<=...<=vec(ind(vec.dim()-1)) (vec may be rectangular, set nondecreasing=false for opposite order)
"""
cb_sortindex(vec::CBMatrix, ind::CBIndexmatrix, nondecreasing::Bool) = @ccall libcb.cb_matrix_sortindex(vec.data::Ptr{Cvoid}, ind.data::Ptr{Cvoid}, nondecreasing::Cint)::Cvoid

@doc raw"""
    cb_display(self::CBMatrix, precision::Integer = 0, width::Integer = 0, screenwidth::Integer = 0)

* @brief displays a matrix in a pretty way for bounded screen widths; for variables of value zero default values are used.
      
"""
cb_display(self::CBMatrix, precision::Integer = 0, width::Integer = 0, screenwidth::Integer = 0) = @ccall libcb.cb_matrix_display(self.data::Ptr{Cvoid}, precision::Cint, width::Cint, screenwidth::Cint)::Cvoid

@doc raw"""
    cb_mfile_output(self::CBMatrix, precision::Integer = 16, width::Integer = 0)

* @brief outputs a matrix A in the format "[ A(0,1) ... A(0,nc-1)\n ... A(nr-1,nc-1)];\n" so that it can be read e.g. by octave as an m-file
     
"""
cb_mfile_output(self::CBMatrix, precision::Integer = 16, width::Integer = 0) = @ccall libcb.cb_matrix_mfile_output(self.data::Ptr{Cvoid}, precision::Cint, width::Cint)::Cvoid

