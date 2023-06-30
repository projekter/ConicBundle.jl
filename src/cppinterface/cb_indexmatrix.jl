@doc raw"""
    CBIndexmatrix()

empty matrix
"""
CBIndexmatrix() = CBIndexmatrix(@ccall libcb.cb_indexmatrix_new()::Ptr{Cvoid})

@doc raw"""
    CBIndexmatrix(A::CBIndexmatrix, d::Integer = 1)

copy constructor, *this=d*A
"""
CBIndexmatrix(A::CBIndexmatrix, d::Integer = 1) = CBIndexmatrix(@ccall libcb.cb_indexmatrix_new2(A.data::Ptr{Cvoid}, d::Cint)::Ptr{Cvoid})

@doc raw"""
    CBIndexmatrix(param0::AbstractRange{<:Integer})

generate a column vector holding the indices of this #CH_Matrix_Classes::Range
"""
CBIndexmatrix(param0::AbstractRange{<:Integer}) = CBIndexmatrix(@ccall libcb.cb_indexmatrix_new3(first(param0)::Cint, last(param0)::Cint, step(param0)::Cint)::Ptr{Cvoid})

@doc raw"""
    CBIndexmatrix(nr::Integer, nc::Integer)

* @brief generate a matrix of size nr x nc but WITHOUT initializing the memory

        If initializing the memory externally and CONICBUNDLE_DEBUG is defined, please use
        set_init() via matrix.set_init(true) in order to avoid warnings concerning improper
        initialization
    
"""
CBIndexmatrix(nr::Integer, nc::Integer) = CBIndexmatrix(@ccall libcb.cb_indexmatrix_new4(nr::Cint, nc::Cint)::Ptr{Cvoid})

@doc raw"""
    CBIndexmatrix(nr::Integer, nc::Integer, d::Integer)

generate a matrix of size nr x nc initializing all elements to the value d
"""
CBIndexmatrix(nr::Integer, nc::Integer, d::Integer) = CBIndexmatrix(@ccall libcb.cb_indexmatrix_new5(nr::Cint, nc::Cint, d::Cint)::Ptr{Cvoid})

@doc raw"""
    CBIndexmatrix(nr::Integer, nc::Integer, dp::Union{<:AbstractVector{Cint},Nothing})

generate a matrix of size nr x nc initializing the elements from the (one dimensional) array dp with increment incr
"""
CBIndexmatrix(nr::Integer, nc::Integer, dp::Union{<:AbstractVector{Cint},Nothing}) = GC.@preserve dp begin
    CBIndexmatrix(@ccall libcb.cb_indexmatrix_new6(nr::Cint, nc::Cint, dp::Ptr{Cint}, stride(dp, 1)::Cint)::Ptr{Cvoid})
end

@doc raw"""
    cb_init!(self::CBIndexmatrix, A::CBIndexmatrix, d::Integer = 1)

initialize to *this=A*d
"""
cb_init!(self::CBIndexmatrix, A::CBIndexmatrix, d::Integer = 1) = (@ccall libcb.cb_indexmatrix_init(self.data::Ptr{Cvoid}, A.data::Ptr{Cvoid}, d::Cint)::Ptr{Cvoid}; return self)

@doc raw"""
    cb_init!(self::CBIndexmatrix, param0::AbstractRange{<:Integer})

initialize *this to a column vector holding the indices of #CH_Matrix_Classes::Range
"""
cb_init!(self::CBIndexmatrix, param0::AbstractRange{<:Integer}) = (@ccall libcb.cb_indexmatrix_init2(self.data::Ptr{Cvoid}, first(param0)::Cint, last(param0)::Cint, step(param0)::Cint)::Ptr{Cvoid}; return self)

@doc raw"""
    cb_init!(self::CBIndexmatrix, nr::Integer, nc::Integer, d::Integer)

intialize *this to a matrix of size nr x nc initializing all elements to the value d
"""
cb_init!(self::CBIndexmatrix, nr::Integer, nc::Integer, d::Integer) = (@ccall libcb.cb_indexmatrix_init3(self.data::Ptr{Cvoid}, nr::Cint, nc::Cint, d::Cint)::Ptr{Cvoid}; return self)

@doc raw"""
    cb_init!(self::CBIndexmatrix, nr::Integer, nc::Integer, dp::Union{<:AbstractVector{Cint},Nothing})

generate a matrix of size nr x nc initializing the elements from the (one dimensional) array dp with increment incr
"""
cb_init!(self::CBIndexmatrix, nr::Integer, nc::Integer, dp::Union{<:AbstractVector{Cint},Nothing}) = GC.@preserve dp begin
    (@ccall libcb.cb_indexmatrix_init4(self.data::Ptr{Cvoid}, nr::Cint, nc::Cint, dp::Ptr{Cint}, stride(dp, 1)::Cint)::Ptr{Cvoid}; return self)
end

@doc raw"""
    cb_newsize!(self::CBIndexmatrix, nr::Integer, nc::Integer)

* @brief resize the matrix to nr x nc elements but WITHOUT initializing the memory

        If initializing the memory externally and CONICBUNDLE_DEBUG is defined, please use
        set_init() via matrix.set_init(true) in order to avoid warnings concerning improper
        initialization
    
"""
cb_newsize!(self::CBIndexmatrix, nr::Integer, nc::Integer) = @ccall libcb.cb_indexmatrix_newsize(self.data::Ptr{Cvoid}, nr::Cint, nc::Cint)::Cvoid

@doc raw"""
    CBIndexmatrix(param0::CBMatrix)

copy with rounding
"""
CBIndexmatrix(param0::CBMatrix) = CBIndexmatrix(@ccall libcb.cb_indexmatrix_new7(param0.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    CBIndexmatrix(param0::CBSparsemat)

copy with rounding
"""
CBIndexmatrix(param0::CBSparsemat) = CBIndexmatrix(@ccall libcb.cb_indexmatrix_new9(param0.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_dim2(self::CBIndexmatrix)

returns the number of rows in _nr and the number of columns in _nc
"""
function cb_dim2(self::CBIndexmatrix)
    _nc = Ref{Int}()
    _nr = Ref{Int}()
    @ccall libcb.cb_indexmatrix_dim(self.data::Ptr{Cvoid}, _nr::Ref{Int}, _nc::Ref{Int})::Cvoid
    return _nr[], _nc[]
end

@doc raw"""
    cb_dim(self::CBIndexmatrix)

returns the dimension rows * columns when the matrix is regarded as a vector
"""
cb_dim(self::CBIndexmatrix) = @ccall libcb.cb_indexmatrix_dim2(self.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_rowdim(self::CBIndexmatrix)

returns the row dimension
"""
cb_rowdim(self::CBIndexmatrix) = @ccall libcb.cb_indexmatrix_rowdim(self.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_coldim(self::CBIndexmatrix)

returns the column dimension
"""
cb_coldim(self::CBIndexmatrix) = @ccall libcb.cb_indexmatrix_coldim(self.data::Ptr{Cvoid})::Cint

@doc raw"""
    Base.setindex!(self::CBIndexmatrix, value::Integer, i::Integer, j::Integer)

returns reference to element (i,j) of the matrix (rowindex i, columnindex j)
"""
Base.setindex!(self::CBIndexmatrix, value::Integer, i::Integer, j::Integer) = @ccall libcb.cb_indexmatrix_set(self.data::Ptr{Cvoid}, i::Cint, j::Cint, value::Cint)::Cvoid

@doc raw"""
    Base.setindex!(self::CBIndexmatrix, value::Integer, i::Integer)

returns reference to element (i) of the matrix if regarded as vector of stacked columns [element (i%rowdim, i/rowdim)]
"""
Base.setindex!(self::CBIndexmatrix, value::Integer, i::Integer) = @ccall libcb.cb_indexmatrix_set2(self.data::Ptr{Cvoid}, i::Cint, value::Cint)::Cvoid

@doc raw"""
    Base.getindex(self::CBIndexmatrix, i::Integer, j::Integer)

returns value of element (i,j) of the matrix (rowindex i, columnindex j)
"""
Base.getindex(self::CBIndexmatrix, i::Integer, j::Integer) = @ccall libcb.cb_indexmatrix_get(self.data::Ptr{Cvoid}, i::Cint, j::Cint)::Cint

@doc raw"""
    Base.getindex(self::CBIndexmatrix, i::Integer)

returns value of element (i) of the matrix if regarded as vector of stacked columns [element (i%rowdim, i/rowdim)]
"""
Base.getindex(self::CBIndexmatrix, i::Integer) = @ccall libcb.cb_indexmatrix_get2(self.data::Ptr{Cvoid}, i::Cint)::Cint

@doc raw"""
    Base.getindex(self::CBIndexmatrix, vecrow::CBIndexmatrix, veccol::CBIndexmatrix)

returns a new submatrix as indexed by vecrow and veccol, A(i,j)=(*this)(vecrow(i),veccol(j)) for 0<=i<vecrow.dim(), 0<=j<veccol.dim()
"""
Base.getindex(self::CBIndexmatrix, vecrow::CBIndexmatrix, veccol::CBIndexmatrix) = CBIndexmatrix(@ccall libcb.cb_indexmatrix_new_get(self.data::Ptr{Cvoid}, vecrow.data::Ptr{Cvoid}, veccol.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    Base.getindex(self::CBIndexmatrix, A::CBIndexmatrix)

returns a new matrix B of the same shape as A with B(i,j)=(*this)(A(i),A(j)) for 0<=i<A.rowdim(), 0<=j<A.coldim()
"""
Base.getindex(self::CBIndexmatrix, A::CBIndexmatrix) = CBIndexmatrix(@ccall libcb.cb_indexmatrix_new_get2(self.data::Ptr{Cvoid}, A.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_col(self::CBIndexmatrix, i::Integer)

returns column i copied to a new matrix
"""
cb_col(self::CBIndexmatrix, i::Integer) = CBIndexmatrix(@ccall libcb.cb_indexmatrix_new_col(self.data::Ptr{Cvoid}, i::Cint)::Ptr{Cvoid})

@doc raw"""
    cb_row(self::CBIndexmatrix, i::Integer)

returns row i copied to a new matrix
"""
cb_row(self::CBIndexmatrix, i::Integer) = CBIndexmatrix(@ccall libcb.cb_indexmatrix_new_row(self.data::Ptr{Cvoid}, i::Cint)::Ptr{Cvoid})

@doc raw"""
    cb_cols(self::CBIndexmatrix, vec::CBIndexmatrix)

returns a matrix of size this->rowdim() x vec.dim(), with column i a copy of column vec(i) of *this
"""
cb_cols(self::CBIndexmatrix, vec::CBIndexmatrix) = CBIndexmatrix(@ccall libcb.cb_indexmatrix_new_cols(self.data::Ptr{Cvoid}, vec.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_rows(self::CBIndexmatrix, vec::CBIndexmatrix)

returns a matrix of size vec.dim() x this->rowdim(), with row i a copy of row vec(i) of *this
"""
cb_rows(self::CBIndexmatrix, vec::CBIndexmatrix) = CBIndexmatrix(@ccall libcb.cb_indexmatrix_new_rows(self.data::Ptr{Cvoid}, vec.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_triu!(self::CBIndexmatrix, d::Integer = 0)

keeps upper triangle starting with diagonal d; set (*this)(i,j)=0 for 0<=i<row dimension, 0<=j<min(i+d,column dimension), returns *this
"""
cb_triu!(self::CBIndexmatrix, d::Integer = 0) = (@ccall libcb.cb_indexmatrix_triu(self.data::Ptr{Cvoid}, d::Cint)::Ptr{Cvoid}; return self)

@doc raw"""
    cb_tril!(self::CBIndexmatrix, d::Integer = 0)

keeps lower triangle starting with diagonal d; set (*this)(i,j)=0 for 0<=i<row dimension, max(0,i+d)<=j<column dimension, returns *this
"""
cb_tril!(self::CBIndexmatrix, d::Integer = 0) = (@ccall libcb.cb_indexmatrix_tril(self.data::Ptr{Cvoid}, d::Cint)::Ptr{Cvoid}; return self)

@doc raw"""
    cb_subassign!(self::CBIndexmatrix, vecrow::CBIndexmatrix, veccol::CBIndexmatrix, A::CBIndexmatrix)

assigns A to a submatrix of *this,  (*this)(vecrow(i),veccol(j))=A(i,j) for 0<=i<vecrow.dim(), 0<=j<veccol.dim()
"""
cb_subassign!(self::CBIndexmatrix, vecrow::CBIndexmatrix, veccol::CBIndexmatrix, A::CBIndexmatrix) = (@ccall libcb.cb_indexmatrix_subassign(self.data::Ptr{Cvoid}, vecrow.data::Ptr{Cvoid}, veccol.data::Ptr{Cvoid}, A.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

@doc raw"""
    cb_subassign!(self::CBIndexmatrix, vec::CBIndexmatrix, A::CBIndexmatrix)

assigns vector A to a subvector of *this,  (*this)(vec(i))=A(i) for 0<=i<vec.dim(), *this, vec, and A may be rectangular matrices, their dimesions are not changed, returns *this
"""
cb_subassign!(self::CBIndexmatrix, vec::CBIndexmatrix, A::CBIndexmatrix) = (@ccall libcb.cb_indexmatrix_subassign2(self.data::Ptr{Cvoid}, vec.data::Ptr{Cvoid}, A.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

@doc raw"""
    cb_delete_rows!(self::CBIndexmatrix, ind::CBIndexmatrix, sorted_increasingly::Bool = false)

all rows indexed by vector ind are deleted, no row should appear twice in ind, remaining rows are moved up keeping their order, returns *this
"""
cb_delete_rows!(self::CBIndexmatrix, ind::CBIndexmatrix, sorted_increasingly::Bool = false) = (@ccall libcb.cb_indexmatrix_delete_rows(self.data::Ptr{Cvoid}, ind.data::Ptr{Cvoid}, sorted_increasingly::Cint)::Ptr{Cvoid}; return self)

@doc raw"""
    cb_delete_cols!(self::CBIndexmatrix, ind::CBIndexmatrix, sorted_increasingly::Bool = false)

all colmuns indexed by vector ind are deleted, no column should appear twice in ind, remaining columns are moved up keeping their order, returns *this
"""
cb_delete_cols!(self::CBIndexmatrix, ind::CBIndexmatrix, sorted_increasingly::Bool = false) = (@ccall libcb.cb_indexmatrix_delete_cols(self.data::Ptr{Cvoid}, ind.data::Ptr{Cvoid}, sorted_increasingly::Cint)::Ptr{Cvoid}; return self)

@doc raw"""
    cb_insert_row!(self::CBIndexmatrix, i::Integer, v::CBIndexmatrix)

insert the row vector v before row i, 0<=i<= row dimension, for i==row dimension the row is appended below; appending to a 0x0 matrix is allowed, returns *this
"""
cb_insert_row!(self::CBIndexmatrix, i::Integer, v::CBIndexmatrix) = (@ccall libcb.cb_indexmatrix_insert_row(self.data::Ptr{Cvoid}, i::Cint, v.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

@doc raw"""
    cb_insert_col!(self::CBIndexmatrix, i::Integer, v::CBIndexmatrix)

insert a column before column i, 0<=i<= column dimension, for i==column dimension the column is appended at the right; appending to a 0x0 matrix is allowed, returns *this
"""
cb_insert_col!(self::CBIndexmatrix, i::Integer, v::CBIndexmatrix) = (@ccall libcb.cb_indexmatrix_insert_col(self.data::Ptr{Cvoid}, i::Cint, v.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

@doc raw"""
    cb_reduce_length!(self::CBIndexmatrix, n::Integer)

(*this) is set to a column vector of length min{max{0,n},dim()}; usually used to truncate a vector, returns *this
"""
cb_reduce_length!(self::CBIndexmatrix, n::Integer) = (@ccall libcb.cb_indexmatrix_reduce_length(self.data::Ptr{Cvoid}, n::Cint)::Ptr{Cvoid}; return self)

@doc raw"""
    cb_concat_right!(self::CBIndexmatrix, A::CBIndexmatrix)

concats matrix A to the right of *this, A or *this may be the 0x0 matrix initally, returns *this
"""
cb_concat_right!(self::CBIndexmatrix, A::CBIndexmatrix) = (@ccall libcb.cb_indexmatrix_concat_right(self.data::Ptr{Cvoid}, A.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

@doc raw"""
    cb_concat_below!(self::CBIndexmatrix, A::CBIndexmatrix)

concats matrix A to the bottom of *this, A or *this may be the 0x0 matrix initally, returns *this
"""
cb_concat_below!(self::CBIndexmatrix, A::CBIndexmatrix) = (@ccall libcb.cb_indexmatrix_concat_below(self.data::Ptr{Cvoid}, A.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

@doc raw"""
    cb_concat_below!(self::CBIndexmatrix, d::Integer)

concat value d at the bottom of *this, *this must be a column vector or the 0x0 matrix, returns *this
"""
cb_concat_below!(self::CBIndexmatrix, d::Integer) = (@ccall libcb.cb_indexmatrix_concat_below2(self.data::Ptr{Cvoid}, d::Cint)::Ptr{Cvoid}; return self)

@doc raw"""
    cb_concat_right!(self::CBIndexmatrix, d::Integer)

concat value d at the right of *this, *this must be a row vector or the 0x0 matrix, returns *this
"""
cb_concat_right!(self::CBIndexmatrix, d::Integer) = (@ccall libcb.cb_indexmatrix_concat_right2(self.data::Ptr{Cvoid}, d::Cint)::Ptr{Cvoid}; return self)

@doc raw"""
    cb_enlarge_right!(self::CBIndexmatrix, addnc::Integer)

enlarge the matrix by addnc>=0 columns without intializaton of the new columns, returns *this (marked as not initialized if nr>0)
"""
cb_enlarge_right!(self::CBIndexmatrix, addnc::Integer) = (@ccall libcb.cb_indexmatrix_enlarge_right(self.data::Ptr{Cvoid}, addnc::Cint)::Ptr{Cvoid}; return self)

@doc raw"""
    cb_enlarge_below!(self::CBIndexmatrix, addnr::Integer)

enlarge the matrix by addnr>=0 rows without intializaton of the new rows, returns *this (marked as not initialized if nc>0)
"""
cb_enlarge_below!(self::CBIndexmatrix, addnr::Integer) = (@ccall libcb.cb_indexmatrix_enlarge_below(self.data::Ptr{Cvoid}, addnr::Cint)::Ptr{Cvoid}; return self)

@doc raw"""
    cb_enlarge_right!(self::CBIndexmatrix, addnc::Integer, d::Integer)

enlarge the matrix by addnc>=0 columns intializing the new columns by value d, returns *this
"""
cb_enlarge_right!(self::CBIndexmatrix, addnc::Integer, d::Integer) = (@ccall libcb.cb_indexmatrix_enlarge_right2(self.data::Ptr{Cvoid}, addnc::Cint, d::Cint)::Ptr{Cvoid}; return self)

@doc raw"""
    cb_enlarge_below!(self::CBIndexmatrix, addnr::Integer, d::Integer)

enlarge the matrix by addnr>=0 rows intializing the new rows by value d, returns *this
"""
cb_enlarge_below!(self::CBIndexmatrix, addnr::Integer, d::Integer) = (@ccall libcb.cb_indexmatrix_enlarge_below2(self.data::Ptr{Cvoid}, addnr::Cint, d::Cint)::Ptr{Cvoid}; return self)

@doc raw"""
    cb_enlarge_right!(self::CBIndexmatrix, addnc::Integer, dp::Union{<:AbstractVector{Cint},Nothing}, d::Integer = 1)

enlarge the matrix by addnc>=0 columns intializing the new columns by the values pointed to by dp times d, returns *this
"""
cb_enlarge_right!(self::CBIndexmatrix, addnc::Integer, dp::Union{<:AbstractVector{Cint},Nothing}, d::Integer = 1) = GC.@preserve dp begin
    (LinearAlgebra.chkstride1(dp); @ccall libcb.cb_indexmatrix_enlarge_right3(self.data::Ptr{Cvoid}, addnc::Cint, dp::Ptr{Cint}, d::Cint)::Ptr{Cvoid}; return self)
end

@doc raw"""
    cb_enlarge_below!(self::CBIndexmatrix, addnr::Integer, dp::Union{<:AbstractVector{Cint},Nothing}, d::Integer = 1)

enlarge the matrix by addnr>=0 rows intializing the new rows by the values pointed to by dp times d, returns *this
"""
cb_enlarge_below!(self::CBIndexmatrix, addnr::Integer, dp::Union{<:AbstractVector{Cint},Nothing}, d::Integer = 1) = GC.@preserve dp begin
    (LinearAlgebra.chkstride1(dp); @ccall libcb.cb_indexmatrix_enlarge_below3(self.data::Ptr{Cvoid}, addnr::Cint, dp::Ptr{Cint}, d::Cint)::Ptr{Cvoid}; return self)
end

@doc raw"""
    cb_get_store(self::CBIndexmatrix)

returns the current address of the internal value array; use cautiously!
"""
cb_get_store(self::CBIndexmatrix) = @ccall libcb.cb_indexmatrix_get_store2(self.data::Ptr{Cvoid})::Ptr{Cint}

@doc raw"""
    cb_diag(A::CBIndexmatrix)

returns a column vector v consisting of the elements v(i)=A(i,i), 0<=i<min(row dimension,column dimension)
"""
cb_diag(A::CBIndexmatrix) = CBIndexmatrix(@ccall libcb.cb_indexmatrix_new_diag(A.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_triu(A::CBIndexmatrix, d::Integer)

retuns a matrix that keeps the upper triangle of A starting with diagonal d, i.e., (i,j)=A(i,j) for 0<=i<row dimension, max(0,i+d)<=j<column dimension, and sets (i,j)=0 otherwise
"""
cb_triu(A::CBIndexmatrix, d::Integer) = CBIndexmatrix(@ccall libcb.cb_indexmatrix_new_triu(A.data::Ptr{Cvoid}, d::Cint)::Ptr{Cvoid})

@doc raw"""
    cb_tril(A::CBIndexmatrix, d::Integer)

retuns a matrix that keeps the lower triangle of A starting with diagonal d, i.e., (i,j)=A(i,j) for 0<=i<row dimension, 0<=j<min(i+d+1,column dimension), and sets (i,j)=0 otherwise
"""
cb_tril(A::CBIndexmatrix, d::Integer) = CBIndexmatrix(@ccall libcb.cb_indexmatrix_new_tril(A.data::Ptr{Cvoid}, d::Cint)::Ptr{Cvoid})

@doc raw"""
    cb_concat_right(A::CBIndexmatrix, B::CBIndexmatrix)

returns a new matrix [A, B], i.e., it concats matrices A and B rowwise; A or B may be a 0x0 matrix
"""
cb_concat_right(A::CBIndexmatrix, B::CBIndexmatrix) = CBIndexmatrix(@ccall libcb.cb_indexmatrix_new_concat_right(A.data::Ptr{Cvoid}, B.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_concat_below(A::CBIndexmatrix, B::CBIndexmatrix)

returns the matrix [A; B], i.e., it concats matrices A and B columnwise; A or B may be a 0x0 matrix
"""
cb_concat_below(A::CBIndexmatrix, B::CBIndexmatrix) = CBIndexmatrix(@ccall libcb.cb_indexmatrix_new_concat_below(A.data::Ptr{Cvoid}, B.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_swap(A::CBIndexmatrix, B::CBIndexmatrix)

swap the content of the two matrices A and B (involves no copying)
"""
cb_swap(A::CBIndexmatrix, B::CBIndexmatrix) = @ccall libcb.cb_indexmatrix_swap(A.data::Ptr{Cvoid}, B.data::Ptr{Cvoid})::Cvoid

@doc raw"""
    cb_xeya!(self::CBIndexmatrix, A::CBIndexmatrix, d::Integer = 1)

sets *this=d*A and returns *this
"""
cb_xeya!(self::CBIndexmatrix, A::CBIndexmatrix, d::Integer = 1) = (@ccall libcb.cb_indexmatrix_xeya(self.data::Ptr{Cvoid}, A.data::Ptr{Cvoid}, d::Cint)::Ptr{Cvoid}; return self)

@doc raw"""
    cb_xpeya!(self::CBIndexmatrix, A::CBIndexmatrix, d::Integer = 1)

sets *this+=d*A and returns *this
"""
cb_xpeya!(self::CBIndexmatrix, A::CBIndexmatrix, d::Integer = 1) = (@ccall libcb.cb_indexmatrix_xpeya(self.data::Ptr{Cvoid}, A.data::Ptr{Cvoid}, d::Cint)::Ptr{Cvoid}; return self)

@doc raw"""
    cb_xbpeya(x::CBIndexmatrix, y::CBIndexmatrix, alpha::Integer, beta::Integer, ytrans::Integer)

returns x= alpha*y+beta*x, where y may be transposed (ytrans=1); if beta==0. then x is initialized to the correct size
"""
cb_xbpeya(x::CBIndexmatrix, y::CBIndexmatrix, alpha::Integer, beta::Integer, ytrans::Integer) = (@ccall libcb.cb_indexmatrix_xbpeya(x.data::Ptr{Cvoid}, y.data::Ptr{Cvoid}, alpha::Cint, beta::Cint, ytrans::Cint)::Ptr{Cvoid}; return self)

@doc raw"""
    cb_xeyapzb(x::CBIndexmatrix, y::CBIndexmatrix, z::CBIndexmatrix, alpha::Integer, beta::Integer)

returns x= alpha*y+beta*z; x is initialized to the correct size, see CH_Matrix_Classes::xeyapzb() for default values of alpha and beta.
"""
cb_xeyapzb(x::CBIndexmatrix, y::CBIndexmatrix, z::CBIndexmatrix, alpha::Integer, beta::Integer) = (@ccall libcb.cb_indexmatrix_xeyapzb(x.data::Ptr{Cvoid}, y.data::Ptr{Cvoid}, z.data::Ptr{Cvoid}, alpha::Cint, beta::Cint)::Ptr{Cvoid}; return self)

@doc raw"""
    cb_genmult(A::CBIndexmatrix, B::CBIndexmatrix, C::CBIndexmatrix, alpha::Integer, beta::Integer, atrans::Integer, btrans::Integer)

returns C=beta*C+alpha*A*B, where A and B may be transposed; C must not be equal to A and B; if beta==0 then C is initialized to the correct size
"""
cb_genmult(A::CBIndexmatrix, B::CBIndexmatrix, C::CBIndexmatrix, alpha::Integer, beta::Integer, atrans::Integer, btrans::Integer) = (@ccall libcb.cb_indexmatrix_genmult(A.data::Ptr{Cvoid}, B.data::Ptr{Cvoid}, C.data::Ptr{Cvoid}, alpha::Cint, beta::Cint, atrans::Cint, btrans::Cint)::Ptr{Cvoid}; return self)

Base.copy!(self::CBIndexmatrix, A::CBIndexmatrix) = (@ccall libcb.cb_indexmatrix_assign(self.data::Ptr{Cvoid}, A.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

MA.operate!(::typeof(*), self::CBIndexmatrix, s::CBIndexmatrix) = (@ccall libcb.cb_indexmatrix_times(self.data::Ptr{Cvoid}, s.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

MA.promote_operation(::typeof(*), ::Type{<:CBIndexmatrix}, ::Type{<:CBIndexmatrix}) = CBIndexmatrix

MA.operate!(::typeof(+), self::CBIndexmatrix, v::CBIndexmatrix) = (@ccall libcb.cb_indexmatrix_plus(self.data::Ptr{Cvoid}, v.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

MA.promote_operation(::typeof(+), ::Type{<:CBIndexmatrix}, ::Type{<:CBIndexmatrix}) = CBIndexmatrix

MA.operate!(::typeof(-), self::CBIndexmatrix, v::CBIndexmatrix) = (@ccall libcb.cb_indexmatrix_minus(self.data::Ptr{Cvoid}, v.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

MA.promote_operation(::typeof(-), ::Type{<:CBIndexmatrix}, ::Type{<:CBIndexmatrix}) = CBIndexmatrix

@doc raw"""
    MA.operate!(::typeof(%), self::CBIndexmatrix, A::CBIndexmatrix)

ATTENTION: this is redefined as the Hadamard product, (*this)(i,j)=(*this)(i,j)*A(i,j) for all i,j
"""
MA.operate!(::typeof(%), self::CBIndexmatrix, A::CBIndexmatrix) = (@ccall libcb.cb_indexmatrix_rem(self.data::Ptr{Cvoid}, A.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

MA.promote_operation(::typeof(%), ::Type{<:CBIndexmatrix}, ::Type{<:CBIndexmatrix}) = CBIndexmatrix

Base.:-(self::CBIndexmatrix) = CBIndexmatrix(@ccall libcb.cb_indexmatrix_new_minus(self.data::Ptr{Cvoid})::Ptr{Cvoid})

MA.operate!(::typeof(*), self::CBIndexmatrix, d::Integer) = (@ccall libcb.cb_indexmatrix_times2(self.data::Ptr{Cvoid}, d::Cint)::Ptr{Cvoid}; return self)

MA.promote_operation(::typeof(*), ::Type{<:CBIndexmatrix}, ::Type{<:Integer}) = CBIndexmatrix

@doc raw"""
    MA.operate!(::typeof(/), self::CBIndexmatrix, d::Integer)

ATTENTION: d is NOT checked for 0
"""
MA.operate!(::typeof(/), self::CBIndexmatrix, d::Integer) = (@ccall libcb.cb_indexmatrix_divide(self.data::Ptr{Cvoid}, d::Cint)::Ptr{Cvoid}; return self)

MA.promote_operation(::typeof(/), ::Type{<:CBIndexmatrix}, ::Type{<:Integer}) = CBIndexmatrix

@doc raw"""
    MA.operate!(::typeof(%), self::CBIndexmatrix, d::Integer)

sets (*this)(i,j)\%=d for all i,j in the modulo meaning of \%, ATTENTION: d is NOT checked for 0
"""
MA.operate!(::typeof(%), self::CBIndexmatrix, d::Integer) = (@ccall libcb.cb_indexmatrix_rem2(self.data::Ptr{Cvoid}, d::Cint)::Ptr{Cvoid}; return self)

MA.promote_operation(::typeof(%), ::Type{<:CBIndexmatrix}, ::Type{<:Integer}) = CBIndexmatrix

@doc raw"""
    MA.operate!(::typeof(+), self::CBIndexmatrix, d::Integer)

sets (*this)(i,j)+=d for all i,j
"""
MA.operate!(::typeof(+), self::CBIndexmatrix, d::Integer) = (@ccall libcb.cb_indexmatrix_plus2(self.data::Ptr{Cvoid}, d::Cint)::Ptr{Cvoid}; return self)

MA.promote_operation(::typeof(+), ::Type{<:CBIndexmatrix}, ::Type{<:Integer}) = CBIndexmatrix

@doc raw"""
    MA.operate!(::typeof(-), self::CBIndexmatrix, d::Integer)

sets (*this)(i,j)-=d for all i,j
"""
MA.operate!(::typeof(-), self::CBIndexmatrix, d::Integer) = (@ccall libcb.cb_indexmatrix_minus2(self.data::Ptr{Cvoid}, d::Cint)::Ptr{Cvoid}; return self)

MA.promote_operation(::typeof(-), ::Type{<:CBIndexmatrix}, ::Type{<:Integer}) = CBIndexmatrix

@doc raw"""
    cb_transpose!(self::CBIndexmatrix)

transposes itself (cheap for vectors, expensive for matrices)
"""
cb_transpose!(self::CBIndexmatrix) = (@ccall libcb.cb_indexmatrix_transpose(self.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

Base.:*(A::CBIndexmatrix, B::CBIndexmatrix) = CBIndexmatrix(@ccall libcb.cb_indexmatrix_new_times(A.data::Ptr{Cvoid}, B.data::Ptr{Cvoid})::Ptr{Cvoid})

Base.:+(A::CBIndexmatrix, B::CBIndexmatrix) = CBIndexmatrix(@ccall libcb.cb_indexmatrix_new_plus(A.data::Ptr{Cvoid}, B.data::Ptr{Cvoid})::Ptr{Cvoid})

Base.:-(A::CBIndexmatrix, B::CBIndexmatrix) = CBIndexmatrix(@ccall libcb.cb_indexmatrix_new_minus2(A.data::Ptr{Cvoid}, B.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    Base.:%(A::CBIndexmatrix, B::CBIndexmatrix)

ATTENTION: this is redefined as the Hadamard product, (*this)(i,j)=(*this)(i,j)*A(i,j) for all i,j
"""
Base.:%(A::CBIndexmatrix, B::CBIndexmatrix) = CBIndexmatrix(@ccall libcb.cb_indexmatrix_new_rem(A.data::Ptr{Cvoid}, B.data::Ptr{Cvoid})::Ptr{Cvoid})

Base.:*(A::CBIndexmatrix, d::Integer) = CBIndexmatrix(@ccall libcb.cb_indexmatrix_new_times2(A.data::Ptr{Cvoid}, d::Cint)::Ptr{Cvoid})

Base.:*(d::Integer, A::CBIndexmatrix) = CBIndexmatrix(@ccall libcb.cb_indexmatrix_new_times3(d::Cint, A.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    Base.:/(A::CBIndexmatrix, d::Integer)

ATTENTION: d is NOT checked for 0
"""
Base.:/(A::CBIndexmatrix, d::Integer) = CBIndexmatrix(@ccall libcb.cb_indexmatrix_new_divide(A.data::Ptr{Cvoid}, d::Cint)::Ptr{Cvoid})

@doc raw"""
    Base.:%(A::CBIndexmatrix, d::Integer)

sets (i,j)=A(i,j)\%d for all i,j in the modulo meaning of \%, ATTENTION: d is NOT checked for 0
"""
Base.:%(A::CBIndexmatrix, d::Integer) = CBIndexmatrix(@ccall libcb.cb_indexmatrix_new_rem2(A.data::Ptr{Cvoid}, d::Cint)::Ptr{Cvoid})

@doc raw"""
    Base.:+(A::CBIndexmatrix, d::Integer)

returns (i,j)=A(i,j)+d for all i,j
"""
Base.:+(A::CBIndexmatrix, d::Integer) = CBIndexmatrix(@ccall libcb.cb_indexmatrix_new_plus2(A.data::Ptr{Cvoid}, d::Cint)::Ptr{Cvoid})

@doc raw"""
    Base.:+(d::Integer, A::CBIndexmatrix)

returns (i,j)=A(i,j)+d for all i,j
"""
Base.:+(d::Integer, A::CBIndexmatrix) = CBIndexmatrix(@ccall libcb.cb_indexmatrix_new_plus3(d::Cint, A.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    Base.:-(A::CBIndexmatrix, d::Integer)

returns (i,j)=A(i,j)-d for all i,j
"""
Base.:-(A::CBIndexmatrix, d::Integer) = CBIndexmatrix(@ccall libcb.cb_indexmatrix_new_minus3(A.data::Ptr{Cvoid}, d::Cint)::Ptr{Cvoid})

@doc raw"""
    Base.:-(d::Integer, A::CBIndexmatrix)

returns (i,j)=d-A(i,j) for all i,j
"""
Base.:-(d::Integer, A::CBIndexmatrix) = CBIndexmatrix(@ccall libcb.cb_indexmatrix_new_minus4(d::Cint, A.data::Ptr{Cvoid})::Ptr{Cvoid})

cb_transpose(A::CBIndexmatrix) = CBIndexmatrix(@ccall libcb.cb_indexmatrix_new_transpose(A.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_rand!(self::CBIndexmatrix, nr::Integer, nc::Integer, lowerb::Integer, upperb::Integer, random_generator::Union{<:CBGB_rand,Nothing} = nothing)

resize *this to an nr x nc matrix and assign to (i,j) a random number uniformly from [lowerb,upperb] for all i,j
"""
cb_rand!(self::CBIndexmatrix, nr::Integer, nc::Integer, lowerb::Integer, upperb::Integer, random_generator::Union{<:CBGB_rand,Nothing} = nothing) = (@ccall libcb.cb_indexmatrix_rand(self.data::Ptr{Cvoid}, nr::Cint, nc::Cint, lowerb::Cint, upperb::Cint, (isnothing(random_generator) ? C_NULL : random_generator.data)::Ptr{Cvoid})::Ptr{Cvoid}; return self)

@doc raw"""
    cb_shuffle!(self::CBIndexmatrix, random_generator::Union{<:CBGB_rand,Nothing} = nothing)

shuffle the elements randomly (does not change dimensions)
"""
cb_shuffle!(self::CBIndexmatrix, random_generator::Union{<:CBGB_rand,Nothing} = nothing) = (@ccall libcb.cb_indexmatrix_shuffle(self.data::Ptr{Cvoid}, (isnothing(random_generator) ? C_NULL : random_generator.data)::Ptr{Cvoid})::Ptr{Cvoid}; return self)

@doc raw"""
    cb_sign!(self::CBIndexmatrix)

using ::sign assign (*this)(i,j)=sign((*this)(i,j)) for all i,j
"""
cb_sign!(self::CBIndexmatrix) = (@ccall libcb.cb_indexmatrix_sign(self.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

@doc raw"""
    cb_abs!(self::CBIndexmatrix)

using ::abs assign (*this)(i,j)=abs((*this)(i,j)) for all i,j
"""
cb_abs!(self::CBIndexmatrix) = (@ccall libcb.cb_indexmatrix_abs(self.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

@doc raw"""
    cb_rand(nr::Integer, nc::Integer, lb::Integer, ub::Integer, random_generator::Union{<:CBGB_rand,Nothing})

return a nr x nc matrix with (i,j) assigned a random number uniformly from [lowerb,upperb] for all i,j
"""
cb_rand(nr::Integer, nc::Integer, lb::Integer, ub::Integer, random_generator::Union{<:CBGB_rand,Nothing}) = CBIndexmatrix(@ccall libcb.cb_indexmatrix_new_rand(nr::Cint, nc::Cint, lb::Cint, ub::Cint, (isnothing(random_generator) ? C_NULL : random_generator.data)::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_sign(A::CBIndexmatrix)

return a matrix of the same size as A with (i,j)=sign(A(i,j)) for all i,j, see also CH_Matrix_Classes::sign()
"""
cb_sign(A::CBIndexmatrix) = CBIndexmatrix(@ccall libcb.cb_indexmatrix_new_sign(A.data::Ptr{Cvoid})::Ptr{Cvoid})

cb_abs(A::CBIndexmatrix) = CBIndexmatrix(@ccall libcb.cb_indexmatrix_new_abs(A.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_trace(A::CBIndexmatrix)

returns the sum of the diagonal elements A(i,i) over all i
"""
cb_trace(A::CBIndexmatrix) = @ccall libcb.cb_indexmatrix_trace(A.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_ip(A::CBIndexmatrix, B::CBIndexmatrix)

returns the usual inner product of A and B, i.e., the sum of A(i,j)*B(i,j) over all i,j
"""
cb_ip(A::CBIndexmatrix, B::CBIndexmatrix) = @ccall libcb.cb_indexmatrix_ip(A.data::Ptr{Cvoid}, B.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_norm2(A::CBIndexmatrix)

returns the Frobenius norm of A, i.e., the square root of the sum of A(i,j)*A(i,j) over all i,j
"""
cb_norm2(A::CBIndexmatrix) = @ccall libcb.cb_indexmatrix_norm2(A.data::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_sumrows(A::CBIndexmatrix)

returns a row vector holding the sum over all rows, i.e., (1 1 ... 1)*A
"""
cb_sumrows(A::CBIndexmatrix) = CBIndexmatrix(@ccall libcb.cb_indexmatrix_new_sumrows(A.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_sumcols(A::CBIndexmatrix)

returns a column vector holding the sum over all columns, i.e., A*(1 1 ... 1)^T
"""
cb_sumcols(A::CBIndexmatrix) = CBIndexmatrix(@ccall libcb.cb_indexmatrix_new_sumcols(A.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_sum(A::CBIndexmatrix)

returns the sum over all elements of A, i.e., (1 1 ... 1)*A*(1 1 ... 1)^T
"""
cb_sum(A::CBIndexmatrix) = @ccall libcb.cb_indexmatrix_sum(A.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_find(self::CBIndexmatrix)

returns an Indexmatrix ind so that (*this)(ind(i)) 0<=i<ind.dim() runs through all nonzero elements
"""
cb_find(self::CBIndexmatrix) = CBIndexmatrix(@ccall libcb.cb_indexmatrix_new_find(self.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_find_number(self::CBIndexmatrix, num::Integer = 0)

returns an Indexmatrix ind so that (*this)(ind(i)) 0<=i<ind.dim() runs through all elements having value num
"""
cb_find_number(self::CBIndexmatrix, num::Integer = 0) = CBIndexmatrix(@ccall libcb.cb_indexmatrix_new_find_number(self.data::Ptr{Cvoid}, num::Cint)::Ptr{Cvoid})

@doc raw"""
    Base.:<(A::CBIndexmatrix, B::CBIndexmatrix)

returns a matrix having elements (i,j)=Integer(A(i,j)<B(i,j)) for all i,j
"""
Base.:<(A::CBIndexmatrix, B::CBIndexmatrix) = CBIndexmatrix(@ccall libcb.cb_indexmatrix_new_less(A.data::Ptr{Cvoid}, B.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    Base.:>(A::CBIndexmatrix, B::CBIndexmatrix)

returns a matrix having elements (i,j)=Integer(A(i,j)>B(i,j)) for all i,j
"""
Base.:>(A::CBIndexmatrix, B::CBIndexmatrix) = CBIndexmatrix(@ccall libcb.cb_indexmatrix_new_greater(A.data::Ptr{Cvoid}, B.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    Base.:<=(A::CBIndexmatrix, B::CBIndexmatrix)

returns a matrix having elements (i,j)=Integer(A(i,j)<=B(i,j)) for all i,j
"""
Base.:<=(A::CBIndexmatrix, B::CBIndexmatrix) = CBIndexmatrix(@ccall libcb.cb_indexmatrix_new_lessequal(A.data::Ptr{Cvoid}, B.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    Base.:>=(A::CBIndexmatrix, B::CBIndexmatrix)

returns a matrix having elements (i,j)=Integer(A(i,j)>=B(i,j)) for all i,j
"""
Base.:>=(A::CBIndexmatrix, B::CBIndexmatrix) = CBIndexmatrix(@ccall libcb.cb_indexmatrix_new_greaterequal(A.data::Ptr{Cvoid}, B.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    Base.:(==)(A::CBIndexmatrix, B::CBIndexmatrix)

returns a matrix having elements (i,j)=Integer(A(i,j)==B(i,j)) for all i,j
"""
Base.:(==)(A::CBIndexmatrix, B::CBIndexmatrix) = CBIndexmatrix(@ccall libcb.cb_indexmatrix_new_equal(A.data::Ptr{Cvoid}, B.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    Base.:<(A::CBIndexmatrix, d::Integer)

returns a matrix having elements (i,j)=Integer(A(i,j)<d) for all i,j
"""
Base.:<(A::CBIndexmatrix, d::Integer) = CBIndexmatrix(@ccall libcb.cb_indexmatrix_new_less2(A.data::Ptr{Cvoid}, d::Cint)::Ptr{Cvoid})

@doc raw"""
    Base.:>(A::CBIndexmatrix, d::Integer)

returns a matrix having elements (i,j)=Integer(A(i,j)>d) for all i,j
"""
Base.:>(A::CBIndexmatrix, d::Integer) = CBIndexmatrix(@ccall libcb.cb_indexmatrix_new_greater2(A.data::Ptr{Cvoid}, d::Cint)::Ptr{Cvoid})

@doc raw"""
    Base.:<=(A::CBIndexmatrix, d::Integer)

returns a matrix having elements (i,j)=Integer(A(i,j)<=d) for all i,j
"""
Base.:<=(A::CBIndexmatrix, d::Integer) = CBIndexmatrix(@ccall libcb.cb_indexmatrix_new_lessequal2(A.data::Ptr{Cvoid}, d::Cint)::Ptr{Cvoid})

@doc raw"""
    Base.:>=(A::CBIndexmatrix, d::Integer)

returns a matrix having elements (i,j)=Integer(A(i,j)>=d) for all i,j
"""
Base.:>=(A::CBIndexmatrix, d::Integer) = CBIndexmatrix(@ccall libcb.cb_indexmatrix_new_greaterequal2(A.data::Ptr{Cvoid}, d::Cint)::Ptr{Cvoid})

@doc raw"""
    Base.:(==)(A::CBIndexmatrix, d::Integer)

returns a matrix having elements (i,j)=Integer(A(i,j)==d) for all i,j
"""
Base.:(==)(A::CBIndexmatrix, d::Integer) = CBIndexmatrix(@ccall libcb.cb_indexmatrix_new_equal2(A.data::Ptr{Cvoid}, d::Cint)::Ptr{Cvoid})

@doc raw"""
    Base.:<(d::Integer, A::CBIndexmatrix)

returns a matrix having elements (i,j)=Integer(d<A(i,j)) for all i,j
"""
Base.:<(d::Integer, A::CBIndexmatrix) = CBIndexmatrix(@ccall libcb.cb_indexmatrix_new_less3(d::Cint, A.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    Base.:>(d::Integer, A::CBIndexmatrix)

returns a matrix having elements (i,j)=Integer(d>A(i,j)) for all i,j
"""
Base.:>(d::Integer, A::CBIndexmatrix) = CBIndexmatrix(@ccall libcb.cb_indexmatrix_new_greater3(d::Cint, A.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    Base.:<=(d::Integer, A::CBIndexmatrix)

returns a matrix having elements (i,j)=Integer(d<=A(i,j)) for all i,j
"""
Base.:<=(d::Integer, A::CBIndexmatrix) = CBIndexmatrix(@ccall libcb.cb_indexmatrix_new_lessequal3(d::Cint, A.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    Base.:>=(d::Integer, A::CBIndexmatrix)

returns a matrix having elements (i,j)=Integer(d>=A(i,j)) for all i,j
"""
Base.:>=(d::Integer, A::CBIndexmatrix) = CBIndexmatrix(@ccall libcb.cb_indexmatrix_new_greaterequal3(d::Cint, A.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    Base.:(==)(d::Integer, A::CBIndexmatrix)

returns a matrix having elements (i,j)=Integer(d==A(i,j)) for all i,j
"""
Base.:(==)(d::Integer, A::CBIndexmatrix) = CBIndexmatrix(@ccall libcb.cb_indexmatrix_new_equal3(d::Cint, A.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_equal(A::CBIndexmatrix, b::CBIndexmatrix)

returns true if both matrices have the same size and the same elements
"""
cb_equal(A::CBIndexmatrix, b::CBIndexmatrix) = Bool(@ccall libcb.cb_indexmatrix_equal(A.data::Ptr{Cvoid}, b.data::Ptr{Cvoid})::Cint)

@doc raw"""
    cb_minrows(A::CBIndexmatrix)

returns a row vector holding in each column the minimum over all rows in this column
"""
cb_minrows(A::CBIndexmatrix) = CBIndexmatrix(@ccall libcb.cb_indexmatrix_new_minrows(A.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_mincols(A::CBIndexmatrix)

returns a column vector holding in each row the minimum over all columns in this row
"""
cb_mincols(A::CBIndexmatrix) = CBIndexmatrix(@ccall libcb.cb_indexmatrix_new_mincols(A.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_min(A::CBIndexmatrix, iindex::Union{<:AbstractVector{Cint},Nothing}, jindex::Union{<:AbstractVector{Cint},Nothing})

returns the minimum value over all elements of the matrix
"""
cb_min(A::CBIndexmatrix, iindex::Union{<:AbstractVector{Cint},Nothing}, jindex::Union{<:AbstractVector{Cint},Nothing}) = GC.@preserve iindex jindex begin
    (LinearAlgebra.chkstride1(jindex); LinearAlgebra.chkstride1(iindex); @ccall libcb.cb_indexmatrix_min(A.data::Ptr{Cvoid}, iindex::Ptr{Cint}, jindex::Ptr{Cint})::Cint)
end

@doc raw"""
    cb_maxrows(A::CBIndexmatrix)

returns a row vector holding in each column the maximum over all rows in this column
"""
cb_maxrows(A::CBIndexmatrix) = CBIndexmatrix(@ccall libcb.cb_indexmatrix_new_maxrows(A.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_maxcols(A::CBIndexmatrix)

returns a column vector holding in each row the maximum over all columns in this row
"""
cb_maxcols(A::CBIndexmatrix) = CBIndexmatrix(@ccall libcb.cb_indexmatrix_new_maxcols(A.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_max(A::CBIndexmatrix, iindex::Union{<:AbstractVector{Cint},Nothing}, jindex::Union{<:AbstractVector{Cint},Nothing})

returns the maximum value over all elements of the matrix
"""
cb_max(A::CBIndexmatrix, iindex::Union{<:AbstractVector{Cint},Nothing}, jindex::Union{<:AbstractVector{Cint},Nothing}) = GC.@preserve iindex jindex begin
    (LinearAlgebra.chkstride1(jindex); LinearAlgebra.chkstride1(iindex); @ccall libcb.cb_indexmatrix_max(A.data::Ptr{Cvoid}, iindex::Ptr{Cint}, jindex::Ptr{Cint})::Cint)
end

@doc raw"""
    cb_sortindex(vec::CBIndexmatrix, nondecreasing::Bool)

returns an Indexmatrix ind so that vec(ind(0))<=vec(ind(1))<=...<=vec(ind(vec.dim()-1)) (vec may be rectangular, set nondecreasing=false for opposite order)
"""
cb_sortindex(vec::CBIndexmatrix, nondecreasing::Bool) = CBIndexmatrix(@ccall libcb.cb_indexmatrix_new_sortindex(vec.data::Ptr{Cvoid}, nondecreasing::Cint)::Ptr{Cvoid})

@doc raw"""
    cb_sortindex(vec::CBIndexmatrix, ind::CBIndexmatrix, nondecreasing::Bool)

sets ind so that vec(ind(0))<=vec(ind(1))<=...<=vec(ind(vec.dim()-1)) (vec may be rectangular, set nondecreasing=false for opposite order)
"""
cb_sortindex(vec::CBIndexmatrix, ind::CBIndexmatrix, nondecreasing::Bool) = @ccall libcb.cb_indexmatrix_sortindex(vec.data::Ptr{Cvoid}, ind.data::Ptr{Cvoid}, nondecreasing::Cint)::Cvoid

@doc raw"""
    cb_display(self::CBIndexmatrix, precision::Integer = 0, width::Integer = 0, screenwidth::Integer = 0)

* @brief displays a matrix in a pretty way for bounded screen widths; for variables of value zero default values are used.
      
"""
cb_display(self::CBIndexmatrix, precision::Integer = 0, width::Integer = 0, screenwidth::Integer = 0) = @ccall libcb.cb_indexmatrix_display(self.data::Ptr{Cvoid}, precision::Cint, width::Cint, screenwidth::Cint)::Cvoid

@doc raw"""
    cb_mfile_output(self::CBIndexmatrix, precision::Integer = 16, width::Integer = 0)

* @brief outputs a matrix A in the format "[ A(0,1) ... A(0,nc-1)\n ... A(nr-1,nc-1)];\n" so that it can be read e.g. by octave as an m-file
     
"""
cb_mfile_output(self::CBIndexmatrix, precision::Integer = 16, width::Integer = 0) = @ccall libcb.cb_indexmatrix_mfile_output(self.data::Ptr{Cvoid}, precision::Cint, width::Cint)::Cvoid

