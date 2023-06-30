@doc raw"""
    CBSparsemat()

empty matrix
"""
CBSparsemat() = CBSparsemat(@ccall libcb.cb_sparsemat_new()::Ptr{Cvoid})

@doc raw"""
    CBSparsemat(A::CBSparsemat, d::Real = 1.)

copy constructor, *this=d*A, abs(values)<tol are removed from the support
"""
CBSparsemat(A::CBSparsemat, d::Real = 1.) = CBSparsemat(@ccall libcb.cb_sparsemat_new2(A.data::Ptr{Cvoid}, d::Cdouble)::Ptr{Cvoid})

@doc raw"""
    CBSparsemat(nr::Integer, nc::Integer)

initialize to zero-matrix of size nr*nc
"""
CBSparsemat(nr::Integer, nc::Integer) = CBSparsemat(@ccall libcb.cb_sparsemat_new3(nr::Cint, nc::Cint)::Ptr{Cvoid})

@doc raw"""
    CBSparsemat(nr::Integer, nc::Integer, nz::Integer, ini::Union{<:AbstractVector{Cint},Nothing}, inj::Union{<:AbstractVector{Cint},Nothing}, va::Union{<:AbstractVector{Cdouble},Nothing})

initialize to size nr*nc and nz nonzeros so that this(ini[i],inj[i])=val[i] for i=0,..,nz-1; multiple elements are summed up.
"""
CBSparsemat(nr::Integer, nc::Integer, nz::Integer, ini::Union{<:AbstractVector{Cint},Nothing}, inj::Union{<:AbstractVector{Cint},Nothing}, va::Union{<:AbstractVector{Cdouble},Nothing}) = GC.@preserve ini inj va begin
    (LinearAlgebra.chkstride1(va); LinearAlgebra.chkstride1(inj); LinearAlgebra.chkstride1(ini); CBSparsemat(@ccall libcb.cb_sparsemat_new4(nr::Cint, nc::Cint, nz::Cint, ini::Ptr{Cint}, inj::Ptr{Cint}, va::Ptr{Cdouble})::Ptr{Cvoid}))
end

@doc raw"""
    CBSparsemat(nr::Integer, nc::Integer, nz::Integer, ini::CBIndexmatrix, inj::CBIndexmatrix, va::CBMatrix)

initialize to size nr*nc and nz nonzeros so that this(ini(i),inj(i))=val(i) for i=0,..,nz-1; multiple elements are summed up.
"""
CBSparsemat(nr::Integer, nc::Integer, nz::Integer, ini::CBIndexmatrix, inj::CBIndexmatrix, va::CBMatrix) = CBSparsemat(@ccall libcb.cb_sparsemat_new5(nr::Cint, nc::Cint, nz::Cint, ini.data::Ptr{Cvoid}, inj.data::Ptr{Cvoid}, va.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_init!(self::CBSparsemat, A::CBSparsemat, d::Real = 1.)

initialize to *this=A*d
"""
cb_init!(self::CBSparsemat, A::CBSparsemat, d::Real = 1.) = (@ccall libcb.cb_sparsemat_init(self.data::Ptr{Cvoid}, A.data::Ptr{Cvoid}, d::Cdouble)::Ptr{Cvoid}; return self)

@doc raw"""
    cb_init!(self::CBSparsemat, A::CBMatrix, d::Real = 1.)

initialize to *this=A*d, abs(values)<tol are removed from the support
"""
cb_init!(self::CBSparsemat, A::CBMatrix, d::Real = 1.) = (@ccall libcb.cb_sparsemat_init2(self.data::Ptr{Cvoid}, A.data::Ptr{Cvoid}, d::Cdouble)::Ptr{Cvoid}; return self)

@doc raw"""
    cb_init!(self::CBSparsemat, A::CBIndexmatrix, d::Real = 1.)

initialize to *this=A*d, zeros are removed from the support
"""
cb_init!(self::CBSparsemat, A::CBIndexmatrix, d::Real = 1.) = (@ccall libcb.cb_sparsemat_init3(self.data::Ptr{Cvoid}, A.data::Ptr{Cvoid}, d::Cdouble)::Ptr{Cvoid}; return self)

@doc raw"""
    cb_init!(self::CBSparsemat, param0::CBSymmatrix, d::Real = 1.)

initialize to *this=A*d, abs(values)<tol are removed from the support
"""
cb_init!(self::CBSparsemat, param0::CBSymmatrix, d::Real = 1.) = (@ccall libcb.cb_sparsemat_init4(self.data::Ptr{Cvoid}, param0.data::Ptr{Cvoid}, d::Cdouble)::Ptr{Cvoid}; return self)

@doc raw"""
    cb_init!(self::CBSparsemat, param0::CBSparsesym, d::Real = 1.)

initialize to *this=A*d
"""
cb_init!(self::CBSparsemat, param0::CBSparsesym, d::Real = 1.) = (@ccall libcb.cb_sparsemat_init5(self.data::Ptr{Cvoid}, param0.data::Ptr{Cvoid}, d::Cdouble)::Ptr{Cvoid}; return self)

@doc raw"""
    cb_init!(self::CBSparsemat, nr::Integer, nc::Integer)

initialize to zero-matrix of size nr*nc
"""
cb_init!(self::CBSparsemat, nr::Integer, nc::Integer) = (@ccall libcb.cb_sparsemat_init6(self.data::Ptr{Cvoid}, nr::Cint, nc::Cint)::Ptr{Cvoid}; return self)

@doc raw"""
    cb_init!(self::CBSparsemat, nr::Integer, nc::Integer, nz::Integer, ini::Union{<:AbstractVector{Cint},Nothing}, inj::Union{<:AbstractVector{Cint},Nothing}, va::Union{<:AbstractVector{Cdouble},Nothing})

initialize to size nr*nc and nz nonzeros so that this(ini[i],inj[i])=val[i] for i=0,..,nz-1; multiple elements are summed up.
"""
cb_init!(self::CBSparsemat, nr::Integer, nc::Integer, nz::Integer, ini::Union{<:AbstractVector{Cint},Nothing}, inj::Union{<:AbstractVector{Cint},Nothing}, va::Union{<:AbstractVector{Cdouble},Nothing}) = GC.@preserve ini inj va begin
    (LinearAlgebra.chkstride1(va); LinearAlgebra.chkstride1(inj); LinearAlgebra.chkstride1(ini); @ccall libcb.cb_sparsemat_init7(self.data::Ptr{Cvoid}, nr::Cint, nc::Cint, nz::Cint, ini::Ptr{Cint}, inj::Ptr{Cint}, va::Ptr{Cdouble})::Ptr{Cvoid}; return self)
end

@doc raw"""
    cb_init!(self::CBSparsemat, nr::Integer, nc::Integer, nz::Integer, ini::CBIndexmatrix, inj::CBIndexmatrix, va::CBMatrix)

initialize to size nr*nc and nz nonzeros so that this(ini(i),inj(i))=val(i) for i=0,..,nz-1; multiple elements are summed up.
"""
cb_init!(self::CBSparsemat, nr::Integer, nc::Integer, nz::Integer, ini::CBIndexmatrix, inj::CBIndexmatrix, va::CBMatrix) = (@ccall libcb.cb_sparsemat_init8(self.data::Ptr{Cvoid}, nr::Cint, nc::Cint, nz::Cint, ini.data::Ptr{Cvoid}, inj.data::Ptr{Cvoid}, va.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

@doc raw"""
    cb_set_tol!(self::CBSparsemat, t::Real)

set tolerance for recognizing zero values to t
"""
cb_set_tol!(self::CBSparsemat, t::Real) = @ccall libcb.cb_sparsemat_set_tol(self.data::Ptr{Cvoid}, t::Cdouble)::Cvoid

@doc raw"""
    CBSparsemat(A::CBMatrix, d::Real = 1.)

initialize to *this=d*A, abs(values)<tol are removed from the support
"""
CBSparsemat(A::CBMatrix, d::Real = 1.) = CBSparsemat(@ccall libcb.cb_sparsemat_new6(A.data::Ptr{Cvoid}, d::Cdouble)::Ptr{Cvoid})

@doc raw"""
    CBSparsemat(A::CBIndexmatrix, d::Real = 1.)

initialize to *this=d*A, zeros are removed from the support
"""
CBSparsemat(A::CBIndexmatrix, d::Real = 1.) = CBSparsemat(@ccall libcb.cb_sparsemat_new7(A.data::Ptr{Cvoid}, d::Cdouble)::Ptr{Cvoid})

@doc raw"""
    CBSparsemat(A::CBSymmatrix, d::Real = 1.)

initialize to *this=d*A, abs(values)<tol are removed from the support
"""
CBSparsemat(A::CBSymmatrix, d::Real = 1.) = CBSparsemat(@ccall libcb.cb_sparsemat_new8(A.data::Ptr{Cvoid}, d::Cdouble)::Ptr{Cvoid})

@doc raw"""
    CBSparsemat(A::CBSparsesym, d::Real = 1.)

initialize to *this=d*A
"""
CBSparsemat(A::CBSparsesym, d::Real = 1.) = CBSparsemat(@ccall libcb.cb_sparsemat_new9(A.data::Ptr{Cvoid}, d::Cdouble)::Ptr{Cvoid})

@doc raw"""
    cb_dim2(self::CBSparsemat)

returns the number of rows in _nr and the number of columns in _nc
"""
function cb_dim2(self::CBSparsemat)
    in_nc = Ref{Int}()
    in_nr = Ref{Int}()
    @ccall libcb.cb_sparsemat_dim(self.data::Ptr{Cvoid}, in_nr::Ref{Int}, in_nc::Ref{Int})::Cvoid
    return in_nr[], in_nc[]
end

@doc raw"""
    cb_dim(self::CBSparsemat)

returns the dimension rows * columns when the matrix is regarded as a vector
"""
cb_dim(self::CBSparsemat) = @ccall libcb.cb_sparsemat_dim2(self.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_rowdim(self::CBSparsemat)

returns the row dimension
"""
cb_rowdim(self::CBSparsemat) = @ccall libcb.cb_sparsemat_rowdim(self.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_coldim(self::CBSparsemat)

returns the column dimension
"""
cb_coldim(self::CBSparsemat) = @ccall libcb.cb_sparsemat_coldim(self.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_nonzeros(self::CBSparsemat)

returns the number of nonzeros
"""
cb_nonzeros(self::CBSparsemat) = @ccall libcb.cb_sparsemat_nonzeros(self.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_col_nonzeros(self::CBSparsemat, i::Integer, startind::Union{<:AbstractVector{Cint},Nothing} = nothing)

returns the number of nonzeros in column i; if nonzeros>0 and startind!=0 then the index of the first nonzero in colindex/colval is stored there
"""
cb_col_nonzeros(self::CBSparsemat, i::Integer, startind::Union{<:AbstractVector{Cint},Nothing} = nothing) = GC.@preserve startind begin
    (LinearAlgebra.chkstride1(startind); @ccall libcb.cb_sparsemat_col_nonzeros(self.data::Ptr{Cvoid}, i::Cint, startind::Ptr{Cint})::Cint)
end

@doc raw"""
    cb_row_nonzeros(self::CBSparsemat, i::Integer, startind::Union{<:AbstractVector{Cint},Nothing} = nothing)

returns the number of nonzeros in row i; if nonzeros>0 and startind!=0 then the index of the first nonzero in rowindex/rowval is stored there
"""
cb_row_nonzeros(self::CBSparsemat, i::Integer, startind::Union{<:AbstractVector{Cint},Nothing} = nothing) = GC.@preserve startind begin
    (LinearAlgebra.chkstride1(startind); @ccall libcb.cb_sparsemat_row_nonzeros(self.data::Ptr{Cvoid}, i::Cint, startind::Ptr{Cint})::Cint)
end

@doc raw"""
    Base.getindex(self::CBSparsemat, i::Integer, j::Integer)

returns value of element (i,j) of the matrix (rowindex i, columnindex j)
"""
Base.getindex(self::CBSparsemat, i::Integer, j::Integer) = @ccall libcb.cb_sparsemat_get(self.data::Ptr{Cvoid}, i::Cint, j::Cint)::Cdouble

@doc raw"""
    Base.getindex(self::CBSparsemat, i::Integer)

returns value of element (i) of the matrix if regarded as vector of stacked columns [element (i%rowdim, i/rowdim)]
"""
Base.getindex(self::CBSparsemat, i::Integer) = @ccall libcb.cb_sparsemat_get2(self.data::Ptr{Cvoid}, i::Cint)::Cdouble

@doc raw"""
    cb_col(self::CBSparsemat, i::Integer)

returns column i copied to a new sparse matrix
"""
cb_col(self::CBSparsemat, i::Integer) = CBSparsemat(@ccall libcb.cb_sparsemat_new_col(self.data::Ptr{Cvoid}, i::Cint)::Ptr{Cvoid})

@doc raw"""
    cb_row(self::CBSparsemat, i::Integer)

returns row i copied to a new sparse matrix
"""
cb_row(self::CBSparsemat, i::Integer) = CBSparsemat(@ccall libcb.cb_sparsemat_new_row(self.data::Ptr{Cvoid}, i::Cint)::Ptr{Cvoid})

@doc raw"""
    cb_cols(self::CBSparsemat, ind::CBIndexmatrix)

returns a sparse matrix of size this->rowdim() x vec.dim(), with column i a copy of column vec(i) of *this
"""
cb_cols(self::CBSparsemat, ind::CBIndexmatrix) = CBSparsemat(@ccall libcb.cb_sparsemat_new_cols(self.data::Ptr{Cvoid}, ind.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_rows(self::CBSparsemat, ind::CBIndexmatrix)

returns a sparse matrix of size vec.dim() x this->rowdim(), with row i a copy of row vec(i) of *this
"""
cb_rows(self::CBSparsemat, ind::CBIndexmatrix) = CBSparsemat(@ccall libcb.cb_sparsemat_new_rows(self.data::Ptr{Cvoid}, ind.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_delete_rows!(self::CBSparsemat, ind::CBIndexmatrix)

all rows indexed by vector ind are deleted, no row should appear twice in ind, remaining rows are moved up keeping their order, returns *this
"""
cb_delete_rows!(self::CBSparsemat, ind::CBIndexmatrix) = (@ccall libcb.cb_sparsemat_delete_rows(self.data::Ptr{Cvoid}, ind.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

@doc raw"""
    cb_delete_cols!(self::CBSparsemat, ind::CBIndexmatrix)

all colmuns indexed by vector ind are deleted, no column should appear twice in ind, remaining columns are moved up keeping their order, returns *this
"""
cb_delete_cols!(self::CBSparsemat, ind::CBIndexmatrix) = (@ccall libcb.cb_sparsemat_delete_cols(self.data::Ptr{Cvoid}, ind.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

@doc raw"""
    cb_insert_row!(self::CBSparsemat, i::Integer, v::CBSparsemat)

insert the row vector v before row i, 0<=i<= row dimension, for i==row dimension the row is appended below; appending to a 0x0 matrix is allowed, returns *this
"""
cb_insert_row!(self::CBSparsemat, i::Integer, v::CBSparsemat) = (@ccall libcb.cb_sparsemat_insert_row(self.data::Ptr{Cvoid}, i::Cint, v.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

@doc raw"""
    cb_insert_col!(self::CBSparsemat, i::Integer, v::CBSparsemat)

insert a column before column i, 0<=i<= column dimension, for i==column dimension the column is appended at the right; appending to a 0x0 matrix is allowed, returns *this
"""
cb_insert_col!(self::CBSparsemat, i::Integer, v::CBSparsemat) = (@ccall libcb.cb_sparsemat_insert_col(self.data::Ptr{Cvoid}, i::Cint, v.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

@doc raw"""
    cb_concat_right!(self::CBSparsemat, A::CBSparsemat)

concats sparse matrix A to the right of *this, A or *this may be the 0x0 matrix initally, returns *this
"""
cb_concat_right!(self::CBSparsemat, A::CBSparsemat) = (@ccall libcb.cb_sparsemat_concat_right(self.data::Ptr{Cvoid}, A.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

@doc raw"""
    cb_concat_below!(self::CBSparsemat, A::CBSparsemat)

concats sparse matrix A to the bottom of *this, A or *this may be the 0x0 matrix initally, returns *this
"""
cb_concat_below!(self::CBSparsemat, A::CBSparsemat) = (@ccall libcb.cb_sparsemat_concat_below(self.data::Ptr{Cvoid}, A.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

@doc raw"""
    cb_get_colinfo(self::CBSparsemat)

returns information on nonzero columns, k by 3, listing: index, %#nonzeros, first index in colindex/colval
"""
cb_get_colinfo(self::CBSparsemat) = (@ccall libcb.cb_sparsemat_get_colinfo(self.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

@doc raw"""
    cb_get_colindex(self::CBSparsemat)

returns the index vector of the column representation holding the row index for each element
"""
cb_get_colindex(self::CBSparsemat) = (@ccall libcb.cb_sparsemat_get_colindex(self.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

@doc raw"""
    cb_get_colval(self::CBSparsemat)

returns the value vector of the column representation holding the value for each element
"""
cb_get_colval(self::CBSparsemat) = (@ccall libcb.cb_sparsemat_get_colval(self.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

@doc raw"""
    cb_get_rowinfo(self::CBSparsemat)

returns information on nonzero rows, k by 3, listing: index, %#nonzeros, first index in rowindex/rowval
"""
cb_get_rowinfo(self::CBSparsemat) = (@ccall libcb.cb_sparsemat_get_rowinfo(self.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

@doc raw"""
    cb_get_rowindex(self::CBSparsemat)

returns the index vector of the row representation holding the column index for each element
"""
cb_get_rowindex(self::CBSparsemat) = (@ccall libcb.cb_sparsemat_get_rowindex(self.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

@doc raw"""
    cb_get_rowval(self::CBSparsemat)

returns the value vector of the row representation holding the value for each element
"""
cb_get_rowval(self::CBSparsemat) = (@ccall libcb.cb_sparsemat_get_rowval(self.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

@doc raw"""
    cb_get_edge_rep(self::CBSparsemat, I::CBIndexmatrix, J::CBIndexmatrix, val::CBMatrix)

stores the nz nonzero values of *this in I,J,val so that this(I(i),J(i))=val(i) for i=0,...,nz-1 and dim(I)=dim(J)=dim(val)=nz (ordered as in row representation)
"""
cb_get_edge_rep(self::CBSparsemat, I::CBIndexmatrix, J::CBIndexmatrix, val::CBMatrix) = @ccall libcb.cb_sparsemat_get_edge_rep(self.data::Ptr{Cvoid}, I.data::Ptr{Cvoid}, J.data::Ptr{Cvoid}, val.data::Ptr{Cvoid})::Cvoid

@doc raw"""
    cb_get_edge(self::CBSparsemat, i::Integer)

stores element i of the get_edge_rep() function (ordered as in row representation); returns 1 if i is out of range, 0 otherwise.
"""
function cb_get_edge(self::CBSparsemat, i::Integer)
    val = Ref{Float64}()
    indj = Ref{Int}()
    indi = Ref{Int}()
    @ccall libcb.cb_sparsemat_get_edge(self.data::Ptr{Cvoid}, i::Cint, indi::Ref{Int}, indj::Ref{Int}, val::Ref{Float64})::Cint
    return indi[], indj[], val[]
end

@doc raw"""
    cb_contains_support(self::CBSparsemat, A::CBSparsemat)

returns 1 if A is of the same dimension and the support of A is contained in the support of *this, 0 otherwise
"""
cb_contains_support(self::CBSparsemat, A::CBSparsemat) = @ccall libcb.cb_sparsemat_contains_support(self.data::Ptr{Cvoid}, A.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_concat_right(A::CBSparsemat, B::CBSparsemat)

returns a new sparse matrix [A, B], i.e., it concats matrices A and B rowwise; A or B may be a 0x0 matrix
"""
cb_concat_right(A::CBSparsemat, B::CBSparsemat) = CBSparsemat(@ccall libcb.cb_sparsemat_new_concat_right(A.data::Ptr{Cvoid}, B.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_concat_below(A::CBSparsemat, B::CBSparsemat)

returns a new sparse matrix [A; B], i.e., it concats matrices A and B columnwise; A or B may be a 0x0 matrix
"""
cb_concat_below(A::CBSparsemat, B::CBSparsemat) = CBSparsemat(@ccall libcb.cb_sparsemat_new_concat_below(A.data::Ptr{Cvoid}, B.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_swap(A::CBSparsemat, B::CBSparsemat)

swap the content of the two sparse matrices A and B (involves no copying)
"""
cb_swap(A::CBSparsemat, B::CBSparsemat) = @ccall libcb.cb_sparsemat_swap(A.data::Ptr{Cvoid}, B.data::Ptr{Cvoid})::Cvoid

@doc raw"""
    cb_xeya!(self::CBSparsemat, A::CBSparsemat, d::Real = 1.)

sets *this=d*A and returns *this
"""
cb_xeya!(self::CBSparsemat, A::CBSparsemat, d::Real = 1.) = (@ccall libcb.cb_sparsemat_xeya(self.data::Ptr{Cvoid}, A.data::Ptr{Cvoid}, d::Cdouble)::Ptr{Cvoid}; return self)

@doc raw"""
    cb_xeya!(self::CBSparsemat, A::CBMatrix, d::Real = 1.)

sets *this=d*A removing abs(values)<tol; returns *this
"""
cb_xeya!(self::CBSparsemat, A::CBMatrix, d::Real = 1.) = (@ccall libcb.cb_sparsemat_xeya2(self.data::Ptr{Cvoid}, A.data::Ptr{Cvoid}, d::Cdouble)::Ptr{Cvoid}; return self)

@doc raw"""
    cb_xeya!(self::CBSparsemat, A::CBIndexmatrix, d::Real = 1.)

sets *this=d*A removing zeros; returns *this
"""
cb_xeya!(self::CBSparsemat, A::CBIndexmatrix, d::Real = 1.) = (@ccall libcb.cb_sparsemat_xeya3(self.data::Ptr{Cvoid}, A.data::Ptr{Cvoid}, d::Cdouble)::Ptr{Cvoid}; return self)

@doc raw"""
    cb_xbpeya(x::CBSparsemat, y::CBSparsemat, alpha::Real, beta::Real, ytrans::Integer)

returns x= alpha*y+beta*x, where y may be transposed (ytrans=1); if beta==0. then x is initialized to the correct size
"""
cb_xbpeya(x::CBSparsemat, y::CBSparsemat, alpha::Real, beta::Real, ytrans::Integer) = (@ccall libcb.cb_sparsemat_xbpeya(x.data::Ptr{Cvoid}, y.data::Ptr{Cvoid}, alpha::Cdouble, beta::Cdouble, ytrans::Cint)::Ptr{Cvoid}; return self)

@doc raw"""
    cb_genmult(A::CBSparsemat, B::CBSparsemat, C::CBMatrix, alpha::Real, beta::Real, atrans::Integer, btrans::Integer)

returns C=beta*C+alpha*A*B, where A and B may be transposed; C must not be equal to A and B; if beta==0. then C is initialized to the correct size
"""
cb_genmult(A::CBSparsemat, B::CBSparsemat, C::CBMatrix, alpha::Real, beta::Real, atrans::Integer, btrans::Integer) = (@ccall libcb.cb_sparsemat_genmult4(A.data::Ptr{Cvoid}, B.data::Ptr{Cvoid}, C.data::Ptr{Cvoid}, alpha::Cdouble, beta::Cdouble, atrans::Cint, btrans::Cint)::Ptr{Cvoid}; return self)

Base.copy!(self::CBSparsemat, A::CBSparsemat) = (@ccall libcb.cb_sparsemat_assign(self.data::Ptr{Cvoid}, A.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

MA.operate!(::typeof(+), self::CBSparsemat, A::CBSparsemat) = (@ccall libcb.cb_sparsemat_plus(self.data::Ptr{Cvoid}, A.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

MA.promote_operation(::typeof(+), ::Type{<:CBSparsemat}, ::Type{<:CBSparsemat}) = CBSparsemat

MA.operate!(::typeof(-), self::CBSparsemat, A::CBSparsemat) = (@ccall libcb.cb_sparsemat_minus(self.data::Ptr{Cvoid}, A.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

MA.promote_operation(::typeof(-), ::Type{<:CBSparsemat}, ::Type{<:CBSparsemat}) = CBSparsemat

Base.:-(self::CBSparsemat) = CBSparsemat(@ccall libcb.cb_sparsemat_new_minus(self.data::Ptr{Cvoid})::Ptr{Cvoid})

MA.operate!(::typeof(*), self::CBSparsemat, d::Real) = (@ccall libcb.cb_sparsemat_times(self.data::Ptr{Cvoid}, d::Cdouble)::Ptr{Cvoid}; return self)

MA.promote_operation(::typeof(*), ::Type{<:CBSparsemat}, ::Type{<:Real}) = CBSparsemat

MA.operate!(::typeof(/), self::CBSparsemat, d::Real) = (@ccall libcb.cb_sparsemat_divide(self.data::Ptr{Cvoid}, d::Cdouble)::Ptr{Cvoid}; return self)

MA.promote_operation(::typeof(/), ::Type{<:CBSparsemat}, ::Type{<:Real}) = CBSparsemat

@doc raw"""
    Base.copy!(self::CBSparsemat, A::CBMatrix)

sets *this=A removing abs(values)<tol; returns *this
"""
Base.copy!(self::CBSparsemat, A::CBMatrix) = (@ccall libcb.cb_sparsemat_assign2(self.data::Ptr{Cvoid}, A.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

@doc raw"""
    Base.copy!(self::CBSparsemat, A::CBIndexmatrix)

sets *this=A removing zeros; returns *this
"""
Base.copy!(self::CBSparsemat, A::CBIndexmatrix) = (@ccall libcb.cb_sparsemat_assign3(self.data::Ptr{Cvoid}, A.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

@doc raw"""
    cb_transpose!(self::CBSparsemat)

transposes itself (swaps row and column representations, thus cheap)
"""
cb_transpose!(self::CBSparsemat) = (@ccall libcb.cb_sparsemat_transpose(self.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

Base.:*(A::CBSparsemat, B::CBSparsemat) = CBSparsemat(@ccall libcb.cb_sparsemat_new_times(A.data::Ptr{Cvoid}, B.data::Ptr{Cvoid})::Ptr{Cvoid})

Base.:+(A::CBSparsemat, B::CBSparsemat) = CBSparsemat(@ccall libcb.cb_sparsemat_new_plus(A.data::Ptr{Cvoid}, B.data::Ptr{Cvoid})::Ptr{Cvoid})

Base.:-(A::CBSparsemat, B::CBSparsemat) = CBSparsemat(@ccall libcb.cb_sparsemat_new_minus2(A.data::Ptr{Cvoid}, B.data::Ptr{Cvoid})::Ptr{Cvoid})

Base.:*(A::CBSparsemat, d::Real) = CBSparsemat(@ccall libcb.cb_sparsemat_new_times2(A.data::Ptr{Cvoid}, d::Cdouble)::Ptr{Cvoid})

Base.:*(d::Real, A::CBSparsemat) = CBSparsemat(@ccall libcb.cb_sparsemat_new_times3(d::Cdouble, A.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    Base.:/(A::CBSparsemat, d::Real)

ATTENTION: d is NOT checked for 0
"""
Base.:/(A::CBSparsemat, d::Real) = CBSparsemat(@ccall libcb.cb_sparsemat_new_divide(A.data::Ptr{Cvoid}, d::Cdouble)::Ptr{Cvoid})

Base.:*(A::CBSparsemat, B::CBMatrix) = CBMatrix(@ccall libcb.cb_sparsemat_new_times4(A.data::Ptr{Cvoid}, B.data::Ptr{Cvoid})::Ptr{Cvoid})

Base.:*(A::CBMatrix, B::CBSparsemat) = CBMatrix(@ccall libcb.cb_sparsemat_new_times5(A.data::Ptr{Cvoid}, B.data::Ptr{Cvoid})::Ptr{Cvoid})

Base.:+(A::CBSparsemat, B::CBMatrix) = CBMatrix(@ccall libcb.cb_sparsemat_new_plus2(A.data::Ptr{Cvoid}, B.data::Ptr{Cvoid})::Ptr{Cvoid})

Base.:+(A::CBMatrix, B::CBSparsemat) = CBMatrix(@ccall libcb.cb_sparsemat_new_plus3(A.data::Ptr{Cvoid}, B.data::Ptr{Cvoid})::Ptr{Cvoid})

Base.:-(A::CBSparsemat, B::CBMatrix) = CBMatrix(@ccall libcb.cb_sparsemat_new_minus3(A.data::Ptr{Cvoid}, B.data::Ptr{Cvoid})::Ptr{Cvoid})

Base.:-(A::CBMatrix, B::CBSparsemat) = CBMatrix(@ccall libcb.cb_sparsemat_new_minus4(A.data::Ptr{Cvoid}, B.data::Ptr{Cvoid})::Ptr{Cvoid})

cb_transpose(A::CBSparsemat) = CBSparsemat(@ccall libcb.cb_sparsemat_new_transpose(A.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_xeya!(self::CBSparsemat, A::CBSparsesym, d::Real = 1.)

sets *this=A*d, abs(values)<tol are removed from the support, and returns *this,
"""
cb_xeya!(self::CBSparsemat, A::CBSparsesym, d::Real = 1.) = (@ccall libcb.cb_sparsemat_xeya4(self.data::Ptr{Cvoid}, A.data::Ptr{Cvoid}, d::Cdouble)::Ptr{Cvoid}; return self)

@doc raw"""
    Base.copy!(self::CBSparsemat, A::CBSymmatrix)

sets *this=A, abs(values)<tol are removed from the support, returns *this
"""
Base.copy!(self::CBSparsemat, A::CBSymmatrix) = (@ccall libcb.cb_sparsemat_assign4(self.data::Ptr{Cvoid}, A.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

@doc raw"""
    Base.copy!(self::CBSparsemat, A::CBSparsesym)

sets *this=A and returns *this
"""
Base.copy!(self::CBSparsemat, A::CBSparsesym) = (@ccall libcb.cb_sparsemat_assign5(self.data::Ptr{Cvoid}, A.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

@doc raw"""
    cb_genmult(A::CBSymmatrix, B::CBSparsemat, C::CBMatrix, alpha::Real, beta::Real, btrans::Integer)

returns C=beta*C+alpha*A*B, where B may be transposed; if beta==0. then C is initialized to the correct size
"""
cb_genmult(A::CBSymmatrix, B::CBSparsemat, C::CBMatrix, alpha::Real, beta::Real, btrans::Integer) = (@ccall libcb.cb_sparsemat_genmult5(A.data::Ptr{Cvoid}, B.data::Ptr{Cvoid}, C.data::Ptr{Cvoid}, alpha::Cdouble, beta::Cdouble, btrans::Cint)::Ptr{Cvoid}; return self)

@doc raw"""
    cb_genmult(A::CBSparsemat, B::CBSymmatrix, C::CBMatrix, alpha::Real, beta::Real, atrans::Integer)

returns C=beta*C+alpha*A*B, where A may be transposed; if beta==0. then C is initialized to the correct size
"""
cb_genmult(A::CBSparsemat, B::CBSymmatrix, C::CBMatrix, alpha::Real, beta::Real, atrans::Integer) = (@ccall libcb.cb_sparsemat_genmult6(A.data::Ptr{Cvoid}, B.data::Ptr{Cvoid}, C.data::Ptr{Cvoid}, alpha::Cdouble, beta::Cdouble, atrans::Cint)::Ptr{Cvoid}; return self)

@doc raw"""
    cb_genmult(A::CBSparsesym, B::CBSparsemat, C::CBMatrix, alpha::Real, beta::Real, btrans::Integer)

returns C=beta*C+alpha*A*B, where B may be transposed; if beta==0. then C is initialized to the correct size
"""
cb_genmult(A::CBSparsesym, B::CBSparsemat, C::CBMatrix, alpha::Real, beta::Real, btrans::Integer) = (@ccall libcb.cb_sparsemat_genmult7(A.data::Ptr{Cvoid}, B.data::Ptr{Cvoid}, C.data::Ptr{Cvoid}, alpha::Cdouble, beta::Cdouble, btrans::Cint)::Ptr{Cvoid}; return self)

@doc raw"""
    cb_genmult(A::CBSparsemat, B::CBSparsesym, C::CBMatrix, alpha::Real, beta::Real, atrans::Integer)

returns C=beta*C+alpha*A*B, where A may be transposed; if beta==0. then C is initialized to the correct size
"""
cb_genmult(A::CBSparsemat, B::CBSparsesym, C::CBMatrix, alpha::Real, beta::Real, atrans::Integer) = (@ccall libcb.cb_sparsemat_genmult8(A.data::Ptr{Cvoid}, B.data::Ptr{Cvoid}, C.data::Ptr{Cvoid}, alpha::Cdouble, beta::Cdouble, atrans::Cint)::Ptr{Cvoid}; return self)

@doc raw"""
    cb_rankadd(A::CBSparsemat, C::CBSymmatrix, alpha::Real, beta::Real, trans::Integer)

returns C=beta*C+alpha* A*A^T, where A may be transposed; if beta==0. then C is initialized to the correct size
"""
cb_rankadd(A::CBSparsemat, C::CBSymmatrix, alpha::Real, beta::Real, trans::Integer) = (@ccall libcb.cb_sparsemat_rankadd(A.data::Ptr{Cvoid}, C.data::Ptr{Cvoid}, alpha::Cdouble, beta::Cdouble, trans::Cint)::Ptr{Cvoid}; return self)

@doc raw"""
    cb_scaledrankadd(A::CBSparsemat, D::CBMatrix, C::CBSymmatrix, alpha::Real, beta::Real, trans::Integer)

returns C=beta*C+alpha* A*D*A^T, where D is a vector representing a diagonal matrix and A may be transposed; if beta==0. then C is initialized to the correct size
"""
cb_scaledrankadd(A::CBSparsemat, D::CBMatrix, C::CBSymmatrix, alpha::Real, beta::Real, trans::Integer) = (@ccall libcb.cb_sparsemat_scaledrankadd(A.data::Ptr{Cvoid}, D.data::Ptr{Cvoid}, C.data::Ptr{Cvoid}, alpha::Cdouble, beta::Cdouble, trans::Cint)::Ptr{Cvoid}; return self)

@doc raw"""
    cb_rank2add(A::CBSparsemat, B::CBMatrix, C::CBSymmatrix, alpha::Real, beta::Real, trans::Integer)

returns C=beta*C+alpha*(A*B^T+B*A^T)/2 [or for transposed (A^T*B+B^T*A)/2]. If beta==0. then C is initiliazed to the correct size.
"""
cb_rank2add(A::CBSparsemat, B::CBMatrix, C::CBSymmatrix, alpha::Real, beta::Real, trans::Integer) = (@ccall libcb.cb_sparsemat_rank2add(A.data::Ptr{Cvoid}, B.data::Ptr{Cvoid}, C.data::Ptr{Cvoid}, alpha::Cdouble, beta::Cdouble, trans::Cint)::Ptr{Cvoid}; return self)

Base.:*(A::CBSymmatrix, B::CBSparsemat) = CBMatrix(@ccall libcb.cb_sparsemat_new_times6(A.data::Ptr{Cvoid}, B.data::Ptr{Cvoid})::Ptr{Cvoid})

Base.:*(A::CBSparsemat, B::CBSymmatrix) = CBMatrix(@ccall libcb.cb_sparsemat_new_times7(A.data::Ptr{Cvoid}, B.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_abs(A::CBSparsemat)

returns a Sparsmat with elements abs((*this)(i,j)) for all i,j
"""
cb_abs(A::CBSparsemat) = CBSparsemat(@ccall libcb.cb_sparsemat_new_abs(A.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_scale_rows!(self::CBSparsemat, vec::CBMatrix)

scales each row i of (*this) by vec(i), i.e., (*this)=diag(vec)*(*this), and returns (*this)
"""
cb_scale_rows!(self::CBSparsemat, vec::CBMatrix) = (@ccall libcb.cb_sparsemat_scale_rows(self.data::Ptr{Cvoid}, vec.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

@doc raw"""
    cb_scale_cols!(self::CBSparsemat, vec::CBMatrix)

scales each column i of (*this) by vec(i), i.e., (*this)=(*this)*diag(vec), and returns (*this)
"""
cb_scale_cols!(self::CBSparsemat, vec::CBMatrix) = (@ccall libcb.cb_sparsemat_scale_cols(self.data::Ptr{Cvoid}, vec.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

@doc raw"""
    cb_trace(A::CBSparsemat)

returns the sum of the diagonal elements A(i,i) over all i
"""
cb_trace(A::CBSparsemat) = @ccall libcb.cb_sparsemat_trace(A.data::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_ip(A::CBSparsemat, B::CBSparsemat)

returns the usual inner product of A and B, i.e., the sum of A(i,j)*B(i,j) over all i,j
"""
cb_ip(A::CBSparsemat, B::CBSparsemat) = @ccall libcb.cb_sparsemat_ip(A.data::Ptr{Cvoid}, B.data::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_ip(A::CBSparsemat, B::CBMatrix)

returns the usual inner product of A and B, i.e., the sum of A(i,j)*B(i,j) over all i,j
"""
cb_ip(A::CBSparsemat, B::CBMatrix) = @ccall libcb.cb_sparsemat_ip2(A.data::Ptr{Cvoid}, B.data::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_ip(A::CBMatrix, B::CBSparsemat)

returns the usual inner product of A and B, i.e., the sum of A(i,j)*B(i,j) over all i,j
"""
cb_ip(A::CBMatrix, B::CBSparsemat) = @ccall libcb.cb_sparsemat_ip3(A.data::Ptr{Cvoid}, B.data::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_colip(A::CBSparsemat, j::Integer, scaling::Union{<:CBMatrix,Nothing})

returns the squared Frobenius norm of col i of A, i.e., the sum of A(i,j)*A(i,j) over all i with possibly (if scaling!=0) each term i multiplied by (*scaling)(i)
"""
cb_colip(A::CBSparsemat, j::Integer, scaling::Union{<:CBMatrix,Nothing}) = @ccall libcb.cb_sparsemat_colip(A.data::Ptr{Cvoid}, j::Cint, (isnothing(scaling) ? C_NULL : scaling.data)::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_rowip(A::CBSparsemat, i::Integer, scaling::Union{<:CBMatrix,Nothing})

returns the squared Frobenius norm of row i of A, i.e., the sum of A(i,j)*A(i,j) over all j with possibly (if scaling!=0) each term j multiplied by (*scaling)(j)
"""
cb_rowip(A::CBSparsemat, i::Integer, scaling::Union{<:CBMatrix,Nothing}) = @ccall libcb.cb_sparsemat_rowip(A.data::Ptr{Cvoid}, i::Cint, (isnothing(scaling) ? C_NULL : scaling.data)::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_colsip(A::CBSparsemat, scaling::Union{<:CBMatrix,Nothing})

returns the column vector of the squared Frobenius norm of all columns j of A, i.e., the sum of A(i,j)*A(i,j) over all i for each j with possibly (if scaling~=0) each term i multiplied by (*scaling)(i)
"""
cb_colsip(A::CBSparsemat, scaling::Union{<:CBMatrix,Nothing}) = CBSparsemat(@ccall libcb.cb_sparsemat_new_colsip(A.data::Ptr{Cvoid}, (isnothing(scaling) ? C_NULL : scaling.data)::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_rowsip(A::CBSparsemat, scaling::Union{<:CBMatrix,Nothing})

returns the row vector of the squared Frobenius norm of all rows i of A, i.e., the sum of A(i,j)*A(i,j) over all j for each i with possibly (if scaling~=0) each term j multiplied by (*scaling)(j)
"""
cb_rowsip(A::CBSparsemat, scaling::Union{<:CBMatrix,Nothing}) = CBSparsemat(@ccall libcb.cb_sparsemat_new_rowsip(A.data::Ptr{Cvoid}, (isnothing(scaling) ? C_NULL : scaling.data)::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_norm2(A::CBSparsemat)

returns the Frobenius norm of A, i.e., the square root of the sum of A(i,j)*A(i,j) over all i,j
"""
cb_norm2(A::CBSparsemat) = @ccall libcb.cb_sparsemat_norm2(A.data::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_sumrows(A::CBSparsemat)

returns a row vector holding the sum over all rows, i.e., (1 1 ... 1)*A
"""
cb_sumrows(A::CBSparsemat) = CBMatrix(@ccall libcb.cb_sparsemat_new_sumrows(A.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_sumcols(A::CBSparsemat)

returns a column vector holding the sum over all columns, i.e., A*(1 1 ... 1)^T
"""
cb_sumcols(A::CBSparsemat) = CBMatrix(@ccall libcb.cb_sparsemat_new_sumcols(A.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_sum(A::CBSparsemat)

returns the sum over all elements of A, i.e., (1 1 ... 1)*A*(1 1 ... 1)^T
"""
cb_sum(A::CBSparsemat) = @ccall libcb.cb_sparsemat_sum(A.data::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_equal(A::CBSparsemat, B::CBSparsemat, eqtol::Real)

returns 1 if both matrices are identical, 0 otherwise
"""
cb_equal(A::CBSparsemat, B::CBSparsemat, eqtol::Real) = @ccall libcb.cb_sparsemat_equal(A.data::Ptr{Cvoid}, B.data::Ptr{Cvoid}, eqtol::Cdouble)::Cint

@doc raw"""
    cb_display(self::CBSparsemat, precision::Integer = 0, width::Integer = 0, screenwidth::Integer = 0)

* @brief displays a matrix in a pretty way for bounded screen widths; for variables of value zero default values are used.
      
"""
cb_display(self::CBSparsemat, precision::Integer = 0, width::Integer = 0, screenwidth::Integer = 0) = @ccall libcb.cb_sparsemat_display(self.data::Ptr{Cvoid}, precision::Cint, width::Cint, screenwidth::Cint)::Cvoid

