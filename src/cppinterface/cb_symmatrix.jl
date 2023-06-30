@doc raw"""
    CBSymmatrix()

empty matrix
"""
CBSymmatrix() = CBSymmatrix(@ccall libcb.cb_symmatrix_new()::Ptr{Cvoid})

@doc raw"""
    CBSymmatrix(A::CBSymmatrix, d::Real = 1.)

copy constructor, *this=d*A
"""
CBSymmatrix(A::CBSymmatrix, d::Real = 1.) = CBSymmatrix(@ccall libcb.cb_symmatrix_new2(A.data::Ptr{Cvoid}, d::Cdouble)::Ptr{Cvoid})

@doc raw"""
    CBSymmatrix(nr::Integer)

* @brief generate a matrix of size nr x nr but WITHOUT initializing the memory

        If initializing the memory externally and CONICBUNDLE_DEBUG is defined, please use
        set_init() via matrix.set_init(true) in order to avoid warnings concerning improper
        initialization
    
"""
CBSymmatrix(nr::Integer) = CBSymmatrix(@ccall libcb.cb_symmatrix_new3(nr::Cint)::Ptr{Cvoid})

@doc raw"""
    CBSymmatrix(nr::Integer, d::Real)

generate a matrix of size nr x nr initializing all elements to the value d
"""
CBSymmatrix(nr::Integer, d::Real) = CBSymmatrix(@ccall libcb.cb_symmatrix_new4(nr::Cint, d::Cdouble)::Ptr{Cvoid})

@doc raw"""
    CBSymmatrix(nr::Integer, dp::Union{<:AbstractVector{Cdouble},Nothing})

generate a matrix of size nr x nr initializing the elements from the (one dimensional) array dp, which must have the elements arranged consecutively in internal order
"""
CBSymmatrix(nr::Integer, dp::Union{<:AbstractVector{Cdouble},Nothing}) = GC.@preserve dp begin
    (LinearAlgebra.chkstride1(dp); CBSymmatrix(@ccall libcb.cb_symmatrix_new5(nr::Cint, dp::Ptr{Cdouble})::Ptr{Cvoid}))
end

@doc raw"""
    cb_init!(self::CBSymmatrix, A::CBSymmatrix, d::Real = 1.)

initialize to *this=A*d
"""
cb_init!(self::CBSymmatrix, A::CBSymmatrix, d::Real = 1.) = (@ccall libcb.cb_symmatrix_init(self.data::Ptr{Cvoid}, A.data::Ptr{Cvoid}, d::Cdouble)::Ptr{Cvoid}; return self)

@doc raw"""
    cb_init!(self::CBSymmatrix, A::CBMatrix, d::Real = 1.)

initialize to *this=d*(A+transpose(A))/2.
"""
cb_init!(self::CBSymmatrix, A::CBMatrix, d::Real = 1.) = (@ccall libcb.cb_symmatrix_init2(self.data::Ptr{Cvoid}, A.data::Ptr{Cvoid}, d::Cdouble)::Ptr{Cvoid}; return self)

@doc raw"""
    cb_init!(self::CBSymmatrix, A::CBIndexmatrix, d::Real = 1.)

initialize to *this=d*(A+transpose(A))/2.
"""
cb_init!(self::CBSymmatrix, A::CBIndexmatrix, d::Real = 1.) = (@ccall libcb.cb_symmatrix_init3(self.data::Ptr{Cvoid}, A.data::Ptr{Cvoid}, d::Cdouble)::Ptr{Cvoid}; return self)

@doc raw"""
    cb_init!(self::CBSymmatrix, A::CBSparsesym, d::Real = 1.)

initialize to *this=A*d
"""
cb_init!(self::CBSymmatrix, A::CBSparsesym, d::Real = 1.) = (@ccall libcb.cb_symmatrix_init4(self.data::Ptr{Cvoid}, A.data::Ptr{Cvoid}, d::Cdouble)::Ptr{Cvoid}; return self)

@doc raw"""
    cb_init!(self::CBSymmatrix, nr::Integer, d::Real)

intialize *this to a matrix of size nr x nr initializing all elements to the value d
"""
cb_init!(self::CBSymmatrix, nr::Integer, d::Real) = (@ccall libcb.cb_symmatrix_init5(self.data::Ptr{Cvoid}, nr::Cint, d::Cdouble)::Ptr{Cvoid}; return self)

@doc raw"""
    cb_init!(self::CBSymmatrix, nr::Integer, dp::Union{<:AbstractVector{Cdouble},Nothing})

generate a matrix of size nr x nc initializing the elements from the (one dimensional) array dp which must have the elements arranged consecutively in internal order
"""
cb_init!(self::CBSymmatrix, nr::Integer, dp::Union{<:AbstractVector{Cdouble},Nothing}) = GC.@preserve dp begin
    (LinearAlgebra.chkstride1(dp); @ccall libcb.cb_symmatrix_init6(self.data::Ptr{Cvoid}, nr::Cint, dp::Ptr{Cdouble})::Ptr{Cvoid}; return self)
end

@doc raw"""
    cb_newsize!(self::CBSymmatrix, n::Integer)

* @brief resize the matrix to nr x nr elements but WITHOUT initializing the memory

        If initializing the memory externally and CONICBUNDLE_DEBUG is defined, please use
        set_init() via matrix.set_init(true) in order to avoid warnings concerning improper
        initialization
    
"""
cb_newsize!(self::CBSymmatrix, n::Integer) = @ccall libcb.cb_symmatrix_newsize(self.data::Ptr{Cvoid}, n::Cint)::Cvoid

@doc raw"""
    CBSymmatrix(param0::CBMatrix, d::Real = 1.)

(*this)=d*(A+transpose(A))/2.
"""
CBSymmatrix(param0::CBMatrix, d::Real = 1.) = CBSymmatrix(@ccall libcb.cb_symmatrix_new6(param0.data::Ptr{Cvoid}, d::Cdouble)::Ptr{Cvoid})

@doc raw"""
    CBSymmatrix(param0::CBIndexmatrix, d::Real = 1.)

(*this)=d*(A+transpose(A))/2.
"""
CBSymmatrix(param0::CBIndexmatrix, d::Real = 1.) = CBSymmatrix(@ccall libcb.cb_symmatrix_new7(param0.data::Ptr{Cvoid}, d::Cdouble)::Ptr{Cvoid})

@doc raw"""
    CBSymmatrix(A::CBSparsesym, d::Real = 1.)

(*this)=d*A
"""
CBSymmatrix(A::CBSparsesym, d::Real = 1.) = CBSymmatrix(@ccall libcb.cb_symmatrix_new8(A.data::Ptr{Cvoid}, d::Cdouble)::Ptr{Cvoid})

@doc raw"""
    cb_dim2(self::CBSymmatrix)

returns the number of rows in _nr and _nc
"""
function cb_dim2(self::CBSymmatrix)
    _nc = Ref{Int}()
    _nr = Ref{Int}()
    @ccall libcb.cb_symmatrix_dim(self.data::Ptr{Cvoid}, _nr::Ref{Int}, _nc::Ref{Int})::Cvoid
    return _nr[], _nc[]
end

@doc raw"""
    cb_dim(self::CBSymmatrix)

returns the dimension rows * columns when the matrix is regarded as a vector
"""
cb_dim(self::CBSymmatrix) = @ccall libcb.cb_symmatrix_dim2(self.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_rowdim(self::CBSymmatrix)

returns the row dimension
"""
cb_rowdim(self::CBSymmatrix) = @ccall libcb.cb_symmatrix_rowdim(self.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_coldim(self::CBSymmatrix)

returns the column dimension
"""
cb_coldim(self::CBSymmatrix) = @ccall libcb.cb_symmatrix_coldim(self.data::Ptr{Cvoid})::Cint

@doc raw"""
    Base.setindex!(self::CBSymmatrix, value::Real, i::Integer, j::Integer)

returns reference to element (i,j) of the matrix (rowindex i, columnindex j)
"""
Base.setindex!(self::CBSymmatrix, value::Real, i::Integer, j::Integer) = @ccall libcb.cb_symmatrix_set(self.data::Ptr{Cvoid}, i::Cint, j::Cint, value::Cdouble)::Cvoid

@doc raw"""
    Base.setindex!(self::CBSymmatrix, value::Real, i::Integer)

returns reference to element (i) of the matrix if regarded as vector of stacked columns [element (i%rowdim, i/rowdim)]
"""
Base.setindex!(self::CBSymmatrix, value::Real, i::Integer) = @ccall libcb.cb_symmatrix_set2(self.data::Ptr{Cvoid}, i::Cint, value::Cdouble)::Cvoid

@doc raw"""
    Base.getindex(self::CBSymmatrix, i::Integer, j::Integer)

returns value of element (i,j) of the matrix (rowindex i, columnindex j)
"""
Base.getindex(self::CBSymmatrix, i::Integer, j::Integer) = @ccall libcb.cb_symmatrix_get(self.data::Ptr{Cvoid}, i::Cint, j::Cint)::Cdouble

@doc raw"""
    Base.getindex(self::CBSymmatrix, i::Integer)

returns value of element (i) of the matrix if regarded as vector of stacked columns [element (i%rowdim, i/rowdim)]
"""
Base.getindex(self::CBSymmatrix, i::Integer) = @ccall libcb.cb_symmatrix_get2(self.data::Ptr{Cvoid}, i::Cint)::Cdouble

@doc raw"""
    cb_col(self::CBSymmatrix, i::Integer)

returns column i copied to a new Matrix
"""
cb_col(self::CBSymmatrix, i::Integer) = CBMatrix(@ccall libcb.cb_symmatrix_new_col(self.data::Ptr{Cvoid}, i::Cint)::Ptr{Cvoid})

@doc raw"""
    cb_row(self::CBSymmatrix, i::Integer)

returns row i copied to a new Matrix
"""
cb_row(self::CBSymmatrix, i::Integer) = CBMatrix(@ccall libcb.cb_symmatrix_new_row(self.data::Ptr{Cvoid}, i::Cint)::Ptr{Cvoid})

@doc raw"""
    cb_cols(self::CBSymmatrix, vec::CBIndexmatrix)

returns a matrix of size this->rowdim() x vec.dim(), with column i a copy of column vec(i) of *this
"""
cb_cols(self::CBSymmatrix, vec::CBIndexmatrix) = CBMatrix(@ccall libcb.cb_symmatrix_new_cols(self.data::Ptr{Cvoid}, vec.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_rows(self::CBSymmatrix, vec::CBIndexmatrix)

returns a matrix of size vec.dim() x this->coldim(), with row i a copy of row vec(i) of *this
"""
cb_rows(self::CBSymmatrix, vec::CBIndexmatrix) = CBMatrix(@ccall libcb.cb_symmatrix_new_rows(self.data::Ptr{Cvoid}, vec.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_swapij!(self::CBSymmatrix, i::Integer, j::Integer)

swaps rows (and columns) i and j
"""
cb_swapij!(self::CBSymmatrix, i::Integer, j::Integer) = (@ccall libcb.cb_symmatrix_swapij(self.data::Ptr{Cvoid}, i::Cint, j::Cint)::Ptr{Cvoid}; return self)

@doc raw"""
    cb_pivot_permute!(self::CBSymmatrix, piv::CBIndexmatrix, inverse::Bool = false)

for i=0 to rowdim row (and column) i of this matrix is swapped with row piv(j); for inverse=true the inverse permutation is generated
"""
cb_pivot_permute!(self::CBSymmatrix, piv::CBIndexmatrix, inverse::Bool = false) = (@ccall libcb.cb_symmatrix_pivot_permute(self.data::Ptr{Cvoid}, piv.data::Ptr{Cvoid}, inverse::Cint)::Ptr{Cvoid}; return self)

@doc raw"""
    cb_principal_submatrix(self::CBSymmatrix, ind::CBIndexmatrix, S::CBSymmatrix)

returns S and in S the principal submatrix indexed by ind (multiple indices are allowed)
"""
cb_principal_submatrix(self::CBSymmatrix, ind::CBIndexmatrix, S::CBSymmatrix) = (@ccall libcb.cb_symmatrix_principal_submatrix(self.data::Ptr{Cvoid}, ind.data::Ptr{Cvoid}, S.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

@doc raw"""
    cb_principal_submatrix(self::CBSymmatrix, ind::CBIndexmatrix)

returns the principal submatrix indexed by ind (multiple indices are allowed)
"""
cb_principal_submatrix(self::CBSymmatrix, ind::CBIndexmatrix) = CBSymmatrix(@ccall libcb.cb_symmatrix_new_principal_submatrix(self.data::Ptr{Cvoid}, ind.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_delete_principal_submatrix!(self::CBSymmatrix, ind::CBIndexmatrix, sorted_increasingly::Bool = false)

returns this afte deleting the principal submatrix indexed by ind (no repetitions!);
"""
cb_delete_principal_submatrix!(self::CBSymmatrix, ind::CBIndexmatrix, sorted_increasingly::Bool = false) = (@ccall libcb.cb_symmatrix_delete_principal_submatrix(self.data::Ptr{Cvoid}, ind.data::Ptr{Cvoid}, sorted_increasingly::Cint)::Ptr{Cvoid}; return self)

@doc raw"""
    cb_enlarge_below!(self::CBSymmatrix, addn::Integer)

increases the order of the matrix by appending storage for further addn rows and columns (marked as not initiliazed if addn>0, no changes if addn<=0)
"""
cb_enlarge_below!(self::CBSymmatrix, addn::Integer) = (@ccall libcb.cb_symmatrix_enlarge_below(self.data::Ptr{Cvoid}, addn::Cint)::Ptr{Cvoid}; return self)

@doc raw"""
    cb_enlarge_below!(self::CBSymmatrix, addn::Integer, d::Real)

increases the order of the matrix by appending storage for further addn rows and columns initialized to d (no changes if addn<=0);
"""
cb_enlarge_below!(self::CBSymmatrix, addn::Integer, d::Real) = (@ccall libcb.cb_symmatrix_enlarge_below2(self.data::Ptr{Cvoid}, addn::Cint, d::Cdouble)::Ptr{Cvoid}; return self)

@doc raw"""
    cb_get_store(self::CBSymmatrix)

returns the current address of the internal value array; use cautiously!
"""
cb_get_store(self::CBSymmatrix) = @ccall libcb.cb_symmatrix_get_store2(self.data::Ptr{Cvoid})::Ptr{Cdouble}

@doc raw"""
    cb_diag(A::CBSymmatrix)

returns a column vector v consisting of the elements v(i)=(*this)(i,i), 0<=i<row dimension
"""
cb_diag(A::CBSymmatrix) = CBMatrix(@ccall libcb.cb_symmatrix_new_diag(A.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_Diag(A::CBMatrix)

returns a symmetric diagonal matrix S of order A.dim() with vec(A) on the diagonal, i.e., S(i,i)=A(i) for all i and S(i,j)=0 for i!=j
"""
cb_Diag(A::CBMatrix) = CBSymmatrix(@ccall libcb.cb_symmatrix_new_diag2(A.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_swap(A::CBSymmatrix, B::CBSymmatrix)

swap the content of the two matrices A and B (involves no copying)
"""
cb_swap(A::CBSymmatrix, B::CBSymmatrix) = @ccall libcb.cb_symmatrix_swap(A.data::Ptr{Cvoid}, B.data::Ptr{Cvoid})::Cvoid

@doc raw"""
    cb_xeya!(self::CBSymmatrix, A::CBSymmatrix, d::Real = 1.)

sets *this=d*A and returns *this
"""
cb_xeya!(self::CBSymmatrix, A::CBSymmatrix, d::Real = 1.) = (@ccall libcb.cb_symmatrix_xeya(self.data::Ptr{Cvoid}, A.data::Ptr{Cvoid}, d::Cdouble)::Ptr{Cvoid}; return self)

@doc raw"""
    cb_xpeya!(self::CBSymmatrix, A::CBSymmatrix, d::Real = 1.)

sets *this+=d*A and returns *this
"""
cb_xpeya!(self::CBSymmatrix, A::CBSymmatrix, d::Real = 1.) = (@ccall libcb.cb_symmatrix_xpeya(self.data::Ptr{Cvoid}, A.data::Ptr{Cvoid}, d::Cdouble)::Ptr{Cvoid}; return self)

@doc raw"""
    cb_rankadd(A::CBMatrix, C::CBSymmatrix, alpha::Real, beta::Real, trans::Integer)

returns C=beta*C+alpha* A*A^T, where A may be transposed. If beta==0. then C is initiliazed to the correct size.
"""
cb_rankadd(A::CBMatrix, C::CBSymmatrix, alpha::Real, beta::Real, trans::Integer) = (@ccall libcb.cb_symmatrix_rankadd(A.data::Ptr{Cvoid}, C.data::Ptr{Cvoid}, alpha::Cdouble, beta::Cdouble, trans::Cint)::Ptr{Cvoid}; return self)

@doc raw"""
    cb_scaledrankadd(A::CBMatrix, D::CBMatrix, C::CBSymmatrix, alpha::Real, beta::Real, trans::Integer)

returns C=beta*C+alpha* A*D*A^T, where D is a vector representing a diagonal matrix and A may be transposed; if beta==0. then C is initialized to the correct size
"""
cb_scaledrankadd(A::CBMatrix, D::CBMatrix, C::CBSymmatrix, alpha::Real, beta::Real, trans::Integer) = (@ccall libcb.cb_symmatrix_scaledrankadd(A.data::Ptr{Cvoid}, D.data::Ptr{Cvoid}, C.data::Ptr{Cvoid}, alpha::Cdouble, beta::Cdouble, trans::Cint)::Ptr{Cvoid}; return self)

@doc raw"""
    cb_rank2add(A::CBMatrix, B::CBMatrix, C::CBSymmatrix, alpha::Real, beta::Real, trans::Integer)

returns C=beta*C+alpha*(A*B^T+B*A^T)/2 [or for transposed (A^T*B+B^T*A)/2]. If beta==0. then C is initiliazed to the correct size.
"""
cb_rank2add(A::CBMatrix, B::CBMatrix, C::CBSymmatrix, alpha::Real, beta::Real, trans::Integer) = (@ccall libcb.cb_symmatrix_rank2add(A.data::Ptr{Cvoid}, B.data::Ptr{Cvoid}, C.data::Ptr{Cvoid}, alpha::Cdouble, beta::Cdouble, trans::Cint)::Ptr{Cvoid}; return self)

@doc raw"""
    cb_xbpeya(x::CBSymmatrix, y::CBSymmatrix, alpha::Real, beta::Real)

returns x= alpha*y+beta*x; if beta==0. then x is initialized to the correct size
"""
cb_xbpeya(x::CBSymmatrix, y::CBSymmatrix, alpha::Real, beta::Real) = (@ccall libcb.cb_symmatrix_xbpeya(x.data::Ptr{Cvoid}, y.data::Ptr{Cvoid}, alpha::Cdouble, beta::Cdouble)::Ptr{Cvoid}; return self)

@doc raw"""
    cb_xeyapzb(x::CBSymmatrix, y::CBSymmatrix, z::CBSymmatrix, alpha::Real, beta::Real)

returns x= alpha*y+beta*z; x is initialized to the correct size
"""
cb_xeyapzb(x::CBSymmatrix, y::CBSymmatrix, z::CBSymmatrix, alpha::Real, beta::Real) = (@ccall libcb.cb_symmatrix_xeyapzb(x.data::Ptr{Cvoid}, y.data::Ptr{Cvoid}, z.data::Ptr{Cvoid}, alpha::Cdouble, beta::Cdouble)::Ptr{Cvoid}; return self)

Base.copy!(self::CBSymmatrix, A::CBSymmatrix) = (@ccall libcb.cb_symmatrix_assign(self.data::Ptr{Cvoid}, A.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

MA.operate!(::typeof(+), self::CBSymmatrix, A::CBSymmatrix) = (@ccall libcb.cb_symmatrix_plus(self.data::Ptr{Cvoid}, A.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

MA.promote_operation(::typeof(+), ::Type{<:CBSymmatrix}, ::Type{<:CBSymmatrix}) = CBSymmatrix

MA.operate!(::typeof(-), self::CBSymmatrix, A::CBSymmatrix) = (@ccall libcb.cb_symmatrix_minus(self.data::Ptr{Cvoid}, A.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

MA.promote_operation(::typeof(-), ::Type{<:CBSymmatrix}, ::Type{<:CBSymmatrix}) = CBSymmatrix

@doc raw"""
    MA.operate!(::typeof(%), self::CBSymmatrix, A::CBSymmatrix)

ATTENTION: this is redefined as the Hadamard product, (*this)(i,j)=(*this)(i,j)*A(i,j) for all i<=j
"""
MA.operate!(::typeof(%), self::CBSymmatrix, A::CBSymmatrix) = (@ccall libcb.cb_symmatrix_rem(self.data::Ptr{Cvoid}, A.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

MA.promote_operation(::typeof(%), ::Type{<:CBSymmatrix}, ::Type{<:CBSymmatrix}) = CBSymmatrix

Base.:-(self::CBSymmatrix) = CBSymmatrix(@ccall libcb.cb_symmatrix_new_minus(self.data::Ptr{Cvoid})::Ptr{Cvoid})

MA.operate!(::typeof(*), self::CBSymmatrix, d::Real) = (@ccall libcb.cb_symmatrix_times(self.data::Ptr{Cvoid}, d::Cdouble)::Ptr{Cvoid}; return self)

MA.promote_operation(::typeof(*), ::Type{<:CBSymmatrix}, ::Type{<:Real}) = CBSymmatrix

@doc raw"""
    MA.operate!(::typeof(/), self::CBSymmatrix, d::Real)

ATTENTION: d is NOT checked for 0
"""
MA.operate!(::typeof(/), self::CBSymmatrix, d::Real) = (@ccall libcb.cb_symmatrix_divide(self.data::Ptr{Cvoid}, d::Cdouble)::Ptr{Cvoid}; return self)

MA.promote_operation(::typeof(/), ::Type{<:CBSymmatrix}, ::Type{<:Real}) = CBSymmatrix

@doc raw"""
    MA.operate!(::typeof(+), self::CBSymmatrix, d::Real)

sets (*this)(i,j)+=d for all i<=j
"""
MA.operate!(::typeof(+), self::CBSymmatrix, d::Real) = (@ccall libcb.cb_symmatrix_plus2(self.data::Ptr{Cvoid}, d::Cdouble)::Ptr{Cvoid}; return self)

MA.promote_operation(::typeof(+), ::Type{<:CBSymmatrix}, ::Type{<:Real}) = CBSymmatrix

@doc raw"""
    MA.operate!(::typeof(-), self::CBSymmatrix, d::Real)

sets (*this)(i,j)-=d for all i<=j
"""
MA.operate!(::typeof(-), self::CBSymmatrix, d::Real) = (@ccall libcb.cb_symmatrix_minus2(self.data::Ptr{Cvoid}, d::Cdouble)::Ptr{Cvoid}; return self)

MA.promote_operation(::typeof(-), ::Type{<:CBSymmatrix}, ::Type{<:Real}) = CBSymmatrix

@doc raw"""
    cb_transpose!(self::CBSymmatrix)

transposes itself (at almost no cost)
"""
cb_transpose!(self::CBSymmatrix) = (@ccall libcb.cb_symmatrix_transpose(self.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

Base.:*(A::CBSymmatrix, B::CBSymmatrix) = CBMatrix(@ccall libcb.cb_symmatrix_new_times(A.data::Ptr{Cvoid}, B.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    Base.:%(A::CBSymmatrix, B::CBSymmatrix)

ATTENTION: this is redefined as the Hadamard product and sets (i,j)=A(i,j)*B(i,j) for all i<=j
"""
Base.:%(A::CBSymmatrix, B::CBSymmatrix) = CBSymmatrix(@ccall libcb.cb_symmatrix_new_rem(A.data::Ptr{Cvoid}, B.data::Ptr{Cvoid})::Ptr{Cvoid})

Base.:+(A::CBSymmatrix, B::CBSymmatrix) = CBSymmatrix(@ccall libcb.cb_symmatrix_new_plus(A.data::Ptr{Cvoid}, B.data::Ptr{Cvoid})::Ptr{Cvoid})

Base.:-(A::CBSymmatrix, B::CBSymmatrix) = CBSymmatrix(@ccall libcb.cb_symmatrix_new_minus2(A.data::Ptr{Cvoid}, B.data::Ptr{Cvoid})::Ptr{Cvoid})

Base.:*(A::CBSymmatrix, B::CBMatrix) = CBMatrix(@ccall libcb.cb_symmatrix_new_times2(A.data::Ptr{Cvoid}, B.data::Ptr{Cvoid})::Ptr{Cvoid})

Base.:*(A::CBMatrix, B::CBSymmatrix) = CBMatrix(@ccall libcb.cb_symmatrix_new_times3(A.data::Ptr{Cvoid}, B.data::Ptr{Cvoid})::Ptr{Cvoid})

Base.:+(A::CBSymmatrix, B::CBMatrix) = CBMatrix(@ccall libcb.cb_symmatrix_new_plus2(A.data::Ptr{Cvoid}, B.data::Ptr{Cvoid})::Ptr{Cvoid})

Base.:+(A::CBMatrix, B::CBSymmatrix) = CBMatrix(@ccall libcb.cb_symmatrix_new_plus3(A.data::Ptr{Cvoid}, B.data::Ptr{Cvoid})::Ptr{Cvoid})

Base.:-(A::CBSymmatrix, B::CBMatrix) = CBMatrix(@ccall libcb.cb_symmatrix_new_minus3(A.data::Ptr{Cvoid}, B.data::Ptr{Cvoid})::Ptr{Cvoid})

Base.:-(A::CBMatrix, B::CBSymmatrix) = CBMatrix(@ccall libcb.cb_symmatrix_new_minus4(A.data::Ptr{Cvoid}, B.data::Ptr{Cvoid})::Ptr{Cvoid})

Base.:*(A::CBSymmatrix, d::Real) = CBSymmatrix(@ccall libcb.cb_symmatrix_new_times4(A.data::Ptr{Cvoid}, d::Cdouble)::Ptr{Cvoid})

Base.:*(d::Real, A::CBSymmatrix) = CBSymmatrix(@ccall libcb.cb_symmatrix_new_times5(d::Cdouble, A.data::Ptr{Cvoid})::Ptr{Cvoid})

Base.:/(A::CBSymmatrix, d::Real) = CBSymmatrix(@ccall libcb.cb_symmatrix_new_divide(A.data::Ptr{Cvoid}, d::Cdouble)::Ptr{Cvoid})

@doc raw"""
    Base.:+(A::CBSymmatrix, d::Real)

returns (i,j)=A(i,j)+d for all i<=j
"""
Base.:+(A::CBSymmatrix, d::Real) = CBSymmatrix(@ccall libcb.cb_symmatrix_new_plus4(A.data::Ptr{Cvoid}, d::Cdouble)::Ptr{Cvoid})

@doc raw"""
    Base.:+(d::Real, A::CBSymmatrix)

returns (i,j)=A(i,j)+d for all i<=j
"""
Base.:+(d::Real, A::CBSymmatrix) = CBSymmatrix(@ccall libcb.cb_symmatrix_new_plus5(d::Cdouble, A.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    Base.:-(A::CBSymmatrix, d::Real)

returns (i,j)=A(i,j)-d for all i<=j
"""
Base.:-(A::CBSymmatrix, d::Real) = CBSymmatrix(@ccall libcb.cb_symmatrix_new_minus5(A.data::Ptr{Cvoid}, d::Cdouble)::Ptr{Cvoid})

@doc raw"""
    Base.:-(d::Real, A::CBSymmatrix)

returns (i,j)=d-A(i,j) for all i<=j
"""
Base.:-(d::Real, A::CBSymmatrix) = CBSymmatrix(@ccall libcb.cb_symmatrix_new_minus6(d::Cdouble, A.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_transpose(A::CBSymmatrix)

(drop it or use a constructor instead)
"""
cb_transpose(A::CBSymmatrix) = CBSymmatrix(@ccall libcb.cb_symmatrix_new_transpose(A.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_xeya!(self::CBSymmatrix, A::CBMatrix, d::Real = 1.)

sets *this=d*(A+transpose(A))/2. and returns *this
"""
cb_xeya!(self::CBSymmatrix, A::CBMatrix, d::Real = 1.) = (@ccall libcb.cb_symmatrix_xeya2(self.data::Ptr{Cvoid}, A.data::Ptr{Cvoid}, d::Cdouble)::Ptr{Cvoid}; return self)

@doc raw"""
    cb_xpeya!(self::CBSymmatrix, A::CBMatrix, d::Real = 1.)

sets *this+=d*(A+transpose(A))/2. and returns *this
"""
cb_xpeya!(self::CBSymmatrix, A::CBMatrix, d::Real = 1.) = (@ccall libcb.cb_symmatrix_xpeya2(self.data::Ptr{Cvoid}, A.data::Ptr{Cvoid}, d::Cdouble)::Ptr{Cvoid}; return self)

@doc raw"""
    cb_xeya!(self::CBSymmatrix, A::CBIndexmatrix, d::Real = 1.)

sets *this=d*(A+transpose(A))/2. and returns *this
"""
cb_xeya!(self::CBSymmatrix, A::CBIndexmatrix, d::Real = 1.) = (@ccall libcb.cb_symmatrix_xeya3(self.data::Ptr{Cvoid}, A.data::Ptr{Cvoid}, d::Cdouble)::Ptr{Cvoid}; return self)

@doc raw"""
    cb_xpeya!(self::CBSymmatrix, A::CBIndexmatrix, d::Real = 1.)

sets *this+=d*(A+transpose(A))/2. and returns *this
"""
cb_xpeya!(self::CBSymmatrix, A::CBIndexmatrix, d::Real = 1.) = (@ccall libcb.cb_symmatrix_xpeya3(self.data::Ptr{Cvoid}, A.data::Ptr{Cvoid}, d::Cdouble)::Ptr{Cvoid}; return self)

@doc raw"""
    cb_xeya!(self::CBSymmatrix, A::CBSparsesym, d::Real = 1.)

sets *this=d*A and returns *this
"""
cb_xeya!(self::CBSymmatrix, A::CBSparsesym, d::Real = 1.) = (@ccall libcb.cb_symmatrix_xeya4(self.data::Ptr{Cvoid}, A.data::Ptr{Cvoid}, d::Cdouble)::Ptr{Cvoid}; return self)

@doc raw"""
    cb_xpeya!(self::CBSymmatrix, A::CBSparsesym, d::Real = 1.)

sets *this+=d*A and returns *this
"""
cb_xpeya!(self::CBSymmatrix, A::CBSparsesym, d::Real = 1.) = (@ccall libcb.cb_symmatrix_xpeya4(self.data::Ptr{Cvoid}, A.data::Ptr{Cvoid}, d::Cdouble)::Ptr{Cvoid}; return self)

@doc raw"""
    cb_xetriu_yza!(self::CBSymmatrix, A::CBMatrix, B::CBMatrix, d::Real = 1.)

sets *this(i,j), i<=j to the upper triangle of the matrix product d*transpose(A)*B
"""
cb_xetriu_yza!(self::CBSymmatrix, A::CBMatrix, B::CBMatrix, d::Real = 1.) = (@ccall libcb.cb_symmatrix_xetriu_yza(self.data::Ptr{Cvoid}, A.data::Ptr{Cvoid}, B.data::Ptr{Cvoid}, d::Cdouble)::Ptr{Cvoid}; return self)

@doc raw"""
    cb_xpetriu_yza!(self::CBSymmatrix, A::CBMatrix, B::CBMatrix, d::Real = 1.)

adds to *this(i,j), i<=j the upper triangle of the matrix product d*transpose(A)*B
"""
cb_xpetriu_yza!(self::CBSymmatrix, A::CBMatrix, B::CBMatrix, d::Real = 1.) = (@ccall libcb.cb_symmatrix_xpetriu_yza(self.data::Ptr{Cvoid}, A.data::Ptr{Cvoid}, B.data::Ptr{Cvoid}, d::Cdouble)::Ptr{Cvoid}; return self)

@doc raw"""
    cb_xetriu_yza!(self::CBSymmatrix, A::CBSparsemat, B::CBMatrix, d::Real = 1.)

sets *this(i,j), i<=j to the upper triangle of the matrix product d*transpose(A)*B
"""
cb_xetriu_yza!(self::CBSymmatrix, A::CBSparsemat, B::CBMatrix, d::Real = 1.) = (@ccall libcb.cb_symmatrix_xetriu_yza2(self.data::Ptr{Cvoid}, A.data::Ptr{Cvoid}, B.data::Ptr{Cvoid}, d::Cdouble)::Ptr{Cvoid}; return self)

@doc raw"""
    cb_xpetriu_yza!(self::CBSymmatrix, A::CBSparsemat, B::CBMatrix, d::Real = 1.)

adds to *this(i,j), i<=j the upper triangle of the matrix product d*transpose(A)*B
"""
cb_xpetriu_yza!(self::CBSymmatrix, A::CBSparsemat, B::CBMatrix, d::Real = 1.) = (@ccall libcb.cb_symmatrix_xpetriu_yza2(self.data::Ptr{Cvoid}, A.data::Ptr{Cvoid}, B.data::Ptr{Cvoid}, d::Cdouble)::Ptr{Cvoid}; return self)

Base.copy!(self::CBSymmatrix, A::CBSparsesym) = (@ccall libcb.cb_symmatrix_assign2(self.data::Ptr{Cvoid}, A.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

MA.operate!(::typeof(+), self::CBSymmatrix, A::CBSparsesym) = (@ccall libcb.cb_symmatrix_plus3(self.data::Ptr{Cvoid}, A.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

MA.promote_operation(::typeof(+), ::Type{<:CBSymmatrix}, ::Type{<:CBSparsesym}) = CBSymmatrix

MA.operate!(::typeof(-), self::CBSymmatrix, A::CBSparsesym) = (@ccall libcb.cb_symmatrix_minus3(self.data::Ptr{Cvoid}, A.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

MA.promote_operation(::typeof(-), ::Type{<:CBSymmatrix}, ::Type{<:CBSparsesym}) = CBSymmatrix

@doc raw"""
    cb_abs(A::CBSymmatrix)

returns a Symmatrix with elements abs(A(i,j))
"""
cb_abs(A::CBSymmatrix) = CBSymmatrix(@ccall libcb.cb_symmatrix_new_abs(A.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_shift_diag!(self::CBSymmatrix, s::Real)

shifts the diagonal by s, i.e., (*this)(i,i)+=s for all i
"""
cb_shift_diag!(self::CBSymmatrix, s::Real) = (@ccall libcb.cb_symmatrix_shift_diag(self.data::Ptr{Cvoid}, s::Cdouble)::Ptr{Cvoid}; return self)

@doc raw"""
    cb_LDLfactor!(self::CBSymmatrix, tol::Real = 1e-10)

computes LDLfactorization (implemented only for positive definite matrices so far, no pivoting), (*this) is overwritten by the factorization; returns 1 if diagonal elements go below tol
"""
cb_LDLfactor!(self::CBSymmatrix, tol::Real = 1e-10) = @ccall libcb.cb_symmatrix_ldlfactor(self.data::Ptr{Cvoid}, tol::Cdouble)::Cint

@doc raw"""
    cb_LDLsolve(self::CBSymmatrix, x::CBMatrix)

computes, after LDLfactor was executed succesfully, the solution to (*old_this)x=rhs; rhs is overwritten by the solution; always returns 0; NOTE: there is NO check against division by zero
"""
cb_LDLsolve(self::CBSymmatrix, x::CBMatrix) = @ccall libcb.cb_symmatrix_ldlsolve(self.data::Ptr{Cvoid}, x.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_LDLinverse(self::CBSymmatrix, S::CBSymmatrix)

computes, after LDLfactor was executed succesfully, the inverse to (*old_this) and stores it in S (numerically not too wise); always returns 0; NOTE: there is NO check against division by zero
"""
cb_LDLinverse(self::CBSymmatrix, S::CBSymmatrix) = @ccall libcb.cb_symmatrix_ldlinverse(self.data::Ptr{Cvoid}, S.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_Chol_factor!(self::CBSymmatrix, tol::Real = 1e-10)

computes the Cholesky factorization, for positive definite matrices only, (*this) is overwritten by the factorization; there is no pivoting; returns 1 if diagonal elements go below tol
"""
cb_Chol_factor!(self::CBSymmatrix, tol::Real = 1e-10) = @ccall libcb.cb_symmatrix_chol_factor(self.data::Ptr{Cvoid}, tol::Cdouble)::Cint

@doc raw"""
    cb_Chol_solve(self::CBSymmatrix, x::CBMatrix)

computes, after Chol_factor was executed succesfully, the solution to (*old_this)x=rhs; rhs is overwritten by the solution; always returns 0; NOTE: there is NO check against division by zero
"""
cb_Chol_solve(self::CBSymmatrix, x::CBMatrix) = @ccall libcb.cb_symmatrix_chol_solve(self.data::Ptr{Cvoid}, x.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_Chol_inverse(self::CBSymmatrix, S::CBSymmatrix)

computes, after Chol_factor was executed succesfully, the inverse to (*old_this) and stores it in S (numerically not too wise); always returns 0; NOTE: there is NO check against division by zero
"""
cb_Chol_inverse(self::CBSymmatrix, S::CBSymmatrix) = @ccall libcb.cb_symmatrix_chol_inverse(self.data::Ptr{Cvoid}, S.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_Chol_Lsolve(self::CBSymmatrix, rhs::CBMatrix)

computes, after Chol_factor into LL^T was executed succesfully, the solution to Lx=rhs; rhs is overwritten by the solution; always returns 0; NOTE: there is NO check against division by zero
"""
cb_Chol_Lsolve(self::CBSymmatrix, rhs::CBMatrix) = @ccall libcb.cb_symmatrix_chol_lsolve(self.data::Ptr{Cvoid}, rhs.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_Chol_Ltsolve(self::CBSymmatrix, rhs::CBMatrix)

computes, after Chol_factor into LL^T was executed succesfully, the solution to L^Tx=rhs; rhs is overwritten by the solution; always returns 0; NOTE: there is NO check against division by zero
"""
cb_Chol_Ltsolve(self::CBSymmatrix, rhs::CBMatrix) = @ccall libcb.cb_symmatrix_chol_ltsolve(self.data::Ptr{Cvoid}, rhs.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_Chol_scaleLi(self::CBSymmatrix, S::CBSymmatrix)

computes, after Chol_factor into LL^T was executed succesfully, L^{-1}SL^{-T} overwriting S
"""
cb_Chol_scaleLi(self::CBSymmatrix, S::CBSymmatrix) = @ccall libcb.cb_symmatrix_chol_scaleli(self.data::Ptr{Cvoid}, S.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_Chol_scaleLt(self::CBSymmatrix, S::CBSymmatrix)

computes, after Chol_factor into LL^T was executed succesfully, L^TSL overwriting S
"""
cb_Chol_scaleLt(self::CBSymmatrix, S::CBSymmatrix) = @ccall libcb.cb_symmatrix_chol_scalelt(self.data::Ptr{Cvoid}, S.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_Chol_Lmult(self::CBSymmatrix, rhs::CBMatrix)

computes, after Chol_factor into LL^T was executed succesfully,  L*rhs, overwriting rhs by the result; always returns 0;
"""
cb_Chol_Lmult(self::CBSymmatrix, rhs::CBMatrix) = @ccall libcb.cb_symmatrix_chol_lmult(self.data::Ptr{Cvoid}, rhs.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_Chol_Ltmult(self::CBSymmatrix, rhs::CBMatrix)

computes, after Chol_factor into LL^T was executed succesfully,  L^Trhs, overwriting rhs by the result; always returns 0;
"""
cb_Chol_Ltmult(self::CBSymmatrix, rhs::CBMatrix) = @ccall libcb.cb_symmatrix_chol_ltmult(self.data::Ptr{Cvoid}, rhs.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_Chol_factor!(self::CBSymmatrix, piv::CBIndexmatrix, tol::Real = 1e-10)

computes the Cholesky factorization with pivoting, for positive semidefinite matrices only, (*this) is overwritten by the factorization; on termination piv.dim() is the number of positive pivots>=tol; returns 1 if negative diagonal element is encountered during computations, 0 otherwise.
"""
cb_Chol_factor!(self::CBSymmatrix, piv::CBIndexmatrix, tol::Real = 1e-10) = @ccall libcb.cb_symmatrix_chol_factor2(self.data::Ptr{Cvoid}, piv.data::Ptr{Cvoid}, tol::Cdouble)::Cint

@doc raw"""
    cb_Chol_solve(self::CBSymmatrix, x::CBMatrix, piv::CBIndexmatrix)

computes, after Chol_factor(Indexmatrix&,Real) with pivoting was executed succesfully, the solution to (*old_this)*x=rhs(piv); rhs is overwritten by the solution arranged in original unpermuted order; always returns 0; NOTE: there is NO check against division by zero
"""
cb_Chol_solve(self::CBSymmatrix, x::CBMatrix, piv::CBIndexmatrix) = @ccall libcb.cb_symmatrix_chol_solve2(self.data::Ptr{Cvoid}, x.data::Ptr{Cvoid}, piv.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_Chol_inverse(self::CBSymmatrix, S::CBSymmatrix, piv::CBIndexmatrix)

computes, after Chol_factor(Indexmatrix&,Real) with pivoting was executed succesfully, the inverse to (*old_this) and stores it in S (the pivoting permutation is undone in S); NOTE: there is NO check against division by zero
"""
cb_Chol_inverse(self::CBSymmatrix, S::CBSymmatrix, piv::CBIndexmatrix) = @ccall libcb.cb_symmatrix_chol_inverse2(self.data::Ptr{Cvoid}, S.data::Ptr{Cvoid}, piv.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_Aasen_factor!(self::CBSymmatrix, piv::CBIndexmatrix)

computes Aasen factorization LTL^T with pivoting, where L is unit lower triangular with first colum e_1 and T is tridiagonal; (*this) is overwritten by the factorization, with column i of L being stored in column i-1 of (*this); always returns 0;
"""
cb_Aasen_factor!(self::CBSymmatrix, piv::CBIndexmatrix) = @ccall libcb.cb_symmatrix_aasen_factor(self.data::Ptr{Cvoid}, piv.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_Aasen_Lsolve(self::CBSymmatrix, x::CBMatrix)

computes, after Aasen_factor into LTL^T was executed, the solution to Lx=rhs; rhs is overwritten by the solution; always returns 0;
"""
cb_Aasen_Lsolve(self::CBSymmatrix, x::CBMatrix) = @ccall libcb.cb_symmatrix_aasen_lsolve(self.data::Ptr{Cvoid}, x.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_Aasen_Ltsolve(self::CBSymmatrix, x::CBMatrix)

computes, after Aasen_factor into LTL^T was executed, the solution to L^Tx=rhs; rhs is overwritten by the solution; always returns 0;
"""
cb_Aasen_Ltsolve(self::CBSymmatrix, x::CBMatrix) = @ccall libcb.cb_symmatrix_aasen_ltsolve(self.data::Ptr{Cvoid}, x.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_Aasen_tridiagsolve(self::CBSymmatrix, x::CBMatrix)

computes, after Aasen_factor into LTL^T was executed, the solution to Tx=rhs; rhs is overwritten by the solution;if the solution fails due to division by zero (=system not solvable) the return value is -(rowindex+1) where this occured in the backsolve
"""
cb_Aasen_tridiagsolve(self::CBSymmatrix, x::CBMatrix) = @ccall libcb.cb_symmatrix_aasen_tridiagsolve(self.data::Ptr{Cvoid}, x.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_Aasen_solve(self::CBSymmatrix, x::CBMatrix, piv::CBIndexmatrix)

computes, after Aasen_factor into LTL^T was executed, the solution to (*old_this)x=rhs; rhs is overwritten by the solution; if the solution fails due to division by zero (=system not solvable) the return value is -(rowindex+1) where this occured in the backsolve
"""
cb_Aasen_solve(self::CBSymmatrix, x::CBMatrix, piv::CBIndexmatrix) = @ccall libcb.cb_symmatrix_aasen_solve(self.data::Ptr{Cvoid}, x.data::Ptr{Cvoid}, piv.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_eig(self::CBSymmatrix, P::CBMatrix, d::CBMatrix, sort_non_decreasingly::Bool = true)

computes an eigenvalue decomposition P*Diag(d)*tranpose(P)=(*this) by symmetric QR; returns 0 on success,
"""
cb_eig(self::CBSymmatrix, P::CBMatrix, d::CBMatrix, sort_non_decreasingly::Bool = true) = @ccall libcb.cb_symmatrix_eig(self.data::Ptr{Cvoid}, P.data::Ptr{Cvoid}, d.data::Ptr{Cvoid}, sort_non_decreasingly::Cint)::Cint

@doc raw"""
    cb_trace(A::CBSymmatrix)

returns the sum of the diagonal elements A(i,i) over all i
"""
cb_trace(A::CBSymmatrix) = @ccall libcb.cb_symmatrix_trace(A.data::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_ip(A::CBSymmatrix, B::CBSymmatrix)

returns the usual inner product of A and B, i.e., the sum of A(i,j)*B(i,j) over all i,j
"""
cb_ip(A::CBSymmatrix, B::CBSymmatrix) = @ccall libcb.cb_symmatrix_ip(A.data::Ptr{Cvoid}, B.data::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_ip(A::CBMatrix, B::CBSymmatrix)

returns the usual inner product of A and B, i.e., the sum of A(i,j)*B(i,j) over all i,j
"""
cb_ip(A::CBMatrix, B::CBSymmatrix) = @ccall libcb.cb_symmatrix_ip2(A.data::Ptr{Cvoid}, B.data::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_ip(A::CBSymmatrix, B::CBMatrix)

returns the usual inner product of A and B, i.e., the sum of A(i,j)*B(i,j) over all i,j
"""
cb_ip(A::CBSymmatrix, B::CBMatrix) = @ccall libcb.cb_symmatrix_ip3(A.data::Ptr{Cvoid}, B.data::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_norm2(A::CBSymmatrix)

returns the Frobenius norm of A, i.e., the square root of the sum of A(i,j)*A(i,j) over all i,j
"""
cb_norm2(A::CBSymmatrix) = @ccall libcb.cb_symmatrix_norm2(A.data::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_sumrows(A::CBSymmatrix)

returns a row vector holding the sum over all rows, i.e., (1 1 ... 1)*A
"""
cb_sumrows(A::CBSymmatrix) = CBMatrix(@ccall libcb.cb_symmatrix_new_sumrows(A.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_sumcols(A::CBSymmatrix)

returns a column vector holding the sum over all columns, i.e., A*(1 1 ... 1)^T
"""
cb_sumcols(A::CBSymmatrix) = CBMatrix(@ccall libcb.cb_symmatrix_new_sumcols(A.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_sum(A::CBSymmatrix)

returns the sum over all elements of A, i.e., (1 1 ... 1)*A*(1 1 ... 1)^T
"""
cb_sum(A::CBSymmatrix) = @ccall libcb.cb_symmatrix_sum(A.data::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_svec(A::CBSymmatrix, v::CBMatrix, a::Real, add::Bool, startindex_vec::Integer, startindex_A::Integer, blockdim::Integer)

the symmetric vec operator stacks the lower triangle of A to a n*(n+1)/2 vector with the same norm2 as A; here it sets svec(A)=[a11,sqrt(2)a12,...,sqrt(2)a1n,a22,...,sqrt(2)a(n-1,n),ann]', multiplies it by a and sets or adds (if add==true) it to v starting from startindex_vec possibly restricted to the subblock of order blockdim (whenever >=0, else blockdim is set to A.rowdim()-startindex_A) starting from startindex_A (must be >=0); if add==false and startindex_vec<0 then vec is also reinitialzed to the appropriate size
"""
cb_svec(A::CBSymmatrix, v::CBMatrix, a::Real, add::Bool, startindex_vec::Integer, startindex_A::Integer, blockdim::Integer) = @ccall libcb.cb_symmatrix_svec(A.data::Ptr{Cvoid}, v.data::Ptr{Cvoid}, a::Cdouble, add::Cint, startindex_vec::Cint, startindex_A::Cint, blockdim::Cint)::Cvoid

@doc raw"""
    cb_svec(A::CBSymmatrix)

the symmetric vec operator, stacks the lower triangle of A to a n*(n+1)/2 vector with the same norm2 as A; i.e., it returns svec(A)=[a11,sqrt(2)a12,...,sqrt(2)a1n,a22,...,sqrt(2)a(n-1,n),ann]'
"""
cb_svec(A::CBSymmatrix) = CBMatrix(@ccall libcb.cb_symmatrix_new_svec(A.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_sveci(v::CBMatrix, A::CBSymmatrix, a::Real, add::Bool, startindex_vec::Integer, startindex_A::Integer, blockdim::Integer)

the inverse operator to svec, extracts from v at startindex_vec (>=0) the symmetric matrix of blockdim adding its mutliple by a into A starting at startindex_A; if add==false and startindex_A<0 A is initialized to the size of blockdim; if the latter is also negative then v.dim()-startindex_vec must match an exact order and  matrix A is initialized to this size. In all other cases the size of the symmetric matrix determines the missing parameters and vec.dim-startindex_vec
"""
cb_sveci(v::CBMatrix, A::CBSymmatrix, a::Real, add::Bool, startindex_vec::Integer, startindex_A::Integer, blockdim::Integer) = @ccall libcb.cb_symmatrix_sveci(v.data::Ptr{Cvoid}, A.data::Ptr{Cvoid}, a::Cdouble, add::Cint, startindex_vec::Cint, startindex_A::Cint, blockdim::Cint)::Cvoid

@doc raw"""
    cb_init_svec!(self::CBSymmatrix, nr::Integer, dp::Union{<:AbstractVector{Cdouble},Nothing}, d::Real = 1.)

initialize from an svec stored in a real array (or matrix)
"""
cb_init_svec!(self::CBSymmatrix, nr::Integer, dp::Union{<:AbstractVector{Cdouble},Nothing}, d::Real = 1.) = GC.@preserve dp begin
    @ccall libcb.cb_symmatrix_init_svec(self.data::Ptr{Cvoid}, nr::Cint, dp::Ptr{Cdouble}, stride(dp, 1)::Cint, d::Cdouble)::Cvoid
end

@doc raw"""
    cb_store_svec(self::CBSymmatrix, dp::Union{<:AbstractVector{Cdouble},Nothing}, d::Real = 1.)

store in the form of an svec in the real array (or matrix)
"""
cb_store_svec(self::CBSymmatrix, dp::Union{<:AbstractVector{Cdouble},Nothing}, d::Real = 1.) = GC.@preserve dp begin
    @ccall libcb.cb_symmatrix_store_svec(self.data::Ptr{Cvoid}, dp::Ptr{Cdouble}, stride(dp, 1)::Cint, d::Cdouble)::Cvoid
end

@doc raw"""
    cb_skron(A::CBSymmatrix, B::CBSymmatrix, alpha::Real, add::Bool, startindex_S::Integer)

the symmetric Kronecker product, defined via (A skron B)svec(C)=(BCA'+ACB')/2; sets or adds (if add==true) the symmetric matrix a*(A skron B) into S starting at startindex_S; if add==false and startindex_S<0, S is initialzed to the correct size
"""
cb_skron(A::CBSymmatrix, B::CBSymmatrix, alpha::Real, add::Bool, startindex_S::Integer) = CBSymmatrix(@ccall libcb.cb_symmatrix_new_skron(A.data::Ptr{Cvoid}, B.data::Ptr{Cvoid}, alpha::Cdouble, add::Cint, startindex_S::Cint)::Ptr{Cvoid})

@doc raw"""
    cb_skron(A::CBSymmatrix, B::CBSymmatrix, S::CBSymmatrix, a::Real, add::Bool, startindex_S::Integer)

def symmetric Kronecker product (A skron B)svec(C)=(BCA'+ACB')/2; sets S=alpha*(A skron B) or S*=... (if add==true) possibly shifted to the block starting at startindex_S;  if add==false and startindex_S<0, S is initialzed to the correct size
"""
cb_skron(A::CBSymmatrix, B::CBSymmatrix, S::CBSymmatrix, a::Real, add::Bool, startindex_S::Integer) = @ccall libcb.cb_symmatrix_skron(A.data::Ptr{Cvoid}, B.data::Ptr{Cvoid}, S.data::Ptr{Cvoid}, a::Cdouble, add::Cint, startindex_S::Cint)::Cvoid

@doc raw"""
    cb_symscale(A::CBSymmatrix, B::CBMatrix, S::CBSymmatrix, a::Real, b::Real, btrans::Integer)

sets S=beta*S+alpha*B'*A*B for symmatrix A and matrix B
"""
cb_symscale(A::CBSymmatrix, B::CBMatrix, S::CBSymmatrix, a::Real, b::Real, btrans::Integer) = @ccall libcb.cb_symmatrix_symscale(A.data::Ptr{Cvoid}, B.data::Ptr{Cvoid}, S.data::Ptr{Cvoid}, a::Cdouble, b::Cdouble, btrans::Cint)::Cvoid

@doc raw"""
    cb_minrows(A::CBSymmatrix)

returns a row vector holding in each column the minimum over all rows in this column
"""
cb_minrows(A::CBSymmatrix) = CBMatrix(@ccall libcb.cb_symmatrix_new_minrows(A.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_mincols(A::CBSymmatrix)

returns a column vector holding in each row the minimum over all columns in this row
"""
cb_mincols(A::CBSymmatrix) = CBMatrix(@ccall libcb.cb_symmatrix_new_mincols(A.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_min(A::CBSymmatrix)

returns the minimum value over all elements of the matrix
"""
cb_min(A::CBSymmatrix) = @ccall libcb.cb_symmatrix_min(A.data::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_maxrows(A::CBSymmatrix)

returns a row vector holding in each column the maximum over all rows in this column
"""
cb_maxrows(A::CBSymmatrix) = CBMatrix(@ccall libcb.cb_symmatrix_new_maxrows(A.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_maxcols(A::CBSymmatrix)

returns a column vector holding in each row the maximum over all columns in this row
"""
cb_maxcols(A::CBSymmatrix) = CBMatrix(@ccall libcb.cb_symmatrix_new_maxcols(A.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_max(A::CBSymmatrix)

returns the maximum value over all elements of the matrix
"""
cb_max(A::CBSymmatrix) = @ccall libcb.cb_symmatrix_max(A.data::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_display(self::CBSymmatrix, precision::Integer = 0, width::Integer = 0, screenwidth::Integer = 0)

* @brief displays a matrix in a pretty way for bounded screen widths; for variables of value zero default values are used.
      
"""
cb_display(self::CBSymmatrix, precision::Integer = 0, width::Integer = 0, screenwidth::Integer = 0) = @ccall libcb.cb_symmatrix_display(self.data::Ptr{Cvoid}, precision::Cint, width::Cint, screenwidth::Cint)::Cvoid

@doc raw"""
    cb_mfile_output(self::CBSymmatrix, precision::Integer = 16, width::Integer = 0)

* @brief outputs a matrix A in the format "[ A(0,1) ... A(0,nc-1)\n ... A(nr-1,nc-1)];\n" so that it can be read e.g. by octave as an m-file
     
"""
cb_mfile_output(self::CBSymmatrix, precision::Integer = 16, width::Integer = 0) = @ccall libcb.cb_symmatrix_mfile_output(self.data::Ptr{Cvoid}, precision::Cint, width::Cint)::Cvoid

