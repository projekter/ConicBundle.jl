@doc raw"""
    CBSparsesym()

empty matrix
"""
CBSparsesym() = CBSparsesym(@ccall libcb.cb_sparsesym_new()::Ptr{Cvoid})

@doc raw"""
    CBSparsesym(A::CBSparsesym, d::Real = 1.)

copy constructor, *this=d*A
"""
CBSparsesym(A::CBSparsesym, d::Real = 1.) = CBSparsesym(@ccall libcb.cb_sparsesym_new2(A.data::Ptr{Cvoid}, d::Cdouble)::Ptr{Cvoid})

@doc raw"""
    CBSparsesym(nr::Integer)

initialize to zero-matrix of size nr*nr
"""
CBSparsesym(nr::Integer) = CBSparsesym(@ccall libcb.cb_sparsesym_new3(nr::Cint)::Ptr{Cvoid})

@doc raw"""
    CBSparsesym(nr::Integer, nz::Integer, ini::Union{<:AbstractVector{Cint},Nothing}, inj::Union{<:AbstractVector{Cint},Nothing}, va::Union{<:AbstractVector{Cdouble},Nothing})

initialize to size nr*nr and nz nonzeros so that this(ini[i],inj[i])=val[i] for i=0,..,nz-1; specify only one of (i,j) and (j,i), multiple elements are summed up.
"""
CBSparsesym(nr::Integer, nz::Integer, ini::Union{<:AbstractVector{Cint},Nothing}, inj::Union{<:AbstractVector{Cint},Nothing}, va::Union{<:AbstractVector{Cdouble},Nothing}) = GC.@preserve ini inj va begin
    (LinearAlgebra.chkstride1(va); LinearAlgebra.chkstride1(inj); LinearAlgebra.chkstride1(ini); CBSparsesym(@ccall libcb.cb_sparsesym_new4(nr::Cint, nz::Cint, ini::Ptr{Cint}, inj::Ptr{Cint}, va::Ptr{Cdouble})::Ptr{Cvoid}))
end

@doc raw"""
    CBSparsesym(nr::Integer, nz::Integer, ini::CBIndexmatrix, inj::CBIndexmatrix, va::CBMatrix)

initialize to size nr*nr and nz nonzeros so that this(ini(i),inj(i))=val[i] for i=0,..,nz-1; specify only one of (i,j) and (j,i), multiple elements are summed up.
"""
CBSparsesym(nr::Integer, nz::Integer, ini::CBIndexmatrix, inj::CBIndexmatrix, va::CBMatrix) = CBSparsesym(@ccall libcb.cb_sparsesym_new5(nr::Cint, nz::Cint, ini.data::Ptr{Cvoid}, inj.data::Ptr{Cvoid}, va.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_init!(self::CBSparsesym, param0::CBSparsesym, d::Real = 1.)

initialize to *this=A*d
"""
cb_init!(self::CBSparsesym, param0::CBSparsesym, d::Real = 1.) = (@ccall libcb.cb_sparsesym_init(self.data::Ptr{Cvoid}, param0.data::Ptr{Cvoid}, d::Cdouble)::Ptr{Cvoid}; return self)

@doc raw"""
    cb_init!(self::CBSparsesym, param0::CBMatrix, d::Real = 1.)

initialize to *this=d*(A+transpose(A))/2., abs(values)<tol are removed from the support
"""
cb_init!(self::CBSparsesym, param0::CBMatrix, d::Real = 1.) = (@ccall libcb.cb_sparsesym_init2(self.data::Ptr{Cvoid}, param0.data::Ptr{Cvoid}, d::Cdouble)::Ptr{Cvoid}; return self)

@doc raw"""
    cb_init!(self::CBSparsesym, param0::CBIndexmatrix, d::Real = 1.)

initialize to *this=d*(A+transpose(A))/2., zeros are removed from the support
"""
cb_init!(self::CBSparsesym, param0::CBIndexmatrix, d::Real = 1.) = (@ccall libcb.cb_sparsesym_init3(self.data::Ptr{Cvoid}, param0.data::Ptr{Cvoid}, d::Cdouble)::Ptr{Cvoid}; return self)

@doc raw"""
    cb_init!(self::CBSparsesym, param0::CBSymmatrix, d::Real = 1.)

initialize to *this=A*d, abs(values)<tol are removed from the support
"""
cb_init!(self::CBSparsesym, param0::CBSymmatrix, d::Real = 1.) = (@ccall libcb.cb_sparsesym_init4(self.data::Ptr{Cvoid}, param0.data::Ptr{Cvoid}, d::Cdouble)::Ptr{Cvoid}; return self)

@doc raw"""
    cb_init!(self::CBSparsesym, param0::CBSparsemat, d::Real = 1.)

initialize to *this=d*(A+transpose(A))/2.
"""
cb_init!(self::CBSparsesym, param0::CBSparsemat, d::Real = 1.) = (@ccall libcb.cb_sparsesym_init5(self.data::Ptr{Cvoid}, param0.data::Ptr{Cvoid}, d::Cdouble)::Ptr{Cvoid}; return self)

@doc raw"""
    cb_init!(self::CBSparsesym, nr::Integer)

initialize to zero-matrix of size nr*nr
"""
cb_init!(self::CBSparsesym, nr::Integer) = (@ccall libcb.cb_sparsesym_init6(self.data::Ptr{Cvoid}, nr::Cint)::Ptr{Cvoid}; return self)

@doc raw"""
    cb_init!(self::CBSparsesym, nr::Integer, nz::Integer, ini::Union{<:AbstractVector{Cint},Nothing}, inj::Union{<:AbstractVector{Cint},Nothing}, va::Union{<:AbstractVector{Cdouble},Nothing})

initialize to size nr*nr and nz nonzeros so that this(ini[i],inj[i])=val[i] for i=0,..,nz-1; specify only one of (i,j) and (j,i), multiple elements are summed up.
"""
cb_init!(self::CBSparsesym, nr::Integer, nz::Integer, ini::Union{<:AbstractVector{Cint},Nothing}, inj::Union{<:AbstractVector{Cint},Nothing}, va::Union{<:AbstractVector{Cdouble},Nothing}) = GC.@preserve ini inj va begin
    (LinearAlgebra.chkstride1(va); LinearAlgebra.chkstride1(inj); LinearAlgebra.chkstride1(ini); @ccall libcb.cb_sparsesym_init7(self.data::Ptr{Cvoid}, nr::Cint, nz::Cint, ini::Ptr{Cint}, inj::Ptr{Cint}, va::Ptr{Cdouble})::Ptr{Cvoid}; return self)
end

@doc raw"""
    cb_init!(self::CBSparsesym, nr::Integer, nz::Integer, ini::CBIndexmatrix, inj::CBIndexmatrix, va::CBMatrix)

initialize to size nr*nr and nz nonzeros so that this(ini(i),inj(i))=val[i] for i=0,..,nz-1; specify only one of (i,j) and (j,i), multiple elements are summed up.
"""
cb_init!(self::CBSparsesym, nr::Integer, nz::Integer, ini::CBIndexmatrix, inj::CBIndexmatrix, va::CBMatrix) = (@ccall libcb.cb_sparsesym_init8(self.data::Ptr{Cvoid}, nr::Cint, nz::Cint, ini.data::Ptr{Cvoid}, inj.data::Ptr{Cvoid}, va.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

@doc raw"""
    cb_init_support!(self::CBSparsesym, A::CBSparsesym, d::Real = 0.)

initialize to the same support as A but with constant value d; the same support will be generated even for d=0.
"""
cb_init_support!(self::CBSparsesym, A::CBSparsesym, d::Real = 0.) = (@ccall libcb.cb_sparsesym_init_support(self.data::Ptr{Cvoid}, A.data::Ptr{Cvoid}, d::Cdouble)::Ptr{Cvoid}; return self)

@doc raw"""
    cb_set_tol!(self::CBSparsesym, t::Real)

set tolerance for recognizing zero values to t
"""
cb_set_tol!(self::CBSparsesym, t::Real) = @ccall libcb.cb_sparsesym_set_tol(self.data::Ptr{Cvoid}, t::Cdouble)::Cvoid

@doc raw"""
    CBSparsesym(param0::CBMatrix, d::Real = 1.)

initialize to *this=d*(A+transpose(A))/2., abs(values)<tol are removed from the support
"""
CBSparsesym(param0::CBMatrix, d::Real = 1.) = CBSparsesym(@ccall libcb.cb_sparsesym_new6(param0.data::Ptr{Cvoid}, d::Cdouble)::Ptr{Cvoid})

@doc raw"""
    CBSparsesym(param0::CBIndexmatrix, d::Real = 1.)

initialize to *this=d*(A+transpose(A))/2., zeros are removed from the support
"""
CBSparsesym(param0::CBIndexmatrix, d::Real = 1.) = CBSparsesym(@ccall libcb.cb_sparsesym_new7(param0.data::Ptr{Cvoid}, d::Cdouble)::Ptr{Cvoid})

@doc raw"""
    CBSparsesym(param0::CBSymmatrix, d::Real = 1.)

initialize to *this=A*d, abs(values)<tol are removed from the support
"""
CBSparsesym(param0::CBSymmatrix, d::Real = 1.) = CBSparsesym(@ccall libcb.cb_sparsesym_new8(param0.data::Ptr{Cvoid}, d::Cdouble)::Ptr{Cvoid})

@doc raw"""
    CBSparsesym(param0::CBSparsemat, d::Real = 1.)

initialize to *this=d*(A+transpose(A))/2.
"""
CBSparsesym(param0::CBSparsemat, d::Real = 1.) = CBSparsesym(@ccall libcb.cb_sparsesym_new9(param0.data::Ptr{Cvoid}, d::Cdouble)::Ptr{Cvoid})

@doc raw"""
    cb_dim2(self::CBSparsesym)

returns the number of rows in _nr and the number of columns in _nc
"""
function cb_dim2(self::CBSparsesym)
    c = Ref{Int}()
    r = Ref{Int}()
    @ccall libcb.cb_sparsesym_dim(self.data::Ptr{Cvoid}, r::Ref{Int}, c::Ref{Int})::Cvoid
    return r[], c[]
end

@doc raw"""
    cb_dim(self::CBSparsesym)

returns the dimension rows * columns when the matrix is regarded as a vector
"""
cb_dim(self::CBSparsesym) = @ccall libcb.cb_sparsesym_dim2(self.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_rowdim(self::CBSparsesym)

returns the row dimension
"""
cb_rowdim(self::CBSparsesym) = @ccall libcb.cb_sparsesym_rowdim(self.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_coldim(self::CBSparsesym)

returns the column dimension
"""
cb_coldim(self::CBSparsesym) = @ccall libcb.cb_sparsesym_coldim(self.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_nonzeros(self::CBSparsesym)

returns the number of nonzeros in the lower triangle (including diagonal)
"""
cb_nonzeros(self::CBSparsesym) = @ccall libcb.cb_sparsesym_nonzeros(self.data::Ptr{Cvoid})::Cint

@doc raw"""
    Base.getindex(self::CBSparsesym, i::Integer, j::Integer)

returns value of element (i,j) of the matrix (rowindex i, columnindex j)
"""
Base.getindex(self::CBSparsesym, i::Integer, j::Integer) = @ccall libcb.cb_sparsesym_get(self.data::Ptr{Cvoid}, i::Cint, j::Cint)::Cdouble

@doc raw"""
    Base.getindex(self::CBSparsesym, i::Integer)

returns value of element (i) of the matrix if regarded as vector of stacked columns [element (i%rowdim, i/rowdim)]
"""
Base.getindex(self::CBSparsesym, i::Integer) = @ccall libcb.cb_sparsesym_get2(self.data::Ptr{Cvoid}, i::Cint)::Cdouble

@doc raw"""
    cb_get_colinfo(self::CBSparsesym)

returns information on nozero diagonal/columns, k by 4, listing: index (<0 for diagonal), # nonzeros, first index in colindex/colval, index in suppport submatrix
"""
cb_get_colinfo(self::CBSparsesym) = (@ccall libcb.cb_sparsesym_get_colinfo(self.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

@doc raw"""
    cb_get_colindex(self::CBSparsesym)

returns the index vector of the column representation holding the row index for each element
"""
cb_get_colindex(self::CBSparsesym) = (@ccall libcb.cb_sparsesym_get_colindex(self.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

@doc raw"""
    cb_get_colval(self::CBSparsesym)

returns the value vector of the column representation holding the value for each element
"""
cb_get_colval(self::CBSparsesym) = (@ccall libcb.cb_sparsesym_get_colval(self.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

@doc raw"""
    cb_get_suppind(self::CBSparsesym)

returns the index vector of the column representation holding the row index w.r.t. the principal support submatrix for each element
"""
cb_get_suppind(self::CBSparsesym) = (@ccall libcb.cb_sparsesym_get_suppind(self.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

@doc raw"""
    cb_get_suppcol(self::CBSparsesym)

returns the vector listing in ascending order the original column indices of the principal support submatrix
"""
cb_get_suppcol(self::CBSparsesym) = (@ccall libcb.cb_sparsesym_get_suppcol(self.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

@doc raw"""
    cb_get_edge_rep(self::CBSparsesym, I::CBIndexmatrix, J::CBIndexmatrix, val::CBMatrix)

stores the nz nonzero values of the lower triangle of *this in I,J,val so that this(I(i),J(i))=val(i) for i=0,...,nz-1 and dim(I)=dim(J)=dim(val)=nz (ordered as in row representation)
"""
cb_get_edge_rep(self::CBSparsesym, I::CBIndexmatrix, J::CBIndexmatrix, val::CBMatrix) = @ccall libcb.cb_sparsesym_get_edge_rep(self.data::Ptr{Cvoid}, I.data::Ptr{Cvoid}, J.data::Ptr{Cvoid}, val.data::Ptr{Cvoid})::Cvoid

@doc raw"""
    cb_contains_support(self::CBSparsesym, A::CBSparsesym)

returns 1 if A is of the same dimension and the support of A is contained in the support of *this, 0 otherwise
"""
cb_contains_support(self::CBSparsesym, A::CBSparsesym) = @ccall libcb.cb_sparsesym_contains_support(self.data::Ptr{Cvoid}, A.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_check_support(self::CBSparsesym, i::Integer, j::Integer)

returns 0 if (i,j) is not in the support, 1 otherwise
"""
cb_check_support(self::CBSparsesym, i::Integer, j::Integer) = @ccall libcb.cb_sparsesym_check_support(self.data::Ptr{Cvoid}, i::Cint, j::Cint)::Cint

@doc raw"""
    cb_diag(A::CBSparsesym)

returns the diagonal of A as a dense Matrix vector
"""
cb_diag(A::CBSparsesym) = CBMatrix(@ccall libcb.cb_sparsesym_new_diag(A.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_sparseDiag(A::CBMatrix, tol::Real)

forms a sparse symmetrix matrix having vector A on its diagonal
"""
cb_sparseDiag(A::CBMatrix, tol::Real) = CBSparsesym(@ccall libcb.cb_sparsesym_new_sparsediag(A.data::Ptr{Cvoid}, tol::Cdouble)::Ptr{Cvoid})

@doc raw"""
    cb_swap(A::CBSparsesym, B::CBSparsesym)

swap the content of the two sparse matrices A and B (involves no copying)
"""
cb_swap(A::CBSparsesym, B::CBSparsesym) = @ccall libcb.cb_sparsesym_swap(A.data::Ptr{Cvoid}, B.data::Ptr{Cvoid})::Cvoid

@doc raw"""
    cb_xeya!(self::CBSparsesym, A::CBSparsesym, d::Real = 1.)

sets *this=d*A and returns *this
"""
cb_xeya!(self::CBSparsesym, A::CBSparsesym, d::Real = 1.) = (@ccall libcb.cb_sparsesym_xeya(self.data::Ptr{Cvoid}, A.data::Ptr{Cvoid}, d::Cdouble)::Ptr{Cvoid}; return self)

@doc raw"""
    cb_xeya!(self::CBSparsesym, A::CBMatrix, d::Real = 1.)

sets and returns *this=d*(A+transpose(A))/2. where abs(values)<tol are removed from the support
"""
cb_xeya!(self::CBSparsesym, A::CBMatrix, d::Real = 1.) = (@ccall libcb.cb_sparsesym_xeya2(self.data::Ptr{Cvoid}, A.data::Ptr{Cvoid}, d::Cdouble)::Ptr{Cvoid}; return self)

@doc raw"""
    cb_xeya!(self::CBSparsesym, A::CBIndexmatrix, d::Real = 1.)

sets and returns *this=d*(A+transpose(A))/2. where zeros are removed from the support
"""
cb_xeya!(self::CBSparsesym, A::CBIndexmatrix, d::Real = 1.) = (@ccall libcb.cb_sparsesym_xeya3(self.data::Ptr{Cvoid}, A.data::Ptr{Cvoid}, d::Cdouble)::Ptr{Cvoid}; return self)

@doc raw"""
    cb_support_xbpeya!(self::CBSparsesym, y::CBSparsesym, alpha::Real = 1., beta::Real = 0.)

returns *this= alpha*y+beta*(*this) restricted to the curent support of *this; if beta==0, then *this is initialized to 0 on its support first
"""
cb_support_xbpeya!(self::CBSparsesym, y::CBSparsesym, alpha::Real = 1., beta::Real = 0.) = (@ccall libcb.cb_sparsesym_support_xbpeya(self.data::Ptr{Cvoid}, y.data::Ptr{Cvoid}, alpha::Cdouble, beta::Cdouble)::Ptr{Cvoid}; return self)

@doc raw"""
    cb_xbpeya(x::CBSparsesym, y::CBSparsesym, alpha::Real, beta::Real)

returns x= alpha*y+beta*x; if beta==0. then x is initialized to the correct size
"""
cb_xbpeya(x::CBSparsesym, y::CBSparsesym, alpha::Real, beta::Real) = (@ccall libcb.cb_sparsesym_xbpeya(x.data::Ptr{Cvoid}, y.data::Ptr{Cvoid}, alpha::Cdouble, beta::Cdouble)::Ptr{Cvoid}; return self)

@doc raw"""
    cb_xeyapzb(x::CBSparsesym, y::CBSparsesym, z::CBSparsesym, alpha::Real, beta::Real)

returns x= alpha*y+beta*z; x is initialized to the correct size
"""
cb_xeyapzb(x::CBSparsesym, y::CBSparsesym, z::CBSparsesym, alpha::Real, beta::Real) = (@ccall libcb.cb_sparsesym_xeyapzb(x.data::Ptr{Cvoid}, y.data::Ptr{Cvoid}, z.data::Ptr{Cvoid}, alpha::Cdouble, beta::Cdouble)::Ptr{Cvoid}; return self)

@doc raw"""
    cb_support_rankadd(A::CBMatrix, C::CBSparsesym, alpha::Real, beta::Real, trans::Integer)

returns C=beta*C+alpha*AA^T (or A^TA), but only on the current support of C
"""
cb_support_rankadd(A::CBMatrix, C::CBSparsesym, alpha::Real, beta::Real, trans::Integer) = (@ccall libcb.cb_sparsesym_support_rankadd(A.data::Ptr{Cvoid}, C.data::Ptr{Cvoid}, alpha::Cdouble, beta::Cdouble, trans::Cint)::Ptr{Cvoid}; return self)

Base.copy!(self::CBSparsesym, A::CBSparsesym) = (@ccall libcb.cb_sparsesym_assign(self.data::Ptr{Cvoid}, A.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

MA.operate!(::typeof(+), self::CBSparsesym, v::CBSparsesym) = (@ccall libcb.cb_sparsesym_plus(self.data::Ptr{Cvoid}, v.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

MA.promote_operation(::typeof(+), ::Type{<:CBSparsesym}, ::Type{<:CBSparsesym}) = CBSparsesym

MA.operate!(::typeof(-), self::CBSparsesym, v::CBSparsesym) = (@ccall libcb.cb_sparsesym_minus(self.data::Ptr{Cvoid}, v.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

MA.promote_operation(::typeof(-), ::Type{<:CBSparsesym}, ::Type{<:CBSparsesym}) = CBSparsesym

Base.:-(self::CBSparsesym) = CBSparsesym(@ccall libcb.cb_sparsesym_new_minus(self.data::Ptr{Cvoid})::Ptr{Cvoid})

MA.operate!(::typeof(*), self::CBSparsesym, d::Real) = (@ccall libcb.cb_sparsesym_times(self.data::Ptr{Cvoid}, d::Cdouble)::Ptr{Cvoid}; return self)

MA.promote_operation(::typeof(*), ::Type{<:CBSparsesym}, ::Type{<:Real}) = CBSparsesym

@doc raw"""
    MA.operate!(::typeof(/), self::CBSparsesym, d::Real)

ATTENTION: d is NOT checked for 0
"""
MA.operate!(::typeof(/), self::CBSparsesym, d::Real) = (@ccall libcb.cb_sparsesym_divide(self.data::Ptr{Cvoid}, d::Cdouble)::Ptr{Cvoid}; return self)

MA.promote_operation(::typeof(/), ::Type{<:CBSparsesym}, ::Type{<:Real}) = CBSparsesym

@doc raw"""
    Base.copy!(self::CBSparsesym, A::CBMatrix)

sets *this=(A+transpose(A))/2. removing abs(values)<tol; returns *this
"""
Base.copy!(self::CBSparsesym, A::CBMatrix) = (@ccall libcb.cb_sparsesym_assign2(self.data::Ptr{Cvoid}, A.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

@doc raw"""
    Base.copy!(self::CBSparsesym, A::CBIndexmatrix)

sets *this=(A+transpose(A))/2. removing zeros; returns *this
"""
Base.copy!(self::CBSparsesym, A::CBIndexmatrix) = (@ccall libcb.cb_sparsesym_assign3(self.data::Ptr{Cvoid}, A.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

@doc raw"""
    cb_transpose!(self::CBSparsesym)

transposes itself (at almost no cost)
"""
cb_transpose!(self::CBSparsesym) = (@ccall libcb.cb_sparsesym_transpose(self.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

Base.:+(A::CBSparsesym, B::CBSparsesym) = CBSparsesym(@ccall libcb.cb_sparsesym_new_plus(A.data::Ptr{Cvoid}, B.data::Ptr{Cvoid})::Ptr{Cvoid})

Base.:-(A::CBSparsesym, B::CBSparsesym) = CBSparsesym(@ccall libcb.cb_sparsesym_new_minus2(A.data::Ptr{Cvoid}, B.data::Ptr{Cvoid})::Ptr{Cvoid})

Base.:*(A::CBSparsesym, d::Real) = CBSparsesym(@ccall libcb.cb_sparsesym_new_times(A.data::Ptr{Cvoid}, d::Cdouble)::Ptr{Cvoid})

Base.:*(d::Real, A::CBSparsesym) = CBSparsesym(@ccall libcb.cb_sparsesym_new_times2(d::Cdouble, A.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    Base.:/(A::CBSparsesym, d::Real)

ATTENTION: d is NOT checked for 0
"""
Base.:/(A::CBSparsesym, d::Real) = CBSparsesym(@ccall libcb.cb_sparsesym_new_divide(A.data::Ptr{Cvoid}, d::Cdouble)::Ptr{Cvoid})

Base.:*(A::CBMatrix, B::CBSparsesym) = CBMatrix(@ccall libcb.cb_sparsesym_new_times3(A.data::Ptr{Cvoid}, B.data::Ptr{Cvoid})::Ptr{Cvoid})

Base.:*(A::CBSparsesym, B::CBMatrix) = CBMatrix(@ccall libcb.cb_sparsesym_new_times4(A.data::Ptr{Cvoid}, B.data::Ptr{Cvoid})::Ptr{Cvoid})

Base.:+(A::CBMatrix, B::CBSparsesym) = CBMatrix(@ccall libcb.cb_sparsesym_new_plus2(A.data::Ptr{Cvoid}, B.data::Ptr{Cvoid})::Ptr{Cvoid})

Base.:+(A::CBSparsesym, B::CBMatrix) = CBMatrix(@ccall libcb.cb_sparsesym_new_plus3(A.data::Ptr{Cvoid}, B.data::Ptr{Cvoid})::Ptr{Cvoid})

Base.:-(A::CBMatrix, B::CBSparsesym) = CBMatrix(@ccall libcb.cb_sparsesym_new_minus3(A.data::Ptr{Cvoid}, B.data::Ptr{Cvoid})::Ptr{Cvoid})

Base.:-(A::CBSparsesym, B::CBMatrix) = CBMatrix(@ccall libcb.cb_sparsesym_new_minus4(A.data::Ptr{Cvoid}, B.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_transpose(A::CBSparsesym)

(drop it or use a constructor instead)
"""
cb_transpose(A::CBSparsesym) = CBSparsesym(@ccall libcb.cb_sparsesym_new_transpose(A.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_xeya!(self::CBSparsesym, A::CBSymmatrix, d::Real = 1.)

sets and returns *this=A*d where abs(values)<tol are removed from the support
"""
cb_xeya!(self::CBSparsesym, A::CBSymmatrix, d::Real = 1.) = (@ccall libcb.cb_sparsesym_xeya4(self.data::Ptr{Cvoid}, A.data::Ptr{Cvoid}, d::Cdouble)::Ptr{Cvoid}; return self)

@doc raw"""
    cb_xeya!(self::CBSparsesym, A::CBSparsemat, d::Real = 1.)

sets and returns *this=d.*(A+transpose(A))/2.
"""
cb_xeya!(self::CBSparsesym, A::CBSparsemat, d::Real = 1.) = (@ccall libcb.cb_sparsesym_xeya5(self.data::Ptr{Cvoid}, A.data::Ptr{Cvoid}, d::Cdouble)::Ptr{Cvoid}; return self)

@doc raw"""
    Base.copy!(self::CBSparsesym, A::CBSymmatrix)

sets and returns *this=A where abs(values)<tol are removed from the support
"""
Base.copy!(self::CBSparsesym, A::CBSymmatrix) = (@ccall libcb.cb_sparsesym_assign4(self.data::Ptr{Cvoid}, A.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

@doc raw"""
    Base.copy!(self::CBSparsesym, A::CBSparsemat)

sets and returns *this=(A+transpose(A))/2.
"""
Base.copy!(self::CBSparsesym, A::CBSparsemat) = (@ccall libcb.cb_sparsesym_assign5(self.data::Ptr{Cvoid}, A.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

@doc raw"""
    cb_sparsemult(self::CBSparsesym, A::CBMatrix)

compute (*this)*A and return the result in a Sparsemat
"""
cb_sparsemult(self::CBSparsesym, A::CBMatrix) = CBSparsemat(@ccall libcb.cb_sparsesym_new_sparsemult(self.data::Ptr{Cvoid}, A.data::Ptr{Cvoid})::Ptr{Cvoid})

Base.:+(A::CBSparsesym, B::CBSymmatrix) = CBSymmatrix(@ccall libcb.cb_sparsesym_new_plus4(A.data::Ptr{Cvoid}, B.data::Ptr{Cvoid})::Ptr{Cvoid})

Base.:+(A::CBSymmatrix, B::CBSparsesym) = CBSymmatrix(@ccall libcb.cb_sparsesym_new_plus5(A.data::Ptr{Cvoid}, B.data::Ptr{Cvoid})::Ptr{Cvoid})

Base.:-(A::CBSparsesym, B::CBSymmatrix) = CBSymmatrix(@ccall libcb.cb_sparsesym_new_minus5(A.data::Ptr{Cvoid}, B.data::Ptr{Cvoid})::Ptr{Cvoid})

Base.:-(A::CBSymmatrix, B::CBSparsesym) = CBSymmatrix(@ccall libcb.cb_sparsesym_new_minus6(A.data::Ptr{Cvoid}, B.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_ip(A::CBSymmatrix, B::CBSparsesym)

returns the usual inner product of A and B, i.e., the sum of A(i,j)*B(i,j) over all i,j
"""
cb_ip(A::CBSymmatrix, B::CBSparsesym) = @ccall libcb.cb_sparsesym_ip(A.data::Ptr{Cvoid}, B.data::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_ip(A::CBSparsesym, B::CBSymmatrix)

returns the usual inner product of A and B, i.e., the sum of A(i,j)*B(i,j) over all i,j
"""
cb_ip(A::CBSparsesym, B::CBSymmatrix) = @ccall libcb.cb_sparsesym_ip2(A.data::Ptr{Cvoid}, B.data::Ptr{Cvoid})::Cdouble

Base.:*(A::CBSparsesym, B::CBSparsemat) = CBMatrix(@ccall libcb.cb_sparsesym_new_times5(A.data::Ptr{Cvoid}, B.data::Ptr{Cvoid})::Ptr{Cvoid})

Base.:*(A::CBSparsemat, B::CBSparsesym) = CBMatrix(@ccall libcb.cb_sparsesym_new_times6(A.data::Ptr{Cvoid}, B.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_abs(A::CBSparsesym)

returns a Sparsesym with elements abs((*this)(i,j)) for all i,j
"""
cb_abs(A::CBSparsesym) = CBSparsesym(@ccall libcb.cb_sparsesym_new_abs(A.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_trace(A::CBSparsesym)

returns the sum of the diagonal elements A(i,i) over all i
"""
cb_trace(A::CBSparsesym) = @ccall libcb.cb_sparsesym_trace(A.data::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_ip(A::CBSparsesym, B::CBSparsesym)

returns the usual inner product of A and B, i.e., the sum of A(i,j)*B(i,j) over all i,j
"""
cb_ip(A::CBSparsesym, B::CBSparsesym) = @ccall libcb.cb_sparsesym_ip3(A.data::Ptr{Cvoid}, B.data::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_ip(A::CBMatrix, B::CBSparsesym)

returns the usual inner product of A and B, i.e., the sum of A(i,j)*B(i,j) over all i,j
"""
cb_ip(A::CBMatrix, B::CBSparsesym) = @ccall libcb.cb_sparsesym_ip4(A.data::Ptr{Cvoid}, B.data::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_ip(A::CBSparsesym, B::CBMatrix)

returns the usual inner product of A and B, i.e., the sum of A(i,j)*B(i,j) over all i,j
"""
cb_ip(A::CBSparsesym, B::CBMatrix) = @ccall libcb.cb_sparsesym_ip5(A.data::Ptr{Cvoid}, B.data::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_norm2(A::CBSparsesym)

returns the Frobenius norm of A, i.e., the square root of the sum of A(i,j)*A(i,j) over all i,j
"""
cb_norm2(A::CBSparsesym) = @ccall libcb.cb_sparsesym_norm2(A.data::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_sumrows(A::CBSparsesym)

returns a row vector holding the sum over all rows, i.e., (1 1 ... 1)*A
"""
cb_sumrows(A::CBSparsesym) = CBMatrix(@ccall libcb.cb_sparsesym_new_sumrows(A.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_sumcols(A::CBSparsesym)

returns a column vector holding the sum over all columns, i.e., A*(1 1 ... 1)^T
"""
cb_sumcols(A::CBSparsesym) = CBMatrix(@ccall libcb.cb_sparsesym_new_sumcols(A.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_sum(A::CBSparsesym)

returns the sum over all elements of A, i.e., (1 1 ... 1)*A*(1 1 ... 1)^T
"""
cb_sum(A::CBSparsesym) = @ccall libcb.cb_sparsesym_sum(A.data::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_equal(A::CBSparsesym, B::CBSparsesym, eqtol::Real)

returns 1 if both matrices are identical, 0 otherwise
"""
cb_equal(A::CBSparsesym, B::CBSparsesym, eqtol::Real) = @ccall libcb.cb_sparsesym_equal(A.data::Ptr{Cvoid}, B.data::Ptr{Cvoid}, eqtol::Cdouble)::Cint

@doc raw"""
    cb_display(self::CBSparsesym, precision::Integer = 0, width::Integer = 0, screenwidth::Integer = 0)

* @brief displays a matrix in a pretty way for bounded screen widths; for variables of value zero default values are used.
      
"""
cb_display(self::CBSparsesym, precision::Integer = 0, width::Integer = 0, screenwidth::Integer = 0) = @ccall libcb.cb_sparsesym_display(self.data::Ptr{Cvoid}, precision::Cint, width::Cint, screenwidth::Cint)::Cvoid

