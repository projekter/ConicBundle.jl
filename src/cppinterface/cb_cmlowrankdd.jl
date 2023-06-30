@doc raw"""
    CBCMlowrankdd(Ain::CBMatrix, Bin::CBMatrix, cip::Union{<:CBCoeffmatInfo,Nothing} = nothing)

copy Ain, Bin and store the user information
"""
CBCMlowrankdd(Ain::CBMatrix, Bin::CBMatrix, cip::Union{<:CBCoeffmatInfo,Nothing} = nothing) = CBCMlowrankdd(@ccall libcb.cb_cmlowrankdd_new(Ain.data::Ptr{Cvoid}, Bin.data::Ptr{Cvoid}, (isnothing(cip) ? C_NULL : cip.data)::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_clone(self::CBCMlowrankdd)

makes an explicit copy of itself and returns a pointer to it
"""
cb_clone(self::CBCMlowrankdd) = CBCoeffmat(@ccall libcb.cb_cmlowrankdd_clone(self.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_dim(self::CBCMlowrankdd)

returns the order of the represented symmetric matrix
"""
cb_dim(self::CBCMlowrankdd) = @ccall libcb.cb_cmlowrankdd_dim(self.data::Ptr{Cvoid})::Cint

@doc raw"""
    Base.getindex(self::CBCMlowrankdd, i::Integer, j::Integer)

returns the value of the matrix element (i,j)
"""
Base.getindex(self::CBCMlowrankdd, i::Integer, j::Integer) = @ccall libcb.cb_cmlowrankdd_get(self.data::Ptr{Cvoid}, i::Cint, j::Cint)::Cdouble

@doc raw"""
    cb_make_symmatrix(self::CBCMlowrankdd, S::CBSymmatrix)

returns a dense symmetric constraint matrix
"""
cb_make_symmatrix(self::CBCMlowrankdd, S::CBSymmatrix) = @ccall libcb.cb_cmlowrankdd_make_symmatrix(self.data::Ptr{Cvoid}, S.data::Ptr{Cvoid})::Cvoid

@doc raw"""
    cb_norm(self::CBCMlowrankdd)

returns the Frobenius norm of the matrix
"""
cb_norm(self::CBCMlowrankdd) = @ccall libcb.cb_cmlowrankdd_norm(self.data::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_subspace(self::CBCMlowrankdd, P::CBMatrix)

delivers a new object on the heap corresponding to the matrix P^T(*this)P, the caller is responsible for deleting the object
"""
cb_subspace(self::CBCMlowrankdd, P::CBMatrix) = CBCoeffmat(@ccall libcb.cb_cmlowrankdd_subspace(self.data::Ptr{Cvoid}, P.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_multiply!(self::CBCMlowrankdd, d::Real)

multiply constraint permanentely by d; this is to allow scaling or sign changes in the constraints
"""
cb_multiply!(self::CBCMlowrankdd, d::Real) = @ccall libcb.cb_cmlowrankdd_multiply(self.data::Ptr{Cvoid}, d::Cdouble)::Cvoid

@doc raw"""
    cb_ip(self::CBCMlowrankdd, S::CBSymmatrix)

returns ip(*this,S)=trace(*this*S), the trace inner product
"""
cb_ip(self::CBCMlowrankdd, S::CBSymmatrix) = @ccall libcb.cb_cmlowrankdd_ip(self.data::Ptr{Cvoid}, S.data::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_gramip(self::CBCMlowrankdd, P::CBMatrix)

returns ip(*this,PP^T)=trace P^T(*this)P
"""
cb_gramip(self::CBCMlowrankdd, P::CBMatrix) = @ccall libcb.cb_cmlowrankdd_gramip(self.data::Ptr{Cvoid}, P.data::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_gramip(self::CBCMlowrankdd, P::CBMatrix, start_row::Integer, Lam::Union{<:CBMatrix,Nothing} = nothing)

returns ip(*this,QQ^T)=trace Q^T(*this)Q for Q=P.rows(start_row,start_row+dim-1)
"""
cb_gramip(self::CBCMlowrankdd, P::CBMatrix, start_row::Integer, Lam::Union{<:CBMatrix,Nothing} = nothing) = @ccall libcb.cb_cmlowrankdd_gramip2(self.data::Ptr{Cvoid}, P.data::Ptr{Cvoid}, start_row::Cint, (isnothing(Lam) ? C_NULL : Lam.data)::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_addmeto(self::CBCMlowrankdd, S::CBSymmatrix, d::Real = 1.)

computes S+=d*(*this);
"""
cb_addmeto(self::CBCMlowrankdd, S::CBSymmatrix, d::Real = 1.) = @ccall libcb.cb_cmlowrankdd_addmeto(self.data::Ptr{Cvoid}, S.data::Ptr{Cvoid}, d::Cdouble)::Cvoid

@doc raw"""
    cb_addprodto(self::CBCMlowrankdd, D::CBMatrix, C::CBMatrix, d::Real = 1.)

computes D+=d*(*this)*C
"""
cb_addprodto(self::CBCMlowrankdd, D::CBMatrix, C::CBMatrix, d::Real = 1.) = @ccall libcb.cb_cmlowrankdd_addprodto(self.data::Ptr{Cvoid}, D.data::Ptr{Cvoid}, C.data::Ptr{Cvoid}, d::Cdouble)::Cvoid

@doc raw"""
    cb_addprodto(self::CBCMlowrankdd, D::CBMatrix, C::CBSparsemat, d::Real = 1.)

computes D+=d*(*this)*C
"""
cb_addprodto(self::CBCMlowrankdd, D::CBMatrix, C::CBSparsemat, d::Real = 1.) = @ccall libcb.cb_cmlowrankdd_addprodto2(self.data::Ptr{Cvoid}, D.data::Ptr{Cvoid}, C.data::Ptr{Cvoid}, d::Cdouble)::Cvoid

@doc raw"""
    cb_left_right_prod(self::CBCMlowrankdd, P::CBMatrix, Q::CBMatrix, R::CBMatrix)

computes R=P^T*(*this)*Q
"""
cb_left_right_prod(self::CBCMlowrankdd, P::CBMatrix, Q::CBMatrix, R::CBMatrix) = @ccall libcb.cb_cmlowrankdd_left_right_prod(self.data::Ptr{Cvoid}, P.data::Ptr{Cvoid}, Q.data::Ptr{Cvoid}, R.data::Ptr{Cvoid})::Cvoid

@doc raw"""
    cb_prodvec_flops(self::CBCMlowrankdd)

returns an estimate of number of flops to compute addprodto for a vector
"""
cb_prodvec_flops(self::CBCMlowrankdd) = @ccall libcb.cb_cmlowrankdd_prodvec_flops(self.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_dense(self::CBCMlowrankdd)

returns 1 if its structure is as bad as its dense symmetric representation, otherwise 0
"""
cb_dense(self::CBCMlowrankdd) = @ccall libcb.cb_cmlowrankdd_dense(self.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_sparse(self::CBCMlowrankdd)

returns 0 if not sparse, otherwise 1
"""
cb_sparse(self::CBCMlowrankdd) = @ccall libcb.cb_cmlowrankdd_sparse(self.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_sparse(self::CBCMlowrankdd, param0::CBIndexmatrix, param1::CBIndexmatrix, param2::CBMatrix, param3::Real)

returns 0 if not sparse. If it is sparse it returns 1 and the nonzero structure in I,J and val, where val is multiplied by d. Only the upper triangle (including diagonal) is delivered
"""
cb_sparse(self::CBCMlowrankdd, param0::CBIndexmatrix, param1::CBIndexmatrix, param2::CBMatrix, param3::Real) = @ccall libcb.cb_cmlowrankdd_sparse2(self.data::Ptr{Cvoid}, param0.data::Ptr{Cvoid}, param1.data::Ptr{Cvoid}, param2.data::Ptr{Cvoid}, param3::Cdouble)::Cint

@doc raw"""
    cb_support_in(self::CBCMlowrankdd, param0::CBSparsesym)

returns 0 if the support of the costraint matrix is not contained in the support of the sparse symmetric matrix S, 1 if it is contained.
"""
cb_support_in(self::CBCMlowrankdd, param0::CBSparsesym) = @ccall libcb.cb_cmlowrankdd_support_in(self.data::Ptr{Cvoid}, param0.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_ip(self::CBCMlowrankdd, S::CBSparsesym)

returns the inner product of the constraint matrix with S
"""
cb_ip(self::CBCMlowrankdd, S::CBSparsesym) = @ccall libcb.cb_cmlowrankdd_ip2(self.data::Ptr{Cvoid}, S.data::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_project(self::CBCMlowrankdd, S::CBSymmatrix, P::CBMatrix)

computes S=P^T*(*this)*P
"""
cb_project(self::CBCMlowrankdd, S::CBSymmatrix, P::CBMatrix) = @ccall libcb.cb_cmlowrankdd_project(self.data::Ptr{Cvoid}, S.data::Ptr{Cvoid}, P.data::Ptr{Cvoid})::Cvoid

@doc raw"""
    cb_add_projection(self::CBCMlowrankdd, S::CBSymmatrix, P::CBMatrix, alpha::Real = 1., start_row::Integer = 0)

computes S+=Q^T(*this)Q for Q=P.rows(start_row,start_row+dim-1)
"""
cb_add_projection(self::CBCMlowrankdd, S::CBSymmatrix, P::CBMatrix, alpha::Real = 1., start_row::Integer = 0) = @ccall libcb.cb_cmlowrankdd_add_projection(self.data::Ptr{Cvoid}, S.data::Ptr{Cvoid}, P.data::Ptr{Cvoid}, alpha::Cdouble, start_row::Cint)::Cvoid

@doc raw"""
    cb_postgenmult(self::CBCMlowrankdd, D::CBMatrix, C::CBMatrix, alpha::Real = 1., beta::Real = 0., dtrans::Integer = 0)

computes C= alpha*(*this)*D^(T if dtrans) + beta*C, C is also returned
"""
cb_postgenmult(self::CBCMlowrankdd, D::CBMatrix, C::CBMatrix, alpha::Real = 1., beta::Real = 0., dtrans::Integer = 0) = (@ccall libcb.cb_cmlowrankdd_postgenmult(self.data::Ptr{Cvoid}, D.data::Ptr{Cvoid}, C.data::Ptr{Cvoid}, alpha::Cdouble, beta::Cdouble, dtrans::Cint)::Ptr{Cvoid}; return self)

@doc raw"""
    cb_pregenmult(self::CBCMlowrankdd, D::CBMatrix, C::CBMatrix, alpha::Real = 1., beta::Real = 0., dtrans::Integer = 0)

computes C= alpha*D^(T if dtrans)*(*this) + beta*C, C is also returned
"""
cb_pregenmult(self::CBCMlowrankdd, D::CBMatrix, C::CBMatrix, alpha::Real = 1., beta::Real = 0., dtrans::Integer = 0) = (@ccall libcb.cb_cmlowrankdd_pregenmult(self.data::Ptr{Cvoid}, D.data::Ptr{Cvoid}, C.data::Ptr{Cvoid}, alpha::Cdouble, beta::Cdouble, dtrans::Cint)::Ptr{Cvoid}; return self)

@doc raw"""
    cb_equal(self::CBCMlowrankdd, p::Union{<:CBCoeffmat,Nothing}, tol::Real = 1e-6)

returns 1, if p is the same derived class and entries differ by less than tol, otherwise zero
"""
cb_equal(self::CBCMlowrankdd, p::Union{<:CBCoeffmat,Nothing}, tol::Real = 1e-6) = @ccall libcb.cb_cmlowrankdd_equal(self.data::Ptr{Cvoid}, (isnothing(p) ? C_NULL : p.data)::Ptr{Cvoid}, tol::Cdouble)::Cint

@doc raw"""
    cb_display(self::CBCMlowrankdd)

display constraint information
"""
cb_display(self::CBCMlowrankdd) = @ccall libcb.cb_cmlowrankdd_display(self.data::Ptr{Cvoid})::Cvoid

@doc raw"""
    cb_out(self::CBCMlowrankdd)

put entire contents onto ostream with the class type in the beginning so that the derived class can be recognized by in().
"""
cb_out(self::CBCMlowrankdd) = @ccall libcb.cb_cmlowrankdd_out(self.data::Ptr{Cvoid})::Cvoid

