@doc raw"""
    CBCMsymdense(Ain::CBSymmatrix, cip::Union{<:CBCoeffmatInfo,Nothing} = nothing)

copy Ain and possibly store the user information
"""
CBCMsymdense(Ain::CBSymmatrix, cip::Union{<:CBCoeffmatInfo,Nothing} = nothing) = CBCMsymdense(@ccall libcb.cb_cmsymdense_new(Ain.data::Ptr{Cvoid}, (isnothing(cip) ? C_NULL : cip.data)::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_clone(self::CBCMsymdense)

makes an explicit copy of itself and returns a pointer to it
"""
cb_clone(self::CBCMsymdense) = CBCoeffmat(@ccall libcb.cb_cmsymdense_clone(self.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_dim(self::CBCMsymdense)

returns the order of the represented symmetric matrix
"""
cb_dim(self::CBCMsymdense) = @ccall libcb.cb_cmsymdense_dim(self.data::Ptr{Cvoid})::Cint

@doc raw"""
    Base.getindex(self::CBCMsymdense, i::Integer, j::Integer)

returns the value of the matrix element (i,j)
"""
Base.getindex(self::CBCMsymdense, i::Integer, j::Integer) = @ccall libcb.cb_cmsymdense_get(self.data::Ptr{Cvoid}, i::Cint, j::Cint)::Cdouble

@doc raw"""
    cb_make_symmatrix(self::CBCMsymdense, S::CBSymmatrix)

returns a dense symmetric constraint matrix
"""
cb_make_symmatrix(self::CBCMsymdense, S::CBSymmatrix) = @ccall libcb.cb_cmsymdense_make_symmatrix(self.data::Ptr{Cvoid}, S.data::Ptr{Cvoid})::Cvoid

@doc raw"""
    cb_norm(self::CBCMsymdense)

returns the Frobenius norm of the matrix
"""
cb_norm(self::CBCMsymdense) = @ccall libcb.cb_cmsymdense_norm(self.data::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_subspace(self::CBCMsymdense, P::CBMatrix)

delivers a new object on the heap corresponding to the matrix P^T(*this)P, the caller is responsible for deleting the object
"""
cb_subspace(self::CBCMsymdense, P::CBMatrix) = CBCoeffmat(@ccall libcb.cb_cmsymdense_subspace(self.data::Ptr{Cvoid}, P.data::Ptr{Cvoid})::Ptr{Cvoid})

@doc raw"""
    cb_multiply!(self::CBCMsymdense, d::Real)

multiply constraint permanentely by d; this is to allow scaling or sign changes in the constraints
"""
cb_multiply!(self::CBCMsymdense, d::Real) = @ccall libcb.cb_cmsymdense_multiply(self.data::Ptr{Cvoid}, d::Cdouble)::Cvoid

@doc raw"""
    cb_ip(self::CBCMsymdense, S::CBSymmatrix)

returns ip(*this,S)=trace(*this*S), the trace inner product
"""
cb_ip(self::CBCMsymdense, S::CBSymmatrix) = @ccall libcb.cb_cmsymdense_ip(self.data::Ptr{Cvoid}, S.data::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_gramip(self::CBCMsymdense, P::CBMatrix)

returns ip(*this,PP^T)=trace P^T(*this)P
"""
cb_gramip(self::CBCMsymdense, P::CBMatrix) = @ccall libcb.cb_cmsymdense_gramip(self.data::Ptr{Cvoid}, P.data::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_gramip(self::CBCMsymdense, P::CBMatrix, start_row::Integer, Lam::Union{<:CBMatrix,Nothing} = nothing)

returns ip(*this,QQ^T)=trace Q^T(*this)Q for Q=P.rows(start_row,start_row+dim-1)
"""
cb_gramip(self::CBCMsymdense, P::CBMatrix, start_row::Integer, Lam::Union{<:CBMatrix,Nothing} = nothing) = @ccall libcb.cb_cmsymdense_gramip2(self.data::Ptr{Cvoid}, P.data::Ptr{Cvoid}, start_row::Cint, (isnothing(Lam) ? C_NULL : Lam.data)::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_addmeto(self::CBCMsymdense, S::CBSymmatrix, d::Real = 1.)

computes S+=d*(*this);
"""
cb_addmeto(self::CBCMsymdense, S::CBSymmatrix, d::Real = 1.) = @ccall libcb.cb_cmsymdense_addmeto(self.data::Ptr{Cvoid}, S.data::Ptr{Cvoid}, d::Cdouble)::Cvoid

@doc raw"""
    cb_addprodto(self::CBCMsymdense, B::CBMatrix, C::CBMatrix, d::Real = 1.)

computes B+=d*(*this)*C
"""
cb_addprodto(self::CBCMsymdense, B::CBMatrix, C::CBMatrix, d::Real = 1.) = @ccall libcb.cb_cmsymdense_addprodto(self.data::Ptr{Cvoid}, B.data::Ptr{Cvoid}, C.data::Ptr{Cvoid}, d::Cdouble)::Cvoid

@doc raw"""
    cb_addprodto(self::CBCMsymdense, B::CBMatrix, C::CBSparsemat, d::Real = 1.)

computes B+=d*(*this)*C
"""
cb_addprodto(self::CBCMsymdense, B::CBMatrix, C::CBSparsemat, d::Real = 1.) = @ccall libcb.cb_cmsymdense_addprodto2(self.data::Ptr{Cvoid}, B.data::Ptr{Cvoid}, C.data::Ptr{Cvoid}, d::Cdouble)::Cvoid

@doc raw"""
    cb_left_right_prod(self::CBCMsymdense, P::CBMatrix, Q::CBMatrix, R::CBMatrix)

computes R=P^T*(*this)*Q
"""
cb_left_right_prod(self::CBCMsymdense, P::CBMatrix, Q::CBMatrix, R::CBMatrix) = @ccall libcb.cb_cmsymdense_left_right_prod(self.data::Ptr{Cvoid}, P.data::Ptr{Cvoid}, Q.data::Ptr{Cvoid}, R.data::Ptr{Cvoid})::Cvoid

@doc raw"""
    cb_prodvec_flops(self::CBCMsymdense)

returns an estimate of number of flops to compute addprodto for a vector
"""
cb_prodvec_flops(self::CBCMsymdense) = @ccall libcb.cb_cmsymdense_prodvec_flops(self.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_dense(self::CBCMsymdense)

returns 1 if its structure is as bad as its dense symmetric representation, otherwise 0
"""
cb_dense(self::CBCMsymdense) = @ccall libcb.cb_cmsymdense_dense(self.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_sparse(self::CBCMsymdense)

returns 0 if not sparse, otherwise 1
"""
cb_sparse(self::CBCMsymdense) = @ccall libcb.cb_cmsymdense_sparse(self.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_sparse(self::CBCMsymdense, param0::CBIndexmatrix, param1::CBIndexmatrix, param2::CBMatrix, param3::Real)

returns 0 if not sparse. If it is sparse it returns 1 and the nonzero structure in I,J and v, where v is multiplied by d. Only the upper triangle (including diagonal) is delivered
"""
cb_sparse(self::CBCMsymdense, param0::CBIndexmatrix, param1::CBIndexmatrix, param2::CBMatrix, param3::Real) = @ccall libcb.cb_cmsymdense_sparse2(self.data::Ptr{Cvoid}, param0.data::Ptr{Cvoid}, param1.data::Ptr{Cvoid}, param2.data::Ptr{Cvoid}, param3::Cdouble)::Cint

@doc raw"""
    cb_support_in(self::CBCMsymdense, param0::CBSparsesym)

returns 0 if the support of the costraint matrix is not contained in the support of the sparse symmetric matrix S, 1 if it is contained.
"""
cb_support_in(self::CBCMsymdense, param0::CBSparsesym) = @ccall libcb.cb_cmsymdense_support_in(self.data::Ptr{Cvoid}, param0.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_ip(self::CBCMsymdense, S::CBSparsesym)

returns the inner product of the constraint matrix with S
"""
cb_ip(self::CBCMsymdense, S::CBSparsesym) = @ccall libcb.cb_cmsymdense_ip2(self.data::Ptr{Cvoid}, S.data::Ptr{Cvoid})::Cdouble

@doc raw"""
    cb_project(self::CBCMsymdense, S::CBSymmatrix, P::CBMatrix)

computes S=P^T*(*this)*P
"""
cb_project(self::CBCMsymdense, S::CBSymmatrix, P::CBMatrix) = @ccall libcb.cb_cmsymdense_project(self.data::Ptr{Cvoid}, S.data::Ptr{Cvoid}, P.data::Ptr{Cvoid})::Cvoid

@doc raw"""
    cb_add_projection(self::CBCMsymdense, S::CBSymmatrix, P::CBMatrix, alpha::Real = 1., start_row::Integer = 0)

computes S+=Q^T(*this)Q for Q=P.rows(start_row,start_row+dim-1)
"""
cb_add_projection(self::CBCMsymdense, S::CBSymmatrix, P::CBMatrix, alpha::Real = 1., start_row::Integer = 0) = @ccall libcb.cb_cmsymdense_add_projection(self.data::Ptr{Cvoid}, S.data::Ptr{Cvoid}, P.data::Ptr{Cvoid}, alpha::Cdouble, start_row::Cint)::Cvoid

@doc raw"""
    cb_postgenmult(self::CBCMsymdense, B::CBMatrix, C::CBMatrix, alpha::Real = 1., beta::Real = 0., btrans::Integer = 0)

computes C= alpha*(*this)*B^(T if btrans) + beta*C, C is also returned
"""
cb_postgenmult(self::CBCMsymdense, B::CBMatrix, C::CBMatrix, alpha::Real = 1., beta::Real = 0., btrans::Integer = 0) = (@ccall libcb.cb_cmsymdense_postgenmult(self.data::Ptr{Cvoid}, B.data::Ptr{Cvoid}, C.data::Ptr{Cvoid}, alpha::Cdouble, beta::Cdouble, btrans::Cint)::Ptr{Cvoid}; return self)

@doc raw"""
    cb_pregenmult(self::CBCMsymdense, B::CBMatrix, C::CBMatrix, alpha::Real = 1., beta::Real = 0., btrans::Integer = 0)

computes C= alpha*B^(T if btrans)*(*this) + beta*C, C is also returned
"""
cb_pregenmult(self::CBCMsymdense, B::CBMatrix, C::CBMatrix, alpha::Real = 1., beta::Real = 0., btrans::Integer = 0) = (@ccall libcb.cb_cmsymdense_pregenmult(self.data::Ptr{Cvoid}, B.data::Ptr{Cvoid}, C.data::Ptr{Cvoid}, alpha::Cdouble, beta::Cdouble, btrans::Cint)::Ptr{Cvoid}; return self)

@doc raw"""
    cb_equal(self::CBCMsymdense, p::Union{<:CBCoeffmat,Nothing}, tol::Real = 1e-6)

returns 1, if p is the same derived class and entries differ by less than tol, otherwise zero
"""
cb_equal(self::CBCMsymdense, p::Union{<:CBCoeffmat,Nothing}, tol::Real = 1e-6) = @ccall libcb.cb_cmsymdense_equal(self.data::Ptr{Cvoid}, (isnothing(p) ? C_NULL : p.data)::Ptr{Cvoid}, tol::Cdouble)::Cint

@doc raw"""
    cb_display(self::CBCMsymdense)

display constraint information
"""
cb_display(self::CBCMsymdense) = @ccall libcb.cb_cmsymdense_display(self.data::Ptr{Cvoid})::Cvoid

@doc raw"""
    cb_out(self::CBCMsymdense)

put entire contents onto outstream with the class type in the beginning so that the derived class can be recognized by in().
"""
cb_out(self::CBCMsymdense) = @ccall libcb.cb_cmsymdense_out(self.data::Ptr{Cvoid})::Cvoid

@doc raw"""
    cb_get_A(self::CBCMsymdense)

returns the const reference to the internal symmetric matrix
"""
cb_get_A(self::CBCMsymdense) = (@ccall libcb.cb_cmsymdense_get_a(self.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

