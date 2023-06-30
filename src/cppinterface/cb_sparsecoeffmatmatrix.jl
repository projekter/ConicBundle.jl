@doc raw"""
    Base.copy!(self::CBSparseCoeffmatMatrix, param0::CBSparseCoeffmatMatrix)

copy
"""
Base.copy!(self::CBSparseCoeffmatMatrix, param0::CBSparseCoeffmatMatrix) = (@ccall libcb.cb_sparsecoeffmatmatrix_assign(self.data::Ptr{Cvoid}, param0.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

@doc raw"""
    cb_clear!(self::CBSparseCoeffmatMatrix)

clears all (empty 0 times 0 matrix)
"""
cb_clear!(self::CBSparseCoeffmatMatrix) = @ccall libcb.cb_sparsecoeffmatmatrix_clear(self.data::Ptr{Cvoid})::Cvoid

@doc raw"""
    cb_init!(self::CBSparseCoeffmatMatrix, block_dim::CBIndexmatrix, col_dim::Integer, block_ind::Union{<:CBIndexmatrix,Nothing} = nothing, col_ind::Union{<:CBIndexmatrix,Nothing} = nothing, coeff_vec::Union{<:CBCoeffmatVector,Nothing} = nothing)

first calls clear() and then it sets the new values (if one of block_ind, col_ind, or coeff_vec is !=NULL, all must be !=NULL and of the same size)
"""
cb_init!(self::CBSparseCoeffmatMatrix, block_dim::CBIndexmatrix, col_dim::Integer, block_ind::Union{<:CBIndexmatrix,Nothing} = nothing, col_ind::Union{<:CBIndexmatrix,Nothing} = nothing, coeff_vec::Union{<:CBCoeffmatVector,Nothing} = nothing) = @ccall libcb.cb_sparsecoeffmatmatrix_init(self.data::Ptr{Cvoid}, block_dim.data::Ptr{Cvoid}, col_dim::Cint, (isnothing(block_ind) ? C_NULL : block_ind.data)::Ptr{Cvoid}, (isnothing(col_ind) ? C_NULL : col_ind.data)::Ptr{Cvoid}, (isnothing(coeff_vec) ? C_NULL : coeff_vec.data)::Ptr{Cvoid})::Cint

@doc raw"""
    CBSparseCoeffmatMatrix(incr::Integer = -1)

set the output and call clear()
"""
CBSparseCoeffmatMatrix(incr::Integer = -1) = CBSparseCoeffmatMatrix(@ccall libcb.cb_sparsecoeffmatmatrix_new(incr::Cint)::Ptr{Cvoid})

@doc raw"""
    CBSparseCoeffmatMatrix(S::CBSparseCoeffmatMatrix, incr::Integer = -1)

set the output and call clear()
"""
CBSparseCoeffmatMatrix(S::CBSparseCoeffmatMatrix, incr::Integer = -1) = CBSparseCoeffmatMatrix(@ccall libcb.cb_sparsecoeffmatmatrix_new2(S.data::Ptr{Cvoid}, incr::Cint)::Ptr{Cvoid})

@doc raw"""
    CBSparseCoeffmatMatrix(in_block_dim::CBIndexmatrix, in_col_dim::Integer, block_ind::Union{<:CBIndexmatrix,Nothing} = nothing, col_ind::Union{<:CBIndexmatrix,Nothing} = nothing, coeff_vec::Union{<:CBCoeffmatVector,Nothing} = nothing, incr::Integer = -1)

set the output and call init() for the given sparse information
"""
CBSparseCoeffmatMatrix(in_block_dim::CBIndexmatrix, in_col_dim::Integer, block_ind::Union{<:CBIndexmatrix,Nothing} = nothing, col_ind::Union{<:CBIndexmatrix,Nothing} = nothing, coeff_vec::Union{<:CBCoeffmatVector,Nothing} = nothing, incr::Integer = -1) = CBSparseCoeffmatMatrix(@ccall libcb.cb_sparsecoeffmatmatrix_new3(in_block_dim.data::Ptr{Cvoid}, in_col_dim::Cint, (isnothing(block_ind) ? C_NULL : block_ind.data)::Ptr{Cvoid}, (isnothing(col_ind) ? C_NULL : col_ind.data)::Ptr{Cvoid}, (isnothing(coeff_vec) ? C_NULL : coeff_vec.data)::Ptr{Cvoid}, incr::Cint)::Ptr{Cvoid})

@doc raw"""
    cb_blockdim(self::CBSparseCoeffmatMatrix)

returns a column vector with row i giving the order of block i
"""
cb_blockdim(self::CBSparseCoeffmatMatrix) = (@ccall libcb.cb_sparsecoeffmatmatrix_blockdim(self.data::Ptr{Cvoid})::Ptr{Cvoid}; return self)

@doc raw"""
    cb_blockdim(self::CBSparseCoeffmatMatrix, i::Integer)

returns the order of block i
"""
cb_blockdim(self::CBSparseCoeffmatMatrix, i::Integer) = @ccall libcb.cb_sparsecoeffmatmatrix_blockdim2(self.data::Ptr{Cvoid}, i::Cint)::Cint

@doc raw"""
    cb_coldim(self::CBSparseCoeffmatMatrix)

returns the number of columns (i.e. of blockdiagonal matrices)
"""
cb_coldim(self::CBSparseCoeffmatMatrix) = @ccall libcb.cb_sparsecoeffmatmatrix_coldim(self.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_rowdim(self::CBSparseCoeffmatMatrix)

returns the number of blocks
"""
cb_rowdim(self::CBSparseCoeffmatMatrix) = @ccall libcb.cb_sparsecoeffmatmatrix_rowdim(self.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_nzcoldim(self::CBSparseCoeffmatMatrix)

returns the number of columns with nonzero coefficient matrices
"""
cb_nzcoldim(self::CBSparseCoeffmatMatrix) = @ccall libcb.cb_sparsecoeffmatmatrix_nzcoldim(self.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_set!(self::CBSparseCoeffmatMatrix, i::Integer, j::Integer, cm::CBCoeffmatPointer)

sets the CoeffmatPointer of block i in column j (blockdiagonal matrix j) (may be empty)
"""
cb_set!(self::CBSparseCoeffmatMatrix, i::Integer, j::Integer, cm::CBCoeffmatPointer) = @ccall libcb.cb_sparsecoeffmatmatrix_set(self.data::Ptr{Cvoid}, i::Cint, j::Cint, cm.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_set!(self::CBSparseCoeffmatMatrix, i::Integer, j::Integer, cm::Union{<:CBCoeffmat,Nothing})

sets the CoeffmatPointer of block i in column j (blockdiagonal matrix j) to point to cm (or deletes it if cm==0)
"""
cb_set!(self::CBSparseCoeffmatMatrix, i::Integer, j::Integer, cm::Union{<:CBCoeffmat,Nothing}) = @ccall libcb.cb_sparsecoeffmatmatrix_set2(self.data::Ptr{Cvoid}, i::Cint, j::Cint, (isnothing(cm) ? C_NULL : cm.data)::Ptr{Cvoid})::Cint

@doc raw"""
    Base.getindex(self::CBSparseCoeffmatMatrix, i::Integer, j::Integer)

returns the CoeffmatPointer of block i in column j (blockdiagonal matrix j) (may be empty)
"""
Base.getindex(self::CBSparseCoeffmatMatrix, i::Integer, j::Integer) = CBCoeffmatPointer(@ccall libcb.cb_sparsecoeffmatmatrix_new_get(self.data::Ptr{Cvoid}, i::Cint, j::Cint)::Ptr{Cvoid})

@doc raw"""
    cb_append_blocks!(self::CBSparseCoeffmatMatrix, append_mat::CBSparseCoeffmatMatrix, blocks::Union{<:CBIndexmatrix,Nothing} = nothing, cols::Union{<:CBIndexmatrix,Nothing} = nothing)

* @brief append append_mat (or its submatrix given by blocks and/or cols) below

        If blocks or cols have negative entries, each of these represent a row of blocks of zeros whose order is the negative of the given value or a column of zeros of appropriate size (the size of the negative value does not matter for columns)
     
"""
cb_append_blocks!(self::CBSparseCoeffmatMatrix, append_mat::CBSparseCoeffmatMatrix, blocks::Union{<:CBIndexmatrix,Nothing} = nothing, cols::Union{<:CBIndexmatrix,Nothing} = nothing) = @ccall libcb.cb_sparsecoeffmatmatrix_append_blocks(self.data::Ptr{Cvoid}, append_mat.data::Ptr{Cvoid}, (isnothing(blocks) ? C_NULL : blocks.data)::Ptr{Cvoid}, (isnothing(cols) ? C_NULL : cols.data)::Ptr{Cvoid})::Cint

@doc raw"""
    cb_append_columns!(self::CBSparseCoeffmatMatrix, append_mat::CBSparseCoeffmatMatrix, blocks::Union{<:CBIndexmatrix,Nothing} = nothing, cols::Union{<:CBIndexmatrix,Nothing} = nothing)

* @brief append append_mat (or its submatrix given by blocks and/or cols) to the right

        If blocks or cols have negative entries, each of these represent a row of blocks of zeros whose order is the negative of the given value (if columns exist already, these need to match the order of the existing blocks) or a column of zeros of appropriate size (the size of the negative value does not matter for columns)
     
"""
cb_append_columns!(self::CBSparseCoeffmatMatrix, append_mat::CBSparseCoeffmatMatrix, blocks::Union{<:CBIndexmatrix,Nothing} = nothing, cols::Union{<:CBIndexmatrix,Nothing} = nothing) = @ccall libcb.cb_sparsecoeffmatmatrix_append_columns(self.data::Ptr{Cvoid}, append_mat.data::Ptr{Cvoid}, (isnothing(blocks) ? C_NULL : blocks.data)::Ptr{Cvoid}, (isnothing(cols) ? C_NULL : cols.data)::Ptr{Cvoid})::Cint

@doc raw"""
    cb_reassign_blocks!(self::CBSparseCoeffmatMatrix, map_to_old::CBIndexmatrix)

afterwards the new block i is the previous block map_to_old(i); no multiple appearances are allowed, but not all have to appear (these are deleted)
"""
cb_reassign_blocks!(self::CBSparseCoeffmatMatrix, map_to_old::CBIndexmatrix) = @ccall libcb.cb_sparsecoeffmatmatrix_reassign_blocks(self.data::Ptr{Cvoid}, map_to_old.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_reassign_columns!(self::CBSparseCoeffmatMatrix, map_to_old::CBIndexmatrix)

afterwards the new column i is the previous column map_to_old(i); no multiple appearances are allowed, but not all have to appear (these are deleted)
"""
cb_reassign_columns!(self::CBSparseCoeffmatMatrix, map_to_old::CBIndexmatrix) = @ccall libcb.cb_sparsecoeffmatmatrix_reassign_columns(self.data::Ptr{Cvoid}, map_to_old.data::Ptr{Cvoid})::Cint

@doc raw"""
    cb_delete_blocks!(self::CBSparseCoeffmatMatrix, delete_indices::CBIndexmatrix, map_to_old::Union{<:CBIndexmatrix,Nothing} = nothing)

generates a map_to_old by deleting the specified indices and calls reassign_blocks; this map_to_old will be returned if the corresponding pointer is not NULL
"""
cb_delete_blocks!(self::CBSparseCoeffmatMatrix, delete_indices::CBIndexmatrix, map_to_old::Union{<:CBIndexmatrix,Nothing} = nothing) = @ccall libcb.cb_sparsecoeffmatmatrix_delete_blocks(self.data::Ptr{Cvoid}, delete_indices.data::Ptr{Cvoid}, (isnothing(map_to_old) ? C_NULL : map_to_old.data)::Ptr{Cvoid})::Cint

@doc raw"""
    cb_delete_columns!(self::CBSparseCoeffmatMatrix, delete_indices::CBIndexmatrix, map_to_old::Union{<:CBIndexmatrix,Nothing} = nothing)

generates a map_to_old by deleting the specified indices and calls reassign_columns; this map_to_old will be returned if the corresponding pointer is not NULL
"""
cb_delete_columns!(self::CBSparseCoeffmatMatrix, delete_indices::CBIndexmatrix, map_to_old::Union{<:CBIndexmatrix,Nothing} = nothing) = @ccall libcb.cb_sparsecoeffmatmatrix_delete_columns(self.data::Ptr{Cvoid}, delete_indices.data::Ptr{Cvoid}, (isnothing(map_to_old) ? C_NULL : map_to_old.data)::Ptr{Cvoid})::Cint

@doc raw"""
    cb_get_dense_cnt(self::CBSparseCoeffmatMatrix, i::Integer)

returns the number of dense matrices in block i
"""
cb_get_dense_cnt(self::CBSparseCoeffmatMatrix, i::Integer) = @ccall libcb.cb_sparsecoeffmatmatrix_get_dense_cnt(self.data::Ptr{Cvoid}, i::Cint)::Cint

@doc raw"""
    Base.:(==)(self::CBSparseCoeffmatMatrix, mat::CBSparseCoeffmatMatrix)

useful for testing purposes; true iff both have the same size and both have the same pointers in the same positions
"""
Base.:(==)(self::CBSparseCoeffmatMatrix, mat::CBSparseCoeffmatMatrix) = Bool(@ccall libcb.cb_sparsecoeffmatmatrix_new_equal(self.data::Ptr{Cvoid}, mat.data::Ptr{Cvoid})::Bool)

@doc raw"""
    cb_Gram_ip(self::CBSparseCoeffmatMatrix, ipvec::CBMatrix, P::CBMatrix, Lam::Union{<:CBMatrix,Nothing} = nothing, ind::Union{<:CBIndexmatrix,Nothing} = nothing)

computes the inner products of (selected) columns (which represent block diagonal symmetric matrices) with the Gram matrix P*P^T (or, if Lam is given, P*Diag(Lam)*P^T) into the column vector ipvec(j)=ip(P*P^T,A.column((*ind)(j)}) (j=0,...,ind->dim()-1); if ind==NULL, use all columns
"""
cb_Gram_ip(self::CBSparseCoeffmatMatrix, ipvec::CBMatrix, P::CBMatrix, Lam::Union{<:CBMatrix,Nothing} = nothing, ind::Union{<:CBIndexmatrix,Nothing} = nothing) = @ccall libcb.cb_sparsecoeffmatmatrix_gram_ip(self.data::Ptr{Cvoid}, ipvec.data::Ptr{Cvoid}, P.data::Ptr{Cvoid}, (isnothing(Lam) ? C_NULL : Lam.data)::Ptr{Cvoid}, (isnothing(ind) ? C_NULL : ind.data)::Ptr{Cvoid})::Cint

@doc raw"""
    cb_Gram_ip(self::CBSparseCoeffmatMatrix, P::CBMatrix, j::Integer)

computes the inner product of the block diagonal symmetric matrix stored in colummn j with the Gram matrix P*P^T into ipval=ip(P*P^T,A.column(j))
"""
function cb_Gram_ip(self::CBSparseCoeffmatMatrix, P::CBMatrix, j::Integer)
    ipval = Ref{Float64}()
    @ccall libcb.cb_sparsecoeffmatmatrix_gram_ip2(self.data::Ptr{Cvoid}, ipval::Ref{Float64}, P.data::Ptr{Cvoid}, j::Cint)::Cint
    return ipval[]
end

@doc raw"""
    cb_primal_ip(self::CBSparseCoeffmatMatrix, ipvec::CBMatrix, primal::Union{<:CBPSCPrimal,Nothing}, ind::Union{<:CBIndexmatrix,Nothing} = nothing)

computes the inner products of (selected) columns (which represent block diagonal symmetric matrices) with the primal into the column vector ipvec(j)=ip(*primal,A.column((*ind)(j)}) (j=0,...,ind->dim()-1); if ind==NULL, use all columns
"""
cb_primal_ip(self::CBSparseCoeffmatMatrix, ipvec::CBMatrix, primal::Union{<:CBPSCPrimal,Nothing}, ind::Union{<:CBIndexmatrix,Nothing} = nothing) = @ccall libcb.cb_sparsecoeffmatmatrix_primal_ip(self.data::Ptr{Cvoid}, ipvec.data::Ptr{Cvoid}, (isnothing(primal) ? C_NULL : primal.data)::Ptr{Cvoid}, (isnothing(ind) ? C_NULL : ind.data)::Ptr{Cvoid})::Cint

@doc raw"""
    cb_primal_ip(self::CBSparseCoeffmatMatrix, primal::Union{<:CBPSCPrimal,Nothing}, j::Integer)

computes the inner product of the block diagonal symmetric matrix stored in colummn j with the primal into ipval=ip(*primal,A.column(j))
"""
function cb_primal_ip(self::CBSparseCoeffmatMatrix, primal::Union{<:CBPSCPrimal,Nothing}, j::Integer)
    value = Ref{Float64}()
    @ccall libcb.cb_sparsecoeffmatmatrix_primal_ip2(self.data::Ptr{Cvoid}, value::Ref{Float64}, (isnothing(primal) ? C_NULL : primal.data)::Ptr{Cvoid}, j::Cint)::Cint
    return value[]
end

@doc raw"""
    cb_project(self::CBSparseCoeffmatMatrix, S::CBSymmatrix, P::CBMatrix, j::Integer)

computes S=P^T*A.column(j)*P
"""
cb_project(self::CBSparseCoeffmatMatrix, S::CBSymmatrix, P::CBMatrix, j::Integer) = @ccall libcb.cb_sparsecoeffmatmatrix_project(self.data::Ptr{Cvoid}, S.data::Ptr{Cvoid}, P.data::Ptr{Cvoid}, j::Cint)::Cint

