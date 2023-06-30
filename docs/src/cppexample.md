# Example for the C++ interface

Here, we outline the [MAXCUT example](https://www-user.tu-chemnitz.de/~helmberg/ConicBundle/Manual/mctutorial.html), although
we don't implement all the custom triangle functionality (though it would certainly be possible).

```@repl
using ConicBundle
import MutableArithmetics as MA

# warning: we have to take care of garbage collection by ourselves!
garbage = [];
add_gc(x) = (push!(garbage, x); x);

STDIN = open("samplegraph.dat", "r")
# read the graph into (indi,indj,val) triples (replace val by 1.)
nnodes = parse(Int, readuntil(STDIN, " "))
medges = parse(Int, readuntil(STDIN, "\n"))

indi = add_gc(CBIndexmatrix(medges, 1, 0))
indj = add_gc(CBIndexmatrix(medges, 1, 0))
val = add_gc(CBMatrix(medges, 1, 1.))

for i in 0:medges-1
    head = parse(Int, readuntil(STDIN, " "))
    tail = parse(Int, readuntil(STDIN, " "))
    d = parse(Float64, readuntil(STDIN, "\n"))
    indi[i] = tail -1
    indj[i] = head -1
    # val[i] = d # edge weight d could be plugged in here
end

close(STDIN)

L = add_gc(CBSparsesym(nnodes, medges, indi, indj, val))
Ldiag = add_gc(CBMatrix(cb_diag(L)))
L = MA.operate!(-, L, add_gc(cb_sparseDiag(Ldiag, 1e-60)))
cb_init!(Ldiag, nnodes, 1, cb_sum(L) / nnodes)
L = MA.operate!(*, L, -1)
L = MA.operate!(+, L, add_gc(cb_sparseDiag(Ldiag, 1e-60)))
L = MA.operate!(/, L, 4.)
bestcut = add_gc(CBMatrix(nnodes, 1, 1.))
bestcut_val = cb_ip(bestcut, L * bestcut)

# form the Affine Matrix Function Oracle
Xdim = add_gc(CBIndexmatrix(1, 1, nnodes))
C = add_gc(CBSparseCoeffmatMatrix(Xdim, 1))

cb_set!(C, 0, 0, CBCMsymsparse(L)) # by default, ownership is managed by the coeffmat
opAt = add_gc(CBSparseCoeffmatMatrix(Xdim, nnodes))
for i in 0:nnodes-1
    cb_set!(opAt, 0, i, CBCMsingleton(nnodes, i, i, -1.)) # by default, the singletons are deleted by the coeffmat
end
mc = add_gc(CBPSCAffineFunction(C, opAt, CBGramSparsePSCPrimal(L))) # ownership of third argument passes to the function
cb_set_out!(mc, 0)

# initialize the solver and the problem
cbsolver = add_gc(CBMatrixCBSolver(1))
lb = add_gc(CBMatrix(nnodes, 1, cb_minus_infinity))
ub = add_gc(CBMatrix(nnodes, 1, cb_plus_infinity))
rhs = add_gc(CBMatrix(nnodes, 1, 1.))
cb_init_problem!(cbsolver, nnodes, lb, ub, nothing, rhs)
cb_add_function!(cbsolver, mc, nnodes, cbft_objective_function, nothing, true)

# call the solver
cb_set_active_bounds_fixing!(cbsolver, true)
cb_solve!(cbsolver) # here we could do all the rounding from the example
cb_print_termination_code(cbsolver)

cb_destroy!(@view(garbage[end:-1:1])...)
```