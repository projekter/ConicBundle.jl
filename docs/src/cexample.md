# Example for the C interface

This example is the direct counterpart to the miniature
[C example for a convex quadratic in two variables](https://www-user.tu-chemnitz.de/~helmberg/ConicBundle/Manual/group__cinterface.html).
```@repl
using ConicBundle

function eval_fun(function_key, x, relprec, max_new_subg, objective_threshold, subgval, subgradient, primal)
    # compute objective
    objective_value = .5*(5x[1]*x[1] + 2x[1]*x[2] + 4x[2]*x[2]) - 12x[1] - 10x[2] + 3;
    new_subg = 1
    subgval[1] = objective_value
    subgradient[1] = (5x[1] + x[2]) - 12;
    subgradient[2] = (x[1]+4*x[2]) - 10;
    return objective_value, new_subg
end

p = CBProblem{Int}() # the type parameter specified the type that our function_key will have - irrelevant here
cb_init_problem!(p, 2) # 2 variables, no bounds
cb_add_function!(p, 0, eval_fun) # the second parameter is our function_key; no extensions, no primal
cb_set_print_level!(p, 1)
cb_set_term_relprec!(p, 1e-8) # set relative precision
cb_solve!(p) # minimize the function up to termination
cb_print_termination_code(p)
x = cb_get_center(p) # retrieve the computed solution
print("x = ", x, ", objval = ", cb_get_objval(p))
# problem is automatically freed
```