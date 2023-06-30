import regex as re, glob

adapt = {"matrix.hxx": ["Matrix"], "indexmat.hxx": ["Indexmatrix"], "sparsmat.hxx": ["Sparsemat"],
         "symmat.hxx": ["Symmatrix"], "sparssym.hxx": ["Sparsesym"],
         "CMgramdense.hxx": ["CMgramdense<:Coeffmat"], "CMgramsparse.hxx": ["CMgramsparse<:Coeffmat"],
         "CMgramsparse_withoutdiag.hxx": ["CMgramsparse_withoutdiag<:Coeffmat"], "CMlowrankdd.hxx": ["CMlowrankdd<:Coeffmat"],
         "CMlowranksd.hxx": ["CMlowranksd<:Coeffmat"], "CMlowrankss.hxx": ["CMlowrankss<:Coeffmat"],
         "CMsingleton.hxx": ["CMsingleton<:Coeffmat"], "CMsymdense.hxx": ["CMsymdense<:Coeffmat"],
         "CMsymsparse.hxx": ["CMsymsparse<:Coeffmat"],
         "SparseCoeffmatMatrix.hxx": ["SparseCoeffmatMatrix"], "gb_rand.hxx": ["GB_rand"], "Coeffmat.hxx": ["CoeffmatInfo"],
         "MatrixCBSolver.hxx": ["PrimalMatrix<:PrimalData", "MatrixCBSolver"],
         "PSCPrimal.hxx": ["BlockPSCPrimal<:PSCPrimal<:PrimalData", "DensePSCPrimal<:PSCPrimal",
                           "GramSparsePSCPrimal<:PSCPrimal", "SparsePSCPrimal<:PSCPrimal"],
         "AFTModification.hxx": ["AFTModification<:OracleModification"],
         "GroundsetModification.hxx": ["GroundsetModification<:OracleModification"],
         "NNCBoxSupportModification.hxx": ["NNCBoxSupportModification<:OracleModification"],
         "PSCAffineModification.hxx": ["PSCAffineModification<:OracleModification"],
         "SOCSupportModification.hxx": ["SOCSupportModification<:OracleModification"],
         "AffineFunctionTransformation.hxx": ["AffineFunctionTransformation"],
         "MinorantPointer.hxx": ["MinorantPointer"], "MinorantUseData.hxx": ["MinorantUseData"],
         "CBSolver.hxx": ["Minorant"],
         "BoxOracle.hxx": ["BoxPrimalExtender<:PrimalExtender", "BoxOracle<:ModifiableOracleObject<:FunctionObject"],
         "PSCOracle.hxx": ["PSCPrimalExtender<:PrimalExtender", "PSCBundleParameters<:BundleParameters"],
         "SOCOracle.hxx": ["SOCPrimalExtender<:PrimalExtender", "SOCBundleParameters<:BundleParameters"],
         "CFunction.hxx": ["CFunctionMinorantExtender<:MinorantExtender",
                           "CFunction<:MatrixFunctionOracle<:ModifiableOracleObject"],
         "NNCBoxSupportFunction.hxx": ["NNCBoxSupportMinorantExtender<:MinorantExtender",
                                       "NNCBoxSupportFunction<:FunctionObject"],
         "PSCAffineFunction.hxx": ["PSCAffineMinorantExtender<:MinorantExtender",
                                   "PSCAffineFunction<:PSCOracle<:ModifiableOracleObject<:FunctionObject"],
         "SOCSupportFunction.hxx": ["SOCSupportMinorantExtender<:MinorantExtender",
                                    "SOCSupportFunction<:SOCOracle<:ModificableOracleObject<:FunctionObject"],
         "BoxModelParameters.hxx": ["BoxModelParameters<:BoxModelParametersObject<:BundleParameters"],
         "NNCModelParameters.hxx": ["NNCModelParameters<:NNCModelParametersObject<:BundleParameters"],
         "PSCModelParameters.hxx": ["PSCModelParameters<:PSCModelParametersObject<:BundleParameters"],
         "SOCModelParameters.hxx": ["SOCModelParameters<:SOCModelParametersObject<:BundleParameters"],
         "SumBundleParameters.hxx": ["SumBundleParameters<:SumBundleParametersObject<:BundleParameters"],
         "AFTData.hxx": ["AFTData<:BundleData<:VariableMetricBundleData"],
         "BoxData.hxx": ["BoxData<:BundleData<:VariableMetricBundleData"],
         "NNCData.hxx": ["NNCData<:BundleData<:VariableMetricBundleData"],
         "PSCData.hxx": ["PSCData<:BundleData<:VariableMetricBundleData"],
         "SOCData.hxx": ["SOCData<:BundleData<:VariableMetricBundleData"],
         "BundleHKWeight.hxx": ["BundleHKWeight<:BundleWeight"],
         "BundleRQBWeight.hxx": ["BundleRQBWeight<:BundleWeight"],
         "BundleDenseTrustRegionProx.hxx": ["BundleDenseTrustRegionProx<:BundleProxObject<:QPSolverProxObject"],
         "BundleDiagonalTrustRegionProx.hxx": ["BundleDiagonalTrustRegionProx<:BundleProxObject"],
         "BundleDLRTrustRegionProx.hxx": ["BundleDLRTrustRegionProx<:BundleProxObject"],
         "BundleIdProx.hxx": ["BundleIdProx<:BundleProxObject"],
         "BundleLowRankTrustRegionProx.hxx": ["BundleLowRankTrustRegionProx<:BundleProxObject"],
         "QPSolverParameters.hxx": ["QPSolverParameters<:QPSolverParametersObject"],
         "QPSolver.hxx": ["QPSolver<:QPSolverObject<:QPModelDataPointer"],
         "UQPSolver.hxx": ["UQPSolver<:QPSolverObject"],
         "LPGroundset.hxx": ["LPGroundset<:Groundset<:VariableMetricModel"],
         "UnconstrainedGroundset.hxx": ["UnconstrainedGroundset<:Groundset"],
         "AFTModel.hxx": ["AFTModel<:SumBlockModel<:BundleModel<:VariableMetricModel"],
         "BoxModel.hxx": ["BoxModel<:ConeModel<:SumBlockModel"],
         "NNCModel.hxx": ["NNCModel<:ConeModel"],
         "PSCModel.hxx": ["PSCModel<:ConeModel"],#, "QPPSCOracleData<:QPPSCOracleDataObject<:QPModelOracleDataObject"],
         "SOCModel.hxx": ["SOCModel<:ConeModel"],
         "SumModel.hxx": ["SumModel<:SumBlockModel"],
         "PSCVariableMetricSelection.hxx": ["PSCVariableMetricSelection<:VariableMetricSelection"],
         "VariableMetricSVDSelection.hxx": ["VariableMetricSVDSelection<:VariableMetricSelection"],
         "QPDirectKKTSolver.hxx": ["QPDirectKKTSolver<:QPKKTSolverObject"],
         "QPIterativeKKTHAeqSolver.hxx": ["QPIterativeKKTHAeqSolver<:QPIterativeKKTSolver<:QPKKTSolverObject"],
         "QPIterativeKKTHASolver.hxx": ["QPIterativeKKTHASolver<:QPIterativeKKTSolver<:QPKKTSolverObject"],
         "QPKKTSolverComparison.hxx": ["QPKKTSolverComparison<:QPKKTSOlverObject"],
         "SumBundleHandler.hxx": ["SumBundleHandler"],
         "QPConeModelBlock.hxx": ["QPConeModelBlock<:QPModelBlock<:QPModelBlockObject"],
         "QPSumModelBlock.hxx": ["QPSumModelBlock<:QPModelBlock"],
         "UQPConeModelBlock.hxx": ["UQPConeModelBlock<:UQPModelBlock<:UQPModelBlockObject"],
         "UQPSumModelBlock.hxx": ["UQPSumModelBlock<:UQPModelBlock"],
         "minres.hxx": ["MinRes<:IterativeSolverObject"],
         "pcg.hxx": ["PCG<:IterativeSolverObject"],
         "psqmr.hxx": ["Psqmr<:IterativeSolverObject"],
         "QPKKTSubspaceHPrecond.hxx": ["QPKKTSubspaceHPrecond<:QPKKTPrecondObject"],
         "SumBundle.hxx": ["SumBundle"],
         "BundleSolver.hxx": ["BundleSolver"],
         "BundleTerminator.hxx": ["BundleTerminator"],
         "clock.hxx": ["Clock", "Microseconds"]}
paramexpr = re.compile(r"(const\s+)?([A-Za-z_:\d]++)\s*([&*])?\s*([A-Za-z_\d]*)(?:\s*=\s*(.*))?")
nativetypes = {"void": "Nothing", "int": "Integer", "float": "Real", "double": "Real", "char": "UInt8", "long": "Integer",
               "short": "Integer", "Integer": "Integer", "Real": "Real", "bool": "Bool", "cb_subgextp": "Ptr{Cvoid}",
               "cb_functionp": "Ptr{Cvoid}"}
enums = {"FunctionTask": "CBFunctionTask", "ModelUpdate": "CBModelUpdate", "Mode": "CBMode", "SumBundle::Mode": "CBMode"}
nativeenums = nativetypes | enums
ignorefuncs = {"set_init", "get_init", "get_mtype", "operator[]", "clone_primal_data"}
operatornames = {"+": "plus", "-": "minus", "*": "times", "/": "divide", "%": "rem", "<": "less", ">": "greater",
                 "<=": "lessequal", ">=": "greaterequal", "==": "equal", "!=": "inequal"}
exclusions = {"cb_matrix_get_store", "cb_matrix_new_find2", "cb_matrix_new_find_number2",
              "cb_indexmatrix_new8", "cb_indexmatrix_new10", "cb_indexmatrix_get_store", "cb_indexmatrix_scale_rows", "cb_indexmatrix_scale_cols", "cb_indexmatrix_new_find2", "cb_indexmatrix_new_find_number2",
              "cb_symmatrix_get_store",
              "cb_sparsecoeffmatmatrix_get_colrep", "cb_sparsecoeffmatmatrix_block", "cb_sparsecoeffmatmatrix_column",
              "cb_matrixcbsolver_new2",
              "cb_sumbundle_eval_model", "cb_sumbundle_lb_model",
              "cb_sumbundlehandler_append_vars", "cb_sumbundlehandler_reassign_vars",
              "cb_pscdata_clear_model_except_bundlevecs",
              "cb_socdata_clear_model_except_bundlevecs",
              "cb_aftmodel_output_bundle_data",
              "cb_nncmodelparameters_new", "cb_pscmodelparameters_new", "cb_socmodelparameters_new",
              "cb_sumbundleparameters_new", "cb_bundlesolver_new"}
juliareserved = {"baremodule", "begin", "break", "catch", "const", "continue", "do", "else", "elseif", "end", "export",
                 "false", "finally", "for", "function", "global", "if", "import", "let", "local", "macro", "module", "quote",
                 "return", "struct", "true", "try", "using", "while"}
defaultreplacements = {"ObjectiveFunction": "cbft_objective_function", "SumBundle::inactive": "cbm_inactive"}
full_qualifier = {"CBProblem", "cb_clear!", "cb_set_default!", "cb_init_problem!", "cb_add_function!", "cb_set_lower_bound!",
                  "cb_set_upper_bound!", "cb_append_variables!", "cb_delete_variables!", "cb_reassign_variables!", "cb_solve!",
                  "cb_termination_code", "cb_print_termination_code", "cb_get_objval", "cb_get_center", "cb_get_center!",
                  "cb_get_sgnorm", "cb_get_subgradient", "cb_get_subgradient!", "cb_get_candidate_value", "cb_get_candidate",
                  "cb_get_candidate!", "cb_set_term_relprec!", "cb_set_new_center_point!", "cb_get_function_status",
                  "cb_get_approximate_slacks", "cb_get_approximate_slacks!", "cb_get_approximate_primal!",
                  "cb_get_center_primal!", "cb_get_candidate_primal!", "cb_set_max_modelsize!", "cb_set_max_bundlesize!",
                  "cb_set_max_new_subgradients!", "cb_get_bundle_parameters", "cb_reinit_function_model!",
                  "cb_get_last_weight", "cb_set_next_weight!", "cb_set_min_weight!", "cb_set_max_weight!",
                  "cb_set_variable_metric!", "cb_get_dim", "cb_get_n_functions", "cb_get_minus_infinity",
                  "cb_get_plus_infinity", "cb_clear_fail_counts!", "cb_set_eval_limit!", "cb_set_inner_update_limit!",
                  "cb_set_active_bounds_fixing!", "cb_get_fixed_active_bounds", "cb_get_fixed_active_bounds!",
                  "cb_set_print_level!", "cb_minus_infinity", "cb_plus_infinity"}

f_total = open("ConicBundle_cpp.jl", "w+")
f_total_c = open("../../ConicBundle/cppinterface/cb_cppinterface.cpp", "w+")
f_classes = open("cb_classes.jl", "w+")
f_doc = open("../../docs/src/cppref.md", "w+")
f_total.write("""function cb_destroy!(obj, rest...)
    cb_destroy!(obj)
    cb_destroy!(rest...)
end

""")#include("cb_classes.jl")\n')
f_total_c.write('#include "cb_cinterface.h"\n')
for fn in adapt:
    f_total_c.write(f'#include "{fn}"\n')
f_total_c.write("""using namespace CH_Matrix_Classes;
using namespace ConicBundle;
using namespace CH_Tools;
#define ModelUpdate BundleModel::ModelUpdate

extern "C" {

""")
f_classes.write(r"""mutable struct StdVector{T}
    start::Ref{T}
    stop::Ref{T}
    max::Ref{T}

    function StdVector(data::AbstractVector{T}) where {T}
        LinearAlgebra.chkstride1(data)
        stop = Ref(data, length(data))
        return new{T}(Ref(data), stop, stop) # no way to get the capacity
    end
end

Base.convert(::Type{<:StdVector}, vec::AbstractVector) = StdVector(vec)

""")
f_doc.write("""
```@meta
CurrentModule = ConicBundle
```
# Reference of C++ interface
The C++ interface was semi-automatically ported to Julia, for which a parser was written. It aims to bring almost the full
functionality of ConicBundle to Julia; however, the way of implementation is neither efficient nor always in a Julian style.
The documentation was extracted verbatim from the CPP code, hence no markdown profits will be available.
This basic interface can be improved in the future for example by exploiting memory alignment of structs in Julia and C++,
which would allow to inline lots of setters and getters instead of having to call external functions.

Note that the C++ interface requires you to manually keep track of memory. Every allocated object must be freed by an
appropriate call to [`cb_destroy!`](@ref). This is in contrast to the automatically managed [`CBProblem`](@ref) from the C
interface.

```@docs
""")
abstracts = set()
julia_signatures = set()
julia_functions = set()

for filename, classes in adapt.items():
    for curfile in glob.glob(f"..\\..\\ConicBundle\\**\\{filename}", recursive=True):
        with open(curfile, 'r') as file:
            datafull = file.read()
        allcomments = []
        def rep(match):
            allcomments.append(match.group(0))
            return f"$COMMENT{len(allcomments)-1}$"
        datafull = re.sub(r"//.*+|/\*(?:.|\n)*?\*/", rep, datafull)
        for curclass in classes:
            curclass_chain = curclass.split("<:")
            curclass = curclass_chain[0]
            if len(curclass_chain) > 1:
                for i in range(len(curclass_chain)-1, 0, -1):
                    if curclass_chain[i] not in abstracts:
                        if i < len(curclass_chain) -1:
                            f_classes.write(f"abstract type CB{curclass_chain[i]} <: CB{curclass_chain[i+1]} end\n\n")
                        else:
                            f_classes.write(f"abstract type CB{curclass_chain[i]} end\n\n")
                        abstracts.add(curclass_chain[i])
            print(curclass)
            nameprefix = f"cb_{curclass.lower()}_"
            usednames = {}
            data = re.search(r"""
\sclass\s+""" + curclass + r"""(?![a-zA-Z0-9_])[^{;]*\{
   (([^{}]*?(?:\{(?2)*\})?)*)
\}""", datafull, re.VERBOSE)
            if not data:
                print(" (no data)")
                continue
            data = data.group(1)
            data = re.sub(r"""
\s(?:class|struct)\s+[^{;]+\{
   ([^{}]*(?:\{(?1)*\})?)*+
\};
""", "", data, 0, re.VERBOSE)
            matches = re.finditer(r"(.*?)\s(public:|private:|protected:|$)+", data, re.DOTALL)
            data = ""
            last = ""
            for match in matches:
                md = match.group(1)
                if md:
                    md = md.strip()
                if last == "public:":
                    data += "\n" + md
                last = match.group(2)
            # while True:
            #     privpos = data.split("private:", 1)
            #     if len(privpos) == 1:
            #         break
            #     privpos[1] = data.split("public:", 1)
            #     if len(privpos[1]) == 1:
            #         data = privpos[0]
            #         break
            #     data = privpos[0] + "\n" + privpos[1][1]
            data = data.replace("CH_Matrix_Classes::", "")
            f_curclass_c = open("../../ConicBundle/cppinterface/" + nameprefix[:-1] + ".cpp", "w+")
            f_curclass_jl = open(nameprefix[:-1] + ".jl", "w+")
            f_total.write('include("' + nameprefix[:-1] + '.jl")\n')
            f_total_c.write('#include "' + nameprefix[:-1] + '.cpp"\n')
            juliatype = "CB" + curclass
            print(f"dll void {nameprefix}destroy({curclass}* self) {{\n  delete self;\n}}\n", file=f_curclass_c)
            julia_functions.add("cb_destroy!")
            print(f"""struct {juliatype}{'' if len(curclass_chain) == 1 else ' <: CB' + curclass_chain[1]}
    data::Ptr{{Cvoid}}
end
MA.mutability(::Type{{{juliatype}}}) = MA.IsMutable()

\"""
    cb_destroy!(obj::{juliatype})
\"""
cb_destroy!(obj::{juliatype}) = @ccall libcb.{nameprefix}destroy(obj.data::Ptr{{Cvoid}})::Cvoid
""", file=f_classes)
            # look for all possible functions
            for match in re.finditer(
                r"""
(?:\$COMMENT(\d+)\$)?
\s*\n\s*
(?:virtual\s+)?
(?:inline\s+)?
(friend\s+)?
(?:inline\s+)?
(?:
   ((?:const\s+)?[A-Za-z_&*:]+)
   \s+
   (operator(?:\(\)|\[\]|[+\-*/%=<>!]=?)|[A-Za-z_][A-Za-z_\d]*)
|
   (~?""" + curclass + """)
)
\s*
\(([^)]*)\)
\s*
(const)?
\s*
(?::[^{]+)?
(\{(?>[^{}]|(?8))*\})?""",
                data, re.VERBOSE | re.MULTILINE):
                comment, friend, returntype, name, conname, parameters, const, code = match.groups()
                origname = name
                if parameters:
                    parameters = re.sub(r"\$COMMENT\d+\$", "", parameters).strip()
                if parameters.find("<") != -1 or name in ignorefuncs or parameters.find("std::istream") != -1:
                    continue # no generic types
                if comment:
                    comment = allcomments[int(comment)]
                    if comment.startswith("/*"):
                        comment = comment[2:-2]
                    else:
                        comment = comment.strip(" /")
                friend = not not friend
                const = not not const
                if returntype == "operator":
                    # conversion operator
                    print("Skipping conversion operator ", origname)
                    continue

                if conname:
                    # this is a constructor/destructor
                    origname = conname
                    if conname[0] == '~':
                        continue
                    else:
                        name = "new"
                        returntype = conname
                name = name.lower()

                if returntype[-1] == "&":
                    returnstyle = "&"
                    returntype = returntype.rstrip("& ")
                elif returntype[-1] == "*":
                    returnstyle = "*"
                    returntype = returntype.rstrip("* ")
                else:
                    returnstyle = ""
                returnbool = returntype == "bool"
                if returntype in {"std::istream", "std::ostream", "FunObjModMap"}:
                    returntype = "void"
                    returnstyle = ""

                julia_parameters = []
                need_promote = False
                if name.startswith("operator"):
                    origname = name
                    if name == "operator()":
                        if const:
                            name = "get"
                            julia_name = "Base.getindex"
                        else:
                            name = "set"
                            julia_name = "Base.setindex" # ! is appended automatically
                            parameters += ", " + returntype + " value"
                            returntype = "void"
                            returnstyle = ""
                    elif name == "operator[]":
                        assert False # filtered
                    elif name == "operator=":
                        name = "assign"
                        julia_name = "Base.copy" # ! is appended automatically
                    elif name[-1] == "=" and name[-2] not in {"<", ">", "=", "!"}:
                        assert len(name) == 10
                        julia_name = "MA.operate" # ! is appended automatically
                        julia_parameters.append("::typeof(" + name[-2] + ")")
                        need_promote = True
                        name = operatornames[name[-2]]
                    else:
                        assert 9 <= len(name) <= 10
                        if name == "operator==":
                            julia_name = "Base.:(==)"
                        else:
                            julia_name = "Base.:" + name[8:]
                        name = "new_" + operatornames[name[8:]]
                    name = nameprefix + name
                else:
                    name = nameprefix + name
                    if returntype == "bool":
                        if returnstyle:
                            print("WARNING: bool reference return value found - assuming sizeof(bool) == 1")
                        else:
                            returntype = "int"
                    if conname:
                        julia_name = "CB" + origname
                    else:
                        julia_name = "cb_" + origname
                if not conname and not friend:
                    julia_parameters.append("self::" + juliatype)
                    julia_bodyparameters = ["self.data::Ptr{Cvoid}"]
                else:
                    julia_bodyparameters = []

                functionheader = returntype
                returntype = returntype.lstrip("const ")
                if returnstyle != "" or returntype not in nativeenums:
                    functionheader += "*"

                julia_body = f"@ccall libcb."
                julia_end = ")::"
                if returntype == "void" and returnstyle == "":
                    functionbody = "  "
                    functionend = ";\n}"
                    julia_end += "Cvoid"
                elif returntype in nativeenums:
                    functionbody = "  return "
                    functionend = ";\n}"
                    if returnstyle == "*":
                        julia_end += r"Ptr{"
                    elif returnstyle == "&":
                        functionbody += "&"
                        if returntype in nativetypes:
                            julia_end += r"Ptr{"
                    if returntype == "Integer":
                        julia_end += "Cint"
                    elif returntype == "Real":
                        julia_end += "Cdouble"
                    elif returntype == "bool":
                        julia_end += "Bool" # warning is printed above, as this might be wrong
                    elif returntype in {"cb_subgextp", "cb_functionp"}:
                        julia_end += "Ptr{Cvoid}"
                    elif returntype in nativetypes:
                        julia_end += f"C{returntype}"
                    elif returntype in enums:
                        julia_end += enums[returntype]
                        functionbody += "(int)"
                        if returnstyle != "":
                            functionheader = "*int"
                        else:
                            functionheader = "int"
                    if returnstyle == "*" or (returnstyle == "&" and returntype in nativetypes):
                        julia_end += "}"
                    elif returnbool:
                        julia_body = "Bool(" + julia_body
                        julia_end += ")"
                elif returnstyle == "*":
                    functionbody = "  return "
                    functionend = ";\n}"
                    julia_body = f"CB{returntype}({julia_body}"
                    julia_end += "Ptr{Cvoid})"
                elif returnstyle == "&":
                    functionbody = "  return &"
                    functionend = ";\n}"
                    julia_body = "(" + julia_body
                    julia_end += r"Ptr{Cvoid}; return self)"
                else:
                    assert returnstyle == ""
                    functionbody = "  return new " + returntype + "("
                    functionend = ");\n}"
                    if not name.startswith(nameprefix + "new"):
                        name = nameprefix + "new_" + name[len(nameprefix):]
                    julia_body = "CB" + returntype + "(" + julia_body
                    julia_end += r"Ptr{Cvoid})"

                if name in usednames:
                    usednames[name] += 1
                    name = name + str(usednames[name])
                else:
                    usednames[name] = 1
                if name in exclusions:
                    # we still do the increment in order to be able to uniquely identify the functions
                    continue
                functionheader += " " + name + "("
                julia_body += name + "("
                julia_return = []

                functionparams = []
                if friend:
                    if not origname.startswith("operator"):
                        functionbody += origname + "("
                        functionend = ")" + functionend
                elif not conname:
                    if const:
                        functionparam = "const "
                    else:
                        functionparam = ""
                        julia_name += "!"
                    functionparam += curclass + "* self"
                    if origname == "operator()":
                        functionbody += "(*self)("
                        if const:
                            functionend = ")" + functionend
                        else:
                            functionend = ") = value" + functionend
                    elif origname == "operator=":
                        functionbody += "(*self = "
                        functionend = ")" + functionend
                    elif origname == "operator-" and parameters == "":
                        functionbody += "-(*self"
                        functionend = ")" + functionend
                    elif origname.startswith("operator"):
                        if origname[-1] != "=":
                            if not returnbool:
                                print(f"WARNING: {origname} ignored, not a friend but also not an assignment")
                                continue
                                # should be corrected in the code, these instances cannot even be resolved by the cpp compiler!
                            assert const
                            functionbody += "self->" + origname + "("
                            functionend = ")" + functionend
                        else:
                            functionbody += "(*self " + origname[8:] + " "
                            functionend = ")" + functionend
                    else:
                        functionbody += "self->" + origname + "("
                        functionend = ")" + functionend
                    functionparams.append(functionparam)

                preserve = []
                if not parameters or parameters == "void":
                    parameters = []
                else:
                    parameters = parameters.split(",")
                    functionbodyparams = []
                    julia_skipparam = 0
                    i = 0
                    while i < len(parameters):
                        try:
                            param = parameters[i]
                            param = param.strip()
                            paramm = re.fullmatch(paramexpr, param)
                            if not paramm:
                                print("ERROR: unsupported pattern in ", name, ": ", param)
                                i = -2
                                break
                            paramconst, paramtype, paramref, paramname, paramdefault = paramm.groups()
                            if paramname in juliareserved:
                                paramname += "_"
                            if not paramname:
                                paramname = f"param{i}"
                            boolparam = False
                            if paramtype == "Range":
                                assert paramref != "*"
                                functionparams.append(f"Integer {paramname}_from, Integer {paramname}_to, Integer {paramname}_step = 1")
                                functionbodyparams.append(f"Range({paramname}_from, {paramname}_to, {paramname}_step)")
                                julia_parameters.append(f"{paramname}::AbstractRange{{<:Integer}}")
                                julia_bodyparameters.append(f"first({paramname})::Cint, last({paramname})::Cint, step({paramname})::Cint")
                                continue
                            elif paramtype == "Realrange":
                                assert paramref != "*"
                                functionparams.append(f"Real {paramname}_from, Real {paramname}_to, Real {paramname}_step = 1, Real {paramname}_tol = 1e-8")
                                functionbodyparams.append(f"Realrange({paramname}_from, {paramname}_to, {paramname}_step, {paramname}_tol)")
                                julia_parameters.append(f"{paramname}::AbstractRange{{<:Real}}, {paramname}_tol::Real = 1e-8")
                                julia_bodyparameters.append(f"first({paramname})::Cdouble, last({paramname})::Cdouble, step({paramname})::Cdouble, {paramname}_tol::Cdouble")
                                continue
                            elif paramtype == "std::ostream":
                                if paramref == "&":
                                    functionbodyparams.append("std::cout")
                                elif paramref == "*":
                                    functionbodyparams.append("&std::cout")
                                else:
                                    assert False
                                continue
                            elif paramtype == "FunObjModMap":
                                if paramref != "*":
                                    print("ERROR: Untranslatable parameter in ", name, ": ", param)
                                    i = -2
                                    break
                                functionbodyparams.append("(FunObjModMap*)0")
                                continue
                            elif paramtype == "CBout":
                                assert paramref == "*"
                                functionbodyparams.append("0")
                                continue
                            paramconst = paramconst == "const "
                            if paramconst:
                                functionparam = "const "
                            else:
                                functionparam = ""
                            julia_bodyparameter = paramname + "::"
                            if paramtype == "bool":
                                boolparam = True
                                paramtype = "int"
                                julia_parameter = f"{paramname}::Bool"
                                if paramref:
                                    julia_bodyparameter += r"Ref{Cint}"
                                    if paramref == "*":
                                        functionbody = f"  *{paramname} = 0;\n{functionbody}"
                                    else:
                                        functionbody = f"  {paramname} = 0;\n{functionbody}"
                                else:
                                    julia_bodyparameter += "Cint"
                            elif paramtype in nativeenums:
                                if paramref == "*":
                                    # in principle StridedVector, but it is enough if the stride interface is supported (may be
                                    # views)
                                    if paramtype == "Integer":
                                        julia_parameter = f"{paramname}::Union{{<:AbstractVector{{Cint}},Nothing}}"
                                    elif paramtype == "Real":
                                        julia_parameter = f"{paramname}::Union{{<:AbstractVector{{Cdouble}},Nothing}}"
                                    else:
                                        julia_parameter = f"{paramname}::Union{{<:AbstractVector{{{nativeenums[paramtype]}}},Nothing}}"
                                    preserve.append(paramname)
                                elif paramref == "&":
                                    # these should be return values
                                    if paramtype == "Integer":
                                        jt = "Int"
                                    elif paramtype == "Real":
                                        jt = "Float64"
                                    else:
                                        jt = nativeenums[paramtype]
                                    julia_body = f"{paramname} = Ref{{{jt}}}()\n{julia_body}"
                                    julia_bodyparameters.append(f"{paramname}::Ref{{{jt}}}")
                                    julia_return.append(f"{paramname}[]")
                                    julia_skipparam = 1
                                    if julia_name == "cb_dim":
                                        julia_name = "cb_dim2"
                                else:
                                    julia_parameter = f"{paramname}::{nativeenums[paramtype]}"
                                if paramref:
                                    julia_bodyparameter += r"Ptr{"
                                if paramtype == "Integer" or paramtype in enums:
                                    julia_bodyparameter += "Cint"
                                elif paramtype == "Real":
                                    julia_bodyparameter += "Cdouble"
                                elif paramtype in {"cb_subgextp", "cb_functionp"}:
                                    julia_bodyparameter += "Ptr{Cvoid}"
                                else:
                                    julia_bodyparameter += f"C{paramtype.rsplit('::', 1)[-1]}"
                                if paramref:
                                    julia_bodyparameter += "}"
                                if paramref == "*":
                                    if i < len(parameters) -1 and parameters[i+1].strip().startswith(("Integer inc", "int inc")):
                                        julia_skipparam = 2
                                        julia_bodyparameter += f", stride({paramname}, 1)::Cint"
                                    elif julia_body[0] == "(":
                                        julia_body = f"(LinearAlgebra.chkstride1({paramname}); {julia_body[1:]}"
                                    else:
                                        julia_body = f"(LinearAlgebra.chkstride1({paramname}); {julia_body}"
                                        julia_end += ")"
                            else:
                                if paramref == "&" or not paramref:
                                    julia_parameter = f"{paramname}::CB{paramtype.rsplit('::', 1)[-1]}"
                                    julia_bodyparameter = f"{paramname}.data::Ptr{{Cvoid}}"
                                else:
                                    assert paramref == "*"
                                    julia_parameter = f"{paramname}::Union{{<:CB{paramtype.rsplit('::', 1)[-1]},Nothing}}"
                                    julia_bodyparameter = f"(isnothing({paramname}) ? C_NULL : {paramname}.data)::Ptr{{Cvoid}}"
                            if paramtype in enums:
                                functionparam += "int"
                            else:
                                functionparam += paramtype
                            if not paramref:
                                if boolparam:
                                    # if there are overloaded variants, we might not be able to rely on automatic conversion
                                    functionbodyparams.append("(bool)" + paramname)
                                elif paramtype in enums:
                                    functionbodyparams.append(f"({paramtype}){paramname}")
                                elif paramtype in nativetypes:
                                    functionbodyparams.append(paramname)
                                else:
                                    functionbodyparams.append("*" + paramname)
                                    paramref = "*"
                            elif paramref == "*":
                                if boolparam:
                                    functionbodyparams.append("(bool*)" + paramname)
                                elif paramtype in enums:
                                    functionbodyparams.append(f"({paramtype}*){paramname}")
                                else:
                                    functionbodyparams.append(paramname)
                            else:
                                assert paramref == "&"
                                if boolparam:
                                    functionbodyparams.append("*(bool*)" + paramname)
                                elif paramtype in enums:
                                    functionbodyparams.append(f"*({paramtype}*){paramname}")
                                else:
                                    functionbodyparams.append("*" + paramname)
                            if paramref:
                                functionparam += "* " + paramname
                            else:
                                functionparam += " " + paramname
                            if paramdefault:
                                functionparam += " = " + paramdefault.replace("false", "0").replace("true", "1")
                                if paramdefault == "0" and paramref == "*":
                                    julia_parameter += " = nothing"
                                else:
                                    for dr, rep in defaultreplacements.items():
                                        paramdefault = paramdefault.replace(dr, rep)
                                    julia_parameter += " = " + paramdefault
                            if julia_skipparam != 1:
                                julia_parameters.append(julia_parameter)
                                julia_bodyparameters.append(julia_bodyparameter)
                            if julia_skipparam > 0:
                                julia_skipparam -= 1
                            functionparams.append(functionparam)
                        finally:
                            i += 1
                    if i < 0:
                        continue
                    if origname == "operator()" and not const:
                        functionbodyparams.pop()
                    if friend and origname.startswith("operator"):
                        functionbody += (" " + origname[8:] + " ").join(functionbodyparams)
                    else:
                        functionbody += ", ".join(functionbodyparams)
                functionheader += ", ".join(functionparams) + ") {\n"
                print(f"dll {functionheader}{functionbody}{functionend}\n", file=f_curclass_c)
                if julia_name == "Base.setindex!":
                    julia_parameters.insert(1, julia_parameters.pop())
                elif julia_name == "Base.:!=":
                    continue # should not be overloaded
                julia_signature = julia_name + "(" + ", ".join(julia_parameters) + ")"
                julia_typesignature = re.sub(r"[a-zA-Z_][a-zA-Z\d_]*::([a-zA-Z_][a-zA-Z\d_]*)(?: = [^),]+)?", r"_::\1", julia_signature)
                if julia_typesignature in julia_signatures:
                    print(f"skipping {julia_signature}, already present in indistinguishable form")
                    continue
                else:
                    julia_signatures.add(julia_typesignature)
                if julia_return:
                    julia_header = f"function {julia_signature}\n    "
                    julia_body = (julia_body.replace("\n", "\n    ") + ", ".join(julia_bodyparameters) + julia_end +
                                  "\n    return " + ", ".join(julia_return) + "\nend")
                    if preserve:
                        julia_body = (f"    GC.@preserve {' '.join(preserve)} begin\n" + julia_body.replace('\n', '\n    ') +
                                      "\n    end")
                else:
                    julia_header = f"{julia_signature} = "
                    julia_body = julia_body + ", ".join(julia_bodyparameters) + julia_end
                    if preserve:
                        julia_body = f"GC.@preserve {' '.join(preserve)} begin\n    {julia_body}\nend"
                if comment:
                    julia_header = ('@doc raw"""\n    ' + julia_name + "(" + ", ".join(julia_parameters) + ")\n\n" +
                        comment + '\n"""\n' + julia_header)
                if julia_name.find(".") == -1:
                    if julia_name not in full_qualifier:
                        julia_functions.add(julia_name)
                    else:
                        for i in range(len(julia_parameters)):
                            julia_parameters[i] = "::" + julia_parameters[i].split("::", 1)[1]
                            julia_parameters[i] = julia_parameters[i].split(" = ", 1)[0]
                        julia_functions.add(f"{julia_name}({', '.join(julia_parameters)})")
                print(f"{julia_header}{julia_body}\n", file=f_curclass_jl)
                if need_promote:
                    # first parameter is the typeof
                    for i in range(1, len(julia_parameters)):
                        julia_parameters[i] = "::Type{<:" + julia_parameters[i].split("::", 1)[1].split(" = ", 1)[0] + "}"
                    print(f"MA.promote_operation({', '.join(julia_parameters)}) = {juliatype}\n",
                          file=f_curclass_jl)
            f_curclass_c.close()
            f_curclass_jl.close()

julia_functions |= {"CBModelUpdate", "cbmu_new_subgradient", "cbmu_descent_step", "cbmu_null_step",
                    "CBMode", "cbm_root", "cbm_child", "cbm_inactive", "cbm_unavailable"}
f_classes.write(r'''const CBCoeffmatVector = StdVector{<:CBCoeffmat}

const CBCoeffmatPointer = Ref{<:CBCoeffmat}

const CBMinorantBundle = StdVector{CBMinorantPointer}

const CBDVector = StdVector{Cdouble}

const CBIVector = StdVector{Cint}

"""
    enum CBModelUpdate

- `cbmu_new_subgradient`
- `cbmu_descent_step`
- `cbmu_null_step`
"""
@enum CBModelUpdate::Cint cbmu_new_subgradient=0 cbmu_descent_step=1 cbmu_null_step=2

"""
    enum CBMode

- `cbm_root`
- `cbm_child`
- `cbm_inactive`
- `cbm_unavailable`
"""
@enum CBMode::Cint cbm_root=0 cbm_child=1 cbm_inactive=2 cbm_unavailable=3

const CBQPSumModelDataObject = Union{CBQPSumModelBlock,CBUQPSumModelBlock}

const CBQPConeModelDataObject = Union{CBQPConeModelBlock,CBUQPSumModelBlock}

const CBQPModelDataObject = Union{<:CBQPSumModelDataObject,<:CBQPConeModelDataObject}

const CBIterativeSystemObject = Union{<:CBQPIterativeKKTSolver}

abstract type CBQPModelOracleDataObject end

abstract type CBQPPSCOracleDataObject <: CBQPModelOracleDataObject end''')
f_classes.close()
julia_functions_notype = sorted({jf.split("(", 1)[0] for jf in julia_functions})
julia_functions = sorted(julia_functions)
f_total.write("export " + ",\n    ".join(julia_functions_notype))
f_doc.write("\n".join(julia_functions) + "\n```")
f_doc.close()
f_total.close()
f_total_c.write("\n}")
f_total_c.close()