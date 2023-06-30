mutable struct StdVector{T}
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

struct CBMatrix
    data::Ptr{Cvoid}
end
MA.mutability(::Type{CBMatrix}) = MA.IsMutable()

"""
    cb_destroy!(obj::CBMatrix)
"""
cb_destroy!(obj::CBMatrix) = @ccall libcb.cb_matrix_destroy(obj.data::Ptr{Cvoid})::Cvoid

struct CBIndexmatrix
    data::Ptr{Cvoid}
end
MA.mutability(::Type{CBIndexmatrix}) = MA.IsMutable()

"""
    cb_destroy!(obj::CBIndexmatrix)
"""
cb_destroy!(obj::CBIndexmatrix) = @ccall libcb.cb_indexmatrix_destroy(obj.data::Ptr{Cvoid})::Cvoid

struct CBSparsemat
    data::Ptr{Cvoid}
end
MA.mutability(::Type{CBSparsemat}) = MA.IsMutable()

"""
    cb_destroy!(obj::CBSparsemat)
"""
cb_destroy!(obj::CBSparsemat) = @ccall libcb.cb_sparsemat_destroy(obj.data::Ptr{Cvoid})::Cvoid

struct CBSymmatrix
    data::Ptr{Cvoid}
end
MA.mutability(::Type{CBSymmatrix}) = MA.IsMutable()

"""
    cb_destroy!(obj::CBSymmatrix)
"""
cb_destroy!(obj::CBSymmatrix) = @ccall libcb.cb_symmatrix_destroy(obj.data::Ptr{Cvoid})::Cvoid

struct CBSparsesym
    data::Ptr{Cvoid}
end
MA.mutability(::Type{CBSparsesym}) = MA.IsMutable()

"""
    cb_destroy!(obj::CBSparsesym)
"""
cb_destroy!(obj::CBSparsesym) = @ccall libcb.cb_sparsesym_destroy(obj.data::Ptr{Cvoid})::Cvoid

abstract type CBCoeffmat end

struct CBCMgramdense <: CBCoeffmat
    data::Ptr{Cvoid}
end
MA.mutability(::Type{CBCMgramdense}) = MA.IsMutable()

"""
    cb_destroy!(obj::CBCMgramdense)
"""
cb_destroy!(obj::CBCMgramdense) = @ccall libcb.cb_cmgramdense_destroy(obj.data::Ptr{Cvoid})::Cvoid

struct CBCMgramsparse <: CBCoeffmat
    data::Ptr{Cvoid}
end
MA.mutability(::Type{CBCMgramsparse}) = MA.IsMutable()

"""
    cb_destroy!(obj::CBCMgramsparse)
"""
cb_destroy!(obj::CBCMgramsparse) = @ccall libcb.cb_cmgramsparse_destroy(obj.data::Ptr{Cvoid})::Cvoid

struct CBCMgramsparse_withoutdiag <: CBCoeffmat
    data::Ptr{Cvoid}
end
MA.mutability(::Type{CBCMgramsparse_withoutdiag}) = MA.IsMutable()

"""
    cb_destroy!(obj::CBCMgramsparse_withoutdiag)
"""
cb_destroy!(obj::CBCMgramsparse_withoutdiag) = @ccall libcb.cb_cmgramsparse_withoutdiag_destroy(obj.data::Ptr{Cvoid})::Cvoid

struct CBCMlowrankdd <: CBCoeffmat
    data::Ptr{Cvoid}
end
MA.mutability(::Type{CBCMlowrankdd}) = MA.IsMutable()

"""
    cb_destroy!(obj::CBCMlowrankdd)
"""
cb_destroy!(obj::CBCMlowrankdd) = @ccall libcb.cb_cmlowrankdd_destroy(obj.data::Ptr{Cvoid})::Cvoid

struct CBCMlowranksd <: CBCoeffmat
    data::Ptr{Cvoid}
end
MA.mutability(::Type{CBCMlowranksd}) = MA.IsMutable()

"""
    cb_destroy!(obj::CBCMlowranksd)
"""
cb_destroy!(obj::CBCMlowranksd) = @ccall libcb.cb_cmlowranksd_destroy(obj.data::Ptr{Cvoid})::Cvoid

struct CBCMlowrankss <: CBCoeffmat
    data::Ptr{Cvoid}
end
MA.mutability(::Type{CBCMlowrankss}) = MA.IsMutable()

"""
    cb_destroy!(obj::CBCMlowrankss)
"""
cb_destroy!(obj::CBCMlowrankss) = @ccall libcb.cb_cmlowrankss_destroy(obj.data::Ptr{Cvoid})::Cvoid

struct CBCMsingleton <: CBCoeffmat
    data::Ptr{Cvoid}
end
MA.mutability(::Type{CBCMsingleton}) = MA.IsMutable()

"""
    cb_destroy!(obj::CBCMsingleton)
"""
cb_destroy!(obj::CBCMsingleton) = @ccall libcb.cb_cmsingleton_destroy(obj.data::Ptr{Cvoid})::Cvoid

struct CBCMsymdense <: CBCoeffmat
    data::Ptr{Cvoid}
end
MA.mutability(::Type{CBCMsymdense}) = MA.IsMutable()

"""
    cb_destroy!(obj::CBCMsymdense)
"""
cb_destroy!(obj::CBCMsymdense) = @ccall libcb.cb_cmsymdense_destroy(obj.data::Ptr{Cvoid})::Cvoid

struct CBCMsymsparse <: CBCoeffmat
    data::Ptr{Cvoid}
end
MA.mutability(::Type{CBCMsymsparse}) = MA.IsMutable()

"""
    cb_destroy!(obj::CBCMsymsparse)
"""
cb_destroy!(obj::CBCMsymsparse) = @ccall libcb.cb_cmsymsparse_destroy(obj.data::Ptr{Cvoid})::Cvoid

struct CBSparseCoeffmatMatrix
    data::Ptr{Cvoid}
end
MA.mutability(::Type{CBSparseCoeffmatMatrix}) = MA.IsMutable()

"""
    cb_destroy!(obj::CBSparseCoeffmatMatrix)
"""
cb_destroy!(obj::CBSparseCoeffmatMatrix) = @ccall libcb.cb_sparsecoeffmatmatrix_destroy(obj.data::Ptr{Cvoid})::Cvoid

struct CBGB_rand
    data::Ptr{Cvoid}
end
MA.mutability(::Type{CBGB_rand}) = MA.IsMutable()

"""
    cb_destroy!(obj::CBGB_rand)
"""
cb_destroy!(obj::CBGB_rand) = @ccall libcb.cb_gb_rand_destroy(obj.data::Ptr{Cvoid})::Cvoid

struct CBCoeffmatInfo
    data::Ptr{Cvoid}
end
MA.mutability(::Type{CBCoeffmatInfo}) = MA.IsMutable()

"""
    cb_destroy!(obj::CBCoeffmatInfo)
"""
cb_destroy!(obj::CBCoeffmatInfo) = @ccall libcb.cb_coeffmatinfo_destroy(obj.data::Ptr{Cvoid})::Cvoid

abstract type CBPrimalData end

struct CBPrimalMatrix <: CBPrimalData
    data::Ptr{Cvoid}
end
MA.mutability(::Type{CBPrimalMatrix}) = MA.IsMutable()

"""
    cb_destroy!(obj::CBPrimalMatrix)
"""
cb_destroy!(obj::CBPrimalMatrix) = @ccall libcb.cb_primalmatrix_destroy(obj.data::Ptr{Cvoid})::Cvoid

struct CBMatrixCBSolver
    data::Ptr{Cvoid}
end
MA.mutability(::Type{CBMatrixCBSolver}) = MA.IsMutable()

"""
    cb_destroy!(obj::CBMatrixCBSolver)
"""
cb_destroy!(obj::CBMatrixCBSolver) = @ccall libcb.cb_matrixcbsolver_destroy(obj.data::Ptr{Cvoid})::Cvoid

abstract type CBPSCPrimal <: CBPrimalData end

struct CBBlockPSCPrimal <: CBPSCPrimal
    data::Ptr{Cvoid}
end
MA.mutability(::Type{CBBlockPSCPrimal}) = MA.IsMutable()

"""
    cb_destroy!(obj::CBBlockPSCPrimal)
"""
cb_destroy!(obj::CBBlockPSCPrimal) = @ccall libcb.cb_blockpscprimal_destroy(obj.data::Ptr{Cvoid})::Cvoid

struct CBDensePSCPrimal <: CBPSCPrimal
    data::Ptr{Cvoid}
end
MA.mutability(::Type{CBDensePSCPrimal}) = MA.IsMutable()

"""
    cb_destroy!(obj::CBDensePSCPrimal)
"""
cb_destroy!(obj::CBDensePSCPrimal) = @ccall libcb.cb_densepscprimal_destroy(obj.data::Ptr{Cvoid})::Cvoid

struct CBGramSparsePSCPrimal <: CBPSCPrimal
    data::Ptr{Cvoid}
end
MA.mutability(::Type{CBGramSparsePSCPrimal}) = MA.IsMutable()

"""
    cb_destroy!(obj::CBGramSparsePSCPrimal)
"""
cb_destroy!(obj::CBGramSparsePSCPrimal) = @ccall libcb.cb_gramsparsepscprimal_destroy(obj.data::Ptr{Cvoid})::Cvoid

struct CBSparsePSCPrimal <: CBPSCPrimal
    data::Ptr{Cvoid}
end
MA.mutability(::Type{CBSparsePSCPrimal}) = MA.IsMutable()

"""
    cb_destroy!(obj::CBSparsePSCPrimal)
"""
cb_destroy!(obj::CBSparsePSCPrimal) = @ccall libcb.cb_sparsepscprimal_destroy(obj.data::Ptr{Cvoid})::Cvoid

abstract type CBOracleModification end

struct CBAFTModification <: CBOracleModification
    data::Ptr{Cvoid}
end
MA.mutability(::Type{CBAFTModification}) = MA.IsMutable()

"""
    cb_destroy!(obj::CBAFTModification)
"""
cb_destroy!(obj::CBAFTModification) = @ccall libcb.cb_aftmodification_destroy(obj.data::Ptr{Cvoid})::Cvoid

struct CBGroundsetModification <: CBOracleModification
    data::Ptr{Cvoid}
end
MA.mutability(::Type{CBGroundsetModification}) = MA.IsMutable()

"""
    cb_destroy!(obj::CBGroundsetModification)
"""
cb_destroy!(obj::CBGroundsetModification) = @ccall libcb.cb_groundsetmodification_destroy(obj.data::Ptr{Cvoid})::Cvoid

struct CBNNCBoxSupportModification <: CBOracleModification
    data::Ptr{Cvoid}
end
MA.mutability(::Type{CBNNCBoxSupportModification}) = MA.IsMutable()

"""
    cb_destroy!(obj::CBNNCBoxSupportModification)
"""
cb_destroy!(obj::CBNNCBoxSupportModification) = @ccall libcb.cb_nncboxsupportmodification_destroy(obj.data::Ptr{Cvoid})::Cvoid

struct CBPSCAffineModification <: CBOracleModification
    data::Ptr{Cvoid}
end
MA.mutability(::Type{CBPSCAffineModification}) = MA.IsMutable()

"""
    cb_destroy!(obj::CBPSCAffineModification)
"""
cb_destroy!(obj::CBPSCAffineModification) = @ccall libcb.cb_pscaffinemodification_destroy(obj.data::Ptr{Cvoid})::Cvoid

struct CBSOCSupportModification <: CBOracleModification
    data::Ptr{Cvoid}
end
MA.mutability(::Type{CBSOCSupportModification}) = MA.IsMutable()

"""
    cb_destroy!(obj::CBSOCSupportModification)
"""
cb_destroy!(obj::CBSOCSupportModification) = @ccall libcb.cb_socsupportmodification_destroy(obj.data::Ptr{Cvoid})::Cvoid

struct CBAffineFunctionTransformation
    data::Ptr{Cvoid}
end
MA.mutability(::Type{CBAffineFunctionTransformation}) = MA.IsMutable()

"""
    cb_destroy!(obj::CBAffineFunctionTransformation)
"""
cb_destroy!(obj::CBAffineFunctionTransformation) = @ccall libcb.cb_affinefunctiontransformation_destroy(obj.data::Ptr{Cvoid})::Cvoid

struct CBMinorantPointer
    data::Ptr{Cvoid}
end
MA.mutability(::Type{CBMinorantPointer}) = MA.IsMutable()

"""
    cb_destroy!(obj::CBMinorantPointer)
"""
cb_destroy!(obj::CBMinorantPointer) = @ccall libcb.cb_minorantpointer_destroy(obj.data::Ptr{Cvoid})::Cvoid

struct CBMinorantUseData
    data::Ptr{Cvoid}
end
MA.mutability(::Type{CBMinorantUseData}) = MA.IsMutable()

"""
    cb_destroy!(obj::CBMinorantUseData)
"""
cb_destroy!(obj::CBMinorantUseData) = @ccall libcb.cb_minorantusedata_destroy(obj.data::Ptr{Cvoid})::Cvoid

struct CBMinorant
    data::Ptr{Cvoid}
end
MA.mutability(::Type{CBMinorant}) = MA.IsMutable()

"""
    cb_destroy!(obj::CBMinorant)
"""
cb_destroy!(obj::CBMinorant) = @ccall libcb.cb_minorant_destroy(obj.data::Ptr{Cvoid})::Cvoid

abstract type CBPrimalExtender end

struct CBBoxPrimalExtender <: CBPrimalExtender
    data::Ptr{Cvoid}
end
MA.mutability(::Type{CBBoxPrimalExtender}) = MA.IsMutable()

"""
    cb_destroy!(obj::CBBoxPrimalExtender)
"""
cb_destroy!(obj::CBBoxPrimalExtender) = @ccall libcb.cb_boxprimalextender_destroy(obj.data::Ptr{Cvoid})::Cvoid

abstract type CBFunctionObject end

abstract type CBModifiableOracleObject <: CBFunctionObject end

struct CBBoxOracle <: CBModifiableOracleObject
    data::Ptr{Cvoid}
end
MA.mutability(::Type{CBBoxOracle}) = MA.IsMutable()

"""
    cb_destroy!(obj::CBBoxOracle)
"""
cb_destroy!(obj::CBBoxOracle) = @ccall libcb.cb_boxoracle_destroy(obj.data::Ptr{Cvoid})::Cvoid

struct CBPSCPrimalExtender <: CBPrimalExtender
    data::Ptr{Cvoid}
end
MA.mutability(::Type{CBPSCPrimalExtender}) = MA.IsMutable()

"""
    cb_destroy!(obj::CBPSCPrimalExtender)
"""
cb_destroy!(obj::CBPSCPrimalExtender) = @ccall libcb.cb_pscprimalextender_destroy(obj.data::Ptr{Cvoid})::Cvoid

abstract type CBBundleParameters end

struct CBPSCBundleParameters <: CBBundleParameters
    data::Ptr{Cvoid}
end
MA.mutability(::Type{CBPSCBundleParameters}) = MA.IsMutable()

"""
    cb_destroy!(obj::CBPSCBundleParameters)
"""
cb_destroy!(obj::CBPSCBundleParameters) = @ccall libcb.cb_pscbundleparameters_destroy(obj.data::Ptr{Cvoid})::Cvoid

struct CBSOCPrimalExtender <: CBPrimalExtender
    data::Ptr{Cvoid}
end
MA.mutability(::Type{CBSOCPrimalExtender}) = MA.IsMutable()

"""
    cb_destroy!(obj::CBSOCPrimalExtender)
"""
cb_destroy!(obj::CBSOCPrimalExtender) = @ccall libcb.cb_socprimalextender_destroy(obj.data::Ptr{Cvoid})::Cvoid

struct CBSOCBundleParameters <: CBBundleParameters
    data::Ptr{Cvoid}
end
MA.mutability(::Type{CBSOCBundleParameters}) = MA.IsMutable()

"""
    cb_destroy!(obj::CBSOCBundleParameters)
"""
cb_destroy!(obj::CBSOCBundleParameters) = @ccall libcb.cb_socbundleparameters_destroy(obj.data::Ptr{Cvoid})::Cvoid

abstract type CBMinorantExtender end

struct CBCFunctionMinorantExtender <: CBMinorantExtender
    data::Ptr{Cvoid}
end
MA.mutability(::Type{CBCFunctionMinorantExtender}) = MA.IsMutable()

"""
    cb_destroy!(obj::CBCFunctionMinorantExtender)
"""
cb_destroy!(obj::CBCFunctionMinorantExtender) = @ccall libcb.cb_cfunctionminorantextender_destroy(obj.data::Ptr{Cvoid})::Cvoid

abstract type CBMatrixFunctionOracle <: CBModifiableOracleObject end

struct CBCFunction <: CBMatrixFunctionOracle
    data::Ptr{Cvoid}
end
MA.mutability(::Type{CBCFunction}) = MA.IsMutable()

"""
    cb_destroy!(obj::CBCFunction)
"""
cb_destroy!(obj::CBCFunction) = @ccall libcb.cb_cfunction_destroy(obj.data::Ptr{Cvoid})::Cvoid

struct CBNNCBoxSupportMinorantExtender <: CBMinorantExtender
    data::Ptr{Cvoid}
end
MA.mutability(::Type{CBNNCBoxSupportMinorantExtender}) = MA.IsMutable()

"""
    cb_destroy!(obj::CBNNCBoxSupportMinorantExtender)
"""
cb_destroy!(obj::CBNNCBoxSupportMinorantExtender) = @ccall libcb.cb_nncboxsupportminorantextender_destroy(obj.data::Ptr{Cvoid})::Cvoid

struct CBNNCBoxSupportFunction <: CBFunctionObject
    data::Ptr{Cvoid}
end
MA.mutability(::Type{CBNNCBoxSupportFunction}) = MA.IsMutable()

"""
    cb_destroy!(obj::CBNNCBoxSupportFunction)
"""
cb_destroy!(obj::CBNNCBoxSupportFunction) = @ccall libcb.cb_nncboxsupportfunction_destroy(obj.data::Ptr{Cvoid})::Cvoid

struct CBPSCAffineMinorantExtender <: CBMinorantExtender
    data::Ptr{Cvoid}
end
MA.mutability(::Type{CBPSCAffineMinorantExtender}) = MA.IsMutable()

"""
    cb_destroy!(obj::CBPSCAffineMinorantExtender)
"""
cb_destroy!(obj::CBPSCAffineMinorantExtender) = @ccall libcb.cb_pscaffineminorantextender_destroy(obj.data::Ptr{Cvoid})::Cvoid

abstract type CBPSCOracle <: CBModifiableOracleObject end

struct CBPSCAffineFunction <: CBPSCOracle
    data::Ptr{Cvoid}
end
MA.mutability(::Type{CBPSCAffineFunction}) = MA.IsMutable()

"""
    cb_destroy!(obj::CBPSCAffineFunction)
"""
cb_destroy!(obj::CBPSCAffineFunction) = @ccall libcb.cb_pscaffinefunction_destroy(obj.data::Ptr{Cvoid})::Cvoid

struct CBSOCSupportMinorantExtender <: CBMinorantExtender
    data::Ptr{Cvoid}
end
MA.mutability(::Type{CBSOCSupportMinorantExtender}) = MA.IsMutable()

"""
    cb_destroy!(obj::CBSOCSupportMinorantExtender)
"""
cb_destroy!(obj::CBSOCSupportMinorantExtender) = @ccall libcb.cb_socsupportminorantextender_destroy(obj.data::Ptr{Cvoid})::Cvoid

abstract type CBModificableOracleObject <: CBFunctionObject end

abstract type CBSOCOracle <: CBModificableOracleObject end

struct CBSOCSupportFunction <: CBSOCOracle
    data::Ptr{Cvoid}
end
MA.mutability(::Type{CBSOCSupportFunction}) = MA.IsMutable()

"""
    cb_destroy!(obj::CBSOCSupportFunction)
"""
cb_destroy!(obj::CBSOCSupportFunction) = @ccall libcb.cb_socsupportfunction_destroy(obj.data::Ptr{Cvoid})::Cvoid

abstract type CBBoxModelParametersObject <: CBBundleParameters end

struct CBBoxModelParameters <: CBBoxModelParametersObject
    data::Ptr{Cvoid}
end
MA.mutability(::Type{CBBoxModelParameters}) = MA.IsMutable()

"""
    cb_destroy!(obj::CBBoxModelParameters)
"""
cb_destroy!(obj::CBBoxModelParameters) = @ccall libcb.cb_boxmodelparameters_destroy(obj.data::Ptr{Cvoid})::Cvoid

abstract type CBNNCModelParametersObject <: CBBundleParameters end

struct CBNNCModelParameters <: CBNNCModelParametersObject
    data::Ptr{Cvoid}
end
MA.mutability(::Type{CBNNCModelParameters}) = MA.IsMutable()

"""
    cb_destroy!(obj::CBNNCModelParameters)
"""
cb_destroy!(obj::CBNNCModelParameters) = @ccall libcb.cb_nncmodelparameters_destroy(obj.data::Ptr{Cvoid})::Cvoid

abstract type CBPSCModelParametersObject <: CBBundleParameters end

struct CBPSCModelParameters <: CBPSCModelParametersObject
    data::Ptr{Cvoid}
end
MA.mutability(::Type{CBPSCModelParameters}) = MA.IsMutable()

"""
    cb_destroy!(obj::CBPSCModelParameters)
"""
cb_destroy!(obj::CBPSCModelParameters) = @ccall libcb.cb_pscmodelparameters_destroy(obj.data::Ptr{Cvoid})::Cvoid

abstract type CBSOCModelParametersObject <: CBBundleParameters end

struct CBSOCModelParameters <: CBSOCModelParametersObject
    data::Ptr{Cvoid}
end
MA.mutability(::Type{CBSOCModelParameters}) = MA.IsMutable()

"""
    cb_destroy!(obj::CBSOCModelParameters)
"""
cb_destroy!(obj::CBSOCModelParameters) = @ccall libcb.cb_socmodelparameters_destroy(obj.data::Ptr{Cvoid})::Cvoid

abstract type CBSumBundleParametersObject <: CBBundleParameters end

struct CBSumBundleParameters <: CBSumBundleParametersObject
    data::Ptr{Cvoid}
end
MA.mutability(::Type{CBSumBundleParameters}) = MA.IsMutable()

"""
    cb_destroy!(obj::CBSumBundleParameters)
"""
cb_destroy!(obj::CBSumBundleParameters) = @ccall libcb.cb_sumbundleparameters_destroy(obj.data::Ptr{Cvoid})::Cvoid

abstract type CBVariableMetricBundleData end

abstract type CBBundleData <: CBVariableMetricBundleData end

struct CBAFTData <: CBBundleData
    data::Ptr{Cvoid}
end
MA.mutability(::Type{CBAFTData}) = MA.IsMutable()

"""
    cb_destroy!(obj::CBAFTData)
"""
cb_destroy!(obj::CBAFTData) = @ccall libcb.cb_aftdata_destroy(obj.data::Ptr{Cvoid})::Cvoid

struct CBBoxData <: CBBundleData
    data::Ptr{Cvoid}
end
MA.mutability(::Type{CBBoxData}) = MA.IsMutable()

"""
    cb_destroy!(obj::CBBoxData)
"""
cb_destroy!(obj::CBBoxData) = @ccall libcb.cb_boxdata_destroy(obj.data::Ptr{Cvoid})::Cvoid

struct CBNNCData <: CBBundleData
    data::Ptr{Cvoid}
end
MA.mutability(::Type{CBNNCData}) = MA.IsMutable()

"""
    cb_destroy!(obj::CBNNCData)
"""
cb_destroy!(obj::CBNNCData) = @ccall libcb.cb_nncdata_destroy(obj.data::Ptr{Cvoid})::Cvoid

struct CBPSCData <: CBBundleData
    data::Ptr{Cvoid}
end
MA.mutability(::Type{CBPSCData}) = MA.IsMutable()

"""
    cb_destroy!(obj::CBPSCData)
"""
cb_destroy!(obj::CBPSCData) = @ccall libcb.cb_pscdata_destroy(obj.data::Ptr{Cvoid})::Cvoid

struct CBSOCData <: CBBundleData
    data::Ptr{Cvoid}
end
MA.mutability(::Type{CBSOCData}) = MA.IsMutable()

"""
    cb_destroy!(obj::CBSOCData)
"""
cb_destroy!(obj::CBSOCData) = @ccall libcb.cb_socdata_destroy(obj.data::Ptr{Cvoid})::Cvoid

abstract type CBBundleWeight end

struct CBBundleHKWeight <: CBBundleWeight
    data::Ptr{Cvoid}
end
MA.mutability(::Type{CBBundleHKWeight}) = MA.IsMutable()

"""
    cb_destroy!(obj::CBBundleHKWeight)
"""
cb_destroy!(obj::CBBundleHKWeight) = @ccall libcb.cb_bundlehkweight_destroy(obj.data::Ptr{Cvoid})::Cvoid

struct CBBundleRQBWeight <: CBBundleWeight
    data::Ptr{Cvoid}
end
MA.mutability(::Type{CBBundleRQBWeight}) = MA.IsMutable()

"""
    cb_destroy!(obj::CBBundleRQBWeight)
"""
cb_destroy!(obj::CBBundleRQBWeight) = @ccall libcb.cb_bundlerqbweight_destroy(obj.data::Ptr{Cvoid})::Cvoid

abstract type CBQPSolverProxObject end

abstract type CBBundleProxObject <: CBQPSolverProxObject end

struct CBBundleDenseTrustRegionProx <: CBBundleProxObject
    data::Ptr{Cvoid}
end
MA.mutability(::Type{CBBundleDenseTrustRegionProx}) = MA.IsMutable()

"""
    cb_destroy!(obj::CBBundleDenseTrustRegionProx)
"""
cb_destroy!(obj::CBBundleDenseTrustRegionProx) = @ccall libcb.cb_bundledensetrustregionprox_destroy(obj.data::Ptr{Cvoid})::Cvoid

struct CBBundleDiagonalTrustRegionProx <: CBBundleProxObject
    data::Ptr{Cvoid}
end
MA.mutability(::Type{CBBundleDiagonalTrustRegionProx}) = MA.IsMutable()

"""
    cb_destroy!(obj::CBBundleDiagonalTrustRegionProx)
"""
cb_destroy!(obj::CBBundleDiagonalTrustRegionProx) = @ccall libcb.cb_bundlediagonaltrustregionprox_destroy(obj.data::Ptr{Cvoid})::Cvoid

struct CBBundleDLRTrustRegionProx <: CBBundleProxObject
    data::Ptr{Cvoid}
end
MA.mutability(::Type{CBBundleDLRTrustRegionProx}) = MA.IsMutable()

"""
    cb_destroy!(obj::CBBundleDLRTrustRegionProx)
"""
cb_destroy!(obj::CBBundleDLRTrustRegionProx) = @ccall libcb.cb_bundledlrtrustregionprox_destroy(obj.data::Ptr{Cvoid})::Cvoid

struct CBBundleIdProx <: CBBundleProxObject
    data::Ptr{Cvoid}
end
MA.mutability(::Type{CBBundleIdProx}) = MA.IsMutable()

"""
    cb_destroy!(obj::CBBundleIdProx)
"""
cb_destroy!(obj::CBBundleIdProx) = @ccall libcb.cb_bundleidprox_destroy(obj.data::Ptr{Cvoid})::Cvoid

struct CBBundleLowRankTrustRegionProx <: CBBundleProxObject
    data::Ptr{Cvoid}
end
MA.mutability(::Type{CBBundleLowRankTrustRegionProx}) = MA.IsMutable()

"""
    cb_destroy!(obj::CBBundleLowRankTrustRegionProx)
"""
cb_destroy!(obj::CBBundleLowRankTrustRegionProx) = @ccall libcb.cb_bundlelowranktrustregionprox_destroy(obj.data::Ptr{Cvoid})::Cvoid

abstract type CBQPSolverParametersObject end

struct CBQPSolverParameters <: CBQPSolverParametersObject
    data::Ptr{Cvoid}
end
MA.mutability(::Type{CBQPSolverParameters}) = MA.IsMutable()

"""
    cb_destroy!(obj::CBQPSolverParameters)
"""
cb_destroy!(obj::CBQPSolverParameters) = @ccall libcb.cb_qpsolverparameters_destroy(obj.data::Ptr{Cvoid})::Cvoid

abstract type CBQPModelDataPointer end

abstract type CBQPSolverObject <: CBQPModelDataPointer end

struct CBQPSolver <: CBQPSolverObject
    data::Ptr{Cvoid}
end
MA.mutability(::Type{CBQPSolver}) = MA.IsMutable()

"""
    cb_destroy!(obj::CBQPSolver)
"""
cb_destroy!(obj::CBQPSolver) = @ccall libcb.cb_qpsolver_destroy(obj.data::Ptr{Cvoid})::Cvoid

struct CBUQPSolver <: CBQPSolverObject
    data::Ptr{Cvoid}
end
MA.mutability(::Type{CBUQPSolver}) = MA.IsMutable()

"""
    cb_destroy!(obj::CBUQPSolver)
"""
cb_destroy!(obj::CBUQPSolver) = @ccall libcb.cb_uqpsolver_destroy(obj.data::Ptr{Cvoid})::Cvoid

abstract type CBVariableMetricModel end

abstract type CBGroundset <: CBVariableMetricModel end

struct CBLPGroundset <: CBGroundset
    data::Ptr{Cvoid}
end
MA.mutability(::Type{CBLPGroundset}) = MA.IsMutable()

"""
    cb_destroy!(obj::CBLPGroundset)
"""
cb_destroy!(obj::CBLPGroundset) = @ccall libcb.cb_lpgroundset_destroy(obj.data::Ptr{Cvoid})::Cvoid

struct CBUnconstrainedGroundset <: CBGroundset
    data::Ptr{Cvoid}
end
MA.mutability(::Type{CBUnconstrainedGroundset}) = MA.IsMutable()

"""
    cb_destroy!(obj::CBUnconstrainedGroundset)
"""
cb_destroy!(obj::CBUnconstrainedGroundset) = @ccall libcb.cb_unconstrainedgroundset_destroy(obj.data::Ptr{Cvoid})::Cvoid

abstract type CBBundleModel <: CBVariableMetricModel end

abstract type CBSumBlockModel <: CBBundleModel end

struct CBAFTModel <: CBSumBlockModel
    data::Ptr{Cvoid}
end
MA.mutability(::Type{CBAFTModel}) = MA.IsMutable()

"""
    cb_destroy!(obj::CBAFTModel)
"""
cb_destroy!(obj::CBAFTModel) = @ccall libcb.cb_aftmodel_destroy(obj.data::Ptr{Cvoid})::Cvoid

abstract type CBConeModel <: CBSumBlockModel end

struct CBBoxModel <: CBConeModel
    data::Ptr{Cvoid}
end
MA.mutability(::Type{CBBoxModel}) = MA.IsMutable()

"""
    cb_destroy!(obj::CBBoxModel)
"""
cb_destroy!(obj::CBBoxModel) = @ccall libcb.cb_boxmodel_destroy(obj.data::Ptr{Cvoid})::Cvoid

struct CBNNCModel <: CBConeModel
    data::Ptr{Cvoid}
end
MA.mutability(::Type{CBNNCModel}) = MA.IsMutable()

"""
    cb_destroy!(obj::CBNNCModel)
"""
cb_destroy!(obj::CBNNCModel) = @ccall libcb.cb_nncmodel_destroy(obj.data::Ptr{Cvoid})::Cvoid

struct CBPSCModel <: CBConeModel
    data::Ptr{Cvoid}
end
MA.mutability(::Type{CBPSCModel}) = MA.IsMutable()

"""
    cb_destroy!(obj::CBPSCModel)
"""
cb_destroy!(obj::CBPSCModel) = @ccall libcb.cb_pscmodel_destroy(obj.data::Ptr{Cvoid})::Cvoid

struct CBSOCModel <: CBConeModel
    data::Ptr{Cvoid}
end
MA.mutability(::Type{CBSOCModel}) = MA.IsMutable()

"""
    cb_destroy!(obj::CBSOCModel)
"""
cb_destroy!(obj::CBSOCModel) = @ccall libcb.cb_socmodel_destroy(obj.data::Ptr{Cvoid})::Cvoid

struct CBSumModel <: CBSumBlockModel
    data::Ptr{Cvoid}
end
MA.mutability(::Type{CBSumModel}) = MA.IsMutable()

"""
    cb_destroy!(obj::CBSumModel)
"""
cb_destroy!(obj::CBSumModel) = @ccall libcb.cb_summodel_destroy(obj.data::Ptr{Cvoid})::Cvoid

abstract type CBVariableMetricSelection end

struct CBPSCVariableMetricSelection <: CBVariableMetricSelection
    data::Ptr{Cvoid}
end
MA.mutability(::Type{CBPSCVariableMetricSelection}) = MA.IsMutable()

"""
    cb_destroy!(obj::CBPSCVariableMetricSelection)
"""
cb_destroy!(obj::CBPSCVariableMetricSelection) = @ccall libcb.cb_pscvariablemetricselection_destroy(obj.data::Ptr{Cvoid})::Cvoid

struct CBVariableMetricSVDSelection <: CBVariableMetricSelection
    data::Ptr{Cvoid}
end
MA.mutability(::Type{CBVariableMetricSVDSelection}) = MA.IsMutable()

"""
    cb_destroy!(obj::CBVariableMetricSVDSelection)
"""
cb_destroy!(obj::CBVariableMetricSVDSelection) = @ccall libcb.cb_variablemetricsvdselection_destroy(obj.data::Ptr{Cvoid})::Cvoid

abstract type CBQPKKTSolverObject end

struct CBQPDirectKKTSolver <: CBQPKKTSolverObject
    data::Ptr{Cvoid}
end
MA.mutability(::Type{CBQPDirectKKTSolver}) = MA.IsMutable()

"""
    cb_destroy!(obj::CBQPDirectKKTSolver)
"""
cb_destroy!(obj::CBQPDirectKKTSolver) = @ccall libcb.cb_qpdirectkktsolver_destroy(obj.data::Ptr{Cvoid})::Cvoid

abstract type CBQPIterativeKKTSolver <: CBQPKKTSolverObject end

struct CBQPIterativeKKTHAeqSolver <: CBQPIterativeKKTSolver
    data::Ptr{Cvoid}
end
MA.mutability(::Type{CBQPIterativeKKTHAeqSolver}) = MA.IsMutable()

"""
    cb_destroy!(obj::CBQPIterativeKKTHAeqSolver)
"""
cb_destroy!(obj::CBQPIterativeKKTHAeqSolver) = @ccall libcb.cb_qpiterativekkthaeqsolver_destroy(obj.data::Ptr{Cvoid})::Cvoid

struct CBQPIterativeKKTHASolver <: CBQPIterativeKKTSolver
    data::Ptr{Cvoid}
end
MA.mutability(::Type{CBQPIterativeKKTHASolver}) = MA.IsMutable()

"""
    cb_destroy!(obj::CBQPIterativeKKTHASolver)
"""
cb_destroy!(obj::CBQPIterativeKKTHASolver) = @ccall libcb.cb_qpiterativekkthasolver_destroy(obj.data::Ptr{Cvoid})::Cvoid

abstract type CBQPKKTSOlverObject end

struct CBQPKKTSolverComparison <: CBQPKKTSOlverObject
    data::Ptr{Cvoid}
end
MA.mutability(::Type{CBQPKKTSolverComparison}) = MA.IsMutable()

"""
    cb_destroy!(obj::CBQPKKTSolverComparison)
"""
cb_destroy!(obj::CBQPKKTSolverComparison) = @ccall libcb.cb_qpkktsolvercomparison_destroy(obj.data::Ptr{Cvoid})::Cvoid

struct CBSumBundleHandler
    data::Ptr{Cvoid}
end
MA.mutability(::Type{CBSumBundleHandler}) = MA.IsMutable()

"""
    cb_destroy!(obj::CBSumBundleHandler)
"""
cb_destroy!(obj::CBSumBundleHandler) = @ccall libcb.cb_sumbundlehandler_destroy(obj.data::Ptr{Cvoid})::Cvoid

abstract type CBQPModelBlockObject end

abstract type CBQPModelBlock <: CBQPModelBlockObject end

struct CBQPConeModelBlock <: CBQPModelBlock
    data::Ptr{Cvoid}
end
MA.mutability(::Type{CBQPConeModelBlock}) = MA.IsMutable()

"""
    cb_destroy!(obj::CBQPConeModelBlock)
"""
cb_destroy!(obj::CBQPConeModelBlock) = @ccall libcb.cb_qpconemodelblock_destroy(obj.data::Ptr{Cvoid})::Cvoid

struct CBQPSumModelBlock <: CBQPModelBlock
    data::Ptr{Cvoid}
end
MA.mutability(::Type{CBQPSumModelBlock}) = MA.IsMutable()

"""
    cb_destroy!(obj::CBQPSumModelBlock)
"""
cb_destroy!(obj::CBQPSumModelBlock) = @ccall libcb.cb_qpsummodelblock_destroy(obj.data::Ptr{Cvoid})::Cvoid

abstract type CBUQPModelBlockObject end

abstract type CBUQPModelBlock <: CBUQPModelBlockObject end

struct CBUQPConeModelBlock <: CBUQPModelBlock
    data::Ptr{Cvoid}
end
MA.mutability(::Type{CBUQPConeModelBlock}) = MA.IsMutable()

"""
    cb_destroy!(obj::CBUQPConeModelBlock)
"""
cb_destroy!(obj::CBUQPConeModelBlock) = @ccall libcb.cb_uqpconemodelblock_destroy(obj.data::Ptr{Cvoid})::Cvoid

struct CBUQPSumModelBlock <: CBUQPModelBlock
    data::Ptr{Cvoid}
end
MA.mutability(::Type{CBUQPSumModelBlock}) = MA.IsMutable()

"""
    cb_destroy!(obj::CBUQPSumModelBlock)
"""
cb_destroy!(obj::CBUQPSumModelBlock) = @ccall libcb.cb_uqpsummodelblock_destroy(obj.data::Ptr{Cvoid})::Cvoid

abstract type CBIterativeSolverObject end

struct CBMinRes <: CBIterativeSolverObject
    data::Ptr{Cvoid}
end
MA.mutability(::Type{CBMinRes}) = MA.IsMutable()

"""
    cb_destroy!(obj::CBMinRes)
"""
cb_destroy!(obj::CBMinRes) = @ccall libcb.cb_minres_destroy(obj.data::Ptr{Cvoid})::Cvoid

struct CBPCG <: CBIterativeSolverObject
    data::Ptr{Cvoid}
end
MA.mutability(::Type{CBPCG}) = MA.IsMutable()

"""
    cb_destroy!(obj::CBPCG)
"""
cb_destroy!(obj::CBPCG) = @ccall libcb.cb_pcg_destroy(obj.data::Ptr{Cvoid})::Cvoid

struct CBPsqmr <: CBIterativeSolverObject
    data::Ptr{Cvoid}
end
MA.mutability(::Type{CBPsqmr}) = MA.IsMutable()

"""
    cb_destroy!(obj::CBPsqmr)
"""
cb_destroy!(obj::CBPsqmr) = @ccall libcb.cb_psqmr_destroy(obj.data::Ptr{Cvoid})::Cvoid

abstract type CBQPKKTPrecondObject end

struct CBQPKKTSubspaceHPrecond <: CBQPKKTPrecondObject
    data::Ptr{Cvoid}
end
MA.mutability(::Type{CBQPKKTSubspaceHPrecond}) = MA.IsMutable()

"""
    cb_destroy!(obj::CBQPKKTSubspaceHPrecond)
"""
cb_destroy!(obj::CBQPKKTSubspaceHPrecond) = @ccall libcb.cb_qpkktsubspacehprecond_destroy(obj.data::Ptr{Cvoid})::Cvoid

struct CBSumBundle
    data::Ptr{Cvoid}
end
MA.mutability(::Type{CBSumBundle}) = MA.IsMutable()

"""
    cb_destroy!(obj::CBSumBundle)
"""
cb_destroy!(obj::CBSumBundle) = @ccall libcb.cb_sumbundle_destroy(obj.data::Ptr{Cvoid})::Cvoid

struct CBBundleSolver
    data::Ptr{Cvoid}
end
MA.mutability(::Type{CBBundleSolver}) = MA.IsMutable()

"""
    cb_destroy!(obj::CBBundleSolver)
"""
cb_destroy!(obj::CBBundleSolver) = @ccall libcb.cb_bundlesolver_destroy(obj.data::Ptr{Cvoid})::Cvoid

struct CBBundleTerminator
    data::Ptr{Cvoid}
end
MA.mutability(::Type{CBBundleTerminator}) = MA.IsMutable()

"""
    cb_destroy!(obj::CBBundleTerminator)
"""
cb_destroy!(obj::CBBundleTerminator) = @ccall libcb.cb_bundleterminator_destroy(obj.data::Ptr{Cvoid})::Cvoid

struct CBClock
    data::Ptr{Cvoid}
end
MA.mutability(::Type{CBClock}) = MA.IsMutable()

"""
    cb_destroy!(obj::CBClock)
"""
cb_destroy!(obj::CBClock) = @ccall libcb.cb_clock_destroy(obj.data::Ptr{Cvoid})::Cvoid

struct CBMicroseconds
    data::Ptr{Cvoid}
end
MA.mutability(::Type{CBMicroseconds}) = MA.IsMutable()

"""
    cb_destroy!(obj::CBMicroseconds)
"""
cb_destroy!(obj::CBMicroseconds) = @ccall libcb.cb_microseconds_destroy(obj.data::Ptr{Cvoid})::Cvoid

const CBCoeffmatVector = StdVector{<:CBCoeffmat}

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

abstract type CBQPPSCOracleDataObject <: CBQPModelOracleDataObject end