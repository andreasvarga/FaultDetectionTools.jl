module FaultDetectionTools

using LinearAlgebra
using MatrixEquations
using MatrixPencils
using DescriptorSystems
using Polynomials
using Random

abstract type AbstractFDDObject end

import LinearAlgebra: BlasFloat, BlasReal, BlasComplex
import DescriptorSystems: DescriptorStateSpace, chess, rcond, gbalmr
import Combinatorics: combinations

export FDIModel, fdimodset, FDFilter, FDFilterIF, FDIFilter, FDIFilterIF, fdIFeval, gbalmr
export fditspec, fdisspec, fdscond, fdif2ngap
export fdigenspec, fdichkspec
export efdsyn, efdisyn, efdbasesel, afdredsyn, afdbasesel
export fdhinfminus, fdhinfmax, binmat2dec, dec2binmat, fditspec_, fdisspec_, fdscond_

const VRS = Union{Vector{Int}, UnitRange{Int}, Int}

include("types/FDIObjects.jl")
include("FDIanalysis.jl")
include("FDIperformance.jl")
include("FDIsynthesis.jl")
include("FDIutils.jl")
end