module FaultDetectionTools

using LinearAlgebra
using MatrixEquations
using MatrixPencils
using DescriptorSystems
using Polynomials
using Random

abstract type AbstractDescriptorStateSpace end

import LinearAlgebra: BlasFloat, BlasReal, BlasComplex
import DescriptorSystems: DescriptorStateSpace, rcond
import Combinatorics: combinations

export FDIModel, fdimodset, FDFilter, FDFilterIF, FDIFilter, FDIFilterIF
export fditspec, fdisspec, fdscond
export efdsyn, efdbasesel
export fdIFeval, fdhinfminus, fdhinfmax

const VRS = Union{Vector{Int}, UnitRange{Int}, Int}

include("types/FDIObjects.jl")
include("FDIperformance.jl")
include("FDIsynthesis.jl")
include("FDIutils.jl")
end