module FaultDetectionTools

using LinearAlgebra
using MatrixEquations
using MatrixPencils
using DescriptorSystems
using Polynomials
using Random

abstract type AbstractFDDObject end

import LinearAlgebra: BlasFloat, BlasReal, BlasComplex
import DescriptorSystems: DescriptorStateSpace, chess, rcond, gbalmr, gminreal, gpole
import Combinatorics: combinations
import Base: +, -, *

export FDIModel, fdimodset, FDFilter, FDFilterIF, FDIFilter, FDIFilterIF, fdIFeval, gbalmr
export fditspec, fdisspec, fdiscond, fdif2ngap, fdimmperf
export fdigenspec, fdichkspec
export efdsyn, efdisyn, efdbasesel, afdsyn, afdredsyn, afdbasesel
export fdhinfminus, fdhinfmax, binmat2dec, dec2binmat, fditspec_, fdisspec_, fdiscond_

const VRS = Union{Vector{Int}, UnitRange{Int}, Int}

include("types/FDIObjects.jl")
include("FDIanalysis.jl")
include("FDIperformance.jl")
include("FDIsynthesis.jl")
include("FDIutils.jl")
end