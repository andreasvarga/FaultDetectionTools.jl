module FaultDetectionTools

using LinearAlgebra
using MatrixEquations
using MatrixPencils
using DescriptorSystems
using Polynomials
using Random

abstract type AbstractFDDObject end

import LinearAlgebra: BlasFloat, BlasReal, BlasComplex
import DescriptorSystems: DescriptorStateSpace, chess, rcond, gbalmr, gminreal, gpole, isstable
import Combinatorics: combinations
import Base: +, -, *

export FDIModel, fdimodset, FDFilter, FDFilterIF, FDIFilter, FDIFilterIF, fdIFeval, gbalmr
export fditspec, fdisspec, fdiscond, fdif2ngap, fdimmperf
export fdigenspec, fdichkspec
export efdsyn, efdisyn, efdbasesel, afdsyn, afdredsyn, afdisyn, afdbasesel, emmsyn, emmbasesel, ammsyn, ammbasesel
export fdhinfminus, fdhinfmax, binmat2dec, dec2binmat, fditspec_, fdisspec_, fdiscond_

export MDMModel, MDModel, mdmodset, MDFilter, MDFilterIF, mdIFeval
export mddist, mddist2c
export mdspec, mdsspec, mdperf, mdmatch, mdgap
export emdsyn, amdsyn, emdbasesel


const VRS = Union{Vector{Int}, UnitRange{Int}, Int}

include("types/FDIObjects.jl")
include("FDIanalysis.jl")
include("FDIperformance.jl")
include("FDIsynthesis.jl")
include("FDIutils.jl")
include("types/MDObjects.jl")
include("MDanalysis.jl")
include("MDperformance.jl")
include("MDsynthesis.jl")
end