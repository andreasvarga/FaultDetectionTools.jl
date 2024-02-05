# Julia scripts for the examples presented in the book "Solving Fault Diagnosis Problems"

## Examples and case studies

This collection contains the Julia scripts which allow to execute all computational examples and case studies presented in Chapters 5-8 of the book:

**Andreas Varga**, [Solving Fault Diagnosis Problems - Linear Synthesis Techniques with Julia Code Examples](https://link.springer.com/book/10.1007/978-3-031-35767-1), vol. 482 of Studies in Systems, Decision and Control, Springer International Publishing, 2024.


| Examples | Description |
| :--- | :--- |
| **Ex2_4**  |         Example 2.4 - Conversion of a LPV model to LTI form |
| **Ex5_4**  |         Example 5.4 - Solution of an exact fault detection problem (EFDP) |
| **Ex5_4c**  |      Compact version of script Ex5_4 using the function efdsyn of FaultDetectionTools |
| **Ex5_6**  |         Example 5.6 - Solution of an approximate fault detection problem (AFDP) |
| **Ex5_6c**  |       Compact version of script Ex5_6 using the function afdsyn of FaultDetectionTools |
| **Ex5_10**  |       Example 5.10 - Solution of an exact fault detection and isolation problem (EFDIP) |
| **Ex5_10c**  |     Compact version of script Ex5_10 using the function efdisyn of FaultDetectionTools |
| **Ex5_11**  |       Example 5.11 - Solution of an approximat fault detection and isolation problem (AFDIP) |
| **Ex5_11c**  |     Compact version of script Ex5_11 using the function afdisyn of FaultDetectionTools |
| **Ex5_12**  |       Example 5.12 - Solution of an exact model-matching problem (EMMP) |
| **Ex5_12c**  |     Compact version of script Ex5_12 using the function emmsyn of  FaultDetectionTools |
| **Ex5_13**  |        Example 5.13 - Solution of an exact model-matching problem (EMMP) |
| **Ex5_13a**  |     Alternative direct approach to solve the EMMP for Example 5.13 |
| **Ex5_13c**  |     Compact version of script Ex5_13 using the function emmsyn of  FaultDetectionTools |
| **Ex5_16**  |       Example 5.16 - Solution of an  H∞  approximate model-matching problem (AMMP) |
| **Ex5_16c**  |    Compact version of script Ex5_16 using the function ammsyn of  FaultDetectionTools |
| **Ex6_1**  |          Example 6.1 - Solution of an exact model detection problem (EMDP) |
| **Ex6_1c**  |        Compact version of script Ex6_1 using the function emdsyn of  FaultDetectionTools |
| **Ex6_2**  |          Example 6.2 - Solution of an approximate model detection problem (AMDP) |
| **Ex6_2c**  |       Compact version of script Ex6_2 using the function amdsyn of  FaultDetectionTools |
| **Ex6_2KF**  |   Kalman-filter based synthesis to address the AMMP for Example 6.2 |
| **Ex7_1**  |         Example 7.1 - Illustrating polynomial root sensitivity |
| **Ex7_3**  |         Example 7.3 - Nullspace-based solution of an EFDP |
| **Ex7_4**  |         Example 7.4 - Nullspace-based least order synthesis to solve an EFDP |

| Case studies | Description |
| :--- | :--- |
| **CS1_1**  |        Case-study 1.1 - Monitoring flight actuator faults (no local surface angle measurements) |
| **CS1_2**  |        Case-study 1.2 - Monitoring flight actuator faults with local surface angle measurements |
| **CS2_1**  |        Case study 2.1 - Monitoring air data sensor faults (robust least order LTI synthesis) |
| **CS2_2**  |        Case study 2.2 - Monitoring air data sensor faults (robust least order LPV synthesis)v

## How to execute the scripts

For the execution of scripts the Julia version 1.8 or higher must be installed, together with the packages **FaultDetectionTools**, **DescriptorSystems**, **MatrixPencils**, **MatrixEquations**, **Polynomials**, **Measurements**, **GenericLinearAlgebra**, **Makie**, **CairoMakie**, **LaTeXStrings**, **JLD2** and **Optim**. 

To execute all examples in Chapters 5, 6 and 7, and all case studies in Chapter 8, execute 

````JULIA
using FaultDetectionTools
cd(joinpath(pkgdir(FaultDetectionTools), "Examples"))
include("runexamples.jl")
include("runcasestudies.jl")
````
To execute a particular example, say Example 5.4 and its compact variant 5.4c, execute 

````JULIA
using FaultDetectionTools
cd(joinpath(pkgdir(FaultDetectionTools), "Examples"))
include("Ex5_4.jl")
include("Ex5_4c.jl")
````
To execute a particular case study, say Case Study 2.1, execute 

````JULIA
using FaultDetectionTools
cd(joinpath(pkgdir(FaultDetectionTools), "Examples"))
include("CS2_1.jl")
````
