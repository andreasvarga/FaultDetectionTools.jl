# Release Notes

## Version 0.6.0 

This minor release concludes the implementation of the functions for solving synthesis problems of _fault detection and isolation_ (FDI) filters. This release provides a new function for the approximate model-matching based synthesis of FDI filters and three new functions for the evaluation of the model-matching performance. The presented examples in Chapters 5 and 7, and the case studies in Chapter 8 of the book 
"[A. Varga, Solving Fault Diagnosis Problems, Linear Synthesis Techniques, Springer 2017](https://www.springer.com/us/book/9783319515588)", can be fully reproduced using the collection of scripts provided in the `src/Examples` directory.  This release relies on the new enhanced version `v1.2` of the [DescriptorSystems.jl](https://github.com/andreasvarga/DescriptorSystems.jl) package. Several tutorials have been prepared to complement the information available in the above book.

## Version 0.5.0 

This minor release provides a new function for the model-matching based synthesis of fault detection and isolation filters and two new functions for the evaluation of the model-matching performance. The series coupling of a descriptor system model and fault detection objects of type `FDIModel` and `FDFilterIF`, as well as for the component-wise series coupling of a vector of descriptor system models and a fault detection object of type `FDIFilterIF` are now supported. 

## Version 0.4.0 

This minor release provides a new function for the approximate synthesis of fault detection and isolation filters and three new functions for the evaluation of the achieved fault-to-noise gap and model-matching performace. The functions of the DescriptorSystems package  `gminreal` and `gpole` have been extended to cover FDIFilter objects as well. 

## Version 0.3.0 

This minor release provides a new function for the approximate synthesis of fault detection filters and two new functions for the evaluation of the achieved fault-to-noise gap and model-matching performace. The functions of the DescriptorSystems package  `gbalmr`, `gminreal` and `gpole` have been extended to cover FDFilter objects as well. 

## Version 0.2.0 

This minor release provides several new functions to define FDI filter objects and associated constructors, functions for the generation of achievable specifications and for the assessment of feasibility of a set of specifications, a synthesis function for the exact solution of synthesis problems of fault detection and isolation filters, and functions for the performance analysis of FDI filters. The functionality of functions `fditspec` and `fdisspec` to evaluate strong structure matrices has been simplified.

## Version 0.1.0

This is the initial release providing prototype implementations of several basic fault detection related objects and their constructors, jointly with functions to perform basic performance analysis and for the exact solution of synthesis problems of fault detection filters.  
