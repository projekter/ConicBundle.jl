```@meta
CurrentModule = ConicBundle
```
# Introduction

`ConicBundle.jl` is the Julia port to the [ConicBundle](https://www-user.tu-chemnitz.de/~helmberg/ConicBundle) library.
The functions that are centered around the [`CBProblem`](@ref) type address the official C interface of ConicBundle; they were
translated by hand, adapted to a Julian way of coding and should work without any problems.
However, the C interface is quite restricted and ConicBundle has a lot more to offer. For this reason, a lot of the
functionality of the C++ interface was also exported for this Julia port and can be accessed by the C++ interface. Do not mix
those two interfaces! The C interface has automatic memory management, the C++ interface does not. Additionally, the latter was
created automatically using a script and lots of regular expressions; it comes without any guarantees that the functions will
even be callable. Furthermore, the way of just exporting everything is not particularly efficient, as Julia would be able to do
a lot of the functionality that is inlined in C++ on its own.

# Overview
```@docs
ConicBundle
```