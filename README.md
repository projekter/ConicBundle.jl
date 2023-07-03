# ConicBundle for Julia

[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://projekter.github.io/ConicBundle.jl/stable)

This repository contains Julia bindings of Christoph Helmberg's
[ConicBundle](https://www-user.tu-chemnitz.de/~helmberg/ConicBundle/Manual/index.html), version 1.a.2.
The source code itself is almost unchanged, but an additional set of export functions have been provided for a lot of the C++
interface.
Note that the C++ part was generated automatically; a manual cleanup would certainly be desirable.
A project file for Visual Studio 2022 was added and the binary for Windows x64 are provided. While it was not tried to compile
the project on Linux, this should be very doable, as this was ConicBundle's original platform.

# Compilation
ConicBundle itself with some modifications is found in the `ConicBundle` subdirectory.
For the compilation of ConicBundle under Linux, just run `make` in this directory, which will create both a static library
`libcb.a` in `ConicBundle/lib` as well as the shared object file `ConicBundle.so` that the Julia interface uses in `bin`.
For more options, see the [the ConicBundle readme](ConicBundle/README).
For Windows, the repository comes with a Visual Studio 2022 project file which you can just open and compile. The output DLL
(in Release mode) should go to the `bin` directory automatically.
Both shared files are prepackaged.
If you want to regenerate the automatically created C++ interface, move your working directory to `src/cppinterface` and run
the `adapter.py` script (requires the `regex` package, which can be easily installed using pip). This will recreate all Julia
files in the `src/cppinterface` directory as well as the corresponding C++ files in `ConicBundle/cppinterface`.

# License
The bindings are put under the same license as ConicBundle itself, which currently is GPL v3.0 or later.