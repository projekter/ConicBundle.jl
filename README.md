# ConicBundle for Julia

[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://projekter.github.io/ConicBundle.jl/dev)

This repository contains Julia bindings of Christoph Helmberg's
[ConicBundle](https://www-user.tu-chemnitz.de/~helmberg/ConicBundle/Manual/index.html), version 1.a.2.
The source code itself is almost unchanged, but an additional set of export functions have been provided for a lot of the C++
interface.
It should be used in conjuction with my [GitHub clone](../ConicBundle), which contains the modifications necessary for these
bindings. The binary artifacs are shipped with the Julia package.

# License
The bindings are put under the same license as ConicBundle itself, which currently is GPL v3.0 or later.