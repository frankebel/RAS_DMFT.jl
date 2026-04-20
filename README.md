# RAS_DMFT

[![Build Status](https://github.com/frankebel/RAS_DMFT.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/frankebel/RAS_DMFT.jl/actions/workflows/CI.yml)
[![codecov](https://codecov.io/gh/frankebel/RAS_DMFT.jl/graph/badge.svg?token=5ACAMMA64E)](https://codecov.io/gh/frankebel/RAS_DMFT.jl)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)
[![code style: runic](https://img.shields.io/badge/code_style-%E1%9A%B1%E1%9A%A2%E1%9A%BE%E1%9B%81%E1%9A%B2-black)](https://github.com/fredrikekre/Runic.jl)

Restricted active space DMFT solver on the real frequency axis.

## Installation

As the package is not inside the [General registry](https://github.com/JuliaRegistries/General),
it needs to be added
[manually](https://pkgdocs.julialang.org/v1/managing-packages/#Adding-unregistered-packages).

```sh
export JULIA_PKG_USE_CLI_GIT="true"
julia --project=path/to/project --eval 'using Pkg; Pkg.add(url="https://github.com/frankebel/RAS_DMFT.jl")'
```

If the package is installed, you can run all tests with

```julia
julia --project=path/to/project --eval 'using Pkg; Pkg.test("RAS_DMFT")'
```

## Documentation

The documentation resides in `docs`.
Currently, it needs to be compiled manually

```sh
julia --project=path/to/package/docs --exec 'using Pkg; Pkg.develop(path=".."); Pkg.instantiate()'
julia --project=path/to/package/docs make.jl
```

It can then be viewed with, e.g. [LiveServer.jl](https://github.com/tlienart/LiveServer.jl)

```julia
using LiveServer
serve(dir="build")
```
