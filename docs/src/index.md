```@meta
CurrentModule = RAS_DMFT
```

# Real frequency DMFT solver

`RAS_DMFT.jl` can be used to calculate a single Hubbard band using DMFT on the real frequency axis.

Given a single impurity Anderson Hamiltonian (SIAM) $H=H_0 + H_\mathrm{int}$ with

```math
H_0
= \sum_\sigma \epsilon_d d_\sigma^\dagger d_\sigma
+ \sum_{k\sigma} \epsilon_{k\sigma} c_{k\sigma}^\dagger c_{k\sigma}
+ \sum_{k\sigma} (V_{k\sigma} d_\sigma^\dagger c_{k\sigma}
  + V_{k\sigma}^* c_{k\sigma}^\dagger d_\sigma),
```

and interaction

```math
H_\mathrm{int} = U d_\uparrow^\dagger d_\downarrow^\dagger d_\downarrow d_\uparrow,
```

it will calculate the retarded Green's function on the real axis

```math
\begin{aligned}
G(t) &= - \mathrm{i} \Theta(t) \langle \{ d_\alpha^\dagger(t), d_\alpha \} \rangle \\
G(ω) &= \lim_{δ→0^+} \int_{-∞}^∞ G(t) \mathrm{e}^{\mathrm{i}(ω+\mathrm{i}δ)t} \mathrm{d}t.
\end{aligned}
```

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

## Modules

The main module is called `RAS_DMFT` and can be put into the namespace by

```julia
using RAS_DMFT
```

There are also submodules

- `RAS_DMFT.Combinatorics`
- `RAS_DMFT.Debug`

which are not necessary to be called in 99 % of all cases.
