# GrobCovBridge.jl

[日本語](README.ja.md)

> **⚠ EXPERIMENTAL — NO MAINTENANCE GUARANTEE**
>
> This package is a temporary bridge and is **not registered** in the Julia General registry.
> It exists solely to expose `grobcov.lib` to OSCAR.jl users until the functionality is
> integrated upstream. Expect breaking changes without notice. Do not depend on it in
> production code.

A Julia bridge that wraps Singular's [`grobcov.lib`](https://www.singular.uni-kl.de/Manual/4-0-2/sing_956.htm)
for use with [OSCAR.jl](https://docs.oscar-system.org/) polynomial ideals.

Internally it uses `Singular.jl`'s `low_level_caller_ring` to call the Singular interpreter
and handles the conversion between Oscar's `MPolyIdeal` and Singular's parametric ring
`Q(params)[main_vars]` required by `grobcov.lib`.

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/xiupos/GrobCovBridge.jl")
```

Or in development mode from a local clone:

```julia
Pkg.develop(path="/path/to/GrobCovBridge.jl")
```

## Quick start

```julia
using Oscar, GrobCovBridge

R, (a, x, y) = polynomial_ring(QQ, ["a", "x", "y"])
I = ideal(R, [a*x - 1, a*y - 1])

# Compute the Groebner cover; [a] are the parameter variables, [x, y] are the main variables.
result = grobcov(I, [a], [x, y])

# result is a Vector{GrobCovSegment}
for seg in result
    println("lpp:     ", seg.lpp)
    println("basis:   ", gens(seg.basis))
    println("segment: ", seg.segment)
    println()
end
```

### Two-parameter example

```julia
R, (a, b, x, y) = polynomial_ring(QQ, ["a", "b", "x", "y"])
I = ideal(R, [a*x^2 + b*y, x - b*y])
result = grobcov(I, [a, b], [x, y])
```

### Comprehensive Groebner system (cgsdr)

```julia
R, (a, x) = polynomial_ring(QQ, ["a", "x"])
I = ideal(R, [a*x - 1])
result = cgsdr(I, [a], [x])

# result is a Vector{CGSSegment}
for seg in result
    println("null:    ", gens(seg.null))
    println("nonnull: ", gens(seg.nonnull))
    println("basis:   ", gens(seg.basis))
    println()
end
```

## API

### `grobcov(I, params, main_vars; options=Int[])`

Compute the Groebner cover of the parametric ideal `I`.

- `I`: an `Oscar.MPolyIdeal` over `QQ`
- `params`: vector of generators of `base_ring(I)` to use as parameters
- `main_vars`: vector of generators to use as polynomial variables
- `options`: optional list of integer flags passed to the Singular procedure

Returns `Vector{GrobCovSegment}`.

#### `GrobCovSegment`

| Field | Type | Description |
|-------|------|-------------|
| `lpp` | `Vector{MPolyRingElem}` | Leading power products of the basis |
| `basis` | `MPolyIdeal` | Reduced Groebner basis on this segment |
| `segment` | `Vector` | P-representation: `[[prime_ideal, [excl_ideal, ...]], ...]` |

The segment describes a constructible set
`∪_i (V(prime_i) \ ∪_j V(excl_ij))` in parameter space.

### `cgsdr(I, params, main_vars; options=Int[])`

Compute a comprehensive Groebner system with disjoint reduced representation.

Returns `Vector{CGSSegment}`.

#### `CGSSegment`

| Field | Type | Description |
|-------|------|-------------|
| `null` | `MPolyIdeal` | Ideal whose variety is the null set of the segment |
| `nonnull` | `MPolyIdeal` | Ideal whose variety must be avoided |
| `basis` | `MPolyIdeal` | Groebner basis on this segment |

## Requirements

- Julia ≥ 1.10
- Oscar.jl ≥ 1.0 (which includes Singular.jl and `grobcov.lib`)

## License

[GPL v3+](LICENSE)
