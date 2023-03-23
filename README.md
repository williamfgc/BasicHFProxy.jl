# BasicHFProxy.jl

Julia implementation of a basic Hartree-Fock proxy application

## Example

```julia
julia> using BasicHFProxy

julia> he4 = joinpath(dirname(pathof(BasicHFProxy)), "../data/he4");

julia> bhfp_sequential(he4);
2e- energy= 4.050176411152184
```
