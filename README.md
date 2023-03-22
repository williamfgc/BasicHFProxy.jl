# BasicHFProxy.jl

Julia implementation of a basic Hartree-Fock proxy application

## Example

```julia
julia> using BasicHFProxy

julia> h4 = joinpath(dirname(pathof(BasicHFProxy)), "../data/he4");

julia> bhfp_sequential(h4);
2e- energy= 4.050176411152184
```
