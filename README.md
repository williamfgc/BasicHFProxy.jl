# BasicHFProxy.jl

Julia implementation of a basic Hartree-Fock proxy application

## Example

```julia
julia> using BasicHFProxy

julia> keys(BasicHFProxy.DATA)
KeySet for a Dict{Symbol, String} with 9 entries. Keys:
  :he128
  :he16
  :he256
  :he4
  :he512
  :he1024
  :he32
  :he8
  :he64

julia> he4 = BasicHFProxy.DATA[:he4];

julia> E = bhfp_sequential(he4)
4.050176411152184

julia> BasicHFProxy.expected_energy(he4) # parsed from data file
4.0501763971342815

julia> E ≈ BasicHFProxy.expected_energy(he4)
true
```
