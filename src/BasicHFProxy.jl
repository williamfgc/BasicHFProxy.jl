module BasicHFProxy

import SpecialFunctions # erf

function datafilter(fname)
    contains(fname, '_') && return false
    endswith(fname, ".inp") && return false
    return true
end

const DATADIR = joinpath(dirname(pathof(BasicHFProxy)), "../data")
const DATA = Dict(Symbol(f) => joinpath(DATADIR, f)
                  for f in filter(datafilter, readdir(DATADIR)))

include("common.jl")
include("sequential.jl")
include("threads.jl")

export bhfp_sequential, bhfp_threads

end
