module BasicHFProxy

import SpecialFunctions # erf
import Atomix
import MPI

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
include("threads_lock.jl")
include("threads_atomix.jl")
include("BasicHFPreferences.jl")
include("jacc.jl")
include("mpi.jl")

export bhfp_sequential,
       bhfp_threads_lock,
       bhfp_threads_atomix,
       bhfp_mpi,
       bhfp_jacc

end
