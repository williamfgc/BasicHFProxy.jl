
import MPI

include("Common.jl")

function main(args)
    MPI.Init()
    myrank = MPI.Comm_rank(MPI.COMM_WORLD)
    nproc = MPI.Comm_size(MPI.COMM_WORLD)

    if myrank == 0
        println("number of ranks: ", nproc)
    end

    input_file = get_input_file(args)

    xpnt = Array{Float64, 1}(undef, ngauss)
    coef = Array{Float64, 1}(undef, ngauss)
    geom = Array{Float64, 2}(undef, 3, natom)
end
