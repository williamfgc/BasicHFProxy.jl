using BasicHFProxy
using MPI

MPI.Init()
comm = MPI.COMM_WORLD
rank = MPI.Comm_rank(comm)

if rank == 0
    printstyled("Input: $(ARGS[1])\n"; bold = true, color = :white)
end
# MPI.Barrier(comm)

t_start = MPI.Wtime()
E = bhfp_mpi()
t_stop = MPI.Wtime()

if rank == 0
    println("E = ", E)
    println("Time: ", round(t_stop - t_start; digits = 5), "sec\n")
end
