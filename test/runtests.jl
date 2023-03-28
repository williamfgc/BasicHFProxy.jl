using Test
using BasicHFProxy
using MPI
@show MPI.identify_implementation()

const TESTDIR = @__DIR__
const MPI_TESTDIR = joinpath(TESTDIR, "mpi")
filterfunc(x) = startswith(x, "mpitest_") && endswith(x, ".jl")
const MPI_TESTFILES = joinpath.(MPI_TESTDIR,
                                sort(filter(filterfunc, readdir(MPI_TESTDIR))))
const MPI_NPROCS = clamp(Sys.CPU_THREADS, 2, 5)

# include all tests_*.jl files from the test/ directory
for f in filter(startswith("tests_"), readdir(@__DIR__))
    # !contains(f, "mpi") && continue
    include(f)
end
