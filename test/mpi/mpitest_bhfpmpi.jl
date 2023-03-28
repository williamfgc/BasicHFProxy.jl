using BasicHFProxy
using MPI
using Test

MPI.Init()
comm = MPI.COMM_WORLD
rank = MPI.Comm_rank(comm)
for s in (:he4, :he8)
    f = BasicHFProxy.DATA[s]
    @testset "$s" begin
        @test bhfp_mpi(f) isa Float64
        # @test bhfp_mpi(f) ≈ bhfp_sequential(f)
        @test bhfp_mpi(f) ≈ BasicHFProxy.expected_energy(f)
    end
end
