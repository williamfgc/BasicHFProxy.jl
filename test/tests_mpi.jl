@testset "MPI" begin
    @testset "$f" for f in MPI_TESTFILES
        mpiexec() do mpirun
            function cmd(n = MPI_NPROCS)
                `$mpirun -n $n $(Base.julia_cmd()) --startup-file=no $f`
            end
            @test success(cmd())
            # run(cmd())
            # @test true
        end
    end
end
