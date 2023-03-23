if Threads.nthreads == 1
    @warn("Running tests_threads.jl with 1 Thread.")
end

@testset "Threads" begin
    for bhfp in (bhfp_threads_lock, bhfp_threads_atomix)
        @testset "$bhfp" begin
            # considering all input files would take too long
            for s in (:he4, :he8)
                f = BasicHFProxy.DATA[s]
                @testset "$s" begin
                    @test bhfp(f) isa Float64
                    @test bhfp(f) ≈ bhfp_sequential(f)
                    @test bhfp(f) ≈ BasicHFProxy.expected_energy(f)
                end
            end
        end
    end
end
