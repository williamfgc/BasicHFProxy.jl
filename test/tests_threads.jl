if Threads.nthreads == 1
    @warn("Running tests_threads.jl with 1 Thread.")
end

@testset "Threads" begin
# considering all input files would take too long
for s in (:he4, :he8)
    f = BasicHFProxy.DATA[s]
    @testset "$s" begin
        @test bhfp_threads(f) isa Float64
        @test bhfp_threads(f) ≈ bhfp_sequential(f)
        @test bhfp_threads(f) ≈ BasicHFProxy.expected_energy(f)
    end
end end
