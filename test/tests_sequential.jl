@testset "Sequential" begin
# considering all input files would take too long
for s in (:he4, :he8)
    f = BasicHFProxy.DATA[s]
    @testset "$s" begin
        @test bhfp_sequential(f) isa Float64
        @test bhfp_sequential(f) ≈ BasicHFProxy.expected_energy(f)
    end
end end
