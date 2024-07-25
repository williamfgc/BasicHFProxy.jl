import JACC
using Test
using BasicHFProxy

@testset "JACC" begin
    for s in (:he4, :he8, :he16, :he32)
        f = BasicHFProxy.DATA[s]
        @testset "$s" begin
            @test bhfp_jacc(f) isa Float64
            @test bhfp_jacc(f) ≈ bhfp_sequential(f)
            @test bhfp_jacc(f) ≈ BasicHFProxy.expected_energy(f)
        end
    end
end
