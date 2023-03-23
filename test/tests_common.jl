const he4 = BasicHFProxy.DATA[:he4]
const he32 = BasicHFProxy.DATA[:he32]

@testset "Input file parsing" begin
    @testset "he4" begin
        input = BasicHFProxy.parse_input_file(he4)
        @test typeof(input) ==
              NamedTuple{(:ngauss, :natom, :xpnt, :coef, :geom),
                         Tuple{Int64, Int64, Vector{Float64}, Vector{Float64},
                               Matrix{Float64}}}
        @test input.ngauss == 3
        @test input.natom == 4
        @test input.xpnt == [6.3624214, 1.1589230, 0.3136498]
        @test input.coef == [0.154328967295, 0.535328142282, 0.444634542185]
        @test input.geom == [0.0 0.0 0.0
                             0.05 0.0 1.0
                             0.1 1.0 0.0
                             1.0 0.2 0.0]'
    end
    @testset "he32" begin
        input = BasicHFProxy.parse_input_file(he32)
        @test typeof(input) ==
              NamedTuple{(:ngauss, :natom, :xpnt, :coef, :geom),
                         Tuple{Int64, Int64, Vector{Float64}, Vector{Float64},
                               Matrix{Float64}}}
        @test input.ngauss == 3
        @test input.natom == 32
        @test input.xpnt == [6.3624214, 1.1589230, 0.3136498]
        @test input.coef == [0.154328967295, 0.535328142282, 0.444634542185]
        @test input.geom ==
              [0.0000 0.0000 0.0000
               0.0000 0.0000 1.4000
               0.0000 0.0000 2.8000
               0.0000 0.0000 4.2000
               0.0000 1.4000 0.0000
               0.0000 1.4000 1.4000
               0.0000 1.4000 2.8000
               0.0000 1.4000 4.2000
               0.0000 2.8000 0.0000
               0.0000 2.8000 1.4000
               0.0000 2.8000 2.8000
               0.0000 2.8000 4.2000
               0.0000 4.2000 0.0000
               0.0000 4.2000 1.4000
               0.0000 4.2000 2.8000
               0.0000 4.2000 4.2000
               1.4000 0.0000 0.0000
               1.4000 0.0000 1.4000
               1.4000 0.0000 2.8000
               1.4000 0.0000 4.2000
               1.4000 1.4000 0.0000
               1.4000 1.4000 1.4000
               1.4000 1.4000 2.8000
               1.4000 1.4000 4.2000
               1.4000 2.8000 0.0000
               1.4000 2.8000 1.4000
               1.4000 2.8000 2.8000
               1.4000 2.8000 4.2000
               1.4000 4.2000 0.0000
               1.4000 4.2000 1.4000
               1.4000 4.2000 2.8000
               1.4000 4.2000 4.2000]'
    end
end
