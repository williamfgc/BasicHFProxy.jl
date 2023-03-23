using Test
using BasicHFProxy

# include all tests_*.jl files from the test/ directory
for f in filter(startswith("tests_"), readdir(@__DIR__))
    include(f)
end
