using BasicHFProxy

printstyled("Warmup run...")
input = BasicHFProxy.DATA[:he4]
E = bhfp_jacc(input)

for k in (:he4, :he8, :he16, :he32, :he64, :he128, :he256)
    printstyled("Input: $(k)\n"; bold = true, color = :white)

    input = BasicHFProxy.DATA[k]
    t = @elapsed begin E = bhfp_threads_atomix(input) end

    println("E = ", E)
    println("Time: ", round(t; digits = 5), "sec\n")
    flush(stdout)
end
