using BasicHFProxy

for s in (:he4, :he8, :he16, :he32, :he64, :he128, :he256)
    f = BasicHFProxy.DATA[s]
    E = bhfp_sequential(f)
    E_expected = BasicHFProxy.expected_energy(f)
    # printing
    printstyled("$s\n"; bold = true, color = :white)
    println("\t", "Energy: ", E)
    println("\t", "Exp. Energy: ", E_expected)
    println("\t", "Absolute diff: ", abs(E - E_expected))
    println("\t", "Relative diff: ",
            2 * abs(E - E_expected) / abs(E + E_expected))
    println()
end
