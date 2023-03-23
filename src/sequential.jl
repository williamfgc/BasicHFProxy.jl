function bhfp_sequential(inputfile = get_input_filename_from_args(); verbose=false)
    (; ngauss, natom, xpnt, coef, geom) = parse_input_file(inputfile)

    dens = Matrix{Float64}(undef, natom, natom)

    # fake ('random') density
    for i in 1:natom
        for j in 1:natom
            dens[i, j] = 0.1
        end
        dens[i, i] = 1.0
    end

    # normalize the primitive GTO weights
    for i in 1:ngauss
        coef[i] = coef[i] * (2.0 * xpnt[i])^0.75
    end

    # scale the geometry to Bohrs for energy calculations in AU
    for i in 1:natom
        geom[1, i] = geom[1, i] * tobohrs
        geom[2, i] = geom[2, i] * tobohrs
        geom[3, i] = geom[3, i] * tobohrs
    end

    fock = Matrix{Float64}(undef, natom, natom)
    for i in 1:natom
        for j in 1:natom
            fock[i, j] = 0.0
        end
    end

    # compute Schwarz Inequality factors for integral screening
    schwarz = Vector{Float64}(undef, ((natom^2) + natom) ÷ 2)

    ij = 0
    for i in 1:natom
        for j in 1:i
            ij = ij + 1
            eri = ssss(i, j, i, j, ngauss, xpnt, coef, geom)
            schwarz[ij] = sqrt(abs(eri))
        end
    end

    # this loop-structure reflects the 8-fold label symmetry of the integrals

    ij = 0
    for i in 1:natom
        for j in 1:i
            ij = ij + 1
            kl = 0
            for k in 1:i
                lx = k
                if (k == i)
                    lx = j
                end
                for l in 1:lx
                    kl = kl + 1
                    if (schwarz[ij] * schwarz[kl] > dtol)
                        eri = ssss(i, j, k, l, ngauss, xpnt, coef, geom)   # main kernel
                        if (i == j)
                            eri = eri * 0.5
                        end
                        if (k == l)
                            eri = eri * 0.5
                        end
                        if (i == k && j == l)
                            eri = eri * 0.5
                        end
                        fock[i, j] = fock[i, j] + dens[k, l] * eri * 4.0
                        fock[k, l] = fock[k, l] + dens[i, j] * eri * 4.0
                        fock[i, k] = fock[i, k] - dens[j, l] * eri
                        fock[i, l] = fock[i, l] - dens[j, k] * eri
                        fock[j, k] = fock[j, k] - dens[i, l] * eri
                        fock[j, l] = fock[j, l] - dens[i, k] * eri
                    end
                end
            end
        end
    end

    # trace Fock with the density, print the 2e- energy
    erep = 0.0
    for i in 1:natom
        for j in 1:natom
            erep = erep + fock[i, j] * dens[i, j]
        end
    end
    E = erep * 0.5
    verbose && println("2e- energy= ", E)
    return E
end
