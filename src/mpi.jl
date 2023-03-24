function bhfp_mpi(inputfile = get_input_filename_from_args(); verbose = false)
    MPI.Init()
    comm = MPI.COMM_WORLD
    myrank = MPI.Comm_rank(comm)
    nproc = MPI.Comm_size(comm)
    if myrank == 0
        println("number of ranks: ", nproc)
    end

    if myrank == 0
        (; ngauss, natom, xpnt, coef, geom) = parse_input_file(inputfile)
    else
        ngauss = -1
        natom = -1
    end

    ngauss = MPI.Bcast(ngauss, 0, comm)
    natom = MPI.Bcast(natom, 0, comm)

    if myrank != 0
        xpnt = Vector{Float64}(undef, ngauss)
        coef = Vector{Float64}(undef, ngauss)
        geom = Matrix{Float64}(undef, 3, natom)
    end

    MPI.Bcast!(xpnt, 0, comm)
    MPI.Bcast!(coef, 0, comm)
    MPI.Bcast!(geom, 0, comm)

    dens = Matrix{Float64}(undef, natom, natom)

    # fake ('random') density
    for i in 1:natom
        for j in 1:natom
            dens[i, j] = 0.1
        end
        dens[i, i] = 1.0
    end

    # Along with the density matrix, above, 3*(N^2) storage is typical
    # of real Hartree-Fock calculations.
    fock = Matrix{Float64}(undef, natom, natom)
    frcv = Matrix{Float64}(undef, natom, natom)

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

    # The following four-loop structure (with permutational symmetry)
    # represents the kernel of Hartree-Fock calculations.
    # Integrals are screened to avoid small terms.
    # At present, tasks are distributed amongst the MPI ranks inside
    # the second loop. This is done primarily for demonstration
    # purposes and does not guarantee optimal load balancing.
    task = 0
    ij = 0
    for i in 1:natom
        for j in 1:i
            task = task + 1
            ij = ij + 1
            if mod(task, nproc) == myrank
                kl = 0
                for k in 1:i
                    lx = k
                    if k == i
                        lx = j
                    end
                    for l in 1:lx
                        kl = kl + 1
                        if schwarz[ij] * schwarz[kl] > dtol
                            eri = ssss(i, j, k, l, ngauss, xpnt, coef, geom)
                            if i == j
                                eri = eri * 0.5
                            end
                            if k == l
                                eri = eri * 0.5
                            end
                            if i == k && j == l
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
    end

    # The following global sum of the Fock matrix is the only significant
    # communication overhead in replicated-data Hartree-Fock.
    # call mpi_allreduce( fock,frcv, natom**2, mpi_real8,mpi_sum,mpi_comm_world,errcon )
    MPI.Allreduce!(fock, frcv, MPI.SUM, comm)

    # trace Fock with the density, print the 2e- energy
    erep = 0.0
    for i in 1:natom
        for j in 1:natom
            erep = erep + frcv[i, j] * dens[i, j]
        end
    end
    E = erep * 0.5
    if verbose && myrank == 0
        println("2e- energy= ", E)
    end
    return E
end
