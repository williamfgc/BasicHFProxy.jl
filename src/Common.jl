
import SpecialFunctions # erf

Base.@kwdef struct Parameters
    pi::Float64 = 3.141592653589793
    sqrpi2::Float64 = (pi^(-0.5)) * 2.0
    dtol::Float64 = 1.0e-10
    rcut::Float64 = 1.0e-12
    tobohrs::Float64 = 1.889725987722
end

function get_input_file(args::Array{String, 1})::String
    if length(args) != 1
        if MPI.Initialized()
            throw(ArgumentError("Usage: mpirun -n <nprocs> julia --project BasicHFProxy_<backend> <input_file>"))
        else
            throw(ArgumentError("Usage: julia --project BasicHFProxy_<backend> <input_file>"))
        end
    end
    return first(args)
end

"""
Main kernel for the calculation of the two-electron integrals
"""
function ssss(i::Int32, j::Int32, k::Int32, l::Int32, ngauss::Int32,
              xpnt::Array{Float64, 1}, coef::Array{Float64, 1},
              geom::Array{Float64, 2}, eri::Float64)
    eri = 0.0

    for ib in 1:ngauss
        for jb in 1:ngauss
            aij = 1.0 / (xpnt[ib] + xpnt[jb])
            dij = coef[ib] * coef[jb] *
                  exp(-xpnt[ib] * xpnt[jb] * aij *
                      ((geom[1, i] - geom[1, j])^2 +
                       (geom[2, i] - geom[2, j])^2 +
                       (geom[3, i] - geom[3, j])^2)) * (aij^1.5)

            if abs(dij) > parameters.dtol
                xij = aij * (xpnt[ib] * geom[1, i] + xpnt[jb] * geom[1, j])
                yij = aij * (xpnt[ib] * geom[2, i] + xpnt[jb] * geom[2, j])
                zij = aij * (xpnt[ib] * geom[3, i] + xpnt[jb] * geom[3, j])

                for kb in 1:ngauss
                    for lb in 1:ngauss
                        akl = 1.0 / (xpnt[kb] + xpnt[lb])
                        dkl = dij * coef[kb] * coef[lb] *
                              exp(-xpnt[kb] * xpnt[lb] * akl *
                                  ((geom[1, k] - geom[1, l])^2 +
                                   (geom[2, k] - geom[2, l])^2 +
                                   (geom[3, k] - geom[3, l])^2)) * (akl^1.5)

                        if abs(dkl) > parameters.dtol
                            aijkl = (xpnt[ib] + xpnt[jb]) *
                                    (xpnt[kb] + xpnt[lb]) /
                                    (xpnt[ib] + xpnt[jb] + xpnt[kb] + xpnt[lb])
                            tt = aijkl *
                                 ((xij -
                                   akl *
                                   (xpnt[kb] * geom[1, k] +
                                    xpnt[lb] * geom[1, l]))^2 +
                                  (yij -
                                   akl *
                                   (xpnt[kb] * geom[2, k] +
                                    xpnt[lb] * geom[2, l]))^2 +
                                  (zij -
                                   akl *
                                   (xpnt[kb] * geom[3, k] +
                                    xpnt[lb] * geom[3, l]))^2)

                            f0t = parameters.sqrpi2
                            if tt > parameters.rcut
                                f0t = (tt^(-0.5)) *
                                      SpecialFunctions.erf(sqrt(tt))
                            end #if
                            eri = eri + dkl * f0t * sqrt(aijkl)
                        end #if
                    end
                end
            end #if
        end
    end
end #ssss
