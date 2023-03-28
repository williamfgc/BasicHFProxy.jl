# constants
const pi = 3.141592653589793
const sqrpi2 = (pi^(-0.5)) * 2.0
const dtol = 1.0e-10
const rcut = 1.0e-12
const tobohrs = 1.889725987722

# input file
function get_input_filename_from_args(args = ARGS)
    ids = collect(keys(BasicHFProxy.DATA))
    if length(args) != 1
        throw(ArgumentError("Please provide an identifier in `ARGS` or an explicit input file." *
                            " Supported identifiers: $(join(ids, ", "))"))
    end
    id = Symbol(only(args))
    if haskey(BasicHFProxy.DATA, id)
        return BasicHFProxy.DATA[id]
    else
        throw(ArgumentError("Unknown input identifier \"$id\"." *
                            " Supported identifiers: $(collect(keys(BasicHFProxy.DATA)))"))
    end
end

"""
Parses given input file and returns a `NamedTuple` with the following keys:
* `ngauss`: number of gaussian-type functions (GTF) per atom
* `natom`: number of atoms
* `xpnt`: vector of expansion exponents
* `coef`: vector of expansion coefficients
* `geom`: geometry matrix
"""
function parse_input_file(filename)
    # read + split entire file
    data = split(read(filename, String))
    # parse first two records as integers
    ngauss = parse(Int, popfirst!(data))
    natom = parse(Int, popfirst!(data))
    @debug("Input data", ngauss, natom)

    # parse next records as exponents and coefficients
    xpnt = Vector{Float64}(undef, ngauss)
    coef = Vector{Float64}(undef, ngauss)
    for i in 1:ngauss
        xpnt[i] = parse(Float64, popfirst!(data))
        coef[i] = parse(Float64, popfirst!(data))
    end
    @debug("exponents", xpnt)
    @debug("coefficients", coef)

    # parse geometry matrix
    geom = Matrix{Float64}(undef, 3, natom)
    for i in 1:natom
        for j in 1:3
            geom[j, i] = parse(Float64, popfirst!(data))
        end
    end
    return (; ngauss, natom, xpnt, coef, geom)
end

"""
Extracts the expected energy from the given input file.
"""
function expected_energy(filename)
    energy_line = only(filter(contains("2e- energy="), readlines(filename)))
    energy_str = strip(replace(energy_line, "2e- energy=" => ""))
    return parse(Float64, energy_str)
end

# computation
"""
Main kernel for the calculation of the two-electron integrals
"""
function ssss(i::Integer, j::Integer, k::Integer, l::Integer, ngauss::Integer,
              xpnt::AbstractVector{<:Float64}, coef::AbstractVector{<:Float64},
              geom::AbstractMatrix{<:Float64})
    eri = 0.0
    for ib in 1:ngauss
        for jb in 1:ngauss
            aij = 1.0 / (xpnt[ib] + xpnt[jb])
            dij = coef[ib] * coef[jb] *
                  exp(-xpnt[ib] * xpnt[jb] * aij *
                      ((geom[1, i] - geom[1, j])^2 +
                       (geom[2, i] - geom[2, j])^2 +
                       (geom[3, i] - geom[3, j])^2)) * (aij^1.5)

            if abs(dij) > dtol
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

                        if abs(dkl) > dtol
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

                            f0t = sqrpi2
                            if tt > rcut
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
    return eri
end #ssss
