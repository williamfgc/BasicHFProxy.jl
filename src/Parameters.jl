
Base.@kwdef struct Parameters
    pi::Float64 = 3.141592653589793
    sqrpi2::Float64 = (pi^(-0.5)) * 2.0
    dtol::Float64 = 1.0e-10
    rcut::Float64 = 1.0e-12
    tobohrs::Float64 = 1.889725987722
end
