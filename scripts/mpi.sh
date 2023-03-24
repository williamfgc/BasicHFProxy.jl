julia --project -e 'using Pkg; Pkg.instantiate();'

for k in he4 he8 he16 he32 he64 he128 he256; do
    mpiexecjl --project -n 6 julia -t1 bhfp_mpi.jl $k
done
