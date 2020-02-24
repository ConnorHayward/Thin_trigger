

function run_sim()
    nthreads = Sys.CPU_THREADS
    ncores = Int(floor(nthreads / 2))

    run_sim(ncores)
end

function run_sim(ncores::Int64)
    cd("build/")
    run(`./PEN -m test.mac -t $ncores`)
    cd("..")
end
