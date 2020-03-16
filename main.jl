include("control.jl")
include("analysis.jl")

abs_files = ["400nm-Fixed-20mm.csv", "450nm-Fixed-20mm.csv"]
light_yields = [3400, 6800]
n_events = 1000000

cd("build/")
for abs_file in abs_files
    for ly in light_yields
    script = "/PEN/det/setDetName $(abs_file[1:end-4])-LY-$ly
       /PEN/det/setAbsFile $abs_file
       /PEN/det/setSigAlpha 0.9
       /PEN/det/setLY $ly
       /run/beamOn $n_events"

       open("test.mac","w") do f
           write(f,script)
       end
        run(`./PEN -m test.mac -t 16`)
    end
end
