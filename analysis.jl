using CSV, DataFrames, StatsBase, Glob, Plots, LsqFit

@. model(x, p) = (1 / sqrt(2*pi*p[2]^2)) * exp(-x-p[2]^2/(2*p[2]^2))

Plots.pyplot()

foil_file = "/home/connor/Documents/Simulations/Thin_Trigger/build/6PMT_Trigger_Foil+Grease_3500_nt_DetectedPhotons.csv"
foil_data = CSV.read(foil_file, comment="#", header = ["N_Trigger","EDep_Trigger", "EDep_PEN", "N_Left", "N_Right", "N_Bottom", "N_Front", "N_Back"])

n_trigger_foil = length(foil_data.EDep_Trigger)
h_edep_foil = fit(Histogram, foil_data.EDep_Trigger, 0:10:1200)
trigger_plt = plot(h_edep_foil, st=:step, label = "Foil", xlabel = "Energy in Trigger [keV]")

h_edep_foil_pen = fit(Histogram, foil_data.EDep_PEN, 0:10:1200)
pen_plt = plot(h_edep_foil_pen, st=:step, label = "Foil", legend = :topleft, xlabel = "Energy in PEN [keV]")

foil_data.N_total = foil_data.N_Left .+ foil_data.N_Right .+ foil_data.N_Bottom .+ foil_data.N_Front .+ foil_data.N_Back
h_nphotons = fit(Histogram, foil_data.N_total, 0:10:2000)
peak_bin = findmax(h_nphotons.weights[30:end])[2]+29
println(peak_bin)
photon_plt = plot(h_nphotons, st=:step, label = "PEN", legend = :topright, xlabel = "Number of Photons", ylabel = "Events / 10 photons")

vline!(photon_plt, [midpoints(0:10:2000)[peak_bin]], label = "$(midpoints(0:10:2000)[peak_bin])")
combi_energy = plot(xlabel = "Energy [keV]", ylabel = "Events / 10 keV")
plot!(combi_energy, h_edep_foil_pen, st=:step, label = "PEN")
plot!(combi_energy, h_edep_foil, st=:step, label = "Trigger")
combined_plt = plot(layout = (2,1), size = (800,800), combi_energy, photon_plt)
savefig(combined_plt, "pen_sim_3500.pdf")
