using FermiCG
using Plots
using JLD2
using LinearAlgebra

@load "/Users/nicole/code/FermiCG-data/he_clusters/07_oct/pes_tpsci.jl.545360.scr/scan_energies_triplet.jld2"

#display(plot([energies_ground.-energies_ground[end], energies_t1.-energies_t1[end], energies_t2.-energies_t2[end], energies_t3.-energies_t3[end], energies_t4.-energies_t4[end], energies_t5.-energies_t5[end], energies_t6.-energies_t6[end], energies_t7.-energies_t7[end]]*627.51, labels = ["ground" "T1" "T2" "T3" "T4" "T5" "T6" "T7"], xlabel="Reaction coordinates",ylabel="Energy (kcal/mol)")
#ylims!((-.5,.5))
