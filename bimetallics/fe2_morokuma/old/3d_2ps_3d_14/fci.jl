using QCBase
using ActiveSpaceSolvers
using InCoreIntegrals
using LinearAlgebra
using Printf
using Test
using Arpack
using NPZ


@load "data_cmf.jld2"

n_elec_a = 9
n_elec_b = 9

norb = size(h1,1)
ansatz = FCIAnsatz(norb, n_elec_a, n_elec_b)

display(ansatz)

solver = SolverSettings(nroots=6, package="arpack")
println(solver)
solution = solve(ints, ansatz, solver)
display(solution)
