# Optimization and Root Finding for Julia
# Sergio Ocampo
# September 2020
cd() # Go to root directory
cd("./Dropbox/Teaching/PhD_Macro_Comp/Julia_Code/Lecture_4_Optimization/")
mkpath("Figures")
using Plots
# using LateXStrings # Pkg.add("LaTeXStrings") # https://github.com/stevengj/LaTeXStrings.jl
using Dierckx # Pkg.add("Dierckx") # https://github.com/kbarbary/Dierckx.jl
using Interpolations # Pkg.add("Interpolations") # https://github.com/JuliaMath/Interpolations.jl
using ForwardDiff # Pkg.add("ForwardDiff") # https://github.com/JuliaDiff/ForwardDiff.jl

println(" ")
println("------------------------")
println("Interpolation in Julia")
println("PWD: ",pwd())
println("This code uses Plots, Interpolations, Dierckx and ForwardDiff")
println("------------------------")
println(" ")


#-----------------------------------------------------------
#-----------------------------------------------------------
# Runge's Function
