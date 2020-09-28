# Sobol numbers in julia
# Sergio Ocampo
# September 2020
# Draw Sobol numbers and compare with uniform random
cd() # Go to root directory
cd("./Dropbox/Teaching/PhD_Macro_Comp/Julia_Code/Lecture_4_Optimization/")
mkpath("Figures/Sobol")
using Plots
using Sobol # Pkg.add("Sobol") # https://github.com/stevengj/Sobol.jl
using Random
println(" ")
println("------------------------")
println("Sobol numbers in Julia")
println("PWD: ",pwd())
println("This code uses Plots, Sobol, Random")
println("------------------------")
println(" ")


#-----------------------------------------------------------
#-----------------------------------------------------------
# Set random seed
Random.seed!(3486);

#-----------------------------------------------------------
#-----------------------------------------------------------
# One Dimension

    # Get Sobol sequence
    s = SobolSeq(1)
    s_vec = Float64[]
    # Get plots comparing numbers in one dimension
    for n=(10,20,50)
        # Get frist n elements of sobol sequence
        global s_vec = [s_vec ; hcat([next!(s) for i = 1:(n-length(s_vec))]...)']
        # Get n random numbers
        r_vec = rand(n,1)
        # Plot results
        gr()
        plt = plot([0;1],[1 2;1 2],linewidth=2,linecolor=RGB(0.6,0.6,0.6),label=nothing,title="Grid Points - n=$n")
        scatter!(r_vec,1*ones(n),marker=(:diamond,6),markercolor=1,label=nothing)
        scatter!(s_vec,2*ones(n),marker=(:star4,6)  ,markercolor=2,label=nothing)
        xlims!(0,1); ylims!(0,3); yticks!([1,2], ["Uniform", "Sobol"])
        savefig("./Figures/Sobol/Sobol_vs_Unifrom_1D_n$n.pdf")
    end

#-----------------------------------------------------------
#-----------------------------------------------------------
# Two Dimensions

    # Get Sobol sequence
    s = SobolSeq(2)
    s_vec = zeros(0,2)
    # Get plots comparing numbers in one dimension
    for n=(10,20,50,1024)
        # Get frist n elements of sobol sequence
        global s_vec = [s_vec ; hcat([next!(s) for i = 1:(n-length(s_vec[:,1]))]...)']
        # Plot results
        gr()
        plt1 = plot(title="Sobol Points - n=$n",legend=:outerright,foreground_color_legend = nothing,background_color_legend = nothing)
        s_scatter = scatter!(s_vec[:,1],s_vec[:,2],marker=(:star4,6),label=nothing,color=2)
        xlims!(0,1); ylims!(0,1);
        savefig("./Figures/Sobol/Sobol_2D_n$n.pdf")
        # Get n random numbers
        r_vec = rand(n,2)
        # Plot results
        gr()
        plt2 = plot(title="Uniform Points - n=$n",legend=:outerright,foreground_color_legend = nothing,background_color_legend = nothing)
        r_scatter = scatter!(r_vec[:,1],r_vec[:,2],marker=(:diamond,6),label=nothing,color=1)
        xlims!(0,1); ylims!(0,1);
        savefig("./Figures/Sobol/Unifrom_2D_n$n.pdf")
        # Joint Figure
        gr()
        plt2 = plot(legend=:outerbottom,foreground_color_legend = nothing,background_color_legend = nothing)
        r_scatter = scatter!(r_vec[:,1],r_vec[:,2],marker=(:diamond,6),label="Uniform",color=1)
        s_scatter = scatter!(s_vec[:,1],s_vec[:,2],marker=(:star4,6),label="Sobol",color=2)
        xlims!(0,1); ylims!(0,1);
        savefig("./Figures/Sobol/Sobol_vs_Unifrom_2D_n$n.pdf")
    end
