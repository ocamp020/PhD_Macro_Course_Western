# VFI Program for Julia
# Sergio Ocampo
# September 2020
# NGM with log utility and full depreciation
# Solve:   V(k) = max{ log(zk^alpha-k') +beta*V(k') }
# Solution by grid search: k_grid = {k_1,...,k_I}
# Graphing
gr()
# Value Function Analytical vs 200
    plot(k_grid_200,V_200a,linewidth=3,label = "Analytical Solution",title="Value Function",legend=(0.75,0.2),foreground_color_legend = nothing,background_color_legend = nothing)
    plot!(k_grid_50,V_50,linetype=:scatter,marker=(:diamond,6),markercolor=RGB(0.5,0.1,0.1),label="VFI - n_k=50")
    plot!(k_grid_200,V_200,linetype=:scatter,marker=(:diamond,3),markercolor=RGB(0.1,0.1,0.1),label = "VFI - n_k=200")
    xlabel!("Capital")
    ylabel!("Value")
    savefig("./Figures/VFI_"*VFI_Type*"_V.pdf")
# Capital Policy Function Analytical vs 200
    plot(k_grid_200,k_grid_200,lw=1,linecolor=RGB(0.6,0.6,0.6),label=nothing)
    plot!(k_grid_200,G_kp_200a,linewidth=3,label = "Analytical Solution",title="Policy Function - K",legend=(0.75,0.2),foreground_color_legend = nothing,background_color_legend = nothing)
    plot!(k_grid_50,k_grid_50[G_kp_50],linetype=:scatter,marker=(:diamond,6),markercolor=RGB(0.5,0.1,0.1),label="VFI - n_k=50")
    plot!(k_grid_200,k_grid_200[G_kp_200],linetype=:scatter,marker=(:diamond,3),markercolor=RGB(0.1,0.1,0.1),label = "VFI - n_k=200")
    xlabel!("Capital")
    ylabel!("Capital")
savefig("./Figures/VFI_"*VFI_Type*"_G_kp.pdf")
    # Euler Error 200
    plot([0,2*k_ss],[0,0],lw=1,linecolor=RGB(0.6,0.6,0.6),label=nothing,title = "Euler Equation Error (%)")
    plot!(k_grid_50,Euler_50,linetype=:scatter,marker=(:diamond,6),markercolor=RGB(0.5,0.1,0.1),label="VFI - n_k=50")
    plot!(k_grid_200,Euler_200,linetype=:scatter,marker=(:diamond,3),markercolor=RGB(0.1,0.1,0.1),label="VFI - n_k=200")
    xlabel!("Capital")
    ylabel!("Percentage Points")
savefig("./Figures/VFI_"*VFI_Type*"_Euler.pdf")

println("\n     Graphs Completed for VFI_$VFI_Type \n")
