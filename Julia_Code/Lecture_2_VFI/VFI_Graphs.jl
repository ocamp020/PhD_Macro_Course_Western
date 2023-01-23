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


 
#-----------------------------------------------------------
#-----------------------------------------------------------
# Caution with Euler Errors 

    # Possible scenarios 
    x_vec = range(-1,1,1000)
    y_vec = x_vec.^3 
    gr() 
    plot([-1,1],[0,0],lw=1,linecolor=RGB(0.6,0.6,0.6),label=nothing)
    plot!(x_vec,y_vec,lw=2,label=nothing)
    ylims!(-1,1)
    xlabel!("Deviations from Policy")
    ylabel!("Euler Error (%)")
    savefig("./Figures/Euler_Error_Warning_1.pdf")
    y_vec = tan.(1.5*x_vec)
    gr() 
    plot([-1,1],[0,0],lw=1,linecolor=RGB(0.6,0.6,0.6),label=nothing)
    plot!(x_vec,y_vec,lw=2,label=nothing)
    ylims!(-1,1)
    xlabel!("Deviations from Policy")
    ylabel!("Euler Error (%)")
    savefig("./Figures/Euler_Error_Warning_2.pdf")
    y_vec = x_vec.^2
    gr() 
    plot([-1,1],[0,0],lw=1,linecolor=RGB(0.6,0.6,0.6),label=nothing)
    plot!(x_vec,y_vec,lw=2,label=nothing)
    ylims!(-1,1)
    xlabel!("Deviations from Policy")
    ylabel!("Euler Error (%)")
    savefig("./Figures/Euler_Error_Warning_3.pdf")
    y_vec = x_vec.^4
    gr() 
    plot([-1,1],[0,0],lw=1,linecolor=RGB(0.6,0.6,0.6),label=nothing)
    plot!(x_vec,y_vec,lw=2,label=nothing)
    ylims!(-1,1)
    xlabel!("Deviations from Policy")
    ylabel!("Euler Error (%)")
    savefig("./Figures/Euler_Error_Warning_4.pdf")

    # Low capital 
    k_0  = 1e-3 
    kp_0,c_0 = G_analytical(k_0,p)
    kp_vec = range(kp_0-1e-3,kp_0+1e-3,1000)
    kpp_vec, cp_vec = G_analytical(kp_vec,p)
    Euler_vec = [Euler_Error(k_0,kp_vec[i],kpp_vec[i],p) for i in 1:1000]
    gr()
    plot([-1e-3,1e-3],[0,0],lw=1,linecolor=RGB(0.6,0.6,0.6),label=nothing,title = "Sensitivity of Euler Error for Low K (%)")
    plot!(kp_vec.-kp_0,Euler_vec,lw=2,label=nothing)
    xlims!(-1e-3,1e-3)
    xlabel!("Deviations from Policy (in levels)")
    ylabel!("Euler Error (%)")
    savefig("./Figures/Euler_Error_Test_Low_K.pdf")
    # Steady State capital 
    k_0  = k_ss
    kp_0,c_0 = G_analytical(k_0,p)
    kp_vec = range(kp_0-1e-3,kp_0+1e-3,1000)
    kpp_vec, cp_vec = G_analytical(kp_vec,p)
    Euler_vec = [Euler_Error(k_0,kp_vec[i],kpp_vec[i],p) for i in 1:1000]
    gr()
    plot([-1e-3,1e-3],[0,0],lw=1,linecolor=RGB(0.6,0.6,0.6),label=nothing,title = "Sensitivity of Euler Error for SS K (%)")
    plot!(kp_vec.-kp_0,Euler_vec,lw=2,label=nothing)
    xlims!(-1e-3,1e-3)
    xlabel!("Deviations from Policy (in levels)")
    ylabel!("Euler Error (%)")
    savefig("./Figures/Euler_Error_Test_Med_K.pdf")
    # High capital 
    k_0  = 1.5*k_ss
    kp_0,c_0 = G_analytical(k_0,p)
    kp_vec = range(kp_0-1e-3,kp_0+1e-3,1000)
    kpp_vec, cp_vec = G_analytical(kp_vec,p)
    Euler_vec = [Euler_Error(k_0,kp_vec[i],kpp_vec[i],p) for i in 1:1000]
    gr()
    plot([-1e-3,1e-3],[0,0],lw=1,linecolor=RGB(0.6,0.6,0.6),label=nothing,title = "Sensitivity of Euler Error for High K (%)")
    plot!(kp_vec.-kp_0,Euler_vec,lw=2,label=nothing)
    xlims!(-1e-3,1e-3)
    xlabel!("Deviations from Policy (in levels)")
    ylabel!("Euler Error (%)")
    savefig("./Figures/Euler_Error_Test_High_K.pdf")


println("\n     Graphs Completed for VFI_$VFI_Type \n")
