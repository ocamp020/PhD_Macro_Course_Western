# EGM in Julia
# Sergio Ocampo
# September 2020
# EGM on the NGM with CRRA utility and AR(1) productivity
# Solve:   V(z,k) = max{ (zk^alpha-k')^(1-γ)/(1-γ) +beta*E[V(z',k')|z] }
#           log(z') = ρ*log(z) + η; η~N(0,σ)
cd() # Go to root directory
cd("./Dropbox/Teaching/PhD_Macro_Comp/Julia_Code/Lecture_6_EGM/")
mkpath("Figures")
using Plots
using Interpolations # Pkg.add("Interpolations") # https://github.com/JuliaMath/Interpolations.jl
using Dierckx # Pkg.add("Dierckx") # https://github.com/kbarbary/Dierckx.jl
using ForwardDiff # Pkg.add("ForwardDiff") # https://github.com/JuliaDiff/ForwardDiff.jl
using Optim # Pkg.add("Optim") # https://julianlsolvers.github.io/Optim.jl/stable/
    using Optim: converged, maximum, maximizer, minimizer, iterations
using Roots # Pkg.add("Roots") # https://github.com/JuliaMath/Roots.jl
using Parameters # Pkg.add("Parameters") # https://github.com/mauro3/Parameters.jl
include("../VFI_Toolbox.jl")
# using .VFI_Toolbox
println(" ")
println("------------------------")
println("EGM in Julia")
println("PWD: ",pwd())
println("This code uses Plots, Interpolations, Dierckx, ForwardDiff, Optim, Roots, Parameters, ScaledInterpolation")
println("EGM in the context of the neoclassical growth model")
println("------------------------")
println(" ")



#-----------------------------------------------------------
#-----------------------------------------------------------
# Paramters and Model Structure
    # Generate structure for parameters using Parameters module
    # We can set default values for our parameters
    @with_kw struct Par
        # Model Parameters
        α::Float64 = 1/3  ; # Production function
        β::Float64 = 0.96 ; # Discount factor
        γ::Float64 = 2.0  ; # Relative risk aversion parameter
        δ::Float64 = 0.05 ; # Depreciation rate
        ρ::Float64 = 0.90 ; # Persistence of productivity process
        σ::Float64 = 0.05 ; # Standard devaiation of productivity innovation
        z_bar::Float64 = 1; # Reference level for productivity
        # VFI Paramters
        max_iter::Int64   = 2000  ; # Maximum number of iterations
        dist_tol::Float64 = 1E-12 ; # Tolerance for distance
        # Howard's Policy Iterations
        H_tol::Float64    = 1E-9  ; # Tolerance for policy function iteration
        N_H::Int64        = 20    ; # Maximum number of policy iterations
        # Minimum consumption for numerical optimization
        c_min::Float64    = 1E-16
    end

    # Allocate paramters to object p for future calling
    p = Par()


    # Generate structure of model objects
    @with_kw struct Model
        # Parameters
        p::Par = Par() # Model paramters in their own structure
        # Steady State Values
        k_ss = (p.β*p.α*p.z_bar/(1-p.β*(1-p.δ)))^(1/(1-p.α))
        # Capital Grid
        θ_k::Float64    = 1.5                        # Curvature of k_grid
        n_k::Int64      = 500                       # Size of k_grid
        n_k_fine::Int64 = 1000                      # Size of fine grid for interpolation
        k_grid          = Make_Grid(n_k     ,θ_k,1E-5,2*k_ss)  # k_grid for model solution
        k_grid_fine     = Make_Grid(n_k_fine,1  ,1E-5,2*k_ss)  # Fine grid for interpolation
        # Productivity process
        n_z       = 5                               # Size of z_grid
        MP_z      = Rouwenhorst95(p.ρ,p.σ,n_z)      # Markov Process for z
        # State matrices
        k_mat     = repeat(k_grid',n_z,1)
        z_mat     = exp.(repeat(MP_z.grid,1,n_k))
        # Y_grid and Marginal Product of Capital
        Y_grid    =     p.z_bar*z_mat.*(k_mat.^(p.α)  ) .+ (1-p.δ).*k_mat
        MPk_mat   = p.α*p.z_bar*z_mat.*(k_mat.^(p.α-1)) .+ (1-p.δ)
        # Value and policy functions
        V         = Array{Float64}(undef,n_z,n_k)       # Value Function
        G_kp      = Array{Float64}(undef,n_z,n_k)       # Policy Function
        G_c       = Array{Float64}(undef,n_z,n_k)       # Policy Function
        V_fine    = Array{Float64}(undef,n_z,n_k_fine)  # Value Function on fine grid
        G_kp_fine = Array{Float64}(undef,n_z,n_k_fine)  # Policy Function on fine grid
        G_c_fine  = Array{Float64}(undef,n_z,n_k_fine)  # Policy Function on fine grid
        # Error in Euler equation
        Euler     = Array{Float64}(undef,n_z,n_k_fine)  # Errors in Euler equation
    end

    # Allocate model to object M for future calling
    M = Model()


#-----------------------------------------------------------
#-----------------------------------------------------------
# Utility function
function utility(z,k,kp,p::Par)
    @unpack z_bar, α, γ, δ, c_min = p
    c = max.(z_bar.*z.*k.^α .+ (1-δ).*k  .- kp,c_min)
    if γ>1
    return (c).^(1-γ)/(1-γ)
    else
    return log.(c)
    end
end

function utility(c,p::Par)
    if p.γ>1
    return (c).^(1-p.γ)/(1-p.γ)
    else
    return log.(c)
    end
end

function d_utility(z,k,kp,p::Par)
    @unpack z_bar, α, γ, δ, c_min = p
    c = max.(z_bar.*z.*k.^α .+ (1-δ).*k  .- kp,c_min)
    return (c).^(-γ)
end

function d_utility(c,p::Par)
    return (c).^(-p.γ)
end

function d_utility_inv(x,p::Par)
    return x.^(-1/p.γ)
end

#-----------------------------------------------------------
#-----------------------------------------------------------
# Steady state values (funciton)
function SS_values(p::Par)
    # This function takes in parameters and provides steady state values
    # Parameters: productivity (z), returns to scale (a) and discount factor (b)
    # Output: values for capital, production, consumption, rental rate, wage
    @unpack z_bar, α, β, δ = p
    k_ss = (β*α*z_bar/(1-β*(1-δ)))^(1/(1-α))
    y_ss = z_bar*k_ss^α
    c_ss = y_ss - δ*k_ss
    r_ss = α*y_ss/k_ss
    w_ss = (1-α)*y_ss
    return k_ss,y_ss,c_ss,r_ss,w_ss
end

    # Test steady state function
    k_ss,y_ss,c_ss,r_ss,w_ss = SS_values(p)
    println(" ")
    println("------------------------")
    println(" Steady State variables")
    println("   Quantities: k = $k_ss; y = $y_ss; c = $c_ss;")
    println("   Prices:     r = $r_ss; w = $w_ss;")
    println("------------------------")
    println(" ")

# Value function for non-stochastic case with full depreciation and log-utility
    function V_analytical(k,p::Par)
        # This function takes in paramters and a value of capital
        # Output is the value of the value function at that level of capital
        # Parameters: productivity (z), returns to scale (a) and discount factor (b)
        @unpack z_bar, α, β = p
        a_1 = α/(1-α*β)
        a_0 = (1/(1-β))*(log(z_bar)+log(1/(1+β*a_1))+β*a_1*log(β*a_1*z_bar/(1+β*a_1*z)))
        V   = a_0 .+ a_1*log.(k)
        return V
    end

#-----------------------------------------------------------
#-----------------------------------------------------------
# Euler Error Function

# G_k interpolation
function G_kp_zk(i_z::Int64,k,M::Model)
    itp = ScaledInterpolations(M.k_grid,M.G_kp[i_z,:], BSpline(Cubic(Line(OnGrid()))))
    return itp(k)
end

# Euler Equation
function Euler_Error(i_z::Int64,k,M::Model)
    # Return percentage error in Euler equation
    @unpack p, MP_z, n_z = M
    @unpack z_bar, α, β, δ = p
    # Iterpolate G_k at current z
    kp  = min(M.k_grid[end],max(M.k_grid[1],G_kp_zk(i_z,k,M)))

    # Compute left hand side of Euler equation
    LHS = d_utility(exp(MP_z.grid[i_z]),k,kp,p)
    # Compute right hand side of Euler equation
        # Marginal product at kp for all z
        MPkp = α*z_bar*exp.(MP_z.grid).*kp^(α-1) .+ (1-δ)
        # Marginal utility at z',k',G_k(z',k')
        up   = zeros(n_z)
        for i_zp=1:n_z
        up[i_zp] = d_utility(exp(MP_z.grid[i_zp]),kp,G_kp_zk(i_zp,kp,M),p)
        end
    RHS = β*(MP_z.Π[i_z,:])'*(MPkp.*up)
    # Return percentage errror in Euler equation
    return (RHS/LHS-1)*100
end

# Euler Equation for Optimization
function Euler_Eq(kp,i_z::Int64,i_k::Int64,Vk,M::Model)
    # Return percentage error in Euler equation
    @unpack p, MP_z, n_z, k_grid = M
    @unpack β = p
    # Compute left hand side of Euler equation
    LHS = d_utility(exp(MP_z.grid[i_z]),k_grid[i_k],kp,p)
    # Compute right hand side of Euler equation
        # Expected value of derivative β*E[Vk[z',kp]|k]
        Vkp  = zeros(n_z)
        for i_zp=1:n_z
            Vk_ip     = ScaledInterpolations(k_grid,Vk[i_zp,:], FritschButlandMonotonicInterpolation())
            Vkp[i_zp] = Vk_ip(kp)
        end
    RHS = β*(MP_z.Π[i_z,:])'*Vkp
    # Return squared errror in Euler equation
    return (RHS/LHS-1)^2
end

#-----------------------------------------------------------
#-----------------------------------------------------------
# VFI Fixed Point
function VFI_Fixed_Point(T::Function,M::Model,V_old=nothing)
    # Unpack model structure
    @unpack p, n_z, n_k, θ_k, k_grid, n_k_fine, k_grid_fine = M
    # VFI paramters
    @unpack max_iter, dist_tol = p
    # Initialize variables for loop
    if V_old==nothing
    # V_old  = zeros(n_z,n_k)     ; # Initialize value function, here I just do 0, not the best option
    V_old  = utility(M.z_mat,M.k_mat,zeros(n_z,n_k),p)  ; # Start at utility with zero savings
    # V_old  = utility((1-p.δ)*(M.k_mat),p)  ; # Start at utility with c=(1-δ)k
    end
    # println("V_0="); display(V_old)
    G_kp_old = copy(M.k_mat)
    V_dist = 1              ; # Initialize distance
    println(" ")
    println("------------------------")
    println("VFI - n_z=$n_z, n_k=$n_k - θ_k=$θ_k")
    for iter=1:max_iter
        # println("\n Iteration $iter")
        # Update value function
        V_new, G_kp, G_c = T(Model(M,V=copy(V_old)))
            # println("T(V) = $V_new")
            # println("  V  = $V_old")
        # Update distance and iterations
        V_dist = maximum(abs.(V_new./V_old.-1))
        # V_dist = maximum(abs.(G_kp./G_kp_old.-1))
        # Update old function
        V_old  = V_new
        # G_kp_old = copy(G_kp)
        # Report progress
        if mod(iter,100)==0
            println("   VFI Loop: iter=$iter, dist=",100*V_dist,"%")
        end
        # Check convergence and return results
        if V_dist<=dist_tol
            println("VFI - n_z=$n_z, n_k=$n_k - θ_k=$θ_k")
            println("Iterations = $iter and Distance = ",100*V_dist,"%")
            println("------------------------")
            println(" ")
            # Interpolate to fine grid
            V_fine    = zeros(n_z,n_k_fine)
            G_kp_fine = zeros(n_z,n_k_fine)
            G_c_fine  = zeros(n_z,n_k_fine)
            for i_z=1:n_z
            V_ip    = ScaledInterpolations(k_grid,V_new[i_z,:], BSpline(Cubic(Line(OnGrid()))))
                V_fine[i_z,:]   .= V_ip.(collect(k_grid_fine))
            G_kp_ip = ScaledInterpolations(k_grid,G_kp[i_z,:] , BSpline(Cubic(Line(OnGrid()))))
                G_kp_fine[i_z,:].= G_kp_ip.(collect(k_grid_fine))
            G_c_ip  = ScaledInterpolations(k_grid,G_c[i_z,:]  , BSpline(Cubic(Line(OnGrid()))))
                G_c_fine[i_z,:] .= G_c_ip.(collect(k_grid_fine))
            end
            # Update model
            M = Model(M; V=V_new,G_kp=G_kp,G_c=G_c,V_fine=V_fine,G_kp_fine=G_kp_fine,G_c_fine=G_c_fine)
            # Euler Equation Errors
            Euler = [Euler_Error(i_z,k_grid_fine[i_k],M) for i_z=1:n_z, i_k in 1:n_k_fine]
            # Update model
            M = Model(M; Euler=Euler)
            # Return results
            return M
        end
    end
    # If loop ends there was no convergence -> Error!
    error("Error in VFI - Solution not found")
end


#-----------------------------------------------------------
#-----------------------------------------------------------
# Graphs
function VFI_Graphs(M::Model,VFI_Type)
    gr()
    # Value Function
        plt = plot(title="Value Function - n_k=$(M.n_k) - θ_k=$(M.θ_k)",legend=:bottomright,foreground_color_legend = nothing,background_color_legend = nothing)
        for i_z=1:M.n_z
        plot!(M.k_grid,M.V[i_z,:],linetype=:scatter,label="V(z_$i_z)")
        plot!(M.k_grid_fine,M.V_fine[i_z,:],linewidth=2.5,linestyle=(:dash),linecolor=RGB(0.4,0.4,0.4),label=nothing)
        end
        xlabel!("Capital")
        ylabel!("Value")
        savefig("./Figures/VFI_"*VFI_Type*"_V_$(M.n_k)_$(M.θ_k).pdf")
    # Capital Policy Function Analytical vs 200
        plt = plot(title="Policy Function - K - n_k=$(M.n_k) - θ_k=$(M.θ_k)",legend=:bottomright,foreground_color_legend = nothing,background_color_legend = nothing)
        plot!(M.k_grid_fine,M.k_grid_fine,lw=1,linecolor=RGB(0.6,0.6,0.6),label=nothing)
        for i_z=1:M.n_z
        plot!(M.k_grid,M.G_kp[i_z,:],linetype=:scatter,label="G_kp(z_$i_z)")
        plot!(M.k_grid_fine,M.G_kp_fine[i_z,:],linewidth=2.5,linestyle=(:dash),linecolor=RGB(0.4,0.4,0.4),label=nothing)
        end
        xlabel!("Capital")
        ylabel!("Capital")
    savefig("./Figures/VFI_"*VFI_Type*"_G_kp_$(M.n_k)_$(M.θ_k).pdf")
        # Euler Error 200
        plt = plot(title="Euler Equation Error (%) - n_k=$(M.n_k) - θ_k=$(M.θ_k)",legend=:bottomright,foreground_color_legend = nothing,background_color_legend = nothing)
        plot!(M.k_grid_fine,zeros(M.n_k_fine),lw=1,linecolor=RGB(0.6,0.6,0.6))
        for i_z=1:M.n_z
        plot!(M.k_grid_fine,M.Euler[i_z,:],linetype=:scatter,label="z_$i_z")
        end
        xlabel!("Capital")
        ylabel!("Percentage Points")
    savefig("./Figures/VFI_"*VFI_Type*"_Euler_$(M.n_k)_$(M.θ_k).pdf")

    println("\n     Graphs Completed for VFI_$VFI_Type - n_k=$(M.n_k) - θ_k=$(M.θ_k)\n")
end


#-----------------------------------------------------------
#-----------------------------------------------------------
# Bellman operator - EGM
function T_EGM(M::Model)
    @unpack p, n_z, MP_z, n_k, k_grid, Y_grid, V, G_kp, G_c = M
    @unpack β = p
    # println("Initial V = $V \n ")
    # Check Monotonicity
    if any( diff(V,dims=2).<0 )
        error("V must be monotone for EGM to work")
    end
    # Define expectation of value function for each (z,k')
    EV = β*MP_z.Π*V  # Rows are present z and columns are tomorrow's k in fixed grid
    # println("Initial EV size=$(size(EV)), and EV=$EV")
    # Check Monotonicity
    if any( diff(EV,dims=2).<0 )
        error("EV must be monotone for EGM to work")
    end
    # Define the derivative of EV for each value
        # I will be lazy and interpolate to then get the derivative with ForwardDiff
    dEV = zeros(size(EV))
    for i_z=1:n_z
        # EV_ip      = ScaledInterpolations(k_grid,EV[i_z,:], BSpline(Cubic(Line(OnGrid()))))
        # dEV_ip(x)  = ForwardDiff.derivative(EV_ip,x)
        # dEV[i_z,:].= dEV_ip.(k_grid)
        # EV_ip      = Spline1D(k_grid,EV[i_z,:])
        # dEV[i_z,:].= derivative(EV_ip, k_grid )
        EV_ip      = ScaledInterpolations(k_grid,EV[i_z,:], FritschButlandMonotonicInterpolation())
        dEV_ip(x)  = ForwardDiff.derivative(EV_ip,x)
        dEV[i_z,:].= dEV_ip.(k_grid)
        # println("\n i_z=$i_z, dEV=$(dEV[i_z,:]) \n")
    end
    # Check Monotonicity
    if any( dEV.<0 )
        for i_z=1:n_z
            println("\n i_z=$i_z, [min,max]=$([minimum(dEV[i_z,:]),maximum(dEV[i_z,:])]) \n")
        end
        error("dEV must be monotone for EGM to work")
    end
        # gr()
        # plot(title="EV Test")
        # plot!(k_grid,EV')
        # plot(title="dEV Test")
        # plot!(k_grid,dEV')
    # Define Consumption from Euler Equation
    C_endo = d_utility_inv(dEV,p)
    # Define endogenous grid on cash on hand
    Y_endo = C_endo .+ M.k_mat
        # Sort Y_endo for interpolation
        for i_z=1:n_z
        sort_ind = sortperm(Y_endo[i_z,:])
        Y_endo[i_z,:] .= Y_endo[i_z,:][sort_ind]
        C_endo[i_z,:] .= C_endo[i_z,:][sort_ind]
        end
    # Define value function on endogenous grid
    V_endo = utility(C_endo,p) .+ EV
    # Interpolate functions on exogenous grid
    for i_z=1:n_z
        # V_ip        = ScaledInterpolations(Y_endo[i_z,:],V_endo[i_z,:], BSpline(Cubic(Line(OnGrid()))))
        V_ip        = Spline1D(Y_endo[i_z,:],V_endo[i_z,:])
        V[i_z,:]   .= V_ip.(Y_grid[i_z,:])
        # C_ip        = ScaledInterpolations(Y_endo[i_z,:],C_endo[i_z,:], BSpline(Cubic(Line(OnGrid()))))
        C_ip        = Spline1D(Y_endo[i_z,:],C_endo[i_z,:])
        G_c[i_z,:] .= C_ip.(Y_grid[i_z,:])
        G_kp[i_z,:].= Y_grid[i_z,:] .- G_c[i_z,:]
    end
    # Solve manually if capital is out of bounds
    for ind = findall(<=(k_grid[1]),G_kp)
        # Minimize square of Euler Error
        k_max      = Y_grid[ind] - p.c_min
        min_result = optimize(x->Euler_Eq(x,ind[1],ind[2],Vk,M),k_grid[1],k_max)
        # Check result
        converged(min_result) || error("Failed to solve Euler Equation in $(iterations(min_result)) iterations")
        # Upddate policy function
        G_kp[ind]  = min_result.minimizer
        G_c[ind]   = Y_grid[ind] - G_kp[ind]
        EV_ip      = ScaledInterpolations(k_grid,EV[ind[1],:], FritschButlandMonotonicInterpolation())
        V[ind]     = utility(G_c[ind],p) + EV_ip(G_kp[ind])
    end
    # Return Results
        # println("T(V) = $V")
        # stop
    return V, G_kp, G_c
end




#-----------------------------------------------------------
#-----------------------------------------------------------
# Bellman operator - ECM
function T_ECM(M::Model)
    @unpack p, n_z, MP_z, n_k, k_grid, Y_grid, V, G_kp, G_c, k_mat, z_mat, MPk_mat = M
    @unpack β, α, δ, z_bar = p
    # println("Initial V = \n "); display(V)
    # Option 1: Get Vk from interpolation
        Vk = zeros(n_z,n_k)
        for i_z=1:n_z
            V_ip      = ScaledInterpolations(k_grid,V[i_z,:], FritschButlandMonotonicInterpolation())
            Vk_ip(x)  = ForwardDiff.derivative(V_ip,x)
            Vk[i_z,:].= Vk_ip.(k_grid)
            # V_ip      = Spline1D(k_grid,V[i_z,:])
            # Vk[i_z,:].= derivative(V_ip, k_grid )
            # println("\n i_z=$i_z, Vk=$(Vk[i_z,:]) \n")
            # V_ip, dV_ip = spline_NR(collect(k_grid),V[i_z,:])
            # Vk[i_z,:]  .= dV_ip.(collect(k_grid))

        end
    # Optioin 2: Get Vk from V (the program actually receives Vk)
        # Vk = V
    # Check Monotonicity
        if any( Vk.<0 )
            for i_z=1:n_z
                println("\n i_z=$i_z, [min,max]=$([minimum(Vk[i_z,:]),maximum(Vk[i_z,:])]) \n")
            end
            error("Vk must be positive for ECM to work")
        end
    # # Check for maximum allowable consumption
    #     Vk .= max.(d_utility(Y_grid.-k_grid[1],p).*MPk_mat,Vk)
    # Define consumption from envelope condition
    G_c  = d_utility_inv(Vk./MPk_mat,p)
    # Define savings from budget  constraint
    G_kp.= Y_grid .- G_c
        # Solve manually if capital is out of bounds
        for ind = findall(<=(k_grid[1]),G_kp)
            # Minimize square of Euler Error
            k_max      = Y_grid[ind] - p.c_min
            min_result = optimize(x->Euler_Eq(x,ind[1],ind[2],Vk,M),k_grid[1],k_max)
            # Check result
            converged(min_result) || error("Failed to solve Euler Equation in $(iterations(min_result)) iterations")
            # Upddate policy function
            G_kp[ind]  = min_result.minimizer
        end
        # Check for non-negativity (actually bounded by min grid point to avoid extrapolations)
        G_kp.= min.(max.(k_grid[1],G_kp),k_grid[end])
        G_c .= Y_grid .- G_kp
        # Vk    = d_utility(G_c,p).*MPk_mat
    # Define expectation of value function for each (z,k')
    EV = zeros(n_z,n_k) # E[V(z',k'(z,k))|z]
    for i_z=1:n_z
        kp = G_kp[i_z,:]
        Vp = zeros(n_z,n_k) # rows are zp and columns are k'(z,k) varying over k, z is fixed
        for i_zp=1:n_z
            # Define interpolated value at (z',G_kp(z,k))
            V_ip       = ScaledInterpolations(k_grid,V[i_zp,:], FritschButlandMonotonicInterpolation())
            Vp[i_zp,:].= V_ip.(kp)
        end
        EV[i_z,:] = MP_z.Π[i_z,:]'*Vp # Rows are present z and columns are future k as in k'=G_kp(z,k)
    end
    # Update value
        # Option 1: Update value function
        V.= utility(G_c,p) .+ β*EV
        # Option 2: Update derivative of value function
        # V = β*MPk_mat.*EV
    # Return Results
        # println("T(V) = "); display(V)
        # stop
    return V, G_kp, G_c
end





#-----------------------------------------------------------
#-----------------------------------------------------------
# Execute VFI and plot graphs

# Execute Numerical VFI - EGM
    println("===============================================\n Solving VFI with EGM")
    @time M_EGM  = VFI_Fixed_Point(T_EGM,Model())
    # Graphs
        VFI_Graphs(M_EGM,"EGM")

# Execute Numerical VFI - ECM
    println("===============================================\n Solving VFI with ECM")
    # Initial value
    s_0  = p.β*p.α
    C_0  = (1-s_0)*M.Y_grid
    # V_0 = d_utility(C_0,p).*M.MPk_mat
    V_0  = (1/(1-p.β))*utility(C_0,p)
    # V_0  = copy(M_EGM.V)
    # V_0 = d_utility(M_EGM.G_c,p).*M.MPk_mat
    # println("V_0="); display(V_0)
    @time M_ECM  = VFI_Fixed_Point(T_ECM,Model(),V_0)
    # Graphs
        VFI_Graphs(M_ECM,"ECM")
