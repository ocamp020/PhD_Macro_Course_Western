# Heterogeneous agent models in Julia
# Sergio Ocampo
# September 2020
# Aiyagari economy with inelastic labor supply
# Solve:   V(ϵ,a) = max{ ((1+r)a+wϵ̄ϵ-a')^(1-γ)/(1-γ) +beta*E[V(ϵ',a')|ϵ] }
#           log(ϵ') = ρ*log(ϵ) + η; η~N(0,σ); ϵ̄=exp(-σ^2/(2(1-ρ^2)))
# The constant ϵ̄ guarantees that E[ϵ]=1 and so aggregate labor L=E[ϵ]=1
cd() # Go to root directory
cd("./Dropbox/Teaching/PhD_Macro_Comp/Julia_Code/Lecture_8_Aiyagari/")
mkpath("Figures")
using SparseArrays
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
println("Aiyagari in Julia")
println("PWD: ",pwd())
println("This code uses Plots, Interpolations, Dierckx, ForwardDiff, Optim, Roots, Parameters, ScaledInterpolation")
println("Solve Aiyagari model with EGM and several implementations of histogram method")
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
        ρ::Float64 = 0.90 ; # Persistence of labor efficiency process
        σ::Float64 = 0.10 ; # Standard devaiation of labor efficiency innovation
        z_bar::Float64 = 1; # Reference level for productivity
        ϵ̄::Float64 = exp(-σ^2/(2*(1-ρ^2))); # Refernce level for labor efficiency
        # Borrowing constraint
        a_min::Float64 = 0; # Borrowing constraint
        # VFI Paramters
        max_iter::Int64   = 100000; # Maximum number of iterations
        dist_tol::Float64 = 1E-6  ; # Tolerance for distance
        # Howard's Policy Iterations
        H_tol::Float64    = 1E-9  ; # Tolerance for policy function iteration
        N_H::Int64        = 20    ; # Maximum number of policy iterations
        # Histogram iteration parameters
        Hist_max_iter     = 10000 ;
        Hist_tol          = 1E-7  ;
        # Histogram iteration parameters
        N_eq              = 1000  ;
        tol_eq            = 1E-7  ;
        η                 = 0.3   ; # Dampen factor for updating capital
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
        a_max::Float64  = 50                         # Max node of a_grid
        θ_a::Float64    = 2.5                        # Curvature of a_grid
        n_a::Int64      = 200                        # Size of a_grid
        n_a_fine::Int64 = 1000                        # Size of fine grid for interpolation and distribution
        a_grid          = Make_Grid(n_a     ,θ_a,p.a_min,a_max,"Poly")  # a_grid for model solution
        a_grid_fine     = Make_Grid(n_a_fine,1  ,p.a_min,a_max,"Poly")  # Fine grid for interpolation
        # Productivity process
        n_ϵ       = 11                               # Size of ϵ_grid
        MP_ϵ      = Rouwenhorst95(p.ρ,p.σ,n_ϵ)       # Markov Process for ϵ
        ϵ_grid    = p.ϵ̄*exp.(MP_ϵ.grid)              # Grid in levels
        # State matrices
        a_mat     = repeat(a_grid',n_ϵ,1)
        a_mat_fine= repeat(a_grid_fine',n_ϵ,1)
        ϵ_mat     = p.ϵ̄*exp.(repeat(MP_ϵ.grid,1,n_a))
        # Prices and aggregates
        r::Float64 = 0.90*(1/p.β - 1) # p.α*p.z_bar*K^(p.α-1) - p.δ
        K::Float64 = (p.α*p.z_bar/(r+p.δ))^(1/(1-p.α)) # k_ss
        Y::Float64 = p.z_bar*K^(p.α)
        w::Float64 = (1-p.α)*p.z_bar*K^(p.α)
        # Value and policy functions
        V         = Array{Float64}(undef,n_ϵ,n_a)       # Value Function
        G_ap      = Array{Float64}(undef,n_ϵ,n_a)       # Policy Function
        G_c       = Array{Float64}(undef,n_ϵ,n_a)       # Policy Function
        V_fine    = Array{Float64}(undef,n_ϵ,n_a_fine)  # Value Function on fine grid
        G_ap_fine = Array{Float64}(undef,n_ϵ,n_a_fine)  # Policy Function on fine grid
        G_c_fine  = Array{Float64}(undef,n_ϵ,n_a_fine)  # Policy Function on fine grid
        # Distribution
        Γ         = 1/(n_ϵ*n_a_fine)*ones(n_ϵ,n_a_fine) # Distribution (initiliazed to uniform)
        # Error in Euler equation
        Euler     = Array{Float64}(undef,n_ϵ,n_a_fine)  # Errors in Euler equation
        # Solver
        Solver    = "PFI"
    end

    # Allocate model to object M for future calling
    M = Model()


#-----------------------------------------------------------
#-----------------------------------------------------------
# Utility function
function utility(c,p::Par)
    if p.γ>1
    return (c).^(1-p.γ)/(1-p.γ)
    else
    return log.(c)
    end
end

function d_utility(c,p::Par)
    return (c).^(-p.γ)
end

function d_utility_inv(x,p::Par)
    return x.^(-1/p.γ)
end


#-----------------------------------------------------------
#-----------------------------------------------------------
# Euler Error Function

# G_ap interpolation
function G_ap_ϵa(i_ϵ::Int64,a,M::Model)
    itp = ScaledInterpolations(M.a_grid,M.G_ap[i_ϵ,:], BSpline(Cubic(Line(OnGrid()))))
    return itp(a)
end

# Euler Equation
function Euler_Error(i_ϵ::Int64,a,M::Model)
    # Return percentage error in Euler equation
    @unpack p, MP_ϵ, n_ϵ, ϵ_grid, r, w = M
    @unpack β, ϵ̄ = p

    # Iterpolate G_ap at current ϵ
    ap  = min(M.a_grid[end],max(M.a_grid[1],G_ap_ϵa(i_ϵ,a,M)))

    # Current consumption
    c   = (1+r)*a + w*ϵ_grid[i_ϵ] - ap

    # Compute left hand side of Euler equation
    LHS = d_utility(c,p)

    # Compute right hand side of Euler equation
        # Marginal utility at ϵ',a',G_ap(ϵ',a')
        up   = zeros(n_ϵ)
        for i_ϵp=1:n_ϵ
        cp       = (1+r)*ap + w*ϵ_grid[i_ϵp] - G_ap_ϵa(i_ϵp,ap,M)
        up[i_ϵp] = d_utility(cp,p)
        end
    RHS = β*(MP_ϵ.Π[i_ϵ,:])'*((1+r).*up)
    # Return percentage errror in Euler equation
    return (RHS/LHS-1)*100
end

# Euler Equation for Optimization
function Euler_Eq(ap,i_ϵ::Int64,i_a::Int64,Va,M::Model)
    # Return percentage error in Euler equation
    @unpack p, MP_ϵ, n_ϵ, a_grid = M
    @unpack β = p
    # Compute left hand side of Euler equation
    c   = (1+r)*a_grid[i_a] + w*ϵ_grid[i_ϵ] - ap
    LHS = d_utility(c,p)
    # Compute right hand side of Euler equation
        # Expected value of derivative β*E[Vk[z',kp]|k]
        Vap  = zeros(n_ϵ)
        for i_ϵp=1:n_ϵ
            Va_ip     = ScaledInterpolations(a_grid,Va[i_ϵp,:], FritschButlandMonotonicInterpolation())
            Vap[i_ϵp] = Va_ip(kp)
        end
    RHS = β*(MP_ϵ.Π[i_ϵ,:])'*Vap
    # Return squared errror in Euler equation
    return (RHS/LHS-1)^2
end

#-----------------------------------------------------------
#-----------------------------------------------------------
# PFI Fixed Point
function PFI_Fixed_Point(T::Function,M::Model,G_ap_old=nothing)
    # Unpack model structure
    @unpack p, n_ϵ, n_a, n_a_fine, θ_a, a_grid, a_grid_fine, r, w = M
    # PFI paramters
    @unpack max_iter, dist_tol = p
    # Initialize variables for loop
    if G_ap_old==nothing
    G_ap_old  = (1+r)*M.a_mat # p.a_min*ones(n_ϵ,n_a) #max.(p.β*((1+r)*M.a_mat+w*M.ϵ_mat),p.a_min) ;
    end
    # println("V_0="); display(V_old)
    G_dist = 1              ; # Initialize distance
    println(" ")
    println("------------------------")
    println("PFI - n_ϵ=$n_ϵ, n_a=$n_a - θ_a=$θ_a - r=$r")
    for iter=1:max_iter
        # println("\n Iteration $iter")
        # Update value function
        G_ap_new, G_c = T(Model(M,G_ap=copy(G_ap_old)))
            # println("T(G_ap) = $G_ap_new")
            # println("  G_ap  = $G_ap_old")
        # Update distance and iterations
        G_dist = sqrt(norm(G_ap_new-G_ap_old,2))
        # Update old function
        G_ap_old  = G_ap_new
        # Report progress
        if mod(iter,250)==0
            println("   PFI Loop: iter=$iter, dist=",G_dist)
        end
        # Check convergence and return results
        if G_dist<=dist_tol
            println("PFI - n_ϵ=$n_ϵ, n_a=$n_a - θ_a=$θ_a - r=$r")
            println("Iterations = $iter and Distance = ",G_dist)
            println("------------------------")
            println(" ")
            # Interpolate to fine grid
            G_ap_fine = zeros(n_ϵ,n_a_fine)
            G_c_fine  = zeros(n_ϵ,n_a_fine)
            for i_ϵ=1:n_ϵ
            G_ap_ip = ScaledInterpolations(a_grid,G_ap_new[i_ϵ,:] , BSpline(Cubic(Line(OnGrid()))))
                G_ap_fine[i_ϵ,:].= G_ap_ip.(collect(a_grid_fine))
            G_c_ip  = ScaledInterpolations(a_grid,G_c[i_ϵ,:]  , BSpline(Cubic(Line(OnGrid()))))
                G_c_fine[i_ϵ,:] .= G_c_ip.(collect(a_grid_fine))
            end
            # Update model
            M = Model(M; G_ap=G_ap_new,G_c=G_c,G_ap_fine=G_ap_fine,G_c_fine=G_c_fine)
            # # Euler Equation Errors
            # Euler = [Euler_Error(i_ϵ,a_grid_fine[i_a],M) for i_ϵ=1:n_ϵ, i_a in 1:n_a_fine]
            # # Update model
            # M = Model(M; Euler=Euler)
            # Return results
            return M
        end
    end
    # If loop ends there was no convergence -> Error!
    error("Error in PFI - Solution not found")
end


# Recover Value Function
function Value_Function(M::Model)
    # Unpack model structure
    @unpack p, n_ϵ, n_a, n_a_fine, a_grid, a_grid_fine, r, w, MP_ϵ, G_ap, G_c = M
    @unpack β, dist_tol = p

    # Compute value function with policy function iteration
    V      = zeros(n_ϵ,n_a)
    V_new  = zeros(n_ϵ,n_a)
    V_dist = 1
    U_mat  = utility(G_c,p)
    while V_dist>dist_tol
        for i_ϵ=1:n_ϵ
            Pr = MP_ϵ.Π[i_ϵ,:]'
        for i_a=1:n_a
            ap = G_ap[i_ϵ,i_a]
            Vp = zeros(n_ϵ,1)
            for i_ϵp=1:n_ϵ
            Vp_ip    = ScaledInterpolations(a_grid,V[i_ϵp,:], BSpline(Cubic(Line(OnGrid()))))
            Vp[i_ϵp] = Vp_ip(ap)
            end
            V_new[i_ϵ,i_a] = U_mat[i_ϵ,i_a] + β*Pr⋅Vp
        end
        end
        V_dist = maximum(abs.(V_new./V.-1))
        V      = V_new
    end
    # Interpolate to fine grid
    V_fine    = zeros(n_ϵ,n_a_fine)
    for i_ϵ=1:n_ϵ
        V_ip    = ScaledInterpolations(a_grid,V[i_ϵ,:], BSpline(Cubic(Line(OnGrid()))))
        V_fine[i_ϵ,:] .= V_ip.(collect(a_grid_fine))
    end
    # Update model
    M = Model(M; V=V,V_fine=V_fine)
    return M
end


#-----------------------------------------------------------
#-----------------------------------------------------------
# VFI Fixed Point
function VFI_Fixed_Point(T::Function,M::Model,V_old=nothing)
    # Unpack model structure
    @unpack p, n_ϵ, n_a, θ_a, a_grid, n_a_fine, a_grid_fine, r, w = M
    # VFI paramters
    @unpack max_iter, dist_tol = p
    # Initialize variables for loop
    if V_old==nothing
    V_old  = utility(M.a_mat+w*M.ϵ_mat,p)  ; # Start at utility with zero savings
    end
    # println("V_0="); display(V_old)

    V_dist = 1              ; # Initialize distance
    println(" ")
    println("------------------------")
    println("VFI - n_ϵ=$n_ϵ, n_a=$n_a - θ_a=$θ_a - r=$r")
    for iter=1:max_iter
        # println("\n Iteration $iter")
        # Update value function
        V_new, G_ap, G_c = T(Model(M,V=copy(V_old)))
            # println("T(V) = $V_new")
            # println("  V  = $V_old")
        # Update distance and iterations
        V_dist = maximum(abs.(V_new./V_old.-1))
        # V_dist = maximum(abs.(G_ap./G_ap_old.-1))
        # Update old function
        V_old  = V_new
        # G_ap_old = copy(G_ap)
        # Report progress
        if mod(iter,100)==0
            println("   VFI Loop: iter=$iter, dist=",100*V_dist,"%")
        end
        # Check convergence and return results
        if V_dist<=dist_tol
            println("VFI - n_ϵ=$n_ϵ, n_a=$n_a - θ_a=$θ_a - r=$r")
            println("Iterations = $iter and Distance = ",100*V_dist,"%")
            println("------------------------")
            println(" ")
            # Interpolate to fine grid
            V_fine    = zeros(n_ϵ,n_a_fine)
            G_ap_fine = zeros(n_ϵ,n_a_fine)
            G_c_fine  = zeros(n_ϵ,n_a_fine)
            for i_z=1:n_ϵ
            V_ip    = ScaledInterpolations(a_grid,V_new[i_ϵ,:], BSpline(Cubic(Line(OnGrid()))))
                V_fine[i_ϵ,:]   .= V_ip.(collect(a_grid_fine))
            G_ap_ip = ScaledInterpolations(a_grid,G_ap[i_ϵ,:] , BSpline(Cubic(Line(OnGrid()))))
                G_ap_fine[i_ϵ,:].= G_ap_ip.(collect(a_grid_fine))
            G_c_ip  = ScaledInterpolations(a_grid,G_c[i_ϵ,:]  , BSpline(Cubic(Line(OnGrid()))))
                G_c_fine[i_ϵ,:] .= G_c_ip.(collect(a_grid_fine))
            end
            # Update model
            M = Model(M; V=V_new,G_ap=G_ap,G_c=G_c,V_fine=V_fine,G_ap_fine=G_ap_fine,G_c_fine=G_c_fine)
            # Return results
            return M
        end
    end
    # If loop ends there was no convergence -> Error!
    error("Error in VFI - Solution not found")
end



#-----------------------------------------------------------
#-----------------------------------------------------------
# Graphs and Stats
function  Aiyagari_Graph(M::Model,H_Type)
    gr()
    # Value Function
        plt = plot(title="Value Function - n_a=$(M.n_a) - θ_a=$(M.θ_a)",legend=:bottomright,foreground_color_legend = nothing,background_color_legend = nothing)
        for i_ϵ=1:M.n_ϵ
        plot!(M.a_grid,M.V[i_ϵ,:],linetype=:scatter,ms=1.5,label="V(ϵ_$i_ϵ)")
        plot!(M.a_grid_fine,M.V_fine[i_ϵ,:],linewidth=1,linestyle=(:dash),linecolor=RGB(0.4,0.4,0.4),label=nothing)
        end
        xlabel!("Capital")
        ylabel!("Value")
        savefig("./Figures/Aiyagari_V_"*H_Type*".pdf")
    # Assets Policy Function
        plt = plot(title="Policy Function - a' - n_a=$(M.n_a) - θ_a=$(M.θ_a)",legend=:bottomright,foreground_color_legend = nothing,background_color_legend = nothing)
        plot!(M.a_grid_fine,M.a_grid_fine,lw=1,linecolor=RGB(0.6,0.6,0.6),label=nothing)
        for i_ϵ=1:M.n_ϵ
        plot!(M.a_grid,M.G_ap[i_ϵ,:],linetype=:scatter,ms=1.5,label="G_ap(ϵ_$i_ϵ)")
        plot!(M.a_grid_fine,M.G_ap_fine[i_ϵ,:],linewidth=1,linestyle=(:dash),linecolor=RGB(0.4,0.4,0.4),label=nothing)
        end
        xlabel!("Assets")
        ylabel!("Assets")
        savefig("./Figures/Aiyagari_G_ap_"*H_Type*".pdf")
    # Disstribution
        Γ_a = sum(M.Γ,dims=1)'
        plt = plot(title="Distribution of Assets - n_a=$(M.n_a) - θ_a=$(M.θ_a)",legend=:bottomright,foreground_color_legend = nothing,background_color_legend = nothing)
        plot!(M.a_grid_fine,Γ_a,lw=2.5,label=nothing)
        xlabel!("Assets")
        ylabel!("Frequencies")
        savefig("./Figures/Aiyagari_Dist_"*H_Type*".pdf")
    # Disstribution
        plt = plot(title="Distribution of Assets by ϵ - n_a=$(M.n_a) - θ_a=$(M.θ_a)",legend=:bottomright,foreground_color_legend = nothing,background_color_legend = nothing)
        plot!(M.a_grid_fine,Γ_a,lw=2.5,linecolor=RGB(0.2,0.2,0.2),label="Γ_a")
        for i_ϵ=1:M.n_ϵ
        plot!(M.a_grid_fine,M.Γ[i_ϵ,:]/sum(M.Γ[i_ϵ,:]),linewidth=1.5,label="Γ_a(ϵ_$i_ϵ)")
        end
        xlabel!("Assets")
        ylabel!("Frequencies")
        savefig("./Figures/Aiyagari_Dist_by_eps_"*H_Type*".pdf")
    println("\n     Graphs Completed for Aiyagari_$Type - n_a=$(M.n_a) - θ_a=$(M.θ_a)\n")
end


#-----------------------------------------------------------
#-----------------------------------------------------------
# Bellman operator - EGM - Iterate on Policy Functions
function T_EGM_G(M::Model)
    @unpack p, n_ϵ, MP_ϵ, n_a, G_ap, r, w = M
    @unpack β, a_min = p
    # Define RHS of Euler equation for each (ϵ,a')
    # Rows are present ϵ and columns are tomorrow's a in fixed grid
    Euler_RHS = β*(1+r)*MP_ϵ.Π*d_utility( (1+r)*M.a_mat + w*M.ϵ_mat - G_ap , p )
    # Check Monotonicity
    if any( Euler_RHS.<0 )
        error("RHS must be monotone for EGM to work")
    end
    # Define consumption from Euler equation
    C_endo = d_utility_inv(Euler_RHS,p)
    # Define endogenous grid on assets
    A_endo = (C_endo .+ M.a_mat - w*M.ϵ_mat)/(1+r)
    # Interpolate functions on exogenous grid
    G_c = Array{Float64}(undef,n_ϵ,n_a)
    for i_ϵ=1:n_ϵ
        # Sort A_endo for interpolation
        sort_ind = sortperm(A_endo[i_ϵ,:])
        A_aux    = A_endo[i_ϵ,:][sort_ind]
        C_aux    = C_endo[i_ϵ,:][sort_ind]
        # Check boundary condition
            # Ap(ϵ,a)=a_min for all a<min(A_aux)
            # Note that in that case C is linear between a_min and min(A_aux)
        if minimum(A_aux)>a_min
            a_vec = M.a_grid[M.a_grid.<minimum(A_aux)]
            A_aux = [a_vec ; A_aux]
            C_aux = [((1+r)*a_vec.+w*M.ϵ_grid[i_ϵ].-a_min) ; C_aux]
        end
        # C_ip        = ScaledInterpolations(A_aux,C_aux, BSpline(Cubic(Line(OnGrid()))))
        C_ip        = Spline1D(A_aux,C_aux)
        G_c[i_ϵ,:] .= C_ip.(M.a_grid)
        Ap_aux      = (1+r)*collect(M.a_grid) .+ w*M.ϵ_grid[i_ϵ] .- G_c[i_ϵ,:]
        # if any(Ap_aux.<a_min)
        # display([M.a_grid G_c[i_ϵ,:] Ap_aux A_endo[i_ϵ,:][sort_ind] C_endo[i_ϵ,:][sort_ind]])
        # println(((1+r)*a_min+w*M.ϵ_grid[i_ϵ]-a_min))
        # end
    end
    # Update policy function
    G_ap .= (1+r)*M.a_mat .+ w*M.ϵ_mat .- G_c
        # Adjust for numerical error
        for ind = findall(<=(1e-10),abs.(G_ap.-a_min))
            G_ap[ind] = a_min
            G_c[ind]  = (1+r)*M.a_mat[ind] + w*M.ϵ_mat[ind] - a_min
        end
        # Check for borrowing constraint
        if any( G_ap.<a_min )
            error("Borrowing Constraint Violated")
        end
    # Return Results
    return G_ap, G_c
end


#-----------------------------------------------------------
#-----------------------------------------------------------
# Bellman operator - EGM - Iterate on Value Function
function T_EGM_V(M::Model)
    @unpack p, n_ϵ, MP_ϵ, n_a, V, G_ap, r, w , a_grid= M
    @unpack β, a_min = p
    # Check Monotonicity
    if any( diff(V,dims=2).<0 )
        error("V must be monotone for EGM to work")
    end
    # Define expectation of value function for each (z,k')
    EV = β*MP_ϵ.Π*V  # Rows are present z and columns are tomorrow's k in fixed grid
    # println("Initial EV size=$(size(EV)), and EV=$EV")
    # Check Monotonicity
    if any( diff(EV,dims=2).<0 )
        error("EV must be monotone for EGM to work")
    end
    # Define the derivative of EV for each value
        # I will be lazy and interpolate to then get the derivative with ForwardDiff
    dEV = zeros(size(EV))
    for i_ϵ=1:n_ϵ
        # EV_ip      = ScaledInterpolations(a_grid,EV[i_ϵ,:], BSpline(Cubic(Line(OnGrid()))))
        # dEV_ip(x)  = ForwardDiff.derivative(EV_ip,x)
        # dEV[i_ϵ,:].= dEV_ip.(a_grid)
        # EV_ip      = Spline1D(a_grid,EV[i_ϵ,:])
        # dEV[i_ϵ,:].= derivative(EV_ip, a_grid )
        EV_ip      = ScaledInterpolations(a_grid,EV[i_ϵ,:], FritschButlandMonotonicInterpolation())
        dEV_ip(x)  = ForwardDiff.derivative(EV_ip,x)
        dEV[i_ϵ,:].= dEV_ip.(a_grid)
        # println("\n i_ϵ=$i_ϵ, dEV=$(dEV[i_ϵ,:]) \n")
    end
    # Check Monotonicity
    if any( dEV.<0 )
        for i_ϵ=1:n_ϵ
            println("\n i_ϵ=$i_ϵ, [min,max]=$([minimum(dEV[i_ϵ,:]),maximum(dEV[i_ϵ,:])]) \n")
        end
        error("dEV must be monotone for EGM to work")
    end
        # gr()
        # plot(title="EV Test")
        # plot!(a_grid,EV')
        # plot(title="dEV Test")
        # plot!(a_grid,dEV')
    # Define Consumption from Euler Equation
    C_endo = d_utility_inv(dEV,p)
    # Define endogenous grid on assets
    A_endo = (C_endo .+ M.a_mat - w*M.ϵ_mat)/(1+r)
    # Define value on endogenous grid
    V_endo = utility(C_endo,p) .+ EV
    # Interpolate functions on exogenous grid
    G_c = Array{Float64}(undef,n_ϵ,n_a)
    for i_ϵ=1:n_ϵ
        # Sort A_endo for interpolation
        sort_ind = sortperm(A_endo[i_ϵ,:])
        A_aux    = A_endo[i_ϵ,:][sort_ind]
        C_aux    = C_endo[i_ϵ,:][sort_ind]
        V_aux    = V_endo[i_ϵ,:][sort_ind]
        # Check boundary condition
            # Ap(ϵ,a)=a_min for all a<min(A_aux)
            # Note that in that case C is linear between a_min and min(A_aux)
        if minimum(A_aux)>a_min
            a_vec = M.a_grid[a_grid.<minimum(A_aux)]
            c_vec = ((1+r)*a_vec.+w*M.ϵ_grid[i_ϵ].-a_min)
            A_aux = [a_vec ; A_aux]
            C_aux = [c_vec ; C_aux]
            V_aux = [utility(c_vec,p).+EV[i_ϵ,1] ; V_aux]
        end
        # C_ip        = ScaledInterpolations(A_aux,C_aux, BSpline(Cubic(Line(OnGrid()))))
        C_ip        = Spline1D(A_aux,C_aux)
        G_c[i_ϵ,:] .= C_ip.(a_grid)
        Ap_aux      = (1+r)*collect(a_grid) .+ w*M.ϵ_grid[i_ϵ] .- G_c[i_ϵ,:]
        # V_ip        = ScaledInterpolations(A_aux,V_endo[i_ϵ,:], BSpline(Cubic(Line(OnGrid()))))
        V_ip        = Spline1D(A_aux,V_aux)
        V[i_ϵ,:]   .= V_ip.(a_grid)
        # if any(Ap_aux.<a_min)
        # display([M.a_grid G_c[i_ϵ,:] Ap_aux A_endo[i_ϵ,:][sort_ind] C_endo[i_ϵ,:][sort_ind]])
        # println(((1+r)*a_min+w*M.ϵ_grid[i_ϵ]-a_min))
        # end
    end
    # Update policy function
    G_ap .= (1+r)*M.a_mat .+ w*M.ϵ_mat .- G_c
        # Adjust for numerical error
        for ind = findall(<=(1e-10),abs.(G_ap.-a_min))
            G_ap[ind] = a_min
            G_c[ind]  = (1+r)*M.a_mat[ind] + w*M.ϵ_mat[ind] - a_min
        end
        # Check for borrowing constraint
        if any( G_ap.<a_min )
            error("Borrowing Constraint Violated")
        end
    # Return Results
    return V, G_ap, G_c
end


#-----------------------------------------------------------
#-----------------------------------------------------------
# Histogram method
function Histogram_Method_Loop(M::Model,N_H=nothing,Γ_0=nothing)
    @unpack p, n_ϵ, MP_ϵ, n_a_fine, a_grid_fine, G_ap_fine = M
    @unpack a_min, Hist_max_iter, Hist_tol = p

    println("\n--------------------------------\nBegining Histogram Method with Loops")

    # Change max iter
    if N_H==nothing
        N_H = Hist_max_iter
    end

    # Initial distribution
    if Γ_0==nothing
        Γ_0 = M.Γ
    end

    # Discretize distribution
    H_ind    = Array{Int64}(undef,n_ϵ,n_a_fine)
    H_weight = Array{Float64}(undef,n_ϵ,n_a_fine)
    a_max    = maximum(a_grid_fine)
    for i_ϵ=1:n_ϵ
    for i_a=1:n_a_fine
        H_ind[i_ϵ,i_a]    = Grid_Inv(G_ap_fine[i_ϵ,i_a],n_a_fine,1,a_min,a_max)
        H_weight[i_ϵ,i_a] = (G_ap_fine[i_ϵ,i_a]-a_grid_fine[H_ind[i_ϵ,i_a]])/(a_max-a_min)
    end
    end
        # Correct corner solutions above
        H_weight[H_ind.==n_a_fine] .= 0
        H_ind[H_ind.==n_a_fine]    .= n_a_fine-1
        # Check bounds for weights
        H_weight = min.(1,max.(0,H_weight))

    # Loop for updating histogram
    H_dist = 1
    for i_H=1:N_H
        # Update histogram
        Γ = zeros(n_ϵ,n_a_fine)
        for i_ϵ=1:n_ϵ # Current ϵ
        for i_a=1:n_a_fine # Current a
            i_ap = H_ind[i_ϵ,i_a]
            ω_ap = H_weight[i_ϵ,i_a]
            for i_ϵp=1:n_ϵ # Future ϵ
                Γ[i_ϵp,i_ap]   = Γ[i_ϵp,i_ap]   +    ω_ap *MP_ϵ.Π[i_ϵ,i_ϵp]*Γ_0[i_ϵ,i_a]
                Γ[i_ϵp,i_ap+1] = Γ[i_ϵp,i_ap+1] + (1-ω_ap)*MP_ϵ.Π[i_ϵ,i_ϵp]*Γ_0[i_ϵ,i_a]
            end
        end
        end
        # Update distance
        H_dist = maximum(abs.(Γ-Γ_0))
        # Update initial distribution
        Γ_0 .= Γ
        # Report progress
        if mod(i_H,250)==0
            println("   Histogram Loop: iter=$i_H, dist=$H_dist")
        end
        # Check convergence
        if H_dist<Hist_tol
            println("Histogram iteartion converged in iteration $i_H. H_dist=$H_dist\n--------------------------------\n")
            M = Model(M; Γ=Γ)
            return M
        end
    end

    # Return Results
    println("Histogram updated for $N_H iteartions. Current H_dist=$H_dist \n--------------------------------\n")
    M = Model(M; Γ=Γ_0)
    return M
end

#-----------------------------------------------------------
#-----------------------------------------------------------
# Histogram method as in Tan (2020)
function Histogram_Method_Tan(M::Model,N_H=nothing,Γ_0=nothing)
    @unpack p, n_ϵ, MP_ϵ, n_a_fine, a_grid_fine, G_ap_fine = M
    @unpack a_min, Hist_max_iter, Hist_tol = p

    println("\n--------------------------------\nBegining Histogram Method with Tan (2020)")

    # Change max iter
    if N_H==nothing
        N_H = Hist_max_iter
    end

    # Initial distribution
    if Γ_0==nothing
        Γ_0 = M.Γ
    end

    # Discretize distribution into transition matrix
    H_ind    = Array{Int64}(undef,n_ϵ*n_a_fine)
    H_weight = Array{Float64}(undef,n_ϵ*n_a_fine)
    a_max    = maximum(a_grid_fine)
    for i_ϵ=1:n_ϵ
    for i_a=1:n_a_fine
        ap = G_ap_fine[i_ϵ,i_a]
        if ap<=a_min
            i_ap = 1
            ω_ap = 1
        elseif ap >=a_max
            i_ap = n_a_fine - 1
            ω_ap = 0
        else
        i_ap = Grid_Inv(ap,n_a_fine,1,a_min,a_max)
        ω_ap = min.(1,max.(0, (ap-a_grid_fine[i_ap])/(a_max-a_min) ))
        end
        # Save index in transition matrix and weight
        H_ind[n_ϵ*(i_a-1)+i_ϵ]    = n_ϵ*(i_ap-1)+i_ϵ
        H_weight[n_ϵ*(i_a-1)+i_ϵ] = ω_ap
    end
    end
    # Assemble transition matrix
    T_a = sparse( repeat(collect(1:n_ϵ*n_a_fine),2) , [H_ind;H_ind.+n_ϵ] , [H_weight; 1 .- H_weight] )

    # Loop for updating histogram
    H_dist = 1
    for i_H=1:N_H
        # Update histogram only for a->a'
        Γ = T_a'*vec(Γ_0)
        # Reshape into n_ϵ by n_a and update ϵ->ϵ'
        Γ = MP_ϵ.Π'*reshape(Γ,n_ϵ,n_a)
        # Update distance
        H_dist = maximum(abs.(Γ-Γ_0))
        # Update initial distribution
        Γ_0 .= Γ
        # Report progress
        if mod(i_H,250)==0
            println("   Histogram Loop: iter=$i_H, dist=$H_dist")
        end
        # Check convergence
        if H_dist<Hist_tol
            println("Histogram iteartion converged in iteration $i_H. H_dist=$H_dist\n--------------------------------\n")
            M = Model(M; Γ=Γ)
            return M
        end
    end

    # Return Results
    println("Histogram updated for $N_H iteartions. Current H_dist=$H_dist \n--------------------------------\n")
    M = Model(M; Γ=Γ_0)
    return M
end

#-----------------------------------------------------------
#-----------------------------------------------------------
# Histogram method
function Histogram_Method(M::Model,Method=nothing,N_H=nothing,Γ_0=nothing)
    # Choose method for solving histogram
    if Method==nothing
        Method = "Loop"
    end

    # Call Histogram
    if Method=="Loop"
        M = Histogram_Method_Loop(M::Model,N_H,Γ_0)
    elseif Method=="Tan"
        M = Histogram_Method_Loop(M::Model,N_H,Γ_0)
    else
        error("No valid method for histogram")
    end

    # Return model
    return M
end

    # # Test Histogram Methods
    # M_aux   = PFI_Fixed_Point(T_EGM_G,Model())
    # @time M_Loop = Histogram_Method(M_aux,"Loop")
    # @time M_Tan  = Histogram_Method(M_aux,"Tan")
    # K_Loop = sum(M_Loop.a_mat_fine.*M_Loop.Γ)
    # K_Tan  = sum(M_Tan.a_mat_fine.*M_Tan.Γ)
    # println("\nHistogram comparisson: K_Loop=$K_Loop, K_Tan=$K_Tan \n")





#-----------------------------------------------------------
#-----------------------------------------------------------
# Excess Demand
function Excess_Demand_K(r,M::Model)
    @unpack p, n_ϵ, n_a_fine, a_mat_fine, Solver = M
    @unpack α, δ, z_bar, N_eq, tol_eq, η = p

    # Capital supply and aggregates at current r
    K_S = ((r+δ)/(α*z_bar))^(-1/(1-α))
    w   = (1-α)*z_bar*K_S^α
    Y   = z_bar*K_S^α

    # Update Model
    M   = Model(M; r=r, w=w, K=K_S, Y=Y)
    # Compute Policy Functions
    if Solver=="PFI"
    M   = PFI_Fixed_Point(T_EGM_G,M)
    elseif Solver=="VFI"
    M   = VFI_Fixed_Point(T_EGM_V,M)
    end
    # Update Distribution
    M   = Histogram_Method(M)
        # println("Distribution check: $(sum(M.Γ)), Total Capital $(sum(a_mat_fine.*M.Γ)), Total Labor $(sum(M.ϵ_grid.*sum(M.Γ,dims=2)))")
    # Caital Demand
    K_D = sum(a_mat_fine.*M.Γ)

    # Return excess demand
    println("\nCapital market at r=$r: K_D=$K_D, K_S=$K_S, error=$(((K_D/K_S)-1)^2) \n")
    return ((K_D/K_S)-1)^2
end

#=
# Graph demand and supply
    r_ss     = 1/p.β-1
    r_vec    = range(0.5*r_ss,0.99*r_ss,length=21)
    r_vec    = range(0.005, 0.04, length = 20)
    K_market = zeros(length(r_vec))
    K_S      = zeros(length(r_vec))
    K_D      = zeros(length(r_vec))
    for i=1:length(r_vec)
        K_S[i]      = ((r_vec[i]+p.δ)/(p.α*p.z_bar))^(-1/(1-p.α))
        # Compute Policy Functions
        M_aux   = PFI_Fixed_Point(T_EGM_G,Model(r=r_vec[i]))
        # Update Distribution
        M_aux   = Histogram_Method(M_aux)
        K_D[i] = sum(M_aux.a_mat_fine.*M_aux.Γ)
        #= Graph (optional)
        # Policy function
            plt = plot(title="Policy Function - a' - n_a=$(M.n_a) - θ_a=$(M.θ_a)",legend=:bottomright,foreground_color_legend = nothing,background_color_legend = nothing)
            plot!(M_aux.a_grid_fine,M_aux.a_grid_fine,lw=1,linecolor=RGB(0.6,0.6,0.6),label=nothing)
            for i_ϵ=1:M.n_ϵ
            plot!(M_aux.a_grid,M_aux.G_ap[i_ϵ,:],linetype=:scatter,ms=1.5,label="G_ap(ϵ_$i_ϵ)")
            plot!(M_aux.a_grid_fine,M_aux.G_ap_fine[i_ϵ,:],linewidth=1,linestyle=(:dash),linecolor=RGB(0.4,0.4,0.4),label=nothing)
            end
            xlabel!("Assets")
            ylabel!("Assets")
        # Distribution
            Γ_a = sum(M_aux.Γ,dims=1)'
            plt = plot(title="Distribution of Assets - n_a=$(M.n_a) - θ_a=$(M.θ_a)",legend=:bottomright,foreground_color_legend = nothing,background_color_legend = nothing)
            plot!(M_aux.a_grid_fine,Γ_a,lw=2.5,label=nothing)
            xlabel!("Assets")
            ylabel!("Frequencies")
        =#
    end
    gr()
    plot(title="Capital Market",legend=:bottomright,foreground_color_legend = nothing,background_color_legend = nothing)
    plot!(K_S,r_vec,label="K Supply",lw=2)
    plot!(K_D,r_vec,label="K Demand",lw=2)
    xlabel!("Capital")
    ylabel!("Interest Rate")
    savefig("./Figures/Aiyagari_K_Market_Loops.pdf")
=#

#-----------------------------------------------------------
#-----------------------------------------------------------
# Aiyagari Equilibrium
function Aiyagari_Equilibrium(M::Model)
    @unpack p, Solver = M
    @unpack β, α, δ, z_bar = p

    # Equilibrium interest rate
    r_ss = 1/β-1
    min_result = optimize(x->Excess_Demand_K(x,M),0.50*r_ss,0.99*r_ss)
    # Check result
    converged(min_result) || error("Failed to clear capital market in $(iterations(min_result)) iterations")
    # Upddate policy function
    r  = min_result.minimizer
    println("Equilibrium found in $(iterations(min_result)) iterations: r=$r")

    # Compute Aggregates
    K   = ((r+δ)/(α*z_bar))^(-1/(1-α))
    w   = (1-α)*z_bar*K^α
    Y   = z_bar*K^α

    # Load aggregates into model
    M = Model(M; r=r, w=w, K=K, Y=Y)

    # Compute Policy Functions
    if Solver=="PFI"
    M   = PFI_Fixed_Point(T_EGM_G,M)
        # Compute Value Function
        M = Value_Function(M)
    elseif Solver=="VFI"
    M   = VFI_Fixed_Point(T_EGM_V,M)
    end

    # Compute Distribution
    M  = Histogram_Method(M)

    # Return Model
    return M

    # No convergence, Display error
    error("No convervence to equilibrium after $N_eq iterations. Current distance of capital: $(100*K_dist)%")
end


#-----------------------------------------------------------
#-----------------------------------------------------------
# Execute Equilibrium and plot graphs

println("===============================================\n Solving Aiyagari with EGM-Histogram(loop)")
@time M_Aiyagari = Aiyagari_Equilibrium(Model(Solver="PFI"))
# Graphs
    Aiyagari_Graph(M_Aiyagari,"Loop")
    # Aiyagari_Stats(M_Aiyagari,"Loop")
