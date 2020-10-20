# Integration for Julia
# Sergio Ocampo
# September 2020
# NGM with log utility and full depreciation
# Solve:   V(k) = max{ log(zk^alpha-k') +beta*V(k') }
cd() # Go to root directory
cd("./Dropbox/Teaching/PhD_Macro_Comp/Julia_Code/Lecture_5_Integration/")
mkpath("Figures")
using Plots
# using LateXStrings # Pkg.add("LaTeXStrings") # https://github.com/stevengj/LaTeXStrings.jl
using Dierckx # Pkg.add("Dierckx") # https://github.com/kbarbary/Dierckx.jl
using Interpolations # Pkg.add("Interpolations") # https://github.com/JuliaMath/Interpolations.jl
using ForwardDiff # Pkg.add("ForwardDiff") # https://github.com/JuliaDiff/ForwardDiff.jl
using Optim # Pkg.add("Optim") # https://julianlsolvers.github.io/Optim.jl/stable/
    using Optim: converged, maximum, maximizer, minimizer, iterations
using Roots # Pkg.add("Roots") # https://github.com/JuliaMath/Roots.jl
using Parameters # Pkg.add("Parameters") # https://github.com/mauro3/Parameters.jl
using Distributions #Pkg.add("Distributions")
using QuadGK # Pkg.add("QuadGK") # https://juliamath.github.io/QuadGK.jl/latest/
using LinearAlgebra
using Random
using Statistics
# Call Scaled Interpolation Functions
    include("../Lecture_3_Interpolation/Scaled_Interpolation_Functions.jl")
println(" ")
println("------------------------")
println("Integration in Julia")
println("PWD: ",pwd())
println("This code uses Plots, Interpolations, Dierckx, ForwardDiff, Optim, ")
println("   Roots, Parameters, ScaledInterpolation, Distributions, QuadGK, LinearAlgebra, Random")
println("Optimization in the context of the neoclassical growth model")
println("------------------------")
println(" ")

#-----------------------------------------------------------
#-----------------------------------------------------------
# Set random seed
Random.seed!(3486);

#-----------------------------------------------------------
# Defien a markov process struct
    # Generate structure for markov processes using Parameters module
    @with_kw struct MP
        # Model Parameters
        N::Int64 # Number of states
        grid     # Grid of discrete markov process
        Π        # Transition matrix
        PDF      # Stationary distribution
        CDF      # Stationary distribution
    end


#-----------------------------------------------------------
#-----------------------------------------------------------
# Tauchen (1986)
    # Objective is to discretize AR(1) process: z'=ρz+η, η~N(0,σ)
    # Code from Kopecky & Suen (2010)
    # Inputs:
        # ρ - Process persisntence
        # σ - Innovation standard deviation
        # N - Size of the grid
        # Ω - Grid expansion in number of standard devaitions (Optional)
    # Outputs:
        # z - Grid of N equally spaced points covering [-Ωσ,Ωσ]
        # Π - Transition matrix, a stochastic matrix (sums to 1 across columns)
        # PDF_z, CDF_z - Stationary distribution of z
function Tauchen86(ρ,σ,N,Ω::Any=3)
    # Create z grid
        z = range(-Ω*σ/sqrt(1-ρ^2),Ω*σ/sqrt(1-ρ^2),length=N)
    # Define intermediate step length
        h = (z[2]-z[1])/2
    # Define auxiliary matrices
        z_0 = repeat(z ,1,N) # Matrix of today's z each row is a value, columns are equal
        z_1 = repeat(z',N,1) # Matrix of tomorrow's z each column is a value, rows are equal
    # Define intervals
        z_lim = zeros(N,N,2) # First matrix is lower bounds. Second matrix is uppor bounds.
        z_lim[:,1      ,1] .= -Inf
        z_lim[:,2:end  ,1] .=  ( z_1[:,2:end  ] - ρ*z_0[:,2:end  ] .- h )./σ
        z_lim[:,1:end-1,2] .=  ( z_1[:,1:end-1] - ρ*z_0[:,1:end-1] .+ h )./σ
        z_lim[:,end    ,2] .=  Inf
    # Define reference distribution
        # This line uses "Distributions"
        F(x) = cdf.(Normal(),x)
    # Fill in transition matrix
        Π_z = F.(z_lim[:,:,2]) - F.(z_lim[:,:,1])
        Π_z = Π_z./repeat(sum(Π_z,dims=2),1,N)
    # Get stationary distribution of markov chain
        PDF_z = real(eigvecs(Π_z')[:,end]); PDF_z = PDF_z/sum(PDF_z) ;
        CDF_z = cumsum(PDF_z)
    # Return
        return MP(N=N,grid=z,Π=Π_z,PDF=PDF_z,CDF=CDF_z)
end


#-----------------------------------------------------------
#-----------------------------------------------------------
# Rouwenhorst (1995)
    # Objective is to discretize AR(1) process: z'=ρz+η, η~N(0,σ)
    # Code from Kopecky & Suen (2010)
    # Inputs:
        # ρ - Process persisntence
        # σ - Innovation standard deviation
        # N - Size of the grid
    # Outputs:
        # z - Grid of N equally spaced points covering [-ψ,ψ]
        # Π - Transition matrix, a stochastic matrix (sums to 1 across columns)
        # PDF_z, CDF_z - Stationary distribution of z
function Rouwenhorst95(ρ,σ,N)
    # Define paramters for Rouwenhorst's approximation
        p = (1+ρ)/2
        q = p                   # Note: I am leaving q here for comparability with source
        ψ = σ*sqrt((N-1)/(1-ρ^2))
        s = (1-q)/(2-(p+q))     # Note: s=0.5, I leave it for comparability with source
    # Fill in transition matrix
    if N==2
        Π_z = [p 1-p ; 1-q q]
    else
        MP_aux = Rouwenhorst95(ρ,σ,N-1)
        o = zeros(N-1)
        Π_z = p*[MP_aux.Π o ; o' 0] + (1-p)*[o MP_aux.Π ; 0 o'] + (1-q)*[o' 0 ; MP_aux.Π o] + q*[0 o' ; o MP_aux.Π]
        # Adjust scale for double counting
        Π_z = Π_z./repeat(sum(Π_z,dims=2),1,N)
    end
    # Distribution
        PDF_z = pdf.(Binomial(N-1,1-s),(0:N-1))
        CDF_z = cumsum(PDF_z)
    # Create z grid
        z    = range(-ψ,ψ,length=N)
    # Return
        return MP(N=N,grid=z,Π=Π_z,PDF=PDF_z,CDF=CDF_z)
end


#-----------------------------------------------------------
#-----------------------------------------------------------
# Simulation of Markov processes

# Simulation function for discrete Markov process
# The result of the simulation is a Markov chain
    # Inputs:
        # Ns - Number of simulated periods
        # z  - Grid with levels of Markov process
        # Π  - Transition matrix of Markov process
        # N_0- Number of periods to throw before reporting simulation (Optional)
    # Output:
        # z_MC - Vector of Ns elemenrts with Markov Chain
function Simulation_MC(Ns,MP::MP,N_0=1000)
    # Compute conditional CDF
    Γ = cumsum(MP.Π,dims=2)
    # Allocate simulated vector
    z_ind    = zeros(Int64,N_0+Ns)
    z_MC     = zeros(N_0+Ns)
    # Starting value for simulation
    z_ind[1] = Int(ceil(length(MP.grid)/2))
    z_MC[1]  = MP.grid[z_ind[1]]
    # Simulate
    for i=2:Ns+N_0
        #= Option 1
        # Draw a uniform random number (r). Compare it with conditional CDF.
        # Conditional CDF given by cumsum of current row of Π, call it Γ
        # Γ[i,j] gives the probability z'<z_j
        # We are looking for j s.t Γ[i,j-1]<r<Γ[i,j]
        # Equivalently, the lowest j s.t. Γ[i,j]-r>0
            z_ind[i] = findmax(sign.(Γ[z_ind[i-1],:] .- rand()))[2]
        =#
        #= Option 2
        # Alternatively we can draw directly from conditional distributional
        # Distributional is categorical with P[z'=z_j] = Π[i,j]
        =#
        z_ind[i] = rand(Categorical(MP.Π[z_ind[i-1],:]))
        z_MC[i]  = MP.grid[z_ind[i]]
    end
    # Throw out first N_0 elements
    z_MC = z_MC[N_0+1:end]
    # Return result
    return z_MC
end

# Moments function for sample
function Moments_MC(z_MC)
    mean_MC      = mean(z_MC)
    std_MC       = std(z_MC)
    auto_corr_MC = cor(z_MC[1:end-1],z_MC[2:end])
    return mean_MC, std_MC, auto_corr_MC
end

#-----------------------------------------------------------
#-----------------------------------------------------------
# Simulation examples
    # Set parameters
    ρ  = 0.95
    σ  = 0.2
    Ns = 10000
    # Get Markov processes with various sizes
    MP_T5  = Tauchen86(ρ,σ,5)
    MP_T11 = Tauchen86(ρ,σ,11)
    MP_T21 = Tauchen86(ρ,σ,21)

    MP_R5  = Rouwenhorst95(ρ,σ,5)
    MP_R11 = Rouwenhorst95(ρ,σ,11)
    MP_R21 = Rouwenhorst95(ρ,σ,21)

    #-------------------------------------------------------
    # Simulate
    T_s_5  = Simulation_MC(Ns,MP_T5 ,0)
    T_s_11 = Simulation_MC(Ns,MP_T11,0)
    T_s_21 = Simulation_MC(Ns,MP_T21,0)

    R_s_5  = Simulation_MC(Ns,MP_R5 ,0)
    R_s_11 = Simulation_MC(Ns,MP_R11,0)
    R_s_21 = Simulation_MC(Ns,MP_R21,0)

    #-------------------------------------------------------
    # Plot paths for simulations, Ns=100
    gr()
    plt = plot(title="Tauchen Simulation Paths",legend=:outerright,foreground_color_legend = nothing,background_color_legend = nothing)
    plot!(1:100,T_s_5[1:100] ,linewidth=3,label="N=5" )
    plot!(1:100,T_s_11[1:100],linewidth=3,label="N=11",linestyle=(:dash))
    plot!(1:100,T_s_21[1:100],linewidth=3,label="N=21",linestyle=(:dot))
    xlabel!("Time")
    savefig("./Figures/Tauchen_Simulation_Paths.pdf")

    gr()
    plt = plot(title="Rouwenhorst Simulation Paths",legend=:outerright,foreground_color_legend = nothing,background_color_legend = nothing)
    plot!(1:100,R_s_5[1:100] ,linewidth=3,label="N=5" )
    plot!(1:100,R_s_11[1:100],linewidth=3,label="N=11",linestyle=(:dash))
    plot!(1:100,R_s_21[1:100],linewidth=3,label="N=21",linestyle=(:dot))
    xlabel!("Time")
    savefig("./Figures/Rouwenhorst_Simulation_Paths.pdf")

    #-------------------------------------------------------
    # Report moments for simulations, Ns=100
    mean_T_5 , std_T_5 , acorr_T_5  = Moments_MC(T_s_5)
    mean_T_11, std_T_11, acorr_T_11 = Moments_MC(T_s_11)
    mean_T_21, std_T_21, acorr_T_21 = Moments_MC(T_s_21)

    mean_R_5 , std_R_5 , acorr_R_5  = Moments_MC(R_s_5)
    mean_R_11, std_R_11, acorr_R_11 = Moments_MC(R_s_11)
    mean_R_21, std_R_21, acorr_R_21 = Moments_MC(R_s_21)

    Moments_Mat = [" " "Data" "T" "T" "T" " " "R" "R" "R";
                   " " " " "N=5" "N=11" "N=21" " " "N=5" "N=11" "N=21";
                   "mean" 0 mean_T_5 mean_T_11 mean_T_21 " " mean_R_5 mean_R_11 mean_R_21;
                   "std"  σ/sqrt(1-ρ^2) std_T_5 std_T_11 std_T_21 " " std_R_5 std_R_11 std_R_21;
                   "acorr" ρ acorr_T_5 acorr_T_11 acorr_T_21 " " acorr_R_5 acorr_R_11 acorr_R_21];

    println("\n Moments from simulated paths: \n")
    display(Moments_Mat)




#-----------------------------------------------------------
#-----------------------------------------------------------
#-----------------------------------------------------------
# Constrained entrepreneur problem - Using integrals to solve VFI
    # V(z,k) = max{c,k'}{ c^(1-γ)/(1-γ) + βE[V(z',k')|z]}
    # s.t.  c+k'=z(k^α)(l^(1-α)) - w*l + (1-δ)k
    #       log(z') = ρlog(z)+η; η~N(0,σ_η)

#-----------------------------------------------------------
# Paramters
    # Generate structure for parameters using Parameters module
    # We can set default values for our parameters
    @with_kw struct Par
        # Model Parameters
        α::Float64 = 1/3  ; # Production function
        β::Float64 = 0.98 ; # Discount factor
        γ::Float64 = 2    ; # Relative risk aversion parameter
        δ::Float64 = 0.05 ; # Depreciation rate
        z_bar::Float64 = 1; # Reference level for productivity
        ρ::Float64 = 0.95 ; # Persistence of productivity process
        σ::Float64 = 0.10 ; # Standard devaiation of productivity innovation
        # Wage level: set as wage in SS of equivalent NGM
        w_bar::Float64 = (1-α)*z_bar*(β*α*z_bar)^(α/(1-α))
        # VFI Paramters
        max_iter::Int64   = 2000  ; # Maximum number of iterations
        dist_tol::Float64 = 1E-9  ; # Tolerance for distance
        # Howard's Policy Iterations
        H_tol::Float64    = 1E-9  ; # Tolerance for policy function iteration
        # Grid points
        N::Int64          = 31    ;
    end
#-----------------------------------------------------------

#-----------------------------------------------------------
# Guess and verify: V(z,k)=v(z)*(k^(1-γ)/(1-γ))
    # v(z) = [ 1 + (βE[v(z')|z])^(1/γ) ]^γ * Γ(z)^(1-γ)
    # Γ(z) = α((1-α)/w_bar)^((1-α)/α)*z^(1/α) + (1-δ)

    #-------------------------------------------------------
    # Define function for profits
    function Γ(z,p::Par=Par())
        @unpack α, δ, w_bar = p
        return α*((1-α)/w_bar)^((1-α)/α)*z.^(1/α) .+ (1-δ)
    end
    #-------------------------------------------------------

    #-------------------------------------------------------
    # Define Bellman operator with Markov Process
    function T_MP(v,MP::MP,p::Par)
        # Check that f has same dimension as the Markov Process
        if length(v)!=MP.N
            error("Expectation requires that dimensiosn agree, v vs MP")
        end
        # Unpack parameters
            @unpack β, γ = p
        # Calculate Bellman operator (done in one step for all states)
            Tv = ( 1 .+ (β* (MP.Π*v) ).^(1/γ) ).^(γ).*(Γ(exp.(MP.grid),p)).^(1-γ)
        # Return
            return Tv
    end
    #-------------------------------------------------------

    #-------------------------------------------------------
    # Define Bellman operator with Gauss-Kronrod Quadrature
    function T_GK(v,log_z,p::Par)
        # Unpack parameters
            @unpack β, γ, ρ, σ = p
        # Define function for interpolation and evaluation out of bounds
            # Note: function evaluated in log(z) bc log(z) is the one with normal innovations
            vp = interpolate(v, BSpline(Cubic(Line(OnGrid()))))
            vp = Interpolations.scale(vp,log_z)
            vp_extrapolate = extrapolate(vp,Line())

        # Define integral limits as 4 std around mean (of zero)
            x_min, x_max = -5*σ, 5*σ
        # Calculate Bellman operator (done in one step for all states)
            Tv = zeros(p.N)
            for i=1:p.N
                # Define integrand
                function F(x)
                    zp = ρ*log_z[i].+x
                    if zp<log_z[1] || zp>log_z[end]
                    return  pdf.(Normal(0,σ),x) .* vp_extrapolate.(zp)
                    else
                    return  pdf.(Normal(0,σ),x) .* vp.(zp)
                    end
                end
            Tv[i] = ( 1 .+ (β* quadgk(F,x_min,x_max)[1] ).^(1/γ) ).^(γ).*(Γ(exp.(log_z[i]),p)).^(1-γ)
            end
        # Return
            return Tv
    end
    #-------------------------------------------------------


    #-------------------------------------------------------
    # Value Function Iteration
    function VFI_Fixed_Point(T::Function,p::Par)
        # VFI paramters
        @unpack max_iter, dist_tol, N = p
        # Initialize variables for loop
        v_old  = zeros(N)     ; # Initialize value function, here I just do 0, not the best option
        V_dist = 1              ; # Initialize distance
        println(" ")
        println("------------------------")
        println("VFI - N=$N ")
        for iter=1:max_iter
            # Update value function
            v_new = T(v_old)
                # println("T(V) = $V_new")
                # println("  V  = $V_old")
            # Update distance and iterations
            V_dist = maximum(abs.(v_new./v_old.-1))
            # Update old function
            v_old  = v_new
            # Report progress
            if mod(iter,100)==0
                println("   VFI Loop: iter=$iter, dist=",100*V_dist,"%")
            end
            # Check convergence and return results
            if V_dist<=dist_tol
                println("VFI - N=$N ")
                println("Iterations = $iter and Distance = ",100*V_dist,"%")
                println("------------------------")
                println(" ")
                # Return results
                return v_old
            end
        end
        # If loop ends there was no convergence -> Error!
        error("Error in VFI - Solution not found")
    end


#-----------------------------------------------------------
#-----------------------------------------------------------
# Expectations on value function

    for rho in [0.5,0.8,0.9,0.95]
    # Set parameters
    p = Par(ρ=rho)
    log_z_grid = range(-4*p.σ,4*p.σ,length=p.N)

    # Solve with Tauchen
    println("\n Solving value with Tauchen approximation \n")
    MP_T = Tauchen86(p.ρ,p.σ,p.N,4)
    @time v_T  = VFI_Fixed_Point(x-> T_MP(x,MP_T,p),p)
    v_T_ip = ScaledInterpolations(MP_T.grid,v_T, BSpline(Cubic(Line(OnGrid())))).(log_z_grid)

    # Solve with Rouwenhorst
    println("\n Solving value with Rouwenhorst approximation \n")
    MP_R = Rouwenhorst95(p.ρ,p.σ,p.N) 
    @time v_R  = VFI_Fixed_Point(x-> T_MP(x,MP_R,p),p)
    v_R_ip = ScaledInterpolations(MP_R.grid,v_R, BSpline(Cubic(Line(OnGrid())))).(log_z_grid)

    # Solve with Quadrature
    println("\n Solving value with Gauss-Kronrod approximation \n")
    @time v_GK  = VFI_Fixed_Point(x-> T_GK(x,log_z_grid,p),p)

    # Plot results
    gr()
    plt = plot(title="Value Function: υ(z)/(1-γ) - ρ=$rho, N=$(p.N)",legend=:bottomright,foreground_color_legend = nothing,background_color_legend = nothing)
    plot!(log_z_grid,v_T_ip/(1-p.γ),linewidth=3,label="Tauchen" )
    plot!(log_z_grid,v_R_ip/(1-p.γ),linewidth=3,label="Rouwenhorst",linestyle=(:dash))
    plot!(log_z_grid,v_GK  /(1-p.γ),linewidth=3,label="Gauss-Kronrod",linestyle=(:dot))
    xlabel!("log(z)")
    savefig("./Figures/Value_Functions_rho_$rho.pdf")

    # Test with ρ=0, all three columns should be the same
    a = ((v_T_ip./((Γ(exp.(log_z_grid),p)).^(1-p.γ))).^(1/p.γ) .- 1).^(p.γ)
    b = ((v_R_ip./((Γ(exp.(log_z_grid),p)).^(1-p.γ))).^(1/p.γ) .- 1).^(p.γ)
    c = ((  v_GK./((Γ(exp.(log_z_grid),p)).^(1-p.γ))).^(1/p.γ) .- 1).^(p.γ)
    gr()
    plt = plot(title="Expectation:  βE[υ(z')|z] - ρ=$rho, N=$(p.N)",legend=:topright,foreground_color_legend = nothing,background_color_legend = nothing)
    plot!(log_z_grid,a,linewidth=3,label="Tauchen" )
    plot!(log_z_grid,b,linewidth=3,label="Rouwenhorst",linestyle=(:dash))
    plot!(log_z_grid,c,linewidth=3,label="Gauss-Kronrod",linestyle=(:dot))
    xlabel!("log(z)")
    savefig("./Figures/Expectations_rho_$rho.pdf")

    end
