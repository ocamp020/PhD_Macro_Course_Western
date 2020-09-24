# Optimization and Root Finding for Julia
# Sergio Ocampo
# September 2020
# NGM with log utility and full depreciation
# Solve:   V(k) = max{ log(zk^alpha-k') +beta*V(k') }
cd() # Go to root directory
cd("./Dropbox/Teaching/PhD_Macro_Comp/Julia_Code/Lecture_4_Optimization/")
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
# Call Scaled Interpolation Functions
    include("../Lecture_3_Interpolation/Scaled_Interpolation_Functions.jl")
println(" ")
println("------------------------")
println("Optimization in Julia")
println("PWD: ",pwd())
println("This code uses Plots, Interpolations, Dierckx, ForwardDiff, Optim, Roots, Parameters, ScaledInterpolation")
println("Optimization in the context of the neoclassical growth model")
println("------------------------")
println(" ")


#-----------------------------------------------------------
#-----------------------------------------------------------
# Minimum Bracketing Function
    # This function comes from Numerical Recipes, Section 10.1
    # Input: numbers a,b that initiliaze bracket and objective function F
        # Optional: Bounds for the brackeet so that a,b,c ∈ [x_min,x_max]
    # Output: numbers (a,b,c) and images (fa,fb,fc) such that a<b<c and fb<fa,fc
    function mnbrak(a,b,F::Function,x_min=-Inf,x_max=Inf)
        # Auxiliariy variables
        Tiny    = 1E-20
        G_limit = 100
        # Define golden ratio
        γ = 1.618034
        # Evaluate function at end points
        fa, fb = F(a), F(b)
        # Swap a and b so that we can go downhill in the direction from a to b.
            # This way we know for certain that fb<fa, we only need to find c such that fb<fc
        if fb>fa
        a , b = b , a
        fa,fb = fb, fa
        end
        # Guess for value c expanding bracket away from b -> (a,b,c) or (c,b,a)
            # Check bounds
        c  = max( min( b + γ*(b-a) , x_max ) , x_min)
        fc = F(c)
        # While fb>fc we need to keep on bracketing
        while fb>fc
            # Compute u (new candidate) by parabolic extrapolation from a, b, c.
                # TINY is used to prevent any possible division by zero.
            r = (b-a)*(fb-fc)
            q = (b-c)*(fb-fa)
            u = min(max(b-((b-c)*q-(b-a)*r)/(2*sign(q-r)*max(abs(q-r),Tiny)),x_min),x_max)
            u_lim = min(max(b+G_limit*(c-b),x_min),x_max)
            # Test various cases for new candidate
            if ((b-u)*(u-c) > 0) # Parabolic u is between b and c
                fu=F(u)
                if (fu < fc) # Got a minimum between b and c so bounds are (b,u,c)
                    a, fa, b, fb = b, fb, u, fu
                    break
                elseif (fu > fb) # Got a minimum between a and u so bounds are (a,b,u)
                    c, fc = u, fu
                    break
                end
                # If you got here is because candidate u failed
                # Get new candidate by expanding interval with golden ratio
                # Check bounds
                u  = max(min( c+γ*(c-b) , x_max ),x_min)
                fu = F(u)
            elseif ((c-u)*(u-u_lim) > 0.0) # Parabolic u is between c and its limit (ulim)
                fu=F(u)
                if (fu < fc) # Drop c and replace it with u, get new u with golden expansion
                    b, c, fb, fc = c, u, fc, fu
                    u  = max(min( c+γ*(c-b) , x_max ),x_min)
                    fu = F(u)
                end
            elseif ((u-u_lim)*(u_lim-c) >= 0.0) # U is above its limit, reign it in!
                u  = u_lim
                fu = F(u)
            else # Nothing worked, use golden expansion
                u  = max(min( c+γ*(c-b) , x_max ),x_min)
                fu = F(u)
            end
            # If no break then forget the oldest point and go onto next iteration
                # This means assigning b->a, c->b, u-> and forgetting about a
            a, b, c, fa, fb, fc = b, c, u, fb, fc, fu
        end
        # Return solution once out of the loop
            # Swap a and c so that a<b<c
            if a>c
            a , c  = c, a
            fa, fc = fc, fa
            end
        # Minimum bracketed in [a,c] with intermediate point b so that fb<fa,fc
        # println("Minimum bracketed in [$a,$c] with intermediate point b=$b \n Function values: F(a)=$fa, F(b)=$fb, F(c)=$fc")
        return a,c,b,fa,fc,fb
    end


#-----------------------------------------------------------
#-----------------------------------------------------------
# Paramters and Model Structure
    # Generate structure for parameters using Parameters module
    # We can set default values for our parameters
    @with_kw struct Par
        # Model Parameters
        z::Float64 = 1    ; # Productivity
        α::Float64 = 1/3  ; # Production function
        β::Float64 = 0.98 ; # Discount factor
        # VFI Paramters
        max_iter::Int64   = 2000  ; # Maximum number of iterations
        dist_tol::Float64 = 1E-9  ; # Tolerance for distance
        # Howard's Policy Iterations
        H_tol::Float64    = 1E-9  ; # Tolerance for policy function iteration
        N_H::Int64        = 20    ; # Maximum number of policy iterations
        # Minimum consumption for numerical optimization
        c_min::Float64    = 1E-16
    end

    # Allocate paramters to object p for future calling
    p = Par()



#-----------------------------------------------------------
#-----------------------------------------------------------
# Utility function
function utility(k,kp,p::Par)
    @unpack z, α, c_min = p
    c = z*k.^α  .- kp
    if c>c_min
    return log(c)
    else
    return log(c_min)
    end
end

function d_utility(k,kp,p::Par)
    @unpack z, α = p
    c = z*k.^α  - kp
    if c>0
    return -1/c
    else
    return -Inf
    end
end

#-----------------------------------------------------------
#-----------------------------------------------------------
# Steady state values (funciton)
function SS_values(p::Par)
    # This function takes in parameters and provides steady state values
    # Parameters: productivity (z), returns to scale (a) and discount factor (b)
    # Output: values for capital, production, consumption, rental rate, wage
    @unpack z, α, β = p
    k_ss = (β*α*z)^(1/(1-α))
    y_ss = z*k_ss^α
    c_ss = y_ss - k_ss
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

#-----------------------------------------------------------
#-----------------------------------------------------------
# Grid
    function Make_K_Grid(n_k,θ_k,p::Par)
        # Get SS
        k_ss,y_ss,c_ss,r_ss,w_ss = SS_values(p)
        # Get k_grid
        if θ_k≠1
        # k_grid = ExpRange(1E-5,2*k_ss;θ=θ_k,N=n_k) ; # Curved grid between 0 and 2*k_ss
        k_grid = PolyRange(1E-5,2*k_ss;θ=θ_k,N=n_k) ; # Curved grid between 0 and 2*k_ss
        else
        k_grid = range(1E-5,2*k_ss,length=n_k)
        end
        # Return
        return k_grid
    end


# Generate structure of model objects
    @with_kw struct Model
        # Parameters
        p::Par = Par() # Model paramters in their own structure
        # Grids
        θ_k::Float64    = 1     # Curvature of k_grid
        n_k::Int64      = 20    # Size of k_grid
        n_k_fine::Int64 = 1000  # Size of fine grid for interpolation
        k_grid          = Make_K_Grid(n_k,θ_k,p)    # k_grid for model solution
        k_grid_fine     = Make_K_Grid(n_k_fine,1,p) # Fine grid for interpolation
        # Value and policy functions
        V         = Array{Float64}(undef,n_k)       # Value Function
        G_kp      = Array{Float64}(undef,n_k)       # Policy Function
        G_c       = Array{Float64}(undef,n_k)       # Policy Function
        V_fine    = Array{Float64}(undef,n_k_fine)  # Value Function on fine grid
        G_kp_fine = Array{Float64}(undef,n_k_fine)  # Policy Function on fine grid
        G_c_fine  = Array{Float64}(undef,n_k_fine)  # Policy Function on fine grid
        # Anaytical Solutions
        V_a       = Array{Float64}(undef,n_k_fine)  # Analytical Value Function on fine grid
        G_kp_a    = Array{Float64}(undef,n_k_fine)  # Analytical Policy Function on fine grid
        G_c_a     = Array{Float64}(undef,n_k_fine)  # Analytical Policy Function on fine grid
        Euler     = Array{Float64}(undef,n_k_fine)  # Errors in Euler equation
    end

    M = Model()


#-----------------------------------------------------------
#-----------------------------------------------------------
# Analytical solution

# Value function
function V_analytical(k,p::Par)
    # This function takes in paramters and a value of capital
    # Output is the value of the value function at that level of capital
    # Parameters: productivity (z), returns to scale (a) and discount factor (b)
    @unpack z, α, β = p
    a_1 = α/(1-α*β)
    a_0 = (1/(1-β))*(log(z)+log(1/(1+β*a_1))+β*a_1*log(β*a_1*z/(1+β*a_1*z)))
    V   = a_0 .+ a_1*log.(k)
    return V
end

# Policy functions
function G_analytical(k,p::Par)
    # This function takes in paramters and a value of capital
    # Output is the a tuple of optimal capital (kp) and consumption (c)
    # Parameters: productivity (z), returns to scale (a) and discount factor (b)
    @unpack z, α, β = p
    a_1 = α/(1-α*β)
    kp  = β*α*z*k.^α
    c   = z*k.^α .- kp
    return kp, c
end

# Euler Equation
function Euler_Error(k,kp,kpp,p::Par)
    # Return percentage error in Euler equation
    @unpack z, α, β = p
    LHS = 1/(z*k^α-kp)
    RHS = β*α*z*kp^(α-1)/(z*kp^α-kpp)
    return (RHS/LHS-1)*100
end

# Analytical Suite
function VFI_Analytical_Results(M::Model)
    @unpack n_k_fine, k_grid_fine, p = M
    # Analytical value function
    V_a = V_analytical(k_grid_fine,p)
    # Analytical policy function
    G_kp_a, G_c_a = G_analytical(k_grid_fine,p)
    # Interpolation of G_kp
    G_kp_ip = ScaledInterpolations(M.k_grid,M.G_kp, BSpline(Cubic(Line(OnGrid()))))
    # Euler error of numerical policy function on grid
    Euler = zeros(n_k_fine)
    for i=1:n_k_fine
        k   = k_grid_fine[i]
        kp  = G_kp_ip(k)
        kpp = G_kp_ip(kp)
        Euler[i] = Euler_Error(k,kp,kpp,p)
    end
    # Return
    return V_a, G_kp_a, G_c_a, Euler
end

#-----------------------------------------------------------
#-----------------------------------------------------------
# VFI Fixed Point
function VFI_Fixed_Point(T::Function,M::Model)
    # Unpack model structure
    @unpack p, n_k, θ_k, k_grid = M
    # VFI paramters
    @unpack max_iter, dist_tol = p
    # Initialize variables for loop
    V_old  = zeros(n_k)     ; # Initialize value function, here I just do 0, not the best option
    # V_old  = utility.(collect(k_grid),zeros(n_k),p) ; # Start at utility with zero savings
    V_dist = 1              ; # Initialize distance
    println(" ")
    println("------------------------")
    println("VFI - n_k=$n_k - θ_k=$θ_k")
    for iter=1:max_iter
        # Update value function
        V_new, G_kp, G_c = T(Model(M,V=copy(V_old)))
            # println("T(V) = $V_new")
            # println("  V  = $V_old")
        # Update distance and iterations
        V_dist = maximum(abs.(V_new./V_old.-1))
        # Update old function
        V_old  = V_new
        # Report progress
        if mod(iter,100)==0
            println("   VFI Loop: iter=$iter, dist=",100*V_dist,"%")
        end
        # Check convergence and return results
        if V_dist<=dist_tol
            println("VFI - n_k=$n_k - θ_k=$θ_k")
            println("Iterations = $iter and Distance = ",100*V_dist,"%")
            println("------------------------")
            println(" ")
            # Interpolate to fine grid
            V_ip = ScaledInterpolations(M.k_grid,V_new, BSpline(Cubic(Line(OnGrid()))))
                V_fine = V_ip.(collect(M.k_grid_fine))
            G_kp_ip = ScaledInterpolations(M.k_grid,G_kp, BSpline(Cubic(Line(OnGrid()))))
                G_kp_fine = G_kp_ip.(collect(M.k_grid_fine))
            G_c_ip = ScaledInterpolations(M.k_grid,G_c, BSpline(Cubic(Line(OnGrid()))))
                G_c_fine = G_c_ip.(collect(M.k_grid_fine))
            # Update model
            M = Model(M; V=V_new,G_kp=G_kp,G_c=G_c,V_fine=V_fine,G_kp_fine=G_kp_fine,G_c_fine=G_c_fine)
            # Analytical results
            V_a, G_kp_a, G_c_a, Euler = VFI_Analytical_Results(M)
            # Update model
            M = Model(M; V_a=V_a,G_kp_a=G_kp_a,G_c_a=G_c_a,Euler=Euler)
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
    # Value Function Analytical vs 200
        plot(M.k_grid_fine,M.V_a,linewidth=4,label = "Analytical Solution",title="Value Function",legend=(0.75,0.2),foreground_color_legend = nothing,background_color_legend = nothing)
        plot!(M.k_grid,M.V,linetype=:scatter,marker=(:diamond,4),markercolor=RGB(0.5,0.1,0.1),label="VFI - n_k=$(M.n_k) - θ_k=$(M.θ_k)")
        plot!(M.k_grid_fine,M.V_fine,linewidth=2.5,linestyle=(:dash),linecolor=RGB(0.4,0.4,0.4),label=nothing)
        xlabel!("Capital")
        ylabel!("Value")
        savefig("./Figures/VFI_"*VFI_Type*"_V_$(M.n_k)_$(M.θ_k).pdf")
    # Capital Policy Function Analytical vs 200
        plot(M.k_grid_fine,M.k_grid_fine,lw=1,linecolor=RGB(0.6,0.6,0.6),label=nothing)
        plot!(M.k_grid_fine,M.G_kp_a,linewidth=3,label = "Analytical Solution",title="Policy Function - K",legend=(0.75,0.2),foreground_color_legend = nothing,background_color_legend = nothing)
        plot!(M.k_grid,M.G_kp,linetype=:scatter,marker=(:diamond,4),markercolor=RGB(0.5,0.1,0.1),label="VFI - n_k=$(M.n_k) - θ_k=$(M.θ_k)")
        plot!(M.k_grid_fine,M.G_kp_fine,linewidth=2.5,linestyle=(:dash),linecolor=RGB(0.4,0.4,0.4),label=nothing)
        xlabel!("Capital")
        ylabel!("Capital")
    savefig("./Figures/VFI_"*VFI_Type*"_G_kp_$(M.n_k)_$(M.θ_k).pdf")
        # Euler Error 200
        plot(M.k_grid_fine,zeros(M.n_k_fine),lw=1,linecolor=RGB(0.6,0.6,0.6),label=nothing,title = "Euler Equation Error (%)",legend=(0.75,0.2),foreground_color_legend = nothing,background_color_legend = nothing)
        plot!(M.k_grid_fine,M.Euler,linetype=:scatter,marker=(:diamond,2),linecolor=RGB(0.5,0.1,0.1),label="VFI - n_k=$(M.n_k) - θ_k=$(M.θ_k)")
        xlabel!("Capital")
        ylabel!("Percentage Points")
    savefig("./Figures/VFI_"*VFI_Type*"_Euler_$(M.n_k)_$(M.θ_k).pdf")

    println("\n     Graphs Completed for VFI_$VFI_Type - n_k=$(M.n_k) - θ_k=$(M.θ_k)\n")
end


#-----------------------------------------------------------
#-----------------------------------------------------------
# Bellman operator - Continuous Choice
function T_cts_max(M::Model)
    @unpack p, n_k, k_grid, V, G_kp, G_c = M
    @unpack z, α, β, c_min = p
    # println("Initial V = $V")
    # Define the interpolation of the current value function (V) for values of Vp
    Vp = ScaledInterpolations(k_grid,V, BSpline(Cubic(Line(OnGrid()))))
        # println("Interpolation test: Vp(k_min)=",Vp(k_grid[1])," Vp(k/2)=",Vp((k_grid[1]+k_grid[end])/2))
    # Define the derivative (might not be used, depends on algorithm)
        # For many methods we don't need to use ForwardDiff, we have direct measures,
        # for example for cubic splines. I am using ForwardDiff as an example
    # dVp(x) = ForwardDiff.derivative(Vp,x)
    for i = 1:n_k
        # Define objective function: Right hand side of Bellman operator
            # Optimizer will minimize so we need the negative of the function
        Obj_Fun = x->(-utility(k_grid[i],x,p) - β*Vp.(x))
            # println("Interpolation test: Obj_Fun(k_min)=",Obj_Fun(k_grid[1])," Obj_Fun(k/2)=",Obj_Fun((k_grid[1]+k_grid[end])/2))
        # Min and max kp given current k
        kp_min = k_grid[1]
        kp_max = z*k_grid[i]^α - c_min
            # Note: we can be smarter here.
            # We know that if k_grid[i]<k_ss then the kp_max is actually k_ss
            # Similarly, if k_grid[i]>k_ss then kp_min = k_ss
        # Check if lower bound binds
            # No need in this problem because without capital there is no production
            # also no consumption. So value would be -Inf.
            # I am leaving two ways of getting the derivative
        # dObj_min = ForwardDiff.derivative(Obj_Fun,kp_min)
        # dObj_min = -d_utility(k_grid[i],kp_min,p) - β*dVp(kp_min)
        dObj_min = -1
        if dObj_min>0
            V[i]    = -Obj_Fun(kp_min)
            G_kp[i] = kp_min
        else
        # Check if upper bound binds
            # This can bind but shouldn't I will leave code to check but comment it
            # I am leaving two ways of getting the derivative
        # dObj_max = ForwardDiff.derivative(Obj_Fun,kp_max)
        # dObj_max = -d_utility(k_grid[i],kp_max,p) - β*dVp(kp_max)
        dObj_max = 1
        if dObj_max<0
            V[i]    = -bj_Fun(kp_max)
            G_kp[i] = kp_max
        else
        # Bracket the solution for k=k_grid[i]
            # In this case this is not necessary, we know the min is bracketed by [kp_min and kp_max]
            # We can still try to use mnbrak to get a better bracket, but in this case the algorithm
            # stops in the initial evaluation, it just verifies that we do have a bracket
            kp_min, kp_max = mnbrak(kp_min,(kp_max+kp_min)/2,Obj_Fun,kp_min,kp_max)
        # Maximize
            min_result = optimize(Obj_Fun,kp_min,kp_max)
        # Check result
            converged(min_result) || error("Failed to solve Bellman max in $(iterations(min_result)) iterations")
        # Record results
            V[i]     = -min_result.minimum
            G_kp[i]  = min_result.minimizer
            # println("   Maximization result k_grid[$i] - kp=",min_result.minimizer," V(k)=",min_result.minimum," in $(iterations(min_result)) iterations")
        end # Upper bound if  
        end # Lower bound if
    end # loop of k_grid[i]
    # Fill in policy for consumption
    G_c = z.*collect(k_grid).^α .- G_kp
    # Return Results
        # println("T(V) = $V")
    return V, G_kp, G_c
end


#-----------------------------------------------------------
#-----------------------------------------------------------
# Bellman operator - Continuous Choice
function T_cts_root(M::Model)
    @unpack p, n_k, k_grid, V, G_kp, G_c = M
    @unpack z, α, β, c_min = p
    # println("Initial V = $V")
    # Define the interpolation of the current value function (V) for values of Vp
    Vp = ScaledInterpolations(k_grid,V, BSpline(Cubic(Line(OnGrid()))))
    # Define the derivative
        # For many methods we don't need to use ForwardDiff, we have direct measures,
        # for example for cubic splines. I am using ForwardDiff as an example
    dVp(x) = ForwardDiff.derivative(Vp,x)
    for i = 1:n_k
        # Define objective function: Right hand side of Bellman operator
            # Optimizer will minimize so we need the negative of the function
        Obj_Fun   = x->( -utility(k_grid[i],x,p) - β*Vp.(x)) # Note the negative sign
        d_Obj_Fun = x->(d_utility(k_grid[i],x,p) + β*dVp.(x)) # Note there is no negative sign
        # d_Obj_Fun(x) = ForwardDiff.derivative(Obj_Fun,x)
            # println("Interpolation test: Obj_Fun(k_min)=",Obj_Fun(k_grid[1])," Obj_Fun(k/2)=",Obj_Fun((k_grid[1]+k_grid[end])/2)," Vp(k_min)=",Vp(k_grid[1]))
            # println("Interpolation test: dObj_Fun(k_min)=",d_Obj_Fun(k_grid[1])," dObj_Fun(k/2)=",d_Obj_Fun((k_grid[1]+k_grid[end])/2)," dVp(k_min)=",dVp(k_grid[1]))
        # Min and max kp given current k
        kp_min = 1.001*k_grid[1]
        kp_max = min( z*k_grid[i]^α - c_min , k_grid[end] )
        # println("\n dObj_Fun=$(d_Obj_Fun.(range(kp_min,kp_max,length=15)))")
            # Note: we can be smarter here.
            # We know that if k_grid[i]<k_ss then the kp_max is actually k_ss
            # Similarly, if k_grid[i]>k_ss then kp_min = k_ss
        # Check if lower bound binds
            # No need in this problem because without capital there is no production
            # also no consumption. So value would be -Inf.
        dObj_min = d_Obj_Fun(kp_min)
        if dObj_min<0
            V[i]    = -Obj_Fun(kp_min)
            G_kp[i] = kp_min
        elseif isnan(dObj_min)
            error("dObj_min=NaN")
        else
        # Check if upper bound binds
            # This can bind but shouldn't I will leave code to check but comment it
        dObj_max = d_Obj_Fun(kp_max)
        if dObj_max>0
            V[i]    = -Obj_Fun(kp_max)
            G_kp[i] = kp_max
        elseif isnan(dObj_max)
            errlr("dObj_max=NaN")
        else
        # Find root (also check for bracketing)
            if dObj_min*dObj_max<0
            # println("kp_min=$kp_min, c=$(z*k_grid[1]^α-kp_min) dObj_min=$dObj_min, dVp=$(dVp(kp_min)), dU=$(d_utility(k_grid[1],kp_min,p))")
            # println("kp_max=$kp_max, c=$(z*k_grid[1]^α-kp_max) dObj_max=$dObj_max, dVp=$(dVp(kp_max)), dU=$(d_utility(k_grid[1],kp_max,p))")
            G_kp[i] = find_zero(d_Obj_Fun,(kp_min,kp_max),Roots.Brent())
            else
            error("Zero not bracketed by kp_min=$kp_min and kp_max=$kp_max, DF_min=$dObj_min, DF_max=$dObj_max")
            end
        # Check result
            if abs(d_Obj_Fun(G_kp[i]))>1E-6
            error("Failed to solve Bellman max")
            end
        # Record results
            V[i]    = -Obj_Fun(G_kp[i])
            # println("   Maximization result k_grid[$i] - kp=",min_result.minimizer," V(k)=",min_result.minimum," in $(iterations(min_result)) iterations")
        end # Upper bound if  
        end # Lower bound if
    end # loop of k_grid[i]
    # Fill in policy for consumption
    G_c = z.*collect(k_grid).^α .- G_kp
    # Return Results
        # println("T(V) = $V")
    return V, G_kp, G_c
end




#-----------------------------------------------------------
#-----------------------------------------------------------
# Execute VFI and plot graphs

# Execute Numerical VFI - equally spaced grid
    @time M_20  = VFI_Fixed_Point(T_cts_max,Model(n_k=20))
    @time M_50  = VFI_Fixed_Point(T_cts_max,Model(n_k=50))
    @time M_100 = VFI_Fixed_Point(T_cts_max,Model(n_k=100))
    # Graphs
        VFI_Graphs(M_20,"cts_max")
        VFI_Graphs(M_50,"cts_max")
        VFI_Graphs(M_100,"cts_max")
# Execute Numerical VFI - curved spaced grid
    @time M_20c  = VFI_Fixed_Point(T_cts_max,Model(n_k=20,θ_k=2.5))
    @time M_50c  = VFI_Fixed_Point(T_cts_max,Model(n_k=50,θ_k=2.5))
    @time M_100c = VFI_Fixed_Point(T_cts_max,Model(n_k=100,θ_k=2.5))
    # Graphs
        VFI_Graphs(M_20c,"cts_max")
        VFI_Graphs(M_50c,"cts_max")
        VFI_Graphs(M_100c,"cts_max")


# Execute Numerical VFI - equally spaced grid
    # @time M_20  = VFI_Fixed_Point(T_cts_root,Model(n_k=20))
    # @time M_50  = VFI_Fixed_Point(T_cts_root,Model(n_k=50))
    @time M_100 = VFI_Fixed_Point(T_cts_root,Model(n_k=100))
     # Graphs
        # VFI_Graphs(M_20,"cts_root")
        # VFI_Graphs(M_50,"cts_root")
        VFI_Graphs(M_100,"cts_root")

# Execute Numerical VFI - curved spaced grid
    @time M_20c  = VFI_Fixed_Point(T_cts_root,Model(n_k=20,θ_k=2.5))
    @time M_50c  = VFI_Fixed_Point(T_cts_root,Model(n_k=50,θ_k=2.5))
    @time M_100c = VFI_Fixed_Point(T_cts_root,Model(n_k=100,θ_k=2.5))
    # Graphs
        VFI_Graphs(M_20c,"cts_root")
        VFI_Graphs(M_50c,"cts_root")
        VFI_Graphs(M_100c,"cts_root")
