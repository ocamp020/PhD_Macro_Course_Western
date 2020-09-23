# VFI Program for Julia
# Sergio Ocampo
# September 2020
# NGM with log utility and full depreciation
# Solve:   V(k) = max{ log(zk^alpha-k') +beta*V(k') }
# Solution by grid search: k_grid = {k_1,...,k_I}
cd() # Go to root directory
cd("./Dropbox/Teaching/PhD_Macro_Comp/Julia_Code/Lecture_2_VFI/")
mkpath("Figures")
using Plots
using Parameters # Pkg.add("Parameters") # https://github.com/mauro3/Parameters.jl
println(" ")
println("------------------------")
println("VFI code for NGM with log-utility and full depreciation")
println("PWD: ",pwd())
println("This code uses Plots and Parameters")
println("------------------------")
println(" ")

#-----------------------------------------------------------
#-----------------------------------------------------------
# Paramters
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
    end

    # Allocate paramters to object p for future calling
    p = Par()

#-----------------------------------------------------------
#-----------------------------------------------------------
# Utility function
function utility(k,kp,p::Par)
    @unpack z, α = p
    c = z*k.^α  - kp
    if c>0
    return log(c)
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
function VFI_Analytical_Results(n_k::Int64,G_kp,p::Par)
    # Get Grid
    k_grid = Make_K_Grid(n_k,p)
    # Analytical value function
    V_a = V_analytical(k_grid,p)
    # Analytical policy function
    G_kp_a, G_c_a = G_analytical(k_grid,p)
    # Euler error of numerical policy function on grid
    Euler = zeros(n_k)
    for i=1:n_k
        k   = k_grid[i]
        kp  = k_grid[G_kp[i]]
        kpp = k_grid[G_kp[G_kp[i]]]
        Euler[i] = Euler_Error(k,kp,kpp,p)
    end
    # Return
    return V_a, G_kp_a, G_c_a, Euler
end

#-----------------------------------------------------------
#-----------------------------------------------------------
# VFI with Grid
    # Grid
    function Make_K_Grid(n_k,p::Par)
        # Get SS
        k_ss,y_ss,c_ss,r_ss,w_ss = SS_values(p)
        # Get k_grid
        k_grid = range(1E-5,2*k_ss;length=n_k) ; # Equally spaced grid between 0 and 2*k_ss
        # Return
        return k_grid
    end

    # Solve VFI with grid search and loops
    function VFI_grid(T::Function,k_grid,p::Par)
        # VFI paramters
        @unpack max_iter, dist_tol = p
        # Initialize variables for loop
        n_k    = length(k_grid) ; # Number of grid nodes
        V_old  = zeros(n_k)     ; # Initial value, a vector of zeros
        V_dist = 1              ; # Initialize distance
        println(" ")
        println("------------------------")
        println("VFI - Grid Search - n_k=$n_k")
        for iter=1:max_iter
            # Update value function
            V_new, G_kp, G_c = T(V_old)
            # Update distance and iterations
            V_dist = maximum(abs.(V_new./V_old.-1))
            iter  += 1
            # Update old function
            V_old  = V_new
            # Report progress
            if mod(iter,100)==0
                println("   VFI Loop: iter=$iter, dist=",100*V_dist,"%")
            end
            # Check convergence and return results
            if V_dist<=dist_tol
                println("VFI - Grid Search - n_k=$n_k")
                println("Iterations = $iter and Distance = ",100*V_dist,"%")
                println("------------------------")
                println(" ")
                return V_new, G_kp, G_c
            end
        end
        # If loop ends there was no convergence -> Error!
        error("Error in VFI - Grid Search - Solution not found")
    end

#-----------------------------------------------------------
#-----------------------------------------------------------
# VFI - Grid Search - Loops

println(" "); println(" "); println(" ")
println("----------------------------------------------------")
println("----------------------------------------------------")
println("Value Function Iteration - Grid Search with Loops")



# Define function for Value update and policy functions
    function T_grid_loop(V_old,k_grid,p::Par)
        @unpack z, α, β = p
        n_k  = length(k_grid)
        V    = zeros(n_k)
        G_kp = fill(0,n_k)
        G_c  = zeros(n_k)
        for i = 1:n_k
            V_aux = zeros(n_k) ; # Empty vector for auxiliary value of V(i,j)
            for j = 1:n_k
                # Evaluate potential value function for combinations of
                # current capital k_i and future capital k_j
                V_aux[j] = utility(k_grid[i],k_grid[j],p) + β*V_old[j]
                #println(V_aux[j]," ",k_grid[i]," ",k_grid[j]," ",utility(k_grid[i],k_grid[j],z,a,b) + b*V_old[j])
            end
            # Choose maximum value given current capital k_i
            V[i], G_kp[i] = findmax(V_aux)
            G_c[i]        = z*k_grid[i]^α - k_grid[G_kp[i]]
        end
        return V, G_kp, G_c
    end


# Solve VFI with grid search and loops
    function Solve_VFI_loop(n_k,p::Par)
        # Get Grid
        k_grid = Make_K_Grid(n_k,p)
        # Solve VFI
        V, G_kp, G_c = VFI_grid(x->T_grid_loop(x,k_grid,p),k_grid,p)
        # Return Solution
        return V,G_kp, G_c, k_grid
    end
    #= I am leaving these lines here as an example of how to do VFI with a single funciton
    function VFI_grid_loop(n_k,z,a,b)
        println(" ")
        println("------------------------")
        println("VFI with grid search and loops - n_k=$n_k")
        # Get SS
        k_ss,y_ss,c_ss,r_ss,w_ss = SS_values(z,a,b)
        # Get k_grid
        k_grid   = range(1E-5,2*k_ss;length=n_k) ; # Equally spaced grid between 0 and 2*k_ss
        # Initialize variables for loop
        V_old  = zeros(n_k) ; # Initial value, a vector of zeros
        iter   = 0          ; # Iteration index
        V_dist = 1          ; # Initialize distance
        while iter<=max_iter && V_dist>dist_tol
            # Update value function
            V_new, G_kp, G_c = T_grid_loop(V_old,k_grid,z,a,b)
            # Update distance and iterations
            V_dist = maximum(abs.(V_new./V_old.-1))
            iter  += 1
            # Update old function
            V_old  = V_new
            # Report progress
            if mod(iter,50)==0
                println("   VFI Loop: iter=$iter, dist=",200*V_dist,"%")
            end
        end
        # Check solution
        if (V_dist<=dist_tol)
            # Recover value and policy functions
            V, G_kp, G_c = T_grid_loop(V_old,k_grid,z,a,b)
            # Return
            println("VFI with grid search and loops - Completed - n_k=$n_k")
            println("Iterations = $iter and Distance = ",200*V_dist,"%")
            println("------------------------")
            println(" ")
            return V, G_kp, G_c, k_grid
        else
            error("Error in VFI with loops - Solution not found")
        end
    end
    =#

    # Execute Numerical VFI
    @time V_20, G_kp_20, G_c_20, k_grid_20 = Solve_VFI_loop(20,p)
    @time V_50, G_kp_50, G_c_50, k_grid_50 = Solve_VFI_loop(50,p)
    @time V_200, G_kp_200, G_c_200, k_grid_200 = Solve_VFI_loop(200,p)
    @time V_1000, G_kp_1000, G_c_1000, k_grid_1000 = Solve_VFI_loop(1000,p)

    # Analytical Solutions
    V_20a, G_kp_20a, G_c_20a, Euler_20 = VFI_Analytical_Results(20,G_kp_20,p)
    V_50a, G_kp_50a, G_c_50a, Euler_50 = VFI_Analytical_Results(50,G_kp_50,p)
    V_200a, G_kp_200a, G_c_200a, Euler_200 = VFI_Analytical_Results(200,G_kp_200,p)
    V_1000a, G_kp_1000a, G_c_1000a, Euler_1000 = VFI_Analytical_Results(1000,G_kp_1000,p)

    # Graphs
    VFI_Type = "grid_loop"
    include("./VFI_Graphs.jl")


println("----------------------------------------------------")
println("----------------------------------------------------")
println(" ")


#-----------------------------------------------------------
#-----------------------------------------------------------
# VFI - Grid Search - Matrices
println(" "); println(" "); println(" ")
println("----------------------------------------------------")
println("----------------------------------------------------")
println("Value Function Iteration - Grid Search with Matrices")


# Define function for Value update and policy functions
    function T_grid_mat(V_old,U_mat,k_grid,p::Par)
        @unpack z, α, β = p
        n_k    = length(V_old)
        V,G_kp = findmax( U_mat .+ β*repeat(V_old',n_k,1) , dims=2 )
        G_kp   = [G_kp[i][2] for i in 1:n_k] # G_kp is a Cartesian index
        G_c    = z*k_grid.^α .- k_grid[G_kp]
        return V, G_kp, G_c
    end


# Solve VFI with grid search and loops
    function Solve_VFI_mat(n_k,p::Par)
        # Get Grid
        k_grid = Make_K_Grid(n_k,p)
        # Utility matrix
        U_mat = [utility(k_grid[i],k_grid[j],p) for i in 1:n_k, j in 1:n_k]
        # Solve VFI
        V, G_kp, G_c = VFI_grid(x->T_grid_mat(x,U_mat,k_grid,p),k_grid,p)
        # Return Solution
        return V,G_kp, G_c, k_grid
    end
    #=
    function VFI_grid_mat(n_k,z,a,b)
        println(" ")
        println("------------------------")
        println("VFI with grid search and matrices - n_k=$n_k")
        # Get SS
        k_ss,y_ss,c_ss,r_ss,w_ss = SS_values(z,a,b)
        # Get k_grid
        k_grid   = range(1E-5,2*k_ss;length=n_k) ; # Equally spaced grid between 0 and 2*k_ss
        # Utility matrix
        U_mat = [utility(k_grid[i],k_grid[j],z,a,b) for i in 1:n_k, j in 1:n_k]
        # Initialize variables for loop
        V_old  = zeros(n_k) ; # Initial value, a vector of zeros
        iter   = 0          ; # Iteration index
        V_dist = 1          ; # Initialize distance
        while iter<=max_iter && V_dist>dist_tol
            # Update value function
            V_new, G_kp, G_c = T_grid_mat(V_old,U_mat,k_grid,z,a,b)
            # Update distance and iterations
            V_dist = maximum(abs.(V_new./V_old.-1))
            iter  += 1
            # Update old function
            V_old  = V_new
            # Report progress
            if mod(iter,50)==0
                println("   VFI Loop: iter=$iter, dist=",200*V_dist,"%")
            end
        end
        # Check solution
        if (V_dist<=dist_tol)
            # Recover value and policy functions
            V, G_kp, G_c = T_grid_loop(V_old,k_grid,z,a,b)
            # Return
            println("VFI with grid search and matrices - Completed - n_k=$n_k")
            println("Iterations = $iter and Distance = ",200*V_dist,"%")
            println("------------------------")
            println(" ")
            return V, G_kp, G_c, k_grid
        else
            error("Error in VFI with matrices - Solution not found")
        end
    end
    =#

    # Execute Numerical VFI
    @time V_20, G_kp_20, G_c_20, k_grid_20 = Solve_VFI_mat(20,p)
    @time V_50, G_kp_50, G_c_50, k_grid_50 = Solve_VFI_mat(50,p)
    @time V_200, G_kp_200, G_c_200, k_grid_200 = Solve_VFI_mat(200,p)
    @time V_1000, G_kp_1000, G_c_1000, k_grid_1000 = Solve_VFI_mat(1000,p)

    # Analytical Solutions
    V_20a, G_kp_20a, G_c_20a, Euler_20 = VFI_Analytical_Results(20,G_kp_20,p)
    V_50a, G_kp_50a, G_c_50a, Euler_50 = VFI_Analytical_Results(50,G_kp_50,p)
    V_200a, G_kp_200a, G_c_200a, Euler_200 = VFI_Analytical_Results(200,G_kp_200,p)
    V_1000a, G_kp_1000a, G_c_1000a, Euler_1000 = VFI_Analytical_Results(1000,G_kp_1000,p)

    # Graphs
    VFI_Type = "grid_mat"
    include("./VFI_Graphs.jl")


println("----------------------------------------------------")
println("----------------------------------------------------")
println(" ")


#-----------------------------------------------------------
#-----------------------------------------------------------
# VFI - Grid Search - Matrices - Howard's Policy Iteration
println(" "); println(" "); println(" ")
println("----------------------------------------------------")
println("----------------------------------------------------")
println("Value Function Iteration - Howard's Policy Iteration")


# Define function for Value update and policy functions
    function HT_grid_mat(V_old,U_mat,k_grid,p::Par,n_H)
        @unpack z, α, β, H_tol = p
        # Get Policy Function
        n_k    = length(V_old)
        V,G_kp = findmax( U_mat .+ β*repeat(V_old',n_k,1) , dims=2 )
        V_old  = V
        # "Optimal" U for Howard's iteration
            U_vec = U_mat[G_kp]
        # Howard's policy iteration
        # G_kp is a Cartesian Index
        for i=1:n_H
            V = U_vec .+ β*repeat(V_old',n_k,1)[G_kp]
            if maximum(abs.(V./V_old.-1))<=H_tol
                break
            end
            V_old = V
        end
        # Recover Policy Functions
        G_kp   = [G_kp[i][2] for i in 1:n_k] # G_kp is a Cartesian index
        G_c    = z*k_grid.^α  .- k_grid[G_kp]
        # Return output
        return V, G_kp, G_c
    end


# Solve VFI with Howard's policy iteration
    function Solve_VFI_HPI(n_H,n_k,p::Par)
        # Get Grid
        k_grid = Make_K_Grid(n_k,p)
        # Utility matrix
        U_mat = [utility(k_grid[i],k_grid[j],p) for i in 1:n_k, j in 1:n_k]
        # Solve VFI
        V, G_kp, G_c = VFI_grid(x->HT_grid_mat(x,U_mat,k_grid,p,n_H),k_grid,p)
        # Return Solution
        return V,G_kp, G_c, k_grid
    end
    #=
    function VFI_HPI_grid_mat(n_H,n_k,z,a,b)
        println(" ")
        println("------------------------")
        println("VFI with Howard's Policy Iteration - n_k=$n_k")
        # Get SS
        k_ss,y_ss,c_ss,r_ss,w_ss = SS_values(z,a,b)
        # Get k_grid
        k_grid   = range(1E-5,2*k_ss;length=n_k) ; # Equally spaced grid between 0 and 2*k_ss
        # Utility matrix
        U_mat = [utility(k_grid[i],k_grid[j],z,a,b) for i in 1:n_k, j in 1:n_k]
        # Initialize variables for loop
        V_old  = zeros(n_k) ; # Initial value, a vector of zeros
        iter   = 0          ; # Iteration index
        V_dist = 1          ; # Initialize distance
        while iter<=max_iter && V_dist>dist_tol
            # Update value function
            V_new, G_kp, G_c = HT_grid_mat(V_old,U_mat,k_grid,z,a,b,n_H)
            # Update distance and iterations
            V_dist = maximum(abs.(V_new./V_old.-1))
            iter  += 1
            # Update old function
            V_old  = V_new
            # Report progress
            if mod(iter,50)==0
                println("   VFI Loop: iter=$iter, dist=",200*V_dist,"%")
            end
        end
        # Check solution
        if (V_dist<=dist_tol)
            # Recover value and policy functions
            V, G_kp, G_c = T_grid_loop(V_old,k_grid,z,a,b)
            # Return
            println("VFI with  Howard's Policy Iteration - Completed - n_k=$n_k")
            println("Iterations = $iter and Distance = ",200*V_dist,"%")
            println("------------------------")
            println(" ")
            return V, G_kp, G_c, k_grid
        else
            error("Error in VFI with Howard's Policy Iteration - Solution not found")
        end
    end
    =#

    # Execute Numerical VFI
    @time V_20, G_kp_20, G_c_20, k_grid_20 = Solve_VFI_HPI(20,20,p)
    @time V_50, G_kp_50, G_c_50, k_grid_50 = Solve_VFI_HPI(20,50,p)
    @time V_200, G_kp_200, G_c_200, k_grid_200 = Solve_VFI_HPI(20,200,p)
    @time V_1000, G_kp_1000, G_c_1000, k_grid_1000 = Solve_VFI_HPI(20,1000,p)


    # Analytical Solutions
    V_20a, G_kp_20a, G_c_20a, Euler_20 = VFI_Analytical_Results(20,G_kp_20,p)
    V_50a, G_kp_50a, G_c_50a, Euler_50 = VFI_Analytical_Results(50,G_kp_50,p)
    V_200a, G_kp_200a, G_c_200a, Euler_200 = VFI_Analytical_Results(200,G_kp_200,p)
    V_1000a, G_kp_1000a, G_c_1000a, Euler_1000 = VFI_Analytical_Results(1000,G_kp_1000,p)

    # Graphs
    VFI_Type = "HPI"
    include("./VFI_Graphs.jl")

println("----------------------------------------------------")
println("----------------------------------------------------")
println(" ")





#-----------------------------------------------------------
#-----------------------------------------------------------
# VFI - Grid Search - Matrices - MacQueen Porteus Bounds
println(" "); println(" "); println(" ")
println("----------------------------------------------------")
println("----------------------------------------------------")
println("Value Function Iteration - MacQueen-Porteus Bounds")

# Fixed point with MPB
    function VFI_grid_MPB(T::Function,k_grid,p::Par)
        @unpack β, max_iter, dist_tol = p
        # Initialize variables for loop
        n_k    = length(k_grid) ; # Number of grid nodes
        V_old  = zeros(n_k)     ; # Initial value, a vector of zeros
        iter   = 0              ; # Iteration index
        V_dist = 1              ; # Initialize distance
        println(" ")
        println("------------------------")
        println("VFI - Grid Search - MPB - n_k=$n_k")
        for iter=1:max_iter
            # Update value function
            V_new, G_kp, G_c = T(V_old)
            # MPB and Distance
            MPB_l  = β/(1-β)*minimum(V_new-V_old)
            MPB_h  = β/(1-β)*maximum(V_new-V_old)
            V_dist = MPB_h - MPB_l
            # Update old function
            V_old  = V_new
            # Report progress
            if mod(iter,100)==0
                println("   VFI Loop: iter=$iter, dist=",100*V_dist,"%")
            end
            # Check Convergence
            if (V_dist<=dist_tol)
                # Recover value and policy functions
                V = V_old .+ (MPB_l+MPB_h)/2
                # Return
                println("VFI - Grid Search - MPB - n_k=$n_k")
                println("Iterations = $iter and Distance = ",100*V_dist)
                println("------------------------")
                println(" ")
                return V, G_kp, G_c
            end
        end
        # Report error for non-convergence
        error("Error in VFI - Grid Search - MPB - Solution not found")
    end



# Solve VFI with MPB
    function Solve_VFI_MPB(n_k,p::Par)
        # Get Grid
        k_grid = Make_K_Grid(n_k,p)
        # Utility matrix
        U_mat = [utility(k_grid[i],k_grid[j],p) for i in 1:n_k, j in 1:n_k]
        # Solve VFI
        V, G_kp, G_c = VFI_grid_MPB(x->T_grid_mat(x,U_mat,k_grid,p),k_grid,p)
        # Return Solution
        return V,G_kp, G_c, k_grid
    end

    # Execute Numerical VFI
    @time V_20, G_kp_20, G_c_20, k_grid_20 = Solve_VFI_MPB(20,p)
    @time V_50, G_kp_50, G_c_50, k_grid_50 = Solve_VFI_MPB(50,p)
    @time V_200, G_kp_200, G_c_200, k_grid_200 = Solve_VFI_MPB(200,p)
    @time V_1000, G_kp_1000, G_c_1000, k_grid_1000 = Solve_VFI_MPB(1000,p)


    # Analytical Solutions
    V_20a, G_kp_20a, G_c_20a, Euler_20 = VFI_Analytical_Results(20,G_kp_20,p)
    V_50a, G_kp_50a, G_c_50a, Euler_50 = VFI_Analytical_Results(50,G_kp_50,p)
    V_200a, G_kp_200a, G_c_200a, Euler_200 = VFI_Analytical_Results(200,G_kp_200,p)
    V_1000a, G_kp_1000a, G_c_1000a, Euler_1000 = VFI_Analytical_Results(1000,G_kp_1000,p)

    # Graphs
    VFI_Type = "MPB"
    include("./VFI_Graphs.jl")
