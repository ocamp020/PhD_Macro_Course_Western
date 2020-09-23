# Interpolation for Julia
# Sergio Ocampo
# September 2020
cd() # Go to root directory
cd("./Dropbox/Teaching/PhD_Macro_Comp/Julia_Code/Lecture_3_Interpolation/")
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
    # Define Runge's function
    f(x)  =1/(1+x^2)
    df(x) = -2*x*(1+x^2)^-2
    # Fine Grid
    xaxis = range(-5,5;length=1000) ;
    # Exact evaluation of Runge's function
    runge            = map(f,xaxis) # Runge's function values
    runge_derivative = map(df,xaxis) # Runge's function derivative values


    #-----------------------------------------------------------
    #-----------------------------------------------------------
    # Polynomial Interpolation: Newton Basis
    # Code from Giray Otken, First semester in numerical analysis with Julia

        function diff(x::Array,y::Array)
            m = length(x) #here m is the number of data points. #the degree of the polynomial is m-1 a=Array{Float64}(undef,m)
            a = zeros(m)
            for i in 1:m
                a[i]=y[i]
            end
            for j in 2:m
                for i in reverse(collect(j:m))
                    a[i]=(a[i]-a[i-1])/(x[i]-x[i-(j-1)])
                end
            end
            return(a)
        end

        function newton(x::Array,y::Array,z)
            m=length(x) #here m is the number of data points, not the degree # of the polynomial
            a=diff(x,y)
            sum=a[1]
            pr=1.0
            for j in 1:(m-1)
                pr=pr*(z-x[j])
                sum=sum+a[j+1]*pr
            end
            return sum
        end

#-----------------------------------------------------------
#-----------------------------------------------------------
# Runge's Example with Newton Polynomial for different n
for n = (4,6,11,21)
    # Grid of nodes for interpolation
    xi = collect(range(-5,5;length=n)) ; # Collect makes it an array instead of a collection
    yi = map(f,xi) # the corresponding y-coordinates
    # Interpolation
    interp=map(z->newton(xi,yi,z),xaxis) # Interpolating poly for the data

    # Plot
    gr()
    plt = plot(title="Interpolation n=$n - Newton Polynomial")
    plot!(xaxis,runge,linewidth=3,label = "Runge's Function",legend=(0.75,0.75),foreground_color_legend = nothing,background_color_legend = nothing)
    plot!(xaxis,interp,linewidth=3,label="Interpolation")
    plot!(xi,yi,linetype=:scatter,marker=(:diamond,9),markercolor=RGB(0.5,0.1,0.1),label = "Data")
    savefig("./Figures/Runge_Newton_n_$n.pdf")
end


#-----------------------------------------------------------
#-----------------------------------------------------------
# Runge's Example with Linear Spline
for n = (4,6,11,21)
    # Grid of nodes for interpolation
    xi = range(-5,5;length=n) ;
    yi = map(f,xi) # the corresponding y-coordinates
    # Interpolation
    interp = interpolate((xi,),yi, Gridded(Linear()))
    # interp = Spline1D(xi, yi; k=1, bc="nearest")
    # Plot
    gr()
    plt = plot(title="Interpolation n=$n - Linear Spline")
    plot!(xaxis,runge,linewidth=3,label = "Runge's Function",legend=(0.75,0.75),foreground_color_legend = nothing,background_color_legend = nothing)
    plot!(xaxis,interp(xaxis),linewidth=3,label="Interpolation")
    plot!(xi,yi,linetype=:scatter,marker=(:diamond,9),markercolor=RGB(0.5,0.1,0.1),label = "Data")
    savefig("./Figures/Runge_LS_n_$n.pdf")
end



#-----------------------------------------------------------
#-----------------------------------------------------------
# Runge's Example with Cubic Spline
for n = (4,6,11,21)
    # Grid of nodes for interpolation
    xi = range(-5,5;length=n) ;
    yi = map(f,xi) # the corresponding y-coordinates
    d_yi = map(df,xi) # The corresponding derivatives
    # Interpolation (interpolates.jl)
    interp = interpolate(yi, BSpline(Cubic(Line(OnGrid()))))
        # interpolate assumes grid is equally spaced in [0,1]
        interp   = scale(interp, xi)
        # Derivative (gradient computes one value at a time)
        deriv_ip(x) = Interpolations.gradient(interp,x)[1]

    # Interpolation (Dierckx.jl)
    # interp = Spline1D(xi, yi; k=3, bc="nearest")
    #     # Derivative
    #     deriv_ip(x) = derivative(interp, x)

    # Derivative using ForwardDiff
        deriv_ip_FD(x) = ForwardDiff.derivative(interp,x)

    # Plot
    gr()
    plt = plot(title="Interpolation n=$n - Cubic Spline")
    plot!(xaxis,runge,linewidth=3,label = "Runge's Function",legend=(0.75,0.75),foreground_color_legend = nothing,background_color_legend = nothing)
    plot!(xaxis,interp(xaxis),linewidth=3,label="Interpolation")
    plot!(xi,yi,linetype=:scatter,marker=(:diamond,9),markercolor=RGB(0.5,0.1,0.1),label = "Data")
    savefig("./Figures/Runge_CS_n_$n.pdf")

    gr()
    plt = plot(title="Derivative Interpolation n=$n - Cubic Spline")
    plot!(xaxis,runge_derivative,linewidth=3,label = "Runge's Derivative",legend=(0.75,0.75),foreground_color_legend = nothing,background_color_legend = nothing)
    plot!(xaxis,deriv_ip.(xaxis),linewidth=3,label="Interpolation")
    plot!(xaxis,deriv_ip_FD.(xaxis),linewidth=2.5,linestyle=(:dash),label="Interpolation FD")
    plot!(xi,d_yi,linetype=:scatter,marker=(:diamond,9),markercolor=RGB(0.5,0.1,0.1),label = "Data")
    savefig("./Figures/Runge_CS_Derivative_n_$n.pdf")
end


#-----------------------------------------------------------
#-----------------------------------------------------------
# Runge's Example with Monotone Perserving Cubic Spline
for n = (4,6,11,21)
    # Grid of nodes for interpolation
    xi = range(-5,5;length=n) ;
        # xi_u = (xi .- xi[1])./(xi[end]-xi[1])
    yi = map(f,xi) # the corresponding y-coordinates
    d_yi = map(df,xi) # The corresponding derivatives
    # Interpolation
    # interp = interpolate(xi, yi, FritschCarlsonMonotonicInterpolation())
     interp = interpolate(xi, yi, FritschButlandMonotonicInterpolation())
        # Monotone interpolation does apply to vectors, only to numbers
    # Derivative (gradient computes one value at a time)
        deriv_ip(x) = Interpolations.gradient(interp,x)[1]
    # Derivative using ForwardDiff
        deriv_ip_FD(x) = ForwardDiff.derivative(interp,x)

    # Plot
    gr()
    plt = plot(title="Interpolation n=$n - Monotone Cubic Spline")
    plot!(xaxis,runge,linewidth=3,label = "Runge's Function",legend=(0.75,0.75),foreground_color_legend = nothing,background_color_legend = nothing)
    plot!(xaxis,interp.(xaxis),linewidth=3,label="Interpolation")
    plot!(xi,yi,linetype=:scatter,marker=(:diamond,9),markercolor=RGB(0.5,0.1,0.1),label = "Data")
    savefig("./Figures/Runge_MS_n_$n.pdf")

    gr()
    plt = plot(title="Derivative Interpolation n=$n - Monotone Cubic Spline")
    plot!(xaxis,runge_derivative,linewidth=3,label = "Runge's Derivative",legend=(0.75,0.75),foreground_color_legend = nothing,background_color_legend = nothing)
    plot!(xaxis,deriv_ip.(xaxis),linewidth=3,label="Interpolation")
    plot!(xaxis,deriv_ip_FD.(xaxis),linewidth=2.5,linestyle=(:dash),label="Interpolation FD")
    plot!(xi,d_yi,linetype=:scatter,marker=(:diamond,9),markercolor=RGB(0.5,0.1,0.1),label = "Data")
    savefig("./Figures/Runge_MS_Derivative_n_$n.pdf")
end



#-----------------------------------------------------------
#-----------------------------------------------------------
# CRRA Example U(c)=c^(1-γ)/(1-γ)
    γ     = 10
    U(x)  = x^(1-γ)/(1-γ)
    dU(x) = x^-γ

    # Fine Grid
    xaxis = range(0.01,0.50;length=1000) ;
    # Exact evaluation of Runge's function
    Utility            = map(U,xaxis) # Utility function values
    Utility_derivative = map(dU,xaxis) # Utility derivative values

#-----------------------------------------------------------
#-----------------------------------------------------------
# CRRA Example U(c)=c^(1-γ)/(1-γ)

for n = (5,10,50,100)
    # Grid of nodes for interpolation
        ci = range(0.01,0.50;length=n) ;
        ui = map(U,ci) # the corresponding y-coordinates
        d_ui = map(dU,ci) # The corresponding derivatives
    # Interpolation Cubic Splines (y(x_0)''=0)
    U_ip_CS_nat = interpolate(ui, BSpline(Cubic(Line(OnGrid()))))
        # interpolate assumes grid is equally spaced in [0,1]
        U_ip_CS_nat   = scale(U_ip_CS_nat, ci)
        # Derivative
        dU_ip_CS_nat(x) = Interpolations.gradient(U_ip_CS_nat,x)[1]
    # Interpolation Cubic Splines (y(x_0)'=0)
    U_ip_CS_flat = interpolate(ui, BSpline(Cubic(Flat(OnGrid()))))
        # interpolate assumes grid is equally spaced in [0,1]
        U_ip_CS_flat   = scale(U_ip_CS_flat, ci)
        # Derivative
        dU_ip_CS_flat(x) = Interpolations.gradient(U_ip_CS_flat,x)[1]
    # Interpolation Monotone Spline
        U_ip_MS = interpolate(ci, ui, FritschButlandMonotonicInterpolation())
        dU_ip_MS(x) = Interpolations.gradient(U_ip_MS,x)[1]
    # Interpolation Newton
        U_ip_NP=map(z->newton(collect(ci),ui,z),collect(xaxis)) # Interpolating poly for the data

    # Plot
    gr()
    plt = plot(title="Interpolation n=$n - Splines")
    plot!(log.(xaxis),Utility,linewidth=3,label = "Utility Function",legend=(0.75,0.75),foreground_color_legend = nothing,background_color_legend = nothing)
    plot!(log.(xaxis),U_ip_CS_nat(xaxis),linewidth=3,label="Interpolation CS-Natural")
    plot!(log.(xaxis),U_ip_CS_flat(xaxis),linewidth=3,linestyle=(:dot),label="Interpolation CS-Flat")
    plot!(log.(xaxis),U_ip_MS.(xaxis),linewidth=3,linestyle=(:dash),label="Interpolation MS")
    # plot!(log.(xaxis),U_ip_NP,linewidth=3,linestyle=(:dashdotdot),label="Interpolation NP")
    plot!(log.(ci),ui,linetype=:scatter,marker=(:diamond,9),markercolor=RGB(0.5,0.1,0.1),label = "Data")
    xlabel!("Log(Consumption)")
    savefig("./Figures/CRRA_n_$n.pdf")

    plot!(title="Interpolation n=$n - Cubic Spline",xlims = (log(0.01),log(0.15)))
    savefig("./Figures/CRRA_zoom_n_$n.pdf")

    gr()
    plt = plot(title="Interpolation n=$n - Splines")
    plot!(log.(xaxis),Utility_derivative,linewidth=3,label = "Derivative of Utility",legend=(0.75,0.75),foreground_color_legend = nothing,background_color_legend = nothing)
    plot!(log.(xaxis),dU_ip_CS_nat.(xaxis),linewidth=3,label="Interpolation CS-Natural")
    plot!(log.(xaxis),dU_ip_CS_flat.(xaxis),linewidth=3,linestyle=(:dot),label="Interpolation CS-Flat")
    plot!(log.(xaxis),dU_ip_MS.(xaxis),linewidth=3,linestyle=(:dash),label="Interpolation MS")
    # plot!(log.(xaxis),U_ip_NP,linewidth=3,linestyle=(:dashdotdot),label="Interpolation NP")
    plot!(log.(ci),d_ui,linetype=:scatter,marker=(:diamond,9),markercolor=RGB(0.5,0.1,0.1),label = "Data")
    xlabel!("Log(Consumption)")
    savefig("./Figures/CRRA_Derivative_n_$n.pdf")

    plot!(title="Interpolation n=$n - Cubic Spline",xlims = (log(0.01),log(0.15)))
    savefig("./Figures/CRRA_Derivative_zoom_n_$n.pdf")
end


#-----------------------------------------------------------
#-----------------------------------------------------------
# CRRA Example U(c)=c^(1-γ)/(1-γ) - Curved Grid

# Call Scaled Interpolation Functions
include("./Scaled_Interpolation_Functions.jl")

for n = (5,10,50,100)
    # Grid of nodes for interpolation
        ci = PolyRange(0.01,0.50;θ=2,N=n)
        ui = map(U,ci) # the corresponding y-coordinates
        d_ui = map(dU,ci) # The corresponding derivatives
    # Interpolation Cubic Splines (y(x_0)''=0)
        U_ip_CS_nat = ScaledInterpolations(ci,ui, BSpline(Cubic(Line(OnGrid()))))
        dU_ip_CS_nat(x) = ForwardDiff.derivative(U_ip_CS_nat,x)
    # Interpolation Cubic Splines (y(x_0)'=0)
        U_ip_CS_flat = ScaledInterpolations(ci,ui, BSpline(Cubic(Flat(OnGrid()))))
        dU_ip_CS_flat(x) = ForwardDiff.derivative(U_ip_CS_flat,x)
    # Interpolation Monotone Splines
        U_ip_MS = ScaledInterpolations(ci, ui, FritschButlandMonotonicInterpolation())
        dU_ip_MS(x) = ForwardDiff.derivative(U_ip_MS,x)
        # # No need for ScaledInterpolations if not using PolyRange
        #     ci = 0.01 .+ (0.50-0.01).*range(0,1;length=n).^2
        #     U_ip_MS = interpolations(ci, ui, FritschButlandMonotonicInterpolation())
    # Interpolation Newton
        U_ip_NP=map(z->newton(collect(ci),ui,z),collect(xaxis)) # Interpolating poly for the data

    # Plot
    gr()
    plt = plot(title="Interpolation n=$n - Splines")
    plot!(log.(xaxis),Utility,linewidth=3,label = "Utility Function",legend=(0.75,0.75),foreground_color_legend = nothing,background_color_legend = nothing)
    plot!(log.(xaxis),U_ip_CS_nat.(xaxis),linewidth=3,label="Interpolation CS-Natural")
    plot!(log.(xaxis),U_ip_CS_flat.(xaxis),linewidth=3,linestyle=(:dot),label="Interpolation CS-Flat")
    plot!(log.(xaxis),U_ip_MS.(xaxis),linewidth=3,linestyle=(:dash),label="Interpolation MS")
    #plot!(log.(xaxis),U_ip_NP,linewidth=3,linestyle=(:dashdotdot),label="Interpolation NP")
    plot!(log.(ci),ui,linetype=:scatter,marker=(:diamond,9),markercolor=RGB(0.5,0.1,0.1),label = "Data")
    xlabel!("Log(Consumption)")
    savefig("./Figures/Curved_CRRA_n_$n.pdf")

    plot!(title="Interpolation n=$n - Cubic Spline",xlims = (log(0.01),log(0.15)))
    savefig("./Figures/Curved_CRRA_Derivative_zoom_n_$n.pdf")

    gr()
    plt = plot(title="Interpolation n=$n - Splines")
    plot!(log.(xaxis),Utility_derivative,linewidth=3,label = "Derivative of Utility",legend=(0.75,0.75),foreground_color_legend = nothing,background_color_legend = nothing)
    plot!(log.(xaxis),dU_ip_CS_nat.(xaxis),linewidth=3,label="Interpolation CS-Natural")
    plot!(log.(xaxis),dU_ip_CS_flat.(xaxis),linewidth=3,linestyle=(:dot),label="Interpolation CS-Flat")
    plot!(log.(xaxis),dU_ip_MS.(xaxis),linewidth=3,linestyle=(:dash),label="Interpolation MS")
    #plot!(log.(xaxis),U_ip_NP,linewidth=3,linestyle=(:dashdotdot),label="Interpolation NP")
    plot!(log.(ci),d_ui,linetype=:scatter,marker=(:diamond,9),markercolor=RGB(0.5,0.1,0.1),label = "Data")
    xlabel!("Log(Consumption)")
    savefig("./Figures/Curved_CRRA_Derivative_n_$n.pdf")

    plot!(title="Interpolation n=$n - Cubic Spline",xlims = (log(0.01),log(0.15)))
    savefig("./Figures/Curved_CRRA_Derivative_zoom_n_$n.pdf")
end

for t = (1,1.5,2,3)
    # Grid size
        n = 10
    # Grid of nodes for interpolation
        ci =PolyRange(0.01,0.50;θ=t,N=n)
        ui = map(U,ci) # the corresponding y-coordinates
        d_ui = map(dU,ci) # The corresponding derivatives
    # Interpolation Cubic Splines (y(x_0)''=0)
        U_ip_CS_nat = ScaledInterpolations(ci,ui, BSpline(Cubic(Line(OnGrid()))))
        dU_ip_CS_nat(x) = ForwardDiff.derivative(U_ip_CS_nat,x)
    # Interpolation Cubic Splines (y(x_0)'=0)
        U_ip_CS_flat = ScaledInterpolations(ci,ui, BSpline(Cubic(Flat(OnGrid()))))
        dU_ip_CS_flat(x) = ForwardDiff.derivative(U_ip_CS_flat,x)
    # Interpolation Monotone Splines
        U_ip_MS = ScaledInterpolations(ci, ui, FritschButlandMonotonicInterpolation())
        dU_ip_MS(x) = ForwardDiff.derivative(U_ip_MS,x)
        # # No need for ScaledInterpolations if not using PolyRange
        #     ci = 0.01 .+ (0.50-0.01).*range(0,1;length=n).^2
        #     U_ip_MS = interpolations(ci, ui, FritschButlandMonotonicInterpolation())
    # Interpolation Newton
        U_ip_NP=map(z->newton(collect(ci),ui,z),collect(xaxis)) # Interpolating poly for the data

    # Plot
    gr()
    plt = plot(title="Interpolation n=$n - θ=$t")
    plot!(log.(xaxis),Utility,linewidth=3,label = "Utility Function",legend=(0.75,0.75),foreground_color_legend = nothing,background_color_legend = nothing)
    plot!(log.(xaxis),U_ip_CS_nat.(xaxis),linewidth=3,label="Interpolation CS-Natural")
    plot!(log.(xaxis),U_ip_CS_flat.(xaxis),linewidth=3,linestyle=(:dot),label="Interpolation CS-Flat")
    plot!(log.(xaxis),U_ip_MS.(xaxis),linewidth=3,linestyle=(:dash),label="Interpolation MS")
    #plot!(log.(xaxis),U_ip_NP,linewidth=3,linestyle=(:dashdotdot),label="Interpolation NP")
    plot!(log.(ci),ui,linetype=:scatter,marker=(:diamond,9),markercolor=RGB(0.5,0.1,0.1),label = "Data")
    xlabel!("Log(Consumption)")
    savefig("./Figures/Curved_CRRA_t_$t.pdf")

    plot!(title="Interpolation n=$n - Cubic Spline",xlims = (log(0.01),log(0.15)))
    savefig("./Figures/Curved_CRRA_zoom_t_$t.pdf")

    gr()
    plt = plot(title="Interpolation n=$n - θ=$t")
    plot!(log.(xaxis),Utility_derivative,linewidth=3,label = "Derivative of Utility",legend=(0.75,0.75),foreground_color_legend = nothing,background_color_legend = nothing)
    plot!(log.(xaxis),dU_ip_CS_nat.(xaxis),linewidth=3,label="Interpolation CS-Natural")
    plot!(log.(xaxis),dU_ip_CS_flat.(xaxis),linewidth=3,linestyle=(:dot),label="Interpolation CS-Flat")
    plot!(log.(xaxis),dU_ip_MS.(xaxis),linewidth=3,linestyle=(:dash),label="Interpolation MS")
    #plot!(log.(xaxis),U_ip_NP,linewidth=3,linestyle=(:dashdotdot),label="Interpolation NP")
    plot!(log.(ci),d_ui,linetype=:scatter,marker=(:diamond,9),markercolor=RGB(0.5,0.1,0.1),label = "Data")
    xlabel!("Log(Consumption)")
    savefig("./Figures/Curved_CRRA_Derivative_t_$t.pdf")

    plot!(title="Interpolation n=$n - Cubic Spline",xlims = (log(0.01),log(0.15)))
    savefig("./Figures/Curved_CRRA_Derivative_zoom_t_$t.pdf")
end

# Plot Grid
marker_list = (:diamond,:circle,:star4)
gr()
p = plot(legend=(0.75,0.95),foreground_color_legend = nothing,background_color_legend = nothing)
for t=(1,2,3)
    n    = 10
    grid = PolyRange(0.0,1.00;θ=t,N=n)
    plot!(grid,t.*ones(n),linewidth=2,label = "Grid, θ=$t",marker=(marker_list[t],7),markercolor=RGB(0.3,0.3,0.3),linecolor=RGB(0.6,0.6,0.6))
end
    xlabel!("Grid"); ylabel!("Curvature θ"); ylims!(0,4)
    savefig("./Figures/Curved_Grid.pdf")
