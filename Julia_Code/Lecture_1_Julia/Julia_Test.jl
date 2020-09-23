# Test Program for Julia
# Sergio Ocampo
# August 2020

cd("./Dropbox/Teaching/PhD_Macro_Comp/Julia_Code/Lecture_1_Julia")

#-----------------------------------------------------------
# Start with a print statement
println(" ")
println("------------------------")
println("Hello World")
println("   I can type latex in Julia: α β ϵ ♢ ±")

# We can also print and evaluate
println("   I will print some additions")
println("   3+2=",3+2)
println("   3+2.1=",3+2.1)
println("   I can merge strings using *")
println("   String number 1 "*"and "*"String number 2")
println("   I can also include evaluations inside strings with dollar sign")
a, b, c = 2, 3, 2+3
println("   Sum $a+$b=$c or 2+3=$(a+b) and also its type $(typeof(a))")
println("   I can also use the '_' to separate digits, it is ignored:")
println("   1_000_123=",1_000_123)
println("------------------------")
println(" ")

#-----------------------------------------------------------
# We can define logical expressions
println(" ")
println("------------------------")
println("Logicals 'true' and 'false' are recognized by Julia")
println("   ",true)
println("   ",false)
println("The '!' operator is a not")
println("   !true  = ",!true)
println("   !false = ",!false)
println("The '&&' operator is an and")
println("   true && true  = ",true && true)
println("   true && false = ",true && false)
println("The '||' operator is an or")
println("   true || true  = ",true || true)
println("   true || false = ",true || false)
println("Comparing numbers and logicals")
println("   2> 1 = ",2>1)
println("   2< 1 = ",2<1)
println("   2==1 = ",2==1)
println("Julia can compare across types")
println("   2.0== 2 = ",2.0==2)
println("   2.0===2 = ",2.0===2)
println("   The '===' is a strict comparisson operator")
println("Careful not to use bitwise and (&) or (|)")
println("------------------------")
println(" ")


#-----------------------------------------------------------
# Defining arrays (matrices)
println(" ")
println("------------------------")
α = 1
β = 2
A = [α,β]
println("Array A=",A)
B = Array{Float64,2}(undef,2,3)
println("Undefined array B:",B)
B = [1 2 3; 4 5.3 6.1]
println("Assigned array B:",B)
println("------------------------")
println(" ")



#-----------------------------------------------------------
# Plots
using Plots
x = 0:0.5:6
y = cos.(x)
gr()
plot(x,y,linewidth=2,marker=(:diamond,9),markercolor=RGB(0.1,0.1,0.1))
plot(x,y,linetype=:scatter)
savefig("./My_first_fig.pdf")
# See list of attributes with plotattr(:Series), plotattr(:Plot), plotattr(:Axis)


#-----------------------------------------------------------
# Call Outside Script
include("./Outside_Script.jl")
println(" ")
println("------------------------")
println("From Ouside_Script.jl: X=",x," and Y=",y)
println("------------------------")
println(" ")


#=
I can also comment multiple lines
All these lines are commented
Look at this! They are all commented!
=#
