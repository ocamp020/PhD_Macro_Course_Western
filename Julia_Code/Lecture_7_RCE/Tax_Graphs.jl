


using Plots
cd("./Dropbox/Teaching/PhD_Macro_Comp/Julia_Code/Lecture_7_RCE/")
mkpath("Figures")


y_bt=range(1,100,length=1000)

 y_at(τ,θ,y_bar)=(1-τ)*y_bt.^(1-θ).+y_bar
dy_at(τ,θ,y_bar)=(1-θ)*(1-τ)*y_bt.^(-θ)

y_at_vec  =  y_at(0.25,0.1,0)
dy_at_vec = dy_at(0.25,0.1,0)

gr()
plot(title="Progressive Taxes: τ=0.25, θ=0.1",legend=:bottomright,foreground_color_legend = nothing,background_color_legend = nothing)
plot!(y_bt,y_bt.-y_at_vec             ,linewidth=2.5,label="Tax Level")
plot!(y_bt,100*(1 .- dy_at_vec)       ,linewidth=2.5,label="Marginal Tax (in pp)")
plot!(y_bt,100*(1 .- y_at_vec./y_bt)  ,linewidth=2.5,label="Average Tax (in pp)")
xlabel!("Before Tax Income")
savefig("./Figures/Progressive_Tax.pdf")

y_at_vec  =  y_at(0.25,-0.05,0)
dy_at_vec = dy_at(0.25,-0.05,0)

gr()
plot(title="Regressive Taxes: τ=0.25, θ=-0.05",legend=:topright,foreground_color_legend = nothing,background_color_legend = nothing)
plot!(y_bt,y_bt.-y_at_vec             ,linewidth=2.5,label="Tax Level")
plot!(y_bt,100*(1 .- dy_at_vec)       ,linewidth=2.5,label="Marginal Tax (in pp)")
plot!(y_bt,100*(1 .- y_at_vec./y_bt)  ,linewidth=2.5,label="Average Tax (in pp)")
xlabel!("Before Tax Income")
savefig("./Figures/Regressive_Tax_Transfer.pdf")


y_at_vec  =  y_at(0.25,0,0)
dy_at_vec = dy_at(0.25,0,0)

gr()
plot(title="Linear Taxes: τ=0.25, θ=0",legend=:topright,foreground_color_legend = nothing,background_color_legend = nothing)
plot!(y_bt,y_bt.-y_at_vec             ,linewidth=2.5,label="Tax Level")
plot!(y_bt,100*(1 .- dy_at_vec)       ,linewidth=2.5,label="Marginal Tax (in pp)")
plot!(y_bt,100*(1 .- y_at_vec./y_bt)  ,linewidth=2.5,label="Average Tax (in pp)")
xlabel!("Before Tax Income")
savefig("./Figures/Linear_Tax_Transfer.pdf")

y_at_vec  =  y_at(0.25,0.1,2)
dy_at_vec = dy_at(0.25,0.1,2)

gr()
plot(title="Progressive Taxes with Transfer: τ=0.25, θ=0.1, y_bar=2",legend=:bottomright,foreground_color_legend = nothing,background_color_legend = nothing)
plot!(y_bt,y_bt.-y_at_vec             ,linewidth=2.5,label="Tax Level")
plot!(y_bt,100*(1 .- dy_at_vec)       ,linewidth=2.5,label="Marginal Tax (in pp)")
plot!(y_bt,100*(1 .- y_at_vec./y_bt)  ,linewidth=2.5,label="Average Tax (in pp)")
xlabel!("Before Tax Income")
savefig("./Figures/Progressive_Tax_Transfer.pdf")
