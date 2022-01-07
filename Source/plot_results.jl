function plot_results(ADMM::Dict,results::Dict,agents::Dict,EOM::Dict)
# iterations = 1:(ADMM["n_iter"])
# p1 = plot(iterations,log10.(ADMM["Residuals"]["Primal"]["EOM"][1:ADMM["n_iter"]]),label="Primal EOM",xlabel="Iterations [-]",ylabel="Log10(Residuals)",legend=:outertopright, linewidth = 2);
# plot!(iterations,log10.(ADMM["Residuals"]["Primal"]["ETS"][1:ADMM["n_iter"]]),label="Primal ETS",legend=:outertopright, linewidth = 2);
# plot!(iterations,log10.(ADMM["Residuals"]["Primal"]["REC"][1:ADMM["n_iter"]]),label="Primal REC",legend=:outertopright, linewidth = 2);
# plot!(iterations,log10.(ADMM["Residuals"]["Primal"]["MSR"][1:ADMM["n_iter"]]),label="Primal MSR",legend=:outertopright, linewidth = 2);
# plot!(iterations,log10.(ADMM["Residuals"]["Dual"]["EOM"][1:ADMM["n_iter"]]),label="Dual EOM",legend=:outertopright, linewidth = 2);
# plot!(iterations,log10.(ADMM["Residuals"]["Dual"]["ETS"][1:ADMM["n_iter"]]),label="Dual ETS",legend=:outertopright, linewidth = 2);
# plot!(iterations,log10.(ADMM["Residuals"]["Dual"]["REC"][1:ADMM["n_iter"]]),label="Dual REC",legend=:outertopright, linewidth = 2);

# p2 = plot(iterations,results["λ_EOM"][1:ADMM["n_iter"],1,1,14],label="EOM",xlabel="Iterations [-]",ylabel="Prices (2031)",legend=:outertopright, linewidth = 2);
# plot!(iterations,results["λ_EUA"][1:ADMM["n_iter"],14],label="ETS",legend=:outertopright, linewidth = 2);
# plot!(iterations,results["λ_REC"][1:ADMM["n_iter"],14],label="REC",legend=:outertopright, linewidth = 2);

# p3 = plot(iterations,ADMM["ρ_EOM"][1:ADMM["n_iter"]],label="EOM",xlabel="Iterations [-]",ylabel="ρ-values",legend=:outertopright, linewidth = 2);
# plot!(iterations,ADMM["ρ_EUA"][1:ADMM["n_iter"]],label="ETS",legend=:outertopright, linewidth = 2);
# plot!(iterations,ADMM["ρ_REC"][1:ADMM["n_iter"]],label="REC",legend=:outertopright, linewidth = 2);

# plot(p1, p2, p3, layout = (3, 1))

# Hours = 1:24
# p1 = plot(Hours,EOM["D"][:,1,1],label="D",xlabel="Hours [-]",ylabel="Generation",legend=:outertopright, linewidth = 2);
# for x = 1:length(agents[:ps])
# plot!(Hours,sum(results["g"][m][ADMM["n_iter"],:,1,1] for m in agents[:ps][1:x]),label=agents[:ps][x],xlabel="Hours [-]",ylabel="Generation",legend=:outertopright, linewidth = 2);
# end
# p2 = plot(Hours,results["λ_EOM"][ADMM["n_iter"],:,1,1],label="EOM",xlabel="Hours [-]",ylabel="Prices",legend=:outertopright, linewidth = 2);
# plot(p1,p2,layout = (2,1))


# plot!(size=(600,1000))


# I = mod_ps.ext[:sets][:I]
# IV = mod_ps.ext[:sets][:IV]
# ID = mod_ps.ext[:sets][:ID]
# JH = mod_ps.ext[:sets][:JH]
# JD = mod_ps.ext[:sets][:JD]
# JY = mod_ps.ext[:sets][:JY]
# JY = mod_ps.ext[:sets][:JY]
# W = mod_ps.ext[:parameters][:W]
# GF = mod_ps.ext[:parameters][:GF] # growth factor of the demand
# D = mod_ps.ext[:timeseries][:D] # demand

# cap = value.(mod_ps.ext[:variables][:cap])
# g = value.(mod_ps.ext[:variables][:g])
# ens = value.(mod_ps.ext[:variables][:ens])
# curt = value.(mod_ps.ext[:expressions][:curt])
# λ = dual.(mod_ps.ext[:constraints][:dem_balance])

# e_ps = [value.(mod_ps.ext[:expressions][:e][jy]) for jy in JY]
# e_ind = mod_ind.ext[:parameters][:e]

# demand_yearly = [sum((1+GF[jy])*D[jh,jd] for jh in JH, jd in JD) for jy in JY]
# g_yearly = Dict(i => [sum(W[jd]*g[i,jh,jd,jy] for jh in JH, jd in JD) for jy in JY] for i in I)
# ens_yearly = [sum(W[jd]*ens[jh,jd,jy] for jh in JH, jd in JD) for jy in JY]
# curt_yearly =  [sum(W[jd]*curt[iv,jh,jd,jy] for jh in JH, jd in JD, iv in IV) for jy in JY]
# fuel_shares_yearly = Dict(i => [g_yearly[i][jy]/demand_yearly[jy] for jy in JY] for i in I)
# cap_yearly = Dict(i => [cap[i,jy] for jy in JY] for i in I)

# using Plots
# # Capacity
# p1= plot(JY,zeros(data["nyears"],1),label="",xlabel="Years [-]",ylabel="New capacity [GW]",legend=:outertopright, linewidth = 2);
# for i in I
# plot!(JY,cap_yearly[i],label=i, xlabel="Years [-]",ylabel="New capacity [GW]",legend=:outertopright, linewidth = 2);
# end
# # Generation
# p2= plot(JY,zeros(data["nyears"],1),label="",xlabel="Years [-]",ylabel="Generation [GWh]",legend=:outertopright, linewidth = 2);
# for i in I
# plot!(JY,fuel_shares_yearly[i],label=i, xlabel="Years [-]",ylabel="Generation [GWh]",legend=:outertopright, linewidth = 2);
# end
# plot(p1, p2, layout = (2, 1))
# plot!(size=(600,800))

# EUA prices
# p0= plot(JY,mod_ps.ext[:parameters][:λ_EUA],label="",xlabel="Years [-]",ylabel="EUA prices [euro/tCO2]",legend=:outertopright, linewidth = 2);
# # Emissions
# p1= plot(JY,e_ind,label="Industry",xlabel="Years [-]",ylabel="Emissions [MtCO2]",legend=:outertopright, linewidth = 2);
# plot!(JY,e_ps,label="Power sector",xlabel="Years [-]",ylabel="Emissions [MtCO2]",legend=:outertopright, linewidth = 2);
# # Supply 
# p2= plot(JY,ETS["S"],label="Supply",xlabel="Years [-]",ylabel="EUAs [MtCO2]",legend=:outertopright, linewidth = 2);
# plot!(JY,ETS["CAP"],label="Cap",xlabel="Years [-]",ylabel="EUAs [MtCO2]",legend=:outertopright, linewidth = 2);
# plot(p0, p1, p2, layout = (3, 1))
# plot!(size=(600,800))

end