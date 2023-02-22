function define_H2_parameters!(H2::Dict,data::Dict,ts::DataFrame,repr_days::DataFrame,H2CN_prod::Dict)

    # H2 demand
    YoY_1= (data["H2_demand_2030"]-data["demand"])/9
    YoY_2= (data["H2_demand_2050"]-data["H2_demand_2030"])/20
    H2["D_y"] = [[data["conv_factor"]*(data["demand"]+YoY_1*(y-1)) for y in 1:10]; [data["conv_factor"]*(data["H2_demand_2030"]+YoY_2*y) for y in 1:20];[data["conv_factor"]*(data["H2_demand_2030"]+YoY_2*20) for y in 1:data["nyears"]-30]]
    H2["D_h"] = [H2["D_y"][y]/8760 for h in 1:data["nTimesteps"], d in 1:data["nReprDays"], y in 1:data["nyears"]] 
    H2["D_d"] = [H2["D_y"][y]/365 for d in 1:data["nReprDays"], y in 1:data["nyears"]] 
    H2["D_m"] = [H2["D_y"][y]/12 for d in 1:data["nMonths"], y in 1:data["nyears"]] 

    # Weighted demand
    H2["W"] = W = Dict(jd => repr_days[!,:weights][jd] for jd=1:data["nReprDays"])
    H2["Dw"] = [W[jd]*H2["D_h"][jh,jd,jy] for jh=1:data["nTimesteps"], jd=1:data["nReprDays"], jy=1:data["nyears"]]
    return H2
end