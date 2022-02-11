function define_ps_parameters!(mod::Model, data::Dict,ts::DataFrame,repr_days::DataFrame,scenario_overview_row::DataFrameRow)
    # Parameters 
    mod.ext[:parameters][:W] = Dict(jd => repr_days[!,:Period_weight][jd] for jd=1:data["nReprDays"]) # weights of each representative day
    mod.ext[:parameters][:VC] = data["fuelCosts"].*[(1+data["YoY_VC"]/100)^(jy-1) for jy in 1:data["nyears"]] # EUR/MWh or MEUR/TWh
    mod.ext[:parameters][:CI] = data["emissions"] # tCO2/MWh or MtCO2/TWh
    mod.ext[:parameters][:IC] = data["OC"].*[(1+data["YoY_OC"]/100)^(jy-1) for jy in 1:data["nyears"]] # EUR/MW or MEUR/TW    
    mod.ext[:parameters][:DELTA_CAP_MAX] = data["max_YoY_new_cap"] # TW
    mod.ext[:parameters][:CAP_SV] =  [maximum([0,1-(data["nyears"]-jy+1)/data["Lifetime"]]) for jy=1:data["nyears"]] 
    mod.ext[:parameters][:LEG_CAP] = data["AF"]*[data["Legcap"]*maximum([0,(data["Legcap_out"]-jy+1)/data["Legcap_out"]]) for jy=1:data["nyears"]]  
    mod.ext[:parameters][:CAP_LT] = zeros(data["nyears"],data["nyears"]) 
    for y=1:data["nyears"]
        if y+data["Leadtime"] < data["nyears"]
            for yy = y+data["Leadtime"]:minimum([y+data["Leadtime"]+data["Lifetime"]-1,data["nyears"]])
                mod.ext[:parameters][:CAP_LT][y,yy] = 1
            end
        end
    end

   # Availability factors
    if data["AF_ts"] != "NA" 
        mod.ext[:timeseries][:AF] = [ts[!,data["AF_ts"]][round(Int,data["nTimesteps"]*repr_days[!,:Period_index][jd]+jh)]/10^3/data["Legcap"] for jh=1:data["nTimesteps"], jd=1:data["nReprDays"]] # scaling: from MW to GW
    else
        mod.ext[:timeseries][:AF] = ones(data["nTimesteps"],data["nReprDays"]) 
    end 
   
    return mod
end