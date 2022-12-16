function define_ps_parameters!(mod::Model, data::Dict,ts::DataFrame,repr_days::DataFrame)
    # Parameters 
    mod.ext[:parameters][:VC] = data["fuelCosts"].*[(1+data["YoY_VC"]/100)^(jy-1) for jy in 1:data["nyears"]] # EUR/MWh or MEUR/TWh
    mod.ext[:parameters][:η] = data["efficiency"] # - 
    mod.ext[:parameters][:CI] = data["emissions"] # tCO2/MWh or MtCO2/TWh
    mod.ext[:parameters][:IC] = data["OC"].*[(1+data["YoY_OC"]/100)^(jy-1) for jy in 1:data["nyears"]] # EUR/MW or MEUR/TW    
    mod.ext[:parameters][:DELTA_CAP_MAX] = data["max_YoY_new_cap"]/100 # fraction
    mod.ext[:parameters][:CAP_SV] =  [maximum([0,1-(data["nyears"]-jy+1)/data["Lifetime"]]) for jy=1:data["nyears"]] 
    mod.ext[:parameters][:LEG_CAP] = zeros(data["nyears"],1)
    if haskey(data,"Legcap_2019")
        mod.ext[:parameters][:LEG_CAP][1] = data["AF"]*data["Legcap_2019"]  
    else
        mod.ext[:parameters][:LEG_CAP][1] = data["AF"]*data["Legcap_2021"]  
    end
    if haskey(data,"Legcap_2020")
        mod.ext[:parameters][:LEG_CAP][2] = data["AF"]*data["Legcap_2020"]  
    else
        mod.ext[:parameters][:LEG_CAP][2] = data["AF"]*data["Legcap_2021"]  
    end
    mod.ext[:parameters][:LEG_CAP][3] = data["AF"]*data["Legcap_2021"]  
    mod.ext[:parameters][:LEG_CAP][4:data["nyears"]] = data["AF"]*[data["Legcap_2021"]*maximum([0,(data["Legcap_out"]-jy+1)/data["Legcap_out"]]) for jy=1:data["nyears"]-3]  
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
        mod.ext[:timeseries][:AF] = [ts[!,data["AF_ts"]][round(Int,data["nTimesteps"]*repr_days[!,:periods][jd]+jh)] for jh=1:data["nTimesteps"], jd=1:data["nReprDays"]]  
    else
        mod.ext[:timeseries][:AF] = ones(data["nTimesteps"],data["nReprDays"]) 
    end 
   
    return mod
end