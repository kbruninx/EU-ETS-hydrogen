function define_H2S_parameters!(mod::Model, data::Dict,ts::DataFrame,repr_days::DataFrame,REC::Dict)
    # Parameters 
    mod.ext[:parameters][:η_E_H2] = data["efficiency_E_H2"] # - 
    mod.ext[:parameters][:η_NG_H2] = data["efficiency_NG_H2"] # -
    mod.ext[:parameters][:CI] = data["emissions"] # tCO2/MWh or MtCO2/TWh
    mod.ext[:parameters][:IC] = data["OC"].*[(1+data["YoY_OC"]/100)^(jy-1) for jy in 1:data["nyears"]] # EUR/MW or MEUR/TW
    mod.ext[:parameters][:DELTA_CAP_MAX] = data["max_YoY_new_cap"]/100 # fraction
    mod.ext[:parameters][:CAP_SV] =  [maximum([0,1-(data["nyears"]-jy+1)/data["Lifetime"]]) for jy=1:data["nyears"]]
    mod.ext[:parameters][:LEG_CAP] = zeros(data["nyears"],1)
    mod.ext[:parameters][:LEG_CAP][1] = data["AF"]*data["Legcap_2021"]  
    mod.ext[:parameters][:LEG_CAP][2:data["nyears"]] = [data["AF"]*data["Legcap_2021"]*maximum([0,(data["Legcap_out"]-jy+1)/data["Legcap_out"]]) for jy=1:data["nyears"]-1] 
    mod.ext[:parameters][:CAP_LT] = zeros(data["nyears"],data["nyears"]) 
    for y=1:data["nyears"]
        if y+data["Leadtime"] < data["nyears"]
            for yy = y+data["Leadtime"]:minimum([y+data["Leadtime"]+data["Lifetime"]- 1, data["nyears"]])
                mod.ext[:parameters][:CAP_LT][y,yy] = 1
            end
        end
    end
    
    if data["Additionality"] == "Yearly"
        mod.ext[:parameters][:ADD_SF] = ones(data["nyears"],1) 
    else
        mod.ext[:parameters][:ADD_SF] = REC["RT"] # Scaling factor for RECs when additionality is not applied = RES target
    end

    return mod
end