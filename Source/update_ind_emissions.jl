function update_ind_emissions!(mod::Model,data::Dict,ETS::Dict)
    # Baseline emissions, corrected for share of industry in emissions 
    E_REF = data["E_ref"]*ones(data["nyears"],1)
    E_REF[4:end] = E_REF[4:end] + data["E_ref_maritime"]*ones(data["nyears"]-3,1)

    for y = 2:data["nyears"] # Emissions are only price-dependent as of 2022 (y=2)
        λ_nom = maximum([0,mod.ext[:parameters][:λ_EUA][y]/(1+data["inflation"])^(y-1)]) # M€/MtCO2, discounted to 2021 values, limited to positive values
        mod.ext[:parameters][:e][y] = minimum([E_REF[y],maximum([0,E_REF[y] - (λ_nom/mod.ext[:parameters][:β])^(1/data["gamma"])])]) # emissions according to MACC

        # Account for maximum price-induced change in emission allowance demand (as percentage of 2021 emissions):
        if mod.ext[:parameters][:e][y-1] > data["max_em_change"]*mod.ext[:parameters][:e][1]/100
            if mod.ext[:parameters][:e][y-1] - mod.ext[:parameters][:e][y] > data["max_em_change"]*mod.ext[:parameters][:e][1]/100
                mod.ext[:parameters][:e][y] = mod.ext[:parameters][:e][y-1]-data["max_em_change"]*mod.ext[:parameters][:e][1]/100
            end
        end
        
        # Compute abatement cost
        mod.ext[:parameters][:AC][y] = mod.ext[:parameters][:β]*(E_REF[y]-mod.ext[:parameters][:e][y])^(data["gamma"]+1)/(data["gamma"]+1)

        # Add covid or overlapping policy
        mod.ext[:parameters][:e][y] = mod.ext[:parameters][:e][y]+ETS["Δe"][y]
    end

    return mod
end