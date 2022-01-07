function update_ind_emissions!(mod::Model,data::Dict,ETS::Dict,scenario_overview_row::DataFrameRow)
    # Baseline emissions, corrected for share of industry in emissions, not accounting for difference in scope (Fit for 55 vs. current rules)
    E_REF = data["E_ref"]*ones(data["nyears"],1)*data["SF_ind"]
  
    for y = 4:data["nyears"] # Emissions are only price-dependent as of 2020 (y=4)
        λ_nom = maximum([0,mod.ext[:parameters][:λ_EUA][y]/(1+data["inflation"])^(y-3)]) # M€/MtCO2, discounted to 2019 values, limited to positive values
        mod.ext[:parameters][:e][y] = minimum([E_REF[y],maximum([0,E_REF[y] - (λ_nom/mod.ext[:parameters][:β])^(1/scenario_overview_row[:gamma])])]) # emissions according to MACC

        # Account for maximum price-induced change in emission allowance demand (as percentage of 2019 emissions):
        if mod.ext[:parameters][:e][y-1] > scenario_overview_row[:max_em_change]*mod.ext[:parameters][:e][3]/100
            if mod.ext[:parameters][:e][y-1] - mod.ext[:parameters][:e][y] > scenario_overview_row[:max_em_change]*mod.ext[:parameters][:e][3]/100
                mod.ext[:parameters][:e][y] = mod.ext[:parameters][:e][y-1]-scenario_overview_row[:max_em_change]*mod.ext[:parameters][:e][3]/100
            end
        end

        # Add covid or overlapping policy
        mod.ext[:parameters][:e][y] = mod.ext[:parameters][:e][y]+ETS["Δe"][y]
    end

    return mod
end