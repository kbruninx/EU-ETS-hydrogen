function define_H2CN_prod_parameters!(H2CN_prod::Dict,data::Dict,ts::DataFrame,repr_days::DataFrame,scenario_overview_row::DataFrameRow)
    # H2 demand
    H2CN_prod["H2CN_PRODT"] = [zeros(5); data["conv_factor"]*scenario_overview_row["CNH2_demand_2024"]*ones(6); data["conv_factor"]*scenario_overview_row["CNH2_demand_2030"]*ones(data["nyears"]-11)]

    return H2CN_prod
end