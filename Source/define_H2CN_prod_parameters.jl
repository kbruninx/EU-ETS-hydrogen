function define_H2CN_prod_parameters!(H2CN_prod::Dict,data::Dict,ts::DataFrame,repr_days::DataFrame,scenario_overview_row::DataFrameRow)
    # Number of agents
    H2CN_prod["nAgents"] = data["nAgents"]

    # H2 demand
    H2CN_prod["H2CN_PRODT"] = [zeros(7); data["H2CN_prod_target_2024"]*ones(6); data["H2CN_prod_target_2030"]*ones(data["nyears"]-13)]

    return H2CN_prod
end