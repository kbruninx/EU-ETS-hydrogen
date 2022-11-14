function define_H2CN_cap_parameters!(H2CN_cap::Dict,data::Dict,ts::DataFrame,repr_days::DataFrame,scenario_overview_row::DataFrameRow)
    # Number of agents
    H2CN_cap["nAgents"] = data["nAgents"]

    # carbon neutral hydrogen production capacity target
    H2CN_cap["H2CN_CAPT"] = [zeros(7); scenario_overview_row["CNH2_cap_target_2024"]*ones(6); scenario_overview_row["CNH2_cap_target_2030"]*ones(data["nyears"]-13)]

    return H2CN_cap
end