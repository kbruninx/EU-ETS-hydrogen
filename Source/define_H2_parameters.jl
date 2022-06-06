function define_H2_parameters!(H2::Dict,data::Dict,ts::DataFrame,repr_days::DataFrame,scenario_overview_row::DataFrameRow)
    # Number of agents
    H2["nAgents"] = data["nAgents"]

    # H2 demand
    H2["D"] = data["demand"]*ones(data["nyears"])
    
    return H2
end