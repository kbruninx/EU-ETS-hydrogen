function define_H2_parameters!(H2::Dict,data::Dict,ts::DataFrame,repr_days::DataFrame,scenario_overview_row::DataFrameRow,H2CN_prod::Dict)
    # Number of agents
    # H2["nAgents"] = data["nAgents"]

    # H2 demand
    H2["D"] = [maximum([data["conv_factor"]*data["demand"],H2CN_prod["H2CN_PRODT"][y]]) for y in 1:data["nyears"]]
    
    return H2
end