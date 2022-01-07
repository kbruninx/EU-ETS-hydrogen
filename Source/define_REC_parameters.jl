function define_REC_parameters!(REC::Dict,data::Dict,ts::DataFrame,repr_days::DataFrame,scenario_overview_row::DataFrameRow)
    # Number of agents
    REC["nAgents"] = data["nAgents"]

    # RES target
    REC["RT"] = [zeros(3); data["RES_target_2020"]*ones(10); scenario_overview_row["RES_target_2030"]*ones(data["nyears"]-13)]
    
    REC["RS_other_2017"] = data["RES_share_other_2017"]

    return REC
end