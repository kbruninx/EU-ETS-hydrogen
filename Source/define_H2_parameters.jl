function define_H2_parameters!(H2::Dict,data::Dict,ts::DataFrame,repr_days::DataFrame,scenario_overview_row::DataFrameRow,H2CN_prod::Dict)

    # H2 demand
    YoY = (scenario_overview_row["H2_demand_2050"]-data["demand"])/20
    H2["D"] = [[data["conv_factor"]*data["demand"] for y in 1:11]; [data["conv_factor"]*(data["demand"]+YoY*y) for y in 1:20];[data["conv_factor"]*(data["demand"]+YoY*20) for y in 1:data["nyears"]-11-20]]

    return H2
end