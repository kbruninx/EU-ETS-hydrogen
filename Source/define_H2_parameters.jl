function define_H2_parameters!(H2::Dict,data::Dict,ts::DataFrame,repr_days::DataFrame,scenario_overview_row::DataFrameRow,H2CN_prod::Dict)

    # H2 demand
    YoY_1= (scenario_overview_row["H2_demand_2030"]-data["demand"])/9
    YoY_2= (scenario_overview_row["H2_demand_2050"]-scenario_overview_row["H2_demand_2030"])/20

    H2["D"] =[data["conv_factor"]*data["demand"];data["conv_factor"]*data["demand"];data["conv_factor"]*data["demand"];[data["conv_factor"]*(data["demand"]+YoY_1*y) for y in 1:8]; [data["conv_factor"]*(scenario_overview_row["H2_demand_2030"]+YoY_2*y) for y in 1:20];[data["conv_factor"]*(scenario_overview_row["H2_demand_2030"]+YoY_2*20) for y in 1:data["nyears"]-31]]

    return H2
end