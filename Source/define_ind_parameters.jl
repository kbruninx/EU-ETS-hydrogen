function define_ind_parameters!(mod::Model, data::Dict, scenario_overview_row::DataFrameRow)
    # Emissions bound to historical values in 2019-2021
    mod.ext[:parameters][:e] = zeros(data["nyears"],1)  
    mod.ext[:parameters][:e][1] = data["E_2019"]*data["SF_ind"] 
    mod.ext[:parameters][:e][2] = data["E_2020"]*data["SF_ind"]
    mod.ext[:parameters][:e][3] = data["E_2021"]*data["SF_ind"]

    # β-value
    if scenario_overview_row[:scen_number] - scenario_overview_row[:ref_scen_number] == 0 # this is a calibration run - provide an initial estimate
        mod.ext[:parameters][:β] = data["P_2021"]/(data["E_ref"]*data["SF_ind"] - data["E_2021"]*data["SF_ind"])^scenario_overview_row[:gamma]
    else # get beta from reference result
        overview_results = CSV.read(joinpath(home_dir,string("overview_results_",data["nReprDays"],"_repr_days.csv")),DataFrame;delim=";")
        overview_results_row = filter(row -> row.scen_number in [scenario_overview_row[:ref_scen_number]], overview_results)
        mod.ext[:parameters][:β] = overview_results_row[!,:Beta][1]
    end

    # Abatement costs
    mod.ext[:parameters][:AC] = zeros(data["nyears"],1)  
    
    return mod
end