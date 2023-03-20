function define_ind_parameters!(mod::Model, data::Dict)
    # Emissions bound to historical values in 2019-2021
    mod.ext[:parameters][:e] = zeros(data["nyears"],1)  
    mod.ext[:parameters][:e][1] = data["E_ind_2021"]

    # β-value
    if data["scen_number"] - data["ref_scen_number"] == 0 && data["sens_number"] == 1 # this is a calibration run - provide an initial estimate
        mod.ext[:parameters][:β] = 1.2 #data["P_calibration"]/(data["E_ref"]*data["SF_ind"] - data["E_ind_2021"])^data["gamma"]
    else # get beta from reference result
        overview_results = CSV.read(joinpath(home_dir,string("overview_results_",data["nReprDays"],"_repr_days.csv")),DataFrame;delim=";")
        overview_results_row = filter(row -> row.scen_number in [data["ref_scen_number"]], overview_results)
        mod.ext[:parameters][:β] = overview_results_row[!,:Beta][1]
    end

    # Abatement costs
    mod.ext[:parameters][:AC] = zeros(data["nyears"],1)  
    
    return mod
end