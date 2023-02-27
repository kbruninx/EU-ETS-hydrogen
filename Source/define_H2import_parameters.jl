function define_H2import_parameters!(mod::Model, data::Dict,ts::DataFrame,repr_days::DataFrame,REC::Dict)
   
    if data["scen_number"] - data["ref_scen_number"] == 0 && data["sens_number"] == 1 # this is a calibration run - provide an initial estimate
        mod.ext[:parameters][:α_H2_import] = 11800
    else
        overview_results = CSV.read(joinpath(home_dir,string("overview_results_",data["nReprDays"],"_repr_days.csv")),DataFrame;delim=";")
        overview_results_row = filter(row -> row.scen_number in [data["ref_scen_number"]], overview_results)
        mod.ext[:parameters][:α_H2_import] = overview_results_row[!,:Alpha][1]
    end
    
    if data["Additionality"] == "Yearly"
        mod.ext[:parameters][:ADD_SF] = ones(data["nyears"],1) 
    else
        mod.ext[:parameters][:ADD_SF] = REC["RT"] # Scaling factor for RECs when additionality is not applied = RES target
    end

    return mod
end