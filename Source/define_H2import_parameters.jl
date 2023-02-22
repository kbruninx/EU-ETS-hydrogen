function define_H2import_parameters!(mod::Model, data::Dict,ts::DataFrame,repr_days::DataFrame,REC::Dict)
   
    mod.ext[:parameters][:α_H2_import] = data["α_H2_import"]
    
    if data["Additionality"] == "Yearly"
        mod.ext[:parameters][:ADD_SF] = ones(data["nyears"],1) 
    else
        mod.ext[:parameters][:ADD_SF] = REC["RT"] # Scaling factor for RECs when additionality is not applied = RES target
    end

    return mod
end