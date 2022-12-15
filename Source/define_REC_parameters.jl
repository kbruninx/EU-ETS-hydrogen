function define_REC_parameters!(REC::Dict,data::Dict,ts::DataFrame,repr_days::DataFrame)
    # RES target
    REC["RT"] = [0; data["RES_target_2020"]*ones(10); data["RES_target_2030"]*ones(data["nyears"]-11)]
    
    # RES share "other technologies" in 2017
    REC["RS_other_2021"] = data["RES_share_other_2021"]

    return REC
end