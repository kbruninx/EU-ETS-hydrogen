function define_ETS_parameters!(ETS::Dict,data::Dict)
    # LRF 
    LRF_tot = data["LRF_stat_2021"]/data["CAP_stat_2021"]*(data["CAP_stat_2021"]+data["CAP_aviation_2021"]+data["CAP_maritime_2021"])
    ETS["LRF"] = zeros(data["nyears"])                                                                                 
    ETS["LRF"][3:12] = LRF_tot*data["LRF_2021"]/0.022*ones(10,1)                            # 2021-2030
    ETS["LRF"][13:end] = LRF_tot*data["LRF_2031"]/0.022*ones(data["nyears"]-12,1)           # 2030 - end ETS
       
    # CAP
    ETS["CAP"] = zeros(data["nyears"])
    ETS["CAP"][1] = data["S_2019"]                                                                          # effective supply in 2019
    ETS["CAP"][2] = data["S_2020"] + data["TNAC_2020"]                                                      # effective supply in 2020 + TNAC at the end of 2020 
    ETS["CAP"][3] = data["CAP_stat_2021"]+data["CAP_aviation_2021"]+data["CAP_maritime_2021"]               # 2021
    for y =4:data["nyears"]                                                                                 # 2022-end
        ETS["CAP"][y]= maximum([ETS["CAP"][y-1]-ETS["LRF"][y-1] 0])
    end

    # Supply 
    ETS["S"] = copy(ETS["CAP"])

    # TNAC 
    ETS["TNAC"] = zeros(data["nyears"]);
    ETS["TNAC"][1] = data["TNAC_2019"]
    ETS["TNAC"][2] = data["TNAC_2020"]

    # MSR 
    ETS["X_MSR"] = zeros(data["nyears"],12);
    ETS["MSR"] = zeros(data["nyears"],12);
    ETS["MSR"][1] = data["MSR_2019"]
    ETS["MSR"][2] = data["MSR_2020"]
    ETS["C"] = zeros(data["nyears"],12);
    ETS["X_MSR_MAX_POS"] = zeros(data["nyears"])
    ETS["X_MSR_MAX_POS"][1:2] = zeros(2,1);
    if data["MSR"] == 2018
        ETS["X_MSR_MAX_POS"][1:5] = data["X_MSR_MAX_POS_2019"]*ones(5);  # 24% until 2023
        ETS["X_MSR_MAX_POS"][6:end] = data["X_MSR_MAX_POS_2023"]*ones(data["nyears"]-5); 
    elseif data["MSR"] == 2021
        ETS["X_MSR_MAX_POS"][1:12] = data["X_MSR_MAX_POS_2019"]*ones(12);  # 24% until 2030
        ETS["X_MSR_MAX_POS"][13:end] = data["X_MSR_MAX_POS_2023"]*ones(data["nyears"]-12); 
    end
    ETS["X_MSR_MAX_NEG"] = zeros(data["nyears"])
    ETS["X_MSR_MAX_NEG"][1:5] = data["X_MSR_MAX_NEG_2019"]*ones(5);  
    ETS["X_MSR_MAX_NEG"][6:end] = data["X_MSR_MAX_NEG_2023"]*ones(data["nyears"]-5); 
    ETS["TNAC_MAX"] = data["TNAC_MAX"]
    ETS["TNAC_MIN"] = data["TNAC_MIN"]
    ETS["TNAC_THRESHOLD"] = data["TNAC_THRESHOLD"]

    # Forced inflow to MSR (backloading, unallocated,...) - does not need to be used here
    ETS["DELTA"] = zeros(data["nyears"]);

    # Compute impact COVID-19 & overlapping policies - not used here
    ETS["Δe"] = zeros(data["nyears"]);
    ETS["Δs"] = zeros(data["nyears"]);

    # Historical prices
    ETS["P_2019"] = data["P_2019"]
    ETS["P_2020"] = data["P_2020"]
    ETS["P_2021"] = data["P_2021"]

    # Historical emissions, including voluntary cancellation
    ETS["E_2019"] = data["E_2019"]
    ETS["E_2020"] = data["E_2020"]
    ETS["E_2021"] = data["E_2021"]

    return ETS
end