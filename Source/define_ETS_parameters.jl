function define_ETS_parameters!(ETS::Dict,data::Dict,scenario_overview_row::DataFrameRow)
    # Number of agents
    ETS["nAgents"] = data["nAgents"]

    # LRF 2017 - 2021, no Brexit, expansion of scope (aviation, maritime) or Green Deal
    ETS["LRF"] = zeros(data["nyears"]);
    ETS["LRF"][1:4] = data["LRF_2017"]*ones(4,1); 
    ETS["LRF"][5:14] = ETS["LRF"][1]*scenario_overview_row[:LRF_2021]/0.0174*ones(10,1);                            # 2021-2030
    ETS["LRF"][15:end] = ETS["LRF"][1]*scenario_overview_row[:LRF_2031]/0.0174*ones(data["nyears"]-14,1);           # 2030 - end ETS
       
    # CAP
    ETS["CAP"] = zeros(data["nyears"]);
    for y =1:data["nyears"]
        ETS["CAP"][y]= maximum([data["CAP_2016"]-sum(ETS["LRF"][1:y]) 0])
    end
    ETS["CAP"][1] = data["S_2017"]+data["EX_2016"] # real supply in 2017 and surplus in market at end 2016

    # Corrections to the cap and LRF if this not a calibration run
    if scenario_overview_row[:ref_scen_number] != scenario_overview_row[:scen_number] 
        if scenario_overview_row[:MSR] == 2018     # Acount for Brexit and the inclusion of aviation - rebasing the cap (see Data-file)
            ETS["LRF"][5:14] = ETS["LRF"][5]*(data["CAP_2021_NoGreenDeal"])/ETS["CAP"][5]*ones(10);                    # 2021-2030
            ETS["LRF"][15:end] = ETS["LRF"][15]*(data["CAP_2021_NoGreenDeal"])/ETS["CAP"][5]*ones(data["nyears"]-14);  # 2030 - end ETS
            ETS["CAP"][5] = data["CAP_2021_NoGreenDeal"]                                                               # 2021 Cap
        elseif scenario_overview_row[:MSR] == 2021  # Account for Green Deal, Brexit and the inclusion of aviation - rebasing the cap (see Data-file)
            ETS["LRF"][5:14] = ETS["LRF"][5]*(data["CAP_2021_GreenDeal"])/ETS["CAP"][5]*ones(10);                      # 2021-2030
            ETS["LRF"][15:end] = ETS["LRF"][15]*(data["CAP_2021_GreenDeal"])/ETS["CAP"][5]*ones(data["nyears"]-14);    # 2030 - end ETS
            ETS["CAP"][5] = data["CAP_2021_GreenDeal"]                                                                 # 2021 Cap
        end

        # Cap 2022-end ETS
        for y=6:data["nyears"]
            ETS["CAP"][y]= maximum([ETS["CAP"][5]-sum(ETS["LRF"][6:y]) 0])
        end
    end

    # Supply 
    ETS["S"] = copy(ETS["CAP"])

    # TNAC 
    ETS["TNAC"] = zeros(data["nyears"]);

    # MSR 
    ETS["X_MSR"] = zeros(data["nyears"],12);
    ETS["MSR"] = zeros(data["nyears"],12);
    ETS["C"] = zeros(data["nyears"],12);
    ETS["X_MSR_MAX_POS"] = zeros(data["nyears"])
    ETS["X_MSR_MAX_POS"][1:2] = zeros(2,1);
    if scenario_overview_row[:MSR] == 2018
        ETS["X_MSR_MAX_POS"][3:7] = data["X_MSR_MAX_POS_2019"]*ones(5);  # 24% until 2023
        ETS["X_MSR_MAX_POS"][8:end] = data["X_MSR_MAX_POS_2023"]*ones(data["nyears"]-7); 
    elseif scenario_overview_row[:MSR] == 2021
        ETS["X_MSR_MAX_POS"][3:14] = data["X_MSR_MAX_POS_2019"]*ones(12);  # 24% until 2030
        ETS["X_MSR_MAX_POS"][15:end] = data["X_MSR_MAX_POS_2023"]*ones(data["nyears"]-14); 
    end
    ETS["X_MSR_MAX_NEG"] = zeros(data["nyears"])
    ETS["X_MSR_MAX_NEG"][1:2] = zeros(2,1);
    ETS["X_MSR_MAX_NEG"][3:7] = data["X_MSR_MAX_NEG_2019"]*ones(5);  
    ETS["X_MSR_MAX_NEG"][8:end] = data["X_MSR_MAX_NEG_2023"]*ones(data["nyears"]-7); 
    ETS["TNAC_MAX"] = data["TNAC_MAX"]
    ETS["TNAC_MIN"] = data["TNAC_MIN"]
    ETS["TNAC_THRESHOLD"] = data["TNAC_THRESHOLD"]

    # Unallocated and backloaded allowances
    ETS["DELTA"] = zeros(data["nyears"]);
    ETS["DELTA"][3]= data["B_2016"]; # Backloaded EUAs are placed in MSR @ start 2019
    ETS["DELTA"][5] = data["UA_2020"]; # Unallocated EUAs from phase 4 @ start of 2021

    # Compute impact COVID-19 & overlapping policies
    ETS["Δe"] = zeros(data["nyears"]);
    ETS["Δs"] = zeros(data["nyears"]);
    for y = 1:data["nyears"]
        if scenario_overview_row[:COVID] == 1 
            if y == 4 # 2020
                ETS["Δe"][y] = - 240
            elseif y == 5 # 2021
                ETS["Δe"][y] = - 240*4/5
            elseif y == 6 # 2022
                ETS["Δe"][y] = - 240*3/5
            elseif y == 7 # 2023
                ETS["Δe"][y] = - 240*2/5
            elseif y == 8 # 2024
                ETS["Δe"][y] = - 240*1/5
            end
        end 

        # Overlapping policy on demand for emission allowances
        if y in range(scenario_overview_row[:start_op]-2016, stop=scenario_overview_row[:stop_op]-2016)
            ETS["Δe"][y] = ETS["Δe"][y] - scenario_overview_row[:op_dem]
        end

        if y in range(scenario_overview_row[:start_op]-2016, stop=scenario_overview_row[:stop_op]-2016)
            ETS["Δs"][y] = scenario_overview_row[:op_supply]
        end
    end

    return ETS
end