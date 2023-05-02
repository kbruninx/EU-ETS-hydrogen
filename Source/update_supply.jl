function update_supply!(e::Array,ETS::Dict,data::Dict)
    # Impact of MSR is only enforced as of 2021. 
    # Note: 
    # TNAC will be shifted by 2 years (i.e., TNAC[y] is the TNAC at the end of year y-2)
    # MSR will be shifted by 1 yeare (i.e., MSR[y,12] is the MSR at the end of year y-1)

    if data["MSR"] == "YES" #   2021-2023 based on 2018 rules, 2024 - end based on 2023 rules
        for y = 1:data["nyears"]
            for m = 1:12
                if y <= 3 # MSR2018 
                    if m <= 8 # For the first 8 months, intake/outflow MSR depends on the TNAC in y-2
                        if ETS["TNAC"][y] >= data["TNAC_MAX"] # If exceeds TNAC_MAX - inflow
                            if min(ETS["CAP"][y],ETS["X_MSR_MAX_POS"][y]*ETS["TNAC"][y]) > ETS["X_MSR_MAX_POS"][y]*data["TNAC_MAX"] # only put EUAs in MSR if total exceeds 100/200 million
                                    ETS["X_MSR"][y,m] = min(ETS["CAP"][y]/8,ETS["X_MSR_MAX_POS"][y]*ETS["TNAC"][y]/12) # one cannot put more in the MSR than what's still left to be auctioned/allocated
                            else
                                    ETS["X_MSR"][y,m] = 0
                            end
                        elseif ETS["TNAC"][y] <= data["TNAC_MIN"] # if below TNAC_MIN - outflow 
                                    ETS["X_MSR"][y,m] = -min(ETS["X_MSR_MAX_NEG"][y]/12,ETS["MSR"][y,12]/8) # outflow limited to what is in MSR
                        else
                                    ETS["X_MSR"][y,m] = 0
                        end
                    else # For the last 4 months, intake/outflow MSR depends on the TNAC in y-1
                        if ETS["TNAC"][y+1] >= data["TNAC_MAX"] # If exceeds TNAC_MAX - inflow
                            if min(ETS["CAP"][y],ETS["X_MSR_MAX_POS"][y]*ETS["TNAC"][y+1])> ETS["X_MSR_MAX_POS"][y]*data["TNAC_MAX"] # only put EUAs in MSR if total exceeds 200/100 million
                                    ETS["X_MSR"][y,m] = min(ETS["CAP"][y]/4,ETS["X_MSR_MAX_POS"][y]*ETS["TNAC"][y+1]/12) # one cannot put more in the MSR than what's still left to be auctioned/allocated
                            else
                                    ETS["X_MSR"][y,m] = 0
                            end
                        elseif ETS["TNAC"][y+1] <= data["TNAC_MIN"] # if below TNAC_MIN - outflow 
                                    ETS["X_MSR"][y,m] = -min(ETS["X_MSR_MAX_NEG"][y]/12,ETS["MSR"][y+1,8]/4) # outflow limited to what is in MSR
                        else
                                    ETS["X_MSR"][y,m] = 0
                        end
                    end

                    # Adapt MSR with backloaded/non-allocated allowances
                    if m == 1
                            ETS["MSR"][y+1,m] = ETS["MSR"][y,12]+ETS["DELTA"][y]+ETS["X_MSR"][y,1]
                    else
                            ETS["MSR"][y+1,m] = ETS["MSR"][y+1,m-1]+ETS["X_MSR"][y,m]
                    end
                    
                    # Cancellation enforced as of 2023
                    if y >= 3
                        if ETS["MSR"][y+1,m] > 0.57*ETS["CAP"][y-1] 
                            ETS["C"][y,m] =  ETS["MSR"][y+1,m]-0.57*ETS["CAP"][y-1]
                            ETS["MSR"][y+1,m] = 0.57*ETS["CAP"][y-1]
                        end
                    else
                            ETS["C"][y,m] = 0
                    end
                else # MSR 2023
                    for m = 1:12
                        if m <= 8 # For the first 8 months, intake/outflow MSR depends on the TNAC in y-2
                            if ETS["TNAC"][y] >= data["TNAC_MAX"] # If exceeds TNAC_MAX - inflow
                                if ETS["TNAC"][y] <= data["TNAC_THRESHOLD"]  # between new threshold and maximum 
                                    ETS["X_MSR"][y,m] = min(ETS["CAP"][y]/8,(ETS["TNAC"][y]-data["TNAC_MAX"])/12) # one cannot put more in the MSR than what's still left to be auctioned/allocated
                                else # regular intake rate 
                                    ETS["X_MSR"][y,m] = min(ETS["CAP"][y]/8,ETS["X_MSR_MAX_POS"][y]*ETS["TNAC"][y]/12) # one cannot put more in the MSR than what's still left to be auctioned/allocated
                                end
                            elseif ETS["TNAC"][y] <= data["TNAC_MIN"] # if below TNAC_MIN - outflow 
                                ETS["X_MSR"][y,m] = -min(ETS["X_MSR_MAX_NEG"][y]/12,ETS["MSR"][y,12]/8) # outflow limited to what is in MSR
                            else
                                ETS["X_MSR"][y,m] = 0
                            end
                        else # For the last 4 months, intake/outflow MSR depends on the TNAC in y-1
                            if ETS["TNAC"][y+1] >= data["TNAC_MAX"] # If exceeds TNAC_MAX - inflow
                                if ETS["TNAC"][y+1] <= data["TNAC_THRESHOLD"] # between new threshold and maximum 
                                    ETS["X_MSR"][y,m] = min(ETS["CAP"][y]/4,(ETS["TNAC"][y+1]-data["TNAC_MAX"])/12) # one cannot put more in the MSR than what's still left to be auctioned/allocated
                                else # regular intake rate  
                                    ETS["X_MSR"][y,m] = min(ETS["CAP"][y]/4,ETS["X_MSR_MAX_POS"][y]*ETS["TNAC"][y+1]/12) # one cannot put more in the MSR than what's still left to be auctioned/allocated
                                end
                            elseif ETS["TNAC"][y+1] <= data["TNAC_MIN"] # if below TNAC_MIN - outflow 
                                ETS["X_MSR"][y,m] = -min(ETS["X_MSR_MAX_NEG"][y]/12,ETS["MSR"][y+1,8]/4) # outflow limited to what is in MSR
                            else
                                ETS["X_MSR"][y,m] = 0
                            end
                        end
        
                        # Adapt MSR with backloaded/non-allocated allowances
                        if m == 1
                                ETS["MSR"][y+1,m] = ETS["MSR"][y,12]+ETS["DELTA"][y]+ETS["X_MSR"][y,1]
                        else
                                ETS["MSR"][y+1,m] = ETS["MSR"][y+1,m-1]+ETS["X_MSR"][y,m]
                        end
        
                        # Cancellation enforced as of 2023
                        if ETS["MSR"][y+1,m] > data["TNAC_MIN"] && y >= 3  
                                ETS["C"][y,m] =  ETS["MSR"][y+1,m]-data["TNAC_MIN"]
                                ETS["MSR"][y+1,m] = data["TNAC_MIN"]
                        else
                                ETS["C"][y,m] = 0
                        end
                    end
                end
            end

            # Corrected supply of EUAS 
            ETS["S"][y]  = maximum([0,ETS["CAP"][y] - sum(ETS["X_MSR"][y,1:12]) + ETS["Δs"][y]]) 
            
            # TNAC
            ETS["TNAC"][y+2] =  ETS["TNAC"][y+1] + ETS["S"][y] - e[y] # recursive to account for "initial" state of TNAC (TNAC is added to supply in preceding year to make EUAs avaialble to actors, otherwise double counting)
        end      
    elseif data["MSR"] == "NO"  # No MSR, supply equals cap except in 2021
        ETS["S"] = copy(ETS["CAP"])
        for y = 1:data["nyears"]
            ETS["TNAC"][y+2] =  ETS["TNAC"][y+1] + ETS["S"][y] - e[y]
        end
    end

    # Correct supply in 2021 with TNAC at end of 2020 (needs to be made avaialble to market parties)
    ETS["S"][1] = ETS["S"][1]+ETS["TNAC"][2]

    return ETS
end
