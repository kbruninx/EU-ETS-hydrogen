function update_supply!(e::Array,ETS::Dict,data::Dict,scenario_overview_row::DataFrameRow)
    # Impact of MSR is only enforced as of 2021. Before 2021, supply is fixed to historical values (see "define_ETS_parameters.jl")
    if scenario_overview_row[:MSR] == 2018 # The MSR according to the 2018 rules
        for y = 3:data["nyears"]
            for m = 1:12
                if m <= 8 # For the first 8 months, intake/outflow MSR depends on the TNAC in y-2
                    if ETS["TNAC"][y-2] >= data["TNAC_MAX"] # If exceeds TNAC_MAX - inflow
                        if min(ETS["CAP"][y],ETS["X_MSR_MAX_POS"][y]*ETS["TNAC"][y-2]) > ETS["X_MSR_MAX_POS"][y]*data["TNAC_MAX"] # only put EUAs in MSR if total exceeds 100/200 million
                                ETS["X_MSR"][y,m] = min(ETS["CAP"][y]/8,ETS["X_MSR_MAX_POS"][y]*ETS["TNAC"][y-2]/12) # one cannot put more in the MSR than what's still left to be auctioned/allocated
                        else
                                ETS["X_MSR"][y,m] = 0
                        end
                    elseif ETS["TNAC"][y-2] <= data["TNAC_MIN"] # if below TNAC_MIN - outflow 
                                ETS["X_MSR"][y,m] = -min(ETS["X_MSR_MAX_NEG"][y]/12,ETS["MSR"][y-1,12]/8) # outflow limited to what is in MSR
                    else
                                ETS["X_MSR"][y,m] = 0
                    end
                else # For the last 4 months, intake/outflow MSR depends on the TNAC in y-1
                    if ETS["TNAC"][y-1] >= data["TNAC_MAX"] # If exceeds TNAC_MAX - inflow
                        if min(ETS["CAP"][y],ETS["X_MSR_MAX_POS"][y]*ETS["TNAC"][y-1])> ETS["X_MSR_MAX_POS"][y]*data["TNAC_MAX"] # only put EUAs in MSR if total exceeds 200/100 million
                                ETS["X_MSR"][y,m] = min(ETS["CAP"][y]/4,ETS["X_MSR_MAX_POS"][y]*ETS["TNAC"][y-1]/12) # one cannot put more in the MSR than what's still left to be auctioned/allocated
                        else
                                ETS["X_MSR"][y,m] = 0
                        end
                    elseif ETS["TNAC"][y-1] <= data["TNAC_MIN"] # if below TNAC_MIN - outflow 
                                ETS["X_MSR"][y,m] = -min(ETS["X_MSR_MAX_NEG"][y]/12,ETS["MSR"][y,8]/4) # outflow limited to what is in MSR
                    else
                                ETS["X_MSR"][y,m] = 0
                    end
                end

                # Adapt MSR with backloaded/non-allocated allowances
                if m == 1
                        ETS["MSR"][y,m] = ETS["MSR"][y-1,12]+ETS["DELTA"][y]+ETS["X_MSR"][y,1]
                else
                        ETS["MSR"][y,m] = ETS["MSR"][y,m-1]+ETS["X_MSR"][y,m]
                end

                # Cancellation enforced as of 2023
                if ETS["MSR"][y,m] > 0.57*ETS["CAP"][y-1] && y >= 5  
                        ETS["C"][y,m] =  ETS["MSR"][y,m]-0.57*ETS["CAP"][y-1]
                        ETS["MSR"][y,m] = 0.57*ETS["CAP"][y-1]
                else
                        ETS["C"][y,m] = 0
                end
            end

            # Corrected supply of EUAS 
            ETS["S"][y]  = maximum([0,ETS["CAP"][y] - sum(ETS["X_MSR"][y,1:12]) + ETS["Δs"][y]]) 
          
            # TNAC
            ETS["TNAC"][y] = sum(ETS["CAP"][1:y]) + sum(ETS["Δs"][1:y]) + sum(ETS["DELTA"][1:y]) - sum(e[1:y]) - sum(ETS["C"][1:y,:]) - ETS["MSR"][y,12]
        end
    elseif scenario_overview_row[:MSR] == 2021 # The MSR as proposed in the Green Deal
        for y = 3:data["nyears"]
            for m = 1:12
                if m <= 8 # For the first 8 months, intake/outflow MSR depends on the TNAC in y-2
                    if ETS["TNAC"][y-2] >= data["TNAC_MAX"] # If exceeds TNAC_MAX - inflow
                        if ETS["TNAC"][y-2] <= data["TNAC_THRESHOLD"]  # between new threshold and maximum 
                            ETS["X_MSR"][y,m] = min(ETS["CAP"][y]/8,(ETS["TNAC"][y-2]-data["TNAC_MAX"])/12) # one cannot put more in the MSR than what's still left to be auctioned/allocated
                        else # regular intake rate 
                            ETS["X_MSR"][y,m] = min(ETS["CAP"][y]/8,ETS["X_MSR_MAX_POS"][y]*ETS["TNAC"][y-2]/12) # one cannot put more in the MSR than what's still left to be auctioned/allocated
                        end
                    elseif ETS["TNAC"][y-2] <= data["TNAC_MIN"] # if below TNAC_MIN - outflow 
                        ETS["X_MSR"][y,m] = -min(ETS["X_MSR_MAX_NEG"][y]/12,ETS["MSR"][y-1,12]/8) # outflow limited to what is in MSR
                    else
                        ETS["X_MSR"][y,m] = 0
                    end
                else # For the last 4 months, intake/outflow MSR depends on the TNAC in y-1
                    if ETS["TNAC"][y-1] >= data["TNAC_MAX"] # If exceeds TNAC_MAX - inflow
                        if ETS["TNAC"][y-1] <= data["TNAC_THRESHOLD"] # between new threshold and maximum 
                            ETS["X_MSR"][y,m] = min(ETS["CAP"][y]/4,(ETS["TNAC"][y-1]-data["TNAC_MAX"])/12) # one cannot put more in the MSR than what's still left to be auctioned/allocated
                        else # regular intake rate  
                            ETS["X_MSR"][y,m] = min(ETS["CAP"][y]/4,ETS["X_MSR_MAX_POS"][y]*ETS["TNAC"][y-1]/12) # one cannot put more in the MSR than what's still left to be auctioned/allocated
                        end
                    elseif ETS["TNAC"][y-1] <= data["TNAC_MIN"] # if below TNAC_MIN - outflow 
                        ETS["X_MSR"][y,m] = -min(ETS["X_MSR_MAX_NEG"][y]/12,ETS["MSR"][y,8]/4) # outflow limited to what is in MSR
                    else
                        ETS["X_MSR"][y,m] = 0
                    end
                end

                # Adapt MSR with backloaded/non-allocated allowances
                if m == 1
                        ETS["MSR"][y,m] = ETS["MSR"][y-1,12]+ETS["DELTA"][y]+ETS["X_MSR"][y,1]
                else
                        ETS["MSR"][y,m] = ETS["MSR"][y,m-1]+ETS["X_MSR"][y,m]
                end

                # Cancellation enforced as of 2023
                if ETS["MSR"][y,m] > data["TNAC_MIN"] && y >= 7  
                        ETS["C"][y,m] =  ETS["MSR"][y,m]-data["TNAC_MIN"]
                        ETS["MSR"][y,m] = data["TNAC_MIN"]
                else
                        ETS["C"][y,m] = 0
                end
            end

            # Corrected supply of EUAS 
            ETS["S"][y]  = maximum([0,ETS["CAP"][y] - sum(ETS["X_MSR"][y,1:12]) + ETS["Δs"][y]]) 
            
            # TNAC
            ETS["TNAC"][y] = sum(ETS["S"][1:y]) - sum(e[1:y]) 
        end
    end

    return ETS
end
