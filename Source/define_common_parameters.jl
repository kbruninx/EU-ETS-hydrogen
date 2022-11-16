function define_common_parameters!(m::String,mod::Model, data::Dict, ts::DataFrame, repr_days::DataFrame, agents::Dict, scenario_overview_row::DataFrameRow)
    # Solver settings
    set_optimizer_attribute(mod, "OutputFlag",0)

    # Define dictonaries for sets, parameters, timeseries, variables, constraints & expressions
    mod.ext[:sets] = Dict()
    mod.ext[:parameters] = Dict()
    mod.ext[:timeseries] = Dict()
    mod.ext[:variables] = Dict()
    mod.ext[:constraints] = Dict()
    mod.ext[:expressions] = Dict()

    # Sets
    mod.ext[:sets][:JY] = 1:data["nyears"]    
    mod.ext[:sets][:JD] = 1:data["nReprDays"]
    mod.ext[:sets][:JH] = 1:data["nTimesteps"]

    # Parameters
    mod.ext[:parameters][:W] = Dict(jd => repr_days[!,:Period_weight][jd] for jd=1:data["nReprDays"])   # weights of each representative day
    mod.ext[:parameters][:A] = ones(data["nyears"],1)                                                   # Discount rate, 2019 as base year due to calibration to 2019 data
    for y in 4:data["nyears"]
        mod.ext[:parameters][:A][y] = 1/(1+data["discount_rate"])^(y-3)
    end

    # Parameters related to the EUA auctions
    mod.ext[:parameters][:λ_EUA] = zeros(data["nyears"],1)       # Price structure
    mod.ext[:parameters][:b_bar] = zeros(data["nyears"],1)       # ADMM penalty term
    mod.ext[:parameters][:ρ_EUA] = data["rho_EUA"]               # ADMM rho value 

    # Parameters related to the EOM
    mod.ext[:parameters][:λ_EOM] = zeros(data["nTimesteps"],data["nReprDays"],data["nyears"])   # Price structure
    mod.ext[:parameters][:g_bar] = zeros(data["nTimesteps"],data["nReprDays"],data["nyears"])   # ADMM penalty term
    mod.ext[:parameters][:ρ_EOM] = data["rho_EOM"]                                              # ADMM rho value 

    # Parameters related to the REC
    mod.ext[:parameters][:λ_y_REC] = zeros(data["nyears"],1)       # Price structure
    mod.ext[:parameters][:r_y_bar] = zeros(data["nyears"],1)       # ADMM penalty term
    if scenario_overview_row["Additionality"] == "Yearly" || scenario_overview_row["Additionality"] == "NA" 
        mod.ext[:parameters][:ρ_y_REC] = data["rho_REC"]            # ADMM rho value 
    else
        mod.ext[:parameters][:ρ_y_REC] = 0                          # ADMM rho value 
    end

    mod.ext[:parameters][:λ_d_REC] = zeros(data["nReprDays"],data["nyears"])       # Price structure
    mod.ext[:parameters][:r_d_bar] = zeros(data["nReprDays"],data["nyears"])       # ADMM penalty term
    if scenario_overview_row["Additionality"] == "Daily" 
        mod.ext[:parameters][:ρ_d_REC] = data["rho_REC"]                            # ADMM rho value 
    else
        mod.ext[:parameters][:ρ_d_REC] = 0                                          # ADMM rho value 
    end
    
    mod.ext[:parameters][:λ_h_REC] = zeros(data["nTimesteps"],data["nReprDays"],data["nyears"])         # Price structure
    mod.ext[:parameters][:r_h_bar] = zeros(data["nTimesteps"],data["nReprDays"],data["nyears"])         # ADMM penalty term
    if scenario_overview_row["Additionality"] == "Hourly" 
        mod.ext[:parameters][:ρ_h_REC] = data["rho_REC"]                                                # ADMM rho value 
    else
        mod.ext[:parameters][:ρ_h_REC] = 0                                                              # ADMM rho value 
    end
    
    # Parameters related to the H2 market
    mod.ext[:parameters][:λ_H2] = zeros(data["nyears"],1)       # Price structure
    mod.ext[:parameters][:gH_bar] = zeros(data["nyears"],1)     # ADMM penalty term
    mod.ext[:parameters][:ρ_H2] = data["rho_H2"]                # ADMM rho value 

    # Parameters related to the carbon-neutral H2 generation subsidy
    mod.ext[:parameters][:λ_H2CN_prod] = zeros(data["nyears"],1)        # Price structure
    mod.ext[:parameters][:gHCN_bar] = zeros(data["nyears"],1)           # ADMM penalty term
    mod.ext[:parameters][:ρ_H2CN_prod] = data["rho_H2CN_prod"]          # ADMM rho value

    # Parameters related to the carbon-neutral H2 production capacity subsidy
    mod.ext[:parameters][:λ_H2CN_cap] = zeros(data["nyears"],1)         # Price structure
    mod.ext[:parameters][:capHCN_bar] = zeros(data["nyears"],1)         # ADMM penalty term
    mod.ext[:parameters][:ρ_H2CN_cap] = data["rho_H2CN_cap"]            # ADMM rho value 

    # Parameters related to the natural gas market
    mod.ext[:parameters][:λ_NG] = zeros(data["nyears"],1)               # Price structure

    # Eligble for RECs?
    if data["REC"] == "YES" 
        mod.ext[:parameters][:REC] = 1
        push!(agents[:rec],m)
    else
        mod.ext[:parameters][:REC] = 0
    end

    # Covered by ETS?
    if data["ETS"] == "YES" 
        mod.ext[:parameters][:ETS] = 1
        push!(agents[:ets],m)
    else
        mod.ext[:parameters][:ETS] = 0
    end
    
    # Covered by EOM?
    if data["EOM"] == "YES" 
        mod.ext[:parameters][:EOM] = 1
        push!(agents[:eom],m)
    else
        mod.ext[:parameters][:EOM] = 0
    end

    # Covered by Hydrogen Market 
    if data["H2"] == "YES" 
        mod.ext[:parameters][:H2] = 1
        push!(agents[:h2],m)
    else
        mod.ext[:parameters][:H2] = 0
    end

    # Covered by incentive scheme for carbon neutral hydrogen?
    if data["H2CN_prod"] == "YES" 
        mod.ext[:parameters][:H2CN_prod] = 1
        push!(agents[:h2cn_prod],m)
    else
        mod.ext[:parameters][:H2CN_prod] = 0
    end

     # Covered by incentive scheme for carbon neutral hydrogen?
     if data["H2CN_cap"] == "YES" 
        mod.ext[:parameters][:H2CN_cap] = 1
        push!(agents[:h2cn_cap],m)
    else
        mod.ext[:parameters][:H2CN_cap] = 0
    end

     # Covered by natural gas market
     if data["NG"] == "YES" 
        mod.ext[:parameters][:NG] = 1
        push!(agents[:ng],m)
    else
        mod.ext[:parameters][:NG] = 0
    end

    return mod, agents
end