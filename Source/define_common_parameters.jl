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
    mod.ext[:parameters][:λ_REC] = zeros(data["nyears"],1)       # Price structure
    mod.ext[:parameters][:r_bar] = zeros(data["nyears"],1)       # ADMM penalty term
    mod.ext[:parameters][:ρ_REC] = data["rho_REC"]               # ADMM rho value 

    
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
    
    # Covered by ETS?
    if data["EOM"] == "YES" 
        mod.ext[:parameters][:EOM] = 1
        push!(agents[:eom],m)
    else
        mod.ext[:parameters][:EOM] = 0
    end

    return mod, agents
end