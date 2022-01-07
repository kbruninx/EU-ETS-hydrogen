## Topic: EU ETS, MSR and overlapping policies
# Author: Kenneth Bruninx
# Last update: November 2021

## 0. Set-up code
# Range of scenarios to be simulated
start_scen = 1
stop_scen = 2

# Include packages 
using JuMP, Gurobi # Optimization packages
using DataFrames, CSV, YAML # dataprocessing
using ProgressBars, Printf # progress bar

# Gurobi environment to suppress output
const GUROBI_ENV = Gurobi.Env()

# Include functions
include(joinpath(@__DIR__,"BUILD_ind_agent.jl"))
include(joinpath(@__DIR__,"SOLVE_ind_agent.jl"))
include(joinpath(@__DIR__,"BUILD_ps_agent.jl"))
include(joinpath(@__DIR__,"SOLVE_ps_agent.jl"))

# Data common to all scenarios data 
data = YAML.load_file(joinpath(@__DIR__,"Input","overview_data.yaml"))
ts = CSV.read(joinpath(@__DIR__,"Input","timeseries_May2018.csv"),delim=",",DataFrame)
repr_days = CSV.read(joinpath(@__DIR__,"Input",string("Period_definitions_",data["nReprDays"],".csv")),delim=",",DataFrame)

# Overview scenarios
scenario_overview = CSV.read(joinpath(@__DIR__,"overview_scenarios.csv"),DataFrame;delim=";")

# Create file with results 
if isfile(joinpath(@__DIR__,"overview_results.csv")) != 1
    CSV.write(joinpath(@__DIR__,"overview_results.csv"),DataFrame(),delim=";",header=["scen_number";"n_iter";"walltime";"PrimalResidual";"DualResidual";"Beta";"EUA_2021";"CumulativeEmissions";"Cancellation"; "WBseal"; "WBL"; "dirWBL"; "indirWBL"])
end

# Create folder for results
if isdir(joinpath(@__DIR__,"Results")) != 1
    mkdir(joinpath(@__DIR__,"Results"))
end

# Scenario number 
# for scen_number in range(start_scen,stop=stop_scen,step=1)
scen_number = 2

println("    ")
println(string("######################                  Scenario ",scen_number,"                 #########################"))
println("    ")
println("Including all required input data: done")
println("   ")

## 1. Read associated input for this simulation
scenario_overview_row = scenario_overview[scen_number,:]

## 2. Initiate model, sets and parameters for representative agents & ETS
mod_ind = Model(() -> Gurobi.Optimizer(GUROBI_ENV))
set_optimizer_attribute(mod_ind, "OutputFlag",0)
mod_ps = Model(() -> Gurobi.Optimizer(GUROBI_ENV))
set_optimizer_attribute(mod_ps, "OutputFlag",0)

# Parameters for representative industry agent
function define_ind_parameters!(mod_ind::Model, data::Dict, scenario_overview_row::DataFrameRow)
    # define sets
    mod_ind.ext[:sets] = Dict()
    mod_ind.ext[:sets][:JY] = 1:data["nyears"]
    
    # define parameters 
    mod_ind.ext[:parameters] = Dict()

    # Emissions representative agents, bound to historical values in 2017-2019
    mod_ind.ext[:parameters][:e] = zeros(data["nyears"],1)
    mod_ind.ext[:parameters][:e][1] = data["E_2017"]
    mod_ind.ext[:parameters][:e][2] = data["E_2018"]
    mod_ind.ext[:parameters][:e][3] = data["E_2019"]

    # β-value
    if scenario_overview_row[:scen_number] - scenario_overview_row[:ref_scen_number] == 0 # this is a calibration run - provide an initial estimate
        mod_ind.ext[:parameters][:β] = data["P_2019"]/(data["MACC"]["E_ref"] - data["E_2019"])^scenario_overview_row[:gamma]
    else # get beta from reference result
        overview_results = CSV.read(joinpath(@__DIR__,"overview_results.csv"),DataFrame;delim=";")
        overview_results_row = filter(row -> row.scen_number in [scenario_overview_row[:ref_scen_number]], overview_results)
        mod_ind.ext[:parameters][:β] = overview_results_row[!,:Beta][1]
    end

    # Parameters related to the EUA auctions -- note these appear in both the power sector and the industry agent
    mod_ind.ext[:parameters][:A] = ones(data["nyears"],1)       # Discount rate, 2019 as base year due to calibration to 2019 data
    for y in 4:data["nyears"]
        mod_ind.ext[:parameters][:A][y] = 1/(1+data["discount_rate"])^(y-3);
    end
    mod_ind.ext[:parameters][:λ] = zeros(data["nyears"],1)      # Price structure
    mod_ind.ext[:parameters][:b_bar] = zeros(data["nyears"],1)  # ADMM penalty term
    mod_ind.ext[:parameters][:ρ] = data["ADMM"]["rho"]          # ADMM rho value 
    
    return mod_ind
end
define_ind_parameters!(mod_ind,data,scenario_overview_row)

# Parameters for power sector agent
function define_ps_parameters!(mod_ps::Model, data::Dict,ts::DataFrame,repr_days::DataFrame,scenario_overview_row::DataFrameRow)
    # define sets
    mod_ps.ext[:sets] = Dict()
    JY = mod_ps.ext[:sets][:JY] = 1:data["nyears"]
    JD = mod_ps.ext[:sets][:JD] = 1:data["nReprDays"]
    JH = mod_ps.ext[:sets][:JH] = 1:data["nTimesteps"]
    ID = mod_ps.ext[:sets][:ID] = [id for id in keys(data["dispatchableGenerators"])] # dispatchable generators
    IV = mod_ps.ext[:sets][:IV] = [iv for iv in keys(data["intermittentGenerators"])] # intermittent generators
    I = mod_ps.ext[:sets][:I] = union(mod_ps.ext[:sets][:ID], mod_ps.ext[:sets][:IV]) # all generators

    # define parameters 
    mod_ps.ext[:parameters] = Dict()
    mod_ps.ext[:parameters][:VOLL] = data["VOLL"] # VOLL
    mod_ps.ext[:parameters][:W] = Dict(jd => repr_days[!,:Period_weight][jd] for jd in JD) # weights of each representative day

    # Time dependency of fuel costs? investment costs? not accounted for. parameters in Yaml or via scenarios?

    d = merge(data["dispatchableGenerators"],data["intermittentGenerators"])
    mod_ps.ext[:parameters][:VC] = Dict(i => d[i]["fuelCosts"] for i in ID) # EUR/MWh
    mod_ps.ext[:parameters][:CI] = Dict(id => d[id]["emissions"] for id in ID) # tCO2/MWh
    mod_ps.ext[:parameters][:IC] = Dict(i => d[i]["OC"] for i in I) # EUR/MW
    mod_ps.ext[:parameters][:DELTA_CAP_MAX] = Dict(i => d[i]["max_YoY_new_cap"] for i in I) # MW

    # Lead time on new capacity, salvage value and availability of legacy capacity
    LeadT = mod_ps.ext[:parameters][:LeadTime] = Dict(i => d[i]["Leadtime"] for i in I) # years
    LifeT = mod_ps.ext[:parameters][:Lifetime] = Dict(i => d[i]["Lifetime"] for i in I) # years
    LegCap2017 = mod_ps.ext[:parameters][:Legcap2017] = Dict(i => d[i]["Legcap"] for i in I) # MW
    LegCap0 = mod_ps.ext[:parameters][:Legcap0] = Dict(i => d[i]["Legcap_out"] for i in I) # years
    mod_ps.ext[:parameters][:CAP_SV] = Dict(i => [maximum([0,1-(data["nyears"]-jy+1)/LifeT[i]]) for jy in JY] for i in I) 
    mod_ps.ext[:parameters][:LEG_CAP] = Dict(i => [d[i]["AF"]*LegCap2017[i]*maximum([0,(LegCap0[i]-jy+1)/LegCap0[i]]) for jy in JY] for i in I) 
    mod_ps.ext[:parameters][:CAP_LT] = Dict(i => zeros(data["nyears"],data["nyears"]) for i in I) 
    for i in I        
        for y in JY
            if y+LeadT[i] < data["nyears"]
                for yy = y+LeadT[i]:minimum([y+LeadT[i]+LifeT[i],data["nyears"]])
                    mod_ps.ext[:parameters][:CAP_LT][i][y,yy]  = 1
                end
            end
        end
    end

    # Growth factor of demand based on 2016 reference scenario
    GF_YoY =  mod_ps.ext[:parameters][:GF_YoY] = [0;  0.001*ones(2,1); 0.0045*ones(11,1); 0.0071*ones(data["nyears"]-14,1)]
    mod_ps.ext[:parameters][:GF] = [sum(GF_YoY[1:jy]) for jy in JY]

    # timeseries
    mod_ps.ext[:timeseries] = Dict()
    mod_ps.ext[:timeseries][:D_tot] = [data["SF_load"]*ts[!,:LOAD][round(Int,data["nTimesteps"]*repr_days[!,:Period_index][jd]+jh)]/1000 for jh in JH, jd in JD] # from MWh to GWh
    mod_ps.ext[:timeseries][:G_other] = [ts[!,:OTHER][round(Int,data["nTimesteps"]*repr_days[!,:Period_index][jd]+jh)]/1000 for jh in JH, jd in JD] # from MWh to GWh
    mod_ps.ext[:timeseries][:D] = mod_ps.ext[:timeseries][:D_tot] - mod_ps.ext[:timeseries][:G_other] 
    mod_ps.ext[:timeseries][:AF] = Dict(iv => zeros(data["nTimesteps"],data["nReprDays"]) for iv in IV) 
    mod_ps.ext[:timeseries][:AF]["Solar"] =  [ts[!,:SOLAR][round(Int,data["nTimesteps"]*repr_days[!,:Period_index][jd]+jh)]/1000/LegCap2017["Solar"] for jh in JH, jd in JD] # from MWh to GWh
    mod_ps.ext[:timeseries][:AF]["WindOffshore"] =  [ts[!,:WIND_OFFSHORE][round(Int,data["nTimesteps"]*repr_days[!,:Period_index][jd]+jh)]/1000/LegCap2017["WindOffshore"] for jh in JH, jd in JD] # from MWh to GWh
    mod_ps.ext[:timeseries][:AF]["WindOnshore"] =  [ts[!,:WIND_ONSHORE][round(Int,data["nTimesteps"]*repr_days[!,:Period_index][jd]+jh)]/1000/LegCap2017["WindOnshore"] for jh in JH, jd in JD] # from MWh to GWh

    # RES target: assumes that "other" RES (hydro, biomass,..) are replaced with "other" RES after they reach end of lifetime
    D_tot_cum = [(1+mod_ps.ext[:parameters][:GF][jy])*sum(mod_ps.ext[:parameters][:W][jd]*mod_ps.ext[:timeseries][:D_tot][jh,jd] for jh in JH, jd in JD) for jy in JY]
    RES_output_2017 = sum(mod_ps.ext[:parameters][:W][jd]*mod_ps.ext[:timeseries][:AF][iv][jh,jd]*LegCap2017[iv] for iv in IV, jh in JH, jd in JD)
    RT_rel = [zeros(3,1); RES_output_2017/D_tot_cum[4]*ones(10,1); data["RES_target_2030"]*ones(data["nyears"]-13,1)]
    mod_ps.ext[:parameters][:RT] = [RT_rel[jy]*D_tot_cum[jy] - data["RES_share_other_2017"]*D_tot_cum[1] for jy in JY]

    # Parameters related to the EUA auctions -- note these appear in both the power sector and the industry agent
    mod_ps.ext[:parameters][:A] = ones(data["nyears"],1)       # Discount rate, 2019 as base year due to calibration to 2019 data
    for y in 4:data["nyears"]
        mod_ps.ext[:parameters][:A][y] = 1/(1+data["discount_rate"])^(y-3);
    end
    mod_ps.ext[:parameters][:λ_EUA] = zeros(data["nyears"],1)      # Price structure
    mod_ps.ext[:parameters][:b_bar] = zeros(data["nyears"],1)  # ADMM penalty term
    mod_ps.ext[:parameters][:ρ_EUA] = data["ADMM"]["rho"]          # ADMM rho value 

    return mod_ps
end
define_ps_parameters!(mod_ps,data,ts,repr_days,scenario_overview_row)

# Parameters/variables ETS 
ETS = Dict()
function define_ETS_parameters!(ETS::Dict,data::Dict,scenario_overview_row::DataFrameRow)
    # LRF 2017 - 2021, no Brexit, expansion of scope (aviation, maritime) or Green Deal
    ETS["LRF"] = zeros(data["nyears"],1);
    ETS["LRF"][1:4] = data["LRF_2017"]*ones(4,1); 
    ETS["LRF"][5:14] = ETS["LRF"][1]*scenario_overview_row[:LRF_2021]/0.0174*ones(10,1);                            # 2021-2030
    ETS["LRF"][15:end] = ETS["LRF"][1]*scenario_overview_row[:LRF_2031]/0.0174*ones(data["nyears"]-14,1);           # 2030 - end ETS
       
    # CAP
    ETS["CAP"] = zeros(data["nyears"],1);
    for y =1:data["nyears"]
        ETS["CAP"][y]= maximum([data["CAP_2016"]-sum(ETS["LRF"][1:y]) 0])
    end
    ETS["CAP"][1] = data["S_2017"]+data["EX_2016"] # real supply in 2017 and surplus in market at end 2016

    # Corrections to the cap and LRF if this not a calibration run
    if scenario_overview_row[:ref_scen_number] != scenario_overview_row[:scen_number] 
        if scenario_overview_row[:MSR] == 2018     # Acount for Brexit and the inclusion of aviation - rebasing the cap (see Data-file)
            ETS["LRF"][5:14] = ETS["LRF"][5]*(data["CAP_2021_NoGreenDeal"])/ETS["CAP"][5]*ones(10,1);                    # 2021-2030
            ETS["LRF"][15:end] = ETS["LRF"][15]*(data["CAP_2021_NoGreenDeal"])/ETS["CAP"][5]*ones(data["nyears"]-14,1);  # 2030 - end ETS
            ETS["CAP"][5] = data["CAP_2021_NoGreenDeal"]                                                                 # 2021 Cap
        elseif scenario_overview_row[:MSR] == 2021  # Account for Green Deal, Brexit and the inclusion of aviation - rebasing the cap (see Data-file)
            ETS["LRF"][5:14] = ETS["LRF"][5]*(data["CAP_2021_GreenDeal"])/ETS["CAP"][5]*ones(10,1);                      # 2021-2030
            ETS["LRF"][15:end] = ETS["LRF"][15]*(data["CAP_2021_GreenDeal"])/ETS["CAP"][5]*ones(data["nyears"]-14,1);    # 2030 - end ETS
            ETS["CAP"][5] = data["CAP_2021_GreenDeal"]                                                                   # 2021 Cap
        end

        # Cap 2022-end ETS
        for y=6:data["nyears"]
            ETS["CAP"][y]= maximum([ETS["CAP"][5]-sum(ETS["LRF"][6:y]) 0])
        end
    end

    # Supply 
    ETS["S"] = copy(ETS["CAP"])

    # TNAC 
    ETS["TNAC"] = zeros(data["nyears"],1);

    # MSR 
    ETS["X_MSR"] = zeros(data["nyears"],12);
    ETS["MSR"] = zeros(data["nyears"],12);
    ETS["C"] = zeros(data["nyears"],12);
    ETS["X_MSR_MAX_POS"] = zeros(data["nyears"],1)
    ETS["X_MSR_MAX_POS"][1:2] = zeros(2,1);
    if scenario_overview_row[:MSR] == 2018
        ETS["X_MSR_MAX_POS"][3:7] = data["X_MSR_MAX_POS_2019"]*ones(5,1);  # 24% until 2023
        ETS["X_MSR_MAX_POS"][8:end] = data["X_MSR_MAX_POS_2023"]*ones(data["nyears"]-7); 
    elseif scenario_overview_row[:MSR] == 2021
        ETS["X_MSR_MAX_POS"][3:14] = data["X_MSR_MAX_POS_2019"]*ones(12,1);  # 24% until 2030
        ETS["X_MSR_MAX_POS"][15:end] = data["X_MSR_MAX_POS_2023"]*ones(data["nyears"]-14); 
    end
    ETS["X_MSR_MAX_NEG"] = zeros(data["nyears"],1)
    ETS["X_MSR_MAX_NEG"][1:2] = zeros(2,1);
    ETS["X_MSR_MAX_NEG"][3:7] = data["X_MSR_MAX_NEG_2019"]*ones(5,1);  
    ETS["X_MSR_MAX_NEG"][8:end] = data["X_MSR_MAX_NEG_2023"]*ones(data["nyears"]-7); 

    # Unallocated and backloaded allowances
    ETS["DELTA"] = zeros(data["nyears"],1);
    ETS["DELTA"][3]= data["B_2016"]; # Backloaded EUAs are placed in MSR @ start 2019
    ETS["DELTA"][5] = data["UA_2020"]; # Unallocated EUAs from phase 4 @ start of 2021

    return ETS
end
define_ETS_parameters!(ETS,data,scenario_overview[scen_number,:])

println("Inititate model, sets and parameters: done")
println("   ")

## 3. Build models
BUILD_ind_agent!(mod_ind)
BUILD_ps_agent!(mod_ps)

println("Build model: done")
println("   ")

## 4. ADMM proces
# Auxiliary function: update emissions representative industry agent @ current emission allowance prices
function update_emissions!(mod_ind::Model,data::Dict,scenario_overview_row::DataFrameRow)
    # Baseline emissions 
    E_REF = data["MACC"]["E_ref"]*ones(data["nyears"],1)
    # Correction depending on the case at hand
    # if scenario_overview_row[:ref_scen_number] != scenario_overview_row[:scen_number] 
    #     if scenario_overview_row[:MSR] == 2018     # Acount for Brexit and the inclusion of aviation - rebasing the baseline emissions based on 2019 emission data
    #         E_REF[5:end] = (data["MACC"]["E_ref"] - 115 + 68)*ones(data["nyears"]-4,1)      
    #     elseif scenario_overview_row[:MSR] == 2021  # Account for Green Deal, Brexit and the inclusion of aviation - rebasing the cap (see Data-file)
    #         E_REF[5:end] = (data["MACC"]["E_ref"] - 115 + 68  + 90)*ones(data["nyears"]-4,1)      
    #     end
    # end

    for y = 4:data["nyears"] # Emissions are only price-dependent as of 2020 (y=4)
        λ_nom = maximum([0,mod_ind.ext[:parameters][:λ][y]/(1+data["inflation"])^(y-3)]) # M€/MtCO2, discounted to 2019 values, limited to positive values

        if y in range(scenario_overview_row[:start_op]-2016, stop=scenario_overview_row[:stop_op]-2016)
            mod_ind.ext[:parameters][:e][y] = maximum([0,E_REF[y] - (λ_nom/mod_ind.ext[:parameters][:β])^(1/scenario_overview_row[:gamma])])  - scenario_overview_row[:op_dem]
        else
            mod_ind.ext[:parameters][:e][y] = maximum([0,E_REF[y] - (λ_nom/mod_ind.ext[:parameters][:β])^(1/scenario_overview_row[:gamma])])
        end # to Mton

        # Account for maximum price-induced change in emission allowance demand (as percentage of 2019 emissions):
        if mod_ind.ext[:parameters][:e][y-1] > 1
            if mod_ind.ext[:parameters][:e][y-1] - mod_ind.ext[:parameters][:e][y] > scenario_overview_row[:max_em_change]*mod_ind.ext[:parameters][:e][3]/100
                mod_ind.ext[:parameters][:e][y] = mod_ind.ext[:parameters][:e][y-1]-scenario_overview_row[:max_em_change]*mod_ind.ext[:parameters][:e][3]/100
            end
        end

        # Correct for impact COVID-19
        if scenario_overview_row[:COVID] == 1
            if y == 4 # 2020
                mod_ind.ext[:parameters][:e][y] = mod_ind.ext[:parameters][:e][y] - 240
            elseif y == 5 # 2021
                mod_ind.ext[:parameters][:e][y] = mod_ind.ext[:parameters][:e][y] - 240*4/5
            elseif y == 6 # 2022
                mod_ind.ext[:parameters][:e][y] = mod_ind.ext[:parameters][:e][y] - 240*3/5
            elseif y == 7 # 2023
                mod_ind.ext[:parameters][:e][y] = mod_ind.ext[:parameters][:e][y] - 240*2/5
            elseif y == 8 # 2024
                mod_ind.ext[:parameters][:e][y] = mod_ind.ext[:parameters][:e][y] - 240*1/5
            end
        end
    
        if mod_ind.ext[:parameters][:e][y] < 0
            mod_ind.ext[:parameters][:e][y] = 0
        end

        # Scale down to reflect industry share: 
        mod_ind.ext[:parameters][:e][y] = mod_ind.ext[:parameters][:e][y]*data["SF_Ind"]
    end

    

    return mod_ind
end

# Auxiliary function: update supply of allowances via MSR
function update_supply!(e::Array,ETS::Dict,data::Dict,scenario_overview_row::DataFrameRow)
    if scenario_overview_row[:MSR] == 2018 # The MSR according to the 2018 rules
        for y = 1:data["nyears"]
            if y >= 3 # MSR only active after 2019
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
                if ETS["MSR"][y,m] > 0.57*ETS["CAP"][y-1] && y >= 7  
                        ETS["C"][y,m] =  ETS["MSR"][y,m]-0.57*ETS["CAP"][y-1]
                        ETS["MSR"][y,m] = 0.57*ETS["CAP"][y-1]
                else
                        ETS["C"][y,m] = 0
                end
            end
            end

            # Corrected supply of EUAS 
            if y in range(scenario_overview_row[:start_op]-2016, stop=scenario_overview_row[:stop_op]-2016) # overlapping policy in play
                ETS["S"][y]  = maximum([0,ETS["CAP"][y] - sum(ETS["X_MSR"][y,1:12]) - scenario_overview_row[:op_supply]]) 
            else
                ETS["S"][y]  = ETS["CAP"][y] - sum(ETS["X_MSR"][y,1:12])
            end

            # TNAC
            if y >= scenario_overview_row[:start_op]-2016 # overlapping policy in play
                ETS["TNAC"][y] = sum(ETS["CAP"][1:y]) - scenario_overview_row[:op_supply]*(minimum([y,scenario_overview_row[:stop_op]-2016])-(scenario_overview_row[:start_op]-2016)+1) + sum(ETS["DELTA"][1:y]) - sum(e[1:y]) - sum(ETS["C"][1:y,:]) - ETS["MSR"][y,12]
            else
                ETS["TNAC"][y] = sum(ETS["CAP"][1:y]) + sum(ETS["DELTA"][1:y]) - sum(e[1:y]) - sum(ETS["C"][1:y,:]) - ETS["MSR"][y,12]
            end
        end
    elseif scenario_overview_row[:MSR] == 2021 # The MSR as proposed in the Green Deal
        for y = 1:data["nyears"]
            if y >= 3 # MSR only active after 2019
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
            end

            # Corrected supply of EUAS 
            if y in range(scenario_overview_row[:start_op]-2016, stop=scenario_overview_row[:stop_op]-2016) # overlapping policy in play
                ETS["S"][y]  = maximum([0,ETS["CAP"][y] - sum(ETS["X_MSR"][y,1:12]) - scenario_overview_row[:op_supply]]) 
            else
                ETS["S"][y]  = ETS["CAP"][y] - sum(ETS["X_MSR"][y,1:12])
            end
            
            # TNAC
            ETS["TNAC"][y] = sum(ETS["S"][1:y]) - sum(e[1:y]) 
        end
    end

    return ETS
end

# Auxiliary function: initialize results -- only these results are retained in every iteration
results = Dict()
ADMM = Dict()
function define_results!(data::Dict,results::Dict,ADMM::Dict)
    results["Ind"] = Dict()
    results["Ind"]["b"] = zeros(data["ADMM"]["max_iter"],data["nyears"])
    results["Ind"]["e"] = zeros(data["ADMM"]["max_iter"],data["nyears"])
    results["PS"] = Dict()
    results["PS"]["b"] = zeros(data["ADMM"]["max_iter"],data["nyears"])
    results["PS"]["e"] = zeros(data["ADMM"]["max_iter"],data["nyears"])
    results["s"] = zeros(data["ADMM"]["max_iter"],data["nyears"])
    results["λ"] = zeros(data["ADMM"]["max_iter"],data["nyears"])

    ADMM["Imbalances"] = Dict()
    ADMM["Imbalances"]["ETS"] = zeros(data["ADMM"]["max_iter"],data["nyears"])
    ADMM["Imbalances"]["MSR"] = zeros(data["ADMM"]["max_iter"],data["nyears"])
    ADMM["Residuals"] = Dict()
    ADMM["Residuals"]["Primal"] = zeros(data["ADMM"]["max_iter"],1)
    ADMM["Residuals"]["Dual"] = zeros(data["ADMM"]["max_iter"],1)
    ADMM["Tolerance"] = data["ADMM"]["epsilon"]*sqrt(2*data["nyears"])  
    ADMM["ρ"] = data["ADMM"]["rho"]*ones(data["ADMM"]["max_iter"],1)
    ADMM["n_iter"] = 1 
    ADMM["walltime"] = 0

    return results, ADMM
end

# ADMM 
function ADMM!(results::Dict,ADMM::Dict,ETS::Dict,mod_ind::Model,scenario_overview_row::DataFrameRow)
    convergence = 0
    iterations = ProgressBar(1:data["ADMM"]["max_iter"]-1)
    for iter in iterations
        if convergence == 0
            # Calculate penalty terms ADMM and update price to most recent value (both agents)
            if iter > 1
                mod_ind.ext[:parameters][:b_bar] = [results["Ind"]["b"][iter-1,jy] + 1/2*ADMM["Imbalances"]["ETS"][iter-1,jy] for jy in mod_ind.ext[:sets][:JY]]
                mod_ps.ext[:parameters][:b_bar] = [results["PS"]["b"][iter-1,jy] + 1/2*ADMM["Imbalances"]["ETS"][iter-1,jy] for jy in mod_ps.ext[:sets][:JY]]
            end
            mod_ps.ext[:parameters][:λ_EUA] = mod_ind.ext[:parameters][:λ] = results["λ"][iter,:] 
            mod_ps.ext[:parameters][:ρ_EUA] = mod_ind.ext[:parameters][:ρ] = ADMM["ρ"][iter]

            # Industry agent  
            update_emissions!(mod_ind,data,scenario_overview_row) # Update emissions 
            results["Ind"]["e"][iter,:] = [mod_ind.ext[:parameters][:e][jy] for jy in mod_ind.ext[:sets][:JY]]
            SOLVE_ind_agent!(mod_ind) # Solve updated decision problem 
            results["Ind"]["b"][iter,:] = [value.(mod_ind.ext[:variables][:b])[jy] for jy in mod_ind.ext[:sets][:JY]]

            # Power sector agent
            SOLVE_ps_agent!(mod_ps) # Solve updated decision problem 
            results["PS"]["b"][iter,:] = [value.(mod_ps.ext[:variables][:b])[jy] for jy in mod_ps.ext[:sets][:JY]]
            results["PS"]["e"][iter,:] = [value.(mod_ps.ext[:expressions][:e])[jy] for jy in mod_ps.ext[:sets][:JY]]

            # Update supply of allowances 
            update_supply!(results["Ind"]["e"][iter,:]+results["PS"]["e"][iter,:],ETS,data,scenario_overview_row)
            results["s"][iter,:] = ETS["S"]

            # Imbalances
            ADMM["Imbalances"]["ETS"][iter,:] = results["s"][iter,:]-results["Ind"]["b"][iter,:]-results["PS"]["b"][iter,:]
            if iter > 1
                ADMM["Imbalances"]["MSR"][iter,:] = results["s"][iter,:]-results["s"][iter-1,:]
            end

            # Primal residuals
            ADMM["Residuals"]["Primal"][iter] = sqrt(sum(ADMM["Imbalances"]["ETS"][iter,:].^2))+sqrt(sum(ADMM["Imbalances"]["MSR"][iter,:].^2))

            # Dual residuals
            if iter > 1  
                ADMM["Residuals"]["Dual"][iter] = sqrt(sum(ADMM["ρ"][iter]*((results["Ind"]["b"][iter,:]-ADMM["Imbalances"]["ETS"][iter,:]/2) - (results["Ind"]["b"][iter-1,:]-ADMM["Imbalances"]["ETS"][iter-1,:]/2)).^2)) 
                                                + sqrt(sum(ADMM["ρ"][iter]*((results["PS"]["b"][iter,:]-ADMM["Imbalances"]["ETS"][iter,:]/2) - (results["PS"]["b"][iter-1,:]-ADMM["Imbalances"]["ETS"][iter-1,:]/2)).^2))
            end

            # Price updates 
            results["λ"][iter+1,:] = results["λ"][iter,:]-mod_ind.ext[:parameters][:ρ]*ADMM["Imbalances"]["ETS"][iter,:]  

            # ρ-updates following Boyd et al. (2010), limited to 0.1 - 0.001 range to avoid divergence
            if ADMM["Residuals"]["Primal"][iter] > 10*ADMM["Residuals"]["Dual"][iter]
                ADMM["ρ"][iter+1] = minimum([10*data["ADMM"]["rho"],2*ADMM["ρ"][iter]])
            elseif ADMM["Residuals"]["Dual"][iter] > 10*ADMM["Residuals"]["Primal"][iter]
                ADMM["ρ"][iter+1] = maximum([0.1*data["ADMM"]["rho"],1/2*ADMM["ρ"][iter]])
            end

            # Progress bar
            set_description(iterations, string(@sprintf("Primal %.3f -- Dual %.3f -- ρ %.3f",  ADMM["Residuals"]["Primal"][iter], ADMM["Residuals"]["Dual"][iter],ADMM["ρ"][iter])))

            # Check convergence: primal and dual satisfy tolerance OR primal remains high and dual near zero. 
            # In case the latter happens, the solution needs to be inspected to ensure this can be interpreted as an equilibrium.
            if (ADMM["Residuals"]["Primal"][iter] <= ADMM["Tolerance"] && ADMM["Residuals"]["Dual"][iter] <= ADMM["Tolerance"]) || (ADMM["Residuals"]["Dual"][iter] <= 0.0001 && iter > 1 && ADMM["Residuals"]["Primal"][iter] <= 100*ADMM["Tolerance"])
                convergence = 1
            end
            ADMM["n_iter"] = copy(iter)
        end
    end
end

# calculate equilibrium 
define_results!(data,results,ADMM) # initialize structure of results 
ADMM["walltime"] = @elapsed ADMM!(results,ADMM,ETS,mod_ind,scenario_overview_row)  

# if the follwoing conditions hold, these runs require calibrating the beta parameter in the MACC
while abs(results["λ"][ADMM["n_iter"],3]-data["P_2019"]) > data["tolerance_calibration"] && scenario_overview_row[:ref_scen_number] == scen_number
    # Calibration β - new estimate:
    println(string("Calibration error 2019 EUA prices :", abs(results["λ"][ADMM["n_iter"],3]-data["P_2019"])," €/tCO2")) 
    println("Update β and retry ...")
    println(string("        "))

    mod_ind.ext[:parameters][:β] = copy(mod_ind.ext[:parameters][:β])*1/(1+(results["λ"][ADMM["n_iter"],3]-data["P_2019"])/data["P_2019"])

    # Calculate equilibrium with new estimate beta
    define_results!(data,results,ADMM) # initialize structure of results 
    ADMM["walltime"] =  @elapsed ADMM!(results,ADMM,ETS,mod_ind,scenario_overview_row)      
end

if ADMM["Residuals"]["Primal"][ADMM["n_iter"]] > ADMM["Tolerance"]
    println(string("Warning: primal residual does not satisfy tolerance. Inspect results to confirm this in an equilibrium solution."))
    println(string("        "))
end
println(string("Required iterations: ",ADMM["n_iter"]))
println(string("Required walltime: ",ADMM["walltime"], " seconds"))
println(string("        "))
println(string("RP: ",  ADMM["Residuals"]["Primal"][ADMM["n_iter"]], " -- Tolerance: ",ADMM["Tolerance"]))
println(string("RD: ",  ADMM["Residuals"]["Dual"][ADMM["n_iter"]], " -- Tolerance: ",ADMM["Tolerance"]))
println(string("        "))

## 5. Postprocessing and save results -> needs to be updated 
# Save results
function SaveResults!(mod_ind::Model,ETS::Dict,data::Dict,ADMM::Dict,results::Dict,scenario_overview_row::DataFrameRow) 
    # Year in which the waterbed is sealed
    Years = range(2017,stop=2017+data["nyears"]-1)
    WBseal = Years[findfirst(ETS["TNAC"] .< data["TNAC_MAX"])]
    # Waterbed leakage 
    if scenario_overview_row[:op_dem] != 0  # overlapping policy that affects emissions
        # Reference results
        ref_results = CSV.read(joinpath(@__DIR__,"Results",string("Scenario_",scenario_overview_row[:ref_scen_number],".csv")),DataFrame;delim=";")
        # WBL 
        WBL =(sum(ETS["C"])-sum(ref_results[:Cancellation]))/(scenario_overview_row[:op_dem]*(scenario_overview_row[:stop_op]-scenario_overview_row[:start_op]+1))
        # Hypothetical emission profile from reference scenario, corrected for overlapping policy
        ref_e = ref_results[:Emissions] 
        for y in range(scenario_overview_row[:start_op]-2016, stop=scenario_overview_row[:stop_op]-2016)
            ref_e[y] = ref_e[y] - scenario_overview_row[:op_dem]
        end        
        # ETS: compute hypothetical actions MSR under hypothetical emission profile 
        ref_ETS = Dict() 
        define_ETS_parameters!(ref_ETS,data,scenario_overview[scen_number,:])       
        update_supply!(ref_e,ref_ETS,data,scenario_overview[scen_number,:]) 
        # Calculate direct waterbed leakage: change in supply of allowances as result of overlapping policy
        dirWBL = (sum(ref_ETS["C"])-sum(ref_results[:Cancellation]))/(scenario_overview_row[:op_dem]*(scenario_overview_row[:stop_op]-scenario_overview_row[:start_op]+1))
        # Calculate indirect waterbed leakage 
        indirWBL = WBL - dirWBL
    elseif scenario_overview_row[:op_supply] != 0 
        # Reference results
        ref_results = CSV.read(joinpath(@__DIR__,"Results",string("Scenario_",scenario_overview_row[:ref_scen_number],".csv")),DataFrame;delim=";")
        # WBL 
        WBL = - (sum(ETS["C"])-sum(ref_results[:Cancellation]))/(scenario_overview_row[:op_supply]*(scenario_overview_row[:stop_op]-scenario_overview_row[:start_op]+1))  
        # ETS: compute hypothetical actions MSR under reference emission profile 
        ref_ETS = Dict() 
        define_ETS_parameters!(ref_ETS,data,scenario_overview[scen_number,:])       
        update_supply!(ref_results[:Emissions],ref_ETS,data,scenario_overview[scen_number,:]) 
        # Calculate direct waterbed leakage: change in supply of allowances as result of overlapping policy
        dirWBL = - (sum(ref_ETS["C"])-sum(ref_results[:Cancellation]))/(scenario_overview_row[:op_supply]*(scenario_overview_row[:stop_op]-scenario_overview_row[:start_op]+1))
        # Calculate indirect waterbed leakage 
        indirWBL = WBL - dirWBL
    else
        WBL = NaN
        dirWBL = NaN
        indirWBL = NaN
    end

    # Vector with aggregate output
    vector_output = [scen_number; ADMM["n_iter"]; ADMM["walltime"]; ADMM["Residuals"]["Primal"][ADMM["n_iter"]]; ADMM["Residuals"]["Dual"][ADMM["n_iter"]]; mod_ind.ext[:parameters][:β]; results["λ"][ADMM["n_iter"],5]; sum(results["e"][ADMM["n_iter"],4:end]); sum(ETS["C"]); WBseal; WBL; dirWBL; indirWBL]
    CSV.write(joinpath(@__DIR__,"overview_results.csv"), DataFrame(reshape(vector_output,1,:),:auto), delim=";",append=true)

    # Equilibrium in ETS year-by-year 
    mat_output = [Years ETS["CAP"] ETS["S"] sum(ETS["C"][:,:],dims=2) ETS["MSR"][:,12] ETS["TNAC"] results["e"][ADMM["n_iter"],:] results["λ"][ADMM["n_iter"],:] results["b"][ADMM["n_iter"],:]]
    CSV.write(joinpath(@__DIR__,"Results",string("Scenario_",scen_number,".csv")), DataFrame(mat_output,:auto), delim=";",header=["Year";"CAP";"Supply";"Cancellation";"MSR";"TNAC";"Emissions";"EUAprice";"EUAs"])
end

SaveResults!(mod_ind,ETS,data,ADMM,results,scenario_overview_row) 

## 6. Plotting (temporary)
I = mod_ps.ext[:sets][:I]
IV = mod_ps.ext[:sets][:IV]
ID = mod_ps.ext[:sets][:ID]
JH = mod_ps.ext[:sets][:JH]
JD = mod_ps.ext[:sets][:JD]
JY = mod_ps.ext[:sets][:JY]
JY = mod_ps.ext[:sets][:JY]
W = mod_ps.ext[:parameters][:W]
GF = mod_ps.ext[:parameters][:GF] # growth factor of the demand
D = mod_ps.ext[:timeseries][:D] # demand

cap = value.(mod_ps.ext[:variables][:cap])
g = value.(mod_ps.ext[:variables][:g])
ens = value.(mod_ps.ext[:variables][:ens])
curt = value.(mod_ps.ext[:expressions][:curt])
λ = dual.(mod_ps.ext[:constraints][:dem_balance])

e_ps = [value.(mod_ps.ext[:expressions][:e][jy]) for jy in JY]
e_ind = mod_ind.ext[:parameters][:e]

demand_yearly = [sum((1+GF[jy])*D[jh,jd] for jh in JH, jd in JD) for jy in JY]
g_yearly = Dict(i => [sum(W[jd]*g[i,jh,jd,jy] for jh in JH, jd in JD) for jy in JY] for i in I)
ens_yearly = [sum(W[jd]*ens[jh,jd,jy] for jh in JH, jd in JD) for jy in JY]
curt_yearly =  [sum(W[jd]*curt[iv,jh,jd,jy] for jh in JH, jd in JD, iv in IV) for jy in JY]
fuel_shares_yearly = Dict(i => [g_yearly[i][jy]/demand_yearly[jy] for jy in JY] for i in I)
cap_yearly = Dict(i => [cap[i,jy] for jy in JY] for i in I)

using Plots
# Capacity
p1= plot(JY,zeros(data["nyears"],1),label="",xlabel="Years [-]",ylabel="New capacity [GW]",legend=:outertopright, linewidth = 2);
for i in I
plot!(JY,cap_yearly[i],label=i, xlabel="Years [-]",ylabel="New capacity [GW]",legend=:outertopright, linewidth = 2);
end
# Generation
p2= plot(JY,zeros(data["nyears"],1),label="",xlabel="Years [-]",ylabel="Generation [GWh]",legend=:outertopright, linewidth = 2);
for i in I
plot!(JY,fuel_shares_yearly[i],label=i, xlabel="Years [-]",ylabel="Generation [GWh]",legend=:outertopright, linewidth = 2);
end
plot(p1, p2, layout = (2, 1))
plot!(size=(600,800))

# EUA prices
p0= plot(JY,mod_ps.ext[:parameters][:λ_EUA],label="",xlabel="Years [-]",ylabel="EUA prices [euro/tCO2]",legend=:outertopright, linewidth = 2);
# Emissions
p1= plot(JY,e_ind,label="Industry",xlabel="Years [-]",ylabel="Emissions [MtCO2]",legend=:outertopright, linewidth = 2);
plot!(JY,e_ps,label="Power sector",xlabel="Years [-]",ylabel="Emissions [MtCO2]",legend=:outertopright, linewidth = 2);
# Supply 
p2= plot(JY,ETS["S"],label="Supply",xlabel="Years [-]",ylabel="EUAs [MtCO2]",legend=:outertopright, linewidth = 2);
plot!(JY,ETS["CAP"],label="Cap",xlabel="Years [-]",ylabel="EUAs [MtCO2]",legend=:outertopright, linewidth = 2);
plot(p0, p1, p2, layout = (3, 1))
plot!(size=(600,800))

println("Postprocessing & save results: done")
println("   ")

# end # end for loop over scenarios

println(string("##############################################################################################"))
