## Topic: EU ETS, MSR and overlapping policies
# Author: Kenneth Bruninx
# Last update: November 2021

## 0. Set-up code
# Range of scenarios to be simulated
start_scen = 3
stop_scen = 4

# Include packages 
using JuMP, Gurobi # Optimization packages
using DataFrames, CSV, YAML, DataStructures # dataprocessing
using ProgressBars, Printf # progress bar
using Plots # visuals
using TimerOutputs # profiling 

# Home directory
const home_dir = @__DIR__

# Gurobi environment to suppress output
const GUROBI_ENV = Gurobi.Env()

# Include functions
include(joinpath(home_dir,"Source","define_common_parameters.jl"))
include(joinpath(home_dir,"Source","define_ps_parameters.jl"))
include(joinpath(home_dir,"Source","define_ind_parameters.jl"))
include(joinpath(home_dir,"Source","define_ETS_parameters.jl"))
include(joinpath(home_dir,"Source","define_EOM_parameters.jl"))
include(joinpath(home_dir,"Source","define_REC_parameters.jl"))
include(joinpath(home_dir,"Source","build_ind_agent.jl"))
include(joinpath(home_dir,"Source","build_ps_agent.jl"))
include(joinpath(home_dir,"Source","define_results.jl"))
include(joinpath(home_dir,"Source","ADMM_v2.jl"))
include(joinpath(home_dir,"Source","update_ind_emissions.jl"))
include(joinpath(home_dir,"Source","solve_ind_agent.jl"))
include(joinpath(home_dir,"Source","solve_ps_agent.jl"))
include(joinpath(home_dir,"Source","update_supply.jl"))
include(joinpath(home_dir,"Source","update_rho.jl"))
include(joinpath(home_dir,"Source","save_results.jl"))
include(joinpath(home_dir,"Source","plot_results.jl"))

# Data common to all scenarios data 
data = YAML.load_file(joinpath(home_dir,"Input","overview_data.yaml"))
ts = CSV.read(joinpath(home_dir,"Input","timeseries_May2018.csv"),delim=",",DataFrame)
repr_days = CSV.read(joinpath(home_dir,"Input",string("Period_definitions_",data["General"]["nReprDays"],".csv")),delim=",",DataFrame)

# Overview scenarios
scenario_overview = CSV.read(joinpath(home_dir,"overview_scenarios.csv"),DataFrame;delim=";")

# Create file with results 
if isfile(joinpath(home_dir,"overview_results.csv")) != 1
    CSV.write(joinpath(home_dir,"overview_results.csv"),DataFrame(),delim=";",header=["scen_number";"n_iter";"walltime";"PrimalResidual_ETS";"PrimalResidual_MSR";"PrimalResidual_EOM";
            "PrimalResidual_REC"; "DualResidual_ETS"; "DualResidual_EOM";"DualResidual_REC";"Beta";"EUA_2021";"CumulativeEmissions";"Cancellation"; "WBseal"; "WBL"; "dirWBL"; "indirWBL"])
end

# Create folder for results
if isdir(joinpath(home_dir,"Results")) != 1
    mkdir(joinpath(home_dir,"Results"))
end

# Scenario number 
for scen_number in range(start_scen,stop=stop_scen,step=1)

println("    ")
println(string("######################                  Scenario ",scen_number,"                 #########################"))
println("    ")
println("Including all required input data: done")
println("   ")

## 1. Read associated input for this simulation
scenario_overview_row = scenario_overview[scen_number,:]

## 2. Initiate models for representative agents 
agents = Dict()
agents[:ps] = [id for id in keys(data["PowerSector"])]
# agents[:hyd] = [id for id in keys(data["HydrogenSector"])]
agents[:ind] = ["Ind"] 
agents[:all] = union(agents[:ps],agents[:ind])  
agents[:eom] = []                  # all agents that participate in the eom
agents[:ets] = []                  # all agents that participate in the eua auctions
agents[:rec] = []                  # all agents that participate in the rec auctions
# agents[:hydm] = []                  # all agents that participate in the rec auctions

mdict = Dict(i => Model(optimizer_with_attributes(() -> Gurobi.Optimizer(GUROBI_ENV))) for i in agents[:all])

## 3. Define parameters for representative agents
for m in agents[:ind]
    define_common_parameters!(m,mdict[m],merge(data["General"],data["ADMM"],data["Industry"]),ts,repr_days,agents,scenario_overview_row)        # Parameters common to all agents
    define_ind_parameters!(mdict[m],merge(data["General"],data["Industry"],data["ETS"]),scenario_overview_row)                                              # Industry
end
for m in agents[:ps]
    define_common_parameters!(m,mdict[m],merge(data["General"],data["ADMM"],data["PowerSector"][m]),ts,repr_days,agents,scenario_overview_row)  # Parameters common to all agents
    define_ps_parameters!(mdict[m],merge(data["General"],data["PowerSector"][m]),ts,repr_days,scenario_overview_row)                            # Power sector
end
# if m in agents[:hyd]
#     define_common_parameters!(m,mdict[m],merge(data["General"],data["ADMM"],data["PowerSector"][m]),ts,repr_days,scenario_overview_row) # Parameters common to all agents
#     define_ps_parameters!(mdict[m],merge(data["General"],data["PowerSector"][m]),ts,repr_days,scenario_overview_row)                    # Power sector
# end

# Parameters/variables ETS 
ETS = Dict()
define_ETS_parameters!(ETS,merge(data["General"],data["ETS"]),scenario_overview[scen_number,:])

# Parameters/variables EOM
EOM = Dict()
define_EOM_parameters!(EOM,merge(data["General"],data["EOM"]),ts,repr_days,scenario_overview[scen_number,:])

# Parameters/variables REC 
REC = Dict()
define_REC_parameters!(REC,merge(data["General"],data["REC"]),ts,repr_days,scenario_overview[scen_number,:])

# # Parameters/variables HydrogenMarket
# REC = Dict()
# define_REC_parameters!(REC,merge(data["General"],data["REC"]),ts,repr_days,scenario_overview[scen_number,:])

println("Inititate model, sets and parameters: done")
println("   ")

## 4. Build models
for m in agents[:all]
    if m in agents[:ind]
        build_ind_agent!(mdict[m])
    end
    if m in agents[:ps]
        build_ps_agent!(mdict[m])
    end
    # if m in agents[:hyd]
    #     build_ps_agent!(mdict[m])
    # end
end

println("Build model: done")
println("   ")

## 5. ADMM proces to calculate equilibrium
results = Dict()
ADMM = Dict()
TO = TimerOutput()
define_results!(merge(data["General"],data["ADMM"]),results,ADMM,agents,ETS,EOM,REC)       # initialize structure of results, only those that will be stored in each iteration
ADMM!(results,ADMM,ETS,EOM,REC,mdict,agents,scenario_overview_row,TO)                      # calculate equilibrium 

# Calibration of industry MACC
while abs(results[ "λ"]["EUA"][end][3]-data["ETS"]["P_2019"]) > data["Industry"]["tolerance_calibration"] && scenario_overview_row[:ref_scen_number] == scen_number
    # Calibration β - new estimate:
    println(string("Calibration error 2019 EUA prices: " , results[ "λ"]["EUA"][end][3]-data["ETS"]["P_2019"]," €/tCO2"))

    mdict["Ind"].ext[:parameters][:β] = copy(mdict["Ind"].ext[:parameters][:β]*1/(1+(results[ "λ"]["EUA"][end][3]-data["ETS"]["P_2019"])/data["ETS"]["P_2019"])^(1/scenario_overview_row[:gamma]))
    
    println(string("New estimate for β: ", mdict["Ind"].ext[:parameters][:β]))
    println(string("        "))

    # Calculate equilibrium with new estimate beta
    define_results!(merge(data["General"],data["ADMM"]),results,ADMM,agents,ETS,EOM,REC)    # initialize structure of results, only those that will be stored in each iteration
    ADMM!(results,ADMM,ETS,EOM,REC,mdict,agents,scenario_overview_row,TO)                   # calculate equilibrium 
end
ADMM["walltime"] =  TimerOutputs.tottime(TO)*10^-9/60                                       # wall time 

println(string("        "))
println(string("Required iterations: ",ADMM["n_iter"]))
println(string("Required walltime: ",ADMM["walltime"], " minutes"))
println(string("        "))
println(string("RP MSR: ",  ADMM["Residuals"]["Primal"]["MSR"][end], " -- Tolerance: ",ADMM["Tolerance"]["ETS"]))
println(string("RP ETS: ",  ADMM["Residuals"]["Primal"]["ETS"][end], " -- Tolerance: ",ADMM["Tolerance"]["ETS"]))
println(string("RD ETS: ",  ADMM["Residuals"]["Dual"]["ETS"][end], " -- Tolerance: ",ADMM["Tolerance"]["ETS"]))
println(string("RP EOM: ",  ADMM["Residuals"]["Primal"]["EOM"][end], " -- Tolerance: ",ADMM["Tolerance"]["EOM"]))
println(string("RD EOM: ",  ADMM["Residuals"]["Dual"]["EOM"][end], " -- Tolerance: ",ADMM["Tolerance"]["EOM"]))
println(string("RP REC: ",  ADMM["Residuals"]["Primal"]["REC"][end], " -- Tolerance: ",ADMM["Tolerance"]["REC"]))
println(string("RD REC: ",  ADMM["Residuals"]["Dual"]["REC"][end], " -- Tolerance: ",ADMM["Tolerance"]["REC"]))
println(string("        "))

## 6. Postprocessing and save results 
save_results(mdict,EOM,ETS,ADMM,results,merge(data["General"],data["ADMM"]),agents,scenario_overview_row) 

println("Postprocessing & save results: done")
println("   ")

end # end for loop over scenarios

println(string("##############################################################################################"))