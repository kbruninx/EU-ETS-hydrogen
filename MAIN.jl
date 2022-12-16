## Topic: EU ETS, MSR and overlapping policies
# Author: Kenneth Bruninx
# Last update: November 2022

## 0. Set-up code
# HPC or not?
HPC = "DelftBlue" # NA, DelftBlue or ThinKing

# Home directory
const home_dir = @__DIR__

if HPC == "DelftBlue"  # only for running this on DelftBlue
    ENV["GRB_LICENSE_FILE"] = "./Hpc/gurobi.lic"
    ENV["GUROBI_HOME"] = "./scratch/kbruninx/gurobi950/linux64"
    println(string("Number of threads: ", Threads.nthreads()))
end

if HPC == "ThinKing"  # only for running this on VSC
    # ENV["GRB_LICENSE_FILE"] = " "
    # ENV["GUROBI_HOME"] = " "
end

# Include packages 
using JuMP, Gurobi # Optimization packages
using DataFrames, CSV, YAML, DataStructures # dataprocessing
using ProgressBars, Printf # progress bar
using TimerOutputs # profiling 
using JLD2
using Base.Threads: @spawn 
using Base: split
using ArgParse # Parsing arguments from the command line

# Gurobi environment to suppress output
println("Define Gurobi environment...")
println("        ")
const GUROBI_ENV = Gurobi.Env()
# set parameters:
GRBsetparam(GUROBI_ENV, "OutputFlag", "0")   
GRBsetparam(GUROBI_ENV, "Threads", "4")   
GRBsetparam(GUROBI_ENV, "TimeLimit", "300")  # will only affect solutions if you're selecting representative days  
println("        ")

# Include functions
include(joinpath(home_dir,"Source","define_common_parameters.jl"))
include(joinpath(home_dir,"Source","define_H2S_parameters.jl"))
include(joinpath(home_dir,"Source","define_ps_parameters.jl"))
include(joinpath(home_dir,"Source","define_ind_parameters.jl"))
include(joinpath(home_dir,"Source","define_ETS_parameters.jl"))
include(joinpath(home_dir,"Source","define_EOM_parameters.jl"))
include(joinpath(home_dir,"Source","define_REC_parameters.jl"))
include(joinpath(home_dir,"Source","define_H2_parameters.jl"))
include(joinpath(home_dir,"Source","define_H2CN_prod_parameters.jl"))
include(joinpath(home_dir,"Source","define_H2CN_cap_parameters.jl"))
include(joinpath(home_dir,"Source","define_NG_parameters.jl"))
include(joinpath(home_dir,"Source","build_ind_agent.jl"))
include(joinpath(home_dir,"Source","build_ps_agent.jl"))
include(joinpath(home_dir,"Source","build_H2S_agent.jl"))
include(joinpath(home_dir,"Source","define_results.jl"))
include(joinpath(home_dir,"Source","ADMM.jl"))
include(joinpath(home_dir,"Source","ADMM_subroutine.jl"))
include(joinpath(home_dir,"Source","update_ind_emissions.jl"))
include(joinpath(home_dir,"Source","solve_ind_agent.jl"))
include(joinpath(home_dir,"Source","solve_ps_agent.jl"))
include(joinpath(home_dir,"Source","solve_H2S_agent.jl"))
include(joinpath(home_dir,"Source","update_supply.jl"))
include(joinpath(home_dir,"Source","update_rho.jl"))
include(joinpath(home_dir,"Source","save_results.jl"))

# Data common to all scenarios data 
temp_data = YAML.load_file(joinpath(home_dir,"Input","overview_data.yaml"))
ts = CSV.read(joinpath(home_dir,"Input","timeseries.csv"),delim=";",DataFrame)
if isfile(joinpath(home_dir,"Input",string("output_",temp_data["General"]["nReprDays"],"_repr_days"),"decision_variables_short.csv"))
    repr_days = rightjoin(CSV.read(joinpath(home_dir,"Input",string("output_",temp_data["General"]["nReprDays"],"_repr_days"),"decision_variables_short.csv"),delim=",",DataFrame), CSV.read(joinpath(home_dir,"Input",string("output_",temp_data["General"]["nReprDays"],"_repr_days"),"weight_day_month.csv"),delim=",",DataFrame),on= :periods)
else    
    using RepresentativePeriodsFinder # representative day finder  - https://ucm.pages.gitlab.kuleuven.be/representativeperiodsfinder.jl/
    # select representative days from the timeseries.csv" file
    println("Selecting representative days (may take up to 300 seconds)... ")
    config_file = joinpath(home_dir,"Input","config_file_repr_days.yaml") # contains general info/specification of representative day finder package
    pf = PeriodsFinder(config_file; populate_entries=true)
    # specific settings to select representative days
    pf.config["method"]["optimization"]["binary_ordering"] = true
    pf.config["method"]["options"]["representative_periods"] = temp_data["General"]["nReprDays"]
    pf.config["results"]["result_dir"] = string("output_",temp_data["General"]["nReprDays"],"_repr_days")
    delete!(pf.config["method"]["optimization"], "time_series_error")
    pf = find_representative_periods(pf, optimizer=optimizer_with_attributes(() -> Gurobi.Optimizer(GUROBI_ENV)))
    repr_days = CSV.read(joinpath(home_dir,"Input",string("output_",temp_data["General"]["nReprDays"],"_repr_days"),"decision_variables_short.csv"),delim=",",DataFrame)
    
    println("        ")
    println("Done!")
    println("         ")
    println("Computing mapping representative days to months...")

    # Create syntethic time series - https://ucm.pages.gitlab.kuleuven.be/representativeperiodsfinder.jl/examples/days_re_ordering/ 
    pf = PeriodsFinder(config_file; populate_entries=true)
    # specific settings to create syntethic time series
    pf.config["method"]["optimization"]["binary_ordering"] = false
    pf.config["method"]["options"]["representative_periods"] = temp_data["General"]["nReprDays"]
    pf.config["results"]["result_dir"] = string("output_",temp_data["General"]["nReprDays"],"_repr_days")
    delete!(pf.config["method"]["optimization"], "duration_curve_error")
    pf.config["method"]["options"]["mandatory_periods"] = repr_days[!,:periods]
    pf.u = zeros(Bool, pf.config["method"]["options"]["total_periods"])
    pf.u[repr_days[!,:periods]] .= 1 
    find_representative_periods(pf, optimizer=optimizer_with_attributes(() -> Gurobi.Optimizer(GUROBI_ENV)), reset=false);
    write_out_synthetic_timeseries(pf)
    # compute mapping day -> month
    n_days_month = [0 31 28 31 30 31 30 31 31 30 31 30 31 0] # to compute weights of day -> month
    weight_day_month = [repr_days[!,:periods] [sum(pf.v[1+sum(n_days_month[1:m]):sum(n_days_month[1:m+1]),jd]) for jd=1:temp_data["General"]["nReprDays"],m=1:12]]
    CSV.write(joinpath(home_dir,"Input",string("output_",temp_data["General"]["nReprDays"],"_repr_days"),"weight_day_month.csv"),DataFrame(weight_day_month,:auto),delim=",",header=[:periods, :January, :February, :March, :April, :May, :June, :July, :August, :September, :October, :November, :December])
    repr_days = rightjoin(repr_days, CSV.read(joinpath(home_dir,"Input",string("output_",temp_data["General"]["nReprDays"],"_repr_days"),"weight_day_month.csv"),delim=",",DataFrame),on= :periods)
    println("        ")
    println("Done!")
    println("        ")
end

# Overview scenarios
scenario_overview = CSV.read(joinpath(home_dir,"overview_scenarios.csv"),DataFrame,delim=";")
sensitivity_overview = CSV.read(joinpath(home_dir,"overview_sensitivity.csv"),DataFrame,delim=";") 

# Create file with results 
# add column for sensitivity analsysis
if isfile(joinpath(home_dir,string("overview_results_",temp_data["General"]["nReprDays"],"_repr_days.csv"))) != 1
    CSV.write(joinpath(home_dir,string("overview_results_",temp_data["General"]["nReprDays"],"_repr_days.csv")),DataFrame(),delim=";",header=["scen_number";"sensitivity";"n_iter";"walltime";"PrimalResidual_ETS";"PrimalResidual_MSR";"PrimalResidual_EOM";"PrimalResidual_REC";"PrimalResidual_H2";"PrimalResidual_H2CN_prod";"PrimalResidual_H2CN_cap"; "DualResidual_ETS"; "DualResidual_EOM";"DualResidual_REC";"DualResidual_H2";"DualResidual_H2CN_prod";"DualResidual_H2CN_cap";"Beta";"EUA_2021";"CumulativeEmissions";"TotalCost"])
end

# Create folder for results
if isdir(joinpath(home_dir,string("Results_",temp_data["General"]["nReprDays"],"_repr_days"))) != 1
    mkdir(joinpath(home_dir,string("Results_",temp_data["General"]["nReprDays"],"_repr_days")))
end

# Scenario number 
if HPC == "DelftBlue" || HPC == "ThinKing"
   function parse_commandline()
       s = ArgParseSettings()
       @add_arg_table! s begin
           "--start_scen"
               help = "Enter the number of the first scenario here"
               arg_type = Int
               default = 1
            "--stop_scen"
               help = "Enter the number of the last scenario here"
               arg_type = Int
               default = 1
       end
       return parse_args(s)
   end
   # Simulation number as argument:
   dict_sim_number =  parse_commandline()
   start_scen = dict_sim_number["start_scen"]
   stop_scen = dict_sim_number["stop_scen"]
else
    # Range of scenarios to be simulated
    start_scen = 1
    stop_scen = 1
end

# scen_number = 1
for scen_number in range(start_scen,stop=stop_scen,step=1)

println("    ")
println(string("######################                  Scenario ",scen_number,"                 #########################"))

## 1. Read associated input for this simulation
scenario_overview_row = Dict(pairs(scenario_overview[scen_number,:])) # create dict from dataframe
scenario_definition = Dict("scenario" => Dict([String(collect(keys(scenario_overview_row))[x]) => collect(values(scenario_overview_row))[x] for x = 1:length(collect(values(scenario_overview_row)))]))  # Keys from Symbol to String
data = YAML.load_file(joinpath(home_dir,"Input","overview_data.yaml")) # reload data to avoid previous sensitivity analysis affected data
data = merge(data,scenario_definition)

if data["scenario"]["Sens_analysis"] == "YES" && data["scenario"][:ref_scen_number] != scen_number
    numb_of_sens = length((sensitivity_overview[!,:Parameter]))
else
    numb_of_sens = 0 
end    
# sens_number = 1 
for sens_number in range(1,stop=numb_of_sens+1,step=1) 
if sens_number >= 2
    println("    ") 
    println(string("#                                  Sensitivity ",sens_number-1,"                                      #"))
    parameter = split(sensitivity_overview[sens_number-1,:Parameter])
    if length(parameter) == 2
        data[parameter[1]][parameter[2]] = sensitivity_overview[sens_number-1,:Scaling]*data[parameter[1]][parameter[2]]
    elseif length(parameter) == 3
        data[parameter[1]][parameter[2]][parameter[3]] = sensitivity_overview[sens_number-1,:Scaling]*data[parameter[1]][parameter[2]][parameter[3]]
    else
        printnl("warning! Sensitivity analysis is not well defined!")
    end
end

# write final set-up to yaml 
if sens_number == 1
    YAML.write_file(joinpath(home_dir,string("Results_",data["General"]["nReprDays"],"_repr_days"),string("Scenario_",data["scenario"]["scen_number"],"_ref.yaml")),data)
else    
    YAML.write_file(joinpath(home_dir,string("Results_",data["General"]["nReprDays"],"_repr_days"),string("Scenario_",data["scenario"]["scen_number"],"_",sensitivity_overview[sens_number-1,:remarks],".yaml")),data)
end

println("    ")
println("Including all required input data: done")
println("   ")

## 2. Initiate models for representative agents 
agents = Dict()
agents[:ps] = [id for id in keys(data["PowerSector"])] 
agents[:h2s] = [id for id in keys(data["HydrogenSector"])]
agents[:ind] = ["Ind"] 
agents[:all] = union(agents[:ps],agents[:h2s],agents[:ind])   
# Different markets - to be completed based on the agents
agents[:eom] = []                  
agents[:ets] = []                  
agents[:rec] = []      
agents[:h2] = []             
agents[:h2cn_prod] = []         
agents[:h2cn_cap] = []   
agents[:ng] = []                    
mdict = Dict(i => Model(optimizer_with_attributes(() -> Gurobi.Optimizer(GUROBI_ENV))) for i in agents[:all])

## 3. Define parameters for markets and representative agents
# Parameters/variables ETS 
ETS = Dict()
define_ETS_parameters!(ETS,merge(data["General"],data["ETS"],data["scenario"]))

# Parameters/variables EOM
EOM = Dict()
define_EOM_parameters!(EOM,merge(data["General"],data["EOM"]),ts,repr_days)

# Parameters/variables REC 
REC = Dict()
define_REC_parameters!(REC,merge(data["General"],data["REC"],data["scenario"]),ts,repr_days)

# Parameters/variables incentive scheme carbon neutral hydrogen 
H2CN_prod = Dict()
define_H2CN_prod_parameters!(H2CN_prod,merge(data["General"],data["H2"],data["scenario"]),ts,repr_days)

# Parameters/variables incentive scheme carbon neutral hydrogen production capacity
H2CN_cap = Dict()
define_H2CN_cap_parameters!(H2CN_cap,merge(data["General"],data["scenario"]),ts,repr_days)

# Parameters/variables Hydrogen Market
H2 = Dict()
define_H2_parameters!(H2,merge(data["General"],data["H2"],data["scenario"]),ts,repr_days,H2CN_prod)

# Parameters/variables natural gas market
NG = Dict()
define_NG_parameters!(NG,merge(data["General"],data["NG"]),ts,repr_days)

# Parameters/variables natural gas market
for m in agents[:ind]
    define_common_parameters!(m,mdict[m],merge(data["General"],data["ADMM"],data["Industry"]),ts,repr_days,agents)           # Parameters common to all agents
    define_ind_parameters!(mdict[m],merge(data["General"],data["Industry"],data["ETS"],data["scenario"]))                    # Industry
end
for m in agents[:ps]
    define_common_parameters!(m,mdict[m],merge(data["General"],data["ADMM"],data["PowerSector"][m]),ts,repr_days,agents)     # Parameters common to all agents
    define_ps_parameters!(mdict[m],merge(data["General"],data["PowerSector"][m]),ts,repr_days)                               # Power sector
end
for m in agents[:h2s]
    define_common_parameters!(m,mdict[m],merge(data["General"],data["ADMM"],data["HydrogenSector"][m]),ts,repr_days,agents)  # Parameters common to all agents
    define_H2S_parameters!(mdict[m],merge(data["General"],data["HydrogenSector"][m],data["scenario"]),ts,repr_days,REC)      # Hydrogen sector
end

# Calculate number of agents in each market
ETS["nAgents"] =  length(agents[:ets])
EOM["nAgents"] = length(agents[:eom])
REC["nAgents"] = length(agents[:rec])
H2["nAgents"] =  length(agents[:h2])
H2CN_prod["nAgents"] = length(agents[:h2cn_prod])
H2CN_cap["nAgents"] = length(agents[:h2cn_cap])
NG["nAgents"] = length(agents[:ng])

println("Inititate model, sets and parameters: done")
println("   ")

## 4. Build models
for m in agents[:ind]
    build_ind_agent!(mdict[m])
end
for m in agents[:ps]
    build_ps_agent!(mdict[m])
end
for m in agents[:h2s]
    build_h2s_agent!(mdict[m])
end

println("Build model: done")
println("   ")

## 5. ADMM proces to calculate equilibrium
println("Find equilibrium solution...")
println("   ")
println("(Progress indicators on primal residuals, relative to tolerance: <1 indicates convergence)")
println("   ")

results = Dict()
ADMM = Dict()
TO = TimerOutput()
define_results!(merge(data["General"],data["ADMM"],data["scenario"]),results,ADMM,agents,ETS,EOM,REC,H2,H2CN_prod,H2CN_cap,NG)            # initialize structure of results, only those that will be stored in each iteration
ADMM!(results,ADMM,ETS,EOM,REC,H2,H2CN_prod,H2CN_cap,NG,mdict,agents,data,TO)                                                             # calculate equilibrium 
ADMM["walltime"] =  TimerOutputs.tottime(TO)*10^-9/60                                                                                     # wall time 

# Calibration of industry MACC
while abs(results[ "λ"]["EUA"][end][3]-data["ETS"]["P_calibration"]) > data["Industry"]["tolerance_calibration"] && data["scenario"]["ref_scen_number"] == scen_number
    # Calibration β - new estimate:
    println(string("Calibration error 2021 EUA prices: " , results[ "λ"]["EUA"][end][3]-data["ETS"]["P_calibration"]," €/tCO2"))

    mdict["Ind"].ext[:parameters][:β] = copy(mdict["Ind"].ext[:parameters][:β]*1/(1+(results[ "λ"]["EUA"][end][3]-data["ETS"]["P_calibration"])/data["ETS"]["P_calibration"])^(1/data["scenario"]["gamma"]))

    println(string("Required iterations: ",ADMM["n_iter"]))
    println(string("Required walltime: ",ADMM["walltime"], " minutes"))
    println(string("New estimate for β: ", mdict["Ind"].ext[:parameters][:β]))
    println(string("        "))

    # Calculate equilibrium with new estimate beta
    define_results!(merge(data["General"],data["ADMM"],data["scenario"]),results,ADMM,agents,ETS,EOM,REC,H2,H2CN_prod,H2CN_cap,NG)      # initialize structure of results, only those that will be stored in each iteration
    ADMM!(results,ADMM,ETS,EOM,REC,H2,H2CN_prod,H2CN_cap,NG,mdict,agents,data,TO)                                                       # calculate equilibrium 
    ADMM["walltime"] =  TimerOutputs.tottime(TO)*10^-9/60                                                                               # wall time 
end

println(string("Done!"))
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
println(string("RP H2: ",  ADMM["Residuals"]["Primal"]["H2"][end], " -- Tolerance: ",ADMM["Tolerance"]["H2"]))
println(string("RD H2: ",  ADMM["Residuals"]["Dual"]["H2"][end], " -- Tolerance: ",ADMM["Tolerance"]["H2"]))
println(string("RP H2CN_prod: ",  ADMM["Residuals"]["Primal"]["H2CN_prod"][end], " -- Tolerance: ",ADMM["Tolerance"]["H2CN_prod"]))
println(string("RD H2CN_prod: ",  ADMM["Residuals"]["Dual"]["H2CN_prod"][end], " -- Tolerance: ",ADMM["Tolerance"]["H2CN_prod"]))
println(string("RP H2CN_cap: ",  ADMM["Residuals"]["Primal"]["H2CN_cap"][end], " -- Tolerance: ",ADMM["Tolerance"]["H2CN_cap"]))
println(string("RD H2CN_cap: ",  ADMM["Residuals"]["Dual"]["H2CN_cap"][end], " -- Tolerance: ",ADMM["Tolerance"]["H2CN_cap"]))
println(string("        "))

## 6. Postprocessing and save results 
if sens_number >= 2
    save_results(mdict,EOM,ETS,ADMM,results,merge(data["General"],data["ADMM"],data["H2"],data["scenario"]),agents,sensitivity_overview[sens_number-1,:remarks]) 
    @save joinpath(home_dir,string("Results_",data["General"]["nReprDays"],"_repr_days"),string("Scenario_",data["scenario"]["scen_number"],"_",sensitivity_overview[sens_number-1,:remarks]))
else
    save_results(mdict,EOM,ETS,ADMM,results,merge(data["General"],data["ADMM"],data["H2"],data["scenario"]),agents,"ref") 
    @save joinpath(home_dir,string("Results_",data["General"]["nReprDays"],"_repr_days"),string("Scenario_",data["scenario"]["scen_number"],"_ref"))
end

println("Postprocessing & save results: done")
println("   ")

end # end loop over sensititivity
end # end for loop over scenarios

println(string("##############################################################################################"))