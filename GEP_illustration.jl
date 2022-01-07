## Generation expansion planning problem
# author: Kenneth Bruninx, Sebastian Gonzato
# last update: October 26, 2020
# description: example of a simple generation expansion plannign problem

## Step 0: Activate environment - ensure consistency accross computers
using Pkg
Pkg.activate(@__DIR__) # @__DIR__ = directory this script is in
Pkg.instantiate() # If a Manifest.toml file exist in the current project, download all the packages declared in that manifest. Else, resolve a set of feasible packages from the Project.toml files and install them.

##  Step 1: input data
using CSV
using DataFrames
using YAML

data = YAML.load_file(joinpath(@__DIR__, "data_gep.yaml"))
ts = CSV.read(joinpath(@__DIR__, "Profiles_12_reprdays.csv"), DataFrame)
repr_days = CSV.read(joinpath(@__DIR__, "Weights_12_reprdays.csv"), DataFrame)

## Step 2: create model & pass data to model
using JuMP
using Gurobi
m = Model(optimizer_with_attributes(Gurobi.Optimizer))

# Step 2a: create sets
function define_sets!(m::Model, data::Dict)
    # create dictionary to store sets
    m.ext[:sets] = Dict()

    # define the sets
    m.ext[:sets][:JH] = 1:data["nTimesteps"] # Timesteps
    m.ext[:sets][:JD] = 1:data["nReprDays"] # Representative days
    m.ext[:sets][:ID] = [id for id in keys(data["dispatchableGenerators"])] # dispatchable generators
    m.ext[:sets][:IV] = [iv for iv in keys(data["variableGenerators"])] # variable generators
    m.ext[:sets][:I] = union(m.ext[:sets][:ID], m.ext[:sets][:IV]) # all generators

    # return model
    return m
end

# Step 2b: add time series
function process_time_series_data!(m::Model, data::Dict, ts::DataFrame)
    # extract the relevant sets
    IV = m.ext[:sets][:IV] # variable generators
    JH = m.ext[:sets][:JH] # Time steps
    JD = m.ext[:sets][:JD] # Days

    # create dictionary to store time series
    m.ext[:timeseries] = Dict()
    m.ext[:timeseries][:AF] = Dict()

    # add time series to dictionary
    m.ext[:timeseries][:D] = [ts.Load[jh+data["nTimesteps"]*(jd-1)] for jh in JH, jd in JD]
    m.ext[:timeseries][:AF][IV[1]] = [ts.LFW[jh+data["nTimesteps"]*(jd-1)] for jh in JH, jd in JD]
    m.ext[:timeseries][:AF][IV[2]] = [ts.LFS[jh+data["nTimesteps"]*(jd-1)] for jh in JH, jd in JD]

    # return model
    return m
end

# step 2c: process input parameters
function process_parameters!(m::Model, data::Dict, repr_days::DataFrame)
    # extract the sets you need
    I = m.ext[:sets][:I]
    ID = m.ext[:sets][:ID]
    IV = m.ext[:sets][:IV]

    # generate a dictonary "parameters"
    m.ext[:parameters] = Dict()

    # input parameters
    αCO2 = m.ext[:parameters][:αCO2] = data["CO2Price"] # EUR/ton
    m.ext[:parameters][:VOLL] = data["VOLL"] # VOLL
    r = m.ext[:parameters][:discountrate] = data["discountrate"] # discount rate
    m.ext[:parameters][:W] = repr_days.Weights # weights of each representative day

    d = merge(data["dispatchableGenerators"],data["variableGenerators"])
    # variable costs
    β = m.ext[:parameters][:β] = Dict(i => d[i]["fuelCosts"] for i in ID) # EUR/MWh
    δ = m.ext[:parameters][:δ] = Dict(i => d[i]["emissions"] for i in ID) # ton/MWh
    m.ext[:parameters][:VC] = Dict(i => β[i]+αCO2*δ[i] for i in ID) # variable costs - EUR/MWh

    # investment costs
    OC = m.ext[:parameters][:OC] = Dict(i => d[i]["OC"] for i in I) # EUR/MW
    LifeTime = m.ext[:parameters][:LT] = Dict(i => d[i]["lifetime"] for i in I) # years
    m.ext[:parameters][:IC] = Dict(i => r*OC[i]/(1-(1+r).^(-LifeTime[i])) for i in I) # EUR/MW/y

    # legacy capacity
    m.ext[:parameters][:LC] = Dict(i => d[i]["legcap"] for i in I) # MW

    # return model
    return m
end

# call functions
define_sets!(m, data)
process_time_series_data!(m, data, ts)
process_parameters!(m, data, repr_days)

## Step 3: construct your model
# Greenfield GEP - single year (Lecture 3 - slide 25, but based on representative days instead of full year)
function build_greenfield_1Y_GEP_model!(m::Model)
    # Clear m.ext entries "variables", "expressions" and "constraints"
    m.ext[:variables] = Dict()
    m.ext[:expressions] = Dict()
    m.ext[:constraints] = Dict()

    # Extract sets
    I = m.ext[:sets][:I]
    ID = m.ext[:sets][:ID]
    IV = m.ext[:sets][:IV]
    JH = m.ext[:sets][:JH]
    JD = m.ext[:sets][:JD]

    # Extract time series data
    D = m.ext[:timeseries][:D] # demand
    AF = m.ext[:timeseries][:AF] # avaiability factors

    # Extract parameters
    VOLL = m.ext[:parameters][:VOLL] # VOLL
    VC = m.ext[:parameters][:VC] # variable costs
    IC = m.ext[:parameters][:IC] # investment costs
    W = m.ext[:parameters][:W] # weights

    # Create variables
    cap = m.ext[:variables][:cap] = @variable(m, [i=I], lower_bound=0, base_name="capacity")
    g = m.ext[:variables][:g] = @variable(m, [i=I,jh=JH,jd=JD], lower_bound=0, base_name="generation")
    ens =  m.ext[:variables][:ens] = @variable(m, [jh=JH,jd=JD], lower_bound=0, base_name="load_shedding")

    # Create affine expressions (= linear combinations of variables)
    curt = m.ext[:expressions][:curt] = @expression(m, [i=IV,jh=JH,jd=JD],
        AF[i][jh,jd] * cap[i] - g[i,jh,jd]
    )

    # Formulate objective 1a
    m.ext[:objective] = @objective(m, Min,
        + sum(IC[i]*cap[i] for i in I)
        + sum(W[jd]*VC[i]*g[i,jh,jd] for i in ID, jh in JH, jd in JD)
        + sum(W[jd]*ens[jh,jd]*VOLL for jh in JH, jd in JD)
    )

    # Formulate constraints
    m.ext[:constraints][:con2a] = @constraint(m, [jh=JH,jd=JD],
        sum(g[i,jh,jd] for i in I) == D[jh,jd] - ens[jh,jd]
    )
    m.ext[:constraints][:con2c] = @constraint(m, [jh=JH,jd=JD],
        ens[jh,jd] <= D[jh,jd]
    )

    # 3a1 - renewables
    m.ext[:constraints][:con3a1res] = @constraint(m, [i=IV,jh=JH,jd=JD],
        g[i,jh,jd] <=  AF[i][jh,jd]*cap[i]
    )

    # 3a1 - conventional
    m.ext[:constraints][:con3a1conv] = @constraint(m, [i=ID,jh=JH,jd=JD],
        g[i,jh,jd] <= cap[i]
    )

    return m
end

# Brownfield GEP - single year
function build_brownfield_1Y_GEP_model!(m::Model)
    # start from the greenfield model
    m = build_greenfield_1Y_GEP_model!(m::Model)

    # extract sets
    ID = m.ext[:sets][:ID]
    IV = m.ext[:sets][:IV]
    JH = m.ext[:sets][:JH]
    JD = m.ext[:sets][:JD]

    # extract parameters
    LC = m.ext[:parameters][:LC]

    # extract time series data
    AF = m.ext[:timeseries][:AF] # avaiability factors

    # extract variables
    g = m.ext[:variables][:g]
    cap = m.ext[:variables][:cap]

    # remove the constraints that need to be changed:
    for iv in IV, jh in JH, jd in JD
        delete(m,m.ext[:constraints][:con3a1res][iv,jh,jd])
    end
    for id in ID, jh in JH, jd in JD
        delete(m,m.ext[:constraints][:con3a1conv][id,jh,jd])
    end

    # define new constraints
    # 3a1 - renewables
    m.ext[:constraints][:con3a1res] = @constraint(m, [i=IV,jh=JH,jd=JD],
        g[i,jh,jd] <=  AF[i][jh,jd]*(cap[i]+LC[i])
    )

    # 3a1 - conventional
    m.ext[:constraints][:con3a1conv] = @constraint(m, [i=ID,jh=JH,jd=JD],
        g[i,jh,jd] <= (cap[i]+LC[i])
    )

    return m
end

# build your model (activate 1 option)
build_greenfield_1Y_GEP_model!(m)
# build_brownfield_1Y_GEP_model!(m)

## Step 4: solve
optimize!(m)

# check termination status
print(
    """

    Termination status: $(termination_status(m))

    """
)


# print some output
@show value(m.ext[:objective])
@show value.(m.ext[:variables][:cap])

## Step 5: interpretation
using Plots
using StatsPlots

# sets
JH = m.ext[:sets][:JH]
JD = m.ext[:sets][:JD]
I = m.ext[:sets][:I]

# parameters
D = m.ext[:timeseries][:D]
W = m.ext[:parameters][:W]
LC = m.ext[:parameters][:LC]

# variables/expressions
cap = value.(m.ext[:variables][:cap])
g = value.(m.ext[:variables][:g])
ens = value.(m.ext[:variables][:g])
curt = value.(m.ext[:expressions][:curt])
λ = dual.(m.ext[:constraints][:con2a])

# create arrays for plotting
λvec = [λ[jh,jd]/W[jd] for jh in JH, jd in JD]
gvec = [g[i,jh,jd] for  i in I, jh in JH, jd in JD]
capvec = [cap[i] for  i in I]

# Select day for which you'd like to plot: 
jd = 2
# electricity price price
p1 = plot(JH,λvec[:,jd], xlabel="Timesteps [-]",ylabel="λ [EUR/MWh]",label="λ [EUR/MWh]",legend=:outertopright);
# dispatch
p2 = groupedbar(transpose(gvec[:,:,jd]),label=["Mid" "Base" "Peak" "Wind" "Solar"],bar_position = :stack,legend=:outertopright);
plot!(JH,D[:,jd], label="Demand", xlabel="Timesteps [-]",ylabel="Generation [MWh]",legend=:outertopright, linewidth = 2);
# capacity
p3 = groupedbar(transpose(capvec), label=["Mid" "Base" "Peak" "Wind" "Solar"],xlabel="Technology [-]",ylabel="New capacity [MW]",bar_position = :stack,legend=:outertopright);
plot(p1, p2, p3, layout = (3, 1))
plot!(size=(600,800))
