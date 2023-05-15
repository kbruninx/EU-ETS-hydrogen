function build_cc_agent!(mod::Model)
    # Extract sets
    JH = mod.ext[:sets][:JH]
    JD = mod.ext[:sets][:JD]
    JM = mod.ext[:sets][:JM]
    JY = mod.ext[:sets][:JY]
       
    # Extract parameters
    W = mod.ext[:parameters][:W] # weight of the representative days
    IC = mod.ext[:parameters][:IC] # overnight investment costs
    LEG_CAP = mod.ext[:parameters][:LEG_CAP] # legacy capacity
    CAP_LT = mod.ext[:parameters][:CAP_LT] # lead time on new capacity
    CAP_SV = mod.ext[:parameters][:CAP_SV] # salvage value of new capacity
    DELTA_CAP_MAX = mod.ext[:parameters][:DELTA_CAP_MAX] # max YoY change in new capacity
    A = mod.ext[:parameters][:A] # discount factors
    η_E_CO2 = mod.ext[:parameters][:η_E_CO2] # efficiency E->H2
    ADD_SF = mod.ext[:parameters][:ADD_SF] # efficiency E->H2

    # Decision variables
    cap = mod.ext[:variables][:cap] = @variable(mod, [jy=JY], lower_bound=0, base_name="capacity")
    cc = mod.ext[:variables][:cc] = @variable(mod, [jh=JH,jd=JD,jy=JY], lower_bound=0, base_name="carbon_captured")
    g = mod.ext[:variables][:g] = @variable(mod, [jh=JH,jd=JD,jy=JY], upper_bound=0, base_name="demand_electricity_hydrogen") # note this is defined as a negative number, consumption
    b = mod.ext[:variables][:b] = @variable(mod, [jy=JY], upper_bound=0, base_name="EUA") 
    r_y = mod.ext[:variables][:r_y] = @variable(mod, [jy=JY], upper_bound=0, base_name="REC_y") 
    r_m = mod.ext[:variables][:r_m] = @variable(mod, [jd=JM,jy=JY], upper_bound=0, base_name="REC_m") # added for consistency with H2/code, not used
    r_d = mod.ext[:variables][:r_d] = @variable(mod, [jd=JD,jy=JY], upper_bound=0, base_name="REC_d") # added for consistency with H2/code, not used
    r_h = mod.ext[:variables][:r_h] = @variable(mod, [jh=JH,jd=JD,jy=JY], upper_bound=0, base_name="REC_h") # added for consistency with H2/code, not used

    # Create affine expressions  
    mod.ext[:expressions][:e] = @expression(mod, [jy=JY],
        -sum(W[jd]*cc[jh,jd,jy] for jh in JH, jd in JD)
    )
    mod.ext[:expressions][:gw] = @expression(mod, [jh=JH,jd=JD,jy=JY],
        W[jd]*g[jh,jd,jy]
    )
    mod.ext[:expressions][:tot_cost] = @expression(mod, 
        sum(A[jy]*(1-CAP_SV[jy])*IC[jy]*cap[jy] for jy in JY)  
    )

    # Definition of the objective function - will be updated 
    mod.ext[:objective] = @objective(mod, Min, 0)
    
    # Constraints
    mod.ext[:constraints][:carbon_capture_limit] = @constraint(mod, [jh=JH,jd=JD,jy=JY],
        cc[jh,jd,jy] <= (sum(CAP_LT[y2,jy]*cap[y2] for y2=1:jy) + LEG_CAP[jy]) # [MtCO2]        
    )

    # Investment limits: YoY investment is limited
    mod.ext[:constraints][:cap_limit] = @constraint(mod,
        cap[1] <= DELTA_CAP_MAX*LEG_CAP[1] # [MtCO2/h]        
    )
    mod.ext[:constraints][:cap_limit] = @constraint(mod, [jy=2:JY[end]],
        cap[jy] <= DELTA_CAP_MAX*(sum(CAP_LT[y2,jy]*cap[y2] for y2=1:jy-1) + LEG_CAP[jy-1]) # [MtCO2/h]        
    )

    # Electricity consumption  
    mod.ext[:constraints][:gen_limit_energy_sources] = @constraint(mod, [jh=JH,jd=JD,jy=JY],
        cc[jh,jd,jy]/η_E_CO2 == -g[jh,jd,jy]  # [TWh]
    )

    # Selling emission allowances
    mod.ext[:constraints][:EUA_balance]  = @constraint(mod, [jy=JY], 
        sum(b[y2] for y2=1:jy) >=  -sum(W[jd]*cc[jh,jd,jy] for jh in JH, jd in JD) # [MtonCO2]
    )

    # REC    
    mod.ext[:constraints][:REC_balance_yearly] = @constraint(mod, [jy=JY],
        r_y[jy] == ADD_SF[jy]*sum(W[jd]*g[jh,jd,jy] for jh in JH, jd in JD)
    )     

    return mod
end