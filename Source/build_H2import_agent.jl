function build_H2import_agent!(mod)
    # Extract sets
    JY = mod.ext[:sets][:JY]
    
    # Extract parameters
    A = mod.ext[:parameters][:A]
    e = mod.ext[:parameters][:e]
    AC = mod.ext[:parameters][:AC]
    λ_EUA = mod.ext[:parameters][:λ_EUA]
    b_bar = mod.ext[:parameters][:b_bar]
    ρ_EUA = mod.ext[:parameters][:ρ_EUA]

    # Extract sets
    JH = mod.ext[:sets][:JH]
    JD = mod.ext[:sets][:JD]
    JM = mod.ext[:sets][:JM]
    JY = mod.ext[:sets][:JY]
       
    # Extract parameters
    W = mod.ext[:parameters][:W] # weight of the representative days
    Wm = mod.ext[:parameters][:Wm] # weight of the representative days
    IC = mod.ext[:parameters][:IC] # overnight investment costs
    CI = mod.ext[:parameters][:CI] # carbon intensity
    LEG_CAP = mod.ext[:parameters][:LEG_CAP] # legacy capacity
    CAP_LT = mod.ext[:parameters][:CAP_LT] # lead time on new capacity
    CAP_SV = mod.ext[:parameters][:CAP_SV] # salvage value of new capacity
    DELTA_CAP_MAX = mod.ext[:parameters][:DELTA_CAP_MAX] # max YoY change in new capacity
    A = mod.ext[:parameters][:A] # discount factors
    η_E_H2 = mod.ext[:parameters][:η_E_H2] # efficiency E->H2
    η_NG_H2 = mod.ext[:parameters][:η_NG_H2] # efficiency NG->H2
    λ_NG = mod.ext[:parameters][:λ_NG] # natural gas price

    # ADMM algorithm parameters
    λ_EUA = mod.ext[:parameters][:λ_EUA] # EUA prices
    b_bar = mod.ext[:parameters][:b_bar] # element in ADMM penalty term related to EUA auctions
    ρ_EUA = mod.ext[:parameters][:ρ_EUA] # rho-value in ADMM related to EUA auctions
    λ_EOM = mod.ext[:parameters][:λ_EOM] # EOM prices
    g_bar = mod.ext[:parameters][:g_bar] # element in ADMM penalty term related to EOM
    ρ_EOM = mod.ext[:parameters][:ρ_EOM] # rho-value in ADMM related to EUA auctions
    λ_H2 = mod.ext[:parameters][:λ_H2] # H2 prices
    gH_bar = mod.ext[:parameters][:gH_bar] # element in ADMM penalty term related to hydrogen market
    ρ_H2 = mod.ext[:parameters][:ρ_H2] # rho-value in ADMM related to H2 market
    λ_H2CN_prod = mod.ext[:parameters][:λ_H2CN_prod] # Carbon neutral H2 generation subsidy
    gHCN_bar = mod.ext[:parameters][:gHCN_bar] # element in ADMM penalty term related to carbon neutral hydrogen generation subsidy
    ρ_H2CN_prod = mod.ext[:parameters][:ρ_H2CN_prod] # rho-value in ADMM related to carbon neutral H2 generation subsidy 
    λ_H2CN_cap = mod.ext[:parameters][:λ_H2CN_cap] # Carbon neutral H2 capacity subsidy
    capHCN_bar = mod.ext[:parameters][:capHCN_bar] # element in ADMM penalty term related to carbon neutral hydrogen capacity subsidy
    ρ_H2CN_cap = mod.ext[:parameters][:ρ_H2CN_cap] # rho-value in ADMM related to carbon neutral H2 capacity subsidy 
    λ_y_REC = mod.ext[:parameters][:λ_y_REC] # REC prices
    r_y_bar = mod.ext[:parameters][:r_y_bar] # element in ADMM penalty term related to REC auctions
    ρ_y_REC = mod.ext[:parameters][:ρ_y_REC] # rho-value in ADMM related to REC auctions
    λ_m_REC = mod.ext[:parameters][:λ_m_REC] # REC prices
    r_m_bar = mod.ext[:parameters][:r_m_bar] # element in ADMM penalty term related to REC auctions
    ρ_m_REC = mod.ext[:parameters][:ρ_m_REC] # rho-value in ADMM related to REC auctions
    λ_d_REC = mod.ext[:parameters][:λ_d_REC] # REC prices
    r_d_bar = mod.ext[:parameters][:r_d_bar] # element in ADMM penalty term related to REC auctions
    ρ_d_REC = mod.ext[:parameters][:ρ_d_REC] # rho-value in ADMM related to REC auctions
    λ_h_REC = mod.ext[:parameters][:λ_h_REC] # REC prices
    r_h_bar = mod.ext[:parameters][:r_h_bar] # element in ADMM penalty term related to REC auctions
    ρ_h_REC = mod.ext[:parameters][:ρ_h_REC] # rho-value in ADMM related to REC auctions
    ADD_SF = mod.ext[:parameters][:ADD_SF] 

    # Decision variables
    capH = mod.ext[:variables][:capH] = @variable(mod, [jy=JY], lower_bound=0, base_name="capacity")
    capHCN = mod.ext[:variables][:capHCN] = @variable(mod, [jy=JY], lower_bound=0, base_name="green_capacity")
    gH = mod.ext[:variables][:gH] = @variable(mod, [jy=JY], lower_bound=0, base_name="generation_hydrogen")
    gHCN = mod.ext[:variables][:gHCN] = @variable(mod, [jy=JY], lower_bound=0, base_name="generation_carbon_neutral_hydrogen")
    g = mod.ext[:variables][:g] = @variable(mod, [jh=JH,jd=JD,jy=JY], upper_bound=0, base_name="demand_electricity_hydrogen") # note this is defined as a negative number, consumption
    dNG = mod.ext[:variables][:dNG] = @variable(mod, [jy=JY], lower_bound=0, base_name="demand_natural_gas_hydrogen")
    b = mod.ext[:variables][:b] = @variable(mod, [jy=JY], lower_bound=0, base_name="EUA") 
    r_y = mod.ext[:variables][:r_y] = @variable(mod, [jy=JY], upper_bound=0, base_name="REC_y") 
    r_m = mod.ext[:variables][:r_m] = @variable(mod, [jd=JM,jy=JY], lower_bound=0, base_name="REC_m") 
    r_d = mod.ext[:variables][:r_d] = @variable(mod, [jd=JD,jy=JY], upper_bound=0, base_name="REC_d") 
    r_h = mod.ext[:variables][:r_h] = @variable(mod, [jh=JH,jd=JD,jy=JY], upper_bound=0, base_name="REC_h") 
    
    # Define variables
    b = mod.ext[:variables][:b] = @variable(mod, [jy=JY], lower_bound = 0, base_name="EUA") 
    
    # Expressions
    mod.ext[:expressions][:tot_cost] = @expression(mod, 
        sum(A[jy]*AC[jy] for jy in JY)
    )
    
    # Objective
    @objective(mod, Min, sum(A[jy]*λ_EUA[jy]*b[jy] for jy in JY) +sum(ρ_EUA/2*(b[jy] - b_bar[jy])^2 for jy in JY))
    
    # Constraints 
    mod.ext[:constraints][:con1]  = @constraint(mod,[jy=JY], 
        sum(b[y2] for y2=1:jy) >= sum(e[y2] for y2=1:jy)
    )
    
    optimize!(mod);
    
    return mod
    end
    