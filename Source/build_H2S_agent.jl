function build_h2s_agent!(mod::Model)
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

    # Create affine expressions - used in postprocessing
    mod.ext[:expressions][:e] = @expression(mod, [jy=JY],
        CI*dNG[jy]
    )
    mod.ext[:expressions][:gw] = @expression(mod, [jh=JH,jd=JD,jy=JY],
        W[jd]*g[jh,jd,jy]
    )
    mod.ext[:expressions][:tot_cost] = @expression(mod, 
        + sum(A[jy]*(1-CAP_SV[jy])*IC[jy]*capH[jy] for jy in JY)  
        + sum(A[jy]*λ_NG[jy]*dNG[jy] for jy in JY) 
    )

    # Definition of the objective function
    mod.ext[:objective] = @objective(mod, Min,
    + sum(A[jy]*(1-CAP_SV[jy])*IC[jy]*capH[jy] for jy in JY) # [MEUR]
    - sum(A[jy]*W[jd]*(λ_EOM[jh,jd,jy])*g[jh,jd,jy] for jh in JH, jd in JD, jy in JY) # [MEUR]
    - sum(A[jy]*ADD_SF[jy]*λ_y_REC[jy]*r_y[jy] for jy in JY)
    - sum(A[jy]*λ_m_REC[jm,jy]*r_m[jm,jy] for jm in JM, jy in JY)
    - sum(A[jy]*W[jd]*λ_d_REC[jd,jy]*r_d[jd,jy] for jd in JD, jy in JY)
    - sum(A[jy]*W[jd]*λ_h_REC[jh,jd,jy]*r_h[jh,jd,jy] for jh in JH, jd in JD, jy in JY)
    + sum(A[jy]*λ_NG[jy]*dNG[jy] for jy in JY) 
    - sum(A[jy]*λ_H2[jy]*gH[jy] for jy in JY)
    - sum(A[jy]*λ_H2CN_prod[jy]*gHCN[jy] for jy in JY) 
    - sum(A[jy]*(1-CAP_SV[jy])*λ_H2CN_cap[jy]*capHCN[jy] for jy in JY) 
    + sum(A[jy]*λ_EUA[jy]*b[jy] for jy in JY) 
    + sum(ρ_EUA/2*(b[jy] - b_bar[jy])^2 for jy in JY)
    + sum(ρ_EOM/2*W[jd]*(g[jh,jd,jy] - g_bar[jh,jd,jy])^2 for jh in JH, jd in JD, jy in JY) # g is electricity
    + sum(ρ_H2/2*(gH[jy] - gH_bar[jy])^2 for jy in JY) 
    + sum(ρ_H2CN_prod/2*(gHCN[jy] - gHCN_bar[jy])^2 for jy in JY)  
    + sum(ρ_H2CN_cap/2*(capHCN[jy] - capHCN_bar[jy])^2 for jy in JY)  
    + sum(ρ_y_REC/2*ADD_SF[jy]*(r_y[jy] - r_y_bar[jy])^2 for jy in JY)
    + sum(ρ_m_REC/2*(r_m[jm,jy] - r_m_bar[jm,jy])^2 for jm in JM, jy in JY)
    + sum(ρ_d_REC/2*W[jd]*(r_d[jd,jy] - r_d_bar[jd,jy])^2 for jd in JD, jy in JY)
    + sum(ρ_h_REC/2*W[jd]*(r_h[jh,jd,jy] - r_h_bar[jh,jd,jy])^2 for jh in JH, jd in JD, jy in JY)
    )
    
    # Constraints
    mod.ext[:constraints][:gen_limit_capacity] = @constraint(mod, [jy=JY],
        gH[jy] <=  8760*((sum(CAP_LT[y2,jy]*capH[y2] for y2=1:jy) + LEG_CAP[jy])/1000) # [TWh]        
    )

    # Investment limits: YoY investment is limited
    mod.ext[:constraints][:cap_limit] = @constraint(mod, [jy=JY],
        capH[jy] <= DELTA_CAP_MAX # [GW]
    )
        
    # Electricity consumption
    mod.ext[:constraints][:elec_consumption] = @constraint(mod, [jh=JH,jd=JD,jy=JY],
        -η_E_H2*g[jh,jd,jy] <= (sum(CAP_LT[y2,jy]*capH[y2] for y2=1:jy) + LEG_CAP[jy])/1000  # [TWh]
    )    

    # Total H2 production
    mod.ext[:constraints][:gen_limit_energy_sources] = @constraint(mod, [jy=JY],
        gH[jy] <= -sum(W[jd]*η_E_H2*g[jh,jd,jy] for jh in JH, jd in JD) + (η_NG_H2*dNG[jy]) # [TWh]
    )
    
    if mod.ext[:parameters][:H2CN_prod] == 1
        mod.ext[:constraints][:gen_limit_carbon_neutral] = @constraint(mod, [jy=JY],
            gHCN[jy] <= gH[jy] # [TWh]
        )
    else
        mod.ext[:constraints][:gen_limit_carbon_neutral] = @constraint(mod, [jy=JY],
            gHCN[jy] == 0 # [TWh]
        )
    end

    if mod.ext[:parameters][:H2CN_cap] == 1
        mod.ext[:constraints][:cap_limit_carbon_neutral]  = @constraint(mod, [jy=JY], 
            capHCN[jy] <= sum(CAP_LT[y2,jy]*capH[y2] for y2=1:jy) + LEG_CAP[jy] # [GW] - cap are capacity additions per year, whereas capHCN needs to be the avaiable capacity
        )
    else
        mod.ext[:constraints][:cap_limit_carbon_neutral]  = @constraint(mod, [jy=JY], 
            capHCN[jy] == 0 # [GW]
        )
    end

    if CI == 0 # Mton CO2/TWh 
        mod.ext[:constraints][:EUA_balance]  = @constraint(mod, [jy=JY], 
            b[jy] == 0 # [MtonCO2]
        )
    else
        mod.ext[:constraints][:EUA_balance]  = @constraint(mod, [jy=JY], 
            sum(b[y2] for y2=1:jy) >=  sum(CI*dNG[y2] for y2=1:jy) # [MtonCO2]
        )
    end

    # recall that g and r are negative numbers (demand for electricity or renewable electricity)
    # ρ_m_REC, ρ_d_REC, ρ_h_REC are only non-zero if additionality is enforced on those time scales
    # in those cases, additionality should not be enforced on an annual basis
    if mod.ext[:parameters][:REC] == 1 && ρ_m_REC == 0 && ρ_d_REC == 0 && ρ_h_REC == 0 
        mod.ext[:constraints][:REC_balance_yearly] = @constraint(mod, [jy=JY],
            -η_E_H2*r_y[jy] >= ADD_SF[jy]*gHCN[jy] 
        )
        mod.ext[:constraints][:REC_balance_yearly_2] = @constraint(mod, [jy=JY],
            r_y[jy] >= ADD_SF[jy]*sum(W[jd]*g[jh,jd,jy] for jh in JH, jd in JD)
        )   
        mod.ext[:constraints][:REC_balance_monthly] = @constraint(mod, [jm=JM,jy=JY],
            r_m[jm,jy] == 0
        )
        mod.ext[:constraints][:REC_balance_daily] = @constraint(mod, [jd=JD,jy=JY],
            r_d[jd,jy] == 0 
        )
        mod.ext[:constraints][:REC_balance_hourly] = @constraint(mod, [jh=JH,jd=JD,jy=JY],
            r_h[jh,jd,jy] == 0 
        )
    elseif mod.ext[:parameters][:REC] == 1 && ρ_m_REC > 0 
        mod.ext[:constraints][:REC_balance_monthly] = @constraint(mod, [jy=JY],
            -η_E_H2*sum(r_m[jm,jy] for jm in JM) >= gHCN[jy] 
        )
        mod.ext[:constraints][:REC_balance_monthly_2] = @constraint(mod, [jm=JM,jy=JY],
            r_m[jm,jy] >= sum(Wm[jd,jm]*g[jh,jd,jy] for jh in JH,jd in JD)
        )
        mod.ext[:constraints][:REC_balance_yearly] = @constraint(mod, [jy=JY],
            r_y[jy] == 0
        )
        mod.ext[:constraints][:REC_balance_daily] = @constraint(mod, [jd=JD,jy=JY],
            r_d[jd,jy] == 0 
        )
        mod.ext[:constraints][:REC_balance_hourly] = @constraint(mod, [jh=JH,jd=JD,jy=JY],
            r_h[jh,jd,jy] == 0 
        )
    elseif mod.ext[:parameters][:REC] == 1 && ρ_d_REC > 0
        mod.ext[:constraints][:REC_balance_daily] = @constraint(mod, [jy=JY],
            -η_E_H2*sum(W[jd]*r_d[jd,jy] for jd in JD)  >= gHCN[jy] 
        )
        mod.ext[:constraints][:REC_balance_daily_2] = @constraint(mod, [jd=JD,jy=JY],
            r_d[jd,jy] >= sum(g[jh,jd,jy] for jh in JH)
        )
        mod.ext[:constraints][:REC_balance_yearly] = @constraint(mod, [jy=JY],
            r_y[jy] == 0
        )
        mod.ext[:constraints][:REC_balance_monthly] = @constraint(mod, [jm=JM,jy=JY],
            r_m[jm,jy] == 0
        )
        mod.ext[:constraints][:REC_balance_hourly] = @constraint(mod, [jh=JH,jd=JD,jy=JY],
            r_h[jh,jd,jy] == 0 
        )
    elseif mod.ext[:parameters][:REC] == 1 && ρ_h_REC > 0
        mod.ext[:constraints][:REC_balance_hourly] = @constraint(mod, [jy=JY],
            -η_E_H2*sum(W[jd]*r_h[jh,jd,jy] for jh in JH, jd in JD) >= gHCN[jy]   
        )
        mod.ext[:constraints][:REC_balance_hourly_2] = @constraint(mod, [jh=JH,jd=JD,jy=JY],
            r_h[jh,jd,jy] >= g[jh,jd,jy] 
        )
        mod.ext[:constraints][:REC_balance_yearly] = @constraint(mod, [jy=JY],
            r_y[jy] == 0
        )
        mod.ext[:constraints][:REC_balance_monthly] = @constraint(mod, [jm=JM,jy=JY],
            r_m[jm,jy] == 0
        )
        mod.ext[:constraints][:REC_balance_daily] = @constraint(mod, [jd=JD,jy=JY],
            r_d[jd,jy] == 0 
        )
    else
        mod.ext[:constraints][:REC_balance_yearly] = @constraint(mod, [jy=JY],
            r_y[jy] == 0
        )
        mod.ext[:constraints][:REC_balance_monthly] = @constraint(mod, [jm=JM,jy=JY],
            r_m[jm,jy] == 0
        )
        mod.ext[:constraints][:REC_balance_daily] = @constraint(mod, [jd=JD,jy=JY],
            r_d[jd,jy] == 0 
        )
        mod.ext[:constraints][:REC_balance_hourly] = @constraint(mod, [jh=JH,jd=JD,jy=JY],
            r_h[jh,jd,jy] == 0 
        )
    end

    return mod

end