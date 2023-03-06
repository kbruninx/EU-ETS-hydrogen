function build_h2import_agent!(mod)
    # Extract sets
    JH = mod.ext[:sets][:JH]
    JD = mod.ext[:sets][:JD]
    JM = mod.ext[:sets][:JM]
    JY = mod.ext[:sets][:JY]
    
    # Extract parameters
    W = mod.ext[:parameters][:W] # weight of the representative days
    Wm = mod.ext[:parameters][:Wm] # weight of the representative days
    A = mod.ext[:parameters][:A]
    λ_h_H2 = mod.ext[:parameters][:λ_h_H2] # H2 prices
    gH_h_bar = mod.ext[:parameters][:gH_h_bar] # element in ADMM penalty term related to hydrogen market
    ρ_h_H2 = mod.ext[:parameters][:ρ_h_H2] # rho-value in ADMM related to H2 market
    λ_d_H2 = mod.ext[:parameters][:λ_d_H2] # H2 prices
    gH_d_bar = mod.ext[:parameters][:gH_d_bar] # element in ADMM penalty term related to hydrogen market
    ρ_d_H2 = mod.ext[:parameters][:ρ_d_H2] # rho-value in ADMM related to H2 market
    λ_m_H2 = mod.ext[:parameters][:λ_m_H2] # H2 prices
    gH_m_bar = mod.ext[:parameters][:gH_m_bar] # element in ADMM penalty term related to hydrogen market
    ρ_m_H2 = mod.ext[:parameters][:ρ_m_H2] # rho-value in ADMM related to H2 market
    λ_y_H2 = mod.ext[:parameters][:λ_y_H2] # H2 prices
    gH_y_bar = mod.ext[:parameters][:gH_y_bar] # element in ADMM penalty term related to hydrogen market
    ρ_y_H2 = mod.ext[:parameters][:ρ_y_H2] # rho-value in ADMM related to H2 market
    λ_H2CN_prod = mod.ext[:parameters][:λ_H2CN_prod] # Carbon neutral H2 generation subsidy
    gHCN_bar = mod.ext[:parameters][:gHCN_bar] # element in ADMM penalty term related to carbon neutral hydrogen generation subsidy
    ρ_H2CN_prod = mod.ext[:parameters][:ρ_H2CN_prod] # rho-value in ADMM related to carbon neutral H2 generation subsidy 
    λ_H2CN_cap = mod.ext[:parameters][:λ_H2CN_cap] # Carbon neutral H2 capacity subsidy
    capHCN_bar = mod.ext[:parameters][:capHCN_bar] # element in ADMM penalty term related to carbon neutral hydrogen capacity subsidy
    ρ_H2CN_cap = mod.ext[:parameters][:ρ_H2CN_cap] # rho-value in ADMM related to carbon neutral H2 capacity subsidy 
    ADD_SF = mod.ext[:parameters][:ADD_SF] 

    α_H2_import = mod.ext[:parameters][:α_H2_import]

    # Decision variables
    gH = mod.ext[:variables][:gH] = @variable(mod, [jh=JH,jd=JD,jy=JY], lower_bound=0, base_name="generation_hydrogen") 
    gHCN = mod.ext[:variables][:gHCN] = @variable(mod, [jy=JY], lower_bound=0, base_name="generation_carbon_neutral_hydrogen")
    #gH_m = mod.ext[:variables][:gH_m] = @variable(mod, [jm=JM,jy=JY],lower_bound=0, base_name="generation_hydrogen_monthly") # needs to be variable to get feasible solution with representative days (combination of days may not allow exact match of montly demand, may be infeasible)
    
    # Create affine expressions 
    gH_y = mod.ext[:expressions][:gH_y] = @expression(mod, [jy=JY],
    sum(W[jd]*gH[jh,jd,jy] for jh in JH, jd in JD)
    )
    gH_d = mod.ext[:expressions][:gH_d] = @expression(mod, [jd=JD,jy=JY],
    sum(gH[jh,jd,jy] for jh in JH)
    )
    gH_h_w = mod.ext[:expressions][:gH_h_w] = @expression(mod, [jh=JH,jd=JD,jy=JY],
    W[jd]*gH[jh,jd,jy] 
    )
    gH_d_w = mod.ext[:expressions][:gH_d_w] = @expression(mod, [jd=JD,jy=JY],
    W[jd]*sum(gH[jh,jd,jy] for jh in JH)
    )
    gH_m = mod.ext[:variables][:gH_m] = @expression(mod, [jm=JM,jy=JY],
    sum(Wm[jd]*sum(gH[jh,jd,jy] for jh in JH) for jd in JD)
    )
    mod.ext[:expressions][:tot_cost] = @expression(mod, 
    sum(A[jy]*(α_H2_import*gH[jh,jd,jy]+ 53)*gH[jh,jd,jy] for jh in JH, jd in JD, jy in JY)
    )
    
    # Objective
    @objective(mod, Min,
    + sum(A[jy]*(α_H2_import*gH[jh,jd,jy]+ 53)*gH[jh,jd,jy]  for jh in JH, jd in JD, jy in JY)
    - sum(A[jy]*gH_y[jy]*λ_y_H2[jy] for jy in JY)
    - sum(A[jy]*λ_H2CN_prod[jy]*gHCN[jy] for jy in JY) 
    + sum(ρ_y_H2/2*(gH_y[jy] - gH_y_bar[jy])^2 for jy in JY)
    + sum(ρ_H2CN_prod/2*(gHCN[jy] - gHCN_bar[jy])^2 for jy in JY)
    )

    #mod.ext[:constraints][:H2_balance] = @constraint(mod, [jh=JH,jd=JD,jy=JY], 
    #gH[jh,jd,jy] == gH[1,1,jy] 
    #)
    if mod.ext[:parameters][:H2CN_prod] == 1
        mod.ext[:constraints][:gen_limit_carbon_neutral] = @constraint(mod, [jy=JY],
            gHCN[jy] <=  sum(W[jd]*gH[jh,jd,jy] for jh in JH, jd in JD) # [TWh]
        )
    else
        mod.ext[:constraints][:gen_limit_carbon_neutral] = @constraint(mod, [jy=JY],
            gHCN[jy] == 0 # [TWh]
        )
    end
    
    optimize!(mod);
    
    return mod
    end
    