function solve_h2import_agent!(mod::Model)
    # Extract sets
    JH = mod.ext[:sets][:JH]
    JD = mod.ext[:sets][:JD]
    JY = mod.ext[:sets][:JY]
       
    # Extract parameters 
    W = mod.ext[:parameters][:W] # weight of the representative days
    Wm = mod.ext[:parameters][:Wm] # weight of the representative days
    A = mod.ext[:parameters][:A]
    λ_y_H2 = mod.ext[:parameters][:λ_y_H2] # H2 prices
    gH_y_bar = mod.ext[:parameters][:gH_y_bar] # element in ADMM penalty term related to hydrogen market
    ρ_y_H2 = mod.ext[:parameters][:ρ_y_H2] # rho-value in ADMM related to H2 market
    α_H2_import = mod.ext[:parameters][:α_H2_import]

    
    # Extract variables and expressions
    ρ_y_H2 = mod.ext[:parameters][:ρ_y_H2] # rho-value in ADMM related to H2 market
    gH = mod.ext[:variables][:gH]
    gH_m = mod.ext[:variables][:gH_m]
    gH_y = mod.ext[:expressions][:gH_y]

    if ρ_y_H2 > 0 
        λ_y_H2 = mod.ext[:parameters][:λ_y_H2] # H2 prices
        gH_y_bar = mod.ext[:parameters][:gH_y_bar] # element in ADMM penalty term related to hydrogen market
        gH_y = mod.ext[:expressions][:gH_y]

        H2_obj = mod.ext[:expressions][:H2_obj] = @expression(mod,
            - sum(A[jy]*λ_y_H2[jy]*gH_y[jy] for jy in JY)
            + sum(ρ_y_H2/2*(gH_y[jy] - gH_y_bar[jy])^2 for jy in JY) 
        )
    end  
 
    if mod.ext[:parameters][:H2CN_prod] == 1
        λ_H2CN_prod = mod.ext[:parameters][:λ_H2CN_prod] # Carbon neutral H2 generation subsidy
        gHCN_bar = mod.ext[:parameters][:gHCN_bar] # element in ADMM penalty term related to carbon neutral hydrogen generation subsidy
        ρ_H2CN_prod = mod.ext[:parameters][:ρ_H2CN_prod] # rho-value in ADMM related to carbon neutral H2 generation subsidy 
        gHCN = mod.ext[:variables][:gHCN] 

        H2CN_obj = mod.ext[:expressions][:H2CN_obj] = @expression(mod,
            - sum(A[jy]*λ_H2CN_prod[jy]*gHCN[jy] for jy in JY) 
            + sum(ρ_H2CN_prod/2*(gHCN[jy] - gHCN_bar[jy])^2 for jy in JY)  
        )
    else
         H2CN_obj = mod.ext[:expressions][:H2CN_obj] = @expression(mod,
            0
         )
    end
    # Update objective 
    mod.ext[:objective] = @objective(mod, Min,
         + α_H2_import*sum(A[jy]*gH[jh,jd,jy]^2 for jh in JH, jd in JD, jy in JY)
         + H2_obj
         + H2CN_obj
    )
     
    optimize!(mod);
 
    return mod
 end