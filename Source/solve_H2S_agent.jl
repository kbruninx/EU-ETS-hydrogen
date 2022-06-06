function solve_h2s_agent!(mod::Model)
   # Extract sets
   JH = mod.ext[:sets][:JH]
   JD = mod.ext[:sets][:JD]
   JY = mod.ext[:sets][:JY]
      
   # Extract parameters
   W = mod.ext[:parameters][:W] # weight of the representative days
   IC = mod.ext[:parameters][:IC] # overnight investment costs
   CAP_SV = mod.ext[:parameters][:CAP_SV] # salvage value of new capacity
   A = mod.ext[:parameters][:A] # discount factors
   λ_NG = mod.ext[:parameters][:λ_NG] # Gas prices
   
   # ADMN algorithm parameters
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

    # Extract variables and expressions
    capH = mod.ext[:variables][:capH] 
    capHCN = mod.ext[:variables][:capHCN] 
    gH = mod.ext[:variables][:gH]
    gHCN = mod.ext[:variables][:gHCN] 
    g = mod.ext[:variables][:g] 
    dNG = mod.ext[:variables][:dNG] 
    b = mod.ext[:variables][:b]    
 
    # Update objective
    if  mod.ext[:parameters][:H2CN_prod] == 1 && mod.ext[:parameters][:ETS] == 0 
        
        mod.ext[:objective] = @objective(mod, Min,
        + sum(A[jy]*(1-CAP_SV[jy])*IC[jy]*capH[jy] for jy in JY) 
        - sum(A[jy]*W[jd]*λ_EOM[jh,jd,jy]*g[jh,jd,jy] for jh in JH, jd in JD, jy in JY)
        + sum(A[jy]*λ_NG[jy]*dNG[jy] for jy in JY) 
        - sum(A[jy]*λ_H2[jy]*gH[jy] for jy in JY) 
        - sum(A[jy]*λ_H2CN_prod[jy]*gHCN[jy] for jy in JY) 
        - sum(A[jy]*(1-CAP_SV[jy])*λ_H2CN_cap[jy]*capHCN[jy] for jy in JY)  
        + sum(ρ_EOM/2*W[jd]*(g[jh,jd,jy] - g_bar[jh,jd,jy])^2 for jh in JH, jd in JD, jy in JY)
        + sum(ρ_H2/2*(gH[jy] - gH_bar[jy])^2 for jy in JY) 
        + sum(ρ_H2CN_prod/2*(gHCN[jy] - gHCN_bar[jy])^2 for jy in JY)
        + sum(ρ_H2CN_cap/2*(capHCN[jy] - capHCN_bar[jy])^2 for jy in JY)
        )            
    
    elseif mod.ext[:parameters][:H2CN_prod] == 0 && mod.ext[:parameters][:ETS] == 1
        mod.ext[:objective] = @objective(mod, Min,
        + sum(A[jy]*(1-CAP_SV[jy])*IC[jy]*capH[jy] for jy in JY)  
        + sum(A[jy]*λ_NG[jy]*dNG[jy] for jy in JY) 
        - sum(A[jy]*λ_H2[jy]*gH[jy] for jy in JY)
        + sum(A[jy]*λ_EUA[jy]*b[jy] for jy in JY) 
        + sum(ρ_EUA/2*(b[jy] - b_bar[jy])^2 for jy in JY)
        + sum(ρ_H2/2*(gH[jy] - gH_bar[jy])^2 for jy in JY) 
        )

    end
    
    optimize!(mod);

    return mod

end