function solve_h2s_agent!(mod::Model)
   # Extract sets
   JH = mod.ext[:sets][:JH]
   JD = mod.ext[:sets][:JD]
   JM = mod.ext[:sets][:JM]
   JY = mod.ext[:sets][:JY]
      
   # Extract parameters
   W = mod.ext[:parameters][:W] # weight of the representative days
   Wm = mod.ext[:parameters][:Wm] # weight of the representative days
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

    # Extract variables and expressions
    capH = mod.ext[:variables][:capH] 
    capHCN = mod.ext[:variables][:capHCN] 
    gH = mod.ext[:variables][:gH]
    gHCN = mod.ext[:variables][:gHCN] 
    g = mod.ext[:variables][:g] 
    dNG = mod.ext[:variables][:dNG] 
    b = mod.ext[:variables][:b]    
    r_y = mod.ext[:variables][:r_y] 
    r_m = mod.ext[:variables][:r_m] 
    r_d = mod.ext[:variables][:r_d]  
    r_h = mod.ext[:variables][:r_h] 
 
    # Update objective
    if  mod.ext[:parameters][:H2CN_prod] == 1 && mod.ext[:parameters][:ETS] == 0 &&  mod.ext[:parameters][:EOM] == 1 # e.g. electrolysis 

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
        + sum(ρ_EOM/2*W[jd]*(g[jh,jd,jy] - g_bar[jh,jd,jy])^2 for jh in JH, jd in JD, jy in JY) # g is electricity
        + sum(ρ_H2/2*(gH[jy] - gH_bar[jy])^2 for jy in JY) 
        + sum(ρ_H2CN_prod/2*(gHCN[jy] - gHCN_bar[jy])^2 for jy in JY)  
        + sum(ρ_H2CN_cap/2*(capHCN[jy] - capHCN_bar[jy])^2 for jy in JY)  
        + sum(ρ_y_REC/2*ADD_SF[jy]*(r_y[jy] - r_y_bar[jy])^2 for jy in JY)
        + sum(ρ_m_REC/2*(r_m[jm,jy] - r_m_bar[jm,jy])^2 for jm in JM, jy in JY)
        + sum(ρ_d_REC/2*W[jd]*(r_d[jd,jy] - r_d_bar[jd,jy])^2 for jd in JD, jy in JY)
        + sum(ρ_h_REC/2*W[jd]*(r_h[jh,jd,jy] - r_h_bar[jh,jd,jy])^2 for jh in JH, jd in JD, jy in JY)
        )

    elseif  mod.ext[:parameters][:H2CN_prod] == 1 && mod.ext[:parameters][:ETS] == 0 &&  mod.ext[:parameters][:EOM] == 0 # e.g. SMR + CCS

        mod.ext[:objective] = @objective(mod, Min,
        + sum(A[jy]*(1-CAP_SV[jy])*IC[jy]*capH[jy] for jy in JY) # [MEUR]
        + sum(A[jy]*λ_NG[jy]*dNG[jy] for jy in JY) 
        - sum(A[jy]*λ_H2[jy]*gH[jy] for jy in JY)
        - sum(A[jy]*λ_H2CN_prod[jy]*gHCN[jy] for jy in JY) 
        - sum(A[jy]*(1-CAP_SV[jy])*λ_H2CN_cap[jy]*capHCN[jy] for jy in JY) 
        + sum(ρ_H2/2*(gH[jy] - gH_bar[jy])^2 for jy in JY) 
        + sum(ρ_H2CN_prod/2*(gHCN[jy] - gHCN_bar[jy])^2 for jy in JY)  
        + sum(ρ_H2CN_cap/2*(capHCN[jy] - capHCN_bar[jy])^2 for jy in JY)  
        )

    elseif mod.ext[:parameters][:H2CN_prod] == 0 && mod.ext[:parameters][:ETS] == 1 && mod.ext[:parameters][:EOM] == 0  # e.g. SMR 

        mod.ext[:objective] = @objective(mod, Min,
        + sum(A[jy]*(1-CAP_SV[jy])*IC[jy]*capH[jy] for jy in JY)  
        + sum(A[jy]*λ_NG[jy]*dNG[jy] for jy in JY) 
        - sum(A[jy]*λ_H2[jy]*gH[jy] for jy in JY)
        + sum(A[jy]*λ_EUA[jy]*b[jy] for jy in JY) 
        + sum(ρ_EUA/2*(b[jy] - b_bar[jy])^2 for jy in JY)
        + sum(ρ_H2/2*(gH[jy] - gH_bar[jy])^2 for jy in JY) 
        )

    else # any other tech

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

    end
    
    optimize!(mod);

    return mod

end