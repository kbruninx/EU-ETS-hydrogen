function solve_h2import_agent!(mod::Model)
    # Extract sets
    JH = mod.ext[:sets][:JH]
    JD = mod.ext[:sets][:JD]
    JM = mod.ext[:sets][:JM]
    JY = mod.ext[:sets][:JY]
       
    # Extract parameters 
    W = mod.ext[:parameters][:W] # weight of the representative days
    Wm = mod.ext[:parameters][:Wm] # weight of the representative days
    A = mod.ext[:parameters][:A]
    λ_y_H2 = mod.ext[:parameters][:λ_y_H2] # H2 prices
    gH_y_bar = mod.ext[:parameters][:gH_y_bar] # element in ADMM penalty term related to hydrogen market

    α_1 = mod.ext[:parameters][:α_1]
    α_2 = mod.ext[:parameters][:α_2]

    
    # Extract variables and expressions    
    ρ_h_H2 = mod.ext[:parameters][:ρ_h_H2] # rho-value in ADMM related to H2 market
    ρ_d_H2 = mod.ext[:parameters][:ρ_d_H2] # rho-value in ADMM related to H2 market
    ρ_m_H2 = mod.ext[:parameters][:ρ_m_H2] # rho-value in ADMM related to H2 market
    ρ_y_H2 = mod.ext[:parameters][:ρ_y_H2] # rho-value in ADMM related to H2 market
    gH = mod.ext[:variables][:gH]
    gH_m = mod.ext[:variables][:gH_m]
    gH_y = mod.ext[:expressions][:gH_y]

    if ρ_h_H2 > 0
            λ_h_H2 = mod.ext[:parameters][:λ_h_H2] # H2 prices
            gH_h_bar = mod.ext[:parameters][:gH_h_bar] # element in ADMM penalty term related to hydrogen market
            gH = mod.ext[:variables][:gH]
    
            H2_obj = mod.ext[:expressions][:H2_obj] = @expression(mod,
                - sum(A[jy]*W[jd]*λ_h_H2[jh,jd,jy]*gH[jh,jd,jy] for jh in JH, jd in JD, jy in JY)
                + sum(ρ_h_H2/2*W[jd]*(gH[jh,jd,jy] - gH_h_bar[jh,jd,jy])^2 for jh in JH, jd in JD, jy in JY) 
            )
    elseif ρ_d_H2 > 0 
        λ_d_H2 = mod.ext[:parameters][:λ_d_H2] # H2 prices
        gH_d_bar = mod.ext[:parameters][:gH_d_bar] # element in ADMM penalty term related to hydrogen market
        gH_d = mod.ext[:expressions][:gH_d]

        H2_obj = mod.ext[:expressions][:H2_obj] = @expression(mod,
            - sum(A[jy]*W[jd]*λ_d_H2[jd,jy]*gH_d[jd,jy] for jd in JD, jy in JY)
            + sum(ρ_d_H2/2*W[jd]*(gH_d[jd,jy] - gH_d_bar[jd,jy])^2 for jd in JD, jy in JY) 
        )
    elseif ρ_m_H2 > 0 
        λ_m_H2 = mod.ext[:parameters][:λ_m_H2] # H2 prices
        gH_m_bar = mod.ext[:parameters][:gH_m_bar] # element in ADMM penalty term related to hydrogen market
        gH_m = mod.ext[:variables][:gH_m]

        H2_obj = mod.ext[:expressions][:H2_obj] = @expression(mod,
            - sum(A[jy]*λ_m_H2[jm,jy]*gH_m[jm,jy] for jm in JM, jy in JY)
            + sum(ρ_m_H2/2*(gH_m[jm,jy] - gH_m_bar[jm,jy])^2 for jm in JM, jy in JY) 
        )
    elseif ρ_y_H2 > 0 
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
         H2CN_obj = mod.ext[:expressions][:H2CN_obj] = @expression(mod, 0)
    end
    # Update objective 
    mod.ext[:objective] = @objective(mod, Min,
         + sum(A[jy]*(α_2*gH[jh,jd,jy]+ α_1)*gH[jh,jd,jy]  for jh in JH, jd in JD, jy in JY)
         + H2_obj
         + H2CN_obj
    )
     
    optimize!(mod);
 
    return mod
 end