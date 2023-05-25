function solve_cc_agent!(mod::Model)
   # Extract sets
   JH = mod.ext[:sets][:JH]
   JD = mod.ext[:sets][:JD]
   JY = mod.ext[:sets][:JY]
      
   # Extract parameters 
   W = mod.ext[:parameters][:W] # weight of the representative days
   IC = mod.ext[:parameters][:IC] # overnight investment costs
   CAP_SV = mod.ext[:parameters][:CAP_SV] # salvage value of new capacity
   A = mod.ext[:parameters][:A] # discount factors
   
   # Extract variables and expressions
   cap = mod.ext[:variables][:cap]  
   cc = mod.ext[:variables][:cc]  

    # Expressions to compute objective
    if mod.ext[:parameters][:EOM] == 1
        λ_EOM = mod.ext[:parameters][:λ_EOM] # EOM prices
        g_bar = mod.ext[:parameters][:g_bar] # element in ADMM penalty term related to EOM
        ρ_EOM = mod.ext[:parameters][:ρ_EOM] # rho-value in ADMM related to EUA auctions
        g = mod.ext[:variables][:g] 

        EOM_obj = mod.ext[:expressions][:EOM_obj] = @expression(mod,
            - sum(A[jy]*W[jd]*(λ_EOM[jh,jd,jy])*g[jh,jd,jy] for jh in JH, jd in JD, jy in JY) # [MEUR]
            + sum(ρ_EOM/2*W[jd]*(g[jh,jd,jy] - g_bar[jh,jd,jy])^2 for jh in JH, jd in JD, jy in JY) 
        )
    else 
        EOM_obj = mod.ext[:expressions][:EOM_obj] = @expression(mod,
            0
        )
    end

    if mod.ext[:parameters][:REC] == 1
        λ_y_REC = mod.ext[:parameters][:λ_y_REC] # REC prices
        r_y_bar = mod.ext[:parameters][:r_y_bar] # element in ADMM penalty term related to REC auctions
        ρ_y_REC = mod.ext[:parameters][:ρ_y_REC] # rho-value in ADMM related to REC auctions
        r_y = mod.ext[:variables][:r_y] 

        REC_obj = mod.ext[:expressions][:REC_obj] = @expression(mod,
            - sum(A[jy]*λ_y_REC[jy]*r_y[jy] for jy in JY)
            + sum(ρ_y_REC/2*(r_y[jy] - r_y_bar[jy])^2 for jy in JY)
        )
    else
        REC_obj = mod.ext[:expressions][:REC_obj] = @expression(mod,
            0
        )
    end 

    if mod.ext[:parameters][:ETS] == 1
        λ_EUA = mod.ext[:parameters][:λ_EUA] # EUA prices
        b_bar = mod.ext[:parameters][:b_bar] # element in ADMM penalty term related to EUA auctions
        ρ_EUA = mod.ext[:parameters][:ρ_EUA] # rho-value in ADMM related to EUA auctions
        b = mod.ext[:variables][:b]    

        ETS_obj = mod.ext[:expressions][:ETS_obj] = @expression(mod,
            + sum(A[jy]*λ_EUA[jy]*b[jy] for jy in JY) 
            + sum(ρ_EUA/2*(b[jy] - b_bar[jy])^2 for jy in JY)
        )
    else
        ETS_obj = mod.ext[:expressions][:ETS_obj] = @expression(mod,
            0
        )
    end
   
    CAPEX_obj = mod.ext[:expressions][:CAPEX_obj] = @expression(mod,
        + sum(A[jy]*(1-CAP_SV[jy])*IC[jy]*cap[jy] for jy in JY) # [MEUR]
    )
       
    OPEX_obj = mod.ext[:expressions][:OPEX_obj] = @expression(mod,
        + sum(A[jy]*(10*10^6*sum(W[jd]*cc[jh,jd,jy] for jh in JH, jd in JD) + 0.001*sum(W[jd]*cc[jh,jd,jy] for jh in JH, jd in JD)^2) for jy in JY) # [MEUR]
    )

    # Update objective 
    mod.ext[:objective] = @objective(mod, Min,
        + CAPEX_obj
        + OPEX_obj
        + EOM_obj
        + REC_obj
        + ETS_obj
    )
    
    optimize!(mod);

    return mod

end