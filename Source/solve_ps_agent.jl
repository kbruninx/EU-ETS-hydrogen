function solve_ps_agent!(mod::Model)
    # Extract sets
    JH = mod.ext[:sets][:JH]
    JD = mod.ext[:sets][:JD]
    JY = mod.ext[:sets][:JY]
   
    # Extract parameters
    W = mod.ext[:parameters][:W] # weight of the representative days
    VC = mod.ext[:parameters][:VC] # variable costs, excluding cost of carbon emissions
    IC = mod.ext[:parameters][:IC] # overnight investment costs
    A = mod.ext[:parameters][:A] # discount factors
    CAP_SV = mod.ext[:parameters][:CAP_SV] # salvage value of new capacity
    DELTA_CAP_MAX = mod.ext[:parameters][:DELTA_CAP_MAX] # max YoY change in new capacity
    λ_EUA = mod.ext[:parameters][:λ_EUA] # EUA prices
    b_bar = mod.ext[:parameters][:b_bar] # element in ADMM penalty term related to EUA auctions
    ρ_EUA = mod.ext[:parameters][:ρ_EUA] # rho-value in ADMM related to EUA auctions
    λ_EOM = mod.ext[:parameters][:λ_EOM] # EOM prices
    g_bar = mod.ext[:parameters][:g_bar] # element in ADMM penalty term related to EOM
    ρ_EOM = mod.ext[:parameters][:ρ_EOM] # rho-value in ADMM related to EUA auctions
    λ_REC = mod.ext[:parameters][:λ_REC] # EOM prices
    r_bar = mod.ext[:parameters][:r_bar] # element in ADMM penalty term related to REC auctions
    ρ_REC = mod.ext[:parameters][:ρ_REC] # rho-value in ADMM related to EUA auctions

    # Extract variables and expressions
    cap = mod.ext[:variables][:cap]  
    g = mod.ext[:variables][:g]  
    b = mod.ext[:variables][:b]  
    r = mod.ext[:variables][:r]

    # Update objective
    if  mod.ext[:parameters][:REC] == 1 && mod.ext[:parameters][:ETS] == 0
        if DELTA_CAP_MAX > 0 
            mod.ext[:objective] = @objective(mod, Min,
            + sum(A[jy]*(1-CAP_SV[jy])*IC[jy]*cap[jy] for jy in JY)
            + sum(A[jy]*W[jd]*VC[jy]*g[jh,jd,jy] for jh in JH, jd in JD, jy in JY)
            - sum(A[jy]*W[jd]*λ_EOM[jh,jd,jy]*g[jh,jd,jy] for jh in JH, jd in JD, jy in JY)
            - sum(A[jy]*λ_REC[jy]*r[jy] for jy in JY)
            + sum(ρ_EOM/2*W[jd]*(g[jh,jd,jy] - g_bar[jh,jd,jy])^2 for jh in JH, jd in JD, jy in JY)
            + sum(ρ_REC/2*(r[jy] - r_bar[jy])^2 for jy in JY)
            )
        else
            mod.ext[:objective] = @objective(mod, Min,
            + sum(A[jy]*W[jd]*VC[jy]*g[jh,jd,jy] for jh in JH, jd in JD, jy in JY)
            - sum(A[jy]*W[jd]*λ_EOM[jh,jd,jy]*g[jh,jd,jy] for jh in JH, jd in JD, jy in JY)
            - sum(A[jy]*λ_REC[jy]*r[jy] for jy in JY)
            + sum(ρ_EOM/2*W[jd]*(g[jh,jd,jy] - g_bar[jh,jd,jy])^2 for jh in JH, jd in JD, jy in JY)
            + sum(ρ_REC/2*(r[jy] - r_bar[jy])^2 for jy in JY)
            )
        end
    elseif mod.ext[:parameters][:REC] == 0 && mod.ext[:parameters][:ETS] == 1
        if DELTA_CAP_MAX > 0 
            mod.ext[:objective] = @objective(mod, Min,
            + sum(A[jy]*(1-CAP_SV[jy])*IC[jy]*cap[jy] for jy in JY)
            + sum(A[jy]*W[jd]*VC[jy]*g[jh,jd,jy] for jh in JH, jd in JD, jy in JY)
            - sum(A[jy]*W[jd]*λ_EOM[jh,jd,jy]*g[jh,jd,jy] for jh in JH, jd in JD, jy in JY)
            + sum(A[jy]*λ_EUA[jy]*b[jy] for jy in JY)
            + sum(ρ_EUA/2*(b[jy] - b_bar[jy])^2 for jy in JY)
            + sum(ρ_EOM/2*W[jd]*(g[jh,jd,jy] - g_bar[jh,jd,jy])^2 for jh in JH, jd in JD, jy in JY)
            )
        else
            mod.ext[:objective] = @objective(mod, Min,
            + sum(A[jy]*W[jd]*VC[jy]*g[jh,jd,jy] for jh in JH, jd in JD, jy in JY)
            - sum(A[jy]*W[jd]*λ_EOM[jh,jd,jy]*g[jh,jd,jy] for jh in JH, jd in JD, jy in JY)
            + sum(A[jy]*λ_EUA[jy]*b[jy] for jy in JY)
            + sum(ρ_EUA/2*(b[jy] - b_bar[jy])^2 for jy in JY)
            + sum(ρ_EOM/2*W[jd]*(g[jh,jd,jy] - g_bar[jh,jd,jy])^2 for jh in JH, jd in JD, jy in JY)
            )
        end
    elseif mod.ext[:parameters][:REC] == 0 && mod.ext[:parameters][:ETS] == 0
        if DELTA_CAP_MAX > 0 
            mod.ext[:objective] = @objective(mod, Min,
            + sum(A[jy]*(1-CAP_SV[jy])*IC[jy]*cap[jy] for jy in JY)
            + sum(A[jy]*W[jd]*VC[jy]*g[jh,jd,jy] for jh in JH, jd in JD, jy in JY)
            - sum(A[jy]*W[jd]*λ_EOM[jh,jd,jy]*g[jh,jd,jy] for jh in JH, jd in JD, jy in JY)
            + sum(ρ_EOM/2*W[jd]*(g[jh,jd,jy] - g_bar[jh,jd,jy])^2 for jh in JH, jd in JD, jy in JY)
            )
        else
            mod.ext[:objective] = @objective(mod, Min,
            + sum(A[jy]*W[jd]*VC[jy]*g[jh,jd,jy] for jh in JH, jd in JD, jy in JY)
            - sum(A[jy]*W[jd]*λ_EOM[jh,jd,jy]*g[jh,jd,jy] for jh in JH, jd in JD, jy in JY)
            + sum(ρ_EOM/2*W[jd]*(g[jh,jd,jy] - g_bar[jh,jd,jy])^2 for jh in JH, jd in JD, jy in JY)
            )
        end
    else
        if DELTA_CAP_MAX > 0 
            mod.ext[:objective] = @objective(mod, Min,
            + sum(A[jy]*(1-CAP_SV[jy])*IC[jy]*cap[jy] for jy in JY)
            + sum(A[jy]*W[jd]*VC[jy]*g[jh,jd,jy] for jh in JH, jd in JD, jy in JY)
            - sum(A[jy]*W[jd]*λ_EOM[jh,jd,jy]*g[jh,jd,jy] for jh in JH, jd in JD, jy in JY)
            - sum(A[jy]*λ_REC[jy]*r[jy] for jy in JY)
            + sum(A[jy]*λ_EUA[jy]*b[jy] for jy in JY)
            + sum(ρ_EUA/2*(b[jy] - b_bar[jy])^2 for jy in JY)
            + sum(ρ_EOM/2*W[jd]*(g[jh,jd,jy] - g_bar[jh,jd,jy])^2 for jh in JH, jd in JD, jy in JY)
            + sum(ρ_REC/2*(r[jy] - r_bar[jy])^2 for jy in JY)
            )
        else
            mod.ext[:objective] = @objective(mod, Min,
            + sum(A[jy]*W[jd]*VC[jy]*g[jh,jd,jy] for jh in JH, jd in JD, jy in JY)
            - sum(A[jy]*W[jd]*λ_EOM[jh,jd,jy]*g[jh,jd,jy] for jh in JH, jd in JD, jy in JY)
            - sum(A[jy]*λ_REC[jy]*r[jy] for jy in JY)
            + sum(A[jy]*λ_EUA[jy]*b[jy] for jy in JY)
            + sum(ρ_EUA/2*(b[jy] - b_bar[jy])^2 for jy in JY)
            + sum(ρ_EOM/2*W[jd]*(g[jh,jd,jy] - g_bar[jh,jd,jy])^2 for jh in JH, jd in JD, jy in JY)
            + sum(ρ_REC/2*(r[jy] - r_bar[jy])^2 for jy in JY)
            )
        end
    end

    optimize!(mod);

    return mod
end