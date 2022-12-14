function build_ps_agent!(mod::Model)
    # Extract sets
    JH = mod.ext[:sets][:JH]
    JD = mod.ext[:sets][:JD]
    JM = mod.ext[:sets][:JM]
    JY = mod.ext[:sets][:JY]

    # Extract time series data
    AF = mod.ext[:timeseries][:AF] # avaiability factors

    # Extract parameters
    W = mod.ext[:parameters][:W] # weight of the representative days
    Wm = mod.ext[:parameters][:Wm] # weight of the representative days
    if mod.ext[:parameters][:NG] == 1
        VC  = mod.ext[:parameters][:VC] = mod.ext[:parameters][:λ_NG]/mod.ext[:parameters][:η]
    else 
        VC  = mod.ext[:parameters][:VC]  
    end
    IC = mod.ext[:parameters][:IC] # overnight investment costs
    CI = mod.ext[:parameters][:CI] # carbon intensity
    LEG_CAP = mod.ext[:parameters][:LEG_CAP] # legacy capacity
    CAP_LT = mod.ext[:parameters][:CAP_LT] # lead time on new capacity
    CAP_SV = mod.ext[:parameters][:CAP_SV] # salvage value of new capacity
    DELTA_CAP_MAX = mod.ext[:parameters][:DELTA_CAP_MAX] # max YoY change in new capacity
    A = mod.ext[:parameters][:A] # discount factors
    λ_EUA = mod.ext[:parameters][:λ_EUA] # EUA prices
    b_bar = mod.ext[:parameters][:b_bar] # element in ADMM penalty term related to EUA auctions
    ρ_EUA = mod.ext[:parameters][:ρ_EUA] # rho-value in ADMM related to EUA auctions
    λ_EOM = mod.ext[:parameters][:λ_EOM] # EOM prices
    g_bar = mod.ext[:parameters][:g_bar] # element in ADMM penalty term related to EOM
    ρ_EOM = mod.ext[:parameters][:ρ_EOM] # rho-value in ADMM related to EOM auctions
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

    # Create variables
    cap = mod.ext[:variables][:cap] = @variable(mod, [jy=JY], lower_bound=0, base_name="capacity")
    g = mod.ext[:variables][:g] = @variable(mod, [jh=JH,jd=JD,jy=JY], lower_bound=0, base_name="generation")
    b = mod.ext[:variables][:b] = @variable(mod, [jy=JY], lower_bound=0, base_name="EUA") 
    r_y = mod.ext[:variables][:r_y] = @variable(mod, [jy=JY], lower_bound=0, base_name="REC_y") 
    r_m = mod.ext[:variables][:r_m] = @variable(mod, [jd=JM,jy=JY], lower_bound=0, base_name="REC_m") 
    r_d = mod.ext[:variables][:r_d] = @variable(mod, [jd=JD,jy=JY], lower_bound=0, base_name="REC_d") 
    r_h = mod.ext[:variables][:r_h] = @variable(mod, [jh=JH,jd=JD,jy=JY], lower_bound=0, base_name="REC_h") 

    # Create affine expressions - used in postprocessing, haven't tested impact on computational cost of having these here computed in every iteration
    mod.ext[:expressions][:curt] = @expression(mod, [jh=JH,jd=JD,jy=JY],
        AF[jh,jd]*(sum(CAP_LT[y2,jy]*cap[jy] for y2=1:jy) + LEG_CAP[jy]) - g[jh,jd,jy]
    )
    mod.ext[:expressions][:e] = @expression(mod, [jy=JY],
        sum(W[jd]*CI*g[jh,jd,jy] for jh in JH, jd in JD)
    )
    mod.ext[:expressions][:gw] = @expression(mod, [jh=JH,jd=JD,jy=JY],
        W[jd]*g[jh,jd,jy]
    )
    mod.ext[:expressions][:tot_cost] = @expression(mod, 
        + sum(A[jy]*(1-CAP_SV[jy])*IC[jy]*cap[jy] for jy in JY)
        + sum(A[jy]*W[jd]*VC[jy]*g[jh,jd,jy] for jh in JH, jd in JD, jy in JY)
    )

    # Objective 
    mod.ext[:objective] = @objective(mod, Min,
        + sum(A[jy]*(1-CAP_SV[jy])*IC[jy]*cap[jy] for jy in JY)
        + sum(A[jy]*W[jd]*VC[jy]*g[jh,jd,jy] for jh in JH, jd in JD, jy in JY)
        - sum(A[jy]*W[jd]*λ_EOM[jh,jd,jy]*g[jh,jd,jy] for jh in JH, jd in JD, jy in JY)
        - sum(A[jy]*λ_y_REC[jy]*r_y[jy] for jy in JY)
        - sum(A[jy]*λ_m_REC[jm,jy]*r_m[jm,jy] for jm in JM, jy in JY)
        - sum(A[jy]*W[jd]*λ_d_REC[jd,jy]*r_d[jd,jy] for jd in JD, jy in JY)
        - sum(A[jy]*W[jd]*λ_h_REC[jh,jd,jy]*r_h[jh,jd,jy] for jh in JH, jd in JD, jy in JY)
        + sum(A[jy]*λ_EUA[jy]*b[jy] for jy in JY)
        + sum(ρ_EUA/2*(b[jy] - b_bar[jy])^2 for jy in JY)
        + sum(ρ_EOM/2*W[jd]*(g[jh,jd,jy] - g_bar[jh,jd,jy])^2 for jh in JH, jd in JD, jy in JY)
        + sum(ρ_y_REC/2*(r_y[jy] - r_y_bar[jy])^2 for jy in JY)
        + sum(ρ_m_REC/2*(r_m[jm,jy] - r_m_bar[jm,jy])^2 for jm in JM, jy in JY)
        + sum(ρ_d_REC/2*W[jd]*(r_d[jd,jy] - r_d_bar[jd,jy])^2 for jd in JD, jy in JY)
        + sum(ρ_h_REC/2*W[jd]*(r_h[jh,jd,jy] - r_h_bar[jh,jd,jy])^2 for jh in JH, jd in JD, jy in JY)
    )

    if DELTA_CAP_MAX > 0 
        # Capacity constraint
        mod.ext[:constraints][:cap_limit] = @constraint(mod, [jh=JH,jd=JD,jy=JY],
            g[jh,jd,jy] <=  AF[jh,jd]*(sum(CAP_LT[y2,jy]*cap[y2] for y2=1:jy) + LEG_CAP[jy])/1000  # scaling factor needed to go from GW -> TWh
        )

        # Investment limits: YoY investment is limited
        mod.ext[:constraints][:invest_limit] = @constraint(mod, [jy=JY],
            cap[jy] <= DELTA_CAP_MAX
        ) 

        # Generation of RES from new capacity that participates in REC auction
        if mod.ext[:parameters][:REC] == 1
            mod.ext[:constraints][:REC_balance_yearly] = @constraint(mod, [jy=JY],
                r_y[jy] <= sum(W[jd]*AF[jh,jd]*(sum(CAP_LT[y2,jy]*cap[y2] for y2=1:jy) + LEG_CAP[jy])  for jh in JH, jd in JD)/1000 - sum(W[jd]*r_h[jh,jd,jy] for jh in JH,jd in JD) - sum(W[jd]*r_d[jd,jy] for jd in JD)  # scaling factor needed to go from GW -> TWh
            )
            mod.ext[:constraints][:REC_balance_monthly] = @constraint(mod, [jm=JM,jy=JY],
                r_m[jm,jy] <= sum(Wm[jd,jm]*AF[jh,jd]*sum(CAP_LT[y2,jy]*cap[y2] for y2=1:jy) for jh in JH, jd in JD)/1000
            )
            mod.ext[:constraints][:REC_balance_daily] = @constraint(mod, [jd=JD,jy=JY],
                r_d[jd,jy] <=  sum(AF[jh,jd]*sum(CAP_LT[y2,jy]*cap[y2] for y2=1:jy) for jh in JH)/1000
            )
            mod.ext[:constraints][:REC_balance_hourly] = @constraint(mod, [jh=JH,jd=JD,jy=JY],
                r_h[jh,jd,jy] <=  AF[jh,jd]*sum(CAP_LT[y2,jy]*cap[y2] for y2=1:jy)/1000 
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

    else
        # Capacity constraint
        mod.ext[:constraints][:cap_limit] = @constraint(mod, [jh=JH,jd=JD,jy=JY],
            g[jh,jd,jy] <=  AF[jh,jd]*LEG_CAP[jy]/1000  # scaling factor needed to go from GW -> TWh
        )

        # Investment limits: YoY investment is limited
        mod.ext[:constraints][:invest_limit] = @constraint(mod, [jy=JY],
            cap[jy] == 0
        ) 

        # Generation of RES from existing capacity that participates in REC auction - only for annual REC
        if mod.ext[:parameters][:REC] == 1
            mod.ext[:constraints][:REC_balance_yearly] = @constraint(mod, [jy=JY],
                r_y[jy] <= sum(W[jd]*AF[jh,jd]*LEG_CAP[jy] for jh in JH, jd in JD)/1000 
            )
        else
            mod.ext[:constraints][:REC_balance_yearly] = @constraint(mod, [jy=JY],
                r_y[jy] == 0
            )
        end
        mod.ext[:constraints][:REC_balance_monthtly] = @constraint(mod, [jm=JM,jy=JY],
            r_m[jm,jy] == 0 
        )
        mod.ext[:constraints][:REC_balance_daily] = @constraint(mod, [jd=JD,jy=JY],
            r_d[jd,jy] == 0 
        )
        mod.ext[:constraints][:REC_balance_hourly] = @constraint(mod, [jh=JH,jd=JD,jy=JY],
            r_h[jh,jd,jy] == 0 
        )
    end

    # EUA balance 
    if mod.ext[:parameters][:ETS] == 1 
        mod.ext[:constraints][:EUA_balance]  = @constraint(mod,[jy=JY], 
            sum(b[y2] for y2=1:jy) >=  sum(W[jd]*CI*g[jh,jd,y2] for jh in JH, jd in JD, y2=1:jy)
        )
    else
        mod.ext[:constraints][:EUA_balance]  = @constraint(mod,[jy=JY], 
            b[jy] == 0 
        )   
    end

    return mod
end