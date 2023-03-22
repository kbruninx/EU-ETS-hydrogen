# Save results
function save_results(mdict::Dict,EOM::Dict,ETS::Dict,H2::Dict,ADMM::Dict,results::Dict,data::Dict,agents::Dict,sens) 
    # note that type of "sens" is not defined as a string stored in a dictionary is of type String31, whereas a "regular" string is of type String. Specifying one or the other may trow errors.
    Years = range(2021,stop=2021+data["nyears"]-1)

    # Aggregate metrics 
    tot_cost = sum(value(mdict[m].ext[:expressions][:tot_cost]) for m in agents[:all])
    tot_em = sum(results["e"][m][end][jy] for m in agents[:ets],jy in mdict[agents[:ps][1]].ext[:sets][:JY]) 
    
    vector_output = [data["scen_number"]; sens; ADMM["n_iter"]; ADMM["walltime"];ADMM["Residuals"]["Primal"]["ETS"][end];ADMM["Residuals"]["Primal"]["MSR"][end]; 
                     ADMM["Residuals"]["Primal"]["EOM"][end];ADMM["Residuals"]["Primal"]["REC_y"][end]+ADMM["Residuals"]["Primal"]["REC_m"][end]+ADMM["Residuals"]["Primal"]["REC_d"][end]+ADMM["Residuals"]["Primal"]["REC_h"][end]; ADMM["Residuals"]["Primal"]["H2_h"][end]+ADMM["Residuals"]["Primal"]["H2_d"][end]+ADMM["Residuals"]["Primal"]["H2_m"][end]+ADMM["Residuals"]["Primal"]["H2_y"][end]; ADMM["Residuals"]["Primal"]["H2CN_prod"][end]; ADMM["Residuals"]["Primal"]["H2CN_cap"][end]; ADMM["Residuals"]["Dual"]["ETS"][end]; ADMM["Residuals"]["Dual"]["EOM"][end]; ADMM["Residuals"]["Dual"]["REC_y"][end]+ADMM["Residuals"]["Dual"]["REC_m"][end]+ADMM["Residuals"]["Dual"]["REC_d"][end]+ADMM["Residuals"]["Dual"]["REC_h"][end]; ADMM["Residuals"]["Dual"]["H2_h"][end]+ADMM["Residuals"]["Dual"]["H2_d"][end]+ADMM["Residuals"]["Dual"]["H2_m"][end]+ADMM["Residuals"]["Dual"]["H2_y"][end];ADMM["Residuals"]["Dual"]["H2CN_prod"][end]; ADMM["Residuals"]["Dual"]["H2CN_cap"][end]; mdict["Ind"].ext[:parameters][:β]; results[ "λ"]["EUA"][end][2]; tot_em; tot_cost]
    CSV.write(joinpath(home_dir,string("overview_results_",data["nReprDays"],"_repr_days.csv")), DataFrame(reshape(vector_output,1,:),:auto), delim=";",append=true);

    # ETS
    # Note: 
    # TNAC will be shifted by 2 years (i.e., TNAC[y] is the TNAC at the end of year y-2)
    # MSR will be shifted by 1 yeare (i.e., MSR[y,12] is the MSR at the end of year y-1)
    mat_output = [Years ETS["CAP"] ETS["S"] sum(ETS["C"][:,:],dims=2) ETS["MSR"][2:end,12] ETS["TNAC"][3:end] sum(results["e"][m][end] for m in agents[:ind]) sum(results["e"][m][end] for m in setdiff(agents[:ets],union(agents[:ind],agents[:ps]))) sum(results["e"][m][end] for m in setdiff(agents[:ets],union(agents[:ind],agents[:h2s]))) results[ "λ"]["EUA"][end] sum(results["b"][m][end] for m in agents[:ind]) sum(results["b"][m][end] for m in setdiff(agents[:ets],agents[:ind]))]
    CSV.write(joinpath(home_dir,string("Results_",data["nReprDays"],"_repr_days"),string("Scenario_",data["scen_number"],"_ETS_",sens,".csv")), DataFrame(mat_output,:auto), delim=";",header=["Year";"CAP";"Supply";"Cancellation";"MSR";"TNAC";"Emissions_Ind";"Emissions_H2S"; "Emissions_PS"; "EUAprice"; "EUAs_Ind"; "EUAs_PS"]);
     
    # Power sector
    fuel_shares = zeros(length(agents[:ps]),data["nyears"])
    available_cap = zeros(length(agents[:ps]),data["nyears"])
    add_cap = zeros(length(agents[:ps]),data["nyears"])
    mm = 1
    for m in agents[:ps]
        gw = value.(mdict[m].ext[:expressions][:gw])
        fuel_shares[mm,:] = sum(gw[jh,jd,:] for jh in mdict[m].ext[:sets][:JH], jd in mdict[m].ext[:sets][:JD])
        CAP_LT = mdict[m].ext[:parameters][:CAP_LT]
        LEG_CAP = mdict[m].ext[:parameters][:LEG_CAP]
        add_cap[mm,:] = cap = value.(mdict[m].ext[:variables][:cap])
        available_cap[mm,:] = [sum(CAP_LT[y2,jy]*cap[y2] for y2=1:jy) + LEG_CAP[jy] for jy in mdict[m].ext[:sets][:JY]]
        mm = mm+1
    end
    gw = Dict()
    for m in agents[:eom] # including electrolysis demand
        gw[m] = value.(mdict[m].ext[:expressions][:gw])
    end
    gw_tot = sum(gw[m] for m in agents[:eom]) # total electricity generation
    gw_res_tot = sum(results["r_y"][m][end][:] for m in agents[:rec]) # total renewable electricty generation participating in the res auctions
    # For the production weighted averages below I assume that REC support for RES from electrolyzers (or electricity costs) are internal transfers - not taken into account:
    λ_EOM_avg = [sum(gw_tot[:,:,jy].*results[ "λ"]["EOM"][end][:,:,jy])./sum(gw_tot[:,:,jy]) for jy in mdict[agents[:ps][1]].ext[:sets][:JY]] # production weighted average electricity price
    λ_REC_avg = [gw_res_tot[jy].*results[ "λ"]["REC_y"][end][jy]./gw_res_tot[jy] for jy in mdict[agents[:ps][1]].ext[:sets][:JY]] # production weighted support for RES  
    mat_output = [Years λ_EOM_avg λ_REC_avg transpose(add_cap) transpose(available_cap) transpose(fuel_shares)]
    CSV.write(joinpath(home_dir,string("Results_",data["nReprDays"],"_repr_days"),string("Scenario_",data["scen_number"],"_PS_",sens,".csv")), DataFrame(mat_output,:auto), delim=";",header=["Year";"EOM_avg";"REC_y";string.("ADD_CAP_",agents[:ps]);string.("CAP_",agents[:ps]);string.("FS_",agents[:ps])]);
    
    # Hydrogen sector 
    h2_cap = zeros(length(agents[:h2s]),data["nyears"])
    h2_prod = zeros(length(agents[:h2s]),data["nyears"])
    h2_add_cap = zeros(length(agents[:h2s]),data["nyears"])

    mm = 1
    for m in agents[:h2s]
        CAP_LT = mdict[m].ext[:parameters][:CAP_LT]
        LEG_CAP = mdict[m].ext[:parameters][:LEG_CAP]
        h2_add_cap[mm,:] = cap = value.(mdict[m].ext[:variables][:capH])
        h2_cap[mm,:] = [sum(CAP_LT[y2,jy]*cap[y2] for y2=1:jy) + LEG_CAP[jy] for jy in mdict[m].ext[:sets][:JY]]    
        h2_prod[mm,:] = value.(mdict[m].ext[:expressions][:gH_y])./data["conv_factor"]
        mm = mm+1
    end
    h2cn_prod = zeros(length(agents[:h2cn_prod]),data["nyears"])
    h2cn_cap = zeros(length(agents[:h2cn_prod]),data["nyears"])
    mm = 1
    for m in agents[:h2cn_prod]
        h2cn_cap[mm,:] = value.(mdict[m].ext[:variables][:capHCN])
        h2cn_prod[mm,:] = value.(mdict[m].ext[:variables][:gHCN])./data["conv_factor"]
        mm = mm+1
    end
    gHw_h = Dict()
    gHw_d = Dict()
    gHw_m = Dict()
    for m in agents[:h2s]  
        gHw_h[m] = value.(mdict[m].ext[:expressions][:gH_h_w])
        gHw_d[m] = value.(mdict[m].ext[:expressions][:gH_d_w])
        gHw_m[m] = value.(mdict[m].ext[:variables][:gH_m])
    end
    gHw_h_tot = sum(gHw_h[m] for m in agents[:h2s]) # total hydrogen production, weighted
    gHw_d_tot = sum(gHw_d[m] for m in agents[:h2s]) # total hydrogen production, weighted
    gHw_m_tot = sum(gHw_m[m] for m in agents[:h2s]) # total hydrogen production 

    if data["H2_balance"] == "Hourly"
        λ_H2_avg = [sum(gHw_h_tot[:,:,jy].*results["λ"]["H2_h"][end][:,:,jy])./sum(gHw_h_tot[:,:,jy])*data["conv_factor"]/1000 for jy in mdict[agents[:h2s][1]].ext[:sets][:JY]]
    elseif data["H2_balance"] == "Daily"
        λ_H2_avg = [sum(gHw_d_tot[:,jy].*results["λ"]["H2_d"][end][:,jy])./sum(gHw_d_tot[:,jy])*data["conv_factor"]/1000 for jy in mdict[agents[:h2s][1]].ext[:sets][:JY]]
    elseif data["H2_balance"] == "Monthly"
        λ_H2_avg = [sum(gHw_m_tot[:,jy].*results["λ"]["H2_m"][end][:,jy])./sum(gHw_m_tot[:,jy])*data["conv_factor"]/1000 for jy in mdict[agents[:h2s][1]].ext[:sets][:JY]]
    elseif data["H2_balance"] == "Yearly"
        λ_H2_avg = results["λ"]["H2_y"][end]*data["conv_factor"]/1000
    end

    mat_output = [Years transpose(h2_add_cap) transpose(h2_cap) transpose(h2_prod) transpose(h2cn_cap) transpose(h2cn_prod) λ_H2_avg results["λ"]["H2CN_prod"][end]*data["conv_factor"]/1000 results["λ"]["H2CN_cap"][end]]
    CSV.write(joinpath(home_dir,string("Results_",data["nReprDays"],"_repr_days"),string("Scenario_",data["scen_number"],"_H2_",sens,".csv")), DataFrame(mat_output,:auto), delim=";",header=["Year";string.("ADD_CAP_",agents[:h2s]);string.("CAP_",agents[:h2s]);string.("PROD_",agents[:h2s]);string.("CN_CAP_",agents[:h2cn_prod]);string.("CN_PROD_",agents[:h2cn_prod]);"PriceH2";"PremiumH2CN_prod";"PremiumH2CN_cap"]);

    # Operational data electricity market
    Timestep = collect(1:data["nTimesteps"]*data["nReprDays"]*data["nyears"])
    el_prod = zeros(length(agents[:eom]), data["nReprDays"] * data["nTimesteps"]*data["nyears"])
    mm = 1
    for m in agents[:eom]
        el_prod[mm, :] = vec(results["g"][m][end][:,:,:])
        mm += 1
    end
    eom_price = vec(results["λ"]["EOM"][end][:,:,:])
    demand = vec(EOM["D"])
    mat_output = [Timestep eom_price -demand transpose(el_prod)]
    CSV.write(
        joinpath(
            home_dir,
            string("Results_", data["nReprDays"], "_repr_days"),
            string("Scenario_", data["scen_number"], "_operational_data_eom_", sens, ".csv")
        ),
        DataFrame(mat_output, :auto),
        delim=";",
        header=[
            "Timestep"; 
            "EOM_price";
            "Demand";
            string.(agents[:eom])
            ]
    )

     # Operational data hydrogen
     h2_prod = zeros(length(agents[:h2s]), data["nReprDays"] * data["nTimesteps"]*data["nyears"])
     mm = 1
     for m in agents[:h2s]
         h2_prod[mm, :] = vec(results["h2_h"][m][end][:,:,:])
         mm += 1
     end
     h2_h_price = vec(results["λ"]["H2_h"][end][:,:,:])
     h2_y_to_h_price = vec([results["λ"]["H2_y"][end][y] for h=1:data["nTimesteps"], d=1:data["nReprDays"], y=1:data["nyears"]])
     h2_d_to_h_price = vec([results["λ"]["H2_d"][end][d,y] for h=1:data["nTimesteps"], d=1:data["nReprDays"], y=1:data["nyears"]])
     h2_m_to_h_price = vec([sum(results["λ"]["H2_m"][end][m,y] for m=1:12)/12 for h=1:data["nTimesteps"], d=1:data["nReprDays"], y=1:data["nyears"]])
     demand = vec(H2["D_h"])
     mat_output = [Timestep h2_h_price h2_d_to_h_price h2_m_to_h_price h2_y_to_h_price demand transpose(h2_prod)]
     CSV.write(
         joinpath(
             home_dir,
             string("Results_", data["nReprDays"], "_repr_days"),
             string("Scenario_", data["scen_number"], "_operational_data_h2m_", sens, ".csv")
         ),
         DataFrame(mat_output, :auto),
         delim=";",
         header=[
             "Timestep"; 
             "H2M_price_h";
             "H2M_price_d";
             "H2M_price_m";
             "H2M_price_y";
             "demand";
             string.(agents[:h2s])
             ]
     )

    # Operational data REC
    rec_prod_h = zeros(length(agents[:rec]), data["nReprDays"] * data["nTimesteps"]*data["nyears"])
    rec_prod_d = zeros(length(agents[:rec]), data["nReprDays"] * data["nTimesteps"]*data["nyears"])
    rec_prod_m = zeros(length(agents[:rec]), data["nReprDays"] * data["nTimesteps"]*data["nyears"])
    rec_prod_y = zeros(length(agents[:rec]), data["nReprDays"] * data["nTimesteps"]*data["nyears"])

    mm = 1
    for m in agents[:rec]
        rec_prod_h[mm, :] = vec(results["r_h"][m][end][:,:,:])
        rec_prod_d[mm,:] = vec([results["r_d"][m][end][d,y] for h=1:data["nTimesteps"], d=1:data["nReprDays"], y=1:data["nyears"]])
        rec_prod_m[mm,:] = vec([sum(results["r_m"][m][end][mm,y] for mm=1:12)/12 for h=1:data["nTimesteps"], d=1:data["nReprDays"], y=1:data["nyears"]])
        rec_prod_y[mm,:] = vec([results["r_y"][m][end][y] for h=1:data["nTimesteps"], d=1:data["nReprDays"], y=1:data["nyears"]])
        mm += 1
    end
    REC_h_price = vec(results["λ"]["REC_h"][end][:,:,:])
    REC_y_to_h_price = vec([results["λ"]["REC_y"][end][y] for h=1:data["nTimesteps"], d=1:data["nReprDays"], y=1:data["nyears"]])
    REC_d_to_h_price = vec([results["λ"]["REC_d"][end][d,y] for h=1:data["nTimesteps"], d=1:data["nReprDays"], y=1:data["nyears"]])
    REC_m_to_h_price = vec([sum(results["λ"]["REC_m"][end][m,y] for m=1:12)/12 for h=1:data["nTimesteps"], d=1:data["nReprDays"], y=1:data["nyears"]])
    
    demand = vec(H2["D_h"])
    mat_output = [Timestep REC_h_price REC_d_to_h_price REC_m_to_h_price REC_y_to_h_price transpose(rec_prod_h) transpose(rec_prod_d) transpose(rec_prod_m) transpose(rec_prod_y)]
    CSV.write(
        joinpath(
            home_dir,
            string("Results_", data["nReprDays"], "_repr_days"),
            string("Scenario_", data["scen_number"], "_operational_data_rec_", sens, ".csv")
        ),
        DataFrame(mat_output, :auto),
        delim=";",
        header=[
            "Timestep"; 
            "REC_price_h";
            "REC_price_d";
            "REC_price_m";
            "REC_price_y";
            string.(agents[:rec],"_h");
            string.(agents[:rec],"_d");
            string.(agents[:rec],"_m");
            string.(agents[:rec],"_y")
            ]
    )
    

end
