# Save results
function save_results(mdict::Dict,EOM::Dict,ETS::Dict,ADMM::Dict,results::Dict,data::Dict,agents::Dict,scenario_overview_row::DataFrameRow,sens) 
    # note that type of "sens" is not defined as a string stored in a dictionary is of type String31, whereas a "regular" string is of type String. Specifying one or the other may trow errors.
    Years = range(2019,stop=2019+data["nyears"]-1)
    Iterations = range(1,stop=data["CircularBufferSize"])

    # Aggregate metrics 
    tot_cost = sum(value(mdict[m].ext[:expressions][:tot_cost]) for m in agents[:all])
    tot_em = sum(results["e"][m][end][jy] for m in agents[:ets],jy in mdict[agents[:ps][1]].ext[:sets][:JY]) 
    
    vector_output = [scenario_overview_row["scen_number"]; sens; ADMM["n_iter"]; ADMM["walltime"];ADMM["Residuals"]["Primal"]["ETS"][end];ADMM["Residuals"]["Primal"]["MSR"][end]; 
                     ADMM["Residuals"]["Primal"]["EOM"][end];ADMM["Residuals"]["Primal"]["REC"][end]; ADMM["Residuals"]["Primal"]["H2"][end]; ADMM["Residuals"]["Primal"]["H2CN_prod"][end]; ADMM["Residuals"]["Primal"]["H2CN_cap"][end]; ADMM["Residuals"]["Dual"]["ETS"][end]; ADMM["Residuals"]["Dual"]["EOM"][end]; ADMM["Residuals"]["Dual"]["REC"][end]; ADMM["Residuals"]["Dual"]["H2"][end];ADMM["Residuals"]["Dual"]["H2CN_prod"][end]; ADMM["Residuals"]["Dual"]["H2CN_cap"][end]; mdict["Ind"].ext[:parameters][:β]; results[ "λ"]["EUA"][end][5]; tot_em; tot_cost]
    CSV.write(joinpath(home_dir,string("overview_results_",data["nReprDays"],"_repr_days.csv")), DataFrame(reshape(vector_output,1,:),:auto), delim=";",append=true);

    # ADMM Convergence 
    mat_output = [Iterations ADMM["Residuals"]["Primal"]["ETS"][1:end] ADMM["Residuals"]["Primal"]["MSR"][1:end] ADMM["Residuals"]["Primal"]["EOM"][1:end] ADMM["Residuals"]["Primal"]["REC"][1:end] ADMM["Residuals"]["Primal"]["H2"][1:end] ADMM["Residuals"]["Primal"]["H2CN_prod"][1:end] ADMM["Residuals"]["Primal"]["H2CN_cap"][1:end] ADMM["Residuals"]["Dual"]["ETS"][1:end] ADMM["Residuals"]["Dual"]["EOM"][1:end]  ADMM["Residuals"]["Dual"]["REC"][1:end] ADMM["Residuals"]["Dual"]["H2"][1:end] ADMM["Residuals"]["Dual"]["H2CN_prod"][1:end] ADMM["Residuals"]["Dual"]["H2CN_cap"][1:end]]
    CSV.write(joinpath(home_dir,string("Results_",data["nReprDays"],"_repr_days"),string("Scenario_",scenario_overview_row["scen_number"],"_convergence_",sens,".csv")), DataFrame(mat_output,:auto), delim=";",header=["Iterations";"PrimalResidual_ETS";"PrimalResidual_MSR";"PrimalResidual_EOM";"PrimalResidual_REC";"PrimalResidual_H2";"PrimalResidual_H2CN_prod";"PrimalResidual_H2CN_cap";"DualResidual_ETS";"DualResidual_EOM";"DualResidual_REC";"DualResidual_H2";"DualResidual_H2CN_prod";"DualResidual_H2CN_cap"])
   
    # ETS
    mat_output = [Years ETS["CAP"] ETS["S"] sum(ETS["C"][:,:],dims=2) ETS["MSR"][:,12] ETS["TNAC"] sum(results["e"][m][end] for m in agents[:ind]) sum(results["e"][m][end] for m in setdiff(agents[:ets],union(agents[:ind],agents[:ps]))) sum(results["e"][m][end] for m in setdiff(agents[:ets],union(agents[:ind],agents[:h2s]))) results[ "λ"]["EUA"][end] sum(results["b"][m][end] for m in agents[:ind]) sum(results["b"][m][end] for m in setdiff(agents[:ets],agents[:ind]))]
    CSV.write(joinpath(home_dir,string("Results_",data["nReprDays"],"_repr_days"),string("Scenario_",scenario_overview_row["scen_number"],"_ETS_",sens,".csv")), DataFrame(mat_output,:auto), delim=";",header=["Year";"CAP";"Supply";"Cancellation";"MSR";"TNAC";"Emissions_Ind";"Emissions_H2S"; "Emissions_PS"; "EUAprice"; "EUAs_Ind"; "EUAs_PS"]);
     
    # Power sector
    fuel_shares = zeros(length(agents[:ps]),data["nyears"])
    available_cap = zeros(length(agents[:ps]),data["nyears"])
    mm = 1
    for m in agents[:ps]
        gw = value.(mdict[m].ext[:expressions][:gw])
        fuel_shares[mm,:] = sum(gw[jh,jd,:] for jh in mdict[m].ext[:sets][:JH], jd in mdict[m].ext[:sets][:JD])
        CAP_LT = mdict[m].ext[:parameters][:CAP_LT]
        LEG_CAP = mdict[m].ext[:parameters][:LEG_CAP]
        cap = value.(mdict[m].ext[:variables][:cap])
        available_cap[mm,:] = [sum(CAP_LT[y2,jy]*cap[y2] for y2=1:jy) + LEG_CAP[jy] for jy in mdict[m].ext[:sets][:JY]]
        mm = mm+1
    end
    gw = Dict()
    for m in agents[:eom]
        gw[m] = value.(mdict[m].ext[:expressions][:gw])
    end
    gw_tot = sum(gw[m] for m in agents[:eom]) # total electricity generation
    gw_res_tot = sum(results["r_y"][m][end][:] for m in agents[:rec]) # total renewable electricty generation participating in the res auctions
    # For the production weighted averages below I assume that REC support for RES from electrolyzers (or electricity costs) are internal transfers - not taken into account:
    λ_EOM_avg = [sum(gw_tot[:,:,jy].*results[ "λ"]["EOM"][end][:,:,jy])./sum(gw_tot[:,:,jy]) for jy in mdict[agents[:ps][1]].ext[:sets][:JY]] # production weighted average electricity price
    λ_REC_avg = [gw_res_tot[jy].*results[ "λ"]["REC_y"][end][jy]./gw_res_tot[jy] for jy in mdict[agents[:ps][1]].ext[:sets][:JY]] # production weighted support for RES  
    mat_output = [Years λ_EOM_avg λ_REC_avg transpose(available_cap) transpose(fuel_shares)]
    CSV.write(joinpath(home_dir,string("Results_",data["nReprDays"],"_repr_days"),string("Scenario_",scenario_overview_row["scen_number"],"_PS_",sens,".csv")), DataFrame(mat_output,:auto), delim=";",header=["Year";"EOM_avg";"REC_y";string.("CAP_",agents[:ps]);string.("FS_",agents[:ps])]);
    
    # Hydrogen sector 
    h2_cap = zeros(length(agents[:h2s]),data["nyears"])
    h2_prod = zeros(length(agents[:h2s]),data["nyears"])
    mm = 1
    for m in agents[:h2s]
        CAP_LT = mdict[m].ext[:parameters][:CAP_LT]
        LEG_CAP = mdict[m].ext[:parameters][:LEG_CAP]
        cap = value.(mdict[m].ext[:variables][:capH])
        h2_cap[mm,:] = [sum(CAP_LT[y2,jy]*cap[y2] for y2=1:jy) + LEG_CAP[jy] for jy in mdict[m].ext[:sets][:JY]]    
        h2_prod[mm,:] = value.(mdict[m].ext[:variables][:gH])./data["conv_factor"]
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
    mat_output = [Years transpose(h2_cap) transpose(h2_prod) transpose(h2cn_cap) transpose(h2cn_prod) results["λ"]["H2"][end]*data["conv_factor"]/1000 results["λ"]["H2CN_prod"][end]*data["conv_factor"]/1000 results["λ"]["H2CN_cap"][end]]
    CSV.write(joinpath(home_dir,string("Results_",data["nReprDays"],"_repr_days"),string("Scenario_",scenario_overview_row["scen_number"],"_H2_",sens,".csv")), DataFrame(mat_output,:auto), delim=";",header=["Year";string.("CAP_",agents[:h2s]);string.("PROD_",agents[:h2s]);string.("CN_CAP_",agents[:h2cn_prod]);string.("CN_PROD_",agents[:h2cn_prod]);"PriceH2";"PremiumH2CN_prod";"PremiumH2CN_cap"]);
end
