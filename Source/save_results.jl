# Save results
function save_results(mdict::Dict,EOM::Dict,ETS::Dict,ADMM::Dict,results::Dict,data::Dict,agents::Dict,scenario_overview_row::DataFrameRow) 
    Years = range(2017,stop=2017+data["nyears"]-1)
    Iterations = range(1,stop=data["CircularBufferSize"])

    # Aggregate metrics 
    WBseal = Years[findfirst(ETS["TNAC"] .< ETS["TNAC_MAX"])]
    if scenario_overview_row[:op_dem] != 0  # overlapping policy that affects emissions
        # to be completed
        # # Reference results
        # ref_results = CSV.read(joinpath(home_dir,"Results",string("Scenario_",scenario_overview_row[:ref_scen_number],".csv")),DataFrame;delim=";")
        # # WBL 
        # WBL =(sum(ETS["C"])-sum(ref_results[:Cancellation]))/(scenario_overview_row[:op_dem]*(scenario_overview_row[:stop_op]-scenario_overview_row[:start_op]+1))
        # # Hypothetical emission profile from reference scenario, corrected for overlapping policy
        # ref_e = ref_results[:Emissions] 
        # for y in range(scenario_overview_row[:start_op]-2016, stop=scenario_overview_row[:stop_op]-2016)
        #     ref_e[y] = ref_e[y] - scenario_overview_row[:op_dem]
        # end        
        # # ETS: compute hypothetical actions MSR under hypothetical emission profile 
        # ref_ETS = Dict() 
        # define_ETS_parameters!(ref_ETS,data,scenario_overview[scen_number,:])       
        # update_supply!(ref_e,ref_ETS,data,scenario_overview[scen_number,:]) 
        # # Calculate direct waterbed leakage: change in supply of allowances as result of overlapping policy
        # dirWBL = (sum(ref_ETS["C"])-sum(ref_results[:Cancellation]))/(scenario_overview_row[:op_dem]*(scenario_overview_row[:stop_op]-scenario_overview_row[:start_op]+1))
        # # Calculate indirect waterbed leakage 
        # indirWBL = WBL - dirWBL
    elseif scenario_overview_row[:op_supply] != 0 
        # to be completed
        # # Reference results
        # ref_results = CSV.read(joinpath(home_dir,"Results",string("Scenario_",scenario_overview_row[:ref_scen_number],".csv")),DataFrame;delim=";")
        # # WBL 
        # WBL = - (sum(ETS["C"])-sum(ref_results[:Cancellation]))/(scenario_overview_row[:op_supply]*(scenario_overview_row[:stop_op]-scenario_overview_row[:start_op]+1))  
        # # ETS: compute hypothetical actions MSR under reference emission profile 
        # ref_ETS = Dict() 
        # define_ETS_parameters!(ref_ETS,data,scenario_overview[scen_number,:])       
        # update_supply!(ref_results[:Emissions],ref_ETS,data,scenario_overview[scen_number,:]) 
        # # Calculate direct waterbed leakage: change in supply of allowances as result of overlapping policy
        # dirWBL = - (sum(ref_ETS["C"])-sum(ref_results[:Cancellation]))/(scenario_overview_row[:op_supply]*(scenario_overview_row[:stop_op]-scenario_overview_row[:start_op]+1))
        # # Calculate indirect waterbed leakage 
        # indirWBL = WBL - dirWBL
    else
        WBL = NaN
        dirWBL = NaN
        indirWBL = NaN
    end
    vector_output = [scenario_overview_row["scen_number"]; ADMM["n_iter"]; ADMM["walltime"];ADMM["Residuals"]["Primal"]["ETS"][end];ADMM["Residuals"]["Primal"]["MSR"][end]; 
                     ADMM["Residuals"]["Primal"]["EOM"][end];ADMM["Residuals"]["Primal"]["REC"][end]; ADMM["Residuals"]["Dual"]["ETS"][end];
                     ADMM["Residuals"]["Dual"]["EOM"][end]; ADMM["Residuals"]["Dual"]["REC"][end]; mdict["Ind"].ext[:parameters][:β];
                     results[ "λ"]["EUA"][end][5]; sum(results["e"][m][end][jy] for m in agents[:ets],jy in mdict[agents[:ps][1]].ext[:sets][:JY]); sum(ETS["C"]); WBseal; WBL; dirWBL; indirWBL]
    CSV.write(joinpath(home_dir,"overview_results.csv"), DataFrame(reshape(vector_output,1,:),:auto), delim=";",append=true);

    # ADMM Convergence
    mat_output = [Iterations ADMM["Residuals"]["Primal"]["ETS"][1:end] ADMM["Residuals"]["Primal"]["MSR"][1:end] ADMM["Residuals"]["Primal"]["EOM"][1:end] ADMM["Residuals"]["Primal"]["REC"][1:end]  ADMM["Residuals"]["Dual"]["ETS"][1:end] ADMM["Residuals"]["Dual"]["EOM"][1:end]  ADMM["Residuals"]["Dual"]["REC"][1:end]]
    CSV.write(joinpath(home_dir,"Results",string("Scenario_",scenario_overview_row["scen_number"],"_convergence.csv")), DataFrame(mat_output,:auto), delim=";",header=["Iterations";"PrimalResidual_ETS";"PrimalResidual_MSR";"PrimalResidual_EOM";
              "PrimalResidual_REC";"DualResidual_ETS";"DualResidual_EOM";"DualResidual_REC"])

     # ETS
     mat_output = [Years ETS["CAP"] ETS["S"] sum(ETS["C"][:,:],dims=2) ETS["MSR"][:,12] ETS["TNAC"] sum(results["e"][m][end] for m in agents[:ind]) sum(results["e"][m][end] for m in setdiff(agents[:ets],agents[:ind])) results[ "λ"]["EUA"][end] sum(results["b"][m][end] for m in agents[:ind]) sum(results["b"][m][end] for m in setdiff(agents[:ets],agents[:ind]))]
     CSV.write(joinpath(home_dir,"Results",string("Scenario_",scenario_overview_row["scen_number"],"_ETS.csv")), DataFrame(mat_output,:auto), delim=";",header=["Year";"CAP";"Supply";"Cancellation";"MSR";"TNAC";"Emissions_Ind"; "Emissions_PS"; "EUAprice"; "EUAs_Ind"; "EUAs_PS"]);
     
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
    λ_EOM_avg = [sum(EOM["Dw"][:,:,jy].*results[ "λ"]["EOM"][end][:,:,jy])./sum(EOM["Dw"][:,:,jy]) for jy in mdict[agents[:ps][1]].ext[:sets][:JY]]
    mat_output = [Years λ_EOM_avg results[ "λ"]["REC"][end] transpose(available_cap) transpose(fuel_shares)]
    CSV.write(joinpath(home_dir,"Results",string("Scenario_",scenario_overview_row["scen_number"],"_PS.csv")), DataFrame(mat_output,:auto), delim=";",header=["Year";"EOM_avg";"REC";string.("CAP_",agents[:ps]);string.("FS_",agents[:ps])]);
    
end