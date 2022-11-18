# ADMM 
function ADMM!(results::Dict,ADMM::Dict,ETS::Dict,EOM::Dict,REC::Dict,H2::Dict,H2CN_prod::Dict,H2CN_cap::Dict,NG::Dict,mdict::Dict,agents::Dict,scenario_overview_row::DataFrameRow,data::Dict,TO::TimerOutput)
    convergence = 0
    iterations = ProgressBar(1:data["ADMM"]["max_iter"])
    for iter in iterations
        if convergence == 0
            # Multi-threaded version
            @sync for m in agents[:all] 
                # created subroutine to allow multi-treading to solve agents' decision problems
                @spawn ADMM_subroutine!(m,results,ADMM,ETS,EOM,REC,H2,H2CN_prod,H2CN_cap,NG,mdict[m],agents,scenario_overview_row,TO)
            end

            # Single-threaded version
            # for m in agents[:all] 
            #     # created subroutine to allow multi-treading to solve agents' decision problems
            #     ADMM_subroutine!(m,results,ADMM,ETS,EOM,REC,H2,H2CN_prod,H2CN_cap,NG,mdict[m],agents,scenario_overview_row,TO)
            # end

            # Update supply of allowances 
            @timeit TO "Update EUA supply" begin
                update_supply!(sum(results["e"][m][end] for m in agents[:ets]),ETS,merge(data["General"],data["ETS"]),scenario_overview_row)
                push!(results["s"],copy(ETS["S"]))
            end

            # Imbalances 
            @timeit TO "Compute imbalances" begin
                push!(ADMM["Imbalances"]["ETS"], results["s"][end]-sum(results["b"][m][end] for m in agents[:ets]))
                push!(ADMM["Imbalances"]["EOM"], sum(results["g"][m][end] for m in agents[:eom]) - EOM["D"][:,:,:])
                if scenario_overview_row["Additionality"] == "Yearly" || scenario_overview_row["Additionality"] == "NA"
                    push!(ADMM["Imbalances"]["REC_y"], sum(results["r_y"][m][end] for m in agents[:rec]) + REC["RS_other_2017"]*EOM["D_cum"][1]*ceil.(REC["RT"]) - REC["RT"][:].*EOM["D_cum"][:])
                elseif scenario_overview_row["Additionality"] == "Daily"
                    push!(ADMM["Imbalances"]["REC_d"], sum(results["r_d"][m][end] for m in agents[:rec]))
                    push!(ADMM["Imbalances"]["REC_y"], sum(results["r_y"][m][end] for m in agents[:rec]) + REC["RS_other_2017"]*EOM["D_cum"][1]*ceil.(REC["RT"]) - REC["RT"][:].*EOM["D_cum"][:])
                elseif scenario_overview_row["Additionality"] == "Hourly"
                    push!(ADMM["Imbalances"]["REC_h"], sum(results["r_h"][m][end] for m in agents[:rec]))
                    push!(ADMM["Imbalances"]["REC_y"], sum(results["r_y"][m][end] for m in agents[:rec]) + REC["RS_other_2017"]*EOM["D_cum"][1]*ceil.(REC["RT"]) - REC["RT"][:].*EOM["D_cum"][:])
                end
                push!(ADMM["Imbalances"]["MSR"], results["s"][end]-results["s"][end-1])
                push!(ADMM["Imbalances"]["H2"], sum(results["h2"][m][end] for m in agents[:h2]) - H2["D"][:])
                push!(ADMM["Imbalances"]["H2CN_prod"], sum(results["h2cn_prod"][m][end] for m in agents[:h2cn_prod]) - H2CN_prod["H2CN_PRODT"][:])
                push!(ADMM["Imbalances"]["H2CN_cap"], sum(results["h2cn_cap"][m][end] for m in agents[:h2cn_cap]) - H2CN_cap["H2CN_CAPT"][:])
            end

            # Note on the EOM residuals: 
            # Penalty factors related to the EOM are weighted in objectives with weight representative days (as balancing constraint appears X times in orignal optimization problem),
            # but not in the calculation of the primal and dual residuals as this adds computational effort (inefficient operations on large arrays).
            # Primal and duals are, however, only used to measure convergence (i.e., tolerance on meeting the balancing constraints, expressed as percentage of demand 
            # on the considered representative days) and to update rho (based on ratio of primal and dual residuals, which doesn't change with the scaling of the representative days)

            # Primal residuals 
            @timeit TO "Compute primal residuals" begin
                push!(ADMM["Residuals"]["Primal"]["ETS"], sqrt(sum(ADMM["Imbalances"]["ETS"][end].^2)))
                push!(ADMM["Residuals"]["Primal"]["MSR"], sqrt(sum(ADMM["Imbalances"]["MSR"][end].^2)))
                push!(ADMM["Residuals"]["Primal"]["EOM"], sqrt(sum(ADMM["Imbalances"]["EOM"][end].^2)))
                push!(ADMM["Residuals"]["Primal"]["REC"], sqrt(sum(ADMM["Imbalances"]["REC_y"][end].^2)) + sqrt(sum(ADMM["Imbalances"]["REC_d"][end].^2)) +  sqrt(sum(ADMM["Imbalances"]["REC_h"][end].^2)))
                push!(ADMM["Residuals"]["Primal"]["H2"], sqrt(sum(ADMM["Imbalances"]["H2"][end].^2)))
                push!(ADMM["Residuals"]["Primal"]["H2CN_prod"], sqrt(sum(ADMM["Imbalances"]["H2CN_prod"][end].^2)))
                push!(ADMM["Residuals"]["Primal"]["H2CN_cap"], sqrt(sum(ADMM["Imbalances"]["H2CN_cap"][end].^2)))
            end

            # Dual residuals
            @timeit TO "Compute dual residuals" begin 
            if iter > 1
                push!(ADMM["Residuals"]["Dual"]["ETS"], sqrt(sum(sum((ADMM["ρ"]["EUA"][end]*((results["b"][m][end]-sum(results["b"][mstar][end] for mstar in agents[:ets])./(ETS["nAgents"]+1)) - (results["b"][m][end-1]-sum(results["b"][mstar][end-1] for mstar in agents[:ets])./(ETS["nAgents"]+1)))).^2 for m in agents[:ets])))) 
                push!(ADMM["Residuals"]["Dual"]["EOM"], sqrt(sum(sum((ADMM["ρ"]["EOM"][end]*((results["g"][m][end]-sum(results["g"][mstar][end] for mstar in agents[:eom])./(EOM["nAgents"]+1)) - (results["g"][m][end-1]-sum(results["g"][mstar][end-1] for mstar in agents[:eom])./(EOM["nAgents"]+1)))).^2 for m in agents[:eom]))))               
                if scenario_overview_row["Additionality"] == "NA" || scenario_overview_row["Additionality"] == "Yearly"
                    push!(ADMM["Residuals"]["Dual"]["REC"], sqrt(sum(sum((ADMM["ρ"]["REC_y"][end]*((results["r_y"][m][end]-sum(results["r_y"][mstar][end] for mstar in agents[:rec])./(REC["nAgents"]+1)) - (results["r_y"][m][end-1]-sum(results["r_y"][mstar][end-1] for mstar in agents[:rec])./(REC["nAgents"]+1)))).^2 for m in agents[:rec]))))    
                elseif scenario_overview_row["Additionality"] == "Daily"
                    push!(ADMM["Residuals"]["Dual"]["REC"], sqrt(sum(sum((ADMM["ρ"]["REC_d"][end]*((results["r_d"][m][end]-sum(results["r_d"][mstar][end] for mstar in agents[:rec])./(REC["nAgents"]+1)) - (results["r_d"][m][end-1]-sum(results["r_d"][mstar][end-1] for mstar in agents[:rec])./(REC["nAgents"]+1)))).^2 for m in agents[:rec]))) 
                    + sqrt(sum(sum((ADMM["ρ"]["REC_y"][end]*((results["r_y"][m][end]-sum(results["r_y"][mstar][end] for mstar in agents[:rec])./(REC["nAgents"]+1)) - (results["r_y"][m][end-1]-sum(results["r_y"][mstar][end-1] for mstar in agents[:rec])./(REC["nAgents"]+1)))).^2 for m in agents[:rec]))))              
                elseif scenario_overview_row["Additionality"] == "Hourly"
                    push!(ADMM["Residuals"]["Dual"]["REC"], sqrt(sum(sum((ADMM["ρ"]["REC_h"][end]*((results["r_h"][m][end]-sum(results["r_h"][mstar][end] for mstar in agents[:rec])./(REC["nAgents"]+1)) - (results["r_h"][m][end-1]-sum(results["r_h"][mstar][end-1] for mstar in agents[:rec])./(REC["nAgents"]+1)))).^2 for m in agents[:rec]))) 
                    + sqrt(sum(sum((ADMM["ρ"]["REC_y"][end]*((results["r_y"][m][end]-sum(results["r_y"][mstar][end] for mstar in agents[:rec])./(REC["nAgents"]+1)) - (results["r_y"][m][end-1]-sum(results["r_y"][mstar][end-1] for mstar in agents[:rec])./(REC["nAgents"]+1)))).^2 for m in agents[:rec]))))             
                end                               
                push!(ADMM["Residuals"]["Dual"]["H2"], sqrt(sum(sum((ADMM["ρ"]["H2"][end]*((results["h2"][m][end]-sum(results["h2"][mstar][end] for mstar in agents[:h2])./(H2["nAgents"]+1)) - (results["h2"][m][end-1]-sum(results["h2"][mstar][end-1] for mstar in agents[:h2])./(H2["nAgents"]+1)))).^2 for m in agents[:h2]))))
                push!(ADMM["Residuals"]["Dual"]["H2CN_prod"], sqrt(sum(sum((ADMM["ρ"]["H2CN_prod"][end]*((results["h2cn_prod"][m][end]-sum(results["h2cn_prod"][mstar][end] for mstar in agents[:h2cn_prod])./(H2CN_prod["nAgents"]+1)) - (results["h2cn_prod"][m][end-1]-sum(results["h2cn_prod"][mstar][end-1] for mstar in agents[:h2cn_prod])./(H2CN_prod["nAgents"]+1)))).^2 for m in agents[:h2cn_prod]))))
                push!(ADMM["Residuals"]["Dual"]["H2CN_cap"], sqrt(sum(sum((ADMM["ρ"]["H2CN_cap"][end]*((results["h2cn_cap"][m][end]-sum(results["h2cn_cap"][mstar][end] for mstar in agents[:h2cn_cap])./(H2CN_cap["nAgents"]+1)) - (results["h2cn_cap"][m][end-1]-sum(results["h2cn_cap"][mstar][end-1] for mstar in agents[:h2cn_cap])./(H2CN_cap["nAgents"]+1)))).^2 for m in agents[:h2cn_cap]))))
            end
            end

            # Price updates 
            # In general, price caps or floors can be imposed, but may slow down convergence (a negative price gives a stronger incentive than a zero price).
            # For hydrogen, a price floor (0) has been imposed to avoid negative hydrogen prices. This may arise when all hydrogen is provided by "carbon-neutral" hydrogen production routes. In that case, the net incentive is the difference between the carbon-neutral hydrogen premium and the hydrogen price. As no other agent is responding to the hydrogen price, the hydrogen prices can drop below zero if the carbon-neutral hydrogen premium is sufficiently high.
            # Such interactions do not occur in other markets.
            @timeit TO "Update prices" begin
                if scenario_overview_row[:ref_scen_number] == scenario_overview_row[:scen_number] # calibration run, 2017-2018 ETS prices fixed to historical values, 2019 to be calibrated
                    push!(results[ "λ"]["EUA"], [ETS["P_2017"]; ETS["P_2018"]; results[ "λ"]["EUA"][end][3:end] - ADMM["ρ"]["EUA"][end]/(100*data["General"]["nReprDays"])*ADMM["Imbalances"]["ETS"][end][3:end]])    
                else # 2017-2019 ETS prices fixed to historical values
                    push!(results[ "λ"]["EUA"], [ETS["P_2017"]; ETS["P_2018"]; ETS["P_2019"]; results[ "λ"]["EUA"][end][4:end] - ADMM["ρ"]["EUA"][end]/(100*data["General"]["nReprDays"])*ADMM["Imbalances"]["ETS"][end][4:end]])    
                end
                push!(results[ "λ"]["EOM"], results[ "λ"]["EOM"][end] - ADMM["ρ"]["EOM"][end]*ADMM["Imbalances"]["EOM"][end])
                push!(results[ "λ"]["REC_y"], results[ "λ"]["REC_y"][end] - ADMM["ρ"]["REC_y"][end]/(100*data["General"]["nReprDays"])*ADMM["Imbalances"]["REC_y"][end])
                push!(results[ "λ"]["REC_d"], results[ "λ"]["REC_d"][end] - ADMM["ρ"]["REC_d"][end]/(10*data["General"]["nReprDays"])*ADMM["Imbalances"]["REC_d"][end])
                push!(results[ "λ"]["REC_h"], results[ "λ"]["REC_h"][end] - ADMM["ρ"]["REC_h"][end]*ADMM["Imbalances"]["REC_h"][end])
                push!(results[ "λ"]["H2"], [maximum([0,results[ "λ"]["H2"][end][jy] - ADMM["ρ"]["H2"][end]/(100*data["General"]["nReprDays"])*data["H2"]["conv_factor"]*ADMM["Imbalances"]["H2"][end][jy]]) for jy in mdict[agents[:all][end]].ext[:sets][:JY]])
                push!(results[ "λ"]["H2CN_prod"], results[ "λ"]["H2CN_prod"][end] - ADMM["ρ"]["H2CN_prod"][end]/(100*data["General"]["nReprDays"])*data["H2"]["conv_factor"]*ADMM["Imbalances"]["H2CN_prod"][end])
                push!(results[ "λ"]["H2CN_cap"], results[ "λ"]["H2CN_cap"][end] - ADMM["ρ"]["H2CN_cap"][end]*ADMM["Imbalances"]["H2CN_cap"][end])
                push!(results[ "λ"]["NG"], NG["λ"])
            end

            # Update ρ-values
            @timeit TO "Update ρ" begin
                 update_rho!(ADMM,iter)
            end

            # Progress bar
            @timeit TO "Progress bar" begin
                set_description(iterations, string(@sprintf("ΔETS: %.3f -- ΔMSR: %.3f -- ΔEOM %.3f -- ΔREC %.3f -- ΔH2 %.3f -- ΔH2CN-PROD %.3f -- ΔH2CN-CAP %.3f ",  ADMM["Residuals"]["Primal"]["ETS"][end]/ADMM["Tolerance"]["ETS"],ADMM["Residuals"]["Primal"]["MSR"][end]/ADMM["Tolerance"]["ETS"], ADMM["Residuals"]["Primal"]["EOM"][end]/ADMM["Tolerance"]["EOM"],ADMM["Residuals"]["Primal"]["REC"][end]/ADMM["Tolerance"]["REC"],ADMM["Residuals"]["Primal"]["H2"][end]/ADMM["Tolerance"]["H2"],ADMM["Residuals"]["Primal"]["H2CN_prod"][end]/ADMM["Tolerance"]["H2CN_prod"],ADMM["Residuals"]["Primal"]["H2CN_cap"][end]/ADMM["Tolerance"]["H2CN_cap"])))
            end

            # Check convergence: primal and dual satisfy tolerance 
            if ADMM["Residuals"]["Primal"]["MSR"][end] <= ADMM["Tolerance"]["ETS"] && ADMM["Residuals"]["Primal"]["ETS"][end] <= ADMM["Tolerance"]["ETS"] && ADMM["Residuals"]["Dual"]["ETS"][end] <= ADMM["Tolerance"]["ETS"] && ADMM["Residuals"]["Primal"]["EOM"][end] <= ADMM["Tolerance"]["EOM"] && ADMM["Residuals"]["Dual"]["EOM"][end] <= ADMM["Tolerance"]["EOM"] && ADMM["Residuals"]["Primal"]["REC"][end] <= ADMM["Tolerance"]["REC"] && ADMM["Residuals"]["Dual"]["REC"][end] <= ADMM["Tolerance"]["REC"] && ADMM["Residuals"]["Primal"]["H2"][end] <= ADMM["Tolerance"]["H2"] && ADMM["Residuals"]["Dual"]["H2"][end] <= ADMM["Tolerance"]["H2"] && ADMM["Residuals"]["Primal"]["H2CN_prod"][end] <= ADMM["Tolerance"]["H2CN_prod"] && ADMM["Residuals"]["Dual"]["H2CN_prod"][end] <= ADMM["Tolerance"]["H2CN_prod"] && ADMM["Residuals"]["Primal"]["H2CN_cap"][end] <= ADMM["Tolerance"]["H2CN_cap"] && ADMM["Residuals"]["Dual"]["H2CN_cap"][end] <= ADMM["Tolerance"]["H2CN_cap"]
                convergence = 1
            end

            # store number of iterations
            ADMM["n_iter"] = copy(iter)
        end
    end
end