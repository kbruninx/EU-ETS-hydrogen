# ADMM 
function ADMM!(results::Dict,ADMM::Dict,ETS::Dict,EOM::Dict,REC::Dict,H2::Dict,H2CN_prod::Dict,H2CN_cap::Dict,NG::Dict,mdict::Dict,agents::Dict,data::Dict,TO::TimerOutput)
    convergence = 0
    iterations = ProgressBar(1:data["ADMM"]["max_iter"])
    for iter in iterations
        if convergence == 0
            # Multi-threaded version
            @sync for m in agents[:all] 
                # created subroutine to allow multi-treading to solve agents' decision problems
                @spawn ADMM_subroutine!(m,data,results,ADMM,ETS,EOM,REC,H2,H2CN_prod,H2CN_cap,NG,mdict[m],agents,TO)
            end
            # Update supply of allowances 
                update_supply!(sum(results["e"][m][end] for m in agents[:ets]),ETS,merge(data["General"],data["ETS"],data["scenario"]))
                push!(results["s"],copy(ETS["S"]))

            # Imbalances 
                push!(ADMM["Imbalances"]["ETS"], results["s"][end]-sum(results["b"][m][end] for m in agents[:ets]))
                push!(ADMM["Imbalances"]["EOM"], sum(results["g"][m][end] for m in agents[:eom]) - EOM["D"][:,:,:])
                push!(ADMM["Imbalances"]["REC_y"], sum(results["r_y"][m][end] for m in agents[:rec]) + REC["RS_other_2021"]*EOM["D_cum"][1]*ceil.(REC["RT"]) - REC["RT"][:].*EOM["D_cum"][:])
                if data["scenario"]["Additionality"] == "Monthly"
                    push!(ADMM["Imbalances"]["REC_m"], sum(results["r_m"][m][end] for m in agents[:rec]))
                elseif data["scenario"]["Additionality"] == "Daily"
                    push!(ADMM["Imbalances"]["REC_d"], sum(results["r_d"][m][end] for m in agents[:rec]))
                elseif data["scenario"]["Additionality"] == "Hourly"
                    push!(ADMM["Imbalances"]["REC_h"], sum(results["r_h"][m][end] for m in agents[:rec]))
                end
                push!(ADMM["Imbalances"]["MSR"], results["s"][end]-results["s"][end-1])
                if data["scenario"]["H2_balance"] == "Hourly"
                    push!(ADMM["Imbalances"]["H2_h"], sum(results["h2_h"][m][end] for m in agents[:h2]) - H2["D_h"][:,:,:])
                elseif data["scenario"]["H2_balance"] == "Daily"
                    push!(ADMM["Imbalances"]["H2_d"], sum(results["h2_d"][m][end] for m in agents[:h2]) - H2["D_d"][:,:])
                elseif data["scenario"]["H2_balance"] == "Monthly"
                    push!(ADMM["Imbalances"]["H2_m"], sum(results["h2_m"][m][end] for m in agents[:h2]) - H2["D_m"][:,:])
                elseif data["scenario"]["H2_balance"] == "Yearly"
                    push!(ADMM["Imbalances"]["H2_y"], sum(results["h2_y"][m][end] for m in agents[:h2]) - H2["D_y"][:])
                end
                push!(ADMM["Imbalances"]["H2CN_prod"], sum(results["h2cn_prod"][m][end] for m in agents[:h2cn_prod]) - H2CN_prod["H2CN_PRODT"][:])
                push!(ADMM["Imbalances"]["H2CN_cap"], sum(results["h2cn_cap"][m][end] for m in agents[:h2cn_cap]) - H2CN_cap["H2CN_CAPT"][:])

            # Primal residuals 
                push!(ADMM["Residuals"]["Primal"]["ETS"], sqrt(sum(ADMM["Imbalances"]["ETS"][end].^2)))
                push!(ADMM["Residuals"]["Primal"]["MSR"], sqrt(sum(ADMM["Imbalances"]["MSR"][end].^2)))
                push!(ADMM["Residuals"]["Primal"]["EOM"], sqrt(sum(ADMM["Imbalances"]["EOM"][end].^2)))
                push!(ADMM["Residuals"]["Primal"]["REC_y"], sqrt(sum(ADMM["Imbalances"]["REC_y"][end].^2))) 
                push!(ADMM["Residuals"]["Primal"]["REC_m"], sqrt(sum(ADMM["Imbalances"]["REC_m"][end].^2))) 
                push!(ADMM["Residuals"]["Primal"]["REC_d"], sqrt(sum(ADMM["Imbalances"]["REC_d"][end].^2))) 
                push!(ADMM["Residuals"]["Primal"]["REC_h"], sqrt(sum(ADMM["Imbalances"]["REC_h"][end].^2)))
                push!(ADMM["Residuals"]["Primal"]["H2_h"], sqrt(sum(ADMM["Imbalances"]["H2_h"][end].^2)))
                push!(ADMM["Residuals"]["Primal"]["H2_d"], sqrt(sum(ADMM["Imbalances"]["H2_d"][end].^2)))
                push!(ADMM["Residuals"]["Primal"]["H2_m"], sqrt(sum(ADMM["Imbalances"]["H2_m"][end]).^2))
                push!(ADMM["Residuals"]["Primal"]["H2_y"], sqrt(sum(ADMM["Imbalances"]["H2_y"][end].^2)))
                push!(ADMM["Residuals"]["Primal"]["H2CN_prod"], sqrt(sum(ADMM["Imbalances"]["H2CN_prod"][end].^2)))
                push!(ADMM["Residuals"]["Primal"]["H2CN_cap"], sqrt(sum(ADMM["Imbalances"]["H2CN_cap"][end].^2)))

            # Dual residuals
            if iter > 1
                push!(ADMM["Residuals"]["Dual"]["ETS"], sqrt(sum(sum((ADMM["ρ"]["EUA"][end]*((results["b"][m][end]-sum(results["b"][mstar][end] for mstar in agents[:ets])./(ETS["nAgents"]+1)) - (results["b"][m][end-1]-sum(results["b"][mstar][end-1] for mstar in agents[:ets])./(ETS["nAgents"]+1)))).^2 for m in agents[:ets])))) 
                push!(ADMM["Residuals"]["Dual"]["EOM"], sqrt(sum(sum((ADMM["ρ"]["EOM"][end]*((results["g"][m][end]-sum(results["g"][mstar][end] for mstar in agents[:eom])./(EOM["nAgents"]+1)) - (results["g"][m][end-1]-sum(results["g"][mstar][end-1] for mstar in agents[:eom])./(EOM["nAgents"]+1)))).^2 for m in agents[:eom]))))               
                push!(ADMM["Residuals"]["Dual"]["REC_y"], sqrt(sum(sum((ADMM["ρ"]["REC_y"][end]*((results["r_y"][m][end]-sum(results["r_y"][mstar][end] for mstar in agents[:rec])./(REC["nAgents"]+1)) - (results["r_y"][m][end-1]-sum(results["r_y"][mstar][end-1] for mstar in agents[:rec])./(REC["nAgents"]+1)))).^2 for m in agents[:rec]))))    
                if data["scenario"]["Additionality"] == "Monthly"
                    push!(ADMM["Residuals"]["Dual"]["REC_m"], sqrt(sum(sum((ADMM["ρ"]["REC_m"][end]*((results["r_m"][m][end]-sum(results["r_m"][mstar][end] for mstar in agents[:rec])./(REC["nAgents"]+1)) - (results["r_m"][m][end-1]-sum(results["r_m"][mstar][end-1] for mstar in agents[:rec])./(REC["nAgents"]+1)))).^2 for m in agents[:rec])))) 
                elseif data["scenario"]["Additionality"] == "Daily"
                    push!(ADMM["Residuals"]["Dual"]["REC_d"], sqrt(sum(sum((ADMM["ρ"]["REC_d"][end]*((results["r_d"][m][end]-sum(results["r_d"][mstar][end] for mstar in agents[:rec])./(REC["nAgents"]+1)) - (results["r_d"][m][end-1]-sum(results["r_d"][mstar][end-1] for mstar in agents[:rec])./(REC["nAgents"]+1)))).^2 for m in agents[:rec]))))              
                elseif data["scenario"]["Additionality"] == "Hourly"
                    push!(ADMM["Residuals"]["Dual"]["REC_h"], sqrt(sum(sum((ADMM["ρ"]["REC_h"][end]*((results["r_h"][m][end]-sum(results["r_h"][mstar][end] for mstar in agents[:rec])./(REC["nAgents"]+1)) - (results["r_h"][m][end-1]-sum(results["r_h"][mstar][end-1] for mstar in agents[:rec])./(REC["nAgents"]+1)))).^2 for m in agents[:rec]))))             
                end       
                if data["scenario"]["H2_balance"] == "Hourly"
                    push!(ADMM["Residuals"]["Dual"]["H2_h"], sqrt(sum(sum((ADMM["ρ"]["H2_h"][end]*((results["h2_h"][m][end]-sum(results["h2_h"][mstar][end] for mstar in agents[:h2])./(H2["nAgents"]+1)) - (results["h2_h"][m][end-1]-sum(results["h2_h"][mstar][end-1] for mstar in agents[:h2])./(H2["nAgents"]+1)))).^2 for m in agents[:h2]))))
                elseif data["scenario"]["H2_balance"] == "Daily"
                    push!(ADMM["Residuals"]["Dual"]["H2_d"], sqrt(sum(sum((ADMM["ρ"]["H2_d"][end]*((results["h2_d"][m][end]-sum(results["h2_d"][mstar][end] for mstar in agents[:h2])./(H2["nAgents"]+1)) - (results["h2_d"][m][end-1]-sum(results["h2_d"][mstar][end-1] for mstar in agents[:h2])./(H2["nAgents"]+1)))).^2 for m in agents[:h2]))))
                elseif data["scenario"]["H2_balance"] == "Monthly"
                    push!(ADMM["Residuals"]["Dual"]["H2_m"], sqrt(sum(sum((ADMM["ρ"]["H2_m"][end]*((results["h2_m"][m][end]-sum(results["h2_m"][mstar][end] for mstar in agents[:h2])./(H2["nAgents"]+1)) - (results["h2_m"][m][end-1]-sum(results["h2_m"][mstar][end-1] for mstar in agents[:h2])./(H2["nAgents"]+1)))).^2 for m in agents[:h2]))))
                elseif data["scenario"]["H2_balance"] == "Yearly"
                    push!(ADMM["Residuals"]["Dual"]["H2_y"], sqrt(sum(sum((ADMM["ρ"]["H2_y"][end]*((results["h2_y"][m][end]-sum(results["h2_y"][mstar][end] for mstar in agents[:h2])./(H2["nAgents"]+1)) - (results["h2_y"][m][end-1]-sum(results["h2_y"][mstar][end-1] for mstar in agents[:h2])./(H2["nAgents"]+1)))).^2 for m in agents[:h2]))))
                end
                push!(ADMM["Residuals"]["Dual"]["H2CN_prod"], sqrt(sum(sum((ADMM["ρ"]["H2CN_prod"][end]*((results["h2cn_prod"][m][end]-sum(results["h2cn_prod"][mstar][end] for mstar in agents[:h2cn_prod])./(H2CN_prod["nAgents"]+1)) - (results["h2cn_prod"][m][end-1]-sum(results["h2cn_prod"][mstar][end-1] for mstar in agents[:h2cn_prod])./(H2CN_prod["nAgents"]+1)))).^2 for m in agents[:h2cn_prod]))))
                push!(ADMM["Residuals"]["Dual"]["H2CN_cap"], sqrt(sum(sum((ADMM["ρ"]["H2CN_cap"][end]*((results["h2cn_cap"][m][end]-sum(results["h2cn_cap"][mstar][end] for mstar in agents[:h2cn_cap])./(H2CN_cap["nAgents"]+1)) - (results["h2cn_cap"][m][end-1]-sum(results["h2cn_cap"][mstar][end-1] for mstar in agents[:h2cn_cap])./(H2CN_cap["nAgents"]+1)))).^2 for m in agents[:h2cn_cap]))))
            end

            # Price updates 
            # In general, price caps or floors can be imposed, but may slow down convergence (a negative price gives a stronger incentive than a zero price).
                if data["scenario"]["ref_scen_number"] == data["scenario"]["scen_number"] && data["scenario"]["sens_number"] == 1 # calibration run, 2021 to be calibrated
                    push!(results[ "λ"]["EUA"], results[ "λ"]["EUA"][end] - ADMM["ρ"]["EUA"][end]/100*ADMM["Imbalances"]["ETS"][end])    
                else # 2019-2021 ETS prices fixed to historical values
                    push!(results[ "λ"]["EUA"], [ETS["P_2021"]; results[ "λ"]["EUA"][end][2:end] - ADMM["ρ"]["EUA"][end]/100*ADMM["Imbalances"]["ETS"][end][2:end]])    
                end
                push!(results[ "λ"]["EOM"], results[ "λ"]["EOM"][end] - ADMM["ρ"]["EOM"][end]*ADMM["Imbalances"]["EOM"][end])
                push!(results[ "λ"]["REC_y"], results[ "λ"]["REC_y"][end] - ADMM["ρ"]["REC_y"][end]/100*ADMM["Imbalances"]["REC_y"][end])
                if data["scenario"]["Additionality"] == "Monthly"
                    push!(results[ "λ"]["REC_m"], results[ "λ"]["REC_m"][end] - ADMM["ρ"]["REC_m"][end]/10*ADMM["Imbalances"]["REC_m"][end])
                elseif data["scenario"]["Additionality"] == "Daily"
                    push!(results[ "λ"]["REC_d"], results[ "λ"]["REC_d"][end] - ADMM["ρ"]["REC_d"][end]/10*ADMM["Imbalances"]["REC_d"][end])
                elseif data["scenario"]["Additionality"] == "Hourly"
                    push!(results[ "λ"]["REC_h"], results[ "λ"]["REC_h"][end] - ADMM["ρ"]["REC_h"][end]*ADMM["Imbalances"]["REC_h"][end])
                end
                if data["scenario"]["H2_balance"] == "Hourly"
                    push!(results[ "λ"]["H2_h"], results[ "λ"]["H2_h"][end] - ADMM["ρ"]["H2_h"][end]*ADMM["Imbalances"]["H2_h"][end])
                elseif data["scenario"]["H2_balance"] == "Daily"
                    push!(results[ "λ"]["H2_d"], results[ "λ"]["H2_d"][end] - ADMM["ρ"]["H2_d"][end]/10*ADMM["Imbalances"]["H2_d"][end])
                elseif data["scenario"]["H2_balance"] == "Monthly"
                    push!(results[ "λ"]["H2_m"], results[ "λ"]["H2_m"][end] - ADMM["ρ"]["H2_m"][end]/10*ADMM["Imbalances"]["H2_m"][end])
                elseif data["scenario"]["H2_balance"] == "Yearly"
                    push!(results[ "λ"]["H2_y"], results[ "λ"]["H2_y"][end] - ADMM["ρ"]["H2_y"][end]/1000*ADMM["Imbalances"]["H2_y"][end])
                end
                push!(results[ "λ"]["H2CN_prod"], results[ "λ"]["H2CN_prod"][end] - ADMM["ρ"]["H2CN_prod"][end]/100*ADMM["Imbalances"]["H2CN_prod"][end])
                push!(results[ "λ"]["H2CN_cap"], results[ "λ"]["H2CN_cap"][end] - ADMM["ρ"]["H2CN_cap"][end]/1000*ADMM["Imbalances"]["H2CN_cap"][end])
                push!(results[ "λ"]["NG"], NG["λ"])

            # Update ρ-values
                update_rho!(ADMM,iter)

            # Check convergence: primal and dual satisfy tolerance 
            if ADMM["Residuals"]["Primal"]["MSR"][end] <= ADMM["Tolerance"]["ETS"] && ADMM["Residuals"]["Primal"]["ETS"][end] <= ADMM["Tolerance"]["ETS"] && ADMM["Residuals"]["Dual"]["ETS"][end] <= ADMM["Tolerance"]["ETS"] && ADMM["Residuals"]["Primal"]["EOM"][end] <= ADMM["Tolerance"]["EOM"] && ADMM["Residuals"]["Dual"]["EOM"][end] <= ADMM["Tolerance"]["EOM"] && ADMM["Residuals"]["Primal"]["REC_y"][end] <= ADMM["Tolerance"]["REC_y"] && ADMM["Residuals"]["Primal"]["REC_m"][end] <= ADMM["Tolerance"]["REC_m"] && ADMM["Residuals"]["Primal"]["REC_d"][end] <= ADMM["Tolerance"]["REC_d"] && ADMM["Residuals"]["Primal"]["REC_h"][end] <= ADMM["Tolerance"]["REC_h"] && ADMM["Residuals"]["Dual"]["REC_y"][end] <= ADMM["Tolerance"]["REC_y"] && ADMM["Residuals"]["Dual"]["REC_m"][end] <= ADMM["Tolerance"]["REC_m"] && ADMM["Residuals"]["Dual"]["REC_d"][end] <= ADMM["Tolerance"]["REC_d"] && ADMM["Residuals"]["Dual"]["REC_h"][end] <= ADMM["Tolerance"]["REC_h"] && ADMM["Residuals"]["Primal"]["H2_h"][end] <= ADMM["Tolerance"]["H2_h"] && ADMM["Residuals"]["Dual"]["H2_h"][end] <= ADMM["Tolerance"]["H2_h"] && ADMM["Residuals"]["Primal"]["H2_d"][end] <= ADMM["Tolerance"]["H2_d"] && ADMM["Residuals"]["Dual"]["H2_d"][end] <= ADMM["Tolerance"]["H2_d"] && ADMM["Residuals"]["Primal"]["H2_m"][end] <= ADMM["Tolerance"]["H2_m"] && ADMM["Residuals"]["Dual"]["H2_m"][end] <= ADMM["Tolerance"]["H2_m"] && ADMM["Residuals"]["Primal"]["H2_y"][end] <= ADMM["Tolerance"]["H2_y"] && ADMM["Residuals"]["Dual"]["H2_y"][end] <= ADMM["Tolerance"]["H2_y"] && ADMM["Residuals"]["Primal"]["H2CN_prod"][end] <= ADMM["Tolerance"]["H2CN_prod"] && ADMM["Residuals"]["Dual"]["H2CN_prod"][end] <= ADMM["Tolerance"]["H2CN_prod"] && ADMM["Residuals"]["Primal"]["H2CN_cap"][end] <= ADMM["Tolerance"]["H2CN_cap"] && ADMM["Residuals"]["Dual"]["H2CN_cap"][end] <= ADMM["Tolerance"]["H2CN_cap"]
                convergence = 1
            end
            # store number of iterations
            ADMM["n_iter"] = copy(iter)
        end
    end
end
