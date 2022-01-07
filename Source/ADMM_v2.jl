# ADMM 
function ADMM!(results::Dict,ADMM::Dict,ETS::Dict,EOM::Dict,REC::Dict,mdict::Dict,agents::Dict,scenario_overview_row::DataFrameRow,TO::TimerOutput)
    convergence = 0
    iterations = ProgressBar(1:data["ADMM"]["max_iter"])
    for iter in iterations
        if convergence == 0
            for m in agents[:all] 
                # Calculate penalty terms ADMM and update price to most recent value 
                @timeit TO "Compute ADMM penalty terms" begin
                    if mdict[m].ext[:parameters][:ETS] == 1
                        mdict[m].ext[:parameters][:b_bar] = results["b"][m][end] + 1/(ETS["nAgents"]+1)*ADMM["Imbalances"]["ETS"][end]
                        mdict[m].ext[:parameters][:λ_EUA] = results[ "λ"]["EUA"][end] 
                        mdict[m].ext[:parameters][:ρ_EUA] = ADMM["ρ"]["EUA"][end]
                    end
                    if mdict[m].ext[:parameters][:EOM] == 1
                        mdict[m].ext[:parameters][:g_bar] = results["g"][m][end] - 1/(EOM["nAgents"]+1)*ADMM["Imbalances"]["EOM"][end]
                        mdict[m].ext[:parameters][:λ_EOM] = results[ "λ"]["EOM"][end] 
                        mdict[m].ext[:parameters][:ρ_EOM] = ADMM["ρ"]["EOM"][end]
                    end
                    if mdict[m].ext[:parameters][:REC] == 1
                        mdict[m].ext[:parameters][:r_bar] = results["r"][m][end] - 1/(REC["nAgents"]+1)*ADMM["Imbalances"]["REC"][end]
                        mdict[m].ext[:parameters][:λ_REC] = results[ "λ"]["REC"][end] 
                        mdict[m].ext[:parameters][:ρ_REC] = ADMM["ρ"]["REC"][end]
                    end
                end

                # Solve agents decision problems:
                if m in agents[:ind]
                    @timeit TO "Solve industry" begin
                        update_ind_emissions!(mdict[m],merge(data["General"],data["Industry"]),ETS,scenario_overview_row) 
                        solve_ind_agent!(mdict[m])  
                    end
                    @timeit TO "Query results" begin
                        push!(results["b"][m], collect(value.(mdict[m].ext[:variables][:b])))
                        push!(results["e"][m], mdict[m].ext[:parameters][:e][:])
                    end
                elseif m in agents[:ps]
                    @timeit TO "Solve power sector" begin
                        solve_ps_agent!(mdict[m])  
                    end
                    @timeit TO "Query results" begin
                        if mdict[m].ext[:parameters][:REC] == 1
                            push!(results["r"][m], collect(value.(mdict[m].ext[:variables][:r])))
                        end
                        if mdict[m].ext[:parameters][:ETS] == 1
                            push!(results["b"][m], collect(value.(mdict[m].ext[:variables][:b])))
                            push!(results["e"][m], collect(value.(mdict[m].ext[:expressions][:e])))
                        end
                        if mdict[m].ext[:parameters][:EOM] == 1
                            push!(results["g"][m], collect(value.(mdict[m].ext[:variables][:g])))
                        end
                    end
                end
            end

            # Update supply of allowances 
            @timeit TO "Update EUA supply" begin
                update_supply!(sum(results["e"][m][end] for m in agents[:ets]),ETS,merge(data["General"],data["ETS"]),scenario_overview_row)
                push!(results["s"],copy(ETS["S"]))
            end

            # Imbalances 
            @timeit TO "Compute imbalances" begin
                push!(ADMM["Imbalances"]["ETS"], results["s"][end]-sum(results["b"][m][end] for m in agents[:ets]))
                push!(ADMM["Imbalances"]["EOM"], sum(results["g"][m][end] for m in agents[:eom]) - EOM["D"][:,:,:])
                push!(ADMM["Imbalances"]["REC"], sum(results["r"][m][end] for m in agents[:rec]) + REC["RS_other_2017"]*EOM["D_cum"][1]*ceil.(REC["RT"]) - REC["RT"][:].*EOM["D_cum"][:])
                push!(ADMM["Imbalances"]["MSR"], results["s"][end]-results["s"][end-1])
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
                push!(ADMM["Residuals"]["Primal"]["REC"], sqrt(sum(ADMM["Imbalances"]["REC"][end].^2)))
            end

            # Dual residuals
            @timeit TO "Compute dual residuals" begin 
            if iter > 1
                push!(ADMM["Residuals"]["Dual"]["ETS"], sqrt(sum(sum((ADMM["ρ"]["EUA"][end]*((results["b"][m][end]-sum(results["b"][mstar][end] for mstar in agents[:ets])./(ETS["nAgents"]+1)) - (results["b"][m][end-1]-sum(results["b"][mstar][end-1] for mstar in agents[:ets])./(ETS["nAgents"]+1)))).^2 for m in agents[:ets])))) 
                push!(ADMM["Residuals"]["Dual"]["EOM"], sqrt(sum(sum((ADMM["ρ"]["EOM"][end]*((results["g"][m][end]-sum(results["g"][mstar][end] for mstar in agents[:eom])./(EOM["nAgents"]+1)) - (results["g"][m][end-1]-sum(results["g"][mstar][end-1] for mstar in agents[:eom])./(EOM["nAgents"]+1)))).^2 for m in agents[:eom]))))               
                push!(ADMM["Residuals"]["Dual"]["REC"], sqrt(sum(sum((ADMM["ρ"]["REC"][end]*((results["r"][m][end]-sum(results["r"][mstar][end] for mstar in agents[:rec])./(REC["nAgents"]+1)) - (results["r"][m][end-1]-sum(results["r"][mstar][end-1] for mstar in agents[:rec])./(REC["nAgents"]+1)))).^2 for m in agents[:rec]))))
            end
            end

            # Price updates 
            @timeit TO "Update prices" begin
                if scenario_overview_row[:ref_scen_number] == scenario_overview_row[:scen_number] # calibration run, 2017-2018 prices fixed, 2019 prices are free variable
                    push!(results[ "λ"]["EUA"], [ETS["P_2017"]; ETS["P_2018"]; results[ "λ"]["EUA"][end][3:end] - ADMM["ρ"]["EUA"][end]/1000*ADMM["Imbalances"]["ETS"][end][3:end]])    
                else # 2017-2019 ETS prices fixed to historical values
                    push!(results[ "λ"]["EUA"], [ETS["P_2017"]; ETS["P_2018"]; ETS["P_2019"]; results[ "λ"]["EUA"][end][4:end] - ADMM["ρ"]["EUA"][end]/1000*ADMM["Imbalances"]["ETS"][end][4:end]])    
                end
                push!(results[ "λ"]["EOM"], results[ "λ"]["EOM"][end] - ADMM["ρ"]["EOM"][end]*ADMM["Imbalances"]["EOM"][end])
                push!(results[ "λ"]["REC"], results[ "λ"]["REC"][end] - ADMM["ρ"]["EUA"][end]/1000*ADMM["Imbalances"]["REC"][end])
            end

            # Update ρ-values
            @timeit TO "Update ρ" begin
                update_rho!(ADMM,iter)
            end

            # Progress bar
            @timeit TO "Progress bar" begin
                set_description(iterations, string(@sprintf("ΔETS: %.3f -- ΔMSR: %.3f -- ΔEOM %.3f -- ΔREC %.3f ",  ADMM["Residuals"]["Primal"]["ETS"][end]/ADMM["Tolerance"]["ETS"],ADMM["Residuals"]["Primal"]["MSR"][end]/ADMM["Tolerance"]["ETS"], ADMM["Residuals"]["Primal"]["EOM"][end]/ADMM["Tolerance"]["EOM"],ADMM["Residuals"]["Primal"]["REC"][end]/ADMM["Tolerance"]["REC"])))
            end

            # Check convergence: primal and dual satisfy tolerance 
            if ADMM["Residuals"]["Primal"]["MSR"][end] <= ADMM["Tolerance"]["ETS"] && ADMM["Residuals"]["Primal"]["ETS"][end] <= ADMM["Tolerance"]["ETS"] && ADMM["Residuals"]["Dual"]["ETS"][end] <= ADMM["Tolerance"]["ETS"] && ADMM["Residuals"]["Primal"]["EOM"][end] <= ADMM["Tolerance"]["EOM"] && ADMM["Residuals"]["Dual"]["EOM"][end] <= ADMM["Tolerance"]["EOM"] && ADMM["Residuals"]["Primal"]["REC"][end] <= ADMM["Tolerance"]["REC"] && ADMM["Residuals"]["Dual"]["REC"][end] <= ADMM["Tolerance"]["REC"]
                convergence = 1
            end

            # store number of iterations
            ADMM["n_iter"] = copy(iter)
        end
    end
end