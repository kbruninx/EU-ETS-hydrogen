# ADMM 
function ADMM!(results::Dict,ADMM::Dict,ETS::Dict,EOM::Dict,REC::Dict,mdict::Dict,agents::Dict,scenario_overview_row::DataFrameRow)
    convergence = 0
    iterations = ProgressBar(1:data["ADMM"]["max_iter"]-1)
    for iter in iterations
        if convergence == 0
            for m in agents[:all] 
                # Calculate penalty terms ADMM and update price to most recent value 
                @timeit TO "Compute ADMM penalty terms" begin
                    if iter > 1
                        mdict[m].ext[:parameters][:b_bar] = results["b"][m][iter-1,:] + 1/(ETS["nAgents"]+1)*ADMM["Imbalances"]["ETS"][iter-1,:]
                        mdict[m].ext[:parameters][:g_bar] = results["g"][m][iter-1,:,:,:] - 1/(EOM["nAgents"]+1)*ADMM["Imbalances"]["EOM"][iter-1,:,:,:]
                        mdict[m].ext[:parameters][:r_bar] = results["r"][m][iter-1,:] - 1/(REC["nAgents"]+1)*ADMM["Imbalances"]["REC"][iter-1,:]
                    end
                    mdict[m].ext[:parameters][:λ_EUA] = results["λ_EUA"][iter,:] 
                    mdict[m].ext[:parameters][:ρ_EUA] = ADMM["ρ_EUA"][iter]
                    mdict[m].ext[:parameters][:λ_EOM] = results["λ_EOM"][iter,:,:,:] 
                    mdict[m].ext[:parameters][:ρ_EOM] = ADMM["ρ_EOM"][iter]
                    mdict[m].ext[:parameters][:λ_REC] = results["λ_REC"][iter,:] 
                    mdict[m].ext[:parameters][:ρ_REC] = ADMM["ρ_REC"][iter]
                end

                # Solve agents decision problems:
                if m in agents[:ind]
                    @timeit TO "Solve industry" begin
                        update_ind_emissions!(mdict[m],merge(data["General"],data["Industry"]),scenario_overview_row) 
                        solve_ind_agent!(mdict[m])  
                    end
                    @timeit TO "Query results" begin
                        results["b"][m][iter,:] = value.(mdict[m].ext[:variables][:b]) 
                        results["e"][m][iter,:] = mdict[m].ext[:parameters][:e][:] 
                    end
                elseif m in agents[:ps]
                    @timeit TO "Solve power sector" begin
                        solve_ps_agent!(mdict[m])  
                    end
                    @timeit TO "Query results" begin
                        if mdict[m].ext[:parameters][:REC] == 1
                            results["r"][m][iter,:] = value.(mdict[m].ext[:variables][:r]) 
                        end
                        if mdict[m].ext[:parameters][:ETS] == 1
                            results["b"][m][iter,:] = value.(mdict[m].ext[:variables][:b]) 
                            results["e"][m][iter,:] = value.(mdict[m].ext[:expressions][:e]) 
                        end
                        results["g"][m][iter,:,:,:] = value.(mdict[m].ext[:variables][:g])
                    end
                end
            end

            # Update supply of allowances 
            @timeit TO "Update EUA supply" begin
                update_supply!(sum(results["e"][m][iter,:] for m in agents[:all]),ETS,merge(data["General"],data["ETS"]),scenario_overview_row)
                results["s"][iter,:] = ETS["S"]
            end
            
            # Imbalances 
            @timeit TO "Compute imbalances" begin
                ADMM["Imbalances"]["ETS"][iter,:] = results["s"][iter,:]-sum(results["b"][m][iter,:] for m in agents[:all])
                ADMM["Imbalances"]["EOM"][iter,:,:,:] = sum(results["g"][m][iter,:,:,:] for m in agents[:all]) - EOM["D"][:,:,:] 
                ADMM["Imbalances"]["REC"][iter,:] = sum(results["r"][m][iter,:] for m in agents[:all]) + REC["RS_other_2017"]*EOM["D_cum"][1]*ceil.(REC["RT"]) - REC["RT"][:].*EOM["D_cum"][:]
                if iter > 1
                    ADMM["Imbalances"]["MSR"][iter,:] = results["s"][iter,:]-results["s"][iter-1,:]
                end
            end

            # Primal residuals 
            @timeit TO "Compute primal residuals" begin
                ADMM["Residuals"]["Primal"]["ETS"][iter] = sqrt(sum(ADMM["Imbalances"]["ETS"][iter,:].^2)) 
                ADMM["Residuals"]["Primal"]["MSR"][iter] = sqrt(sum(ADMM["Imbalances"]["MSR"][iter,:].^2))
                ADMM["Residuals"]["Primal"]["EOM"][iter] = sqrt(sum(ADMM["Imbalances"]["EOM"][iter,:,:,:].^2))
                ADMM["Residuals"]["Primal"]["REC"][iter] = sqrt(sum(ADMM["Imbalances"]["REC"][iter,:].^2))
            end

            # Dual residuals
            @timeit TO "Compute dual residuals" begin 
            if iter > 1  
                ADMM["Residuals"]["Dual"]["ETS"][iter] = sqrt(sum(sum((ADMM["ρ_EUA"][iter]*((results["b"][m][iter,:]-sum(results["b"][mstar][iter,:] for mstar in agents[:all])./(ETS["nAgents"]+1)) - (results["b"][m][iter-1,:]-sum(results["b"][mstar][iter-1,:] for mstar in agents[:all])./(ETS["nAgents"]+1)))).^2 for m in agents[:all]))) 
                ADMM["Residuals"]["Dual"]["EOM"][iter] = sqrt(sum(sum((ADMM["ρ_EOM"][iter]*((results["g"][m][iter,:,:,:]-sum(results["g"][mstar][iter,:,:,:] for mstar in agents[:all])./(EOM["nAgents"]+1))-(results["g"][m][iter-1,:,:,:]-sum(results["g"][mstar][iter-1,:,:,:] for mstar in agents[:all])./(EOM["nAgents"]+1)))).^2 for m in agents[:all])))               
                ADMM["Residuals"]["Dual"]["REC"][iter] = sqrt(sum(sum((ADMM["ρ_REC"][iter]*((results["r"][m][iter,:]-sum(results["r"][mstar][iter,:] for mstar in agents[:all])./(REC["nAgents"]+1)) - (results["r"][m][iter-1,:]-sum(results["r"][mstar][iter-1,:] for mstar in agents[:all])./(REC["nAgents"]+1)))).^2 for m in agents[:all]))) 
            end
            end

            # Price updates 
            @timeit TO "Update prices" begin
                results["λ_EUA"][iter+1,:] = results["λ_EUA"][iter,:]-ADMM["ρ_EUA"][iter]/365*ADMM["Imbalances"]["ETS"][iter,:]  
                results["λ_EOM"][iter+1,:,:,:] = results["λ_EOM"][iter,:,:,:]- ADMM["ρ_EOM"][iter]*ADMM["Imbalances"]["EOM"][iter,:,:,:]
                results["λ_REC"][iter+1,:] = results["λ_REC"][iter,:]-ADMM["ρ_EUA"][iter]/365*ADMM["Imbalances"]["REC"][iter,:]
            end

            # Update ρ-values
            @timeit TO "Update ρ" begin
                update_rho!(ADMM,iter)
            end

            # Progress bar
            @timeit TO "Progress bar" begin
                set_description(iterations, string(@sprintf("ΔETS: %.3f MtCO2 -- ΔEOM %.3f TWh -- ΔREC %.3f TWh",  ADMM["Residuals"]["Primal"]["ETS"][iter], ADMM["Residuals"]["Primal"]["EOM"][iter],ADMM["Residuals"]["Primal"]["REC"][iter])))
            end

            # Check convergence: primal and dual satisfy tolerance 
            if ADMM["Residuals"]["Primal"]["MSR"][iter] <= ADMM["Tolerance"] && ADMM["Residuals"]["Primal"]["ETS"][iter] <= ADMM["Tolerance"] && ADMM["Residuals"]["Dual"]["ETS"][iter] <= ADMM["Tolerance"] && ADMM["Residuals"]["Primal"]["EOM"][iter] <= ADMM["Tolerance"] && ADMM["Residuals"]["Dual"]["EOM"][iter] <= ADMM["Tolerance"] && ADMM["Residuals"]["Primal"]["REC"][iter] <= ADMM["Tolerance"] && ADMM["Residuals"]["Dual"]["REC"][iter] <= ADMM["Tolerance"]
                convergence = 1
            end

            # store number of iterations
            ADMM["n_iter"] = copy(iter)
        end
    end
end