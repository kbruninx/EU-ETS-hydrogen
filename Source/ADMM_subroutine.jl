function ADMM_subroutine!(m::String,results::Dict,ADMM::Dict,ETS::Dict,EOM::Dict,REC::Dict,H2::Dict,H2CN_prod::Dict,H2CN_cap::Dict,NG::Dict,mdict::Dict,agents::Dict,scenario_overview_row::DataFrameRow,TO::TimerOutput)
TO_local = TimerOutput()
# Calculate penalty terms ADMM and update price to most recent value 
@timeit TO_local "Compute ADMM penalty terms" begin
    if mdict[m].ext[:parameters][:ETS] == 1
        mdict[m].ext[:parameters][:b_bar] = results["b"][m][end] + 1/(ETS["nAgents"]+1)*ADMM["Imbalances"]["ETS"][end]
        mdict[m].ext[:parameters][:λ_EUA] = results["λ"]["EUA"][end] 
        mdict[m].ext[:parameters][:ρ_EUA] = ADMM["ρ"]["EUA"][end]
    end
    if mdict[m].ext[:parameters][:EOM] == 1
        mdict[m].ext[:parameters][:g_bar] = results["g"][m][end] - 1/(EOM["nAgents"]+1)*ADMM["Imbalances"]["EOM"][end]
        mdict[m].ext[:parameters][:λ_EOM] = results["λ"]["EOM"][end] 
        mdict[m].ext[:parameters][:ρ_EOM] = ADMM["ρ"]["EOM"][end]
    end
    if mdict[m].ext[:parameters][:REC] == 1  
        # Yearly
        mdict[m].ext[:parameters][:r_y_bar] = results["r_y"][m][end] - 1/(REC["nAgents"]+1)*ADMM["Imbalances"]["REC_y"][end]
        mdict[m].ext[:parameters][:λ_y_REC] = results["λ"]["REC_y"][end] 
        mdict[m].ext[:parameters][:ρ_y_REC] = ADMM["ρ"]["REC_y"][end]
        # Daily
        mdict[m].ext[:parameters][:r_d_bar] = results["r_d"][m][end] - 1/(REC["nAgents"]+1)*ADMM["Imbalances"]["REC_d"][end]
        mdict[m].ext[:parameters][:λ_d_REC] = results["λ"]["REC_d"][end] 
        mdict[m].ext[:parameters][:ρ_d_REC] = ADMM["ρ"]["REC_d"][end]
        # Hourly
        mdict[m].ext[:parameters][:r_h_bar] = results["r_h"][m][end] - 1/(REC["nAgents"]+1)*ADMM["Imbalances"]["REC_h"][end]
        mdict[m].ext[:parameters][:λ_h_REC] = results["λ"]["REC_h"][end] 
        mdict[m].ext[:parameters][:ρ_h_REC] = ADMM["ρ"]["REC_h"][end]
    end
    if mdict[m].ext[:parameters][:H2] == 1
        mdict[m].ext[:parameters][:gH_bar] = results["h2"][m][end] - 1/(H2["nAgents"]+1)*ADMM["Imbalances"]["H2"][end]
        mdict[m].ext[:parameters][:λ_H2] = results["λ"]["H2"][end] 
        mdict[m].ext[:parameters][:ρ_H2] = ADMM["ρ"]["H2"][end]
    end
    if mdict[m].ext[:parameters][:H2CN_prod] == 1
        mdict[m].ext[:parameters][:gHCN_bar] = results["h2cn_prod"][m][end] - 1/(H2CN_prod["nAgents"]+1)*ADMM["Imbalances"]["H2CN_prod"][end]
        mdict[m].ext[:parameters][:λ_H2CN_prod] = results["λ"]["H2CN_prod"][end] 
        mdict[m].ext[:parameters][:ρ_H2CN_prod] = ADMM["ρ"]["H2CN_prod"][end]
        # Yearly
        mdict[m].ext[:parameters][:r_y_bar] = results["g_y"][m][end] - 1/(REC["nAgents"]+1)*ADMM["Imbalances"]["REC_y"][end]
        mdict[m].ext[:parameters][:λ_y_REC] = results["λ"]["REC_y"][end] 
        mdict[m].ext[:parameters][:ρ_y_REC] = ADMM["ρ"]["REC_y"][end]
        # Daily
        mdict[m].ext[:parameters][:r_d_bar] = results["g_d"][m][end] - 1/(REC["nAgents"]+1)*ADMM["Imbalances"]["REC_d"][end]
        mdict[m].ext[:parameters][:λ_d_REC] = results["λ"]["REC_d"][end] 
        mdict[m].ext[:parameters][:ρ_d_REC] = ADMM["ρ"]["REC_d"][end]
        # Hourly
        mdict[m].ext[:parameters][:r_h_bar] = results["g"][m][end] - 1/(REC["nAgents"]+1)*ADMM["Imbalances"]["REC_h"][end]
        mdict[m].ext[:parameters][:λ_h_REC] = results["λ"]["REC_h"][end] 
        mdict[m].ext[:parameters][:ρ_h_REC] = ADMM["ρ"]["REC_h"][end]
    end
    if mdict[m].ext[:parameters][:H2CN_cap] == 1
        mdict[m].ext[:parameters][:capHCN_bar] = results["h2cn_cap"][m][end] - 1/(H2CN_cap["nAgents"]+1)*ADMM["Imbalances"]["H2CN_cap"][end]
        mdict[m].ext[:parameters][:λ_H2CN_cap] = results[ "λ"]["H2CN_cap"][end] 
        mdict[m].ext[:parameters][:ρ_H2CN_cap] = ADMM["ρ"]["H2CN_cap"][end]
    end
    if mdict[m].ext[:parameters][:NG] == 1 # NG is not yet responsive to changes in demand - could be done at later stage
        mdict[m].ext[:parameters][:λ_NG] =  results[ "λ"]["NG"][end] 
    end
end

# Solve agents decision problems:
if m in agents[:ind]
    @timeit TO_local "Solve industry" begin
        update_ind_emissions!(mdict[m],merge(data["General"],data["Industry"]),ETS,scenario_overview_row) 
        solve_ind_agent!(mdict[m])  
    end
elseif m in agents[:ps]
    @timeit TO_local "Solve power sector" begin
        solve_ps_agent!(mdict[m])  
    end
elseif m in agents[:h2s]
    @timeit TO_local "Solve hydrogen sector" begin
        solve_h2s_agent!(mdict[m])  
    end
end

# Query results
@timeit TO_local "Query results" begin
    if mdict[m].ext[:parameters][:ETS] == 1
        push!(results["b"][m], collect(value.(mdict[m].ext[:variables][:b])))
        if m in agents[:ind]  
            push!(results["e"][m], mdict[m].ext[:parameters][:e][:])
        else
            push!(results["e"][m], collect(value.(mdict[m].ext[:expressions][:e])))
        end
    end
    if mdict[m].ext[:parameters][:EOM] == 1
        push!(results["g"][m], collect(value.(mdict[m].ext[:variables][:g])))
        push!(results["g_y"][m], collect(value.(mdict[m].ext[:expressions][:g_y])))
        push!(results["g_d"][m], collect(value.(mdict[m].ext[:expressions][:g_d])))
    end
    if mdict[m].ext[:parameters][:REC] == 1
        push!(results["r_y"][m], collect(value.(mdict[m].ext[:variables][:r_y])))
        push!(results["r_d"][m], collect(value.(mdict[m].ext[:variables][:r_d])))
        push!(results["r_h"][m], collect(value.(mdict[m].ext[:variables][:r_h])))
    end
    if mdict[m].ext[:parameters][:H2] == 1
        push!(results["h2"][m], collect(value.(mdict[m].ext[:variables][:gH])))
    end                     
    if mdict[m].ext[:parameters][:H2CN_prod] == 1
        push!(results["h2cn_prod"][m], collect(value.(mdict[m].ext[:variables][:gHCN])))
    end
    if mdict[m].ext[:parameters][:H2CN_cap] == 1
        push!(results["h2cn_cap"][m], collect(value.(mdict[m].ext[:variables][:capHCN])))
    end
end

# Merge local TO with TO:
merge!(TO,TO_local)
end