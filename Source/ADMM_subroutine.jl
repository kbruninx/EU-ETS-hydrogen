function ADMM_subroutine!(m::String,data::Dict,results::Dict,ADMM::Dict,ETS::Dict,EOM::Dict,REC::Dict,H2::Dict,H2CN_prod::Dict,H2CN_cap::Dict,NG::Dict,mod::Model,agents::Dict,TO::TimerOutput)
TO_local = TimerOutput()
# Calculate penalty terms ADMM and update price to most recent value 
@timeit TO_local "Compute ADMM penalty terms" begin
    if mod.ext[:parameters][:ETS] == 1
        mod.ext[:parameters][:b_bar] = results["b"][m][end] + 1/(ETS["nAgents"]+1)*ADMM["Imbalances"]["ETS"][end]
        mod.ext[:parameters][:λ_EUA] = results["λ"]["EUA"][end] 
        mod.ext[:parameters][:ρ_EUA] = ADMM["ρ"]["EUA"][end]
    end
    if mod.ext[:parameters][:EOM] == 1
        mod.ext[:parameters][:g_bar] = results["g"][m][end] - 1/(EOM["nAgents"]+1)*ADMM["Imbalances"]["EOM"][end]
        mod.ext[:parameters][:λ_EOM] = results["λ"]["EOM"][end] 
        mod.ext[:parameters][:ρ_EOM] = ADMM["ρ"]["EOM"][end]
    end
    if mod.ext[:parameters][:REC] == 1  
        # Yearly
        mod.ext[:parameters][:r_y_bar] = results["r_y"][m][end] - 1/(REC["nAgents"]+H2CN_prod["nAgents"]+1)*ADMM["Imbalances"]["REC_y"][end]
        mod.ext[:parameters][:λ_y_REC] = results["λ"]["REC_y"][end] 
        mod.ext[:parameters][:ρ_y_REC] = ADMM["ρ"]["REC_y"][end]
        mod.ext[:parameters][:ρ_y_REC_pre2030] = ADMM["ρ"]["REC_y_pre2030"][end]
        mod.ext[:parameters][:ρ_y_REC_post2030] = ADMM["ρ"]["REC_y_post2030"][end]
        # Monthly
        mod.ext[:parameters][:r_m_bar] = results["r_m"][m][end] - 1/(REC["nAgents"]+H2CN_prod["nAgents"]+1)*ADMM["Imbalances"]["REC_m"][end]
        mod.ext[:parameters][:λ_m_REC] = results["λ"]["REC_m"][end] 
        mod.ext[:parameters][:ρ_m_REC_pre2030] = ADMM["ρ"]["REC_m_pre2030"][end]
        mod.ext[:parameters][:ρ_m_REC_post2030] = ADMM["ρ"]["REC_m_post2030"][end]
        # Daily
        mod.ext[:parameters][:r_d_bar] = results["r_d"][m][end] - 1/(REC["nAgents"]+H2CN_prod["nAgents"]+1)*ADMM["Imbalances"]["REC_d"][end]
        mod.ext[:parameters][:λ_d_REC] = results["λ"]["REC_d"][end] 
        mod.ext[:parameters][:ρ_d_REC_pre2030] = ADMM["ρ"]["REC_d_pre2030"][end]
        mod.ext[:parameters][:ρ_d_REC_post2030] = ADMM["ρ"]["REC_d_post2030"][end]
        # Hourly
        mod.ext[:parameters][:r_h_bar] = results["r_h"][m][end] - 1/(REC["nAgents"]+H2CN_prod["nAgents"]+1)*ADMM["Imbalances"]["REC_h"][end]
        mod.ext[:parameters][:λ_h_REC] = results["λ"]["REC_h"][end] 
        mod.ext[:parameters][:ρ_h_REC_pre2030] = ADMM["ρ"]["REC_h_pre2030"][end]
        mod.ext[:parameters][:ρ_h_REC_post2030] = ADMM["ρ"]["REC_h_post2030"][end]
    end
    if mod.ext[:parameters][:H2] == 1
        mod.ext[:parameters][:gH_h_bar] = results["h2_h"][m][end] - 1/(H2["nAgents"]+1)*ADMM["Imbalances"]["H2_h"][end]
        mod.ext[:parameters][:λ_h_H2] = results["λ"]["H2_h"][end] 
        mod.ext[:parameters][:ρ_h_H2] = ADMM["ρ"]["H2_h"][end]
        mod.ext[:parameters][:gH_d_bar] = results["h2_d"][m][end] - 1/(H2["nAgents"]+1)*ADMM["Imbalances"]["H2_d"][end]
        mod.ext[:parameters][:λ_d_H2] = results["λ"]["H2_d"][end] 
        mod.ext[:parameters][:ρ_d_H2] = ADMM["ρ"]["H2_d"][end]
        mod.ext[:parameters][:gH_m_bar] = results["h2_m"][m][end] - 1/(H2["nAgents"]+1)*ADMM["Imbalances"]["H2_m"][end]
        mod.ext[:parameters][:λ_m_H2] = results["λ"]["H2_m"][end] 
        mod.ext[:parameters][:ρ_m_H2] = ADMM["ρ"]["H2_m"][end]
        mod.ext[:parameters][:gH_y_bar] = results["h2_y"][m][end] - 1/(H2["nAgents"]+1)*ADMM["Imbalances"]["H2_y"][end]
        mod.ext[:parameters][:λ_y_H2] = results["λ"]["H2_y"][end] 
        mod.ext[:parameters][:ρ_y_H2] = ADMM["ρ"]["H2_y"][end]
    end
    if mod.ext[:parameters][:H2CN_prod] == 1
        mod.ext[:parameters][:gHCN_bar] = results["h2cn_prod"][m][end] - 1/(H2CN_prod["nAgents"]+1)*ADMM["Imbalances"]["H2CN_prod"][end]
        mod.ext[:parameters][:λ_H2CN_prod] = results["λ"]["H2CN_prod"][end] 
        mod.ext[:parameters][:ρ_H2CN_prod] = ADMM["ρ"]["H2CN_prod"][end]
    end
    if mod.ext[:parameters][:H2CN_cap] == 1
        mod.ext[:parameters][:capHCN_bar] = results["h2cn_cap"][m][end] - 1/(H2CN_cap["nAgents"]+1)*ADMM["Imbalances"]["H2CN_cap"][end]
        mod.ext[:parameters][:λ_H2CN_cap] = results["λ"]["H2CN_cap"][end] 
        mod.ext[:parameters][:ρ_H2CN_cap] = ADMM["ρ"]["H2CN_cap"][end]
    end
    if mod.ext[:parameters][:NG] == 1 # NG is not yet responsive to changes in demand - could be done at later stage
        mod.ext[:parameters][:λ_NG] =  results[ "λ"]["NG"][end] 
    end
end

# Solve agents decision problems:
if m in agents[:ind]
    @timeit TO_local "Solve industry" begin
        update_ind_emissions!(mod,merge(data["General"],data["Industry"],data["scenario"]),ETS) 
        solve_ind_agent!(mod)  
    end
elseif m in agents[:ps]
    @timeit TO_local "Solve power sector" begin
        solve_ps_agent!(mod)  
    end
elseif m in agents[:h2s]
    @timeit TO_local "Solve hydrogen sector" begin
        solve_h2s_agent!(mod)  
    end
end

# Query results
@timeit TO_local "Query results" begin
    if mod.ext[:parameters][:ETS] == 1
        push!(results["b"][m], collect(value.(mod.ext[:variables][:b])))
        if m in agents[:ind]  
            push!(results["e"][m], mod.ext[:parameters][:e][:])
        else
            push!(results["e"][m], collect(value.(mod.ext[:expressions][:e])))
        end
    end
    if mod.ext[:parameters][:EOM] == 1
        push!(results["g"][m], collect(value.(mod.ext[:variables][:g])))
    end
    if mod.ext[:parameters][:REC] == 1
        push!(results["r_y"][m], collect(value.(mod.ext[:variables][:r_y])))
        push!(results["r_m"][m], collect(value.(mod.ext[:variables][:r_m])))
        push!(results["r_d"][m], collect(value.(mod.ext[:variables][:r_d])))
        push!(results["r_h"][m], collect(value.(mod.ext[:variables][:r_h])))
    end
    if mod.ext[:parameters][:H2] == 1
        push!(results["h2_h"][m], collect(value.(mod.ext[:variables][:gH])))
        push!(results["h2_d"][m], collect(value.(mod.ext[:expressions][:gH_d])))
        push!(results["h2_m"][m], collect(value.(mod.ext[:variables][:gH_m])))
        push!(results["h2_y"][m], collect(value.(mod.ext[:expressions][:gH_y])))
    end                     
    if mod.ext[:parameters][:H2CN_prod] == 1
        push!(results["h2cn_prod"][m], collect(value.(mod.ext[:variables][:gHCN])))
    end
    if mod.ext[:parameters][:H2CN_cap] == 1
        push!(results["h2cn_cap"][m], collect(value.(mod.ext[:variables][:capHCN])))
    end
end

# Merge local TO with TO:
merge!(TO,TO_local)
end