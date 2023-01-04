function update_rho!(ADMM::Dict, iter::Int64)
    if mod(iter,1) == 0
        # ρ-updates following Boyd et al.  
        if ADMM["Residuals"]["Primal"]["ETS"][end]+ADMM["Residuals"]["Primal"]["MSR"][end]> 2*ADMM["Residuals"]["Dual"]["ETS"][end]
            push!(ADMM["ρ"]["EUA"], minimum([1000,1.1*ADMM["ρ"]["EUA"][end]]))
        elseif ADMM["Residuals"]["Dual"]["ETS"][end]+ADMM["Residuals"]["Primal"]["MSR"][end] > 2*ADMM["Residuals"]["Primal"]["ETS"][end]
            push!(ADMM["ρ"]["EUA"], 1/1.1*ADMM["ρ"]["EUA"][end])
        end

        if ADMM["Residuals"]["Primal"]["EOM"][end] > 2*ADMM["Residuals"]["Dual"]["EOM"][end]
            push!(ADMM["ρ"]["EOM"], minimum([1000,1.1*ADMM["ρ"]["EOM"][end]]))
        elseif ADMM["Residuals"]["Dual"]["EOM"][end] > 2*ADMM["Residuals"]["Primal"]["EOM"][end]
            push!(ADMM["ρ"]["EOM"], 1/1.1*ADMM["ρ"]["EOM"][end])
        end

        if ADMM["Residuals"]["Primal"]["REC"][end] > 2*ADMM["Residuals"]["Dual"]["REC"][end]
            push!(ADMM["ρ"]["REC_y"], minimum([1000,1.1*ADMM["ρ"]["REC_y"][end]]))
            push!(ADMM["ρ"]["REC_m"], minimum([1000,1.1*ADMM["ρ"]["REC_m"][end]]))
            push!(ADMM["ρ"]["REC_d"], minimum([1000,1.1*ADMM["ρ"]["REC_d"][end]]))
            push!(ADMM["ρ"]["REC_h"], minimum([1000,1.1*ADMM["ρ"]["REC_h"][end]]))
        elseif ADMM["Residuals"]["Dual"]["EOM"][end] > 2*ADMM["Residuals"]["Primal"]["REC"][end]
            push!(ADMM["ρ"]["REC_y"], 1/1.1*ADMM["ρ"]["REC_y"][end])
            push!(ADMM["ρ"]["REC_m"], 1/1.1*ADMM["ρ"]["REC_m"][end])
            push!(ADMM["ρ"]["REC_d"], 1/1.1*ADMM["ρ"]["REC_d"][end])
            push!(ADMM["ρ"]["REC_h"], 1/1.1*ADMM["ρ"]["REC_h"][end])
        end

        if ADMM["Residuals"]["Primal"]["H2"][end] > 2*ADMM["Residuals"]["Dual"]["H2"][end]
            push!(ADMM["ρ"]["H2"], minimum([1000,1.1*ADMM["ρ"]["H2"][end]]))
        elseif ADMM["Residuals"]["Dual"]["H2"][end] > 2*ADMM["Residuals"]["Primal"]["H2"][end]
            push!(ADMM["ρ"]["H2"], 1/1.1*ADMM["ρ"]["H2"][end])
        end

        if ADMM["Residuals"]["Primal"]["H2CN_prod"][end] > 2*ADMM["Residuals"]["Dual"]["H2CN_prod"][end]
            push!(ADMM["ρ"]["H2CN_prod"], minimum([1000,1.1*ADMM["ρ"]["H2CN_prod"][end]]))
        elseif ADMM["Residuals"]["Dual"]["H2CN_prod"][end] > 2*ADMM["Residuals"]["Primal"]["H2CN_prod"][end]
            push!(ADMM["ρ"]["H2CN_prod"], 1/1.1*ADMM["ρ"]["H2CN_prod"][end])
        end

        if ADMM["Residuals"]["Primal"]["H2CN_cap"][end] > 2*ADMM["Residuals"]["Dual"]["H2CN_cap"][end]
            push!(ADMM["ρ"]["H2CN_cap"], minimum([1000,1.1*ADMM["ρ"]["H2CN_cap"][end]]))
        elseif ADMM["Residuals"]["Dual"]["H2CN_cap"][end] > 2*ADMM["Residuals"]["Primal"]["H2CN_cap"][end]
            push!(ADMM["ρ"]["H2CN_cap"], 1/1.1*ADMM["ρ"]["H2CN_cap"][end])
        end
    end
end