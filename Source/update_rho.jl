function update_rho!(ADMM::Dict, iter::Int64)
    if mod(iter,1) == 0
        # ρ-updates following Boyd et al. (2010)
        if ADMM["Residuals"]["Primal"]["ETS"][end]> 2*ADMM["Residuals"]["Dual"]["ETS"][end]
            push!(ADMM["ρ"]["EUA"], 1.1*ADMM["ρ"]["EUA"][end])
        elseif ADMM["Residuals"]["Dual"]["ETS"][end] > 2*ADMM["Residuals"]["Primal"]["ETS"][end]
            push!(ADMM["ρ"]["EUA"], 1/1.1*ADMM["ρ"]["EUA"][end])
        end

        if ADMM["Residuals"]["Primal"]["EOM"][end] > 2*ADMM["Residuals"]["Dual"]["EOM"][end]
            push!(ADMM["ρ"]["EOM"], 1.1*ADMM["ρ"]["EOM"][end])
        elseif ADMM["Residuals"]["Dual"]["EOM"][end] > 2*ADMM["Residuals"]["Primal"]["EOM"][end]
            push!(ADMM["ρ"]["EOM"], 1/1.1*ADMM["ρ"]["EOM"][end])
        end

        if ADMM["Residuals"]["Primal"]["REC"][end] > 2*ADMM["Residuals"]["Dual"]["REC"][end]
            push!(ADMM["ρ"]["REC"], 1.1*ADMM["ρ"]["REC"][end])
        elseif ADMM["Residuals"]["Dual"]["EOM"][end] > 2*ADMM["Residuals"]["Primal"]["REC"][end]
            push!(ADMM["ρ"]["REC"], 1/1.1*ADMM["ρ"]["REC"][end])
        end
    end
end