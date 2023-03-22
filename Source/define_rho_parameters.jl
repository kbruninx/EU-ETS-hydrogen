function define_rho_parameters!(data)
    # Additionality 
    if data["scenario"]["Additionality_pre_2030"] == "Yearly" 
        data["ADMM"]["rho_REC_y_pre2030"] = data["ADMM"]["rho_REC"]
    else
        data["ADMM"]["rho_REC_y_pre2030"] = 0
    end 
    if data["scenario"]["Additionality_pre_2030"] == "Monthly" 
        data["ADMM"]["rho_REC_m_pre2030"] = data["ADMM"]["rho_REC"]
    else
        data["ADMM"]["rho_REC_m_pre2030"] = 0
    end 
    if data["scenario"]["Additionality_pre_2030"] == "Daily" 
        data["ADMM"]["rho_REC_d_pre2030"] = data["ADMM"]["rho_REC"]
    else
        data["ADMM"]["rho_REC_d_pre2030"] = 0
    end 
    if data["scenario"]["Additionality_pre_2030"] == "Hourly" 
        data["ADMM"]["rho_REC_h_pre2030"] = data["ADMM"]["rho_REC"]
    else
        data["ADMM"]["rho_REC_h_pre2030"] = 0
    end 

    if data["scenario"]["Additionality_post_2030"] == "Yearly" 
        data["ADMM"]["rho_REC_y_post2030"] = data["ADMM"]["rho_REC"]
    else
        data["ADMM"]["rho_REC_y_post2030"] = 0
    end 
    if data["scenario"]["Additionality_post_2030"] == "Monthly" 
        data["ADMM"]["rho_REC_m_post2030"] = data["ADMM"]["rho_REC"]
    else
        data["ADMM"]["rho_REC_m_post2030"] = 0
    end 
    if data["scenario"]["Additionality_post_2030"] == "Daily" 
        data["ADMM"]["rho_REC_d_post2030"] = data["ADMM"]["rho_REC"]
    else
        data["ADMM"]["rho_REC_d_post2030"] = 0
    end 
    if data["scenario"]["Additionality_post_2030"] == "Hourly" 
        data["ADMM"]["rho_REC_h_post2030"] = data["ADMM"]["rho_REC"]
    else
        data["ADMM"]["rho_REC_h_post2030"] = 0
    end 
    
    data["ADMM"]["rho_REC_y"] = data["ADMM"]["rho_REC"]

    # Hydrogen demand   
    if data["scenario"]["H2_balance"] == "Hourly" 
        data["ADMM"]["rho_H2_h"] = data["ADMM"]["rho_H2"]
    else
        data["ADMM"]["rho_H2_h"] = 0
    end
    if data["scenario"]["H2_balance"] == "Daily" 
        data["ADMM"]["rho_H2_d"] = data["ADMM"]["rho_H2"]
    else
        data["ADMM"]["rho_H2_d"] = 0
    end
    if data["scenario"]["H2_balance"] == "Monthly" 
        data["ADMM"]["rho_H2_m"] = data["ADMM"]["rho_H2"]
    else
        data["ADMM"]["rho_H2_m"] = 0
    end
    if data["scenario"]["H2_balance"] == "Yearly" 
        data["ADMM"]["rho_H2_y"] = data["ADMM"]["rho_H2"]
    else
        data["ADMM"]["rho_H2_y"] = 0
    end

    return data
end
