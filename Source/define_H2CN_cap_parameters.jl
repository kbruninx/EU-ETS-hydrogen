function define_H2CN_cap_parameters!(H2CN_cap::Dict,data::Dict,ts::DataFrame,repr_days::DataFrame,scenario_overview_row::DataFrameRow)
    # carbon neutral hydrogen production capacity target
    H2CN_cap["H2CN_CAPT"] = [zeros(5); scenario_overview_row["CNH2_cap_target_2024"]*ones(6); scenario_overview_row["CNH2_cap_target_2030"]*ones(data["nyears"]-11)]

    return H2CN_cap
end