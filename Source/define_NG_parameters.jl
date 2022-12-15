function define_NG_parameters!(NG::Dict,data::Dict,ts::DataFrame,repr_days::DataFrame)
    YoY = data["YoY"]/100*data["P_2020"]

    NG["λ"] = [data["P_2019"]; data["P_2020"]; data["P_2021"]; data["P_2022"]; data["P_2023"];data["P_2024"];data["P_2025"];data["P_2026"];data["P_2027"];data["P_2028"];data["P_2029"]; [data["P_2030"]+YoY*(jy-1) for jy in 1:data["nyears"]-11]]
    
    return NG
end