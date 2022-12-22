function define_NG_parameters!(NG::Dict,data::Dict,ts::DataFrame,repr_days::DataFrame)
    YoY = data["SF"]*data["YoY"]/100*data["P_2020"]

    NG["λ"] = [data["P_2021"]; data["P_2022"]; data["SF"]*data["P_2023"]; data["SF"]*data["P_2024"]; data["SF"]*data["P_2025"]; data["SF"]*data["P_2026"]; data["SF"]*data["P_2027"]; data["SF"]*data["P_2028"]; data["SF"]*data["P_2029"]; [data["SF"]*data["P_2030"]+YoY*(jy-1) for jy in 1:data["nyears"]-9]]
    
    return NG
end