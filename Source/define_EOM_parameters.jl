function define_EOM_parameters!(EOM::Dict,data::Dict,ts::DataFrame,repr_days::DataFrame)
    # Growth factor of demand based on 2016 reference scenario
    EOM["GF_YoY"] = data["SF"].*[0; 0.014*ones(10); 0.007*ones(data["nyears"]-11)]  # growth rate according to EU reference scenario 2020
    EOM["GF"] = [sum(EOM["GF_YoY"][1:jy]) for jy=1:data["nyears"]]

    # timeseries
    EOM["D_tot"] = [ts[!,:LOAD][round(Int,data["nTimesteps"]*repr_days[!,:periods][jd]+jh)]/10^6 for jh=1:data["nTimesteps"], jd=1:data["nReprDays"]] # from MWh to TWh
    EOM["G_other"] = [ts[!,:OTHER][round(Int,data["nTimesteps"]*repr_days[!,:periods][jd]+jh)]/10^6 for jh=1:data["nTimesteps"], jd=1:data["nReprDays"]] # from MWh to TWh
    EOM["D"] = [(1+EOM["GF"][jy])*EOM["D_tot"][jh,jd] - EOM["G_other"][jh,jd] for jh=1:data["nTimesteps"], jd=1:data["nReprDays"], jy=1:data["nyears"]]

    # Total and weighted
    EOM["W"] = W = Dict(jd => repr_days[!,:weights][jd] for jd=1:data["nReprDays"])
    EOM["D_cum"] = [sum(W[jd]*EOM["D"][jh,jd,jy] for jh=1:data["nTimesteps"],jd=1:data["nReprDays"]) for jy=1:data["nyears"]]
    EOM["Dw"] = [W[jd]*EOM["D"][jh,jd,jy] for jh=1:data["nTimesteps"], jd=1:data["nReprDays"], jy=1:data["nyears"]]

    return EOM
end