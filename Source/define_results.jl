function define_results!(data::Dict,results::Dict,ADMM::Dict,agents::Dict,ETS::Dict,EOM::Dict,REC::Dict) 
    results["e"] = Dict()
    results["g"] = Dict()
    results["r"] = Dict()
    results["b"] = Dict()
    for m in agents[:ets]
        results["b"][m] = CircularBuffer{Array{Float64,1}}(data["CircularBufferSize"])  
        push!(results["b"][m],zeros(data["nyears"]))
        results["e"][m] = CircularBuffer{Array{Float64,1}}(data["CircularBufferSize"])  
        push!(results["e"][m],zeros(data["nyears"]))
    end
    for m in agents[:rec]
        results["r"][m] = CircularBuffer{Array{Float64,1}}(data["CircularBufferSize"])  
        push!(results["r"][m],zeros(data["nyears"]))
    end
    for m in agents[:eom]
        results["g"][m] = CircularBuffer{Array{Float64,3}}(data["CircularBufferSize"]) 
        push!(results["g"][m],zeros(data["nTimesteps"],data["nReprDays"],data["nyears"]))
    end
    results["s"] = CircularBuffer{Array{Float64,1}}(data["CircularBufferSize"])  
    push!(results["s"],zeros(data["nyears"]))

    results["λ"] = Dict()

    results[ "λ"]["EUA"] = CircularBuffer{Array{Float64,1}}(data["CircularBufferSize"]) 
    push!(results[ "λ"]["EUA"],[5; 5; zeros(data["nyears"]-2)])

    results[ "λ"]["REC"] = CircularBuffer{Array{Float64,1}}(data["CircularBufferSize"])  
    push!(results[ "λ"]["REC"],zeros(data["nyears"]))
    results[ "λ"]["EOM"] = CircularBuffer{Array{Float64,3}}(data["CircularBufferSize"]) 
    push!(results[ "λ"]["EOM"],zeros(data["nTimesteps"],data["nReprDays"],data["nyears"]))

    ADMM["Imbalances"] = Dict()
    ADMM["Imbalances"]["ETS"] = CircularBuffer{Array{Float64,1}}(data["CircularBufferSize"])  
    push!(ADMM["Imbalances"]["ETS"],zeros(data["nyears"]))
    ADMM["Imbalances"]["MSR"] = CircularBuffer{Array{Float64,1}}(data["CircularBufferSize"])
    push!(ADMM["Imbalances"]["MSR"],zeros(data["nyears"]))
    ADMM["Imbalances"]["EOM"] = CircularBuffer{Array{Float64,3}}(data["CircularBufferSize"])
    push!(ADMM["Imbalances"]["EOM"],zeros(data["nTimesteps"],data["nReprDays"],data["nyears"]))
    ADMM["Imbalances"]["REC"] = CircularBuffer{Array{Float64,1}}(data["CircularBufferSize"])
    push!(ADMM["Imbalances"]["REC"],zeros(data["nyears"]))

    ADMM["Residuals"] = Dict()
    ADMM["Residuals"]["Primal"] = Dict()
    ADMM["Residuals"]["Primal"]["ETS"] = CircularBuffer{Float64}(data["CircularBufferSize"])
    push!(ADMM["Residuals"]["Primal"]["ETS"],0)
    ADMM["Residuals"]["Primal"]["MSR"] = CircularBuffer{Float64}(data["CircularBufferSize"])
    push!(ADMM["Residuals"]["Primal"]["MSR"],0)
    ADMM["Residuals"]["Primal"]["REC"] = CircularBuffer{Float64}(data["CircularBufferSize"])
    push!(ADMM["Residuals"]["Primal"]["REC"],0)
    ADMM["Residuals"]["Primal"]["EOM"] = CircularBuffer{Float64}(data["CircularBufferSize"])
    push!(ADMM["Residuals"]["Primal"]["EOM"],0)
    ADMM["Residuals"]["Dual"] = Dict()
    ADMM["Residuals"]["Dual"]["ETS"] = CircularBuffer{Float64}(data["CircularBufferSize"])
    push!(ADMM["Residuals"]["Dual"]["ETS"],0)
    ADMM["Residuals"]["Dual"]["REC"] = CircularBuffer{Float64}(data["CircularBufferSize"])
    push!(ADMM["Residuals"]["Dual"]["REC"],0)
    ADMM["Residuals"]["Dual"]["EOM"] = CircularBuffer{Float64}(data["CircularBufferSize"])
    push!(ADMM["Residuals"]["Dual"]["EOM"],0)

    ADMM["Tolerance"] = Dict()
    ADMM["Tolerance"]["ETS"] = data["epsilon"]/100*sum(ETS["CAP"])
    ADMM["Tolerance"]["EOM"] = data["epsilon"]/100*sum(EOM["D"])
    ADMM["Tolerance"]["REC"] = data["epsilon"]/100*sum(REC["RT"].*EOM["D_cum"])
    
    ADMM["ρ"] = Dict()
    ADMM["ρ"]["EUA"] = CircularBuffer{Float64}(data["CircularBufferSize"]) 
    push!(ADMM["ρ"]["EUA"],data["rho_EUA"])
    ADMM["ρ"]["EOM"] = CircularBuffer{Float64}(data["CircularBufferSize"]) 
    push!(ADMM["ρ"]["EOM"],data["rho_EOM"])
    ADMM["ρ"]["REC"] = CircularBuffer{Float64}(data["CircularBufferSize"]) 
    push!(ADMM["ρ"]["REC"],data["rho_REC"])
    
    ADMM["n_iter"] = 1 
    ADMM["walltime"] = 0
    
    return results, ADMM
end