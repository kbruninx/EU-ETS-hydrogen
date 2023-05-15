# Packages 
using DataFrames, CSV, YAML, DataStructures # dataprocessing

# Constants
const home_dir = @__DIR__

# Functions 
# Remove duplicate rows, keep last entry, modified from: https://discourse.julialang.org/t/delete-duplicate-rows-in-a-dataframe/65487/10
function dropduplicates(df, cols; keep = :first) 
    keep in [:first, :last] || error("keep parameter should be :first or :last")
    combine(groupby(df, cols)) do sdf
        if nrow(sdf) == 1 
            DataFrame(
                sdf
            )
        else
            DataFrame(
              filter(
                r->rownumber(r)==(keep == :first ? 1 : nrow(sdf)), 
                eachrow(sdf)
              )
            )
        end
    end
end

# Read overview results file and clean up 
df = CSV.read(joinpath(home_dir,"overview_results_8_repr_days.csv"),delim=";",DataFrame)
df_filtered = dropduplicates(df, [:scen_number, :sensitivity], keep = :last)
sort!(df_filtered, [:scen_number])

# # Limit on number of iterations hit
df_walltime_hit = filter(row -> row.n_iter == 10000, df_filtered)
CSV.write(joinpath(home_dir,"overview_results_8_repr_days_ITERATION_LIMIT.csv"),delim=";",df_walltime_hit)

# Additional results added to overview results
Years= range(2021,2060,step=1)
waterbed_sealed = zeros(nrow(df_filtered))
CumulativeEmissions_PS= zeros(nrow(df_filtered))
CumulativeEmissions_H2S= zeros(nrow(df_filtered))
CumulativeEmissions_IND= zeros(nrow(df_filtered))
nodenames = ["scen_number";"sens";string.(Years)]
RES = DataFrame([[] for _ = nodenames] , nodenames)
PriceElectricity = DataFrame([[] for _ = nodenames] , nodenames)

for n = 1:nrow(df_filtered)
    ets_temp = CSV.read(joinpath(home_dir,string("Results_8_repr_days"),string("Scenario_",df_filtered[n,:scen_number],"_ETS_",df_filtered[n,:sensitivity],".csv")), delim=";",DataFrame)
    waterbed_sealed[n] = Years[findfirst(x -> x<833, ets_temp[:,:TNAC])]
    CumulativeEmissions_PS[n] = sum(ets_temp[:,:Emissions_PS])
    CumulativeEmissions_IND[n] = sum(ets_temp[:,:Emissions_Ind])
    CumulativeEmissions_H2S[n] = sum(ets_temp[:,:Emissions_H2S])

    H2_temp = CSV.read(joinpath(home_dir,string("Results_8_repr_days"),string("Scenario_",df_filtered[n,:scen_number],"_H2_",df_filtered[n,:sensitivity],".csv")), delim=";",DataFrame)
    # which aggregate metrics do we need to look at?

    PS_temp = CSV.read(joinpath(home_dir,string("Results_8_repr_days"),string("Scenario_",df_filtered[n,:scen_number],"_PS_",df_filtered[n,:sensitivity],".csv")), delim=";",DataFrame)
    FS_RES = transpose(PS_temp[:,:FS_Biomass]+ PS_temp[:,:FS_Solar] + PS_temp[:,:FS_WindOnshore]+PS_temp[:,:FS_WindOffshore])
    append!(RES,DataFrame([df_filtered[n,:scen_number] df_filtered[n,:sensitivity] FS_RES],nodenames))
    append!(RES,DataFrame([df_filtered[n,:scen_number] df_filtered[n,:sensitivity] transpose(PS_temp[:,:EOM_avg])],nodenames))

end
df_filtered.waterbed_sealed = waterbed_sealed;
df_filtered.CumulativeEmissions_PS = CumulativeEmissions_PS;
df_filtered.CumulativeEmissions_IND = CumulativeEmissions_IND;
df_filtered.CumulativeEmissions_H2S = CumulativeEmissions_H2S;

# Write out results
CSV.write(joinpath(home_dir,"overview_results_8_repr_days_CLEAN.csv"),delim=";",df_filtered)
CSV.write(joinpath(home_dir,"FS_RES_8_repr_days.csv"),delim=";",RES)


