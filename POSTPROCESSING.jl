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
df_filtered = dropduplicates!(df, [:scen_number, :sensitivity], keep = :last)
sort!(df_filtered, [:scen_number])
CSV.write(joinpath(home_dir,"overview_results_8_repr_days_CLEAN.csv"),delim=";",df_filtered)


# Limit on number of iterations hit
df_walltime_hit = filter!(row -> row.n_iter == 10000, df_filtered)
CSV.write(joinpath(home_dir,"overview_results_8_repr_days_ITERATION_LIMIT.csv"),delim=";",df_filtered)



