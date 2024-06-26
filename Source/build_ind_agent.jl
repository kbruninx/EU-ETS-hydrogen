function build_ind_agent!(mod)
# Extract sets
JY = mod.ext[:sets][:JY]

# Extract parameters
A = mod.ext[:parameters][:A]
e = mod.ext[:parameters][:e]
AC = mod.ext[:parameters][:AC]
λ_EUA = mod.ext[:parameters][:λ_EUA]
b_bar = mod.ext[:parameters][:b_bar]
ρ_EUA = mod.ext[:parameters][:ρ_EUA]

# Define variables
b = mod.ext[:variables][:b] = @variable(mod, [jy=JY], lower_bound = 0, base_name="EUA") 

# Expressions
mod.ext[:expressions][:tot_cost] = @expression(mod, 
    sum(A[jy]*AC[jy] for jy in JY)
)

# Objective
@objective(mod, Min, sum(A[jy]*λ_EUA[jy]*b[jy] for jy in JY) +sum(ρ_EUA/2*(b[jy] - b_bar[jy])^2 for jy in JY))

# Constraints 
mod.ext[:constraints][:con1]  = @constraint(mod,[jy=JY], 
    sum(b[y2] for y2=1:jy) >= sum(e[y2] for y2=1:jy)
)

optimize!(mod);

return mod
end
